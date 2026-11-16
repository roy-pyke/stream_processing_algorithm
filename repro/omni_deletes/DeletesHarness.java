// Driver for the journal version's code (github.com/wiegerpunter/omni-deletes, commit 17b969e) on
// update streams with deletions.
//
// The upstream sources are compiled UNMODIFIED (scripts/setup_omni_deletes.sh); this class only
// replaces their I/O harness (CSV stream files under hard-coded folders, DuckDB mixing) so that the
// upstream synopses see the update streams and query workloads of our C++ benchmark. The synopses
// are built exactly as in their experiment runner (RunExperiments.omniExperiment, aSHExperiment):
//   s2   OmniSketch with batched deletions (OmniSketchBuilder.QFirstSampleLaterArrayWithBufferOptBatchDeletes,
//        their setting expOmniSketchVLDBArrayWithBufferOptBatchDeletes), d = 3, w = 28, b = 31 and
//        B = (ram / (d w p) - 32) / b ("GridSearch"), ingestion and deletion buffers of
//        8e6 / (d w p) / 32 / 2 hashes each, seed = repetition, beta = 1 (their getBeta() returns 1
//        when its table paramTable/bufferMinwiseTable.csv is absent; it is not in the repository)
//   ash  adaptive Sample-and-Hold (omni.synopses.baselines.aSH), parameters from
//        RamToPar.getParamsAdapSampling(ram, 0.1), useNoCast = true, withDeletes = true
// A record is passed as long[]{id, a_1, ..., a_p} (id = row index + 1), a query as long[p] with -1
// for attributes without a predicate, as in their runQueriesNoCast.
//
// Usage:
//   java omni.DeletesHarness gen <out.bin> <numRecords> <numAttributes> <alpha>
//        the generator of their data pipeline (dataset-pipeline, SyntheticDatasetRef + ZipfGenerator):
//        commons-math ZipfDistribution over {1..10000} driven by Well19937c(0), p values per record
//   java omni.DeletesHarness run s2|ash <base.bin> <ops.bin> <workload.wl> <out.csv> <name> <repetition> <mem_mb>...
//        ops.bin from bin/mkdelstream; one line per query appended to out.csv:
//        name,mem_mb,repetition,query_index,truth,estimate
package omni;

import omni.Experiments.parameterSetting.RamToPar;
import omni.Experiments.utils.QueryInfo;
import omni.synopses.SynopsisRefactor;
import omni.synopses.baselines.aSH.aSH;
import omni.synopses.omniFactory.OmniSketch;
import omni.synopses.omniFactory.OmniSketchBuilder;
import org.apache.commons.math3.distribution.ZipfDistribution;
import org.apache.commons.math3.random.Well19937c;

import java.io.*;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.channels.FileChannel;
import java.util.ArrayList;
import java.util.List;

public class DeletesHarness {
    static int p;
    static long n;
    static int[] data;

    static void loadData(String path) throws IOException {
        try (DataInputStream in = new DataInputStream(new BufferedInputStream(new FileInputStream(path), 1 << 20))) {
            byte[] hdr = new byte[16];
            in.readFully(hdr);
            ByteBuffer bb = ByteBuffer.wrap(hdr).order(ByteOrder.LITTLE_ENDIAN);
            if (hdr[0] != 'S' || hdr[1] != 'T' || hdr[2] != 'R' || hdr[3] != 'M') throw new IOException("bad magic");
            p = bb.getInt(4);
            n = bb.getLong(8);
            long total = n * p;
            if (total > Integer.MAX_VALUE - 8) throw new IOException("dataset too large for one int[]");
            data = new int[(int) total];
            byte[] buf = new byte[1 << 22];
            int pos = 0;
            while (pos < total) {
                int want = (int) Math.min(buf.length, (total - pos) * 4);
                in.readFully(buf, 0, want);
                ByteBuffer b2 = ByteBuffer.wrap(buf, 0, want).order(ByteOrder.LITTLE_ENDIAN);
                for (int i = 0; i < want / 4; i++) data[pos++] = b2.getInt();
            }
        }
    }

    static int[] loadOps(String path) throws IOException {
        try (FileChannel ch = new FileInputStream(path).getChannel()) {
            ByteBuffer hdr = ByteBuffer.allocate(12).order(ByteOrder.LITTLE_ENDIAN);
            while (hdr.hasRemaining()) ch.read(hdr);
            hdr.flip();
            if (hdr.get(0) != 'U' || hdr.get(1) != 'P' || hdr.get(2) != 'D' || hdr.get(3) != 'S') throw new IOException("bad magic");
            long m = hdr.getLong(4);
            int[] ops = new int[(int) m];
            ByteBuffer b = ByteBuffer.allocate(1 << 22).order(ByteOrder.LITTLE_ENDIAN);
            int pos = 0;
            while (pos < m) {
                b.clear();
                b.limit((int) Math.min(b.capacity(), (m - pos) * 4));
                while (b.hasRemaining()) if (ch.read(b) < 0) throw new IOException("short update stream");
                b.flip();
                while (b.hasRemaining()) ops[pos++] = b.getInt();
            }
            return ops;
        }
    }

    static class WQ {
        long truth;
        int np;
        long[] q;
    }

    static ArrayList<WQ> loadWorkload(String path) throws IOException {
        ArrayList<WQ> qs = new ArrayList<>();
        try (BufferedReader r = new BufferedReader(new FileReader(path))) {
            String line;
            while ((line = r.readLine()) != null) {
                String[] t = line.trim().split("\\s+");
                if (t.length < 4) continue;
                WQ w = new WQ();
                w.np = Integer.parseInt(t[1]);
                w.truth = Long.parseLong(t[2]);
                w.q = new long[p];
                java.util.Arrays.fill(w.q, -1);
                for (int i = 0; i < w.np; i++) w.q[Integer.parseInt(t[4 + 2 * i])] = Long.parseLong(t[5 + 2 * i]);
                qs.add(w);
            }
        }
        return qs;
    }

    public static void main(String[] args) throws Exception {
        if (args[0].equals("gen")) {
            gen(args);
            return;
        }
        String kind = args[1];
        loadData(args[2]);
        int[] ops = loadOps(args[3]);
        ArrayList<WQ> wl = loadWorkload(args[4]);
        String outPath = args[5], name = args[6];
        int rep = Integer.parseInt(args[7]);
        System.err.println("loaded " + n + " rows, " + p + " attributes, " + ops.length + " updates, " + wl.size() + " queries");
        try (PrintWriter out = new PrintWriter(new FileWriter(outPath, true))) {
            for (int ai = 8; ai < args.length; ai++) {
                long ram = (long) (Double.parseDouble(args[ai]) * 8e6);
                SynopsisRefactor s;
                String cfg;
                if (kind.equals("s2")) {
                    int d = 3, w = 28, b = 31;
                    int B = (int) ((ram / (d * w * p) - 32) / b);
                    int bufferSize = 8000000 / (d * w * p) / 32;
                    s = OmniSketchBuilder.QFirstSampleLaterArrayWithBufferOptBatchDeletes()
                            .setRam(ram)
                            .setIngestBufferSize(bufferSize / 2)
                            .setDeleteBufferSize(bufferSize / 2)
                            .setNumStoredAttributes(p)
                            .setParams(new int[]{d, w, B, b})
                            .setSeed(rep)
                            .setBufferBeta(1.0)
                            .build();
                    cfg = "d=" + d + ";w=" + w + ";B=" + B + ";b=" + b + ";buf=" + bufferSize / 2;
                } else if (kind.equals("ash")) {
                    RamToPar rtp = new RamToPar(p, new ArrayList<>(List.of(ram)));
                    int[] params = rtp.getParamsAdapSampling(ram, 0.1);
                    s = new aSH(ram, true, p, params, rep, true, 0.1);
                    cfg = "sample=" + params[0] + ";buffer=" + params[1];
                } else {
                    throw new IllegalArgumentException("unknown synopsis " + kind);
                }
                long t0 = System.nanoTime();
                long live = 0, dels = 0;
                for (int op : ops) {
                    int i = (op > 0 ? op : -op) - 1;
                    long[] rec = new long[p + 1];
                    rec[0] = i + 1;
                    int base = i * p;
                    for (int a = 0; a < p; a++) rec[a + 1] = data[base + a] & 0xFFFFFFFFL;
                    if (op > 0) {
                        s.add(rec);
                        live++;
                    } else {
                        s.delete(rec);
                        live--;
                        dels++;
                    }
                }
                double tIngest = (System.nanoTime() - t0) / 1e9;
                long memBits = s.getMemoryUsage();
                long t1 = System.nanoTime();
                StringBuilder sb = new StringBuilder();
                for (int qi = 0; qi < wl.size(); qi++) {
                    WQ w = wl.get(qi);
                    int est = s.query(w.q.clone(), w.np, new QueryInfo());
                    sb.append(name).append(',').append(args[ai]).append(',').append(rep).append(',').append(qi).append(',')
                            .append(w.truth).append(',').append(est).append('\n');
                }
                double tQuery = (System.nanoTime() - t1) / 1e9;
                out.print(sb);
                out.flush();
                System.out.printf("JDRUN,%s,%s,%s,%d,%d,%d,%d,%.3f,%.2f,%.2f,%.3f,%s%n", name, kind, args[ai], rep, p, ops.length, live,
                        tIngest, tIngest / ops.length * 1e9, tQuery / wl.size() * 1e6, memBits / 8e6, cfg + ";deletes=" + dels);
                s.reset();
            }
        }
    }

    static void gen(String[] args) throws IOException {
        long nRec = Long.parseLong(args[2]);
        int nAttr = Integer.parseInt(args[3]);
        double alpha = Double.parseDouble(args[4]);
        ZipfDistribution zipf = new ZipfDistribution(new Well19937c(0), 10000, alpha);
        try (BufferedOutputStream bo = new BufferedOutputStream(new FileOutputStream(args[1]), 1 << 22)) {
            ByteBuffer hdr = ByteBuffer.allocate(16).order(ByteOrder.LITTLE_ENDIAN);
            hdr.put((byte) 'S').put((byte) 'T').put((byte) 'R').put((byte) 'M').putInt(nAttr).putLong(nRec);
            bo.write(hdr.array());
            ByteBuffer row = ByteBuffer.allocate(4 * nAttr).order(ByteOrder.LITTLE_ENDIAN);
            for (long i = 0; i < nRec; i++) {
                row.clear();
                for (int a = 0; a < nAttr; a++) row.putInt(zipf.sample());
                bo.write(row.array());
            }
        }
    }
}
