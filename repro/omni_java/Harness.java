// Driver for the authors' OmniSketch implementation (github.com/wiegerpunter/omnisketch).
//
// The upstream sources are compiled UNMODIFIED; this class only replaces their I/O
// harness (hard-coded paths, CSV dataset files, CAIDA/SNMP parsers) so that the
// upstream sketch can be fed the binary datasets and query workloads used by our
// C++ benchmark. Everything that determines the estimates is upstream code:
//   - parameter selection:  omniVLDB.RamToPar / omniVLDB.DetermineB
//   - sketch + estimator:   omniVLDB.omni.OmniSketch (add / query), Kmin, CountMin
// It lives in package omniVLDB because RecordSynth's constructor is package-private.
//
// Usage:
//   java omniVLDB.Harness run <data.bin> <workload.wl> <out.csv> <name> <repetition> <mem_mb>...
//   java omniVLDB.Harness hashes <repetition> <b> <width> <depth>     (prints test vectors)
//
// Output (one line per query, appended to out.csv):
//   name,mem_mb,repetition,query_index,truth,estimate,case2Estimate,intersectSize,thrm33Case2
// where "estimate" is what OmniSketch.query() returns (Theorem 3.3: the bound when the
// intersection is small) and "case2Estimate" is ceil(|S_cap| * n_max / B) (Corollary 3.4).
package omniVLDB;

import omniVLDB.omni.Kmin;
import omniVLDB.omni.OmniSketch;

import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.BufferedInputStream;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.ArrayList;
import java.util.Random;

public class Harness {
    static int p;
    static long n;
    static int[] data; // row-major, n * p unsigned 32-bit values

    static void loadData(String path) throws Exception {
        try (DataInputStream in = new DataInputStream(new BufferedInputStream(new FileInputStream(path), 1 << 20))) {
            byte[] hdr = new byte[16];
            in.readFully(hdr);
            ByteBuffer bb = ByteBuffer.wrap(hdr).order(ByteOrder.LITTLE_ENDIAN);
            if (hdr[0] != 'S' || hdr[1] != 'T' || hdr[2] != 'R' || hdr[3] != 'M') throw new RuntimeException("bad magic");
            p = bb.getInt(4);
            n = bb.getLong(8);
            long total = n * p;
            if (total > Integer.MAX_VALUE - 8) throw new RuntimeException("dataset too large for one int[]");
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

    static class WQ {
        long truth;
        int[] attrs;
        long[] vals;
    }

    // Text format, one query per line: kind npred truth minfreq a1 v1 a2 v2 ...
    static ArrayList<WQ> loadWorkload(String path) throws Exception {
        ArrayList<WQ> qs = new ArrayList<>();
        try (BufferedReader r = new BufferedReader(new FileReader(path))) {
            String line;
            while ((line = r.readLine()) != null) {
                String[] t = line.trim().split("\\s+");
                if (t.length < 4) continue;
                int np = Integer.parseInt(t[1]);
                WQ q = new WQ();
                q.truth = Long.parseLong(t[2]);
                q.attrs = new int[np];
                q.vals = new long[np];
                for (int i = 0; i < np; i++) {
                    q.attrs[i] = Integer.parseInt(t[4 + 2 * i]);
                    q.vals[i] = Long.parseLong(t[5 + 2 * i]);
                }
                qs.add(q);
            }
        }
        return qs;
    }

    public static void main(String[] args) throws Exception {
        if (args[0].equals("hashes")) {
            hashes(args);
            return;
        }
        if (args[0].equals("gen")) {
            gen(args);
            return;
        }
        String dataPath = args[1], wlPath = args[2], outPath = args[3], name = args[4];
        Main.repetition = Integer.parseInt(args[5]);
        Main.runOnODC = false;
        Main.rangeQueries = false;
        Main.useDS = false;
        Main.sanityBound = false;
        Main.eps = 0.1;
        Main.delta = 0.1;

        loadData(dataPath);
        ArrayList<WQ> wl = loadWorkload(wlPath);
        Main.numAttributes = p;
        Main.streamSize = (int) n;
        System.err.println("loaded " + n + " records, " + p + " attributes, " + wl.size() + " queries");

        boolean[] hasPred = new boolean[p];
        for (int i = 0; i < p; i++) hasPred[i] = true;

        long[] rams = new long[args.length - 6];
        for (int i = 6; i < args.length; i++) rams[i - 6] = (long) (Double.parseDouble(args[i]) * 8e6);
        RamToPar rtp = new RamToPar(p, rams); // upstream parameter selection (eps = 0.1)

        try (PrintWriter out = new PrintWriter(new FileWriter(outPath, true))) {
            for (int ri = 0; ri < rams.length; ri++) {
                long ram = rams[ri];
                double[] params = rtp.getParamsSketch(ram);
                Main.depth = (int) params[0];
                Main.width = (int) params[1];
                Main.maxSize = (int) params[2];
                Main.b = (int) params[3];
                if (Main.maxSize < 1) {
                    System.err.println("budget too small: " + args[6 + ri]);
                    continue;
                }
                OmniSketch s = new OmniSketch(hasPred, ram);
                long t0 = System.nanoTime();
                for (int i = 0; i < n; i++) {
                    long[] rec = new long[p];
                    int base = i * p;
                    for (int a = 0; a < p; a++) rec[a] = data[base + a] & 0xFFFFFFFFL;
                    s.add(new RecordSynth(i, rec));
                }
                double tIngest = (System.nanoTime() - t0) / 1e9;
                long memBits = s.getMemoryUsage();
                long t1 = System.nanoTime();
                int qi = 0;
                StringBuilder sb = new StringBuilder();
                for (WQ w : wl) {
                    Query q = new Query(qi);
                    for (int k = 0; k < w.attrs.length; k++) q.add(w.attrs[k], w.vals[k]);
                    q.exactAnswer = w.truth;
                    int est = s.query(q);
                    sb.append(name).append(',').append(args[6 + ri]).append(',').append(Main.repetition).append(',')
                      .append(qi).append(',').append(w.truth).append(',').append(est).append(',')
                      .append((long) q.case2Estimate).append(',').append((long) q.intersectSize).append(',')
                      .append(q.thrm33Case2 ? 1 : 0).append('\n');
                    qi++;
                }
                double tQuery = (System.nanoTime() - t1) / 1e9;
                out.print(sb);
                out.flush();
                System.out.printf("JRUN,%s,%s,%d,%d,%d,%.3f,%.2f,%.2f,%.3f,w=%d;d=%d;B=%d;b=%d%n", name, args[6 + ri],
                        Main.repetition, p, n, tIngest, tIngest / n * 1e9, tQuery / wl.size() * 1e6, memBits / 8e6,
                        Main.width, Main.depth, Main.maxSize, Main.b);
                s.reset();
            }
        }
    }

    // gen <out.bin> <numRecords> <numAttributes> <alpha>: synthetic stream produced by the
    // upstream generator (omniVLDB.ZipfGenerator: commons-math ZipfDistribution over
    // {1..10000}; value k has probability proportional to k^-alpha, i.e. the frequent values
    // are the small consecutive integers), written in our binary format.
    static void gen(String[] args) throws Exception {
        int nRec = Integer.parseInt(args[2]), nAttr = Integer.parseInt(args[3]);
        double alpha = Double.parseDouble(args[4]);
        long[][] z = ZipfGenerator.zipfData(nRec, nAttr, 10000, alpha);
        try (java.io.FileOutputStream fo = new java.io.FileOutputStream(args[1]);
             java.io.BufferedOutputStream bo = new java.io.BufferedOutputStream(fo, 1 << 20)) {
            ByteBuffer hdr = ByteBuffer.allocate(16).order(ByteOrder.LITTLE_ENDIAN);
            hdr.put((byte) 'S').put((byte) 'T').put((byte) 'R').put((byte) 'M').putInt(nAttr).putLong(nRec);
            bo.write(hdr.array());
            ByteBuffer row = ByteBuffer.allocate(4 * nAttr).order(ByteOrder.LITTLE_ENDIAN);
            for (int i = 0; i < nRec; i++) {
                row.clear();
                for (int a = 0; a < nAttr; a++) row.putInt((int) z[i][a]);
                bo.write(row.array());
            }
        }
    }

    // Test vectors for the bit-exact C++ port of the upstream hash functions:
    //   record hash  g(id)   = Kmin.hash(id)            (java.util.Random seeded with id + repetition)
    //   cell hashes  h_j(v)  = Random(v + repetition).nextInt(width), j = 1..depth
    static void hashes(String[] args) {
        Main.repetition = Integer.parseInt(args[1]);
        Main.b = Integer.parseInt(args[2]);
        Main.width = Integer.parseInt(args[3]);
        Main.depth = Integer.parseInt(args[4]);
        Main.maxSize = 16;
        Kmin k = new Kmin();
        long[] ids = {0, 1, 2, 3, 12345, 4999999, 11000000, 2147483000L};
        for (long id : ids) {
            Random rn = new Random();
            rn.setSeed(id + Main.repetition);
            long direct = rn.nextLong(Math.min((int) Math.pow(2, Main.b), Integer.MAX_VALUE));
            System.out.println("G " + id + " " + direct);
        }
        long[] vals = {0, 1, 2, 9999, 65535, 3232235777L, 4294967295L};
        for (long v : vals) {
            Random rn = new Random();
            rn.setSeed(v + Main.repetition);
            StringBuilder sb = new StringBuilder("H " + v);
            for (int j = 0; j < Main.depth; j++) sb.append(' ').append(rn.nextInt(Main.width));
            System.out.println(sb);
        }
    }
}
