// Update stream with bounded deletions, built as in the journal version's data pipeline
// (github.com/wiegerpunter/dataset-pipeline: splitDeletesSynth.py, MixInMemorysynthFromDiskDB.java).
//   mkdelstream <base.bin> <n_res> <perc> <seed> <out.ops> [<res.bin>]
// Rows [0, n_res) of base.bin are the residue (the records that stay); the next floor(perc * n_res)
// rows are the noise records, each inserted and deleted again. The list of residue ids plus every
// noise id twice is shuffled; the first occurrence of a noise id is its insertion and the second its
// deletion, so insertions and deletions are spread over the whole stream and every deletion follows
// its insertion (their spreadOutDeletes setting). Output: "UPDS" | u64 m | m x int32, +(i+1) inserts
// row i, -(i+1) deletes it. res.bin (optional) receives the residue as a dataset, on which the query
// workloads and their true answers are computed.
#include "common.hpp"

int main(int argc, char** argv) {
    if (argc < 6) {
        fprintf(stderr, "usage: mkdelstream <base.bin> <n_res> <perc> <seed> <out.ops> [<res.bin>]\n");
        return 1;
    }
    Dataset d;
    if (!d.load(argv[1])) { fprintf(stderr, "cannot load %s\n", argv[1]); return 1; }
    uint64_t nres = strtoull(argv[2], 0, 10);
    double perc = atof(argv[3]);
    uint64_t seed = strtoull(argv[4], 0, 10);
    uint64_t nnoise = (uint64_t)(perc * (double)nres);
    if (nres + nnoise > d.n || nres + nnoise > (uint64_t)INT32_MAX) {
        fprintf(stderr, "base has %llu rows, need %llu\n", (unsigned long long)d.n, (unsigned long long)(nres + nnoise));
        return 1;
    }
    std::vector<int32_t> ids;
    ids.reserve(nres + 2 * nnoise);
    for (uint64_t i = 0; i < nres; i++) ids.push_back((int32_t)i);
    for (int rep = 0; rep < 2; rep++)
        for (uint64_t i = nres; i < nres + nnoise; i++) ids.push_back((int32_t)i);
    Rng r(seed);
    for (uint64_t i = ids.size(); i > 1; i--) std::swap(ids[i - 1], ids[r.below(i)]);
    std::vector<uint8_t> seen(nnoise, 0);
    uint64_t live = 0, peak = 0;
    for (auto& x : ids) {
        if ((uint64_t)x < nres) {
            x = x + 1;
            live++;
        } else if (!seen[x - nres]) {
            seen[x - nres] = 1;
            x = x + 1;
            live++;
        } else {
            x = -(x + 1);
            live--;
        }
        peak = std::max(peak, live);
    }
    FILE* f = fopen(argv[5], "wb");
    uint64_t m = ids.size();
    fwrite("UPDS", 1, 4, f);
    fwrite(&m, 8, 1, f);
    fwrite(ids.data(), 4, m, f);
    fclose(f);
    if (argc > 6) {
        FILE* g = fopen(argv[6], "wb");
        fwrite("STRM", 1, 4, g);
        fwrite(&d.p, 4, 1, g);
        fwrite(&nres, 8, 1, g);
        fwrite(d.v.data(), 4, nres * d.p, g);
        fclose(g);
    }
    printf("updates=%llu inserts=%llu deletes=%llu residue=%llu peak_live=%llu\n", (unsigned long long)m,
           (unsigned long long)(nres + nnoise), (unsigned long long)nnoise, (unsigned long long)nres, (unsigned long long)peak);
    return 0;
}
