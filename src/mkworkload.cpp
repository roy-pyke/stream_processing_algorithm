// Workload generator with exact ground truth.
//   mkworkload <data.bin> <out.txt> <queries_per_k> <kmax> <seed>
// Two families of conjunctive equality queries (all non-empty by construction):
//   kind 0  "record-sampled": pick a record uniformly at random, pick k random
//           attributes, use the record's values. This is the workload of the
//           OmniSketch paper; it favours frequent values.
//   kind 1  "frequency-stratified anchor": pick an attribute, pick a frequency
//           decade [10^j, 10^(j+1)) uniformly among the non-empty ones (j >= 1),
//           pick a value of that attribute uniformly inside the decade, pick a
//           random record carrying it, and add k-1 further predicates taken
//           from that record. This covers rare and frequent anchors evenly
//           (the "drill-down on one entity" query shape).
#include "common.hpp"
#include <unordered_map>
#include <unordered_set>

struct VInfo {
    uint32_t count = 0;
    uint32_t rec = 0;
};

int main(int argc, char** argv) {
    if (argc < 6) {
        fprintf(stderr, "usage: mkworkload data out queries_per_k kmax seed [prefix_p]\n"
                        "       mkworkload data out pvldb|vldbj num_records seed [prefix_p]\n");
        return 1;
    }
    Dataset d;
    if (!d.load(argv[1])) {
        fprintf(stderr, "cannot load %s\n", argv[1]);
        return 1;
    }
    if (argc > 6) d.project_prefix((uint32_t)atoi(argv[6]));
    std::string family = argv[3];
    int per_k = atoi(argv[3]);
    uint32_t kmax = std::min<uint32_t>((uint32_t)atoi(argv[4]), d.p);
    Rng rng(strtoull(argv[5], 0, 10));
    const uint32_t p = d.p;

    // per-attribute value statistics with a size-1 reservoir of carrying records
    std::vector<std::unordered_map<uint32_t, VInfo>> stats(p);
    for (uint64_t i = 0; i < d.n; i++) {
        const uint32_t* r = d.row(i);
        for (uint32_t a = 0; a < p; a++) {
            VInfo& vi = stats[a][r[a]];
            vi.count++;
            if (rng.below(vi.count) == 0) vi.rec = (uint32_t)i;
        }
    }
    // decade lists
    std::vector<std::vector<std::vector<uint32_t>>> dec(p); // dec[a][j] = values with count in [10^j,10^(j+1))
    std::vector<uint32_t> eligible;
    for (uint32_t a = 0; a < p; a++) {
        dec[a].resize(10);
        for (auto& kv : stats[a]) {
            int j = (int)std::floor(std::log10((double)kv.second.count));
            dec[a][j].push_back(kv.first);
        }
        for (auto& l : dec[a]) std::sort(l.begin(), l.end());
        fprintf(stderr, "attr %u: %zu distinct\n", a, stats[a].size());
        if (stats[a].size() >= 20) eligible.push_back(a);
    }

    std::vector<Query> qs;
    std::unordered_set<std::string> seen;
    auto add_query = [&](int kind, uint64_t rec, uint32_t k, int anchor) {
        std::vector<uint32_t> attrs;
        if (anchor >= 0) attrs.push_back((uint32_t)anchor);
        while (attrs.size() < k) {
            uint32_t a = (uint32_t)rng.below(p);
            if (std::find(attrs.begin(), attrs.end(), a) == attrs.end()) attrs.push_back(a);
        }
        std::vector<uint32_t> sorted = attrs;
        std::sort(sorted.begin(), sorted.end());
        std::string key;
        const uint32_t* r = d.row(rec);
        for (uint32_t a : sorted) key += std::to_string(a) + "=" + std::to_string(r[a]) + ";";
        if (!seen.insert(std::to_string(kind) + "|" + key).second) return false;
        Query q;
        q.kind = kind;
        for (uint32_t a : attrs) q.attrs.push_back(a), q.vals.push_back(r[a]);
        qs.push_back(std::move(q));
        return true;
    };
    // Workload procedures of the OmniSketch papers (selected with "pvldb" / "vldbj" in place of
    // queries_per_k; the next argument is then the number of sampled records):
    //   kind 2 "pvldb" (Dataset.synthWorkload in the authors' omnisketch repository): for every
    //           sampled record, one query per subset of at least two attributes.
    //   kind 3 "vldbj" (SyntheticDataset.generateQueries in their omni-deletes repository): for
    //           every sampled record and every k = 1..p, one query on k random attributes.
    auto add_fixed = [&](int kind, uint64_t rec, std::vector<uint32_t> attrs) {
        std::sort(attrs.begin(), attrs.end());
        std::string key;
        const uint32_t* r = d.row(rec);
        for (uint32_t a : attrs) key += std::to_string(a) + "=" + std::to_string(r[a]) + ";";
        if (!seen.insert(std::to_string(kind) + "|" + key).second) return false;
        Query q;
        q.kind = kind;
        for (uint32_t a : attrs) q.attrs.push_back(a), q.vals.push_back(r[a]);
        qs.push_back(std::move(q));
        return true;
    };
    if (family == "pvldb" || family == "vldbj") {
        int nrec = atoi(argv[4]);
        for (int i = 0; i < nrec; i++) {
            uint64_t rec = rng.below(d.n);
            if (family == "pvldb") {
                if (p > 16) { fprintf(stderr, "pvldb workload needs p <= 16\n"); return 1; }
                for (uint32_t mask = 1; mask < (1u << p); mask++) {
                    if (__builtin_popcount(mask) < 2) continue;
                    std::vector<uint32_t> attrs;
                    for (uint32_t a = 0; a < p; a++)
                        if (mask >> a & 1) attrs.push_back(a);
                    add_fixed(2, rec, attrs);
                }
            } else {
                for (uint32_t k = 1; k <= p; k++) {
                    std::vector<uint32_t> attrs;
                    while (attrs.size() < k) {
                        uint32_t a = (uint32_t)rng.below(p);
                        if (std::find(attrs.begin(), attrs.end(), a) == attrs.end()) attrs.push_back(a);
                    }
                    add_fixed(3, rec, attrs);
                }
            }
        }
        kmax = 0;
    }
    for (uint32_t k = 1; k <= kmax; k++) {
        int made = 0, tries = 0;
        while (made < per_k && tries++ < per_k * 50) made += add_query(0, rng.below(d.n), k, -1);
        if (eligible.empty()) continue;
        made = 0, tries = 0;
        while (made < per_k && tries++ < per_k * 50) {
            uint32_t a = eligible[rng.below(eligible.size())];
            std::vector<int> js;
            for (int j = 1; j < 10; j++)
                if (!dec[a][j].empty()) js.push_back(j);
            if (js.empty()) continue;
            int j = js[rng.below(js.size())];
            uint32_t v = dec[a][j][rng.below(dec[a][j].size())];
            made += add_query(1, stats[a][v].rec, k, (int)a);
        }
    }
    // exact ground truth: index queries by their first predicate
    std::unordered_map<uint64_t, std::vector<uint32_t>> index;
    std::vector<char> used(p, 0);
    for (uint32_t i = 0; i < qs.size(); i++) {
        index[((uint64_t)qs[i].attrs[0] << 32) | qs[i].vals[0]].push_back(i);
        used[qs[i].attrs[0]] = 1;
    }
    for (uint64_t i = 0; i < d.n; i++) {
        const uint32_t* r = d.row(i);
        for (uint32_t a = 0; a < p; a++) {
            if (!used[a]) continue;
            auto it = index.find(((uint64_t)a << 32) | r[a]);
            if (it == index.end()) continue;
            for (uint32_t qi : it->second)
                if (qs[qi].matches(r)) qs[qi].truth++;
        }
    }
    FILE* f = fopen(argv[2], "w");
    for (auto& q : qs) {
        uint64_t mf = UINT64_MAX;
        for (size_t i = 0; i < q.attrs.size(); i++) mf = std::min<uint64_t>(mf, stats[q.attrs[i]][q.vals[i]].count);
        fprintf(f, "%d %zu %llu %llu", q.kind, q.attrs.size(), (unsigned long long)q.truth, (unsigned long long)mf);
        for (size_t i = 0; i < q.attrs.size(); i++) fprintf(f, " %u %u", q.attrs[i], q.vals[i]);
        fprintf(f, "\n");
    }
    fclose(f);
    // number of distinct full tuples (relevant for tuple-aggregating synopses)
    std::unordered_set<uint64_t> tup;
    for (uint64_t i = 0; i < d.n; i++) {
        const uint32_t* r = d.row(i);
        uint64_t h = 1469598103934665603ULL;
        for (uint32_t a = 0; a < p; a++) h = mix64(h ^ r[a]);
        tup.insert(h);
    }
    fprintf(stderr, "wrote %zu queries to %s; n=%llu p=%u distinct_tuples=%zu\n", qs.size(), argv[2],
            (unsigned long long)d.n, p, tup.size());
    return 0;
}
