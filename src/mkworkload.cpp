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

// kind 4 "record-sampled range": pick a record r uniformly at random and k = 1..4 predicates.
// The first predicate is a range on a random range attribute a, of width W drawn log-uniformly
// from [2, max_a + 1] and placed uniformly among the positions that contain r_a (clipped to the
// domain); for k >= 2 the second predicate is a range on another range attribute with probability
// 1/2; the remaining predicates are equalities on further random attributes, with r's values.
// Every query is non-empty by construction. Each predicate is written as  a lo hi.
static int range_workload(const Dataset& d, const char* out, int per_k, uint64_t seed, const std::string& alist) {
    std::vector<uint32_t> ra;
    for (size_t s = 0; s < alist.size();) {
        size_t e = alist.find(',', s);
        ra.push_back((uint32_t)atoi(alist.substr(s, e == std::string::npos ? e : e - s).c_str()));
        s = e == std::string::npos ? alist.size() : e + 1;
    }
    if (ra.empty()) { fprintf(stderr, "range: give the range attributes as a1,a2,...\n"); return 1; }
    const uint32_t p = d.p;
    std::vector<uint32_t> mx(p, 0);
    for (uint64_t i = 0; i < d.n; i++)
        for (uint32_t a = 0; a < p; a++) mx[a] = std::max(mx[a], d.row(i)[a]);
    Rng rng(seed);
    std::vector<Query> qs;
    std::unordered_set<std::string> seen;
    auto range_pred = [&](Query& q, uint32_t a, uint32_t x) {
        double W = std::exp(std::log(2.0) + rng.unit() * (std::log((double)mx[a] + 1.0) - std::log(2.0)));
        uint64_t w = std::max<uint64_t>(2, (uint64_t)W);
        uint64_t off = rng.below(w); // position of x inside the range
        int64_t lo = (int64_t)x - (int64_t)off, hi = lo + (int64_t)w - 1;
        if (lo < 0) lo = 0;
        if (hi > (int64_t)mx[a]) hi = mx[a];
        if (lo == hi) { if (hi < (int64_t)mx[a]) hi++; else if (lo > 0) lo--; }
        q.attrs.push_back(a), q.vals.push_back((uint32_t)lo), q.hi.push_back((uint32_t)hi);
    };
    for (uint32_t k = 1; k <= std::min<uint32_t>(4, p); k++) {
        int made = 0, tries = 0;
        while (made < per_k && tries++ < per_k * 50) {
            const uint32_t* r = d.row(rng.below(d.n));
            Query q;
            q.kind = 4;
            uint32_t a0 = ra[rng.below(ra.size())];
            range_pred(q, a0, r[a0]);
            if (k >= 2 && ra.size() > 1 && rng.below(2) == 0) {
                uint32_t a1;
                do a1 = ra[rng.below(ra.size())]; while (a1 == a0);
                range_pred(q, a1, r[a1]);
            }
            while (q.attrs.size() < k) {
                uint32_t a = (uint32_t)rng.below(p);
                if (std::find(q.attrs.begin(), q.attrs.end(), a) != q.attrs.end()) continue;
                q.attrs.push_back(a), q.vals.push_back(r[a]), q.hi.push_back(r[a]);
            }
            std::string key;
            std::vector<size_t> ord(q.attrs.size());
            for (size_t i = 0; i < ord.size(); i++) ord[i] = i;
            std::sort(ord.begin(), ord.end(), [&](size_t x, size_t y) { return q.attrs[x] < q.attrs[y]; });
            for (size_t i : ord) key += std::to_string(q.attrs[i]) + "=" + std::to_string(q.vals[i]) + ".." + std::to_string(q.hi[i]) + ";";
            if (!seen.insert(key).second) continue;
            qs.push_back(std::move(q));
            made++;
        }
    }
    // exact answers and per-predicate frequencies by a scan
    std::vector<std::vector<uint64_t>> pf(qs.size());
    for (size_t i = 0; i < qs.size(); i++) pf[i].assign(qs[i].attrs.size(), 0);
    for (uint64_t i = 0; i < d.n; i++) {
        const uint32_t* r = d.row(i);
        for (size_t j = 0; j < qs.size(); j++) {
            const Query& q = qs[j];
            bool all = true;
            for (size_t t = 0; t < q.attrs.size(); t++) {
                uint32_t x = r[q.attrs[t]];
                bool ok = x >= q.vals[t] && x <= q.hi[t];
                pf[j][t] += ok;
                all &= ok;
            }
            qs[j].truth += all;
        }
    }
    FILE* f = fopen(out, "w");
    for (size_t j = 0; j < qs.size(); j++) {
        const Query& q = qs[j];
        uint64_t mf = *std::min_element(pf[j].begin(), pf[j].end());
        fprintf(f, "%d %zu %llu %llu", q.kind, q.attrs.size(), (unsigned long long)q.truth, (unsigned long long)mf);
        for (size_t t = 0; t < q.attrs.size(); t++) fprintf(f, " %u %u %u", q.attrs[t], q.vals[t], q.hi[t]);
        fprintf(f, "\n");
    }
    fclose(f);
    fprintf(stderr, "wrote %zu range queries to %s\n", qs.size(), out);
    return 0;
}

int main(int argc, char** argv) {
    if (argc < 6) {
        fprintf(stderr, "usage: mkworkload data out queries_per_k kmax seed [prefix_p]\n"
                        "       mkworkload data out pvldb|vldbj num_records seed [prefix_p]\n"
                        "       mkworkload data out range queries_per_k seed a1,a2,...\n");
        return 1;
    }
    Dataset d;
    if (!d.load(argv[1])) {
        fprintf(stderr, "cannot load %s\n", argv[1]);
        return 1;
    }
    std::string family = argv[3];
    if (family == "range") return range_workload(d, argv[2], atoi(argv[4]), strtoull(argv[5], 0, 10), argc > 6 ? argv[6] : "");
    if (argc > 6) d.project_prefix((uint32_t)atoi(argv[6]));
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
    // MKW_TEMPLATES=T: record-sampled queries (kind 0 only) whose attributes are drawn from one of
    // T fixed random attribute sets of size kmax ("query templates"), instead of from all p
    // attributes; models a workload whose attribute combinations are known in advance.
    std::vector<std::vector<uint32_t>> tmpl;
    if (const char* e = getenv("MKW_TEMPLATES")) {
        Rng tr(strtoull(argv[5], 0, 10) ^ 0x7E3D1B5A9C2F4E61ULL);
        for (int t = 0; t < atoi(e); t++) {
            std::vector<uint32_t> s;
            while (s.size() < kmax) {
                uint32_t a = (uint32_t)tr.below(p);
                if (std::find(s.begin(), s.end(), a) == s.end()) s.push_back(a);
            }
            tmpl.push_back(s);
        }
        eligible.clear(); // no anchored family
    }
    auto add_query = [&](int kind, uint64_t rec, uint32_t k, int anchor) {
        std::vector<uint32_t> attrs;
        if (anchor >= 0) attrs.push_back((uint32_t)anchor);
        const std::vector<uint32_t>* from = tmpl.empty() ? nullptr : &tmpl[rng.below(tmpl.size())];
        while (attrs.size() < k) {
            uint32_t a = from ? (*from)[rng.below(from->size())] : (uint32_t)rng.below(p);
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
