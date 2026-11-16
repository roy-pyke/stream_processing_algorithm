// Benchmark driver.
//   bench acc   <data.bin> <workload.txt> <out.csv> <name> <mem_mb> <seed> [-p k] <spec>...
//   bench time  <data.bin> <name> <mem_mb> <seed> <reps> [-p k] [-w workload.txt] <spec>...
//               (-w supplies range attributes for dy=1 / omnid and the workload for psel=wl)
//   bench merge <data.bin> <workload.txt> <out.csv> <name> <mem_mb> <seed> <parts> [-p k] <spec>...
// Method specs (all sized to the same memory budget of mem_mb * 8e6 bits):
//   res            uniform row reservoir
//   uss            Unbiased Space-Saving on full tuples
//   coco           CocoSketch (Zhang et al., SIGCOMM 2021) on full tuples, d = 2, 32p + 32 bits/bucket
//   omnip          OmniSketch, paper parameters (eps=delta=0.1) and paper memory accounting
//   omnil          OmniSketch, paper (w,d), lean accounting (32 bits/sample, dense)
//   omni:w:d       OmniSketch with given (w,d); B chosen by an oracle so that the
//                  *actually stored* samples fill the budget (most favourable accounting)
//   omnid:w:d      OmniSketch with dyadic range support (Sec. 3.3 of the paper) on the attributes
//                  that carry range predicates in the workload; oracle B, 32 bits per stored hash
//                  and per non-empty cell; "est" holds the reference code's estimator, "ht" the
//                  threshold estimator (see OmniDyad)
//   omniv          OmniSketch as configured in the journal version's code (w=28, d=3, 31 bits/hash)
//   ash / ashl     the journal version's sampling baseline (adaptive Sample-and-Hold on full
//                  tuples) with its accounting (32p+384 bits/entry) / lean accounting (32p+64)
//   resj           reservoir sized as in the journal version's code (64p+32 bits per row)
//   omnijs         as omnij, but "est" holds the always-scaled estimate and "ht" the Theorem 3.3 one
//   omnij          bit-exact port of the authors' Java implementation (hashes, parameter search,
//                  accounting, Theorem 3.3 estimator); the seed is used as their "repetition"
//   cs:w:lf:gamma[:phi2[:phi3[:rake[:index[:cap]]]]]  CoStrata with w buckets/attribute, house
//                  weight lf, exponent gamma; phi2 / phi3 > 0 add pair / triple counter tables
//                  using that fraction of the budget; rake = 1 calibrates to all anchors of a query
//                  (IPF) instead of the lightest one; index = 1 builds a query-time bucket index
//                  (not charged); cap > 0 limits the average number of table increments per record
//                  by updating the tables for a hashed fraction of the records only (sampled tables)
//   csh:w:lf[:phi2[:phi3]]  CoStrata-H: CoStrata over distinct tuples with sample-and-hold counts
//   ussc / cococ:w[:phi2[:phi3]]  USS / CocoSketch with C1 (and C2 / C3) calibration
// Options, appended to a spec as /key=value (e.g. cs:0:1:1:0.25/psel=mi/np=136):
//   ix=2           (cs) build the query-time index and charge it: 32p bits more per row and
//                  32p(w+1) bits of bucket offsets come out of the budget
//   dy=1           (cs) dyadic counters of the range attributes as further anchors (charged)
//   psel=wl|mi|sc|rnd, np=m   (cs, csh, ussc, cococ) keep m pair tables only, chosen by: the
//                  number of workload queries using the pair (the workload is known), the mutual
//                  information on a stream prefix, the anchor reduction on a prefix (see
//                  select_pairs), or at random; the pair tables then use their budget fraction
//                  with m tables instead of C(p,2)
// Output columns of "acc": name,spec,mem,seed,query,est,ht,sd,cm[,m,wmax]. For OmniSketch specs
// the last three hold: the Theorem 3.3 estimate (bound when the intersection is small), |S_cap|,
// n_max; for omnij "est" is the Theorem 3.3 estimate and "ht" the always-scaled one. CoStrata and
// CoStrata-H add the number of matching rows (held count) and the largest calibrated weight a
// matching record can have, used for the small-count confidence intervals.
// "acc" appends one line per (method, query) to out.csv and prints a RUN line. "merge" splits the
// stream into <parts> contiguous pieces, builds one synopsis per piece (same seed, the full budget
// each), merges them in order and answers the workload from the merged synopsis (cs, omni:w:d).
#include "synopses.hpp"
#include "omni_java_compat.hpp"
#include <map>
#include <memory>
#include <unordered_map>

struct Spec {
    std::string raw, kind;
    std::vector<double> a;
    std::map<std::string, std::string> opt;
    std::string get(const std::string& k, const std::string& def = "") const {
        auto it = opt.find(k);
        return it == opt.end() ? def : it->second;
    }
    double num(const std::string& k, double def = 0) const {
        auto it = opt.find(k);
        return it == opt.end() ? def : atof(it->second.c_str());
    }
};
static Spec parse_spec(const std::string& s) {
    Spec sp;
    sp.raw = s;
    std::string body = s.substr(0, s.find('/'));
    for (size_t q = s.find('/'); q != std::string::npos;) {
        size_t nx = s.find('/', q + 1);
        std::string kv = s.substr(q + 1, nx == std::string::npos ? nx : nx - q - 1);
        size_t eq = kv.find('=');
        sp.opt[kv.substr(0, eq)] = eq == std::string::npos ? "1" : kv.substr(eq + 1);
        q = nx;
    }
    size_t pos = body.find(':');
    sp.kind = body.substr(0, pos);
    while (pos != std::string::npos) {
        size_t nx = body.find(':', pos + 1);
        sp.a.push_back(atof(body.substr(pos + 1, nx == std::string::npos ? nx : nx - pos - 1).c_str()));
        pos = nx;
    }
    return sp;
}

template <class S>
static double ingest(S& syn, const Dataset& d, uint64_t from = 0, uint64_t to = UINT64_MAX) {
    double t0 = now_sec();
    to = std::min(to, d.n);
    for (uint64_t i = from; i < to; i++) syn.insert(i, d.row(i));
    return now_sec() - t0;
}

// w = 0 selects the default: the largest power of two such that the counters use <= 5% of the budget.
static uint32_t auto_w(double wspec, uint64_t bits, uint32_t p) {
    if (wspec > 0) return (uint32_t)wspec;
    uint32_t w = 64;
    while (w < 65536 && 32.0 * p * (2.0 * w) <= 0.05 * (double)bits) w *= 2;
    return w;
}

// ---------------------------------------------------------------- pair selection (wide streams)
// Ranks the C(p,2) attribute pairs and returns the m best:
//   wl   number of workload queries whose predicates include the pair (workload known in advance)
//   mi   mutual information of the two attributes on a stream prefix
//   sc   anchor reduction on a prefix: the mean over prefix records r of 1 - n_ab / min(n_a, n_b),
//        the fraction of the lighter C1 cell that the pair cell removes for the query (r_a, r_b)
//   rnd  uniformly random
// Ties are broken at random. Prefix statistics are exact counts over the first `prefix` records.
static std::vector<std::pair<uint32_t, uint32_t>> select_pairs(const std::string& mode, uint64_t m, const Dataset& d,
                                                               const std::vector<Query>& qs, uint64_t seed,
                                                               uint64_t prefix = 50000) {
    const uint32_t p = d.p;
    std::vector<std::pair<uint32_t, uint32_t>> all;
    for (uint32_t a = 0; a < p; a++)
        for (uint32_t b = a + 1; b < p; b++) all.push_back({a, b});
    Rng rng(seed ^ 0xA5A5A5A55A5A5A5AULL);
    std::vector<double> score(all.size(), 0.0), tie(all.size());
    for (auto& t : tie) t = rng.unit();
    if (mode == "wl") {
        std::vector<double> c((size_t)p * p, 0.0);
        for (auto& q : qs)
            for (size_t i = 0; i < q.attrs.size(); i++)
                for (size_t j = i + 1; j < q.attrs.size(); j++) {
                    uint32_t a = std::min(q.attrs[i], q.attrs[j]), b = std::max(q.attrs[i], q.attrs[j]);
                    c[(size_t)a * p + b] += 1;
                }
        for (size_t t = 0; t < all.size(); t++) score[t] = c[(size_t)all[t].first * p + all[t].second];
    } else if (mode == "mi" || mode == "sc") {
        uint64_t n0 = std::min<uint64_t>(prefix, d.n);
        std::vector<std::unordered_map<uint32_t, uint32_t>> one(p);
        for (uint64_t i = 0; i < n0; i++)
            for (uint32_t a = 0; a < p; a++) one[a][d.row(i)[a]]++;
        std::unordered_map<uint64_t, uint32_t> two;
        for (size_t t = 0; t < all.size(); t++) {
            uint32_t a = all[t].first, b = all[t].second;
            two.clear();
            for (uint64_t i = 0; i < n0; i++) two[((uint64_t)d.row(i)[a] << 32) | d.row(i)[b]]++;
            double s = 0;
            for (auto& kv : two) {
                double nab = kv.second, na = one[a][(uint32_t)(kv.first >> 32)], nb = one[b][(uint32_t)kv.first];
                if (mode == "mi") s += nab / (double)n0 * std::log(nab * (double)n0 / (na * nb));
                else s += nab / (double)n0 * (1.0 - nab / std::min(na, nb));
            }
            score[t] = s;
        }
    } // rnd: all scores 0, order decided by the random tie-breaker
    std::vector<size_t> ord(all.size());
    for (size_t t = 0; t < ord.size(); t++) ord[t] = t;
    std::sort(ord.begin(), ord.end(), [&](size_t x, size_t y) { return score[x] != score[y] ? score[x] > score[y] : tie[x] < tie[y]; });
    std::vector<std::pair<uint32_t, uint32_t>> sel;
    for (size_t t = 0; t < ord.size() && sel.size() < m; t++) sel.push_back(all[ord[t]]);
    return sel;
}

// Range attributes: those carrying a (non-degenerate) range predicate in the workload. levels[a] =
// number of dyadic levels (bits of the largest value + 1), 1 for all other attributes.
struct Ranges {
    std::vector<uint32_t> levels, maxv;
    bool any = false;
    void init(const Dataset& d, const std::vector<Query>& qs) {
        levels.assign(d.p, 1);
        maxv.assign(d.p, 0);
        std::vector<char> is((size_t)d.p, 0);
        for (auto& q : qs)
            for (size_t i = 0; i < q.attrs.size(); i++)
                if (!q.is_eq(i) && q.attrs[i] < d.p) is[q.attrs[i]] = 1;
        for (uint64_t i = 0; i < d.n; i++)
            for (uint32_t a = 0; a < d.p; a++)
                if (is[a]) maxv[a] = std::max(maxv[a], d.row(i)[a]);
        for (uint32_t a = 0; a < d.p; a++)
            if (is[a]) {
                uint32_t bits = 1;
                while (bits < 32 && (maxv[a] >> bits)) bits++;
                levels[a] = bits + 1;
                any = true;
            }
    }
};

struct Ctx { // what spec construction may need besides the data
    const std::vector<Query>* qs = nullptr;
    const Ranges* rg = nullptr;
};

// Widths of the pair / triple calibration tables from the budget fractions given as the
// spec arguments number `at` (pairs) and `at + 1` (triples); 0 disables a level. With the
// option np=m only m pair tables are kept (chosen by psel) and share the pair budget.
static void marg_widths(const Spec& sp, size_t at, uint64_t bits, uint32_t p, uint32_t& w2, uint32_t& w3, uint64_t& npair) {
    double phi2 = sp.a.size() > at ? sp.a[at] : 0.0, phi3 = sp.a.size() > at + 1 ? sp.a[at + 1] : 0.0;
    uint64_t np = (uint64_t)sp.num("np", 0);
    npair = np > 0 ? std::min<uint64_t>(np, MargCal::npairs(p)) : 0;
    w2 = phi2 > 0 ? MargCal::width(phi2 * (double)bits, npair ? npair : MargCal::npairs(p)) : 0;
    w3 = phi3 > 0 && p <= 64 ? MargCal::width(phi3 * (double)bits, MargCal::ntriples(p)) : 0;
}
static std::vector<std::pair<uint32_t, uint32_t>> pairs_for(const Spec& sp, uint64_t npair, const Dataset& d, const Ctx& cx, uint64_t seed) {
    if (!npair) return {};
    static const std::vector<Query> none;
    return select_pairs(sp.get("psel", "rnd"), npair, d, cx.qs ? *cx.qs : none, seed);
}

static OmniSketch* make_omni(const Spec& sp, const Dataset& d, uint64_t bits, uint64_t seed, std::string& info) {
    uint32_t w, dd;
    uint64_t B = 1;
    if (sp.kind == "omnip" || sp.kind == "omnil") {
        OmniSketch::paper_wd(0.1, 0.1, w, dd);
        if (sp.kind == "omnip") {
            // exactly the B that the authors' code selects (RamToPar/DetermineB)
            int32_t jb = OmniJava::searchB((int64_t)bits, (int32_t)dd, (int32_t)w, d.p, 0.1);
            B = jb < 1 ? 1 : (uint64_t)jb;
        } else {
            uint64_t lo = 1, hi = 1ull << 34;
            while (lo < hi) {
                uint64_t mid = (lo + hi + 1) / 2;
                uint64_t need = (uint64_t)d.p * dd * w * (32 * mid + 32);
                if (need <= bits) lo = mid; else hi = mid - 1;
            }
            B = lo;
        }
    } else if (sp.kind == "omniv") {
        // configuration of the journal version's code: d = 3, w = 28, 31 bits per stored hash
        // (Formulas.ramOmniKmin: d * w * p * (31 B + 32) bits)
        w = 28;
        dd = 3;
        double per_cell = (double)bits / ((double)d.p * dd * w);
        B = per_cell > 63 ? (uint64_t)((per_cell - 32) / 31) : 1;
    } else {
        w = (uint32_t)sp.a[0];
        dd = (uint32_t)sp.a[1];
        OmniSketch probe(d.p, dd, w, 1, seed);
        std::vector<uint32_t> c((size_t)d.p * dd * w, 0);
        for (uint64_t i = 0; i < d.n; i++) {
            const uint32_t* r = d.row(i);
            for (uint32_t a = 0; a < d.p; a++)
                for (uint32_t j = 0; j < dd; j++) c[probe.cell(a, j, r[a])]++;
        }
        uint64_t fixed = 32ull * c.size();
        uint64_t budget = bits > fixed ? (bits - fixed) / 32 : 0;
        uint64_t lo = 1, hi = d.n;
        while (lo < hi) {
            uint64_t mid = (lo + hi + 1) / 2, tot = 0;
            for (uint32_t x : c) tot += std::min<uint64_t>(x, mid);
            if (tot <= budget) lo = mid; else hi = mid - 1;
        }
        B = lo;
    }
    char buf[128];
    snprintf(buf, sizeof buf, "w=%u;d=%u;B=%llu", w, dd, (unsigned long long)B);
    info = buf;
    return new OmniSketch(d.p, dd, w, B, seed);
}

static OmniDyad* make_omnid(const Spec& sp, const Dataset& d, const Ranges& rg, uint64_t bits, uint64_t seed, std::string& info) {
    uint32_t w = (uint32_t)sp.a[0], dd = (uint32_t)sp.a[1];
    OmniDyad probe(d.p, dd, w, 1, rg.levels, seed);
    std::vector<uint32_t> c = probe.loads(d);
    uint64_t ne = 0;
    for (uint32_t x : c) ne += x > 0;
    uint64_t budget = bits > 32 * ne ? (bits - 32 * ne) / 32 : 0;
    uint64_t lo = 1, hi = d.n;
    while (lo < hi) {
        uint64_t mid = (lo + hi + 1) / 2, tot = 0;
        for (uint32_t x : c) tot += std::min<uint64_t>(x, mid);
        if (tot <= budget) lo = mid; else hi = mid - 1;
    }
    char buf[128];
    snprintf(buf, sizeof buf, "w=%u;d=%u;B=%llu;cells=%zu", w, dd, (unsigned long long)lo, c.size());
    info = buf;
    return new OmniDyad(d.p, dd, w, lo, rg.levels, seed);
}

// CoStrata from its spec (see the header); `info` describes the configuration.
static CoStrata* make_cs(const Spec& sp, const Dataset& d, uint64_t bits, uint64_t seed, const Ctx& cx, std::string& info) {
    uint32_t w = auto_w(sp.a.empty() ? 0 : sp.a[0], bits, d.p);
    double lf = sp.a.size() > 1 ? sp.a[1] : 1.0, gamma = sp.a.size() > 2 ? sp.a[2] : 1.0;
    uint32_t w2, w3;
    uint64_t npair;
    marg_widths(sp, 3, bits, d.p, w2, w3, npair);
    uint64_t cal = MargCal::bits_for(d.p, w2, w3, npair), extra = 0;
    bool dy = sp.num("dy") > 0 && cx.rg && cx.rg->any;
    if (dy) cal += DyadCnt::bits_for(cx.rg->levels, cx.rg->maxv, w);
    if (sp.num("ix") >= 2) cal += 32ull * d.p * (w + 1), extra = 32ull * d.p; // charged index
    uint64_t M = CoStrata::rows_for_budget(bits, d.p, w, cal, extra);
    if (M < 16) return nullptr;
    CoStrata* s = new CoStrata(d.p, w, M, lf, gamma, seed);
    double cap = sp.a.size() > 7 ? sp.a[7] : 0.0; // max. table increments per record (0: no limit)
    double beta = MargCal::beta_for_cap(d.p, w2, w3, cap);
    auto sel = pairs_for(sp, npair, d, cx, seed);
    if (w2 || w3) s->enable_marginals(w2, w3, beta, sel.empty() ? nullptr : &sel);
    if (dy) s->enable_ranges(cx.rg->levels, cx.rg->maxv);
    char buf[200];
    snprintf(buf, sizeof buf, "w=%u;w2=%u;w3=%u;M=%llu;npair=%u;beta=%.4g;dy_mb=%.3f", w, w2, w3, (unsigned long long)M, s->mc.n2, beta,
             (double)s->dy.mem_bits() / 8e6);
    info = buf;
    return s;
}

template <class TS>
static TupleCal<TS>* make_tcal(const Spec& sp, const Dataset& d, uint64_t bits, uint64_t seed, const Ctx& cx, std::string& info) {
    uint32_t w = auto_w(sp.a.empty() ? 0 : sp.a[0], bits, d.p);
    uint32_t w2, w3;
    uint64_t npair;
    marg_widths(sp, 1, bits, d.p, w2, w3, npair);
    auto* s = new TupleCal<TS>(d.p, w, std::max<uint64_t>(1, TupleCal<TS>::rows_for_budget(bits, d.p, w, MargCal::bits_for(d.p, w2, w3, npair))), seed);
    auto sel = pairs_for(sp, npair, d, cx, seed);
    if (w2 || w3) s->enable_marginals(w2, w3, 1.0, sel.empty() ? nullptr : &sel);
    info = "m=" + std::to_string(s->uss.entries()) + ";w=" + std::to_string(w) + ";w2=" + std::to_string(w2) + ";w3=" + std::to_string(w3) +
           ";npair=" + std::to_string(s->mc.n2);
    return s;
}

static CoStrataH* make_csh(const Spec& sp, const Dataset& d, uint64_t bits, uint64_t seed, const Ctx& cx, std::string& info) {
    uint32_t w = auto_w(sp.a.empty() ? 0 : sp.a[0], bits, d.p);
    double lf = sp.a.size() > 1 ? sp.a[1] : 1.0;
    uint32_t w2, w3;
    uint64_t npair;
    marg_widths(sp, 2, bits, d.p, w2, w3, npair);
    uint64_t M = CoStrataH::rows_for_budget(bits, d.p, w, MargCal::bits_for(d.p, w2, w3, npair));
    if (M < 16) return nullptr;
    CoStrataH* s = new CoStrataH(d.p, w, M, lf, seed);
    auto sel = pairs_for(sp, npair, d, cx, seed);
    if (w2 || w3) s->mc.init(d.p, w2, w3, s->useed ^ 0xC2B2AE3D27D4EB4FULL, 1.0, sel.empty() ? nullptr : &sel);
    char buf[64];
    snprintf(buf, sizeof buf, "w=%u;M=%llu", w, (unsigned long long)M);
    info = buf;
    return s;
}

struct Out { // per-query outputs of one method
    std::vector<double> est, ht, sd, cm, m, wmax;
    void init(size_t n, bool full, bool diag) {
        est.assign(n, 0);
        if (full) ht.assign(n, 0), sd.assign(n, 0), cm.assign(n, 0);
        if (diag) m.assign(n, 0), wmax.assign(n, 0);
    }
};

int main(int argc, char** argv) {
    if (argc < 2) return 1;
    std::string mode = argv[1];
    Dataset d;
    if (mode == "acc" || mode == "merge") {
        bool mg = mode == "merge";
        if (argc < (mg ? 10 : 9)) { fprintf(stderr, "bad args\n"); return 1; }
        if (!d.load(argv[2])) { fprintf(stderr, "cannot load data\n"); return 1; }
        std::vector<Query> qs;
        if (!load_workload(argv[3], qs)) { fprintf(stderr, "cannot load workload\n"); return 1; }
        FILE* out = fopen(argv[4], "a");
        std::string name = argv[5];
        double mem_mb = atof(argv[6]);
        uint64_t seed = strtoull(argv[7], 0, 10);
        uint64_t bits = (uint64_t)(mem_mb * 8e6);
        int ai = 8;
        uint32_t parts = 1;
        if (mg) parts = (uint32_t)atoi(argv[ai++]);
        if (std::string(argv[ai]) == "-p") { d.project_prefix((uint32_t)atoi(argv[ai + 1])); ai += 2; }
        Ranges rg;
        rg.init(d, qs);
        Ctx cx;
        cx.qs = &qs;
        cx.rg = &rg;
        for (; ai < argc; ai++) {
            Spec sp = parse_spec(argv[ai]);
            std::string info;
            double t_ing = 0, t_q = 0;
            uint64_t used_bits = 0;
            Out o;
            o.init(qs.size(), false, false);
            if (mg) {
                // ---------------------------------------------------- merge experiment
                auto piece = [&](uint32_t j, uint64_t& from, uint64_t& to) {
                    from = d.n * j / parts;
                    to = d.n * (j + 1) / parts;
                };
                uint64_t from, to;
                if (sp.kind == "cs") {
                    CoStrata* acc = make_cs(sp, d, bits, seed * 7919 + 4, cx, info);
                    if (!acc) { fprintf(stderr, "budget too small for %s\n", sp.raw.c_str()); continue; }
                    piece(0, from, to);
                    t_ing = ingest(*acc, d, from, to);
                    for (uint32_t j = 1; j < parts; j++) {
                        std::string inf2;
                        CoStrata* s = make_cs(sp, d, bits, seed * 7919 + 4, cx, inf2);
                        piece(j, from, to);
                        t_ing += ingest(*s, d, from, to);
                        acc->merge(*s);
                        delete s;
                    }
                    acc->finalize();
                    o.init(qs.size(), true, true);
                    double t0 = now_sec();
                    for (size_t i = 0; i < qs.size(); i++) {
                        CoStrata::Est e = acc->estimate(qs[i]);
                        o.est[i] = e.hajek; o.ht[i] = e.ht; o.sd[i] = e.sd_hajek; o.cm[i] = e.cm; o.m[i] = e.m; o.wmax[i] = e.wmax;
                    }
                    t_q = now_sec() - t0;
                    used_bits = acc->mem_bits();
                    info += ";held=" + std::to_string(acc->size) + ";compactions=" + std::to_string(acc->compactions);
                    delete acc;
                } else if (sp.kind.rfind("omni", 0) == 0 && sp.kind != "omnid") {
                    OmniSketch* acc = make_omni(sp, d, bits, seed * 7919 + 3, info);
                    piece(0, from, to);
                    t_ing = ingest(*acc, d, from, to);
                    for (uint32_t j = 1; j < parts; j++) {
                        OmniSketch s(acc->p, acc->d, acc->w, acc->B, seed * 7919 + 3);
                        piece(j, from, to);
                        t_ing += ingest(s, d, from, to);
                        acc->merge(s);
                    }
                    acc->finalize();
                    o.init(qs.size(), true, false);
                    double t0 = now_sec();
                    for (size_t i = 0; i < qs.size(); i++) {
                        OmniSketch::Est e = acc->estimate_full(qs[i]);
                        o.est[i] = e.scaled; o.ht[i] = e.thm33; o.sd[i] = e.inter; o.cm[i] = e.nmax;
                    }
                    t_q = now_sec() - t0;
                    used_bits = acc->mem_bits_lean();
                    delete acc;
                } else {
                    fprintf(stderr, "merge: unsupported spec %s\n", sp.raw.c_str());
                    continue;
                }
                sp.raw += "@" + std::to_string(parts);
            } else if (sp.kind == "res" || sp.kind == "resj") {
                // resj: reservoir sized as in the journal version's code (64 p + 32 bits per row)
                uint64_t row_bits = sp.kind == "res" ? Reservoir::bits_per_row(d.p) : 64ull * d.p + 32;
                Reservoir s(d.p, std::max<uint64_t>(1, bits / row_bits), seed * 7919 + 1);
                t_ing = ingest(s, d);
                s.finalize();
                double t0 = now_sec();
                for (size_t i = 0; i < qs.size(); i++) o.est[i] = s.estimate(qs[i]);
                t_q = now_sec() - t0;
                used_bits = s.M * row_bits;
                info = "M=" + std::to_string(s.M);
            } else if (sp.kind == "ash" || sp.kind == "ashl") {
                // ash: journal accounting (32 p + 384 bits per entry); ashl: lean (32 p + 64)
                uint64_t eb = sp.kind == "ash" ? ASH::bits_per_entry_journal(d.p) : ASH::bits_per_entry_lean(d.p);
                ASH s(d.p, std::max<uint64_t>(4, bits / eb), 0.1, seed * 7919 + 5);
                t_ing = ingest(s, d);
                s.finalize();
                double t0 = now_sec();
                for (size_t i = 0; i < qs.size(); i++) o.est[i] = s.estimate(qs[i]);
                t_q = now_sec() - t0;
                used_bits = s.maxSize * eb;
                info = "max=" + std::to_string(s.maxSize) + ";held=" + std::to_string(s.size) + ";ejects=" + std::to_string(s.ejects);
            } else if (sp.kind == "uss") {
                USS s(d.p, std::max<uint64_t>(1, bits / USS::bits_per_row(d.p)), seed * 7919 + 2);
                t_ing = ingest(s, d);
                double t0 = now_sec();
                for (size_t i = 0; i < qs.size(); i++) o.est[i] = s.estimate(qs[i]);
                t_q = now_sec() - t0;
                used_bits = s.mem_bits();
                info = "m=" + std::to_string(s.m) + ";held=" + std::to_string(s.size);
            } else if (sp.kind == "coco") {
                Coco s(d.p, std::max<uint64_t>(2, bits / Coco::bits_per_row(d.p)), seed * 7919 + 7);
                t_ing = ingest(s, d);
                double t0 = now_sec();
                for (size_t i = 0; i < qs.size(); i++) o.est[i] = s.estimate(qs[i]);
                t_q = now_sec() - t0;
                used_bits = s.mem_bits();
                info = "buckets=" + std::to_string(s.entries()) + ";held=" + std::to_string(s.held());
            } else if (sp.kind == "ussc" || sp.kind == "cococ") {
                auto run = [&](auto* s) {
                    t_ing = ingest(*s, d);
                    s->finalize();
                    double t0 = now_sec();
                    for (size_t i = 0; i < qs.size(); i++) o.est[i] = s->estimate(qs[i]);
                    t_q = now_sec() - t0;
                    used_bits = s->mem_bits();
                    delete s;
                };
                if (sp.kind == "ussc") run(make_tcal<USS>(sp, d, bits, seed * 7919 + 2, cx, info));
                else run(make_tcal<Coco>(sp, d, bits, seed * 7919 + 7, cx, info));
            } else if (sp.kind == "omnij" || sp.kind == "omnijs") {
                bool scaled_first = sp.kind == "omnijs";
                OmniJava s(d.p, (int64_t)bits, (int32_t)seed);
                if (s.B < 1) { fprintf(stderr, "budget too small for omnij\n"); continue; }
                t_ing = ingest(s, d);
                o.init(qs.size(), true, false);
                double t0 = now_sec();
                for (size_t i = 0; i < qs.size(); i++) {
                    OmniJava::Est e = s.estimate(qs[i]);
                    o.est[i] = scaled_first ? e.case2 : e.thm33;
                    o.ht[i] = scaled_first ? e.thm33 : e.case2;
                    o.sd[i] = e.inter; o.cm[i] = e.nmax;
                }
                t_q = now_sec() - t0;
                used_bits = s.mem_bits();
                char buf[128];
                snprintf(buf, sizeof buf, "w=%d;d=%d;B=%d;b=%d", s.width, s.depth, s.B, s.b);
                info = buf;
            } else if (sp.kind == "omnid") {
                if (!rg.any) { fprintf(stderr, "omnid: no range predicates in the workload\n"); continue; }
                OmniDyad* s = make_omnid(sp, d, rg, bits, seed * 7919 + 3, info);
                t_ing = ingest(*s, d);
                s->finalize();
                o.init(qs.size(), true, false);
                double t0 = now_sec();
                for (size_t i = 0; i < qs.size(); i++) {
                    OmniDyad::Est e = s->estimate(qs[i]);
                    o.est[i] = e.code; o.ht[i] = e.thr; o.sd[i] = e.inter; o.cm[i] = 0;
                }
                t_q = now_sec() - t0;
                used_bits = s->mem_bits_lean();
                delete s;
            } else if (sp.kind.rfind("omni", 0) == 0) {
                OmniSketch* s = make_omni(sp, d, bits, seed * 7919 + 3, info);
                t_ing = ingest(*s, d);
                s->finalize();
                o.init(qs.size(), true, false);
                double t0 = now_sec();
                for (size_t i = 0; i < qs.size(); i++) {
                    OmniSketch::Est e = s->estimate_full(qs[i]);
                    o.est[i] = e.scaled; o.ht[i] = e.thm33; o.sd[i] = e.inter; o.cm[i] = e.nmax;
                }
                t_q = now_sec() - t0;
                used_bits = sp.kind == "omnip" ? OmniSketch::paper_bits(s->p, s->d, s->w, s->B)
                          : sp.kind == "omniv" ? (uint64_t)s->p * s->d * s->w * (31 * s->B + 32)
                          : s->mem_bits_lean();
                delete s;
            } else if (sp.kind == "csh") {
                CoStrataH* s = make_csh(sp, d, bits, seed * 7919 + 6, cx, info);
                if (!s) { fprintf(stderr, "budget too small for %s\n", sp.raw.c_str()); continue; }
                t_ing = ingest(*s, d);
                s->finalize();
                o.init(qs.size(), true, true);
                double t0 = now_sec();
                for (size_t i = 0; i < qs.size(); i++) {
                    CoStrataH::Est e = s->estimate(qs[i]);
                    o.est[i] = e.hajek; o.ht[i] = e.ht; o.sd[i] = e.sd_hajek; o.cm[i] = e.cm; o.m[i] = e.m; o.wmax[i] = e.wmax;
                }
                t_q = now_sec() - t0;
                used_bits = s->mem_bits();
                char buf[200];
                snprintf(buf, sizeof buf, ";held=%llu;B=%.4g;compactions=%llu;rounds=%llu;E=%llu", (unsigned long long)s->size, s->Bcur,
                         (unsigned long long)s->compactions, (unsigned long long)s->rounds, (unsigned long long)s->E);
                info += buf;
                delete s;
            } else if (sp.kind == "cs") {
                CoStrata* s = make_cs(sp, d, bits, seed * 7919 + 4, cx, info);
                if (!s) { fprintf(stderr, "budget too small for %s\n", sp.raw.c_str()); continue; }
                t_ing = ingest(*s, d);
                s->finalize();
                double ixo = sp.num("ix", sp.a.size() > 6 ? sp.a[6] : 0);
                if (ixo > 0) s->build_index(); // query-time bucket index
                o.init(qs.size(), true, true);
                double t0 = now_sec();
                bool rake = sp.a.size() > 5 && sp.a[5] > 0; // calibrate to all anchors (raking)
                for (size_t i = 0; i < qs.size(); i++) {
                    CoStrata::Est e = s->estimate(qs[i]);
                    o.est[i] = rake ? s->estimate_rake(qs[i]) : e.hajek;
                    o.ht[i] = e.ht; o.sd[i] = e.sd_hajek; o.cm[i] = e.cm; o.m[i] = e.m; o.wmax[i] = e.wmax;
                }
                t_q = now_sec() - t0;
                used_bits = s->mem_bits() + (ixo >= 2 ? s->index_bits() : 0);
                char buf[200];
                snprintf(buf, sizeof buf, ";held=%llu;B=%.4g;compactions=%llu;E=%llu;index_mb=%.2f", (unsigned long long)s->size, s->Bcur,
                         (unsigned long long)s->compactions, (unsigned long long)s->E, (double)s->index_bits() / 8e6);
                info += buf;
                delete s;
            } else {
                fprintf(stderr, "unknown spec %s\n", sp.raw.c_str());
                continue;
            }
            for (size_t i = 0; i < qs.size(); i++) {
                if (o.ht.empty()) fprintf(out, "%s,%s,%g,%llu,%zu,%.12g,,,\n", name.c_str(), sp.raw.c_str(), mem_mb, (unsigned long long)seed, i, o.est[i]);
                else if (o.m.empty())
                    fprintf(out, "%s,%s,%g,%llu,%zu,%.12g,%.12g,%.12g,%.12g\n", name.c_str(), sp.raw.c_str(), mem_mb, (unsigned long long)seed, i,
                            o.est[i], o.ht[i], o.sd[i], o.cm[i]);
                else
                    fprintf(out, "%s,%s,%g,%llu,%zu,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g\n", name.c_str(), sp.raw.c_str(), mem_mb,
                            (unsigned long long)seed, i, o.est[i], o.ht[i], o.sd[i], o.cm[i], o.m[i], o.wmax[i]);
            }
            printf("RUN,%s,%s,%g,%llu,%u,%llu,%.3f,%.2f,%.2f,%.3f,%s\n", name.c_str(), sp.raw.c_str(), mem_mb,
                   (unsigned long long)seed, d.p, (unsigned long long)d.n, t_ing, t_ing / (double)d.n * 1e9,
                   t_q / (double)qs.size() * 1e6, (double)used_bits / 8e6, info.c_str());
            fflush(stdout);
        }
        fclose(out);
    } else if (mode == "time") {
        if (argc < 8) { fprintf(stderr, "bad args\n"); return 1; }
        if (!d.load(argv[2])) { fprintf(stderr, "cannot load data\n"); return 1; }
        std::string name = argv[3];
        double mem_mb = atof(argv[4]);
        uint64_t seed = strtoull(argv[5], 0, 10);
        int reps = atoi(argv[6]);
        uint64_t bits = (uint64_t)(mem_mb * 8e6);
        int ai = 7;
        if (ai < argc && std::string(argv[ai]) == "-p") { d.project_prefix((uint32_t)atoi(argv[ai + 1])); ai += 2; }
        // optional workload (-w file): needed for range attributes (dy=1, omnid) and psel=wl;
        // without it pair selection by workload falls back to an empty ranking
        std::vector<Query> qs;
        Ranges rg;
        Ctx cx;
        if (ai + 1 < argc && std::string(argv[ai]) == "-w") {
            if (!load_workload(argv[ai + 1], qs)) { fprintf(stderr, "cannot load workload\n"); return 1; }
            rg.init(d, qs);
            cx.qs = &qs;
            cx.rg = &rg;
            ai += 2;
        }
        for (; ai < argc; ai++) {
            Spec sp = parse_spec(argv[ai]);
            std::vector<double> ts;
            std::string info;
            for (int r = 0; r < reps; r++) {
                if (sp.kind == "res") {
                    Reservoir s(d.p, std::max<uint64_t>(1, bits / Reservoir::bits_per_row(d.p)), seed + r);
                    ts.push_back(ingest(s, d));
                } else if (sp.kind == "uss") {
                    USS s(d.p, std::max<uint64_t>(1, bits / USS::bits_per_row(d.p)), seed + r);
                    ts.push_back(ingest(s, d));
                } else if (sp.kind == "coco") {
                    Coco s(d.p, std::max<uint64_t>(2, bits / Coco::bits_per_row(d.p)), seed + r);
                    ts.push_back(ingest(s, d));
                } else if (sp.kind == "omnij" || sp.kind == "omnijs") {
                    OmniJava s(d.p, (int64_t)bits, (int32_t)(seed + r));
                    ts.push_back(ingest(s, d));
                } else if (sp.kind == "omnid") {
                    if (!rg.any) { fprintf(stderr, "omnid: needs a workload with range predicates (-w)\n"); break; }
                    OmniDyad* s = make_omnid(sp, d, rg, bits, seed + r, info);
                    ts.push_back(ingest(*s, d));
                    delete s;
                } else if (sp.kind.rfind("omni", 0) == 0) {
                    OmniSketch* s = make_omni(sp, d, bits, seed + r, info);
                    ts.push_back(ingest(*s, d));
                    delete s;
                } else if (sp.kind == "ash" || sp.kind == "ashl") {
                    uint64_t eb = sp.kind == "ash" ? ASH::bits_per_entry_journal(d.p) : ASH::bits_per_entry_lean(d.p);
                    ASH s(d.p, std::max<uint64_t>(4, bits / eb), 0.1, seed + r);
                    ts.push_back(ingest(s, d));
                } else if (sp.kind == "ussc") {
                    auto* s = make_tcal<USS>(sp, d, bits, seed + r, cx, info);
                    ts.push_back(ingest(*s, d));
                    delete s;
                } else if (sp.kind == "cococ") {
                    auto* s = make_tcal<Coco>(sp, d, bits, seed + r, cx, info);
                    ts.push_back(ingest(*s, d));
                    delete s;
                } else if (sp.kind == "csh") {
                    CoStrataH* s = make_csh(sp, d, bits, seed + r, cx, info);
                    if (!s) break;
                    ts.push_back(ingest(*s, d));
                    info = "compactions=" + std::to_string(s->compactions);
                    delete s;
                } else if (sp.kind == "cs") {
                    CoStrata* s = make_cs(sp, d, bits, seed + r, cx, info);
                    if (!s) break;
                    ts.push_back(ingest(*s, d));
                    info = "compactions=" + std::to_string(s->compactions) + ";npair=" + std::to_string(s->mc.n2);
                    delete s;
                }
            }
            if (ts.empty()) continue;
            std::sort(ts.begin(), ts.end());
            double med = ts[ts.size() / 2];
            printf("TIME,%s,%s,%g,%u,%llu,%.2f,%.3f,%s\n", name.c_str(), sp.raw.c_str(), mem_mb, d.p,
                   (unsigned long long)d.n, med / (double)d.n * 1e9, (double)d.n / med / 1e6, info.c_str());
            fflush(stdout);
        }
    }
    return 0;
}
