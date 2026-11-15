// Benchmark driver.
//   bench acc  <data.bin> <workload.txt> <out.csv> <name> <mem_mb> <seed> [-p k] <spec>...
//   bench time <data.bin> <name> <mem_mb> <seed> <reps> [-p k] <spec>...
// Method specs (all sized to the same memory budget of mem_mb * 8e6 bits):
//   res            uniform row reservoir
//   uss            Unbiased Space-Saving on full tuples
//   omnip          OmniSketch, paper parameters (eps=delta=0.1) and paper memory accounting
//   omnil          OmniSketch, paper (w,d), lean accounting (32 bits/sample, dense)
//   omni:w:d       OmniSketch with given (w,d); B chosen by an oracle so that the
//                  *actually stored* samples fill the budget (most favourable accounting)
//   omniv          OmniSketch as configured in the journal version's code (w=28, d=3, 31 bits/hash)
//   ash / ashl     the journal version's sampling baseline (adaptive Sample-and-Hold on full
//                  tuples) with its accounting (32p+384 bits/entry) / lean accounting (32p+64)
//   resj           reservoir sized as in the journal version's code (64p+32 bits per row)
//   omnijs         as omnij, but "est" holds the always-scaled estimate and "ht" the Theorem 3.3 one
//   omnij          bit-exact port of the authors' Java implementation (hashes, parameter search,
//                  accounting, Theorem 3.3 estimator); the seed is used as their "repetition"
//   cs:w:lf:gamma[:phi2[:phi3[:rake]]]  CoStrata with w buckets/attribute, house weight lf,
//                  exponent gamma; phi2 / phi3 > 0 add pair / triple counter tables using that
//                  fraction of the budget; rake = 1 calibrates to all anchors of a query (IPF)
//                  instead of the lightest one; index = 1 builds a query-time bucket index
//                  (spec: cs:w:lf:gamma:phi2:phi3:rake:index)
//   csh:w:lf       CoStrata-H: CoStrata over distinct tuples with sample-and-hold counts
// Output columns of "acc": name,spec,mem,seed,query,est,ht,sd,cm. For OmniSketch specs the last
// three hold: the Theorem 3.3 estimate (bound when the intersection is small), |S_cap|, n_max;
// for omnij "est" is the Theorem 3.3 estimate and "ht" the always-scaled one.
// "acc" appends one line per (method, query) to out.csv and prints a RUN line.
#include "synopses.hpp"
#include "omni_java_compat.hpp"

struct Spec {
    std::string raw, kind;
    std::vector<double> a;
};
static Spec parse_spec(const std::string& s) {
    Spec sp;
    sp.raw = s;
    size_t pos = s.find(':');
    sp.kind = s.substr(0, pos);
    while (pos != std::string::npos) {
        size_t nx = s.find(':', pos + 1);
        sp.a.push_back(atof(s.substr(pos + 1, nx == std::string::npos ? nx : nx - pos - 1).c_str()));
        pos = nx;
    }
    return sp;
}

template <class S>
static double ingest(S& syn, const Dataset& d) {
    double t0 = now_sec();
    for (uint64_t i = 0; i < d.n; i++) syn.insert(i, d.row(i));
    return now_sec() - t0;
}

// w = 0 selects the default: the largest power of two such that the counters use <= 5% of the budget.
static uint32_t auto_w(double wspec, uint64_t bits, uint32_t p) {
    if (wspec > 0) return (uint32_t)wspec;
    uint32_t w = 64;
    while (w < 65536 && 32.0 * p * (2.0 * w) <= 0.05 * (double)bits) w *= 2;
    return w;
}

// Widths of the pair / triple calibration tables from the budget fractions given as the
// spec arguments number `at` (pairs) and `at + 1` (triples); 0 disables a level.
static void marg_widths(const Spec& sp, size_t at, uint64_t bits, uint32_t p, uint32_t& w2, uint32_t& w3) {
    double phi2 = sp.a.size() > at ? sp.a[at] : 0.0, phi3 = sp.a.size() > at + 1 ? sp.a[at + 1] : 0.0;
    w2 = phi2 > 0 ? MargCal::width(phi2 * (double)bits, MargCal::npairs(p)) : 0;
    w3 = phi3 > 0 && p <= 64 ? MargCal::width(phi3 * (double)bits, MargCal::ntriples(p)) : 0;
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

int main(int argc, char** argv) {
    if (argc < 2) return 1;
    std::string mode = argv[1];
    Dataset d;
    if (mode == "acc") {
        if (argc < 9) { fprintf(stderr, "bad args\n"); return 1; }
        if (!d.load(argv[2])) { fprintf(stderr, "cannot load data\n"); return 1; }
        std::vector<Query> qs;
        if (!load_workload(argv[3], qs)) { fprintf(stderr, "cannot load workload\n"); return 1; }
        FILE* out = fopen(argv[4], "a");
        std::string name = argv[5];
        double mem_mb = atof(argv[6]);
        uint64_t seed = strtoull(argv[7], 0, 10);
        uint64_t bits = (uint64_t)(mem_mb * 8e6);
        int ai = 8;
        if (std::string(argv[ai]) == "-p") { d.project_prefix((uint32_t)atoi(argv[ai + 1])); ai += 2; }
        for (; ai < argc; ai++) {
            Spec sp = parse_spec(argv[ai]);
            std::string info;
            double t_ing = 0, t_q = 0;
            uint64_t used_bits = 0;
            std::vector<double> est(qs.size()), ht, sd, cm;
            if (sp.kind == "res" || sp.kind == "resj") {
                // resj: reservoir sized as in the journal version's code (64 p + 32 bits per row)
                uint64_t row_bits = sp.kind == "res" ? Reservoir::bits_per_row(d.p) : 64ull * d.p + 32;
                Reservoir s(d.p, std::max<uint64_t>(1, bits / row_bits), seed * 7919 + 1);
                t_ing = ingest(s, d);
                s.finalize();
                double t0 = now_sec();
                for (size_t i = 0; i < qs.size(); i++) est[i] = s.estimate(qs[i]);
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
                for (size_t i = 0; i < qs.size(); i++) est[i] = s.estimate(qs[i]);
                t_q = now_sec() - t0;
                used_bits = s.maxSize * eb;
                info = "max=" + std::to_string(s.maxSize) + ";held=" + std::to_string(s.size) + ";ejects=" + std::to_string(s.ejects);
            } else if (sp.kind == "uss") {
                USS s(d.p, std::max<uint64_t>(1, bits / USS::bits_per_row(d.p)), seed * 7919 + 2);
                t_ing = ingest(s, d);
                double t0 = now_sec();
                for (size_t i = 0; i < qs.size(); i++) est[i] = s.estimate(qs[i]);
                t_q = now_sec() - t0;
                used_bits = s.mem_bits();
                info = "m=" + std::to_string(s.m) + ";held=" + std::to_string(s.size);
            } else if (sp.kind == "ussc") {
                uint32_t w = auto_w(sp.a.empty() ? 0 : sp.a[0], bits, d.p);
                uint32_t w2, w3;
                marg_widths(sp, 1, bits, d.p, w2, w3);
                USSC s(d.p, w, std::max<uint64_t>(1, USSC::rows_for_budget(bits, d.p, w, MargCal::bits_for(d.p, w2, w3))),
                       seed * 7919 + 2);
                if (w2 || w3) s.enable_marginals(w2, w3);
                t_ing = ingest(s, d);
                s.finalize();
                double t0 = now_sec();
                for (size_t i = 0; i < qs.size(); i++) est[i] = s.estimate(qs[i]);
                t_q = now_sec() - t0;
                used_bits = s.mem_bits();
                info = "m=" + std::to_string(s.uss.m) + ";w=" + std::to_string(w) + ";w2=" + std::to_string(w2) +
                       ";w3=" + std::to_string(w3);
            } else if (sp.kind == "omnij" || sp.kind == "omnijs") {
                bool scaled_first = sp.kind == "omnijs";
                OmniJava s(d.p, (int64_t)bits, (int32_t)seed);
                if (s.B < 1) { fprintf(stderr, "budget too small for omnij\n"); continue; }
                t_ing = ingest(s, d);
                ht.resize(qs.size()); sd.resize(qs.size()); cm.resize(qs.size());
                double t0 = now_sec();
                for (size_t i = 0; i < qs.size(); i++) {
                    OmniJava::Est e = s.estimate(qs[i]);
                    est[i] = scaled_first ? e.case2 : e.thm33;
                    ht[i] = scaled_first ? e.thm33 : e.case2;
                    sd[i] = e.inter; cm[i] = e.nmax;
                }
                t_q = now_sec() - t0;
                used_bits = s.mem_bits();
                char buf[128];
                snprintf(buf, sizeof buf, "w=%d;d=%d;B=%d;b=%d", s.width, s.depth, s.B, s.b);
                info = buf;
            } else if (sp.kind.rfind("omni", 0) == 0) {
                OmniSketch* s = make_omni(sp, d, bits, seed * 7919 + 3, info);
                t_ing = ingest(*s, d);
                s->finalize();
                ht.resize(qs.size()); sd.resize(qs.size()); cm.resize(qs.size());
                double t0 = now_sec();
                for (size_t i = 0; i < qs.size(); i++) {
                    OmniSketch::Est e = s->estimate_full(qs[i]);
                    est[i] = e.scaled; ht[i] = e.thm33; sd[i] = e.inter; cm[i] = e.nmax;
                }
                t_q = now_sec() - t0;
                used_bits = sp.kind == "omnip" ? OmniSketch::paper_bits(s->p, s->d, s->w, s->B)
                          : sp.kind == "omniv" ? (uint64_t)s->p * s->d * s->w * (31 * s->B + 32)
                          : s->mem_bits_lean();
                delete s;
            } else if (sp.kind == "csh") {
                uint32_t w = auto_w(sp.a.empty() ? 0 : sp.a[0], bits, d.p);
                double lf = sp.a.size() > 1 ? sp.a[1] : 1.0;
                uint32_t w2, w3;
                marg_widths(sp, 2, bits, d.p, w2, w3);
                uint64_t M = CoStrataH::rows_for_budget(bits, d.p, w, MargCal::bits_for(d.p, w2, w3));
                if (M < 16) { fprintf(stderr, "budget too small for %s\n", sp.raw.c_str()); continue; }
                CoStrataH s(d.p, w, M, lf, seed * 7919 + 6);
                if (w2 || w3) s.enable_marginals(w2, w3);
                t_ing = ingest(s, d);
                s.finalize();
                ht.resize(qs.size()); sd.resize(qs.size()); cm.resize(qs.size());
                double t0 = now_sec();
                for (size_t i = 0; i < qs.size(); i++) {
                    CoStrataH::Est e = s.estimate(qs[i]);
                    est[i] = e.hajek; ht[i] = e.ht; sd[i] = e.sd_hajek; cm[i] = e.cm;
                }
                t_q = now_sec() - t0;
                used_bits = s.mem_bits();
                char buf[200];
                snprintf(buf, sizeof buf, "w=%u;M=%llu;held=%llu;B=%.4g;compactions=%llu;rounds=%llu;E=%llu", w,
                         (unsigned long long)M, (unsigned long long)s.size, s.Bcur, (unsigned long long)s.compactions,
                         (unsigned long long)s.rounds, (unsigned long long)s.E);
                info = buf;
            } else if (sp.kind == "cs") {
                uint32_t w = auto_w(sp.a.empty() ? 0 : sp.a[0], bits, d.p);
                double lf = sp.a.size() > 1 ? sp.a[1] : 1.0, gamma = sp.a.size() > 2 ? sp.a[2] : 1.0;
                uint32_t w2, w3;
                marg_widths(sp, 3, bits, d.p, w2, w3);
                uint64_t M = CoStrata::rows_for_budget(bits, d.p, w, MargCal::bits_for(d.p, w2, w3));
                if (M < 16) { fprintf(stderr, "budget too small for %s\n", sp.raw.c_str()); continue; }
                CoStrata s(d.p, w, M, lf, gamma, seed * 7919 + 4);
                if (w2 || w3) s.enable_marginals(w2, w3);
                t_ing = ingest(s, d);
                s.finalize();
                bool use_index = sp.a.size() > 6 && sp.a[6] > 0; // query-time bucket index
                if (use_index) s.build_index();
                ht.resize(qs.size()); sd.resize(qs.size()); cm.resize(qs.size());
                double t0 = now_sec();
                bool rake = sp.a.size() > 5 && sp.a[5] > 0; // calibrate to all anchors (raking)
                for (size_t i = 0; i < qs.size(); i++) {
                    CoStrata::Est e = s.estimate(qs[i]);
                    est[i] = rake ? s.estimate_rake(qs[i]) : e.hajek;
                    ht[i] = e.ht; sd[i] = e.sd_hajek; cm[i] = e.cm;
                }
                t_q = now_sec() - t0;
                used_bits = s.mem_bits();
                char buf[160];
                snprintf(buf, sizeof buf, "w=%u;w2=%u;w3=%u;M=%llu;held=%llu;B=%.4g;compactions=%llu;E=%llu;index_mb=%.2f", w, w2, w3,
                         (unsigned long long)M, (unsigned long long)s.size, s.Bcur, (unsigned long long)s.compactions,
                         (unsigned long long)s.E, (double)s.index_bits() / 8e6);
                info = buf;
            } else {
                fprintf(stderr, "unknown spec %s\n", sp.raw.c_str());
                continue;
            }
            for (size_t i = 0; i < qs.size(); i++) {
                if (ht.empty()) fprintf(out, "%s,%s,%g,%llu,%zu,%.12g,,,\n", name.c_str(), sp.raw.c_str(), mem_mb, (unsigned long long)seed, i, est[i]);
                else fprintf(out, "%s,%s,%g,%llu,%zu,%.12g,%.12g,%.12g,%.12g\n", name.c_str(), sp.raw.c_str(), mem_mb, (unsigned long long)seed, i, est[i], ht[i], sd[i], cm[i]);
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
        if (std::string(argv[ai]) == "-p") { d.project_prefix((uint32_t)atoi(argv[ai + 1])); ai += 2; }
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
                } else if (sp.kind == "omnij" || sp.kind == "omnijs") {
                    OmniJava s(d.p, (int64_t)bits, (int32_t)(seed + r));
                    ts.push_back(ingest(s, d));
                } else if (sp.kind.rfind("omni", 0) == 0) {
                    OmniSketch* s = make_omni(sp, d, bits, seed + r, info);
                    ts.push_back(ingest(*s, d));
                    delete s;
                } else if (sp.kind == "ash" || sp.kind == "ashl") {
                    uint64_t eb = sp.kind == "ash" ? ASH::bits_per_entry_journal(d.p) : ASH::bits_per_entry_lean(d.p);
                    ASH s(d.p, std::max<uint64_t>(4, bits / eb), 0.1, seed + r);
                    ts.push_back(ingest(s, d));
                } else if (sp.kind == "ussc") {
                    uint32_t w = auto_w(sp.a.empty() ? 0 : sp.a[0], bits, d.p);
                    uint32_t w2, w3;
                    marg_widths(sp, 1, bits, d.p, w2, w3);
                    USSC s(d.p, w, std::max<uint64_t>(1, USSC::rows_for_budget(bits, d.p, w, MargCal::bits_for(d.p, w2, w3))), seed + r);
                    if (w2 || w3) s.enable_marginals(w2, w3);
                    ts.push_back(ingest(s, d));
                } else if (sp.kind == "csh") {
                    uint32_t w = auto_w(sp.a.empty() ? 0 : sp.a[0], bits, d.p);
                    double lf = sp.a.size() > 1 ? sp.a[1] : 1.0;
                    uint32_t w2, w3;
                    marg_widths(sp, 2, bits, d.p, w2, w3);
                    CoStrataH s(d.p, w, CoStrataH::rows_for_budget(bits, d.p, w, MargCal::bits_for(d.p, w2, w3)), lf, seed + r);
                    if (w2 || w3) s.enable_marginals(w2, w3);
                    ts.push_back(ingest(s, d));
                    info = "compactions=" + std::to_string(s.compactions);
                } else if (sp.kind == "cs") {
                    uint32_t w = auto_w(sp.a.empty() ? 0 : sp.a[0], bits, d.p);
                    double lf = sp.a.size() > 1 ? sp.a[1] : 1.0, gamma = sp.a.size() > 2 ? sp.a[2] : 1.0;
                    uint32_t w2, w3;
                    marg_widths(sp, 3, bits, d.p, w2, w3);
                    CoStrata s(d.p, w, CoStrata::rows_for_budget(bits, d.p, w, MargCal::bits_for(d.p, w2, w3)), lf, gamma, seed + r);
                    if (w2 || w3) s.enable_marginals(w2, w3);
                    ts.push_back(ingest(s, d));
                    info = "compactions=" + std::to_string(s.compactions);
                }
            }
            std::sort(ts.begin(), ts.end());
            double med = ts[ts.size() / 2];
            printf("TIME,%s,%s,%g,%u,%llu,%.2f,%.3f,%s\n", name.c_str(), sp.raw.c_str(), mem_mb, d.p,
                   (unsigned long long)d.n, med / (double)d.n * 1e9, (double)d.n / med / 1e6, info.c_str());
            fflush(stdout);
        }
    }
    return 0;
}
