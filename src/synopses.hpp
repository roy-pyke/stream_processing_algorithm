// Stream synopses for COUNT queries with ad-hoc conjunctive predicates.
//   Reservoir  : uniform row sample (Vitter's Algorithm L)
//   OmniSketch : Punter, Papapetrou, Garofalakis (PVLDB 17(3), 2023), re-implemented
//   USS        : Unbiased Space-Saving over full attribute tuples (Ting, SIGMOD 2018)
//   CoStrata   : coordinated, hash-stratified row sample with count-driven thresholds (this work)
#pragma once
#include "common.hpp"

// =============================================================== Reservoir
struct Reservoir {
    uint32_t p;
    uint64_t M, N = 0, size = 0, next;
    double W;
    Rng rng;
    std::vector<std::vector<uint32_t>> col;

    Reservoir(uint32_t p_, uint64_t M_, uint64_t seed) : p(p_), M(M_), rng(seed), col(p_) {
        for (auto& c : col) c.resize(M);
        W = std::exp(std::log(rng.unit() + 1e-300) / (double)M);
        next = M + skip();
    }
    uint64_t skip() { return (uint64_t)std::floor(std::log(rng.unit() + 1e-300) / std::log1p(-W)) + 1; }
    static uint64_t bits_per_row(uint32_t p) { return 32ull * p; }
    uint64_t mem_bits() const { return M * bits_per_row(p); }

    inline void insert(uint64_t /*rid*/, const uint32_t* rec) {
        uint64_t i = N++;
        if (i < M) {
            for (uint32_t a = 0; a < p; a++) col[a][i] = rec[a];
            size = i + 1;
        } else if (i + 1 == next) {
            uint64_t s = rng.below(M);
            for (uint32_t a = 0; a < p; a++) col[a][s] = rec[a];
            W *= std::exp(std::log(rng.unit() + 1e-300) / (double)M);
            next += skip();
        }
    }
    void finalize() {}
    double estimate(const Query& q) const {
        if (size == 0) return 0;
        const uint32_t* c0 = col[q.attrs[0]].data();
        uint32_t v0 = q.vals[0];
        uint64_t hits = 0;
        for (uint64_t s = 0; s < size; s++) {
            if (c0[s] != v0) continue;
            bool ok = true;
            for (size_t k = 1; k < q.attrs.size(); k++)
                if (col[q.attrs[k]][s] != q.vals[k]) { ok = false; break; }
            hits += ok;
        }
        return (double)hits * (double)N / (double)size;
    }
};

// =============================================================== OmniSketch
// p attribute sketches, each d x w cells; a cell keeps a counter and the B
// smallest hashes g(rid) of the records hashed into it (bottom-B / K-minwise).
// Estimator (Eq. 4 of the paper): n_max / B * | intersection of the p_q * d samples |.
// Engineering beyond the reference Java code: bounded binary max-heaps in flat
// vectors instead of red-black trees, and a dense threshold array so that the
// common "reject" path costs one comparison per cell.
struct OmniSketch {
    uint32_t p, d, w;
    uint64_t B;
    uint64_t N = 0, gseed;
    std::vector<uint64_t> hseed;             // p*d
    std::vector<uint32_t> cnt;               // p*d*w
    std::vector<uint64_t> thr;               // p*d*w : current max hash once full, else UINT64_MAX
    std::vector<std::vector<uint64_t>> smp;  // p*d*w : heap during ingest, sorted after finalize()
    bool sorted = false;

    OmniSketch(uint32_t p_, uint32_t d_, uint32_t w_, uint64_t B_, uint64_t seed)
        : p(p_), d(d_), w(w_), B(B_), hseed((size_t)p_ * d_), cnt((size_t)p_ * d_ * w_, 0),
          thr((size_t)p_ * d_ * w_, UINT64_MAX), smp((size_t)p_ * d_ * w_) {
        Rng r(seed);
        gseed = r.next();
        for (auto& s : hseed) s = r.next();
    }
    inline size_t cell(uint32_t a, uint32_t j, uint32_t v) const {
        return ((size_t)a * d + j) * w + bucket_of(v, hseed[(size_t)a * d + j], w);
    }
    // "Lean" accounting: 32 bits per stored hash + one 32-bit counter per cell.
    uint64_t mem_bits_lean() const {
        uint64_t e = 0;
        for (auto& s : smp) e += s.size();
        return 32ull * e + 32ull * cnt.size();
    }
    // Accounting used in the paper / reference code: B*(b + 3*32 + 1) + 32 bits per cell.
    static uint64_t paper_bits(uint32_t p, uint32_t d, uint32_t w, uint64_t B, double delta = 0.1) {
        double b = std::ceil(std::log2(4.0 * std::pow((double)B, 2.5) / delta));
        return (uint64_t)((double)p * d * w * ((double)B * (b + 97.0) + 32.0));
    }
    // Parameters prescribed by the paper for (eps, delta).
    static void paper_wd(double eps, double delta, uint32_t& w, uint32_t& d) {
        d = (uint32_t)std::ceil(std::log(2.0 / delta));
        w = 1 + (uint32_t)std::ceil(std::exp(1.0) * std::pow((1.0 + eps) / eps, 1.0 / d));
    }

    inline void insert(uint64_t rid, const uint32_t* rec) {
        N++;
        uint64_t g = mix64(rid ^ gseed) >> 1;
        for (uint32_t a = 0; a < p; a++) {
            for (uint32_t j = 0; j < d; j++) {
                size_t c = cell(a, j, rec[a]);
                cnt[c]++;
                if (g >= thr[c]) continue;
                auto& h = smp[c];
                if (h.size() < B) {
                    h.push_back(g);
                    std::push_heap(h.begin(), h.end());
                    if (h.size() == B) thr[c] = h.front();
                } else {
                    std::pop_heap(h.begin(), h.end());
                    h.back() = g;
                    std::push_heap(h.begin(), h.end());
                    thr[c] = h.front();
                }
            }
        }
    }
    void finalize() {
        for (auto& h : smp) std::sort(h.begin(), h.end());
        sorted = true;
    }
    double estimate(const Query& q) const {
        std::vector<size_t> cells;
        uint32_t nmax = 0;
        for (size_t k = 0; k < q.attrs.size(); k++)
            for (uint32_t j = 0; j < d; j++) {
                size_t c = cell(q.attrs[k], j, q.vals[k]);
                cells.push_back(c);
                nmax = std::max(nmax, cnt[c]);
            }
        std::sort(cells.begin(), cells.end(), [&](size_t x, size_t y) { return smp[x].size() < smp[y].size(); });
        std::vector<uint64_t> cand = smp[cells[0]], tmp;
        for (size_t i = 1; i < cells.size() && !cand.empty(); i++) {
            const auto& s = smp[cells[i]];
            tmp.clear();
            std::set_intersection(cand.begin(), cand.end(), s.begin(), s.end(), std::back_inserter(tmp));
            cand.swap(tmp);
        }
        double scale = std::max(1.0, (double)nmax / (double)B);
        return scale * (double)cand.size();
    }
};

// =============================================================== USS
// Unbiased Space-Saving on the full p-attribute tuple. Subset sums (any
// predicate over the tuple) are estimated by summing the counters of the
// matching monitored tuples.
struct USS {
    uint32_t p;
    uint64_t m, size = 0, N = 0, tmask;
    Rng rng;
    std::vector<uint32_t> keys, cnt, heap, hpos;
    std::vector<int32_t> table;
    uint64_t hseed;

    USS(uint32_t p_, uint64_t m_, uint64_t seed) : p(p_), m(m_), rng(seed) {
        keys.resize(m * p);
        cnt.assign(m, 0);
        heap.resize(m);
        hpos.resize(m);
        uint64_t cap = 16;
        while (cap < 2 * m) cap <<= 1;
        table.assign(cap, -1);
        tmask = cap - 1;
        hseed = rng.next();
    }
    static uint64_t bits_per_row(uint32_t p) { return 32ull * p + 64; }
    uint64_t mem_bits() const { return m * bits_per_row(p); }
    inline uint64_t hkey(const uint32_t* r) const {
        uint64_t h = hseed;
        for (uint32_t a = 0; a < p; a++) h = (h ^ r[a]) * 0x9E3779B97F4A7C15ULL, h ^= h >> 29;
        return mix64(h);
    }
    inline bool eq(uint32_t s, const uint32_t* r) const { return memcmp(&keys[(size_t)s * p], r, 4 * p) == 0; }
    void sift_down(uint64_t i) {
        uint32_t s = heap[i];
        uint32_t c = cnt[s];
        for (;;) {
            uint64_t l = 2 * i + 1;
            if (l >= size) break;
            if (l + 1 < size && cnt[heap[l + 1]] < cnt[heap[l]]) l++;
            if (cnt[heap[l]] >= c) break;
            heap[i] = heap[l];
            hpos[heap[i]] = (uint32_t)i;
            i = l;
        }
        heap[i] = s;
        hpos[s] = (uint32_t)i;
    }
    void sift_up(uint64_t i) {
        uint32_t s = heap[i];
        uint32_t c = cnt[s];
        while (i > 0) {
            uint64_t par = (i - 1) / 2;
            if (cnt[heap[par]] <= c) break;
            heap[i] = heap[par];
            hpos[heap[i]] = (uint32_t)i;
            i = par;
        }
        heap[i] = s;
        hpos[s] = (uint32_t)i;
    }
    void tab_insert(uint32_t s) {
        uint64_t i = hkey(&keys[(size_t)s * p]) & tmask;
        while (table[i] >= 0) i = (i + 1) & tmask;
        table[i] = (int32_t)s;
    }
    void tab_erase(uint32_t s) {
        uint64_t i = hkey(&keys[(size_t)s * p]) & tmask;
        while (table[i] != (int32_t)s) i = (i + 1) & tmask;
        // backward-shift deletion for linear probing
        uint64_t j = i;
        for (;;) {
            j = (j + 1) & tmask;
            if (table[j] < 0) break;
            uint64_t k = hkey(&keys[(size_t)table[j] * p]) & tmask;
            bool movable = (i <= j) ? (k <= i || k > j) : (k <= i && k > j);
            if (movable) {
                table[i] = table[j];
                i = j;
            }
        }
        table[i] = -1;
    }
    inline void insert(uint64_t /*rid*/, const uint32_t* rec) {
        N++;
        uint64_t i = hkey(rec) & tmask;
        while (table[i] >= 0) {
            if (eq((uint32_t)table[i], rec)) {
                uint32_t s = (uint32_t)table[i];
                cnt[s]++;
                sift_down(hpos[s]);
                return;
            }
            i = (i + 1) & tmask;
        }
        if (size < m) {
            uint32_t s = (uint32_t)size++;
            memcpy(&keys[(size_t)s * p], rec, 4 * p);
            cnt[s] = 1;
            heap[s] = s;
            hpos[s] = s;
            sift_up(s);
            table[i] = (int32_t)s;
            return;
        }
        uint32_t s = heap[0];
        cnt[s]++;
        if (rng.unit() * (double)cnt[s] < 1.0) {
            tab_erase(s);
            memcpy(&keys[(size_t)s * p], rec, 4 * p);
            tab_insert(s);
        }
        sift_down(0);
    }
    void finalize() {}
    double estimate(const Query& q) const {
        uint64_t sum = 0;
        for (uint64_t s = 0; s < size; s++)
            if (q.matches(&keys[s * p])) sum += cnt[s];
        return (double)sum;
    }
};

// =============================================================== CoStrata
// One Count-Min-style row of w exact counters per attribute defines p hash
// stratifications of the stream (plus the trivial stratification {all}).
// A record r with uniform priority u_r in [0,1) is retained iff
//        u_r * nstar_r < Bcur ,   nstar_r = min( G , min_a s(n_{a,h_a(r_a)}) ),
// i.e. it is sampled with probability pi_r = min(1, Bcur / nstar_r): the
// *maximum* over the sampling rates of all strata it belongs to. Thresholds
// are functions of counters only (no per-stratum heaps); the single global
// parameter Bcur is lowered by selection whenever the row store is full.
struct CoStrata {
    uint32_t p, w;
    uint64_t M;          // row capacity
    double lf;           // weight of the global ("house") stratum relative to all hash strata; 0 disables
    double gamma;        // allocation exponent: s(n) = n^gamma (1 = equal allocation per stratum)
    double slack;        // fraction of the store freed by a compaction
    uint64_t N = 0, size = 0, E = 0, compactions = 0;
    double Bcur = std::numeric_limits<double>::infinity();
    double G = 0;        // running max of N / (lf * E)
    uint64_t useed;
    std::vector<uint64_t> hseed;                // p
    std::vector<uint32_t> cnt;                  // p*w
    std::vector<std::vector<uint32_t>> col;     // p columns of M
    std::vector<double> u;                      // M priorities
    // query-time working state (built by finalize)
    std::vector<double> invpi, Wb, Vb;

    CoStrata(uint32_t p_, uint32_t w_, uint64_t M_, double lf_, double gamma_, uint64_t seed, double slack_ = 0.1)
        : p(p_), w(w_), M(M_), lf(lf_), gamma(gamma_), slack(slack_), hseed(p_), cnt((size_t)p_ * w_, 0), col(p_), u(M_) {
        Rng r(seed);
        useed = r.next();
        for (auto& s : hseed) s = r.next();
        for (auto& c : col) c.resize(M);
    }
    static uint64_t bits_per_row(uint32_t p) { return 32ull * p + 64; }
    uint64_t mem_bits() const { return M * bits_per_row(p) + 32ull * cnt.size(); }
    static uint64_t rows_for_budget(uint64_t bits, uint32_t p, uint32_t w) {
        uint64_t fixed = 32ull * p * w;
        return bits > fixed ? (bits - fixed) / bits_per_row(p) : 0;
    }
    inline double sfun(double n) const { return gamma == 1.0 ? n : std::pow(n, gamma); }
    inline double nstar_slot(uint64_t s) const {
        uint32_t nmin = UINT32_MAX;
        for (uint32_t a = 0; a < p; a++) nmin = std::min(nmin, cnt[(size_t)a * w + bucket_of(col[a][s], hseed[a], w)]);
        double v = sfun((double)nmin);
        return (lf > 0 && G < v) ? G : v;
    }
    inline void insert(uint64_t rid, const uint32_t* rec) {
        N++;
        uint32_t nmin = UINT32_MAX;
        for (uint32_t a = 0; a < p; a++) {
            uint32_t& c = cnt[(size_t)a * w + bucket_of(rec[a], hseed[a], w)];
            if (c == 0) E++;
            c++;
            nmin = std::min(nmin, c);
        }
        double ns = sfun((double)nmin);
        if (lf > 0) {
            double g = (double)N / (lf * (double)E);
            if (g > G) G = g;
            if (G < ns) ns = G;
        }
        double ur = to_unit(mix64(rid ^ useed));
        if (ur * ns < Bcur) {
            uint64_t s = size++;
            for (uint32_t a = 0; a < p; a++) col[a][s] = rec[a];
            u[s] = ur;
            if (size == M) compact();
        }
    }
    void compact() {
        compactions++;
        std::vector<double> z(size), zs;
        for (uint64_t s = 0; s < size; s++) z[s] = u[s] * nstar_slot(s);
        uint64_t K = (uint64_t)((1.0 - slack) * (double)M);
        if (K >= size) K = size - 1;
        zs = z;
        std::nth_element(zs.begin(), zs.begin() + K, zs.end());
        if (zs[K] < Bcur) Bcur = zs[K];
        uint64_t o = 0;
        for (uint64_t s = 0; s < size; s++) {
            if (z[s] < Bcur) {
                if (o != s) {
                    for (uint32_t a = 0; a < p; a++) col[a][o] = col[a][s];
                    u[o] = u[s];
                }
                o++;
            }
        }
        size = o;
    }
    // Drop stale rows, compute Horvitz-Thompson weights and per-stratum aggregates.
    void finalize() {
        uint64_t o = 0;
        invpi.assign(size, 0.0);
        Wb.assign((size_t)p * w, 0.0);
        Vb.assign((size_t)p * w, 0.0);
        for (uint64_t s = 0; s < size; s++) {
            double ns = nstar_slot(s);
            if (!(u[s] * ns < Bcur)) continue;
            double pi = std::min(1.0, Bcur / ns);
            if (o != s) {
                for (uint32_t a = 0; a < p; a++) col[a][o] = col[a][s];
                u[o] = u[s];
            }
            invpi[o] = 1.0 / pi;
            for (uint32_t a = 0; a < p; a++) {
                size_t b = (size_t)a * w + bucket_of(col[a][o], hseed[a], w);
                Wb[b] += 1.0 / pi;
                Vb[b] += (1.0 - pi) / (pi * pi);
            }
            o++;
        }
        size = o;
        invpi.resize(size);
    }
    struct Est {
        double hajek, ht, sd_hajek, sd_ht, cm; // cm = Count-Min style upper bound from the counters
    };
    Est estimate(const Query& q) const {
        const uint32_t* c0 = col[q.attrs[0]].data();
        uint32_t v0 = q.vals[0];
        double num = 0, A = 0;
        for (uint64_t s = 0; s < size; s++) {
            if (c0[s] != v0) continue;
            bool ok = true;
            for (size_t k = 1; k < q.attrs.size(); k++)
                if (col[q.attrs[k]][s] != q.vals[k]) { ok = false; break; }
            if (ok) {
                double ip = invpi[s];
                num += ip;
                A += ip * ip - ip; // (1-pi)/pi^2
            }
        }
        // lightest stratum touched by the query
        size_t best = 0;
        uint32_t nb = UINT32_MAX;
        for (size_t k = 0; k < q.attrs.size(); k++) {
            size_t b = (size_t)q.attrs[k] * w + bucket_of(q.vals[k], hseed[q.attrs[k]], w);
            if (cnt[b] < nb) nb = cnt[b], best = b;
        }
        Est e;
        e.ht = num;
        e.sd_ht = std::sqrt(A);
        e.cm = (double)nb;
        if (Wb[best] > 0) {
            double R = num / Wb[best];
            double scale = (double)nb / Wb[best];
            e.hajek = (double)nb * R;
            double var = scale * scale * ((1 - R) * (1 - R) * A + R * R * std::max(0.0, Vb[best] - A));
            e.sd_hajek = std::sqrt(std::max(0.0, var));
        } else {
            e.hajek = 0;
            e.sd_hajek = 0;
        }
        return e;
    }
};

// =============================================================== USS + calibration
// USS over full tuples plus one row of w exact counters per attribute; the
// subset-sum estimate is ratio-calibrated to the exact count of the lightest
// hash stratum touched by the query (same calibration as CoStrata).
struct USSC {
    USS uss;
    uint32_t p, w;
    std::vector<uint64_t> hseed;
    std::vector<uint32_t> cnt;
    std::vector<double> Wb;
    USSC(uint32_t p_, uint32_t w_, uint64_t m, uint64_t seed) : uss(p_, m, seed), p(p_), w(w_), hseed(p_), cnt((size_t)p_ * w_, 0) {
        Rng r(seed ^ 0x5bd1e995);
        for (auto& x : hseed) x = r.next();
    }
    static uint64_t rows_for_budget(uint64_t bits, uint32_t p, uint32_t w) {
        uint64_t fixed = 32ull * p * w;
        return bits > fixed ? (bits - fixed) / USS::bits_per_row(p) : 0;
    }
    uint64_t mem_bits() const { return uss.mem_bits() + 32ull * cnt.size(); }
    inline void insert(uint64_t rid, const uint32_t* rec) {
        for (uint32_t a = 0; a < p; a++) cnt[(size_t)a * w + bucket_of(rec[a], hseed[a], w)]++;
        uss.insert(rid, rec);
    }
    void finalize() {
        Wb.assign((size_t)p * w, 0.0);
        for (uint64_t s = 0; s < uss.size; s++)
            for (uint32_t a = 0; a < p; a++)
                Wb[(size_t)a * w + bucket_of(uss.keys[s * p + a], hseed[a], w)] += uss.cnt[s];
    }
    double estimate(const Query& q) const {
        double num = uss.estimate(q);
        size_t best = 0;
        uint32_t nb = UINT32_MAX;
        for (size_t k = 0; k < q.attrs.size(); k++) {
            size_t b = (size_t)q.attrs[k] * w + bucket_of(q.vals[k], hseed[q.attrs[k]], w);
            if (cnt[b] < nb) nb = cnt[b], best = b;
        }
        return Wb[best] > 0 ? (double)nb * num / Wb[best] : 0.0;
    }
};
