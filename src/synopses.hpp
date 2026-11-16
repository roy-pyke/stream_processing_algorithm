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
        uint64_t hits = 0;
        if (q.has_range()) {
            for (uint64_t s = 0; s < size; s++) {
                bool ok = true;
                for (size_t k = 0; k < q.attrs.size() && ok; k++) {
                    uint32_t x = col[q.attrs[k]][s];
                    ok = x >= q.vals[k] && x <= q.hi[k];
                }
                hits += ok;
            }
            return (double)hits * (double)N / (double)size;
        }
        const uint32_t* c0 = col[q.attrs[0]].data();
        uint32_t v0 = q.vals[0];
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
    // Accounting of the reference code: B*(b + 3*32 + 1) + 32 bits per cell, with
    // b = ceil(ln(4 B^2.5 / delta)) (the Java code uses the natural logarithm).
    static uint64_t paper_bits(uint32_t p, uint32_t d, uint32_t w, uint64_t B, double delta = 0.1) {
        double b = std::ceil(std::log(4.0 * std::pow((double)B, 2.5) / delta));
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
    // Merge a sketch with the same parameters and seeds built over a disjoint part of the stream
    // (before finalize): counters add and every cell keeps the B smallest hashes of the union of the
    // two samples. This is exact: the result is the sketch of the concatenated stream.
    void merge(const OmniSketch& o) {
        N += o.N;
        std::vector<uint64_t> u;
        for (size_t c = 0; c < cnt.size(); c++) {
            cnt[c] += o.cnt[c];
            if (o.smp[c].empty()) continue;
            u.assign(smp[c].begin(), smp[c].end());
            u.insert(u.end(), o.smp[c].begin(), o.smp[c].end());
            std::sort(u.begin(), u.end());
            u.erase(std::unique(u.begin(), u.end()), u.end());
            if (u.size() > B) u.resize(B);
            smp[c].assign(u.begin(), u.end());
            std::make_heap(smp[c].begin(), smp[c].end());
            thr[c] = smp[c].size() == B ? smp[c].front() : UINT64_MAX;
        }
    }
    // scaled: Eq. 4 / Corollary 3.4 (always n_max/B * |S_cap|); thm33: Theorem 3.3, which
    // returns the bound 2 n_max L / (B eps^2), L = ln(4 p_q d sqrt(B)/delta), whenever
    // |S_cap| < 3 L / eps^2 -- this is what the authors' Java query() returns.
    struct Est {
        double scaled, thm33, inter, nmax;
    };
    Est estimate_full(const Query& q, double eps = 0.1, double delta = 0.1) const {
        Est e;
        e.scaled = estimate(q, &e.inter, &e.nmax);
        double L = std::log(4.0 * (double)q.attrs.size() * d * std::sqrt((double)B) / delta);
        e.thm33 = e.inter < 3 * L / (eps * eps) ? std::ceil(2 * e.nmax * L / ((double)B * eps * eps)) : e.scaled;
        return e;
    }
    double estimate(const Query& q, double* inter = nullptr, double* nmax_out = nullptr) const {
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
        if (inter) *inter = (double)cand.size();
        if (nmax_out) *nmax_out = (double)nmax;
        return scale * (double)cand.size();
    }
};

// =============================================================== OmniSketch with dyadic ranges
// Range support of OmniSketch (Section 3.3 of the PVLDB paper; class CountMinDyad of the authors'
// code): every range attribute a has L_a = bits_a + 1 levels, level l being a d x w array of cells
// over the dyadic intervals of length 2^l (interval index v >> l); every other attribute has level 0
// only. Each record is inserted at every level. A predicate lo <= r_a <= hi is decomposed into its
// canonical dyadic cover; for every sketch row the samples of the covering cells are united, and the
// united samples of all predicates and rows are intersected. Two estimators:
//   code : |S_cap| * n_max / B_virtual as in the reference code (OmniSketch.rangeQuery, case 2):
//          n_max is the largest summed count of the covering cells of a predicate in one row, over
//          all predicates and rows, and B_virtual = B times the number of covering pieces of the
//          last predicate of the query;
//   thr  : the united sample of every (predicate, row) is a coordinated sample of the union of its
//          cells up to the smallest threshold among them; the intersection, truncated to the smallest
//          of these thresholds tau*, is a uniform sample of R_cap at rate tau*, scaled by 1/tau*.
// Memory is charged with the lean accounting of `omni:w:d` (32 bits per stored hash and per non-empty
// cell) and B is chosen by the oracle from the final cell loads.
struct OmniDyad {
    uint32_t p, d, w;
    uint64_t B, N = 0, gseed;
    std::vector<uint32_t> L, off;            // levels per attribute, first level id of attribute a
    std::vector<uint64_t> hseed;             // per (level id, row)
    std::vector<uint32_t> cnt;
    std::vector<uint64_t> thr;
    std::vector<std::vector<uint64_t>> smp;
    static constexpr uint64_t G63 = 1ull << 63;

    OmniDyad(uint32_t p_, uint32_t d_, uint32_t w_, uint64_t B_, const std::vector<uint32_t>& levels, uint64_t seed)
        : p(p_), d(d_), w(w_), B(B_), L(levels), off(p_) {
        uint32_t T = 0;
        for (uint32_t a = 0; a < p; a++) off[a] = T, T += L[a];
        Rng r(seed);
        gseed = r.next();
        hseed.resize((size_t)T * d);
        for (auto& s : hseed) s = r.next();
        cnt.assign((size_t)T * d * w, 0);
        thr.assign(cnt.size(), UINT64_MAX);
        smp.resize(cnt.size());
    }
    inline size_t cell(uint32_t a, uint32_t l, uint32_t j, uint64_t x) const {
        size_t lv = (size_t)(off[a] + l) * d + j;
        return lv * w + bucket_of((uint32_t)x, hseed[lv], w);
    }
    // cell loads of a stream, for the oracle choice of B (same hash functions as the sketch)
    std::vector<uint32_t> loads(const Dataset& ds) const {
        std::vector<uint32_t> c(cnt.size(), 0);
        for (uint64_t i = 0; i < ds.n; i++) {
            const uint32_t* r = ds.row(i);
            for (uint32_t a = 0; a < p; a++)
                for (uint32_t l = 0; l < L[a]; l++)
                    for (uint32_t j = 0; j < d; j++) c[cell(a, l, j, (uint64_t)r[a] >> l)]++;
        }
        return c;
    }
    uint64_t mem_bits_lean() const {
        uint64_t e = 0, ne = 0;
        for (size_t c = 0; c < cnt.size(); c++) e += smp[c].size(), ne += cnt[c] > 0;
        return 32ull * (e + ne);
    }
    inline void insert(uint64_t rid, const uint32_t* rec) {
        N++;
        uint64_t g = mix64(rid ^ gseed) >> 1;
        for (uint32_t a = 0; a < p; a++)
            for (uint32_t l = 0; l < L[a]; l++)
                for (uint32_t j = 0; j < d; j++) {
                    size_t c = cell(a, l, j, (uint64_t)rec[a] >> l);
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
    void finalize() {
        for (auto& h : smp) std::sort(h.begin(), h.end());
    }
    struct Est {
        double code, thr, inter;
    };
    Est estimate(const Query& q) const {
        const size_t k = q.attrs.size();
        std::vector<std::vector<uint64_t>> U(k * d);
        std::vector<double> n(k * d, 0.0), pieces(k, 1.0);
        std::vector<uint64_t> tau(k * d, UINT64_MAX);
        std::vector<std::pair<uint32_t, uint64_t>> cov;
        for (size_t i = 0; i < k; i++) {
            uint32_t a = q.attrs[i];
            uint64_t lo = q.vals[i], hi = q.is_eq(i) ? q.vals[i] : q.hi[i];
            if (L[a] > 1) hi = std::min<uint64_t>(hi, (1ull << (L[a] - 1)) - 1);
            if (lo == hi) cov.assign(1, {0u, lo});
            else if (lo < hi) dyadic_cover(lo, hi, cov); // L[a] > 1 for every attribute with range predicates
            else cov.clear();
            pieces[i] = std::max<double>(1.0, (double)cov.size());
            for (uint32_t j = 0; j < d; j++) {
                auto& u = U[i * d + j];
                for (auto& lx : cov) {
                    size_t c = cell(a, lx.first, j, lx.second);
                    n[i * d + j] += cnt[c];
                    tau[i * d + j] = std::min(tau[i * d + j], thr[c]);
                    u.insert(u.end(), smp[c].begin(), smp[c].end());
                }
                if (cov.size() > 1) {
                    std::sort(u.begin(), u.end());
                    u.erase(std::unique(u.begin(), u.end()), u.end());
                }
            }
        }
        std::vector<size_t> ord(k * d);
        for (size_t i = 0; i < ord.size(); i++) ord[i] = i;
        std::sort(ord.begin(), ord.end(), [&](size_t x, size_t y) { return U[x].size() < U[y].size(); });
        std::vector<uint64_t> cand = U[ord[0]], tmp;
        for (size_t i = 1; i < ord.size() && !cand.empty(); i++) {
            tmp.clear();
            std::set_intersection(cand.begin(), cand.end(), U[ord[i]].begin(), U[ord[i]].end(), std::back_inserter(tmp));
            cand.swap(tmp);
        }
        Est e;
        e.inter = (double)cand.size();
        // as OmniSketch.rangeQuery (case-2 estimate): n_max over all predicates and rows, and
        // B_virtual = (number of covering pieces of the last predicate) * B
        double nmax = 0;
        for (double x : n) nmax = std::max(nmax, x);
        e.code = e.inter * nmax / (pieces[k - 1] * (double)B);
        uint64_t ts = *std::min_element(tau.begin(), tau.end());
        if (ts == UINT64_MAX) e.thr = e.inter; // every cell holds all its records: exact
        else {
            // records of R_cap with hash below tau* are all in the intersection; their number over
            // P(hash < tau*) is the bottom-k estimate with the k-th value excluded
            double m = (double)(std::lower_bound(cand.begin(), cand.end(), ts) - cand.begin());
            e.thr = ts > 0 ? m * (double)G63 / (double)ts : 0.0;
        }
        return e;
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
    template <class F> void for_each(F f) const {
        for (uint64_t s = 0; s < size; s++) f(&keys[s * p], (double)cnt[s]);
    }
    uint64_t entries() const { return m; }
    uint64_t held() const { return size; }
};

// =============================================================== CocoSketch
// CocoSketch (Zhang et al., SIGCOMM 2021), software version, following the authors' code
// (github.com/yindazhang/CocoSketch, CPU/Multiple/OurSoft.h): d = 2 arrays of buckets, each holding a
// full key -- here the p-attribute tuple -- and a count. An update hashes the key to one bucket per
// array; if one of them holds the key it is incremented, otherwise the bucket with the smallest count
// is incremented and takes over the key with probability 1/count ("stochastic variance
// minimisation"). A partial-key query -- here any predicate on the tuple -- sums the counts of the
// buckets whose key matches. Charged 32p + 32 bits per bucket.
struct Coco {
    uint32_t p, d;
    uint64_t L, N = 0;
    Rng rng;
    std::vector<uint32_t> keys, cnt; // d*L*p, d*L
    std::vector<uint64_t> hs;
    Coco(uint32_t p_, uint64_t buckets, uint64_t seed, uint32_t d_ = 2) : p(p_), d(d_), rng(seed) {
        L = std::max<uint64_t>(1, buckets / d);
        keys.assign((size_t)d * L * p, 0);
        cnt.assign((size_t)d * L, 0);
        hs.resize(d);
        for (auto& x : hs) x = rng.next();
    }
    static uint64_t bits_per_row(uint32_t p) { return 32ull * p + 32; }
    uint64_t mem_bits() const { return (uint64_t)d * L * bits_per_row(p); }
    inline uint64_t hkey(const uint32_t* r, uint64_t seed) const {
        uint64_t h = seed;
        for (uint32_t a = 0; a < p; a++) h = (h ^ r[a]) * 0x9E3779B97F4A7C15ULL, h ^= h >> 29;
        return mix64(h);
    }
    inline void insert(uint64_t /*rid*/, const uint32_t* rec) {
        N++;
        uint32_t minc = UINT32_MAX;
        size_t minb = 0;
        for (uint32_t i = 0; i < d; i++) {
            size_t b = (size_t)i * L + (size_t)(((unsigned __int128)hkey(rec, hs[i]) * L) >> 64);
            if (cnt[b] > 0 && memcmp(&keys[b * p], rec, 4 * p) == 0) {
                cnt[b]++;
                return;
            }
            if (cnt[b] < minc) minc = cnt[b], minb = b;
        }
        uint32_t c = ++cnt[minb];
        if (rng.below(c) == 0) memcpy(&keys[minb * p], rec, 4 * p);
    }
    void finalize() {}
    double estimate(const Query& q) const {
        uint64_t sum = 0;
        for (size_t b = 0; b < cnt.size(); b++)
            if (cnt[b] && q.matches(&keys[b * p])) sum += cnt[b];
        return (double)sum;
    }
    template <class F> void for_each(F f) const {
        for (size_t b = 0; b < cnt.size(); b++)
            if (cnt[b]) f(&keys[b * p], (double)cnt[b]);
    }
    uint64_t entries() const { return (uint64_t)d * L; }
    uint64_t held() const {
        uint64_t h = 0;
        for (uint32_t c : cnt) h += c > 0;
        return h;
    }
};

// =============================================================== aSH (journal-version baseline)
// Port of the sampling baseline of the OmniSketch journal version (class
// omni.synopses.baselines.aSH.aSH in the authors' omni-deletes repository): adaptive
// Sample-and-Hold of Cohen, Cormode and Duffield (SIGMETRICS 2012) keyed on the *full
// attribute tuple*. An entry holds a count c and a threshold tau; when the table reaches
// `maxSize` entries, it is cut back to `sketchSize` by a batched EjectOne step that draws
// fresh u ~ U(0,1), z ~ Exp(1) per entry, evicts the entries with the smallest
// T = max(tau/u, c/z), and adjusts the survivors (c -= tstar*z if tstar*u > tau; tau = tstar).
// The estimate of a query is the truncated sum of c + tau over the matching entries.
// The port follows the Java control flow; random numbers are not bit-identical.
// Insert-only streams are assumed (the Java class also handles deletions).
struct ASH {
    uint32_t p;
    uint64_t maxSize, sketchSize, N = 0, tmask, ejects = 0;
    Rng rn, rn1;
    std::vector<uint32_t> keys; // maxSize * p
    std::vector<double> c, tau;
    std::vector<int32_t> table;
    uint64_t size = 0, hseed;

    // sampleSize entries in total, of which a fraction `buffer` is the ingest buffer
    ASH(uint32_t p_, uint64_t sampleSize, double buffer, uint64_t seed)
        : p(p_), rn(seed * 31 + 12), rn1(seed * 31 + 10) {
        uint64_t bufferSize = std::max<uint64_t>((uint64_t)(buffer * (double)sampleSize), 1);
        sketchSize = sampleSize > bufferSize ? sampleSize - bufferSize : 1;
        maxSize = sketchSize + bufferSize;
        keys.resize(maxSize * p);
        c.assign(maxSize, 0);
        tau.assign(maxSize, 0);
        uint64_t cap = 16;
        while (cap < 2 * maxSize) cap <<= 1;
        table.assign(cap, -1);
        tmask = cap - 1;
        hseed = rn.next();
    }
    // accounting of the journal code (Formulas.ramASH): 32 p + 5*64 + 2*32 bits per entry
    static uint64_t bits_per_entry_journal(uint32_t p) { return 32ull * p + 5 * 64 + 2 * 32; }
    // lean accounting: the tuple plus two 32-bit floats (c, tau)
    static uint64_t bits_per_entry_lean(uint32_t p) { return 32ull * p + 64; }
    inline uint64_t hkey(const uint32_t* r) const {
        uint64_t h = hseed;
        for (uint32_t a = 0; a < p; a++) h = (h ^ r[a]) * 0x9E3779B97F4A7C15ULL, h ^= h >> 29;
        return mix64(h);
    }
    inline void insert(uint64_t /*rid*/, const uint32_t* rec) {
        N++;
        uint64_t i = hkey(rec) & tmask;
        while (table[i] >= 0) {
            if (memcmp(&keys[(size_t)table[i] * p], rec, 4 * p) == 0) {
                c[table[i]] += 1;
                return;
            }
            i = (i + 1) & tmask;
        }
        uint32_t s = (uint32_t)size++;
        memcpy(&keys[(size_t)s * p], rec, 4 * p);
        c[s] = 1;
        tau[s] = 0;
        table[i] = (int32_t)s;
        if (size >= maxSize) eject();
    }
    void eject() {
        ejects++;
        uint64_t numToEject = size - sketchSize;
        if (numToEject == 0) return;
        std::vector<double> T(size), u(size), z(size), tmp;
        for (uint64_t s = 0; s < size; s++) {
            u[s] = rn.unit();
            z[s] = std::log(rn1.unit() + 1e-300); // = -Exp(1)
            T[s] = std::max(tau[s] / std::max(u[s], 1e-300), c[s] / (-z[s]));
        }
        tmp = T;
        std::nth_element(tmp.begin(), tmp.begin() + (numToEject - 1), tmp.end());
        double tstar = tmp[numToEject - 1];
        uint64_t o = 0;
        for (uint64_t s = 0; s < size; s++) {
            if (T[s] <= tstar) continue;
            double cs = c[s], ts = tau[s];
            if (ts <= tstar) {
                if (tstar * u[s] > ts) cs += tstar * z[s];
                ts = tstar;
            }
            if (o != s) memcpy(&keys[(size_t)o * p], &keys[(size_t)s * p], 4 * p);
            c[o] = cs;
            tau[o] = ts;
            o++;
        }
        size = o;
        std::fill(table.begin(), table.end(), -1);
        for (uint64_t s = 0; s < size; s++) {
            uint64_t i = hkey(&keys[(size_t)s * p]) & tmask;
            while (table[i] >= 0) i = (i + 1) & tmask;
            table[i] = (int32_t)s;
        }
    }
    void finalize() {
        if (size > sketchSize) eject();
    }
    double estimate(const Query& q) const {
        double sum = 0;
        for (uint64_t s = 0; s < size; s++)
            if (q.matches(&keys[s * p])) sum += c[s] + tau[s];
        return std::floor(sum);
    }
};

// =============================================================== low-order marginal tables
// Calibration totals beyond the one-dimensional counter rows: for every attribute pair (and,
// optionally, triple) one table of exact counters indexed by a hash of the value combination,
// i.e. a depth-one Count-Min sketch of each two- (three-) dimensional marginal. A query is
// calibrated to the lightest cell among all cells addressed by subsets of its predicates.
// Updating costs C(p,2) (+ C(p,3)) counter increments per record.
struct MargCal {
    uint32_t p = 0, w2 = 0, w3 = 0, n2 = 0, n3 = 0;
    // Pair tables kept: all C(p,2) pairs by default, or a selected subset (wide streams). pidx maps
    // (a,b), a<b, to the table index or -1; with all pairs, the order (and hence the seeds and the
    // cells) is that of the canonical enumeration a<b, so that results are unchanged.
    std::vector<std::pair<uint32_t, uint32_t>> plist;
    std::vector<int32_t> pidx;
    // Sampled tables: only a hash-selected fraction `beta` of the records updates the tables and
    // counts are scaled by 1/beta. This caps the table work per record for wide streams; the
    // anchors are then estimates from a Bernoulli(beta) sample of the stream instead of exact counts.
    double beta = 1.0, inv_beta = 1.0;
    uint32_t bthr = UINT32_MAX;
    uint64_t bseed = 0;
    std::vector<uint64_t> seed2, seed3;
    std::vector<uint32_t> c2, c3;       // n2*w2, n3*w3 counters
    std::vector<double> W2, V2, W3, V3; // query-time: weight and variance mass of the sample per cell
    std::vector<uint32_t> tri;          // p^3 lookup of triple indexes (only if w3 > 0)

    bool on() const { return w2 > 0 || w3 > 0; }
    static uint64_t npairs(uint32_t p) { return (uint64_t)p * (p - 1) / 2; }
    static uint64_t ntriples(uint32_t p) { return (uint64_t)p * (p - 1) * (p - 2) / 6; }
    // Largest power of two width such that `ntab` tables use at most `bits` bits (0 if < 16).
    static uint32_t width(double bits, uint64_t ntab) {
        if (ntab == 0 || 32.0 * (double)ntab * 16 > bits) return 0;
        uint32_t wd = 16;
        while (wd < (1u << 26) && 32.0 * (double)ntab * (2.0 * wd) <= bits) wd *= 2;
        return wd;
    }
    static uint64_t bits_for(uint32_t p, uint32_t w2, uint32_t w3, uint64_t npair_sel = 0) {
        return 32ull * w2 * (npair_sel ? npair_sel : npairs(p)) + 32ull * w3 * ntriples(p);
    }
    uint64_t mem_bits() const { return 32ull * (c2.size() + c3.size()); }
    // Fraction of records that may update the tables if at most `cap` table increments per record
    // are allowed on average (cap <= 0: no limit).
    static double beta_for_cap(uint32_t p, uint32_t w2, uint32_t w3, double cap) {
        double per = (w2 ? (double)npairs(p) : 0.0) + (w3 ? (double)ntriples(p) : 0.0);
        return cap > 0 && per > cap ? cap / per : 1.0;
    }
    // `sel`: the pairs to keep (nullptr or empty: all pairs)
    void init(uint32_t p_, uint32_t w2_, uint32_t w3_, uint64_t seed, double beta_ = 1.0,
              const std::vector<std::pair<uint32_t, uint32_t>>* sel = nullptr) { // before the first add()
        beta = beta_;
        inv_beta = 1.0 / beta_;
        bthr = beta_ >= 1.0 ? UINT32_MAX : (uint32_t)(beta_ * 4294967296.0);
        bseed = mix64(seed ^ 0x51ED270B0A1D3F27ULL);
        p = p_;
        w2 = w2_;
        w3 = w3_;
        n3 = w3 ? (uint32_t)ntriples(p) : 0;
        Rng r(seed);
        plist.clear();
        pidx.assign((size_t)p * p, -1);
        if (w2) {
            // seeds are drawn for all pairs in canonical order, so a selected pair has the seed it would have
            std::vector<uint64_t> all((size_t)npairs(p));
            for (auto& x : all) x = r.next();
            std::vector<char> keep((size_t)p * p, sel == nullptr || sel->empty());
            if (sel)
                for (auto& ab : *sel) keep[(size_t)std::min(ab.first, ab.second) * p + std::max(ab.first, ab.second)] = 1;
            seed2.clear();
            uint32_t pi = 0;
            for (uint32_t a = 0; a < p; a++)
                for (uint32_t b = a + 1; b < p; b++, pi++)
                    if (keep[(size_t)a * p + b]) {
                        pidx[(size_t)a * p + b] = (int32_t)plist.size();
                        plist.push_back({a, b});
                        seed2.push_back(all[pi]);
                    }
        }
        n2 = (uint32_t)plist.size();
        seed3.resize(n3);
        for (auto& x : seed3) x = r.next();
        c2.assign((size_t)n2 * w2, 0);
        c3.assign((size_t)n3 * w3, 0);
        if (w3) {
            tri.assign((size_t)p * p * p, 0);
            uint32_t t = 0;
            for (uint32_t a = 0; a < p; a++)
                for (uint32_t b = a + 1; b < p; b++)
                    for (uint32_t c = b + 1; c < p; c++) tri[((size_t)a * p + b) * p + c] = t++;
        }
    }
    inline size_t cell2t(uint32_t t, uint32_t va, uint32_t vb) const { // table t
        uint64_t h = mix64((((uint64_t)va << 32) | vb) * 0x9E3779B97F4A7C15ULL + seed2[t]);
        return (size_t)t * w2 + (size_t)(((h >> 32) * (uint64_t)w2) >> 32);
    }
    inline bool has_pair(uint32_t a, uint32_t b) const { return w2 && pidx[(size_t)a * p + b] >= 0; } // a < b
    inline size_t cell2(uint32_t a, uint32_t b, uint32_t va, uint32_t vb) const { // a < b, pair kept
        return cell2t((uint32_t)pidx[(size_t)a * p + b], va, vb);
    }
    inline size_t cell3(uint32_t a, uint32_t b, uint32_t c, uint32_t va, uint32_t vb, uint32_t vc) const { // a < b < c
        uint32_t ti = tri[((size_t)a * p + b) * p + c];
        uint64_t h = mix64((((uint64_t)va << 32) | vb) * 0x9E3779B97F4A7C15ULL + seed3[ti]);
        h = mix64(h ^ ((uint64_t)vc * 0xC2B2AE3D27D4EB4FULL));
        return (size_t)ti * w3 + (size_t)(((h >> 32) * (uint64_t)w3) >> 32);
    }
    inline void add(uint64_t rid, const uint32_t* rec) {
        if (beta < 1.0 && (uint32_t)(mix64(rid ^ bseed) >> 32) >= bthr) return;
        for (uint32_t t = 0; t < n2; t++) c2[cell2t(t, rec[plist[t].first], rec[plist[t].second])]++;
        if (w3)
            for (uint32_t a = 0; a + 2 < p; a++)
                for (uint32_t b = a + 1; b + 1 < p; b++)
                    for (uint32_t c = b + 1; c < p; c++) c3[cell3(a, b, c, rec[a], rec[b], rec[c])]++;
    }
    // deletion of a record added before (the tables are linear: counts of the live records)
    inline void sub(uint64_t rid, const uint32_t* rec) {
        if (beta < 1.0 && (uint32_t)(mix64(rid ^ bseed) >> 32) >= bthr) return;
        for (uint32_t t = 0; t < n2; t++) c2[cell2t(t, rec[plist[t].first], rec[plist[t].second])]--;
        if (w3)
            for (uint32_t a = 0; a + 2 < p; a++)
                for (uint32_t b = a + 1; b + 1 < p; b++)
                    for (uint32_t c = b + 1; c < p; c++) c3[cell3(a, b, c, rec[a], rec[b], rec[c])]--;
    }
    // tables built with the same seeds over a disjoint part of the stream
    void merge(const MargCal& o) {
        for (size_t i = 0; i < c2.size(); i++) c2[i] += o.c2[i];
        for (size_t i = 0; i < c3.size(); i++) c3[i] += o.c3[i];
    }
    void begin_weights() {
        W2.assign(c2.size(), 0.0);
        V2.assign(c2.size(), 0.0);
        W3.assign(c3.size(), 0.0);
        V3.assign(c3.size(), 0.0);
    }
    inline void add_weight(const uint32_t* rec, double wgt, double var) {
        for (uint32_t t = 0; t < n2; t++) {
            size_t c = cell2t(t, rec[plist[t].first], rec[plist[t].second]);
            W2[c] += wgt;
            V2[c] += var;
        }
        if (w3)
            for (uint32_t a = 0; a + 2 < p; a++)
                for (uint32_t b = a + 1; b + 1 < p; b++)
                    for (uint32_t c = b + 1; c < p; c++) {
                        size_t x = cell3(a, b, c, rec[a], rec[b], rec[c]);
                        W3[x] += wgt;
                        V3[x] += var;
                    }
    }
    // Exact count and sample weight of the pair (triple) cell addressed by predicates i, j (, l) of q.
    inline void pair_stat(const Query& q, size_t i, size_t j, double& n, double& W) const {
        if (q.attrs[i] > q.attrs[j]) std::swap(i, j);
        if (!has_pair(q.attrs[i], q.attrs[j])) { n = 0; W = 0; return; } // not kept: no constraint
        size_t c = cell2(q.attrs[i], q.attrs[j], q.vals[i], q.vals[j]);
        n = c2[c] * inv_beta;
        W = W2[c];
    }
    inline void triple_stat(const Query& q, size_t i, size_t j, size_t l, double& n, double& W) const {
        size_t x[3] = {i, j, l};
        std::sort(x, x + 3, [&](size_t a, size_t b) { return q.attrs[a] < q.attrs[b]; });
        size_t c = cell3(q.attrs[x[0]], q.attrs[x[1]], q.attrs[x[2]], q.vals[x[0]], q.vals[x[1]], q.vals[x[2]]);
        n = c3[c] * inv_beta;
        W = W3[c];
    }
    // If a pair/triple cell addressed by the query is lighter than nb, make it the anchor.
    inline void anchor(const Query& q, double& nb, double& Wa, double& Va) const {
        size_t k = q.attrs.size();
        if (!on() || k < 2 || k > 16) return;
        uint32_t idx[16];
        for (size_t i = 0; i < k; i++) idx[i] = (uint32_t)i;
        std::sort(idx, idx + k, [&](uint32_t x, uint32_t y) { return q.attrs[x] < q.attrs[y]; });
        if (w2)
            for (size_t i = 0; i + 1 < k; i++)
                for (size_t j = i + 1; j < k; j++) {
                    if (!has_pair(q.attrs[idx[i]], q.attrs[idx[j]])) continue;
                    size_t c = cell2(q.attrs[idx[i]], q.attrs[idx[j]], q.vals[idx[i]], q.vals[idx[j]]);
                    if (c2[c] * inv_beta < nb) nb = c2[c] * inv_beta, Wa = W2[c], Va = V2[c];
                }
        if (w3 && k >= 3)
            for (size_t i = 0; i + 2 < k; i++)
                for (size_t j = i + 1; j + 1 < k; j++)
                    for (size_t l = j + 1; l < k; l++) {
                        size_t c = cell3(q.attrs[idx[i]], q.attrs[idx[j]], q.attrs[idx[l]], q.vals[idx[i]],
                                         q.vals[idx[j]], q.vals[idx[l]]);
                        if (c3[c] * inv_beta < nb) nb = c3[c] * inv_beta, Wa = W3[c], Va = V3[c];
                    }
    }
};

// =============================================================== raking the query table
// Calibration with *all* anchors of a query instead of the lightest one. T[b], b in [0, 2^k), is
// the sample's weight of the rows whose set of satisfied predicates is exactly b (bit i = predicate
// i). For every subset U of at most three predicates whose cell statistics are available, the
// number of records satisfying U is estimated by the single-anchor estimator n_c * W_U / W_c
// (exact when the cell has no collisions); the table is then fitted to these marginals and to the
// stream size by iterative proportional fitting (Deming-Stephan raking), and the cell "all
// predicates hold" is returned. For k = |U| this is the single-anchor estimate; for larger k it is
// the classical estimator of a contingency table with known lower-order marginals, whose variance
// is governed by the smallest cell of the table rather than by the answer.
// `single(i, n, W)` returns the exact count and sample weight of the C1 bucket of predicate i.
template <class Single>
static double rake_query(const Query& q, std::vector<double>& T, double N, const MargCal& mc, Single single) {
    const size_t k = q.attrs.size(), full = ((size_t)1 << k) - 1;
    if (T[full] <= 0) return 0.0;
    double tot = 0;
    for (double x : T) tot += x;
    if (tot <= 0) return 0.0;
    struct Con {
        size_t mask;
        double target;
    };
    std::vector<Con> cons;
    auto mass = [&](size_t mask) {
        double m = 0;
        for (size_t b = 0; b <= full; b++)
            if ((b & mask) == mask) m += T[b];
        return m;
    };
    auto add = [&](size_t mask, double n, double W) {
        if (W <= 0) return; // no sampled row in the cell: nothing known about its composition
        double t = n * mass(mask) / W;
        cons.push_back({mask, std::min(std::max(t, 0.0), N)});
    };
    for (size_t i = 0; i < k; i++) {
        double n, W;
        single(i, n, W);
        add((size_t)1 << i, n, W);
    }
    if (mc.w2)
        for (size_t i = 0; i + 1 < k; i++)
            for (size_t j = i + 1; j < k; j++) {
                double n, W;
                mc.pair_stat(q, i, j, n, W);
                add(((size_t)1 << i) | ((size_t)1 << j), n, W);
            }
    if (mc.w3)
        for (size_t i = 0; i + 2 < k; i++)
            for (size_t j = i + 1; j + 1 < k; j++)
                for (size_t l = j + 1; l < k; l++) {
                    double n, W;
                    mc.triple_stat(q, i, j, l, n, W);
                    add(((size_t)1 << i) | ((size_t)1 << j) | ((size_t)1 << l), n, W);
                }
    for (double& x : T) x *= N / tot;
    for (int it = 0; it < 200; it++) {
        double dev = 0;
        for (const Con& c : cons) {
            double in = 0;
            for (size_t b = 0; b <= full; b++)
                if ((b & c.mask) == c.mask) in += T[b];
            double out = N - in;
            if (in <= 0 || out <= 0) continue;
            double fi = c.target / in, fo = (N - c.target) / out;
            dev = std::max(dev, std::fabs(fi - 1));
            for (size_t b = 0; b <= full; b++) T[b] *= ((b & c.mask) == c.mask) ? fi : fo;
        }
        if (dev < 1e-7) break;
    }
    return T[full];
}

// =============================================================== dyadic range anchors
// Exact counters of the dyadic intervals of the range attributes, levels 1 .. L_a - 1 (level 0,
// the single values, is the attribute's C1 row). Level l of attribute a counts the records by
// the interval index r_a >> l, directly when there are at most `wr` intervals at that level and
// hashed into wr counters otherwise. A range predicate is anchored to the summed counts of the
// cells of its canonical dyadic cover; the sample's weight in the same cells (with multiplicity,
// as a record may collide into several of them) is computed at query time.
struct DyadCnt {
    uint32_t wr = 0;
    std::vector<int32_t> ridx;              // attribute -> index among the range attributes, or -1
    std::vector<uint32_t> attrs, L;         // range attributes and their numbers of levels
    std::vector<size_t> off;                // per (range attribute, level): first counter
    std::vector<uint32_t> width;            // per (range attribute, level): table width
    std::vector<uint8_t> direct;            // per (range attribute, level)
    std::vector<uint64_t> seed;             // per (range attribute, level)
    std::vector<uint32_t> c;
    bool on() const { return !attrs.empty(); }
    uint64_t mem_bits() const { return 32ull * c.size(); }
    // levels[a] = number of dyadic levels of attribute a (0 or 1: not a range attribute)
    static uint64_t bits_for(const std::vector<uint32_t>& levels, const std::vector<uint32_t>& maxv, uint32_t wr) {
        uint64_t t = 0;
        for (size_t a = 0; a < levels.size(); a++)
            for (uint32_t l = 1; l < levels[a]; l++) t += std::min<uint64_t>(wr, ((uint64_t)maxv[a] >> l) + 1);
        return 32ull * t;
    }
    void init(uint32_t p, const std::vector<uint32_t>& levels, const std::vector<uint32_t>& maxv, uint32_t wr_, uint64_t s0) {
        wr = wr_;
        ridx.assign(p, -1);
        Rng r(s0);
        size_t t = 0;
        for (uint32_t a = 0; a < p; a++) {
            if (levels[a] < 2) continue;
            ridx[a] = (int32_t)attrs.size();
            attrs.push_back(a);
            L.push_back(levels[a]);
            for (uint32_t l = 0; l < levels[a]; l++) {
                uint64_t n = l == 0 ? 0 : std::min<uint64_t>(wr, ((uint64_t)maxv[a] >> l) + 1);
                off.push_back(t);
                width.push_back((uint32_t)n);
                direct.push_back(l > 0 && ((uint64_t)maxv[a] >> l) + 1 <= wr);
                seed.push_back(r.next());
                t += n;
            }
        }
        c.assign(t, 0);
    }
    inline size_t tab(uint32_t i, uint32_t l) const { // index of (range attribute i, level l)
        size_t k = 0;
        for (uint32_t j = 0; j < i; j++) k += L[j];
        return k + l;
    }
    inline size_t cellt(size_t t, uint64_t x) const {
        return off[t] + (direct[t] ? (size_t)std::min<uint64_t>(x, width[t] - 1) : (size_t)bucket_of((uint32_t)x, seed[t], width[t]));
    }
    inline void add(const uint32_t* rec) {
        size_t t = 0;
        for (size_t i = 0; i < attrs.size(); i++) {
            t++; // level 0 is the C1 row
            for (uint32_t l = 1; l < L[i]; l++, t++) c[cellt(t, (uint64_t)rec[attrs[i]] >> l)]++;
        }
    }
    void merge(const DyadCnt& o) {
        for (size_t i = 0; i < c.size(); i++) c[i] += o.c[i];
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
    std::vector<uint32_t> u;                    // M priorities, 32-bit fixed point: (u + 0.5) / 2^32
    MargCal mc;                                 // optional pair/triple tables (calibration only)
    DyadCnt dy;                                 // optional dyadic counters of range attributes (calibration only)
    // query-time working state (built by finalize)
    std::vector<double> invpi, Wb, Vb;
    double Wtot = 0, Vtot = 0;                  // sample weight / variance mass of all rows (anchor N)
    uint32_t nmax_c1 = 0;                       // min over attributes of the largest C1 count
    // Optional query-time index: for every attribute, the stored rows grouped by C1 bucket
    // (CSR layout, p * (w + 1 + size) 32-bit words). A query then scans only the rows of the
    // bucket with the fewest stored rows among those addressed by its predicates.
    bool indexed = false;
    std::vector<uint32_t> ix_start, ix_rows;    // p*(w+1), p*size
    void build_index() { // after finalize()
        ix_start.assign((size_t)p * (w + 1), 0);
        ix_rows.resize((size_t)p * size);
        for (uint32_t a = 0; a < p; a++) {
            uint32_t* st = &ix_start[(size_t)a * (w + 1)];
            for (uint64_t s = 0; s < size; s++) st[bucket_of(col[a][s], hseed[a], w) + 1]++;
            for (uint32_t b = 0; b < w; b++) st[b + 1] += st[b];
            std::vector<uint32_t> pos(st, st + w);
            for (uint64_t s = 0; s < size; s++) ix_rows[(size_t)a * size + pos[bucket_of(col[a][s], hseed[a], w)]++] = (uint32_t)s;
        }
        indexed = true;
    }
    uint64_t index_bits() const { return 32ull * (ix_start.size() + ix_rows.size()); }

    CoStrata(uint32_t p_, uint32_t w_, uint64_t M_, double lf_, double gamma_, uint64_t seed, double slack_ = 0.1)
        : p(p_), w(w_), M(M_), lf(lf_), gamma(gamma_), slack(slack_), hseed(p_), cnt((size_t)p_ * w_, 0), col(p_), u(M_) {
        Rng r(seed);
        useed = r.next();
        for (auto& s : hseed) s = r.next();
        for (auto& c : col) c.resize(M);
    }
    static uint64_t bits_per_row(uint32_t p) { return 32ull * p + 32; }
    static inline double prio(uint32_t x) { return ((double)x + 0.5) * (1.0 / 4294967296.0); }
    uint64_t mem_bits() const { return M * bits_per_row(p) + 32ull * cnt.size() + mc.mem_bits() + dy.mem_bits() + delete_bits(pend_cap); }
    // row_extra: further bits charged per row (e.g. 32 p for the query-time index)
    static uint64_t rows_for_budget(uint64_t bits, uint32_t p, uint32_t w, uint64_t cal_bits = 0, uint64_t row_extra = 0) {
        uint64_t fixed = 32ull * p * w + cal_bits;
        return bits > fixed ? (bits - fixed) / (bits_per_row(p) + row_extra) : 0;
    }
    void enable_marginals(uint32_t w2, uint32_t w3, double beta = 1.0, const std::vector<std::pair<uint32_t, uint32_t>>* sel = nullptr) {
        mc.init(p, w2, w3, useed ^ 0xC2B2AE3D27D4EB4FULL, beta, sel); // before the first insert
    }
    void enable_ranges(const std::vector<uint32_t>& levels, const std::vector<uint32_t>& maxv) {
        dy.init(p, levels, maxv, w, useed ^ 0x2545F4914F6CDD1DULL);
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
        if (mc.on()) mc.add(rid, rec);
        if (dy.on()) dy.add(rec);
        double ns = sfun((double)nmin);
        if (lf > 0) {
            double g = (double)N / (lf * (double)E);
            if (g > G) G = g;
            if (G < ns) ns = G;
        }
        uint32_t ui = (uint32_t)(mix64(rid ^ useed) >> 32);
        double ur = prio(ui);
        if (ur * ns < Bcur) {
            uint64_t s = size++;
            for (uint32_t a = 0; a < p; a++) col[a][s] = rec[a];
            u[s] = ui;
            if (size == M && !pend.empty()) flush_deletes();
            if (size == M) compact();
        }
    }
    // ------------------------------------------------------------ deletions (bounded-deletion streams)
    // Only for the uniform design (house weight lf so large that n* = G for every row; RS, RS+C1/C2/C3):
    // a row is then stored iff u_r G < B, and since G only grows and B only shrinks, a live record
    // passing the test at the end was admitted and never evicted. A deletion decrements N, the C1
    // counters and the pair/triple tables (all linear), and removes the record's row if it is stored,
    // so that the store holds exactly the live records with u_r G < B: a Bernoulli sample of the live
    // data with a known rate, which the anchors calibrate as before. Rows are found in batches: the
    // keys (priority, tuple hash) of deleted records that pass the test are buffered, and one pass over
    // the store removes them when the buffer is full, before a compaction and in finalize(). The
    // buffer (64 bits per key) is charged to the budget. (The stratified design would need thresholds
    // from insertion counts, which deletions do not lower; not implemented.)
    std::vector<uint64_t> pend;
    uint64_t pend_cap = 0, deletes = 0, pend_flushes = 0, rows_deleted = 0;
    void enable_deletes(uint64_t cap) { pend_cap = cap; pend.reserve(cap); }
    static uint64_t delete_bits(uint64_t cap) { return 64ull * cap; }
    inline uint64_t del_key(uint32_t ui, const uint32_t* rec) const {
        uint64_t h = 0x8E2B7F1C3A5D9E47ULL;
        for (uint32_t a = 0; a < p; a++) h = mix64(h ^ ((uint64_t)rec[a] + 0x9E3779B97F4A7C15ULL * (a + 1)));
        return ((uint64_t)ui << 32) | (h >> 32);
    }
    inline void remove(uint64_t rid, const uint32_t* rec) {
        N--;
        deletes++;
        uint32_t nmin = UINT32_MAX;
        for (uint32_t a = 0; a < p; a++) {
            uint32_t& c = cnt[(size_t)a * w + bucket_of(rec[a], hseed[a], w)];
            c--;
            if (c == 0) E--;
            nmin = std::min(nmin, c + 1); // the record's own count included, as when it was stored
        }
        if (mc.on()) mc.sub(rid, rec);
        if (dy.on()) { fprintf(stderr, "CoStrata::remove: dyadic counters do not support deletions\n"); exit(1); }
        uint32_t ui = (uint32_t)(mix64(rid ^ useed) >> 32);
        double ns = sfun((double)nmin);
        if (lf > 0 && G < ns) ns = G;
        if (prio(ui) * ns < Bcur) { // otherwise the row is not stored, or stale (dropped at finalize)
            pend.push_back(del_key(ui, rec));
            if (pend.size() >= pend_cap) flush_deletes();
        }
    }
    void flush_deletes() {
        if (pend.empty()) return;
        pend_flushes++;
        std::sort(pend.begin(), pend.end());
        std::vector<uint32_t> rowbuf(p);
        for (uint64_t s = 0; s < size;) {
            for (uint32_t a = 0; a < p; a++) rowbuf[a] = col[a][s];
            uint64_t k = del_key(u[s], rowbuf.data());
            auto it = std::lower_bound(pend.begin(), pend.end(), k);
            if (it != pend.end() && *it == k) {
                pend.erase(it);
                size--;
                rows_deleted++;
                if (s != size) {
                    for (uint32_t a = 0; a < p; a++) col[a][s] = col[a][size];
                    u[s] = u[size];
                }
                if (pend.empty()) break;
                continue;
            }
            s++;
        }
        pend.clear(); // keys of rows evicted before their deletion
    }
    void compact() {
        compactions++;
        std::vector<double> z(size), zs;
        for (uint64_t s = 0; s < size; s++) z[s] = prio(u[s]) * nstar_slot(s);
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
    // Merge a CoStrata built with the same parameters and seed over a disjoint part of the stream
    // (both before finalize). Counters add, the threshold becomes the smaller one, and a row of
    // either store is kept iff u_r * nstar_r < B under the merged counters. Since nstar_r can only
    // grow and B only shrink, every row passing this test had passed the test of its own part, so the
    // merged store holds exactly the rows of either part that a synopsis with the merged counters and
    // threshold B selects; it is then compacted to capacity as during ingestion.
    void merge(const CoStrata& o) {
        N += o.N;
        for (size_t i = 0; i < cnt.size(); i++) cnt[i] += o.cnt[i];
        E = 0;
        for (uint32_t c : cnt) E += c > 0;
        if (lf > 0) G = std::max(std::max(G, o.G), (double)N / (lf * (double)E));
        Bcur = std::min(Bcur, o.Bcur);
        if (mc.on()) mc.merge(o.mc);
        if (dy.on()) dy.merge(o.dy);
        uint64_t tot = size + o.size;
        if (tot > M) {
            for (auto& c : col) c.resize(tot);
            u.resize(tot);
        }
        for (uint64_t s = 0; s < o.size; s++) {
            for (uint32_t a = 0; a < p; a++) col[a][size + s] = o.col[a][s];
            u[size + s] = o.u[s];
        }
        size = tot;
        uint64_t k = 0;
        for (uint64_t s = 0; s < size; s++) {
            if (!(prio(u[s]) * nstar_slot(s) < Bcur)) continue;
            if (k != s) {
                for (uint32_t a = 0; a < p; a++) col[a][k] = col[a][s];
                u[k] = u[s];
            }
            k++;
        }
        size = k;
        while (size >= M) compact();
        for (auto& c : col) c.resize(M);
        u.resize(M);
    }
    // Drop stale rows, compute Horvitz-Thompson weights and per-stratum aggregates.
    void finalize() {
        flush_deletes();
        uint64_t o = 0;
        invpi.assign(size, 0.0);
        Wb.assign((size_t)p * w, 0.0);
        Vb.assign((size_t)p * w, 0.0);
        Wtot = Vtot = 0;
        nmax_c1 = UINT32_MAX;
        for (uint32_t a = 0; a < p; a++) {
            uint32_t mx = 0;
            for (uint32_t b = 0; b < w; b++) mx = std::max(mx, cnt[(size_t)a * w + b]);
            nmax_c1 = std::min(nmax_c1, mx);
        }
        mc.begin_weights();
        std::vector<uint32_t> rowbuf(p);
        for (uint64_t s = 0; s < size; s++) {
            double ns = nstar_slot(s);
            if (!(prio(u[s]) * ns < Bcur)) continue;
            double pi = std::min(1.0, Bcur / ns);
            if (o != s) {
                for (uint32_t a = 0; a < p; a++) col[a][o] = col[a][s];
                u[o] = u[s];
            }
            invpi[o] = 1.0 / pi;
            Wtot += 1.0 / pi;
            Vtot += (1.0 - pi) / (pi * pi);
            for (uint32_t a = 0; a < p; a++) {
                size_t b = (size_t)a * w + bucket_of(col[a][o], hseed[a], w);
                Wb[b] += 1.0 / pi;
                Vb[b] += (1.0 - pi) / (pi * pi);
            }
            if (mc.on()) {
                for (uint32_t a = 0; a < p; a++) rowbuf[a] = col[a][o];
                mc.add_weight(rowbuf.data(), 1.0 / pi, (1.0 - pi) / (pi * pi));
            }
            o++;
        }
        size = o;
        invpi.resize(size);
    }
    struct Est {
        double hajek, ht, sd_hajek, sd_ht, cm; // cm = Count-Min style upper bound from the counters
        // for confidence intervals: number of matching rows (a lower bound on the answer), and
        // wmax = 1 / pi_q, where pi_q is the smallest inclusion probability that any matching record
        // can have (from the lightest C1 cell of the query); the calibration factor is hajek / ht
        double m = 0, wmax = 0;
        // the anchor cell: sample weight W_a and variance mass V_a of its stored rows, and
        // wanc = 1 / pi_a, where pi_a is the smallest inclusion probability that any record of the
        // anchor cell can have (that of its C1 bucket for a C1 anchor; otherwise the global floor)
        double wa = 0, va = 0, wanc = 0;
    };
    inline bool row_matches(const Query& q, uint64_t s) const {
        for (size_t i = 0; i < q.attrs.size(); i++) {
            uint32_t x = col[q.attrs[i]][s];
            if (q.hi.empty() ? x != q.vals[i] : (x < q.vals[i] || x > q.hi[i])) return false;
        }
        return true;
    }
    // Count, weight and variance terms of the dyadic anchor of range predicate i of q. A stored row in
    // m >= 1 cover cells (m >= 2 only through hash collisions) enters W with weight m / pi; for the
    // linearised variance sum (1/pi^2 - 1/pi) (y - R m)^2 of the ratio estimator, V = sum c m^2 and
    // Vy = sum c y m, with c = 1/pi^2 - 1/pi and y = [row matches q].
    void range_anchor(const Query& q, size_t i, double& n, double& W, double& V, double& Vy) const {
        uint32_t a = q.attrs[i];
        int32_t ri = dy.ridx[a];
        uint32_t La = dy.L[ri];
        uint64_t hi = std::min<uint64_t>(q.hi[i], (1ull << (La - 1)) - 1); // domain [0, 2^(La-1))
        std::vector<std::pair<uint32_t, uint64_t>> cov;
        if (q.vals[i] <= hi) dyadic_cover(q.vals[i], hi, cov);
        size_t t0 = dy.tab((uint32_t)ri, 0);
        std::vector<size_t> pc(cov.size());
        n = 0;
        for (size_t j = 0; j < cov.size(); j++) {
            uint32_t l = cov[j].first;
            if (l >= La) l = La - 1; // cannot happen for hi < 2^La
            if (l == 0) {
                pc[j] = (size_t)a * w + bucket_of((uint32_t)cov[j].second, hseed[a], w);
                n += cnt[pc[j]];
            } else {
                pc[j] = dy.cellt(t0 + l, cov[j].second);
                n += dy.c[pc[j]];
            }
        }
        W = V = Vy = 0;
        const uint32_t* ca = col[a].data();
        for (uint64_t s = 0; s < size; s++) {
            uint32_t m = 0;
            for (size_t j = 0; j < cov.size(); j++) {
                uint32_t l = cov[j].first;
                size_t c = l == 0 ? (size_t)a * w + bucket_of(ca[s], hseed[a], w) : dy.cellt(t0 + l, (uint64_t)ca[s] >> l);
                m += c == pc[j];
            }
            if (m) {
                double ip = invpi[s], c = ip * ip - ip;
                W += m * ip;
                V += (double)m * m * c;
                if (row_matches(q, s)) Vy += m * c;
            }
        }
    }
    // Estimate calibrated to all anchors of the query (see rake_query).
    double estimate_rake(const Query& q) const {
        const size_t k = q.attrs.size();
        if (k > 16 || q.has_range()) return estimate(q).hajek;
        std::vector<double> T((size_t)1 << k, 0.0);
        std::vector<const uint32_t*> cp(k);
        for (size_t i = 0; i < k; i++) cp[i] = col[q.attrs[i]].data();
        for (uint64_t s = 0; s < size; s++) {
            size_t b = 0;
            for (size_t i = 0; i < k; i++) b |= (size_t)(cp[i][s] == q.vals[i]) << i;
            T[b] += invpi[s];
        }
        return rake_query(q, T, (double)N, mc, [&](size_t i, double& n, double& W) {
            size_t b = (size_t)q.attrs[i] * w + bucket_of(q.vals[i], hseed[q.attrs[i]], w);
            n = cnt[b];
            W = Wb[b];
        });
    }
    Est estimate(const Query& q) const {
        const uint32_t* c0 = col[q.attrs[0]].data();
        uint32_t v0 = q.vals[0];
        double num = 0, A = 0, mcount = 0;
        bool any_eq = false;
        for (size_t k = 0; k < q.attrs.size(); k++) any_eq |= q.is_eq(k);
        if (indexed && any_eq) {
            const uint32_t *lo = nullptr, *hi = nullptr;
            for (size_t k = 0; k < q.attrs.size(); k++) {
                if (!q.is_eq(k)) continue;
                uint32_t a = q.attrs[k], b = bucket_of(q.vals[k], hseed[a], w);
                const uint32_t* st = &ix_start[(size_t)a * (w + 1)];
                const uint32_t* base = &ix_rows[(size_t)a * size];
                if (!lo || st[b + 1] - st[b] < (uint32_t)(hi - lo)) lo = base + st[b], hi = base + st[b + 1];
            }
            for (const uint32_t* it = lo; it != hi; ++it) {
                uint32_t s = *it;
                if (row_matches(q, s)) {
                    double ip = invpi[s];
                    num += ip;
                    A += ip * ip - ip;
                    mcount++;
                }
            }
        } else if (q.has_range()) {
            for (uint64_t s = 0; s < size; s++)
                if (row_matches(q, s)) {
                    double ip = invpi[s];
                    num += ip;
                    A += ip * ip - ip;
                    mcount++;
                }
        } else
        for (uint64_t s = 0; s < size; s++) {
            if (c0[s] != v0) continue;
            bool ok = true;
            for (size_t k = 1; k < q.attrs.size(); k++)
                if (col[q.attrs[k]][s] != q.vals[k]) { ok = false; break; }
            if (ok) {
                double ip = invpi[s];
                num += ip;
                A += ip * ip - ip; // (1-pi)/pi^2
                mcount++;
            }
        }
        // lightest stratum touched by the query
        size_t best = SIZE_MAX;
        double nb = 1e300;
        for (size_t k = 0; k < q.attrs.size(); k++) {
            if (!q.is_eq(k)) continue;
            size_t b = (size_t)q.attrs[k] * w + bucket_of(q.vals[k], hseed[q.attrs[k]], w);
            if (cnt[b] < nb) nb = cnt[b], best = b;
        }
        double n1 = nb; // lightest C1 cell: bounds the inclusion probability of matching records
        double Wa = best == SIZE_MAX ? 0.0 : Wb[best], Va = best == SIZE_MAX ? 0.0 : Vb[best];
        bool dyadic = false; // the anchor is a dyadic cover (a multiset of records)
        double Vya = 0;
        if (q.has_range()) {
            mc.anchor(q.eq_part(), nb, Wa, Va);
            if (dy.on())
                for (size_t k = 0; k < q.attrs.size(); k++) {
                    if (q.is_eq(k) || dy.ridx[q.attrs[k]] < 0) continue;
                    double n, W, V, Vy;
                    range_anchor(q, k, n, W, V, Vy);
                    if (n < nb) nb = n, Wa = W, Va = V, Vya = Vy, dyadic = true;
                }
            if (nb >= (double)N) nb = (double)N, Wa = Wtot, Va = Vtot, dyadic = false; // the stream size is an anchor too
        } else
            mc.anchor(q, nb, Wa, Va); // lighter pair/triple cell, if any
        Est e;
        e.ht = num;
        e.sd_ht = std::sqrt(A);
        e.cm = (double)nb;
        e.m = mcount;
        double nsq = sfun(best == SIZE_MAX ? (double)nmax_c1 : n1);
        if (lf > 0 && G < nsq) nsq = G;
        double piq = std::min(1.0, Bcur / nsq);
        e.wmax = 1.0 / piq;
        e.wa = Wa;
        e.va = Va;
        if (dyadic) e.wanc = 0; // no retention floor for the non-matching part of a multiset anchor
        else if (best != SIZE_MAX && nb == n1 && Wa == Wb[best]) e.wanc = e.wmax; // the C1 bucket is the anchor
        else {
            double nsf = sfun((double)nmax_c1); // no record has a larger n*: the floor of all inclusion probabilities
            if (lf > 0 && G < nsf) nsf = G;
            e.wanc = 1.0 / std::min(1.0, Bcur / nsf);
        }
        if (Wa > 0) {
            double R = num / Wa;
            double scale = (double)nb / Wa;
            e.hajek = (double)nb * R;
            double var = dyadic ? scale * scale * (A - 2 * R * Vya + R * R * Va)
                                : scale * scale * ((1 - R) * (1 - R) * A + R * R * std::max(0.0, Va - A));
            e.sd_hajek = std::sqrt(std::max(0.0, var));
        } else {
            e.hajek = 0;
            e.sd_hajek = 0;
        }
        return e;
    }
};

// =============================================================== tuple sketches + calibration
// A tuple-aggregating sketch (USS or CocoSketch) over full tuples plus one row of w exact counters
// per attribute (and optional pair/triple tables); the subset-sum estimate is ratio-calibrated to
// the exact count of the lightest cell addressed by the query (same calibration as CoStrata).
template <class TS>
struct TupleCal {
    TS uss; // the tuple sketch (named after the first instance, USS)
    uint32_t p, w;
    std::vector<uint64_t> hseed;
    std::vector<uint32_t> cnt;
    std::vector<double> Wb;
    MargCal mc;
    uint64_t mseed;
    TupleCal(uint32_t p_, uint32_t w_, uint64_t m, uint64_t seed) : uss(p_, m, seed), p(p_), w(w_), hseed(p_), cnt((size_t)p_ * w_, 0) {
        Rng r(seed ^ 0x5bd1e995);
        for (auto& x : hseed) x = r.next();
        mseed = r.next();
    }
    static uint64_t rows_for_budget(uint64_t bits, uint32_t p, uint32_t w, uint64_t cal_bits = 0) {
        uint64_t fixed = 32ull * p * w + cal_bits;
        return bits > fixed ? (bits - fixed) / TS::bits_per_row(p) : 0;
    }
    void enable_marginals(uint32_t w2, uint32_t w3, double beta = 1.0, const std::vector<std::pair<uint32_t, uint32_t>>* sel = nullptr) {
        mc.init(p, w2, w3, mseed, beta, sel);
    }
    uint64_t mem_bits() const { return uss.mem_bits() + 32ull * cnt.size() + mc.mem_bits(); }
    inline void insert(uint64_t rid, const uint32_t* rec) {
        for (uint32_t a = 0; a < p; a++) cnt[(size_t)a * w + bucket_of(rec[a], hseed[a], w)]++;
        if (mc.on()) mc.add(rid, rec);
        uss.insert(rid, rec);
    }
    void finalize() {
        Wb.assign((size_t)p * w, 0.0);
        mc.begin_weights();
        uss.for_each([&](const uint32_t* key, double c) {
            for (uint32_t a = 0; a < p; a++) Wb[(size_t)a * w + bucket_of(key[a], hseed[a], w)] += c;
            if (mc.on()) mc.add_weight(key, c, 0.0);
        });
    }
    double estimate(const Query& q) const {
        double num = uss.estimate(q);
        size_t best = 0;
        double nb = 1e300;
        for (size_t k = 0; k < q.attrs.size(); k++) {
            if (!q.is_eq(k)) continue;
            size_t b = (size_t)q.attrs[k] * w + bucket_of(q.vals[k], hseed[q.attrs[k]], w);
            if (cnt[b] < nb) nb = cnt[b], best = b;
        }
        if (nb == 1e300) return num; // no equality predicate: the tuple counts already sum to N
        double Wa = Wb[best], Va = 0;
        if (q.has_range()) mc.anchor(q.eq_part(), nb, Wa, Va);
        else mc.anchor(q, nb, Wa, Va);
        return Wa > 0 ? (double)nb * num / Wa : 0.0;
    }
};
using USSC = TupleCal<USS>;
using CocoC = TupleCal<Coco>;

// =============================================================== CoStrata-H (hold variant)
// CoStrata over *distinct tuples* with sample-and-hold counting. The thresholds are those of
// CoStrata (tau_t = min(1, Bcur / nstar_t), nstar_t computed from the counter rows), but the
// store keeps one entry per sampled tuple: the entry is opened by the first row r of the tuple
// with priority u_r < tau_t and from then on counts every later row of the tuple exactly.
// When a threshold has dropped below the priority of an entry's start row, the entry is
// re-anchored at the next row with priority below the new threshold; since the rows counted
// after the start row were never tested, their priorities are still unconditioned, and the
// number of rows to skip is Geometric(tau) (principle of deferred decisions). Re-anchoring is
// done lazily at compactions and before querying. With final threshold tau the entry state is
// that of plain sample-and-hold run at rate tau, so that  count - 1 + 1/tau  estimates the
// multiplicity of the tuple. For streams without repeated tuples this is CoStrata with a
// larger per-entry overhead (count and hash-index slot); with repeated tuples, heavy tuples
// are counted exactly as in (Unbiased) Space-Saving.
struct CoStrataH {
    uint32_t p, w;
    uint64_t M;
    double lf, slack;
    uint64_t N = 0, size = 0, E = 0, compactions = 0, rounds = 0, tmask;
    double Bcur = std::numeric_limits<double>::infinity();
    double G = 0;
    uint64_t useed, tseed;
    Rng rng;
    std::vector<uint64_t> hseed;   // p
    std::vector<uint32_t> cnt;     // p*w
    std::vector<uint32_t> keys;    // M*p
    std::vector<uint32_t> c, u;    // M counts, M start priorities (32-bit fixed point)
    std::vector<int32_t> table;
    std::vector<double> est, evar, Wb, Vb; // query-time working state (built by finalize)
    MargCal mc;                            // optional pair/triple tables (calibration only)

    CoStrataH(uint32_t p_, uint32_t w_, uint64_t M_, double lf_, uint64_t seed, double slack_ = 0.1)
        : p(p_), w(w_), M(M_), lf(lf_), slack(slack_), rng(seed ^ 0x9E3779B97F4A7C15ULL), hseed(p_),
          cnt((size_t)p_ * w_, 0), keys((size_t)M_ * p_), c(M_), u(M_) {
        Rng r(seed);
        useed = r.next();
        tseed = r.next();
        for (auto& s : hseed) s = r.next();
        uint64_t cap = 16;
        while (cap < 2 * M) cap <<= 1;
        table.assign(cap, -1);
        tmask = cap - 1;
    }
    // tuple, count, start priority, hash-index slot
    static uint64_t bits_per_entry(uint32_t p) { return 32ull * p + 96; }
    uint64_t mem_bits() const { return M * bits_per_entry(p) + 32ull * cnt.size() + mc.mem_bits(); }
    static uint64_t rows_for_budget(uint64_t bits, uint32_t p, uint32_t w, uint64_t cal_bits = 0) {
        uint64_t fixed = 32ull * p * w + cal_bits;
        return bits > fixed ? (bits - fixed) / bits_per_entry(p) : 0;
    }
    void enable_marginals(uint32_t w2, uint32_t w3, double beta = 1.0) { mc.init(p, w2, w3, useed ^ 0xC2B2AE3D27D4EB4FULL, beta); }
    static inline double prio(uint32_t x) { return ((double)x + 0.5) * (1.0 / 4294967296.0); }
    inline uint64_t hkey(const uint32_t* r) const {
        uint64_t h = tseed;
        for (uint32_t a = 0; a < p; a++) h = (h ^ r[a]) * 0x9E3779B97F4A7C15ULL, h ^= h >> 29;
        return mix64(h);
    }
    inline double nstar_slot(uint64_t s) const {
        uint32_t nmin = UINT32_MAX;
        const uint32_t* r = &keys[(size_t)s * p];
        for (uint32_t a = 0; a < p; a++) nmin = std::min(nmin, cnt[(size_t)a * w + bucket_of(r[a], hseed[a], w)]);
        double v = (double)nmin;
        return (lf > 0 && G < v) ? G : v;
    }
    inline void insert(uint64_t rid, const uint32_t* rec) {
        N++;
        uint32_t nmin = UINT32_MAX;
        for (uint32_t a = 0; a < p; a++) {
            uint32_t& x = cnt[(size_t)a * w + bucket_of(rec[a], hseed[a], w)];
            if (x == 0) E++;
            x++;
            nmin = std::min(nmin, x);
        }
        if (mc.on()) mc.add(rid, rec);
        double ns = (double)nmin;
        if (lf > 0) {
            double g = (double)N / (lf * (double)E);
            if (g > G) G = g;
            if (G < ns) ns = G;
        }
        uint64_t i = hkey(rec) & tmask;
        while (table[i] >= 0) {
            if (memcmp(&keys[(size_t)table[i] * p], rec, 4 * p) == 0) {
                c[table[i]]++; // hold
                return;
            }
            i = (i + 1) & tmask;
        }
        uint32_t ui = (uint32_t)(mix64(rid ^ useed) >> 32);
        if (prio(ui) * ns < Bcur) {
            uint32_t s = (uint32_t)size++;
            memcpy(&keys[(size_t)s * p], rec, 4 * p);
            c[s] = 1;
            u[s] = ui;
            table[i] = (int32_t)s;
            if (size == M) compact();
        }
    }
    // Re-anchor every entry whose start priority is not below its current threshold.
    void adjust_all() {
        uint64_t o = 0;
        for (uint64_t s = 0; s < size; s++) {
            double ns = nstar_slot(s);
            double tau = std::min(1.0, Bcur / ns);
            uint32_t cs = c[s], us = u[s];
            if (!(prio(us) < tau)) {
                // skip to the next counted row with priority < tau
                double g = 1.0 + std::floor(std::log(rng.unit() + 1e-300) / std::log1p(-tau));
                if (!(g < (double)cs)) continue; // no such row among the counted ones: drop
                cs -= (uint32_t)g;
                double x = std::floor(tau * rng.unit() * 4294967296.0);
                us = (uint32_t)std::min(x, 4294967295.0);
                if (!(prio(us) < tau)) us = us > 0 ? us - 1 : 0;
            }
            if (o != s) memcpy(&keys[(size_t)o * p], &keys[(size_t)s * p], 4 * p);
            c[o] = cs;
            u[o] = us;
            o++;
        }
        size = o;
        std::fill(table.begin(), table.end(), -1);
        for (uint64_t s = 0; s < size; s++) {
            uint64_t i = hkey(&keys[(size_t)s * p]) & tmask;
            while (table[i] >= 0) i = (i + 1) & tmask;
            table[i] = (int32_t)s;
        }
    }
    void compact() {
        compactions++;
        uint64_t K = (uint64_t)((1.0 - slack) * (double)M);
        adjust_all();
        std::vector<double> z;
        while (size > K) {
            rounds++;
            z.resize(size);
            for (uint64_t s = 0; s < size; s++) z[s] = prio(u[s]) * nstar_slot(s);
            std::nth_element(z.begin(), z.begin() + K, z.end());
            if (z[K] < Bcur) Bcur = z[K];
            else Bcur = std::nextafter(Bcur, 0.0);
            adjust_all();
        }
    }
    void finalize() {
        adjust_all();
        est.assign(size, 0.0);
        evar.assign(size, 0.0);
        Wb.assign((size_t)p * w, 0.0);
        Vb.assign((size_t)p * w, 0.0);
        mc.begin_weights();
        for (uint64_t s = 0; s < size; s++) {
            double ns = nstar_slot(s);
            double tau = std::min(1.0, Bcur / ns);
            double e = tau < 1.0 ? (double)c[s] - 1.0 + 1.0 / tau : (double)c[s];
            double v = tau < 1.0 ? (1.0 - tau) / (tau * tau) : 0.0;
            est[s] = e;
            evar[s] = v;
            const uint32_t* r = &keys[(size_t)s * p];
            for (uint32_t a = 0; a < p; a++) {
                size_t b = (size_t)a * w + bucket_of(r[a], hseed[a], w);
                Wb[b] += e;
                Vb[b] += v;
            }
            if (mc.on()) mc.add_weight(r, e, v);
        }
    }
    struct Est {
        double hajek, ht, sd_hajek, cm;
        double m = 0, wmax = 0; // held count of the matching entries (a lower bound); see CoStrata::Est
        double wa = 0, va = 0, wanc = 0;
    };
    uint32_t nmax_c1() const { // min over attributes of the largest C1 count
        uint32_t r = UINT32_MAX;
        for (uint32_t a = 0; a < p; a++) {
            uint32_t mx = 0;
            for (uint32_t b = 0; b < w; b++) mx = std::max(mx, cnt[(size_t)a * w + b]);
            r = std::min(r, mx);
        }
        return r;
    }
    Est estimate(const Query& q) const {
        double num = 0, A = 0, held = 0;
        for (uint64_t s = 0; s < size; s++) {
            if (!q.matches(&keys[(size_t)s * p])) continue;
            num += est[s];
            A += evar[s];
            held += c[s];
        }
        size_t best = SIZE_MAX;
        double nb = 1e300;
        for (size_t k = 0; k < q.attrs.size(); k++) {
            if (!q.is_eq(k)) continue;
            size_t b = (size_t)q.attrs[k] * w + bucket_of(q.vals[k], hseed[q.attrs[k]], w);
            if (cnt[b] < nb) nb = cnt[b], best = b;
        }
        double n1 = nb;
        double Wa = 0, Va = 0;
        if (best != SIZE_MAX) Wa = Wb[best], Va = Vb[best];
        if (q.has_range()) mc.anchor(q.eq_part(), nb, Wa, Va);
        else mc.anchor(q, nb, Wa, Va);
        if (nb >= (double)N) {
            nb = (double)N;
            Wa = Va = 0;
            for (uint64_t s = 0; s < size; s++) Wa += est[s], Va += evar[s];
        }
        Est e;
        e.ht = num;
        e.cm = (double)nb;
        e.m = held;
        double nsq = best == SIZE_MAX ? (double)N : n1;
        if (lf > 0 && G < nsq) nsq = G;
        double tq = std::min(1.0, Bcur / nsq);
        e.wmax = 1.0 / tq;
        e.wa = Wa;
        e.va = Va;
        if (best != SIZE_MAX && nb == n1 && Wa == Wb[best]) e.wanc = e.wmax;
        else {
            double nsf = (double)nmax_c1();
            if (lf > 0 && G < nsf) nsf = G;
            e.wanc = 1.0 / std::min(1.0, Bcur / nsf);
        }
        if (Wa > 0) {
            double R = num / Wa;
            double scale = (double)nb / Wa;
            e.hajek = (double)nb * R;
            double var = scale * scale * ((1 - R) * (1 - R) * A + R * R * std::max(0.0, Va - A));
            e.sd_hajek = std::sqrt(std::max(0.0, var));
        } else {
            e.hajek = 0;
            e.sd_hajek = 0;
        }
        return e;
    }
};
