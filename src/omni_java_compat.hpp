// Bit-exact C++ port of the authors' Java OmniSketch (github.com/wiegerpunter/omnisketch,
// commit 8dc67b5; classes omniVLDB.omni.{OmniSketch,CountMin,Kmin}, omniVLDB.{RamToPar,DetermineB}).
//
// Purpose: validation only. Fed the same records and queries, this port must return the
// same integers as the Java code (scripts/validate_omni_java.py checks this per query).
// It deliberately reproduces every implementation detail that influences the estimates:
//   * hash functions built from java.util.Random re-seeded per call
//       cell   h_j(v) = Random(v + repetition).nextInt(width), j-th draw     (same for all attributes)
//       record g(id)  = Random(id + repetition).nextLong(min((int)2^b, Integer.MAX_VALUE))
//   * b = ceil(ln(4 B^2.5 / delta))  (natural logarithm), so g has at most 2^b distinct values
//   * Kmin.add(): a TreeSet of *distinct* hash values; `curSampleSize` is incremented even
//     when the value was already present, and a duplicate arriving at a full cell evicts
//     the maximum without adding anything, so a cell can hold fewer than B values for good
//   * the estimator of Theorem 3.3 (query()), which returns the bound
//     2 n_max ln(4 p d sqrt(B)/delta) / (B eps^2) when |S_cap| < 3 ln(4 p d sqrt(B)/delta)/eps^2,
//     and the always-scaled value ceil(|S_cap| n_max / B) that the code logs as "case2Estimate"
//   * the memory accounting B (b + 3*32 + 1) + 32 bits per cell and the search for B
#pragma once
#include "common.hpp"
#include <set>

struct JavaRandom {
    uint64_t s = 0;
    static constexpr uint64_t MASK = (1ULL << 48) - 1;
    void setSeed(int64_t seed) { s = ((uint64_t)seed ^ 0x5DEECE66DULL) & MASK; }
    int32_t next(int bits) {
        s = (s * 0x5DEECE66DULL + 0xBULL) & MASK;
        return (int32_t)(uint32_t)(s >> (48 - bits));
    }
    int32_t nextInt(int32_t bound) {
        int32_t r = next(31);
        int32_t m = bound - 1;
        if ((bound & m) == 0) {
            r = (int32_t)(((int64_t)bound * (int64_t)r) >> 31);
        } else {
            for (int32_t u = r;; u = next(31)) {
                r = u % bound;
                if ((int32_t)((uint32_t)u - (uint32_t)r + (uint32_t)m) >= 0) break;
            }
        }
        return r;
    }
    int64_t nextLong() {
        int64_t hi = (int64_t)next(32);
        int64_t lo = (int64_t)next(32);
        return (int64_t)(((uint64_t)hi << 32) + (uint64_t)lo);
    }
    // RandomGenerator.nextLong(bound) default (jdk.internal.util.random.RandomSupport.boundedNextLong)
    int64_t nextLong(int64_t bound) {
        const int64_t m = bound - 1;
        int64_t r = nextLong();
        if ((bound & m) == 0) {
            r &= m;
        } else {
            for (int64_t u = (int64_t)((uint64_t)r >> 1);; u = (int64_t)((uint64_t)nextLong() >> 1)) {
                r = u % bound;
                if ((int64_t)((uint64_t)u + (uint64_t)m - (uint64_t)r) >= 0) break;
            }
        }
        return r;
    }
};

struct OmniJava {
    struct Cell {
        std::set<int64_t> sketch; // distinct hash values (Java: TreeSet in reverse order)
        int32_t n = 0, curSampleSize = 0;
        int64_t curTreeRoot = INT64_MAX;
    };
    uint32_t p;
    int32_t depth, width, B, b, repetition;
    double eps, delta;
    int64_t maxHash;
    std::vector<Cell> cells; // p * depth * width

    // ---- parameter selection: RamToPar.getInfoViaB + DetermineB
    static double compRam(int32_t depth, int32_t width, uint32_t p, double delta, int32_t B) {
        int32_t smallb = (int32_t)std::ceil(std::log(4 * std::pow((double)B, 2.5) / delta));
        return (double)depth * width * (double)(B * (smallb + 3 * 32 + 1) + 32) * p;
    }
    static int32_t searchB(int64_t ram, int32_t depth, int32_t width, uint32_t p, double delta) {
        int32_t prevB = 0, Bc = 2, highB = INT32_MAX;
        for (int iters = 0;;) {
            if (Bc <= 1) return Bc;
            if (iters > 100) return compRam(depth, width, p, delta, prevB) < (double)ram ? prevB : -1;
            iters++;
            double used = compRam(depth, width, p, delta, Bc);
            if ((double)ram * 0.99 <= used && used <= (double)ram) return Bc;
            if (used < (double)ram * 0.99) {
                int32_t nb = highB == INT32_MAX ? Bc * 2 : std::min(highB, Bc * 2);
                prevB = Bc;
                Bc = nb;
            } else {
                highB = Bc;
                Bc = (int32_t)std::ceil((double)(Bc + prevB) / 2);
            }
        }
    }

    OmniJava(uint32_t p_, int64_t ram_bits, int32_t repetition_, double eps_ = 0.1, double delta_ = 0.1)
        : p(p_), repetition(repetition_), eps(eps_), delta(delta_) {
        double deltaCM = delta / 2;
        depth = (int32_t)std::ceil(std::log(1 / deltaCM) / std::log(std::exp(1.0)));
        double epsCM = std::pow(eps / (1 + eps), 1.0 / depth);
        width = 1 + (int32_t)std::ceil(std::exp(1.0) / epsCM);
        B = searchB(ram_bits, depth, width, p, delta);
        b = (int32_t)std::ceil(std::log(4 * std::pow((double)B, 2.5) / delta));
        double ph = std::pow(2.0, b);
        int32_t asInt = ph >= 2147483647.0 ? INT32_MAX : (int32_t)ph; // Java (int) cast saturates
        maxHash = std::min<int64_t>(asInt, INT32_MAX);
        cells.resize((size_t)p * depth * width);
    }
    uint64_t paper_bits_nominal() const { return (uint64_t)compRam(depth, width, p, delta, B); }
    uint64_t mem_bits() const { // OmniSketch.getMemoryUsage()
        uint64_t m = 0;
        for (auto& c : cells) m += (uint64_t)c.sketch.size() * (uint64_t)(b + 3 * 32 + 1) + 32;
        return m;
    }
    inline int64_t ghash(uint64_t rid) const {
        JavaRandom rn;
        rn.setSeed((int64_t)((int32_t)(uint32_t)rid + repetition));
        return rn.nextLong(maxHash);
    }
    inline void insert(uint64_t rid, const uint32_t* rec) {
        int64_t hx = ghash(rid);
        JavaRandom rn;
        for (uint32_t a = 0; a < p; a++) {
            rn.setSeed((int64_t)rec[a] + repetition);
            for (int32_t j = 0; j < depth; j++) {
                int32_t k = rn.nextInt(width);
                Cell& c = cells[((size_t)a * depth + j) * width + k];
                c.n++;
                if (c.curSampleSize < B) {
                    c.sketch.insert(hx);
                    c.curSampleSize++;
                } else if (hx < c.curTreeRoot) {
                    c.sketch.erase(std::prev(c.sketch.end()));
                    c.sketch.insert(hx);
                    c.curTreeRoot = *c.sketch.rbegin();
                }
            }
        }
    }
    void finalize() {}
    struct Est {
        double thm33, case2, inter, nmax;
        bool case2_branch;
    };
    Est estimate(const Query& q) const {
        std::vector<const Cell*> cs;
        int32_t nmax = 0;
        JavaRandom rn;
        for (size_t i = 0; i < q.attrs.size(); i++) {
            rn.setSeed((int64_t)q.vals[i] + repetition);
            for (int32_t j = 0; j < depth; j++) {
                int32_t k = rn.nextInt(width);
                const Cell* c = &cells[((size_t)q.attrs[i] * depth + j) * width + k];
                cs.push_back(c);
                nmax = std::max(nmax, c->n);
            }
        }
        size_t best = 0;
        for (size_t i = 1; i < cs.size(); i++)
            if (cs[i]->sketch.size() < cs[best]->sketch.size()) best = i;
        int64_t S = 0;
        for (int64_t h : cs[best]->sketch) {
            bool ok = true;
            for (size_t i = 0; i < cs.size() && ok; i++)
                if (i != best && !cs[i]->sketch.count(h)) ok = false;
            S += ok;
        }
        int32_t pq = (int32_t)q.attrs.size();
        double L = std::log((4 * pq * depth * std::sqrt((double)B)) / delta);
        double constraint = 3 * L / (eps * eps);
        Est e;
        e.inter = (double)S;
        e.nmax = nmax;
        e.case2 = std::ceil((double)S * nmax / B);
        e.case2_branch = !((double)S < constraint);
        e.thm33 = e.case2_branch ? e.case2 : std::ceil(2 * nmax * L / (B * eps * eps));
        return e;
    }
};
