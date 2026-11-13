// Common utilities: hashing, dataset I/O, query representation.
#pragma once
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <string>
#include <vector>

// ---------------------------------------------------------------- hashing
static inline uint64_t mix64(uint64_t x) {
    x ^= x >> 33;
    x *= 0xff51afd7ed558ccdULL;
    x ^= x >> 33;
    x *= 0xc4ceb9fe1a85ec53ULL;
    x ^= x >> 33;
    return x;
}
// Hash a 32-bit attribute value into [0, w).
static inline uint32_t bucket_of(uint32_t v, uint64_t seed, uint32_t w) {
    uint64_t h = mix64((uint64_t)v * 0x9E3779B97F4A7C15ULL + seed);
    return (uint32_t)(((h >> 32) * (uint64_t)w) >> 32);
}
// Uniform double in [0,1) from a 64-bit hash.
static inline double to_unit(uint64_t h) { return (double)(h >> 11) * (1.0 / 9007199254740992.0); }

struct Rng { // splitmix64
    uint64_t s;
    explicit Rng(uint64_t seed) : s(seed) {}
    uint64_t next() {
        uint64_t z = (s += 0x9E3779B97F4A7C15ULL);
        z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
        z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
        return z ^ (z >> 31);
    }
    double unit() { return to_unit(next()); }
    uint64_t below(uint64_t n) { return (uint64_t)(((unsigned __int128)next() * n) >> 64); }
};

// ---------------------------------------------------------------- dataset
// Binary layout: "STRM" | u32 p | u64 n | n*p u32 (row-major).
struct Dataset {
    uint64_t n = 0;
    uint32_t p = 0;
    std::vector<uint32_t> v;
    const uint32_t* row(uint64_t i) const { return v.data() + i * p; }
    bool load(const std::string& path) {
        FILE* f = fopen(path.c_str(), "rb");
        if (!f) return false;
        char magic[4];
        if (fread(magic, 1, 4, f) != 4 || memcmp(magic, "STRM", 4) != 0) { fclose(f); return false; }
        if (fread(&p, 4, 1, f) != 1 || fread(&n, 8, 1, f) != 1) { fclose(f); return false; }
        v.resize(n * p);
        size_t got = fread(v.data(), 4, n * p, f);
        fclose(f);
        return got == n * p;
    }
    bool save(const std::string& path) const {
        FILE* f = fopen(path.c_str(), "wb");
        if (!f) return false;
        fwrite("STRM", 1, 4, f);
        fwrite(&p, 4, 1, f);
        fwrite(&n, 8, 1, f);
        fwrite(v.data(), 4, n * p, f);
        fclose(f);
        return true;
    }
    // Keep only the first k attributes (used for the attribute-scaling experiment).
    void project_prefix(uint32_t k) {
        if (k >= p) return;
        std::vector<uint32_t> nv(n * k);
        for (uint64_t i = 0; i < n; i++)
            for (uint32_t j = 0; j < k; j++) nv[i * k + j] = v[i * p + j];
        v.swap(nv);
        p = k;
    }
};

// ---------------------------------------------------------------- queries
// Conjunction of equality predicates: attrs[i] == vals[i].
struct Query {
    std::vector<uint32_t> attrs, vals;
    uint64_t truth = 0;
    uint64_t minfreq = 0; // frequency of the rarest single predicate
    int kind = 0; // workload family
    bool matches(const uint32_t* rec) const {
        for (size_t i = 0; i < attrs.size(); i++)
            if (rec[attrs[i]] != vals[i]) return false;
        return true;
    }
};

// Text format, one query per line: kind npred truth minfreq a1 v1 a2 v2 ...
static inline bool load_workload(const std::string& path, std::vector<Query>& qs) {
    FILE* f = fopen(path.c_str(), "r");
    if (!f) return false;
    int kind, np;
    unsigned long long truth, minfreq;
    while (fscanf(f, "%d %d %llu %llu", &kind, &np, &truth, &minfreq) == 4) {
        Query q;
        q.kind = kind;
        q.truth = truth;
        q.minfreq = minfreq;
        for (int i = 0; i < np; i++) {
            unsigned a, v;
            if (fscanf(f, "%u %u", &a, &v) != 2) { fclose(f); return false; }
            q.attrs.push_back(a);
            q.vals.push_back(v);
        }
        qs.push_back(std::move(q));
    }
    fclose(f);
    return true;
}

static inline double now_sec() {
    return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count();
}
