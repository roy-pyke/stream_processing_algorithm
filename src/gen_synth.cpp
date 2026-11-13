// Synthetic stream generator.
//   gen_synth <out.bin> <n> <p> <alpha> <domain> <seed>
// Every attribute is drawn independently from Zipf(alpha) over {1..domain}
// (alpha = 0 gives the uniform distribution). This mirrors the synthetic
// streams of the OmniSketch paper (p = 5, domain = 10000, n = 5M).
#include "common.hpp"

int main(int argc, char** argv) {
    if (argc < 7) {
        fprintf(stderr, "usage: gen_synth out n p alpha domain seed\n");
        return 1;
    }
    Dataset d;
    d.n = strtoull(argv[2], 0, 10);
    d.p = (uint32_t)atoi(argv[3]);
    double alpha = atof(argv[4]);
    uint32_t dom = (uint32_t)atoi(argv[5]);
    Rng rng(strtoull(argv[6], 0, 10));
    std::vector<double> cdf(dom);
    double s = 0;
    for (uint32_t k = 1; k <= dom; k++) s += std::pow((double)k, -alpha), cdf[k - 1] = s;
    for (auto& c : cdf) c /= s;
    // Per-attribute random relabeling so that value ids carry no rank information.
    std::vector<std::vector<uint32_t>> perm(d.p, std::vector<uint32_t>(dom));
    for (uint32_t a = 0; a < d.p; a++) {
        for (uint32_t k = 0; k < dom; k++) perm[a][k] = k + 1;
        for (uint32_t k = dom - 1; k > 0; k--) std::swap(perm[a][k], perm[a][rng.below(k + 1)]);
    }
    d.v.resize(d.n * d.p);
    for (uint64_t i = 0; i < d.n * d.p; i++) {
        double x = rng.unit();
        uint32_t k = (uint32_t)(std::lower_bound(cdf.begin(), cdf.end(), x) - cdf.begin());
        if (k >= dom) k = dom - 1;
        d.v[i] = perm[i % d.p][k];
    }
    if (!d.save(argv[1])) {
        fprintf(stderr, "cannot write %s\n", argv[1]);
        return 1;
    }
    fprintf(stderr, "wrote %s: n=%llu p=%u alpha=%.2f domain=%u\n", argv[1], (unsigned long long)d.n, d.p, alpha, dom);
    return 0;
}
