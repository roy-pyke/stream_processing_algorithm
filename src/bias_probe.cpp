// Unbiasedness probe for the threshold samplers on a small synthetic stream.
//   bias_probe <variant> <N> <p> <M> <w> <lf> <seeds> [domain] [alpha]
// For every seed (hash functions and priorities) the stream is summarised and the
// Horvitz-Thompson estimate of (a) the stream size N and (b) the count of a fixed frequent
// value and of a fixed rare-ish value of attribute 0 is compared with the truth. Prints the
// mean relative deviation and its standard error over the seeds. A small stream and a small
// store make compactions frequent and allow many thousands of seeds.
//   variant: cs (CoStrata), csh (CoStrata-H)
#include "synopses.hpp"

int main(int argc, char** argv) {
    if (argc < 8) return 1;
    std::string variant = argv[1];
    uint64_t N = strtoull(argv[2], 0, 10);
    uint32_t p = (uint32_t)atoi(argv[3]);
    uint64_t M = strtoull(argv[4], 0, 10);
    uint32_t w = (uint32_t)atoi(argv[5]);
    double lf = atof(argv[6]);
    int seeds = atoi(argv[7]);
    uint32_t dom = argc > 8 ? (uint32_t)atoi(argv[8]) : 1000;
    double alpha = argc > 9 ? atof(argv[9]) : 1.0;
    // fixed stream: Zipf(alpha) over `dom` values per attribute, independent attributes
    std::vector<double> cdf(dom);
    double z = 0;
    for (uint32_t k = 0; k < dom; k++) z += std::pow(k + 1.0, -alpha), cdf[k] = z;
    Rng g(12345);
    std::vector<uint32_t> data(N * p);
    for (auto& x : data) {
        double u = g.unit() * z;
        x = (uint32_t)(std::lower_bound(cdf.begin(), cdf.end(), u) - cdf.begin());
    }
    Query qa, qb;
    qa.attrs = {0}; qa.vals = {0};
    qb.attrs = {0}; qb.vals = {dom / 10};
    double ta = 0, tb = 0;
    for (uint64_t i = 0; i < N; i++) ta += data[i * p] == 0, tb += data[i * p] == dom / 10;
    double s[3] = {0, 0, 0}, s2[3] = {0, 0, 0};
    for (int sd = 1; sd <= seeds; sd++) {
        double e[3];
        if (variant == "cs") {
            CoStrata c(p, w, M, lf, 1.0, (uint64_t)sd * 7919 + 4);
            for (uint64_t i = 0; i < N; i++) c.insert(i, &data[i * p]);
            c.finalize();
            double tot = 0;
            for (double x : c.invpi) tot += x;
            e[0] = tot / (double)N - 1;
            e[1] = c.estimate(qa).ht / ta - 1;
            e[2] = c.estimate(qb).ht / tb - 1;
        } else {
            CoStrataH c(p, w, M, lf, (uint64_t)sd * 7919 + 6);
            for (uint64_t i = 0; i < N; i++) c.insert(i, &data[i * p]);
            c.finalize();
            double tot = 0;
            for (double x : c.est) tot += x;
            e[0] = tot / (double)N - 1;
            e[1] = c.estimate(qa).ht / ta - 1;
            e[2] = c.estimate(qb).ht / tb - 1;
        }
        for (int j = 0; j < 3; j++) s[j] += e[j], s2[j] += e[j] * e[j];
    }
    const char* nm[3] = {"total N", "frequent value", "rarer value"};
    printf("%s N=%llu p=%u M=%llu w=%u lf=%g seeds=%d (truths: %.0f, %.0f)\n", variant.c_str(), (unsigned long long)N, p,
           (unsigned long long)M, w, lf, seeds, ta, tb);
    for (int j = 0; j < 3; j++) {
        double mean = s[j] / seeds, var = s2[j] / seeds - mean * mean;
        printf("  %-15s mean rel. deviation %+.5f%%  +- %.5f%%  (z = %+.2f)\n", nm[j], 100 * mean, 100 * std::sqrt(var / seeds),
               mean / std::sqrt(var / seeds));
    }
    return 0;
}
