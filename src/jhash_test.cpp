// Prints test vectors of the ported java.util.Random-based hash functions in the same
// format as `java omniVLDB.Harness hashes <rep> <b> <width> <depth>`; the two outputs must
// be identical (checked by scripts/validate_omni_java.sh).
#include "omni_java_compat.hpp"

int main(int argc, char** argv) {
    if (argc < 5) return 1;
    int rep = atoi(argv[1]), b = atoi(argv[2]), width = atoi(argv[3]), depth = atoi(argv[4]);
    double ph = std::pow(2.0, b);
    int64_t maxHash = ph >= 2147483647.0 ? INT32_MAX : (int32_t)ph;
    long long ids[] = {0, 1, 2, 3, 12345, 4999999, 11000000, 2147483000LL};
    for (long long id : ids) {
        JavaRandom r;
        r.setSeed((int64_t)((int32_t)id + rep));
        printf("G %lld %lld\n", id, (long long)r.nextLong(maxHash));
    }
    long long vals[] = {0, 1, 2, 9999, 65535, 3232235777LL, 4294967295LL};
    for (long long v : vals) {
        JavaRandom r;
        r.setSeed(v + rep);
        printf("H %lld", v);
        for (int j = 0; j < depth; j++) printf(" %d", r.nextInt(width));
        printf("\n");
    }
    return 0;
}
