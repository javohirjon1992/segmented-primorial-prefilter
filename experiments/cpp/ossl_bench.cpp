// Comparison with OpenSSL: prime generation time and distribution of residues.
//   ossl_bench gen <bits> <K>     -> bits,ms
//   ossl_bench resid <bits> <N>   -> m3,m5,m7,m11,m13
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <string>
extern "C" {
typedef struct bignum_st BIGNUM;
BIGNUM* BN_new(void);
void BN_free(BIGNUM*);
int BN_generate_prime_ex(BIGNUM*, int, int, const BIGNUM*, const BIGNUM*, void*);
unsigned long long BN_mod_word(const BIGNUM*, unsigned long long);
const char* OpenSSL_version(int);
}
int main(int argc, char** argv) {
    if (argc < 4) { std::fprintf(stderr, "ossl_bench gen|resid bits count\n"); return 1; }
    std::string cmd = argv[1]; int bits = std::atoi(argv[2]), n = std::atoi(argv[3]);
    std::fprintf(stderr, "%s\n", OpenSSL_version(0));
    BIGNUM* p = BN_new();
    if (cmd == "gen") {
        std::printf("bits,ms\n");
        for (int i = 0; i < n; ++i) {
            auto a = std::chrono::steady_clock::now();
            if (BN_generate_prime_ex(p, bits, 0, nullptr, nullptr, nullptr) != 1) return 2;
            double ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - a).count();
            std::printf("%d,%.3f\n", bits, ms); std::fflush(stdout);
        }
    } else {
        std::printf("m3,m5,m7,m11,m13\n");
        for (int i = 0; i < n; ++i) {
            if (BN_generate_prime_ex(p, bits, 0, nullptr, nullptr, nullptr) != 1) return 2;
            std::printf("%llu,%llu,%llu,%llu,%llu\n", BN_mod_word(p, 3), BN_mod_word(p, 5), BN_mod_word(p, 7), BN_mod_word(p, 11), BN_mod_word(p, 13));
        }
    }
    BN_free(p);
    return 0;
}
