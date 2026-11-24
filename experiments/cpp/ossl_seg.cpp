// OpenSSL's prime generation algorithm and variants with a segmented primorial prefilter (V0-V3).
// Everything uses OpenSSL BN arithmetic only (no GMP).
//   ossl_seg <bits> <K> <B1,B2,...>   -> CSV
//   ossl_seg <bits> <K> <B> confirm4 [warmup=20]   -> confirmatory V0-V3 experiment at 2048/3072 bits
//        (docs/protocol_V0-V3_2048_3072_EN.txt): in every block the four variants run in the order of row (k mod 4) of a
//        Williams Latin square; the first `warmup` blocks are not recorded.
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
extern "C" {
typedef struct bignum_st BIGNUM; typedef struct bignum_ctx BN_CTX; typedef struct bn_gencb_st BN_GENCB;
BN_CTX* BN_CTX_new(void); void BN_CTX_free(BN_CTX*);
BIGNUM* BN_new(void); void BN_free(BIGNUM*); BIGNUM* BN_copy(BIGNUM*, const BIGNUM*);
int BN_priv_rand_ex(BIGNUM*, int, int, int, unsigned int, BN_CTX*);
unsigned long long BN_mod_word(const BIGNUM*, unsigned long long);
int BN_add_word(BIGNUM*, unsigned long long); int BN_mul_word(BIGNUM*, unsigned long long); int BN_set_word(BIGNUM*, unsigned long long);
int BN_num_bits(const BIGNUM*); int BN_gcd(BIGNUM*, const BIGNUM*, const BIGNUM*, BN_CTX*); int BN_is_one(const BIGNUM*);
int BN_mul(BIGNUM*, const BIGNUM*, const BIGNUM*, BN_CTX*); int BN_mod(BIGNUM*, const BIGNUM*, const BIGNUM*, BN_CTX*) ;
int BN_div(BIGNUM*, BIGNUM*, const BIGNUM*, const BIGNUM*, BN_CTX*);
int BN_is_prime_fasttest_ex(const BIGNUM*, int, BN_CTX*, int, BN_GENCB*);
int BN_generate_prime_ex(BIGNUM*, int, int, const BIGNUM*, const BIGNUM*, BN_GENCB*);
BN_GENCB* BN_GENCB_new(void); void BN_GENCB_free(BN_GENCB*); void BN_GENCB_set(BN_GENCB*, int (*)(int, int, BN_GENCB*), void*);
const char* OpenSSL_version(int);
}
#define BN_mod(r, a, m, ctx) BN_div(NULL, (r), (a), (m), (ctx))
static double now() { return std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now().time_since_epoch()).count(); }
static std::vector<uint64_t> primes_upto(uint64_t n) {          // starting from 2 (as in OpenSSL's primes[] table)
    std::vector<char> s(n + 1, 1); std::vector<uint64_t> p; s[0] = s[1] = 0;
    for (uint64_t i = 2; i <= n; ++i) if (s[i]) { p.push_back(i); for (uint64_t j = i * i; j <= n; j += i) s[j] = 0; }
    return p;
}
static uint64_t gcd64(uint64_t a, uint64_t b) { while (b) { uint64_t t = a % b; a = b; b = t; } return a; }
static int calc_trial_divisions(int bits) { return bits <= 512 ? 64 : bits <= 1024 ? 128 : bits <= 2048 ? 384 : bits <= 4096 ? 1024 : 2048; }

// Segmented primorial prefilter on OpenSSL BN. Primes [from_idx, B]; S0 = 5, r = 2.
// Word-sized segment -> BN_mod_word + word GCD; larger segment -> BN_gcd.
struct SegPF {
    std::vector<uint64_t> word; std::vector<BIGNUM*> big; int nwordsegs = 0;
    SegPF(const std::vector<uint64_t>& P, size_t from, uint64_t B) {
        size_t i = from, size = 5;
        while (i < P.size() && P[i] <= B) {
            size_t j = i; unsigned __int128 prod = 1; bool fits = true;
            for (; j < P.size() && P[j] <= B && j - i < size; ++j) { prod *= P[j]; if (prod >> 64) fits = false; }
            if (fits) word.push_back((uint64_t)prod);
            else { BIGNUM* M = BN_new(); BN_set_word(M, 1); for (size_t t = i; t < j; ++t) BN_mul_word(M, P[t]); big.push_back(M); }
            i = j; size *= 2;
        }
    }
    bool passes(const BIGNUM* c, BIGNUM* g, BN_CTX* ctx) const {
        for (uint64_t m : word) if (gcd64(BN_mod_word(c, m), m) != 1) return false;
        for (BIGNUM* M : big) {
            if (BN_num_bits(M) > BN_num_bits(c)) { BN_div(NULL, g, M, c, ctx); BN_gcd(g, g, c, ctx); }   // M mod c, then GCD
            else BN_gcd(g, c, M, ctx);
            if (!BN_is_one(g)) return false;
        }
        return true;
    }
};
// Word primorials: consecutive primes are multiplied while the product fits into 64 bits; BN_mod_word + word GCD, early exit
struct WordPF {
    std::vector<uint64_t> word;
    WordPF(const std::vector<uint64_t>& P, size_t from, uint64_t B) {
        unsigned __int128 prod = 1;
        for (size_t i = from; i < P.size() && P[i] <= B; ++i) {
            if ((prod * P[i]) >> 64) { word.push_back((uint64_t)prod); prod = 1; }
            prod *= P[i];
        }
        if (prod > 1) word.push_back((uint64_t)prod);
    }
    bool passes(const BIGNUM* c, BIGNUM*, BN_CTX*) const {
        for (uint64_t m : word) if (gcd64(BN_mod_word(c, m), m) != 1) return false;
        return true;
    }
};
static long g_cb_count = 0;
static int cb_count(int a, int, BN_GENCB*) { if (a == 0) ++g_cb_count; return 1; }

int main(int argc, char** argv) {
    int bits = std::atoi(argv[1]), K = std::atoi(argv[2]); std::vector<uint64_t> Bs;
    for (char* t = std::strtok(argv[3], ","); t; t = std::strtok(nullptr, ",")) Bs.push_back(std::strtoull(t, nullptr, 10));
    auto P = primes_upto(20000000); int td = calc_trial_divisions(bits);
    BN_CTX* ctx = BN_CTX_new(); BIGNUM *rnd = BN_new(), *g = BN_new();
    std::fprintf(stderr, "%s\n", OpenSSL_version(0));
    std::vector<WordPF*> hyb, ind; for (auto B : Bs) { hyb.push_back(new WordPF(P, td, B)); ind.push_back(new WordPF(P, 1, B)); }
    std::vector<unsigned> mods(td);
    // V1/V2: replica of OpenSSL's probable_prime; if pf != nullptr, the segmented prefilter runs before MR
    auto hybrid = [&](const WordPF* pf, long& mr, long& rej) {
        for (;;) {
            BN_priv_rand_ex(rnd, bits, 1 /*TOP_TWO*/, 1 /*BOTTOM_ODD*/, 0, ctx);
            for (int i = 1; i < td; ++i) mods[i] = (unsigned)BN_mod_word(rnd, P[i]);
            uint64_t delta = 0;
        loop:
            for (int i = 1; i < td; ++i) if ((mods[i] + delta) % P[i] == 0) { delta += 2; goto loop; }
            BN_add_word(rnd, delta);
            if (BN_num_bits(rnd) != bits) continue;
            if (pf && !pf->passes(rnd, g, ctx)) { ++rej; continue; }
            ++mr; if (BN_is_prime_fasttest_ex(rnd, 0, ctx, 0, nullptr) == 1) return;
        }
    };
    auto independent = [&](const WordPF* pf, long& mr, long& rej) {
        for (;;) {
            BN_priv_rand_ex(rnd, bits, 1, 1, 0, ctx);
            if (!pf->passes(rnd, g, ctx)) { ++rej; continue; }
            ++mr; if (BN_is_prime_fasttest_ex(rnd, 0, ctx, 0, nullptr) == 1) return;
        }
    };
    if (argc > 4 && std::string(argv[4]) == "confirm") {        // confirmatory experiment: V0 and V3 only (B = Bs[0]), alternating order
        std::printf("bits,k,method,B,ms,mr_calls,prefilter_rejects\n");
        BN_GENCB* cb0 = BN_GENCB_new(); BN_GENCB_set(cb0, cb_count, nullptr);
        for (int k = 0; k < K; ++k) for (int t = 0; t < 2; ++t) {
            int m = (k % 2 == 0) ? t : 1 - t; long mr = 0, rej = 0; double a = now();
            if (m == 0) { g_cb_count = 0; BN_generate_prime_ex(rnd, bits, 0, nullptr, nullptr, cb0); mr = g_cb_count; }
            else independent(ind[0], mr, rej);
            std::printf("%d,%d,%s,%llu,%.3f,%ld,%ld\n", bits, k, m == 0 ? "V0_OpenSSL" : "V3_independent_WORD",
                        (unsigned long long)(m == 0 ? P[td - 1] : Bs[0]), now() - a, mr, rej); std::fflush(stdout);
        }
        return 0;
    }
    if (argc > 4 && std::string(argv[4]) == "confirm4") {       // confirmatory V0-V3 experiment (protocol of 2026-11-23)
        static const int W[4][4] = {{0, 1, 3, 2}, {1, 2, 0, 3}, {2, 3, 1, 0}, {3, 0, 2, 1}};   // 4x4 Williams square (balanced for first-order carry-over effects)
        static const char* NM[4] = {"V0_OpenSSL", "V1_replica", "V2_hybrid_WORD", "V3_independent_WORD"};
        int warm = argc > 5 ? std::atoi(argv[5]) : 20;
        BN_GENCB* cb4 = BN_GENCB_new(); BN_GENCB_set(cb4, cb_count, nullptr);
        std::printf("bits,k,pos,method,B,ms,mr_calls,prefilter_rejects\n");
        for (int k = -warm; k < K; ++k) for (int pos = 0; pos < 4; ++pos) {
            int m = W[((k % 4) + 4) % 4][pos]; long mr = 0, rej = 0; double a = now(); uint64_t B = P[td - 1];
            if (m == 0) { g_cb_count = 0; BN_generate_prime_ex(rnd, bits, 0, nullptr, nullptr, cb4); mr = g_cb_count; }
            else if (m == 1) hybrid(nullptr, mr, rej);
            else if (m == 2) { hybrid(hyb[0], mr, rej); B = Bs[0]; }
            else { independent(ind[0], mr, rej); B = Bs[0]; }
            double ms = now() - a;
            if (BN_num_bits(rnd) != bits) { std::fprintf(stderr, "ERROR: bit length\n"); return 2; }
            if (k >= 0) { std::printf("%d,%d,%d,%s,%llu,%.3f,%ld,%ld\n", bits, k, pos, NM[m], (unsigned long long)B, ms, mr, rej); std::fflush(stdout); }
        }
        return 0;
    }
    if (argc > 4 && std::string(argv[4]) == "resid") {          // residues: V0 (stock) and V3 (independent sampling + word prefilter)
        std::printf("method,m3,m5,m7,m11,m13\n"); long a = 0, b = 0;
        for (int k = 0; k < K; ++k) for (int m = 0; m < 2; ++m) {
            if (m == 0) BN_generate_prime_ex(rnd, bits, 0, nullptr, nullptr, nullptr); else independent(ind[0], a, b);
            std::printf("%s,%llu,%llu,%llu,%llu,%llu\n", m == 0 ? "V0_OpenSSL" : "V3_independent_WORD", BN_mod_word(rnd, 3), BN_mod_word(rnd, 5), BN_mod_word(rnd, 7), BN_mod_word(rnd, 11), BN_mod_word(rnd, 13));
        }
        return 0;
    }
    BN_GENCB* cb = BN_GENCB_new(); BN_GENCB_set(cb, cb_count, nullptr);
    std::printf("bits,k,method,B,ms,mr_calls,prefilter_rejects\n");
    int nm = 2 + 2 * (int)Bs.size();
    for (int k = 0; k < K; ++k) for (int m0 = 0; m0 < nm; ++m0) {
        int m = (m0 + k) % nm; long mr = 0, rej = 0; double a = now(); std::string name; uint64_t B = 0;
        if (m == 0) { g_cb_count = 0; BN_generate_prime_ex(rnd, bits, 0, nullptr, nullptr, cb); mr = g_cb_count; name = "V0_OpenSSL"; B = P[td - 1]; }
        else if (m == 1) { hybrid(nullptr, mr, rej); name = "V1_replica"; B = P[td - 1]; }
        else { int j = (m - 2) / 2; B = Bs[j];
            if ((m - 2) % 2 == 0) { hybrid(hyb[j], mr, rej); name = "V2_hybrid_WORD"; } else { independent(ind[j], mr, rej); name = "V3_independent_WORD"; } }
        std::printf("%d,%d,%s,%llu,%.2f,%ld,%ld\n", bits, k, name.c_str(), (unsigned long long)B, now() - a, mr, rej); std::fflush(stdout);
    }
}
