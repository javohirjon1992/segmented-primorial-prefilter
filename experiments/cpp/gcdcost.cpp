// Cost of BN_gcd and BN_mod_word in OpenSSL for 1024- and 1536-bit operands (usage: ./gcdcost).
#include <initializer_list>
#include <chrono>
#include <cstdio>
extern "C" { typedef struct bignum_st BIGNUM; typedef struct bignum_ctx BN_CTX;
BN_CTX* BN_CTX_new(void); BIGNUM* BN_new(void); int BN_priv_rand_ex(BIGNUM*, int, int, int, unsigned int, BN_CTX*);
int BN_gcd(BIGNUM*, const BIGNUM*, const BIGNUM*, BN_CTX*); int BN_div(BIGNUM*, BIGNUM*, const BIGNUM*, const BIGNUM*, BN_CTX*);
unsigned long long BN_mod_word(const BIGNUM*, unsigned long long); }
int main() { BN_CTX* ctx = BN_CTX_new(); BIGNUM *a = BN_new(), *b = BN_new(), *g = BN_new(), *M = BN_new();
  for (int bits : {1024, 1536}) { BN_priv_rand_ex(a, bits, 1, 1, 0, ctx); BN_priv_rand_ex(b, bits - 3, 1, 1, 0, ctx);
    auto t0 = std::chrono::steady_clock::now(); for (int i = 0; i < 2000; ++i) BN_gcd(g, a, b, ctx);
    double us = std::chrono::duration<double, std::micro>(std::chrono::steady_clock::now() - t0).count() / 2000;
    t0 = std::chrono::steady_clock::now(); volatile unsigned long long s = 0; for (int i = 0; i < 200000; ++i) s += BN_mod_word(a, 0xfffffffffffffc5ULL);
    double mw = std::chrono::duration<double, std::nano>(std::chrono::steady_clock::now() - t0).count() / 200000;
    BN_priv_rand_ex(M, 100000, 1, 1, 0, ctx); t0 = std::chrono::steady_clock::now(); for (int i = 0; i < 200; ++i) BN_div(NULL, g, M, a, ctx);
    double dv = std::chrono::duration<double, std::micro>(std::chrono::steady_clock::now() - t0).count() / 200;
    std::printf("%d bit: BN_gcd %.1f us, BN_mod_word %.0f ns, 100k-bit mod %.1f us\n", bits, us, mw, dv); } }
