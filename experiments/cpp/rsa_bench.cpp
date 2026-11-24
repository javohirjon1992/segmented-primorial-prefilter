// Complete RSA key generation: spp (segmented prefilter at B*, trial division at B = 1000) and OpenSSL RSA_generate_key_ex.
#include "spp/spp.hpp"
#include <chrono>
#include <cstdio>
#include <cstdlib>
extern "C" {
typedef struct rsa_st RSA; typedef struct bignum_st BIGNUM;
RSA* RSA_new(void); void RSA_free(RSA*); int RSA_generate_key_ex(RSA*, int, BIGNUM*, void*);
BIGNUM* BN_new(void); void BN_free(BIGNUM*); int BN_set_word(BIGNUM*, unsigned long long);
}
using namespace spp;
static double now() { return std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now().time_since_epoch()).count(); }
// p: bits/2 bits, until gcd(e, p-1) = 1
static void gen_p(mpz_ptr p, unsigned hb, const SegmentedPrefilter* pf, const TrialDivision* td, RandomSource& r, int t, GenStats& st) {
    BigInt pm1, g;
    for (;;) {
        if (pf) generate_prime(p, hb, *pf, r, t, &st);
        else { for (;;) { random_odd(p, hb, r); ++st.candidates; if (!td->passes(p)) continue; ++st.mr_calls; auto m = miller_rabin(p, t); st.mr_rounds += m.rounds; if (m.probable_prime) break; } }
        mpz_sub_ui(pm1.get(), p, 1); if (mpz_fdiv_ui(pm1.get(), 65537) != 0) return;
    }
}
int main(int argc, char** argv) {
    unsigned bits = std::atoi(argv[1]); int K = std::atoi(argv[2]); unsigned hb = bits / 2; int t = hb > 1024 ? 128 : 64;   // round counts of OpenSSL 3
    SegmentedPrefilter pf(recommended_bound(hb)); TrialDivision td(1000);
    std::printf("bits,k,method,ms,candidates,mr_calls\n");
    for (int k = 0; k < K; ++k) {
        for (int m = 0; m < 3; ++m) {
            int mm = (m + k) % 3; double a = now(); GenStats st;
            if (mm < 2) { SeededRandom r(4000 + bits * 10 + k); BigInt p, q, n, e(65537UL), l, d, pm1, qm1, g;
                gen_p(p.get(), hb, mm == 0 ? &pf : nullptr, &td, r, t, st); do gen_p(q.get(), hb, mm == 0 ? &pf : nullptr, &td, r, t, st); while (!mpz_cmp(p.get(), q.get()));
                mpz_mul(n.get(), p.get(), q.get()); mpz_sub_ui(pm1.get(), p.get(), 1); mpz_sub_ui(qm1.get(), q.get(), 1);
                mpz_gcd(g.get(), pm1.get(), qm1.get()); mpz_mul(l.get(), pm1.get(), qm1.get()); mpz_tdiv_q(l.get(), l.get(), g.get());
                mpz_invert(d.get(), e.get(), l.get());
            } else { RSA* rsa = RSA_new(); BIGNUM* e = BN_new(); BN_set_word(e, 65537); RSA_generate_key_ex(rsa, bits, e, nullptr); BN_free(e); RSA_free(rsa); }
            std::printf("%u,%d,%s,%.2f,%llu,%llu\n", bits, k, mm == 0 ? "spp_GEOM" : mm == 1 ? "spp_TD" : "OpenSSL", now() - a,
                        (unsigned long long)st.candidates, (unsigned long long)st.mr_calls); std::fflush(stdout);
        }
    }
}
