// Complete generators (built on spp):
//   prov <bits> <K> <B>    provable prime (Maurer-type recursion + Pocklington), B>0: segmented prefilter, B=0: trial division (1000)
//   prob <bits> <K> <t>    probable prime: independent sampling (S1) + segmented prefilter (B*) + t MR rounds with INDEPENDENT
//                          RANDOM bases (t = 0: number of rounds from mr_rounds_for_error(bits, -128))
//   rsa  <nlen> <K>        pair of RSA primes with the FIPS 186-5 conditions: set S3, eps/2 = 2^-129 for each prime, L(S3) = 2.045,
//                          random bases (same rule as generate_rsa_primes_fips of the library)
// The bases of prob and rsa come from a separate random stream, so the candidate streams (and the numbers of candidates,
// Miller-Rabin calls, and rounds) are identical to those of a run with candidate-derived bases.
#include "spp/spp.hpp"
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
using namespace spp;
static double now() { return std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now().time_since_epoch()).count(); }
struct PStats { long cand = 0, prefilter_rej = 0, modexp = 0, levels = 0; };

// small prime (<= 32 bits) found by complete trial division
static void small_prime(mpz_ptr p, unsigned bits, RandomSource& r) {
    for (;;) {
        random_odd(p, bits, r); unsigned long n = mpz_get_ui(p); bool ok = n > 2;
        for (unsigned long d = 3; d * d <= n && ok; d += 2) if (n % d == 0) ok = false;
        if (ok) return;
    }
}
// random number in [lo, hi]
static void rand_range(mpz_ptr x, mpz_srcptr lo, mpz_srcptr hi, RandomSource& r) {
    BigInt w, t; mpz_sub(w.get(), hi, lo); mpz_add_ui(w.get(), w.get(), 1);
    random_odd(t.get(), (unsigned)mpz_sizeinbase(w.get(), 2) + 64, r); mpz_mod(t.get(), t.get(), w.get()); mpz_add(x, t.get(), lo);
}
// Provable generation: p = 2Rq + 1, q a proven prime of floor(b/2)+1 bits, hence 2q + 1 > sqrt(p) for every b (q > sqrt(p) may
// fail for odd b). Lemma 2: if a^(p-1) = 1 and gcd(a^((p-1)/q) - 1, p) = 1, every prime divisor of p is 1 (mod 2q), i.e.
// >= 2q + 1 > sqrt(p), so p is prime.
static void prov_prime(mpz_ptr p, unsigned b, const SegmentedPrefilter* pf, const TrialDivision* td, RandomSource& r, PStats& st) {
    if (b <= 32) { small_prime(p, b, r); return; }
    ++st.levels;
    BigInt q; prov_prime(q.get(), b / 2 + 1, pf, td, r, st);
    BigInt lo, hi, twoq, t, rem, R, x, y, g, a(2UL), e;
    mpz_mul_2exp(twoq.get(), q.get(), 1);
    mpz_set_ui(t.get(), 1); mpz_mul_2exp(t.get(), t.get(), b - 1); mpz_sub_ui(t.get(), t.get(), 1);      // 2^(b-1) - 1
    mpz_tdiv_q(lo.get(), t.get(), twoq.get()); mpz_tdiv_r(rem.get(), t.get(), twoq.get()); if (mpz_sgn(rem.get())) mpz_add_ui(lo.get(), lo.get(), 1);
    mpz_set_ui(t.get(), 1); mpz_mul_2exp(t.get(), t.get(), b); mpz_sub_ui(t.get(), t.get(), 2);          // 2^b - 2
    mpz_tdiv_q(hi.get(), t.get(), twoq.get());
    for (;;) {
        rand_range(R.get(), lo.get(), hi.get(), r); mpz_mul(p, twoq.get(), R.get()); mpz_add_ui(p, p, 1); ++st.cand;
        if (pf ? !pf->passes(p) : !td->passes(p)) { ++st.prefilter_rej; continue; }
        mpz_mul_2exp(e.get(), R.get(), 1);                       // (p-1)/q = 2R
        mpz_powm(x.get(), a.get(), e.get(), p); ++st.modexp;      // x = a^(2R)
        mpz_powm(y.get(), x.get(), q.get(), p); ++st.modexp;      // y = a^(p-1)
        if (mpz_cmp_ui(y.get(), 1) != 0) continue;                // Fermat condition fails: composite
        mpz_sub_ui(x.get(), x.get(), 1); mpz_gcd(g.get(), x.get(), p);
        if (mpz_cmp_ui(g.get(), 1) == 0) return;                  // Lemma 2 (Pocklington): p is a proven prime
    }
}
int main(int argc, char** argv) {
    std::string cmd = argv[1]; unsigned bits = std::atoi(argv[2]); int K = std::atoi(argv[3]);
    if (cmd == "prov") {
        unsigned long B = std::strtoul(argv[4], nullptr, 10); SegmentedPrefilter* pf = B ? new SegmentedPrefilter(B) : nullptr; TrialDivision td(1000);
        std::printf("bits,k,method,ms,candidates_total,prefilter_rej,modexp,levels\n");
        for (int k = 0; k < K; ++k) {
            SeededRandom r(70000 + bits * 13 + k); PStats st; BigInt p; double a = now();
            prov_prime(p.get(), bits, pf, &td, r, st); double ms = now() - a;
            if (mpz_sizeinbase(p.get(), 2) != bits || !miller_rabin(p.get(), 20).probable_prime) { std::fprintf(stderr, "ERROR\n"); return 2; }
            std::printf("%u,%d,%s,%.2f,%ld,%ld,%ld,%ld\n", bits, k, B ? "prov_GEOM" : "prov_TD1000", ms, st.cand, st.prefilter_rej, st.modexp, st.levels); std::fflush(stdout);
        }
    } else if (cmd == "prob") {
        int t = std::atoi(argv[4]); if (t == 0) t = mr_rounds_for_error(bits, -128.0); SegmentedPrefilter pf(recommended_bound(bits));
        std::printf("bits,k,method,ms,candidates,mr_calls,mr_rounds\n");
        for (int k = 0; k < K; ++k) {
            SeededRandom r(70000 + bits * 13 + k), rb(170000 + bits * 13 + k); GenStats st;   // rb: separate stream for the bases (the candidate stream does not depend on t)
            BigInt p; double a = now();
            for (;;) {      // same candidate stream as generate_prime, but random bases
                random_odd(p.get(), bits, r); ++st.candidates; if (!pf.passes(p.get())) { ++st.prefilter_rejected; continue; }
                ++st.mr_calls; auto m = miller_rabin(p.get(), t, BaseMode::Random, &rb); st.mr_rounds += m.rounds; if (m.probable_prime) break; }
            std::printf("%u,%d,prob_t%d,%.2f,%llu,%llu,%llu\n", bits, k, t, now() - a, (unsigned long long)st.candidates, (unsigned long long)st.mr_calls, (unsigned long long)st.mr_rounds); std::fflush(stdout);
        }
    } else if (cmd == "rsa") {      // probable primes in the style of FIPS 186-5 A.1.3 (with its checks)
        unsigned hb = bits / 2; int t = mr_rounds_for_error(hb, -129.0, loss_bound(CandidateSet::S3_RsaCoprime)); SegmentedPrefilter pf(recommended_bound(hb));
        std::printf("nlen,k,t,ms,gcd_rej,pq_rej,top_ok\n");
        BigInt sq2;                    // instead of comparing with floor(sqrt(2) * 2^(hb-1)), the two top bits are set
        for (int k = 0; k < K; ++k) {
            SeededRandom r(80000 + bits * 7 + k), rb(180000 + bits * 7 + k); BigInt p, q, d, pm1, g, e(65537UL), lim; long gr = 0, pr = 0; double a = now();
            auto one = [&](mpz_ptr x) { for (;;) {
                random_odd(x, hb, r); mpz_setbit(x, hb - 2);                       // two top bits = 11 (FIPS 186-5: p >= sqrt(2)*2^(hb-1))
                mpz_sub_ui(pm1.get(), x, 1); mpz_gcd(g.get(), pm1.get(), e.get());
                if (mpz_cmp_ui(g.get(), 1) != 0) { ++gr; continue; }               // condition gcd(e, p-1) = 1 checked before the test
                if (!pf.passes(x)) continue;
                if (miller_rabin(x, t, BaseMode::Random, &rb).probable_prime) return; } };
            one(p.get());
            mpz_set_ui(lim.get(), 1); mpz_mul_2exp(lim.get(), lim.get(), hb - 100);
            for (;;) { one(q.get()); mpz_sub(d.get(), p.get(), q.get()); if (mpz_sgn(d.get()) < 0) mpz_neg(d.get(), d.get()); if (mpz_cmp(d.get(), lim.get()) > 0) break; ++pr; }
            // two top bits: p >= 1.5 * 2^(hb-1) > sqrt(2) * 2^(hb-1)
            bool top = mpz_tstbit(p.get(), hb - 1) && mpz_tstbit(p.get(), hb - 2) && mpz_tstbit(q.get(), hb - 1) && mpz_tstbit(q.get(), hb - 2);
            std::printf("%u,%d,%d,%.2f,%ld,%ld,%d\n", bits, k, t, now() - a, gr, pr, top ? 1 : 0); std::fflush(stdout);
        }
    }
}
