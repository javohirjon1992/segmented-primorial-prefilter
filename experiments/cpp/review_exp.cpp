// Replication of seed-dependent counts and check of the provable generator (spp 1.1.1).
//   provodd <K> <bits...>        provable primes at odd bit lengths; each is checked with GMP's mpz_probab_prime_p
//   genrep <bits> <K> <B> <t> [plain=1]   GEOM(B) and TD(1000) [and MR only] on the same candidate stream: Miller-Rabin calls
//   saferep <bits> <K> <B>       safe primes: td_pair and geom_product on the same stream: modular exponentiations
// Seeds: genrep 31000 + 7*bits + k (as diss_bench gen), saferep 9000 + 100*bits + k (as diss_bench safe),
// provodd 5000000 + 10000*bits + k. The counts (candidates, Miller-Rabin calls, exponentiations) do not depend on the platform.
#include "spp/spp.hpp"
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>
using namespace spp;
static double now() { return std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now().time_since_epoch()).count(); }
int main(int argc, char** argv) {
    std::string cmd = argv[1];
    if (cmd == "provodd") {
        int K = std::atoi(argv[2]);
        std::printf("bits,k,ms,candidates,prefilter_rej,modexps,levels,levels_q_below_sqrt,gmp_check\n");
        for (int i = 3; i < argc; ++i) {
            unsigned bits = std::atoi(argv[i]); SegmentedPrefilter pf(bits < 256 ? 1000 : recommended_bound(bits));
            for (int k = 0; k < K; ++k) {
                SeededRandom r(5000000ULL + 10000ULL * bits + k); ProvableStats st; BigInt p; double a = now();
                generate_provable_prime(p.get(), bits, pf, r, &st); double ms = now() - a;
                int chk = mpz_probab_prime_p(p.get(), 25);
                if (mpz_sizeinbase(p.get(), 2) != bits || chk == 0) { std::fprintf(stderr, "ERROR bits=%u k=%d\n", bits, k); return 2; }
                std::printf("%u,%d,%.3f,%lu,%lu,%lu,%lu,%lu,%d\n", bits, k, ms, st.candidates, st.prefilter_rejected, st.modexps, st.levels, st.q_below_sqrt, chk);
            }
            std::fflush(stdout);
        }
    } else if (cmd == "genrep") {
        unsigned bits = std::atoi(argv[2]); int K = std::atoi(argv[3]); unsigned long B = std::strtoul(argv[4], nullptr, 10); int t = std::atoi(argv[5]);
        bool plain = argc > 6 ? std::atoi(argv[6]) != 0 : true;
        SegmentedPrefilter pf(B); TrialDivision td(1000);
        std::printf("bits,k,method,candidates,mr_calls,mr_rounds,composite_calls,composite_rounds_beyond_first,ms\n");
        auto out = [&](int k, const char* m, const GenStats& s, double ms) {
            unsigned long long cc = s.mr_calls - 1, beyond = s.mr_rounds - t - cc;
            std::printf("%u,%d,%s,%llu,%llu,%llu,%llu,%llu,%.3f\n", bits, k, m, (unsigned long long)s.candidates, (unsigned long long)s.mr_calls,
                        (unsigned long long)s.mr_rounds, cc, beyond, ms); };
        for (int k = 0; k < K; ++k) {
            BigInt a, b, c; GenStats sa, sb, sc; SeededRandom r1(31000 + bits * 7 + k), r3(31000 + bits * 7 + k);
            double t0 = now(); generate_prime(a.get(), bits, pf, r1, t, &sa); double t1 = now();
            for (;;) { random_odd(c.get(), bits, r3); ++sc.candidates; if (!td.passes(c.get())) continue; ++sc.mr_calls; auto m = miller_rabin(c.get(), t); sc.mr_rounds += m.rounds; if (m.probable_prime) break; }
            double t2 = now();
            if (mpz_cmp(a.get(), c.get())) { std::fprintf(stderr, "MISMATCH\n"); return 2; }
            out(k, "GEOM", sa, t1 - t0); out(k, "TD", sc, t2 - t1);
            if (plain) { SeededRandom r2(31000 + bits * 7 + k); double t3 = now(); generate_prime_plain(b.get(), bits, r2, t, &sb);
                if (mpz_cmp(a.get(), b.get())) { std::fprintf(stderr, "MISMATCH\n"); return 2; } out(k, "MR", sb, now() - t3); }
            std::fflush(stdout);
        }
    } else if (cmd == "saferep") {
        unsigned bits = std::atoi(argv[2]); int K = std::atoi(argv[3]); unsigned long Bg = std::strtoul(argv[4], nullptr, 10);
        SegmentedPrefilter pf(Bg); auto P1 = odd_primes_upto(1000);
        std::printf("bits,k,strategy,candidates,rejected,modexp,ms\n");
        const char* S[] = {"td_pair", "geom_product"};
        for (int k = 0; k < K; ++k) { BigInt first;
          for (int s = 0; s < 2; ++s) {
            SeededRandom r(9000 + 100 * bits + k); BigInt q, p, qp, two(2UL), e, x; long cand = 0, rej = 0, mexp = 0; double a = now();
            for (;;) {
                random_odd(q.get(), bits - 1, r); ++cand;
                mpz_mul_2exp(p.get(), q.get(), 1); mpz_add_ui(p.get(), p.get(), 1);
                bool ok;
                if (s == 0) { ok = true; for (auto pr : P1) { unsigned long rr = mpz_fdiv_ui(q.get(), pr); if (rr == 0 || rr == (pr - 1) / 2) { ok = false; break; } } }
                else { mpz_mul(qp.get(), q.get(), p.get()); ok = pf.passes(qp.get()); }
                if (!ok) { ++rej; continue; }
                auto m1 = miller_rabin(q.get(), 1); mexp += m1.rounds; if (!m1.probable_prime) continue;
                mpz_sub_ui(e.get(), p.get(), 1); mpz_powm(x.get(), two.get(), e.get(), p.get()); ++mexp;
                if (mpz_cmp_ui(x.get(), 1) != 0) continue;
                auto m2 = miller_rabin(q.get(), 40); mexp += m2.rounds; if (m2.probable_prime) break;
            }
            double ms = now() - a;
            if (s == 0) mpz_set(first.get(), p.get()); else if (mpz_cmp(first.get(), p.get())) { std::fprintf(stderr, "MISMATCH\n"); return 2; }
            if (!mpz_probab_prime_p(p.get(), 25) || !mpz_probab_prime_p(q.get(), 25)) { std::fprintf(stderr, "ERROR\n"); return 2; }
            std::printf("%u,%d,%s,%ld,%ld,%ld,%.2f\n", bits, k, S[s], cand, rej, mexp, ms); std::fflush(stdout);
          }
        }
    }
    return 0;
}
