// Experiment driver for the C++ measurements of the paper (built on the spp library).
//   diss_bench consts|sweep|strat|s0|safe|safesweep|gen|gaps|resid  [parameters]  -> CSV stdout
//   e.g. diss_bench gen <bits> <K> <B> <t> [k0];  diss_bench safe <bits> <K> <B> [k0]   (seeds k0..k0+K-1)
#include "spp/spp.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
using namespace spp;
static double now() { return std::chrono::duration<double, std::nano>(std::chrono::steady_clock::now().time_since_epoch()).count(); }
static double median(std::vector<double> v) { std::sort(v.begin(), v.end()); return v[v.size() / 2]; }
static std::vector<BigInt> cands(unsigned bits, size_t n, uint64_t seed) {
    SeededRandom r(seed); std::vector<BigInt> c(n); for (auto& x : c) random_odd(x.get(), bits, r); return c;
}
template <class F> static double per(const F& f, const std::vector<BigInt>& c, int reps) {
    std::vector<double> t; volatile unsigned long s = 0;
    for (int k = 0; k < reps; ++k) { double a = now(); for (auto& n : c) s += f(n.get()); t.push_back((now() - a) / c.size()); }
    return median(t);
}
static double mr_cost(unsigned bits, int k, uint64_t seed) {
    auto c = cands(bits, k, seed); std::vector<double> t;
    for (auto& n : c) { double a = now(); miller_rabin(n.get(), 1); t.push_back(now() - a); }
    return median(t);
}
static std::vector<unsigned long> grid() {
    std::vector<unsigned long> g; for (double e = 2.5; e <= 6.51; e += 0.25) g.push_back((unsigned long)std::llround(std::pow(10, e))); return g;
}
static double survival(unsigned long B) { double s = 1; for (auto p : odd_primes_upto(B)) s *= 1.0 - 1.0 / p; return s; }

int main(int argc, char** argv) {
    std::string cmd = argv[1];
    if (cmd == "consts") {              // bits c_MR c_mod beta alpha
        std::printf("bits,c_mr_ns,c_mod_ns,beta_ns_per_bit,alpha_ns\n");
        auto P = odd_primes_upto(2000000);
        for (int i = 2; i < argc; ++i) {
            unsigned bits = std::atoi(argv[i]);
            double cmr = mr_cost(bits, bits > 2500 ? 15 : 40, 101 + bits);
            auto c = cands(bits, 400, 7 + bits);
            TrialDivision td(997);
            // c_mod: cost of one residue modulo a small prime, measured on candidates without small divisors
            std::vector<BigInt> surv; for (auto& n : c) if (td.passes(n.get())) surv.push_back(n);
            double tdc = per([&](mpz_srcptr n) { unsigned long s = 0; for (auto p : odd_primes_upto(997)) s += mpz_divisible_ui_p(n, p); return s; }, surv, 5) / 167.0;
            std::vector<double> xs, ys; BigInt M;
            for (int k : {256, 512, 1024, 2048, 4096, 8192, 16384}) {
                mpz_set_ui(M.get(), 1); for (int j = 0; j < k; ++j) mpz_mul_ui(M.get(), M.get(), P[j]);
                BigInt g; double t = per([&](mpz_srcptr n) { mpz_gcd(g.get(), n, M.get()); return 1; }, c, 3);
                xs.push_back((double)mpz_sizeinbase(M.get(), 2)); ys.push_back(t);
            }
            double mx = 0, my = 0; for (size_t j = 0; j < xs.size(); ++j) { mx += xs[j]; my += ys[j]; } mx /= xs.size(); my /= ys.size();
            double sxy = 0, sxx = 0; for (size_t j = 0; j < xs.size(); ++j) { sxy += (xs[j] - mx) * (ys[j] - my); sxx += (xs[j] - mx) * (xs[j] - mx); }
            double beta = sxy / sxx, alpha = my - beta * mx;
            std::printf("%u,%.1f,%.2f,%.4f,%.1f\n", bits, cmr, tdc, beta, alpha); std::fflush(stdout);
        }
    } else if (cmd == "sweep") {         // bits B segs pass geom_ns td_ns pi
        unsigned bits = std::atoi(argv[2]); size_t n = std::atoi(argv[3]);
        auto c = cands(bits, n, 2026 + bits);
        std::printf("bits,B,segments,pass,geom_ns,td_ns,pi\n");
        for (auto B : grid()) {
            SegmentedPrefilter pf(B); size_t pass = 0; for (auto& x : c) pass += pf.passes(x.get());
            double g = per([&](mpz_srcptr x) { return pf.passes(x); }, c, 5), t = -1;
            if (B <= 100000) { TrialDivision td(B); std::vector<BigInt> sub(c.begin(), c.begin() + std::min<size_t>(n, 1500)); t = per([&](mpz_srcptr x) { return td.passes(x); }, sub, 3); }
            std::printf("%u,%lu,%zu,%.5f,%.1f,%.1f,%.6f\n", bits, B, pf.segment_count(), (double)pass / n, g, t, survival(B)); std::fflush(stdout);
        }
    } else if (cmd == "strat") {         // B=1000 strategies, interleaved reps
        unsigned bits = std::atoi(argv[2]); auto c = cands(bits, 5000, 77 + bits);
        TrialDivision td(1000);
        SegmentedPrefilter sp(1000, 100000, 1.0), f16(1000, 16, 1.0), g3(1000, 3, 2.0), g5(1000, 5, 2.0), g8(1000, 8, 2.0), g5r3(1000, 5, 3.0), g5r15(1000, 5, 1.5);
        const char* nm[] = {"TD", "SP", "SEG-fixed16", "GEOM S0=3", "GEOM S0=5", "GEOM S0=8", "GEOM S0=5 r=1.5", "GEOM S0=5 r=3"};
        std::vector<std::vector<double>> T(8);
        for (int rep = 0; rep < 21; ++rep) for (int m = 0; m < 8; ++m) {
            int k = (m + rep) % 8; volatile unsigned long s = 0; double a = now();
            for (auto& x : c) { mpz_srcptr n = x.get();
                switch (k) { case 0: s += td.passes(n); break; case 1: s += sp.passes(n); break; case 2: s += f16.passes(n); break; case 3: s += g3.passes(n); break;
                             case 4: s += g5.passes(n); break; case 5: s += g8.passes(n); break; case 6: s += g5r15.passes(n); break; default: s += g5r3.passes(n); } }
            T[k].push_back((now() - a) / c.size());
        }
        std::printf("bits,method,median_ns,q1_ns,q3_ns\n");
        for (int k = 0; k < 8; ++k) { auto v = T[k]; std::sort(v.begin(), v.end()); std::printf("%u,%s,%.1f,%.1f,%.1f\n", bits, nm[k], v[10], v[5], v[15]); }
    } else if (cmd == "s0") {
        unsigned bits = std::atoi(argv[2]); auto c = cands(bits, 3000, 99 + bits);
        std::printf("bits,ratio,S0,first_bits,median_ns\n");
        for (double r : {1.5, 2.0, 3.0}) for (int s0 = 1; s0 <= 28; ++s0) {
            SegmentedPrefilter pf(1000, s0, r);
            std::printf("%u,%.1f,%d,%zu,%.1f\n", bits, r, s0, pf.segment_bits()[0], per([&](mpz_srcptr x) { return pf.passes(x); }, c, 9)); std::fflush(stdout);
        }
    } else if (cmd == "safesweep") {    // pair pass fraction + cost per q-candidate
        unsigned bits = std::atoi(argv[2]); size_t n = std::atoi(argv[3]);
        SeededRandom r(555 + bits); std::vector<BigInt> Q(n), QP(n);
        for (size_t i = 0; i < n; ++i) { random_odd(Q[i].get(), bits - 1, r); BigInt p; mpz_mul_2exp(p.get(), Q[i].get(), 1); mpz_add_ui(p.get(), p.get(), 1); mpz_mul(QP[i].get(), Q[i].get(), p.get()); }
        double cmr = mr_cost(bits - 1, 30, 3 + bits);
        std::printf("bits,B,method,pre_ns,pass,pi2,pi1,cmr_ns\n");
        for (auto B : grid()) {
            auto P = odd_primes_upto(B); double pi1 = 1, pi2 = 1; for (auto p : P) { pi1 *= 1 - 1.0 / p; pi2 *= 1 - 2.0 / p; }
            SegmentedPrefilter pf(B); size_t pass = 0; for (auto& x : QP) pass += pf.passes(x.get());
            std::vector<BigInt> sub(QP.begin(), QP.begin() + std::min<size_t>(n, 3000));
            std::printf("%u,%lu,geom_product,%.1f,%.6f,%.6f,%.6f,%.1f\n", bits, B, per([&](mpz_srcptr x) { return pf.passes(x); }, sub, 3), (double)pass / n, pi2, pi1, cmr);
            if (B <= 30000) {
                std::vector<BigInt> sq(Q.begin(), Q.begin() + std::min<size_t>(n, 3000));
                double t = per([&](mpz_srcptr q) { for (auto p : P) { unsigned long rr = mpz_fdiv_ui(q, p); if (rr == 0 || rr == (p - 1) / 2) return 0; } return 1; }, sq, 3);
                std::printf("%u,%lu,td_pair,%.1f,%.6f,%.6f,%.6f,%.1f\n", bits, B, t, (double)pass / n, pi2, pi1, cmr);
            }
            std::fflush(stdout);
        }
    } else if (cmd == "safe") {          // generate safe primes with strategies
        unsigned bits = std::atoi(argv[2]); int K = std::atoi(argv[3]); unsigned long Bg = std::strtoul(argv[4], nullptr, 10);
        int k0 = argc > 5 ? std::atoi(argv[5]) : 0;   // optional first seed index (for splitting replicates)
        TrialDivision td(1000); SegmentedPrefilter pf(Bg); auto P1 = odd_primes_upto(1000);
        std::printf("bits,k,strategy,candidates,rejected,modexp,ms\n");
        const char* S[] = {"naive", "td_pair", "geom_product"};
        for (int k = k0; k < k0 + K; ++k) for (int s = 0; s < 3; ++s) {
            SeededRandom r(9000 + 100 * bits + k); BigInt q, p, qp, two(2UL), e, x; long cand = 0, rej = 0, mexp = 0; double a = now();
            for (;;) {
                random_odd(q.get(), bits - 1, r); ++cand;
                mpz_mul_2exp(p.get(), q.get(), 1); mpz_add_ui(p.get(), p.get(), 1);
                bool ok;
                if (s == 0) ok = td.passes(q.get());
                else if (s == 1) { ok = true; for (auto pr : P1) { unsigned long rr = mpz_fdiv_ui(q.get(), pr); if (rr == 0 || rr == (pr - 1) / 2) { ok = false; break; } } }
                else { mpz_mul(qp.get(), q.get(), p.get()); ok = pf.passes(qp.get()); }
                if (!ok) { ++rej; continue; }
                auto m1 = miller_rabin(q.get(), 1); mexp += m1.rounds; if (!m1.probable_prime) continue;
                mpz_sub_ui(e.get(), p.get(), 1); mpz_powm(x.get(), two.get(), e.get(), p.get()); ++mexp;
                if (mpz_cmp_ui(x.get(), 1) != 0) continue;
                auto m2 = miller_rabin(q.get(), 40); mexp += m2.rounds; if (m2.probable_prime) break;
            }
            std::printf("%u,%d,%s,%ld,%ld,%ld,%.2f\n", bits, k, S[s], cand, rej, mexp, (now() - a) / 1e6); std::fflush(stdout);
        }
    } else if (cmd == "gen") {           // independent generation: prefilter vs plain on same streams
        unsigned bits = std::atoi(argv[2]); int K = std::atoi(argv[3]); unsigned long B = std::strtoul(argv[4], nullptr, 10); int t = std::atoi(argv[5]);
        int k0 = argc > 6 ? std::atoi(argv[6]) : 0;   // optional first seed index (for splitting replicates)
        SegmentedPrefilter pf(B), p1(1000); TrialDivision td(1000);
        std::printf("bits,k,method,candidates,mr_calls,mr_rounds,ms\n");
        for (int k = k0; k < k0 + K; ++k) {
            BigInt a, b; GenStats sa, sb; SeededRandom r1(31000 + bits * 7 + k), r2(31000 + bits * 7 + k);
            double t0 = now(); generate_prime(a.get(), bits, pf, r1, t, &sa); double t1 = now(); generate_prime_plain(b.get(), bits, r2, t, &sb); double t2 = now();
            // TD(1000) on the same stream
            SeededRandom r3(31000 + bits * 7 + k); GenStats sc; BigInt c; double t3 = now();
            for (;;) { random_odd(c.get(), bits, r3); ++sc.candidates; if (!td.passes(c.get())) continue; ++sc.mr_calls; auto m = miller_rabin(c.get(), t); sc.mr_rounds += m.rounds; if (m.probable_prime) break; }
            double t4 = now();
            if (mpz_cmp(a.get(), b.get()) || mpz_cmp(a.get(), c.get())) { std::fprintf(stderr, "MISMATCH\n"); return 2; }
            std::printf("%u,%d,GEOM,%llu,%llu,%llu,%.3f\n", bits, k, (unsigned long long)sa.candidates, (unsigned long long)sa.mr_calls, (unsigned long long)sa.mr_rounds, (t1 - t0) / 1e6);
            std::printf("%u,%d,MR,%llu,%llu,%llu,%.3f\n", bits, k, (unsigned long long)sb.candidates, (unsigned long long)sb.mr_calls, (unsigned long long)sb.mr_rounds, (t2 - t1) / 1e6);
            std::printf("%u,%d,TD,%llu,%llu,%llu,%.3f\n", bits, k, (unsigned long long)sc.candidates, (unsigned long long)sc.mr_calls, (unsigned long long)sc.mr_rounds, (t4 - t3) / 1e6);
            std::fflush(stdout);
        }
    } else if (cmd == "gaps") {          // size-biased gap samples
        unsigned bits = std::atoi(argv[2]); int S = std::atoi(argv[3]);
        SegmentedPrefilter pf(recommended_bound(bits)); SeededRandom r(4242 + bits);
        std::printf("bits,gap\n");
        auto isp = [&](mpz_srcptr n) { return pf.passes(n) && miller_rabin(n, 20).probable_prime; };
        for (int i = 0; i < S; ++i) {
            BigInt x, a, b; random_odd(x.get(), bits, r);
            mpz_set(a.get(), x.get()); while (!isp(a.get())) mpz_add_ui(a.get(), a.get(), 2);
            mpz_sub_ui(b.get(), x.get(), 2); while (!isp(b.get())) mpz_sub_ui(b.get(), b.get(), 2);
            BigInt g; mpz_sub(g.get(), a.get(), b.get()); std::printf("%u,%lu\n", bits, mpz_get_ui(g.get()));
        }
    } else if (cmd == "resid") {         // residues of primes: independent and incremental
        unsigned bits = std::atoi(argv[2]); int N = std::atoi(argv[3]);
        SegmentedPrefilter pf(recommended_bound(bits)); SeededRandom r1(6060 + bits), r2(7070 + bits);
        std::printf("method,m3,m5,m7,m11,m13\n");
        BigInt p;
        for (int i = 0; i < N; ++i) { generate_prime(p.get(), bits, pf, r1, 20);
            std::printf("independent,%lu,%lu,%lu,%lu,%lu\n", mpz_fdiv_ui(p.get(), 3), mpz_fdiv_ui(p.get(), 5), mpz_fdiv_ui(p.get(), 7), mpz_fdiv_ui(p.get(), 11), mpz_fdiv_ui(p.get(), 13)); }
        for (int i = 0; i < N; ++i) { generate_prime_incremental(p.get(), bits, pf, r2, 20);
            std::printf("incremental,%lu,%lu,%lu,%lu,%lu\n", mpz_fdiv_ui(p.get(), 3), mpz_fdiv_ui(p.get(), 5), mpz_fdiv_ui(p.get(), 7), mpz_fdiv_ui(p.get(), 11), mpz_fdiv_ui(p.get(), 13)); }
    }
    return 0;
}
