// spp_gen — generates k primes of a given bit length with two methods and compares their times.
// Both methods test EXACTLY THE SAME candidate stream (seed + i for the i-th prime).
//
//   spp_gen 1024            one 1024-bit prime
//   spp_gen 2048 10         ten 2048-bit primes
//   spp_gen 1024 1000 12345 with a seed (reproducible)
#include "common.hpp"
#include "spp/spp.hpp"

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <random>

using namespace spp;
using tools::fmt_time;
using tools::group;

int main(int argc, char** argv) {
    if (argc < 2 || std::string(argv[1]) == "-h" || std::string(argv[1]) == "--help") {
        std::printf("Usage: spp_gen <bit length> [k] [seed] [--no-file]\n");
        return argc < 2 ? 1 : 0;
    }
    unsigned bits = static_cast<unsigned>(std::strtoul(argv[1], nullptr, 10));
    unsigned long k = argc > 2 ? std::strtoul(argv[2], nullptr, 10) : 1;
    std::uint64_t seed;
    bool save = true;
    if (argc > 3 && std::string(argv[3]) != "--no-file") seed = std::strtoull(argv[3], nullptr, 10);
    else { SystemRandom sr; seed = sr.next_u64() >> 16; }
    for (int i = 1; i < argc; ++i) if (std::string(argv[i]) == "--no-file") save = false;
    if (bits < 16 || k < 1) { std::fprintf(stderr, "the bit length must be >= 16 and k >= 1\n"); return 1; }

    const int rounds = kDefaultRounds;
    SegmentedPrefilter pf(recommended_bound(bits));
    std::printf("================================================================================================\n");
    std::printf("%lu primes of %u bits | GMP %s | spp %s\n", k, bits, gmp_version, version());
    std::printf("Method of the paper: B = %lu (%zu odd primes, %zu segments, S0 = %u, r = %g); setup: %s\n",
                pf.bound(), pf.prime_count(), pf.segment_count(), kDefaultS0, kDefaultRatio,
                fmt_time(pf.build_seconds() * 1e9).c_str());
    std::printf("Miller-Rabin: %d rounds | seed = %llu\n", rounds, static_cast<unsigned long long>(seed));
    std::printf("------------------------------------------------------------------------------------------------\n");
    bool verbose = k <= 50;
    if (verbose)
        std::printf("%5s %8s %12s %18s %16s %12s %12s %9s\n", "#", "cand.", "prefilter rej.", "MR calls (A/B)",
                    "MR rounds (A/B)", "A: paper", "B: plain MR", "speedup");

    GenStats SA, SB;
    std::vector<double> TA, TB;
    std::vector<std::string> primes;
    bool all_same = true;
    BigInt pa, pb;
    for (unsigned long i = 0; i < k; ++i) {
        GenStats sa, sb;
        double ta, tb;
        SeededRandom ra(seed + i), rb(seed + i);
        auto runA = [&] { double t0 = tools::now_ns(); generate_prime(pa.get(), bits, pf, ra, rounds, &sa); ta = tools::now_ns() - t0; };
        auto runB = [&] { double t0 = tools::now_ns(); generate_prime_plain(pb.get(), bits, rb, rounds, &sb); tb = tools::now_ns() - t0; };
        if (i % 2 == 0) { runA(); runB(); } else { runB(); runA(); }   // alternate the order
        bool same = mpz_cmp(pa.get(), pb.get()) == 0 && sa.candidates == sb.candidates;
        all_same = all_same && same;
        SA.add(sa); SB.add(sb); TA.push_back(ta); TB.push_back(tb);
        primes.push_back(pa.str());
        if (verbose)
            std::printf("%5lu %8llu %12llu %8llu / %-8llu %7llu / %-7llu %12s %12s %8.2fx%s\n", i + 1,
                        (unsigned long long)sa.candidates, (unsigned long long)sa.prefilter_rejected,
                        (unsigned long long)sa.mr_calls, (unsigned long long)sb.mr_calls,
                        (unsigned long long)sa.mr_rounds, (unsigned long long)sb.mr_rounds,
                        fmt_time(ta).c_str(), fmt_time(tb).c_str(), tb / ta, same ? "" : "  <-- MISMATCH!");
        else if ((i + 1) % std::max<unsigned long>(1, k / 20) == 0)
            std::printf("  ... %lu / %lu done\n", i + 1, k), std::fflush(stdout);
    }
    double sumA = 0, sumB = 0;
    for (double t : TA) sumA += t;
    for (double t : TB) sumB += t;
    std::printf("------------------------------------------------------------------------------------------------\n");
    std::printf("SUMMARY\n");
    std::printf("  Both methods tested the same candidates and found the same primes: %s\n", all_same ? "YES" : "NO");
    std::printf("  Candidates tested: %s (%.1f per prime; expected ~ %.1f)\n",
                group(SA.candidates).c_str(), double(SA.candidates) / k, bits * std::log(2.0) / 2);
    std::printf("  Rejected by the prefilter: %s (%.1f%%)\n", group(SA.prefilter_rejected).c_str(), 100.0 * SA.prefilter_rejected / SA.candidates);
    std::printf("  Miller-Rabin calls : paper %s | plain %s (%.1f%% fewer)\n", group(SA.mr_calls).c_str(),
                group(SB.mr_calls).c_str(), 100.0 * (1 - double(SA.mr_calls) / SB.mr_calls));
    std::printf("  Miller-Rabin rounds: paper %s | plain %s (%.1f%% fewer)\n", group(SA.mr_rounds).c_str(),
                group(SB.mr_rounds).c_str(), 100.0 * (1 - double(SA.mr_rounds) / SB.mr_rounds));
    std::printf("\n  %-24s %14s %22s %14s\n", "", "total time", "per prime (mean)", "median");
    std::printf("  %-24s %14s %22s %14s\n", "A: method of the paper", fmt_time(sumA).c_str(), fmt_time(sumA / k).c_str(), fmt_time(tools::median(TA)).c_str());
    std::printf("  %-24s %14s %22s %14s\n", "B: plain Miller-Rabin", fmt_time(sumB).c_str(), fmt_time(sumB / k).c_str(), fmt_time(tools::median(TB)).c_str());
    std::printf("  Speedup (B / A): %.2fx\n", sumB / sumA);
    if (k == 1) std::printf("------------------------------------------------------------------------------------------------\nPrime:\n%s\n", primes[0].c_str());
    if (save) {
        std::string fname = "primes_" + std::to_string(bits) + "bit_" + std::to_string(k) + "_seed" + std::to_string(seed) + ".txt";
        std::ofstream f(fname);
        f << "# " << k << " primes of " << bits << " bits, seed = " << seed << ", MR rounds = " << rounds << "\n";
        for (const auto& p : primes) f << p << "\n";
        std::printf("Primes saved to: %s\n", fname.c_str());
    }
    std::printf("================================================================================================\n");
    return all_same ? 0 : 2;
}
