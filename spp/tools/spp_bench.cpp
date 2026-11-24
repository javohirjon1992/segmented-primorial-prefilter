// spp_bench — repeats the prefilter measurements of the paper and calibrates B* for the current machine.
//
//   spp_bench 1024              (10000 candidates)
//   spp_bench 2048 4000 --sweep (sweep over B: find B*)
#include "common.hpp"
#include "spp/spp.hpp"

#include <cmath>
#include <cstdlib>
#include <cstring>

using namespace spp;
using tools::fmt_time;

template <class F>
static double per_candidate_ns(const F& filter, const std::vector<BigInt>& c, int reps) {
    std::vector<double> t;
    volatile unsigned long sink = 0;
    for (int r = 0; r < reps; ++r) {
        double t0 = tools::now_ns();
        for (const auto& n : c) sink += filter.passes(n.get());
        t.push_back((tools::now_ns() - t0) / c.size());
    }
    (void)sink;
    return tools::median(t);
}

static double survival(const std::vector<unsigned long>& P) {
    double s = 1;
    for (unsigned long p : P) s *= 1.0 - 1.0 / p;
    return s;
}

int main(int argc, char** argv) {
    if (argc < 2 || !std::strcmp(argv[1], "-h") || !std::strcmp(argv[1], "--help")) {
        std::printf("Usage: spp_bench <bit length> [number of candidates=10000] [--sweep] [--seed S] [--reps R]\n");
        return argc < 2 ? 1 : 0;
    }
    unsigned bits = static_cast<unsigned>(std::strtoul(argv[1], nullptr, 10));
    std::size_t count = 10000;
    bool sweep = false;
    std::uint64_t seed = 20261101;
    int reps = 11;
    for (int i = 2; i < argc; ++i) {
        if (!std::strcmp(argv[i], "--sweep")) sweep = true;
        else if (!std::strcmp(argv[i], "--seed") && i + 1 < argc) seed = std::strtoull(argv[++i], nullptr, 10);
        else if (!std::strcmp(argv[i], "--reps") && i + 1 < argc) reps = std::atoi(argv[++i]);
        else count = std::strtoull(argv[i], nullptr, 10);
    }
    const int rounds = kDefaultRounds;
    SeededRandom rng(seed);
    std::vector<BigInt> cand(count);
    for (auto& n : cand) random_odd(n.get(), bits, rng);
    unsigned long Bstar = recommended_bound(bits);

    std::printf("================================================================================\n");
    std::printf("spp_bench: %zu random odd %u-bit candidates | GMP %s | spp %s | seed %llu\n", count, bits,
                gmp_version, version(), static_cast<unsigned long long>(seed));

    // --- the Miller-Rabin time of every candidate is measured once (per-candidate measurement)
    std::vector<double> mr_ns(count);
    std::vector<int> mr_rounds(count);
    std::vector<char> is_prime(count);
    double comp_sum = 0; std::size_t comp_n = 0;
    for (std::size_t i = 0; i < count; ++i) {
        double t0 = tools::now_ns();
        MRResult m = miller_rabin(cand[i].get(), rounds);
        mr_ns[i] = tools::now_ns() - t0;
        mr_rounds[i] = m.rounds; is_prime[i] = m.probable_prime;
        if (!m.probable_prime) { comp_sum += mr_ns[i]; ++comp_n; }
    }
    double c_comp = comp_n ? comp_sum / comp_n : 0;         // composite MR call (usually one round)
    std::printf("MR per composite candidate: %s\n", fmt_time(c_comp).c_str());

    // --- prefilter cost (median, ns per candidate)
    TrialDivision td1(1000), tdS(Bstar);
    SegmentedPrefilter g1(1000), gS(Bstar);
    struct Row { const char* name; double pre; std::vector<char> pass; unsigned long B; };
    auto passes_of = [&](const auto& f) { std::vector<char> v(count); for (std::size_t i = 0; i < count; ++i) v[i] = f.passes(cand[i].get()); return v; };
    std::vector<Row> rows = {
        {"TD (B = 1000)", per_candidate_ns(td1, cand, reps), passes_of(td1), 1000},
        {"GEOM (B = 1000)", per_candidate_ns(g1, cand, reps), passes_of(g1), 1000},
        {"TD (B = B*)", per_candidate_ns(tdS, cand, reps), passes_of(tdS), Bstar},
        {"GEOM (B = B*)", per_candidate_ns(gS, cand, reps), passes_of(gS), Bstar},
    };
    bool equal1 = rows[0].pass == rows[1].pass, equalS = rows[2].pass == rows[3].pass;
    std::printf("Identical prefilter decisions: at B = 1000 %s, at B* = %lu %s\n", equal1 ? "yes" : "NO", Bstar, equalS ? "yes" : "NO");
    std::printf("--------------------------------------------------------------------------------\n");
    std::printf("%-18s %12s %12s %12s %14s %10s %10s\n", "Method", "prefilter", "passed", "W(B)", "per candidate", "vs. MR", "vs. TD");
    double mr_only = 0; long long mr_only_rounds = 0;
    for (std::size_t i = 0; i < count; ++i) { mr_only += mr_ns[i]; mr_only_rounds += mr_rounds[i]; }
    mr_only /= count;
    std::printf("%-18s %12s %12s %12s %14s %9.2fx %10s\n", "MR only", "-", "100%", "-", fmt_time(mr_only).c_str(), 1.0, "");
    double td_total = 0;
    for (std::size_t r = 0; r < rows.size(); ++r) {
        double mr = 0; std::size_t passed = 0;
        for (std::size_t i = 0; i < count; ++i) if (rows[r].pass[i]) { mr += mr_ns[i]; ++passed; }
        double total = rows[r].pre + mr / count;
        if (r == 0) td_total = total;
        std::printf("%-18s %12s %11.2f%% %12.4f %14s %9.2fx %9.2fx\n", rows[r].name, fmt_time(rows[r].pre).c_str(),
                    100.0 * passed / count, survival(odd_primes_upto(rows[r].B)), fmt_time(total).c_str(),
                    mr_only / total, td_total / total);
    }
    std::printf("(B* = %lu; speedups refer to this set of candidates; primes in the sample: %d)\n", Bstar,
                static_cast<int>(std::count(is_prime.begin(), is_prime.end(), 1)));

    if (sweep) {
        std::printf("--------------------------------------------------------------------------------\n");
        std::printf("Sweep over B (segmented, S0 = 5, r = 2): expected cost = prefilter + pass fraction x MR\n");
        std::printf("%10s %9s %12s %12s %18s\n", "B", "segments", "prefilter", "passed", "expected cost");
        const unsigned long grid[] = {300, 1000, 3000, 10000, 30000, 100000, 300000, 1000000, 3000000};
        double best = 1e300; unsigned long bestB = 0;
        for (unsigned long B : grid) {
            SegmentedPrefilter pf(B);
            double pre = per_candidate_ns(pf, cand, std::max(3, reps / 2));
            std::size_t passed = 0;
            for (const auto& n : cand) passed += pf.passes(n.get());
            double ratio = double(passed) / count;
            double cost = pre + ratio * c_comp;
            if (cost < best) { best = cost; bestB = B; }
            std::printf("%10lu %9zu %12s %11.2f%% %18s\n", B, pf.segment_count(), fmt_time(pre).c_str(), 100 * ratio, fmt_time(cost).c_str());
        }
        std::printf("B* for this machine = %lu (recommended_bound() returns %lu)\n", bestB, Bstar);
    }
    std::printf("================================================================================\n");
    return (equal1 && equalS) ? 0 : 2;
}
