// spp_check — tests a number with the method of the paper (segmented primorial prefilter + Miller–Rabin)
// and with the plain Miller–Rabin test, and reports the time taken.
//
//   spp_check 2^521-1
//   spp_check "(2^127-1)*(2^89-1)" --rounds 40
//   spp_check                       (interactive mode)
#include "common.hpp"
#include "spp/spp.hpp"

#include <cstdlib>
#include <cstring>
#include <iostream>
#include <map>
#include <memory>

using namespace spp;
using tools::fmt_time;

static const char* stage_name(Stage s) {
    switch (s) {
        case Stage::Trivial: return "trivial case (n < 2, n = 2, or even n)";
        case Stage::SmallSieve: return "n <= B: from the sieve table";
        case Stage::Prefilter: return "prefilter: small divisor found";
        case Stage::MillerRabin: return "passed the prefilter -> Miller-Rabin";
    }
    return "";
}

static void check(const std::string& text, int rounds) {
    BigInt n = parse_expression(text);
    unsigned bits = static_cast<unsigned>(std::max<std::size_t>(1, n.bits()));
    unsigned long B = recommended_bound(bits);
    static std::map<unsigned long, std::unique_ptr<SegmentedPrefilter>> cache;
    bool fresh = !cache.count(B);
    if (fresh) cache[B].reset(new SegmentedPrefilter(B));
    const SegmentedPrefilter& pf = *cache[B];

    std::printf("==============================================================================\n");
    std::printf("n = %s\n", tools::shorten(n.str()).c_str());
    std::printf("Bit length: %u   |   GMP %s   |   spp %s\n", bits, gmp_version, version());
    if (fresh)
        std::printf("Segments prepared: B = %lu, %zu primes, %zu segments (%s, once)\n",
                    B, pf.prime_count(), pf.segment_count(), fmt_time(pf.build_seconds() * 1e9).c_str());
    std::printf("------------------------------------------------------------------------------\n");

    // alternating measurements: A, B, A, B, ...
    std::vector<double> ta, tb;
    TestReport ra, rb;
    double total = 0;
    while ((total < 6e8 || ta.size() < 5) && ta.size() < 20000) {
        double t0 = tools::now_ns(); ra = is_probable_prime(n.get(), pf, rounds);
        double t1 = tools::now_ns(); rb = is_probable_prime_plain(n.get(), rounds);
        double t2 = tools::now_ns();
        ta.push_back(t1 - t0); tb.push_back(t2 - t1); total += t2 - t0;
        if (ta.size() >= 5 && total > 3e9) break;
    }
    double ma = tools::median(ta), mb = tools::median(tb);
    double pre = 0;
    if (ra.stage == Stage::Prefilter || ra.stage == Stage::MillerRabin) {
        const int reps = 2000;
        double t0 = tools::now_ns();
        for (int i = 0; i < reps; ++i) (void)pf.passes(n.get());
        pre = (tools::now_ns() - t0) / reps;
    }
    auto verdict = [](bool p) { return p ? "PROBABLE PRIME" : "COMPOSITE"; };
    std::printf("1) Method of the paper: segmented primorial prefilter + Miller-Rabin\n");
    std::printf("   Stage   : %s", stage_name(ra.stage));
    if (ra.segment > 0) std::printf(" (segment %d)", ra.segment);
    std::printf("\n");
    if (ra.factor) std::printf("   Divisor : %lu\n", ra.factor);
    std::printf("   MR rounds: %d\n   Result   : %s\n", ra.mr_rounds, verdict(ra.probable_prime));
    std::printf("   Time     : %s  (median of %zu runs)\n", fmt_time(ma).c_str(), ta.size());
    if (pre > 0) std::printf("              of which prefilter: %s\n", fmt_time(pre).c_str());
    std::printf("\n2) Plain Miller-Rabin test (up to %d rounds)\n", rounds);
    std::printf("   MR rounds: %d\n   Result   : %s\n", rb.mr_rounds, verdict(rb.probable_prime));
    std::printf("   Time     : %s  (median of %zu runs)\n", fmt_time(mb).c_str(), tb.size());
    std::printf("------------------------------------------------------------------------------\n");
    if (ra.probable_prime != rb.probable_prime) {
        std::printf("WARNING: the results disagree!\n");
    } else {
        double r = mb / ma;
        const char* who = r > 1.05 ? "method of the paper faster" : (r < 0.95 ? "plain MR faster" : "about equal");
        std::printf("Results agree. Speedup: %.2fx (%s)\n", r, who);
        if (ra.stage == Stage::MillerRabin && !ra.probable_prime)
            std::printf("Note: this number has no divisor up to %lu, so the prefilter only adds cost.\n", B);
    }
    std::printf("==============================================================================\n");
}

int main(int argc, char** argv) {
    int rounds = kDefaultRounds;
    std::string expr;
    for (int i = 1; i < argc; ++i) {
        if (!std::strcmp(argv[i], "--rounds") && i + 1 < argc) rounds = std::atoi(argv[++i]);
        else if (!std::strcmp(argv[i], "-h") || !std::strcmp(argv[i], "--help")) {
            std::printf("Usage: spp_check [number or expression] [--rounds t]\nExample: spp_check 2^1279-1\n");
            return 0;
        } else { if (!expr.empty()) expr += ' '; expr += argv[i]; }
    }
    try {
        if (!expr.empty()) { check(expr, rounds); return 0; }
        std::printf("Primality check (method of the paper and plain Miller-Rabin). Quit: q\n");
        std::string line;
        while (std::printf("\nEnter a number: "), std::fflush(stdout), std::getline(std::cin, line)) {
            if (line == "q" || line == "quit" || line == "exit") break;
            if (line.find_first_not_of(" \t\r") == std::string::npos) continue;
            try { check(line, rounds); } catch (const std::exception& e) { std::printf("Error: %s\n", e.what()); }
        }
    } catch (const std::exception& e) {
        std::fprintf(stderr, "Error: %s\n", e.what());
        return 1;
    }
    return 0;
}
