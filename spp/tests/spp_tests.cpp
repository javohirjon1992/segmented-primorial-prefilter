// Unit tests for spp. Exit code 0 = all tests passed.
#include "spp/spp.hpp"

#include <algorithm>
#include <cstdio>
#include <string>
#include <vector>

static int failures = 0, checks = 0;
#define CHECK(cond) do { ++checks; if (!(cond)) { ++failures; std::printf("FAIL %s:%d: %s\n", __FILE__, __LINE__, #cond); } } while (0)

using namespace spp;

static void test_primes() {
    CHECK(odd_primes_upto(1000).size() == 167);
    CHECK(odd_primes_upto(10000).size() == 1228);
    CHECK(odd_primes_upto(2).empty());
    SegmentedPrefilter pf(1000);
    CHECK(pf.prime_count() == 167);
    CHECK(pf.largest_prime() == 997);
    // geometric sizes 5,10,20,40,80,12 and bit lengths of the segment-survival table of the paper
    std::vector<std::size_t> sz = pf.segment_sizes(), bt = pf.segment_bits();
    CHECK((sz == std::vector<std::size_t>{5, 10, 20, 40, 80, 12}));
    CHECK((bt == std::vector<std::size_t>{14, 50, 133, 321, 744, 119}));
}

static void test_equivalence() {
    SeededRandom rng(20261101);
    BigInt n;
    const unsigned long bounds[] = {300, 1000, 10000, 100000};
    const unsigned s0s[] = {1, 5, 16};
    const double ratios[] = {1.0, 1.5, 2.0, 4.0};
    const unsigned bitlens[] = {64, 512, 1024, 2048};
    for (unsigned long B : bounds) {
        TrialDivision td(B);
        for (unsigned s0 : s0s)
            for (double r : ratios) {
                SegmentedPrefilter pf(B, s0, r);
                for (unsigned bits : bitlens)
                    for (int i = 0; i < 300; ++i) {
                        random_odd(n.get(), bits, rng);
                        bool a = pf.passes(n.get()), b = td.passes(n.get());
                        CHECK(a == b);
                        FilterResult f = pf.check(n.get());
                        CHECK(f.passed == a);
                        if (!f.passed) CHECK(f.factor > 2 && mpz_divisible_ui_p(n.get(), f.factor));
                    }
            }
    }
}

static void test_miller_rabin() {
    const char* composites[] = {"561", "1105", "1729", "2465", "2821", "6601", "8911", "25326001", "3215031751",
                                "2152302898747", "3474749660383", "341550071728321", "3825123056546413051",
                                "318665857834031151167461",
                                // RSA-100 (product of two 50-digit primes)
                                "1522605027922533360535618378132637429718068114961380688657908494580122963258952897654000350692006139"};
    for (const char* c : composites) {
        BigInt n(c);
        CHECK(!miller_rabin(n.get()).probable_prime);
        CHECK(!is_probable_prime(n.get()).probable_prime);
        SystemRandom sr;
        CHECK(!miller_rabin(n.get(), 40, BaseMode::Random, &sr).probable_prime);
    }
    const char* primes[] = {"2^89-1", "2^107-1", "2^127-1", "2^521-1", "2^607-1", "2^1279-1",
                            "37975227936943673922808872755445627854565536638199"};
    for (const char* p : primes) {
        BigInt n = parse_expression(p);
        CHECK(miller_rabin(n.get()).probable_prime);
        CHECK(is_probable_prime(n.get()).probable_prime);
        CHECK(is_probable_prime_plain(n.get()).probable_prime);
    }
    // small numbers against the sieve (covers the n <= p_k branch of Lemma 1)
    std::vector<unsigned long> P = odd_primes_upto(5000);
    SegmentedPrefilter pf(1000);
    BigInt n;
    for (unsigned long v = 0; v <= 5000; ++v) {
        bool expect = (v == 2) || std::binary_search(P.begin(), P.end(), v);
        mpz_set_ui(n.get(), v);
        CHECK(is_probable_prime(n.get(), pf).probable_prime == expect);
        CHECK(is_probable_prime_plain(n.get()).probable_prime == expect);
    }
}

static void test_generation() {
    SegmentedPrefilter pf(recommended_bound(512));
    BigInt a, b, c;
    for (std::uint64_t seed = 1; seed <= 5; ++seed) {
        SeededRandom r1(seed), r2(seed);
        GenStats s1, s2;
        generate_prime(a.get(), 512, pf, r1, 40, &s1);
        generate_prime_plain(b.get(), 512, r2, 40, &s2);
        CHECK(mpz_cmp(a.get(), b.get()) == 0);          // identical candidate streams -> identical prime
        CHECK(s1.candidates == s2.candidates);
        CHECK(s1.mr_calls + s1.prefilter_rejected == s1.candidates);
        CHECK(a.bits() == 512);
        CHECK(miller_rabin(a.get(), 64).probable_prime);
        SeededRandom r3(seed);
        generate_prime_incremental(c.get(), 512, pf, r3);
        CHECK(c.bits() == 512 && miller_rabin(c.get(), 64).probable_prime);
    }
    SystemRandom sys;
    generate_prime(a.get(), 256, SegmentedPrefilter(1000), sys);
    CHECK(a.bits() == 256);
}

static void test_parser() {
    CHECK(parse_expression("2^10").str() == "1024");
    CHECK(parse_expression("2**10 - 24").str() == "1000");
    CHECK(parse_expression("(2^3+1)*(10-3)").str() == "63");
    CHECK(parse_expression("0xff").str() == "255");
    CHECK(parse_expression("1_000_007").str() == "1000007");
    bool threw = false;
    try { parse_expression("2^^3"); } catch (const std::exception&) { threw = true; }
    CHECK(threw);
}

int main() {
    test_primes();
    test_equivalence();
    test_miller_rabin();
    test_generation();
    test_parser();
    // ---- tests of version 1.1
    {
        CHECK(spp::mr_rounds_for_error(1024) == 6); CHECK(spp::mr_rounds_for_error(1536) == 4); CHECK(spp::mr_rounds_for_error(2048) == 3);
        // table of the paper: RSA pair (S3, eps/2) needs 8, 6, 5, 4, 3, 3 rounds at 768..3072 bits; 13 at 512 bits because of L(S)
        {   const unsigned bb[] = {768, 1024, 1280, 1536, 2048, 3072}; const int tt[] = {8, 6, 5, 4, 3, 3};
            for (int i = 0; i < 6; ++i) CHECK(spp::mr_rounds_for_error(bb[i], -129.0, spp::loss_bound(spp::CandidateSet::S3_RsaCoprime)) == tt[i]);
            CHECK(spp::mr_rounds_for_error(512, -129.0, 2.045) == 13); CHECK(spp::mr_rounds_for_error(512, -128.0, 1.0) == 12); }
        for (unsigned bits : {33u, 35u, 47u, 65u, 129u}) for (int k = 0; k < 4; ++k) {      // odd b: condition 2q + 1 > sqrt(p) of Lemma 2
            spp::SeededRandom r(500 + bits + k); spp::BigInt p; spp::generate_provable_prime(p.get(), bits, spp::SegmentedPrefilter(1000), r);
            CHECK(mpz_sizeinbase(p.get(), 2) == bits); CHECK(mpz_probab_prime_p(p.get(), 30) > 0); }
        for (auto s : {spp::CandidateSet::S1_Odd, spp::CandidateSet::S2_TopTwo, spp::CandidateSet::S3_RsaCoprime}) for (int k = 0; k < 3; ++k) {
            spp::SeededRandom r(700 + k); spp::BigInt p; spp::gen_prime(p.get(), 768, s, -128.0, spp::SegmentedPrefilter(30000), r);
            CHECK(mpz_sizeinbase(p.get(), 2) == 768); CHECK(mpz_probab_prime_p(p.get(), 30) > 0);
            if (s != spp::CandidateSet::S1_Odd) CHECK(mpz_tstbit(p.get(), 766)); }
        spp::SegmentedPrefilter pf(100000);
        for (unsigned bits : {40u, 100u, 256u, 512u, 1024u}) for (int k = 0; k < 4; ++k) {
            spp::SeededRandom r(900 + bits + k); spp::BigInt p; spp::ProvableStats st;
            spp::generate_provable_prime(p.get(), bits, pf, r, &st);
            CHECK(mpz_sizeinbase(p.get(), 2) == bits); CHECK(spp::miller_rabin(p.get(), 30).probable_prime);
        }
        for (unsigned nlen : {1024u, 2048u}) for (int k = 0; k < 3; ++k) {
            spp::SeededRandom r(1900 + nlen + k); spp::BigInt p, q, d, g, pm1, e(65537UL);
            spp::generate_rsa_primes_fips(p.get(), q.get(), nlen, 65537, pf, r);
            unsigned hb = nlen / 2;
            CHECK(mpz_sizeinbase(p.get(), 2) == hb && mpz_tstbit(p.get(), hb - 2)); CHECK(mpz_sizeinbase(q.get(), 2) == hb && mpz_tstbit(q.get(), hb - 2));
            mpz_sub_ui(pm1.get(), p.get(), 1); mpz_gcd(g.get(), pm1.get(), e.get()); CHECK(mpz_cmp_ui(g.get(), 1) == 0);
            CHECK(spp::miller_rabin(p.get(), 30).probable_prime && spp::miller_rabin(q.get(), 30).probable_prime);
        }
    }
    std::printf("spp %s: %d checks, %d failures\n", version(), checks, failures);
    
    return failures == 0 ? 0 : 1;
}
