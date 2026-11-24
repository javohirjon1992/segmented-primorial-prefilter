// Minimal example of the spp API.
#include <spp/spp.hpp>
#include <cstdio>

int main() {
    // 1) primality test with the method of the paper (B chosen automatically from the bit length)
    spp::BigInt n = spp::parse_expression("2^607-1");
    spp::TestReport r = spp::is_probable_prime(n.get());
    std::printf("2^607-1: %s (stage %d, %d MR rounds)\n", r.probable_prime ? "probable prime" : "composite",
                static_cast<int>(r.stage), r.mr_rounds);

    // 2) numbers received from others: use random Miller-Rabin bases from the OS CSPRNG
    spp::SystemRandom os_rng;
    spp::BigInt m("1522605027922533360535618378132637429718068114961380688657908494580122963258952897654000350692006139");
    r = spp::is_probable_prime(m.get(), 40, spp::BaseMode::Random, &os_rng);
    std::printf("RSA-100: %s\n", r.probable_prime ? "probable prime" : "composite");

    // 3) GenPrime: a 2048-bit prime from S2 with error probability below 2^-128 (independent candidates, bases
    //    drawn from the OS CSPRNG); the prefilter is built once and can be reused (and shared between threads)
    spp::SegmentedPrefilter pf(spp::recommended_bound(2048));
    spp::BigInt p;
    spp::GenStats st;
    spp::gen_prime(p.get(), 2048, spp::CandidateSet::S2_TopTwo, -128.0, pf, os_rng, 65537, &st);
    std::printf("2048-bit prime found after %llu candidates (%llu rejected by the prefilter, %llu Miller-Rabin calls)\n",
                (unsigned long long)st.candidates, (unsigned long long)st.prefilter_rejected, (unsigned long long)st.mr_calls);
    std::printf("p = %s...\n", p.str().substr(0, 40).c_str());
    return 0;
}
