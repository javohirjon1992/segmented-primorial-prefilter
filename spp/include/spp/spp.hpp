// spp — Segmented Primorial Prefilter for Miller–Rabin primality testing.
//
// C++17 library built on GMP. It implements the method of the paper
// "Segmented Primorial Prefiltering and Cost-Optimal Sieving Bounds for Prime Generation":
//   * odd primes 3..B are split into segments whose sizes grow geometrically (S0, S0*r, S0*r^2, ...);
//   * every segment is stored as its primorial M_i;
//   * a candidate n is rejected as soon as gcd(n, M_i) > 1 for some segment;
//   * survivors go to the Miller–Rabin test.
// For n > p_k (largest sieving prime) the prefilter rejects exactly the candidates rejected by trial
// division over the same primes (Lemma 1 of the paper).
//
// Thread safety: SegmentedPrefilter and TrialDivision are immutable after construction; their const
// member functions may be called concurrently from several threads.

#ifndef SPP_SPP_HPP
#define SPP_SPP_HPP

#include <gmp.h>

#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

namespace spp {

constexpr int kDefaultRounds = 40;     // Miller–Rabin rounds (error probability <= 4^-40)
constexpr unsigned kDefaultS0 = 5;     // initial segment size (Section IV-C of the paper)
constexpr double kDefaultRatio = 2.0;  // geometric growth ratio r

// ----------------------------------------------------------------------------- BigInt
// Small RAII wrapper around mpz_t. Use .get() to pass it to GMP or to spp functions.
class BigInt {
public:
    BigInt();
    explicit BigInt(unsigned long v);
    explicit BigInt(const std::string& decimal_or_hex);  // "123", "0x7f", base detected automatically
    BigInt(const BigInt& o);
    BigInt(BigInt&& o) noexcept;
    BigInt& operator=(const BigInt& o);
    BigInt& operator=(BigInt&& o) noexcept;
    ~BigInt();

    mpz_ptr get() { return v_; }
    mpz_srcptr get() const { return v_; }
    std::size_t bits() const;
    std::string str(int base = 10) const;

private:
    mpz_t v_;
};

// ----------------------------------------------------------------------------- random sources
// Source of random bytes for candidates and Miller–Rabin bases.
class RandomSource {
public:
    virtual ~RandomSource() = default;
    virtual void fill(std::uint8_t* out, std::size_t len) = 0;
    std::uint64_t next_u64();
};

// Operating-system CSPRNG (BCryptGenRandom on Windows, getrandom()/dev/urandom elsewhere).
// Use this for real keys.
class SystemRandom final : public RandomSource {
public:
    void fill(std::uint8_t* out, std::size_t len) override;
};

// Deterministic, reproducible generator (xoshiro256**). NOT cryptographically secure:
// intended for benchmarks and experiments where the same candidate sequence must be repeated.
class SeededRandom final : public RandomSource {
public:
    explicit SeededRandom(std::uint64_t seed);
    void fill(std::uint8_t* out, std::size_t len) override;

private:
    std::uint64_t s_[4];
    std::uint64_t next();
};

// ----------------------------------------------------------------------------- prefilters
// Default sieving bound for a candidate bit length (near the best bounds of the paper; calibrate per platform with spp_bench --sweep).
unsigned long recommended_bound(unsigned bits);

// Odd primes 3..B (sieve of Eratosthenes).
std::vector<unsigned long> odd_primes_upto(unsigned long bound);

struct FilterResult {
    bool passed = true;          // true: no small prime factor found
    int segment = -1;            // 1-based index of the segment that rejected n (or -1)
    unsigned long factor = 0;    // a small prime factor of n (0 if passed)
};

class SegmentedPrefilter {
public:
    // bound: sieving bound B (odd primes 3..B are used); s0, ratio: geometric segmentation parameters.
    explicit SegmentedPrefilter(unsigned long bound, unsigned s0 = kDefaultS0, double ratio = kDefaultRatio);
    ~SegmentedPrefilter();
    SegmentedPrefilter(const SegmentedPrefilter&) = delete;
    SegmentedPrefilter& operator=(const SegmentedPrefilter&) = delete;
    SegmentedPrefilter(SegmentedPrefilter&&) noexcept;
    SegmentedPrefilter& operator=(SegmentedPrefilter&&) noexcept;

    // Fast path: true if n has no prime factor in 3..B. Precondition: n odd and n > largest_prime().
    bool passes(mpz_srcptr n) const;
    // Detailed version (also reports the segment and a small factor).
    FilterResult check(mpz_srcptr n) const;

    unsigned long bound() const;
    unsigned long largest_prime() const;
    std::size_t prime_count() const;
    std::size_t segment_count() const;
    // sizes (number of primes) and bit lengths of the segment primorials
    std::vector<std::size_t> segment_sizes() const;
    std::vector<std::size_t> segment_bits() const;
    const std::vector<unsigned long>& primes() const;
    double build_seconds() const;   // time spent constructing the segments

private:
    struct Impl;
    std::unique_ptr<Impl> p_;
};

// Classical trial division over the same primes (baseline for comparisons).
class TrialDivision {
public:
    explicit TrialDivision(unsigned long bound);
    bool passes(mpz_srcptr n) const;
    FilterResult check(mpz_srcptr n) const;
    unsigned long largest_prime() const;
    std::size_t prime_count() const;

private:
    std::vector<unsigned long> primes_;
};

// ----------------------------------------------------------------------------- Miller–Rabin
enum class BaseMode {
    Random,          // bases from a RandomSource (recommended for numbers supplied by others)
    FromCandidate    // bases derived deterministically from n: identical work for identical n
                     // (reproducible comparisons; do not use on adversarial inputs, cf. Albrecht et al., CCS 2018)
};

struct MRResult {
    bool probable_prime = false;
    int rounds = 0;               // Miller–Rabin rounds actually executed
};

// Miller–Rabin test with t rounds. Handles n < 2, n = 2, 3 and even n.
MRResult miller_rabin(mpz_srcptr n, int rounds = kDefaultRounds, BaseMode mode = BaseMode::FromCandidate,
                      RandomSource* rng = nullptr);

// ----------------------------------------------------------------------------- full test
enum class Stage { Trivial, SmallSieve, Prefilter, MillerRabin };

struct TestReport {
    bool probable_prime = false;
    Stage stage = Stage::Trivial;  // where the decision was made
    int segment = -1;              // rejecting segment (Stage::Prefilter)
    unsigned long factor = 0;      // small factor found (if any)
    int mr_rounds = 0;
};

// Method of the paper: segmented primorial prefilter followed by Miller–Rabin.
TestReport is_probable_prime(mpz_srcptr n, const SegmentedPrefilter& pf, int rounds = kDefaultRounds,
                             BaseMode mode = BaseMode::FromCandidate, RandomSource* rng = nullptr);
// Convenience overload: chooses B from the bit length and caches the prefilter internally.
TestReport is_probable_prime(mpz_srcptr n, int rounds = kDefaultRounds,
                             BaseMode mode = BaseMode::FromCandidate, RandomSource* rng = nullptr);
// Baseline: Miller–Rabin alone.
TestReport is_probable_prime_plain(mpz_srcptr n, int rounds = kDefaultRounds,
                                   BaseMode mode = BaseMode::FromCandidate, RandomSource* rng = nullptr);

// ----------------------------------------------------------------------------- prime generation
struct GenStats {
    std::uint64_t candidates = 0;          // candidates examined
    std::uint64_t prefilter_rejected = 0;  // rejected without modular exponentiation
    std::uint64_t mr_calls = 0;            // Miller–Rabin tests started
    std::uint64_t mr_rounds = 0;           // Miller–Rabin rounds executed
    void add(const GenStats& o);
};

// Uniformly random odd integer with exactly `bits` bits (top bit and lowest bit set).
void random_odd(mpz_ptr out, unsigned bits, RandomSource& rng);

// Random probable prime of exactly `bits` bits using the segmented prefilter.
void generate_prime(mpz_ptr out, unsigned bits, const SegmentedPrefilter& pf, RandomSource& rng,
                    int rounds = kDefaultRounds, GenStats* stats = nullptr);
// Same candidate stream, Miller–Rabin only (baseline). With equal RNG state both functions return the same prime.
void generate_prime_plain(mpz_ptr out, unsigned bits, RandomSource& rng,
                          int rounds = kDefaultRounds, GenStats* stats = nullptr);
// Incremental search n0, n0+2, ... from a random odd start, using the prefilter.
void generate_prime_incremental(mpz_ptr out, unsigned bits, const SegmentedPrefilter& pf, RandomSource& rng,
                                int rounds = kDefaultRounds, GenStats* stats = nullptr);

// Parses decimal/hex integers and simple expressions with + - * ^ ( ), e.g. "2^521-1", "(2^127-1)*(2^89-1)".
BigInt parse_expression(const std::string& text);

const char* version();


// ---- 1.1: complete generators --------------------------------------------------
// Admissible candidate sets (Proposition 2 of the paper): S1 — all odd b-bit integers (L = 1); S2 — the two top bits set
// (L <= 2.004); S3 — S2 with gcd(e, n-1) = 1 for a prime e >= 65537 (L <= 2.045, 512 <= b <= 4096). No guarantee is given
// for other restrictions.
enum class CandidateSet { S1_Odd, S2_TopTwo, S3_RsaCoprime };
double loss_bound(CandidateSet s);

// Smallest t >= 3 (and t <= bits/9) with L * p_{bits,t} / (1 - p_{bits,t}) <= 2^log2_eps, where p_{bits,t} is the
// Damgard-Landrock-Pomerance bound and loss = L(S); for each prime of an RSA pair pass log2_eps - 1 (union bound).
int mr_rounds_for_error(unsigned bits, double log2_eps = -128.0, double loss = 1.0);

// GenPrime (Algorithm 2): independent candidates from S, the segmented prefilter, and t Miller-Rabin rounds with
// INDEPENDENT RANDOM bases drawn from rng. Probability of a composite output <= 2^log2_eps (average case, Proposition 2).
// Use SystemRandom for real keys.
void gen_prime(mpz_ptr out, unsigned bits, CandidateSet s, double log2_eps, const SegmentedPrefilter& pf, RandomSource& rng,
               unsigned long e = 65537, GenStats* stats = nullptr);

struct ProvableStats { unsigned long candidates = 0, prefilter_rejected = 0, modexps = 0, levels = 0;
                       unsigned long q_below_sqrt = 0;   // levels with q^2 < p, where only 2q + 1 > sqrt(p) holds
};
// Provable prime: p = 2Rq + 1, q a recursively proven prime of floor(b/2)+1 bits, so 2q + 1 > sqrt(p) for every b;
// hence 2^(p-1) = 1 (mod p) and gcd(2^(2R) - 1, p) = 1 imply that p is prime (Pocklington, Lemma 2 of the paper).
// Error probability 0; the prefilter is applied at every level.
void generate_provable_prime(mpz_ptr out, unsigned bits, const SegmentedPrefilter& pf, RandomSource& rng, ProvableStats* stats = nullptr);

// Pair of RSA primes with the FIPS 186-5 conditions: each prime from gen_prime(S3, log2_eps - 1), |p - q| > 2^(nlen/2 - 100).
// Error probability of the pair <= 2^log2_eps (union bound).
void generate_rsa_primes_fips(mpz_ptr p, mpz_ptr q, unsigned nlen, unsigned long e, const SegmentedPrefilter& pf, RandomSource& rng,
                              double log2_eps = -128.0);
}  // namespace spp

#endif  // SPP_SPP_HPP
