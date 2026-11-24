// spp — Segmented Primorial Prefilter for Miller–Rabin primality testing (implementation).
#include "spp/spp.hpp"

#include <algorithm>
#include <cctype>
#include <chrono>
#include <climits>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <limits>
#include <map>
#include <mutex>
#include <stdexcept>

#if defined(_WIN32)
#  ifndef NOMINMAX
#    define NOMINMAX
#  endif
#  include <windows.h>
#  include <bcrypt.h>
#elif defined(__linux__)
#  include <sys/random.h>
#  include <cerrno>
#elif defined(__APPLE__) || defined(__FreeBSD__) || defined(__OpenBSD__)
#  include <unistd.h>
#endif

namespace spp {

const char* version() { return "1.1.1"; }

// ============================================================================ BigInt
BigInt::BigInt() { mpz_init(v_); }
BigInt::BigInt(unsigned long v) { mpz_init_set_ui(v_, v); }
BigInt::BigInt(const std::string& s) {
    mpz_init(v_);
    if (mpz_set_str(v_, s.c_str(), 0) != 0) {
        mpz_clear(v_);
        throw std::invalid_argument("spp::BigInt: not an integer: " + s);
    }
}
BigInt::BigInt(const BigInt& o) { mpz_init(v_); mpz_set(v_, o.v_); }
BigInt::BigInt(BigInt&& o) noexcept { mpz_init(v_); std::swap(v_[0], o.v_[0]); }
BigInt& BigInt::operator=(const BigInt& o) { if (this != &o) mpz_set(v_, o.v_); return *this; }
BigInt& BigInt::operator=(BigInt&& o) noexcept { std::swap(v_[0], o.v_[0]); return *this; }
BigInt::~BigInt() { mpz_clear(v_); }
std::size_t BigInt::bits() const { return mpz_cmp_ui(v_, 0) == 0 ? 0 : mpz_sizeinbase(v_, 2); }
std::string BigInt::str(int base) const {
    std::string buf(mpz_sizeinbase(v_, base) + 2, '\0');
    mpz_get_str(&buf[0], base, v_);
    buf.resize(std::strlen(buf.c_str()));
    return buf;
}

// ============================================================================ random sources
std::uint64_t RandomSource::next_u64() {
    std::uint8_t b[8];
    fill(b, 8);
    std::uint64_t v = 0;
    for (int i = 0; i < 8; ++i) v = (v << 8) | b[i];
    return v;
}

void SystemRandom::fill(std::uint8_t* out, std::size_t len) {
#if defined(_WIN32)
    if (BCryptGenRandom(nullptr, out, static_cast<ULONG>(len), BCRYPT_USE_SYSTEM_PREFERRED_RNG) != 0)
        throw std::runtime_error("spp::SystemRandom: BCryptGenRandom failed");
#elif defined(__linux__)
    std::size_t done = 0;
    while (done < len) {
        ssize_t r = getrandom(out + done, len - done, 0);
        if (r < 0) {
            if (errno == EINTR) continue;
            throw std::runtime_error("spp::SystemRandom: getrandom failed");
        }
        done += static_cast<std::size_t>(r);
    }
#elif defined(__APPLE__) || defined(__FreeBSD__) || defined(__OpenBSD__)
    for (std::size_t done = 0; done < len; done += 256) {
        std::size_t n = std::min<std::size_t>(256, len - done);
        if (getentropy(out + done, n) != 0) throw std::runtime_error("spp::SystemRandom: getentropy failed");
    }
#else
    std::FILE* f = std::fopen("/dev/urandom", "rb");
    if (!f || std::fread(out, 1, len, f) != len) {
        if (f) std::fclose(f);
        throw std::runtime_error("spp::SystemRandom: cannot read /dev/urandom");
    }
    std::fclose(f);
#endif
}

static std::uint64_t splitmix64(std::uint64_t& x) {
    std::uint64_t z = (x += 0x9E3779B97F4A7C15ULL);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}
static inline std::uint64_t rotl(std::uint64_t x, int k) { return (x << k) | (x >> (64 - k)); }

SeededRandom::SeededRandom(std::uint64_t seed) {
    for (auto& s : s_) s = splitmix64(seed);
}
std::uint64_t SeededRandom::next() {   // xoshiro256**
    const std::uint64_t result = rotl(s_[1] * 5, 7) * 9;
    const std::uint64_t t = s_[1] << 17;
    s_[2] ^= s_[0]; s_[3] ^= s_[1]; s_[1] ^= s_[2]; s_[0] ^= s_[3];
    s_[2] ^= t; s_[3] = rotl(s_[3], 45);
    return result;
}
void SeededRandom::fill(std::uint8_t* out, std::size_t len) {
    while (len) {
        std::uint64_t v = next();
        for (int i = 0; i < 8 && len; ++i, --len) { *out++ = static_cast<std::uint8_t>(v); v >>= 8; }
    }
}

// ============================================================================ helpers
namespace {

// per-thread scratch integers (avoids allocations in the hot path)
struct Scratch {
    mpz_t g, a, d, x, nm1, span;
    Scratch() { mpz_init(g); mpz_init(a); mpz_init(d); mpz_init(x); mpz_init(nm1); mpz_init(span); }
    ~Scratch() { mpz_clear(g); mpz_clear(a); mpz_clear(d); mpz_clear(x); mpz_clear(nm1); mpz_clear(span); }
};
Scratch& scratch() { thread_local Scratch s; return s; }

std::uint64_t low64(mpz_srcptr n) {
    std::uint64_t lo = mpz_get_ui(n);
    if (sizeof(unsigned long) < 8) {          // Windows: unsigned long is 32 bits
        mpz_t t; mpz_init(t);
        mpz_tdiv_q_2exp(t, n, 32);
        lo = (lo & 0xFFFFFFFFULL) | (static_cast<std::uint64_t>(mpz_get_ui(t) & 0xFFFFFFFFUL) << 32);
        mpz_clear(t);
    }
    return lo;
}

void random_below(mpz_ptr out, mpz_srcptr bound, RandomSource& rng) {  // uniform-ish in [0, bound)
    std::size_t bytes = (mpz_sizeinbase(bound, 2) + 7) / 8 + 8;     // 64 extra bits: negligible bias
    std::vector<std::uint8_t> buf(bytes);
    rng.fill(buf.data(), bytes);
    mpz_import(out, bytes, 1, 1, 0, 0, buf.data());
    mpz_mod(out, out, bound);
}

bool is_even(mpz_srcptr n) { return mpz_tstbit(n, 0) == 0; }

}  // namespace

// ============================================================================ primes
std::vector<unsigned long> odd_primes_upto(unsigned long bound) {
    std::vector<unsigned long> out;
    if (bound < 3) return out;
    std::vector<bool> comp(bound + 1, false);
    for (unsigned long i = 3; i * i <= bound; i += 2)
        if (!comp[i])
            for (unsigned long j = i * i; j <= bound; j += 2 * i) comp[j] = true;
    for (unsigned long i = 3; i <= bound; i += 2)
        if (!comp[i]) out.push_back(i);
    return out;
}

unsigned long recommended_bound(unsigned bits) {
    if (bits <= 64) return 300;
    if (bits <= 256) return 1000;
    if (bits <= 600) return 10000;
    if (bits <= 1200) return 100000;
    if (bits <= 2500) return 300000;
    return 1000000;
}

// ============================================================================ SegmentedPrefilter
struct SegmentedPrefilter::Impl {
    struct Segment {
        std::size_t first = 0, count = 0;   // indices into primes
        bool word = false;                  // product fits into unsigned long -> mpz_gcd_ui fast path
        unsigned long prod_ul = 0;
        mpz_t prod;
        std::size_t bits = 0;
    };
    unsigned long bound = 0;
    std::vector<unsigned long> primes;
    std::vector<Segment> segs;
    double build_s = 0;

    ~Impl() { for (auto& s : segs) mpz_clear(s.prod); }

    unsigned long factor_in(const Segment& s, mpz_srcptr g) const {
        for (std::size_t i = s.first; i < s.first + s.count; ++i)
            if (mpz_divisible_ui_p(g, primes[i])) return primes[i];
        return 0;
    }
};

SegmentedPrefilter::SegmentedPrefilter(unsigned long bound, unsigned s0, double ratio) : p_(new Impl) {
    if (s0 < 1) throw std::invalid_argument("spp::SegmentedPrefilter: s0 must be >= 1");
    if (ratio < 1.0) throw std::invalid_argument("spp::SegmentedPrefilter: ratio must be >= 1");
    auto t0 = std::chrono::steady_clock::now();
    p_->bound = bound;
    p_->primes = odd_primes_upto(bound);
    const auto& P = p_->primes;
    std::size_t i = 0;
    for (int k = 0; i < P.size(); ++k) {
        std::size_t size = std::max<std::size_t>(1, static_cast<std::size_t>(std::llround(s0 * std::pow(ratio, k))));
        size = std::min(size, P.size() - i);
        p_->segs.emplace_back();
        Impl::Segment& s = p_->segs.back();
        s.first = i; s.count = size;
        mpz_init_set_ui(s.prod, 1);
        for (std::size_t j = i; j < i + size; ++j) mpz_mul_ui(s.prod, s.prod, P[j]);
        s.bits = mpz_sizeinbase(s.prod, 2);
        s.word = mpz_fits_ulong_p(s.prod) != 0;
        if (s.word) s.prod_ul = mpz_get_ui(s.prod);
        i += size;
    }
    p_->build_s = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
}
SegmentedPrefilter::~SegmentedPrefilter() = default;
SegmentedPrefilter::SegmentedPrefilter(SegmentedPrefilter&&) noexcept = default;
SegmentedPrefilter& SegmentedPrefilter::operator=(SegmentedPrefilter&&) noexcept = default;

bool SegmentedPrefilter::passes(mpz_srcptr n) const {
    Scratch& sc = scratch();
    for (const auto& s : p_->segs) {
        if (s.word) {
            if (mpz_gcd_ui(nullptr, n, s.prod_ul) != 1) return false;
        } else {
            mpz_gcd(sc.g, n, s.prod);
            if (mpz_cmp_ui(sc.g, 1) != 0) return false;
        }
    }
    return true;
}

FilterResult SegmentedPrefilter::check(mpz_srcptr n) const {
    Scratch& sc = scratch();
    FilterResult r;
    for (std::size_t k = 0; k < p_->segs.size(); ++k) {
        const auto& s = p_->segs[k];
        bool hit;
        if (s.word) {
            unsigned long g = mpz_gcd_ui(nullptr, n, s.prod_ul);
            hit = g != 1;
            if (hit) mpz_set_ui(sc.g, g);
        } else {
            mpz_gcd(sc.g, n, s.prod);
            hit = mpz_cmp_ui(sc.g, 1) != 0;
        }
        if (hit) {
            r.passed = false;
            r.segment = static_cast<int>(k) + 1;
            r.factor = p_->factor_in(s, sc.g);
            return r;
        }
    }
    return r;
}

unsigned long SegmentedPrefilter::bound() const { return p_->bound; }
unsigned long SegmentedPrefilter::largest_prime() const { return p_->primes.empty() ? 2 : p_->primes.back(); }
std::size_t SegmentedPrefilter::prime_count() const { return p_->primes.size(); }
std::size_t SegmentedPrefilter::segment_count() const { return p_->segs.size(); }
std::vector<std::size_t> SegmentedPrefilter::segment_sizes() const {
    std::vector<std::size_t> v; for (const auto& s : p_->segs) v.push_back(s.count); return v;
}
std::vector<std::size_t> SegmentedPrefilter::segment_bits() const {
    std::vector<std::size_t> v; for (const auto& s : p_->segs) v.push_back(s.bits); return v;
}
const std::vector<unsigned long>& SegmentedPrefilter::primes() const { return p_->primes; }
double SegmentedPrefilter::build_seconds() const { return p_->build_s; }

// ============================================================================ TrialDivision
TrialDivision::TrialDivision(unsigned long bound) : primes_(odd_primes_upto(bound)) {}
bool TrialDivision::passes(mpz_srcptr n) const {
    for (unsigned long p : primes_)
        if (mpz_divisible_ui_p(n, p)) return false;
    return true;
}
FilterResult TrialDivision::check(mpz_srcptr n) const {
    FilterResult r;
    for (std::size_t i = 0; i < primes_.size(); ++i)
        if (mpz_divisible_ui_p(n, primes_[i])) { r.passed = false; r.segment = static_cast<int>(i) + 1; r.factor = primes_[i]; return r; }
    return r;
}
unsigned long TrialDivision::largest_prime() const { return primes_.empty() ? 2 : primes_.back(); }
std::size_t TrialDivision::prime_count() const { return primes_.size(); }

// ============================================================================ Miller–Rabin
MRResult miller_rabin(mpz_srcptr n, int rounds, BaseMode mode, RandomSource* rng) {
    MRResult res;
    if (mpz_cmp_ui(n, 2) < 0) return res;
    if (mpz_cmp_ui(n, 3) <= 0) { res.probable_prime = true; return res; }
    if (is_even(n)) return res;

    Scratch& sc = scratch();
    mpz_sub_ui(sc.nm1, n, 1);
    const mp_bitcnt_t s = mpz_scan1(sc.nm1, 0);
    mpz_tdiv_q_2exp(sc.d, sc.nm1, s);
    mpz_sub_ui(sc.span, n, 3);                 // bases a in [2, n-2]

    SeededRandom derived(low64(n) ^ 0x5350505F4D522D31ULL);
    SystemRandom sys;
    RandomSource* src = (mode == BaseMode::FromCandidate) ? static_cast<RandomSource*>(&derived)
                                                          : (rng ? rng : static_cast<RandomSource*>(&sys));
    for (int r = 1; r <= rounds; ++r) {
        res.rounds = r;
        random_below(sc.a, sc.span, *src);
        mpz_add_ui(sc.a, sc.a, 2);
        mpz_powm(sc.x, sc.a, sc.d, n);
        if (mpz_cmp_ui(sc.x, 1) == 0 || mpz_cmp(sc.x, sc.nm1) == 0) continue;
        bool witness = true;
        for (mp_bitcnt_t j = 1; j < s; ++j) {
            mpz_mul(sc.x, sc.x, sc.x);
            mpz_mod(sc.x, sc.x, n);
            if (mpz_cmp(sc.x, sc.nm1) == 0) { witness = false; break; }
        }
        if (witness) return res;               // composite
    }
    res.probable_prime = true;
    return res;
}

// ============================================================================ full test
TestReport is_probable_prime(mpz_srcptr n, const SegmentedPrefilter& pf, int rounds, BaseMode mode, RandomSource* rng) {
    TestReport rep;
    if (mpz_cmp_ui(n, 2) < 0) return rep;
    if (mpz_cmp_ui(n, 2) == 0) { rep.probable_prime = true; return rep; }
    if (is_even(n)) { rep.factor = 2; return rep; }
    if (mpz_cmp_ui(n, pf.largest_prime()) <= 0) {          // Lemma 1 requires n > p_k
        rep.stage = Stage::SmallSieve;
        const auto& P = pf.primes();
        unsigned long v = mpz_get_ui(n);
        rep.probable_prime = std::binary_search(P.begin(), P.end(), v);
        return rep;
    }
    FilterResult f = pf.check(n);
    if (!f.passed) {
        rep.stage = Stage::Prefilter; rep.segment = f.segment; rep.factor = f.factor;
        return rep;
    }
    rep.stage = Stage::MillerRabin;
    MRResult m = miller_rabin(n, rounds, mode, rng);
    rep.probable_prime = m.probable_prime;
    rep.mr_rounds = m.rounds;
    return rep;
}

TestReport is_probable_prime(mpz_srcptr n, int rounds, BaseMode mode, RandomSource* rng) {
    static std::mutex mu;
    static std::map<unsigned long, std::unique_ptr<SegmentedPrefilter>> cache;
    unsigned bits = mpz_cmp_ui(n, 0) == 0 ? 1 : static_cast<unsigned>(mpz_sizeinbase(n, 2));
    unsigned long B = recommended_bound(bits);
    const SegmentedPrefilter* pf;
    {
        std::lock_guard<std::mutex> lk(mu);
        auto& slot = cache[B];
        if (!slot) slot.reset(new SegmentedPrefilter(B));
        pf = slot.get();
    }
    return is_probable_prime(n, *pf, rounds, mode, rng);
}

TestReport is_probable_prime_plain(mpz_srcptr n, int rounds, BaseMode mode, RandomSource* rng) {
    TestReport rep;
    rep.stage = Stage::MillerRabin;
    MRResult m = miller_rabin(n, rounds, mode, rng);
    rep.probable_prime = m.probable_prime;
    rep.mr_rounds = m.rounds;
    if (mpz_cmp_ui(n, 4) < 0 || is_even(n)) rep.stage = Stage::Trivial;
    return rep;
}

// ============================================================================ generation
void GenStats::add(const GenStats& o) {
    candidates += o.candidates; prefilter_rejected += o.prefilter_rejected;
    mr_calls += o.mr_calls; mr_rounds += o.mr_rounds;
}

void random_odd(mpz_ptr out, unsigned bits, RandomSource& rng) {
    if (bits < 2) throw std::invalid_argument("spp::random_odd: bits must be >= 2");
    std::size_t bytes = (bits + 7) / 8;
    std::vector<std::uint8_t> buf(bytes);
    rng.fill(buf.data(), bytes);
    unsigned top = bits % 8;
    if (top) buf[0] &= static_cast<std::uint8_t>((1u << top) - 1);   // keep exactly `bits` bits
    mpz_import(out, bytes, 1, 1, 0, 0, buf.data());
    mpz_setbit(out, bits - 1);
    mpz_setbit(out, 0);
}

static bool small_or_filtered(mpz_srcptr c, const SegmentedPrefilter& pf, bool& is_small_prime) {
    // returns true if the candidate is decided without Miller–Rabin
    if (mpz_cmp_ui(c, pf.largest_prime()) <= 0) {
        const auto& P = pf.primes();
        is_small_prime = std::binary_search(P.begin(), P.end(), mpz_get_ui(c));
        return true;
    }
    is_small_prime = false;
    return !pf.passes(c);
}

void generate_prime(mpz_ptr out, unsigned bits, const SegmentedPrefilter& pf, RandomSource& rng, int rounds, GenStats* st) {
    GenStats local;
    for (;;) {
        random_odd(out, bits, rng);
        ++local.candidates;
        bool small_prime;
        if (small_or_filtered(out, pf, small_prime)) {
            if (small_prime) break;
            ++local.prefilter_rejected;
            continue;
        }
        ++local.mr_calls;
        MRResult m = miller_rabin(out, rounds, BaseMode::FromCandidate);
        local.mr_rounds += m.rounds;
        if (m.probable_prime) break;
    }
    if (st) st->add(local);
}

void generate_prime_plain(mpz_ptr out, unsigned bits, RandomSource& rng, int rounds, GenStats* st) {
    GenStats local;
    for (;;) {
        random_odd(out, bits, rng);
        ++local.candidates;
        ++local.mr_calls;
        MRResult m = miller_rabin(out, rounds, BaseMode::FromCandidate);
        local.mr_rounds += m.rounds;
        if (m.probable_prime) break;
    }
    if (st) st->add(local);
}

void generate_prime_incremental(mpz_ptr out, unsigned bits, const SegmentedPrefilter& pf, RandomSource& rng,
                                int rounds, GenStats* st) {
    GenStats local;
    random_odd(out, bits, rng);
    for (;;) {
        if (mpz_sizeinbase(out, 2) > bits) random_odd(out, bits, rng);   // ran past 2^bits: restart
        ++local.candidates;
        bool small_prime;
        if (small_or_filtered(out, pf, small_prime)) {
            if (small_prime) break;
            ++local.prefilter_rejected;
        } else {
            ++local.mr_calls;
            MRResult m = miller_rabin(out, rounds, BaseMode::FromCandidate);
            local.mr_rounds += m.rounds;
            if (m.probable_prime) break;
        }
        mpz_add_ui(out, out, 2);
    }
    if (st) st->add(local);
}

// ============================================================================ expression parser
namespace {
class Parser {
public:
    explicit Parser(const std::string& s) : s_(s) {}
    BigInt parse() {
        BigInt v = expr();
        skip();
        if (i_ != s_.size()) fail("unexpected character");
        return v;
    }

private:
    const std::string& s_;
    std::size_t i_ = 0;
    [[noreturn]] void fail(const char* m) { throw std::invalid_argument(std::string("spp::parse_expression: ") + m); }
    void skip() { while (i_ < s_.size() && (s_[i_] == ' ' || s_[i_] == '\t')) ++i_; }
    bool eat(char c) { skip(); if (i_ < s_.size() && s_[i_] == c) { ++i_; return true; } return false; }
    BigInt expr() {
        BigInt v = term();
        for (;;) {
            if (eat('+')) { BigInt t = term(); mpz_add(v.get(), v.get(), t.get()); }
            else if (eat('-')) { BigInt t = term(); mpz_sub(v.get(), v.get(), t.get()); }
            else return v;
        }
    }
    BigInt term() {
        BigInt v = power();
        while (eat('*')) { BigInt t = power(); mpz_mul(v.get(), v.get(), t.get()); }
        return v;
    }
    BigInt power() {
        BigInt base = unary();
        skip();
        bool pow_op = eat('^');
        if (!pow_op && i_ + 1 < s_.size() && s_[i_] == '*' && s_[i_ + 1] == '*') { i_ += 2; pow_op = true; }
        if (pow_op) {
            BigInt e = power();
            if (!mpz_fits_ulong_p(e.get()) || mpz_cmp_ui(e.get(), 10000000) > 0) fail("exponent too large");
            mpz_pow_ui(base.get(), base.get(), mpz_get_ui(e.get()));
        }
        return base;
    }
    BigInt unary() {
        if (eat('-')) { BigInt v = unary(); mpz_neg(v.get(), v.get()); return v; }
        if (eat('(')) { BigInt v = expr(); if (!eat(')')) fail("missing ')'"); return v; }
        skip();
        std::size_t start = i_;
        if (i_ + 1 < s_.size() && s_[i_] == '0' && (s_[i_ + 1] == 'x' || s_[i_ + 1] == 'X')) {
            i_ += 2;
            while (i_ < s_.size() && std::isxdigit(static_cast<unsigned char>(s_[i_]))) ++i_;
        } else {
            while (i_ < s_.size() && (std::isdigit(static_cast<unsigned char>(s_[i_])) || s_[i_] == '_')) ++i_;
        }
        if (start == i_) fail("number expected");
        std::string digits;
        for (std::size_t k = start; k < i_; ++k) if (s_[k] != '_') digits += s_[k];
        return BigInt(digits);
    }
};
}  // namespace

BigInt parse_expression(const std::string& text) { return Parser(text).parse(); }


// ---- 1.1: complete generators --------------------------------------------------
double loss_bound(CandidateSet s) { return s == CandidateSet::S1_Odd ? 1.0 : s == CandidateSet::S2_TopTwo ? 2.004 : 2.045; }
int mr_rounds_for_error(unsigned k, double log2_eps, double loss) {
    for (int t = 3; t <= 128; ++t) {
        double l = 1.5 * std::log2((double)k) + t - 0.5 * std::log2((double)t) + 2.0 * (2.0 - std::sqrt((double)t * k));
        double bound = std::log2(loss) + l - std::log2(1.0 - std::exp2(l));
        if (k >= 21 && t <= (int)(k / 9) && bound <= log2_eps) return t;
    }
    return 64;
}
void gen_prime(mpz_ptr out, unsigned bits, CandidateSet s, double log2_eps, const SegmentedPrefilter& pf, RandomSource& rng,
               unsigned long e, GenStats* st) {
    const int t = mr_rounds_for_error(bits, log2_eps, loss_bound(s));
    GenStats local; BigInt nm1, g, ee(e);
    for (;;) {
        random_odd(out, bits, rng); ++local.candidates;
        if (s != CandidateSet::S1_Odd) mpz_setbit(out, bits - 2);
        if (s == CandidateSet::S3_RsaCoprime) {
            mpz_sub_ui(nm1.get(), out, 1); mpz_gcd(g.get(), nm1.get(), ee.get());
            if (mpz_cmp_ui(g.get(), 1) != 0) continue;
        }
        if (!pf.passes(out)) { ++local.prefilter_rejected; continue; }
        ++local.mr_calls;
        MRResult m1 = miller_rabin(out, 1, BaseMode::Random, &rng); local.mr_rounds += m1.rounds;
        if (!m1.probable_prime) continue;
        MRResult m2 = miller_rabin(out, t - 1, BaseMode::Random, &rng); local.mr_rounds += m2.rounds;
        if (m2.probable_prime) break;
    }
    if (st) st->add(local);
}
namespace {
void small_prime_td(mpz_ptr p, unsigned bits, RandomSource& r) {
    for (;;) {
        random_odd(p, bits, r); unsigned long n = mpz_get_ui(p); bool ok = n > 2;
        for (unsigned long d = 3; d * d <= n && ok; d += 2) if (n % d == 0) ok = false;
        if (ok) return;
    }
}
void rand_range(mpz_ptr x, mpz_srcptr lo, mpz_srcptr hi, RandomSource& r) {
    BigInt w, t; mpz_sub(w.get(), hi, lo); mpz_add_ui(w.get(), w.get(), 1);
    random_odd(t.get(), (unsigned)mpz_sizeinbase(w.get(), 2) + 64, r); mpz_mod(t.get(), t.get(), w.get()); mpz_add(x, t.get(), lo);
}
}  // namespace
void generate_provable_prime(mpz_ptr p, unsigned b, const SegmentedPrefilter& pf, RandomSource& r, ProvableStats* st) {
    ProvableStats dummy; ProvableStats& s = st ? *st : dummy;
    if (b <= 32) { small_prime_td(p, b < 2 ? 2 : b, r); return; }
    ++s.levels;
    BigInt q; generate_provable_prime(q.get(), b / 2 + 1, pf, r, &s);
    BigInt lo, hi, twoq, t, rem, R, x, y, g, a(2UL), e;
    mpz_mul_2exp(twoq.get(), q.get(), 1);
    mpz_set_ui(t.get(), 1); mpz_mul_2exp(t.get(), t.get(), b - 1); mpz_sub_ui(t.get(), t.get(), 1);
    mpz_tdiv_q(lo.get(), t.get(), twoq.get()); mpz_tdiv_r(rem.get(), t.get(), twoq.get()); if (mpz_sgn(rem.get())) mpz_add_ui(lo.get(), lo.get(), 1);
    mpz_set_ui(t.get(), 1); mpz_mul_2exp(t.get(), t.get(), b); mpz_sub_ui(t.get(), t.get(), 2); mpz_tdiv_q(hi.get(), t.get(), twoq.get());
    for (;;) {
        rand_range(R.get(), lo.get(), hi.get(), r); mpz_mul(p, twoq.get(), R.get()); mpz_add_ui(p, p, 1); ++s.candidates;
        if (!pf.passes(p)) { ++s.prefilter_rejected; continue; }
        mpz_mul_2exp(e.get(), R.get(), 1);
        mpz_powm(x.get(), a.get(), e.get(), p); mpz_powm(y.get(), x.get(), q.get(), p); s.modexps += 2;
        if (mpz_cmp_ui(y.get(), 1) != 0) continue;
        mpz_sub_ui(x.get(), x.get(), 1); mpz_gcd(g.get(), x.get(), p);
        if (mpz_cmp_ui(g.get(), 1) == 0) { mpz_mul(x.get(), q.get(), q.get()); if (mpz_cmp(x.get(), p) < 0) ++s.q_below_sqrt; return; }  // Pocklington (Lemma 2): 2q + 1 > sqrt(p) for every b -> p is prime
    }
}
void generate_rsa_primes_fips(mpz_ptr p, mpz_ptr q, unsigned nlen, unsigned long e, const SegmentedPrefilter& pf, RandomSource& r, double log2_eps) {
    const unsigned hb = nlen / 2; BigInt d, lim;
    gen_prime(p, hb, CandidateSet::S3_RsaCoprime, log2_eps - 1.0, pf, r, e);
    mpz_set_ui(lim.get(), 1); mpz_mul_2exp(lim.get(), lim.get(), hb - 100);
    for (;;) {
        gen_prime(q, hb, CandidateSet::S3_RsaCoprime, log2_eps - 1.0, pf, r, e);
        mpz_sub(d.get(), p, q); if (mpz_sgn(d.get()) < 0) mpz_neg(d.get(), d.get());
        if (mpz_cmp(d.get(), lim.get()) > 0) return;
    }
}

}  // namespace spp
