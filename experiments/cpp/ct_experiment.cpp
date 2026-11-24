// Timing-leakage experiment for prefilter variants, including the constant-time prefilter CT-TREE.
//
// Method: dudect (Reparaz, Balasch, Verbauwhede, 2017). Two classes of inputs (e.g. a fixed prime and random
// primes) are measured in random order, and Welch's t statistic is computed. |t| > 4.5 indicates that the
// running time depends on the input, i.e. that information leaks.
//
// Variants:
//   TD_early   — trial division, stops at the first divisor (customary implementation)
//   SEG_early  — segmented prefilter of the paper (stops at the first GCD > 1)
//   SEG_full   — all segments are evaluated (no early exit), but mpz_gcd is data dependent
//   CT_TD      — constant-time trial division: mpn_sec_div_r for every p, no early exit
//   CT_TREE    — constant-time hierarchical primorial prefilter: descends a product tree of the small primes
//                with mpn_sec_div_r and finishes with branch-free word-level Barrett reductions.
//                All operand sizes depend only on B and on the bit length of n.
//
// Tests:
//   T1  fixed prime versus random primes (does the value of the final prime leak?)
//   T2  numbers divisible by 3 versus numbers divisible only by the largest sieving prime
//       (leakage through rejected candidates, e.g. in incremental search)
//   T3  modular exponentiation inside Miller–Rabin: mpz_powm versus mpz_powm_sec
//
// Usage:  ct_experiment [bits=1024] [B=1000] [measurements_per_class=200000] [output.csv]
#include "spp/spp.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <functional>
#include <map>
#include <memory>
#include <string>
#include <vector>

#if defined(__x86_64__) || defined(_M_X64) || defined(__i386__)
#  if defined(_MSC_VER)
#    include <intrin.h>
#  else
#    include <x86intrin.h>
#  endif
static inline std::uint64_t ticks() { _mm_lfence(); std::uint64_t t = __rdtsc(); _mm_lfence(); return t; }
static const char* TICK_UNIT = "CPU cycles (rdtsc)";
#else
static inline std::uint64_t ticks() {
    return std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now().time_since_epoch()).count();
}
static const char* TICK_UNIT = "nanoseconds";
#endif

using namespace spp;

// ----------------------------------------------------------------------------- constant-time helpers
// The ctgrind analysis showed that Clang 18.1 turns the masked subtraction of the Barrett reduction into a conditional branch.
// The masks are therefore hidden from the optimizer by an empty inline-assembly barrier (disable with -DCT_NO_BARRIER; the
// timing measurements of the paper were made with GCC builds without the barrier).
#if defined(CT_NO_BARRIER) || defined(_MSC_VER)
#  define CT_BARRIER(x) ((void)0)
#else
#  define CT_BARRIER(x) __asm__ volatile("" : "+r"(x))
#endif
static inline std::uint64_t ct_is_zero(std::uint64_t x) { return ((x | (0 - x)) >> 63) ^ 1; }

struct Barrett {                      // r < 2^32, p < 2^32: r mod p by multiplication (no division)
    std::uint64_t p, m;
    explicit Barrett(std::uint64_t p_) : p(p_) {
#if defined(__SIZEOF_INT128__)
        m = static_cast<std::uint64_t>((static_cast<unsigned __int128>(1) << 64) / p);
#else
        m = (~0ULL) / p;   // approximate for MSVC (still exact for r < 2^32 after the final correction)
#endif
    }
    inline std::uint64_t mod(std::uint64_t r) const {
#if defined(__SIZEOF_INT128__)
        std::uint64_t q = static_cast<std::uint64_t>((static_cast<unsigned __int128>(r) * m) >> 64);
#else
        std::uint64_t q = __umulh(r, m);
#endif
        std::uint64_t rem = r - q * p;
        std::uint64_t m1 = 0 - static_cast<std::uint64_t>(rem >= p); CT_BARRIER(m1); rem -= p & m1;
        std::uint64_t m2 = 0 - static_cast<std::uint64_t>(rem >= p); CT_BARRIER(m2); rem -= p & m2;
        return rem;
    }
};

// CT_TD: constant-time residue for every small prime
class CtTrialDivision {
public:
    CtTrialDivision(unsigned long B, std::size_t nlimbs) : primes_(odd_primes_upto(B)), nn_(nlimbs) {
        tmp_.resize(nn_); scratch_.resize(mpn_sec_div_r_itch(nn_, 1) + 4);
    }
    bool passes(mpz_srcptr n) {
        const mp_limb_t* src = mpz_limbs_read(n);
        std::uint64_t acc = 0;
        for (unsigned long p : primes_) {
            std::memcpy(tmp_.data(), src, nn_ * sizeof(mp_limb_t));
            mp_limb_t d = p;
            mpn_sec_div_r(tmp_.data(), nn_, &d, 1, scratch_.data());
            acc |= ct_is_zero(tmp_[0]);
        }
        return acc == 0;
    }
private:
    std::vector<unsigned long> primes_;
    std::size_t nn_;
    std::vector<mp_limb_t> tmp_, scratch_;
};

// CT_TREE: constant-time hierarchical primorial prefilter
class CtTree {
public:
    CtTree(unsigned long B, std::size_t nlimbs) : nn_(nlimbs) {
        std::vector<unsigned long> P = odd_primes_upto(B);
        // leaves: groups of consecutive primes whose product is below 2^32
        std::vector<std::vector<unsigned long>> groups;
        std::uint64_t prod = 1; std::vector<unsigned long> cur;
        for (unsigned long p : P) {
            if (!cur.empty() && prod * p >= (1ULL << 32)) { groups.push_back(cur); cur.clear(); prod = 1; }
            cur.push_back(p); prod *= p;
        }
        if (!cur.empty()) groups.push_back(cur);
        for (auto& g : groups) {
            std::vector<Barrett> bs; for (auto p : g) bs.emplace_back(p);
            leaves_.push_back(bs);
        }
        // tree levels (0 = leaves)
        std::vector<Node> lvl;
        for (auto& g : groups) {
            Node nd; mpz_t v; mpz_init_set_ui(v, 1);
            for (auto p : g) mpz_mul_ui(v, v, p);
            nd.limbs.assign(mpz_limbs_read(v), mpz_limbs_read(v) + mpz_size(v)); mpz_clear(v);
            lvl.push_back(nd);
        }
        levels_.push_back(lvl);
        while (levels_.back().size() > 1) {
            const auto& L = levels_.back(); std::vector<Node> up;
            for (std::size_t i = 0; i < L.size(); i += 2) {
                if (i + 1 < L.size()) up.push_back(mul(L[i], L[i + 1])); else up.push_back(L[i]);
            }
            levels_.push_back(up);
        }
        // starting level: the highest level whose nodes are all not longer than n
        top_ = 0;
        for (std::size_t l = 0; l < levels_.size(); ++l) {
            std::size_t mx = 0; for (auto& nd : levels_[l]) mx = std::max(mx, nd.limbs.size());
            if (mx <= nn_) top_ = l;
        }
        // storage for the residues
        res_.resize(levels_.size());
        for (std::size_t l = 0; l <= top_; ++l) {
            res_[l].resize(levels_[l].size());
            for (std::size_t i = 0; i < levels_[l].size(); ++i) res_[l][i].resize(levels_[l][i].limbs.size());
        }
        std::size_t maxn = nn_; tmp_.resize(maxn + 4);
        scratch_.resize(mpn_sec_div_r_itch(maxn, maxn) + 16);
        groups_ = groups.size();
    }
    bool passes(mpz_srcptr n) {
        const mp_limb_t* src = mpz_limbs_read(n);
        // every node of the starting level: n mod node
        for (std::size_t i = 0; i < levels_[top_].size(); ++i) reduce(src, nn_, levels_[top_][i], res_[top_][i]);
        // descend the tree
        for (std::size_t l = top_; l-- > 0;) {
            for (std::size_t i = 0; i < levels_[l].size(); ++i) {
                const auto& par = res_[l + 1][i / 2];
                reduce(par.data(), par.size(), levels_[l][i], res_[l][i]);
            }
        }
        // leaves: word-level Barrett reduction
        std::uint64_t acc = 0;
        for (std::size_t g = 0; g < leaves_.size(); ++g) {
            std::uint64_t r = res_[0][g][0];
            for (const auto& b : leaves_[g]) acc |= ct_is_zero(b.mod(r));
        }
        return acc == 0;
    }
    std::size_t groups() const { return groups_; }
    std::size_t depth() const { return top_ + 1; }
private:
    struct Node { std::vector<mp_limb_t> limbs; };
    std::size_t nn_, top_ = 0, groups_ = 0;
    std::vector<std::vector<Barrett>> leaves_;
    std::vector<std::vector<Node>> levels_;
    std::vector<std::vector<std::vector<mp_limb_t>>> res_;
    std::vector<mp_limb_t> tmp_, scratch_;
    static Node mul(const Node& a, const Node& b) {
        mpz_t x, y, z; mpz_init(x); mpz_init(y); mpz_init(z);
        mpz_import(x, a.limbs.size(), -1, sizeof(mp_limb_t), 0, 0, a.limbs.data());
        mpz_import(y, b.limbs.size(), -1, sizeof(mp_limb_t), 0, 0, b.limbs.data());
        mpz_mul(z, x, y);
        Node nd; nd.limbs.assign(mpz_limbs_read(z), mpz_limbs_read(z) + mpz_size(z));
        mpz_clear(x); mpz_clear(y); mpz_clear(z);
        return nd;
    }
    void reduce(const mp_limb_t* a, std::size_t an, const Node& d, std::vector<mp_limb_t>& out) {
        std::memcpy(tmp_.data(), a, an * sizeof(mp_limb_t));
        mpn_sec_div_r(tmp_.data(), an, d.limbs.data(), d.limbs.size(), scratch_.data());
        std::memcpy(out.data(), tmp_.data(), d.limbs.size() * sizeof(mp_limb_t));
    }
};

// ----------------------------------------------------------------------------- statistika
struct Welch {
    double t, mean0, mean1;
};
static Welch welch(const std::vector<double>& a, const std::vector<double>& b) {
    auto mv = [](const std::vector<double>& v, double& m, double& s2) {
        m = 0; for (double x : v) m += x; m /= v.size();
        s2 = 0; for (double x : v) s2 += (x - m) * (x - m); s2 /= (v.size() - 1);
    };
    double m0, v0, m1, v1; mv(a, m0, v0); mv(b, m1, v1);
    double den = std::sqrt(v0 / a.size() + v1 / b.size());
    return {den > 0 ? (m0 - m1) / den : 0.0, m0, m1};
}
// as in dudect: largest |t| over several percentile crops
static double max_abs_t(const std::vector<double>& a, const std::vector<double>& b, double& m0, double& m1) {
    std::vector<double> all(a); all.insert(all.end(), b.begin(), b.end());
    std::sort(all.begin(), all.end());
    double best = 0; m0 = m1 = 0;
    for (double pct : {0.50, 0.75, 0.90, 0.95, 0.99, 1.0}) {
        double thr = all[std::min(all.size() - 1, static_cast<std::size_t>(pct * (all.size() - 1)))];
        std::vector<double> x, y;
        for (double v : a) if (v <= thr) x.push_back(v);
        for (double v : b) if (v <= thr) y.push_back(v);
        if (x.size() < 100 || y.size() < 100) continue;
        Welch w = welch(x, y);
        if (pct == 1.0) { m0 = w.mean0; m1 = w.mean1; }
        best = std::max(best, std::fabs(w.t));
    }
    return best;
}

// measure the two classes in random order
static void measure(const std::function<void(mpz_srcptr)>& f, const std::vector<BigInt>& c0, const std::vector<BigInt>& c1,
                    std::size_t per_class, std::vector<double>& t0, std::vector<double>& t1, std::uint64_t seed) {
    SeededRandom rng(seed);
    t0.clear(); t1.clear(); t0.reserve(per_class); t1.reserve(per_class);
    // Every input is copied into one fixed buffer (outside the timed region), so that both classes are read
    // from the same memory address and cache effects cannot cause spurious "leakage".
    BigInt work;
    mpz_set(work.get(), c1[0].get());
    for (int w = 0; w < 1000; ++w) { mpz_set(work.get(), c0[w % c0.size()].get()); f(work.get()); }   // warm-up
    while (t0.size() < per_class || t1.size() < per_class) {
        bool cls = (rng.next_u64() & 1) != 0;
        if (!cls && t0.size() >= per_class) cls = true;
        if (cls && t1.size() >= per_class) cls = false;
        const auto& pool = cls ? c1 : c0;
        mpz_set(work.get(), pool[rng.next_u64() % pool.size()].get());
        std::uint64_t a = ticks(); f(work.get()); std::uint64_t b = ticks();
        (cls ? t1 : t0).push_back(static_cast<double>(b - a));
    }
}

static double median(std::vector<double> v) { std::sort(v.begin(), v.end()); return v[v.size() / 2]; }

int main(int argc, char** argv) {
    unsigned bits = argc > 1 ? std::atoi(argv[1]) : 1024;
    unsigned long B = argc > 2 ? std::strtoul(argv[2], nullptr, 10) : 1000;
    std::size_t N = argc > 3 ? std::strtoull(argv[3], nullptr, 10) : 200000;
    std::string out = argc > 4 ? argv[4] : ("ct_results_" + std::to_string(bits) + "_" + std::to_string(B) + ".csv");
    std::size_t nlimbs = (bits + GMP_LIMB_BITS - 1) / GMP_LIMB_BITS;

    std::printf("=== Timing-leakage experiment | %u bits, B = %lu, %zu measurements per class | unit: %s\n",
                bits, B, N, TICK_UNIT);
    std::printf("GMP %s, spp %s\n", gmp_version, version());

    SeededRandom rng(20261101);
    SegmentedPrefilter seg(B);
    TrialDivision td(B);
    CtTrialDivision cttd(B, nlimbs);
    CtTree cttree(B, nlimbs);
    std::printf("CT_TREE: %zu leaf groups, tree depth %zu\n", cttree.groups(), cttree.depth());

    // --- input data
    std::vector<BigInt> primes(512), fixed(1), div3(512), divmax(512), random_odd_v(2000);
    for (auto& p : primes) generate_prime(p.get(), bits, seg, rng);
    fixed[0] = primes[0];
    unsigned long pmax = seg.largest_prime();
    for (auto& x : div3) {                                  // odd number of exactly `bits` bits, divisible by 3
        do { random_odd(x.get(), bits, rng); mpz_sub_ui(x.get(), x.get(), mpz_fdiv_ui(x.get(), 3)); if (mpz_tstbit(x.get(), 0) == 0) mpz_add_ui(x.get(), x.get(), 3); }
        while (mpz_sizeinbase(x.get(), 2) != bits);
    }
    for (auto& x : divmax) {                                // divisible by pmax only (no other small divisor)
        BigInt q;
        do { generate_prime(q.get(), bits - static_cast<unsigned>(std::log2(pmax)) - 1, seg, rng); mpz_mul_ui(x.get(), q.get(), pmax); }
        while (mpz_sizeinbase(x.get(), 2) != bits);
    }
    for (auto& x : random_odd_v) random_odd(x.get(), bits, rng);

    // --- correctness: all variants make the same decision
    std::size_t mism = 0;
    for (auto& x : random_odd_v) {
        bool a = td.passes(x.get()), b = seg.passes(x.get()), c = cttd.passes(x.get()), d = cttree.passes(x.get());
        mism += (a != b) + (a != c) + (a != d);
    }
    for (auto& x : primes) mism += !cttree.passes(x.get()) + !cttd.passes(x.get());
    for (auto& x : divmax) mism += cttree.passes(x.get()) + cttd.passes(x.get()) + seg.passes(x.get());
    std::printf("Correctness: number of disagreements between variants = %zu (must be 0)\n\n", mism);

    struct Variant { const char* name; std::function<void(mpz_srcptr)> f; };
    volatile int sink = 0;
    std::vector<Variant> V = {
        {"TD_early", [&](mpz_srcptr n) { sink += td.passes(n); }},
        {"SEG_early", [&](mpz_srcptr n) { sink += seg.passes(n); }},
        {"SEG_full", [&](mpz_srcptr n) { (void)n; }},       // replaced below
        {"CT_TD", [&](mpz_srcptr n) { sink += cttd.passes(n); }},
        {"CT_TREE", [&](mpz_srcptr n) { sink += cttree.passes(n); }},
    };
    // SEG_full: all segments, no early exit
    {
        std::vector<BigInt> mods;
        std::size_t i = 0; const auto& P = seg.primes();
        for (std::size_t sz : seg.segment_sizes()) {
            BigInt m(1UL); for (std::size_t j = i; j < i + sz; ++j) mpz_mul_ui(m.get(), m.get(), P[j]);
            mods.push_back(m); i += sz;
        }
        auto g = std::make_shared<BigInt>();
        V[2].f = [mods, g, &sink](mpz_srcptr n) { int acc = 0; for (auto& m : mods) { mpz_gcd(g->get(), n, m.get()); acc |= mpz_cmp_ui(g->get(), 1) != 0; } sink += acc; };
    }

    std::FILE* csv = std::fopen(out.c_str(), "w");
    std::fprintf(csv, "bits,B,variant,test,median_class0,median_class1,mean_class0,mean_class1,max_abs_t,leak\n");
    std::printf("%-10s %-34s %12s %12s %9s %8s\n", "variant", "test", "median A", "median B", "max|t|", "leak");
    double base_cost = 0;
    std::vector<double> a, b;
    for (auto& v : V) {
        struct T { const char* name; const std::vector<BigInt>* c0; const std::vector<BigInt>* c1; };
        T tests[] = {{"T1 fixed prime vs random primes", &fixed, &primes}, {"T2 divisible by 3 vs by pmax", &div3, &divmax}};
        for (auto& t : tests) {
            measure(v.f, *t.c0, *t.c1, N, a, b, 777);
            double m0, m1, mt = max_abs_t(a, b, m0, m1);
            double md0 = median(a), md1 = median(b);
            if (!std::strcmp(v.name, "TD_early") && t.c0 == &fixed) base_cost = md1;
            bool leak = mt > 4.5;
            std::printf("%-10s %-34s %12.0f %12.0f %9.1f %8s\n", v.name, t.name, md0, md1, mt, leak ? "YES" : "no");
            std::fprintf(csv, "%u,%lu,%s,%s,%.1f,%.1f,%.1f,%.1f,%.2f,%d\n", bits, B, v.name, t.name, md0, md1, m0, m1, mt, leak ? 1 : 0);
        }
    }
    // cost: median on random primes, relative to TD_early
    std::printf("\nCost (median on random primes, TD_early = 1):\n");
    for (auto& v : V) {
        measure(v.f, primes, primes, N / 4, a, b, 99);
        double md = median(b);
        std::printf("  %-10s %10.0f  (%.2fx)\n", v.name, md, md / base_cost);
        std::fprintf(csv, "%u,%lu,%s,cost (random primes),%.1f,%.1f,,,,\n", bits, B, v.name, md, md / base_cost);
    }

    // T3: Miller–Rabin exponentiation
    std::printf("\nT3: exponentiation inside Miller-Rabin (base 2, exponent (n-1)/2^s)\n");
    std::map<mpz_srcptr, std::shared_ptr<BigInt>> dexp;     // n -> (n-1)/2^s
    for (const auto* pool : {&fixed, &primes})
        for (const auto& p : *pool) {
            auto e = std::make_shared<BigInt>();
            mpz_sub_ui(e->get(), p.get(), 1); mpz_tdiv_q_2exp(e->get(), e->get(), mpz_scan1(e->get(), 0));
            dexp[p.get()] = e;
        }
    // exponent matching the n held in the work buffer: the lowest limb of n is used as key (before timing)
    std::map<mp_limb_t, mpz_srcptr> by_low;
    for (auto& kv : dexp) by_low[mpz_limbs_read(kv.first)[0]] = kv.second->get();
    auto lookup = [&](mpz_srcptr n) { return by_low[mpz_limbs_read(n)[0]]; };
    auto r = std::make_shared<BigInt>(); BigInt two(2UL);
    struct PV { const char* name; bool sec; };
    for (PV pv : {PV{"mpz_powm", false}, PV{"mpz_powm_sec", true}}) {
        auto f = [&, pv](mpz_srcptr n) {
            mpz_srcptr e = lookup(n);
            if (pv.sec) mpz_powm_sec(r->get(), two.get(), e, n); else mpz_powm(r->get(), two.get(), e, n);
        };
        std::size_t n3 = std::max<std::size_t>(2000, N / 50);
        measure(f, fixed, primes, n3, a, b, 555);
        double m0, m1, mt = max_abs_t(a, b, m0, m1);
        std::printf("  %-14s median A %10.0f  median B %10.0f  max|t| %6.1f  leak: %s\n", pv.name, median(a), median(b), mt, mt > 4.5 ? "YES" : "no");
        std::fprintf(csv, "%u,%lu,%s,T3 exponentiation,%.1f,%.1f,%.1f,%.1f,%.2f,%d\n", bits, B, pv.name, median(a), median(b), m0, m1, mt, mt > 4.5);
    }
    std::fclose(csv);
    std::printf("\nResults written to: %s\n", out.c_str());
    return mism == 0 ? 0 : 2;
}
