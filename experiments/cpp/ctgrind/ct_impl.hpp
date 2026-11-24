// Constant-time prefilters of ct_experiment.cpp (CT_TD, CT_TREE), extracted for the ctgrind analysis.
// -DCT_ASM_BARRIER: the masks of the Barrett reduction are hidden from the optimizer by an empty inline-assembly barrier.
#pragma once
#include "spp/spp.hpp"
#include <algorithm>
#include <cstdint>
#include <cstring>
#include <vector>
#if defined(CT_ASM_BARRIER)
#  define CT_BARRIER(x) __asm__ volatile("" : "+r"(x))
#else
#  define CT_BARRIER(x) ((void)0)
#endif
using namespace spp;
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

