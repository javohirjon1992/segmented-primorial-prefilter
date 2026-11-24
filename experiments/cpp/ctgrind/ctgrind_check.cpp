// ctgrind analysis (A. Langley's method): the limbs of the candidate are marked as undefined for Valgrind's Memcheck, which
// then reports every conditional branch (or cmov) and every memory address that depends on them as an error.
// The decision (accept/reject) is declassified after the call (VALGRIND_MAKE_MEM_DEFINED), because a generator discards
// rejected candidates anyway.
// Ishlatish: valgrind --error-limit=no ./ctgrind_check <variant> <bits> <B> <N>
//   variant: TD_early | SEG_early | SEG_full | CT_TD | CT_TREE
#include "ct_impl.hpp"
#include <valgrind/memcheck.h>
#include <cstdio>
#include <cstdlib>
#include <string>
int main(int argc, char** argv) {
    std::string v = argv[1]; unsigned bits = std::atoi(argv[2]); unsigned long B = std::strtoul(argv[3], nullptr, 10); int N = std::atoi(argv[4]);
    std::size_t nl = (bits + 63) / 64;
    SegmentedPrefilter seg(B); TrialDivision td(B); CtTrialDivision cttd(B, nl); CtTree tree(B, nl);
    std::vector<BigInt> mods; { std::size_t i = 0; const auto& P = seg.primes();
        for (std::size_t sz : seg.segment_sizes()) { BigInt m(1UL); for (std::size_t j = i; j < i + sz; ++j) mpz_mul_ui(m.get(), m.get(), P[j]); mods.push_back(m); i += sz; } }
    SeededRandom r(424242 + bits); BigInt g; int accepted = 0;
    for (int k = 0; k < N; ++k) {
        BigInt n;
        if (k % 3 == 0) generate_prime(n.get(), bits, seg, r, 20);          // primes (accepted candidates)
        else random_odd(n.get(), bits, r);                                    // random odd candidates (mostly rejected)
        mp_limb_t* d = n.get()->_mp_d; std::size_t sz = (std::size_t)n.get()->_mp_size;   // bit uzunligi ochiq, limblar maxfiy
        VALGRIND_MAKE_MEM_UNDEFINED(d, sz * sizeof(mp_limb_t));
        int res;
        if (v == "TD_early") res = td.passes(n.get());
        else if (v == "SEG_early") res = seg.passes(n.get());
        else if (v == "SEG_full") { int acc = 0; for (auto& m : mods) { mpz_gcd(g.get(), n.get(), m.get()); acc |= mpz_cmp_ui(g.get(), 1) != 0; } res = !acc; }
        else if (v == "CT_TD") res = cttd.passes(n.get());
        else res = tree.passes(n.get());
        VALGRIND_MAKE_MEM_DEFINED(&res, sizeof res);                          // qarorni ochish
        VALGRIND_MAKE_MEM_DEFINED(d, sz * sizeof(mp_limb_t));
        accepted += res;
    }
    std::printf("variant=%s bits=%u B=%lu N=%d accepted=%d tree_groups=%zu depth=%zu\n", v.c_str(), bits, B, N, accepted, tree.groups(), tree.depth());
}
