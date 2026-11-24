// Safe primes: segmented versus batch (remainder-tree) screening of the products q(2q+1).
#include "spp/spp.hpp"
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <algorithm>
using namespace spp;
static double now() { return std::chrono::duration<double, std::nano>(std::chrono::steady_clock::now().time_since_epoch()).count(); }
static void prod_range(mpz_ptr out, const std::vector<unsigned long>& P, size_t a, size_t b) {   // balanced product
    if (b - a <= 16) { mpz_set_ui(out, 1); for (size_t i = a; i < b; ++i) mpz_mul_ui(out, out, P[i]); return; }
    BigInt L, R; size_t m = (a + b) / 2; prod_range(L.get(), P, a, m); prod_range(R.get(), P, m, b); mpz_mul(out, L.get(), R.get());
}
// batch screening: product tree + remainder tree
static size_t batch_pass(const std::vector<BigInt>& X, size_t a, size_t b, mpz_srcptr M, std::vector<char>& ok) {
    std::vector<std::vector<BigInt>> T; T.emplace_back(X.begin() + a, X.begin() + b);
    while (T.back().size() > 1) { auto& L = T.back(); std::vector<BigInt> U((L.size() + 1) / 2);
        for (size_t i = 0; i < U.size(); ++i) { if (2 * i + 1 < L.size()) mpz_mul(U[i].get(), L[2 * i].get(), L[2 * i + 1].get()); else mpz_set(U[i].get(), L[2 * i].get()); }
        T.push_back(std::move(U)); }
    std::vector<BigInt> rem(1); mpz_mod(rem[0].get(), M, T.back()[0].get());
    for (int l = (int)T.size() - 2; l >= 0; --l) { std::vector<BigInt> nr(T[l].size());
        for (size_t i = 0; i < nr.size(); ++i) mpz_mod(nr[i].get(), rem[i / 2].get(), T[l][i].get()); rem.swap(nr); }
    size_t pass = 0; BigInt g;
    for (size_t i = 0; i < rem.size(); ++i) { mpz_gcd(g.get(), rem[i].get(), T[0][i].get()); ok[a + i] = mpz_cmp_ui(g.get(), 1) == 0; pass += ok[a + i]; }
    return pass;
}
int main(int argc, char** argv) {
    unsigned bits = std::atoi(argv[1]); size_t N = std::atoi(argv[2]); unsigned long maxB = argc > 3 ? std::strtoul(argv[3], nullptr, 10) : 30000000UL;
    SeededRandom r(777 + bits); std::vector<BigInt> Q(N), X(N);
    for (size_t i = 0; i < N; ++i) { random_odd(Q[i].get(), bits - 1, r); BigInt p; mpz_mul_2exp(p.get(), Q[i].get(), 1); mpz_add_ui(p.get(), p.get(), 1); mpz_mul(X[i].get(), Q[i].get(), p.get()); }
    std::vector<double> mt; for (int i = 0; i < 25; ++i) { double a = now(); miller_rabin(Q[i].get(), 1); mt.push_back(now() - a); }
    std::sort(mt.begin(), mt.end()); double cmr = mt[12];
    std::printf("bits,B,method,batch,pre_ns,pass,cmr_ns,build_s\n");
    for (unsigned long B : {100000UL, 1000000UL, 3000000UL, 10000000UL, 30000000UL}) {
        if (B > maxB) break;
        auto P = odd_primes_upto(B); double t0 = now(); BigInt M; prod_range(M.get(), P, 0, P.size()); double build = (now() - t0) / 1e9;
        std::vector<char> ok(N);
        for (size_t bs : {256UL, 1024UL, 4096UL}) {
            if (bs > N) continue; double a = now(); size_t pass = 0;
            for (size_t s = 0; s < N; s += bs) pass += batch_pass(X, s, std::min(N, s + bs), M.get(), ok);
            std::printf("%u,%lu,batch,%zu,%.1f,%.6f,%.1f,%.2f\n", bits, B, bs, (now() - a) / N, (double)pass / N, cmr, build); std::fflush(stdout);
        }
        if (B <= 10000000UL) {
            SegmentedPrefilter pf(B); size_t n = std::min<size_t>(N, 1000), pass = 0; double a = now();
            for (size_t i = 0; i < n; ++i) pass += pf.passes(X[i].get());
            std::printf("%u,%lu,segmented,1,%.1f,%.6f,%.1f,%.2f\n", bits, B, (now() - a) / n, (double)pass / n, cmr, pf.build_seconds()); std::fflush(stdout);
        }
    }
}
