# spp — segmented primorial prefilter for Miller–Rabin (C++17 library)

`spp` is the GMP-based C++ implementation of the method of the paper
"Segmented Primorial Prefiltering and Cost-Optimal Sieving Bounds for Prime Generation".

1. The odd primes from 3 to a bound *B* are split into segments whose sizes grow geometrically
   (*s*<sub>0</sub>, *s*<sub>0</sub>*r*, *s*<sub>0</sub>*r*<sup>2</sup>, ...; defaults *s*<sub>0</sub> = 5, *r* = 2).
2. The product *M*<sub>*i*</sub> of every segment is computed once.
3. A candidate *n* is rejected as soon as gcd(*n*, *M*<sub>*i*</sub>) > 1 for some segment.
4. Candidates that pass the prefilter go to the Miller–Rabin test.

For *n* larger than the largest sieving prime, the prefilter rejects exactly the candidates that trial division over
the same primes rejects (Lemma 1 of the paper), so it changes the speed of a generator but not its output.

## Features

- `SegmentedPrefilter` — the prefilter of the paper. Segments whose product fits into a machine word use the
  single-word GCD (`mpz_gcd_ui`), following the word-size rule for the first segment.
- `TrialDivision` — classical trial division over the same primes (baseline).
- `miller_rabin()` — Miller–Rabin test with random bases (`BaseMode::Random`) or bases derived deterministically
  from the candidate (`BaseMode::FromCandidate`, for reproducible comparisons).
- `is_probable_prime()` — prefilter followed by Miller–Rabin; *B* is chosen from the bit length.
- `gen_prime(out, bits, S, log2_eps, pf, rng, e)` — GenPrime (Algorithm 2 of the paper): independent candidates
  from an admissible set `CandidateSet::S1_Odd`, `S2_TopTwo` (two top bits set), or `S3_RsaCoprime`
  (S2 with gcd(*e*, *n* − 1) = 1), the prefilter, and Miller–Rabin rounds with bases drawn from `rng`,
  independently of the candidate. The number of rounds is the smallest *t* ≥ 3 with
  *L*(*S*)·*p*<sub>*b,t*</sub>/(1 − *p*<sub>*b,t*</sub>) ≤ 2<sup>log2_eps</sup> (average-case error bound).
- `generate_provable_prime()` — provable primes by a Maurer-type recursion with Pocklington's criterion
  (Lemma 2 of the paper: *p* = 2*Rq* + 1 with 2*q* + 1 > √*p*), with the prefilter applied at every level;
  `ProvableStats::q_below_sqrt` counts the levels at which *q*² < *p*.
- `generate_rsa_primes_fips()` — pairs of RSA primes with the conditions of FIPS 186-5 (two top bits,
  gcd(*e*, *p* − 1) = 1, |*p* − *q*| > 2<sup>*n*/2−100</sup>); each prime comes from `gen_prime(S3, log2_eps − 1)`,
  so the error probability of the pair is at most 2<sup>log2_eps</sup>.
- `generate_prime()`, `generate_prime_plain()`, `generate_prime_incremental()` — benchmarking generators. With the
  same seed, `generate_prime()` and `generate_prime_plain()` test exactly the same candidates and return the same
  prime; they use candidate-derived bases and are not intended for key generation.
- `recommended_bound(bits)` — default bound: 10<sup>4</sup> up to 600 bits, 10<sup>5</sup> up to 1200 bits,
  3·10<sup>5</sup> up to 2500 bits, 10<sup>6</sup> above (near the best bounds measured in the paper; the cost
  curves are flat, but the bound should be calibrated on the target platform with `spp_bench --sweep`).
- `mr_rounds_for_error(bits, log2_eps, loss)` and `loss_bound(CandidateSet)` — round counts of Table 2 of the paper
  (`loss_bound`: 1, 2.004, 2.045 for S1, S2, S3).
- `parse_expression()` — reads expressions such as `2^1279-1`, `(2^127-1)*(2^89-1)`, `0xff`, `1_000_007`.
- `SystemRandom` (operating-system CSPRNG) and `SeededRandom` (xoshiro256**, reproducible, for experiments only).
- Prefilter objects are immutable after construction and can be shared between threads.

## Building

Requirements: a C++17 compiler, GMP, and `make` or CMake.

```bash
make && make test                                   # build/libspp.a, the tools, and the unit tests
cmake -B build && cmake --build build && ctest --test-dir build   # alternatively with CMake
```

Linux: `sudo apt install g++ libgmp-dev cmake`; macOS: `brew install gmp cmake`;
Windows: MSYS2 UCRT64 (`pacman -S mingw-w64-ucrt-x86_64-gcc mingw-w64-ucrt-x86_64-gmp mingw-w64-ucrt-x86_64-cmake make`)
or Visual Studio with `vcpkg install gmp:x64-windows` and
`-DCMAKE_TOOLCHAIN_FILE=<vcpkg>/scripts/buildsystems/vcpkg.cmake`.
For a non-standard GMP location use `make GMP_INC=... GMP_LIB=...` or
`cmake -B build -DGMP_INCLUDE_DIR=... -DGMP_LIBRARY=...`.

## Command-line tools

| Tool | Purpose |
|---|---|
| `spp_check` | tests a number with the prefilter + Miller–Rabin and with Miller–Rabin alone, and compares the times |
| `spp_gen` | generates *k* primes with both methods on the same candidate stream |
| `spp_bench` | measures the prefilters on random candidates; `--sweep` finds the best bound for the current machine |
| `spp_example` | short example of the API |

```bash
build/spp_check 2^1279-1
build/spp_check                      # interactive mode
build/spp_gen 1024 1000              # 1000 primes of 1024 bits
build/spp_gen 2048 20 12345          # with a seed
build/spp_bench 1024 10000 --sweep   # calibrate B
```

## Using the API

```cpp
#include <spp/spp.hpp>

// primality test (B chosen automatically); random bases for numbers supplied by others
spp::SystemRandom rng;
spp::BigInt n = spp::parse_expression("2^607-1");
spp::TestReport r = spp::is_probable_prime(n.get(), 40, spp::BaseMode::Random, &rng);
// r.probable_prime, r.stage (Trivial/SmallSieve/Prefilter/MillerRabin), r.factor, r.mr_rounds

// key generation: the prefilter is built once and reused
spp::SegmentedPrefilter pf(spp::recommended_bound(2048));
spp::BigInt p, q;
spp::gen_prime(p.get(), 2048, spp::CandidateSet::S2_TopTwo, -128.0, pf, rng);   // error probability < 2^-128
spp::generate_rsa_primes_fips(p.get(), q.get(), 3072, 65537, pf, rng);
```

`BigInt` is a thin RAII wrapper around `mpz_t`; `.get()` passes it to any GMP or `spp` function.
To use the library from another CMake project:

```cmake
find_package(spp REQUIRED)          # or: add_subdirectory(spp)
target_link_libraries(my_app PRIVATE spp::spp)
```

## Security notes

- Generate key material only with `SystemRandom`; `SeededRandom` is reproducible but not cryptographically secure.
- Use `BaseMode::Random` for numbers supplied by others. `BaseMode::FromCandidate` is meant for reproducible
  comparisons and is not suitable for adversarially chosen inputs.
- The error bounds of `gen_prime()` and `generate_rsa_primes_fips()` are average-case bounds for their own randomly
  drawn candidates; other candidate restrictions require their own bound on *L*(*S*) (Section VI of the paper), and
  validated implementations must follow the minimum round counts of the applicable standard.
- The prefilter is not constant-time. When the timing of key generation can be observed, use a constant-time
  variant such as CT-TREE (`experiments/cpp/ct_experiment.cpp`) and analyze the compiled code for the target compiler
  (`experiments/cpp/ctgrind/`).

## Relation to the paper

| Paper | Library |
|---|---|
| Geometric segmentation, word-size rule (Section IV) | `SegmentedPrefilter(B, s0, r)`, `mpz_gcd_ui` fast path |
| Lemma 1, condition *n* > *p*<sub>*k*</sub> | `Stage::SmallSieve` branch of `is_probable_prime()` |
| Cost-optimal bound | `recommended_bound()`, `spp_bench --sweep` |
| Admissible sets, round counts (Proposition 2, Table 2), GenPrime (Algorithm 2) | `CandidateSet`, `loss_bound()`, `mr_rounds_for_error()`, `gen_prime()` |
| Lemma 2, provable generator | `generate_provable_prime()` |
| RSA primes | `generate_rsa_primes_fips()` |

The unit tests (`tests/spp_tests.cpp`, 174,662 checks) cover the equivalence of the prefilter with trial division
(four bounds, three initial segment sizes, four growth ratios, four bit lengths), the segment sizes and bit lengths
for *B* = 1000 (5, 10, 20, 40, 80, 12 primes; 14, 50, 133, 321, 744, 119 bits), rejection of Carmichael numbers,
strong pseudoprimes and RSA-100, acceptance of Mersenne primes, correct classification of all integers below 5000,
identical outputs of the paired generators, and the complete generators, including provable primes of odd bit
lengths.

## Files

```
include/spp/spp.hpp     public API
src/spp.cpp             implementation
tools/                  spp_check, spp_gen, spp_bench
examples/example.cpp    short example
tests/spp_tests.cpp     unit tests
CMakeLists.txt          CMake build (supports install and find_package)
Makefile                build without CMake
```

## License

MIT (see `../LICENSE`).
