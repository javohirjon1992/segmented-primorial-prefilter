# Reproducing the tables and figures

All commands are run from `experiments/cpp/` after `make` (C++ programs) or from the repository root (Python).
Every C++ experiment uses fixed seeds; the seed formulas below are taken from the source code. The counts of
candidates, Miller–Rabin calls and rounds, and modular exponentiations are deterministic functions of the seeds
and are identical on every platform (they were bit-identical on the server and the laptop of the paper);
timings depend on the platform. `./run_all.sh` runs steps 1–15 in order and writes the raw per-run files to
`experiments/cpp/results/`. Table and figure numbers refer to the paper.

## Platform parameters (Table 5)

```bash
./system_info.sh > results/system.txt     # CPU, cores, caches, memory, OS, compilers, GMP, OpenSSL, Python
```
`system_info.sh` is the "platform_info script" named in the protocol of the planned V0–V3 experiment.

## Tables

| Table | Command(s) | Seeds | Output / aggregated table |
|---|---|---|---|
| 2 (round counts), loss factors | `python3 analysis/derive_paper_numbers.py` (Sections 1–2); exact per-set counts: `spp::mr_rounds_for_error(bits, log2_eps, spp::loss_bound(S))` | — | printed |
| 3 (status of guarantees) | — (analytical) | — | — |
| 4 (compiled-code analysis) | `./ctgrind/run_ctgrind.sh` (Valgrind, g++, clang++); `python3 analysis/ctgrind_summary.py` | candidates `424242+bits` | `data/raw/ctgrind/`, `data/tables/ctgrind_code_locations.csv` |
| 6 (platform constants) | `./diss_bench consts 768 1280 1536 2560 3584` | candidates `7+bits` | `data/tables/platform_constants.csv` |
| 7 (survival per segment) | `python experiments/python/prototype_experiments.py` (experiment E2) | prototype `SEED+2` | `data/tables/segment_survival_1024_prototype.csv` |
| 8 (organizations at B = 1000) | `./diss_bench strat 768`, `./diss_bench strat 1536` | candidates `77+bits` | `data/tables/strategies_B1000.csv` |
| 9 (ablation) | `./diss_bench sweep <bits> <n>` for 768/1280 (n = 2000), 2560 (1000), 3584 (600); then `derive_paper_numbers.py` (Section 3) | candidates `2026+bits` | `data/tables/bound_sweep.csv` |
| 10, Fig. 5 (best bounds) | as Table 9; `python3 analysis/optimal_bound_formulas.py`; `derive_paper_numbers.py` (Sections 4, 5, 8); prototype: `experiments/python/03_bound_theory.ipynb` (MODE = "full") | as Table 9; prototype `20261101` | `data/tables/best_bounds.csv`, `best_bounds_prototype.csv` |
| 11 (prime generation) | `./diss_bench gen 1280 60 100000 40`, `./diss_bench gen 2560 12 300000 40` | `31000 + 7·bits + k` | `data/tables/prime_generation.csv` |
| 12 (safe primes) | `./diss_bench safe 768 10 100000`, `./diss_bench safe 1280 3 100000` | `9000 + 100·bits + k` | `data/tables/safe_prime_generation.csv` |
| 11–12, replicated counts | `./review_exp genrep 1280 200 100000 40 1`, `./review_exp genrep 2560 60 300000 40 1`, `./review_exp saferep 768 40 100000`, `./review_exp saferep 1280 40 100000`; `python3 analysis/replication_analysis.py` | same formulas, k = 0..199 / 0..59 / 0..39 | `data/raw/replicates/` |
| 13, Fig. 7 (confirmatory OpenSSL) | `./ossl_seg 1024 1000 10000 confirm`, `./ossl_seg 1536 1000 10000 confirm` on each platform; `python3 analysis/confirmatory_analysis.py <4 files>`; pooling and sensitivity: `derive_paper_numbers.py` (Section 6) | OpenSSL DRBG | `data/tables/openssl_confirmatory.csv` |
| 14 (complete generators) | `./fullgen prob <bits> <K> 0` and `... 64`, `./fullgen prov <bits> <K> 100000` and `... 0`, `./fullgen rsa 2048 30`, `./fullgen rsa 3072 10` (K = 60/40/20 at 1024/1536/2048 bits; `./run_all.sh replicate`); `python3 analysis/replication_analysis.py` | candidates `70000 + 13·bits + k` (prob, prov), `80000 + 7·nlen + k` (rsa); bases `170000 + 13·bits + k`, `180000 + 7·nlen + k` | `data/raw/complete_generators/`, `data/tables/complete_generators.csv` |
| Section VIII-I, 22,860 provable primes | `./review_exp provodd 1270 33 35 39 47 63 65 97 127 129 193 255 385 511 769 1023 1025 1279 1535`; `python3 analysis/replication_analysis.py` | `5000000 + 10000·bits + k` | `data/raw/replicates/provodd.csv` |
| 15 (timing leakage) | `./ct_experiment 768 1000 30000 <csv>`, `./ct_experiment 1280 3000 4000 <csv>` (build with `-DCT_NO_BARRIER` for the version that was timed) | `20261101` | `data/tables/timing_leakage.csv` |

Text-only results: OpenSSL versus the library and RSA keys (`./ossl_bench gen <bits> <K>`; `./rsa_bench <nlen> <K>`, whose
library part uses the seeds `4000 + 10·nlen + k`; `openssl_vs_spp_generation.csv`, `rsa_key_generation.csv`); exploratory OpenSSL
variants V0–V3 (`./ossl_seg <bits> <K> 3000,10000,30000`; `openssl_variants_exploratory.csv`); safe-prime sweep and
batch screening (`./diss_bench safesweep <bits> 6000`, seed `555+bits`; `./safe_batch 1280 4096 30000000`, seed
`777+bits`); residues (`./diss_bench resid 448 10000`, `./ossl_bench resid 448 10000`, `./ossl_seg 448 5000 3000 resid`;
`residue_uniformity.csv`). A run of the complete generators on the server with candidate-derived bases (identical
candidate streams and counts) is kept in `data/tables/complete_generators_server_derived_bases.csv`.

The barrier in the Barrett reduction of CT-TREE (`ct_experiment.cpp`, `ctgrind/ct_impl.hpp`) is on by default in
`ct_experiment.cpp` and is selected with `-DCT_ASM_BARRIER` in the ctgrind build; Clang 18.1 compiles the version
without it into secret-dependent branches, whereas GCC 13.3 produces branch-free code for both versions.

## Figures

```bash
python3 figures/make_figures.py            # Figs. 2-7 from data/tables -> figures/out/*.pdf
pdflatex figures/workflow_diagram.tex      # Fig. 1
```

## Confirmatory V0–V3 experiment at 2048 and 3072 bits (protocol fixed, data not yet collected)

Protocol: `docs/protocol_V0-V3_2048_3072_EN.txt` (original: `..._UZ.txt`). On each platform:
`PLATFORM=server ./run_all.sh plan` (or `./ossl_seg 2048 2000 30000 confirm4 20` and `./ossl_seg 3072 2000 60000 confirm4 20`);
analysis: `python3 analysis/confirm4_analysis.py conf4_<platform>_<bits>.csv ...`.

## Software versions

| Component | Server | Laptop |
|---|---|---|
| Compiler | GCC 13.3 | GCC 16.2 (MSYS2 UCRT64) |
| GMP | 6.3.0 | MSYS2 package |
| OpenSSL | 3.0.13 | 3.6.4 |
| Python (prototype) | — | 3.11.3, gmpy2 2.3.2 |

The replicated counts, the 22,860 provable primes, the complete generators of Table 14, and the compiled-code
analysis were produced in a cloud container (Intel Xeon, GCC 13.3, Clang 18.1, Valgrind 3.22, GMP 6.3.0,
OpenSSL 3.0.13; shared load). The analysis scripts were checked with Python 3.13, NumPy 2.5, SciPy 1.18, pandas 3.0,
and Matplotlib 3.11.
