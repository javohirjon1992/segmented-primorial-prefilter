# Raw data

* `replicates/` — raw per-seed outputs of `review_exp` (`./run_all.sh replicate`): `genrep_*.csv` (ordinary primes,
  Miller–Rabin calls and rounds), `saferep_*.csv` (safe primes, modular exponentiations), and `provodd.csv`
  (22,860 provable primes at 18 odd bit lengths, recursion levels with q² < p, check with GMP's
  `mpz_probab_prime_p`). The counts are deterministic functions of the seeds; the `ms` columns are platform-dependent
  timings from a cloud container under shared load. Analysis: `analysis/replication_analysis.py`
  (writes `replication_summary.txt`).
* `complete_generators/` — raw outputs of `fullgen` for Table 14 (`prob_<bits>_tdlp` with t from Table 2,
  `prob_<bits>_t64`, `prov_<bits>_geom` and `prov_<bits>_td`, `rsa_<nlen>`).
* `ctgrind/` — results and Valgrind logs of the compiled-code analysis (`experiments/cpp/ctgrind/run_ctgrind.sh`);
  `analysis/ctgrind_summary.py` unpacks the logs and counts the code locations of Table 4.

The raw per-run files of the timing experiments on the server and the laptop are written by
`experiments/cpp/run_all.sh` to `experiments/cpp/results/` (one directory per platform) and summarized by
`analysis/aggregate_raw_results.py`; their aggregates are the tables in `data/tables/`.
