#!/usr/bin/env python3
"""Analysis of the replicated seed-dependent counts, of the 22,860 provable primes, and of the complete generators
(Table "complete" of the paper). Bootstrap over seeds, 10,000 resamples, seed 20261123.

    python3 analysis/replication_analysis.py [data/raw/replicates] [data/raw/complete_generators] [OLD]

OLD (optional): directory with diss_bench gen_*.csv / safe_*.csv files to compare the counts of the original seeds.
Writes replication_summary.txt into the first directory."""
import csv, math, random, sys, os
from statistics import mean, stdev
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
R = sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, "data", "raw", "replicates")
CG = sys.argv[2] if len(sys.argv) > 2 else os.path.join(ROOT, "data", "raw", "complete_generators")
OLD = sys.argv[3] if len(sys.argv) > 3 else None
rd = lambda f: list(csv.DictReader(open(os.path.join(R, f))))
def boot(stat, groups, B=10000, seed=20261123):
    r = random.Random(seed); n = len(groups); v = sorted(stat([groups[r.randrange(n)] for _ in range(n)]) for _ in range(B))
    return v[int(0.025 * B)], v[int(0.975 * B) - 1]
out = []
P = lambda s: (print(s), out.append(s))
# ---- 1. provable primes
rows = rd("provodd.csv"); bl = sorted({int(r["bits"]) for r in rows})
P(f"[provodd] {len(rows)} primes, {len(bl)} bit lengths: {bl}")
P(f"  GMP mpz_probab_prime_p: all confirmed = {all(int(r['gmp_check']) > 0 for r in rows)}; 2 (definitely prime) = {sum(int(r['gmp_check']) == 2 for r in rows)}, 1 (probably prime) = {sum(int(r['gmp_check']) == 1 for r in rows)}")
L = sum(int(r["levels"]) for r in rows); Lq = sum(int(r["levels_q_below_sqrt"]) for r in rows)
P(f"  recursion levels {L}; levels with q^2 < p (condition q > sqrt(p) violated) {Lq} ({100*Lq/L:.1f} %)")
top = [r for r in rows]
for b in bl:
    rr = [r for r in rows if int(r["bits"]) == b]
    P(f"  {b} bits: n={len(rr)}, mean {mean(float(r['ms']) for r in rr):.2f} ms, levels {rr[0]['levels']}, share of levels with q^2<p {sum(int(r['levels_q_below_sqrt']) for r in rr)/sum(int(r['levels']) for r in rr):.3f}")
# ---- 2. ordinary primes: GEOM and TD (and MR only)
for f, b in (("genrep_1280.csv", 1280), ("genrep_2560.csv", 2560)):
    if not os.path.exists(os.path.join(R, f)): continue
    rows = rd(f); ks = sorted({int(r["k"]) for r in rows})
    by = {(int(r["k"]), r["method"]): r for r in rows}
    ks = [k for k in ks if (k, "GEOM") in by and (k, "TD") in by and (k, "MR") in by]
    g = [(int(by[k, "GEOM"]["mr_calls"]), int(by[k, "TD"]["mr_calls"])) for k in ks]
    red = lambda G: 1 - sum(a for a, _ in G) / sum(c for _, c in G)
    lo, hi = boot(red, g)
    comp = sum(int(by[k, m]["composite_calls"]) for k in ks for m in ("GEOM", "TD", "MR")); compMR = sum(int(by[k, "MR"]["composite_calls"]) for k in ks)
    beyond = sum(int(by[k, m]["composite_rounds_beyond_first"]) for k in ks for m in ("GEOM", "TD", "MR"))
    P(f"[genrep {b}] {len(ks)} seeds: reduction of Miller-Rabin calls by GEOM relative to TD(1000) {100*red(g):.1f} % (95 % bootstrap {100*lo:.1f}-{100*hi:.1f} %)")
    P(f"  Miller-Rabin calls on composites {comp} (of which in the MR-only baseline {compMR}); rounds beyond the first {beyond}")
    P(f"  mean candidates per prime {mean(int(by[k,'GEOM']['candidates']) for k in ks):.1f}; mean Miller-Rabin calls GEOM/TD/MR {mean(a for a,_ in g):.1f} / {mean(c for _,c in g):.1f} / {mean(int(by[k,'MR']['mr_calls']) for k in ks):.1f}")
    if OLD and os.path.exists(os.path.join(OLD, f.replace("genrep", "gen"))):
        old = {(r["k"], r["method"]): r for r in csv.DictReader(open(os.path.join(OLD, f.replace("genrep", "gen"))))}
        new = {(r["k"], r["method"]): r for r in rows}; n = mm = 0
        for key, r in old.items():
            if key in new: n += 1; mm += sum(r[c] != new[key][c] for c in ("candidates", "mr_calls", "mr_rounds"))
        P(f"  comparison with the original seeds: {n} records, differences {mm}")
# ---- 3. safe primes
for f, b, hl in (("saferep_768.csv", 768, None), ("saferep_1280.csv", 1280, None)):
    if not os.path.exists(os.path.join(R, f)): continue
    rows = rd(f); by = {(int(r["k"]), r["strategy"]): r for r in rows}
    ks = sorted(k for k in {int(r["k"]) for r in rows} if (k, "td_pair") in by and (k, "geom_product") in by)
    g = [(int(by[k, "td_pair"]["modexp"]), int(by[k, "geom_product"]["modexp"]), int(by[k, "geom_product"]["candidates"])) for k in ks]
    rat = lambda G: sum(a for a, _, _ in G) / sum(c for _, c, _ in G); lo, hi = boot(rat, g)
    cm = lambda G: mean(c for _, _, c in G); clo, chi = boot(cm, g)
    x = (b - 1.5) * math.log(2); C2 = 0.6601618158
    P(f"[saferep {b}] {len(ks)} seeds: ratio of modular exponentiations td_pair/geom_product {rat(g):.3f} (95 % bootstrap {lo:.3f}-{hi:.3f})")
    P(f"  mean candidates per safe prime {cm(g):,.0f} ({clo:,.0f}-{chi:,.0f}); Hardy-Littlewood estimate {x*(x+math.log(2))/(4*C2):,.0f}")
    if OLD and os.path.exists(os.path.join(OLD, f.replace("saferep", "safe"))):
        old = {(r["k"], r["strategy"]): r for r in csv.DictReader(open(os.path.join(OLD, f.replace("saferep", "safe"))))}
        new = {(r["k"], r["strategy"]): r for r in rows}; n = mm = 0
        for key, r in old.items():
            if key in new: n += 1; mm += sum(r[c] != new[key][c] for c in ("candidates", "rejected", "modexp"))
        P(f"  comparison with the original seeds: {n} records, differences {mm}")
# ---- 4. complete generators (mean time per prime +- standard error, ms)
if os.path.isdir(CG):
    P("[complete generators] mean time per prime +- standard error (ms); timings are platform dependent")
    for b in (1024, 1536, 2048):
        T = {}
        for f in (f"prob_{b}_tdlp", f"prob_{b}_t64", f"prov_{b}_geom", f"prov_{b}_td"):
            x = [float(r["ms"]) for r in csv.DictReader(open(os.path.join(CG, f + ".csv")))]
            T[f] = (mean(x), stdev(x) / math.sqrt(len(x)), len(x))
        P(f"  {b} bits: " + ", ".join(f"{f} {m:.1f}+-{s:.1f} (n={n})" for f, (m, s, n) in T.items()))
        P(f"    t from Table 2 vs t = 64: {T[f'prob_{b}_t64'][0] / T[f'prob_{b}_tdlp'][0]:.2f}x; provable vs 64 rounds: "
          f"{T[f'prob_{b}_t64'][0] / T[f'prov_{b}_geom'][0]:.2f}x; segmented vs TD within provable: {T[f'prov_{b}_td'][0] / T[f'prov_{b}_geom'][0]:.2f}x")
    for n in (2048, 3072):
        rr = list(csv.DictReader(open(os.path.join(CG, f"rsa_{n}.csv"))))
        P(f"  RSA-{n}: {len(rr)} pairs, t = {rr[0]['t']}, mean {mean(float(r['ms']) for r in rr):.1f} ms, regenerations of q: "
          f"{sum(int(r['pq_rej']) for r in rr)}, top bits ok: {all(r['top_ok'] == '1' for r in rr)}")
open(os.path.join(R, "replication_summary.txt"), "w").write("\n".join(out) + "\n")
