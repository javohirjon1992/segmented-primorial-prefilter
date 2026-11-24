"""Closed-form (eq. closed) and refined (eq. refined) approximations of the best bound, evaluated with the
measured platform constants and compared with the measured best bounds (Table "bstar" of the paper).

    python analysis/optimal_bound_formulas.py
Expects results/server/{consts.csv,results.json} and results/laptop/{consts.csv,results.json}, i.e. the raw
results of run_all.sh and the output of aggregate_raw_results.py for each platform.
"""
import numpy as np, pandas as pd, json, math
N = 5_000_000
s = np.ones(N + 1, bool); s[:2] = False; s[4::2] = False
for i in range(3, int(N ** .5) + 1, 2):
    if s[i]: s[i*i::2*i] = False
P = np.nonzero(s)[0]; P = P[P > 2]
logS = np.concatenate([[0.0], np.cumsum(np.log(1 - 1.0 / P))])      # logS[k] = ln prod_{i<k}(1-1/p_i)
sizes, starts, i, j = [], [], 0, 0
while i < len(P):
    z = round(5 * 2 ** j); sizes.append(z); starts.append(i); i += z; j += 1
starts = np.array(starts); sizes = np.array(sizes)
def seg_info(B):
    k = int(np.searchsorted(P, B, side="right"))            # primes <= B
    g = int(np.searchsorted(starts, k - 1, side="right") - 1)  # segment of the k-th prime
    R = math.exp(logS[k] - logS[starts[g]]) ** -1           # S_start / S_k
    return sizes[g], R
def solve(f, lo=50, hi=4e7):
    for _ in range(200):
        m = math.sqrt(lo * hi)
        if f(m) > 0: hi = m
        else: lo = m
    return math.sqrt(lo * hi)
rows = []
for plat, cf, bf in (("server", "results/server/consts.csv", "results/server/results.json"),
                     ("laptop", "results/laptop/consts.csv", "results/laptop/results.json")):
    c = pd.read_csv(cf).set_index("bits"); bs = pd.DataFrame(json.load(open(bf))["bstar"]).set_index("bits")
    for b in (768, 1280, 2560, 3584):
        cmr, beta, alpha, cmod = c.loc[b, "c_mr_ns"], c.loc[b, "beta_ns_per_bit"], c.loc[b, "alpha_ns"], c.loc[b, "c_mod_ns"]
        simple = solve(lambda B: B * beta * math.log2(B) - cmr)
        refined = solve(lambda B: B * (beta * math.log2(B) + alpha / seg_info(B)[0]) * seg_info(B)[1] - cmr)
        sB, RB = seg_info(refined)
        rows.append(dict(platform=plat, bits=b, measured=bs.loc[b, "Bstar"], simple=simple, refined=refined,
                         ratio_simple=simple / bs.loc[b, "Bstar"], ratio_refined=refined / bs.loc[b, "Bstar"],
                         R=RB, seg=sB, td_measured=bs.loc[b, "Btd"], td_formula=cmr / cmod))
d = pd.DataFrame(rows); pd.set_option("display.width", 200); print(d.round(3).to_string())
print("simple ratio range", d.ratio_simple.min().round(2), d.ratio_simple.max().round(2), "| refined", d.ratio_refined.min().round(2), d.ratio_refined.max().round(2))
print("mean |log ratio| simple", np.mean(np.abs(np.log(d.ratio_simple))).round(3), "refined", np.mean(np.abs(np.log(d.ratio_refined))).round(3))
d.to_csv("theory_bounds.csv", index=False)
