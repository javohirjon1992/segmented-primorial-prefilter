"""Aggregates the raw per-run CSV files written by experiments/cpp/run_all.sh into one JSON summary.

    python analysis/aggregate_raw_results.py experiments/cpp/results      -> experiments/cpp/results/results.json

The aggregated values of both platforms are the basis of the tables in data/tables. The prime-gap and exact
small-size distribution summaries at the end are exploratory and are not reported in the paper.
"""
import json, math, os, sys
import numpy as np, pandas as pd
from scipy import stats
D = sys.argv[1] if len(sys.argv) > 1 else "results"
R = {}
C2 = 0.6601618158; CR = (1 - np.euler_gamma) / np.log(2)
c = pd.read_csv(f"{D}/consts.csv"); R["consts"] = c.round(4).to_dict("records")
k_mr = np.polyfit(np.log(c.bits), np.log(c.c_mr_ns), 1)[0]; k_b = np.polyfit(np.log(c.bits), np.log(c.beta_ns_per_bit), 1)[0]
# ---------- B sweeps and B*
rows = []; curves = {}
def pmin(B, y):
    i = int(np.argmin(y))
    if 0 < i < len(y) - 1:
        a, b_, _ = np.polyfit(np.log10(B[i-1:i+2]), y[i-1:i+2], 2)
        if a > 0: return 10 ** (-b_ / (2 * a))
    return B[i]
def closed(rhs):
    B = rhs / 20
    for _ in range(60): B = rhs / math.log2(B)
    return B
for b in (768, 1280, 2560, 3584):
    s = pd.read_csv(f"{D}/sweep_{b}.csv"); cc = c[c.bits == b].iloc[0]
    s["cost"] = s.geom_ns + s["pass"] * cc.c_mr_ns
    s["td_cost"] = np.where(s.td_ns > 0, s.td_ns + s["pass"] * cc.c_mr_ns, np.nan)
    curves[b] = s
    st = s.dropna(subset=["td_cost"])
    Bm = pmin(s.B.values, s.cost.values); Btd = pmin(st.B.values, st.td_cost.values)
    c1000 = s.loc[(s.B - 1000).abs().idxmin()]
    best = s.loc[s.cost.idxmin()]
    rows.append(dict(bits=b, Bstar=Bm, Bclosed=closed(cc.c_mr_ns / cc.beta_ns_per_bit), Btd=Btd, Btd_closed=cc.c_mr_ns / cc.c_mod_ns,
                     cost1000=c1000.cost, costbest=best.cost, gain=100 * (1 - best.cost / c1000.cost),
                     mrsave=100 * (1 - best["pass"] / c1000["pass"]), tdbest=st.td_cost.min(), geom_vs_td=100 * (1 - best.cost / st.td_cost.min()),
                     edge=bool(best.B == s.B.max())))
bs = pd.DataFrame(rows); R["bstar"] = bs.round(3).to_dict("records")
k_bs = np.polyfit(np.log(bs.bits), np.log(bs.Bstar), 1)[0]
R["scaling"] = dict(k_mr=k_mr, k_beta=k_b, k_Bstar=k_bs, pred=k_mr - k_b)
# ---------- strategies, S0
R["strat"] = pd.concat([pd.read_csv(f"{D}/strat_768.csv"), pd.read_csv(f"{D}/strat_1536.csv")]).to_dict("records")
s0 = pd.read_csv(f"{D}/s0_1280.csv"); R["s0"] = s0.to_dict("records")
# ---------- safe primes
ss = pd.concat([pd.read_csv(f"{D}/safesweep_768.csv"), pd.read_csv(f"{D}/safesweep_1280.csv")])
ss["cost"] = ss.pre_ns + ss["pass"] * ss.cmr_ns
sb = []
for b in (768, 1280):
    for m in ("geom_product", "td_pair"):
        x = ss[(ss.bits == b) & (ss.method == m)]
        sb.append(dict(bits=b, method=m, B=int(x.loc[x.cost.idxmin()].B), costbest=x.cost.min(), cost1000=float(x[(x.B - 1000).abs() < 1].cost.iloc[0]) if (x.B == 1000).any() else float(x.loc[(x.B - 1000).abs().idxmin()].cost)))
R["safe_bstar"] = sb
R["safe_ratio"] = ss[ss.method == "geom_product"].assign(r=lambda d: d.pi2 / d.pi1 ** 2)[["bits", "B", "pass", "pi2", "r"]].round(5).to_dict("records")
sg = pd.concat([pd.read_csv(f"{D}/safe_768.csv"), pd.read_csv(f"{D}/safe_1280.csv")])
R["safe_gen"] = sg.groupby(["bits", "strategy"]).agg(n=("k", "count"), cand=("candidates", "mean"), rej=("rejected", "sum"), candsum=("candidates", "sum"), modexp=("modexp", "mean"), ms=("ms", "mean")).reset_index().round(1).to_dict("records")
# ---------- generation
g = {}
for f in ("gen_1280", "gen_2560", "gen64_1280", "gen128_2560"):
    d = pd.read_csv(f"{D}/{f}.csv"); g[f] = d.groupby("method").agg(n=("k", "count"), cand=("candidates", "sum"), mr=("mr_calls", "sum"), rounds=("mr_rounds", "sum"), ms=("ms", "sum"), med=("ms", "median")).round(2).reset_index().to_dict("records")
R["gen"] = g
o = pd.read_csv(f"{D}/ossl_gen.csv"); R["ossl_gen"] = o.groupby("bits").ms.agg(["count", "mean", "median"]).round(1).reset_index().to_dict("records")
# ---------- distribution (exploratory)
def exact(b):
    lo, hi = 1 << (b - 1), 1 << b
    s = np.ones(hi, bool); s[:2] = False; s[4::2] = False
    for i in range(3, int(hi ** .5) + 1, 2):
        if s[i]: s[i*i::2*i] = False
    P = np.nonzero(s[lo:hi])[0] + lo; m = len(P)
    w = np.empty(m); w[1:] = (P[1:] - P[:-1]) // 2; w[0] = (P[0] - lo - 1) // 2 + 1 + (hi - 1 - P[-1]) // 2
    pw = w / w.sum(); H = -(pw * np.log2(pw)).sum()
    return dict(bits=b, m=m, loss=math.log2(m) - H, ratio=float(w.max() / w.min()), inc_mod3_2=float(w[P % 3 == 2].sum() / w.sum()), ind_mod3_2=float(np.mean(P % 3 == 2)))
R["exact"] = [exact(b) for b in (18, 21, 24)]
ent = []
rng = np.random.default_rng(1)
for b in (768, 1536):
    gg = pd.read_csv(f"{D}/gaps_{b}.csv").gap.values.astype(float)
    f = lambda x: np.mean(np.log2(x)) - np.log2(1 / np.mean(1 / x))
    bt = [f(gg[rng.integers(0, len(gg), len(gg))]) for _ in range(1000)]
    ent.append(dict(bits=b, n=len(gg), gbar=1 / np.mean(1 / gg), bln2=b * math.log(2), loss=f(gg), lo=np.percentile(bt, 2.5), hi=np.percentile(bt, 97.5)))
R["entropy"] = ent
res = pd.read_csv(f"{D}/resid_448.csv"); os3 = pd.read_csv(f"{D}/ossl_resid_448.csv"); os3["method"] = "OpenSSL 3.0.13"
allr = pd.concat([res, os3]); rr = []
for m, d in allr.groupby("method"):
    for s_, col in ((3, "m3"), (5, "m5"), (7, "m7"), (11, "m11"), (13, "m13")):
        cnt = np.array([(d[col] == k).sum() for k in range(1, s_)]); ch = stats.chisquare(cnt)
        rr.append(dict(method=m, s=s_, N=len(d), chi2=ch.statistic, p=ch.pvalue, r1=cnt[0] / len(d), r2=cnt[1] / len(d), r3=(cnt[2] / len(d) if s_ > 3 else None)))
R["resid"] = rr
# ---------- CT
ct = []
for f in ("ct_768_1000", "ct_1280_3000"):
    d = pd.read_csv(f"{D}/{f}.csv"); ct.append(d.to_dict("records"))
R["ct"] = ct
json.dump(R, open(os.path.join(D, "results.json"), "w"), indent=1, default=float)
print(json.dumps(R["scaling"], default=float)); print(bs.round(2).to_string())
print(pd.DataFrame(R["safe_bstar"]).round(2).to_string()); print(pd.DataFrame(R["exact"]).round(4).to_string()); print(pd.DataFrame(ent).round(3).to_string())
print(pd.DataFrame(rr).round(4).to_string())
