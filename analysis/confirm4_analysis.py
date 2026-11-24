#!/usr/bin/env python3
"""PRE-SPECIFIED analysis of the confirmatory V0-V3 experiment at 2048/3072 bits (protocol of 2026-11-23,
docs/protocol_V0-V3_2048_3072_EN.txt). Only messages and comments were translated from the original Uzbek version;
the computations are unchanged.
Usage:  python3 confirm4_analysis.py conf4_<platform>_<bits>.csv [...]
The platform and the bit length in the file name define the setting. Python standard library only."""
import csv, math, random, sys, os
from statistics import mean, variance
def norm_ppf(p):                       # Acklam's approximation
    a = [-39.6968302866538, 220.946098424521, -275.928510446969, 138.357751867269, -30.6647980661472, 2.50662827745924]
    b = [-54.4760987982241, 161.585836858041, -155.698979859887, 66.8013118877197, -13.2806815528857]
    c = [-7.78489400243029e-03, -0.322396458041136, -2.40075827716184, -2.54973253934373, 4.37466414146497, 2.93816398269878]
    d = [7.78469570904146e-03, 0.32246712907004, 2.445134137143, 3.75440866190742]
    if p < 0.02425:
        q = math.sqrt(-2 * math.log(p)); return (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) / ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1)
    if p > 1 - 0.02425: return -norm_ppf(1 - p)
    q = p - 0.5; r = q * q
    return (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q / (((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1)
pnorm2 = lambda z: math.erfc(abs(z) / math.sqrt(2))          # two-sided p
CONTRASTS = [("H1 complete pipeline V0/V3", "V0_OpenSSL", "V3_independent_WORD"),
             ("H2 prefilter V1/V2", "V1_replica", "V2_hybrid_WORD"),
             ("H3 sampling method V2/V3", "V2_hybrid_WORD", "V3_independent_WORD")]
MARGIN = math.log(1.05)                                       # equivalence margin for V0 and V1 (+-5 %)
def load(paths):
    S = {}
    for p in paths:
        base = os.path.basename(p).replace(".csv", "").split("_")      # conf4_<platform>_<bits>
        key = (base[1], int(base[2]))
        rows = list(csv.DictReader(open(p)))
        S[key] = {m: [float(r["ms"]) for r in rows if r["method"] == m] for m in {r["method"] for r in rows}}
        S[key]["_mr"] = {m: mean(float(r["mr_calls"]) for r in rows if r["method"] == m) for m in {r["method"] for r in rows}}
    return S
def logratio(x, y):
    mx, my = mean(x), mean(y)
    se = math.sqrt(variance(x) / (len(x) * mx * mx) + variance(y) / (len(y) * my * my))
    return math.log(mx / my), se
def boot(x, y, B=int(os.environ.get("BOOT_B", "10000")), seed=20261123):
    r = random.Random(seed); v = []
    for _ in range(B):
        v.append(mean(r.choices(x, k=len(x))) / mean(r.choices(y, k=len(y))))
    v.sort(); return v[int(0.025 * B)], v[int(0.975 * B) - 1]
def welch_p(x, y):
    vx, vy = variance(x) / len(x), variance(y) / len(y); t = (mean(x) - mean(y)) / math.sqrt(vx + vy)
    return pnorm2(t)                                           # normal approximation is adequate for n >= 1000
def fixed(est):
    w = [1 / s ** 2 for _, s in est]; m = sum(wi * e for wi, (e, _) in zip(w, est)) / sum(w); se = math.sqrt(1 / sum(w))
    Q = sum(wi * (e - m) ** 2 for wi, (e, _) in zip(w, est)); df = len(est) - 1
    I2 = max(0.0, (Q - df) / Q) if Q > 0 else 0.0
    return m, se, Q, df, I2
def dersimonian(est):
    m, se, Q, df, _ = fixed(est); w = [1 / s ** 2 for _, s in est]
    tau2 = max(0.0, (Q - df) / (sum(w) - sum(x * x for x in w) / sum(w))) if df > 0 else 0.0
    ws = [1 / (s ** 2 + tau2) for _, s in est]; mr = sum(wi * e for wi, (e, _) in zip(ws, est)) / sum(ws)
    return mr, math.sqrt(1 / sum(ws)), tau2
def main(paths):
    S = load(paths); bits_list = sorted({b for _, b in S}); fam = []
    print("=== Per setting (secondary analysis)")
    for key in sorted(S):
        d = S[key]
        print(f"{key[0]} {key[1]} bit: n = " + ", ".join(f"{m}={len(d[m])}" for m in sorted(k for k in d if not k.startswith('_'))))
        for name, a, b in CONTRASTS:
            e, s = logratio(d[a], d[b]); lo, hi = boot(d[a], d[b])
            print(f"  {name}: {math.exp(e):.4f} [{math.exp(e-1.96*s):.4f}, {math.exp(e+1.96*s):.4f}], Welch p = {welch_p(d[a], d[b]):.4g}, bootstrap [{lo:.4f}, {hi:.4f}]")
        print("  Miller-Rabin calls per prime: " + ", ".join(f"{m}={v:.1f}" for m, v in sorted(d['_mr'].items())))
    print("\n=== Primary endpoints: platforms pooled per bit length (fixed effect), Holm, familywise alpha = 0.05")
    for bits in bits_list:
        keys = [k for k in S if k[1] == bits]
        for name, a, b in CONTRASTS:
            est = [logratio(S[k][a], S[k][b]) for k in keys]
            m, se, Q, df, I2 = fixed(est); p = pnorm2(m / se); mr, ser, tau2 = dersimonian(est)
            fam.append(dict(name=name, bits=bits, R=math.exp(m), lo=math.exp(m - 1.96 * se), hi=math.exp(m + 1.96 * se), p=p, Q=Q, df=df, I2=I2,
                            Rr=math.exp(mr), lor=math.exp(mr - 1.96 * ser), hir=math.exp(mr + 1.96 * ser)))
    order = sorted(range(len(fam)), key=lambda i: fam[i]["p"]); k = len(fam); stop = False
    for rank, i in enumerate(order):
        thr = 0.05 / (k - rank); fam[i]["holm_thr"] = thr
        fam[i]["reject"] = (not stop) and fam[i]["p"] <= thr
        if not fam[i]["reject"]: stop = True
    for f in fam:
        flag = "REJECTED (difference)" if f["reject"] else "not rejected"
        warn = "  [I^2 > 50 %: pooled value reported but interpreted per platform]" if f["I2"] > 0.5 else ""
        print(f"{f['bits']} bit, {f['name']}: R = {f['R']:.4f} [{f['lo']:.4f}, {f['hi']:.4f}], p = {f['p']:.4g}, Holm threshold {f['holm_thr']:.4g} -> {flag};"
              f" Q = {f['Q']:.2f} (df = {f['df']}), I^2 = {100*f['I2']:.0f} %; DL: {f['Rr']:.4f} [{f['lor']:.4f}, {f['hir']:.4f}]{warn}")
    print("\n=== Check of the replica V1: equivalence of V0 and V1 (TOST, margin +-5 %, alpha = 0.05; the 90 % CI must lie within the margin)")
    for key in sorted(S):
        e, s = logratio(S[key]["V0_OpenSSL"], S[key]["V1_replica"]); z = norm_ppf(0.95)
        lo, hi = e - z * s, e + z * s; ok = lo > -MARGIN and hi < MARGIN
        print(f"{key[0]} {key[1]} bit: V0/V1 = {math.exp(e):.4f}, 90 % CI [{math.exp(lo):.4f}, {math.exp(hi):.4f}] -> {'equivalent' if ok else 'equivalence not shown'}")
if __name__ == "__main__":
    if len(sys.argv) < 2: sys.exit(__doc__)
    main(sys.argv[1:])
