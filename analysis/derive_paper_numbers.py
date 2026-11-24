"""Recomputes the derived quantities reported in the paper from the CSV tables in data/tables.

    python analysis/derive_paper_numbers.py        (run from the repository root)

Sections
  1. Damgard-Landrock-Pomerance bound p_{b,t} and the round counts of Table "rounds" (eq. tchoice)
  2. Loss factors L(S2) and L(S3) of the admissible candidate sets (Proposition 2)
  3. Ablation: expected cost per candidate relative to trial division at B = 1000
  4. Accuracy of the closed-form and refined approximations of the best bound
  5. Power-law scaling of the best bound with the bit length (95% confidence intervals)
  6. Pooled analysis of the confirmatory OpenSSL experiment (primary endpoint, Cochran's Q, I^2) and its
     sensitivity analysis (DerSimonian-Laird random effects, leave-one-out, pooling by bit length)
  7. Bound (eq. rc) on the mean number of Miller-Rabin rounds r_C of composites that pass the prefilter
  8. Trial division beyond the evaluated grid: extrapolation of eq. (ctd) to its model optimum
Only the Python standard library and NumPy are required.
"""
import csv, math, os
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
T = os.path.join(ROOT, 'data', 'tables')


def load(name):
    with open(os.path.join(T, name)) as f:
        rows = list(csv.DictReader(f))
    for r in rows:
        for k, v in r.items():
            try:
                r[k] = float(v) if v not in ('', None) else None
            except ValueError:
                pass
    return rows


def t_quantile_975(df):
    """Two-sided 95% quantile of Student's t distribution (table values, df <= 30)."""
    tab = {1: 12.706, 2: 4.303, 3: 3.182, 4: 2.776, 5: 2.571, 6: 2.447, 7: 2.365, 8: 2.306,
           9: 2.262, 10: 2.228, 15: 2.131, 20: 2.086, 30: 2.042}
    return tab[df]


def norm_sf(z):
    return 0.5 * math.erfc(z / math.sqrt(2))


def chi2_sf(x, df):
    """Upper tail of the chi-square distribution for integer df (closed-form series)."""
    if df % 2 == 0:
        s, term = 0.0, 1.0
        for i in range(df // 2):
            if i > 0:
                term *= (x / 2) / i
            s += term
        return math.exp(-x / 2) * s
    s, term = 2 * norm_sf(math.sqrt(x)), math.sqrt(2 * x / math.pi) * math.exp(-x / 2)
    for i in range(3, df + 1, 2):
        s += term
        term *= x / i
    return s


print('=' * 78)
print('1. DLP bound log2 p_{b,t} and chosen number of rounds (eps = 2^-128)')
print('=' * 78)


def dlp_log2(b, t):
    """log2 of the Damgard-Landrock-Pomerance bound p_{b,t} (valid for t >= 3, b >= 21, t <= b/9)."""
    return 1.5 * math.log2(b) + t - 0.5 * math.log2(t) + 2 * (2 - math.sqrt(t * b))


L_MAX = 2.05          # upper bound on L(S) for S1, S2, S3 (Section 2 below)
EPS = 2.0 ** -128
print(f"{'b':>5} {'t=3':>8} {'t=4':>8} {'t=5':>8} {'t=6':>8} {'t':>3} {'log2 p':>8}  single  RSA pair")
for b in (768, 1024, 1280, 1536, 2048, 3072):
    t = 3
    # eq. (tchoice) with eps' = eps/2 (RSA pair) and the largest loss factor; this also covers single primes
    while L_MAX * 2 ** dlp_log2(b, t) / (1 - 2 ** dlp_log2(b, t)) > EPS / 2:
        t += 1
    p = 2 ** dlp_log2(b, t)
    ok_single = L_MAX * p / (1 - p) <= EPS
    ok_pair = L_MAX * p / (1 - p) <= EPS / 2
    print(f"{b:>5} " + ' '.join(f"{dlp_log2(b, s):8.1f}" for s in (3, 4, 5, 6)) +
          f" {t:>3} {dlp_log2(b, t):8.1f}  {ok_single!s:>6}  {ok_pair!s:>8}")

print()
print('=' * 78)
print('2. Loss factors (Dusart: x/(ln x - 1) <= pi(x) <= x/(ln x - 1.1))')
print('=' * 78)


def pi_lo(lnx):      # lower bound on pi(x)/x, x >= 5393
    return 1.0 / (lnx - 1.0)


def pi_hi(lnx):      # upper bound on pi(x)/x, x >= 60184
    return 1.0 / (lnx - 1.1)


E = 65537
print(f"{'b':>5} {'L(S2) <=':>10} {'L(S3) <=':>10}")
worst2 = worst3 = 0.0
for b in range(512, 4097, 64):
    lx = b * math.log(2)                     # all quantities in units of 2^b
    all_hi = pi_hi(lx) - 0.5 * pi_lo(lx + math.log(0.5))          # Pi(O_b)  <= ...
    top_lo = pi_lo(lx) - 0.75 * pi_hi(lx + math.log(0.75))        # Pi(S2)   >= ...
    removed = 0.25 / (2 * E) + 2.0 ** -b                          # n = 1 (mod 2e) in [0.75*2^b, 2^b)
    L2, L3 = all_hi / top_lo, all_hi / (top_lo - removed)
    worst2, worst3 = max(worst2, L2), max(worst3, L3)
    if b in (512, 768, 1024, 2048, 3072, 4096):
        print(f"{b:>5} {L2:10.4f} {L3:10.4f}")
print(f"max over 512..4096: L(S2) <= {worst2:.4f}, L(S3) <= {worst3:.4f}")

print()
print('=' * 78)
print('3. Ablation: reduction of expected cost (%) relative to trial division at B = 1000')
print('=' * 78)
sw = load('bound_sweep.csv')
print(f"{'platform':>8} {'bits':>5} {'segm B=1000':>12} {'TD best':>16} {'segm best':>16}")
for plat in ('server', 'laptop'):
    for b in (768, 1280, 2560, 3584):
        rows = [r for r in sw if r['platform'] == plat and r['bits'] == b]
        cmr = rows[0]['c_MR_us']
        B = np.array([r['B'] for r in rows])
        P = np.array([r['pass_fraction_measured'] for r in rows])
        G = np.array([r['segmented_prefilter_us'] for r in rows])
        Td = np.array([r['trial_division_us'] if r['trial_division_us'] is not None else np.nan for r in rows])
        cg, ct = G + P * cmr, Td + P * cmr
        i0 = list(B).index(1000)
        ig, it = int(np.argmin(cg)), int(np.nanargmin(ct))
        base = ct[i0]
        edge = ' (edge)' if it == int(np.max(np.where(~np.isnan(ct)))) else ''
        print(f"{plat:>8} {b:>5} {100 * (base - cg[i0]) / base:12.1f} "
              f"{100 * (base - ct[it]) / base:6.1f} ({B[it]:.1e}){edge:<7} {100 * (base - cg[ig]) / base:6.1f} ({B[ig]:.1e})")

print()
print('=' * 78)
print('4. Accuracy of the approximations of the best bound (Table "bstar")')
print('=' * 78)
ob = load('best_bounds.csv')
for r in ob:
    print(f"{r['platform']:>8} {int(r['bits']):>5}  B*={r['Bstar_measured']:>9.0f}  B*_TD={r['Bstar_TD_measured']:>7.0f}  "
          f"ratio={r['Bstar_measured'] / r['Bstar_TD_measured']:5.1f}  closed/meas={r['closed_over_measured']:.2f}  "
          f"refined/meas={r['refined_over_measured']:.2f}")
lc = [abs(math.log(r['closed_over_measured'])) for r in ob]
lr = [abs(math.log(r['refined_over_measured'])) for r in ob]
print(f"mean |log error|: closed {np.mean(lc):.2f} (factor {math.exp(np.mean(lc)):.2f}), "
      f"refined {np.mean(lr):.2f} (factor {math.exp(np.mean(lr)):.2f})")
td = [r['Bstar_TD_formula'] / r['Bstar_TD_measured'] for r in ob if r['Bstar_TD_measured'] < 99999]
print(f"c_MR/c_mod over measured B*_TD (excluding grid edge): {min(td):.2f}-{max(td):.2f}")

print()
print('=' * 78)
print('5. Power-law fits B* ~ b^k with 95% confidence intervals')
print('=' * 78)


def fit(X, y):
    X = np.asarray(X, float)
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    res = y - X @ beta
    df = len(y) - X.shape[1]
    s2 = res @ res / df
    cov = s2 * np.linalg.inv(X.T @ X)
    return beta, np.sqrt(np.diag(cov)), df, res


for plat in ('server', 'laptop'):
    rows = [r for r in ob if r['platform'] == plat]
    x = np.log([r['bits'] for r in rows]); y = np.log([r['Bstar_measured'] for r in rows])
    beta, se, df, res = fit(np.c_[np.ones_like(x), x], y)
    q = t_quantile_975(df)
    print(f"{plat:>8}: k = {beta[1]:.2f} (95% CI {beta[1] - q * se[1]:.2f}-{beta[1] + q * se[1]:.2f}), "
          f"max residual factor {math.exp(np.max(np.abs(res))):.2f}")
x = np.log([r['bits'] for r in ob]); y = np.log([r['Bstar_measured'] for r in ob])
d = np.array([1.0 if r['platform'] == 'laptop' else 0.0 for r in ob])
beta, se, df, res = fit(np.c_[np.ones_like(x), d, x], y)
q = t_quantile_975(df)
print(f"   joint: k = {beta[2]:.2f} (95% CI {beta[2] - q * se[2]:.2f}-{beta[2] + q * se[2]:.2f}), "
      f"platform-specific intercepts, max residual factor {math.exp(np.max(np.abs(res))):.2f}")
op = load('best_bounds_prototype.csv')
x = np.log([r['bits'] for r in op]); y = np.log([r['Bstar_measured'] for r in op])
beta, se, df, res = fit(np.c_[np.ones_like(x), x], y)
q = t_quantile_975(df)
print(f"prototype: k = {beta[1]:.2f} (95% CI {beta[1] - q * se[1]:.2f}-{beta[1] + q * se[1]:.2f}), "
      f"{len(op)} bit lengths, max residual factor {math.exp(np.max(np.abs(res))):.2f}")
# The exponent 2.31 quoted for the prototype (c_MR ~ b^3.02, kappa ~ b^0.71) is computed from the
# prototype's own cost constants, which experiments/python/prototype_experiments.py measures and writes
# to its summary file (keys E5_*); the C++ constants are listed in data/tables/platform_constants.csv.

print()
print('=' * 78)
print('6. Confirmatory OpenSSL experiment: fixed-effect pooling of log(V0/V3)')
print('=' * 78)
cf = [r for r in load('openssl_confirmatory.csv') if r['primes_per_variant']]
lr = np.array([math.log(r['ratio_V0_over_V3']) for r in cf])
se = np.array([r['se_log_ratio'] for r in cf])
w = 1 / se ** 2
m = (w * lr).sum() / w.sum(); sm = 1 / math.sqrt(w.sum())
p = 2 * norm_sf(abs(m / sm))
Q = (w * (lr - m) ** 2).sum()
for r, l, s in zip(cf, lr, se):
    print(f"{r['setting']:<32} ratio {math.exp(l):.3f} [{math.exp(l - 1.96 * s):.3f}, {math.exp(l + 1.96 * s):.3f}]  "
          f"Welch p = {r['welch_p']:.4f}  fewer MR calls {100 * (1 - r['V3_mr_calls'] / r['V0_mr_calls']):.1f}%")
print(f"pooled ratio {math.exp(m):.3f} [{math.exp(m - 1.96 * sm):.3f}, {math.exp(m + 1.96 * sm):.3f}], p = {p:.3f}; "
      f"Cochran Q = {Q:.1f} (df {len(cf) - 1}), p = {chi2_sf(Q, len(cf) - 1):.3f}")
I2 = max(0.0, (Q - (len(cf) - 1)) / Q)
tau2 = max(0.0, (Q - (len(cf) - 1)) / (w.sum() - (w ** 2).sum() / w.sum()))
wr = 1 / (se ** 2 + tau2); mr = (wr * lr).sum() / wr.sum(); sr = 1 / math.sqrt(wr.sum())
print(f"I^2 = {100 * I2:.0f}%; random effects (DerSimonian-Laird): {math.exp(mr):.3f} "
      f"[{math.exp(mr - 1.96 * sr):.3f}, {math.exp(mr + 1.96 * sr):.3f}], p = {2 * norm_sf(abs(mr / sr)):.2f}")
loo = []
for i in range(len(cf)):
    j = [k for k in range(len(cf)) if k != i]
    loo.append(math.exp((w[j] * lr[j]).sum() / w[j].sum()))
print(f"leave-one-out fixed-effect estimates: {min(loo):.3f}-{max(loo):.3f}")
for bits in ('1024', '1536'):
    j = [k for k, r in enumerate(cf) if bits in r['setting']]
    mb = (w[j] * lr[j]).sum() / w[j].sum(); sb = 1 / math.sqrt(w[j].sum()); Qb = (w[j] * (lr[j] - mb) ** 2).sum()
    print(f"{bits} bits pooled: {math.exp(mb):.3f} [{math.exp(mb - 1.96 * sb):.3f}, {math.exp(mb + 1.96 * sb):.3f}], "
          f"p = {2 * norm_sf(abs(mb / sb)):.2g}; Q = {Qb:.2f} (df 1), p = {chi2_sf(Qb, 1):.2f}")

print()
print('=' * 78)
print('7. Bound on r_C - 1 (eq. rc): (4/3) p_{b,1}/(1-p_{b,1}) rho_b/(W(B)-rho_b), W(B) = 2 e^-gamma / ln B')
print('=' * 78)
GAMMA = 0.5772156649015329
for b in (512, 768, 1024, 2048):
    rho = 2 / (b * math.log(2)); p1 = 2 ** (2 * math.log2(b) + 2 * (2 - math.sqrt(b)))
    ex = [(4 / 3) * p1 / (1 - p1) * rho / (2 * math.exp(-GAMMA) / math.log(B) - rho) for B in (1e3, 1e5, 1e7)]
    print(f"b = {b:4d}: r_C - 1 <= " + ', '.join(f"{e:.2e} (B=1e{k})" for e, k in zip(ex, (3, 5, 7))))

print()
print('=' * 78)
print('8. Trial division beyond the grid (model extrapolation of eq. ctd with the measured c_mod)')
print('=' * 78)
N = 2_000_000
sv = np.ones(N + 1, bool); sv[:2] = False; sv[4::2] = False
for i in range(3, int(N ** .5) + 1, 2):
    if sv[i]:
        sv[i * i::2 * i] = False
PR = np.nonzero(sv)[0]; PR = PR[PR > 2]
pc = load('platform_constants.csv')
for plat in ('server', 'laptop'):
    for b in (768, 1280, 2560, 3584):
        rows = [r for r in sw if r['platform'] == plat and r['bits'] == b]
        cmr = rows[0]['c_MR_us']; cmod = [r for r in pc if r['bits'] == b][0][f'{plat}_c_mod_ns'] / 1000
        B = np.array([r['B'] for r in rows]); Pf = np.array([r['pass_fraction_measured'] for r in rows])
        G = np.array([r['segmented_prefilter_us'] for r in rows])
        Td = np.array([r['trial_division_us'] if r['trial_division_us'] is not None else np.nan for r in rows])
        ct = Td + Pf * cmr; last = int(np.max(np.where(~np.isnan(ct))))
        wk, cost, best, bestB = Pf[last], Td[last], ct[last], B[last]
        for p in PR[PR > B[last]]:
            cost += wk * cmod; wk *= 1 - 1 / p
            if cost + wk * cmr < best:
                best, bestB = cost + wk * cmr, p
        gbest = (G + Pf * cmr).min()
        if bestB > B[last]:
            print(f"{plat:>8} {b:>5}: grid optimum of TD at the edge; extrapolated optimum {bestB:.3g} "
                  f"(c_MR/c_mod = {cmr / cmod:.3g}); gain of the segmented prefilter {100 * (1 - gbest / np.nanmin(ct)):.1f}% "
                  f"-> {100 * (1 - gbest / best):.1f}%")

