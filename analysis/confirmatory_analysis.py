"""Analysis of the confirmatory OpenSSL experiment as specified in docs/confirmatory_protocol_EN.txt
(written before data collection). Input: the raw CSV files written by `ossl_seg <bits> 1000 10000 confirm`.
   python3 analysis/confirmatory_analysis.py server_1024.csv server_1536.csv laptop_1024.csv laptop_1536.csv"""
import sys, numpy as np, pandas as pd
from scipy import stats
rng = np.random.default_rng(20261105)
names = ["server 1024", "server 1536", "laptop 1024", "laptop 1536"]
L, W, rows = [], [], []
for name, f in zip(names, sys.argv[1:5]):
    d = pd.read_csv(f); a = d[d.method == "V0_OpenSSL"].ms.values; x = d[d.method == "V3_independent_WORD"].ms.values
    R = a.mean() / x.mean()
    se = np.sqrt((a.std(ddof=1) / np.sqrt(len(a)) / a.mean()) ** 2 + (x.std(ddof=1) / np.sqrt(len(x)) / x.mean()) ** 2)
    bs = [rng.choice(a, len(a)).mean() / rng.choice(x, len(x)).mean() for _ in range(10000)]
    p = stats.ttest_ind(a, x, equal_var=False).pvalue
    mr0 = d[d.method == "V0_OpenSSL"].mr_calls.mean(); mr3 = d[d.method == "V3_independent_WORD"].mr_calls.mean()
    rows.append(dict(setting=name, n0=len(a), n3=len(x), V0_ms=a.mean(), V3_ms=x.mean(), ratio=R, ci_low=np.percentile(bs, 2.5), ci_high=np.percentile(bs, 97.5), welch_p=p, MR_V0=mr0, MR_V3=mr3))
    L.append(np.log(R)); W.append(1 / se ** 2)
L, W = np.array(L), np.array(W); est = (L * W).sum() / W.sum(); se = 1 / np.sqrt(W.sum()); z = est / se
pd.set_option("display.width", 200); print(pd.DataFrame(rows).round(4).to_string(index=False))
Q = ((L - est) ** 2 * W).sum(); print(f"\nPRIMARY ENDPOINT: pooled ratio V0/V3 = {np.exp(est):.4f} (95% CI {np.exp(est - 1.96 * se):.4f}-{np.exp(est + 1.96 * se):.4f}), p = {2 * stats.norm.sf(abs(z)):.4g}")
print(f"Heterogeneity (Cochran Q) = {Q:.2f}, df = 3, p = {stats.chi2.sf(Q, 3):.3f}")
