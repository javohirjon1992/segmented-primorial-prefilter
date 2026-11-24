# %% [markdown]
# # Segmented primorial prefilter: prototype experiments (Python, CPython integers and gmpy2/GMP)
#
# Experiments of the prototype used in the paper for correctness checks and for validating the
# probabilistic and cost models over wide parameter grids:
#
# | # | Experiment |
# |---|---|
# | E0 | Environment (CPU, Python, libraries) |
# | E1 | Correctness: identical decisions of all prefilters, reference Miller-Rabin test, condition n > p_k |
# | E2 | Survival through the segments: measurement versus Mertens products |
# | E3 | Prefilter strategies with bootstrap confidence intervals and Wilcoxon tests |
# | E4 | Grid over the initial segment size S0 and the growth ratio r |
# | E5 | Cost model: GCD cost as a function of the product length |
# | E6 | Optimization of the sieving bound B |
# | E7 | Complete primality-testing pipeline, 512-4096 bits |
# | E8 | Timing noise in end-to-end measurements |
# | E9 | Prime generation: independent, incremental, sieve-based |
#
# Usage: choose MODE in the first cell (`smoke`: 1-2 min, `quick`: 15-40 min, `full`: settings of the
# paper, several hours), then run all cells. Results are written to `results/` and `results.zip`.
# Prototype timings include interpreter overhead and are not used for performance claims.

# %%
# 0. Installation (once). Uncomment if needed:
# %pip install gmpy2 sympy scipy pandas matplotlib numpy

# %%
# 1. CONFIGURATION
MODE = "quick"          # "smoke" | "quick" | "full"
SEED = 20261101
OUTDIR = "results"
USE_GMPY2 = True        # used if installed
MR_ROUNDS = 40          # Miller-Rabin rounds of the pipeline experiments
MR_ROUNDS_SENS = [1, 5, 10, 20, 40, 64]   # for the theoretical speedup bound

PRESETS = {
    "smoke": dict(
        R=3, N_PREF=400, BITS_PREF=[1024],
        S0_GRID=[2, 4, 5, 8], RATIO_GRID=[2.0], R_GRID=3,
        B_GRID=[1000, 10000], TD_MAX_B=10000, BITS_B=[1024], N_B=300, R_B=3,
        BITS_PIPE=[512, 1024], N_PIPE={512: 400, 1024: 300},
        N_E2E=200, R_E2E=2, PRIMEGEN={512: 3}, BOOT=500,
    ),
    "quick": dict(
        R=15, N_PREF=3000, BITS_PREF=[1024, 2048],
        S0_GRID=list(range(1, 33)), RATIO_GRID=[1.5, 2.0, 3.0, 4.0], R_GRID=9,
        B_GRID=[300, 1000, 3000, 10000, 30000, 100000, 300000], TD_MAX_B=100000,
        BITS_B=[512, 1024, 2048, 4096], N_B=1500, R_B=5,
        BITS_PIPE=[512, 1024, 2048, 3072, 4096],
        N_PIPE={512: 5000, 1024: 3000, 2048: 1000, 3072: 500, 4096: 300},
        N_E2E=1000, R_E2E=5, PRIMEGEN={1024: 10, 2048: 4}, BOOT=2000,
    ),
    "full": dict(
        R=31, N_PREF=10000, BITS_PREF=[512, 1024, 2048, 4096],
        S0_GRID=list(range(1, 33)), RATIO_GRID=[1.5, 2.0, 3.0, 4.0], R_GRID=21,
        B_GRID=[300, 1000, 3000, 10000, 30000, 100000, 300000, 1000000], TD_MAX_B=300000,
        BITS_B=[512, 1024, 2048, 3072, 4096], N_B=4000, R_B=11,
        BITS_PIPE=[512, 1024, 2048, 3072, 4096],
        N_PIPE={512: 20000, 1024: 10000, 2048: 4000, 3072: 2000, 4096: 1000},
        N_E2E=3000, R_E2E=11, PRIMEGEN={1024: 50, 2048: 20, 3072: 6}, BOOT=5000,
    ),
}
CFG = PRESETS[MODE]
print("MODE =", MODE)
for k, v in CFG.items():
    print(f"  {k} = {v}")

# %%
# 2. IMPORTS AND HELPER FUNCTIONS
import os, sys, gc, json, math, time, random, platform, zipfile, datetime, itertools, warnings
from time import perf_counter_ns
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Polygon

warnings.filterwarnings("ignore")
try:
    from scipy import stats as sstats
    HAVE_SCIPY = True
except Exception:
    HAVE_SCIPY = False
try:
    import sympy
    HAVE_SYMPY = True
except Exception:
    HAVE_SYMPY = False
try:
    import gmpy2
    HAVE_GMPY2 = True and USE_GMPY2
except Exception:
    HAVE_GMPY2 = False

os.makedirs(f"{OUTDIR}/tables", exist_ok=True)
os.makedirs(f"{OUTDIR}/figures", exist_ok=True)
plt.rcParams.update({"figure.dpi": 110, "savefig.dpi": 300, "font.size": 10,
                     "axes.grid": True, "grid.alpha": 0.3, "font.family": "DejaVu Serif"})
RNG = np.random.default_rng(SEED)
REPORT = []     # (kind, name, caption, payload)
SUMMARY = {"mode": MODE, "seed": SEED, "mr_rounds": MR_ROUNDS}
T_START = time.time()

def log(msg):
    print(f"[{time.time() - T_START:7.1f}s] {msg}", flush=True)

def fmt(x, sig=4):
    if isinstance(x, (int, np.integer)):
        return f"{x}"
    if isinstance(x, (float, np.floating)):
        if not np.isfinite(x):
            return str(x)
        if x == 0:
            return "0"
        if abs(x) >= 1e5 or abs(x) < 1e-3:
            return f"{x:.{sig-1}e}"
        return f"{x:.{sig}g}" if abs(x) >= 1 else f"{x:.{sig}f}".rstrip("0")
    return str(x)

def df_to_md(df):
    cols = [str(c) for c in df.columns]
    out = ["| " + " | ".join(cols) + " |", "|" + "|".join(["---"] * len(cols)) + "|"]
    for _, row in df.iterrows():
        out.append("| " + " | ".join(fmt(v) for v in row.values) + " |")
    return "\n".join(out)

def df_to_latex(df, caption, label):
    esc = lambda s: str(s).replace("_", r"\_").replace("%", r"\%").replace("&", r"\&").replace("#", r"\#")
    cols = "l" + "r" * (len(df.columns) - 1)
    lines = [r"\begin{table}[ht]", r"\centering", rf"\caption{{{esc(caption)}}}", rf"\label{{{label}}}",
             rf"\begin{{tabular}}{{{cols}}}", r"\hline",
             " & ".join(esc(c) for c in df.columns) + r" \\", r"\hline"]
    for _, row in df.iterrows():
        lines.append(" & ".join(esc(fmt(v)) for v in row.values) + r" \\")
    lines += [r"\hline", r"\end{tabular}", r"\end{table}"]
    return "\n".join(lines)

def save_table(df, name, caption):
    df.to_csv(f"{OUTDIR}/tables/{name}.csv", index=False)
    with open(f"{OUTDIR}/tables/{name}.tex", "w", encoding="utf-8") as f:
        f.write(df_to_latex(df, caption, f"tab:{name}"))
    with open(f"{OUTDIR}/tables/{name}.md", "w", encoding="utf-8") as f:
        f.write(f"**{caption}**\n\n" + df_to_md(df) + "\n")
    REPORT.append(("table", name, caption, df_to_md(df)))
    try:
        from IPython.display import display
        print(caption); display(df)
    except Exception:
        print(caption); print(df.to_string())

def save_fig(fig, name, caption):
    fig.tight_layout()
    fig.savefig(f"{OUTDIR}/figures/{name}.png", bbox_inches="tight")
    fig.savefig(f"{OUTDIR}/figures/{name}.pdf", bbox_inches="tight")
    REPORT.append(("figure", name, caption, f"figures/{name}.png"))
    plt.show()
    plt.close(fig)

# --- statistics ---
def boot_ci_median(x, nb=None):
    nb = nb or CFG["BOOT"]
    x = np.asarray(x, float)
    idx = RNG.integers(0, len(x), (nb, len(x)))
    m = np.median(x[idx], axis=1)
    return float(np.percentile(m, 2.5)), float(np.percentile(m, 97.5))

def boot_ratio(a, b, nb=None):
    """CI of median(a)/median(b), a and b resampled independently."""
    nb = nb or CFG["BOOT"]
    a, b = np.asarray(a, float), np.asarray(b, float)
    ia = RNG.integers(0, len(a), (nb, len(a))); ib = RNG.integers(0, len(b), (nb, len(b)))
    r = np.median(a[ia], axis=1) / np.median(b[ib], axis=1)
    return float(np.median(a) / np.median(b)), float(np.percentile(r, 2.5)), float(np.percentile(r, 97.5))

def wilcoxon_p(a, b):
    if not HAVE_SCIPY or len(a) < 5:
        return float("nan")
    try:
        return float(sstats.wilcoxon(a, b).pvalue)
    except Exception:
        return float("nan")

# %%
# 3. ARITHMETIC BACKENDS, PRIMES, CANDIDATES
class Backend:
    def __init__(self, name):
        self.name = name
        if name == "gmpy2":
            self.mpz, self.gcd, self.powmod = gmpy2.mpz, gmpy2.gcd, gmpy2.powmod
            self.version = gmpy2.version() + " / GMP " + gmpy2.mp_version()
        else:
            self.mpz, self.gcd, self.powmod = int, math.gcd, pow
            self.version = "CPython int " + platform.python_version()

BACKENDS = [Backend("python")] + ([Backend("gmpy2")] if HAVE_GMPY2 else [])
MRB = BACKENDS[-1]          # Miller-Rabin uses the fastest available backend
log(f"Prefilter backends: {[b.name for b in BACKENDS]};  MR backend: {MRB.name}")

def odd_primes_upto(B):
    """Odd primes 3..B (2 excluded: candidates are odd)."""
    s = np.ones(B + 1, dtype=bool); s[:2] = False
    for i in range(2, int(B ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = False
    return [int(p) for p in np.nonzero(s)[0] if p > 2]

_PCACHE = {}
def primes(B):
    if B not in _PCACHE:
        _PCACHE[B] = odd_primes_upto(B)
    return _PCACHE[B]

def survival_odd(ps):
    """Mertens-type density of odd integers with no factor in ps: prod(1-1/p)."""
    return float(np.prod([1 - 1 / p for p in ps])) if ps else 1.0

def prime_density_odd(bits):
    """Probability that a random odd bits-bit integer is prime ~ 2/ln(2^bits)."""
    return 2.0 / (bits * math.log(2))

def candidates(bits, n, seed):
    r = random.Random(seed)
    return [r.getrandbits(bits) | (1 << (bits - 1)) | 1 for _ in range(n)]

# %%
# 4. PREFILTERS (all over the same set of primes P = {3..B}; precondition: n > B)
def sizes_fixed(P, k):
    return [k] * math.ceil(len(P) / k)

def sizes_bits(P, target_bits):
    sizes, cur, prod = [], 0, 1
    for p in P:
        if cur and (prod * p).bit_length() > target_bits:
            sizes.append(cur); cur, prod = 0, 1
        cur += 1; prod *= p
    if cur:
        sizes.append(cur)
    return sizes

def sizes_geom(P, s0, ratio=2.0):
    sizes, tot, i = [], 0, 0
    while tot < len(P):
        s = max(1, int(round(s0 * ratio ** i)))
        sizes.append(s); tot += s; i += 1
    return sizes

def make_blocks(P, sizes):
    blocks, i = [], 0
    for s in sizes:
        b = P[i:i + s]
        if not b:
            break
        blocks.append(b); i += s
    return blocks

def make_segmented(P, sizes, be):
    mods = [be.mpz(math.prod(b)) for b in make_blocks(P, sizes)]
    g = be.gcd
    def f(n):
        for M in mods:
            if g(n, M) != 1:
                return False
        return True
    f.mods = mods
    return f

def make_trial(P, be):
    def f(n):
        for p in P:
            if n % p == 0:
                return False
        return True
    return f

def make_single(P, be):
    M = be.mpz(math.prod(P)); g = be.gcd
    def f(n):
        return g(n, M) == 1
    return f

def make_psm(P, be, word_bits=64):
    """GMP-style: reduce n modulo word-sized primorial chunks, then test small residue."""
    chunks, cur, prod = [], [], 1
    for p in P:
        if prod * p >= (1 << word_bits):
            chunks.append((be.mpz(prod), cur)); cur, prod = [], 1
        cur.append(p); prod *= p
    if cur:
        chunks.append((be.mpz(prod), cur))
    def f(n):
        for M, ps in chunks:
            r = int(n % M)
            for p in ps:
                if r % p == 0:
                    return False
        return True
    return f

def make_null(P, be):
    def f(n):
        return True
    return f

def block_of_rejection(n, blocks_int):
    """Index of first block sharing a factor with n, or -1 if n passes (instrumented, not timed)."""
    for i, M in enumerate(blocks_int):
        if math.gcd(n, M) != 1:
            return i
    return -1

def method_suite(B, be, s0_list=(4, 5)):
    P = primes(B)
    d = {
        "TD (trial division)": make_trial(P, be),
        "SP (single primorial)": make_single(P, be),
        "PSM (word-primorial + small mod, GMP-style)": make_psm(P, be),
        "SEG-fixed16": make_segmented(P, sizes_fixed(P, 16), be),
        "SEG-bits256": make_segmented(P, sizes_bits(P, 256), be),
    }
    for s0 in s0_list:
        d[f"SEG-geom(S0={s0},r=2)"] = make_segmented(P, sizes_geom(P, s0, 2.0), be)
    return d

# %%
# 5. MILLER-RABIN AND TIMING
def make_mr(be, t=MR_ROUNDS, seed=SEED):
    pm = be.powmod
    def mr(n):
        """Returns (is_probable_prime, rounds_executed). Bases are deterministic per n."""
        m = int(n) - 1
        s = (m & -m).bit_length() - 1
        d = be.mpz(m >> s)
        N = be.mpz(n); Nm1 = N - 1
        rnd = random.Random(seed ^ (int(n) & ((1 << 61) - 1)))
        rounds = 0
        for _ in range(t):
            a = be.mpz(rnd.randrange(2, int(n) - 1))
            x = pm(a, d, N); rounds += 1
            if x == 1 or x == Nm1:
                continue
            for _ in range(s - 1):
                x = x * x % N
                if x == Nm1:
                    break
            else:
                return False, rounds
        return True, rounds
    return mr

def time_once(f, cands):
    t0 = perf_counter_ns()
    for n in cands:
        f(n)
    return (perf_counter_ns() - t0) / len(cands)

def bench_interleaved(methods, cands, R, warmup=1):
    """ns/candidate for each method; method order shuffled each repetition (paired design)."""
    names = list(methods)
    res = {k: [] for k in names}
    for _ in range(warmup):
        for k in names:
            time_once(methods[k], cands)
    order_rng = random.Random(SEED + 7)
    gc_was = gc.isenabled(); gc.disable()
    try:
        for _ in range(R):
            order_rng.shuffle(names)
            for k in names:
                res[k].append(time_once(methods[k], cands))
    finally:
        if gc_was:
            gc.enable()
    return {k: np.array(v) for k, v in res.items()}

def summarize_times(res, baseline):
    rows = []
    base = res[baseline]
    for k, v in res.items():
        lo, hi = boot_ci_median(v)
        sp, slo, shi = boot_ratio(base, v)
        rows.append({"Method": k, "Median (ns/cand)": float(np.median(v)), "95% CI low": lo, "95% CI high": hi,
                     "CV (%)": float(100 * np.std(v) / np.mean(v)),
                     f"Speedup vs {baseline.split(' ')[0]}": sp, "Speedup CI low": slo, "Speedup CI high": shi,
                     "Wilcoxon p": wilcoxon_p(base, v) if k != baseline else float("nan")})
    return pd.DataFrame(rows)

# %% [markdown]
# ## E0. Environment

# %%
def cpu_name():
    try:
        with open("/proc/cpuinfo") as f:
            for line in f:
                if "model name" in line:
                    return line.split(":", 1)[1].strip()
    except Exception:
        pass
    return platform.processor() or platform.machine()

env = {
    "Date": datetime.datetime.now().isoformat(timespec="seconds"),
    "OS": platform.platform(), "CPU": cpu_name(), "Logical cores": os.cpu_count(),
    "Python": sys.version.split()[0], "NumPy": np.__version__, "pandas": pd.__version__,
    "matplotlib": matplotlib.__version__,
    "SciPy": (sstats and __import__("scipy").__version__) if HAVE_SCIPY else "n/a",
    "SymPy": sympy.__version__ if HAVE_SYMPY else "n/a",
    "gmpy2/GMP": BACKENDS[-1].version if HAVE_GMPY2 else "not installed",
    "MR backend": MRB.name, "MR rounds t": MR_ROUNDS, "Mode": MODE, "Seed": SEED,
}
SUMMARY["environment"] = env
save_table(pd.DataFrame({"Item": list(env), "Value": [str(v) for v in env.values()]}),
           "T00_environment", "Experimental environment and configuration")

# %% [markdown]
# ## E1. Correctness
# (a) All prefilters make identical decisions. (b) For `n <= B`, small primes would be rejected incorrectly,
# which is why Lemma 1 requires `n > p_k`. (c) The Miller-Rabin implementation is compared with a reference test.

# %%
log("E1 correctness")
rows = []
for bits in sorted(set(CFG["BITS_PREF"] + [512])):
    C = candidates(bits, CFG["N_PREF"], SEED + bits)
    for B in [1000, 10000]:
        for be in BACKENDS:
            ms = method_suite(B, be)
            Cb = [be.mpz(n) for n in C]
            dec = {k: [f(n) for n in Cb] for k, f in ms.items()}
            ref = dec["TD (trial division)"]
            same = all(v == ref for v in dec.values())
            rows.append({"Bits": bits, "B": B, "Backend": be.name, "Candidates": len(C),
                         "Passed": sum(ref), "Pass ratio": sum(ref) / len(C),
                         "Mertens prediction": survival_odd(primes(B)),
                         "Methods compared": len(ms), "All decisions identical": "yes" if same else "NO"})
            assert same, "Prefilter decisions differ!"
df_corr = pd.DataFrame(rows)
save_table(df_corr, "T01_prefilter_equivalence", "Prefilter equivalence: all methods give identical decisions")

# (b) precondition n > B
P = primes(1000); f = make_segmented(P, sizes_geom(P, 5), BACKENDS[0])
small_rejected = [p for p in P if not f(p)]
pre = pd.DataFrame([{"Test": "Odd primes p <= B=1000 fed to prefilter", "Count": len(P),
                     "Wrongly rejected": len(small_rejected),
                     "Conclusion": "Theorem 1 requires n > p_k (always true for cryptographic sizes)"}])
save_table(pre, "T02_precondition", "Precondition check: candidates must exceed the largest sieving prime")

# (c) Miller-Rabin vs reference
mr = make_mr(MRB)
known_composites = [561, 1105, 1729, 2465, 2821, 6601, 8911, 25326001, 3215031751, 2152302898747,
                    3474749660383, 341550071728321, 3825123056546413051,
                    318665857834031151167461]   # Carmichael numbers and strong pseudoprimes to several bases
known_primes = [2**89 - 1, 2**107 - 1, 2**127 - 1, 2**521 - 1, 2**607 - 1, 2**1279 - 1]
mr_rows = [{"Set": "Carmichael / strong pseudoprimes", "Count": len(known_composites),
            "MR says prime": sum(mr(n)[0] for n in known_composites), "Expected": 0},
           {"Set": "Mersenne primes (89..1279 bits)", "Count": len(known_primes),
            "MR says prime": sum(mr(n)[0] for n in known_primes), "Expected": len(known_primes)}]
if HAVE_SYMPY or HAVE_GMPY2:
    ref_isprime = (lambda n: bool(gmpy2.is_prime(n, 50))) if HAVE_GMPY2 else sympy.isprime
    C = candidates(512, min(CFG["N_PREF"], 3000), SEED + 99)
    passed = [n for n in C if make_trial(primes(1000), BACKENDS[0])(n)]
    agree = sum(mr(n)[0] == ref_isprime(n) for n in passed)
    mr_rows.append({"Set": "512-bit candidates passing prefilter vs reference",
                    "Count": len(passed), "MR says prime": sum(mr(n)[0] for n in passed),
                    "Expected": sum(ref_isprime(n) for n in passed)})
    mr_rows[-1]["Agreement"] = f"{agree}/{len(passed)}"
save_table(pd.DataFrame(mr_rows), "T03_mr_correctness", "Miller-Rabin implementation correctness")

# %% [markdown]
# ## E2. Survival through the segments: measured versus predicted (Mertens products)

# %%
log("E2 block survival")
bits = 1024; B = 1000; P = primes(B)
C = candidates(bits, max(CFG["N_PREF"], 2000), SEED + 2)
surv_rows, configs = [], {"geom S0=4 (paper)": sizes_geom(P, 4), "geom S0=5": sizes_geom(P, 5),
                           "fixed16": sizes_fixed(P, 16), "bits256": sizes_bits(P, 256)}
for cname, sizes in configs.items():
    blocks = make_blocks(P, sizes)
    ints = [math.prod(b) for b in blocks]
    rej = np.array([block_of_rejection(n, ints) for n in C])
    before = []
    for i, b in enumerate(blocks):
        reach = int(np.sum((rej == -1) | (rej >= i)))
        rejected = int(np.sum(rej == i))
        prior = [p for bb in blocks[:i] for p in bb]
        pred_reach = survival_odd(prior) * len(C)
        pred_rej = (survival_odd(prior) - survival_odd(prior + b)) * len(C)
        surv_rows.append({"Config": cname, "Block": i + 1, "Primes": len(b), "First prime": b[0], "Last prime": b[-1],
                          "Bits": ints[i].bit_length(), "Reaching (measured)": reach, "Reaching (Mertens)": pred_reach,
                          "Rejected (measured)": rejected, "Rejected (Mertens)": pred_rej})
df_surv = pd.DataFrame(surv_rows)
save_table(df_surv, "T04_block_survival", f"Candidates reaching/rejected by each block, measured vs Mertens prediction ({bits}-bit, B={B}, N={len(C)})")
# max relative deviation of reaching counts
dev = np.max(np.abs(df_surv["Reaching (measured)"] - df_surv["Reaching (Mertens)"]) / df_surv["Reaching (Mertens)"])
SUMMARY["E2_max_rel_dev_reaching"] = float(dev)

fig, ax = plt.subplots(1, 2, figsize=(11, 4))
for cname in configs:
    d = df_surv[df_surv.Config == cname]
    ax[0].plot(d.Block, d["Reaching (measured)"] / len(C) * 100, "o-", label=f"{cname} measured")
    ax[0].plot(d.Block, d["Reaching (Mertens)"] / len(C) * 100, "k:", lw=1)
ax[0].set_xlabel("Block index"); ax[0].set_ylabel("Candidates reaching block (%)")
ax[0].set_title("Survival through blocks (dotted: Mertens)"); ax[0].legend(fontsize=7)
d = df_surv[df_surv.Config == "geom S0=5"]
x = np.arange(len(d))
ax[1].bar(x - 0.2, d["Rejected (measured)"], 0.4, label="measured")
ax[1].bar(x + 0.2, d["Rejected (Mertens)"], 0.4, label="Mertens")
ax[1].set_xticks(x); ax[1].set_xticklabels([f"{r.Block}\n({r.Primes}p, {r.Bits}b)" for r in d.itertuples()], fontsize=8)
ax[1].set_ylabel("Rejected candidates"); ax[1].set_title("Geometric S0=5: rejections per block"); ax[1].legend()
save_fig(fig, "F03_block_survival", "Measured block survival agrees with the Mertens-product prediction")

# %% [markdown]
# ## E3. Prefilter strategies: statistical comparison (CPython integers and gmpy2)

# %%
log("E3 strategies")
B = 1000
strat_rows, strat_raw = [], {}
for bits in CFG["BITS_PREF"]:
    C = candidates(bits, CFG["N_PREF"], SEED + 3 + bits)
    for be in BACKENDS:
        ms = {"NULL (loop overhead)": make_null(primes(B), be)}
        ms.update(method_suite(B, be))
        Cb = [be.mpz(n) for n in C]
        res = bench_interleaved(ms, Cb, CFG["R"])
        df = summarize_times(res, "TD (trial division)")
        df.insert(0, "Backend", be.name); df.insert(0, "Bits", bits)
        strat_rows.append(df); strat_raw[(bits, be.name)] = res
        log(f"  {bits}-bit {be.name} done")
df_strat = pd.concat(strat_rows, ignore_index=True)
save_table(df_strat, "T05_strategy_comparison",
           f"Prefilter strategies (B={B}, odd primes 3..997, {CFG['R']} interleaved repetitions, median and bootstrap 95% CI)")

fig, axes = plt.subplots(1, len(BACKENDS), figsize=(6 * len(BACKENDS), 4.5), squeeze=False)
for j, be in enumerate(BACKENDS):
    ax = axes[0, j]
    d = df_strat[(df_strat.Backend == be.name) & (df_strat.Bits == 1024 if 1024 in CFG["BITS_PREF"] else df_strat.Bits == CFG["BITS_PREF"][0])]
    d = d[d.Method != "NULL (loop overhead)"].sort_values("Median (ns/cand)")
    y = np.arange(len(d))
    ax.barh(y, d["Median (ns/cand)"], xerr=[d["Median (ns/cand)"] - d["95% CI low"], d["95% CI high"] - d["Median (ns/cand)"]],
            color=["tab:green" if "geom" in m else "tab:gray" for m in d.Method], capsize=3)
    ax.set_yticks(y); ax.set_yticklabels(d.Method, fontsize=8); ax.invert_yaxis()
    ax.set_xlabel("Median time per candidate (ns)"); ax.set_title(f"Backend: {be.name}")
save_fig(fig, "F04_strategy_comparison", "Prefilter cost per candidate with 95% bootstrap confidence intervals")

# summary
for be in BACKENDS:
    d = df_strat[(df_strat.Backend == be.name) & (df_strat.Bits == CFG["BITS_PREF"][0])].set_index("Method")
    SUMMARY[f"E3_{be.name}_speedup_geom5_vs_TD"] = float(d.loc["SEG-geom(S0=5,r=2)", "Speedup vs TD"])
    SUMMARY[f"E3_{be.name}_best_method"] = d.drop("NULL (loop overhead)")["Median (ns/cand)"].idxmin()

# %% [markdown]
# ## E5. Cost model: GCD cost as a function of the product length (run before E4, which uses the model)

# %%
log("E5 cost model")
bits = 1024
C = candidates(bits, min(CFG["N_PREF"], 2000), SEED + 5)
P_big = primes(100000)
cost_rows, COST = [], {}
for be in BACKENDS:
    Cb = [be.mpz(n) for n in C]
    meas = []
    for k in [1, 2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64, 96, 128, 192, 256, 384, 512, 768, 1024]:
        M = be.mpz(math.prod(P_big[:k])); g = be.gcd
        f = (lambda M_: (lambda n: g(n, M_) != 1))(M)
        ts = [time_once(f, Cb) for _ in range(max(3, CFG["R"] // 3))]
        meas.append((int(M).bit_length(), float(np.median(ts))))
    mod_f = lambda n: n % 997 == 0
    c_mod = float(np.median([time_once(mod_f, Cb) for _ in range(max(3, CFG["R"] // 3))]))
    xb = np.array([m[0] for m in meas], float); yc = np.array([m[1] for m in meas])
    beta, alpha = np.polyfit(xb, yc, 1)
    COST[be.name] = dict(bits=xb, cost=yc, alpha=alpha, beta=beta, c_mod=c_mod)
    for xbits, c in meas:
        cost_rows.append({"Backend": be.name, "Primorial bits": xbits, "gcd cost (ns)": c,
                          "Linear fit (ns)": alpha + beta * xbits})
    SUMMARY[f"E5_{be.name}_gcd_alpha_ns"] = float(alpha); SUMMARY[f"E5_{be.name}_gcd_beta_ns_per_bit"] = float(beta)
    SUMMARY[f"E5_{be.name}_mod_small_ns"] = c_mod
save_table(pd.DataFrame(cost_rows), "T06_gcd_cost", f"Cost of gcd(n, M) for {bits}-bit n versus primorial bit length")

def gcd_cost(be_name, b):
    c = COST[be_name]
    if b <= c["bits"][-1]:
        return float(np.interp(b, c["bits"], c["cost"]))
    return float(c["alpha"] + c["beta"] * b)

def model_cost_segmented(P, sizes, be_name):
    """E[cost] = sum_i c(bits(M_i)) * Pr[reach block i]."""
    tot, reach = 0.0, 1.0
    for b in make_blocks(P, sizes):
        tot += reach * gcd_cost(be_name, math.prod(b).bit_length())
        reach *= survival_odd(b)
    return tot

def model_cost_trial(P, be_name):
    tot, reach = 0.0, 1.0
    for p in P:
        tot += reach; reach *= (1 - 1 / p)
    return tot * COST[be_name]["c_mod"]

fig, ax = plt.subplots(figsize=(6, 4))
for be in BACKENDS:
    c = COST[be.name]
    ax.plot(c["bits"], c["cost"], "o", label=f"{be.name} measured")
    xs = np.linspace(0, c["bits"].max(), 100); ax.plot(xs, c["alpha"] + c["beta"] * xs, "--", label=f"{be.name} fit: {c['alpha']:.0f} + {c['beta']:.3f}·bits")
ax.axvline(bits, color="k", lw=0.8, ls=":"); ax.text(bits, ax.get_ylim()[1] * 0.9, " n bits", fontsize=8)
ax.set_xlabel("Primorial segment bit length"); ax.set_ylabel("ns per gcd"); ax.legend(fontsize=8)
ax.set_title(f"gcd(n, M) cost, {bits}-bit n")
save_fig(fig, "F07_gcd_cost_model", "Cost of a single gcd versus segment size (basis of the cost model)")

# %% [markdown]
# ## E4. Initial segment size S0 and growth ratio r: measured versus model

# %%
log("E4 S0 x ratio grid")
bits, B = 1024, 1000
P = primes(B)
C = candidates(bits, min(CFG["N_PREF"], 2000), SEED + 4)
grid_rows = []
for be in BACKENDS:
    Cb = [be.mpz(n) for n in C]
    ms = {}
    for r in CFG["RATIO_GRID"]:
        for s0 in CFG["S0_GRID"]:
            ms[(s0, r)] = make_segmented(P, sizes_geom(P, s0, r), be)
    ms_named = {f"{k[0]}|{k[1]}": v for k, v in ms.items()}
    res = bench_interleaved(ms_named, Cb, CFG["R_GRID"])
    for (s0, r), fobj in ms.items():
        v = res[f"{s0}|{r}"]; lo, hi = boot_ci_median(v)
        grid_rows.append({"Backend": be.name, "S0": s0, "Ratio r": r, "Blocks": len(fobj.mods),
                          "Median (ns/cand)": float(np.median(v)), "95% CI low": lo, "95% CI high": hi,
                          "Model (ns/cand)": model_cost_segmented(P, sizes_geom(P, s0, r), be.name)})
    log(f"  grid {be.name} done")
df_grid = pd.DataFrame(grid_rows)
save_table(df_grid, "T07_S0_ratio_grid", f"Geometric segmentation: initial size S0 and growth ratio r ({bits}-bit, B={B}, {CFG['R_GRID']} repetitions)")

opt_rows = []
for be in BACKENDS:
    d = df_grid[df_grid.Backend == be.name]
    best = d.loc[d["Median (ns/cand)"].idxmin()]
    # configurations statistically indistinguishable from the best (CI overlap)
    tie = d[d["95% CI low"] <= best["95% CI high"]]
    mbest = d.loc[d["Model (ns/cand)"].idxmin()]
    rho = float(sstats.spearmanr(d["Median (ns/cand)"], d["Model (ns/cand)"]).correlation) if HAVE_SCIPY else float("nan")
    mape = float(np.mean(np.abs(d["Model (ns/cand)"] - d["Median (ns/cand)"]) / d["Median (ns/cand)"]) * 100)
    d2 = d[d["Ratio r"] == 2.0]
    best2 = d2.loc[d2["Median (ns/cand)"].idxmin()] if len(d2) else best
    opt_rows.append({"Backend": be.name, "Best S0": int(best.S0), "Best r": best["Ratio r"],
                     "Best median (ns)": best["Median (ns/cand)"],
                     "Best S0 at r=2": int(best2.S0), "Configs within CI of best": len(tie),
                     "S0 range within CI": f"{int(tie.S0.min())}-{int(tie.S0.max())}",
                     "Model-optimal S0": int(mbest.S0), "Model-optimal r": mbest["Ratio r"],
                     "Spearman(model, measured)": rho, "Model MAPE (%)": mape})
df_opt = pd.DataFrame(opt_rows)
save_table(df_opt, "T08_S0_optimum", "Optimal geometric configuration and statistically equivalent configurations")
SUMMARY["E4_optimum"] = df_opt.to_dict(orient="records")

fig, axes = plt.subplots(1, len(BACKENDS), figsize=(6 * len(BACKENDS), 4), squeeze=False)
for j, be in enumerate(BACKENDS):
    ax = axes[0, j]
    for r in CFG["RATIO_GRID"]:
        d = df_grid[(df_grid.Backend == be.name) & (df_grid["Ratio r"] == r)].sort_values("S0")
        l, = ax.plot(d.S0, d["Median (ns/cand)"], "o-", ms=3, label=f"r={r} measured")
        ax.fill_between(d.S0, d["95% CI low"], d["95% CI high"], color=l.get_color(), alpha=0.15)
        ax.plot(d.S0, d["Model (ns/cand)"], ":", color=l.get_color())
    ax.set_xlabel("Initial segment size S0 (primes)"); ax.set_ylabel("ns per candidate")
    ax.set_title(f"{be.name}: measured (solid, 95% CI) vs model (dotted)"); ax.legend(fontsize=7)
save_fig(fig, "F05_S0_sensitivity", "Sensitivity to the initial segment size S0 and growth ratio r")

if len(CFG["RATIO_GRID"]) > 1 and len(CFG["S0_GRID"]) > 1:
    fig, axes = plt.subplots(1, len(BACKENDS), figsize=(6 * len(BACKENDS), 3.2), squeeze=False)
    for j, be in enumerate(BACKENDS):
        d = df_grid[df_grid.Backend == be.name].pivot(index="Ratio r", columns="S0", values="Median (ns/cand)")
        im = axes[0, j].imshow(d.values, aspect="auto", cmap="viridis_r", origin="lower")
        axes[0, j].set_yticks(range(len(d.index))); axes[0, j].set_yticklabels(d.index)
        axes[0, j].set_xticks(range(0, len(d.columns), max(1, len(d.columns) // 8)))
        axes[0, j].set_xticklabels(d.columns[::max(1, len(d.columns) // 8)])
        axes[0, j].set_xlabel("S0"); axes[0, j].set_ylabel("ratio r"); axes[0, j].grid(False)
        axes[0, j].set_title(f"{be.name}: ns/candidate"); fig.colorbar(im, ax=axes[0, j])
    save_fig(fig, "F06_S0_ratio_heatmap", "Prefilter cost over the (S0, r) grid")

fig, ax = plt.subplots(figsize=(5, 4))
for be in BACKENDS:
    d = df_grid[df_grid.Backend == be.name]
    ax.scatter(d["Model (ns/cand)"], d["Median (ns/cand)"], s=10, label=be.name)
lim = [0, df_grid[["Model (ns/cand)", "Median (ns/cand)"]].values.max() * 1.05]
ax.plot(lim, lim, "k--", lw=0.8); ax.set_xlabel("Model prediction (ns)"); ax.set_ylabel("Measured (ns)")
ax.set_title("Cost model validation"); ax.legend()
save_fig(fig, "F08_model_vs_measured", "Cost-model predictions versus measured prefilter cost over all (S0, r) configurations")

# %% [markdown]
# ## E6. Optimization of the sieving bound B
# For every B: measured prefilter cost, pass fraction, and cost of one MR round.
# Expected cost per candidate: `C(B) = c_pre(B) + π(B)·c_round`.

# %%
log("E6 sieve bound optimisation")
def mr_round_cost(bits, k=None):
    """Median time of one MR round (one modexp + squarings) on random odd composites."""
    k = k or (40 if bits <= 1024 else 15 if bits <= 2048 else 6)
    C = candidates(bits, k, SEED + 600 + bits)
    mr1 = make_mr(MRB, t=1)
    ts = []
    for n in C:
        t0 = perf_counter_ns(); mr1(MRB.mpz(n)); ts.append(perf_counter_ns() - t0)
    return float(np.median(ts))

C_ROUND = {}
for b in sorted(set(CFG["BITS_B"] + CFG["BITS_PIPE"])):
    C_ROUND[b] = mr_round_cost(b)
    log(f"  MR round cost {b}-bit: {C_ROUND[b]/1e6:.3f} ms ({MRB.name})")
SUMMARY["mr_round_cost_ns"] = C_ROUND

brows = []
for bits in CFG["BITS_B"]:
    C = candidates(bits, CFG["N_B"], SEED + 6 + bits)
    for be in BACKENDS:
        Cb = [be.mpz(n) for n in C]
        for B in CFG["B_GRID"]:
            P = primes(B)
            ms = {"GEOM": make_segmented(P, sizes_geom(P, 5, 2.0), be)}
            if B <= CFG["TD_MAX_B"]:
                ms["TD"] = make_trial(P, be)
            res = bench_interleaved(ms, Cb, CFG["R_B"])
            passed = sum(ms["GEOM"](n) for n in Cb) / len(Cb)
            row = {"Bits": bits, "Backend": be.name, "B": B, "Primes": len(P), "Blocks": len(ms["GEOM"].mods),
                   "Pass ratio": passed, "Mertens pass": survival_odd(P),
                   "GEOM prefilter (us)": float(np.median(res["GEOM"])) / 1e3,
                   "TD prefilter (us)": float(np.median(res["TD"])) / 1e3 if "TD" in res else float("nan"),
                   "MR round (us)": C_ROUND[bits] / 1e3}
            row["Expected cost/cand GEOM (us)"] = row["GEOM prefilter (us)"] + passed * row["MR round (us)"]
            row["Expected cost/cand TD (us)"] = row["TD prefilter (us)"] + passed * row["MR round (us)"]
            brows.append(row)
        log(f"  B-sweep {bits}-bit {be.name} done")
df_B = pd.DataFrame(brows)
save_table(df_B, "T09_B_sweep", f"Sieve bound sweep: prefilter cost, pass ratio and expected cost per candidate (S0=5, r=2; MR backend {MRB.name})")

bopt_rows, BSTAR = [], {}
for (bits, ben), d in df_B.groupby(["Bits", "Backend"]):
    d = d.sort_values("B")
    i = d["Expected cost/cand GEOM (us)"].idxmin(); best = d.loc[i]
    ref = d[d.B == 1000].iloc[0] if (d.B == 1000).any() else d.iloc[0]
    dt = d.dropna(subset=["Expected cost/cand TD (us)"])
    bt = dt.loc[dt["Expected cost/cand TD (us)"].idxmin()] if len(dt) else None
    if ben == MRB.name:
        BSTAR[bits] = int(best.B)
    bopt_rows.append({"Bits": bits, "Backend": ben, "B* (GEOM)": int(best.B), "Pass ratio at B*": best["Pass ratio"],
                      "Cost/cand at B=1000 (us)": ref["Expected cost/cand GEOM (us)"],
                      "Cost/cand at B* (us)": best["Expected cost/cand GEOM (us)"],
                      "Gain vs B=1000 (%)": 100 * (1 - best["Expected cost/cand GEOM (us)"] / ref["Expected cost/cand GEOM (us)"]),
                      "MR calls saved vs B=1000 (%)": 100 * (1 - best["Pass ratio"] / ref["Pass ratio"]),
                      "B* (TD, tested range)": int(bt.B) if bt is not None else float("nan"),
                      "Cost/cand TD at its B* (us)": bt["Expected cost/cand TD (us)"] if bt is not None else float("nan")})
df_bopt = pd.DataFrame(bopt_rows)
save_table(df_bopt, "T10_B_optimum", "Optimal sieve bound B* for the geometric prefilter versus the conventional B=1000")
SUMMARY["E6_Bstar"] = df_bopt.to_dict(orient="records")
# pipeline bit lengths not covered by the B-sweep use the nearest tested bit length
_tested = sorted(BSTAR)
for b in CFG["BITS_PIPE"] + list(CFG["PRIMEGEN"]):
    if b not in BSTAR:
        BSTAR[b] = BSTAR[min(_tested, key=lambda x: abs(x - b))] if _tested else 1000
log(f"  B* used in pipeline: {BSTAR}")

fig, axes = plt.subplots(1, 2, figsize=(12, 4.2))
for be in BACKENDS:
    d = df_B[(df_B.Backend == be.name) & (df_B.Bits == CFG["BITS_B"][min(1, len(CFG["BITS_B"]) - 1)])].sort_values("B")
    axes[0].loglog(d.B, d["GEOM prefilter (us)"], "o-", label=f"GEOM, {be.name}")
    dd = d.dropna(subset=["TD prefilter (us)"])
    axes[0].loglog(dd.B, dd["TD prefilter (us)"], "s--", label=f"TD, {be.name}")
axes[0].set_xlabel("Sieve bound B"); axes[0].set_ylabel("Prefilter time per candidate (us)")
axes[0].set_title(f"Prefilter cost vs B ({CFG['BITS_B'][min(1, len(CFG['BITS_B']) - 1)]}-bit)"); axes[0].legend(fontsize=8)
for bits in CFG["BITS_B"]:
    d = df_B[(df_B.Backend == MRB.name) & (df_B.Bits == bits)].sort_values("B")
    if not len(d):
        continue
    rel = d["Expected cost/cand GEOM (us)"] / d[d.B == d.B.min()]["Expected cost/cand GEOM (us)"].iloc[0] if not (d.B == 1000).any() \
        else d["Expected cost/cand GEOM (us)"] / d[d.B == 1000]["Expected cost/cand GEOM (us)"].iloc[0]
    l, = axes[1].semilogx(d.B, rel, "o-", label=f"{bits}-bit")
    j = rel.values.argmin(); axes[1].plot(d.B.values[j], rel.values[j], "*", ms=14, color=l.get_color())
axes[1].axhline(1, color="k", lw=0.8, ls=":")
axes[1].set_xlabel("Sieve bound B"); axes[1].set_ylabel("Expected cost per candidate (relative to B=1000)")
axes[1].set_title(f"Total expected cost vs B (star = B*), backend {MRB.name}"); axes[1].legend(fontsize=8)
save_fig(fig, "F09_B_optimisation", "Effect of the sieve bound B: prefilter cost (left) and total expected cost per candidate (right)")

# %% [markdown]
# ## E7. Complete primality-testing pipeline, 512-4096 bits
# The MR time of every candidate is measured **once** and stored (the number of MR rounds is deterministic),
# so that differences between methods are due to the prefilter alone and MR timing noise is removed.
# The number of MR rounds is reported as a hardware-independent measure. Because the sample contains few primes
# (each adds t rounds), the "expected" column normalizes the speedup with the theoretical prime density.

# %%
log("E7 full pipeline")
def theoretical_speedup(bits, B, t):
    rho = prime_density_odd(bits); pi = survival_odd(primes(B))
    return ((1 - rho) + rho * t) / ((pi - rho) + rho * t)

pipe_rows = []
for bits in CFG["BITS_PIPE"]:
    N = CFG["N_PIPE"][bits]
    C = candidates(bits, N, SEED + 7 + bits)
    Cm = [MRB.mpz(n) for n in C]
    mrt = make_mr(MRB)
    mr_time, mr_rounds, is_p = np.zeros(N), np.zeros(N, int), np.zeros(N, bool)
    for i, n in enumerate(Cm):
        t0 = perf_counter_ns(); ok, rds = mrt(n); mr_time[i] = perf_counter_ns() - t0
        mr_rounds[i], is_p[i] = rds, ok
    log(f"  {bits}-bit: MR on {N} candidates done, primes={is_p.sum()}")
    be = MRB
    Bs = sorted({1000, BSTAR.get(bits, 1000)})
    methods = {"MR only": (None, None)}
    methods["TD(B=1000) + MR"] = (make_trial(primes(1000), be), 1000)
    methods["SP(B=1000) + MR"] = (make_single(primes(1000), be), 1000)
    methods["GEOM(S0=5,B=1000) + MR"] = (make_segmented(primes(1000), sizes_geom(primes(1000), 5), be), 1000)
    if BSTAR.get(bits, 1000) != 1000:
        Bx = BSTAR[bits]
        methods[f"GEOM(S0=5,B*={Bx}) + MR"] = (make_segmented(primes(Bx), sizes_geom(primes(Bx), 5), be), Bx)
    pref = {k: v[0] for k, v in methods.items() if v[0] is not None}
    Csub = Cm[:min(N, 3000)]
    tres = bench_interleaved(pref, Csub, max(5, CFG["R"] // 2))
    base_rows = {}
    for k, (f, B) in methods.items():
        if f is None:
            mask = np.ones(N, bool); pre_s = 0.0
        else:
            mask = np.array([f(n) for n in Cm]); pre_s = float(np.median(tres[k])) * N / 1e9
        mr_s = float(mr_time[mask].sum()) / 1e9
        row = {"Bits": bits, "Method": k, "N": N, "Passed to MR": int(mask.sum()), "Primes found": int(is_p[mask].sum()),
               "MR rounds": int(mr_rounds[mask].sum()), "Prefilter (s)": pre_s, "MR (s)": mr_s, "Total (s)": pre_s + mr_s,
               "Theoretical bound": theoretical_speedup(bits, B, MR_ROUNDS) if B else 1.0}
        base_rows[k] = row
    # prime-density-normalised expectation (removes sampling noise from the few primes in the sample)
    rho = prime_density_odd(bits)
    c_comp = float(mr_time[~is_p].mean()) / 1e9                      # MR time of a composite candidate (s)
    c_rnd = float((mr_time[~is_p] / np.maximum(mr_rounds[~is_p], 1)).mean()) / 1e9
    e_none = (1 - rho) * c_comp + rho * MR_ROUNDS * c_rnd
    for k, (f, B) in methods.items():
        row = base_rows[k]
        pi_m = row["Passed to MR"] / N
        e = (row["Prefilter (s)"] / N) + max(pi_m - rho, 0) * c_comp + rho * MR_ROUNDS * c_rnd
        row["Expected time/cand (ms)"] = e * 1e3
        row["Speedup vs MR only (expected)"] = e_none / e
    t_none = base_rows["MR only"]["Total (s)"]; r_none = base_rows["MR only"]["MR rounds"]
    t_td = base_rows["TD(B=1000) + MR"]["Total (s)"]
    for k, row in base_rows.items():
        row["Speedup vs MR only (time)"] = t_none / row["Total (s)"]
        row["Speedup vs MR only (rounds)"] = r_none / row["MR rounds"]
        row["Speedup vs TD+MR (time)"] = t_td / row["Total (s)"]
        pipe_rows.append(row)
df_pipe = pd.DataFrame(pipe_rows)
save_table(df_pipe, "T11_pipeline", f"Complete primality-testing pipeline (t={MR_ROUNDS} MR rounds, MR backend {MRB.name}; MR cost measured once per candidate)")
SUMMARY["E7_pipeline"] = df_pipe[["Bits", "Method", "Speedup vs MR only (expected)", "Speedup vs MR only (time)", "Speedup vs MR only (rounds)",
                                  "Speedup vs TD+MR (time)", "Theoretical bound"]].to_dict(orient="records")

# theoretical bound table for t sensitivity
tb = []
for bits in [256, 512, 1024, 2048, 3072, 4096]:
    for B in [1000, 10000, 100000]:
        r = {"Bits": bits, "B": B, "Pass ratio (Mertens)": survival_odd(primes(B)), "Prime density (odd)": prime_density_odd(bits)}
        for t in MR_ROUNDS_SENS:
            r[f"t={t}"] = theoretical_speedup(bits, B, t)
        tb.append(r)
df_tb = pd.DataFrame(tb)
save_table(df_tb, "T12_theoretical_bound", "Theoretical maximum speedup of prefilter+MR over MR alone (in MR rounds)")

fig, ax = plt.subplots(figsize=(7.5, 4.5))
bl = sorted(df_pipe.Bits.unique())
for k, mk in [("TD(B=1000) + MR", "s"), ("GEOM(S0=5,B=1000) + MR", "o")]:
    d = df_pipe[df_pipe.Method == k].set_index("Bits").reindex(bl)
    l, = ax.plot(bl, d["Speedup vs MR only (expected)"], mk + "-", label=f"{k} (expected)")
    ax.plot(bl, d["Speedup vs MR only (time)"], mk, mfc="none", color=l.get_color(), alpha=0.6, label=f"{k} (this sample)")
dstar = df_pipe[df_pipe.Method.str.contains("B\\*")]
if len(dstar):
    ax.plot(dstar.Bits, dstar["Speedup vs MR only (expected)"], "D-", label="GEOM(S0=5, B*) + MR (expected)")
xs = np.array(bl)
ax.plot(xs, [theoretical_speedup(b, 1000, MR_ROUNDS) for b in xs], "k--", lw=1, label=f"Theory B=1000, t={MR_ROUNDS}")
ax.plot(xs, [theoretical_speedup(b, 1000, 1) for b in xs], color="gray", ls="-.", lw=1, label="Theory B=1000, t=1 (upper bound)")
if len(dstar):
    ax.plot(xs, [theoretical_speedup(b, BSTAR.get(b, 1000), MR_ROUNDS) for b in xs], color="tab:green", ls="--", lw=0.8, label=f"Theory B*, t={MR_ROUNDS}")
ax.set_xscale("log", base=2); ax.set_xticks(xs); ax.set_xticklabels(xs)
ax.set_xlabel("Candidate bit length"); ax.set_ylabel("Speedup over Miller–Rabin alone")
ax.set_title("Pipeline speedup: measured vs theoretical bound"); ax.legend(fontsize=7)
save_fig(fig, "F11_pipeline_speedup", "Measured pipeline speedup with theoretical bounds")

fig, ax = plt.subplots(figsize=(7.5, 4.2))
for k in df_pipe.Method.unique():
    d = df_pipe[df_pipe.Method == k]
    ax.plot(d.Bits, d["Total (s)"] / d.N * 1e3, "o-", label=k)
ax.set_xscale("log", base=2); ax.set_yscale("log"); ax.set_xticks(bl); ax.set_xticklabels(bl)
ax.set_xlabel("Candidate bit length"); ax.set_ylabel("Average time per candidate (ms)")
ax.set_title("Average pipeline time per candidate"); ax.legend(fontsize=7)
save_fig(fig, "F12_pipeline_time", "Average time per candidate for the complete pipeline")

# %% [markdown]
# ## E8. End-to-end timing noise compared with prefilter differences

# %%
log("E8 end-to-end noise")
bits = 1024
C = [MRB.mpz(n) for n in candidates(bits, CFG["N_E2E"], SEED + 8)]
mrt = make_mr(MRB); P = primes(1000)
pf = {"TD": make_trial(P, MRB), "SP": make_single(P, MRB), "GEOM(S0=5)": make_segmented(P, sizes_geom(P, 5), MRB)}
tot = {k: [] for k in pf}; prt = {k: [] for k in pf}
ordr = random.Random(SEED + 80)
for rep in range(CFG["R_E2E"]):
    names = list(pf); ordr.shuffle(names)
    for k in names:
        f = pf[k]
        t0 = perf_counter_ns(); passed = [n for n in C if f(n)]; t1 = perf_counter_ns()
        for n in passed:
            mrt(n)
        t2 = perf_counter_ns()
        tot[k].append((t2 - t0) / 1e9); prt[k].append((t1 - t0) / 1e9)
e2e = []
for k in pf:
    e2e.append({"Method": k, "Runs": len(tot[k]), "Total mean (s)": np.mean(tot[k]), "Total std (s)": np.std(tot[k], ddof=1) if len(tot[k]) > 1 else float("nan"),
                "Prefilter mean (s)": np.mean(prt[k]), "MR part mean (s)": np.mean(np.array(tot[k]) - np.array(prt[k])),
                "MR part std (s)": np.std(np.array(tot[k]) - np.array(prt[k]), ddof=1) if len(tot[k]) > 1 else float("nan")})
df_e2e = pd.DataFrame(e2e)
pre_diff = df_e2e["Prefilter mean (s)"].max() - df_e2e["Prefilter mean (s)"].min()
df_e2e["Prefilter difference (max-min, s)"] = pre_diff
df_e2e["Noise / prefilter difference"] = df_e2e["MR part std (s)"] / pre_diff if pre_diff > 0 else float("nan")
save_table(df_e2e, "T14_e2e_noise", f"End-to-end runs ({bits}-bit, N={len(C)}, {CFG['R_E2E']} interleaved runs): MR timing noise versus prefilter differences")
SUMMARY["E8_noise_ratio_mean"] = float(np.nanmean(df_e2e["Noise / prefilter difference"]))

def _bp(a, data, labels):
    try:
        a.boxplot(data, tick_labels=labels)
    except TypeError:
        a.boxplot(data, labels=labels)
fig, ax = plt.subplots(1, 2, figsize=(10, 3.8))
_bp(ax[0], [tot[k] for k in pf], list(pf)); ax[0].set_ylabel("Total time (s)"); ax[0].set_title("End-to-end total (MR dominates)")
_bp(ax[1], [prt[k] for k in pf], list(pf)); ax[1].set_ylabel("Prefilter time (s)"); ax[1].set_title("Prefilter part only")
save_fig(fig, "F13_e2e_noise", "End-to-end timing: variability of the Miller–Rabin part exceeds prefilter differences")

# %% [markdown]
# ## E9. Prime generation
# (A) independent random candidates + TD(B=1000), (B) independent + GEOM(B*),
# (C) incremental search n0, n0+2, ... + GEOM(B*), (D) incremental sieve (OpenSSL style, B*).
# (C) and (D) must find the same prime from the same starting point (correctness check).

# %%
log("E9 prime generation")
def gen_independent(bits, seed, f, mr):
    r = random.Random(seed); tested = rounds = 0
    while True:
        n = MRB.mpz(r.getrandbits(bits) | (1 << (bits - 1)) | 1); tested += 1
        if f(n):
            ok, rd = mr(n); rounds += rd
            if ok:
                return int(n), tested, rounds

def gen_incremental(bits, seed, f, mr):
    r = random.Random(seed); n = MRB.mpz(r.getrandbits(bits) | (1 << (bits - 1)) | 1); tested = rounds = 0
    while True:
        tested += 1
        if f(n):
            ok, rd = mr(n); rounds += rd
            if ok:
                return int(n), tested, rounds
        n += 2

def gen_sieve(bits, seed, P, mr, W=None):
    W = W or max(1024, 2 * bits)
    r = random.Random(seed); n0 = MRB.mpz(r.getrandbits(bits) | (1 << (bits - 1)) | 1); tested = rounds = 0
    while True:
        mark = bytearray(W)
        for p in P:
            j = ((p - int(n0 % p)) * ((p + 1) // 2)) % p
            if j < W:
                mark[j::p] = b"\x01" * ((W - 1 - j) // p + 1)
        for j in range(W):
            tested += 1
            if not mark[j]:
                n = n0 + 2 * j
                ok, rd = mr(n); rounds += rd
                if ok:
                    return int(n), tested, rounds
        n0 += 2 * W

pg_rows, pg_raw = [], {}
for bits, K in CFG["PRIMEGEN"].items():
    Bx = BSTAR.get(bits, 1000)
    mr = make_mr(MRB); P1, Px = primes(1000), primes(Bx)
    fTD = make_trial(P1, MRB); fG = make_segmented(Px, sizes_geom(Px, 5), MRB)
    strategies = {"A: random + TD(B=1000)": lambda s: gen_independent(bits, s, fTD, mr),
                  f"B: random + GEOM(B*={Bx})": lambda s: gen_independent(bits, s, fG, mr),
                  f"C: incremental + GEOM(B*={Bx})": lambda s: gen_incremental(bits, s, fG, mr),
                  f"D: incremental sieve (B*={Bx})": lambda s: gen_sieve(bits, s, Px, mr)}
    found = {}
    for name, g in strategies.items():
        ts, tested, rounds, ps = [], [], [], []
        for k in range(K):
            t0 = perf_counter_ns(); p, nt, nr = g(SEED * 1000 + bits * 10 + k); ts.append((perf_counter_ns() - t0) / 1e6)
            tested.append(nt); rounds.append(nr); ps.append(p)
        found[name] = ps; pg_raw[(bits, name)] = ts
        lo, hi = boot_ci_median(ts)
        pg_rows.append({"Bits": bits, "Strategy": name, "Primes": K, "Median time/prime (ms)": float(np.median(ts)),
                        "95% CI low": lo, "95% CI high": hi, "Mean time/prime (ms)": float(np.mean(ts)),
                        "Candidates/prime": float(np.mean(tested)), "MR rounds/prime": float(np.mean(rounds))})
        log(f"  {bits}-bit {name}: {np.median(ts):.1f} ms/prime")
    keys = list(found)
    same_AB = found[keys[0]] == found[keys[1]]; same_CD = found[keys[2]] == found[keys[3]]
    for r in pg_rows[-4:]:
        r["Same primes as paired strategy"] = "yes" if (same_AB if r["Strategy"][0] in "AB" else same_CD) else "NO"
df_pg = pd.DataFrame(pg_rows)
save_table(df_pg, "T15_prime_generation", f"Prime generation: time, candidates and MR rounds per prime (t={MR_ROUNDS}, MR backend {MRB.name})")
SUMMARY["E9_primegen"] = df_pg.to_dict(orient="records")

fig, axes = plt.subplots(1, len(CFG["PRIMEGEN"]), figsize=(6 * len(CFG["PRIMEGEN"]), 4), squeeze=False)
for j, bits in enumerate(CFG["PRIMEGEN"]):
    d = df_pg[df_pg.Bits == bits]
    y = np.arange(len(d))
    axes[0, j].barh(y, d["Median time/prime (ms)"], xerr=[d["Median time/prime (ms)"] - d["95% CI low"], d["95% CI high"] - d["Median time/prime (ms)"]], capsize=3)
    axes[0, j].set_yticks(y); axes[0, j].set_yticklabels(d.Strategy, fontsize=8); axes[0, j].invert_yaxis()
    axes[0, j].set_xlabel("Median time per prime (ms)"); axes[0, j].set_title(f"{bits}-bit prime generation")
save_fig(fig, "F14_prime_generation", "Prime generation strategies (median and 95% CI)")

# %% [markdown]
# ## Diagrams: workflow and segmentation schematic

# %%
log("Diagrams")
def box(ax, x, y, w, h, text, fc, ec="k", shape="rect", fs=8.5):
    if shape == "diamond":
        ax.add_patch(Polygon([[x, y + h / 2], [x + w / 2, y], [x, y - h / 2], [x - w / 2, y]], closed=True, fc=fc, ec=ec, lw=1))
    else:
        ax.add_patch(FancyBboxPatch((x - w / 2, y - h / 2), w, h, boxstyle="round,pad=0.02,rounding_size=0.08", fc=fc, ec=ec, lw=1))
    ax.text(x, y, text, ha="center", va="center", fontsize=fs)

def arrow(ax, x1, y1, x2, y2, text=None, color="k", tx=0.08, ty=0.0):
    ax.annotate("", xy=(x2, y2), xytext=(x1, y1), arrowprops=dict(arrowstyle="-|>", color=color, lw=1.1))
    if text:
        ax.text((x1 + x2) / 2 + tx, (y1 + y2) / 2 + ty, text, fontsize=8, color=color)

fig, ax = plt.subplots(figsize=(7.5, 9.5)); ax.set_xlim(0, 10); ax.set_ylim(0, 13); ax.axis("off")
G, Y, R_, B_ = "#e8f4e8", "#fff4d6", "#fde2e2", "#e3edf9"
box(ax, 5, 12.3, 5.4, 0.8, "Random odd candidate $n$  ($n > p_k$)", B_)
box(ax, 5, 11.0, 6.4, 1.0, "Precomputed once: $P=\\{p_1,\\dots,p_k\\}$ (odd primes $\\leq B$)\n"
    "partitioned into $P^{(1)},\\dots,P^{(m)}$, $|P^{(i)}|=\\mathrm{round}(S_0 r^{i-1})$;  $M_i=\\prod_{p\\in P^{(i)}}p$", "#f0f0f0")
box(ax, 5, 9.7, 2.2, 0.6, "$i \\leftarrow 1$", G)
box(ax, 5, 8.6, 3.4, 0.7, "$g_i=\\gcd(n, M_i)$", G)
box(ax, 5, 7.2, 2.6, 1.2, "$g_i>1$ ?", Y, shape="diamond")
box(ax, 8.4, 7.2, 2.6, 0.9, "Reject: composite\n(small prime factor)", R_)
box(ax, 5, 5.6, 2.6, 1.2, "$i=m$ ?", Y, shape="diamond")
box(ax, 1.8, 5.6, 2.0, 0.6, "$i \\leftarrow i+1$", G)
box(ax, 5, 4.1, 4.0, 0.8, "Miller–Rabin, $t$ rounds", B_)
box(ax, 5, 2.6, 2.8, 1.2, "probable\nprime ?", Y, shape="diamond")
box(ax, 8.4, 2.6, 2.4, 0.8, "Reject: composite", R_)
box(ax, 5, 1.0, 3.6, 0.8, "Accept: probable prime", G)
arrow(ax, 5, 11.9, 5, 11.5); arrow(ax, 5, 10.5, 5, 10.0); arrow(ax, 5, 9.4, 5, 8.95); arrow(ax, 5, 8.25, 5, 7.8)
arrow(ax, 6.3, 7.2, 7.1, 7.2, "yes", "tab:red", -0.3, 0.12); arrow(ax, 5, 6.6, 5, 6.2, "no", "tab:green", 0.08)
arrow(ax, 3.7, 5.6, 2.8, 5.6, "no", "tab:green", -0.2, 0.12)
ax.annotate("", xy=(3.3, 8.6), xytext=(1.8, 5.9), arrowprops=dict(arrowstyle="-|>", lw=1.1, connectionstyle="angle,angleA=90,angleB=180"))
arrow(ax, 5, 5.0, 5, 4.5, "yes", "tab:green", 0.12, -0.12); arrow(ax, 5, 3.7, 5, 3.2)
arrow(ax, 6.4, 2.6, 7.2, 2.6, "no", "tab:red", -0.3, 0.12); arrow(ax, 5, 2.0, 5, 1.4, "yes", "tab:green", 0.08)
ax.text(5, 6.35, "", fontsize=7)
ax.add_patch(FancyBboxPatch((0.4, 4.85), 9.2, 5.25, boxstyle="round,pad=0.05", fc="none", ec="purple", ls="--", lw=1))
ax.text(0.6, 9.9, "Segmented primorial prefilter", color="purple", fontsize=9, style="italic")
save_fig(fig, "F01_workflow", "Workflow of the segmented primorial prefilter (notation consistent with Section 3)")

P = primes(1000); sz = sizes_geom(P, 5); blocks = make_blocks(P, sz)
d = df_surv[df_surv.Config == "geom S0=5"]
fig, ax = plt.subplots(figsize=(10, 3.0))
cols = plt.cm.Greens(np.linspace(0.25, 0.75, len(blocks)))
reach0 = d["Reaching (measured)"].iloc[0]
for i, b in enumerate(blocks):
    x = i * 1.7
    ax.add_patch(FancyBboxPatch((x, 0), 1.4, 1.2, boxstyle="round,pad=0.02", fc=cols[i], ec="k"))
    ax.text(x + 0.7, 0.6, f"$M_{{{i+1}}}$\n{len(b)} primes\n{b[0]}–{b[-1]}\n{math.prod(b).bit_length()} bits",
            ha="center", va="center", fontsize=8)
    reach = d.iloc[i]["Reaching (measured)"] / reach0 * 100
    ax.text(x + 0.7, 1.42, f"reach: {reach:.1f}%", ha="center", fontsize=8, color="navy")
    if i < len(blocks) - 1:
        ax.annotate("", xy=(x + 1.68, 0.6), xytext=(x + 1.42, 0.6), arrowprops=dict(arrowstyle="-|>", lw=1))
ax.text(len(blocks) * 1.7 - 0.2, 0.6, "→ MR", fontsize=9, va="center")
ax.set_xlim(-0.2, len(blocks) * 1.7 + 0.8); ax.set_ylim(-0.15, 1.7); ax.axis("off")
ax.set_title("Geometric segmentation (S0=5, r=2, B=1000): size, bit length and share of candidates reaching each segment", fontsize=9)
save_fig(fig, "F02_segmentation_schematic", "Geometric segmentation: segment sizes, bit lengths and fraction of candidates reaching each segment")

# %% [markdown]
# ## Final report and archive

# %%
SUMMARY["runtime_minutes"] = (time.time() - T_START) / 60
with open(f"{OUTDIR}/summary.json", "w", encoding="utf-8") as f:
    json.dump(SUMMARY, f, indent=2, default=lambda o: o.item() if hasattr(o, "item") else str(o))

lines = [f"# Segmented Primorial Prefilter — experiment report", "",
         f"Mode: **{MODE}**, seed {SEED}, MR rounds t={MR_ROUNDS}, runtime {SUMMARY['runtime_minutes']:.1f} min", ""]
for kind, name, cap, payload in REPORT:
    if kind == "table":
        lines += [f"## {name}: {cap}", "", payload, ""]
    else:
        lines += [f"## {name}: {cap}", "", f"![{name}]({payload})", ""]
lines += ["## Summary (JSON)", "", "```json", json.dumps(SUMMARY, indent=1, default=str)[:20000], "```"]
with open(f"{OUTDIR}/report.md", "w", encoding="utf-8") as f:
    f.write("\n".join(lines))

with zipfile.ZipFile("results.zip", "w", zipfile.ZIP_DEFLATED) as z:
    for root, _, files in os.walk(OUTDIR):
        for fn in files:
            z.write(os.path.join(root, fn))
log(f"DONE. Results archive: results.zip  ({os.path.getsize('results.zip')/1e6:.1f} MB)")
