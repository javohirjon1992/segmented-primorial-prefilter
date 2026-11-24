"""
dtlib — shared helper module of the Python prototype experiments.

The notebooks import this module; keep it in the same folder as the notebooks. On Windows, multiprocessing
worker functions must live in a module rather than in a notebook, so they are defined here as well.
"""
import ctypes, ctypes.util, datetime, glob, json, math, os, platform, random, statistics, sys, time, zipfile
from time import perf_counter_ns

import numpy as np

# ----------------------------------------------------------------------------- arifmetika
try:
    import gmpy2
    mpz, gcd, powmod = gmpy2.mpz, gmpy2.gcd, gmpy2.powmod
    BACKEND = f"GMP {gmpy2.mp_version()} (gmpy2 {gmpy2.version()})"
    HAVE_GMPY2 = True

    def ref_is_prime(n, reps=30):
        return bool(gmpy2.is_prime(n, reps))

    def ref_next_prime(n):
        """Smallest prime >= n."""
        return int(gmpy2.next_prime(int(n) - 1))

    def ref_prev_prime(n):
        """Largest prime < n."""
        return int(gmpy2.prev_prime(int(n))) if hasattr(gmpy2, "prev_prime") else _prev_prime_fallback(n)
except ImportError:                                  # fallback for testing (slow)
    mpz, gcd, powmod = int, math.gcd, pow
    BACKEND = f"CPython int (Python {sys.version.split()[0]}) — gmpy2 tavsiya etiladi"
    HAVE_GMPY2 = False
    import sympy

    def ref_is_prime(n, reps=30):
        return bool(sympy.isprime(int(n)))

    def ref_next_prime(n):
        return int(sympy.nextprime(int(n) - 1))

    def ref_prev_prime(n):
        return int(sympy.prevprime(int(n)))


def _prev_prime_fallback(n):
    n = int(n) - 1
    if n % 2 == 0:
        n -= 1
    while not ref_is_prime(n):
        n -= 2
    return n


# ----------------------------------------------------------------------------- primes
_PCACHE = {}


def odd_primes_upto(B):
    if B in _PCACHE:
        return _PCACHE[B]
    s = np.ones(B + 1, dtype=bool); s[:2] = False
    s[4::2] = False
    for i in range(3, int(B ** 0.5) + 1, 2):
        if s[i]:
            s[i * i::2 * i] = False
    P = [int(p) for p in np.nonzero(s)[0] if p > 2]
    _PCACHE[B] = P
    return P


def survival(P, k=1):
    """prod (1 - k/p): k = 1 for an ordinary candidate, k = 2 for a safe-prime pair (q, 2q+1)."""
    return float(np.prod([1.0 - k / p for p in P])) if P else 1.0


def geom_sizes(nprimes, s0=5, r=2.0):
    sizes, tot, i = [], 0, 0
    while tot < nprimes:
        s = max(1, int(round(s0 * r ** i)))
        s = min(s, nprimes - tot)
        sizes.append(s); tot += s; i += 1
    return sizes


def prod(xs):
    """Balanced product (much faster than sequential multiplication for large primorials)."""
    xs = [mpz(x) for x in xs]
    if not xs:
        return mpz(1)
    while len(xs) > 1:
        xs = [xs[i] * xs[i + 1] if i + 1 < len(xs) else xs[i] for i in range(0, len(xs), 2)]
    return xs[0]


class SegPrefilter:
    """Segmented primorial prefilter of the paper (B, S0, r)."""

    def __init__(self, B, s0=5, r=2.0):
        self.B = B
        self.P = odd_primes_upto(B)
        self.pset = set(self.P)
        self.blocks, i = [], 0
        for s in geom_sizes(len(self.P), s0, r):
            self.blocks.append(self.P[i:i + s]); i += s
        self.mods = [prod(b) for b in self.blocks]
        self.bits = [int(m).bit_length() for m in self.mods]

    def passes(self, n):
        for M in self.mods:
            if gcd(n, M) != 1:
                return False
        return True

    def passes_pair(self, q, p):
        """For safe primes: neither q nor p = 2q+1 has a small divisor."""
        for M in self.mods:
            if gcd(q, M) != 1 or gcd(p, M) != 1:
                return False
        return True

    def passes_product(self, qp):
        """One GCD with the product q*(2q+1) per segment."""
        for M in self.mods:
            if gcd(qp, M) != 1:
                return False
        return True


class TrialDiv:
    def __init__(self, B):
        self.B = B
        self.P = odd_primes_upto(B)

    def passes(self, n):
        for p in self.P:
            if n % p == 0:
                return False
        return True

    def passes_pair(self, q):
        """OpenSSL-style double sieve: q mod s != 0 and 2q+1 mod s != 0, i.e. r != 0 and r != (s-1)/2."""
        for p in self.P:
            r = q % p
            if r == 0 or r == (p - 1) // 2:
                return False
        return True


# ----------------------------------------------------------------------------- Miller–Rabin
def miller_rabin(n, t=40, seed=None):
    """Returns (probable prime?, rounds executed). Bases are derived deterministically from n (or from a seed)."""
    n_int = int(n)
    if n_int < 2:
        return False, 0
    if n_int in (2, 3):
        return True, 0
    if n_int % 2 == 0:
        return False, 0
    m = n_int - 1
    s = (m & -m).bit_length() - 1
    d = mpz(m >> s); N = mpz(n_int); Nm1 = N - 1
    rnd = random.Random((n_int & ((1 << 61) - 1)) if seed is None else seed)
    for r in range(1, t + 1):
        x = powmod(mpz(rnd.randrange(2, n_int - 1)), d, N)
        if x == 1 or x == Nm1:
            continue
        for _ in range(s - 1):
            x = x * x % N
            if x == Nm1:
                break
        else:
            return False, r
    return True, t


def fermat2(n):
    """Fermat test to base 2 (one modular exponentiation)."""
    N = mpz(n)
    return powmod(mpz(2), N - 1, N) == 1


def random_odd(bits, rnd):
    return rnd.getrandbits(bits) | (1 << (bits - 1)) | 1


def recommended_bound(bits):
    if bits <= 64: return 300
    if bits <= 256: return 1000
    if bits <= 600: return 10_000
    if bits <= 1200: return 100_000
    if bits <= 2500: return 300_000
    return 1_000_000


# ----------------------------------------------------------------------------- safe primes
def find_safe_prime(bits, seed, strategy, B=1000, t=40, pf=None, td=None):
    """Finds a safe prime p = 2q+1 of `bits` bits. All strategies see the same stream of candidates q.
    The confirmation stage is identical for all: MR(q, 1) -> Pocklington 2^(p-1) mod p == 1 -> MR(q, t-1 more rounds).
    Qaytaradi: (p, statistika)."""
    rnd = random.Random(seed)
    st = dict(candidates=0, sieve_rejected=0, modexp=0)
    top = 1 << (bits - 2)
    while True:
        q = rnd.getrandbits(bits - 1) | top | 1          # q: odd, bits-1 bits; p = 2q+1: bits bits
        st["candidates"] += 1
        p = 2 * q + 1
        if strategy == "naive":                           # trial division of q only
            ok = td.passes(q)
        elif strategy == "td_pair":                       # OpenSSL-style double sieve
            ok = td.passes_pair(q)
        elif strategy == "geom_pair":                     # separate GCDs for q and p
            ok = pf.passes_pair(mpz(q), mpz(p))
        elif strategy == "geom_product":                  # one GCD with the product q*(2q+1)
            ok = pf.passes_product(mpz(q) * mpz(p))
        elif strategy == "mr_only":
            ok = True
        else:
            raise ValueError(strategy)
        if not ok:
            st["sieve_rejected"] += 1
            continue
        # confirmation stage
        ok_q, r = miller_rabin(q, 1); st["modexp"] += r
        if not ok_q:
            continue
        st["modexp"] += 1
        if not fermat2(p):                                # Pocklington: if q is prime, p is prime <=> 2^(p-1) = 1 (mod p)
            continue
        ok_q, r = miller_rabin(q, t, seed=int(q) ^ 0xABCDEF); st["modexp"] += r
        if ok_q:
            return p, st


# ----------------------------------------------------------------------------- batch prefilter
def product_tree(xs):
    tree = [list(xs)]
    while len(tree[-1]) > 1:
        lvl = tree[-1]
        tree.append([lvl[i] * lvl[i + 1] if i + 1 < len(lvl) else lvl[i] for i in range(0, len(lvl), 2)])
    return tree


def batch_passes(cands, P_mpz):
    """Bernstein-style remainder tree: P mod n_j for all candidates, then gcd(n_j, P mod n_j)."""
    tree = product_tree([mpz(c) for c in cands])
    rems = [P_mpz % tree[-1][0]]
    for lvl in reversed(tree[:-1]):
        rems = [rems[i // 2] % lvl[i] for i in range(len(lvl))]
    return [gcd(r, n) == 1 for r, n in zip(rems, tree[0])]


# ----------------------------------------------------------------------------- parallel ishchilar
def worker_generate(args):
    """Generates k primes in one process (for multiprocessing)."""
    bits, k, seed, B = args
    pf = SegPrefilter(B)
    rnd = random.Random(seed)
    t0 = perf_counter_ns(); cand = 0; primes = []
    while len(primes) < k:
        n = random_odd(bits, rnd); cand += 1
        if not pf.passes(mpz(n)):
            continue
        if miller_rabin(n)[0]:
            primes.append(n)
    return dict(seconds=(perf_counter_ns() - t0) / 1e9, candidates=cand, primes=[str(p) for p in primes])


# ----------------------------------------------------------------------------- OpenSSL (ctypes)
class OpenSSL:
    """Access to libcrypto through ctypes: BN_generate_prime_ex and BN_check_prime."""

    def __init__(self):
        self.lib, self.path = self._load()
        if self.lib is None:
            return
        L = self.lib
        L.BN_new.restype = ctypes.c_void_p
        L.BN_free.argtypes = [ctypes.c_void_p]
        L.BN_CTX_new.restype = ctypes.c_void_p
        L.BN_CTX_free.argtypes = [ctypes.c_void_p]
        L.BN_generate_prime_ex.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_int, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
        L.BN_generate_prime_ex.restype = ctypes.c_int
        L.BN_bn2hex.argtypes = [ctypes.c_void_p]; L.BN_bn2hex.restype = ctypes.c_void_p
        L.BN_hex2bn.argtypes = [ctypes.POINTER(ctypes.c_void_p), ctypes.c_char_p]; L.BN_hex2bn.restype = ctypes.c_int
        L.BN_num_bits.argtypes = [ctypes.c_void_p]; L.BN_num_bits.restype = ctypes.c_int
        self.has_check = hasattr(L, "BN_check_prime")
        if self.has_check:
            L.BN_check_prime.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
            L.BN_check_prime.restype = ctypes.c_int
        free = getattr(L, "CRYPTO_free", None)
        if free is not None:
            free.argtypes = [ctypes.c_void_p, ctypes.c_char_p, ctypes.c_int]
        self._free = free
        try:
            L.OpenSSL_version.restype = ctypes.c_char_p
            self.version = L.OpenSSL_version(0).decode()
        except Exception:
            self.version = "noma'lum"
        self.ctx = L.BN_CTX_new()

    @staticmethod
    def _load():
        names = []
        found = ctypes.util.find_library("crypto")
        if found:
            names.append(found)
        names += ["libcrypto-3-x64.dll", "libcrypto-3.dll", "libcrypto-1_1-x64.dll", "libcrypto-1_1.dll",
                  "libcrypto.so.3", "libcrypto.so.1.1", "libcrypto.dylib"]
        dirs = [os.path.join(sys.base_prefix, "DLLs"), os.path.join(sys.base_prefix, "Library", "bin"),
                os.path.join(sys.prefix, "DLLs"), os.path.join(sys.prefix, "Library", "bin"), sys.base_prefix]
        for d in dirs:
            names += glob.glob(os.path.join(d, "libcrypto*.dll"))
        for nm in names:
            try:
                return ctypes.CDLL(nm), nm
            except OSError:
                continue
        return None, None

    @property
    def ok(self):
        return self.lib is not None

    def generate_prime(self, bits, safe=False):
        bn = self.lib.BN_new()
        try:
            if self.lib.BN_generate_prime_ex(bn, bits, 1 if safe else 0, None, None, None) != 1:
                raise RuntimeError("BN_generate_prime_ex xatosi")
            h = self.lib.BN_bn2hex(bn)
            val = int(ctypes.string_at(h).decode(), 16)
            if self._free:
                self._free(h, b"dtlib", 0)
            return val
        finally:
            self.lib.BN_free(bn)

    def check_prime(self, n):
        bn = ctypes.c_void_p(None)
        self.lib.BN_hex2bn(ctypes.byref(bn), format(int(n), "x").encode())
        try:
            return self.lib.BN_check_prime(bn, self.ctx, None) == 1
        finally:
            self.lib.BN_free(bn)


# ----------------------------------------------------------------------------- statistika
RNG = np.random.default_rng(20261101)


def boot_ci(x, nb=2000, stat=np.median):
    x = np.asarray(x, float)
    if len(x) < 2:
        return float(stat(x)), float(stat(x))
    idx = RNG.integers(0, len(x), (nb, len(x)))
    v = stat(x[idx], axis=1)
    return float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))


def time_per_item(fn, items, reps=5):
    ts = []
    for _ in range(reps):
        t0 = perf_counter_ns()
        for it in items:
            fn(it)
        ts.append((perf_counter_ns() - t0) / len(items))
    return float(np.median(ts)), ts


def mr_round_cost(bits, k=30, seed=7):
    rnd = random.Random(seed); ts = []
    for _ in range(k):
        n = random_odd(bits, rnd)
        t0 = perf_counter_ns(); miller_rabin(n, 1); ts.append(perf_counter_ns() - t0)
    return float(np.median(ts))


# ----------------------------------------------------------------------------- report
class Report:
    """Writes tables, figures and a summary to results_XX/ and zips the folder at the end."""

    def __init__(self, name, mode, seed):
        self.dir = f"results_{name}"
        os.makedirs(f"{self.dir}/tables", exist_ok=True)
        os.makedirs(f"{self.dir}/figures", exist_ok=True)
        self.name, self.items, self.t0 = name, [], time.time()
        self.summary = {"experiment": name, "mode": mode, "seed": seed, "backend": BACKEND}
        self.summary["environment"] = env_info()

    def log(self, msg):
        print(f"[{time.time() - self.t0:8.1f}s] {msg}", flush=True)

    def table(self, df, key, caption):
        df.to_csv(f"{self.dir}/tables/{key}.csv", index=False)
        md = _df_md(df)
        with open(f"{self.dir}/tables/{key}.md", "w", encoding="utf-8") as f:
            f.write(f"**{caption}**\n\n{md}\n")
        self.items.append(("table", key, caption, md))
        try:
            from IPython.display import display
            print(caption); display(df)
        except Exception:
            print(caption); print(df.to_string())

    def figure(self, fig, key, caption):
        fig.tight_layout()
        fig.savefig(f"{self.dir}/figures/{key}.png", dpi=200, bbox_inches="tight")
        fig.savefig(f"{self.dir}/figures/{key}.pdf", bbox_inches="tight")
        self.items.append(("figure", key, caption, f"figures/{key}.png"))
        import matplotlib.pyplot as plt
        plt.show(); plt.close(fig)

    def finish(self):
        self.summary["runtime_minutes"] = (time.time() - self.t0) / 60
        with open(f"{self.dir}/summary.json", "w", encoding="utf-8") as f:
            json.dump(self.summary, f, indent=2, ensure_ascii=False, default=str)
        lines = [f"# {self.name}: results", "", f"Mode: {self.summary['mode']}, runtime {self.summary['runtime_minutes']:.1f} min", ""]
        for kind, key, cap, payload in self.items:
            lines += [f"## {key}: {cap}", ""]
            lines += [payload, ""] if kind == "table" else [f"![{key}]({payload})", ""]
        lines += ["## summary.json", "", "```json", json.dumps(self.summary, indent=1, ensure_ascii=False, default=str)[:15000], "```"]
        with open(f"{self.dir}/report.md", "w", encoding="utf-8") as f:
            f.write("\n".join(lines))
        zname = f"{self.dir}.zip"
        with zipfile.ZipFile(zname, "w", zipfile.ZIP_DEFLATED) as z:
            for root, _, files in os.walk(self.dir):
                for fn in files:
                    z.write(os.path.join(root, fn))
        self.log(f"DONE. Results archive: {zname} ({os.path.getsize(zname) / 1e6:.1f} MB)")


def _fmt(v):
    if isinstance(v, (float, np.floating)):
        if v != v: return "–"
        if v == 0: return "0"
        if abs(v) >= 1e5 or abs(v) < 1e-3: return f"{v:.3e}"
        return f"{v:.4g}"
    return str(v)


def _df_md(df):
    cols = [str(c) for c in df.columns]
    out = ["| " + " | ".join(cols) + " |", "|" + "|".join(["---"] * len(cols)) + "|"]
    for _, row in df.iterrows():
        out.append("| " + " | ".join(_fmt(v) for v in row.values) + " |")
    return "\n".join(out)


def env_info():
    cpu = platform.processor() or platform.machine()
    try:
        with open("/proc/cpuinfo") as f:
            for line in f:
                if "model name" in line:
                    cpu = line.split(":", 1)[1].strip(); break
    except Exception:
        pass
    return {"date": datetime.datetime.now().isoformat(timespec="seconds"), "os": platform.platform(), "cpu": cpu,
            "cores": os.cpu_count(), "python": sys.version.split()[0], "backend": BACKEND}
