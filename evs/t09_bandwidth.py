# -*- coding: utf-8 -*-
r"""t09_bandwidth.py - EVS family BW | the RBF bandwidth median heuristic, executed VERBATIM

WHAT IS TESTED (no production file is edited; the REAL tree stays the measuring target)
================================================================================
The single number that sets the SCALE of every MMD value in the paper.  Two LIVE copies exist and
neither had ever been executed by a test (t01/t07 construct their own bw2 and pass it in):

    p6b_01_corrected_pipeline.py:267-272   ->  evs_load.BLOCK_BW_P6B
    p10a_prepare_cache.py:224-228          ->  evs_load.BLOCK_BW_P10A

Both blocks are extracted by line range and exec'd as-is (L.exec_block), with the environment the
block declares it needs: SEED, MMD_M, Dtr (raw source-train design matrix), sc (a StandardScaler
fitted on that same Dtr), np, pdist.  A name the block does not bind comes back as None - never a
substitute value.  `ib` / `Dtr_std` are requested as OPTIONAL names: a correct rewrite may name its
row selection anything, so an absent name is "not observable" (callers fall back to a derived pair
set), NOT an execution failure.  Only a missing / non-finite `bw2` counts as an exec failure.

DOCUMENTED RULE - the source of every expectation (quoted verbatim, never paraphrased from code)
================================================================================
Methods_v1.md sec.9 (L1 input-space shift signals, "Maximum mean discrepancy."), lines 223-229,
quoted verbatim (markdown bold markers included as they appear in the source):
    "MMD² was estimated with an RBF kernel[4] on source-training-standardized features, with
     bandwidth set by the median heuristic applied to source **training** pairwise squared
     distances and held fixed thereafter — the bandwidth was not scanned or selected against any
     outcome. The estimator is the unbiased U-statistic, and the source reference was fixed at
     m = 1,000 observations with K = 5 draws."

Methods_v1.md sec.15 (corrections made while the target labels remained sealed), lines 425-428:
    "a source-model preprocessing mismatch, in which some preprocessing parameters had been derived
     from pooled rather than source training data, and an MMD implementation in which the inputs
     were not standardized and the kernel bandwidth was derived from pooled source and target data."

The RBF kernel form k(x,y) = exp(-||x-y||^2 / bw2) itself is Gretton, Borgwardt, Rasch, Schoelkopf &
Smola, "A Kernel Two-Sample Test", JMLR 13 (2012) 723-773, eq. (3); it is used here only to say what
a fallback bw2 does to the kernel - the kernel estimator is t01/t07's item, not re-tested here.

METHODOLOGICAL SPLIT: (a) mathematics / frozen value vs (b) project convention
================================================================================
(a) MATHEMATICS / FROZEN VALUE, hard gates.  Median of pairwise SQUARED Euclidean distances between
    standardized source-training observations.  TWO readings are admissible for the consumed row
    set, and a correct implementation must match AT LEAST ONE of them:
      (i)  the declared draw min(MMD_M, n_train) from default_rng(SEED + 3) - the corrected spec
           keeps the 1,000-row bandwidth subsample, and it is the convention the frozen p10_cache
           bandwidths were produced under; and
      (ii) the exact full-pair median over ALL source-train rows - Methods sec.9's stated
           median-heuristic rule, and the remedy BW-07 measures.
    BW-01/01b hold the value to that identity (independent double loop, no scipy); BW-01c is a
    control on the reference itself; BW-04 pins median-of-SQUARED against three alternatives;
    BW-06 holds the two live copies bit-identical on a shared fixture.
(b) CONVENTION LAYER, WARN-class where a correct implementation may legitimately differ: the single
    SEED+3 draw and the exact subsample size (instrumented in BW-02; a draw-free full-pair
    implementation cannot confirm the draw convention -> WARN, never FAIL), the 1.0 degenerate
    fallback VALUE (BW-05b, warn_only), and the outcome-independence claim of sec.9 (BW-02c).
    The guard-FIRING itself is NOT convention: BW-05 hard-checks that the raw statistic really hit
    the degenerate branch and that the returned value IS the block's own fallback.
    A mismatch against (a) is a defect.  A mismatch against (b) means either the code or the
    documented convention is wrong - both sides are reported, with the exact source sentence.

SCOPE / NON-DUPLICATION
================================================================================
p6_02_label_free_signals.py:230 holds a THIRD copy that computes the bandwidth INSIDE the K-draw
loop from the POOLED subsample (its own v1.0 header documents that convention); it is already
measured by t07 (FORK-04b) and is not re-tested here - it is quoted once in BW-02b as the live
contrast to "held fixed thereafter".  p10b:88-92 consumes the p10a-cached scalar A["bw2"] only.
The synthetic fixtures below are deterministic; BW-07 additionally scans the frozen p10_cache
READ-ONLY (202 tasks / 101 unique source models; per pickle: bw2 + Dtr_std) to measure the REAL
subsample gap - no cache file is ever written and no other test reads it.

HONEST LIMIT OF BW-03/03b (stated so the coverage claim is not overstated)
================================================================================
The block receives `sc` and `Dtr` FROM its caller, so a behavioural test that supplies its own
environment cannot observe how the caller fitted the scaler.  The pooled-scaler counterfactual
therefore measures what the documented defect WOULD do to the value, while the caller's own
provenance is pinned separately by BW-03c (every scaler fit found in each file, at least one per
file, quoted by line).  A mutant that changes p6b:262 / p10a:221 is caught by BW-03c, not by BW-03
(verified by mutation testing during EVS v1.1 development).

EXPECTED STATUS MATRIX (REAL unmodified tree -> PATCHED tree)
================================================================================
  BW-00    PASS -> PASS   HARD infra: constants PARSED from each file and gated equal to the
                          documented SEED=20261005 / MMD_M=1000 (mutating either number in either
                          file is a hard FAIL, not a free-text change); both blocks extracted and
                          their interface census satisfied (SEED/MMD_M are NOT required tokens - a
                          full-pair implementation legitimately does not consume MMD_M)
  BW-01    PASS -> PASS   HARD value identity, n_train > MMD_M: the value must equal the median over
                          at least one admissible pair set (the declared SEED+3 / min(MMD_M,n)
                          draw, or the exact full-pair median); a bound `ib` is reported, an absent
                          one is not an error (the declared draw is reconstructed)
  BW-01b   PASS -> PASS   HARD value identity, n_train < MMD_M: value == all-pairs median exactly
                          (no Monte Carlo, no draw dependence); "consumed set was a full
                          permutation" is reported as a characterisation, not gated
  BW-01c   PASS -> PASS   control on the reference itself (chunked vs pure-Python scalar double loop)
  BW-02    PASS -> PASS   HARD: bw2 bit-reproducible, invariant to ambient global RNG state and to
                          MMD_K.  WARN-class convention (never FAIL): exactly one default_rng(SEED+3)
                          + one choice(n, min(MMD_M,n)) + one pdist call, bound `ib` equal to that
                          draw - a draw-free full-pair implementation cannot confirm the draw
                          convention (WARN)
  BW-02b   PASS -> PASS   "held fixed thereafter": the scalar is produced once per task, outside the
                          K loop, and no consumer reassigns it (quoted call sites)
  BW-02c   PASS -> PASS   "not selected against any outcome": no outcome/target name is in scope
                          (compile-time scan) and bw2 is unchanged when y / Dte vary in the env
  BW-03    PASS -> PASS   HARD: production is on the SOURCE-TRAIN-only standardization side and its
                          value is reproduced independently over the consumed pair set; the pooled-fit
                          scaler counterfactual differs (quantified)
  BW-03b   PASS -> PASS   HARD: further counterfactuals: unstandardized distances; pooled source+
                          target rows (both quantified, both distinguishable from production)
  BW-03c   PASS -> PASS   static conformance of the caller context: every scaler fit found in each
                          file is fitted on Dtr (the source-train split, quoted; `X=Dtr` accepted);
                          at least one such fit per file - supporting evidence for BW-03/03b
  BW-04    PASS -> PASS   HARD: production's position = median of SQUARED distances over the consumed
                          pair set; the other three alternatives quantified; monotone sqrt reported
  BW-05    PASS -> PASS   HARD: the degenerate guard FIRES on both triggers - on the all-identical
                          fixture the raw median is exactly 0 and on the NaN-majority fixture the
                          raw np.median is non-finite, and BOTH returned values equal the block's own
                          fallback literal parsed from its text; a NaN-skipping median (raw NaN,
                          returned = finite nan-median != fallback) must go RED
  BW-05b   PASS -> PASS   warn_only characterisation: the fallback constant is exactly 1.0 and is
                          NOT scale-matched (the kernel collapses towards the zero kernel) - a future
                          mismatch degrades to WARN, never FAIL, because the value is a convention
  BW-06    PASS -> PASS   p6b vs p10a parity on a shared fixture, both sizes; index vectors compared
                          only when both copies bind them
  BW-07    WARN -> WARN   warn_only characterisation on the REAL frozen corpus (p10_cache): the
                          min(MMD_M, n_train) subsample median deviates from the EXACT full-pair
                          median by up to +52.50% (12/101 source models beyond 10%, 3/101 beyond
                          25%) - convention layer, never FAIL; the subsample is DELIBERATELY KEPT
                          by the correction spec (CS-06), and the WARN is NUMERIC (5% negligible
                          criterion) so disclosing the subsample cannot clear it - advisory only

v1.3 acceptance property (the two manufactured FAILUREs this round repaired): a full-pair-median
implementation (BW-07's own recommended remedy) and a bit-identical rewrite that renames the
block's internal index vector BOTH run green - no hard gate may reject them.

Exit code is always 0; every test is wrapped by L.guard / guard_soft / guard_warn so an exception
becomes one FAIL line.
OUTPUT CONTRACT: RESULT|test_id|PASS/FAIL/WARN|expected (source)|observed|note
"""
import os
import pickle
import re
import sys
import traceback
import warnings

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
from scipy.spatial.distance import pdist
from sklearn.preprocessing import StandardScaler

import evs_load as L

FAM = "BW"
TOL_REL = 1e-9          # relative tolerance for an O(n^2) sum re-ordered vs scipy.pdist
BW01_FULLPAIR_MAX_PAIRS = 8_000_000   # compute the exact full-pair candidate below this many pairs
                                      # (the n=1200 fixture needs 719,400)

DEF_MED = ('Methods_v1.md sec.9 (lines 224-225): "bandwidth set by the median heuristic applied to '
           'source **training** pairwise squared distances and held fixed thereafter — the bandwidth '
           'was not scanned or selected against any outcome"')
DEF_FIX = ('Methods_v1.md sec.15: corrected defect - "an MMD implementation in which the inputs were '
           'not standardized and the kernel bandwidth was derived from pooled source and target data"')
DEF_KER = ('Gretton et al. JMLR 13 (2012) 723-773 eq.(3): the RBF kernel k(x,y)=exp(-norm(x-y)^2/bw2) is '
           'defined only for a finite positive bandwidth bw2')

L.banner(f"{FAM} - t09_bandwidth.py | the median-heuristic bandwidth, {L.block_source(L.BLOCK_BW_P6B)}"
         f" + {L.block_source(L.BLOCK_BW_P10A)}")


# ============================================================ production text (read-only)
def read_lines(path):
    with open(path, encoding="utf-8") as fh:
        return fh.read().split("\n")


_RE_SEED = re.compile(r"^SEED\s*=\s*(\d+)\s*$")
_RE_MMD = re.compile(r"^MMD_M\s*,\s*MMD_K\s*,\s*NBIN\s*=\s*(\d+)\s*,\s*(\d+)\s*,\s*(\d+)\s*$")


def read_constants(path):
    """Parse the module-level constants OUT of the production text (never trust a hard-coded copy)."""
    d = {}
    for ln in read_lines(path)[:80]:
        m = _RE_SEED.match(ln)
        if m:
            d["SEED"] = int(m.group(1))
        m = _RE_MMD.match(ln)
        if m:
            d["MMD_M"], d["MMD_K"], d["NBIN"] = int(m.group(1)), int(m.group(2)), int(m.group(3))
    return d


def find_line(path, pattern):
    """First line matching `pattern`; returns (lineno, stripped text) or (None, None)."""
    rx = re.compile(pattern)
    for i, ln in enumerate(read_lines(path), 1):
        if rx.search(ln):
            return i, ln.strip()
    return None, None


try:                                     # the module must still reach its RESULT lines (exit 0) if a
    C_P6B = read_constants(L.P6B)        # production file was moved or corrupted, instead of dying
    C_P10A = read_constants(L.P10A)      # with a traceback before any test has run
    _SRC_ERR = None
except Exception as _e:
    C_P6B, C_P10A, _SRC_ERR = {}, {}, f"{type(_e).__name__}: {_e}"
SEED_DOC, MMD_M_DOC = 20261005, 1000     # documented constants; BW-00 checks the files still carry them
SEED = C_P6B.get("SEED", SEED_DOC)
MMD_M = C_P6B.get("MMD_M", MMD_M_DOC)


# ============================================================ deterministic synthetic fixtures
def source_fixture(n, seed):
    """Synthetic source-TRAIN design matrix, shaped like production `design(Xs, Hs)`: 2 continuous
    columns with very different scales (age in years; a skewed positive lab) and 5 one-hot dummies
    (4 t-bins, 3 past-AST states, reference level dropped) -> n x 7."""
    rng = np.random.default_rng(seed)
    age = rng.normal(58.0, 17.0, n)
    lab = np.exp(rng.normal(0.0, 0.8, n))
    tb = rng.integers(0, 4, n)
    sa = rng.integers(0, 3, n)
    cols = [age, lab]
    cols += [(tb == k).astype(float) for k in (1, 2, 3)]
    cols += [(sa == k).astype(float) for k in (1, 2)]
    return np.column_stack(cols)


def target_fixture(n, seed):
    """Synthetic TARGET design matrix: same 7 columns, shifted means/scales and a different case
    mix - the domain shift that a pooled (source+target) standardization would absorb."""
    rng = np.random.default_rng(seed)
    age = rng.normal(50.0, 20.5, n)
    lab = np.exp(rng.normal(0.60, 1.05, n))
    tb = rng.choice(4, n, p=[0.30, 0.32, 0.22, 0.16])
    sa = rng.choice(3, n, p=[0.46, 0.32, 0.22])
    cols = [age, lab]
    cols += [(tb == k).astype(float) for k in (1, 2, 3)]
    cols += [(sa == k).astype(float) for k in (1, 2)]
    return np.column_stack(cols)


N_BIG, N_SMALL = 1200, 450
DTR_BIG = source_fixture(N_BIG, 90001)       # n_train > MMD_M -> subsample = MMD_M rows
DTR_SMALL = source_fixture(N_SMALL, 90002)   # n_train < MMD_M -> subsample = all rows
DTGT = target_fixture(1500, 90003)           # used by the pooled counterfactuals only


# ============================================================ independent distance references
def pair_sqdist_chunked(Z):
    """Condensed squared Euclidean distances over UNORDERED pairs i<j, from the definition
    ||x-y||^2 = sum_k (x_k - y_k)^2, by an explicit double loop (outer loop over i, inner loop
    vectorised over j).  scipy is never used by this reference."""
    n = len(Z)
    out = np.empty(n * (n - 1) // 2, dtype=float)
    k = 0
    for i in range(n - 1):
        d = ((Z[i + 1:] - Z[i]) ** 2).sum(axis=1)
        out[k:k + len(d)] = d
        k += len(d)
    return out


def pair_sqdist_scalar(Z):
    """Pure-Python scalar double loop over the same unordered pairs - the control on the reference."""
    out = []
    for i in range(len(Z)):
        zi = Z[i]
        for j in range(i + 1, len(Z)):
            zj = Z[j]
            s = 0.0
            for k in range(len(zi)):
                d = zi[k] - zj[k]
                s += d * d
            out.append(s)
    return np.asarray(out)


def med_rel(a, b):
    return abs(a - b) / max(abs(b), 1e-300)


# ============================================================ instrumentation proxies
class _CountingPdist:
    def __init__(self, fn):
        self.fn, self.calls = fn, 0

    def __call__(self, *a, **kw):
        self.calls += 1
        return self.fn(*a, **kw)


class _GenProxy:
    def __init__(self, gen, log):
        self._gen, self._log = gen, log

    def choice(self, *a, **kw):
        self._log["choice_calls"] += 1
        self._log["choice_args"] = (a[0] if len(a) > 0 else kw.get("a"),
                                    a[1] if len(a) > 1 else kw.get("size"))
        return self._gen.choice(*a, **kw)

    def __getattr__(self, name):
        return getattr(self._gen, name)


class _RandomProxy:
    def __init__(self, real, log):
        self._real, self._log = real, log

    def default_rng(self, *a, **kw):
        self._log["default_rng_calls"] += 1
        self._log["default_rng_args"] = (tuple(a), dict(kw))
        return _GenProxy(self._real.default_rng(*a, **kw), self._log)

    def __getattr__(self, name):
        return getattr(self._real, name)


class _NpProxy:
    """Forwards everything to numpy but counts np.random.default_rng() and the choice() calls on the
    generator it returns.  The block source stays VERBATIM; only the bound name `np` differs, which
    is part of the environment the block itself declares it needs."""

    def __init__(self, real, log):
        self._real, self._log = real, log
        self.random = _RandomProxy(real.random, log)
        self.isfinite = real.isfinite
        self.median = real.median

    def __getattr__(self, name):
        return getattr(self._real, name)


def blank_log():
    return {"default_rng_calls": 0, "choice_calls": 0, "default_rng_args": None, "choice_args": None}


def run_block(block, Dtr, sc, seed=None, mmd_m=None, pdist_fn=pdist, np_mod=np, extra=None):
    """exec the production bandwidth block VERBATIM with the environment it declares.
    `bw2` is the only hard requirement; `ib` / `Dtr_std` are OPTIONAL (absent -> None)."""
    env = dict(SEED=SEED if seed is None else seed,
               MMD_M=MMD_M if mmd_m is None else mmd_m,
               Dtr=Dtr, sc=sc, np=np_mod, pdist=pdist_fn)
    if extra:
        env.update(extra)
    got, src = L.exec_block(block, env, names=["bw2", "ib", "Dtr_std"])
    return got, src


def sc_of(Dtr):
    return StandardScaler().fit(Dtr)


def try_run(block, D, sc):
    """Run the block on a fixture; on exec failure return a sentinel so EVERY test reports its own
    FAIL line instead of the module dying before the rest of the family has run."""
    try:
        return run_block(block, D, sc)
    except Exception as e:
        return ({"bw2": float("nan"), "ib": None, "Dtr_std": None},
                f"<EXEC FAILED: {type(e).__name__}: {e}>")


def exec_failed(got):
    """True when the block itself failed: the exec raised (sentinel), or no finite `bw2` came out.
    An absent / renamed `ib` / `Dtr_std` is NOT an execution failure (v1.3 repair): a correct
    implementation may name its row selection anything, and callers fall back to derived
    quantities instead of emitting a phantom 'block exec failed ... see BW-00' while BW-00 is green."""
    v = got.get("bw2")
    try:
        return v is None or not np.isfinite(v)
    except Exception:
        return True


def sel_of(got, n):
    """The block's own row-index vector when it is BOUND and usable, else None.
    Never raises and never treats a missing/renamed name as an error."""
    ib = got.get("ib")
    if ib is None:
        return None
    try:
        ib = np.asarray(ib)
        if ib.ndim != 1 or ib.dtype.kind not in "iu" or len(ib) < 2:
            return None
        if ib.min() < 0 or ib.max() >= n:
            return None
        return ib
    except Exception:
        return None


def declared_draw(n):
    """The declared draw convention: min(MMD_M, n) rows from np.random.default_rng(SEED + 3)."""
    return np.random.default_rng(SEED + 3).choice(n, min(MMD_M, n), replace=False)


def pair_sel(got, n):
    """(indices, how) - the row set over which to recompute the block's value:
    the block's own bound `ib` when usable; otherwise the declared draw reconstructed from
    default_rng(SEED + 3) (which is the full row set whenever min(MMD_M, n) == n)."""
    sel = sel_of(got, n)
    if sel is not None:
        return sel, "block's own bound index vector"
    return declared_draw(n), "declared draw reconstructed (SEED+3, min(MMD_M, n)); `ib` not bound"


_RE_FB = re.compile(r"bw2\s*=\s*([0-9]+(?:\.[0-9]+)?(?:[eE][+-]?\d+)?)\s*$")


def fallback_literal(block):
    """The degenerate-guard fallback VALUE, parsed OUT of the block's own text (never assumed):
    the `bw2 = <number>` assignment on the non-finite/<=0 guard line.  None when not parseable."""
    for ln in L.grab(*block).split("\n"):
        if "isfinite" in ln and re.search(r"bw2\s*=", ln):
            m = _RE_FB.search(ln)
            if m:
                return float(m.group(1))
    return None


def guard_warn(tid, fn, expected, note=""):
    """Like L.guard, but the outcome is a WARN at worst: used where the check pins a deliberately
    documented CONVENTION / current state (rule 5: such a check is a characterisation, not a test)."""
    try:
        ok, observed = fn()
    except Exception as e:
        return L.emit(tid, "WARN", expected,
                      f"EXCEPTION {type(e).__name__}: {e}",
                      (note + " ; " + traceback.format_exc().splitlines()[-1]).strip(" ;"))
    return L.check(tid, ok, expected, observed, note, warn_only=True)


def guard_soft(tid, fn, expected, note=""):
    """Three-state guard: hard failure -> FAIL; hard ok but a documented convention unconfirmed ->
    WARN; otherwise PASS.  `fn` returns (hard_ok, soft_ok, observed).  Used where an assertion has a
    mathematical core (hard) AND a project-convention shell that a correct alternative
    implementation need not satisfy (rule 4: no manufactured FAIL against a correct tree)."""
    try:
        hard_ok, soft_ok, observed = fn()
    except Exception as e:
        return L.emit(tid, "FAIL", expected, f"EXCEPTION {type(e).__name__}: {e}",
                      (note + " ; " + traceback.format_exc().splitlines()[-1]).strip(" ;"))
    if not hard_ok:
        return L.emit(tid, "FAIL", expected, observed, note)
    return L.check(tid, soft_ok, expected, observed, note, warn_only=True)


# ============================================================ frozen-corpus measurement (BW-07)
CACHE_DIR = os.path.join(L.OUT, "p10_cache")
BW07_EXACT_MAX_PATTERN_PAIRS = 40_000_000   # multiplicity-expanded pairs; above this use the sampler
BW07_SAMPLE_PAIRS = 20_000_000              # unbiased unordered-pair sample on the fallback path
BW07_NEGLIGIBLE_PCT = 5.0                   # abs(gap) at or below this = declared negligible
BW07_BRUTE_MAX_ROWS = 4000                  # largest n whose FULL pair set is brute-forced as a control
BW07_CTRL_SAMPLE_PAIRS = 2_000_000
BW07_FRESH_DRAWS = 30


class _PassThroughScaler:
    """Identity transform, used ONLY to feed the cached standardized matrix back into the live block:
    the pickle stores the exact array the frozen pdist consumed, so the scaler step must not alter a
    single value (re-fitting a StandardScaler here would perturb it at the 1e-16 level)."""

    def transform(self, X):
        return np.asarray(X, dtype=float)


def sample_pair_median(Z, n_pairs, seed):
    """Median of `n_pairs` UNORDERED row pairs drawn UNIFORMLY with replacement from the FULL row set
    ((u, v) uniform, u == v rejected, i = min(u, v), j = max(u, v)); the unbiased-sampling fallback
    of BW-07 for a source model too large for the exact multiplicity expansion."""
    n = len(Z)
    rng = np.random.default_rng(seed)
    out = np.empty(n_pairs, dtype=float)
    got, chunk = 0, min(2_000_000, n_pairs)
    while got < n_pairs:
        m = min(chunk, n_pairs - got)
        u = rng.integers(0, n, m)
        v = rng.integers(0, n, m)
        keep = u != v
        u, v = u[keep], v[keep]
        i, j = np.minimum(u, v), np.maximum(u, v)
        d2 = ((Z[i] - Z[j]) ** 2).sum(axis=1)
        out[got:got + len(d2)] = d2
        got += len(d2)
    return float(np.median(out))


def full_pair_median_exact(Z):
    """EXACT median of the squared distances over ALL n(n-1)/2 unordered pairs of Z, from the
    unique-row multiplicity expansion: a pair of distinct patterns (i, j) contributes w_i*w_j pairs
    at distance ||u_i - u_j||^2, and a repeated pattern contributes w_i*(w_i-1)/2 pairs at distance 0.
    The ordered pair multiset is therefore reproduced exactly, and positions (N-1)//2 and N//2
    reproduce np.median's even-count averaging convention.  Returns (median, total_pairs, m)."""
    U, w = np.unique(Z, axis=0, return_counts=True)
    m = len(U)
    vs, ws = [], []
    for i in range(m - 1):
        vs.append(((U[i + 1:] - U[i]) ** 2).sum(axis=1))
        ws.append(w[i] * w[i + 1:].astype(np.float64))
    vals = np.concatenate(vs) if vs else np.zeros(0)
    mult = np.concatenate(ws) if ws else np.zeros(0)
    zero = float((w.astype(np.int64) * (w.astype(np.int64) - 1) // 2).sum())
    order = np.argsort(vals, kind="mergesort")
    vals, mult = vals[order], mult[order]
    allv = np.concatenate(([0.0], vals))
    allw = np.concatenate(([zero], mult))
    keep = allw > 0
    allv, allw = allv[keep], allw[keep]
    cum = np.cumsum(allw)
    N = int(cum[-1])
    k1, k2 = (N - 1) // 2, N // 2
    return (0.5 * (float(allv[np.searchsorted(cum, k1, side="right")])
                   + float(allv[np.searchsorted(cum, k2, side="right")])), N, m)


def scan_frozen_corpus():
    """READ-ONLY scan of the frozen p10_cache.  For every task: recompute the full-pair median
    (exact multiplicity expansion; sampled fallback above the pattern-pair budget) and reproduce the
    frozen bw2 through the LIVE production block fed with the cached Dtr_std (pass-through scaler),
    recording gap% = 100*(frozen - full)/full.  A pathological live block (e.g. a full-pair rewrite
    whose pdist cost explodes on a 60k-row cached matrix) must not kill the scan: its live value is
    recorded as non-reproducible and the WARN characterisation carries on.  Identical source models
    (same organism/drug/source site) are matched on the fingerprint (n, Z.sum(), Z[0,0], bw2) and
    are NOT re-executed.  Returns (files, models), models = {key: [record, ...]}."""
    files = sorted(f for f in os.listdir(CACHE_DIR) if f.endswith(".pkl"))
    models = {}
    for fn in files:
        with open(os.path.join(CACHE_DIR, fn), "rb") as fh:
            d = pickle.load(fh)
        Z = np.asarray(d["Dtr_std"], dtype=float)
        n, bw = len(Z), float(d["bw2"])
        key = (str(d["organism"]), str(d["drug"]), str(d["source_site"]))
        fp = (n, float(Z.sum()), float(Z[0, 0]))
        prev = models.get(key)
        if prev and prev[0]["fp"] == fp and prev[0]["bw"] == bw:
            prev.append(dict(prev[0], task=str(d["task_id"]), fn=fn, dup=True))
            del d
            continue
        m_uniq = len(np.unique(Z, axis=0))
        if m_uniq * (m_uniq - 1) // 2 <= BW07_EXACT_MAX_PATTERN_PAIRS:
            full, N, _ = full_pair_median_exact(Z)
            method = "exact"
        else:
            full, N, method = sample_pair_median(Z, BW07_SAMPLE_PAIRS, 20261008), \
                n * (n - 1) // 2, "sampled"
        live, live_err = float("nan"), None
        try:
            got, _ = L.exec_block(L.BLOCK_BW_P6B, dict(SEED=SEED, MMD_M=MMD_M, Dtr=Z,
                                                       sc=_PassThroughScaler(), np=np, pdist=pdist),
                                  names=["bw2"])
            if got["bw2"] is not None:
                live = float(got["bw2"])
        except Exception as e:
            live_err = f"{type(e).__name__}: {e}"
        rel_live = (abs(live - bw) / max(abs(bw), 1e-300)) if np.isfinite(live) else float("inf")
        models.setdefault(key, []).append(dict(
            key=key, task=str(d["task_id"]), fn=fn, n=n, m=m_uniq, bw=bw, live=live,
            live_err=live_err, rel_live=rel_live, full=full, pairs=N, method=method,
            gap=100.0 * (bw - full) / full, fp=fp, dup=False))
        del d
    return files, models


def main():
    # ------------------------------------------------------------------ BW-00 infra
    def bw00():
        if _SRC_ERR:
            return False, f"production text could not be read/parsed: {_SRC_ERR}"
        parsed = all(C_P6B.get(k) is not None for k in ("SEED", "MMD_M")) and \
                 all(C_P10A.get(k) is not None for k in ("SEED", "MMD_M"))
        same = (C_P6B.get("SEED") == C_P10A.get("SEED")
                and C_P6B.get("MMD_M") == C_P10A.get("MMD_M"))
        # the documented constants are GATED BY VALUE, not by text presence: SEED drives every draw
        # and Methods sec.9 documents the reference size m = 1,000, so a mutation of either number
        # in EITHER file (or a silent divergence between the two) silently rescales the whole family;
        # BW-01..BW-06 then read the drifted values as their own expectation, which is exactly why
        # this gate must be a hard FAIL and not a free-text note
        doc_ok = (SEED == SEED_DOC and MMD_M == MMD_M_DOC)
        src_ok, census = {}, {}
        for blk, tag in ((L.BLOCK_BW_P6B, "p6b"), (L.BLOCK_BW_P10A, "p10a")):
            s = L.grab(*blk)                       # the census always reads the live text
            try:
                got, _ = L.exec_block(blk, dict(SEED=SEED, MMD_M=MMD_M,
                                                Dtr=np.zeros((12, 3)),
                                                sc=StandardScaler().fit(np.zeros((12, 3))),
                                                np=np, pdist=pdist), names=["bw2"])
                raised = None
            except Exception as e:                 # a defect block may not even run: report it
                got, raised = {}, f"{type(e).__name__}: {e}"
            # interface conformance only (the ALGORITHM is judged behaviourally, never by text):
            # the extracted region must define bw2 from Dtr / sc.transform / pdist.  SEED / MMD_M
            # are NOT required tokens - a full-pair-median implementation legitimately never
            # consumes MMD_M - but their presence is still reported in the census.
            census[tag] = {tok: (tok in s) for tok in ("bw2 =", "Dtr", "sc.transform", "pdist",
                                                       "SEED", "MMD_M")}
            if raised:
                census[tag]["block raised on the degenerate fixture"] = raised
                src_ok[tag] = False
            else:
                src_ok[tag] = (all(census[tag][t] for t in ("bw2 =", "Dtr", "sc.transform", "pdist"))
                               and got["bw2"] is not None and np.isfinite(got["bw2"]))
        return (parsed and same and all(src_ok.values()) and doc_ok,
                f"parsed constants: p6b SEED={C_P6B.get('SEED')} MMD_M={C_P6B.get('MMD_M')}; "
                f"p10a SEED={C_P10A.get('SEED')} MMD_M={C_P10A.get('MMD_M')}; both files agree: {same}; "
                f"documented SEED=={SEED_DOC} and MMD_M=={MMD_M_DOC}: {doc_ok} "
                f"({'still current' if doc_ok else 'CHANGED vs doc -> hard FAIL'}); "
                f"interface census {census}; required interface (bw2<-Dtr/sc.transform/pdist) present "
                f"and finite on a degenerate fixture: {src_ok}")

    L.guard("BW-00", bw00,
            "environment contract (HARD): the block's constants and inputs come from the production "
            "text (a name the block does not bind comes back None - never substituted) AND the "
            "constants PARSED from both files still equal the documented SEED=20261005 (drives every "
            "draw) and MMD_M=1000 (Methods sec.9: reference m = 1,000): mutating either number in "
            "either file, or letting the two files diverge, is a hard FAIL - not a change of free "
            "text; plus the block must define bw2 from Dtr / sc.transform / pdist",
            "census: two live median-heuristic copies (p6b:267-272, p10a:224-228); p6_02:230 is a THIRD, "
            "pooled-bandwidth copy already measured by t07 FORK-04b; the documented values are the "
            "source of BW-01..BW-06's own expectations, so this gate protects the whole family")

    sc_big = sc_of(DTR_BIG)
    sc_small = sc_of(DTR_SMALL)
    got_big, src_big = try_run(L.BLOCK_BW_P6B, DTR_BIG, sc_big)
    got_small, _ = try_run(L.BLOCK_BW_P6B, DTR_SMALL, sc_small)
    bw_prod = got_big["bw2"]
    SEL_BIG, SEL_HOW = pair_sel(got_big, len(DTR_BIG))   # consumed set: bound ib, else declared draw
    Z_BIG = sc_big.transform(DTR_BIG)
    N_PAIRS_BIG = len(DTR_BIG) * (len(DTR_BIG) - 1) // 2

    # ------------------------------------------------------------------ BW-01 value identity (n > MMD_M)
    def bw01():
        if exec_failed(got_big):
            return False, ("the block produced no finite bw2 on the n=1200 fixture - see BW-00 for "
                           "the infrastructure reason")
        n = len(DTR_BIG)
        k_doc = min(MMD_M, n)
        sel = sel_of(got_big, n)
        ref_draw = declared_draw(n)
        cands = [(f"declared draw default_rng(SEED+3).choice({n}, min(MMD_M,{n})={k_doc})",
                  float(np.median(pair_sqdist_chunked(Z_BIG[ref_draw]))))]
        if N_PAIRS_BIG <= BW01_FULLPAIR_MAX_PAIRS:
            cands.append(("EXACT full-pair median over all source-train rows (Methods sec.9's rule)",
                          float(np.median(pair_sqdist_chunked(Z_BIG)))))
        rels = [(lab, v, med_rel(v, bw_prod)) for lab, v in cands]
        best_lab, best_v, best_rel = min(rels, key=lambda t: t[2])
        ok = best_rel <= TOL_REL
        if sel is None:
            ihow = ("index vector `ib` not bound (a correct rewrite may name it anything) - the "
                    "reconstructed declared draw listed below is the admissible reading it must match")
        else:
            ihow = (f"index vector bound: len={len(sel)}, == declared draw: "
                    f"{np.array_equal(np.asarray(sel), ref_draw)}, median over it = "
                    f"{float(np.median(pair_sqdist_chunked(Z_BIG[sel])))!r}")
        return (ok, f"n_train={n} min(MMD_M,n_train)={k_doc} full_pairs={N_PAIRS_BIG}; {ihow}; "
                    f"prod_bw2={bw_prod!r}; candidates: "
                    + "; ".join(f"{lab}: value={v!r} rel={rel:.2e}" for lab, v, rel in rels)
                    + f"; best={best_lab!r} rel={best_rel:.2e}")

    L.guard("BW-01", bw01,
            f"{DEF_MED} --value identity (HARD): the stored bw2 must equal the median of the squared "
            f"distances over AT LEAST ONE admissible pair set, recomputed by an independent path "
            f"(sc.transform + explicit double loop + np.median, no scipy): either the declared "
            f"min(MMD_M, n_train) draw from default_rng(SEED+3) or the exact full-pair median over all "
            f"source-train rows; an implementation that consumes every row (full-pair median) and one "
            f"that renames its internal index vector are both correct and must not go red; the bound "
            f"`ib`, when present, and the draw identity are reported",
            "the two admissible readings are (i) the declared draw kept by the corrected spec and "
            "(ii) Methods sec.9's full-pair median heuristic, which BW-07 measures and recommends; "
            "matching neither is a hard defect (e.g. a different RNG stream changes every frozen scale)")

    # ------------------------------------------------------------------ BW-01b value identity (n < MMD_M)
    def bw01b():
        if exec_failed(got_small):
            return False, ("the block produced no finite bw2 on the n=450 fixture (e.g. the "
                           "min(MMD_M, len(Dtr)) small-n branch is gone) - see BW-00")
        n = len(DTR_SMALL)
        Z = sc_small.transform(DTR_SMALL)
        ref_all = float(np.median(pair_sqdist_chunked(Z)))
        sel = sel_of(got_small, n)
        if sel is None:
            draw_txt = "consumed set: not observable (`ib` not bound) - value held to the full-row median"
            perm = None
        else:
            perm = bool(np.array_equal(np.sort(np.asarray(sel)), np.arange(n)))
            draw_txt = (f"consumed set: len={len(sel)} == min(MMD_M,n_train)={min(MMD_M, n)}; "
                        f"full permutation of all rows: {perm}")
        rel = med_rel(got_small["bw2"], ref_all)
        return (rel <= TOL_REL,
                f"n_train={n} {draw_txt}; prod_bw2={got_small['bw2']!r} all-pairs median={ref_all!r} "
                f"rel={rel:.3e} pairs={n * (n - 1) // 2}")

    L.guard("BW-01b", bw01b,
            f"{DEF_MED} --when n_train <= MMD_M the min() branch consumes every row (or the readonly "
            f"full-pair rule does outright), so the production value must equal the median over ALL "
            f"source-train pairs exactly (no Monte Carlo, no draw dependence): the HARD gate is that "
            f"value identity; 'the consumed set is a full permutation' is a characterisation (an "
            f"absent `ib` is not an error)",
            "small-n branch of the same verbatim block")

    # ------------------------------------------------------------------ BW-01c control on the reference
    def bw01c():
        Zc = sc_small.transform(DTR_SMALL[:120])
        a = pair_sqdist_chunked(Zc)
        b = pair_sqdist_scalar(Zc)
        rel = float(np.max(np.abs(a - b) / np.maximum(np.abs(b), 1e-300)))
        return (len(a) == len(b) == 120 * 119 // 2 and rel <= 1e-12,
                f"pairs={len(a)} max_rel_diff={rel:.3e} chunked_median={float(np.median(a))!r} "
                f"scalar_double_loop_median={float(np.median(b))!r}")

    L.guard("BW-01c", bw01c,
            "reference self-check: the chunked double-loop reference must equal a pure-Python scalar "
            "double loop on the same 120-row slice (control on the reference itself, not on production)",
            "guards against a wrong reference silently passing BW-01")

    # ------------------------------------------------------------------ BW-02 determinism / draw identity
    def bw02():
        log = blank_log()
        pd_ct = _CountingPdist(pdist)
        got, _ = run_block(L.BLOCK_BW_P6B, DTR_BIG, sc_big, np_mod=_NpProxy(np, log), pdist_fn=pd_ct)
        n, k_doc = len(DTR_BIG), min(MMD_M, len(DTR_BIG))
        ib_ref = declared_draw(n)
        np.random.seed(1234)
        np.random.random(97)                                    # perturb the ambient global stream
        got_amb, _ = run_block(L.BLOCK_BW_P6B, DTR_BIG, sc_big, extra=dict(MMD_K=100000))
        np.random.seed(4321)
        np.random.random(5)
        got_amb2, _ = run_block(L.BLOCK_BW_P6B, DTR_BIG, sc_big)
        # ---- HARD: the frozen scalar is deterministic and invariant to ambient state / MMD_K -----
        vals = [got["bw2"], got_amb["bw2"], got_amb2["bw2"]]
        hard_ok = (all(v is not None and np.isfinite(v) for v in vals) and vals[0] == vals[1] == vals[2])
        ibs = [got.get("ib"), got_amb.get("ib"), got_amb2.get("ib")]
        if all(x is not None for x in ibs):
            ib_det = bool(np.array_equal(np.asarray(ibs[0]), np.asarray(ibs[1]))
                          and np.array_equal(np.asarray(ibs[0]), np.asarray(ibs[2])))
            hard_ok = hard_ok and ib_det
            ib_cmp_txt = f"bound index vector identical across the three execs: {ib_det}"
        else:
            ib_cmp_txt = "bound index vector: not observable (name not bound / not bound everywhere)"
        # ---- WARN-class convention: the ONE declared draw (SEED+3, min(MMD_M,n)) ------------------
        seed_arg = log["default_rng_args"][0] if log["default_rng_args"] else ()
        ib = got.get("ib")
        if ib is None:
            ib_claim, ib_ok = "not observable (name not bound)", True
        else:
            ib_ok = bool(np.array_equal(np.asarray(ib), ib_ref))
            ib_claim = str(ib_ok)
        conv_ok = (log["default_rng_calls"] == 1 and log["choice_calls"] == 1 and pd_ct.calls == 1
                   and seed_arg == (SEED + 3,) and log["choice_args"] == (n, k_doc) and ib_ok)
        return (hard_ok, conv_ok,
                f"HARD: bw2 identical across three execs under advanced ambient global RNG state and "
                f"MMD_K variation: {vals[0] == vals[1] == vals[2]} ({vals[0]!r}); {ib_cmp_txt}; "
                f"WARN-class convention (unconfirmed -> WARN, never FAIL): "
                f"np.random.default_rng calls={log['default_rng_calls']} args={seed_arg} "
                f"(SEED+3={SEED + 3}); Generator.choice calls={log['choice_calls']} "
                f"args={log['choice_args']} (declared: ({n}, {k_doc})); pdist calls={pd_ct.calls}; "
                f"ib == default_rng(SEED+3) draw: {ib_claim}; convention verdict: {conv_ok}")

    guard_soft("BW-02", bw02,
               f"{DEF_MED} --HARD part: the stored scalar is bit-reproducible and independent of the "
               f"ambient global RNG state and of MMD_K. WARN-class convention (never FAIL): the "
               f"subsample is drawn ONCE per task by np.random.default_rng(SEED + 3), one "
               f"Generator.choice(min(MMD_M, n_train)) call and one pdist call inside the block",
               "the SEED+3 single-draw identity is a project convention: a draw-free full-pair "
               "implementation cannot confirm it (WARN), and a block drawing from another stream is "
               "reported as WARN here while BW-01's value identity fails it (its value matches "
               "neither admissible reading).  Implication of the hard part: bw2 is a deterministic "
               "function of (the block's stream, n_train, the fitted scaler, the source-train rows), "
               "so no per-run randomness enters any MMD value.  v1.3 repair: a renamed index vector "
               "must not print 'ib == draw: False' when it is simply unobservable")

    # ------------------------------------------------------------------ BW-02b held fixed across the K draws
    def bw02b():
        mmd_fn = L.load_def(L.P6B, "mmd_rbf")
        body = mmd_fn.__evs_src__
        sig = re.search(r"^def\s+mmd_rbf\(([^)]*)\)", body, re.M)
        takes = sig is not None and "bw2" in sig.group(1)
        reassigned = re.search(r"\bbw2\s*=(?!=)", body) is not None
        i_call, call_site = find_line(L.P6B, r"mmd_rbf\(Ds_std, Dt_std, M\[.bw2.\]")
        i_k, k_line = find_line(L.P10A, r"Ka = np\.exp\(-pdist\(a, \"sqeuclidean\"\) / M\[.bw2.\]\)")
        i_store, store_line = find_line(L.P10A, r"bw2=M\[.bw2.\]")
        i_redraw, redraw_line = find_line(L.P10A, r"ia = r_ref\.choice\(len\(Ds_std\)")
        ok = (takes and not reassigned and i_call and i_k and i_store and i_redraw)
        return (ok, f"mmd_rbf takes bw2 as a parameter and never reassigns it (reassignment found: {reassigned}); "
                    f"p6b:{i_call} = {call_site!r}; p10a:{i_k} (inside the K loop) = {k_line!r}; "
                    f"p10a:{i_store} stores {store_line!r}; the subset that IS redrawn per draw is the MMD "
                    f"reference: p10a:{i_redraw} = {redraw_line!r} (p6b:220-221 inside mmd_rbf), matching "
                    f"Methods sec.9 'K = 5 draws' - the SCALE is not redrawn")

    L.guard("BW-02b", bw02b,
            f"{DEF_MED} --'held fixed thereafter': the bandwidth scalar is produced once per task OUTSIDE "
            f"the K-draw loop and every consumer reads the cached scalar, so all K MMD draws share one scale",
            "live contrast (t07 FORK-04b, not re-tested): p6_02:230 recomputes the bandwidth INSIDE its K "
            "loop from pooled a+b and redraws it per draw - a different documented convention")

    # ------------------------------------------------------------------ BW-02c no outcome dependence
    def bw02c():
        code = compile(L.grab(*L.BLOCK_BW_P6B), "<bw-block>", "exec")
        names = sorted(set(code.co_names))
        forbidden = {"y", "yte", "ytr", "yv", "target", "outcome", "Dt", "Dte", "Dt_raw", "Dall",
                     "ps", "auroc", "delta", "harmful"}
        leak = sorted(set(names) & forbidden)
        got_a, _ = run_block(L.BLOCK_BW_P6B, DTR_BIG, sc_big,
                             extra=dict(y=np.zeros(len(DTR_BIG)), Dte=np.zeros((50, 7)), auroc=0.1))
        got_b, _ = run_block(L.BLOCK_BW_P6B, DTR_BIG, sc_big,
                             extra=dict(y=np.ones(len(DTR_BIG)), Dte=np.ones((50, 7)), auroc=0.9))
        bw_ok = (all(g.get("bw2") is not None and np.isfinite(g["bw2"]) for g in (got_a, got_b))
                 and got_a["bw2"] == got_b["bw2"])
        ib_same = None
        if got_a.get("ib") is not None and got_b.get("ib") is not None:
            ib_same = bool(np.array_equal(np.asarray(got_a["ib"]), np.asarray(got_b["ib"])))
        ok = (not leak and bw_ok and ib_same is not False)
        return (ok, f"compile-time names in scope = {names}; outcome/target names present = {leak or 'none'}; "
                    f"bw2 unchanged when y / Dte / auroc vary in the env: {bw_ok}; bound index vector "
                    f"unchanged: {'not observable (name not bound)' if ib_same is None else ib_same}")

    L.guard("BW-02c", bw02c,
            'Methods_v1.md sec.9: "the bandwidth was not scanned or selected against any outcome" -- '
            "compile-time scan of the verbatim block for outcome/target names, plus a behavioural probe "
            "that supplies two different y / Dte / auroc settings in the environment",
            "a bandwidth that read the outcome or the target sample would change every MMD value and "
            "leak target information into the kernel scale")

    # ------------------------------------------------------------------ BW-03 / 03b standardization source
    try:                                   # guarded: a failed fixture exec must not kill the family
        sc_pooled = sc_of(np.vstack([DTR_BIG, DTGT]))
        bw_pooled_scaler = float(np.median(pdist(sc_pooled.transform(DTR_BIG)[SEL_BIG], "sqeuclidean")))
        pooled_rows = np.vstack([DTR_BIG, DTGT])
        ib_pooled = declared_draw(len(pooled_rows))
        bw_pooled_rows = float(np.median(pdist(sc_big.transform(pooled_rows)[ib_pooled], "sqeuclidean")))
        bw_raw = float(np.median(pdist(DTR_BIG[SEL_BIG], "sqeuclidean")))
    except Exception:
        bw_pooled_scaler = bw_pooled_rows = bw_raw = float("nan")

    def bw03():
        if exec_failed(got_big):
            return False, ("the block produced no finite bw2 on the n=1200 fixture - see BW-00 for "
                           "the infrastructure reason")
        ref = float(np.median(pair_sqdist_chunked(Z_BIG[SEL_BIG])))
        d_pool = 100.0 * (bw_pooled_scaler - bw_prod) / bw_prod
        ok = med_rel(ref, bw_prod) <= TOL_REL and med_rel(bw_pooled_scaler, bw_prod) > 1e-3
        return (ok, f"prod={bw_prod!r} on source-TRAIN-standardized rows (independent ref={ref!r} over "
                    f"the consumed pair set: {SEL_HOW}); "
                    f"pooled-scaler counterfactual={bw_pooled_scaler!r} ({d_pool:+.1f}%); "
                    f"production is on the SOURCE-TRAIN-only side")

    L.guard("BW-03", bw03,
            f"{DEF_FIX} --the bandwidth is computed on data standardized by a scaler fitted on the "
            f"SOURCE-TRAIN rows only; the counterfactual fits the scaler on pooled source+target (the "
            f"documented defect) and must be clearly distinguishable, else the test cannot discriminate",
            "Methods sec.9 for the correct side: 'on source-training-standardized features'")

    def bw03b():
        d_rows = 100.0 * (bw_pooled_rows - bw_prod) / bw_prod
        d_raw = 100.0 * (bw_raw - bw_prod) / bw_prod
        ok = med_rel(bw_pooled_rows, bw_prod) > 1e-3 and med_rel(bw_raw, bw_prod) > 1e-3
        return (ok, f"pooled source+target ROWS with the source-train scaler={bw_pooled_rows!r} ({d_rows:+.1f}%); "
                    f"UNSTANDARDIZED raw Dtr={bw_raw!r} ({d_raw:+.1f}%); both halves of the documented "
                    f"defect differ from production ({bw_prod!r})")

    L.guard("BW-03b", bw03b,
            f"{DEF_FIX} --second half of the documented defect, split out: (i) distances taken on pooled "
            f"source+target rows and (ii) distances taken WITHOUT standardization must both differ from "
            f"production - so a partial repair (scaler fixed but rows still pooled, or the reverse) is visible",
            "decomposition of the pooled/standardization counterfactuals")

    def bw03c():
        hits = []
        for path, tag in ((L.P6B, "p6b"), (L.P10A, "p10a")):
            for i, ln in enumerate(read_lines(path), 1):
                if re.search(r"\bsc\s*=\s*StandardScaler\(\)\.fit\(", ln):
                    hits.append((tag, i, ln.strip()))
        # lenient by design (v1.3): every scaler fit FOUND must be fitted on Dtr (accepting the
        # `X=Dtr` spelling) and EACH file must contain at least one; the behavioural gates are
        # BW-03/03b, this is quoted supporting evidence, so >=1 per file replaces ==2 exactly.
        dtr_fit = {i: bool(re.search(r"\.fit\(\s*(?:X\s*=\s*)?Dtr\s*\)", h[2])) for i, h in enumerate(hits)}
        per_file = {t: sum(1 for h in hits if h[0] == t) for t in ("p6b", "p10a")}
        fits_ok = bool(hits) and all(dtr_fit.values()) and all(per_file[t] >= 1 for t in ("p6b", "p10a"))
        i6, l6 = find_line(L.P6B, r"Dtr,\s*Dva,\s*Dte\s*=\s*Dall\[i_tr\]")
        i10, l10 = find_line(L.P10A, r"Dtr,\s*Dva,\s*Dte\s*=\s*Dall\[i_tr\]")
        dtr_ok = i6 is not None and i10 is not None
        dtr_std = got_big.get("Dtr_std")
        std_txt = ("not observable (`Dtr_std` not bound)" if dtr_std is None
                   else str(bool(np.array_equal(dtr_std, sc_big.transform(DTR_BIG)))))
        return (fits_ok and dtr_ok,
                f"scaler fits found ({len(hits)}; per file {per_file}; every one fitted on Dtr: "
                f"{all(dtr_fit.values())}): {[(h[0], h[1], h[2]) for h in hits]}; Dtr is the train "
                f"split: p6b:{i6} = {l6!r}, p10a:{i10} = {l10!r}; the verbatim block applies that "
                f"scaler (Dtr_std == sc.transform(Dtr): {std_txt})")

    L.guard("BW-03c", bw03c,
            f"{DEF_MED} --static conformance of the caller context: every scaler fit FOUND in each file "
            f"(at least one per file; `X=Dtr` spelling accepted) is fitted on Dtr, the source-train "
            f"split, quoted by line; the verbatim block then transforms Dtr with that scaler",
            "supporting evidence for BW-03/03b, not a substitute for them; de-brittled in v1.3 "
            "(0 hits in a file, or a helper-moved fit, is reported by the quotes - the behavioural "
            "gates remain BW-03/03b)")

    # ------------------------------------------------------------------ BW-04 four-way discrimination
    def bw04():
        if exec_failed(got_big):
            return False, ("the block produced no finite bw2 on the n=1200 fixture - see BW-00 for "
                           "the infrastructure reason")
        d2 = pair_sqdist_chunked(Z_BIG[SEL_BIG])
        n_pairs = len(d2)
        med_sq = float(np.median(d2))
        mean_sq = float(d2.mean())
        med_eu = float(np.median(np.sqrt(d2)))
        mean_eu = float(np.sqrt(d2).mean())
        sqrt_med = float(np.sqrt(med_sq))
        gap = abs(med_eu - sqrt_med) / sqrt_med
        ds = pair_sqdist_chunked(sc_small.transform(DTR_SMALL))          # odd pair count
        exact_odd = abs(float(np.median(np.sqrt(ds))) - float(np.sqrt(np.median(ds))))
        ok = (med_rel(med_sq, bw_prod) <= TOL_REL
              and med_rel(mean_sq, med_sq) > 1e-6
              and med_rel(med_eu, med_sq) > 1e-6
              and exact_odd <= 1e-12)
        return (ok, f"production position = MEDIAN of SQUARED = {med_sq!r} (prod={bw_prod!r}); "
                    f"mean of squared={mean_sq!r} ({100 * (mean_sq / med_sq - 1):+.1f}%); "
                    f"median unsquared={med_eu!r} and mean unsquared={mean_eu!r} (sqrt scale: "
                    f"{100 * (med_eu / med_sq - 1):+.1f}% vs the squared median, NOT a like-for-like kernel "
                    f"argument); monotone relation median(sqrt(x)) vs sqrt(median(x)): {gap:.2e} relative on "
                    f"{n_pairs} pairs (even count, second-order), exactly 0 on the odd-count n={N_SMALL} "
                    f"fixture: {exact_odd:.1e}")

    L.guard("BW-04", bw04,
            f'{DEF_MED} --"median heuristic applied to source **training** pairwise squared distances": the '
            f"production value must sit at median-of-SQUARED and be distinguishable from mean-of-squared, "
            f"median-of-Euclidean and mean-of-Euclidean on the same fixture",
            "median-of-squared is the documented choice; its position is reported, not assumed")

    # ------------------------------------------------------------------ BW-05 degenerate guard fires
    D_IDENT = np.tile(np.array([[60.0, 1.5, 0.0, 1.0, 0.0, 1.0, 0.0]]), (300, 1))

    def bw05():
        fb = fallback_literal(L.BLOCK_BW_P6B)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            sc_i = sc_of(D_IDENT)
            got_i, _ = run_block(L.BLOCK_BW_P6B, D_IDENT, sc_i)
            rows_i = pair_sel(got_i, len(D_IDENT))[0]
            raw_i = float(np.median(pdist(sc_i.transform(D_IDENT)[rows_i], "sqeuclidean")))
            Dnan = source_fixture(400, 90004)
            Dnan[:300, :] = np.nan                     # 75% of rows NaN -> a majority of pairs non-finite
            sc_n = sc_of(Dnan)
            got_n, _ = run_block(L.BLOCK_BW_P6B, Dnan, sc_n)
            rows_n = pair_sel(got_n, len(Dnan))[0]
            d2_n = pdist(sc_n.transform(Dnan)[rows_n], "sqeuclidean")
            raw_n = float(np.median(d2_n))
            nanmed_n = float(np.nanmedian(d2_n))       # what a NaN-skipping median would have returned
        if exec_failed(got_i) or exec_failed(got_n):
            return False, (f"block produced no finite value on a degenerate fixture "
                           f"(identical rows -> {got_i.get('bw2')!r}; NaN majority -> {got_n.get('bw2')!r})")
        pos = lambda v: bool(v is not None and np.isfinite(v) and v > 0)
        if fb is None:                                 # no parseable guard line: fall back to the
            ok_i = raw_i == 0.0 and pos(got_i["bw2"])  # weaker finite-positive invariant
            ok_n = (not np.isfinite(raw_n)) and pos(got_n["bw2"])
        else:
            # the guard FIRED iff the raw statistic hit the degenerate branch (exactly 0 on the
            # identical-rows fixture; non-finite on the NaN-majority fixture) AND the returned value
            # is the block's own fallback literal - not some NaN-skipping statistic that never took
            # the guard path.  This is the v1.3 repair of the manufactured PASS where np.nanmedian
            # returned 13.66 while the block's own guard was never executed.
            ok_i = (raw_i == 0.0) and (got_i["bw2"] == fb)
            ok_n = (not np.isfinite(raw_n)) and (got_n["bw2"] == fb)
        return (bool(ok_i and ok_n),
                f"identical-rows fixture: raw median={raw_i!r} -> returned bw2={got_i['bw2']!r} "
                f"(guard fired: {ok_i}); NaN-majority fixture (300/400 rows NaN): raw median={raw_n!r} "
                f"(a NaN-skipping median would have returned {nanmed_n!r}) -> returned bw2={got_n['bw2']!r} "
                f"(guard fired: {ok_n}); block's own fallback literal parsed from its text: {fb!r}")

    L.guard("BW-05", bw05,
            f"{DEF_KER} --a bandwidth must remain finite and positive for any input, otherwise the kernel "
            f"is undefined or degenerates: the degenerate guard must FIRE on both the zero-distance branch "
            f"(all-identical rows, raw median exactly 0) and the non-finite branch (NaN-majority input, raw "
            f"np.median non-finite) - i.e. the RETURNED value must equal the block's own fallback literal "
            f"parsed from its text; a NaN-skipping median that never takes the guard path is a defect",
            "the fallback VALUE itself is a convention and is characterised separately in BW-05b; two "
            "independent triggers constructed; an infinity column cannot reach the guard because "
            "sklearn's scaler raises before the block runs")

    def bw05b():
        if exec_failed(got_big):
            return False, ("the block produced no finite bw2 on the n=1200 fixture - see BW-00 for "
                           "the infrastructure reason")
        src = L.grab(*L.BLOCK_BW_P6B)
        only_one = len(re.findall(r"(?<![\d.])1\.0(?![\d])", src)) == 1 and "bw2 = 1.0" in src
        d2 = pair_sqdist_chunked(Z_BIG[SEL_BIG])
        k_fb = np.exp(-d2 / 1.0)
        k_mh = np.exp(-d2 / bw_prod)
        return (only_one,
                f"fallback constant is exactly the literal 1.0 (the only 1.0 in the block: {only_one}); "
                f"median kernel entry {float(np.median(k_mh)):.6f} at the median-heuristic scale (exactly 1/e, "
                f"half the pairs have d2 <= bw2) vs {float(np.median(k_fb)):.3e} at the 1.0 fallback; "
                f"fraction of kernel entries > 1e-3: {float((k_mh > 1e-3).mean()):.4f} vs "
                f"{float((k_fb > 1e-3).mean()):.4f}; max fallback entry {float(k_fb.max()):.3e} - the "
                f"fallback makes k ~ 0 for essentially every pair, so MMD2 collapses towards 0 and would "
                f"read as 'no shift' rather than surfacing the degeneracy; the fallback is NOT scale-matched "
                f"(for d standardized features the typical squared distance is ~2d; here the median-heuristic "
                f"scale is {bw_prod:.2f}, i.e. {bw_prod:.1f}x the fallback)")

    guard_warn("BW-05b", bw05b,
               f'{DEF_MED} --characterisation (warn_only, NOT a correctness gate): the fallback VALUE is a '
               f'convention, not mathematics - the check records which value the code uses and quantifies '
               f'what it does to k(x,y)=exp(-norm(x-y)^2/bw2) under source-train standardization',
               "the fallback keeps the pipeline alive but is not a distance scale - reported for the author; "
               "BW-05 (hard) separately requires that the guard actually FIRED on both trigger fixtures")

    # ------------------------------------------------------------------ BW-06 p6b vs p10a parity
    def bw06():
        out, ok = [], True
        for tag, D, sc in (("n=1200", DTR_BIG, sc_big), ("n=450", DTR_SMALL, sc_small)):
            g6, _ = run_block(L.BLOCK_BW_P6B, D, sc)
            g10, _ = run_block(L.BLOCK_BW_P10A, D, sc)
            same = (g6["bw2"] is not None and g10["bw2"] is not None
                    and np.isfinite(g6["bw2"]) and np.isfinite(g10["bw2"]) and g6["bw2"] == g10["bw2"])
            ib_cmp = None
            if g6.get("ib") is not None and g10.get("ib") is not None:
                ib_cmp = bool(np.array_equal(np.asarray(g6["ib"]), np.asarray(g10["ib"])))
                same = same and ib_cmp
            ok = ok and same
            out.append(f"{tag}: p6b={g6['bw2']!r} p10a={g10['bw2']!r} identical={same} "
                       f"(index vectors identical: "
                       f"{'not observable (name not bound)' if ib_cmp is None else ib_cmp})")
        return (ok, "; ".join(out))

    L.guard("BW-06", bw06,
            "parity: the two live copies of a frozen quantity must agree EXACTLY (bit-for-bit) on a shared "
            "fixture - they are separate rewrites of the same rule in two files that both feed results; a "
            "bound index vector must agree too, an unbound one is reported as not observable",
            "t07 covers fork parity for the MMD kernel; the bandwidth copies are covered here")

    # ------------------------------------------------------------------ BW-07 subsample vs full-pair median (REAL frozen corpus)
    def bw07():
        files, models = scan_frozen_corpus()
        n_files, n_models = len(files), len(models)
        if n_models == 0:
            return False, f"no cached tasks in {CACHE_DIR} - characterisation cannot run"
        first = [rs[0] for rs in models.values()]
        allrec = [r for rs in models.values() for r in rs]
        gaps = np.array([r["gap"] for r in first], dtype=float)
        finite = bool(np.isfinite(gaps).all()
                      and all(np.isfinite(r["full"]) and np.isfinite(r["bw"]) for r in first))
        n_exact = sum(1 for r in first if r["method"] == "exact")
        n_samp = n_models - n_exact
        n_bind = sum(1 for r in allrec if r["n"] > MMD_M)
        n_min, n_max = min(r["n"] for r in allrec), max(r["n"] for r in allrec)
        n_dup = sum(len(rs) - 1 for rs in models.values())
        dup_disagree = sum(1 for rs in models.values() for r in rs[1:] if not r.get("dup"))
        live_rel = [r["rel_live"] for r in allrec]
        n_live_fail = sum(1 for v in live_rel if not np.isfinite(v))
        max_live = max((v for v in live_rel if np.isfinite(v)), default=float("nan"))

        def _load(r):
            with open(os.path.join(CACHE_DIR, r["fn"]), "rb") as fh:
                return np.asarray(pickle.load(fh)["Dtr_std"], dtype=float)

        # ---- controls: the exact full-pair median must equal an independent FULL enumeration on
        # small n before it is trusted on the big ones; the sampler is controlled on the same task.
        small = min(first, key=lambda r: r["pairs"])
        Zs = _load(small)
        ctrl_small = float(np.median(pdist(Zs, "sqeuclidean")))
        rel_small = abs(ctrl_small - small["full"]) / max(abs(ctrl_small), 1e-300)
        samp_ctrl = sample_pair_median(Zs, BW07_CTRL_SAMPLE_PAIRS, 4242)
        rel_samp = abs(samp_ctrl - ctrl_small) / max(abs(ctrl_small), 1e-300)
        mid = max((r for r in first if r["n"] <= BW07_BRUTE_MAX_ROWS), key=lambda r: r["n"], default=None)
        if mid is not None:
            ctrl_mid = float(np.median(pdist(_load(mid), "sqeuclidean")))
            rel_mid = abs(ctrl_mid - mid["full"]) / max(abs(ctrl_mid), 1e-300)
        else:
            rel_mid = float("nan")
        # ---- the worst model, characterised one step further: 30 fresh 1000-row draws under the
        # SAME convention - if the frozen draw sits at an extreme, the convention itself is spread.
        wr = max(first, key=lambda r: abs(r["gap"]))
        Zw = _load(wr)
        nw = len(Zw)
        draws = np.array([float(np.median(pdist(Zw[np.random.default_rng(777000 + i).choice(
            nw, min(MMD_M, nw), replace=False)], "sqeuclidean"))) for i in range(BW07_FRESH_DRAWS)])
        pct = 100.0 * float((draws < wr["bw"]).mean())

        ok_meas = (finite and n_files > 0 and n_live_fail == 0 and max_live <= 1e-9
                   and dup_disagree == 0
                   and rel_small <= 1e-9 and (mid is None or rel_mid <= 1e-9) and rel_samp <= 1e-2)
        material = finite and float(np.abs(gaps).max()) > BW07_NEGLIGIBLE_PCT
        ok = bool(ok_meas and not material)

        method_s = (f"EXACT (unique-row multiplicity expansion) for {n_exact}/{n_models}"
                    + (f", sampled ({BW07_SAMPLE_PAIRS / 1e6:.0f}M unbiased pairs per model) for {n_samp}"
                       if n_samp else " - the sampled fallback was not needed here"))
        ctrl_s = (f"brute-force control vs the exact median: n={small['n']} rel={rel_small:.1e}"
                  + (f", n={mid['n']} rel={rel_mid:.1e}" if mid is not None else "")
                  + f"; sampler control rel={rel_samp:.1e}")
        top = sorted(first, key=lambda r: -abs(r["gap"]))[:3]
        top_s = "; ".join(f"{r['key'][0]}/{r['key'][1]}/{r['key'][2]} {r['gap']:+.2f}% "
                          f"(frozen {r['bw']:.4f} vs full {r['full']:.4f})" for r in top)
        live_s = (f"frozen bw2 reproduced through the LIVE block on the cached Dtr_std: max rel {max_live:.1e}"
                  + (f" - live reproduction FAILED on {n_live_fail}/{n_files} tasks (the live block's own "
                     f"cost/behaviour, reported not hidden)" if n_live_fail else ""))
        observed = (
            f"frozen p10_cache scanned READ-ONLY: {n_files} tasks / {n_models} unique source models "
            f"(duplicate tasks matched on the fingerprint (n, Z.sum(), Z[0,0], bw2) and NOT re-executed: "
            f"{n_dup - dup_disagree}/{n_dup} identical; n_train {n_min}-{n_max}; {n_bind}/{n_files} tasks "
            f"with n_train > MMD_M={MMD_M}, so the fixed {MMD_M}-row subsample binds for every reported "
            f"bandwidth); {live_s}; full-pair median {method_s}; {ctrl_s}; "
            f"gap_pct = 100*(frozen - full)/full over the {n_models} models: min {np.min(gaps):+.2f} "
            f"p25 {np.percentile(gaps, 25):+.2f} median {np.median(gaps):+.2f} p75 {np.percentile(gaps, 75):+.2f} "
            f"max {np.max(gaps):+.2f}; abs(gap)>10pct: {int((np.abs(gaps) > 10).sum())}/{n_models}, "
            f"abs(gap)>25pct: {int((np.abs(gaps) > 25).sum())}/{n_models}; largest: {top_s}; worst model "
            f"fresh {BW07_FRESH_DRAWS} draws span {draws.min():.2f}-{draws.max():.2f} "
            f"(median {np.median(draws):.2f}), frozen draw at percentile {pct:.0f}"
        )
        return ok, observed

    guard_warn("BW-07", bw07,
               f'{DEF_MED} --characterisation (warn_only, NOT a correctness gate) on the REAL frozen corpus: '
               f'the code uses min(MMD_M, n_train) rows, Methods sec.9 documents the median heuristic on the '
               f'source-training pairwise squared distances and does NOT document a bandwidth subsample; '
               f'the reported quantity is the frozen bw2 against the EXACT full-pair median of the same '
               f'source-train matrix (gap_pct), with abs(gap) <= {BW07_NEGLIGIBLE_PCT}% declared negligible',
               "WARN-class characterisation only: the subsample of min(MMD_M, n_train) rows is DELIBERATELY "
               "KEPT by the correction spec (CS-06) and is undisclosed in Methods sec.9 (which documents "
               "m = 1,000 / K = 5 only for the MMD reference), so the deviation is reported, never failed. "
               "The WARN is NUMERIC: the 5% negligible criterion is a fixed number, so documenting/disclosing "
               "the subsample cannot clear it - only changing the estimator (e.g. to the full-pair median, "
               "gap -> 0) or re-defining the criterion would. On the frozen corpus the measured gap is NOT "
               "numerically negligible (up to +52.5%): advisory for the author, not a to-do that flips PASS")

    print("# t09_bandwidth.py finished", flush=True)


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        L.emit("BW-00", "FAIL", "every BW test must be guarded by L.guard",
               f"EXCEPTION {type(e).__name__}: {e}",
               traceback.format_exc().splitlines()[-1])
    sys.exit(0)
