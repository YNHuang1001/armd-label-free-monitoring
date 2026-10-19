# -*- coding: utf-8 -*-
r"""t01_mmd.py — family MMD | estimator validation suite (red-light baseline)

WHAT IS TESTED
  The production kernel-MMD estimator `mmd_rbf` (p6b_01_corrected_pipeline.py:216-228) and its two
  split re-implementations (p10a_prepare_cache.py:273-280 source within-term `Kaa`;
  p10b_subsample_curves.py:85-93 target within-term + cross-term + combination), loaded VERBATIM
  through evs_load (L.load_def / L.exec_block).  No production file is edited.

MATHEMATICAL DEFINITION IMPLEMENTED BY THE REFERENCE (written FROM this text, never adapted from
the production code):

  Unbiased U-statistic maximum mean discrepancy -- Gretton, Borgwardt, Rasch, Schoelkopf & Smola,
  "A Kernel Two-Sample Test", JMLR 13 (2012) 723-773, eq. (3):

      MMD_u^2(A,B) = 1/(n(n-1)) * sum_{i != j} k(a_i,a_j)
                   + 1/(m(m-1)) * sum_{i != j} k(b_i,b_j)
                   - 2/(n*m)   * sum_{i,j}   k(a_i,b_j)

  with the Gaussian/RBF kernel k(x,y) = exp(-||x-y||^2 / bw2), k(x,x)=1, built on the FULL n x n and
  n x m matrices by an explicit O(n*m) row-chunked double loop (numpy broadcasting; scipy is NOT
  used by the reference) plus, for the small-sample case, a pure-Python scalar double loop.
  bw2 = median heuristic = median of the squared Euclidean distances over UNORDERED source-train
  pairs (p6b header lines 9-11 and the task spec), computed here from the full squared-distance
  matrix, not from scipy.pdist.

DERIVED IDENTITY for A == B (used in MMD-01; derivation from the definition above):
  Write C = (1/n^2) * sum_{i,j} k(a_i,a_j) (FULL matrix, diagonal INCLUDED, k(x,x)=1).  Because
  sum_{i != j} k = n^2 C - n,
        A_ref = B_ref = (n^2 C - n) / (n(n-1)) = (nC - 1)/(n-1),
  and for A == B the cross matrix is the same n x n matrix, so C_ref = C and
        MMD_u^2(A,A) = 2(nC-1)/(n-1) - 2C = 2(C - 1)/(n - 1)      (exact, deterministic).
  NOTE: the naive claim "A == B => MMD^2 = 0" is FALSE for the realised statistic; only the
  expectation over repeated independent draws from one distribution is 0 (Gretton 2012 sec. 5).
  The production estimator divides the pdist UPPER TRIANGLE (n(n-1)/2 unordered pairs, whose sum is
  exactly half of sum_{i != j}) by n(n-1), so both within terms come out at 0.5 x A_ref and
        MMD_prod(A,A) = A_ref - 2C = ((2-n)C - 1)/(n-1)   -> too negative by exactly A_ref.

CONFIRMED DEFECT UNDER TEST
  (a) pdist returns only the upper triangle while the normalisation n(n-1) counts ordered pairs =>
      both within-sample kernel terms are HALF of the unbiased U-statistic.

SCOPE NOTE (tests that are NOT defect probes, and why): MMD-02c (cross term), MMD-06 (positive
  under separation) and the *R reference lower-bound controls are RED-LIGHT CONTROLS: they pass on
  the current code by construction because the pdist defect does not touch the cross term and is
  symmetric in A and B.  They are included to prove the harness discriminates, and they guard
  against a future fix that overshoots.  MMD-05 keeps its swap-symmetry control role but is ALSO
  a value probe since the equal-size-blind-spot hardening (see next paragraph).

EQUAL-SIZE BLIND SPOT (closed by the hardening in this revision): every exact-value comparison
  used to be built on n == m fixtures (1000/1000 or 40/40), while the one unequal-size test
  (MMD-05, n=40 m=60) asserted only A<->B symmetry -- and a symmetric halving of both within terms
  is invisible to a symmetry check.  An adversarial verifier showed that a WRONG patch applying
  the factor 2 ONLY when len(a) == len(b) leaves the pre-hardening suite with ZERO FAILs (the
  equal-size fixtures all get the factor 2; the unequal-size check only sees a symmetric error).
  MMD-09 / MMD-10 close this by brute-forcing the exact U-statistic VALUE on small unequal pairs
  (n=30,m=50 and the reverse n=50,m=30; both sizes < MMD_M so every row is consumed, the K draws
  are the same permutation-invariant statistic and the production value is deterministic, not
  Monte Carlo), and MMD-05 now carries the same value identity at n=40, m=60.

OUTPUT CONTRACT: RESULT|test_id|PASS/FAIL/WARN|expected (source)|observed|note
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
from scipy.spatial.distance import pdist, cdist

import evs_load as L

FAM = "MMD"
SEED_TASK = 20261005          # p6b:30 / p10a:32 / p10b:27
MMD_M, MMD_K = 1000, 5        # p6b:31 / p10a:33 / p10b:28
TOL = 1e-8                    # float-summation tolerance for O(n^2) sums of O(1) terms

DEF_SRC = "Gretton et al. 2012 JMLR 13:723-773 eq.(3) unbiased U-statistic, k(x,y)=exp(-||x-y||^2/bw2)"
BOUND_SRC = "PSD kernel with k(x,x)=1 => unbiased MMD_u^2 >= -1/(n-1) - 1/(m-1)"

L.banner(f"{FAM} — t01_mmd.py | production mmd_rbf (p6b:216-228) + p10a:273-280 + p10b:85-93")


# ============================================================ reference (from the definition only;
# scipy is used by the production's own calls, never by the reference) =====================
def sqdist_full(A, B, chunk=200):
    """Full |A| x |B| matrix of squared Euclidean distances, explicit row-chunked double loop."""
    out = np.empty((len(A), len(B)), dtype=float)
    for i0 in range(0, len(A), chunk):
        a = A[i0:i0 + chunk]
        out[i0:i0 + chunk] = ((a[:, None, :] - B[None, :, :]) ** 2).sum(-1)
    return out


def kernel_full(A, B, bw2, chunk=200):
    """Full |A| x |B| Gaussian kernel matrix k(x,y)=exp(-||x-y||^2/bw2); diagonal is 1 when A is B."""
    return np.exp(-sqdist_full(A, B, chunk) / bw2)


def ref_terms(A, B, bw2):
    """Gretton 2012 eq.(3) terms from the definition, on FULL kernel matrices."""
    n, m = len(A), len(B)
    Kaa = kernel_full(A, A, bw2)
    Kbb = kernel_full(B, B, bw2)
    Kab = kernel_full(A, B, bw2)
    np.fill_diagonal(Kaa, 0.0)          # i != j
    np.fill_diagonal(Kbb, 0.0)          # i != j
    a_within = Kaa.sum() / (n * (n - 1))
    b_within = Kbb.sum() / (m * (m - 1))
    cross = Kab.sum() / (n * m)
    return dict(a_within=float(a_within), b_within=float(b_within), cross=float(cross),
                mmd2=float(a_within + b_within - 2.0 * cross))


def ref_terms_scalar(A, B, bw2):
    """Same definition, pure-Python scalar double loop (independent of numpy matrix tricks)."""
    n, m = len(A), len(B)
    sa = 0.0
    for i in range(n):
        for j in range(n):
            if i != j:
                sa += float(np.exp(-np.sum((A[i] - A[j]) ** 2) / bw2))
    sb = 0.0
    for i in range(m):
        for j in range(m):
            if i != j:
                sb += float(np.exp(-np.sum((B[i] - B[j]) ** 2) / bw2))
    sc = 0.0
    for i in range(n):
        for j in range(m):
            sc += float(np.exp(-np.sum((A[i] - B[j]) ** 2) / bw2))
    a_within, b_within, cross = sa / (n * (n - 1)), sb / (m * (m - 1)), sc / (n * m)
    return dict(a_within=a_within, b_within=b_within, cross=cross,
                mmd2=a_within + b_within - 2.0 * cross)


def bw2_median(X, chunk=200):
    """Median heuristic: median of squared distances over unordered pairs (p6b header lines 9-11)."""
    D2 = sqdist_full(X, X, chunk)
    iu = np.triu_indices(len(X), 1)
    return float(np.median(D2[iu]))


# ------------------------------------------------------- production loading (verbatim, instrumented)
def instrumented():
    """pdist/cdist wrappers that RECORD what the production passes them; returned arrays are the
    genuine scipy outputs, so the production's arithmetic is bit-for-bit unaffected."""
    calls = []

    def _pdist(a, metric="euclidean", **kw):
        v = pdist(a, metric, **kw)
        calls.append(("pdist", a, None, v))
        return v

    def _cdist(a, b, metric="euclidean", **kw):
        v = cdist(a, b, metric, **kw)
        calls.append(("cdist", a, b, v))
        return v

    return calls, _pdist, _cdist


def load_prod_p6b():
    calls, _p, _c = instrumented()
    ns = {"MMD_M": MMD_M, "MMD_K": MMD_K, "pdist": _p, "cdist": _c, "np": np}
    return L.load_def(L.P6B, "mmd_rbf", extra_ns=ns), calls


def parse_p6b_calls(calls, bw2):
    """Group the recorded (pdist a, pdist b, cdist a-b) triples into per-replicate term records.
    The production writes np.exp(-pdist(...)/bw2), so the recorded raw distances are exponentiated
    here with the SAME bw2 in order to reconstruct its kernel terms."""
    reps = []
    if len(calls) % 3 != 0:
        raise ValueError(f"pdist/cdist call count {len(calls)} is not a multiple of 3")
    for i in range(0, len(calls), 3):
        t0, a0, _, v0 = calls[i]
        t1, a1, _, v1 = calls[i + 1]
        t2, a2, b2, v2 = calls[i + 2]
        if (t0, t1, t2) != ("pdist", "pdist", "cdist"):
            raise ValueError(f"unexpected call pattern {(t0, t1, t2)}")
        n, m = len(a0), len(a1)
        reps.append(dict(n=n, m=m,
                         a_within=float(np.exp(-v0 / bw2).sum() / (n * (n - 1))),
                         b_within=float(np.exp(-v1 / bw2).sum() / (m * (m - 1))),
                         cross=float(np.exp(-v2 / bw2).mean())))
    return reps


def make_block(rng, n, shift=0.0, d=8, binary=(4, 5, 6, 7)):
    """Raw feature block: d continuous N(shift,1) + 4 Bernoulli(0.3) columns."""
    Z = rng.normal(shift, 1.0, size=(n, d))
    for j in binary:
        Z[:, j] = (rng.random(n) < 0.3).astype(float)
    return Z


def scaler_of(X):
    sd = X.std(0)
    sd[sd == 0] = 1.0
    return X.mean(0), sd


def standardize(X, mu, sd):
    return (X - mu) / sd


# ======================================================================================== tests
def _main():
    # ---- shared 1000 x 1000 identical-sample fixture (A == B exactly) -------------------------
    rng = np.random.default_rng(SEED_TASK)
    Z = make_block(rng, 1000)
    mu, sd = scaler_of(Z)
    A = standardize(Z, mu, sd)                 # mimics StandardScaler fitted on source train
    B = A.copy()                               # exactly the same rows
    bw2 = bw2_median(A)                        # median heuristic on source-train pairs
    n = len(A)
    Kfull = kernel_full(A, A, bw2)             # full matrix, diagonal INCLUDED (k(x,x)=1)
    C = float(Kfull.mean())
    A_ref = (n * C - 1.0) / (n - 1)            # exactly sum_{i != j} K / (n(n-1))
    ref_aa = 2.0 * (C - 1.0) / (n - 1)         # derived closed form for A == B

    # ------------------------------------------------------------------------------- MMD-01
    def mmd01():
        fn, calls = load_prod_p6b()
        prod, prod_sd = fn(A, B, bw2, np.random.default_rng(SEED_TASK + 2))
        rt = ref_terms(A, A, bw2)              # explicit full-matrix control of the closed form
        L.check("MMD-01b", abs(ref_aa - rt["mmd2"]) <= 1e-9,
                "t01_mmd.py docstring derivation for A==B: MMD_u^2 = 2(C-1)/(n-1); explicit full-matrix U-statistic is the control",
                f"closed_form={ref_aa:.10f} explicit_full_matrix={rt['mmd2']:.10f} delta={ref_aa-rt['mmd2']:.2e}",
                "control on the reference: closed form vs explicit O(n^2) sum")
        return (abs(prod - ref_aa) <= TOL,
                f"prod={prod:.8f} ref={ref_aa:.8f} absdiff={abs(prod-ref_aa):.3e} tol={TOL:g} C={C:.6f} A_ref={A_ref:.6f} prod_sd={prod_sd:.2e}")

    L.guard("MMD-01", mmd01,
            f"{DEF_SRC} | A==B closed form 2(C-1)/(n-1) derived in the t01 docstring (n=m=1000)",
            "identical samples: production is too negative by exactly A_ref because both within terms are halved")

    # ------------------------------------------------------------------------------- MMD-02
    def mmd02():
        rng2 = np.random.default_rng(SEED_TASK + 7)
        Za, Zb = make_block(rng2, 40), make_block(rng2, 40, shift=0.4)
        mu2, sd2 = scaler_of(Za)               # scaler fit on SOURCE only, as production does
        A2, B2 = standardize(Za, mu2, sd2), standardize(Zb, mu2, sd2)
        bw22 = bw2_median(A2)
        fn, calls = load_prod_p6b()
        prod, _ = fn(A2, B2, bw22, np.random.default_rng(SEED_TASK + 8))
        reps = parse_p6b_calls(calls, bw22)
        r0 = reps[0]
        ref, ref_sc = ref_terms(A2, B2, bw22), ref_terms_scalar(A2, B2, bw22)
        L.check("MMD-02c", abs(ref["cross"] - r0["cross"]) <= 1e-10,
                f"{DEF_SRC} | cross term = (1/(n*m)) sum over all i,j = plain mean of the FULL n x m kernel",
                f"prod_cross={r0['cross']:.10f} ref_cross={ref['cross']:.10f} scalar_double_loop={ref_sc['cross']:.10f}",
                "control: the cross term uses cdist (full matrix) and is NOT affected by the pdist defect")
        recon = float(np.mean([(r["a_within"] + r["b_within"] - 2 * r["cross"]) for r in reps]))
        L.check("MMD-02b", abs(recon - prod) <= 1e-12,
                "characterisation / current-state pin: the test-side capture with the sum/(n(n-1)) normalisation mirrors the production return value ONLY while production keeps that same unordered-pair normalisation; reads FAIL once the production within-term normalisation is repaired",
                f"recon={recon:.12f} prod={prod:.12f} delta={abs(recon-prod):.2e}",
                "diagnostic only, never a correctness check", warn_only=True)
        L.check("MMD-02d", abs(ref_sc["mmd2"] - ref["mmd2"]) <= 1e-9,
                "reference self-check: row-chunked full matrix vs pure-Python scalar double loop",
                f"broadcast={ref['mmd2']:.10f} scalar_double_loop={ref_sc['mmd2']:.10f}",
                "control on the reference itself")
        spread = max(abs(r["cross"] - r0["cross"]) for r in reps) + max(
            abs((r["a_within"] + r["b_within"]) - (r0["a_within"] + r0["b_within"])) for r in reps)
        assert spread <= 1e-12, f"K={len(reps)} reps are not identical (spread={spread:.3e}); the inversion assumes n=m<=MMD_M so every draw is a full-sample permutation"
        w_prod = prod + 2.0 * r0["cross"]     # production's own two within terms, recovered by inverting its return value
        ratio = w_prod / (ref["a_within"] + ref["b_within"])
        ra, rb = r0["a_within"] / ref["a_within"], r0["b_within"] / ref["b_within"]   # test-side reconstruction, pinned to the current normalisation (diagnostic only)
        return (abs(ratio - 1.0) <= 1e-9,
                f"ratio={ratio:.10f} W_prod={w_prod:.10f} ref_within={ref['a_within']+ref['b_within']:.10f}; test-side reconstruction (pinned to the current normalisation): ratio_A={ra:.10f} ratio_B={rb:.10f} prod_A={r0['a_within']:.6f} ref_A={ref['a_within']:.6f} prod_B={r0['b_within']:.6f} ref_B={ref['b_within']:.6f}")

    L.guard("MMD-02", mmd02,
            f"{DEF_SRC} | the production's own two within terms, recovered by inverting its return value (W_prod = mmd_return + 2*cross), must equal 1/(n(n-1)) * sum over the ORDERED pairs i != j + 1/(m(m-1)) * sum over the ORDERED pairs i != j -> ratio = 1.0",
            "n=m=40<=MMD_M: every K draw is the full sample permuted, so the inversion is exact; ratio 0.5 is the money shot (both production within terms exactly halved)")

    # ---------------------------------------------------------------------- MMD-03 / MMD-04
    def bound_case(tid, A3, B3, bw23, label):
        fn, _ = load_prod_p6b()
        prod, _ = fn(A3, B3, bw23, np.random.default_rng(SEED_TASK + 31))
        ref = ref_terms(A3, B3, bw23)
        n3, m3 = len(A3), len(B3)
        bound = -1.0 / (n3 - 1) - 1.0 / (m3 - 1)
        L.check(tid + "R", ref["mmd2"] >= bound - 1e-12, BOUND_SRC,
                f"ref={ref['mmd2']:.8f} bound={bound:.8f} margin={ref['mmd2']-bound:.3e}",
                f"reference adheres to the bound ({label})")
        return (prod >= bound - 1e-12,
                f"prod={prod:.8f} bound={bound:.8f} violation={prod-bound:.3e} ref={ref['mmd2']:.8f} ({label})")

    def mmd03():
        rng3 = np.random.default_rng(SEED_TASK + 17)
        Za, Zb = make_block(rng3, 40), make_block(rng3, 40, shift=0.4)
        mu3, sd3 = scaler_of(Za)
        A3, B3 = standardize(Za, mu3, sd3), standardize(Zb, mu3, sd3)
        return bound_case("MMD-03", A3, B3, bw2_median(A3), "small sample n=m=40")

    L.guard("MMD-03", mmd03, BOUND_SRC + " | production estimator, small sample",
            "the halved within terms push the production below the theoretical floor")

    L.guard("MMD-04", lambda: bound_case("MMD-04", A, B, bw2, "realistic 1000 x 1000, A==B"),
            BOUND_SRC + " | production estimator, realistic 1000x1000",
            "A==B is the adversarial worst case for the bound")

    # ------------------------------------------------------------------------------- MMD-05
    def mmd05():
        rng5 = np.random.default_rng(SEED_TASK + 23)
        Za, Zb = make_block(rng5, 40), make_block(rng5, 60, shift=0.3)
        mu5, sd5 = scaler_of(Za)
        A5, B5 = standardize(Za, mu5, sd5), standardize(Zb, mu5, sd5)
        bw25 = bw2_median(np.vstack([A5, B5]))
        fn, _ = load_prod_p6b()
        ab, _ = fn(A5, B5, bw25, np.random.default_rng(SEED_TASK + 24))
        ba, _ = fn(B5, A5, bw25, np.random.default_rng(SEED_TASK + 25))
        rab, rba = ref_terms(A5, B5, bw25), ref_terms(B5, A5, bw25)
        sym = abs(ab - ba) <= 1e-10 and abs(rab["mmd2"] - rba["mmd2"]) <= 1e-10
        val = abs(ab - rab["mmd2"]) <= 1e-10 and abs(ba - rba["mmd2"]) <= 1e-10   # added: exact VALUE identity
        ok = sym and val
        return (ok,
                f"symmetry: prod_AB={ab:.10f} prod_BA={ba:.10f} absdiff={abs(ab-ba):.2e} ref_AB={rab['mmd2']:.10f} ref_BA={rba['mmd2']:.10f}; value identity (n=40,m=60 both fully consumed): d_AB={abs(ab-rab['mmd2']):.3e} d_BA={abs(ba-rba['mmd2']):.3e} tol=1e-10")

    L.guard("MMD-05", mmd05,
            f"{DEF_SRC} | n=40, m=60, both sizes <= MMD_M so every row is consumed and all K draws are the same statistic: the production value must EQUAL the full-matrix U-statistic to 1e-10 in BOTH orientations AND be A<->B invariant",
            "unequal sizes; the value identity is what kills a wrong 'factor 2 only when len(a)==len(b)' conditional patch (a symmetric halving was invisible to the symmetry check alone)")

    # ------------------------------------------------------------------------------- MMD-06
    def mmd06():
        rng6 = np.random.default_rng(SEED_TASK + 41)
        Za, Zb = make_block(rng6, 1000), make_block(rng6, 1000, shift=5.0)
        mu6, sd6 = scaler_of(Za)
        A6, B6 = standardize(Za, mu6, sd6), standardize(Zb, mu6, sd6)
        bw26 = bw2_median(A6)
        fn, _ = load_prod_p6b()
        prod, _ = fn(A6, B6, bw26, np.random.default_rng(SEED_TASK + 42))
        ref = ref_terms(A6, B6, bw26)
        ratio = prod / ref["mmd2"]
        L.check("MMD-06b", abs(ratio - 1.0) <= 0.02,
                "with within terms ~1 and cross ~0 the production estimator must equal the U-statistic -> prod/ref = 1.0",
                f"prod/ref={ratio:.6f} prod={prod:.6f} ref={ref['mmd2']:.6f}",
                "diagnostic only: ~0.5 is the halved-within-term signature", warn_only=True)
        return ((prod > 0) and (ref["mmd2"] > 0),
                f"prod={prod:.6f} ref={ref['mmd2']:.6f} cross={ref['cross']:.3e} cluster_shift=5.0")

    L.guard("MMD-06", mmd06,
            "disjoint well-separated clusters => MMD^2 clearly > 0 (definition: within ~1, cross ~0)",
            "sign sanity / red-light control: passes by construction, it is not a defect probe")

    # ---------------------------------------------------------------------- MMD-07 / MMD-08
    def mmd07():
        calls, _p, _c = instrumented()
        env = dict(np=np, MMD_K=MMD_K, MMD_M=MMD_M, pdist=_p, cdist=_c,
                   Ds_std=A, M={"bw2": bw2}, r_ref=np.random.default_rng(SEED_TASK + 51))
        got, _src = L.exec_block((L.P10A, 291, 298), env, names=["Kaa", "mmd_src"])
        Kaa = got["Kaa"]
        if Kaa is None or len(Kaa) != MMD_K:
            raise ValueError("p10a block did not produce Kaa")
        ref = ref_terms(A, A, bw2)
        ra = Kaa[0] / ref["a_within"]
        return (abs(ra - 1.0) <= 1e-9,
                f"ratio={ra:.10f} Kaa[0]={Kaa[0]:.6f} ref_A={ref['a_within']:.6f} n={len(A)}")

    L.guard("MMD-07", mmd07,
            f"{DEF_SRC} | p10a source within term = 1/(n(n-1)) * sum over ORDERED pairs i != j",
            "verbatim block p10a_prepare_cache.py:273-280; same factor-2 defect as p6b, so a fix must touch this file too")

    def mmd08():
        calls, _p, _c = instrumented()
        src_env = dict(np=np, MMD_K=MMD_K, MMD_M=MMD_M, pdist=_p, cdist=_c,
                       Ds_std=A, M={"bw2": bw2}, r_ref=np.random.default_rng(SEED_TASK + 52))
        got, _ = L.exec_block((L.P10A, 291, 298), src_env, names=["Kaa", "mmd_src"])
        n_before = len(calls)
        rec = {}
        tgt_env = dict(np=np, MMD_K=MMD_K, MMD_M=MMD_M, pdist=_p, cdist=_c,
                       Dt_std_s=B, n=len(B), rng=np.random.default_rng(SEED_TASK + 53),
                       A={"Kaa": got["Kaa"], "mmd_src": got["mmd_src"], "bw2": bw2}, rec=rec)
        got2, _ = L.exec_block((L.P10B, 106, 114), tgt_env, names=["rec"])
        rec = got2["rec"]
        tail = calls[n_before:]
        if len(tail) != 2 * MMD_K:
            raise ValueError(f"expected {2*MMD_K} pdist/cdist calls from the p10b block, got {len(tail)}")
        m_ = len(B)
        Kbb0 = np.exp(-tail[0][3] / bw2)
        Kab0 = np.exp(-tail[1][3] / bw2)
        ref = ref_terms(A, B, bw2)
        ref_within = ref["a_within"] + ref["b_within"]
        vals_rep = [got["Kaa"][k] + float(np.exp(-tail[2 * k][3] / bw2).sum()) / (m_ * (m_ - 1))
                    - 2.0 * float(np.exp(-tail[2 * k + 1][3] / bw2).mean()) for k in range(MMD_K)]
        spread = max(vals_rep) - min(vals_rep)
        assert spread <= 1e-12, f"K={MMD_K} reps are not identical (spread={spread:.3e}); the inversion assumes m<=MMD_M so every draw is a full-sample permutation"
        w_prod = rec["mmd"] + 2.0 * float(Kab0.mean())    # p10b's two within terms, recovered by inverting its return value
        rw = w_prod / ref_within
        old_ts = got["Kaa"][0] + float(Kbb0.sum() / (m_ * (m_ - 1)))   # test-side reconstruction, pinned to the current normalisation (diagnostic only)
        L.check("MMD-08b", abs(rec["mmd"] - ref["mmd2"]) <= TOL,
                f"{DEF_SRC} | combined value of the p10a+p10b split implementation",
                f"p10b_mmd={rec['mmd']:.8f} ref={ref['mmd2']:.8f} absdiff={abs(rec['mmd']-ref['mmd2']):.3e} prod_cross={float(Kab0.mean()):.6f}",
                "cross term is correct, so the whole deviation comes from the two halved within terms")
        return (abs(rw - 1.0) <= 1e-9,
                f"ratio_within={rw:.10f} W_prod={w_prod:.6f} ref_within={ref_within:.6f} Kaa[0]={got['Kaa'][0]:.6f}; test-side reconstruction (pinned to the current normalisation): ratio={old_ts/ref_within:.10f}")

    L.guard("MMD-08", mmd08,
            f"{DEF_SRC} | both within terms, recovered by inverting the p10b return value (W_prod = mmd + 2*cross), are ordered-pair means over n(n-1) and m(m-1)",
            "verbatim chain p10a:273-280 + p10b:85-93, A==B 1000x1000; per-file copy 3/3")

    # ---------------------------------------------------------------------- MMD-09 / MMD-10
    def exact_unequal(tid, na, nb, seed_off):
        """Exact-value probe on a SMALL unequal pair: na, nb < MMD_M, so every row is consumed,
        the K draws are the same permutation-invariant statistic and the production return value
        is exact -- compared against the full-matrix U-statistic at 1e-10."""
        rng9 = np.random.default_rng(SEED_TASK + seed_off)
        Za, Zb = make_block(rng9, na), make_block(rng9, nb, shift=0.3)
        mu9, sd9 = scaler_of(Za)
        A9, B9 = standardize(Za, mu9, sd9), standardize(Zb, mu9, sd9)
        bw29 = bw2_median(A9)
        fn, _ = load_prod_p6b()
        prod, prod_sd = fn(A9, B9, bw29, np.random.default_rng(SEED_TASK + seed_off + 1))
        ref = ref_terms(A9, B9, bw29)                     # definition, full n x n / m x m / n x m matrices
        ref_sc = ref_terms_scalar(A9, B9, bw29)           # independent pure-Python control on the reference
        L.check(tid + "b", abs(ref_sc["mmd2"] - ref["mmd2"]) <= 1e-9,
                "reference self-check on the small unequal-size fixture: row-chunked full matrix vs pure-Python scalar double loop",
                f"broadcast={ref['mmd2']:.12f} scalar_double_loop={ref_sc['mmd2']:.12f}",
                "control on the reference itself")
        d = abs(prod - ref["mmd2"])
        return (d <= 1e-10,
                f"prod={prod:.12f} ref={ref['mmd2']:.12f} absdiff={d:.3e} tol=1e-10; n={na} m={nb} prod_sd={prod_sd:.2e} "
                f"ref_within={ref['a_within']+ref['b_within']:.6f} ref_cross={ref['cross']:.6f}")

    L.guard("MMD-09", lambda: exact_unequal("MMD-09", 30, 50, 61),
            f"{DEF_SRC} | exact value identity at n != m: n=30, m=50, both < MMD_M=1000 so every row is consumed and the production value is exact (no Monte Carlo) -- production return must equal the full-matrix U-statistic to 1e-10",
            "closes the unequal-sample-size blind spot: a wrong 'factor 2 only when len(a)==len(b)' conditional patch leaves both within terms halved here and is caught")

    L.guard("MMD-10", lambda: exact_unequal("MMD-10", 50, 30, 62),
            f"{DEF_SRC} | exact value identity at n != m, REVERSE orientation: n=50, m=30 -- same identity, neither side privileged",
            "reverse of MMD-09: a one-sided or equal-size-conditional repair cannot pass both orientations")

    print("# t01_mmd.py finished", flush=True)


if __name__ == "__main__":
    try:
        _main()
    except Exception as e:                     # never exit non-zero: the orchestrator aggregates
        import traceback
        L.emit("MMD-00", "FAIL", "every MMD test must be guarded by L.guard",
               f"EXCEPTION {type(e).__name__}: {e}",
               traceback.format_exc().splitlines()[-1])
    sys.exit(0)
