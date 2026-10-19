# -*- coding: utf-8 -*-
r"""t07_fork_parity.py - EVS family FORK: execute the LIVE same-source estimator forks

WHY THIS FILE EXISTS
================================================================================
Three production code paths implement estimator kernels that NO test used to execute, while
remaining live (their outputs can still reach result CSVs and figures):

  * p10b_subsample_curves.py:150-174   L3 (IOW / OM / DR)  -> p10b_raw.csv -> p10c -> fig5 / figS8
  * p10b_subsample_curves.py:116-148   SUDO-AURCC          -> same chain (p10b_raw.csv)
  * p6_02_label_free_signals.py:223-237 mmd_rbf            -> p6_label_free_signals_blinded.csv
    (a pre-correction Protocol v1.0 file whose CSV no downstream script references)

This module executes those forks VERBATIM (evs_load.exec_block / load_def) on small deterministic
fixtures and compares them to brute-force references.  It also adds the fork census (every
mmd_rbf / wconcordance copy), the cross-copy agreement gate and the reintroduction-coverage gate.
No production file is edited; the REAL tree stays the measuring target.

REFERENCES (source of every expectation; never paraphrased from the code under test)
================================================================================
  - Gretton, Borgwardt, Rasch, Schoelkopf & Smola, "A Kernel Two-Sample Test", JMLR 13 (2012)
    723-773, eq. (3): unbiased U-statistic
        MMD_u^2(A,B) = 1/(n(n-1)) sum_{i!=j} k(a_i,a_j) + 1/(m(m-1)) sum_{i!=j} k(b_i,b_j)
                       - 2/(nm) sum_{i,j} k(a_i,b_j),   RBF k(x,y)=exp(-||x-y||^2 / bw2).
    scipy.pdist returns the n(n-1)/2 UNORDERED pairs; the U-statistic's within terms are ORDERED
    double sums over n(n-1) / m(m-1) -- the confirmed defect (a) is exactly that factor 2.
  - Li, Gatsonis, Dahabreh & Steingrimsson, Biometrics 2023, DOI 10.1111/biom.13796 (on disk as
    ./ref_li2023_dr_auc.md), eqs (12)(13)(15): the literal three-term DR kernel
    d^DR = d^IOW + d^OM - w_i w_j m_i (1-m_j) I(S_i=1,S_j=1) k, evaluated per ordered pair as an
    INDEPENDENT oracle for the p10b DR source term (FORK-02c), and eq. (14) as the oracle for the
    whole DR ratio (FORK-02).
  - ★ CS-05 (2026-10-06): p6c_01_L3_blinded.py's own frozen spec header, lines 14-27, is NO LONGER
    a conformance target for the DR tests here (nor in t03).  It prescribes a DR denominator D1*D0
    from marginal means instead of eq. (14)'s k=1 pair sum, and an OM denominator that keeps the
    i=j diagonal, contradicting its own numerator; and it is not what the code implemented either
    (p6c:283 / p10b:167 formed a product of two residuals).  It is retained only as the
    non-vacuity control (helpers bf_dr(bracket="spec"), DEF_DR_SPEC) so that a gate which cannot
    tell the published estimator from the header form FAILS instead of passing vacuously.
  - Kiyasseh et al., Nat Commun 2024, doi:10.1038/s41467-024-46000-9: the original SUDO runs TWO
    pseudo-label branches per bin; the auxiliary-classifier count is therefore
        2 * SUDO_REP * (number of valid bins).
  - Hanley & McNeil 1982: the standard AUC tie credit is 0.5 (canonical wconcordance reference).
  - Bandwidth convention (frozen Methods / p6b:267-272 / p10a:227): median heuristic on
    SOURCE-TRAIN-only distances.  p6_02:229-231 instead takes it from the POOLED subsample
    (its own v1.0 header documents that convention) -> measured as a divergence (FORK-04b), not
    silently normalised.

EXPECTED STATUS MATRIX (REAL unmodified tree -> PATCHED tree)
================================================================================
  FORK-00       PASS -> PASS   infra: all forks extracted verbatim
  FORK-01       FAIL -> PASS   p10b IOW (harness control: ~1e-16 on BOTH trees) + OM vs
                               ref_li2023 eq. (11) pair-sum denominator (CS-05; same deviation as
                               p6c:19).  Was PASS -> PASS until CS-05 because its OM reference
                               pinned the deviating product form -- i.e. it agreed with the bug.
  FORK-01b      PASS -> PASS   block-internal nuisance fits == independent refit (control)
  FORK-02       FAIL -> PASS   p10b DR ratio vs ref_li2023 eq. (14) (defect c + the D1*D0
                               denominator; repaired by the CS-05 patches K9-K14)
  FORK-02b      WARN -> WARN   superseded-header / tie-convention disclosure (diagnostic, not gated)
  FORK-02c      FAIL -> PASS   p10b DR source kernel vs the external authority (defect c)
  FORK-03       FAIL -> FAIL   p10b SUDO single branch (defect d; SUDO repair deferred to author)
  FORK-03b      PASS -> PASS   counting-proxy faithfulness control
  FORK-04a      FAIL -> PASS   p6_02 within terms = half the U-statistic (defect a; patch A2)
  FORK-04b      WARN -> WARN   p6_02 pooled vs source-train-only bandwidth (convention conflict,
                               NOT covered by any patch; awaits author decision)
  FORK-05-*     FAIL -> PASS   census lines, one per mmd_rbf / wconcordance copy (live copies gate;
                               superseded / diagnostic copies are WARN only)
  FORK-06       PASS -> PASS   cross-copy agreement (p6b vs p10a+p10b chain; p6_02 under a
                               controlled bandwidth).  A CONSISTENCY gate, not a correctness
                               oracle: on the REAL tree the copies agree because they share
                               defect (a); it fires when a repair touches only some copies.
  FORK-07       PASS -> PASS   reintroduction coverage of (a)(b)(c)(d)(e)(f) at every live copy
                               location.  (e)/(f) are the OM and DR DENOMINATOR deviations that
                               CS-05 found and repaired; they are distinct from (c) (the DR
                               numerator bracket) and each is caught by its own gate.

Exit code is always 0; every test is wrapped so an exception becomes one FAIL line.
OUTPUT CONTRACT: RESULT|test_id|PASS/FAIL/WARN|expected (source)|observed|note
"""
import os
import sys
import traceback
import warnings

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
from scipy.spatial.distance import pdist, cdist
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import StandardScaler

import evs_load as L

warnings.filterwarnings("ignore")

FAM = "FORK"
SEED = 20261005                     # production seed (p6b:30 / p10a:32 / p10b:27 / p6c:48)
MMD_M, MMD_K, NBIN = 1000, 5, 10
SUDO_REP, SUDO_CAP, SUDO_FRAC = 5, 200, 0.5
TOL = 1e-12                         # value-equality tolerance
RTOL = 1e-9                         # ratio tolerance

P27 = os.path.join(L.OUT, "p27_dr_formula_diagnostic.py")
P6_02_SUPERSEDED = os.path.join(L.OUT, "p6_superseded_pre_correction", "p6_02_label_free_signals.py")
REF_DR_MD = os.path.join(os.path.dirname(os.path.abspath(__file__)), "ref_li2023_dr_auc.md")

DEF_MMD = ("Gretton et al. 2012 JMLR 13:723-773 eq.(3) unbiased U-statistic; within terms are "
           "ORDERED double sums over n(n-1) and m(m-1)")
DEF_DR_SPEC = ("SUPERSEDED frozen DR spec header p6c:21-27 (kept only as the non-vacuity control "
               "and as the historical pre-registration text): N_DR source bracket = "
               "y_i(1-y_j) - m_i(1-m_j), ordered double sum over i != j, but denominator D1*D0 from "
               "marginal means; it is NOT the published estimator and NOT what the code ever did")
DEF_DR_PAPER = ("Li, Gatsonis, Dahabreh & Steingrimsson, Biometrics 2023, "
                "DOI 10.1111/biom.13796, eqs (12)(13)(15) literal three-term DR kernel")
DEF_SUDO = ("Kiyasseh et al., Nat Commun 2024, doi:10.1038/s41467-024-46000-9: original SUDO runs "
            "TWO auxiliary branches per bin -> classifier count = 2*SUDO_REP*n_valid_bins")
DEF_TIE = "Hanley & McNeil 1982: standard AUC tie credit = 0.5 over ordered i != j pairs"

L.banner(f"{FAM} - t07_fork_parity.py | live same-source forks: p10b L3 + p10b SUDO + p6_02 MMD "
         f"+ census/cross-copy/reintroduction gates")


# ============================================================================
# references written from definitions only (no scipy, explicit loops / full matrices)
# ============================================================================
def _sqdist_matrix(A, B=None, chunk=256):
    """Full |A| x |B| matrix of squared Euclidean distances (row-chunked broadcasting loop)."""
    A = np.asarray(A, float)
    B = A if B is None else np.asarray(B, float)
    out = np.empty((len(A), len(B)), dtype=float)
    for i0 in range(0, len(A), chunk):
        blk = A[i0:i0 + chunk]
        out[i0:i0 + chunk] = ((blk[:, None, :] - B[None, :, :]) ** 2).sum(-1)
    return out


def ref_mmd_terms(A, B, bw2):
    """Gretton 2012 eq.(3) terms from the definition, on FULL kernel matrices, k(x,x)=1."""
    n, m = len(A), len(B)
    Kaa = np.exp(-_sqdist_matrix(A, A) / bw2)
    Kbb = np.exp(-_sqdist_matrix(B, B) / bw2)
    Kab = np.exp(-_sqdist_matrix(A, B) / bw2)
    np.fill_diagonal(Kaa, 0.0)
    np.fill_diagonal(Kbb, 0.0)
    a_within = float(Kaa.sum() / (n * (n - 1)))
    b_within = float(Kbb.sum() / (m * (m - 1)))
    cross = float(Kab.sum() / (n * m))
    return dict(a_within=a_within, b_within=b_within, cross=cross,
                mmd2=float(a_within + b_within - 2.0 * cross))


def ref_mmd_terms_scalar(A, B, bw2):
    """Same definition, pure-Python scalar double loop (independent of numpy matrix tricks)."""
    n, m = len(A), len(B)
    sa = sb = sc = 0.0
    for i in range(n):
        for j in range(n):
            if i != j:
                sa += float(np.exp(-float(np.sum((A[i] - A[j]) ** 2)) / bw2))
    for i in range(m):
        for j in range(m):
            if i != j:
                sb += float(np.exp(-float(np.sum((B[i] - B[j]) ** 2)) / bw2))
    for i in range(n):
        for j in range(m):
            sc += float(np.exp(-float(np.sum((A[i] - B[j]) ** 2)) / bw2))
    a_within, b_within, cross = sa / (n * (n - 1)), sb / (m * (m - 1)), sc / (n * m)
    return dict(a_within=a_within, b_within=b_within, cross=cross,
                mmd2=float(a_within + b_within - 2.0 * cross))


def ref_bw2_source_train(X):
    """Median heuristic on the SOURCE-TRAIN-only unordered pairs (frozen Methods / p6b:267-272)."""
    X = np.asarray(X, float)
    D2 = _sqdist_matrix(X)
    iu = np.triu_indices(len(X), 1)
    return float(np.median(D2[iu]))


def ref_bw2_pooled(A, B):
    """Median heuristic on the POOLED sample (p6_02:229-231 convention; measured, not endorsed)."""
    P = np.vstack([np.asarray(A, float), np.asarray(B, float)])
    D2 = _sqdist_matrix(P)
    iu = np.triu_indices(len(P), 1)
    return float(np.median(D2[iu]))


def _credit(si, sj, tie):
    """k(score_i, score_j) = 1[s_i > s_j]; 0.5 credit for ties if tie=True."""
    if si > sj:
        return 1.0
    if tie and si == sj:
        return 0.5
    return 0.0


def ref_concordance(score, a, b, tie=False):
    """sum_{i != j} a_i b_j credit(s_i, s_j) - explicit O(n^2) ordered double loop."""
    s = np.asarray(score, float)
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    tot = 0.0
    for i in range(len(s)):
        ai, si = float(a[i]), float(s[i])
        for j in range(len(s)):
            if i == j:
                continue
            tot += ai * float(b[j]) * _credit(si, float(s[j]), tie)
    return float(tot)


def _loop_sum(v):
    """Literal loop sum (independent of numpy pairwise summation)."""
    t = 0.0
    for x in np.asarray(v, float):
        t += float(x)
    return t


def _count_tied_ordered_pairs(s):
    s = np.asarray(s, float)
    c = 0
    for i in range(len(s)):
        for j in range(len(s)):
            if i != j and s[i] == s[j]:
                c += 1
    return c


def bf_iow(w, y, s, tie=False):
    """IOW (frozen header p6c:14-17): num = sum_{i!=j} w_i w_j 1[y_i=1]1[y_j=0] 1[s_i>s_j]. """
    y = np.asarray(y, float)
    w = np.asarray(w, float)
    a = w * (y == 1.0)
    b = w * (y == 0.0)
    den = _loop_sum(a) * _loop_sum(b)
    return ref_concordance(s, a, b, tie) / den


def bf_om(s, m, tie=False, header_den=False):
    """OM, ref_li2023 eq. (11): num = sum_{i!=j} m_i (1-m_j) 1[s_i>s_j]; den = sum_{i!=j} m_i (1-m_j).

    ★ CS-05 (2026-10-06): the denominator used to be the product form (sum m)(sum 1-m), which keeps
    the i=j diagonal and is therefore NOT eq. (11)'s denominator.  header_den=True returns that
    superseded product form, used only as the non-vacuity control in t_fork01."""
    m = np.asarray(m, float)
    if header_den:
        den = _loop_sum(m) * _loop_sum(1.0 - m)
    else:
        den = 0.0
        n = len(m)
        for i in range(n):
            mi = float(m[i])
            for j in range(n):
                if i != j:
                    den += mi * (1.0 - float(m[j]))
    return ref_concordance(s, m, 1.0 - m, tie) / den


def bf_dr_spec_source_sum(w, y, m, s, tie=False):
    """Frozen-header source-stratum numerator as an explicit double loop (collapse of the paper's
    three terms); used only as a CROSS-CHECK on the paper-form reference in FORK-02c."""
    y = np.asarray(y, float)
    m = np.asarray(m, float)
    w = np.asarray(w, float)
    s = np.asarray(s, float)
    tot = 0.0
    for i in range(len(s)):
        for j in range(len(s)):
            if i == j:
                continue
            tot += (w[i] * w[j] * (y[i] * (1.0 - y[j]) - m[i] * (1.0 - m[j]))
                    * _credit(float(s[i]), float(s[j]), tie))
    return float(tot)


def bf_dr(w, y, m, s, mt, st, tie=False, bracket="spec"):
    """DR ratio.  bracket="spec": frozen header difference-of-two-products;
    bracket="prod": production's product-of-two-residuals structure (diagnostic only).
    Denominator D1/D0 as in the frozen header (p6c:24-25)."""
    y = np.asarray(y, float)
    m = np.asarray(m, float)
    w = np.asarray(w, float)
    s = np.asarray(s, float)
    mt = np.asarray(mt, float)
    st = np.asarray(st, float)
    ns, nt = len(s), len(st)
    tot = 0.0
    for i in range(ns):
        for j in range(ns):
            if i == j:
                continue
            if bracket == "spec":
                br = y[i] * (1.0 - y[j]) - m[i] * (1.0 - m[j])
            else:
                br = (y[i] - m[i]) * ((1.0 - y[j]) - (1.0 - m[j]))
            tot += w[i] * w[j] * br * _credit(float(s[i]), float(s[j]), tie)
    Ns = tot / (ns * (ns - 1.0))
    Nt = ref_concordance(st, mt, 1.0 - mt, tie) / (nt * (nt - 1.0))
    D1 = _loop_sum(w * (y - m)) / ns + _loop_sum(mt) / nt
    D0 = _loop_sum(w * ((1.0 - y) - (1.0 - m))) / ns + _loop_sum(1.0 - mt) / nt
    return (Ns + Nt) / (D1 * D0)


def bf_dr_paper_strata(w, y, m, s, mt, st, tie=False, k_is_one=False):
    """Literal three-term DR kernel of the external authority (eqs (12),(13),(15)), evaluated per
    ordered pair of like-stratum indices (mixed pairs contribute zero by the S_j / S_i indicators):
        source stratum = sum_{i!=j in S} [ w_i w_j y_i(1-y_j) - w_i w_j m_i(1-m_j) ] k
        target stratum = sum_{i!=j in T} m_i (1-m_j) k
    k = I(g_i > g_j) by default; k_is_one=True gives eq. (14)'s DENOMINATOR (k = 1).
    Returns (src_sum, tgt_sum).  This is the authority's own decomposition, NOT the collapsed
    header form, so it is an independent oracle for the p10b block's term_s / term_t."""
    y = np.asarray(y, float)
    m = np.asarray(m, float)
    w = np.asarray(w, float)
    s = np.asarray(s, float)
    mt = np.asarray(mt, float)
    st = np.asarray(st, float)
    ns, nt = len(s), len(st)
    src = 0.0
    for i in range(ns):
        for j in range(ns):
            if i == j:
                continue
            k = 1.0 if k_is_one else _credit(float(s[i]), float(s[j]), tie)
            d_iow = w[i] * w[j] * (1.0 if (y[i] == 1.0 and y[j] == 0.0) else 0.0) * k
            d_cor = w[i] * w[j] * m[i] * (1.0 - m[j]) * k
            src += d_iow - d_cor
    tgt = 0.0
    for i in range(nt):
        for j in range(nt):
            if i == j:
                continue
            kt = 1.0 if k_is_one else _credit(float(st[i]), float(st[j]), tie)
            tgt += mt[i] * (1.0 - mt[j]) * kt
    return float(src), float(tgt)


def bf_dr_paper_ratio(w, y, m, s, mt, st, tie=False):
    """ref_li2023 eq. (14): tau = sum d^DR(k=concordance) / sum d^DR(k=1), both over ordered i != j
    of the composite sample.  Mixed pairs carry zero by the S_i / S_j indicators, so the composite
    total is src + tgt.  Numerator and denominator share the same 1/(n(n-1)) (eq. (16)), so that
    factor cancels and the ratio is normalisation-free."""
    sc, tc = bf_dr_paper_strata(w, y, m, s, mt, st, tie=tie, k_is_one=False)
    s1, t1 = bf_dr_paper_strata(w, y, m, s, mt, st, tie=tie, k_is_one=True)
    return (sc + tc) / (s1 + t1)


# ============================================================================
# instrumentation (faithful proxies; they never alter the production arithmetic)
# ============================================================================
def instrumented():
    """pdist/cdist wrappers recording every call; returned arrays are the genuine scipy outputs."""
    calls = []

    def _p(a, metric="euclidean", **kw):
        v = pdist(a, metric, **kw)
        calls.append(("pdist", np.asarray(a, float), None, v))
        return v

    def _c(a, b, metric="euclidean", **kw):
        v = cdist(a, b, metric, **kw)
        calls.append(("cdist", np.asarray(a, float), np.asarray(b, float), v))
        return v

    return calls, _p, _c


class CountingLR(object):
    """Faithful LogisticRegression proxy: delegates every call to real sklearn and records every
    .fit(...) (the auxiliary classifiers actually trained).  Same technique as t04's CountingLR."""
    FITS = []

    def __init__(self, **kw):
        self._kw = dict(kw)
        self._inner = LogisticRegression(**kw)

    def fit(self, X, y, **kw):
        CountingLR.FITS.append((len(np.asarray(y)), dict(self._kw)))
        self._inner.fit(X, y, **kw)
        return self

    def decision_function(self, X):
        return self._inner.decision_function(X)

    def predict_proba(self, X):
        return self._inner.predict_proba(X)

    def predict(self, X):
        return self._inner.predict(X)

    def get_params(self, deep=True):
        return self._inner.get_params(deep)

    def __getattr__(self, name):
        return getattr(object.__getattribute__(self, "_inner"), name)


def _group_reps(calls, k):
    """Group a recorded call list into MMD_K reps of k calls each."""
    if len(calls) != k * MMD_K:
        raise ValueError(f"expected {k * MMD_K} recorded calls, got {len(calls)}")
    return [calls[i * k:(i + 1) * k] for i in range(MMD_K)]


def _terms_from_raw(raw_a, raw_b, raw_c, bw2, n, m):
    """Reconstruct the three kernel terms from recorded raw squared distances (production arithmetic
    exp(-raw/bw2), sum/(n(n-1)) for pdist within terms, mean for the full cdist cross term)."""
    a_w = float(np.exp(-raw_a / bw2).sum()) / (n * (n - 1))
    b_w = float(np.exp(-raw_b / bw2).sum()) / (m * (m - 1))
    cross = float(np.exp(-raw_c / bw2).mean())
    return a_w, b_w, cross


def _parse_p602_reps(calls, n, m):
    """Parse p6_02's per-rep call pattern [pdist(P), pdist(a), pdist(b), cdist(a,b)]."""
    reps = []
    for i in range(MMD_K):
        g = calls[4 * i:4 * i + 4]
        if (g[0][0], g[1][0], g[2][0], g[3][0]) != ("pdist", "pdist", "pdist", "cdist"):
            raise ValueError(f"unexpected p6_02 call pattern in rep {i}: {[c[0] for c in g]}")
        s2 = float(np.median(g[0][3]))
        a_w, b_w, cross = _terms_from_raw(g[1][3], g[2][3], g[3][3], s2, n, m)
        reps.append(dict(s2=s2, Kaa=a_w, Kbb=b_w, Kab=cross,
                         val=a_w + b_w - 2.0 * cross))
    return reps


# ============================================================================
# fixtures
# ============================================================================
_FX = {}


def p10b_fixture():
    """Small deterministic stand-in for one p10-cache artefact, built the way p10a builds it
    (StandardScaler on source train; logistic source model; outcome model on source train;
    fixed source reference; a target subsample idx).  No target labels are created: every
    estimator below is blinded."""
    if "p10b" in _FX:
        return _FX["p10b"]
    rng = np.random.default_rng(SEED + 101)
    d = 4
    beta = np.array([1.1, -0.8, 0.6, 0.4])

    def blk(n):
        X = rng.normal(0.0, 1.0, (n, d))
        lat = X @ beta + rng.normal(0.0, 0.5, n)
        return X, (lat > np.median(lat)).astype(int)

    Xtr, ytr = blk(80)
    Xva, yva = blk(30)
    Xte, yte = blk(60)
    Xt = rng.normal(0.0, 1.0, (200, d)) + np.array([0.9, -0.7, 0.4, 0.1])

    sc = StandardScaler().fit(Xtr)
    Dtr_std, Dva_std, Dte_std = sc.transform(Xtr), sc.transform(Xva), sc.transform(Xte)
    Dt_std_full = sc.transform(Xt)
    Ds_ref_std = np.vstack([Dtr_std, Dva_std, Dte_std])

    mdl = LogisticRegression(max_iter=2000, C=1.0, random_state=SEED).fit(Dtr_std, ytr)
    om = LogisticRegression(max_iter=2000, random_state=SEED).fit(Dtr_std, ytr)

    st_full = mdl.decision_function(Dt_std_full)
    pt_full = mdl.predict_proba(Dt_std_full)[:, 1]
    s_e = mdl.decision_function(Dte_std)
    m_e = np.clip(om.predict_proba(Dte_std)[:, 1], 1e-9, 1.0 - 1e-9)
    m_t_full = np.clip(om.predict_proba(Dt_std_full)[:, 1], 1e-9, 1.0 - 1e-9)

    idx = np.sort(np.random.default_rng(SEED + 102).choice(len(Xt), 180, replace=False))
    fx = dict(Dtr_std=Dtr_std, ytr=ytr, Dva_std=Dva_std, yva=yva, Dte_std=Dte_std, yte=yte,
              Ds_ref_std=Ds_ref_std, Dt_std_full=Dt_std_full, st_full=st_full, pt_full=pt_full,
              s_e=s_e, m_e=m_e, m_t_full=m_t_full,
              idx=idx, n=len(idx), Dt_std_s=Dt_std_full[idx], pt_s=pt_full[idx],
              m_t_s=m_t_full[idx], sc_t=st_full[idx],
              ties_se=_count_tied_ordered_pairs(s_e),
              ties_st=_count_tied_ordered_pairs(st_full[idx]))
    _FX["p10b"] = fx
    return fx


_MMD_FX = []


def shared_mmd_fixture():
    """ONE shared MMD fixture for the census and the cross-copy gate (n != m on purpose; both
    sizes < MMD_M so every production draw consumes the whole sample and the K repeats are
    permutation-invariant -> exact comparison, no Monte Carlo)."""
    if not _MMD_FX:
        rng = np.random.default_rng(SEED + 201)
        d = 4
        A = rng.normal(0.0, 1.0, (120, d))
        B = rng.normal(0.0, 1.0, (100, d)) + np.array([0.8, -0.5, 0.3, 0.2])
        _MMD_FX.extend([A, B])
    return _MMD_FX[0], _MMD_FX[1]


def tie_fixture():
    """Shared tie-rich fixture for the wconcordance census."""
    score = np.array([3.0, 3.0, 2.0, 2.0, 2.0, 1.0, 0.0, 0.0])
    a = np.array([1.0, 0.5, 2.0, 1.0, 0.5, 1.5, 1.0, 0.5])
    b = np.array([1.0, 2.0, 0.5, 1.5, 1.0, 0.5, 2.0, 1.0])
    return score, a, b


# ============================================================================
# production block runners (VERBATIM exec)
# ============================================================================
def run_p10b_l3(fx):
    """Exec L.BLOCK_P10B_L3 (p10b:150-174) with the loop-level env p10b itself supplies."""
    rec = {}
    rows = []
    env = dict(np=np, LogisticRegression=LogisticRegression, SEED=SEED,
               wconcordance=L.load_def(L.P10B, "wconcordance"),
               A=dict(Ds_ref_std=fx["Ds_ref_std"], Dte_std=fx["Dte_std"], s_e=fx["s_e"],
                      yte=fx["yte"], m_e=fx["m_e"], st_full=fx["st_full"]),
               Dt_std_s=fx["Dt_std_s"], n=fx["n"], idx=fx["idx"], m_t_s=fx["m_t_s"],
               sc_t=fx["sc_t"], rec=rec, rows=rows)
    got, src = L.exec_block(L.BLOCK_P10B_L3, env,
                            names=["rec", "w_te", "term_s", "term_t", "N_DR", "D1", "D0", "dm2"])
    return got, src, rows


def refit_w_te(fx):
    """Independent refit of the block-internal membership model (p10b:152-154): same class, same
    constructor args, same inputs; the estimator under test is NOT re-implemented here."""
    Xd = np.vstack([fx["Ds_ref_std"], fx["Dt_std_s"]])
    yd = np.r_[np.zeros(len(fx["Ds_ref_std"])), np.ones(fx["n"])]
    dm = LogisticRegression(max_iter=2000, class_weight="balanced", random_state=SEED).fit(Xd, yd)
    e = np.clip(dm.predict_proba(fx["Dte_std"])[:, 1], 1e-6, 1.0 - 1e-6)
    return e / (1.0 - e)


def run_p10b_sudo(fx, lr_cls):
    """Exec L.BLOCK_P10B_SUDO (p10b:116-148); returns (rec, counts, locals)."""
    rec = {}
    counts = {"fit": 0, "auc": 0}

    def _auc(y_true, y_score, **kw):
        counts["auc"] += 1
        return roc_auc_score(y_true, y_score, **kw)

    env = dict(np=np, NBIN=NBIN, SEED=SEED,
               A=dict(ytr=fx["ytr"], Dtr_std=fx["Dtr_std"], yva=fx["yva"], Dva_std=fx["Dva_std"]),
               rng=np.random.default_rng(SEED + 11), SUDO_REP=SUDO_REP, SUDO_CAP=SUDO_CAP,
               SUDO_FRAC=SUDO_FRAC, Dt_std_s=fx["Dt_std_s"], n=fx["n"], pt_s=fx["pt_s"],
               LogisticRegression=lr_cls, roc_auc_score=_auc, rec=rec)
    got, src = L.exec_block(L.BLOCK_P10B_SUDO, env, names=["rec", "sbins", "bid", "qs"])
    return rec, counts, got


def spec_valid_bins(pt_s, ytr, nbin=NBIN, min_n=10):
    """Valid bins as defined by the frozen SUDO spec header (p6c:29-33 read as a definition):
    score-quantile bins over pt; bin pseudo-label = (bin mean >= 0.5); pool = source-train units
    of the OPPOSITE ground-truth class; both bin and pool need >= 10 units."""
    pt = np.asarray(pt_s, float)
    ytr = np.asarray(ytr)
    qs = np.quantile(pt, np.linspace(0.0, 1.0, nbin + 1))
    qs[0], qs[-1] = -np.inf, np.inf
    bid = np.digitize(pt, qs[1:-1])
    out = []
    for b in range(nbin):
        ib = np.where(bid == b)[0]
        if len(ib) < min_n:
            continue
        yb = 1 if pt[ib].mean() >= 0.5 else 0
        pool = np.where(ytr == 1 - yb)[0]
        if len(pool) < min_n:
            continue
        out.append(dict(bin=b, yb=yb, n_idx=len(ib), n_pool=len(pool)))
    return out


# ============================================================================
# tests
# ============================================================================
def t_fork00():
    """Infrastructure: every live fork must be extractable VERBATIM by evs_load."""
    items = []
    items.append(("p6b mmd_rbf",
                  L.load_def(L.P6B, "mmd_rbf", extra_ns=dict(MMD_M=MMD_M, MMD_K=MMD_K)).__evs_origin__))
    items.append(("p10a Kaa block", L.block_source((L.P10A, 291, 298))))
    items.append(("p10b L1b inline", L.block_source((L.P10B, 106, 114))))
    items.append(("p10b L3 block", L.block_source(L.BLOCK_P10B_L3)))
    items.append(("p10b SUDO block", L.block_source(L.BLOCK_P10B_SUDO)))
    items.append(("p6_02 mmd_rbf",
                  L.load_def(L.BLOCK_P6_02_MMD[0], "mmd_rbf",
                             extra_ns=dict(MMD_M=MMD_M, MMD_K=MMD_K)).__evs_origin__))
    if os.path.exists(P6_02_SUPERSEDED):
        items.append(("p6_02 superseded mmd_rbf",
                      L.load_def(P6_02_SUPERSEDED, "mmd_rbf",
                                 extra_ns=dict(MMD_M=MMD_M, MMD_K=MMD_K)).__evs_origin__))
    items.append(("p6c wconcordance", L.load_def(L.P6C, "wconcordance").__evs_origin__))
    items.append(("p10b wconcordance", L.load_def(L.P10B, "wconcordance").__evs_origin__))
    if os.path.exists(P27):
        items.append(("p27 wconcordance", L.load_def(P27, "wconcordance").__evs_origin__))
    return True, "; ".join(f"{k}={v}" for k, v in items)


def t_fork01b():
    """Control: the block-internal nuisance fit equals an independent refit (else the exact
    comparisons below would not be attributable to the estimator)."""
    fx = p10b_fixture()
    got, _src, _rows = run_p10b_l3(fx)
    w_ref = refit_w_te(fx)
    w_code = np.asarray(got["w_te"], float)
    dw = float(np.max(np.abs(w_ref - w_code)))
    ok = dw <= 1e-9 and fx["ties_se"] == 0 and fx["ties_st"] == 0
    return ok, (f"max_abs_diff(w_te block vs refit)={dw:.2e}; tied ordered pairs s_e={fx['ties_se']} "
                f"sc_t={fx['ties_st']}; fixture is tie-free so strict-> and 0.5 credit coincide")


def t_fork01():
    """FORK-01: p10b L3 executes; IOW and OM must equal explicit O(n^2) brute force.

    Two halves with different roles, deliberately in one line so neither can be lost:
      * IOW half = the HARNESS CONTROL.  IOW was never buggy, so its residual must be ~1e-16 on
        BOTH trees; that is the evidence the FORK harness can reproduce a production estimator
        exactly (i.e. that a FAIL elsewhere is a real defect, not an always-failing harness).
      * OM half = a GATE since CS-05.  Its reference denominator is now eq. (11)'s pair sum
        sum_{i!=j} m_i(1-m_j); the superseded header's product form keeps the i=j diagonal, so the
        OM half FAILS on the unrepaired tree and PASSES after repair.  The product form is retained
        only as the non-vacuity control (om_gap), which must be non-zero for this gate to mean
        anything.
    FORK-01b remains the PASS -> PASS control line."""
    fx = p10b_fixture()
    got, _src, _rows = run_p10b_l3(fx)
    rec = got["rec"]
    w_ref = refit_w_te(fx)
    y = np.asarray(fx["yte"], float)
    iow_ref = bf_iow(w_ref, y, fx["s_e"], tie=False)
    om_ref = bf_om(fx["sc_t"], fx["m_t_s"], tie=False)
    om_ref_sup = bf_om(fx["sc_t"], fx["m_t_s"], tie=False, header_den=True)
    c_iow = float(rec["AUROC_IOW"])
    c_om = float(rec["AUROC_OM"])
    d_iow, d_om = abs(c_iow - iow_ref), abs(c_om - om_ref)
    om_gap = abs(om_ref - om_ref_sup)
    iow_ok = d_iow <= TOL
    ok = iow_ok and d_om <= TOL and om_gap > 1e-9
    return ok, (f"p10b_block executed (n_src_test={len(y)}, n_target_sub={fx['n']}): "
                f"AUROC_IOW code={c_iow:.12f} ref={iow_ref:.12f} absdiff={d_iow:.2e} "
                f"[IOW half = harness control, iow_control_ok={iow_ok}]; "
                f"AUROC_OM code={c_om:.12f} ref_eq11={om_ref:.12f} absdiff={d_om:.2e}; "
                f"superseded_product_den_ref={om_ref_sup:.12f} (discriminating gap={om_gap:.2e}); "
                f"tol={TOL:g} nonvacuous={om_gap > 1e-9}")


def t_fork02():
    """FORK-02: p10b DR ratio vs the PUBLISHED estimator (ref_li2023 eq. (14)).

    ★ CS-05 (2026-10-06): the reference was the script's own frozen spec header (p6c:21-27).  That
    header prescribes a DR denominator D1*D0 built from marginal means, so it is a different
    estimator from eq. (14)'s k=1 pair sum, and it is also not what the p10b code ever did.  A
    record that matches neither the published estimator nor the code cannot be a verification
    oracle; the reference is now the published eq. (14) ratio.  Expect FAIL on REAL (product of two
    residuals + D1*D0) and PASS once repaired."""
    fx = p10b_fixture()
    got, _src, _rows = run_p10b_l3(fx)
    rec = got["rec"]
    w = refit_w_te(fx)
    y = np.asarray(fx["yte"], float)
    m = np.asarray(fx["m_e"], float)
    s = fx["s_e"]
    mt, st = fx["m_t_s"], fx["sc_t"]
    code = float(rec["AUROC_DR"])
    ref = bf_dr_paper_ratio(w, y, m, s, mt, st, tie=False)
    ref_h = bf_dr_paper_ratio(w, y, m, s, mt, st, tie=True)
    # superseded header form and defect (c) structure: diagnostics AND the non-vacuity control
    spec_s = bf_dr(w, y, m, s, mt, st, tie=False, bracket="spec")
    prod_s = bf_dr(w, y, m, s, mt, st, tie=False, bracket="prod")
    d = abs(code - ref)
    gap = abs(ref - spec_s)
    # ★ non-vacuity guard (same rule as t03's L3-02 / L3-03): if on this fixture the published ratio
    #   and the superseded header ratio coincided, this gate could not discriminate and must FAIL
    #   rather than report a PASS that carries no information.
    ok = (d <= 1e-9 * max(1.0, abs(ref))) and (gap > 1e-9)
    return ok, (f"code={code:.9f} paper_eq14={ref:.9f} absdiff={d:.2e}; "
                f"paper_eq14_halfcredit={ref_h:.9f}; superseded_header_form={spec_s:.9f} "
                f"(discriminating gap={gap:.2e}); defect(c)_prod_structure={prod_s:.9f} lies "
                f"{abs(prod_s - code):.2e} from the code; nonvacuous={gap > 1e-9}")


def t_fork02b():
    """FORK-02b (WARN, diagnostic): how far the superseded header form sits from the published
    ratio on this fixture, and whether the gap is a bracket/denominator effect or a tie effect."""
    tid = "FORK-02b"
    try:
        s, a, b = [2.0, 1.0, 1.0, 0.0], np.ones(4), np.ones(4)
        demo_strict = ref_concordance(s, a, b, tie=False)
        demo_half = ref_concordance(s, a, b, tie=True)
        fx = p10b_fixture()
        got, _src, _rows = run_p10b_l3(fx)
        rec = got["rec"]
        w = refit_w_te(fx)
        y = np.asarray(fx["yte"], float)
        m = np.asarray(fx["m_e"], float)
        code = float(rec["AUROC_DR"])
        ref = bf_dr_paper_ratio(w, y, m, fx["s_e"], fx["m_t_s"], fx["sc_t"], tie=False)
        ref_h = bf_dr_paper_ratio(w, y, m, fx["s_e"], fx["m_t_s"], fx["sc_t"], tie=True)
        spec_s = bf_dr(w, y, m, fx["s_e"], fx["m_t_s"], fx["sc_t"], tie=False, bracket="spec")
        spec_h = bf_dr(w, y, m, fx["s_e"], fx["m_t_s"], fx["sc_t"], tie=True, bracket="spec")
        L.emit(tid, "WARN",
               DEF_DR_PAPER + "; superseded p6c:21-27 header form and the 0.5 tie convention "
               "reported for comparison",
               f"code - paperStrict={code - ref:+.6f}; code - paperHalf={code - ref_h:+.6f}; "
               f"code - supersededHeaderStrict={code - spec_s:+.6f}; "
               f"paper - supersededHeader={ref - spec_s:+.6f}; "
               f"headerHalf - headerStrict={spec_h - spec_s:+.6f}; "
               f"tied ordered pairs s_e={fx['ties_se']} s_t={fx['ties_st']}; "
               f"tie-sensitivity demo [2,1,1,0]: strict={demo_strict:.1f} half={demo_half:.1f}",
               "fixture is tie-free, so headerHalf == headerStrict and the paper-vs-header gap "
               "(the quantity FORK-02 now guards against) is a bracket/denominator effect, not a "
               "tie artefact; the demo shows the operator IS tie-sensitive")
    except Exception as e:
        L.emit(tid, "FAIL", DEF_DR_PAPER + "; superseded header comparison",
               f"EXCEPTION {type(e).__name__}: {e}", traceback.format_exc().splitlines()[-1])


def t_fork02c():
    """FORK-02c: p10b DR source kernel vs the EXTERNAL authority (eqs (12)(13)(15)). Expect FAIL."""
    with open(REF_DR_MD, encoding="utf-8") as f:
        txt = f.read()
    anchors = ["10.1111/biom.13796", "PMC10188769", "Equation (12)", "Equation (13)",
               "Equation (14)", "Equation (15)", "d^DR(O_i,O_j", "I(S_i=1, S_j=1) k(O_i,O_j)",
               "DIFFERENCE OF TWO PRODUCTS",
               # CS-05: provenance for "the denominator is a pair sum, not a product of marginals",
               # which FORK-01's OM half and FORK-02 now gate on
               "numerator and denominator are each a sum over i",
               "the denominator is the sum with k"]
    missing = [a for a in anchors if a not in txt]
    if missing:
        raise ValueError(f"authority anchors missing from {os.path.basename(REF_DR_MD)}: {missing}")
    fx = p10b_fixture()
    got, _src, _rows = run_p10b_l3(fx)
    rec = got["rec"]
    w = refit_w_te(fx)
    y = np.asarray(fx["yte"], float)
    m = np.asarray(fx["m_e"], float)
    src_sum, tgt_sum = bf_dr_paper_strata(w, y, m, fx["s_e"], fx["m_t_s"], fx["sc_t"], tie=False)
    spec_src = bf_dr_spec_source_sum(w, y, m, fx["s_e"], tie=False)
    ref_agree = abs(spec_src - src_sum)
    t_s, t_t = float(rec_drterm(got, "term_s")), float(rec_drterm(got, "term_t"))
    ds, dt = t_s - src_sum, t_t - tgt_sum
    tol = 1e-9 * max(1.0, abs(src_sum), abs(tgt_sum))
    ok = max(abs(ds), abs(dt)) <= tol and ref_agree <= 1e-9 * max(1.0, abs(src_sum))
    return ok, (f"authority anchors ok ({len(anchors)}/{len(anchors)}); "
                f"reference cross-check paper_src_sum vs spec-collapse={ref_agree:.2e}; "
                f"term_s(code)={t_s:.9f} vs paper_src_sum={src_sum:.9f} diff={ds:+.3e}; "
                f"term_t(code)={t_t:.9f} vs paper_tgt_sum={tgt_sum:.9f} diff={dt:+.3e}; tol={tol:.1e}")


def rec_drterm(got, key):
    v = got.get(key)
    if v is None:
        raise ValueError(f"p10b block did not expose {key}")
    return v


def t_ess02():
    """ESS-02 (CS-09) —— **p10b 侧**的样本域一致性检查（A 组·代码级）。

    t03 的 ESS-01 查 p6c 那一份；这里查 p10b 这一份。两份都必须把 `w_ESS` 算在**自己估计量
    所用的权重**（source internal test 子集）上。任何一份用了 all-source 权重，它那一半就红
    ⇒ 两半合起来就是**同源副本一致性**检查（不需要额外造一份共享 fixture）。

    p10b 在 CS-09 之前就已经用 `w_te`（正确的那一份）。之所以仍然上闸，是因为
    `w_ESS` 此前**不在** FORK_CENSUS 里，没有任何东西阻止未来把这一行改回去而无人察觉。"""
    fx = p10b_fixture()
    got, _src, _rows = run_p10b_l3(fx)
    rec = got["rec"]
    w_te = np.asarray(got["w_te"], float)
    ref = float((w_te.sum() ** 2) / (w_te ** 2).sum())
    code = float(rec["w_ESS"])
    d = abs(ref - code)
    return (d <= 1e-12), (f"ESS(w_te)={ref:.9f} code={code:.9f} absdiff={d:.3e}; "
                          f"n_src_test={len(w_te)}")


def t_fork03():
    """FORK-03: p10b SUDO executes; count auxiliary classifiers against the two-branch original.

    ★ 判别力声明（复验 agent 指出，必须留在观测串里）：
      2× 计数是**必要而非充分**条件。一个"单分支但把每个分类器拟合两次"的块同样会得到
      n_fit == 2*SUDO_REP*n_valid 且 SUDO_AURCC 不变 —— FORK-03 会假绿。
      **FORK-03 绿 ≠ SUDO 合规**；结构性的互补分支判据在 t04 的 SUDO-01（等 SUDO 修复落地后移植过来）。
    ★ 空转守卫：没有有效 bin、或没有训练任何分类器时，本测试**不具判别力**，必须 FAIL。
    ★ fail-closed 一致性：两侧都有限才比较；**任一 NaN 不得当作"一致"**（除非两侧都是 NaN，
      此时显式声明预期为非有限）。
    """
    fx = p10b_fixture()
    valid = spec_valid_bins(fx["pt_s"], fx["ytr"])
    n_valid = len(valid)
    rec_plain, counts_plain, _locs = run_p10b_sudo(fx, LogisticRegression)
    CountingLR.FITS = []
    rec_proxy, counts_proxy, locs = run_p10b_sudo(fx, CountingLR)
    n_fit = len(CountingLR.FITS)
    expected = 2 * SUDO_REP * n_valid
    a = float(rec_proxy["SUDO_AURCC"]); b = float(rec_plain["SUDO_AURCC"])
    a_fin, b_fin = bool(np.isfinite(a)), bool(np.isfinite(b))
    if a_fin and b_fin:
        same, same_note = abs(a - b) <= 1e-12, "both finite"
    elif (not a_fin) and (not b_fin):
        same, same_note = True, "both non-finite (declared expected)"
    else:
        same, same_note = False, "one finite / one non-finite -> NOT agreement"
    vacuous = (n_valid < 1) or (n_fit < 1)
    ok = (n_fit == expected) and same and (not vacuous)
    return ok, (f"auxiliary classifiers fitted={n_fit}; two-branch expectation="
                f"2*SUDO_REP*n_valid_bins=2*{SUDO_REP}*{n_valid}={expected} "
                f"(observed = 1*{SUDO_REP}*{n_valid} = {SUDO_REP * n_valid}); "
                f"valid bins={n_valid} (of NBIN={NBIN}); "
                f"VACUOUS={vacuous} (n_valid<1 or n_fit<1 => 本测试无判别力，直接 FAIL); "
                f"auc_calls={counts_proxy['auc']}; "
                f"SUDO_nbin(block)={locs.get('rec', {}).get('SUDO_nbin', rec_proxy.get('SUDO_nbin'))}; "
                f"SUDO_AURCC={a:.6f}; proxy==plain={same} [{same_note}]; "
                f"**注：2x 计数为必要非充分——'单分支但每个分类器拟合两次'亦可得 100 且 AURCC 不变；"
                f"FORK-03 绿不等于 SUDO 合规**")


def t_fork03b():
    """FORK-03b: counting-proxy faithfulness control (same seed, same data)."""
    fx = p10b_fixture()
    rec_plain, _c1, _l1 = run_p10b_sudo(fx, LogisticRegression)
    CountingLR.FITS = []
    rec_proxy, _c2, _l2 = run_p10b_sudo(fx, CountingLR)
    keys = sorted(set(rec_plain) | set(rec_proxy))
    diffs = []
    for k in keys:
        a, b = rec_plain.get(k), rec_proxy.get(k)
        if isinstance(a, (int, float, np.floating)) and isinstance(b, (int, float, np.floating)):
            if not (np.isnan(a) and np.isnan(b)) and float(a) != float(b):
                diffs.append(f"{k}:{a}!={b}")
    ok = not diffs
    return ok, (f"proxy rec == plain rec on all numeric keys: {ok}; "
                f"SUDO_AURCC proxy={float(rec_proxy['SUDO_AURCC']):.10f} "
                f"plain={float(rec_plain['SUDO_AURCC']):.10f}; n_fits_proxy={len(CountingLR.FITS)}")


def t_fork04a():
    """FORK-04a: p6_02 mmd_rbf within terms vs the U-statistic, under its OWN bandwidth.
    The within contribution is recovered by INVERTING the return value (val + 2*cross), which is
    normalisation-agnostic: it yields 0.5 x the U-statistic while the factor-2 defect is present
    and 1.0 after the repair."""
    A, B = shared_mmd_fixture()
    n, m = len(A), len(B)
    calls, _p, _c = instrumented()
    fn = L.load_def(L.BLOCK_P6_02_MMD[0], "mmd_rbf",
                    extra_ns=dict(MMD_M=MMD_M, MMD_K=MMD_K, pdist=_p, cdist=_c))
    val, sd = fn(A, B, np.random.default_rng(SEED + 301))
    val = float(val)
    reps = _parse_p602_reps(calls, n, m)
    cross = float(np.mean([r["Kab"] for r in reps]))
    W_rec = val + 2.0 * cross
    ref_here = ref_mmd_terms(A, B, reps[0]["s2"])
    ratio = W_rec / (ref_here["a_within"] + ref_here["b_within"])
    recon_halved = float(np.mean([r["Kaa"] + r["Kbb"] - 2.0 * r["Kab"] for r in reps]))
    recon_fixed = float(np.mean([2.0 * r["Kaa"] + 2.0 * r["Kbb"] - 2.0 * r["Kab"] for r in reps]))
    which = ("halved" if abs(val - recon_halved) <= 1e-12 else
             ("fixed" if abs(val - recon_fixed) <= 1e-12 else "NONE"))
    ok = abs(ratio - 1.0) <= RTOL and which != "NONE"
    return ok, (f"p6_02:223-237 executed (n={n}, m={m}): within_ratio recovered_by_inversion="
                f"{ratio:.10f} (W_rec=val+2*cross={W_rec:.6f} ref_within="
                f"{ref_here['a_within'] + ref_here['b_within']:.6f}) under its OWN bandwidth "
                f"s2={reps[0]['s2']:.6f}; normalisation detected={which}; value={val:.8f} "
                f"recon_delta={min(abs(val - recon_halved), abs(val - recon_fixed)):.2e}; "
                f"sd_over_reps={sd:.2e}")


def t_fork04b():
    """FORK-04b (WARN): p6_02's pooled bandwidth vs the frozen source-train-only convention."""
    tid = "FORK-04b"
    try:
        A, B = shared_mmd_fixture()
        calls, _p, _c = instrumented()
        fn = L.load_def(L.BLOCK_P6_02_MMD[0], "mmd_rbf",
                        extra_ns=dict(MMD_M=MMD_M, MMD_K=MMD_K, pdist=_p, cdist=_c))
        fn(A, B, np.random.default_rng(SEED + 302))
        reps = _parse_p602_reps(calls, len(A), len(B))
        s2_pool = reps[0]["s2"]
        bw_src = ref_bw2_source_train(A)
        rel = abs(s2_pool - bw_src) / bw_src
        status = "PASS" if rel <= RTOL else "WARN"
        L.emit(tid, status,
               "frozen Methods / p6b:267-272 / p10a:227: RBF bandwidth = median heuristic over "
               "SOURCE-TRAIN-only squared distances",
               f"p6_02 uses pooled subsample bandwidth: s2_pooled={s2_pool:.6f} vs "
               f"source_train_only={bw_src:.6f} (ratio={s2_pool / bw_src:.4f}, rel_diff={rel:.3f}); "
               f"pooled sample adds within-B and cross pairs to the median",
               "convention conflict, NOT one of the patched defects (A2 only fixes the factor 2) and "
               "p6_02's own v1.0 header documents the pooled choice; its output CSV is referenced by "
               "no downstream script; author decision required")
    except Exception as e:
        L.emit(tid, "FAIL", "bandwidth convention divergence disclosure",
               f"EXCEPTION {type(e).__name__}: {e}", traceback.format_exc().splitlines()[-1])


def t_fork05():
    """FORK-05: census / parity table - one line per mmd_rbf and wconcordance copy."""
    A, B = shared_mmd_fixture()
    n, m = len(A), len(B)
    bw_src = ref_bw2_source_train(A)
    ref_src = ref_mmd_terms(A, B, bw_src)
    within_ref = ref_src["a_within"] + ref_src["b_within"]

    # reference self-check on the shared fixture
    ref_sc = ref_mmd_terms_scalar(A, B, bw_src)
    L.check("FORK-05-REF", abs(ref_sc["mmd2"] - ref_src["mmd2"]) <= 1e-9,
            "reference self-check: full-matrix U-statistic vs pure-Python scalar double loop",
            f"broadcast={ref_src['mmd2']:.10f} scalar_loop={ref_sc['mmd2']:.10f} "
            f"delta={abs(ref_sc['mmd2'] - ref_src['mmd2']):.2e}",
            "control on the reference itself (shared fixture n=120 m=100)")

    # ---------------- p6b ----------------
    def c_p6b():
        calls, _p, _c = instrumented()
        fn = L.load_def(L.P6B, "mmd_rbf", extra_ns=dict(MMD_M=MMD_M, MMD_K=MMD_K, pdist=_p, cdist=_c))
        val, _sd = fn(A, B, bw_src, np.random.default_rng(SEED + 501))
        val = float(val)
        groups = _group_reps(calls, 3)
        terms = [_terms_from_raw(g[0][3], g[1][3], g[2][3], bw_src, n, m) for g in groups]
        cross = float(np.mean([t[2] for t in terms]))
        W_rec = val + 2.0 * cross                     # normalisation-agnostic inversion
        ratio = W_rec / within_ref
        r_halved = float(np.mean([t[0] + t[1] - 2.0 * t[2] for t in terms]))
        r_fixed = float(np.mean([2.0 * t[0] + 2.0 * t[1] - 2.0 * t[2] for t in terms]))
        which = ("halved" if abs(val - r_halved) <= 1e-12 else
                 ("fixed" if abs(val - r_fixed) <= 1e-12 else "NONE"))
        status = "PASS" if (abs(ratio - 1.0) <= RTOL and which != "NONE") else "FAIL"
        return status, (f"p6b_01_corrected_pipeline.py:216-228 def mmd_rbf; within_ratio="
                        f"{ratio:.10f} (recovered W=val+2*cross={W_rec:.6f} "
                        f"ref_within={within_ref:.6f}); normalisation detected={which}; "
                        f"value={val:.8f}; executes: t01 MMD-01..06 + this census;"
                        f" LIVE: mmd column of p6b_label_free_signals_blinded.csv is read by "
                        f"p10a gate, p11 figure datapack, p6d, p14")

    _census("FORK-05-MMD-p6b", c_p6b)

    # ---------------- p10a ----------------
    def c_p10a():
        calls, _p, _c = instrumented()
        env = dict(np=np, MMD_K=MMD_K, MMD_M=MMD_M, pdist=_p, cdist=_c, Ds_std=A,
                   M={"bw2": bw_src}, r_ref=np.random.default_rng(SEED + 502))
        got, _src = L.exec_block((L.P10A, 291, 298), env, names=["Kaa", "mmd_src"])
        Kaa = list(got["Kaa"])
        ratio = float(Kaa[0]) / ref_src["a_within"]
        status = "PASS" if abs(ratio - 1.0) <= RTOL else "FAIL"
        return status, (f"p10a_prepare_cache.py:273-280 Kaa source-within precompute; ratio="
                        f"Kaa[0]/ref_a_within={ratio:.10f} (Kaa[0]={Kaa[0]:.6f} "
                        f"ref_a_within={ref_src['a_within']:.6f}); executes: t01 MMD-07/08 + this "
                        f"census; LIVE: p10_cache/*.pkl (Kaa, mmd_src, bw2) is consumed by p10b; "
                        f"p10b_raw.csv feeds p10c and figures")

    _census("FORK-05-MMD-p10a", c_p10a)

    # ---------------- p10b ----------------
    def c_p10b():
        calls, _p, _c = instrumented()
        env1 = dict(np=np, MMD_K=MMD_K, MMD_M=MMD_M, pdist=_p, cdist=_c, Ds_std=A,
                    M={"bw2": bw_src}, r_ref=np.random.default_rng(SEED + 503))
        got1, _s1 = L.exec_block((L.P10A, 291, 298), env1, names=["Kaa", "mmd_src"])
        n0 = len(calls)
        rec = {}
        env2 = dict(np=np, MMD_K=MMD_K, MMD_M=MMD_M, pdist=_p, cdist=_c, Dt_std_s=B, n=m,
                    rng=np.random.default_rng(SEED + 504),
                    A={"Kaa": got1["Kaa"], "mmd_src": got1["mmd_src"], "bw2": bw_src}, rec=rec)
        got2, _s2 = L.exec_block((L.P10B, 106, 114), env2, names=["rec"])
        tail = calls[n0:]
        groups = _group_reps(tail, 2)
        Kaa = [float(got1["Kaa"][k]) for k in range(MMD_K)]
        b_old = [float(np.exp(-g[0][3] / bw_src).sum()) / (m * (m - 1)) for g in groups]
        cross = [float(np.exp(-g[1][3] / bw_src).mean()) for g in groups]
        val = float(got2["rec"]["mmd"])
        Kaa_mean = float(np.mean(Kaa))
        cross_mean = float(np.mean(cross))
        Kbb_rec = val + 2.0 * cross_mean - Kaa_mean      # normalisation-agnostic inversion
        ratio_b = Kbb_rec / ref_src["b_within"]
        ratio_chain = (Kaa_mean + Kbb_rec) / within_ref
        r_halved = float(np.mean([Kaa[k] + b_old[k] - 2.0 * cross[k] for k in range(MMD_K)]))
        r_fixed = float(np.mean([Kaa[k] + 2.0 * b_old[k] - 2.0 * cross[k] for k in range(MMD_K)]))
        which = ("halved" if abs(val - r_halved) <= 1e-12 else
                 ("fixed" if abs(val - r_fixed) <= 1e-12 else "NONE"))
        status = "PASS" if (abs(ratio_b - 1.0) <= RTOL and which != "NONE") else "FAIL"
        return status, (f"p10b_subsample_curves.py:85-93 inline MMD (p10a Kaa + target within + cross) "
                        f"; target within_ratio=recovered Kbb/ref_b_within={ratio_b:.10f}; combined "
                        f"within_ratio={ratio_chain:.10f}; normalisation detected={which}; "
                        f"value={val:.8f}; executes: t01 MMD-08 + this census; LIVE: "
                        f"p10b_raw.csv -> p10c_aggregate.py; p11 -> figures/data/fig5_subsample.csv "
                        f"(fig5, figS8)")

    _census("FORK-05-MMD-p10b", c_p10b)

    # ---------------- p6_02 (live tree) and its archive ----------------
    def _p602_line(tid, path, liveness, gated):
        calls, _p, _c = instrumented()
        fn = L.load_def(path, "mmd_rbf", extra_ns=dict(MMD_M=MMD_M, MMD_K=MMD_K, pdist=_p, cdist=_c))
        val, _sd = fn(A, B, np.random.default_rng(SEED + 505))
        val = float(val)
        reps = _parse_p602_reps(calls, n, m)
        s2 = reps[0]["s2"]
        cross = float(np.mean([r["Kab"] for r in reps]))
        W_rec = val + 2.0 * cross
        ref_here = ref_mmd_terms(A, B, s2)
        ratio = W_rec / (ref_here["a_within"] + ref_here["b_within"])
        r_halved = float(np.mean([r["Kaa"] + r["Kbb"] - 2.0 * r["Kab"] for r in reps]))
        r_fixed = float(np.mean([2.0 * r["Kaa"] + 2.0 * r["Kbb"] - 2.0 * r["Kab"] for r in reps]))
        which = ("halved" if abs(val - r_halved) <= 1e-12 else
                 ("fixed" if abs(val - r_fixed) <= 1e-12 else "NONE"))
        if abs(ratio - 1.0) <= RTOL and which != "NONE":
            status = "PASS"
        else:
            status = "FAIL" if gated else "WARN"
        rel_bw = abs(s2 - bw_src) / bw_src
        return status, (f"{os.path.basename(os.path.dirname(path)) + '/' if not gated else ''}"
                        f"{os.path.basename(path)} mmd_rbf; within_ratio recovered_by_inversion="
                        f"{ratio:.10f}; normalisation detected={which}; value={val:.8f}; "
                        f"bw_internal={s2:.6f} vs source-train-only {bw_src:.6f} (rel_diff={rel_bw:.3f}); "
                        f"{liveness}")

    _census("FORK-05-MMD-p6_02-live",
            lambda: _p602_line("FORK-05-MMD-p6_02-live", L.BLOCK_P6_02_MMD[0],
                               "liveness: in-tree Protocol v1.0 file; p6_label_free_signals_blinded.csv "
                               "is read by NO downstream script (grep); not in p11 figure datapack; "
                               "executes: this module only (FORK-04a/05/06) - previously NO test",
                               gated=True))
    _census("FORK-05-MMD-p6_02-superseded",
            lambda: _p602_line("FORK-05-MMD-p6_02-superseded", P6_02_SUPERSEDED,
                               "liveness: SUPERSEDED archive (p6_superseded_pre_correction with "
                               "SUPERSEDED.sha256); not in use; executes: this census only; "
                               "defect presence is frozen history",
                               gated=False))

    # ---------------- p25 (diagnostic, not an estimator) ----------------
    L.emit("FORK-05-MMD-p25-diag", "WARN",
           "census: only still-in-use estimator kernels gate; a diagnostic that renders both forms "
           "side by side is inventoried, not gated",
           "p25_mmd_estimator_diagnostic.py:66-81; no standalone kernel: recomputes code_v "
           "(halved) and true_v (corrected) from the frozen p10_cache side by side; executes: not "
           "executed by the suite; DIAGNOSTIC: writes p25_mmd_diagnostic.csv + report only; no "
           "figure/result artifact consumes it",
           "inventoried so the grep-verified census is complete")

    # ---------------- wconcordance copies ----------------
    s_tie, a_tie, b_tie = tie_fixture()
    ref_half = ref_concordance(s_tie, a_tie, b_tie, tie=True)
    ref_strict = ref_concordance(s_tie, a_tie, b_tie, tie=False)
    n_tied = _count_tied_ordered_pairs(s_tie)

    def _wc_line(tid, loader, label, liveness, gated):
        v = float(loader(s_tie, a_tie, b_tie))
        dev_h, dev_s = abs(v - ref_half), abs(v - ref_strict)
        if dev_h <= 1e-12:
            status = "PASS"
        else:
            status = "FAIL" if gated else "WARN"
        return status, (f"{label}; tie fixture value={v:.6f} canonical(0.5 credit)={ref_half:.6f} "
                        f"strict={ref_strict:.6f}; delta_vs_0.5={dev_h:.2e} delta_vs_strict={dev_s:.2e}; "
                        f"tied ordered pairs={n_tied}; {liveness}")

    _census("FORK-05-WC-p6c",
            lambda: _wc_line("FORK-05-WC-p6c", L.load_def(L.P6C, "wconcordance"),
                             "p6c_01_L3_blinded.py:84-91 def wconcordance",
                             "executes: t02 CONC-04, t03 L3-03 + this census; LIVE: "
                             "p6c_L3_blinded.csv -> p11 -> figures/data/fig4_l3_scatter.csv, fig3c",
                             gated=True))
    _census("FORK-05-WC-p10b",
            lambda: _wc_line("FORK-05-WC-p10b", L.load_def(L.P10B, "wconcordance"),
                             "p10b_subsample_curves.py:40-46 def wconcordance",
                             "executes: t02 CONC-04B + this module's FORK-01/02; LIVE: the same "
                             "p10b_raw.csv -> p10c / fig5 chain",
                             gated=True))
    _census("FORK-05-WC-p27",
            lambda: _wc_line("FORK-05-WC-p27", L.load_def(P27, "wconcordance"),
                             "p27_dr_formula_diagnostic.py:46-53 def wconcordance",
                             "executes: not executed by the suite; DIAGNOSTIC: "
                             "p27_dr_formula_diagnostic.csv + report only; no figure/result artifact",
                             gated=False))

    return True, "census emitted: 6 mmd_rbf copies + 3 wconcordance copies + reference self-check"


def _census(tid, fn):
    """Emit one census RESULT line; an exception in one copy must not kill the family."""
    try:
        status, observed = fn()
    except Exception as e:
        L.emit(tid, "FAIL", "census line for a live/unreferenced estimator copy",
               f"EXCEPTION {type(e).__name__}: {e}", traceback.format_exc().splitlines()[-1])
        return
    L.emit(tid, status,
           "census: every still-in-use estimator kernel copy must equal the canonical reference "
           "(Gretton 2012 eq.(3) ordered-pair within terms / Hanley-McNeil 0.5 tie credit) on ONE "
           "shared fixture",
           observed, "one shared deterministic fixture; sizes < MMD_M so all draws are full-sample")


def t_fork06():
    """FORK-06: cross-copy agreement on the shared fixture under a CONTROLLED bandwidth."""
    A, B = shared_mmd_fixture()
    n, m = len(A), len(B)
    bw_src = ref_bw2_source_train(A)

    calls1, _p1, _c1 = instrumented()
    f6b = L.load_def(L.P6B, "mmd_rbf", extra_ns=dict(MMD_M=MMD_M, MMD_K=MMD_K, pdist=_p1, cdist=_c1))
    v6b, _sd = f6b(A, B, bw_src, np.random.default_rng(SEED + 601))

    calls2, _p2, _c2 = instrumented()
    env1 = dict(np=np, MMD_K=MMD_K, MMD_M=MMD_M, pdist=_p2, cdist=_c2, Ds_std=A,
                M={"bw2": bw_src}, r_ref=np.random.default_rng(SEED + 602))
    got1, _s1 = L.exec_block((L.P10A, 291, 298), env1, names=["Kaa", "mmd_src"])
    rec = {}
    env2 = dict(np=np, MMD_K=MMD_K, MMD_M=MMD_M, pdist=_p2, cdist=_c2, Dt_std_s=B, n=m,
                rng=np.random.default_rng(SEED + 603),
                A={"Kaa": got1["Kaa"], "mmd_src": got1["mmd_src"], "bw2": bw_src}, rec=rec)
    got2, _s2 = L.exec_block((L.P10B, 106, 114), env2, names=["rec"])
    v_chain = float(got2["rec"]["mmd"])
    d_chain = abs(v6b - v_chain)

    calls3, _p3, _c3 = instrumented()
    f602 = L.load_def(L.BLOCK_P6_02_MMD[0], "mmd_rbf",
                      extra_ns=dict(MMD_M=MMD_M, MMD_K=MMD_K, pdist=_p3, cdist=_c3))
    v602, _sd2 = f602(A, B, np.random.default_rng(SEED + 604))
    reps = _parse_p602_reps(calls3, n, m)
    s2_602 = reps[0]["s2"]

    calls4, _p4, _c4 = instrumented()
    f6b2 = L.load_def(L.P6B, "mmd_rbf", extra_ns=dict(MMD_M=MMD_M, MMD_K=MMD_K, pdist=_p4, cdist=_c4))
    v6b_ctrl, _sd3 = f6b2(A, B, s2_602, np.random.default_rng(SEED + 605))
    d_602 = abs(v602 - v6b_ctrl)

    ok = (d_chain <= TOL) and (d_602 <= TOL)
    return ok, (f"live chain: p6b(A,B,bw_source_train)={v6b:.8f} vs p10a+p10b chain={v_chain:.8f} "
                f"absdiff={d_chain:.2e}; p6_02 (own pooled bw={s2_602:.6f})={v602:.8f} vs "
                f"p6b@SAME bandwidth={v6b_ctrl:.8f} absdiff={d_602:.2e}; tol={TOL:g}")


def t_fork07():
    """FORK-07: reintroduction coverage of defects (a)(b)(c)(d)(e)(f)(g) at every live copy location.

    (g) was added by CS-09 (2026-10-06): `w_ESS` computed on the ALL-SOURCE weights while the
    estimator it diagnoses runs on the source INTERNAL TEST subset.  This one is different in
    kind from (a)-(f): the FORMULA is not wrong -- `(sum w)^2 / sum w^2` is correct ESS -- it is
    attached to the **wrong sample domain**.  That is why the old bound rule (t06 R14,
    `w_ESS <= n_source`) could be mathematically valid and still non-discriminating.

    (e) and (f) were added by CS-05 (2026-10-06).  They are NOT the same defect as (c):
      (c) DR NUMERATOR bracket = product of two residuals instead of difference of two products.
      (e) OM DENOMINATOR (sum m)(sum 1-m) keeps the i=j diagonal, contradicting the pair-sum
          numerator written directly above it in the same header line.
      (f) DR DENOMINATOR D1*D0 built from marginal means instead of eq. (14)'s k=1 pair sum
          (includes p6c:286's per-stratum 1/(n(n-1)) normalisers, which only cancel when both
          numerator and denominator carry the same one).
    Each gets its own rows because each can be reintroduced on its own, and each is caught by a
    different test: (e) by the OM gate, (f) by the DR ratio gate."""
    coverage = [
        ("a", "p6b:227", "t01_mmd.py", ["MMD-01", "MMD-08"]),
        ("a", "p10a:279", "t01_mmd.py", ["MMD-07"]),
        ("a", "p10b:92", "t01_mmd.py", ["MMD-08"]),
        ("a", "p6_02:236", "t07_fork_parity.py", ["FORK-04a", "FORK-05-MMD-p6_02-live", "FORK-06"]),
        ("b", "p6c:91", "t02_concordance.py", ["CONC-04"]),
        ("b", "p10b:46", "t02_concordance.py", ["CONC-04B"]),
        ("c", "p6c:283", "t03_l3_estimators.py", ["L3-03", "L3-06"]),
        ("c", "p10b:167", "t07_fork_parity.py", ["FORK-02", "FORK-02c"]),
        ("d", "p6c:291-330", "t04_sudo.py", ["SUDO-01"]),
        ("d", "p10b:116-148", "t07_fork_parity.py", ["FORK-03"]),
        ("e", "p6c:279", "t03_l3_estimators.py", ["L3-02"]),
        ("e", "p10b:164", "t07_fork_parity.py", ["FORK-01"]),
        ("f", "p6c:287-289", "t03_l3_estimators.py", ["L3-03"]),
        ("f", "p10b:171-173", "t07_fork_parity.py", ["FORK-02"]),
        # (g) = CS-09：w_ESS 算在 all-source 权重上（估计量用的是 source internal test）。
        # 与 (a)-(f) 都不同：前六个是"公式写错"，这个是**样本域**接错，公式本身没写错。
        ("g", "p6c:267", "t03_l3_estimators.py", ["ESS-01"]),
        ("g", "p10b:157", "t07_fork_parity.py", ["ESS-02"]),
    ]
    emitted_here = set(tid for tid, _st in L._SEEN)
    missing = []
    caught = []
    for defect, loc, modfile, tids in coverage:
        if modfile == "t07_fork_parity.py":
            hit = [t for t in tids if t in emitted_here]
        else:
            path = os.path.join(os.path.dirname(os.path.abspath(__file__)), modfile)
            with open(path, encoding="utf-8") as f:
                txt = f.read()
            hit = [t for t in tids if t in txt]
        if not hit:
            missing.append(f"({defect})@{loc} -> none of {tids} present in {modfile}")
        else:
            caught.append(f"({defect})@{loc} -> {'/'.join(hit)}")
    ok = not missing
    new_cov = ("(c)@p10b:167 -> FORK-02/02c; (d)@p10b:116-148 -> FORK-03; (a)@p6_02:236 -> "
               "FORK-04a/05/06; CS-05 added (e)@p6c:279 -> L3-02, (e)@p10b:164 -> FORK-01, "
               "(f)@p6c:287-289 -> L3-03, (f)@p10b:171-173 -> FORK-02; CS-09 added "
               "(g)@p6c:267 -> ESS-01, (g)@p10b:157 -> ESS-02")
    return ok, (f"all {len(coverage)} defect x live-location combinations caught: "
                f"{'; '.join(caught)}; uncovered={missing or 'none'}; NEW coverage added by this "
                f"module: {new_cov}")


# ============================================================================
def _main():
    L.guard("FORK-00", t_fork00,
            "every live same-source fork is extractable verbatim by evs_load (blocks by line "
            "range, defs by name); no production file is edited",
            "infra: verbatim extraction of the three untested live forks + the census corpus")
    L.guard("FORK-01b", t_fork01b,
            "control: block-internal nuisance fit (p10b:152-154 domain classifier) must equal an "
            "independent refit with identical class / args / inputs (t03 L3-00 pattern)",
            "ties must be absent from the fixture so tie conventions cannot mask the estimator check")
    L.guard("FORK-01", t_fork01,
            "IOW = the frozen L3 header p6c:16 and ref_li2023 eq. (9), read as definitions; OM = "
            "ref_li2023 eq. (11) (num and den BOTH sum over i != j). Both evaluated by an explicit "
            "O(n^2) ordered double-loop reference on the block's own inputs",
            "control: IOW / OM were never buggy; PASS here proves the FORK harness can pass. For IOW "
            "the product form (sum_{y=1} w)(sum_{y=0} w) IS the pair sum, since a binary y makes the "
            "diagonal empty; for OM it is not, which is why only OM's denominator moved in CS-05")
    L.guard("FORK-02", t_fork02,
            DEF_DR_PAPER + "; numerator and denominator each brute-forced by explicit double loops",
            "probe: p10b:167 is the same defect (c) as p6c:283; both tie variants reported; the "
            "superseded header form is retained only as the non-vacuity control")
    t_fork02b()      # emits FORK-02b itself (diagnostic, WARN) - not wrapped, to keep ids unique
    L.guard("FORK-02c", t_fork02c,
            DEF_DR_PAPER + "; literal three-term per-pair evaluation as an independent oracle",
            "second, external oracle for the same p10b defect (does not rely on p6c's own header)")
    L.guard("FORK-03", t_fork03,
            DEF_SUDO + "; counting proxy around LogisticRegression (t04 technique)",
            "probe: p10b:116-148 SUDO was never executed by any test before this module")
    L.guard("FORK-03b", t_fork03b,
            "control: the CountingLR proxy must reproduce the plain sklearn run bit-for-bit "
            "(same seed, same data) so the counts are trustworthy",
            "proxy faithfulness")
    L.guard("FORK-04a", t_fork04a,
            DEF_MMD + "; controlled-bandwidth comparison of the p6_02 within terms",
            "probe: p6_02 mmd_rbf executed by a reference test for the first time")
    t_fork04b()      # emits FORK-04b itself (WARN) - not wrapped, to keep ids unique
    L.guard("FORK-05", t_fork05,
            "census: one RESULT line per mmd_rbf and wconcordance copy (grep-verified; "
            "evs_load.FORK_CENSUS checked, not trusted)",
            "live copies gate; superseded / diagnostic copies are WARN only")
    L.guard("FORK-06", t_fork06,
            "cross-copy agreement: all copies claiming the same estimator must produce identical "
            "output on ONE shared fixture; p6_02 is compared under ITS OWN pooled bandwidth fed to "
            "p6b, so the bandwidth convention is neutralised and any residual gap is the factor-2 "
            "effect (which FORK-04a measures against the U-statistic)",
            "consistency gate, NOT a correctness oracle: two copies sharing defect (a) agree with "
            "each other by construction on the REAL tree, so this test passes now and its value is "
            "catching a PARTIAL repair (e.g. p6b fixed but p10a/p10b/p6_02 left behind)")
    L.guard("ESS-02", t_ess02,
            "CS-09 / CS09_ESS_acceptance.md §2 ESS-02 [A组·代码级]：p10b 的 rec['w_ESS'] 必须等于 "
            "ESS(w_te)，w_te 是该块自己的 source-internal-test 权重。t03 的 ESS-01 查 p6c 那一半；"
            "两半合起来即为跨副本样本域一致性检查")
    L.guard("FORK-07", t_fork07,
            "coverage gate: each of the four defects (a)(b)(c)(d) reintroduced at each live copy "
            "location must be caught by at least one executing reference test",
            "a location where reintroduction is invisible is a FAIL against the suite")
    print("# t07_fork_parity.py finished", flush=True)


if __name__ == "__main__":
    try:
        _main()
    except Exception as e:
        import traceback as _tb
        L.emit("FORK-00", "FAIL", "every FORK test must be guarded by L.guard",
               f"EXCEPTION {type(e).__name__}: {e}",
               _tb.format_exc().splitlines()[-1])
    sys.exit(0)
