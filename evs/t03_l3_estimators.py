# -*- coding: utf-8 -*-
r"""t03_l3_estimators.py - EVS family L3: transported-AUC estimators (IOW / OM / DR)

Production under test (READ-ONLY, executed verbatim by evs_load):
    L.BLOCK_L3B = p6c_01_L3_blinded.py lines 257-289  ("L3-B: transported AUC")
    executed with L.exec_block(L.BLOCK_L3B, env, names=["rec", ...]) on the synthetic fixture.

================================================================================
WHAT THE BRUTE-FORCE REFERENCE IMPLEMENTS, AND FROM WHERE
================================================================================

IOW  (source: harness task spec, identical to protocol v1.2 §11'' step 3 as transcribed
      in p6c:16):
        num = sum_{i != j in S_test} w_i w_j 1[y_i=1] 1[y_j=0] 1[s_i > s_j]
        den = (sum_{S_test} w_i 1[y_i=1]) * (sum_{S_test} w_j 1[y_j=0])
        Ahat_IOW = num / den

OM   (source: the EXTERNAL AUTHORITY — ./ref_li2023_dr_auc.md eq. (11)):
        num = sum_{i != j in T} m_i (1 - m_j) 1[s_i > s_j]
        den = sum_{i != j in T} m_i (1 - m_j)     <- PAIR SUM; eq. (11)'s own denominator
        Ahat_OM = num / den
      Equivalently den = (sum_T m)(sum_T (1-m)) - sum_T m_i(1-m_i).  The two forms differ by the
      diagonal sum sum_i m_i(1-m_i), which vanishes only when every m is 0 or 1; on calibrated
      probabilities it is strictly positive (this fixture: 4.8109 = 1.07% of the product form).
      L3-02 gates the pair form and additionally requires the fixture to separate the two forms,
      so a fixture on which they coincide FAILS instead of passing vacuously.

DR   (source: the EXTERNAL AUTHORITY — ./ref_li2023_dr_auc.md eqs (12)-(15), plus that file's own
      statement that "numerator and denominator are each a sum over i≠j"):

     External authority: ./ref_li2023_dr_auc.md, the on-disk transcription (equations 1-29) of
     Li, Gatsonis, Dahabreh & Steingrimsson, Biometrics 2023, DOI 10.1111/biom.13796
     (PMC10188769).  Its DR kernel — equations (12)+(13)+(15), summed over ordered i != j of
     the composite sample with k(O_i,O_j) = I(g_i > g_j) (eq. (1)/(14)) — is evaluated
     literally by _bf_dr_paper and gated by test L3-06; L3-05 reads the file from disk.
     Collapsed over strata (ref file section 6, marked there as inference and re-derived here
     by the independent three-term loop):
         source pairs:  w_i w_j [ y_i(1-y_j) - m_i(1-m_j) ] k    <- DIFFERENCE of two products
         target pairs:  m_i (1-m_j) k
         mixed pairs:   0
     L3-03 gates the FULL ratio sum d^DR(k=conc)/sum d^DR(k=1) against the authority, and L3-06
     gates the numerator kernel separately.  Those two together already imply the denominator:
     if N_code = N_paper and N_code/D_code = N_paper/D_paper with N_paper != 0, then
     D_code = D_paper.  A separate denominator gate would add no information, and would have to
     name internal variables (d_s/d_t), which is what made this module brittle before.

     ★ CS-05 (2026-10-06) — THE FROZEN SPEC HEADER IS NO LONGER A CONFORMANCE TARGET HERE.
     p6c lines 19-27 prescribes (a) an OM denominator (sum m)(sum (1-m)) that keeps the i=j
     diagonal, so the header's own numerator (written sum_{i != j}) and its own denominator
     contradict each other; and (b) a DR denominator D1*D0 built from marginal means instead of
     the k=1 pair sum.  It is also not what the code ever implemented: p6c:283 formed a product
     of two residuals, which is not the header's difference of two products either (defect (c)).
     A record that agrees with neither the published estimator nor the code it claims to specify
     cannot serve as a verification oracle.  It is retained below as the SUPERSEDED
     pre-registration text, and used ONLY as a diagnostic (helper _bf_dr_spec, and the L3-03
     non-vacuity gap and the L3-07 disclosure) to quantify how far it deviated.  Superseded
     header, p6c lines 21-27, transcribed VERBATIM:

[DR]   augmented-IPW（分子与分母各自双重稳健）
       N_DR = (1/(n_s(n_s−1))) Σ_{i≠j∈S_test} w_i w_j [ y_i(1−y_j) − m_i(1−m_j) ] 1[s_i>s_j]
              + (1/(n_t(n_t−1))) Σ_{i≠j∈T} m_i (1−m_j) 1[s_i>s_j]
       D1_DR = (1/n_s) Σ_{S_test} w_i (y_i − m_i) + (1/n_t) Σ_T m_i
       D0_DR = (1/n_s) Σ_{S_test} w_i ((1−y_i) − (1−m_i)) + (1/n_t) Σ_T (1−m_i)
       Ahat_DR  = N_DR / (D1_DR · D0_DR)
       一致性：outcome model 或 density ratio 任一正确即可（双重稳健）。

      The header is written with 1[s_i>s_j] (strict); it does not state a tie convention,
      so the reference is evaluated BOTH with strict credit and with 0.5 credit for ties
      (Hanley–McNeil 1982 standard AUC convention) and both numbers are reported.
      "Σ_{i≠j}" is read as the ordered double sum over all i≠j (each ordered pair once),
      which is the same reading production's wconcordance implements.

================================================================================
WHY THIS REFERENCE IS NOT A COPY OF THE PRODUCTION ESTIMATOR
================================================================================
The reference is three explicit O(n^2) Python double loops over ordered index pairs
(helper `_bf_pair`, `_bf_dr_spec`, `_bf_dr_prod`, `_bf_dr_paper` below).  It shares NO data structure,
NO sorting/bincount trick and NO code with production's `wconcordance` (which groups by
unique score, uses argsort/bincount/cumsum).  Agreement is therefore evidence about the
estimator arithmetic, not a tautology.  The production code itself is never edited and
never re-implemented: it is extracted as text and exec'd by evs_load.

NUISANCE MODELS ARE RE-FITTED IN THIS TEST, ON PURPOSE (and it is NOT a copy of the
estimator): the production block fits its own domain classifier `dm` (p6c:260-261) and
outcome model `om` (p6c:275) INTERNALLY, so a brute force that only received (w, m, s)
would be unable to claim it saw the same weights.  This module therefore re-fits both
models with the identical class, identical constructor arguments, identical
random_state=SEED and identical inputs (Dtr / Ds_std / Dt / Dt_std), then feeds its OWN
fitted w and m to the O(n^2) reference.  Test L3-00 verifies elementwise that the refit
nuisance quantities equal the ones production computed inside the block; only then are
the comparisons exact.  What is copied is the nuisance-fitting step (deterministic
nuisance estimation, sanctioned by the task), never the estimator under test.

KNOWN CONFIRMED DEFECT UNDER TEST (defect (c) of the harness background)
-----------------------------------------------------------------------
Production DR line p6c:283 forms a PRODUCT OF TWO RESIDUALS inside ONE bilinear form,
    term_s = wconcordance(s_e, w_te*(y - m_e), w_te*((1-y) - (1-m_e)))
           = sum_{i!=j} w_i w_j (y_i - m_i)((1-y_j) - (1-m_j)) 1[s_i>s_j],
whereas the published estimator requires the DIFFERENCE OF TWO PRODUCTS,
sum w_i w_j [ y_i(1-y_j) - m_i(1-m_j) ] 1[s_i>s_j].  L3-03 and L3-06 must therefore FAIL on the
unrepaired tree and PASS only once the kernel is repaired (defect (c) caught twice, by two
independent oracles: the eq. (14) ratio and the numerator kernel respectively).

ALGEBRAIC PROOF FOR L3-04 (a := y - m throughout)
    [y_i(1-y_j) - m_i(1-m_j)] - [(y_i-m_i)((1-y_j)-(1-m_j))]
  = [y_i - y_i*y_j - m_i + m_i*m_j] - [(y_i-m_i)(m_j-y_j)]
  = y_i - y_i*y_j - m_i + m_i*m_j - y_i*m_j + y_i*y_j + m_i*m_j - m_i*y_j
  = y_i - m_i - y_i*m_j - m_i*y_j + 2*m_i*m_j
  and  a_i(1-m_j) - m_i*a_j
  = (y_i-m_i)(1-m_j) - m_i(y_j-m_j)
  = y_i - y_i*m_j - m_i + m_i*m_j - m_i*y_j + m_i*m_j
  = y_i - m_i - y_i*m_j - m_i*y_j + 2*m_i*m_j.  QED
Consequence: production's second factor is (1-y_j)-(1-m_j) = -(y_j-m_j) = -a_j, so the
production bracket is -a_i*a_j, which is NOT the spec bracket in general; the gap is
exactly a_i(1-m_j) - m_i*a_j.

EXTERNAL AUTHORITY WIRED IN (see L3-05/L3-06): Li, Gatsonis, Dahabreh & Steingrimsson,
Biometrics 2023, DOI 10.1111/biom.13796, transcribed on disk as ./ref_li2023_dr_auc.md.
L3-05 checks the transcription is present, complete, and still carries the two sentences this
module's denominator gates rest on (that each estimator is a ratio of two sum-over-i≠j terms,
and that the denominator is the k=1 sum).  L3-02, L3-03 and L3-06 gate production against it.
This module now claims paper-equality of the OM ratio and of the whole DR ratio, not merely
kernel-level conformance.

Exit code is always 0; every test is wrapped so an exception becomes one FAIL line.
"""
import sys, os, traceback, warnings

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np

import evs_load as L
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score

warnings.filterwarnings("ignore")

SEED = 20261005


# ============================================================================
# brute-force reference: explicit O(n^2) double loops, independent of wconcordance
# ============================================================================
def _credit(si, sj, tie):
    """1[s_i > s_j] with strict convention; optionally +0.5 for s_i == s_j."""
    if si > sj:
        return 1.0
    if tie and si == sj:
        return 0.5
    return 0.0


def _bf_pair(s, a, b, tie=False):
    """sum_{i != j} a_i b_j credit(s_i, s_j), literal double loop (ordered pairs)."""
    n = len(s)
    tot = 0.0
    for i in range(n):
        ai = float(a[i])
        si = float(s[i])
        for j in range(n):
            if i == j:
                continue
            tot += ai * float(b[j]) * _credit(si, float(s[j]), tie)
    return tot


def _loop_sum(v):
    """Literal loop sum (independent of numpy's pairwise summation)."""
    t = 0.0
    for x in v:
        t += float(x)
    return t


def _bf_dr_spec(w, y, m, s, mt, st, tie=False):
    """Frozen spec header (p6c:21-27), difference-of-two-products bracket."""
    ns, nt = len(s), len(st)
    Ns = 0.0
    for i in range(ns):
        wi, yi, mi, si = float(w[i]), float(y[i]), float(m[i]), float(s[i])
        for j in range(ns):
            if i == j:
                continue
            yj, mj = float(y[j]), float(m[j])
            Ns += wi * float(w[j]) * (yi * (1.0 - yj) - mi * (1.0 - mj)) * _credit(si, float(s[j]), tie)
    Nt = 0.0
    for i in range(nt):
        mti, sti = float(mt[i]), float(st[i])
        for j in range(nt):
            if i == j:
                continue
            Nt += mti * (1.0 - float(mt[j])) * _credit(sti, float(st[j]), tie)
    N = Ns / (ns * (ns - 1.0)) + Nt / (nt * (nt - 1.0))
    D1 = _loop_sum(w * (y - m)) / ns + _loop_sum(mt) / nt
    D0 = _loop_sum(w * ((1.0 - y) - (1.0 - m))) / ns + _loop_sum(1.0 - mt) / nt
    return N / (D1 * D0)


def _bf_dr_prod(w, y, m, s, mt, st, tie=False):
    """Production's structure (p6c:283): product of two residuals in one bilinear form,
    with D1/D0 as in the header.  Used only as a diagnostic to separate the tie effect
    from the bracket-structure effect."""
    ns, nt = len(s), len(st)
    Ns = 0.0
    for i in range(ns):
        ai = float(w[i]) * (float(y[i]) - float(m[i]))
        si = float(s[i])
        for j in range(ns):
            if i == j:
                continue
            bj = float(w[j]) * ((1.0 - float(y[j])) - (1.0 - float(m[j])))
            Ns += ai * bj * _credit(si, float(s[j]), tie)
    Nt = 0.0
    for i in range(nt):
        mti, sti = float(mt[i]), float(st[i])
        for j in range(nt):
            if i == j:
                continue
            Nt += mti * (1.0 - float(mt[j])) * _credit(sti, float(st[j]), tie)
    N = Ns / (ns * (ns - 1.0)) + Nt / (nt * (nt - 1.0))
    D1 = _loop_sum(w * (y - m)) / ns + _loop_sum(mt) / nt
    D0 = _loop_sum(w * ((1.0 - y) - (1.0 - m))) / ns + _loop_sum(1.0 - mt) / nt
    return N / (D1 * D0)


def _bf_dr_paper(w_s, y_s, m_s, s_s, m_t, s_t, tie=False, k_is_one=False):
    """Literal per-pair evaluation of the DR kernel of the EXTERNAL AUTHORITY
    (ref_li2023_dr_auc.md; Li, Gatsonis, Dahabreh & Steingrimsson, Biometrics 2023,
    equations (12), (13), (15)) on the composite sample (source test: S=1; target: S=0):

        d^DR(O_i,O_j) = d^IOW + d^OM - w(X_i,X_j) m(X_i,X_j) I(S_i=1,S_j=1) k(O_i,O_j)   (15)
        d^IOW = w(X_i,X_j) I(Y_i=1,Y_j=0,S_i=1,S_j=1) k(O_i,O_j)                        (12)
        d^OM  = m(X_i,X_j) I(S_i=0,S_j=0) k(O_i,O_j)                                    (13)

    with w(X_i,X_j) = w_i w_j (eq. (10) factored), m(X_i,X_j) = q_i (1-q_j) (eq. (6)/(11))
    and k(O_i,O_j) = I(g_i > g_j) strict (eq. (1)/(14); the ref file section 5 records that
    the paper prints no tie convention), or k = 1 when k_is_one=True (eq. (14) denominator,
    (16)-(17)).  Returns (total, source-stratum, target-stratum) ordered double sums.
    Target-observation w is never needed: every term carrying w has I(S_i=1,S_j=1).
    This is the authority's OWN three-term form, not the collapsed p6c header form."""
    ns, nt = len(s_s), len(s_t)
    S = np.r_[np.ones(ns), np.zeros(nt)]
    G = np.r_[np.asarray(s_s, float), np.asarray(s_t, float)]
    Y = np.r_[np.asarray(y_s, float), np.zeros(nt)]
    Q = np.r_[np.asarray(m_s, float), np.asarray(m_t, float)]
    W = np.r_[np.asarray(w_s, float), np.ones(nt)]
    n = ns + nt
    tot = src = tgt = 0.0
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            kf = 1.0 if k_is_one else _credit(G[i], G[j], tie)
            both_s = (S[i] == 1.0 and S[j] == 1.0)
            both_t = (S[i] == 0.0 and S[j] == 0.0)
            d_iow = W[i] * W[j] * (1.0 if (both_s and Y[i] == 1.0 and Y[j] == 0.0) else 0.0) * kf
            d_om = Q[i] * (1.0 - Q[j]) * (1.0 if both_t else 0.0) * kf
            d_cor = W[i] * W[j] * Q[i] * (1.0 - Q[j]) * (1.0 if both_s else 0.0) * kf
            d = d_iow + d_om - d_cor
            tot += d
            if both_s:
                src += d
            if both_t:
                tgt += d
    return tot, src, tgt


def _count_tied_ordered_pairs(s):
    n = len(s)
    c = 0
    for i in range(n):
        for j in range(n):
            if i != j and s[i] == s[j]:
                c += 1
    return c


# ============================================================================
# fixture
# ============================================================================
def _fixture_data():
    """Small hand-checkable synthetic data: n_source_test = 40, n_target = 60, d = 4.
    Source covariates N(0, I); target covariates shifted (covariate shift) so the domain
    classifier is non-trivial.  Source labels are balanced by construction (median split
    of a latent score that depends on X).  Target Y is NEVER used (blinded estimators)."""
    rng = np.random.default_rng(20261005)
    d = 4

    def src(n):
        X = rng.normal(0.0, 1.0, size=(n, d))
        lat = (1.2 * X[:, 0] - 0.9 * X[:, 1] + 0.7 * X[:, 2] + 0.6 * X[:, 3]
               + rng.normal(0.0, 0.6, size=n))
        y = (lat > np.median(lat)).astype(int)
        return X, y

    Xtr, ytr = src(60)
    Xva, yva = src(20)
    Xte, yte = src(40)
    Xt = rng.normal(0.0, 1.0, size=(60, d)) + np.array([1.1, -0.9, 0.5, 0.2])
    return dict(Xtr=Xtr, ytr=ytr, Xva=Xva, yva=yva, Xte=Xte, yte=yte, Xt=Xt)


def build_fixture():
    dat = _fixture_data()
    Dtr, ytr = dat["Xtr"], dat["ytr"]
    Dva, yva = dat["Xva"], dat["yva"]
    Dte, yte = dat["Xte"], dat["yte"]
    Dt = dat["Xt"]

    # ---- nuisance fits replicated EXACTLY as p6c:224-225 (class, args, random_state, inputs)
    sc = StandardScaler().fit(Dtr)
    mdl = LogisticRegression(max_iter=2000, C=1.0, random_state=SEED).fit(sc.transform(Dtr), ytr)
    M = dict(Dtr=Dtr, Dva=Dva, Dte=Dte, ytr=ytr, yva=yva, yte=yte, sc=sc, mdl=mdl, af=0.0)

    Ds_std = sc.transform(np.vstack([Dtr, Dva, Dte]))
    Dt_std = sc.transform(Dt)
    i_te = np.arange(len(yte)) + len(ytr) + len(yva)
    # p6c:255 lives OUTSIDE the block and is read by it at line 279
    sc_t = mdl.decision_function(sc.transform(Dt))

    rec = {}
    env = dict(np=np, LogisticRegression=LogisticRegression, roc_auc_score=roc_auc_score,
               SEED=SEED, M=M, Dt=Dt, Ds_std=Ds_std, Dt_std=Dt_std, i_te=i_te,
               sc_t=sc_t, wconcordance=L.load_def(L.P6C, "wconcordance"), rec=rec)
    out, src = L.exec_block(L.BLOCK_L3B, env,
                            names=["rec", "w_all", "s_e", "y_e", "m_t", "m_e",
                                   "term_s", "term_t", "N_DR", "D1", "D0"])

    # ---- refit the two nuisance models here (see docstring: sanctioned, not the estimator)
    Xd = np.vstack([Ds_std, Dt_std])
    yd = np.r_[np.zeros(len(Ds_std)), np.ones(len(Dt_std))]
    dm2 = LogisticRegression(max_iter=2000, class_weight="balanced", random_state=SEED).fit(Xd, yd)
    e_s2 = np.clip(dm2.predict_proba(Ds_std)[:, 1], 1e-6, 1.0 - 1e-6)
    w2 = e_s2 / (1.0 - e_s2)
    w_te2 = w2[i_te]

    om2 = LogisticRegression(max_iter=2000, random_state=SEED).fit(sc.transform(Dtr), ytr)
    m_t2 = np.clip(om2.predict_proba(sc.transform(Dt))[:, 1], 1e-9, 1.0 - 1e-9)
    m_e2 = np.clip(om2.predict_proba(sc.transform(Dte))[:, 1], 1e-9, 1.0 - 1e-9)

    return dict(dat=dat, M=M, rec=out["rec"], src=src,
                s_e=mdl.decision_function(sc.transform(Dte)), sc_t=sc_t,
                w_all=w2, w_te=w_te2, y_e=yte.astype(float),
                m_e=m_e2, m_t=m_t2, i_te=i_te,
                prod_w=out["w_all"], prod_m_t=out["m_t"], prod_m_e=out["m_e"],
                prod_term_s=out["term_s"], prod_term_t=out["term_t"],
                # N_DR exists on BOTH trees (REAL: per-stratum normalised; PATCHED: raw pair sum).
                # D1/D0 exist ONLY on REAL -- CS-05 removed those names from the block, so anything
                # reading them must tolerate None or it dies on the patched tree.
                prod_N_DR=out["N_DR"], D1=out["D1"], D0=out["D0"])


# ============================================================================
# tests
# ============================================================================
def _test_00(fx):
    dw = float(np.max(np.abs(fx["w_all"] - fx["prod_w"])))
    dmt = float(np.max(np.abs(fx["m_t"] - fx["prod_m_t"])))
    dme = float(np.max(np.abs(fx["m_e"] - fx["prod_m_e"])))
    ok = max(dw, dmt, dme) < 1e-9
    return ok, f"max_abs_dw={dw:.2e} max_abs_dm_t={dmt:.2e} max_abs_dm_e={dme:.2e}"


def _test_01(fx):
    s = fx["s_e"]
    y = fx["y_e"]
    w = fx["w_te"]
    a = np.where(y == 1.0, w, 0.0)
    b = np.where(y == 0.0, w, 0.0)
    num = _bf_pair(s, a, b, tie=False)
    den = _loop_sum(a) * _loop_sum(b)
    ref = num / den
    code = float(fx["rec"]["AUROC_IOW"])
    d = abs(ref - code)
    return (d <= 1e-12), f"brute={ref:.12f} code={code:.12f} absdiff={d:.2e}"


def _bf_om_den_pairs(m):
    """sum_{i != j} m_i (1 - m_j) -- literal ordered double loop, == ref_li2023 eq. (11) denominator.
    Shares no code with production's wconcordance (which groups by unique score)."""
    tot = 0.0
    n = len(m)
    for i in range(n):
        mi = float(m[i])
        for j in range(n):
            if i != j:
                tot += mi * (1.0 - float(m[j]))
    return tot


def _test_02(fx):
    s = fx["sc_t"]
    m = fx["m_t"]
    num = _bf_pair(s, m, 1.0 - m, tie=False)
    den = _bf_om_den_pairs(m)                            # eq. (11) denominator, i != j
    den_sup = _loop_sum(m) * _loop_sum(1.0 - m)          # superseded header denominator, keeps i=j
    ref, ref_sup = num / den, num / den_sup
    code = float(fx["rec"]["AUROC_OM"])
    d = abs(ref - code)
    gap = abs(ref - ref_sup)
    # ★ non-vacuity guard: on a fixture where the two denominators coincide this test cannot tell
    #   them apart, so a PASS would carry no information.  Must FAIL rather than pass vacuously.
    #   (A binary m vector is exactly such a fixture: sum_i m_i(1-m_i) = 0 makes the forms equal.)
    ok = (d <= 1e-12) and (gap > 1e-9)
    diag = den_sup - den
    return ok, (f"brute_pairs={ref:.12f} code={code:.12f} absdiff={d:.2e}; "
                f"superseded_header_form={ref_sup:.12f} (discriminating gap={gap:.2e}); "
                f"den_pairs={den:.4f} den_product={den_sup:.4f} diagonal={diag:.4f} "
                f"({100.0 * diag / den_sup:.2f}% of the product form); nonvacuous={gap > 1e-9}")


def _bf_ess(v):
    """ESS = (sum w)^2 / sum w^2, literal loop sums (Cauchy-Schwarz: 1 <= ESS <= N for w >= 0).
    Shares no code with production, which uses numpy .sum()/.sum()."""
    tot = 0.0
    sq = 0.0
    for x in np.asarray(v, float):
        xf = float(x)
        tot += xf
        sq += xf * xf
    return (tot * tot) / sq


def _test_ess01(fx):
    """ESS-01 (CS-09) -- p6c half of the cross-copy domain check.  ★ Phase A 规范
    (CS09_ESS_acceptance.md §2 ESS-01, A 组·代码级).

    Requirement: p6c's rec["w_ESS"] must be the ESS of the SAME weight vector its own
    transported-AUC estimators use -- w_te = w_all[i_te], the source INTERNAL TEST subset --
    not of w_all (all source).  Domain identity: the diagnostic's sample must equal the
    estimator's sample."""
    w_all = np.asarray(fx["w_all"], float)
    w_te = w_all[np.asarray(fx["i_te"], int)]
    ref = _bf_ess(w_te)
    code = float(fx["rec"]["w_ESS"])
    d = abs(ref - code)
    return (d <= 1e-12), (f"ESS(w_te)={ref:.9f} code={code:.9f} absdiff={d:.3e}; "
                          f"n_src_test={len(w_te)} n_source={len(w_all)}")


def _test_ess03(fx):
    """ESS-03 (CS-09) -- non-vacuity guard for ESS-01 (A 组·代码级).

    On this fixture the all-source ESS and the source-test ESS MUST actually differ.  If they
    coincided, ESS-01 could not tell the two domains apart and a PASS would carry no
    information -- so this must FAIL rather than report non-vacuity it cannot deliver.
    Same rule as L3-02 / L3-03: a criterion that cannot discriminate is not a criterion."""
    w_all = np.asarray(fx["w_all"], float)
    w_te = w_all[np.asarray(fx["i_te"], int)]
    a, b = _bf_ess(w_all), _bf_ess(w_te)
    gap = abs(a - b)
    return (gap > 1e-9), (f"ESS(w_all)={a:.6f} ESS(w_te)={b:.6f} gap={gap:.6f} "
                          f"ratio={a / b:.4f}; nonvacuous={gap > 1e-9}")


def _bf_dr_paper_ratio(w, y, m, s, mt, st, tie=False):
    """ref_li2023 eq. (14): tau = sum d^DR(k = concordance) / sum d^DR(k = 1), both sums over
    ordered i != j of the composite sample.  Numerator and denominator carry the SAME
    normalisation, so the paper's 1/(n(n-1)) (eq. (16)) cancels and the ratio is
    normalisation-free -- which is why dropping p6c's per-stratum 1/(n_s(n_s-1)) and
    1/(n_t(n_t-1)) from the repaired numerator does not change the estimand."""
    num, _, _ = _bf_dr_paper(w, y, m, s, mt, st, tie=tie, k_is_one=False)
    den, _, _ = _bf_dr_paper(w, y, m, s, mt, st, tie=tie, k_is_one=True)
    return num / den


def _test_03(fx):
    code = float(fx["rec"]["AUROC_DR"])
    w, y, m, s = fx["w_te"], fx["y_e"], fx["m_e"], fx["s_e"]
    mt, st = fx["m_t"], fx["sc_t"]
    ref = _bf_dr_paper_ratio(w, y, m, s, mt, st, tie=False)
    ref_half = _bf_dr_paper_ratio(w, y, m, s, mt, st, tie=True)
    sp_s = _bf_dr_spec(w, y, m, s, mt, st, tie=False)   # superseded header: diagnostic only
    pr_s = _bf_dr_prod(w, y, m, s, mt, st, tie=False)   # defect (c) structure: diagnostic only
    d_paper = abs(code - ref)
    gap = abs(ref - sp_s)
    # ★ non-vacuity guard (same rule as L3-02): the fixture must separate the published ratio from
    #   the superseded header ratio, otherwise this gate cannot discriminate and must FAIL rather
    #   than pass vacuously.
    ok = (d_paper <= 1e-12) and (gap > 1e-9)
    L.check("L3-03", ok,
            "ref_li2023_dr_auc.md eq. (14): tau_DR = sum_{i!=j} d^DR(k=I[g_i>g_j]) / "
            "sum_{i!=j} d^DR(k=1), three-term kernel eqs (12),(13),(15); the superseded "
            "p6c:21-27 header is NOT the reference (see module docstring, CS-05)",
            f"code={code:.9f} paper_eq14={ref:.9f} absdiff={d_paper:.2e}; "
            f"superseded_header_form={sp_s:.9f} (discriminating gap={gap:.2e}); "
            f"defect(c)_prod_structure={pr_s:.9f} lies {abs(pr_s - code):.2e} from the code; "
            f"paper_eq14_in_unit_interval={0.0 <= ref <= 1.0} nonvacuous={gap > 1e-9}",
            "the denominator is IMPLIED rather than separately gated: L3-06 pins the numerator "
            "kernel and this line pins the ratio, so N_code=N_paper and N_code/D_code="
            "N_paper/D_paper give D_code=D_paper directly")
    ties_se = _count_tied_ordered_pairs(s)
    ties_st = _count_tied_ordered_pairs(st)
    demo_s = np.array([2.0, 1.0, 1.0, 0.0])
    demo_on = np.ones(4)
    d_strict_demo = _bf_pair(demo_s, demo_on, demo_on, tie=False)
    d_half_demo = _bf_pair(demo_s, demo_on, demo_on, tie=True)
    L.emit("L3-03b", "WARN",
           "tie convention alternative for the published eq. (14): 0.5 credit when s_i = s_j "
           "(Hanley-McNeil 1982 AUC convention); ref file section 5 records that the paper prints "
           "no tie convention",
           f"code-paperStrict={code - ref:+.6f} code-paperHalf={code - ref_half:+.6f}; "
           f"code-supersededHeader={code - sp_s:+.6f}; tied ordered pairs s_e={ties_se} s_t={ties_st}",
           "no ties in these scores, so the two tie variants coincide; the operator IS tie-sensitive "
           "(artificial vector [2,1,1,0]: strict={:.1f} half={:.1f}), so L3-03's residual is a kernel "
           "difference and not a tie artefact; tie path exercised by L3-03c".format(
               d_strict_demo, d_half_demo))
    return True, "emitted L3-03 / L3-03b"


def _test_04(fx):
    # exhaustive grid (exact): y_i,y_j in {0,1}; m_i,m_j in {0,.25,.5,.75,1}
    grid = [0.0, 0.25, 0.5, 0.75, 1.0]
    worst = 0.0
    for yi in (0.0, 1.0):
        for yj in (0.0, 1.0):
            for mi in grid:
                for mj in grid:
                    lhs = (yi * (1 - yj) - mi * (1 - mj)) - (yi - mi) * ((1 - yj) - (1 - mj))
                    rhs = (yi - mi) * (1 - mj) - mi * (yj - mj)
                    worst = max(worst, abs(lhs - rhs))
    rng = np.random.default_rng(7)
    N = 200000
    yi = rng.integers(0, 2, N).astype(float)
    yj = rng.integers(0, 2, N).astype(float)
    mi = rng.random(N)
    mj = rng.random(N)
    lhs = (yi * (1 - yj) - mi * (1 - mj)) - (yi - mi) * ((1 - yj) - (1 - mj))
    rhs = (yi - mi) * (1 - mj) - mi * (yj - mj)
    worst_r = float(np.max(np.abs(lhs - rhs)))
    ok = max(worst, worst_r) < 1e-12
    return ok, f"maxabs_grid={worst:.2e} (100 exact quadruples) maxabs_random={worst_r:.2e} (N={N})"


def _authority_path():
    return os.path.join(os.path.dirname(os.path.abspath(__file__)), "ref_li2023_dr_auc.md")


def _authority_anchors():
    """Substrings that must survive in the transcription for it to still be the authority
    this module's OM and DR oracles are derived from: equations (9),(11),(12)-(15), the strict-tie
    note, and -- added by CS-05 -- the two sentences the DENOMINATOR gates now rest on.  Those two
    are load-bearing: if the transcription were edited to drop them, L3-02's pair-sum denominator
    and L3-03's ratio denominator would lose their provenance and L3-05 must go red."""
    return [
        "10.1111/biom.13796",
        "PMC10188769",
        "Equation (9)",
        "Equation (11)",
        "Equation (14)",
        "Equation (15)",
        "d^DR(O_i,O_j",
        "I(S_i=1, S_j=1) k(O_i,O_j)",
        "DIFFERENCE OF TWO PRODUCTS",
        "STRICT inequality",
        "product-of-two-residuals",
        # CS-05: provenance for "the denominator is a pair sum, not a product of marginals"
        "numerator and denominator are each a sum over i",
        "the denominator is the sum with k",
    ]


def _test_05():
    """The external authority is read from disk, not quoted from memory."""
    path = _authority_path()
    with open(path, encoding="utf-8") as f:
        txt = f.read()
    anchors = _authority_anchors()
    missing = [a for a in anchors if a not in txt]
    return (not missing), (
        f"{os.path.basename(path)}: {len(txt)} chars, {len(anchors) - len(missing)}/{len(anchors)} "
        f"anchors present" + (f"; MISSING={missing}" if missing else ""))


def _test_06(fx):
    """Production DR numerator kernel vs the EXTERNAL authority (eqs (12),(13),(15))."""
    w, y, m, s = fx["w_te"], fx["y_e"], fx["m_e"], fx["s_e"]
    mt, st = fx["m_t"], fx["sc_t"]
    ns, nt = len(s), len(st)
    _, src_sum, tgt_sum = _bf_dr_paper(w, y, m, s, mt, st, tie=False)
    p_ts, p_tt = fx["prod_term_s"], fx["prod_term_t"]
    if p_ts is not None and p_tt is not None:
        ds = float(p_ts) - src_sum
        dt = float(p_tt) - tgt_sum
        d = max(abs(ds), abs(dt))
        tol = 1e-9 * max(1.0, abs(src_sum), abs(tgt_sum))
        return (d <= tol), (
            f"term_s(code)={float(p_ts):.9f} vs paper_src_sum={src_sum:.9f} diff={ds:+.3e}; "
            f"term_t(code)={float(p_tt):.9f} vs paper_tgt_sum={tgt_sum:.9f} diff={dt:+.3e}; "
            f"tol={tol:.1e}")
    # production names changed: fall back to the ratio's implied numerator
    N_paper = src_sum / (ns * (ns - 1.0)) + tgt_sum / (nt * (nt - 1.0))
    N_code = float(fx["rec"]["AUROC_DR"]) * fx["D1"] * fx["D0"]
    d = N_code - N_paper
    return (abs(d) <= 1e-9 * max(1.0, abs(N_paper))), (
        f"term_s/term_t not exposed by the block; fallback: code*D1*D0={N_code:.9f} vs paper "
        f"kernel normalised as the header={N_paper:.9f} diff={d:+.3e}")


def _test_07(fx):
    """Numbers for the L3-07 disclosure (no assertion)."""
    w, y, m, s = fx["w_te"], fx["y_e"], fx["m_e"], fx["s_e"]
    mt, st = fx["m_t"], fx["sc_t"]
    num, _, _ = _bf_dr_paper(w, y, m, s, mt, st, tie=False)
    den, _, _ = _bf_dr_paper(w, y, m, s, mt, st, tie=False, k_is_one=True)
    auc = float(fx["rec"]["AUROC_DR"])
    # ★ The denominator the code actually used, recovered WITHOUT naming an internal variable.
    #   The block's own identity is AUROC_DR = N_DR / denominator on BOTH trees (on REAL N_DR is
    #   the per-stratum-normalised numerator and the denominator is D1*D0; on PATCHED both are raw
    #   pair sums), so N_DR/tau recovers whichever denominator that tree used.  Reading fx["D1"] /
    #   fx["D0"] instead is wrong on PATCHED -- CS-05 deleted those names and this line emitted
    #   "TypeError: unsupported operand type(s) for *: 'float' and 'NoneType'" until this fix.
    N_code = fx["prod_N_DR"]
    D_code = (float(N_code) / auc) if (N_code is not None and auc != 0.0) else float("nan")
    om_prod = float(mt.sum() * (1.0 - mt).sum())
    om_diag = float((mt * (1.0 - mt)).sum())
    return (f"paper(eq14) tau={num / den:.6f}; code tau={auc:.6f}; "
            f"paper den(sum k=1)={den:.9f} vs code den(N_DR/tau)={D_code:.9f} "
            f"diff={D_code - den:+.3e}; "
            f"OM den product={om_prod:.4f} vs i!=j pair sum={om_prod - om_diag:.4f}; "
            f"code tau within unit interval={0.0 <= auc <= 1.0} "
            f"(reported, NOT gated: with mixed-sign d^DR the ratio is not provably bounded, so a "
            f"unit-interval assertion would be a false guarantee)")


def _test_03c(fx):
    """Characterisation: exercise the production operator's tie path on a tied vector.
    The main fixture has zero tied ordered pairs, so without this a gross tie regression in
    wconcordance would pass L3-01/02/03 untouched."""
    wc = L.load_def(L.P6C, "wconcordance")
    rng = np.random.default_rng(3)
    n = 64
    s = np.round(rng.normal(0.0, 1.0, n), 1)          # many exact ties
    a = rng.random(n)
    b = rng.random(n)
    op = float(wc(s, a, b))
    strict = _bf_pair(s, a, b, tie=False)             # 0.0 credit (frozen header / ref_li2023 sec. 5)
    half = _bf_pair(s, a, b, tie=True)                # 0.5 credit (Hanley-McNeil / family CONC target)
    ok = abs(op - strict) <= 1e-9 or abs(op - half) <= 1e-9
    return ok, (f"tied vector n={n} tied_ordered_pairs={_count_tied_ordered_pairs(s)}: "
                f"operator={op:.6f} brute_strict0={strict:.6f} brute_half0.5={half:.6f} "
                f"op-strict={op - strict:+.3e} op-half={op - half:+.3e}")


def _main():
    L.banner("L3 - transported AUC estimators (IOW / OM / DR), p6c lines 257-289")
    try:
        fx = build_fixture()
    except Exception as e:
        L.emit("L3-00", "FAIL", "fixture construction: synthetic data + verbatim execution of "
               "L.BLOCK_L3B (p6c:257-289)",
               f"EXCEPTION {type(e).__name__}: {e}", "all L3 tests depend on this fixture")
        L.emit("L3-00b", "WARN", "coverage contract: a fixture failure must not look like a clean run",
               "L3-01/02/03/03b/03c/04/05/06/07 NOT RUN (early return)",
               "family is still RED via L3-00; this line makes the dropped coverage explicit")
        return

    L.guard("L3-00", lambda: _test_00(fx),
            "nuisance-model refit identity: same class / args / random_state / inputs as "
            "p6c:260-261 (domain) and p6c:275 (outcome)",
            "brute force must see the same w and m as production's internal fits")
    L.guard("L3-01", lambda: _test_01(fx),
            "IOW definition (harness task spec, p6c:16): num = sum_{i!=j in S_test} w_i w_j "
            "1[y_i=1] 1[y_j=0] 1[s_i>s_j]; den = (sum_S w 1[y=1])(sum_S w 1[y=0])")
    L.guard("L3-02", lambda: _test_02(fx),
            "OM definition, ref_li2023_dr_auc.md eq. (11): num = sum_{i!=j in T} m_i (1-m_j) "
            "1[s_i>s_j]; den = sum_{i!=j in T} m_i (1-m_j) -- a PAIR SUM. The superseded p6c:19 "
            "header's den = (sum_T m)(sum_T (1-m)) keeps the i=j diagonal and is a diagnostic, "
            "not the reference; this test also FAILS if the fixture cannot tell the two apart")
    try:
        _test_03(fx)
    except Exception as e:
        L.emit("L3-03", "FAIL",
               "ref_li2023_dr_auc.md eq. (14) DR ratio sum d^DR(k=conc)/sum d^DR(k=1) over "
               "ordered i!=j of the composite sample",
               f"EXCEPTION {type(e).__name__}: {e}", traceback.format_exc().splitlines()[-1][:180])
    try:
        ok_3c, obs_3c = _test_03c(fx)
    except Exception as e:
        ok_3c, obs_3c = False, f"EXCEPTION {type(e).__name__}: {e}"
    L.check("L3-03c", ok_3c,
            "characterisation (not a correctness gate): credit on tied scores must equal EITHER the "
            "strict convention (0; ref_li2023 section 5 records the paper prints no tie convention) "
            "OR the 0.5 Hanley-McNeil convention that family CONC targets; this line only flags "
            "incoherent tie credit",
            obs_3c,
            "the main fixture has zero tied pairs, so the operator's tie branch is exercised only here",
            warn_only=True)
    L.guard("L3-04", lambda: _test_04(fx),
            "algebraic identity [y_i(1-y_j)-m_i(1-m_j)] - [(y_i-m_i)((1-y_j)-(1-m_j))] "
            "= a_i(1-m_j) - m_i a_j with a := y-m")
    L.guard("L3-05", _test_05,
            "external authority on disk: ref_li2023_dr_auc.md (Li, Gatsonis, Dahabreh & Steingrimsson, "
            "Biometrics 2023, DOI 10.1111/biom.13796, PMC10188769) must exist and contain equations "
            "(9), (11), (12)-(15) and the DR three-term kernel L3-06's oracle is derived from",
            "read from disk at test time; an INDEPENDENT transcription, not a paraphrase of p6c")
    L.guard("L3-06", lambda: _test_06(fx),
            "external authority ref_li2023_dr_auc.md eqs (12),(13),(15), composite-sample ordered "
            "i!=j: source-pair DR kernel = w_i w_j [y_i(1-y_j) - m_i(1-m_j)] 1[g_i>g_j], target-pair "
            "kernel = m_i(1-m_j) 1[g_i>g_j] (paper's own difference-of-two-products ordering)",
            "oracle _bf_dr_paper is the literal three-term loop; denominator conventions are out of "
            "scope for this test (disclosed in L3-07)")
    try:
        obs_07 = _test_07(fx)
    except Exception as e:
        obs_07 = f"EXCEPTION {type(e).__name__}: {e}"
    L.emit("L3-07", "WARN",
           "characterisation (disclosure, not a gate): the denominator the code ACTUALLY used, "
           "recovered as N_DR/tau, versus the published denominator (ref_li2023 eq (14) k=1 pair "
           "sum; eq (11) for OM), plus the unit-interval status of tau. The unit interval is "
           "REPORTED ONLY: with mixed-sign per-pair weights d^DR the ratio is not provably bounded, "
           "so gating on it would assert something that is not a theorem",
           obs_07,
           "denominator locality diagnostic -- L3-03 gates the ratio and L3-06 the numerator, which "
           "together imply the denominator; this line says WHICH part moved when L3-03 goes red")

    # ---- CS-09 ESS domain identity（Phase A 规范见 FREEZE_.../CS09_ESS_acceptance.md）----
    # ESS-01 = p6c 侧；ESS-02 = p10b 侧（在 t07，那里已有 p10b 的 fixture）。
    # 两半合起来构成同源副本一致性检查：若一侧用 w_all、另一侧用 w_te，两半不可能同时绿。
    L.guard("ESS-01", lambda: _test_ess01(fx),
            "CS-09 / CS09_ESS_acceptance.md §2 ESS-01 [A组·代码级]：p6c 的 rec['w_ESS'] 必须等于 "
            "ESS(w_te)，其中 w_te = w_all[i_te] 是 AUROC_IOW / AUROC_DR **实际使用**的同一个权重"
            "向量（样本域一致性）；不得是 all-source 的 ESS(w_all)")
    L.guard("ESS-03", lambda: _test_ess03(fx),
            "CS-09 §2 ESS-03 [A组·代码级]：ESS-01 的非空转守卫 —— 本 fixture 上 ESS(w_all) 与 "
            "ESS(w_te) 必须真的不同，否则 ESS-01 无法区分两个样本域，必须判红而不是白给一个 PASS")


if __name__ == "__main__":
    try:
        _main()
    except Exception:
        L.emit("L3-99", "FAIL", "module must not raise an uncaught exception",
               traceback.format_exc().splitlines()[-1][:200], "")
