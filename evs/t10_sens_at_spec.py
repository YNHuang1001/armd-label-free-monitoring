# -*- coding: utf-8 -*-
r"""t10_sens_at_spec.py - EVS family SENS : the LIVE `sens_at_spec` copies

WHY THIS FILE EXISTS
================================================================================
`sens_at_spec` is one of the estimator kernels that NO test used to execute while staying live:

  * p6d_01_unseal_evaluation.py:263-270  ->  p6d_signal_eval.csv  column `sens_at_spec80`
        ->  p11_figure_datapack.py:147 copies it into fig2_forest.csv (Figure 2 datapack)
        ->  the manuscript promises exactly this statistic: Methods_v1.md §14
            "Harmful-transfer discrimination was summarized by AUROC, with sensitivity at a
             prespecified fixed specificity of 0.80."
  * p10b_subsample_curves.py:256-263 (nested def) -> p10b_subsample_curves.csv `sens80_mean`
  * p10c_aggregate.py:64-73 (top-level def)       -> p10c_subsample_curves.csv `sens80_mean`

COPY MEMBERSHIP — governed by ACTIVE_PATH_MANIFEST.yaml (path_id SENS_AT_SPEC_80), NOT by runtime
discovery.  Every copy below is executed VERBATIM from the measured tree via evs_load.
================================================================================
  P10B  p10b_subsample_curves.py:256-263    NESTED def inside run_task()  -> L.exec_block(BLOCK_P10B_SENS)
  P10C  p10c_aggregate.py:64-73             top-level def                 -> L.load_def(L.P10C, ...)
  P6D   p6d_01_unseal_evaluation.py:263-270 top-level def                 -> L.load_def(L.P6D, ...)
The third copy is the one that writes the manuscript column, so it is executed and compared here
instead of being left as an untested fork.
★ v1.4: the heuristic discovery line **SENS-07 was DELETED** (runtime scan of the measured tree for
`def sens_at_spec`).  Signature discovery failed in both directions (v1.3 adversarial round: a
renamed fork was invisible, and a copy under a colliding basename in a subdirectory was exempted),
and a fourth copy must be caught by **registering it in the manifest** (t14 MC-01/MC-02), not by
growing the pattern.  An absent named copy is now declared by t14, not by this module.

REFERENCES (source of every expectation; never paraphrased from the code under test)
================================================================================
  (a) DEFINITION (CORRECTION_SPEC CS-07 v1.3, admits one correct value).
      "Among all ACHIEVABLE thresholds, choose the one maximising sensitivity subject to achieved
      specificity >= s0; deterministic tie-break if tied."  The reference implements exactly that
      as an exhaustive search:
          achieved specificity at threshold t = #{neg : s < t} / n_neg   (positive iff s >= t)
          sensitivity                        = #{pos : s >= t} / n_pos
          candidate thresholds = every observed score value (the ROC operating points) PLUS the
              boundary just above max(s) (the "predict nothing positive" point; the same operating
              point sklearn.metrics.roc_curve reports as its inf threshold)
          tie-break (predefined, deterministic): among the spec-compliant candidates attaining the
              maximum sensitivity take the LARGEST achieved specificity; if several still tie, the
              largest threshold (== the last such candidate in an ascending scan)
      The candidate set is EXHAUSTIVE over achievable operating points: any real threshold t has
      the same (achieved_spec, sens) pair as the smallest observed score >= t (or as the boundary
      point when t > max(s)), because no observed value lies strictly between them -- so no
      operating point can be missed and the maximum is exact.  Equivalently the selected
      threshold sits just above the k-th smallest negative score, k = ceil(s0 * n_neg) (the
      smallest spec-compliant real threshold).  Standard ROC operating-point definition (Fawcett,
      Pattern Recognition Letters 27 (2006) 861-874).  A score-value-only search WITHOUT the
      boundary point would be incomplete whenever no observed value reaches the requested
      specificity (e.g. all-tied scores, or one negative at the top): measured -- nan instead of
      the defined value 0.0 on 2/300 fuzz datasets; sklearn.metrics.roc_curve keeps that boundary
      point as its thresholds[0].  SENS-09 cross-checks this reference against
      sklearn.metrics.roc_curve's own ROC point selection (CS-07 acceptance criterion ③).
  (b) DOCUMENTED CONVENTION.
      Protocol v1.0 FROZEN (LabelFree_AMR_Deployment_Benchmark_Protocol_v1.0_FROZEN.md) §7.2,
      line 248: "sensitivity at fixed specificity（**specificity 固定为 0.80**，预先写死）".
      Methods_v1.md §14: "sensitivity at a prespecified fixed specificity of 0.80".
      => the delivered operating point must have ACHIEVED specificity >= 0.80.  Nothing in the
      protocol or Methods defines a quantile-interpolation rule, so the interpolation/tie layer is
      an UNDOCUMENTED code convention that this module identifies by measurement.
  (c) The sign of the threshold search is fixed by the two call-site orientations, not assumed:
      p6d builds `r_failure` so that HIGHER = more likely to fail (its own header: "越大 = 位移越大
      = 越可能失败"; SUDO/transported-AUC are negated for that reason), and the copies predict
      positive as `s >= thr` with tp = 1[pred & y == 1].  Hence a higher threshold = fewer positives.

WHAT THIS MODULE MEASURES (never silently adopts the production convention as the definition)
================================================================================
  SENS-01  both/three copies execute on a shared fixture                      EXPECT PASS
  SENS-02  cross-copy agreement, exact (fixtures + 300-dataset fuzz)          EXPECT PASS
  SENS-03  definitional correctness by explicit counting + convention probe   EXPECT FAIL (S-1)
  SENS-04  boundary: ties, inclusivity, attainability at delivered spec 0.80  EXPECT FAIL (S-1)
  SENS-05  label/score alignment under NaN filtering                          EXPECT PASS
  SENS-06  degenerate inputs: per-copy contract, single-class / no-negative   EXPECT FAIL (S-2)
  SENS-08  impact on the frozen `sens_at_spec80` column (202 tasks)           EXPECT FAIL (S-1)
  SENS-08b characterisation: frozen column reproduced from frozen inputs      PASS -> WARN (pin)
  SENS-09  green-ability self-test + reference cross-check vs sklearn         EXPECT PASS

  Green-ability (rule: a test that cannot go GREEN on a correct implementation is a pin, not a
  test).  Every FAIL above is produced by values the COPIES return, so a corrected copy set turns
  them green by construction: SENS-03 compares against the definitional operating point taken over
  the ACHIEVABLE thresholds (CS-07 v1.3 target definition: observed score values + the boundary
  above max(s)), SENS-04 checks that the delivered value is attainable at specificity >= 0.80
  (never exceeding the definitional maximum; a non-finite value fails closed), SENS-06 checks
  cross-copy agreement (plus no fabricated value where the statistic is undefined, and a defined
  non-nan value where the statistic IS defined), SENS-08 applies the same attainability check to
  the frozen corpus.  SENS-09 does not take that on faith: it runs a definition-only synthetic
  implementation through the same predicates (fixtures + fuzz corpus + frozen corpus) and requires
  them to go green, while the one-token `1-spec -> spec` fix must still be rejected -- otherwise
  the criterion would have no power; it also cross-checks the reference itself against
  sklearn.metrics.roc_curve's ROC operating-point selection.  SENS-08b is the one deliberate
  current-state pin (frozen-column reproduction) and is emitted warn_only=True with the word
  "characterisation" in its expected field.  Acceptance self-test: point EVS_PROD_DIR at a
  corrected copy set -> all lines except SENS-08b (PASS -> WARN) go green.  A corrected copy must
  implement the CS-07 target definition -- equivalently, place the threshold just above the
  k-th smallest negative score (k = ceil(spec * n_neg)); nudging the interpolated spec-quantile by
  one ulp is NOT sufficient on tie-heavy data (measured in SENS-04 / SENS-09: delivered
  specificity < 0.80 on 132/300 fuzz datasets, value above the attainable maximum on 62/300).

  S-1 = threshold taken at the (1 - spec) quantile of the negative scores, i.e. the complementary
        operating point (achieved specificity ~0.20 instead of 0.80).  NOT one of the six confirmed
        defects and NOT covered by any existing patch -> this module is its acceptance criterion.
  S-2 = latent cross-copy divergence on degenerate inputs (p6d fabricates 0.0 / raises IndexError
        where p10b and p10c return nan; p10b raises where p10c returns nan).  Unreachable with the
        current binary `harmful` labels, so it does NOT change the frozen column today.
"""
import os
import re
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import evs_load as L

SPEC = 0.80
QMETHODS = ("linear", "lower", "higher", "nearest", "midpoint")

REF_P10B = "p10b_subsample_curves.py:256-263"
REF_P10C = "p10c_aggregate.py:64-73"
REF_P6D = "p6d_01_unseal_evaluation.py:263-270"
LIVE = (REF_P10B, REF_P10C, REF_P6D)
SHORT = {REF_P10B: "p10b", REF_P10C: "p10c", REF_P6D: "p6d"}

COPY_NOTE = ("registered live copies (ACTIVE_PATH_MANIFEST.yaml :: SENS_AT_SPEC_80): "
             "p10b:256-263 (nested def, exec_block), p10c:64-73 (top-level def), p6d:263-270 "
             "(top-level def, writes p6d_signal_eval.csv sens_at_spec80 -> fig2_forest.csv)")


# ------------------------------------------------------------------ verbatim loading
def _load_copies():
    """Load every live copy VERBATIM from the measured tree.  A copy whose file is absent from the
    tree being measured (e.g. a trimmed patched tree) is reported in `absent`, never invented;
    the manifest (t14 MC-01/MC-02) is what declares an absent copy, not this module.  The verbatim
    source of each copy is kept so that the caller can check the loaded text really carries the
    `def sens_at_spec` signature (anchor-drift guard: a shifted block must not silently execute a
    neighbouring computation)."""
    copies, absent, srcs = {}, {}, {}
    try:
        ns, src = L.exec_block(L.BLOCK_P10B_SENS, {"np": np}, names=["sens_at_spec"])
        if ns.get("sens_at_spec") is None:
            absent[REF_P10B] = "block executed but no sens_at_spec in it"
        else:
            copies[REF_P10B] = ns["sens_at_spec"]
            srcs[REF_P10B] = src
    except Exception as e:
        absent[REF_P10B] = "%s: %s" % (type(e).__name__, e)
    for ref, attr in ((REF_P10C, "P10C"), (REF_P6D, "P6D")):
        path = getattr(L, attr, None)
        if path is None:
            absent[ref] = "evs_load has no anchor %s" % attr
            continue
        try:
            fn = L.load_def(path, "sens_at_spec")
            copies[ref] = fn
            srcs[ref] = getattr(fn, "__evs_src__", "") or ""
        except Exception as e:
            absent[ref] = "%s: %s" % (type(e).__name__, e)
    return copies, absent, srcs


COPIES, ABSENT, SRCS = _load_copies()
PRESENT = [r for r in LIVE if r in COPIES]


# ------------------------------------------------------------------ fixtures
# MAIN: 10 negatives / 5 positives, no ties, no NaN.  neg[8] is the definitional threshold.
Y_MAIN = np.array([0] * 10 + [1] * 5, dtype=float)
S_MAIN = np.array([0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50,
                   0.30, 0.45, 0.60, 0.75, 0.90])

# BOUND: 6 negatives; the definitional threshold sits just above neg[4] = 0.70, giving achieved
# specificity 5/6 = 0.8333 and sensitivity 2/4 = 0.50.  This fixture probes the one-token fix
# (1-spec -> spec with the threshold kept at the interpolated quantile = 0.70): with the copies'
# comparison s >= thr the negative AT 0.70 is predicted positive, so the DELIVERED achieved
# specificity is 4/6 = 0.6667 < 0.80 -- a specificity shortfall, not a sensitivity overshoot --
# while its returned value 0.50 happens to COINCIDE with the definitional value, so the value
# alone cannot expose it here (the shortfall is counted directly in SENS-04).  The fixture also
# does not discriminate the >= vs > rule (SENS-04 identifies it on the fuzz corpus instead).
Y_BOUND = np.array([0] * 6 + [1] * 4, dtype=float)
S_BOUND = np.array([0.10, 0.40, 0.50, 0.60, 0.70, 0.80, 0.40, 0.45, 0.75, 0.90])

# TIE: heavy ties inside the negative class.
Y_TIE = np.array([0, 0, 0, 0, 0, 1, 1, 1], dtype=float)
S_TIE = np.array([0.10, 0.20, 0.20, 0.20, 0.30, 0.20, 0.30, 0.40])

# NAN: two non-finite scores; the two dropped rows carry DIFFERENT labels, so a mask applied to
# one array only is detectable.
Y_NAN = np.array([0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1], dtype=float)
S_NAN = np.array([0.90, 0.40, 0.60, 0.50, np.nan, 0.20, 0.80, np.nan, 0.30, 0.70, 0.50, 0.10])

FIXTURES = (("MAIN", Y_MAIN, S_MAIN), ("BOUND", Y_BOUND, S_BOUND),
            ("TIE", Y_TIE, S_TIE), ("NAN", Y_NAN, S_NAN))

_FUZZ_CACHE = {}


def fuzz(n=300, seed=20261005):
    """Deterministic negative-control corpus: binary y, both classes, heavy score ties."""
    key = (n, seed)
    if key not in _FUZZ_CACHE:
        rng = np.random.default_rng(seed)
        out = []
        while len(out) < n:
            k = int(rng.integers(7, 120))
            y = rng.integers(0, 2, k)
            if len(np.unique(y)) < 2:
                continue
            s = np.round(rng.normal(size=k), int(rng.integers(1, 3)))
            out.append((y.astype(float), s))
        _FUZZ_CACHE[key] = out
    return _FUZZ_CACHE[key]


# ------------------------------------------------------------------ references
def mask_pairs(y, s):
    """The mask that the copies apply: isfinite(s), applied to BOTH arrays."""
    y = np.asarray(y)
    s = np.asarray(s, dtype=float)
    m = np.isfinite(s)
    return y[m], s[m]


def nquant(a, q, method):
    try:
        return float(np.quantile(a, q, method=method))
    except TypeError:                     # numpy < 1.22
        return float(np.quantile(a, q))


def sens_at_thr(y, s, thr, ge=True):
    pos = np.asarray(y) == 1
    hit = (np.asarray(s, dtype=float) >= thr) if ge else (np.asarray(s, dtype=float) > thr)
    return float((hit & pos).sum()) / int(pos.sum()) if int(pos.sum()) else float("nan")


def spec_at_thr(y, s, thr):
    neg = np.asarray(y) == 0
    if not int(neg.sum()):
        return float("nan")
    return float(((np.asarray(s, dtype=float) < thr) & neg).sum()) / int(neg.sum())


def ref_operating_point(y, s, spec=SPEC):
    """(a) the definition (CORRECTION_SPEC CS-07 v1.3), evaluated EXHAUSTIVELY over the achievable
    thresholds.  Returns (thr, achieved_specificity, sensitivity), or nan when the statistic is
    undefined (single-class y / no negatives) or when no threshold reaches the requested
    specificity.

    Target definition: among all ACHIEVABLE thresholds choose the one maximising sensitivity
    subject to achieved specificity >= spec; deterministic tie-break if tied.  Implemented here
    as the exhaustive search the definition names:

      * candidate thresholds = every observed score value (the ROC operating points at which the
        confusion counts can change) PLUS the boundary just above max(s) (the "predict nothing
        positive" point; sklearn.metrics.roc_curve reports the same point as its inf threshold);
      * achieved specificity at t = #{neg : s < t} / n_neg (predict positive iff s >= t);
      * among candidates with achieved specificity >= spec, take the maximum sensitivity;
      * tie-break (predefined, deterministic): largest achieved specificity, then largest
        threshold (== the last such candidate in an ascending scan).

    The candidate set is exhaustive over achievable operating points: any real threshold t has
    the same (achieved_spec, sens) pair as the smallest observed score >= t (or as the boundary
    point when t > max(s)), because no observed value lies strictly between them.  Hence the
    maximum found here is the exact definitional value; equivalently the selected threshold sits
    just above the k-th smallest negative score, k = ceil(spec * n_neg) (the smallest
    spec-compliant real threshold).  A score-value-only search WITHOUT the boundary point would
    be incomplete whenever no observed value reaches the requested specificity (e.g. all-tied
    scores, or one negative at the top): measured -- nan instead of the defined value 0.0 on
    2/300 fuzz datasets.  SENS-09 cross-checks this reference against
    sklearn.metrics.roc_curve's own ROC point selection.
    """
    y, s = mask_pairs(y, s)
    if len(np.unique(y)) < 2:
        return float("nan"), float("nan"), float("nan")
    neg = (y == 0)
    pos = (y == 1)
    n_neg, n_pos = int(neg.sum()), int(pos.sum())
    if n_neg == 0 or n_pos == 0:
        return float("nan"), float("nan"), float("nan")
    cands = [float(v) for v in np.unique(s)]
    cands.append(float(np.nextafter(cands[-1], np.inf)))   # boundary just above max(s)
    if spec <= 0:                                          # -inf predicts every row positive
        cands.insert(0, -np.inf)
    best = None                                            # (tp, tn), threshold
    for t in cands:
        tn = int(((s < t) & neg).sum())
        if tn / n_neg < spec - 1e-12:
            continue
        tp = int(((s >= t) & pos).sum())
        key = (tp, tn)
        if best is None or key >= best[0]:                 # ascending scan -> ties keep the last
            best = (key, t)
    if best is None:                                       # no threshold reaches this specificity
        return float("nan"), float("nan"), float("nan")
    return float(best[1]), best[0][1] / n_neg, best[0][0] / n_pos


def candidate(y, s, side, method, ge, spec=SPEC):
    """One point of the (side x interpolation-method x comparison) convention grid."""
    y, s = mask_pairs(y, s)
    neg = s[y == 0]
    if neg.size == 0:
        return float("nan"), float("nan")
    thr = nquant(neg, (1.0 - spec) if side == "1-spec" else spec, method)
    return sens_at_thr(y, s, thr, ge), thr


def one_token_spec(y, s, method="linear", nudge=False, spec=SPEC):
    """The one-token fix probe: threshold kept at the spec-quantile of the negatives (the same
    interpolation machinery the copies use), optionally nudged one ulp up.  Returns
    (delivered_achieved_specificity, returned_value) -- NOT a reference, only a probe."""
    y, s = mask_pairs(y, s)
    neg = s[y == 0]
    if neg.size == 0:
        return float("nan"), float("nan")
    thr = nquant(neg, spec, method)
    if nudge:
        thr = float(np.nextafter(thr, np.inf))
    return spec_at_thr(y, s, thr), sens_at_thr(y, s, thr)


def call(ref, y, s, spec=SPEC):
    """Execute a copy for DISPLAY; never raises.  6 decimals -> never used for comparison."""
    try:
        v = float(COPIES[ref](y, s, spec))
        return "nan" if v != v else "%.6f" % v
    except Exception as e:
        return type(e).__name__


def value(ref, y, s, spec=SPEC):
    """Execute a copy and return the exact float; an exception becomes nan."""
    try:
        return float(COPIES[ref](y, s, spec))
    except Exception:
        return float("nan")


def key_of(ref, y, s, spec=SPEC):
    """Exact comparison key: repr(float) (round-trip exact), 'nan', or the exception type."""
    try:
        v = float(COPIES[ref](y, s, spec))
    except Exception as e:
        return type(e).__name__
    return "nan" if v != v else repr(v)


_CELLS = None


def convention_cells():
    """Exact-match census of the (side x interpolation x comparison) grid against the copies on
    the fuzz corpus: {cell: n_datasets_reproduced}.  Cached.  This is how the implemented
    convention is IDENTIFIED by measurement instead of being assumed from the production text."""
    global _CELLS
    if _CELLS is None:
        _CELLS = {(side, m, ge): 0 for side in ("1-spec", "spec")
                  for m in QMETHODS for ge in (True, False)}
        if PRESENT:
            for y, s in fuzz():
                ref_v = key_of(PRESENT[0], y, s)
                if any(key_of(r, y, s) != ref_v for r in PRESENT):
                    continue                     # cross-copy divergence is SENS-02's finding
                v = value(PRESENT[0], y, s)
                for key in _CELLS:
                    vv, _ = candidate(y, s, key[0], key[1], key[2])
                    if abs(vv - v) < 1e-12:
                        _CELLS[key] += 1
    return _CELLS


def unique_cell():
    """The single grid cell that reproduces every copy on every fuzz dataset, else None."""
    hits = convention_cells()
    full = [k for k, c in hits.items() if c == len(fuzz())]
    return full[0] if len(full) == 1 else None


def delivered(y, s, spec=SPEC):
    """The operating point the copies deliver, as identified (never assumed) by unique_cell():
    (achieved_specificity, sensitivity, threshold), or None when their convention lies outside
    the quantile grid (e.g. after a correct fix to an exact-threshold implementation)."""
    cell = unique_cell()
    if cell is None:
        return None
    sens, thr = candidate(y, s, cell[0], cell[1], cell[2], spec)
    return spec_at_thr(y, s, thr), sens, thr


# ------------------------------------------------------------------ SENS-01
def t01():
    obs, ok = [], True
    for ref in LIVE:
        if ref in ABSENT:
            obs.append("%s ABSENT(%s)" % (ref, ABSENT[ref]))
            continue
        v = call(ref, Y_MAIN, S_MAIN)
        obs.append("%s -> %s" % (ref, v))
        try:
            vv = float(COPIES[ref](Y_MAIN, S_MAIN, SPEC))
            ok = ok and bool(np.isfinite(vv) and 0.0 <= vv <= 1.0)
        except Exception:
            ok = False
    if not PRESENT:
        ok = False
        obs.append("no live copy could be loaded")
    if ABSENT:
        obs.append("absent-in-this-tree: " + ", ".join(sorted(ABSENT)))
    return ok, "; ".join(obs)


L.guard("SENS-01", t01,
        "every live sens_at_spec copy is extracted VERBATIM from the measured tree and returns a "
        "finite value in [0,1] on the shared fixture (MAIN: 10 negatives, 5 positives, spec=0.80)",
        COPY_NOTE)


# ------------------------------------------------------------------ SENS-02
def t02():
    vals = {ref: call(ref, Y_MAIN, S_MAIN) for ref in PRESENT}
    spread = len(set(vals.values()))
    fz = fuzz()
    ndiv, exemplar = 0, ""
    for i, (y, s) in enumerate(fz):
        v = {ref: key_of(ref, y, s) for ref in PRESENT}
        if len(set(v.values())) != 1:
            ndiv += 1
            if not exemplar:
                exemplar = " dataset#%d %s" % (i, v)
    ok = (spread == 1) and (ndiv == 0)
    obs = ("MAIN copies: " + ", ".join("%s=%s" % (SHORT[k], v) for k, v in vals.items())
           + " (distinct values = %d)" % spread
           + "; fuzz %d datasets (binary y, ties, both classes): divergent = %d%s"
           % (len(fz), ndiv, exemplar))
    return ok, obs


L.guard("SENS-02", t02,
        "cross-copy agreement, EXACT: the same statistic must not have two answers depending on "
        "which script computed it (identical inputs -> bit-identical output on every live copy)",
        "parity gate, not a correctness oracle: agreement does not make a wrong value right, it "
        "only proves the three forks are the same implementation")


# ------------------------------------------------------------------ SENS-03
def t03():
    rows, ok = [], True
    for nm, y, s in (("MAIN", Y_MAIN, S_MAIN), ("BOUND", Y_BOUND, S_BOUND), ("TIE", Y_TIE, S_TIE)):
        thr_d, sp_d, sn_d = ref_operating_point(y, s)
        v_prod = value(PRESENT[0], y, s)
        prod_ok = all(key_of(r, y, s) == key_of(PRESENT[0], y, s) for r in PRESENT)
        ok = ok and prod_ok and (abs(v_prod - sn_d) < 1e-12)
        dl = delivered(y, s)
        dlt = ("delivered thr=%.6f achieved_spec=%.4f" % (dl[2], dl[0])) if dl else \
              "delivered operating point outside the quantile grid"
        rows.append("%s prod=%s def=%s (thr_def=%.6f achieved_spec_def=%.4f; %s)"
                    % (nm, call(PRESENT[0], y, s),
                       "nan" if sn_d != sn_d else "%.6f" % sn_d, thr_d, sp_d, dlt))
    # definitional equality is also counted over the whole fuzz corpus: a copy that matches the
    # definition on three small fixtures but not on tie-heavy data must not be certified
    fz = fuzz()
    n_eq_fuzz = 0
    for y, s in fz:
        _, _, sn_d = ref_operating_point(y, s)
        v = value(PRESENT[0], y, s)
        if np.isfinite(sn_d) and np.isfinite(v) and abs(v - sn_d) < 1e-12:
            n_eq_fuzz += 1
    ok = ok and (n_eq_fuzz == len(fz))
    # convention probe: which grid cell reproduces the copies on every fuzz dataset
    hits = convention_cells()
    n = len(fuzz())
    allmatch = [k for k, c in hits.items() if c == n]
    best_opposite = max([c for k, c in hits.items() if k[0] == "spec"] or [0])
    obs = ("; ".join(rows)
           + "; copy == definitional value on %d/%d fuzz datasets" % (n_eq_fuzz, n)
           + "; convention probe over %d fuzz datasets: cells reproducing production on ALL = %s"
           % (n, ", ".join("%s/%s/%s" % k for k in allmatch) or "none")
           + "; best (1-spec)-vs-spec opposite-side cell = %d/%d" % (best_opposite, n))
    return ok, obs


L.guard("SENS-03", t03,
        "definition (a), CORRECTION_SPEC CS-07 v1.3: sensitivity at fixed specificity 0.80 = TPR "
        "at the point chosen by 'among all ACHIEVABLE thresholds, the one maximising sensitivity "
        "subject to achieved specificity >= 0.80; deterministic tie-break if tied'.  The "
        "reference implements it as an exhaustive search over the achievable thresholds (observed "
        "score values + the boundary above max(s)) with achieved specificity #{neg: s < t}/n_neg "
        "and sensitivity #{pos: s >= t}/n_pos [standard ROC operating point, Fawcett 2006]; "
        "equivalently the threshold just above the k-th smallest negative score, "
        "k = ceil(0.80 * n_neg).  Methods_v1.md §14: 'Harmful-transfer discrimination was "
        "summarized by AUROC, with sensitivity at a prespecified fixed specificity of 0.80.'",
        "the production convention is NEVER adopted as the definition: the (1-spec) quantile is "
        "the COMPLEMENTARY operating point (specificity ~0.20), and the probe identifies which "
        "grid cell (side x interpolation x comparison) the code actually implements.  The "
        "reference search is exhaustive over achievable operating points (a pure score-value "
        "search without the boundary point would be incomplete whenever no observed value "
        "reaches the specificity, e.g. all-tied scores) and is cross-checked against "
        "sklearn.metrics.roc_curve's point selection in SENS-09")


# ------------------------------------------------------------------ SENS-04
def t04():
    # (i) inclusivity is MEASURED, not assumed: on datasets where the identified cell's '>=' and
    #     '>' variants differ, check which one the copies reproduce.
    idc = unique_cell()
    fz = fuzz()
    disc, inc_hit = 0, 0
    if idc is not None:
        for y, s in fz:
            v = value(PRESENT[0], y, s)
            veq, _ = candidate(y, s, idc[0], idc[1], True)
            vgt, _ = candidate(y, s, idc[0], idc[1], False)
            if abs(veq - vgt) > 1e-12:
                disc += 1
                if abs(veq - v) < 1e-12:
                    inc_hit += 1
    # (ii) black-box boundary contract: whatever operating point a copy delivers, its returned
    #      value cannot exceed the maximum sensitivity attainable at achieved specificity >= SPEC
    #      (the definitional value).  The assertion depends ONLY on copy-produced values, so any
    #      spec-compliant implementation passes by construction and the test can go GREEN.
    #      The same loop also measures the one-token fix (iii) so the fuzz corpus is walked once.
    rows, viol = [], []
    for nm, y, s in (("MAIN", Y_MAIN, S_MAIN), ("BOUND", Y_BOUND, S_BOUND), ("TIE", Y_TIE, S_TIE)):
        _, sp_d, sn_d = ref_operating_point(y, s)
        v_max = max(value(r, y, s) for r in PRESENT)
        over = v_max - sn_d
        # a NaN must not silence the boundary check: these fixtures have both classes and finite
        # scores, so a non-finite copy value is a violation in its own right
        if (not np.isfinite(v_max)) or v_max > sn_d + 1e-12:
            viol.append(nm)
        rows.append("%s copy=%s max_attainable@spec>=%.2f=%.6f over_max=%s"
                    % (nm, ("%.6f" % v_max) if np.isfinite(v_max) else "nan", SPEC, sn_d,
                       ("%+.4f" % over) if np.isfinite(over) else "n/a (non-finite copy)"))
    n_viol_fuzz = n_naive_short = n_nudge_short = n_naive_valbad = n_nudge_valbad = 0
    n_nudge_overshoot = 0
    for y, s in fz:
        _, _, sn_d = ref_operating_point(y, s)
        vmax = max(value(r, y, s) for r in PRESENT)
        # fail closed: the fuzz datasets all have both classes and finite scores, so a correct
        # copy always returns a finite value -- NaN is a violation, not an absence of violation
        if (not np.isfinite(vmax)) or vmax > sn_d + 1e-12:
            n_viol_fuzz += 1
        sp_n, v_n = one_token_spec(y, s, "linear", nudge=False)
        sp_u, v_u = one_token_spec(y, s, "linear", nudge=True)
        if sp_n < SPEC - 1e-12:
            n_naive_short += 1
        if sp_u < SPEC - 1e-12:
            n_nudge_short += 1
        if not (np.isfinite(v_n) and np.isfinite(sn_d) and abs(v_n - sn_d) < 1e-12):
            n_naive_valbad += 1
        if not (np.isfinite(v_u) and np.isfinite(sn_d) and abs(v_u - sn_d) < 1e-12):
            n_nudge_valbad += 1
        if np.isfinite(v_u) and np.isfinite(sn_d) and v_u > sn_d + 1e-12:
            n_nudge_overshoot += 1
    # (iii) the one-token fix (1-spec -> spec, threshold kept at the interpolated spec-quantile) is
    #      measured, not assumed: on BOUND the quantile lands exactly ON the negative score 0.70
    #      and, with the copies' comparison s >= thr, that negative is predicted positive, so the
    #      DELIVERED achieved specificity is 4/6 = 0.6667 < 0.80 -- a specificity SHORTFALL, not a
    #      sensitivity overshoot.  The value it returns (0.50) coincides with the definitional
    #      value on this fixture, so only the delivered-specificity probe exposes it here; the
    #      shortfall is therefore counted directly over the fuzz corpus (naive and nextafter-
    #      nudged variants).  A corrected copy must move the threshold just ABOVE the k-th
    #      smallest negative score, k = ceil(spec * n_neg).
    naive = {}
    for m in QMETHODS:
        _, t_nv = candidate(Y_BOUND, S_BOUND, "spec", m, True)
        naive[m] = (spec_at_thr(Y_BOUND, S_BOUND, t_nv), sens_at_thr(Y_BOUND, S_BOUND, t_nv))
    n_below = sum(1 for m in QMETHODS if naive[m][0] < SPEC - 1e-12)
    dl = delivered(Y_BOUND, S_BOUND)
    dl_txt = ("BOUND delivered operating point: thr=%.6f achieved_spec=%.4f sens=%.4f"
              % (dl[2], dl[0], dl[1])) if dl else \
             "BOUND delivered operating point: outside the quantile grid"
    inc_txt = ("inclusive s>=thr measured against cell %s/%s/%s: on %d/%d fuzz datasets the >= "
               "and > variants differ and the copies matched >= on %d/%d"
               % (idc[0], idc[1], idc[2], disc, len(fz), inc_hit, disc)) if idc is not None else \
              "inclusivity not measurable: no grid cell reproduces the copies on the fuzz corpus"
    ok = (not viol) and (n_viol_fuzz == 0)
    obs = ("%s; %s; %s; fuzz boundary violation (copy value non-finite or above the max "
           "attainable at spec>=0.80) on %d/%d datasets; fixture violations: %s; one-token fix (spec-quantile "
           "threshold, same comparison) DELIVERS achieved specificity < 0.80 on BOUND for %d/%d "
           "interpolation methods [linear: delivered_spec=%.4f value=%.4f == definitional %.4f, so "
           "this fixture is a specificity SHORTFALL, not a sensitivity overshoot], and over the "
           "fuzz corpus on %d/%d datasets (its nextafter-nudged variant still on %d/%d, its value "
           "above the attainable maximum on %d/%d and differing from the definitional value on "
           "%d/%d vs %d/%d for the naive variant)"
           % ("; ".join(rows), inc_txt, dl_txt, n_viol_fuzz, len(fz),
              ", ".join(viol) or "none", n_below, len(QMETHODS),
              naive["linear"][0], naive["linear"][1], ref_operating_point(Y_BOUND, S_BOUND)[2],
              n_naive_short, len(fz), n_nudge_short, len(fz), n_nudge_overshoot, len(fz),
              n_nudge_valbad, len(fz), n_naive_valbad, len(fz)))
    return ok, obs


L.guard("SENS-04", t04,
        "convention layer (b) + boundary: a value labelled sensitivity-at-specificity-0.80 must "
        "be attainable at a delivered specificity >= 0.80; no spec-compliant operating point can "
        "exceed the definitional value (CS-07 v1.3: the maximum sensitivity over the ACHIEVABLE "
        "thresholds -- observed scores + the boundary above max(s) -- with achieved specificity "
        "#{neg: s < thr}/n_neg >= 0.80), so a copy-produced value above that maximum proves a "
        "sub-spec operating point, and a non-finite copy value is a violation in itself. "
        "Protocol v1.0 FROZEN §7.2 line 248: 'sensitivity at fixed specificity（specificity 固定为 "
        "0.80，预先写死）'; Methods_v1.md line 412 §14: 'sensitivity at a prespecified fixed "
        "specificity of 0.80'",
        "the assertion consumes ONLY the copies' returned values (green-able on a corrected "
        "implementation; a NaN copy value fails closed on BOTH arms, the fixtures and the fuzz "
        "corpus); inclusivity and the delivered operating point are identified by measurement, "
        "never assumed; the one-token fix 1-spec -> spec is probed and shown to leave a "
        "DELIVERED-specificity shortfall (0.6667 < 0.80 on BOUND; measured over the fuzz corpus "
        "for the naive and nextafter-nudged variants) -- the corrected implementation must "
        "satisfy the CS-07 target definition (equivalently: threshold just above the k-th "
        "smallest negative score, k = ceil(spec * n_neg))")


# ------------------------------------------------------------------ SENS-05
def t05():
    yv, sv = mask_pairs(Y_NAN, S_NAN)
    rows, ok = [], True
    for ref in PRESENT:
        v_raw = call(ref, Y_NAN, S_NAN)
        v_msk = call(ref, yv, sv)
        same = (key_of(ref, Y_NAN, S_NAN) == key_of(ref, yv, sv))
        ok = ok and same
        rows.append("%s raw=%s masked=%s %s" % (SHORT[ref], v_raw, v_msk,
                                                "OK" if same else "MISALIGNED"))
    thr_d, sp_d, sn_d = ref_operating_point(yv, sv)
    _, _, sn_roll = ref_operating_point(np.roll(yv, 1), sv)
    power = abs(sn_d - sn_roll) > 1e-12
    ok = ok and power
    obs = ("%d rows -> %d after the shared mask; %s; definitional reference on the paired mask = "
           "%.6f (thr=%.6f, achieved_spec=%.4f) vs %.6f for a one-step rolled label vector -> "
           "misalignment would change the answer (power=%s); y-NaN note: a NaN in y is NOT filtered "
           "(row is silently dropped from both numerator and denominator)"
           % (len(Y_NAN), len(sv), "; ".join(rows), sn_d, thr_d, sp_d, sn_roll, power))
    return ok, obs


L.guard("SENS-05", t05,
        "the mask m = isfinite(s) must subset BOTH y and s (paired filtering); filtering one array "
        "only would silently re-pair labels with the wrong scores and inflate sensitivity",
        "the NaN fixture drops two rows carrying DIFFERENT labels, and the rolled-label probe "
        "proves the reference can see a mis-pairing (it changes the value)")


# ------------------------------------------------------------------ SENS-06
def t06():
    deg = (
        ("single_class_all0", np.array([0, 0, 0, 0], float), np.array([0.1, 0.2, 0.3, 0.4])),
        ("single_class_all1", np.array([1, 1, 1], float), np.array([0.1, 0.2, 0.3])),
        ("two_labels_no_zero", np.array([1, 2, 2, 1], float), np.array([0.1, 0.9, 0.5, 0.3])),
        ("all_tied", np.array([0, 1, 0, 1, 1], float), np.array([0.5] * 5)),
        ("empty", np.array([], float), np.array([], float)),
        ("all_nan_scores", np.array([0, 1, 0, 1], float), np.array([np.nan] * 4)),
    )
    rows, ok = [], True
    for nm, y, s in deg:
        v = {SHORT[ref]: key_of(ref, y, s) for ref in PRESENT}
        agree = len(set(v.values())) == 1
        # a returned 'nan' is NOT agreement where the statistic is DEFINED: all_tied carries both
        # classes and negatives, so its value is defined (0.0); only genuine divergence or a nan
        # VALUE fails this row (same-exception outcomes stay permitted by the documented contract)
        nan_defined = nm == "all_tied" and any(vv == "nan" for vv in v.values())
        ok = ok and agree and not nan_defined
        extra = ""
        if nm == "all_tied":
            dl = delivered(y, s)
            extra = " delivered achieved_spec=%.4f" % spec_at_thr(y, s, dl[2] if dl else 0.5)
            if nan_defined:
                extra += " [statistic is DEFINED here: nan is not agreement]"
        rows.append("%s -> %s%s%s" % (nm, v, extra,
                                      "" if (agree and not nan_defined) else " DIVERGE"))
    # a zero-positive-class input must not fabricate a number: sensitivity is undefined
    v0 = {SHORT[ref]: key_of(ref, np.array([0, 0, 0, 0], float), np.array([0.1, 0.2, 0.3, 0.4]))
          for ref in PRESENT}
    fab = [k for k, v in v0.items() if v != "nan"]
    ok = ok and not fab
    obs = ("; ".join(rows) + "; fabricated values on a zero-positive-class input (must be undefined "
           "nan, not 0.0): %s; nan-is-not-agreement enforced where the statistic is defined "
           "(all_tied must return the defined finite value 0.0)"
           % (", ".join("%s=%s" % (k, v0[k]) for k in fab) or "none"))
    return ok, obs


L.guard("SENS-06", t06,
        "degenerate-input contract: every live copy must return the same outcome (same value, or "
        "the same exception type) on every degenerate input, and a single-class y must not yield a "
        "fabricated sensitivity (with zero positive cases the statistic is undefined, not 0.0)",
        "single-class y / empty negatives / all-tied scores are outside the frozen corpus (65 "
        "harmful of 202) but reachable from any re-aggregation of it; DIVERGE rows are latent "
        "cross-copy divergence, not a value error on the frozen column")


# ------------------------------------------------------------------ SENS-08
EXP08 = ("definition (a) on the FROZEN corpus: a value reported as sens_at_spec80 must be "
         "attainable at a delivered specificity >= 0.80, i.e. it cannot exceed the definitional "
         "sensitivity at specificity 0.80 (CS-07 v1.3: the maximum sensitivity over the "
         "ACHIEVABLE thresholds -- observed scores + the boundary above max(s) -- with achieved "
         "specificity >= 0.80); a non-finite recon fails closed. Methods_v1.md line 412 §14: "
         "'Harmful-transfer discrimination was summarized by AUROC, with sensitivity at a "
         "prespecified fixed specificity of 0.80.'; Protocol v1.0 FROZEN "
         "§7.2 line 248: 'sensitivity at fixed specificity（specificity 固定为 0.80，预先写死）'")

EXP08B = ("characterisation (current-state pin, warn_only by design): the frozen "
          "p6d_signal_eval.csv column sens_at_spec80 (copied at p11_figure_datapack.py:147 into "
          "fig2_forest.csv) is reproduced from the frozen CSV inputs by the live copy to < 1e-9, "
          "which pins the measured object to the manuscript column.  A correct fix changes the "
          "copies, so the OLD frozen column is no longer reproduced until the pipeline is re-run "
          "-> this line is expected PASS -> WARN in PATCHED mode; it is a pin, not a test.")


def _impact():
    need = ["p6d_outcome.csv", "p6b_label_free_signals_blinded.csv", "p6c_L3_blinded.csv",
            "p6d_signal_eval.csv"]
    miss = [f for f in need if not os.path.exists(os.path.join(L.OUT, f))]
    if miss:
        return None
    O = pd.read_csv(os.path.join(L.OUT, "p6d_outcome.csv"))
    SIG = pd.read_csv(os.path.join(L.OUT, "p6b_label_free_signals_blinded.csv"))
    L3 = pd.read_csv(os.path.join(L.OUT, "p6c_L3_blinded.csv"))
    E = O.merge(SIG[["task_id", "domain_auc", "mmd", "psi_mean", "psi_max", "js_mean", "js_max",
                     "pred_mean_shift", "entropy_mean_shift", "entropy_js", "conf_dist_js",
                     "pred_pos_prev_shift"]], on="task_id", how="left") \
         .merge(L3[["task_id", "SUDO_AURCC", "AUROC_IOW", "AUROC_OM", "AUROC_DR"]],
                on="task_id", how="left")
    METHODS = [("domain_auc", lambda x: x["domain_auc"]), ("mmd", lambda x: x["mmd"]),
               ("psi_mean", lambda x: x["psi_mean"]), ("psi_max", lambda x: x["psi_max"]),
               ("js_mean", lambda x: x["js_mean"]), ("js_max", lambda x: x["js_max"]),
               ("pred_mean_shift", lambda x: x["pred_mean_shift"].abs()),
               ("entropy_mean_shift", lambda x: x["entropy_mean_shift"].abs()),
               ("entropy_js", lambda x: x["entropy_js"]),
               ("conf_dist_js", lambda x: x["conf_dist_js"]),
               ("pred_pos_prev_shift", lambda x: x["pred_pos_prev_shift"].abs()),
               ("SUDO_AURCC", lambda x: -x["SUDO_AURCC"]),
               ("AUROC_DR", lambda x: -x["AUROC_DR"]), ("AUROC_IOW", lambda x: -x["AUROC_IOW"]),
               ("AUROC_OM", lambda x: -x["AUROC_OM"])]
    frozen = pd.read_csv(os.path.join(L.OUT, "p6d_signal_eval.csv")).set_index("method")
    ref = REF_P6D if REF_P6D in COPIES else PRESENT[0]
    y = E["harmful"].values.astype(int)
    rec, maxdiff, maxover, n_bad = [], 0.0, -np.inf, 0
    specs, rep, dfn, pairs = [], [], [], []
    for name, fn in METHODS:
        r = fn(E).values.astype(float)
        okm = np.isfinite(r)
        yv, sv = y[okm], r[okm]
        fz = float(frozen.loc[name, "sens_at_spec80"])
        v = float(COPIES[ref](yv, sv, SPEC))
        _, _, sn_d = ref_operating_point(yv, sv)
        # fail closed on a non-finite recon: a NaN must never read as 'reproduced to 0.0e+00'
        # (maxdiff) or 'max overshoot = -inf' (maxover) -- Python max keeps the running value
        # on a False nan comparison, so the update must be gated, not merely folded in
        fin = bool(np.isfinite(v) and np.isfinite(fz) and np.isfinite(sn_d))
        if fin:
            maxdiff = max(maxdiff, abs(v - fz))
            maxover = max(maxover, v - sn_d)
        else:
            n_bad += 1
        dl = delivered(yv, sv)
        specs.append(dl[0] if dl else float("nan"))
        rep.append(v); dfn.append(sn_d)
        pairs.append((name, yv, sv))
        rec.append("%s: frozen=%.4f recon=%.4f delivered_spec=%s def_sens=%.4f over_max=%s"
                   % (name, fz, v, ("%.4f" % dl[0]) if dl else "grid-n/a", sn_d,
                      ("%+.4f" % (v - sn_d)) if fin else "n/a (non-finite)"))
    return dict(E=E, y=y, nneg=int((y == 0).sum()), rec=rec, maxdiff=maxdiff, maxover=maxover,
                n_bad=n_bad, specs=specs, rep=rep, dfn=dfn, nm=len(METHODS), pairs=pairs)


R = None
try:
    R = _impact()
    if R is None:
        L.emit("SENS-08", "WARN", EXP08, "input absent from %s" % L.OUT,
               "frozen inputs unavailable: value not assertable")
        L.emit("SENS-08b", "WARN", EXP08B, "input absent from %s" % L.OUT,
               "characterisation not computable")
    else:
        s_fin = [x for x in R["specs"] if x == x]
        sp_txt = ("%.4f..%.4f" % (min(s_fin), max(s_fin))) if s_fin else "grid-n/a"
        rep_fin = [x for x in R["rep"] if np.isfinite(x)]
        dfn_fin = [x for x in R["dfn"] if np.isfinite(x)]
        rep_txt = ("%.4f..%.4f" % (min(rep_fin), max(rep_fin))) if rep_fin else "all non-finite"
        dfn_txt = ("%.4f..%.4f" % (min(dfn_fin), max(dfn_fin))) if dfn_fin else "all non-finite"
        if R["n_bad"] == 0:
            over_txt = "max overshoot = %+.4f" % R["maxover"]
        elif R["n_bad"] == R["nm"]:
            over_txt = "max overshoot = n/a: all %d method reconstructions non-finite" % R["nm"]
        else:
            over_txt = ("max overshoot over the %d finite method(s) = %+.4f; NON-FINITE recon "
                        "or definition on %d/%d methods"
                        % (R["nm"] - R["n_bad"], R["maxover"], R["n_bad"], R["nm"]))
        rep_note = ("the live copy reproduces the old frozen column to max abs diff %.1e (pin "
                    "line SENS-08b)" % R["maxdiff"]) if R["n_bad"] == 0 else \
                   ("old-column reproduction is NOT computable: %d/%d method recon(s) non-finite "
                    "(pin line SENS-08b)" % (R["n_bad"], R["nm"]))
        obs8 = ("E reconstructed as p6d does (%d tasks, %d harmful / %d non-harmful); %s; "
                "copy values vs the definitional maximum attainable at specificity >= %.2f: %s; "
                "reported sens_at_spec80 range %s vs definitional %s; delivered specificity under "
                "the identified grid convention = %s; per-method: %s"
                % (len(R["E"]), int(R["y"].sum()), R["nneg"], rep_note, SPEC, over_txt,
                   rep_txt, dfn_txt, sp_txt, "; ".join(R["rec"])))
        L.check("SENS-08", R["maxover"] <= 1e-12 and R["n_bad"] == 0, EXP08, obs8,
                "the assertion consumes ONLY copy-produced values, so a spec-compliant "
                "implementation passes by construction, and a non-finite recon fails closed (a "
                "NaN cannot read as 'max overshoot = -inf'); the frozen artefact is regenerated "
                "only by re-running p6d, so the pin (bit-reproduction of the old column) is split "
                "out into SENS-08b as a warn_only characterisation")
        b_txt = ("max abs(frozen - recon) = %.1e over the %d finite method recon(s); non-finite "
                 "on %d/%d" % (R["maxdiff"], R["nm"] - R["n_bad"], R["n_bad"], R["nm"])) \
            if R["n_bad"] else ("max abs(frozen - recon) = %.1e over %d methods"
                                % (R["maxdiff"], R["nm"]))
        L.check("SENS-08b", R["maxdiff"] < 1e-9 and R["n_bad"] == 0, EXP08B, b_txt,
                "warn_only=True by design: a current-state characterisation (pin); PASS -> WARN "
                "on a corrected tree is expected behaviour, not a defect; a non-finite recon "
                "cannot claim a reproduction", warn_only=True)
except Exception as _e:
    L.emit("SENS-08", "FAIL", EXP08, "EXCEPTION %s: %s" % (type(_e).__name__, _e),
           "impact test could not be computed")
    L.emit("SENS-08b", "FAIL", EXP08B, "EXCEPTION %s: %s" % (type(_e).__name__, _e),
           "characterisation could not be computed")


# ------------------------------------------------------------------ SENS-09
def corrected_sens(y, s, spec=SPEC):
    """Definition-only SYNTHETIC implementation, used ONLY by SENS-09 to demonstrate that the
    acceptance criterion of this module is satisfiable (green-ability self-test).  It is NOT a
    reference for production values (the reference is ref_operating_point, which evaluates the
    CS-07 target definition by exhaustive enumeration over achievable thresholds) and it is not a
    production copy.

    CS-07 v1.3 target definition: among all achievable thresholds choose the one maximising
    sensitivity subject to achieved specificity >= spec.  Equivalently (same value, same
    operating point): the smallest real threshold with achieved specificity >= spec -- i.e. just
    above the k-th smallest negative score, k = ceil(spec * n_neg); predict positive as s >= thr.
    """
    y = np.asarray(y)
    s = np.asarray(s, dtype=float)
    m = np.isfinite(s)
    y, s = y[m], s[m]
    if len(np.unique(y)) < 2:
        return float("nan")
    neg = np.sort(s[y == 0])
    if neg.size == 0:
        return float("nan")
    k = int(np.ceil(spec * neg.size - 1e-12))
    if k > neg.size:                       # no threshold can reach this specificity
        return float("nan")
    thr = -np.inf if k < 1 else float(np.nextafter(neg[k - 1], np.inf))
    pos = y == 1
    return float(((s >= thr) & pos).sum()) / int(pos.sum())


def sklearn_operating_point(y, s, spec=SPEC):
    """External cross-check oracle (NOT the reference): sklearn.metrics.roc_curve's own ROC
    operating points, evaluated with the CS-07 target rule and the same deterministic tie-break
    as ref_operating_point (maximise sensitivity subject to achieved specificity >= spec; ties ->
    largest achieved specificity, then largest threshold).  Returns (thr, achieved_spec,
    sensitivity), a nan triple when the statistic is undefined, or None when sklearn is not
    importable.  roc_curve's threshold set is the observed score values plus its boundary point
    (thresholds[0] = inf), so this is the external implementation of the same achievable set."""
    try:
        from sklearn.metrics import roc_curve
    except Exception:
        return None
    y, s = mask_pairs(y, s)
    if len(np.unique(y)) < 2:
        return float("nan"), float("nan"), float("nan")
    fpr, tpr, thr = roc_curve(y, s, drop_intermediate=False)
    spec_arr = 1.0 - fpr
    compl = spec_arr >= spec - 1e-12
    if not compl.any():
        return float("nan"), float("nan"), float("nan")
    best = None                                # key, threshold
    for i in np.where(compl)[0]:
        key = (float(tpr[i]), float(spec_arr[i]), float(thr[i]))
        if best is None or key >= best[0]:
            best = (key, float(thr[i]))
    return best[1], best[0][1], best[0][0]


EXP09 =("green-ability self-test (acceptance-criterion satisfiability) + reference audit: the "
         "predicates this module asserts on production -- SENS-03 equality with the definitional "
         "reference and SENS-04 never above the definitional maximum -- must be satisfiable by a "
         "mathematically correct implementation of the CS-07 v1.3 target definition 'among all "
         "achievable thresholds, the one maximising sensitivity subject to achieved specificity "
         ">= spec; deterministic tie-break if tied' (equivalently thr = just above the k-th "
         "smallest negative score, k = ceil(spec * n_neg); sensitivity = #{pos: s >= thr}/n_pos).  "
         "The reference search itself is cross-checked against sklearn.metrics.roc_curve's ROC "
         "operating-point selection.  The synthetic implementation used here probes the "
         "acceptance predicate only; it is never a reference and never a production copy")


def t09():
    cases = list(FIXTURES) + [
        ("single_class", np.array([1.0, 1.0, 1.0]), np.array([0.1, 0.2, 0.3])),
        ("all_tied", np.array([0.0, 1.0, 0.0, 1.0, 1.0]), np.array([0.5] * 5)),
    ]
    rows, bad = [], []
    for nm, y, s in cases:
        _, _, sn_d = ref_operating_point(y, s)
        v = corrected_sens(y, s)
        eq = (v != v and sn_d != sn_d) or (abs(v - sn_d) < 1e-12)
        if not eq:
            bad.append(nm)
        rows.append("%s def=%.6f synthetic=%.6f %s" % (nm, sn_d, v, "OK" if eq else "MISMATCH"))
    fz = fuzz()
    n_fuzz_bad = n_naive_bad = 0
    for y, s in fz:
        _, _, sn_d = ref_operating_point(y, s)
        c = corrected_sens(y, s)
        # fail closed: fuzz datasets always have both classes and a spec-compliant point, so a
        # correct synthetic value is finite -- a non-finite value is NOT equality
        if not (np.isfinite(c) and np.isfinite(sn_d) and abs(c - sn_d) < 1e-12):
            n_fuzz_bad += 1
        _, v_n = one_token_spec(y, s, "linear", nudge=False)
        if not (np.isfinite(v_n) and np.isfinite(sn_d) and abs(v_n - sn_d) < 1e-12):
            n_naive_bad += 1
    frozen_txt = "frozen corpus not available: frozen-corpus arm skipped"
    n_frozen_bad = None
    if R is not None and R.get("pairs"):
        n_frozen_bad, worst = 0, 0.0
        for _nmf, yv, sv in R["pairs"]:
            _, _, sn_d = ref_operating_point(yv, sv)
            c = corrected_sens(yv, sv)
            fin = bool(np.isfinite(c) and np.isfinite(sn_d))
            if fin:
                worst = max(worst, abs(c - sn_d))
            if not (fin and abs(c - sn_d) <= 1e-12):
                n_frozen_bad += 1
        frozen_txt = ("frozen corpus: the synthetic implementation reproduces the definitional "
                      "value on %d/%d method vectors (max abs diff = %.1e)"
                      % (len(R["pairs"]) - n_frozen_bad, len(R["pairs"]), worst))
    # reference audit (CS-07 acceptance criterion ③): the exhaustive search over achievable
    # thresholds must reproduce sklearn.metrics.roc_curve's own ROC operating-point selection on
    # every non-degenerate case; when sklearn is not importable the cross-check is disclosed as
    # skipped instead of silently passed
    sk_pool = list(FIXTURES) + [("fuzz#%d" % i, y, s) for i, (y, s) in enumerate(fz)] \
        + (list(R["pairs"]) if R is not None and R.get("pairs") else [])
    n_sk_cmp = n_sk_bad = 0
    sk_ok, sk_txt = True, "sklearn not importable: cross-check skipped (disclosed, not hidden)"
    if sklearn_operating_point(Y_MAIN, S_MAIN) is not None:
        for nm2, y2, s2 in sk_pool:
            t_k, sp_k, sn_k = sklearn_operating_point(y2, s2)
            _, sp_d, sn_d = ref_operating_point(y2, s2)
            n_sk_cmp += 1
            same = ((sp_k != sp_k and sp_d != sp_d) or
                    (np.isfinite(sp_k) and np.isfinite(sp_d) and
                     abs(sp_k - sp_d) < 1e-12 and abs(sn_k - sn_d) < 1e-12))
            if not same:
                n_sk_bad += 1
        sk_ok = (n_sk_bad == 0)
        sk_txt = ("sklearn.metrics.roc_curve point-selection cross-check over %d cases "
                  "(fixtures + fuzz + frozen method vectors): disagreements with the reference "
                  "= %d" % (n_sk_cmp, n_sk_bad))
    ok = (not bad) and n_fuzz_bad == 0 and n_naive_bad > 0 \
        and (n_frozen_bad is None or n_frozen_bad == 0) and sk_ok
    obs = ("%s; fuzz corpus (%d datasets): definitional equality mismatches = %d; the one-token "
           "fix is still rejected by the same predicate on %d/%d datasets (power check); %s; %s"
           % ("; ".join(rows), len(fz), n_fuzz_bad, n_naive_bad, len(fz), frozen_txt, sk_txt))
    return ok, obs


L.guard("SENS-09", t09, EXP09,
        "the synthetic implementation is definition-only and is compared against the reference on "
        "fixtures + fuzz + frozen corpus; the one-token 1-spec -> spec fix must FAIL the same "
        "predicate (otherwise the criterion would have no power); the reference is cross-checked "
        "against sklearn.metrics.roc_curve's own point selection and a wrong synthetic value/NaN "
        "fails closed; a corrected copy must implement the CS-07 target definition (equivalently, "
        "the threshold just above the k-th smallest negative score) -- nextafter(np.quantile(neg, "
        "spec)) is not sufficient on tie-heavy data (measured in SENS-04: delivered specificity "
        "< 0.80 on 132/300 fuzz datasets, value above the attainable maximum on 62/300)")

print("\n[t10_sens_at_spec done]", flush=True)
