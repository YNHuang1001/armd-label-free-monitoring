# -*- coding: utf-8 -*-
r"""t13_eval_core.py — EVS family EVALCORE: the un-closed core of the evaluation layer

WHY THIS FILE EXISTS
================================================================================
EVS v1.2 closed the estimator-kernel layer and part of the evaluation layer
(t12: p6f settlement, bootstrap mechanics, p12 calib, p6f LOO section).  A
completeness critique found four number-producing paths still executed by NO
test:

  (A) p11_figure_datapack.py — the WHOLE file: the summary statistics behind
      most figures (fig1c quantiles / P90-P10, fig1d per-direction medians,
      fig1e aligned medians + harmful counts, fig1b counts, fig3a category
      labels, fig3c decomposition, fig4d MAE + calibration + rho_ess, fig5/5b
      curves, figS6 counts, tableS1/S2 merges).
  (B) p6d_01_unseal_evaluation.py:76-219 — the first layer: 3-site model
      rebuild, auroc_target, d_auroc, harmful/stable -> p6d_outcome.csv, the
      backbone of Fig 1C/1D and of every downstream count (65/32.2%, 49).
  (C) p6d_01_unseal_evaluation.py:243-304 — the per-method loop: harmful_auc,
      the organism-drug cluster bootstrap CI, Spearman, AUPRC, grouped_cv_oof,
      in-file LOO, coverage -> Fig 2A/2B and tableS2.
  (D) p6e_01_lodo_fold_qc.py:93-145 — fit_predict / run / fold composition /
      macro_mean / fold quantiles; plus the documented inconsistency that
      p6e's macro_lo/hi is a SIZE-WEIGHTED bootstrap CI labelled against an
      UNWEIGHTED macro_mean.

HOW THE REFERENCES ARE BUILT (rule: definition or external source only)
--------------------------------------------------------------------------------
  (a) DEFINITION: d_auroc = auroc_target - auroc_source_internal;
      harmful = I(d_auroc <= -0.05); stable = I(|d_auroc| <= 0.02)
      (task brief; protocol v1.2 section 11; p6d's own frozen spec header).
  (b) DEFINITION: AUROC = Mann-Whitney rank statistic with MIDRANKS on ties;
      AUPRC = step-wise average precision sum_n (R_n - R_{n-1}) * P_n over
      distinct descending score thresholds (sklearn.metrics used as the
      external implementation cross-check).
  (c) DEFINITION: cluster bootstrap = resample CLUSTERS (organism-drug
      groups) WITH replacement, statistic recomputed over the resample
      RETAINING multiplicity (Efron & Tibshirani 1993, ch. 6; CS-08).
  (d) DEFINITION: leave-one-group-out = fit the logistic model on the other
      groups ONLY, predict the held-out group; macro = UNWEIGHTED mean over
      held-out groups (LODO definition; p6d:243-255 / p6e:93-112 semantics).
  (e) DEFINITION: quantile = linear interpolation on sorted values
      (Hyndman-Fan type 7, the numpy/pandas default); sample sd = ddof=1;
      MAE = mean|est - true| over finite estimates; calibration = group means
      of the actual value within the fixed 0.05-wide bins of the predicted
      delta (the documented calibration convention).
  (f) EXTERNAL: frozen artifacts (p6d_outcome.csv, p6d_signal_eval.csv,
      p6e_lodo_*.csv, figures/data/*).  Reconciliation against a frozen CSV is
      labelled as such; it is NOT the same as executing the code.
  The production scripts are NOT references; they are the objects under test,
  loaded VERBATIM from the tree under EVS_PROD_DIR (evs_load.grab / load_def)
  or box-extracted from the file text by content anchors.  Code is executed
  verbatim; only logging statements (p(...)) and the final to_csv writes to
  the PRODUCTION directory are replaced/redirected so the module stays
  read-only on the production tree.  Every exclusion is listed by block label
  in the module's census line (EC-15).

WHAT THIS MODULE MEASURES (16 RESULT lines)
--------------------------------------------------------------------------------
  EC-01  (B) p6d first layer EXECUTED VERBATIM (3-site STORE + models + loop)
             reproduces frozen p6d_outcome.csv                          PASS/FAIL
  EC-02  (B) frozen p6d_outcome.csv columns vs DEFINITIONS
             (identity, flags, one-model-per-source invariant)          PASS/FAIL
  EC-03  (C) p6d per-method loop EXECUTED VERBATIM reproduces frozen
             p6d_signal_eval.csv (all 12 columns)                       PASS/FAIL
  EC-04  (C) independent harmful-AUROC / AUPRC reference (own midrank
             implementation) vs frozen columns                          PASS/FAIL
  EC-05  (C) independent cluster-bootstrap replay (clusters, multiplicity)
             reproduces frozen sp_lo/sp_hi/auc_lo/auc_hi; set-dedup control
             shown to move the interval                                 PASS/FAIL
  EC-06  (A) fig3a: executed block == frozen fig3a_lodo.csv AND frozen
             categories == CI-vs-0.5 rule (15/15)                       PASS/FAIL
  EC-06b (A) WARN characterisation: 2/15 fig3a categories flip if p6e's
             size-weighted CI convention is used instead                WARN
  EC-07  (A) fig3c decomposition: executed == frozen; subset/aggregation
             identities; cross-artifact identity with p6e macro_mean     PASS/FAIL
  EC-08  (A) fig4d: executed == frozen; independent numpy MAE, scipy
             Spearman, calibration bin means, rho_ess                    PASS/FAIL
  EC-09  (A) fig1c summary: executed == frozen; independent np.quantile /
             mean / sd(ddof=1) / counts                                PASS/FAIL
  EC-10  (A) fig5/5b: executed == frozen; independent per-rep AUC replay
             (own implementation, all 15 methods x 6 level-rows)        PASS/FAIL
  EC-11  (A) remaining panels executed == frozen (fig1b, fig1d, fig1e,
             fig2, fig4e, figS6, tableS1, tableS2) + merge integrity;
             independent fig1e recomputation                           PASS/FAIL
  EC-12  (D) p6e fold composition + aggregation loop EXECUTED VERBATIM
             reproduces the three frozen p6e CSVs to <1e-15             PASS/FAIL
  EC-13  (D) p6e fit_predict / run semantics on oracle fixtures
             (single-class fold exclusion, unweighted macro, sd=0 and NaN
             handling, no y_test leakage)                              PASS/FAIL
  EC-14  (D) WARN characterisation: macro_lo/hi is a size-weighted bootstrap
             labelled against an UNWEIGHTED macro_mean; numbers           WARN
  EC-15  census: every named path executed (or UNTESTED with the reason) PASS/WARN

GREEN-ABILITY (a test that cannot go green on a correct implementation is a
PIN, not a test)
--------------------------------------------------------------------------------
  All EC-01..EC-13 checks are equality/identity checks against definitions or
  frozen artifacts that the verified correct tree matches to <=1.1e-16
  (prototype measurements on the REAL tree).  No check asserts current bugs;
  WARN lines (EC-06b/EC-14) are characterisations.  Failure modes hunted:
  NaN-hostile comparisons (non-finite -> inf, never silently skipped), no
  exception swallowing (every exec is guarded and turns FAIL), no constants
  injected from the file under test into any reference.
  FREEZE DISCIPLINE: the equality gates compare the EXECUTED code path against
  the frozen artifact (or the frozen artifact against a definition).  If the
  code is changed deliberately, the corresponding artifact must be regenerated
  together; a code change with a stale artifact then reports FAIL, which is the
  intended staleness signal, not a false red.  Anchor lines are content-pinned:
  a moved/renamed block fails loudly with "-> re-pin" instead of silently
  executing whatever now sits at fixed line numbers.
"""
import os
import re
import sys
import time
import traceback
import warnings

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score, average_precision_score
from scipy.stats import spearmanr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import evs_load as L

OUT = L.OUT
PROD = L.PROD_DIR
HERE = os.path.dirname(os.path.abspath(__file__))
FIGD = os.path.join(L.PROJ, "figures", "data")
P6D = L.P6D
P6E = os.path.join(PROD, "p6e_01_lodo_fold_qc.py")
P11 = os.path.join(PROD, "p11_figure_datapack.py")

SEED = 20261005
TOL = 1e-9


# ===================================================================== helpers
def _lines(path):
    with open(path, encoding="utf-8") as f:
        return f.read().split("\n")


# 生产脚本的可移植路径常量行（2026-10-19 起）。锚点必须与源码**逐字**一致，
# 故由 re.escape 生成，避免手工转义漂移。
_BASE_LINE = re.escape(
    "BASE = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))); "
    "OUT = os.path.dirname(os.path.abspath(__file__))")

def _anchor(path, sp, ep, include_end=True):
    """Content anchor: first line whose STRIPPED text fully matches sp, through
    the (unique) line after it fully matching ep.  Ambiguity is an error, so a
    file that moved or renamed a block fails loudly instead of silently
    executing unrelated lines."""
    src = _lines(path)
    si = [i for i, ln in enumerate(src) if re.fullmatch(sp, ln.strip())]
    if len(si) != 1:
        raise RuntimeError("anchor START hit %dx in %s: %r -> re-pin"
                           % (len(si), os.path.basename(path), sp))
    ei = [i for i, ln in enumerate(src) if i >= si[0] and re.fullmatch(ep, ln.strip())]
    if len(ei) != 1:
        raise RuntimeError("anchor END hit %dx in %s: %r -> re-pin"
                           % (len(ei), os.path.basename(path), ep))
    a, b = si[0] + 1, ei[0] + 1
    if not include_end:
        b -= 1
    return a, b


def _exec_text(ns, text, tag):
    exec(compile(text, "<verbatim:%s>" % tag, "exec"), ns)


def _exec_anchored(ns, blocks, path, key, sp, ep, include_end=True, required=()):
    a, b = _anchor(path, sp, ep, include_end)
    label = "%s:%d-%d" % (os.path.basename(path), a, b)
    _exec_text(ns, "\n".join(_lines(path)[a - 1:b]), label)
    missing = [n for n in required if n not in ns]
    if missing:
        raise RuntimeError("block %s did not bind %s" % (label, missing))
    blocks.append((key, label))
    return label


def _load_def_inj(path, name, **extra):
    fn = L.load_def(path, name)
    fn.__globals__.update(extra)
    return fn


# ------------------------------------------------------------- comparisons
def _cmp_series(a, b, nan_ok=True):
    """(worst diff, note).  Float: |a-b|; a non-finite pattern that differs OR
    non-finite where the other side is finite is +inf; NaN vs NaN counts 0
    only when nan_ok.  Non-float: mismatch count."""
    a = np.asarray(a)
    b = np.asarray(b)
    if a.shape != b.shape:
        return float("inf"), "shape %s vs %s" % (a.shape, b.shape)
    if a.dtype.kind in "fc" or b.dtype.kind in "fc":
        aa = a.astype(float)
        bb = b.astype(float)
        na = ~np.isfinite(aa)
        nb = ~np.isfinite(bb)
        if (na != nb).any():
            return float("inf"), "non-finite pattern differs in %d/%d rows" % (
                int((na != nb).sum()), len(a))
        both = na & nb
        if both.any() and not nan_ok:
            return float("inf"), "non-finite in %d rows" % int(both.sum())
        d = np.abs(np.where(both, 0.0, aa) - np.where(both, 0.0, bb))
        m = float(d.max()) if len(d) else 0.0
        return m, ("" if np.isfinite(m) else "non-finite diff")
    ne = int((a.astype(str) != b.astype(str)).sum())
    return float(ne), ("%d string mismatches" % ne if ne else "")


def _cmp_frames(new, fro, nan_ok=True):
    if list(new.columns) != list(fro.columns):
        return float("inf"), "columns differ: new=%s frozen=%s" % (
            list(new.columns), list(fro.columns))
    if len(new) != len(fro):
        return float("inf"), "rows %d vs frozen %d" % (len(new), len(fro))
    worst = 0.0
    notes = []
    for c in new.columns:
        w, n = _cmp_series(new[c].values, fro[c].values, nan_ok=nan_ok)
        if n:
            notes.append("%s: %s" % (c, n))
        worst = max(worst, w)
    return worst, "; ".join(notes)


def _scalar_diff(x, y):
    x = float(x)
    y = float(y)
    if not (np.isfinite(x) and np.isfinite(y)):
        return float("inf")
    return abs(x - y)


def _read(path):
    return pd.read_csv(path)


# ======================================================= independent references
def _midrank(x):
    """Average ranks with tie handling (midranks); own implementation."""
    x = np.asarray(x, float)
    u, inv, cnt = np.unique(x, return_inverse=True, return_counts=True)
    starts = np.r_[0, np.cumsum(cnt)[:-1]]
    avg = starts + (cnt + 1) / 2.0
    return avg[inv]


def ref_auc(y, s):
    """Mann-Whitney AUROC with midranks (definition (b)); own implementation."""
    y = np.asarray(y, float)
    s = np.asarray(s, float)
    m = np.isfinite(s)
    y, s = y[m], s[m]
    if len(np.unique(y)) < 2:
        return np.nan
    r = _midrank(s)
    n1 = float(y.sum())
    n0 = float(len(y) - n1)
    if n1 == 0 or n0 == 0:
        return np.nan
    return float((r[y == 1].sum() - n1 * (n1 + 1) / 2.0) / (n1 * n0))


def ref_ap(y, s):
    """Step-wise average precision (definition (b)); own implementation."""
    y = np.asarray(y, float)
    s = np.asarray(s, float)
    m = np.isfinite(s)
    y, s = y[m], s[m]
    if len(np.unique(y)) < 2:
        return np.nan
    o = np.argsort(-s, kind="mergesort")
    ys, ss = y[o], s[o]
    npos = ys.sum()
    ap = 0.0
    tp = 0.0
    fp = 0.0
    prev_r = 0.0
    i = 0
    while i < len(ss):
        j = i
        while j + 1 < len(ss) and ss[j + 1] == ss[i]:
            j += 1
        block = ys[i:j + 1]
        tp += block.sum()
        fp += (j - i + 1) - block.sum()
        rec = tp / npos
        prec = tp / (tp + fp)
        ap += (rec - prev_r) * prec
        prev_r = rec
        i = j + 1
    return float(ap)


def ref_spearman(a, b):
    """Spearman = Pearson correlation of midranks; own implementation."""
    ra = _midrank(np.asarray(a, float))
    rb = _midrank(np.asarray(b, float))
    ra = ra - ra.mean()
    rb = rb - rb.mean()
    den = np.sqrt((ra ** 2).sum() * (rb ** 2).sum())
    return float((ra * rb).sum() / den) if den > 0 else np.nan


def _ref_selftest():
    """Controls for the reference implementations AND the comparison helpers:
    if either were broken (constant reference / NaN-blind comparator), every
    check in this module would be vacuous."""
    rng = np.random.default_rng(7)
    r = rng.integers(0, 3, size=400).astype(float)          # heavy ties
    y = (rng.random(400) < 0.4).astype(int)
    a = rng.integers(0, 4, 300).astype(float)
    b = 0.5 * a + rng.integers(0, 3, 300).astype(float)     # ties on both sides
    checks = {
        "auc-vs-sklearn": abs(ref_auc(y, r) - roc_auc_score(y, r)) < 1e-12,
        "ap-vs-sklearn": abs(ref_ap(y, r) - average_precision_score(y, r)) < 1e-12,
        "spearman-vs-scipy": abs(ref_spearman(a, b) - float(spearmanr(a, b).statistic)) < 1e-12,
        "auc-reversal": abs(ref_auc(y, -r) - (1.0 - ref_auc(y, r))) < 1e-12,
        "all-tied=0.5": abs(ref_auc(np.array([0, 0, 1, 1, 1.0]), np.ones(5)) - 0.5) < 1e-12,
        "cmp-nan-vs-number": _cmp_series([1.0, np.nan], [1.0, 2.0])[0] == float("inf"),
        "cmp-nan-vs-nan": _cmp_series([1.0, np.nan], [1.0, np.nan])[0] == 0.0,
        "cmp-detects-diff": _cmp_series([1.0, 2.0], [1.0, 2.001])[0] > 1e-6,
        "cmp-shape": _cmp_series([1.0], [1.0, 2.0])[0] == float("inf"),
        "cmp-string": _cmp_series(np.array(["a", "b"]), np.array(["a", "c"]))[0] == 1.0,
    }
    bad = [k for k, v in checks.items() if not v]
    return (not bad), "references+comparators selftest: %s" % (
        "all 10 OK" if not bad else "BROKEN: " + ", ".join(bad))


# r_failure direction spec (which column, which sign) — the protocol's direction
# table, NOT read from any production file: L1 shift = +; L2 = |.|; SUDO = -AURCC;
# transported = -AUROC.
_METHOD_SPECS = [
    ("domain_auc", "L1", lambda x: x["domain_auc"]),
    ("mmd", "L1", lambda x: x["mmd"]),
    ("psi_mean", "L1", lambda x: x["psi_mean"]),
    ("psi_max", "L1", lambda x: x["psi_max"]),
    ("js_mean", "L1", lambda x: x["js_mean"]),
    ("js_max", "L1", lambda x: x["js_max"]),
    ("pred_mean_shift", "L2", lambda x: x["pred_mean_shift"].abs()),
    ("entropy_mean_shift", "L2", lambda x: x["entropy_mean_shift"].abs()),
    ("entropy_js", "L2", lambda x: x["entropy_js"]),
    ("conf_dist_js", "L2", lambda x: x["conf_dist_js"]),
    ("pred_pos_prev_shift", "L2", lambda x: x["pred_pos_prev_shift"].abs()),
    ("SUDO_AURCC", "L2", lambda x: -x["SUDO_AURCC"]),
    ("AUROC_DR", "L3-primary", lambda x: -x["AUROC_DR"]),
    ("AUROC_IOW", "L3-sens", lambda x: -x["AUROC_IOW"]),
    ("AUROC_OM", "L3-sens", lambda x: -x["AUROC_OM"]),
]

_SIG_COLS = ["task_id", "domain_auc", "mmd", "psi_mean", "psi_max", "js_mean", "js_max",
             "pred_mean_shift", "entropy_mean_shift", "entropy_js", "conf_dist_js",
             "pred_pos_prev_shift"]
_L3_COLS = ["task_id", "SUDO_AURCC", "AUROC_IOW", "AUROC_OM", "AUROC_DR", "w_max", "w_ESS"]


def build_E(O=None):
    """Merge the frozen (or re-executed) outcome with the frozen signal tables.
    This is the input-construction step, not a quantity under test."""
    if O is None:
        O = _read(os.path.join(OUT, "p6d_outcome.csv"))
    SIG = _read(os.path.join(OUT, "p6b_label_free_signals_blinded.csv"))
    L3 = _read(os.path.join(OUT, "p6c_L3_blinded.csv"))
    E = O.merge(SIG[_SIG_COLS], on="task_id", how="left").merge(
        L3[_L3_COLS], on="task_id", how="left")
    E["group"] = E["organism"] + "|" + E["drug"]
    return E


_STATE = {}


# ================================================================== (B)/(C)
def _p6d_state():
    """Execute the p6d first and second layers VERBATIM (content-anchored
    boxes).  Blocks exec'd: constants/DRUGMAP, STORE build, helper defs (name
    loaded), MC build, first-layer loop; then METHODS, second-layer defs, the
    per-method loop.  Excluded: the logging p(...) statements (outside every
    boxed span) and the two to_csv writes to the PRODUCTION dir (p6d:205,
    p6d:304); the artifacts are compared in memory instead."""
    if "p6d" in _STATE:
        return _STATE["p6d"]
    st = dict(ok=False, err=None, O=None, R=None, E=None, blocks=[], secs=None)
    _STATE["p6d"] = st
    if not os.path.exists(P6D):
        st["err"] = "production file absent: %s" % P6D
        return st
    try:
        t0 = time.time()
        blocks = st["blocks"]
        ns = dict(os=os, time=time, warnings=warnings, np=np, pd=pd,
                  LogisticRegression=LogisticRegression, StandardScaler=StandardScaler,
                  roc_auc_score=roc_auc_score, average_precision_score=average_precision_score,
                  spearmanr=spearmanr,
                  __file__=os.path.join(OUT, "p6d_01_unseal_evaluation.py"))
        _exec_anchored(ns, blocks, P6D, "p6d-const+DRUGMAP",
                       _BASE_LINE,
                       r'def _canon\(n\):', include_end=False)
        for nm in ("_canon", "norm_spec", "to_day", "age_lo"):
            fn = _load_def_inj(P6D, nm, pd=pd, DRUGMAP=ns["DRUGMAP"], FOCAL=ns["FOCAL"])
            ns[nm] = fn
            blocks.append(("p6d-def:" + nm, fn.__evs_origin__))
        _exec_anchored(ns, blocks, P6D, "p6d-STORE", r"STORE = \{\}", r"del d, wd, dem")
        fn_bh = _load_def_inj(P6D, "build_hist", pd=pd)
        fn_de = _load_def_inj(P6D, "design", pd=pd)
        fn_pr = _load_def_inj(P6D, "prep", pd=pd, build_hist=fn_bh)   # cross-ref: prep -> build_hist
        for nm, fn in (("build_hist", fn_bh), ("design", fn_de), ("prep", fn_pr)):
            ns[nm] = fn
            blocks.append(("p6d-def:" + nm, fn.__evs_origin__))
        _exec_anchored(ns, blocks, P6D, "p6d-MC",
                       r"SIG = pd\.read_csv\(os\.path\.join\(OUT, \"p6b_label_free_signals_blinded\.csv\"\)\)",
                       r"MC\[\(sym, drug, src\)\] = dict\(sc=sc, mdl=mdl, yte=y\[i_te\], Dte=Dall\[i_te\], af=af\)")
        _exec_anchored(ns, blocks, P6D, "p6d-first-layer", r"rows = \[\]", r"O = pd\.DataFrame\(rows\)")
        _exec_anchored(ns, blocks, P6D, "p6d-METHODS", r"METHODS = \[", r"\]")
        for nm in ("grouped_cv_oof", "safe_auc", "safe_ap", "sens_at_spec"):
            fn = _load_def_inj(P6D, nm, pd=pd, LogisticRegression=LogisticRegression,
                               roc_auc_score=roc_auc_score,
                               average_precision_score=average_precision_score)
            ns[nm] = fn
            blocks.append(("p6d-def:" + nm, fn.__evs_origin__))
        Oex = ns["O"]
        E = build_E(Oex)
        ns["E"] = E
        _exec_anchored(ns, blocks, P6D, "p6d-method-loop", r"res = \[\]", r"R = pd\.DataFrame\(res\)")
        st["O"] = ns["O"]
        st["R"] = ns["R"]
        st["E"] = E
        st["secs"] = time.time() - t0
        st["ok"] = True
    except Exception as e:
        st["err"] = "%s: %s" % (type(e).__name__, e)
    return st


def test_EC01():
    st = _p6d_state()
    if st["err"]:
        return False, "p6d first-layer execution failed: " + st["err"]
    fro = _read(os.path.join(OUT, "p6d_outcome.csv"))
    cols = ["task_id", "organism", "drug", "source_site", "target_site", "n_target",
            "auroc_source_internal", "auroc_target", "d_auroc", "harmful", "stable"]
    new = st["O"][cols].reset_index(drop=True)
    fr = fro[cols].reset_index(drop=True)
    w, note = _cmp_frames(new, fr, nan_ok=False)
    ok = (w <= TOL) and list(new.task_id) == list(fr.task_id)
    obs = ("EXECUTED (%.0fs) first layer, blocks %s; vs frozen p6d_outcome.csv: rows=%d, "
           "task order identical=%s, max_diff=%.3e%s; harmful=%d stable=%d"
           % (st["secs"] or -1, "; ".join(l for _k, l in st["blocks"]),
              len(new), list(new.task_id) == list(fr.task_id), w,
              ("; " + note) if note else "", int(fr.harmful.sum()), int(fr.stable.sum())))
    return ok, obs


def test_EC01D():
    """★ EC-01 的**定义维度**一半（2026-10-05 新增，post-freeze amendment）。

    与 EC-01 的区别：EC-01 把执行结果与**冻结产物**逐项比对（产物验收，重算前必然不符）；
    本判据只断言**执行出来的输出**满足定义/恒等式 —— 与取值无关，故修正前后都应通过。

    依据：GPT6.1「同一路径可能同时需要代码验收和产物验收……不能把整条路径后移，
    连代码验证一起跳过」；EC-03 已按此拆分，EC-01 同理。
    """
    st = _p6d_state()
    if st["err"]:
        return False, "p6d first-layer execution failed: " + st["err"]
    O = st["O"].reset_index(drop=True)
    ident = float(np.abs(O.d_auroc.values
                         - (O.auroc_target.values - O.auroc_source_internal.values)).max())
    haz = (O.d_auroc.values <= -0.05).astype(int)
    sta = (np.abs(O.d_auroc.values) <= 0.02).astype(int)
    flag_bad = int((haz != O.harmful.values).sum()) + int((sta != O.stable.values).sum())
    grp = O.groupby(["organism", "drug", "source_site"]).auroc_source_internal.nunique()
    nonconst = int((grp > 1).sum())
    tid_rule = (O.organism + "_" + O.drug + "_" + O.source_site + "2" + O.target_site).values
    tid_ok = int((O.task_id.values == tid_rule).sum())
    in_range = int(((O.auroc_target.values >= 0) & (O.auroc_target.values <= 1)).sum())
    ok = (ident <= TOL and flag_bad == 0 and nonconst == 0
          and tid_ok == len(O) and in_range == len(O) and len(O) > 0)
    return ok, ("EXECUTED (%.0fs) first layer, blocks %s; **定义核对（不依赖冻结取值）**: "
                "d_auroc identity max_diff=%.2e; harmful/stable flag mismatches=%d; "
                "source-AUROC non-constant over (organism,drug,source)=%d; task_id rule %d/%d; "
                "auroc_target in [0,1] %d/%d; rows=%d (counts are a characterisation, NOT pinned)"
                % (st["secs"] or -1, "; ".join(l for _k, l in st["blocks"]),
                   ident, flag_bad, nonconst, tid_ok, len(O), in_range, len(O), len(O)))


def test_EC02():
    O = _read(os.path.join(OUT, "p6d_outcome.csv"))
    ident = float(np.abs(O.d_auroc.values - (O.auroc_target.values - O.auroc_source_internal.values)).max())
    haz = (O.d_auroc.values <= -0.05).astype(int)
    sta = (np.abs(O.d_auroc.values) <= 0.02).astype(int)
    flag_bad = int((haz != O.harmful.values).sum()) + int((sta != O.stable.values).sum())
    grp = O.groupby(["organism", "drug", "source_site"]).auroc_source_internal.nunique()
    nonconst = int((grp > 1).sum())
    tid_rule = (O.organism + "_" + O.drug + "_" + O.source_site + "2" + O.target_site).values
    tid_ok = int((O.task_id.values == tid_rule).sum())
    q = O.d_auroc.values
    p10, p50, p90 = (float(np.quantile(q, x)) for x in (0.10, 0.50, 0.90))
    # The identities/flags/keys are DEFINITIONAL and must hold for any regenerated
    # artifact.  The counts (202/65/49) are a characterisation of the CURRENT
    # frozen artifact, reported but NOT asserted: the CS-01/02/03 corrections will
    # legitimately change them when p6d is re-run (a pin here would false-RED the
    # corrected tree).
    ok = (ident <= TOL and flag_bad == 0 and nonconst == 0 and tid_ok == len(O))
    obs = ("frozen reconciliation vs DEFINITIONS: identity max_diff=%.2e; flag mismatches=%d; "
           "source-AUROC non-constant over (organism,drug,source)=%d; task_id rule %d/%d; "
           "current-artifact counts (characterisation, not pinned): rows=%d harmful=%d(%.1f%%) "
           "stable=%d; P10=%+.4f P50=%+.4f P90=%+.4f P90-P10=%.4f"
           % (ident, flag_bad, nonconst, tid_ok, len(O), len(O), int(O.harmful.sum()),
              100.0 * O.harmful.mean(), int(O.stable.sum()), p10, p50, p90, p90 - p10))
    return ok, obs


def test_EC03():
    st = _p6d_state()
    if st["err"]:
        return False, "p6d second-layer execution failed: " + st["err"]
    fro = _read(os.path.join(OUT, "p6d_signal_eval.csv"))
    new = st["R"]
    w, note = _cmp_frames(new, fro, nan_ok=False)
    ok = (w <= TOL) and list(new.method) == list(fro.method)
    obs = ("EXECUTED per-method loop (285-317, 1000-draw cluster bootstrap, LODO, in-file LOO) vs "
           "frozen p6d_signal_eval.csv: 15 methods, columns=%d, order identical=%s, max_diff=%.3e%s; "
           "exemplars domain_auc harmful_auc=%.6f [%.6f,%.6f], mmd=%.6f"
           % (len(new.columns), list(new.method) == list(fro.method), w,
              ("; " + note) if note else "",
              float(fro.loc[fro.method == "domain_auc", "harmful_auc"].iat[0]),
              float(fro.loc[fro.method == "domain_auc", "auc_lo"].iat[0]),
              float(fro.loc[fro.method == "domain_auc", "auc_hi"].iat[0]),
              float(fro.loc[fro.method == "mmd", "harmful_auc"].iat[0])))
    return ok, obs


def test_EC04():
    sel_ok, sel_obs = _ref_selftest()
    E = build_E()
    fro = _read(os.path.join(OUT, "p6d_signal_eval.csv")).set_index("method")
    wa = 0.0
    wp = 0.0
    for name, _role, fn in _METHOD_SPECS:
        r = fn(E).values.astype(float)
        y = E["harmful"].values
        okm = np.isfinite(r)
        wa = max(wa, _scalar_diff(ref_auc(y[okm], r[okm]), fro.loc[name, "harmful_auc"]))
        wp = max(wp, _scalar_diff(ref_ap(y[okm], r[okm]), fro.loc[name, "harmful_auprc"]))
    ok = sel_ok and wa <= TOL and wp <= TOL
    obs = ("15 methods; max abs(own-midrank AUROC - frozen harmful_auc)=%.3e; "
           "max abs(own step-AP - frozen harmful_auprc)=%.3e; %s" % (wa, wp, sel_obs))
    return ok, obs


def test_EC05():
    E = build_E()
    fro = _read(os.path.join(OUT, "p6d_signal_eval.csv")).set_index("method")
    g = E["group"].values
    ugs = np.unique(g)
    dd = E["d_auroc"].values
    worst = 0.0
    for name, _role, fn in _METHOD_SPECS:
        r = fn(E).values.astype(float)
        y = E["harmful"].values
        rng = np.random.default_rng(SEED + 7)     # per-method stream (artifact provenance)
        sps = []
        aucs = []
        for _ in range(1000):
            pick = rng.choice(len(ugs), len(ugs), replace=True)
            idx = np.concatenate([np.where(g == ugs[i])[0] for i in pick])
            if len(np.unique(y[idx])) < 2:
                continue
            sps.append(ref_spearman(r[idx], dd[idx]))
            aucs.append(ref_auc(y[idx], r[idx]))
        for got, key in ((np.percentile(sps, 2.5), "sp_lo"), (np.percentile(sps, 97.5), "sp_hi"),
                         (np.percentile(aucs, 2.5), "auc_lo"), (np.percentile(aucs, 97.5), "auc_hi")):
            worst = max(worst, _scalar_diff(got, fro.loc[name, key]))
    mmd_fn = [f for n, _r, f in _METHOD_SPECS if n == "mmd"][0]
    r = mmd_fn(E).values.astype(float)
    y = E["harmful"].values
    rngc = np.random.default_rng(SEED + 7)
    ca = []
    for _ in range(1000):
        pick = rngc.choice(len(ugs), len(ugs), replace=True)
        idx = np.concatenate([np.where(g == ugs[i])[0] for i in np.unique(pick)])
        if len(np.unique(y[idx])) < 2:
            continue
        ca.append(ref_auc(y[idx], r[idx]))
    lo_c, hi_c = (float(x) for x in np.percentile(ca, [2.5, 97.5]))
    shift = max(abs(lo_c - fro.loc["mmd", "auc_lo"]), abs(hi_c - fro.loc["mmd", "auc_hi"]))
    ok = worst <= TOL and shift > 1e-3
    obs = ("own replay (clusters = organism-drug, WITH multiplicity, fresh rng(SEED+7) per method) vs "
           "frozen 15x4 CI columns: max_diff=%.3e; set-dedup CONTROL on mmd -> CI [%.6f,%.6f] vs frozen "
           "[%.6f,%.6f] (shift %.4f: the check goes red if multiplicity were dropped)"
           % (worst, lo_c, hi_c, float(fro.loc["mmd", "auc_lo"]), float(fro.loc["mmd", "auc_hi"]), shift))
    return ok, obs


# ====================================================================== (A)
def _p11_state():
    """Execute p11_figure_datapack.py's blocks VERBATIM with FIG redirected to a
    scratch dir (so the module is read-only on the production tree) and p()
    replaced by a file sink.  Line 17 (FIG = production path) and lines 21-23
    (open README / def p) are the only lines not exec'd; every other line of
    every named panel block is executed verbatim."""
    if "p11" in _STATE:
        return _STATE["p11"]
    st = dict(ok=False, err=None, scratch=None, blocks=[], secs=None)
    _STATE["p11"] = st
    if not os.path.exists(P11):
        st["err"] = "production file absent: %s" % P11
        return st
    try:
        import tempfile
        t0 = time.time()
        scratch = tempfile.mkdtemp(prefix="evs_t13_p11_")
        st["scratch"] = scratch
        blocks = st["blocks"]
        ns = dict(os=os, np=np, pd=pd, roc_auc_score=roc_auc_score,
                  OUT=OUT, FIG=scratch, SEED=SEED, B=4000,
                  __file__=os.path.join(OUT, "p11_figure_datapack.py"))
        ns["rng"] = np.random.default_rng(SEED)
        logf = open(os.path.join(scratch, "p11_log.txt"), "w", encoding="utf-8")

        def _sink(*a):
            logf.write(" ".join(str(x) for x in a) + "\n")

        ns["p"] = _sink
        base = _BASE_LINE
        _exec_anchored(ns, blocks, P11, "p11-BASE(read-only paths)", base, base)
        mk = r"os\.makedirs\(FIG, exist_ok=True\)"
        _exec_anchored(ns, blocks, P11, "p11-mkdir(scratch FIG)", mk, mk)
        sd = r"SEED = 20261005; rng = np\.random\.default_rng\(SEED\); B = 4000"
        _exec_anchored(ns, blocks, P11, "p11-SEED", sd, sd)
        _exec_anchored(ns, blocks, P11, "p11-inputs",
                       r'O = pd\.read_csv\(os\.path\.join\(OUT, "p6d_outcome\.csv"\)\)',
                       r'E\["direction"\] = E\.source_site \+ "->" \+ E\.target_site')
        _exec_anchored(ns, blocks, P11, "p11-fig1b",
                       r"cand = len\(UNI\); elig = int\(UNI\.eligible\.sum\(\)\)",
                       r'\]\)\.to_csv\(os\.path\.join\(FIG, "fig1b_universe_long\.csv"\), index=False\)')
        _exec_anchored(ns, blocks, P11, "p11-fig1c",
                       r'f1c = E\[\["task_id", "organism", "drug", "source_site", "target_site", "direction",',
                       r'\)\.to_csv\(os\.path\.join\(FIG, "fig1c_summary\.csv"\), index=False\)')
        _exec_anchored(ns, blocks, P11, "p11-fig1d",
                       r'f1d = E\.groupby\("direction"\)\["d_auroc"\]\.describe\(\)\[\["count", "mean", "std", "min", "25%", "50%", "75%", "max"\]\]',
                       r'f1d\.to_csv\(os\.path\.join\(FIG, "fig1d_by_direction\.csv"\)\)')
        _exec_anchored(ns, blocks, P11, "p11-fig1e",
                       r'CFG12 = pd\.read_csv\(os\.path\.join\(OUT, "p12_configs_corrected\.csv"\)\)',
                       r'os\.path\.join\(FIG, "fig1e_configs_aligned\.csv"\), index=False\)')
        _exec_anchored(ns, blocks, P11, "p11-fig3a",
                       r'f3a = SET\[\["method", "role", "n_eval_fold", "C_macro_eval", "C_lo", "C_hi"\]\]\.rename\(',
                       r'f3a\.to_csv\(os\.path\.join\(FIG, "fig3a_lodo\.csv"\), index=False\)')
        _exec_anchored(ns, blocks, P11, "p11-fig3c",
                       r'f3c = SET\[SET\.method == "SUDO_AURCC"\]\[\["A_pooled_all202", "B_pooled_eval", "C_macro_eval",',
                       r'f3c\.to_csv\(os\.path\.join\(FIG, "fig3c_sudo_decomposition\.csv"\), index=False\)')
        _exec_anchored(ns, blocks, P11, "p11-METHODS",
                       r'METHODS = \[\("domain_auc", "L1", lambda x: x\["domain_auc"\]\), \("mmd", "L1", lambda x: x\["mmd"\]\),',
                       r'\("AUROC_OM", "L3-sens", lambda x: -x\["AUROC_OM"\]\)\]')
        _exec_anchored(ns, blocks, P11, "p11-fig4",
                       r'f4 = E\[\["task_id", "organism", "drug", "source_site", "target_site",',
                       r'f4\.to_csv\(os\.path\.join\(FIG, "fig4_l3_scatter\.csv"\), index=False\)')
        _exec_anchored(ns, blocks, P11, "p11-fig4d",
                       r"rows4 = \[\]",
                       r'f4d\.to_csv\(os\.path\.join\(FIG, "fig4d_mae\.csv"\), index=False\)')
        _exec_anchored(ns, blocks, P11, "p11-fig4e+rho_ess",
                       r'f4e = f4\[f4\.AUROC_DR > 1\]\[\["task_id", "organism", "drug", "source_site", "target_site",',
                       r'rho_ess = pd\.Series\(f4\.w_ESS\)\.corr\(pd\.Series\(np\.abs\(f4\.AUROC_DR - f4\.auroc_target\)\), method="spearman"\)')
        _exec_anchored(ns, blocks, P11, "p11-fig2",
                       r'f2 = SEV\[\["method", "role", "coverage", "harmful_auc", "auc_lo", "auc_hi",',
                       r'f2\.to_csv\(os\.path\.join\(FIG, "fig2_forest\.csv"\), index=False\)')
        _exec_anchored(ns, blocks, P11, "p11-fig5",
                       r'O2 = O\[\["task_id", "harmful"\]\]',
                       r'f5\.to_csv\(os\.path\.join\(FIG, "fig5_subsample\.csv"\), index=False\)')
        _exec_anchored(ns, blocks, P11, "p11-fig5b",
                       r'_lv = f5\[f5\.level\.isin\(\["100", "250", "500", "1000", "full"\]\)\]',
                       r'pd\.DataFrame\(rows5b\)\.to_csv\(os\.path\.join\(FIG, "fig5b_across_n_range\.csv"\), index=False\)')
        _exec_anchored(ns, blocks, P11, "p11-figS6",
                       r'FC\.to_csv\(os\.path\.join\(FIG, "figS6_fold_composition\.csv"\), index=False\)',
                       r'S6S\.to_csv\(os\.path\.join\(FIG, "figS6_fold_summary\.csv"\), index=False\)')
        _exec_anchored(ns, blocks, P11, "p11-tableS1",
                       r'tS1 = E\[\["task_id", "organism", "drug", "source_site", "target_site", "n_target",',
                       r'tS1\.to_csv\(os\.path\.join\(FIG, "tableS1_tasks\.csv"\), index=False\)')
        _exec_anchored(ns, blocks, P11, "p11-tableS2",
                       r'tS2 = SEV\.merge\(LAGG\[\["method", "pooled_OOF", "macro_mean", "n_evaluable_folds",',
                       r'tS2\.to_csv\(os\.path\.join\(FIG, "tableS2_metrics\.csv"\), index=False\)')
        logf.close()
        st["ns"] = ns
        st["secs"] = time.time() - t0
        st["ok"] = True
    except Exception as e:
        st["err"] = "%s: %s" % (type(e).__name__, e)
    return st


def _p11_cmp(fname, nan_ok=True):
    st = _p11_state()
    if st["err"]:
        raise RuntimeError("p11 execution failed: " + st["err"])
    new = _read(os.path.join(st["scratch"], fname))
    fro = _read(os.path.join(FIGD, fname))
    return _cmp_frames(new, fro, nan_ok=nan_ok)


def _cat_rule(lo, hi):
    return np.where(np.asarray(lo) > 0.5, "above_chance",
                    np.where(np.asarray(hi) < 0.5, "below_chance", "indistinguishable"))


def test_EC06():
    st = _p11_state()
    if st["err"]:
        return False, "p11 execution failed: " + st["err"]
    w, note = _p11_cmp("fig3a_lodo.csv")
    fro = _read(os.path.join(FIGD, "fig3a_lodo.csv"))
    new = _read(os.path.join(st["scratch"], "fig3a_lodo.csv"))
    rule_fro = _cat_rule(fro.lo_95.values, fro.hi_95.values)
    rule_new = _cat_rule(new.lo_95.values, new.hi_95.values)
    bad = int((rule_fro != fro.category.values).sum()) + int((rule_new != new.category.values).sum())
    FC = _read(os.path.join(OUT, "p6e_lodo_fold_composition.csv"))
    n_fold_def = int(FC.auroc_definable.sum())
    n_task_eval = int(FC.loc[FC.auroc_definable == 1, "n"].sum())
    folds_ok = bool((new.n_eval_fold == n_fold_def).all()) and \
        bool((new.n_tasks_eval == n_task_eval).all())
    counts = fro.category.value_counts().to_dict()
    ok = (w <= TOL) and bad == 0 and folds_ok
    obs = ("exec block vs frozen: 15 rows, max_diff=%.3e%s; category rule (CI vs 0.5) mismatches=%d/15 "
           "on frozen AND on executed frame; evaluable folds=%d and tasks=%d per definition "
           "(p6e_lodo_fold_composition: definable folds / their n-sum); categories: %s"
           % (w, ("; " + note) if note else "", bad, n_fold_def, n_task_eval, counts))
    return ok, obs


def _warn_EC06b():
    fro = _read(os.path.join(FIGD, "fig3a_lodo.csv")).set_index("method")
    LAGG = _read(os.path.join(OUT, "p6e_lodo_aggregation.csv")).set_index("method")
    lo_e = LAGG.macro_lo.reindex(fro.index).values
    hi_e = LAGG.macro_hi.reindex(fro.index).values
    cat_f = _cat_rule(fro.lo_95.values, fro.hi_95.values)
    cat_e = _cat_rule(lo_e, hi_e)
    flips = []
    for i, m in enumerate(fro.index):
        if cat_f[i] != cat_e[i]:
            flips.append("%s %s->%s (p6f CI [%.3f,%.3f] vs p6e weighted CI [%.3f,%.3f])"
                         % (m, cat_f[i], cat_e[i], fro.lo_95.values[i], fro.hi_95.values[i],
                            lo_e[i], hi_e[i]))
    c_f = {k: int((cat_f == k).sum()) for k in ("above_chance", "below_chance", "indistinguishable")}
    c_e = {k: int((cat_e == k).sum()) for k in ("above_chance", "below_chance", "indistinguishable")}
    return ("executed frame (p6f C_lo/C_hi, unweighted) vs the p6e size-weighted CI convention",
            "%d/15 categories flip: %s; counts p6f/unweighted %s vs p6e/weighted %s"
            % (len(flips), "; ".join(flips) if flips else "none", c_f, c_e))


def test_EC07():
    st = _p11_state()
    if st["err"]:
        return False, "p11 execution failed: " + st["err"]
    w, note = _p11_cmp("fig3c_sudo_decomposition.csv")
    fro = _read(os.path.join(FIGD, "fig3c_sudo_decomposition.csv")).iloc[0]
    SET = _read(os.path.join(OUT, "p6f_lodo_settlement.csv"))
    row = SET[SET.method == "SUDO_AURCC"].iloc[0]
    LAGG = _read(os.path.join(OUT, "p6e_lodo_aggregation.csv"))
    p6d = _read(os.path.join(OUT, "p6d_signal_eval.csv"))
    # decomposition identity (stage difference, earlier - later; p6f:146 AtB=A-Bv, BtC=Bv-C):
    id_sub = abs(fro.task_subset_effect - (fro.pooled_202task - fro.pooled_138task))
    id_agg = abs(fro.aggregation_effect - (fro.pooled_138task - fro.macro_138task))
    # the same identity must hold on all 15 rows of the p6f source artifact
    id_all = max(float(np.abs(SET.AtB - (SET.A_pooled_all202 - SET.B_pooled_eval)).max()),
                 float(np.abs(SET.BtC - (SET.B_pooled_eval - SET.C_macro_eval)).max()))
    x_set = max(abs(fro.pooled_202task - row.A_pooled_all202),
                abs(fro.pooled_138task - row.B_pooled_eval))
    x_p6e = _scalar_diff(fro.macro_138task,
                         LAGG.loc[LAGG.method == "SUDO_AURCC", "macro_mean"].iat[0])
    x_p6d = _scalar_diff(fro.pooled_202task,
                         p6d.loc[p6d.method == "SUDO_AURCC", "grouped_cv_auc"].iat[0])
    ok = (w <= TOL) and max(id_sub, id_agg, id_all, x_set, x_p6e, x_p6d) <= TOL
    obs = ("exec == frozen (max_diff=%.3e%s); identity subset=A-B (%.2e), aggregation=B-C (%.2e), "
           "same over all 15 p6f rows (%.2e); cross-artifact SUDO_AURCC: pooled_202 == p6f A (%.2e), "
           "pooled_138 == p6f B (%.2e), macro_138 == p6e macro_mean (%.2e), "
           "pooled_202 == p6d grouped_cv_auc (%.2e)"
           % (w, ("; " + note) if note else "", id_sub, id_agg, id_all,
              abs(fro.pooled_202task - row.A_pooled_all202),
              abs(fro.pooled_138task - row.B_pooled_eval), x_p6e, x_p6d))
    return ok, obs


def test_EC08():
    st = _p11_state()
    if st["err"]:
        return False, "p11 execution failed: " + st["err"]
    w, note = _p11_cmp("fig4d_mae.csv")
    wc = 0.0
    notes = []
    for lab in ("DR", "IOW", "OM"):
        ww, nn = _p11_cmp("fig4d_calibration_%s.csv" % lab)
        wc = max(wc, ww)
        if nn:
            notes.append("%s: %s" % (lab, nn))
    we, _n = _p11_cmp("fig4e_dr_out_of_range.csv")
    fro = _read(os.path.join(FIGD, "fig4d_mae.csv")).set_index("estimator")
    scatter = _read(os.path.join(FIGD, "fig4_l3_scatter.csv"))
    worst_mae = 0.0
    worst_rho = 0.0
    cal_bad = 0.0
    EDGES = [-np.inf, -0.10, -0.05, 0.0, 0.05, np.inf]
    LABELS = ["<-0.10", "-0.10~-0.05", "-0.05~0", "0~0.05", ">0.05"]
    for lab, col in (("DR", "AUROC_DR"), ("IOW", "AUROC_IOW"), ("OM", "AUROC_OM")):
        v = scatter[col].values.astype(float)
        t = scatter["auroc_target"].values.astype(float)
        m = np.isfinite(v)
        worst_mae = max(worst_mae, abs(float(np.mean(np.abs(v[m] - t[m]))) - fro.loc[lab, "mae"]))
        worst_rho = max(worst_rho, _scalar_diff(
            float(spearmanr(v[m], t[m]).statistic), fro.loc[lab, "rho_pred_true"]))
        pred_d = v - scatter["auroc_source_internal"].values
        b = pd.cut(pred_d[m], bins=EDGES, labels=LABELS, right=False)
        g = pd.DataFrame({"b": b, "pred": pred_d[m], "act": scatter["d_auroc"].values[m]}) \
            .groupby("b", observed=False).agg(n=("pred", "size"), pm=("pred", "mean"), am=("act", "mean"))
        cf = _read(os.path.join(FIGD, "fig4d_calibration_%s.csv" % lab))
        d1 = float(np.abs(g["n"].values - cf["n"].values).max())
        d2 = float(np.abs(g["pm"].values - cf["pred_d_mean"].values).max())
        d3 = float(np.abs(g["am"].values - cf["actual_d_mean"].values).max())
        cal_bad = max(cal_bad, d1, d2, d3)
    ns = st["ns"]
    rho_ess_exec = float(ns["rho_ess"])
    w_ess = scatter["w_ESS"].values.astype(float)
    resid = np.abs(scatter["AUROC_DR"].values.astype(float) - scatter["auroc_target"].values.astype(float))
    rho_ess_ref = float(spearmanr(w_ess, resid).statistic)
    ok = (w <= TOL) and (wc <= TOL) and (we <= TOL) and worst_mae <= TOL and worst_rho <= TOL \
        and cal_bad <= TOL and _scalar_diff(rho_ess_exec, rho_ess_ref) <= TOL
    obs = ("exec == frozen: MAE %.3e, calibration max_diff=%.3e, fig4e rows %.3e%s; independent refs: "
           "numpy MAE diff %.2e, scipy Spearman diff %.2e, calibration group means diff %.2e; "
           "rho_ess exec=%+.6f vs own scipy ref=%+.6f; frozen MAE DR=%.4f IOW=%.4f OM=%.4f"
           % (w, max(wc, cal_bad), we, ("; " + note) if note else "", worst_mae, worst_rho,
              cal_bad, rho_ess_exec, rho_ess_ref,
              fro.loc["DR", "mae"], fro.loc["IOW", "mae"], fro.loc["OM", "mae"]))
    return ok, obs


def test_EC09():
    st = _p11_state()
    if st["err"]:
        return False, "p11 execution failed: " + st["err"]
    w, note = _p11_cmp("fig1c_summary.csv")
    fro = _read(os.path.join(FIGD, "fig1c_summary.csv")).iloc[0]
    O = _read(os.path.join(OUT, "p6d_outcome.csv"))
    q = O.d_auroc.values
    refs = {
        "p05": np.quantile(q, .05), "p10": np.quantile(q, .10), "p25": np.quantile(q, .25),
        "p50": np.quantile(q, .50), "p75": np.quantile(q, .75), "p90": np.quantile(q, .90),
        "p95": np.quantile(q, .95), "mean": np.mean(q), "sd": np.std(q, ddof=1),
        "p90_minus_p10": np.quantile(q, .90) - np.quantile(q, .10),
        "dmin": np.min(q), "dmax": np.max(q),
    }
    worst = max(_scalar_diff(refs[k], fro[k]) for k in refs)
    # n_tasks/harmful/stable are characterisation of the current artifact, not
    # asserted (they change legitimately once the corrected benchmark is re-frozen).
    counts_note = ("current-artifact counts (characterisation): n/harmful/stable = %d/%d/%d"
                   % (int(fro.n_tasks), int(fro.n_harmful), int(fro.n_stable)))
    ok = (w <= TOL) and worst <= TOL
    obs = ("exec == frozen (max_diff=%.3e%s); independent np.quantile(linear)/mean/sd(ddof=1) vs frozen "
           "columns: max_diff=%.3e; P10=%+.4f P50=%+.4f P90=%+.4f P90-P10=%.4f; %s"
           % (w, ("; " + note) if note else "", worst, refs["p10"], refs["p50"], refs["p90"],
              refs["p90_minus_p10"], counts_note))
    return ok, obs


def test_EC10():
    st = _p11_state()
    if st["err"]:
        return False, "p11 execution failed: " + st["err"]
    w, note = _p11_cmp("fig5_subsample.csv")
    wb, nb = _p11_cmp("fig5b_across_n_range.csv")
    fro = _read(os.path.join(FIGD, "fig5_subsample.csv"))
    RW = _read(os.path.join(OUT, "p10b_raw.csv"))
    O = _read(os.path.join(OUT, "p6d_outcome.csv"))
    RW = RW.merge(O[["task_id", "harmful"]], on="task_id", how="left")
    levels = ["100", "250", "500", "1000", "full"]
    worst = 0.0
    n_paired = 0
    for name, _role, fn in _METHOD_SPECS:
        per_level = {}
        for lv in levels:
            sub = RW[RW.level == lv]
            vals = {}
            for rr, gg in sub.groupby("rep"):
                gg = gg.sort_values("task_id")
                s = fn(gg).values.astype(float)
                yy = gg["harmful"].values.astype(int)
                mm = np.isfinite(s)
                if len(np.unique(yy[mm])) < 2:
                    continue
                vals[rr] = ref_auc(yy[mm], s[mm])
            v = pd.Series(vals)
            per_level[lv] = v
            frow = fro[(fro.method == name) & (fro.level == lv)].iloc[0]
            if len(v) != int(frow.n_repeat):
                worst = float("inf")
            worst = max(worst,
                        _scalar_diff(np.mean(v.values), frow.auroc_mean),
                        _scalar_diff(np.percentile(v.values, 2.5), frow.auroc_lo),
                        _scalar_diff(np.percentile(v.values, 97.5), frow.auroc_hi))
        d = (per_level["full"] - per_level["100"]).dropna()
        n_paired += 1
        frow = fro[(fro.method == name) & (fro.level == "PAIRED_full_minus_100")].iloc[0]
        if len(d) != int(frow.n_repeat):
            worst = float("inf")
        worst = max(worst,
                    _scalar_diff(np.mean(d.values), frow.auroc_mean),
                    _scalar_diff(np.percentile(d.values, 2.5), frow.auroc_lo),
                    _scalar_diff(np.percentile(d.values, 97.5), frow.auroc_hi))
    ok = (w <= TOL) and (wb <= TOL) and worst <= TOL
    obs = ("exec == frozen (fig5 max_diff=%.3e%s; fig5b max_diff=%.3e%s); independent replay of all "
           "15 methods x 5 levels x 200 reps + %d paired rows with own midrank AUC vs frozen "
           "mean/2.5%%/97.5%%: max_diff=%.3e; frozen domain_auc full mean=%.6f"
           % (w, ("; " + note) if note else "", wb, ("; " + nb) if nb else "",
              n_paired, worst,
              float(fro[(fro.method == "domain_auc") & (fro.level == "full")].auroc_mean.iat[0])))
    return ok, obs


def test_EC11():
    st = _p11_state()
    if st["err"]:
        return False, "p11 execution failed: " + st["err"]
    worst = 0.0
    bad = []
    for fname in ["fig1b_universe.csv", "fig1b_universe_long.csv", "fig1c_delta_distribution.csv",
                  "fig1d_by_direction.csv", "fig1e_configs.csv", "fig1e_configs_aligned.csv",
                  "fig1e_configs_full.csv", "fig2_forest.csv", "fig4_l3_scatter.csv",
                  "figS6_fold_composition.csv", "figS6_fold_summary.csv", "tableS1_tasks.csv",
                  "tableS2_metrics.csv"]:
        ww, nn = _p11_cmp(fname)
        worst = max(worst, ww)
        if nn:
            bad.append("%s(%s)" % (fname, nn))
    CFG12 = _read(os.path.join(OUT, "p12_configs_corrected.csv"))
    fro = _read(os.path.join(FIGD, "fig1e_configs.csv"))
    P1set = set(CFG12[(CFG12.config == "P1 History+LR") & CFG12.eligible].task_id)
    worst_fe = 0.0
    for cfg in ["P1 History+LR", "S1 Core+LR", "S2 History+LGBM", "S3 Core+LGBM"]:
        own = CFG12[CFG12.config == cfg]
        s = own[own.task_id.isin(P1set)]
        d = s.d_auroc.values.astype(float)
        frow = fro[fro.config == cfg].iloc[0]
        ref = dict(median=np.median(d), p10=np.quantile(d, .1), p90=np.quantile(d, .9),
                   pm=float(np.quantile(d, .9) - np.quantile(d, .1)),
                   nh=int((d <= -0.05).sum()), pct=100.0 * float((d <= -0.05).mean()),
                   n=len(d), nelig=int(own.eligible.sum()))
        worst_fe = max(worst_fe,
                       _scalar_diff(ref["median"], frow.median_aligned),
                       _scalar_diff(ref["p10"], frow.p10_aligned),
                       _scalar_diff(ref["p90"], frow.p90_aligned),
                       _scalar_diff(ref["pm"], frow.p90_minus_p10_aligned),
                       abs(ref["nh"] - int(frow.n_harmful_aligned)),
                       _scalar_diff(ref["pct"], frow.pct_harmful_aligned),
                       abs(ref["n"] - int(frow.n_aligned)),
                       abs(ref["nelig"] - int(frow.n_eligible_own)))
    t2 = _read(os.path.join(FIGD, "tableS2_metrics.csv"))
    integ = (len(t2) == 15 and t2.method.nunique() == 15 and int(t2.isna().sum().sum()) == 0)
    ok = (worst <= TOL) and (worst_fe <= TOL) and integ
    obs = ("exec == frozen for 13 panel files (max_diff=%.3e%s); independent fig1e recompute over the "
           "P1-202 aligned set (median/P10/P90/harmful counts/pct/eligible) max_diff=%.3e; tableS2 "
           "merge integrity rows=15 unique methods=15 nulls=%d"
           % (worst, ("; " + "; ".join(bad)) if bad else "", worst_fe,
              int(t2.isna().sum().sum())))
    return ok, obs


# ====================================================================== (D)
def _p6e_state():
    """Execute p6e's fold-composition block and aggregation loop VERBATIM.
    Excluded: the two to_csv writes to the production dir (p6e:68, 144-145);
    the frames are rebuilt from the exec'd lists and compared in memory."""
    if "p6e" in _STATE:
        return _STATE["p6e"]
    st = dict(ok=False, err=None, FC=None, R=None, FT=None, blocks=[], secs=None)
    _STATE["p6e"] = st
    if not os.path.exists(P6E):
        st["err"] = "production file absent: %s" % P6E
        return st
    try:
        t0 = time.time()
        blocks = st["blocks"]
        E = build_E()
        y_all = E["harmful"].values.astype(int)
        g_all = E["group"].values
        ns = dict(np=np, pd=pd, os=os, OUT=OUT, E=E, y_all=y_all, g_all=g_all,
                  METHODS=_METHOD_SPECS)
        _exec_anchored(ns, blocks, P6E, "p6e-fold-composition", r"rows = \[\]",
                       r'FC = pd\.DataFrame\(rows\)\.sort_values\(\["organism", "drug"\]\)\.reset_index\(drop=True\)')
        fit_predict = _load_def_inj(P6E, "fit_predict", pd=pd, y_all=y_all, g_all=g_all,
                                    LogisticRegression=LogisticRegression)
        safe_auc = _load_def_inj(P6E, "safe_auc", pd=pd, y_all=y_all, g_all=g_all,
                                 roc_auc_score=roc_auc_score)
        run = _load_def_inj(P6E, "run", pd=pd, y_all=y_all, g_all=g_all,
                            fit_predict=fit_predict, safe_auc=safe_auc)
        boot_macro = _load_def_inj(P6E, "boot_macro", pd=pd)
        rng = np.random.default_rng(SEED)
        boot_macro.__globals__["rng"] = rng
        for nm, fn in (("fit_predict", fit_predict), ("safe_auc", safe_auc),
                       ("run", run), ("boot_macro", boot_macro)):
            blocks.append(("p6e-def:" + nm, fn.__evs_origin__))
        ns.update(run=run, boot_macro=boot_macro, rng=rng)
        _exec_anchored(ns, blocks, P6E, "p6e-aggregation-loop", r"res = \[\]; foldtab = \[\]",
                       r"foldtab\.append\(dict\(method=name, group=k, fold_auc=a, fold_n=nn\)\)")
        st["FC"] = ns["FC"]
        st["R"] = pd.DataFrame(ns["res"])
        st["FT"] = pd.DataFrame(ns["foldtab"])
        st["secs"] = time.time() - t0
        st["ok"] = True
    except Exception as e:
        st["err"] = "%s: %s" % (type(e).__name__, e)
    return st


def test_EC12():
    st = _p6e_state()
    if st["err"]:
        return False, "p6e execution failed: " + st["err"]
    wC, nC = _cmp_frames(st["FC"], _read(os.path.join(OUT, "p6e_lodo_fold_composition.csv")),
                         nan_ok=False)
    wR, nR = _cmp_frames(st["R"], _read(os.path.join(OUT, "p6e_lodo_aggregation.csv")),
                         nan_ok=False)
    wF, nF = _cmp_frames(st["FT"], _read(os.path.join(OUT, "p6e_lodo_fold_auc_long.csv")),
                         nan_ok=False)
    ok = max(wC, wR, wF) <= TOL
    r0 = st["R"]
    obs = ("EXECUTED fold composition + aggregation loop (%.0fs): composition 38 folds max_diff=%.3e; "
           "aggregation 15 methods x 12 cols max_diff=%.3e%s; fold-long 375 rows max_diff=%.3e%s; "
           "macro_mean max_diff=%.3e, macro_lo/hi max_diff=%.3e/%.3e"
           % (st["secs"] or -1, wC, wR, ("; " + nR) if nR else "", wF, ("; " + nF) if nF else "",
              float(np.abs(r0.macro_mean.values -
                           _read(os.path.join(OUT, "p6e_lodo_aggregation.csv")).macro_mean.values).max()),
              float(np.abs(r0.macro_lo.values -
                           _read(os.path.join(OUT, "p6e_lodo_aggregation.csv")).macro_lo.values).max()),
              float(np.abs(r0.macro_hi.values -
                           _read(os.path.join(OUT, "p6e_lodo_aggregation.csv")).macro_hi.values).max())))
    return ok, obs


def test_EC13():
    if not os.path.exists(P6E):
        raise RuntimeError("p6e absent")
    y_all = np.array([1, 1, 1, 1, 0, 1, 0, 0, 0], int)
    g_all = np.array(["A"] * 3 + ["B"] * 3 + ["C"] * 3)
    fit_predict = _load_def_inj(P6E, "fit_predict", pd=pd, y_all=y_all, g_all=g_all,
                                LogisticRegression=LogisticRegression)
    safe_auc = _load_def_inj(P6E, "safe_auc", pd=pd, y_all=y_all, g_all=g_all,
                             roc_auc_score=roc_auc_score)
    run = _load_def_inj(P6E, "run", pd=pd, y_all=y_all, g_all=g_all,
                        fit_predict=fit_predict, safe_auc=safe_auc)
    r = np.array([0.5, 0.4, 0.6, 0.9, 0.2, 0.7, 0.3, 0.1, 0.2])
    pooled, fold_map, fa = run(r)
    keys_ok = list(fold_map.keys()) == ["B"]
    fold_ok = (abs(fold_map["B"][0] - 1.0) < 1e-12) and (fold_map["B"][1] == 3)
    macro_ok = abs(np.mean(fa) - 1.0) < 1e-12
    # own OOF reference (definition (d)): fit on the other groups FOR THIS GROUP,
    # predict the held-out group; pooled = AUROC over the pooled OOF predictions.
    oof_ref = np.full(len(y_all), np.nan)
    for g_ in np.unique(g_all):
        te_ = g_all == g_
        tr_ = ~te_
        if len(np.unique(y_all[tr_])) < 2 or te_.sum() == 0:
            continue
        mu_ = np.nanmean(r[tr_])
        sd_ = np.nanstd(r[tr_])
        # ★ 2026-10-06：原写 `... or 1.0`，不防 nan（nan 为真 ⇒ 求值仍是 nan）——
        #   整列 NaN 时会崩在 sklearn 的通用报错里。此处显式判，报出是哪一组。
        if not (np.isfinite(mu_) and np.isfinite(sd_)) or sd_ == 0.0:
            raise ValueError("reference OOF 前置统计量非有限：group=%r nanmean=%s nanstd=%s"
                             "（该列在此训练子集上可能整列 NaN，常见于产物代次不一致）"
                             % (g_, mu_, sd_))
        mm_ = LogisticRegression(max_iter=1000).fit(((r[tr_] - mu_) / sd_).reshape(-1, 1), y_all[tr_])
        oof_ref[te_] = mm_.predict_proba(((r[te_] - mu_) / sd_).reshape(-1, 1))[:, 1]
    pooled_ref = ref_auc(y_all, oof_ref)
    pooled_ok = abs(pooled - pooled_ref) < 1e-9
    # fit_predict semantics fixture: sd>0; held-out only.
    # Separate def copies so the fixture labels (length 40) bind to ITS y_all.
    # (r is kept NaN-free: np.nanmean/np.nanstd are plumbed, but NaN in TRAIN
    # scores reaches sklearn validation and raises; the frozen 15 methods are
    # all finite, so that path is not exercised by this fixture.)
    rng = np.random.default_rng(11)
    rr = rng.normal(size=40).round(3)
    yy = (rng.random(40) < 0.5).astype(int)
    tr = np.arange(30)
    te = np.arange(30, 40)
    fp_a = _load_def_inj(P6E, "fit_predict", pd=pd, y_all=yy, g_all=np.zeros(40, int),
                         LogisticRegression=LogisticRegression)
    yy_perm = yy.copy()
    yy_perm[te] = 1 - yy_perm[te]           # labels of the HELD-OUT rows flipped
    fp_b = _load_def_inj(P6E, "fit_predict", pd=pd, y_all=yy_perm, g_all=np.zeros(40, int),
                         LogisticRegression=LogisticRegression)
    got = fp_a(tr, te, rr)
    mu = np.nanmean(rr[tr])
    sd = np.nanstd(rr[tr])
    # ★ 2026-10-06：同 t12/t13 的另两处，原写 `... or 1.0` 不防 nan。
    if not (np.isfinite(mu) and np.isfinite(sd)) or sd == 0.0:
        raise ValueError("独立参照的前置统计量非有限：nanmean=%s nanstd=%s"
                         "（该列在此训练子集上可能整列 NaN，常见于产物代次不一致）" % (mu, sd))
    m = LogisticRegression(max_iter=1000).fit(((rr[tr] - mu) / sd).reshape(-1, 1), yy[tr])
    refp = m.predict_proba(((rr[te] - mu) / sd).reshape(-1, 1))[:, 1]
    d_sem = float(np.abs(got - refp).max())
    got2 = fp_b(tr, te, rr)          # same labels, but y[te] is not an input to fit_predict
    leak = float(np.abs(got - got2).max())
    # sd==0 fixture
    r0 = rr.copy()
    r0[tr] = 5.0
    got0 = fp_a(tr, te, r0)
    m0 = LogisticRegression(max_iter=1000).fit(np.zeros((30, 1)), yy[tr])
    ref0 = m0.predict_proba(np.zeros((10, 1)))[:, 1]
    d_zero = float(np.abs(got0 - ref0).max())
    ok = (keys_ok and fold_ok and pooled_ok and macro_ok
          and d_sem <= 1e-9 and leak == 0.0 and d_zero <= 1e-9)
    obs = ("oracle fixture {A all-harmful, B mixed, C all-non-harmful}: fold_map keys=%s (single-class "
           "held-out folds excluded; A/C predictions still enter the POOLED OOF), fold AUC=%.6f n=%d, "
           "pooled=%.6f vs own OOF reference=%.6f, unweighted macro=%.6f; standardization: "
           "abs(exec - manual mean/std(ddof=0)->1.0 pipeline)=%.2e; y[test] permutation leak=%.2e "
           "(0 = no test-label leakage); sd=0 fallback diff=%.2e"
           % (list(fold_map.keys()),
              fold_map["B"][0] if "B" in fold_map else float("nan"),
              fold_map["B"][1] if "B" in fold_map else -1, pooled, pooled_ref, float(np.mean(fa)),
              d_sem, leak, d_zero))
    return ok, obs


def _warn_EC14():
    st = _p6e_state()
    if st["err"]:
        raise RuntimeError("p6e execution failed: " + st["err"])
    R = st["R"]
    FT = st["FT"]
    lo_u = []
    hi_u = []
    cross_w = []
    cross_u = []
    rng = np.random.default_rng(SEED)
    for name in R.method:
        fm = {g_: (a, n) for g_, a, n in FT[FT.method == name][["group", "fold_auc", "fold_n"]].values}
        keys = list(fm)
        us = []
        ws = []
        for _ in range(2000):
            idx = rng.integers(0, len(keys), len(keys))
            vals = [fm[keys[i]][0] for i in idx]
            w = np.array([fm[keys[i]][1] for i in idx], float)
            us.append(np.mean(vals))
            ws.append(np.average(vals, weights=w))
        lu, hu = (float(x) for x in np.percentile(us, [2.5, 97.5]))
        lw, hw = (float(x) for x in np.percentile(ws, [2.5, 97.5]))
        lo_u.append(lu)
        hi_u.append(hu)
        cross_w.append((lw <= 0.5) and (0.5 <= hw))
        cross_u.append((lu <= 0.5) and (0.5 <= hu))
    R2 = R.copy()
    R2["lo_u"] = lo_u
    R2["hi_u"] = hi_u
    gap = np.abs(R2.macro_mean - R2.size_weighted_macro)
    i_max = int(gap.idxmax())
    outside = int(((R2.macro_mean < R2.macro_lo) | (R2.macro_mean > R2.macro_hi)).sum())
    flips = R2[(np.array(cross_w) != np.array(cross_u))]
    ev_shift = float(max(np.abs(R2.macro_lo.values - R2.lo_u.values).max(),
                         np.abs(R2.macro_hi.values - R2.hi_u.values).max()))
    return ("p6e aggregation: macro_lo/hi = size-WEIGHTED bootstrap CI vs UNWEIGHTED macro_mean (column "
            "labels); numbers requested for the freeze decision",
            "max abs(macro_mean - size_weighted_macro) = %.4f (%s); methods with gap>0.005: %d/15; "
            "unweighted macro_mean outside its weighted CI: %d/15; CI covers 0.5: weighted %d/15 vs "
            "unweighted-%d/15 (flip: %s); max weighted-vs-unweighted CI endpoint shift %.4f"
            % (float(gap.max()), R2.method.iat[i_max], int((gap > 0.005).sum()), outside,
               int(sum(cross_w)), int(sum(cross_u)),
               ", ".join(flips.method.values) if len(flips) else "none", ev_shift))


def test_EC15():
    st6 = _p6d_state()
    st11 = _p11_state()
    st6e = _p6e_state()
    parts = []
    untested = []
    for tag, st in (("p6d", st6), ("p11", st11), ("p6e", st6e)):
        if st["ok"]:
            parts.append("%s EXECUTED [%s]" % (tag, "; ".join(l for _k, l in st["blocks"])))
        else:
            msg = "%s NOT EXECUTED (%s)" % (tag, st["err"])
            parts.append(msg)
            if L.IS_PATCHED and st["err"] and "absent" in st["err"]:
                untested.append(msg)
    ok = st6["ok"] and st11["ok"] and st6e["ok"]
    if not ok and L.IS_PATCHED and untested:
        L.emit("EC-15", "WARN",
               "census: every named path executed, or UNTESTED with the reason (never silent)",
               "; ".join(parts),
               "UNTESTED entries are files sim_fix does not copy into _patched; no PASS is claimed for them")
        return None
    return ok, "; ".join(parts)


# ====================================================================== runner
def _run(tid, expected, fn, note=""):
    try:
        r = fn()
    except Exception as e:
        return L.emit(tid, "FAIL", expected, "EXCEPTION %s: %s" % (type(e).__name__, e),
                      (note + " / " + traceback.format_exc().splitlines()[-1]).strip(" /"))
    if r is None:
        return None
    ok, obs = r
    return L.check(tid, ok, expected, obs, note)


def _run_warn(tid, expected, fn, note=""):
    try:
        expected2, obs = fn()
    except Exception as e:
        return L.emit(tid, "FAIL", expected, "EXCEPTION %s: %s" % (type(e).__name__, e),
                      (note + " / " + traceback.format_exc().splitlines()[-1]).strip(" /"))
    return L.emit(tid, "WARN", expected2, obs, note)


def _skip_if_absent(path, tid, what):
    if os.path.exists(path):
        return False
    if L.IS_PATCHED:
        L.emit(tid, "WARN", "UNTESTED: %s absent from the PATCHED tree" % what,
               "%s does not exist under EVS_PROD_DIR=%s; neither the code path nor the "
               "executed-vs-frozen comparison can run" % (path, PROD),
               "sim_fix LOADER_FILES does not copy this file -- UNTESTED is reported instead "
               "of a silent green")
        return True
    raise RuntimeError("named production file absent from the REAL tree: %s" % path)


def main():
    L.banner("EVALCORE - t13: p11 datapack / p6d first+second layer / p6e LQC")
    print("  mode            : %s" % ("PATCHED" if L.IS_PATCHED else "REAL"))
    print("  production tree : %s" % PROD)
    print("  frozen figures  : %s" % FIGD)
    print("  note: EC-01/EC-03/EC-12 execute production boxes verbatim; the 3-site model rebuild "
          "takes minutes.  All to_csv target writes are excluded or redirected to scratch, so the "
          "production tree is never written.")

    _run("EC-01",
         "EXECUTED p6d first layer verbatim (3-site STORE + model rebuild + per-task loop) == "
         "frozen p6d_outcome.csv (execution, not mere reconciliation; tol 1e-9)",
         test_EC01,
         "columns task_id..stable compared row-ordered; n_target exact, AUROCs tol 1e-9")

    _run("EC-01D",
         "EXECUTED p6d first layer verbatim -> the EXECUTED output satisfies the DEFINITIONS "
         "(d_auroc = auroc_target - auroc_source_internal; harmful = I(d<=-0.05); "
         "stable = I(abs(d)<=0.02); one source model per (organism,drug,source); "
         "task_id rule; auroc_target in [0,1]).  Value-independent: passes before AND after "
         "the corrected rerun.  EC-01 keeps the frozen-artifact comparison and is deferred.",
         test_EC01D,
         "code-level half of EC-01; the product-level half (== frozen p6d_outcome.csv) is a "
         "deferred gate, so this path is NOT moved wholesale to post_corrected_rerun")

    _run("EC-02",
         "frozen p6d_outcome.csv vs DEFINITIONS: d_auroc = auroc_target - auroc_source_internal; "
         "harmful = I(d<=-0.05); stable = I(abs(d)<=0.02); one source model per (organism,drug,source); "
         "task_id = org_drug_src2tgt; 202 tasks / 65 harmful / 49 stable",
         test_EC02,
         "reconciliation against the frozen CSV, not execution")

    _run("EC-03",
         "EXECUTED p6d per-method loop (285-317) verbatim == frozen p6d_signal_eval.csv, all 12 "
         "quantity columns incl. cluster-bootstrap CIs, LODO and in-file LOO",
         test_EC03,
         "the loop writes nothing (p6d:304 excluded); frame compared in memory")

    _run("EC-04",
         "independent midrank AUROC and step-AP (own implementations, sklearn/scipy cross-checked "
         "on tied fixtures) vs frozen harmful_auc / harmful_auprc over all 15 methods",
         test_EC04,
         "reference self-control in observed; reconciliation against the frozen CSV")

    _run("EC-05",
         "independent cluster bootstrap: resample organism-drug CLUSTERS with replacement, keep "
         "multiplicity (CS-08), statistic recomputed over the resample == frozen sp_lo/sp_hi and "
         "auc_lo/auc_hi (15 methods); a set-dedup variant must move the interval",
         test_EC05,
         "own Spearman/midrank-AUC implementations; per-method rng stream per the artifact provenance")

    if not _skip_if_absent(P11, "EC-06", "p11_figure_datapack.py"):
        _run("EC-06",
             "EXECUTED p11 fig3a block == frozen figures/data/fig3a_lodo.csv; category labels == "
             "CI-vs-0.5 rule on frozen AND executed frame; 25 evaluable folds / 138 tasks per the "
             "fold-composition definition",
             test_EC06,
             "fig3a is the categorical claim that flips with the aggregation convention (see EC-06b)")
    _run_warn("EC-06b",
              "characterisation: fig3a category labels under the p6f (unweighted) CI convention vs "
              "the p6e size-weighted CI convention -- how many of the 15 labels move",
              _warn_EC06b,
              "diagnostic only; p11 deliberately uses p6f's C_lo/C_hi")
    if not _skip_if_absent(P11, "EC-07", "p11_figure_datapack.py"):
        _run("EC-07",
             "EXECUTED p11 fig3c block == frozen fig3c_sudo_decomposition.csv; decomposition "
             "identities subset=B-A and aggregation=C-B; cross-artifact SUDO_AURCC pooled/macro "
             "agreement across p6f / p6e / p6d artifacts",
             test_EC07,
             "reconciliation + execution")
    if not _skip_if_absent(P11, "EC-08", "p11_figure_datapack.py"):
        _run("EC-08",
             "EXECUTED p11 fig4d block == frozen fig4d_mae.csv + 3 calibration files + fig4e; "
             "independent numpy MAE, scipy Spearman, bin group means, rho_ess",
             test_EC08,
             "MAE = mean(abs(est-true)) over finite estimates; calibration bins -inf/-0.10/-0.05/0/0.05/+inf")
    if not _skip_if_absent(P11, "EC-09", "p11_figure_datapack.py"):
        _run("EC-09",
             "EXECUTED p11 fig1c summary block == frozen fig1c_summary.csv; independent "
             "np.quantile(linear) / mean / sd(ddof=1) / P90-P10 / counts",
             test_EC09,
             "202 tasks 65 harmful 49 stable; P10/P50/P90 quoted in the brief reproduced")
    if not _skip_if_absent(P11, "EC-10", "p11_figure_datapack.py"):
        _run("EC-10",
             "EXECUTED p11 fig5/fig5b blocks == frozen fig5_subsample.csv + fig5b_across_n_range.csv; "
             "independent per-rep AUC replay (own midrank implementation) of all 15 methods x "
             "5 levels x 200 reps + paired full-100 rows",
             test_EC10,
             "RAW10 = frozen p10b_raw.csv; harm labels from frozen p6d_outcome.csv")
    if not _skip_if_absent(P11, "EC-11", "p11_figure_datapack.py"):
        _run("EC-11",
             "EXECUTED p11 remaining panels (fig1b, fig1d, fig1e, fig2, fig4e, figS6, tableS1, "
             "tableS2) == frozen files; independent fig1e recompute; tableS2 merge integrity",
             test_EC11,
             "13 panel files compared")

    if not _skip_if_absent(P6E, "EC-12", "p6e_01_lodo_fold_qc.py"):
        _run("EC-12",
             "EXECUTED p6e fold-composition block + aggregation loop verbatim == frozen "
             "p6e_lodo_fold_composition.csv / p6e_lodo_aggregation.csv / p6e_lodo_fold_auc_long.csv",
             test_EC12,
             "the two to_csv writes to the production dir are excluded; frames compared in memory")
    if not _skip_if_absent(P6E, "EC-13", "p6e_01_lodo_fold_qc.py"):
        _run("EC-13",
             "p6e run/fit_predict semantics on oracle fixtures: single-class held-out folds excluded "
             "from AUROC while their OOF predictions still enter the pooled AUROC, macro = UNWEIGHTED "
             "mean, standardization from TRAIN only, sd=0 -> 1.0 fallback, no y_test leakage",
             test_EC13,
             "fixture oracle: {A all-harmful, B mixed, C all-non-harmful} -> fold_map={B}, macro=1, "
             "pooled OOF == own reference.  np.nanmean/nanstd are plumbed but NaN in TRAIN scores "
             "reaches sklearn validation and raises; the frozen 15 methods are all finite")
    if not _skip_if_absent(P6E, "EC-14", "p6e_01_lodo_fold_qc.py"):
        _run_warn("EC-14",
                  "characterisation: p6e macro_lo/hi is a size-weighted bootstrap CI while macro_mean "
                  "(and the SUDO decomposition user) is unweighted -- numbers for the freeze decision",
                  _warn_EC14,
                  "documented inconsistency; p11 already refuses to use p6e's CI for Fig 3A")
    _run("EC-15",
         "census: every named path either EXECUTED verbatim (with block labels) or reported "
         "UNTESTED with the reason; never a silent green",
         test_EC15,
         "under PATCHED, files not copied by sim_fix are UNTESTED; under REAL they must exist")

    st11 = _STATE.get("p11")
    if st11 and st11.get("scratch"):
        import shutil
        shutil.rmtree(st11["scratch"], ignore_errors=True)


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        L.emit("EVALCORE-FATAL", "FAIL", "the module must run to completion",
               "EXCEPTION %s: %s" % (type(e).__name__, e),
               traceback.format_exc().splitlines()[-1])
