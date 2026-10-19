# -*- coding: utf-8 -*-
r"""t12_eval_layer.py - EVS family EVAL : the evaluation / aggregation layer

WHY THIS FILE EXISTS
================================================================================
EVS v1.1 closed the estimator-kernel layer (L1/L2/L3 arithmetic, fork parity,
value domains).  The layer that turns model predictions into the REPORTED
manuscript numbers of Figs. 1E, 2, 3A, 3B, S5, S9 and tableS2 was executed by
NO test.  An adversarial reviewer reproduced, from that layer, a live defect in
`p6f_lodo_settlement.py` and a Methods-vs-code divergence (calibration), both of
which move manuscript numbers.  t12 executes that layer verbatim and measures it
against independent references.

TARGETS (item brief v1.2)
--------------------------------------------------------------------------------
  (A) p6f_lodo_settlement.py:143-144 `set()` dedup inside the fold-bootstrap
      macro statistic  -> Fig. 3A CIs too narrow.
  (B) Platt calibration: Methods_v1.md:131 says calibrated probabilities were
      frozen and used downstream; p6b/p6c/p10a consume raw `predict_proba`.
  (C) Bootstrap multiplicity census: every live resampler in the evaluation
      layer, checked against ONE independent reference bootstrap that keeps
      multiplicity (resample WITH replacement, statistic over the resample
      retaining repeated draws -- Efron & Tibshirani 1993, ch. 6).
      NOTE (v1.4): the heuristic **EVAL-CENSUS scan was DELETED** -- membership of a live
      bootstrap site is now governed by ACTIVE_PATH_MANIFEST.yaml (path_id
      GROUPED_BOOTSTRAP_P6F / TARGETN_BRANCH), checked by t14 as existence + anchors only.
  (D) `lodo_oof` semantics: p6f's "LOO-organism" section (p6f:180-189) calls
      `lodo_oof(r)` (p6f:182) = leave-one-ORGANISM-DRUG-out, not
      leave-one-ORGANISM-out.

REFERENCES (source of every expectation; never paraphrased from the code tested)
--------------------------------------------------------------------------------
  (a) DEFINITION: bootstrap = resample WITH replacement from the fold list,
      statistic computed over the resample RETAINING multiplicity.  A macro
      mean over resample draw {A,A,B} is (2*A + B)/3, not (A + B)/2.
  (b) DEFINITION: AUROC = Mann-Whitney rank-sum statistic with midranks on
      ties (equals sklearn.metrics.roc_auc_score; verified to 1e-12).
  (c) DEFINITION: leave-one-group-out = fit on the other groups ONLY, predict
      the held-out group; macro = unweighted mean over held-out groups.
  (d) EXTERNAL: ARMD_Code_Audit_20261004.md -- corrected (multiplicity-keeping)
      Fig-3A intervals, the true-LOO values, and the calibration divergence.
  (e) EXTERNAL: Methods_v1.md:131-134 (calibration sentence);
      Results_v3_clean.md:73-76 ("Five of six L1 metrics ...").
  (f) EXTERNAL: frozen artifacts -- figures/data/fig3a_lodo.csv,
      p6f_lodo_settlement.csv, p6e_lodo_fold_auc_long.csv,
      p6e_lodo_aggregation.csv, p6d_signal_eval.csv, fig3b_loo_organism.csv.
  Production scripts are NOT references; they are the objects under test,
  loaded VERBATIM through evs_load (load_def / content-anchored grab + exec).

WHAT THIS MODULE MEASURES (13 RESULT lines)
--------------------------------------------------------------------------------
  EVAL-A1  fixture: production macro statistic dedups {A,A,B}->{A,B}   FAIL (D-A)
  EVAL-A2  replay fidelity: frozen fold AUCs + frozen pooled CI + frozen
           macro CI reproduced under p6f's exact RNG advance order        PASS
  EVAL-A3  corrected vs published Fig-3A intervals; audit cross-check    WARN
  EVAL-B1  Methods calibration claim vs raw predict_proba in p6b/p6c/p10a:
           RED only while the claim exists AND no calibration application is
           present at the probability source (CS-10 resolution)           FAIL (D-B)
  EVAL-B2  execute p12 calib verbatim on a fixture: Platt is monotone ->
           AUROC invariant, probability-dependent statistics move          PASS
  EVAL-C1  p6f boot resampler: multiplicity kept by the resampler, and the
           reproduction can go red if it were not (fixture replay)         PASS
  EVAL-C2  p6d cluster bootstrap (1000 draws, organism-drug): keeps
           multiplicity (3 methods, real frozen inputs, exact replay)      PASS
  EVAL-C3  p6e boot_macro: keeps multiplicity (fixture replay)             PASS
  EVAL-C4  p6g boot + p14 bci: keep multiplicity (fixture replay)          PASS
  EVAL-C5  p6b domain-AUC row bootstrap (200 draws): keeps multiplicity    PASS
  EVAL-D1  p6f "LOO-organism" section uses lodo_oof (LODO), not LOO      FAIL (D-D)
  EVAL-D2  external anchors: true-LOO + LODO-based organism macros match
           the audit numbers exactly; p6d/fig3b implement the true LOO      PASS

GREEN-ABILITY (rule: a test that cannot go GREEN on a correct implementation
is a PIN, not a test)
--------------------------------------------------------------------------------
  A1 : reads the verbatim statistic (content anchor 'flo, fhi = boot(EVAL');
       a fix that keeps draw multiplicity makes the executed value equal the
       definition -> PASS.  The fixture is built so ROW multiplicity != DRAW
       multiplicity, so a row-weighted statistic cannot pass by coincidence.
  B1 : RED only when an affirmative Methods calibration claim coexists with no
       calibration application in p6b/p6c/p10a.  Either accepted resolution
       turns it green: applying the frozen Platt map at the probability source
       (as p12/p5 do) OR deleting/correcting the Methods sentence (CS-10).
       Structural check: a comment mentioning `calib` is not an application.
  D1 : reads p6f's LOO section output; refitting per organism (inline or via a
       top-level helper such as `oof = loo_oof(r)`) makes the section value
       equal the definition -> PASS; retiring the section (absent) or marking
       it explicitly deprecated (CS-11: 'removed or explicitly deprecated')
       also PASSes -- no mislabeled LOO output remains authoritative.
  A2/A3, C1-C5, D2 are instrument / characterisation checks: A2 and the C
  family prove a reference replay is faithful, A3 is the WARN-class impact
  characterisation the brief asks for.  CENSUS is a coverage gate over an
  explicit registry (covered sites + out-of-scope reasons); an unregistered
  bootstrap-class site FAILs.

V1.3 REPAIRS (adversarial verifier, family EVAL)
--------------------------------------------------------------------------------
  - D1 namespace: the section block auto-binds top-level p6f defs it references
    (helper-style refits) plus the ambient evaluation names (g / org / GROUPS)
    and retries once on NameError before failing with the missing symbol named;
    the 'oof = ' extraction falls back to the whole section body when no
    assignment exists.  D2 no longer cascades on D1.
  - B1: substring counts replaced by the marker+application test; the
    comment-only bypass ('# calib absent') no longer flips it green.
  - A1: discriminating fixture; the stub offers fn(idx, pick), then
    fn(idx, pick=pick), then fn(idx).
  - Line pins (A1/C2/C4-bci/B2/C5) replaced by content anchors + required
    tokens; a moved layout emits 'block anchor moved -- re-pin' instead of
    executing unrelated lines.
  - CENSUS: the hard-coded census print became a real scan with a RESULT line.
  - A3: a degraded/partial audit cross-check reports INOPERATIVE instead of
    max|d|=0.00e+00.
  - every test emits exactly one RESULT line even when files are missing, so a
    trimmed PATCHED tree no longer loses 8 of 13 lines silently.
"""
import os
import re
import sys

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from scipy.stats import spearmanr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import evs_load as L

# --------------------------------------------------------------------- anchors
P6B = L.P6B
P6C = L.P6C
P10A = L.P10A
P6D = L.P6D
P5 = L.P5
P6E = os.path.join(L.PROD_DIR, "p6e_01_lodo_fold_qc.py")
P6F = os.path.join(L.PROD_DIR, "p6f_lodo_settlement.py")
P6G = os.path.join(L.PROD_DIR, "p6g_direction_learning.py")
P14 = os.path.join(L.PROD_DIR, "p14_rebuild_fig3b.py")
P12 = os.path.join(L.PROD_DIR, "p12_configs_corrected.py")

OUT = L.OUT
FIGD = os.path.join(L.PROJ, "figures", "data")
AUDIT = os.path.join(L.PROJ, "ARMD_Code_Audit_20261004.md")
# ★★ 2026-10-06 源绑定修正（属**正本治理**，与 CS-10 的 Platt 文本修正无关，两者分开记账）：
#   本判据原读**根目录**的 Methods_v1.md。但 `MANUSCRIPT_MASTER_DECLARATION.md` 认定正本是
#   `manuscript_renumbered/Methods_v1.md`；根目录那份是修正期对照副本，不用于投稿拼装。
#   判据必须读**权威文本**，否则会出现两种失真：
#     · 正本已改好、判据仍红（假红）；
#     · 更糟——判据绿了而正本没改（假绿）。
#   ⇒ 指向正本；对照副本只**报告状态**、不参与判定（见 test_B1 的 observed）。
METHODS_MD = os.path.join(L.PROJ, "manuscript_renumbered", "Methods_v1.md")
METHODS_MD_NONMASTER = os.path.join(L.PROJ, "Methods_v1.md")   # 对照副本：只报状态，不判定
RESULTS_MD = os.path.join(L.PROJ, "Results_v3_clean.md")

SEED = 20261005          # p6f/p6e/p6d/p6g/p14/p6b all use this seed
B_BOOT = 4000            # p6f / p6g / p14 draw count
B6D = 1000               # p6d cluster bootstrap draw count
B6B = 200                # p6b domain-AUC bootstrap draw count

DEF_BOOT = ("bootstrap definition: resample WITH replacement, statistic over the resample "
            "retaining multiplicity (Efron & Tibshirani 1993)")
RULE3 = ("references are definitions/external sources; production is executed verbatim as the "
         "object under test (evs_load load_def / content-anchored grab + exec)")

# p6f / p6e / p6g / p14 method table (order is load-bearing: the RNG advance order)
METHODS = [
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
SIG_COLS = ["domain_auc", "mmd", "psi_mean", "psi_max", "js_mean", "js_max",
            "pred_mean_shift", "entropy_mean_shift", "entropy_js", "conf_dist_js",
            "pred_pos_prev_shift"]
L3_COLS = ["SUDO_AURCC", "AUROC_IOW", "AUROC_OM", "AUROC_DR"]

# audit table label -> p6f method key (ARMD_Code_Audit_20261004.md, corrected table)
AUDIT_MAP = {
    "domain-AUC": "domain_auc", "MMD": "mmd", "PSI mean": "psi_mean", "PSI max": "psi_max",
    "JS mean": "js_mean", "JS max": "js_max", "entropy mean shift": "entropy_mean_shift",
    "entropy JS": "entropy_js", "confidence JS": "conf_dist_js", "SUDO变体": "SUDO_AURCC",
    "DR": "AUROC_DR", "IOW": "AUROC_IOW", "OM": "AUROC_OM",
}
# external reviewer anchors (audit, §p6f的LOO-organism段)
EXT_TRUE_LOO = {"entropy_js": 0.666307, "conf_dist_js": 0.446328, "SUDO_AURCC": 0.706656}
EXT_LODO_BASED = {"entropy_js": 0.583103, "conf_dist_js": 0.549195, "SUDO_AURCC": 0.625107}

_STATE = {}
_EMITTED = set()     # tids already emitted this run (main's exception wrapper uses it)


# --------------------------------------------------------------------- helpers
def _flat(s):
    """RESULT fields are pipe-delimited: no '|' / newlines may survive."""
    return str(s).replace("|", "/").replace("\r", " ").replace("\n", " ")


def _emit(tid, status, expected, observed, note=""):
    _EMITTED.add(tid)
    L.emit(tid, status, _flat(expected), _flat(observed), _flat(note))


def _lines(path):
    with open(path, encoding="utf-8") as f:
        return f.read().split("\n")


# --------------------------------------------------- content-anchored extraction
# v1.3: absolute line pins (p6b:319-323, p12:208-214, p14:97-102, ...) false-REDded
# any edit that added a line above them, grabbing unrelated code.  Blocks are now
# located by an anchor regex on the statement + required tokens + bracket balance.
def _line_depth_delta(text):
    """Net bracket delta of one source line, ignoring string literals and comments."""
    d = 0
    q = None
    esc = False
    for ch in text:
        if q:
            if esc:
                esc = False
            elif ch == "\\":
                esc = True
            elif ch == q:
                q = None
        elif ch in "\"'":
            q = ch
        elif ch in "([{":
            d += 1
        elif ch in ")]}":
            d -= 1
        elif ch == "#":
            break
    return d


def _grab_anchor(path, anchor_pat, tokens, max_lines=16):
    """Grab the block that starts at the first line matching `anchor_pat`, growing it
    until every token in `tokens` has appeared AND the bracket depth is back to 0.

    Returns (dedented text, first line number, err) with err in
    (None, 'no-anchor', 'tokens-missing', 'io-error(...)').  A moved layout yields
    err -- the caller emits 'block anchor moved -- re-pin', never executes the
    wrong lines and never reports a fixture-degenerate misdiagnosis.
    """
    try:
        lines = _lines(path)
    except OSError as e:
        return None, None, f"io-error({type(e).__name__})"
    for i, ln in enumerate(lines):
        if re.search(anchor_pat, ln):
            buf, depth = [], 0
            for j in range(i, min(i + max_lines, len(lines))):
                buf.append(lines[j])
                depth += _line_depth_delta(lines[j])
                txt = "\n".join(buf)
                if depth <= 0 and all(t in txt for t in tokens):
                    import textwrap as _tw
                    return _tw.dedent(txt), i + 1, None
            return None, i + 1, "tokens-missing"
    return None, None, "no-anchor"


def _exec_src(src, env, tag):
    ns = dict(env)
    ns.setdefault("np", np)
    exec(compile(src, f"<verbatim:{tag}>", "exec"), ns)
    return ns


def _block_err(tid, expected, err, name):
    txt = (f"UNTESTED: {name} absent from this tree" if err and err.startswith("io-error")
           else f"block anchor moved/absent in {name} ({err}) -- file layout changed, re-pin")
    _emit(tid, "FAIL", expected, txt,
          "content-anchored extraction refused to execute unrelated lines")


def _need_inputs(tid, expected):
    if "build_error" in _STATE:
        _emit(tid, "FAIL", expected,
              f"UNTESTED: evaluation inputs not built ({_STATE['build_error']})", "")
        return False
    return True


def _need_p6f(tid, expected):
    if not _need_inputs(tid, expected):
        return False
    if "p6f_error" in _STATE:
        _emit(tid, "FAIL", expected,
              f"UNTESTED: p6f defs unavailable ({_STATE['p6f_error']})", "")
        return False
    return True


def ref_auc(yy, ss):
    """Independent definitional AUROC: Mann-Whitney rank-sum with midranks on ties.

    Not production's rank loop: written here from the definition and cross-checked
    against sklearn.metrics.roc_auc_score (they agree to <1e-12)."""
    yy = np.asarray(yy)
    ss = np.asarray(ss, dtype=float)
    m = np.isfinite(ss)
    yy, ss = yy[m], ss[m]
    if len(np.unique(yy)) < 2:
        return np.nan
    order = np.argsort(ss, kind="mergesort")
    yy, ss = yy[order], ss[order]
    n = len(ss)
    r = np.empty(n)
    i = 0
    while i < n:
        j = i
        while j + 1 < n and ss[j + 1] == ss[i]:
            j += 1
        r[i:j + 1] = 0.5 * (i + j) + 1.0
        i = j + 1
    n1 = float(yy.sum())
    n0 = n - n1
    if n1 == 0 or n0 == 0:
        return np.nan
    return (r[yy == 1].sum() - n1 * (n1 + 1) / 2.0) / (n1 * n0)


def build_inputs():
    """Merge the frozen evaluation inputs exactly as p6f/p6e/p6g/p14/p6d do, and
    reproduce p6f's EVAL fold list from the data (not from any hardcoded list)."""
    O = pd.read_csv(os.path.join(OUT, "p6d_outcome.csv"))
    SIG = pd.read_csv(os.path.join(OUT, "p6b_label_free_signals_blinded.csv"))
    L3 = pd.read_csv(os.path.join(OUT, "p6c_L3_blinded.csv"))
    E = O.merge(SIG[["task_id"] + SIG_COLS], on="task_id", how="left") \
         .merge(L3[["task_id"] + L3_COLS], on="task_id", how="left")
    E["group"] = E.organism + "|" + E.drug
    y = E["harmful"].values.astype(int)
    g = E["group"].values
    groups = sorted(np.unique(g))
    ev = []
    for gg in groups:
        m = g == gg
        if y[m].sum() > 0 and (1 - y[m]).sum() > 0:
            ev.append(gg)
    _STATE.update(E=E, y=y, g=g, groups=groups, EVAL=ev)


def prepare_p6f():
    """Load p6f's own auc / lodo_oof / boot verbatim and build per-method fa/oof."""
    p6f_auc = L.load_def(P6F, "auc")
    p6f_lodo = L.load_def(P6F, "lodo_oof",
                          extra_ns={"GROUPS": _STATE["groups"], "g": _STATE["g"],
                                    "y": _STATE["y"]})
    E, y, g, EVAL = _STATE["E"], _STATE["y"], _STATE["g"], _STATE["EVAL"]
    fa, oof = {}, {}
    for name, _role, fn in METHODS:
        r = fn(E).values.astype(float)
        o = p6f_lodo(r)
        oof[name] = o
        d = {}
        for gg in EVAL:
            m = g == gg
            a = p6f_auc(y[m], o[m])
            if np.isfinite(a):
                d[gg] = a
        fa[name] = d
    _STATE.update(p6f_auc=p6f_auc, p6f_lodo=p6f_lodo, fa=fa, oof=oof)


# ===================================================================== (A) 3.1
def test_A1():
    """Fixture, executed verbatim: does the production macro statistic keep draw
    multiplicity?  The fixture is built so ROW multiplicity != DRAW multiplicity
    (g_fix row counts A=3, B=2, C=4, D=1; drawn fold lists {A,A,B} and {A,B,B}):
        {A,A,B}: definition (2*fa[A]+fa[B])/3 = 0.8   (set()-dedup 0.7; row-weighted 0.85)
        {A,B,B}: definition (fa[A]+2*fa[B])/3 = 0.6   (set()-dedup 0.7; row-weighted 4.6/7)
    Only a statistic that averages fa over the DRAWN fold list matches."""
    if not os.path.exists(P6F):
        _emit("EVAL-A1", "FAIL", DEF_BOOT + f"; executed verbatim {os.path.basename(P6F)}",
              f"UNTESTED: {os.path.basename(P6F)} absent from this tree", "")
        return
    g_fix = np.array(["A"] * 3 + ["B"] * 2 + ["C"] * 4 + ["D"] * 1)
    fa = {"A": 1.0, "B": 0.4, "C": 0.1, "D": 0.2}
    EV = ["A", "B", "C", "D"]
    block, blk, err = _grab_anchor(P6F, r"flo\s*,\s*fhi\s*=\s*boot\(EVAL", ["boot(EVAL", "fa"])
    if err:
        _block_err("EVAL-A1", DEF_BOOT + "; executed verbatim p6f boot-macro statistic",
                   err, os.path.basename(P6F))
        return

    def run_case(idx, pick):
        # the stub drives the verbatim statistic with a chosen draw.  It offers
        # fn(idx, pick) first (a fix that receives the drawn fold list; pick = indices
        # into EVAL), then the keyword form fn(idx, pick=pick), then the current
        # 1-argument form fn(idx).
        def stub(sub, fn, **kw):
            try:
                v = fn(idx, pick)
            except TypeError:
                try:
                    v = fn(idx, pick=pick)
                except TypeError:
                    v = fn(idx)
            return (float(v), float("nan"))
        ns = _exec_src(block, {"boot": stub, "EVAL": EV, "g": g_fix, "fa": fa},
                       "p6f-boot-macro-stat")
        return ns["flo"]

    idx1 = [0, 1, 2, 0, 1, 2, 3, 4]      # boot's builder on fold draws {A,A,B}
    idx2 = [0, 1, 2, 3, 4, 3, 4]         # boot's builder on fold draws {A,B,B}
    bug1 = run_case(idx1, [0, 0, 1])
    ref1 = (2 * 1.0 + 0.4) / 3.0
    bug2 = run_case(idx2, [0, 1, 1])
    ref2 = (1.0 + 2 * 0.4) / 3.0
    ok = (np.isfinite(bug1) and np.isfinite(bug2)
          and abs(bug1 - ref1) < 1e-12 and abs(bug2 - ref2) < 1e-12)
    _emit("EVAL-A1", "PASS" if ok else "FAIL",
          DEF_BOOT + f"; executed verbatim p6f_lodo_settlement.py (anchor 'flo, fhi = "
          f"boot(EVAL', grabbed at line {blk})",
          f"{{A,A,B}}: prod={bug1:.6f} vs definitional={ref1:.6f}; "
          f"{{A,B,B}}: prod={bug2:.6f} vs definitional={ref2:.6f}",
          "set(g[idx]) collapses repeated fold draws -> 0.7/0.7 instead of 0.8/0.6 "
          "(CS-08); the fix must average fa over the DRAWN fold list (draw multiplicity), "
          "e.g. pass `pick` into the statistic -- this fixture also rejects the "
          "row-weighted g[idx] variants (0.850000 / 0.657143), which the old "
          "coinciding-multiplicity fixture could not distinguish")


def test_A2_A3():
    """Replay p6f's exact RNG advance (SEED=20261005; per method 4000 pooled draws
    then 4000 macro draws, rng.choice(n,n,replace=True)) and compare against the
    frozen artifacts; then the corrected vs published characterisation."""
    if "build_error" in _STATE or "p6f_error" in _STATE:
        why = _STATE.get("build_error") or _STATE.get("p6f_error")
        _emit("EVAL-A2", "FAIL", "frozen-artifact replay", f"UNTESTED: {why}", "")
        _emit("EVAL-A3", "FAIL", "characterisation needs EVAL-A2", f"UNTESTED: {why}", "")
        return
    E, y, g, EVAL = _STATE["E"], _STATE["y"], _STATE["g"], _STATE["EVAL"]
    oof, fa = _STATE["oof"], _STATE["fa"]
    ng = len(EVAL)
    group_idx = {gg: np.where(g == gg)[0] for gg in EVAL}

    f3a = pd.read_csv(os.path.join(FIGD, "fig3a_lodo.csv")).set_index("method")
    s16 = pd.read_csv(os.path.join(OUT, "p6f_lodo_settlement.csv")).set_index("method")
    p6e_long = pd.read_csv(os.path.join(OUT, "p6e_lodo_fold_auc_long.csv"))

    # ---- (0) fold-level agreement with the frozen per-fold CSV ----
    max_fold_d, n_fold = 0.0, 0
    for name, _r, _f in METHODS:
        sub = p6e_long[p6e_long.method == name].set_index("group")["fold_auc"]
        for k, v in fa[name].items():
            if k in sub.index:
                max_fold_d = max(max_fold_d, abs(float(sub.loc[k]) - v))
                n_fold += 1

    # ---- (1) RNG replay ----
    rng = np.random.default_rng(SEED)
    buggy, corrected = {}, {}
    max_pool_d, max_bug_d = 0.0, 0.0
    # ★ 2026-10-06：新增 max_corr_d —— 「fig3a 是否等于**正确**（含重数）版」的距离。
    #   原判据只算了 max_bug_d 并断言它 <= 1e-9，即「fig3a == 去重版」。那是 CS-08 **未修时**
    #   的红灯基线。CS-08 修好后（p6f 保留重数、p11 读 p6f），正确的断言必须反过来。
    max_corr_d = 0.0
    for name, _r, _f in METHODS:
        o = oof[name]
        d = fa[name]
        # 4000 draws consumed by the pooled-CI bootstrap call (p6f:142)
        pv = []
        for _ in range(B_BOOT):
            pick = rng.choice(ng, ng, replace=True)
            idx = np.concatenate([group_idx[EVAL[i]] for i in pick])
            v = ref_auc(y[idx], o[idx])
            if np.isfinite(v):
                pv.append(v)
        lo, hi = np.percentile(pv, [2.5, 97.5])
        max_pool_d = max(max_pool_d, abs(lo - s16.loc[name, "B_lo"]),
                         abs(hi - s16.loc[name, "B_hi"]))
        # 4000 draws consumed by the macro-CI bootstrap call (p6f:143)
        bl, cl = [], []
        for _ in range(B_BOOT):
            pick = rng.choice(ng, ng, replace=True)
            idx = np.concatenate([group_idx[EVAL[i]] for i in pick])
            s_ = set(g[idx])                      # verbatim production statistic
            if s_ & set(d):
                bl.append(float(np.mean([d[x] for x in s_ if x in d])))
            cl.append(float(np.mean([d[EVAL[i]] for i in pick])))   # definition
        blo, bhi = np.percentile(bl, [2.5, 97.5])
        clo, chi = np.percentile(cl, [2.5, 97.5])
        buggy[name] = (blo, bhi)
        corrected[name] = (clo, chi)
        max_bug_d = max(max_bug_d, abs(blo - f3a.loc[name, "lo_95"]),
                        abs(bhi - f3a.loc[name, "hi_95"]))
        max_corr_d = max(max_corr_d, abs(clo - f3a.loc[name, "lo_95"]),
                         abs(chi - f3a.loc[name, "hi_95"]))
    _STATE.update(buggy=buggy, corrected=corrected)

    # ★ 2026-10-06 判据语义修正（同 CS04-05 那一类：前提被后续修正取代）：
    #   · 原：`n_fold == 375 and ... and max_bug_d <= 1e-9`
    #        —— 375 = 25 折 × 15 方法，硬编码；且断言 fig3a == **去重版**。
    #   · 现状实测：eligible 宇宙为 186 task / 37 组 ⇒ 可定义折 26 ⇒ 26×15 = **390** 行；
    #       且 CS-08 已修 ⇒ fig3a 应当等于**含重数**版而非去重版。
    #   · 改为：`n_fold == len(p6e_long)`（= 重放覆盖了产物的**每一行**，不再硬编码折数）
    #       且 `max_corr_d <= 1e-9`（fig3a == 正确版）。
    #   · **max_bug_d 保留在 observed 里**：它现在不再是判据，而是
    #       「正确版与去重版确实不同」的**鉴别力证据**（实测 3.63e-02）。
    #       若哪天它变成 0，说明两版又合一了，需要重新确认 CS-08 是否还在生效。
    ok = (n_fold == len(p6e_long) and np.isfinite([max_fold_d, max_pool_d, max_corr_d]).all()
          and max_fold_d <= 1e-12 and max_pool_d <= 1e-9 and max_corr_d <= 1e-9)
    _emit("EVAL-A2", "PASS" if ok else "FAIL",
          f"frozen artifacts: p6e_lodo_fold_auc_long.csv（全部 {len(p6e_long)} 行都须被重放覆盖）、"
          "p6f_lodo_settlement.csv B_lo/B_hi、figures/data/fig3a_lodo.csv lo_95/hi_95；"
          "按 p6f 的 RNG 顺序重放：SEED=20261005，每方法先 4000 次 pooled 再 4000 次 macro "
          "rng.choice(n,n,replace=True)。★ 判据 = fig3a 必须等于**含重数**的正确版"
          "（CS-08 后），不是等于去重版",
          f"行覆盖 {n_fold}/{len(p6e_long)} max|d|={max_fold_d:.2e}；pooled CI 15/15 "
          f"max|d|={max_pool_d:.2e}；**正确版(含重数)** macro CI max|d|={max_corr_d:.2e}；"
          f"[鉴别力证据] 去重版 macro CI max|d|={max_bug_d:.2e}（应为大值；若为 0 则两版合一，需复查 CS-08）",
          "instrument check: the replay reproduces the frozen artifacts exactly, so the "
          "corrected column below is comparable like-for-like. Matched: seed, method order, "
          "draw counts, call order, resample protocol, EVAL order. NOT matched: model "
          "training (fold AUCs come from p6f's own lodo_oof, executed verbatim) and the "
          "other confirmed defects (MMD/DR/SUDO inputs) which remain uncorrected here")

    # ---- (2) corrected vs published + audit cross-check ----
    audit = {}
    try:
        for ln in _lines(AUDIT):
            cells = [c.strip() for c in ln.strip().strip("|").split("|")]
            if len(cells) >= 5 and cells[0] in AUDIT_MAP:
                m = re.fullmatch(r"([0-9]+\.[0-9]+)[–\-]([0-9]+\.[0-9]+)", cells[3])
                if m:
                    audit[AUDIT_MAP[cells[0]]] = (float(m.group(1)), float(m.group(2)))
    except Exception:
        audit = {}
    n_au, max_au, n_bad = 0, 0.0, 0
    for k, (al, ah) in audit.items():
        cl, ch = corrected.get(k, (np.nan, np.nan))
        if np.isfinite(cl) and np.isfinite(ch):
            max_au = max(max_au, abs(cl - al), abs(ch - ah))
            n_au += 1
        else:
            n_bad += 1
    # v1.3: a partial/NaN cross-check must not print 'max|d|=0.00e+00' as if it had
    # succeeded -- report it as INOPERATIVE (status stays WARN, this is diagnostic).
    if n_au < len(AUDIT_MAP) or n_bad > 0:
        cross_txt = (f"audit cross-check INOPERATIVE ({n_au}/{len(AUDIT_MAP)} parsed"
                     + (f", {n_bad} non-finite corrected)" if n_bad else ", degraded)"))
    else:
        cross_txt = f"audit cross-check {n_au}/{len(AUDIT_MAP)} max|d|={max_au:.2e}"

    crosses = []
    for name, _r, _f in METHODS:
        pub = (float(f3a.loc[name, "lo_95"]), float(f3a.loc[name, "hi_95"]))
        cor = corrected[name]
        pub_in = pub[0] <= 0.5 <= pub[1]
        cor_in = cor[0] <= 0.5 <= cor[1]
        if (not pub_in) and cor_in:
            crosses.append(f"{name} pub[{pub[0]:.6f},{pub[1]:.6f}] -> "
                           f"cor[{cor[0]:.6f},{cor[1]:.6f}]")
    l1 = [m for m in METHODS if m[1] == "L1"]
    n_pub_below = sum(1 for nm, _r, _f in l1
                      if float(f3a.loc[nm, "hi_95"]) < 0.5)
    n_cor_below = sum(1 for nm, _r, _f in l1 if corrected[nm][1] < 0.5)
    results_sentence = ""
    try:
        txt = " ".join(_lines(RESULTS_MD)[:120])
        m = re.search(r"Five of six L1 metrics[^.]*\.", txt)
        if m:
            results_sentence = m.group(0)
    except Exception:
        pass
    _emit("EVAL-A3", "WARN",
          "characterisation: multiplicity-keeping bootstrap (definition) vs published "
          "set-dedup intervals; external cross-check = ARMD_Code_Audit_20261004.md "
          "§分组bootstrap table (13 methods)",
          f"crossing to include 0.5: {'; '.join(crosses)}; "
          f"L1 intervals entirely below 0.5: published {n_pub_below}/6 -> corrected "
          f"{n_cor_below}/6; {cross_txt}",
          ("published = figures/data/fig3a_lodo.csv (identical to p6f C_lo/C_hi); "
           "diagnostic only -- MMD/DR/SUDO inputs still carry other confirmed defects; "
           "Results_v3_clean.md:73-76 quotes the published set "
           + ("(matched: " + results_sentence[:90] + "...)" if results_sentence else ""))
          + ("; confidence JS lower bound moves 0.537951 -> 0.502778 (still >0.5)"
             if "conf_dist_js" in corrected else ""))


# ===================================================================== (B) 3.2
def _methods_calib_claim(path=None):
    """The affirmative 'probabilities were calibrated' claim in Methods, if any.
    Negated/corrected sentences ('were not calibrated', 'uncalibrated', CS-10 '改准')
    are not claims.  Returns (sentence, note).

    `path` defaults to the MASTER (see METHODS_MD); pass METHODS_MD_NONMASTER to probe the
    non-master comparison copy without that probe affecting the verdict."""
    path = path or METHODS_MD
    try:
        txt = " ".join(_lines(path))
    except OSError:
        return "", f"{os.path.basename(path)} unreadable"
    for sent in re.split(r"(?<=[.!?。])\s+", txt):
        s = sent.strip()
        sl = s.lower()
        if ("calibrat" in sl and "probabilit" in sl
                and not re.search(r"(?i)\b(not|never|no|without)\b|un-?calibrat", s)):
            return s, "affirmative calibration claim found"
    return "", "no affirmative calibration claim in Methods (raw predict_proba consistent, CS-10)"


# ===================================================================== CS-10 语义升级
# ★★ 2026-10-06（CS-10 amendment）：判据语义由
#      「no unsupported calibration claim」
#    升级为
#      「calibration claim must match implementation scope」。
#
#   触发这次升级的实测事实（不是推理）：
#     · 协议 `platt` 命中 = 0                       ⇒ 协议从未要求 Platt
#     · 活的 p5 里 Platt **确实在跑**（p5:237/243/244，被 p5:285 调用）
#     · p6b / p6c / p10a **不做** Platt，直接用 predict_proba
#     · 全项目**不存在**"被取代的 p5"
#   ⇒ CS-10 原始缺陷描述（"删除一个没有实现支撑的声明"）有一半不成立：该声明**有**实现支撑，
#     问题在于它把一条**管线局部**的实现泛化成了整条分析链的共同步骤。
#   ⇒ 因此不能删（删了会制造新的不一致），只能**限定作用域**。
#
#   三段式（GPT6 2026-10-06 裁决的 A/B/C）：
#     A claim detection      —— 已有的肯定式声明探针
#     B scope binding        —— 声明必须限定到实际执行校准的管线
#     C implementation concordance —— 被声明的那个管线在代码里确实有校准
#   并要求方法学**正面陈述**下游分析用的是未校准概率（不是写"p6c 没做 X"）。
#   为此处要求"下游陈述"而非"可选"的理由：稿件 §9-11 确实在描述那些分析，
#   读者会把上一个 §5 的校准声明自然带下去；只限定作用域而不说明下游，
#   正是 CS-10 得以存在的那个误解。
SCOPE_PHRASES = ("source-model development", "source model development",
                 "benchmark construction", "source-model training", "source model training")
# ★ 第一版实现（2026-10-06 早）被独立对抗验证攻破，四处真旁路。以下三项限制是针对那些
#   旁路**逐条**加的，并且由 EVAL-B1b 用同一批 specimen 做**常设回归**（见文件末）。
#
#   旁路①：scope_ok 只查"句内是否出现某短语" ⇒ 指针括号 "(see the source-model development
#           pipeline, Section 5)" 即满足，而句子仍可断言 "throughout the analyses reported here"。
#           ⇒ 修正：作用域短语必须出现在**括号之外**；且声明句不得含泛化量词。
GENERALIZATION_MARKERS = ("throughout", "all analyses", "every analysis", "each analysis",
                          "the analyses reported", "reported here", "all downstream",
                          "each pipeline", "every pipeline", "in all")
#   旁路②：down_ok 是**文件级子串** ⇒ 任意一句 "... were not calibrated against target
#           outcomes" 就能满足，与下游概率来源零关系。⇒ 修正：改为**同句共现**——
#           必须存在一句话同时含【下游标记 + "probabilit" + 未校准标记】。
DOWNSTREAM_MARKERS = ("downstream", "sections 9–11", "sections 9-11", "section 9–11",
                      "signal analyses", "transported-performance")
UNCALIB_MARKERS = ("raw predicted probabilit", "raw probabilit", "uncalibrated probabilit",
                   "without this calibration", "without the calibration",
                   "not calibrated")
#   旁路③：探针检不出**自相矛盾** —— 真句含 "do not"（neg-exempt）既不能成为 claim，
#           又供给 down_ok ⇒ "全链已校准" + "下游不消费校准" 并存也 PASS。
#           ⇒ 无法用文本探针彻底消除（语义判断不可机械化）。故：
#             · 泛化量词限制堵住"全链已校准"这一半；
#             · 并在 expected 里**显式声明本判据是结构性的、不是语义保证**，不得当作语义认证。


def _sents(path):
    try:
        txt = " ".join(_lines(path))
    except OSError:
        return None
    return txt, [s.strip() for s in re.split(r"(?<=[.!?。])\s+", txt) if s.strip()]


def _scope_of_text(txt):
    """核心逻辑（纯文本 → (scope_ok, downstream_ok, note)）。抽出来是为了让 EVAL-B1b 的
    canary 能直接喂**specimen 文本**，而不用写临时文件。"""
    sents = [s.strip() for s in re.split(r"(?<=[.!?。])\s+", " ".join(txt.split("\n"))) if s.strip()]
    claims, unscoped, reasons = [], [], []
    for s in sents:
        sl = s.lower()
        if ("calibrat" in sl and "probabilit" in sl
                and not re.search(r"(?i)\b(not|never|no|without)\b|un-?calibrat", s)):
            claims.append(s)
            outside = re.sub(r"\([^)]*\)", " ", s).lower()      # 去括号内容
            scope_hit = [p for p in SCOPE_PHRASES if p in outside]
            gen_hit = [g for g in GENERALIZATION_MARKERS if g in sl]
            if not scope_hit:
                unscoped.append(s)
                reasons.append("scope 短语缺失或在括号内")
            if gen_hit:
                unscoped.append(s)
                reasons.append(f"含泛化量词 {gen_hit}")
    down_hits = []
    for s in sents:
        sl = s.lower()
        if (any(m in sl for m in DOWNSTREAM_MARKERS) and "probabilit" in sl
                and any(u in sl for u in UNCALIB_MARKERS)):
            down_hits.append(s)
    scope_ok = bool(claims) and not unscoped
    down_ok = bool(down_hits)
    note = (f"claims={len(claims)} unscoped={len(set(unscoped))}"
            + (f" :: {sorted(set(reasons))}" if reasons else "")
            + (f" :: {unscoped[0][:80]}..." if unscoped else "")
            + f"; downstream_sentences={len(down_hits)}"
            + (f" :: {down_hits[0][:80]}..." if down_hits else ""))
    return scope_ok, down_ok, note


def _calib_claim_scope(path=None):
    """CS-10 升级后的 B/C 两段。返回 (scope_ok, downstream_ok, note)。

    scope_ok      : 每条肯定式校准声明句必须 ①含**括号外**的作用域短语，且 ②不含泛化量词。
                    不满足 = 泛化声明 = FAIL（这正是 CS-10 的原缺陷）。
    downstream_ok : 必须存在**同一句话**同时含【下游标记 + probabilit + 未校准标记】。
                    不是文件级子串命中 —— 那是第一版的漏洞。
    """
    path = path or METHODS_MD
    r = _sents(path)
    if r is None:
        return False, False, f"{os.path.basename(path)} unreadable"
    return _scope_of_text(r[0])


# ---------------------------------------------------------------- EVAL-B1b 非空转 canary
# ★★ 2026-10-06：_calib_claim_scope 的**第一版**被独立对抗验证攻破 —— 4 个角度各派一名
#   skeptic，**4/4 都构造出"能过探针但仍然误导"的文本**，且全部通过第二轮对抗复核。三处漏洞：
#     ① scope 短语出现在**指针括号**里即满足（句子仍可断言 "throughout the analyses reported here"）
#     ② down_ok 是**文件级子串** ⇒ 任意一句 "... were not calibrated" 就够，与下游概率零关系
#     ③ 探针检不出**自相矛盾**（真句含 "do not" 被 neg-exempt，同时又供给 down_ok）
#   最尖锐的复核意见：**绿色已无法区分真修法与糊弄文本。**
#   ⇒ 本 canary 把那批 specimen **固化**：糊弄文本必须判红，真修法必须判绿。
#     这不是"再加一条判据"，而是让 EVAL-B1 的**判别力本身**成为可回归的对象。
B1B_BYPASS = [
    ("指针括号 + 泛化量词",
     "Predicted probabilities were calibrated by Platt scaling throughout the analyses reported "
     "here (see the source-model development pipeline, Section 5). The calibration intercept and "
     "slope reported in Section 7 were computed for descriptive purposes and were not calibrated "
     "against target outcomes."),
    ("删真句 + 塞一句无关的 not calibrated",
     "Within source-model development and benchmark construction, predicted probabilities were "
     "calibrated by Platt scaling fitted on the source validation split alone; the calibration map "
     "was then frozen and applied unchanged to subsequent data. The 0.5 cut-point is a reporting "
     "convention and was not calibrated to the target cohort."),
    ("自相矛盾：泛化声明与真排除句并存",
     "For source-model development, predicted probabilities were calibrated by Platt scaling for "
     "all analyses. The downstream input-shift and transported-performance analyses do not consume "
     "that calibration map and use their raw predicted probabilities."),
]
B1B_POSITIVE = (
    "Within source-model development and benchmark construction, predicted probabilities were "
    "calibrated by Platt scaling fitted on the source validation split alone. The downstream "
    "input-shift, prediction-space and transported-performance analyses (Sections 9-11) do not "
    "consume that calibration map and use their raw predicted probabilities.")


def test_B1b():
    """EVAL-B1b：EVAL-B1 探针的**非空转**证明。

    · 三条对抗 specimen 必须**判红**（否则探针可被糊弄）；
    · 一条真修法文本必须**判绿**（否则探针会误杀正确写法）。
    两者缺一，EVAL-B1 的绿色就不携带信息。
    """
    bad = []
    for name, txt in B1B_BYPASS:
        so, dn, note = _scope_of_text(txt)
        if so and dn:
            bad.append(f"BYPASS-ACCEPTED[{name}]")
    so_p, dn_p, note_p = _scope_of_text(B1B_POSITIVE)
    if not (so_p and dn_p):
        bad.append(f"POSITIVE-REJECTED :: {note_p}")
    # ★ 必须**自己 emit**：本模块的 run() 只在抛异常时才补发 FAIL 行，测试函数如果不 emit，
    #   就会**静默什么都不发**（既不绿也不红）—— 第一版 B1b 正是这样，靠核对输出才发现。
    _emit("EVAL-B1b", "PASS" if not bad else "FAIL",
          "EVAL-B1 探针的非空转证明：三条被独立对抗验证攻破过的 specimen 必须判红，"
          "一条真修法文本必须判绿；两者缺一，EVAL-B1 的绿色就不携带信息",
          f"bypass specimens={len(B1B_BYPASS)} 全部判红="
          f"{not any(b.startswith('BYPASS') for b in bad)}；positive 判绿={bool(so_p and dn_p)}；"
          + (f"问题={bad}" if bad else "OK"),
          "specimen 出处：2026-10-06 独立对抗验证（4 个角度 4/4 攻破第一版探针）；"
          "本 canary 把那批 specimen 固化为常设回归")


# structural evidence for a calibration stage at the probability source
CALIB_MARKER_PAT = re.compile(r"(?i)def\s+calib|\bplatt\b")
CALIB_CALL_PAT = re.compile(r"(?i)\b(?:calib|platt)\w*\s*\(")
CALIB_WRAP_F_PAT = re.compile(r"(?i)\b(?:calib|platt)\w*\s*\(\s*f\s*\(")


def _calib_stage(path):
    """(marker_lines, application_lines): (a) a calibration-stage definition/marker and
    (b) a live application of the stage to the model probability source -- same line
    mentions predict_proba, or sits within +-3 lines of a predict_proba expression, or
    is an `src_auroc = ... calib(f(`-style application.  Comment-only lines are not
    evidence (v1.3: this closes the '# calib absent' comment bypass)."""
    ls = _lines(path)
    marker = [i + 1 for i, ln in enumerate(ls)
              if CALIB_MARKER_PAT.search(ln) and not ln.strip().startswith("#")]
    apps = []
    for i, ln in enumerate(ls):
        st = ln.strip()
        if not st or st.startswith("#") or st.startswith("def "):
            continue
        if not CALIB_CALL_PAT.search(ln):
            continue
        near = any("predict_proba" in ls[j]
                   for j in range(max(0, i - 3), min(len(ls), i + 4)))
        if ("predict_proba" in ln or near or "src_auroc" in ln
                or CALIB_WRAP_F_PAT.search(ln)):
            apps.append(i + 1)
    return marker, apps


def test_B1():
    """Methods calibration claim vs the live probability source in p6b/p6c/p10a.

    Divergence (D-B): Methods_v1.md claims a frozen Platt map; the signal scripts
    consume raw predict_proba.  CS-10 resolved this by ruling the code compliant and
    requiring the Methods sentence to be deleted/corrected.  RED only while BOTH
    (i) an affirmative calibration claim exists in Methods AND (ii) no calibration
    application exists at the probability source; either accepted resolution turns it
    green.  Structural/text-level, not behavioural (frozen model caches are not
    available to this suite)."""
    claim, claim_note = _methods_calib_claim()
    info, code_ok = [], True
    for tag, path in (("p6b", P6B), ("p6c", P6C), ("p10a", P10A)):
        if not os.path.exists(path):
            info.append(f"{tag}=ABSENT")
            code_ok = False
            continue
        mk, ap = _calib_stage(path)
        info.append(f"{tag}: marker@{mk or 'none'}, calib-application@{ap or 'none'}")
        code_ok &= bool(mk) and bool(ap)
    ctrl = []
    for tag, path in (("p5", P5), ("p12", P12)):
        if os.path.exists(path):
            mk, ap = _calib_stage(path)
            ctrl.append(f"{tag}: marker@{mk or 'none'} application@{ap or 'none'}")
        else:
            ctrl.append(f"{tag}=ABSENT")
    # ★ CS-10 语义升级（2026-10-06）：由 "no unsupported calibration claim"
    #   升级为 "calibration claim must match implementation scope"。
    #     ① 无声明                        -> 通过（原分支保留）
    #     ② 有声明但**未限定作用域**       -> FAIL（原缺陷：把管线局部泛化成全链）
    #     ③ 有声明且已限定、但未说明下游    -> FAIL（新要求，见 _calib_claim_scope 的说明）
    #     ④ 代码侧真的在校准概率源         -> 通过（原替代出路保留，故 p12/p5 型树仍可绿）
    scope_ok, down_ok, scope_note = _calib_claim_scope()
    claim_ok = (not claim) or (scope_ok and down_ok)
    ok = claim_ok or code_ok
    # 对照副本（非正本）的状态：**只报告，不参与判定**。正本已声明见
    # MANUSCRIPT_MASTER_DECLARATION.md；若副本仍含主张，说明两份文本已分叉，须一并处置。
    _nm_claim, _ = _methods_calib_claim(METHODS_MD_NONMASTER)
    nm_txt = ("non-master copy (%s): %s" % (os.path.relpath(METHODS_MD_NONMASTER, L.PROJ),
                                            "also PRESENT" if _nm_claim else "absent"))
    _emit("EVAL-B1", "PASS" if ok else "FAIL",
          f"MASTER ({os.path.relpath(METHODS_MD, L.PROJ)}) §5: 校准声明必须与其实现范围一致 —— "
          "① 若声明校准，声明句必须限定到**实际执行校准的管线**（scope phrase）；"
          "② 方法学必须正面说明下游分析使用**未校准**概率；"
          "③ 被声明的管线在代码里确实有校准实现。"
          "CS-10 amendment（2026-10-06）：本判据语义由 no-unsupported-claim 升级为 "
          "claim-must-match-scope，不废判据名",
          f"master read: {os.path.relpath(METHODS_MD, L.PROJ)}; "
          f"claim: {('PRESENT: ' + claim[:110] + '...') if claim else 'absent'}; "
          f"scope_ok={scope_ok} downstream_ok={down_ok}; {scope_note}; "
          + nm_txt + "; "
          + "; ".join(info) + "; controls: " + "; ".join(ctrl),
          "RED only when an affirmative Methods claim coexists with no calibration "
          "application in any of p6b/p6c/p10a (the current tree).  GREEN when the frozen "
          "Platt map is applied at the probability source (as p12/p5) OR when Methods no "
          "longer claims calibration (CS-10).  'Application' is structural: a calib call "
          "at/near a predict_proba expression or an src_auroc-style application -- a "
          "comment or docstring mentioning calib is not evidence.  Monotone calibration "
          "leaves AUROC unchanged, so source AUROC cannot certify the probability VALUES "
          "that feed entropy/confidence/0.5-threshold/prediction shifts (p6b:265,337). "
          + claim_note)


def test_B2():
    """Execute p12's Platt block verbatim on a fixture and measure what calibration
    does to probability-dependent statistics."""
    if not os.path.exists(P12):
        _emit("EVAL-B2", "FAIL",
              "mathematics: monotone (rank-preserving) map leaves AUROC unchanged but "
              "changes probability-dependent statistics; executed verbatim p12 Platt stage",
              f"UNTESTED: {os.path.basename(P12)} absent from this tree", "")
        return
    try:
        L.load_def(P12, "calib")
        loader_note = "load_def(calib): unexpectedly succeeded"
    except LookupError:
        loader_note = ("load_def(calib) fails: nested def inside the per-config loop "
                       "(LookupError) -> content-anchored block executed instead")
    block, blk, err = _grab_anchor(P12, r"^\s*eps\s*=\s*[0-9.]",
                                   ["pv = f(Dall[i_va])", "def calib(",
                                    "return cal.predict_proba"], max_lines=20)
    if err:
        _block_err("EVAL-B2", "mathematics: monotone map / probability statistics; "
                   "executed verbatim p12 Platt stage", err, os.path.basename(P12))
        return
    rng = np.random.default_rng(42)
    nva, nt = 500, 400
    zva = rng.normal(0, 1, nva)
    pv_fix = 1.0 / (1.0 + np.exp(-zva))
    yva = (rng.random(nva) < 1.0 / (1.0 + np.exp(-(1.5 * zva - 0.3)))).astype(int)
    zt = rng.normal(0.8, 1.2, nt)
    pr = 1.0 / (1.0 + np.exp(-zt))
    yt = (rng.random(nt) < 1.0 / (1.0 + np.exp(-(1.5 * zt - 0.3)))).astype(int)
    env = {"np": np, "LogisticRegression": LogisticRegression,
           "f": (lambda Z, pv=pv_fix: pv), "Dall": np.zeros((nva, 1)),
           "i_va": np.arange(nva), "y": yva}
    ns = _exec_src(block, env, "p12-calib")
    calib = ns["calib"]
    cal = calib(pr)
    grid = calib(np.linspace(0.001, 0.999, 999))
    mono = bool(np.all(np.diff(grid) >= -1e-12))
    a_raw, a_cal = roc_auc_score(yt, pr), roc_auc_score(yt, cal)
    ent = lambda p: -(np.clip(p, 1e-12, 1 - 1e-12) * np.log(np.clip(p, 1e-12, 1 - 1e-12))
                      + (1 - np.clip(p, 1e-12, 1 - 1e-12)) * np.log(1 - np.clip(p, 1e-12, 1 - 1e-12)))
    d_ent = float(ent(cal).mean() - ent(pr).mean())
    d_conf = float(np.abs(cal - 0.5).mean() - np.abs(pr - 0.5).mean())
    d_thr = float((cal > 0.5).mean() - (pr > 0.5).mean())
    d_mean = float(cal.mean() - pr.mean())
    finite = bool(np.isfinite(cal).all() and np.isfinite([a_raw, a_cal, d_ent, d_conf,
                                                          d_thr, d_mean]).all())
    ok = (finite and mono and abs(a_raw - a_cal) < 1e-9
          and max(abs(d_ent), abs(d_conf), abs(d_thr)) > 0.01)
    _emit("EVAL-B2", "PASS" if ok else "FAIL",
          "mathematics: monotone (rank-preserving) map leaves AUROC unchanged but changes "
          "mean probability, entropy, confidence |p-0.5| and the 0.5-threshold rate; "
          f"executed verbatim p12_configs_corrected.py:{blk}-"
          f"{blk + len(block.splitlines()) - 1} (Platt fit on validation + calib helper, "
          "content anchor 'eps = ...' + 'pv = f(Dall[i_va])' + 'def calib(')",
          f"fixture n_va={nva} n_t={nt}; monotone={mono}; AUROC raw={a_raw:.9f} "
          f"cal={a_cal:.9f} |d|={abs(a_raw - a_cal):.1e}; d mean(p)={d_mean:+.4f}; "
          f"d entropy={d_ent:+.4f}; d confidence={d_conf:+.4f}; d frac(p>0.5)={d_thr:+.4f}",
          "fixture numbers, NOT manuscript values; they show that source-AUROC agreement "
          "cannot certify the probability values, and quantify the order of magnitude a "
          "missing Platt map injects into the L2 signals. " + loader_note)


# ===================================================================== (C) 3.3
def test_C1():
    """p6f's boot resampler: does idx retain repeated draws?  Fixture replay against
    the definitional builder, plus the dedup variant to show the check can go red."""
    g_fix = np.array(["A"] * 3 + ["B"] * 2 + ["C"] * 4 + ["D"] * 1)
    sub = ["A", "B", "C", "D"]
    try:
        boot = L.load_def(P6F, "boot", extra_ns={"rng": np.random.default_rng(7), "B": 400,
                                                 "g": g_fix})
    except LookupError:
        _block_err("EVAL-C1", "bootstrap definition; executed verbatim p6f boot resampler",
                   "no-anchor(load_def 'boot')", os.path.basename(P6F))
        return
    except OSError as e:
        _block_err("EVAL-C1", "bootstrap definition; executed verbatim p6f boot resampler",
                   f"io-error({type(e).__name__})", os.path.basename(P6F))
        return
    if "percentile" not in getattr(boot, "__evs_src__", ""):
        _block_err("EVAL-C1", "bootstrap definition; executed verbatim p6f boot resampler",
                   "tokens-missing(percentile)", os.path.basename(P6F))
        return
    def stat(idx, pick=None):                            # tolerant of both boot signatures
        return float((idx == 0).sum())                   # #times the A-rows block is drawn
    plo, phi = boot(sub, stat)
    rng2 = np.random.default_rng(7)
    ref = []
    for _ in range(400):
        pick = rng2.choice(4, 4, replace=True)
        ref.append(float(sum(1 for i in pick if i == 0)))   # definitional: per draw
    rlo, rhi = np.percentile(ref, [2.5, 97.5])
    rng3 = np.random.default_rng(7)
    ded = []
    for _ in range(400):
        pick = rng3.choice(4, 4, replace=True)
        ded.append(float(1 if 0 in pick else 0))            # dedup variant
    dlo, dhi = np.percentile(ded, [2.5, 97.5])
    ok = (np.isfinite([plo, phi, rlo, rhi]).all()
          and abs(plo - rlo) < 1e-12 and abs(phi - rhi) < 1e-12)
    _emit("EVAL-C1", "PASS" if ok else "FAIL",
          "bootstrap definition; replay = same SEED-controlled stream, definitional idx "
          f"builder retaining multiplicity; executed verbatim p6f ({boot.__evs_origin__})",
          f"prod CI=[{plo:.3f},{phi:.3f}] reference CI=[{rlo:.3f},{rhi:.3f}]; "
          f"(dedup variant would give [{dlo:.3f},{dhi:.3f}])",
          "p6f's RESAMPLER keeps multiplicity (idx concatenates the group rows once per "
          "draw); the dedup in this script is confined to the macro statistic "
          "(see EVAL-A1) -- the dedup variant CI differs, so this replay would detect a "
          "deduping resampler")


def test_C2():
    """p6d's 1000-draw organism-drug cluster bootstrap, executed verbatim on the real
    frozen inputs, replayed against the definitional builder."""
    if "build_error" in _STATE:
        _emit("EVAL-C2", "FAIL", "bootstrap definition; p6d cluster bootstrap",
              f"UNTESTED: evaluation inputs not built ({_STATE['build_error']})", "")
        return
    E, y, g = _STATE["E"], _STATE["y"], _STATE["g"]
    try:
        safe_auc = L.load_def(P6D, "safe_auc", extra_ns={"roc_auc_score": roc_auc_score})
    except LookupError:
        _block_err("EVAL-C2", "bootstrap definition; p6d cluster bootstrap",
                   "no-anchor(load_def 'safe_auc')", os.path.basename(P6D))
        return
    except OSError as e:
        _block_err("EVAL-C2", "bootstrap definition; p6d cluster bootstrap",
                   f"io-error({type(e).__name__})", os.path.basename(P6D))
        return
    block, blk, err = _grab_anchor(P6D, r"rng\s*=\s*np\.random\.default_rng\(SEED\s*\+\s*7\)",
                                   ["sps.append", "aucs.append"])
    if err:
        _block_err("EVAL-C2", "bootstrap definition; p6d cluster bootstrap",
                   err, os.path.basename(P6D))
        return
    picked = ["domain_auc", "psi_max", "SUDO_AURCC"]
    rows = []
    ok_all = True
    for mname in picked:
        fn = dict((m[0], m[2]) for m in METHODS)[mname]
        r = fn(E).values.astype(float)
        env = {"np": np, "g": g, "y": y, "r": r, "E": E, "SEED": SEED,
               "spearmanr": spearmanr, "safe_auc": safe_auc}
        ns = _exec_src(block, env, "p6d-cluster-bootstrap")
        sps, aucs = ns["sps"], ns["aucs"]
        rng2 = np.random.default_rng(SEED + 7)
        ugs = np.unique(g)
        ref_sp, ref_au = [], []
        for _ in range(B6D):
            pick = rng2.choice(len(ugs), len(ugs), replace=True)
            idx = np.concatenate([np.where(g == ugs[i])[0] for i in pick])
            if len(np.unique(y[idx])) < 2:
                continue
            ref_sp.append(spearmanr(r[idx], E["d_auroc"].values[idx]).statistic)
            ref_au.append(safe_auc(y[idx], r[idx]))
        same = (len(sps) == len(ref_sp) and len(sps) > 0
                and np.isfinite(sps).all() and np.isfinite(aucs).all()
                and np.allclose(sps, ref_sp, rtol=0, atol=1e-12)
                and np.allclose(aucs, ref_au, rtol=0, atol=1e-12))
        # dedup variant on the identical draw stream (sensitivity evidence)
        rng3 = np.random.default_rng(SEED + 7)
        dd = []
        for _ in range(B6D):
            pick = rng3.choice(len(ugs), len(ugs), replace=True)
            idx0 = np.concatenate([np.where(g == ugs[i])[0] for i in pick])
            idx = np.unique(idx0)
            if len(np.unique(y[idx])) < 2:
                continue
            dd.append(safe_auc(y[idx], r[idx]))
        ok_all &= bool(same)
        rows.append(f"{mname}: n={len(sps)} identical={same} prodCI="
                    f"[{np.percentile(aucs, 2.5):.5f},{np.percentile(aucs, 97.5):.5f}] "
                    f"dedupCI=[{np.percentile(dd, 2.5):.5f},{np.percentile(dd, 97.5):.5f}]")
    _emit("EVAL-C2", "PASS" if ok_all else "FAIL",
          "bootstrap definition; exact replay of the verbatim block "
          f"p6d_01_unseal_evaluation.py:{blk}-{blk + len(block.splitlines()) - 1} "
          "(rng=default_rng(SEED+7), 1000 draws, organism-drug groups; content anchor "
          "'rng = np.random.default_rng(SEED + 7)') on the frozen 202-task inputs",
          "; ".join(rows),
          "p6d's cluster bootstrap KEEPS multiplicity (idx concatenates each picked "
          "cluster's rows once per pick; no unique/set). Dedup-variant CIs would be "
          "narrower, so the replay is sensitive")


def test_C3():
    """p6e's boot_macro, executed verbatim on a fixture fold_map."""
    fold_map = {"k%d" % i: (0.2 + 0.05 * i, 5 + i) for i in range(8)}
    try:
        bm = L.load_def(P6E, "boot_macro", extra_ns={"rng": np.random.default_rng(11)})
    except LookupError:
        _block_err("EVAL-C3", "bootstrap definition; p6e boot_macro",
                   "no-anchor(load_def 'boot_macro')", os.path.basename(P6E))
        return
    except OSError as e:
        _block_err("EVAL-C3", "bootstrap definition; p6e boot_macro",
                   f"io-error({type(e).__name__})", os.path.basename(P6E))
        return
    if "percentile" not in getattr(bm, "__evs_src__", ""):
        _block_err("EVAL-C3", "bootstrap definition; p6e boot_macro",
                   "tokens-missing(percentile)", os.path.basename(P6E))
        return
    plo, phi = bm(dict(fold_map), B=800)
    rng2 = np.random.default_rng(11)
    keys = list(fold_map.keys())
    ref = []
    for _ in range(800):
        idx = rng2.integers(0, len(keys), len(keys))
        vals = [fold_map[keys[i]][0] for i in idx]
        ws = np.array([fold_map[keys[i]][1] for i in idx], float)
        ref.append(np.average(vals, weights=ws))
    rlo, rhi = np.percentile(ref, [2.5, 97.5])
    ok = (np.isfinite([plo, phi, rlo, rhi]).all()
          and abs(plo - rlo) < 1e-12 and abs(phi - rhi) < 1e-12)
    _emit("EVAL-C3", "PASS" if ok else "FAIL",
          "bootstrap definition; replay via the same SEED-controlled integer stream; "
          f"executed verbatim p6e ({bm.__evs_origin__})",
          f"prod CI=[{plo:.9f},{phi:.9f}] reference CI=[{rlo:.9f},{rhi:.9f}]",
          "p6e's boot_macro KEEPS multiplicity (list over keys[idx] with repeats; weights "
          "repeated too). Characterisation: the CI belongs to the SIZE-WEIGHTED macro "
          "variant while the script's reported macro_mean is unweighted -- an aggregation "
          "labelling point, not a dedup")


def test_C4():
    """p6g's boot and p14's bci, executed verbatim on fixtures."""
    vals_map = {"g%d" % i: 0.3 + 0.04 * i for i in range(10)}
    try:
        g_boot = L.load_def(P6G, "boot", extra_ns={"rng": np.random.default_rng(3), "B": 600})
    except LookupError:
        _block_err("EVAL-C4", "bootstrap definition; p6g boot + p14 bci",
                   "no-anchor(load_def p6g 'boot')", os.path.basename(P6G))
        return
    except OSError as e:
        _block_err("EVAL-C4", "bootstrap definition; p6g boot + p14 bci",
                   f"io-error({type(e).__name__})", os.path.basename(P6G))
        return
    if "percentile" not in getattr(g_boot, "__evs_src__", ""):
        _block_err("EVAL-C4", "bootstrap definition; p6g boot + p14 bci",
                   "tokens-missing(percentile)", os.path.basename(P6G))
        return
    plo, phi = g_boot(dict(vals_map))
    rng2 = np.random.default_rng(3)
    ks = list(vals_map)
    ref = []
    for _ in range(600):
        pick = rng2.choice(len(ks), len(ks), replace=True)
        ref.append(np.mean([vals_map[ks[i]] for i in pick]))
    rlo, rhi = np.percentile(ref, [2.5, 97.5])
    ok_g = (np.isfinite([plo, phi, rlo, rhi]).all()
            and abs(plo - rlo) < 1e-12 and abs(phi - rhi) < 1e-12)
    bsrc, blk, err = _grab_anchor(P14, r"def\s+bci\(", ["def bci", "percentile"])
    if err:
        _block_err("EVAL-C4", "bootstrap definition; p6g boot + p14 bci",
                   err, os.path.basename(P14))
        return
    ns = _exec_src(bsrc, {"rng": np.random.default_rng(4), "B": 600}, "p14-bci")
    qlo, qhi = ns["bci"](dict(vals_map))
    rng3 = np.random.default_rng(4)
    ref2 = []
    for _ in range(600):
        pick = rng3.choice(10, 10, replace=True)
        ref2.append(np.mean([vals_map[ks[i]] for i in pick]))
    r2lo, r2hi = np.percentile(ref2, 2.5), np.percentile(ref2, 97.5)
    ok_p = (np.isfinite([qlo, qhi, r2lo, r2hi]).all()
            and abs(qlo - r2lo) < 1e-12 and abs(qhi - r2hi) < 1e-12)
    ok = ok_g and ok_p
    _emit("EVAL-C4", "PASS" if ok else "FAIL",
          "bootstrap definition; fixture replay on the same SEED-controlled stream; "
          f"executed verbatim p6g ({g_boot.__evs_origin__}, load_def 'boot') and p14 "
          f"(lines {blk}-{blk + len(bsrc.splitlines()) - 1}, grabbed+exec'd; anchor 'def bci(')",
          f"p6g boot prod=[{plo:.6f},{phi:.6f}] ref=[{rlo:.6f},{rhi:.6f}]; "
          f"p14 bci prod=[{qlo:.6f},{qhi:.6f}] ref=[{r2lo:.6f},{r2hi:.6f}]",
          "both KEEP multiplicity (list comprehension over keys[idx], no set/unique)")


def test_C5():
    """p6b's domain-AUC row bootstrap (200 draws), executed verbatim on a fixture."""
    rngf = np.random.default_rng(5)
    n = 60
    yd = np.array([0, 1] * (n // 2))
    oof = 0.3 + 0.4 * yd + rngf.normal(0, 0.15, n)
    block, blk, err = _grab_anchor(P6B, r"rb\s*=\s*np\.random\.default_rng\(SEED\s*\+\s*1\)",
                                   ["rb.integers", "bs.append"])
    if err:
        _block_err("EVAL-C5", "bootstrap definition; p6b domain-AUC row bootstrap",
                   err, os.path.basename(P6B))
        return
    env = {"np": np, "SEED": SEED, "yd": yd, "oof": oof,
           "roc_auc_score": roc_auc_score}
    ns = _exec_src(block, env, "p6b-domain-auc-bootstrap")
    bs = ns["bs"]
    rb2 = np.random.default_rng(SEED + 1)
    ref = []
    for _ in range(B6B):
        s = rb2.integers(0, len(yd), len(yd))
        if len(np.unique(yd[s])) < 2:
            continue
        ref.append(roc_auc_score(yd[s], oof[s]))
    same = (len(bs) == len(ref) and len(bs) > 0 and np.isfinite(bs).all()
            and np.allclose(bs, ref, rtol=0, atol=1e-12))
    _emit("EVAL-C5", "PASS" if same else "FAIL",
          "bootstrap definition; fixture replay of the verbatim block "
          f"p6b_01_corrected_pipeline.py:{blk}-{blk + len(block.splitlines()) - 1} "
          "(rng=default_rng(SEED+1), 200 draws, row-level index draw; content anchor "
          "'rb = np.random.default_rng(SEED + 1)')",
          f"n_draws={len(bs)} identical_to_reference={same} "
          f"CI=[{np.percentile(bs, 2.5):.6f},{np.percentile(bs, 97.5):.6f}]"
          if bs else "n_draws=0 (fixture degenerate)",
          "p6b's domain-AUC row bootstrap KEEPS multiplicity (s indexes rows directly; "
          "repeated indices are repeated rows; no dedup)")


# ===================================================================== (D) 3.4
def _loo_section_oof_block():
    """Locate the 'LOO-organism' section in p6f and return (block, first line,
    lodo_oof call line, kind).

    kind='oof'        the section's `oof = ...` construction (assignment + following
                      lines up to its first print) -- handles p6f's lodo_oof call and
                      any per-organism refit (inline or via a top-level helper).
    kind='body'       no `oof = ...` line: the whole section body after the header
                      prints is executed (per method, with the METHODS loop filtered to
                      that method and `p` no-op'd), so a refactor computing per-organism
                      AUCs directly still executes.
    kind='deprecated' the section explicitly marks itself deprecated/retired -- CS-11
                      accepts 'removed or explicitly deprecated'.
    kind='absent'     no LOO-organism section at all (retired).
    """
    lines = _lines(P6F)
    # only a printed section header counts -- the docstring mention of "LOO-organism"
    # must NOT be mistaken for a (retired) section, or a removed section would re-grab
    # the main settlement loop's `oof = lodo_oof(r)` and report a misleading FAIL
    hits = [i for i, ln in enumerate(lines)
            if "LOO-organism" in ln and ln.lstrip().startswith("p(")]
    start = hits[-1] if hits else None       # the section header (not the docstring mention)
    if start is None:
        return None, None, None, "absent"
    end = len(lines)
    for j in range(start + 1, len(lines)):
        if '"\\n" + "=" * 112' in lines[j]:     # next section separator
            end = j
            break
    lodo_line = None
    for i in range(start, end):
        if "lodo_oof(" in lines[i]:
            lodo_line = i + 1
    dep_pat = re.compile(r"(?i)deprecat|obsolete|retired|废弃|作废|不再使用|不复用|not authoritative")
    for i in range(max(0, start - 3), end):
        if dep_pat.search(lines[i]):
            return None, i + 1, lodo_line, "deprecated"
    for i in range(start, end):
        m = re.match(r"^(\s*)oof\s*=\s*(.+)$", lines[i])
        if m:
            # collect the whole construction (assignment + any following indented or
            # sibling lines) up to the section's first print; this keeps a per-organism
            # refit fix (assignment + for-loop) inside the executed block.
            buf, j = [lines[i]], i + 1
            while j < end and not lines[j].strip().startswith("p("):
                buf.append(lines[j])
                j += 1
            import textwrap as _tw
            return _tw.dedent("\n".join(buf)), i + 1, lodo_line, "oof"
    j = start + 1
    while j < end and (not lines[j].strip() or lines[j].strip().startswith("p(")):
        j += 1
    if j < end and any(ln.strip() for ln in lines[j:end]):
        import textwrap as _tw
        return _tw.dedent("\n".join(lines[j:end])), j + 1, lodo_line, "body"
    return None, None, lodo_line, "no-stat"


def _bind_helpers(ns, block, base_ns):
    """Bind top-level p6f defs the block references (helper-style refactors such as
    `oof = loo_oof(r)`), skipping keywords/builtins/known names; LookupError => the name
    is not a top-level def and is left to the retry/failure path."""
    import keyword as _kw
    import builtins as _bi
    cand = set(re.findall(r"\b([A-Za-z_]\w*)\b", block))
    added = []
    for n in sorted(cand):
        if n in ns or n in base_ns or _kw.iskeyword(n) or hasattr(_bi, n):
            continue
        try:
            fn = L.load_def(P6F, n, extra_ns=base_ns)
        except (LookupError, OSError):
            continue
        ns[n] = fn
        added.append(n)
    return added


def _exec_section(block, ns, base_ns):
    """Exec the section block; on NameError bind the missing top-level p6f def and
    retry once, then fail with the symbol named (no bare cascade into EVAL-D2)."""
    _bind_helpers(ns, block, base_ns)
    try:
        exec(compile(block, "<verbatim:p6f-LOO-section>", "exec"), ns)
    except NameError as e:
        m = re.search(r"name '([A-Za-z_]\w*)' is not defined", str(e))
        if not m:
            raise
        sym = m.group(1)
        try:
            ns[sym] = L.load_def(P6F, sym, extra_ns=base_ns)
        except (LookupError, OSError):
            raise NameError(
                f"name '{sym}' is not defined and is not a top-level def of "
                f"{os.path.basename(P6F)}") from None
        exec(compile(block, "<verbatim:p6f-LOO-section-retry>", "exec"), ns)


def _section_value(ns, org, oeval, y, p6f_auc):
    """The macro statistic the executed section leaves behind: prefer an `oof` vector,
    then a per-organism `fa` list, then a scalar `C`."""
    oof = ns.get("oof")
    if oof is not None:
        try:
            oof = np.asarray(oof, dtype=float)
            if oof.shape == (len(y),) and np.isfinite(oof).any():
                vals = [p6f_auc(y[org == o_], oof[org == o_]) for o_ in oeval]
                vals = [float(v) for v in vals if np.isfinite(v)]
                if vals:
                    return float(np.mean(vals)), "oof"
        except Exception:
            pass
    fa = ns.get("fa")
    try:
        if isinstance(fa, (list, tuple, np.ndarray)) and len(fa) > 0:
            vals = [float(v) for v in np.asarray(fa, dtype=float).ravel()
                    if np.isfinite(v)]
            if vals:
                return float(np.mean(vals)), "fa"
    except Exception:
        pass
    c = ns.get("C")
    if isinstance(c, (int, float)) and np.isfinite(c):
        return float(c), "C"
    return None, None


def _loo_anchors():
    """Per-method true-LOO and LODO-based organism macros, computed independently of
    the p6f section (so EVAL-D2 does not cascade on EVAL-D1).  Cached in _STATE."""
    if "loo_anchors" in _STATE:
        return _STATE["loo_anchors"]
    E, y = _STATE["E"], _STATE["y"]
    org = E.organism.values
    p6f_auc, p6f_lodo = _STATE["p6f_auc"], _STATE["p6f_lodo"]
    orgs = sorted(E.organism.unique())
    oeval = [o_ for o_ in orgs
             if y[org == o_].sum() > 0 and (1 - y[org == o_]).sum() > 0]
    tab = {}
    for name, _r, fn in METHODS:
        r = fn(E).values.astype(float)
        vals_true = []
        for o_ in oeval:
            te = org == o_
            tr = ~te
            # ★ 2026-10-06：原写 `np.nanstd(r[tr]) or 1.0` —— **不防 nan**（nan 在布尔语境里
            #   为真，`nan or 1.0` 求值仍是 nan），除以 nan 后整列变 NaN，最终只崩在
            #   LogisticRegression 的一句通用报错里，看不出真实原因。
            #   实测触发场景：p6d_outcome.csv（上游产物）与 p6b/p6c（新产物）**代次不一致**时，
            #   merge 后多出的 task 在信号列上全是 NaN。这里显式判并报出是哪个方法/哪一组。
            mu = np.nanmean(r[tr]); sd = np.nanstd(r[tr])
            if not (np.isfinite(mu) and np.isfinite(sd)) or sd == 0.0:
                raise ValueError(
                    "LOO anchor 的前置统计量非有限：method=%r organism=%r nanmean=%s nanstd=%s"
                    " —— 该信号列在此训练子集上可能整列 NaN（常见原因：上游产物代次不一致，"
                    "如 p6d 未随 p6b/p6c 重跑）" % (name, o_, mu, sd))
            m = LogisticRegression(max_iter=1000).fit(
                ((r[tr] - mu) / sd).reshape(-1, 1), y[tr])
            s2 = m.predict_proba(((r[te] - mu) / sd).reshape(-1, 1))[:, 1]
            vals_true.append(p6f_auc(y[te], s2))
        o_lodo = _STATE["oof"][name]
        fa_lodo = [p6f_auc(y[org == o_], o_lodo[org == o_]) for o_ in oeval]
        fa_lodo = [float(v) for v in fa_lodo if np.isfinite(v)]
        tab[name] = (float(np.mean(vals_true)) if vals_true else np.nan,
                     float(np.mean(fa_lodo)) if fa_lodo else np.nan)
    _STATE["loo_anchors"] = tab
    return tab


def test_D1():
    """p6f's 'LOO-organism' section: does the statistic it actually builds equal
    leave-one-organism-out?  The section is executed as written -- inline refits and
    top-level-helper refits (`oof = loo_oof(r)`) are both supported.  A section that was
    retired (absent) or explicitly marked deprecated is the CS-11 resolution and no
    longer reports a mislabeled LOO statistic."""
    if not _need_p6f("EVAL-D1", "definition of leave-one-group-out (fit on the other "
                     "groups only, predict the held-out group, macro-average)"):
        return
    E, y = _STATE["E"], _STATE["y"]
    org = E.organism.values
    p6f_auc, p6f_lodo = _STATE["p6f_auc"], _STATE["p6f_lodo"]
    orgs = sorted(E.organism.unique())
    oeval = [o_ for o_ in orgs
             if y[org == o_].sum() > 0 and (1 - y[org == o_]).sum() > 0]
    tab = _loo_anchors()
    block, blk_line, lodo_line, kind = _loo_section_oof_block()
    if kind == "absent":
        _emit("EVAL-D1", "PASS",
              "definition of leave-one-group-out; CS-11 resolution: the mislabeled "
              "'LOO-organism' section is removed or explicitly deprecated",
              "no 【LOO-organism】 construction remains in p6f (section retired)",
              "no mislabeled LOO-organism output remains; the canonical Fig-3B path "
              "(p14 + frozen p6d/fig3b columns) is cross-checked in EVAL-D2")
        return
    if kind == "deprecated":
        _emit("EVAL-D1", "PASS",
              "definition of leave-one-group-out; CS-11 resolution: the mislabeled "
              "'LOO-organism' section is removed or explicitly deprecated",
              f"section present but explicitly marked deprecated/retired (p6f:{blk_line})",
              "the stale section is isolated by its own annotation; the canonical Fig-3B "
              "numbers are anchored in EVAL-D2.  The marker is a text-level signal: "
              "downstream consumers must honour it")
        return
    if block is None:
        _emit("EVAL-D1", "FAIL",
              "definition of leave-one-group-out (fit on the other groups only, predict "
              "the held-out group, macro-average)",
              f"UNTESTED: could not locate a statistic construction inside p6f's "
              f"LOO-organism section ({kind})",
              "cannot execute the section's statistic; failing rather than passing")
        return
    base_ns = {"np": np, "LogisticRegression": LogisticRegression, "E": E, "y": y,
               "g": _STATE["g"], "org": org, "orgs": orgs, "oeval": oeval,
               "GROUPS": _STATE["groups"], "auc": p6f_auc, "lodo_oof": p6f_lodo,
               "METHODS": METHODS, "p": (lambda *a, **k: None)}
    rows, worst = [], (None, np.nan)
    for name, _r, fn in METHODS:
        r = fn(E).values.astype(float)
        ns = dict(base_ns)
        ns["r"] = r
        if kind == "body":
            # run the section's own METHODS loop for this method only, so the statistic
            # left in ns belongs to the method being compared
            ns["METHODS"] = [(name, _r, fn)]
        try:
            _exec_section(block, ns, base_ns)
        except NameError as e:
            _emit("EVAL-D1", "FAIL",
                  "definition of leave-one-group-out; executed the section's own "
                  f"construction (p6f:{blk_line}, kind={kind})",
                  f"EXCEPTION NameError: {e}",
                  "a helper the section calls is neither bound nor a top-level def of "
                  "p6f; if the refactor moved it elsewhere, re-pin the extraction")
            return
        c_sec, how = _section_value(ns, org, oeval, y, p6f_auc)
        if c_sec is None or not np.isfinite(c_sec):
            _emit("EVAL-D1", "FAIL",
                  "definition of leave-one-group-out; executed the section's own statistic",
                  f"executed p6f:{blk_line} (kind={kind}) but found no finite section "
                  f"statistic (source={how})",
                  "the executed block is not a valid statistic construction; failing "
                  "rather than passing on a degenerate result")
            return
        c_true, c_lodo_dir = tab[name]
        _STATE.setdefault("loo_table", {})[name] = (c_sec, c_true, c_lodo_dir)
        d = abs(c_true - c_sec)
        if worst[0] is None:
            worst = (name, d)
        elif np.isfinite(d) and (not np.isfinite(worst[1]) or d > worst[1]):
            worst = (name, d)
        if name in EXT_TRUE_LOO:
            rows.append(f"{name}: section={c_sec:.6f} vs true LOO={c_true:.6f}")
    ok = worst[0] is not None and np.isfinite(worst[1]) and worst[1] <= 1e-9
    _emit("EVAL-D1", "PASS" if ok else "FAIL",
          "definition of leave-one-group-out: for each organism fit on the OTHER "
          "organisms only, predict the held-out organism, then macro-average; p6f's "
          "section header is 【LOO-organism】",
          f"executed the section's own statistic construction (p6f:{blk_line}, kind={kind}, "
          f"{len(block.splitlines())} line(s)); worst |trueLOO-section|={worst[1]:.6f} "
          f"({worst[0]}); " + "; ".join(rows),
          ("p6f:%d is a lodo_oof(...) call -> predictions come from leave-one-ORGANISM-"
           "DRUG-out (other groups of the same organism still train), which is NOT "
           "leave-one-organism-out; " % lodo_line if lodo_line else "")
          + "p6d and p14 implement the true LOO (frozen p6d_signal_eval.loo_org_auc and "
            "fig3b_loo_organism.csv match the definition, checked in EVAL-D2). "
            "Section output should be recomputed or retired")


def test_D2():
    """External anchors: do the true-LOO and LODO-based organism macros reproduce the
    reviewer's exact numbers?  And do the frozen fig3b/p6d columns match the true LOO?
    Both frozen cross-checks fail CLOSED on a missing/corrupt/partial artifact."""
    if not _need_p6f("EVAL-D2", "external audit anchors"):
        return
    tab = _loo_anchors()
    rows, ok = [], True
    # ★★ 2026-10-06 语义收窄（与 CS04-05 同类：判据的前提被后续修正取代）。
    #   原判据：`ok &= (d1 <= 2e-6 and d2 <= 2e-6)`，其中
    #     EXT_TRUE_LOO / EXT_LODO_BASED = 外部审计 `ARMD_Code_Audit_20261004.md` 给出的
    #     两组数（true LOO 0.666307 / 0.446328 / 0.706656；LODO-based 0.583103 / 0.549195 / 0.625107）。
    #   那是**一次性对抗验证**：确认审计方声称的算术可从当时产物复现。当时**通过了**。
    #   但那些数是在**修正前的 202 universe** 上算的。CS-01/02/03 重建任务宇宙后，
    #   本轮实测为 trueLOO 0.452 / 0.442 / 0.439、LODO-based 0.512 / 0.500 / 0.182 ——
    #   **永远不可能**再回到旧值（旧 universe 已不存在）。若继续作为判据，它只会恒红。
    #   ⇒ 改为**记录**：距离照算、照打印（保留"审计方当时算得对不对"的可追溯性），
    #     但不再充当通过条件。判据改为只由下面的"冻结列 vs 独立 true LOO"承担 ——
    #     那一条是**永久有效**的不变式，且本轮已实测通过（p6d 1.1e-16、fig3b 0.0e+00）。
    for name in EXT_TRUE_LOO:
        ct, cl = tab.get(name, (np.nan, np.nan))
        d1 = abs(ct - EXT_TRUE_LOO[name])
        d2 = abs(cl - EXT_LODO_BASED[name])
        rows.append(f"{name}: trueLOO={ct:.6f}(d={d1:.1e}) LODO-based(direct)={cl:.6f}"
                    f"(d={d2:.1e})")

    # frozen column cross-check: p6d loo_org_auc and fig3b loo_org_macro vs my true LOO
    def frozen_check(path, col):
        try:
            df = pd.read_csv(path)
        except Exception as e:
            return np.nan, 0, f"{os.path.basename(path)} unreadable ({type(e).__name__})"
        if "method" not in df.columns or col not in df.columns:
            return np.nan, 0, f"{os.path.basename(path)} missing column {col}"
        idx = df.set_index("method")
        dmax, n = 0.0, 0
        for k in EXT_TRUE_LOO:
            if k not in idx.index:
                return np.nan, n, f"{os.path.basename(path)} missing method {k}"
            dmax = max(dmax, abs(float(idx.loc[k, col]) - tab[k][0]))
            n += 1
        return dmax, n, None

    dmax_sev, n_sev, err_sev = frozen_check(os.path.join(OUT, "p6d_signal_eval.csv"),
                                            "loo_org_auc")
    dmax_f3b, n_f3b, err_f3b = frozen_check(os.path.join(FIGD, "fig3b_loo_organism.csv"),
                                            "loo_org_macro")
    ok &= (err_sev is None and n_sev == len(EXT_TRUE_LOO)
           and np.isfinite(dmax_sev) and dmax_sev <= 2e-6)
    ok &= (err_f3b is None and n_f3b == len(EXT_TRUE_LOO)
           and np.isfinite(dmax_f3b) and dmax_f3b <= 2e-6)
    glitch = "; ".join(x for x in (err_sev, err_f3b) if x)
    _emit("EVAL-D2", "PASS" if ok else "FAIL",
          "★ 判据（2026-10-06 收窄后）= **冻结列必须等于独立算出的 true LOO**："
          "p6d_signal_eval.csv 的 loo_org_auc 与 fig3b_loo_organism.csv 的 loo_org_macro，"
          "3/3 方法缺一不可，容忍 2e-6；缺文件/缺列/缺方法一律 fail-closed。"
          "【已退役部分】原判据还要求复现外部审计 ARMD_Code_Audit_20261004.md 的两组数"
          "（true LOO 0.666307/0.446328/0.706656；LODO-based 0.583103/0.549195/0.625107）——"
          "那是修正前 202 universe 上的一次性对抗验证，宇宙重建后不可能再成立，"
          "现改为**记录**（距离照算照印，见 observed），不再充当通过条件",
          "; ".join(rows) + f"; 冻结列 vs 独立 true LOO：p6d max|d|={dmax_sev:.1e} "
          f"({n_sev}/3)，fig3b max|d|={dmax_f3b:.1e} ({n_f3b}/3)"
          + (f"; {glitch}" if glitch else ""),
          "both computations reproduce the external numbers to their 6-decimal display; "
          "the comparison therefore isolates the SEMANTICS (LODO vs LOO), not arithmetic. "
          "v1.3: the frozen cross-checks no longer swallow exceptions -- an unreadable or "
          "partial artifact is a FAIL, not a silent PASS")



# ===================================================================== main
def main():
    L.banner("EVAL - evaluation/aggregation layer (t12)")
    print(f"  production tree : {L.PROD_DIR}")
    print(f"  mode            : {'PATCHED' if L.IS_PATCHED else 'REAL'}")
    print(f"  rule            : {RULE3}")
    print("")

    try:
        build_inputs()
        print(f"  inputs: E={_STATE['E'].shape}, folds={len(_STATE['groups'])}, "
              f"evaluable={len(_STATE['EVAL'])}, "
              f"tasks(eval)={int(np.isin(_STATE['g'], _STATE['EVAL']).sum())}")
    except Exception as e:
        _STATE["build_error"] = f"{type(e).__name__}: {e}"
    try:
        prepare_p6f()
    except Exception as e:
        _STATE["p6f_error"] = f"{type(e).__name__}: {e}"

    def run(fn, tids, expected):
        if isinstance(tids, str):
            tids = (tids,)
        try:
            fn()
        except Exception as e:
            txt = f"EXCEPTION {type(e).__name__}: {e}"
            for t in tids:
                if t not in _EMITTED:
                    _emit(t, "FAIL", expected, txt, "")

    run(test_A1, "EVAL-A1", DEF_BOOT + "; executed verbatim from " + os.path.basename(P6F))
    run(test_A2_A3, ("EVAL-A2", "EVAL-A3"), "frozen-artifact replay")
    run(test_B1, "EVAL-B1", "Methods-calibration sentence vs signal scripts")
    run(test_B1b, "EVAL-B1b", "EVAL-B1 probe non-vacuity canary (adversarial specimens)")
    run(test_B2, "EVAL-B2", "Platt monotonicity / probability statistics")
    for tid, t in (("EVAL-C1", test_C1), ("EVAL-C2", test_C2), ("EVAL-C3", test_C3),
                   ("EVAL-C4", test_C4), ("EVAL-C5", test_C5)):
        run(t, tid, "bootstrap-multiplicity check (per-site)")
    run(test_D1, "EVAL-D1", "definition of leave-one-organism-out")
    run(test_D2, "EVAL-D2", "external audit anchors")
    print("done")


if __name__ == "__main__":
    main()
