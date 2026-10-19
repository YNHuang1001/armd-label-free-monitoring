# -*- coding: utf-8 -*-
r"""LODO fold composition QC —— 不改任何结果，只查清 aggregation ambiguity

回答三个问题（GPT6 指定）：
 (1) 全部 organism-drug fold 中，多少个同时含 harmful 与 non-harmful（AUROC 可定义）？[原写死 "38 个"；fold 数随宇宙变化，本期 37，由本脚本实测报告]
 (2) 每个可定义 fold 的 n+ / n- 是多少？
 (3) per-fold AUROC 的分布（并给 cluster bootstrap CI on macro mean）

同时并列报告 pooled-OOF AUROC 与 evaluable-fold macro mean AUROC。
不做任何 retroactive primary 指定。

输出：report_p6e_lodo_qc.txt / p6e_lodo_fold_composition.csv / p6e_lodo_aggregation.csv
"""
import os, warnings
import numpy as np, pandas as pd
warnings.filterwarnings("ignore")
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score

BASE = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))); OUT = os.path.dirname(os.path.abspath(__file__))
SEED = 20261005; rng = np.random.default_rng(SEED)

rep = open(os.path.join(OUT, "report_p6e_lodo_qc.txt"), "w", encoding="utf-8")
def p(*a):
    line = " ".join(str(x) for x in a); print(line, flush=True); rep.write(line + "\n"); rep.flush()

O = pd.read_csv(os.path.join(OUT, "p6d_outcome.csv"))
SIG = pd.read_csv(os.path.join(OUT, "p6b_label_free_signals_blinded.csv"))
L3 = pd.read_csv(os.path.join(OUT, "p6c_L3_blinded.csv"))

E = O.merge(SIG[["task_id", "domain_auc", "mmd", "psi_mean", "psi_max", "js_mean", "js_max",
                 "pred_mean_shift", "entropy_mean_shift", "entropy_js", "conf_dist_js",
                 "pred_pos_prev_shift"]], on="task_id", how="left") \
     .merge(L3[["task_id", "SUDO_AURCC", "AUROC_IOW", "AUROC_OM", "AUROC_DR", "w_max", "w_ESS"]],
            on="task_id", how="left")
E["group"] = E["organism"] + "|" + E["drug"]

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

y_all = E["harmful"].values.astype(int)
g_all = E["group"].values

# ============ (1)(2) fold composition ============
rows = []
for g in sorted(np.unique(g_all)):
    m = g_all == g
    n = int(m.sum()); npos = int(y_all[m].sum()); nneg = n - npos
    org, drug = g.split("|")
    rows.append(dict(group=g, organism=org, drug=drug, n=n, n_harmful=npos,
                     n_nonharmful=nneg, auroc_definable=int(npos > 0 and nneg > 0)))
FC = pd.DataFrame(rows).sort_values(["organism", "drug"]).reset_index(drop=True)
FC.to_csv(os.path.join(OUT, "p6e_lodo_fold_composition.csv"), index=False)

n_fold = len(FC); n_def = int(FC.auroc_definable.sum())
p("=" * 100)
p("[QC-1/2] LODO fold composition  (organism-drug group = fold)")
p("=" * 100)
p(f"  总 fold 数                       = {n_fold}")
p(f"  AUROC 可定义 fold（n+>0 且 n->0）= {n_def}  ({n_def/n_fold*100:.1f}%)")
p(f"  不可定义 fold                    = {n_fold - n_def}"
  f"  [全 harmful = {int((FC.n_nonharmful == 0).sum())}, 全 non-harmful = {int((FC.n_harmful == 0).sum())}]")
p(f"  fold 规模 n:  min={FC.n.min()}  P25={FC.n.quantile(.25):.0f}  median={FC.n.median():.0f}"
  f"  P75={FC.n.quantile(.75):.0f}  max={FC.n.max()}")
p(f"  可定义 fold 的 n+ : min={FC[FC.auroc_definable==1].n_harmful.min()}"
  f"  median={FC[FC.auroc_definable==1].n_harmful.median():.0f}"
  f"  max={FC[FC.auroc_definable==1].n_harmful.max()}")
p(f"  可定义 fold 的 n- : min={FC[FC.auroc_definable==1].n_nonharmful.min()}"
  f"  median={FC[FC.auroc_definable==1].n_nonharmful.median():.0f}"
  f"  max={FC[FC.auroc_definable==1].n_nonharmful.max()}")
p("")
p("  全部 fold 明细 (group, n, n+, n-, 可定义):")
for _, r in FC.iterrows():
    flag = "OK " if r.auroc_definable else "NA "
    p(f"    {flag}{r.group:<14} n={r.n:>2}  n+={r.n_harmful:>2}  n-={r.n_nonharmful:>2}")

# ============ (3) aggregation ============
def fit_predict(tr, te, r):
    mu, sd = np.nanmean(r[tr]), np.nanstd(r[tr]) or 1.0
    m = LogisticRegression(max_iter=1000).fit(((r[tr] - mu) / sd).reshape(-1, 1), y_all[tr])
    return m.predict_proba(((r[te] - mu) / sd).reshape(-1, 1))[:, 1]

def safe_auc(y, s):
    m = np.isfinite(s)
    return roc_auc_score(y[m], s[m]) if len(np.unique(y[m])) > 1 else np.nan

def run(r):
    oof = np.full(len(r), np.nan); fold = {}
    for g in np.unique(g_all):
        te = g_all == g; tr = ~te
        if len(np.unique(y_all[tr])) < 2 or te.sum() == 0: continue
        oof[te] = fit_predict(tr, te, r)
        a = safe_auc(y_all[te], oof[te])
        if np.isfinite(a): fold[g] = (a, int(te.sum()))
    pooled = safe_auc(y_all, oof)
    fa = np.array([v[0] for v in fold.values()])
    return pooled, fold, fa

def boot_macro(fold_map, B=2000):
    """cluster bootstrap over evaluable folds, resample folds with replacement"""
    keys = list(fold_map.keys()); out = []
    for _ in range(B):
        idx = rng.integers(0, len(keys), len(keys))
        vals = [fold_map[keys[i]][0] for i in idx]
        # weight by fold size for a size-weighted macro variant
        ws = np.array([fold_map[keys[i]][1] for i in idx], float)
        out.append(np.average(vals, weights=ws))
    return np.nanpercentile(out, [2.5, 97.5])

res = []; foldtab = []
for name, role, fn in METHODS:
    r = fn(E).values.astype(float)
    pooled, fold_map, fa = run(r)
    macro = np.mean(fa) if len(fa) else np.nan
    wmacro = np.average(fa, weights=[fold_map[k][1] for k in fold_map]) if len(fa) else np.nan
    lo, hi = boot_macro(fold_map) if len(fa) else (np.nan, np.nan)
    res.append(dict(method=name, role=role, n_evaluable_folds=len(fold_map),
                    pooled_OOF=pooled, macro_mean=macro, size_weighted_macro=wmacro,
                    macro_lo=lo, macro_hi=hi,
                    fold_min=np.min(fa) if len(fa) else np.nan,
                    fold_p25=np.percentile(fa, 25) if len(fa) else np.nan,
                    fold_median=np.median(fa) if len(fa) else np.nan,
                    fold_p75=np.percentile(fa, 75) if len(fa) else np.nan,
                    fold_max=np.max(fa) if len(fa) else np.nan,
                    n_fold_above_0p5=int((fa > 0.5).sum()) if len(fa) else 0))
    for k, (a, nn) in fold_map.items():
        foldtab.append(dict(method=name, group=k, fold_auc=a, fold_n=nn))

R = pd.DataFrame(res); R.to_csv(os.path.join(OUT, "p6e_lodo_aggregation.csv"), index=False)
pd.DataFrame(foldtab).to_csv(os.path.join(OUT, "p6e_lodo_fold_auc_long.csv"), index=False)

p("")
p("=" * 100)
p("[QC-3] 两种聚合口径并列（不做 retroactive primary 指定）")
p("=" * 100)
p(f"{'method':<22}{'role':<12}{'kFold':>6}{'pooledOOF':>11}{'macroMean':>11}"
  f"{'sizeWtd':>9}{'macro_95%CI':>18}{'F>0.5':>7}")
for _, r in R.iterrows():
    ci = f"[{r.macro_lo:.3f},{r.macro_hi:.3f}]"
    p(f"{r.method:<22}{r.role:<12}{r.n_evaluable_folds:>6}{r.pooled_OOF:>11.3f}"
      f"{r.macro_mean:>11.3f}{r.size_weighted_macro:>9.3f}{ci:>18}{r.n_fold_above_0p5:>7}")

p("")
p("  逐 fold AUROC 分布（仅可定义 fold）")
p(f"{'method':<22}{'min':>7}{'P25':>7}{'median':>8}{'P75':>7}{'max':>7}{'IQR':>7}")
for _, r in R.iterrows():
    p(f"{r.method:<22}{r.fold_min:>7.3f}{r.fold_p25:>7.3f}{r.fold_median:>8.3f}"
      f"{r.fold_p75:>7.3f}{r.fold_max:>7.3f}{r.fold_p75-r.fold_p25:>7.3f}")

p("")
p("  最大分歧项 (pooledOOF - macroMean):")
R["gap"] = R.pooled_OOF - R.macro_mean
for _, r in R.reindex(R.gap.abs().sort_values(ascending=False).index).head(6).iterrows():
    p(f"    {r.method:<22} pooled={r.pooled_OOF:.3f}  macro={r.macro_mean:.3f}  gap={r.gap:+.3f}")

p("")
p("  NOTE: aggregation strategy was NOT prespecified in Protocol v1.3; both are reported.")
rep.close()
print("done")
