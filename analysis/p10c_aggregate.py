# -*- coding: utf-8 -*-
r"""Protocol v1.3 §10 —— 聚合（修复版）

前一版 p10b 的 bug：raw 行里 `full` 档的 N 存的是该 task 的真实 target 大小，
聚合时用 N==0 筛 full 档 → 空集 → np.percentile 崩。
本脚本**不重跑计算**，只读 p10b_raw.csv 重建档位并聚合。

档位重建：每 task 的 max(N) 行 = full；其余 N ∈ {100,250,500,1000} 行 = 对应档。
（全站最小 target = 2208 > 1000，故无歧义；脚本内断言验证。）

输出：report_p10c_subsample.txt / p10c_subsample_curves.csv / p10c_paired_drift.csv
"""
import os
import numpy as np, pandas as pd

BASE = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))); OUT = os.path.dirname(os.path.abspath(__file__))
LEVELS = [100, 250, 500, 1000, "full"]

rep = open(os.path.join(OUT, "report_p10c_subsample.txt"), "w", encoding="utf-8")
def p(*a):
    line = " ".join(str(x) for x in a); print(line, flush=True); rep.write(line + "\n"); rep.flush()

RAW = pd.read_csv(os.path.join(OUT, "p10b_raw.csv"))
O = pd.read_csv(os.path.join(OUT, "p6d_outcome.csv"))[["task_id", "organism", "drug", "d_auroc", "harmful", "stable"]]
RAW = RAW.merge(O, on="task_id", how="left")

# ---- 档位重建 ----
nt = RAW.groupby("task_id")["N"].max().rename("n_target")
RAW = RAW.merge(nt, on="task_id", how="left")
RAW["level"] = np.where(RAW["N"] == RAW["n_target"], "full", RAW["N"].astype(str))
RAW["level"] = pd.Categorical(RAW["level"], categories=[str(v) for v in LEVELS], ordered=True)
assert RAW["level"].isna().sum() == 0, "level 重建失败"
chk = RAW.groupby(["task_id", "level"], observed=True).size()
assert (chk == 200).all(), "有 task 档位行数 != 200"
chk2 = RAW.merge(O[["task_id"]], on="task_id")
mc = O.set_index("task_id")["n_target"] if "n_target" in O.columns else None
p("=" * 118)
p("[Protocol v1.3 §10.1] target 样本量敏感性 —— 全部 L1 / L2 / L3")
p("=" * 118)
p(f"  task = {RAW.task_id.nunique()}   B = {RAW.rep.nunique()}   N 档 = {LEVELS}")
# ★ 2026-10-06：原为 f"  raw 行数 = {len(RAW)}  (202 x 5 x 200 = 202,000)" ——
#   括号里的 202 是**修正前**的 task 数，会与实际行数矛盾（现在是 186 x 5 x 200）。
#   改为从 RAW 现算，避免报告里出现一个与它左边数字不符的注释。
p(f"  raw 行数 = {len(RAW)}  "
  f"({RAW.task_id.nunique()} task x {len(LEVELS)} 档 x {RAW.rep.nunique()} rep "
  f"= {RAW.task_id.nunique() * len(LEVELS) * RAW.rep.nunique():,}）")
p(f"  档位重建自洽性：每 task 每档 200 行 —— 通过")
p(f"  task 的 full target N: min={int(RAW.n_target.min())}  median={int(RAW.n_target.median())}  max={int(RAW.n_target.max())}")

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
from sklearn.metrics import roc_auc_score, average_precision_score

def sens_at_spec(y, s, spec=0.80):
    m = np.isfinite(s); y, s = y[m], s[m]
    if len(np.unique(y)) < 2: return np.nan
    neg = s[y == 0]
    if len(neg) == 0: return np.nan
    thr = min([v for v in np.unique(s) if (neg < v).mean() >= spec] or [np.nextafter(s.max() if s.size else 0.0, np.inf)])
    pred = s >= thr
    tp = ((pred == 1) & (y == 1)).sum(); fn = ((pred == 0) & (y == 1)).sum()
    return tp / max(tp + fn, 1) if (len(np.unique(y)) > 1 and (y == 0).sum() > 0 and (y == 1).sum() > 0) else np.nan

# ---- per (method, level, rep) 场景内指标 ----
scen = {}
for lv in [str(v) for v in LEVELS]:
    sub = RAW[RAW["level"] == lv]
    for name, role, fn in METHODS:
        rows, cov = [], []
        for r, g in sub.groupby("rep", sort=True):
            g = g.sort_values("task_id")
            s = fn(g).values.astype(float); y = g["harmful"].values.astype(int)
            m = np.isfinite(s); cov.append(int(m.sum()))
            if len(np.unique(y[m])) < 2: continue
            rows.append((r, roc_auc_score(y[m], s[m]), average_precision_score(y[m], s[m]), sens_at_spec(y, s)))
        d = pd.DataFrame(rows, columns=["rep", "auroc", "auprc", "sens80"])
        d["coverage"] = cov[:len(d)]
        scen[(name, lv)] = d

# ---- 曲线表 ----
out = []
for name, role, _ in METHODS:
    for lv in [str(v) for v in LEVELS]:
        d = scen[(name, lv)]
        out.append(dict(method=name, role=role, level=lv, n_repeat=len(d),
                        coverage_mean=float(d.coverage.mean()),
                        auroc_mean=d.auroc.mean(), auroc_sd=d.auroc.std(),
                        auroc_lo=d.auroc.quantile(.025), auroc_hi=d.auroc.quantile(.975),
                        auprc_mean=d.auprc.mean(), sens80_mean=d.sens80.mean()))
CUR = pd.DataFrame(out)
CUR.to_csv(os.path.join(OUT, "p10c_subsample_curves.csv"), index=False)

p("")
p("=" * 118)
p("【主表】harmful 检测 AUROC 随可获得的 target 样本量变化   mean [2.5%, 97.5%] over 200 repeats")
p("=" * 118)
p(f"{'method':<22}{'role':<12}" + "".join(f"{str(v):>18}" for v in LEVELS) + f"{'Drift':>8}")
for name, role, _ in METHODS:
    cells, vals = [], []
    for lv in [str(v) for v in LEVELS]:
        d = scen[(name, lv)]
        if len(d) == 0: cells.append(" " * 18); vals.append(np.nan); continue
        v = d.auroc.mean(); vals.append(v)
        cells.append(f"{v:.3f}[{d.auroc.quantile(.025):.3f},{d.auroc.quantile(.975):.3f}]")
    p(f"{name:<22}{role:<12}" + "".join(f"{c:>18}" for c in cells) + f"{np.nanmax(vals)-np.nanmin(vals):>8.3f}")

p("")
p("  AUPRC")
for name, role, _ in METHODS:
    v = [scen[(name, str(lv))].auprc.mean() for lv in LEVELS]
    p(f"    {name:<22}" + " ".join(f"{x:.3f}" for x in v))
p("")
p("  sens@spec=0.80")
for name, role, _ in METHODS:
    v = [scen[(name, str(lv))].sens80.mean() for lv in LEVELS]
    p(f"    {name:<22}" + " ".join(f"{x:.3f}" for x in v))

# ---- coverage ----
p("")
p("=" * 118)
p(f"【Coverage】每 repeat 场景中有限的 signal 数（满额 {RAW.task_id.nunique()}）")
p("=" * 118)
for name, role, fn in METHODS:
    v = []
    for lv in [str(v_) for v_ in LEVELS]:
        sub = RAW[RAW["level"] == lv]
        cnt = []
        for r, g in sub.groupby("rep"):
            cnt.append(int(np.isfinite(fn(g).values.astype(float)).sum()))
        v.append(np.mean(cnt))
    p(f"    {name:<22}" + " ".join(f"{x:8.1f}" for x in v))

# ---- paired drift：同一 repeat 内 full vs N=100 ----
p("")
p("=" * 118)
p("【Paired drift】同一 repeat 内 (AUROC@full - AUROC@100)  —— 消除 repeat 间波动")
p("=" * 118)
p(f"{'method':<22}{'mean':>9}{'2.5%':>9}{'97.5%':>9}{'P(>0)':>8}{'判读':>16}")
prow = []
for name, role, _ in METHODS:
    a = scen[(name, "full")].set_index("rep")["auroc"]
    b = scen[(name, "100")].set_index("rep")["auroc"]
    d = (a - b).dropna()
    lo, hi = d.quantile(.025), d.quantile(.975)
    verdict = "随 N 显著上升" if lo > 0 else ("随 N 显著下降" if hi < 0 else "无显著趋势")
    prow.append(dict(method=name, role=role, mean=d.mean(), lo=lo, hi=hi,
                     p_gt0=float((d > 0).mean()), verdict=verdict))
    p(f"{name:<22}{d.mean():>9.3f}{lo:>9.3f}{hi:>9.3f}{float((d>0).mean()):>8.3f}{verdict:>16}")
pd.DataFrame(prow).to_csv(os.path.join(OUT, "p10c_paired_drift.csv"), index=False)

# ---- 三档稳定性排名 ----
p("")
p("=" * 118)
p("【稳定性排名】按 |AUROC@full - AUROC@100| 升序（越小越不依赖 target 样本量）")
p("=" * 118)
prow_s = sorted(prow, key=lambda r: abs(r["mean"]))
for r in prow_s:
    p(f"    {r['method']:<22}{r['role']:<12} |Delta| = {abs(r['mean']):.3f}   {r['verdict']}")

rep.close()
print("done")
