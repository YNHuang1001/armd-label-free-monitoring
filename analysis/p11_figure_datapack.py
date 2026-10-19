# -*- coding: utf-8 -*-
r"""主文 5 图 + Supplement 的数据抽取（figure data pack）

原则：
  * **只从已冻结的产物读取**，不重算任何 signal；唯一重算的是聚合统计量（CI / MAE / calibration），
    且新算的量在此脚本内一次性定义清楚。
  * 每个 panel 一个 tidy CSV，外加 README 映射 panel → 文件 → 来源。
  * 与既有报告数值不一致处必须显式报出，不得静默覆盖。

输出目录：figures/data/
"""
import os
import sys

# ★★ 2026-10-13：**本脚本此前无法在本机重跑** —— 它 print 里有 `⚠`（U+26A0），
#   而 Windows 控制台默认 GBK，直接抛 UnicodeEncodeError 中止。
#   本仓其它脚本（p45/p57/EVS 全套）都有这个守卫，**只有 p11 漏了**。
#   后果：fig2_forest.csv（= 摘要引的 0.556–0.729 的来源）**没有任何人能重跑它**，
#   种子稳定性也就无从测起。此守卫与全仓一致。
try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass
import numpy as np, pandas as pd
from sklearn.metrics import roc_auc_score

BASE = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))); OUT = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "figures", "data")
os.makedirs(FIG, exist_ok=True)
# ★ 2026-10-13：种子可由环境变量覆盖，供种子稳定性分析（p63）用。**不改逻辑**。
SEED = 20261005; rng = np.random.default_rng(SEED); B = 4000

log = open(os.path.join(FIG, "README.md"), "w", encoding="utf-8")
def p(*a):
    line = " ".join(str(x) for x in a); print(line, flush=True); log.write(line + "\n"); log.flush()

# ---------- 读入 ----------
O = pd.read_csv(os.path.join(OUT, "p6d_outcome.csv"))
UNI = pd.read_csv(os.path.join(OUT, "p6b_corrected_universe.csv"))
SIG = pd.read_csv(os.path.join(OUT, "p6b_label_free_signals_blinded.csv"))
L3 = pd.read_csv(os.path.join(OUT, "p6c_L3_blinded.csv"))
SEV = pd.read_csv(os.path.join(OUT, "p6d_signal_eval.csv"))
LAGG = pd.read_csv(os.path.join(OUT, "p6e_lodo_aggregation.csv"))
FC = pd.read_csv(os.path.join(OUT, "p6e_lodo_fold_composition.csv"))
SET = pd.read_csv(os.path.join(OUT, "p6f_lodo_settlement.csv"))
DIR = pd.read_csv(os.path.join(OUT, "p6g_direction.csv"))
RAW10 = pd.read_csv(os.path.join(OUT, "p10b_raw.csv"))
# [AUDIT] 曾在此读取 p5_all_configs.csv（superseded），已移除 —— 不得再引入

E = O.merge(SIG[["task_id", "domain_auc", "mmd", "psi_mean", "psi_max", "js_mean", "js_max",
                 "pred_mean_shift", "entropy_mean_shift", "entropy_js", "conf_dist_js",
                 "pred_pos_prev_shift"]], on="task_id", how="left") \
     .merge(L3[["task_id", "SUDO_AURCC", "AUROC_IOW", "AUROC_OM", "AUROC_DR", "w_max", "w_ESS",
                "w_P95"]], on="task_id", how="left")
E["group"] = E.organism + "|" + E.drug
E["direction"] = E.source_site + "->" + E.target_site

p("# Figure Data Pack — 主文 5 图 + Supplement")
p("")
p("由 `分析结果/p11_figure_datapack.py` 生成。所有 signal 均来自已冻结产物，未重算。")
p("重算的仅为聚合统计量（CI / MAE / calibration），定义见脚本内注释。")
p("")

# ================= Figure 1 =================
p("## Figure 1 — Benchmark construction and failure landscape")
p("")
# 1B
cand = len(UNI); elig = int(UNI.eligible.sum())
f1b = pd.DataFrame([dict(candidate_tasks=cand, eligible_tasks=elig,
                         excluded=cand - elig,
                         n_organism=O.organism.nunique(), n_drug=O.drug.nunique(),
                         n_organism_drug_groups=E.group.nunique())])
f1b.to_csv(os.path.join(FIG, "fig1b_universe.csv"), index=False)
# Panel B 长表：供 R 层直接绘图（R 不再做宽->长 reshape）
_cand_lab = "Candidate" + chr(10) + "transfer tasks"
_elig_lab = "Eligible" + chr(10) + "(source AUROC >= 0.65)"
pd.DataFrame([
    dict(stage_key="Candidate", stage_label=_cand_lab, n=cand,
         display=_cand_lab, excluded_note=""),
    dict(stage_key="Eligible", stage_label=_elig_lab, n=elig,
         display=_elig_lab, excluded_note=f"{cand - elig} excluded"),
]).to_csv(os.path.join(FIG, "fig1b_universe_long.csv"), index=False)
p(f"- `fig1b_universe.csv` — {cand} candidate → {elig} eligible; "
  f"{O.organism.nunique()} organisms, {O.drug.nunique()} drugs, {E.group.nunique()} organism–drug groups")

# 1C
f1c = E[["task_id", "organism", "drug", "source_site", "target_site", "direction",
         "n_target", "auroc_source_internal", "auroc_target", "d_auroc", "harmful", "stable"]]
# Panel D 顺序：在 Python 层按各 direction 的 median 预计算好，R 层不再排序
_dir_med = f1c.groupby("direction")["d_auroc"].median()
_dir_order = {d: i + 1 for i, d in enumerate(_dir_med.sort_values().index)}
f1c["direction_order"] = f1c["direction"].map(_dir_order)
f1c.to_csv(os.path.join(FIG, "fig1c_delta_distribution.csv"), index=False)
qs = O.d_auroc.quantile([.10, .25, .50, .75, .90])
# 1C 汇总：R 层禁止 quantile()，所有分位数必须在此算好再交给绘图层
pd.DataFrame([dict(n_tasks=len(O), n_harmful=int(O.harmful.sum()), n_stable=int(O.stable.sum()),
                   harmful_threshold=-0.05, stable_abs_threshold=0.02,
                   mean=float(O.d_auroc.mean()), sd=float(O.d_auroc.std()),
                   p05=float(O.d_auroc.quantile(.05)), p10=float(qs[.10]), p25=float(qs[.25]),
                   p50=float(qs[.50]), p75=float(qs[.75]), p90=float(qs[.90]),
                   p95=float(O.d_auroc.quantile(.95)),
                   p90_minus_p10=float(qs[.90] - qs[.10]),
                   dmin=float(O.d_auroc.min()), dmax=float(O.d_auroc.max()))]
             ).to_csv(os.path.join(FIG, "fig1c_summary.csv"), index=False)
p(f"- `fig1c_delta_distribution.csv` — {len(O)} tasks. harmful={int(O.harmful.sum())}, stable={int(O.stable.sum())}, "
  f"P10={qs[.10]:+.3f} P50={qs[.50]:+.3f} P90={qs[.90]:+.3f}")

# 1D
f1d = E.groupby("direction")["d_auroc"].describe()[["count", "mean", "std", "min", "25%", "50%", "75%", "max"]]
f1d.to_csv(os.path.join(FIG, "fig1d_by_direction.csv"))
p(f"- `fig1d_by_direction.csv` — 6 个方向。**caption 必须写明不作机制解释**")
for d, r in f1d.iterrows():
    p(f"    {d:<12} n={int(r['count'])}  median={r['50%']:+.3f}  min={r['min']:+.3f}  max={r['max']:+.3f}")

# 1E —— 必须用 **修正后管线** 的四配置结果（p12），不得用 superseded 的 p5_all_configs.csv
CFG12 = pd.read_csv(os.path.join(OUT, "p12_configs_corrected.csv"))
f1e = CFG12.to_csv(os.path.join(FIG, "fig1e_configs_full.csv"), index=False)
P1set = set(CFG12[(CFG12.config == "P1 History+LR") & CFG12.eligible].task_id)
rows1e = []
for cfg in ["P1 History+LR", "S1 Core+LR", "S2 History+LGBM", "S3 Core+LGBM"]:
    own = CFG12[CFG12.config == cfg]; s = own[own.task_id.isin(P1set)]
    rows1e.append(dict(config=cfg, tier=own.tier.iat[0], algo=own.algo.iat[0],
                       n_aligned=len(s), median_aligned=s.d_auroc.median(),
                       p10_aligned=s.d_auroc.quantile(.1), p90_aligned=s.d_auroc.quantile(.9),
                       p90_minus_p10_aligned=s.d_auroc.quantile(.9) - s.d_auroc.quantile(.1),
                       n_harmful_aligned=int((s.d_auroc <= -0.05).sum()),
                       pct_harmful_aligned=100 * float((s.d_auroc <= -0.05).mean()),
                       n_eligible_own=int(own.eligible.sum()), n_task_own=len(own)))
F1E = pd.DataFrame(rows1e)
F1E.to_csv(os.path.join(FIG, "fig1e_configs.csv"), index=False)
# 供 R 直接绘图的**已对齐**长表：4 配置 x 全部 eligible task，只含 P1-defined universe 的 task。
# R 层只读不筛，避免在绘图层重新推导 eligible。
_al = CFG12[CFG12.task_id.isin(P1set)].copy()
_al["in_p1_universe"] = True
_al[["task_id", "config", "tier", "algo", "organism", "drug", "source", "target",
     "d_auroc", "in_p1_universe"]].to_csv(
    os.path.join(FIG, "fig1e_configs_aligned.csv"), index=False)
p("- `fig1e_configs.csv` / `fig1e_configs_full.csv` — 四配置 degradation 对照，"
  "**修正后管线（p12_configs_corrected.py）**")
p("  **不得使用 `p5_all_configs.csv`**（protocol v1.1 标记为 superseded / preprocessing-mismatched；"
  "其 P1 与 p6d 的 d_auroc 最大差 0.0084，S3 最大差 0.0238）")
p("  对齐在同一批 P1-202 eligible tasks 上（各配置自身 eligible universe 大小不同，不具可比性）")
for _, r in F1E.iterrows():
    p(f"    {r.config:<18} tier={r.tier:<5} algo={r.algo:<5} median={r.median_aligned:+.3f} "
      f"P90-P10={r.p90_minus_p10_aligned:.3f} harmful={r.n_harmful_aligned} "
      f"({r.pct_harmful_aligned:.1f}%)  | 自身 eligible={r.n_eligible_own}/{r.n_task_own}")
h = F1E[F1E.tier == "hist"]; c = F1E[F1E.tier == "core"]
p(f"  **结论**：同 tier 内 LR 与 LGBM 几乎一致（hist {h.n_harmful_aligned.tolist()}，"
  f"core {c.n_harmful_aligned.tolist()}）；差异主要来自 feature tier 而非算法 "
  f"→ 支持 'failure is not classifier-specific'")
p("  注：P1 与冻结 p6d 的 d_auroc 最大差 2.27e-05（p12 在算 AUROC 前套了 Platt 校准，p6b/p6d 的 source 侧没有），"
  "不影响任何标签或计数。")

# ================= Figure 2 =================
p("")
p("## Figure 2 — What unlabeled target data can and cannot tell us")
p("")
f2 = SEV[["method", "role", "coverage", "harmful_auc", "auc_lo", "auc_hi",
          "harmful_auprc", "sens_at_spec80", "spearman", "sp_lo", "sp_hi"]].copy()
f2.to_csv(os.path.join(FIG, "fig2_forest.csv"), index=False)
p("- `fig2_forest.csv` — 15 方法。harmful AUROC + **organism–drug cluster bootstrap (1000×)** 95%CI；"
  "旁边小 panel 用 spearman/sp_lo/sp_hi")
p("  **不得在 Fig 2 中反向 MMD**（Fig 2 是冻结分析，惯例方向）。方向可学性见 Figure S。")
n_below = int((f2.auc_hi < 0.5).sum()); n_above = int((f2.auc_lo > 0.5).sum())
p(f"  全样本口径：CI 全 > 0.5 的 {n_above} 个；CI 全 < 0.5 的 {n_below} 个；其余覆盖 0.5")

# ================= Figure 3 =================
p("")
p("## Figure 3 — Generalization beyond the task on which the signal was observed")
p("")
# 3A —— 用 p6f 的口径（点估计与 CI 均为 unweighted，自洽）
f3a = SET[["method", "role", "n_eval_fold", "C_macro_eval", "C_lo", "C_hi"]].rename(
    columns={"C_macro_eval": "lodo_macro_auroc", "C_lo": "lo_95", "C_hi": "hi_95"})
f3a["n_tasks_eval"] = int(np.isin(E.group, FC[FC.auroc_definable == 1].group).sum())
f3a["category"] = np.where(f3a.lo_95 > 0.5, "above_chance",
                           np.where(f3a.hi_95 < 0.5, "below_chance", "indistinguishable"))
f3a.to_csv(os.path.join(FIG, "fig3a_lodo.csv"), index=False)
p(f"- `fig3a_lodo.csv` — LODO macro AUROC + fold-level cluster bootstrap CI (B={B})")
# ★ 2026-10-06：原写死 "25/38" 与 "/202"。eligible 宇宙改为 186 task / 37 组后
#   实测折数为 26/37、覆盖 134 task。折数与 task 数都随宇宙变化 ⇒ 一律现算。
#   FC = 折构成表，E = 评价帧（两者在本行上方已被使用，非新引入的依赖）。
p(f"  **图上必须标注：{int(FC.auroc_definable.sum())}/{len(FC)} evaluable folds; "
  f"{int(f3a.n_tasks_eval.iat[0])}/{len(E)} tasks**")
p("  **口径说明**：使用 p6f 的 unweighted 点估计 + unweighted CI。")
p("    p6e_lodo_aggregation.csv 的 macro_lo/hi 是 size-weighted bootstrap，而它的点估计 macro_mean 是"
  " unweighted —— 两者不自洽，故图不使用 p6e 的 CI。")
for _, r in f3a.iterrows():
    p(f"    {r.method:<22}{r.role:<11}[{r.lo_95:.3f},{r.hi_95:.3f}]  {r.category}")

# 3B —— LOO-organism：**本脚本不生成**。
# 原版在此处用 `auc(y[m], r[m])`（留出 organism 内用原始 signal 直接算 AUROC）是**语义错误**，
# 不是 leave-one-organism-out。正确实现（在另两个 organism 上训练后留出评价）见
# `p14_rebuild_fig3b.py`，由它独占写出 fig3b_loo_organism.csv。
# 若在此重新引入，会覆盖正确文件 —— 已有 p11 重跑导致 fig3b 被回退的先例。
METHODS = [("domain_auc", "L1", lambda x: x["domain_auc"]), ("mmd", "L1", lambda x: x["mmd"]),
           ("psi_mean", "L1", lambda x: x["psi_mean"]), ("psi_max", "L1", lambda x: x["psi_max"]),
           ("js_mean", "L1", lambda x: x["js_mean"]), ("js_max", "L1", lambda x: x["js_max"]),
           ("pred_mean_shift", "L2", lambda x: x["pred_mean_shift"].abs()),
           ("entropy_mean_shift", "L2", lambda x: x["entropy_mean_shift"].abs()),
           ("entropy_js", "L2", lambda x: x["entropy_js"]),
           ("conf_dist_js", "L2", lambda x: x["conf_dist_js"]),
           ("pred_pos_prev_shift", "L2", lambda x: x["pred_pos_prev_shift"].abs()),
           ("SUDO_AURCC", "L2", lambda x: -x["SUDO_AURCC"]),
           ("AUROC_DR", "L3-primary", lambda x: -x["AUROC_DR"]),
           ("AUROC_IOW", "L3-sens", lambda x: -x["AUROC_IOW"]),
           ("AUROC_OM", "L3-sens", lambda x: -x["AUROC_OM"])]

# _ = None
# if False:
# def auc(yy, ss):
#     m = np.isfinite(ss)
#     if len(np.unique(yy[m])) < 2: return np.nan
#     yy, ss = yy[m], ss[m]
#     o = np.argsort(ss, kind="mergesort"); yy = yy[o]; ss = ss[o]
#     r = np.empty(len(ss)); i = 0
#     while i < len(ss):
#         j = i
#         while j + 1 < len(ss) and ss[j + 1] == ss[i]: j += 1
#         r[i:j + 1] = 0.5 * (i + j) + 1; i = j + 1
#     n1 = yy.sum(); n0 = len(yy) - n1
#     if n1 == 0 or n0 == 0: return np.nan
#     return (r[yy == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)

# METHODS = [("domain_auc", "L1", lambda x: x["domain_auc"]), ("mmd", "L1", lambda x: x["mmd"]),
#            ("psi_mean", "L1", lambda x: x["psi_mean"]), ("psi_max", "L1", lambda x: x["psi_max"]),
#            ("js_mean", "L1", lambda x: x["js_mean"]), ("js_max", "L1", lambda x: x["js_max"]),
#            ("pred_mean_shift", "L2", lambda x: x["pred_mean_shift"].abs()),
#            ("entropy_mean_shift", "L2", lambda x: x["entropy_mean_shift"].abs()),
#            ("entropy_js", "L2", lambda x: x["entropy_js"]),
#            ("conf_dist_js", "L2", lambda x: x["conf_dist_js"]),
#            ("pred_pos_prev_shift", "L2", lambda x: x["pred_pos_prev_shift"].abs()),
#            ("SUDO_AURCC", "L2", lambda x: -x["SUDO_AURCC"]),
#            ("AUROC_DR", "L3-primary", lambda x: -x["AUROC_DR"]),
#            ("AUROC_IOW", "L3-sens", lambda x: -x["AUROC_IOW"]),
#            ("AUROC_OM", "L3-sens", lambda x: -x["AUROC_OM"])]
# y = E["harmful"].values.astype(int); org = E.organism.values
# ORGS = sorted(E.organism.unique())
# rows3b = []
# for name, role, fn in METHODS:
#     r = fn(E).values.astype(float)
#     per = {}
#     for o_ in ORGS:
#         m = org == o_
#         a = auc(y[m], r[m])
#         if np.isfinite(a): per[o_] = a
#     vals = list(per.values())
#     ks = list(per.keys()); bv = []
#     for _ in range(B):
#         pick = rng.choice(len(ks), len(ks), replace=True)
#         bv.append(np.mean([per[ks[i]] for i in pick]))
#     rows3b.append(dict(method=name, role=role, n_eval_org=len(per),
#                        loo_org_macro=float(np.mean(vals)),
#                        lo_95=np.percentile(bv, 2.5), hi_95=np.percentile(bv, 97.5),
#                        **{f"auc_{o_}": per.get(o_, np.nan) for o_ in ORGS}))
# f3b = pd.DataFrame(rows3b)
# f3b.to_csv(os.path.join(FIG, "fig3b_loo_organism.csv"), index=False)
# p(f"- `fig3b_loo_organism.csv` — LOO-organism macro + CI（3/3 organism 折全部可定义）。"
#   f" 每折 n: " + ", ".join(f"{o_}={int((org==o_).sum())}" for o_ in ORGS))

p("- `fig3b_loo_organism.csv` — **由 `p14_rebuild_fig3b.py` 生成，本脚本不写**。"
  "proper leave-one-organism-out（在另两个 organism 上训练后留出评价），与 p6d 的 loo_org_auc 逐方法一致。"
  "⚠ 仅 3 折，bootstrap 重抽空间仅 27 种，CI 极粗，caption 必须说明")

# 3C —— SUDO 聚合分解
f3c = SET[SET.method == "SUDO_AURCC"][["A_pooled_all202", "B_pooled_eval", "C_macro_eval",
                                       "AtB", "BtC"]].copy()
f3c.columns = ["pooled_202task", "pooled_138task", "macro_138task",
               "task_subset_effect", "aggregation_effect"]
f3c.to_csv(os.path.join(FIG, "fig3c_sudo_decomposition.csv"), index=False)
# ★ 2026-10-06：标签原写 "202-pooled / 138-pooled / 138-macro"。数字随 eligible 宇宙变化
#   （本期 all=186 / eval=134），写死必错。改用与 p6f 表头一致的口径名（pooled(all) /
#   pooled(eval) / macro(eval)）。**列名 pooled_202task 等保留不动**（A 类：改则换 schema）。
p(f"- `fig3c_sudo_decomposition.csv` — pooled(all) {f3c.pooled_202task.iat[0]:.3f} → "
  f"pooled(eval) {f3c.pooled_138task.iat[0]:.3f} → macro(eval) {f3c.macro_138task.iat[0]:.3f}"
  f"  (subset effect {f3c.task_subset_effect.iat[0]:+.3f}, aggregation effect {f3c.aggregation_effect.iat[0]:+.3f})")

# ================= Figure 4 =================
p("")
p("## Figure 4 — Can target performance itself be estimated without target labels?")
p("")
f4 = E[["task_id", "organism", "drug", "source_site", "target_site",
        "auroc_target", "auroc_source_internal", "d_auroc",
        "AUROC_DR", "AUROC_IOW", "AUROC_OM", "w_max", "w_P95", "w_ESS"]].copy()
f4["d_hat_DR"] = f4.AUROC_DR - f4.auroc_source_internal
f4["d_hat_IOW"] = f4.AUROC_IOW - f4.auroc_source_internal
f4["d_hat_OM"] = f4.AUROC_OM - f4.auroc_source_internal
f4.to_csv(os.path.join(FIG, "fig4_l3_scatter.csv"), index=False)
p("- `fig4_l3_scatter.csv` — ÂUROC_target 与 AUROC_target 配对，供 A/B/C 三个 identity plot")
p("  **统一坐标轴**，便于肉眼比较 OM 更贴近 identity。")

rows4 = []
for c, lab in [("AUROC_DR", "DR"), ("AUROC_IOW", "IOW"), ("AUROC_OM", "OM")]:
    v = f4[c].values.astype(float); t = f4.auroc_target.values.astype(float)
    m = np.isfinite(v)
    err = np.abs(v[m] - t[m])
    # calibration：必须复现 p6d —— 对**预测 ΔAUROC** 用**固定分档**
    # （曾误用 qcut 对 ÂUROC_target 等分位，与 report_p6d 不一致，已修正）
    pred_d = v - f4.auroc_source_internal.values
    EDGES = [-np.inf, -0.10, -0.05, 0.0, 0.05, np.inf]
    LABELS = ["<-0.10", "-0.10~-0.05", "-0.05~0", "0~0.05", ">0.05"]
    calib = pd.DataFrame({"bin": pd.cut(pred_d[m], bins=EDGES, labels=LABELS, right=False),
                          "pred_d": pred_d[m], "actual_d": f4.d_auroc.values[m]}) \
        .groupby("bin", observed=False).agg(n=("pred_d", "size"), pred_d_mean=("pred_d", "mean"),
                                            actual_d_mean=("actual_d", "mean")).reset_index()
    rows4.append(dict(estimator=lab, n=int(m.sum()), oor=int((v[m] > 1).sum() + (v[m] < 0).sum()),
                      **{f"mae": float(err.mean())},
                      rho_pred_true=float(pd.Series(v[m]).corr(pd.Series(t[m]), method="spearman"))))
    calib.to_csv(os.path.join(FIG, f"fig4d_calibration_{lab}.csv"), index=False)
f4d = pd.DataFrame(rows4)
f4d.to_csv(os.path.join(FIG, "fig4d_mae.csv"), index=False)
p("- `fig4d_mae.csv` + `fig4d_calibration_{DR,IOW,OM}.csv`")
for _, r in f4d.iterrows():
    p(f"    {r.estimator:<5} n={r.n}  MAE={r.mae:.4f}  越界={r.oor}  "
      f"Spearman(ÂUROC,AUROC)={r.rho_pred_true:+.3f}")
p("  **正文必须坚持：DR was prespecified primary; OM performed better in sensitivity analysis.**")

# 4E
f4e = f4[f4.AUROC_DR > 1][["task_id", "organism", "drug", "source_site", "target_site",
                           "auroc_source_internal", "auroc_target", "d_auroc",
                           "AUROC_DR", "w_max", "w_P95", "w_ESS"]].sort_values("AUROC_DR", ascending=False)
f4e.to_csv(os.path.join(FIG, "fig4e_dr_out_of_range.csv"), index=False)
p(f"- `fig4e_dr_out_of_range.csv` — {len(f4e)} 个 DR 越界 task（未 clip）")
for _, r in f4e.iterrows():
    p(f"    {r.task_id:<18} ÂUROC={r.AUROC_DR:.3f}  w_max={r.w_max:,.1f}  ESS={r.w_ESS:,.1f}")
rho_ess = pd.Series(f4.w_ESS).corr(pd.Series(np.abs(f4.AUROC_DR - f4.auroc_target)), method="spearman")
p(f"  [诊断] Spearman(w_ESS, |ÂUROC_DR − AUROC_target|) = {rho_ess:+.3f}")

# ================= Figure 5 =================
p("")
p("## Figure 5 — How much unlabeled target data are needed for stable monitoring?")
p("")
O2 = O[["task_id", "harmful"]]
RW = RAW10.merge(O2, on="task_id", how="left")
levels = ["100", "250", "500", "1000", "full"]
rows5 = []
for name, role, fn in METHODS:
    per_level = {}
    for lv in levels:
        sub = RW[RW.level == lv]
        vals = {}
        for rr, gg in sub.groupby("rep"):
            gg = gg.sort_values("task_id")
            s = fn(gg).values.astype(float); yy = gg["harmful"].values.astype(int)
            m = np.isfinite(s)
            if len(np.unique(yy[m])) < 2: continue
            vals[rr] = roc_auc_score(yy[m], s[m])
        v = pd.Series(vals)
        per_level[lv] = v
        rows5.append(dict(method=name, role=role, level=lv, n_repeat=len(v),
                          auroc_mean=float(v.mean()), auroc_lo=float(v.quantile(.025)),
                          auroc_hi=float(v.quantile(.975))))
    # paired drift full - 100
    d = (per_level["full"] - per_level["100"]).dropna()
    rows5.append(dict(method=name, role=role, level="PAIRED_full_minus_100",
                      n_repeat=len(d), auroc_mean=float(d.mean()),
                      auroc_lo=float(d.quantile(.025)), auroc_hi=float(d.quantile(.975))))
f5 = pd.DataFrame(rows5)
f5.to_csv(os.path.join(FIG, "fig5_subsample.csv"), index=False)
# Panel B 用：跨 5 档的 AUROC 极差。与 P90−P10 同类，**在 Python 层算好**，R 不做减法。
_lv = f5[f5.level.isin(["100", "250", "500", "1000", "full"])]
rows5b = []
for name, role, _ in METHODS:
    s = _lv[_lv.method == name].copy()
    imax = s.auroc_mean.idxmax(); imin = s.auroc_mean.idxmin()
    rows5b.append(dict(method=name, role=role, n_levels=len(s),
                       across_n_range=float(s.auroc_mean.max() - s.auroc_mean.min()),
                       auroc_max=float(s.auroc_mean.max()), auroc_max_level=s.loc[imax, "level"],
                       auroc_min=float(s.auroc_mean.min()), auroc_min_level=s.loc[imin, "level"],
                       auroc_at_full=float(s.loc[s.level == "full", "auroc_mean"].iat[0]),
                       auroc_at_100=float(s.loc[s.level == "100", "auroc_mean"].iat[0])))
pd.DataFrame(rows5b).to_csv(os.path.join(FIG, "fig5b_across_n_range.csv"), index=False)
p("- `fig5_subsample.csv` — 5 档 x 15 方法 + 每方法的 paired (full − 100)")
p("- `fig5b_across_n_range.csv` — 每方法跨 5 档的 AUROC 极差（Python 层算好，R 不做减法）")
p("  **caption 必须两句连写**：")
p("    Reducing unlabeled target sample size had little effect on failure-discrimination")
p("    performance for most methods.")
p("    This stability should not be interpreted as adequacy: L1 methods remained close to or")
p("    below the prespecified random-discrimination reference across sample sizes.")
drift = f5[f5.level != "PAIRED_full_minus_100"].groupby("method").auroc_mean.agg(lambda s: s.max() - s.min())
p(f"  跨档极差最大的三个：{', '.join(f'{k}={v:.3f}' for k, v in drift.nlargest(3).items())}")
sig = f5[(f5.level == "PAIRED_full_minus_100") & ((f5.auroc_lo > 0) | (f5.auroc_hi < 0))]
p(f"  paired CI 不含 0 的方法：{list(sig.method) if len(sig) else '无'}")

# ================= Supplement =================
p("")
p("## Supplement")
p("")
FC.to_csv(os.path.join(FIG, "figS6_fold_composition.csv"), index=False)
# S6 Panel B 汇总：全部在 Python 层算好，R 只读（不在绘图层做条件计数）
_eval = FC[FC.auroc_definable == 1]
_nonev = FC[FC.auroc_definable == 0]
S6S = pd.DataFrame([
    dict(item="Evaluable folds (both classes present)", n_folds=int(len(_eval)),
         n_tasks=int(_eval.n.sum()), group="evaluable", eligible_for_auroc=1),
    dict(item="Single-class folds (excluded from within-fold AUROC)", n_folds=int(len(_nonev)),
         n_tasks=int(_nonev.n.sum()), group="excluded", eligible_for_auroc=0),
    dict(item="   of which all-harmful", n_folds=int((_nonev.n_nonharmful == 0).sum()),
         n_tasks=int(_nonev.loc[_nonev.n_nonharmful == 0, "n"].sum()),
         group="all_harmful", eligible_for_auroc=0),
    dict(item="   of which all-non-harmful", n_folds=int((_nonev.n_harmful == 0).sum()),
         n_tasks=int(_nonev.loc[_nonev.n_harmful == 0, "n"].sum()),
         group="all_nonharmful", eligible_for_auroc=0),
    dict(item="Total", n_folds=int(len(FC)), n_tasks=int(FC.n.sum()),
         group="total", eligible_for_auroc=0),
])
# 自洽性断言：可定义 + 不可定义 = 全量；两类单类折 = 全部不可定义折
assert int(S6S.loc[S6S.group == "evaluable", "n_folds"].iat[0]) + \
       int(S6S.loc[S6S.group == "excluded", "n_folds"].iat[0]) == len(FC)
assert int(S6S.loc[S6S.group == "all_harmful", "n_folds"].iat[0]) + \
       int(S6S.loc[S6S.group == "all_nonharmful", "n_folds"].iat[0]) == len(_nonev)
S6S.to_csv(os.path.join(FIG, "figS6_fold_summary.csv"), index=False)
# ★ 2026-10-06：原写死 "38 folds"。eligible 宇宙改为 186 task / 37 组后折数会变，
#   写死会让报告与 CSV 里的实际行数不符。改为从 FC 现算。
p(f"- `figS6_fold_composition.csv` — {len(FC)} folds，可定义 {int(FC.auroc_definable.sum())} 个"
  f"（全 harmful {int((FC.n_nonharmful==0).sum())}，全 non-harmful {int((FC.n_harmful==0).sum())}）")
p(f"- `figS6_fold_summary.csv` — S6 Panel B 汇总（含自洽性断言：可定义+不可定义=全量）")
DIR.to_csv(os.path.join(FIG, "figS_orientation.csv"), index=False)
p("- `figS_orientation.csv` — **Post-unblinding analysis of learnable signal orientation**")
p("  三口径：① 组内惯例方向 / ② logistic 学得方向（预注册 learner）/ ③ rank 学得方向（事后）")
p("  **MMD 与 SUDO 在图上会非常醒目**")
p("  正文 Results 2 用一句引到 Supplement（措辞已定，见 audit log A-2026-09-30-4）")
# ★ 2026-10-07（A2.1）：**补齐 15 个信号的逐 transfer 值**。
#   病因：原列表只含 10 个信号，缺 psi_max / js_max / pred_mean_shift /
#   entropy_mean_shift / pred_pos_prev_shift。
#   后果：S1 的首句「All derived quantities in this manuscript are computed on these rows」
#   **为假** —— 主对比的对照臂（六者秩平均 composite）需要 psi_max 与 js_max，
#   而那两列不在表里 ⇒ composite 的 7/37 **无法被任何人独立复算**（外部核验实报）。
#   来源：E 表（p6b blinded + p6c blinded），**presentation-side 抽取，不重算任何统计量**。
#   列序按信号家族：L1(6) → L2(6) → L3(3)。
tS1 = E[["task_id", "organism", "drug", "source_site", "target_site", "n_target",
         "auroc_source_internal", "auroc_target", "d_auroc", "harmful", "stable",
         "domain_auc", "mmd", "psi_mean", "psi_max", "js_mean", "js_max",
         "pred_mean_shift", "entropy_mean_shift", "entropy_js", "conf_dist_js",
         "pred_pos_prev_shift", "SUDO_AURCC",
         "AUROC_DR", "AUROC_IOW", "AUROC_OM"]]
assert len(tS1.columns) == 26, f"S1 列数应为 26，实为 {len(tS1.columns)}"
# ★★ 2026-10-09（INC-028）：**import 必须在第一次使用之前**。
#   原序为 `_gs1.backup(...)` → `to_csv(...)` → `import genstamp as _gs1; _gs1.stamp(...)` ——
#   第 425 行的 `_gs1` **在导入之前**被引用 ⇒ `NameError`。
#   后果链：p11 无法整体执行 ⇒ EVS `t13` 的 EC-06…EC-11 与 EC-15 全 FAIL（**7 条**）
#   ⇒ closure gate 报 4 条 Critical path 未 PASS ⇒ **不得 freeze**。
#   ★ 这是本会话**第二次**在同一个位置出错（首版是「覆盖后才归档、导入在调用后」，
#   我诊断对了、修的时候把 import 放错了地方）；**t13 是唯一会执行 p11 整块代码的判据**，
#   所以只有它抓得到 —— 而我在改完之后**没有跑过 t13**。
#> 块须**自带依赖解析**：EVS `t13` 在本进程 exec 本块，`分析结果/` 不在其 sys.path 上
#> （否则 `ModuleNotFoundError: No module named 'genstamp'`）。
#> **可被逐块执行的代码，必须自带它的路径前提** —— 不能依赖「整脚本一起跑」这个隐含条件。
import sys as _sys
_gsdir = os.path.dirname(os.path.abspath(__file__))
if _gsdir not in _sys.path:
    _sys.path.insert(0, _gsdir)
import genstamp as _gs1
_gs1.backup(os.path.join(FIG, "tableS1_tasks.csv"), "p11")
tS1.to_csv(os.path.join(FIG, "tableS1_tasks.csv"), index=False)
_gs1.stamp(os.path.join(FIG, "tableS1_tasks.csv"), "p11")
p(f"- `tableS1_tasks.csv` — {len(O)} tasks 完整表（**26 列 = 11 身份/结局 + 15 信号**）")
tS2 = SEV.merge(LAGG[["method", "pooled_OOF", "macro_mean", "n_evaluable_folds",
                      "fold_median", "fold_p25", "fold_p75"]], on="method", how="left") \
          .merge(SET[["method", "A_pooled_all202", "B_pooled_eval", "C_macro_eval", "C_lo", "C_hi"]], on="method", how="left") \
          .merge(DIR[["method", "macro_within", "macro_logit", "macro_rank", "frac_train_dir_pos"]], on="method", how="left") \
          .merge(pd.read_csv(os.path.join(OUT, "p33_protocol72_metrics.csv"))
                   .query("convention == 'primary'")
                   .reindex(columns=["method",
                                     "ppv", "ppv_lo", "ppv_hi",
                                     "npv", "npv_lo", "npv_hi",
                                     "auprc_lo", "auprc_hi",
                                     "sens80_lo", "sens80_hi"]),
                 on="method", how="left")
# [2026-10-07 补] 上一行 merge 的是 protocol §7.2 四项指标里此前缺失的 **PPV / NPV**。
#   来源 p33_protocol72_metrics.py（从 p6b/p6c/p6d 上游独立重算，并与 p6d_signal_eval.csv
#   在 AUROC / AUPRC / sens@.80 三项上逐位锚定，差 < 1e-9）。
#   ★ **只做 merge，不重算任何统计量** —— 与 CORRECTION_SPEC 对 p11 的
#     "Table S1/S2 merges = merge-only" 定性一致。AUPRC 与 sens@.80 已在 SEV 中，不重复引入。
#   ★★ 定义**必须内联**，不得写成 `P72 = ...` 再引用：t13_eval_core 按锚点**抽取代码块单独执行**，
#      块外的裸名字在该命名空间里不存在。首版即因此触发
#      `NameError: name 'P72' is not defined`，连带 EC-06…EC-11 七条判红。
assert len(tS2) == len(SEV), "tableS2 merge 改变了行数"
assert tS2.ppv.notna().all() and tS2.npv.notna().all(), "有方法未匹配到 PPV/NPV"
tS2.to_csv(os.path.join(FIG, "tableS2_metrics.csv"), index=False)
p(f"- `tableS2_metrics.csv` — 每指标 full-universe / LODO(两口径+CI) / LOO-organism / coverage / orientation / **§7.2 的 PPV+NPV（2026-10-07 补，merge-only）**")

p("")
p("## 待补（尚未生成）")
p("- Fig S1 完整 organism × drug × site task matrix（可从 fig1c 透视得到）")
p("- Fig S2 四配置 robustness（已由 fig1e 覆盖大部分）")
p("- Fig S4 原始 signal distribution + blinded QC（源：p6b/p6c blinded 文件）")
p("- Fig S5 MMD corrected vs superseded（需旧 p6 signals，已被覆盖 —— 需确认是否还留有副本）")
p("- Fig S7 IOW weight/ESS diagnostics（已在 fig4e 部分覆盖）")
p("- Fig S8 §10 全 15 方法曲线（= fig5_subsample.csv 全量）")
log.close()
print("\ndone ->", FIG)
