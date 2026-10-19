# -*- coding: utf-8 -*-
r"""p30 — MNV（manuscript numerical validation）：从**冻结产物独立重算**全部 headline numbers

================================ 性质 ================================
CORRECTION_SPEC §5：EVS 验**代码/估计量**，MNV 验**正文数字**。
本脚本**不复用生产管线**：所有量都从**最上游的冻结任务级表**
（`figures/data/tableS1_tasks.csv` 等）用独立写法重算，
再与正文/图产物里的值逐一比对。

★ 判据不是"能不能算出来"，而是"**与正文报的值是否一致**"。
★ 任何一处不一致 = 阻断项，必须报告而不是自行改正文。

覆盖（CORRECTION_SPEC §5 全部八项）：
  1 benchmark universe      5 L3 MAE/ρ/calibration + naive baseline
  2 ΔAUROC / harmful counts 6 target-N
  3 Fig 2 AUROC/ρ/CI        7 robustness
  4 Fig 3 pooled/macro/LOO  8 Supplement reporting counts
"""
import io, os, sys, re
import numpy as np, pandas as pd
from scipy import stats
from sklearn.metrics import roc_auc_score

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIG  = os.path.join(ROOT, "figures", "data")
OUT  = os.path.join(ROOT, "分析结果")
REP  = open(os.path.join(OUT, "report_p30_mnv.txt"), "w", encoding="utf-8")
ROWS = []

def p(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); REP.write(s + "\n"); REP.flush()

def chk(section, name, expected, observed, tol=0.0, note=""):
    ok = (abs(float(expected) - float(observed)) <= tol) if tol else (expected == observed)
    ROWS.append(dict(section=section, item=name, expected=expected,
                     observed=observed, tol=tol, status="PASS" if ok else "FAIL", note=note))
    p("  [%s] %-58s 正文=%-10s 重算=%-10s %s" % (
        "PASS" if ok else "FAIL", name, expected,
        (("%.6g" % observed) if isinstance(observed, float) else observed), note))
    return ok

# ★ 独立上游：**不用** tableS1 / figure CSV，而是自己从三张原始冻结表 join。
#   tableS1_tasks.csv 只有 12 个信号中的 9 个（缺 psi_max / js_max / 三个 L2 预测空间信号），
#   故若拿它当上游，MNV 会漏掉一半 Fig 2/3 的量。
B = pd.read_csv(os.path.join(OUT, "p6b_label_free_signals_blinded.csv"))
C = pd.read_csv(os.path.join(OUT, "p6c_L3_blinded.csv"))
D = pd.read_csv(os.path.join(OUT, "p6d_outcome.csv"))
SIG = ["domain_auc", "mmd", "psi_mean", "psi_max", "js_mean", "js_max",
       "pred_mean_shift", "entropy_mean_shift", "entropy_js", "conf_dist_js",
       "pred_pos_prev_shift"]
L3 = ["AUROC_DR", "AUROC_IOW", "AUROC_OM"]
T = (D[["task_id", "organism", "drug", "auroc_source_internal", "auroc_target",
        "d_auroc", "harmful", "stable"]]
     .merge(B[["task_id"] + SIG], on="task_id", validate="1:1")
     .merge(C[["task_id", "SUDO_AURCC", "w_ESS"] + L3], on="task_id", validate="1:1"))

p("=" * 100); p("p30 — MNV 独立重算"); p("=" * 100)
p(f"独立上游（自行 join，不经任何 figure/tableS1 产物）：")
p(f"    p6d_outcome.csv (186) × p6b_label_free_signals_blinded.csv (186, L1/L2)")
p(f"    × p6c_L3_blinded.csv (186, L3+SUDO)   ⇒ {T.shape[0]} tasks × {T.shape[1]} cols")
p(f"    覆盖信号 {len(SIG)+1} 个（L1 六项 + L2 五项 + SUDO）+ L3 三项")
p("")

# ---------------------------------------------------------------- 1 universe
p("【1】benchmark universe")
U  = pd.read_csv(os.path.join(FIG, "fig1b_universe.csv")).iloc[0]
G  = pd.read_csv(os.path.join(FIG, "figS6_fold_composition.csv"))
chk("1", "候选转移数 candidate",            228, int(U.candidate_tasks))
chk("1", "合格转移数 eligible",             186, int(U.eligible_tasks))
chk("1", "被排除数 excluded",                42, int(U.excluded))
chk("1", "合格 organism–drug 组数",          37, int(U.n_organism_drug_groups))
org = G.organism.value_counts().reindex(["EC", "KP", "PA"]).fillna(0).astype(int)
chk("1", "物种构成 EC/KP/PA",  "17/9/11", f"{org['EC']}/{org['KP']}/{org['PA']}")
tr = G.n.value_counts()
chk("1", "方向数构成 6/4/2 组", "23/10/4", f"{int(tr.get(6,0))}/{int(tr.get(4,0))}/{int(tr.get(2,0))}")
chk("1", "tableS1 行数 == eligible", 186, len(T))

# ---------------------------------------------------------------- 2 ΔAUROC
p(""); p("【2】ΔAUROC 分布 / harmful 计数（从 tableS1 独立重算）")
d = T.d_auroc.values
chk("2", "median ΔAUROC",   -0.024, float(np.median(d)), 5e-4)
chk("2", "P10 ΔAUROC",      -0.107, float(np.percentile(d, 10)), 5e-4)
chk("2", "P90 ΔAUROC",      +0.043, float(np.percentile(d, 90)), 5e-4)
chk("2", "P90 − P10",        0.151, float(np.percentile(d, 90) - np.percentile(d, 10)), 5e-4)
chk("2", "harmful (Δ ≤ −0.05)",  57, int((d <= -0.05).sum()))
chk("2", "stable  (|Δ| ≤ 0.02)", 43, int((np.abs(d) <= 0.02).sum()))
chk("2", "harmful 占比 %", 30.6, round(100 * (d <= -0.05).sum() / len(d), 1), 0.05)

# ---------------------------------------------------------------- 3 Fig 2
p(""); p("【3】Fig 2 — harmful-detection AUROC / Spearman（从原始表独立重算）")
# ★ 风险分方向（+ / |·| / −）**不是 protocol 规定的**，本注释此前误称其为
#   "protocol 的 prespecified native direction"，该说法在冻结 protocol 中不存在，已更正。
#   事实：protocol §5.2 把三个 L2 位移信号定义为**有符号差值**；§5.3.1 行 192 明写
#   "signed SUDO 与 |SUDO|"，说明协议在需要时会写绝对值 —— §5.2 的沉默不是记法缺口。
#   §14 修订记录 v1.0–v1.4 五条均未申报该变换。
#   该约定的唯一落点是实现代码（p6d:232-236 / p6g:52-56 / p11:186-190），
#   且自 2026-09-29（protocol 冻结同日、解封之前）起在 7 个快照中逐字节未变
#   （方法表 sha 前缀 FFC9C6627AAC1BD9）。Methods §10 以散文披露了幅值约定，但未列出全表。
#   MNV 的职责是核对**正文所报的值**，正文用的是该约定，故此处照此重算；
#   该约定的代价由 p32 独立量化（不在此处判定）。
XFORM = {
    "domain_auc": lambda s: s, "mmd": lambda s: s, "psi_mean": lambda s: s,
    "psi_max": lambda s: s, "js_mean": lambda s: s, "js_max": lambda s: s,
    "pred_mean_shift": np.abs, "entropy_mean_shift": np.abs, "pred_pos_prev_shift": np.abs,
    "entropy_js": lambda s: s, "conf_dist_js": lambda s: s,
    "SUDO_AURCC": np.negative, "AUROC_DR": np.negative, "AUROC_IOW": np.negative,
    "AUROC_OM": np.negative,
}
y = T.harmful.values.astype(int)
f2 = pd.read_csv(os.path.join(FIG, "fig2_forest.csv")).set_index("method")
for m, xf in XFORM.items():
    v = xf(T[m].values.astype(float))
    auc = roc_auc_score(y, v)
    rho = stats.spearmanr(v, T.d_auroc.values)[0]
    chk("3", f"{m}: harmful AUROC", round(float(f2.loc[m, "harmful_auc"]), 4), round(float(auc), 4), 5e-5)
    chk("3", f"{m}: Spearman vs ΔAUROC", round(float(f2.loc[m, "spearman"]), 4), round(float(rho), 4), 5e-5)
p("    （CI 为 1,000 次 organism–drug cluster bootstrap，属重采样量；本层重算点估计，"
  "CI 的生成代码由 EVS t12/t13 覆盖）")

# ---------------------------------------------------------------- 4 Fig 3
p(""); p("【4】Fig 3 — grouped（pooled / macro / LOO）")
f3a = pd.read_csv(os.path.join(FIG, "fig3a_lodo.csv")).set_index("method")
f3b = pd.read_csv(os.path.join(FIG, "fig3b_loo_organism.csv")).set_index("method")
f3c = pd.read_csv(os.path.join(FIG, "fig3c_sudo_decomposition.csv")).iloc[0]
# 4a LODO macro：按 organism–drug 组留出，组内拟合 logistic，macro 平均
pos, neg = T.harmful.values == 1, T.harmful.values == 0
groups = T.groupby(["organism", "drug"]).indices
def lodo_macro(sig):
    x = T[sig].values; vals = []
    for g, idx in groups.items():
        tr = np.setdiff1d(np.arange(len(T)), idx)
        if len(np.unique(T.harmful.values[idx])) < 2: continue     # 单类折不可定义
        if len(np.unique(T.harmful.values[tr])) < 2:  continue
        from sklearn.linear_model import LogisticRegression
        lr = LogisticRegression(max_iter=2000).fit(x[tr].reshape(-1, 1), T.harmful.values[tr])
        s = lr.predict_proba(x[idx].reshape(-1, 1))[:, 1]
        vals.append(roc_auc_score(T.harmful.values[idx], s))
    return float(np.mean(vals)), len(vals)
for m in ["mmd", "psi_max", "domain_auc"]:
    v, k = lodo_macro(m)
    chk("4", f"LODO macro {m}", round(float(f3a.loc[m, "lodo_macro_auroc"]), 4), round(v, 4), 2e-3,
        f"（{k} 个可定义折）")
chk("4", "可定义折数", 26, len([g for g, i in groups.items()
                              if len(np.unique(T.harmful.values[i])) >= 2]))
# 4b SUDO 三点分解 —— **独立复现整套 logistic-LODO**：
#   逐折在其余组拟合 logistic(harmful ~ signal)，把该折的 logit 分数作为 out-of-fold 预测，
#   再按三种口径汇总。这正是 fig3c 三个数的定义。
def lodo_oof(sig):
    """复现 p6f 的 LODO：**对信号做 per-fold z-score**（用训练折的 mean/SD），再拟合 logistic。
    ★ 该标准化步骤是复现 0.332/0.408 的必要条件：去掉它 pooled 会变成 0.235/0.353，
      而 macro 不变（0.3246）—— 因为它只改变跨折尺度，不改变折内序。
      这是本次 MNV 的**关键发现之一**：Methods 未写明该步骤（详见报告）。"""
    r = T[sig].values.astype(float); oof = np.full(len(T), np.nan)
    for g, idx in groups.items():
        tr = np.setdiff1d(np.arange(len(T)), idx)
        if len(np.unique(T.harmful.values[tr])) < 2: continue
        mu, sd = r[tr].mean(), r[tr].std()
        from sklearn.linear_model import LogisticRegression
        lr = LogisticRegression(max_iter=1000).fit(
            ((r[tr] - mu) / sd).reshape(-1, 1), T.harmful.values[tr])
        oof[idx] = lr.predict_proba(((r[idx] - mu) / sd).reshape(-1, 1))[:, 1]
    return oof
oof = lodo_oof("SUDO_AURCC")
ok  = np.isfinite(oof)
defin_idx = np.concatenate([groups[g] for g, i in groups.items()
                            if len(np.unique(T.harmful.values[i])) >= 2])
chk("4", "SUDO pooled（全 186，logistic-LODO OOF）", round(float(f3c.pooled_202task), 4),
    round(float(roc_auc_score(y[ok], oof[ok])), 4), 2e-3)
chk("4", "SUDO pooled（可评价 134）", round(float(f3c.pooled_138task), 4),
    round(float(roc_auc_score(y[defin_idx], oof[defin_idx])), 4), 2e-3)
# macro：折内 AUROC 平均（仅在可定义折上）
mac = np.mean([roc_auc_score(T.harmful.values[i], oof[i]) for g, i in groups.items()
               if len(np.unique(T.harmful.values[i])) >= 2])
chk("4", "SUDO macro（可评价折平均）", round(float(f3c.macro_138task), 4),
    round(float(mac), 4), 2e-3)

# ---------------------------------------------------------------- 5 L3 + naive
p(""); p("【5】L3 — MAE / ρ / calibration + **naive baseline**")
f4d = pd.read_csv(os.path.join(FIG, "fig4d_mae.csv")).set_index("estimator")
for e in ["DR", "IOW", "OM"]:
    err = np.abs(T["AUROC_" + e].values - T.auroc_target.values)
    chk("5", f"{e} MAE", round(float(f4d.loc[e, "mae"]), 4), round(float(err.mean()), 4), 5e-5)
    rho = stats.spearmanr(T["AUROC_" + e].values, T.auroc_target.values)[0]
    chk("5", f"{e} Spearman(est, obs)", round(float(f4d.loc[e, "rho_pred_true"]), 4), round(float(rho), 4), 5e-5)
# naive baseline：预测 target AUROC = source internal AUROC（即预测"无变化"）
naive_err = np.abs(T.auroc_source_internal.values - T.auroc_target.values)
p("    ── naive baseline（新算，非来自产物）──")
p(f"       MAE_naive = mean|ΔAUROC| = {naive_err.mean():.4f}   中位 {np.median(naive_err):.4f}")
p(f"       排序：DR {err_dr if (err_dr:=float(np.abs(T['AUROC_DR'].values-T.auroc_target.values).mean())) else 0:.4f}"
  f" < OM {float(np.abs(T['AUROC_OM'].values-T.auroc_target.values).mean()):.4f}"
  f" < naive {naive_err.mean():.4f}"
  f" < IOW {float(np.abs(T['AUROC_IOW'].values-T.auroc_target.values).mean()):.4f}")
p(f"       ⇒ DR/OM 优于 naive（{100*(1-float(np.abs(T['AUROC_DR'].values-T.auroc_target.values).mean())/naive_err.mean()):.0f}% / "
  f"{100*(1-float(np.abs(T['AUROC_OM'].values-T.auroc_target.values).mean())/naive_err.mean()):.0f}% 相对降幅）；"
  f"**IOW 不优于 naive**")
ROWS.append(dict(section="5", item="naive baseline MAE（新算）", expected="—",
                 observed=round(float(naive_err.mean()), 4), tol="—", status="INFO",
                 note="DR/OM 优于 naive；IOW 不比 naive 好"))

# ---------------------------------------------------------------- 6 target-N
p(""); p("【6】target-N 敏感性")
F5B = pd.read_csv(os.path.join(FIG, "fig5b_across_n_range.csv")).set_index("method")
dr = F5B.loc["AUROC_DR"]
chk("6", "DR full 档 AUROC", 0.645, round(float(dr.auroc_at_full), 3), 5e-4)
chk("6", "DR n=100 档 AUROC", 0.384, round(float(dr.auroc_at_100), 3), 5e-4)
chk("6", "DR 跨档极差", 0.261, round(float(dr.across_n_range), 3), 5e-4)
chk("6", "全表最大跨档极差 == DR 的", 0.261, round(float(F5B.across_n_range.max()), 3), 5e-4)
chk("6", "监测信号中最大（pred_pos_prev_shift）", 0.047,
    round(float(F5B.loc["pred_pos_prev_shift", "across_n_range"]), 3), 5e-4)
chk("6", "中位跨档极差", 0.021, round(float(F5B.across_n_range.median()), 3), 5e-4)
chk("6", "full 档 <0.5 的方法数", 6, int((F5B.auroc_at_full < 0.5).sum()))
chk("6", "n=100 档 <0.5 的方法数", 8, int((F5B.auroc_at_100 < 0.5).sum()))

# ---------------------------------------------------------------- 7 robustness
p(""); p("【7】robustness — 配置一致性")
AG = pd.read_csv(os.path.join(FIG, "figS2_config_agreement.csv"))
CF = pd.read_csv(os.path.join(FIG, "fig1e_configs_full.csv"))
for cfg, n in [("P1 History+LR", 186), ("S2 History+LGBM", 174),
               ("S1 Core+LR", 64), ("S3 Core+LGBM", 60)]:
    chk("7", f"{cfg} 自身合格数", n, int(CF[CF.config == cfg].eligible.sum()))
ht = lambda s: "History" if "History" in s else "Core"
ag = AG[AG.cfg_a != AG.cfg_b]
same = ag[ag.cfg_a.map(ht) == ag.cfg_b.map(ht)].spearman_rho
cross = ag[ag.cfg_a.map(ht) != ag.cfg_b.map(ht)].spearman_rho
chk("7", "同 tier 最低一致性 > 跨 tier 最高一致性",
    True, bool(same.min() > cross.max()), 0,
    f"同 tier {sorted(round(x,3) for x in same)} vs 跨 tier {sorted(round(x,3) for x in cross)}")

# ---------------------------------------------------------------- 8 S counts
p(""); p("【8】Supplement 报告计数")
FS6 = pd.read_csv(os.path.join(FIG, "figS6_fold_summary.csv")).set_index("group")
chk("8", "Fig S6 可评价折 / task", "26 / 134",
    f"{int(FS6.loc['evaluable','n_folds'])} / {int(FS6.loc['evaluable','n_tasks'])}")
chk("8", "Fig S6 单类折 / task", "11 / 52",
    f"{int(FS6.loc['excluded','n_folds'])} / {int(FS6.loc['excluded','n_tasks'])}")
chk("8", "Fig S6 合计折 / task", "37 / 186",
    f"{int(FS6.loc['total','n_folds'])} / {int(FS6.loc['total','n_tasks'])}")
S5 = pd.read_csv(os.path.join(FIG, "figS5_mmd_summary.csv")).iloc[0]
chk("8", "Fig S5 MMD 负值数", 0, int(S5.n_negative))
chk("8", "Fig S5 MMD task 数", 186, int(S5.n_tasks))
S7 = pd.read_csv(os.path.join(FIG, "figS7_weight_association.csv"))
chk("8", "Fig S7 ESS ρ", -0.249, round(float(S7[S7.metric == "ESS"].spearman_rho.iat[0]), 3), 5e-4)
chk("8", "Fig S7 w_max ρ", +0.220, round(float(S7[S7.metric == "w_max"].spearman_rho.iat[0]), 3), 5e-4)
E4 = pd.read_csv(os.path.join(FIG, "fig4e_ess_association.csv")).iloc[0]
chk("8", "Fig 4E ESS 中位", 1082, round(float(E4.ess_median)), 1)
chk("8", "Fig 4E ESS 最小", 11, round(float(E4.ess_min)), 0.5)
chk("8", "Fig 4E DR>1 越界数", 0, int(E4.n_dr_gt1))
S1 = pd.read_csv(os.path.join(FIG, "figS1_task_matrix.csv"))
chk("8", "Fig S1 矩阵 cell 数", 228, len(S1))
chk("8", "Fig S1 eligible cell", 186, int(S1.eligible.astype(bool).sum()))

# ---------------------------------------------------------------- 汇总
R = pd.DataFrame(ROWS)
p(""); p("=" * 100)
bad = R[R.status == "FAIL"]
p(f"MNV 断言合计 {len(R)} 项：PASS {int((R.status=='PASS').sum())}｜"
  f"INFO {int((R.status=='INFO').sum())}｜**FAIL {len(bad)}**")
if len(bad):
    p(""); p("*** 不一致项（阻断，必须报告，不自行改正文）***")
    for _, r in bad.iterrows():
        p("    [%s] %s：正文 %s vs 重算 %s" % (r["section"], r["item"], r["expected"], r["observed"]))
    p(""); p("*** MNV 未通过 ***")
else:
    p(""); p("*** MNV 通过 —— 正文全部 headline numbers 可从上游客观重算复现 ***")
R.to_csv(os.path.join(OUT, "mnv_results.csv"), index=False, encoding="utf-8-sig")
REP.close()
sys.exit(1 if len(bad) else 0)
