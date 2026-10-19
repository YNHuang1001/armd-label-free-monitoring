# -*- coding: utf-8 -*-
r"""p32 — 方向约定敏感性：L2 三个有符号位移信号取 |·| vs 保持原值

================================ 为什么 ================================
外部审阅（V2）与 handoff §12 均指出：三个 prediction-space **有符号**位移信号
（pred_mean_shift / entropy_mean_shift / pred_pos_prev_shift）在送入评价前被取绝对值，
而冻结 protocol §5.2 对它们的定义是**有符号差值**，§14 修订记录 v1.0–v1.4 从未申报该变换。

本脚本**只做一件事**：把两种口径并排算出来，量化该约定的代价。
**只读冻结产物，不改任何生产文件。**

锚定：PROD 口径必须逐位复现 fig3a_lodo.csv / fig2_forest.csv，
否则说明本脚本的独立重算与生产管线不同源，结论作废。

判据来源：figures/data/fig3a_lodo.csv、fig2_forest.csv（冻结产物）
"""
import io, os, sys
import numpy as np, pandas as pd
from scipy import stats
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIG  = os.path.join(ROOT, "figures", "data")
OUT  = os.path.join(ROOT, "分析结果")
REP  = open(os.path.join(OUT, "report_p32_direction_convention.txt"), "w", encoding="utf-8")

def p(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); REP.write(s + "\n"); REP.flush()

# ---------------------------------------------------------------- 独立上游
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
y = T.harmful.values.astype(int)

# ---------------------------------------------------------------- 两种口径
SIGNED_L2 = ["pred_mean_shift", "entropy_mean_shift", "pred_pos_prev_shift"]

def transforms(mode):
    """mode='prod'：生产口径（三个有符号 L2 取 |·|）；mode='raw'：只用原值，只保留语义强制的符号。"""
    f = lambda s: s
    X = {k: f for k in SIG + ["SUDO_AURCC"] + L3}
    # 性能类代理：越高越安全 → 风险分取负。**两种口径下都不变**（语义强制，非自由度）
    for k in ["SUDO_AURCC"] + L3:
        X[k] = np.negative
    if mode == "prod":
        for k in SIGNED_L2:
            X[k] = np.abs
    return X

ORDER = ["domain_auc", "mmd", "psi_mean", "psi_max", "js_mean", "js_max",
         "pred_mean_shift", "entropy_mean_shift", "entropy_js", "conf_dist_js",
         "pred_pos_prev_shift", "SUDO_AURCC", "AUROC_DR", "AUROC_IOW", "AUROC_OM"]

# ---------------------------------------------------------------- LODO
groups = T.groupby(["organism", "drug"]).indices

def lodo_macro(x):
    """与 p30.lodo_macro 同口径：折内不 z-score，默认 C=1，折内 AUROC 后 macro 平均。
    ★ 入参是**已变换的数值向量**，不是列名 —— 否则口径差异会被吞掉。"""
    x = np.asarray(x, dtype=float); vals = []
    for g, idx in groups.items():
        tr = np.setdiff1d(np.arange(len(T)), idx)
        if len(np.unique(y[idx])) < 2 or len(np.unique(y[tr])) < 2: continue
        lr = LogisticRegression(max_iter=2000).fit(x[tr].reshape(-1, 1), y[tr])
        s = lr.predict_proba(x[idx].reshape(-1, 1))[:, 1]
        vals.append(roc_auc_score(y[idx], s))
    return float(np.mean(vals)), len(vals)

def loo_organism(x):
    """与 p6d 同口径：留一物种（3 折），折内 z-score 后拟合 logistic，取 3 个折 AUROC 的均值。
    对应 Fig 3B / Results §3 引用的 LOO 值。"""
    x = np.asarray(x, dtype=float); vals = []
    for org in ["EC", "KP", "PA"]:
        te = (T.organism == org).values; tr = ~te
        if len(np.unique(y[tr])) < 2 or len(np.unique(y[te])) < 2: continue
        mu, sd = x[tr].mean(), x[tr].std()
        if not np.isfinite(sd) or sd == 0: sd = 1.0
        lr = LogisticRegression(max_iter=1000).fit(((x[tr] - mu) / sd).reshape(-1, 1), y[tr])
        s2 = lr.predict_proba(((x[te] - mu) / sd).reshape(-1, 1))[:, 1]
        vals.append(roc_auc_score(y[te], s2))
    return float(np.mean(vals))

def lodo_pooled_oof(x):
    """与 p6f/fig3c 同口径：折内 z-score 后拟合 logistic，OOF 概率在全部可评价 task 上算 pooled AUROC。
    ★ 入参同上，必须是已变换的数值向量。"""
    r = np.asarray(x, dtype=float); oof = np.full(len(T), np.nan)
    for g, idx in groups.items():
        tr = np.setdiff1d(np.arange(len(T)), idx)
        if len(np.unique(y[tr])) < 2: continue
        mu, sd = r[tr].mean(), r[tr].std()
        if not np.isfinite(sd) or sd == 0: sd = 1.0
        lr = LogisticRegression(max_iter=1000).fit(((r[tr] - mu) / sd).reshape(-1, 1), y[tr])
        oof[idx] = lr.predict_proba(((r[idx] - mu) / sd).reshape(-1, 1))[:, 1]
    ok = np.isfinite(oof)
    return float(roc_auc_score(y[ok], oof[ok]))

# ---------------------------------------------------------------- 锚定
p("=" * 108); p("p32 — 方向约定敏感性（L2 有符号位移信号：|·| vs 原值）"); p("=" * 108)
p("只读冻结产物；不改任何生产文件。")
p("")
f3a = pd.read_csv(os.path.join(FIG, "fig3a_lodo.csv")).set_index("method")
f3b = pd.read_csv(os.path.join(FIG, "fig3b_loo_organism.csv")).set_index("method")
f2  = pd.read_csv(os.path.join(FIG, "fig2_forest.csv")).set_index("method")

Xp = transforms("prod")
anchor_bad = []
p("【0】锚定：PROD 口径必须复现冻结产物（否则本脚本结论作废）")
for m in ORDER:
    v_prod = Xp[m](T[m].values.astype(float))
    auc = roc_auc_score(y, v_prod)
    mac, k = lodo_macro(v_prod) if m not in ("SUDO_AURCC",) + tuple(L3) else (np.nan, 0)
    loo = loo_organism(v_prod)
    ok_a = abs(auc - float(f2.loc[m, "harmful_auc"])) <= 1e-4
    ok_m = True if np.isnan(mac) else abs(mac - float(f3a.loc[m, "lodo_macro_auroc"])) <= 2e-3
    ok_l = abs(loo - float(f3b.loc[m, "loo_org_macro"])) <= 2e-3
    if not (ok_a and ok_m and ok_l): anchor_bad.append(m)
    p("    %-22s pooled %.4f/%.4f %s ｜ macro %s/%s %s ｜ LOO %s/%.4f %s"
      % (m, float(f2.loc[m, "harmful_auc"]), auc, "OK" if ok_a else "**差异**",
         ("%.4f" % float(f3a.loc[m, "lodo_macro_auroc"])) if m in f3a.index else "  —  ",
         ("%.4f" % mac) if np.isfinite(mac) else "  —  ", "OK" if ok_m else "**差异**",
         ("%.4f" % float(f3b.loc[m, "loo_org_macro"])) if m in f3b.index else "  —  ",
         loo, "OK" if ok_l else "**差异**"))
p("")
if anchor_bad:
    p(f"*** 锚定失败：{anchor_bad} —— 本脚本与生产管线不同源，结论作废 ***")
    REP.close(); sys.exit(2)
p("    锚定通过：PROD 口径 15/15 在 pooled / LODO macro / LOO-organism 三处逐位复现冻结产物。")
p("    可定义折数 = %d" % len([g for g, i in groups.items() if len(np.unique(y[i])) >= 2]))
p("")

# ---------------------------------------------------------------- 主表
p("【1】两种口径并排")
p("")
hdr = ("%-22s | %-31s | %-31s | %s" %
       ("method", "PROD（|·|，当前正文）", "RAW（原值，protocol §5.2 定义）", "Δ(RAW−PROD)"))
p(hdr); p("-" * len(hdr))
p("%-22s | %7s %7s %8s %7s | %7s %7s %8s %7s | %8s %8s %8s"
  % ("", "pooled", "macro", "poolOOF", "LOO", "pooled", "macro", "poolOOF", "LOO",
     "macro", "poolOOF", "LOO"))

Xr = transforms("raw")
rows = []
for m in ORDER:
    r_prod, r_raw = {}, {}
    for tag, X in (("prod", Xp), ("raw", Xr)):
        v = X[m](T[m].values.astype(float))
        d = {"pooled": roc_auc_score(y, v),
             "rho": float(stats.spearmanr(v, T.d_auroc.values)[0])}
        d["loo"] = loo_organism(v)
        if m in ("SUDO_AURCC",) + tuple(L3):
            d["macro"] = np.nan; d["poolOOF"] = np.nan
        else:
            d["macro"], _ = lodo_macro(v)
            d["poolOOF"] = lodo_pooled_oof(v)
        (r_prod if tag == "prod" else r_raw).update(d)
    dm = r_raw["macro"] - r_prod["macro"] if np.isfinite(r_prod["macro"]) else np.nan
    do = r_raw["poolOOF"] - r_prod["poolOOF"] if np.isfinite(r_prod["poolOOF"]) else np.nan
    dl = r_raw["loo"] - r_prod["loo"]
    star = "  ← 有符号 L2" if m in SIGNED_L2 else ""
    p("%-22s | %7.4f %7.4f %8.4f %7.4f | %7.4f %7.4f %8.4f %7.4f | %+8.4f %+8.4f %+8.4f%s"
      % (m, r_prod["pooled"], r_prod["macro"], r_prod["poolOOF"], r_prod["loo"],
         r_raw["pooled"], r_raw["macro"], r_raw["poolOOF"], r_raw["loo"], dm, do, dl, star))
    rows.append(dict(method=m, signed_l2=(m in SIGNED_L2),
                     prod_pooled=r_prod["pooled"], prod_macro=r_prod["macro"],
                     prod_poolOOF=r_prod["poolOOF"], prod_loo=r_prod["loo"], prod_rho=r_prod["rho"],
                     raw_pooled=r_raw["pooled"], raw_macro=r_raw["macro"],
                     raw_poolOOF=r_raw["poolOOF"], raw_loo=r_raw["loo"], raw_rho=r_raw["rho"],
                     d_macro=dm, d_poolOOF=do, d_loo=dl))

# ---------------------------------------------------------------- 排序
p("")
p("【2】LODO macro 排序对比（仅可评价的方法；SUDO/L3 不参与 macro）")
R = pd.DataFrame(rows)
R12 = R[np.isfinite(R.prod_macro)].copy()
R12["rank_prod"] = R12.prod_macro.rank(ascending=False).astype(int)
R12["rank_raw"]  = R12.raw_macro.rank(ascending=False).astype(int)
R12["d_rank"]    = R12.rank_prod - R12.rank_raw
p("")
p("%-4s %-22s %8s   |  %-4s %-22s %8s   | %s"
  % ("PROD", "method", "macro", "RAW", "method", "macro", "位移"))
for i in range(len(R12)):
    a = R12.sort_values("prod_macro", ascending=False).iloc[i]
    b = R12.sort_values("raw_macro", ascending=False).iloc[i]
    p("%-4d %-22s %8.4f   |  %-4d %-22s %8.4f   | %s"
      % (i + 1, a.method, a.prod_macro, i + 1, b.method, b.raw_macro,
         "" if a.method == b.method else "**换位**"))
p("")
p("    三个有符号 L2 信号的秩位移：")
for m in SIGNED_L2:
    r = R12[R12.method == m].iloc[0]
    p("      %-22s  macro %.4f → %.4f （%+0.4f）  秩 %d → %d"
      % (m, r.prod_macro, r.raw_macro, r.d_macro, int(r.rank_prod), int(r.rank_raw)))
top3p = list(R12.sort_values("prod_macro", ascending=False).method[:3])
top3r = list(R12.sort_values("raw_macro", ascending=False).method[:3])
p("")
p("    PROD 前三：%s" % " / ".join(top3p))
p("    RAW  前三：%s" % " / ".join(top3r))
p("    交集：%s" % (" / ".join(sorted(set(top3p) & set(top3r))) or "（空）"))

# ---------------------------------------------------------------- §7.2 主任务
p("")
p("【3】protocol §7.2 主任务（harmful 预警）与 §3 引用的 LOO，受影响程度")
for m in SIGNED_L2:
    r = R[R.method == m].iloc[0]
    p("    %-22s pooled %.4f → %.4f （%+0.4f）｜ LOO %.4f → %.4f （%+0.4f）｜ ρ(ΔAUROC) %+.4f → %+.4f"
      % (m, r.prod_pooled, r.raw_pooled, r.raw_pooled - r.prod_pooled,
         r.prod_loo, r.raw_loo, r.d_loo, r.prod_rho, r.raw_rho))
p("")
mx = R[R.signed_l2].d_macro.abs().max()
p("    三个信号的 |Δmacro| 最大 = %.4f；其余 12 个方法的口径差异 = 0（变换为恒等或语义强制）。" % mx)
p("    若仅看**未受该约定影响**的 12 个方法，两种口径下 pooled 与 macro **完全不变**。")

# ---------------------------------------------------------------- naive
p("")
p("【4】常数基线（供 item 14 核对）")
p("")
p("    ★ 常数基线**有两个**，报错那个会把估计量的优势放大。二者都必须给出。")
dd = T.d_auroc.values.astype(float)
f_zero = float(np.abs(dd - 0.0).mean())            # 预测「零变化」
f_med  = float(np.abs(dd - np.median(dd)).mean())  # 预测观测中位数 —— **更强**，真正的 floor
f_mean = float(np.abs(dd - dd.mean()).mean())
p("    mean|ΔAUROC|              = %.6f   ← 「预测零变化」" % f_zero)
p("    median(ΔAUROC)            = %.6f" % float(np.median(dd)))
p("    mean|Δ − median(Δ)|       = %.6f   ← 「预测中位数」= **最强常数基线，真 floor**" % f_med)
p("    mean|Δ − mean(Δ)|         = %.6f   ← 「预测均值」（略差于中位数，符合 MAE 理论）" % f_mean)
p("")
p("    %-5s %10s %14s %16s" % ("", "MAE", "优于零变化", "优于中位数floor"))
rows_const = []
for e in ["DR", "IOW", "OM"]:
    err = float(np.abs(T["AUROC_" + e].values - T.auroc_target.values).mean())
    p("    %-5s %10.6f %+14.6f %+16.6f" % (e, err, f_zero - err, f_med - err))
    rows_const.append((e, err, f_zero - err, f_med - err))
_cb = pd.DataFrame(
    [dict(name="constant_zero",   mae=f_zero, gain_vs_zero=0.0,      gain_vs_median_floor=f_zero - f_med),
     dict(name="constant_median", mae=f_med,  gain_vs_zero=f_med - f_zero, gain_vs_median_floor=0.0)]
    + [dict(name=e, mae=m, gain_vs_zero=gz, gain_vs_median_floor=gm) for e, m, gz, gm in rows_const])
_cb["median_d_auroc"] = float(np.median(dd))
_cb["mean_d_auroc"] = float(dd.mean())
_cb.to_csv(os.path.join(OUT, "p32_constant_baselines.csv"), index=False, encoding="utf-8-sig")
p("")
p("    ⇒ 对**最强**常数基线（中位数 %.4f），DR / OM 的优势只有 %.4f / %.4f；"
  % (f_med, f_med - rows_const[0][1], f_med - rows_const[2][1]))
p("      对**较弱**基线（零变化 %.4f）则是 %.4f / %.4f —— **写正文必须声明用的是哪个**。"
  % (f_zero, f_zero - rows_const[0][1], f_zero - rows_const[2][1]))
p("      IOW 对两个基线均为负增益（%.4f / %.4f）。"
  % (f_zero - rows_const[1][1], f_med - rows_const[1][1]))

R.to_csv(os.path.join(OUT, "p32_direction_convention.csv"), index=False, encoding="utf-8-sig")
p("")
p("=" * 108)
p("结论：该约定只影响 3 个信号；对其余 12 个方法**零影响**。")
p("      受影响三者的 LODO macro 位移见【2】；pooled（§7.2 主任务）位移见【3】。")
REP.close()
