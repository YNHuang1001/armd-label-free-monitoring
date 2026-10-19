# -*- coding: utf-8 -*-
r"""L1 direction-learning sensitivity（v2）—— post-unblinding explanatory analysis

**为什么做**：把 L1 的 LODO AUROC < 0.5 直接解读成"X-space 指标无信息"不够严谨。
必须在三种互斥解释中区分：
  (a) 无关联                        → 所有口径都 ≈ 0.5
  (b) 关联方向与惯例相反、且**方向可迁移** → 学到反向符号后应 > 0.5
  (c) 方向**不可迁移**（组间不一致 / 组内样本太小 / 非线性）→ 学符号后仍 ≈ 0.5 或更低

**v1 的错误（已修正）**：v1 用 Pearson 相关的符号作为"训练集符号"诊断。
但 L1 的线性相关本质为零（mmd: point-biserial r = −0.0040, p = 0.955），
其符号是噪声，会造出"符号不稳定"的假象。改用与 AUROC 同口径的**秩方向**。

**本版三个口径**：
  ① macro_within_raw  : 原始 signal 在每个组内直接算 AUROC 再平均（不学任何东西；惯例方向）
  ② macro_learned_logit: 每折在其余组上拟合 logistic harmful~signal，符号+斜率全由训练组学得
  ③ macro_learned_rank : 每折用训练组的 AUROC 定符号（sign = sign(AUC_train − 0.5)），
                        再套到留出组 —— 与 AUROC 同口径、无线性假设

**定性**：看完全部解封结果后设计，属 post-unblinding explanatory analysis，不得冒充预注册 primary。

输出：report_p6g_direction.txt / p6g_direction.csv / p6g_fold_signs.csv
"""
import os
import numpy as np, pandas as pd
from sklearn.linear_model import LogisticRegression
from scipy.stats import spearmanr

BASE = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))); OUT = os.path.dirname(os.path.abspath(__file__))
SEED = 20261005; rng = np.random.default_rng(SEED); B = 4000

rep = open(os.path.join(OUT, "report_p6g_direction.txt"), "w", encoding="utf-8")
def p(*a):
    line = " ".join(str(x) for x in a); print(line, flush=True); rep.write(line + "\n"); rep.flush()

O = pd.read_csv(os.path.join(OUT, "p6d_outcome.csv"))
SIG = pd.read_csv(os.path.join(OUT, "p6b_label_free_signals_blinded.csv"))
L3 = pd.read_csv(os.path.join(OUT, "p6c_L3_blinded.csv"))
E = O.merge(SIG[["task_id", "domain_auc", "mmd", "psi_mean", "psi_max", "js_mean", "js_max",
                 "pred_mean_shift", "entropy_mean_shift", "entropy_js", "conf_dist_js",
                 "pred_pos_prev_shift"]], on="task_id", how="left") \
     .merge(L3[["task_id", "SUDO_AURCC", "AUROC_IOW", "AUROC_OM", "AUROC_DR"]], on="task_id", how="left")
E["group"] = E.organism + "|" + E.drug

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

y = E["harmful"].values.astype(int)
g = E["group"].values
GROUPS = sorted(np.unique(g))

def auc(yy, ss):
    m = np.isfinite(ss)
    if len(np.unique(yy[m])) < 2: return np.nan
    yy, ss = yy[m], ss[m]
    o = np.argsort(ss, kind="mergesort"); yy = yy[o]; ss = ss[o]
    r = np.empty(len(ss)); i = 0
    while i < len(ss):
        j = i
        while j + 1 < len(ss) and ss[j + 1] == ss[i]: j += 1
        r[i:j + 1] = 0.5 * (i + j) + 1
        i = j + 1
    n1 = yy.sum(); n0 = len(yy) - n1
    if n1 == 0 or n0 == 0: return np.nan
    return (r[yy == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)

def boot(vals_map):
    ks = list(vals_map.keys()); out = []
    for _ in range(B):
        pick = rng.choice(len(ks), len(ks), replace=True)
        out.append(np.mean([vals_map[ks[i]] for i in pick]))
    return np.percentile(out, [2.5, 97.5])

p("=" * 122)
p("L1 direction-learning sensitivity（v2）—— 区分「无信息」/「可迁移的反向」/「方向不可迁移」")
p("=" * 122)
p("  ① macro_within_raw   : 组内原始方向（惯例：越大越可能 failure），不学任何东西")
p("  ② macro_learned_logit: 每折在其余组拟合 logistic harmful~signal，符号+斜率全由训练组学")
p("  ③ macro_learned_rank : 每折用训练组 AUROC 定符号（与 AUROC 同口径，无线性假设）")
p("  注：post-unblinding explanatory analysis，非预注册 primary。")
p("  注：v1 曾用 Pearson 相关定符号，因 L1 线性相关≈0 而产生假象，已弃用。")

rows, fold_rows = [], []
for name, role, fn in METHODS:
    r = fn(E).values.astype(float)
    oof = np.full(len(r), np.nan)
    within, rank_dir, logit_sign, sp_sign, tr_auc = {}, {}, {}, {}, {}
    for gg in GROUPS:
        te = g == gg; tr = ~te
        if te.sum() == 0: continue
        a_w = auc(y[te], r[te])
        if np.isfinite(a_w): within[gg] = a_w
        if len(np.unique(y[tr])) == 2:
            a_tr = auc(y[tr], r[tr])
            if np.isfinite(a_tr):
                tr_auc[gg] = a_tr
                rank_dir[gg] = 1.0 if a_tr > 0.5 else (-1.0 if a_tr < 0.5 else 0.0)
            mu, sd = np.nanmean(r[tr]), np.nanstd(r[tr]) or 1.0
            try:
                m = LogisticRegression(max_iter=1000).fit(((r[tr] - mu) / sd).reshape(-1, 1), y[tr])
                oof[te] = m.predict_proba(((r[te] - mu) / sd).reshape(-1, 1))[:, 1]
                logit_sign[gg] = float(np.sign(m.coef_[0, 0]))
            except Exception:
                pass
            if np.nanstd(r[tr]) > 0 and len(np.unique(y[tr])) == 2:
                sp_sign[gg] = float(np.sign(spearmanr(r[tr], y[tr]).statistic))
        fold_rows.append(dict(method=name, group=gg, n_fold=int(te.sum()),
                              n_pos=int(y[te].sum()), n_neg=int((1 - y[te]).sum()),
                              within_auc=a_w, train_auc=tr_auc.get(gg, np.nan),
                              rank_dir=rank_dir.get(gg, np.nan),
                              logit_sign=logit_sign.get(gg, np.nan),
                              spearman_sign=sp_sign.get(gg, np.nan),
                              oof_auc=auc(y[te], oof[te])))

    ev = sorted(within.keys())
    mac_within = float(np.mean([within[k] for k in ev]))
    lo_w, hi_w = boot(within)
    # ③ rank-direction learning：留出组得分 = rank_dir * r，AUC = AUC if dir>0 else 1-AUC
    rk = {}
    for k in ev:
        if k not in rank_dir or not np.isfinite(within[k]): continue
        rk[k] = within[k] if rank_dir[k] > 0 else (1 - within[k]) if rank_dir[k] < 0 else np.nan
    rk = {k: v for k, v in rk.items() if np.isfinite(v)}
    mac_rank = float(np.mean(list(rk.values()))) if rk else np.nan
    lo_r, hi_r = boot(rk) if rk else (np.nan, np.nan)
    # ② logistic
    oo = {k: auc(y[g == k], oof[g == k]) for k in ev}
    oo = {k: v for k, v in oo.items() if np.isfinite(v)}
    mac_logit = float(np.mean(list(oo.values()))) if oo else np.nan
    lo_l, hi_l = boot(oo) if oo else (np.nan, np.nan)
    # 方向一致率
    agree = np.mean([rank_dir[k] == logit_sign.get(k, np.nan) for k in ev if k in rank_dir and k in logit_sign]) \
        if any(k in rank_dir and k in logit_sign for k in ev) else np.nan
    rows.append(dict(method=name, role=role, n_eval_fold=len(ev),
                     macro_within=mac_within, within_lo=lo_w, within_hi=hi_w,
                     macro_logit=mac_logit, logit_lo=lo_l, logit_hi=hi_l,
                     macro_rank=mac_rank, rank_lo=lo_r, rank_hi=hi_r,
                     pooled_raw_all202=auc(y, r),
                     frac_train_dir_pos=float(np.mean([rank_dir[k] > 0 for k in ev if k in rank_dir])),
                     logit_vs_rank_agreement=agree))
D = pd.DataFrame(rows)
D.to_csv(os.path.join(OUT, "p6g_direction.csv"), index=False)
pd.DataFrame(fold_rows).to_csv(os.path.join(OUT, "p6g_fold_signs.csv"), index=False)

p("\n" + "=" * 122)
p(f"{'method':<22}{'role':<11}{'k':>4}{'① within_raw':>16}{'② learn_logit':>16}"
  f"{'③ learn_rank':>16}{'pooled(all)':>11}{'训练方向为正':>12}")
for _, r_ in D.iterrows():
    p(f"{r_.method:<22}{r_.role:<11}{r_.n_eval_fold:>4}"
      f"{r_.macro_within:>8.3f}[{r_.within_lo:.2f},{r_.within_hi:.2f}]"
      f"{r_.macro_logit:>8.3f}[{r_.logit_lo:.2f},{r_.logit_hi:.2f}]"
      f"{r_.macro_rank:>8.3f}[{r_.rank_lo:.2f},{r_.rank_hi:.2f}]"
      f"{r_.pooled_raw_all202:>11.3f}{r_.frac_train_dir_pos:>12.2f}")

p("\n" + "=" * 122)
p("【判读】")
p("=" * 122)
p("  (a) 无信息          : ① 与 ③ 的 CI 都含 0.5")
p("  (b) 可迁移·惯例方向 : 训练方向 100% 为正 且 ③ 的 CI 下界 > 0.5")
p("  (b') 可迁移·反向    : 训练方向 100% 为负 且 ③ 的 CI 下界 > 0.5 ← GPT6 预判的情形")
p("  (c) 方向不可迁移    : ③ CI 含 0.5，但训练方向在折间不一致（0 < 比例 < 1）")
p("  (d) 学方向反而更差  : ③ 的 CI 上界 < 0.5")
p("")
for _, r_ in D.iterrows():
    rk_lo, rk_hi = r_.rank_lo, r_.rank_hi
    f = r_.frac_train_dir_pos
    if rk_hi < 0.5:
        v = "(d) 学方向反而更差 —— ③ 显著 < 0.5"
    elif rk_lo > 0.5 and f >= 0.99:
        v = "(b) 可迁移·惯例方向"
    elif rk_lo > 0.5 and f <= 0.01:
        v = "(b') 可迁移·反向 ← 方向可学但与惯例相反"
    elif rk_lo > 0.5:
        v = "③ CI 下界 > 0.5 但训练方向不稳定，需谨慎"
    elif 0.01 < f < 0.99:
        v = "(c) 方向不可迁移 —— 训练方向折间不一致"
    else:
        v = "(a) 无信息 —— 全程 ≈ 0.5"
    p(f"  {r_.method:<22}{r_.role:<11} ①[{r_.within_lo:.3f},{r_.within_hi:.3f}] "
      f"③[{rk_lo:.3f},{rk_hi:.3f}]  训练方向为正比例 {f:.2f}  → {v}")

p("\n" + "=" * 122)
p("【L1 专项：对主结论第二条的直接影响】")
p("=" * 122)
p(f"{'method':<16}{'线性相关(点二列)':>17}{'秩相关(Spearman)':>18}{'①组内':>9}{'③学方向':>10}{'惯例方向可用?':>16}")
for _, r_ in D[D.role == "L1"].iterrows():
    r_vals = fn = None
    rr = METHODS[[m[0] for m in METHODS].index(r_.method)][2](E).values.astype(float)
    from scipy.stats import pointbiserialr
    pb = pointbiserialr(y, rr).statistic
    sp = spearmanr(rr, y).statistic
    verdict = "可用" if r_.rank_lo > 0.5 else ("不可用（方向学不到）" if r_.rank_hi >= 0.5 else "不可用（方向相反）")
    p(f"{r_.method:<16}{pb:>+17.4f}{sp:>+18.4f}{r_.macro_within:>9.3f}{r_.macro_rank:>10.3f}{verdict:>16}")

p("\n  事后全局翻转（in-sample，部署不可得，仅作对照）:")
for _, r_ in D[D.role == "L1"].iterrows():
    p(f"    {r_.method:<16} 1 - ①组内 = {1 - r_.macro_within:.3f}   "
      f"1 - pooled(all) = {1 - r_.pooled_raw_all202:.3f}")
rep.close()
print("done")
