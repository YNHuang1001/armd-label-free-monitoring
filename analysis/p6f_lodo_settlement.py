# -*- coding: utf-8 -*-
r"""LODO aggregation ambiguity —— 彻底结清

背景：protocol v1.3 §8.1.1 记录了一个未预注册的口径歧义。
      pooled-OOF 与 macro mean 在本数据上给出实质不同的估计（SUDO 0.511 vs 0.756）。

**此前未被指出的问题**：这两个数字不只是在"聚合方式"上不同，它们用的**task 子集也不同**。
  - pooled-OOF  = 全部 202 个 task 的 OOF 分数放在一起算
  - macro mean  = 只有 25 个「可定义折」（n+>0 且 n->0）对应的 task 才进入计算
因此直接的 pooled-vs-macro 对比是 apples-to-oranges。

本脚本把差异分解为两个来源：
  A) 全部 202 task  的 pooled-OOF            （= 之前报告的 "pooledOOF"）
  B) 可定义折 task  的 pooled-OOF            （与 macro 同一 task 子集）
  C) 可定义折       的 macro mean            （= 之前报告的 "macro mean"）
  A vs B  → 「task 子集」效应
  B vs C  → 「聚合方式」效应

两者各配 fold-level cluster bootstrap（重抽 25 个可定义折）。
★ 2026-10-06（CS-11）：本节**不含** LOO-organism —— 原段误调 lodo_oof（把 LODO 标成 LOO），已删除；canonical Fig 3B 的 LOO 路径是 p14_rebuild_fig3b.py，不在本脚本内。

不改动任何既有结果，只做分解与 CI。
输出：report_p6f_lodo_settlement.txt / p6f_lodo_settlement.csv
"""
import os
import numpy as np, pandas as pd

BASE = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))); OUT = os.path.dirname(os.path.abspath(__file__))
SEED = 20261005; rng = np.random.default_rng(SEED); B = 4000

rep = open(os.path.join(OUT, "report_p6f_lodo_settlement.txt"), "w", encoding="utf-8")
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
    # 秩和法（精确处理 ties）
    r = np.empty(len(ss)); i = 0
    while i < len(ss):
        j = i
        while j + 1 < len(ss) and ss[j + 1] == ss[i]: j += 1
        r[i:j + 1] = 0.5 * (i + j) + 1
        i = j + 1
    n1 = yy.sum(); n0 = len(yy) - n1
    if n1 == 0 or n0 == 0: return np.nan
    return (r[yy == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)

# ---- fold 构成 ----
fold_info = {}
for gg in GROUPS:
    m = g == gg
    fold_info[gg] = dict(n=int(m.sum()), npos=int(y[m].sum()), nneg=int((1 - y[m]).sum()))
EVAL = [gg for gg in GROUPS if fold_info[gg]["npos"] > 0 and fold_info[gg]["nneg"] > 0]
NONEVAL = [gg for gg in GROUPS if gg not in EVAL]
tmask_eval = np.isin(g, EVAL)
p("=" * 112)
p("LODO aggregation ambiguity —— 分解为「task 子集」与「聚合方式」两个来源")
p("=" * 112)
p(f"  全部 fold = {len(GROUPS)}   可定义 AUROC 的 fold = {len(EVAL)}   不可定义 = {len(NONEVAL)}")
p(f"  全部 task = {len(E)}   可定义折覆盖的 task = {int(tmask_eval.sum())}   被排除 = {int((~tmask_eval).sum())}")
p(f"  被排除的 task 中 harmful = {int(y[~tmask_eval].sum())} / {int((~tmask_eval).sum())}")
p(f"\n  不可定义 fold（{len(NONEVAL)} 个）:")
for gg in NONEVAL:
    p(f"    {gg:<12} n={fold_info[gg]['n']}  n+={fold_info[gg]['npos']}  n-={fold_info[gg]['nneg']}")

def lodo_oof(r):
    oof = np.full(len(r), np.nan)
    from sklearn.linear_model import LogisticRegression
    for gg in GROUPS:
        te = g == gg; tr = ~te
        if len(np.unique(y[tr])) < 2 or te.sum() == 0: continue
        mu, sd = np.nanmean(r[tr]), np.nanstd(r[tr]) or 1.0
        try:
            m = LogisticRegression(max_iter=1000).fit(((r[tr] - mu) / sd).reshape(-1, 1), y[tr])
            oof[te] = m.predict_proba(((r[te] - mu) / sd).reshape(-1, 1))[:, 1]
        except Exception:
            continue
    return oof

def boot(sub_groups, fn):
    """fold-level cluster bootstrap：重抽 fold，对每次重抽算 fn(选中 task 的并集)"""
    vals = []
    for _ in range(B):
        pick = rng.choice(len(sub_groups), len(sub_groups), replace=True)
        idx = np.concatenate([np.where(g == sub_groups[i])[0] for i in pick])
        v = fn(idx, pick) if fn.__code__.co_argcount >= 2 else fn(idx)
        if np.isfinite(v): vals.append(v)
    return (np.percentile(vals, [2.5, 97.5]) if vals else (np.nan, np.nan))

rows = []
p("\n" + "=" * 112)
p("【分解】A = 全部 task 的 pooled-OOF ｜ B = 可定义折 task 的 pooled-OOF ｜ C = 可定义折 macro mean")
p("=" * 112)
p(f"{'method':<22}{'role':<11}{'A: pooled(all)':>15}{'B: pooled(eval)':>16}{'C: macro(eval)':>15}"
  f"{'A-B':>8}{'B-C':>8}{'C 95%CI':>18}")
for name, role, fn in METHODS:
    r = fn(E).values.astype(float)
    oof = lodo_oof(r)
    A = auc(y, oof)
    Bv = auc(y[tmask_eval], oof[tmask_eval])
    fa = {}
    for gg in EVAL:
        m = g == gg
        a = auc(y[m], oof[m])
        if np.isfinite(a): fa[gg] = a
    C = float(np.mean(list(fa.values())))
    lo, hi = boot(EVAL, lambda idx: auc(y[idx], oof[idx]))   # 与 B 同口径的 CI
    flo, fhi = boot(EVAL, lambda idx, pick: float(np.mean([fa[EVAL[i]] for i in pick if EVAL[i] in fa]))
                    if any(EVAL[i] in fa for i in pick) else np.nan)
    rows.append(dict(method=name, role=role, A_pooled_all202=A, B_pooled_eval=Bv, C_macro_eval=C,
                     AtB=A - Bv, BtC=Bv - C, C_lo=flo, C_hi=fhi, B_lo=lo, B_hi=hi,
                     n_eval_fold=len(fa)))
    p(f"{name:<22}{role:<11}{A:>15.3f}{Bv:>16.3f}{C:>15.3f}{A-Bv:>+8.3f}{Bv-C:>+8.3f}"
      f"{'[' + format(flo, '.3f') + ',' + format(fhi, '.3f') + ']':>18}")
pd.DataFrame(rows).to_csv(os.path.join(OUT, "p6f_lodo_settlement.csv"), index=False)

D = pd.DataFrame(rows)
p("\n  分解汇总（跨 15 个方法）:")
p(f"    |A-B| 中位（task 子集效应）= {D.AtB.abs().median():.3f}   max = {D.AtB.abs().max():.3f}")
p(f"    |B-C| 中位（聚合方式效应）= {D.BtC.abs().median():.3f}   max = {D.BtC.abs().max():.3f}")
p(f"    |A-C| 中位（原先报告的差距）= {(D.A_pooled_all202-D.C_macro_eval).abs().median():.3f}")

p("\n  口径 C 的 95%CI 是否含 0.5（即该指标在留出分组上是否与随机不可区分）:")
for _, r in D.iterrows():
    v = "CI 含 0.5（与随机不可区分）" if (r.C_lo <= 0.5 <= r.C_hi) else \
        ("CI 全 > 0.5（显著优于随机）" if r.C_lo > 0.5 else "CI 全 < 0.5（显著差于随机）")
    p(f"    {r.method:<22} [{r.C_lo:.3f},{r.C_hi:.3f}]  {v}")

# ★ 2026-10-06（CS-11）：原此处为 LOO-organism 段，已删除 —— 它误调 lodo_oof。
#   （原删除在此留下一条重复的 "=" 分隔线，是残迹，一并去掉。）
p("【结论 —— 供正文引用】")
p("=" * 112)
p("  1. 此前报告的 pooled-OOF 与 macro mean 之差，同时包含两种来源，二者必须分开陈述。")
p("  2. 若把 task 子集统一（都只用可定义折覆盖的 task），剩余差距即为纯粹的聚合方式差异。")
p(f"     实测：task 子集效应 |A-B| 中位 {D.AtB.abs().median():.3f}，聚合方式效应 |B-C| 中位 {D.BtC.abs().median():.3f}。")
p("  3. 本 protocol 冻结了 leave-one-organism-drug-out，但未冻结 fold 如何聚合；")
p("     两种口径均并列报告，不指定 retroactive primary。")
p("  4. 正文须写明：the aggregation strategy was not prespecified and produced materially different estimates。")
# ★ 2026-10-06：原写死 "25/38 折可定义；被排除的 13 折中 10 折为「全 non-harmful」"
#   与 "38 个分组"。eligible 宇宙改为 186 task / 37 组后全部失效（实测 26/37、11 折、
#   其中 10 折全 non-harmful、1 折全 harmful）。改为从 EVAL/NONEVAL/fold_info 现算，
#   这些量随宇宙变化，写死必错。
_p0 = sum(1 for gg in NONEVAL if fold_info[gg]["npos"] == 0)
_n0 = sum(1 for gg in NONEVAL if fold_info[gg]["nneg"] == 0)
p(f"  5. {len(EVAL)}/{len(GROUPS)} 折可定义；被排除的 {len(NONEVAL)} 折中 {_p0} 折为「全 non-harmful」"
  f"、{_n0} 折为「全 harmful」，故 macro mean 建立在非随机、更异质的子集上，")
p(f"     不得表述为 {len(GROUPS)} 个分组的代表值。")
rep.close()
print("done")
