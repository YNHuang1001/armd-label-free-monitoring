# -*- coding: utf-8 -*-
r"""t20_identity_gate.py — ★ BENCHMARK IDENTITY GATE（CORRECTION_SPEC §2.1）

在 CS-01/02/03 已进生产之后，回答 spec §2.1 的三个问题，并（按 GPT 2026-10-05 的要求）
输出四层对照：

  Layer 1  organism-drug 组数        38 -> ?
  Layer 2  candidate transfer / eligible task   228 -> ? / 202 -> ?
  Layer 3  **掉了哪些、为什么**（逐组逐 task，含原因归因）
  Layer 4  性能变化预估（样本量、train/test 构成、source AUROC 筛选的变化）

数据来源（**均为实测产物**，不是推演）：
  修正前 = FREEZE_pre_correction_20261005/PRODUCTION_outputs_backup/分析结果/
  修正后 = 分析结果/（2026-10-06 用 CS-01/02/03 修后的生产代码重跑 p5 得到）

★ spec §2.1 的硬约束：若三个数字有变，后面任何"旧结果 vs 新结果"的比较
  **只能作为 correction impact analysis**，**不得强行维持原来的 202 universe**。
"""
import os, io, sys, json
import numpy as np
import pandas as pd

import evs_load as L

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

PROJ = L.PROJ
OUT = os.path.join(PROJ, "分析结果")
BAK = os.path.join(PROJ, "FREEZE_pre_correction_20261005",
                   "PRODUCTION_outputs_backup", "分析结果")
REPORT = os.path.join(PROJ, "BENCHMARK_IDENTITY_GATE.md")
AUC_THR = 0.65                    # 预先写死的 source 合格线（p5 的 availability screen）

rep = io.open(REPORT, "w", encoding="utf-8")


def p(*a):
    line = " ".join(str(x) for x in a)
    print(line, flush=True)
    rep.write(line + "\n")
    rep.flush()


def rd(base, f):
    q = os.path.join(base, f)
    return pd.read_csv(q) if os.path.exists(q) else None


p("# BENCHMARK IDENTITY GATE")
p("")
p(f"**执行时点**：2026-10-06 ｜ **触发条件**：CS-01（药物映射，含 **v2 别名归属修复**）、"
  f"CS-02（episode 单元）、CS-03（S/R 冲突裁决）已进生产并各自验证通过；"
  f"p4/p5 已在修正后重跑")
p(f"**依据**：CORRECTION_SPEC §2.1 ｜ 对照基准 = 冻结快照 "
  f"`PRODUCTION_outputs_backup/`（修正前产物）")
p("")
p("---")
p("")

old_t = rd(BAK, "p5_tasks_primary.csv")
new_t = rd(OUT, "p5_tasks_primary.csv")
old_g = rd(BAK, "p5_gate.csv")
new_g = rd(OUT, "p5_gate.csv")
# ★ 2026-10-06：原在这里读入 old_u / new_u（p6b_corrected_universe.csv），
#   但**全文从未使用**（死变量）。已删除，理由是它掩盖了一件事实：
#   本 gate **不**检查 p6b 的产物，因此它的判定不依赖下游是否已重跑。
#   保留死变量会让读者以为"这份 gate 验过 p6b"，而实际上那份 CSV 至今
#   仍是 2026-09-29 的修正前副本。

assert old_t is not None and new_t is not None, "缺 p5_tasks_primary.csv（前/后）"


def key(df):
    """★ 合并键必须是**完整 task 身份**：organism|drug|source|target。
    第一版只用了 organism|drug —— 而一个组有 6 个有向迁移，于是 merge 变成组内笛卡尔积，
    Layer 3/4 的 n_target「前→后」是**错配对**的数字（同一行会重复出现）。已修正。"""
    return (df.organism.astype(str) + "|" + df.drug.astype(str) + "|"
            + df.source.astype(str) + "|" + df.target.astype(str))


# ---------------------------------------------------------------- Layer 1
grp = lambda d: sorted(set(d.organism.astype(str) + "|" + d.drug.astype(str)))
old_groups = grp(old_t)
new_groups = grp(new_t)
lost_groups = sorted(set(old_groups) - set(new_groups))
new_groups_only = sorted(set(new_groups) - set(old_groups))

p("## Layer 1 — organism-drug 组数")
p("")
p("| | 修正前 | 修正后 |")
p("|---|---|---|")
p(f"| organism-drug 组数（**candidate**，即 228 迁移覆盖的组） | **{len(old_groups)}** | **{len(new_groups)}** |")
p(f"| 菌数 | {old_t.organism.nunique()} | {new_t.organism.nunique()} |")
p(f"| 药数 | {old_t.drug.nunique()} | {new_t.drug.nunique()} |")

# ★★ 2026-10-06 补：**eligible 覆盖的组数**。
#   为什么必须有这一行：Layer 1 上行报的是 **candidate** 集合的组数，它对 CS-01/02/03
#   **不敏感**（228 个迁移一个没少）。若只报那一行，读者会把"38 组不变"读成
#   **"被分析的就是 38 组"** —— 而实测：eligible 只覆盖 **37** 组，
#   `KP|NIT` 的 6 个迁移在修正后 src_auroc 全部落到 0.65 以下（0.6217/0.6458/0.6384），
#   整组出局。这是**仅看 candidate 计数就会漏掉**的实质变化。
def egrp(d):
    return sorted(set(d.loc[d.eligible.astype(bool), "organism"].astype(str) + "|"
                      + d.loc[d.eligible.astype(bool), "drug"].astype(str)))
_eo, _en = egrp(old_t), egrp(new_t)
_lost_g = sorted(set(_eo) - set(_en))
p(f"| organism-drug 组数（**eligible 覆盖**，即真正进入下游分析的组） | **{len(_eo)}** | **{len(_en)}** |")
if _lost_g:
    p("")
    p(f"★ **eligible 覆盖的组少了 {len(_lost_g)} 个**：{', '.join(_lost_g)} —— "
      "该组全部 6 个迁移都掉出合格集。")
    p("★ 表述要求：稿件若写「38 个 organism-drug 组」，**必须限定为 candidate 集合**；"
      "凡描述**被分析**的宇宙，须用 eligible 覆盖的组数（本期为 "
      f"{len(_en)}）。两者不是同一个量。")
else:
    p("")
    p("（eligible 覆盖的组数与 candidate 相同。）")
p("")
if lost_groups:
    p(f"**candidate 里消失的组（{len(lost_groups)}）**：" + ", ".join(lost_groups))
else:
    p(f"**candidate 集合没有组消失**：{len(old_groups)} 个 organism-drug 组在 228 迁移层面全部保留。"
      f"★ 注意与上表区的区别 —— eligible 覆盖的组少 1 个**不是**「组消失」，"
      f"而是「该组的迁移**全部不合格**」：两者含义不同，不得混用措辞。")
if new_groups_only:
    p(f"**新增的组（{len(new_groups_only)}）**：" + ", ".join(new_groups_only))
p("")

# ---------------------------------------------------------------- Layer 2
n_cand_o, n_cand_n = len(old_t), len(new_t)
n_el_o = int(old_g.n_eligible.iat[0]) if old_g is not None else int(old_t.eligible.sum())
n_el_n = int(new_g.n_eligible.iat[0]) if new_g is not None else int(new_t.eligible.sum())
p("## Layer 2 — candidate transfer / eligible task")
p("")
p("| | 修正前 | 修正后 | 变化 |")
p("|---|---|---|---|")
p(f"| candidate transfer | {n_cand_o} | {n_cand_n} | {n_cand_n - n_cand_o:+d} |")
p(f"| eligible task（source AUROC ≥ {AUC_THR}） | **{n_el_o}** | **{n_el_n}** | "
  f"**{n_el_n - n_el_o:+d}** |")
p(f"| excluded | {n_cand_o - n_el_o} | {n_cand_n - n_el_n} | {(n_cand_n - n_el_n) - (n_cand_o - n_el_o):+d} |")
if old_g is not None and new_g is not None:
    p(f"| N_harmful (Δ ≤ −0.05) | {int(old_g.n_harmful.iat[0])} | **{int(new_g.n_harmful.iat[0])}** | "
      f"{int(new_g.n_harmful.iat[0]) - int(old_g.n_harmful.iat[0]):+d} |")
    p(f"| N_stable (\\|Δ\\| ≤ 0.02) | {int(old_g.n_stable.iat[0])} | **{int(new_g.n_stable.iat[0])}** | "
      f"{int(new_g.n_stable.iat[0]) - int(old_g.n_stable.iat[0]):+d} |")
    p("")
    p("★ **这两个绝对计数的变化不得表述为「模型性能下降」**——它们改变是因为**被评价的 "
      "task universe 被重建**（分母变了），不是同一个 universe 上的性能变差。正确表述：")
    p("")
    p("> The absolute counts changed because the evaluated task universe was reconstructed.")
    p(f"| P90−P10(ΔAUROC) | {old_g.spread_p90_p10.iat[0]:.4f} | {new_g.spread_p90_p10.iat[0]:.4f} | "
      f"{new_g.spread_p90_p10.iat[0] - old_g.spread_p90_p10.iat[0]:+.4f} |")
p("")

# ---------------------------------------------------------------- Layer 3
old_t2 = old_t.copy(); old_t2["grp"] = key(old_t)
new_t2 = new_t.copy(); new_t2["grp"] = key(new_t)
m = old_t2.merge(new_t2, on="grp", suffixes=("_o", "_n"))
assert len(m) == len(old_t2) == len(new_t2), f"合并后行数 {len(m)} != {len(old_t2)}（键不唯一）"
m["flip"] = m.eligible_o.astype(bool) & (~m.eligible_n.astype(bool))
flip = m[m.flip].copy()
p("## Layer 3 — 掉了哪些 task、为什么")
p("")
# ★ 必须**同时**给出 removed / newly / net —— 只列 removed 会让读者问
#   「为什么列 24 却只掉 16」（净变化 = removed − newly）。
_newly = n_el_n - (n_el_o - len(flip))
p("| | 数量 |")
p("|---|---|")
p(f"| **Removed eligible tasks**（eligible → not eligible） | **{len(flip)}** |")
p(f"| **Newly qualified tasks**（not eligible → eligible） | **{_newly}** |")
p(f"| **Net eligible change** | **{n_el_n - n_el_o:+d}**（{n_el_o} → {n_el_n}） |")
p("")
p(f"（{len(flip)} − {_newly} = {len(flip) - _newly}，与 Layer 2 的净变化一致）")
p("")
if len(flip):
    cols = ["grp", "organism_o", "drug_o", "source_o", "target_o", "src_auroc_o", "src_auroc_n",
            "n_src_test_o", "n_src_test_n", "n_target_o", "n_target_n", "d_auroc_o", "d_auroc_n"]
    cols = [c for c in cols if c in flip.columns]
    f2 = flip[cols].sort_values("src_auroc_n")
    # ★ 拆成独立列（外部审计意见）：原来把完整 task 身份放在一格、又单列一列 source→target，
    #   两处重复，读者会误以为 `UT→SU` 就是 task 身份本身。现按 organism / drug / source / target 分列。
    p("| organism | drug | source site | target site | AUROC 前→后 | n_src_test 前→后 | n_target 前→后 | 归因 |")
    p("|---|---|---|---|---|---|---|---|")
    for _, r in f2.iterrows():
        why = []
        if r.src_auroc_n < AUC_THR <= r.src_auroc_o:
            why.append(f"source AUROC 掉到 {r.src_auroc_n:.3f} < {AUC_THR}")
        if r.n_target_n < r.n_target_o:
            why.append(f"n_target −{int(r.n_target_o - r.n_target_n):,}")
        if r.n_src_test_n < r.n_src_test_o:
            why.append(f"n_src_test −{int(r.n_src_test_o - r.n_src_test_n):,}")
        p(f"| {r.organism_o} | {r.drug_o} | {r.source_o} | {r.target_o} | "
          f"{r.src_auroc_o:.4f} → **{r.src_auroc_n:.4f}** | "
          f"{int(r.n_src_test_o):,} → {int(r.n_src_test_n):,} | "
          f"{int(r.n_target_o):,} → {int(r.n_target_n):,} | {'；'.join(why) or '—'} |")
    p("")
    p("### 归因（区分「原因」与「结果」）")
    p("")
    p("**因果链**：")
    p("")
    p("```")
    p("CS-01/02/03 修正  →  训练/评价样本被重建  →  source 模型性能被**重新计算**")
    p("                  →  AUROC 落在预设合格线 0.65 之下  →  task 出局")
    p("```")
    p("")
    p("★ **门槛规则本身没有改变**（仍是预先写死的 source 内部 AUROC ≥ 0.65）；"
      "**被纠正的是评价对象**——episode 单元与 S/R 结局定义。")
    p("")
    p("因此正确的表述是：")
    p("")
    p("> All removed tasks failed the **predefined** source-AUROC eligibility criterion after")
    p("> CS-01/02/03-induced reconstruction of the benchmark unit and outcome labels.")
    p("")
    p("**不得**写成「因为 AUROC 降了所以掉」——AUROC 下降是**结果**，不是原因。")
    p("")
    p("| 结果层的数字 | 计数 |")
    p("|---|---|")
    p(f"| removed task 的 AUROC 现在低于门槛 | {int((f2.src_auroc_n < AUC_THR).sum())} / {len(f2)} |")
    p(f"| 其中 n_target 同时减少 | {int((f2.n_target_n < f2.n_target_o).sum())} |")
    p(f"| 其中 n_src_test 同时减少 | {int((f2.n_src_test_n < f2.n_src_test_o).sum())} |")
else:
    p("**没有 task 掉出合格集。**")
p("")

# ---------------------------------------------------------------- Layer 4
p("## Layer 4 — 性能变化预估（不重训下游模型即可给出）")
p("")
d_src = m.n_src_test_n - m.n_src_test_o
d_tgt = m.n_target_n - m.n_target_o
d_auc = m.src_auroc_n - m.src_auroc_o
p("| 量 | 修正前（中位/合计） | 修正后 | 变化 |")
p("|---|---|---|---|")
p(f"| n_src_test | {m.n_src_test_o.median():,.0f} / {int(m.n_src_test_o.sum()):,} | "
  f"{m.n_src_test_n.median():,.0f} / {int(m.n_src_test_n.sum()):,} | "
  f"中位 {d_src.median():+,.0f}；合计 {int(d_src.sum()):+,} |")
p(f"| n_target | {m.n_target_o.median():,.0f} / {int(m.n_target_o.sum()):,} | "
  f"{m.n_target_n.median():,.0f} / {int(m.n_target_n.sum()):,} | "
  f"中位 {d_tgt.median():+,.0f}；合计 {int(d_tgt.sum()):+,} |")
_eo = m[m.eligible_o.astype(bool)]
_en = m[m.eligible_n.astype(bool)]
_keep = m[m.eligible_o.astype(bool) & m.eligible_n.astype(bool)]
p(f"| src_auroc 中位（**全部 228 candidate**） | {m.src_auroc_o.median():.4f} | "
  f"{m.src_auroc_n.median():.4f} | {d_auc.median():+.4f} |")
p(f"| src_auroc 中位（**eligible 子集**） | {_eo.src_auroc_o.median():.4f}（n={len(_eo)}） | "
  f"{_en.src_auroc_n.median():.4f}（n={len(_en)}） | "
  f"{_en.src_auroc_n.median() - _eo.src_auroc_o.median():+.4f} |")
p(f"| src_auroc 中位（**前后都合格**的 {len(_keep)} 个） | {_keep.src_auroc_o.median():.4f} | "
  f"{_keep.src_auroc_n.median():.4f} | "
  f"{_keep.src_auroc_n.median() - _keep.src_auroc_o.median():+.4f} |")
p("")
p(f"- ★ **存活任务的质量**：eligible 子集的 src_auroc 中位 "
  f"{_eo.src_auroc_o.median():.4f} → {_en.src_auroc_n.median():.4f}"
  f"（{_en.src_auroc_n.median() - _eo.src_auroc_o.median():+.4f}）——**掉出去的是边缘 task，不是随机掉**。")
p("")
p("★ ★ **不得过度包装**（外部审计意见）：**模型没有重新训练**，所以"
  "「存活任务的中位 AUROC 更高」**不等于**「修正提升了模型性能」。禁止写 "
  "*correction improved model performance*。正确表述二选一：")
p("")
p("> The correction reduced the evaluation universe by excluding transfers that no longer "
  "satisfied predefined source-quality criteria.")
p("")
p("> The retained transfers showed higher median source discrimination, reflecting removal "
  "of marginal-quality tasks.")
p("")
p(f"- **样本量整体变化**：source 侧合计 {int(d_src.sum()):+,}（{100.0*d_src.sum()/m.n_src_test_o.sum():+.1f}%），"
  f"target 侧合计 {int(d_tgt.sum()):+,}（{100.0*d_tgt.sum()/m.n_target_o.sum():+.1f}%）")
p(f"- **受影响最重的 5 个 task**（按 target 侧缩减比例）：")
worst = (m.assign(rel=(m.n_target_o - m.n_target_n) / m.n_target_o)
          .sort_values("rel", ascending=False).head(5))
for _, r in worst.iterrows():
    p(f"  - {r.grp} {r.source_o}→{r.target_o}：n_target {int(r.n_target_o):,} → {int(r.n_target_n):,}"
      f"（−{100*r.rel:.1f}%）")
p("")

# ---------------------------------------------------------------- 判读
p("---")
p("")
p("## 判读（CORRECTION_SPEC §2.1）")
p("")
q1 = "**不变**" if len(lost_groups) == 0 else f"**减少 {len(lost_groups)}**"
q2 = "**仍是 228**" if n_cand_n == n_cand_o else f"**{n_cand_o} → {n_cand_n}**"
q3 = "**仍是 202**" if n_el_n == n_el_o else f"**{n_el_o} → {n_el_n}（{n_el_n - n_el_o:+d}）**"
p(f"| 问题 | 结果 |")
p("|---|---|")
p(f"| **Q1** 38 个 organism–drug 组还剩多少？ | {len(new_groups)}（{q1}） |")
p(f"| **Q2** 228 个 candidate transfer 是否还是 228？ | {q2} |")
p(f"| **Q3** source AUROC ≥ {AUC_THR} 之后是否还是 202？ | {q3} |")
p("")
if n_el_n != n_el_o or len(lost_groups):
    p("★ **三个数字有变** ⇒ 按 §2.1 的硬约束：**后面任何「旧结果 vs 新结果」的比较只能作为 "
      "correction impact analysis，不得强行维持原来的 202 universe。**")
else:
    p("★ 三个数字均未变 ⇒ 旧结果与新结果可直接比较。")
p("")
p("**边界（2026-10-06 更正）**：")
p("")
p("· 本 gate 判定的是 **CS-01 / CS-01 v2 / CS-02 / CS-03** 的影响 —— 即改变"
  "**任务身份与结局定义**的那几项。")
p("· 原文写「CS-04…CS-12 尚未修，最终数字仍会再动一次」，**该句已作废**："
  "CS-05…CS-12 已于 2026-10-06 批量进生产（存档 "
  "`FREEZE_pre_correction_20261005/APPLIED_CORRECTIONS_20261006.txt`）。")
p("· **它们不会改动本节任何数字**，依据是实测而非推断：该批 37 处补丁落在 9 个文件上 ——")
p("  `p10a / p10b / p10c / p21a / p6_02 / p6b / p6c / p6d / p6f`，**无一触及 p4 或 p5**。")
p("  故本 gate 给出的 38 / 228 / 186 对本项目全部修正项而言**已是终值**。")
p("· 唯一提到 `p5` 的修正是 CS-10（Platt），但它是**文本条款**，且 Platt 是单调映射 ⇒ "
  "不改变 AUROC 类量（`src_auroc` / `d_auroc` / eligibility）。")
p("· 仍未发生的是**下游重跑**（p6b 起）：那是 `post_corrected_rerun` 阶段的进度，"
  "不是本 gate 的未决项。")
rep.close()
print("\ndone")
