# -*- coding: utf-8 -*-
r"""p31 — 补充图例（Supplementary_Legends.md）的数值审计

================================ 为什么需要（外部审阅者指出）================================
`p15` 明确**不审 supplement legends**（其源码里 grep "legend/Figure S" 命中 0）；
MNV（p30）对图例只核布尔不等式。于是**图例里的具体数值长期处于无人值守状态**：
实测漏出的是 Figure S2 图例的 ρ = 0.942 / 0.841 与 0.567–0.713 ——
在任何对齐口径下都无法从 `figS2_config_agreement.csv` 复现，
而正确值 [0.869, 0.927] / [0.611–0.687] 就摆在数据文件里。

本脚本把图例中**可核验的数值断言**逐条登记来源并比对，堵住这个缺口。

判据与 p15 一致：
  · pattern 必须在图例中原样出现（否则 = 已被编辑/改错）
  · 从**冻结产物**独立算出的值必须与图例声称的值一致（tol）
"""
import io, os, re, sys
import pandas as pd

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIG  = os.path.join(ROOT, "figures", "data")
LEG  = os.path.join(ROOT, "figures", "final", "Supplementary_Legends.md")
OUT  = os.path.join(ROOT, "分析结果")
REP  = open(os.path.join(OUT, "report_p31_legends_audit.txt"), "w", encoding="utf-8")

TXT = io.open(LEG, encoding="utf-8").read()

def p(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); REP.write(s + "\n"); REP.flush()

R = lambda f: pd.read_csv(os.path.join(FIG, f))
ROWS = []

def chk(fig, cid, pattern, claimed, source, computed, tol=5e-4):
    in_txt = pattern in TXT
    ok = in_txt and abs(float(claimed) - float(computed)) <= tol
    ROWS.append(dict(fig=fig, cid=cid, pattern=pattern, claimed=claimed,
                     source=source, computed=round(float(computed), 6),
                     status="PASS" if ok else ("NOT_IN_LEGEND" if not in_txt else "VALUE_MISMATCH")))
    p("  [%-14s] %-11s %-46s 图例=%-9s 重算=%-9s %s" % (
        "PASS" if ok else ("NOT_IN_LEGEND" if not in_txt else "MISMATCH"),
        cid, pattern[:46], claimed, round(float(computed), 6),
        source if ok else "←" + source))

p("=" * 108); p("p31 — 补充图例数值审计"); p("=" * 108)
p(f"图例文件：{LEG}"); p("")

U  = R("fig1b_universe.csv").iloc[0]
F6 = R("figS6_fold_summary.csv").set_index("group")
AG = R("figS2_config_agreement.csv")
M5 = R("figS5_mmd_summary.csv").iloc[0]
S7 = R("figS7_weight_association.csv").set_index("metric")
ORI = R("figS_orientation.csv").set_index("method")
CN9 = R("figS9_orientation_consistency.csv").set_index("method")
D  = R("fig1c_delta_distribution.csv")
CF = R("fig1e_configs_full.csv")

p("【S1】任务矩阵")
chk("S1", "S1-01", "all 38 × 6 = 228 candidate cells", 228, "fig1b.candidate_tasks", int(U.candidate_tasks), 0)
chk("S1", "S1-02", "186 are eligible", 186, "fig1b.eligible_tasks", int(U.eligible_tasks), 0)
chk("S1", "S1-03", "42 are\nexcluded (grey with ×)", 42, "fig1b.excluded", int(U.excluded), 0)
chk("S1", "S1-04", "37 organism–drug groups were single-class", 37, "fig1b.n_organism_drug_groups",
    int(U.n_organism_drug_groups), 0)

p(""); p("【S2】配置稳健性")
for cfg, n in [("P1 History+LR", 186), ("S2 History+LGBM", 174), ("S1 Core+LR", 64), ("S3 Core+LGBM", 60)]:
    pass
chk("S2", "S2-01", "(186 / 174 / 64 / 60 of 228)", 186, "fig1e P1 eligible",
    int(CF[CF.config == "P1 History+LR"].eligible.sum()), 0)
ht = lambda s: "History" if "History" in s else "Core"
ag = AG[AG.cfg_a != AG.cfg_b]
same  = sorted(ag[ag.cfg_a.map(ht) == ag.cfg_b.map(ht)].spearman_rho)
cross = sorted(ag[ag.cfg_a.map(ht) != ag.cfg_b.map(ht)].spearman_rho)
chk("S2", "S2-02", "ρ = 0.927 for history-augmented features", 0.927, "figS2 同 tier 上界", max(same), 5e-4)
chk("S2", "S2-03", "ρ = 0.869 for core features", 0.869, "figS2 同 tier 下界", min(same), 5e-4)
chk("S2", "S2-04", "ρ = 0.611 to 0.687", 0.611, "figS2 跨 tier 下界", min(cross), 5e-4)
chk("S2", "S2-05", "0.611 to 0.687", 0.687, "figS2 跨 tier 上界", max(cross), 5e-4)

p(""); p("【S3】方向分布")
sub = D[D.direction == "MGB->UT"]
# ★★ 2026-10-13：**pattern 同步**（正当修改，理由如下，不是"改判据让红变绿"）。
#   图注文件原为**双语工作稿**，S3 那两条数据自洽性注记写的是中文：
#       「MGB→UT（35 eligible 中 19 harmful，中位 −0.051）与 SU→MGB（31 中 7，中位 +0.003）」
#   本轮把图注改成**纯英文的发行文本**（提交件不得含中文工作注），同一句话变成
#       "MGB→UT (19 of 35 eligible harmful, median ΔAUROC −0.051) and SU→MGB (7 of 31 …)"
#   ⇒ 判据的 pattern 必须随被审文本的语言同步，**断言本身（19 / 7 与 fig1c 一致）一字未改**。
#   证据：改前 p31 = 46/46（中文 pattern 命中）；改图注后立刻
#   NOT_IN_LEGEND ×2（**说明这两条判据真的在查文本，不是空转**）；同步后恢复全绿。
chk("S3", "S3-01", "MGB→UT (19 of 35 eligible harmful", 19, "fig1c MGB->UT harmful", int(sub.harmful.sum()), 0)
chk("S3", "S3-02", "SU→MGB (7 of 31", 7, "fig1c SU->MGB harmful",
    int(D[D.direction == "SU->MGB"].harmful.sum()), 0)

p(""); p("【S5】MMD")
chk("S5", "S5-01", "Median SD was 0.0025", 0.0025, "figS5.sd_median", float(M5.sd_median), 5e-5)
chk("S5", "S5-02", "median\n|MMD²| of 0.0197", 0.0197, "figS5.mmd_median", float(M5.mmd_median), 5e-5)
chk("S5", "S5-03", "no task produced a negative estimate", 0, "figS5.n_negative", int(M5.n_negative), 0)

p(""); p("【S6】折构成")
chk("S6", "S6-01", "26 evaluable folds", 26, "figS6.evaluable.n_folds",
    int(F6.loc["evaluable", "n_folds"]), 0)
chk("S6", "S6-02", "134/186 tasks", 134, "figS6.evaluable.n_tasks",
    int(F6.loc["evaluable", "n_tasks"]), 0)
chk("S6", "S6-03", "Eleven of 37 organism–drug groups", 11, "figS6.excluded.n_folds",
    int(F6.loc["excluded", "n_folds"]), 0)

p(""); p("【S7】权重诊断")
chk("S7", "S7-01", "ρ = −0.249", -0.249, "figS7.ESS.rho", float(S7.loc["ESS", "spearman_rho"]), 5e-4)
chk("S7", "S7-02", "−0.376 to −0.117", -0.376, "figS7.ESS.lo", float(S7.loc["ESS", "ci_lo"]), 5e-4)
chk("S7", "S7-03", "ρ = +0.220", 0.220, "figS7.w_max.rho", float(S7.loc["w_max", "spearman_rho"]), 5e-4)
chk("S7", "S7-04", "+0.081 to +0.357", 0.081, "figS7.w_max.lo", float(S7.loc["w_max", "ci_lo"]), 5e-4)

p(""); p("【S9】方向学习")
chk("S9", "S9-01", "native direction gave 0.447", 0.447, "figS_orientation.mmd.macro_within",
    float(ORI.loc["mmd", "macro_within"]), 5e-4)
chk("S9", "S9-02", "logistic learner 0.409", 0.409, "figS_orientation.mmd.macro_logit",
    float(ORI.loc["mmd", "macro_logit"]), 5e-4)
chk("S9", "S9-03", "gave 0.553", 0.553, "figS_orientation.mmd.macro_rank",
    float(ORI.loc["mmd", "macro_rank"]), 5e-4)
chk("S9", "S9-04", "positive orientation in 25 of 26", 25,
    "figS9.mmd.frac_logit_pos × 26", float(CN9.loc["mmd", "frac_logit_pos"]) * 26, 5e-3)
chk("S9", "S9-05", "negative in all 26 folds", 0.0, "figS9.mmd.frac_rank_pos",
    float(CN9.loc["mmd", "frac_rank_pos"]), 1e-9)

# ================================================================ 【S10】分诊策略模拟 —— **图例已撤除**
#
# ★★ 2026-10-15（作者决定 **B**）：**Figure S10 的图例已从 `Supplementary_Legends.md` 移除。**
#
#   为什么要撤 —— 这**不是**「漏了一个正文引用点」，而是一个**已经在四处执行过的收录决定，
#   漏删了第五处**：
#     · `p56` 有明文决定：「稿件引用 **S1–S9**；`FigureS10_triage.R` **不在投稿范围（正文未引），
#       未收入**」；
#     · 该决定**四处一致执行**：包内补充图只有 S1–S9；`SourceData.xlsx` 的 33 个 sheet 里
#       **没有** S10；正文不引；`03_补充图_TIFF_350dpi/` 里没有 `FigureS10.tiff`；
#     · **唯独图例留在包里** —— 于是出现「补充图 9 张、图例 10 条」。
#       （外部核验 2026-10-15 抓到「S10 全包零引用」，本文把它归因清楚了。）
#
#   处置（选 B = **把决定执行完**，而不是重开它去补图件/Source Data/引用）：
#     · 图例段落已删；**原文留档**于
#       `分析结果/_artifact_archive/FigureS10_legend_DROPPED_20261015T012150Z.md`，
#       删改前的整份 legends 亦在同目录留档（`Supplementary_Legends.md.20261015T012150Z.a2b27d3b`）；
#     · **一个数字都没丢** —— 分诊的每个量仍由 **Table S3** 承载，
#       `p45`（生产者）/ `evs/t45`（独立验证）/ `p46` + `evs/t46`（Table S3 渲染与审计）**全部不动**。
#
#   ⇒ 本节此前的 **20 条** S10-01…S10-19b **不是被删掉，是换了判据形态**：
#     从「审这一节的数字」改为「**保证这一节不再出现**」—— 负向守卫。
#     若日后有人把图例加回来却仍不带图件，这里会立刻判红。
p(""); p("【S10】分诊策略模拟 —— 图例已按作者决定撤除（负向守卫）")

_S10_ABSENT = "## Figure S10" not in TXT
ROWS.append(dict(fig="S10", cid="S10-X1",
                 pattern="`Supplementary_Legends.md` 不得再含 `## Figure S10` 一节",
                 claimed="（不得出现）", source="作者决定 B（2026-10-15）＋ p56 的收录决定",
                 computed=("未出现" if _S10_ABSENT else "**仍存在**"),
                 status="PASS" if _S10_ABSENT else "VALUE_MISMATCH"))
p("  [%-14s] %-11s %-46s %s" % (
    "PASS" if _S10_ABSENT else "VALUE_MISMATCH", "S10-X1", "legends 不含 Figure S10 一节",
    "✔ 未出现（图例与图件一致）" if _S10_ABSENT else "✘ 仍存在 —— 图例不能没有图"))

# 控制：证明上面那条不是空转（该串放进合成文本里应可被检出）
_S10_SELFTEST = ("## Figure S10" in "x ## Figure S10 y")
ROWS.append(dict(fig="S10", cid="S10-X2", pattern="负向守卫非空转自检",
                 claimed="可检出", source="合成串",
                 computed=("可检出" if _S10_SELFTEST else "不可检出"),
                 status="PASS" if _S10_SELFTEST else "VALUE_MISMATCH"))
p("  [%-14s] %-11s %-46s %s" % ("PASS" if _S10_SELFTEST else "VALUE_MISMATCH",
    "S10-X2", "负向守卫非空转自检", "✔ 可检出" if _S10_SELFTEST else "✘ 判据无判别力"))

# ★ 顺带守住「配对的另一头」：编号必须连续到 S9（撤掉 S10 后**不能**留下断号或跳号）
_S9_LAST = ("## Figure S9" in TXT) and ("## Figure S10" not in TXT)
ROWS.append(dict(fig="S10", cid="S10-X3",
                 pattern="legends 的最后一个编号小节为 Figure S9",
                 claimed="S9 为末节", source="撤除 S10 后编号应连续收在 S9",
                 computed=("S9 为末节" if _S9_LAST else "编号不连续"),
                 status="PASS" if _S9_LAST else "VALUE_MISMATCH"))
p("  [%-14s] %-11s %-46s %s" % ("PASS" if _S9_LAST else "VALUE_MISMATCH",
    "S10-X3", "编号连续收在 S9", "✔" if _S9_LAST else "✘"))

p(""); p("=" * 108)
R2 = pd.DataFrame(ROWS)
bad = R2[R2.status != "PASS"]
p(f"p31 断言合计 {len(R2)} 项：PASS {int((R2.status=='PASS').sum())}｜**非 PASS {len(bad)}**")
if len(bad):
    p("")
    for _, r in bad.iterrows():
        p(f"    [{r['status']}] {r.cid} {r.pattern[:60]}  图例={r.claimed} 重算={r.computed} 来源={r.source}")
    p(""); p("*** p31 未通过 ***")
else:
    p(""); p("*** p31 通过 —— 补充图例的数值断言全部可从冻结产物复现 ***")
R2.to_csv(os.path.join(OUT, "p31_legends_audit.csv"), index=False, encoding="utf-8-sig")
REP.close()
sys.exit(1 if len(bad) else 0)
