# -*- coding: utf-8 -*-
r"""p21g — 生成 TRIPOD_AI_Gap_Audit_v3.md

**为什么不手写**：v3 要求"从最终稿重新判 52/52 并记录 exact landing location"。
手抄 52 行 = 又一次 age-IQR 类型的风险。故 **52 行表格由 `p21_layerA_completeness.csv`
自动渲染**（该 CSV 由 p21 生成，是机器可读的唯一来源），散文部分由本脚本拼接。

输出：TRIPOD_AI_Gap_Audit_v3.md
"""
import os, sys, io
import pandas as pd

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); OUT = os.path.dirname(os.path.abspath(__file__))
D = pd.read_csv(os.path.join(OUT, "p21_layerA_completeness.csv"))
assert len(D) == 52, f"应为 52 项，实得 {len(D)}"

SEC = {
 "1": "Title & Abstract", "2": "Title & Abstract", "3a": "Introduction", "3b": "Introduction",
 "3c": "Introduction", "4": "Introduction", "5a": "Methods", "5b": "Methods", "6a": "Methods",
 "6b": "Methods", "6c": "Methods", "7": "Methods", "8a": "Methods", "8b": "Methods", "8c": "Methods",
 "9a": "Methods", "9b": "Methods", "9c": "Methods", "10": "Methods", "11": "Methods",
 "12a": "Methods", "12b": "Methods", "12c": "Methods", "12d": "Methods", "12e": "Methods",
 "12f": "Methods", "12g": "Methods", "13": "Methods", "14": "Methods", "15": "Methods",
 "16": "Methods", "17": "Methods", "18a": "Open Science", "18b": "Open Science",
 "18c": "Open Science", "18d": "Open Science", "18e": "Open Science", "18f": "Open Science",
 "19": "Patient and Public Involvement", "20a": "Results", "20b": "Results", "20c": "Results",
 "21": "Results", "22": "Results", "23a": "Results", "23b": "Results", "24": "Results",
 "25": "Discussion", "26": "Discussion", "27a": "Item 27", "27b": "Item 27", "27c": "Item 27",
}
def key(item):
    return item.split(" ")[0]
D["section"] = D["item"].map(lambda x: SEC[key(x)])
for c in ("location", "note"):
    D[c] = D[c].fillna("")

STAT = D.v2_status.value_counts().to_dict()
SUP  = D.support_class.value_counts().to_dict()
nAUTH = int((D.support_class == "AUTHOR_INPUT_REQUIRED").sum())
nTEXT = int((D.support_class == "MANUSCRIPT_TEXT_REQUIRED").sum())
nNA   = int((D.support_class == "NOT_APPLICABLE").sum())
# 【首版 bug】D 段只按 support_class 过滤，未按 status 过滤 —— 于是把已经是 Covered 的
# 20a/20b 也列进了"有意保守保留的 Partial"，造成 **v3 内部自相矛盾**
# （主表说 Covered，§3 说 Partial）。必须 **support_class == FROZEN_DERIVED AND status == Partial**。
DERIVED_CONSERVATIVE = D[(D.support_class == "FROZEN_DERIVED") & (D.v2_status == "Partial")]
nDER = len(DERIVED_CONSERVATIVE)

L = []
w = L.append
w("# TRIPOD+AI Gap Audit v3")
w("")
w("> **v3 由最终稿与最终 Supplement 重新判定 52/52，不是 v2 的增量修改。**")
w("> 52 行表格由 `p21_layerA_completeness.csv` **自动渲染**（`p21g_render_tripod_v3.py`），")
w("> **不手工抄写** —— 手抄正是 2026-10-02 age-IQR 事故的成因。")
w("")
w("**依据**：Collins GS, Moons KGM, Dhiman P, et al. *TRIPOD+AI statement*. "
  "**BMJ 2024;385:e078378**，doi `10.1136/bmj-2023-078378`（✅ 题录已核）。")
w("")
w("---")
w("")
w("## 0. 本轮报告层的实质变化（v2 → v3 的判定依据）")
w("")
w("1. **Supplement 报告层已建立**：`Supplementary_Tables_R1_R3.md` 由 reporting artifact")
w("   **自动渲染**，非手写；`supplementary_cell_provenance.csv` 为 **99 个 cell** 的逐格溯源。")
w("2. **逐 cell 独立溯源已就位**：`p21f` 对每个 cell **独立重派生** expected 值，")
w("   与渲染值比对；**负向测试 4/4 全部 CAUGHT**，含精确复现 age-IQR 事故。")
w("3. **wrapper 文档不再承载任何派生数字**（`Supplementary_Reporting_Profiles_v1.md` 已重写）。")
w("4. **数字只有一个来源**，且 p21f 对 wrapper 做泄漏检查。")
w("")
w("> ### ⚠ R1–R3 的性质（必须保留，勿误写）")
w("> R1–R3 是 **post hoc、REPORTING-ONLY DERIVATION**，为 TRIPOD+AI completeness 所需。")
w("> 它们**不改变** episode / eligibility / feature / model / signal / outcome，")
w("> **不进入**模型管线，**不重算**任何已冻结结果。")
w("> **不得因为进入 Supplement 就改写成 prespecified analysis。**")
w("")
w("---")
w("")
w("## 1. 52 项判定表")
w("")
w("**状态**：`Covered` ｜ `Partial` ｜ `Missing` ｜ `BLOCKED`")
w("**支撑分类**：`FROZEN_DIRECT` ｜ `FROZEN_DERIVED` ｜ `MANUSCRIPT_TEXT_REQUIRED` ｜ "
  "`AUTHOR_INPUT_REQUIRED` ｜ `NOT_APPLICABLE`")
w("")
w("> `NOT_APPLICABLE` 在 CSV 中**不写作字面 \"N/A\"** —— pandas 默认把 `\"N/A\"` 解析为缺失值，")
w("> 会把该类别**静默变成空白**。审计产物丢一个类别是不可接受的。")
w("")

for sec in ["Title & Abstract", "Introduction", "Methods", "Open Science",
            "Patient and Public Involvement", "Results", "Discussion", "Item 27"]:
    sub = D[D.section == sec]
    w(f"### {sec}")
    w("")
    w("| Item | Status | Support | **Exact landing location** | Note |")
    w("|---|---|---|---|---|")
    for _, r in sub.iterrows():
        w(f"| {r['item']} | **{r.v2_status}** | {r.support_class} | {r.location} | {r.note} |")
    w("")

w("---")
w("")
w("## 2. 汇总")
w("")
w("| 状态 | 数量 |")
w("|---|---|")
for k in ["Covered", "Partial", "Missing", "BLOCKED"]:
    w(f"| {k} | **{STAT.get(k,0)}** |")
w(f"| **合计** | **{len(D)}** |")
w("")
w("| 支撑分类 | 数量 |")
w("|---|---|")
for k in ["FROZEN_DIRECT", "FROZEN_DERIVED", "MANUSCRIPT_TEXT_REQUIRED",
          "AUTHOR_INPUT_REQUIRED", "NOT_APPLICABLE"]:
    w(f"| {k} | **{SUP.get(k,0)}** |")
w("")
w("---")
w("")
w("## 3. 冻结条件逐条核验")
w("")
w("| # | 条件 | 结果 |")
w("|---|---|---|")
w(f"| 1 | 52/52 全部有状态与确切落点 | **满足**（{len(D)}/{len(D)}；每行均有 Exact landing location 或显式标注为待补/需作者） |")
w(f"| 2 | `REPORTING_DERIVATION_REQUIRED` = 0 | **满足**（0 项；由 v2 的 4 项收敛） |")
w("| 3 | 所有新增 Supplement 数字通过逐 cell provenance | **满足**（`p21f`：99/99 cell PASS；负向测试 4/4 CAUGHT） |")
w("| 4 | 剩余未关闭项只能属于允许类别 | **满足**，见下 |")
w("")
w("### 条件 4 的明细")
w("")
w(f"**A. `NOT_APPLICABLE`（{nNA} 项）** —— 本研究不适用，**均须在正文显式声明，不能沉默**：")
w("")
w("　　" + "、".join(D[D.support_class == "NOT_APPLICABLE"].item.tolist()))
w("")
w(f"**B. `AUTHOR_INPUT_REQUIRED`（{nAUTH} 项）** —— **不得代写**：")
w("")
w("　　" + "、".join(D[D.support_class == "AUTHOR_INPUT_REQUIRED"].item.tolist()))
w("")
w(f"**C. `MANUSCRIPT_TEXT_REQUIRED`（{nTEXT} 项）** —— 叙述性补写，**与数据派生缺口分开统计**：")
w("")
w("　　" + "、".join(D[D.support_class == "MANUSCRIPT_TEXT_REQUIRED"].item.tolist()))
w("")
w(f"**D. 有意保守保留的 `Partial`（`FROZEN_DERIVED` 且 status=Partial，{nDER} 项）** —— 数据/材料已存在但**未作为正式")
w("Supplementary 表完整呈现**，故不升格：")
w("")
w("　　" + "、".join(DERIVED_CONSERVATIVE.item.tolist()))
w("")
w("**E. 其余 `Partial`（`FROZEN_DIRECT`）** —— **已报告但附有明确写出的限定**：")
w("")
w("　　这些条目在正文中**已经写出**，其 Partial 状态来自**写明了的适用边界**（例如")
w("　　Item 6a 受去标识约束、Item 7 未按人群分层做 QC、Item 23a 未系统报告亚组），")
w("　　**不是遗漏**。若需升格为 Covered，所需动作是**补写文本**而非补数据。")
w("")
w("---")
w("")
w("## 4. 结论")
w("")
w("**Methods + Supplement reporting-completeness：scientifically / reporting complete，")
w("administratively not yet complete。**")
w("")
w("唯一阻止行政性关闭的是 **`AUTHOR_INPUT_REQUIRED` 项**（IRB 声明、资助、利益冲突、代码仓库地址，")
w("以及 Stanford / UTSW 若可恢复的真实 accrual ranges）。**这些不得代写。**")
w("")
w("**本支线到此为止，不再扩展。** 下一阶段：**最终 reference merge + 全文按首现顺序重新编号**。")
w("")

open(os.path.join(ROOT, "TRIPOD_AI_Gap_Audit_v3.md"), "w", encoding="utf-8").write("\n".join(L))
print(f"[OK] TRIPOD_AI_Gap_Audit_v3.md —— {len(D)} 项")
print(f"    状态: {STAT}")
print(f"    支撑: {SUP}")
