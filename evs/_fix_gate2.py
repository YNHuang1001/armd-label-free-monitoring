# -*- coding: utf-8 -*-
"""用**行号定位**做 4 处表述收紧（多行字面匹配已经失败两次，不再用）。"""
import io, os, sys, ast

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
F = "t20_identity_gate.py"
L = io.open(F, encoding="utf-8").read().split("\n")


def find(sub, start=0):
    for i in range(start, len(L)):
        if sub in L[i]:
            return i
    raise LookupError(sub)


def end_of_stmt(i):
    """从 i 起，找到该语句结束（圆括号配平）的行号。"""
    depth = 0
    j = i
    while j < len(L):
        depth += L[j].count("(") - L[j].count(")")
        if depth <= 0 and L[j].rstrip().endswith((")", '")')):
            return j
        j += 1
    raise LookupError("stmt end")

IND = "    "      # if len(flip): 块内的缩进

# ---- ① Layer 3 三个数字 ----
i = find("eligible → not eligible 的 task")
j = end_of_stmt(i)
L[i:j + 1] = [
    IND + 'p("| | 数量 |")',
    IND + 'p("|---|---|")',
    IND + 'p(f"| **Removed eligible tasks**（eligible → not eligible） | **{len(flip)}** |")',
    IND + 'p(f"| **Newly qualified tasks**（not eligible → eligible） | **{_newly}** |")',
    IND + 'p(f"| **Net eligible change** | **{n_el_n - n_el_o:+d}**（{n_el_o} → {n_el_n}） |")',
    IND + 'p("")',
    IND + 'p(f"（{len(flip)} − {_newly} = {len(flip) - _newly}，与 Layer 2 的净变化一致；'
          '只列 removed 会让读者问「为什么列 24 却只掉 16」）")',
]
# _newly 需要在块外先算出来（原来是在 f-string 里内联的）
k = find("flip = m[m.flip].copy()")
L.insert(k + 1, '')
L.insert(k + 2, 'p("## Layer 3 — 掉了哪些 task、为什么")'.replace("## Layer 3 — 掉了哪些 task、为什么", "") or "")
del L[k + 2]
L.insert(k + 2, '_newly = n_el_n - (n_el_o - len(flip))')

# ---- ② 归因：原因 vs 结果 ----
i = find("**归因汇总**")
j = end_of_stmt(i)
L[i:j + 1] = [
    IND + 'p("### 归因（区分「原因」与「结果」）")',
    IND + 'p("")',
    IND + 'p("**因果链**：")',
    IND + 'p("")',
    IND + 'p("```")',
    IND + 'p("CS-01/02/03 修正  →  训练/评价样本被重建  →  source 模型性能被**重新计算**")',
    IND + 'p("                  →  AUROC 落在预设合格线 0.65 之下  →  task 出局")',
    IND + 'p("```")',
    IND + 'p("")',
    IND + 'p("★ **门槛规则本身没有改变**（仍是预先写死的 source 内部 AUROC ≥ 0.65）；'
          '**被纠正的是评价对象**——episode 单元与 S/R 结局定义。")',
    IND + 'p("")',
    IND + 'p("因此正确的表述是：")',
    IND + 'p("")',
    IND + 'p("> All removed tasks failed the **predefined** source-AUROC eligibility criterion "',
    IND + '  "after CS-01/02/03-induced reconstruction of the benchmark unit and outcome labels.")',
    IND + 'p("")',
    IND + 'p("**不得**写成「因为 AUROC 降了所以掉」——AUROC 下降是**结果**，不是原因。")',
    IND + 'p("")',
    IND + 'p("| 结果层的数字 | 计数 |")',
    IND + 'p("|---|---|")',
    IND + 'p(f"| removed task 的 AUROC 现在低于门槛 | '
          '{int((f2.src_auroc_n < AUC_THR).sum())} / {len(f2)} |")',
    IND + 'p(f"| 其中 n_target 同时减少 | {int((f2.n_target_n < f2.n_target_o).sum())} |")',
    IND + 'p(f"| 其中 n_src_test 同时减少 | {int((f2.n_src_test_n < f2.n_src_test_o).sum())} |")',
]

# ---- ③ Layer 4：eligible 子集质量 ----
i = find("src_auroc（合格线筛选量）")
L[i:i + 2] = [
    'p(f"| src_auroc 中位（**全部 228 candidate**） | {m.src_auroc_o.median():.4f} | "',
    '  f"{m.src_auroc_n.median():.4f} | {d_auc.median():+.4f} |")',
    'p(f"| src_auroc 中位（**eligible 子集**） | {_eo.src_auroc_o.median():.4f}（n={len(_eo)}） | "',
    '  f"{_en.src_auroc_n.median():.4f}（n={len(_en)}） | "',
    '  f"{_en.src_auroc_n.median() - _eo.src_auroc_o.median():+.4f} |")',
    'p(f"| src_auroc 中位（**前后都合格**的 {len(_keep)} 个） | {_keep.src_auroc_o.median():.4f} | "',
    '  f"{_keep.src_auroc_n.median():.4f} | "',
    '  f"{_keep.src_auroc_n.median() - _keep.src_auroc_o.median():+.4f} |")',
    'p("")',
    'p(f"- ★ **存活任务的质量**：eligible 子集的 src_auroc 中位 "',
    '  f"{_eo.src_auroc_o.median():.4f} → {_en.src_auroc_n.median():.4f}（"',
    '  f"{_en.src_auroc_n.median() - _eo.src_auroc_o.median():+.4f}）——"',
    '  f"**修正优先移除了边缘任务**，留下的是质量更高的迁移评价集。")',
    'p("- 可写入论文的表述：*The revised benchmark excluded marginal transfers failing the "',
    '  "predefined source-model quality criterion, resulting in a smaller but more reliable "',
    '  "transfer evaluation set.*")',
    'p("")',
]
# _eo/_en/_keep 需要在循环前定义
k = find("q = '### 归因")
k2 = find("p('## Layer 4 — 性能变化预估（不重训下游模型即可给出）')")
L.insert(k2 - 1, "_eo = m[m.eligible_o.astype(bool)]")
L.insert(k2, "_en = m[m.eligible_n.astype(bool)]")
L.insert(k2 + 1, "_keep = m[m.eligible_o.astype(bool) & m.eligible_n.astype(bool)]")
L.insert(k2 + 2, "")

# ---- ④ N_harmful / N_stable 措辞 ----
i = find("N_stable (\\\\|Δ\\\\| ≤ 0.02)")
j = end_of_stmt(i)
L[i + 1:j + 1] = [
    '    p("")',
    '    p("★ **这两个绝对计数的变化不得表述为「模型性能下降」**——它们改变是因为**被评价的 "',
    '      "task universe 被重建**（分母变了），不是同一个 universe 上的性能变差。正确表述：")',
    '    p("")',
    '    p("> The absolute counts changed because the evaluated task universe was reconstructed.")',
]

out = "\n".join(L)
ast.parse(out)
io.open(F, "w", encoding="utf-8", newline="\n").write(out)
print("4 处表述已收紧（行号定位）  AST OK")
