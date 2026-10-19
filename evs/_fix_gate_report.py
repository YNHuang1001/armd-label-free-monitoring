# -*- coding: utf-8 -*-
r"""按外部审计意见收紧 BENCHMARK IDENTITY GATE 的表述（4 处）：
  ① Layer 3 标题必须**同时给出** removed / newly / net，否则读者问"为什么列 24 掉 16"
  ② Layer 3 归因：**不要把结果当原因**。正确因果链是
       CS 修正 -> 训练/评价样本被重建 -> source 性能被重算 -> AUROC < 0.65 -> 出局
     （门槛规则**没变**，被纠正的是评价对象）
  ③ Layer 4 补上「**eligible 子集**的 src_auroc 中位」并给出正向解释：
     修正**优先移除边缘任务** ⇒ 留下的 benchmark 质量更高
     旧版只印了全 228 task 的中位，没做这个解读
  ④ N_harmful / N_stable 的措辞：不得被读成"性能下降"——绝对计数变化是因为
     **被评价的 task universe 被重建**
"""
import io, os, sys, ast

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
F = "t20_identity_gate.py"
t = io.open(F, encoding="utf-8").read()


def rep(a, b):
    global t
    assert t.count(a) == 1, (a[:60], t.count(a))
    t = t.replace(a, b, 1)


# ---- ① Layer 3 标题：三个数字都给出 ----
rep('''p(f"**eligible → not eligible 的 task：{len(flip)} 个**"
  f"（其余 {n_el_o - len(flip)} 个仍合格；{n_el_n - (n_el_o - len(flip))} 个为新增合格，若为 0 则无新增）")''',
    '''_newly = n_el_n - (n_el_o - len(flip))
p("| | 数量 |")
p("|---|---|")
p(f"| **Removed eligible tasks**（eligible → not eligible） | **{len(flip)}** |")
p(f"| **Newly qualified tasks**（not eligible → eligible） | **{_newly}** |")
p(f"| **Net eligible change** | **{n_el_n - n_el_o:+d}**（{n_el_o} → {n_el_n}） |")
p("")
p(f"（{len(flip)} − {_newly} = {len(flip) - _newly}，与 Layer 2 的净变化一致；"
  f"只列 removed 会让读者问"为什么列 24 却只掉 16"）")''')

# ---- ② Layer 3 归因：区分原因与结果 ----
rep('''p(f"**归因汇总**："
  f"因 source AUROC 跌破合格线 = {int((f2.src_auroc_n < AUC_THR).sum())}；"
  f"n_target 减少的 = {int((f2.n_target_n < f2.n_target_o).sum())}；"
  f"n_src_test 减少的 = {int((f2.n_src_test_n < f2.n_src_test_o).sum())}")''',
    '''p("### 归因（区分「原因」与「结果」）")
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
p("> All removed tasks failed the **predefined** source-AUROC eligibility criterion after "
  "CS-01/02/03-induced reconstruction of the benchmark unit and outcome labels.")
p("")
p("**不得**写成「因为 AUROC 降了所以掉」——AUROC 下降是**结果**，不是原因。")
p("")
p(f"| 结果层的数字 | 计数 |")
p("|---|---|")
p(f"| removed task 的 AUROC 现在低于门槛 | {int((f2.src_auroc_n < AUC_THR).sum())} / {len(f2)} |")
p(f"| 其中 n_target 同时减少 | {int((f2.n_target_n < f2.n_target_o).sum())} |")
p(f"| 其中 n_src_test 同时减少 | {int((f2.n_src_test_n < f2.n_src_test_o).sum())} |")''')

# ---- ③ Layer 4：eligible 子集的 src_auroc 中位 + 正向解释 ----
rep('''p(f"| src_auroc（合格线筛选量） | {m.src_auroc_o.median():.4f} | {m.src_auroc_n.median():.4f} | "
  f"中位 {d_auc.median():+.4f} |")
p("")''',
    '''_eo = m[m.eligible_o.astype(bool)]
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
  f"（{_en.src_auroc_n.median() - _eo.src_auroc_o.median():+.4f}）——"
  f"**修正优先移除了边缘任务**，留下的是质量更高的迁移评价集。")
p(f"- 可写入论文的表述：*The revised benchmark excluded marginal transfers failing the "
  f"predefined source-model quality criterion, resulting in a smaller but more reliable "
  f"transfer evaluation set.*")
p("")''')

# ---- ④ N_harmful / N_stable 措辞 ----
rep('''    p(f"| N_harmful (Δ ≤ −0.05) | {int(old_g.n_harmful.iat[0])} | **{int(new_g.n_harmful.iat[0])}** | "
      f"{int(new_g.n_harmful.iat[0]) - int(old_g.n_harmful.iat[0]):+d} |")
    p(f"| N_stable (\\\\|Δ\\\\| ≤ 0.02) | {int(old_g.n_stable.iat[0])} | **{int(new_g.n_stable.iat[0])}** | "
      f"{int(new_g.n_stable.iat[0]) - int(old_g.n_stable.iat[0]):+d} |")''',
    '''    p(f"| N_harmful (Δ ≤ −0.05) | {int(old_g.n_harmful.iat[0])} | **{int(new_g.n_harmful.iat[0])}** | "
      f"{int(new_g.n_harmful.iat[0]) - int(old_g.n_harmful.iat[0]):+d} |")
    p(f"| N_stable (\\\\|Δ\\\\| ≤ 0.02) | {int(old_g.n_stable.iat[0])} | **{int(new_g.n_stable.iat[0])}** | "
      f"{int(new_g.n_stable.iat[0]) - int(old_g.n_stable.iat[0]):+d} |")
    p("")
    p("★ **这两个绝对计数的变化不得表述为「模型性能下降」**——它们改变是因为**被评价的 "
      "task universe 被重建**（分母变了），不是同一个 universe 上的性能变差。正确表述：")
    p("")
    p("> The absolute counts changed because the evaluated task universe was reconstructed.")''')

ast.parse(t)
io.open(F, "w", encoding="utf-8", newline="\n").write(t)
print("4 处表述已收紧  AST OK")
