# -*- coding: utf-8 -*-
r"""p21e — 由冻结产物**自动渲染** Supplementary Table S2

================================ 为什么需要 ================================
正文 Results §7.2 小节引用了 "Supplementary Table S2"，但该表的**排版件不存在** ——
`figures/data/tableS2_metrics.csv` 只是数据件。正文引一个不存在的排版件是**引文完整性问题**，
不是可选项。

本脚本仿 `p21d`（Table R1–R3 的渲染器），规矩一致：
  · 从**冻结 CSV** 渲染，不手工录入；
  · **每个显示出来的 cell 一行 provenance**（来源文件 / 列 / 变换规则 / 显示值）；
  · 渲染件**不得包含**任何未登记的数值。

============================ 本脚本自带的对账 ============================
外部建议的两条便宜检查，都实现在此：

  C1 **行数 / 列头一致**：渲染出的行集与列头必须与冻结表逐项对应，多一个少一个即判红。
  C2 **逐 cell 回读**：把渲染出来的每个数字**解析回 float**，与冻结列比对（容差按显示精度）。
     这条能抓住格式化/舍入/列错位 —— 只比字符串是抓不住列错位的。

输出：
  figures/data/reporting_profiles/Supplementary_Table_S2.md
  figures/data/reporting_profiles/supplementary_tableS2_provenance.csv
用法：python p21e_render_tableS2.py
"""
import io, os, re, sys
import pandas as pd

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIG  = os.path.join(ROOT, "figures", "data")
D    = os.path.join(FIG, "reporting_profiles")
VER  = "p21e v1.0 (2026-10-07)"
os.makedirs(D, exist_ok=True)

T = pd.read_csv(os.path.join(FIG, "tableS2_metrics.csv"))
# ★ 2026-10-09（A5.10）：Panel C 的源 —— `p33_protocol72_metrics.csv`，与 Results §7.2 同源。
#   选它而不是 `p32_direction_convention.csv`：p33 在该口径下**四项协议指标齐全**（含区间），
#   而 p32 只有 pooled 一个量。两者在该量上逐位一致（0.362437 等），互为旁证。
_P33F = os.path.join(ROOT, "分析结果", "p33_protocol72_metrics.csv")
P72 = pd.read_csv(_P33F) if os.path.exists(_P33F) else None
ROWS = []
def rec(cid, item, ok, note=""):
    ROWS.append(dict(cid=cid, item=item, status="PASS" if ok else "FAIL", note=note))
    print("  [%-4s] %-56s %s" % ("PASS" if ok else "FAIL", item[:56], note), flush=True)
    return ok

cells = []       # 每显示一个数就记一行
def cell(table, row, col, src_col, xform, shown, numeric):
    cells.append(dict(table=table, method=row, column=col,
                      artifact_file="tableS2_metrics.csv", artifact_column=src_col,
                      transform_rule=xform, displayed_value=str(shown),
                      numeric_value=float(numeric)))
    return shown

# ---------------------------------------------------------------- 显示格式
def f3(x):
    return "" if pd.isna(x) else "%.3f" % float(x)
def f3ci(pt, lo, hi):
    if pd.isna(pt): return ""
    return "%.3f (%.3f–%.3f)" % (float(pt), float(lo), float(hi))
def f2s(x):
    return "" if pd.isna(x) else "%+.2f" % float(x)
def f2ci(pt, lo, hi):
    if pd.isna(pt): return ""
    return "%+.2f (%+.2f to %+.2f)" % (float(pt), float(lo), float(hi))
def f0(x):
    return "" if pd.isna(x) else "%d" % int(round(float(x)))

ORDER = ["domain_auc", "mmd", "psi_mean", "psi_max", "js_mean", "js_max",
         "pred_mean_shift", "entropy_mean_shift", "entropy_js", "conf_dist_js",
         "pred_pos_prev_shift", "SUDO_AURCC", "AUROC_DR", "AUROC_IOW", "AUROC_OM"]
FAMLAB = {"L1": "L1 input-space", "L2": "L2 prediction-space",
          "L3-primary": "L3 transported (primary)", "L3-sens": "L3 transported (sensitivity)"}

print("=" * 100); print("p21e — Supplementary Table S2 渲染"); print("=" * 100)

# ---------------------------------------------------------------- C1 冻结表身份
missing = [m for m in ORDER if m not in set(T.method)]
rec("C1-1", "冻结表中 15 个方法齐全（不多不少）",
    len(T) == 15 and not missing, f"行数={len(T)}；缺={missing or '无'}")

NEED_A = ["harmful_auc", "auc_lo", "auc_hi", "harmful_auprc", "auprc_lo", "auprc_hi",
          "sens_at_spec80", "sens80_lo", "sens80_hi",
          "ppv", "ppv_lo", "ppv_hi", "npv", "npv_lo", "npv_hi"]
NEED_B = ["role", "coverage", "spearman", "sp_lo", "sp_hi",
          "macro_logit", "n_evaluable_folds", "loo_org_auc"]
miss_cols = [c for c in NEED_A + NEED_B if c not in T.columns]
rec("C1-2", "渲染所需列全部存在于冻结表", not miss_cols, f"缺列={miss_cols or '无'}")

T = T.set_index("method").loc[ORDER].reset_index()
rec("C1-3", "渲染行序与登记顺序一致（protocol order，不按性能排）",
    list(T.method) == ORDER, "L1 六项 → L2 六项 → L3 三项")

# ---------------------------------------------------------------- 渲染
prev = None
md = ["# Supplementary Table S2 — Protocol-specified warning-task metrics and grouped summaries",
      "",
      "> **Task definition (protocol §7.2).** A transfer is *harmful* when ΔAUROC ≤ −0.05; the "
      "risk score *r*_failure is oriented so that a larger value denotes greater risk. Sensitivity "
      "is reported at the prespecified operating point of **0.80** specificity, implemented as the "
      "**smallest threshold whose realised specificity is at least 0.80** (Methods Section 14); "
      "with 129 non-failing transfers this is 104/129 = **0.806** for every signal, so every row "
      "carries the same FP = 25 and TN = 104 and only the true-positive count varies. Predictive "
      "values are reported at the observed harmful prevalence (57 of 186 transfers, 30.6%). "
      "Intervals are 95% organism–drug cluster bootstrap intervals.",
      "",
      "## Panel A — Warning-task discrimination (all 186 transfers)",
      "",
      "| Signal | Family | n | AUROC (95% CI) | AUPRC (95% CI) | Sens@0.80 (95% CI) | PPV (95% CI) | NPV (95% CI) |",
      "|---|---|---|---|---|---|---|---|"]
for _, r in T.iterrows():
    m = r["method"]; fam = FAMLAB.get(r["role"], r["role"])
    if m.startswith("AUROC_") and prev != "L3":
        pass
    a  = cell("S2A", m, "AUROC", "harmful_auc", "point + percentile CI", f3ci(r.harmful_auc, r.auc_lo, r.auc_hi), r.harmful_auc)
    ap = cell("S2A", m, "AUPRC", "harmful_auprc", "point + percentile CI", f3ci(r.harmful_auprc, r.auprc_lo, r.auprc_hi), r.harmful_auprc)
    sn = cell("S2A", m, "Sens@0.80", "sens_at_spec80", "point + percentile CI", f3ci(r.sens_at_spec80, r.sens80_lo, r.sens80_hi), r.sens_at_spec80)
    pv = cell("S2A", m, "PPV", "ppv", "point + percentile CI", f3ci(r.ppv, r.ppv_lo, r.ppv_hi), r.ppv)
    nv = cell("S2A", m, "NPV", "npv", "point + percentile CI", f3ci(r.npv, r.npv_lo, r.npv_hi), r.npv)
    cell("S2A", m, "coverage", "coverage", "as is", f0(r.coverage), r.coverage)
    md.append("| `%s` | %s | %s | %s | %s | %s | %s | %s |" % (m, fam, f0(r.coverage), a, ap, sn, pv, nv))

md += ["",
       "> **The two entropy- and confidence-distribution JS rows are near-identical by construction, not",
       "> by error.** Entropy and confidence are both deterministic functions of the same predicted",
       "> probability, so their source-to-target divergences track each other closely: Pearson",
       "> *r* = 0.993 (95% interval 0.972 to 0.997) across the 186 transfers (Table S2 rows 9–10).",
       # ★★ 2026-10-09 更正（A5.13/A5.14）：原注作「the 129 non-harmful transfers quantise to",
       #   **about 26** false positives, so two vectors with r ≈ 0.99 **can readily land** on the same",
       #   operating point」—— **两处都不准，且方向是把「必然」写成了「可能」。**
       #   实测（`p33_protocol72_metrics.csv`，15 primary ＋ 3 signed 共 18 行）：
       #   **FP 恒为 25、TN 恒为 104，无一例外；只有 TP 变**（6,7,8,9,10,13,17,21,22）。
       #   原因是纯算术，与信号排序无关：
       #     规则 = 满足 spec ≥ 0.80 的**最小**阈值；n_neg = 129
       #     spec = (129−FP)/129 ≥ 0.80  ⇒  FP ≤ 129×0.20 = 25.8  ⇒  FP ≤ 25（整数）
       #     灵敏度最大 ⟺ 阈值最低 ⟺ FP 取允许的最大值  ⇒  **FP ≡ 25，TN ≡ 104**
       #   ⇒ 这三个量的一切共同点**都是被规则强迫的**，不是巧合。唯一能变的是 TP。
       "> The identical Sens@0.80 / PPV / NPV entries are **forced by the operating-point rule, not",
       "> coincidental**: at 80% specificity the 129 non-harmful transfers admit at most",
       "> 129 × 0.20 = 25.8 false positives, and the rule takes the smallest threshold satisfying the",
       "> constraint, so every signal lands on exactly **FP = 25, TN = 104** regardless of its ranking.",
       "> Only TP can differ between signals — here 9 for both rows. Their AUROCs differ",
       "> (0.587 versus 0.585), which also rules out one column being read twice.",
       "",
       "## Panel B — Rank association and grouped summaries",
       "",
       "| Signal | Family | Spearman ρ vs ΔAUROC (95% CI) | LODO macro AUROC | LOO-organism AUROC | Evaluable folds |",
       "|---|---|---|---|---|---|"]
for _, r in T.iterrows():
    m = r["method"]; fam = FAMLAB.get(r["role"], r["role"])
    sp = cell("S2B", m, "Spearman", "spearman", "point + percentile CI", f2ci(r.spearman, r.sp_lo, r.sp_hi), r.spearman)
    # ★ 2026-10-09 更正（INC-019）：本列原渲染 `macro_mean`（来自 p6e_lodo_aggregation.csv），
    #   而 **Fig 3A 画的是 `macro_logit`**（= `C_macro_eval`，同族里唯一与 fig3a_lodo 逐位相同的列，
    #   15/15 bit-identical），**正文 §3 引的也是它**。两者本是同一个量（LODO macro，unweighted）
    #   被**两处代码各算一次**：9/15 行只差 1 ULP，**entropy_mean_shift 恰好跨在三位小数的进位边界上**
    #   ⇒ 表印 0.613、正文与图标 0.612，**读者一比即见**。
    #   改为读 `macro_logit`：图、表、正文从此同源。
    #   为什么选 logit 一侧：`fig3a` 是**点估计与 CI 同一口径**（p6f，unweighted）的自洽对，
    #   而 p6e 的 `macro_lo/hi` 是 size-weighted、与它自己的 unweighted 点估计口径不一致
    #   （EVS t13 EC-14 已把该不一致记为 WARN）。故图侧是应当被跟随的那一个。
    mc = cell("S2B", m, "LODO macro", "macro_logit", "as is", f3(r.macro_logit), r.macro_logit) if not pd.isna(r.macro_logit) else "—"
    lo = cell("S2B", m, "LOO-organism", "loo_org_auc", "as is", f3(r.loo_org_auc), r.loo_org_auc) if not pd.isna(r.loo_org_auc) else "—"
    nf = cell("S2B", m, "n folds", "n_evaluable_folds", "as is", f0(r.n_evaluable_folds), r.n_evaluable_folds) if not pd.isna(r.n_evaluable_folds) else "—"
    md.append("| `%s` | %s | %s | %s | %s | %s |" % (m, fam, sp, mc, lo, nf))

# ---------------------------------------------------------------- ★ Panel C（2026-10-09，A5.10）
#   病历：§3 的方向约定披露句引了三个 **signed 口径**的值（0.362 / 0.457 / 0.342），
#   而它们**在稿内没有任何读者可见的出处** —— 不在任何图、不在任何补充表、Table 1 里也没有。
#   由 `p54`（指针-值判据）发现，登记为 OPEN_GAP 挂在清单 A5.10 上。
#
#   裁决（外部核验 2026-10-09，理由采纳）：**扩 Table S2，不开新表。**
#     · 这三个是「同三个方法、**另一口径**的同一个量」—— 而 S2 Panel A 正是每方法数值的家；
#     · 为三个数立一张新表太重；且 S2 已在重打包同步清单里，扩面板不加件；
#     · **Table 1 不能承接**：它是**约定表**，内容靠定义成立、不需锚定；
#       而这三个是**结果数**，需要锚定。混排等于让一张表同时背两种担保等级。
#   数据源：`p33_protocol72_metrics.csv` 的 `convention == "signed"` 三行（AUROC 及其余四项齐全）。
SIG = P72[P72.convention == "signed"] if P72 is not None and "convention" in P72.columns else None
if SIG is not None and len(SIG):
    md += ["",
           "## Panel C — Direction-convention sensitivity (signed convention)",
           "",
           "> The magnitude convention applied to the three signed prediction-space shifts is a modelling",
           "> choice, not a definitional consequence (Methods Section 11). This panel reports the same",
           "> protocol §7.2 quantities for those three signals under the **untransformed (signed)**",
           "> values, so that the convention sensitivity stated in Results can be checked here.",
           "",
           "| Signal | Family | n | AUROC (95% CI) | AUPRC (95% CI) | Sens@0.80 (95% CI) | PPV (95% CI) | NPV (95% CI) |",
           "|---|---|---|---|---|---|---|---|"]
    for _, r in SIG.iterrows():
        m = r["method"]; fam = FAMLAB.get(r["role"], r["role"])
        md.append("| `%s` | %s | %s | %s | %s | %s | %s | %s |" % (
            m, fam, f0(r["n"]),
            f3ci(r.auroc, r.auroc_lo, r.auroc_hi),
            f3ci(r.auprc, r.auprc_lo, r.auprc_hi),
            f3ci(r.sens80, r.sens80_lo, r.sens80_hi),
            f3ci(r.ppv, r.ppv_lo, r.ppv_hi),
            f3ci(r.npv, r.npv_lo, r.npv_hi)))
    for _, r in SIG.iterrows():
        for _c in ("auroc", "auprc", "sens80", "ppv", "npv"):
            cell("S2C", r["method"], _c, _c, "signed convention", f3(r[_c]), r[_c])
    # ★ 2026-10-09（外部核验建议）：三行共用同一组敏感度/PPV/NPV —— 审稿人必问「为什么一模一样」。
    #   预先拆掉：**因为它们是同一个混淆矩阵**。计数已从 p33 的 tp/fp/fn/tn 列实测核对：
    #   TP=9, FN=48, FP=25, TN=104，合计 186 ⇒ sens=9/57、PPV=9/34、NPV=104/152，三行必然相等。
    #   只陈述**事实与算术**，不替它给机理解释（那是尚待验证的主张）。
    _q = SIG.iloc[0]
    md += ["",
           "> **The three rows share a single confusion matrix, and this is forced rather than",
           "> coincidental.** The prespecified rule is the smallest threshold achieving specificity",
           "> ≥ 0.80; with 129 non-failing transfers that requires FP ≤ 129 × 0.20 = 25.8, so",
           "> **FP = 25 and TN = 104 for every signal, independently of its ranking**. The counts here are",
           "> TP %d, FN %d, FP %d, TN %d, summing to 186, giving sensitivity 9/57, PPV 9/34 and" % (
               int(_q.tp), int(_q.fn), int(_q.fp), int(_q.tn)),
           "> NPV 104/152. Only TP is free to vary between signals, which is why AUROC and AUPRC — both",
           "> rank-based — still differ across the three rows while the three threshold-dependent",
           "> quantities do not. Counts are read from the `tp`/`fp`/`fn`/`tn` columns of",
           "> the frozen protocol-metrics table on the `signed` rows."]
    md += ["",
           "> All five quantities for each row, and the 95% intervals, come from "
           "the frozen protocol-metrics table, `signed` convention. The corresponding",
           "> primary-convention values for the same three signals are in Panel A; comparing the two",
           "> panels shows the direction convention changes the ordering, which is the point of the panel.",
           ""]

md += ["",
       "> **Negative predictive values are high for every method, including those that discriminate no",
       "> better than chance**, because 69.4% of transfers are not harmful; a high NPV should not be read",
       "> as evidence that a signal works. Sensitivity at the prespecified operating point is at most",
       "> 0.386, so even the best method identifies a minority of failing transfers.",
       ""]

io.open(os.path.join(D, "Supplementary_Table_S2.md"), "w", encoding="utf-8").write("\n".join(md))
pdf = pd.DataFrame(cells)
pdf.to_csv(os.path.join(D, "supplementary_tableS2_provenance.csv"), index=False, encoding="utf-8-sig")

# ---------------------------------------------------------------- C2 逐 cell 回读
print()
# Panel A 记录 6 个数值列（AUROC/AUPRC/Sens/PPV/NPV/coverage）；
#   Panel B 记录 4 个（Spearman/LODO macro/LOO-organism/n folds）。
#   family 与行标签是文字，不记 cell —— 计数要按**实际记录的数值列**算。
rec("C2-1", "Panel A 6 个数值列 × 15 行", len([c for c in cells if c["table"] == "S2A"]) == 15 * 6,
    f"S2A cell 数 = {len([c for c in cells if c['table']=='S2A'])}（期望 15×6=90）")
rec("C2-2b", "Panel C 5 个数值列 × 3 行（signed 口径）", len([c for c in cells if c["table"] == "S2C"]) == 3 * 5,
    f"S2C cell 数 = {len([c for c in cells if c['table']=='S2C'])}（期望 3×5=15）")
rec("C2-2", "Panel B 4 个数值列 × 15 行", len([c for c in cells if c["table"] == "S2B"]) == 15 * 4,
    f"S2B cell 数 = {len([c for c in cells if c['table']=='S2B'])}（期望 15×4=60）")
rec("C2-3", "provenance 行数 == 显示 cell 数", len(pdf) == len(cells), f"{len(pdf)} 行")
# ★★ 2026-10-09（外部核验指出）：Panel C 的**表注是一段含 9/48/25/104/186 的散文**，
#   而 `p54` 的句扫描只覆盖正文四节＋摘要 ⇒ **这段散文的算术此前只享受了人工核对**。
#   「人工扫描可以是发现手段，不能是判据」（INC-021）。故在此把它**逐值对源**：
#   回读渲染件，要求表注里那四个计数与 `p33` signed 行的 tp/fp/fn/tn **逐值一致**。
_S2TXT = io.open(os.path.join(D, "Supplementary_Table_S2.md"), encoding="utf-8").read()
if SIG is not None and len(SIG):
    _q = SIG.iloc[0]
    _need = "TP %d, FN %d, FP %d, TN %d" % (int(_q.tp), int(_q.fn), int(_q.fp), int(_q.tn))
    # ★ 编号更正：本断言初版误用了 **C2-5** —— 该 id 已被「列头与登记列集合一致」占用。
    #   **同一张审计表里两个不同主张共用一个 cid**，会让「按 cid 引用某条断言」变得有歧义。
    #   故改为 C2-6，并在下方加 **cid 唯一性断言**（本仓第一次出现该形态，属预防）。
    # ★ 2026-10-09：把「FP ≡ 25、TN ≡ 104 与排序无关」**做成断言**（原为表注里的散文主张）。
    #   规则推演：FP ≤ 129×0.20 = 25.8 ⇒ FP ≤ 25；取最小阈值 ⇒ FP 取最大允许值 = 25。
    _all = P72.copy()
    _bad = _all[(_all.fp != 25) | (_all.tn != 104)]
    rec("C2-8", "操作点规则的结构性后果：**所有**行 FP=25、TN=104（与信号排序无关）",
        len(_bad) == 0,
        f"共 {len(_all)} 行（primary+signed）；违例 {len(_bad)} —— 规则：FP ≤ 129×0.20 = 25.8 ⇒ FP=25")
    # ★ 2026-10-09（A5.13/A5.14b）：**操作点规则**现已写进表注。断言渲染件真的印了规则与实现值，
    #   免得「规则写进了生成器、渲染件却没跟上」——这正是 INC-021 的形状（判据的输入与结论分家）。
    rec("C2-9", "任务定义表注印出操作点规则与实现特异度（104/129 = 0.806）",
        ("smallest threshold whose realised specificity is at least 0.80" in _S2TXT
         and "104/129" in _S2TXT),
        "Methods §14 与表注须同述一条规则；实现值 104/129 由 C2-8 结构性保证")
    rec("C2-6", "Panel C 表注的混淆矩阵计数：渲染件回读 == p33 signed 行",
        _need in _S2TXT, f"要求渲染件含「{_need}」（源 p33_protocol72_metrics.csv 的 tp/fp/fn/tn）")

# 回读：把显示串里的**第一个数字**解析回来，与**按显示精度舍入后**的冻结值比。
#   ★ 容差必须由显示串的小数位数推导，不能用固定值 —— 首版用 5e-4，
#     于是 "+0.13"（两位小数）对 0.12628 被判不符，14 项假阳性。
bad = []
for c in cells:
    if c["displayed_value"] == "—":
        continue
    mm = re.search(r"[-+]?\d*\.?\d+", c["displayed_value"])
    if not mm:
        bad.append((c["method"], c["column"], c["displayed_value"], "无可解析数字"))
        continue
    tok = mm.group()
    dec = len(tok.split(".")[1]) if "." in tok else 0
    if abs(float(tok) - round(c["numeric_value"], dec)) > 1e-9:
        bad.append((c["method"], c["column"], tok, c["numeric_value"]))
rec("C2-4", "每个显示值回读后与冻结数值一致（按显示精度）", not bad,
    f"不符 {len(bad)} 项" + (f"：{bad[:3]}" if bad else ""))

# 列头一致性：渲染件里出现的列名必须都在登记列集合内
hdr_ok = all(any(h in ln for ln in md) for h in
             ["AUROC (95% CI)", "AUPRC (95% CI)", "Sens@0.80 (95% CI)", "PPV (95% CI)",
              "NPV (95% CI)", "LODO macro AUROC", "LOO-organism AUROC"])
rec("C2-5", "渲染件的列头与登记列集合一致", hdr_ok, "")

print()
# ★ 2026-10-09：**cid 唯一性**。本文件新增 C2-6 时曾误用已被占用的 C2-5，
#   造成同一张审计表里两个不同主张共用一个编号 —— 「按 cid 引用某条断言」即失效。
#   属 INC-021 族在**判据自身编号空间**上的形态。此处加断言预防。
#   注意：必须在 `pd.DataFrame(ROWS)` **之前**调用 —— 否则新行进 ROWS 却进不了已建好的 Dd。
_dupcid = pd.Series([r["cid"] for r in ROWS])
_d = _dupcid[_dupcid.duplicated()].unique().tolist()
rec("C2-7", "审计表内 cid 唯一（无两条主张共用编号）", not _d,
    "重复 cid = " + (", ".join(_d) if _d else "无"))
Dd = pd.DataFrame(ROWS)
nbad = int((Dd.status != "PASS").sum())
print("=" * 100)
print(f"p21e 断言合计 {len(Dd)} 项：PASS {int((Dd.status=='PASS').sum())}｜**非 PASS {nbad}**")
if nbad:
    print("*** p21e 未通过 ***")
else:
    print("*** p21e 通过 —— Table S2 已渲染，且行/列头与逐 cell 回读全部对账一致 ***")
Dd.to_csv(os.path.join(D, "supplementary_tableS2_render_audit.csv"), index=False, encoding="utf-8-sig")

# ================================================================ Table S1
# 2026-10-07 增补（外部建议 + p35 清点）：论文的**分析单位总体**是这 186 个迁移，
# 而正文此前无处指向它们的清单。渲染 S1 并在 Methods §6 引用一次。
print()
print("-" * 100)
print("p21e — Supplementary Table S1（完整任务清单）渲染")
print("-" * 100)
T1 = pd.read_csv(os.path.join(FIG, "tableS1_tasks.csv"))
S1_ROWS = []
def rec1(cid, item, ok, note=""):
    S1_ROWS.append(dict(cid=cid, item=item, status="PASS" if ok else "FAIL", note=note))
    print("  [%-4s] %-56s %s" % ("PASS" if ok else "FAIL", item[:56], note), flush=True)

U = pd.read_csv(os.path.join(FIG, "fig1b_universe.csv")).iloc[0]
rec1("T1-1", "行数 == 合格宇宙（fig1b.eligible_tasks）", len(T1) == int(U.eligible_tasks),
     f"{len(T1)} vs {int(U.eligible_tasks)}")
rec1("T1-2", "task_id 唯一", T1.task_id.is_unique, f"{T1.task_id.nunique()} 个唯一 id")
rec1("T1-3", "harmful 计数 == fig1c.n_harmful",
     int(T1.harmful.sum()) == int(pd.read_csv(os.path.join(FIG, "fig1c_summary.csv")).n_harmful.iat[0]),
     f"{int(T1.harmful.sum())}")
rec1("T1-4", "ΔAUROC 恒等式 d_auroc == auroc_target − auroc_source_internal",
     bool((T1.d_auroc - (T1.auroc_target - T1.auroc_source_internal)).abs().max() < 1e-12),
     "逐行核对")
rec1("T1-5", "harmful 标志 == (d_auroc <= -0.05)", bool((T1.harmful.astype(bool) == (T1.d_auroc <= -0.05)).all()))
rec1("T1-6", "无缺失值", int(T1.isna().sum().sum()) == 0, f"{int(T1.isna().sum().sum())} 个 NaN")

# ★ 渲染**全部列**（21 列），不做编辑性收窄。
#   首版只渲染身份+结局 11 列，把 12 个信号列留在数据文件里 —— 那是我做的编辑判断，
#   会让「完整清单」名不副实，且读者无法从渲染件核对信号取值。
S1_cols = list(T1.columns)
# ★★ 2026-10-07：渲染件是**投稿材料**，正文与表注必须**全英文**。
#   首版把内部工作注记（"自动渲染，不手工录入"、脚本文件名、版本号）一并写进了表头，
#   并沿用了中文表注 —— 这些都会原样进到投稿的 Word 里。
#   处置：期刊件的表头改为英文陈述；**内部溯源不写在表里**（它已经在
#   `supplementary_tableS1_render_audit.csv` 与 `supplementary_tableS2_provenance.csv` 里）。
# ★ 2026-10-07（A2.1）：列数 21 → **26**（补齐 5 个信号）。断言随之更新。
HDR = {"task_id": "task_id", "organism": "organism", "drug": "drug",
       "source_site": "source", "target_site": "target", "n_target": "n target",
       "auroc_source_internal": "AUROC source", "auroc_target": "AUROC target",
       "d_auroc": "ΔAUROC", "harmful": "harmful", "stable": "stable",
       "domain_auc": "domain_auc", "mmd": "MMD²", "psi_mean": "PSI mean",
       "psi_max": "PSI max", "js_mean": "JS mean", "js_max": "JS max",
       "pred_mean_shift": "pred. mean shift", "entropy_mean_shift": "entropy mean shift",
       "entropy_js": "entropy JS", "conf_dist_js": "confidence JS",
       "pred_pos_prev_shift": "pred. prevalence shift", "SUDO_AURCC": "SUDO AURCC",
       "AUROC_DR": "ÂUROC DR", "AUROC_IOW": "ÂUROC IOW", "AUROC_OM": "ÂUROC OM"}
assert len(S1_cols) == 26, "tableS1 列数变了（%d），HDR 需同步" % len(S1_cols)
sig_cols = [c for c in T1.columns if c not in
            ("task_id", "organism", "drug", "source_site", "target_site", "n_target",
             "auroc_source_internal", "auroc_target", "d_auroc", "harmful", "stable")]
md1 = ["# Supplementary Table S1 — The 186 eligible directed cross-system transfers",
       "",
       # ★ 2026-10-07（A2.2）：首句「All derived quantities … computed on these rows」
       #   在补齐 5 个信号列**之前为假**（composite 的成员信号不在表里）。
       #   现已为真，并**指出机械担保**：p49 从这 26 列复算 composite 的 7/37 与 S2 Panel A 五行。
       "> **This table is the analysis-unit population of the study.** It lists all 186 eligible "
       "directed cross-system transfers — the complete benchmark universe — with the source and "
       "target discrimination of each. All derived quantities in this manuscript are computed on "
       "these rows, and that claim is mechanically checked: the triage contrast of the "
       "retrospective policy simulation and the five signal rows of Supplementary Table S2 "
       "Panel A are reproduced bit-for-bit "
       "from the columns below.",
       # ★ 2026-10-08：**期刊件不得含脚本文件名**（p41 V7-3；投稿材料只描述可核验性，不报内部实现）。校验器的名字留在仓库，不进补充表。
       "",
       f"> `task_id` has the form `<organism>_<drug>_<source>2<target>`. **harmful** denotes "
       f"ΔAUROC ≤ −0.05, the prespecified deployment endpoint; **stable** denotes |ΔAUROC| ≤ 0.02, "
       f"a descriptive near-zero band that is not a second endpoint.",
       # ★ 2026-10-07（A2.2）：原句两处不实，均已改正 ——
       #   ① 「All derived quantities … computed on these rows」：补齐 5 个信号列**之后**才为真；
       #      此前 composite（主对比的对照臂）所需的 psi_max 与 js_max 不在表里。
       #   ② 「in the orientation in which each signal enters the evaluation」：实测**不成立** ——
       #      表内各列一律是**原始值**（pn 未取绝对值、AURCC 与 L3 未取负）。
       #      评价时的方向约定在 Methods §11 与正文 Table 1，读者须自行应用。
       f"> All {len(S1_cols)} columns are shown, including the per-transfer value of each of the "
       f"{len(sig_cols)} label-free monitoring signals. **The signal columns hold the raw "
       f"computed values, not the oriented risk scores**: the per-signal orientation applied at "
       f"evaluation (identity, absolute magnitude, or negation) is specified in Methods §11 and "
       f"main-text Table 1, and must be applied by the reader before reproducing any aggregate. "
       f"Aggregated signal performance is reported in Supplementary Table S2.",
       "",
       "| " + " | ".join(HDR[c] for c in S1_cols) + " |",
       "|" + "---|" * len(S1_cols)]
for _, r in T1.iterrows():
    cells1 = []
    for c in S1_cols:
        v = r[c]
        if c in ("auroc_source_internal", "auroc_target", "d_auroc"):
            cells1.append("%.3f" % float(v))
        elif c in ("harmful", "stable"):
            cells1.append("yes" if bool(v) else "")
        elif c == "n_target":
            cells1.append("%d" % int(v))
        elif isinstance(v, (int, float)) and not isinstance(v, bool):
            cells1.append("%.4f" % float(v))          # 全部信号：统一 4 位小数
        else:
            cells1.append(str(v))
    md1.append("| " + " | ".join(cells1) + " |")
md1 += [""]
io.open(os.path.join(D, "Supplementary_Table_S1.md"), "w", encoding="utf-8").write("\n".join(md1))

rnd1 = io.open(os.path.join(D, "Supplementary_Table_S1.md"), encoding="utf-8").read()
data_rows = [l for l in rnd1.splitlines() if l.startswith("| ") and "---" not in l][1:]
rec1("T1-7", "渲染行数 == 冻结表行数", len(data_rows) == len(T1), f"{len(data_rows)} vs {len(T1)}")
rec1("T1-8", "渲染列数 == 登记列数", rnd1.splitlines()[0].count("|") == 0 and
     all(h in rnd1 for h in HDR.values()), "列头齐全")
bad1 = []
for _, r in T1.iterrows():
    ln = [l for l in data_rows if l.startswith("| %s |" % r["task_id"])]
    if not ln:
        bad1.append(r["task_id"]); continue
    cells = [c.strip() for c in ln[0].strip("|").split("|")]
    if len(cells) != len(S1_cols):
        bad1.append("%s:列数%d" % (r["task_id"], len(cells))); continue
    for c, shown in zip(S1_cols, cells):
        v = r[c]
        if c in ("harmful", "stable"):
            exp = "yes" if bool(v) else ""
        elif c == "n_target":
            exp = "%d" % int(v)
        elif c in ("auroc_source_internal", "auroc_target", "d_auroc"):
            exp = "%.3f" % float(v)
        elif isinstance(v, (int, float)) and not isinstance(v, bool):
            exp = "%.4f" % float(v)
        else:
            exp = str(v)
        if shown != exp:
            bad1.append("%s:%s(%s≠%s)" % (r["task_id"], c, shown, exp))
rec1("T1-9", f"逐格回读：全部 {len(S1_cols)} 列 × 186 行与冻结表一致", not bad1,
     f"不符 {len(bad1)} 项" + (f"：{bad1[:3]}" if bad1 else ""))

print()
D1 = pd.DataFrame(S1_ROWS); n1 = int((D1.status != "PASS").sum())
print(f"p21e-S1 断言合计 {len(D1)} 项：PASS {int((D1.status=='PASS').sum())}｜**非 PASS {n1}**")
D1.to_csv(os.path.join(D, "supplementary_tableS1_render_audit.csv"), index=False, encoding="utf-8-sig")

sys.exit(1 if (nbad or n1) else 0)
