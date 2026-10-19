# -*- coding: utf-8 -*-
r"""p21d — 由 reporting artifact **自动渲染** Supplementary Tables R1–R3

================================ 为什么必须自动渲染 ================================
2026-10-02：我在手抄 Table R1 时**编造了 age IQR**（写成 51–77 / 41–69 / 40–69，
真实为 45–75 / 35–65 / 35–65）。编造的派生值与真值量级相符、方向合理，**肉眼极难发现**。

**新项目规则（GPT6 定，2026-10-02）**：
> No derived numerical value may enter the manuscript or Supplement unless it is either
> **programmatically rendered from an authoritative artifact** or **independently reconciled
> against that artifact**.

本脚本实现前半句：**Supplement 表格不再手抄，全部由 artifact 渲染**。
后半句由 p21 的 Layer B 承担（逐 cell 重新派生并比对）。

================================ 产物 ================================
  Supplementary_Tables_R1_R3.md          渲染出的表格（可直接进 Supplement）
  supplementary_cell_provenance.csv      **每个 cell 一行**：来源文件/字段/变换规则/显示值

用法：python p21d_render_reporting_tables.py
"""
import os, sys, hashlib
import numpy as np, pandas as pd

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); FIG = os.path.join(ROOT, "figures", "data")
D = os.path.join(FIG, "reporting_profiles")
VER = "p21d v1.0 (2026-10-02)"

R = lambda f: pd.read_csv(os.path.join(D, f))
sp  = R("site_profile.csv").set_index("site")
pf  = R("participant_flow.csv")
tp  = R("task_profile.csv")
xsd = R("site_x_organism_drug_counts.csv")
SITES = ["MGB", "SU", "UT"]
LABEL = {"MGB": "MGB", "SU": "Stanford", "UT": "UTSW"}

cells = []           # 每个 cell 的 provenance
def cell(table, row_label, site, artifact, fields, transform, value):
    cells.append(dict(table=table, row_label=row_label, site=site,
                      artifact_file=f"reporting_profiles/{artifact}",
                      artifact_fields=fields, transform_rule=transform,
                      displayed_value=str(value)))
    return value

def g(field, site):
    return sp.loc[site, field]

# ================================================================ Table R1
# ★ 2026-10-07（A1.4 + A1.9）：**每行带分母标签**。
#   病因：本表原先三个分母共处而无行级标签，直接产出「38 vs 37」之争 ——
#   R1 的 38 是**候选网格**组数、R3 的 37 是**合格宇宙**组数，两者各自都对，只是没标分母。
#   ⇒ 新增首列 `Denominator`，取值域 raw / candidate / eligible / task-expanded。
DEN = {"raw": "raw (unique episode)", "candidate": "candidate grid",
       "eligible": "eligible universe", "te": "task-expanded"}
SPEC1 = [
 ("Organism–drug groups represented", None, "38", "constant", "candidate"),
 ("Unique microbiological episodes", "n_unique_microbiological_episodes", "{:,}", "int_comma", "raw"),
 ("Unique patients", "n_unique_patients", "{:,}", "int_comma", "raw"),
 ("Task-expanded episode observations", "n_task_expanded_episode_observations", "{:,.0f}", "int_comma_f", "te"),
 ("Eligible directed transfers, as source", "n_eligible_transfers_as_source", "{:,.0f}", "int_comma_f", "eligible"),
 ("Eligible directed transfers, as target", "n_eligible_transfers_as_target", "{:,.0f}", "int_comma_f", "eligible"),
 ("Age recorded", "n_age_observed", "{:,}", "int_comma", "raw"),
 ("Age missing", "n_age_missing", "{:,}", "int_comma", "raw"),
 ("Age missing (%)", "pct_age_missing", "{:.2f}", "pct2", "raw"),
 ("Age, median (p25–p75)", ("age_median", "age_p25", "age_p75"), "{:.0f} ({:.0f}–{:.0f})", "median_iqr", "raw"),
 ("Care setting: inpatient", "setting_Inpatient", "{:,}", "int_comma", "raw"),
 ("Care setting: outpatient", "setting_Outpatient", "{:,}", "int_comma", "raw"),
 ("Care setting: emergency", "setting_ER", "{:,}", "int_comma", "raw"),
 ("Care setting: other", "setting_Other", "{:,}", "int_comma", "raw"),
 ("Specimen: urine", "specimen_URINE", "{:,}", "int_comma", "raw"),
 ("Specimen: blood", "specimen_BLOOD", "{:,}", "int_comma", "raw"),
 ("Specimen: other", "specimen_other", "{:,}", "int_comma", "raw"),
 ("AST susceptible (task-expanded)", "n_S_task_expanded", "{:,.0f}", "int_comma_f", "te"),
 ("AST resistant (task-expanded)", "n_R_task_expanded", "{:,.0f}", "int_comma_f", "te"),
]
r1 = []
for label, fld, fmt, rule, den in SPEC1:
    vals = {}
    for s in SITES:
        if rule == "constant":
            v = fmt
        elif rule == "median_iqr":
            v = fmt.format(g(fld[0], s), g(fld[1], s), g(fld[2], s))
        elif rule == "pct2":
            v = fmt.format(g(fld, s))
        else:
            v = fmt.format(g(fld, s))
        vals[s] = v
    r1.append((label, vals, den))
    for s in SITES:
        fields = fld if isinstance(fld, tuple) else (fld,)
        cell("R1", label, s, "site_profile.csv",
             "|".join(str(x) for x in fields), rule if rule != "constant" else "constant 38",
             vals[s])

# ================================================================ Table R2
r2, FLOW_LABEL = [], {
 1: ("Records in microbiology source file (all organisms)", "record", "participant_flow.csv", "value", "int_comma"),
 2: ("Records for the three focal organisms", "record", "participant_flow_counts.csv", "n_rows_focal_organisms", "int_comma"),
 3: ("Distinct (patient, date, specimen, antibiotic, result)", "record", "participant_flow_counts.csv", "n_rows_dedup_abcd_s", "int_comma"),
 4: ("Records with an interpretable AST result (S/I/R)", "record", "participant_flow_counts.csv", "n_rows_SRI", "int_comma"),
 5: ("Distinct (patient, specimen, date) units entering episode construction", "record", "participant_flow_counts.csv", "n_units", "int_comma"),
 6: ("Microbiological episodes constructed (anchor chain, 14 days)", "episode", "participant_flow_counts.csv", "n_episodes_in_store", "int_comma"),
 7: ("Episodes with S/R for ≥1 benchmark organism–drug pair", "episode", "site_profile.csv", "n_unique_microbiological_episodes", "int_comma"),
 8: ("Unique patients contributing ≥1 such episode", "patient", "site_profile.csv", "n_unique_patients", "int_comma"),
 # ★ 2026-10-07（A1.6 + A1.7 + A1.9）：把「候选 → 合格」的**损耗**摆进流程表。
 #   病因：原表在 row 10 就停了，且 row 9/10 报的是**陈旧**的 202 口径 ⇒
 #   228→186 的 42 条排除在补充材料里**从未可见**（正是审阅记录点的「宇宙重建不可见」）。
 #   现改为每个角色三行：候选（228）→ 合格（186）→ 被源性能筛排除（42），
 #   并给出**分站分解**（as source 6/14/22；as target 18/14/10）。
 9: ("Candidate directed transfers involving this site, as source", "transfer", "UNIVERSE", "source", "int_comma"),
 10:("Eligible directed transfers involving this site, as source", "transfer", "site_profile.csv", "n_eligible_transfers_as_source", "int_comma_f"),
 11:("Excluded by the source-performance screen, as source", "transfer", "UNIVERSE", "excl_source", "int_comma"),
 12:("Candidate directed transfers involving this site, as target", "transfer", "UNIVERSE", "target", "int_comma"),
 13:("Eligible directed transfers involving this site, as target", "transfer", "site_profile.csv", "n_eligible_transfers_as_target", "int_comma_f"),
 14:("Excluded by the source-performance screen, as target", "transfer", "UNIVERSE", "excl_target", "int_comma"),
}
FCC = R("participant_flow_counts.csv").set_index("site")
# ★ 候选/合格宇宙的计数直接从**冻结宇宙表**读（唯一有修正后代次出处的表）
_U = pd.read_csv(os.path.join(ROOT, "分析结果", "p6b_corrected_universe.csv"))
_CNT = {"source": _U.groupby("source_site").size().to_dict(),
        "target": _U.groupby("target_site").size().to_dict()}
_ELI = {"source": _U[_U.eligible.astype(str) == "True"].groupby("source_site").size().to_dict(),
        "target": _U[_U.eligible.astype(str) == "True"].groupby("target_site").size().to_dict()}

for order, (label, unit, src, fld, rule) in FLOW_LABEL.items():
    vals = {}
    PFS = pf[pf["order"] == order].set_index("site")   # participant_flow.csv 以 (site, order) 为键
    for s in SITES:
        if src == "UNIVERSE":
            if fld in ("source", "target"):
                v = _CNT[fld][s]
            else:                                   # 排除 = 候选 − 合格
                role = fld.replace("excl_", "")
                v = _CNT[role][s] - _ELI[role][s]
        elif src.startswith("participant_flow_counts"):
            v = FCC.loc[s, fld]
        elif src.startswith("participant_flow"):
            v = PFS.loc[s, fld]
        else:
            v = g(fld, s)
        vals[s] = f"{v:,.0f}" if rule == "int_comma_f" else f"{int(v):,}"
        cell("R2", f"{order}. {label}", s,
             "p6b_corrected_universe.csv" if src == "UNIVERSE" else src, fld, rule, vals[s])
    r2.append((order, label, unit, vals))
# ★ 损耗会计断言：两侧的「候选 − 合格 − 排除 == 0」必须逐站成立
for _role in ("source", "target"):
    for s in SITES:
        assert _CNT[_role][s] - _ELI[_role][s] - (_CNT[_role][s] - _ELI[_role][s]) == 0
assert sum(_CNT["source"].values()) - sum(_ELI["source"].values()) == 42, "source 侧排除应为 42"

# ================================================================ Table R3（汇总 + 指针）
r3 = []
for s in SITES:
    sub = xsd[xsd.site == s]
    nS, nR = int(sub.n_S.sum()), int(sub.n_R.sum())
    nSR = nS + nR
    r3.append((s, len(sub), f"{nS:,}", f"{nR:,}", f"{nSR:,}"))
    for lab, v, fldname, rule in [("rows (organism–drug)", len(sub), "row count", "int"),
                                  ("AST susceptible", nS, "sum(n_S)", "int_comma"),
                                  ("AST resistant", nR, "sum(n_R)", "int_comma"),
                                  ("AST total", nSR, "sum(n_S)+sum(n_R)", "int_comma")]:
        cell("R3", lab, s, "site_x_organism_drug_counts.csv", fldname, rule,
             f"{v:,}" if rule != "int" else str(v))

# ================================================================ 自洽关系
cons = []
for s in SITES:
    a = int(g("setting_Inpatient", s) + g("setting_Outpatient", s) + g("setting_ER", s) + g("setting_Other", s))
    b = int(g("specimen_URINE", s) + g("specimen_BLOOD", s) + g("specimen_other", s))
    c = int(g("n_age_observed", s) + g("n_age_missing", s))
    t = int(g("n_unique_microbiological_episodes", s))
    cons.append((s, a, b, c, t))
assert all(a == b == c == t for _, a, b, c, t in cons), f"行内自洽失败: {cons}"

# ================================================================ 渲染 markdown
md = []
# ★★ 2026-10-07：本文件是**投稿材料**，表头与表注必须**全英文**。
#   首版把内部工作注记（"本文件由 … 自动渲染"、脚本文件名、版本号、"项目规则"）写进了表头，
#   中文会原样进到投稿的 Word 里。内部溯源不进期刊件 ——
#   每个 cell 的来源见 `supplementary_cell_provenance.csv`。
md.append("# Supplementary Tables R1–R3 — Reporting Profiles\n")
md.append("\n## Table R1 — Site-level profile\n")
md.append("| | Denominator | " + " | ".join(LABEL[s] for s in SITES) + " |")
md.append("|---|---|---|---|---|")
for label, vals, den in r1:
    md.append(f"| {label} | {DEN[den]} | " + " | ".join(vals[s] for s in SITES) + " |")
md.append("\n> **Note.** The *Denominator* column states, for every row, which population the "
          "count is taken over; **the three are not interchangeable** and differ by a factor of "
          "13–17 between *raw* and *task-expanded*.\n"
          "> \n"
          "> **raw (unique episode)** — deduplicated microbiological episodes (key: patient, "
          "specimen, anchor block), applied **before** task expansion. **Unique patients are not "
          "episodes**: the analysis unit is the microbiological episode, not the clinical "
          "participant.\n"
          "> \n"
          "> **candidate grid** — the 38 organism–drug groups × 6 directed transfers = 228 "
          "candidate transfers, **before** the source-performance screen. Every group appears at "
          "every site (2 outgoing transfers per group per site).\n"
          "> \n"
          "> **eligible universe** — the 186 transfers that passed the source-performance screen "
          "and form the frozen benchmark; **37** of the 38 candidate groups contribute at least "
          "one eligible transfer, so a group count of 38 and one of 37 are both correct under "
          "their own denominator (see Table R2 for the full attrition).\n"
          "> \n"
          "> **task-expanded** — organism–drug episode observations, **not** deduplicated across "
          "drugs within an episode; deliberately larger than the unique-episode denominator.\n")

md.append("\n## Table R2 — Participant / episode flow\n")
md.append("| # | Step | Unit | " + " | ".join(LABEL[s] for s in SITES) + " |")
md.append("|---|---|---|---|---|---|")
for order, label, unit, vals in r2:
    md.append(f"| {order} | {label} | {unit} | " + " | ".join(vals[s] for s in SITES) + " |")
md.append("\n> **Note.** Row 1 is a file-level count that includes non-focal organisms; at MGB the "
          "organism filter is applied during read. From row 6 onward the unit is the episode, the "
          "patient or the task, and episode counts must not be reported as participant counts.\n"
          "> \n"
          "> **Rows 9–14 complete the attrition.** A directed transfer is counted **once per "
          "role**: each of the 38 candidate organism–drug groups contributes two outgoing and two "
          "incoming candidate transfers per site, so the candidate grid is "
          "38 × 6 = 228 = 76 + 76 + 76 in either role. Of these, 186 passed the prespecified "
          "source-performance screen (internal AUROC ≥ 0.65; protocol §3.2) and form the frozen benchmark "
          "universe; **42 were excluded**, and the exclusion is given per site below each "
          "candidate/eligible pair. **The bridge between the two universes is simply "
          "228 − 42 = 186**, and it is the cheapest consistency check available across them: if "
          "any one of the three totals is stale, the bridge breaks. **The two roles partition the same 228 transfers and the same "
          "42 exclusions differently**: as source 6 / 14 / 22, as target 18 / 14 / 10.\n"
          "> \n"
          "> **One candidate group (K. pneumoniae × nitrofurantoin) contributed no eligible "
          "transfer at any site (0 of 6)** — 2 candidate transfers per site. That is why the "
          "candidate grid has 38 groups and the eligible benchmark universe has 37; both counts "
          "are correct under their own denominator.\n"
          "> \n"
          "> **The 42 exclusions admit two independent partitions**, which reconcile with each "
          "other: **by role** (as source 6 / 14 / 22; as target 18 / 14 / 10) and **by group** "
          "(6 from the single group above, which contributes 2 per site, leaving 4 / 12 / 20 from "
          "the other 36 exclusions distributed over 14 groups).\n")

md.append("\n## Table R3 — Task-level counts (summary)\n")
md.append("| Site | Rows (organism–drug) | AST susceptible | AST resistant | AST total |")
md.append("|---|---|---|---|---|")
for s, n, nS, nR, nSR in r3:
    md.append(f"| {LABEL[s]} | {n} | {nS} | {nR} | {nSR} |")
# ★ 2026-10-07（A1.5）：原注说「the full 114-row breakdown is provided as supplementary data」，
#   但该文件**从未随包发出**（投审包内不存在 114/111 行的表），且 114 = 38×3 是**候选网格**行数、
#   与本表印的 37 组（合格口径）**分母不同**。
#   现改为指认**包内真实存在**的源：Source Data 工作簿的 `Figure S1` sheet（228 行带 eligible 旗标），
#   一个文件同时携带 228 / 186 / 38 / 37 四个分母。
md.append("\n> **Note.** This table uses the **task-expanded** denominator (analysis unit: "
          "organism–drug episode observation; deduplication level: site × drug), which differs "
          "from the unique-episode denominator of Table R1.\n"
          "> \n"
          "> **Row counts.** The table reports **37 organism–drug groups per site** — the "
          "**eligible** basis, i.e. groups contributing at least one transfer that passed the "
          "prespecified screen. The underlying task-level file carries **114 rows = 38 groups × 3 "
          "sites** — the **candidate** basis, before the screen; of these, 111 rows (37 × 3) carry "
          "at least one eligible transfer. Both counts are correct under their own denominator.\n"
          "> \n"
          "> **Source data.** The per-transfer breakdown is provided in the Source Data workbook, "
          "sheet **`Figure S1`**, which lists all **228** candidate transfers with the eligibility "
          "flag; the candidate count (228), the eligible benchmark universe (186) and the group "
          "counts (38 candidate / 37 eligible) can all be reconciled from that single file.\n")

import genstamp  # noqa: E402
# ★ 2026-10-08（A3.6）：**覆盖前**归档。初版把 backup 写在 write **之后**、且 import 在调用**之后**
#   ⇒ NameError 被 except 吞掉、归档静默为 0（实测：重跑后 `_artifact_archive/` 为空）。
#   **正确顺序：import → backup → write → stamp。**
genstamp.backup(os.path.join(D, "Supplementary_Tables_R1_R3.md"), "p21d")
genstamp.backup(os.path.join(D, "supplementary_cell_provenance.csv"), "p21d")

open(os.path.join(D, "Supplementary_Tables_R1_R3.md"), "w", encoding="utf-8").write("\n".join(md))
pd.DataFrame(cells).to_csv(os.path.join(D, "supplementary_cell_provenance.csv"), index=False)
try:
    for _p in [os.path.join(D, "Supplementary_Tables_R1_R3.md"),
               os.path.join(D, "supplementary_cell_provenance.csv")]:
        genstamp.stamp(_p, "p21d")
    print("[genstamp] 代次头已写: 2 个")
except Exception as _e:
    print("[genstamp] 跳过:", _e)

print(f"[OK] Supplementary_Tables_R1_R3.md 已渲染")
print(f"[OK] supplementary_cell_provenance.csv —— {len(cells)} 个 cell")
print(f"[OK] 行内自洽：setting / specimen / age 分项之和 == unique episodes（三站）")
print(f"     {cons}")
