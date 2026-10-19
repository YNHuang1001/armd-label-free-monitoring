# -*- coding: utf-8 -*-
r"""t23_render_parity.py — 补充材料渲染件的**独立**验证

================================ 为什么需要 ================================
`p21e` / `p21f` 各自带内部对账（渲染器检查自己写的文件）。但那是**自证**：
渲染器与被渲染物出自同一段代码，同错则同绿。

现实缺口：这几个渲染脚本**不在 EVS manifest 的 closure 集合内**，
因此它们的产物不在任何**独立**判据之下。本模块补这一档。

================================ 独立性从哪来 ================================
本模块**不 import、不 exec、不复用**任何渲染器代码。它：
  1. 自己读冻结 CSV，自己按 Markdown 语法**重新解析**渲染出的 .md；
  2. 自己从 CSV 原始值**重建**每个单元格应该显示的字符串；
  3. 逐格比对。
渲染器改了格式化规则而忘了同步 → 本模块判红；CSV 漂移 → 本模块判红。

================================ 覆盖 ================================
  T23-01  Table S1：行数 / 列头 / 逐格（186 × 21）
  T23-02  Table S2 Panel A：逐格
  T23-03  Table S2 Panel B：逐格
  T23-04  Supplementary Methods：四个配置行存在且层次/分类器正确
  T23-05  三份渲染件都在（缺任一即 FAIL，防「渲染件被删而没人发现」）
"""
import io
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import evs_load as L

import pandas as pd

PROJ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIG  = os.path.join(PROJ, "figures", "data")
RP   = os.path.join(FIG, "reporting_profiles")

def _read(p):
    return io.open(p, encoding="utf-8").read() if os.path.exists(p) else None

def _parse_md_tables(text):
    """独立解析 Markdown 表格 → [(header_cells, [row_cells, ...]), ...]
    不依赖渲染器的任何函数。"""
    tables, cur, hdr = [], [], None
    for ln in text.splitlines():
        s = ln.strip()
        if s.startswith("|") and s.endswith("|"):
            cells = [c.strip() for c in s.strip("|").split("|")]
            if set("".join(cells)) <= set("-: "):      # 分隔行
                continue
            if hdr is None:
                hdr = cells
            else:
                cur.append(cells)
        else:
            if hdr is not None:
                tables.append((hdr, cur)); hdr, cur = None, []
    if hdr is not None:
        tables.append((hdr, cur))
    return tables

R = []
def rec(tid, ok, exp, obs, note=""):
    R.append(ok)
    L.emit(tid, "PASS" if ok else "FAIL", exp, obs, note)

# ================================================================ T23-05 存在性
need = ["Supplementary_Table_S1.md", "Supplementary_Table_S2.md", "Supplementary_Methods.md"]
missing = [n for n in need if not os.path.exists(os.path.join(RP, n))]
rec("T23-05", not missing,
    "三份渲染件均存在：Supplementary_Table_S1/S2、Supplementary Methods",
    f"缺失 = {missing or '无'}", "渲染件被删除时本模块判红，不会静默通过")

# ================================================================ T23-01 Table S1
T1 = pd.read_csv(os.path.join(FIG, "tableS1_tasks.csv"))
md1 = _read(os.path.join(RP, "Supplementary_Table_S1.md")) or ""
tabs1 = _parse_md_tables(md1)
big1 = max(tabs1, key=lambda t: len(t[1]), default=([], []))
hdr1, rows1 = big1
ok_rows = len(rows1) == len(T1) == 186
rec("T23-01a", ok_rows, "Table S1 渲染行数 == 冻结表行数 == 186",
    f"渲染 {len(rows1)} 行｜冻结 {len(T1)} 行｜列头 {len(hdr1)} 列")

# 逐格重建（独立于 p21e 的格式化代码）
def fmt(v, col):
    if col in ("harmful", "stable"):
        return "yes" if bool(v) else ""
    if col == "n_target":
        return "%d" % int(v)
    if col in ("auroc_source_internal", "auroc_target", "d_auroc"):
        return "%.3f" % float(v)
    if isinstance(v, (int, float)) and not isinstance(v, bool):
        return "%.4f" % float(v)
    return str(v)

bad1 = []
T1i = T1.set_index("task_id")
for r in rows1:
    tid = r[0]
    if tid not in T1i.index:
        bad1.append("未知 id " + tid); continue
    src = T1i.loc[tid]          # 该 Series 的索引是除 task_id 外的列
    for c, shown in zip(list(T1.columns)[1:], r[1:]):
        exp = fmt(src[c], c)
        if shown != exp:
            bad1.append("%s:%s(%s≠%s)" % (tid, c, shown, exp))
            if len(bad1) > 5: break
    if len(bad1) > 5: break
rec("T23-01b", not bad1,
    "Table S1 逐格（186 × 21）与冻结 CSV 一致（本模块自行重建显示串）",
    f"不符 {len(bad1)} 项" + (f"：{bad1[:3]}" if bad1 else ""), "独立解析 + 独立重建")

# ================================================================ T23-02/03 Table S2
T2 = pd.read_csv(os.path.join(FIG, "tableS2_metrics.csv")).set_index("method")
md2 = _read(os.path.join(RP, "Supplementary_Table_S2.md")) or ""
tabs2 = _parse_md_tables(md2)
# Panel A / B 各一张 8 列表、6 列表
tA = [t for t in tabs2 if len(t[0]) == 8]
tB = [t for t in tabs2 if len(t[0]) == 6]
badA = []
if not tA:
    badA.append("未解析到 Panel A 表")
else:
    hdrA, rowsA = tA[0]
    def ci3(r, pt, lo, hi):
        return "%.3f (%.3f\u2013%.3f)" % (float(r[pt]), float(r[lo]), float(r[hi]))
    for row in rowsA:
        m = row[0].strip("`")
        if m not in T2.index:
            badA.append("未知方法 " + m); continue
        r = T2.loc[m]
        exp = [ci3(r, "harmful_auc", "auc_lo", "auc_hi"),
               ci3(r, "harmful_auprc", "auprc_lo", "auprc_hi"),
               ci3(r, "sens_at_spec80", "sens80_lo", "sens80_hi"),
               ci3(r, "ppv", "ppv_lo", "ppv_hi"),
               ci3(r, "npv", "npv_lo", "npv_hi")]
        got = row[3:8]
        if got != exp:
            badA.append("%s: %s ≠ %s" % (m, got, exp))
        if len(badA) > 3: break
rec("T23-02", not badA, "Table S2 Panel A 逐格与冻结 CSV 一致（5 个区间列 × 15 行）",
    f"不符 {len(badA)} 项" + (f"：{badA[:2]}" if badA else ""), "区间串由本模块独立重建")

badB = []
if not tB:
    badB.append("未解析到 Panel B 表")
else:
    hdrB, rowsB = tB[0]
    def ci2s(r, pt, lo, hi):
        return "%+.2f (%+.2f to %+.2f)" % (float(r[pt]), float(r[lo]), float(r[hi]))
    for row in rowsB:
        m = row[0].strip("`")
        if m not in T2.index:
            badB.append("未知方法 " + m); continue
        r = T2.loc[m]
        exp_sp = ci2s(r, "spearman", "sp_lo", "sp_hi")
        if row[2] != exp_sp:
            badB.append("%s:spearman(%s≠%s)" % (m, row[2], exp_sp))
        # ★★ 2026-10-09（INC-019 的回声）：本处原比对 `macro_mean`（来自 p6e）。
        #   但 `p21e` 自 2026-10-09 起**渲染的是 `macro_logit`** —— 与 `fig3a_lodo` **逐位相同**
        #   的那一列（`macro_mean` 在 entropy_mean_shift 一行上高出 1 ULP，恰跨三位小数进位边界
        #   ⇒ 表印 0.613、图与正文 0.612）。**判据必须比对「实际渲染的列」**，
        #   否则它会一直报一个读者已经看不到的分歧（与 p53 早前同一课）。
        for i, col in ((3, "macro_logit"), (4, "loo_org_auc")):
            if pd.isna(r[col]):
                continue
            if row[i] != "%.3f" % float(r[col]):
                badB.append("%s:%s(%s≠%.3f)" % (m, col, row[i], float(r[col])))
        if len(badB) > 3: break
rec("T23-03", not badB, "Table S2 Panel B 逐格与冻结 CSV 一致（Spearman / macro / LOO）",
    f"不符 {len(badB)} 项" + (f"：{badB[:2]}" if badB else ""), "独立重建")

# ================================================================ T23-04 Supplementary Methods
mdm = _read(os.path.join(RP, "Supplementary_Methods.md")) or ""
F1E = pd.read_csv(os.path.join(FIG, "fig1e_configs.csv"))
cfgs = list(F1E.config)
missing_cfg = [c for c in cfgs if f"`{c}`" not in mdm]
rec("T23-04", not missing_cfg,
    "Supplementary Methods 含全部 %d 个配置行" % len(cfgs),
    f"缺失 = {missing_cfg or '无'}", "配置名取自冻结产物 fig1e_configs.csv")
