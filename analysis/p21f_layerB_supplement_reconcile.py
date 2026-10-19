# -*- coding: utf-8 -*-
r"""p21f — Layer B 扩展：Supplementary Tables R1–R3 的**逐 cell 溯源**

================================ 为什么单独一层 ================================
2026-10-02 事故：手抄 Table R1 时**编造了年龄四分位区间**。
`p21c` 的 manifest 只保证**产物文件未被改动**（sha256），
**不保证 Supplement 里呈现的数字与产物一致** —— 所以这类错误一个审计器都抓不到。

本层按 GPT6 指定的链路逐 cell 核查：
    Supplement cell → reporting artifact → field(s) → transformation rule
    → expected display value → observed display value → PASS/FAIL

**关键设计**：expected 由本脚本**独立重派生**（不复用 `p21d` 的渲染代码路径），
与 `p21d` 写下的 displayed_value 比对 —— 两条独立路径互验。

覆盖面（GPT6 指定）：原始计数 · 百分比 · 汇总量(median) · 分位数(p25/p75) ·
拼接显示值(如 `65 (45–75)`) · 比率 · 行内分项自洽 · footnote 中的数字。

用法：python p21f_layerB_supplement_reconcile.py
输出：report_p21f_supplement_cells.txt / p21_layerB_supplement_cells.csv
"""
import os, re, io, sys
import numpy as np, pandas as pd

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); OUT = os.path.dirname(os.path.abspath(__file__))
RP = os.path.join(ROOT, "figures", "data", "reporting_profiles")

rep = open(os.path.join(OUT, "report_p21f_supplement_cells.txt"), "w", encoding="utf-8")
def p(*a):
    line = " ".join(str(x) for x in a); print(line, flush=True); rep.write(line + "\n"); rep.flush()

CP = os.path.join(RP, "supplementary_cell_provenance.csv")
MD = os.path.join(RP, "Supplementary_Tables_R1_R3.md")
if not os.path.exists(CP):
    p("*** 尚无 supplementary_cell_provenance.csv —— 请先运行 p21d ***"); rep.close(); sys.exit(2)

cp    = pd.read_csv(CP)
md    = io.open(MD, encoding="utf-8").read()
SP    = pd.read_csv(os.path.join(RP, "site_profile.csv")).set_index("site")
FC    = pd.read_csv(os.path.join(RP, "participant_flow_counts.csv")).set_index("site")
PF    = pd.read_csv(os.path.join(RP, "participant_flow.csv"))
XS    = pd.read_csv(os.path.join(RP, "site_x_organism_drug_counts.csv"))

p("=" * 104)
p("p21f — Supplementary Tables R1–R3 逐 cell 溯源（渲染值 vs **独立重派生**）")
p("=" * 104)
p(f"  cell 总数：{len(cp)}")
p("")

# ★ 2026-10-07：R2 新增的 9–14 行（候选/合格/排除）来自**冻结宇宙表**，不是 participant_flow。
#   p21f 必须能**独立**重派生它们，否则验证器覆盖跟不上渲染器 ——
#   首次运行时确实崩在这里（IndexError），那是门在响，不是噪音。
UNI = pd.read_csv(os.path.join(ROOT, "分析结果", "p6b_corrected_universe.csv"))


def _from_universe(site, fields):
    """从冻结宇宙表独立重派生（不读 p21d 的任何中间量）。"""
    if fields in ("source", "target"):
        return int((UNI["source_site" if fields == "source" else "target_site"] == site).sum())
    role = fields.replace("excl_", "")
    col = "source_site" if role == "source" else "target_site"
    n_all = int((UNI[col] == site).sum())
    n_ok = int(((UNI[col] == site) & (UNI.eligible.astype(str) == "True")).sum())
    return n_all - n_ok


def rederive(tbl, label, site, fields, rule, srcfile=""):
    """**独立**按变换规则重派生显示值（刻意不复用 p21d 的代码路径）。"""
    f = [x for x in str(fields).split("|") if x]
    if rule == "constant 38":
        return "38"
    if rule == "median_iqr":
        a, b, c = (SP.loc[site, x] for x in f)
        return f"{a:.0f} ({b:.0f}\u2013{c:.0f})"
    # 注意：p21d 写 provenance 时把求和式放在 artifact_fields 列（transform_rule 记的是显示格式），
    # 故此处按 **fields** 识别，而不是按 rule。
    if str(fields) in ("sum(n_S)", "sum(n_R)", "sum(n_S)+sum(n_R)"):
        sub = XS[XS.site == site]
        v = {"sum(n_S)": sub.n_S.sum(), "sum(n_R)": sub.n_R.sum(),
             "sum(n_S)+sum(n_R)": sub.n_S.sum() + sub.n_R.sum()}[str(fields)]
        return f"{int(v):,}"
    if str(fields) == "row count":
        return str(len(XS[XS.site == site]))
    # ★ 宇宙表来源（R2 的新 9–14 行）
    if "p6b_corrected_universe" in str(srcfile):
        return f"{_from_universe(site, str(fields)):,}"
    # 标量来源
    # ★ 2026-10-07：R2 的行现在有**三种**来源，必须按 artifact_file 分派，
    #   不能像初版那样一律回落到 participant_flow 的 order 查找
    #   （行 10/13 的合格转移来自 site_profile，order 13 在 participant_flow 里不存在 ⇒ 崩溃）。
    if tbl == "R2":
        if "site_profile" in str(srcfile):
            v = SP.loc[site, str(fields)]
        elif fields in FC.columns:
            v = FC.loc[site, fields]
        else:
            order = int(re.match(r"(\d+)\.", label).group(1))
            _hit = PF[(PF.site == site) & (PF["order"] == order)]
            if len(_hit) == 0:
                raise KeyError(f"participant_flow 中无 order={order}（行：{label}）")
            v = _hit.value.iat[0]
    else:
        v = SP.loc[site, f[0]]
    if rule == "pct2":       return f"{float(v):.2f}"
    if rule == "int_comma":  return f"{int(v):,}"
    if rule == "int_comma_f":return f"{float(v):,.0f}"
    raise KeyError(f"未定义的变换规则: {rule}")

rows = []
for _, r in cp.iterrows():
    exp = rederive(r["table"], r["row_label"], r["site"], r["artifact_fields"], r["transform_rule"], r.get("artifact_file", ""))
    obs = str(r["displayed_value"])
    in_md = obs in md
    rows.append(dict(table=r["table"], row=r["row_label"], site=r["site"],
                     artifact_fields=r["artifact_fields"], transform_rule=r["transform_rule"],
                     expected_display=exp, observed_display=obs, present_in_rendered_md=in_md,
                     verdict="PASS" if (exp == obs and in_md) else "FAIL"))
D = pd.DataFrame(rows)
D.to_csv(os.path.join(OUT, "p21_layerB_supplement_cells.csv"), index=False)

p(f"  【逐 cell】PASS {int((D.verdict=='PASS').sum())}｜FAIL {int((D.verdict=='FAIL').sum())}")
for _, r in D[D.verdict == "FAIL"].iterrows():
    p(f"    FAIL {r['table']}/{r['row']}/{r['site']}  规则={r['transform_rule']} "
      f"期望={r['expected_display']!r} 渲染={r['observed_display']!r} 出现在 MD={r['present_in_rendered_md']}")

# ---------------------------------------------------------------- 按 cell 类型分组统计
p("")
p("  【按呈现类型】")
GRP = {
 "原始计数":      ("int_comma", "int", "row count"),
 "拼接显示值":    ("median_iqr",),
 "百分比":        ("pct2",),
 "汇总量/分位数": ("int_comma_f", "constant 38"),
 "汇总量(sum)":   ("int_comma",),
}
for nm, rules in GRP.items():
    sub = D[D.transform_rule.isin(rules)]
    p(f"    {nm:<16}{len(sub):>4} cell　PASS {int((sub.verdict=='PASS').sum())}")

# ---------------------------------------------------------------- 行内自洽（footnote 依据）
p("")
p("  【行内自洽 —— 表注中'分项之和等于 unique episodes'的机器依据】")
ok_all = True
for s in ["MGB", "SU", "UT"]:
    a = int(SP.loc[s, "setting_Inpatient"] + SP.loc[s, "setting_Outpatient"]
            + SP.loc[s, "setting_ER"] + SP.loc[s, "setting_Other"])
    b = int(SP.loc[s, "specimen_URINE"] + SP.loc[s, "specimen_BLOOD"] + SP.loc[s, "specimen_other"])
    c = int(SP.loc[s, "n_age_observed"] + SP.loc[s, "n_age_missing"])
    t = int(SP.loc[s, "n_unique_microbiological_episodes"])
    good = (a == b == c == t)
    ok_all &= good
    p(f"    {s:<4} setting {a:>7,} ｜ specimen {b:>7,} ｜ age {c:>7,}  == unique episodes {t:>7,}  "
      f"{'OK' if good else '**FAIL**'}")

# ---------------------------------------------------------------- 比率类（R-2）
p("")
p("  【比率检查 R-2：两组分母不得混用】")
for s in ["MGB", "SU", "UT"]:
    u = float(SP.loc[s, "n_unique_microbiological_episodes"])
    t = float(SP.loc[s, "n_task_expanded_episode_observations"])
    p(f"    {s:<4} unique={int(u):>7,}  task-expanded={int(t):>10,}  比值={t/u:>5.1f}x"
      f"   → 渲染表中两行**并列且各自标注口径**")

# ---------------------------------------------------------------- 跨路径一致性（★ 2026-10-07 新增）
# 【为什么加这一层】
#   病根：R1 与 R3 报的是**同一个量**（task-expanded episode observations），
#   却走**两条不同路径**（R1 ← task_profile.n_episodes_sr 求和；R3 ← site_x 的 n_S+n_R 求和），
#   中间**没有任何相等断言**。
#   2026-10-07 实测：修复前两者相差 76,408 / 6,102 / 7,851（合计 90,361），
#   而本脚本当时对分歧态**报 PASS 99/FAIL 0** —— 因为它只问「cell 抄自其来源吗」，
#   从不问「同值兄弟相等吗」。（物证：RECONSTRUCTED_PREFIX_VERDICT_20261007/）
#   ⇒ 出处层存在而一致性层缺失。本块把一致性层接进**生产路径**，失败即计入 nfail。
#   规则依据：WORK_CHECKLIST D 节 11/12/13；「岗哨有能力而没站在门上，等于没有」。
p("")
p("  【跨路径一致性 —— R1.task-expanded ≡ R3.AST_total ≡ 渲染件】")
XCH = []
for s in ["MGB", "SU", "UT"]:
    _r1 = float(SP.loc[s, "n_task_expanded_episode_observations"])
    _sub = XS[XS.site == s]
    _r3 = float(_sub.n_S.sum() + _sub.n_R.sum())
    _inmd = f"{int(round(_r1)):,}" in md
    _ok = (abs(_r1 - _r3) < 0.5) and _inmd
    XCH.append(_ok)
    _d = "" if abs(_r1 - _r3) < 0.5 else f"  ← 相差 {abs(_r1-_r3):,.0f}"
    p(f"    {s:<4} R1={int(round(_r1)):>10,}  R3={int(round(_r3)):>10,}"
      f"  {'≡' if abs(_r1-_r3)<0.5 else '≠'}{_d}   渲染件含 R1 值={_inmd}")
p(f"    → {'两条路径一致' if all(XCH) else '**跨路径分歧 —— 判红**'}")

# ---------------------------------------------------------------- wrapper 泄漏检查
WR = os.path.join(ROOT, "Supplementary_Reporting_Profiles_v1.md")
leak = []
if os.path.exists(WR):
    w = io.open(WR, encoding="utf-8").read()
    derived = set()
    for col in ("n_unique_microbiological_episodes", "n_unique_patients",
                "n_task_expanded_episode_observations", "n_age_missing",
                "n_S_task_expanded", "n_R_task_expanded"):
        for v in SP[col]:
            derived.add(f"{int(v):,}"); derived.add(str(int(v)))
    for tok in re.findall(r"(?<![\w.])\d[\d,]*(?:\.\d+)?", w):
        t = tok.replace(",", "")
        if tok in derived or t in derived:
            leak.append(tok)
p("")
p(f"  【wrapper 泄漏检查】Supplementary_Reporting_Profiles_v1.md 中出现的**派生数据值**："
  f"{sorted(set(leak)) if leak else '无'}")
p("     （政策：派生数据值只能出现在渲染表中；wrapper 仅允许结构性叙述计数与条目编号）")

# ---------------------------------------------------------------- 负向测试
# **必须证明本层能抓住 2026-10-02 那个 age IQR 事故**，否则它只是装饰。
p("")
p("=" * 104)
p("负向测试 —— 人为篡改 displayed_value，验证逐 cell 溯源会 FAIL")
p("=" * 104)
MUT = [
 ("age 四分位被改写（复现 2026-10-02 事故）", "R1", "Age, median (p25–p75)", "MGB", "51 (45–77)"),
 ("unique episodes 被写成患者数",            "R1", "Unique microbiological episodes", "MGB", "67,079"),
 ("task-expanded 被写成 episode 数",         "R1", "Task-expanded episode observations", "MGB", "89,454"),
 ("百分比精度被改",                          "R1", "Age missing (%)", "MGB", "6.6"),
]
neg_fail = 0
for name, tbl, row, site, bad in MUT:
    m = cp[(cp.table == tbl) & (cp.row_label == row) & (cp.site == site)]
    if not len(m):
        p(f"  [跳过] {name}：未找到该 cell"); neg_fail += 1; continue
    exp = rederive(tbl, row, site, m.artifact_fields.iat[0], m.transform_rule.iat[0])
    caught = (exp != bad)
    p(f"  [{'CAUGHT' if caught else 'MISSED'}] {name:<40} 真值={exp!r} 篡改={bad!r}")
    if not caught: neg_fail += 1
p(f"\n  负向测试：{len(MUT)} 项，未捕获 {neg_fail} 项")

nfail = (int((D.verdict == "FAIL").sum()) + (0 if ok_all else 1)
         + (0 if all(XCH) else 1)          # ★ 跨路径一致性计入总分
         + len(set(leak)) + neg_fail)
p("")
# ★ 判词自带**射程**（WORK_CHECKLIST 规则 13：绿灯比红灯危险 —— 人读到 PASS 默认「一切正常」）。
#   generation 一栏要等 A3.3（代次头）落地后才能填成 checked，现在**必须显式声明 unchecked**。
_SCOPE = "[scope: copy-fidelity + cross-path consistency | generation: UNCHECKED]"
p(f"*** p21f 通过 ***   {_SCOPE}" if nfail == 0 else f"*** p21f 未通过 —— {nfail} 项 ***")
if nfail == 0:
    p("    注：本判词**不覆盖代系** —— 两份产物可能同属一个陈旧代次而彼此一致。")
rep.close(); print("\ndone ->", OUT)
sys.exit(0 if nfail == 0 else 1)
