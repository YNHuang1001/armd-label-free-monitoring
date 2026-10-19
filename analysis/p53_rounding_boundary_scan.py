# -*- coding: utf-8 -*-
r"""p53 — **同一量、多生产者**：显示精度下的分歧扫描

================================ 为什么需要 ================================
INC-015（DR AUROC 区间两个种子）与 INC-016（p34 看不见小数）之外，2026-10-09 又撞上第三例：

    `entropy_mean_shift` 的 LODO macro 真值**恰好是 0.6125**，
    而两个产物各存了一个 1 ULP 的近似：
        `fig3a_lodo.csv`     -> 0.61249999999999993  （**低于** 0.6125）-> 印 0.612
        `tableS2_metrics.csv`-> 0.61250000000000004  （**高于** 0.6125）-> 印 0.613
    15 行里 9 行两文件相差 1 ULP，**只有这一行恰好跨在 3 位小数的进位边界上**，
    于是 1 ULP 的分歧变成了**读者可见的两个字符串**。

⇒ **ULP 级分歧通常不可见，因此通常无害；但可见性不是随机分布的 —— 真值越靠近进位边界，
同一个 ULP 越可能被放大成一位数字的差。**

================================ 本模块查什么 ================================
把**已知的「同一量」对应关系**显式登记（PAIRS），逐对比较：
  · `n_ulp`  两文件值不相等的格数（**允许** —— 不同求和顺序必然产生 ULP 噪声）
  · `n_3dp`  按**显示精度**渲染后字符串不同的格数（**报警** —— 读者会同时看到两者）

★ 为什么必须**显式登记**对应关系，而不是按列名自动配对：
  自动配对会把 `fig2_forest.auc_hi`（Fig 2A 的 harmful AUROC 区间上界）与
  `fig3a_lodo.hi_95`（Fig 3A 的 LODO macro 区间上界）配成一对 —— **列名族相同、量不同**，
  首版扫描因此报出 15 个假警报。**同名 ≠ 同量**，这是本仓已记录过的 114 同形异指的近亲。

================================ 判据 ================================
`n_3dp == 0`。任一格非零即退出码 1。
"""
import io
import itertools
import os
import sys

import numpy as np
import pandas as pd

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIG = os.path.join(ROOT, "figures", "data")
AN = os.path.join(ROOT, "分析结果")
REP = io.open(os.path.join(AN, "report_p53_rounding_scan.txt"), "w", encoding="utf-8")


def p(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    REP.write(s + "\n")
    REP.flush()


def rd(path, conv=None, key="method"):
    d = pd.read_csv(path)
    if conv:
        d = d[d.convention == conv]
    return d.set_index(key)


import rephead
# ★ 2026-10-09（外部核验）：输入声明必须**覆盖所有被读的文件**。
#   第 16 对（p32.mae vs fig4d.mae）读的这两份原先**不在声明里** ⇒
#   「3 0 0」那一行要么来自未声明输入（指纹不覆盖），要么是转抄 —— 两种都不该留。
#   现补入 `fig4d_mae.csv` 与 `p32_constant_baselines.csv`；**该行由脚本现算**（见 ②）。
for _l in rephead.header(__file__, [os.path.join(FIG, f) for f in
        ("fig2_forest.csv", "tableS2_metrics.csv", "fig3a_lodo.csv", "fig3b_loo_organism.csv",
         "fig4d_mae.csv")] + [os.path.join(AN, f) for f in
        ("p33_protocol72_metrics.csv", "p32_constant_baselines.csv", "p6e_lodo_fold_auc_long.csv")]):
    p(_l)
p("=" * 100)
p("p53 — 同一量多生产者：显示精度分歧扫描")
p("=" * 100)

# ---------------------------------------------------------------- 产物
F2 = rd(os.path.join(FIG, "fig2_forest.csv"))
TS2 = rd(os.path.join(FIG, "tableS2_metrics.csv"))
G3A = rd(os.path.join(FIG, "fig3a_lodo.csv"))
G3B = rd(os.path.join(FIG, "fig3b_loo_organism.csv"))
P33 = rd(os.path.join(AN, "p33_protocol72_metrics.csv"), conv="primary")
# ★ 2026-10-09（外部核验发现，第 16 对）：
#   **MAE 有两个生产者**：`p32_constant_baselines.csv`（→ 经 `p52` 喂 §4 的点估计与区间）
#   与 `fig4d_mae.csv`（→ 经 `p19` 喂 Discussion 的 0.0416/0.0417）。
#   **两条链在 Discussion 相遇。** 今日两者**逐位相同**（已实测 DR/IOW/OM 三个估计量），
#   故登记后不报警 —— 但登记的意义正在于此：**INC-019 事故前的形状就是「15 行里 14 行一致」**。
#   这是对 INC-019 的**补防**，不是对现有分歧的响应。
P32 = pd.read_csv(os.path.join(AN, "p32_constant_baselines.csv")).set_index("name")
F4D = pd.read_csv(os.path.join(FIG, "fig4d_mae.csv")).set_index("estimator")

# ---------------------------------------------------------------- ★ 显式登记的「同一量」对应
#   每一行 = (左标签, 右标签, 左表, 右表, 左列, 右列, 显示精度)
#   **只登记确知是同一个量的对**。不确定的不登记 —— 宁可漏，不可假。
PAIRS = [
    ("fig2_forest.harmful_auc",   "tableS2.harmful_auc",   F2, TS2, "harmful_auc", "harmful_auc", 3),
    ("fig2_forest.auc_lo",        "tableS2.auc_lo",        F2, TS2, "auc_lo", "auc_lo", 3),
    ("fig2_forest.auc_hi",        "tableS2.auc_hi",        F2, TS2, "auc_hi", "auc_hi", 3),
    ("fig2_forest.harmful_auprc", "tableS2.harmful_auprc", F2, TS2, "harmful_auprc", "harmful_auprc", 3),
    ("fig2_forest.sens_at_spec80", "tableS2.sens_at_spec80", F2, TS2, "sens_at_spec80", "sens_at_spec80", 3),
    ("fig2_forest.harmful_auc",   "p33.auroc",             F2, P33, "harmful_auc", "auroc", 3),
    ("fig2_forest.auc_lo",        "p33.auroc_lo",          F2, P33, "auc_lo", "auroc_lo", 3),
    ("fig2_forest.auc_hi",        "p33.auroc_hi",          F2, P33, "auc_hi", "auroc_hi", 3),
    # ↑ 这两对是 **INC-015 的已登记分歧**（两个 bootstrap 种子）。它们的点估计逐位相同、
    #   仅百分位界不同。**当前不存在读者可见的后果**：Fig 2 与 Table S2 的 AUROC 列同读
    #   `fig2_forest`；§7.2 的正文亦已于 2026-10-08 改绑 fig2。故 p33 的 AUROC 界
    #   **当前未被任何渲染件使用** ⇒ 归为 `latent`（潜伏），不报警，但**不删** ——
    #   删掉等于把「哪天有人渲染了 p33 的 AUROC」这个风险也一并删掉。
    ("fig2_forest.harmful_auprc", "p33.auprc",             F2, P33, "harmful_auprc", "auprc", 3),
    ("fig2_forest.sens_at_spec80", "p33.sens80",           F2, P33, "sens_at_spec80", "sens80", 3),
    # ★ 2026-10-09 更正（INC-019）：本对原写 `tableS2.macro_mean` —— 那是**渲染前**的候选列，
    #   而 Table S2 **实际渲染的是 `macro_logit`**（2026-10-09 起；此前渲染 `macro_mean`，
    #   在 entropy_mean_shift 一行与图/正文差一位）。**比对应以「实际渲染的列」为准**，
    #   否则扫描器会一直报一个读者已经看不到的分歧，而真正的渲染分歧无人看守。
    ("fig3a_lodo.lodo_macro_auroc", "tableS2.macro_logit", G3A, TS2, "lodo_macro_auroc", "macro_logit", 3),
    ("fig3a_lodo.lo_95",          "tableS2.C_lo",          G3A, TS2, "lo_95", "C_lo", 3),
    ("fig3a_lodo.hi_95",          "tableS2.C_hi",          G3A, TS2, "hi_95", "C_hi", 3),
    # ★ `fig3b_loo_organism.csv` **自带两列同量**：`loo_org_macro`（重算）与
    #   `p6d_loo_org_auc`（p6d 原值）。两列只差 1 ULP —— 即**同一个文件里也有两个生产者**。
    #   故两列分别与 Table S2 比对，看正文该跟随哪一个。
    # ★ 第 16 对（2026-10-09，外部核验发现）：**MAE 的两个生产者**。
    #   `p32` 的键含常数基线（`constant_zero/median`），故与 `fig4d` 的键交集恰为
    #   {DR, IOW, OM} 三行 —— **一对登记即覆盖三个估计量**。
    ("p32.mae（§4 链）",           "fig4d.mae（Discussion 链）", P32, F4D, "mae", "mae", 4),
    ("fig3b_loo.loo_org_macro",   "tableS2.loo_org_auc",   G3B, TS2, "loo_org_macro", "loo_org_auc", 3),
    ("fig3b_loo.p6d_loo_org_auc", "tableS2.loo_org_auc",   G3B, TS2, "p6d_loo_org_auc", "loo_org_auc", 3),
]

# ---------------------------------------------------------------- 非空自检
#   若因列名/键变动导致某一对**没有可比格**，必须响 —— 否则「0 分歧」与「没比」不可区分。
# ★ 2026-10-09：**不改门槛，改声明**。
#   新增的 MAE 对只有 3 个可比格（DR/IOW/OM），低于通用门槛 5。
#   把门槛降到 3 是「因为新案子过不去就放松判据」—— 本仓明令不做。
#   故改为**逐对声明期望格数**：声明了的按**精确值**判（比「≥5」更严），
#   未声明的仍守 ≥5 的下限。
EXPECT = {("p32.mae（§4 链）", "fig4d.mae（Discussion 链）"): 3}
probe = [(a, b, int(len(DA.index.intersection(DB.index)))) for a, b, DA, DB, _, _, _ in PAIRS]
empty = []
for _a, _b, _n in probe:
    _exp = EXPECT.get((_a, _b))
    if (_n != _exp) if _exp is not None else (_n < 5):
        empty.append((_a, _b, _n, _exp if _exp is not None else ">=5"))
assert not empty, "自检失败：以下对应关系的可比格数与声明不符或低于下限：%s" % empty
p("① 非空自检：%d 对对应关系 —— 已声明期望格数的按**精确值**判（%s），其余守 ≥5 下限  ✔"
  % (len(PAIRS), "、".join("%s=%d" % (k[0], v) for k, v in EXPECT.items()) or "无"))
p("")

# ---------------------------------------------------------------- 扫描
p("② 逐对比较（**ULP 差允许；显示精度差报警**）")
p(f"    {'同一量的两个生产者':<52}{'可比格':>7}{'ULP 差':>8}{'显示差':>8}")
LATENT = {("fig2_forest.auc_lo", "p33.auroc_lo"), ("fig2_forest.auc_hi", "p33.auroc_hi")}
ALARM = []
WATCH = []
for a, b, DA, DB, ca, cb, nd in PAIRS:
    idx = DA.index.intersection(DB.index)
    x = pd.to_numeric(DA.loc[idx, ca], errors="coerce")
    y = pd.to_numeric(DB.loc[idx, cb], errors="coerce")
    m = x.notna() & y.notna()
    x, y = x[m], y[m]
    nbit = int((x.values != y.values).sum())
    bad = [(i, float(u), float(v)) for i, u, v in zip(x.index, x.values, y.values)
           if ("%.*f" % (nd, u)) != ("%.*f" % (nd, v))]
    (WATCH if (a, b) in LATENT else ALARM).extend([(i, a, b, u, v, nd) for i, u, v in bad])
    tagx = "   **报警**" if bad and (a, b) not in LATENT else ("   (latent·已登记 INC-015)" if bad else "")
    p(f"    {a + '  vs  ' + b:<52}{len(x):>7}{nbit:>8}{len(bad):>8}{tagx}")
p("")

# ---------------------------------------------------------------- ★ 正典生产者 + 可判定的正典性
# GLM 建议（2026-10-09）：给每对登记「**渲染应跟随哪一侧**」，下次跨边界就不必由人裁决。
# 但**登记一个主张还不够** —— 主张本身也应尽量**机械化**。凡是能从上游精确算出来的，
# 就在这里**算出来**，而不是把「跟随图侧」当成一个约定写死了事。
# ★ 2026-10-09（外部核验要求）：每条正典主张必须标注**等级** ——
#   **可证**（有有理根基，能精确算出该跟随谁）｜**约定**（无「真值」，只能由生产者声明裁决）。
#   让将来读报告的人一眼看出：这条 canon 是**定理**还是**条约**。
#   等级**不写在这里** —— 由 `grade()` 算出来（见下）。本表只声明**跟随谁**及理由。
CANON = {
    ("fig2_forest.auc_lo", "p33.auroc_lo"):
        ("fig2", "INC-015 裁定 fig2 为主生产者；两者是两个 bootstrap 种子"),
    ("fig2_forest.auc_hi", "p33.auroc_hi"):
        ("fig2", "同上"),
    ("fig3a_lodo.lodo_macro_auroc", "tableS2.macro_logit"):
        ("图（同源）", "fig3a 与 macro_logit **逐位相同 15/15**；分歧发生在它与 macro_mean 之间"),
    ("fig3b_loo.loo_org_macro", "tableS2.loo_org_auc"):
        ("（潜伏）", "Table S2 现读 `p6d_loo_org_auc` 侧（与 S2 0 ULP）；`loo_org_macro` 3 行差 1 ULP"),
}
# ★ 2026-10-09（外部核验要求，但**做法改了**）：
#   外部建议「给每条 canon 标注 可证/约定」。**标签不是判据** —— 故本模块不写标签，
#   而是**尝试证明**：有冻结输入且能一次舍入的，就用 `Fraction` 算出精确值，
#   再看哪一侧的 double 等于 `float(exact)`。
#
#   两条轴**正交**，故分两列：
#     **等级** = 可证（有有理根基，算过）｜ 约定（无真值，只能声明）
#     **状态** = 活跃（两个消费者都在渲染 → 分歧读者可见）｜ 潜伏（仅一个在渲染）
#   实测提醒：**潜伏的对也可以是可证的**，反之亦然。
PROOF = {
    ("fig3a_lodo.lodo_macro_auroc", "tableS2.macro_logit"): dict(
        src=os.path.join(AN, "p6e_lodo_fold_auc_long.csv"), key="method", col="fold_auc",
        note="LODO macro = 可评价折 AUROC 的**未加权均值**"),
    ("fig3b_loo.loo_org_macro", "tableS2.loo_org_auc"): dict(
        src=os.path.join(FIG, "fig3b_loo_organism.csv"), key="method",
        cols=["loo_auc_EC", "loo_auc_KP", "loo_auc_PA"],
        note="LOO-organism macro = 3 个物种留一 AUROC 的均值"),
}


def grade(pair, sideA_vals, sideB_vals, nd):
    """返回 (等级, 说明)。能证明的返回「可证」，否则「约定」并说明为什么证不出。"""
    pr = PROOF.get(pair)
    if pr is None:
        return "约定", "无冻结输入可做精确算术（分位/重采样量无「真值」）"
    try:
        d = pd.read_csv(pr["src"])
        from fractions import Fraction
        diff = [m for m in sideA_vals if ("%.*f" % (nd, sideA_vals[m])) != ("%.*f" % (nd, sideB_vals.get(m, sideA_vals[m])))]
        if not diff:
            return "可证", "两侧已逐位相同，无分歧可裁"
        msgs = []
        for m in diff:
            r = d[d[pr["key"]] == m]
            if len(r) == 0:
                return "约定", f"{m} 在上游无行，证不出"
            vals = ([r[pr["col"]].iat[0]] if "col" in pr else [r[c].iat[0] for c in pr["cols"]])
            exact = sum(Fraction(str(v)) for v in vals) / len(vals)
            rounded = float(exact)
            hit = [s for s, vv in (("A", sideA_vals[m]), ("B", sideB_vals[m])) if repr(vv) == repr(rounded)]
            if not hit:
                # ── 失败的证法也是结果：说明该量**不是**冻结输入的一次舍入聚合
                return "约定", (f"{m}: 精确均值 {float(exact)!r} 与两侧 double **都不相等** "
                                f"⇒ 该量非一次舍入聚合（管线做了多步舍入），无法从冻结件复现其舍入路径")
            msgs.append(f"{m}: exact={float(exact)!r} == {hit[0]} 侧")
        return "可证", "；".join(msgs)
    except Exception as e:
        return "约定", f"证明过程异常：{e}"


srcs_view = {}
for _a, _b, _DA, _DB, _ca, _cb, _nd in PAIRS:
    srcs_view[_a] = pd.to_numeric(_DA[_ca], errors="coerce")
    srcs_view[_b] = pd.to_numeric(_DB[_cb], errors="coerce")

p("③ **正典生产者登记**（每个量：渲染跟随谁 ｜ 状态）")
p(f"    {'量（两个生产者）':<52}{'正典':<12}{'状态':<6}理由")
for pair, (side, why) in sorted(CANON.items()):
    active = pair not in LATENT
    p(f"    {pair[0] + ' vs ' + pair[1]:<52}{side:<12}{'活跃' if active else '潜伏':<6}{why}")
p("    ⇒ **状态**：活跃 = 两侧都在被渲染（分歧读者可见）｜潜伏 = 仅一侧在渲染。")
p("")

# ================================================================ ⑥ 正典性：**尝试证明**，而不是贴标签
# ★ 2026-10-09（外部核验建议「标注等级」，**做法改了**）：标签不是判据。
#   本模块对**真正分歧的那一对**做精确算术，能证则「可证」，证不出则**说明为什么证不出** ——
#   失败本身也是结果：它说明该量**不是**冻结输入的一次舍入聚合。
PROOF = [
    dict(pair="LODO macro：fig3a_lodo vs p6e/tableS2.macro_mean",
         A=("fig3a_lodo.lodo_macro_auroc", float(G3A.loc["entropy_mean_shift", "lodo_macro_auroc"])),
         B=("tableS2.macro_mean", float(TS2.loc["entropy_mean_shift", "macro_mean"])),
         src=os.path.join(AN, "p6e_lodo_fold_auc_long.csv"), key="method", m="entropy_mean_shift",
         cols=["fold_auc"], note="未加权均值，26 个可评价折"),
    dict(pair="LOO-organism：fig3b.loo_org_macro vs tableS2.loo_org_auc",
         A=("fig3b.loo_org_macro", None), B=("tableS2.loo_org_auc", None),
         src=os.path.join(FIG, "fig3b_loo_organism.csv"), key="method", m=None,
         cols=["loo_auc_EC", "loo_auc_KP", "loo_auc_PA"], note="3 个物种留一的均值"),
]
from fractions import Fraction

p("④ **正典性 —— 逐个尝试证明**（不是贴标签）")
for spec in PROOF:
    src = pd.read_csv(spec["src"])
    diffs = []
    if spec["m"]:                                   # 指定某一 method 的**全部行**（如 26 个折）
        r = src[src[spec["key"]] == spec["m"]]
        assert len(r) > 1, "method=%s 只有 %d 行 —— 折级精确均值无意义" % (spec["m"], len(r))
        exact = sum(Fraction(str(v)) for c in spec["cols"] for v in r[c]) / (len(r) * len(spec["cols"]))
        diffs.append((spec["m"], spec["A"][1], spec["B"][1]))
    else:                                           # 逐行：只取**两个生产者真的不同**的行
        t = TS2
        for _, r in src.iterrows():
            m = r[spec["key"]]
            if m not in t.index or m not in G3B.index:
                continue
            a = float(r["loo_org_macro"]); b = float(t.loc[m, "loo_org_auc"])
            if a != b:
                row = src[src[spec["key"]] == m]
                ex = sum(Fraction(str(row[c].iat[0])) for c in spec["cols"]) / len(spec["cols"])
                diffs.append((m, a, b, ex))
        exact = None
    p(f"    **{spec['pair']}**")
    p(f"      输入：{os.path.basename(spec['src'])} ｜ 聚合：{spec['note']}")
    if spec["m"]:
        got = float(exact)
        hitA = repr(spec["A"][1]) == repr(got)
        hitB = repr(spec["B"][1]) == repr(got)
        p(f"      exact = {exact}")
        p(f"      float(exact) = {got!r}")
        p(f"      A {spec['A'][0]} = {spec['A'][1]!r}  {'← 等于精确值' if hitA else ''}")
        p(f"      B {spec['B'][0]} = {spec['B'][1]!r}  {'← 等于精确值' if hitB else ''}")
        if hitA or hitB:
            p(f"      ⇒ **可证**：正确舍入在 {'A' if hitA else 'B'} 侧，另一侧是多步舍入的产物。")
        else:
            p("      ⇒ **约定**：精确值与两侧都不等 —— 该量非一次舍入聚合，无法从冻结件复现。")
    else:
        n_hit = sum(1 for m, a, b, ex in diffs if repr(a) == repr(float(ex)) or repr(b) == repr(float(ex)))
        p(f"      两侧逐位不同的行：{len(diffs)} ｜ 其中精确值命中某一侧：{n_hit}")
        for m, a, b, ex in diffs[:3]:
            p(f"        {m:<22} A={a!r}  B={b!r}  float(exact)={float(ex)!r}")
        if len(diffs) and n_hit == len(diffs):
            p("      ⇒ **可证**（每个分歧行都有唯一正确侧）")
        else:
            p("      ⇒ **约定（证不出）**：至少一行上，精确均值与两侧 double **都不相等** ——")
            p("         该量由**多步舍入**产生（先加后除、逐步舍入），其舍入路径无法从冻结件复现。")
            p("         **失败本身是结果**：它把「可证」的适用条件收窄为「冻结输入 ＋ **一次**舍入聚合」。")
    p("")

# （原「④ LODO 精确判定」段已由上方**通用化的 ④ 证明段**取代 —— 不再硬编码单个量的结论。）

# ---------------------------------------------------------------- 判词
if WATCH:
    p("⑤ **latent（已登记、当前无读者可见后果）**：%d 格 —— %s" %
      (len(WATCH), "; ".join(sorted({x[1] + " vs " + x[2] for x in WATCH}))))
    p("")
p("=" * 100)
if ALARM:
    p(f"*** p53 **未通过** —— {len(ALARM)} 格：**同一量、两个读者可见的渲染、显示精度不同** ***")
    p("")
    for i, a, b, u, v, nd in ALARM:
        p(f"  · {i}  （行 {i}）")
        p(f"      {a:<34} {u!r}  ->  {('%.*f' % (nd, u))}")
        p(f"      {b:<34} {v!r}  ->  {('%.*f' % (nd, v))}")
    p("")
    p("★ 形态：真值越靠近**进位边界**，同一个 ULP 越容易被放大成一位数字的差。")
    p("  处置方向（择一，须记入事故日志）：① 让两个消费者**读同一个生产者**；")
    p("  ② 让下游在**全精度**上再舍入一次（禁止用已舍入的中间值二次舍入）。")
else:
    p("*** p53 通过 —— 登记的同一量对中，无任何一对在显示精度上分歧 ***")
REP.close()
sys.exit(1 if ALARM else 0)
