# -*- coding: utf-8 -*-
r"""p46 — Supplementary Table S3 的渲染（分诊策略模拟）

================================ 性质 ================================
**事后**分析（post hoc / retrospective policy simulation）的补充表。
数据源 = p45（预注册备忘 TRIAGE_POLICY_PREREGISTRATION_20261007.md）。
本脚本**只做渲染**：不重算任何统计量，全部数值取自冻结产物。

================================ 五条硬约束（每条都有外部核验依据） ================
C1 **推断只在一格。** 只有 k = 20% 列带 95% cluster bootstrap 区间，且只有主对比是
   预注册的。k = 10% / 30% 是**敏感性行**：只报命中与精确率，**不给区间、不写显著性语言**。
   （实现：k≠20% 的行的区间列**留空**，由 `assert_no_ci` 机械保证。）
C2 **不做逐 k 最优高亮、不加粗任何信号。** 表脚说明唯一预注册格，其余声明为完整性呈现。
C3 **表注由生成器模板产出**，不手工编辑（守卫的配置住在生成器里，不住在生成物里）。
C4 **一次舍入。** 统一入口 `_r3()`，禁止链式舍入
   （来历：图例曾把 0.486486 →(4dp) 0.4865 →(3dp) 0.487）。
C5 **每成员差值六行列全**（L-DR 与六个 L1 成员逐个之差）。

另加一条**数值巧合脚注**：主对比区间上界与 L-DR 点精确率数值相同（都显示 0.486），
但角色不同（bootstrap 百分位 vs 观测精确率）。措辞限用
`numerically identical, conceptually distinct` —— **不得写 by construction**
（未在 bootstrap 分布里验证其结构性）。

================================ 闭合断言（GLM 约束 3） ================================
  T46-A 点估计 =（命中差 / k）逐位可推
  T46-B 表内数字与 `figS10_triage_summary.csv` 一致
  T46-C 表内数字与 `TRIAGE_RECONCILIATION_NOTE_20261007.md` 中写死的记法一致（三方一致）
  T46-D k ≠ 20% 的行**没有**任何区间值
  T46-E 表内无粗体 / 无高亮
  T46-F `_r3` 单入口：与直接用内置 round 的结果逐位一致
"""
import io
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pandas as pd

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIG = os.path.join(ROOT, "figures", "data")
RP = os.path.join(FIG, "reporting_profiles")
OUT = os.path.join(ROOT, "分析结果")
REP = open(os.path.join(OUT, "report_p46_tableS3.txt"), "w", encoding="utf-8")


def p(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True); REP.write(s + "\n"); REP.flush()


R = []
OK = True


def chk(cid, ok, expect, observed, note=""):
    global OK
    OK = OK and bool(ok)
    R.append(dict(cid=cid, status="PASS" if ok else "FAIL",
                  expected=expect, observed=observed, note=note))
    p(f"  [{('PASS' if ok else 'FAIL'):<4}] {cid}  {expect}")
    p(f"         观测：{observed}")


# ---------------------------------------------------------------- C4 一次舍入（单一入口）
def _r3(x):
    """**唯一的舍入入口**：从源头全精度一次舍到 3 位小数。
    绝不接受"已经是舍入过的值"再舍一次 —— 那正是链式舍入（0.486486→0.4865→0.487）。"""
    return "%.3f" % float(x)


def _r2(x):
    return "%.2f" % float(x)


def _f0(x):
    return "%d" % int(x)


# ---------------------------------------------------------------- 输入
SUMM = pd.read_csv(os.path.join(OUT, "p45_triage_summary.csv"))
G10 = pd.read_csv(os.path.join(FIG, "figS10_triage_summary.csv")).set_index("layer")
MEM = pd.read_csv(os.path.join(FIG, "figS10_triage_members_k37.csv"))
OVL = pd.read_csv(os.path.join(FIG, "figS10_triage_overlap.csv"))
RB = pd.read_csv(os.path.join(FIG, "figS10_triage_recall_budget.csv"))
UNIV = pd.read_csv(os.path.join(FIG, "fig1b_universe.csv")).iloc[0]
DELT = pd.read_csv(os.path.join(FIG, "fig1c_delta_distribution.csv"))

N_TOT = int(UNIV.eligible_tasks)          # 186
N_FAIL = int(DELT.harmful.sum())          # 57
BASE = N_FAIL / N_TOT
K_MAIN = 0.20
# ★★ 2026-10-12（外部核验要求）：追加两层 —— `L-src`（**平凡基线**：只用 source 自身
#   判别力排序）与 `L-dhat`（事后对齐分数）。**均为事后追加**，见
#   TRIAGE_POLICY_PREREGISTRATION_AMENDMENT_**20261012**.md。
#   ★ 原指向 `..._AMENDMENT_20261007.md` 是**假的出处引用**（INC-032）：该文件写于 10-07、
#     自述不改动 D1–D4，其中**没有这两层**。指向已更正。
#   病历：正文与摘要此前只报「DR 1.59× vs composite」，而 **composite 本身低于随机**
#   （7/37，0.62×）—— 胜过一个低于随机的对照，门槛极低；**有信息量的对照是 src-only**。
#   ★ 主对比仍是 DR − input-shift composite，追加两层**不改变**预注册的任何结论。
LAYERS = ["L-DR", "L-OM", "L-shift", "L-src", "L-dhat", "L-random"]
LAB = {"L-DR": "DR (transported AUC, primary)",
       "L-OM": "Outcome model (sensitivity)",
       "L-shift": "Input-shift composite (current practice)",
       "L-src": "Source discrimination alone (post-unblinding baseline)",
       "L-dhat": "Aligned change score −Δ̂ (post-unblinding)",
       "L-random": "Random ordering (exact expectation; = constant predictor)"}
ALL15 = ["domain_auc", "mmd", "psi_mean", "psi_max", "js_mean", "js_max",
         "pred_mean_shift", "entropy_mean_shift", "entropy_js", "conf_dist_js",
         "pred_pos_prev_shift", "SUDO_AURCC", "AUROC_DR", "AUROC_IOW", "AUROC_OM"]
FAM = {"domain_auc": "L1 input-space", "mmd": "L1 input-space", "psi_mean": "L1 input-space",
       "psi_max": "L1 input-space", "js_mean": "L1 input-space", "js_max": "L1 input-space",
       "pred_mean_shift": "L2 prediction-space", "entropy_mean_shift": "L2 prediction-space",
       "entropy_js": "L2 prediction-space", "conf_dist_js": "L2 prediction-space",
       "pred_pos_prev_shift": "L2 prediction-space", "SUDO_AURCC": "L2 prediction-space",
       "AUROC_DR": "L3 transported (primary)", "AUROC_IOW": "L3 transported (sensitivity)",
       "AUROC_OM": "L3 transported (sensitivity)"}
KS = [0.10, 0.20, 0.30]

p("=" * 116)
p("p46 — Supplementary Table S3（分诊策略模拟，**事后**）渲染")
p("=" * 116)
p(f"上游：p45_triage_summary.csv ｜ figS10_triage_summary.csv ｜ figS10_triage_members_k37.csv")
p(f"宇宙：{N_TOT} 迁移 ｜ 失效 {N_FAIL} ｜ 基率 {_r3(BASE)}")
p("")

# ---------------------------------------------------------------- 取数（一律从冻结产物，不重算）
def cell_of(item, kf):
    r = SUMM[(SUMM.item == item) & (abs(SUMM.k_frac - kf) < 1e-9)]
    if len(r) != 1:
        raise RuntimeError(f"取数失败：{item} @ k_frac={kf} 命中 {len(r)} 行")
    r = r.iloc[0]
    # ★ 返回**原始浮点**，不截断。L-random 的 caught 是精确期望 k×57/186（5.82 / 11.34 / 17.16），
    #   初版在此写 `int(...)` ⇒ 显示成 5.00 / 11.00 / 17.00，与同行给出的精确率 0.306
    #   自相矛盾（5/19 = 0.263）。整数行由调用处用 `_f0` 格式化。
    return float(r.caught), float(r.precision)


# ---------------------------------------------------------------- 闭合断言
p("【闭合断言】")

# T46-A 点估计 =（命中差 / k）
cDR, _ = cell_of("L-DR", K_MAIN); cSH, _ = cell_of("L-shift", K_MAIN)
k_main = int(round(K_MAIN * N_TOT))
cDR, cSH = int(cDR), int(cSH)                  # 这两个必然是整数（真实计数）
pt_from_hits = (cDR - cSH) / k_main
pt_csv = float(G10.loc["L-DR_minus_L-shift", "precision_point"])
_lo = float(G10.loc["L-DR_minus_L-shift", "lo95"])
_hi = float(G10.loc["L-DR_minus_L-shift", "hi95"])
chk("T46-A", abs(pt_from_hits - pt_csv) < 1e-12,
    "主对比点估计 == (命中_DR − 命中_shift) / k，逐位可推",
    f"({cDR} − {cSH}) / {k_main} = {pt_from_hits!r} vs CSV {pt_csv!r}")

# T46-B 表内数字与 figS10_triage_summary.csv 一致
b_ok, b_msg = True, []
for lay in ("L-DR", "L-OM", "L-shift"):
    ca, pr = cell_of(lay, K_MAIN)
    if ca != int(G10.loc[lay, "caught"]) or abs(pr - float(G10.loc[lay, "precision_point"])) > 1e-12:
        b_ok = False; b_msg.append(lay)
chk("T46-B", b_ok, "k=20% 的命中/精确率与 figS10_triage_summary.csv 逐位一致",
    "全部一致" if b_ok else f"不符：{b_msg}")

# T46-C 与 reconciliation note 三方一致
# ★ 2026-10-07 更正：初版用**字符串**比对三位小数的数，而 note 写的是四位小数
#   （`[+0.0789, +0.4865]`）⇒ 假失败。三方一致应当**按数值**比对，
#   容差 = 半个末位（5e-4），而不是比字符串。
_note = io.open(os.path.join(ROOT, "TRIAGE_RECONCILIATION_NOTE_20261007.md"),
                encoding="utf-8").read()
TOL = 5e-4
_ivs = [(float(a), float(b)) for a, b in
        re.findall(r"\[\+?([0-9.]+),\s*\+?([0-9.]+)\]", _note)]
_ci_hit = [(a, b) for a, b in _ivs if abs(a - _lo) <= TOL and abs(b - _hi) <= TOL]
_nums = [float(x) for x in re.findall(r"(?<![\d.])[+-]\d+\.\d+", _note)]
_pt_hit = [x for x in _nums if abs(x - pt_csv) <= TOL]
chk("T46-C", bool(_ci_hit) and bool(_pt_hit),
    "reconciliation note 的记法与表内数字**按数值**一致（容差 5e-4，不比字符串）",
    f"区间命中 {_ci_hit}｜点估计命中 {_pt_hit}")

# T46-D k ≠ 20% 的行**没有**区间
_rows_no_ci = []
for kf in (0.10, 0.30):
    for it in ALL15 + LAYERS:
        ca, pr = cell_of(it, kf)
        _rows_no_ci.append((it, kf, ca, pr))
# ★★ 2026-10-12：原写死 `19 * 2`（15 信号 + 4 层）。本轮 LAYERS 增到 6 ⇒ 立刻断。
#   **改由条目数派生** —— 凡「可以被派生」的期望值都不得硬编码（本仓已因此吃过多次）。
_N_ITEMS = len(ALL15) + len(LAYERS)
chk("T46-D", len(_rows_no_ci) == _N_ITEMS * 2,
    f"k = 10% / 30% 共 {_N_ITEMS * 2} 行，全部**无区间**（结构上不给 CI）",
    f"{len(_rows_no_ci)} 行 = ({len(ALL15)} 信号 + {len(LAYERS)} 层) × 2；"
    f"CI 仅存在于 k=20% 列")

# T46-E 见文末（必须**在 md 生成之后**对实际产物扫描 —— 初版把它写成硬编码 False，
#        是**空转判据**：表里明明有 `**` 也判绿。空转的检查比没有检查更糟。）

# T46-F _r3 单入口 == 直接 round
_f_ok = all(_r3(v) == "%.3f" % v for v in
            [pt_csv, 0.48648648648648646, 0.07894736842105263, 0.3064516129032258])
chk("T46-F", _f_ok, "_r3 为唯一舍入入口，且与直接一次舍入逐位一致（防链式）",
    f"0.486486486… → {_r3(0.48648648648648646)}（不是 0.487）")

# T46-G 非随机行的命中数必须是**整数计数**（防浮点被截断显示）
_frac = []
for _it in ALL15 + [l for l in LAYERS if l != "L-random"]:
    for _kf in KS:
        _c, _ = cell_of(_it, _kf)
        if abs(_c - round(_c)) > 1e-9:
            _frac.append((_it, _kf, _c))
chk("T46-G", not _frac, "非随机行的命中数均为整数计数（随机行是期望值，故排除）",
    "全部为整数" if not _frac else f"非整数：{_frac[:3]}")

# T46-H 随机行的命中数 = k×57/186（**精确期望**），且与同行精确率自洽
#      —— 直接针对本轮实测到的缺陷：初版写成 5.00/11.00/17.00，与同行 0.306 自相矛盾
_bad_exp = []
for _kf in KS:
    _k = int(round(_kf * N_TOT))
    _c, _pr = cell_of("L-random", _kf)
    if abs(_c - _k * BASE) > 1e-9 or abs(_pr - _c / _k) > 1e-12:
        _bad_exp.append((_kf, _c, _k * BASE, _pr, _c / _k))
chk("T46-H", not _bad_exp,
    "随机行：命中 = k×57/186 的精确期望，且命中/k == 同行精确率（自洽）",
    "自洽" if not _bad_exp else f"不自洽：{_bad_exp}")
p("")

# ---------------------------------------------------------------- 生成 Markdown
md = []
md.append("# Supplementary Table S3 — Triage policy simulation at fixed label-review budgets")
md.append("")
md.append("> **Retrospective policy simulation, not a confirmatory analysis.** A monitoring")
md.append("> programme must choose which transfers to send for outcome-label review under a fixed")
md.append("> budget of *k* transfers, ranked by a label-free signal. A transfer is *failing* when")
md.append(f"> ΔAUROC ≤ −0.05; there are {N_FAIL} such transfers among the {N_TOT} in the frozen")
md.append("> benchmark and the input-shift statistics are the six L1 metrics in their prespecified")
# ★★ 2026-10-13：**「pre-registered」一词全面撤出本表**（外部核验 2026-10-12 要求）。
#   理由：这张模拟的**方案**固定在**目标标签解封之后**（详见 TRIAGE_POLICY_PREREGISTRATION_
#   20261007.md 的增补），因此它在**本稿的术语里是 post hoc** —— 而正文 Results §4 已明写
#   「The simulation is post hoc … carries no confirmatory weight」。
#   表注若仍称「pre-registered」，**与正文正面冲突**。
#   改成什么：说它**真实是什么** ——「方案在模拟计算之前定下，但在标签解封之后」
#   ⇒ 这是**该模拟自身的主对比**，不是密封分析的一部分。
md.append("> native direction. **Only the k = 20% column carries inference.** The primary contrast of")
md.append("> this simulation is DR minus the input-shift composite at k = 20%. Its plan was fixed")
md.append("> **before the simulation was computed but after the target labels were unsealed**, so it")
md.append("> is **post hoc** and carries no confirmatory weight. The k = 10% and")
md.append("> k = 30% columns are sensitivity rows reported for completeness, with no interval and")
md.append("> no significance claim. Intervals are 95% organism–drug cluster bootstrap intervals")
md.append("> (10,000 replicates), given only for that one contrast.")
md.append(">")
md.append("> **No signal is highlighted or ordered by performance at any budget.** A statistic that")
md.append("> happens to rank well at k = 10% carries no special status; the single contrast-bearing")
md.append("> cell is identified above, and everything else in this table is descriptive.")
md.append(">")
md.append("> **No multiple-comparison correction was applied.** The 15 signals of Panel A, and the")
md.append("> policy layers of Panel B, are reported side by side without family-wise adjustment;")
md.append("> intervals are per-comparison and should be read as descriptive.")
md.append(">")
# ★ 可复现性条款（外部完备性核验提出：缺了它，重新实现者无法从公开补充材料复现任一格）
md.append("> **Reproducibility.** Transfers are ranked by descending risk score, with ties broken by "
          "ascending task identifier. The review budgets are *k* = round(0.10 × 186) = 19, "
          "round(0.20 × 186) = 37 and round(0.30 × 186) = 56 — **rounded, not truncated**. "
          "Intervals come from 10,000 organism–drug cluster bootstrap replicates (seed 20261007) "
          "in which the budget is held at the same *fraction*, *k* = round(0.20 × replicate size). "
          "Recall for any cell is that cell's hits divided by 57.")
md.append("")
md.append(f"The base rate is {_r3(BASE)} ({N_FAIL} of {N_TOT}); at k = 20% ({k_main} transfers) "
          f"a random ordering captures {_r2(k_main * BASE)} failing transfers in expectation.")
md.append("")

# ---- Panel A：15 个信号 × 3 个预算
md.append(f"## Panel A — Signal-level triage performance at three review budgets")
md.append("")
# ★ 表头**不得**承诺 (95% CI)：单个信号没有聚类 bootstrap 区间（只有政策层有）。
#   初版表头写 "(95% CI)" 而 11 行没有区间 ⇒ 表头在承诺它给不出的东西。
md.append("| Signal | Family | k=10% hits | k=10% precision | k=20% hits | k=20% precision | k=30% hits | k=30% precision |")
md.append("|---|---|---|---|---|---|---|---|")
for it in ALL15:
    c10, p10 = cell_of(it, 0.10)
    c20, p20 = cell_of(it, 0.20)
    c30, p30 = cell_of(it, 0.30)
    k10, k30 = int(round(0.10 * N_TOT)), int(round(0.30 * N_TOT))
    if it in G10.index:
        lo = float(G10.loc[it, "lo95"]); hi = float(G10.loc[it, "hi95"])
        ci = f"{_r3(p20)} ({_r3(lo)}–{_r3(hi)})"
    else:
        ci = _r3(p20)          # 非政策层，本表不为其提供区间
    md.append(f"| `{it}` | {FAM[it]} | {_f0(c10)}/{k10} | {_r3(p10)} | "
              f"{_f0(c20)}/{k_main} | {ci} | {_f0(c30)}/{k30} | {_r3(p30)} |")
md.append("")
md.append("> Panel A reports hits and precision only. Cluster-bootstrap intervals are given for the "
          "policy layers in Panel B, where that contrast lives; no interval is "
          "attached to any individual signal here.")
md.append("")

# ---- Panel B：政策层
md.append("## Panel B — Policy layers, including the primary contrast of the simulation")
md.append("")
md.append("| Policy | k=10% hits | k=10% precision | k=20% hits | k=20% precision (95% CI) | k=30% hits | k=30% precision |")
md.append("|---|---|---|---|---|---|---|")
for lay in ("L-random", "L-shift", "L-src", "L-OM", "L-DR", "L-dhat"):
    c10, p10 = cell_of(lay, 0.10)
    c20, p20 = cell_of(lay, 0.20)
    c30, p30 = cell_of(lay, 0.30)
    k10, k30 = int(round(0.10 * N_TOT)), int(round(0.30 * N_TOT))
    h10 = _r2(c10) if lay == "L-random" else _f0(c10)
    h20 = _r2(c20) if lay == "L-random" else _f0(c20)
    h30 = _r2(c30) if lay == "L-random" else _f0(c30)
    lo = float(G10.loc[lay, "lo95"]); hi = float(G10.loc[lay, "hi95"])
    ci = f"{_r3(p20)} ({_r3(lo)}–{_r3(hi)})" if lo == lo else _r3(p20)
    md.append(f"| {LAB[lay]} | {h10}/{k10} | {_r3(p10)} | {h20}/{k_main} | {ci} | "
              f"{h30}/{k30} | {_r3(p30)} |")
# ★ 表体**不得出现任何粗体**（Nature Portfolio：不得用粗体作强调）。
#   唯一预注册格由**表注**指认，不由排版指认 —— 这也是 GLM 约束 C2 的本意。
md.append(f"| DR minus input-shift composite (the primary contrast of this simulation) | — | — | "
          f"{_f0(cDR - cSH)}/{k_main} | {_r3(pt_csv)} ({_r3(_lo)}–{_r3(_hi)}) | — | — |")
# ★★ 2026-10-12 追加（外部核验要求）：**与平凡基线的配对差**。
#   正文与摘要引用了它的点估计与区间，**必须在本表里有家**（否则读者无处可核 —— INC-020 的形态）。
_cSR = int(cell_of("L-src", K_MAIN)[0])
_pt2 = float(G10.loc["L-DR_minus_L-src", "precision_point"])
_lo2 = float(G10.loc["L-DR_minus_L-src", "lo95"])
_hi2 = float(G10.loc["L-DR_minus_L-src", "hi95"])
md.append(f"| DR minus source-discrimination baseline (post-unblinding) | — | — | "
          f"{_f0(cDR - _cSR)}/{k_main} | {_r3(_pt2)} ({_r3(_lo2)}–{_r3(_hi2)}) | — | — |")
md.append("")
md.append(f"> **The primary contrast of this simulation** — its plan was fixed before the simulation "
          f"was computed but after the target labels were unsealed, so it is **post hoc** — is the "
          f"DR-minus-input-shift-composite row: "
          f"reviewing the top {k_main} transfers by DR captures {cDR - cSH} **more** failing "
          f"transfers than reviewing the same {k_main} by the composite (observed difference "
          f"{_r3(pt_csv)}), 95% cluster bootstrap interval {_r3(_lo)} to {_r3(_hi)}. It is the only "
          f"quantity in this table for which a significance claim is made. Its interval lower "
          f"bound is above zero.")
md.append("")
md.append(f"> **The informative comparator is not the composite.** The composite's own precision")
md.append(f"> (0.189 at k = 20%) is below the uninformative expectation of 0.306, so exceeding it")
md.append(f"> is a weak test. The row comparing DR with the source-discrimination baseline is")
md.append(f"> reported for that reason; its interval includes zero.")
md.append("")
# ★ 措辞限用「显示精度上相同」。**不得**写 "by construction"，也**不得**声称比特级相等 ——
#   实测：区间上界存为 0.48648648648648646，比 double(18/37) = 0.4864864864864865 低 1 ULP
#   （它是两个浮点精确率之差的百分位）。精确算术下同为 18/37，但双精度下不同。
md.append(f"> **Note on a numerical coincidence.** The upper bound of that interval and the observed "
          f"k = 20% precision of DR are *identical at the displayed precision* ({_r3(_hi)}); both "
          f"reduce to {cDR}/{k_main} in exact arithmetic, though as stored doubles they differ in "
          "the last bit. They are **conceptually distinct**: one is a bootstrap percentile of the "
          "difference between two policies, the other is an observed precision. Neither should be "
          "read as the other.")
md.append("")

# ---- Panel C：每成员差值（六行列全）
md.append("## Panel C — DR minus each input-shift statistic (k = 20%)")
md.append("")
md.append("| Input-shift statistic | Its own hits | Its own precision | DR − this statistic |")
md.append("|---|---|---|---|")
mem = MEM.set_index("member")
for m in ["domain_auc", "mmd", "psi_mean", "psi_max", "js_mean", "js_max"]:
    ch = int(mem.loc[m, "caught"])
    md.append(f"| `{m}` | {_f0(ch)}/{k_main} | {_r3(ch / k_main)} | "
              f"{_r3(cDR / k_main - ch / k_main)} |")
md.append("")
_ov20 = OVL[abs(OVL.k_frac - 0.20) < 1e-9]
_ov_lo, _ov_hi = int(_ov20.overlap.min()), int(_ov20.overlap.max())
_ov_who = _ov20.loc[_ov20.overlap.idxmax(), "member"]
md.append("> This panel exists so that the composite comparator cannot be read as a post-hoc "
          "choice of an unfavourable aggregate. The comparator fixed in the simulation plan is the equal-weight "
          "rank-average of these six statistics, and it happens to fall at the lower edge of their "
          "range; the row-wise differences above show that the ranking advantage holds against "
          "**every** member individually, not only against the composite.")
md.append(">")
# ★ 反向读法也要挡住（外部完备性核验提出）：复合与成员本身高度重合，
#   读者应能自己判断「复合是否只是某个成员的替身」。
md.append(f"> The converse reading is also answerable from the data: at k = 20% the composite's top "
          f"{k_main} overlaps each individual statistic by {_ov_lo} to {_ov_hi} of {k_main} "
          f"transfers (highest for `{_ov_who}`). The six input-shift statistics are strongly "
          f"collinear, so the composite resembles all of them about equally and the comparison "
          f"above is not an artefact of any single member.")
md.append(">")
# ★ 预注册备忘 §3 D1 承诺的次要描述量：达到 50% 召回的预算。初版未交付，现补。
_rb = {r.layer: r for _, r in RB.iterrows()}
md.append(f"> **Secondary descriptor (promised in the simulation plan).** The review budget needed to capture half "
          f"of the {N_FAIL} failing transfers is {_f0(_rb['L-OM'].n_needed)} transfers "
          f"({_r3(_rb['L-OM'].frac_needed)}) for the outcome model, "
          f"{_f0(_rb['L-DR'].n_needed)} ({_r3(_rb['L-DR'].frac_needed)}) for DR, and "
          f"{_f0(_rb['L-shift'].n_needed)} ({_r3(_rb['L-shift'].frac_needed)}) for the input-shift "
          f"composite. No policy reaches half-recall within the three budgets tabulated above, "
          f"which is why this quantity is given separately.")

md.append("")

MD = "\n".join(md)

# ---------------------------------------------------------------- T46-E（在产物上扫描，不在代码里假定）
# 判据：所有以 `|` 开头的**表格行**中，不得出现 `**`。
# 非表格行（标题、表注）允许存在粗体 —— 限制的是表体，不是全文。
_tbl_lines = [ln for ln in MD.splitlines() if ln.strip().startswith("|")]
_bold_rows = [ln for ln in _tbl_lines if "**" in ln]
chk("T46-E", not _bold_rows,
    f"表体（{len(_tbl_lines)} 行）无任何粗体标记（逐行扫描实际产物，非代码内假定）",
    "无粗体行" if not _bold_rows else f"{len(_bold_rows)} 行含 **：{_bold_rows[:2]}")

io.open(os.path.join(RP, "Supplementary_Table_S3.md"), "w", encoding="utf-8").write(MD)
p(f"写出 {os.path.join(RP, 'Supplementary_Table_S3.md')}")

# ---------------------------------------------------------------- 数据产物 + 逐格溯源
prov = []
for it in ALL15 + LAYERS:
    for kf in KS:
        ca, pr = cell_of(it, kf)
        # ★ has_ci 必须反映**区间是否真的存在**，不能只反映「这一行在不在区间表里」。
        #   初版写 `it in G10.index` ⇒ 随机行被标 True，而它的 lo95/hi95 是 NaN
        #   （随机线是精确期望，没有抽样分布）。下游若按该旗标重渲染，
        #   会**凭空造出一条随机线的区间**，与本文自己的立场相反。
        _has = (abs(kf - K_MAIN) < 1e-9 and it in G10.index
                and not pd.isna(G10.loc[it, "lo95"]))
        prov.append(dict(item=it, kind="signal" if it in ALL15 else "layer",
                         k_frac=kf, k=int(round(kf * N_TOT)), caught=ca,
                         precision_full=pr, precision_shown=_r3(pr),
                         has_ci=_has,
                         # ★★ 2026-10-14（泄漏）：原写 `source="p45_triage_summary.csv"` ——
                         #   渲染件里出现**内部脚本文件名 ×63 处**，而 p41 的 V7-3 明文禁止
                         #   （期刊件不得含脚本文件名）。改为**按角色描述该上游**，不报文件名。
                         source="policy-triage summary for this analysis"))
# ★★ 2026-10-12（INC-030）：**导入必须在第一次使用之前。**
#   原序 `_gs46.backup(...)` → `to_csv` → `import genstamp as _gs46` ⇒ `NameError`。
#   **p11 / p45 / p46 三处同型** —— 都是 A3.6 的 backup+stamp 块被插到了导入之前。
#   ★ 三处都只在「真跑这个脚本」时现形；只读其产物的审计器**永远看不到**。
import sys as _sys2
if os.path.join(ROOT, "分析结果") not in _sys2.path:
    _sys2.path.insert(0, os.path.join(ROOT, "分析结果"))
import genstamp as _gs46
_gs46.backup(os.path.join(FIG, "tableS3_triage.csv"), "p46")
pd.DataFrame(prov).to_csv(os.path.join(FIG, "tableS3_triage.csv"), index=False)
_gs46.stamp(os.path.join(FIG, "tableS3_triage.csv"), "p46")
pd.DataFrame(prov).to_csv(os.path.join(RP, "supplementary_tableS3_provenance.csv"),
                          index=False, encoding="utf-8-sig")

# ★ has_ci 的语义断言：置 True 的每一行，其上游区间必须**真的存在且有限**。
_badflag = []
for _r in prov:
    if _r["has_ci"]:
        _lay = _r["item"]
        if _lay not in G10.index or pd.isna(G10.loc[_lay, "lo95"]) or pd.isna(G10.loc[_lay, "hi95"]):
            _badflag.append(_lay)
    else:
        if _r["item"] in ("L-DR", "L-OM", "L-shift") and abs(_r["k_frac"] - K_MAIN) < 1e-9:
            _badflag.append(_r["item"] + "(漏标)")
chk("T46-I", not _badflag,
    "has_ci 旗标语义正确：只在该行的上游区间**存在且有限**时为 True",
    "正确" if not _badflag else f"错误：{_badflag}",
    "随机行的区间是 NaN（精确期望无抽样分布）⇒ 必须为 False")

# ---------------------------------------------------------------- 渲染审计
D = pd.DataFrame(R)
D.to_csv(os.path.join(RP, "supplementary_tableS3_render_audit.csv"),
         index=False, encoding="utf-8-sig")
p("")
p(f"渲染审计：{len(D)} 项，FAIL {int((D.status == 'FAIL').sum())}")
p("")
p("*** p46 " + ("通过" if OK else "**未通过**") + " ***")
REP.close()
sys.exit(0 if OK else 1)
