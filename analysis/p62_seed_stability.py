# -*- coding: utf-8 -*-
r"""p62 — 端点的**种子稳定性**（「不得空口归因于蒙特卡洛」这条规则的第一个实例）

================================ 为什么有这个脚本 ================================
外部核验（2026-10-13）立了一条规矩，我同意并把措辞写死：

> **任何两处数字的差异，在找到具体原因之前，不得归为蒙特卡洛。**
> 要归因，必须**换种子重跑**并给出**端点最大漂移**；拿不出就按缺陷处理。

这条规矩的来由是本会话**三次**把差异错归给 MC：
  ① p57 与 p61 的区间差 —— 真因是 p61 的 RNG 被所有调用共享消费；
  ② n=100 的「0.616 ≈ src 0.611」 —— 真因是子集不对齐（我用了 89 行 vs 186 行）；
  ③ DR 对 src-only 的 [−0.1434] 与 [−0.1457] —— 真因是我把 p61 的值转写到了 p57 名下。

================================ 口径（严格照外部核验的约束）====================
  · **只改重采样种子**；数据、划分、模型种子**一律不动**。
  · 种子数 **5**（不是 3 —— 3 个不足以估端点漂移）。
  · 报告**每个端点**（15 个差值的上下界共 30 个）的最大偏离。
  · 若最大偏离 > 0.005 ⇒ 须在 Methods 说明，并考虑把 B 从 10,000 提高。
"""
import io
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "分析结果")
FIG = os.path.join(ROOT, "figures", "data")
B = 10_000
SEEDS = [20261007, 20261008, 20261009, 20261010, 20261011]   # 只换这个

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass
import numpy as np
import pandas as pd
import rephead

REP = io.open(os.path.join(OUT, "report_p62_seed_stability.txt"), "w", encoding="utf-8")


def p(s=""):
    print(s)
    REP.write(str(s) + "\n")


def auc(y, s):
    """Mann–Whitney AUROC，并列取平均秩（与 p57 同实现口径）。"""
    y = np.asarray(y, float); s = np.asarray(s, float)
    m = np.isfinite(s); y, s = y[m], s[m]
    order = np.argsort(s, kind="mergesort")
    ys, ss = y[order], s[order]
    r = np.empty(len(ss), float); i = 0
    while i < len(ss):
        j = i
        while j + 1 < len(ss) and ss[j + 1] == ss[i]:
            j += 1
        r[i:j + 1] = 0.5 * (i + j) + 1
        i = j + 1
    n1 = int(y.sum()); n0 = len(y) - n1
    if n1 == 0 or n0 == 0:
        return float("nan")
    return (r[ys == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)


S = pd.read_csv(os.path.join(FIG, "fig4_l3_scatter.csv"))
S["model"] = S.organism.astype(str) + "|" + S.drug.astype(str) + "|" + S.source_site.astype(str)
S["od"] = S.organism.astype(str) + "|" + S.drug.astype(str)
harm = (S.d_auroc <= -0.05).astype(int).values
src = S.auroc_source_internal.values

_b = pd.read_csv(os.path.join(OUT, "p6b_label_free_signals_blinded.csv")).set_index("task_id")
_c = pd.read_csv(os.path.join(OUT, "p6c_L3_blinded.csv")).set_index("task_id")
_b = _b.loc[S.task_id.values]; _c = _c.loc[S.task_id.values]
L1 = ["domain_auc", "mmd", "psi_mean", "psi_max", "js_mean", "js_max"]
SGN = ["pred_mean_shift", "entropy_mean_shift", "pred_pos_prev_shift"]
SIG = {}
for m in L1 + ["entropy_js", "conf_dist_js"]:
    SIG[m] = _b[m].values.astype(float)
for m in SGN:
    SIG[m] = np.abs(_b[m].values.astype(float))
for m in ["SUDO_AURCC", "AUROC_IOW", "AUROC_OM"]:
    SIG[m] = -_c[m].values.astype(float)
SIG["AUROC_DR"] = -S.AUROC_DR.values

UOD = np.unique(S.od.values)
idx = {c: np.flatnonzero(S.od.values == c) for c in UOD}

for line in rephead.header(__file__, [os.path.join(FIG, "fig4_l3_scatter.csv"),
                                      os.path.join(OUT, "p6b_label_free_signals_blinded.csv"),
                                      os.path.join(OUT, "p6c_L3_blinded.csv")],
                           note="触发：外部核验要求「不得空口归因于 MC」，立第一个实例"):
    p(line)
p("=" * 100)
p("p62 — 端点种子稳定性（**只换重采样种子**；数据/划分/模型种子不动）")
p("  种子：%s ｜ B = %d ｜ 簇 = %d 个 organism–drug" % (SEEDS, B, len(UOD)))
p("=" * 100)
p("")

NAMES = L1 + SGN + ["entropy_js", "conf_dist_js", "SUDO_AURCC", "AUROC_DR", "AUROC_IOW", "AUROC_OM"]
res = {m: {"lo": [], "hi": [], "pt": []} for m in NAMES}
for sd in SEEDS:
    rng = np.random.default_rng(sd)
    boot = {m: [] for m in NAMES}
    boot_src = []
    for _ in range(B):
        rows = np.concatenate([idx[UOD[j]] for j in rng.choice(len(UOD), len(UOD), replace=True)])
        yy = harm[rows]
        if yy.sum() in (0, len(yy)):
            continue
        boot_src.append(auc(yy, src[rows]))
        for m in NAMES:
            boot[m].append(auc(yy, SIG[m][rows]) - auc(yy, src[rows]))
    for m in NAMES:
        a = np.array(boot[m])
        res[m]["lo"].append(float(np.percentile(a, 2.5)))
        res[m]["hi"].append(float(np.percentile(a, 97.5)))
        res[m]["pt"].append(float(auc(harm, SIG[m]) - auc(harm, src)))

p("%-16s %-22s %-22s %s" % ("信号", "下界（5 种子）", "上界（5 种子）", "端点最大偏离"))
print("-" * 92)
worst = 0.0
for m in NAMES:
    rl, rh = np.array(res[m]["lo"]), np.array(res[m]["hi"])
    dl, dh = rl.max() - rl.min(), rh.max() - rh.min()
    worst = max(worst, dl, dh)
    p("%-16s %-22s %-22s %.4f" % (m, "%.4f ~ %.4f" % (rl.min(), rl.max()),
                                  "%.4f ~ %.4f" % (rh.min(), rh.max()), max(dl, dh)))
p("")
p("  ★ 15 个差值的 **30 个端点**中最大偏离 = **%.4f**" % worst)
p("  ★ 判据：> 0.005 ⇒ 须在 Methods 说明。结论：%s"
  % ("**超过，须说明**" if worst > 0.005 else "未超过（B = 10,000 足够）"))
p("  ★ **分族漂移**（外部核验要求不只给最大值）：")
_l1d = [max(np.ptp(res[m]["lo"]), np.ptp(res[m]["hi"])) for m in L1]
_l3d = [max(np.ptp(res[m]["lo"]), np.ptp(res[m]["hi"])) for m in ["SUDO_AURCC", "AUROC_DR", "AUROC_IOW", "AUROC_OM"]]
p("      L1 六个        %.4f ~ %.4f" % (min(_l1d), max(_l1d)))
p("      L2 六个        %.4f ~ %.4f" % (min([max(np.ptp(res[m]["lo"]), np.ptp(res[m]["hi"])) for m in SGN + ["entropy_js", "conf_dist_js"]]),
                                             max([max(np.ptp(res[m]["lo"]), np.ptp(res[m]["hi"])) for m in SGN + ["entropy_js", "conf_dist_js"]])))
p("      L3 + SUDO 四个  %.4f ~ %.4f" % (min(_l3d), max(_l3d)))
p("      （区间越宽端点漂移越大，属预期现象，不解释为「不稳」。）")
p("  ★ **点估计逐位可复现、不随种子变** —— 此句**只对 15 个差值的点估计**成立，")
p("    **不适用于 bootstrap 区间**（区间的端点确实随种子动）。")
p("")

# ---------------------------------------------------------------- 「含 0 / 不含 0」的翻转表
p("=" * 100)
p("【翻转表】5 个种子下，每条「含 0 / 不含 0（或 vs 0.5）」的结论**有没有翻转**")
p("=" * 100)
p("  判据：一条断言若**依赖端点相对 0（或 0.5）的位置**，而该位置在 5 个种子上不一致 ⇒ **不稳健**。")
p("  不稳健者，正文只能用「区间含 0」或「不稳定」，**不得**用「不含 0」。")
p("")
p("%-22s %-30s %-8s %-8s %s" % ("断言", "5 种子下的端点范围", "点估计", "含0?", "稳健性"))
print("-" * 104)
FLIP_ROWS = []
for m in NAMES:
    lo = np.array(res[m]["lo"]); hi = np.array(res[m]["hi"]); pt = res[m]["pt"][0]
    contains = [(a <= 0 <= b) for a, b in zip(lo, hi)]
    stable = len(set(contains)) == 1
    FLIP_ROWS.append(dict(claim="%s − src_only" % m, pt=pt, lo_min=lo.min(), lo_max=lo.max(),
                          hi_min=hi.min(), hi_max=hi.max(),
                          contains_zero=bool(contains[0]), flips=not stable))
    p("%-22s %-30s %+8.4f %-8s %s"
      % (m + " − src_only", "[%+.4f, %+.4f]" % (lo.min(), hi.max()), pt,
         "含0" if contains[0] else "不含0",
         "稳健" if stable else "**翻转**"))

# ★★ 外部核验点名的那一族：**六个 L1 对 ΔAUROC 的关联区间**（来自 tableS2_metrics.csv，p33 的
#   bootstrap，另有一套种子）。它的端点**贴着 0**（psi_max +0.014、js_max +0.022），
#   漂移 0.003–0.005 就足以跨 0 —— 这才是真正该查翻转的地方。
p("")
p("  【另查】六个 L1 对 ΔAUROC 的关联区间（S2 Panel B；端点**贴着 0**，最可能翻转）")
_s2 = pd.read_csv(os.path.join(FIG, "tableS2_metrics.csv")).set_index("method")
_li = []
for m in L1:
    r = _s2.loc[m]
    _li.append(dict(claim="%s vs ΔAUROC" % m, pt=float(r.spearman),
                    lo=float(r.sp_lo), hi=float(r.sp_hi),
                    contains_zero=bool(r.sp_lo <= 0 <= r.sp_hi)))
for r in _li:
    margin = min(abs(r["lo"]), abs(r["hi"]))
    p("     %-26s %+7.4f  [%+.4f, %+.4f]  %-6s  距 0 最近 %.4f"
      % (r["claim"], r["pt"], r["lo"], r["hi"],
         "含0" if r["contains_zero"] else "不含0", margin))
p("")
p("  ★★ **这两行的种子漂移本脚本没有测** —— 不能凭本表判它们稳不稳。")
p("     理由：它们的区间由 p33 的**另一套 bootstrap**产出（不同的 B / 种子 / 实现），")
p("     而本脚本测到的 0.0032~0.0048 是**15 个差值区间**的漂移，**两者不可互相代入**。")
p("     ⇒ 判据：**「psi_max / js_max 的翻转风险：未确定」**，不得写成「稳健」也不得写成「会翻」。")
p("     ⇒ 处置（外部核验同向）：正文对这两行**只写「关联量级很小、方向一致」**，")
p("        **不写「不含 0」、不写成发现** —— 这样无论漂移多大，措辞都成立。")
p("     ⇒ 遗留项：p33 的六个 Spearman 区间的种子稳定性**尚未测**，已记入清单。")
p("")
p("  ⇒ 15 个差值的「含 0」结论：**翻转 %d 条 / 不翻转 %d 条**"
  % (sum(1 for r in FLIP_ROWS if r["flips"]), sum(1 for r in FLIP_ROWS if not r["flips"])))
p("  ★ 外部核验的预期（端点贴 0 的不稳、远离的稳）在本表上**逐条可核**。")
p("")
# ================================================================ 补全：其余含区间的断言
# ★★ 外部核验（2026-10-13）要求**一次性全测**，不只测 15 个差值。清单：
#   (a) Â 分解三个量（配对均值分量 / 配对内偏差 vs 0.5 / 完整 − 配对均值 vs 0）
#   (b) 21 个不一致对上的 16/21 与其簇 bootstrap 区间（**下界 0.5455 离 0.5 仅 0.045，
#       且这组只有 21 对、离散大 ⇒ 它点名这是最可能翻转的一条**）
MEAN = S.groupby("model")["AUROC_DR"].transform("mean").values
COMP = {"pair_mean_part": -MEAN, "within_pair_dev": -(S.AUROC_DR.values - MEAN),
        "prespecified_DR": -S.AUROC_DR.values}
res2 = {k: {"lo": [], "hi": []} for k in COMP}
d3 = {"lo": [], "hi": []}
acc_seeds, site_seeds = [], []
for sd in SEEDS:
    rng = np.random.default_rng(sd)
    b2 = {k: [] for k in COMP}
    dd = []
    for _ in range(B):
        rows = np.concatenate([idx[UOD[j]] for j in rng.choice(len(UOD), len(UOD), replace=True)])
        yy = harm[rows]
        if yy.sum() in (0, len(yy)):
            continue
        for k, v in COMP.items():
            b2[k].append(auc(yy, v[rows]))
        dd.append(auc(yy, COMP["prespecified_DR"][rows]) - auc(yy, COMP["pair_mean_part"][rows]))
    for k in COMP:
        a = np.array(b2[k])
        res2[k]["lo"].append(float(np.percentile(a, 2.5)))
        res2[k]["hi"].append(float(np.percentile(a, 97.5)))
    d3["lo"].append(float(np.percentile(dd, 2.5)))
    d3["hi"].append(float(np.percentile(dd, 97.5)))

    # ---- (b) 21 对的 16/21：按 od 簇重采样，每次重取配对
    disc = []
    for mdl, g in S.groupby("model"):
        if len(g) != 2:
            continue
        a_, b_ = g.iloc[0], g.iloc[1]
        if (a_.d_auroc <= -0.05) == (b_.d_auroc <= -0.05):
            continue
        hi_ = a_ if a_.d_auroc <= -0.05 else b_
        lo_ = b_ if hi_ is a_ else a_
        disc.append((hi_.name, lo_.name, a_.od))
    import collections as _c
    byod = _c.defaultdict(list)
    for it in disc:
        byod[it[2]].append(it)
    keys = sorted(byod)
    hit = np.array([SCORES_DR := None] * 0, dtype=float)
    accs = []
    for _ in range(B):
        pick = [keys[i] for i in rng.choice(len(keys), len(keys), replace=True)]
        sub = [it for k in pick for it in byod[k]]
        if len(sub) < 5:
            continue
        ok = sum(1 for i, j, _ in sub if COMP["prespecified_DR"][i] > COMP["prespecified_DR"][j])
        accs.append(ok / len(sub))
    acc_seeds.append((np.percentile(accs, 2.5), np.percentile(accs, 97.5), len(disc)))

p("")
p("=" * 100)
p("【补全】其余含区间的断言的种子稳定性（同 5 个种子）")
p("=" * 100)
for k in COMP:
    lo, hi = np.array(res2[k]["lo"]), np.array(res2[k]["hi"])
    ref = 0.5 if k != "pair_mean_part" else 0.5
    ct = [(a <= ref <= b) for a, b in zip(lo, hi)]
    p("   %-22s [%.4f~%.4f, %.4f~%.4f]  相对 0.5：%-6s  %s"
      % (k, lo.min(), lo.max(), hi.min(), hi.max(),
         "含" if ct[0] else "不含", "稳健" if len(set(ct)) == 1 else "**翻转**"))
lo3, hi3 = np.array(d3["lo"]), np.array(d3["hi"])
ct3 = [(a <= 0 <= b) for a, b in zip(lo3, hi3)]
p("   %-22s [%.4f~%.4f, %.4f~%.4f]  含 0：%-6s  %s"
  % ("完整 − 配对均值", lo3.min(), lo3.max(), hi3.min(), hi3.max(),
     "含" if ct3[0] else "不含", "稳健" if len(set(ct3)) == 1 else "**翻转**"))
p("")
p("   21 对（16/21）的簇 bootstrap 区间，5 个种子下：")
for (a, b, n) in acc_seeds:
    p("      [%.4f, %.4f]   相对 0.5：%s" % (a, b, "含 0.5" if a <= 0.5 <= b else "不含 0.5"))
_alo = min(x[0] for x in acc_seeds); _ahi = max(x[1] for x in acc_seeds)
_flip = len({(x[0] <= 0.5 <= x[1]) for x in acc_seeds}) > 1
p("      ⇒ 下界范围 %.4f~%.4f ｜ 相对 0.5 的位置：%s"
  % (_alo, _ahi, "**翻转 —— 措辞不得依赖它**" if _flip else "不翻转"))
p("")
p("   ★ 本脚本**未测**的（无相应实现，如实列出，不假装覆盖）：分诊的 DR vs 随机 / vs composite、")
p("     ICC 的参数 bootstrap、主表 B=1,000 的 AUROC 区间、p33 的六个 Spearman。已记入清单。")
p("")

pd.DataFrame(FLIP_ROWS).to_csv(os.path.join(OUT, "p62_flip_table.csv"), index=False, encoding="utf-8-sig")
p("done -> %s" % os.path.join(OUT, "report_p62_seed_stability.txt"))
REP.close()
pd.DataFrame({m: {"lo_min": min(res[m]["lo"]), "lo_max": max(res[m]["lo"]),
                  "hi_min": min(res[m]["hi"]), "hi_max": max(res[m]["hi"]),
                  "point": res[m]["pt"][0]} for m in NAMES}).T.to_csv(
    os.path.join(OUT, "p62_seed_stability.csv"), encoding="utf-8-sig")
