# -*- coding: utf-8 -*-
r"""p52 — L3 传输估计量的 MAE 对「最优常数预测器」的**配对**增益：区间

================================ 为什么需要 ================================
Results §4 现在写着：「The transported estimators therefore carry a **small but consistent
advantage** over predicting no change at all, and a much smaller one over the best constant
available.」

「**consistent**」是一个**需要区间支撑**的词 —— 而全项目**没有任何产物**给出
DR/OM 的 MAE 与常数基线之差的区间（`p32_constant_baselines.csv` 只有点估计：
gain_vs_median_floor = 0.003501 / 0.003452 / −0.008885）。

相邻句子已经建立先例：DR 与 OM 的配对差**报了区间**（spanning zero）。同一节里
一个差报了区间、另一个没有，读者会默认后者也经过了同样的检验 —— 这正是
本仓反复记录的「**标签承诺的保证等级高于交付**」。

故本脚本补这一个区间。**不新增任何数据**：上游仍是 `fig4_l3_scatter.csv`。

================================ 配对的定义（写死）================================
· 常数预测器**不是**在每次重采样里重新拟合的 —— 它预测**全样本观测中位数**
  （−0.0242587），文中口径即「predicting the observed median change (−0.024)」。
  重采样只改变**评估样本**，不改变这个常数。
· 两个 MAE 在**同一次重采样内同时计算**，其差即配对统计量。
  （INC 教训：曾把配对量写成两个独立累加器，见 p45 的病历。）
· cluster bootstrap：按 (organism, drug) 重采样 37 个簇，B = 10,000，seed = 20261007
  （与 p33 / p45 / Table S3 同种子）。

================================ 判据 ================================
先**逐位锚定**三个 MAE 点估计与 `p32_constant_baselines.csv` 一致（1e-12），
不一致即中止 —— 本脚本不得成为第二个 MAE 生产者。
"""
import io
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
OUT = os.path.join(ROOT, "分析结果")
REP = open(os.path.join(OUT, "report_p52_mae_vs_constant.txt"), "w", encoding="utf-8")


def p(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    REP.write(s + "\n")
    REP.flush()


B = 10_000
SEED = 20261007
EST = ["DR", "IOW", "OM"]

import rephead
for _l in rephead.header(__file__, [os.path.join(FIG, "fig4_l3_scatter.csv"), os.path.join(OUT, "p32_constant_baselines.csv")]):
    p(_l)
p("=" * 100)
p("p52 — MAE 对最优常数预测器的配对增益（organism–drug cluster bootstrap）")
p("=" * 100)

# ---------------------------------------------------------------- 上游
S = pd.read_csv(os.path.join(FIG, "fig4_l3_scatter.csv"))
need = ["task_id", "organism", "drug", "d_auroc"] + ["d_hat_" + e for e in EST]
miss = [c for c in need if c not in S.columns]
assert not miss, "fig4_l3_scatter.csv 缺列：%s" % miss
S = S[need].copy()
S["cluster"] = S["organism"].astype(str) + "|" + S["drug"].astype(str)
p(f"上游：fig4_l3_scatter.csv  {len(S)} 行 ｜ 簇 {S.cluster.nunique()} 个")

CONST = float(S.d_auroc.median())
p(f"常数预测器所预测的值 = 全样本观测中位数 ΔAUROC = {CONST!r}")
p("")

# ---------------------------------------------------------------- 点估计 + 逐位锚定
def mae(delta):
    """预测值 = 常数预测器的固定值 + delta（delta=0 即常数本身）。"""
    return float(np.abs(S.d_auroc - (CONST + delta)).mean())


CB = pd.read_csv(os.path.join(OUT, "p32_constant_baselines.csv")).set_index("name")
POINT = {e: float(np.abs(S.d_auroc - S["d_hat_" + e]).mean()) for e in EST}
POINT["constant_median"] = mae(0.0)
POINT["constant_zero"] = float(S.d_auroc.abs().mean())

p("① 点估计锚定（与 p32_constant_baselines.csv 逐位比对，tol=1e-12）")
ok = True
for k, col in [("constant_zero", "mae"), ("constant_median", "mae")] + [(e, "mae") for e in EST]:
    ref = float(CB.loc[k, col])
    d = abs(POINT[k] - ref)
    good = d < 1e-12
    ok &= good
    p(f"    {k:<18}{POINT[k]:.12f}  vs p32 {ref:.12f}   {'✔' if good else '✘ 差 %.2e' % d}")
assert ok, "点估计与 p32_constant_baselines.csv 不一致 —— 中止（不得成为第二个生产者）"
p("    ⇒ 本脚本与 p32 同源，非第二个 MAE 生产者")
p("")

# ---------------------------------------------------------------- 簇 bootstrap
clusters = S.cluster.unique()
idx_by_cluster = {c: np.flatnonzero(S.cluster.values == c) for c in clusters}
rng = np.random.default_rng(SEED)

d_obs = S.d_auroc.values
d_hat = {e: S["d_hat_" + e].values for e in EST}

acc = {e: [] for e in EST}          # MAE_est − MAE_median（配对）
acc0 = {e: [] for e in EST}         # MAE_est − MAE_zero（配对）
for _ in range(B):
    pick = rng.choice(len(clusters), size=len(clusters), replace=True)
    rows = np.concatenate([idx_by_cluster[clusters[j]] for j in pick])
    y = d_obs[rows]
    base = np.abs(y - CONST).mean()
    base0 = np.abs(y).mean()
    for e in EST:
        m = np.abs(y - d_hat[e][rows]).mean()
        acc[e].append(m - base)
        acc0[e].append(m - base0)

p(f"② 配对差，B = {B:,}，seed = {SEED}")
p(f"    {'估计量':<8}{'口径':<12}{'点估计':>12}{'bootstrap 均值':>16}{'95% 区间':>26}   是否含 0")
RES = []
for e in EST:
    row = dict(estimator=e, n=len(S), mae=POINT[e],
               mae_constant_median=POINT["constant_median"], mae_constant_zero=POINT["constant_zero"])
    for tag, a_, pt in (("median", acc[e], POINT[e] - POINT["constant_median"]),
                        ("zero", acc0[e], POINT[e] - POINT["constant_zero"])):
        a = np.array(a_)
        lo, hi = np.percentile(a, [2.5, 97.5])
        spans = bool(lo <= 0.0 <= hi)
        row["gain_vs_%s_point" % tag] = pt
        row["gain_vs_%s_lo95" % tag] = float(lo)
        row["gain_vs_%s_hi95" % tag] = float(hi)
        row["gain_vs_%s_spans_zero" % tag] = spans
        p(f"    {e:<8}{'vs ' + tag:<12}{pt:>+12.6f}{a.mean():>+16.6f}   [{lo:+.6f}, {hi:+.6f}]   "
          f"{'**含 0**' if spans else '不含 0'}")
    RES.append(row)
p("")

R = pd.DataFrame(RES)
R.to_csv(os.path.join(OUT, "p52_mae_vs_constant.csv"), index=False, encoding="utf-8-sig")

# ---------------------------------------------------------------- 判词
p("=" * 100)
p("判词")
p("")
for _, r in R.iterrows():
    for tag, lab in (("median", "最优常数（观测中位数）"), ("zero", "预测零变化")):
        pt = r["gain_vs_%s_point" % tag]
        lo = r["gain_vs_%s_lo95" % tag]
        hi = r["gain_vs_%s_hi95" % tag]
        spans = r["gain_vs_%s_spans_zero" % tag]
        if spans:
            p(f"· **{r.estimator}** vs {lab}：MAE 增益 {pt:+.4f}，95% 区间 [{lo:+.4f}, {hi:+.4f}]"
              f" —— **含 0**。")
            p("  ⇒ 该增益**不可与「无增益」区分**。§4 的 “advantage” 对这一个估计量、这一个基线**没有区间支撑**。")
        else:
            p(f"· **{r.estimator}** vs {lab}：MAE 增益 {pt:+.4f}，95% 区间 [{lo:+.4f}, {hi:+.4f}] —— 不含 0。")
p("")
p("★ 口径提醒：这里的对比对象是**全样本观测中位数**这一个固定常数，不是"
  "「每个重采样里重新取中位数」的oracle常数。后者会更难超越。")
p("★ 与 §4 现有文字的关系：本脚本**不改变**任何点估计，只补一个此前缺失的区间。")
REP.close()
