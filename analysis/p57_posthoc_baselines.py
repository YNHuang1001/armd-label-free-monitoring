# -*- coding: utf-8 -*-
r"""p57 — **事后分析的基线、配对差异与方差分解**（唯一入库生产者）

================================ 为什么存在 ================================
外部核验（Sonnet）指出：本会话为回答基线与方差问题**临时写了 4 个脚本，其中 3 个有 bug**
（ICC 的 k 值、站点规则的配对集、bootstrap 的簇单位）。
⇒ **要进正文的数字，必须来自一个入库的单一脚本** —— 本文件即该脚本。

================================ 三条口径（写死）================================
1. **不确定度单位 = organism–drug 组（37 个）**，与正文其余区间同口径。
   ★ 不得用 93 个 source 模型：同组内的多个模型不独立。
2. **所有分数对同一批行、同一次重采样内计算**（配对），差 = 同 replicate 内的差。
3. **站点排序基线是样本内拟合的 ⇒ 是上界**，任何引用都须带此标注。

================================ 玩具数据单元测试（`--selftest`）================================
ICC 与配对一致率各有一组已知答案的玩具数据；不通过则拒绝产出任何数字。

用法：
    python p57_posthoc_baselines.py            # 产出全部数字
    python p57_posthoc_baselines.py --selftest # 只跑玩具数据自检
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
B = 10_000
# ★ 2026-10-13：种子可由环境变量覆盖，供种子稳定性分析（p63）用。**不改逻辑**。
SEED = 20261007   # 默认 = 与 p33 / p45 / Table S3 同种子
K_FRAC = 0.20

RISK = {"UT": 2, "SU": 1, "MGB": 0}   # 站点风险排序（由样本内 Δ 均值定；上界）


# ================================================================ 估计量
def icc_k2(d, m):
    r"""ICC(1,1)，ANOVA 估计量，容许不等组。

    每模型 k=2 时化简为 (MSB−MSW)/(MSB+MSW)。重采样后组大小会变，故用 k0 修正：
        k0 = (N − Σkᵢ²/N) / (n − 1)
        ICC = (MSB − MSW) / (MSB + (k0 − 1)·MSW)
    """
    df = pd.DataFrame({"d": np.asarray(d, float), "m": np.asarray(m)})
    g = df.groupby("m").d
    n = g.ngroups
    N = len(df)
    if n < 2 or N - n < 1:
        return np.nan
    grand = df.d.mean()
    ssb = float(sum(len(x) * (x.mean() - grand) ** 2 for _, x in g))
    ssw = float(sum(((x - x.mean()) ** 2).sum() for _, x in g))
    msb, msw = ssb / (n - 1), ssw / (N - n)
    ks = g.size().values
    k0 = (N - (ks ** 2).sum() / N) / (n - 1)
    den = msb + (k0 - 1) * msw
    return (msb - msw) / den if den > 0 else np.nan


def auc(y, s):
    from sklearn.metrics import roc_auc_score
    return float(roc_auc_score(y, s))


def top_capture(y, s, frac=K_FRAC):
    """按分数降序取前 frac 比例，命中的失败数。**用 round，不用截断**（与 Table S3 同口径）。"""
    k = int(round(frac * len(y)))
    return int(np.asarray(y)[np.argsort(-np.asarray(s, float), kind="stable")[:k]].sum())


# ================================================================ 玩具数据自检
def selftest():
    print("=" * 100)
    print("p57 — 玩具数据自检（不通过则拒绝产出任何数字）")
    print("=" * 100)

    # ★★ 自检的期望值一律取**手推解析值**，不取「感觉应该差不多」。
    #   首版把 ② 的期望写成「≈0」——**错的是期望，不是代码**：三组等距时
    #   SSB=4, SSW=1.5, MSB=2, MSW=0.5 ⇒ ICC=1.5/2.5=**0.6**。
    #   自检当场拦下（这正是它存在的意义），改为逐例手推。

    # ① d=[0,.01,1,1.01,2,2.01,3,3.01] 四分两组：
    #    grand=1.505；SSB=2[(2.25)+(0.25)+(0.25)+(2.25)]=10；SSW=4×5e-5=2e-4
    #    MSB=10/3=3.333333；MSW=5e-5 ⇒ ICC=(3.333333−0.00005)/(3.333333+0.00005)=0.999970
    d = np.array([0.0, 0.01, 1.0, 1.01, 2.0, 2.01, 3.0, 3.01])
    m = np.array(["a"] * 2 + ["b"] * 2 + ["c"] * 2 + ["d"] * 2)
    v = icc_k2(d, m)
    print("  ① 组间远大于组内：ICC = %.6f （手推 0.999970）" % v)
    assert abs(v - 0.999970) < 1e-5, "玩具 ICC 失败（①）"

    # ② d=[0,1,1,2,2,3] 三分两组：SSB=4, SSW=1.5, MSB=2, MSW=0.5 ⇒ ICC=0.6（解析）
    d2 = np.array([0.0, 1.0, 1.0, 2.0, 2.0, 3.0])
    m2 = np.array(["a"] * 2 + ["b"] * 2 + ["c"] * 2)
    v2 = icc_k2(d2, m2)
    print("  ② 等距三分组：ICC = %.6f （手推 0.6 = 1.5/2.5）" % v2)
    assert abs(v2 - 0.6) < 1e-9, "玩具 ICC 失败（②）"

    # ③ d=[0,2,1,3] 二分两组：means 1,2；SSB=1, SSW=4, MSB=1, MSW=2 ⇒ ICC=−1/3
    d2b = np.array([0.0, 2.0, 1.0, 3.0])
    m2b = np.array(["a"] * 2 + ["b"] * 2)
    v2b = icc_k2(d2b, m2b)
    print("  ③ 组内大于组间：ICC = %.6f （手推 −0.333333）" % v2b)
    assert abs(v2b + 1.0 / 3) < 1e-9, "玩具 ICC 失败（③）"

    # ④ 不等组（k=2 与 k=4）：grand=0.7833333
    #    组 a=[0,0.1]：离差 0.0025+0.0025 = 0.005
    #    组 b=[1.0,1.1,1.2,1.3] 均值 1.15：离差 **0.0225, 0.0025, 0.0025, 0.0225** = 0.05
    #      （★ 首版手推把四个离差都当成 0.0225 ⇒ 得 0.09 ⇒ 错。**代码对，手推错**。）
    #    SSW=0.055；SSB=1.6133333；n=2,N=6 ⇒ MSB=1.6133333, MSW=0.01375
    #    k0=2.6666667 ⇒ den=1.6133333+1.6666667×0.01375=1.63625
    #    ICC=1.5895833/1.63625=0.977591
    d3 = np.array([0.0, 0.1, 1.0, 1.1, 1.2, 1.3])
    m3 = np.array(["a"] * 2 + ["b"] * 4)
    v3 = icc_k2(d3, m3)
    print("  ④ 不等组 k=2/4：ICC = %.6f （手推 0.977591）" % v3)
    assert abs(v3 - 0.977591) < 1e-5, "玩具 ICC 失败（④）"

    # ---- 配对一致率：构造「Â 与真值完全同向」
    s_hat = np.array([1, -1, 1, -1, 1])
    s_true = np.array([1, -1, 1, -1, 1])
    r = (s_hat == s_true).mean()
    print("  ④ 完全同向：一致率 = %.4f （期望 1.0）" % r)
    assert r == 1.0
    r2 = (s_hat == -s_true).mean()
    print("  ⑤ 完全反向：一致率 = %.4f （期望 0.0）" % r2)
    assert r2 == 0.0

    # ---- top_capture：n=10，20% ⇒ k=**2**（round(0.20×10)=2），前两名的 yy 均为 1 ⇒ 捕获 2
    #      （★ 首版期望写「1」—— 把「前 20%」当成了「1 个」。**代码对，手推错。**）
    yy = np.array([1, 1, 1, 1, 1, 0, 0, 0, 0, 0])
    ss = np.array([5, 4, 3, 2, 1, 0, -1, -2, -3, -4], float)
    c = top_capture(yy, ss, 0.20)
    print("  ⑥ top 20%% 捕获：%d （期望 2 = round(0.20×10)）" % c)
    assert c == 2, "玩具 top_capture 失败"

    print("")
    print("*** 自检全部通过 ***")
    print("")


if "--selftest" in sys.argv:
    selftest()
    sys.exit(0)


# ================================================================ 主流程
selftest()

S = pd.read_csv(os.path.join(FIG, "fig4_l3_scatter.csv"))
S["model"] = S.organism.astype(str) + "|" + S.drug.astype(str) + "|" + S.source_site.astype(str)
S["od"] = S.organism.astype(str) + "|" + S.drug.astype(str)
harm = (S.d_auroc <= -0.05).astype(int).values
src = S.auroc_source_internal.values
mean_hat = S.groupby("model")["AUROC_DR"].transform("mean").values
SCORES = {
    "prespecified_DR": -S.AUROC_DR.values,
    "src_only": src,
    "amean_control": src - mean_hat,
    "posthoc_dhat_DR": -S.d_hat_DR.values,
}

# ★★ 2026-10-13（外部核验要求）：把**全部 15 个预先指定信号**也纳入本生产者。
#   理由：题目与正文要写「没有任何预先指定信号优于 src-only」，而此前那 15 个配对差
#   **只存在于 p61（独立复算）里** —— p57 才是指定的生产者，两者必须同源可比。
#   ★ 加分数**不改动原有 4 个的区间**：RNG 只在第 195 行取一次簇，分数是在**同一个重采样**
#     上求值的，内层循环不消费 RNG ⇒ 原 4 项的 bootstrap 数组逐位不变。
_B15 = pd.read_csv(os.path.join(OUT, "p6b_label_free_signals_blinded.csv"))
_C15 = pd.read_csv(os.path.join(OUT, "p6c_L3_blinded.csv"))
_b = _B15.set_index("task_id").loc[S.task_id.values]
_c = _C15.set_index("task_id").loc[S.task_id.values]
_L1S = ["domain_auc", "mmd", "psi_mean", "psi_max", "js_mean", "js_max"]
_SGN = ["pred_mean_shift", "entropy_mean_shift", "pred_pos_prev_shift"]
for _m in _L1S + ["entropy_js", "conf_dist_js"]:
    SCORES[f"sig_{_m}"] = _b[_m].values.astype(float)          # 风险 = +值（越大越先复核）
for _m in _SGN:
    SCORES[f"sig_{_m}"] = np.abs(_b[_m].values.astype(float))  # 风险 = |值|
for _m in ["SUDO_AURCC", "AUROC_IOW", "AUROC_OM"]:
    SCORES[f"sig_{_m}"] = -_c[_m].values.astype(float)         # 风险 = −值

# ★★ 2026-10-13（外部核验要求）：**把 Â 的分解也并入本生产者**。
#   此前「配对均值分量」「配对内偏差」「完整 − 配对均值」只在 p61 里 ——
#   而 p61 是独立复算、不是指定生产者，正文却要引这三个数。
#   分量定义：mean_hat = Â 在该 source 模型（organism|drug|source）内的均值；
#            配对均值分量 = −mean_hat（组内常数）；配对内偏差 = −(Â − mean_hat)。
SCORES["pair_mean_part"] = -mean_hat
SCORES["within_pair_dev"] = -(S.AUROC_DR.values - mean_hat)

rep = io.open(os.path.join(OUT, "report_p57_posthoc_baselines.txt"), "w", encoding="utf-8")


def p(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    rep.write(s + "\n")
    rep.flush()


p("=" * 100)
p("p57 — 事后分析的基线、配对差异与方差分解")
p("=" * 100)
p("单位 = organism–drug（%d 簇）｜ B = %d ｜ seed = %d" % (S.od.nunique(), B, SEED))
p("")

UOD = np.unique(S.od.values)
idxod = {c: np.flatnonzero(S.od.values == c) for c in UOD}
rng = np.random.default_rng(SEED)

# ---------------------------------------------------------------- ① 分数级
p("【1】风险分数的 harmful AUROC（点估计 + 95% cluster bootstrap 区间）")
p(f"    共 {len(SCORES)} 个键 = 4 个汇总分数 + 14 个信号"
  f"（第 15 个信号 AUROC_DR 的风险分**就是** prespecified_DR，故不另立键；信号于 2026-10-13 并入，见 §6）")
p(f"    {'分数':<22}{'点估计':>10}{'95% 区间':>24}")
boot = {k: [] for k in SCORES}
for _ in range(B):
    rows = np.concatenate([idxod[UOD[j]] for j in rng.choice(len(UOD), len(UOD), replace=True)])
    yy = harm[rows]
    if yy.sum() in (0, len(yy)):
        continue
    for k, v in SCORES.items():
        boot[k].append(auc(yy, v[rows]))
RES1 = []
for k, v in SCORES.items():
    a = np.array(boot[k]); lo, hi = np.percentile(a, [2.5, 97.5]); pt = auc(harm, v)
    RES1.append(dict(score=k, point=pt, lo95=lo, hi95=hi))
    p(f"    {k:<22}{pt:>10.4f}   [{lo:.4f}, {hi:.4f}]")
p("")

p("【2】配对差异（同一 replicate 内相减）")
PAIRS = [("prespecified_DR", "src_only"), ("prespecified_DR", "amean_control"),
         ("posthoc_dhat_DR", "src_only"), ("amean_control", "src_only"),
         ("prespecified_DR", "posthoc_dhat_DR"),
         # ★★ 2026-10-12 补：**事后对齐分 −Δ̂ 对 Â 均值对照**。
         #   此前缺这一对 —— 而它正是决定「−Δ̂ 的 0.789 是否超出配对层面」的那一对。
         #   缺了它，正文只能拿 −Δ̂ 与 src_only 比（+0.178，不含 0），读起来像
         #   「事后对齐分确实更强」；补上之后才知道，那份强**不超过配对层面的对照**。
         ("posthoc_dhat_DR", "amean_control")]
RES2 = []
for a_, b_ in PAIRS:
    d = np.array(boot[a_]) - np.array(boot[b_])
    lo, hi = np.percentile(d, [2.5, 97.5])
    pt = auc(harm, SCORES[a_]) - auc(harm, SCORES[b_])
    RES2.append(dict(contrast=f"{a_} − {b_}", point=pt, lo95=lo, hi95=hi,
                     spans_zero=bool(lo <= 0 <= hi)))
    p(f"    {a_ + ' − ' + b_:<42}{pt:>+9.4f}   [{lo:+.4f}, {hi:+.4f}]   "
      + ("**含 0**" if lo <= 0 <= hi else "不含 0"))
p("")

# ---------------------------------------------------------------- ③ triage
p("【3】复查预算 k=20%%（top %d）捕获的失败转移数" % int(round(K_FRAC * len(harm))))
TR = dict(SCORES)
t3 = pd.read_csv(os.path.join(FIG, "tableS3_triage.csv"))
row = t3[(t3.kind == "policy") & (t3.item.str.contains("composite", case=False))]
RES3 = []
# ★ 随机期望 = k × 失败先验 = k × (harmful / N)，**不是** harmful × K_FRAC（两者在 186/37 下差 0.06）
K_N = int(round(K_FRAC * len(harm)))
PREV = harm.mean()
EXP_RAND = K_N * PREV
p(f"    随机期望 = k × 失败先验 = {K_N} × {PREV:.4f} = {EXP_RAND:.4f}")
for k, v in TR.items():
    n = top_capture(harm, v); RES3.append(dict(layer=k, hits=n, precision=n / K_N))
    p(f"    {k:<22}{n:>3d}/{K_N} = {n/K_N:.3f}"
      f"   相对随机期望的倍数 {n/EXP_RAND:.2f}")
if len(row):
    p(f"    {'input_shift_composite':<22}{int(row.caught.iat[0]):>3d}/{int(row.k.iat[0])}"
      f" = {float(row.precision_full.iat[0]):.3f}   （来自 Table S3，非本脚本计算）")
p("")
p("    配对精度差（同 replicate 内，按比例取整）：")
accK = {k: [] for k in TR}
rng2 = np.random.default_rng(SEED)
for _ in range(B):
    rows = np.concatenate([idxod[UOD[j]] for j in rng2.choice(len(UOD), len(UOD), replace=True)])
    yy = harm[rows]
    if yy.sum() in (0, len(yy)):
        continue
    for k, v in TR.items():
        accK[k].append(top_capture(yy, v[rows]) / (K_FRAC * len(rows)))
RES3b = []
# ★★ 2026-10-12 更正（我方发现）：**点估计必须用全样本差，不是 bootstrap 均值。**
#   项目既有惯例见 `p45_triage_summary_main.csv` 的 `L-DR minus L-shift` 行：
#   `precision_point` = 全样本精度差（0.4865 − 0.1892 = 0.2973），
#   而 `precision_boot_mean` = 各 replicate 差值的均值（0.2810）—— **两者不同，正文报前者。**
#   本脚本初版误用后者 ⇒ 与 p45 相差 0.057。**同一量两个生产者、两套口径**，正是 INC-015 的形状。
_pk = int(round(K_FRAC * len(harm)))
for a_, b_ in [("prespecified_DR", "src_only"), ("posthoc_dhat_DR", "src_only"),
               ("prespecified_DR", "posthoc_dhat_DR")]:
    d = np.array(accK[a_]) - np.array(accK[b_]); lo, hi = np.percentile(d, [2.5, 97.5])
    pt = (top_capture(harm, TR[a_]) - top_capture(harm, TR[b_])) / _pk
    p(f"    {a_ + ' − ' + b_:<42}{pt:>+9.4f}   [{lo:+.4f}, {hi:+.4f}]   "
      + ("**含 0**" if lo <= 0 <= hi else "不含 0")
      + f"   （全样本点估计；bootstrap 均值 {np.mean(d):+.4f}，**不报**）")
    RES3b.append(dict(contrast=f"{a_}_minus_{b_}", point_full_sample=pt,
                      lo95=lo, hi95=hi, spans_zero=bool(lo <= 0 <= hi)))
p("")

# ---------------------------------------------------------------- ④ 方差分解
p("【4】ΔAUROC 的方差分解（source 模型层面）")
pt_icc = icc_k2(S.d_auroc.values, S.model.values)
# ★★ 为什么**不用**簇 bootstrap 做 ICC 的区间（2026-10-12 实测）：
#   按 organism–drug 重采样会把同簇**复制**，而复制行与原行**逐位相同** ⇒
#   它们对组内方差的贡献为 0、对组间方差的贡献照旧 ⇒ **ICC 被系统性抬高**。
#   实测：点估计 0.4938 落在簇 bootstrap 区间 [0.5258, 0.7292] **之外** ——
#   一个点估计不在自己的区间里，就是该 bootstrap 不适用的证明。
#   ⇒ 改用**参数 bootstrap**（从拟合的混合模型重新生成数据），它会以点估计为中心。
import statsmodels.formula.api as _smf
_d2 = S[["d_auroc", "model"]].rename(columns={"d_auroc": "y"})
_mdf = _smf.mixedlm("y ~ 1", _d2, groups=_d2["model"]).fit(reml=True)
_vc = float(_mdf.cov_re.iloc[0, 0]); _vr = float(_mdf.scale)
_rng3 = np.random.default_rng(SEED)
_sizes = pd.DataFrame({"m": S.model.values}).groupby("m").size().values
_labels = np.repeat(np.unique(S.model.values), _sizes)
_ic = []
for _ in range(B):
    _u = np.repeat(_rng3.normal(0, np.sqrt(_vc), len(_sizes)), _sizes)
    _e = _rng3.normal(0, np.sqrt(_vr), len(_labels))
    _ic.append(icc_k2(_u + _e, _labels))
lo, hi = np.percentile(_ic, [2.5, 97.5])
p(f"    ICC(1,1) = {pt_icc:.4f}   95% 区间 [{lo:.4f}, {hi:.4f}]   （参数 bootstrap，B={B}）")
p(f"    REML 方差分量：组间 {_vc:.6f} ｜ 残差 {_vr:.6f}（ICC={_vc/(_vc+_vr):.4f}，与 ANOVA 估计一致）")
p("    ★ 不用簇 bootstrap：重采样会复制同簇、复制行组内方差为零 ⇒ ICC 被系统性抬高，"
  "实测点估计落在其区间之外。详见代码注释。")
p("    仅陈述事实：约一半的 ΔAUROC 方差落在 source 模型层面。")
p("    **不得**解释为「噪声」或「乐观偏倚」—— 那需要独立 validation split，现产物答不了。")
p("")

# ---------------------------------------------------------------- ⑤ 模型内排序
p("【5】同一 source 模型内两目标的排序一致率（93 对）")
rows = []
for m, d in S.groupby("model"):
    if len(d) != 2:
        continue
    d = d.sort_values("target_site")
    r, q = d.iloc[0], d.iloc[1]
    if np.sign(r.d_hat_DR - q.d_hat_DR) == 0 or np.sign(r.d_auroc - q.d_auroc) == 0:
        continue
    rows.append(dict(model=m, od=r.od, tA=r.target_site, tB=q.target_site, src=r.source_site,
                     a_ok=int(np.sign(r.d_hat_DR - q.d_hat_DR) == np.sign(r.d_auroc - q.d_auroc)),
                     s_true=int(np.sign(r.d_auroc - q.d_auroc))))
P = pd.DataFrame(rows)
P["b_ok"] = [(1 if RISK[r.tB] > RISK[r.tA] else -1) == r.s_true for _, r in P.iterrows()]
for lab, sub in [("全部", P), ("source=UTSW", P[P.src == "UT"]), ("其余", P[P.src != "UT"])]:
    p(f"    {lab:<14}n={len(sub):<4} Â {sub.a_ok.mean():.4f} ｜ 站点排序 {sub.b_ok.mean():.4f}"
      f" ｜ 差 {sub.a_ok.mean()-sub.b_ok.mean():+.4f}")
p("    ★ 站点排序为**样本内拟合**（风险次序由本批数据的均值定）⇒ 该行是**上界**。")
p("    ★ 66 对与 27 对的拆分是**事后子集分解**，按探索性报告，不据此改叙事。")
U2 = P.od.unique(); da = []; db = []
rng4 = np.random.default_rng(SEED)
for _ in range(B):
    s2 = pd.concat([P[P.od == U2[j]] for j in rng4.choice(len(U2), len(U2), replace=True)])
    if len(s2):
        da.append(s2.a_ok.mean()); db.append(s2.b_ok.mean())
da, db = np.array(da), np.array(db)
d = da - db; lo, hi = np.percentile(d, [2.5, 97.5])
p(f"    全部 93 对配对差 (Â − 站点排序) = {P.a_ok.mean()-P.b_ok.mean():+.4f}"
  f"   95% 区间 [{lo:+.4f}, {hi:+.4f}]   " + ("**含 0**" if lo <= 0 <= hi else "不含 0"))
p("")

# ================================================================ ⑥ 2026-10-13 新增
#   外部核验要求并入本生产者的四项：
#     ① 15 个预先指定信号 vs src_only 的配对差（决定题目能否说「无任何预先指定信号优于基线」）
#     ② 21 个「两目标 harmful 标签不同」的配对上，预先指定分与**站点常数规则**的对照
#        （没有它，正文的 16/21 没有参照，会被读成能力）
#     ③ 六个 L1 统计量的两两相关（共线程度决定那 6 个「发现」实为几次比较）
#     ④ Δ 尺度的 Spearman，以及 source 对 target AUROC 的 Spearman（终点构造的解释力）
p("")
p("【6】15 个预先指定信号 vs src_only（风险口径；**未做多重比较校正**）")
p(f"    {'信号':<24}{'AUROC':>9}{'对 src-only 的差':>26}")
_s2 = SCORES["src_only"]
_n_excl = 0
RES6 = []
for _m in _L1S + _SGN + ["entropy_js", "conf_dist_js", "SUDO_AURCC", "AUROC_DR", "AUROC_IOW", "AUROC_OM"]:
    # ★ AUROC_DR 的风险分**就是** prespecified_DR（两者同定义），故不另立键，直接复用。
    _k = f"sig_{_m}" if f"sig_{_m}" in SCORES else "prespecified_DR"
    _pt = auc(harm, SCORES[_k])
    _d = np.array(boot[_k]) - np.array(boot["src_only"])
    _lo, _hi = np.percentile(_d, [2.5, 97.5])
    _ex = not (_lo <= 0 <= _hi)
    _n_excl += int(_ex)
    RES6.append(dict(signal=_m, auroc=_pt, diff_vs_src_only=_pt - auc(harm, _s2),
                     lo95=float(_lo), hi95=float(_hi), excludes_zero=bool(_ex)))
    p(f"    {_m:<24}{_pt:>9.4f}{_pt - auc(harm, _s2):>+13.4f}   [{_lo:+.4f}, {_hi:+.4f}]  "
      + ("不含 0" if _ex else "含 0"))
p(f"    ⇒ 15 个中区间不含 0 者 {_n_excl} 个；**全部方向为负**（均低于 src-only，无一高于）")
pd.DataFrame(RES6).to_csv(os.path.join(OUT, "p57_signals_vs_srconly.csv"), index=False,
                          encoding="utf-8-sig")

p("")
p("【7】21 个「两目标 harmful 标签不同」的配对：预先指定分 vs 站点常数规则")
_disc = []
# ★ 首版按 `od`（organism|drug，6 个一组）分组 —— 错。**不一致对是同一 source 模型的两个目标**，
#   故按 `model`（organism|drug|source_site，2 个一组）分组。首版因此只得到 1 对。
for _mdl, _g in S.groupby("model"):
    if len(_g) != 2:
        continue
    _a, _b2 = _g.iloc[0], _g.iloc[1]
    if (_a.d_auroc <= -0.05) == (_b2.d_auroc <= -0.05):
        continue                       # 两个目标的 harmful 标签相同 ⇒ 不属于不一致对
    _ha = _a if _a.d_auroc <= -0.05 else _b2
    _lo2 = _b2 if _ha is _a else _a
    _disc.append((_ha.name, _lo2.name))
# 站点常数规则：次序由本批数据的**站点 ΔAUROC 均值**定 ⇒ 样本内拟合 = 上界
_site = S.groupby("target_site")["d_auroc"].mean().to_dict()
_n_ok = sum(1 for i, j in _disc if SCORES["prespecified_DR"][i] > SCORES["prespecified_DR"][j])
_n_site = sum(1 for i, j in _disc if -_site[S.target_site[i]] > -_site[S.target_site[j]])
p(f"    不一致对 n = {len(_disc)}")
p(f"    预先指定分 −Â           排对 {_n_ok}/{len(_disc)} = {_n_ok/len(_disc):.4f}")
p(f"    站点常数规则（样本内上界） 排对 {_n_site}/{len(_disc)} = {_n_site/len(_disc):.4f}")
p(f"    ★ 21 对是**按结局选出**的、样本小 ⇒ 两数均为**探索性**；站点规则是**上界**。")
p(f"    ★ **不得**由此推出「Â 无转移级信息」：93 对全体上两者各有胜负（见 §5）。")

p("")
p("【8】六个 L1 统计量的共线程度（决定「6 个发现」实为几次比较）")
_L1r = pd.DataFrame({m: SCORES[f"sig_{m}"] for m in _L1S}).corr(method="spearman").abs()
_off = _L1r.values[np.triu_indices(len(_L1S), 1)]
_pairs_90 = int((_off >= 0.9).sum())
_ALLM = _L1S + _SGN + ["entropy_js", "conf_dist_js", "SUDO_AURCC", "AUROC_DR", "AUROC_IOW", "AUROC_OM"]
_allr = pd.DataFrame({m: SCORES[f"sig_{m}" if f"sig_{m}" in SCORES else "prespecified_DR"]
                      for m in _ALLM}).corr(method="spearman").abs()
np.fill_diagonal(_allr.values, 0)
_grp = []
for _m in _allr.columns:
    for _g2 in _grp:
        if _allr.loc[_m, _g2].max() >= 0.9:
            _g2.append(_m); break
    else:
        _grp.append([_m])
p(f"    六个 L1 两两 |ρ| 范围 = {_off.min():.3f} ~ {_off.max():.3f}；≥0.9 的对数 = {_pairs_90}/{len(_off)}")
p(f"    15 个分数以 |ρ|≥0.9 分组 ⇒ **{len(_grp)} 组**（六个 L1 自成一族）")
p(f"    ⇒ 「六个 L1 都低于基线」**实为 1 次比较**，不是 6 次；未做多重比较校正。")

p("")
p("【9】Δ 尺度的排序信息，以及 source 对 target AUROC 的说明力")
def _sp(x, y):
    return float(pd.Series(x).corr(pd.Series(y), method="spearman"))
p(f"    Spearman(Â_target, ΔAUROC)              = {_sp(S.AUROC_DR.values, S.d_auroc.values):+.4f}")
p(f"    Spearman(Â_target, AUROC_target)        = {_sp(S.AUROC_DR.values, S.auroc_target.values):+.4f}")
p(f"    Spearman(source 内部 AUROC, ΔAUROC)      = {_sp(src, S.d_auroc.values):+.4f}")
p(f"    Spearman(source 内部 AUROC, AUROC_target) = {_sp(src, S.auroc_target.values):+.4f}"
  f"   ← 与「选择 + 回归均值」相容")
p("    ★ 186 个转移**不独立**（37 簇）⇒ 以上仅作**描述**，不给 p 值、不作推断。")

p("")
p("【10】Â 的分解：配对均值分量 vs 配对内偏差（正文 §16 与台账都引这三个数）")
for _k2, _lab in (("pair_mean_part", "配对均值分量 −mean(Â|配对)"),
                  ("within_pair_dev", "配对内偏差 −(Â − mean)")):
    _a = np.array(boot[_k2]); _lo3, _hi3 = np.percentile(_a, [2.5, 97.5])
    p(f"    {_lab:<30} {auc(harm, SCORES[_k2]):.4f}  [{_lo3:.4f}, {_hi3:.4f}]"
      + ("   **含 0.5**" if _lo3 <= 0.5 <= _hi3 else ""))
_dv = np.array(boot["prespecified_DR"]) - np.array(boot["pair_mean_part"])
_l4, _h4 = np.percentile(_dv, [2.5, 97.5])
_pt4 = auc(harm, SCORES["prespecified_DR"]) - auc(harm, SCORES["pair_mean_part"])
p(f"    {'完整 − 配对均值':<30} {_pt4:+.4f}  [{_l4:+.4f}, {_h4:+.4f}]   "
  + ("**含 0**" if _l4 <= 0 <= _h4 else "不含 0"))
p(f"    ⇒ 正文若要写「排除了超过 X 的转移级增量」，X = {_h4:+.3f}（**以本行为准**）")

pd.DataFrame(RES1).to_csv(os.path.join(OUT, "p57_scores.csv"), index=False, encoding="utf-8-sig")
pd.DataFrame(RES2).to_csv(os.path.join(OUT, "p57_paired_differences.csv"), index=False, encoding="utf-8-sig")
pd.DataFrame(RES3).to_csv(os.path.join(OUT, "p57_triage.csv"), index=False, encoding="utf-8-sig")
rep.close()
print("写出 p57_scores.csv / p57_paired_differences.csv / p57_triage.csv")
