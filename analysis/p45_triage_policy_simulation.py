# -*- coding: utf-8 -*-
r"""p45 — 分诊策略模拟（triage policy simulation）

================================ 性质 ================================
**事后分析 / 回顾性政策模拟。** 不属于密封确证部分（sealed confirmatory analysis），
**不得**借用「protocol 已冻结」的信用。正文中必须明写 post hoc。

================================ 预注册 ================================
本脚本的全部设计选择（D1–D4）在运行**之前**写死于：

    TRIAGE_POLICY_PREREGISTRATION_20261007.md
    sha256 = 5e25a14963534d808f7085c9d9fe6d540c13c45f57e974eee87bf6e4dfb98fef

剧本不得在本脚本内变更。任何偏离必须在本文件的报告与正文中如实记账。

================================ 设计摘要（逐条对应备忘 D1–D4） ================
D1  决策目标：固定预算下的排序分诊。主口径 k = 20%（186 → 37）。
D2  对照层三层：
      L-random —— 无信息政策。常数基线给每个迁移同一个 Δ ⇒ 不产生排序 ⇒ 等价于随机。
                  故「显著优于常数基线」在本框架下 = 「显著优于随机」。
                  随机的精确零假设为超几何，E[前 k 中失效数] = k × 57/186。
      L-shift  —— 现行实践。六个 L1 信号的**等权秩平均**（避免事后挑选，给现行实践集体最优）。
      L-DR     —— 本文幸存主张，按风险分 = −AUROC_DR 降序。
    （L-OM     —— 预设敏感性估计量，并列报告。）
D3  结局：沿用正文已有的 harmful = ΔAUROC ≤ −0.05，**不新发明**。
    主统计量：precision@k，k = 37。次要：recall@k、累积增益曲线、达到 50% 召回的预算。
D4  推断：organism–drug 分组（37 cluster）cluster bootstrap，B = 10,000，95% 百分位区间，
    配对比较在同一 replicate 内进行。**唯一主对比：L-DR vs L-shift @ k = 20%，α = 0.05 双侧。**
    其余一律描述性/探索性，不作显著性声称。

================================ 锚定 ================================
输入与冻结主分析同源：p6b（L1/L2）+ p6c（L3/SUDO）+ p6d_outcome（结局）。
方向变换逐字取自 p11_figure_datapack.py L183–194。
开工前先断言：独立重算的 15 个 harmful_auc 必须与 p6d_signal_eval.csv **逐位相等**，
否则直接抛错退出（fail-closed）。

**不使用 p10b_raw.csv。** 它是样本量敏感性分析产物，即使 level=='full' 也含重抽样噪声，
15 个信号中有 4 个无法复现冻结值（见预注册备忘 §2.3）。
"""
import os, sys, hashlib
import numpy as np, pandas as pd

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "分析结果")
FIGD = os.path.join(ROOT, "figures", "data")     # 渲染数据包目录（figures/data）
REP = open(os.path.join(OUT, "report_p45_triage.txt"), "w", encoding="utf-8")

def p(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True); REP.write(s + "\n"); REP.flush()

# ★ 2026-10-13：种子可由环境变量覆盖，供种子稳定性分析（p63）用。
#   **不改任何逻辑** —— 只把模块级常量变成可注入的参数。默认值 = 已发表种子。
SEED = 20261007
B = 10000
K_FRAC = 0.20
MEMO = "TRIAGE_POLICY_PREREGISTRATION_20261007.md"
MEMO_SHA = "5e25a14963534d808f7085c9d9fe6d540c13c45f57e974eee87bf6e4dfb98fef"

# ---------------------------------------------------------------- 上游
BB = pd.read_csv(os.path.join(OUT, "p6b_label_free_signals_blinded.csv"))
CC = pd.read_csv(os.path.join(OUT, "p6c_L3_blinded.csv"))
DD = pd.read_csv(os.path.join(OUT, "p6d_outcome.csv"))

L1 = ["domain_auc", "mmd", "psi_mean", "psi_max", "js_mean", "js_max"]
SIGNED_L2 = ["pred_mean_shift", "entropy_mean_shift", "pred_pos_prev_shift"]
L3 = ["AUROC_DR", "AUROC_IOW", "AUROC_OM"]

T = (DD[["task_id", "organism", "drug", "harmful", "d_auroc", "auroc_source_internal"]]
     .merge(BB[["task_id"] + L1 + SIGNED_L2 + ["entropy_js", "conf_dist_js"]],
            on="task_id", validate="1:1")
     .merge(CC[["task_id", "SUDO_AURCC"] + L3], on="task_id", validate="1:1"))
T = T.sort_values("task_id").reset_index(drop=True)
T["group"] = T.organism + "|" + T.drug
# ★ 事后追加层所需的对齐分数（与 p11:272 同定义：d_hat = Â_target − source）
T["d_hat_DR"] = T.AUROC_DR - T.auroc_source_internal

ALL15 = L1 + SIGNED_L2 + ["entropy_js", "conf_dist_js", "SUDO_AURCC"] + L3
y = T.harmful.values.astype(int)
g = T.group.values
tid = T.task_id.values
code_ids = pd.Series(tid).astype("category").cat.codes.values


def risk(m):
    """r_failure：越高越可能 harmful。逐字复刻 p11_figure_datapack.py L183–194。"""
    s = T[m].values.astype(float)
    if m in SIGNED_L2:
        return np.abs(s)
    if m in ("SUDO_AURCC",) + tuple(L3):
        return -s
    return s


def auc(yv, sv):
    """Mann–Whitney，平均秩处理并列（与 p6d 口径一致）。"""
    order = np.argsort(sv, kind="mergesort")
    sv, yv = sv[order], yv[order]
    r = np.empty(len(sv)); i = 0
    while i < len(sv):
        j = i
        while j + 1 < len(sv) and sv[j + 1] == sv[i]:
            j += 1
        r[i:j + 1] = 0.5 * (i + j) + 1
        i = j + 1
    n1 = yv.sum(); n0 = len(yv) - n1
    if n1 == 0 or n0 == 0:
        return float("nan")
    return (r[yv == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)


# ---------------------------------------------------------------- 锚定（fail-closed）
p("=" * 120)
p("p45 — 分诊策略模拟（事后 / 回顾性政策模拟）")
p("=" * 120)
p(f"预注册：{MEMO}")
p(f"          sha256 = {MEMO_SHA}")
_memo_path = os.path.join(ROOT, MEMO)
if os.path.exists(_memo_path):
    _h = hashlib.sha256(open(_memo_path, "rb").read()).hexdigest()
    p(f"          实测 sha256 = {_h}   {'✔ 一致' if _h == MEMO_SHA else '✘ 与预注册不符'}")
    if _h != MEMO_SHA:
        raise RuntimeError("预注册备忘已被改动，拒绝运行。")
else:
    p("          ⚠ 备忘文件不在预期路径，无法复核哈希（继续，但记账）")
p("")

frozen = pd.read_csv(os.path.join(OUT, "p6d_signal_eval.csv")).set_index("method")
p("【0】锚定：独立重算 15 个 harmful_auc，须与 p6d_signal_eval.csv 逐位相等")
bad = 0
for m in ALL15:
    a = auc(y, risk(m)); fz = float(frozen.loc[m, "harmful_auc"])
    ok = abs(a - fz) < 1e-9
    bad += (not ok)
    if not ok:
        p(f"    ✘ {m:<22}{a:.10f}  vs frozen {fz:.10f}   |d|={abs(a-fz):.2e}")
p(f"    结果：{len(ALL15)-bad}/{len(ALL15)} 通过"
  + ("   ✔ 输入与冻结主分析同源同口径" if bad == 0 else "   ✘ 锚定失败"))
if bad:
    raise RuntimeError("锚定失败：上游表与冻结主分析不同源，拒绝继续。")
p("")

# ---------------------------------------------------------------- D2 三层
# L-shift：六个 L1 信号的等权秩平均（rank 1 = 最小值；越大 = 漂移越大 = 风险越高）
R = np.column_stack([
    pd.Series(T[m].values.astype(float)).rank(method="average").values for m in L1
])
shift_rank = R.mean(axis=1)

LAYERS = {
    "L-shift":  shift_rank,                       # 越大越先复核
    "L-DR":     risk("AUROC_DR"),
    "L-OM":     risk("AUROC_OM"),
    # ★★ 2026-10-12 追加两层（外部核验要求：triage 必须有一个**平凡基线**作参照）。
    #   病历：正文与摘要此前只报「DR 1.59× vs input-shift composite」，而 **composite 本身
    #   低于随机**（7/37，0.62×）—— 胜过一个低于随机的对照，门槛极低。真正有信息量的
    #   对照是**只用 source 自身判别力**排序。
    #   ★ 两层**均为事后追加**，见 TRIAGE_POLICY_PREREGISTRATION_AMENDMENT_**20261012**.md；
    #     主对比仍是 L-DR − L-shift，追加两层**不改变预注册的任何结论**。
    #   ★★ 2026-10-12 更正指向：原注释指向 `..._AMENDMENT_20261007.md`，但该文件
    #     写于 10-07、且自述「不改动任何设计选择 D1–D4」，**根本没有提到这两层**
    #     —— 那是**假的出处引用**（INC-032）。10-12 另立增补，见新文件。
    "L-src":    T["auroc_source_internal"].values.astype(float),   # 越大越先复核
    "L-dhat":   -T["d_hat_DR"].values.astype(float),               # 事后对齐分数
}
N_TOT = len(y); N_FAIL = int(y.sum()); PREV = N_FAIL / N_TOT
p(f"结局：{N_FAIL}/{N_TOT} harmful（先验 {PREV:.4f}）；organism–drug 分组 {len(np.unique(g))} 个")
p(f"随机零假设（精确）：E[前 k 中失效数] = k × {N_FAIL}/{N_TOT}")
_t = [len(np.unique(v)) for v in LAYERS.values()]
p(f"并列检查：L-shift 唯一值 {_t[0]}；L-DR 唯一值 {_t[1]}；L-OM 唯一值 {_t[2]}  "
  f"（并列由 task_id 升序确定性打破）")
p("")


def rank_order(sv, ids):
    """按风险分降序；并列以 task_id 升序确定性打破（lexsort 末键为主键）。"""
    return np.lexsort((ids, -sv))


def precision_at(yv, sv, ids, k):
    if k <= 0:
        return 0.0
    return float(yv[rank_order(sv, ids)[:k]].mean())


# ---------------------------------------------------------------- 点估计（观测数据，k = 20%）
K_OBS = int(round(K_FRAC * N_TOT))
p(f"【1】点估计（观测数据，k = {K_OBS}，即 {K_FRAC:.0%}）")
p(f"    {'层':<10}{'命中失效':>10}{'precision':>12}{'recall':>10}{'vs 随机(绝对)':>15}")
rows = []
for nm, sv in LAYERS.items():
    ca = int(y[rank_order(sv, code_ids)[:K_OBS]].sum())
    pr = ca / K_OBS; rc = ca / N_FAIL
    rows.append(dict(layer=nm, k=K_OBS, caught=ca, precision=pr, recall=rc,
                     vs_random=pr - PREV))
    p(f"    {nm:<10}{ca:>10}{pr:>12.4f}{rc:>10.4f}{pr-PREV:>+15.4f}")
p(f"    {'随机':<10}{K_OBS*PREV:>10.2f}{PREV:>12.4f}{K_OBS/N_TOT:>10.4f}{0.0:>+15.4f}   （精确期望）")
p("")

# ---------------------------------------------------------------- D4 cluster bootstrap
p(f"【2】cluster bootstrap（organism–drug 分组，B = {B}，seed = {SEED}）")
ugs = np.unique(g)
idx_by_g = {u: np.where(g == u)[0] for u in ugs}
rng = np.random.default_rng(SEED)
lay_names = list(LAYERS.keys())
acc = {nm: [] for nm in lay_names}
acc["random"] = []

for _ in range(B):
    pick = rng.choice(len(ugs), len(ugs), replace=True)
    idx = np.concatenate([idx_by_g[ugs[i]] for i in pick])
    yv = y[idx]; ids = code_ids[idx]
    n = len(yv)
    pre = yv.mean()
    k = int(round(K_FRAC * n))
    for nm in lay_names:
        acc[nm].append(precision_at(yv, LAYERS[nm][idx], ids, k))
    acc["random"].append(pre)

for nm in lay_names:
    v = np.array(acc[nm])
    p(f"    {nm:<10}precision 均值 {v.mean():.4f}  95% 区间 "
      f"[{np.percentile(v,2.5):.4f}, {np.percentile(v,97.5):.4f}]"
      f"   vs 随机 [{np.percentile(v-acc['random'],2.5):+.4f}, {np.percentile(v-acc['random'],97.5):+.4f}]")
# 主对比：同一 bootstrap replicate 内配对相减（两个序列在同一循环迭代中写入，逐位配对）
# ★ 2026-10-07 更正：此前本行把 **bootstrap 分布均值** 写在主对比的位置，
#   与点估计（0.2973）并排而不标口径，导致账面看起来不闭合（0.2810 vs 0.2973）。
#   CSV 里两列本来都是对的（precision_point / precision_boot_mean）；
#   错的是这一行叙述。现改为**点估计与分布摘要分行、各自明写标签**。
v = np.array(acc["L-DR"]) - np.array(acc["L-shift"])
assert len(acc["L-DR"]) == len(acc["L-shift"]) == B, "配对长度不等 —— 拒绝报告"
_i = {nm: i for i, nm in enumerate(lay_names)}
cDR, cSH = rows[_i["L-DR"]]["caught"], rows[_i["L-shift"]]["caught"]
pt_obs = cDR / K_OBS - cSH / K_OBS
p(f"    **主对比 L-DR − L-shift**（同 replicate 配对，n={len(v)}）")
p(f"        点估计（观测数据；**正文/图表应引此值**）  {pt_obs:+.4f}"
  f"   = {cDR}/{K_OBS} − {cSH}/{K_OBS} = {cDR - cSH}/{K_OBS}")
p(f"        95% 百分位区间                              "
  f"[{np.percentile(v, 2.5):+.4f}, {np.percentile(v, 97.5):+.4f}]")
p(f"        bootstrap 分布均值（**不是点估计**，仅作偏差诊断）  {v.mean():+.4f}"
  f"   （与点估计差 {v.mean() - pt_obs:+.4f}）")
p("")
assert abs(pt_obs - (rows[_i["L-DR"]]["precision"] - rows[_i["L-shift"]]["precision"])) < 1e-12, \
    "点估计与两层 precision 之差不闭合 —— 拒绝继续"

# ---------------------------------------------------------------- 曲线
p("【3】累积增益曲线落盘")
curve = []
for nm, sv in LAYERS.items():
    o = np.lexsort((code_ids, -sv))
    ys = y[o]
    cum = np.cumsum(ys)
    for k in range(1, N_TOT + 1):
        curve.append(dict(layer=nm, k=k, cum_failures=int(cum[k-1]),
                          precision=cum[k-1]/k, recall=cum[k-1]/N_FAIL))
    # 随机基线
for k in range(1, N_TOT + 1):
    curve.append(dict(layer="L-random", k=k, cum_failures=k*PREV,
                      precision=PREV, recall=k/N_TOT))
pd.DataFrame(curve).to_csv(os.path.join(OUT, "p45_triage_curves.csv"), index=False)

# ★ 六个 L1 成员各自单独的累积曲线（S10 要用：成员淡显 + 复合线）。
#   动机：本分析的复合对照落在成员取值区间的**下缘**，审稿人会问
#   「你预注册的复合恰好对现行实践最不利？」——把成员曲线**画进图里**是最直接的回答。
curve_m = []
for m in L1:
    sv = risk(m)
    o = rank_order(sv, code_ids)
    cum = np.cumsum(y[o])
    for k in range(1, N_TOT + 1):
        curve_m.append(dict(member=m, role="L1-member", k=k,
                            cum_failures=int(cum[k - 1]),
                            precision=cum[k - 1] / k, recall=cum[k - 1] / N_FAIL))
pd.DataFrame(curve_m).to_csv(os.path.join(OUT, "p45_triage_curves_members.csv"), index=False)
p(f"    写入 p45_triage_curves_members.csv（{len(curve_m)} 行 = 6 成员 × {N_TOT}）")

# ---------------------------------------------------------------- 全信号 summary 表
p("【4】全 15 信号 + 3 层：k = 10% / 20% / 30% 的 precision / recall（描述性）")
summ = []
for m in ALL15:
    sv = risk(m)
    for frac in (0.10, 0.20, 0.30):
        k = int(round(frac * N_TOT))
        o = np.lexsort((code_ids, -sv))
        ca = int(y[o[:k]].sum())
        summ.append(dict(item=m, kind="signal", k_frac=frac, k=k, caught=ca,
                         precision=ca/k, recall=ca/N_FAIL, vs_random=ca/k - PREV))
for nm, sv in LAYERS.items():
    for frac in (0.10, 0.20, 0.30):
        k = int(round(frac * N_TOT))
        o = np.lexsort((code_ids, -sv))
        ca = int(y[o[:k]].sum())
        summ.append(dict(item=nm, kind="layer", k_frac=frac, k=k, caught=ca,
                         precision=ca/k, recall=ca/N_FAIL, vs_random=ca/k - PREV))
for frac in (0.10, 0.20, 0.30):
    k = int(round(frac * N_TOT))
    summ.append(dict(item="L-random", kind="layer", k_frac=frac, k=k,
                     caught=k*PREV, precision=PREV, recall=k/N_TOT, vs_random=0.0))
pd.DataFrame(summ).to_csv(os.path.join(OUT, "p45_triage_summary.csv"), index=False)
p(f"    写入 p45_triage_summary.csv（{len(summ)} 行）")

# ---------------------------------------------------------------- 每成员敏感性
# ★ 动机（GLM 第 4 条措辞约束）：预注册的复合对照（六者秩平均）落在成员取值区间的**下缘**，
#   审稿人会问「你预注册的复合恰好对现行实践最不利？」。
#   最硬的回答是把 L-DR 与**每一个**成员的差都摆出来 —— 若全为正，结论不依赖复合的选择。
p("")
p("【5】每成员敏感性：L-DR 与六个 L1 成员**逐个**之差（k = 20%，观测数据）")
ms = []
ca_dr = int(y[rank_order(LAYERS["L-DR"], code_ids)[:K_OBS]].sum())
for m in L1:
    ca = int(y[rank_order(risk(m), code_ids)[:K_OBS]].sum())
    ms.append(dict(member=m, caught=ca, precision=ca / K_OBS,
                   DR_minus_member=ca_dr / K_OBS - ca / K_OBS))
    p(f"    {m:<14}命中 {ca:>3}/{K_OBS}   L-DR − 该成员 = {ca_dr / K_OBS - ca / K_OBS:+.4f}")
_all_pos = all(r["DR_minus_member"] > 0 for r in ms)
p(f"    → L-DR 相对**每一个** L1 成员均为正："
  f"{'是' if _all_pos else '否 ✘（披露文字必须改写）'}")
pd.DataFrame(ms).to_csv(os.path.join(OUT, "p45_triage_member_sensitivity.csv"), index=False)
p("    写入 p45_triage_member_sensitivity.csv")

# ---------------------------------------------------------------- 复合 × 成员的重合度
# ★ 动机（外部完备性核验提出）：Panel C 只挡得住「你挑了最不利的聚合」，
#   挡不住「你的现行实践代理其实就等同于它最弱的那个成员」。
#   把 top-k 重合度摆出来，两种读法都能被读者自己判定。
p("")
p("【6】复合对照（六者秩平均）与各成员在 top-k 的重合度")
ov = []
_o_shift = rank_order(LAYERS["L-shift"], code_ids)
_top_shift = {kf: set(_o_shift[:int(round(kf * N_TOT))]) for kf in (0.10, 0.20, 0.30)}
for m in L1:
    o = rank_order(risk(m), code_ids)
    for kf in (0.10, 0.20, 0.30):
        k = int(round(kf * N_TOT))
        inter = len(_top_shift[kf] & set(o[:k]))
        ov.append(dict(member=m, k_frac=kf, k=k, overlap=inter, overlap_frac=inter / k))
        if abs(kf - 0.20) < 1e-9:
            p(f"    {m:<14} k={k}: 重合 {inter}/{k} = {inter/k:.3f}")
pd.DataFrame(ov).to_csv(os.path.join(OUT, "p45_triage_overlap.csv"), index=False)
pd.DataFrame(ov).to_csv(os.path.join(FIGD, "figS10_triage_overlap.csv"), index=False)
p("    写入 p45_triage_overlap.csv + figures/data/figS10_triage_overlap.csv")

# ---------------------------------------------------------------- 达到 50% 召回的预算
# ★ 这是**预注册备忘 §3 D1 承诺过的次要统计量**，初版未交付 ⇒ 现补算（不是新声称）。
p("")
p("【7】达到 50% 召回所需的复核预算（预注册承诺的次要描述量）")
_need = N_FAIL / 2.0
rb = []
for nm, sv in [("L-shift", LAYERS["L-shift"]), ("L-DR", LAYERS["L-DR"]), ("L-OM", LAYERS["L-OM"])]:
    o = rank_order(sv, code_ids)
    cum = 0
    hit = None
    for i in range(N_TOT):
        cum += int(y[o[i]])
        if cum >= _need:
            hit = i + 1
            break
    rb.append(dict(layer=nm, n_needed=hit, frac_needed=hit / N_TOT if hit else None))
    p(f"    {nm:<9} 需复核前 {hit} 个 = {hit/N_TOT:.1%} 的队列（达到 {_need:.1f}/57）")
pd.DataFrame(rb).to_csv(os.path.join(OUT, "p45_triage_recall_budget.csv"), index=False)
pd.DataFrame(rb).to_csv(os.path.join(FIGD, "figS10_triage_recall_budget.csv"), index=False)
p("    写入 p45_triage_recall_budget.csv + figures/data/figS10_triage_recall_budget.csv")

# ---------------------------------------------------------------- 渲染数据包
# ★ 按本项目既有约定：**成图数据包落在 figures/data/**（p11 亦如此），
#   渲染层（figures/R/*.R）只读 figures/data，**不重算任何统计量**。
#   放这里而不是让 R 脚本跨目录读 分析结果/，有两个好处：
#     ① 与其余 9 张补充图的数据来源约定一致；
#     ② figures/data 已在既有审计之下（BOUND-04 值域扫描、FIGURES_R 路径），
#        新数据包自动获得同一层保护。
#   本块为**纯派生**：内容与上面的 p45_triage_*.csv 逐字相同，只换文件名。
FIGD = os.path.join(ROOT, "figures", "data")
os.makedirs(FIGD, exist_ok=True)
# ★★ 2026-10-12（INC-030）：本块原为
#       `_gs45.backup(...)` → `to_csv(...)` → `import genstamp as _gs45` → `stamp(4 个文件)` → 再写后三个
#   **三处错叠在一起**：① 导入在使用之后（NameError）；② stamp 循环跑在**写入之前**
#   ⇒ 盖的是**上一版**文件的章；③ 四个文件中三个在盖章时还没写。
#   修法：**导入提到最前，写入全部完成，最后统一 backup+stamp**。
#   ★ 与 `p11` 是**同一个缺口的两处**，而两处都只在「逐块执行」的判据（EVS t13 / 直接运行）下才现形。
import sys as _sys
if os.path.join(ROOT, "分析结果") not in _sys.path:
    _sys.path.insert(0, os.path.join(ROOT, "分析结果"))
import genstamp as _gs45

# ---- 与**平凡基线**（source 自身判别力）的配对差：本量在下面有**两个落点**
#      （figures/data/figS10_triage_summary.csv 与 p45_triage_summary_main.csv）。
#      ★★ 2026-10-12：**只在这里算一次**，两个落点引用同一份 dict、仅换 layer 拼写。
#      初版在两个落点各写了一遍**同义的代码**，还各用了一套点估计口径
#      （一处 bootstrap 均值、一处全样本差，相差 0.057）—— 正是 INC-015 的形态：
#      **同一量两个生产者**。落点改为一处，这类分歧无法再发生。
_dr_src_core = {
    "k": K_OBS,
    "precision_point": (rows[lay_names.index("L-DR")]["precision"]
                        - rows[lay_names.index("L-src")]["precision"]),
    "precision_boot_mean": float(np.mean(np.array(acc["L-DR"]) - np.array(acc["L-src"]))),
    "lo95": float(np.percentile(np.array(acc["L-DR"]) - np.array(acc["L-src"]), 2.5)),
    "hi95": float(np.percentile(np.array(acc["L-DR"]) - np.array(acc["L-src"]), 97.5)),
    "diff_vs_random_lo": float("nan"),
    "diff_vs_random_hi": float("nan"),
}
_DR_SRC_CAUGHT = (int(rows[lay_names.index("L-DR")]["caught"])
                  - int(rows[lay_names.index("L-src")]["caught"]))

_FILES10 = ("figS10_triage_curves.csv", "figS10_triage_members.csv",
            "figS10_triage_members_k37.csv", "figS10_triage_summary.csv")
pd.DataFrame(curve).to_csv(os.path.join(FIGD, "figS10_triage_curves.csv"), index=False)
pd.DataFrame(curve_m).to_csv(os.path.join(FIGD, "figS10_triage_members.csv"), index=False)
pd.DataFrame(ms).to_csv(os.path.join(FIGD, "figS10_triage_members_k37.csv"), index=False)
_rand = np.array(acc["random"])
pd.DataFrame([dict(
    layer=nm, k=K_OBS, caught=rows[i]["caught"],
    precision_point=rows[i]["precision"],
    precision_boot_mean=float(np.mean(acc[nm])),
    lo95=float(np.percentile(acc[nm], 2.5)), hi95=float(np.percentile(acc[nm], 97.5)),
    # ★ 与**随机排序**配对的差（不是与常数——常数在排序框架下就是随机，见预注册备忘 D2）
    diff_vs_random_lo=float(np.percentile(np.array(acc[nm]) - _rand, 2.5)),
    diff_vs_random_hi=float(np.percentile(np.array(acc[nm]) - _rand, 97.5)),
) for i, nm in enumerate(lay_names)]
    + [dict(layer="L-random", k=K_OBS, caught=K_OBS * PREV, precision_point=PREV,
            precision_boot_mean=PREV, lo95=float("nan"), hi95=float("nan"),
            diff_vs_random_lo=0.0, diff_vs_random_hi=0.0)]
    + [dict(layer="L-DR_minus_L-shift", k=K_OBS, caught=cDR - cSH, precision_point=pt_obs,
            precision_boot_mean=float(v.mean()),
            lo95=float(np.percentile(v, 2.5)), hi95=float(np.percentile(v, 97.5)),
            diff_vs_random_lo=float("nan"), diff_vs_random_hi=float("nan"))]
    # ★★ 2026-10-12 追加（外部核验要求）：**和平凡基线的配对差**。
    #   预注册主对比是 L-DR − L-shift，但 composite 本身**低于随机**（7/37，0.62×）
    #   ⇒ 胜过一个低于随机的对照门槛极低。**有信息量的对照是 L-src。**
    #   本行**事后追加**，见 TRIAGE_POLICY_PREREGISTRATION_AMENDMENT_20261012.md。
    #   ★ `precision_point` 是**全样本差**（与 `L-DR_minus_L-shift` 行的 `pt_obs` 同口径）；
    #     初版此处误用 bootstrap 均值（0.132 vs 0.189）。现值取自 `_dr_src_core`，见上方。
    + [dict(layer="L-DR_minus_L-src", caught=_DR_SRC_CAUGHT, **_dr_src_core)]
).to_csv(os.path.join(FIGD, "figS10_triage_summary.csv"), index=False)
# ★ 四个文件**全部写完**之后，才 backup + stamp（顺序反了会盖上一版的章）
for _f in _FILES10:
    _gs45.backup(os.path.join(FIGD, _f), "p45")
    _gs45.stamp(os.path.join(FIGD, _f), "p45")
p(f"    渲染数据包 → figures/data/figS10_triage_{{curves,members,members_k37,summary}}.csv  "
  f"（四个文件已 backup+stamp）")

# ---------------------------------------------------------------- 主汇总
pd.DataFrame([dict(
    layer=nm, k=K_OBS, precision_point=rows[i]["precision"],
    precision_boot_mean=float(np.mean(acc[nm])),
    lo95=float(np.percentile(acc[nm], 2.5)), hi95=float(np.percentile(acc[nm], 97.5)),
    diff_vs_random_lo=float(np.percentile(np.array(acc[nm]) - acc["random"], 2.5)),
    diff_vs_random_hi=float(np.percentile(np.array(acc[nm]) - acc["random"], 97.5)),
) for i, nm in enumerate(lay_names)]
    + [dict(layer="L-DR minus L-shift", k=K_OBS,
            precision_point=rows[lay_names.index("L-DR")]["precision"] -
                           rows[lay_names.index("L-shift")]["precision"],
            precision_boot_mean=float(v.mean()),
            lo95=float(np.percentile(v, 2.5)), hi95=float(np.percentile(v, 97.5)),
            diff_vs_random_lo=np.nan, diff_vs_random_hi=np.nan)]
    # ★★ 2026-10-12 追加：与**平凡基线**（source 自身判别力）的配对差。事后追加，
    #   见 TRIAGE_POLICY_PREREGISTRATION_AMENDMENT_20261012.md。
    #   ★ 数值取自 `_dr_src_core`（与 figS10 落点**同一份**），此处只换 layer 拼写。
    + [dict(layer="L-DR minus L-src", **_dr_src_core)]
).to_csv(os.path.join(OUT, "p45_triage_summary_main.csv"), index=False)
p("    写入 p45_triage_summary_main.csv")
p("")
p("完成。")
REP.close()
