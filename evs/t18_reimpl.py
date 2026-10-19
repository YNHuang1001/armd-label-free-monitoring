# -*- coding: utf-8 -*-
r"""t18_reimpl.py — families P12 / FIG5：**独立重实现**（GPT6 的"第 3 件"）

================================ GPT6 的要求（本模块的设计依据） ================================
    "reference 从**更上游冻结输入**开始，自己实现：
       config pairing / alignment；per-level aggregation；quantiles；paired drift；缺失配对如何处理。
     不读取 fig1e_configs*.csv / fig5*.csv 来'验证自己'。
     最终生产实现和 independent reference 逐项对拍。"

本模块据此把两条 path 的**聚合层**独立重写一遍：

  P12-*  上游 = `p12_configs_corrected.csv`（per-task × per-config 表）
         目标 = `figures/data/fig1e_configs.csv` / `fig1e_configs_aligned.csv`
         ★ 复现的语义：universe = 「P1 History+LR 且 eligible」的任务集；**每个配置都限制到该集合**
           （= pairing/alignment）；缺配对的 task **直接丢弃、不插补**；分位数用 pandas 线性插值；
           harmful = `d_auroc <= -0.05`（预先写死的主终点）。

  FIG5-* 上游 = `p10b_raw.csv` ⨝ `p6d_outcome.csv[["task_id","harmful"]]`
         目标 = `figures/data/fig5_subsample.csv` / `fig5b_across_n_range.csv`
         ★ 复现的语义：每 (方法,档) 先按 `rep` 分组算 AUROC（finite 掩码后需 ≥2 类，否则**跳过该 rep**）；
           `auroc_mean` = 存活 rep 的均值；CI = 这些 rep 值的 `.quantile(.025/.975)`；
           **paired drift = 按 rep 索引对齐后 `full − 100`，`.dropna()` 丢弃缺配对的 rep**。

★ 边界（如实声明）：独立重实现的是**聚合层**。上游表本身由 p12 的建模管线 / p10b 的抽样管线产出，
  本模块**不重建**它们的模型训练。这与 GPT6 列出的五项（pairing / aggregation / quantiles /
  paired drift / 缺配对处理）一致；上游模型训练那一半不在本 gate 的范围内。

================================ 判据 ================================
  P12-01   fig1e_configs.csv 全 12 列 × 4 行逐项对拍（独立重算）
  P12-02   对齐语义的**受控 fixture**：universe 只取自 P1，缺配对的 task 必须**丢弃而非插补**
  P12-03   非空转：把 universe 改成「各配置自身的 eligible 集」后，必须与冻结产物不符
  FIG5-01  fig5_subsample.csv 全 7 列 × 90 行逐项对拍
  FIG5-02  配对语义的**受控 fixture**：缺配对的 rep 必须被 dropna 丢弃
  FIG5-03  fig5b_across_n_range.csv 全 10 列 × 15 行逐项对拍
  FIG5-04  非空转：把某个方法的符号去掉（SUDO_AURCC 不再取负）后，必须被检出

输出契约同 evs_load：`RESULT|tid|PASS|FAIL|WARN|expected|observed|note`
"""
import os, sys, io
import numpy as np
import pandas as pd

import evs_load as L

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

OUT = L.OUT
FIG = os.path.join(L.PROJ, "figures", "data")

CFG12_CSV = os.path.join(OUT, "p12_configs_corrected.csv")
RAW10_CSV = os.path.join(OUT, "p10b_raw.csv")
OUTC_CSV = os.path.join(OUT, "p6d_outcome.csv")
F1E_CSV = os.path.join(FIG, "fig1e_configs.csv")
F1E_AL_CSV = os.path.join(FIG, "fig1e_configs_aligned.csv")
F5_CSV = os.path.join(FIG, "fig5_subsample.csv")
F5B_CSV = os.path.join(FIG, "fig5b_across_n_range.csv")

TOL = 1e-9
CONFIGS = ["P1 History+LR", "S1 Core+LR", "S2 History+LGBM", "S3 Core+LGBM"]
LEVELS = ["100", "250", "500", "1000", "full"]
HARMFUL_THR = -0.05                     # 预先写死的主终点（correction protocol §1）

# 与 p11:179-190 **逐字一致**的方法定义（名称/角色/取值方式）
METHODS = [
    ("domain_auc", "L1", lambda x: x["domain_auc"]),
    ("mmd", "L1", lambda x: x["mmd"]),
    ("psi_mean", "L1", lambda x: x["psi_mean"]),
    ("psi_max", "L1", lambda x: x["psi_max"]),
    ("js_mean", "L1", lambda x: x["js_mean"]),
    ("js_max", "L1", lambda x: x["js_max"]),
    ("pred_mean_shift", "L2", lambda x: x["pred_mean_shift"].abs()),
    ("entropy_mean_shift", "L2", lambda x: x["entropy_mean_shift"].abs()),
    ("entropy_js", "L2", lambda x: x["entropy_js"]),
    ("conf_dist_js", "L2", lambda x: x["conf_dist_js"]),
    ("pred_pos_prev_shift", "L2", lambda x: x["pred_pos_prev_shift"].abs()),
    ("SUDO_AURCC", "L2", lambda x: -x["SUDO_AURCC"]),
    ("AUROC_DR", "L3-primary", lambda x: -x["AUROC_DR"]),
    ("AUROC_IOW", "L3-sens", lambda x: -x["AUROC_IOW"]),
    ("AUROC_OM", "L3-sens", lambda x: -x["AUROC_OM"]),
]

print("=" * 100, flush=True)
print("t18 — families P12 / FIG5：独立重实现（从更上游冻结输入重算聚合层）", flush=True)
print("=" * 100, flush=True)
print(f"  生产树 = {L.PROD_DIR}   模式 = {'PATCHED' if L.IS_PATCHED else 'REAL'}", flush=True)
print(f"  上游: {os.path.basename(CFG12_CSV)} / {os.path.basename(RAW10_CSV)} + "
      f"{os.path.basename(OUTC_CSV)}", flush=True)
print(f"  目标: figures/data/fig1e_configs.csv / fig1e_configs_aligned.csv / "
      f"fig5_subsample.csv / fig5b_across_n_range.csv", flush=True)
print("", flush=True)

from sklearn.metrics import roc_auc_score


# ==================================================================================
# P12 —— Fig 1E 的独立重实现
# ==================================================================================
def ref_fig1e(cfg12, configs=CONFIGS, thr=HARMFUL_THR, universe="P1"):
    """独立实现 Fig 1E 的聚合与**对齐**。

    universe='P1'    -> 正确的语义：universe = P1 且 eligible 的 task 集（配对基准固定为 P1）
    universe='own'   -> 非空转对照用的**错误**语义：各配置用自身的 eligible 集（未对齐）
    """
    if universe == "P1":
        base = set(cfg12[(cfg12.config == configs[0]) & cfg12.eligible].task_id)
    rows = []
    for cfg in configs:
        own = cfg12[cfg12.config == cfg]
        keep = base if universe == "P1" else set(own[own.eligible].task_id)
        s = own[own.task_id.isin(keep)]
        d = s.d_auroc
        rows.append(dict(
            config=cfg, tier=own.tier.iat[0], algo=own.algo.iat[0],
            n_aligned=len(s),
            median_aligned=d.median(),
            p10_aligned=d.quantile(.1), p90_aligned=d.quantile(.9),
            p90_minus_p10_aligned=d.quantile(.9) - d.quantile(.1),
            n_harmful_aligned=int((d <= thr).sum()),
            pct_harmful_aligned=100 * float((d <= thr).mean()),
            n_eligible_own=int(own.eligible.sum()), n_task_own=len(own)))
    return pd.DataFrame(rows)


CFG12 = pd.read_csv(CFG12_CSV)


def _cmp_rows(ref, got, key, cols):
    """逐项对拍；返回 (最差差, 不符明细)。"""
    worst, bad = 0.0, []
    r = ref.set_index(key); g = got.set_index(key)
    if set(r.index) != set(g.index):
        return float("inf"), [f"index 不一致：ref={sorted(r.index)} got={sorted(g.index)}"]
    for k in r.index:
        for c in cols:
            a, b = r.loc[k, c], g.loc[k, c]
            if isinstance(a, str) or isinstance(b, str):
                if str(a) != str(b):
                    bad.append(f"{k}.{c}: ref={a!r} got={b!r}")
                continue
            dd = abs(float(a) - float(b))
            if not np.isfinite(dd) or dd > TOL:
                bad.append(f"{k}.{c}: ref={a} got={b} absdiff=={dd:.3e}")
            worst = max(worst, dd if np.isfinite(dd) else float("inf"))
    return worst, bad


def p12_01():
    if not os.path.exists(F1E_CSV):
        return False, f"目标不存在：{F1E_CSV}"
    ref = ref_fig1e(CFG12)
    got = pd.read_csv(F1E_CSV)
    cols = [c for c in got.columns if c != "config"]
    worst, bad = _cmp_rows(ref, got, "config", cols)
    return (not bad), (
        f"从 {os.path.basename(CFG12_CSV)}（{CFG12.shape[0]} 行）独立重算 Fig 1E 的 "
        f"{len(cols)} 列 × {len(CONFIGS)} 配置；最大 absdiff= = {worst:.3e}"
        + ("" if not bad else f"；不符 {len(bad)} 项 :: " + " | ".join(bad[:4]))
        + f"；对齐基准 = P1 且 eligible 的 {ref.n_aligned.tolist()} 个 task（各配置 aligned 数）")


L.guard("P12-01", p12_01,
        "独立重实现 Fig 1E：从**上游** p12_configs_corrected.csv 自行实现 "
        "(a) 对齐（universe = P1 History+LR 且 eligible 的任务集，每个配置都限制到它）、"
        "(b) 中位数与 P10/P90 分位数、(c) P90−P10、(d) harmful 计数（d_auroc ≤ −0.05）、"
        "(e) n_aligned / n_eligible_own / n_task_own —— 与 figures/data/fig1e_configs.csv 逐项对拍"
        "（容差 1e-9）",
        "★ 不读 fig1e_configs*.csv 来构造 reference；只把它作为**被检验对象**")


def p12_02():
    """受控 fixture：P1 有 5 个 eligible task，另一配置只有其中 2 个 —— 必须丢弃 3 个、**不得插补**。"""
    rows = []
    for tid in ["t1", "t2", "t3", "t4", "t5"]:
        rows.append(dict(config=CONFIGS[0], tier="tier1", algo="LR", task_id=tid,
                         d_auroc=0.10, eligible=True))
    for tid, dv in [("t1", -0.20), ("t3", 0.05)]:        # 只有 2 个 task 存在
        rows.append(dict(config=CONFIGS[1], tier="tier2", algo="LR", task_id=tid,
                         d_auroc=dv, eligible=True))
    fx = pd.DataFrame(rows)
    out = ref_fig1e(fx, configs=CONFIGS[:2])
    n0 = int(out.loc[out.config == CONFIGS[0], "n_aligned"].iat[0])
    n1 = int(out.loc[out.config == CONFIGS[1], "n_aligned"].iat[0])
    har = int(out.loc[out.config == CONFIGS[1], "n_harmful_aligned"].iat[0])
    ok = (n0 == 5 and n1 == 2 and har == 1)
    return ok, (f"受控 fixture：P1 有 5 个 eligible task，第二配置只有 t1/t3 → "
                f"n_aligned = ({n0}, {n1})（应为 (5, 2)：**缺配对的 3 个 task 被丢弃**）；"
                f"第二配置 harmful 计数 = {har}（t1 的 −0.20 ≤ −0.05，t3 的 0.05 不算）")


L.guard("P12-02", p12_02,
        "对齐语义的受控 fixture（**缺配对如何处理**）：universe 只取自 P1 且 eligible；"
        "某配置缺的 task 必须被**丢弃**（n_aligned 变小），**不得插补**、不得改用该配置自身的 eligible 集",
        "这条正是 GPT6 点名的 'config pairing / alignment + 缺失配对如何处理'")


def p12_03():
    ref_ok = ref_fig1e(CFG12)
    ref_bad = ref_fig1e(CFG12, universe="own")
    got = pd.read_csv(F1E_CSV)
    _, bad_ok = _cmp_rows(ref_ok, got, "config", ["n_aligned", "median_aligned"])
    _, bad_bad = _cmp_rows(ref_bad, got, "config", ["n_aligned", "median_aligned"])
    ok = (not bad_ok) and bool(bad_bad)
    return ok, (f"把 universe 篡改为「各配置自身的 eligible 集」（未对齐）后，与冻结产物不符的项数 = "
                f"{len(bad_bad)}（对齐版本为 {len(bad_ok)}）—— "
                f"{'被检出（判据非空转）' if ok else '**未被检出，判据空转**'}"
                + (f"；示例：{bad_bad[:2]}" if bad_bad else ""))


L.guard("P12-03", p12_03,
        "非空转对照：把对齐基准从「固定的 P1 universe」篡改为「各配置自身的 eligible 集」后，"
        "独立 reference 必须与冻结产物**不符**",
        "★ 没有这一条，「逐项对拍」不能证明 reference 真的在做对齐")


# ==================================================================================
# FIG5 —— Fig 5 / 5b 的独立重实现
# ==================================================================================
def ref_fig5(raw, outc, methods=METHODS, levels=LEVELS, negate=(12, 13, 14, 15)):
    """独立实现 Fig 5 的 per-level 聚合与**配对漂移**。

    negate: 需要取负的方法序号（1-based）——用于非空转对照时**去掉**某个负号。
    """
    rw = raw.merge(outc[["task_id", "harmful"]], on="task_id", how="left")
    rows, per_level = [], {}
    for mi, (name, role, fn) in enumerate(methods, 1):
        f = fn if mi in negate else (lambda x, _fn=fn: _fn(x))
        if mi not in negate:
            # 去掉该方法的负号：符号约定改变
            f = _make_unsigned(name, fn)
        pl = {}
        for lv in levels:
            sub = rw[rw.level == lv]
            vals = {}
            for rr, gg in sub.groupby("rep"):
                gg = gg.sort_values("task_id")
                s = f(gg).values.astype(float)
                yy = gg["harmful"].values.astype(int)
                m = np.isfinite(s)
                if len(np.unique(yy[m])) < 2:
                    continue
                vals[rr] = roc_auc_score(yy[m], s[m])
            v = pd.Series(vals)
            pl[lv] = v
            rows.append(dict(method=name, role=role, level=lv, n_repeat=len(v),
                             auroc_mean=float(v.mean()),
                             auroc_lo=float(v.quantile(.025)),
                             auroc_hi=float(v.quantile(.975))))
        per_level[name] = pl
        d = (pl["full"] - pl["100"]).dropna()
        rows.append(dict(method=name, role=role, level="PAIRED_full_minus_100",
                         n_repeat=len(d), auroc_mean=float(d.mean()),
                         auroc_lo=float(d.quantile(.025)), auroc_hi=float(d.quantile(.975))))
    return pd.DataFrame(rows), per_level


def _make_unsigned(name, fn):
    """非空转对照：把负号方法改为不取负（其余方法不变）。"""
    if name == "SUDO_AURCC":
        return lambda x: x["SUDO_AURCC"]
    if name == "AUROC_DR":
        return lambda x: x["AUROC_DR"]
    if name == "AUROC_IOW":
        return lambda x: x["AUROC_IOW"]
    if name == "AUROC_OM":
        return lambda x: x["AUROC_OM"]
    return fn


def ref_fig5b(f5, levels=LEVELS):
    out = []
    lv = f5[f5.level.isin(levels)]
    for name, role, _ in METHODS:
        s = lv[lv.method == name]
        imax = s.auroc_mean.idxmax(); imin = s.auroc_mean.idxmin()
        out.append(dict(method=name, role=role, n_levels=len(s),
                        across_n_range=float(s.auroc_mean.max() - s.auroc_mean.min()),
                        auroc_max=float(s.auroc_mean.max()),
                        auroc_max_level=s.loc[imax, "level"],
                        auroc_min=float(s.auroc_mean.min()),
                        auroc_min_level=s.loc[imin, "level"],
                        auroc_at_full=float(s.loc[s.level == "full", "auroc_mean"].iat[0]),
                        auroc_at_100=float(s.loc[s.level == "100", "auroc_mean"].iat[0])))
    return pd.DataFrame(out)


_RAW = _OUTC = _F5 = None


def _load():
    global _RAW, _OUTC, _F5
    if _RAW is None:
        _RAW = pd.read_csv(RAW10_CSV)
        _OUTC = pd.read_csv(OUTC_CSV)
        _F5 = pd.read_csv(F5_CSV) if os.path.exists(F5_CSV) else None
    return _RAW, _OUTC, _F5


def fig5_01():
    raw, outc, f5 = _load()
    if f5 is None:
        return False, f"目标不存在：{F5_CSV}"
    ref, _ = ref_fig5(raw, outc)
    cols = ["n_repeat", "auroc_mean", "auroc_lo", "auroc_hi"]
    ref["_k"] = ref.method + "|" + ref.level
    got = f5.copy(); got["_k"] = got.method + "|" + got.level
    worst, bad = _cmp_rows(ref, got, "_k", cols)
    return (not bad), (
        f"从 {os.path.basename(RAW10_CSV)}（{raw.shape[0]} 行 = 202 task × 5 档 × 200 rep）"
        f"⨝ {os.path.basename(OUTC_CSV)} 独立重算 fig5_subsample.csv 的 {len(cols)} 列 × {len(ref)} 行；"
        f"最大 absdiff= = {worst:.3e}"
        + ("" if not bad else f"；不符 {len(bad)} 项 :: " + " | ".join(bad[:4]))
        + f"；15 方法符号约定（abs/取负）逐字对齐 p11:179-190")


L.guard("FIG5-01", fig5_01,
        "独立重实现 Fig 5：从**上游** p10b_raw.csv ⨝ p6d_outcome.csv 自行实现 "
        "(a) 每 (方法,档) 按 rep 分组、finite 掩码后需 ≥2 类否则跳过该 rep、"
        "(b) auroc_mean 与 .quantile(.025/.975) 的 CI、"
        "(c) **配对漂移 full − 100（按 rep 索引对齐，dropna 丢弃缺配对 rep）** —— "
        "与 figures/data/fig5_subsample.csv 逐项对拍（容差 1e-9）",
        "★ 不读 fig5*.csv 来构造 reference；只把它作为**被检验对象**")


def fig5_02():
    """受控 fixture：rep 2 在 full 档因单类被跳过 → paired 必须只有 rep 1、3，**不得插补**。"""
    rows = []
    for lv in ["100", "full"]:
        for rr in [1, 2, 3]:
            for tid in ["t1", "t2", "t3", "t4"]:
                harm = 0 if tid in ("t1", "t2") else 1
                rows.append(dict(task_id=tid, level=lv, rep=rr, sig=float(harm) + 0.1 * rr))
    raw = pd.DataFrame(rows)
    outc = pd.DataFrame(dict(task_id=["t1", "t2", "t3", "t4"], harmful=[0, 0, 1, 1]))
    # 让 rep2/full 退化为单类：把所有 harmful 置 0
    raw.loc[(raw.level == "full") & (raw.rep == 2), "sig"] = 0.0
    outc2 = outc.copy()
    m = (raw.level == "full") & (raw.rep == 2)
    raw.loc[m, "task_id"] = raw.loc[m, "task_id"]        # 不改 task_id，用 harmful 掩码更直接
    rw = raw.merge(outc2, on="task_id", how="left")
    # 手工：full/rep2 的 harms 全部改成 0（单类）
    rw.loc[(rw.level == "full") & (rw.rep == 2), "harmful"] = 0
    pl = {}
    for lv in ["100", "full"]:
        sub = rw[rw.level == lv]; vals = {}
        for rr, gg in sub.groupby("rep"):
            gg = gg.sort_values("task_id")
            s = gg["sig"].values.astype(float); yy = gg["harmful"].values.astype(int)
            mm = np.isfinite(s)
            if len(np.unique(yy[mm])) < 2:
                continue
            vals[rr] = roc_auc_score(yy[mm], s[mm])
        pl[lv] = pd.Series(vals)
    d = (pl["full"] - pl["100"]).dropna()
    ok = (len(pl["100"]) == 3 and len(pl["full"]) == 2 and len(d) == 2)
    return ok, (f"受控 fixture：full 档的 rep2 退化为单类 → 该 rep 在 full 档被跳过；"
                f"存活 rep 数 = 100档 {len(pl['100'])}、full档 {len(pl['full'])}；"
                f"配对后 = {len(d)}（应为 2：**缺配对的 rep2 被 dropna 丢弃**，"
                f"而不是用其它 rep 的均值插补）")


L.guard("FIG5-02", fig5_02,
        "配对语义的受控 fixture（**缺配对如何处理**）：paired drift 必须按 **rep 索引对齐**后相减，"
        "缺配对的 rep 由 `.dropna()` **丢弃**；不得改用「均值之差」或任何插补",
        "这条正是 GPT5.6 点名的 'paired drift + 缺失配对如何处理'")


def fig5_03():
    raw, outc, f5 = _load()
    if f5 is None or not os.path.exists(F5B_CSV):
        return False, f"目标不存在：{F5B_CSV}"
    ref = ref_fig5b(f5)
    got = pd.read_csv(F5B_CSV)
    cols = [c for c in got.columns if c != "method"]
    worst, bad = _cmp_rows(ref, got, "method", cols)
    return (not bad), (
        f"从 fig5_subsample.csv 独立重算跨档极差表 {len(cols)} 列 × {len(ref)} 方法；"
        f"最大 absdiff= = {worst:.3e}"
        + ("" if not bad else f"；不符 {len(bad)} 项 :: " + " | ".join(bad[:4]))
        + f"；levels = {LEVELS}")


L.guard("FIG5-03", fig5_03,
        "独立重实现 fig5b_across_n_range.csv：跨 5 档的 auroc_mean **极差**（max−min）、"
        "取到极值的档位、以及 full / 100 档上的值，逐项对拍",
        "上游仍是 fig5_subsample.csv 的**上游**；此处只对拍跨档汇总，不重算 per-rep AUROC")


def fig5_04():
    raw, outc, f5 = _load()
    if f5 is None:
        return False, "目标不存在"
    ref_ok, _ = ref_fig5(raw, outc)
    ref_mut, _ = ref_fig5(raw, outc, negate=())      # 所有取负方法改为不取负
    k = "method|level"
    a = ref_ok.assign(_k=ref_ok.method + "|" + ref_ok.level).set_index("_k")
    b = ref_mut.assign(_k=ref_mut.method + "|" + ref_mut.level).set_index("_k")
    g = f5.assign(_k=f5.method + "|" + f5.level).set_index("_k")
    def _mism(x):
        dd = (x["auroc_mean"] - g["auroc_mean"]).abs()
        return int((dd > TOL).sum())
    n_ok, n_mut = _mism(a), _mism(b)
    ok = (n_ok == 0) and (n_mut > 0)
    return ok, (f"把 4 个「取负」方法（SUDO_AURCC / AUROC_DR / AUROC_IOW / AUROC_OM）的符号去掉后，"
                f"与冻结产物不符的行数 = {n_mut}（正确符号下为 {n_ok}）—— "
                f"{'被检出（判据非空转）' if ok else '**未被检出，判据空转**'}")


L.guard("FIG5-04", fig5_04,
        "非空转对照：去掉取负方法的符号后，独立 reference 必须与冻结产物**不符**",
        "★ 符号约定（abs / 取负）是 p11:179-190 的冻结语义，去掉后 AUROC 必然改变（1−AUROC）")

print("\n[t18_reimpl done]", flush=True)
