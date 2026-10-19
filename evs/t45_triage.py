# -*- coding: utf-8 -*-
r"""t45_triage.py — 分诊策略模拟（p45）的**独立**验证

================================ 为什么需要 ================================
`p45_triage_policy_simulation.py` 是**事后**分析，但它产出的数字要进稿件
（Supplementary Figure S10 / Table S3）。它是**唯一**决定这些数字的地方，
因此必须落在**独立**判据之下，而不是靠 p45 自己的内部对账（那是自证）。

本模块**不 import、不 exec、不 subprocess** p45 的任何代码。

================================ 独立性从哪来 ================================
  1. 自己从**更上游的冻结表**（p6b / p6c / p6d_outcome）重建排名与命中数；
  2. 自己实现平均秩、排序、并列处理（纯 Python，不用 numpy/pandas 的排序）；
  3. 与 p45 落盘的 CSV 逐项比对。
p45 改了排名规则而 CSV 没同步 → 判红；CSV 漂移 → 判红。

================================ 覆盖 ================================
  T45-01  预注册备忘哈希未变（防「看完结果再改备忘」）
  T45-02  锚定：独立重算 15 个 harmful_auc == p6d_signal_eval.csv（逐位）
  T45-02b source 内部 AUROC 的两表一致性（p6b `source_auroc_internal`
          vs p6d `auroc_source_internal`，逐位）＋ 值域 [0,1]
  T45-03  k=20% 的命中数逐项（21 项：15 信号 + 5 层 + 随机）
  T45-04  累积增益曲线内部一致（单调 / 终值 / precision 定义）
  T45-05  主对比自洽（点估计 = 两层之差；bootstrap 区间含点估计）
  T45-06  L-shift 秩平均定义的独立重建
  T45-07  六个 L1 成员逐个命中数与 summary 一致（增补文件所依赖的稳健性披露）
  T45-08  互相闭合：点估计与 bootstrap 均值在报告文本中分列且各标各的数
  T45-09  **事后**对比（DR − src-only）自洽 ＋ 两个落点逐位一致
          （同一量在 figS10 与 summary_main 各有落点 —— 防 INC-015 形态）

★ 2026-10-12：L-src / L-dhat 两层是本轮事后追加的（原判据不认识它们，T45-03 报
  「未知条目」）。**处置不是把它们加进白名单，而是独立重建它们**：两层都从
  p6b 的 `source_auroc_internal` 出发，本模块不 import p45 的任何代码。
"""
import csv
import hashlib
import io
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import evs_load as L

PROJ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(PROJ, "分析结果")

MEMO = os.path.join(PROJ, "TRIAGE_POLICY_PREREGISTRATION_20261007.md")
MEMO_SHA = "5e25a14963534d808f7085c9d9fe6d540c13c45f57e974eee87bf6e4dfb98fef"

L1 = ["domain_auc", "mmd", "psi_mean", "psi_max", "js_mean", "js_max"]
SIGNED_L2 = ["pred_mean_shift", "entropy_mean_shift", "pred_pos_prev_shift"]
L3 = ["AUROC_DR", "AUROC_IOW", "AUROC_OM"]
ALL15 = L1 + SIGNED_L2 + ["entropy_js", "conf_dist_js", "SUDO_AURCC"] + L3
K_FRAC = 0.20


def rd(fn):
    p = os.path.join(OUT, fn)
    if not os.path.exists(p):
        return None
    with io.open(p, encoding="utf-8") as f:
        return {r["task_id"]: r for r in csv.DictReader(f)}


B_ = rd("p6b_label_free_signals_blinded.csv")
C_ = rd("p6c_L3_blinded.csv")
D_ = rd("p6d_outcome.csv")
CURVE = None
SUMM = []
MAIN = []

_p = os.path.join(OUT, "p45_triage_curves.csv")
if os.path.exists(_p):
    with io.open(_p, encoding="utf-8") as f:
        CURVE = list(csv.DictReader(f))
_p = os.path.join(OUT, "p45_triage_summary.csv")
if os.path.exists(_p):
    with io.open(_p, encoding="utf-8") as f:
        SUMM = list(csv.DictReader(f))
_p = os.path.join(OUT, "p45_triage_summary_main.csv")
if os.path.exists(_p):
    with io.open(_p, encoding="utf-8") as f:
        MAIN = list(csv.DictReader(f))


# ---------------------------------------------------------------- 独立实现
def val(m, t):
    v = float(B_[t][m]) if m in B_[t] else float(C_[t][m])
    if m in SIGNED_L2:
        return abs(v)
    if m in ("SUDO_AURCC",) + tuple(L3):
        return -v
    return v


def avg_rank(vals):
    """平均秩（1 = 最小）。纯 Python，独立于 numpy/pandas。"""
    order = sorted(range(len(vals)), key=lambda i: vals[i])
    r = [0.0] * len(vals)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and vals[order[j + 1]] == vals[order[i]]:
            j += 1
        a = 0.5 * (i + j) + 1
        for k in range(i, j + 1):
            r[order[k]] = a
        i = j + 1
    return r


def caught(score, tasks, harm, k, tie="up"):
    idx = list(range(len(tasks)))
    idx.sort(key=lambda i: (-score[i], tasks[i] if tie == "up" else "\uffff" + tasks[i]))
    return sum(harm[tasks[i]] for i in idx[:k])


def wmw_auc(yv, sv):
    """Mann–Whitney，平均秩处理并列。"""
    pairs = sorted(zip(sv, yv))
    n1 = sum(yv); n0 = len(yv) - n1
    if n1 == 0 or n0 == 0:
        return float("nan")
    r = [0.0] * len(pairs); i = 0
    while i < len(pairs):
        j = i
        while j + 1 < len(pairs) and pairs[j + 1][0] == pairs[i][0]:
            j += 1
        a = 0.5 * (i + j) + 1
        for k in range(i, j + 1):
            r[k] = a
        i = j + 1
    return (sum(r[k] for k in range(len(pairs)) if pairs[k][1] == 1) - n1 * (n1 + 1) / 2) / (n1 * n0)


# ---------------------------------------------------------------- T45-01 备忘哈希
if os.path.exists(MEMO):
    h = hashlib.sha256(open(MEMO, "rb").read()).hexdigest()
    L.check("T45-01", h == MEMO_SHA,
            "预注册备忘 sha256 == %s" % MEMO_SHA[:16],
            "实测 %s" % h[:16],
            "哈希不符即判红：说明预注册在被计算之前就已不是当初写死的那份")
else:
    L.check("T45-01", False, "预注册备忘存在", "文件缺失", "缺失即无法证明事前性")

if not (B_ and C_ and D_ and CURVE and SUMM and MAIN):
    L.check("T45-02..07", False, "上游冻结表与 p45 产物齐备",
            "缺失：%s" % [n for n, v in [("p6b", B_), ("p6c", C_), ("p6d", D_),
                                          ("curves", CURVE), ("summary", SUMM),
                                          ("main", MAIN)] if not v],
            "缺件即判红，不静默跳过")
else:
    tasks = sorted(D_)
    harm = {t: int(D_[t]["harmful"]) for t in tasks}
    N = len(tasks)
    K = int(round(K_FRAC * N))
    NFAIL = sum(harm.values())
    PREV = NFAIL / N

    # ------------------------------------------------------------ T45-02 锚定
    fr = {}
    with io.open(os.path.join(OUT, "p6d_signal_eval.csv"), encoding="utf-8") as f:
        for r in csv.DictReader(f):
            fr[r["method"]] = float(r["harmful_auc"])
    bad = []
    for m in ALL15:
        a = wmw_auc([harm[t] for t in tasks], [val(m, t) for t in tasks])
        if abs(a - fr[m]) > 1e-9:
            bad.append("%s(%.6f≠%.6f)" % (m, a, fr[m]))
    L.check("T45-02", not bad,
            "独立重算 15 个 harmful_auc == p6d_signal_eval.csv（逐位）",
            "不符 %d 项%s" % (len(bad), ("：" + "; ".join(bad[:3])) if bad else ""),
            "锚定失败即判红：p45 的输入将与冻结主分析不同源")

    # ------------------------------------------------- T45-02b source 内部 AUROC
    # ★ 2026-10-12 新增。两层事后基线（L-src / L-dhat）**共同**依赖 source 自身的内部
    #   AUROC。该量出现在两张冻结表里，拼写不同（p6b `source_auroc_internal`、
    #   p6d `auroc_source_internal`）—— 同一量两个名字即两个漂移机会，必须逐位对齐。
    #   值域 [0,1] 是 AUROC 的定义域。
    bad2b = []
    for t in tasks:
        bv = float(B_[t]["source_auroc_internal"])
        dv = float(D_[t]["auroc_source_internal"])
        if abs(bv - dv) > 1e-12:
            bad2b.append("%s(p6b %.9g ≠ p6d %.9g)" % (t, bv, dv))
        if not (0.0 <= bv <= 1.0):
            bad2b.append("%s(越界 %.6f)" % (t, bv))
    L.check("T45-02b", not bad2b,
            "p6b.source_auroc_internal == p6d.auroc_source_internal 且均落在 [0,1]",
            "不符 %d 项%s" % (len(bad2b), ("：" + "; ".join(bad2b[:3])) if bad2b else ""),
            "L-src / L-dhat 两层都建立在 source 内部 AUROC 上；两表若漂移，事后基线即失锚")

    # ------------------------------------------------------------ T45-06 秩平均定义
    R6 = [avg_rank([float(B_[t][m]) for t in tasks]) for m in L1]
    shift = [sum(R6[j][i] for j in range(6)) / 6.0 for i in range(N)]
    st = sorted([(shift[i], tasks[i]) for i in range(N)], key=lambda x: (-x[0], x[1]))
    top_ids = [t for _, t in st[:K]]
    # 与 p45 summary 中 L-shift 的命中数比对
    ls_row = [r for r in SUMM if r["item"] == "L-shift" and abs(float(r["k_frac"]) - K_FRAC) < 1e-9]
    ok6 = bool(ls_row) and abs(sum(harm[t] for t in top_ids) - float(ls_row[0]["caught"])) < 1e-9
    L.check("T45-06", ok6,
            "L-shift 秩平均（六者平均秩）独立重建的命中数 == summary 中的值",
            "独立 %d｜summary %s" % (sum(harm[t] for t in top_ids),
                                     ls_row[0]["caught"] if ls_row else "缺行"),
            "秩平均的定义是主对照的一半，定义错了整条结论失效")

    # ------------------------------------------------------------ T45-03 逐项 precision
    # ★ 2026-10-12 追加两层（p45 事后新增，见 TRIAGE_POLICY_PREREGISTRATION_AMENDMENT_20261012.md）。
    #   独立重建，**不 import p45**：
    #     L-src  —— 风险分就是 source 自身内部 AUROC（高者先复核，与 p45 同向）。
    #                方向由数据裁定：以该列作 harm 分，AUROC = 0.6107（取负则为 0.3893），
    #                ⇒ 高分=高危 是唯一使 src-only AUROC > 0.5 的方向。
    #     L-dhat —— 事后对齐分 −(Â_target − AUROC_source)。本模块的 `val("AUROC_DR", t)`
    #                已经返回 −Â_target ⇒ 直接加回 source 内部 AUROC 即得。
    src_int = [float(B_[t]["source_auroc_internal"]) for t in tasks]
    layer = {"L-shift": shift,
             "L-DR": [val("AUROC_DR", t) for t in tasks],
             "L-OM": [val("AUROC_OM", t) for t in tasks],
             "L-src": src_int,
             "L-dhat": [val("AUROC_DR", t) + src_int[i] for i, t in enumerate(tasks)]}
    exp_map = {}
    for m in ALL15:
        exp_map[m] = caught([val(m, t) for t in tasks], tasks, harm, K)
    for nm, s in layer.items():
        exp_map[nm] = caught(s, tasks, harm, K)
    bad3 = []
    seen = set()
    for r in SUMM:
        if abs(float(r["k_frac"]) - K_FRAC) > 1e-9:
            continue
        it = r["item"]
        if it in ("L-random",):
            if abs(float(r["precision"]) - PREV) > 1e-9:
                bad3.append("L-random precision")
            continue
        if it not in exp_map:
            bad3.append("未知条目 " + it); continue
        seen.add(it)
        if abs(float(r["caught"]) - exp_map[it]) > 1e-9:
            bad3.append("%s(%s≠%d)" % (it, r["caught"], exp_map[it]))
    miss = set(exp_map) - seen
    if miss:
        bad3.append("summary 缺条目 " + ",".join(sorted(miss)))
    L.check("T45-03", not bad3,
            "k=%d 的命中数逐项（%d 项）== 独立重算" % (K, len(exp_map)),
            "不符 %d 项%s" % (len(bad3), ("：" + "; ".join(bad3[:3])) if bad3 else ""),
            "这是要进 S10 / Table S3 的数字")

    # ------------------------------------------------------------ T45-04 曲线自洽
    bad4 = []
    for lay in set(r["layer"] for r in CURVE):
        rows = sorted([r for r in CURVE if r["layer"] == lay], key=lambda r: int(r["k"]))
        if [int(r["k"]) for r in rows] != list(range(1, N + 1)):
            bad4.append("%s:k 序列不完整" % lay); continue
        if abs(float(rows[-1]["cum_failures"]) - (NFAIL if lay != "L-random" else NFAIL)) > 1e-6:
            bad4.append("%s:终值 %s" % (lay, rows[-1]["cum_failures"]))
        prev = -1
        for r in rows:
            c, k = float(r["cum_failures"]), int(r["k"])
            if c < prev - 1e-9:
                bad4.append("%s:k=%d 非单调" % (lay, k)); break
            prev = c
            if abs(float(r["precision"]) - c / k) > 1e-9:
                bad4.append("%s:k=%d precision 定义不符" % (lay, k)); break
            if abs(float(r["recall"]) - c / NFAIL) > 1e-9:
                bad4.append("%s:k=%d recall 定义不符" % (lay, k)); break
        if bad4:
            break
    L.check("T45-04", not bad4,
            "累积增益曲线自洽（k 完整 1..%d、cum 单调、precision=cum/k、recall=cum/57）" % N,
            "不符 %d 项%s" % (len(bad4), ("：" + "; ".join(bad4[:3])) if bad4 else ""),
            "曲线是 S10 的唯一数据源")

    # ------------------------------------------------------------ T45-05 主对比自洽
    bad5 = []
    mrow = [r for r in MAIN if r["layer"].startswith("L-DR minus")]
    dr = [r for r in MAIN if r["layer"] == "L-DR"]
    ls = [r for r in MAIN if r["layer"] == "L-shift"]
    if not (mrow and dr and ls):
        bad5.append("summary_main 缺行")
    else:
        pt = float(dr[0]["precision_point"]) - float(ls[0]["precision_point"])
        if abs(pt - float(mrow[0]["precision_point"])) > 1e-9:
            bad5.append("点估计 %.6f ≠ 两层之差 %.6f" % (float(mrow[0]["precision_point"]), pt))
        lo, hi = float(mrow[0]["lo95"]), float(mrow[0]["hi95"])
        if not (lo < float(mrow[0]["precision_boot_mean"]) < hi):
            bad5.append("bootstrap 区间不含均值")
        if not (lo <= pt <= hi):
            bad5.append("点估计落在区间外（%.4f not in [%.4f,%.4f]）" % (pt, lo, hi))
        if not (0 < float(mrow[0]["lo95"])):
            L.emit("T45-05b", "WARN", "主对比 95% 下界 > 0（预注册的主对比）",
                   "下界 = %.4f" % lo, "下界 ≤ 0 时**不得**在正文声称主对比显著")
        else:
            L.emit("T45-05b", "PASS", "主对比 95% 下界 > 0", "下界 = %.4f" % lo,
                   "预注册主对比 L-DR − L-shift")
    L.check("T45-05", not bad5, "主对比自洽（点估计 = 两层之差；bootstrap 区间含点估计与均值）",
            "不符 %d 项%s" % (len(bad5), ("：" + "; ".join(bad5[:2])) if bad5 else ""),
            "这是全文唯一作显著性声称的对比")

    # ------------------------------------------------------------ T45-07 六个 L1 成员
    bad7 = []
    got = {r["item"]: float(r["caught"]) for r in SUMM
           if abs(float(r["k_frac"]) - K_FRAC) < 1e-9}
    for m in L1:
        if m not in got:
            bad7.append("缺 " + m); continue
        if abs(got[m] - exp_map[m]) > 1e-9:
            bad7.append("%s(%s≠%d)" % (m, got[m], exp_map[m]))
    over = [m for m in L1 if got.get(m, 0) > K * PREV]
    L.emit("T45-07b", "PASS" if not over else "WARN",
           "六个 L1 成员无一超过随机期望 %.2f" % (K * PREV),
           "超过者 = %s｜各自命中 = %s" % (over or "无", {m: int(got.get(m, -1)) for m in L1}),
           "增补文件 §3 的稳健性陈述依赖此条；若某成员超过随机期望，该陈述必须改写")
    L.check("T45-07", not bad7, "六个 L1 成员逐个命中数 == 独立重算",
            "不符 %d 项%s" % (len(bad7), ("：" + "; ".join(bad7[:3])) if bad7 else ""),
            "防「只报秩平均、隐藏成员」")

    # ------------------------------------------------------------ T45-08 互相闭合
    # ★ 2026-10-07 新增（外部核验建议的「第 13 项」）。
    #   缺口：逐点独立重算防的是**算错**，防不了**「各自都对、彼此不闭环」**。
    #   实际发生过的缺陷正是这一类：报告把 bootstrap 分布均值写在主对比的位置，
    #   与点估计并排而不标口径 ⇒ 0.4865−0.1892=0.2973 与 +0.2810 看起来矛盾。
    #   CSV 里两列本来都对；错的是**叙述**。故本检查必须读**报告文本**，不能只读 CSV。
    bad8 = []
    mrow = [r for r in MAIN if r["layer"].startswith("L-DR minus")]
    dr_r = [r for r in MAIN if r["layer"] == "L-DR"]
    sh_r = [r for r in MAIN if r["layer"] == "L-shift"]
    if not (mrow and dr_r and sh_r):
        bad8.append("summary_main 缺行")
    else:
        pt_csv = float(mrow[0]["precision_point"])
        bm_csv = float(mrow[0]["precision_boot_mean"])
        # (1) 点估计必须逐位等于两层之差
        if abs(pt_csv - (float(dr_r[0]["precision_point"]) - float(sh_r[0]["precision_point"]))) > 1e-12:
            bad8.append("点估计 ≠ 两层 precision_point 之差")
        # (2)(3)(4) 报告文本必须分列且各自带标签
        rep_path = os.path.join(OUT, "report_p45_triage.txt")
        if not os.path.exists(rep_path):
            bad8.append("报告文本缺失")
        else:
            txt = io.open(rep_path, encoding="utf-8").read()
            def _grab(label):
                for ln in txt.splitlines():
                    if label in ln:
                        mm = re.search(r"([+-]\d+\.\d{4})", ln)
                        if mm:
                            return float(mm.group(1)), ln
                return None, None
            # (a) 点估计：**必须**存在、带标签、且等于 CSV —— 这一条单独就能抓住本次真实缺陷
            pt_txt, pt_ln = _grab("点估计（观测数据")
            if pt_txt is None:
                bad8.append("报告未标「点估计（观测数据」")
            elif abs(pt_txt - pt_csv) > 5e-5:
                bad8.append("报告点估计 %.4f ≠ CSV %.4f" % (pt_txt, pt_csv))
            # (b) bootstrap 均值：**仅当该数值确实出现在文本里**才要求标签。
            #     若作者选择只报点估计（这是推荐做法），不得判红 —— 否则是过度约束
            #     （负向测试第 4 例暴露过这个假阳性）。
            if ("0.2810" in txt) or ("%.4f" % bm_csv) in txt:
                bm_txt, bm_ln = _grab("bootstrap 分布均值")
                if bm_txt is None:
                    bad8.append("文本出现 bootstrap 均值数值，但未标「bootstrap 分布均值」")
                elif abs(bm_txt - bm_csv) > 5e-5:
                    bad8.append("报告 bootstrap 均值 %.4f ≠ CSV %.4f" % (bm_txt, bm_csv))
                if pt_ln is not None and bm_ln is not None and pt_ln == bm_ln:
                    bad8.append("点估计与 bootstrap 均值出现在同一行（口径会被误读）")
    L.check("T45-08", not bad8,
            "互相闭合：点估计 = 两层之差；报告文本把「点估计」与「bootstrap 分布均值」分列且各标各的数",
            "不符 %d 项%s" % (len(bad8), ("：" + "; ".join(bad8[:3])) if bad8 else ""),
            "本项针对的缺陷类型：**每个数单独都对，但并排之后不闭环**")

    # ------------------------------------------------------ T45-09 事后对比自洽
    # ★ 2026-10-12 新增。事后对比 DR − src-only 的数字**要进稿件**（正文与 Table S3），
    #   因此必须落判据。三重检查：
    #     ① 点估计逐位 == 两层 precision_point 之差（**必须**是全样本差，不是 bootstrap 均值）；
    #     ② 区间含点估计与 bootstrap 均值；
    #     ③ 该量在**两个落点**（figures/data/figS10_triage_summary.csv 与
    #        p45_triage_summary_main.csv）逐位一致 —— 防 INC-015 形态：
    #        **同一量两个生产者**，各写一遍同义的代码，一次改动只动一处。
    FIGD = os.path.join(PROJ, "figures", "data")
    bad9 = []
    post_main = [r for r in MAIN if r["layer"] == "L-DR minus L-src"]
    dr_m = [r for r in MAIN if r["layer"] == "L-DR"]
    sr_m = [r for r in MAIN if r["layer"] == "L-src"]
    if not (post_main and dr_m and sr_m):
        bad9.append("summary_main 缺行（L-DR / L-src / 事后对比）")
    else:
        pt9 = float(post_main[0]["precision_point"])
        exp_pt9 = float(dr_m[0]["precision_point"]) - float(sr_m[0]["precision_point"])
        if abs(pt9 - exp_pt9) > 1e-12:
            bad9.append("点估计 %.9g ≠ 两层之差 %.9g" % (pt9, exp_pt9))
        lo9, hi9 = float(post_main[0]["lo95"]), float(post_main[0]["hi95"])
        bm9 = float(post_main[0]["precision_boot_mean"])
        if not (lo9 < bm9 < hi9):
            bad9.append("bootstrap 区间不含分布均值")
        if not (lo9 <= pt9 <= hi9):
            bad9.append("点估计落在区间外（%.6f not in [%.6f, %.6f]）" % (pt9, lo9, hi9))
        if not (-1.0 <= lo9 <= hi9 <= 1.0):
            bad9.append("差值区间越出 [-1,1]（[%.6f, %.6f]）" % (lo9, hi9))

        # ③ 两个落点逐位一致（layer 拼写不同：下划线 vs 空格）
        p_fig = os.path.join(FIGD, "figS10_triage_summary.csv")
        if not os.path.exists(p_fig):
            bad9.append("figS10_triage_summary.csv 缺失")
        else:
            with io.open(p_fig, encoding="utf-8") as f:
                fig_rows = [r for r in csv.DictReader(f) if r["layer"] == "L-DR_minus_L-src"]
            if not fig_rows:
                bad9.append("figS10 缺 L-DR_minus_L-src 行")
            else:
                for col in ("precision_point", "precision_boot_mean", "lo95", "hi95"):
                    a, b = float(fig_rows[0][col]), float(post_main[0][col])
                    if abs(a - b) > 1e-12:
                        bad9.append("落点不一致 %s（figS10 %.9g vs main %.9g）" % (col, a, b))
    L.check("T45-09", not bad9,
            "事后对比 DR − src-only 自洽（点估计 = 全样本差；区间含点估计；两个落点逐位一致）",
            "不符 %d 项%s" % (len(bad9), ("：" + "; ".join(bad9[:3])) if bad9 else ""),
            "该数字要进正文与 Table S3；同一量在两个文件各有一个落点，须防漂移")

    # 事后对比的区间是否跨 0 —— 决定正文允不允许作显著性声称
    if post_main:
        lo9, hi9 = float(post_main[0]["lo95"]), float(post_main[0]["hi95"])
        L.emit("T45-09b", "PASS",
               "事后对比的 95% 区间跨 0 时，正文**不得**声称该对比显著（本条只报事实）",
               "区间 = [%+.4f, %+.4f] ⇒ %s" % (lo9, hi9, "跨 0" if lo9 < 0 < hi9 else "不跨 0"),
               "事后对比 DR − src-only；主对比仍是预注册的 DR − L-shift（T45-05b 管）")
