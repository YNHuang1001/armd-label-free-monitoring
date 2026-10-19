# -*- coding: utf-8 -*-
r"""p61 — 转移级 vs 配对层面的分解（**取向注册表** ＋ 已知答案测试 ＋ 区间）

================================ 为什么有这个脚本 ================================
外部核验（两轮）要求把三件事落成**入库的单一生产者**，而不是停在诊断脚本里：

  ① **取向注册表**：本会话连续 **四次** 在「风险取向没对齐」上出错
     （−Δ̂ 的配对内一致率、镜像的朝向、对照组写成 −(src−mean)、③ 表里 src 的朝向）。
     根因是**各脚本各自写负号**。⇒ 一切风险分**只能从本表的 `RISK` 取**，
     任何脚本里不得再出现临时的 `-x`。
     ★ 光有恒等式 `一致率(s)+一致率(−s)=1` 没用 —— 它是**恒真**的，
       抓不到「对照组写成 −(src−mean)」这一类的**组合错误**。
       真正管用的是**已知答案测试**：钉在 p57 的冻结值上，不符即红。

  ② **分解要带区间**：0.6382 / 0.5535 / 0.504 / 0.3761 此前全是点估计。

  ③ **配对内的一致率（0.785，连续 Δ 的排序）与去均值分量的池化 AUROC（0.5535，二分类）
     是两个量**，读者会问为什么一强一弱。⇒ 补一个**直接的**配对内 AUROC：
     只取**两个目标 harmful 标签恰好不同**的配对，看风险分能否把 harmful 那个排在前面。

  ④ **n 档的覆盖面与可用性**：`有限 ≠ 可用`。Â 在 n=100 能到 −86。
     并把 **同子集的 full 档**一并算出，否则「0.645 → 0.384」混合了覆盖面丢失与区分度下降。

================================ 风险取向（本文件唯一的口径来源）========================
    **分数越高 ⇒ 越可能 harmful（观测 ΔAUROC 越低）**
    `src_only` 的风险取向是 **+src**（p45 的 L-src 即此方向），其 harmful-AUROC = 0.6107。
"""
import io
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "分析结果")
B = 10000
SEED = 20261007

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

import numpy as np
import pandas as pd
import genstamp
import rephead


# ================================================================ 取向注册表
def _registry(frame, mean_hat):
    """**唯一**定义风险分的地方。键 = 正文里用的名字，值 = (风险分向量, 取向说明, p57 冻结值)。"""
    src = frame["auroc_source_internal"].values.astype(float)
    A = frame["AUROC_DR"].values.astype(float)
    mm = mean_hat
    return {
        "prespecified_DR":   (-A,            "风险分 = −Â_target（预先指定）",            0.6436828505371956),
        "src_only":          (src,           "风险分 = +src（越大越先复核）",              0.6107031499651972),
        "amean_control":     (src - mm,      "风险分 = +(src − mean(Â|配对))",           0.7880458442916880),
        "posthoc_dhat_DR":   (src - A,       "风险分 = +(src − Â) = −Δ̂（事后对齐）",      0.7886582580926395),
        "pair_mean_part":    (-mm,           "分解：−mean(Â|配对)（配对层面）",            None),
        "within_pair_dev":   (-(A - mm),     "分解：−(Â − mean)（配对内偏差）",             None),
    }


def wmw(ys, ss):
    """Mann–Whitney AUROC，并列取平均秩。**取向由调用方负责。**"""
    ys = np.asarray(ys, float); ss = np.asarray(ss, float)
    m = np.isfinite(ss); ys, ss = ys[m], ss[m]
    if len(ys) == 0 or len(set(ys)) < 2:
        return float("nan")
    p = sorted(zip(ss, ys)); n1 = int(sum(ys)); n0 = len(ys) - n1
    if n1 == 0 or n0 == 0:
        return float("nan")
    rk = [0.0] * len(p); i = 0
    while i < len(p):
        j = i
        while j + 1 < len(p) and p[j + 1][0] == p[i][0]:
            j += 1
        a = 0.5 * (i + j) + 1
        for k in range(i, j + 1):
            rk[k] = a
        i = j + 1
    return (sum(rk[k] for k in range(len(p)) if p[k][1] == 1) - n1 * (n1 + 1) / 2) / (n1 * n0)


def main():
    REP = io.open(os.path.join(OUT, "report_p61_transfer_level.txt"), "w", encoding="utf-8")

    def p(s=""):
        print(s)
        REP.write(str(s) + "\n")

    D = pd.read_csv(os.path.join(OUT, "p6d_outcome.csv"))
    # ★ 整表读入：⑥ 还要用 AUROC_IOW / AUROC_OM / SUDO_AURCC（首版只留了 task_id+AUROC_DR，撞了 KeyError）
    CFULL = pd.read_csv(os.path.join(OUT, "p6c_L3_blinded.csv"))
    S = D.merge(CFULL[["task_id", "AUROC_DR"]], on="task_id", validate="1:1") \
         .sort_values("task_id").reset_index(drop=True)
    S["harm"] = (S["d_auroc"] <= -0.05).astype(int)
    S["cluster"] = S["organism"] + "|" + S["drug"]
    S["pair"] = S["organism"] + "|" + S["drug"] + "|" + S["source_site"]
    mean_hat = S.groupby("pair")["AUROC_DR"].transform("mean").values
    RISK = _registry(S, mean_hat)
    y = S["harm"].values

    for line in rephead.header(__file__, [os.path.join(OUT, "p6d_outcome.csv"),
                                          os.path.join(OUT, "p6c_L3_blinded.csv"),
                                          os.path.join(OUT, "p10b_raw.csv")],
                               note="触发：外部核验要求把转移级/配对层面的分解落成入库生产者"):
        p(line)
    p("=" * 100)
    p("p61 — 转移级 vs 配对层面：取向注册表 ＋ 已知答案测试 ＋ 带区间的分解")
    p("=" * 100)
    p("风险取向（本文件唯一定义）：**分数越高 ⇒ 越可能 harmful**")
    for k, (_, note, exp) in RISK.items():
        p("    %-18s %s%s" % (k, note, "" if exp is None else "   冻结值 %.4f" % exp))
    p("")

    # ---------------------------------------------------------------- ① 已知答案测试
    p("【① 已知答案测试】钉在 p57 的冻结值上；任一不符即判红（这是抓「组合写错」的唯一手段）")
    nfail = 0
    for k, (sc, note, exp) in RISK.items():
        if exp is None:
            continue
        got = wmw(y, sc)
        ok = abs(got - exp) < 5e-5
        if not ok:
            nfail += 1
        p("    %-18s 实测 %.6f ｜ 冻结 %.6f ｜ %s" % (k, got, exp, "OK" if ok else "**FAIL**"))
    p("")

    # ---------------------------------------------------------------- ② 分解 + 区间
    cl = sorted(S["cluster"].unique())
    groups = {c: np.where(S["cluster"].values == c)[0] for c in cl}
    rng = np.random.default_rng(SEED)

    def boot(fn, B=B, seed=SEED):
        # ★★ 2026-10-13 修：**每次调用用独立的生成器**。
        #   原先 `rng` 在 `boot` 外只建一次、被所有调用**共享消费** ⇒ 只要给 p61 增删任何一处
        #   `boot()` 调用，其后**所有**区间都会漂移。实测：`AUROC_DR − src_only` 的区间
        #   在两次编辑之间由 `[−0.1457, +0.1937]` 变成 `[−0.1400, +0.1937]`。
        #   一个「独立复算」若会因无关编辑而改变自己的输出，就不能当稳定的对照。
        rng = np.random.default_rng(seed)
        out = []
        for _ in range(B):
            pick = rng.choice(len(cl), len(cl), replace=True)
            idx = np.concatenate([groups[cl[i]] for i in pick])
            out.append(fn(idx))
        return np.array([v for v in out if np.isfinite(v)])

    p("【② 分解（全 186），带 37 簇 bootstrap 区间】")
    rows = []
    for k in ("prespecified_DR", "pair_mean_part", "within_pair_dev", "src_only",
              "amean_control", "posthoc_dhat_DR"):
        sc = RISK[k][0]
        pt = wmw(y, sc)
        bs = boot(lambda idx, _sc=sc: wmw(y[idx], _sc[idx]))
        lo, hi = np.percentile(bs, [2.5, 97.5])
        rows.append(dict(score=k, point=pt, lo95=float(lo), hi95=float(hi)))
        p("    %-18s %-22s %.4f  [%.4f, %.4f]"
          % (k, RISK[k][1][:20], pt, lo, hi))
    pd.DataFrame(rows).to_csv(os.path.join(OUT, "p61_decomposition.csv"), index=False)

    # 差：完整分数 − 仅用配对均值的分数
    scA = RISK["prespecified_DR"][0]
    scM = RISK["pair_mean_part"][0]
    pt = wmw(y, scA) - wmw(y, scM)
    bs = boot(lambda idx: wmw(y[idx], scA[idx]) - wmw(y[idx], scM[idx]))
    lo, hi = np.percentile(bs, [2.5, 97.5])
    p("    %-18s %.4f  [%.4f, %.4f]  %s"
      % ("完整 − 配对均值", pt, lo, hi, "含 0" if lo < 0 < hi else "不含 0"))

    # ★★ 2026-10-12：**事后对齐分 − 配对层面控制** —— 全文最关键的那一对（p57 已补）。
    #   缺了它，正文只能拿 −Δ̂ 与 src_only 比（+0.178，不含 0），读起来像「事后分确实更强」。
    scH = RISK["posthoc_dhat_DR"][0]
    scC = RISK["amean_control"][0]
    pt2 = wmw(y, scH) - wmw(y, scC)
    bs2 = boot(lambda idx: wmw(y[idx], scH[idx]) - wmw(y[idx], scC[idx]))
    lo2, hi2 = np.percentile(bs2, [2.5, 97.5])
    p("    %-18s %.4f  [%.4f, %.4f]  %s   ← **决定 −Δ̂ 是否超出配对层面的那一对**"
      % ("−Δ̂ − Â均值对照", pt2, lo2, hi2, "含 0" if lo2 < 0 < hi2 else "不含 0"))
    # 已知答案：p57 冻结值（点估计逐位可复现；区间受 RNG 流影响，仅比点估计）
    KA = {"完整 − 配对均值": 0.0055, "−Δ̂ − Â均值对照": 0.0006}
    for nm, got, exp in (("完整 − 配对均值", pt, KA["完整 − 配对均值"]),
                         ("−Δ̂ − Â均值对照", pt2, KA["−Δ̂ − Â均值对照"])):
        okk = abs(got - exp) < 5e-5
        if not okk:
            nfail += 1
        p("       已知答案 %-16s 实测 %+.6f ｜ p57 冻结 %+.6f ｜ %s"
          % (nm, got, exp, "OK" if okk else "**FAIL**"))
    p("    ★ 该「均值」用了**两个目标**的协变量 ⇒ 只能说**配对层面**，不得说 source 模型层面")
    p("")

    # ---------------------------------------------------------------- ③ 配对内直接 AUROC
    p("【③ 配对内直接 AUROC】只取两个目标 harmful 标签**恰好不同**的配对")
    discord = []
    for _, gg in S.groupby("pair"):
        if len(gg) != 2:
            continue
        a, b = gg.iloc[0], gg.iloc[1]
        if a["harm"] == b["harm"]:
            continue
        hi_risk = a if (a["harm"] == 1) else b        # 应当被排在前面的那个
        lo_risk = b if hi_risk is a else a
        discord.append((hi_risk.name, lo_risk.name, hi_risk["cluster"]))
    p("    不一致对数量 = %d（占 %d 个配对）" % (len(discord), S["pair"].nunique()))
    hits = 0
    for i, j, _ in discord:
        if RISK["prespecified_DR"][0][i] > RISK["prespecified_DR"][0][j]:
            hits += 1
    p("    Â（−Â 风险取向）把 harmful 那个排前面的比例 = %.4f  (%d/%d)"
      % (hits / len(discord), hits, len(discord)))
    hits2 = sum(1 for i, j, _ in discord if RISK["posthoc_dhat_DR"][0][i] > RISK["posthoc_dhat_DR"][0][j])
    p("    −Δ̂ 同上                                     = %.4f  (%d/%d)"
      % (hits2 / len(discord), hits2, len(discord)))
    # 簇 bootstrap
    if discord:
        arr = np.array([(RISK["prespecified_DR"][0][i] > RISK["prespecified_DR"][0][j])
                        for i, j, _ in discord], dtype=float)
        cls = np.array([c for _, _, c in discord])
        ucl = sorted(set(cls))
        bs = []
        r2 = np.random.default_rng(SEED)
        for _ in range(B):
            pick = r2.choice(len(ucl), len(ucl), replace=True)
            sel = np.concatenate([np.where(cls == ucl[i])[0] for i in pick])
            if len(sel):
                bs.append(arr[sel].mean())
        bs = np.array(bs)
        p("    其 95%% 簇 bootstrap 区间 = [%.4f, %.4f]%s"
          % (np.percentile(bs, 2.5), np.percentile(bs, 97.5),
             "  ⇒ 含 0.5" if np.percentile(bs, 2.5) < 0.5 < np.percentile(bs, 97.5) else "  ⇒ 不含 0.5"))
    p("")

    # ---------------------------------------------------------------- ④ n 档覆盖面与可用性
    p("【④ n 档：有限 ≠ 可用】Â 超出 [0,1] 即不可用；越界的**不得**计入「可用」")
    raw = pd.read_csv(os.path.join(OUT, "p10b_raw.csv"))
    cov = []
    for lv in ["100", "250", "500", "1000", "full"]:
        d = raw[raw.level == lv]
        fin = np.isfinite(d["AUROC_DR"].values)
        inr = fin & (d["AUROC_DR"].values >= 0) & (d["AUROC_DR"].values <= 1)
        per_rep = d.assign(fin=fin, inr=inr).groupby("rep")[["fin", "inr"]].sum()
        cov.append(dict(level=lv, cells=len(d), finite_frac=float(fin.mean()),
                        inrange_frac=float(inr.mean()),
                        finite_per_rep=float(per_rep["fin"].mean()),
                        inrange_per_rep=float(per_rep["inr"].mean()),
                        oor=int((fin & ~inr).sum())))
        p("    level=%-5s 有限 %.3f ｜ 落在[0,1]内 %.3f ｜ 有限但越界 %6d 格 ｜ 每重复中位：有限 %.1f / 可用 %.1f（共 186）"
          % (lv, fin.mean(), inr.mean(), int((fin & ~inr).sum()),
             per_rep["fin"].median(), per_rep["inr"].median()))
    pd.DataFrame(cov).to_csv(os.path.join(OUT, "p61_coverage_by_level.csv"), index=False)
    p("")
    p("【⑤ 同子集的 full 档】每重复用**该重复在 n=100 有限的那批转移**，同批上算 full 档")
    d100 = raw[raw.level == "100"]
    dfull = raw[raw.level == "full"]
    srcmap = S.set_index("task_id")["auroc_source_internal"]
    hm = S.set_index("task_id")["harm"]
    diffs = []
    for rep, g100 in d100.groupby("rep"):
        m = np.isfinite(g100["AUROC_DR"].values)
        ids = g100["task_id"].values[m]
        gf = dfull[(dfull.rep == rep) & (dfull.task_id.isin(ids))]
        if len(gf) < 30:
            continue
        yy = hm.loc[gf["task_id"].values].values.astype(int)
        a100 = wmw(yy, -g100.set_index("task_id").loc[gf["task_id"].values, "AUROC_DR"].values)
        afull = wmw(yy, -gf["AUROC_DR"].values)
        diffs.append(afull - a100)
    diffs = np.array([d for d in diffs if np.isfinite(d)])
    p("    掩码=「n=100 有限」        full − n=100 = %+.4f  (sd %.4f, n=%d 重复)"
      % (diffs.mean(), diffs.std(), len(diffs)))
    # ★ 更好的掩码是「**可用**」（在 [0,1] 内）—— Â 能到 −86，「有限」远不等于「可用」。
    diffs2, ns2 = [], []
    for rep, g100 in d100.groupby("rep"):
        v = g100["AUROC_DR"].values
        m = np.isfinite(v) & (v >= 0) & (v <= 1)
        if m.sum() < 20:
            continue
        ids = g100["task_id"].values[m]
        gf = dfull[(dfull.rep == rep) & (dfull.task_id.isin(ids))]
        if len(gf) < 20:
            continue
        yy = hm.loc[gf["task_id"].values].values.astype(int)
        ser = pd.Series(v[m], index=g100["task_id"].values[m])
        a100 = wmw(yy, -ser.loc[gf["task_id"].values].values)
        afull = wmw(yy, -gf["AUROC_DR"].values)
        if np.isfinite(a100) and np.isfinite(afull):
            diffs2.append(afull - a100); ns2.append(int(m.sum()))
    diffs2 = np.array(diffs2)
    p("    掩码=「n=100 **可用**(∈[0,1])」 full − n=100 = %+.4f  (sd %.4f, n=%d 重复, 子集中位 %d 个转移)"
      % (diffs2.mean(), diffs2.std(), len(diffs2), int(np.median(ns2))))
    p("    ⇒ 朴素的「0.645 − 0.384 = +0.261」把**覆盖面丢失**与**区分度下降**混在了一起；")
    p("      上面两行是**同批转移**上的配对差，才是可比的那个数。")
    p("")

    # ---------------------------------------------------------------- ⑥ 15 个预先指定信号 vs src-only
    # ★ 外部核验（2026-10-12）：题目若要写「No **prespecified** label-free signal …」，
    #   就必须对**全部 15 个**预先指定的信号都给出 (信号 − src_only) 的配对差异，
    #   不能只给 DR。**任一信号的区间不含 0，题目就得改。**
    p("")
    p("【⑥ 15 个预先指定信号 vs src-only 的配对差异】（决定题目能否用 prespecified 这一支）")
    p("    ★ **未做多重比较校正** —— 15 个区间同时看，属探索性并列表述，不据此下单一结论。")
    B15 = pd.read_csv(os.path.join(OUT, "p6b_label_free_signals_blinded.csv"))
    L1 = ["domain_auc", "mmd", "psi_mean", "psi_max", "js_mean", "js_max"]
    SG = ["pred_mean_shift", "entropy_mean_shift", "pred_pos_prev_shift"]
    OTHER = ["entropy_js", "conf_dist_js", "SUDO_AURCC"]
    L3 = ["AUROC_DR", "AUROC_IOW", "AUROC_OM"]
    b = B15.set_index("task_id")
    c = CFULL.set_index("task_id")
    SIG = {}
    # 出处逐列核过：L1 / SG / entropy_js / conf_dist_js 在 p6b；SUDO_AURCC 与 L3 在 p6c。
    for m in L1 + ["entropy_js", "conf_dist_js"]:
        SIG[m] = b.loc[S["task_id"].values, m].values.astype(float)
    for m in SG:
        SIG[m] = np.abs(b.loc[S["task_id"].values, m].values.astype(float))
    for m in L3 + ["SUDO_AURCC"]:
        SIG[m] = -c.loc[S["task_id"].values, m].values.astype(float)
    frozen = {r["method"]: float(r["harmful_auc"])
              for r in pd.read_csv(os.path.join(OUT, "p6d_signal_eval.csv")).to_dict("records")}

    p("    %-22s %-10s %-22s %s" % ("信号", "AUROC", "对 src-only 的配对差", "冻结值核对"))
    sc_src = RISK["src_only"][0]
    n_excl = 0
    tab = []
    for m in L1 + SG + OTHER + L3:
        pt = wmw(y, SIG[m])
        bs = boot(lambda idx, _s=SIG[m]: wmw(y[idx], _s[idx]) - wmw(y[idx], sc_src[idx]))
        lo, hi = np.percentile(bs, [2.5, 97.5])
        excl = not (lo < 0 < hi)
        n_excl += int(excl)
        ok = abs(pt - frozen[m]) < 5e-5
        if not ok:
            nfail += 1
        tab.append(dict(signal=m, auroc=pt, d_vs_src_only=pt - wmw(y, sc_src),
                        lo95=float(lo), hi95=float(hi), excludes_zero=bool(excl)))
        p("    %-22s %.4f     %+.4f [%+.4f, %+.4f]  %-4s  %s"
          % (m, pt, pt - wmw(y, sc_src), lo, hi, "不含0" if excl else "含0",
             "OK" if ok else "**与冻结值不符**"))
    pd.DataFrame(tab).to_csv(os.path.join(OUT, "p61_signals_vs_srconly.csv"), index=False)
    p("")
    p("    src-only 本身 = %.4f [%.4f, %.4f]"
      % (wmw(y, sc_src), *np.percentile(boot(lambda idx: wmw(y[idx], sc_src[idx])), [2.5, 97.5])))
    p("    15 个信号中，对 src-only 的配对差**区间不含 0** 者：%d 个" % n_excl)
    if n_excl:
        p("    ⇒ **题目不得声称「全部预先指定信号都未显示优于 src-only」** —— 见上表。")
    else:
        p("    ⇒ 没有任何预先指定信号被证明优于 src-only（**不等于**「无增益」：")
        p("      每个区间给出了各自能排除的上界，须逐个照实写）。")

    p("")
    p("【汇总】已知答案测试 FAIL %d" % nfail)
    p("done -> %s" % os.path.join(OUT, "report_p61_transfer_level.txt"))
    REP.close()
    return 1 if nfail else 0


if __name__ == "__main__":
    sys.exit(main())
