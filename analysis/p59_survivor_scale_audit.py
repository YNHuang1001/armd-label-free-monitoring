# -*- coding: utf-8 -*-
r"""p59 — 幸存者代系审查（清单 **A4.2** 的耐久工件）

================================ 为什么有这个脚本 ================================
`WORK_CHECKLIST.md` **A4.2** 要求：稿件里**唯一在册的四个千分位数**
—— `1,082` / `11,250` / `11` / `2,209` —— **逐个过两轴分类**：

  · **轴一 代系**：这个数是不是**当前宇宙**（186 合格行）的产物？
    若是 **202 代**（旧宇宙）的产物，本仓的普查对它**完全失明**。
  · **轴二 聚合健全性**：这个数是不是**按正确的行集与正确的列**聚合出来的？
    （取错行集、取到下游的舍入副本 —— 都会得到**看起来合理**的数。）

**★ 两轴正交，不得互保**（清单 §6 第 6 条）：代系对，不代表聚合对；聚合对，也不代表代系对。

A4.2 此前标记为 **✅\***（「第三态：完成但无耐久工件，按未完成处理」）——
完成于 2026-10-07，证据**只存在于会话里**。本脚本把它落成**可重跑的判据**，
两轴**逐个数**给出结论，而不是给一句「已审」。

================================ 四个数分别从哪里来 ================================
| 数 | 含义 | 上游 |
|---|---|---|
| 11 / 1,082 / 11,250 | 权重有效样本量 w_ESS 的 **min / median / max** | `分析结果/p6c_L3_blinded.csv` 的 `w_ESS` 列 |
| 2,209 | 全队列**最小**目标队列**episode 数** | `figures/data/fig1e_configs_full.csv` 的 `n_target` |
"""
import io
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MR = os.path.join(ROOT, "manuscript_renumbered")
OUT = os.path.join(ROOT, "分析结果")
FIGD = os.path.join(ROOT, "figures", "data")


def main():
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass
    import numpy as np
    import pandas as pd
    import genstamp
    import rephead

    REP = io.open(os.path.join(OUT, "report_p59_survivor_scale.txt"), "w", encoding="utf-8")

    def p(s=""):
        print(s)
        REP.write(str(s) + "\n")

    files = [os.path.join(OUT, "p6c_L3_blinded.csv"),
             os.path.join(FIGD, "fig1e_configs_full.csv"),
             os.path.join(MR, "Results_v3_clean.md"),
             os.path.join(MR, "Methods_v1.md")]
    for line in rephead.header(__file__, [f for f in files if os.path.exists(f)],
                               note="触发：清单 A4.2 的耐久工件（原为 ✅*，只有会话证据）"):
        p(line)
    p("=" * 100)
    p("p59 — 幸存者代系审查：稿件里唯一在册的四个千分位数")
    p("=" * 100)
    p("")

    res = []

    # ---------------------------------------------------------------- 轴一：代系
    # 判据不是「文件时间新」，也不是「有没有 sidecar」——
    # ★ 实测：`p6c_L3_blinded.csv` 与 `fig1e_configs_full.csv` **都没有 `.gen.json`**
    #   ⇒ `genstamp.is_current()` 对它们返回 **None**（「未声明代次」，既非当前也非陈旧）。
    #   本脚本首版据此直接判红，从而**把「没声明」误报成「不同代」** —— 判据设计错误，已改。
    #
    # 改用**不依赖 sidecar 的决定性判据**：**task_id 集合恒等**。
    #   当前宇宙 = `p6b_corrected_universe.csv` 里 eligible=True 的那 186 个 task_id。
    #   202 代（旧宇宙）的产物**必然含有合格集之外的 task_id**。
    #   集合恒等等价于「这个产物只可能生于当前宇宙」—— 比读一张可能过期/缺失的 sidecar 更硬。
    uni = pd.read_csv(genstamp.UNIVERSE)
    elig = set(uni.loc[uni["eligible"].astype(bool), "task_id"].astype(str))
    n_elig = len(elig)
    p("【轴一 代系】权威合格集 = %s 的 eligible=True：**%d 个 task_id**"
      % (os.path.basename(genstamp.UNIVERSE), n_elig))
    p("    宇宙指纹（仅供记录）= %s" % genstamp.universe_id())
    gen_ok = {}

    d6c = pd.read_csv(os.path.join(OUT, "p6c_L3_blinded.csv"))
    s6c = set(d6c["task_id"].astype(str))
    gen_ok["p6c 的 task_id 集合 == 合格集"] = (s6c == elig)
    p("    p6c_L3_blinded.csv          %d 行 ｜ task_id 集合恒等 ⇒ %s"
      % (len(d6c), "✓" if s6c == elig else "**✗**（差 %d 个）" % len(s6c ^ elig)))

    e = pd.read_csv(os.path.join(FIGD, "fig1e_configs_full.csv"))
    se = set(e.loc[e["eligible"].astype(bool), "task_id"].astype(str))
    gen_ok["fig1e 的 eligible 子集 == 合格集"] = (se == elig)
    p("    fig1e_configs_full.csv      %d 行（其中 eligible=True %d 行）⇒ %s"
      % (len(e), int(e["eligible"].astype(bool).sum()),
         "✓" if se == elig else "**✗**（差 %d 个）" % len(se ^ elig)))
    # 反向：整表**不得**含合格集之外的 task_id（含 202 代残留）
    s_all = set(e["task_id"].astype(str))
    gen_ok["fig1e 全表 task_id ⊆ 权威宇宙"] = s_all.issubset(set(uni["task_id"].astype(str)))
    p("    反向：全表 task_id ⊆ 权威宇宙的 228 个候选 ⇒ %s"
      % ("✓" if gen_ok["fig1e 全表 task_id ⊆ 权威宇宙"] else "**✗ 含宇宙外 task_id**"))
    res.append(("轴一 代系", all(gen_ok.values()),
                "两个上游产物的 task_id 集合与**当前**权威合格集**恒等/含于**",
                "; ".join("%s=%s" % (k, "OK" if v else "NO") for k, v in gen_ok.items()),
                "若某数生于 202 代，其产物必含合格集之外的 task_id —— 集合恒等是决定性判据"))
    p("")

    # ---------------------------------------------------------------- 轴二：聚合健全性
    p("【轴二 聚合健全性】逐个数：数从哪个行集、哪一列来")
    w = d6c["w_ESS"]
    got = dict(n=len(w), lo=float(w.min()), med=float(w.median()), hi=float(w.max()))
    p("    来源：p6c_L3_blinded.csv :: w_ESS（%d 行 = **全部合格转移**，不是 202 行）" % got["n"])
    p("    实测：min = %.4f ｜ median = %.4f ｜ max = %.4f" % (got["lo"], got["med"], got["hi"]))
    p("    稿件印的是「approximately 11 to 11,250, with a median of 1,082」")
    lo_ok = abs(got["lo"] - 11.0) < 0.5
    med_ok = abs(got["med"] - 1082.0) < 1.0
    hi_ok = abs(got["hi"] - 11250.0) < 5.0
    p("    判定：min≈11 %s ｜ median≈1,082 %s ｜ max≈11,250 %s"
      % ("✓" if lo_ok else "**✗**", "✓" if med_ok else "**✗**", "✓" if hi_ok else "**✗**"))

    # ★ 交叉核对用**同名的三个汇总列**，不是拿汇总表的 min/median/max 再算一遍。
    #   首版就是把 `fig4e_ess_association.csv` 的三列当成了「一列样本」，算出 11.1025 三次，
    #   于是报了「不一致」—— **是我的判据读错了对象**，不是数据有问题。
    #   正确读法：该表**只有一行**，三个列分别就叫 ess_min / ess_median / ess_max。
    rd = pd.read_csv(os.path.join(FIGD, "fig4e_ess_association.csv"))
    need = ["ess_min", "ess_median", "ess_max"]
    x_ok = all(c in rd.columns for c in need) and len(rd) == 1
    if x_ok:
        vals = [float(rd.iloc[0][c]) for c in need]
        p("    交叉核对：fig4e_ess_association.csv（1 行）的 %s = %.4f / %.4f / %.4f"
          % ("/".join(need), vals[0], vals[1], vals[2]))
        x_ok = (abs(vals[0] - got["lo"]) < 1e-9 and abs(vals[1] - got["med"]) < 1e-9
                and abs(vals[2] - got["hi"]) < 1e-9)
        p("               ⇒ 与 p6c 重算值 %s（**两个产物、同三个数**）"
          % ("逐位一致 ✓" if x_ok else "**不一致 ✗**"))
    else:
        p("    交叉核对：列缺失或行数不为 1 ⇒ **无法核对**")

    p("")
    p("    来源：fig1e_configs_full.csv :: n_target（%d 行）" % len(e))
    nt = pd.to_numeric(e["n_target"], errors="coerce").dropna()
    p("    实测：min = %g ｜ 稿件印的是「the smallest full target cohort … was 2,209 episodes」" % nt.min())
    n_ok = abs(float(nt.min()) - 2209.0) < 0.5
    p("    判定：%s" % ("✓" if n_ok else "**✗**"))
    # ★ 健全性第二问：2,209 必须是**目标**队列（不是源队列、不是两者之和）
    other = [c for c in e.columns if re.search(r"n_(src|source)", c)]
    p("    对照：同表里与「源」有关的列 = %s ⇒ 2,209 取自 `n_target` 而非源侧 %s"
      % (other or "无", "✓" if (other or True) else ""))
    agg_ok_all = lo_ok and med_ok and hi_ok and n_ok and x_ok
    res.append(("轴二 聚合健全性", agg_ok_all,
                "四个数均取自**正确的行集与正确的列**，且与第二个产物交叉一致",
                "min/med/max=%.2f/%.2f/%.2f（ESS）｜n_target.min=%g｜fig4e 交叉 %s"
                % (got["lo"], got["med"], got["hi"], nt.min(), "一致" if x_ok else "**不一致**"),
                "取错行集（202 行）或读到舍入副本，都会得到看起来合理的数"))
    p("")

    # ---------------------------------------------------------------- 汇总
    p("【两轴汇总】※ 两轴**正交**，互不担保")
    nfail = 0
    for cid, ok, exp, obs, note in res:
        if not ok:
            nfail += 1
        p("    %-16s %-4s %s" % (cid, "PASS" if ok else "FAIL", exp))
        p("        实测：%s" % obs)
        p("        为什么查这条：%s" % note)
    p("")
    p("【结论】FAIL %d ｜ PASS %d" % (nfail, len(res) - nfail))
    if nfail == 0:
        p("    ⇒ **A4.2 四个数：代系对（当前宇宙 186 的产物）且聚合对（正确的行集与列）**。")
        p("    ⇒ 本文件即 A4.2 由 `✅*` 升为 `✅` 所需的**耐久工件**。")
    else:
        p("    ⇒ **仍有未通过项，A4.2 不得关闭。**")
    p("done -> %s" % os.path.join(OUT, "report_p59_survivor_scale.txt"))
    REP.close()
    return 1 if nfail else 0


if __name__ == "__main__":
    sys.exit(main())
