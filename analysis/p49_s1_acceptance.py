# -*- coding: utf-8 -*-
r"""p49 — Table S1 的可复现性验收（A2.3）

================================ 验收条款（外部核验给定，完整版）================
**从 S1 的打印列出发**（不看任何其它冻结产物），必须逐位复现：
  (a) composite（六个 L1 的等权秩平均）在 k=37 的命中数 = **7**
  (b) **S2 Panel A 那五行**（此前缺列的五个信号）的 harmful-detection AUROC

若任一不成立 ⇒ **停线**，不得进文字轮。

================================ 为什么这层必须存在 ================================
S1 的首句声称「All derived quantities in this manuscript are computed on these rows」。
在本轮补齐 5 个信号列之前，**composite 的 7/37 无法被任何人从发布的材料复算** ——
外部核验实测报出。本脚本就是那句首句的**机械担保**。

================================ 口径（一律取自 Methods §11 / 冻结产物，不新发明）========
方向变换：+s 用于六个 L1 与 entropy_js / conf_dist_js；
          |s| 用于 pred_mean_shift / entropy_mean_shift / pred_pos_prev_shift；
          −s 用于 SUDO_AURCC 与三个 L3。
并列按 task_id 升序打破。
"""
import csv
import io
import os
import sys

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIG = os.path.join(ROOT, "figures", "data")
OUT = os.path.join(ROOT, "分析结果")
REP = open(os.path.join(OUT, "report_p49_s1_acceptance.txt"), "w", encoding="utf-8")


def p(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True); REP.write(s + "\n"); REP.flush()


L1 = ["domain_auc", "mmd", "psi_mean", "psi_max", "js_mean", "js_max"]
SIGNED = ["pred_mean_shift", "entropy_mean_shift", "pred_pos_prev_shift"]
L3 = ["AUROC_DR", "AUROC_IOW", "AUROC_OM"]
NEG = ["SUDO_AURCC"] + L3

S1 = list(csv.DictReader(io.open(os.path.join(FIG, "tableS1_tasks.csv"), encoding="utf-8")))
T = [r["task_id"] for r in S1]
p("=" * 104)
p("p49 — Table S1 可复现性验收（A2.3）")
p("=" * 104)
p(f"S1 行数 {len(S1)} ｜ 列数 {len(S1[0])} ｜ **仅用 S1 自身的列**，不读其它冻结产物")


def risk(m):
    v = [float(r[m]) for r in S1]
    return [-x for x in v] if m in NEG else ([abs(x) for x in v] if m in SIGNED else v)


def avg_rank(v):
    o = sorted(range(len(v)), key=lambda i: v[i])
    r = [0.0] * len(v); i = 0
    while i < len(o):
        j = i
        while j + 1 < len(o) and v[o[j + 1]] == v[o[i]]:
            j += 1
        a = 0.5 * (i + j) + 1
        for k in range(i, j + 1):
            r[o[k]] = a
        i = j + 1
    return r


def top(s, k):
    return sorted(range(len(T)), key=lambda i: (-s[i], T[i]))[:k]


def wmw(y, s):
    pr = sorted(zip(s, y)); n1 = sum(y); n0 = len(y) - n1
    r = [0.0] * len(pr); i = 0
    while i < len(pr):
        j = i
        while j + 1 < len(pr) and pr[j + 1][0] == pr[i][0]:
            j += 1
        a = 0.5 * (i + j) + 1
        for k in range(i, j + 1):
            r[k] = a
        i = j + 1
    return (sum(r[k] for k in range(len(pr)) if pr[k][1] == 1) - n1 * (n1 + 1) / 2) / (n1 * n0)


y = [int(r["harmful"]) for r in S1]
NF = sum(y)
p(f"失效 {NF} / {len(y)}")

# ---------------------------------------------------------------- (a) composite @ k=37
R6 = [avg_rank([float(r[m]) for r in S1]) for m in L1]
comp = [sum(R6[j][i] for j in range(6)) / 6.0 for i in range(len(S1))]
caught = sum(y[i] for i in top(comp, 37))
p("")
p("(a) composite（六个 L1 的等权秩平均）@ k=37")
p(f"    命中 = **{caught} / 37**   冻结值 = **7**   ⇒ {'✔ 逐位复现' if caught == 7 else '✘ 复现失败'}")

# ---------------------------------------------------------------- (b) S2 Panel A 五行
S2 = {r["method"]: r for r in csv.DictReader(
    io.open(os.path.join(FIG, "tableS2_metrics.csv"), encoding="utf-8"))}
p("")
p("(b) S2 Panel A 那五行（此前缺列的五个信号）的 harmful-detection AUROC")
bad = []
for m in ["psi_max", "js_max", "pred_mean_shift", "entropy_mean_shift", "pred_pos_prev_shift"]:
    a = wmw(y, risk(m)); fz = float(S2[m]["harmful_auc"])
    ok = abs(a - fz) < 1e-9
    if not ok:
        bad.append(m)
    p(f"    {m:<22} 从 S1 复算 {a:.6f}   冻结 {fz:.6f}   {'✔' if ok else '✘'}")

p("")
ok_all = (caught == 7) and not bad
p("=" * 104)
p("*** p49 " + ("通过 —— S1 的首句在打印列上成立 ***" if ok_all else f"**未通过：{bad}** ***"))
REP.close()
sys.exit(0 if ok_all else 1)
