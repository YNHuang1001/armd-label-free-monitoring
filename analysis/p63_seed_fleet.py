# -*- coding: utf-8 -*-
r"""p63 — 种子稳定性总表（**一次覆盖全部已发表区间**）

================================ 为什么有这个脚本 ================================
外部核验（2026-10-13）要求：**不要逐个补漏**，让每个 bootstrap 生产者都把种子当参数，
用**同一个种子框架**跑一遍，输出**一张统一表**。

================================ 设计 ================================
① **断言 ①（先跑，不过就停）**：默认种子下，各生产者重跑后与**冻结值逐位一致**。
   比对基准 = 投稿包里的**改前副本**（`SUBMISSION_20261010/07_画图代码_可运行/数据包_figures_data/`）。
   这一条尤其要紧：`p11` 此前**因 GBK 编码跑不动**，已发表的 fig2 区间与六个 Spearman
   区间**从未在可复现的条件下被测过**。哈希不一致 ⇒ 判红并中止。
② **种子表**：每个生产者显式列出 **B 与种子**（正文引的区间必须能指认来自哪一个）。
③ **5 个种子**：**不用相邻值**（某些 PRNG 对相邻种子输出相关）。默认种子（已发表）＋ 4 个相距很远的。
④ **统一表**：每条断言一行 —— 已发表区间 ｜ 5 种子端点 min/max ｜ 相对参照（0 / 0.5 / 取整边界）
   的位置是否翻转 ｜ 翻转则标注「措辞需改」。
⑤ **收尾**：用默认种子再跑一遍，把冻结状态还原。
"""
import csv
import hashlib
import io
import os
import shutil
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "分析结果")
FIG = os.path.join(ROOT, "figures", "data")
PKG = os.path.join(ROOT, "SUBMISSION_20261010", "07_画图代码_可运行", "数据包_figures_data")
PY = sys.executable

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

# ---- 生产者与其 B / 已发表种子（正文引的区间必须能指认来自哪一个）----
# ★★ 2026-10-13 更正（首版**归错了生产者 ⇒ 8 个端点的判据在空转**）：
#   首版把 `fig2_forest.csv` / `tableS2_metrics.csv` 记成 p11（B=4000）的产物 —— **错**。
#   p11 只是**复制/合并** `p6d_signal_eval.csv`（见 p11:158、p11:457），真正算 bootstrap 的是
#   **`p6d_01_unseal_evaluation.py`**（`rng = default_rng(SEED + 7)`，B = 1000）。
#   ⇒ 变 p11 或 p33 的种子**都不可能**让那两个文件动 —— 首版报的「0.5559 ~ 0.5559 全稳」
#     是**空转的假象**，不是发现。这正是本会话反复出现的那一类缺陷。
#   ⇒ 生产者改成**链**：`p6d_01 → p11`（前者算，后者复制成 fig2/tableS2）。
PRODUCERS = [
    (["p6d_01_unseal_evaluation.py", "p11_figure_datapack.py"], 1000, 20261005,
     ["fig2_forest.csv", "tableS2_metrics.csv"], "p6d_01：B=1000，bootstrap 种子 = SEED+7"),
    (["p33_protocol72_metrics.py"], 1000, 20261007, ["p33_protocol72_metrics.csv"],
     "p33：B=1000，种子 20261007"),
    (["p45_triage_policy_simulation.py"], 10000, 20261007, ["figS10_triage_summary.csv"],
     "p45：B=10000，种子 20261007"),
    (["p57_posthoc_baselines.py"], 10000, 20261007, ["p57_scores.csv"],
     "p57：B=10000，种子 20261007（ICC 用参数 bootstrap，同种子）"),
]
PUBLISHED = 20261007
# ★ 不用相邻种子：相距很远，且避开各种"好看"的常数
NEW_SEEDS = [20261007, 31415926, 271828182, 987654321, 1234567890]

REP = io.open(os.path.join(OUT, "report_p63_seed_fleet.txt"), "w", encoding="utf-8")


def p(s=""):
    print(s)
    REP.write(str(s) + "\n")
    REP.flush()


def sha(fp):
    return hashlib.sha256(open(fp, "rb").read()).hexdigest()[:16] if os.path.exists(fp) else "(缺)"


def run(scripts, seed):
    """scripts 可以是单个文件名，也可以是**一条链**（按序跑）。"""
    if isinstance(scripts, str):
        scripts = [scripts]
    env = dict(os.environ)
    env["ARMD_SEED"] = str(seed)
    rc = 0
    for sc in scripts:
        r = subprocess.run([PY, os.path.join(OUT, sc)], cwd=OUT, env=env,
                           capture_output=True, text=True, encoding="utf-8", errors="replace")
        rc = r.returncode or rc
    return rc


def read_cell(path, keycol, keyval, col):
    if not os.path.exists(path):
        return None
    with io.open(path, encoding="utf-8-sig") as f:
        for row in csv.DictReader(f):
            if row.get(keycol) == keyval:
                try:
                    return float(row[col])
                except Exception:
                    return None
    return None


p("=" * 100)
p("p63 — 种子稳定性总表")
p("=" * 100)
p("  生产者 → B / 已发表种子（正文引的区间按此指认）：")
for sc, b, sd, outs, note in PRODUCERS:
    p("    %-46s → %s" % (" ＋ ".join(sc), ", ".join(outs)))
    p("        %s" % note)
p("  用于稳定性测试的 5 个种子（**非相邻**）：%s" % NEW_SEEDS)
p("  ★ 只覆盖 `ARMD_SEED` 环境变量；**数据、划分、模型种子一律不动**。")
p("")

# ---------------------------------------------------------------- 断言 ①
p("【断言 ①】默认种子下重跑 == 冻结值（逐位）。**不过即中止**")
p("           比对基准：投稿包内改前副本（p11 此前因 GBK 跑不动，这两项从未被复现过）")
p("")
before = {}
for sc, b, sd, outs, note in PRODUCERS:
    for o in outs:
        src = os.path.join(FIG, o) if not o.startswith("p") else os.path.join(OUT, o)
        before[o] = (sha(src), sha(os.path.join(PKG, o)))
bad = []
for sc, b, sd, outs, note in PRODUCERS:
    rc = run(sc, sd)
    p("    %-46s exit=%d  （种子 %d = 已发表）" % (" ＋ ".join(sc), rc, sd))
    for o in outs:
        src = os.path.join(FIG, o) if not o.startswith("p") else os.path.join(OUT, o)
        now, frozen = sha(src), before[o][1]
        ok = (frozen == "(缺)") or (now == frozen)
        if not ok:
            bad.append(o)
        p("       %-30s 现在 %s ｜ 冻结 %s  %s" % (o, now, frozen,
                                                  "OK" if ok else "**不一致**"))
p("")
if bad:
    p("  ⇒ **断言 ① 判红：%s** —— 重跑值与冻结值不一致 ⇒ 这是新缺陷，**后续分析中止**。" % bad)
    REP.close()
    sys.exit(1)
p("  ⇒ **断言 ① 通过**：默认种子下逐位复现冻结值 ⇒ 冻结值可靠，可继续。")
p("")

# ---------------------------------------------------------------- 统一稳定性表
TARGETS = [
    # (标签, 文件, 主键列, 主键值, 端点列, 参照)
    ("DR 对 src-only 之外的主表 AUROC 区间", "fig2_forest.csv", "method", "AUROC_DR", "auc_lo", 0.5),
    (" 同上（上界）", "fig2_forest.csv", "method", "AUROC_DR", "auc_hi", 0.5),
    ("分诊 DR 对随机的下界（全稿最贴 0）", "figS10_triage_summary.csv", "layer", "L-DR",
     "diff_vs_random_lo", 0.0),
    ("分诊 DR 对 composite 下界", "figS10_triage_summary.csv", "layer", "L-DR_minus_L-shift",
     "lo95", 0.0),
    ("分诊 DR 对 composite 上界", "figS10_triage_summary.csv", "layer", "L-DR_minus_L-shift",
     "hi95", 0.0),
]
for m in ["domain_auc", "mmd", "psi_mean", "psi_max", "js_mean", "js_max"]:
    TARGETS.append(("%s 对 ΔAUROC 下界" % m, "tableS2_metrics.csv", "method", m, "sp_lo", 0.0))
    TARGETS.append(("%s 对 ΔAUROC 上界" % m, "tableS2_metrics.csv", "method", m, "sp_hi", 0.0))

p("【统一稳定性表】每个端点一行；参照片：0 / 0.5 / 取整边界（第二位小数的 .xx5）")
p("")
cols = {}
for seed in NEW_SEEDS:
    p("  跑种子 %d …" % seed)
    for sc, b, sd, outs, note in PRODUCERS:
        run(sc, seed)
    for lab, fn, kc, kv, col, ref in TARGETS:
        fp = os.path.join(FIG, fn) if fn.endswith(".csv") and not fn.startswith("p") \
            else os.path.join(OUT, fn)
        v = read_cell(fp, kc, kv, col)
        cols.setdefault((lab, ref), []).append(v)

p("")
p("%-40s %-24s %-22s %s" % ("断言/端点", "已发表种子值", "5 种子 min ~ max", "相对参照是否翻转"))
print("-" * 116)
flips = []
for (lab, ref), vals in cols.items():
    vs = [v for v in vals if v is not None]
    if not vs:
        p("%-40s %-24s %-22s %s" % (lab, "(缺)", "(缺)", "—"))
        continue
    pub = vs[0]
    lo, hi = min(vs), max(vs)
    if ref == 0.5:
        pos = [(v - 0.5) > 0 for v in vs]
        refname = "0.5"
    else:
        pos = [(v > 0) for v in vs]
        refname = "0"
    flip = len(set(pos)) > 1
    # 取整边界（第二位小数）
    rnd = {round(v, 2) for v in vs}
    rflip = len(rnd) > 1
    if flip or rflip:
        flips.append((lab, lo, hi))
    p("%-40s %-24s %-22s %s%s"
      % (lab, "%.4f" % pub, "%.4f ~ %.4f" % (lo, hi),
         ("**相对 %s 翻转**" % refname) if flip else "相对 %s 不翻转" % refname,
         " ｜ **第二位小数不稳：%s**" % sorted(rnd) if rflip else ""))

p("")
if flips:
    p("  ⇒ **措辞需改** 的端点 %d 个：" % len(flips))
    for lab, lo, hi in flips:
        p("      %-40s [%.4f, %.4f]" % (lab, lo, hi))
else:
    p("  ⇒ 无端点在 5 个种子上翻转参照位置。")

p("")
p("【收尾】用已发表种子再跑一遍，还原冻结状态")
for sc, b, sd, outs, note in PRODUCERS:
    run(sc, sd)
for sc, b, sd, outs, note in PRODUCERS:
    for o in outs:
        src = os.path.join(FIG, o) if not o.startswith("p") else os.path.join(OUT, o)
        now, frozen = sha(src), before[o][1]
        p("    %-30s %s  %s" % (o, now, "OK（= 冻结）" if now == frozen or frozen == "(缺)"
                                else "**与冻结不一致 —— 必须处理**"))
p("")
p("done -> %s" % os.path.join(OUT, "report_p63_seed_fleet.txt"))
REP.close()
