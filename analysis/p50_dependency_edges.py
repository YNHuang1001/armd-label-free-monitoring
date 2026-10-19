# -*- coding: utf-8 -*-
r"""p50 — **依赖边登记**（A3.4）：把「谁住在谁下游」写成可机读的边，而非默契

================================ 为什么需要（规则 9 的机械化）================================
A3.3 的代次头回答了「**这个产物是不是当前宇宙的**」——
但那只在被问到的时候才回答。**缺的是「谁该被问」**：
`task_profile.csv` 陈旧三轮而无人发现，不是因为查不出来，是因为**没有任何东西声明
「p21 链在 p6b 下游」**。宇宙再冻结时，没有任何机制把下游标脏。

⇒ 本模块把边写成**声明**（`EDGES`，手写、单一出处），并据此逐产物判定：
  · `current`   —— 代次头与当前宇宙一致
  · `STALE`     —— 代次头存在但与当前宇宙**不一致** ⇒ 必须重跑
  · `UNSTAMPED` —— **无代次头** ⇒ 无法判断，按未声明处理（规则 15 精神）

**本清单直接喂 A5.1（重建投稿包）的检查表**：任何 STALE / UNSTAMPED 都不得进包。

用法：python p50_dependency_edges.py   （exit 0 = 全部 current）
"""
import io
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import genstamp as gs

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
AN = os.path.join(ROOT, "分析结果")
# ★ 2026-10-08：**root 覆盖** —— 供**负向测试在副本上跑**（外部核验的三合一修法：
#   副本上跑 ＋ tee 落盘 ＋ 引用时自带围栏）。默认仍是真 root；测试传 P50_DATA_ROOT。
_DR = os.environ.get("P50_DATA_ROOT")
RP = os.path.join(_DR, "reporting_profiles") if _DR else os.path.join(ROOT, "figures", "data", "reporting_profiles")
FG = _DR if _DR else os.path.join(ROOT, "figures", "data")

# ================================================================ 依赖边（手写，单一出处）
# 每边：上游 → 下游产物（＋产出它的脚本）。**加新产物必须在此登记**，否则等于没有边。
# ★ 2026-10-08（Q3 裁决）：三元组由 (path, 生产者) 升级为 **(path, 生产者, expected_final)——
#   「本产物的合法终态应由谁收尾」**。
#   病因：p21a 与 p21c **都写** `site_profile.csv`（p21a 写原始口径版；p21c 读它、合并 te/rs 后
#   **覆盖**回真名）。修复「p21a 也盖章」之后，**若 p21a 收尾**，sidecar 记 p21a ＋ 原版 sha，
#   **自洽**，p50 两道检查全绿 —— 而文件里躺的是**未合并、缺 5 列**的原版。
#   ⇒ sha 校验抓「未盖章的覆盖」，**看不见「序」**。这是「绿但错版」，
#     与让 202 活三轮的形状同一，只上移了一层。
#   边界（如实记账）：p21a 版**不是**合法终态（p21d 会因缺列而崩，是**响亮**失败）；
#   但 p50 会为它开绿灯，且**没有任何东西检查 schema**。
EDGES = [
    dict(upstream="p6b_corrected_universe.csv",
         artifacts=[(os.path.join(RP, "site_profile.pending.csv"), "p21b", "p21b"),
                    (os.path.join(RP, "site_profile.csv"), "p21c", "p21c"),
                    (os.path.join(RP, "Supplementary_Tables_R1_R3.md"), "p21d", "p21d"),
                    (os.path.join(RP, "supplementary_cell_provenance.csv"), "p21d", "p21d")]),
    dict(upstream="p6b_label_free_signals_blinded.csv + p6c_L3_blinded.csv",
         artifacts=[(os.path.join(FG, "tableS1_tasks.csv"), "p11", "p11"),
                    (os.path.join(FG, "figS10_triage_curves.csv"), "p45", "p45")]),
    dict(upstream="p6d_outcome.csv",
         artifacts=[(os.path.join(FG, "tableS3_triage.csv"), "p46", "p46")]),
]

rep = open(os.path.join(AN, "report_p50_dependency_edges.txt"), "w", encoding="utf-8")


def p(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True); rep.write(s + "\n"); rep.flush()


p("=" * 104)
p("p50 — 依赖边登记与代次判定（A3.4）")
p("=" * 104)
p(f"当前宇宙：{gs.universe_id()}")
p("")
stale, unstamped, current = [], [], []
for e in EDGES:
    p(f"【上游】{e['upstream']}")
    for path, prod, expected in e["artifacts"]:
        name = os.path.basename(path)
        if not os.path.exists(path):
            p(f"    {name:<40} ✘ **产物缺失**（生产者 {prod}）")
            stale.append(name); continue
        r = gs.is_current(path)
        g = gs.read(path) or {}
        # ★ 2026-10-08 补洞：**sidecar 的 artifact_sha 必须与实际文件相符** ——
        #   否则「盖章后被别的脚本覆盖」会静默通过（p21a 覆盖 p21c 写过的 site_profile.csv 即为该形态）。
        import hashlib as _hl
        _real = _hl.sha256(open(path, "rb").read()).hexdigest()[:16]
        if g.get("artifact_sha") and g["artifact_sha"] != _real:
            p(f"    {name:<40} ✘ **工件与代次头不符**（头 {g['artifact_sha']} ≠ 实 {_real}）"
              f"⇒ 盖章后被覆盖；生产者 {prod}")
            stale.append(name); continue
        # ★ Q3：**序**检查 —— sidecar 记的产出者必须是声明的收尾者
        if g.get("producer") and not str(g["producer"]).startswith(expected):
            stale.append(name)
            p(f"    {name:<40} ✘ **错序**（sidecar 产出者 {g['producer']} ≠ 声明的收尾者 {expected}）"
              f"⇒ **被中断的链**；文件可能是非终态版本")
            continue
        if r is True:
            current.append(name)
            p(f"    {name:<40} ✔ current     （{prod}｜{g.get('universe_sha')}｜收尾 {g.get('producer')}）")
        elif r is False:
            stale.append(name)
            p(f"    {name:<40} ✘ **STALE**     （{prod}｜头 {g.get('universe_sha')} "
              f"≠ 当前 {gs.universe_id()['universe_sha']}）⇒ **必须重跑**")
        else:
            unstamped.append(name)
            p(f"    {name:<40} ⚠ **UNSTAMPED** （无代次头；生产者 {prod}）⇒ 无法判断")
    p("")

p("=" * 104)
p(f"判定汇总：current {len(current)} ｜ STALE {len(stale)} ｜ UNSTAMPED {len(unstamped)}")
if stale:
    p(f"  ⚠ STALE：{stale} —— 这些产物与当前宇宙不同代，**不得进包**（A5.1 检查表）")
if unstamped:
    p(f"  ⚠ UNSTAMPED：{unstamped} —— 生产者未接 genstamp，代次不可判")
ok = not stale and not unstamped
p("")
p("*** p50 " + ("通过 —— 全部下游产物与当前宇宙同代 ***" if ok
                else "**未通过 —— 存在 STALE 或 UNSTAMPED** ***"))
rep.close()
sys.exit(0 if ok else 1)
