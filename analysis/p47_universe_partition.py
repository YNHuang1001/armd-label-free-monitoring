# -*- coding: utf-8 -*-
r"""p47 — 候选/合格宇宙的**分区审计**（只读，出耐久工件）

================================ 为什么单独一层 ================================
2026-10-07：清单里一句「其余 36（15 个组）」被外部核验指出是错的 ——
36 条非 KP_NIT 排除住在 **14** 个组，而 15 = 38−23 是**含 KP_NIT** 的未全合格组数。
该错误经由消息转述扩散（外部核验据它推出了 22/15/1 的分区，被 groupby 证伪）。
⇒ 分区断言必须有**耐久工件**可引，不能住在消息里（WORK_CHECKLIST 规则 15）。

================================ 本脚本产出 ================================
  · `figures/data/universe_partition_audit.csv` —— (group, site) 级计数
  · 下行四个分区的**逐组**明细
  · 报告 `report_p47_universe_partition.txt`

================================ 两个验证模式（外部核验提出，采纳）=============
  · **千分位数**（1,372,185 / 2,820,657 之流）→ 走**指纹扫描**（文本中可区分）
  · **小整数**（23 / 38 / 42 / 76 之流）→ 只能走**重算**（扫了等于没扫）
  本脚本就是小整数的重算器。
"""
import csv
import io
import os
import sys
from collections import Counter, defaultdict

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "分析结果")
FIG = os.path.join(ROOT, "figures", "data")
UNI = os.path.join(OUT, "p6b_corrected_universe.csv")
REP = open(os.path.join(OUT, "report_p47_universe_partition.txt"), "w", encoding="utf-8")


def p(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True); REP.write(s + "\n"); REP.flush()


U = list(csv.DictReader(io.open(UNI, encoding="utf-8")))
p("=" * 100)
p("p47 — 候选 / 合格宇宙分区审计（只读，出耐久工件）")
p("=" * 100)
p(f"上游：{UNI}")
p(f"总行数 {len(U)}")

GROUPS = sorted(set((r["organism"], r["drug"]) for r in U))
SITES = ["MGB", "SU", "UT"]
p(f"候选组 {len(GROUPS)} ｜ 站点 {len(SITES)} ｜ 方向/组 {len(U)//len(GROUPS)}")

# ---------------------------------------------------------------- (group, site) 级
cell = Counter((r["organism"] + "|" + r["drug"], r["source_site"]) for r in U)
rows = []
for g in GROUPS:
    for s in SITES:
        n = cell[(g[0] + "|" + g[1], s)]
        rows.append(dict(group=g[0] + "|" + g[1], site=s, n_candidate=n))
with io.open(os.path.join(FIG, "universe_partition_audit.csv"), "w", encoding="utf-8", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["group", "site", "n_candidate"])
    w.writeheader(); w.writerows(rows)
p(f"写出 {os.path.join(FIG, 'universe_partition_audit.csv')}（{len(rows)} 行 = {len(GROUPS)} 组 × {len(SITES)} 站）")

# ---------------------------------------------------------------- ① 均匀性
dist = Counter(cell.values())
p("")
p("① 组×站 均匀性")
p(f"   (group, site) 计数分布 = {dict(dist)}")
_uniform = set(cell.values()) == {2}
p(f"   ⇒ {'每组每站恰 2 条 ⇒ 228 = 38×2×3 是**设计不变量**' if _uniform else '✘ 不均匀'}")

# ---------------------------------------------------------------- ② 按组合格类
per = defaultdict(int)
for r in U:
    if str(r["eligible"]) == "True":
        per[(r["organism"], r["drug"])] += 1
for g in GROUPS:
    per[g] += 0                      # ★ 显式补 0：0 合格组否则永远不进字典（本轮实测过的 bug）
cls = Counter(per[g] for g in GROUPS)
p("")
p("② 按组「合格转移数」分布（★ 已显式补 0，防 0 合格组漏计）")
for k in sorted(cls, reverse=True):
    p(f"   {k}/6 合格: {cls[k]} 组")
_full, _zero = cls.get(6, 0), cls.get(0, 0)
_part = sum(v for k, v in cls.items() if 0 < k < 6)
p(f"   ⇒ 全合格 {_full} ｜ 部分 {_part} ｜ 零合格 {_zero} ｜ 合计 {_full + _part + _zero} 组")
p(f"   ⇒ 转移 {sum(k * v for k, v in cls.items())}")

# ---------------------------------------------------------------- ③ 排除 42 的两个正交分区
exc = [(r["organism"], r["drug"], r["source_site"]) for r in U if str(r["eligible"]) != "True"]
p("")
p(f"③ 排除 {len(exc)} 条的两个正交分区")
by_site = Counter(s for _, _, s in exc)
by_grp = Counter((o, d) for o, d, _ in exc)
p(f"   按站: {dict(by_site)}  合计 {sum(by_site.values())}")
p(f"   按组: {len(by_grp)} 个组贡献排除")
_kn = by_grp[("KP", "NIT")]
p(f"        KP_NIT {_kn} 条（该组 6 条候选全灭）")
p(f"        其余 {sum(by_grp.values()) - _kn} 条分布于 **{len(by_grp) - 1} 个组**")
net = Counter(s for o, d, s in exc if (o, d) != ("KP", "NIT"))
p(f"   净去 KP_NIT 的站点分解: {dict(net)}")
p("")
p(f"   ★ 措辞纪律：正确说法是「42 条分布于 {len(by_grp)} 个组」或"
  f"「其余 {sum(by_grp.values()) - _kn} 条分布于 {len(by_grp) - 1} 个组」")
p(f"     来自**含 KP_NIT** 的未全合格组数 = {len(GROUPS)} − {_full} = {len(GROUPS) - _full}")
p(f"     它不是「其余 36 条」的组数 —— 那是 {len(by_grp) - 1}")

# ---------------------------------------------------------------- ④ 正文组级句子的锚
p("")
p("④ 与正文组级句子的对账（这是小整数的**文本锚**，与重算互为双重出处）")
p(f"   正文：23 组贡献全部六条、10 组四条、4 组两条")
p(f"   实测：{_full} 组 6/6 ｜ {cls.get(4, 0)} 组 4/6 ｜ {cls.get(2, 0)} 组 2/6 ｜ {_zero} 组 0/6")
_ok = (_full == 23 and cls.get(4, 0) == 10 and cls.get(2, 0) == 4)
p(f"   ⇒ {'✔ 逐位吻合（正文组级句子通过**当前代**验证）' if _ok else '✘ 不符'}")

# ---------------------------------------------------------------- 断言
p("")
p("=" * 100)
_bad = []
if not _uniform: _bad.append("均匀性")
if _full + _part + _zero != len(GROUPS): _bad.append("组数不闭合")
# ★ 初版拿**合格**转移数（186）去对**候选**总数（len(U)=228）—— 断言写错了，被自己判红。
#   正确对照物是 eligible==True 的行数。
_N_ELIG = sum(1 for r in U if str(r["eligible"]) == "True")
if sum(k * v for k, v in cls.items()) != _N_ELIG:
    _bad.append(f"合格转移数不闭合（{sum(k * v for k, v in cls.items())} ≠ {_N_ELIG}）")
if not _ok: _bad.append("正文组级句子")
p("*** p47 " + ("通过" if not _bad else f"**未通过：{_bad}**") + " ***")
REP.close()
sys.exit(0 if not _bad else 1)
