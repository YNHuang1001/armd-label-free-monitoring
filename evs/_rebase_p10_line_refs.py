# -*- coding: utf-8 -*-
r"""_rebase_p10_line_refs.py —— 把 p10a / p10b 的**行号引用**平移到插入代码之后的位置

★ 事故（2026-10-06，本人造成）：为加"缓存产出代码签名 + resume 依赖检查"，给
  `p10a_prepare_cache.py` 插入了 32 行、给 `p10b_subsample_curves.py` 插入了 49 行。
  这两个文件里**有大量硬编码行号引用**（evs_load 的 BLOCK_*、t01/t07 的 exec_block、
  manifest 的 anchors），插入后**全部错位** —— 若不改，FORK 系列会去测**错误的代码块**。

★ 位移图（用 difflib 逐段算出，并与 7 个**按内容**定位到的块交叉验证，见下方 verify）：
    p10a : 旧 >=35 -> +18 ；>=333 -> +19 ；>=344 -> +32
    p10b : 旧 31..61 -> +13 ；62..201 -> +21 ；202..203 -> +30 ；204..206 -> +46 ；207..211 -> +47 ；>=212 -> +49

★ 脚本只改**行号数字**，不改任何逻辑；改完逐条**用内容断言**验证。
"""
import io
import os
import re
import sys

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

HERE = os.path.dirname(os.path.abspath(__file__))
PROJ = os.path.dirname(HERE)
OUT = os.path.join(PROJ, "分析结果")
PAT = os.path.join(HERE, "_patched")          # 合并后、我动手前的状态 = 行号的"旧"基准
DRY = "--dry-run" in sys.argv


def shift_p10a(n):
    if n >= 344:
        return n + 32
    if n >= 333:
        return n + 19
    if n >= 35:
        return n + 18
    return n


def shift_p10b(n):
    if n >= 212:
        return n + 49
    if n >= 207:
        return n + 47
    if n >= 204:
        return n + 46
    if n >= 202:
        return n + 30
    if n >= 62:
        return n + 21
    if n >= 31:
        return n + 13
    return n


# 需要改的 (文件, 正则, 用哪张位移表, 组号是行号的组)
TARGETS = [
    ("evs_load.py",        r"BLOCK_P10B_SUDO = \(P10B, (\d+), (\d+)\)", "b"),
    ("evs_load.py",        r"BLOCK_P10B_L3 = \(P10B, (\d+), (\d+)\)",   "b"),
    ("evs_load.py",        r"BLOCK_P10B_L2 = \(P10B, (\d+), (\d+)\)",   "b"),
    ("evs_load.py",        r"BLOCK_P10B_SENS = \(P10B, (\d+), (\d+)\)", "b"),
    ("evs_load.py",        r"BLOCK_P10A_MMD = \(P10A, (\d+), (\d+)\)",  "a"),
    ("evs_load.py",        r"BLOCK_BW_P10A = \(P10A, (\d+), (\d+)\)",   "a"),
    ("t01_mmd.py",         r"\(L\.P10A, (\d+), (\d+)\)",                 "a"),
    ("t01_mmd.py",         r"\(L\.P10B, (\d+), (\d+)\)",                 "b"),
    ("t07_fork_parity.py", r"\(L\.P10A, (\d+), (\d+)\)",                 "a"),
    ("t07_fork_parity.py", r"\(L\.P10B, (\d+), (\d+)\)",                 "b"),
]
SH = {"a": shift_p10a, "b": shift_p10b}

print("=" * 100)
print("平移 p10a / p10b 的行号引用  %s" % ("【DRY-RUN】" if DRY else "【真实写入】"))
print("=" * 100)

total = 0
for fn, pat, which in TARGETS:
    p = os.path.join(HERE, fn)
    src = io.open(p, encoding="utf-8").read()
    f = SH[which]
    hits = []

    def sub(m):
        a, b = int(m.group(1)), int(m.group(2))
        na, nb = f(a), f(b)
        if (na, nb) != (a, b):
            hits.append((a, b, na, nb))
        return m.group(0).replace("%d, %d" % (a, b), "%d, %d" % (na, nb))

    new = re.sub(pat, sub, src)
    if hits:
        total += len(hits)
        print("  %-22s %-46s %d 处" % (fn, pat[:44], len(hits)))
        for a, b, na, nb in hits:
            print("        %d-%d  ->  %d-%d" % (a, b, na, nb))
        if not DRY:
            io.open(p, "w", encoding="utf-8", newline="\n").write(new)
print("\n  合计 %d 处" % total)
print()

# ---------------------------------------------------------------- 验证
print("[验证] 按新行号从生产抽块，与 _patched 按旧行号抽块**逐字节比对**")
import importlib
sys.path.insert(0, HERE)
import evs_load as L
importlib.reload(L)
CHECKS = [
    ("BLOCK_BW_P10A",  (L.P10A, 224, 228), (L.P10A, 242, 246)),
    ("BLOCK_P10A_MMD", (L.P10A, 273, 280), (L.P10A, 291, 298)),
    ("p10b L1b",       (L.P10B, 85, 93),   (L.P10B, 106, 114)),
    ("BLOCK_P10B_L2",  (L.P10B, 107, 114), (L.P10B, 128, 135)),
    ("BLOCK_P10B_SUDO", (L.P10B, 116, 148), (L.P10B, 137, 169)),
    ("BLOCK_P10B_L3",  (L.P10B, 150, 174), (L.P10B, 171, 195)),
    ("BLOCK_P10B_SENS", (L.P10B, 256, 263), (L.P10B, 305, 312)),
]
allok = True
for name, old, new in CHECKS:
    fo = os.path.join(PAT, os.path.basename(old[0]))
    fn_ = os.path.join(OUT, os.path.basename(new[0]))
    o = io.open(fo, encoding="utf-8").read().split("\n")[old[1] - 1:old[2]]
    n = io.open(fn_, encoding="utf-8").read().split("\n")[new[1] - 1:new[2]]
    ok = (o == n)
    allok &= ok
    print("    %-18s 旧 %-11s -> 新 %-11s  相同=%s" % (name, "%d-%d" % old[1:], "%d-%d" % new[1:], ok))
print()
print("== %s ==" % ("全部通过" if allok else "**有不同**"))
sys.exit(0 if allok else 3)
