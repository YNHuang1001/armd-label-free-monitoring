# -*- coding: utf-8 -*-
r"""_apply_batch_to_production.py —— 把 sim_fix 的 PATCHES + DELETIONS **批量应用到生产树**

★ 这是本项目至今不可逆性最高的一步。设计原则：

  1. **先备份**：9 个受影响文件先整份复制到 FREEZE_.../PRODUCTION_backup/CS_BATCH_20261006/，
     并落 sha256（备份后立即逐字节复核）。
  2. **逐行字节级断言**：每一处补丁在写入**之前**断言 `old` 精确命中（含 CRLF），
     不命中立即中止，不写任何文件 —— fail-closed。
  3. **与 sim_fix.apply_patch 逐字同规则**（含 CRLF 判定），不另立一套。
  4. **顺序**：PATCHES 按表内顺序（有 D→S 顺序替换链，顺序不可换）→ 再 DELETIONS。
  5. **收尾对拍**：应用后用 sha256 与 `evs/_patched/` **逐文件比对**（17 个）。
     两条独立路径（本脚本逐行替换 vs sim_fix 的补丁机制）必须给出同一结果。
     不一致 = 有东西没按预期落地，立即报告，不掩盖。
  6. **行数断言**：PATCHES 后每个文件行数必须不变；DELETIONS 后 p6f 必须恰好少 25 行。

用法：
    python _apply_batch_to_production.py            # 真正应用
    python _apply_batch_to_production.py --dry-run  # 只核验命中，不写盘
"""
import io
import os
import re
import sys
import hashlib
import shutil

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

HERE = os.path.dirname(os.path.abspath(__file__))
PROJ = os.path.dirname(HERE)
OUT = os.path.join(PROJ, "分析结果")
PATCHED = os.path.join(HERE, "_patched")
BACKUP = os.path.join(PROJ, "FREEZE_pre_correction_20261005",
                      "PRODUCTION_backup", "CS_BATCH_20261006")
DRY = "--dry-run" in sys.argv

# ---------------------------------------------------------------- 读表（不 import：sim_fix 是脚本）
_t = io.open(os.path.join(HERE, "sim_fix.py"), encoding="utf-8").read()
_ns = {}
for _name in ("PATCHES", "DELETIONS", "LOADER_FILES"):
    _m = re.search(r"^%s = \[(.*?)^\]" % _name, _t, re.S | re.M)
    exec("%s = [" % _name + _m.group(1) + "]", _ns)
PATCHES, DELETIONS, LOADER_FILES = _ns["PATCHES"], _ns["DELETIONS"], _ns["LOADER_FILES"]

print("=" * 100)
print("批量进生产  %s" % ("【DRY-RUN：只核验，不写盘】" if DRY else "【真实写入】"))
print("=" * 100)
print("  生产树 : %s" % OUT)
print("  备份到 : %s" % BACKUP)
print("  补丁 %d 处 / 删除 %d 段 / 生产文件 %d 个" % (len(PATCHES), len(DELETIONS), len(LOADER_FILES)))
print()


def sha(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()


def read_bytes(p):
    with open(p, "rb") as f:
        return f.read()


# ---------------------------------------------------------------- [1] 备份
touched = sorted({p[1] for p in PATCHES} | {d[1] for d in DELETIONS})
before = {fn: sha(os.path.join(OUT, fn)) for fn in touched}
before_lines = {fn: len(read_bytes(os.path.join(OUT, fn)).split(b"\n")) for fn in touched}

print("[1] 备份 %d 个受影响文件" % len(touched))
if not DRY:
    if os.path.isdir(BACKUP):
        shutil.rmtree(BACKUP)
    os.makedirs(BACKUP)
    for fn in touched:
        shutil.copy2(os.path.join(OUT, fn), os.path.join(BACKUP, fn))
        assert sha(os.path.join(BACKUP, fn)) == before[fn], "备份后 sha256 不一致: " + fn
for fn in touched:
    print("    %-42s %s  %d 行" % (fn, before[fn][:16], before_lines[fn]))
print("    ✓ 全部备份并逐字节复核" if not DRY else "    （dry-run：跳过）")
print()

# ---------------------------------------------------------------- [2] 应用补丁
print("[2] 应用 PATCHES（表内顺序；D→S 是顺序替换链，顺序不可换）")
cur, napplied, miss = {}, 0, []
for tag, fn, ln, old, new, why in PATCHES:
    p = os.path.join(OUT, fn)
    if fn not in cur:
        raw = read_bytes(p)
        cur[fn] = {"crlf": b"\r\n" in raw, "lines": raw.split(b"\n")}
    c = cur[fn]
    nl = b"\r" if c["crlf"] else b""
    A = c["lines"]
    if not (1 <= ln <= len(A)):
        miss.append((tag, fn, ln, "行号越界"))
        continue
    old_b, new_b = old.encode("utf-8") + nl, new.encode("utf-8") + nl
    if A[ln - 1] != old_b:
        miss.append((tag, fn, ln, "old 未命中\n      expect %r\n      actual %r"
                     % (old_b[:110], A[ln - 1][:110])))
        continue
    A[ln - 1] = new_b
    napplied += 1
    print("    [%s] %-38s :%-4d  OK   %s" % (tag, fn, ln, why[:58]))

print("    已应用 %d / %d" % (napplied, len(PATCHES)))
if miss:
    print()
    print("  !! 有 %d 处未命中 —— **中止，不写任何文件**：" % len(miss))
    for m in miss:
        print("     [%s] %s:%s  %s" % m)
    sys.exit(2)
print()

# ---------------------------------------------------------------- [3] 应用删除
print("[3] 应用 DELETIONS（在 PATCHES 之后；行数会变）")
for tag, fn, a, b, why in DELETIONS:
    p = os.path.join(OUT, fn)
    if fn not in cur:
        raw = read_bytes(p)
        cur[fn] = {"crlf": b"\r\n" in raw, "lines": raw.split(b"\n")}
    A = cur[fn]["lines"]
    seg = A[a - 1:b]
    head = seg[0].decode("utf-8", "replace")[:60]
    tail = seg[-1].decode("utf-8", "replace")[:60]
    print("    [%s] %s 删 %d-%d（%d 行）" % (tag, fn, a, b, len(seg)))
    print("         首行: %s" % head)
    print("         末行: %s" % tail)
    del A[a - 1:b]
    print("         -> 剩余 %d 行（原 %d，应少 %d）"
          % (len(A), len(A) + len(seg), len(seg)))
print()

# ---------------------------------------------------------------- [4] 写盘 + 行数断言
print("[4] 写盘并断言行数")
if not DRY:
    for fn, c in cur.items():
        p = os.path.join(OUT, fn)
        with open(p, "wb") as f:
            f.write(b"\n".join(c["lines"]))
    for fn, c in cur.items():
        now = len(read_bytes(os.path.join(OUT, fn)).split(b"\n"))
        exp = before_lines[fn] + sum(-(d[3] - d[2] + 1) for d in DELETIONS if d[1] == fn)
        status = "OK" if now == exp else "行数不符!"
        print("    %-42s %d -> %d  (期望 %d)  %s" % (fn, before_lines[fn], now, exp, status))
        assert now == exp, "行数不符: " + fn
else:
    print("    （dry-run：跳过写盘）")
print()

# ---------------------------------------------------------------- [5] 与 _patched 对拍
print("[5] 与 evs/_patched/ 逐文件对拍（两条独立路径必须一致）")
if DRY:
    print("    （dry-run：生产未写盘，与 _patched 必然不同 —— 此步跳过。"
          "dry-run 的用途只是核验 37 处 old 是否全部命中。）")
    print()
    print("=" * 100)
    print("dry-run 结束：%d/%d 处补丁命中；未写任何文件。" % (napplied, len(PATCHES)))
    print("=" * 100)
    sys.exit(0)
bad = []
for fn in LOADER_FILES:
    a = os.path.join(OUT, fn)
    b = os.path.join(PATCHED, fn)
    if not os.path.exists(b):
        print("    %-42s _patched 无此文件" % fn)
        continue
    sa, sb = sha(a), sha(b)
    ok = sa == sb
    if not ok:
        bad.append(fn)
    print("    %-42s %s  %s" % (fn, sa[:16], "OK" if ok else
          "DIFF (patched %s)" % sb[:16]))
print()
if bad:
    print("  !! %d 个文件与 _patched 不一致：%s" % (len(bad), bad))
    print("     这不代表应用失败 —— 也可能 _patched 本身不是最新。**必须人工查清后再决定**。")
    sys.exit(3)

print("=" * 100)
print("完成：%d 处补丁 + %d 段删除已%s；17 个生产文件与 _patched **逐字节一致**。"
      % (len(PATCHES), len(DELETIONS), "（dry-run 未写盘）" if DRY else "写入生产"))
print("备份：%s" % BACKUP)
print("=" * 100)
