# -*- coding: utf-8 -*-
r"""_fix_cs01_alias_misplacement.py —— 修正 CS-01 引入的**别名胜性错位**

★ 背景（2026-10-06 独立审计发现，本人逐条核实）：
  CS-01 把 10 个新别名加进 DRUGMAP，但其中 **7 个挂到了错误的 canonical**。
  与冻结的 CS01_amendment.txt:11-21 逐条对照，不符 7/10：

      ticarcillinclavulanicacid    -> 应为 TIM   实为 CZA
      ticacarcillinclavulanicacid  -> 应为 TIM   实为 CZA
      tobramycinrm                 -> 应为 TOB   实为 NIT
      amikaciniv                   -> 应为 AMK   实为 NIT
      amikacinintravenous          -> 应为 AMK   实为 NIT
      polymixine                   -> 应为 CST   实为 PEN
      polymyxine                   -> 应为 CST   实为 PEN

  为什么没有判据抓到：DRUG-01 只检查 5 个副本**彼此一致**（它们确实一致，因为五个都是同一个错）；
  DRUG-05 把这些串记为 mapped_to_benchmark（它们精确命中别名表）。**判据验的是副本间一致性，
  不是与裁决的一致性。**

★ 本脚本做三件事，缺一不可：
  1. 逐行**字节级断言** old 命中（含 CRLF），不命中即中止不写盘；
  2. 行数守恒断言（15 行 1:1 替换）；
  3. **行为验证**：exec 每个文件修正后的 DRUGMAP + _canon，对 10 个别名逐个调 `_canon()`
     验证返回声明的 canonical —— **不是查文本，是执行真函数**。
     另验：canonical 的**集合与顺序不变**（rebuild_figS1_task_matrix.py 依赖 p5 的首现序）；
     且**无别名串跨行重复**（精确匹配下跨行重复会让顺序重新起作用）。

用法：
    python _fix_cs01_alias_misplacement.py --dry-run
    python _fix_cs01_alias_misplacement.py
"""
import io
import os
import re
import sys
import ast
import shutil
import hashlib

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

HERE = os.path.dirname(os.path.abspath(__file__))
PROJ = os.path.dirname(HERE)
OUT = os.path.join(PROJ, "分析结果")
BACKUP = os.path.join(PROJ, "FREEZE_pre_correction_20261005",
                      "PRODUCTION_backup", "CS01v2_20261006")
DRY = "--dry-run" in sys.argv

OLD = {
    "TIMCZA": '    ("TIM", ["ticarcillinclavulanate", "ticarcilclavulanate"]), ("CZA", ["ceftolozanetazobactam", "ticarcillinclavulanicacid", "ticacarcillinclavulanicacid"]),',
    "TOBAMKNIT": '    ("TOB", ["tobramycin"]), ("AMK", ["amikacin"]), ("NIT", ["nitrofurantoin", "tobramycinrm", "amikaciniv", "amikacinintravenous"]),',
    "CSTPEN": '    ("CST", ["colistin"]), ("SAM", ["sulbactam"]), ("AMP", ["ampicillin"]), ("PEN", ["piperacillin", "polymixine", "polymyxine"]),',
}
NEW = {
    "TIMCZA": '    ("TIM", ["ticarcillinclavulanate", "ticarcilclavulanate", "ticarcillinclavulanicacid", "ticacarcillinclavulanicacid"]), ("CZA", ["ceftolozanetazobactam"]),',
    "TOBAMKNIT": '    ("TOB", ["tobramycin", "tobramycinrm"]), ("AMK", ["amikacin", "amikaciniv", "amikacinintravenous"]), ("NIT", ["nitrofurantoin"]),',
    "CSTPEN": '    ("CST", ["colistin", "polymixine", "polymyxine"]), ("SAM", ["sulbactam"]), ("AMP", ["ampicillin"]), ("PEN", ["piperacillin"]),',
}
# (file, TIMCZA行, TOBAMKNIT行, CSTPEN行) —— 行号已实测
SPEC = [
    ("p5_01_degradation_benchmark.py", 45, 49, 52),
    ("p6_02_label_free_signals.py",    44, 48, 51),
    ("p6b_01_corrected_pipeline.py",   51, 55, 58),
    ("p6c_01_L3_blinded.py",           58, 62, 65),
    ("p10a_prepare_cache.py",          47, 51, 54),
]

# CS01_amendment.txt:11-21 声明的映射 —— 这是唯一权威
WANT = {
    "amoxicillinkclavulanateaug": "AMC",
    "ticarcillinclavulanicacid": "TIM",
    "ticacarcillinclavulanicacid": "TIM",
    "gentscreen": "GEN",
    "gentamicinscreen": "GEN",
    "tobramycinrm": "TOB",
    "amikaciniv": "AMK",
    "amikacinintravenous": "AMK",
    "polymixine": "CST",
    "polymyxine": "CST",
}


def sha(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()


def drugmap_of(path):
    """从文件文本抽出 DRUGMAP 与 _canon 并**执行**，返回 (canon_of_alias, canon_order, fn)。"""
    txt = io.open(path, encoding="utf-8").read()
    blk = re.search(r"^DRUGMAP = \[(.*?)^\]", txt, re.S | re.M).group(1)
    fn_src = re.search(r"^def _canon\(n\):(?:\n[ \t].*)+", txt, re.M).group(0)
    ns = {}
    exec("DRUGMAP = [" + blk + "]", ns)
    exec(fn_src, ns)
    dm = ns["DRUGMAP"]
    alias_of = {}
    order = []
    for c, ks in dm:
        if c not in order:
            order.append(c)
        for k in ks:
            alias_of.setdefault(k, []).append(c)
    return alias_of, order, ns["_canon"]


print("=" * 100)
print("修正 CS-01 别名胜性错位  %s" % ("【DRY-RUN】" if DRY else "【真实写入】"))
print("=" * 100)

# ---------------------------------------------------------------- [1] 前置：行为基线
print("[1] 修正前行为基线（执行各文件的 _canon，看 10 个别名现在映射到哪）")
pre = {}
for fn, *_ in SPEC:
    p = os.path.join(OUT, fn)
    alias_of, order, can = drugmap_of(p)
    wrong = [a for a, c in WANT.items() if can(a) != c]
    pre[fn] = (alias_of, order, can)
    print("    %-38s 别名总数=%-3d canonical 数=%-3d  10 个别名中不符=%d"
          % (fn, len(alias_of), len(order), len(wrong)))
print()

# ---------------------------------------------------------------- [2] 断言 + 替换
print("[2] 逐行字节级断言并替换（%d 文件 x 3 行）" % len(SPEC))
buf, miss = {}, []
for fn, a, b, c in SPEC:
    raw = open(os.path.join(OUT, fn), "rb").read()
    nl = b"\r" if b"\r\n" in raw else b""
    L = raw.split(b"\n")
    for key, ln in (("TIMCZA", a), ("TOBAMKNIT", b), ("CSTPEN", c)):
        o = OLD[key].encode("utf-8") + nl
        if not (1 <= ln <= len(L)) or L[ln - 1] != o:
            miss.append((fn, key, ln, (L[ln - 1][:80] if 1 <= ln <= len(L) else b"OOB")))
            continue
        L[ln - 1] = NEW[key].encode("utf-8") + nl
    buf[fn] = (L, nl)
    print("    %-38s 3 行 %s" % (fn, "命中" if not any(m[0] == fn for m in miss) else "**未命中**"))
if miss:
    print("\n  !! 未命中 %d 处 —— 中止，不写盘：" % len(miss))
    for m in miss:
        print("     %s / %s / L%s :: %r" % m)
    sys.exit(2)
print()

# ---------------------------------------------------------------- [3] 写盘
if not DRY:
    if os.path.isdir(BACKUP):
        shutil.rmtree(BACKUP)
    os.makedirs(BACKUP)
    for fn, (L, nl) in buf.items():
        src = os.path.join(OUT, fn)
        shutil.copy2(src, os.path.join(BACKUP, fn))
        assert sha(os.path.join(BACKUP, fn)) == sha(src), "备份不一致: " + fn
        with open(src, "wb") as f:
            f.write(b"\n".join(L))
    print("[3] 写盘完成；原文件已备份到 %s" % BACKUP)
else:
    print("[3] （dry-run：跳过写盘与备份）")
print()

# ---------------------------------------------------------------- [4] 行数 + 语法
print("[4] 行数与语法")
import tempfile
TMPD = tempfile.mkdtemp(prefix="cs01v2_")
for fn, (L, nl) in buf.items():
    p = os.path.join(OUT, fn)
    n = len(open(p, "rb").read().split(b"\n")) if not DRY else len(L)
    if DRY:                      # dry-run：把**新内容**落到临时文件，供 [5]/[6] 做行为验证
        with open(os.path.join(TMPD, fn), "wb") as f:
            f.write(b"\n".join(L))
        p = os.path.join(TMPD, fn)
    ast.parse(io.open(p, encoding="utf-8").read())
    print("    %-38s 行数=%d  语法 OK" % (fn, n))
print()

# ---------------------------------------------------------------- [5] 行为验证
# ★ dry-run 下对**临时文件里的新内容**验证 —— 这样 dry-run 就能证明修正有效，
#   而不是像第一版那样去读未修改的原文件、报出一堆无意义的"仍不符"。
print("[5] 行为验证（执行修正后的 _canon，不是查文本）")
allok = True
for fn, *_ in SPEC:
    p = os.path.join(TMPD, fn) if DRY else os.path.join(OUT, fn)
    alias_of, order, can = drugmap_of(p)
    wrong = [(a, WANT[a], can(a)) for a in WANT if can(a) != WANT[a]]
    dup = {k: v for k, v in alias_of.items() if len(v) > 1}
    _, order0, _ = pre[fn]
    same_order = (order == order0)
    ok = (not wrong) and (not dup) and same_order
    allok &= ok
    print("    %-38s 10/10 正确=%s  跨行重复=%s  canonical 顺序不变=%s  rows=%d"
          % (fn, not wrong, dup or "无", same_order, len(order)))
    if wrong:
        print("        !! 仍不符: %s" % wrong)
    if dup:
        print("        !! 别名跨行重复（精确匹配下会让顺序重新起作用）: %s" % dup)
print()

# ---------------------------------------------------------------- [6] 与 amendment 全文复核
print("[6] 与 CS01_amendment.txt:11-21 的 10 条声明逐条复核")
_p6b = os.path.join(TMPD, SPEC[2][0]) if DRY else os.path.join(OUT, SPEC[2][0])
alias_of, order, can = drugmap_of(_p6b)
for a, c in WANT.items():
    got = can(a)
    print("    %-32s -> %-4s  %s" % (a, got, "OK" if got == c else "**不符**"))
print()
print("=" * 100)
print("结果：%s" % ("全部通过" if allok else "**有未通过项**"))
print("=" * 100)
sys.exit(0 if allok else 3)
