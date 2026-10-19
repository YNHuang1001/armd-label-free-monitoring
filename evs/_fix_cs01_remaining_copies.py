# -*- coding: utf-8 -*-
r"""_fix_cs01_remaining_copies.py —— 对**另外 5 份**药物映射副本执行 CS-01

★ 背景：CS-01 的「位置」栏列的是当时普查到的 5 个文件；全项目穷尽清点实为 12 份。
  另外 5 份仍是**子串匹配 + CS-01 之前的 42 别名表**（与已修的 5 份不一致）：
      p6d_01_unseal_evaluation.py        ← **在结局路径上**：p6d:85 d["abx"] = ...map(_canon)
      p12_configs_corrected.py           ← 在生产 17 个文件内（Fig 1E）
      p21a_reporting_profile_derivation.py ← 在生产 17 个文件内（reporting profiles）
      p26_analysis_unit_diagnostic.py    ← 生产树外，无引用
      p28_drug_name_audit.py             ← 生产树外，无引用
  **不纳入**：p4_01_task_availability.py（表是**第三支系**，含 LEX/CEP，与标准表差 13 行 —— 需单独裁定）

★ 本脚本对每份做两件事（= CS-01 对那 5 份做过的同一件事）：
    ① 5 行别名表 -> 标准表（含 10 个新别名，位置正确）
    ② `any(k in up for k in ks)` -> `up in ks`（废掉 substring-first-match）
  共 5 文件 x 6 行 = 30 处，全部行数守恒。

★ 验证（**行为验证 + 跨副本等价**，不是查文本）：
    · exec 每份修正后的 DRUGMAP + _canon
    · 对 **标准表 p6b 的全部别名** 逐一比对映射结果 —— 必须**完全一致**
    · 10 条 amendment 声明逐条正确
    · 无别名串跨行重复；canonical 顺序 == p6b 的顺序

用法：--dry-run / 直接跑
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
                      "PRODUCTION_backup", "CS01v2_remaining_20261006")
CANON = "p6b_01_corrected_pipeline.py"      # 标准表（已修 v2）
DRY = "--dry-run" in sys.argv

# (file, AMC行, TIMCZA行, GEN行, TOBAMKNIT行, CSTPEN行, _canon行) —— 全部实测
SPEC = [
    ("p6d_01_unseal_evaluation.py",            31, 35, 38, 39, 42, 47),
    ("p12_configs_corrected.py",               33, 37, 40, 41, 44, 49),
    ("p21a_reporting_profile_derivation.py",   68, 72, 75, 76, 79, 85),
    ("p26_analysis_unit_diagnostic.py",        47, 51, 54, 55, 58, 63),
    ("p28_drug_name_audit.py",                 40, 44, 47, 48, 51, 56),
]

OLD = {
    "AMC":       '    ("AMC", ["amoxicillinclavulanate", "amoxicillinclavulanicacid", "amoxclavulanate", "amoxicillinclav"]),',
    "TIMCZA":    '    ("TIM", ["ticarcillinclavulanate", "ticarcilclavulanate"]), ("CZA", ["ceftolozanetazobactam"]),',
    "GEN":       '    ("CIP", ["ciprofloxacin"]), ("LVX", ["levofloxacin"]), ("GEN", ["gentamicin"]),',
    "TOBAMKNIT": '    ("TOB", ["tobramycin"]), ("AMK", ["amikacin"]), ("NIT", ["nitrofurantoin"]),',
    "CSTPEN":    '    ("CST", ["colistin"]), ("SAM", ["sulbactam"]), ("AMP", ["ampicillin"]), ("PEN", ["piperacillin"]),',
    "CANON":     '        if any(k in up for k in ks): return c',
}
NEW = {
    "AMC":       '    ("AMC", ["amoxicillinclavulanate", "amoxicillinclavulanicacid", "amoxclavulanate", "amoxicillinclav", "amoxicillinkclavulanateaug"]),',
    "TIMCZA":    '    ("TIM", ["ticarcillinclavulanate", "ticarcilclavulanate", "ticarcillinclavulanicacid", "ticacarcillinclavulanicacid"]), ("CZA", ["ceftolozanetazobactam"]),',
    "GEN":       '    ("CIP", ["ciprofloxacin"]), ("LVX", ["levofloxacin"]), ("GEN", ["gentamicin", "gentscreen", "gentamicinscreen"]),',
    "TOBAMKNIT": '    ("TOB", ["tobramycin", "tobramycinrm"]), ("AMK", ["amikacin", "amikaciniv", "amikacinintravenous"]), ("NIT", ["nitrofurantoin"]),',
    "CSTPEN":    '    ("CST", ["colistin", "polymixine", "polymyxine"]), ("SAM", ["sulbactam"]), ("AMP", ["ampicillin"]), ("PEN", ["piperacillin"]),',
    "CANON":     '        if up in ks: return c',
}
ORDER = ["AMC", "TIMCZA", "GEN", "TOBAMKNIT", "CSTPEN", "CANON"]

WANT = {
    "amoxicillinkclavulanateaug": "AMC", "ticarcillinclavulanicacid": "TIM",
    "ticacarcillinclavulanicacid": "TIM", "gentscreen": "GEN", "gentamicinscreen": "GEN",
    "tobramycinrm": "TOB", "amikaciniv": "AMK", "amikacinintravenous": "AMK",
    "polymixine": "CST", "polymyxine": "CST",
}


def sha(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()


def dm_of(path):
    txt = io.open(path, encoding="utf-8").read()
    blk = re.search(r"^DRUGMAP = \[(.*?)^\]", txt, re.S | re.M).group(1)
    fn_src = re.search(r"^def _canon\(n\):(?:\n[ \t].*)+", txt, re.M).group(0)
    ns = {}
    exec("DRUGMAP = [" + blk + "]", ns)
    exec(fn_src, ns)
    alias_of, order = {}, []
    for c, ks in ns["DRUGMAP"]:
        if c not in order:
            order.append(c)
        for k in ks:
            alias_of.setdefault(k, []).append(c)
    return alias_of, order, ns["_canon"]


print("=" * 100)
print("对另外 5 份执行 CS-01  %s" % ("【DRY-RUN】" if DRY else "【真实写入】"))
print("=" * 100)

canon_alias, canon_order, canon_fn = dm_of(os.path.join(OUT, CANON))
print("标准表 %s：别名 %d 个 / canonical %d 个" % (CANON, len(canon_alias), len(canon_order)))
print()

print("[1] 修正前：各副本与标准表的映射分歧")
for fn, *_ in SPEC:
    a, o, f = dm_of(os.path.join(OUT, fn))
    diff = [k for k in set(list(a) + list(canon_alias)) if f(k) != canon_fn(k)]
    print("    %-40s 别名=%-3d 与标准表分歧别名=%d" % (fn, len(a), len(diff)))
print()

print("[2] 逐行字节级断言并替换（%d 文件 x 6 行）" % len(SPEC))
buf, miss = {}, []
for fn, *lns in SPEC:
    raw = open(os.path.join(OUT, fn), "rb").read()
    nl = b"\r" if b"\r\n" in raw else b""
    L = raw.split(b"\n")
    for key, ln in zip(ORDER, lns):
        o = OLD[key].encode("utf-8") + nl
        if not (1 <= ln <= len(L)) or L[ln - 1] != o:
            miss.append((fn, key, ln, (L[ln - 1][:85] if 1 <= ln <= len(L) else b"OOB")))
            continue
        L[ln - 1] = NEW[key].encode("utf-8") + nl
    buf[fn] = (L, nl)
    print("    %-40s 6 行 %s" % (fn, "命中" if not any(m[0] == fn for m in miss) else "**未命中**"))
if miss:
    print("\n  !! 未命中 %d 处 —— 中止，不写盘：" % len(miss))
    for m in miss:
        print("     %s / %s / L%s :: %r" % m)
    sys.exit(2)
print()

import tempfile
TMPD = tempfile.mkdtemp(prefix="cs01rem_")
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
    print("[3] 写盘完成；备份到 %s" % BACKUP)
else:
    for fn, (L, nl) in buf.items():
        with open(os.path.join(TMPD, fn), "wb") as f:
            f.write(b"\n".join(L))
    print("[3] （dry-run：新内容写到临时目录供验证，未动生产）")
print()

print("[4] 行数与语法")
for fn, *lns in SPEC:
    p = os.path.join(TMPD, fn) if DRY else os.path.join(OUT, fn)
    n = len(open(p, "rb").read().split(b"\n"))
    ast.parse(io.open(p, encoding="utf-8").read())
    print("    %-40s 行数=%d  语法 OK" % (fn, n))
print()

print("[5] 行为验证 + 跨副本等价（执行真 _canon；对**标准表全部别名**比对）")
allok = True
for fn, *lns in SPEC:
    p = os.path.join(TMPD, fn) if DRY else os.path.join(OUT, fn)
    a, o, f = dm_of(p)
    wrong10 = [(k, WANT[k], f(k)) for k in WANT if f(k) != WANT[k]]
    disagree = [k for k in set(list(a) + list(canon_alias)) if f(k) != canon_fn(k)]
    dup = {k: v for k, v in a.items() if len(v) > 1}
    ok = (not wrong10) and (not disagree) and (not dup) and (o == canon_order)
    allok &= ok
    print("    %-40s 10/10=%s  与标准表分歧=%d  跨行重复=%s  顺序同标准=%s"
          % (fn, not wrong10, len(disagree), dup or "无", o == canon_order))
    if wrong10:
        print("        !! 10 条声明仍不符: %s" % wrong10)
    if disagree:
        print("        !! 与标准表分歧的别名（前 8）: %s" % disagree[:8])
print()
print("=" * 100)
print("结果：%s" % ("全部通过 —— 10 份副本的药物映射语义完全一致" if allok else "**有未通过项**"))
print("=" * 100)
sys.exit(0 if allok else 3)
