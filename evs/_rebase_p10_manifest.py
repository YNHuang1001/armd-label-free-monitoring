# -*- coding: utf-8 -*-
r"""_rebase_p10_manifest.py —— 把 ACTIVE_PATH_MANIFEST.yaml 里 p10a / p10b 的**行号锚点**平移

与 `_rebase_p10_line_refs.py` 同一张位移表（由 difflib 逐段算出、并与 7 个按内容定位的块交叉验证）：
    p10a : 旧 >=35 -> +18 ；>=333 -> +19 ；>=344 -> +32
    p10b : 旧 31..61 -> +13 ；62..201 -> +21 ；202..203 -> +30 ；204..206 -> +46 ；207..211 -> +47 ；>=212 -> +49

★ 只改锚点里的行号与 `p10a:NNN` / `p10b:NNN` 形式的文字引用；
  **不碰** p5/p6b/p6c/p6d/p10c 等其它文件的引用（它们的行号没变）。
★ 改完用 t14 的 MC-02（锚点可按文本解析）复核；锚点还带 `contains` 正则，错位会被它抓出来。
"""
import io
import os
import re
import sys

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

HERE = os.path.dirname(os.path.abspath(__file__))
PROJ = os.path.dirname(HERE)
MAN = os.path.join(PROJ, "ACTIVE_PATH_MANIFEST.yaml")
DRY = "--dry-run" in sys.argv


def sa(n):
    return n + 32 if n >= 344 else (n + 19 if n >= 333 else (n + 18 if n >= 35 else n))


def sb(n):
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


SH = {"p10a": sa, "p10b": sb}
NUMS = {"p10a": "p10a_prepare_cache.py", "p10b": "p10b_subsample_curves.py"}

src = io.open(MAN, encoding="utf-8").read()
orig = src
log = []

# ---- ① 锚点：{file: p10X_*.py, kind: lines, value: 'A-B', ...}
def fix_anchor(m):
    fn, a, b = m.group(1), int(m.group(2)), int(m.group(3))
    key = "p10a" if fn.startswith("p10a") else ("p10b" if fn.startswith("p10b") else None)
    if key is None:
        return m.group(0)
    na, nb = SH[key](a), SH[key](b)
    if (na, nb) == (a, b):
        return m.group(0)
    log.append(("anchor", fn, "%d-%d" % (a, b), "%d-%d" % (na, nb)))
    return m.group(0).replace("value: '%d-%d'" % (a, b), "value: '%d-%d'" % (na, nb))


src = re.sub(r"\{file: (p10[ab]_[a-z_]+\.py), kind: lines, value: '(\d+)-(\d+)'", fix_anchor, src)

# ---- ② 文字引用：p10a:NNN  /  p10b:NNN-MMM
def fix_text(m):
    key, a, b = m.group(1), int(m.group(2)), (int(m.group(3)) if m.group(3) else None)
    na = SH[key](a)
    nb = SH[key](b) if b is not None else None
    if na == a and (b is None or nb == b):
        return m.group(0)
    log.append(("text", key, m.group(0), "%s:%d%s" % (key, na, ("-%d" % nb) if nb else "")))
    return "%s:%d%s" % (key, na, ("-%d" % nb) if nb else "")


src = re.sub(r"\b(p10[ab]):(\d+)(?:-(\d+))?", fix_text, src)

print("=" * 100)
print("平移 manifest 里 p10a/p10b 的行号引用  %s" % ("【DRY-RUN】" if DRY else "【真实写入】"))
print("=" * 100)
for kind, fn, old, new in log:
    print("  [%-6s] %-28s %-12s -> %s" % (kind, fn, old, new))
print("\n  合计 %d 处" % len(log))
if not DRY:
    io.open(MAN, "w", encoding="utf-8", newline="\n").write(src)
    print("  已写入 %s" % MAN)
