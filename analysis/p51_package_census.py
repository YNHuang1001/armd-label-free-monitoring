# -*- coding: utf-8 -*-
r"""p51 — 包级**双向**普查 ＋ **闭包**检查（A4.3）

================================ 为什么是「双向」================================
**单向普查看不见「该补的行没补上」**（外部核验）。故两向都查：
  · **旧指纹全缺** —— 陈旧值不得残留在任何期刊件里
  · **新值按位就位** —— 重建后的值必须真的在

================================ 两种验证模式（外部核验给定，采纳）================
  · **千分位数**（1,372,185 之流）→ **指纹扫描**（文本中可区分）
  · **小整数**（23/38/42/76/114 之流）→ **只能重算**（扫了等于没扫）⇒ 归 p47

================================ ★ 闭包检查（GLM 结构注记，2026-10-08）============
> **p44 是 F 节漂移的形状** —— 一份**手写名单**，与世界静默分叉。
> F 节换血是它，代码包缺件也是它。**修法是补录不是机制** ——
> 下一个会话第 114 个脚本落地时，同样的静默排除会再来一次，
> **除非有东西在「名单外的文件」上开火。**

故本模块的核心不是指纹，是**闭包**：
  · **树上每个文件必须被账到**（收录 or 显式排除 ＋ 理由）
  · **名单上每项必须存在**（悬空指针）
**114 撞击与这个闭包，是同一扇门的两面。**

================================ ★ 114 的同形异指（停放观察 O2）================
C 节把 `114` 列为**陈旧行数**（→111），但 `114` 同时是 **(组 × 站) 组合数**——
一个**活不变量**（`{2:114}`、`universe_partition_audit.csv` 114 行）。
指纹扫描必撞，故本模块**按上下文分类**而非一律判红。
"""
import csv
import io
import os
import re
import sys

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
AN = os.path.join(ROOT, "分析结果")
PKG = os.path.join(ROOT, "SUBMISSION_PACKAGE_20261007")
REP = open(os.path.join(AN, "report_p51_package_census.txt"), "w", encoding="utf-8")


def p(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True); REP.write(s + "\n"); REP.flush()


# ---------------------------------------------------------------- 扫描目标（期刊件）
def journal_files():
    out = []
    for base, _, fs in os.walk(PKG):
        for f in fs:
            if f.lower().endswith((".md", ".csv", ".txt", ".json")):
                out.append(os.path.join(base, f))
    return sorted(out)


# ---------------------------------------------------------------- 旧指纹 / 新值
OLD = {
    # ★ 2026-10-08 更正（本模块首跑即自证）：**裸小整数不做文本扫描** ——
    #   GLM 的规则「小整数扫了等于没扫」，我在文件头写了，却在代码里违反：
    #   `74/68/60/64/67/71` 在 SHA 串、图注小数、坐标里到处都是 ⇒ 首跑 7 类全假阳性。
    #   它们的正确检查方式是**锚定到 R 表的行结构**（见 RLOW 块），不是全文扫描。
    # **202 亦非陈旧指纹**：命中处一是「已被取代的 MMD 实现…200 of its 202 tasks…identified as
    #   invalid and superseded（历史披露，**必须保留**）」，二是冻结列名 `pooled_202task` /
    #   `pooled_raw_all202`（改列名会改冻结产物）。故移出通缉名单，改由下方 SUPERSEDED 白名单确认。
    "1,372,185（task-expanded 旧值）": r"1,372,185",
    "818,466（task-expanded 旧值）":   r"818,466",
    "630,006（task-expanded 旧值）":   r"630,006",
    "2,820,657（旧合计）":             r"2,820,657",
}
# 小整数旧指纹**不作为文本扫描对象**（扫了等于没扫）——见文件头「两种验证模式」
SMALL_INT_STALE = ["202", "13（旧空行数）", "16（幻影差额）", "21（S1 旧列数）"]

NEW = {
    "70 / 62 / 54（as source）": r"70.*62.*54|70 \| 62 \| 54",
    "58 / 62 / 66（as target）": r"58 \| 62 \| 66",
    "6 / 14 / 22（排除·按站）":   r"6 / 14 / 22",
    "4 / 12 / 20（净去 KP_NIT）": r"4 / 12 / 20",
    "1,295,777（task-expanded）": r"1,295,777",
    "228 − 42 = 186（桥）":        r"228 − 42 = 186",
    "candidate grid（分母标签）":  r"candidate grid",
    "26（S1 列数上限）":           r"26 columns",
}

# ---------------------------------------------------------------- R 表行结构锚（替代裸小整数扫描）
# 小整数的正确判据：**该行是否还印着旧值**。行结构固定，故可精确锚定。
RLOW = [
    ("R1/R2 as source 应已改 202-口径", r"as source\s*\|?\s*\*{0,2}\s*\|?\s*7[04]",
     r"^\| Eligible directed transfers[^|]*as source \| .*(70|62|54)"),
    ("R1/R2 as target 应已改 202-口径", None,
     r"^\| Eligible directed transfers[^|]*as target \| .*(58|62|66)"),
]
STALE_ROWS = [
    # R 表行：`| ... as source | 74 | ...` 里的旧站值（新值应为 70/62/54）
    "as source[^\n]*\b74\b",
    "as target[^\n]*\b64\b",
    # 陈旧行数表述（新表述为 111 行；`114 rows = 38 groups x 3 sites` 属不变量义，见 ④）
    "full 114-row breakdown",
]
# 「已被取代」的白名单：出现这些词即属**依法披露**，不是陈旧残留
SUPERSEDED = ["superseded", "earlier MMD implementation", "was identified as invalid",
              "pooled_202task", "pooled_raw_all202", "p6_superseded_pre_correction"]

p("=" * 108)
p("p51 — 包级双向普查 ＋ 闭包检查（A4.3）")
p("=" * 108)
FILES = journal_files()
p(f"扫描目标：{PKG} 下 {len(FILES)} 个期刊件（.md/.csv/.txt/.json）")
p("")

# ---------------------------------------------------------------- ① 旧指纹（应全缺）
p("①a **R 表行结构锚定**（小整数的正确判据：该行是否还印着旧值）")
for pat in STALE_ROWS:
    hits = []
    for f in FILES:
        txt = io.open(f, encoding="utf-8", errors="ignore").read()
        for m in re.finditer(pat, txt, re.M):
            seg = txt[max(0, m.start() - 80):m.end() + 80]
            if any(k in seg for k in SUPERSEDED):
                continue
            hits.append((os.path.relpath(f, PKG), seg.strip()[:80]))
    p(f"    {pat[:34]:<36}{'✔ 无陈旧行' if not hits else '✘ ' + str(hits[:2])}")
    if hits:
        bad_old.append((pat, hits))

p("")
p("①b 旧指纹（**仅千分位数**做文本扫描；小整数见 ①a）—— **应全缺**")
bad_old = []
for name, pat in OLD.items():
    hits = []
    for f in FILES:
        try:
            t = io.open(f, encoding="utf-8", errors="ignore").read()
        except Exception:
            continue
        if re.search(pat, t):
            hits.append(os.path.relpath(f, PKG))
    if hits:
        bad_old.append((name, hits))
    p(f"    {name:<34}{'✘ **命中 ' + str(len(hits)) + ' 个文件**' if hits else '✔ 全缺'}"
      + (f"  → {hits[:2]}" if hits else ""))

# ---------------------------------------------------------------- ② 新值（应按位就位）
p("")
p("② 新值 —— **应按位就位**")
bad_new = []
for name, pat in NEW.items():
    hits = [os.path.relpath(f, PKG) for f in FILES
            if re.search(pat, io.open(f, encoding="utf-8", errors="ignore").read())]
    if not hits:
        bad_new.append(name)
    p(f"    {name:<34}{'✔ 就位（' + str(len(hits)) + ' 件）' if hits else '✘ **缺席**'}")

# ---------------------------------------------------------------- ③ ★ 闭包检查
p("")
p("③ ★ **闭包检查**（GLM 结构注记：树上每个文件必须被账到；名单上每项必须存在）")
bad_closure = []

# 3a. p44 的收录 ∪ 排除 是否**划分**了 分析结果 下的全部 .py
# ★ 2026-10-08 更正：初版 `import p44_package_code` —— 那是**脚本不是模块**，
#   import 会**执行整个打包**（副作用：打印 + 重写代码包）。改为**文本解析**名单。
_SRC44 = io.open(os.path.join(AN, "p44_package_code.py"), encoding="utf-8").read()


def _names(listname):
    m = re.search(r"^%s = \[(.*?)^\]" % listname, _SRC44, re.S | re.M)
    if not m:
        return set()
    return set(re.findall(r'\("([^"]+\.py)"', m.group(1)))


LISTED = _names("CHAIN") | _names("AUDIT") | _names("RENDER")
THROWAWAY_RE = re.compile(r"^_")
_OTHER_M = re.search(r'^OTHER_PROJECT = r"([^"]+)"', _SRC44, re.M)
OTHER_RE = re.compile(_OTHER_M.group(1) + r"|^litq\.py$|^ms_master\.py$|^style_map\.py$") if _OTHER_M else None
ALL_PY = {f for f in os.listdir(AN) if f.endswith(".py")}
UNLISTED = sorted(ALL_PY - LISTED)
p(f"    3a. `分析结果/` 下 .py 共 {len(ALL_PY)} 个；p44 名单 {len(LISTED)} 项")
p(f"        **名单外（必须各有排除理由）**：{len(UNLISTED)} 个")
for f in UNLISTED:
    why = ("一次性修补（_ 前缀）" if THROWAWAY_RE.match(f)
           else "另一课题" if (OTHER_RE and OTHER_RE.match(f))
           else "**⚠ 无理由**")
    if why == "**⚠ 无理由**":
        bad_closure.append(("名单外且无理由", f))
    p(f"          {f:<46}{why}")

# 3b. 名单上每项必须存在（悬空指针）
p("    3b. 名单**悬空**检查（列了但文件不存在）")
dangling = sorted(fn for fn in LISTED if not os.path.exists(os.path.join(AN, fn)))
for f in dangling:
    bad_closure.append(("悬空", f))
p(f"        {'✔ 无悬空' if not dangling else '✘ ' + str(dangling)}")

# 3c. 包内 MANIFEST 每项必须存在（同上，在包侧）
p("    3c. 包内 `08_code/MANIFEST.csv` 的每项必须存在")
MAN = os.path.join(PKG, "08_code", "MANIFEST.csv")
if os.path.exists(MAN):
    rows = list(csv.DictReader(io.open(MAN, encoding="utf-8")))
    key = "rel" if "rel" in (rows[0] if rows else {}) else list(rows[0])[0]
    miss = [r[key] for r in rows if not os.path.exists(os.path.join(PKG, "08_code", r[key]))]
    for m in miss:
        bad_closure.append(("MANIFEST 悬空", m))
    p(f"        MANIFEST {len(rows)} 项 ｜ 缺失 {len(miss)} 项"
      + (f" → {miss[:3]}" if miss else "  ✔"))
else:
    p("        ⚠ 无 MANIFEST.csv")

# ---------------------------------------------------------------- ④ 114 的同形异指
p("")
p("④ **114 的同形异指**（停放观察 O2：指纹扫描必撞活不变量）")
# ★ 2026-10-08 更正：`114 rows = 38 groups × 3 sites` 是**不变量义**，初版被误判为行数义。
#   行数义的关键词收紧为「陈旧行数」的明确表述。
CTX_ROW = re.compile(r"full 114|114[- ]row breakdown|114 行明细")
CTX_INV = re.compile(r"\{2: ?114\}|114 行 =|38 × 3|38 组 × 3|组合数")
n_row = n_inv = 0
for f in FILES:
    t = io.open(f, encoding="utf-8", errors="ignore").read()
    for m in re.finditer(r"(?<!\d)114(?!\d)", t):
        seg = t[max(0, m.start() - 60):m.end() + 60]
        if CTX_ROW.search(seg):
            n_row += 1
        elif CTX_INV.search(seg):
            n_inv += 1
p(f"    以「行数」义出现：**{n_row}**（陈旧指纹，应为 0 —— 已改为 111）")
p(f"    以「不变量」义出现（组 × 站）：{n_inv}（**活值，不得判红**）")
p("    ⇒ 判据：**按上下文分类，不按字面一律判红** —— 这是 O2 预登记的撞击点")

# ---------------------------------------------------------------- 判词
p("")
p("=" * 108)
n_bad = len(bad_old) + len(bad_new) + len(bad_closure) + (1 if n_row else 0)
p(f"① 旧指纹命中 {len(bad_old)} 类 ｜ ② 新值缺席 {len(bad_new)} 项 ｜ "
  f"③ 闭包缺口 {len(bad_closure)} 项 ｜ ④ 114 行数义 {n_row} 处")
if bad_closure:
    for k, f in bad_closure[:8]:
        p(f"    闭包：{k} → {f}")
p("")
p("*** p51 " + ("通过 —— 包级双向普查与闭包检查全部成立 ***" if n_bad == 0
                else f"**未通过（{n_bad} 处）** ***"))
REP.close()
sys.exit(0 if n_bad == 0 else 1)
