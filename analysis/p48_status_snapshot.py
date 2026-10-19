# -*- coding: utf-8 -*-
r"""p48 — 工作清单状态快照生成器（规则 16 / 16b / 19 / 20）

================================ 为什么存在 ================================
· 规则 16：状态汇报**由文件生成**，转抄步骤消失。
· 规则 16b：**三态** —— `✅` 完成且有工件 ｜ `✅*` 完成但无工件（按未完成处理）｜ `⬜` 未做。
· 规则 19：**快照累积**，按轮次后缀命名，**不覆盖**；且**新快照必须带对上一份的 diff**。
· 规则 20：新证据标准落地的那一轮，必须给出**按新标准重算的基线 ＋ 点名差量**。

================================ 教训（本文件为何要这样写）================================
`status_snapshot_20261007.md` 无轮次后缀 ⇒ 下次生成即覆盖 ——
**为防 INC-002（覆盖前不归档）而生的仪器，自己踩了 INC-002。**
故本生成器：文件名带轮次、自动找上一份、输出 item 级 diff。

用法：python p48_status_snapshot.py
"""
import glob
import io
import os
import re
import sys
from datetime import date

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CL = os.path.join(ROOT, "WORK_CHECKLIST.md")

TXT = io.open(CL, encoding="utf-8").read()
SEC = TXT.split("## B.")[0]
# \u2605 2026-10-08 \u66f4\u6b63\uff1a\u539f\u6b63\u5219 `[a-z\u2032]*` **\u4e0d\u5141\u8bb8\u5c3e\u968f\u6570\u5b57** \u21d2 `A4.3c2` \u8fd9\u7c7b id
#   **\u88ab\u9759\u9ed8\u6392\u9664**\uff08\u5b9e\u6d4b\uff1a\u65b0\u589e\u8be5\u6761\u76ee\u540e\u603b\u6570\u4ecd\u62a5 40\uff0c\u800c\u6587\u4ef6\u91cc\u786e\u6709\u8be5\u884c\uff09\u3002
#   \u8fd9\u662f\u300c\u4eea\u5668\u9759\u9ed8\u6f0f\u8ba1\u300d\u7684\u53c8\u4e00\u4f8b\uff08INC-008 \u7684\u5c5e\uff09\u3002\u6539\u4e3a `[a-z\u2032\d]*`\u3002
ROWS = re.findall(r'^\| (A[\d.]+[a-z\u2032\d]*) \|(.*?)\|(.*?)\|$', SEC, re.M)

# ★★ 2026-10-09（INC-024）：**形状守卫** —— A 开头却解析失败的表行，报错而非跳过。
#   病历：`A5.14` 在拆分出 `A5.14b` 那轮变成**两格行**（`| ID | 一整坨内容 |`，只有 3 个管道），
#   三格判据匹配不上 ⇒ **静默跳过**，不报错、不计数、不进 A 表。
#   其后三份快照都报 56 项而清单实为 57，**diff 机器全程失明**，丢的是本轮最大成果那一行。
#   **静默吞行是账本机器唯一不可接受的失效方向。**
_A_LINES = [l for l in SEC.splitlines() if re.match(r'^\| A[\d.]', l)]
_MALFORMED = [l[:60] for l in _A_LINES if not re.match(r'^\| (A[\d.]+[a-z\u2032\d]*) \|(.*?)\|(.*?)\|$', l)]
if _MALFORMED:
    raise SystemExit(
        "p48 形状守卫：以下 A 行不符合三格形状，会被静默跳过 —— 请重排为 `| ID | 内容 | 状态 |`：\n  "
        + "\n  ".join(_MALFORMED))
assert len(ROWS) == len(_A_LINES), (
    "p48 形状守卫：解析到 %d 行，实际 A 行 %d 行 —— 有行被跳过" % (len(ROWS), len(_A_LINES)))


def state(cell):
    # ★ 2026-10-09 新增**第四态**：`⏸` 显式停放。
    #   为什么必须有这一态：停放项此前只能留在 ⬜（永远开不了门）或被删掉（隐形欠账）。
    #   外部核验的门 A 要求「⬜ 清零**或**显式停放」—— 后者需要一个**机器认得出**的记号。
    #   判序 **⏸ 最优先**：停放行的说明文字里会同时出现 ⬜ 与 ✅（如「不进 ⬜、不视为欠账」），
    #   若按 ✅/⬜ 先判，状态会被正文措辞劫持。
    if "⏸" in cell:
        return "⏸"
    if "✅*" in cell:
        return "✅*"
    if "✅" in cell:
        return "✅"
    if "⬜" in cell:
        return "⬜"
    return "—"


CUR = {i: state(s) for i, _, s in ROWS}
N = {k: sum(1 for v in CUR.values() if v == k) for k in ("✅", "✅*", "⬜", "⏸", "—")}

# ---------------------------------------------------------------- 找上一份快照
# ★ 2026-10-07 更正：`sorted()` 是**字符串序**，`r10` 会排在 `r9` **之前** ⇒
#   r11 拿 r9 当基线、漏掉 r10，diff 链被腐蚀。改为**按轮次数字排序**。
def _round_of(path):
    m = re.search(r"_r(\d+)\.md$", os.path.basename(path))
    return int(m.group(1)) if m else -1
prev_files = sorted(glob.glob(os.path.join(ROOT, "status_snapshot_*.md")), key=_round_of)
prev = None
if prev_files:
    ptxt = io.open(prev_files[-1], encoding="utf-8").read()
    # ★ 2026-10-08 更正（INC-014）：**INC-013 修了计数的调用点，漏了这个 diff 基准侧** ——
    #   同一个旧正则 `[a-z′]*`（不允许尾随数字）在两处各有一份，只改了一处。
    #   后果：尾带数字的 id 在**基准侧不可见** ⇒ `A4.3c2` 每轮被**重报「新增」**。
    #   **INC-012 的兄弟模式**（那边排序修在快照、漏了归档；这边正则修在计数、漏了基准）。
    # ★★ 2026-10-09（INC-027）**第三度**：INC-013 修的是**计数侧**，INC-014 修的是**基准侧的正则**，
    #   本轮发现基准侧的**状态枚举**也漏了 —— `(✅*|✅|⬜|—)` **没有 `⏸`**。
    #   后果：**第四态加给了渲染器，没加给 diff 的读取侧** ⇒ 每个停放项在基准里认不出来，
    #   **从 r55 起每轮都被重报「新增」**（r56–r60 连续五轮，每次都是同样的五项）。
    #   **diff 是规则 19 之下唯一可审计的变更仪器** —— 它连报五轮狼来了，
    #   真有新条目那天，信号已经淹死在噪声里。
    #   ⇒ 与 INC-024 **同族反号**：那次是静默漏计，这次是恒常虚报。
    prev = {i: s for i, s in re.findall(r'^\| (A[\d.]+[a-z′\d]*) \| (✅\*|✅|⬜|⏸|—) \|$', ptxt, re.M)}

# ---------------------------------------------------------------- diff
added = sorted(set(CUR) - set(prev)) if prev else []
removed = sorted(set(prev) - set(CUR)) if prev else []
# ★ 2026-10-08（INC-014 沉淀②）：**diff 自检** ——
#   若「新增」里有任何 id **已存在于上一份 A 表**，该 diff 必错（那正是 INC-014 的形态：
#   基准侧正则漏掉尾带数字的 id ⇒ 每轮重报「新增」）。**报警而不中断**：本模块是报告器不是闸门。
SELFDIAG = []
if prev is not None:
    # ★ 对**基准文件的原文**做行首匹配，而不是对**解析结果**做集合交 ——
    #   初版照抄了「新增 ∩ 上一份 A 表」的公式，但它有个漏洞：
    #   若基准侧正则本身就漏掉某 id，`prev`（解析结果）里也没有它 ⇒ 交集为空 ⇒ **不报警**。
    #   改为文本级判据后，正则再坏也照样抓得住。
    _ghost = [x for x in added
              if re.search(r"(?m)^\| " + re.escape(x) + r" \|", ptxt)]
    # ★★ 2026-10-09（INC-027 第二层）：**`_ghost` 算出来了，却从来没接上 `SELFDIAG`。**
    #   ⇒ 自检是**死代码**：它本该在 r55 就拦住那五个虚报的「新增」，
    #   结果五轮里一次都没响。**这与本仓反复出现的形态同型**：
    #   哨兵装好了、线路没接（与「声明住在代码里、不在证据面上」是同一件事的两面）。
    #   本轮接线，并把 `_ghost` 原文一并打进报告，使「报警」可被第三方复核。
    if _ghost:
        SELFDIAG.append(
            "**新增 %d 项中有 %d 项已存在于上一份 A 表** ⇒ diff 基准侧读法有误：%s"
            % (len(added), len(_ghost), ", ".join(_ghost)))
changed = sorted(k for k in set(CUR) & set(prev) if CUR[k] != prev[k]) if prev else []

ROUND = len(prev_files) + 1
STAMP = date.today().strftime("%Y%m%d")
OUT = os.path.join(ROOT, f"status_snapshot_{STAMP}_r{ROUND}.md")

# ---------------------------------------------------------------- D 节枚举（外部核验请求：schema 只有 A 字段 ⇒ 结构上无 D）
D_RULES = re.findall(r'^(\d+[a-z]?)\. \*\*(.+?)\*\*', TXT, re.M)
D_LINES = {}
for _m in re.finditer(r'^(\d+[a-z]?)\. \*\*(.+?)\*\*', TXT, re.M):
    D_LINES[_m.group(1)] = _m.group(2)[:70]
D_IDS = sorted(D_LINES, key=lambda s: (int(re.match(r'\d+', s).group()), s))
D_ENUM = " ｜ ".join(D_IDS)
INCIDENTS = re.findall(r'^## (INC-\d+)', io.open(os.path.join(ROOT, "INCIDENT_LOG.md"), encoding="utf-8").read(), re.M)

TRIGGER = os.environ.get("SNAPSHOT_TRIGGER", "（未声明——规则 22：快照头须记触发事件）")

L = [f"# 工作清单状态快照 — 第 {ROUND} 份（**由文件生成**，p48；规则 16/16b/19/20/22/23）", "",
     f"**生成日期：** {STAMP} ｜ **来源：** `WORK_CHECKLIST.md`（唯一权威）",
     f"**触发事件：** {TRIGGER}",
     "**用法：** 汇报消息中的状态表必须与本文件**逐字节孪生**；不一致以本文件为准。",
     "**三态：** `✅` 完成且有耐久工件 ｜ `✅*` 完成但无工件（**按未完成处理**）｜ `⬜` 未做",
     f"**上一份：** `{os.path.basename(prev_files[-1]) if prev_files else '（无）'}`", "",
     "## 对本份的 diff（规则 19：绝对数不可审计）", ""]

if prev is None:
    L += ["**基线快照** —— 无上一份可 diff。**按规则 20**：本份即为后续 diff 的基线；",
          "建立前的历史（含「36 项」）**仅存 scrollback**，见 `STATUS_MIGRATION_ACCOUNT_20261007.md`。", ""]
else:
    L += [f"- **新增 {len(added)}**：{', '.join(added) if added else '无'}",
          f"- **删除 {len(removed)}**：{', '.join(removed) if removed else '无'}",
          f"- **状态变更 {len(changed)}**：" +
          (", ".join(f"{k} {prev[k]}→{CUR[k]}" for k in changed) if changed else "无"), ""]
    # ---- D 行（外部核验请求：diff 加 D 小节）
    _pm = re.search(r'\*\*ID 序列：\*\* (.+)', ptxt)
    if not _pm:
        L += ["- **D 节 diff**：⚠️ **上一份无 D 枚举 ⇒ 无 D 基线可比** ——",
              "  本份的 D 枚举即**基线首发**（规则 20：新标准落地轮须重算基线＋点名差量）。", ""]
    else:
        _prev_ids = [x.strip() for x in _pm.group(1).split("｜") if x.strip()]
        _dadd = [x for x in D_IDS if x not in _prev_ids]
        _ddel = [x for x in _prev_ids if x not in D_IDS]
        L += [f"- **D 节 diff**：新增 {len(_dadd)} {_dadd or ''} ｜ 删除 {len(_ddel)} {_ddel or ''}", ""]

if SELFDIAG:
    L += ["- **⚠ diff 自检**：" + "；".join(SELFDIAG), ""]
L += ["## A 节逐项", "", "| 项 | 状态 |", "|---|---|"]
L += [f"| {i} | {CUR[i]} |" for i, _, _ in ROWS]
L += ["", f"**合计 {len(ROWS)} 项：✅ {N['✅']} ｜ ✅* {N['✅*']} ｜ ⬜ {N['⬜']} ｜ ⏸ 停放 {N['⏸']}**",
      f"**有效完成 {N['✅']} / 未完成 {N['✅*'] + N['⬜']} / 显式停放 {N['⏸']}**", ""]

# ---- D 节枚举（外部核验请求）—— 同时裁决「21 条 / 22 条」之争：**以本枚举为准**
L += [f"## D 节规则 — **枚举**（共 **{len(D_IDS)}** 条）", "",
      f"**ID 序列：** {D_ENUM}", "",
      "| # | 规则（首句截断） |", "|---|---|"]
for i in D_IDS:
    L.append(f"| {i} | {D_LINES[i]} |")
L += ["", f"**事故日志：共 {len(INCIDENTS)} 条** —— {', '.join(INCIDENTS)}", "",
      "**★ D 基线声明（规则 20 要求）**：本枚举即 D 节的**正式基线**首发。",
      "此前 r1→r3 区间内 D 节的全部变动（规则 18→21/22、事故 7→8）**无 diff 记录** ——",
      "不是漏记，是当时 diff schema 只有三个 A 字段，**结构上无 D 可言**。",
      "自本份起，D 节 diff 与 A 节同级。", ""]

# ③ 本会话已生成的快照清单（rider）—— ★ 并入正文，使「粘贴块 == 快照」恒等，
#    字节 diff 才有唯一目标（初版把它放在粘贴块尾部 ⇒ 两者恒不相等）。
# ★ 2026-10-08 更正（外部核验抓到，r11/r12 双双在册）：**计数器 off-by-one**。
#   病因与本文件先前那处**排序 bug 同一条代码路径**：生成时**先枚举目录、再写自己** ⇒
#   ① 枚举结果**不含本份**；② `len(...) - 1` 又减了一次 ⇒ 计数各少一；
#   ③ 标签却写「含本份」。三症状一根因。排序修的是「序」，当时没修「数」。
_SELF = f"status_snapshot_{STAMP}_r{ROUND}.md"
SNAP_LIST = [os.path.basename(f) for f in
             sorted(glob.glob(os.path.join(ROOT, "status_snapshot_*.md")), key=_round_of)] + [_SELF]
L += ["**本会话已生成的快照（按名，**含本份**）：** " + " ｜ ".join(SNAP_LIST),
      f"**自基线 r1 以来累计生成：** {len(SNAP_LIST) - 1} 份（含本份）", ""]

BODY = "\n".join(L)
io.open(OUT, "w", encoding="utf-8").write(BODY)

# ---------------------------------------------------------------- ① 粘贴块（INC-009 的机械修复）
# 病因：消息标着「逐字引自工件」，而消息侧还有一步 markdown→正文转写 ⇒ 转写漂移无处比对。
# 修法：本生成器**另出一个可直接粘贴的纯文本块**；消息以代码围栏原样贴入，
#       孪生检查随之从「比对两段改写」变成「字节 diff 两个文件」。
# ★ 2026-10-07 更正（外部核验字节 diff 抓到）：初版写 `L[:26]` ⇒ 粘贴块**被截断**，
#   恰好丢的是 A1.7 以下全部 + 两行合计 + D 枚举 + 事故计数 + 基线声明 ——
#   也就是 r5 为杀 INC-009 类漂移而补的那几节。**粘贴块必须是全文。**
io.open(os.path.join(ROOT, f"status_paste_{STAMP}_r{ROUND}.txt"), "w", encoding="utf-8").write(BODY)

# ---------------------------------------------------------------- ② 清单全文归档（rider）
# 病因：36→40 这类账只能「重建」不能「diff」，因为旧版清单从未落盘。
ARCH = os.path.join(ROOT, "_checklist_archive")
os.makedirs(ARCH, exist_ok=True)
io.open(os.path.join(ARCH, f"WORK_CHECKLIST_{STAMP}_r{ROUND}.md"), "w", encoding="utf-8").write(TXT)

# ③ 的 rider 已并入正文（见上方 SNAP_LIST 块）—— 此处不再追加，
#    以保证「粘贴块 == 快照」**逐字节恒等**，字节 diff 有唯一目标。

# ---------------------------------------------------------------- ④ 节选头（INC-009 第三次发作的机械修复）
# 病因：消息标「逐字引自」而实际是**手裁节选**（r12：105 行里贴了 17 行）。
#   与前两次同种：标签承诺的保证等级高于交付。**修法与粘贴块同型：让节选头也由本生成器产出**，
#   自带 sha256 与省略清单 ⇒ 消息侧不再有手写空间，读者可 `cmp` 校验。
import hashlib
_sha = hashlib.sha256(io.open(OUT, "rb").read()).hexdigest()[:16]
EXC = [
    f"**节选自 `{_SELF}`**（sha256 前 16 位 `{_sha}`；全文同目录 `status_paste_{STAMP}_r{ROUND}.txt`，可 `cmp` 校验）",
    "",
    f"**已省略**：D 节枚举表、事故计数与编号列表、D 基线声明、A 节 {len(ROWS)} 行明细表、快照清单行。",
    "",
    "**本份关键数**（均可在全文中逐位对上）：",
    f"- A 节 {len(ROWS)} 项：✅ {N['✅']} ｜ ✅* {N['✅*']} ｜ ⬜ {N['⬜']}",
    f"- 有效完成 {N['✅']} / 未完成 {N['✅*'] + N['⬜']}",
    f"- diff：新增 {len(added)} ｜ 删除 {len(removed)} ｜ 状态变更 {len(changed)}"
    + (f"（{', '.join(f'{k} {prev[k]}→{CUR[k]}' for k in changed)}）" if changed else ""),
    f"- D 节 ID：1–{len(D_IDS)} 共 {len(D_IDS)} 个 ｜ 事故日志 {len(INCIDENTS)} 条",
]
io.open(os.path.join(ROOT, f"status_excerpt_{STAMP}_r{ROUND}.md"), "w", encoding="utf-8").write(
    "\n".join(EXC) + "\n")

print(f"写出 {OUT}")
print(f"  粘贴块 status_paste_{STAMP}_r{ROUND}.txt ｜ 节选头 status_excerpt_{STAMP}_r{ROUND}.md"
      f" ｜ 清单归档 _checklist_archive/WORK_CHECKLIST_{STAMP}_r{ROUND}.md")
print(f"  第 {ROUND} 份 ｜ 上一份：{os.path.basename(prev_files[-1]) if prev_files else '（无）'}")
print(f"  diff：新增 {len(added)} ｜ 删除 {len(removed)} ｜ 变更 {len(changed)}")
print(f"  合计 {len(ROWS)} 项：✅ {N['✅']} ｜ ✅* {N['✅*']} ｜ ⬜ {N['⬜']}")
