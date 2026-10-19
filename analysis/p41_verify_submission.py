# -*- coding: utf-8 -*-
r"""p41 — 投稿包验收（独立于 p40 的回读校验）

================================ 查什么 ================================
`p40` 生成时已做一次回读；本脚本是**第二道**，查的是 p40 **自己查不到**的东西：

  V1 **内容零丢失**：源 Markdown 里的每个**非空正文行**，其文字必须出现在 docx 里。
     这不是抽查 —— 是逐行核对（归一化后做子串包含）。
     首版风险：md_to_doc 万一漏掉某类行（如无空行分隔的列表、缩进段落），
     生成时会静默通过，只有逐行核对才发现。
  V2 **表格数一致**：源里的 Markdown 表个数 == docx 里的表个数（按段分别核）。
  V3 **无内部产物泄漏**：批注标记、References 子表、审计词、Markdown 语法一律不得出现。
  V4 **数字守恒**：源文本里的数字 token **多重集**必须与 docx 里的完全一致
     （防止转换过程吞掉或改写任何数字）。
  V5 **包内实物齐全**：图件、验证证据、README 全在。
  V6 **Word 可打开**：docx 能被解析，且非空。
"""
import io, os, re, sys
from collections import Counter
from docx import Document

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MR   = os.path.join(ROOT, "manuscript_renumbered")
PKG  = os.path.join(ROOT, "SUBMISSION_PACKAGE_20261007")
R = []
def rec(cid, desc, ok, note=""):
    R.append(ok)
    print("  [%-4s] %-58s %s" % ("PASS" if ok else "FAIL", desc[:58], note), flush=True)

def read(p):
    return io.open(p, encoding="utf-8").read()

def norm(s):
    """归一化：统一破折号/引号/空白，便于跨格式比对。

    ★ 还要**吃掉 `X_yz` 里的下划线**：p40 现在把数学符号排成真下标，
       文字层面 `k_ij` 变成了 `kij`。两侧同规则归一后才能比对；
       代码字段名（`auroc_target`，等宽体保留了下划线）同样被吃掉，两侧一致，不影响判定。"""
    s = re.sub(r"(?<=[A-Za-z])_(?=[A-Za-z0-9])", "", s)
    s = (s.replace("\u2013", "-").replace("\u2014", "-").replace("\u2212", "-")
          .replace("\u2019", "'").replace("\u2018", "'")
          .replace("\u201c", '"').replace("\u201d", '"')
          .replace("\u00a0", " "))
    return " ".join(s.split())

def docx_all(path):
    d = Document(path)
    parts = [p.text for p in d.paragraphs]
    ntab = len(d.tables)
    for t in d.tables:
        for row in t.rows:
            for c in row.cells:
                parts.append(c.text)
    return "\n".join(parts), ntab, d

def md_tables(lines):
    n, i = 0, 0
    while i < len(lines):
        if lines[i].strip().startswith("|") and lines[i].strip().endswith("|"):
            n += 1
            while i < len(lines) and lines[i].strip().startswith("|"):
                i += 1
        else:
            i += 1
    return n

# ★ 只剥 p40 **真正会剥掉**的标记：加粗/斜体/代码。
#   **不剥 `|`** —— 段落里的竖线是正文内容（如 `|ΔAUROC| ≤ 0.02`、`n_s = |S|`），
#   p40 会原样保留。首版把 `|` 一并剥掉，导致 4 个段落被误判为「缺失」。
MD_MARK = re.compile(r"\*\*|\*|`")
def blocks(md):
    """按**空行**切段，去掉标题/表格/批注行，再做 markdown 去标记 + 空白归一。

    ★ 必须按**段落块**比对，不能按物理行 —— 稿件正文是折行写的，
       一个段落在源里占好几行，单行是它的片段，永远不可能整段出现在 docx 里。
       首版按物理行比对，产生大批假阳性，也因此**掩盖了真问题**
       （跨行的 `**…**` 因两端配不上对而原样落进了 Word）。"""
    out, buf = [], []
    def flush():
        if not buf: return
        t = " ".join(buf)
        t = MD_MARK.sub("", t)
        t = norm(t)
        if len(t) > 15:
            out.append(t)
        buf.clear()
    for ln in md.splitlines():
        s = ln.strip()
        if not s or s == "---":
            flush(); continue
        if s.startswith("#") or s.startswith(">") or s.startswith("|"):
            flush(); continue
        if re.match(r"^(\d+\.|[-*])\s", s):
            flush(); buf.append(s); continue
        buf.append(s)
    flush()
    return out

print("=" * 100); print("p41 — 投稿包验收"); print("=" * 100); print("")

# ================================================================ 正文 docx
DOC1 = os.path.join(PKG, "01_manuscript_text", "Manuscript_main_text.docx")
rec("V5-0", "正文 docx 存在", os.path.exists(DOC1), DOC1)
T1, NT1, D1 = docx_all(DOC1)
N1 = norm(T1)

SRC1 = [
    ("Abstract",      read(os.path.join(ROOT, "Abstract_v1.md"))),
    ("Introduction",  read(os.path.join(MR, "Introduction_v1.md"))),
    ("Results",       read(os.path.join(MR, "Results_v3_clean.md"))),
    ("Methods",       read(os.path.join(MR, "Methods_v1.md"))),
    ("Discussion",    read(os.path.join(MR, "Discussion_v1.md"))),
    ("References",    read(os.path.join(MR, "References_FINAL.md"))),
    # 2026-10-07 新增：Data / Code availability 两节（置于 Methods 与 References 之间）
    ("Availability",  read(os.path.join(MR, "DataCodeAvailability_v1.md"))),
]

# --- V1 内容零丢失（逐段核对）
print("【V1】内容零丢失（逐段核对；源按空行切段、去 markdown、空白归一）")
tot_lines = 0
for name, md in SRC1:
    # 只取该文件在正文里的**那一段**（与 p40 的边界一致）
    if name == "Abstract":
        seg = md.split("## Abstract", 1)[1].split("\n---", 1)[0]
    elif name == "Introduction":
        seg = md.split("## Introduction", 1)[1].split("## References cited", 1)[0]
    elif name == "Results":
        seg = md.split("## A heterogeneous", 1)[1]
    elif name == "Methods":
        seg = md.split("# Methods", 1)[1].split("## References cited", 1)[0]
    elif name == "Discussion":
        seg = md.split("## Discussion", 1)[1].split("## References used", 1)[0]
    elif name == "Availability":
        seg = md.split("## Data availability", 1)[1]
    else:
        # References：边界与 p40 一致 —— **从表头行起**，表头之前的说明文字属内部产物
        seg = "| # |" + md.split("| # |", 1)[1]
    bl = blocks(seg)
    miss = [b for b in bl if b not in N1]
    tot_lines += len(bl)
    rec("V1-" + name[:4], f"{name}: {len(bl)} 个段落全部出现在 docx 中",
        not miss, f"缺 {len(miss)} 段" + (f"：{miss[0][:70]!r}" if miss else ""))
print(f"     逐段核对合计 {tot_lines} 段")
print("")

# --- V2 表格数一致
print("【V2】表格数一致")
src_tab = 0
for name, md in SRC1:
    if name == "Methods":
        seg = md.split("# Methods", 1)[1].split("## References cited", 1)[0]
        src_tab += md_tables(seg.splitlines())
    # ★ References **不再计入源侧表数**：p40 现在把它渲染为**编号列表**（期刊格式），
    #   源文件的三列表是内部对账格式，本就不该出现在正文里。
    #   首版仍按表格计数，于是「源 2 表 vs docx 1 表」被误判为不一致。
# Results 也可能有表
res_seg = SRC1[2][1].split("## A heterogeneous", 1)[1]
src_tab += md_tables(res_seg.splitlines())
rec("V2-1", f"源 Markdown 表数 == docx 表数（源 {src_tab}）", NT1 == src_tab,
    f"docx {NT1} 表 vs 源 {src_tab} 表")
print("")

# --- V3 无内部产物泄漏
print("【V3】无内部产物泄漏")
FORBID = {
    "Markdown 加粗": r"\*\*", "Markdown 表格分隔": r"\|\s*-{3,}",
    "Markdown 引用块": r"(?m)^\s*>\s",
    # 标题：要求 `#` 后有**实际内容**（References 表的列头就是一个孤立的 "#"，不是标题）
    "Markdown 标题": r"(?m)^#{1,6}\s+\S",
    "References 子表": r"References cited in this|References used in this",
    # 内部批注：只列**确定属于内部产物**的标题词。
    # ★ 不列 "审计" / "audit log" —— 正文里有一句合法的
    #   "a dated analysis audit log records every analysis decision"，首版把它误判为泄漏。
    "内部批注": r"内部记录更正|为什么是这一版|被否的候选|## Abstract|规范约束",
    "反引号": r"`",
}
# ★ 逐**单元**扫描，不把段落与表格单元格拼成一整块再扫 ——
#   拼接会在边界处造出假模式（如 `... #` + `\n` + `F...` 被当成行首标题）。
#   首版即因此误报一次。
CELLS1 = [c.text for t in D1.tables for row in t.rows for c in row.cells]
UNITS1 = [p.text for p in D1.paragraphs] + CELLS1
for label, pat in FORBID.items():
    hits = [u for u in UNITS1 if re.search(pat, u)]
    rec("V3-" + label[:6], f"正文 docx 无「{label}」", not hits,
        f"命中 {len(hits)}" + (f"：{hits[0][:50]!r}" if hits else ""))
print("")

# --- V4 数字守恒
print("【V4】数字守恒（源 ↔ docx 的数字多重集）")
NUM = re.compile(r"(?<![\w.])\d[\d,]*(?:\.\d+)?")
def all_nums(md):
    """**口径上无可争议的**数字抽取：源段里的**全部**数字，只去掉行内标记符号
       （`**` `*` 反引号），**不做任何行过滤**。

       ★ 为什么不再按行过滤：首版按「会变成段落的行」过滤，必须逐条镜像 p40 的渲染规则
         （`#`→标题段、`>`→段落、`|`→表格、列表行→独立段……）。镜像永远差一条，
         实测就为了一个 `0.5` 反复对不上 —— 而每次都可能是**真丢失**被掩盖成口径问题。
         改成两边都不过滤后，比对的是「源里有的数字 docx 里都有，反之亦然」，
         标题与表格也一并纳入，比原先覆盖更全。
       docx 侧对应取 **段落 + 表格单元格**。"""
    return MD_MARK.sub("", md)

_WORD = re.compile(r"[A-Za-zÀ-ɏ][A-Za-zÀ-ɏ'\-]{2,}")
src_nums = Counter()
src_words = set()                      # ★ INC-036：词级守恒，见下方 V4-3
for name, md in SRC1:
    if name == "Abstract":
        seg = md.split("## Abstract", 1)[1].split("\n---", 1)[0]
    elif name == "Introduction":
        seg = md.split("## Introduction", 1)[1].split("## References cited", 1)[0]
    elif name == "Results":
        seg = md.split("## A heterogeneous", 1)[1]
    elif name == "Methods":
        seg = md.split("# Methods", 1)[1].split("## References cited", 1)[0]
    elif name == "Discussion":
        seg = md.split("## Discussion", 1)[1].split("## References used", 1)[0]
    elif name == "Availability":
        # 与 p40 同边界：从 "## Data availability" 起，到「给作者的填写说明」止
        seg = md.split("## Data availability", 1)[1].split("> **给作者的填写说明", 1)[0]
    else:
        # References：与 p40 边界一致 —— 从表头行起（表前的内部说明文字不进 Word）
        seg = "| # |" + md.split("| # |", 1)[1]
    src_nums.update(NUM.findall(all_nums(seg)))
    # ★ V4-3 只统计**散文各节**。References 那一节 p40 是按设计**重建**的
    #   （三列表 → 编号列表），其内部表头文字本就不该出现在 Word 里
    #   —— 若把 References 计入，会得到一条**永远为红**的判据，那种判据等于没有。
    if name != "References":
        src_words.update(w.lower() for w in _WORD.findall(all_nums(seg)))
# docx 侧：同口径 —— 段落（含标题段）+ 全部表格单元格，都不过滤。
doc_nums = Counter(NUM.findall(T1))
only_src = src_nums - doc_nums
extra = doc_nums - src_nums
rec("V4-1", "docx 未丢失任何数字 token", not only_src,
    f"仅源有 {sum(only_src.values())} 个" + (f"：{list(only_src)[:6]}" if only_src else ""))
rec("V4-2", "docx 未凭空多出数字 token", not extra,
    f"仅 docx 有 {sum(extra.values())} 个" + (f"：{list(extra)[:6]}" if extra else ""))

# ★★ 2026-10-12（INC-036）：**只看数字是不够的**。
#   病历：Methods §1 的一行 `#STU-2023-0583, and the MGB data are distributed de-identified
#   by PhysioNet under credentialed` 因行首是 `#` 被当成 1 级标题、**整行被丢弃**。
#   它能被抓到，**纯粹因为那行含数字**。⇒ 一行不含数字的正文若同样被丢，
#   V4-1/V4-2 全绿，交付的 Word **悄悄少一段而没有任何判据会响**。
#   ⇒ 补一条**词级**守恒：源里的每个词形（小写、去标点）必须至少在 docx 里出现一次。
#   用**词形集合**而非计数：docx 侧会因折行合并、连字符重排而改变词序/切分，
#   集合口径对排版差异稳健，对**整行消失**敏感 —— 正是本条要抓的那一类。
# 判定用**子串**而非**整词**：docx 侧会把 `AUROC_aux` 这类下标渲染成多个 run，
# `p.text` 拼回来可能变成 `AUROCaux` —— 整词口径会把它误报成「丢失」。
# 子串口径对**渲染造成的切分差异**稳健，对**整行消失**同样敏感（行没了，词也没了）。
_doc_flat = " ".join(T1.lower().split())
_missing_words = sorted(w for w in src_words if w not in _doc_flat)
rec("V4-3", "docx 未丢失任何**词形**（不依赖该行是否含数字；用子串口径）", not _missing_words,
    "仅源有 %d 个词形%s" % (len(_missing_words),
                            ("：" + ", ".join(_missing_words[:12])) if _missing_words else "")
    + "｜未纳入本节者：References（p40 有意把内部三列表**重建**为编号列表，"
      "表头 `Full citation` / `DOI` 两列标题按设计不进入 Word）")
print("")

# ================================================================ 补充 docx
print("【V1-V3】补充表 docx")
DOC2 = os.path.join(PKG, "02_supplementary_tables", "Supplementary_Tables.docx")
rec("V5-0b", "补充 docx 存在", os.path.exists(DOC2), DOC2)
T2, NT2, D2 = docx_all(DOC2)
UNITS2 = ([p.text for p in D2.paragraphs]
          + [c.text for t in D2.tables for row in t.rows for c in row.cells])
for label, pat in {"Markdown 表格分隔": r"\|\s*-{3,}",
                   "Markdown 标题": r"^#{1,6}\s+\S",   # 逐单元扫，「#」列头不会误伤
                   "反引号": r"`", "Markdown 加粗": r"\*\*"}.items():
    hits = [u for u in UNITS2 if re.search(pat, u)]
    rec("V3b-" + label[:6], f"补充 docx 无「{label}」", not hits,
        f"命中 {len(hits)}" + (f"：{hits[0][:50]!r}" if hits else ""))
RP = os.path.join(ROOT, "figures", "data", "reporting_profiles")
for f, sec in [("Supplementary_Table_S1.md", "Supplementary Table S1"),
               ("Supplementary_Table_S2.md", "Supplementary Table S2"),
               ("Supplementary_Tables_R1_R3.md", "Table R1"),
               ("Supplementary_Methods.md", "Supplementary Methods")]:
    src = read(os.path.join(RP, f))
    n_src = md_tables(src.splitlines())
    rec("V2b-" + f[:22], f"{f}: 表数一致（源 {n_src}）", n_src <= NT2, f"docx 共 {NT2} 表")
    rec("V1b-" + f[:22], f"{f}: 章节标题在 docx 中", sec in T2)
print("")

# ================================================================ V5 实物齐全
print("【V5】包内实物齐全")
CHECKS = [
    ("README_投稿包说明.md", os.path.join(PKG, "README_投稿包说明.md")),
    ("_build_log.txt",       os.path.join(PKG, "_build_log.txt")),
]
for i in range(1, 6):
    for e in ("png", "svg", "tiff"):
        CHECKS.append((f"Figure{i}.{e}", os.path.join(PKG, "03_figures_main", f"Figure{i}.{e}")))
CHECKS.append(("Figure_Captions.md", os.path.join(PKG, "03_figures_main", "Figure_Captions.md")))
for i in range(1, 10):
    for e in ("png", "svg", "tiff"):
        CHECKS.append((f"FigureS{i}.{e}", os.path.join(PKG, "04_figures_supplementary", f"FigureS{i}.{e}")))
CHECKS.append(("Supplementary_Legends.md",
               os.path.join(PKG, "04_figures_supplementary", "Supplementary_Legends.md")))
for n in ("report_evs.txt", "coverage_table.csv", "manifest_sha256.txt", "SHA256SUMS.txt"):
    CHECKS.append((n, os.path.join(PKG, "05_verification", n)))
missing = [n for n, p in CHECKS if not os.path.exists(p)]
rec("V5-1", f"包内 {len(CHECKS)} 项实物全部存在", not missing, f"缺 {missing or '无'}")
print("")

# ================================================================ V6 可打开
print("【V6】Word 可打开且非空")
for nm, p in [("正文", DOC1), ("补充表", DOC2)]:
    d = Document(p)
    nch = sum(len(q.text) for q in d.paragraphs) + sum(
        len(c.text) for t in d.tables for row in t.rows for c in row.cells)
    rec("V6-" + nm, f"{nm} docx 可解析且非空（{nch} 字符）", nch > 5000)
print("")

# ================================================================ V7 稿件形态
# ★★ 2026-10-07 补：**这一类检查此前完全不存在**，后果是 33/33 的「通过」放行了一份
#    不能投稿的文件 —— 三列参考文献表、中文表注、字面 `k_ij`，全部是使用者用眼睛发现的。
#    V1–V6 查的是「内容有没有丢、数字对不对」；**V7 查的是「这份文件像不像一篇稿件」**。
#    两类必须都有：只查前者，产物可以是「内容完整但根本不能投」。
print("【V7】稿件形态（此前缺失的一类）")
CJK = re.compile(r"[　-〿㐀-䶿一-鿿＀-￯]")

for nm, path in [("正文", DOC1), ("补充表", DOC2)]:
    d = Document(path)
    units = [p.text for p in d.paragraphs] + \
            [c.text for t in d.tables for row in t.rows for c in row.cells]
    hits = sorted({u.strip()[:70] for u in units if CJK.search(u)})
    rec("V7-1" + nm[:2], f"{nm} docx **不含任何中文/全角字符**", not hits,
        f"命中 {len(hits)} 段" + (f"：{hits[0]!r}" if hits else ""))

# 内部产物字符串
INTERNAL = {
    "脚本文件名": r"\bp\d{2}[a-z]?_[A-Za-z0-9_]*\.py",
    "渲染器版本号": r"\bp2\d[a-z] v[\d.]+",
    "内部路径": r"分析结果/|figures/data/|figures/final/|reporting_profiles/",
    "自述性注记": r"自动渲染|不手工录入|REPORTING-ONLY DERIVATION|provenance\.csv",
}
for nm, path in [("正文", DOC1), ("补充表", DOC2)]:
    d = Document(path)
    units = [p.text for p in d.paragraphs] + \
            [c.text for t in d.tables for row in t.rows for c in row.cells]
    for label, pat in INTERNAL.items():
        hits = sorted({u.strip()[:70] for u in units if re.search(pat, u)})
        rec("V7-2" + nm[:2] + label[:4], f"{nm} docx 不含「{label}」", not hits,
            f"命中 {len(hits)}" + (f"：{hits[0]!r}" if hits else ""))

# 参考文献必须是编号列表，不是表格
d = Document(DOC1)
# ★★ 必须**限定在 References 标题之后**再数。
#   首版直接在全文数「以 `N. ` 开头的段落」，把 **Methods §1–§15 的节标题**也算了进去 ——
#   于是即使参考文献整段丢失，这条检查照样通过（又一个假绿）。
_paras = [p for p in d.paragraphs]
_ri = next((i for i, p in enumerate(_paras)
            if p.text.strip() == "References" and p.style.name.startswith("Heading")), None)
rec("V7-3a", "正文中存在 References 标题段", _ri is not None,
    f"位置 段{_ri}" if _ri is not None else "找不到")
ref_paras = [p.text for p in (_paras[_ri + 1:] if _ri is not None else [])
             if re.match(r"^\s*\d+\.\s+\S", p.text)]
tab_like = [t for t in d.tables
            if t.rows and any(re.search(r"full citation|^#$", c.text, re.I)
                              for row in t.rows[:1] for c in row.cells)]
rec("V7-3", "参考文献是**编号列表**而非三列表格",
    len(ref_paras) >= 20 and not tab_like,
    f"编号段落 {len(ref_paras)} 条；疑似参考文献表格 {len(tab_like)} 个")
rec("V7-3b", "正文仅剩的表格是内容表（不该再有参考文献表）", len(d.tables) <= 3,
    f"正文表格数 {len(d.tables)}")

# 数学下标：正文里遗留的 `X_yz` 必须都是**代码字段名**，不能是数学符号
CODE_TOKENS = set()
for _, md in SRC1:
    CODE_TOKENS |= set(re.findall(r"`([A-Za-z][A-Za-z0-9_]*_[A-Za-z0-9_]+)`", md))
doc_ids = set(re.findall(r"\b[A-Za-z]{1,8}_[A-Za-z0-9_]{1,24}\b", T1))
stray = sorted(doc_ids - CODE_TOKENS)
rec("V7-4", "正文里没有**未按下标排版**的数学符号（`k_ij` 这类）", not stray,
    f"遗留 {len(stray)} 种" + (f"：{stray[:8]}" if stray else "") +
    f"（允许的代码字段名 {len(CODE_TOKENS)} 种）")
print("")

# ================================================================ V8 主表规范
# ★★ 2026-10-07 补：使用者发现正文里有一张**无编号、无引用、嵌在 Methods 中间**的表格。
#    Nature Communications 的四条硬规定一条都没满足：
#      · 主表必须连续编号（Table 1, 2, …）；· 必须在正文被引用；
#      · 必须**集中在正文文件末尾**；· 表格内**不得用粗体**（除非表注中定义）。
#    这一类此前同样没有任何检查 —— 与 V7 一样，只有「内容/数字」层面的判据。
print("【V8】主表规范")
_d = Document(DOC1)
_ps = [p for p in _d.paragraphs]
_idx_ref = next((i for i, p in enumerate(_ps) if p.text.strip() == "References"
                 and p.style.name.startswith("Heading")), None)
_capt = [p.text.strip() for p in _ps if re.match(r"^Table \d+ \|", p.text.strip())]
_txt = "\n".join(p.text for p in _ps)

rec("V8-1", f"正文表格数 == 编号主表数（docx {len(_d.tables)} 表 / 标题 {len(_capt)} 条）",
    len(_d.tables) == len(_capt), f"{len(_d.tables)} vs {len(_capt)}")
_nums = [int(re.match(r"^Table (\d+) \|", c).group(1)) for c in _capt]
rec("V8-2", "主表编号从 1 连续", _nums == list(range(1, len(_nums) + 1)), f"{_nums}")
rec("V8-3", "每条主表都有概括性标题（`Table N | …`）",
    all(len(c.split("|", 1)[1].strip()) > 10 for c in _capt),
    f"最短标题 {min((len(c.split('|',1)[1].strip()) for c in _capt), default=0)} 字符")
# 位置：表必须在 References 之后（集中在正文文件末尾）
_ti = next((i for i, p in enumerate(_ps) if re.match(r"^Table \d+ \|", p.text.strip())), None)
rec("V8-4", "主表集中在正文文件**末尾**（References 之后）",
    _idx_ref is not None and _ti is not None and _ti > _idx_ref,
    f"References 段{_idx_ref} → Table 1 段{_ti}")
# 引用：正文里必须出现 "Table N"
for n in _nums:
    rec(f"V8-5-{n}", f"Table {n} 在正文中被引用",
        len(re.findall(rf"\bTable {n}\b", _txt)) >= 2,
        f"全文出现 {len(re.findall(rf'Table {n}', _txt))} 次（含标题 1 次）")
# 表体不得有粗体
_bold_body = []
for t in _d.tables:
    for ri, row in enumerate(t.rows):
        if ri == 0:
            continue
        for c in row.cells:
            if any(r.bold for p in c.paragraphs for r in p.runs if r.text.strip()):
                _bold_body.append(c.text[:30])
rec("V8-6", "表体无粗体（期刊要求，除非在表注中定义）", not _bold_body,
    f"命中 {len(_bold_body)}" + (f"：{_bold_body[:3]}" if _bold_body else ""))

# ── 加粗只允许出现在三处（Nature Portfolio 明确许可）：Methods 的段落引导词、参考文献卷号、表头。
#    「Please do not use italics or bold for emphasis」—— 正文强调性粗体在生产环节会被直接删掉。
_bad_bold, _runin, _vol = [], 0, 0
_first_text_seen = False
for p in _ps:
    t = p.text.strip()
    if not t or p.style.name.startswith("Heading"):
        continue
    # 论文**标题**：粗体合规（且通常居中）
    if not _first_text_seen:
        _first_text_seen = True
        continue
    # **表题**前缀 `Table N |` 粗体：Nature Communications 的表题格式
    if re.match(r"^Table \d+ \|", t):
        continue
    if re.match(r"^\d+\.\s", t):                        # 参考文献条目 → 卷号粗体合规
        _vol += 1
        continue
    bolds = [r.text.strip() for r in p.runs if r.bold and r.text.strip()]
    if not bolds:
        continue
    # run-in heading：段首加粗、且以句号结尾
    first = bolds[0]
    if t.startswith(first) and first.endswith(".") and len(bolds) == 1:
        _runin += 1
        continue
    _bad_bold.append((t[:26], bolds[:3]))
rec("V8-7", "正文无**强调性**粗体（只允许 Methods 段首引导词与参考文献卷号）", not _bad_bold,
    f"违规 {len(_bad_bold)} 段（合规：段首引导词 {_runin}、参考文献 {_vol}）"
    + (f"：{_bad_bold[:2]}" if _bad_bold else ""))
print("")

nbad = sum(1 for x in R if not x)
print("=" * 100)
print(f"p41 验收合计 {len(R)} 项：PASS {len(R)-nbad}｜**非 PASS {nbad}**")
print("*** p41 通过 —— 投稿包逐行核对无内容丢失、无泄漏、数字守恒 ***" if not nbad
      else "*** p41 未通过 ***")
sys.exit(1 if nbad else 0)
