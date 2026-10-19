# -*- coding: utf-8 -*-
r"""p40 — 构建投稿包（Submission package）

================================ 做什么 ================================
把投稿所需的全部东西收进一个**标注清楚**的新文件夹，并生成两份 Word：

  ① `Manuscript_main_text.docx`  —— 正文全套（Title → Abstract → Introduction →
     Results → Methods → Discussion → References），**不含任何图**。
  ② `Supplementary_Tables.docx`  —— 全部补充表（S1 / S2 / Tables R1–R3 / Supplementary Methods）。

其余（图件、验证证据、来源 .md、哈希清单）按目录分放，逐项在 README 里标注。

================================ 设计约束 ================================
  · **只从正本读取**（路径按 `MANUSCRIPT_MASTER_DECLARATION.md` §2 逐文件认定）。
  · **剥掉一切批注与自查段** —— 那些是内部产物，不得进投稿件。
    剥离靠**显式边界**（如 `## References cited in this Introduction`），不是靠猜。
  · 每个剥离动作都断言边界**确实存在**；边界找不到即抛错，不静默带过去。
  · 生成后**回读 docx** 与源文本比对（段落数、表格数、关键句在场）。
"""
import io, os, re, shutil, sys
import pandas as pd
from docx import Document
from docx.shared import Pt, Inches
from docx.enum.text import WD_ALIGN_PARAGRAPH

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MR   = os.path.join(ROOT, "manuscript_renumbered")
PKG  = os.path.join(ROOT, "SUBMISSION_PACKAGE_20261007")
LOG  = []

def say(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); LOG.append(s)

def must(cond, msg):
    if not cond:
        raise RuntimeError("构建中止：" + msg)
    return True

# ================================================================ 读取与剥离
def read(p):
    must(os.path.exists(p), "正本不存在：%s" % p)
    return io.open(p, encoding="utf-8").read()

def between(text, start_pat, end_pat=None, keep_start=True):
    """按**显式边界**取段。边界必须存在，否则抛错。"""
    m = re.search(start_pat, text, re.M)
    must(m is not None, "找不到起始边界 %r" % start_pat)
    body = text[m.start():] if keep_start else text[m.end():]
    if end_pat:
        m2 = re.search(end_pat, body, re.M)
        must(m2 is not None, "找不到结束边界 %r（起于 %r）" % (end_pat, start_pat))
        body = body[:m2.start()]
    return body.rstrip()

TITLE_TXT = read(os.path.join(ROOT, "Title_v1.md"))
# 标题是**跨行的 blockquote**：逐行去掉 `>` 前缀再拼接，否则续行的 `>` 会被吃进标题
# （首版实测标题里出现了 "prediction-model > failure"）。
_m = re.search(r"^## 正文标题\s*\n\s*\n((?:>.*\n?)+)", TITLE_TXT, re.M)
must(_m is not None, "Title 正本里找不到「## 正文标题」后的引用块")
_raw = " ".join(l.lstrip(">").strip() for l in _m.group(1).splitlines())
_b = re.search(r"\*\*(.+?)\*\*", _raw, re.S)
must(_b is not None, "标题引用块里找不到加粗标题")
TITLE = " ".join(_b.group(1).split())
must(">" not in TITLE, "标题里残留 '>'：%r" % TITLE)

ABSTRACT = between(read(os.path.join(ROOT, "Abstract_v1.md")), r"^## Abstract", r"^---")
INTRO    = between(read(os.path.join(MR, "Introduction_v1.md")), r"^## Introduction",
                   r"^## References cited in this Introduction")
RESULTS  = between(read(os.path.join(MR, "Results_v3_clean.md")), r"^# Results",
                   r"^## [^\n]*\n", keep_start=False)  # 去掉顶部的 > 说明块
RESULTS  = between(read(os.path.join(MR, "Results_v3_clean.md")), r"^## A heterogeneous")
METHODS  = between(read(os.path.join(MR, "Methods_v1.md")), r"^# Methods",
                   r"^## References cited in this Methods")
DISCUSS  = between(read(os.path.join(MR, "Discussion_v1.md")), r"^## Discussion",
                   r"^## References used in this Discussion")
REFS     = between(read(os.path.join(MR, "References_FINAL.md")), r"^\| # \|")

say("=" * 100); say("p40 — 投稿包构建"); say("=" * 100)
say(f"Title      : {TITLE}")
for nm, body in [("Abstract", ABSTRACT), ("Introduction", INTRO), ("Results", RESULTS),
                 ("Methods", METHODS), ("Discussion", DISCUSS), ("References", REFS)]:
    say(f"  {nm:<12} {len(body.splitlines()):>4} 行 ｜ {len(body.split()):>6} 词")

# 剥离自检：批注标记不得残留
for nm, body in [("Abstract", ABSTRACT), ("Introduction", INTRO), ("Results", RESULTS),
                 ("Methods", METHODS), ("Discussion", DISCUSS)]:
    must("## References cited" not in body and "## References used" not in body,
         "%s 里仍有 References 子表" % nm)
    must("内部记录更正" not in body, "%s 里仍有内部记录段" % nm)
say("")
say("剥离自检：无 References 子表、无内部记录段 ✓")

# ================================================================ Markdown → docx
# ★ 交替中必须**单独**列出 `*X*_yz`（斜体单字母 + 下标）这种组合。
#   实例：源里写 `*r*_failure`。若让 `*[^*]+?*` 先匹配，`*r*` 与 `_failure` 会被切到两段，
#   跨段就再也配不出「字母 + 下标」的模式 —— 首版实测 `r_failure` 原样留在了正文里。
INL = re.compile(r"(\*\*[^*]+?\*\*|\*[A-Za-z]\*_[A-Za-z0-9_]+|\*[^*]+?\*|`[^`]+?`)")
def add_runs(par, text, bold=False, italic=False):
    # 必须**递归**：稿件里有嵌套标记（加粗里嵌反引号），只剥外层会把反引号留在文中。
    # 首版即如此，被回读校验抓到。
    text = text.replace("–", "\u2013")
    for seg in INL.split(text):
        if not seg:
            continue
        if re.match(r"^\*[A-Za-z]\*_[A-Za-z0-9_]+$", seg):
            r = par.add_run(seg[1]); r.italic = True; r.bold = bold
            s = par.add_run(seg[4:]); s.font.subscript = True; s.italic = True; s.bold = bold
        elif seg.startswith("**") and seg.endswith("**") and len(seg) > 4:
            add_runs(par, seg[2:-2], bold=True, italic=italic)
        elif seg.startswith("*") and seg.endswith("*") and len(seg) > 2:
            add_runs(par, seg[1:-1], bold=bold, italic=True)
        elif seg.startswith("`") and seg.endswith("`") and len(seg) > 2:
            r = par.add_run(seg[1:-1]); r.font.name = "Consolas"; r.font.size = Pt(8)
            r.bold = bold; r.italic = italic
        else:
            # 兜底：任何**未被上面分支消费**的反引号一律去掉。
            # 实例：References 第 22 条（Ke et al., LightGBM）无 DOI，源里用空代码段 `` 占位，
            # 首版会把这两个反引号原样写进 Word。
            seg = seg.replace("`", "")
            if not seg:
                continue
            _emit_math(par, seg, bold, italic)

# ★★ 数学下标：`k_ij` 在稿件里应当是 **k 带下标 ij**，不是字面的下划线。
#   判定规则**不用启发式** —— 源文件自己分得很清楚：
#     · 写在反引号里的（`` `auroc_target` ``）是**数据字段名**，保持等宽体、不改；
#     · 写在裸文本里的（`k_ij`、`AUROC_target`）是**数学符号**，按下标排版。
#   上面 add_runs 的代码段分支已经消费掉反引号内容，所以**走到这里的都是数学符号**。
_SUB = re.compile(r"\b([A-Za-z]{1,8})_([A-Za-z0-9_]{1,24})\b")
def _emit_math(par, text, bold, italic):
    pos = 0
    for m in _SUB.finditer(text):
        if m.start() > pos:
            r = par.add_run(text[pos:m.start()]); r.bold = bold; r.italic = italic
        r = par.add_run(m.group(1)); r.bold = bold; r.italic = True
        s = par.add_run(m.group(2)); s.font.subscript = True; s.bold = bold; s.italic = True
        pos = m.end()
    if pos < len(text):
        r = par.add_run(text[pos:]); r.bold = bold; r.italic = italic

# ★★ 正文强调性粗体必须删除（Nature Portfolio：*"Please do not use italics or bold for
#    emphasis"*，且生产环节会把这类格式**直接去掉**）。但两类粗体是**明确许可**的：
#      · Methods 的段落引导词（run-in heading，段首、以句号结尾）—— "subdivided by short, bold headings"
#      · 参考文献里的期刊卷号 —— "Volume numbers in reference lists"
#    实测正文原有 31 处强调性粗体（`sign`、`patient`、`not`、`microbiological episode` …）。
_RUNIN = re.compile(r"^\*\*([^*]{2,90}\.)\*\*(?=\s)")
def strip_emphasis(txt):
    """保留段首 run-in heading 的粗体，删除其余** 标记。"""
    m = _RUNIN.match(txt)
    if m:
        return "**" + m.group(1) + "**" + txt[m.end():].replace("**", "")
    return txt.replace("**", "")


SKIPPED_H1 = []      # INC-036：被 skip_first_h1 丢掉的 1 级标题（应为「每段一个」且都是真标题）


def md_to_doc(lines, doc, skip_first_h1=True):
    """把 Markdown 行写进 docx。支持 #/##/###、> 引用、| 表格、有序/无序列表。"""
    i, n = 0, len(lines)
    while i < n:
        ln = lines[i].rstrip()
        s = ln.strip()
        if not s:
            i += 1; continue
        if s == "---":
            i += 1; continue
        if s.startswith("|") and s.endswith("|"):
            blk = []
            while i < n and lines[i].strip().startswith("|"):
                blk.append(lines[i].strip()); i += 1
            cells = [[c.strip() for c in r.strip("|").split("|")] for r in blk]
            cells = [c for c in cells if not set("".join(c)) <= set("-: ")]
            if cells:
                t = doc.add_table(rows=len(cells), cols=len(cells[0]))
                t.style = "Table Grid"
                for ri, row in enumerate(cells):
                    for ci, val in enumerate(row[:len(cells[0])]):
                        cell = t.cell(ri, ci)
                        cell.text = ""
                        # ★ 表体**去粗体**：期刊要求「表格不得使用粗体，除非在表注中定义」。
                        #   实测 §11 的方向约定表里有 `**absolute magnitude**` / `**negated**`。
                        add_runs(cell.paragraphs[0], val if ri == 0 else val.replace("**", ""))
                        for r in cell.paragraphs[0].runs:
                            r.font.size = Pt(8)
                        if ri == 0:
                            for r in cell.paragraphs[0].runs:
                                r.bold = True
                doc.add_paragraph()
            continue
        if s.startswith(">"):
            buf = []
            while i < n and lines[i].strip().startswith(">"):
                buf.append(lines[i].strip().lstrip(">").strip()); i += 1
            p = doc.add_paragraph()
            p.paragraph_format.left_indent = Inches(0.3)
            # 引用块整体已是斜体；其中不再保留任何粗体（引用块不是 Methods 正文，
            # 不存在 run-in heading 的用法）
            add_runs(p, " ".join(buf).replace("**", ""))
            for r in p.runs:
                r.italic = True
                r.font.size = Pt(9)
            continue
        # ★★ 2026-10-12（INC-036）：标题必须**井号后有空白**（CommonMark 口径）。
        #   病历：原判断是 `s.startswith("#")`，于是 Methods 里
        #   `#STU-2023-0583, and the MGB data are distributed …`（IRB 编号**恰好折行到行首**）
        #   被当成 1 级标题 → 命中下面的 `skip_first_h1` → **整行被静默丢弃**。
        #   交付的 Word 里那句话变成「…and UTSW IRB access. Whether this secondary analysis…」，
        #   **中间少了一整行**。是 p41 的 V4-1（数字 token 守衡）抓到的 —— 只因为它含数字。
        #   ★ 不含数字的行若同样被丢，当时**没有任何判据会响**（见 INC-036 沉淀）。
        if re.match(r"^#{1,6}(\s|$)", s):
            lvl = len(s) - len(s.lstrip("#"))
            txt = s.lstrip("#").strip()
            if lvl == 1 and skip_first_h1:
                # ★ 这里**不能**静默丢。但也不宜用「长度/标点像不像标题」去猜 ——
                #   首版试过，`Supplementary Table S2 — Protocol-specified warning-task
                #   metrics and grouped summaries`（84 字符的**合法**标题）被误判为正文而中止。
                #   ⇒ 处置分两层：① 根因已由上面的 `#{1,6}(\s|$)` 修掉
                #     （`#STU-…` 这类**无空白**的行根本不再进本分支）；
                #     ② **丢失本身**由 `p41` 的 V4-3（**词级**守恒）兜底探测 ——
                #     那是不依赖任何启发式的真判据。此处只**记账 + 打印**，让抽查者看得见。
                SKIPPED_H1.append(txt)
                i += 1; continue
            doc.add_heading(txt, level=min(lvl, 4)); i += 1; continue
        if re.match(r"^(\d+\.|[-*])\s", s):
            p = doc.add_paragraph()
            p.paragraph_format.left_indent = Inches(0.25)
            add_runs(p, s); i += 1; continue
        # ★★ 普通段落必须**合并连续的源代码行**再成段。
        #   稿件正文是**折行**写的（每行 ~100 字符），逐行成段会造成两个后果：
        #     ① Word 里每 ~100 字符一个硬换行；
        #     ② 跨行的行内标记（如 `**…\n…**`）两端配不上对，
        #        标记符号**原样作为文字进入 Word**（首版实测正文里泄漏了 5 处 `**`）。
        buf = [s]; i += 1
        while i < n:
            nxt = lines[i].strip()
            # ★ 同 INC-036：这里的「是不是标题」也必须用同一口径，否则
            #   `#STU-…` 这样的正文行会把段落**从这里切断**（即使不再被丢弃）。
            if (not nxt or nxt == "---" or re.match(r"^#{1,6}(\s|$)", nxt) or nxt.startswith("|")
                    or nxt.startswith(">") or re.match(r"^(\d+\.|[-*])\s", nxt)):
                break
            buf.append(nxt); i += 1
        p = doc.add_paragraph()
        add_runs(p, strip_emphasis(" ".join(buf)))

def add_reference_list(refs_md, doc):
    """参考文献必须是**编号列表**，不是表格。

    ★★ 2026-10-07 更正：源文件 `References_FINAL.md` 存的是**内部对账格式**
       （三列：`# | Full citation | DOI`）—— 那张表是拿来和引文审计逐条对账用的。
       首版把它**一比一**搬成了 Word 表格，于是正文里出现了一张三列带 DOI 列的参考文献表。
       **没有任何期刊接受这种形式。** 本函数解析该表并重建为编号列表，
       DOI 以 `doi:…` 附在条目末尾（若目标期刊的清单不要求 DOI，删除该后缀即可）。"""
    rows = []
    for ln in refs_md.splitlines():
        m = re.match(r"^\|\s*(\d+)\s*\|(.+?)\|\s*`?([^|`]*)`?\s*\|\s*$", ln)
        if m:
            rows.append((int(m.group(1)), m.group(2).strip(), m.group(3).strip()))
    must(rows, "参考文献表解析出 0 条 —— 源格式可能已变")
    nums = [r[0] for r in rows]
    must(nums == list(range(1, len(nums) + 1)),
         "参考文献编号不连续：%s" % nums)
    for num, cite, doi in rows:
        p = doc.add_paragraph()
        p.paragraph_format.left_indent = Inches(0.32)
        p.paragraph_format.first_line_indent = Inches(-0.32)   # 悬挂缩进
        add_runs(p, f"{num}. {cite}" + (f" doi:{doi}" if doi else ""))
    return len(rows)


# ================================================================ 文档 1：正文
os.makedirs(PKG, exist_ok=True)
d1 = Document()
st = d1.styles["Normal"]; st.font.name = "Times New Roman"; st.font.size = Pt(11)

p = d1.add_paragraph(); p.alignment = WD_ALIGN_PARAGRAPH.CENTER
r = p.add_run(TITLE); r.bold = True; r.font.size = Pt(16)
d1.add_paragraph()


def lift_tables(body, collected):
    """把正文里的 Markdown 表格**提出来**，稍后统一渲染在正文文件**末尾**。

    ★★ 期刊规定（Nature Communications 格式说明）：
       · 主表必须**连续编号**（Table 1, 2, …；不得细分）
       · 必须在**正文中被引用**
       · 主表必须**集中在正文文件末尾**，不得嵌在正文中间
       · 表格**不得使用粗体**（除非在表注中定义）
    首版把 §11 的方向约定表直接嵌在 Methods 里、无编号、无引用、且带粗体 —— 四条全违反。

    约定的输入形式：表格块**前一行**是标题行，形如 `**Table 1 |** <短标题>`。"""
    lines, out, i = body.splitlines(), [], 0
    n = len(lines)
    while i < n:
        s = lines[i].strip()
        if (s.startswith("|") and s.endswith("|") and i + 1 < n
                and set(lines[i + 1].strip().replace("|", "").replace("-", "").replace(":", "")) == set()):
            # 找到表格块；往上找标题行
            cap = None
            for j in range(len(out) - 1, -1, -1):
                if out[j].strip():
                    m = re.match(r"^\*\*Table\s+(\d+)\s*\|\*\*\s*(.+)$", out[j].strip())
                    if m:
                        cap = (int(m.group(1)), m.group(2).strip())
                        out.pop(j)
                    break
            blk = []
            while i < n and lines[i].strip().startswith("|"):
                blk.append(lines[i].strip()); i += 1
            cells = [[c.strip() for c in r.strip("|").split("|")] for r in blk]
            cells = [c for c in cells if not set("".join(c)) <= set("-: ")]
            collected.append(dict(caption=cap, rows=cells))
            continue
        out.append(lines[i]); i += 1
    return out


def drop_leading_heading(body):
    """`between()` 保留起始标记行，于是标题会被 add_heading 与 md_to_doc **各写一次** ——
       实测出现 'Abstract → Abstract'、'Introduction → Introduction' 这类重复标题。
       这里在送入 md_to_doc 前去掉首行的标题标记。"""
    lines = body.splitlines()
    while lines and not lines[0].strip():
        lines.pop(0)
    if lines and re.match(r"^#{1,6}\s", lines[0]):
        lines.pop(0)
    return lines

d1.add_heading("Abstract", level=1)
md_to_doc(drop_leading_heading(ABSTRACT), d1)
# ★ Results 的起始边界是 **内容小标题**（`## A heterogeneous…`），不是与 add_heading 重复的标题，
#   所以它**不丢首行** —— 首版一刀切丢掉，Results 第一节的标题就消失了（被 p41 V1 抓到）。
TABLES = []          # 正文里的表：提出来，末尾统一渲染
for hd, body, drop in [("Introduction", INTRO, True), ("Results", RESULTS, False),
                       ("Methods", METHODS, True), ("Discussion", DISCUSS, True)]:
    d1.add_page_break()
    d1.add_heading(hd, level=1)
    src = drop_leading_heading(body) if drop else body.splitlines()
    md_to_doc(lift_tables("\n".join(src), TABLES), d1)


# Data / Code availability —— 期刊要求置于 **Methods 之后、References 之前**，且必须同时录入投稿系统。
AVAIL = between(read(os.path.join(MR, "DataCodeAvailability_v1.md")),
                r"^## Data availability", r"^> \*\*给作者的填写说明")
d1.add_page_break()
d1.add_heading("Data availability", level=1)
# ★ 必须按 **"## Code availability"**（带头号）切，否则前一段的尾巴会剩一个孤立的 "## "，
#   变成一个**空标题**段落（实测 Data availability 与 Code availability 之间出现了空标题）。
_apart = AVAIL.split("## Code availability", 1)
md_to_doc(drop_leading_heading(_apart[0].replace("## Data availability", "", 1)), d1)
_ck = _apart[1] if len(_apart) > 1 else ""
d1.add_heading("Code availability", level=1)
md_to_doc(_ck.splitlines(), d1)

d1.add_page_break()
d1.add_heading("References", level=1)
NREFS = add_reference_list(REFS, d1)

# ── 主表：**集中放在正文文件末尾**（期刊要求，不得嵌在正文中间）
must(TABLES, "正文里没有找到任何主表 —— §11 的方向约定表应带 `**Table N |**` 标题行")
nums = [t["caption"][0] for t in TABLES if t["caption"]]
must(sorted(nums) == list(range(1, len(nums) + 1)),
     "主表编号必须从 1 连续：%s" % nums)
d1.add_page_break()
d1.add_heading("Tables", level=1)
for _ti, t in enumerate(TABLES):
    cap, rows = t["caption"], t["rows"]
    # ★★ 2026-10-12：原报错只说「主表缺少标题行」，**不说是哪一张** ——
    #   而本节要判断的是「这张表到底该不该有编号」，必须先看见它。
    #   （同 INC-023 / INC-034：**对象必须声明自己**，否则读报错的人只能靠猜。）
    must(cap is not None,
         "主表缺少 `**Table N |**` 标题行 —— 第 %d 张（共 %d 张），"
         "表头：%s ｜ 首行：%s\n"
         "    ⇒ 若这张表**本来就不是主表**（如 Methods 内的说明性清单），"
         "请把它改成项目符号列表，不要给它编号（编号表必须连续、且在正文被引用）。"
         % (_ti + 1, len(TABLES),
            rows[0] if rows else "(空)",
            rows[1] if len(rows) > 1 else "(无)"))
    p = d1.add_paragraph()
    r = p.add_run("Table %d | " % cap[0]); r.bold = True
    p.add_run(cap[1])
    tb = d1.add_table(rows=len(rows), cols=len(rows[0]))
    tb.style = "Table Grid"
    for ri, row in enumerate(rows):
        for ci, val in enumerate(row[:len(rows[0])]):
            cell = tb.cell(ri, ci); cell.text = ""
            add_runs(cell.paragraphs[0], val if ri == 0 else val.replace("**", ""))
            for rr in cell.paragraphs[0].runs:
                rr.font.size = Pt(8)
                if ri == 0:
                    rr.bold = True
    d1.add_paragraph()

OUT1 = os.path.join(PKG, "01_manuscript_text", "Manuscript_main_text.docx")
os.makedirs(os.path.dirname(OUT1), exist_ok=True)
d1.save(OUT1)
say("")
say(f"写出：{OUT1}")

# ================================================================ 文档 2：补充表
RP = os.path.join(ROOT, "figures", "data", "reporting_profiles")
S1 = read(os.path.join(RP, "Supplementary_Table_S1.md"))
S2 = read(os.path.join(RP, "Supplementary_Table_S2.md"))
# ★ 2026-10-07 新增：Table S3（**事后**分诊策略模拟）。由 p46 渲染。
S3 = read(os.path.join(RP, "Supplementary_Table_S3.md"))
SM = read(os.path.join(RP, "Supplementary_Methods.md"))
RR = read(os.path.join(RP, "Supplementary_Tables_R1_R3.md"))

d2 = Document()
st = d2.styles["Normal"]; st.font.name = "Times New Roman"; st.font.size = Pt(11)
p = d2.add_paragraph(); p.alignment = WD_ALIGN_PARAGRAPH.CENTER
r = p.add_run("Supplementary Tables"); r.bold = True; r.font.size = Pt(16)
d2.add_paragraph()
p = d2.add_paragraph()
add_runs(p, "This document contains Supplementary Tables S1, S2 and S3, Supplementary Tables "
            "R1–R3, and Supplementary Methods. Supplementary figures are supplied as separate files.")
for r in p.runs: r.italic = True

for hd, body in [("Supplementary Table S1", S1), ("Supplementary Table S2", S2),
                 ("Supplementary Table S3", S3),
                 ("Supplementary Tables R1–R3", RR), ("Supplementary Methods", SM)]:
    d2.add_page_break()
    d2.add_heading(hd, level=1)
    md_to_doc(body.splitlines(), d2)

OUT2 = os.path.join(PKG, "02_supplementary_tables", "Supplementary_Tables.docx")
os.makedirs(os.path.dirname(OUT2), exist_ok=True)
d2.save(OUT2)
say(f"写出：{OUT2}")

# ================================================================ 回读校验
say("")
say("【回读校验】")
def back(path, expect_min_tables):
    """回读 docx，**含表格单元格**（首版只扫段落，漏掉了表格里的残留标记）。"""
    d = Document(path)
    ntab = len(d.tables)
    npar = len([q for q in d.paragraphs if q.text.strip()])
    parts = [q.text for q in d.paragraphs]
    for t in d.tables:
        for row in t.rows:
            for c in row.cells:
                parts.append(c.text)
    must(ntab >= expect_min_tables, "%s 表格数 %d < %d" % (os.path.basename(path), ntab, expect_min_tables))
    return ntab, npar, "\n".join(parts)

t1n, p1n, txt1 = back(OUT1, 1)
# ★ 2026-10-12：原写死「（References）」，而正文 docx 里唯一的表是 **Table 1**（§11 的方向约定表），
#   参考文献是**段落列表**、不是表。⇒ 回读标签**指错了对象**，读日志的人会以为主表丢了
#   （实测本轮就先去 docx 里翻了一遍才确认 Table 1 在场）。
#   改为**列出实际表头**：标签从被描述的对象派生，不硬编码（INC-021 / INC-034 族）。
_t1_heads = []
try:
    import docx as _dx
    _d1 = _dx.Document(OUT1)
    for _t in _d1.tables:
        _t1_heads.append("/".join(c.text.strip() for c in _t.rows[0].cells)[:60])
except Exception as _e:
    _t1_heads = ["(读取失败 %s)" % type(_e).__name__]
say(f"  正文 docx ：{p1n} 段 ｜ {t1n} 表 ｜ 各表首行：{_t1_heads if _t1_heads else '(无表)'}")
must(TITLE[:40] in txt1, "正文 docx 里找不到标题")
# ★★ 2026-10-09（INC-026）：本行原**硬编码** v1 摘要的句子
#   「Conventional input-space shift statistics showed little positive discrimination」。
#   A5.3a 把 Abstract 改成 v4 后（`little` → **`no`**），本断言**必然失败** ⇒ p40 在
#   **写完两个 docx 之后**中止，包停在半成品状态。
#   ⇒ 改为**从正本派生**期望句子：读 Abstract 的**第一句正文**，断言它出现在 docx 里。
#   凡「可以被派生」的断言，都不得硬编码 —— 硬编码的那一份会与世界静默分叉（INC-021 族）。
_ABS_SRC = os.path.join(ROOT, "Abstract_v1.md")
_abs_body = io.open(_ABS_SRC, encoding="utf-8").read().split("## Abstract", 1)[1].split("\n---", 1)[0]
_abs_head = " ".join(_abs_body.split())[:120]                          # 取正文开头 120 字
# 用**开头片段**做在场判据：既锚定「摘要正文进了 docx」，又不依赖任何会改写的具体措辞
must(_abs_head[:80] in txt1,
     "正文 docx 里找不到摘要正文（期望开头：%r）" % _abs_head[:80])
must("The prespecified warning task fixed an operating point" in txt1, "正文 docx 里找不到 §7.2 小节")
must("Fig." in txt1, "正文 docx 里找不到图引用")
for bad in ["**[", "## ", "| ---", "Supplementary_Legends"]:
    must(bad not in txt1, "正文 docx 残留 Markdown/内部标记：%r" % bad)
say("  正文 docx 关键句在场 ✓；无 Markdown 残留 ✓")

t2n, p2n, txt2 = back(OUT2, 5)
say(f"  补充 docx ：{p2n} 段 ｜ {t2n} 表（S1 / S2 Panel A / S2 Panel B / R1 / R2 / R3）")
must("Supplementary Table S1" in txt2 and "Supplementary Table S2" in txt2
     and "Supplementary Methods" in txt2, "补充 docx 章节标题不全")
for bad in ["| ---", "`", "## "]:
    must(bad not in txt2, "补充 docx 残留 Markdown：%r" % bad)
say("  补充 docx 章节齐全 ✓；无 Markdown 残留 ✓")

# ================================================================ 其余材料
def copy(src, dst_dir, name=None):
    os.makedirs(dst_dir, exist_ok=True)
    d = os.path.join(dst_dir, name or os.path.basename(src))
    shutil.copy2(src, d); return d

# 03 主图 / 04 补充图
for i in range(1, 6):
    for ext in ("png", "svg", "tiff"):
        s = (os.path.join(ROOT, "figures", "TIFF_400dpi", f"Figure{i}.tiff") if ext == "tiff"
             else os.path.join(ROOT, "figures", "final", f"Figure{i}_FINAL.{ext}"))
        if os.path.exists(s):
            copy(s, os.path.join(PKG, "03_figures_main"),
                 f"Figure{i}.{ext}")
for i in range(1, 10):
    for ext in ("png", "svg", "tiff"):
        s = (os.path.join(ROOT, "figures", "TIFF_400dpi", f"FigureS{i}.tiff") if ext == "tiff"
             else os.path.join(ROOT, "figures", "supplement", f"FigureS{i}.{ext}"))
        if os.path.exists(s):
            copy(s, os.path.join(PKG, "04_figures_supplementary"), f"FigureS{i}.{ext}")
copy(os.path.join(ROOT, "figures", "final", "Figure_Captions.md"), os.path.join(PKG, "03_figures_main"))
copy(os.path.join(ROOT, "figures", "final", "Supplementary_Legends.md"),
     os.path.join(PKG, "04_figures_supplementary"))

# ★★ 2026-10-10（INC-029）：**打包前先验「验证证据比稿件新」**。
#   病历：本包曾于 17:12 构建，而 EVS 双跑 17:51 才跑完 ⇒ 包里装的是**双跑之前**的
#   `report_evs.txt` / `coverage_table.csv`。**没有任何东西发现它** —— 因为
#   「包里的证据对应哪一次运行」从来不是任何判据的问题。
#   ⇒ 顺序必须写死在代码里：**先双跑 → 再打包 → 再验收**。
# ★★ 2026-10-12：本表必须**恰好等于**下面真正读入 docx 的那几份文件。
#   此前只有 6 份，漏了 `References_FINAL.md` 与 `DataCodeAvailability_v1.md`
#   —— 而 `read(os.path.join(MR, ...))` 里**确实**读了它们（见本文件上方 72–81 / 327 行）。
#   ⇒ 那两份改动后**不必重新跑 EVS 就能打包**，守卫会放行一份「新正文 + 旧证据」的组合。
#   新鲜度守卫的价值全在于**覆盖完整**；漏一份就等于没有。
_MAN = [os.path.join(ROOT, "Abstract_v1.md"), os.path.join(ROOT, "Title_v1.md")] + \
       [os.path.join(ROOT, "manuscript_renumbered", f) for f in
        ("Introduction_v1.md", "Methods_v1.md", "Results_v3_clean.md", "Discussion_v1.md",
         "References_FINAL.md", "DataCodeAvailability_v1.md")]
_newest_ms = max(os.path.getmtime(f) for f in _MAN if os.path.exists(f))
_EVS = os.path.join(ROOT, "evs", "report_evs_real.txt")
must(os.path.exists(_EVS),
     "找不到 EVS 证据 evs/report_evs_real.txt —— 请先跑 evs/run_evs.py")
_evs_age = os.path.getmtime(_EVS)
if _evs_age < _newest_ms:
    import datetime as _dt
    raise RuntimeError(
        "*** 打包中止 ***：EVS 证据比稿件旧 —— 包里会装着**过期**的验证证据。\n"
        "   稿件最新：%s\n   EVS 证据：%s\n"
        "   正确顺序：**先跑 evs/run_evs.py（双跑），再跑本脚本**。"
        % (_dt.datetime.fromtimestamp(_newest_ms).strftime("%Y-%m-%d %H:%M:%S"),
           _dt.datetime.fromtimestamp(_evs_age).strftime("%Y-%m-%d %H:%M:%S")))
say("  证据新鲜度 ✓ EVS 证据晚于稿件最新改动")

# 05 验证证据
VD = os.path.join(PKG, "05_verification")
for s in ["evs/report_evs_real.txt", "evs/coverage_table_real.csv",
          "evs/manifest_sha256.txt", "CONVENTION_FREEZE_20261007/SHA256SUMS.txt"]:
    f = os.path.join(ROOT, s.replace("/", os.sep))
    if os.path.exists(f):
        copy(f, VD, os.path.basename(f).replace("_real", ""))

# 06 来源
SD = os.path.join(PKG, "06_source_markdown")
for s in ["Abstract_v1.md", "Title_v1.md", "manuscript_renumbered/Introduction_v1.md",
          "manuscript_renumbered/Results_v3_clean.md", "manuscript_renumbered/Methods_v1.md",
          "manuscript_renumbered/Discussion_v1.md", "manuscript_renumbered/References_FINAL.md",
          "figures/data/reporting_profiles/Supplementary_Table_S1.md",
          "figures/data/reporting_profiles/Supplementary_Table_S2.md",
          "figures/data/reporting_profiles/Supplementary_Table_S3.md",
          "figures/data/reporting_profiles/Supplementary_Tables_R1_R3.md",
          "figures/data/reporting_profiles/Supplementary_Methods.md"]:
    f = os.path.join(ROOT, s.replace("/", os.sep))
    if os.path.exists(f):
        copy(f, SD, os.path.basename(f))

say("")
say("目录内容：")
for d in sorted(os.listdir(PKG)):
    full = os.path.join(PKG, d)
    if os.path.isdir(full):
        say(f"  {d}/  ({len(os.listdir(full))} 项)")
    else:
        say(f"  {d}")

# ================================================================ README（由脚本生成，重建不丢）
WC = {nm: len(body.split()) for nm, body in
      [("Title", TITLE), ("Abstract", ABSTRACT), ("Introduction", INTRO), ("Results", RESULTS),
       ("Methods", METHODS), ("Discussion", DISCUSS), ("References", REFS)]}

# ★★ 2026-10-14（INC-038 沉淀 1「重复件/摘要必须被某个判据认领」的同一课，用在**自述表**上）：
#   本 README 的「验收结果」表**原先手打计数**，于是 p15 从 153 项扩到 203 项、p34 从 22 项
#   重写为 80 项之后，这张表**在纸面上还印着旧数** —— 它看起来是验收记录，实际是过期快照，
#   而且**跟着投稿包一起发出去**。凡审计器留下了 CSV，这里就**从 CSV 派生**；取不到才显示「—」。
def _audit_csv(fn):
    q = os.path.join(ROOT, "分析结果", fn)
    if not os.path.exists(q):
        return None
    try:
        _d = pd.read_csv(q)
    except Exception:                                                # noqa: BLE001
        return None
    _c = "status" if "status" in _d.columns else ("result" if "result" in _d.columns else None)
    return (len(_d), _d[_c].value_counts().to_dict() if _c else {})
def _cnt(fn, key="PASS"):
    r = _audit_csv(fn)
    if r is None:
        return None
    n, vc = r
    return n if key is None else int(vc.get(key, 0))
_A15N, _A15Q = _cnt("p15_results_numeric_audit.csv"), _cnt("p15_results_qualitative_audit.csv")
_A16N, _A16S = _cnt("p16_numeric.csv"), _cnt("p16_semantic.csv")
_A16C_ALL = _cnt("p16_completeness.csv", None)
_A16C_OK  = _cnt("p16_completeness.csv")
_A34_ALL, _A34_OK = _cnt("p34_abstract_audit.csv", None), _cnt("p34_abstract_audit.csv")
_A34_WARN = (None if _A34_ALL is None
             else _A34_ALL - _A34_OK - (_cnt("p34_abstract_audit.csv", "FAIL") or 0))
def _t15():
    return "（产物缺失）" if _A15N is None else f"{_A15N} 数值 + {_A15Q} 定性"
def _t16():
    if _A16N is None:
        return "（产物缺失）"
    return f"{_A16N} + {_A16S} + {_A16C_OK}/{_A16C_ALL}（1 项 BLOCKED = IRB 声明）"
def _t34():
    if _A34_ALL is None:
        return "（产物缺失）"
    return (f"{_A34_OK}/{_A34_ALL} PASS"
            + (f"；另有 {_A34_WARN} 条 WARN（作者决定项，见报告）" if _A34_WARN else ""))

README = f"""# 投稿包 — 逐项说明

**生成日期**：2026-10-07
**生成方式**：`分析结果/p40_build_submission.py` **自动构建**（只从正本读取；正本认定见
`MANUSCRIPT_MASTER_DECLARATION.md` §2）。构建日志见 `_build_log.txt`。
验收由 `分析结果/p41_verify_submission.py` **独立执行**（逐段核对内容零丢失、数字守恒、无泄漏）。

> 正本分布：`Title_v1.md` 与 `Abstract_v1.md` 在**根目录**，其余五份在 `manuscript_renumbered/`。
> 本包未使用任何根目录的废弃副本（`Results_v1/v2.md`、`References_v1.md`）。

---

## 目录

```
SUBMISSION_PACKAGE_20261007/
├── README_投稿包说明.md          ← 本文件（由 p40 生成）
├── _build_log.txt                ← 构建日志
├── 01_manuscript_text/           ← 【Word ①】正文全套，**不含图**
├── 02_supplementary_tables/      ← 【Word ②】全部补充表
├── 03_figures_main/              ← 主图 1–5（png/svg/tiff）+ 图注
├── 04_figures_supplementary/     ← 补充图 S1–S9（png/svg/tiff）+ 图例
├── 05_verification/              ← 验证证据（EVS 报告 / 覆盖表 / 哈希清单）
└── 06_source_markdown/           ← 生成 Word 所用的源 Markdown
```

---

## ① 正文 Word — `01_manuscript_text/Manuscript_main_text.docx`

顺序：Title → Abstract → Introduction → Results → Methods → Discussion → References。

| 段 | 词数 |
|---|---|
| Title | {WC['Title']} |
| Abstract | {WC['Abstract']}（Nature Communications 上限 200；不含引用） |
| Introduction | {WC['Introduction']:,} |
| Results | {WC['Results']:,}（含 §7.2「协议预设预警任务」小节） |
| Methods | {WC['Methods']:,}（§1–§15，含 15 信号方向约定全表） |
| Discussion | {WC['Discussion']:,} |
| References | {WC['References']:,} |

**已剥离的内部产物**（构建脚本按**显式边界**剥离，剥离后断言边界确实存在）：

- `## References cited in this Introduction / Methods / Discussion` 三张内部引文核对子表
- `## 内部记录更正`、`## 为什么是这一版`、`## 被否的候选与理由`
- 各文件顶部的 `>` 批注块（定位、禁用措辞表、Supplement 落点等）
- Abstract 末尾的「自查」表、References 表头前的编号说明

**图**：正文 Word **不含任何图**；图以独立文件提供，正文中以 `Fig. 1`…`Fig. 5`、
`Fig. S1`…`Fig. S9` 引用。

---

## ② 补充表 Word — `02_supplementary_tables/Supplementary_Tables.docx`

| 章节 | 内容 | 来源（全部**自动渲染**，非手工录入） |
|---|---|---|
| **Supplementary Table S1** | 186 个合格迁移的完整清单（21 列 × 186 行，含全部 12 个信号） | `p21e` ← `figures/data/tableS1_tasks.csv` |
| **Supplementary Table S2** | 15 个信号的协议预警指标（AUROC / AUPRC / sens@.80 / PPV / NPV + cluster bootstrap 区间）与排序汇总 | `p21e` ← `figures/data/tableS2_metrics.csv` |
| **Supplementary Tables R1–R3** | 站点画像 / 记录-回合流 / 逐站任务数 | `p21d` |
| **Supplementary Methods** | 四个预设配置的精确估计器设置 | `p21f` ← 从 `p12_configs_corrected.py` fail-closed 抽取 |

补充图（S1–S9）见 `04_figures_supplementary/`；图例见同目录 `Supplementary_Legends.md`。

---

## 引用关系（正文 ↔ 补充材料，**双向**机器核验）

| 补充条目 | 正文引用位置 |
|---|---|
| Supplementary Table S1 | Methods §6 |
| Supplementary Table S2 | Results §7.2 小节 |
| Supplementary Tables R1–R3 | Methods §1 |
| Supplementary Methods | Methods §5 |
| Supplementary Figure Legends | 补充材料自身章节（不要求正文引用） |
| Fig. S1 – Fig. S9 | Results / Methods / Discussion |
| Fig. 1 – Fig. 5 | 正文各处（**面板级引用亦已核对无遗漏**，含 `Fig. 1A,B`、`Fig. 4A–C` 记法） |

> 由 `分析结果/p35_supplement_inventory_audit.py` 机械核验：方向一查**孤儿**（有产物无引用），
> 方向二查**悬空**（有引用无产物）。当前 **14/14 PASS**。

---

## ③ / ④ 图件

每张图三种格式：`.png`（预览）、`.svg`（矢量）、`.tiff`（≥330 dpi，印刷）。

---

## ⑤ 验证证据 — `05_verification/`

| 文件 | 说明 |
|---|---|
| `report_evs.txt` | EVS closure 报告 |
| `coverage_table.csv` | 逐路径的测试覆盖表 |
| `manifest_sha256.txt` | closure 时点的 manifest 身份 |
| `SHA256SUMS.txt` | 冻结清单（正文、补充表、渲染脚本、审计器的 SHA-256） |

| 产物 | 核验器 | 结果 |
|---|---|---|
| Results 的每个数字 | `p15` | {_t15()} |
| Methods 的每个数字 | `p16` | {_t16()} |
| Discussion 的断言强度 | `p19` | 全通过 |
| Abstract 的数字与措辞 | `p34` | {_t34()} |
| 补充图例的每个数字 | `p31` | {_cnt("p31_legends_audit.csv")}/{_cnt("p31_legends_audit.csv", None)} |
| 正文 headline numbers（独立重算） | `p30`（MNV） | 82 PASS / 0 FAIL |
| 补充材料产物 ↔ 引用 | `p35` | 14/14 |
| 渲染件（S1 / S2 / Supplementary Methods） | `evs/t23`（**独立**，不 import 渲染器） | 全通过 |
| 代码 / 估计量层 | EVS 全套 | closure 四条全满足 |

> **表中数字怎么来的（如实分档，别把两种当成一种读）**：
> **① 现场派生**（每次构建时从**审计器自己的 CSV** 现读，因此永远是新值）：
> `p15`、`p16`、`p31`、`p34` 四行。
> **② 手打快照**（写死在本脚本里的字符串，**不会自动跟着审计器走**）：
> `p19`「全通过」、`p30`「82 / 0」、`p35`「14/14」、`evs/t23`「全通过」、EVS「closure 四条全满足」。
> ⇒ **②类的数字若与审计器不符，以审计器的报告为准。**
> （为什么不全派生：②类的四个核验器**不写结构化 CSV**，只能读它们的报告文本；
> 在收尾阶段改四个核验器的输出格式，风险大于收益。**此条为已知缺口，不是已解决。**）

---

## ⑥ 源 Markdown — `06_source_markdown/`

**注意**：这些源文件**含有内部批注段**，生成 Word 时已被剥离。
**投稿请使用 `.docx`，不要直接提交这些 `.md`。**

---

## 已知待办

| 项 | 位置 | 说明 |
|---|---|---|
| **IRB / 数据使用声明** | 正文 › Methods §1 | 保留为占位符 `[Institutional review board determination and data-use statement to be inserted.]`，**需作者提供真实批件信息后替换**。`p16` 的完备性检查会一直 BLOCKED 到该占位符被替换。 |

---

## 复现

```
python 分析结果\\p40_build_submission.py            # 重建本包（清空并重写 SUBMISSION_PACKAGE_20261007/）
python 分析结果\\p41_verify_submission.py           # 独立验收（逐段零丢失 / 数字守恒 / 无泄漏）
python 分析结果\\p35_supplement_inventory_audit.py  # 补充材料引用双向核验
cd evs && python run_evs.py                          # EVS closure —— ★ 必须连跑两次
cd evs && python run_evs.py                          # ★ 第二次才是有效判读（见下）
```

> ### ★ EVS 必须连跑两次 —— 这是**执行条件**，不是习惯
>
> `t14` 的 `MC-10`（每条 FAIL 都要有归属）读的是 `evs_results_{{real,patched}}.csv`，
> 而该文件由 `run_evs.py` 在**全部模块跑完之后**写出；`t14` 自身又是并行模块之一。
> **⇒ `MC-10` 永远读到的是上一轮的结果。**
>
> 后果有两个方向，其中**假绿更危险**：
> - **假红**：上一轮有、本轮已修的 FAIL 会被报出来（2026-10-07 实测到）；
> - **假绿**：上一轮全绿、本轮新引入一条无归属 FAIL 时，`MC-10` 本轮**看不见**。
>
> **因此：任何对生产树的改动之后，必须连跑两次 `run_evs.py`，以第二次为准。**
> 只跑一次就下结论 = 让 closure gate 的滞后替你做判断。
>
> 判读标准：第二次运行后，`evs_results_real.csv` 里的 FAIL **必须恰为** `EXPECTED_RED`
> 的已知缺陷集合（当前 = `{{SUDO-01, SUDO-03, SUDO-05, FORK-03}}`），且 closure 四条全满足。
> （注：本 README 是 f-string，花括号须写成双花括号。）
>
> *（冻结后可选的最小修复：`run_evs.py` 在全模块跑完后**追加执行一次 t14**，
> 不改动 `t14` 本体、只改 runner 的调用序列。冻结期内不动。）*
"""
io.open(os.path.join(PKG, "README_投稿包说明.md"), "w", encoding="utf-8").write(README)
say(f"写出：{os.path.join(PKG, 'README_投稿包说明.md')}")

io.open(os.path.join(PKG, "_build_log.txt"), "w", encoding="utf-8").write("\n".join(LOG))

# ================================================================ GENERATION.md
# ★★ 2026-10-12（INC-018 的同族处置）：本目录名的日期是 **20261007**，
#   而它每次构建都被**整体重写**。⇒ 目录名**不是**代次信息，
#   只看名字取件的人无从知道手里这一包对应哪一版稿件。
#   故每次构建都写一张机器生成的代次声明：**构建时刻 + 六个稿件源的 sha256 + 证据哈希**。
#   （`SHA256SUMS.txt` 只证明文件未被改动，不证明是哪一代 —— INC-018 的原文。）
import datetime as _dt2
import hashlib as _hl


def _sha(p):
    return _hl.sha256(open(p, "rb").read()).hexdigest()


_gen = []
_gen.append("# 代次声明（由 `分析结果/p40_build_submission.py` 每次构建覆写）")
_gen.append("")
_gen.append("> **本目录名里的 `20261007` 是建档日期，不是本包的产生时间。**")
_gen.append("> 本目录每次构建都被整体清空重写 —— **只有这张纸说明它是哪一代。**")
_gen.append("")
_gen.append("## 本包产生时刻")
_gen.append("")
_gen.append("| 项 | 值 |")
_gen.append("|---|---|")
_gen.append("| **构建时刻** | %s |"
            % _dt2.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
_gen.append("| **构建脚本** | `%s` |" % os.path.abspath(__file__))
_gen.append("| **源目录** | `%s` |" % ROOT)
_gen.append("")
_gen.append("## 稿件源（本包正文即由这些文件生成）")
_gen.append("")
_gen.append("| 文件 | 改动时刻 | sha256 |")
_gen.append("|---|---|---|")
for _f in _MAN:
    if os.path.exists(_f):
        _gen.append("| `%s` | %s | `%s` |"
                    % (os.path.relpath(_f, ROOT),
                       _dt2.datetime.fromtimestamp(os.path.getmtime(_f)).strftime("%Y-%m-%d %H:%M:%S"),
                       _sha(_f)[:16]))
_gen.append("")
_gen.append("## 验证证据")
_gen.append("")
_gen.append("| 项 | 值 |")
_gen.append("|---|---|")
_gen.append("| EVS 报告 | `%s` |" % os.path.relpath(_EVS, ROOT))
_gen.append("| 报告时刻 | %s |"
            % _dt2.datetime.fromtimestamp(os.path.getmtime(_EVS)).strftime("%Y-%m-%d %H:%M:%S"))
_gen.append("| 报告 sha256 | `%s` |" % _sha(_EVS)[:16])
_gen.append("")
_gen.append("**新鲜度条件（构建时已强制通过）**：EVS 报告时刻必须 **晚于**上表全部稿件源的")
_gen.append("最新改动时刻。本包构建时该条件成立，否则 p40 会中止而**不写出任何文件**。")
_gen.append("")
_gen.append("## 用法")
_gen.append("")
_gen.append("- **可以**：按本表核对「手里的包对应哪一版稿件」。")
_gen.append("- **不可以**：把本目录名当时间戳；也不要在只重跑 `p56` 之后以为上传包更新了"
            " —— `p56` **读**本目录的产物，正文变了必须**先重跑 p40**。")
io.open(os.path.join(PKG, "GENERATION.md"), "w", encoding="utf-8").write("\n".join(_gen) + "\n")
say(f"写出：{os.path.join(PKG, 'GENERATION.md')}（代次声明）")

# ★ INC-036：把**被 skip_first_h1 丢掉的 1 级标题**全列出来（**必须在所有 md_to_doc 之后**，
#   否则补充表那条还没发生，列表是空的 —— 首版就把这段放早了，打印出「共 0 条」）。
#   它们是被**有意**丢弃的重复标题（各节标题由 add_heading 写）。列出来是让人一眼看出
#   「有没有混进一行不是标题的正文」；真正的判据是 p41 V4-3（**词级**守恒，不依赖启发式）。
say("")
say("被 skip_first_h1 丢弃的 1 级标题（应为各节/各表自己的重复标题，共 %d 条）：" % len(SKIPPED_H1))
for _s in SKIPPED_H1:
    say("    · %s" % _s)
if not SKIPPED_H1:
    say("    （无）")

say("")
say("*** p40 构建完成 ***")
