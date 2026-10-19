# -*- coding: utf-8 -*-
r"""p60 — 投稿包**泄漏扫描闸门**（recursive leak gate）

================================ 为什么有这个脚本 ================================
外部核验（2026-10-12）指出一个反复出现的失效模式，并且我同意：
**「查了一个文件集就当成全查。」** 本会话里它已经犯了至少三次 ——
`pre-registered` 只扫了 `manuscript_renumbered/`（图注与 S1/S3 漏网）、
`§16` 的指针没打开被指向的 S2、把自己文件里的哈希行说成是别人转储带来的。

⇒ 需要一道**递归、按格式、含文件名**的闸门，而不是再靠人记得扫哪里。

================================ 扫什么 ================================
  ① **文件名本身**（`新建 文本文档.txt`、`p15_*.py` 这类）
  ② **纯文本**（.md/.txt/.csv/.json/.R）
  ③ **docx** —— 解包读 `word/document.xml`（批注与页眉页脚也读）
  ④ **xlsx** —— 解包读 `xl/sharedStrings.xml` + 工作表名
  ⑤ **SVG** —— 是 XML，按文本读
  （TIFF/PNG 为二进制，只入文件名检查；其内容由 p31 的数值审计负责）

================================ 扫描项（零命中才算过）========================
    CJK ｜ sha256/Locked./.gen.json ｜ audit log ｜ 脚本名 p\d+_ ｜ 绝对路径
    pre-registered ｜ superseded ｜ corrected ｜ 占位符 ｜ "statement to make"
    GPT6/GPT-6 ｜ 内部目录名（分析结果/、figures/data/、manuscript_renumbered/）

================================ 非空转的证据（本脚本自带的负向测试）============
`--negtest <dir>` 对**指定的一个目录**跑一遍并打印命中。用本仓 `_artifact_archive/` 里
那些**修改前**的包即可证明闸门有判别力 —— 判据能对真实的脏版本判红，才不是空转。
"""
import io
import os
import re
import sys
import hashlib          # ★ 2026-10-15：「包 vs 源」逐位比对用（首版补丁漏了导入 —— 补丁里那句
                       #   `s.replace("import io, os, re, sys", …)` 根本没匹配上，因为本文件的
                       #   import 是**一行一个**。**静默 no-op 的补丁**等于没打。）
import zipfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEFAULT = os.path.join(ROOT, "SUBMISSION_20261010")

PATTERNS = [
    ("CJK",              r"[\u4e00-\u9fff]"),
    ("sha256/Locked",    r"(?i)\bsha256\b|\bLocked\.|\.gen\.json"),
    # ★ 2026-10-13 收紧：Methods §15 **正面描述**本仓的审计做法（"a dated analysis audit log
    #   records every analysis decision…"）——那是**方法上的披露，不是内部引用**。
    #   要抓的是**指向某条目**的引用（`见 audit log A-2026-…`），不是这个短语本身。
    ("audit log 引用",    r"(?i)audit[ _]log\s+A-\d{4}|见\s*audit log|见审计日志"),
    # ★ 2026-10-13 收紧：首版写成 `\bp\d{1,3}[a-z]?_[a-z]`，把**列名**也当成脚本名
    #   （`p90_minus_p10`、`p10_aligned`、`p6d_loo_org_auc` 全是数据列）——
    #   与早先「先用后导」扫描同一类的假阳性。**要求带文件扩展名**才判。
    ("脚本名 p\\d+_",     r"\bp\d{1,3}[a-z]?_[A-Za-z0-9_\-]+\.(?:py|csv|R|r|txt|json)"),
    ("绝对路径",          r"[A-Za-z]:\\\\|[A-Za-z]:\\|分析结果/|figures/data/|manuscript_renumbered/"),
    ("pre-registered",   r"(?i)pre-?registered"),
    # ★ 2026-10-13 收紧：`superseded` / `corrected` **当裸词用太宽** —— Methods §15 正面写着
    #   「a dated analysis audit log records … including use of a superseded configuration
    #   source」，那是**方法披露**，不是残留。要抓的是**图注里的修订残留**（具体搭配）。
    ("修订残留",          r"Corrected MMD|superseded denominator|after correction"),
    ("占位符",            r"TO BE COMPLETED BY THE AUTHORS|to be inserted by the authors|\[TODO|待填"),
    ("statement to make", r"(?i)statement to make"),
    ("GPT6",             r"(?i)\bGPT[- ]?6\b"),
]

# ★ 2026-10-14：**Word 自身的排版部件不算内容。**
#   `word/fontTable.xml` 与 `word/theme/theme1.xml` 里的中文是**字体名**（宋体 / ＭＳ 明朝），
#   由 Word 在保存时自动写入，**与作者写了什么无关** —— 无法通过编辑正文消除。
#   ⇒ 按**成员名**豁免（不是按文件豁免：docx 的正文 `word/document.xml` **照扫**）。
_WORD_CHROME = ("word/fontTable.xml", "word/theme/theme1.xml")

# ★ 2026-10-14：**IRB / 数据使用声明只能由作者提供**，正文里留的是**明写的占位符**。
#   它是**已知且被别的判据盯着的**项（`p16` 的 C-01 = BLOCKED），不是意外泄漏
#   ⇒ 内容级豁免**仅限这一句原文**，并把它单列进 `DISCLOSE`，
#     保证「豁免」不等于「看不见」（DISCLOSE 会打印在报告里，证明判据查过它）。
#   ★ 2026-10-14 首版按**记忆**写成 `[IRB determination and data-use statement…]` —— 与正文
#     实际措辞不符（正文是 `[Institutional review board determination and data-use statement
#     to be inserted by the authors.]`），于是**豁免没生效**、命中照旧。
#     「凭记忆写 pattern」是 INC-038 沉淀 2 的同族 ⇒ 改为**对正文实际两种写法都覆盖**。
AUTHOR_FILL = [
    (r"\[(?:IRB|Institutional review board)[^\]]{0,80}to be inserted by the authors\.\]",
     "IRB / 数据使用声明（`p16` C-01 = BLOCKED）"),
    (r"\[TO BE COMPLETED BY THE AUTHORS[^\]]*\]",
     "代码可用性：仓库 URL / Zenodo DOI / 许可"),
    (r"\[TODO[^\]]*\]", "作者自留 TODO"),
    # ★ 2026-10-14：IRB 段里还有**一句同批的散文**（不是占位符）：
    #   「…is the authors' statement to make.」—— 它说明**该判定的措辞属作者**，
    #   与上面的 `[Institutional review board…]` 占位符是**同一件事**，只能由作者落实。
    #   `PATTERNS` 的 `statement to make` 命中它属**假阳性**（那条抓的是内部指令的形状），
    #   故移入本桶：既不判红，也在「投稿前必须为 0」的清单里可见。
    (r"(?i)is the authors['’]{0,2} statement to make\.",
     "IRB 段：判定措辞归作者（与上面的占位符同批）"),
]
IRB_PLACEHOLDER = re.compile(AUTHOR_FILL[0][0])

# 允许的例外：出现在这里的内容**不是泄漏**，且理由必须写明（不得用通配掩盖）
DISCLOSE = []   # 「查过且干净」的披露（证明判据不是空转）

ALLOW = [
    # 上传包里 `00_作者待填_★先看这里/` 是**给作者的空骨架**，本来就含中文说明与占位符。
    # 但它**不该随投稿 zip 上传** —— 故本脚本把它单列为一条「须人工处置」，而不是静默放行。
    (r"^00_作者待填_★先看这里/", "作者待填骨架（含中文说明与占位符）——上传 zip 时应排除"),
    (r"^README_投稿包说明\.md$", "给作者看的包说明（非投稿件）"),
    (r"^GENERATION\.md$", "代次声明（本地记账，非投稿件）"),
    # ★★ 2026-10-14（**显式决定，不是静默放过**）：**代码归档内的中文**（`08_代码归档/`
    #   的脚本注释与 `MANIFEST.csv` 的 `role` 列说明）**判为允许**。
    #   · 理由一：它**不是隐私或身份泄漏** —— 不含姓名、机构、路径、患者信息；
    #   · 理由二：外部核验对同类（脚本注释）的判定是「不算隐私泄漏…可选处理」，已照此；
    #   · 理由三：**路径**那一项已单独修掉（`source` 列改为英文相对布局，无绝对路径）。
    #   ★ 代价如实记账：开放仓库的读者会看到中英混排，观感上是「半成品」。
    #     若日后要清，做法是**翻译注释**（改动代码文件 ⇒ 须重跑已知答案测试）。
    (r"^08_代码归档/", "代码归档：中文注释与 role 说明 —— 允许（非隐私泄漏；路径项已单独修）"),
    # ★ 2026-10-14：**可运行画图脚本**同理 —— 脚本注释是中文，非隐私泄漏。
    #   （首版只豁免了文件名扫描，没豁免内容 ⇒ 报了 32 条；这是判据的缺口，已补。）
    (r"^07_画图代码_可运行/", "画图脚本：中文注释 —— 允许（与 08_代码归档 同例）"),
    # ★ 2026-10-14：`09_验证证据_可选/` 是**投不投由作者定**的可选目录，其内容（EVS 报告、
    #   覆盖表）天然含测试编号与脚本名。**它不该进上传 zip** —— 故内容不扫、
    #   但提升为「须排除」的显式判定（见下面的清单检查与 ⚠ 行）。
    #   ★★ 上面这段注释从写下起就**没有对应的条目** —— 注释说「内容不扫」，实际照扫，
    #     于是每次都以 9 条命中把闸门判红。**注释写了、代码没写**（与 INC-039 同族：
    #     「写下规矩 ≠ 规矩被执行」）。现补齐，并保留「须排除」的显式提醒。
    (r"^09_验证证据_可选/", "可选验证证据（EVS 报告/覆盖表，含测试编号）——**不上传 zip**"),
]


def _read_members(zf, names):
    out = []
    for n in names:
        try:
            out.append((n, zf.read(n).decode("utf-8", "ignore")))
        except Exception:
            pass
    return out


def _office_meta(path):
    """★ 2026-10-13（外部核验要求）：**文档属性、修订痕迹、批注、隐藏工作表** ——
       这些靠正文文本正则扫不到，但它们携带作者名、公司名、本机路径与改动历史。
       返回 [(标签, 文本)]。"""
    out = []
    low = path.lower()
    try:
        with zipfile.ZipFile(path) as z:
            names = z.namelist()
            if low.endswith((".docx", ".xlsx", ".xlsm", ".pptx")):
                # ① 文档属性（作者 / 最后保存者 / 公司 / 模板路径）
                for n in names:
                    if n.startswith("docProps/") and n.endswith(".xml"):
                        out += _read_members(z, [n])
                        # ★ 2026-10-13：把**查到的属性字段**一并披露 —— 只报「0 命中」
                        #   看不出这条判据到底看没看（本会话反复出现的「空转」形态）。
                        if n.endswith("core.xml"):
                            t = z.read(n).decode("utf-8", "ignore")
                            f = {}
                            for tag in ("creator", "lastModifiedBy", "title", "subject",
                                        "keywords", "description"):
                                # 用 [^<]* 而非 .*? —— 后者会跨标签吞掉下一个字段（实测把
                                # title 读成了 "<dc:subject/><dc:creator>python-docx"）。
                                m = re.search(r"<[a-z:]*%s[^>]*>([^<]*)</" % tag, t)
                                f[tag] = (m.group(1).strip() if m else "")
                            DISCLOSE.append("%s :: %s" % (os.path.basename(path),
                                                          " ｜ ".join("%s=%r" % (k, v)
                                                                      for k, v in f.items())))
                # ② 修订痕迹与批注
                for n in names:
                    if n.endswith(("comments.xml", "people.xml")) or "revision" in n.lower():
                        out += _read_members(z, [n])
                # ③ docx 的 w:ins / w:del（tracked changes）
                if low.endswith(".docx") and "word/document.xml" in names:
                    x = z.read("word/document.xml").decode("utf-8", "ignore")
                    nrev = x.count("<w:ins ") + x.count("<w:del ")
                    if nrev:
                        out.append(("修订痕迹 tracked-changes", "<w:ins>/<w:del> ×%d" % nrev))
                # ④ xlsx 的隐藏工作表
                if low.endswith((".xlsx", ".xlsm")) and "xl/workbook.xml" in names:
                    wb = z.read("xl/workbook.xml").decode("utf-8", "ignore")
                    hid = re.findall(r'<sheet[^>]*name="([^"]+)"[^>]*state="hidden"', wb)
                    hid += re.findall(r'<sheet[^>]*state="hidden"[^>]*name="([^"]+)"', wb)
                    if hid:
                        out.append(("隐藏工作表", "、".join(hid)))
    except Exception:
        pass
    return out


def texts_of(path):
    """返回 [(标签, 文本)]；能解包的按成员展开，其余按文本读。"""
    low = path.lower()
    out = [(os.path.basename(path), "")]
    try:
        if low.endswith(".docx"):
            with zipfile.ZipFile(path) as z:
                mem = [n for n in z.namelist()
                       if n.startswith("word/") and n.endswith(".xml")]
                out = _read_members(z, mem) or out
                out += _office_meta(path)
        elif low.endswith((".xlsx", ".xlsm")):
            with zipfile.ZipFile(path) as z:
                mem = [n for n in z.namelist()
                       if n.endswith(".xml") and ("sharedStrings" in n or "worksheets/" in n
                                                  or "workbook" in n)]
                out = _read_members(z, mem) or [(n, "") for n in z.namelist()]
                out += [(os.path.basename(path) + "::sheetnames",
                         "\n".join(z.namelist()))]
                out += _office_meta(path)
        elif low.endswith((".md", ".txt", ".csv", ".json", ".r", ".svg", ".tsv", ".yml", ".yaml")):
            with io.open(path, encoding="utf-8", errors="ignore") as f:
                out = [(os.path.basename(path), f.read())]
    except Exception as e:
        out = [(os.path.basename(path), "<UNREADABLE %s>" % type(e).__name__)]
    return out


def scan(scan_root, max_show=6):
    hits = []
    author_hits = []          # 「作者必须填」的占位符：不是泄漏，是待办
    n_files = 0
    for dirpath, dirnames, filenames in os.walk(scan_root):
        rel_dir = os.path.relpath(dirpath, scan_root).replace("\\", "/")
        rel_dir = "" if rel_dir == "." else rel_dir + "/"
        for fn in filenames:
            n_files += 1
            full = os.path.join(dirpath, fn)
            rel = rel_dir + fn
            allowed = next((why for pat, why in ALLOW if re.search(pat, rel)), None)
            # ★ 2026-10-13：**代码目录里的脚本文件名不算泄漏** —— 那两个目录存在的意义
            #   就是随包交付脚本（Zenodo 归档 / 可运行 R 脚本）。文件名检查对它们豁免；
            #   **内容检查照旧**（数据包里若嵌了脚本名，那才是泄漏）。
            _code_dir = bool(re.match(r"^(08_代码归档|07_画图代码_可运行/R_scripts)/", rel))
            # ① 文件名
            # ★ 2026-10-14 修：ALLOW 此前**只挡内容扫描、没挡文件名扫描** ——
            #   于是 `00_作者待填/01_Title_page_模板.md` 这类**被判允许的路径**，
            #   其**中文文件名**仍被逐条报出（实测 50 条里大半是这么来的，属判据自身的缺口）。
            #   语义上：一条路径被判为「允许」，指的是**这条路径整体不再审**。
            if not _code_dir and not allowed:
                for name, pat in PATTERNS:
                    if re.search(pat, fn):
                        hits.append((rel, "文件名", name, fn[:70]))
            # ② 内容
            if not allowed:
                for lab, txt in texts_of(full):
                    if not txt:
                        continue
                    if lab in _WORD_CHROME:
                        continue                      # Word 自带字体表/主题：不是作者内容
                    # ★★ 「必须由作者填」的占位符：**既不算泄漏、也不静默豁免**。
                    #   它们是**投稿前的待办**（漏填就直接投出去了）⇒ 单列一个桶：
                    #     · 不进 `hits`（闸门只对**泄漏**判红；否则红灯永远亮着，没人再看）；
                    #     · 逐条打印成**交付清单**，末尾给计数 —— **0 才算真能投**。
                    for _p, _why in AUTHOR_FILL:
                        for _m in re.finditer(_p, txt):
                            author_hits.append((rel, lab, _why, _m.group(0)[:88]))
                        txt = re.sub(_p, " ", txt)
                    for name, pat in PATTERNS:
                        ms = re.findall(pat, txt)
                        if ms:
                            m = re.search(pat, txt)
                            seg = txt[max(0, m.start() - 40):m.end() + 40].replace("\n", " ")
                            hits.append((rel, lab, name, "%d 处｜…%s…" % (len(ms), seg[:90])))
                            break
    return n_files, hits, author_hits


# ================================================================ 【包 vs 源】新鲜度与同一性
# ★★ 2026-10-15（外部核验后补）：此前**没有任何判据**回答这个问题 ——
#   「现在包里的这一份，是不是**现在这批源**做出来的？」
#     · `p41` 的 `CHECKS` 对图件 / 图注 / 图例**只查「文件存在」**，不比对内容；
#     · `p40` 的新鲜度守卫只盯 **EVS 证据**，不盯包本身。
#   ⇒ 后果：**改了源却没重建包**、或**只重建了两个包中的一个**，都不会有红灯。
#     而作者接下来必然要做的事（填 4 处占位符 → 改正文 → 重建两个包）**正好触发这个场景**；
#     `DataCodeAvailability` 里那句「只重跑 p56 会打出旧的正文」此前**只靠人记住**。
#
#   判据是三条**事实**，不是时间戳代理：
#     A 拷贝件   —— 包内那份与源文件**逐位相同**（内容相同 ⇒ 无所谓新旧）；
#     B 派生物   —— docx 由 md 渲染，无法逐位比 ⇒ 要求**源不得比它新**；
#     C 两包之间 —— 同一件交付物在 `SUBMISSION_PACKAGE_20261007/` 与 `SUBMISSION_20261010/`
#                   必须**逐位相同**（专抓「只重建了一个包」）。
#   D 临时/锁文件（`~$*` / `.DS_Store` / `Thumbs.db`）：**只列不判** ——
#     那是 Word / 系统生成的瞬时残渣，不是包的内容缺陷，但**上传前必须不在**。
#
#   ★ 映射是**机械的**（由 `p40` 的拷贝语句直接读出，不是手抄的文件清单），
#     且每条都已在本机**实测逐位相同**后才写进来。
_PKG = os.path.join(ROOT, "SUBMISSION_PACKAGE_20261007")
_UP  = os.path.join(ROOT, "SUBMISSION_20261010")
_PROF = os.path.join(ROOT, "figures", "data", "reporting_profiles")

_SRC_OF = {
    ".tiff": lambda base, ext: os.path.join(ROOT, "figures", "TIFF_400dpi", base + ext),
    ".png":  lambda base, ext: (os.path.join(ROOT, "figures", "final", base + "_FINAL" + ext)
                                if base.startswith("Figure") and not base.startswith("FigureS")
                                else os.path.join(ROOT, "figures", "supplement", base + ext)),
    ".svg":  lambda base, ext: (os.path.join(ROOT, "figures", "final", base + "_FINAL" + ext)
                                if base.startswith("Figure") and not base.startswith("FigureS")
                                else os.path.join(ROOT, "figures", "supplement", base + ext)),
}
_MD_SRC = {
    "Figure_Captions.md":       os.path.join(ROOT, "figures", "final", "Figure_Captions.md"),
    "Supplementary_Legends.md": os.path.join(ROOT, "figures", "final", "Supplementary_Legends.md"),
}

# ★★ 2026-10-15（修我自己的错）：`06_source_markdown/` 的源**必须按显式表给**。
#   首版写成「先试根目录、再试 `manuscript_renumbered/`」的**猜测式**查找 —— 而
#   **根目录恰好有四个旧的「非正本副本」**（`Discussion_v1.md` / `Introduction_v1.md` /
#   `Methods_v1.md` / `Results_v3_clean.md`，见 `MANUSCRIPT_MASTER_DECLARATION.md` §2），
#   于是它比了**错的对象**，报出 4 条假红。`p40` 是**按显式表**拷的。
#   ⇒ 沉淀：**「猜源」在一个存在多代次副本的仓里必然出错**；源必须是**登记过的**，
#     猜不到就**报未登记**，绝不静默挑一个能打开的。
_SM_SRC = {
    "Abstract_v1.md":                "Abstract_v1.md",
    "Title_v1.md":                   "Title_v1.md",
    "Introduction_v1.md":            "manuscript_renumbered/Introduction_v1.md",
    "Methods_v1.md":                 "manuscript_renumbered/Methods_v1.md",
    "Results_v3_clean.md":           "manuscript_renumbered/Results_v3_clean.md",
    "Discussion_v1.md":              "manuscript_renumbered/Discussion_v1.md",
    "References_FINAL.md":           "manuscript_renumbered/References_FINAL.md",
    "Supplementary_Table_S1.md":     "figures/data/reporting_profiles/Supplementary_Table_S1.md",
    "Supplementary_Table_S2.md":     "figures/data/reporting_profiles/Supplementary_Table_S2.md",
    "Supplementary_Table_S3.md":     "figures/data/reporting_profiles/Supplementary_Table_S3.md",
    "Supplementary_Tables_R1_R3.md": "figures/data/reporting_profiles/Supplementary_Tables_R1_R3.md",
    "Supplementary_Methods.md":      "figures/data/reporting_profiles/Supplementary_Methods.md",
}


def _sha(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()


def _pkg_pairs():
    """产出 (包内文件, 源文件, 说明) —— 逐位相同类。"""
    out = []
    for d, tag in ((os.path.join(_PKG, "03_figures_main"), "主图"),
                   (os.path.join(_PKG, "04_figures_supplementary"), "补充图")):
        if not os.path.isdir(d):
            continue
        for fn in sorted(os.listdir(d)):
            if fn in _MD_SRC:
                out.append((os.path.join(d, fn), _MD_SRC[fn], tag + "/图注"))
                continue
            base, ext = os.path.splitext(fn)
            f = _SRC_OF.get(ext.lower())
            if f:
                out.append((os.path.join(d, fn), f(base, ext), tag))
    for d, src in ((os.path.join(_UP, "05_补充材料"), _PROF),):
        if os.path.isdir(d):
            for fn in sorted(os.listdir(d)):
                if fn.endswith(".md"):
                    out.append((os.path.join(d, fn), os.path.join(src, fn), "补充材料"))
    d6 = os.path.join(_PKG, "06_source_markdown")
    if os.path.isdir(d6):
        for fn in sorted(os.listdir(d6)):
            rel = _SM_SRC.get(fn)
            if rel is None:
                # ★ 未登记的源 ⇒ **不猜**（见下），交给上层报「找不到源」
                out.append((os.path.join(d6, fn), os.path.join(ROOT, "__UNMAPPED__", fn),
                            "源 markdown（未登记）"))
            else:
                out.append((os.path.join(d6, fn), os.path.join(ROOT, rel.replace("/", os.sep)),
                            "源 markdown"))
    return out


_MAN_SRC = [os.path.join(ROOT, x.replace("/", os.sep)) for x in
            ("Abstract_v1.md", "Title_v1.md",
             "manuscript_renumbered/Introduction_v1.md", "manuscript_renumbered/Methods_v1.md",
             "manuscript_renumbered/Results_v3_clean.md", "manuscript_renumbered/Discussion_v1.md",
             "manuscript_renumbered/References_FINAL.md",
             "manuscript_renumbered/DataCodeAvailability_v1.md")]
_SUPP_SRC = [os.path.join(_PROF, x) for x in
             ("Supplementary_Table_S1.md", "Supplementary_Table_S2.md",
              "Supplementary_Table_S3.md", "Supplementary_Tables_R1_R3.md",
              "Supplementary_Methods.md")]


def check_pkg_fresh():
    """返回 (违规列表, 提示列表)。违规进 hits（判红），提示只列。"""
    bad, notes = [], []
    pairs = _pkg_pairs()
    n_ok = 0
    for dst, src, tag in pairs:
        if not os.path.exists(dst):
            bad.append(("(包 vs 源)", os.path.relpath(dst, ROOT), "缺件",
                        "包内缺 %s" % tag))
        elif not os.path.exists(src):
            bad.append(("(包 vs 源)", os.path.relpath(dst, ROOT), "找不到源",
                        "包里有、源里没有 ⇒ 映射已失效，须人工核对"))
        elif _sha(dst) != _sha(src):
            bad.append(("(包 vs 源)", os.path.basename(dst), "与源不符",
                        "** %s 与源逐位不符 —— 改了源但没重建包 **" % tag))
        else:
            n_ok += 1
    # B 派生物：docx 必须晚于它的源
    for dx, srcs, nm in (
            (os.path.join(_PKG, "01_manuscript_text", "Manuscript_main_text.docx"), _MAN_SRC, "正文"),
            (os.path.join(_PKG, "02_supplementary_tables", "Supplementary_Tables.docx"), _SUPP_SRC, "补充表"),
            (os.path.join(_UP, "01_正文", "Manuscript_main_text.docx"), _MAN_SRC, "上传包正文"),
            (os.path.join(_UP, "05_补充材料", "Supplementary_Tables.docx"), _SUPP_SRC, "上传包补充表")):
        if not os.path.exists(dx):
            bad.append(("(包 vs 源)", "%s docx" % nm, "缺件", dx))
            continue
        newer = [p for p in srcs if os.path.exists(p) and os.path.getmtime(p) > os.path.getmtime(dx)]
        if newer:
            bad.append(("(包 vs 源)", "%s docx" % nm, "落后于源",
                        "** 源比它新：%s ⇒ 重建 docx **"
                        % "、".join(os.path.basename(p) for p in newer[:4])))
    # C 两包之间
    for a, b, nm in (
            (os.path.join(_PKG, "01_manuscript_text", "Manuscript_main_text.docx"),
             os.path.join(_UP, "01_正文", "Manuscript_main_text.docx"), "正文 docx"),
            (os.path.join(_PKG, "02_supplementary_tables", "Supplementary_Tables.docx"),
             os.path.join(_UP, "05_补充材料", "Supplementary_Tables.docx"), "补充表 docx"),
            (os.path.join(_PKG, "04_figures_supplementary", "Supplementary_Legends.md"),
             os.path.join(_UP, "04_图注", "Supplementary_Legends.md"), "补充图例")):
        if os.path.exists(a) and os.path.exists(b) and _sha(a) != _sha(b):
            bad.append(("(包 vs 源)", nm, "两包不一致",
                        "** 两个包里的 %s 不同 ⇒ 只重建了其中一个（先 p40 再 p56）**" % nm))
    # D 临时/锁文件：只列不判
    for r, _, fs in os.walk(_UP):
        for f in fs:
            if f.startswith("~$") or f in (".DS_Store", "Thumbs.db"):
                notes.append(os.path.relpath(os.path.join(r, f), _UP))
    return bad, notes, n_ok, len(pairs)


# ================================================================ 全称声明专扫
# ★★ 2026-10-13（外部核验要求）：本稿**出错最多的一类**是**全称声明** ——
#   即「断言的范围大于其证据的范围」。本轮已抓获的实例：
#     · "for all four L1 metrics tested"（L1 实际有六个）
#     · "no clipping was applied anywhere"（只在全队列成立；各档 14–28% 越界）
#     · "the sign reflects the native scale of each metric"（15 行里 7 行不成立）
#     · "only eight groups survived at |r| ≥ 0.90"（实测 Spearman 7 组、Pearson 9 组）
#   这一类**不是靠正则判对错**（对错要看上下文），而是**靠把它们全部列出来人工逐条核对**。
#   故本扫描**只列不判**，产出清单供通读时逐条过。
UNIVERSAL = [
    ("all", r"\ball\b"), ("both", r"\bboth\b"), ("neither", r"\bneither\b"),
    ("every", r"\bevery\b"), ("none", r"\bnone\b"), ("no … fell", r"\bno\b[^.]{0,40}\bfell\b"),
    ("anywhere", r"\banywhere\b"), ("entirely", r"\bentirely\b"),
    ("always", r"\balways\b"), ("never", r"\bnever\b"),
    ("completely", r"\bcompletely\b"), ("uniformly", r"\buniformly\b"),
    ("consistently", r"\bconsistent(?:ly)?\b"), ("in no case", r"\bin no case\b"),
    ("invariably", r"\binvariably\b"), ("exclusively", r"\bexclusively\b"),
    ("universal", r"\buniversal(?:ly)?\b"), ("across all", r"\bacross all\b"),
]

# ★★ 「危险形状」的**判据**（必须可核，不能是感觉）：
#    量词（all/both/neither/every/none/only/any of the/each of the）**后紧跟**
#      ① 阿拉伯数字或英文数词（two…ten）—— 核对对象是**计数**；或
#      ② 一个**封闭类名词**（folds / signals / systems / strata / layers / statistics /
#         variants / categories / tiers / classes / levels / metrics）—— 核对对象是**集合的完备性**。
#    本轮的四个真实实例全落在这两种形状里：
#      `all four L1 metrics`（②：L1 实为六个）｜`only eight groups`（①：实为 7/9）
#      `the native scale of each metric`（②：15 行里 7 行不成立）｜`no clipping … anywhere`（量词单独出现）
CLOSED_CLASS = (r"folds?|signals?|systems?|strata|layers?|statistics|metrics?|variants?|"
                r"categories|tiers?|classes|levels?|methods?|estimators?|organisms?|drugs?|"
                r"transfers?|tasks?|groups?|directions?|panels?|conventions?|learners?|"
                r"diagnostics?|coefficients?|correlations?|intervals?|associations?")
NUMWORD = r"\d+|two|three|four|five|six|seven|eight|nine|ten|eleven|twelve|thirteen|fourteen|fifteen"
RISKY_RE = re.compile(
    r"(?i)\b(all|both|neither|every|none|only|any of the|each of the)\s+"
    r"(?:of\s+(?:the\s+)?)?(?:%s)\s+(?:%s)\b" % (NUMWORD, CLOSED_CLASS))

# 结构计数表 —— 核「all 186」这类断言时的**基准**（全部来自冻结产物或设计）
STRUCT_COUNTS = {
    "candidate transfers": 228, "eligible transfers": 186, "excluded candidates": 42,
    "harmful transfers": 57, "stable transfers": 43, "organism–drug groups (eligible)": 37,
    "candidate groups": 38, "within-model pairs": 93, "evaluable folds": 26,
    "single-class folds": 11, "evaluable tasks": 134, "signals": 15,
    "L1 statistics": 6, "L2 signals": 6, "L3 estimators": 3, "health systems": 3,
    "drugs": 20, "sample-size levels": 5, "resampling repeats": 200,
    "organism–drug clusters": 37,
}


def scan_universal(scan_root, max_each=40):
    """只**列**全称声明，不判对错 —— 清单供通读时逐条核对覆盖范围。

    ★★ 2026-10-14 修（**覆盖面缺陷**）：本函数原先只读 `.md` / `.txt`
      （`exts = (".md", ".txt")`）—— **看不见 docx**。而 docx 正是**要投出去的那份正文**。
      此前它之所以「看起来」覆盖了正文，是因为包里还带着 `01_正文/源_markdown/*.md`
      的副本（旧布局）。那些副本一撤，**全称声明专扫就完全不再覆盖正文** ——
      计数从 235/36 掉到 96/12，看数字像是「清理了重复件」，实际是**判据瞎了**。
      ⇒ 改为与 `scan()` 同一套取文机制（`texts_of`，可读 docx / xlsx 内部 XML），
        并对**任何**扩展名都尝试，读不出就跳过。**覆盖面的缺陷，数字上看不出来。**
    """
    out = {}
    for dp, _, fs in os.walk(scan_root):
        for f in fs:
            rel = os.path.relpath(os.path.join(dp, f), scan_root).replace("\\", "/")
            if re.match(r"^(00_作者待填|08_代码归档|07_画图代码_可运行|09_验证证据)", rel) \
                    or rel in ("README_投稿包说明.md", "GENERATION.md") \
                    or "新建 文本文档" in rel:
                continue                       # 非投稿件（含作者转储），不列
            for _lab, _txt in texts_of(os.path.join(dp, f)):
                if not _txt:
                    continue
                _emit(rel, _plain(_txt), out)
    return out


# ★ 2026-10-14：`texts_of` 对 docx/xlsx 返回的是**内部 XML**，不是纯文本。
#   直接扫会把**排版属性**当成正文的全称词 —— 实测 `w:val="none"`（边框/底纹设置）
#   把 `none` 从 3 处抬到 **303 处**，整张清单立刻失去意义。
#   ⇒ 只有在**确认是 Office XML** 时才剥标签：不能对普通文本无条件 `re.sub(r"<[^>]+>")`，
#     那会把正文里的 `predicted ΔAUROC < −0.10, and … > …` 当成一个标签整段吃掉。
_XMLY = re.compile(r"<\?xml|<w:|</w:|<a:|</a:")
def _plain(txt):
    return re.sub(r"<[^>]*>", " ", txt) if _XMLY.search(txt) else txt


def _emit(rel, txt, out):
    """把一段文本里的全称声明逐个记进 `out`（按词分组）。"""
    for name, pat in UNIVERSAL:
        for m in re.finditer(pat, txt, re.I):
            seg = txt[max(0, m.start() - 60):m.end() + 70].replace("\n", " ")
            seg = re.sub(r"\s+", " ", seg).strip()[:150]
            # ★ 「危险形状」优先：量词落在**计数**或**封闭类**上 ——
            #   本轮的四个实例全是这个形状（all four L1 / eight groups /
            #   native scale of each metric / no clipping … anywhere）。
            risky = bool(RISKY_RE.search(seg))
            out.setdefault(name, []).append((rel, seg, risky))
    return out


def main():
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass
    argv = sys.argv[1:]
    negtest = None
    if "--negtest" in argv:
        negtest = argv[argv.index("--negtest") + 1]
    target = negtest or DEFAULT

    # ★ 2026-10-13（外部核验要求）：**扫成品 zip，不扫工作目录**。
    #   理由：p41 的 V7-3 明文禁止脚本名，而 S2 里仍带着 —— 说明此前扫的不是最终交付物。
    tmp = None
    zip_listing = None
    if target.lower().endswith(".zip"):
        import tempfile
        import zipfile as _zf
        tmp = tempfile.mkdtemp(prefix="p60zip_")
        with _zf.ZipFile(target) as z:
            zip_listing = z.namelist()
            z.extractall(tmp)
        print("=" * 100)
        print("p60 — 投稿包泄漏扫描闸门（**扫成品 zip 解压后的内容**）")
        print("  zip：%s（%d 个成员）" % (target, len(zip_listing)))
        target = tmp

    print("=" * 100)
    print("p60 — 投稿包泄漏扫描闸门")
    print("  扫描根：%s" % target)
    print("  模式：%s" % ("**负向测试**（预期判红）" if negtest else "正式（零命中才算过）"))
    print("=" * 100)
    n, hits, author_hits = scan(target)
    # ★ zip 清单检查必须**在打印命中表之前**追加，否则它进了计数却不在表里（自相矛盾的输出）。
    if zip_listing is not None:
        forbidden = sorted({m.split("/")[0] for m in zip_listing
                            if m.startswith(("00_作者待填", "09_验证证据"))
                            or "新建 文本文档" in m})   # ← 三类：作者骨架 / 可选证据 / 你的转储
        if forbidden:
            hits.append(("(zip 清单)", "不应上传的条目", "排除清单", "、".join(forbidden)))
    print("  文件数 %d ｜ 命中 %d 条" % (n, len(hits)))
    if DISCLOSE:
        print("")
        print("  【查过并披露（非命中）】文档属性：")
        for d in DISCLOSE:
            print("      " + d)
    if hits:
        by_kind = {}
        for rel, lab, kind, ev in hits:
            by_kind.setdefault(kind, []).append((rel, lab, ev))
        for kind, rows in sorted(by_kind.items(), key=lambda kv: -len(kv[1])):
            print("\n  【%s】%d 条" % (kind, len(rows)))
            for rel, lab, ev in rows[:6]:
                print("      %-58s %-22s %s" % (rel[:58], lab[:22], ev[:96]))
            if len(rows) > 6:
                print("      …（其余 %d 条同类）" % (len(rows) - 6))
    print()
    if negtest:
        print("  负向测试结论：%s" % ("**判红 —— 闸门有判别力**" if hits else "*** 零命中：闸门对这份样本无判别力，须查 ***"))
        return 0 if hits else 1
    for _dn, _why in (("00_作者待填_★先看这里", "含中文说明与占位符，是给作者的骨架，不是投稿件"),
                      ("09_验证证据_可选", "投不投由你定；含测试编号与脚本名，不是投稿件")):
        if os.path.isdir(os.path.join(target, _dn)):
            print("  ⚠ 需人工处置：`%s/` 在包内 —— **打上传 zip 时应排除**" % _dn)
            print("     （%s）" % _why)
    print()
    # ---- 包 vs 源（2026-10-15 补：此前无任何判据回答"包是不是现在这批源做的"）----
    _pf_bad, _pf_notes, _pf_ok, _pf_n = check_pkg_fresh()
    print()
    print("  【包 vs 源】拷贝件逐位比对 %d 对：相同 %d ｜ **不符 %d**"
          % (_pf_n, _pf_ok, len([x for x in _pf_bad if x[2] == "与源不符"])))
    if _pf_bad:
        print("      ⚠ 下面这些说明**包不是现在这批源做的** —— 重建：先 p40 再 p56")
        for rel, lab, kind, ev in _pf_bad[:12]:
            print("      %-46s %-12s %s" % (lab[:46], kind, ev[:78]))
        hits.extend(_pf_bad)
    else:
        print("      ✔ 包内每一件都与源一致，且两个包互相一致")
    if _pf_notes:
        print("      ⚠ 瞬时/锁文件（Word 生成，**上传前必须不在**）：%s" % "、".join(_pf_notes))
        print("         → 关掉 Word 后重跑 p56 即可清掉")

    # ---- 全称声明专扫（只列不判；清单供通读逐条核对）----
    uni = scan_universal(target)
    all_rows = [(n, r, s, k) for n, v in uni.items() for (r, s, k) in v]
    n_uni = len(all_rows)
    risky = [(n, r, s) for (n, r, s, k) in all_rows if k]
    print()
    print("=" * 100)
    print("【全称声明专扫】只**列**不判 —— 本稿出错最多的一类（范围大于证据）")
    print("  共 %d 处；其中**危险形状**（量词落在计数/封闭类上）%d 处。" % (n_uni, len(risky)))
    print("=" * 100)
    print("\n  ★★ 危险形状优先核（本轮的四个实例全在这一类）：")
    for n, rel, seg in risky[:24]:
        print("      [%s] %-30s %s" % (n, rel[:30], seg[:96]))
    if len(risky) > 24:
        print("      …（其余 %d 处同形）" % (len(risky) - 24))
    print("\n  —— 其余按词列出（大多无害，通读时抽核）——")
    for name, rows in sorted(uni.items(), key=lambda kv: -len(kv[1])):
        print("  ▸ %s（%d 处，危险 %d）" % (name, len(rows), sum(1 for r in rows if r[2])))

    # ★★ 「带数字的全称」的核对基准表 —— 量词的核对对象是**计数**，
    #   故把结构计数一并打印，供逐条比对（这是可核的部分，不是清单）。
    print()
    print("  【核对基准】结构计数（全称若带数字，必须落在这些值上）：")
    for k, v in STRUCT_COUNTS.items():
        print("      %-34s %d" % (k, v))

    # ★★ 反闸门（外部核验要求）：**台账里没有的数字不得出现在稿里**。
    #   本函数只做**候选提取**（正文数字 → 待查表），台账侧的对账在 p15。
    print()
    print("  【数字候选】正文里出现的全部数字 token（供反向对账：每一个都应有台账行或结构声明）")
    # ★ 2026-10-14 修：此处原写 `scan_root` —— 那是 `scan()` / `scan_universal()` 的**形参名**，
    #   在 `main()` 里并不存在（`main()` 里的根目录叫 `target`）。⇒ 该段**从未成功执行过**：
    #   每次都在这里 `NameError` 崩掉，于是**最后两行结论（泄漏闸门判定）永远不打印**，
    #   而上游的命中表已经打完了 ⇒ 看日志的人只看到一张表，**看不到判定**。
    #   这是 INC-021 族的又一例：**报表打了大半，判据没跑到**。修法是改用 `target`。
    nums = {}
    for dp, _, fs in os.walk(target):
        for f in fs:
            if not f.lower().endswith((".md", ".txt")):
                continue
            rel = os.path.relpath(os.path.join(dp, f), target).replace("\\", "/")
            if re.match(r"^(00_作者待填|08_代码归档|07_画图代码_可运行|09_验证证据)", rel) \
                    or rel in ("README_投稿包说明.md", "GENERATION.md") or "新建 文本文档" in rel:
                continue
            try:
                txt = io.open(os.path.join(dp, f), encoding="utf-8", errors="ignore").read()
            except Exception:
                continue
            for m in re.finditer(r"(?<![\w.])(\d+(?:\.\d+)?)", txt):
                nums.setdefault(m.group(1), 0)
                nums[m.group(1)] += 1
    print("      共 %d 个**不同的**数字 token：%s"
          % (len(nums), ", ".join(sorted(nums, key=lambda x: -nums[x])[:30])))

    # ---- 作者待办清单（**不判红**，但**必须为 0 才能投**）----
    print()
    if author_hits:
        print("  【作者待填占位符】%d 处 —— **不是泄漏**，是投稿前待办" % len(author_hits))
        for rel, lab, why, ev in author_hits[:10]:
            print("      %-52s %-22s %s" % (rel[:52], why[:22], ev[:70]))
        if len(author_hits) > 10:
            print("      …（其余 %d 处）" % (len(author_hits) - 10))
    else:
        print("  【作者待填占位符】0 处 —— 作者项已全部落实")

    ok = not hits
    print()
    print("*** %s ***" % ("泄漏闸门通过：零命中" if ok else "泄漏闸门未通过：见上表"))
    print("    ★ 全称声明那 %d 处**不影响**本闸门的通过与否，但必须在通读时逐条核对。" % n_uni)
    print("    ★ 作者待填占位符 %d 处 —— 同样不影响本闸门，但**投稿前必须为 0**。" % len(author_hits))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
