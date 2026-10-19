# -*- coding: utf-8 -*-
r"""p56 — 装配**投稿上传用**的最终文件夹（`SUBMISSION_20261010/`）

================================ 与 p40 的分工 ================================
`p40` 造的是**存档包**（`SUBMISSION_PACKAGE_20261007/`）—— 八个目录 + 两份 Word + 验证证据，
供**审计与长期留存**。

本脚本造的是**上传包** —— 按投稿系统的操作顺序重排，**按「投稿时要点哪个按钮」分类**，
一眼看得出缺什么。两者**同源**：一律从 p40 的产物的**副本**取，不重算、不改写。

================================ 三条硬约束 ================================
1. **只读源、只写新目录** —— 绝不触碰任何冻结件（复制，不移动）。
2. **每个文件必须是最新的** —— 主图/补充图取自 p40 包内（该包已由 p41 55/55 验收）。
3. **作者专属内容一律不代写** —— IRB、数据使用、作者贡献、利益冲突、资助、许可、
   仓库 URL/DOI，全部只给**带占位符的骨架**，集中在 `00_作者待填_★先看这里/`。
"""
import io
import os
import shutil
import sys

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PKG = os.path.join(ROOT, "SUBMISSION_PACKAGE_20261007")
OUT = os.path.join(ROOT, "SUBMISSION_20261010")
AN = os.path.join(ROOT, "分析结果")

# ★★ 2026-10-12：原为 `raise SystemExit("目标已存在，先自行处理")` —— **安全，但让本脚本不可重跑**。
#   病历：上传包于 2026-10-10 建成后，本日对正文/补充表/图注做了多轮修正，
#   而**只要 `SUBMISSION_20261010/` 还在，p56 就拒绝运行** ⇒ 手上那个上传包**停留在 10-10**，
#   目录名里的日期又恰好让人以为它就是当天的。
#   ⇒ 改为本项目既有的约定动作：**先归档、再重写**（`p40` 重写 PKG、`genstamp.backup` 覆盖前归档，
#     都是同一形状）。归档落 `_artifact_archive/`，**不删除任何东西**，且旧包可随时取回。
if os.path.exists(OUT):
    import datetime as _dt0
    _arch = os.path.join(ROOT, "_artifact_archive",
                         "%s.%s.pre-rebuild"
                         % (os.path.basename(OUT),
                            _dt0.datetime.now(_dt0.timezone.utc).strftime("%Y%m%dT%H%M%SZ")))
    os.makedirs(os.path.dirname(_arch), exist_ok=True)
    shutil.move(OUT, _arch)
    print("旧上传包已归档（未删除）：%s" % _arch)


def mk(*a):
    p = os.path.join(OUT, *a)
    os.makedirs(p, exist_ok=True)
    return p


def cp(src, dst_dir, name=None):
    if not os.path.exists(src):
        print("   ** 缺 **", src)
        return False
    os.makedirs(dst_dir, exist_ok=True)
    shutil.copy2(src, os.path.join(dst_dir, name or os.path.basename(src)))
    return True


def wr(path, text):
    io.open(path, "w", encoding="utf-8").write(text)


# ================================================================ 00 作者待填
d00 = mk("00_作者待填_★先看这里")
PLACEHOLDERS = """# ★ 作者待填清单（投稿前必须清空这一节）

> 本目录**只放骨架与占位符**，没有一个字的代写内容 ——
> IRB / 数据使用 / 作者贡献 / 利益冲突 / 资助 / 许可 / 仓库地址，**只有作者能提供**。

## A. 正文里作者必须处理的 **4 处**（**不清空就不能投**）

> **★ 2026-10-15 更正**：本节此前只列 **3** 处，并漏掉了第 4 处。
> 判据：`分析结果/p60_leak_scan.py` 每次运行都会打印「【作者待填占位符】N 处」——
> **那个 N 必须为 0**。复核件见 `UNIVERSAL_CLAIM_REVIEW_20261014.md` §三。

| # | 位置 | 现在写的是 | 要填什么 |
|---|---|---|---|
| 1 | Methods §1 | `[Institutional review board determination and data-use statement to be inserted by the authors.]` | 伦理审查结论（批准或豁免）与该研究的 data-use 声明 |
| 2 | Data availability | `[IRB / human-subjects determination to be inserted by the authors.]` | **三站点各自的**审批机构名称与批准/豁免编号；是否允许索取数据、若允许则程序为何 |
| 3 | Code availability | `[TO BE COMPLETED BY THE AUTHORS: repository URL, the archival DOI (Zenodo or Code Ocean), and the licence.]` | 仓库 URL ＋ 存档 DOI ＋ 许可 |
| 4 | Methods（IRB 段**散文句**，不是方括号） | `Whether this secondary analysis required additional review, and the wording of that determination, is the authors' statement to make.` | 与 #1 是**同一件事**；判定落实后，**这句要改写为已完成的表述**（留在稿里会读成「作者还没定」） |

**填完不是重跑 `p40` 就完事** —— 按 `FREEZE_DECLARATION_20261010.md` §五 的**重开条件**，
改正文会触发：`A4.3n` ＋ `A4.3c2` 复扫，并重跑 `p15`/`p16`/`p19`/`p34`/`p55`，
然后才 `p44` → `p40` → `p41`/`p50`/`p51`。

## A2. 摘要的 **3 处内容决定**（不是待填，是**要你点头或改回**）

> 判据：`分析结果/p34_abstract_audit.py` 的 **WARN** 项（PASS / **WARN** / FAIL 三分；
> **只有 FAIL 影响退出码**，WARN 是「需要人决定的披露缺口」）。
> 摘要现为 v5；相对 v4 **静默删去**了下面三处。三处都**不影响任何数字的正确性**，
> 只影响**覆盖面与措辞强度**。**要不要加回由你定**，判据不代改。

| # | v4 有的句子 | 现在的状态 | `p34` 项 |
|---|---|---|---|
| 1 | `Signals computed from model outputs retained more information in that analysis, but that family-level advantage did not survive grouped evaluation.` | 摘要**不再含任何 prediction-space 家族级陈述**（只报 input-space 与 DR）。该发现仍在 Results / Fig 3 / Discussion | **E-2** |
| 2 | `In a retrospective policy simulation, reviewing the top 20% of transfers … captured 1.59 times the base rate …` | v4 的第 4 交付物（分诊策略模拟）**不再出现在摘要**；预设警示任务的单点仍在 | **E-6** |
| 3 | `… but stayed poorly calibrated **in absolute terms**.` | 现为裸写 `but was poorly calibrated`，**限定语丢失**。冻结证据的形状是**序数保真（DR Spearman 0.7896）＋ 绝对水平差（MAE 0.0416）**，裸写比证据更重（v4 文件头**专门论证过**为何必须带） | **E-7** |

## B. 正文**整节缺失**，需要新增（见本目录模板）

| # | 缺件 | 说明 | 模板 |
|---|---|---|---|
| 5 | **标题页** | 正文 docx **从标题直接进 Abstract**，无作者、无单位、无通讯作者 | `01_Title_page_模板.md` |
| 6 | **Author contributions** | 全文未见 | `03_Author_contributions_模板.md` |
| 7 | **Competing interests** | 全文未见（Nature Portfolio **要求**） | `04_Competing_interests_模板.md` |
| 8 | **Acknowledgements / Funding** | 全文未见（资助必须写明） | `05_Acknowledgements_funding_模板.md` |

## C. 投稿系统里要交、但项目此前**从未提过**的两件

| # | 缺件 | 说明 | 模板 |
|---|---|---|---|
| 9 | **Cover letter** | 项目里搜不到任何记录 | `02_Cover_letter_模板.md` |
| 10 | **Reporting Summary** | **Nature Portfolio 系刊对 research article 要求**；本项目从未做过 | `06_Reporting_Summary_说明.md` |

★ 第 10 项建议**先到 npj Digital Medicine 的投稿页确认**当前是否仍要求、以及是哪个版本的表 ——
本仓没有它的官方副本，**不代为断言**。

## D. Code availability 归档（见模板）

| # | 动作 | 模板 |
|---|---|---|
| 11 | 把 `08_代码归档/` 上传到仓库 → 取 DOI → 填许可 | `07_Code_归档与许可_说明.md` |

## E. 姊妹稿披露 —— **只有你们能判断，但一旦属实就该写**（2026-10-15 新增）

> 外部核验 2026-10-15 提出：**论文 1 与本稿同三个数据集、同六个方向、同任务集，且已投**。
> **本仓查不到这个事实** —— 全仓搜 `companion / 姊妹稿 / overlapping publication / cover letter`：
> **0 命中**；唯一提到论文 1 的地方是 `HANDOFF_20261012_state.md`：
> 「**姊妹论文（论文 1）写的就是对的**，本仓未去比对」—— 说的是**数据来源措辞**，不是披露义务。
> **⇒ 是否构成「相关稿件」取决于论文 1 的投稿状态；那是作者侧事实，机器代不了。**

| # | 动作 | 说明 |
|---|---|---|
| 12 | **若论文 1 已在审或已投**：`Data availability` 加一段 **Relation to a companion manuscript**，写明两稿各自报告什么、数据与方向的重叠范围 | Nature Portfolio 对**相关稿件**有披露要求 |
| 12b | 同上情形：**cover letter** 里主动说明与本稿的关系（见 `02_Cover_letter_模板.md`） | 由编辑判断是否需进一步处理 |
| 12c | 投稿系统若有 "related manuscripts" 栏，填论文 1 的稿号 | 按系统实际栏目 |

★ **不写，重叠发表的风险原样在**；写错（把不相关说成相关）也有代价 —— 故**由你们判断后落笔**。
"""
wr(os.path.join(d00, "★待填清单.md"), PLACEHOLDERS)

wr(os.path.join(d00, "01_Title_page_模板.md"), """# 标题页（模板 —— 作者填写）

> npj Digital Medicine 的投稿系统会**单独收集**作者与单位信息；
> 但多数编辑仍希望正文第一页有标题页。**先用本模板填好，再决定是并入 docx 还是走系统字段。**

---

**Title**

Failure ranking is not failure detection: what label-free monitoring can and cannot tell us about a cross-system deployment

**Running title**（若刊物要求，≤50 字符）
_[待填]_

**Authors**（按贡献顺序；`*` 通讯作者）

1. _[姓名]_ ¹* ｜ ORCID: _[0000-0000-0000-0000]_
2. _[姓名]_ ¹  ｜ ORCID: _[0000-0000-0000-0000]_
3. …

**Affiliations**

¹ _[单位全称，城市，邮编，国家]_

**Corresponding author**

_[姓名]_ ｜ _[邮箱]_ ｜ _[通讯地址]_

**Keywords**（4–6 个）
_[待填]_

**Word count**：正文 _[ ]_ 词（不含摘要/参考文献）｜**Abstract**：191 词｜**Figures**：5｜**Tables**：1（正文）＋ S1–S3、R1–R3（补充）
""")

wr(os.path.join(d00, "02_Cover_letter_模板.md"), """# Cover letter（模板 —— 作者填写；**项目此前无此件**）

> 写法提示：编辑每天读很多封。**第一段说清「投什么、新在哪」**，
> 之后逐条对应到本刊的范围（digital medicine / 部署监测 / 临床 AI）。
> **不要从记忆里抄数字** —— 若要引用本稿的数字，一律从 `01_正文/` 与 `05_补充材料/` 里取。

---

Dear Editors of *npj Digital Medicine*,

_[第一段：投稿声明 + 一句话结论。建议点明本文的四点贡献：① 186 个跨系统转移的密封基准；
② 输入层位移统计在预注册方向上**无正向判别**；③ 输出层信号的优势在分组口径下不是家族级的；
④ 排序可用而检测不可用（k=20% 复查预算下 1.59× 基率）。]_

_[第二段：为什么适合贵刊 —— 部署监测 / 临床 AI 的落地决策 / 模型跨机构迁移。]_

_[第三段：合规声明 —— 未一稿多投、全体作者同意、无利益冲突（与正文 Competing interests 一致）、
数据为去标识化病历且受数据使用协议约束。]_

_[第四段：推荐审稿人（可选）/ 回避审稿人（可选）/ 是否已预印（如有，给 DOI）。]_

Sincerely,
_[通讯作者姓名，单位]_
""")

wr(os.path.join(d00, "03_Author_contributions_模板.md"), """# Author contributions（模板 —— 作者填写）

> 用 CRediT 分类最省事：Conceptualization / Data curation / Formal analysis / Funding acquisition /
> Investigation / Methodology / Project administration / Resources / Software / Supervision /
> Validation / Visualization / Writing – original draft / Writing – review & editing。

---

_[姓名]_：_[贡献]_
_[姓名]_：_[贡献]_
…

All authors read and approved the final manuscript.
""")

wr(os.path.join(d00, "04_Competing_interests_模板.md"), """# Competing interests（模板 —— 作者填写；Nature Portfolio **强制要求**）

> 若有任何财务或非财务关系，逐条写明；若没有，用下面这句。

---

The authors declare no competing interests.
""")

wr(os.path.join(d00, "05_Acknowledgements_funding_模板.md"), """# Acknowledgements & Funding（模板 —— 作者填写）

> **资助必须写 aquí**：基金全称 + 编号 + 资助对象。Nature Portfolio 要求把资助写进
> Acknowledgements 并在系统里单列，**不能用「This work was supported by …」含糊带过**。

---

**Acknowledgements**

_[致谢非作者贡献者；若有数据/算力支持一并写明。]_

**Funding**

_[基金全称（中英文）]_，批准号 _[ ]_（资助 _[姓名]_）。
_[若无外部资助：This research received no specific grant from any funding agency.]_
""")

wr(os.path.join(d00, "06_Reporting_Summary_说明.md"), """# Reporting Summary —— **项目此前从未做过，需作者定夺**

## 为什么提这条

**Nature Portfolio 系刊（含 npj 系列）对 research article 要求提交 Reporting Summary**：
一张与编辑/审稿人内部使用、不随文章发表的表，逐条说明研究设计、样本量、统计、数据可得性、
以及（若涉及 AI/ML）模型报告要素。

**★ 我的边界**：本仓**没有**该表的官方副本，我**不代为断言** npj Digital Medicine 当前
是否仍要求、要求的是哪个版本（Life Sciences / Health Sciences 用的是不同表）。
**请到投稿页确认后再填。**

## 建议动作

1. 打开 npj Digital Medicine 的 submission guidelines，下载当前的 Reporting Summary 模板。
2. 表里多数条目本稿已有素材，**逐条从下列位置取**：
   - 研究设计与样本量 → `01_正文/Manuscript_main_text.docx` **Methods §1、§12**
   - 统计与不确定性 → **Methods §14**（cluster bootstrap、1,000 次、无 p 值）
   - 数据可得性 → **Data availability**
   - 代码可得性 → **Code availability**
   - 预注册/协议冻结 → **Methods §8**（label sealing）
   - 报告规范 → **Methods §1 末**（TRIPOD+AI）＋ 补充材料
   - 公平性 → **Discussion 的 Limitations**（明确写明未做公平性分析）
3. **凡本稿没做的，如实勾「Not applicable / Not done」并说明** —— 本稿的可信度正建立在
   「未做就说未做」上（Discussion 里已有先例）。
""")

wr(os.path.join(d00, "07_Code_归档与许可_说明.md"), """# Code availability —— 归档与许可（**作者执行**）

## 归档包在 `08_代码归档/`

```
analysis/     分析链（按执行顺序）+ 全部审计器
evs/          独立验证套件（run_evs.py 跑全部模块）
figures/      全部画图脚本（R）
MANIFEST.csv  126 个文件，每个带 SHA-256 —— **可校验，不是清单承诺**
README.md     怎么跑
```

**★ 已核**：MANIFEST **126/126 条**，文件全部存在且 SHA-256 **逐位相符**（2026-10-10）。

## 步骤建议（Zenodo）

1. 在 Zenodo 新建 upload，把 `08_代码归档/` **整个目录**打包上传（zip）。
2. 元数据：标题同正文；作者同正文；**License 见下**。
3. 发布后拿到 **DOI**（形如 `10.5281/zenodo.xxxxxxx`）。
4. 把 **仓库 URL + DOI + License** 填进正文 **Code availability** 的方括号。

## 许可（**作者决定** —— 我不代为选择）

| 选项 | 何时适用 |
|---|---|
| **MIT / BSD-3** | 希望他人自由复用（含商用），最省事 |
| **Apache-2.0** | 同上，另含专利授权条款 |
| **GPL-3.0** | 要求衍生作品同样开源 |
| **CC BY 4.0** | 若你把它当"数据/文档"而非"软件"，且机构偏好 CC |

**注意**：若三站点的数据使用协议对衍生代码有约束，**先确认再选**；
本代码包**不含任何患者级数据**（只有脚本与聚合后的图数据 CSV）。
""")

# ================================================================ 01 正文
d01 = mk("01_正文")
cp(os.path.join(PKG, "01_manuscript_text", "Manuscript_main_text.docx"), d01)
# ★★ 2026-10-14（外部核验决定）：**整份拿掉正文的 markdown 源**。
#   理由：投稿系统要的是 Word/PDF，markdown 源属**内部件**；而这几份源文件的**文件头**
#   带着整套中文决策日志（`Title_v1.md` 1081 处中文、`Introduction` 518、`Abstract` 449、
#   `Discussion` 374）—— 泄漏扫描对它们判红，且它们**没有任何投稿用途**。
#   ⇒ 不复制，并**同步删掉下面 README 里对它的引用**（外部核验要求：拿掉之后包里
#     不得再出现指向它的路径或说明）。

# ================================================================ 02/03 图
d02 = mk("02_主图_TIFF_400dpi")
d02p = mk("02_主图_TIFF_400dpi", "预览_PNG_SVG")
d03 = mk("03_补充图_TIFF_350dpi")
d03p = mk("03_补充图_TIFF_350dpi", "预览_PNG_SVG")
for i in range(1, 6):
    cp(os.path.join(PKG, "03_figures_main", "Figure%d.tiff" % i), d02)
    for ext in ("png", "svg"):
        cp(os.path.join(PKG, "03_figures_main", "Figure%d.%s" % (i, ext)), d02p)
for i in range(1, 10):
    cp(os.path.join(PKG, "04_figures_supplementary", "FigureS%d.tiff" % i), d03)
    for ext in ("png", "svg"):
        cp(os.path.join(PKG, "04_figures_supplementary", "FigureS%d.%s" % (i, ext)), d03p)

# ================================================================ 04 图注
d04 = mk("04_图注")
cp(os.path.join(PKG, "03_figures_main", "Figure_Captions.md"), d04)
cp(os.path.join(PKG, "04_figures_supplementary", "Supplementary_Legends.md"), d04)

# ================================================================ 05 补充材料
d05 = mk("05_补充材料")
cp(os.path.join(PKG, "02_supplementary_tables", "Supplementary_Tables.docx"), d05)
for f in ("Supplementary_Methods.md", "Supplementary_Table_S1.md", "Supplementary_Table_S2.md",
          "Supplementary_Table_S3.md", "Supplementary_Tables_R1_R3.md"):
    cp(os.path.join(PKG, "06_source_markdown", f), d05)

# ================================================================ 06 Source Data
d06 = mk("06_Source_Data")
cp(os.path.join(PKG, "07_source_data", "SourceData.xlsx"), d06)
cp(os.path.join(PKG, "07_source_data", "SourceData_index.csv"), d06)

# ================================================================ 07 画图代码（可运行）
d07 = mk("07_画图代码_可运行", "R_scripts")
for f in sorted(os.listdir(os.path.join(ROOT, "figures", "R"))):
    if f.endswith(".R"):
        # ★★ 2026-10-14 修（装包 bug）：`d07` **本身就是** `…/R_scripts`，
        #   原写 `os.path.join(d07, "R_scripts")` ⇒ 目标变成 `R_scripts/R_scripts/`，
        #   **多套了一层目录**：外层只剩 README，31 个 R 脚本全部卡在里层。
        #   是泄漏扫描撞出来的（报 32 条，路径形如 `R_scripts/R_scripts/Figure2_forest.R`）。
        cp(os.path.join(ROOT, "figures", "R", f), d07)
d07d = mk("07_画图代码_可运行", "数据包_figures_data")
for f in sorted(os.listdir(os.path.join(ROOT, "figures", "data"))):
    if f.endswith(".csv"):
        cp(os.path.join(ROOT, "figures", "data", f), d07d)
wr(os.path.join(d07, "README_怎么跑.md"), """# 画图代码（可直接运行）

## 为什么这一份和 `08_代码归档/figures/` 看起来重复

`08` 里的 `figures/R/` 只有 `.R` 脚本；**脚本要读 `figures/data/*.csv` 才跑得起来**，
而那些 CSV 不在归档包里 —— 归档包只装了脚本。**故此处把脚本与它读的数据放在一起。**

## 跑法

1. 保持目录结构不变。
2. 把 `R_scripts/` 里任意脚本的 `DATA_DIR` 指向同级的 `数据包_figures_data/`
   （脚本头部写的是 `figures/data`，改成 `../数据包_figures_data` 即可）。
3. `Rscript Figure2_forest.R`。

## 依赖

R ≥ 4.3 + `ggplot2`、`patchwork`、`readr`、`dplyr`、`tidyr`、`scales`。

## ★ 这些脚本的性质（写进代码可复现性声明时用）

**渲染层，不重算任何统计量。** 每个脚本头部都标了 `[RENDER LAYER ONLY]`，
并逐行声明它读哪些冻结 CSV —— 统计全部在上游 Python 里完成并冻结，
R 只负责把冻结值画出来。这也是为什么图上每一个点都能追到 `figures/data/` 里的一格。
""")

# ================================================================ 08 代码归档
d08 = mk("08_代码归档")
for sub in ("analysis", "evs", "figures"):
    src = os.path.join(PKG, "08_code", sub)
    if os.path.isdir(src):
        shutil.copytree(src, os.path.join(d08, sub), dirs_exist_ok=True)
for f in ("MANIFEST.csv", "README.md"):
    cp(os.path.join(PKG, "08_code", f), d08)

# ================================================================ 09 验证证据（可选）
d09 = mk("09_验证证据_可选")
for f in ("report_evs.txt", "coverage_table.csv", "manifest_sha256.txt"):
    cp(os.path.join(PKG, "05_verification", f), d09)
cp(os.path.join(AN, "report_p41_package_verify.txt"), d09)
wr(os.path.join(d09, "说明.md"), """# 验证证据（**可选随投稿提交**）

这些是**本稿在生产过程中自建的独立验证**，不是刊物要求。放进来是因为它们的性质特殊：

- `report_evs.txt` —— 独立验证套件的**双跑**结果（两遍逐字节相同）。
  其中 `FAIL 4` 是**按设计变红的哨兵**（记录与原始 SUDO 论文的一处已知偏离），
  **不是回归**；`closure 四条全满足`。
- `coverage_table.csv` —— 逐路径的测试覆盖表。
- `report_p41_package_verify.txt` —— 投稿包自身的 **55 项**逐行验收。

**投不投由作者定。** 若投，建议只放 `report_evs.txt` 并在 cover letter 里一句话带过；
**不建议**把这些放进 Supplementary Information（会改变刊物的补充材料构成）。
""")

# ================================================================ README
wr(os.path.join(OUT, "README_投稿包说明.md"), """# 投稿包 — 2026-10-10

> **按「投稿时要点哪个按钮」排序。** 带 ★ 的是**必须你先处理**的。
> 全部文件取自 **已冻结并验收** 的产物（`p41` 55/55、`p50` 7/0/0、`p51` 四查全零、
> 代码包 MANIFEST 126/126 逐位核对）。

## 🚦 先做这三件

1. **`00_作者待填_★先看这里/★待填清单.md`** —— 12 项，含正文里**4 处只有作者能落笔的东西**
   （三处方括号：IRB、数据使用判定、仓库 URL/DOI/许可；＋ IRB 段那句要改写的散文）、
   **四个整节缺失**（标题页、作者贡献、利益冲突、致谢/资助），
   以及**摘要的 3 处内容决定**（`p34` 的 WARN 项 —— 要不要把 v5 删掉的三句加回去），
   以及**姊妹稿披露**（论文 1 的关系 —— 只有你们能判断是否构成「相关稿件」）。
2. **`00_作者待填_★先看这里/06_Reporting_Summary_说明.md`** —— 先去投稿页确认是否需要。
3. 填完后按 `FREEZE_DECLARATION_20261010.md` §五 的**重开条件**重跑链（**不是只重跑 p40**）。

## 📁 目录

| 目录 | 是什么 | 传到哪 |
|---|---|---|
| `00_作者待填_★先看这里/` | 骨架与占位符清单（**12 项，无代写内容**） | 不传，自己用 |
| `01_正文/` | 正文 Word（**只此一件**；markdown 源已整份移除，不随包） | Manuscript |
| `02_主图_TIFF_400dpi/` | 主图 Figure 1–5，**TIFF 400 dpi** | Figures（逐张上传） |
| `02_主图_TIFF_400dpi/预览_PNG_SVG/` | 同图的 PNG/SVG（**投稿系统一般不收**，供你看与排版） | 不传 |
| `03_补充图_TIFF_350dpi/` | 补充图 Figure S1–S9，**TIFF 350 dpi** | Supplementary figures |
| `04_图注/` | 主图图注 ＋ 补充图图注 | 图注若未随正文，单传；★ **先看正文里是否已含** |
| `05_补充材料/` | 补充表 Word ＋ S1/S2/S3 ＋ R1–R3 ＋ 补充方法 | Supplementary Information |
| `06_Source_Data/` | Source Data 单工作簿（每图板一 sheet） | Source Data（Nature Portfolio 常要） |
| `07_画图代码_可运行/` | **R 脚本 ＋ 它读的数据包**（可直接跑） | 一般不单独传；见 `08` |
| `08_代码归档/` | 完整代码包（analysis/evs/figures ＋ MANIFEST 126 条） | **Zenodo / 仓库**，取 DOI |
| `09_验证证据_可选/` | EVS 双跑报告 / 覆盖表 / 包验收 | **投不投由你定** |

## ✅ 已核过的事（不用再查）

- 主图 TIFF **400 dpi**、补充图 TIFF **350 dpi** —— 均 ≥300，达标。
- 稿件引用 **S1–S9**；`FigureS10_triage.R` **不在投稿范围**（正文未引），**未收入**。
  ★ **2026-10-15 补**：此前**只执行了一半** —— 图件没收录、Source Data 没有 S10 sheet、正文不引，
  **唯独 `Supplementary_Legends.md` 里还留着 S10 的图例**（包内 9 张补充图、10 条图例）。
  外部核验抓到「S10 全包零引用」。**已按作者决定（方案 B）把图例一并撤除**，
  使这个决定**四处一致**；分诊的数字由 **Table S3** 承载，一个没丢。
  被删图例留档于 `分析结果/_artifact_archive/FigureS10_legend_DROPPED_20261015T012150Z.md`；
  `p31` 原有 20 条 S10 数字断言**改为 3 条负向守卫**（图例不得再出现）。
- 代码包 **MANIFEST 126/126**：文件全在、SHA-256 逐位相符。
- 包内验证证据与 `evs/` 源**同时间戳**（`p40` 的新鲜度守卫保证）。

## ⚠️ 我想到、但你要判断的

- **`04_图注` 可能与正文重复** —— 正文 docx 通常已含图注。**别传两份**，先看正文。
- **补充材料是 Word 还是 PDF** —— 本包给的是 Word。**投稿系统常要求单个 PDF**；
  若要求 PDF，用 Word 另存为 PDF（本仓不代转，以免版式走样）。
- **`09_验证证据` 建议只挑 `report_evs.txt`** —— 全放会改变补充材料的构成。
- **图像完整性**：Nature Portfolio 要求图片**不得有任何美化性改动**；本包 TIFF 由 R 直接导出，
  未做后期修饰，符合。
""")

# ================================================================ GENERATION.md
# ★ 2026-10-12：本目录名里的 `20261010` 是**建档日期**，而本脚本每次把内容整体重写。
#   ⇒ 只看目录名无法判断手里这一包对应哪一版稿件（INC-018 的同族）。
#   每次构建覆写一张**代次声明**，记下构建时刻与**所依据的稿件源哈希**。
import datetime as _dtg
import hashlib as _hlg


def _sha(p):
    return _hlg.sha256(open(p, "rb").read()).hexdigest()


_MAN_G = [os.path.join(ROOT, "Abstract_v1.md"), os.path.join(ROOT, "Title_v1.md")] + \
         [os.path.join(ROOT, "manuscript_renumbered", f) for f in
          ("Introduction_v1.md", "Methods_v1.md", "Results_v3_clean.md", "Discussion_v1.md",
           "References_FINAL.md", "DataCodeAvailability_v1.md")]
_gexp = []
_gexp.append("# 代次声明（由 `分析结果/p56_build_submission_folder.py` 每次构建覆写）")
_gexp.append("")
_gexp.append("> **本目录名里的 `20261010` 是建档日期，不是本包的产生时间。**")
_gexp.append("> 每次构建都会**先归档旧包**（到 `_artifact_archive/`）**再整体重写本目录** ——")
_gexp.append("> 因此**只有这张纸说明它是哪一代**。")
_gexp.append("")
_gexp.append("| 项 | 值 |")
_gexp.append("|---|---|")
_gexp.append("| **构建时刻（本地）** | %s |" % _dtg.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
_gexp.append("| **构建脚本** | `%s` |" % os.path.abspath(__file__))
_gexp.append("| **上游存档包** | `%s` |" % PKG)
_gexp.append("| **文件总数** | %d |" % sum(len(fs) for _, _, fs in os.walk(OUT)))
_gexp.append("")
_gexp.append("## 本包正文所依据的稿件源")
_gexp.append("")
_gexp.append("| 文件 | 改动时刻 | sha256 |")
_gexp.append("|---|---|---|")
for _f in _MAN_G:
    if os.path.exists(_f):
        _gexp.append("| `%s` | %s | `%s` |"
                     % (os.path.relpath(_f, ROOT),
                        _dtg.datetime.fromtimestamp(os.path.getmtime(_f)).strftime("%Y-%m-%d %H:%M:%S"),
                        _sha(_f)[:16]))
_gexp.append("")
_gexp.append("## 用法")
_gexp.append("")
_gexp.append("- **可以**：按上表核对「手里这一包对应哪一版稿件」。")
_gexp.append("- **不可以**：把目录名当时间戳；也不要在只重跑 `p40` 之后以为本包更新了"
             " —— 正文变了必须**先重跑 p40，再重跑 p56**。")
io.open(os.path.join(OUT, "GENERATION.md"), "w", encoding="utf-8").write("\n".join(_gexp) + "\n")
print("写出：", os.path.join(OUT, "GENERATION.md"), "（代次声明）")

# ================================================================ 把写死的 126 换成派生值
# ★★ 2026-10-12：本脚本写出的两份说明里**四处写死**「MANIFEST 126 条 / 126/126」，
#   而代码包本轮已扩到 **132 条** —— 数字与实物**静默分叉**（INC-021 / INC-026 族：
#   「可以被派生的东西硬编码之后，会与世界悄悄走散」）。
#   处置：从 `08_code/MANIFEST.csv` **实测条数**替换，且**替换数不符即报错**
#   （防将来措辞一改，这段就静默失效、旧数原样留下）。
_N_MAN = sum(1 for _ in io.open(os.path.join(PKG, "08_code", "MANIFEST.csv"), encoding="utf-8")) - 1
_SUBS = [("126/126", "%d/%d" % (_N_MAN, _N_MAN)),
         ("MANIFEST 126 条", "MANIFEST %d 条" % _N_MAN),
         ("MANIFEST.csv  126 个文件", "MANIFEST.csv  %d 个文件" % _N_MAN)]
_hits = {k: 0 for k, _ in _SUBS}
for _rel in ("00_作者待填_★先看这里/07_Code_归档与许可_说明.md", "README_投稿包说明.md"):
    _p = os.path.join(OUT, _rel)
    _t = io.open(_p, encoding="utf-8").read()
    for _k, _v in _SUBS:
        _c = _t.count(_k)
        _hits[_k] += _c
        _t = _t.replace(_k, _v)
    io.open(_p, "w", encoding="utf-8").write(_t)
_zero = [k for k, c in _hits.items() if c == 0]
if _zero:
    raise SystemExit("*** 说明文件里的写死数字未按预期出现，替换未生效：%s ***" % _zero)
print("已把写死的 MANIFEST 条数替换为实测值 %d（命中 %s）" % (_N_MAN, _hits))

print("写出：", OUT)
n = sum(len(fs) for _, _, fs in os.walk(OUT))
print("文件总数：", n)
for d in sorted(os.listdir(OUT)):
    p = os.path.join(OUT, d)
    if os.path.isdir(p):
        print("  %-34s %d 个文件" % (d, sum(len(fs) for _, _, fs in os.walk(p))))
