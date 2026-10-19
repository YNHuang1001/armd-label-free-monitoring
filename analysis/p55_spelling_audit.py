# -*- coding: utf-8 -*-
r"""p55 — 拼写 / 术语一致性判据（**全稿**，不只 Discussion）

================================ 为什么存在（INC-021 未了项，A5.8）================================
全稿唯一的拼写判据原先长在 `p19` 里，而 **`p19` 只审 Discussion** —— 它覆盖 1/4 的正文明文。
后果（2026-10-09 实测）：`Results_v3_clean.md` 两处 `programme`、`Methods_v1.md` 三处
`artefact` 与一处 `modelling` **长期在册，从未被发现**；发现它们的是**人眼**，不是判据。

**人工扫描可以是发现手段，不能是判据。**

================================ ★ 覆盖范围（规则 26：覆盖要点名）================================
本判据覆盖 **8 个文件**，逐条列出，未覆盖的显式写明：

    manuscript_renumbered/Introduction_v1.md
    manuscript_renumbered/Methods_v1.md
    manuscript_renumbered/Results_v3_clean.md
    manuscript_renumbered/Discussion_v1.md
    Abstract_v1.md
    Title_v1.md
    figures/final/Figure_Captions.md          ← 主图注
    figures/final/Supplementary_Legends.md    ← 补充图例

**未覆盖（显式声明）**：`Supplementary_Table_S*.md`（三张补充表；其正文由 `p21e` 渲染，
拼写受生成器控制，改生成器即改全部）；`References_FINAL.md`（**逐字引用，不得改** ——
英式拼写在文献题录里是合法的）。
"""
import io
import os
import re
import sys

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "分析结果")

COVERED = [
    os.path.join(ROOT, "manuscript_renumbered", f) for f in
    ("Introduction_v1.md", "Methods_v1.md", "Results_v3_clean.md", "Discussion_v1.md")
] + [os.path.join(ROOT, "Abstract_v1.md"),
     os.path.join(ROOT, "Title_v1.md"),
     os.path.join(ROOT, "figures", "final", "Figure_Captions.md"),
     os.path.join(ROOT, "figures", "final", "Supplementary_Legends.md")]

NOT_COVERED = {
    "Supplementary_Table_S1/S2/S3.md": "由 p21e 渲染，拼写受生成器控制（改生成器即改全部）",
    "References_FINAL.md": "**逐字引用，不得改** —— 英式拼写在文献题录里合法",
}

# 美式拼写为全稿规范。每条 = (英式模式, 说明)。**只收真正的英式变体**，
# 且必须避开**词内子串**误伤（首版把 `organism` 匹配给了 `\borganis`，55 处假阳性）。
UK = {
    "programme": r"\bprogramme\b",
    "artefact": r"\bartefact\b",
    "modelling": r"\bmodelling\b",
    # ★ 2026-10-09 自纠：初版作 `\banalyse[ds]?\b`，把 **`analyses`** 也抓了进来。
    #   而 `analyses` 是 **analysis 的复数名词**，**美式拼写同样是 analyses** ——
    #   英美之别只在**动词**（analyse/analyzed）。首跑因此报 5 处假阳性（全在 Methods/Results）。
    #   ⇒ 只取动词形：analyse / analysed / analysing（**排除 analyses**）。
    "analyse(verb)": r"\banalys(?:e|ed|ing)\b",
    "behaviour": r"\bbehaviour",
    "colour": r"\bcolour",
    "whilst": r"\bwhilst\b",
    "fulfil": r"\bfulfil\b",
    "judgement": r"\bjudgement\b",
    "ageing": r"\bageing\b",
    "catalogue": r"\bcatalogue\b",
    "dialogue": r"\bdialogue\b",
    "centre (prose only)": r"\bcentres?\b",
    "labelling": r"\blabelling\b",
    "cancelled": r"\bcancelled\b",
}
# 允许清单：出现也**不算违规**（逐字引用的文献题录、专有名词）
ALLOW = re.compile(r"multi-centre|Multi-centre")

REP = io.open(os.path.join(OUT, "report_p55_spelling.txt"), "w", encoding="utf-8")


def p(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    REP.write(s + "\n")
    REP.flush()


p("=" * 100)
p("p55 — 拼写 / 术语一致性判据（全稿）")
p("=" * 100)
p(f"覆盖 {len(COVERED)} 个文件：")
for f in COVERED:
    p(f"    {'✔' if os.path.exists(f) else '✘ 缺'}  {os.path.relpath(f, ROOT)}")
p("未覆盖（显式声明）：")
for k, v in NOT_COVERED.items():
    p(f"    —  {k} —— {v}")
p("")

hits = []
for f in COVERED:
    if not os.path.exists(f):
        hits.append((os.path.relpath(f, ROOT), 0, "**文件不存在**", ""))
        continue
    for ln, line in enumerate(io.open(f, encoding="utf-8"), 1):
        for name, pat in UK.items():
            for m in re.finditer(pat, line, re.I):
                seg = line[max(0, m.start() - 50):m.end() + 50].replace("\n", " ")
                if ALLOW.search(seg):
                    continue
                hits.append((os.path.relpath(f, ROOT), ln, name, seg.strip()))

p(f"命中 {len(hits)} 处：")
for f, ln, name, seg in hits:
    p(f"    [WRONG] {f}:{ln}  {name}")
    p(f"            …{seg[:130]}…")
p("")

p("=" * 100)
if hits:
    p(f"*** p55 **未通过** —— {len(hits)} 处英式拼写残留（全稿规范为**美式**）***")
else:
    p("*** p55 通过 —— 覆盖范围内的 8 个文件无英式拼写残留 ***")
REP.close()
sys.exit(1 if hits else 0)
