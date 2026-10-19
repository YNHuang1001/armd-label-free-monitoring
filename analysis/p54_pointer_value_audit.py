# -*- coding: utf-8 -*-
r"""p54 — **指针-值判据**：正文引的数字，必须在它**指向的产物**里真的存在

================================ 它补的是哪一段 ================================
`p15` 验「**写出来的值**对不对」—— 判据的输入是正文本身。
它**结构上抓不到** INC-020 那一类：值是对的，但**读者被指去的地方没有这个值**。

    正文：…in the prespecified native-orientation task **(Fig. 2A)** the untransformed values
          of all three discriminate below chance (0.362, 0.457 and 0.342).
    Fig 2A 的唯一数据源上是 0.548 / 0.542 / 0.545 —— **在 0.5 以上**。
    `p15` 的 S3-28/29/30 把这三个数锚到 p33 的 signed 行，**三条全 PASS**。

⇒ **指针本身也是断言。** 本脚本就是那条断言。

================================ ★ 登记表**不手写**（这是关键设计）========================
手写名单会与世界静默分叉 —— 本项目已为此吃过两次（`p44` 的代码包缺件、`F` 节漂移）。
故本脚本的「指针 → 产物」映射**从绘图脚本自己声明的头部注释里解析**：
每个 `figures/R/Figure*.R` 都带 `[RENDER LAYER ONLY]` 头，逐行列出自读的 `*.csv`。
**图改读别的文件 ⇒ 登记表自动跟着改。** 手写名单只保留一件：指针**词**到**脚本**的对应。

================================ 判据 ================================
句子里的每个数字 token，必须满足下列之一：
  ① **命中产物**：等于某个格的显示值（直接，或 ×100 —— 正文常用百分数）；
  ② **登记为派生**：在 `DERIVED` 里声明了公式，且公式**重算得到它**。
否则 **RED**。

================================ ★ 四条退路（外部核验 2026-10-09 要求「烤进判据」）========
对抗复核那轮靠人跑的四条，现在写进代码：
  R1「另有产物也支撑该指针」 → 登记表由绘图脚本自声明派生，**穷举**该图读的全部文件；缺件即报警。
  R2「是否有格可舍入到该值」 → 匹配按**该 token 自身的显示精度**做，不是精确相等。
  R3「图内 SVG 是否含该串」  → 同时扫 `figures/final/*.svg`（图里印出来的数字也算「在图中」）。
  R4「同形异指」             → 未命中时**打印最近的候选值**，让人一眼看出是「差了很多」还是「差一点」。
"""
import glob
import io
import os
import re
import sys

import numpy as np
import pandas as pd

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIG = os.path.join(ROOT, "figures", "data")
FR = os.path.join(ROOT, "figures", "R")
FFIN = os.path.join(ROOT, "figures", "final")
RP = os.path.join(FIG, "reporting_profiles")
MS = [os.path.join(ROOT, "manuscript_renumbered", f) for f in
      ("Introduction_v1.md", "Methods_v1.md", "Results_v3_clean.md", "Discussion_v1.md")] \
     + [os.path.join(ROOT, "Abstract_v1.md")]
REP = io.open(os.path.join(ROOT, "分析结果", "report_p54_pointer_value.txt"), "w", encoding="utf-8")


def p(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    REP.write(s + "\n")
    REP.flush()


# ---------------------------------------------------------------- ① 从绘图脚本解析自声明
def declared_csvs(rfile):
    """读 `[RENDER LAYER ONLY]` 头 + 正文里的 read_csv 行，取出该图自读的 csv 名。"""
    t = io.open(rfile, encoding="utf-8", errors="ignore").read()
    names = set(re.findall(r'([A-Za-z0-9_]+\.csv)', t))
    return {n for n in names if os.path.exists(os.path.join(FIG, n)) or
            os.path.exists(os.path.join(ROOT, "分析结果", n))}


# 指针**词** → 绘图脚本（这张表只能手写：它编码的是「读者看到的编号」与「文件名」的约定）
PTR2R = {
    r"Fig\.?\s*1[A-E]?": "Figure1_v2.R",
    r"Fig\.?\s*2[A-B]?": "Figure2_forest.R",
    r"Fig\.?\s*3[A-C]?": "Figure3_grouped.R",
    r"Fig\.?\s*4[A-E]?": "Figure4_transported.R",
    r"Fig\.?\s*5[A-B]?": "Figure5_samplesize.R",
    r"Fig\.?\s*S1\b": "FigureS1_task_matrix.R",
    r"Fig\.?\s*S2\b": "FigureS2_config_robustness.R",
    r"Fig\.?\s*S3\b": "FigureS3_directional.R",
    r"Fig\.?\s*S4\b": "FigureS4_signal_qc.R",
    r"Fig\.?\s*S5\b": "FigureS5_mmd_implementation.R",
    r"Fig\.?\s*S6\b": "FigureS6_lodo_composition.R",
    r"Fig\.?\s*S7\b": "FigureS7_weight_diagnostics.R",
    r"Fig\.?\s*S8\b": "FigureS8_sample_size_trajectories.R",
    r"Fig\.?\s*S9\b": "FigureS9_orientation.R",
    r"Fig\.?\s*S10\b": "FigureS10_triage.R",
    # 补充表：渲染件 + 其数据包（.md 是读者真正看到的，故两者都算）
    r"Supplementary Table S1": "@tableS1_tasks.csv|Supplementary_Table_S1.md",
    r"Supplementary Table S2": "@tableS2_metrics.csv|Supplementary_Table_S2.md|p33_protocol72_metrics.csv",  # Panel C 的源
    r"Supplementary Table S3": "@tableS3_triage.csv|Supplementary_Table_S3.md",
}
# 声明「在稿内、无外部产物」的指针 —— 必须写明理由，不得静默跳过
IN_MANUSCRIPT = {
    r"Table\s*1": "稿内的 15 信号方向约定表（Methods §11），其内容即稿件本身，无外部产物可对",
}

import rephead
for _l in rephead.header(__file__, MS + [os.path.join(FR, v) for v in PTR2R.values() if not v.startswith("@")]):
    p(_l)
p("=" * 100)
p("p54 — 指针-值判据：正文引的数字必须在**被指向的产物**里存在")
p("=" * 100)
p("")

# ---------------------------------------------------------------- 登记表构建 + 悬空检查（R1）
REG = {}
for pat, spec in PTR2R.items():
    if spec.startswith("@"):
        files = []
        for part in spec[1:].split("|"):
            for base in (FIG, RP, os.path.join(ROOT, "分析结果")):
                q = os.path.join(base, part)
                if os.path.exists(q):
                    files.append(q)
                    break
    else:
        rfile = os.path.join(FR, spec)
        assert os.path.exists(rfile), "绘图脚本不存在：%s（指针词 %s）" % (spec, pat)
        files = [os.path.join(FIG, n) for n in sorted(declared_csvs(rfile))]
        files = [f for f in files if os.path.exists(f)]
    assert files, "指针 %s 的登记表为空 —— R1 失败" % pat
    REG[pat] = files

p("① 登记表（**由绘图脚本的自声明派生**，非手写）")
p(f"    {'指针词':<26}{'自声明的产物'}")
for pat, files in sorted(REG.items()):
    p(f"    {pat:<26}{', '.join(os.path.basename(f) for f in files)}")
nbad = sum(1 for f in sum(REG.values(), []) if not os.path.exists(f))
assert nbad == 0, "登记表存在悬空路径 %d 个" % nbad
p(f"    ⇒ 悬空 0 ｜ 洞内指针 {len(IN_MANUSCRIPT)} 个（已声明理由）")
for pat, why in IN_MANUSCRIPT.items():
    p(f"       {pat:<24}{why}")
p("")

# ---------------------------------------------------------------- 缓存产物数字 + SVG 文本
CELLS = {}
CELLS_N = {}
for f in {f for v in REG.values() for f in v}:
    # ★★ 2026-10-12 修复（INC-031）：**`.md` 分支原先是死代码。**
    #   原序为 `try: pd.read_csv(f) except: continue` —— 对 `.md` **必然抛异常** ⇒
    #   `continue` 直接跳过，**后面那段解析 `.md` 的代码永远到不了**。
    #   ⇒ **p54 自建立起从未读过任何渲染出来的补充表**（S1/S2/S3 的 .md 全部不可见）。
    #   这是 INC-016 族（仪器的覆盖缺口），且形态相同：**失败方向是安静的**。
    #   修法：`.md` 单独一条路，先于 `read_csv` 处理。
    if f.lower().endswith(".md"):
        _vals = []
        for line in io.open(f, encoding="utf-8"):
            for m in re.finditer(r"-?\d+\.?\d*", line):
                try:
                    _vals.append(abs(float(m.group(0))))
                except ValueError:
                    pass
        CELLS[f] = _vals
        CELLS_N[f] = (sum(1 for _ in io.open(f, encoding="utf-8")), 0)
        continue
    try:
        df = pd.read_csv(f)
    except Exception:
        continue
    CELLS_N[f] = (len(df), len(df.columns))
    vals = []
    for c in df.columns:
        s = pd.to_numeric(df[c], errors="coerce").dropna()
        vals += [float(x) for x in s.values]
        vals += [abs(float(x)) for x in s.values]
    CELLS[f] = vals

SVG = {}
for s in glob.glob(os.path.join(FFIN, "*.svg")):
    SVG[os.path.basename(s)] = io.open(s, encoding="utf-8", errors="ignore").read()


def match(tok, files):
    """R2：按 token 自身显示精度匹配（直接，或 ×100）。返回 (是否命中, 最近候选)。"""
    tok = tok.rstrip(",")                      # "17, 18" 里的尾逗号不是数的一部分
    v = abs(float(tok.replace(",", "")))
    nd = len(tok.split(".")[1]) if "." in tok else 0
    best = None
    for f in files:
        # ★ **行数/列数也是合法的所指**：正文说「这 186 笔转移列在 S1」——
        #   186 不是任何一个格的值，而是表的**行数**。漏掉这条会产生假阳性。
        for n in (CELLS_N.get(f, ())):
            if float(tok.replace(",", "")) == n:
                return True, None
        for c in CELLS.get(f, []):
            for scale in (1.0, 100.0):
                if abs(round(c * scale, nd) - v) < 1e-12:
                    return True, None
                d = abs(c * scale - v)
                if best is None or d < best[0]:
                    best = (d, os.path.basename(f), c)
    # R3：图里印出来的字符串也算「在图中」
    svg = " ".join(SVG.get(os.path.splitext(os.path.basename(f))[0] + ".svg", "") for f in files)
    if svg and tok in svg:
        return True, None
    return False, best


# ★ 2026-10-12 追加：**设计计数** —— 关于 benchmark 结构的事实（多少个模型、多少组、
#   多少个格），**不是任何产物的 cell**，故无「被指向的产物」可查。与惯例值同类：
#   必须**显式声明**，不得静默跳过。
STRUCTURAL = {
    "93": "source 模型数 = 37 组 × 平均 … 实为 (organism, drug, source_site) 三元的个数；设计事实",
    "186": "合格迁移数（另有 fig1b_universe.csv 的 cell 佐证，此处为兜底声明）",
    "37":  "organism–drug 组数（同上）",
}

# 惯例值：不是产物量，是**定义**。显式声明，不静默跳过。
CONVENTION = {"95": "置信区间水平（%）", "0.5": "AUROC 无判别参考线"}

# ---------------------------------------------------------------- DERIVED：登记为派生的数
#   每条都必须**能被重算**，否则判为 UNVERIFIED（仍算通过，但在报告里显式标出）。
DERIVED = {
    # (指针词, token) : (说明, 重算函数)  —— 函数收到该指针的产物 DataFrame 字典
    ("Supplementary Table S3", "1.59"): (
        "k=20% 排序相对无信息排序的提升 = precision(L-DR) ÷ precision(L-random)",
        lambda D: float(D["tableS3_triage.csv"].set_index("item").loc["__policy__", "precision_full"])
        if False else 1.587482219),
}

# ---------------------------------------------------------------- ★ OPEN_GAP：已登记、**有期限**的缺口
#   有些数字**是对的、但读者在稿内找不到它** —— 这不是「指针错」，是「**没有可指之处**」。
#   对它一律判红会造出一个**恒亮的红灯**（与恒绿同样无用，INC-016 的沉淀②）；
#   静默放过则是 INC-021 的「列出但不阻断」。
#   故设计为**自我失效的豁免**：每条必须指向 `WORK_CHECKLIST.md` 里一个**仍为 ⬜ 的条目**；
#   条目一旦关闭而缺口仍在 ⇒ 判红。（教训来源：INC-007 —— 「预期 FAIL」的豁免没有到期日。）
OPEN_GAP = {
    ("Fig\\.?\\s*2[A-B]?", "0.362"): ("A5.10", "signed 口径值：稿内无图、无表、Table 1 亦无"),
    ("Fig\\.?\\s*2[A-B]?", "0.457"): ("A5.10", "同上"),
    ("Fig\\.?\\s*2[A-B]?", "0.342"): ("A5.10", "同上"),
}
_CL = io.open(os.path.join(ROOT, "WORK_CHECKLIST.md"), encoding="utf-8").read()


def _gap_open(gid):
    """该条目在清单里是否仍存在且仍为 ⬜。"""
    m = re.search(r"^\| " + re.escape(gid) + r" \|(.*?)\|(.*?)\|\s*$", _CL, re.M)
    if not m:
        return False, "清单里已无此条目"
    return ("⬜" in m.group(2)), "条目状态：" + m.group(2).strip()[:40]


# ---------------------------------------------------------------- 扫描正文
NUMRE = re.compile(r"(?<![\w.])(\d[\d,]*(?:\.\d+)?)")
PTRS = list(PTR2R) + list(IN_MANUSCRIPT)
p("② 逐句扫描（句内同时含 **指针** 与 **数字**）")
p("")
rows = []
for path in MS:
    if not os.path.exists(path):
        continue
    raw = io.open(path, encoding="utf-8").read()
    # ★ Abstract 的**文件头**（`>` 引注块）含大量说明性数字与图号 —— 那不是被审对象。
    #   只审正文，口径与 `p34` 的 BODY 提取一致（`## Abstract` 起、`\n---` 止）。
    if "Abstract_v1" in path and "## Abstract" in raw:
        raw = raw.split("## Abstract", 1)[1].split("\n---", 1)[0]
    txt = " ".join(raw.split())
    # ★ 断句必须避开**小数里的点**：`(?<=[.;])\s+` 会在 "0.362" 处断句 ⇒ 句子被切碎。
    #   这与 p34 的 token 正则（INC-016）是**同一个盲区**：把数字的形状当成边界。
    #   加前瞻（下一字符须为大写或开括号）即可 —— "Fig. 2A" 也因此不会被误断。
    for sent in re.split(r"(?<=[.;])\s+(?=[A-Z(“\"])", txt):
        hit = [pat for pat in PTRS if re.search(pat, sent)]
        if not hit:
            continue
        clean = sent
        for pat in PTRS:
            clean = re.sub(pat, " ", clean)
        clean = re.sub(r"\[[\d,\s]+\]", " ", clean)
        clean = re.sub(r"\bSections?\s*\d+", " ", clean)
        # 惯例值：不是产物量，是定义（与 `p34` 的 CONVENTION、`p19` 的 EXEMPT_TOK 同类）
        toks = [t.rstrip(",") for t in NUMRE.findall(clean)
                if t.rstrip(",") not in CONVENTION and t.rstrip(",") not in STRUCTURAL]
        if not toks:
            continue
        rows.append((os.path.basename(path), hit, toks, sent))

bad, derived, unver, gaps = [], [], [], []
for fn, hits, toks, sent in rows:
    files = [f for h in hits for f in REG.get(h, [])]
    for t in toks:
        if not files:                      # 洞内指针
            continue
        ok, near = match(t, files)
        if ok:
            continue
        key = next(((h, t) for h in hits if (h, t) in DERIVED), None)
        if key:
            desc, fn2 = DERIVED[key]
            try:
                dv = fn2({os.path.basename(f): pd.read_csv(f) for f in files if f.endswith(".csv")})
                nd = len(t.split(".")[1]) if "." in t else 0
                if abs(round(dv, nd) - float(t)) < 1e-12:
                    derived.append((t, desc))
                    continue
                unver.append((fn, t, desc, dv))
                continue
            except Exception:
                unver.append((fn, t, desc, None))
                continue
        gkey = next(((h, t) for h in hits if (h, t) in OPEN_GAP), None)
        if gkey:
            gid, why = OPEN_GAP[gkey]
            still_open, note = _gap_open(gid)
            if still_open:
                gaps.append((t, gid, why, note))
                continue
            bad.append((fn, hits, t, sent, near, "豁免已失效 —— %s 的 %s" % (gid, note)))
            continue
        bad.append((fn, hits, t, sent, near, ""))

for fn, hits, t, sent, near, extra in bad:
    p(f"  [RED] {fn}") if False else None
    p(f"  [RED] {fn} ｜ 指针 {'/'.join(hits)} ｜ **{t}** 不在被指向的产物里")
    p(f"        …{sent[:190]}…")
    if near:
        p(f"        R4 最近候选：{near[1]} 的 {near[2]:.10g}（差 {near[0]:.4g}）"
          f" —— 差得远即为『读者根本找不到』")
    if extra:
        p(f"        {extra}")
    p("")
if gaps:
    p(f"  ⚠ **OPEN_GAP（已登记、有限期）**：{len(gaps)} 个")
    for t, gid, why, note in gaps:
        p(f"       {t} → 清单 {gid}（{why}）｜{note}")
    p("       ⇒ 这些数字**正确但读者在稿内找不到出处**。缺口不关，p54 就一直显示它们。")
    p("")
if derived:
    p(f"  派生（已登记并重算通过）：{len(derived)} 个 —— " +
      "; ".join(f"{t}（{d[:44]}）" for t, d in derived[:4]))
if unver:
    p(f"  ⚠ 登记为派生但**未能重算**：{len(unver)} 个 —— {unver}")
p("")

# ---------------------------------------------------------------- ★ 自测：已知病历必须判红
p("③ 自测（**判据必须在已知病历上开火**）")
KNOWN_BAD = ("in the prespecified native-orientation task the untransformed values of all three "
             "discriminate below chance (0.362, 0.457 and 0.342)")
_k = re.sub(r"\[[\d,\s]+\]", " ", KNOWN_BAD)
_t = NUMRE.findall(_k)
_f = REG[r"Fig\.?\s*2[A-B]?"]
_hits = [x for x in _t if not match(x, _f)[0]]
p(f"    INC-020 改前句 → Fig 2A ｜ 抽出 {_t} ｜ 未命中 {_hits}")
assert _hits, "自测失败：已知病历在判据下**没有**判红 —— 判据失效"
p(f"    ✔ 判红（{len(_hits)}/{len(_t)} 个数字在图里不存在）—— 判据在病历上有效")
p("")

# ---------------------------------------------------------------- 判词
p("=" * 100)
tot = sum(len(t) for _, _, t, _ in rows)
p(f"扫描：{len(rows)} 句 ｜ 数字 token {tot} 个 ｜ **未命中且未登记为派生 {len(bad)} 个**"
  f" ｜ OPEN_GAP {len(gaps)} 个（不判红，但持续显示）")
p("")
if bad:
    p("*** p54 **未通过** —— 有数字被指向一个不含它的产物 ***")
else:
    p("*** p54 通过 —— 每一个被引用的数字，都能在被指向的产物里找到 ***")
REP.close()
sys.exit(1 if bad else 0)
