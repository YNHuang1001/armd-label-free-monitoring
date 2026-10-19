# -*- coding: utf-8 -*-
r"""p34 — Abstract 审计器

================================ 为什么需要 ================================
`p15` 审 Results、`p16` 审 Methods、`p19` 审 Discussion、`p31` 审补充图例 ——
**Abstract 落在全部审计器之外**（实测：全项目 grep `Abstract_v1.md` 只有
`ms_master.py` 的路径助手与两个引用它的报告脚本，无任何审计器）。

这是全文被读最多的两百词，却比正文少一层锚定。

============================= 本审计器查什么 =============================
作者 2026-10-07 收敛为**五条可机械判定**的（不采纳外部建议中「必需声明在场」里
包含 "DR 三项最优" 的那半句 —— 该框架已被裁定不成立：DR 与 OM 在五项上
bootstrap 区间全部重叠，正文用词是 indistinguishable）：

  A. **数字溯源** —— Abstract 正文里**每一个**数字 token 都必须登记来源，
     且与冻结产物一致。未登记的数字 = FAIL（这是 p15 的同款逻辑；
     重点抓的是「摘要中独有、正文查无此数」的雷）。
     ★ 2026-10-08（INC-016）：**本条曾长期是空转的** —— 旧 token 正则
       `(?<![\w.])(\d[\d,]*)` 看不见小数（见 A 节内的详细病历）。修好后
       **覆盖面**才等于「每一个数字」；并加 **A-ST 非空自检**，正则再退化即响。
     ★ 判等按**显示精度**（`_val_ok`）：0.556 对 0.5558628841607566 判等 ——
       「显示精度相同 ≠ 精确算术相等 ≠ 比特级相等」，三档不可混。
     ★★ 2026-10-14（**本表第二次结构性重写**）：摘要换成两位小数版本后，
       登记表暴露**两个真缺陷**，都是「旧值恰好没踩到」才一直没响：
         (1) **同一 token 多次出现、含义不同**时，旧表以 token 为键、只存一条
             ⇒ 只查得到其中一处。新摘要里 `0.03` 出现 **3 次**（DR−src 的点
             估计一处，控制对照的下界与上界各一处），旧结构对此**结构性失明**。
             改为 token → **条目列表**，并加「出现次数 == 登记条数」断言。
         (2) **值在源码里手打**，不是从冻结产物读的 —— 违反本仓的出处纪律。
             改为**全部现场读取**（p57 两个 CSV + p57 报告的 fail-closed 抽取）。
  B. **词数** —— Nature Communications 硬上限 200；项目目标带 175–195。
  C. **禁用词** —— 复用 `p19` 的 BAN 表（first / novel / proves / causes / due to …）。
  D. **零引用** —— NC Article abstract 不得含 references。
  E. **家族断言带限定** —— 抽象的家族级比较必须带口径限定语
     （2026-10-07 新增；此前该句无任何限定，而家族均值在 grouped 口径下反转）。
     ★ 2026-10-14：新摘要**不再包含任何家族级比较句**。E-1 改为**条件式**判据
       （出现了家族级优势句就必须带限定；没有则不触发），E-2 改为 **WARN**，
       如实记录「L2 家族在 grouped 口径下失效」这一披露**当前不在摘要里** ——
       这是作者的内容决定，审计器**不代改**，但**必须让它在报告里可见**。

============================= 状态三分 ==================================
PASS / WARN / FAIL。**只有 FAIL 影响退出码**；WARN 是「需要人决定的披露缺口」，
必须与静默通过区分开（否则审计器会把「没查」和「查了没问题」印成同一行）。
"""
import io, os, re, sys
import pandas as pd

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIG  = os.path.join(ROOT, "figures", "data")
OUT  = os.path.join(ROOT, "分析结果")
# ★ 正本认定见 MANUSCRIPT_MASTER_DECLARATION.md §2：Abstract 的正本在**根目录**
# 可选 argv[1] = 待审文件路径（供负向测试指向改动副本，不碰正本）
ABS  = sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, "Abstract_v1.md")
# ★ 2026-10-14（INC-039）：**给了 argv[1] 就是负向测试，产物必须落 TEST 后缀。**
#   同族病历见 `p15_results_numeric_audit.py`：负向测试覆写了生产报告，在证据面上留下**假红**。
_TEST = len(sys.argv) > 1
_SUF  = "_TEST" if _TEST else ""
if _TEST:
    print("*** 负向测试模式（argv[1]=%s）—— 产物写 %s 后缀，不碰生产报告 ***" % (ABS, _SUF))
REP  = open(os.path.join(OUT, "report_p34_abstract%s.txt" % _SUF), "w", encoding="utf-8")

def p(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); REP.write(s + "\n"); REP.flush()

TXT_ALL = io.open(ABS, encoding="utf-8").read()
BODY = TXT_ALL.split("## Abstract", 1)[1].split("\n---", 1)[0]
ROWS = []
def rec(cid, kind, item, ok, note="", status=None):
    """status=None ⇒ 由 ok 推 PASS/FAIL；显式传 "WARN" 表示「不是缺陷，但需人决定」。"""
    st = status or ("PASS" if ok else "FAIL")
    ROWS.append(dict(cid=cid, kind=kind, item=item, status=st, note=note))
    p("  [%-4s] %-4s %-52s %s" % (st, cid, item[:52], note))
    return st == "PASS"

R = lambda f: pd.read_csv(os.path.join(FIG, f))

# ★ 2026-10-14：p57 的三个产物加入出处清单 —— 新摘要的 20 个数字里有 7 个由 p57 生产，
#   而它们此前**不在本审计器的出处清单内**（旧摘要没有这些数，所以一直没暴露）。
P57D = os.path.join(OUT, "p57_paired_differences.csv")
P57S = os.path.join(OUT, "p57_scores.csv")
P57R = os.path.join(OUT, "report_p57_posthoc_baselines.txt")

import rephead
for _l in rephead.header(__file__, [ABS, os.path.join(FIG, "fig1b_universe.csv"),
        os.path.join(FIG, "fig2_forest.csv"), os.path.join(FIG, "fig5_subsample.csv"),
        P57D, P57S, P57R]):
    p(_l)
p("=" * 104); p("p34 — Abstract 审计器"); p("=" * 104)
p(f"正本：{ABS}")
p(f"（MANUSCRIPT_MASTER_DECLARATION.md §2：Abstract 正本在根目录，manuscript_renumbered/ 无对应文件）")
p("")

# ================================================================ A 数字溯源
p("【A】数字溯源 —— 正文每个数字 token 必须登记来源")
U  = R("fig1b_universe.csv").iloc[0]
F5 = R("fig5_subsample.csv")
LEVELS = sorted({str(v) for v in F5.level.unique()})
F2  = R("fig2_forest.csv").set_index("method")


def g2(m, c):
    return float(F2.loc[m, c])


# ★★ 2026-10-08 更正（INC-016）：**旧 token 正则在结构上看不见任何小数。**
#   旧式 `(?<![\w.])(\d[\d,]*)`：在 "0.439" 上先命中 `0`（负向断言通过），
#   随后 `439` 的**前导字符是 "."** ⇒ 负向断言失败 ⇒ 整段被丢弃。
#   实测后果：v4 摘要里六个小数被并成一个「出现 6 次的 token 0」，
#   `1.59` 被报成 token `1` —— **判据从未看过摘要里任何一个小数**，
#   而本文摘要的数字**全部**落在第三位小数上。
#   与 INC-008/012/013/014 同族（仪器静默漏计），但**量级不同**：
#   那几次漏的是**一条**，这次漏的是**一整类 token**；且失败形态是**红灯**
#   （「未登记的 0」），看的人只会以为「有个数字没登记」，
#   **不会想到「小数根本没进扫描」**。⇒ 沉淀：**覆盖面缺陷要配非空自检**。
TOK_RE = re.compile(r"(?<![\w.])(\d[\d,]*(?:\.\d+)?)")

# ---- 非空自检：证明扫描器确实看得见小数（正则再被改坏时立刻响）----
_SELFTEST = {"0.439": "(AUROC was 0.439", "1.59": "(1.59 times", "186": "(186 met"}
for _tok, _frag in _SELFTEST.items():
    _got = TOK_RE.findall(_frag)
    assert _tok in _got, f"p34 token 扫描自检失败：{_frag!r} 未取到 {_tok!r}，实得 {_got}"
rec("A-ST", "A", "token 扫描器非空自检（小数可见）", True,
    " ／ ".join(f"{k}←{v}" for k, v in _SELFTEST.items()))

# ---------------------------------------------------------------- 冻结产物现场读取
# ★ 出处纪律：**登记表里的值一律现场从冻结产物读**，不手打。手打的值与产物脱节时，
#   本审计器会变成「用自己的记忆核对别人的记忆」，红灯绿灯都不再有意义。
_PD = pd.read_csv(P57D).set_index("contrast")
_PS = pd.read_csv(P57S).set_index("score")
def p57d(c, col):
    return float(_PD.loc[c, col])
def p57s(c, col):
    return float(_PS.loc[c, col])

# ICC(1,1) 目前**只有报告文本、没有 CSV 产物** ⇒ 从报告 fail-closed 抽取
# （抽不到即抛错，绝不回退到默认值 —— 同 p21f 的 grab 纪律）。
_REP57 = io.open(P57R, encoding="utf-8").read()
def grab57(pat, what):
    m = re.search(pat, _REP57)
    if not m:
        raise RuntimeError("p57 报告抽取失败（%s）：pattern=%r —— 报告格式可能已改，"
                           "须人工核对后更新 pattern，禁止回退默认值" % (what, pat))
    return float(m.group(1))

ICC_PT = grab57(r"ICC\(1,1\)\s*=\s*([\d.]+)", "ICC 点估计")
ICC_LO = grab57(r"ICC\(1,1\)\s*=\s*[\d.]+\s*95%\s*区间\s*\[([\d.]+)", "ICC 下界")
ICC_HI = grab57(r"ICC\(1,1\)\s*=\s*[\d.]+\s*95%\s*区间\s*\[[\d.]+,\s*([\d.]+)\]", "ICC 上界")

# ---------------------------------------------------------------- 登记表
# 结构：token → **条目列表**（每条形如 (来源, 值, 出处句)）。
# ★ 2026-10-14：旧表以 token 为键、只容一条 ⇒ 同一 token 多次出现且含义不同时
#   **只查得到一处**。摘要里 `0.03` 出现 3 次（DR−src 点估计 1 处、控制对照两界 2 处），
#   旧结构对此完全失明。改为列表后，另加「出现次数 == 登记条数」断言兜住漏登记。
REG = {
    "186": [("fig1b_universe.csv:eligible_tasks", float(U.eligible_tasks),
             "benchmark of 186 directed cross-system transfers")],
    "37":  [("fig1b_universe.csv:n_organism_drug_groups", float(U.n_organism_drug_groups),
             "37 organism–drug susceptibility tasks")],
    "0.439": [("fig2_forest.csv:js_max.harmful_auc", g2("js_max", "harmful_auc"),
               "(AUROC 0.439 to 0.478")],
    "0.478": [("fig2_forest.csv:domain_auc.harmful_auc", g2("domain_auc", "harmful_auc"),
               "(AUROC 0.439 to 0.478")],
    "0.644": [("fig2_forest.csv:AUROC_DR.harmful_auc", g2("AUROC_DR", "harmful_auc"),
               "reached an AUROC of 0.644 (0.56 to 0.73)")],
    # ★ 区间正本 = fig2_forest.csv（Fig 2 所绘 / Table S2 所印 / §2 所引同源）。
    #   曾一度写成 0.557–0.729，那是 **lo 取 §7.2、hi 取图表** 的拼配值，不是任何单一产物
    #   的值；两个生产者的种子不同（见 INC-015）。
    #   2026-10-14：按时序稳定性（5 种子端点最大漂移 0.0081 > 0.005）改为**两位小数**显示
    #   ⇒ 判等仍按**显示精度**，故 0.56 对 0.555863、0.73 对 0.729476 判等。
    "0.56": [("fig2_forest.csv:AUROC_DR.auc_lo", g2("AUROC_DR", "auc_lo"),
              "reached an AUROC of 0.644 (0.56 to 0.73)")],
    "0.73": [("fig2_forest.csv:AUROC_DR.auc_hi", g2("AUROC_DR", "auc_hi"),
              "reached an AUROC of 0.644 (0.56 to 0.73)")],
    "0.611": [("p57_scores.csv:src_only.point", p57s("src_only", "point"),
               "baseline using source discrimination alone (0.611")],
    "0.03": [
        ("p57_paired_differences.csv:prespecified_DR − src_only.point",
         p57d("prespecified_DR − src_only", "point"),
         "difference +0.03, −0.14 to about 0.2"),
        ("p57_paired_differences.csv:posthoc_dhat_DR − amean_control.lo95",
         p57d("posthoc_dhat_DR − amean_control", "lo95"),
         "(+0.00; −0.03 to +0.03)"),
        ("p57_paired_differences.csv:posthoc_dhat_DR − amean_control.hi95",
         p57d("posthoc_dhat_DR − amean_control", "hi95"),
         "(+0.00; −0.03 to +0.03)"),
    ],
    "0.14": [("p57_paired_differences.csv:prespecified_DR − src_only.lo95",
              p57d("prespecified_DR − src_only", "lo95"),
              "−0.14 to about 0.2")],
    # ★ 上界写「约 0.2」而非 0.1915 的四舍五入值：5 个重采样种子上该端点最高到 +0.1984
    #   ⇒ 第二位小数会翻（见 SEED_STABILITY_TABLE_20261014.md）。按显示精度 0.2 判等。
    "0.2": [("p57_paired_differences.csv:prespecified_DR − src_only.hi95",
             p57d("prespecified_DR − src_only", "hi95"),
             "to about 0.2)")],
    "0.789": [("p57_scores.csv:posthoc_dhat_DR.point", p57s("posthoc_dhat_DR", "point"),
               "aligned with the change endpoint (0.789)")],
    "0.00": [("p57_paired_differences.csv:posthoc_dhat_DR − amean_control.point",
              p57d("posthoc_dhat_DR − amean_control", "point"),
              "(+0.00; −0.03 to +0.03)")],
    "0.49": [("p57 报告:ICC(1,1) 点估计", ICC_PT, "(ICC 0.49, 0.32 to 0.63)")],
    "0.32": [("p57 报告:ICC(1,1) 区间下界", ICC_LO, "(ICC 0.49, 0.32 to 0.63)")],
    # ★ ICC 上界 0.6346 距 0.63/0.64 的取整边界仅 0.0004 ⇒ 换一个重采样种子就可能显示成
    #   0.64。摘要按两位小数显示 0.63 是**当前种子下的正确显示**；不稳定性已在 Methods
    #   声明（第二位小数不可依赖），故此处不改措辞、只如实登记。
    "0.63": [("p57 报告:ICC(1,1) 区间上界", ICC_HI, "(ICC 0.49, 0.32 to 0.63)")],
    # ★ sens_at_spec80 是**比例**，摘要写的是**百分数** ⇒ 这里乘 100 后再按整数精度判等。
    #   0.385965 × 100 = 38.5965 → 39。该值同时是 15 个信号里的最大值（次高 AUROC_OM 0.368421），
    #   故摘要的 "the best signal" 成立（见 E-5）。
    "39": [("fig2_forest.csv:AUROC_DR.sens_at_spec80 × 100",
            g2("AUROC_DR", "sens_at_spec80") * 100,
            "the best signal identified 39% of failing transfers")],
}
# 惯例值：不是产物量，是**定义**。显式声明，免得「未登记」的红灯掩盖真问题。
CONVENTION = {
    "0.5": "AUROC 无判别参考线（区间是否跨过 0.5）",
    "0.05": "ΔAUROC ≤ −0.05 —— 预注册的 harmful 端点定义，非估计量",
    "80": "预设操作点 0.80 特异度（方案定义的操作点，非估计量）",
}


def _nd(tok):
    """token 自身的小数位数 —— 判据钉在**正文写的那个精度**上。

    规则：显示精度相同 ≠ 精确算术相等 ≠ 比特级相等。故 0.556 应对
    `fig2_forest` 的 0.5558628841607566 **按三位小数判等**，而不是逐位判等。
    """
    return len(tok.split(".")[1]) if "." in tok else 0


def _val_ok(tok, val):
    """★ 2026-10-14 修：**token 抽取不含符号**（`TOK_RE` 只捕数字），
    ⇒ 本表按**绝对值**判等；**符号由 `REG` 的「出处句」承担**
    （该句必须原样出现在摘要里，含 U+2212 负号）。"""
    return abs(abs(round(float(val), _nd(tok))) - abs(float(tok.replace(",", "")))) < 1e-12


toks = {}
for m in TOK_RE.finditer(BODY):
    toks[m.group(1)] = toks.get(m.group(1), 0) + 1
FLAT_A = " ".join(BODY.split())
for t, n in sorted(toks.items()):
    if t in CONVENTION:
        rec("A-" + t, "A", f"数字 {t}（出现 {n} 次）· 惯例值", True, CONVENTION[t])
    elif t in REG:
        ents = REG[t]
        # ★ 出现次数 == 登记条数：漏登记一处即响（旧表结构性看不见这一层）
        rec("A-" + t + "-n", "A", f"数字 {t} 出现 {n} 次 == 登记 {len(ents)} 条", n == len(ents),
            "" if n == len(ents) else "少登记或多登记 ⇒ 有含义未登记的用法")
        for i, (src, val, pat) in enumerate(ents, 1):
            tag = "" if len(ents) == 1 else f"·{i}"
            rec("A-" + t + tag, "A", f"数字 {t}{tag} 出处句在场", pat in FLAT_A, f"来源 {src}")
            rec("A-" + t + tag + "v", "A", f"数字 {t}{tag} 与冻结产物一致（按显示精度）",
                _val_ok(t, val), f"正文 {t} ｜ 全精度 {val:.10g} ｜ 舍入 {round(val, _nd(t)):g}")
    else:
        rec("A-" + t, "A", f"**未登记的数字 {t}**（出现 {n} 次）", False,
            "Abstract 中出现的数字必须在 REG 里登记来源，或声明为惯例值")
p(f"  共发现数字 token {len(toks)} 种：{', '.join(sorted(toks))}")
p(f"  已登记产物 {len(REG)} 种（共 {sum(len(v) for v in REG.values())} 条）｜ "
  f"声明惯例 {len(CONVENTION)} 种 ｜ fig5_subsample.csv 的 level 档位 = {LEVELS}")
p(f"  现场读取：p57 点估计/区间 {len(_PD)} 行、p57 汇总得分 {len(_PS)} 行、"
  f"ICC(1,1) = {ICC_PT} [{ICC_LO}, {ICC_HI}]（报告文本 fail-closed 抽取）")
p("")

# ================================================================ B 词数
p("【B】词数")
# ★ 2026-10-08 更正（INC-016 第二处，与 A 节同一根因）：旧式以**字母数字串**为单位，
#   于是 ① 每个小数被拆成两「词」（"0.439" → "0" + "439"）；
#   ② `organism–drug` 的连接号是 **U+2013**（不在字符类内）⇒ 被拆成两「词」。
#   实测：真实 191 词的正文被报成 **199**，**逼近 200 上限** ⇒ 会逼出无谓的压缩。
#   根因与 A 节同：**把「数字/标点的形状」当成了词边界**。改为空白切分。
w = len(BODY.split())
rec("B-0", "B", "词数口径 = 空白切分（与投稿系统一致，写死）", True,
    "小数、连字符、en dash 均计 1 词")
rec("B-1", "B", "词数 <= 200（目标刊 npj Digital Medicine；沿用 Nature Portfolio 200 词口径）",
    w <= 200, f"实际 {w}")
rec("B-2", "B", "词数落在项目目标带 175-195", 175 <= w <= 195, f"实际 {w}")
p("")

# ================================================================ C 禁用词
p("【C】禁用词（复用 p19 的 BAN 表）")
BAN = {
 "first":             r"\bfirst\b",
 "novel":             r"\bnovel\b|\bnovelty\b",
 "proves":            r"\bproves?\b|\bproven\b",
 "demonstrates that": r"demonstrat\w+ that",
 "causes":            r"\bcauses?\b|\bcaused by\b",
 "due to":            r"\bdue to\b",
 "drives":            r"\bdrives?\b",
 "as a result of":    r"\bas a result of\b",
 "we are the first":  r"we are the first",
 "unprecedented":     r"\bunprecedented\b",
}
for name, pat in BAN.items():
    hits = re.findall(pat, BODY, re.I)
    rec("C", "C", f"禁用词 {name!r}", not hits, f"命中 {len(hits)} 次" if hits else "")
p("")

# ================================================================ D 零引用
p("【D】零引用")
cits = re.findall(r"\[\d+(?:\s*[,\u2013-]\s*\d+)*\]", BODY)
rec("D-1", "D", "不含文献引用标记 [n]", not cits, f"命中 {cits}" if cits else "")
p("")

# ================================================================ E 家族断言带限定
p("【E】家族断言带限定")
FLAT = " ".join(BODY.split())

# E-1：**条件式**判据。危险形态是「无口径限定的家族级优势句」；
#   摘要若根本不作家族级比较，则该危险不存在，判据不触发。
#   ★ 不写成「必须含某句」——那是把某版的措辞钉死成判据，换一版措辞就误报
#     （2026-10-14 实测：旧 E-1/E-2 正是这么红了 26 条的）。
FAM_CLAIM = re.compile(
    r"(?:input|prediction)[\s-]space[\w\s-]{0,40}?"
    r"(?:better|outperform\w*|superior|more information|retained more)", re.I)
fam_hits = FAM_CLAIM.findall(FLAT)
QUAL = ("in that analysis", "grouped evaluation", "indistinguishable")
fam_ok = (not fam_hits) or any(q in FLAT for q in QUAL)
rec("E-1", "E", "家族级优势句若在场，必带口径限定（in that analysis / grouped evaluation）",
    fam_ok, f"命中家族优势句 {len(fam_hits)} 处" if fam_hits else "摘要无家族级优势句 ⇒ 不触发")

# E-2：**披露缺口，作者决定项** ⇒ WARN，不进退出码，但必须在报告里可见。
#   本项记录的是：Results 的家族级发现（L2 在 pooled 口径下信息更多、在 grouped 口径下
#   不成立）**当前是否出现在摘要里**。摘要不含该句时，读者看不到「测过 prediction-space
#   家族」这件事；这是内容决定，审计器不代改。
L2_DISC = "did not survive grouped evaluation" in FLAT
rec("E-2", "E", "L2 家族在 grouped 口径下失效 —— 该披露在摘要中是否在场",
    L2_DISC, "在场" if L2_DISC else
    "**不在场**：摘要现不含任何 prediction-space 家族级陈述（仅报告 input-space 与 DR）。"
    "结果与限制仍在 Results/Fig 3；是否补入摘要属作者的内容决定。",
    status="PASS" if L2_DISC else "WARN")

rec("E-3", "E", "未出现「DR 三项最优」类无区间支撑的排序声称",
    not re.search(r"\bthree of (?:the )?(?:four|five)\b|\bbest on three\b", FLAT, re.I),
    "DR 与 OM 五项 bootstrap 区间全部重叠，正文用词是 indistinguishable")

# E-4：摘要里的 input-space 句必须带**范围词**，否则会被读成「所有 label-free 信号」
rec("E-4", "E", "input-space 结论带范围词（限定在 input-space，不外推到全部信号）",
    "Input-space shift statistics showed no positive discrimination" in FLAT,
    "摘要句首即限定 Input-space")

# E-5：摘要称 "the best signal"，须核对该信号确为 15 个信号中该指标的最优者
_sens = F2["sens_at_spec80"].astype(float)
_best = F2.index[int(_sens.values.argmax())]
rec("E-5", "E", "「the best signal」所指方法确为 15 个信号中 sens@80%spec 最优者",
    _best == "AUROC_DR",
    f"实测最优 = {_best}（{_sens.max():.6f}），摘要所述量出自该方法")

# ---- E-6 / E-7：v5（2026-10-14 外部核验稿）**删去的两处 v4 内容** ----
# 两条都**不影响任何数字的正确性**，只影响**覆盖面与措辞强度** ⇒ WARN，不阻断。
# 记录它们是因为：这两处在 v4 阶段**已经被论证过**（见 `_abs_draft_formatted.md` 的归档件、
# WORK_CHECKLIST A5.3a），现在被静默删掉；若审计器不响，删除就**不留痕**。

# E-6：v4 的「第 4 交付物」（分诊策略模拟）是否仍在摘要里
_policy = ("policy simulation" in FLAT) or ("reviewing the top" in FLAT)
rec("E-6", "E", "v4 的分诊策略模拟句是否仍在摘要中",
    _policy, "在场" if _policy else
    "**不在场**：v4 的 `In a retrospective policy simulation, reviewing the top 20% …` 已删。"
    "预设警示任务的**单点**仍在（39% of failing transfers），但**分诊那条对比不在摘要**。"
    "该对比仍在 Results / Figure S10 / Table S3。是否补回属作者的内容决定。",
    status="PASS" if _policy else "WARN")

# E-7：`poorly calibrated` 是否带 `in absolute terms` 限定。
#   v4 的文件头**专门论证过**为何必须带：冻结证据的形状是**序数保真（DR Spearman 0.7896）
#   ＋ 绝对水平差（MAE 0.0416）**，裸写 `poorly calibrated` 比证据更重（两个方向都不准）。
#   v5 把限定语丢了 ⇒ 这里是**回归**，故单独设一项盯着。
_calib = "poorly calibrated" in FLAT
_qual  = "in absolute terms" in FLAT
rec("E-7", "E", "校准措辞若为 `poorly calibrated`，须带 `in absolute terms` 限定",
    (not _calib) or _qual,
    ("带限定" if _qual else
     "**限定语丢失**：现为裸写 `poorly calibrated`。冻结证据是序数保真（DR Spearman 0.7896）"
     "＋绝对水平差（MAE 0.0416）⇒ 裸写比证据更重（v4 文件头已论证过一次）。措辞属作者决定。"),
    status="PASS" if ((not _calib) or _qual) else "WARN")
p("")

# ================================================================ 汇总
p("=" * 104)
D = pd.DataFrame(ROWS)
bad  = D[D.status == "FAIL"]
warn = D[D.status == "WARN"]
p(f"p34 断言合计 {len(D)} 项：PASS {int((D.status=='PASS').sum())}｜"
  f"**WARN {len(warn)}**｜**FAIL {len(bad)}**")
if len(warn):
    p("")
    p("  ── WARN（不是缺陷，是**需要人决定的披露缺口**）──")
    for _, r in warn.iterrows():
        p(f"    [WARN] {r['cid']} {r['item']}")
        p(f"           {r['note']}")
if len(bad):
    p("")
    for _, r in bad.iterrows():
        p(f"    [FAIL] {r['cid']} {r['item']}  {r['note']}")
    p(""); p("*** p34 未通过 ***")
else:
    p("")
    p("*** p34 通过 —— Abstract 的每个数字可溯源，且五类检查全部成立 ***")
    if len(warn):
        p(f"*** 但含 {len(warn)} 条 WARN —— 见上，属作者决定项，不阻断 ***")
D.to_csv(os.path.join(OUT, "p34_abstract_audit%s.csv" % _SUF), index=False, encoding="utf-8-sig")
REP.close()
sys.exit(1 if len(bad) else 0)
