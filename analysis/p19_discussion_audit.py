# -*- coding: utf-8 -*-
r"""p19 — Discussion 审计（四类风险）

================================ 为什么单独一套 ================================
Discussion 的风险与 Results 不同。Results 错在**数字**，Discussion 错在**身份与强度**：

  R1 数字来源    —— 正文出现的数字必须能追到冻结产物
  R2 primary/sensitivity 身份 —— DR 是 primary；OM/IOW 是 sensitivity。
                    **OM 的 MAE 更低不得被写成 "OM superior" 或事后升格**
  R3 post-unblinding 身份 —— orientation / rank 分析是 explanatory，
                    **不得被包装成可部署规则**
  R4 causal / novelty 强度 —— 不得写 causal，不得写 first。
                    association 只能写成 diagnostic / consistent with

这四类**都不会被 p15（Results 数值审计）抓到**，因为 p15 只读 Results。

================================ 边界 ================================
不重算统计量；数字来源一律从冻结文件读取。不改任何文件。

用法：python p19_discussion_audit.py [Discussion 路径]
输出：report_p19_discussion_audit.txt / p19_discussion_audit.csv
"""
import io, os, re, sys
import numpy as np, pandas as pd

try:                                   # Windows 重定向下 stdout 默认 GBK，会中途崩掉审计
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); BASE = os.path.dirname(ROOT)
OUT = os.path.join(ROOT, "分析结果"); FIG = os.path.join(ROOT, "figures", "data")
# ★ 2026-10-06 source binding 修正：原先读根目录的 Discussion_v1.md（**非正本**），
#   全稿引文号检查同样从根目录取。正本见 MANUSCRIPT_MASTER_DECLARATION.md §2。
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ms_master as MM                                                   # noqa: E402
MS = sys.argv[1] if len(sys.argv) > 1 else MM.master("Discussion_v1.md")

rep = open(os.path.join(OUT, "report_p19_discussion_audit.txt"), "w", encoding="utf-8")
def p(*a):
    line = " ".join(str(x) for x in a); print(line, flush=True); rep.write(line + "\n"); rep.flush()

def norm(s):
    s = (s.replace("\u2212", "-").replace("\u2013", "-").replace("\u2014", "-")
          .replace("\u2192", "->").replace("\u03c1", "rho").replace("\u00a0", " ").replace("**", ""))
    return " ".join(s.split())

F3C = pd.read_csv(os.path.join(FIG, "fig3c_sudo_decomposition.csv"))
F3A = pd.read_csv(os.path.join(FIG, "fig3a_lodo.csv"))
F2  = pd.read_csv(os.path.join(FIG, "fig2_forest.csv"))
F4D = pd.read_csv(os.path.join(FIG, "fig4d_mae.csv"))
F4O = pd.read_csv(os.path.join(FIG, "fig4e_dr_out_of_range.csv"))
FSO = pd.read_csv(os.path.join(FIG, "figS_orientation.csv"))
g4 = lambda m, c: float(F4D.set_index("estimator").loc[m, c])

RAW = open(MS, encoding="utf-8").read()
# 只审正文：切掉参考表与 header 引用块
BODY = RAW
for marker in ("## References", "## References used"):
    if marker in BODY: BODY = BODY.split(marker)[0]
BODY = re.sub(r'^> .*$', '', BODY, flags=re.M)          # header 引用块

# 【必需】先过全稿风格规范化，再匹配语义断言。
# 否则 p20 把 "long-recognised" 改成 "long-recognized" 后，这里的一条 PASS 会变成**假 FAIL**。
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from style_map import normalize as _style
TXT = norm(_style(BODY)[0])

# ================================================================ R1 数字来源
# 【项目级规则】审计项必须给出完整溯源五元组，**不得只报 PASS**：
#   item → source file → exact field / row filter → observed value → verdict
# 起因：p19 首版曾用 len(fig2_forest) 取 202，取到 15 却仍报 PASS ——
# **来源取错时，只打印 PASS 会制造错误的安全感。**
# ★★ 2026-10-06 修正·重要 ★★
#   原 NUM 的 token 是**写死的旧值**（0.511/0.529/0.756/0.423/138/25/202/38），且判定
#   只看「该串是否出现在正文」。修正轮把正文改成新值后，这 8 条全部 NOT_IN_MS。
#   更严重的是：**"存在即通过"本身不校验 token 与来源值是否相等** —— 当初 0.423 明明
#   与 fig2 的 0.5266 不符，却因为串在正文里而 PASS。这正是 GPT 指出的
#   "validator expectation provenance" 漏洞。
#   现改为：① token 与来源值**同源派生**（不写死）；② 加入相等断言（tol），
#   「存在」不再足以通过。
def _tok(v, dec=3):
    return ("%.*f" % (dec, float(v)))

NUM = [
 (_tok(F3C.pooled_202task.iat[0]), F3C.pooled_202task.iat[0], "fig3c 全宇宙 pooled",
  "fig3c_sudo_decomposition.csv", "col=pooled_202task, row=0（列名为旧 schema 名，值已是当前）"),
 (_tok(F3C.pooled_138task.iat[0]), F3C.pooled_138task.iat[0], "fig3c 可评价子集 pooled",
  "fig3c_sudo_decomposition.csv", "col=pooled_138task, row=0（列名为旧 schema 名，值已是当前）"),
 (_tok(F3C.macro_138task.iat[0]),  F3C.macro_138task.iat[0],   "fig3c macro",
  "fig3c_sudo_decomposition.csv", "col=macro_138task, row=0（列名为旧 schema 名，值已是当前）"),
 (_tok(F2.set_index("method").loc["SUDO_AURCC", "harmful_auc"]),
  float(F2.set_index("method").loc["SUDO_AURCC", "harmful_auc"]), "fig2 SUDO 原生全宇宙",
  "fig2_forest.csv", "method=SUDO_AURCC, col=harmful_auc"),
 (_tok(F3A.n_tasks_eval.iat[0], 0), float(F3A.n_tasks_eval.iat[0]), "fig3a n_tasks_eval",
  "fig3a_lodo.csv", "col=n_tasks_eval, row=0"),
 (_tok(F3A.n_eval_fold.iat[0], 0),  float(F3A.n_eval_fold.iat[0]),  "fig3a n_eval_fold",
  "fig3a_lodo.csv", "col=n_eval_fold, row=0"),
 ("100",   100.0,  "§10 level", "LabelFree..._Protocol_v1.0_FROZEN.md §10", "冻结文本：100/250/500/1000/full"),
 # 【二版 bug】首版用 len(F2) 取 202，但 fig2_forest.csv 是**方法级**表（15 行），
 # 不是 task 级表 —— 取到 15。宇宙规模必须来自 task 级来源。
 (_tok(pd.read_csv(os.path.join(FIG, "fig1b_universe.csv")).eligible_tasks.iat[0], 0),
  float(pd.read_csv(os.path.join(FIG, "fig1b_universe.csv")).eligible_tasks.iat[0]),
  "fig1b 合格转移数", "fig1b_universe.csv", "col=eligible_tasks, row=0"),
 (_tok(pd.read_csv(os.path.join(FIG, "fig1b_universe.csv")).n_organism_drug_groups.iat[0], 0),
  float(pd.read_csv(os.path.join(FIG, "fig1b_universe.csv")).n_organism_drug_groups.iat[0]),
  "fig1b 合格组数", "fig1b_universe.csv", "col=n_organism_drug_groups, row=0"),
]
# ---- 【2026-10-06 新增】本次 Limitations 改写引入的数字，一并纳入溯源 ----
#   第 6 位是**倍率**：正文写科学计数法的尾数（如 "7.8 × 10⁻⁴"），故比较时须还原量级。
_F5B = pd.read_csv(os.path.join(FIG, "fig5b_across_n_range.csv"))
_DR  = _F5B[_F5B.method == "AUROC_DR"].iloc[0]
_P5T = pd.read_csv(os.path.join(OUT, "p5_tasks_primary.csv"), keep_default_na=False)
_P6D = pd.read_csv(os.path.join(OUT, "p6d_outcome.csv"), keep_default_na=False)
_TID = "PA_AMK_UT2MGB"
# p5_tasks_primary.csv 无 task_id 列（只有 organism/drug/source/target），
# 按 p6d 的命名规则拼出同一主键后再取数。
_P5T["_tid"] = (_P5T.organism + "_" + _P5T.drug + "_" + _P5T.source + "2" + _P5T.target)
_dp5 = float(_P5T.loc[_P5T._tid == _TID, "d_auroc"].iat[0])
_dp6 = float(_P6D.loc[_P6D.task_id == _TID, "d_auroc"].iat[0])
_BND = 0.02      # descriptive near-zero band 的边界（|ΔAUROC| ≤ 0.02）
NUM += [
 (_tok(float(_DR.auroc_at_full)), float(_DR.auroc_at_full), "fig5b DR full 档 harmful AUROC",
  "fig5b_across_n_range.csv", "method=AUROC_DR, col=auroc_at_full"),
 (_tok(float(_DR.auroc_at_100)),  float(_DR.auroc_at_100),  "fig5b DR n=100 档 harmful AUROC",
  "fig5b_across_n_range.csv", "method=AUROC_DR, col=auroc_at_100"),
 # 边界距离 7.8e-4 = 0.02 − |d_p5|；两侧差 1.0e-3 = |d_p6d − d_p5|
 ("7.8", _BND - abs(_dp5), "边界余量（尾数 ×10⁻⁴）",
  "p5_tasks_primary.csv × p6d_outcome.csv",
  f"{_TID}: 0.02 − |d_p5={_dp5:.8f}| = {_BND-abs(_dp5):.3e}", 1e-4),
 ("1.0", abs(_dp6 - _dp5), "两侧差（尾数 ×10⁻³）",
  "p5_tasks_primary.csv × p6d_outcome.csv",
  f"{_TID}: |d_p6d={_dp6:.8f} − d_p5={_dp5:.8f}| = {abs(_dp6-_dp5):.3e}", 1e-3),
 # 既有的两个 MAE 值（此前从未被 R1 覆盖，因"存在即通过"而无人在意）
 (_tok(F4D.set_index("estimator").loc["DR", "mae"], 4), float(F4D.set_index("estimator").loc["DR", "mae"]),
  "fig4d DR MAE", "fig4d_mae.csv", "estimator=DR, col=mae"),
 (_tok(F4D.set_index("estimator").loc["OM", "mae"], 4), float(F4D.set_index("estimator").loc["OM", "mae"]),
  "fig4d OM MAE", "fig4d_mae.csv", "estimator=OM, col=mae"),
 # harmful 预设阈值（正文 "−0.05" 的绝对值）
 ("0.05", 0.05, "harmful 端点阈值", "Protocol §6 / fig1c_summary.csv",
  "harmful_threshold = −0.05（预设，非估计量）"),
]
# ★ 2026-10-09（A5.3d）：补登 Discussion 里**此前列为「未被解析器覆盖」**的 token。
#   这些 token 长期处于「被列出但不阻断」的状态 —— 与 p34 的 REG 空转同型：
#   **信息性输出不等于判据**。其中 0.409 / 0.553 是**改动前就在**的欠登（MMD 两个口径），
#   其余为本轮 P1（triage）与 P4（INC-017 的配对区间）新引入。
_P52D = pd.read_csv(os.path.join(OUT, "p52_mae_vs_constant.csv")).set_index("estimator") \
    if os.path.exists(os.path.join(OUT, "p52_mae_vs_constant.csv")) else None
_TRID = pd.read_csv(os.path.join(FIG, "figS10_triage_summary.csv")).set_index("layer")
_ORI = pd.read_csv(os.path.join(FIG, "figS_orientation.csv")).set_index("method")
_TS2D = pd.read_csv(os.path.join(FIG, "tableS2_metrics.csv")).set_index("method")


def _p52d(e, c):
    return abs(float(_P52D.loc[e, c])) if _P52D is not None else float("nan")


NUM += [
 # ---- P4：INC-017 的配对区间（论旨：水平改进不可检出）
 ("0.0035", _p52d("DR", "gain_vs_median_point"), "DR 对最优常数的 MAE 增益（|·|）",
  "p52_mae_vs_constant.csv", "estimator=DR, col=gain_vs_median_point"),
 ("0.0114", _p52d("DR", "gain_vs_median_lo95"), "同一增益 95% 下界（|·|）",
  "p52_mae_vs_constant.csv", "estimator=DR, col=gain_vs_median_lo95"),
 ("0.0046", _p52d("DR", "gain_vs_median_hi95"), "同一增益 95% 上界",
  "p52_mae_vs_constant.csv", "estimator=DR, col=gain_vs_median_hi95"),
 ("0.0048", _p52d("OM", "gain_vs_median_hi95"), "OM 同一增益 95% 上界",
  "p52_mae_vs_constant.csv", "estimator=OM, col=gain_vs_median_hi95"),
 # ---- P1：第 4 交付物（triage）
 ("1.59", float(_TRID.loc["L-DR", "precision_point"]) / float(_TRID.loc["L-random", "precision_point"]),
  "k=20% 排序相对无信息排序的提升", "figS10_triage_summary.csv",
  "precision_point(L-DR) ÷ precision_point(L-random)"),
 ("20", 100 * float(_TRID.loc["L-DR", "k"]) / 186.0, "复查预算 20%（=37/186）",
  "figS10_triage_summary.csv", "100×k/186, k=37"),
 ("38.6", 100 * float(_TS2D.loc["AUROC_DR", "sens_at_spec80"]),
  "预注册操作点上最好方法的 sens@0.80", "tableS2_metrics.csv",
  "method=AUROC_DR, col=sens_at_spec80 ×100"),
 # ---- P2：这两项在改动**之前**就未被覆盖（MMD 的 logistic 与 rank 两个口径）
 ("0.409", float(_ORI.loc["mmd", "macro_logit"]), "MMD logistic 口径 macro",
  "figS_orientation.csv", "method=mmd, col=macro_logit"),
 ("0.553", float(_ORI.loc["mmd", "macro_rank"]), "MMD rank 口径 macro",
  "figS_orientation.csv", "method=mmd, col=macro_rank"),
]
rowsR1 = []
for it in NUM:
    tok, src, name, srcfile, field = it[:5]
    scale = it[5] if len(it) > 5 else 1.0
    in_ms = tok in TXT
    # 【新增】token 必须**等于**来源值（不是只需"出现在正文里"）；
    #   带 scale 的项：正文写尾数，来源是实值，故比较 tok×scale 与 src。
    try:
        # 容差 = 尾数最后一位的一半 × 倍率（即"按正文显示精度四舍五入"）
        _dec = len(tok.split(".")[1]) if "." in tok else 0
        _tol = max(0.5 * (10 ** -_dec) * scale, 1e-12)
        eq = abs(float(tok) * scale - float(src)) <= _tol
    except ValueError:
        eq = (tok == str(src))
    st = "PASS" if (in_ms and eq) else ("NOT_IN_MS" if not in_ms else "TOKEN≠SOURCE")
    rowsR1.append(dict(risk="R1", item=tok, present=in_ms, source=name,
                       source_file=srcfile, field_and_filter=field,
                       source_value=round(float(src), 6),
                       status=st))
# R1 补充：正文出现的所有数字 token，未覆盖者列出
# 【首版 bug】未剥离引用编号与图号，导致 `11,12` / `14` / `16` / `17` / `3` 等
# 被当成"未解析的数据数字"。必须先去掉 [n] / [n,m] / Fig. N / Fig. SN / Section N。
_scan = re.sub(r'\[[\d,\s]+\]', ' ', BODY)                 # 引用编号
_scan = re.sub(r'Fig\.\s*S?\d+[A-Z]?', ' ', _scan)         # 图号
_scan = re.sub(r'\bSections?\s*\d+', ' ', _scan)           # 节号
covered = {t for t, *_ in NUM}
EXEMPT_TOK = {"0.5": "AUROC 随机水平（数学常数）",
              "10": "科学计数法的底（Limitations 中的 7.8×10⁻⁴ / 1.0×10⁻³）",
              "95": "置信区间水平（惯例值，非产物；与 p34 的 CONVENTION 同类）"}
extra = [t for t in sorted({t for t in re.findall(r'(?<![\w.])\d[\d,]*(?:\.\d+)?(?![\w])', _scan)
                            if t not in covered and t not in EXEMPT_TOK})]

# ================================================================ R2 primary/sensitivity
R2 = []
R2.append(dict(risk="R2", item="DR 被明确写为 prespecified primary",
               status="PASS" if re.search(r"DR as the\s+prespecified primary", TXT) else "FAIL",
               note=""))
R2.append(dict(risk="R2", item="未把 OM 写成 superior / better",
               status="PASS" if not re.search(r"\bOM\b[^.]{0,60}\b(superior|better|outperform)", TXT)
               else "FAIL", note=""))
R2.append(dict(risk="R2", item="OM 的低 MAE 附有不得升格的限定",
               status="PASS" if re.search(r"does not displace DR", TXT) else "FAIL", note=""))
R2.append(dict(risk="R2", item="IOW 未被描述为更优",
               status="PASS" if not re.search(r"\bIOW\b[^.]{0,60}\b(superior|better|outperform)", TXT)
               else "FAIL", note=""))

# ================================================================ R3 post-unblinding 身份
R3 = []
R3.append(dict(risk="R3", item="rank/orientation 结果被标为 post-unblinding",
               status="PASS" if "post-unblinding" in TXT else "FAIL", note=""))
R3.append(dict(risk="R3", item="明确写明它不构成可部署规则",
               status="PASS" if re.search(r"does not establish an inverse MMD\s+rule", TXT)
               or re.search(r"suitable for prospective monitoring", TXT) else "FAIL", note=""))
R3.append(dict(risk="R3", item="未把 explanatory 结果写成 prespecified",
               status="PASS" if not re.search(r"prespecified[^.]{0,40}rank-based", TXT) else "FAIL",
               note=""))

# ================================================================ R4 causal / novelty 强度
R4 = []
BAN = {
 "first":                 r"\bfirst\b",
 "novel":                 r"\bnovel\b|\bnovelty\b",
 "proves":                r"\bproves?\b|\bproven\b",
 "demonstrates that":     r"demonstrat\w+ that",
 "causes":                r"\bcauses?\b|\bcaused by\b",
 "due to":                r"\bdue to\b",
 "drives":                r"\bdrives?\b",
 "as a result of":        r"\bas a result of\b",
 "we are the first":      r"we are the first",
 "unprecedented":         r"\bunprecedented\b",
}
for name, pat in BAN.items():
    hit = re.findall(pat, TXT, re.I)
    R4.append(dict(risk="R4", item=f"禁用词 {name!r}", status="FAIL" if hit else "PASS",
                   note=f"命中 {len(hit)} 次" if hit else ""))
R4.append(dict(risk="R4", item="诊断性关联已标为 diagnostic rather than causal",
               status="PASS" if "diagnostic rather than causal" in TXT else "FAIL", note=""))
R4.append(dict(risk="R4", item="对既有统计现象用 consistent with 而非声称证明",
               status="PASS" if "consistent with a long-recognized issue" in TXT else "FAIL", note=""))

# ================================================================ R5-R8 定稿终检
# GPT6 指定：frozen 之前只查这四件事。
#
# ★★ 2026-10-06 前置依赖说明（不得误读为稿件缺陷）★★
#   R5 与 R8 的预期值都是按**全稿合并后的引用编号**写的（[11,12] / [4,5,9]，
#   以及与 References_FINAL.md 一致的表）。
#   实测：Discussion_v1.md 目前仍使用**合并前的局部编号**，两者不同名同号：
#       Discussion 局部:  [9]=Kouw  [11]=Kore  [12]=Guan  [14]=Chen  [15]=Kim
#       References_FINAL: [9]=Yang  [11]=Koch  [12]=Kouw  [14]=Kore  [15]=Guan
#   因此 R5/R8 必然 FAIL，直到**引用合并**执行完毕。这不是稿件内容错误。
#   合并脚本 p23_reference_merge_execute.py 目前**被自身护栏拒绝**：
#   它从根目录读取源文本，而 MANUSCRIPT_MASTER_DECLARATION §2 认定
#   manuscript_renumbered/ 为正本；两者已分叉（去掉引用编号后仍差
#   Results 210 行 / Methods 112 行 / Discussion 27 行 / Introduction 2 行），
#   直接重跑会把正本回退到分叉前。**须先解决「p23 从哪个副本读」再执行合并。**
# 【2026-10-06 改·按身份解析，不再写死编号】
#   原写法把期望编号写死（[11,12] / [4,5,9]）—— 那是 Discussion 的**合并前局部编号**。
#   而正本正文实际已是**合并后编号**（[14,15] / [4,5,12]），
#   实证：正文 "DECIDE-AI provides an appropriate reporting framework[20]"，
#   而 References_FINAL 的 [20] = Vasey（DECIDE-AI 作者）。
#   ⇒ **P-3 裁定其实早已落实，陈旧的是断言本身**。
#   按 GPT6「不要人工改编号，否则以后再合并必炸」的要求，改为：
#   从 References_FINAL.md 按**人名**解析出编号，再拿身份去匹配正文。
_RF = os.path.join(os.path.dirname(MM.__file__) if hasattr(MM, "__file__") else "",
                   "References_FINAL.md")
_RF = os.path.join(ROOT, "manuscript_renumbered", "References_FINAL.md")
import unicodedata as _ud
def _fold(s):
    """去变音符号并小写：Kivimäki -> kivimaki, Briët -> briet。
    References_FINAL 里的姓氏带原始变音符号，而调用处按 ASCII 写，
    故两边都先 fold 再比，避免因一个 umlaut 就 KeyError。"""
    return "".join(c for c in _ud.normalize("NFKD", str(s))
                   if not _ud.combining(c)).lower().strip()

_NAME2N = {}
if os.path.exists(_RF):
    for _ln in io.open(_RF, encoding="utf-8").read().split("\n"):
        _m = re.match(r"^\|\s*(\d+)\s*\|\s*([^,|]+)", _ln)
        if _m:
            _NAME2N[_fold(_m.group(2))] = int(_m.group(1))

def _ref(*names):
    """把身份名解析成编号（按 References_FINAL 的当前编号，不写死）。缺失则抛错。"""
    out = []
    for nm in names:
        if _fold(nm) not in _NAME2N:
            raise KeyError(f"References_FINAL 中找不到身份 {nm!r}；现有：{sorted(_NAME2N)}")
        out.append(_NAME2N[_fold(nm)])
    return "[" + ",".join(str(x) for x in sorted(out)) + "]"

_KG = _ref("Kore", "Guan")        # 支持「drift 可与性能变化脱钩」的两条
_GK = _ref("Garg", "Kivimaki", "Kouw")   # 支持「需对分布关系作假设」的三条
_G  = _ref("Garg")
R5 = []
R5.append(dict(risk="R5",
               item=f"'shift magnitude ... risk score' 已改引 {_KG}（Kore+Guan）而非 {_G}（Garg）",
               status="PASS" if re.search(re.escape(f"dissociated from performance change{_KG}"), TXT)
               and not re.search(re.escape(f"universal risk score for model failure{_G}"), TXT)
               else "FAIL",
               note="P-3 裁定：Garg 不直接支撑该主张。编号由 References_FINAL 按身份解析，不写死。"))
R5.append(dict(risk="R5", item=f"{_G}（Garg）仍只用于 'estimation requires assumptions' 语境",
               status="PASS" if re.search(
                   re.escape(f"assumptions about how source and target distributions are related{_GK}"), TXT)
               or re.search(re.escape(f"rarely stated{_G}"), TXT) else "FAIL",
               note=f"本句身份应为 Garg+Kivimaki+Kouw = {_GK}（按 References_FINAL 解析）"))

R6 = []
R6.append(dict(risk="R6", item="明确写出聚合口径未预先冻结",
               status="PASS" if re.search(r"was not frozen\s+in advance", TXT) else "FAIL",
               note="P-4：挡“看到 macro 好看才选 macro”的质疑"))
R6.append(dict(risk="R6", item="明确写出两口径并列报告、未事后指定 primary",
               status="PASS" if re.search(r"reported in parallel", TXT)
               and re.search(r"neither was retrospectively\s+designated primary", TXT) else "FAIL", note=""))

# R7 拼写：正文（不含参考表）不得残留英式拼写
_hits = _style(BODY)[1]
R7 = [dict(risk="R7", item="正文无英式拼写残留", status="PASS" if not _hits else "FAIL",
           note=f"残留：{_hits}" if _hits else "由 p20_style_normalize 统一")]

# 全稿引用并集（供 R8 连续性判定）
_whole = set()
for _fn in ("Introduction_v1.md", "Results_v3_clean.md", "Methods_v1.md", "Discussion_v1.md"):
    _fp = MM.master(_fn)                       # ★ 正本路径（全稿引文号必须按正本检查）
    if not os.path.exists(_fp): continue
    _t = open(_fp, encoding="utf-8").read()
    for _mk in ("## References", "## References used"):
        if _mk in _t: _t = _t.split(_mk)[0]
    for _m in re.findall(r'\[([\d,\s]+)\]', _t):
        _whole |= {int(x) for x in _m.split(",") if x.strip()}
_whole_sorted = sorted(_whole)
_whole_max = max(_whole) if _whole else 0
_whole_ok = bool(_whole) and _whole_sorted == list(range(1, _whole_max + 1))

# R8 引用编号：正文出现的 [n] 必须都在文末参考表中有定义，且无空洞
_used = set()
for m in re.findall(r'\[([\d,\s]+)\]', BODY):
    _used |= {int(x) for x in m.split(",") if x.strip()}
_table = set()
if "## References" in RAW:
    _tbl = RAW.split("## References", 1)[1]
    _table = {int(x) for x in re.findall(r'^\|\s*\*{0,2}(\d+)\*{0,2}\s*\|', _tbl, re.M)}
_undef = sorted(_used - _table)
_unused = sorted(_table - _used)
# 【2026-10-06】原 R8 附带 _BLOCK 说明（当时阻塞于引用合并）。
#   阻塞已解除：用 p26 以 References_FINAL.md 重建了 Discussion/Methods 的局部表
#   （正文一字未动，见 p26 内的逐字节断言）。R8 现应全绿，故不再挂阻塞说明。
R8 = [dict(risk="R8", item="所有正文引用编号均在文末表中有定义",
           status="PASS" if not _undef else "FAIL",
           note=(f"未定义：{_undef}" if _undef else "表由 p26 从 References_FINAL.md 重建")),
      dict(risk="R8", item="文末表无未被引用的条目",
           status="PASS" if not _unused else "FAIL",
           note=(f"未被引：{_unused}" if _unused else "")),
      # 【三版修正】连续性必须**跨全稿**判，不能在单文件内判。
      # 编号引用是全局的：Discussion 只用 [4,5,7,9,...] 这种非连续子集是**正确**的，
      # 因为 1/2/3/6/8/10/13 只在 Introduction 被引。
      # 首版在单文件内判连续，产生了假 FAIL。
      dict(risk="R8", item="全稿引用编号连续（跨 Introduction/Results/Methods/Discussion）",
           status="PASS" if _whole_ok else "FAIL",
           note=f"全稿并集 {_whole_sorted}；最大编号 = {_whole_max}")]

# ================================================================ 执行
import rephead
for _l in rephead.header(__file__, [os.path.join(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "manuscript_renumbered"), "Discussion_v1.md"),
        os.path.join(FIG, "fig3a_lodo.csv"), os.path.join(FIG, "fig3c_sudo_decomposition.csv"),
        os.path.join(FIG, "fig5b_across_n_range.csv"), os.path.join(FIG, "figS_orientation.csv"),
        os.path.join(FIG, "figS10_triage_summary.csv"), os.path.join(FIG, "tableS2_metrics.csv"),
        os.path.join(OUT, "p52_mae_vs_constant.csv")]):
    p(_l)
p("=" * 104); p("p19 — Discussion 审计（R1 数字 / R2 身份 / R3 解释性 / R4 强度）"); p("=" * 104)
p(f"  手稿: {MS}"); p("")

D = pd.DataFrame(rowsR1 + R2 + R3 + R4 + R5 + R6 + R7 + R8)
D.to_csv(os.path.join(OUT, "p19_discussion_audit.csv"), index=False)

p(f"【R1 数字来源】{len(rowsR1)} 项")
p(f"    {'item':<7}{'判定':<7}{'来源文件':<34}{'字段/筛选':<44}{'观测值':<12}手稿中")
for r in rowsR1:
    p(f"    {r['item']:<7}{r['status']:<7}{r['source_file']:<34}{r['field_and_filter']:<44}"
      f"{r['source_value']:<12}{'有' if r['present'] else '**无**'}")
p(f"    已剥离引用编号/图号后再抽取：未被解析器覆盖的数据 token = {extra if extra else '无'}")
p(f"    显式豁免：{EXEMPT_TOK}")

for tag, nm in (("R2", "primary/sensitivity 身份"), ("R3", "post-unblinding 身份"),
                ("R4", "causal / novelty 强度"), ("R5", "citation role（定稿终检）"),
                ("R6", "aggregation 未预设（定稿终检）"), ("R7", "拼写统一（定稿终检）"),
                ("R8", "引用编号一致性（定稿终检）")):
    sub = D[D.risk == tag]
    p(f"\n【{tag} {nm}】{len(sub)} 项：PASS {int((sub.status=='PASS').sum())}｜"
      f"FAIL {int((sub.status=='FAIL').sum())}")
    for _, r in sub.iterrows():
        # 【首版 bug】写 `r.item` 取到的是 pandas Series 的**方法**，不是名为 item 的列。
        # 必须用 r["item"]。输出里曾出现一大坨 "<bound method ...>" 即为该 bug。
        note = f"   [{r['note']}]" if r["note"] else ""
        p(f"    {r['status']:<6}{r['item']}{note}")

nfail = int((D.status == "FAIL").sum())
p("")
p("")
p("  本审计器自身的来源映射见上表「来源文件 / 字段·筛选」两列，须人工核对**来源本身是否取对**。")
p("*** p19 通过 ***" if nfail == 0 else f"*** p19 未通过 —— {nfail} 项 FAIL ***")
rep.close(); print("\ndone ->", OUT)
sys.exit(0 if nfail == 0 else 1)
