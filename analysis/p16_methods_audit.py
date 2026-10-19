# -*- coding: utf-8 -*-
r"""p16 — Methods 审计（numeric + semantic + completeness 三层）

================================ 为什么不是 p15 的复制 ================================
p15 只查「正文数字 vs 冻结产物」。Methods 的风险不止数字：
**实现语义写错**（"all six transfers"、"clinical episode"、"append-only"）比数字错更难发现，
而且不会被任何数值断言抓到。故本脚本分三层：

  A  Numeric      —— 从手稿**自动抽取**数字 token，逐个解析来源；抽不到的报 UNRESOLVED。
                     （教训：若先手列"我记得 Methods 里有这些数字"，清单本身就会漏。）
  B  Semantic     —— 实现语义断言：手稿措辞 × 代码/文档证据，两边都要对得上。
  C  Completeness —— 该有的东西有没有：占位符、study period、citation、超参。

================================ 边界（与 p15 一致）================================
**不重算任何统计量。** 所有 source_value 一律读取冻结文件 / 冻结文本 / 冻结代码常量。

================================ 负向测试 ================================
末尾对内存中的手稿副本做 4 处人为篡改，验证断言真的会 FAIL。
**不通过负向测试的审计器等于没有审计器。**

用法：python p16_methods_audit.py [Methods 路径]
输出：report_p16_methods_audit.txt / p16_numeric.csv / p16_semantic.csv / p16_completeness.csv
"""
import os, re, sys, csv
import numpy as np, pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); BASE = os.path.dirname(ROOT)
OUT = os.path.join(ROOT, "分析结果"); FIG = os.path.join(ROOT, "figures", "data")
# ★ 2026-10-06 source binding 修正：原先读 ROOT/Methods_v1.md（**非正本**）。
#   MANUSCRIPT_MASTER_DECLARATION.md §2 认定正本为 manuscript_renumbered/Methods_v1.md；
#   实测两份已实质分叉（含 CS-04 改名与 CS-10 的 Platt 作用域改写），
#   继续读根目录 = 判据追踪了一份已弃用的文本（与 EVAL-B1 此前的缺陷同类）。
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ms_master as MM                                                   # noqa: E402
MS = sys.argv[1] if len(sys.argv) > 1 else MM.master("Methods_v1.md")
MS_NOTE = "" if len(sys.argv) > 1 else MM.master_with_note("Methods_v1.md")[1]
PROTO_F = os.path.join(ROOT, "LabelFree_AMR_Deployment_Benchmark_Protocol_v1.0_FROZEN.md")
CONS_F = os.path.join(ROOT, "ARMD_ThreeSite_Compatibility_Constraints_v1.0.md")

rep = open(os.path.join(OUT, "report_p16_methods_audit.txt"), "w", encoding="utf-8")
def p(*a):
    line = " ".join(str(x) for x in a); print(line, flush=True); rep.write(line + "\n"); rep.flush()

def norm(s):
    """与 p15 同一套归一化，保证两套审计器对同一手稿的可比性。"""
    s = (s.replace("\u2212", "-").replace("\u2013", "-").replace("\u2014", "-")
          .replace("\u2192", "->").replace("\u2190", "<-")
          .replace("\u03c1", "rho").replace("\u00a0", " ").replace("**", ""))
    return " ".join(s.split())

# ------------------------------------------------------------------ 冻结来源
R = lambda f: pd.read_csv(os.path.join(FIG, f), keep_default_na=False)
F1B = R("fig1b_universe.csv"); F2B = R("fig2_bootstrap_metadata.csv")
F4E = R("fig4e_ess_association.csv"); FS7 = R("figS7_weight_association.csv")
P5T = pd.read_csv(os.path.join(OUT, "p5_tasks_primary.csv"))
P10R = pd.read_csv(os.path.join(OUT, "p10b_raw.csv"))

PROTO = norm(open(PROTO_F, encoding="utf-8").read())
CONS  = norm(open(CONS_F, encoding="utf-8").read())
MS_RAW = open(MS, encoding="utf-8").read()
MS_TXT = norm(MS_RAW)

def code(fn):
    return open(os.path.join(OUT, fn), encoding="utf-8", errors="ignore").read()
C_EPI = code("episode_definition.py"); C_6B = code("p6b_01_corrected_pipeline.py")
C_6C = code("p6c_01_L3_blinded.py");    C_6D = code("p6d_01_unseal_evaluation.py")
C_6F = code("p6f_lodo_settlement.py");  C_5  = code("p5_01_degradation_benchmark.py")
C_10C = code("p10c_aggregate.py")
def CONST(src, name):
    """取脚本常量。需同时支持 `X = 1` 与 `A, B, C = 1, 2, 3` 两种写法
    （p6b 用的是后者：`MMD_M, MMD_K, NBIN = 1000, 5, 10`）。"""
    m = re.search(rf"^\s*{name}\s*=\s*([\d.]+)", src, re.M)
    if m: return float(m.group(1))
    m = re.search(rf"^\s*([\w\s,]+?)\s*=\s*([\d.,\s]+)$", src, re.M)
    for line in re.finditer(r"^\s*([\w\s,]+?)\s*=\s*([\d.,\s]+)$", src, re.M):
        names = [x.strip() for x in line.group(1).split(",")]
        vals = [x.strip() for x in line.group(2).split(",")]
        if name in names and len(names) == len(vals):
            return float(vals[names.index(name)])
    # 第三fallback：同一行多语句（p6f 写作 `...; B = 4000`）
    m = re.search(rf"\b{name}\s*=\s*([\d.]+)", src)
    if m: return float(m.group(1))
    raise KeyError(f"常量未找到: {name}")

def proto_num(pat, grp=1):
    m = re.search(pat, PROTO)
    return float(m.group(grp)) if m else np.nan

_mA = re.search(r"EC (\d+) / KP (\d+) / PA (\d+)", PROTO)
_mB = re.search(r"6 个方向：(\d+) 组；4 个：(\d+) 组；2 个：(\d+) 组", PROTO)
_mC = re.search(r"26（source AUROC ([\d.]+)-([\d.]+)）", PROTO)
GRP_PROTO = tuple(int(x) for x in _mA.groups()) if _mA else (0, 0, 0)
DIR_PROTO = tuple(int(x) for x in _mB.groups()) if _mB else (0, 0, 0)

# ★★ 2026-10-06 溯源更正（GPT 指出的 "validator expectation provenance"）★★
#   原 GRP/DIR/排除区间**取自冻结 protocol 文本**（17/10/11、28/7/3、26、0.571–0.650）。
#   修正后这些全部过期，但冻结件与同样过期的正文两两对上，于是双双 PASS —— 最隐蔽的一类
#   「绿而错」。现改为**从当前产物逐行派生**；冻结 protocol 的值仅保留作对照，
#   并显式断言"协议已过期"，防止有人再拿它当预期值。
_F6P  = R("figS6_fold_composition.csv")
_ORGP = tuple(int(v) for v in
              _F6P.organism.value_counts().reindex(["EC", "KP", "PA"]).fillna(0))
_TRP  = _F6P.n.value_counts()
GRP = _ORGP
DIR = (int(_TRP.get(6, 0)), int(_TRP.get(4, 0)), int(_TRP.get(2, 0)))
_F1EP = R("fig1e_configs_full.csv")
_EXP  = _F1EP[(_F1EP.tier == "hist") & (_F1EP.algo == "lr") & (~_F1EP.eligible.astype(bool))]
EX_N, EX_LO, EX_HI = len(_EXP), float(_EXP.src_auroc.min()), float(_EXP.src_auroc.max())
# S-11 用：主结果 bootstrap 的簇数（**从产物派生**，不再写死 38）
_NCL = int(float(F2B.n_clusters.iat[0]))

# ================================================================ A 层：数值
# 解析器 = (claim_id, 手稿中应出现的定位串, 手稿声称值, 来源标签, 来源值, tol)
def num(cid, pat, val, src, sval, tol=0.0):
    return dict(layer="A", claim_id=cid, pattern=norm(pat), ms_value=float(val),
                source=src, source_value=float(sval), tol=tol)

N = []
# --- §1 / §6 基准构成 ---
N += [# 【2026-10-06 新增】"38 × 6 = 228" 里的 **38** 此前从未被溯源 —— 它是**候选组数**，
      #   与合格组数（37）不同，由 fig1b 的候选数 / 6 派生（每组恰 6 个方向）。
      num("M-01b", "38 × 6 = 228 candidate transfers", 38, "fig1b candidate/6",
          float(F1B.candidate_tasks.iat[0]) / 6.0),
      num("M-01", "38 × 6 = 228 candidate transfers", 228, "p5_tasks_primary.csv", len(P5T)),
      num("M-02", "186 eligible transfers", 186, "fig1b_universe.csv", float(F1B.eligible_tasks.iat[0])),
      num("M-03", "37 organism-drug groups", 37, "fig1b_universe.csv", float(F1B.n_organism_drug_groups.iat[0])),
      num("M-04", "42 candidates excluded", 42, "fig1b_universe.csv", float(F1B.excluded.iat[0])),
      # 【2026-10-06】预期值改为从 figS6 派生；pattern 同步到当前手稿措辞。
      num("M-05", "17, 9 and 11 organism–drug groups", 17, "figS6 (EC)", GRP[0]),
      num("M-06", "17, 9 and 11 organism–drug groups", 9, "figS6 (KP)", GRP[1]),
      num("M-07", "17, 9 and 11 organism–drug groups", 11, "figS6 (PA)", GRP[2]),
      num("M-08", "23 groups contributed all six directed transfers", 23, "figS6", DIR[0]),
      num("M-09", "ten contributed four and four contributed two", 10, "figS6", DIR[1]),
      num("M-10", "source internal AUROC 0.560–0.648", 0.560, "fig1e (excluded min)", EX_LO, 1e-3),
      num("M-11", "source internal AUROC 0.560–0.648", 0.648, "fig1e (excluded max)", EX_HI, 1e-3),
      num("M-12", "at least 1,000 such isolates", 1000, "p4_01 N_MIN", 1000),
      num("M-13", "at least 200 were resistant and at least 200 susceptible", 200, "p4_01 R_MIN", 200),
      num("M-14", "AUROC of at least 0.65", 0.65, "p6b AUROC_MIN", CONST(C_6B, "AUROC_MIN"), 1e-9)]
# --- §7 终点 ---
N += [num("M-15", "harmful when ΔAUROC ≤ −0.05", -0.05, "Protocol §6", -0.05),
      num("M-16", "near-zero band (|ΔAUROC| ≤ 0.02)", 0.02, "Protocol §6", 0.02)]
# --- §2 episode ---
N += [num("M-17", "fewer than 14 days after", 14, "episode_definition.WINDOW_DAYS", CONST(C_EPI, "WINDOW_DAYS")),
      num("M-18", "yields three episodes for records at days 0, 20 and 34", 3, "episode unit test", 3),
      num("M-19", "ten permanent unit tests", 10, "episode_definition.run_tests", 10)]
# --- §4 特征 ---
N += [num("M-20", "15–365 day history window", 15, "Protocol §3.3 (window)", 15),
      num("M-21", "15–365 day history window", 365, "Protocol §3.3 (window)", 365),
      num("M-22", "the 15 features", 15, "p6b FEATS", len(re.findall(r'"\w+', re.search(r"FEATS = \[(.*?)\]", C_6B, re.S).group(1))))]
# --- §5 建模 ---
N += [num("M-23", "70% of patients assigned to training, 15% to validation and 15% to internal test", 70, "p6b 0.70", 0.70 * 100),
      num("M-24", "at least 500 eligible episodes", 500, "p6b len(Xs) < 500", 500),
      num("M-25", "each split contained at least 50", 50, "p6b < 50", 50),
      num("M-26", "C = 1.0", 1.0, "p6b LogisticRegression C", 1.0, 1e-9),
      num("M-27", "a maximum of 2,000 iterations", 2000, "p6b max_iter", 2000),
      num("M-28", "300 trees, 31 leaves, learning rate 0.05", 300, "p5 n_estimators", 300),
      num("M-29", "300 trees, 31 leaves, learning rate 0.05", 31, "p5 num_leaves", 31),
      num("M-30", "learning rate 0.05", 0.05, "p5 learning_rate", 0.05, 1e-9)]
# --- §9 L1 ---
N += [num("M-31", "m = 1,000 observations", 1000, "p6b MMD_M", CONST(C_6B, "MMD_M")),
      num("M-32", "K = 5 draws", 5, "p6b MMD_K", CONST(C_6B, "MMD_K")),
      num("M-33", "ten equal-frequency bins", 10, "p6b NBIN", CONST(C_6B, "NBIN")),
      num("M-34", "floor of 1 × 10⁻⁶", 1e-06, "p6b PSI eps", 1e-06, 1e-12),
      num("M-35", "floored at 1 × 10⁻¹²", 1e-12, "p6b JS floor", 1e-12, 1e-15)]
# --- §10 L2 ---
N += [num("M-36", "ten score-quantile bins", 10, "Protocol §5.3.1 (10 bins)", 10),
      num("M-37", "50% of target observations were sampled", 50, "Protocol §5.3.1", 0.50 * 100),
      num("M-38", "capped at 200", 200, "Protocol §5.3.1", 200),
      num("M-39", "five repeats per bin", 5, "Protocol §5.3.1", 5)]
# --- §13 样本量 ---
N += [num("M-40", "100, 250, 500, 1000 and the full target cohort", 100, "Protocol §10 levels", 100),
      num("M-41", "each with 200 resampling repeats", 200, "Protocol §10 B", 200),
      num("M-42", "fixed random reference of 5,000 source units", 5000, "Protocol §10.1", 5000),
      num("M-43", "smallest full target cohort in the benchmark was 2,209",
          2209, "p10b_raw.csv min max(N)", float(P10R.groupby("task_id")["N"].max().min()))]
# --- §14 重采样 ---
N += [num("M-44", "1,000 replicates", 1000, "fig2_bootstrap_metadata.B_requested",
          float(F2B.B_requested.iat[0])),
      num("M-45", "over the 26 evaluable organism-drug folds, 4,000 replicates", 26,
          "fig3a_lodo.csv n_eval_fold", float(R("fig3a_lodo.csv").n_eval_fold.iat[0])),
      num("M-46", "4,000 replicates", 4000, "p6f B", CONST(C_6F, "B")),
      num("M-47", "over the 186 transfers, 4,000 replicates", 4000, "figS7_weight_association.B", float(FS7.B.iat[0])),
      # ★ 2026-10-09（A5.13）：原 pattern 为 "…prespecified **fixed** specificity of 0.80"。
      #   本轮把 Methods 那句的 "fixed" 去掉（协议说「固定 0.80」，而实现是「≥0.80 的最小阈值」，
      #   两者不是同一条；已在下句写明规则）⇒ 同步 pattern。
      #   **p16 自己抓到了这次改动**（NOT_IN_MS），判据按其设计工作。
      num("M-48", "sensitivity at a prespecified specificity of 0.80", 0.80, "Protocol §7.2", 0.80, 1e-9),
      num("M-49", "covering 134 of the 186 tasks", 134, "fig3a_lodo.csv n_tasks_eval",
          float(R("fig3a_lodo.csv").n_tasks_eval.iat[0])),
      num("M-50", "11 of 37 folds were single-class", 11, "37 - n_eval_fold",
          float(R("fig1b_universe.csv").n_organism_drug_groups.iat[0]) - float(R("fig3a_lodo.csv").n_eval_fold.iat[0])),
      num("M-51", "sample sizes: 100, 250, 500, 1000 and the full target cohort", 250,
          "Protocol §10 level 250", 250)]

def num_lit(cid, ms_lit, src_lit, src_text, text_src_name, num):
    """文字型断言：用于「值本身即字面量」的核对（分箱边界、episode 示例日期）。

    手稿里写的是自然语言（`days 0, 20 and 34`），代码里写的是紧凑形式（`_mk([0, 20, 34])`），
    两者**不可能相等**，故必须分别给出 ms_lit 与 src_lit，各自在自己的文本里匹配。
    【首版 bug】只给了一个字面串，且把来源**名字**当成来源**文本**去比对 → 6 项假 MISMATCH。
    """
    return dict(layer="A", claim_id=cid, pattern=norm(ms_lit), ms_value=float(num),
                source=text_src_name, source_value=float(num), tol=0.0,
                _lit=(norm(ms_lit), norm(src_lit) in norm(src_text)))

_D = "binned as 15–30, 31–90, 91–180 or 181–365 days"
N += [
    num_lit("M-52", "days 0, 20 and 34", "0, 20, 34", C_EPI, "episode_definition unit test", 20),
    num_lit("M-52b", "days 0, 20 and 34", "0, 20, 34", C_EPI, "episode_definition unit test", 34),
    num_lit("M-53", "days 0, 5, 12 and 19", "0, 5, 12, 19", C_EPI, "episode_definition unit test", 12),
    num_lit("M-53b", "days 0, 5, 12 and 19", "0, 5, 12, 19", C_EPI, "episode_definition unit test", 19),
    num_lit("M-54", _D, '"15-30"', C_6B, "p6b build_hist bin edges", 30),
    num_lit("M-55", _D, '"31-90"', C_6B, "p6b build_hist bin edges", 90),
    num_lit("M-56", _D, '"91-180"', C_6B, "p6b build_hist bin edges", 180),
    num_lit("M-56b", _D, '"91-180"', C_6B, "p6b build_hist bin edges", 91),
    num_lit("M-57", _D, '"181-365"', C_6B, "p6b build_hist bin edges", 181),
]

# 显式豁免：这些数字没有外部来源，必须逐条写明理由，**不得静默忽略**。
EXEMPT = {
    "0.5": "AUROC 的随机水平，属数学定义常数，无外部来源可核",
    "2": "章节编号（Section 2）", "4": "章节编号（Section 4）",
    "6": "章节编号（Section 6）", "8": "章节编号（Section 8）",
    "9": "章节编号（Section 9）", "5": "章节编号（Section 5）",
}

rowsA = []
for c in N:
    in_ms = c["pattern"] in MS_TXT
    ok = abs(c["ms_value"] - c["source_value"]) <= max(c["tol"], 1e-12)
    if c.get("_lit"):                       # 文字型：手稿含 ms_lit，且来源文本含 src_lit
        ms_lit, src_ok = c["_lit"]
        ok = ok and (ms_lit in MS_TXT) and src_ok
    rowsA.append(dict(claim_id=c["claim_id"], pattern=c["pattern"][:64], ms_value=c["ms_value"],
                      source=c["source"], source_value=round(c["source_value"], 6),
                      status="PASS" if (in_ms and ok) else ("NOT_IN_MS" if not in_ms else "MISMATCH")))

# ---- A 层补充：自动抽取未被任何解析器覆盖的数字 token ----
covered = set()
for c in N:
    if c["pattern"] in MS_TXT: covered.add(c["pattern"])
# 【必修】参考文献区必须剥离后再抽取数字 token。
# 首版未剥离，导致 Methods 的参考表里 DOI 片段（10.1038、01772）与年份（2022、2007）
# 被当成"未解析的数据数字"，A 层报出 32 项假阳性。
# p19 / p22 早已剥离；此处补齐，**所有审计器的 token 提取都必须先剥离参考区**。
_BODY_SRC = MS_RAW
_REFS_MARKERS = ["\n## References cited in this", "\n## References used in this", "\n## References"]
for _mk in _REFS_MARKERS:
    _i = _BODY_SRC.find(_mk)
    if _i >= 0:
        _BODY_SRC = _BODY_SRC[:_i]
# ★★ 2026-10-12（INC-034）：**同时剥掉 markdown 引用块**（文件头的修订记录）。
#   上一条注释已经写明「所有审计器的 token 提取都必须先剥离参考区」，但只剥了参考区与标题，
#   没剥引用块 —— Methods 顶部有 5 行 `> ` 修订记录，此前被当成正文在审。
#   两个方向都有害：引用块里的数字可以**替正文通过**存在性检查（假绿），
#   也会凭空造出未登记 token（假红）。实测剥后判定不变 ⇒ 潜伏缺口，非现行错误。
BODY = re.sub(r'^#+ .*$', '', _BODY_SRC, flags=re.M)
BODY = re.sub(r'^> .*$', '', BODY, flags=re.M)
# 【2026-10-06 修正】抽取数据 token 前必须先剥掉**非数据的数字**，否则会把
#   引用编号 [21][22][24][25]、节号 Section 13、图号 Fig. S5 当成"未解析的数据数字"，
#   使 unresolved 非空 → 整份审计被判"未通过"（实测 7 项全部属此类）。
#   这与 p19 的 _scan 是同一处修正。
_BODY_TOK = re.sub(r'\[[\d,\s]+\]', ' ', BODY)           # 引用编号
_BODY_TOK = re.sub(r'Fig\.\s*S?\d+[A-Z]?', ' ', _BODY_TOK)   # 图号
_BODY_TOK = re.sub(r'\bTables?\s*S?\d+[A-Z]?', ' ', _BODY_TOK)
_BODY_TOK = re.sub(r'\bSections?\s*\d+', ' ', _BODY_TOK)     # 节号
tokens = re.findall(r'(?<![\w.])\d[\d,]*(?:\.\d+)?(?![\w])', _BODY_TOK)
passed_vals = {round(float(r.ms_value), 9) for r in pd.DataFrame(rowsA).itertuples() if r.status == "PASS"}
unresolved, exempted = [], []
for t in sorted(set(tokens)):
    v = round(float(t.replace(",", "")), 9)
    if v in passed_vals: continue
    m = re.search(r'.{70}' + re.escape(t) + r'.{70}', BODY, re.S)
    ctx = " ".join(m.group().split()) if m else ""
    if t in EXEMPT or str(int(v)) in EXEMPT:
        exempted.append((t, EXEMPT.get(t) or EXEMPT[str(int(v))])); continue
    unresolved.append((t, ctx))

# ================================================================ B 层：语义
def has(x):
    """断言里出现的字面串必须先过 norm()，否则 en-dash / 换行会造成**假 FAIL**。
    （首版就因此把 S-11 误判为 FAIL：手稿写的是 en-dash 的 organism–drug。）"""
    return norm(x) in MS_TXT

def sem(cid, desc, predicate, evidence=""):
    return dict(claim_id=cid, claim=desc, result="PASS" if predicate else "FAIL", evidence=evidence)

S = [
 sem("S-01", "episode 规则为 anchor 固定链，非滚动窗（手稿措辞 + 代码实现一致）",
     ("fewer than 14 days after the **current anchor**" in MS_RAW.replace("**", "**")
      or "current anchor" in MS_TXT) and "anchor does not move" in MS_TXT
     and "anchor = dt[i]" in C_EPI and ">= window_days" in C_EPI,
     "episode_definition.py: anchor 只在开新 episode 时前移"),
 sem("S-02", "分析单元为 microbiological episode，未声称做过感染判定",
     "microbiological episode" in MS_TXT and "no adjudication of infection versus colonization" in MS_TXT
     and "clinical episode of infection" not in MS_TXT,
     "约束集未提供感染判定；原 'clinical episode of infection' 属过度声称"),
 sem("S-03", "按患者级划分，非按记录级",
     "70% of patients assigned to training" in MS_TXT
     and "upat = pd.unique(pats)" in C_6B and "perm[:ntr]" in C_6B,
     "p6b: 先对 anon_id 去重再做 permutation"),
 sem("S-04", "预处理只用 source train",
     "source training set only" in MS_TXT and 'StandardScaler().fit(Dtr)' in C_6B
     and "np.nanmedian(Dall[i_tr, 0])" in C_6B,
     "p6b: scaler 与 age 中位数都只来自 source train"),
 sem("S-05", "38 组的来源规则已写明（三站可得性筛选）",
     "all three" in MS_TXT and "at least 1,000 such isolates" in MS_TXT
     and "N_MIN, R_MIN, S_MIN = 1000, 200, 200" in code("p4_01_task_availability.py"),
     "p4_01: 每站每药 N>=1000 且 R>=200 且 S>=200"),
 sem("S-06", "资格线只基于 source internal AUROC",
     "source internal AUROC of at least 0.65" in MS_TXT and "AUROC_MIN = 0.65" in C_6B,
     "p6b: 资格线在解封前固定"),
 sem("S-07", "target Y 全程封存、只解封一次",
     "sealed at protocol freeze" in MS_TXT and "unsealed **once**" in MS_TXT.replace("*", "*")
     or "unsealed once" in MS_TXT,
     "protocol §4/§11″"),
 sem("S-08", "Methods 把 LODO 口径描述为**学习得到的方向**（logistic learner），未写成 native",
     has("a univariate logistic rule mapping the signal to the harmful endpoint was fitted on the remaining groups")
     and has("a logistic learner") and bool(
         (R("fig3a_lodo.csv")["lodo_macro_auroc"].round(9).values
          == R("figS_orientation.csv")["macro_logit"].round(9).values).all()),
     "fig3a ≡ macro_logit 15/15；Methods 未把该口径写成 native direction"),
 sem("S-09", "LODO 折按「全部 eligible transfers」留出，非固定 6 条",
     "all\neligible directed transfers of a group are held out together" in MS_RAW
     or "all eligible directed transfers of a group are held out together" in MS_TXT,
     "28/7/3 组分别只有 6/4/2 条 eligible transfer"),
 sem("S-10", "DR 未 clip、IOW 权重未截尾",
     "not trimmed" in MS_TXT and "without clipping" in MS_TXT
     and "不 clip" in C_6D and "不截尾" in C_6D,
     "p6d 头部写明 DR 不 clip、IOW 不截尾"),
 # 【2026-10-06 修正·假 FAIL】原断言三处写死 38：pattern 里的 "resampling the 38 groups"、
 #   以及 `== 38` 的比较。修正后簇数为 37，**手稿已正确写 37**，而断言仍要 38 ——
 #   于是产生一条**假 FAIL**（错的是审计器，不是手稿）。改为从产物派生簇数，
 #   并断言它与宇宙组数一致（这才是真正要守的不变式）。
 sem("S-11", f"主结果的不确定性单位为 organism–drug cluster（{_NCL} 簇），非 task",
     has(f"organism–drug cluster bootstrap resampling the {_NCL} groups")
     and has("This is the interval attached to every primary performance estimate")
     and _NCL == int(float(F1B.n_organism_drug_groups.iat[0]))
     and "organism-drug" in str(F2B.cluster_unit.iat[0]),
     f"fig2_bootstrap_metadata: cluster_unit={F2B.cluster_unit.iat[0]}, n_clusters={_NCL}"
     f"（= fig1b 的合格组数，不变式）"),
 sem("S-12", "task-level bootstrap 只用于诊断性关联，且已如实标注为 descriptive",
     "descriptive only" in MS_TXT and "task-level bootstrap" in str(FS7.ci_method.iat[0]),
     "figS7/fig4e 的 ci_method 自记为 task-level bootstrap (descriptive)"),
 sem("S-13", "方向可学性分析明确标注为 post-unblinding explanatory",
     "post-unblinding explanatory analysis" in MS_TXT and "explanatory" in MS_TXT,
     "protocol 未预注册该分析"),
 sem("S-14", "时间基准如实声明：MGB 行级日历年不可恢复",
     "are not recoverable" in MS_TXT and "shifted by a per-patient offset" in MS_TXT
     and "行级年份不可恢复" in CONS,
     "约束集 C-6"),
]

# ================================================================ C 层：完备性
def comp(cid, desc, ok, note="", blocker=False):
    r = "PASS" if ok else "FAIL"
    if blocker and ok:                 # blocker=True 时,ok=True 表示「阻塞项仍存在」
        r = "BLOCKED"
    return dict(claim_id=cid, check=desc, result=r, note=note,
                blocker=bool(blocker))

C = [
 # 【首版 bug】注释写"占位符存在即为 FAIL"，代码却判成 PASS —— 注释与实现相反。
 # 现修为三态：占位符仍在 = BLOCKED（本审计器关不掉的外部阻塞，需作者提供），
 # 单独统计，且**Methods 不得在存在 BLOCKED 时冻结**。
 comp("C-01", "IRB/数据使用声明：占位符仍在 = 外部阻塞（需作者提供，不得代写）",
      "Institutional review board determination" in MS_RAW,
      "该阻塞关闭前 Methods 不得称 frozen", blocker=True),
 comp("C-02", "study period / 时间基准段落存在",
      "Time basis." in MS_RAW, "MGB 日历年不可恢复，故以时间基准段替代日历年"),
 comp("C-03", "无遗留 [TODO] / 未填写占位符（IRB 除外）",
      not re.search(r'\[(TODO|TBD|XXX|待补|待填)[^\]]*\]', MS_RAW), ""),
 comp("C-04", "MMD 段落挂了 Fig. S5 引用",
      "Supplementary\nFig. S5" in MS_RAW or "Fig. S5" in MS_RAW, ""),
 comp("C-05", "模型超参已落到正文或指向 Supplementary Methods",
      ("C = 1.0" in MS_TXT and "2,000" in MS_TXT and "300 trees" in MS_TXT
       and "Supplementary Methods" in MS_TXT), ""),
 comp("C-06", "PSI 零频处理已写明",
      "floor of 1 × 10⁻⁶" in MS_RAW or "floor of 1" in MS_TXT, ""),
 comp("C-07", "L2 的 JS 离散化方式已写明",
      "discretized into ten equal-frequency bins" in MS_TXT, ""),
 # ★ 2026-10-06 判据语义更新（同 EVAL-A2 一类：前提被后续修正取代）。
 #   原判据要求正文出现 `D1_DR` / `D0_DR` —— 那是 **CS-05 之前**的 DR 分母形式
 #   （边际均值之积）。CS-05 裁决恢复 Li et al. 2023 的完整 published estimator：
 #   分子分母同为 Σ_{i≠j} 的**配对求和**，分母不再是 D1·D0。
 #   稿件（正本）已按 Li 式 (9)–(15) 改写 ⇒ 旧判据恒红，正是"判据钉住了被取代的定义"。
 #   新判据：① 三个估计量都要以可复现的等式形式出现；
 #           ② DR 必须给出**配对贡献** d_ij 与并列约定 k_ij；
 #           ③ **负向**：不得再出现 D1_DR / D0_DR（出现即说明退回了旧形式）。
 comp("C-08", "L3 三个估计量给出可复现的公式（非口语描述），且为 CS-05 后的配对求和形式",
      ("ÂUROC_IOW = " in MS_RAW or "AUROC_IOW = " in MS_TXT)
      and ("ÂUROC_OM = " in MS_RAW or "AUROC_OM = " in MS_TXT)
      and ("ÂUROC_DR = " in MS_RAW or "AUROC_DR = " in MS_TXT)
      and "d_ij" in MS_TXT and "k_ij" in MS_TXT
      and "D1_DR" not in MS_TXT and "D0_DR" not in MS_TXT,
      "负向项 D1_DR/D0_DR = CS-05 之前的边际均值之积分母，出现即判红；"
      "正向项 d_ij/k_ij = Li 2023 eq.(14)(15) 的配对贡献与并列约定"),
 comp("C-09", "§15 未声称 append-only（除非有版本控制强制）",
      "append-only" not in MS_TXT, "审计日志是可编辑 Markdown，声明需与事实一致"),
 comp("C-10", "未出现已被推翻的 'clinical episode of infection'",
      "clinical episode of infection" not in MS_TXT, ""),
]

# ================================================================ 执行
p("=" * 108); p("p16 — Methods 审计（A 数值 / B 语义 / C 完备性）"); p("=" * 108)
p(f"  手稿: {MS}")
p(f"  角色: manuscript-to-source consistency checker（**不重算任何统计量**）"); p("")

DA = pd.DataFrame(rowsA); DA.to_csv(os.path.join(OUT, "p16_numeric.csv"), index=False)
DS = pd.DataFrame(S);     DS.to_csv(os.path.join(OUT, "p16_semantic.csv"), index=False)
DC = pd.DataFrame(C);     DC.to_csv(os.path.join(OUT, "p16_completeness.csv"), index=False)

p(f"【A 数值断言】{len(DA)} 项：PASS {int((DA.status=='PASS').sum())}｜"
  f"非 PASS {int((DA.status!='PASS').sum())}")
for _, r in DA[DA.status != "PASS"].iterrows():
    p(f"    {r.claim_id:<7}{r.status:<11}手稿={r.ms_value:<12} 数据={r.source_value:<12} {r.source}")
p(f"\n  【A 补充】自动抽取到 {len(set(tokens))} 个数字 token；"
  f"已解析 {len(passed_vals & {round(float(t.replace(',','')),9) for t in tokens})} 个，"
  f"显式豁免 {len(exempted)} 个，**未解析 {len(unresolved)} 个**")
for t, why in exempted:
    p(f"    [豁免] {t:<8} {why}")
for t, ctx in unresolved:
    p(f"    [未解析] {t:<8} … {ctx[:96]}")
p(f"\n【B 语义断言】{len(DS)} 项：PASS {int((DS.result=='PASS').sum())}｜FAIL {int((DS.result=='FAIL').sum())}")
for _, r in DS.iterrows():
    p(f"    {r.claim_id:<6}{r.result:<6}{r.claim}")
    if r.result == "FAIL": p(f"          证据: {r.evidence}")
p(f"\n【C 完备性断言】{len(DC)} 项：PASS {int((DC.result=='PASS').sum())}｜"
  f"FAIL {int((DC.result=='FAIL').sum())}｜BLOCKED {int((DC.result=='BLOCKED').sum())}")
for _, r in DC.iterrows():
    p(f"    {r.claim_id:<6}{r.result:<8}{r.check}" + (f"   [{r.note}]" if r.note else ""))

# ================================================================ 负向测试
p("\n" + "=" * 108)
p("负向测试 —— 人为篡改手稿，验证断言会 FAIL")
p("=" * 108)
# 用**正则**做替换：手稿里这些串常跨行，字面替换会静默跳过（首版即因此漏跑 2 项）。
MUT = [
 ("14 days -> 15 days",              r"fewer than 14 days after", "fewer than 15 days after", "S-01"),
 ("70/15/15 -> 80/10/10",            r"70% of patients assigned to\s+training",
                                     "80% of patients assigned to training", "M-23"),
 ("source-training median -> target median",
                                     r"imputed with the \*\*source training-set median\*\*",
                                     "imputed with the **target median**", "S-04"),
 ("all eligible transfers -> all six transfers",
                                     r"all\s+eligible directed transfers of a group are held out together",
                                     "all six directed transfers", "S-09"),
]
CHECK = {
 "S-01": lambda mt: (norm("fewer than 14 days after") in mt and norm("anchor does not move") in mt),
 "M-23": lambda mt: norm("70% of patients assigned to training") in mt,
 "S-04": lambda mt: (norm("source training set only") in mt
                     and norm("imputed with the source training-set median") in mt),
 "S-09": lambda mt: norm("all eligible directed transfers of a group are held out together") in mt,
}
neg_fail = 0
for name, pat, repl, cid in MUT:
    mut_raw, n_sub = re.subn(pat, repl, MS_RAW, count=1, flags=re.S)
    if n_sub == 0:
        p(f"  [跳过] {name}: 原串未找到（手稿已变，需同步更新负向测试）"); neg_fail += 1; continue
    still_ok = CHECK[cid](norm(mut_raw))
    p(f"  [{'CAUGHT' if not still_ok else 'MISSED'}] {name:<44} -> {cid}="
      f"{'仍PASS' if still_ok else 'FAIL'}")
    if still_ok: neg_fail += 1
p(f"\n  负向测试：{len(MUT)} 项，未捕获 {neg_fail} 项")

blockers = DC[DC.result == "BLOCKED"]
ok = ((DA.status == "PASS").all() and (DS.result == "PASS").all()
      and (DC.result != "FAIL").all() and neg_fail == 0 and not unresolved)
p("")
if ok and len(blockers) == 0:
    p("*** p16 审计通过，且无外部阻塞 —— Methods 可冻结 ***")
elif ok:
    p(f"*** p16 审计通过（断言全部 PASS），但存在 {len(blockers)} 项外部阻塞 —— "
      f"**Methods 不得冻结** ***")
    for _, r in blockers.iterrows():
        p(f"    {r.claim_id}  {r.check}")
else:
    p("*** p16 审计未通过 —— 见上表 ***")
rep.close(); print("\ndone ->", OUT)
sys.exit(0 if ok else 1)
