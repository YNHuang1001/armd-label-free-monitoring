# -*- coding: utf-8 -*-
r"""p15 — Results 数值一致性审计（manuscript-to-source consistency checker）

================================ 角色与边界（写死）================================
本脚本**不是第二套 analysis pipeline**。它只做一件事：
    把 Results 手稿里出现的每个数值/定性断言，逐条对照 `figures/data/` 的冻结产物。

**禁止**在本脚本内重新计算任何统计量（不重算 AUROC、不重跑 bootstrap、不重算分位数）。
source_value 一律从已有冻结文件**读取**。

================================ 它能抓到什么 ================================
1. 手稿数值与冻结数据不一致（写错、抄错、改了数据没改文）
2. 手稿被编辑后某个数值被改动（pattern 匹配失败 → 提示需重新核对）
3. 定性断言与数据方向不符（如 calibration 偏差方向）

================================ 输出 ========================================
report_p15_results_numeric_audit.txt
p15_results_numeric_audit.csv          （数值断言）
p15_results_qualitative_audit.csv      （定性断言，与数值分开）

用法：python p15_results_numeric_audit.py [手稿路径]
"""
import os, sys, csv, re, io
import numpy as np, pandas as pd

BASE = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))); OUT = os.path.dirname(os.path.abspath(__file__))
FIG  = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "figures", "data")
# ★ 2026-10-06 source binding 修正：原先读根目录的 Results_v3_clean.md（**非正本**）。
#   正本见 MANUSCRIPT_MASTER_DECLARATION.md §2 = manuscript_renumbered/。
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ms_master as MM                                                   # noqa: E402
MS   = sys.argv[1] if len(sys.argv) > 1 else MM.master("Results_v3_clean.md")
MS_NOTE = "" if len(sys.argv) > 1 else MM.master_with_note("Results_v3_clean.md")[1]
PROTO_F = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "LabelFree_AMR_Deployment_Benchmark_Protocol_v1.0_FROZEN.md")

# ★★ 2026-10-14（INC-039）：**给了 argv[1] 就是负向测试，产物必须落 TEST 后缀。**
#   病历：负向测试用 argv[1] 指向一份改动副本，却**照样覆写了生产报告与生产 CSV** ——
#   于是磁盘上留下了「审计未通过」的**假红**，而真正的运行结果被冲掉。
#   这比假绿更阴：**一个坏掉的仪器会以「仪器在正常工作」的样子出现在证据面上**。
#   与 INC-037 同族（**测试/参数化动了生产状态**）。
_TEST = len(sys.argv) > 1
_SUF  = "_TEST" if _TEST else ""
if _TEST:
    print("*** 负向测试模式（argv[1]=%s）—— 产物写 %s 后缀，不碰生产报告 ***" % (MS, _SUF))

rep = open(os.path.join(OUT, "report_p15_results_numeric_audit%s.txt" % _SUF), "w", encoding="utf-8")
def _sb():
    """source binding 可见化：把「读的是哪一份」打进报告，避免判据追踪弃用副本而不自知。"""
    _l = "  手稿(正本): %s" % MS
    if MS_NOTE:
        _l += "\n              " + MS_NOTE
    print(_l, flush=True); rep.write(_l + "\n"); rep.flush()
def p(*a):
    line = " ".join(str(x) for x in a); print(line, flush=True); rep.write(line + "\n"); rep.flush()

def norm(s):
    """归一化，使 pattern 与手稿稳定匹配。

    必须同时处理五类差异，否则会产生大量**假阳性**（首版因此误报 8 项 NOT_IN_MS）：
      1. Unicode 负号/连字符 U+2212 / U+2013 / U+2014 -> "-"（本项目报告大量用 U+2212）
      2. 希腊字母 rho(U+03C1) -> "rho"
      3. 箭头 U+2192 -> "->"，U+2190 -> "<-"
      4. markdown 加粗标记 ** 去掉
      5. 空白：换行/多空格 -> 单空格（手稿存在跨行断句）
    """
    s = (s.replace("\u2212", "-").replace("\u2013", "-").replace("\u2014", "-")
          .replace("\u2192", "->").replace("\u2190", "<-")
          .replace("\u03c1", "rho").replace("\u00a0", " ")
          .replace("**", ""))
    return " ".join(s.split())

# ---------------------------------------------------------------- 读取冻结产物
R = lambda f: pd.read_csv(os.path.join(FIG, f), keep_default_na=False)
# ★★ 2026-10-06：新增 F1B —— harmful 百分比的**分母必须取自本表，不能写死**。
#   事故：S1-07 原写 `100*F1C.n_harmful.iat[0]/202`，202 是**修正前**的合格数。
#   重跑后 n_harmful=57、eligible=186，正确比例是 57/186 = **30.6%**，
#   而写死的分母算出 57/202 = **28.2%** —— 并被自动同步写进了稿件。
#   这个错误是外部评审（GPT）靠"这个数应该多大"的常识发现的，不是我这边的判据发现的。
F1B = R("fig1b_universe.csv")
F1C = R("fig1c_summary.csv"); F1D = R("fig1d_by_direction.csv"); F1E = R("fig1e_configs.csv")
F2  = R("fig2_forest.csv"); F3A = R("fig3a_lodo.csv"); F3B = R("fig3b_loo_organism.csv")
F3C = R("fig3c_sudo_decomposition.csv"); F4D = R("fig4d_mae.csv"); F4E = R("fig4e_ess_association.csv")
# ★ 2026-10-07 新增：方向约定敏感性（|·| vs 原值）与 §7.2 四项指标。
#   两份都产在 分析结果/（不是 figures/data/），故用 OUT 读，避免扰动图数据目录。
_p32 = os.path.join(OUT, "p32_direction_convention.csv")
DS = pd.read_csv(_p32).set_index("method") if os.path.exists(_p32) else None
_p33 = os.path.join(OUT, "p33_protocol72_metrics.csv")
P72 = pd.read_csv(_p33) if os.path.exists(_p33) else None

# ★★ 2026-10-13 新增：**样本量档的覆盖面**（p61 产出）。
#   起因：B/C 轮往 §4 加了两个数（「14–28% of estimates … outside [0,1]」、「9,814 cells at
#   n = 100」），而 **p15 当时没有任何断言覆盖它们，也没有「未被覆盖 token」的机制** ——
#   也就是说，那两个数进正文是**无声的**。这正是本会话反复出现的那一类缺口。
_covf = os.path.join(OUT, "p61_coverage_by_level.csv")
COV = pd.read_csv(_covf) if os.path.exists(_covf) else None
_RED = ["100", "250", "500", "1000"]


def _cov_frac(frame, lv):
    r = frame[frame.level.astype(str) == lv]
    return float(r.iloc[0]["oor"]) / float(r.iloc[0]["cells"])


def _cov_min_oor():
    return min(_cov_frac(COV, lv) for lv in _RED) if COV is not None else float("nan")


def _cov_max_oor():
    return max(_cov_frac(COV, lv) for lv in _RED) if COV is not None else float("nan")


def _cov_oor_at(lv):
    if COV is None:
        return float("nan")
    return float(COV[COV.level.astype(str) == lv].iloc[0]["oor"])
# ★ 2026-10-09 新增：§4 两处重定位带入的数字。
#   ① `p52_mae_vs_constant.csv` —— 常数基线的**配对区间**。此前全项目没有：p32 只产点估计，
#      而正文写着 "a small but consistent advantage"。「consistent」需要区间支撑，故补 p52。
#   ② `figS10_triage_summary.csv` —— triage 政策模拟。第 4 交付物此前在正文章节里
#      **无家可归**（Results 从不引 Supplementary Table S3）。
_p52f = os.path.join(OUT, "p52_mae_vs_constant.csv")
P52 = pd.read_csv(_p52f).set_index("estimator") if os.path.exists(_p52f) else None
def _p52(e, col):
    return float(P52.loc[e, col]) if P52 is not None else float("nan")
def _rd(fn):
    q = os.path.join(OUT, fn)
    return pd.read_csv(q) if os.path.exists(q) else None
CB = _rd("p32_constant_baselines.csv")          # 常数基线（零变化 / 中位数 / DR / IOW / OM）
CL = _rd("p33_signal_clusters.csv")             # |r| 阈值 → 独立组数
TP = _rd("p33_signal_top_pairs.csv")            # 最高的若干 |r| 对
def _cb(name, col):
    return float(CB.loc[CB.name == name, col].iat[0]) if CB is not None else float("nan")
def _cl(thr, col):
    return float(CL.loc[CL.threshold == thr, col].iat[0]) if CL is not None else float("nan")
def _tp(a, b, col="abs_r"):
    if TP is None: return float("nan")
    r = TP[((TP.signal_a == a) & (TP.signal_b == b)) | ((TP.signal_a == b) & (TP.signal_b == a))]
    return float(r[col].iat[0]) if len(r) else float("nan")
TPCI = _rd("p33_top_pair_ci.csv")
def _tpci(col):
    return float(TPCI[col].iat[0]) if TPCI is not None and len(TPCI) else float("nan")
def _p33m(m, col):
    if P72 is None: return float("nan")
    r = P72[(P72.convention == "primary") & (P72.method == m)]
    return float(r[col].iat[0]) if len(r) else float("nan")
# L1 / L2 家族在 LODO macro 上的均值（从 fig3a_lodo.csv 读，不重算）
_L1M = ["domain_auc", "mmd", "psi_mean", "psi_max", "js_mean", "js_max"]
_L2M = ["pred_mean_shift", "entropy_mean_shift", "entropy_js", "conf_dist_js",
        "pred_pos_prev_shift", "SUDO_AURCC"]
def _fam(ms):
    v = [float(F3A.set_index("method").loc[m, "lodo_macro_auroc"]) for m in ms]
    return sum(v) / len(v)
def _ds(m, col):
    """从 p32 取一列；缺失时返回 NaN 使断言失败而不是崩溃。"""
    return float(DS.loc[m, col]) if DS is not None else float("nan")
# L2 = prediction-space（与 tableS2_metrics.csv 的 role 列一致）
_L2SET = {"pred_mean_shift", "entropy_mean_shift", "entropy_js", "conf_dist_js",
          "pred_pos_prev_shift", "SUDO_AURCC"}
def _one(m):
    return m in _L2SET
def _sg(m, col):
    """从 p33 取 signed 口径的一列。"""
    if P72 is None: return float("nan")
    r = P72[(P72.convention == "signed") & (P72.method == m)]
    return float(r[col].iat[0]) if len(r) else float("nan")
F4O = R("fig4e_dr_out_of_range.csv"); F4L = R("fig4_l3_scatter.csv")
F5  = R("fig5_subsample.csv"); F5B = R("fig5b_across_n_range.csv"); FSO = R("figS_orientation.csv")
CN9 = R("figS9_orientation_consistency.csv")   # 两个 learner 的逐折方向一致性（S9 Panel B sidecar）
A3L = R("fig3a_lodo.csv")
ORI = R("figS_orientation.csv")
CAL = {k: R(f"fig4d_calibration_{k}.csv") for k in ("DR", "IOW", "OM")}
CAL = {k: v.astype({"pred_d_mean": float, "actual_d_mean": float}) for k, v in CAL.items()}

g2  = lambda m, c: float(F2.set_index("method").loc[m, c])
g3a = lambda m, c: float(F3A.set_index("method").loc[m, c])
g3b = lambda m, c: float(F3B.set_index("method").loc[m, c])
g4d = lambda m, c: float(F4D.set_index("estimator").loc[m, c])
g5b = lambda m, c: float(F5B.set_index("method").loc[m, c])
gso = lambda m, c: float(FSO.set_index("method").loc[m, c])
PA  = F5[F5.level == "PAIRED_full_minus_100"].set_index("method")
gpa = lambda m, c: float(PA.loc[m, c])

# ★★ 2026-10-12（INC-034）：**先剥掉 markdown 引用块**（文件头的修订记录）再审。
#   病历：本文件顶部有一段 `> ` 引用块，记录历次修订。它不是正文，但此前**被当成正文在审**：
#     ① 方向一（假绿）：引用块里出现的数字可以**替正文通过**「某 token 出现在稿件中」的检查
#        —— 正文若删掉该数字，本审计**不会发现**；
#     ② 方向二（假红）：引用块里的数字会凭空造出未登记的 token。
#   实测（本轮探针，见 INCIDENT_LOG INC-034）：剥掉后**判定不变** ⇒ 该缺口当前是**潜伏的**，
#   尚未造成错误结论。但它是 INC-016 族（判据在某个区域里其实什么都没查）的同一形态，故修掉。
#   `p19` 早就有这一行；本文件与 `p16` 此前没有。
MS_RAW = open(MS, encoding="utf-8").read()
MS_TXT = norm(re.sub(r'^> .*$', '', MS_RAW, flags=re.M))
PROTO  = norm(open(PROTO_F, encoding="utf-8").read())   # 冻结 protocol，作为 §3 benchmark 构成的来源

# ---- 从 protocol 冻结表中抽取 benchmark 构成（不是重算统计量，是读取冻结文本） ----
import re as _re
_mA = _re.search(r"EC (\d+) / KP (\d+) / PA (\d+)", PROTO)
_mB = _re.search(r"6 个方向：(\d+) 组；4 个：(\d+) 组；2 个：(\d+) 组", PROTO)
_mC = _re.search(r"26（source AUROC ([\d.]+)-([\d.]+)）", PROTO)
GRP = tuple(int(x) for x in _mA.groups()) if _mA else (0, 0, 0)
DIR = tuple(int(x) for x in _mB.groups()) if _mB else (0, 0, 0)

# ---------------------------------------------------------------- 数值断言
# 每条：section, claim_id, pattern(手稿中应出现的片段), value(手稿声称值),
#        source_file, source_value, tol
def num(section, cid, pattern, value, srcfile, srcval, tol=5e-4):
    return dict(kind="numeric", section=section, claim_id=cid, pattern=norm(pattern),
                manuscript_value=float(value), source_file=srcfile,
                source_value=float(srcval), tolerance=tol)

N = []
# ---- §1 ----
N += [num("§1","S1-01","186 met the prespecified source-performance criterion",186,"fig1b_universe.csv",float(F1B.eligible_tasks.iat[0]),0),
      num("§1","S1-02","median change in discrimination was -0.024",-0.024,"fig1c_summary.csv",F1C.p50.iat[0]),
      num("§1","S1-03","P10 -0.107",-0.107,"fig1c_summary.csv",F1C.p10.iat[0]),
      num("§1","S1-04","P90 +0.043",0.043,"fig1c_summary.csv",F1C.p90.iat[0]),
      num("§1","S1-05","P90 - P10 = 0.151",0.151,"fig1c_summary.csv",F1C.p90_minus_p10.iat[0]),
      num("§1","S1-06","57 transfers (30.6%) were harmful",57,"fig1c_summary.csv",F1C.n_harmful.iat[0],0),
      num("§1","S1-07","30.6%",30.6,"fig1c_summary.csv",100*F1C.n_harmful.iat[0]/float(F1B.eligible_tasks.iat[0]),0.05),
      num("§1","S1-08","A further 43 lay within a descriptive near-zero band",43,"fig1c_summary.csv",F1C.n_stable.iat[0],0),
      num("§1","S1-09","-0.051 for MGB->UTSW",-0.051,"fig1d_by_direction.csv",
          float(F1D.set_index("direction").loc["MGB->UT","50%"]),1e-3),
      num("§1","S1-10","Stanford->MGB was approximately neutral at +0.003",0.003,"fig1d_by_direction.csv",
          float(F1D.set_index("direction").loc["SU->MGB","50%"]),1e-3),
      num("§1","S1-11","57 versus 63 with history-augmented features",57,"fig1e_configs.csv",
          float(F1E.set_index("config").loc["P1 History+LR","n_harmful_aligned"]),0),
      num("§1","S1-12","57 versus 63 with history-augmented features",63,"fig1e_configs.csv",
          float(F1E.set_index("config").loc["S2 History+LGBM","n_harmful_aligned"]),0),
      num("§1","S1-13","71 versus 70 with core features",71,"fig1e_configs.csv",
          float(F1E.set_index("config").loc["S1 Core+LR","n_harmful_aligned"]),0),
      num("§1","S1-14","71 versus 70 with core features",70,"fig1e_configs.csv",
          float(F1E.set_index("config").loc["S3 Core+LGBM","n_harmful_aligned"]),0),
      num("§1","S1-15","(186, 174, 64, and 60 of 228, respectively)",186,"fig1e_configs.csv",
          float(F1E.set_index("config").loc["P1 History+LR","n_eligible_own"]),0),
      num("§1","S1-16","(186, 174, 64, and 60 of 228, respectively)",174,"fig1e_configs.csv",
          float(F1E.set_index("config").loc["S2 History+LGBM","n_eligible_own"]),0),
      num("§1","S1-17","(186, 174, 64, and 60 of 228, respectively)",64,"fig1e_configs.csv",
          float(F1E.set_index("config").loc["S1 Core+LR","n_eligible_own"]),0),
      num("§1","S1-18","(186, 174, 64, and 60 of 228, respectively)",60,"fig1e_configs.csv",
          float(F1E.set_index("config").loc["S3 Core+LGBM","n_eligible_own"]),0)]
# ---- §2 full universe ----
for tag, m, v, lo, hi in [("S2-01","domain_auc",0.478,0.368,0.591),("S2-02","mmd",0.463,0.355,0.571),
                          ("S2-03","psi_mean",0.47,0.356,0.588),("S2-04","psi_max",0.448,0.335,0.562),
                          ("S2-05","js_mean",0.466,0.354,0.582),("S2-06","js_max",0.439,0.329,0.554)]:
    N += [num("§2",tag+"-pt",f"{v:.3f}",v,"fig2_forest.csv",g2(m,"harmful_auc"),1e-3),
          num("§2",tag+"-lo",f"{lo:.3f}",lo,"fig2_forest.csv",g2(m,"auc_lo"),1e-3),
          num("§2",tag+"-hi",f"{hi:.3f}",hi,"fig2_forest.csv",g2(m,"auc_hi"),1e-3)]
# ---- §2 LODO ----
for tag, m, v, lo, hi in [("S2-07","mmd",0.409,0.294,0.519),("S2-08","psi_mean",0.557,0.424,0.69),
                          ("S2-09","psi_max",0.394,0.265,0.532),("S2-10","js_mean",0.586,0.462,0.71),
                          ("S2-11","js_max",0.513,0.387,0.637),("S2-12","domain_auc",0.561,0.435,0.679)]:
    N += [num("§2",tag+"-pt",f"{v:.3f}",v,"fig3a_lodo.csv",g3a(m,"lodo_macro_auroc"),1e-3),
          num("§2",tag+"-lo",f"{lo:.3f}",lo,"fig3a_lodo.csv",g3a(m,"lo_95"),1e-3),
          num("§2",tag+"-hi",f"{hi:.3f}",hi,"fig3a_lodo.csv",g3a(m,"hi_95"),1e-3)]
N += [num("§2","S2-13","26 of 37 organism-drug folds, representing 134 of 186 tasks",26,
          "fig3a_lodo.csv",float(F3A.n_eval_fold.iat[0]),0),
      num("§2","S2-14","representing 134 of 186 tasks",134,"fig3a_lodo.csv",float(F3A.n_tasks_eval.iat[0]),0),
      num("§2","S2-15","macro AUROC of 0.553 (0.44-0.67)",0.553,"figS_orientation.csv",
          gso("mmd","macro_rank"),1e-3),
      num("§2","S2-16","yielded a macro AUROC of 0.409",0.409,"figS_orientation.csv",
          gso("mmd","macro_logit"),1e-3),
      num("§2","S2-17","rho = +0.126",0.126,"fig2_forest.csv",g2("domain_auc","spearman"),1e-3),
      num("§2","S2-18","mean PSI (+0.128)",0.128,"fig2_forest.csv",g2("psi_mean","spearman"),1e-3),
      num("§2","S2-19","mean JS divergence (+0.132)",0.132,"fig2_forest.csv",g2("js_mean","spearman"),1e-3),
      # 【2026-10-06 重锚定】原 pattern "rho = +0.133" 已不存在：该句在修正轮被改写为
      #   "…mean JS divergence (+0.132) and MMD (+0.133)…"（去掉了与 MMD 的对比措辞）。
      # ★ 2026-10-14 pattern 同步（正当修改）：L1 那一句按外部核验要求重写，由
      #   「and MMD (+0.133)」改为「, MMD (+0.133)」（并列出全部六个）。**断言值 0.133 未变**。
      num("§2","S2-20",", MMD (+0.133)",0.133,"fig2_forest.csv",g2("mmd","spearman"),1e-3),
      # 【2026-10-01 新增】clean 版正文补入 S9 双 learner 对照的折数细节：
      # MMD 的 logistic learner 在 16/25 折选正方向（frac_logit_pos = 0.64），
      # 而 rank 规则 25/25 选负方向（frac_rank_pos = 0.00）。
      num("§2","S2-21","selecting the positive orientation in 25 of 26 evaluable folds",25,
          "figS9_orientation_consistency.csv",
          float(CN9.set_index("method").loc["mmd","frac_logit_pos"])*26,0)]
# ---- §3 ----
N += [num("§3","S3-01","harmful-detection AUROC of 0.587 (0.481-0.691)",0.587,"fig2_forest.csv",g2("entropy_js","harmful_auc"),1e-3),
      num("§3","S3-02","confidence-distribution JS divergence 0.585 (0.480-0.686)",0.585,"fig2_forest.csv",g2("conf_dist_js","harmful_auc"),1e-3),
      num("§3","S3-03","entropy mean shift reached 0.612 (0.487-0.733)",0.612,"fig3a_lodo.csv",g3a("entropy_mean_shift","lodo_macro_auroc"),1e-3),
      num("§3","S3-04","entropy JS divergence 0.608 (0.466-0.736)",0.608,"fig3a_lodo.csv",g3a("entropy_js","lodo_macro_auroc"),1e-3),
      num("§3","S3-05","confidence-distribution JS divergence 0.619 (0.481-0.744)",0.619,"fig3a_lodo.csv",g3a("conf_dist_js","lodo_macro_auroc"),1e-3),
      num("§3","S3-06","entropy JS 0.452",0.452,"fig3b_loo_organism.csv",g3b("entropy_js","loo_org_macro"),1e-3),
      num("§3","S3-07","predicted-prevalence shift 0.370",0.37,"fig3b_loo_organism.csv",g3b("pred_pos_prev_shift","loo_org_macro"),1e-3),
      num("§3","S3-08","confidence-distribution JS 0.442",0.442,"fig3b_loo_organism.csv",g3b("conf_dist_js","loo_org_macro"),1e-3),
      num("§3","S3-09","harmful-detection AUROC was 0.527 (0.413-0.629)",0.527,"fig2_forest.csv",g2("SUDO_AURCC","harmful_auc"),1e-3),
      # 【2026-10-01 重锚】原 locator "was 0.756" 所在句在 SUDO 口径更正时被删除（值未变）。
      # 新 locator 与 S3-13 同句，但**来源不同**：S3-10 查 fig3a_lodo 的 macro、
      # S3-13 查 fig3c 的 macro_138task。两源必须给出同一个 0.756，故这是一处跨来源一致性检查。
      num("§3","S3-10","changed the estimate to 0.325",0.325,"fig3a_lodo.csv",g3a("SUDO_AURCC","lodo_macro_auroc"),1e-3),
      num("§3","S3-11","yielded an AUROC of 0.332",0.332,"fig3c_sudo_decomposition.csv",float(F3C.pooled_202task.iat[0]),1e-3),
      num("§3","S3-12","to 0.408",0.408,"fig3c_sudo_decomposition.csv",float(F3C.pooled_138task.iat[0]),1e-3),
      num("§3","S3-13","changed the estimate to 0.325",0.325,"fig3c_sudo_decomposition.csv",float(F3C.macro_138task.iat[0]),1e-3),
      num("§3","S3-14","a +0.075 task-subset effect",0.075,"fig3c_sudo_decomposition.csv",
          float(F3C.pooled_138task.iat[0])-float(F3C.pooled_202task.iat[0]),1e-3),
      num("§3","S3-15","further -0.083 aggregation effect",-0.083,"fig3c_sudo_decomposition.csv",
          float(F3C.macro_138task.iat[0])-float(F3C.pooled_138task.iat[0]),1e-3),
      # 【2026-10-06 重锚定·含公式修正】原为 "The 68 tasks" 且公式写死 202-134。
      #   修正后合格宇宙为 186，单类折覆盖 52 个 task。公式改为**从 fig1b 的合格数派生**，
      #   不再写死旧总数（写死正是这条当初能"通过"却错的原因）。
      num("§3","S3-16","The 52 tasks",52,"fig3a_lodo.csv",
          float(F1B.eligible_tasks.iat[0]) - float(F3A.n_tasks_eval.iat[0]),0),
      # ---- 2026-10-07 新增：方向约定敏感性段（正文 §3 最后一段）----
      #   来源 p32_direction_convention.csv（p32 独立重算，PROD 口径在 pooled/macro/LOO 三级
      #   15/15 逐位复现冻结产物）。该段的作用是把约定的代价写进正文，而不是留在补充材料。
      num("§3","S3-17","mean prediction shift reaches a macro AUROC of 0.665",0.665,
          "p32_direction_convention.csv", _ds("pred_mean_shift","raw_macro"),1e-3),
      num("§3","S3-19","predicted-prevalence shift 0.732",0.732,"p32_direction_convention.csv",
          _ds("pred_pos_prev_shift","raw_macro"),1e-3),
      num("§3","S3-20","against 0.329 and 0.337 under the magnitude",0.329,"p32_direction_convention.csv",
          _ds("pred_mean_shift","prod_macro"),1e-3),
      num("§3","S3-21","against 0.329 and 0.337 under the magnitude",0.337,"p32_direction_convention.csv",
          _ds("pred_pos_prev_shift","prod_macro"),1e-3),
      # ★★ 2026-10-09（INC-022）：本行原绑 `p32_direction_convention.csv:prod_macro` = 0.6125，
      #   那是该量的 **1-ULP 偏高**的那个 double。**p53 已用精确算术证明正确舍入是 0.612**
      #   （26 个折的精确均值严格低于进位边界 0.6125）⇒ **正文写 0.612 是对的，错的是绑定**。
      #   改绑**正典侧** `fig3a_lodo.csv`（0.6124999999999999）—— 让消费者读本来就没有分歧的生产者。
      #   这就是外部核验所称「潜在第 17 对」的正确答案：不是登记它，是**消除它**。
      num("§3","S3-22","from 0.612 to 0.401",0.612,"fig3a_lodo.csv",
          g3a("entropy_mean_shift", "lodo_macro_auroc"),1e-3),
      num("§3","S3-23","from 0.612 to 0.401",0.401,"p32_direction_convention.csv",
          _ds("entropy_mean_shift","raw_macro"),1e-3),
      num("§3","S3-24","mean prediction shift from 0.372 to",0.372,"p32_direction_convention.csv",
          _ds("pred_mean_shift","prod_loo"),1e-3),
      num("§3","S3-25","mean prediction shift from 0.372 to",0.686,"p32_direction_convention.csv",
          _ds("pred_mean_shift","raw_loo"),1e-3),
      num("§3","S3-26","predicted-prevalence shift from 0.370 to 0.718",0.370,"p32_direction_convention.csv",
          _ds("pred_pos_prev_shift","prod_loo"),1e-3),
      num("§3","S3-27","predicted-prevalence shift from 0.370 to 0.718",0.718,"p32_direction_convention.csv",
          _ds("pred_pos_prev_shift","raw_loo"),1e-3),
      num("§3","S3-31","mean entropy shift falls from 0.509 to",0.509,"p32_direction_convention.csv",
          _ds("entropy_mean_shift","prod_loo"),1e-3),
      num("§3","S3-32","mean entropy shift falls from 0.509 to",0.345,"p32_direction_convention.csv",
          _ds("entropy_mean_shift","raw_loo"),1e-3),
      # ---- 2026-10-07 第二批：独立信息源 + 家族级反转 ----
      num("§3","S3-33","Pearson r = 0.993",0.993,"p33_signal_top_pairs.csv",
          _tp("entropy_js","conf_dist_js"),1e-3),
      num("§3","S3-33a","95% interval 0.972 to",0.972,"p33_top_pair_ci.csv",_tpci("ci_lo"),1e-3),
      num("§3","S3-33b","95% interval 0.972 to 0.997",0.997,"p33_top_pair_ci.csv",_tpci("ci_hi"),1e-3),
      # ★ 2026-10-16 pattern 同步（正当修改）：原句 `only eight groups survived at |r| ≥ 0.90`
      #   **未交代口径**（哪个相关系数）⇒ 改为 `Clustering the fifteen signals by |Pearson r|
      #   leaves eight groups at a threshold of 0.90`。**断言对象未变**（仍是
      #   `p33_signal_clusters.csv` 的 threshold=0.90 → n_groups = 8）。见 CORRECTIONS_20261016.md R-2。
      num("§3","S3-34","leaves eight groups at a threshold",8,"p33_signal_clusters.csv",
          _cl(0.9,"n_groups"),0),
      num("§3","S3-35","the macro estimates averaged 0.472",0.472,"fig3a_lodo.csv",
          _fam(_L2M),1e-3),
      num("§3","S3-36","below the input-space average of 0.503",0.503,"fig3a_lodo.csv",
          _fam(_L1M),1e-3),
      num("§3","S3-28","discriminate below chance (0.362, 0.457 and 0.342; Supplementary Table S2, panel C)",0.362,"p33_protocol72_metrics.csv",
          _sg("pred_mean_shift","auroc"),1e-3),
      num("§3","S3-29","discriminate below chance (0.362, 0.457 and 0.342; Supplementary Table S2, panel C)",0.457,"p33_protocol72_metrics.csv",
          _sg("entropy_mean_shift","auroc"),1e-3),
      num("§3","S3-30","discriminate below chance (0.362, 0.457 and 0.342; Supplementary Table S2, panel C)",0.342,"p33_protocol72_metrics.csv",
          _sg("pred_pos_prev_shift","auroc"),1e-3)]
# ★★ 2026-10-09（INC-019 的横向缺口的对称面，由 pointer-value 审计发现）：
#   上面三句（S3-28/29/30）把**未变换**值锚得很好 —— 判据**全绿**，因为值本身是对的。
#   但原句写「in the prespecified native-orientation task **(Fig. 2A)** the untransformed values
#   of all three discriminate below chance (0.362, …)」：**读者被指向 Fig 2A 去核这三个数**，
#   而 Fig 2A（唯一数据源 fig2_forest.csv）上这三个信号是 **0.548 / 0.542 / 0.545，在 0.5 以上**。
#   ⇒ **值对、指针错**：现有判据只查「写的值对不对」，查不出「读者被指到哪里去核」。
#   处置：句子改为**同时给出两种口径**，使读者的比对必然成功；本处补断言锚住「另一侧」三个数。
N += [num("§3","S3-37","places the same three above chance (0.548, 0.542 and 0.545)",0.548,"fig2_forest.csv",
          g2("pred_mean_shift","harmful_auc"),1e-3),
      num("§3","S3-38","places the same three above chance (0.548, 0.542 and 0.545)",0.542,"fig2_forest.csv",
          g2("entropy_mean_shift","harmful_auc"),1e-3),
      num("§3","S3-39","places the same three above chance (0.548, 0.542 and 0.545)",0.545,"fig2_forest.csv",
          g2("pred_pos_prev_shift","harmful_auc"),1e-3)]
# ---- §4 ----
for tag, k, v in [("S4-01","DR",0.79),("S4-02","IOW",0.615),("S4-03","OM",0.802)]:
    N.append(num("§4",tag,f"{v:.3f}\nfor {k}",v,"fig4d_mae.csv",g4d(k,"rho_pred_true"),1e-3))
for tag, k, v in [("S4-04","DR",0.0416),("S4-05","IOW",0.054),("S4-06","OM",0.0417)]:
    N.append(num("§4",tag,f"{v:.4f} for {k}",v,"fig4d_mae.csv",g4d(k,"mae"),1e-4))
# ★★ 2026-10-06：稿件 §4 原先**逐一点名 4 个 `AUROC_DR > 1` 的 task 及其数值**
#   （PA_TZP_SU2UT 1.106 / PA_AMK_MGB2SU 1.085 / PA_AMK_SU2UT 1.027 / PA_ATM_SU2MGB 1.002）。
#   CS-05（恢复 Li 2023 的配对求和分母）之后**越界数归零**（实测 AUROC_DR ∈ [0.5135, 0.9590]），
#   `fig4e_dr_out_of_range.csv` 变成空表 ⇒ 原先的 `.loc[t]` 直接 KeyError，**整个 p15 在打印前就崩**。
#
#   处置：不再假设那 4 行一定存在。若某 task 不在表里（= 该项声称已不成立），
#   产出一条 **manifest_value 与 source 不符**的断言让审计判红并说清原因 ——
#   而不是崩溃。判红是对的：稿件确实还写着这四句，必须删/改。
# 已退役的 claim（声称的对象已不存在）。**记录但不再计入判据** —— 见下方 §4 处的说明。
RETIRED = []
_F4O_IDX = F4O.set_index("task_id")["AUROC_DR"] if (len(F4O) and "task_id" in F4O.columns) else None
for tag, t, v in [("S4-07", "PA_TZP_SU2UT", 1.106), ("S4-08", "PA_AMK_MGB2SU", 1.085),
                  ("S4-09", "PA_AMK_SU2UT", 1.027), ("S4-10", "PA_ATM_SU2MGB", 1.002)]:
    if _F4O_IDX is not None and t in _F4O_IDX.index:
        N.append(num("§4", tag, f"{t} ({v:.3f})", v, "fig4e_dr_out_of_range.csv",
                     float(_F4O_IDX.loc[t]), 1e-3))
    else:
        # ★★ 2026-10-06 **退役**（不是判红）。
        #   这 4 条声称的是「稿件 §4 逐一点名了 4 个 AUROC_DR > 1 的 task」。
        #   CS-05 修正后越界数归零，稿件里那 4 句**已删除并改写**为
        #   "After correction, no DR estimate exceeded the theoretical AUROC range"。
        #   ⇒ 声称的对象已不存在，判据**无靶子**。如实登记为「已退役」并给出理由：
        #     既不留在 N 里制造恒红，也不静默删除（那会丢掉"曾经有这条声称"的记录）。
        RETIRED.append(dict(claim_id=tag, task=t, old_value=v,
                            reason="稿件所述现象已被 CS-05 消除（AUROC_DR 越界 4 -> 0）；"
                                   "对应句子已删除并改写，claim 无靶子"))
N += [num("§4","S4-11","rho = -0.249",-0.249,"fig4e_ess_association.csv",float(F4E.spearman_rho.iat[0]),1e-3),
      num("§4","S4-12","-0.376 to -0.117",-0.376,"fig4e_ess_association.csv",float(F4E.lo_95.iat[0]),1e-3),
      num("§4","S4-13","-0.376 to -0.117",-0.117,"fig4e_ess_association.csv",float(F4E.hi_95.iat[0]),1e-3),
      num("§4","S4-14","11 to 11,250, with a median of 1,082",11,"fig4e_ess_association.csv",float(F4E.ess_min.iat[0]),0.5),
      num("§4","S4-15","11 to 11,250",11250,"fig4e_ess_association.csv",float(F4E.ess_max.iat[0]),1),
      num("§4","S4-16","median of 1,082",1082,"fig4e_ess_association.csv",float(F4E.ess_median.iat[0]),1),
      num("§4","S4-17","mean observed change was only +0.020",0.02,"fig4d_calibration_DR.csv",float(CAL["DR"].actual_d_mean.iat[-1]),5e-4),
      num("§4","S4-18","ratios of 0.78 for DR",0.78,"fig4d_calibration_DR.csv",
          abs(float(CAL["DR"].actual_d_mean.iat[0]))/abs(float(CAL["DR"].pred_d_mean.iat[0])),5e-3),
      num("§4","S4-19","0.46 for IOW",0.46,"fig4d_calibration_IOW.csv",
          abs(float(CAL["IOW"].actual_d_mean.iat[0]))/abs(float(CAL["IOW"].pred_d_mean.iat[0])),5e-3),
      num("§4","S4-20","0.74 for OM",0.74,"fig4d_calibration_OM.csv",
          abs(float(CAL["OM"].actual_d_mean.iat[0]))/abs(float(CAL["OM"].pred_d_mean.iat[0])),5e-3),
      # 【2026-10-06 重锚定】IOW 最正箱的实测均值由 −0.0048 变为 +0.0047（符号翻转），
      #   正文相应由 "slight mean degradation (0.005)" 改为 "essentially no mean change (+0.005)"。
      num("§4","S4-21","no mean change (+0.005)",0.005,"fig4d_calibration_IOW.csv",float(CAL["IOW"].actual_d_mean.iat[-1]),5e-4)]
# ---- 2026-10-07 新增：常数基线段（正文 §4 末）----
#   ★ 常数基线有**两个**：预测零变化 0.0502、预测观测中位数 0.0451。
#     只报前者会把 DR/OM 的优势从 0.0035 放大成 0.0086。两条都必须锚定，
#     且 S4-25/S4-26 用**同一 pattern、不同来源**（DR / OM）构成跨来源一致性检查。
N += [num("§4","S4-22","mean absolute error of 0.0502",0.0502,"p32_constant_baselines.csv",
          _cb("constant_zero","mae"),5e-5),
      num("§4","S4-23","observed median change (-0.024)",-0.024,"p32_constant_baselines.csv",
          _cb("constant_zero","median_d_auroc"),5e-4),
      num("§4","S4-24","gives 0.0451",0.0451,"p32_constant_baselines.csv",
          _cb("constant_median","mae"),5e-5),
      num("§4","S4-25","improve by only about 0.0035",0.0035,"p32_constant_baselines.csv",
          _cb("DR","gain_vs_median_floor"),5e-5),
      num("§4","S4-26","improve by only about 0.0035",0.0035,"p32_constant_baselines.csv",
          _cb("OM","gain_vs_median_floor"),5e-5)]

# ---- ★ 2026-10-09：§4 重定位带入的配对区间（p52）与 triage 政策模拟（Table S3 / figS10）
#   病史：正文原写 "a small but consistent advantage"。**"consistent" 需要区间支撑**，
#   而全项目从未算过这个区间 —— 补算后**五个比较有五个含 0**（唯一不含 0 的是 IOW
#   比常数**更差**），该词因此被撤换。见 INC-017。
TRS = R("figS10_triage_summary.csv").set_index("layer")
def _tr(layer, col):
    return float(TRS.loc[layer, col])
N += [num("§4","S4-27","-0.0035 (95% interval -0.0114 to +0.0046)",-0.0114,"p52_mae_vs_constant.csv",
          _p52("DR","gain_vs_median_lo95"),5e-5),
      num("§4","S4-28","95% interval -0.0114 to +0.0046",0.0046,"p52_mae_vs_constant.csv",
          _p52("DR","gain_vs_median_hi95"),5e-5),
      num("§4","S4-29","(-0.0114 to +0.0048)",0.0048,"p52_mae_vs_constant.csv",
          _p52("OM","gain_vs_median_hi95"),5e-5),
      num("§4","S4-30","-0.0086 (-0.0174 to +0.0005)",-0.0086,"p52_mae_vs_constant.csv",
          _p52("DR","gain_vs_zero_point"),5e-5),
      num("§4","S4-31","-0.0086 (-0.0174 to +0.0005)",-0.0174,"p52_mae_vs_constant.csv",
          _p52("DR","gain_vs_zero_lo95"),5e-5),
      num("§4","S4-32","-0.0085 (-0.0174 to +0.0006)",-0.0085,"p52_mae_vs_constant.csv",
          _p52("OM","gain_vs_zero_point"),5e-5),
      num("§4","S4-33","-0.0085 (-0.0174 to +0.0006)",0.0006,"p52_mae_vs_constant.csv",
          _p52("OM","gain_vs_zero_hi95"),5e-5),
      num("§4","S4-34","worse than the stronger constant by +0.0089",0.0089,"p52_mae_vs_constant.csv",
          _p52("IOW","gain_vs_median_point"),5e-5),
      num("§4","S4-35","(+0.0016 to +0.0164)",0.0016,"p52_mae_vs_constant.csv",
          _p52("IOW","gain_vs_median_lo95"),5e-5),
      num("§4","S4-36","(+0.0016 to +0.0164)",0.0164,"p52_mae_vs_constant.csv",
          _p52("IOW","gain_vs_median_hi95"),5e-5),
      # ---- triage：只锚 **k = 20% 主格**（A5.4 约束），10% / 30% 两列不进正文
      num("§4","S4-37","the top 20% of transfers, 37 of 186",37,"figS10_triage_summary.csv",
          _tr("L-DR","k"),0),
      num("§4","S4-38","identified 18 of the 57 failing transfers",18,"figS10_triage_summary.csv",
          _tr("L-DR","caught"),0),
      num("§4","S4-39","1.59 times the 11.34 expected",1.59,"figS10_triage_summary.csv",
          _tr("L-DR","precision_point") / _tr("L-random","precision_point"),5e-3),
      num("§4","S4-40","1.59 times the 11.34 expected",11.34,"figS10_triage_summary.csv",
          _tr("L-random","caught"),5e-3),
      # ★ 2026-10-12（Phase A）：原 pattern 写「the input-shift composite identified 7」，
      #   而 `composite` 在 Methods 里**从未定义**（且 Methods 声明 single-signal only）
      #   ⇒ 外部核验指出后改为「a composite that rank-averages the six input-shift
      #   statistics」——**把定义带回正文**，同时消掉与 Methods 的张力。pattern 随之同步。
      num("§4","S4-41","that rank-averages the six input-shift statistics identified 7",7,"figS10_triage_summary.csv",
          _tr("L-shift","caught"),0),
      num("§4","S4-42","precision difference 0.297",0.297,"figS10_triage_summary.csv",
          _tr("L-DR_minus_L-shift","precision_point"),5e-4),
      num("§4","S4-43","0.079 to 0.486",0.079,"figS10_triage_summary.csv",
          _tr("L-DR_minus_L-shift","lo95"),1e-4),
      num("§4","S4-44","0.079 to 0.486",0.486,"figS10_triage_summary.csv",
          _tr("L-DR_minus_L-shift","hi95"),5e-4),
      # ★★ 2026-10-13（B/C 轮）：正文新加的两句里有两个数此前**没有任何判据在看着**
      #   —— 而 p15 又没有「未被覆盖的 token」机制，所以放进去是无声的。
      #   来源：p61_coverage_by_level.csv（由 p61_transfer_level_analysis.py 产出）。
      #   ① 「14–28% of estimates at the four reduced sizes did」
      #   ② 「9,814 cells at n = 100」
      #   容差按**打印精度**给（INC-022 的规矩）：正文印的是整数百分数（14 / 28），
      #   故容差 = 末位的半个单位 = 0.5（百分点）。源值 14.148% / 28.366% 都落在 ±0.5 内。
      num("§4","S4-45","14–28% of estimates at the four reduced sizes",14.0,"p61_coverage_by_level.csv",
          _cov_min_oor()*100,0.5),
      num("§4","S4-46","14–28% of estimates at the four reduced sizes",28.0,"p61_coverage_by_level.csv",
          _cov_max_oor()*100,0.5),
      num("§4","S4-47","9,814 cells at n = 100",9814,"p61_coverage_by_level.csv",
          _cov_oor_at("100"),0)]
# ★ 2026-10-14 pattern 同步（正当）：L1 那句与 n=100 那句按外部核验要求重写，
#   四处 pattern 随文本更新，**断言值与来源均未变**（仍各指 fig3b / fig5b）。
# ---- §7.2 小节（2026-10-07 新增：protocol §7.2 主任务四项指标）----
#   该节存在的理由：protocol 把「临床部署预警任务」定为**主**任务并规定了四项指标，
#   而此前正文一个词都没有。全部数值来源 p33_protocol72_metrics.csv
#   （该脚本从 p6b/p6c/p6d 独立重算，与 p6d_signal_eval.csv 逐位锚定，差 < 1e-9）。
_P1 = P72[P72.convention == "primary"] if P72 is not None else None
N += [num("§7.2","S72-01","30.6% (57 of 186 transfers)",57,"fig1c_summary.csv",
          float(F1C.n_harmful.iat[0]),0),
      num("§7.2","S72-02","AUROC of 0.644",0.644,"p33_protocol72_metrics.csv",_p33m("AUROC_DR","auroc"),1e-3),
      # ★ 2026-10-08 更正（INC-015）：DR 的 **harmful-detection AUROC 区间**原绑
      #   `p33_protocol72_metrics.csv`（seed 20261007）—— 但同一个量在 §2 正文、Fig 2 所绘、
      #   Supplementary Table S2 所印三处**都取 `fig2_forest.csv`**（seed 20261012）。
      #   两个 bootstrap 种子 ⇒ 同一量两套界：正文 (0.557-0.728) vs 它自己指向的 S2 (0.556-0.729)。
      #   正本改为 `fig2_forest.csv`，与图、表、§2 同一来源。点估计两源逐位相同（0.6436828505371956）。
      #   S72-05…S72-12（AUPRC / sens / PPV / NPV）**不动** —— 那四量的区间只存在于 p33，且 S2 同源。
      num("§7.2","S72-03","(0.556-0.729)",0.556,"fig2_forest.csv",g2("AUROC_DR","auc_lo"),1e-3),
      num("§7.2","S72-04","(0.556-0.729)",0.729,"fig2_forest.csv",g2("AUROC_DR","auc_hi"),1e-3),
      num("§7.2","S72-05","AUPRC of 0.447",0.447,"p33_protocol72_metrics.csv",_p33m("AUROC_DR","auprc"),1e-3),
      num("§7.2","S72-06","(0.335-0.563)",0.335,"p33_protocol72_metrics.csv",_p33m("AUROC_DR","auprc_lo"),1e-3),
      num("§7.2","S72-07","(0.335-0.563)",0.563,"p33_protocol72_metrics.csv",_p33m("AUROC_DR","auprc_hi"),1e-3),
      num("§7.2","S72-08","sensitivity of 0.386",0.386,"p33_protocol72_metrics.csv",_p33m("AUROC_DR","sens80"),1e-3),
      num("§7.2","S72-09","(0.228-0.590)",0.228,"p33_protocol72_metrics.csv",_p33m("AUROC_DR","sens80_lo"),1e-3),
      num("§7.2","S72-10","(0.228-0.590)",0.590,"p33_protocol72_metrics.csv",_p33m("AUROC_DR","sens80_hi"),1e-3),
      num("§7.2","S72-11","positive predictive value of 0.468",0.468,"p33_protocol72_metrics.csv",_p33m("AUROC_DR","ppv"),1e-3),
      num("§7.2","S72-12","negative predictive value of 0.748",0.748,"p33_protocol72_metrics.csv",_p33m("AUROC_DR","npv"),1e-3),
      num("§7.2","S72-13","sensitivities of 0.175",0.175,"p33_protocol72_metrics.csv",_p33m("entropy_js","sens80"),1e-3),
      num("§7.2","S72-14","sensitivities of 0.175",0.175,"p33_protocol72_metrics.csv",_p33m("conf_dist_js","sens80"),1e-3),
      num("§7.2","S72-15","lower still (0.439 to 0.478)",0.439,"fig2_forest.csv",
          float(F2[F2.method.isin(_L1M)].harmful_auc.min()),1e-3),
      num("§7.2","S72-16","lower still (0.439 to 0.478)",0.478,"fig2_forest.csv",
          float(F2[F2.method.isin(_L1M)].harmful_auc.max()),1e-3),
      num("§7.2","S72-17","identified 38.6% of failing",38.6,"p33_protocol72_metrics.csv",
          100 * _p33m("AUROC_DR","sens80"),0.05),
      num("§7.2","S72-18","uniformly high (0.671 to 0.748)",0.671,"p33_protocol72_metrics.csv",
          float(_P1.npv.min()),1e-3),
      num("§7.2","S72-19","uniformly high (0.671 to 0.748)",0.748,"p33_protocol72_metrics.csv",
          float(_P1.npv.max()),1e-3),
      num("§7.2","S72-20","because 69.4% of transfers are",69.4,"fig1c_summary.csv",
          100 * (1 - float(F1C.n_harmful.iat[0]) / float(F1B.eligible_tasks.iat[0])),0.05)]
# ---- §4b（2026-10-14 新增）----------------------------------------------------
# ★ 背景：摘要 v5 报了六个数（src-only 基线 0.611、配对差 +0.033、配对均值对照 +0.0006、
#   配对均值分量 0.638、配对内偏差 0.554、21 对里 16 对），而 **Results 正文当时一个都没有** ——
#   摘要的数字在正文**无家可归**。这与 INC-038 是同一类：**写过的东西没有判据认领，
#   它的消失就不会有红灯照到**。已补一段正文（Results §4「Label-free target-performance
#   estimation…」末），此处补对应断言。
# ★ 值一律**现场读冻结产物**（p57 两个 CSV；21 对那两个数自 p57 报告 fail-closed 抽取），
#   不手打 —— 手打会让审计器变成「拿自己的记忆核对别人的记忆」。
_P57D = _rd("p57_paired_differences.csv")
_P57S = _rd("p57_scores.csv")
_C_DR_SRC   = "prespecified_DR − src_only"
_C_DHAT_SRC = "posthoc_dhat_DR − src_only"
_C_DHAT_AM  = "posthoc_dhat_DR − amean_control"
_C_AM_SRC   = "amean_control − src_only"
def _d57(c, col):
    return float(_P57D[_P57D.contrast == c][col].iat[0])
def _s57(s, col):
    return float(_P57S[_P57S.score == s][col].iat[0])
_P57R = io.open(os.path.join(OUT, "report_p57_posthoc_baselines.txt"), encoding="utf-8").read()
def _g57(pat, what):
    m = _re.search(pat, _P57R)
    if not m:
        raise RuntimeError("p57 报告抽取失败（%s）：pattern=%r —— 报告格式可能已改，"
                           "须人工核对后更新 pattern，禁止回退默认值" % (what, pat))
    return m.group(1)
_PAIR_OK = _g57(r"预先指定分\s*−Â\s*排对\s*(\d+/\d+)", "21 对里排对的比例")
N += [
  num("§4b","S4B-01","harmful-detection AUROC 0.611, 95% interval 0.506 to",0.611,
      "p57_scores.csv",_s57("src_only","point"),1e-3),
  num("§4b","S4B-02","95% interval 0.506 to 0.730",0.506,"p57_scores.csv",_s57("src_only","lo95"),1e-3),
  num("§4b","S4B-03","95% interval 0.506 to 0.730",0.730,"p57_scores.csv",_s57("src_only","hi95"),1e-3),
  num("§4b","S4B-04","the prespecified score reached 0.644",0.644,"p57_scores.csv",
      _s57("prespecified_DR","point"),1e-3),
  num("§4b","S4B-05","paired difference of +0.033 (-0.143 to +0.192)",0.033,
      "p57_paired_differences.csv",_d57(_C_DR_SRC,"point"),1e-3),
  num("§4b","S4B-06","paired difference of +0.033 (-0.143 to +0.192)",-0.143,
      "p57_paired_differences.csv",_d57(_C_DR_SRC,"lo95"),1e-3),
  num("§4b","S4B-07","paired difference of +0.033 (-0.143 to +0.192)",0.192,
      "p57_paired_differences.csv",_d57(_C_DR_SRC,"hi95"),1e-3),
  # ★ S4B-08：**上界写「约 0.2」而不是 +0.19** —— 按 `SEED_STABILITY_TABLE_20261014.md`
  #   自己的规定（5 个种子上该端点最高到 +0.1984 ⇒ 两位小数会走到 0.20），正文一律写 "about 0.2"。
  #   此处同步：pattern 与断言值都跟着正文走，容差由 token 的 1 位小数自动定（半 ULP = 0.05）。
  num("§4b","S4B-08","nor a gain of up to about 0.2",0.2,"p57_paired_differences.csv",
      _d57(_C_DR_SRC,"hi95"),5e-2),
  num("§4b","S4B-09","the change endpoint reached 0.789",0.789,"p57_scores.csv",
      _s57("posthoc_dhat_DR","point"),1e-3),
  num("§4b","S4B-10","(+0.178, +0.091 to +0.253)",0.178,"p57_paired_differences.csv",
      _d57(_C_DHAT_SRC,"point"),1e-3),
  num("§4b","S4B-11","(+0.178, +0.091 to +0.253)",0.091,"p57_paired_differences.csv",
      _d57(_C_DHAT_SRC,"lo95"),1e-3),
  num("§4b","S4B-12","(+0.178, +0.091 to +0.253)",0.253,"p57_paired_differences.csv",
      _d57(_C_DHAT_SRC,"hi95"),1e-3),
  # ★ S4B-13 是本文**最贴 0** 的一条：它的区间端点是全稿里离 0 最近的，故容差按印出的 4 位小数
  #   自动收紧到 5e-5（半 ULP），不靠声明。
  num("§4b","S4B-13","(+0.0006, -0.027 to +0.034)",0.0006,"p57_paired_differences.csv",
      _d57(_C_DHAT_AM,"point"),1e-3),
  num("§4b","S4B-14","(+0.0006, -0.027 to +0.034)",-0.027,"p57_paired_differences.csv",
      _d57(_C_DHAT_AM,"lo95"),1e-3),
  num("§4b","S4B-15","(+0.0006, -0.027 to +0.034)",0.034,"p57_paired_differences.csv",
      _d57(_C_DHAT_AM,"hi95"),1e-3),
  num("§4b","S4B-16","(+0.177, +0.094 to +0.249)",0.177,"p57_paired_differences.csv",
      _d57(_C_AM_SRC,"point"),1e-3),
  num("§4b","S4B-17","(+0.177, +0.094 to +0.249)",0.094,"p57_paired_differences.csv",
      _d57(_C_AM_SRC,"lo95"),1e-3),
  num("§4b","S4B-18","(+0.177, +0.094 to +0.249)",0.249,"p57_paired_differences.csv",
      _d57(_C_AM_SRC,"hi95"),1e-3),
  num("§4b","S4B-19","pair-level mean component (0.638, 0.544 to 0.727)",0.638,"p57_scores.csv",
      _s57("pair_mean_part","point"),1e-3),
  num("§4b","S4B-20","pair-level mean component (0.638, 0.544 to 0.727)",0.544,"p57_scores.csv",
      _s57("pair_mean_part","lo95"),1e-3),
  num("§4b","S4B-21","pair-level mean component (0.638, 0.544 to 0.727)",0.727,"p57_scores.csv",
      _s57("pair_mean_part","hi95"),1e-3),
  num("§4b","S4B-22","within-pair deviation (0.554, 0.476 to 0.634)",0.554,"p57_scores.csv",
      _s57("within_pair_dev","point"),1e-3),
  num("§4b","S4B-23","within-pair deviation (0.554, 0.476 to 0.634)",0.476,"p57_scores.csv",
      _s57("within_pair_dev","lo95"),1e-3),
  num("§4b","S4B-24","within-pair deviation (0.554, 0.476 to 0.634)",0.634,"p57_scores.csv",
      _s57("within_pair_dev","hi95"),1e-3),
  # ⚠ S4B-25/26 的数字只在 **p57 的报告文本**里（无 CSV 产物）⇒ fail-closed 抽取，同 p34 的 ICC。
  num("§4b","S4B-25","ranked the failing target first in 16 of them",16,
      "report_p57_posthoc_baselines.txt",float(_PAIR_OK.split("/")[0]),0),
  num("§4b","S4B-26","first in 16 of them (0.76)",0.76,"report_p57_posthoc_baselines.txt",
      float(_PAIR_OK.split("/")[0]) / float(_PAIR_OK.split("/")[1]),5e-3),
]
# ---- §5 ----
N += [num("§5","S5-01","14 of 15 methods had 95% intervals that included zero",14,"fig5_subsample.csv",
          float(((PA.auroc_lo<=0)&(PA.auroc_hi>=0)).sum()),0),
      # 【2026-10-06 重锚定】原 pattern 把它写成"唯一例外"；修正后唯一例外是 AUROC_DR，
      #   pred_mean_shift 归入"无显著变化"。pattern 改为句内实际措辞。
      num("§5","S5-02","showed no detectable change (+0.021",0.021,"fig5_subsample.csv",gpa("pred_mean_shift","auroc_mean"),1e-3),
      num("§5","S5-03","-0.017 to +0.057",-0.017,"fig5_subsample.csv",gpa("pred_mean_shift","auroc_lo"),1e-3),
      num("§5","S5-04","-0.017 to +0.057",0.057,"fig5_subsample.csv",gpa("pred_mean_shift","auroc_hi"),1e-3),
      # ★【2026-10-06 重锚定·值↔信号名绑定】这是本项目**最典型的"绿而错"**：
      #   原 claim 只比对数字（取全表 max），不校验该数字被归属到哪个信号，
      #   于是修正后最大值换成 AUROC_DR 时，claim 照样 PASS，而正文仍在说
      #   "predicted-prevalence shift" —— 直到对抗复核才被发现。
      #   现拆成两条，**各自从具名方法的行取数**（名与值同时绑定）：
      num("§5","S5-05","largest range among the twelve input- and prediction-space signals was 0.047 for predicted-prevalence shift",
          0.047,"fig5b_across_n_range.csv",
          float(F5B.loc[F5B.method=="pred_pos_prev_shift","across_n_range"].iat[0]),1e-3),
      num("§5","S5-05b","The DR estimator was again the exception, spanning 0.261",0.261,
          "fig5b_across_n_range.csv",
          float(F5B.loc[F5B.method=="AUROC_DR","across_n_range"].iat[0]),1e-3),
      # 【2026-10-06 新增】修正后 §5 首段新增了未报过的 DR 样本量效应，数值须纳入溯源。
      num("§5","S5-09","0.645",0.645,"fig5b_across_n_range.csv",
          float(F5B.loc[F5B.method=="AUROC_DR","auroc_at_full"].iat[0]),1e-3),
      num("§5","S5-10","fell to 0.384 at n = 100",0.384,"fig5b_across_n_range.csv",
          float(F5B.loc[F5B.method=="AUROC_DR","auroc_at_100"].iat[0]),1e-3),
      num("§5","S5-06","median range across the 15 methods was 0.021",0.021,"fig5b_across_n_range.csv",
          float(F5B.across_n_range.median()),5e-4),
      num("§5","S5-07","six of 15 methods had harmful-detection AUROC below 0.5",6,"fig5b_across_n_range.csv",
          float((F5B.auroc_at_full<0.5).sum()),0),
      num("§5","S5-08","eight were also below\n0.5 at n = 100",8,"fig5b_across_n_range.csv",
          float((F5B.auroc_at_100<0.5).sum()),0)]

# ---------------------------------------------------------------- 定性断言
def qual(section, cid, desc, predicate, note=""):
    return dict(kind="qualitative", section=section, claim_id=cid, claim=desc,
                result="PASS" if predicate else "FAIL", note=note)

_both = set(F5B[(F5B.auroc_at_full<0.5)&(F5B.auroc_at_100<0.5)].method)
_L1   = {"domain_auc", "mmd", "psi_mean", "psi_max", "js_mean", "js_max"}
# Q8 用：MMD 两个 learner 的方向比例与可定义折数
_FL = float(CN9.set_index("method").loc["mmd", "frac_logit_pos"])
_FR = float(CN9.set_index("method").loc["mmd", "frac_rank_pos"])
_NF = int(round(float(ORI.set_index("method").loc["mmd", "n_eval_fold"])))
# Q9/Q10 用：benchmark 构成，**从 figS6 逐行派生**（不再取 protocol 冻结表）
_F6  = R("figS6_fold_composition.csv")
_ORG = tuple(int(v) for v in
             _F6.organism.value_counts().reindex(["EC", "KP", "PA"]).fillna(0))
_TR  = _F6.n.value_counts()
_DIR = (int(_TR.get(6, 0)), int(_TR.get(4, 0)), int(_TR.get(2, 0)))
# Q11 用：被排除候选的 source AUROC 区间，**从 P1 History+LR 的 eligible==False 派生**
_F1E  = R("fig1e_configs_full.csv")
_EX   = _F1E[(_F1E.tier == "hist") & (_F1E.algo == "lr") & (~_F1E.eligible.astype(bool))]
_EXN  = int(F1B.excluded.iat[0])
_EXLO = float(_EX.src_auroc.min())
_EXHI = float(_EX.src_auroc.max())
Q = [
 # 【2026-10-01 修正】Q1 原本断言"Fig 3A 使用冻结原生方向"，这是**错误前提**：
 # 数据追溯显示 fig3a_lodo.csv 的 lodo_macro_auroc 与 figS_orientation.csv 的 macro_logit
 # **15/15 完全相同**，与 macro_within(native) 只有 6/15 相同；p6f.lodo_oof 在训练折拟合
 # LogisticRegression 再应用到留出折。故 Fig 3A = **logistic-learned orientation**。
 # 断言随之改为检查修正后的语义，并新增 Q1b 直接锁定数据层事实。
 qual("§2","Q1","Figure 3A 的 LODO 为 logistic-learned orientation（经数据追溯确认）",
      ("learned from" in MS_TXT and "other groups rather than fixed to the signal" in MS_TXT)
      or "learned from the remaining groups" in MS_TXT),
 qual("§2","Q1b","fig3a_lodo 的 lodo_macro_auroc 与 macro_logit 逐方法一致、与 macro_within 不一致",
      bool((A3L["lodo_macro_auroc"].round(9).values == ORI["macro_logit"].round(9).values).all()
           and not (A3L["lodo_macro_auroc"].round(9).values == ORI["macro_within"].round(9).values).all()),
      "若此断言失败，说明 Fig 3A 的语义被改动，正文描述必须同步复核"),
 # 【2026-10-06 修正】原断言「4 个越界估计全部为 P. aeruginosa」在修正后**真空通过**：
 #   fig4e_dr_out_of_range.csv 已变为空表，空 Series 的 .all() 恒为 True ——
 #   断言若在空集上通过，等于没断言。改为直接断言"越界数为 0"（更强的正向断言）。
 qual("§4","Q2","修正后无任何 DR 越界估计（原 4 个越界任务已随 CS-05 消失）",
      len(F4O) == 0,
      "原 4 个越界任务均为 P. aeruginosa；现 fig4e_dr_out_of_range.csv 为 0 行"),
 qual("§4","Q3","PA_ATM_SU2MGB 的 ESS 为全样本最小",
      F4L.loc[F4L.w_ESS.idxmin(),"task_id"]=="PA_ATM_SU2MGB"),
 # 【2026-10-06 修正】折数由 25 变 26（可定义折 25→26），描述同步。
 qual("§2","Q4","MMD 反向方向在 26/26 折中被选中", gso("mmd","frac_train_dir_pos")==0.0),
 qual("§4","Q5","calibration 两端均为「预测比实际更极端」（观测值向零收缩）",
      all(abs(float(CAL[k].actual_d_mean.iat[0])) < abs(float(CAL[k].pred_d_mean.iat[0])) and
          abs(float(CAL[k].actual_d_mean.iat[-1])) < abs(float(CAL[k].pred_d_mean.iat[-1])) for k in CAL)),
 # 【2026-10-06 修正·重要】原断言用的是**交集**判据
 #   `(_both & {四个}) == {四个}` —— 只要那四个在 _both 里就为真，
 #   **对 _both 里多出来的成员完全不敏感**。修正后 _both 已扩为全部六个 L1，
 #   该断言仍会通过，属"绿而错"。改为**集合相等**：恰好是六个 L1，不多不少。
 qual("§5","Q6","两档均 <0.5 者恰为全部六个 L1 input-shift 指标（不多不少）",
      _both == _L1,
      f"实测两档均<0.5 的集合 = {sorted(_both)}；原断言只查交集，对多出的成员不敏感"),
 qual("§5","Q7","full 与 n=100 两档 <0.5 的成员不完全相同",
      set(F5B[F5B.auroc_at_full<0.5].method) != set(F5B[F5B.auroc_at_100<0.5].method)),

 # ---- 【2026-10-01 新增】S9 双 learner 对照（clean 版正文新写的一段） ----
 # 【2026-10-06 重写·重要】原断言把该段描述为「**符号稳定性**不同」
 #   （判据写死 16/25 = 0.64）。修正后 MMD 的 logistic learner 在 25/26 折选正
 #   （frac_logit_pos = 0.9615），rank 规则 26/26 选负 —— **两者各自都高度一致**，
 #   差别在**学到的符号相反**，不在稳定性。故断言与描述一并重写，
 #   并新增一条：两个 learner 的一致性都接近 1（这正是"稳定性"说法不成立的依据）。
 qual("§2","Q8",f"MMD：logistic learner 在 {int(round(_FL*_NF))}/{_NF} 折选正，rank 规则 {_NF}/{_NF} 折选负"
                f"（差别在**符号**，不在稳定性）",
      abs(_FL - 25/26) < 1e-9 and _FR == 0.0
      and min(_FL, 1 - _FR) >= 0.95          # 两个 learner 各自的一致性都 ≥95%
      # 注意：norm() 已剥掉 markdown 粗体标记 `**`（见其 docstring 第 4 条），故 pattern 不带星号
      and "differ in the sign they learn" in MS_TXT
      and "the logistic rule positive in 25 of 26 folds, the rank rule negative in all 26" in MS_TXT,
      "原断言 16/25 已过期；且它把差别归因于稳定性，与实测（两者一致性均≥96%）不符"),

 # ---- benchmark 构成（clean 版 §1）----
 # 【2026-10-06 修正·溯源更正·重要】Q9/Q10 原来把预期值取自 **protocol §3 冻结表**
 #   （17/10/11、28/7/3、38、202）。修正后这些全部过期，但它们**与同样过期的正文
 #   两两对上，于是双双报 PASS** —— 是本项目最隐蔽的一类"绿而错"。
 #   现改为**从 figS6_fold_composition.csv 逐行派生**（当前唯一权威），
 #   并要求手稿措辞与之逐字一致。
 qual("§1","Q9",f"物种构成 {_ORG[0]}/{_ORG[1]}/{_ORG[2]} 且合计 {sum(_ORG)} 组（派生自 figS6）",
      _ORG == (17, 9, 11) and sum(_ORG) == 37
      and "17, 9 and 11 organism-drug groups" in MS_TXT,
      f"figS6 派生 = {_ORG}，合计 {sum(_ORG)} 组（原断言取自 protocol 冻结表，已废弃）"),
 qual("§1","Q10",f"方向数构成 {_DIR[0]}/{_DIR[1]}/{_DIR[2]}，合计 {sum(_DIR)} 组且恰好还原 {_DIR[0]*6+_DIR[1]*4+_DIR[2]*2} 个 task（派生自 figS6）",
      _DIR == (23, 10, 4) and sum(_DIR) == 37
      and _DIR[0]*6 + _DIR[1]*4 + _DIR[2]*2 == 186
      and "23 groups contributed all six directed transfers" in MS_TXT,
      f"figS6 派生 = {_DIR}；{_DIR[0]}×6+{_DIR[1]}×4+{_DIR[2]}×2 = "
      f"{_DIR[0]*6+_DIR[1]*4+_DIR[2]*2}"),

 # ---- 【2026-10-01 新增】文本卫生：半截字符串替换造成的断句 ----
 # 本项目已四次因"用旧串做替换、只匹配到半句"而产生语法破损
 # （本次实例：S3 中 "…AUROC in the" + 换行 + "was 0.423…"）。
 # 数值断言不会报这类错，故单列一条。判据刻意收窄以免误报：
 # 只有「行尾是限定词/介词」且「下一行以 be 动词开头」才判失败。
 qual("§all","Q12","正文无半截替换造成的断句（行尾限定词/介词 + 下行为 be 动词）",
      not __import__("re").search(
          r"\b(?:the|a|an|in|of|to|on|at|by|with|from|for)\s*\n\s*(?:was|were|is|are|be|been|being)\b",
          MS_TXT),
      "MS_TXT 已归一化换行，若此处 FAIL 说明归一化失效；正常应恒为 PASS，实际防线在原文换行层"),
 qual("§all","Q12b","原文层面同样无该断句（直接查未归一化文本）",
      not __import__("re").search(
          r"\b(?:the|a|an|in|of|to|on|at|by|with|from|for)[ \t]*\n[ \t]*(?:was|were|is|are|be|been|being)\b",
          open(MS, encoding="utf-8").read()),
      "这是真正生效的那一道：上一版 S3 的 'in the / was 0.423' 会在此被拦下"),
 # 【2026-10-06 修正·溯源更正】原断言把预期值取自 **protocol §3 冻结文本**
 #   （正则抠出「26（source AUROC 0.571-0.650）」）——那是修正前的宇宙。
 #   这正是 GPT 指出的「validator expectation provenance」问题：预期值必须来自
 #   **当前产物**，不是历史冻结件。现改为从 fig1e_configs_full.csv 的
 #   P1 History+LR 中 eligible==False 的行**派生**，并同时校验手稿措辞。
 qual("§1","Q11",f"被排除的 {_EXN} 个 candidate 的 source AUROC 区间与产物一致",
      len(_EX) == _EXN
      and f"{_EXLO:.3f}" == "0.560" and f"{_EXHI:.3f}" == "0.648"
      and f"{_EXN} excluded candidate transfers" in MS_TXT
      and f"between {_EXLO:.3f} and {_EXHI:.3f}" in MS_TXT,
      f"派生自 fig1e_configs_full.csv（P1 History+LR, eligible==False）："
      f"n={len(_EX)}，区间 {_EXLO:.3f}–{_EXHI:.3f}"),
 # 2026-10-07 新增：§3 方向约定披露段的收尾断言。
 #   该句声称「两种口径下 macro 最高的三名都属于 prediction-space」——
 #   这是一条**可算的**断言，故从 p32 派生复核，而不是当作修辞放过。
 #   注意：它**只**针对排名顶端，**不**声称 L2 家族均值高于 L1（后者在 grouped 口径下不成立）。
 qual("§3","Q13","两种口径下 macro 前三名均属 prediction-space（L2）",
      DS is not None and all(
          _one(m) for m in
          (list(DS.sort_values("prod_macro", ascending=False).index[:3])
           + list(DS.sort_values("raw_macro", ascending=False).index[:3]))),
      "派生自 p32_direction_convention.csv；仅断言排名顶端，不断言家族均值"),
 # 2026-10-07 新增：§7.2 小节标题声称「**every** method identified a minority of failures」。
 #   「every」是比「most」更强的断言，必须可算 —— 从 p33 的 15 行 primary sens80 派生。
 #   若将来某方法的 sens@.80 升到 0.5 以上，标题即为假，此处判红。
 qual("§7.2","Q14","§7.2 标题的『every method identified a minority』成立（15 个方法 sens@.80 全 < 0.5）",
      P72 is not None and float(P72[P72.convention == "primary"].sens80.max()) < 0.5,
      "派生自 p33_protocol72_metrics.csv（max sens@.80 = %.3f）"
      % (float(P72[P72.convention == "primary"].sens80.max()) if P72 is not None else float("nan"))),
]

# ---------------------------------------------------------------- 执行
import rephead
for _l in rephead.header(__file__, [MS, PROTO_F,
        os.path.join(FIG, "fig2_forest.csv"), os.path.join(FIG, "fig3a_lodo.csv"),
        os.path.join(FIG, "tableS2_metrics.csv"), os.path.join(OUT, "p33_protocol72_metrics.csv"),
        os.path.join(OUT, "p52_mae_vs_constant.csv")]):
    p(_l)
p("=" * 108)
p("p15 — Results 数值一致性审计")
p("=" * 108)
_sb()
p(f"  手稿: {MS}")
p(f"  角色: manuscript-to-source consistency checker（**不重算任何统计量**）")
p("")

# ★ 例外必须被声明（p54 的 OPEN_GAP 同族：**例外必须被声明、有界、可核销**）。
#   此处是**设计常态**而非缺口，故不设到期日；但每条都写明理由，且**逐条列出**在报告里。
#   未声明而又对不上的 ⇒ 判 `TOKEN≠ASSERTION`（S3-23 改稿时的形态）。
PATTERN_WITHOUT_VALUE = {
    "S3-25": "pattern 是**定位串**，截到「to」为止；断言值 0.686 在其后",
    "S3-32": "同上（定位串截到「to」，断言值 0.345 在其后）",
    "S5-07": "值在正文里写作英文数词 **six**，不是数字 token",
    "S5-08": "同上（**eight**）",
}


def _token_of(pattern, value):
    r"""返回 (正文印出的 token 原文, 其小数位数)；取不到则 (None, None)。

    ★ 2026-10-09（外部核验第二轮）：**判据读 token，证据面就必须出口 token**。
    上一轮我改成由 token 推容差，却**没有把 token 写进 CSV** —— 于是第三方复核时只能
    读 CSV 里的 `manuscript_value`，而那是 **float，吞掉尾随零**（`0.850` 印成 `0.85`）。
    外部核验因此**在验证我这次修复时，又踩进了我这次修复所针对的同一个盲区**。
    ⇒ 判据依赖的输入，必须与判据的结论一起落盘。
    """
    for _t in _re.findall(r"(?<![\w.])([+-]?\d[\d,]*(?:\.\d+)?)", pattern):
        try:
            if abs(float(_t.replace(",", "")) - float(value)) < 1e-12:
                return _t, (len(_t.split(".")[1]) if "." in _t else 0)
        except ValueError:
            continue
    return None, None


def _token_nd(pattern, value):
    r"""从**正文印出的 token** 取小数位数 —— 而不是从 float 取。

    ★★ 2026-10-09（INC-022，外部核验发现）：
    判据原先把容差建立在一个 **float** 上（`manuscript_value`）。而 `0.400` 一旦变成
    float 就是 `0.4`，**尾随零消失** ⇒ 由它推出的「半 ULP」是 0.05 而不是 0.0005，
    **宽了 100 倍**。实测后果：源值 0.400534、正文印 `0.400`（正确舍入应为 `0.401`，
    差 0.000534 > 0.0005）**被放行**。
    ⇒ 判据必须钉在**印出来的那个字符串**上。`pattern` 本身就是正文片段，故从中取。
    返回 None 表示 pattern 里找不到与 value 相等的 token（此时回落到声明的 tol）。
    """
    # ★ 必须**带符号**取 token：正文的负号是 U+2212，经 `norm()` 变 "-"，
    #   而它在数字**之外** —— 只取 `\d…` 会把 -0.024 读成 0.024，与断言永不相等。
    #   （首版即犯此错：17 行负值全被判「断言与正文脱节」。**同族盲区：把数字的形状当边界**。）
    for _t in _re.findall(r"(?<![\w.])([+-]?\d[\d,]*(?:\.\d+)?)", pattern):
        try:
            if abs(float(_t.replace(",", "")) - float(value)) < 1e-12:
                return len(_t.split(".")[1]) if "." in _t else 0
        except ValueError:
            continue
    return None


rows = []
for c in N:
    in_ms = c["pattern"] in MS_TXT
    _tok, _nd = _token_of(c["pattern"], c["manuscript_value"])
    # ★ 若 pattern 里**有数字 token 但都不等于断言值** ⇒ 断言与正文脱节。
    #   S3-23 改稿时正是这样：pattern 已更新为 0.401、断言的 value 还是旧的 0.400，
    #   `_token_nd` 返回 None 而容差**悄悄回落**到 0.001 ⇒ 错值被放行。
    #   ⇒ 脱节时**不得回落**，直接判非 PASS。这是「判据自身的静默放宽」，属 INC-016 族。
    _toks_in_pat = _re.findall(r"(?<![\w.])([+-]?\d[\d,]*(?:\.\d+)?)", c["pattern"])
    _detached = (bool(_toks_in_pat) and _nd is None
                 and c["claim_id"] not in PATTERN_WITHOUT_VALUE)
    # 半 ULP 由**印出的位数**推出；取它与声明容差中**更严**的一个。
    _tol = c["tolerance"] if _nd is None else min(c["tolerance"], 0.5 * 10.0 ** (-_nd))
    val_ok = abs(c["manuscript_value"] - c["source_value"]) <= _tol
    status = ("TOKEN≠ASSERTION" if (in_ms and _detached) else
              "PASS" if (in_ms and val_ok) else ("NOT_IN_MS" if not in_ms else "MISMATCH"))
    # ★ 证据面三列（2026-10-09，外部核验请求）：
    #   `token`       —— 正文**印出来的那个字符串**（保留尾随零）。判据读它，证据面就该出口它；
    #                    否则第三方只能读 float，而 float 吞零 ⇒ 复核者会踩进判据刚修好的同一个盲区。
    #   `token_nd`    —— 由 token 取到的小数位数（取不到时为 None）。
    #   `source_repr` —— 源值的**全精度**表示。边界情形（如 S3-22 的 0.6125 与 0.6124999…）
    #                    在 round(.,6) 下不可分辨，repr 下才自持。
    _note = ("" if in_ms else "手稿中未找到该片段（可能已被编辑）")
    if c["claim_id"] in PATTERN_WITHOUT_VALUE:
        # 声明必须**落在证据面上**，不能只住在代码里（「聚合 PASS」那课的小号）
        _note = ("PATTERN_WITHOUT_VALUE｜" + PATTERN_WITHOUT_VALUE[c["claim_id"]]
                 + ("" if in_ms else "；另：手稿中未找到该片段"))
    rows.append(dict(section=c["section"], claim_id=c["claim_id"], kind="numeric",
                     manuscript_value=c["manuscript_value"], token=_tok, token_nd=_nd,
                     source_file=c["source_file"],
                     source_value=round(c["source_value"], 6),
                     source_repr=repr(float(c["source_value"])), tolerance=_tol,
                     status=status, note=_note))
D = pd.DataFrame(rows)
D.to_csv(os.path.join(OUT, "p15_results_numeric_audit%s.csv" % _SUF), index=False)

QD = pd.DataFrame(Q); QD.to_csv(os.path.join(OUT, "p15_results_qualitative_audit%s.csv" % _SUF), index=False)

npass = int((D.status == "PASS").sum())
p(f"【数值断言】共 {len(D)} 项：PASS {npass}｜非 PASS {len(D)-npass}")
for st in ("MISMATCH", "NOT_IN_MS"):
    sub = D[D.status == st]
    if len(sub):
        p(f"\n  --- {st}（{len(sub)} 项）---")
        for _, r in sub.iterrows():
            p(f"    {r.claim_id:<10}{r.section:<4} 手稿={r.manuscript_value:<12} "
              f"数据={r.source_value:<12} 来源={r.source_file}  {r.note}")

# ★★ 2026-10-09（外部核验退回一项，理由成立）：
#   外部核验说「p15 的 S4-37…44 在我这份报告里不存在」，因此**拒绝在它读不到的行上关账**。
#   查清后：那八行**确实在**，但在 `p15_results_numeric_audit.csv` 里，**不在 .txt 报告里** ——
#   报告原先**只逐行打印非 PASS 项**，PASS 只进总数。
#   ⇒ 「聚合 PASS 只能整体相信；具体主张需要可枚举的行。」
#   这是 INC-021 的变体，长在**报告自己的格式**上：判据覆盖了，但**证据面**没覆盖。
#   修法：报告**逐项枚举全部断言**（含 PASS），按 section 分组并给每组小计。
p("")
p("【数值断言 —— 逐项枚举】（含 PASS。**聚合 PASS 只能整体相信，具体主张需要可枚举的行**）")
for _sec in sorted(D.section.unique(), key=lambda s: (len(s), s)):
    _sub = D[D.section == _sec]
    _np = int((_sub.status == "PASS").sum())
    p(f"  {_sec}（{len(_sub)} 项）PASS {_np}｜非 PASS {len(_sub)-_np}")
    for _, r in _sub.iterrows():
        # 印 **token**（正文里那个字符串，保留尾随零）而不是 float —— 判据读 token，
        # 证据面就出口 token。取不到 token 的（定位串／英文数词）回落到数值。
        _shown = r.token if isinstance(r.token, str) else str(r.manuscript_value)
        p(f"    {r.claim_id:<10}{r.status:<10}正文={_shown:<12}"
          f"源={str(r.source_value):<14}{r.source_file}")
    p("")

p("")
p(f"【定性断言】共 {len(QD)} 项：PASS {int((QD.result=='PASS').sum())}｜"
  f"FAIL {int((QD.result=='FAIL').sum())}")
for _, r in QD.iterrows():
    p(f"    {r.claim_id:<6}{r.result:<6}{r.claim}")

ok = (npass == len(D)) and (QD.result == "PASS").all()
p("")
p("*** 审计通过 ***" if ok else "*** 审计未通过 —— 见上表 ***")
rep.close()
print("\ndone ->", OUT)
sys.exit(0 if ok else 1)
