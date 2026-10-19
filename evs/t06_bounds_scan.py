# -*- coding: utf-8 -*-
r"""t06_bounds_scan — family BOUND: definitional-range scan of every frozen result CSV.

================================================================================
WHAT THIS MODULE IMPLEMENTS, AND FROM WHERE (independence statement)
================================================================================
This is a pure *data-file* bounds auditor.  It does not re-implement, import or
execute any production estimator; the only production artifact used is the shared
loader ``evs_load`` (path constants + the RESULT output contract).  Every bound
below is taken from a mathematical definition, a published formula, or a frozen
spec header that is quoted verbatim in this docstring — never from the production
implementation of the same statistic.  The scanned inputs are frozen outputs
(CSVs) only.

THE BOUND TABLE (rule id -> bound -> source of the expectation)
--------------------------------------------------------------
R01 AUROC_DR in [0,1]           DIAGNOSTIC.  A finite-sample augmented / doubly-robust
    estimator is not a probability by construction; t06 brief classifies AUROC_DR as
    diagnostic-only ("unconstrained augmented estimator").  Out-of-range => WARN.
R02 d_hat_DR in [-1,1]          DIAGNOSTIC (inherits the DR class).
R03 d_hat_IOW / d_hat_OM        HARD [-1,1]: difference of two AUROC-scale values in [0,1].
R04 d_auroc in [-1,1]           HARD: d_auroc = AUROC_target - AUROC_source_internal,
    difference of two [0,1] rank statistics.
R05 SUDO_AURCC in [0,1]         HARD: frozen formula header of p6c_01_L3_blinded.py:
    "AURCC = trapezoid area of (completeness in [0,1], cumulative |SUDO| in [0,1])".
R06 signed SUDO in [-1,1]       HARD: same frozen header "signed SUDO = 2*(AUC_aux - 0.5)".
R07 MMD2 >= -1/(m-1) - 1/(n-1)  HARD: lower bound of the UNBIASED U-statistic MMD^2 for a
    PSD kernel with k(x,x)=1 (Gretton et al. 2012 JMLR; the project's own p25 header states
    the same bound).  m,n = minor(1000, n_source/n_target) because p6b subsamples
    min(MMD_M=1000, n); figS5 carries m_fixed=1000 => floor = -0.002002002...  The bound is
    applied to the estimator columns themselves and to any (order-preserving) summary of
    them (min/median/p05/p95), since a statistic of values >= floor is >= floor.
R08 figS4 rows metric=='mmd'    same MMD floor, applied row-wise to median/p05/p95/vmin/vmax.
R09 mmd_sd, sd_median/min/max   HARD >= 0 (standard deviation).
R10 JS-type divergences         HARD >= 0 (Lin 1991: JS divergence is symmetric and >= 0);
    diagnostic upper check <= ln 2 = 0.6931 for natural-log JS (Cover & Thomas) -> WARN only
    (base-2 JS would be <= 1, so the upper bound is deliberately not a hard FAIL).
R11 psi_* >= 0                  HARD: PSI = sum (p-q) ln(p/q) >= 0.
R12 w_P50/P95/P99/max >= 0      HARD: density-ratio weight w = e/(1-e) >= 0.
R13 ESS-type >= 1               HARD: Cauchy-Schwarz, ESS = (sum w)^2 / sum w^2 in [1,N].
R14 w_ESS <= n_src_test (p6c)   HARD: Cauchy-Schwarz identity, N = the sample the weights
    are actually built on.  ★ 2026-10-06 ESS-04 收紧：原写 n_source。CS-09 后 p6c 的 ESS
    建在 w_te（source internal test）上，界必须同步收紧 —— 否则判据静默空转
    （修正后 max(ESS/n_source)=0.147，旧界永不触发）。双向实测见 _ess_upper_array。
R15 AUROC / AUPRC family [0,1]  HARD: rank-statistic definition
    AUROC = P(score_pos > score_neg) + 0.5 P(tie), AUPRC = area under a PR curve.
    (excludes AUROC_DR / d_auroc / d_hat_* / *definable, which have their own rules)
R15b fig5_subsample auroc_*     HARD [0,1] for absolute levels, HARD [-1,1] for the
    "PAIRED_full_minus_100" rows, which by the file's own level label are paired
    DIFFERENCES of AUROCs rather than AUROCs.
R16 indicator columns {0,1}     HARD (harmful/stable/eligible/... are 0/1 flags).
R17 sign columns [-1,1]         HARD (sign statistics in {-1,0,1}).
R18 counts >= 0 and integral    HARD (n_*, B, seed, bin, order, ...).
R19 pct_* in [0,100]            HARD (percentages; fig1e/figS3 store percent units).
R20 fractions in [0,1]          HARD (frac_*, sens_at_spec80).
R21 correlations in [-1,1]      HARD (Spearman rho / CI bounds); file-specific rules for
    ci_lo/ci_hi (figS2/figS7), lo_95/hi_95 (fig4e = bootstrap CI of a Spearman rho).
R22 AUROC CIs in [0,1]          HARD: CI of a macro AUROC (fig3a/fig3b/figS_orientation/
    p6g_direction: lo_95/hi_95/within_*/logit_*/rank_*).
R23 macro-* in [0,1]            HARD: macro-average of AUROC-scale members.
R24 pooled AUROC-scale [0,1]    HARD (pooled_OOF, pooled_raw_all202, A/B/C pooled, C_lo...).
R25 fold quantiles in [0,1]     HARD (fold_min/p25/median/p75/max are fold-AUROC quantiles).
R26 AtB/BtC/effects [-1,1]      HARD: differences of two AUROC-scale quantities.
R27 fig1c/fig1d/figS3 stats     HARD [-1,1]: col is a statistic (mean/median/quantile/
    min/max/dmin/dmax) of d_auroc, itself in [-1,1].
R28 sd/std >= 0                 HARD.
R29 *_minus_* in [-2,2]         HARD: spread between two quantiles of a [-1,1] variable.
R30 fig1e quantiles [-1,1]      HARD (aligned-d_auroc quantiles).
R31 pred_mean_shift / pred_pos_prev_shift in [-1,1]  HARD (differences of probabilities/rates).
R32 entropy_mean_shift in [-ln2, ln2]  HARD (difference of two binary entropies, natural log).
R33 pred_d_mean / actual_d_mean in [-1,1]  HARD (mean of d_auroc per bin).
R34 mae >= 0                    HARD.
R35 range_min/median/max >= 0   HARD (range of density-ratio weights).
R36 tolerance >= 0              HARD.
R37 ratio_median >= 0           HARD (sd / |mmd|).
R38 confidence_level in [0,1]   HARD (proportion).
R39 coverage >= 0, integral     HARD (in every scanned file "coverage" is a task COUNT, 202).
R40 logit_vs_rank_agreement     DIAGNOSTIC [0,1] (agreement statistic whose exact
    definitional range - fraction in [0,1] vs correlation in [-1,1] - is ambiguous).
R41 p-value columns in [0,1]    HARD (definition); no column in the scanned set matches
    this pattern - reported so the non-coverage is explicit.
S01..S04                        explicit SKIP entries: columns deliberately not scanned,
    each with its reason (mixed-unit audit columns, design thresholds, mixed summary stats,
    free-text annotation parsed as float).

Emit contract: RESULT|<test_id>|<PASS|FAIL|WARN>|<expected: source>|<observed>|<note>
HARD violations => FAIL; DR-class and other diagnostics => WARN.  The module always
exits 0 so the orchestrator can aggregate.

Run:  python t06_bounds_scan.py
"""
import glob
import os
import re
import sys
import traceback

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)
import evs_load as L

TOL = 1e-9
LN2 = float(np.log(2.0))

FIGDIR = os.path.join(L.PROJ, "figures", "data")
FIG_FILES = sorted(os.path.basename(p) for p in glob.glob(os.path.join(FIGDIR, "*.csv")))
PARENT_FILES = [
    "p6b_label_free_signals_blinded.csv",
    "p6c_L3_blinded.csv",
    "p6d_outcome.csv",
    "p6d_signal_eval.csv",
    "p6e_lodo_aggregation.csv",
    "p6f_lodo_settlement.csv",
    "p6g_direction.csv",
    "p6g_fold_signs.csv",
    "p6e_lodo_fold_auc_long.csv",
    "p15_results_numeric_audit.csv",          # optional ("if present")
]
PARENT_REQUIRED = [f for f in PARENT_FILES if f != "p15_results_numeric_audit.csv"]
# Declared manifest of the t06 brief's figure/table input list.  The BOUND-02a precondition
# must check THIS list - a glob of the files that happen to be on disk can never report a
# missing file (fig_missing would be provably empty).
FIG_REQUIRED = [
    "fig1c_summary.csv",
    "fig1d_by_direction.csv",
    "fig1e_configs.csv",
    "fig3a_lodo.csv",
    "fig3b_loo_organism.csv",
    "fig4e_ess_association.csv",
    "fig5_subsample.csv",
    "figS2_config_agreement.csv",
    "figS3_direction_summary.csv",
    "figS4_signal_summary.csv",
    "figS5_mmd_summary.csv",
    "figS7_weight_association.csv",
    "figS_orientation.csv",
    "tableS1_tasks.csv",
    "tableS2_metrics.csv",
    # ★ 2026-10-07：Supplementary Figure S10（**事后**分诊策略模拟）的渲染数据包。
    #   列入本表后，"文件缺失"会被 BOUND-02a 报出，而不是静默通过。
    "figS10_triage_curves.csv",
    "figS10_triage_members.csv",
    "figS10_triage_members_k37.csv",
    "figS10_triage_summary.csv",
    "tableS3_triage.csv",
    "figS10_triage_overlap.csv",
    "figS10_triage_recall_budget.csv",
    "universe_partition_audit.csv",
]

FIG_S4 = "figS4_signal_summary.csv"
FIG_S5 = "figS5_mmd_summary.csv"
FIG_1C = "fig1c_summary.csv"
FIG_1D = "fig1d_by_direction.csv"
FIG_1E = "fig1e_configs.csv"
FIG_3A = "fig3a_lodo.csv"
FIG_3B = "fig3b_loo_organism.csv"
FIG_4E = "fig4e_ess_association.csv"
FIG_S2 = "figS2_config_agreement.csv"
FIG_S3 = "figS3_direction_summary.csv"
FIG_S7 = "figS7_weight_association.csv"
FIG_SO = "figS_orientation.csv"
FIG_5S = "fig5_subsample.csv"
FIG_T2 = "tableS2_metrics.csv"
# ★ 2026-10-07：Supplementary Figure S10（事后分诊策略模拟）的渲染数据包
FIG_S10_C = "figS10_triage_curves.csv"
FIG_S10_M = "figS10_triage_members.csv"
FIG_S10_K = "figS10_triage_members_k37.csv"
FIG_S10_S = "figS10_triage_summary.csv"
# ★ Supplementary Table S3（事后分诊）的逐格溯源表
FIG_T3 = "tableS3_triage.csv"
# ★ 复合×成员重合度 / 50% 召回预算（2026-10-07 新增）
FIG_OVL = "figS10_triage_overlap.csv"
FIG_RB = "figS10_triage_recall_budget.csv"
FIG_UPA = "universe_partition_audit.csv"   # A4.2h
P6C = "p6c_L3_blinded.csv"
P6G_DIR = "p6g_direction.csv"
P15 = "p15_results_numeric_audit.csv"

FLAGS = {"harmful", "stable", "eligible", "in_p1_universe", "auroc_definable",
         "eligible_for_auroc"}
SIGNS = {"learned_sign_EC", "learned_sign_KP", "learned_sign_PA",
         "rank_dir", "logit_sign", "spearman_sign"}


def mmd_floor_array(df):
    """Unbiased MMD^2 U-statistic floor -1/(m-1) - 1/(n-1) (Gretton et al. 2012).

    p6b subsamples min(MMD_M=1000, n) points from each sample, so m,n are read from
    n_source / n_target when present; m_fixed (figS5) pins the value at 1000.
    """
    nrow = len(df)
    m = np.full(nrow, 1000.0)
    n = np.full(nrow, 1000.0)

    def num(col):
        if col in df.columns:
            return pd.to_numeric(df[col], errors="coerce").to_numpy(dtype=float)
        return None

    ca, cb = num("n_source"), num("n_target")
    if ca is not None:
        m = np.minimum(m, np.where(np.isfinite(ca), ca, 1000.0))
    if cb is not None:
        n = np.minimum(n, np.where(np.isfinite(cb), cb, 1000.0))
    cm = num("m_fixed")
    if cm is not None:
        m = np.minimum(m, np.where(np.isfinite(cm), cm, 1000.0))
        n = np.minimum(n, np.where(np.isfinite(cm), cm, 1000.0))
    m = np.maximum(m, 2.0)
    n = np.maximum(n, 2.0)
    return -1.0 / (m - 1.0) - 1.0 / (n - 1.0)


def _ess_upper_array(df):
    """ESS 的上界 = **该 ESS 所依据的那组权重**的样本量（Cauchy-Schwarz）。

    ★ 2026-10-06（ESS-04 收紧，CS09_ESS_acceptance.md §3 B 组）：
      原实现用 `n_source`。CS-09 之前 p6c 的 ESS 确实建在 `w_all`（全 source）上 ——
      界与实现**自洽，但两者都错**：估计量真正使用的权重是 `w_te`（source internal test）。
      于是 170/202 个 task 的 ESS **超过了 n_src_test**，最大 6.44 倍 ——
      一个建在 N 个权重上的 ESS 不可能超过 N，这是定义性违反。

      CS-09 把实现改到 `w_te` 后，界必须同步收紧，否则判据**静默空转**：
      实测修正后 max(ESS / n_source) = 0.147，旧界永远触发不了。

      双向实测（判据有鉴别力，非空转）：
        修正前产物  170/202 违约，max(ESS/n_src_test)=6.436  ->  判红
        重跑后产物    0/186 违约，max(ESS/n_src_test)=0.941  ->  判绿
    """
    for col in ("n_src_test", "n_source"):
        if col in df.columns:
            v = pd.to_numeric(df[col], errors="coerce").to_numpy(dtype=float)
            return np.where(np.isfinite(v), v, np.inf)
    return np.full(len(df), np.inf)


# ------------------------------------------------------------------ rule table

# ================================================================ 域声明（A3.5，2026-10-08）
# 外部核验给定 schema：BOUND/EVS 条目从「值」升级为「**值 ＋ 域 ＋ 出处**」。
#   **域** ∈ raw / candidate228 / eligible186 / task-expanded / composite / config
#   **出处** = 冻结产物 ＋ 列位（已在各规则的 files/cols 里）
#
# **为什么需要域**：宇宙再冻结时，失效传播有**两个粒度** ——
#   顺依赖边找到**脏文件**（p50 已做），顺**域字段**找到**脏整数**（本表）。
#   一个数若属 `eligible186` 或 `task-expanded`，宇宙一变它**必须**变；
#   属 `raw` 或 `config` 的则不受影响 —— 这正是「两轴分类」的机械形态。
#
# **手写、单一出处**（与 p50 的 EDGES 同法）：键 = 规则 id 的**前缀**，最长匹配优先。
DOMAIN_OF = [
    # ---- config：与宇宙无关的设计常量 ----
    ("R36_tol", "config"), ("R16_flag", "config"), ("S02_thresholds", "config"),
    ("S04", "config"),
    # ---- candidate228：候选网格口径（38 组 × 6 向 = 228）----
    ("R46_s10_k", "candidate228"), ("R55_t3_k", "candidate228"),
    ("R54_t3_kfrac", "candidate228"), ("R59_ovl_k", "candidate228"),
    ("R60_ovl_kfrac", "candidate228"), ("R61_ovl_n", "candidate228"),
    ("R61b_ovl_frac", "candidate228"), ("R62_rb_n", "candidate228"),
    ("R63_rb_frac", "candidate228"),
    # ---- raw：原始站点数据口径（血统与宇宙无关）----
    ("R36", "raw"), ("R18_counts", "raw"),
    # ---- task-expanded ----
    ("R47_s10_cumfail", "task-expanded"), ("R56_t3_caught", "task-expanded"),
    ("R50_s10_caught", "task-expanded"),
    # ---- composite：复合对照（六个 L1 的秩平均）----
    ("R48_s10_precision", "composite"), ("R49_s10_recall", "composite"),
    ("R51_s10_ci", "composite"), ("R52_s10_ci_order", "composite"),
    ("R53_s10_diffrand", "composite"), ("R57_t3_prec", "composite"),
    ("R58_t3_hasci", "composite"),
]
_DEFAULT_DOMAIN = "eligible186"   # 其余一律按合格宇宙口径（per-transfer / 每方法汇总）


def domain_of(rule_id):
    """最长前缀匹配；未命中则回落到 eligible186（**默认值本身是一条声明**）。"""
    best = None
    for pre, dom in DOMAIN_OF:
        if rule_id.startswith(pre) and (best is None or len(pre) > len(best[0])):
            best = (pre, dom)
    return best[1] if best else _DEFAULT_DOMAIN

# keys: id, kind (HARD/DIAG/SKIP), pat/cols/files/files_exclude/not_pat,
#       lo/hi (scalar), lo_fn/hi_fn (array builder), integral, flag01,
#       cond (row predicate), expect (source of the expectation), justify.
RULES = [
    # ---- DR class: diagnostic only -----------------------------------------
    dict(id="R01_auroc_dr", kind="DIAG", cols={"AUROC_DR"}, lo=0.0, hi=1.0,
         expect="DIAGNOSTIC: unconstrained augmented (DR) estimator - may leave [0,1]",
         justify="t06 brief: AUROC_DR is DR-class (unconstrained augmented estimator); "
                 "out-of-range is reported, never called an implementation error"),
    dict(id="R02_dhat_dr", kind="DIAG", cols={"d_hat_DR"}, lo=-1.0, hi=1.0,
         expect="DIAGNOSTIC: difference of an unconstrained DR estimator and an AUROC",
         justify="inherits the DR class"),
    # ---- hard ranges --------------------------------------------------------
    dict(id="R03_dhat_aux", kind="HARD", cols={"d_hat_IOW", "d_hat_OM"}, lo=-1.0, hi=1.0,
         expect="definition: difference of two AUROC-scale values in [0,1]",
         justify="IOW/OM are weighted concordance ratios, bounded in [0,1]"),
    dict(id="R04_d_auroc", kind="HARD", cols={"d_auroc"}, lo=-1.0, hi=1.0,
         expect="definition: d_auroc = AUROC_target - AUROC_source_internal, both in [0,1]",
         justify="difference of two rank statistics"),
    dict(id="R05_sudo_aurcc", kind="HARD", cols={"SUDO_AURCC", "SUDO_absmean"}, lo=0.0, hi=1.0,
         expect="frozen spec p6c header: AURCC = area under (completeness, |SUDO|) curve; "
                "both axes in [0,1]",
         justify="p6c_01_L3_blinded.py frozen formula header (L3-A block)"),
    dict(id="R06_sudo_signed", kind="HARD", cols={"SUDO_signed_min", "SUDO_signed_max"},
         lo=-1.0, hi=1.0,
         expect="frozen spec p6c header: signed SUDO = 2*(AUC_aux - 0.5) in [-1,1]",
         justify="p6c_01_L3_blinded.py frozen formula header (L3-A block)"),
    dict(id="R07_mmd_floor", kind="HARD",
         cols={"mmd", "mmd_min", "mmd_max", "mmd_median", "mmd_p05", "mmd_p95"},
         lo=None, lo_fn=mmd_floor_array, hi=None,
         expect="theory: unbiased MMD^2 U-statistic >= -1/(m-1) - 1/(n-1) for a PSD kernel "
                "with k(x,x)=1 (Gretton 2012; p25 header); m=n=1000 => >= -0.002002",
         justify="a statistic of values >= floor is >= floor; m,n from n_source/n_target "
                 "(p6b subsamples min(1000, n)), m_fixed=1000 in figS5; independent numerical "
                 "check of the same estimator family (RBF, optimized configurations, m=n=2..8) "
                 "attains -2/m and never goes below the cited floor"),
    dict(id="R08_mmd_floor_figs4", kind="HARD", files={FIG_S4},
         cols={"median", "p05", "p95", "vmin", "vmax"},
         cond=lambda r: str(r.get("metric")).strip().lower() == "mmd",
         cond_desc="metric=='mmd'",
         lo=None, lo_fn=mmd_floor_array, hi=None,
         expect="theory: summary stats of unbiased MMD^2 values must respect the MMD^2 floor "
                "(Gretton 2012; p25 header)",
         justify="row-wise rule on the metric=='mmd' row of figS4"),
    dict(id="R08b_auroc_dr_figs4", kind="DIAG", files={FIG_S4},
         cols={"median", "p05", "p95", "vmin", "vmax"},
         cond=lambda r: str(r.get("metric")).strip().lower() == "auroc_dr",
         cond_desc="metric=='AUROC_DR'",
         lo=0.0, hi=1.0,
         expect="DIAGNOSTIC: summary stats of the unconstrained DR estimator",
         justify="same DR-class justification as R01"),
    dict(id="R09_sd", kind="HARD", cols={"mmd_sd", "sd_median", "sd_min", "sd_max", "sd", "std"},
         pat=None, lo=0.0, hi=None,
         expect="definition: standard deviation >= 0",
         justify="sd/std are dispersion measures"),
    dict(id="R10_js_lower", kind="HARD", pat=r"^js_|_js$", lo=0.0, hi=None,
         expect="definition: Jensen-Shannon divergence >= 0 (Lin 1991)",
         justify="symmetric divergence of two distributions"),
    dict(id="R10b_js_upper", kind="DIAG", pat=r"^js_|_js$", lo=None, hi=LN2,
         expect="diagnostic: JS divergence with natural log <= ln 2 = 0.6931 (Cover & Thomas)",
         justify="upper bound is base-dependent (base-2 JS <= 1) -> WARN only"),
    dict(id="R11_psi", kind="HARD", pat=r"^psi_", lo=0.0, hi=None,
         expect="definition: PSI = sum (p-q) ln(p/q) >= 0",
         justify="population stability index is a f-divergence"),
    dict(id="R12_weights", kind="HARD", cols={"w_P50", "w_P95", "w_P99", "w_max"},
         lo=0.0, hi=None,
         expect="definition: density-ratio weight w = e/(1-e) >= 0",
         justify="weights are non-negative by construction"),
    dict(id="R13_ess_lower", kind="HARD", cols={"w_ESS", "ess_min", "ess_median", "ess_max"},
         lo=1.0, hi=None,
         expect="Cauchy-Schwarz: ESS = (sum w)^2 / sum w^2 in [1, N]",
         justify="ESS is at least 1 for non-negative weights"),
    dict(id="R14_ess_upper", kind="HARD", files={P6C}, cols={"w_ESS"},
         lo=None, hi=None, hi_fn=_ess_upper_array,
         expect="Cauchy-Schwarz: ESS <= N, where N = n_src_test —— **该 ESS 所依据的权重域**"
                "（CS-09 后 p6c 的 ESS 建在 w_te = source internal test 上）。"
                "★ 2026-10-06 ESS-04 收紧：原界写 n_source，与修正前的实现自洽但两者都错",
         justify="upper bound of the ESS identity; 双向实测：修正前 170/202 违约"
                 "（max 6.436）判红，重跑后 0/186（max 0.941）判绿 —— 判据有鉴别力"),
    dict(id="R15_auroc", kind="HARD", pat=r"(?i)auroc|_auc|^auc|auprc",
         not_pat=r"(?i)auroc_dr|definable|^d_auroc|^d_hat_", files_exclude={FIG_5S},
         lo=0.0, hi=1.0,
         expect="definition: AUROC = P(score_pos > score_neg) + 0.5 P(tie) in [0,1]; "
                "AUPRC = area under a PR curve in [0,1]",
         justify="rank-statistic definition"),
    dict(id="R15b_fig5_cond", kind="HARD", files={FIG_5S},
         cols={"auroc_mean", "auroc_lo", "auroc_hi"},
         lo=None, hi=None,
         lo_fn=lambda df: np.where(df["level"].astype(str).str.contains("minus").to_numpy(),
                                   -1.0, 0.0),
         hi_fn=lambda df: np.full(len(df), 1.0),
         expect="definition: AUROC level rows in [0,1]; the file's own 'PAIRED_full_minus_100' "
                "level rows are paired DIFFERENCES of AUROCs => [-1,1]",
         justify="level column is part of the frozen file's semantics"),
    dict(id="R16_flag", kind="HARD", cols=FLAGS, flag01=True,
         expect="definition: indicator column must be in {0,1}",
         justify="binary design/outcome flags"),
    dict(id="R17_sign", kind="HARD", cols=SIGNS, lo=-1.0, hi=1.0,
         expect="definition: sign statistic in {-1,0,1} subset of [-1,1]",
         justify="sign columns"),
    dict(id="R18_counts", kind="HARD",
         pat=r"^(n(_|$)|B$|B_requested$|B_valid_class_check$|bin$|count$|order$|group_rank$|"
             r"dir_rank$|direction_order$|oor$|m_fixed$|SUDO_nbin$|events$|SR$|RS$|seed$|"
             r"n_dr_gt1$|fold_n$)",
         cols={"candidate_tasks", "eligible_tasks", "excluded"},
         lo=0.0, hi=None, integral=True,
         expect="definition: a count is a non-negative integer",
         justify="counts/indices"),
    dict(id="R19_pct", kind="HARD", pat=r"^pct_|_pct$", lo=0.0, hi=100.0,
         expect="definition: percentage in [0,100]",
         justify="percent-unit columns (pct_harmful_aligned, coverage_pct, SR_pct)"),
    dict(id="R20_frac", kind="HARD", pat=r"^frac_|^sens_at_spec", lo=0.0, hi=1.0,
         expect="definition: fraction / sensitivity in [0,1]",
         justify="fractions of folds or cases"),
    dict(id="R21_corr", kind="HARD", pat=r"spearman|rho|sp_lo|sp_hi", lo=-1.0, hi=1.0,
         expect="definition: correlation coefficient in [-1,1]",
         justify="Spearman/Pearson rho and its CI bounds"),
    dict(id="R21b_ci", kind="HARD", files={FIG_S2, FIG_S7}, cols={"ci_lo", "ci_hi"},
         lo=-1.0, hi=1.0,
         expect="definition: CI of a correlation coefficient in [-1,1]",
         justify="ci_lo/ci_hi in these two files bound a Spearman rho"),
    dict(id="R21c_ci4e", kind="HARD", files={FIG_4E}, cols={"lo_95", "hi_95"},
         lo=-1.0, hi=1.0,
         expect="definition: bootstrap CI of a Spearman rho in [-1,1]",
         justify="fig4e_ess_association reports rho of ESS vs DR error"),
    dict(id="R22_auroc_ci", kind="HARD", files={FIG_3A, FIG_3B, FIG_SO, P6G_DIR},
         cols={"lo_95", "hi_95", "within_lo", "within_hi", "logit_lo", "logit_hi",
               "rank_lo", "rank_hi"},
         lo=0.0, hi=1.0,
         expect="definition: CI of a macro AUROC in [0,1]",
         justify="bootstrap CI bounds of AUROC-scale macro averages"),
    dict(id="R23_macro", kind="HARD", pat=r"macro", not_pat=r"auroc", lo=0.0, hi=1.0,
         expect="definition: macro-average of AUROC-scale values in [0,1]",
         justify="macro_*/size_weighted_macro/org-macro columns"),
    dict(id="R24_pooled", kind="HARD",
         cols={"pooled_OOF", "pooled_raw_all202", "pooled_202task", "pooled_138task",
               "A_pooled_all202", "B_pooled_eval", "C_lo", "C_hi", "B_lo", "B_hi"},
         lo=0.0, hi=1.0,
         expect="definition: pooled AUROC-scale statistic in [0,1]",
         justify="pooled/within/organism-level AUROC aggregates"),
    dict(id="R25_fold_q", kind="HARD", pat=r"^fold_(min|p25|median|p75|max)",
         lo=0.0, hi=1.0,
         expect="definition: fold-level AUROC quantile in [0,1]",
         justify="fold AUC summaries"),
    dict(id="R26_deltas", kind="HARD",
         pat=r"^AtB$|^BtC$|^task_subset_effect$|^aggregation_effect$", lo=-1.0, hi=1.0,
         expect="definition: difference of two AUROC-scale quantities in [-1,1]",
         justify="settlement/decomposition differences"),
    dict(id="R27_dstat", kind="HARD", files={FIG_1C, FIG_1D, FIG_S3},
         cols={"mean", "median", "min", "max", "dmin", "dmax", "p05", "p10", "p25", "p50",
               "p75", "p90", "p95", "25%", "50%", "75%"},
         lo=-1.0, hi=1.0,
         expect="definition: statistic of d_auroc, itself in [-1,1]",
         justify="summary statistics preserve the range of the underlying variable"),
    dict(id="R28_std", kind="HARD", cols={"sd", "std"}, lo=0.0, hi=None,
         expect="definition: standard deviation >= 0",
         justify="dispersion is non-negative"),
    dict(id="R29_spread", kind="HARD", pat=r"minus", lo=-2.0, hi=2.0,
         expect="definition: spread of two quantiles of a [-1,1] variable in [-2,2]",
         justify="p90_minus_p10(_aligned)"),
    dict(id="R30_1e", kind="HARD", files={FIG_1E},
         cols={"median_aligned", "p10_aligned", "p90_aligned"}, lo=-1.0, hi=1.0,
         expect="definition: quantile of aligned d_auroc in [-1,1]",
         justify="aligned delta-AUROC quantiles"),
    dict(id="R31_pshift", kind="HARD",
         cols={"pred_mean_shift", "pred_pos_prev_shift"}, lo=-1.0, hi=1.0,
         expect="definition: difference of two probabilities/rates in [-1,1]",
         justify="prediction-space shift statistics"),
    dict(id="R32_entropy_shift", kind="HARD", cols={"entropy_mean_shift"},
         lo=-LN2, hi=LN2,
         expect="definition: difference of two binary entropies (natural log) in [-ln2, ln2]",
         justify="p6b ent(p) = -(p ln p + (1-p) ln(1-p)) in [0, ln2]"),
    dict(id="R33_calib", kind="HARD", cols={"pred_d_mean", "actual_d_mean"}, lo=-1.0, hi=1.0,
         expect="definition: mean d_auroc per calibration bin in [-1,1]",
         justify="mean of a [-1,1] variable"),
    dict(id="R34_mae", kind="HARD", cols={"mae"}, lo=0.0, hi=None,
         expect="definition: mean absolute error >= 0",
         justify="absolute error"),
    dict(id="R35_wrange", kind="HARD", cols={"range_min", "range_median", "range_max"},
         lo=0.0, hi=None,
         expect="definition: range of density-ratio weights (w = e/(1-e) >= 0)",
         justify="figS7 weight ranges"),
    dict(id="R36_tol", kind="HARD", cols={"tolerance"}, lo=0.0, hi=None,
         expect="definition: audit tolerance >= 0",
         justify="p15 numeric audit"),
    dict(id="R37_ratio", kind="HARD", cols={"ratio_median"}, lo=0.0, hi=None,
         expect="definition: sd / |mmd| ratio >= 0",
         justify="figS5 dispersion ratio"),
    dict(id="R38_conf", kind="HARD", cols={"confidence_level"}, lo=0.0, hi=1.0,
         expect="definition: confidence level is a proportion in [0,1]",
         justify="fig2 bootstrap metadata"),
    dict(id="R39_coverage", kind="HARD", cols={"coverage"}, lo=0.0, hi=None, integral=True,
         expect="definition: in every scanned file 'coverage' is a task COUNT (202) => >= 0, "
                "integral",
         justify="fig2/tableS2/p6d_signal_eval coverage columns"),
    dict(id="R40_agree", kind="DIAG", cols={"logit_vs_rank_agreement"}, lo=0.0, hi=1.0,
         expect="DIAGNOSTIC: agreement statistic (fraction in [0,1] or rho in [-1,1]) - "
                "definitional range ambiguous",
         justify="figS_orientation; WARN only"),
    dict(id="R41_pvalue", kind="HARD", pat=r"pval|p_value|p_adj|_p$|qval|fdr", lo=0.0, hi=1.0,
         expect="definition: p-value in [0,1]",
         justify="no column in the required file set matches this pattern (non-coverage is "
                 "reported explicitly)"),
    dict(id="R42_across_n", kind="HARD", cols={"across_n_range"}, lo=0.0, hi=1.0,
         expect="definition: max-min of AUROC values in [0,1]",
         justify="fig5b across-n spread"),
    dict(id="R43_predictive_values", kind="HARD",
         cols={"ppv", "ppv_lo", "ppv_hi", "npv", "npv_lo", "npv_hi"},
         lo=0.0, hi=1.0,
         expect="definition: positive/negative predictive value and their bootstrap CI bounds "
                "are proportions in [0,1]",
         justify="protocol §7.2 warning-task metrics (PPV/NPV at the observed harmful "
                 "prevalence); merged into tableS2_metrics.csv on 2026-10-07 from "
                 "p33_protocol72_metrics.csv (merge-only, no statistic recomputed). "
                 "Registered rather than skipped: a proportion has a hard definitional range."),
    # ---- explicit skips -----------------------------------------------------
    dict(id="S01_p15_mixed", kind="SKIP", files={P15},
         cols={"manuscript_value", "source_value"},
         skip="mixed-unit audit column (counts and signed d_auroc in one column) - no single "
              "definitional range"),
    dict(id="S02_thresholds", kind="SKIP",
         cols={"harmful_threshold", "stable_abs_threshold"},
         skip="design thresholds (signed d_auroc cutoffs) - not derived statistics"),
    dict(id="S03_figs4_mixed", kind="SKIP", files={FIG_S4},
         cols={"median", "p05", "p95", "vmin", "vmax"},
         skip="summary stats of a mixed-unit signal column; row-wise bounds applied only for "
              "metric=='mmd' (R08) and metric=='AUROC_DR' (R08b)"),
    dict(id="S04_p15_note", kind="SKIP", files={P15}, cols={"note"},
         skip="free-text annotation parsed as float (empty)"),
    # ★ 2026-10-09（INC-022 的证据面改动）：`p15` 新增两列，均由**判据的输入**直接派生，
    #   不是统计量，故 SKIP —— 但**必须显式声明**（规则 26：覆盖要点名），否则就是覆盖缺口。
    dict(id="S05_p15_evidence_face", kind="SKIP", files={P15},
         cols={"token_nd", "source_repr"},
         skip="evidence-face columns added 2026-10-09: `token_nd` = decimal count of the printed "
              "token (a display index, 0-6); `source_repr` = full-precision repr of the source "
              "value, carrying the same definitional domain as `source_value`, which S01 covers. "
              "Neither is a derived statistic, so neither has a definitional bound."),
]

# ================================================================ R44 / R45
# CI 成对关系（2026-10-07 新增）
#
# 起因：PPV / NPV 的 bootstrap 区间经 p11 的 merge-only 扩展进入 tableS2 后，
#   「下界 <= 上界」与「点估计落在区间内」两条关系**无人检查** —— 原先只有逐列的
#   值域规则（R15_auroc / R43），而 pivot 或列顺序错误恰好藏在**跨列**关系里。
#
# ★ 分级依据（**与外部建议不同**）：外部建议把包络 lo <= 点估计 <= hi 列为 HARD。
#   不采纳 —— t06 的约定是 HARD = 数学/定义性恒等式。percentile bootstrap 的
#   lo/hi 是**重采样分布**的 2.5/97.5 分位数，点估计来自**原样本**，
#   数学上**不保证**点估计落在区间内 ⇒ 包络只能是 DIAG。
#   真正定义性的是 lo <= hi（分位数天然有序）⇒ 那一条才是 HARD。
#   注：lo <= hi 可由两条包络传递推出，但仍单列 HARD，
#       因为互换 lo/hi 是包络能捕获的最典型故障，需要**硬失败**而非 WARN。
_CI_TRIPLES = [("auc_lo", "auc_hi", "harmful_auc"),
               ("auprc_lo", "auprc_hi", "harmful_auprc"),
               ("sens80_lo", "sens80_hi", "sens_at_spec80"),
               ("ppv_lo", "ppv_hi", "ppv"),
               ("npv_lo", "npv_hi", "npv")]
def _col(colname):
    def f(df):
        return df[colname].to_numpy(dtype=float)
    return f


# ★ 2026-10-12 新增行类判据。
#   缺陷：分诊汇总表 figS10_triage_summary.csv 的 lo95/hi95 两列**不是同质的**。
#     ① 策略层的行 —— 区间是 **precision 本身**（比例），定义域 [0,1]；
#     ② 差值行（layer 名含 "minus"）—— 区间是 **两个 precision 之差**，定义域 [-1,1]。
#   本轮给 p45 追加 L-DR_minus_L-src 一行后，该行下界 -0.138889 被 R51「比例必在 [0,1]」
#   判红。**这是判据的作用域声明错了，不是数据错了**：独立重抽样（4000 次，另一条 RNG）
#   给出下界 -0.1351，且 12.5% 的 replicate 差值为负 —— 差值的下界**本来就可以为负**。
#   处置（不是放宽）：按**行类**分担范围规则，次序规则对两类同等适用。
def _is_diff_row(r):
    return "minus" in str(r.get("layer", "")).lower()
for _lo, _hi, _pt in _CI_TRIPLES:
    RULES += [
        dict(id="R44_ci_order_" + _lo, kind="HARD", files={FIG_T2}, cols={_lo},
             lo=0.0, hi_fn=_col(_hi),
             expect="definition: a lower bootstrap bound cannot exceed the upper bound",
             justify="percentile 2.5 <= percentile 97.5 on the same resample distribution; "
                     "a swap here is a pivot/ordering defect, so it must fail hard"),
        dict(id="R45_ci_env_" + _lo, kind="DIAG", files={FIG_T2}, cols={_lo},
             lo=0.0, hi_fn=_col(_pt),
             expect="observed (NOT definitional): the lower bound lies below the point estimate",
             justify="the point estimate comes from the original sample while lo/hi come from the "
                     "resample distribution, so containment is not guaranteed for a percentile "
                     "bootstrap; WARN only"),
        dict(id="R45_ci_env_" + _hi, kind="DIAG", files={FIG_T2}, cols={_hi},
             lo_fn=_col(_pt), hi=1.0,
             expect="observed (NOT definitional): the upper bound lies above the point estimate",
             justify="same as R45_ci_env_*_lo: containment is empirical, not definitional; WARN only"),
    ]

# ================================================================ R46..R53
# Supplementary Figure S10（**事后**分诊策略模拟）数据包的值域归属（2026-10-07 新增）
#
# 为什么需要：p45 落盘了 4 个新的 figures/data CSV。BOUND-04 是**覆盖审计**——
# 每个数值列必须被某条规则扫到，或被显式 SKIP 并写明理由；未被归属 = 覆盖缺口。
# 新数据进 figures/data 而不登记规则，会让审计出现盲区（数据在，判据不在）。
RULES += [
    dict(id="R46_s10_k", kind="HARD",
         files={FIG_S10_C, FIG_S10_M, FIG_S10_K, FIG_S10_S},
         cols={"k"}, lo=0.0, hi=186.0, integral=True,
         expect="definition: k is a count of reviewed transfers out of the 186-transfer benchmark",
         justify="预注册备忘 D1：k = round(0.20 × n)，观测数据 n = 186；上界即 benchmark 规模"),
    dict(id="R47_s10_cumfail", kind="HARD", files={FIG_S10_C, FIG_S10_M},
         cols={"cum_failures"}, lo=0.0, hi=57.0,
         expect="definition: cumulative number of the 57 failing transfers captured",
         justify="不得超过冻结宇宙中的失效数（fig1c.harmful 合计 = 57）；"
                 "L-random 行是浮点期望值，故不设 integral"),
    dict(id="R48_s10_precision", kind="HARD",
         files={FIG_S10_C, FIG_S10_M, FIG_S10_K, FIG_S10_S},
         cols={"precision", "precision_point", "precision_boot_mean"}, lo=0.0, hi=1.0,
         expect="definition: a precision is a proportion of the reviewed queue",
         justify="比值，定义域 [0,1]"),
    dict(id="R49_s10_recall", kind="HARD", files={FIG_S10_C, FIG_S10_M},
         cols={"recall"}, lo=0.0, hi=1.0,
         expect="definition: recall is a proportion of the 57 failing transfers",
         justify="比值，定义域 [0,1]"),
    dict(id="R50_s10_caught", kind="HARD", files={FIG_S10_K, FIG_S10_S},
         cols={"caught"}, lo=0.0, hi=57.0,
         expect="definition: number of failing transfers captured",
         justify="不得超过全部失效数 57；随机政策那一行是 37×57/186 = 11.34（非整数，故不设 integral）"),
    dict(id="R51_s10_ci", kind="HARD", files={FIG_S10_S}, cols={"lo95", "hi95"},
         cond=lambda r: not _is_diff_row(r),
         cond_desc="layer 不含 'minus'（该行是 precision 本身，不是两个 precision 之差）",
         lo=0.0, hi=1.0,
         expect="definition: percentile bootstrap bounds on a proportion",
         justify="策略层的区间是 precision 分布的分位数，必在 [0,1]（差值行由 R51b 另管）"),
    dict(id="R51b_s10_ci_diff", kind="HARD", files={FIG_S10_S}, cols={"lo95", "hi95"},
         cond=_is_diff_row,
         cond_desc="layer 含 'minus'（两个 precision 之差的 bootstrap 区间）",
         lo=-1.0, hi=1.0,
         expect="definition: percentile bootstrap bounds on a difference of two proportions",
         justify="两个 [0,1] 量之差的区间，定义域 [-1,1]，与 R53 对 diff_vs_random_* 的口径一致。"
                 "★ 本条不是对 R51 的放宽，而是**行类分工**：差值行下界可以为负"
                 "（实测 12.5% 的重抽样 replicate 差值为负）"),
    dict(id="R52_s10_ci_order", kind="HARD", files={FIG_S10_S}, cols={"lo95"},
         hi_fn=_col("hi95"),
         expect="definition: a lower bootstrap bound cannot exceed the upper bound",
         justify="同一次重抽样分布上 2.5% 分位 ≤ 97.5% 分位；顺序颠倒属支点/排序缺陷，必须硬失败。"
                 "★ 2026-10-12：本条**不再同时校验下界 ≥ 0** —— 它是**次序**恒等式，"
                 "对比例行与差值行同等成立；下界的**范围**已由 R51 / R51b 按行类分担，"
                 "去掉重复的后一条不会漏检（比例行的负下界仍被 R51 判红）"),
    dict(id="R53_s10_diffrand", kind="HARD", files={FIG_S10_S},
         cols={"diff_vs_random_lo", "diff_vs_random_hi"}, lo=-1.0, hi=1.0,
         expect="definition: difference between two proportions",
         justify="两个 [0,1] 量之差，定义域 [-1,1]"),
    # ---- Supplementary Table S3（事后分诊）的逐格溯源表（2026-10-07 新增）----
    dict(id="R54_t3_kfrac", kind="HARD", files={FIG_T3}, cols={"k_frac"}, lo=0.0, hi=1.0,
         expect="definition: the review budget as a fraction of the benchmark",
         justify="预注册备忘 D1：k_frac ∈ {0.10, 0.20, 0.30}，是比例不是计数"),
    dict(id="R55_t3_k", kind="HARD", files={FIG_T3}, cols={"k"}, lo=0.0, hi=186.0,
         integral=True,
         expect="definition: k = round(k_frac × 186), a count of reviewed transfers",
         justify="上界 = benchmark 规模 186"),
    dict(id="R56_t3_caught", kind="HARD", files={FIG_T3}, cols={"caught"}, lo=0.0, hi=57.0,
         expect="definition: failing transfers captured among the top k",
         justify="不得超过全部失效数 57；L-random 行是 k×57/186 的**期望值**（非整数，故不设 integral）"),
    dict(id="R57_t3_prec", kind="HARD", files={FIG_T3},
         cols={"precision_full", "precision_shown"}, lo=0.0, hi=1.0,
         expect="definition: precision is a proportion of the reviewed queue",
         justify="precision_full 是全精度值，precision_shown 是它的三位渲染；两者同在 [0,1]"),
    dict(id="R58_t3_hasci", kind="HARD", files={FIG_T3}, cols={"has_ci"}, flag01=True,
         expect="definition: a boolean flag (does this cell carry a bootstrap interval?)",
         justify="布尔量，取值只能是 0 / 1"),
    # ---- 复合×成员重合度、50% 召回预算（2026-10-07 新增，外部完备性核验要求）----
    dict(id="R59_ovl_k", kind="HARD", files={FIG_OVL}, cols={"k"}, lo=0.0, hi=186.0,
         integral=True,
         expect="definition: k = round(k_frac × 186)",
         justify="同 R55"),
    dict(id="R60_ovl_kfrac", kind="HARD", files={FIG_OVL}, cols={"k_frac"}, lo=0.0, hi=1.0,
         expect="definition: review budget as a fraction", justify="同 R54"),
    # ★ 初版把 overlap 与 overlap_frac 合成一条、同给 [0,1] —— **错**：
    #   overlap 是**计数**（实测 15–54），不是比值。BOUND-01 当场判红 18/18 项。
    #   现拆成两条：计数归计数，比值归比值。
    dict(id="R61_ovl_n", kind="HARD", files={FIG_OVL}, cols={"overlap"}, lo=0.0, hi=186.0,
         integral=True,
         expect="definition: size of the intersection of two top-k sets (a count)",
         justify="前 k 名集合的交集不超过 k（k ≤ 186）"),
    dict(id="R61b_ovl_frac", kind="HARD", files={FIG_OVL}, cols={"overlap_frac"}, lo=0.0, hi=1.0,
         expect="definition: overlap / k, a proportion",
         justify="两个 top-k 集合的交集占比，必在 [0,1]"),
    dict(id="R62_rb_n", kind="HARD", files={FIG_RB}, cols={"n_needed"}, lo=0.0, hi=186.0,
         integral=True,
         expect="definition: the number of transfers that must be reviewed to reach 50% recall",
         justify="不超过队列长度 186；整数计数"),
    dict(id="R63_rb_frac", kind="HARD", files={FIG_RB}, cols={"frac_needed"}, lo=0.0, hi=1.0,
         expect="definition: n_needed / 186", justify="比值，定义域 [0,1]"),
    # ---- A4.2h：均匀性升格为**设计不变量**（2026-10-08）----
    #   判据：每个 (organism–drug, site) 组合**恰好** 2 条候选转移 ⇒ 228 = 38×2×3。
    #   这不是「观察到的范围」，是**设计恒等式**：lo == hi == 2。
    #   工件：figures/data/universe_partition_audit.csv（114 行 = 38 组 × 3 站）
    dict(id="R64_uniformity", kind="HARD", files={FIG_UPA}, cols={"n_candidate"},
         lo=2, hi=2, integral=True,
         expect="design identity: every (organism–drug group, site) cell carries exactly 2 "
                "candidate transfers, so 228 = 38 groups x 2 directions x 3 sites",
         justify="A4.2h（2026-10-08）：外部核验预言「每站均匀 ⇒ KP_NIT 三站各 2」并实测命中，"
                 "本规则把该均匀性从**推断**升格为**受断言保护的设计不变量**；"
                 "工件 p47_universe_partition.py 输出（114 行，计数分布 {2:114}）"),
]

ID_COLS = ["task_id", "claim_id", "group", "method", "metric", "item", "config",
           "estimator", "direction", "stage_key", "x",
           # ★ 2026-10-12 追加：分诊产物（p45）与成图数据包的行标识是 `layer`。
           #   此前 ID_COLS 里没有它，判红行只能报 "row#7"，读报告的人无法知道是哪一层
           #   —— 本次定位 L-DR_minus_L-src 的下界缺陷时就先被这个占位绕了一圈。
           #   放在末尾：只在前面所有列都不存在时才生效，不改变任何既有文件的行标签。
           "layer"]


def _matches(rule, fname, col):
    if rule.get("files") and fname not in rule["files"]:
        return False
    if rule.get("files_exclude") and fname in rule["files_exclude"]:
        return False
    if rule.get("cols") is not None and col in rule["cols"]:
        return True
    if not rule.get("pat"):
        return False
    if not re.search(rule["pat"], col):
        return False
    if rule.get("not_pat") and re.search(rule["not_pat"], col):
        return False
    return True


def _san(s):
    return str(s).replace("|", "/").replace("\n", " ")


def _row_label(df, i):
    for c in ID_COLS:
        if c in df.columns:
            return f"{c}={df[c].iat[i]}"
    return f"row#{i}"


def scan_file(path):
    """Return dict with violations, scanned map, skip map, unattributed numerics."""
    fname = os.path.basename(path)
    df = pd.read_csv(path)
    out = dict(file=fname, path=path, rows=len(df), n_hard=0, n_diag=0, viol=[],
               scanned=[], skipped=[], unattributed=[], numeric=[],
               rule_hits=set(), zero_used={})
    for col in df.columns:
        s_raw = df[col]
        if s_raw.dtype.kind not in "ifb":
            out["skipped"].append((col, "non-numeric (string/identifier/bool column)"))
            continue
        out["numeric"].append(col)
        s = pd.to_numeric(s_raw, errors="coerce")
        matched = [r for r in RULES if _matches(r, fname, col)]
        active = [r for r in matched if r["kind"] != "SKIP"]
        skips = [r for r in matched if r["kind"] == "SKIP"]
        if not matched:
            out["unattributed"].append(col)
            continue
        if skips and not active:
            out["skipped"].append((col, skips[0].get("skip", "skipped")))
            continue
        out["scanned"].append(col)
        out["rule_hits"].update(r["id"] for r in active)
        mask0 = s.notna().to_numpy()
        for rule in active:
            lo = rule.get("lo", None)
            hi = rule.get("hi", None)
            if rule["kind"] == "SKIP":
                continue
            if rule.get("lo_fn"):
                lo = rule["lo_fn"](df)
            if rule.get("hi_fn"):
                hi = rule["hi_fn"](df)
            mask = mask0.copy()
            if rule.get("cond"):
                mask &= np.array([bool(rule["cond"](r)) for _, r in df.iterrows()])
            n_used = int(mask.sum())
            if n_used == 0:
                # Fail-closed: a rule that silently evaluates no value must not be reported
                # as coverage (e.g. a relabelled row disables a row-filtered rule, or an
                # all-NaN column is "scanned" but never checked).
                out["zero_used"]["%s::%s" % (rule["id"], col)] = (
                    rule["kind"], rule.get("cond_desc", ""), rule["expect"], col)
                continue
            v = s.to_numpy(dtype=float)
            viol = np.zeros(len(df), dtype=bool)
            lo_arr = np.full(len(df), -np.inf) if lo is None else np.asarray(lo, dtype=float)
            hi_arr = np.full(len(df), np.inf) if hi is None else np.asarray(hi, dtype=float)
            lo_b = np.broadcast_to(lo_arr, (len(df),)) if lo_arr.ndim == 0 else lo_arr
            hi_b = np.broadcast_to(hi_arr, (len(df),)) if hi_arr.ndim == 0 else hi_arr
            viol |= mask & np.isfinite(lo_b) & (v < lo_b - TOL)
            viol |= mask & np.isfinite(hi_b) & (v > hi_b + TOL)
            if rule.get("integral"):
                viol |= mask & (np.abs(v - np.round(v)) > 1e-9)
            if rule.get("flag01"):
                viol |= mask & ~np.isin(v, (0.0, 1.0))
            if not viol.any():
                continue
            idx = np.where(viol)[0]
            order = idx[np.argsort(v[idx])]
            ids = [_row_label(df, int(i)) for i in order[:4]]
            vals = [f"{v[int(i)]:.6g}" for i in order[:4]]
            if len(idx) > 4:
                ids = [_row_label(df, int(i)) for i in order[:2]] + \
                      [_row_label(df, int(i)) for i in order[-2:]]
                vals = [f"{v[int(i)]:.6g}" for i in order[:2]] + \
                       [f"{v[int(i)]:.6g}" for i in order[-2:]]
            lo_s = "-inf" if lo is None else (
                f"min={np.nanmin(lo_b[np.isfinite(lo_b)]):.6g}" if np.ndim(lo)
                else f"{float(lo):.6g}")
            hi_s = "+inf" if hi is None else (
                f"max={np.nanmax(hi_b[np.isfinite(hi_b)]):.6g}" if np.ndim(hi)
                else f"{float(hi):.6g}")
            bound = "bound=[" + lo_s + ", " + hi_s + "]" + \
                    ("+integral" if rule.get("integral") else "") + \
                    ("+flag01" if rule.get("flag01") else "")
            out["viol"].append(dict(
                rule=rule["id"], col=col, kind=rule["kind"], bound=bound,
                expect=rule["expect"], justify=rule["justify"],
                n_viol=int(len(idx)), n_used=n_used,
                vmin=float(np.nanmin(v[idx])), vmax=float(np.nanmax(v[idx])),
                ids=ids, vals=vals,
                cond=rule.get("cond_desc", "")))
            if rule["kind"] == "HARD":
                out["n_hard"] += 1
            else:
                out["n_diag"] += 1
    return out


def main():
    L.banner("t06 bounds scan (family BOUND) — definitional ranges of frozen result CSVs")
    print("# rule table: %d rules, %d hard / %d diag / %d explicit-skip"
          % (len(RULES),
             sum(1 for r in RULES if r["kind"] == "HARD"),
             sum(1 for r in RULES if r["kind"] == "DIAG"),
             sum(1 for r in RULES if r["kind"] == "SKIP")))
    print("# scanned file set: %d figures/data/*.csv + %d parent result CSVs"
          % (len(FIG_FILES), len(PARENT_FILES)))

    # ---------------- BOUND-00: the auditor is audited -----------------------
    def _rule_table_audit():
        ids = [r["id"] for r in RULES]
        dup = sorted({i for i in ids if ids.count(i) > 1})
        nosrc = [r["id"] for r in RULES if r["kind"] != "SKIP" and not r.get("expect")]
        noskip = [r["id"] for r in RULES if r["kind"] == "SKIP" and not r.get("skip")]
        ok = (not dup) and (not nosrc) and (not noskip)
        return ok, ("rules=%d; duplicate_ids=%s; bound_without_source=%s; skip_without_reason=%s"
                    % (len(ids), dup or "none", nosrc or "none", noskip or "none"))
    L.guard("BOUND-00", _rule_table_audit,
            "auditor self-audit: unique rule ids; every bound names the source of its expectation; "
            "every explicit skip names its reason")

    # ---------------- BOUND-02 precondition: required files present ----------
    fig_missing = [f for f in FIG_REQUIRED
                   if not os.path.exists(os.path.join(FIGDIR, f))]
    par_missing = [f for f in PARENT_REQUIRED
                   if not os.path.exists(os.path.join(L.OUT, f))]
    missing = fig_missing + par_missing
    n_fig_found = len(FIG_REQUIRED) - len(fig_missing)
    L.check("BOUND-02a", len(missing) == 0 and len(FIG_FILES) >= len(FIG_REQUIRED),
            "required frozen result files present (declared manifest FIG_REQUIRED + "
            "PARENT_REQUIRED; t06 brief file list)",
            f"figures/data: {n_fig_found}/{len(FIG_REQUIRED)} declared files present "
            f"({len(FIG_FILES)} *.csv on disk); parent: {len(PARENT_REQUIRED)} required "
            f"files; missing={missing if missing else 'none'}",
            "p15 optional: " + ("present" if os.path.exists(os.path.join(L.OUT, P15))
                                else "absent"))

    results = []
    for f in FIG_FILES:
        p = os.path.join(FIGDIR, f)
        try:
            results.append(scan_file(p))
        except Exception as e:
            L.emit("BOUND-01::%s::__file__" % f, "FAIL",
                   "file must parse (frozen result CSV)",
                   f"EXCEPTION {type(e).__name__}: {e}",
                   traceback.format_exc().splitlines()[-1])
    for f in PARENT_FILES:
        p = os.path.join(L.OUT, f)
        if not os.path.exists(p):
            L.emit("BOUND-01::%s::__file__" % f, "WARN",
                   "optional file present (t06 brief: 'if present')", "absent",
                   "skipped - not part of the red-light baseline")
            continue
        try:
            results.append(scan_file(p))
        except Exception as e:
            L.emit("BOUND-01::%s::__file__" % f, "FAIL",
                   "file must parse (frozen result CSV)",
                   f"EXCEPTION {type(e).__name__}: {e}",
                   traceback.format_exc().splitlines()[-1])

    # ---------------- BOUND-01: violations -----------------------------------
    n_hard = n_diag = 0
    for res in results:
        for v in res["viol"]:
            if v["kind"] == "HARD":
                n_hard += 1
                status = "FAIL"
            else:
                n_diag += 1
                status = "WARN"
            tid = "BOUND-01::%s::%s::%s" % (res["file"], v["col"], v["rule"])
            expected = "%s %s | source: %s" % (v["kind"], v["bound"], v["expect"])
            observed = ("n_viol=%d/%d; viol_range=[%.6g, %.6g]; extremal: %s -> %s"
                        % (v["n_viol"], v["n_used"], v["vmin"], v["vmax"],
                           ", ".join(v["ids"]), ", ".join(v["vals"])))
            note = v["justify"] + (" | row filter: %s" % v["cond"] if v["cond"] else "")
            L.emit(tid, status, _san(expected), _san(observed), _san(note))

    # fail-closed coverage: a rule that evaluated ZERO values is not a check
    for res in results:
        for key in sorted(res.get("zero_used", {})):
            kind, cdesc, expect, col = res["zero_used"][key]
            rid = key.split("::")[0]
            status = "FAIL" if kind == "HARD" else "WARN"
            why = ("row predicate %s selected 0 rows" % cdesc if cdesc
                   else "column has no non-NaN values")
            L.emit("BOUND-01::%s::%s::%s" % (res["file"], col, rid), status,
                   _san("%s rule must evaluate at least one value | source: %s"
                        % (kind, expect)),
                   "%s in %s" % (why, res["file"]),
                   "fail-closed: otherwise this rule silently checks nothing while the column "
                   "is reported as scanned")

    L.check("BOUND-01-AGG", n_hard == 0,
            "aggregate: zero HARD definitional-range violations in the scanned set "
            "(expectation: each HARD bound is a mathematical/definitional identity)",
            f"hard_violations={n_hard}; diagnostic_violations={n_diag}; "
            f"files={len(results)}",
            "a FAIL here is the intended red-light signal for the frozen flawed outputs",
            warn_only=False)

    # ---------------- BOUND-02: summary + scanned column list ----------------
    scanned_list = []
    n_numeric = n_scanned = 0
    for res in results:
        n_numeric += len(res["numeric"])
        n_scanned += len(res["scanned"])
        if res["scanned"]:
            scanned_list.append("%s[%s]" % (res["file"], ",".join(res["scanned"])))
    L.emit("BOUND-02", "PASS" if not missing else "FAIL",
           "summary is defined by the t06 brief: counts of hard/diagnostic violations and the "
           "exact list of scanned columns (t06 brief BOUND-02)",
           "files_scanned=%d; numeric_columns=%d; columns_scanned=%d; hard_violations=%d; "
           "diagnostic_violations=%d" % (len(results), n_numeric, n_scanned, n_hard, n_diag),
           "SCANNED COLUMNS: " + " ; ".join(scanned_list))

    # ---------------- BOUND-03: the DR column, explicitly --------------------
    dr_rows = []
    for res in results:
        fname = res["file"]
        p = res["path"]
        try:
            df = pd.read_csv(p)
        except Exception:
            continue
        for col in ("AUROC_DR",):
            if col not in df.columns:
                continue
            s = pd.to_numeric(df[col], errors="coerce")
            gt1 = np.where((s > 1.0 + TOL).to_numpy())[0]
            lt0 = np.where((s < -TOL).to_numpy())[0]
            for i in list(gt1) + list(lt0):
                dr_rows.append((fname, col, _row_label(df, int(i)), float(s.iat[int(i)]),
                                ">1" if i in gt1 else "<0"))
    n_gt1 = sum(1 for r in dr_rows if r[4] == ">1")
    mx = max((r[3] for r in dr_rows), default=float("nan"))
    uniq = {}
    for fname, col, lab, val, _how in dr_rows:
        uniq.setdefault(lab.split("=", 1)[-1], []).append((fname, val))
    detail = "; ".join("%s: max=%.6f (reported in %d files)"
                       % (t, max(v for _, v in vs), len(vs))
                       for t, vs in sorted(uniq.items(), key=lambda kv: -max(v for _, v in kv[1])))
    if dr_rows:
        L.emit("BOUND-03", "WARN",
               "DR-class diagnostic (t06 brief BOUND-03): count of AUROC_DR values outside "
               "[0,1], offending tasks and max",
               "n_gt1_mentions=%d; n_lt0=%d; n_unique_tasks_gt1=%d; max=%.6f; offenders: %s"
               % (n_gt1, len(dr_rows) - n_gt1, len(uniq), mx, detail),
               "justification: unconstrained augmented estimator - a finite-sample DR estimate "
               "may legitimately leave [0,1], so this is WARN, not an implementation error; "
               "the same 4 tasks are restated in fig4_l3_scatter/fig4e/tableS1/p6c")
    else:
        L.emit("BOUND-03", "PASS",
               "DR-class diagnostic (t06 brief BOUND-03): count of AUROC_DR values outside [0,1]",
               "n_gt1=0; n_lt0=0", "justification: unconstrained augmented estimator")

    # ---------------- BOUND-04: coverage of the scan per file ----------------
    for res in results:
        fname = res["file"]
        n_num = len(res["numeric"])
        n_sc = len(res["scanned"])
        n_sk = len(res["skipped"])
        n_un = len(res["unattributed"])
        skip_txt = "; ".join("%s(%s)" % (c, why) for c, why in res["skipped"][:6])
        if n_un:
            L.emit("BOUND-04::%s" % fname, "FAIL",
                   "coverage audit (t06 brief BOUND-04): every numeric column is either "
                   "scanned by a bound rule or explicitly skipped with a reason",
                   "numeric=%d; scanned=%d; skipped=%d; UNATTRIBUTED=%d: %s"
                   % (n_num, n_sc, n_sk, n_un, ",".join(res["unattributed"])),
                   "unattributed numeric columns are a coverage gap - add a bound or a "
                   "justified skip; skipped reasons (first 6): " + (skip_txt or "none"))
        else:
            # 本文件所涉规则的**域集合**（A3.5）
            _doms = sorted({domain_of(r["id"]) for r in RULES
                            if r.get("files") and fname in r["files"]}) or [_DEFAULT_DOMAIN]
            _dom_txt = "+".join(_doms)
            L.emit("BOUND-04::%s" % fname, "PASS",
                   ("coverage audit (t06 brief BOUND-04) [domain: %s]: every numeric column "
                    "is either scanned by a bound rule or explicitly skipped with a reason"
                    % _dom_txt),
                   "numeric=%d; scanned=%d; skipped=%d; UNATTRIBUTED=0"
                   % (n_num, n_sc, n_sk),
                   "skipped columns (first 6): " + (skip_txt or "none"))

    # ---------------- BOUND-05: MMD corroboration (defect (a)) ---------------
    try:
        mmd_rows = []
        for p in [os.path.join(L.OUT, "p6b_label_free_signals_blinded.csv"),
                  os.path.join(FIGDIR, "tableS1_tasks.csv")]:
            if os.path.exists(p):
                df = pd.read_csv(p)
                if "mmd" in df.columns:
                    mmd_rows.append((os.path.basename(p), pd.to_numeric(df["mmd"],
                                                                        errors="coerce")))
        if mmd_rows and len(mmd_rows[0][1].dropna()) > 0:
            name, s = mmd_rows[0]
            s = s.dropna()
            neg = int((s < 0).sum())
            floor = -1.0 / 999.0 - 1.0 / 999.0
            observed = ("%s: n=%d; negative=%d (%.1f%%); median=%.6f; min=%.6f; floor(m=n=1000)"
                        "=%.6f; median/floor=%.1fx"
                        % (name, len(s), neg, 100.0 * neg / len(s), float(s.median()),
                           float(s.min()), floor, abs(float(s.median()) / floor)))
            L.check("BOUND-05", neg <= len(s) / 2,
                    "diagnostic corroboration of defect (a): a valid unbiased MMD^2 estimator "
                    "centred near 0 should not be negative in the large majority of tasks "
                    "(p25 judgment criterion)",
                    observed,
                    "consistent with p6b mmd_rbf halving BOTH within-sample kernel terms "
                    "(pdist upper triangle divided by n(n-1))", warn_only=True)
        else:
            L.emit("BOUND-05", "WARN",
                   "diagnostic corroboration of defect (a): the MMD floor scan input must "
                   "contain at least one finite 'mmd' value (p6b_label_free_signals_blinded.csv "
                   "or tableS1_tasks.csv)",
                   "no usable 'mmd' column found (files_with_mmd=%d; finite_values=%d)"
                   % (len(mmd_rows), len(mmd_rows[0][1].dropna()) if mmd_rows else 0),
                   "non-coverage reported explicitly (this branch previously emitted no RESULT "
                   "line at all)")
    except Exception as e:
        L.emit("BOUND-05", "FAIL", "MMD corroboration diagnostic must run",
               "EXCEPTION %s: %s" % (type(e).__name__, e), "")

    # ---------------- ESS-04: ESS 上界锚定（CS-09 产物级，B 组） ----------------
    # 依据：FREEZE_pre_correction_20261005/CS09_ESS_acceptance.md §3 B 组
    #   「上界锚定：w_ESS <= n_src_test；现状 170/202 违约」。
    # ★ 为什么不只靠 R14：R14 是**只报违规**的规则（通过时不发射 RESULT 行），
    #   所以它不能充当 manifest 里的 test_id —— 会造出 manifest 自己警告过的
    #   「登记了但从不发射」的空条目。本判据**恒发射**，把「ESS 不得超过其权重域的
    #   样本量」做成一条可登记、可销账的门。
    # ★ 非空转：双向实测 —— 修正前产物 170/202 违约（max 6.436）判红；
    #   重跑后 0/186（max 0.941）判绿。若哪天它恒绿，说明产物已换掉，需重新取 RED 基线。
    try:
        _p6c = os.path.join(L.OUT, "p6c_L3_blinded.csv")
        _exp = ("Cauchy-Schwarz 定义性上界：建在 N 个权重上的 ESS 不可能超过 N。"
                "CS-09 后 p6c 的 ESS 建在 w_te（source internal test）上 ⇒ N = n_src_test")
        if not os.path.exists(_p6c):
            L.emit("ESS-04", "FAIL", _exp, "p6c_L3_blinded.csv 不存在 —— 判据无靶子",
                   "缺文件判红而不是跳过：跳过会让本路径静默失去 ESS 覆盖")
        else:
            _d = pd.read_csv(_p6c)
            if "w_ESS" not in _d.columns or "n_src_test" not in _d.columns:
                L.emit("ESS-04", "FAIL", _exp,
                       "p6c_L3_blinded.csv 缺列 w_ESS / n_src_test（有 w_ESS=%s, n_src_test=%s）"
                       % ("w_ESS" in _d.columns, "n_src_test" in _d.columns),
                       "列缺失判红，避免判据在无靶子时恒绿")
            else:
                _e = pd.to_numeric(_d["w_ESS"], errors="coerce").to_numpy(dtype=float)
                _n = pd.to_numeric(_d["n_src_test"], errors="coerce").to_numpy(dtype=float)
                _ok = np.isfinite(_e) & np.isfinite(_n)
                _viol = int(np.sum(_ok & (_e > _n + 1e-9)))
                _ratio = float(np.max(_e[_ok] / _n[_ok])) if _ok.any() else float("nan")
                L.check("ESS-04", _viol == 0, _exp,
                        "n=%d（可判 %d）w_ESS > n_src_test 的项 = %d；max(ESS/n_src_test) = %.3f"
                        % (len(_d), int(_ok.sum()), _viol, _ratio),
                        "双向实测：修正前 170/202（max 6.436）判红 → 重跑后 0/186（max 0.941）判绿"
                        " ⇒ 判据有鉴别力，非空转")
    except Exception as e:                                              # noqa: BLE001
        L.emit("ESS-04", "FAIL",
               "Cauchy-Schwarz 定义性上界 w_ESS <= n_src_test（CS-09 产物级，B 组）",
               "EXCEPTION %s: %s" % (type(e).__name__, e),
               "判据必须运行；异常不得被当成通过")

    # ---------------- ESS-05: 下游依赖对账（CS-09 产物级，B 组） ----------------
    # 依据：FREEZE_pre_correction_20261005/CS09_ESS_acceptance.md §3 B 组「下游依赖检查」：
    #   · fig4e_dr_out_of_range.csv / fig4e_ess_association.csv 必须在 CS-05 **与** CS-09
    #     都落地后**重新生成**；
    #   · fig4e_ess_association.csv 的 `n_dr_gt1` 必须等于**同期** p6c_L3_blinded.csv 中
    #     `AUROC_DR > 1` 的行数；
    #   · `ess_min/median/max` 必须等于同期 `w_ESS` 列的 min/median/max。
    # ★ 为什么这条重要：它把「sidecar 与它声称依据的那期 p6c」绑在一起。
    #   2026-10-06 实测前，fig4e_ess_association.csv 停在修正前的一期
    #   （n_dr_gt1=4、ess_min=96.65），而同期 p6c 已是 AUROC_DR 无越界 —— 正是本条要抓的。
    # ★ 非空转：判据在**旧 sidecar + 新 p6c** 的组合上必然失败（实测差 4 与 85.5）。
    try:
        _f4e = os.path.join(L.PROJ, "figures", "data", "fig4e_ess_association.csv")
        _f4o = os.path.join(L.PROJ, "figures", "data", "fig4e_dr_out_of_range.csv")
        _p6c2 = os.path.join(L.OUT, "p6c_L3_blinded.csv")
        _exp5 = ("CS-09 §3 B 组下游依赖：fig4e_ess_association.csv 的 n_dr_gt1 与 "
                 "ess_min/median/max 必须等于**同期** p6c_L3_blinded.csv 的对应实测值")
        _miss = [p for p in (_f4e, _f4o, _p6c2) if not os.path.exists(p)]
        if _miss:
            L.emit("ESS-05", "FAIL", _exp5,
                   "缺文件：" + ", ".join(os.path.basename(p) for p in _miss),
                   "缺文件判红而不是跳过：跳过会让本判据静默失去靶子")
        else:
            _e = pd.read_csv(_f4e).iloc[0]
            _l = pd.read_csv(_p6c2)
            _o = pd.read_csv(_f4o)
            _n_dr = int((pd.to_numeric(_l["AUROC_DR"], errors="coerce") > 1).sum())
            _w = pd.to_numeric(_l["w_ESS"], errors="coerce").astype(float)
            _chk = [
                ("n_dr_gt1", float(_e["n_dr_gt1"]), float(_n_dr)),
                ("ess_min", float(_e["ess_min"]), float(_w.min())),
                ("ess_median", float(_e["ess_median"]), float(_w.median())),
                ("ess_max", float(_e["ess_max"]), float(_w.max())),
                ("dr_out_of_range 行数", float(len(_o)), float(_n_dr)),
            ]
            _bad = [(k, a, b) for k, a, b in _chk if abs(a - b) > 1e-6]
            L.check("ESS-05", not _bad, _exp5,
                    "；".join("%s csv=%g vs 同期实测=%g" % (k, a, b) for k, a, b in _chk)
                    + ("" if not _bad else "  ★ 不符：" + str(_bad)),
                    "双向实测：旧 sidecar（n_dr_gt1=4、ess_min=96.65）配新 p6c 必然 FAIL；"
                    "sidecar 与 p6c 同期后 5/5 项相符 ⇒ PASS")
    except Exception as e:                                              # noqa: BLE001
        L.emit("ESS-05", "FAIL",
               "CS-09 §3 B 组下游依赖：fig4e 两个 sidecar 必须与同期 p6c 对账",
               "EXCEPTION %s: %s" % (type(e).__name__, e),
               "判据必须运行；异常不得被当成通过")

    # ---------------- BOUND-06: rule non-coverage, explicitly ----------------
    rule_hits = set()
    for res in results:
        rule_hits.update(res.get("rule_hits", ()))
    non_covered = [r["id"] for r in RULES
                   if r["kind"] != "SKIP" and r["id"] not in rule_hits]
    L.check("BOUND-06", not non_covered,
            "docstring contract: rule non-coverage is reported explicitly - every rule id "
            "that matched zero numeric columns in the scanned set is listed",
            "rules_with_zero_matched_columns=%d: %s"
            % (len(non_covered), ", ".join(non_covered) if non_covered else "none"),
            "zero matched columns means the rule checked nothing in this run; the expected "
            "steady-state list is exactly R41_pvalue - no p-value column exists in the frozen set",
            warn_only=True)

    print("\n[t06 bounds scan done]  hard=%d diag=%d; exit 0 by contract" % (n_hard, n_diag))


try:
    main()
except Exception:
    L.emit("BOUND-00", "FAIL", "module must run to completion (t06 contract)",
           "EXCEPTION " + traceback.format_exc().splitlines()[-1],
           "unhandled top-level error")
sys.exit(0)
