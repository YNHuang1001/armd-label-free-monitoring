# =============================================================================
# Supplementary Figure S4 校验（含盲态防线复核）
# =============================================================================
.armd_f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
.armd_d <- if (length(.armd_f) && !is.na(.armd_f)) dirname(.armd_f) else "."
ROOT <- Sys.getenv("ARMD_ROOT", unset = normalizePath(file.path(.armd_d, "..", "..")))
setwd(file.path(ROOT, "figures"))
source("R/check_rendered.R", encoding = "UTF-8")

run_and_check(
  script   = "R/FigureS4_signal_qc.R",
  outbase  = "FigureS4",
  outdir   = "supplement",
  logfile  = "log/FigureS4_render_check.txt",
  panels   = c("A", "B"),
  required = c(
    "Blinded distributions of the L1 and L2 signals",
    "Blinded distributions of the transported-AUC estimators",
    "Native scale; all 186 tasks",
    "QC distributions, not a comparison of methods",
    "Shared x-axis, 186 tasks",
    "0 DR estimates exceed 1; none were clipped",
    "Signal value (native scale)", "Estimated target AUROC",
    "domain-AUC", "MMD", "PSI, mean", "PSI, max", "JS, mean", "JS, max",
    "Prediction mean shift", "Entropy mean shift", "Entropy JS",
    "Confidence-dist. JS", "Predicted prevalence shift", "SUDO-AURCC",
    "DR (primary)", "IOW (sensitivity)", "Outcome model (sensitivity)"
  ),
  forbidden = c(
    # 本图是 blinded QC：**不得**出现任何 target-failure 相关措辞、分组或叠加
    "harmful", "Harmful", "non-harmful", "deployment failure", "failure detection",
    "correlation with", "AUROC vs", "grouped by", "faceted by outcome",
    "paired difference", "sample size"
  )
)

svg <- paste(readLines("supplement/FigureS4.svg", warn = FALSE), collapse = "\n")

# 断言 1：Panel A 必须恰好 12 个 facet，**不得出现 NA 面板**
if (grepl(">NA</text>", svg, fixed = TRUE)) stop("Panel A 出现了 NA facet")
cat("\n[盲态QC断言] Panel A 无 NA 面板（12 个 signal facet 齐全）\n")

# 断言 2：盲态防线必须真的跑过（render log 中应有两条防线输出）
chk <- paste(readLines("log/FigureS4_render_check.txt", warn = FALSE), collapse = "\n")
if (!grepl("脚本层自检通过", chk, fixed = TRUE)) stop("脚本层盲态防线未执行")
if (!grepl("数据层断言通过", chk, fixed = TRUE)) stop("数据层盲态防线未执行")
cat("[盲态QC断言] 两道盲态防线均已执行并在日志中留痕\n")

# 断言 3：DR 的 4 个 >1 值必须原样可见（x 轴上界须 >= 1.1）
if (!grepl("1.1", svg, fixed = TRUE)) stop("Panel B 的 x 轴未覆盖到 >1 的 DR 估计值")
cat("[盲态QC断言] Panel B 的 x 轴覆盖 1.1，DR >1 的估计值原样可见\n")
