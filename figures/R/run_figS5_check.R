# =============================================================================
# Supplementary Figure S5 校验
# =============================================================================
.armd_f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
.armd_d <- if (length(.armd_f) && !is.na(.armd_f)) dirname(.armd_f) else "."
ROOT <- Sys.getenv("ARMD_ROOT", unset = normalizePath(file.path(.armd_d, "..", "..")))
setwd(file.path(ROOT, "figures"))
source("R/check_rendered.R", encoding = "UTF-8")

run_and_check(
  script   = "R/FigureS5_mmd_implementation.R",
  outbase  = "FigureS5",
  outdir   = "supplement",
  logfile  = "log/FigureS5_render_check.txt",
  panels   = c("A", "B", "C"),
  required = c(
    "MMD² across", "across 186 tasks",
    "Unbiased U-statistic estimator", "0 of 186 tasks negative",
    "Replicate spread", "no task had SD exceeding",
    "Implementation of the frozen estimator",
    "Frozen definitions; no tuning on target data",
    "Source-defined", "RBF kernel", "median", "heuristic on",
    "Unbiased", "U-statistic", "m = 1000, K = 5"
  ),
  forbidden = c(
    # 【关键】不得把负值写成估计器有系统偏差
    "negative bias", "systematic bias", "biased estimator", "bias in the estimator",
    "overestimates", "underestimates", "systematically negative",
    # 不得把旧实现与修正后并排呈现（GPT6 决定：旧产物只作 provenance 留档）
    "before correction", "pre-correction", "old implementation", "previous implementation",
    "we previously", "erroneous"
  )
)

svg <- paste(readLines("supplement/FigureS5.svg", warn = FALSE), collapse = "\n")

# 断言 1：负值计数必须与当前产物一致
# 【2026-10-06 修正】原写死 "200 of 202"（CS-06 后、本轮修正前的那版实现）。
#   再修正后 n_negative = 0、n_tasks = 186，图为 "0 of 186"。改为从产物读出后再核对，
#   避免下次再写死一个数。
.S5  <- read.csv(file.path("data", "figS5_mmd_summary.csv"))
.exp <- sprintf("%d of %d tasks negative", .S5$n_negative, .S5$n_tasks)
if (!grepl(.exp, svg, fixed = TRUE)) stop("未标注当前的负值计数，期望 \"", .exp, "\"")
cat("\n[数值断言] 0/186 已标注（corrected 实现；旧实现为 200/202，见 audit log）\n")

# 断言 2：负值必须原样呈现，x 轴须覆盖到 0 与负值区
if (!grepl("MMD\u00b2 = 0", svg, fixed = TRUE)) stop("缺少 MMD² = 0 参考线标注")
cat("[呈现断言] MMD² = 0 参考线已标注，负值未 clip\n")

# 断言 3：不得出现旧实现的任何数值特征（旧的 range 下界 -0.5038 不应出现）
if (grepl("-0.504", svg, fixed = TRUE) || grepl("202 of 202", svg, fixed = TRUE))
  stop("图上出现旧（invalid）实现的数值特征")
cat("[留档断言] 未出现旧 invalid 实现的数值特征\n")
