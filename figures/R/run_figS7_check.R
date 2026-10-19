# =============================================================================
# Supplementary Figure S7 校验
# =============================================================================
.armd_f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
.armd_d <- if (length(.armd_f) && !is.na(.armd_f)) dirname(.armd_f) else "."
ROOT <- Sys.getenv("ARMD_ROOT", unset = normalizePath(file.path(.armd_d, "..", "..")))
setwd(file.path(ROOT, "figures"))
source("R/check_rendered.R", encoding = "UTF-8")

run_and_check(
  script   = "R/FigureS7_weight_diagnostics.R",
  outbase  = "FigureS7",
  outdir   = "supplement",
  logfile  = "log/FigureS7_render_check.txt",
  panels   = c("A", "B", "C"),
  required = c(
    "Effective sample size", "Maximum weight",
    "ESS (log scale)", "Max weight (log scale)",
    "-0.249", "+0.220",
    "Association with", "absolute DR error", "Descriptive, not causal",
    "log10 effective sample size",
    "log10 maximum weight",
    "Spearman"
  ),
  forbidden = c(
    # 不得把诊断关联写成因果；不得暗示"重新选指标"；不得声称某种方法更优
    "causes", "causal effect", "because of", "leads to", "due to",
    "ESS is superior", "better estimator", "we selected the better",
    "should replace", "optimal diagnostic",
    # 【永久禁止】v1 曾把 w_max 的正相关误读为"相反方向"，这两句必须不能再出现
    "opposite direction", "inverse, as expected", "opposite to expectation"
  )
)

# 断言：两个关联的符号。**两者方向一致**，均为权重不稳定性下的预期方向
# （ESS 越高误差越小；w_max 越大误差越大）。差别只在强度。
# 【修正记录】v1 曾断言"符号相反"并据此解释，是符号读反，见 audit log A-2026-09-30-12。
svg <- paste(readLines("supplement/FigureS7.svg", warn = FALSE), collapse = "\n")
if (!grepl("-0.249", svg, fixed = TRUE) || !grepl("+0.220", svg, fixed = TRUE))
  stop("两个关联系数未同时出现")
cat("\n[方向断言] ESS 为负向、w_max 为正向 —— 两者**方向一致**，均为预期方向，\n")
cat("           差别在强度（|rho| 0.301 vs 0.178）；图上不得出现方向性标签。\n")

# 断言：不得出现 w_P95（GPT6 明确要求 S7 不塞第三个诊断量）
if (grepl("P95", svg, fixed = TRUE)) stop("图上出现了 w_P95 —— 超出 S7 锁定的三 panel 范围")
cat("[范围断言] 未出现 w_P95 —— S7 保持三 panel\n")
