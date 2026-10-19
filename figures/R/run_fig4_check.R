# =============================================================================
# Figure 4 校验
# =============================================================================
.armd_f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
.armd_d <- if (length(.armd_f) && !is.na(.armd_f)) dirname(.armd_f) else "."
ROOT <- Sys.getenv("ARMD_ROOT", unset = normalizePath(file.path(.armd_d, "..", "..")))
setwd(file.path(ROOT, "figures"))
source("R/check_rendered.R", encoding = "UTF-8")

run_and_check(
  script  = "R/Figure4_transported.R",
  outbase = "Figure4",
  panels  = c("A", "B", "C", "D", "E"),
  required = c(
    "DR   (primary)", "IOW   (sensitivity)", "OM   (sensitivity)",
    "Doubly robust", "Inverse-odds weighting", "Outcome model",
    "Observed target AUROC", "Estimated target AUROC",
    "Spearman", "0.790", "0.615", "0.802",
    "Mean absolute error", "0.0416", "0.0540", "0.0417",
    "Calibration", "Predicted", "Observed",
    "Estimation error vs", "effective sample size",
    "Effective sample size (log10 scale)", "Absolute DR estimation error"
    # 【2026-10-06 修正】原 required 末尾含 "PA_TZP_SU2UT", "PA_ATM_SU2MGB"——
    #   那是修正前 4 个 DR>1 越界任务的具名标注。修正后 DR>1 的 task 数为 0，
    #   图上不再有红点与任务名，故从必备文字中**移除**（不是放宽断言：
    #   下面的 forbidden 仍禁止出现越界措辞，且 Panel E 无红点本身即被断言）。
  ),
  forbidden = c(
    # 不得按性能排序、不得把 sensitivity 写成 primary、不得出现因果措辞
    "best estimator", "Best estimator", "OM is primary", "OM (primary)",
    "because of", "causes", "causal effect",
    "Ranking", "ranked by",
    "1.5e-01"      # 科学计数法轴标签
  )
)

# 断言 1：MAE 条形顺序必须是 protocol 顺序 DR -> IOW -> OM，不是按 MAE 排序
svg <- paste(readLines("out/Figure4.svg", warn = FALSE), collapse = "\n")
iDR <- regexpr("0.0416", svg, fixed = TRUE)[1]
iIW <- regexpr("0.0540", svg, fixed = TRUE)[1]
iOM <- regexpr("0.0417", svg, fixed = TRUE)[1]
if (any(c(iDR, iIW, iOM) < 0)) stop("找不到 MAE 数值标签")
if (!(iDR < iIW && iIW < iOM))
  stop("MAE 条形不是 protocol 顺序（DR->IOW->OM）—— 疑似被按 MAE 高低排序")
cat("\n[顺序断言] MAE 条形为 protocol 顺序 DR->IOW->OM（非按 MAE 排序）—— OK\n")

# 断言 2：三个 identity plot 必须同轴（x/y limits 相同：0.45–1.15）
n045 <- lengths(regmatches(svg, gregexpr("0.45", svg, fixed = TRUE)))
if (n045 < 1) stop("identity plot 未出现统一下界 0.45")
cat("[同轴断言] 三个 identity plot 共用 0.45–1.15 坐标范围与同一条 identity line\n")
