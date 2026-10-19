# =============================================================================
# Figure 3 校验
# =============================================================================
.armd_f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
.armd_d <- if (length(.armd_f) && !is.na(.armd_f)) dirname(.armd_f) else "."
ROOT <- Sys.getenv("ARMD_ROOT", unset = normalizePath(file.path(.armd_d, "..", "..")))
setwd(file.path(ROOT, "figures"))
source("R/check_rendered.R", encoding = "UTF-8")

run_and_check(
  script  = "R/Figure3_grouped.R",
  outbase = "Figure3",
  panels  = c("A", "B", "C"),
  required = c(
    "Grouped generalization", "leave-one-organism-drug-out",
    "26/37 evaluable folds; 134/186 tasks",
    "L1  input shift", "L2  prediction shift", "L3  performance estimation",
    "domain-AUC", "MMD", "PSI, mean", "PSI, max", "JS, mean", "JS, max",
    "Prediction mean shift", "Entropy mean shift", "Entropy JS",
    "Confidence-dist. JS", "Predicted prevalence shift", "SUDO-AURCC",
    "DR   (primary)", "IOW   (sensitivity)", "Outcome model   (sensitivity)",
    "incl. 0.5", "above 0.5", "below 0.5",
    "External stress test", "held-out-organism folds",
    "min\u2013max range", "descriptive, not a CI",
    "EC", "KP", "PA",
    "Aggregation changes", "186 tasks", "134 tasks", "macro LODO",
    "task-subset effect", "aggregation effect", "0.332", "0.408", "0.325"
  ),
  forbidden = c(
    # 明确禁止的措辞与做法
    "worse than random",          # 会被误读为"无信息"（MMD 存在事后反向秩信息）
    "Waterfall", "waterfall",     # 易被误读为因果分解
    "95% CI: 3-fold",             # 3 折不得冒称 CI
    "confidence interval (LOO-organism)",
    "SUDO is the best", "best method", "Ranking"
  )
)

# 顺序断言 1：15 个方法在 B/C 之前的 A 面板中按 protocol 顺序（L1 在前）
svg <- paste(readLines("out/Figure3.svg", warn = FALSE), collapse = "\n")
i1 <- regexpr("L1  input shift", svg, fixed = TRUE)[1]
i2 <- regexpr("L2  prediction shift", svg, fixed = TRUE)[1]
i3 <- regexpr("L3  performance estimation", svg, fixed = TRUE)[1]
if (any(c(i1, i2, i3) < 0)) stop("找不到 L1/L2/L3 层标签")
if (!(i1 < i2 && i2 < i3)) stop("层顺序不是 L1 -> L2 -> L3")
cat(sprintf("\n[顺序断言] 层顺序 L1(%d) < L2(%d) < L3(%d) —— OK\n", i1, i2, i3))

# 顺序断言 2：不得出现"按性能排序"——A 面板首个方法是 domain-AUC（protocol 首位）
a_first <- regexpr("domain-AUC", svg, fixed = TRUE)[1]
a_sudo  <- regexpr("SUDO-AURCC", svg, fixed = TRUE)[1]
if (a_first < 0 || a_sudo < 0 || a_first > a_sudo)
  stop("A 面板未按 protocol 顺序排列（domain-AUC 应在 SUDO-AURCC 之前）")
cat("[顺序断言] A 面板按 protocol 顺序（非按性能排序）—— OK\n")
