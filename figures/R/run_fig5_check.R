# =============================================================================
# Figure 5 校验
# =============================================================================
.armd_f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
.armd_d <- if (length(.armd_f) && !is.na(.armd_f)) dirname(.armd_f) else "."
ROOT <- Sys.getenv("ARMD_ROOT", unset = normalizePath(file.path(.armd_d, "..", "..")))
setwd(file.path(ROOT, "figures"))
source("R/check_rendered.R", encoding = "UTF-8")

run_and_check(
  script  = "R/Figure5_samplesize.R",
  outbase = "Figure5",
  panels  = c("A", "B"),
  required = c(
    "Change from full cohort to n = 100",
    "Variation across all five target sizes",
    "Paired difference within each repeat",
    "Range of harmful-detection AUROC over 100 / 250 / 500 / 1000 / full",
    "AUROC range (max", "across target sizes)",
    "L1  input shift", "L2  prediction shift", "L3  performance estimation",
    "domain-AUC", "MMD", "PSI, mean", "PSI, max", "JS, mean", "JS, max",
    "Prediction mean shift", "Entropy mean shift", "Entropy JS",
    "Confidence-dist. JS", "Predicted prevalence shift", "SUDO-AURCC",
    "DR   (primary)", "IOW   (sensitivity)", "Outcome model   (sensitivity)"
  ),
  forbidden = c(
    # 不得出现"稳定即够用"的措辞；不得把唯一名义显著做成主故事
    "sufficient", "Sufficient", "adequate", "Adequate", "enough",
    "sample size is not a problem",
    "significant*", "p < 0.05", "**",
    "Ranking", "ranked by"
  )
)

# 断言 1：层顺序 L1 -> L2 -> L3
svg <- paste(readLines("out/Figure5.svg", warn = FALSE), collapse = "\n")
i1 <- regexpr("L1  input shift", svg, fixed = TRUE)[1]
i2 <- regexpr("L2  prediction shift", svg, fixed = TRUE)[1]
i3 <- regexpr("L3  performance estimation", svg, fixed = TRUE)[1]
if (any(c(i1, i2, i3) < 0)) stop("找不到 L1/L2/L3 层标签")
if (!(i1 < i2 && i2 < i3)) stop("层顺序不是 L1 -> L2 -> L3")
cat(sprintf("\n[顺序断言] 层顺序 L1(%d) < L2(%d) < L3(%d) —— OK\n", i1, i2, i3))

# 断言 2：Panel A 必须按 protocol 顺序（domain-AUC 在 SUDO-AURCC 之前）
a_first <- regexpr("domain-AUC", svg, fixed = TRUE)[1]
a_sudo  <- regexpr("SUDO-AURCC", svg, fixed = TRUE)[1]
if (a_first < 0 || a_sudo < 0 || a_first > a_sudo)
  stop("未按 protocol 顺序排列")
cat("[顺序断言] 按 protocol 顺序（非按效应量排序）—— OK\n")

# 断言 3：across-N 极差最大值必须为 0.261（与 Python 层一致）
# 【2026-10-06 修正】原为 0.038（旧宇宙；当时最大者是 predicted-prevalence shift）。
#   修正后 AUROC_DR 的跨档极差为 0.2608，为全表最大，也是本图现在要传达的主要发现。
if (!grepl("0.261", svg, fixed = TRUE))
  stop("Panel B 未出现预期的最大极差 0.261 —— 数据或标注可能不符")
cat("[数值断言] Panel B 最大 across-N 极差 0.261 —— OK\n")
