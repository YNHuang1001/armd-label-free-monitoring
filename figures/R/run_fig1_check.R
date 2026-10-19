# =============================================================================
# Figure 1 回归测试：用 check_rendered.R 校验
# =============================================================================
.armd_f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
.armd_d <- if (length(.armd_f) && !is.na(.armd_f)) dirname(.armd_f) else "."
ROOT <- Sys.getenv("ARMD_ROOT", unset = normalizePath(file.path(.armd_d, "..", "..")))
setwd(file.path(ROOT, "figures"))
source("R/check_rendered.R", encoding = "UTF-8")

run_and_check(
  script  = "R/Figure1_v2.R",
  outbase = "Figure1_v2",
  panels  = c("A", "B", "C", "D", "E"),
  required = c(
    # Panel A
    "SOURCE", "TARGET",
    "All 6 directed cross-system pairs",
    "37 organism-drug groups",
    "input shift", "prediction shift", "performance estimation",
    # Panel B
    "Task universe", "Candidate", "Eligible",
    # Panel C
    "Cross-system degradation", "P10", "P50", "P90",
    "harmful", "57", "43", "0.151",
    # Panel D
    "Directional heterogeneity", "Ordered by precomputed median",
    # Panel E
    "Degradation patterns are similar", "across classifiers",
    "Same P1-defined 186-task benchmark universe",
    "History", "LGBM", "Core",
    "57 harmful", "63 harmful", "71 harmful", "70 harmful",
    # 三层正序（L1 必须出现在 L3 之前）
    "L1", "L2", "L3"
  ),
  forbidden = c(
    # v1 的旧标题/旧措辞，不得残留
    "Failure is not classifier-specific",
    "6 directed transfers per organism-drug",
    "ordered by median",          # v1 的小写措辞
    "med ",                       # v1 的缩写注释
    "No mechanistic"              # 大小写变体，v2 用小写 no
  )
)

# 额外的顺序断言：L1 的 y 位置必须高于 L3（在 SVG 中即先出现）
svg <- readLines("out/Figure1_v2.svg", warn = FALSE)
s <- paste(svg, collapse = "\n")
i1 <- regexpr(">L1</text>", s, fixed = TRUE)[1]
i3 <- regexpr(">L3</text>", s, fixed = TRUE)[1]
if (i1 < 0 || i3 < 0) stop("找不到 L1/L3 标签")
if (i1 > i3) stop("L1 在 SVG 中出现在 L3 之后 —— 三层顺序仍为反序")
cat(sprintf("\n[顺序断言] L1@%d < L3@%d  ->  L1..L3 正序 OK\n", i1, i3))
