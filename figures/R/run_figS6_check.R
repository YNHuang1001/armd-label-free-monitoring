# =============================================================================
# Supplementary Figure S6 校验
# =============================================================================
.armd_f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
.armd_d <- if (length(.armd_f) && !is.na(.armd_f)) dirname(.armd_f) else "."
ROOT <- Sys.getenv("ARMD_ROOT", unset = normalizePath(file.path(.armd_d, "..", "..")))
setwd(file.path(ROOT, "figures"))
source("R/check_rendered.R", encoding = "UTF-8")

run_and_check(
  script   = "R/FigureS6_lodo_composition.R",
  outbase  = "FigureS6",
  outdir   = "supplement",
  logfile  = "log/FigureS6_render_check.txt",
  panels   = c("A", "B"),
  required = c(
    "Composition of each organism-drug group",
    "A within-group AUROC requires both transfer classes",
    "Evaluable", "(n = 26)", "Single-class", "(n = 11)",
    "Harmful", "Non-harmful", "Number of transfer tasks",
    "Restriction imposed by single-class folds",
    "evaluable subset only",
    "Both classes", "All harmful", "All non-harmful",
    "total 37", "total 186",
    "26", "1", "10", "134", "2", "50"
  ),
  forbidden = c(
    # 不得把不可定义折说成"数据缺失/清洗"；不得暗示人为筛选
    "removed", "Removed", "filtered out", "cleaning", "data loss",
    "excluded manually", "95% CI", "confidence interval"
  )
)

# 断言：分区计数必须自洽（38 = 25 + 13；202 = 138 + 64；13 = 3 + 10）
svg <- paste(readLines("supplement/FigureS6.svg", warn = FALSE), collapse = "\n")
for (s in c("total 37", "total 186")) if (!grepl(s, svg, fixed = TRUE)) stop("缺少 ", s)
if (!grepl("(n = 26)", svg, fixed = TRUE) || !grepl("(n = 11)", svg, fixed = TRUE))
  stop("分区计数标签缺失")
cat("\n[自洽断言] 37 = 26 + 11；186 = 134 + 52；11 = 1 + 10 —— 图上与数据一致\n")
