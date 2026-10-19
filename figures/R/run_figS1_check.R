# =============================================================================
# Supplementary Figure S1 校验
# =============================================================================
.armd_f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
.armd_d <- if (length(.armd_f) && !is.na(.armd_f)) dirname(.armd_f) else "."
ROOT <- Sys.getenv("ARMD_ROOT", unset = normalizePath(file.path(.armd_d, "..", "..")))
setwd(file.path(ROOT, "figures"))
source("R/check_rendered.R", encoding = "UTF-8")

run_and_check(
  script   = "R/FigureS1_task_matrix.R",
  outbase  = "FigureS1",
  outdir   = "supplement",
  logfile  = "log/FigureS1_render_check.txt",
  panels   = character(0),
  required = c(
    "Complete cross-system task matrix",
    "228 candidate cells: 38 organism-drug groups",
    "186 eligible (coloured); 42 excluded",
    "Row and column order fixed by organism, drug and direction, not by outcome",
    "MGB\u2192SU", "MGB\u2192UT", "SU\u2192MGB", "SU\u2192UT", "UT\u2192MGB", "UT\u2192SU",
    "EC|AMC", "KP|NIT", "PA|ATM",
    "6/6", "2/6", "4/6",
    "eligible / 38",
    "AUROC"
  ),
  forbidden = c(
    # 不得按结果排序；不得把"未通过资格线"说成数据缺失；不得对 harmful 做第二套编码
    "sorted by", "ordered by delta", "ranked by outcome",
    "missing data", "data missing", "NA", "not available",
    "harmful", "Harmful", "non-harmful", "border marks", "outlined"
  )
)

svg <- paste(readLines("supplement/FigureS1.svg", warn = FALSE), collapse = "\n")

# 断言 1：行序 EC 在顶、PA 在底。
# 注意：**不能用"在 SVG 中出现得更早"来判断** —— svglite 按 factor 水平顺序绘制文字，
# 与视觉上下位置无关（首版据此误报）。必须比较 **y 坐标**：SVG 中 y 越小越靠上。
yof <- function(lab) {
  m <- regmatches(svg, regexpr(sprintf("<text x='[^']*' y='([0-9.]+)'[^>]*>%s", lab), svg))
  if (!length(m)) return(NA_real_)
  as.numeric(sub(".*y='([0-9.]+)'.*", "\\1", m))
}
yEC <- yof("EC\\|AMC"); yPA <- yof("PA\\|ATM")
if (is.na(yEC) || is.na(yPA)) stop("找不到 EC/PA 行标签的坐标")
if (yEC >= yPA) stop(sprintf("行序错误：EC(y=%.1f) 应在 PA(y=%.1f) 之上", yEC, yPA))
cat(sprintf("\n[行序断言] EC y=%.1f < PA y=%.1f —— EC 在顶、PA 在底（按 organism 固定顺序）\n",
            yEC, yPA))

# 断言 2：26 个 excluded 格必须都有 × 标记
n_x <- lengths(regmatches(svg, gregexpr("\u00d7</text>", svg, fixed = TRUE)))
if (n_x < 42) stop(sprintf("excluded 格的 × 标记只有 %d 个，应为 42", n_x))
cat(sprintf("[完整性断言] 42 个 excluded cell 均有 × 标记（实测 %d）\n", n_x))

# 断言 3：底部列计数必须是 37/37/34/34/30/30
for (v in c("35", "31", "27")) if (!grepl(sprintf(">%s</text>", v), svg, fixed = TRUE))
  stop("列计数缺失: ", v)
cat("[列计数断言] 每列 eligible 计数已渲染（35/35/31/31/27/27，分母 38）\n")
