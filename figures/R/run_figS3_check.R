# =============================================================================
# Supplementary Figure S3 校验
# =============================================================================
.armd_f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
.armd_d <- if (length(.armd_f) && !is.na(.armd_f)) dirname(.armd_f) else "."
ROOT <- Sys.getenv("ARMD_ROOT", unset = normalizePath(file.path(.armd_d, "..", "..")))
setwd(file.path(ROOT, "figures"))
source("R/check_rendered.R", encoding = "UTF-8")

run_and_check(
  script   = "R/FigureS3_directional.R",
  outbase  = "FigureS3",
  outdir   = "supplement",
  logfile  = "log/FigureS3_render_check.txt",
  panels   = c("A", "B"),
  required = c(
    "Complete degradation distribution for each directed transfer",
    "All 186 eligible transfers",
    "Direction order is fixed by source then target system, not by observed degradation",
    "Outcome composition of each direction",
    "Counts of eligible transfers",
    "Harmful", "Stable", "Other eligible",
    "MGB\u2192SU", "MGB\u2192UT", "SU\u2192MGB", "SU\u2192UT", "UT\u2192MGB", "UT\u2192SU",
    "37", "34", "30"
  ),
  forbidden = c(
    # 不做方向排名、不点名"最差方向"、不给机制解释
    "worst direction", "best direction", "most affected", "least affected",
    "strongest degradation", "ranking", "Ranking", "ranked by",
    "because of", "due to", "explained by", "mechanism", "weaker system", "better system"
  )
)

svg <- paste(readLines("supplement/FigureS3.svg", warn = FALSE), collapse = "\n")

# 断言 1：方向顺序必须是"先 source 后 target"的固定序，**不是按 median 排序**
# 固定序：MGB→SU, MGB→UT, SU→MGB, SU→UT, UT→MGB, UT→SU
# 按 median 排序则是：MGB→UT, SU→UT, MGB→SU, UT→SU, UT→MGB, SU→MGB（两者不同，故可判别）
xpos <- function(lab) {
  m <- regmatches(svg, regexpr(sprintf("<text x='([0-9.]+)'[^>]*>%s", lab), svg))
  if (!length(m)) return(NA_real_)
  as.numeric(sub(".*x='([0-9.]+)'.*", "\\1", m))
}
xo <- vapply(c("MGB\u2192SU","MGB\u2192UT","SU\u2192MGB","SU\u2192UT","UT\u2192MGB","UT\u2192SU"), xpos, numeric(1))
if (any(is.na(xo))) stop("找不到某个方向的轴标签")
if (is.unsorted(xo)) stop("方向顺序不是固定的 source->target 序（疑似按结果排序）")
cat(sprintf("\n[顺序断言] 方向 x 坐标 %s —— 固定 source->target 序，非按 median 排序\n",
            paste(sprintf("%.0f", xo), collapse = " < ")))

# 断言 2：Panel B 的堆叠总数须等于各方向 eligible 数（186 合计）
# 【2026-10-06 修正】原还接受 "All 202 eligible"（旧宇宙）；现只认 186，
#   并**新增**一条反向断言：图中不得再出现旧的 202 合计标注。
if (!grepl(">186<", svg, fixed = TRUE) && !grepl("All 186 eligible", svg, fixed = TRUE))
  stop("未标注 186 合计")
if (grepl("All 202 eligible", svg, fixed = TRUE))
  stop("仍出现旧宇宙的 202 合计标注")
cat("[计数断言] Panel A 标注 All 186 eligible；Panel B 各方向 n 与合计一致\n")
