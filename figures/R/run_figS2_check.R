# =============================================================================
# Supplementary Figure S2 校验
# =============================================================================
.armd_f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
.armd_d <- if (length(.armd_f) && !is.na(.armd_f)) dirname(.armd_f) else "."
ROOT <- Sys.getenv("ARMD_ROOT", unset = normalizePath(file.path(.armd_d, "..", "..")))
setwd(file.path(ROOT, "figures"))
source("R/check_rendered.R", encoding = "UTF-8")

run_and_check(
  script   = "R/FigureS2_config_robustness.R",
  outbase  = "FigureS2",
  outdir   = "supplement",
  logfile  = "log/FigureS2_render_check.txt",
  panels   = c("A", "B"),
  required = c(
    "Degradation by transfer direction under four analysis specifications",
    "Same aligned 186-task universe for all four configurations",
    "Configurations in protocol order, not ordered by performance",
    "Agreement between configurations",
    "Task-level rank agreement",
    "History + LR", "History + LGBM", "Core + LR", "Core + LGBM",
    "MGB\u2192SU", "MGB\u2192UT", "SU\u2192MGB", "SU\u2192UT", "UT\u2192MGB", "UT\u2192SU",
    "0.927", "0.869", "0.687", "0.611", "0.660"
  ),
  forbidden = c(
    # 本图问的是"结论是否依赖配置"，不是"哪个配置更好"
    "better", "Better", "best", "Best", "superior", "outperforms",
    "worse", "Worse", "ranking", "Ranking", "ranked by", "winner",
    "preferred configuration", "recommend"
  )
)

svg <- paste(readLines("supplement/FigureS2.svg", warn = FALSE), collapse = "\n")

# 断言 1：Panel A 的 4 个配置顺序必须是 protocol 顺序（History+LR 在最上）
yof <- function(lab) {
  m <- regmatches(svg, regexpr(sprintf("<text x='[^']*' y='([0-9.]+)'[^>]*>%s", lab), svg))
  if (!length(m)) return(NA_real_)
  as.numeric(sub(".*y='([0-9.]+)'.*", "\\1", m))
}
y1 <- yof("History \\+ LR"); y2 <- yof("History \\+ LGBM")
y3 <- yof("Core \\+ LR");    y4 <- yof("Core \\+ LGBM")
if (any(is.na(c(y1, y2, y3, y4)))) stop("找不到某个配置的 facet 标签")
if (!(y1 < y2 && y2 < y3 && y3 < y4))
  stop("配置顺序不是 protocol 顺序（应 History+LR -> History+LGBM -> Core+LR -> Core+LGBM，自顶向下）")
cat(sprintf("\n[顺序断言] 配置自顶向下 y = %.0f/%.0f/%.0f/%.0f，为 protocol 顺序 —— OK\n", y1, y2, y3, y4))

# 断言 2：同一 feature tier 内的一致性应高于跨 tier（这正是"对分类器稳健、对特征层敏感"的证据）
# 【2026-10-06 修正】原为写死的旧 rho（同 tier 0.927/0.869、跨 tier 0.668/0.713/0.697/0.567）。
#   改为从 figS2_config_agreement.csv 派生：同 tier = 只差 classifier 的两对，
#   跨 tier = 其余四对。断言本身（组内 > 跨组）保持不变。
AG2 <- read.csv("data/figS2_config_agreement.csv", fileEncoding = "UTF-8")
# 去自身配对：**不能**用 rho < 1 —— 自身配对的 rho 是 1.0 减去浮点 epsilon，
# 实测 `1.0000000 < 1` 为 TRUE，过滤不掉。改用配置名是否相同。
AG2 <- AG2[as.character(AG2$cfg_a) != as.character(AG2$cfg_b), ]
# tier 指 History / Core（不是 P/S 前缀）：P1 与 S2 为 History，S1 与 S3 为 Core。
tf <- function(x) ifelse(grepl("History", as.character(x)), "History", "Core")
tier <- tf(AG2$cfg_a); tier_b <- tf(AG2$cfg_b)
same_tier  <- AG2$spearman_rho[tier == tier_b]
cross_tier <- AG2$spearman_rho[tier != tier_b]
stopifnot(length(same_tier) == 2, length(cross_tier) == 4)
if (min(same_tier) <= max(cross_tier)) stop("同 tier 内一致性未高于跨 tier，与结论不符")
cat(sprintf("[一致性断言] 同 tier 内 rho = %s ；跨 tier rho = %s —— 组内一致性更高\n",
            paste(sprintf("%.3f", same_tier), collapse = "/"),
            paste(sprintf("%.3f", cross_tier), collapse = "/")))
