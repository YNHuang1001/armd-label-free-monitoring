# =============================================================================
# Supplementary Figure S7 - Overlap / weight diagnostics for the DR estimator
# =============================================================================
# [RENDER LAYER ONLY] 只读 figures/data/*.csv，不重算任何统计量。
#   fig4_l3_scatter.csv            : 202 task 的 w_ESS / w_max / DR 估计与观测
#   figS7_weight_association.csv   : 两个诊断量与 |DR error| 的关联（Python 层算好）
#
# 目的：主文 Figure 4E 已选用 ESS 而非 w_max 作横轴。本图**为该选择提供可核查依据**，
#       不是重新挑指标，也不主张任何因果/机制解释。
#
# 规格：宽 6.5 in ｜ PNG 350 dpi + SVG ｜ 白底 ｜ Arial ｜ panel label 粗体
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2); library(patchwork); library(readr)
  library(dplyr); library(scales); library(svglite)
})

if (.Platform$OS.type == "windows") {
  try(grDevices::windowsFonts(Arial = grDevices::windowsFont("Arial")), silent = TRUE)
}

.armd_f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
.armd_d <- if (length(.armd_f) && !is.na(.armd_f)) dirname(.armd_f) else "."
ROOT <- Sys.getenv("ARMD_ROOT", unset = normalizePath(file.path(.armd_d, "..", "..")))
DATA_DIR <- file.path(ROOT, "figures", "data")
SUP_DIR  <- file.path(ROOT, "figures", "supplement")
dir.create(SUP_DIR, recursive = TRUE, showWarnings = FALSE)

D  <- read_csv(file.path(DATA_DIR, "fig4_l3_scatter.csv"), show_col_types = FALSE)
AS <- read_csv(file.path(DATA_DIR, "figS7_weight_association.csv"), show_col_types = FALSE)
# 【2026-10-06 修正】原为硬编码 202（旧宇宙）。改为从 fig1b 派生。
#   注：修正后 D$oor（AUROC_DR > 1）恒为空，下游高亮层对 0 例安全。
UB <- read_csv(file.path(DATA_DIR, "fig1b_universe.csv"), show_col_types = FALSE)
stopifnot(nrow(D) == UB$eligible_tasks, nrow(AS) == 2,
          AS$n_tasks[1] == UB$eligible_tasks)
D$abs_err <- abs(D$AUROC_DR - D$auroc_target)
D$oor <- D$AUROC_DR > 1

PT <- "#2B2B2B"; REF <- "#B22222"; OOR <- "#C1272D"

BASE <- theme_classic(base_size = 9, base_family = "Arial") +
  theme(panel.grid.major = element_line(linewidth = 0.20, colour = "grey93"),
        axis.line  = element_line(linewidth = 0.35, colour = "#666666"),
        axis.ticks = element_line(linewidth = 0.35, colour = "#666666"),
        axis.title = element_text(size = 9, colour = "#333333"),
        axis.text  = element_text(size = 8, colour = "#444444"),
        plot.title = element_text(size = 9.5, face = "bold", colour = "#222222"),
        plot.subtitle = element_text(size = 8, colour = "#555555"),
        panel.border = element_blank(),
        plot.margin = margin(5, 8, 5, 6))

mk_hist <- function(col, ttl, xlab, row, brk, lbl) {
  v <- D[[col]]
  rho <- AS$spearman_rho[AS$metric == row]
  ggplot(D, aes(x = v)) +
    geom_histogram(aes(fill = oor), bins = 30, colour = "white", linewidth = 0.20) +
    scale_fill_manual(values = c(`FALSE` = "#8FB0CF", `TRUE` = OOR), guide = "none") +
    scale_x_log10(breaks = brk, labels = lbl) +
    labs(title = ttl,
         subtitle = sprintf("\u03c1 vs |DR error| = %+.3f", rho),
         x = xlab, y = "Transfer tasks") +
    BASE
}
pA <- mk_hist("w_ESS", "Effective sample size", "ESS (log scale)", "ESS",
              c(100, 10000), c("100", "10,000"))
pB <- mk_hist("w_max", "Maximum weight", "Max weight (log scale)", "w_max",
              c(1, 10, 100, 1000), c("1", "10", "100", "1,000"))

# ---- Panel C：两个诊断量的关联对照 ----
AS$lab <- factor(c("Spearman \u03c1 vs |DR error|")[rep(1, 2)], levels = "Spearman \u03c1 vs |DR error|")
AS$metric <- factor(AS$metric, levels = c("ESS", "w_max"))
# 行标签只写中性名称，不在标签内做方向判断。
# 【修正记录】原版在此标了 "(inverse, as expected)" 与 "(opposite direction)"，
# 是把 w_max 的 rho = +0.178 符号读反了：正相关即"权重越大、误差越大"，属**预期方向**。
# 两个诊断量方向一致，差别只在强度（|rho| 0.301 vs 0.178）。
AS$row <- c(ESS = "log10 effective sample size",
            w_max = "log10 maximum weight")[AS$metric]
AS$row <- factor(AS$row, levels = rev(unique(AS$row[order(AS$metric)])))
AS$sign <- ifelse(AS$spearman_rho < 0, "neg", "pos")

pC <- ggplot(AS, aes(y = row)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.6, colour = REF) +
  # 不用 geom_errorbarh()：它在 ggplot2 4.0 已弃用，产生的 deprecation 警告
  # 会命中 check_rendered 的致命模式（"Warning message"）。
  geom_errorbar(aes(xmin = ci_lo, xmax = ci_hi), orientation = "y",
                width = 0, linewidth = 0.6, colour = PT) +
  geom_point(aes(x = spearman_rho, shape = sign), size = 2.2, colour = PT) +
  geom_text(aes(x = spearman_rho, label = sprintf("%+.3f", spearman_rho)),
            vjust = -1.5, size = 2.8, colour = "#333333") +
  scale_shape_manual(values = c(neg = 16, pos = 17), guide = "none") +
  scale_x_continuous(limits = c(-0.50, 0.40), breaks = c(-0.4, -0.2, 0, 0.2, 0.4)) +
  # \u6bcf\u4e2a panel \u4ec5\u7ea6 2 in \u5bbd\uff0c\u6807\u9898/\u526f\u6807\u9898/x \u6807\u7b7e\u5fc5\u987b\u77ed\uff0c\u5426\u5219\u6ea2\u51fa\u5230\u76f8\u90bb panel
  # 折两行：Panel C 的 y 轴标签很长（"log10 effective sample size"），
  # 绘图区仅约 1.6 in，单行标题需 2.4 in 会被右边界截断
  labs(title = "Association with\nabsolute DR error",
       subtitle = "Descriptive, not causal",
       x = "Spearman \u03c1 (95% CI)", y = NULL) +
  BASE +
  theme(legend.position = "none",
        axis.text.y = element_text(size = 7.4, lineheight = 0.95))

figS7 <- (pA | pB | pC) + plot_layout(widths = c(1, 1, 1.15)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(family = "Arial", face = "bold", size = 12, colour = "#111111"))

ggsave(file.path(SUP_DIR, "FigureS7.png"), figS7, width = 6.5, height = 4.0,
       units = "in", dpi = 350, bg = "white")
ggsave(file.path(SUP_DIR, "FigureS7.svg"), figS7, width = 6.5, height = 4.0,
       units = "in", device = svglite::svglite, bg = "white")

cat("WROTE FigureS7.{png,svg}\n")
cat(sprintf("rho(ESS)=%+.3f [%.3f,%.3f]  rho(w_max)=%+.3f [%.3f,%.3f]\n",
            AS$spearman_rho[1], AS$ci_lo[1], AS$ci_hi[1],
            AS$spearman_rho[2], AS$ci_lo[2], AS$ci_hi[2]))
