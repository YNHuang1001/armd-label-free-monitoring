# =============================================================================
# Supplementary Figure S2 - Robustness to model and feature specification
# =============================================================================
# [RENDER LAYER ONLY] 只读 figures/data/*.csv，不重算任何统计量。
#   fig1e_configs_aligned.csv   : 4 配置 x **对齐的 202 task**（同一批 task，可直比）
#   fig1e_configs.csv           : 各配置自身 eligible 计数（202/184/64/66）
#   figS2_config_agreement.csv  : 配置间 task-level Spearman（Python 层算好，含 bootstrap 区间）
#
# 定位（GPT6 明确）：本图的科学问题是 **robustness to analysis specification**，
#   **不是** performance comparison。因此：
#   * 按 protocol 顺序排列（History+LR, History+LGBM, Core+LR, Core+LGBM），**不按 ΔAUROC 高低排序**
#   * Panel B 展示**配置间一致性**而非排名
#   * 不在图上做任何"哪个配置更好"的标注
#   * Panel A 使用**对齐后的同一批 202 task**（否则各配置 eligible universe 不同，不可直比）
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2); library(patchwork); library(readr); library(svglite)
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

AL <- read_csv(file.path(DATA_DIR, "fig1e_configs_aligned.csv"), show_col_types = FALSE)
CS <- read_csv(file.path(DATA_DIR, "fig1e_configs.csv"), show_col_types = FALSE)
AG <- read_csv(file.path(DATA_DIR, "figS2_config_agreement.csv"), show_col_types = FALSE)
# 【2026-10-06 修正】原为硬编码 808（= 4 配置 × 202，旧宇宙）。改为从 fig1b 派生。
UB <- read_csv(file.path(DATA_DIR, "fig1b_universe.csv"), show_col_types = FALSE)
stopifnot(nrow(AL) == 4L * UB$eligible_tasks, nrow(CS) == 4, nrow(AG) == 10)

CFG <- c("P1 History+LR", "S2 History+LGBM", "S1 Core+LR", "S3 Core+LGBM")
LAB <- c("P1 History+LR" = "History + LR", "S2 History+LGBM" = "History + LGBM",
         "S1 Core+LR" = "Core + LR", "S3 Core+LGBM" = "Core + LGBM")
AL$cfg <- factor(LAB[AL$config], levels = LAB[CFG])
AG$cfg_a <- factor(LAB[AG$cfg_a], levels = LAB[CFG])
AG$cfg_b <- factor(LAB[AG$cfg_b], levels = LAB[CFG])

DIR_ORDER <- c("MGB->SU","MGB->UT","SU->MGB","SU->UT","UT->MGB","UT->SU")
AL$dir <- factor(gsub("->", "\u2192", AL$target), levels = gsub("->", "\u2192", DIR_ORDER))
AL$dir <- factor(paste0(AL$source, "\u2192", AL$target), levels = gsub("->", "\u2192", DIR_ORDER))
stopifnot(!any(is.na(AL$dir)))

PT <- "#2B2B2B"; REF <- "#B22222"
C_HIST <- "#5B8DB8"; C_CORE <- "#B8B8B8"

BASE <- theme_classic(base_size = 9, base_family = "Arial") +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(linewidth = 0.20, colour = "grey93"),
        axis.line  = element_line(linewidth = 0.35, colour = "#666666"),
        axis.ticks = element_line(linewidth = 0.35, colour = "#666666"),
        axis.title = element_text(size = 8.5, colour = "#333333"),
        axis.text  = element_text(size = 7.5, colour = "#444444"),
        plot.title = element_text(size = 9.5, face = "bold", colour = "#222222"),
        plot.subtitle = element_text(size = 7.6, colour = "#555555"),
        panel.border = element_blank(),
        strip.background = element_blank(),
        strip.text.y.left = element_text(angle = 0, hjust = 1, size = 8,
                                         face = "bold", colour = "#555555"),
        strip.placement = "outside",
        panel.spacing.y = unit(5, "pt"),
        plot.margin = margin(5, 8, 5, 5))

# =============================================================================
# Panel A - ΔAUROC by direction, faceted by configuration（对齐的同一批 202 task）
# =============================================================================
pA <- ggplot(AL, aes(x = dir, y = d_auroc)) +
  geom_hline(yintercept = 0, linewidth = 0.30, colour = "#8C8C8C") +
  geom_hline(yintercept = -0.05, linetype = "dashed", linewidth = 0.40, colour = REF) +
  geom_boxplot(width = 0.56, outlier.shape = NA, fill = "#E8EDF2",
               colour = "#444444", linewidth = 0.38) +
  geom_point(position = position_jitter(width = 0.14, seed = 1),
             size = 0.6, alpha = 0.45, colour = PT, stroke = 0) +
  facet_grid(cfg ~ ., switch = "y") +
  scale_x_discrete(labels = function(x) gsub("\u2192", "\u2192", x)) +
  scale_y_continuous(limits = c(-0.26, 0.19), breaks = seq(-0.2, 0.1, 0.1)) +
  labs(title = "Degradation by transfer direction under four analysis specifications",
       # \u6298\u6210\u4e09\u884c\u4e14\u6bcf\u884c\u90fd\u77ed\uff1a\u5355\u884c\u8fc7\u957f\u4f1a\u88ab\u753b\u5e03\u53f3\u8fb9\u754c\u622a\u65ad
       # 【2026-10-06 修正】原为硬编码 202（旧宇宙），改为从 fig1b 派生。
       subtitle = paste0("Same aligned ", UB$eligible_tasks,
                         "-task universe for all four configurations.\n",
                         "Dashed line = harmful threshold.\n",
                         "Configurations in protocol order, not ordered by performance."),
       x = NULL, y = expression(Delta*AUROC)) +
  BASE + theme(axis.text.x = element_text(angle = 28, hjust = 1, size = 7.2))

# =============================================================================
# Panel B - 配置间 task-level 一致性（Spearman）
# =============================================================================
M <- AG[AG$cfg_a != AG$cfg_b, ]
M$cfg_a <- factor(M$cfg_a, levels = rev(LAB[CFG]))
M$cfg_b <- factor(M$cfg_b, levels = LAB[CFG])

pB <- ggplot(M, aes(x = cfg_b, y = cfg_a, fill = spearman_rho)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  geom_text(aes(label = sprintf("%.3f", spearman_rho)), size = 2.9, colour = "#1A1A1A") +
  scale_fill_gradient(low = "#F0F4F8", high = "#3C6E9F", limits = c(0.5, 1),
                      name = "Spearman \u03c1", breaks = c(0.6, 0.7, 0.8, 0.9, 1.0)) +
  scale_x_discrete(position = "top") +
  labs(title = "Agreement between configurations",
       subtitle = "Task-level rank agreement of \u0394AUROC; upper triangle only",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 9, base_family = "Arial") +
  theme(axis.text.x.top = element_text(size = 7.6, colour = "#333333"),
        axis.text.y = element_text(size = 7.6, colour = "#333333"),
        panel.grid = element_blank(),
        plot.title = element_text(size = 9.5, face = "bold", colour = "#222222"),
        plot.subtitle = element_text(size = 7.6, colour = "#555555"),
        legend.key.height = unit(24, "pt"), legend.title = element_text(size = 8),
        legend.text = element_text(size = 7.5),
        plot.margin = margin(5, 8, 5, 5))

figS2 <- pA / pB + plot_layout(heights = c(1.85, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(family = "Arial", face = "bold", size = 12, colour = "#111111"))

ggsave(file.path(SUP_DIR, "FigureS2.png"), figS2, width = 6.5, height = 7.8,
       units = "in", dpi = 350, bg = "white")
ggsave(file.path(SUP_DIR, "FigureS2.svg"), figS2, width = 6.5, height = 7.8,
       units = "in", device = svglite::svglite, bg = "white")

cat("WROTE FigureS2.{png,svg}\n")
cat("aligned rows:", nrow(AL), " configs:", nlevels(AL$cfg), " directions:", nlevels(AL$dir), "\n")
cat("agreement range:", sprintf("%.3f - %.3f", min(M$spearman_rho), max(M$spearman_rho)), "\n")
