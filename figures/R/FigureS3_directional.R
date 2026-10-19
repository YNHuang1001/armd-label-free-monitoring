# =============================================================================
# Supplementary Figure S3 - Directional distribution of degradation
# =============================================================================
# [RENDER LAYER ONLY] 只读 figures/data/*.csv，不重算任何统计量。
#   fig1c_delta_distribution.csv  : 202 task 的 ΔAUROC 与 direction
#   figS3_direction_summary.csv   : 每方向的 eligible / harmful / stable 计数（Python 层算好）
#
# 定位（GPT6 明确）：回答"Fig 1D 展示方向异质性，那六个方向的完整分布是什么样"。
#   * 方向顺序**固定为先 source 后 target**（与 S1 的列一致），**不按结果排序**
#   * **不做 direction ranking、不标 "worst direction"、不解释机制**
#     （正文 Fig 1D 已明确 descriptive only; no mechanistic interpretation，此处继续保持）
#   * Panel B 的增量在于：degradation 是否集中在某些 transfer direction
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

D <- read_csv(file.path(DATA_DIR, "fig1c_delta_distribution.csv"), show_col_types = FALSE)
S <- read_csv(file.path(DATA_DIR, "figS3_direction_summary.csv"), show_col_types = FALSE)
# 【2026-10-06 修正】原为硬编码 202/65/49（旧宇宙）。改为从 D 自身与 fig1b 派生：
#   方向计数表 S 必须与任务级表 D 逐项自洽，且总数等于合格宇宙。
UB <- read_csv(file.path(DATA_DIR, "fig1b_universe.csv"), show_col_types = FALSE)
stopifnot(nrow(D) == UB$eligible_tasks, nrow(S) == 6, sum(S$n_eligible) == nrow(D),
          sum(S$n_harmful) == sum(D$harmful), sum(S$n_stable) == sum(D$stable))

DL <- gsub("->", "\u2192", S$direction[order(S$order)])
D$dir <- factor(gsub("->", "\u2192", D$direction), levels = DL)
S$dir <- factor(gsub("->", "\u2192", S$direction), levels = DL)

PT <- "#2B2B2B"; REF <- "#B22222"

BASE <- theme_classic(base_size = 9, base_family = "Arial") +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(linewidth = 0.20, colour = "grey93"),
        axis.line  = element_line(linewidth = 0.35, colour = "#666666"),
        axis.ticks = element_line(linewidth = 0.35, colour = "#666666"),
        axis.title = element_text(size = 8.5, colour = "#333333"),
        axis.text  = element_text(size = 7.6, colour = "#444444"),
        plot.title = element_text(size = 9.5, face = "bold", colour = "#222222"),
        plot.subtitle = element_text(size = 7.6, colour = "#555555"),
        panel.border = element_blank(),
        plot.margin = margin(5, 8, 5, 5))

# =============================================================================
# Panel A - 六个方向的完整 ΔAUROC 分布
# =============================================================================
pA <- ggplot(D, aes(x = dir, y = d_auroc)) +
  geom_hline(yintercept = 0, linewidth = 0.30, colour = "#8C8C8C") +
  geom_hline(yintercept = -0.05, linetype = "dashed", linewidth = 0.55, colour = REF) +
  geom_violin(width = 0.82, trim = FALSE, fill = "#DCE7EF",
              colour = "#7A93A8", linewidth = 0.35) +
  geom_boxplot(width = 0.16, outlier.shape = NA, fill = "white",
               colour = "#3A3A3A", linewidth = 0.38) +
  geom_point(position = position_jitter(width = 0.10, seed = 1),
             size = 0.75, alpha = 0.45, colour = PT, stroke = 0) +
  geom_text(data = S, aes(x = dir, y = 0.175,
                          label = sprintf("n = %d", n_eligible)),
            inherit.aes = FALSE, size = 2.5, colour = "#555555") +
  scale_y_continuous(breaks = seq(-0.2, 0.1, 0.1)) +
  # 必须用 coord_cartesian：violin 的核密度会超出数据范围，
  # scale limits 会把这些行静默丢弃（曾报 "Removed 23 rows (geom_violin())"）
  coord_cartesian(ylim = c(-0.27, 0.20)) +
  labs(title = "Complete degradation distribution for each directed transfer",
       # 【2026-10-06 修正】原为硬编码 202（旧宇宙），改为从 D 派生。
       subtitle = paste0("All ", nrow(D), " eligible transfers. Dashed line = harmful threshold ",
                         "(\u0394AUROC \u2264 \u22120.05).\nDirection order is fixed by source then ",
                         "target system, not by observed degradation."),
       x = NULL, y = expression(Delta*AUROC)) +
  BASE + theme(axis.text.x = element_text(size = 8.4))

# =============================================================================
# Panel B - 每方向的 eligible / harmful / stable 构成
# =============================================================================
long <- rbind(
  data.frame(dir = S$dir, cls = "Harmful", n = S$n_harmful),
  data.frame(dir = S$dir, cls = "Stable", n = S$n_stable),
  data.frame(dir = S$dir, cls = "Other eligible", n = S$n_eligible - S$n_harmful - S$n_stable))
long$cls <- factor(long$cls, levels = c("Harmful", "Stable", "Other eligible"))

pB <- ggplot(long, aes(x = dir, y = n, fill = cls)) +
  geom_col(width = 0.62, colour = "white", linewidth = 0.35,
           position = position_stack(reverse = TRUE)) +
  geom_text(data = S, aes(x = dir, y = n_eligible, label = n_eligible),
            inherit.aes = FALSE, vjust = -0.6, size = 2.7, colour = "#444444") +
  scale_fill_manual(values = c(Harmful = "#C1272D", Stable = "#7BA7C7",
                               `Other eligible` = "#D9D9D9"), name = NULL) +
  scale_y_continuous(limits = c(0, 44), breaks = seq(0, 40, 10)) +
  labs(title = "Outcome composition of each direction",
       subtitle = "Counts of eligible transfers; harmful = \u0394AUROC \u2264 \u22120.05, stable = |\u0394AUROC| \u2264 0.02",
       x = NULL, y = "Eligible transfers") +
  BASE +
  theme(legend.position = "right", legend.key.size = unit(10, "pt"),
        legend.text = element_text(size = 7.6),
        axis.text.x = element_text(size = 8.4))

figS3 <- pA / pB + plot_layout(heights = c(1.35, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(family = "Arial", face = "bold", size = 12, colour = "#111111"))

ggsave(file.path(SUP_DIR, "FigureS3.png"), figS3, width = 6.5, height = 6.2,
       units = "in", dpi = 350, bg = "white")
ggsave(file.path(SUP_DIR, "FigureS3.svg"), figS3, width = 6.5, height = 6.2,
       units = "in", device = svglite::svglite, bg = "white")

cat("WROTE FigureS3.{png,svg}\n")
cat("directions:", nlevels(D$dir), " eligible:", sum(S$n_eligible),
    " harmful:", sum(S$n_harmful), " stable:", sum(S$n_stable), "\n")
