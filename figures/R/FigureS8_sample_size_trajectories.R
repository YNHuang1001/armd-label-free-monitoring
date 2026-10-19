# =============================================================================
# Supplementary Figure S8 - Failure-discrimination across target sample sizes
# =============================================================================
# [RENDER LAYER ONLY] 只读 figures/data/fig5_subsample.csv，不重算任何统计量。
#
# 设计约束（依 GPT6 定稿意见）：
#   * 三层 facet **共享 y 轴**（0.25-0.80），**不用 free_y**。
#     理由：本图要证明"随 target n 变化 AUROC 基本稳定"。自由 y 轴会把 0.01-0.03 的
#     小变化视觉放大，制造出本不存在的波动感，反而破坏结论。
#   * 水平 0.5 虚线：让"有些曲线稳定地低于 0.5"直接可见 → stability is not validity。
#   * 95% 区间用细 error bar，**不用 ribbon**（6 条 ribbon 叠在一起会糊）。
#   * 方法沿用 protocol 顺序，**不按数值排序**；标题中性，不写"performance improves"。
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

D <- read_csv(file.path(DATA_DIR, "fig5_subsample.csv"), show_col_types = FALSE)
D <- D[D$level != "PAIRED_full_minus_100", ]          # 剔除 paired 行
stopifnot(nrow(D) == 75, length(unique(D$method)) == 15)

ORDER <- c("domain_auc", "mmd", "psi_mean", "psi_max", "js_mean", "js_max",
           "pred_mean_shift", "entropy_mean_shift", "entropy_js", "conf_dist_js",
           "pred_pos_prev_shift", "SUDO_AURCC",
           "AUROC_DR", "AUROC_IOW", "AUROC_OM")
LAB <- c(domain_auc = "domain-AUC", mmd = "MMD",
         psi_mean = "PSI, mean", psi_max = "PSI, max",
         js_mean = "JS, mean", js_max = "JS, max",
         pred_mean_shift = "Prediction mean shift",
         entropy_mean_shift = "Entropy mean shift",
         entropy_js = "Entropy JS", conf_dist_js = "Confidence-dist. JS",
         pred_pos_prev_shift = "Predicted prevalence shift",
         SUDO_AURCC = "SUDO-AURCC",
         AUROC_DR = "DR (primary)", AUROC_IOW = "IOW (sensitivity)",
         AUROC_OM = "Outcome model (sensitivity)")
map_role <- function(x) dplyr::case_when(
  x == "L1" ~ "L1  input shift", x == "L2" ~ "L2  prediction shift",
  startsWith(x, "L3") ~ "L3  performance estimation", TRUE ~ NA_character_)
ROLE_LV <- c("L1  input shift", "L2  prediction shift", "L3  performance estimation")

# 每层内 6/6/3 个方法各配固定颜色（层间色系分开，层内明度梯度）。
# **必须用显示名（LAB）命名**：colour 映射的是 lab 而不是 method 代码，
# 用代码命名会得到 "No shared levels found" 并导致所有曲线变成默认灰色。
PAL <- setNames(
  c("#08306B", "#2171B5", "#6BAED6", "#74C476", "#238B45", "#00441B",   # L1
    "#7F2704", "#D94801", "#FDAE6B", "#CB181D", "#FC9272", "#67000D",   # L2
    "#3F007D", "#6A51A3", "#9E9AC8"),                                   # L3
  LAB[ORDER])

D$layer <- factor(map_role(D$role), levels = ROLE_LV)
D$lab   <- factor(LAB[D$method], levels = LAB[ORDER])
D$lev   <- factor(D$level, levels = c("100", "250", "500", "1000", "full"))
D$levx  <- as.numeric(D$lev)
REF <- "#B22222"

BASE <- theme_classic(base_size = 9, base_family = "Arial") +
  theme(panel.grid.major.y = element_line(linewidth = 0.20, colour = "grey93"),
        panel.grid.major.x = element_blank(),
        axis.line  = element_line(linewidth = 0.35, colour = "#666666"),
        axis.ticks = element_line(linewidth = 0.35, colour = "#666666"),
        axis.title = element_text(size = 9, colour = "#333333"),
        axis.text  = element_text(size = 8, colour = "#444444"),
        plot.title = element_text(size = 9.5, face = "bold", colour = "#222222"),
        plot.subtitle = element_text(size = 8, colour = "#555555"),
        panel.border = element_blank(),
        strip.background = element_blank(),
        strip.text.y.left = element_text(angle = 0, hjust = 1, size = 8.5,
                                         face = "bold", colour = "#555555"),
        strip.placement = "outside",
        panel.spacing.y = unit(5, "pt"),
        legend.key.size = unit(9, "pt"), legend.text = element_text(size = 7.2),
        plot.margin = margin(5, 8, 5, 5))

p <- ggplot(D, aes(x = levx, y = auroc_mean, colour = lab, group = lab)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", linewidth = 0.6, colour = REF) +
  geom_errorbar(aes(ymin = auroc_lo, ymax = auroc_hi), width = 0.10,
                linewidth = 0.35, alpha = 0.75) +
  geom_line(linewidth = 0.55) +
  geom_point(size = 1.3) +
  facet_grid(layer ~ ., scales = "fixed", space = "fixed", switch = "y") +
  scale_colour_manual(values = PAL, name = NULL,
                      breaks = LAB[ORDER], labels = LAB[ORDER]) +
  scale_x_continuous(breaks = 1:5,
                     labels = c("100", "250", "500", "1000", "full"),
                     limits = c(0.82, 5.18), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.25, 0.80), breaks = seq(0.3, 0.8, 0.1)) +
  labs(title = "Failure-discrimination across target sample sizes",
       # 副标题必须短：过长会被画布右边界截断
       subtitle = paste0("Mean and 95% interval over 200 resampling repeats.\n",
                         "Shared y-axis across layers; methods in protocol order."),
       x = "Available unlabeled target observations", y = "Harmful-detection AUROC") +
  BASE +
  # 图例必须放底部且限制列数：15 个方法标签很长，右侧 + nrow=3 会展开成 5 列，
  # 总宽超过画布，**把绘图区挤成一条窄缝**（曾表现为"图里完全没有数据"）。
  guides(colour = guide_legend(nrow = 5, byrow = TRUE)) +
  theme(legend.position = "bottom")

ggsave(file.path(SUP_DIR, "FigureS8.png"), p, width = 6.5, height = 7.0,
       units = "in", dpi = 350, bg = "white")
ggsave(file.path(SUP_DIR, "FigureS8.svg"), p, width = 6.5, height = 7.0,
       units = "in", device = svglite::svglite, bg = "white")

cat("WROTE FigureS8.{png,svg}\n")
cat(sprintf("mean range %.4f-%.4f ; CI range %.4f-%.4f\n",
            min(D$auroc_mean), max(D$auroc_mean), min(D$auroc_lo), max(D$auroc_hi)))
cat("methods:", length(unique(D$method)), " rows:", nrow(D), "\n")
