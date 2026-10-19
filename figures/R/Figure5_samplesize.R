# =============================================================================
# Figure 5 - How much unlabeled target data are needed for stable monitoring?
# =============================================================================
# [RENDER LAYER ONLY] 只读 figures/data/*.csv，不重算任何统计量。
#   fig5_subsample.csv        : 5 档 AUROC + 每方法 paired (full − n=100)
#   fig5b_across_n_range.csv  : 每方法跨 5 档的 AUROC 极差（Python 层已算好）
#
# 设计约束（依 GPT6 定稿意见）：
#   主文**不画**十五条近乎平行的五档折线（spaghetti 图）——完整轨迹放 Figure S8。
#   A  端点变化：paired (full − n=100) 的均值与 95% 区间，0 为参考线。
#      **不特别标注/高亮唯一名义显著的方法**，避免把单个 nominal finding 做成主故事。
#   B  整个 sample-size range 内漂了多少：跨 100/250/500/1000/full 的 AUROC 极差。
#      两 panel 均按 L1 -> L2 -> L3 分层、层内 protocol 顺序，**不排序、不排名**。
#
# caption 必须保留的解释边界：stable ≠ good。
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
OUT_DIR  <- file.path(ROOT, "figures", "out")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

S5 <- read_csv(file.path(DATA_DIR, "fig5_subsample.csv"), show_col_types = FALSE)
RG <- read_csv(file.path(DATA_DIR, "fig5b_across_n_range.csv"), show_col_types = FALSE)

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
         AUROC_DR = "DR   (primary)", AUROC_IOW = "IOW   (sensitivity)",
         AUROC_OM = "Outcome model   (sensitivity)")
map_role <- function(x) dplyr::case_when(
  x == "L1" ~ "L1", x == "L2" ~ "L2", startsWith(x, "L3") ~ "L3", TRUE ~ NA_character_)
ROLE_LAB <- c(L1 = "L1  input shift", L2 = "L2  prediction shift",
              L3 = "L3  performance estimation")

PA <- S5[S5$level == "PAIRED_full_minus_100", ]
stopifnot(nrow(PA) == 15, nrow(RG) == 15)
PA$layer  <- factor(ROLE_LAB[map_role(PA$role)], levels = ROLE_LAB)
RG$layer  <- factor(ROLE_LAB[map_role(RG$role)], levels = ROLE_LAB)
PA$lab <- factor(LAB[PA$method], levels = rev(LAB[ORDER]))
RG$lab <- factor(LAB[RG$method], levels = rev(LAB[ORDER]))
YALL <- rev(LAB[ORDER])

PT <- "#2B2B2B"; REF <- "#B22222"

BASE <- theme_classic(base_size = 8, base_family = "Arial") +
  theme(panel.grid.major.x = element_line(linewidth = 0.20, colour = "grey92"),
        panel.grid.major.y = element_blank(),
        axis.line  = element_line(linewidth = 0.35, colour = "#666666"),
        axis.ticks = element_line(linewidth = 0.35, colour = "#666666"),
        axis.title = element_text(size = 7.6, colour = "#333333"),
        axis.text  = element_text(size = 6.9, colour = "#444444"),
        plot.title = element_text(size = 8.6, face = "bold", colour = "#222222"),
        plot.subtitle = element_text(size = 6.6, colour = "#555555"),
        panel.border = element_blank(),
        strip.background = element_blank(),
        strip.text.y.left = element_text(angle = 0, hjust = 1, size = 6.8,
                                         face = "bold", colour = "#555555"),
        strip.placement = "outside",
        panel.spacing.y = unit(4, "pt"),
        plot.margin = margin(6, 6, 5, 4))

# =============================================================================
# Panel A - paired (full - n=100)
# =============================================================================
# 【2026-10-06 修正·重要】原 limits 硬编码 c(-0.095, 0.085)，是为修正前的数据设的。
#   修正后 AUROC_DR 的配对差为 +0.261 (95% 区间 +0.170 ~ +0.341)，**整条落在此范围之外**，
#   于是 DR 这一个点被静默丢弃（渲染日志 "Removed 1 row ... (`geom_point()`)"）——
#   而它恰是本图要传达的主要发现。改为由数据派生（含区间端点与点估计的极值）。
XA_R <- range(c(PA$auroc_lo, PA$auroc_hi, PA$auroc_mean))

pA <- ggplot(PA, aes(y = lab)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.6, colour = REF) +
  geom_errorbar(aes(xmin = auroc_lo, xmax = auroc_hi), orientation = "y",
                width = 0, linewidth = 0.5, colour = PT) +
  geom_point(aes(x = auroc_mean), size = 1.6, colour = PT) +
  facet_grid(layer ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_y_discrete(drop = TRUE) +
  # 【2026-10-06 修正·重要】原 limits 硬编码 c(-0.095, 0.085)，是为修正前的数据设的。
  #   修正后 AUROC_DR 的配对差为 +0.261 (95% 区间 +0.170 ~ +0.341)，**整条落在此范围之外**，
  #   于是 DR 这一个点被静默丢弃（渲染日志 "Removed 1 row ... (`geom_point()`)"）——
  #   而它恰是本图要传达的主要发现。改为由数据派生（含区间端点的极值）。
  scale_x_continuous(limits = c(XA_R[1] - 0.05 * diff(XA_R), XA_R[2] + 0.05 * diff(XA_R)),
                     breaks = scales::pretty_breaks(6), expand = c(0, 0)) +
  labs(title = "Change from full cohort to n = 100",
       subtitle = paste0("Paired difference within each repeat ",
                         "(AUROC_full \u2212 AUROC_n=100) \u00b7 95% interval"),
       x = expression(AUROC[full] - AUROC[n == 100]), y = NULL) +
  BASE

# =============================================================================
# Panel B - across-N range
# =============================================================================
pB <- ggplot(RG, aes(y = lab)) +
  geom_segment(aes(x = 0, xend = across_n_range, y = lab, yend = lab),
               linewidth = 0.5, colour = "#8FA3B8") +
  geom_point(aes(x = across_n_range), size = 1.7, colour = PT) +
  geom_text(aes(x = across_n_range, label = sprintf("%.3f", across_n_range)),
            hjust = -0.22, size = 2.2, colour = "#555555") +
  facet_grid(layer ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_y_discrete(drop = TRUE) +
  # 【2026-10-06 修正·重要】原 limits 硬编码 c(0, 0.052)。修正后 AUROC_DR 的跨档极差为
  #   0.261，远超 0.052，被静默丢弃（同 Panel A）。改为由数据派生，并留出数值标签的余量。
  scale_x_continuous(limits = c(0, max(RG$across_n_range) * 1.18),
                     breaks = scales::pretty_breaks(6), expand = c(0, 0)) +
  labs(title = "Variation across all five target sizes",
       subtitle = "Range of harmful-detection AUROC over 100 / 250 / 500 / 1000 / full",
       x = "AUROC range (max \u2212 min across target sizes)", y = NULL) +
  # 两个 panel 上下堆叠，各自带方法标签，读者不必跨 panel 数行
  BASE

# =============================================================================
# 拼版
# =============================================================================
# 上下堆叠而非左右并排：左右并排会把 A 的绘图区压到约 1.8 in，
# 而标题宽约 2.2 in，必然溢出到 B。堆叠后两 panel 各占满宽度，标题从容。
# 两 panel 行结构完全一致，堆叠后行对齐，读起来也像一张双行 forest。
fig5 <- pA / pB + plot_layout(heights = c(1, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(family = "Arial", face = "bold", size = 12, colour = "#111111"))

ggsave(file.path(OUT_DIR, "Figure5.png"), fig5, width = 7.2, height = 7.4,
       units = "in", dpi = 400, bg = "white")
ggsave(file.path(OUT_DIR, "Figure5.svg"), fig5, width = 7.2, height = 7.4,
       units = "in", device = svglite::svglite, bg = "white")

sig <- PA[PA$auroc_lo > 0 | PA$auroc_hi < 0, ]
cat("WROTE Figure5.{png,svg}\n")
cat("paired CI 不含 0 的方法:", if (nrow(sig)) paste(sig$method, collapse = ", ") else "无", "\n")
cat(sprintf("across-N range: max=%.4f median=%.4f min=%.4f\n",
            max(RG$across_n_range), median(RG$across_n_range), min(RG$across_n_range)))
