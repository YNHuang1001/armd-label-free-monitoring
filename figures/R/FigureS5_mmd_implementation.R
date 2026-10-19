# =============================================================================
# Supplementary Figure S5 - MMD implementation and negative U-statistic
# =============================================================================
# [RENDER LAYER ONLY] 只读 p6b（corrected, blinded）与 figS5_mmd_summary.csv
#
# 目的：**解释最终有效实现**，不是展示"修正前 vs 修正后"。
#   旧 MMD 实现（未按 source-train 标准化、bandwidth 取自 pooled source+target）
#   已判定为 invalid，其原始产物已保全于 分析结果/p6_superseded_pre_correction/
#   （含 sha256），但**不作为科学结果呈现**，除非审稿人要求。
#
# 措辞边界（GPT6 明确）：
#   * 不得写 "systematic negative bias"。负值只是 unbiased U-statistic 估计器在
#     当前 sample / kernel / bandwidth 配置下的有限样本实现结果。
#   * 200/202 为负是**该配置下的实现事实**，不等于估计器有系统偏差。
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
RES_DIR  <- file.path(ROOT, "\u5206\u6790\u7ed3\u679c")
DATA_DIR <- file.path(ROOT, "figures", "data")
SUP_DIR  <- file.path(ROOT, "figures", "supplement")
dir.create(SUP_DIR, recursive = TRUE, showWarnings = FALSE)

B <- read_csv(file.path(RES_DIR, "p6b_label_free_signals_blinded.csv"), show_col_types = FALSE)
S <- read_csv(file.path(DATA_DIR, "figS5_mmd_summary.csv"), show_col_types = FALSE)

# 盲态/洁净断言：本图只涉及 X-space 量，不得含任何 target-Y 派生字段
FORB <- c("auroc_target", "d_auroc", "delta_auroc", "harmful", "stable")
hit <- names(B)[vapply(names(B), function(c)
  any(vapply(FORB, function(x) grepl(x, tolower(c), fixed = TRUE), logical(1))), logical(1))]
if (length(hit)) stop("[断言] p6b 含 target-Y 派生字段: ", paste(hit, collapse = ", "))
# 【2026-10-06 修正】原为硬编码 202（旧宇宙）。改为从 fig1b 派生。
UB <- read_csv(file.path(DATA_DIR, "fig1b_universe.csv"), show_col_types = FALSE)
stopifnot(nrow(B) == UB$eligible_tasks, S$n_tasks == UB$eligible_tasks)

PT <- "#2B2B2B"; REF <- "#B22222"; C_HIST <- "#8FB0CF"

BASE <- theme_classic(base_size = 9, base_family = "Arial") +
  theme(axis.line  = element_line(linewidth = 0.35, colour = "#666666"),
        axis.ticks = element_line(linewidth = 0.35, colour = "#666666"),
        axis.title = element_text(size = 8.5, colour = "#333333"),
        axis.text  = element_text(size = 7.5, colour = "#444444"),
        plot.title = element_text(size = 9.5, face = "bold", colour = "#222222"),
        plot.subtitle = element_text(size = 7.8, colour = "#555555"),
        panel.border = element_blank(),
        plot.margin = margin(5, 8, 5, 5))

# =============================================================================
# Panel A - corrected MMD² across all eligible tasks（负值原样保留）
# =============================================================================
dA <- data.frame(mmd = B$mmd)
pA <- ggplot(dA, aes(x = mmd)) +
  geom_histogram(bins = 30, fill = C_HIST, colour = "white", linewidth = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.6, colour = REF) +
  # \u30102026-10-06 \u4fee\u6b63\u3011CS-06 \u540e\u5206\u5e03\u7531\u300c200/202 \u4e3a\u8d1f\u3001\u4e2d\u4f4d \u22120.359\u300d\u53d8\u4e3a\u300c\u5168\u90e8\u975e\u8d1f\u3001
  #   \u4e2d\u4f4d +0.020\u300d\uff0c\u6700\u9ad8\u67f1\u56e0\u6b64\u7d27\u8d34\u9762\u677f\u9876\u90e8\uff08\u4ec5\u4f59 ~8pt\uff09\u3002\u539f vjust 1.6 / 3.8 \u4f7f
  #   "median \u2026" \u88ab\u6700\u9ad8\u67f1**\u7a7f\u5b57**\u3001"MMD\u00b2 = 0" \u8d34\u9876\u3002\u6539\u4e3a\uff1a\u2460 \u9876\u90e8\u52a0 20% \u7559\u767d\uff1b
  #   \u2461 \u4e24\u6761\u6807\u6ce8\u90fd\u653e\u8fdb\u8be5\u7559\u767d\u5185\u5e76\u4e0a\u4e0b\u9519\u5f00\u3002
  #   \u2462 \u518d\u4fee\u6b63\uff1avjust 0.9 \u65f6\u5b57\u5f62\u76d2\u9876 = 54.19\uff0c\u800c panel \u9876 = 54.97 \u2014\u2014 "MMD\u00b2" \u7684
  #      \u4e0a\u6807**\u88ab\u88c1\u6389 0.78pt**\u3002\u4e24\u6761\u5404\u81ea\u4e0b\u79fb\uff1a0.9 -> 1.3\u30012.6 -> 4.0\uff08\u4e92\u4e0d\u76f8\u78b0\u3001
  #      \u4e5f\u90fd\u4ecd\u5728\u6700\u9ad8\u67f1\u9876 84.36 \u4e4b\u4e0a\uff09\u3002
  annotate("text", x = 0.01, y = Inf, vjust = 1.3, hjust = -0.05,
           size = 2.6, colour = REF,
           label = sprintf("MMD\u00b2 = 0")) +
  annotate("text", x = S$mmd_median, y = Inf, vjust = 4.0, hjust = -0.06,
           size = 2.6, colour = "#333333",
           label = sprintf("median %.3f", S$mmd_median)) +
  # 标题/副标题必须短：Panel A/B 并排时过长会互相压住
  # \u30102026-10-06 \u4fee\u6b63\u3011\u6807\u9898\u539f\u4e3a\u786c\u7f16\u7801 202\uff1b\u526f\u6807\u9898\u539f\u5199 "negatives retained"\uff0c
  #   \u4f46\u4fee\u6b63\u540e n_negative = 0\uff0c\u6545\u6539\u4e3a\u4e2d\u6027\u7684\u5b9e\u73b0\u63cf\u8ff0\uff08\u6570\u503c\u4ecd\u5168\u90e8\u4ece S \u6d3e\u751f\uff09\u3002
  labs(title = sprintf("MMD\u00b2 across %d tasks", S$n_tasks),
       subtitle = sprintf("Unbiased U-statistic estimator.\n%d of %d tasks negative.",
                          S$n_negative, S$n_tasks),
       x = "MMD\u00b2 (native units, source-standardised feature space)",
       y = "Transfer tasks") +
  # 顶部留白：为上面两条标注腾出空间（见其注释）
  scale_y_continuous(expand = expansion(mult = c(0, 0.20))) +
  BASE

# =============================================================================
# Panel B - replicate spread：|MMD²| 与 5 次 draw 的 SD
# =============================================================================
dB <- data.frame(abs_mmd = abs(B$mmd), sd = B$mmd_sd)
pB <- ggplot(dB, aes(x = abs_mmd, y = sd)) +
  # 不画 SD = |MMD2| 对角线：所有点离它极远，该参考线不提供任何信息，只是噪声
  geom_point(size = 0.95, colour = PT, alpha = 0.55, stroke = 0) +
  annotate("text", x = 0.05, y = 0.125, hjust = 0, size = 2.5, colour = "#444444",
           label = sprintf("median SD = %.4f\nmedian |MMD\u00b2| = %.3f", S$sd_median, abs(S$mmd_median))) +
  labs(title = "Replicate spread",
       subtitle = sprintf("SD across %d fixed subsample draws;\nno task had SD exceeding its |MMD\u00b2|.",
                          S$n_draws),
       x = "|MMD\u00b2|", y = "SD across draws") +
  scale_x_continuous(limits = c(0, 0.52)) + scale_y_continuous(limits = c(0, 0.135)) +
  BASE

# =============================================================================
# Panel C - 实现流程（文字化 schematic）
# =============================================================================
flow <- data.frame(
  x = c(1, 2.35, 3.7, 5.05, 6.4),
  lab = c("Source-defined\nstandardisation",
          "RBF kernel\non 15 features",
          "Bandwidth: median\nheuristic on\nsource-train pairs",
          "Unbiased\nU-statistic\nMMD\u00b2",
          "Fixed reference:\nm = 1000, K = 5\ndraws"))

pC <- ggplot(flow) +
  annotate("rect", xmin = flow$x - 0.62, xmax = flow$x + 0.62,
           ymin = 0.25, ymax = 0.95, fill = "#F2F5F8",
           colour = "#8FA3B8", linewidth = 0.4) +
  geom_text(aes(x = x, y = 0.60, label = lab), size = 2.35, lineheight = 0.95,
            colour = "#26456B") +
  annotate("segment", x = flow$x[-5] + 0.62, xend = flow$x[-1] - 0.62,
           y = 0.60, yend = 0.60, linewidth = 0.4, colour = "#8FA3B8",
           arrow = grid::arrow(length = grid::unit(0.055, "inches"), type = "closed")) +
  coord_cartesian(xlim = c(0.2, 7.2), ylim = c(0, 1.35), clip = "off") +
  labs(title = "Implementation of the frozen estimator",
       subtitle = paste0("Frozen definitions; no tuning on target data. ",
                         "Reproducibility notes for the MMD component only.")) +
  theme_void(base_size = 9, base_family = "Arial") +
  theme(plot.title = element_text(size = 9.5, face = "bold", colour = "#222222"),
        plot.subtitle = element_text(size = 7.8, colour = "#555555"),
        plot.title.position = "plot",
        plot.margin = margin(4, 8, 4, 4))

figS5 <- (pA | pB) / pC + plot_layout(heights = c(1.55, 0.72)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(family = "Arial", face = "bold", size = 12, colour = "#111111"))

ggsave(file.path(SUP_DIR, "FigureS5.png"), figS5, width = 6.5, height = 5.4,
       units = "in", dpi = 350, bg = "white")
ggsave(file.path(SUP_DIR, "FigureS5.svg"), figS5, width = 6.5, height = 5.4,
       units = "in", device = svglite::svglite, bg = "white")

# ★★ 2026-10-14：曾在此处加过 TIFF 导出 —— **已撤**。
#   理由：本仓 TIFF 的**唯一正规生产者是根目录的 `_make_tiff.py`**（统一 LZW/RGB/dpi 标签，
#   并做「S 图不得硬写 400」的防伪检查）。在渲染脚本里另开一条 TIFF 通道 = **第二生产者**。
cat(sprintf("negative %d/%d ; median MMD2 %.4f ; median SD %.4f\n",
            S$n_negative, S$n_tasks, S$mmd_median, S$sd_median))
