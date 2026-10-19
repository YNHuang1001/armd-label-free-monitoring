# =============================================================================
# Figure 2 - What can and cannot be inferred about deployment failure without labels
# =============================================================================
# [RENDER LAYER ONLY] 只读 figures/data/*.csv，不重算任何统计量。
#   fig2_forest.csv              : 点估计与 CI（AUROC + Spearman）
#   fig2_bootstrap_metadata.csv  : CI 生成元信息（cluster unit / B / seed / method）
# 不重新 bootstrap、不重算 AUROC / Spearman、不排序、不做 post-hoc inversion。
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

F  <- read_csv(file.path(DATA_DIR, "fig2_forest.csv"), show_col_types = FALSE)
MT <- read_csv(file.path(DATA_DIR, "fig2_bootstrap_metadata.csv"), show_col_types = FALSE)

# 【2026-10-06 修正】原为硬编码 202（旧宇宙）。改为从元数据派生：
#   覆盖任务数必须一致，且等于元数据自报的 coverage。
stopifnot(nrow(F) == 15L,
          length(unique(F$coverage)) == 1L,
          unique(F$coverage) == as.integer(MT$coverage))

# ---- protocol 固定顺序（**不按 AUROC 高低排序**）----
ORDER <- c("domain_auc", "mmd", "psi_mean", "psi_max", "js_mean", "js_max",
           "pred_mean_shift", "entropy_mean_shift", "entropy_js", "conf_dist_js",
           "pred_pos_prev_shift", "SUDO_AURCC",
           "AUROC_DR", "AUROC_IOW", "AUROC_OM")
stopifnot(setequal(ORDER, F$method))
F <- F[match(ORDER, F$method), ]

# ★★ 2026-10-14：右侧 panel 的**取向标签**（逐行）。
#   `spearman` 列在族间混用取向，故必须逐行标明「ρ 的正负相对于什么方向」。
ORIENT <- c(
  native    = "native",      # 原生方向：位移越大 = ρ 越正
  magnitude = "|value|",     # 有符号 L2 取绝对值
  risk      = "risk-score"   # 风险分（−值）：越正 = 越像失效
)
F$ORIENT <- ifelse(F$method %in% c("SUDO_AURCC", "AUROC_DR", "AUROC_IOW", "AUROC_OM"),
                   ORIENT[["risk"]],
            ifelse(F$method %in% c("pred_mean_shift", "entropy_mean_shift",
                                   "pred_pos_prev_shift"),
                   ORIENT[["magnitude"]], ORIENT[["native"]]))

# 分组间留白（L1 / L2 / L3 各占一段）
YPOS <- c(17, 16, 15, 14, 13, 12,     # L1
          10,  9,  8,  7,  6,  5,     # L2
           3,  2,  1)                 # L3
F$ypos <- YPOS

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
F$lab <- factor(LAB[F$method], levels = rev(LAB[ORDER]))

# 分组底纹（不用颜色编码"好坏"）
BANDS <- data.frame(
  ymin = c(11.5, 4.5, 0.4), ymax = c(17.5, 10.5, 3.5),
  grp  = c("L1  input shift", "L2  prediction shift", "L3  performance estimation"),
  fill = c("#F6F6F6", "#EEF3F8", "#FAF2EC"))

theme_f2 <- function() {
  theme_classic(base_size = 8, base_family = "Arial") +
    theme(
      panel.grid.major.x = element_line(linewidth = 0.20, colour = "grey92"),
      panel.grid.major.y = element_blank(),
      axis.line  = element_line(linewidth = 0.35, colour = "#666666"),
      axis.ticks = element_line(linewidth = 0.35, colour = "#666666"),
      axis.title = element_text(size = 8, colour = "#333333"),
      axis.text  = element_text(size = 7, colour = "#444444"),
      plot.title = element_text(size = 9.5, face = "bold", colour = "#222222"),
      plot.subtitle = element_text(size = 7.2, colour = "#555555"),
      panel.border = element_blank(),
      plot.margin = margin(6, 8, 6, 6)
    )
}

PT <- "#2B2B2B"; REF <- "#B22222"

# =============================================================================
# 主 panel：harmful AUROC + cluster-bootstrap 95% CI，0.5 参考线
# =============================================================================
pA <- ggplot(F, aes(y = ypos)) +
  geom_rect(data = BANDS, aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax),
            fill = BANDS$fill, inherit.aes = FALSE) +
  geom_vline(xintercept = 0.5, linetype = "dashed", linewidth = 0.7, colour = REF) +
  geom_errorbar(aes(xmin = auc_lo, xmax = auc_hi), orientation = "y",
                width = 0, linewidth = 0.5, colour = PT) +
  geom_point(aes(x = harmful_auc), size = 1.7, colour = PT) +
  # 上限必须留出分组标题的空间，否则 geom_text 会被静默丢弃（曾发生：L1 标题 y=17.75 > 17.6）
  scale_y_continuous(breaks = F$ypos, labels = F$lab,
                     limits = c(0.4, 18.4), expand = c(0, 0)) +
  scale_x_continuous(limits = c(0.25, 0.78), breaks = seq(0.3, 0.7, 0.1),
                     expand = c(0, 0)) +
  # \u6807\u9898/\u526f\u6807\u9898\u5fc5\u987b\u6536\u5728 A \u9762\u677f\u5bbd\u5ea6\u5185\uff0c\u5426\u5219\u4f1a\u6ea2\u51fa\u538b\u4f4f B \u9762\u677f\uff08\u66fe\u53d1\u751f\uff09
  labs(title = "Failure discrimination without labels",
       subtitle = paste0("harmful AUROC (", MT$coverage, " tasks; ", MT$n_clusters,
                         " organism-drug clusters)",
                         "\n95% CI: organism-drug cluster bootstrap, ",
                         MT$B_requested, " resamples"),
       x = "AUROC for detecting harmful transfer (\u0394AUROC \u2264 \u22120.05)",
       y = NULL) +
  theme_f2() +
  theme(axis.text.y = element_text(size = 7.2, colour = "#333333"))

# 分组标题（置于每段顶部，左对齐）
pA <- pA +
  annotate("text", x = 0.252, y = c(17.90, 10.90, 3.90),
           label = BANDS$grp, hjust = 0, vjust = 0.5,
           size = 2.45, fontface = "bold", colour = "#666666")

# =============================================================================
# 右侧窄 panel：Spearman vs \u0394AUROC，0 参考线，共用方法顺序
# =============================================================================
pB <- ggplot(F, aes(y = ypos)) +
  geom_rect(data = BANDS, aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax),
            fill = BANDS$fill, inherit.aes = FALSE) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.7, colour = REF) +
  geom_errorbar(aes(xmin = sp_lo, xmax = sp_hi), orientation = "y",
                width = 0, linewidth = 0.5, colour = PT) +
  geom_point(aes(x = spearman), size = 1.7, colour = PT) +
  # ★★ 2026-10-14（外部核验要求）：**取向必须在图里逐行标出**，不能只靠图注 ——
  #   因为 `spearman` 这一列在**族与族之间混用取向**：L1 与 entropy/conf-dist 用原生方向、
  #   三个有符号 L2 用 |值|、SUDO 与三个 L3 用风险方向。不标出来，读者会把
  #   `domain_auc +0.13` 与 `AUROC_DR −0.37` 并排读成「两者反向」，而它们其实同号。
  geom_text(aes(x = 0.52, label = ORIENT), hjust = 0, size = 1.7, colour = "#777777") +
  scale_y_continuous(limits = c(0.4, 18.4), expand = c(0, 0)) +
  scale_x_continuous(limits = c(-0.55, 0.98), breaks = c(-0.4, -0.2, 0, 0.2, 0.4),
                     expand = c(0, 0)) +
  # \u30102026-10-06 \u4fee\u6b63\u3011\u539f\u5355\u884c\u6807\u9898\u53f3\u7aef 522.66 > \u753b\u5e03 518.40\uff08\u8d8a\u754c 4.26pt\uff09\u3002\u6539\u4e3a\u4e24\u884c\u3002
  labs(title = "Rank association\nwith \u0394AUROC", subtitle = "Spearman \u03c1 (95% CI)",
       x = "Spearman \u03c1", y = NULL) +
  theme_f2() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        plot.title = element_text(size = 8.2), plot.subtitle = element_text(size = 7),
        plot.margin = margin(6, 6, 6, 2))

# =============================================================================
# 拼版
# =============================================================================
fig2 <- (pA + pB + plot_layout(widths = c(3.05, 1.25))) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(family = "Arial", face = "bold", size = 12, colour = "#111111"))

ggsave(file.path(OUT_DIR, "Figure2.png"), fig2, width = 7.2, height = 4.9,
       units = "in", dpi = 400, bg = "white")
ggsave(file.path(OUT_DIR, "Figure2.svg"), fig2, width = 7.2, height = 4.9,
       units = "in", device = svglite::svglite, bg = "white")

cat("WROTE Figure2.{png,svg}\n")
cat("methods:", nrow(F), " clusters:", MT$n_clusters, " B:", MT$B_requested, "\n")
