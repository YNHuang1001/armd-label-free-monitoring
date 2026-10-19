# =============================================================================
# Figure 4 - Can target performance be estimated without target labels?
# =============================================================================
# [RENDER LAYER ONLY] 只读 figures/data/*.csv，不重算任何统计量。
#   fig4_l3_scatter.csv          : 202 task 的 observed / estimated AUROC 与权重诊断
#   fig4d_mae.csv                : 各 estimator 的 MAE / 越界数 / 预测-真值 Spearman
#   fig4d_calibration_{DR,IOW,OM}.csv : 预测 ΔAUROC 固定分档下的实际 ΔAUROC
#   fig4e_ess_association.csv    : Panel E 的关联统计（已预计算，含描述性区间）
#
# 设计约束（依 GPT6 定稿意见）：
#   A/B/C  三个 identity plot **必须同轴**（0.45–1.15），并共用同一条 y=x。
#          不统一会让某个 estimator 显得"更贴线"。
#          【2026-10-06 修正】原注释称该上界是为容纳 DR 的 1.106 而设；CORRECTION 后
#          DR 最大值降为 0.959（observed 0.487–0.934，IOW/OM 上限 0.961），故上界现已
#          宽于数据需要。**几何范围保持不变**（仍完整包含全部取值，改则改变 A/B/C 的相对观感）。
#   D      MAE 条形**按 protocol 顺序 DR -> IOW -> OM**，不按 MAE 高低排序
#          （primary/sensitivity 是设计身份，不是性能名次）。
#   E      全部 N_task 个 task；x = log10(ESS)；DR>1 的 task 高亮（最多标 2 个名字）。
#          【2026-10-06 修正】CORRECTION 后 DR>1 的 task 数为 **0**，故无高亮、无标签；
#          代码对 0 例情形安全（sum(oor)==0 -> tag_tasks 为空）。
#          只写 diagnostic association，不写因果。
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

D  <- read_csv(file.path(DATA_DIR, "fig4_l3_scatter.csv"), show_col_types = FALSE)
MAE<- read_csv(file.path(DATA_DIR, "fig4d_mae.csv"), show_col_types = FALSE)
AS <- read_csv(file.path(DATA_DIR, "fig4e_ess_association.csv"), show_col_types = FALSE)
CAL<- lapply(c(DR = "DR", IOW = "IOW", OM = "OM"), function(k)
  read_csv(file.path(DATA_DIR, paste0("fig4d_calibration_", k, ".csv")), show_col_types = FALSE))

# 【2026-10-06 修正】原为硬编码 202（旧宇宙）。改为从 fig1b_universe.csv 派生。
.uB <- read_csv(file.path(DATA_DIR, "fig1b_universe.csv"), show_col_types = FALSE)
N_ELIG <- as.integer(.uB$eligible_tasks)
stopifnot(nrow(D) == N_ELIG, nrow(MAE) == 3, all(!is.na(D$AUROC_DR)))
D$abs_err_DR <- abs(D$AUROC_DR - D$auroc_target)

# protocol 顺序（**不按 MAE 排序**）
EST   <- c("DR", "IOW", "OM")
ESTLAB<- c(DR = "DR   (primary)", IOW = "IOW   (sensitivity)", OM = "OM   (sensitivity)")
ESTCOL<- c(DR = "#26456B", IOW = "#6E8FAE", OM = "#A9BCCD")
PT <- "#2B2B2B"; REF <- "#B22222"
AXLO <- 0.45; AXHI <- 1.15

BASE <- theme_classic(base_size = 8, base_family = "Arial") +
  theme(panel.grid.major = element_line(linewidth = 0.20, colour = "grey93"),
        axis.line  = element_line(linewidth = 0.35, colour = "#666666"),
        axis.ticks = element_line(linewidth = 0.35, colour = "#666666"),
        axis.title = element_text(size = 7.6, colour = "#333333"),
        axis.text  = element_text(size = 6.8, colour = "#444444"),
        plot.title = element_text(size = 8.6, face = "bold", colour = "#222222"),
        plot.subtitle = element_text(size = 6.6, colour = "#555555"),
        panel.border = element_blank(),
        plot.margin = margin(6, 7, 5, 5))

# =============================================================================
# A / B / C - identity plots（同轴 + 同一 identity line）
# =============================================================================
id_plot <- function(col, lab, sub) {
  dd <- data.frame(obs = D$auroc_target, est = D[[col]])
  corr <- MAE$rho_pred_true[MAE$estimator == lab]
  ggplot(dd, aes(x = obs, y = est)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                linewidth = 0.55, colour = REF) +
    geom_point(size = 0.95, colour = ESTCOL[[lab]], alpha = 0.62, stroke = 0) +
    annotate("text", x = AXLO + 0.03, y = AXHI - 0.05,
             label = sprintf("Spearman \u03c1 = %.3f", corr),
             hjust = 0, vjust = 1, size = 2.4, colour = "#444444") +
    scale_x_continuous(limits = c(AXLO, AXHI), breaks = seq(0.5, 1.1, 0.2)) +
    scale_y_continuous(limits = c(AXLO, AXHI), breaks = seq(0.5, 1.1, 0.2)) +
    coord_fixed(ratio = 1) +
    labs(title = ESTLAB[[lab]], subtitle = sub,
         x = "Observed target AUROC", y = "Estimated target AUROC") +
    BASE
}
pA <- id_plot("AUROC_DR", "DR", "Doubly robust (augmented IPW)")
pB <- id_plot("AUROC_IOW", "IOW", "Inverse-odds weighting")
pC <- id_plot("AUROC_OM", "OM", "Outcome model")

# =============================================================================
# D1 - MAE bars（protocol 顺序）
# =============================================================================
MAE$est <- factor(MAE$estimator, levels = EST)
pD1 <- ggplot(MAE, aes(x = mae, y = est, fill = est)) +
  geom_col(width = 0.58, colour = "#3A3A3A", linewidth = 0.35, show.legend = FALSE) +
  geom_text(aes(label = sprintf("%.4f", mae)), hjust = -0.12, size = 2.5, fontface = "bold") +
  scale_fill_manual(values = ESTCOL) +
  scale_y_discrete(limits = rev(EST), labels = rev(unname(ESTLAB[EST]))) +
  scale_x_continuous(limits = c(0, 0.075), breaks = seq(0, 0.06, 0.02),
                     expand = c(0, 0)) +
  labs(title = "Mean absolute error",
       subtitle = "|estimated \u2212 observed target AUROC|",
       x = "MAE", y = NULL) +
  BASE + theme(axis.text.y = element_text(size = 6.8))

# =============================================================================
# D2 - calibration（同样顺序、固定配色）
# =============================================================================
CALD <- bind_rows(lapply(EST, function(k) {
  cc <- CAL[[k]]; cc$est <- k; cc
}))
CALD$est <- factor(CALD$est, levels = EST)
pD2 <- ggplot(CALD, aes(x = pred_d_mean, y = actual_d_mean,
                        colour = est, shape = est)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              linewidth = 0.55, colour = REF) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 1.9) +
  scale_colour_manual(values = ESTCOL, name = NULL, breaks = EST, labels = ESTLAB[EST]) +
  scale_shape_manual(values = c(DR = 16, IOW = 17, OM = 15), name = NULL,
                     breaks = EST, labels = ESTLAB[EST]) +
  # breaks 显式给标签：seq() 的浮点误差会让 ggplot 退回科学计数法（-1.5e-01）。
  # 该子图绘图区窄，刻度必须减到 3 个，否则标签互相贴住（0.05 间隔放不下 6 个标签）。
  scale_x_continuous(limits = c(-0.17, 0.14),
                     breaks = c(-0.15, -0.05, 0.05),
                     labels = c("-0.15", "-0.05", "0.05")) +
  scale_y_continuous(limits = c(-0.13, 0.05),
                     breaks = c(-0.10, 0.00),
                     labels = c("-0.10", "0.00")) +
  labs(title = "Calibration", subtitle = "predicted \u0394AUROC bin vs observed mean",
       x = "Predicted \u0394AUROC", y = "Observed \u0394AUROC") +
  BASE +
  # 图例放右下空白区：放左上会被红色 identity line 穿过
  theme(legend.position = c(0.98, 0.02), legend.justification = c(1, 0),
        legend.background = element_rect(fill = alpha("white", 0.80), colour = NA),
        legend.key.size = unit(8, "pt"), legend.text = element_text(size = 6.2))

# =============================================================================
# E - DR 误差 vs effective sample size
# =============================================================================
D$oor <- D$AUROC_DR > 1
LAB_MAX <- 2
D$tag <- NA_character_
ord <- order(-D$abs_err_DR[D$oor])
tag_tasks <- D$task_id[D$oor][ord][seq_len(min(LAB_MAX, sum(D$oor)))]
D$tag[D$task_id %in% tag_tasks] <- D$task_id[D$task_id %in% tag_tasks]

pE <- ggplot(D, aes(x = w_ESS, y = abs_err_DR)) +
  geom_point(data = D[!D$oor, ], size = 0.95, colour = PT, alpha = 0.45, stroke = 0) +
  geom_point(data = D[D$oor, ], size = 1.9, colour = REF, stroke = 0) +
  geom_text(data = D[!is.na(D$tag), ], aes(label = tag), hjust = -0.10, vjust = 0.6,
            size = 2.15, colour = REF) +
  scale_x_log10(breaks = c(100, 1000, 10000),
                labels = c("100", "1,000", "10,000")) +
  scale_y_continuous(limits = c(0, 0.45), breaks = seq(0, 0.4, 0.1)) +
  # \u6807\u9898/\u526f\u6807\u9898\u987b\u6536\u5728 panel \u5185\uff08\u66fe\u56e0\u8fc7\u957f\u88ab\u53f3\u8fb9\u7f18\u622a\u65ad\uff09
  # \u30102026-10-06 \u4fee\u6b63\u3011\u539f\u5355\u884c\u6807\u9898\u53f3\u7aef 560.39 > \u753b\u5e03 518.40\uff08\u8d8a\u754c 41.99pt\uff09\uff0c\u6539\u4e3a\u4e24\u884c\u3002
  labs(title = "Estimation error vs\neffective sample size",
       subtitle = sprintf("Spearman \u03c1 = %.3f [%.3f, %.3f]",
                          AS$spearman_rho, AS$lo_95, AS$hi_95),
       x = "Effective sample size (log10 scale)",
       y = "Absolute DR estimation error") +
  BASE

# =============================================================================
# 拼版：上行 3 个 identity plot，下行 D(两子图) + E
# =============================================================================
# 扁平 3x2：嵌套 patchwork 放进 `|` 组合会触发 wrap_dims 报错，故不嵌套。
# D 占两个子图但只标一次："D" 给 MAE bar，" " 给 calibration。
pA  <- pA  + labs(tag = "A")
pB  <- pB  + labs(tag = "B")
pC  <- pC  + labs(tag = "C")
pD1 <- pD1 + labs(tag = "D")
pD2 <- pD2 + labs(tag = " ")
pE  <- pE  + labs(tag = "E")

fig4 <- (pA | pB | pC) / (pD1 | pD2 | pE) +
  plot_layout(heights = c(1, 0.95)) &
  theme(plot.tag = element_text(family = "Arial", face = "bold", size = 12, colour = "#111111"))

ggsave(file.path(OUT_DIR, "Figure4.png"), fig4, width = 7.2, height = 5.4,
       units = "in", dpi = 400, bg = "white")
ggsave(file.path(OUT_DIR, "Figure4.svg"), fig4, width = 7.2, height = 5.4,
       units = "in", device = svglite::svglite, bg = "white")

cat("WROTE Figure4.{png,svg}\n")
cat("MAE:", paste(sprintf("%s=%.4f", MAE$estimator, MAE$mae), collapse = "  "), "\n")
cat("out-of-range DR>1:", sum(D$oor), " | labeled:", paste(tag_tasks, collapse = ", "), "\n")
