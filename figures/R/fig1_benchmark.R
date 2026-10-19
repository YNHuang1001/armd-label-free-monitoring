# =============================================================================
# Figure 1 - Benchmark construction and failure landscape
# =============================================================================
# [RENDER LAYER ONLY] 只读取 figures/data/*.csv 出图。
# 严禁在本脚本内重新计算任何统计量：不重算 AUROC / 分位数 / bootstrap /
# harmful-stable 判定 / median / P90-P10 / eligible 筛选。
# 所有数值（P10/P50/P90、harmful 计数、各方向 median、各配置 median）均来自
# Python 层已算好的 CSV。
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2); library(patchwork); library(readr)
  library(dplyr); library(scales); library(svglite)
})

ROOT <- Sys.getenv("ARMD_ROOT", unset = normalizePath(file.path(.armd_d, "..", "..")))
DATA <- file.path(ROOT, "figures", "data")
.armd_f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
.armd_d <- if (length(.armd_f) && !is.na(.armd_f)) dirname(.armd_f) else "."
ROOT <- Sys.getenv("ARMD_ROOT", unset = normalizePath(file.path(.armd_d, "..", "..")))
OUT <- file.path(ROOT, "figures", "out")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# ---- 读取（只读，不派生任何统计量）----
UNI <- read_csv(file.path(DATA, "fig1b_universe.csv"),         show_col_types = FALSE)
TAS <- read_csv(file.path(DATA, "fig1c_delta_distribution.csv"), show_col_types = FALSE)
SUM <- read_csv(file.path(DATA, "fig1c_summary.csv"),          show_col_types = FALSE)
DIR <- read_csv(file.path(DATA, "fig1d_by_direction.csv"),     show_col_types = FALSE)
CFA <- read_csv(file.path(DATA, "fig1e_configs_aligned.csv"),  show_col_types = FALSE)
CFS <- read_csv(file.path(DATA, "fig1e_configs.csv"),          show_col_types = FALSE)

# ---- 配色（Okabe-Ito，色盲安全）----
C_MGB <- "#0072B2"; C_SU <- "#D55E00"; C_UT <- "#009E73"
C_HARM <- "#C1272D"; C_STAB <- "#7BA7C7"; C_GREY <- "#4D4D4D"; C_LIGHT <- "#D9D9D9"
TH <- SUM$harmful_threshold       # 从数据包读，不硬编码
SA <- SUM$stable_abs_threshold

BASE <- theme_bw(base_size = 9, base_family = "sans") +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linewidth = 0.20, colour = "grey92"),
        panel.border     = element_rect(linewidth = 0.45, colour = "grey30"),
        axis.ticks       = element_line(linewidth = 0.35),
        plot.title       = element_text(size = 9.5, face = "bold", hjust = 0),
        plot.subtitle    = element_text(size = 7.6, colour = "grey30", hjust = 0),
        legend.key.size  = unit(9, "pt"))

# =============================================================================
# Panel A - three-site source -> target design and the three signal layers
# =============================================================================
sites <- data.frame(lab = c("MGB", "Stanford", "UTSW"),
                    x = c(1.7, 5.0, 8.3), y = 4.30,
                    col = c(C_MGB, C_SU, C_UT), stringsAsFactors = FALSE)

# 注意：geom_curve 的 curvature 只接受标量，传向量会报
# "the condition has length > 1" 并杀死绘图设备 —— 故按符号拆成两个 layer。
edges <- expand.grid(from = 1:3, to = 1:3) |> filter(from != to) |>
  mutate(x = sites$x[from], xend = sites$x[to], y = 3.85, yend = 3.85,
         grp = ifelse(to > from, "fwd", "bwd"))
E_FWD <- filter(edges, grp == "fwd"); E_BWD <- filter(edges, grp == "bwd")

lay <- data.frame(
  y0   = c(0.15, 0.68, 1.21), y1 = c(0.62, 1.15, 1.68),
  fill = c("#F4F4F4", "#E8EEF4", "#F8EEE7"),
  tag  = c("L1", "L2", "L3"),
  txt  = c("input shift  \u2014  domain-AUC \u00b7 MMD \u00b7 PSI \u00b7 JS",
           "prediction shift  \u2014  pred mean \u00b7 entropy \u00b7 entropy-JS \u00b7 conf-JS \u00b7 prevalence \u00b7 SUDO",
           "performance estimation  \u2014  transported AUROC: DR (primary) \u00b7 IOW \u00b7 outcome-model"))

# 双列布局：左 = source（有标签），右 = target（X 可见 / Y sealed），6 条有向箭头
# 比三节点环形布局清晰得多（环形会把 6 条弧线叠成一团）
SRC_X <- 2.35; TGT_X <- 7.65
src_nodes <- data.frame(lab = c("MGB", "Stanford", "UTSW"),
                        y = c(5.30, 4.55, 3.80), col = c(C_MGB, C_SU, C_UT))
tgt_nodes <- data.frame(lab = c("MGB", "Stanford", "UTSW"),
                        y = c(5.30, 4.55, 3.80), col = c(C_MGB, C_SU, C_UT))
arrow_df <- expand.grid(si = 1:3, ti = 1:3) |>
  mutate(x = SRC_X + 1.00, xend = TGT_X - 1.00,
         y = src_nodes$y[si], yend = tgt_nodes$y[ti],
         curve = ifelse(si == ti, 0, ifelse(ti > si, -0.22, 0.22)),
         grp = ifelse(si == ti, "self", ifelse(ti > si, "up", "dn")))
A_UP <- filter(arrow_df, grp == "up"); A_DN <- filter(arrow_df, grp == "dn")
A_SELF <- filter(arrow_df, grp == "self")

pA <- ggplot() +
  # 两列节点
  geom_rect(data = src_nodes, aes(xmin = SRC_X - 1.00, xmax = SRC_X + 1.00,
                                  ymin = y - 0.26, ymax = y + 0.26),
            fill = "white", colour = "grey25", linewidth = 0.45) +
  geom_rect(data = src_nodes, aes(xmin = SRC_X - 1.00, xmax = SRC_X - 0.85,
                                  ymin = y - 0.26, ymax = y + 0.26),
            fill = src_nodes$col, colour = NA) +
  geom_text(data = src_nodes, aes(x = SRC_X + 0.07, y = y, label = lab),
            size = 2.9, fontface = "bold") +
  geom_rect(data = tgt_nodes, aes(xmin = TGT_X - 1.00, xmax = TGT_X + 1.00,
                                  ymin = y - 0.26, ymax = y + 0.26),
            fill = "grey97", colour = "grey25", linewidth = 0.45) +
  geom_rect(data = tgt_nodes, aes(xmin = TGT_X - 1.00, xmax = TGT_X - 0.85,
                                  ymin = y - 0.26, ymax = y + 0.26),
            fill = tgt_nodes$col, colour = NA) +
  geom_text(data = tgt_nodes, aes(x = TGT_X + 0.07, y = y, label = lab),
            size = 2.9, fontface = "bold") +
  # 6 条有向迁移（curvature 为标量，分三层绘制）
  geom_curve(data = A_UP, aes(x = x, y = y, xend = xend, yend = yend), curvature = -0.22,
             arrow = arrow(length = unit(3.8, "pt"), type = "closed"),
             linewidth = 0.40, colour = "grey45") +
  geom_curve(data = A_DN, aes(x = x, y = y, xend = xend, yend = yend), curvature = 0.22,
             arrow = arrow(length = unit(3.8, "pt"), type = "closed"),
             linewidth = 0.40, colour = "grey45") +
  geom_segment(data = A_SELF, aes(x = x, y = y, xend = xend, yend = yend),
               arrow = arrow(length = unit(3.8, "pt"), type = "closed"),
               linewidth = 0.40, colour = "grey45") +
  annotate("text", x = SRC_X, y = 6.42, label = "SOURCE", size = 3.0, fontface = "bold") +
  annotate("text", x = SRC_X, y = 6.14, label = "trained on labeled data", size = 2.5, colour = "grey35") +
  annotate("text", x = TGT_X, y = 6.42, label = "TARGET", size = 3.0, fontface = "bold") +
  annotate("text", x = TGT_X, y = 6.14, label = "X visible  |  Y sealed", size = 2.5, colour = "grey35") +
  annotate("text", x = 5, y = 3.05,
           label = "6 directed transfers per organism-drug  ×  38 organism-drug groups",
           size = 2.7, colour = "grey30") +
  # 三层 signal
  geom_rect(data = lay, aes(xmin = 0.15, xmax = 9.85, ymin = y0, ymax = y1),
            fill = lay$fill, colour = "grey62", linewidth = 0.35) +
  geom_text(data = lay, aes(x = 0.42, y = (y0 + y1) / 2, label = tag),
            hjust = 0, size = 3.0, fontface = "bold") +
  geom_text(data = lay, aes(x = 1.30, y = (y0 + y1) / 2, label = txt),
            hjust = 0, size = 2.55, colour = "grey15") +
  coord_cartesian(xlim = c(0, 10), ylim = c(0, 6.7), clip = "off") +
  theme_void(base_size = 9, base_family = "sans") +
  theme(plot.margin = margin(2, 6, 2, 2))

# =============================================================================
# Panel B - 228 candidate -> 202 eligible
# =============================================================================
bdf <- data.frame(
  step = factor(c("Candidate\ntransfer tasks", "Eligible\n(source AUROC \u2265 0.65)"),
                levels = c("Candidate\ntransfer tasks", "Eligible\n(source AUROC \u2265 0.65)")),
  n = c(UNI$candidate_tasks, UNI$eligible_tasks), fill = c("a", "b"))

pB <- ggplot(bdf, aes(x = step, y = n)) +
  geom_col(aes(fill = fill), width = 0.55, colour = "grey20", linewidth = 0.4) +
  geom_text(aes(label = n), vjust = -0.45, size = 3.5, fontface = "bold") +
  scale_fill_manual(values = c(a = C_LIGHT, b = C_MGB), guide = "none") +
  scale_y_continuous(limits = c(0, max(bdf$n) * 1.24), expand = c(0, 0)) +
  labs(title = "Task universe", x = NULL, y = "number of transfer tasks",
       subtitle = sprintf("%d organisms  \u00b7  %d drugs  \u00b7  %d organism-drug groups",
                          UNI$n_organism, UNI$n_drug, UNI$n_organism_drug_groups)) +
  BASE

# =============================================================================
# Panel C - corrected dAUROC distribution
# =============================================================================
ylim_hi <- max(ggplot_build(ggplot(TAS, aes(x = d_auroc)) +
                              geom_histogram(binwidth = 0.015, boundary = 0))$data[[1]]$count) * 1.34

pC <- ggplot(TAS, aes(x = d_auroc)) +
  geom_histogram(aes(fill = d_auroc <= TH), binwidth = 0.015,
                 colour = "white", linewidth = 0.25, boundary = 0) +
  geom_vline(xintercept = c(SUM$p10, SUM$p90), linetype = "dashed",
             linewidth = 0.42, colour = C_GREY) +
  geom_vline(xintercept = SUM$p50, linewidth = 0.5, colour = "black") +
  geom_vline(xintercept = TH, linetype = "22", linewidth = 0.55, colour = C_HARM) +
  # P10 / P90 标在顶部并向外对齐；P50 标在底部，彻底避开碰撞
  # \u4e09\u4e2a\u5206\u4f4d\u6807\u7b7e\u540c\u9ad8\u3001\u5404\u81ea\u5411\u5916\u5c55\u5f00\uff1aP10 \u5411\u5de6\uff0cP50/P90 \u5411\u53f3\uff0c\u6c34\u5e73\u4e0d\u91cd\u53e0
  annotate("text", x = SUM$p10, y = ylim_hi * 0.95, hjust = 1.08, size = 2.5, colour = C_GREY,
           label = sprintf("P10 %+.3f", SUM$p10)) +
  annotate("text", x = SUM$p50, y = ylim_hi * 0.83, hjust = 1.10, size = 2.6, fontface = "bold",
           label = sprintf("P50 %+.3f", SUM$p50)) +
  annotate("text", x = SUM$p90, y = ylim_hi * 0.95, hjust = -0.08, size = 2.5, colour = C_GREY,
           label = sprintf("P90 %+.3f", SUM$p90)) +
  # harmful \u6807\u6ce8\u653e\u5728\u5de6\u4fa7\u4f4e\u67f1\u533a\u4e0a\u65b9\uff0c\u907f\u5f00\u67f1\u4f53
  annotate("text", x = -0.155, y = ylim_hi * 0.46, hjust = 0.5, size = 2.5,
           colour = C_HARM, lineheight = 0.92,
           label = sprintf("harmful\n\u0394 \u2264 %.2f", TH)) +
  scale_fill_manual(values = c(`TRUE` = C_HARM, `FALSE` = C_STAB), guide = "none") +
  scale_y_continuous(limits = c(0, ylim_hi), expand = c(0, 0)) +
  labs(title = "Cross-system degradation",
       x = expression(Delta*AUROC~"(target \u2212 source internal)"), y = "transfer tasks",
       subtitle = sprintf("harmful = %d   \u00b7   stable (|\u0394| \u2264 %.2f) = %d\nP90 \u2212 P10 = %.3f",
                          SUM$n_harmful, SA, SUM$n_stable, SUM$p90_minus_p10)) +
  BASE

# =============================================================================
# Panel D - six transfer directions (distribution, not just median)
# =============================================================================
ord <- DIR$direction[order(DIR$`50%`)]
TAS$direction <- factor(TAS$direction, levels = ord)

pD <- ggplot(TAS, aes(x = direction, y = d_auroc)) +
  geom_hline(yintercept = TH, linetype = "22", linewidth = 0.45, colour = C_HARM) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey55") +
  geom_boxplot(width = 0.52, outlier.size = 0.45, outlier.colour = "grey50",
               fill = "grey95", colour = "grey30", linewidth = 0.38) +
  geom_point(aes(colour = source_site), position = position_jitter(width = 0.15, seed = 1),
             size = 0.8, alpha = 0.70, stroke = 0) +
  scale_colour_manual(values = c(MGB = C_MGB, SU = C_SU, UT = C_UT), name = "source") +
  labs(title = "Direction of transfer", x = NULL, y = expression(Delta*AUROC),
       subtitle = "ordered by median \u00b7 no mechanistic interpretation") +
  BASE +
  theme(plot.subtitle = element_text(face = "italic"),
        legend.position = "right", legend.title = element_text(size = 7.4),
        legend.text = element_text(size = 7),
        axis.text.x = element_text(angle = 30, hjust = 1, size = 7.6))

# =============================================================================
# Panel E - model/feature class (SAME P1-defined 202-task universe)
# =============================================================================
CFA$config <- factor(CFA$config, levels = CFS$config)
# 轴标签只放短名；median / harmful 作为面板内注释，避免轴文字互撞
# 轴标签拆成两行（"History" / "+ LR"），否则 4 个长标签会互相压住
cfg_lab <- sub("^([A-Za-z]+)\\+(.+)$", "\\1\n+ \\2", sub("^[^ ]+ ", "", CFS$config))

pE <- ggplot(CFA, aes(x = config, y = d_auroc, fill = tier)) +
  geom_hline(yintercept = TH, linetype = "22", linewidth = 0.45, colour = C_HARM) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey55") +
  geom_violin(width = 0.84, linewidth = 0.35, colour = "grey35",
              scale = "width", alpha = 0.55) +
  geom_boxplot(width = 0.14, outlier.shape = NA, fill = "white",
               colour = "grey25", linewidth = 0.38) +
  geom_text(data = data.frame(config = factor(CFS$config, levels = CFS$config)),
            aes(x = config, y = -0.315, label = sprintf("med %+.3f\nharm %d", CFS$median_aligned, CFS$n_harmful_aligned)),
            inherit.aes = FALSE, size = 2.35, lineheight = 0.92, colour = "grey20") +
  scale_fill_manual(values = c(hist = C_MGB, core = "#AEB8C2"), name = "features") +
  scale_x_discrete(labels = cfg_lab) +
  # 范围必须完整容纳 aligned 集（min -0.2421 / max 0.1706），否则 ggplot 会静默丢点
  scale_y_continuous(limits = c(-0.36, 0.20)) +
  labs(title = "Failure is not classifier-specific", x = NULL, y = expression(Delta*AUROC),
       subtitle = "Same P1-defined 202-task benchmark universe") +
  BASE +
  theme(plot.subtitle = element_text(face = "bold"),
        axis.text.x = element_text(size = 7.2),
        legend.position = "right", legend.title = element_text(size = 7.4),
        legend.text = element_text(size = 7))

# =============================================================================
# 拼版：A 占满整行，B|C 与 D|E 各半
# =============================================================================
fig1 <- pA / (pB | pC) / (pD | pE) +
  plot_layout(heights = c(0.88, 1, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 12, face = "bold"))

W <- 7.2; H <- 8.8
ggsave(file.path(OUT, "Figure1_draft_v1.png"), fig1, width = W, height = H, dpi = 400)
ggsave(file.path(OUT, "Figure1_draft_v1.svg"), fig1, width = W, height = H,
       device = svglite::svglite)
cat("WROTE Figure1_draft_v1.{png,svg}\n")
cat("harmful:", SUM$n_harmful, " stable:", SUM$n_stable,
    " P90-P10:", round(SUM$p90_minus_p10, 3), "\n")
