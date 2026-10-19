# =============================================================================
# Figure 1 v2 - Label-free AMR deployment benchmark
# =============================================================================
# [RENDER LAYER ONLY]
#   - 只读取冻结的 figures/data/*.csv
#   - 不重算 AUROC / bootstrap / 分位数 / harmful 判定 / eligible 筛选 / median
#   - 不覆盖 v1（figures/R/fig1_benchmark.R）
#
# v2 相对 v1 的改动（依 GPT6 定稿意见）：
#   A  删除 3 条 same-site 箭头（MGB->MGB / SU->SU / UT->UT）；L1->L2->L3 改为正序
#   E  四组重排为 History+LR / History+LGBM / Core+LR / Core+LGBM；标题收敛；去 legend
#   B-E 去掉四周黑色 panel 边框
#   C  改为单色直方图 + 淡红 harmful 背景 + 红色阈值虚线（不再在绘图层按 ΔAUROC 派生颜色）
#
# 与 GPT6 提供脚本的差异（因列名/文件名不符，否则会在 4 处 stop()）：
#   B  读 fig1b_universe_long.csv（原脚本假设的长表；宽表无 stage/n 列）
#   D  读 fig1c_delta_distribution.csv（原脚本 ^fig1d.* 解析到的是 describe() 汇总表，无 ΔAUROC 列）
#      并用 **Python 层预计算的 direction_order**（原脚本 fallback 到 unique() 首次出现顺序，
#      与副标题 "Ordered by precomputed median" 不符）
#   E  汇总列名是 median_aligned / n_harmful_aligned
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2); library(patchwork); library(readr)
  library(dplyr); library(scales); library(svglite)
})

# Arial 位于 Windows 字体目录 %WINDIR%/Fonts/arial.ttf，但未注册进 R 的 Windows 字体库；
# 不注册会导致 PNG 设备静默 fallback（"font family not found"），排版度量与所见不符。
if (.Platform$OS.type == "windows") {
  try(grDevices::windowsFonts(Arial = grDevices::windowsFont("Arial")), silent = TRUE)
}

.armd_f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
.armd_d <- if (length(.armd_f) && !is.na(.armd_f)) dirname(.armd_f) else "."
ROOT <- Sys.getenv("ARMD_ROOT", unset = normalizePath(file.path(.armd_d, "..", "..")))
DATA_DIR <- file.path(ROOT, "figures", "data")
OUT_DIR  <- file.path(ROOT, "figures", "out")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

rd <- function(f) read_csv(file.path(DATA_DIR, f), show_col_types = FALSE)

dB  <- rd("fig1b_universe_long.csv")     # stage_key / n / excluded_note
uB  <- rd("fig1b_universe.csv")           # 宇宙计数（唯一权威来源）
dC  <- rd("fig1c_delta_distribution.csv") # 任务级：d_auroc / direction / direction_order / source_site
sC  <- rd("fig1c_summary.csv")            # P10/P50/P90 / harmful 阈值 / 计数（全部预计算）
dE  <- rd("fig1e_configs_aligned.csv")    # 4 x N_task 已对齐
sE  <- rd("fig1e_configs.csv")            # 各配置 median_aligned / n_harmful_aligned

# 【2026-10-06 修正】原为硬编码 202 / 38 / 808（旧宇宙）。改为从 fig1b_universe.csv 派生，
#   以免下次修正后再次静默画错。
N_ELIG   <- as.integer(uB$eligible_tasks)
N_GROUPS <- as.integer(uB$n_organism_drug_groups)
N_ORG    <- as.integer(uB$n_organism)
N_DRUG   <- as.integer(uB$n_drug)
stopifnot(N_ELIG == nrow(dC), N_GROUPS == dplyr::n_distinct(paste(dC$organism, dC$drug)))

# ---- 冻结基准值（全部来自数据包，不硬编码）----
TH       <- sC$harmful_threshold
SA       <- sC$stable_abs_threshold
P10      <- sC$p10
P50      <- sC$p50
P90      <- sC$p90
N_HARM   <- sC$n_harmful
N_STABLE <- sC$n_stable
SPREAD   <- sC$p90_minus_p10

# dE 是「配置 × 任务」的对齐长表（4 个配置），故行数必为任务数的 4 倍。
stopifnot(nrow(dE) == 4L * nrow(dC), nrow(dC) == dplyr::n_distinct(dC$task_id))

fmt_s <- function(x, d = 3) sprintf(paste0("%+.", d, "f"), x)

# ---- 配色 ----
C_MGB <- "#0072B2"; C_SU <- "#D55E00"; C_UT <- "#009E73"
C_RED <- "#C1272D"; C_GREY <- "#737373"; C_LIGHTGREY <- "#D9D9D9"
C_HIST <- "#91B5D0"; C_CORE <- "#B8B8B8"

theme_pub <- function() {
  theme_classic(base_size = 8, base_family = "Arial") +
    theme(
      plot.title    = element_text(size = 9.5, face = "bold", colour = "#222222",
                                   margin = margin(b = 3)),
      plot.subtitle = element_text(size = 7.2, colour = "#555555",
                                   margin = margin(b = 5)),
      axis.title    = element_text(size = 8, colour = "#333333"),
      axis.text     = element_text(size = 7, colour = "#444444"),
      axis.line     = element_line(linewidth = 0.35, colour = "#666666"),
      axis.ticks    = element_line(linewidth = 0.35, colour = "#666666"),
      panel.border  = element_blank(),
      legend.title  = element_text(size = 7),
      legend.text   = element_text(size = 7),
      plot.margin   = margin(6, 7, 6, 7)
    )
}

# =============================================================================
# Panel A - six cross-system directed transfers + L1 -> L2 -> L3
# =============================================================================
pA <- ggplot() + coord_cartesian(xlim = c(0, 10), ylim = c(0, 10), clip = "off") +
  theme_void(base_family = "Arial") + theme(plot.margin = margin(4, 10, 2, 10))

pA <- pA +
  # 【2026-10-06 修正】原 SOURCE 与下一行基线间距仅 4.59pt（8.54px 字），两行几乎相贴。
  #   本 panel 1 数据单位 ≈ 9.77pt，故把间距放到 ~0.97 单位 ≈ 9.5pt（≈1.12 × majorant 字号）。
  annotate("text", x = 2.75, y = 9.90, label = "SOURCE", size = 3.0, fontface = "bold") +
  annotate("text", x = 2.75, y = 8.93, label = "trained on labeled data",
           size = 2.35, colour = "#595959") +
  annotate("text", x = 7.25, y = 9.90, label = "TARGET", size = 3.0, fontface = "bold") +
  annotate("text", x = 7.25, y = 8.93, label = "X visible  |  Y sealed",
           size = 2.35, colour = "#595959")

nodes <- data.frame(site = c("MGB", "Stanford", "UTSW"), y = c(8.00, 6.80, 5.60),
                    colour = c(C_MGB, C_SU, C_UT), stringsAsFactors = FALSE)

for (i in seq_len(nrow(nodes))) {
  pA <- pA +
    annotate("rect", xmin = 1.75, xmax = 3.75, ymin = nodes$y[i] - 0.36,
             ymax = nodes$y[i] + 0.36, fill = "white", colour = "#404040", linewidth = 0.45) +
    annotate("rect", xmin = 1.75, xmax = 1.90, ymin = nodes$y[i] - 0.36,
             ymax = nodes$y[i] + 0.36, fill = nodes$colour[i], colour = NA) +
    annotate("text", x = 2.80, y = nodes$y[i], label = nodes$site[i],
             size = 2.7, fontface = "bold") +
    annotate("rect", xmin = 6.25, xmax = 8.25, ymin = nodes$y[i] - 0.36,
             ymax = nodes$y[i] + 0.36, fill = "#F7F7F7", colour = "#404040", linewidth = 0.45) +
    annotate("rect", xmin = 6.25, xmax = 6.40, ymin = nodes$y[i] - 0.36,
             ymax = nodes$y[i] + 0.36, fill = nodes$colour[i], colour = NA) +
    annotate("text", x = 7.30, y = nodes$y[i], label = nodes$site[i],
             size = 2.7, fontface = "bold")
}

# ---- 恰好 6 条跨院有向迁移；**无 same-site 箭头** ----
# \u4e3a\u4ec0\u4e48\u4e0d\u518d\u9010\u6761\u753b 6 \u6761\u5f27\u7ebf\uff1a
#   \u5728\u5faa\u73af\u91cc\u5199 aes(y = arcs$y[i]) \u4f1a\u8e29 lazy-eval \u9677\u9631 \u2014\u2014 aes() \u5b58\u7684\u662f**\u8868\u8fbe\u5f0f**\u800c\u975e\u503c\uff0c
#   \u5230 ggsave() build \u65f6\u624d\u6c42\u503c\uff0c\u6b64\u65f6 i \u5df2\u662f\u5faa\u73af\u672b\u503c\uff0c\u516d\u6761\u7ebf\u4f1a\u5168\u90e8\u53e0\u6210\u540c\u4e00\u6761
#   \uff08\u8868\u73b0\u4e3a\u5168\u90e8\u4ece UTSW \u51fa\u53d1\u3001\u5168\u90e8\u6c47\u5165 Stanford\uff09\u3002\u672c panel \u53ea\u9700\u8bf4\u660e benchmark \u8bbe\u8ba1\uff0c
#   6 \u4e2a\u5177\u4f53\u65b9\u5411\u5df2\u7531 Panel D \u5b8c\u6574\u5448\u73b0\uff0c\u6545\u6539\u5355\u4e00\u7c97\u7bad\u5934 + \u6587\u5b57\uff0c\u65e0\u4fe1\u606f\u635f\u5931\u3002
# 【2026-10-06 改·箭头】原 length=0.16in + linewidth=1.5：实测箭头三角实际视觉高
#   约 14.7pt（11.5pt 多边形 + 3.2pt 同色描边外扩），而本图节点框仅 8.04pt 高 ——
#   箭头几乎是节点框的两倍粗，与全图 0.45–0.55 的线宽语言不一致，显得笨重。
#   改为 length=0.11in + linewidth=0.85（描边外扩随之减半）。
ARW_BIG <- grid::arrow(length = grid::unit(0.11, "inches"), type = "closed")
pA <- pA +
  geom_segment(aes(x = 4.02, y = 6.80, xend = 5.98, yend = 6.80),
               linewidth = 0.85, colour = C_GREY, lineend = "butt",
               arrow = ARW_BIG, inherit.aes = FALSE) +
  # 标签须整体落在最下面一排节点框（UTSW，框底 SVG y=79.35）与 L1 条框（顶 y=99.24）
  # 之间的 19.9pt 空隙内。三项实测约束（SVG 坐标）：
  #   UTSW 框底 79.35 | L1 条框顶 99.24 | 该空隙仅 19.89pt
  # 【2026-10-06 修正】y=4.68/3.80 时："37 organism-drug groups" 字形底 99.6
  #   **压在 L1 条框顶 99.24 上**。按上缘留 1.5pt、行间距 3.0pt、下缘留 1.4pt 重排：
  #   两行基线分别落在 SVG 86.5 / 96.4，对应数据坐标 4.83 / 3.94。
  annotate("text", x = 5, y = 4.83, label = "All 6 directed cross-system pairs",
           size = 2.55, fontface = "bold", colour = "#404040") +
  annotate("text", x = 5, y = 3.94, label = paste0(N_GROUPS, " organism-drug groups"),
           size = 2.45, colour = "#595959")

# ---- 三层 signal：L1 在上、L3 在下（与全文正序一致）----
layer_df <- data.frame(
  layer = c("L1", "L2", "L3"),
  y     = c(3.15, 2.22, 1.29),
  fill  = c("#F4F4F4", "#E8EEF4", "#F8EEE7"),
  text  = c("input shift  \u2014  domain-AUC \u00b7 MMD \u00b7 PSI \u00b7 JS",
            "prediction shift  \u2014  pred mean \u00b7 entropy \u00b7 entropy-JS \u00b7 conf-JS \u00b7 prevalence \u00b7 SUDO",
            "performance estimation  \u2014  transported AUROC: DR (primary) \u00b7 IOW \u00b7 outcome-model"),
  stringsAsFactors = FALSE)
pA <- pA +
  # 【2026-10-06 修正】半高原为 0.31 单位（条高 6.86pt）。说明文字含降部（g/p/y），
  #   且 ggplot 的文字居中**不含降部补偿**，实测文字字形盒比条框低 0.5pt —— 降部
  #   穿出条框下沿 0.32pt（自动检测 L2d 报出）。加高到 0.36（条高 7.97pt）后容纳。
  #   条间距 0.93 单位 ≈ 10.29pt，加高后条间仍余 2.3pt，不会相接。
  geom_rect(data = layer_df, aes(xmin = 1.05, xmax = 8.95,
                                 ymin = y - 0.36, ymax = y + 0.36),
            fill = layer_df$fill, colour = "#9E9E9E", linewidth = 0.34) +
  geom_text(data = layer_df, aes(x = 1.32, y = y, label = layer),
            hjust = 0, size = 2.70, fontface = "bold") +
  geom_text(data = layer_df, aes(x = 2.05, y = y, label = text),
            hjust = 0, size = 2.35, colour = "#262626")

# =============================================================================
# Panel B - task universe（读长表；R 只做标签映射，不做任何统计）
# =============================================================================
b_lab <- c(Candidate = "Candidate\ntransfer tasks",
           Eligible  = "Eligible\n(source AUROC \u2265 0.65)")
dB$stage_key <- factor(dB$stage_key, levels = c("Candidate", "Eligible"))
b_excl <- dB$excluded_note[dB$stage_key == "Eligible"][1]

pB <- ggplot(dB, aes(x = stage_key, y = n)) +
  geom_col(aes(fill = stage_key), width = 0.55, colour = "#333333",
           linewidth = 0.40, show.legend = FALSE) +
  geom_text(aes(label = n), vjust = -0.40, fontface = "bold", size = 3.0) +
  scale_fill_manual(values = c(Candidate = C_LIGHTGREY, Eligible = C_MGB)) +
  scale_x_discrete(labels = b_lab) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(title = "Task universe",
       subtitle = paste0(N_ORG, " organisms  \u00b7  ", N_DRUG, " drugs  \u00b7  ",
                         N_GROUPS, " organism-drug groups"),
       x = NULL, y = "Number of transfer tasks") +
  theme_pub()

# =============================================================================
# Panel C - cross-system degradation（单色柱 + 淡红 harmful 背景）
# =============================================================================
pC <- ggplot(dC, aes(x = d_auroc)) +
  annotate("rect", xmin = -Inf, xmax = TH, ymin = -Inf, ymax = Inf,
           fill = C_RED, alpha = 0.055) +
  geom_histogram(bins = 22, fill = "#7BA7C7", colour = "white", linewidth = 0.25) +
  geom_vline(xintercept = TH, colour = C_RED, linetype = "dashed", linewidth = 0.55) +
  geom_vline(xintercept = P10, colour = "#555555", linetype = "dashed", linewidth = 0.40) +
  geom_vline(xintercept = P50, colour = "#222222", linewidth = 0.55) +
  geom_vline(xintercept = P90, colour = "#555555", linetype = "dashed", linewidth = 0.40) +
  annotate("text", x = P10, y = Inf, label = paste0("P10 ", fmt_s(P10)),
           hjust = 1.05, vjust = 1.25, size = 2.30, colour = "#4D4D4D") +
  annotate("text", x = P90, y = Inf, label = paste0("P90 ", fmt_s(P90)),
           hjust = -0.05, vjust = 1.25, size = 2.30, colour = "#4D4D4D") +
  # 【2026-10-06 修正·其二】上面一步只解决了"压柱子"，但 hjust = 0.5 使标注**居中于
  #   P50 竖线**（实测线在 x=414.57，标注 398.0–431.1）—— 黑线从字中间穿过；
  #   同时红色 harmful 阈值线在 x=398.83，也落在标注跨度内。**两条线都穿字**。
  #   另：P50 原为 fontface="bold" 且 size 2.38，而 P10/P90 为 2.30 非粗体 ——
  #   三条本是同一类分位标注，线型（实线 vs 虚线）已承担层级，粗体属多余且不一致。
  #   现改为：与 P90 同法**右置**（hjust = -0.05），字号/字重与 P10/P90 统一，
  #   颜色对应各自线色（P50 实线 #222222；P10/P90 虚线灰 #555555→字 #4D4D4D）。
  annotate("text", x = P50, y = Inf, label = paste0("P50 ", fmt_s(P50)),
           hjust = -0.05, vjust = 1.25, size = 2.30, colour = "#222222") +
  annotate("text", x = TH, y = Inf, label = paste0("harmful\n\u0394 \u2264 ", sprintf("%.2f", TH)),
           hjust = 1.10, vjust = 8.0, size = 2.25, colour = C_RED) +
  labs(title = "Cross-system degradation",
       subtitle = paste0("harmful = ", N_HARM, "   \u00b7   stable (|\u0394| \u2264 ",
                         sprintf("%.2f", SA), ") = ", N_STABLE,
                         "\nP90 \u2212 P10 = ", sprintf("%.3f", SPREAD)),
       x = "Delta AUROC (target - source internal)", y = "Transfer tasks") +
  scale_x_continuous(expand = expansion(mult = c(0.04, 0.05))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  coord_cartesian(clip = "off") +
  theme_pub()

# =============================================================================
# Panel D - six transfer directions（保留 v1 的 source 配色；顺序用预计算的 direction_order）
# =============================================================================
levD <- dC %>% distinct(direction, direction_order) %>% arrange(direction_order) %>% pull(direction)
dC$direction <- factor(dC$direction, levels = levD)

pD <- ggplot(dC, aes(x = direction, y = d_auroc)) +
  geom_hline(yintercept = TH, colour = C_RED, linetype = "dashed", linewidth = 0.45) +
  geom_hline(yintercept = 0, linewidth = 0.30, colour = "#8C8C8C") +
  geom_boxplot(width = 0.55, outlier.shape = NA, fill = "#DCE7EF",
               colour = "#444444", linewidth = 0.42) +
  geom_point(aes(colour = source_site),
             position = position_jitter(width = 0.13, seed = 1),
             size = 0.85, alpha = 0.75, stroke = 0) +
  scale_colour_manual(values = c(MGB = C_MGB, SU = C_SU, UT = C_UT), name = "source") +
  labs(title = "Directional heterogeneity",
       subtitle = "Ordered by precomputed median  \u00b7  descriptive only; no mechanistic interpretation",
       x = NULL, y = "Delta AUROC") +
  scale_y_continuous(expand = expansion(mult = c(0.07, 0.07))) +
  theme_pub() +
  theme(plot.subtitle = element_text(face = "italic"),
        legend.position = "right",
        axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1, size = 6.6))

# =============================================================================
# Panel E - model / feature class（四组重排；去 legend；标题收敛）
# =============================================================================
CONFIG_LEVELS <- c("P1 History+LR", "S2 History+LGBM", "S1 Core+LR", "S3 Core+LGBM")
dE$config <- factor(dE$config, levels = CONFIG_LEVELS)
sE$.__c__ <- factor(sE$config, levels = CONFIG_LEVELS)
esum <- sE %>% filter(!is.na(.__c__)) %>% arrange(.__c__)

# 轴标签：短名 + 换行（纯标签映射，不涉及统计）
cfg_lab <- c("P1 History+LR" = "History\nLR", "S2 History+LGBM" = "History\nLGBM",
             "S1 Core+LR" = "Core\nLR", "S3 Core+LGBM" = "Core\nLGBM")

# 注释位置为绘图坐标；数值本身全部读自冻结汇总表
ann_y <- max(dE$d_auroc, na.rm = TRUE) +
  0.085 * diff(range(dE$d_auroc, na.rm = TRUE))
esum$ann_y <- ann_y
# 紧凑两行：第一行 = median ΔAUROC（caption 中说明），第二行 = harmful 计数
esum$ann_label <- paste0(fmt_s(esum$median_aligned),
                         "\n", esum$n_harmful_aligned, " harmful")

pE <- ggplot(dE, aes(x = config, y = d_auroc)) +
  geom_hline(yintercept = TH, colour = C_RED, linetype = "dashed", linewidth = 0.45) +
  geom_violin(aes(fill = tier), width = 0.88, trim = FALSE, alpha = 0.50,
              colour = NA, show.legend = FALSE) +
  geom_boxplot(width = 0.22, outlier.shape = NA, fill = "white",
               colour = "#404040", linewidth = 0.43) +
  # 【2026-10-06 修正】原缺 seed —— 抖动每次渲染重抽，导致本图 **PNG/SVG 逐字节不可复现**，
  #   故原先记录的 sha256 对本图不成立。改为与同图另一处（width=0.13）相同的显式
  #   position_jitter(seed = 1) 写法；`geom_jitter(seed=)` 经实测**不生效**（仍不可复现）。
  geom_point(position = position_jitter(width = 0.07, height = 0, seed = 1),
             size = 0.55, alpha = 0.25, colour = "#555555") +
  geom_text(data = esum, aes(x = .__c__, y = ann_y, label = ann_label),
            inherit.aes = FALSE, size = 2.15, lineheight = 0.94, colour = "#333333") +
  scale_fill_manual(values = c(hist = C_HIST, core = C_CORE)) +
  scale_x_discrete(labels = cfg_lab) +
  scale_y_continuous(expand = expansion(mult = c(0.07, 0.20))) +
  labs(title = paste0("Degradation patterns are similar", "\n", "across classifiers"),
       subtitle = paste0("Same P1-defined ", N_ELIG, "-task benchmark universe"),
       x = NULL, y = "Delta AUROC") +
  theme_pub() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 6.8, lineheight = 0.90))

# =============================================================================
# 拼版
# =============================================================================
fig1 <- pA / (pB + pC + plot_layout(widths = c(0.92, 1.08))) /
  (pD + pE + plot_layout(widths = c(1, 1))) +
  plot_layout(heights = c(0.92, 1.12, 1.28)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(family = "Arial", face = "bold", size = 12, colour = "#111111"))

ggsave(file.path(OUT_DIR, "Figure1_v2.png"), fig1, width = 7.2, height = 8.8,
       units = "in", dpi = 400, bg = "white")
ggsave(file.path(OUT_DIR, "Figure1_v2.svg"), fig1, width = 7.2, height = 8.8,
       units = "in", device = svglite::svglite, bg = "white")

cat("WROTE Figure1_v2.{png,svg}\n")
cat("Panel C tasks:", nrow(dC), " Panel E aligned:", nrow(dE), "\n")
cat("harmful:", N_HARM, " stable:", N_STABLE, " P90-P10:", round(SPREAD, 3), "\n")
