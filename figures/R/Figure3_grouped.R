# =============================================================================
# Figure 3 - Generalization beyond the task on which the signal was observed
# =============================================================================
# [RENDER LAYER ONLY] 只读 figures/data/*.csv，不重算任何统计量。
#
# 设计约束（依 GPT6 定稿意见）：
#   A  按 L1/L2/L3 分层，**层内沿用 protocol 顺序**；归类用轻标记表达，**不用排序表达**。
#      图上必须标注 25/38 evaluable folds; 138/202 tasks。
#      below-0.5 的措辞用 "CI entirely below 0.5"，**不得写 "worse than random"**
#      （MMD 存在 post-unblinding 的反向 rank 信息，"worse than random" 会被误读为"无信息"）。
#   B  **不画 3-fold bootstrap CI**（27 种重抽会制造假精确感）。
#      改为：3 个 organism 各自一个点 + 黑色菱形（3 折均值）+ min–max 横线（描述性范围）。
#   C  **不用 waterfall**（易被误读为因果分解）。改为三点连线 + 两段差值标注。
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

A3 <- read_csv(file.path(DATA_DIR, "fig3a_lodo.csv"), show_col_types = FALSE)
B3 <- read_csv(file.path(DATA_DIR, "fig3b_loo_organism.csv"), show_col_types = FALSE)
C3 <- read_csv(file.path(DATA_DIR, "fig3c_sudo_decomposition.csv"), show_col_types = FALSE)

# 【2026-10-06 修正】原为硬编码 25 / 138（旧宇宙）。改为从 figS6_fold_summary.csv 派生。
#   注：fig3c 的列名 pooled_202task / pooled_138task / macro_138task 是**schema 名**，
#   按 STALE_TEXT_INVENTORY_20261006.md 的决定保留不改；其值已是当前宇宙。
FS <- read_csv(file.path(DATA_DIR, "figS6_fold_summary.csv"), show_col_types = FALSE)
N_EVAL_FOLD <- FS$n_folds[FS$group == "evaluable"]
N_EVAL_TASK <- FS$n_tasks[FS$group == "evaluable"]
N_FOLD_TOT  <- FS$n_folds[FS$group == "total"]
N_TASK_TOT  <- FS$n_tasks[FS$group == "total"]

stopifnot(nrow(A3) == 15, nrow(B3) == 15, nrow(C3) == 1)
stopifnot(all(A3$n_eval_fold == N_EVAL_FOLD), all(A3$n_tasks_eval == N_EVAL_TASK))
stopifnot(N_EVAL_FOLD + FS$n_folds[FS$group == "excluded"] == N_FOLD_TOT)

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
# 注意：role 的实际取值是 L1 / L2 / L3-primary / L3-sens（不是裸 L3），
# 直接按名字映射会产生 NA facet。必须先归并到三个层次。
map_role <- function(x) dplyr::case_when(
  x == "L1" ~ "L1", x == "L2" ~ "L2",
  startsWith(x, "L3") ~ "L3", TRUE ~ NA_character_)
ROLE_LAB <- c(L1 = "L1  input shift", L2 = "L2  prediction shift",
              L3 = "L3  performance estimation")
CAT_LAB  <- c(above_chance = "above 0.5", indistinguishable = "incl. 0.5",
              below_chance = "below 0.5")
CAT_COL  <- c(above_chance = "#3C6E9F", indistinguishable = "#8A8A8A",
              below_chance = "#B07A3C")

BASE <- theme_classic(base_size = 8, base_family = "Arial") +
  theme(panel.grid.major.x = element_line(linewidth = 0.20, colour = "grey92"),
        panel.grid.major.y = element_blank(),
        axis.line = element_line(linewidth = 0.35, colour = "#666666"),
        axis.ticks = element_line(linewidth = 0.35, colour = "#666666"),
        axis.title = element_text(size = 8, colour = "#333333"),
        axis.text  = element_text(size = 7, colour = "#444444"),
        plot.title = element_text(size = 9.2, face = "bold", colour = "#222222"),
        plot.subtitle = element_text(size = 7.1, colour = "#555555"),
        panel.border = element_blank(),
        plot.margin = margin(6, 8, 6, 6))
PT <- "#2B2B2B"; REF <- "#B22222"

# =============================================================================
# Panel A - LODO macro AUROC（facet by layer；层内 protocol 顺序）
# =============================================================================
A3 <- A3[match(ORDER, A3$method), ]
# facet 内只保留本层的方法：给 lab 一个「按层嵌套」的因子，free_y 才能正确丢层
A3$layer <- map_role(A3$role)
stopifnot(!any(is.na(A3$layer)))
A3$lab   <- factor(LAB[A3$method], levels = rev(LAB[ORDER]))
A3$rolef <- factor(ROLE_LAB[A3$layer], levels = ROLE_LAB)
A3$catlab <- CAT_LAB[A3$category]

# 右侧归类栏：把 x 轴延伸到 [0.15, 1.06]，数据区在 [0.15, 0.88]，栏位在 0.91
XMAX_DATA <- 0.88
pA <- ggplot(A3, aes(y = lab)) +
  geom_vline(xintercept = 0.5, linetype = "dashed", linewidth = 0.6, colour = REF) +
  geom_vline(xintercept = XMAX_DATA, linewidth = 0.3, colour = "grey85") +
  geom_errorbar(aes(xmin = lo_95, xmax = hi_95), orientation = "y",
                width = 0, linewidth = 0.5, colour = PT) +
  geom_point(aes(x = lodo_macro_auroc), size = 1.6, colour = PT) +
  geom_text(aes(x = 0.905, label = catlab), hjust = 0, size = 2.15,
            colour = CAT_COL[A3$category]) +
  facet_grid(rolef ~ ., scales = "free_y", space = "free_y", switch = "y") +
  # 不要显式写 limits —— 那会强制所有 facet 共用同一组水平，架空 free_y 的丢弃，
  # 表现为每个 facet 都列出全部 15 个方法。顺序由 lab 的 factor levels 决定。
  scale_y_discrete(drop = TRUE) +
  scale_x_continuous(limits = c(0.15, 1.06), breaks = seq(0.2, 0.8, 0.2),
                     expand = c(0, 0)) +
  labs(title = "Grouped generalization (leave-one-organism-drug-out)",
       subtitle = paste0("Macro AUROC over evaluable folds; ",
                         N_EVAL_FOLD, "/", N_FOLD_TOT, " evaluable folds; ",
                         N_EVAL_TASK, "/", N_TASK_TOT, " tasks",
                         "\n95% CI: fold-level cluster bootstrap (descriptive); ",
                         "methods in protocol order within each layer"),
       x = "AUROC (macro mean over evaluable folds)", y = NULL) +
  BASE +
  theme(strip.background = element_blank(),
        strip.text.y.left = element_text(angle = 0, hjust = 1, size = 7,
                                         face = "bold", colour = "#555555"),
        strip.placement = "outside",
        axis.text.y = element_text(size = 7),
        panel.spacing.y = unit(4, "pt"),
        plot.subtitle = element_text(size = 6.8))

# =============================================================================
# Panel B - LOO-organism：3 个真实 held-out organism，**不画 bootstrap CI**
# =============================================================================
B3 <- B3[match(ORDER, B3$method), ]
long <- rbind(
  data.frame(method = B3$method, org = "EC", auc = B3$loo_auc_EC),
  data.frame(method = B3$method, org = "KP", auc = B3$loo_auc_KP),
  data.frame(method = B3$method, org = "PA", auc = B3$loo_auc_PA))
long$lab <- factor(LAB[long$method], levels = rev(LAB[ORDER]))
stopifnot(!any(is.na(long$auc)))
rng <- long %>% group_by(lab) %>% summarise(mn = min(auc), mx = max(auc), .groups = "drop")
B3$lab <- factor(LAB[B3$method], levels = rev(LAB[ORDER]))
rng$lab <- factor(rng$lab, levels = levels(B3$lab))
B3$lab   <- factor(B3$lab, levels = levels(rng$lab))

pB <- ggplot() +
  geom_vline(xintercept = 0.5, linetype = "dashed", linewidth = 0.6, colour = REF) +
  geom_segment(data = rng, aes(y = lab, yend = lab, x = mn, xend = mx),
               linewidth = 0.5, colour = "#9A9A9A") +
  # 【2026-10-06 改】原为实心 16/17/15。Panel B 每行有 3 个 organism × 3 折 = 9 个点，
  #   实测重叠严重（实心互相遮盖，看不出到底几个点）。改为空心 1/2/0：
  #   重叠处能看到轮廓相切，且背后的 min–max 灰线可透出，不干扰读数。
  #   3 折均值仍是**实心**菱形（shape = 18），保持视觉锚点。
  geom_point(data = long, aes(x = auc, y = lab, shape = org),
             size = 1.9, colour = PT, stroke = 0.8) +
  # 【2026-10-06 改】实心菱形 = 3 折均值。原先只写在 Panel B 副标题里，**图例键中缺席**，
  #   读者在图上找不到对应符号。现把它作为第 4 个 shape level 显式加进图例。
  geom_point(data = B3, aes(x = loo_org_macro, y = lab, shape = "3-fold mean"),
             size = 2.9, colour = "black") +
  scale_shape_manual(name = "held-out\norganism",
                     values = c(EC = 1, KP = 2, PA = 0, "3-fold mean" = 18),
                     breaks = c("EC", "KP", "PA", "3-fold mean")) +
  # 必须 unname()：LAB[ORDER] 是带名字的向量，直接把它交给离散 scale 的 limits，
  # ggplot 会用 names() 当轴标签（= 原始方法名）而不是 values（= 显示名）。
  # Panel A 未传 limits（走 factor levels），故不受影响。
  scale_y_discrete(limits = unname(rev(LAB[ORDER])), drop = FALSE) +
  scale_x_continuous(limits = c(0.15, 0.95), breaks = seq(0.2, 0.8, 0.2),
                     expand = c(0, 0)) +
  labs(title = "External stress test\n(leave-one-organism-out)",
       subtitle = paste0("3 held-out-organism folds; black diamond = 3-fold mean,",
                         "\ngrey line = min\u2013max range (descriptive, not a CI)"),
       x = "AUROC on the fully held-out organism", y = NULL) +
  BASE +
  theme(legend.position = "right", legend.title = element_text(size = 6.8),
        legend.text = element_text(size = 6.6), legend.key.size = unit(9, "pt"),
        axis.text.y = element_text(size = 7),
        plot.subtitle = element_text(size = 6.6))

# =============================================================================
# Panel C - SUDO 聚合分解：三点连线 + 两段差值（**不用 waterfall**）
# =============================================================================
pts <- data.frame(
  x = 1:3,
  y = c(C3$pooled_202task, C3$pooled_138task, C3$macro_138task),
  # 轴标签只保留一眼可读的三步；完整术语放 caption（括号挤进轴上会互相重叠）
  lab = c(paste0(N_TASK_TOT, " tasks\npooled OOF"),
          paste0(N_EVAL_TASK, " tasks\npooled OOF"),
          paste0(N_EVAL_TASK, " tasks\nmacro LODO")),
  # 【2026-10-06 修正】数值标签的避让。
  #   实测（SVG 坐标）：折线 P1(412.63,437.58)->P2(451.10,399.87)->P3(489.56,441.42)。
  #   原 vjust 全为 -1.0（标签在点上方）时：
  #     "0.332" 字形盒 x 403.4–421.9 / y 426.7–433.9 —— 上升段在 x=421.9 处 y=428.5，
  #             **穿盒**；
  #     "0.325" 字形盒 x 480.3–498.8 / y 430.6–437.8 —— 下降段在 x=486.2–489.6 处穿盒。
  #   故第 1、3 点的标签改放**点的下方**（折线自 P1 向上、至 P3 结束，下方恒空），
  #   中间峰值点仍在上方并左移。轴在 y=468.1，两个下置标签底约 450–454，仍有余量。
  vhj = c(0.5, 1.15, 0.5),
  vvj = c(2.0, -1.0, 2.0))
# 差值标注：显示相邻两点之间的**实际增量**（即后来的值 − 先前的值）。
# 注意 fig3c 中 task_subset_effect / aggregation_effect 存的是"前 − 后"的相反数，
# 此处显示的是读者在图上看到的方向，并在方法中注明约定。
d12 <- pts$y[2] - pts$y[1]
d23 <- pts$y[3] - pts$y[2]
# \u6807\u6ce8\u4f4d\u7f6e\u987b\u907f\u5f00\u6570\u503c\u6807\u7b7e\u4e0e\u8fde\u7ebf\uff1a\u7b2c\u4e00\u6bb5\u7f6e\u4e8e\u7ebf\u4e0a\u65b9\uff1b\u7b2c\u4e8c\u6bb5\u53f3\u5bf9\u9f50\u3001\u8d34\u5728\u659c\u7ebf\u5de6\u4fa7
# \u30102026-10-06 \u4fee\u6b63\u00b7\u91cd\u8981\u3011\u539f y \u8f74 limits \u786c\u7f16\u7801 0.44\u20130.84\u3001\u4e24\u4e2a\u5dee\u503c\u6807\u6ce8\u653e\u5728 y=0.610/0.700\uff0c
#   \u90fd\u662f\u4e3a\u65e7\u503c (0.511 / 0.529 / 0.756) \u8bbe\u7684\u3002\u4fee\u6b63\u540e\u4e09\u70b9\u4e3a 0.332 / 0.408 / 0.325\uff0c
#   **\u5168\u90e8\u843d\u5728 0.44 \u4ee5\u4e0b**\uff0c\u4e8e\u662f Panel C \u7684\u753b\u7ebf/\u70b9/\u6570\u503c\u6807\u7b7e\u88ab\u6574\u5e45\u4e22\u5f03
#   \uff08\u6e32\u67d3\u65e5\u5fd7 "Removed 3 rows containing missing values or values outside the scale range"\uff09\uff0c
#   \u5373\u8be5\u9762\u677f\u6b64\u524d\u662f**\u7a7a\u767d**\u7684\u3002\u6b64\u5904\u6539\u4e3a\u7531 pts$y \u6d3e\u751f\u8f74\u8303\u56f4\u4e0e\u6807\u6ce8\u4f4d\u7f6e\uff0c\u675c\u7edd\u518d\u6b21\u53d1\u751f\u3002
YR   <- range(pts$y)
PADY <- diff(YR)
# 上界必须 ≥ 0.5，否则 geom_hline(yintercept = 0.5)（随机参考线）会被丢弃
YLIM <- c(YR[1] - 0.50 * PADY, max(YR[2] + 0.95 * PADY, 0.52))
# \u4f4d\u7f6e\u7ea6\u675f\uff08\u4f9d SVG \u5b9e\u6d4b\u5750\u6807\uff09\uff1a\u4e09\u70b9\u4e3a P1(412.6,437.6) P2(451.1,399.9) P3(489.6,441.4)\uff0c
#   \u6570\u503c\u6807\u7b7e 0.408 \u5360 x 429.8\u2013448.3 / y 389\u2013395\u3002\u539f\u4e24\u6bb5\u6807\u6ce8 x \u8de8\u5ea6 407\u2013457 \u4e0e 412\u2013463
#   **\u51e0\u4e4e\u91cd\u5408**\uff0c\u4e14\u540e\u8005\u7684 y\u2248389\u2013397 \u6b63\u538b\u5728 0.408 \u4e0a\u3002
# \u30102026-10-06 \u518d\u4fee\u6b63\u3011\u4e0a\u4e00\u7248\u628a\u4e24\u6761\u6807\u6ce8\u5e76\u6392\u7f6e\u4e8e\u540c\u4e00\u9ad8\u5ea6\uff0c\u4f46\u90a3**\u653e\u4e0d\u4e0b**\uff1a
#   \u5b9e\u6d4b Panel C \u7ed8\u56fe\u533a\u4ec5 x 408.02\u2013504.95\uff08\u5bbd 96.9pt\uff09\uff0c\u800c "task-subset effect" 50pt\u3001
#   "aggregation effect" 51pt \u2014\u2014 \u5e76\u6392\u9700 >101pt\uff0c\u4e8e\u662f\u4e24\u6761**\u5de6\u53f3\u4e24\u7aef\u90fd\u88ab\u9762\u677f\u88c1\u5207**
#   \uff08\u56fe\u4e0a\u8bfb\u6210 "sk-subset effect" / "aggregation effec"\uff09\u3002
#   \u6539\u4e3a\uff1a**\u4e0a\u4e0b\u5806\u53e0\u3001\u7edf\u4e00\u5de6\u5bf9\u9f50**\u5230\u7ed8\u56fe\u533a\u5de6\u7f18\u9644\u8fd1\uff08x = 0.90 \u6570\u636e\u5355\u4f4d\uff09\uff0c
#   \u4e24\u6761\u5404\u5360\u4e00\u884c\u5757\uff0c\u5747\u843d\u5728\u6298\u7ebf\u5cf0\u9876 (y\u2248399.9) \u4e4b\u4e0a\uff0c\u4e92\u4e0d\u76f8\u78b0\u3002
segs <- data.frame(x = c(0.90, 0.90),
                   y = c(YR[2] + 0.91 * PADY, YR[2] + 0.44 * PADY), hj = c(0, 0),
                   lab = c(sprintf("task-subset effect\n\u0394 = %+.3f", d12),
                           sprintf("aggregation effect\n\u0394 = %+.3f", d23)))

pC <- ggplot(pts, aes(x = x, y = y)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", linewidth = 0.6, colour = REF) +
  geom_line(linewidth = 0.7, colour = PT) +
  geom_point(size = 2.3, colour = PT) +
  geom_text(aes(label = sprintf("%.3f", y), hjust = vhj, vjust = vvj),
            size = 2.6, fontface = "bold") +
  geom_text(data = segs, aes(x = x, y = y, label = lab, hjust = hj), inherit.aes = FALSE,
            size = 2.25, colour = "#555555", lineheight = 0.95, vjust = 0.5) +
  scale_x_continuous(breaks = 1:3, labels = pts$lab, limits = c(0.6, 3.4),
                     expand = c(0, 0)) +
  scale_y_continuous(limits = YLIM, breaks = seq(round(YLIM[1], 1), round(YLIM[2], 1), by = 0.05)) +
  labs(title = "Aggregation changes\nthe estimate",
       subtitle = "SUDO-AURCC: same signal,\nthree reporting conventions",
       x = NULL, y = "AUROC") +
  BASE +
  # 轴标签为两行多词文本，Panel C 绘图区窄，字号必须压到能容下最宽一行
  # （曾因 6.8pt 下 "pooled OOF" ≈ 0.47in > 刻度间距 0.45in 而相邻标签相撞）
  theme(axis.text.x = element_text(size = 5.7, lineheight = 0.90),
        plot.subtitle = element_text(size = 6.9))

# =============================================================================
# 拼版
# =============================================================================
# B 的绘图区本身很空（0.2–0.9 范围、点集中在 0.3–0.8），把宽度匀给 C 以容纳轴标签
fig3 <- pA / (pB + pC + plot_layout(widths = c(1.05, 1.25))) +
  plot_layout(heights = c(1.55, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(family = "Arial", face = "bold", size = 12, colour = "#111111"))

ggsave(file.path(OUT_DIR, "Figure3.png"), fig3, width = 7.2, height = 7.0,
       units = "in", dpi = 400, bg = "white")
ggsave(file.path(OUT_DIR, "Figure3.svg"), fig3, width = 7.2, height = 7.0,
       units = "in", device = svglite::svglite, bg = "white")

cat("WROTE Figure3.{png,svg}\n")
cat(sprintf("SUDO decomposition: %.4f -> %.4f -> %.4f  (d12=%+.4f d23=%+.4f)\n",
            pts$y[1], pts$y[2], pts$y[3], d12, d23))
