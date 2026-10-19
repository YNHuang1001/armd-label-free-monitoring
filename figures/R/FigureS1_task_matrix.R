# =============================================================================
# Supplementary Figure S1 - Complete cross-system task matrix
# =============================================================================
# [RENDER LAYER ONLY] 只读 figures/data/figS1_task_matrix.csv（Python 层已装配并断言）
#
# 定位：228 个 candidate transfer 的完整 task landscape，
#       以及最终 202 个 benchmark task 落在哪里、各自 ΔAUROC 如何。
#
# 关键设计（依 GPT6 定稿意见）：
#   * 行 = 38 个 organism–drug group，按 EC -> KP -> PA、组内按药物固定顺序，
#     **不按 ΔAUROC 排序**（排序方式本身即信息，不做事后择优排序）。
#   * 列 = 6 个 directed transfer，固定顺序。
#   * fill = ΔAUROC，发散色标**中心固定 0**（对称 limits，避免两侧强度不可比）。
#   * 不 eligible 的 26 格：灰色 + ×，**不是空白** —— 让 228->202 的原因一目了然。
#   * 不额外用重颜色编码 harmful（二分类已由图 1C 及后续全部分析承担，
#     同一 tile 同时承担连续量与二分类会过杂）。
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2); library(readr); library(svglite)
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

M <- read_csv(file.path(DATA_DIR, "figS1_task_matrix.csv"), show_col_types = FALSE)

# ---- 复核 Python 层的自洽性断言（R 层不得放宽）----
# 【2026-10-06 修正】原为硬编码 228/202/26/38（旧宇宙）。改为从 fig1b_universe.csv 派生。
#   「候选组数」用 nrow(M)/6（每组恰 6 个方向），因为本图矩阵画的是**候选**网格（38 组），
#   而 fig1b 的 n_organism_drug_groups 是**合格**组数（37）——两者不是同一个量。
UB <- read_csv(file.path(DATA_DIR, "fig1b_universe.csv"), show_col_types = FALSE)
stopifnot(nrow(M) == UB$candidate_tasks,
          sum(M$eligible) == UB$eligible_tasks,
          sum(!M$eligible) == UB$excluded)
stopifnot(sum(M$eligible & is.na(M$d_auroc)) == 0)
stopifnot(sum(!M$eligible & !is.na(M$d_auroc)) == 0)
stopifnot(length(unique(M$group)) == nrow(M) / 6, length(unique(M$direction)) == 6)

# ggplot 把 factor 第 1 水平画在**底部**；要实现"EC 在顶、PA 在底"必须 reverse。
M$grp_f <- factor(M$group, levels = rev(unique(M$group[order(M$group_rank)])))
M$dir_f <- factor(M$direction, levels = unique(M$direction[order(M$dir_rank)]))
M$fillv <- ifelse(M$eligible, M$d_auroc, NA_real_)
# 行标签附 eligible 计数（轻量 summary，用轴标签承载，不额外占版面）
rowlab <- tapply(M$eligible, M$group, sum)
levels(M$grp_f) <- sprintf("%-8s %d/6", levels(M$grp_f), as.integer(rowlab[levels(M$grp_f)]))

LR <- max(abs(M$d_auroc), na.rm = TRUE)          # 对称 limits，保证 0 居中且两侧强度可比
NEG <- "#B2182B"; MID <- "#F7F7F7"; POS <- "#2166AC"

p <- ggplot(M, aes(x = dir_f, y = grp_f)) +
  geom_tile(aes(fill = fillv), colour = "white", linewidth = 0.30) +
  geom_text(data = M[!M$eligible, ], aes(label = "\u00d7"),
            size = 2.9, colour = "#5A5A5A") +
  scale_fill_gradient2(low = NEG, mid = MID, high = POS, midpoint = 0,
                       limits = c(-LR, LR), na.value = "#D9D9D9",
                       name = expression(Delta*AUROC),
                       breaks = c(-0.2, -0.1, 0, 0.1, 0.2)) +
  scale_x_discrete(position = "top") +
  labs(title = "Complete cross-system task matrix",
       # \u6298\u6210\u4e09\u884c\uff1a\u5355\u884c\u8fc7\u957f\u4f1a\u88ab\u753b\u5e03\u53f3\u8fb9\u754c\u622a\u65ad
       # \u30102026-10-06 \u4fee\u6b63\u3011\u539f\u4e3a\u786c\u7f16\u7801 228/38/202/26\uff08\u65e7\u5b87\u5b99\uff09\uff0c\u6539\u4e3a\u4ece M \u4e0e fig1b \u6d3e\u751f\u3002
       subtitle = paste0(
         nrow(M), " candidate cells: ", nrow(M) / 6, " organism-drug groups \u00d7 6 directed transfers.\n",
         sum(M$eligible), " eligible (coloured); ", sum(!M$eligible), " excluded (grey \u00d7). ",
         "Fill = observed \u0394AUROC (target \u2212 source internal).\n",
         "Row and column order fixed by organism, drug and direction, not by outcome.")) +
  theme_minimal(base_size = 9, base_family = "Arial") +
  theme(axis.title = element_blank(),
        axis.text.x.top = element_text(size = 8, colour = "#333333"),
        axis.text.y = element_text(size = 6.4, colour = "#333333", family = "mono"),
        panel.grid = element_blank(),
        plot.title = element_text(size = 10, face = "bold", colour = "#222222"),
        plot.subtitle = element_text(size = 7.4, colour = "#555555"),
        legend.position = "right", legend.key.height = unit(26, "pt"),
        legend.title = element_text(size = 8), legend.text = element_text(size = 7.5),
        plot.margin = margin(6, 8, 20, 6))

# 底部每列 eligible 计数（轻量）
colcnt <- tapply(M$eligible, M$direction, sum)
ymin_num <- 0.45
p <- p +
  annotate("text", x = levels(M$dir_f), y = ymin_num,
           label = sprintf("%d", as.integer(colcnt[levels(M$dir_f)])),
           size = 2.6, colour = "#666666", vjust = 1) +
  annotate("text", x = levels(M$dir_f), y = ymin_num - 0.55,
           label = "eligible / 38", size = 2.2, colour = "#888888", vjust = 1) +
  coord_cartesian(ylim = c(0, 38.5), clip = "off")

ggsave(file.path(SUP_DIR, "FigureS1.png"), p, width = 6.5, height = 7.6,
       units = "in", dpi = 350, bg = "white")
ggsave(file.path(SUP_DIR, "FigureS1.svg"), p, width = 6.5, height = 7.6,
       units = "in", device = svglite::svglite, bg = "white")

cat("WROTE FigureS1.{png,svg}\n")
cat("cells:", nrow(M), " eligible:", sum(M$eligible), " ineligible:", sum(!M$eligible), "\n")
cat("symmetric fill limits: +/-", round(LR, 3), "\n")
