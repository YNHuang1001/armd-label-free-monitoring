# =============================================================================
# Supplementary Figure S4 - Pre-unblinding QC of label-free signals
# =============================================================================
# [RENDER LAYER ONLY]
#
# ======================= 代码级盲态防线（本图特有）=========================
# 本脚本**只允许**读取两份 blinded 产物：
#     p6b_label_free_signals_blinded.csv   （L1/L2 signal）
#     p6c_L3_blinded.csv                   （L3 + SUDO）
# 并在运行时执行两道断言：
#   (1) 数据层：两份文件**不得**含任何 target-Y 派生字段
#       （auroc_target / d_auroc / delta_auroc / harmful / stable）
#   (2) 脚本层：本脚本源码**不得**出现任何 unblinded 文件引用
#       （p6d_* / unseal / fig2_forest / fig4_* / fig5_*）
# 这样"这是 blinded QC"不只靠声明，而是运行时可验证。
#
# 本图**不做**任何与 target failure 的关联、分组、着色或 harmful overlay。
# 那些属于 unblinded evaluation，见主文 Figure 2-4。
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2); library(patchwork); library(readr)
  library(dplyr); library(svglite)
})

if (.Platform$OS.type == "windows") {
  try(grDevices::windowsFonts(Arial = grDevices::windowsFont("Arial")), silent = TRUE)
}

.armd_f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
.armd_d <- if (length(.armd_f) && !is.na(.armd_f)) dirname(.armd_f) else "."
ROOT <- Sys.getenv("ARMD_ROOT", unset = normalizePath(file.path(.armd_d, "..", "..")))
RES_DIR <- file.path(ROOT, "\u5206\u6790\u7ed3\u679c")
DATA_DIR <- file.path(ROOT, "figures", "data")
SUP_DIR  <- file.path(ROOT, "figures", "supplement")
dir.create(SUP_DIR, recursive = TRUE, showWarnings = FALSE)

# ---------------- 防线 (2)：脚本层自检 ----------------
# 实现要点：
#   a) 被扫描的源码必须**排除本检查块自身**（否则禁用词表里的字面量会自我触发），
#      用 SELFCHECK-BEGIN/END 哨兵标记；
#   b) 注释行一并排除（注释里说明性地提到这些文件名是允许的）；
#   c) fixed = TRUE 与 ignore.case 不能同用，故改用 tolower() 归一化比较。
SELFCHECK_BEGIN <- "# >>> SELFCHECK-BEGIN"
SELFCHECK_END   <- "# <<< SELFCHECK-END"
# >>> SELFCHECK-BEGIN
# 【2026-10-06 修正】原为相对路径 "R/FigureS4_signal_qc.R"，**只在 cwd = figures/ 时成立**，
#   换任何工作目录都会 readLines 失败并中止。改为从 Rscript 自身取路径，与 cwd 解耦。
.self <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
.self <- if (length(.self) == 1L && nzchar(.self)) .self else "figures/R/FigureS4_signal_qc.R"
.sf <- readLines(.self, warn = FALSE)
# 必须用**整行精确匹配**：赋值行 `SELFCHECK_BEGIN <- "# >>> SELFCHECK-BEGIN"` 自身
# 含有哨兵子串，用 grep(substring) 会把区间错切成 [赋值行, 赋值行]，导致真正的
# 检查块未被排除（曾因此误报 "脚本引用 p6d"）。
.b <- which(trimws(.sf) == SELFCHECK_BEGIN)
.e <- which(trimws(.sf) == SELFCHECK_END)
# 顺序很重要：必须**先按原始行号切除哨兵区间，再剔除注释行**。
# 反过来的话，注释行被删后行号整体前移，哨兵区间会切错位置。
if (length(.b) && length(.e)) .sf <- .sf[seq_along(.sf) < .b[1] | seq_along(.sf) > .e[1]]
.sf <- .sf[!grepl("^\\s*#", .sf)]
.low <- tolower(paste(.sf, collapse = "\n"))
for (f in c("p6d", "unseal", "fig2_forest", "fig4_l3_scatter", "fig4d_", "fig5_subsample")) {
  if (grepl(f, .low, fixed = TRUE)) stop("[盲态防线] 脚本源码中出现被禁文件引用: ", f)
}
rm(.sf, .b, .e, .low)
cat("[盲态防线] 脚本层自检通过：未引用任何 unblinded 文件\n")
# <<< SELFCHECK-END

# ---------------- 防线 (1)：数据层断言 ----------------
B6 <- read_csv(file.path(RES_DIR, "p6b_label_free_signals_blinded.csv"), show_col_types = FALSE)
B3 <- read_csv(file.path(RES_DIR, "p6c_L3_blinded.csv"), show_col_types = FALSE)
FORBIDDEN_COLS <- c("auroc_target", "d_auroc", "delta_auroc", "harmful", "stable")
for (nm in c("B6", "B3")) {
  df <- get(nm)
  hit <- names(df)[vapply(names(df), function(c)
    any(vapply(FORBIDDEN_COLS, function(x) grepl(x, tolower(c), fixed = TRUE), logical(1))), logical(1))]
  if (length(hit)) stop("[盲态防线] ", nm, " 含 target-Y 派生字段: ", paste(hit, collapse = ", "))
}
# 【2026-10-06 修正】原为硬编码 202（旧宇宙）。改为从 fig1b 派生。
UB <- read_csv(file.path(DATA_DIR, "fig1b_universe.csv"), show_col_types = FALSE)
stopifnot(nrow(B6) == UB$eligible_tasks, nrow(B3) == UB$eligible_tasks)
cat("[盲态防线] 数据层断言通过：两份 blinded 文件均无 target-Y 派生字段；各",
    UB$eligible_tasks, "行\n")

SUM <- read_csv(file.path(DATA_DIR, "figS4_signal_summary.csv"), show_col_types = FALSE)
stopifnot(nrow(SUM) == 15, all(SUM$n_finite == UB$eligible_tasks))

# ---------------- 装配长表 ----------------
SPEC <- list(
  c("domain_auc", "L1"), c("mmd", "L1"), c("psi_mean", "L1"),
  c("psi_max", "L1"), c("js_mean", "L1"), c("js_max", "L1"),
  c("pred_mean_shift", "L2"), c("entropy_mean_shift", "L2"), c("entropy_js", "L2"),
  c("conf_dist_js", "L2"), c("pred_pos_prev_shift", "L2"), c("SUDO_AURCC", "L2"))

long <- bind_rows(lapply(SPEC, function(s) {
  src <- if (s[1] %in% names(B6)) B6 else B3
  data.frame(metric = s[1], layer = s[2], value = as.numeric(src[[s[1]]]))
}))
long$layer <- factor(long$layer, levels = c("L1", "L2"))
LAYER_LAB <- c(L1 = "L1  input shift", L2 = "L2  prediction shift")
long$layerf <- factor(LAYER_LAB[long$layer], levels = LAYER_LAB)

LAB <- c(domain_auc = "domain-AUC", mmd = "MMD", psi_mean = "PSI, mean", psi_max = "PSI, max",
         js_mean = "JS, mean", js_max = "JS, max",
         pred_mean_shift = "Prediction mean shift", entropy_mean_shift = "Entropy mean shift",
         entropy_js = "Entropy JS", conf_dist_js = "Confidence-dist. JS",
         pred_pos_prev_shift = "Predicted prevalence shift", SUDO_AURCC = "SUDO-AURCC")
ORD <- names(LAB)
long$lab <- factor(LAB[long$metric], levels = LAB[ORD])
SUM$lab  <- factor(LAB[SUM$metric], levels = LAB[ORD])

# 参考线：仅在同尺度上真正有自然参考点时才画（domain-AUC -> 0.5；shift 类 -> 0）
REFV <- c(domain_auc = 0.5, mmd = 0, psi_mean = 0, psi_max = 0, js_mean = 0, js_max = 0,
          pred_mean_shift = 0, entropy_mean_shift = 0, entropy_js = 0, conf_dist_js = 0,
          pred_pos_prev_shift = 0, SUDO_AURCC = NA_real_)
SUM$ref <- REFV[SUM$metric]
# Panel A 专用的 L1/L2 子集。若把 L3 行也喂给 facet_wrap 的图层，其 lab 为 NA，
# ggplot 会**多开一个 "NA" 面板**（曾出现一格空白的 NA facet）。
SUM_A <- SUM[SUM$layer != "L3", ]

PT <- "#2B2B2B"; REF <- "#B22222"; C_HIST <- "#8FB0CF"

BASE <- theme_classic(base_size = 9, base_family = "Arial") +
  theme(axis.line  = element_line(linewidth = 0.35, colour = "#666666"),
        axis.ticks = element_line(linewidth = 0.30, colour = "#666666"),
        axis.title = element_text(size = 8, colour = "#333333"),
        axis.text  = element_text(size = 7, colour = "#444444"),
        plot.title = element_text(size = 9.5, face = "bold", colour = "#222222"),
        plot.subtitle = element_text(size = 8, colour = "#555555"),
        panel.border = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 7.8, face = "bold", colour = "#444444"),
        panel.spacing = unit(6, "pt"),
        plot.margin = margin(5, 8, 5, 5))

# =============================================================================
# Panel A - 12 个 L1/L2 signal 的原始分布（small multiples，free x）
# =============================================================================
pA <- ggplot(long, aes(x = value)) +
  geom_histogram(bins = 24, fill = C_HIST, colour = "white", linewidth = 0.18) +
  geom_vline(data = SUM_A[!is.na(SUM_A$ref), ], aes(xintercept = ref),
             linetype = "dashed", linewidth = 0.5, colour = REF) +
  geom_vline(data = SUM_A, aes(xintercept = median), linewidth = 0.5, colour = PT) +
  facet_wrap(~ lab, ncol = 4, scales = "free_x") +
  labs(title = "Blinded distributions of the L1 and L2 signals",
       # 【2026-10-06 修正】原为硬编码 202（旧宇宙），改为从 fig1b 派生。
       subtitle = paste0("Native scale; all ", UB$eligible_tasks,
                         " tasks. Black = median, red dashed = natural reference.\n",
                         "Free x-axis: QC distributions, not a comparison of methods."),
       x = "Signal value (native scale)", y = "Transfer tasks") +
  BASE

# =============================================================================
# Panel B - 3 个 transported-AUC estimator 的 blinded 分布（共享 x）
# =============================================================================
L3 <- bind_rows(lapply(c("AUROC_DR", "AUROC_IOW", "AUROC_OM"), function(m)
  data.frame(metric = m, value = as.numeric(B3[[m]]))))
L3$lab <- factor(c(AUROC_DR = "DR (primary)", AUROC_IOW = "IOW (sensitivity)",
                   AUROC_OM = "Outcome model (sensitivity)")[L3$metric],
                 levels = c("DR (primary)", "IOW (sensitivity)", "Outcome model (sensitivity)"))
S3 <- SUM[SUM$layer == "L3", ]; S3$lab <- factor(
  c(AUROC_DR = "DR (primary)", AUROC_IOW = "IOW (sensitivity)",
    AUROC_OM = "Outcome model (sensitivity)")[S3$metric], levels = levels(L3$lab))

pB <- ggplot(L3, aes(x = value)) +
  geom_histogram(bins = 26, fill = C_HIST, colour = "white", linewidth = 0.18) +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.55, colour = REF) +
  geom_vline(data = S3, aes(xintercept = median), linewidth = 0.5, colour = PT) +
  facet_wrap(~ lab, ncol = 3) +
  scale_x_continuous(breaks = c(0.6, 0.7, 0.8, 0.9, 1.0, 1.1)) +
  # 用 coord_cartesian 固定共享轴：scale limits 会把直方图边缘分箱静默丢弃
  # （曾报 "Removed 6 rows (geom_bar())"），coord_cartesian 只裁剪视野、不丢数据。
  coord_cartesian(xlim = c(0.50, 1.12)) +
  labs(title = "Blinded distributions of the transported-AUC estimators",
       # 【2026-10-06 修正】原写 "202 tasks" 且 "Four DR estimates exceed 1"（旧宇宙）。
       #   修正后 DR>1 的 task 数为 0，故改为由数据导出的计数。
       subtitle = paste0("Shared x-axis, ", UB$eligible_tasks,
                         " tasks. Red dashed = AUROC 1.\n",
                         sum(B3$AUROC_DR > 1), " DR estimates exceed 1; none were clipped."),
       x = "Estimated target AUROC", y = "Transfer tasks") +
  BASE

figS4 <- pA / pB + plot_layout(heights = c(1.5, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(family = "Arial", face = "bold", size = 12, colour = "#111111"))

ggsave(file.path(SUP_DIR, "FigureS4.png"), figS4, width = 6.5, height = 7.4,
       units = "in", dpi = 350, bg = "white")
ggsave(file.path(SUP_DIR, "FigureS4.svg"), figS4, width = 6.5, height = 7.4,
       units = "in", device = svglite::svglite, bg = "white")

cat("WROTE FigureS4.{png,svg}\n")
cat("L1/L2 signals:", length(ORD), " L3:", nlevels(L3$lab), "\n")
cat("DR > 1:", sum(L3$value[L3$metric == "AUROC_DR"] > 1, na.rm = TRUE), "\n")
