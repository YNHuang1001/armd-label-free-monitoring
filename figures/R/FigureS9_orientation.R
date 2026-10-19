# =============================================================================
# Supplementary Figure S9 - Post-unblinding analysis of learnable signal orientation
# =============================================================================
# [RENDER LAYER ONLY] 只读 figures/data/figS_orientation.csv（列已齐备，无需新增分析）
#
# 定位：**post-unblinding explanatory analysis**，不是预注册 primary 或 sensitivity。
# 该分析承担的关键限定：AUROC < 0.5（在其预设方向上）**不等于**该信号不含信息 ——
# native 读法可能失效，但某些信号（尤其 MMD）含有**可迁移的反向秩信息**。
#
# 设计（GPT6 定稿）：
#   A 15 方法 × 三种方向口径（native / logistic-learned / rank-learned）的 grouped discrimination，
#     按 L1->L2->L3、protocol 顺序，不按结果排序；0.5 参考线。
#   B orientation consistency：每个可定义折中选择**正方向**的比例
#     （0 = 25/25 选负方向；1 = 全部选正方向）。
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
DATA_DIR <- file.path(ROOT, "figures", "data")
SUP_DIR  <- file.path(ROOT, "figures", "supplement")
dir.create(SUP_DIR, recursive = TRUE, showWarnings = FALSE)

D <- read_csv(file.path(DATA_DIR, "figS_orientation.csv"), show_col_types = FALSE)
CN <- read_csv(file.path(DATA_DIR, "figS9_orientation_consistency.csv"), show_col_types = FALSE)
# 【2026-10-06 修正】原为硬编码 25（旧宇宙）。改为从折汇总表派生；
#   并断言同折数、两个 learner 的方向比例表与 AUC 表方法集合一致。
SM <- read_csv(file.path(DATA_DIR, "figS6_fold_summary.csv"), show_col_types = FALSE)
stopifnot(nrow(D) == 15,
          all(D$n_eval_fold == SM$n_folds[SM$group == "evaluable"]),
          nrow(CN) == nrow(D), all(CN$n_folds == D$n_eval_fold),
          setequal(CN$method, D$method))

ORDER <- c("domain_auc", "mmd", "psi_mean", "psi_max", "js_mean", "js_max",
           "pred_mean_shift", "entropy_mean_shift", "entropy_js", "conf_dist_js",
           "pred_pos_prev_shift", "SUDO_AURCC",
           "AUROC_DR", "AUROC_IOW", "AUROC_OM")
LAB <- c(domain_auc = "domain-AUC", mmd = "MMD", psi_mean = "PSI, mean", psi_max = "PSI, max",
         js_mean = "JS, mean", js_max = "JS, max",
         pred_mean_shift = "Prediction mean shift", entropy_mean_shift = "Entropy mean shift",
         entropy_js = "Entropy JS", conf_dist_js = "Confidence-dist. JS",
         pred_pos_prev_shift = "Predicted prevalence shift", SUDO_AURCC = "SUDO-AURCC",
         AUROC_DR = "DR (primary)", AUROC_IOW = "IOW (sensitivity)",
         AUROC_OM = "Outcome model (sensitivity)")
map_role <- function(x) dplyr::case_when(
  x == "L1" ~ "L1  input shift", x == "L2" ~ "L2  prediction shift",
  startsWith(x, "L3") ~ "L3  performance estimation", TRUE ~ NA_character_)

D <- D[match(ORDER, D$method), ]
D$layer <- factor(map_role(D$role),
                  levels = c("L1  input shift", "L2  prediction shift",
                             "L3  performance estimation"))
# 层内序号 -> 数值 y（便于同一方法下并列三个口径）
D$idx <- ave(seq_len(nrow(D)), D$layer, FUN = seq_along)
D$ypos <- D$idx

PT <- "#2B2B2B"; REF <- "#B22222"
C_NAT <- "#3A3A3A"; C_LOG <- "#C98A2B"; C_RNK <- "#2E6DA4"

BASE <- theme_classic(base_size = 9, base_family = "Arial") +
  theme(panel.grid.major.x = element_line(linewidth = 0.20, colour = "grey92"),
        panel.grid.major.y = element_blank(),
        axis.line  = element_line(linewidth = 0.35, colour = "#666666"),
        axis.ticks = element_line(linewidth = 0.35, colour = "#666666"),
        axis.title = element_text(size = 8.2, colour = "#333333"),
        axis.text  = element_text(size = 7.2, colour = "#444444"),
        plot.title = element_text(size = 9.5, face = "bold", colour = "#222222"),
        plot.subtitle = element_text(size = 7.4, colour = "#555555"),
        panel.border = element_blank(),
        strip.background = element_blank(),
        strip.text.y.left = element_text(angle = 0, hjust = 1, size = 7.8,
                                         face = "bold", colour = "#555555"),
        strip.placement = "outside",
        panel.spacing.y = unit(5, "pt"),
        legend.key.size = unit(9, "pt"), legend.text = element_text(size = 7.2),
        plot.margin = margin(5, 8, 5, 5))

# =============================================================================
# Panel A - 三种方向口径下的 grouped discrimination
# =============================================================================
# 【修正记录】初版用**连续 y 轴 + 显式 breaks**，结果 15 个标签被画进**每一个 facet**
# （连续 scale 的 breaks 是全局的，facet 无法各自过滤），整图不可读。
# 改为**离散因子 y + position_dodge**：facet 的 free_y 才能正常丢弃层外方法。
A <- rbind(
  data.frame(method = D$method, layer = D$layer, lo = D$within_lo,
             est = D$macro_within, hi = D$within_hi, ori = "Native direction"),
  data.frame(method = D$method, layer = D$layer, lo = D$logit_lo,
             est = D$macro_logit, hi = D$logit_hi, ori = "Logistic-learned"),
  data.frame(method = D$method, layer = D$layer, lo = D$rank_lo,
             est = D$macro_rank, hi = D$rank_hi, ori = "Rank-learned"))
A$ori <- factor(A$ori, levels = c("Native direction", "Logistic-learned", "Rank-learned"))
A$lab <- factor(LAB[A$method], levels = rev(LAB[ORDER]))

pA <- ggplot(A, aes(y = lab, x = est, colour = ori, shape = ori, group = ori)) +
  geom_vline(xintercept = 0.5, linetype = "dashed", linewidth = 0.6, colour = REF) +
  geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y",
                width = 0, linewidth = 0.42, position = position_dodge(width = 0.62)) +
  geom_point(size = 1.6, position = position_dodge(width = 0.62)) +
  facet_grid(layer ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_colour_manual(values = c(`Native direction` = C_NAT,
                                 `Logistic-learned` = C_LOG,
                                 `Rank-learned` = C_RNK), name = NULL) +
  scale_shape_manual(values = c(`Native direction` = 16,
                                `Logistic-learned` = 17,
                                `Rank-learned` = 15), name = NULL) +
  scale_x_continuous(limits = c(0.15, 0.90), breaks = seq(0.2, 0.8, 0.2)) +
  labs(title = "Macro discrimination under\nthree orientation conventions",
       # 措辞必须精确：native 那一列**不做任何学习**，把整图统称 LODO 会误导。
       # 只有 logistic / rank 两列是"在其余组学方向、在留出组评价"。
       # 三行：两行会被右边界截断（宽度受左侧长方法标签挤压）
       # 【2026-10-06 修正】原副标题硬编码 "25 evaluable folds"（旧宇宙），改为从折汇总表派生。
       subtitle = paste0("Native = the signal's own direction within the held-out group;\n",
                         "learned = a univariate rule fitted on the remaining groups.\n",
                         SM$n_folds[SM$group == "evaluable"], " evaluable folds; dashed line = 0.5."),
       x = "Macro AUROC over evaluable folds", y = NULL) +
  BASE +
  guides(colour = guide_legend(nrow = 1), shape = guide_legend(nrow = 1)) +
  theme(legend.position = "bottom", legend.box.margin = margin(0, 0, 0, 0))

# =============================================================================
# Panel B - orientation consistency（**两个 learner 并列**）
# =============================================================================
# 【修正记录】初版只画 frac_train_dir_pos，它来自 **rank 规则**（mean(rank_dir > 0)），
# 不含 logistic learner；而 Panel A 同时有 logistic 与 rank 两个 learner，
# 只画一个会让读者以为该列同时代表两者。现并列两个 learner。
B <- rbind(
  data.frame(method = CN$method, layer = D$layer[match(CN$method, D$method)],
             frac = CN$frac_rank_pos, learner = "Rank-based rule"),
  data.frame(method = CN$method, layer = D$layer[match(CN$method, D$method)],
             frac = CN$frac_logit_pos, learner = "Logistic learner"))
B$learner <- factor(B$learner, levels = c("Logistic learner", "Rank-based rule"))
B$lab <- factor(LAB[B$method], levels = rev(LAB[ORDER]))

pB <- ggplot(B, aes(y = lab, x = frac, colour = learner, shape = learner, group = learner)) +
  geom_vline(xintercept = 0.5, linetype = "dashed", linewidth = 0.5, colour = "#9A9A9A") +
  geom_line(aes(group = lab), linewidth = 0.4, colour = "#C8D2DA") +
  geom_point(size = 1.9, position = position_dodge(width = 0.5)) +
  facet_grid(layer ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_colour_manual(values = c(`Logistic learner` = C_LOG,
                                 `Rank-based rule` = C_RNK), name = NULL) +
  scale_shape_manual(values = c(`Logistic learner` = 17, `Rank-based rule` = 15), name = NULL) +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1),
                     labels = c("0\n(all negative)", "0.5", "1\n(all positive)")) +
  labs(title = "Orientation selected by the training folds",
       subtitle = paste0("Proportion of evaluable folds selecting the positive orientation.\n",
                         "0.5 = no consistent orientation across folds."),
       x = "Proportion of folds selecting positive orientation", y = NULL) +
  BASE +
  guides(colour = guide_legend(nrow = 1), shape = guide_legend(nrow = 1)) +
  theme(legend.position = "bottom", axis.text.y = element_text(size = 7.2))

figS9 <- pA / pB + plot_layout(heights = c(1.55, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(family = "Arial", face = "bold", size = 12, colour = "#111111"))

ggsave(file.path(SUP_DIR, "FigureS9.png"), figS9, width = 6.5, height = 8.6,
       units = "in", dpi = 350, bg = "white")
ggsave(file.path(SUP_DIR, "FigureS9.svg"), figS9, width = 6.5, height = 8.6,
       units = "in", device = svglite::svglite, bg = "white")

cat("WROTE FigureS9.{png,svg}\n")
cat("methods:", nrow(D), " folds:", unique(D$n_eval_fold), "\n")
cat("MMD: native", D$macro_within[D$method == "mmd"],
    " logistic", D$macro_logit[D$method == "mmd"],
    " rank", D$macro_rank[D$method == "mmd"], "\n")
