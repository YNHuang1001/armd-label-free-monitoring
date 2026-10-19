# =============================================================================
# Supplementary Figure S6 - LODO fold composition
# =============================================================================
# [RENDER LAYER ONLY] 只读 figures/data/*.csv，不重算任何统计量。
#   figS6_fold_composition.csv : 38 个 organism-drug group 的 harmful / non-harmful 构成
#   figS6_fold_summary.csv     : 汇总计数（含 Python 层自洽性断言）
#
# 目的：让读者看清**为什么只有 25/38 folds 能算 within-group AUROC**，
#       以及被排除的 13 folds 到底是什么 —— 这个 restriction 不是人为筛选。
#
# 规格（GPT6 锁定）：宽 6.5 in ｜ PNG 350 dpi + SVG ｜ 白底 ｜ Arial ｜
#                   panel label 粗体 ｜ axis title ~9pt ｜ panel title 9-10pt
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
SUP_DIR  <- file.path(ROOT, "figures", "supplement")
dir.create(SUP_DIR, recursive = TRUE, showWarnings = FALSE)

FC <- read_csv(file.path(DATA_DIR, "figS6_fold_composition.csv"), show_col_types = FALSE)
SM <- read_csv(file.path(DATA_DIR, "figS6_fold_summary.csv"), show_col_types = FALSE)
# 【2026-10-06 修正】原为硬编码 38（旧宇宙）。改为从 fig1b 与汇总表派生；
#   并断言逐组表与汇总表自洽（可定义 26 + 单类 11 = 37 折，覆盖 134 + 52 = 186 task）。
UB <- read_csv(file.path(DATA_DIR, "fig1b_universe.csv"), show_col_types = FALSE)
stopifnot(nrow(FC) == UB$n_organism_drug_groups,
          nrow(FC) == SM$n_folds[SM$group == "total"],
          sum(FC$n) == UB$eligible_tasks,
          sum(FC$auroc_definable) == SM$n_folds[SM$group == "evaluable"],
          sum(FC$n[FC$auroc_definable == 1]) == SM$n_tasks[SM$group == "evaluable"])

# ---- 分区：可定义折在前，单类折在后 ----
FC$zone <- ifelse(FC$auroc_definable == 1, "Evaluable folds", "Single-class folds")
FC$zone <- factor(FC$zone, levels = c("Evaluable folds", "Single-class folds"))
FC$orgf <- factor(FC$organism, levels = c("EC", "KP", "PA"))
FC <- FC[order(FC$zone, FC$orgf, FC$drug), ]
FC$lab <- factor(FC$group, levels = rev(FC$group))

# ---- 单类折细分（用于给排除组上色）----
FC$cls <- ifelse(FC$auroc_definable == 1, "Both classes",
          ifelse(FC$n_harmful == 0, "All non-harmful", "All harmful"))
FC$cls <- factor(FC$cls, levels = c("Both classes", "All harmful", "All non-harmful"))

C_BOTH <- "#4C72B0"; C_HARM <- "#C1272D"; C_NONH <- "#D9A441"; C_GREY <- "#5A5A5A"

BASE <- theme_classic(base_size = 9, base_family = "Arial") +
  theme(axis.line  = element_line(linewidth = 0.35, colour = "#666666"),
        axis.ticks = element_line(linewidth = 0.35, colour = "#666666"),
        axis.title = element_text(size = 9, colour = "#333333"),
        axis.text  = element_text(size = 8, colour = "#444444"),
        plot.title = element_text(size = 9.5, face = "bold", colour = "#222222"),
        plot.subtitle = element_text(size = 8, colour = "#555555"),
        panel.border = element_blank(),
        legend.title = element_text(size = 8), legend.text = element_text(size = 8),
        strip.background = element_blank(),
        strip.text.y.left = element_text(angle = 0, hjust = 1, size = 8.5,
                                         face = "bold", colour = "#444444"),
        strip.placement = "outside",
        panel.spacing.y = unit(6, "pt"),
        plot.margin = margin(5, 8, 5, 5))

# =============================================================================
# Panel A - 每个 group 的 harmful / non-harmful 构成
# =============================================================================
# facet 标签必须短：行数多，长标签会把绘图区从左侧压扁（首版即因此把柱体压成一条）
# 【2026-10-06 修正】原为硬编码 "(n = 25)" / "(n = 13)"（旧宇宙）。改为从 FC 派生：
#   先按 zone 取唯一值，再映射成标签，避免再次静默画错。
.ZN <- tapply(FC$zone, FC$zone, length)
ZONE_LAB <- c("Evaluable folds"     = paste0("Evaluable\n(n = ",     .ZN[["Evaluable folds"]], ")"),
              "Single-class folds"  = paste0("Single-class\n(n = ",  .ZN[["Single-class folds"]], ")"))

long <- rbind(
  data.frame(group = FC$group, zone = FC$zone, lab = FC$lab, cls = FC$cls,
             kind = "Harmful", n = FC$n_harmful),
  data.frame(group = FC$group, zone = FC$zone, lab = FC$lab, cls = FC$cls,
             kind = "Non-harmful", n = FC$n_nonharmful))
long$kind <- factor(long$kind, levels = c("Harmful", "Non-harmful"))
tot <- FC %>% group_by(group) %>% summarise(n = sum(n), .groups = "drop")
tot$lab <- factor(tot$group, levels = levels(FC$lab))
tot$zone <- FC$zone[match(tot$group, FC$group)]

pA <- ggplot(long, aes(y = lab, x = n, fill = kind)) +
  geom_col(width = 0.72, colour = "white", linewidth = 0.30) +
  geom_text(data = tot, aes(y = lab, x = n, label = n), inherit.aes = FALSE,
            hjust = -0.35, size = 2.9, colour = "#555555") +
  facet_grid(zone ~ ., scales = "free_y", space = "free_y", switch = "y",
             labeller = labeller(zone = ZONE_LAB)) +
  scale_fill_manual(values = c(Harmful = C_HARM, `Non-harmful` = "#8FB0CF"), name = NULL) +
  scale_x_continuous(limits = c(0, 7.2), breaks = 0:6, expand = c(0, 0)) +
  labs(title = "Composition of each organism-drug group",
       subtitle = "A within-group AUROC requires both transfer classes",
       x = "Number of transfer tasks", y = NULL) +
  BASE +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        legend.key.size = unit(10, "pt"),
        axis.text.y = element_text(size = 6.8))

# =============================================================================
# Panel B - 汇总（folds 与 tasks 两个尺度）
# =============================================================================
gv <- function(grp, col) SM[[col]][match(grp, SM$group)]
seg <- data.frame(
  scale = factor(c(rep("Folds", 3), rep("Tasks", 3)),
                 levels = c("Folds", "Tasks")),
  cls   = factor(c("Both classes", "All harmful", "All non-harmful",
                   "Both classes", "All harmful", "All non-harmful"),
                 levels = c("Both classes", "All harmful", "All non-harmful")),
  n     = c(gv("evaluable","n_folds"), gv("all_harmful","n_folds"), gv("all_nonharmful","n_folds"),
            gv("evaluable","n_tasks"), gv("all_harmful","n_tasks"), gv("all_nonharmful","n_tasks")))
tot_lab <- c(Folds = sum(seg$n[seg$scale == "Folds"]), Tasks = sum(seg$n[seg$scale == "Tasks"]))
lab_seg <- seg %>% filter(n > 0) %>% mutate(txt = sprintf("%d", n))

# folds 与 tasks 量纲相差 5 倍以上，**必须 free_x 分面**。
# 首版把两者画在同一 0–250 轴上，folds(38) 被压成一条线、标签糊在一起。
# 两个坑：
#  1) 默认 position_stack 按 **因子逆序** 堆叠，与图例顺序相反 → 必须 reverse = TRUE
#  2) geom_text 若不带 fill 映射，position_stack 无法按同一分组堆叠，标签会错位糊在一起
pB <- ggplot(seg, aes(y = scale, x = n, fill = cls)) +
  geom_col(width = 0.46, colour = "white", linewidth = 0.4,
           position = position_stack(reverse = TRUE)) +
  # 窄段（如 Tasks 的 All harmful = 6/202）放不下白色数字，故全部标签置于柱体上方。
  # 这里的 cumsum 只是**标签的绘图坐标**，不涉及任何统计量。
  geom_text(data = seg %>%
              group_by(scale) %>%
              mutate(center = cumsum(n) - n / 2) %>%
              ungroup() %>%
              filter(n > 0),
            aes(y = scale, x = center, label = n), inherit.aes = FALSE,
            # 必须 vjust = 0：默认 vjust = 0.5 会让文字下半部分落进柱体上沿，
            # 净空仅约 0.05 数据单位，视觉上就是"数字贴在柱子上"。
            position = position_nudge(y = 0.36), vjust = 0, size = 2.9, colour = "#333333") +
  geom_text(data = data.frame(scale = factor(names(tot_lab), levels = levels(seg$scale)),
                              n = as.numeric(tot_lab),
                              txt = sprintf("total %d", as.numeric(tot_lab))),
            aes(y = scale, x = n, label = txt), inherit.aes = FALSE,
            hjust = -0.06, size = 3.0, colour = "#444444") +
  facet_wrap(~ scale, ncol = 1, scales = "free_x") +
  scale_fill_manual(values = c(`Both classes` = C_BOTH, `All harmful` = C_HARM,
                               `All non-harmful` = C_NONH), name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.16))) +
  # 上方留白：标签用 vjust=0 向上延伸，默认离散 y 范围会把 "138/6/58" 裁掉
  scale_y_discrete(expand = expansion(add = c(0.55, 1.05))) +
  labs(title = "Restriction imposed by single-class folds",
       subtitle = "The macro estimate is computed on the evaluable subset only",
       x = NULL, y = NULL) +
  BASE +
  theme(legend.position = "right", legend.key.size = unit(11, "pt"),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 9, face = "bold", hjust = 0, colour = "#444444"),
        panel.spacing.y = unit(3, "pt"))

# =============================================================================
# 拼版
# =============================================================================
figS6 <- pA / pB + plot_layout(heights = c(4.1, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(family = "Arial", face = "bold", size = 12, colour = "#111111"))

ggsave(file.path(SUP_DIR, "FigureS6.png"), figS6, width = 6.5, height = 9.2,
       units = "in", dpi = 350, bg = "white")
ggsave(file.path(SUP_DIR, "FigureS6.svg"), figS6, width = 6.5, height = 9.2,
       units = "in", device = svglite::svglite, bg = "white")

cat("WROTE FigureS6.{png,svg}\n")
cat("evaluable folds:", gv("evaluable","n_folds"), "tasks:", gv("evaluable","n_tasks"),
    "| excluded folds:", gv("excluded","n_folds"), "tasks:", gv("excluded","n_tasks"), "\n")
