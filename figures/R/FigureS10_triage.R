# =============================================================================
# Supplementary Figure S10 - Triage policy simulation
# =============================================================================
# [RENDER LAYER ONLY] 只读 figures/data/figS10_triage_*.csv，**不重算任何统计量**。
# 数据源 = p45_triage_policy_simulation.py（事后分析；预注册备忘见项目根目录）。
#
# 设计约束（每条都有理由，不要"顺手改漂亮"）：
#   * **随机线不带置信带。** 曾拟画超几何 95% 带，**已撤回**：
#     固定标签集下随机抽取的命中数服从超几何，其 95% 带 ≈ [7,16]，
#     会让 L-DR(18) 看起来"轻松出带"（+2.65 SD），而**预注册的**聚类 bootstrap
#     下界只有 +0.029。同一画布放两个不确定性模型 = 图替统计说话，
#     且说的比统计多 —— 那正是本文批评的「评估约定悄悄决定结论」的图形版。
#     ⇒ 随机只画**精确期望**（一条线），全部推断只在 k=37 处用预注册 CI 标注。
#   * **六个 L1 成员全部淡显画出来。** 预注册的复合对照落在成员取值区间的**下缘**
#     （7 vs 成员 7–10）。只画复合线会被读成"挑了最不利的聚合"；
#     把成员都摆出来才看得出：无一超过随机期望。
#   * **57 失效天花板**：累积曲线的物理上界，防止读者把"接近随机"看成"接近完美"。
#   * **k=37 竖线 + 明确标注**：推断预注册于此点，其余区间仅描述。
#   * 语言用**富集**（1.59× 基率）而非"检出"：0.487 的精确率意味着
#     被标记队列里**过半不是失效**，必须让读者直接看到这一点。
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2); library(patchwork); library(readr)
  library(dplyr); library(scales)
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

CU  <- read_csv(file.path(DATA_DIR, "figS10_triage_curves.csv"),  show_col_types = FALSE)
MB  <- read_csv(file.path(DATA_DIR, "figS10_triage_members.csv"), show_col_types = FALSE)
SM  <- read_csv(file.path(DATA_DIR, "figS10_triage_summary.csv"), show_col_types = FALSE)

N_TOT <- 186L; N_FAIL <- 57L; K_PREREG <- 37L; BASE <- N_FAIL / N_TOT

stopifnot(nrow(CU) == 4 * N_TOT, nrow(MB) == 6 * N_TOT)
stopifnot(setequal(unique(MB$member),
                   c("domain_auc","mmd","psi_mean","psi_max","js_mean","js_max")))
stopifnot(abs(BASE - 57/186) < 1e-12)
# 天花板与 k=37 处的点数必须与冻结汇总一致
chk <- function(lay) SM$caught[SM$layer == lay]
stopifnot(abs(chk("L-DR") - 18) < 1e-9, abs(chk("L-shift") - 7) < 1e-9,
          abs(chk("L-OM") - 17) < 1e-9)
stopifnot(abs(SM$precision_point[SM$layer == "L-DR_minus_L-shift"] - 11/37) < 1e-12)

LAB <- c("L-DR" = "DR (primary)", "L-OM" = "Outcome model",
         "L-shift" = "Input-shift composite", "L-random" = "Random (exact expectation)")
LAB_M <- c(domain_auc = "domain-AUC", mmd = "MMD", psi_mean = "PSI, mean",
           psi_max = "PSI, max", js_mean = "JS, mean", js_max = "JS, max")

C_DR <- "#3F007D"; C_OM <- "#6A51A3"; C_SH <- "#D94801"
C_MEM <- "#9ECAE1"; C_REF <- "#B22222"; C_CEIL <- "#888888"
# ★ 名字必须**带引号**写成 "L-shift" 形式：`c(L.shift = ...)` 会生成 `L.shift`，
#   与 `$layer` 返回的 "L-shift" 对不上，查找得 NA ⇒ 该层全部行被静默丢弃（曾发生）。
LW <- c("L-shift" = 0.6, "L-OM" = 0.9, "L-DR" = 1.25)

BASE_TH <- theme_classic(base_size = 9, base_family = "Arial") +
  theme(panel.grid.major.y = element_line(linewidth = 0.20, colour = "grey93"),
        panel.grid.major.x = element_blank(),
        axis.line  = element_line(linewidth = 0.35, colour = "#666666"),
        axis.ticks = element_line(linewidth = 0.35, colour = "#666666"),
        axis.title = element_text(size = 9, colour = "#333333"),
        axis.text  = element_text(size = 8, colour = "#444444"),
        plot.title = element_text(size = 9.5, face = "bold", colour = "#222222"),
        plot.subtitle = element_text(size = 8, colour = "#555555"),
        legend.key.size = unit(9, "pt"), legend.text = element_text(size = 7.2),
        legend.title = element_blank(),
        plot.margin = margin(5, 8, 5, 5))

# ---------------------------------------------------------------- Panel A 曲线
CU$lab <- factor(LAB[CU$layer], levels = LAB[c("L-shift", "L-random", "L-OM", "L-DR")])
MB$mlab <- factor(LAB_M[MB$member], levels = LAB_M)

pA <- ggplot() +
  # 六个 L1 成员：淡、细，画在底层
  geom_line(data = MB, aes(x = k, y = cum_failures, group = mlab),
            colour = C_MEM, linewidth = 0.35, alpha = 0.85) +
  # 天花板
  geom_hline(yintercept = N_FAIL, linetype = "dotted", linewidth = 0.5, colour = C_CEIL) +
  # ★ 天花板标注**右对齐**放在右上角：曾放在左侧 (x=4, hjust=0)，
  #   文字横向伸到 x≈62，**与 k=37 的竖线相撞**。右对齐后完全避开竖线，
  #   且仍在全部曲线（最高 57）之上。
  annotate("text", x = N_TOT - 1, y = N_FAIL + 3.4, hjust = 1, size = 2.6,
           colour = "#666666", label = "ceiling: all 57 failing transfers") +
  # 随机：精确期望，**不带带**
  geom_line(data = subset(CU, layer == "L-random"),
            aes(x = k, y = cum_failures, colour = lab),
            linetype = "dashed", linewidth = 0.55) +
  # 三条实体曲线
  geom_line(data = subset(CU, layer %in% c("L-shift", "L-OM", "L-DR")),
            aes(x = k, y = cum_failures, colour = lab),
            linewidth = LW[subset(CU, layer %in% c("L-shift", "L-OM", "L-DR"))$layer]) +
  # 预注册推断点
  geom_vline(xintercept = K_PREREG, linetype = "solid",
             linewidth = 0.45, colour = "#333333") +
  # 标注放**上方空白区**：x=42 处 DR 曲线仅到 y≈21，y=51 一带是空的
  # （曾放在 y=6，与输入漂移复合线在 k=37 的点及曲线相撞）。
  annotate("text", x = K_PREREG + 5, y = 51, hjust = 0, vjust = 1, size = 2.7,
           colour = "#333333",
           label = paste0("k = 37 (20% budget)\npreregistered inference here;\n",
                          "the rest of the range is descriptive")) +
  geom_point(data = subset(CU, k == K_PREREG & layer != "L-random"),
             aes(x = k, y = cum_failures, colour = lab), size = 1.9) +
  scale_colour_manual(values = c("DR (primary)" = C_DR, "Outcome model" = C_OM,
                                 "Input-shift composite" = C_SH,
                                 "Random (exact expectation)" = C_REF),
                      breaks = c("DR (primary)", "Outcome model",
                                 "Input-shift composite", "Random (exact expectation)")) +
  scale_x_continuous(limits = c(0, N_TOT), breaks = seq(0, 180, 30),
                     expand = c(0.01, 0)) +
  scale_y_continuous(limits = c(0, N_FAIL + 8), breaks = seq(0, 50, 10)) +
  labs(title = "Triage: failures captured as the label-review budget grows",
       subtitle = paste0("Six input-shift statistics are shown faintly; ",
                         "none exceeds the random expectation at k = 37."),
       x = "Transfers reviewed, ranked by signal (k)",
       y = "Failing transfers captured") +
  BASE_TH +
  # ★ 图例放**面板外底部**。曾放在面板内 (0.62, 0.28)——那里正是三条曲线
  #   在 k≈100–160 经过的区间，图例把数据盖住了。置于面板外则**不可能**遮挡数据。
  theme(legend.position = "bottom",
        legend.margin = margin(t = -2, b = 2),
        plot.margin = margin(5, 8, 2, 5))

# ---------------------------------------------------------------- Panel B k=37
# ★ 全部 x 必须是**离散水平**。本面板最初把六个成员点画在数值 x=1.62，
#   与因子 x 混用 ⇒ ggplot 把 x 判成连续标度后拒收离散值，整张图渲染失败。
#   现改为把「六个成员」本身作为第一个离散水平，成员点在水平内做确定性抖动。
LEVB <- c("Input-shift metrics\n(individual, n = 6)", "Input-shift\ncomposite",
          "Outcome\nmodel", "DR\n(primary)")
B <- SM[SM$layer %in% c("L-shift", "L-OM", "L-DR"), ]
B$xl <- factor(LEVB[c(2L, 3L, 4L)][match(B$layer, c("L-shift", "L-OM", "L-DR"))],
               levels = LEVB)
B$enr <- B$precision_point / BASE

# 每成员在 k=37 的精确率（淡点，与 A 面板呼应）
MK <- read_csv(file.path(DATA_DIR, "figS10_triage_members_k37.csv"), show_col_types = FALSE)
stopifnot(all(MK$member %in% names(LAB_M)), nrow(MK) == 6)
MK$prec <- MK$caught / K_PREREG
MK$xl <- factor(rep(LEVB[1], nrow(MK)), levels = LEVB)
BASE_TXT <- data.frame(xl = factor(LEVB[1], levels = LEVB), y = BASE + 0.030,
                       lb = sprintf("base rate %.3f", BASE))

pB <- ggplot() +
  geom_hline(yintercept = BASE, linetype = "dashed", linewidth = 0.55, colour = C_REF) +
  geom_text(data = BASE_TXT, aes(x = xl, y = y, label = lb),
            hjust = 0, size = 2.6, colour = C_REF) +
  # ★ 抖动幅度 0.16 → 0.34：六个成员中有四个同为 8/37，原幅度下它们叠在一起，
  #   只看得见 5 个点，与副标题「the six input-shift statistics」不符。
  geom_point(data = MK, aes(x = xl, y = prec), colour = C_MEM, size = 1.7, alpha = 0.95,
             position = position_jitter(width = 0.34, height = 0, seed = 20261007)) +
  geom_errorbar(data = B, aes(x = xl, ymin = lo95, ymax = hi95),
                width = 0.10, linewidth = 0.5, colour = "#333333") +
  geom_point(data = B, aes(x = xl, y = precision_point, colour = xl), size = 2.4) +
  geom_text(data = B, aes(x = xl, y = hi95 + 0.030, colour = xl,
                          label = sprintf("%.3f  (%.2f\u00d7 base)", precision_point, enr)),
            size = 2.6, vjust = 0) +
  scale_colour_manual(values = c("DR\n(primary)" = C_DR, "Outcome\nmodel" = C_OM,
                                 "Input-shift\ncomposite" = C_SH), guide = "none") +
  scale_x_discrete(limits = LEVB, drop = FALSE) +
  scale_y_continuous(limits = c(0, 0.80), breaks = seq(0, 0.8, 0.2)) +
  labs(title = "Precision at the preregistered operating point (k = 37)",
       subtitle = paste0("Points: observed. Bars: 95% organism\u2013drug cluster bootstrap ",
                         "(10,000 replicates).\nFaint points at left: the six ",
                         "input-shift statistics individually \u2014 none exceeds the base rate. ",
                         "Even the best policy\nflags a queue in which more than half ",
                         "the transfers are not failing."),
       # ★ 纵轴标题**折成三行**。原来是一整行 49 字符，旋转后纵向长度超过面板高度，
       #   底端垂到横坐标以下。折行后最长一行 25 字符，纵向占用减半。
       x = NULL, y = "Precision\n(failing transfers in the\nreviewed queue)") +
  BASE_TH +
  theme(axis.text.x = element_text(size = 8.2, colour = "#333333"))

# ★ 分栏标签 A / B：与 S2 / S3 / S6 的既有约定一致（均为 tag_levels = "A"，大写）。
#   正文引用只到图号（`Fig. S10`），不带分栏字母 —— 与其余 9 张补充图的引用方式一致。
P <- pA / pB + plot_layout(heights = c(1.55, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 11, face = "bold", family = "Arial"))

ggsave(file.path(SUP_DIR, "FigureS10.png"), P, width = 6.8, height = 7.4,
       units = "in", dpi = 350, bg = "white")
ggsave(file.path(SUP_DIR, "FigureS10.svg"), P, width = 6.8, height = 7.4,
       units = "in", device = svglite::svglite, bg = "white")

cat("WROTE FigureS10.{png,svg}\n")
cat(sprintf("k=37 caught: DR %d, OM %d, shift %d, members %s\n",
            chk("L-DR"), chk("L-OM"), chk("L-shift"),
            paste(sort(MK$caught), collapse = "/")))
cat(sprintf("enrichment: DR %.3fx, OM %.3fx, shift %.3fx (base %.4f)\n",
            B$enr[B$layer == "L-DR"], B$enr[B$layer == "L-OM"],
            B$enr[B$layer == "L-shift"], BASE))
