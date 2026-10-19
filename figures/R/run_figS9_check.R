# =============================================================================
# Supplementary Figure S9 校验
# =============================================================================
.armd_f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
.armd_d <- if (length(.armd_f) && !is.na(.armd_f)) dirname(.armd_f) else "."
ROOT <- Sys.getenv("ARMD_ROOT", unset = normalizePath(file.path(.armd_d, "..", "..")))
setwd(file.path(ROOT, "figures"))
source("R/check_rendered.R", encoding = "UTF-8")

run_and_check(
  script   = "R/FigureS9_orientation.R",
  outbase  = "FigureS9",
  outdir   = "supplement",
  logfile  = "log/FigureS9_render_check.txt",
  panels   = c("A", "B"),
  required = c(
    "Macro discrimination under", "three orientation conventions",
    "Native = the signal's own direction within the held-out group",
    "learned = a univariate rule fitted on the remaining groups.",
    "26 evaluable folds; dashed line = 0.5.",
    "Native direction", "Logistic-learned", "Rank-learned",
    "Orientation selected by the training folds",
    "Proportion of evaluable folds selecting the positive orientation.",
    "Logistic learner", "Rank-based rule",
    "no consistent orientation",
    "domain-AUC", "MMD", "PSI, mean", "PSI, max", "JS, mean", "JS, max",
    "Prediction mean shift", "Entropy mean shift", "Entropy JS",
    "Confidence-dist. JS", "Predicted prevalence shift", "SUDO-AURCC",
    "DR (primary)", "IOW (sensitivity)", "Outcome model (sensitivity)"
  ),
  forbidden = c(
    # 本图是 post-unblinding exploratory：不得写成预注册；不得写成"翻转救回"
    "prespecified primary", "prespecified sensitivity", "pre-registered",
    "rescued", "recovered the signal", "salvage",
    "flipping the sign", "sign correction improved",
    "best orientation", "optimal orientation", "chosen orientation"
  )
)

svg <- paste(readLines("supplement/FigureS9.svg", warn = FALSE), collapse = "\n")

# 断言 1：MMD 的三种口径。
# 注意：图上**不标数值**（数值只是点的坐标、不是文字），故不能在 SVG 里找具体数字。
# 改为核对数据层，并确认三种口径的图例与 15 个方法标签齐全。
# 【2026-10-06 修正】原把三个数值写死为 0.282 / 0.302 / 0.718、折数写死 25（旧宇宙）。
#   修正后为 native 0.447 / logistic 0.409 / rank 0.553，折数 26。改为只断言**内部关系**
#   （rank 规则必然把 MMD 的符号定为负 -> frac_train_dir_pos = 0），不再写死任何数值。
D <- read.csv("data/figS_orientation.csv", fileEncoding = "UTF-8")
SM <- read.csv("data/figS6_fold_summary.csv", fileEncoding = "UTF-8")
mmd <- D[D$method == "mmd", ]
stopifnot(nrow(mmd) == 1,
          mmd$frac_train_dir_pos == 0,                 # rank 规则对 MMD 恒选负方向
          all(is.finite(c(mmd$macro_within, mmd$macro_logit, mmd$macro_rank))))
stopifnot(nrow(D) == 15, all(D$n_eval_fold == SM$n_folds[SM$group == "evaluable"]))
cat(sprintf("\n[数值断言] MMD 三口径 = native %.3f / logistic %.3f / rank %.3f；正方向比例 %.2f（rank 规则全选负方向）\n",
            mmd$macro_within, mmd$macro_logit, mmd$macro_rank, mmd$frac_train_dir_pos))

# 断言 2：层顺序 L1 -> L2 -> L3（用 y 坐标，不用元素顺序）
yof <- function(lab) {
  m <- regmatches(svg, regexpr(sprintf("<text x='[^']*' y='([0-9.]+)'[^>]*>%s", lab), svg))
  if (!length(m)) return(NA_real_)
  as.numeric(sub(".*y='([0-9.]+)'.*", "\\1", m))
}
y1 <- yof("L1  input shift"); y2 <- yof("L2  prediction shift"); y3 <- yof("L3  performance estimation")
if (any(is.na(c(y1, y2, y3)))) stop("找不到层标签")
if (!(y1 < y2 && y2 < y3)) stop("层顺序不是 L1 -> L2 -> L3")
cat(sprintf("[顺序断言] 层顺序 L1(y=%.0f) < L2(y=%.0f) < L3(y=%.0f) —— OK\n", y1, y2, y3))

# 断言 3：Panel B 必须出现 0 与 1 两端刻度（0 = 全部选负方向）
if (!grepl("all negative", svg, fixed = TRUE) || !grepl("all positive", svg, fixed = TRUE))
  stop("Panel B 的 0/1 端点标签缺失")
cat("[一致度断言] Panel B 两端标签（all negative / all positive）齐全\n")
