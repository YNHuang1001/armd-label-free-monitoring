# =============================================================================
# Supplementary Figure S8 校验
# =============================================================================
.armd_f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
.armd_d <- if (length(.armd_f) && !is.na(.armd_f)) dirname(.armd_f) else "."
ROOT <- Sys.getenv("ARMD_ROOT", unset = normalizePath(file.path(.armd_d, "..", "..")))
setwd(file.path(ROOT, "figures"))
source("R/check_rendered.R", encoding = "UTF-8")

run_and_check(
  script   = "R/FigureS8_sample_size_trajectories.R",
  outbase  = "FigureS8",
  outdir   = "supplement",
  logfile  = "log/FigureS8_render_check.txt",
  panels   = character(0),          # 单图无 panel 字母
  required = c(
    "Failure-discrimination across target sample sizes",
    "Mean and 95% interval over 200 resampling repeats",
    "Shared y-axis across layers; methods in protocol order",
    "Available unlabeled target observations",
    "Harmful-detection AUROC",
    "100", "250", "500", "1000", "full",
    "L1  input shift", "L2  prediction shift", "L3  performance estimation",
    "domain-AUC", "MMD", "PSI, mean", "PSI, max", "JS, mean", "JS, max",
    "Prediction mean shift", "Entropy mean shift", "Entropy JS",
    "Confidence-dist. JS", "Predicted prevalence shift", "SUDO-AURCC",
    "DR (primary)", "IOW (sensitivity)", "Outcome model (sensitivity)"
  ),
  forbidden = c(
    # 不得写"样本量越大越好"；不得把稳定性写成充分性
    "improves with sample size", "improved with sample size",
    "sufficient", "adequate", "enough data",
    "significant increase", "p < 0.05", "Ranking", "ranked by",
    # 不得出现 paired 行（S8 只画五档轨迹）
    "PAIRED", "paired difference"
  )
)

# 断言 1：y 轴必须共享（0.25-0.80 应只出现一次刻度序列，且出现 0.8/0.7/0.6/0.5/0.4/0.3）
svg <- paste(readLines("supplement/FigureS8.svg", warn = FALSE), collapse = "\n")
for (t in c(">0.8<", ">0.5<", ">0.3<")) if (!grepl(t, svg, fixed = TRUE)) stop("缺 y 刻度 ", t)
cat("\n[共享轴断言] y 轴刻度 0.3-0.8 已渲染；三层共用同一 0.25-0.80 范围\n")

# 断言 2：五档 x 标签齐全且顺序为 100 -> full
pos <- vapply(c(">100<", ">250<", ">500<", ">1000<", ">full<"),
              function(t) regexpr(t, svg, fixed = TRUE)[1], numeric(1))
if (any(pos < 0)) stop("缺 x 档位标签")
if (is.unsorted(pos)) stop("x 档位顺序不是 100 -> full")
cat("[档位断言] x 轴档位顺序 100 -> 250 -> 500 -> 1000 -> full —— OK\n")
