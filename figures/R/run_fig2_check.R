# =============================================================================
# Figure 2 校验
# =============================================================================
.armd_f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
.armd_d <- if (length(.armd_f) && !is.na(.armd_f)) dirname(.armd_f) else "."
ROOT <- Sys.getenv("ARMD_ROOT", unset = normalizePath(file.path(.armd_d, "..", "..")))
setwd(file.path(ROOT, "figures"))
source("R/check_rendered.R", encoding = "UTF-8")

run_and_check(
  script  = "R/Figure2_forest.R",
  outbase = "Figure2",
  panels  = c("A", "B"),
  required = c(
    "Failure discrimination without labels",
    "Rank association", "with ΔAUROC", "Spearman",
    "L1  input shift", "L2  prediction shift", "L3  performance estimation",
    "domain-AUC", "MMD", "PSI, mean", "PSI, max", "JS, mean", "JS, max",
    "Prediction mean shift", "Entropy mean shift", "Entropy JS",
    "Confidence-dist. JS", "Predicted prevalence shift", "SUDO-AURCC",
    "DR   (primary)", "IOW   (sensitivity)", "Outcome model   (sensitivity)",
    "organism-drug cluster bootstrap", "1000", "37",
    "AUROC for detecting harmful transfer"
  ),
  forbidden = c(
    # 禁止"排行榜"式呈现与 post-hoc inversion 措辞。
    # 注意：不可禁用裸词 "Rank" —— Panel B 标题 "Rank association with ΔAUROC"
    # 是合法且必要的措辞（这里是 Spearman 秩相关）。曾因规则过糙而误报。
    "Ranking", "ranking", "ranked by", "Rank of methods",
    "best method", "Best method", "inverted", "Inverted", "reversed MMD"
  )
)

# 顺序断言：15 个方法在 SVG 中必须按 protocol 顺序出现（即 L1 -> L2 -> L3）
svg <- paste(readLines("out/Figure2.svg", warn = FALSE), collapse = "\n")
ORDER_LAB <- c("domain-AUC", "MMD", "PSI, mean", "PSI, max", "JS, mean", "JS, max",
               "Prediction mean shift", "Entropy mean shift", "Entropy JS",
               "Confidence-dist. JS", "Predicted prevalence shift", "SUDO-AURCC",
               "DR   (primary)", "IOW   (sensitivity)", "Outcome model   (sensitivity)")
pos <- vapply(ORDER_LAB, function(x) regexpr(x, svg, fixed = TRUE)[1], numeric(1))
if (any(pos < 0)) stop("顺序断言：找不到标签 ", paste(ORDER_LAB[pos < 0], collapse = ", "))
if (is.unsorted(pos)) {
  stop("顺序断言失败：方法在 SVG 中的出现顺序与 protocol 顺序不一致\n",
       paste(sprintf("  %-32s @%d", ORDER_LAB, pos), collapse = "\n"))
}
cat("\n[顺序断言] 15 个方法按 protocol 顺序出现（L1->L2->L3），非按 AUROC 排序 —— OK\n")
