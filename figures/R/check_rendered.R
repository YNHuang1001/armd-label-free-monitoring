# =============================================================================
# check_rendered.R - 通用渲染校验（Figure 2-5 复用）
# =============================================================================
# 存在理由：本轮 Figure 1 反复出现"代码跑通、图是错的"：
#   - geom_curve(curvature=<vector>) 直接杀死设备 -> 半空图
#   - 循环内 aes(y = arcs$y[i]) 踩 lazy-eval -> 6 条线叠成 1 条
#   - Arial 未注册 -> PNG 静默 fallback，预览与产物不一致
#   - scale_y_continuous(limits=) 静默丢点
# 这些都不会让 Rscript 返回非零退出码。
#
# 三层校验：
#   L1 SVG 内容检查（文件存在/非空、必备文字、禁用旧文字、panel 字母齐全）
#   L2 render log 检查（致命模式 -> stop；允许白名单）
#   L3 人工成图检查（不可自动化，脚本只打印清单）
#
# 用法：
#   source("R/check_rendered.R")
#   run_and_check(
#     script   = "R/Figure2_forest.R",
#     outbase  = "Figure2",
#     required = c("..."), forbidden = c("..."),
#     panels   = c("A","B","C"))
# =============================================================================

# --- 致命模式（默认）---
FATAL_DEFAULT <- c(
  "Removed [0-9]+ rows",
  "values outside the scale range",
  "font family not found",
  "Killing locked device",
  "infinite values",
  "non-finite",
  "Error in ",
  "Execution halted",
  "Warning message"
)

# --- 允许白名单（正则，匹配到的行不算致命）---
ALLOW_DEFAULT <- c(
  "^Loading required package",
  "^Attaching package",
  "masked from",
  "^Warning: package .* was built under"
)

.say <- function(...) cat(..., "\n", sep = "")

#' 运行 R 脚本并执行三层校验
#'
#' @param script   相对 figures/ 的脚本路径
#' @param outbase  输出文件名前缀（不含扩展名），如 "Figure2"
#' @param required SVG 中必须出现的字符串向量
#' @param forbidden SVG 中禁止出现的字符串向量（旧标题/旧术语）
#' @param panels   期望存在的 panel 字母，如 c("A","B","C")
#' @param outdir   输出目录（默认 figures/out）
#' @param logfile  render log 路径（默认 out/<outbase>_render.log）
#' @param fatal    额外致命模式
#' @param allow    额外白名单
#' @param rscript  Rscript 可执行文件路径
#' @return 不可见地返回 TRUE；任一层失败则 stop()
run_and_check <- function(script, outbase,
                          required = character(0), forbidden = character(0),
                          panels = character(0),
                          outdir = "out", logfile = NULL,
                          fatal = character(0), allow = character(0),
                          rscript = "Rscript") {

  figdir <- normalizePath(".", mustWork = TRUE)
  svg <- file.path(figdir, outdir, paste0(outbase, ".svg"))
  png <- file.path(figdir, outdir, paste0(outbase, ".png"))
  if (is.null(logfile)) logfile <- file.path(figdir, outdir, paste0(outbase, "_render.log"))

  fatal <- unique(c(FATAL_DEFAULT, fatal))
  allow <- unique(c(ALLOW_DEFAULT, allow))

  fails <- character(0)

  # ---------------- L2 渲染 + log 检查 ----------------
  .say("\n== [L2] render: ", script, " ==")
  res <- system2(rscript, shQuote(file.path(figdir, script)),
                 stdout = logfile, stderr = logfile)
  logtxt <- if (file.exists(logfile)) readLines(logfile, warn = FALSE) else character(0)
  cat(paste(logtxt, collapse = "\n"), "\n")

  if (!identical(as.integer(res), 0L))
    fails <- c(fails, sprintf("脚本退出码 = %s（非 0）", res))

  # 逐行判定：命中致命模式且未被白名单豁免
  hit_fatal <- function(ln) {
    if (any(vapply(allow, function(a) grepl(a, ln, perl = TRUE), logical(1)))) return(FALSE)
    any(vapply(fatal, function(f) grepl(f, ln, perl = TRUE), logical(1)))
  }
  bad <- logtxt[vapply(logtxt, hit_fatal, logical(1))]
  if (length(bad)) {
    fails <- c(fails, sprintf("render log 命中致命模式 %d 行", length(bad)))
    .say("  致命行：")
    for (b in utils::head(bad, 8)) .say("    | ", b)
  } else {
    .say("  无致命模式。")
  }

  # ---------------- L1 文件 + SVG 内容 ----------------
  .say("\n== [L1] 输出与 SVG 内容 ==")
  for (f in c(svg, png)) {
    if (!file.exists(f)) { fails <- c(fails, paste("缺失文件:", basename(f))); next }
    sz <- file.info(f)$size
    if (is.na(sz) || sz < 1000) fails <- c(fails, sprintf("%s 过小(%s bytes)", basename(f), sz))
    else .say(sprintf("  %-28s %8.1f KB", basename(f), sz / 1024))
  }

  if (file.exists(svg)) {
    s <- paste(readLines(svg, warn = FALSE), collapse = "\n")
    miss <- required[!vapply(required, function(x) grepl(x, s, fixed = TRUE), logical(1))]
    if (length(miss)) {
      fails <- c(fails, sprintf("SVG 缺少必备文字 %d 项", length(miss)))
      for (m in miss) .say("    缺: ", m)
    } else .say(sprintf("  必备文字 %d 项全部存在。", length(required)))

    pres <- forbidden[vapply(forbidden, function(x) grepl(x, s, fixed = TRUE), logical(1))]
    if (length(pres)) {
      fails <- c(fails, sprintf("SVG 含禁用(旧)文字 %d 项", length(pres)))
      for (m in pres) .say("    禁用却出现: ", m)
    } else if (length(forbidden)) .say(sprintf("  禁用文字 %d 项均未出现。", length(forbidden)))

    if (length(panels)) {
      # patchwork 的 tag 渲染为 <text ...>A</text>
      lack <- panels[!vapply(panels, function(p)
        grepl(paste0(">", p, "</text>"), s, fixed = TRUE), logical(1))]
      if (length(lack)) {
        fails <- c(fails, sprintf("SVG 缺少 panel 字母: %s", paste(lack, collapse = ",")))
      } else .say(sprintf("  panel 字母 %s 齐全。", paste(panels, collapse = "/")))
    }

    # ---------------- L2b 文字越界检测（2026-10-06 新增，自动）----------------
    # 起因：本轮修正后 Figure 2 Panel B 标题右端 522.66 > 画布 518.40（+4.26），
    #       Figure 4 Panel E 标题右端 560.39 > 518.40（+41.99）——两处都被画布**裁切**，
    #       而 L1（文字存在）与 L2（无致命模式）**全部通过**。故把「无裁切」从
    #       纯人工 L3 提升为可自动检测的 L2b。
    # 依据：svglite 为每个 <text> 写出 textLength，anchor 为 start/middle/end 三态，
    #       故右边界 = x + textLength（start）/ x + len/2（middle）/ x（end），可精确求得。
    #
    # 【字形盒口径】降部只在字符串**确实含降部字母**时才计（g j p q y 及 , ; ( ) [ ] { } |）。
    #   否则 "MGB"/"L1"/"37 organism-drug groups" 这类会被平白加上 0.21em 的下沿，
    #   造成假报（首版即因此把节点框内的 "MGB" 误判为与框相交）。
    # 不用正则：R 的 TRE 对 [ ] { } 在字符类内的转义挑剔（实测直接报错）。
    # 逐字符 %in% 判断更稳。
    .DESCCH <- c("g", "j", "p", "q", "y", ",", ";", "(", ")", "[", "]", "{", "}",
                 "|", "/", "@", "Q", "J")
    .isdesc <- function(txt) {
      if (is.na(txt) || !nzchar(txt)) return(FALSE)
      any(strsplit(txt, "")[[1]] %in% .DESCCH)
    }
    .glyph <- function(tg, body) {
      x  <- suppressWarnings(as.numeric(sub(".*[^a-z-]x='([-0-9.]+)'.*", "\\1", tg)))
      y  <- suppressWarnings(as.numeric(sub(".*[^a-z-]y='([-0-9.]+)'.*", "\\1", tg)))
      ln <- suppressWarnings(as.numeric(sub(".*textLength='([0-9.]+)px'.*", "\\1", tg)))
      fs <- suppressWarnings(as.numeric(sub(".*font-size:\\s*([0-9.]+)px.*", "\\1", tg)))
      if (any(is.na(c(x, y, ln)))) return(NULL)
      if (is.na(fs)) fs <- 8
      anc <- if (grepl("text-anchor='middle'", tg, fixed = TRUE)) "middle" else
             if (grepl("text-anchor='end'",    tg, fixed = TRUE)) "end"    else "start"
      L <- switch(anc, end = x - ln, middle = x - ln / 2, start = x)
      d <- if (.isdesc(body)) 0.21 else 0.02
      c(L, L + ln, y - 0.78 * fs, y + d * fs)
    }
    .textel <- "<text[^>]*>[^<]*</text>"
    svgtag <- regmatches(s, regexpr("<svg[^>]*>", s))
    cv <- if (length(svgtag))
      suppressWarnings(as.numeric(sub(".*width='([0-9.]+)pt'.*", "\\1", svgtag))) else NA_real_
    if (!is.na(cv)) {
      over <- character(0)
      for (tg in regmatches(s, gregexpr(.textel, s))[[1]]) {
        gb <- .glyph(sub(">.*", ">", tg), sub(".*?>([^<]*)</text>", "\\1", tg))
        if (is.null(gb)) next
        if (gb[2] > cv + 0.01)
          over <- c(over, sprintf("右端 %.2f > 画布 %.2f (+%.2f)", gb[2], cv, gb[2] - cv))
      }
      if (length(over)) {
        fails <- c(fails, sprintf("SVG 有 %d 处文字越出画布右边界（会被裁切）", length(over)))
        for (o in over) .say("    越界: ", o)
      } else .say("  文字均在画布内（越界检测 L2b）。")
    }

    # ---------------- L2c 面板裁切检测（2026-10-06 新增，自动）----------------
    # 起因：Figure 3 Panel C 的两条标注并排同高，合计宽 >101pt，而该面板绘图区仅 96.9pt，
    #       于是**左右两端被面板裁切**（图上读成 "sk-subset effect" / "aggregation effec"），
    #       而 L1 / L2 / L2b **全部通过** —— L2b 只查画布外沿，查不到更内层的面板边界。
    # 依据：svglite 为每个 panel 写一个 <clipPath id='cp...'><rect .../></clipPath>，
    #       并把内容包在 <g clip-path='url(#cp...)'> 中；按 <g>/</g> 维护栈，
    #       即可得到每个 <text> 所属的**最内层**裁剪矩形，从而判定是否被裁。
    .cps <- list()
    for (blk in regmatches(s, gregexpr("<clipPath id='[^']+'>[^<]*<rect[^>]*>", s))[[1]]) {
      id <- sub(".*id='([^']+)'.*", "\\1", blk)
      v  <- function(k) suppressWarnings(as.numeric(sub(paste0(".*", k, "='([-0-9.]+)'.*"), "\\1", blk)))
      x <- v("x"); y <- v("y"); w <- v("width"); h <- v("height")
      if (!any(is.na(c(x, y, w, h)))) .cps[[id]] <- c(x, x + w, y, y + h)
    }
    if (length(.cps)) {
      .toks <- regmatches(s, gregexpr(paste0("<g[^>]*>|</g>|", .textel), s))[[1]]
      .st <- list(); .clip <- character(0)
      for (tg in .toks) {
        if (startsWith(tg, "</g>")) { if (length(.st)) .st[[length(.st)]] <- NULL; next }
        if (startsWith(tg, "<g")) {
          id <- sub(".*clip-path='url\\(#([^)]+)\\)'.*", "\\1", tg)
          .st[[length(.st) + 1L]] <- if (!identical(id, tg)) .cps[[id]] else NULL
          next
        }
        r <- NULL
        for (k in rev(seq_along(.st))) if (!is.null(.st[[k]])) { r <- .st[[k]]; break }
        if (is.null(r)) next
        body <- sub(".*?>([^<]*)</text>", "\\1", tg)
        gb <- .glyph(sub(">.*", ">", tg), body)
        if (is.null(gb)) next
        L <- gb[1]; R <- gb[2]; top <- gb[3]; bot <- gb[4]
        if (L < r[1] - 0.05 || R > r[2] + 0.05 || top < r[3] - 0.05 || bot > r[4] + 0.05) {
          .clip <- c(.clip, sprintf("%s  [L%.1f R%.1f] vs 面板 [%.1f,%.1f]",
                                    substr(body, 1, 34), L, R, r[1], r[2]))
        }
      }
      if (length(.clip)) {
        fails <- c(fails, sprintf("SVG 有 %d 处文字被所属面板裁切", length(.clip)))
        for (o in .clip) .say("    面板裁切: ", o)
      } else .say("  文字均在其面板内（面板裁切检测 L2c）。")
    }

    # ---------------- L2d 文字×图元 部分重叠检测（2026-10-06 新增，自动）----------------
    # 起因：Figure 1 Panel A 的 "N organism-drug groups" 字形底 99.6 压在 L1 条框顶
    #       99.24 上（重叠 0.36pt）。这是**文字 vs 矩形**的碰撞，L2b（画布）与
    #       L2c（面板裁剪框）都查不到，只能靠人眼。
    # 判据：文字**完全落在**矩形内 = 有意为之（如 "L1" 在色条内、方程式文字在框内）；
    #       文字与矩形**部分相交** = 缺陷。背景大矩形（面积 ≥ 面板 70%）不参与比较。
    .rects <- list()
    for (tg in regmatches(s, gregexpr("<rect[^>]*>", s))[[1]]) {
      x <- suppressWarnings(as.numeric(sub(".*[^a-z-]x='([-0-9.]+)'.*", "\\1", tg)))
      y <- suppressWarnings(as.numeric(sub(".*[^a-z-]y='([-0-9.]+)'.*", "\\1", tg)))
      w <- suppressWarnings(as.numeric(sub(".*[^a-z-]width='([-0-9.]+)'.*", "\\1", tg)))
      h <- suppressWarnings(as.numeric(sub(".*[^a-z-]height='([-0-9.]+)'.*", "\\1", tg)))
      if (any(is.na(c(x, y, w, h))) || w <= 0 || h <= 0) next
      .rects[[length(.rects) + 1L]] <- c(x, x + w, y, y + h)
    }
    # 【收窄】只与"内容框"比较：面板背景、画布底等**大矩形**一律排除，
    #   否则轴标题/副标题贴着面板边缘会被大量误报（首版 7/14 图误报，已弃）。
    #   判据：高 ≤ 45pt 且面积 ≤ 画布 5% —— 如 Figure 1 的 L1/L2/L3 色条、
    #   节点小框；而面板背景动辄高 > 100pt。
    .panelA <- 0
    for (r in .rects) .panelA <- max(.panelA, (r[2] - r[1]) * (r[4] - r[3]))
    .rects <- Filter(function(r) (r[4] - r[3]) <= 45 &&
                              (r[2] - r[1]) * (r[4] - r[3]) <= 0.05 * .panelA, .rects)
    .col <- character(0)
    for (tg in regmatches(s, gregexpr(.textel, s))[[1]]) {
      gb <- .glyph(sub(">.*", ">", tg), sub(".*?>([^<]*)</text>", "\\1", tg))
      if (is.null(gb)) next
      L <- gb[1]; tb <- gb
      for (r in .rects) {
        ox <- min(tb[2], r[2]) - max(tb[1], r[1])
        oy <- min(tb[4], r[4]) - max(tb[3], r[3])
        if (ox <= 0.15 || oy <= 0.15) next                              # 不相交
        inside <- tb[1] >= r[1] - 0.15 && tb[2] <= r[2] + 0.15 &&
                  tb[3] >= r[3] - 0.15 && tb[4] <= r[4] + 0.15
        if (!inside) {
          .col <- c(.col, sprintf("文字[x %.1f–%.1f, y %.1f–%.1f]  与图元[x %.1f–%.1f, y %.1f–%.1f] 部分相交 %.1fx%.1f",
                                  tb[1], tb[2], tb[3], tb[4], r[1], r[2], r[3], r[4], ox, oy))
        }
      }
    }
    if (length(.col)) {
      fails <- c(fails, sprintf("SVG 有 %d 处文字与图元部分重叠", length(.col)))
      for (o in .col) .say("    重叠: ", o)
    } else .say("  文字与图元无部分重叠（L2d）。")
  }

  # ---------------- L3 人工清单 ----------------
  .say("\n== [L3] 人工成图检查（不可自动化，必须逐条看图确认）==")
  .say("  !! L1/L2 通过 **不能** 替代 L3。已发生实例：Figure 3 Panel C 的 x 轴标签")
  .say("     因绘图区过窄而相邻相撞，L1（文字存在）与 L2（无致命模式）均完全通过。")
  for (q in c("无裁切（标题/轴标签/注释完整）",
              "无文本碰撞",
              "图例与数据一致",
              "每个 panel 都非空",
              "箭头/连线方向正确",
              "颜色与符号未因 factor level 变化而错配",
              "参考线/阈值线位置正确")) .say("  [ ] ", q)

  # ---------------- 结论 ----------------
  if (length(fails)) {
    .say("\n*** 校验失败 ***")
    for (f in fails) .say("  - ", f)
    stop("check_rendered: ", length(fails), " 项未通过", call. = FALSE)
  }
  .say("\n*** 校验通过（L1+L2）。L3 仍需人工确认。***")
  invisible(TRUE)
}
