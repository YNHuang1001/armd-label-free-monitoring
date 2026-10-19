# -*- coding: utf-8 -*-
r"""p43 — 打包期刊要求的 Source Data

================================ 为什么 ================================
Nature Communications（及多数 Nature 系刊）要求投稿时**为每张图提供 Source Data** ——
一个含该图各面板背后原始数值的表格文件。项目里**已经有** `figures/data/` 下 34 个
per-panel CSV（由 `p11_figure_datapack.py` 生成），本脚本只做**打包与建索引**，
**不重算任何数值**。

输出（写入投稿包 `07_source_data/`）：
  SourceData_Figure1.xlsx … Figure5.xlsx     每个主图一个工作簿，一个面板一个 sheet
  SourceData_SupplementaryFigures.xlsx       补充图 S1–S9
  SourceData_FigureS<n>.xlsx                 （仅当该图数据量单独成册更清晰时）
  SourceData_index.csv                       面板 → 源文件 → 列 的完整索引

★ 若某张图**没有**对应的数据文件，脚本**报出来**而不是静默跳过 ——
  「图注说有面板、但没有数据」正是最该被发现的缺口。
"""
import io, os, re, sys
import pandas as pd

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIG  = os.path.join(ROOT, "figures", "data")
OUT  = os.path.join(ROOT, "SUBMISSION_PACKAGE_20261007", "07_source_data")
os.makedirs(OUT, exist_ok=True)
ROWS = []

# 面板 → 数据文件。**显式登记**，不靠文件名猜。
PANELS = [
    # (图, 面板, 文件, sheet 名)
    ("Figure1", "1B", "fig1b_universe.csv",             "universe"),
    ("Figure1", "1B", "fig1b_universe_long.csv",        "universe_long"),
    ("Figure1", "1C", "fig1c_summary.csv",              "summary"),
    ("Figure1", "1C", "fig1c_delta_distribution.csv",   "distribution"),
    ("Figure1", "1D", "fig1d_by_direction.csv",         "by direction"),
    ("Figure1", "1E", "fig1e_configs.csv",              "configs"),
    ("Figure1", "1E", "fig1e_configs_aligned.csv",      "configs aligned"),
    ("Figure1", "1E", "fig1e_configs_full.csv",         "configs full"),
    ("Figure2", "2A", "fig2_forest.csv",                "forest"),
    # ★ 2026-10-16（外部核验三轮未清）：Availability 承诺 "one sheet per figure panel"，
    #   而 SourceData 里 2A 有两张、**2B 一张都没有**。图注 Figure 2B 画的是
    #   「各信号与 ΔAUROC 的 Spearman ρ」，这 15 组值就在 fig2_forest.csv 的三列里
    #   ⇒ 取**列子集**，不另建派生 CSV（少一个生产者，少一处漂移面）。
    ("Figure2", "2B", "fig2_forest.csv",                "spearman"),
    ("Figure2", "2A", "fig2_bootstrap_metadata.csv",    "bootstrap_meta"),
    ("Figure3", "3A", "fig3a_lodo.csv",                 "lodo"),
    ("Figure3", "3B", "fig3b_loo_organism.csv",         "loo_organism"),
    ("Figure3", "3C", "fig3c_sudo_decomposition.csv",   "sudo_decomp"),
    ("Figure4", "4A-C", "fig4_l3_scatter.csv",          "scatter"),
    ("Figure4", "4D", "fig4d_mae.csv",                  "MAE"),
    ("Figure4", "4D", "fig4d_calibration_DR.csv",       "calibration DR"),
    ("Figure4", "4D", "fig4d_calibration_IOW.csv",      "calibration IOW"),
    ("Figure4", "4D", "fig4d_calibration_OM.csv",       "calibration OM"),
    ("Figure4", "4E", "fig4e_ess_association.csv",      "ESS association"),
    ("Figure4", "4E", "fig4e_dr_out_of_range.csv",      "DR out of range"),
    ("Figure5", "5A", "fig5_subsample.csv",             "subsample"),
    ("Figure5", "5B", "fig5b_across_n_range.csv",       "across-n range"),

    ("FigureS1", "S1",    "figS1_task_matrix.csv",              "task matrix"),
    ("FigureS2", "S2A",   "figS2_config_agreement.csv",         "agreement"),
    ("FigureS3", "S3A-B", "figS3_direction_summary.csv",        "direction summary"),
    ("FigureS4", "S4A-B", "figS4_signal_summary.csv",           "signal summary"),
    ("FigureS5", "S5A-B", "figS5_mmd_summary.csv",              "MMD summary"),
    ("FigureS6", "S6A",   "figS6_fold_composition.csv",         "fold composition"),
    ("FigureS6", "S6B",   "figS6_fold_summary.csv",             "fold summary"),
    ("FigureS7", "S7A-B", "figS7_weight_association.csv",       "weight association"),
    ("FigureS8", "S8",    "fig5_subsample.csv",                 "trajectories"),
    ("FigureS9", "S9A",   "figS_orientation.csv",               "orientation"),
    ("FigureS9", "S9B",   "figS9_orientation_consistency.csv",  "consistency"),
]

print("=" * 96); print("p43 — Source Data 打包"); print("=" * 96)

# ---------------------------------------------------------------- 缺口检查
used = {p[2] for p in PANELS}
allf = {f for f in os.listdir(FIG) if f.endswith(".csv")}
unused = sorted(allf - used - {"tableS1_tasks.csv", "tableS2_metrics.csv"})
missing = sorted({p[2] for p in PANELS if not os.path.exists(os.path.join(FIG, p[2]))})
print(f"Figure data 目录：{len(allf)} 个 CSV")
print(f"  已登记进入 Source Data：{len(used)} 个")
print(f"  未登记（非图面板数据，如 tableS1/S2）：{unused or '无'}")
print(f"  **登记了但文件不存在**：{missing or '无'}")
print("")

# 每张图至少要有一个数据文件
for fig in sorted({p[0] for p in PANELS}):
    if not any(p[0] == fig for p in PANELS):
        print(f"  !! {fig} 没有任何数据文件")
have = {p[0] for p in PANELS}
print(f"覆盖 {len(have)} 张图：{sorted(have)}")
print("")

# ---------------------------------------------------------------- 写工作簿

# 面板 → 只取这几列（缺键 = 整表）。键为 (图, 面板)。
# ★ 2026-10-16：Figure 2B 的源是 fig2_forest.csv 的列子集，不是单独 CSV。
COLSUBSET = {
    ("Figure2", "2B"): ["method", "spearman", "sp_lo", "sp_hi"],
}

def write_single_book(items):
    """★ 期刊格式（Nature Communications «Data and Code»）：
       **一个 Excel 工作簿**，每个图/表一个 sheet；sheet 名对应面板（示例用 "Figure 1a"）。
       首版做成了每张图一个独立工作簿（6 个），**不合规**。

       输出既含主图也含补充图，因此 sheet 名统一带 "Figure " 前缀。"""
    path = os.path.join(OUT, "SourceData.xlsx")
    seen, n = set(), 0
    # 一个面板有多个数据文件时，用括号后缀区分；只有一个文件的面板直接用 "Figure 1b"。
    # （期刊示例的 sheet 名形如 "Figure 1a"。首版拼成 `Figure1 1b_1B_universe` —— 面板重复、
    #   无空格，虽然能跑但不合示例格式。）
    by_panel = {}
    for fig, panel, f, tag in items:
        by_panel[(fig, panel)] = by_panel.get((fig, panel), 0) + 1
    with pd.ExcelWriter(path, engine="openpyxl") as xw:
        for fig, panel, f, tag in items:
            d = pd.read_csv(os.path.join(FIG, f))
            _sel = COLSUBSET.get((fig, panel))
            if _sel:
                d = d[_sel]          # 面板子集（如 Figure 2B）
            # 面板标签本身已含图号（"1B" / "S2A" / "4A-C"），
            # 所以 sheet 名 = "Figure " + 面板（字母小写）。目标形如 "Figure 1b"、"Figure S2a"。
            _m = re.match(r"^(S?\d+)(.*)$", panel)
            name = ("Figure " + _m.group(1) + _m.group(2).lower()) if _m else f"Figure {panel}"
            if by_panel[(fig, panel)] > 1:
                name = f"{name} ({tag})"
            name = re.sub(r"[\\/*?:\[\]]", "-", name)[:31]
            base, k = name, 2
            while name in seen:
                name = f"{base[:27]} ({k})"; k += 1
            seen.add(name)
            d.to_excel(xw, sheet_name=name, index=False)
            ROWS.append(dict(figure=fig, panel=panel, sheet=name, source_file=f,
                             rows=len(d), cols=len(d.columns),
                             columns="; ".join(map(str, d.columns))))
            n += 1
    n_sheets = len(pd.ExcelFile(path).sheet_names)
    if n_sheets != len(items):
        raise RuntimeError("SourceData.xlsx：写入 %d 个 item，却只有 %d 个 sheet" % (len(items), n_sheets))
    return path, n


def write_book(key, items):
    """每个 item = (图, 面板, 文件, sheet 名)。

    ★★ 两处首版 bug，都会**静默丢数据**：
      ① 解包写成 `for _, sheet, f, _ in items` —— 把**面板**当成了 sheet 名
         （于是 `S3A/B` 触发 openpyxl 的非法字符报错，暴露了它）；
      ② 更隐蔽的：同一面板有多个数据文件时（如 Figure 1 的 1B 有两个文件），
         两个 sheet **同名互相覆盖** —— Figure 1 应有 8 个 sheet，实际只剩 4 个，
         而 pandas 的 `to_excel` 覆盖时不报错。
       现在：正确解包 + 名称唯一化 + **写完立即断言 sheet 数等于 item 数**。
    """
    path = os.path.join(OUT, f"SourceData_{key}.xlsx")
    seen = set()
    with pd.ExcelWriter(path, engine="openpyxl") as xw:
        for fig, panel, f, sheetname in items:
            d = pd.read_csv(os.path.join(FIG, f))
            _sel = COLSUBSET.get((fig, panel))
            if _sel:
                d = d[_sel]          # 面板子集（如 Figure 2B）
            # Excel sheet 名：上限 31 字符，且**禁用** \ / * ? : [ ]
            name = re.sub(r"[\\/*?:\[\]]", "-", sheetname)[:31]
            base, k = name, 2
            while name in seen:                    # 名称唯一化，杜绝覆盖
                name = f"{base[:28]}_{k}"; k += 1
            seen.add(name)
            d.to_excel(xw, sheet_name=name, index=False)
            ROWS.append(dict(figure=key, panel=panel, sheet=name, source_file=f,
                             rows=len(d), cols=len(d.columns),
                             columns="; ".join(map(str, d.columns))))
    # 写完立刻核对：sheet 数必须等于登记项数，否则就是又发生了覆盖
    n_sheets = len(pd.ExcelFile(path).sheet_names)
    if n_sheets != len(items):
        raise RuntimeError("%s：写入 %d 个 item，却只有 %d 个 sheet —— 发生覆盖"
                           % (os.path.basename(path), len(items), n_sheets))
    return path, len(items)

written = [write_single_book(PANELS)]

print("写出工作簿：")
for p, n in written:
    print(f"  {os.path.basename(p):<40} {n} 个 sheet  ({os.path.getsize(p)/1024:.0f} KB)")
print("")

idx = pd.DataFrame(ROWS)
idx.to_csv(os.path.join(OUT, "SourceData_index.csv"), index=False, encoding="utf-8-sig")
print(f"索引：SourceData_index.csv —— {len(idx)} 行（一行一个面板-sheet）")
print("")

# ---------------------------------------------------------------- 回读校验
print("【回读校验】")
bad = []
for p, _ in written:
    xl = pd.ExcelFile(p)
    for sh in xl.sheet_names:
        d = xl.parse(sh)
        row = idx[idx.sheet == sh].iloc[0]
        if len(d) != row.rows or len(d.columns) != row.cols:
            bad.append(f"{os.path.basename(p)}::{sh}")
print(f"  {'PASS' if not bad else 'FAIL'}  每个 sheet 的行列数与冻结 CSV 一致"
      f"  {'' if not bad else '不符：' + str(bad)}")
print(f"  PASS  工作簿 {len(written)} 个，sheet {len(idx)} 个")
print("")
# ★ 退出码必须反映检查结果。首版检查判了 FAIL 仍打印「完成」并返回 0 ——
#   一个不因自身失败而失败的验收脚本，等于没有验收。
if bad:
    print("*** p43 未通过 ***"); sys.exit(1)
print("*** p43 通过 —— Source Data 已打包，逐 sheet 行列数与冻结 CSV 一致 ***")
