# -*- coding: utf-8 -*-
r"""p44 — 打包投稿用的代码（Code Availability 的可交付件）

================================ 选什么、不选什么 ================================
`分析结果/` 下有 73 个 .py，但**不是都属于这篇论文** —— 其中混有另一个课题的脚本
（`p2_*` 三站点可行性、`p3_01` 多菌种扫描、`tt_01_02` 界标阳性率）。
本脚本按**显式名单**甄别，并把**每一个不收录的文件连同理由**打印出来，
不做静默过滤（静默过滤会让「漏了关键脚本」变得不可见）。

输出到投稿包 `08_code/`：
  analysis/   —— 论文的分析链（按执行顺序编号说明）
  figures/    —— 图渲染 R 脚本
  evs/        —— 估计量验证套件
  README.md   —— 环境、执行顺序、路径说明、以及**未收录清单与理由**
  MANIFEST.csv—— 每个文件的来源路径 + SHA-256
"""
import hashlib, io, os, shutil, sys

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC  = os.path.join(ROOT, "分析结果")
EVS  = os.path.join(ROOT, "evs")
RF   = os.path.join(ROOT, "figures", "R")
OUT  = os.path.join(ROOT, "SUBMISSION_PACKAGE_20261007", "08_code")

# ---- 收录：论文分析链（按执行顺序）
CHAIN = [
    ("episode_definition.py",            "权威 episode 定义（anchor 链），10 项单元测试"),
    ("p4_01_task_availability.py",       "organism–drug 可得性筛选（≥1000 isolates / ≥200 R / ≥200 S）"),
    ("p5_01_degradation_benchmark.py",   "benchmark 构建 + source 模型 + degradation 分布"),
    ("p6b_01_corrected_pipeline.py",     "修正后 signal 流水线（L1/L2）"),
    ("p6c_01_L3_blinded.py",             "L3 估计量（IOW / OM / DR）+ SUDO-AURCC"),
    ("p6d_01_unseal_evaluation.py",      "解封评价（主结局、§7.2 四项指标、bootstrap）"),
    ("p6e_01_lodo_fold_qc.py",           "LODO 折构成 QC"),
    ("p6f_lodo_settlement.py",           "LODO 聚合口径裁定"),
    ("p6g_direction_learning.py",        "方向学习（post-unblinding explanatory）"),
    ("p10a_prepare_cache.py",            "§10 样本量敏感性：缓存准备"),
    ("p10b_subsample_curves.py",         "§10 样本量敏感性：曲线"),
    ("p10c_aggregate.py",                "§10 样本量敏感性：汇总"),
    ("p12_configs_corrected.py",         "四配置（feature tier × classifier）对照"),
    ("p11_figure_datapack.py",           "图数据包（图 1–5 / S1–S9 + Table S1/S2）"),
    ("p6_02_label_free_signals.py",      "**已被 p6b 取代**的早期 signal 流水线；仅作溯源保留"),
]
def _cnt_rows(fn):
    """从审计器自己的 CSV 取断言条数 —— **不手打**。

    ★ 2026-10-15：本条 role 原先写死「153 数值 + 16 定性」，p15 扩到 203 项之后就一直是旧的。
      与 `p40` 的验收表同族：**手工计数的快照会静默过期，而且会跟着代码包一起发出去**。
    """
    import os as _o, csv as _c
    q = _o.path.join(_o.path.dirname(_o.path.abspath(__file__)), fn)
    if not _o.path.exists(q):
        return "?"
    # ★ `csv.reader` **不是上下文管理器**（首版写成 `with _c.reader(...)` ⇒ TypeError）。
    _f = open(q, encoding="utf-8-sig")
    try:
        return str(sum(1 for _ in _c.reader(_f)) - 1)     # 减掉表头
    finally:
        _f.close()


# ---- 收录：稿件审计器
AUDIT = [
    ("p15_results_numeric_audit.py",     "Results 数值一致性（%s 数值 + %s 定性）"
     % (_cnt_rows("p15_results_numeric_audit.csv"), _cnt_rows("p15_results_qualitative_audit.csv"))),
    ("p16_methods_audit.py",             "Methods 数值/语义/完备性（61+14+9）"),
    ("p19_discussion_audit.py",          "Discussion 断言强度"),
    ("p30_mnv_independent_recompute.py", "MNV：正文 headline numbers 独立重算"),
    ("p31_legends_numeric_audit.py",     "补充图例数值审计"),
    ("p32_direction_convention_sensitivity.py", "方向约定敏感性"),
    # ★ 2026-10-09：本会话新增的四件 —— 由 `p51` 的闭包检查发现不在名单里
    #   （「名单外且无理由」4 项，正是 p44 此前静默缺件的同一个缺口复发）。
    ("p52_mae_vs_constant.py",           "MAE 对常数基线的**配对区间**（补 INC-017 缺的区间）"),
    ("p53_rounding_boundary_scan.py",    "同一量多生产者·显示精度扫描（INC-019 发现的仪器）"),
    ("p54_pointer_value_audit.py",       "指针-值判据（INC-020 发现的仪器）"),
    ("rephead.py",                       "审计报告的自我声明头（被 p15/p19/p34/p52/p53/p54 共用）"),
    ("p55_spelling_audit.py",            "拼写/术语一致性（**全稿**，覆盖 8 文件并显式声明未覆盖者）"),
    ("p33_protocol72_metrics.py",        "protocol §7.2 四项指标（含 PPV/NPV）"),
    ("p34_abstract_audit.py",            "Abstract 审计器"),
    ("p35_supplement_inventory_audit.py","补充材料产物↔引用双向清点"),
    ("p21a_reporting_profile_derivation.py", "Reporting profile 推导"),
    ("p21b_reporting_profiles.py",       "Reporting profile 计算"),
    ("p21c_reporting_profiles_finalize.py", "Reporting profile 定稿"),
    ("p21d_render_reporting_tables.py",  "渲染 Supplementary Tables R1–R3"),
    # ★ 2026-10-08 补录：本会话新增的生产/验证脚本此前**不在任何名单** ⇒ 落入「未列入收录名单」，
    #   即代码包**静默陈旧**（投稿要求 Source Code 可复现，缺件即不可复现）。
    ("p21e_render_tableS2.py",           "渲染 Supplementary Tables S1 / S2"),
    ("p21f_layerB_supplement_reconcile.py", "逐 cell 溯源 ＋ 跨路径一致性（含负向测试）"),
    ("p45_triage_policy_simulation.py",  "分诊策略模拟（事后；预注册）"),
    ("p46_render_tableS3.py",            "渲染 Supplementary Table S3"),
    ("p47_universe_partition.py",        "候选/合格宇宙分区审计"),
    ("p48_status_snapshot.py",           "工作清单状态快照（由文件生成）"),
    ("p49_s1_acceptance.py",             "Table S1 可复现性验收（7/37 ＋ S2 Panel A 五行）"),
    ("p50_dependency_edges.py",          "依赖边登记与代次判定（同代/未覆盖/不错序）"),
    ("p50_negtest.py",                   "p50 负向测试（在副本上跑，输出落盘）"),
    ("genstamp.py",                      "内部产物代次头 ＋ 覆盖前自动归档"),
    ("p51_package_census.py",            "包级双向普查 ＋ 闭包检查"),
    ("p44_package_code.py",              "代码打包（本脚本自身）"),
    ("p42_reorganize_freeze.py",         "冻结层重组织"),
    ("p6_00_make_manifest.py",           "构造 p6 清单"),
    ("p6_01_anchor_rerun.py",            "anchor 链重跑"),
    ("p10_preflight_gate.py",            "p10 前置门"),
    ("p21g_render_tripod_v3.py",         "TRIPOD 渲染 v3"),
    ("run_all_audits.py",                "审计批跑驱动"),
    # ★★ 2026-10-12：**又是一次同族复发**（第三次）。
    #   `p57_posthoc_baselines.py` 是事后分析的**唯一入库生产者** ——
    #   正文与摘要里「source-only 基线」的每一个数字都由它产生，
    #   而它此前**不在任何名单** ⇒ 代码包缺件，投稿方**无法复现**这些数。
    #   与 2026-10-08 / 2026-10-09 两次是同一形态：**新脚本写完就忘了登记**。
    #   ⇒ 凡新增生产者/判据脚本，**落盘同一次就要登记**，否则 p51 下一轮才报。
    ("p57_posthoc_baselines.py",
     "事后分析的**唯一入库生产者**（source-only 基线 / 配对差 / 方差分解 / 分诊复核）"),
    # ★ 2026-10-12：INC-035 发现的仪器 —— 参考文献**首现序**。
    #   此前全库有 p15/p16/p19/p22/p24/p35 一批引用相关审计，查的都是「号是否存在、
    #   有无条目、DOI 是否核实」，**没有一条查出现的先后**，而该规则明写在 References 表头里。
    ("p58_reference_order_audit.py",
     "参考文献首现序审计（INC-035 发现的仪器；含 R-1..R-4 与三例判据自测）"),
    # ★ 2026-10-12：p57 的**独立复算**三件（自写实现，未读 p57）。
    #   它们不是"再算一遍"，而是「**换一个实现**再算一遍」—— 是正文里
    #   source-only 基线那批数字**唯一**的独立证据。收进代码包，
    #   否则投稿方拿到的只是 p57 的自证（p57 自己算的数，自己当然对得上）。
    ("indep_verify_p57_auroc_bootstrap.py",
     "**独立**（自写 Mann–Whitney ＋ 自写重采样）复算 p57 的四个 AUROC 与配对差"),
    ("indep_verify_p57_engine_check.py",
     "**独立**复算引擎的一致性检查"),
    ("indep_verify_p57_mcnoise_check.py",
     "**独立**度量蒙特卡洛散布（3 seed × B=10k）＋ B=200k 收敛参考 ＋ 列恒等式"),
    # ★ 2026-10-12：清单 A4.2 的**耐久工件**（该项原为 `✅*`：完成但只有会话证据，按未完成处理）。
    ("p59_survivor_scale_audit.py",
     "幸存者代系审查（A4.2 工件）：四个千分位数过**代系 × 聚合健全性**两轴"),
    # ★★ 2026-10-14：**第四次同型**（p57 / p58 / p59 之后）—— 新脚本写完又忘了登记，
    #   每次都靠 p51 的闭包检查报「名单外且无理由」才发现。四件一起补：
    ("p60_leak_scan.py",
     "投稿包**泄漏闸门**（递归扫文件名/文本/docx-xlsx 内部 XML/文档属性/修订痕迹/隐藏表）"
     "＋ **全称声明专扫**（只列不判，供通读逐条核对覆盖范围）"),
    ("p61_transfer_level_analysis.py",
     "转移级 vs 配对层面：**取向注册表** ＋ 已知答案测试 ＋ 带区间的分解 ＋ 不一致对 ＋ 覆盖/可用性"),
    ("p62_seed_stability.py",
     "端点**种子稳定性**（5 个非相邻种子）；「不得空口归因于蒙特卡洛」规则的第一个实例"),
    ("p63_seed_fleet.py",
     "种子稳定性**总表** —— ★ 状态：**第二遍作废**（INC-037：覆盖了 p6d_01 的种子，"
     "连带改了患者级划分）。第一遍的分诊两条有效。**保留在包内是为了留下这道事故的仪器。**"),
]
# ---- 收录：渲染与交付
RENDER = [
    ("p21e_render_tableS2.py",           "渲染 Supplementary Table S1 / S2"),
    ("p21f_render_supplementary_methods.py", "渲染 Supplementary Methods（估计器设置）"),
    ("p40_build_submission.py",          "投稿包构建（Markdown → Word）"),
    ("p41_verify_submission.py",         "投稿包验收（内容零丢失 / 数字守恒 / 稿件形态）"),
    ("p43_package_source_data.py",       "Source Data 打包"),
    # ★ 2026-10-12：装配**上传用**的 SUBMISSION_20261010/（读 p40 的产物）。
    #   与 p40/p41/p43 同类（交付层），故收在 RENDER 而非 AUDIT。
    ("p56_build_submission_folder.py",   "投稿上传文件夹装配（读 p40 产物，含新鲜度与缺件守卫）"),
]
# ---- 不收录：属于**其他课题**
OTHER_PROJECT = r"^p2_\d|^p3_01_|^tt_01_02_|^p13_|^p14_|^p17_|^p18_|^p20_|^p21_reporting|^p22_|^p23_|^p24_|^p25_|^p26_|^p27_|^p28_|^p29_|^litq\.py$|^ms_master\.py$|^style_map\.py$"
# ★ 2026-10-08（A4.3 闭包检查）：这三个此前落「未列入收录名单」（**无理由**）——
#   闭包检查要求树上每个文件**要么收录、要么带理由排除**，不留无理由者。
# （已并入上一行的正则，避免文本解析器只读首行）
# ---- 不收录：一次性修补脚本
THROWAWAY = r"^_"

def sha(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()

print("=" * 96); print("p44 — 代码打包"); print("=" * 96)

# ★ 2026-10-19：先清空输出目录再写。此前直接覆盖，把 figures/ 下已迁到 figures/R/
#   的扁平旧副本留了下来 —— 归档里混入两份同名脚本（违反「不得混入旧版本」）。
assert os.path.basename(OUT) == "08_code", OUT
shutil.rmtree(OUT, ignore_errors=True)
os.makedirs(OUT, exist_ok=True)
manifest, excluded = [], []

def take(src_path, subdir, note):
    if not os.path.exists(src_path):
        excluded.append((os.path.basename(src_path), "**文件不存在**"))
        return False
    d = os.path.join(OUT, subdir); os.makedirs(d, exist_ok=True)
    dst = os.path.join(d, os.path.basename(src_path))
    shutil.copy2(src_path, dst)
    # ★★ 2026-10-14（泄漏）：`source` 原写 `os.path.relpath(src_path, ROOT)` ⇒
    #   值是 `分析结果/p15_...py` 这类**内部中文目录路径**，会随 `MANIFEST.csv` 进开放仓库。
    #   外部泄漏扫描判红，且它暴露的是本机目录结构。改为**英文的仓库内相对布局**，
    #   与包内的 `analysis/`、`evs/`、`figures/` 三棵子树对齐。
    _rel = os.path.relpath(src_path, ROOT).replace(os.sep, "/")
    for _zh, _en in (("分析结果/", "src/"), ("evs/", "src/evs/"),
                     ("figures/", "figures/"), ("/", "/")):
        if _rel.startswith(_zh):
            _rel = _en + _rel[len(_zh):]
            break
    _key = f"{subdir}/{os.path.basename(src_path)}"
    if any(m["rel"] == _key for m in manifest):
        return True                       # 已在册：不重复登记（p21e 曾同时列在两个名单）
    manifest.append(dict(rel=_key,
                         source=_rel,
                         sha256=sha(src_path), bytes=os.path.getsize(src_path), role=note))
    return True

for fn, note in CHAIN:
    take(os.path.join(SRC, fn), "analysis", note)
for fn, note in AUDIT:
    take(os.path.join(SRC, fn), "analysis", note)
for fn, note in RENDER:
    take(os.path.join(SRC, fn), "analysis", note)

n_r = 0
for fn in sorted(os.listdir(RF)):
    if fn.endswith(".R") and take(os.path.join(RF, fn), "figures/R", "图渲染 / 渲染校验"):
        n_r += 1

n_e = 0
for fn in sorted(os.listdir(EVS)):
    if fn.endswith(".py") and take(os.path.join(EVS, fn), "evs", "估计量验证套件"):
        n_e += 1

# ---- 未收录清单（显式、带理由）
import re as _re
listed = {m["source"].split("/")[-1] for m in manifest}
for fn in sorted(os.listdir(SRC)):
    if not fn.endswith(".py") or fn in listed:
        continue
    if _re.search(THROWAWAY, fn):
        excluded.append((fn, "一次性修补脚本（_ 前缀），非分析链组成部分"))
    elif _re.search(OTHER_PROJECT, fn):
        excluded.append((fn, "属于**另一课题**（三站点可行性 / 报告规范等），与本论文无关"))
    else:
        excluded.append((fn, "未列入收录名单"))

print(f"收录 {len(manifest)} 个文件：analysis {len([m for m in manifest if m['rel'].startswith('analysis')])}"
      f" ｜ figures {n_r} ｜ evs {n_e}")
print("")
print(f"未收录 {len(excluded)} 个（{os.path.join('分析结果')} 下的 .py）：")
for fn, why in excluded:
    print(f"   · {fn:<42} {why}")
print("")

# ---- README
readme = f"""# Analysis code

Code for *Label-free monitoring did not outperform a source-discrimination baseline in a
prespecified cross-system deployment benchmark*.

## Contents

| Directory | What it is |
|---|---|
| `analysis/` | The analysis chain, in the order listed below, plus the manuscript auditors |
| `figures/R/` | R scripts that render Figures 1–5 and Supplementary Figures S1–S9 from the frozen data pack |
| `evs/`      | Estimator validation suite: independent checks of every estimator and evaluation path |

## Execution order (`analysis/`)

| # | Script | Role |
|---|---|---|
""" + "\n".join(f"| {i+1} | `{fn}` | {note} |" for i, (fn, note) in
                enumerate(CHAIN + AUDIT + RENDER)) + """

`p6_02_label_free_signals.py` is retained for provenance only: it is the earlier signal
pipeline that `p6b_01_corrected_pipeline.py` supersedes. The manuscript reports results
from the corrected pipeline.

## Environment

- Python 3.13 with `pandas`, `numpy`, `scipy`, `scikit-learn`, `lightgbm`, `openpyxl`, `python-docx`
- R 4.6 with `ggplot2`, `svglite`, `patchwork`
- Random seeds are fixed and stated in the scripts; every resampling procedure is seeded.

## A note on paths

Scripts contain no machine-specific drive letters: each derives its location from
`__file__`, so `ROOT` / `BASE` / `OUT` follow the script wherever the tree is checked
out. `analysis/` holds the pipeline and audit scripts, `evs/` the estimator-verification
suite, and `figures/R/` the figure scripts. The R scripts honour `ARMD_ROOT` if set.
Three raw microbiology tables are not redistributable; their directory is read from
`ARMD_RAW`.

## Data sources and access

**No data are shipped in this repository.** The three cohorts are obtained from their own
repositories, under the terms each imposes:

| Cohort | Repository | Identifier | Access terms |
|---|---|---|---|
| ARMD-MGB | PhysioNet | `10.13026/2r5-b955` | Credentialed Health Data Licence and signed data use agreement |
| ARMD-Stanford | Dryad | `10.5061/dryad.jq2bvq8kp` | CC0 |
| ARMD-UTSW | Dryad | `10.5061/dryad.0rxwdbsd5` | CC0 |

`BASE` — the directory holding the cohort folders (`ARMD-MGB/`, `ARMD-Stanford/`,
`ARMD-UTSW/`) — is derived as three levels above each script; `evs/t15_drugmap.py`
takes its parent from `ARMD_RAW` when that is set. The pipeline also consumes cached
intermediates, which are not part of this archive.

The task-level analysis dataset is Supplementary Table S1; the figure source data are
provided as Source Data files. Raw electronic health record data cannot be shared
(see Data Availability in the manuscript).
"""
io.open(os.path.join(OUT, "README.md"), "w", encoding="utf-8").write(readme)

# ★ 2026-10-19：随包写 .gitignore。数据一律不进仓库；MANIFEST.csv 例外。
GITIGNORE = """# Data — never committed: raw records are governed by their own agreements
*.csv
!MANIFEST.csv
*.tsv
*.xlsx
*.parquet
*.feather
*.pkl
*.npy
*.npz
*.rds
*.RData
*.sqlite
*.db
data/
raw/
figures/data/
figures/out/

# Environment and secrets
.env
.env.*
*.pem
*.key

# Python
__pycache__/
*.py[cod]

# R
.Rhistory
.Rproj.user/
*.Rproj

# OS / temporary
.DS_Store
Thumbs.db
~$*
*.tmp
"""
io.open(os.path.join(OUT, ".gitignore"), "w", encoding="utf-8", newline="\n").write(GITIGNORE)
# ★ 2026-10-19：锁定字节。MANIFEST.csv 记的是文件字节的 SHA-256；
#   Git for Windows 的 core.autocrlf=true 会在入库时改 CRLF，导致 clone 后哈希全不对。
#   `* -text` 关闭一切换行转换：入库字节 == 磁盘字节。
GITATTR = """# Keep deposited bytes exactly as recorded in MANIFEST.csv (per-file SHA-256).
# Overrides core.autocrlf. Do not edit without regenerating MANIFEST.csv.
* -text
"""
io.open(os.path.join(OUT, ".gitattributes"), "w", encoding="utf-8", newline="\n").write(GITATTR)

import pandas as pd
pd.DataFrame(manifest).to_csv(os.path.join(OUT, "MANIFEST.csv"), index=False, encoding="utf-8-sig")
print(f"写出 README.md 与 MANIFEST.csv（{len(manifest)} 条，含 SHA-256）")
print("")
print("*** p44 完成 ***")
