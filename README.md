# Analysis code

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
| 1 | `episode_definition.py` | 权威 episode 定义（anchor 链），10 项单元测试 |
| 2 | `p4_01_task_availability.py` | organism–drug 可得性筛选（≥1000 isolates / ≥200 R / ≥200 S） |
| 3 | `p5_01_degradation_benchmark.py` | benchmark 构建 + source 模型 + degradation 分布 |
| 4 | `p6b_01_corrected_pipeline.py` | 修正后 signal 流水线（L1/L2） |
| 5 | `p6c_01_L3_blinded.py` | L3 估计量（IOW / OM / DR）+ SUDO-AURCC |
| 6 | `p6d_01_unseal_evaluation.py` | 解封评价（主结局、§7.2 四项指标、bootstrap） |
| 7 | `p6e_01_lodo_fold_qc.py` | LODO 折构成 QC |
| 8 | `p6f_lodo_settlement.py` | LODO 聚合口径裁定 |
| 9 | `p6g_direction_learning.py` | 方向学习（post-unblinding explanatory） |
| 10 | `p10a_prepare_cache.py` | §10 样本量敏感性：缓存准备 |
| 11 | `p10b_subsample_curves.py` | §10 样本量敏感性：曲线 |
| 12 | `p10c_aggregate.py` | §10 样本量敏感性：汇总 |
| 13 | `p12_configs_corrected.py` | 四配置（feature tier × classifier）对照 |
| 14 | `p11_figure_datapack.py` | 图数据包（图 1–5 / S1–S9 + Table S1/S2） |
| 15 | `p6_02_label_free_signals.py` | **已被 p6b 取代**的早期 signal 流水线；仅作溯源保留 |
| 16 | `p15_results_numeric_audit.py` | Results 数值一致性（203 数值 + 16 定性） |
| 17 | `p16_methods_audit.py` | Methods 数值/语义/完备性（61+14+9） |
| 18 | `p19_discussion_audit.py` | Discussion 断言强度 |
| 19 | `p30_mnv_independent_recompute.py` | MNV：正文 headline numbers 独立重算 |
| 20 | `p31_legends_numeric_audit.py` | 补充图例数值审计 |
| 21 | `p32_direction_convention_sensitivity.py` | 方向约定敏感性 |
| 22 | `p52_mae_vs_constant.py` | MAE 对常数基线的**配对区间**（补 INC-017 缺的区间） |
| 23 | `p53_rounding_boundary_scan.py` | 同一量多生产者·显示精度扫描（INC-019 发现的仪器） |
| 24 | `p54_pointer_value_audit.py` | 指针-值判据（INC-020 发现的仪器） |
| 25 | `rephead.py` | 审计报告的自我声明头（被 p15/p19/p34/p52/p53/p54 共用） |
| 26 | `p55_spelling_audit.py` | 拼写/术语一致性（**全稿**，覆盖 8 文件并显式声明未覆盖者） |
| 27 | `p33_protocol72_metrics.py` | protocol §7.2 四项指标（含 PPV/NPV） |
| 28 | `p34_abstract_audit.py` | Abstract 审计器 |
| 29 | `p35_supplement_inventory_audit.py` | 补充材料产物↔引用双向清点 |
| 30 | `p21a_reporting_profile_derivation.py` | Reporting profile 推导 |
| 31 | `p21b_reporting_profiles.py` | Reporting profile 计算 |
| 32 | `p21c_reporting_profiles_finalize.py` | Reporting profile 定稿 |
| 33 | `p21d_render_reporting_tables.py` | 渲染 Supplementary Tables R1–R3 |
| 34 | `p21e_render_tableS2.py` | 渲染 Supplementary Tables S1 / S2 |
| 35 | `p21f_layerB_supplement_reconcile.py` | 逐 cell 溯源 ＋ 跨路径一致性（含负向测试） |
| 36 | `p45_triage_policy_simulation.py` | 分诊策略模拟（事后；预注册） |
| 37 | `p46_render_tableS3.py` | 渲染 Supplementary Table S3 |
| 38 | `p47_universe_partition.py` | 候选/合格宇宙分区审计 |
| 39 | `p48_status_snapshot.py` | 工作清单状态快照（由文件生成） |
| 40 | `p49_s1_acceptance.py` | Table S1 可复现性验收（7/37 ＋ S2 Panel A 五行） |
| 41 | `p50_dependency_edges.py` | 依赖边登记与代次判定（同代/未覆盖/不错序） |
| 42 | `p50_negtest.py` | p50 负向测试（在副本上跑，输出落盘） |
| 43 | `genstamp.py` | 内部产物代次头 ＋ 覆盖前自动归档 |
| 44 | `p51_package_census.py` | 包级双向普查 ＋ 闭包检查 |
| 45 | `p44_package_code.py` | 代码打包（本脚本自身） |
| 46 | `p42_reorganize_freeze.py` | 冻结层重组织 |
| 47 | `p6_00_make_manifest.py` | 构造 p6 清单 |
| 48 | `p6_01_anchor_rerun.py` | anchor 链重跑 |
| 49 | `p10_preflight_gate.py` | p10 前置门 |
| 50 | `p21g_render_tripod_v3.py` | TRIPOD 渲染 v3 |
| 51 | `run_all_audits.py` | 审计批跑驱动 |
| 52 | `p57_posthoc_baselines.py` | 事后分析的**唯一入库生产者**（source-only 基线 / 配对差 / 方差分解 / 分诊复核） |
| 53 | `p58_reference_order_audit.py` | 参考文献首现序审计（INC-035 发现的仪器；含 R-1..R-4 与三例判据自测） |
| 54 | `indep_verify_p57_auroc_bootstrap.py` | **独立**（自写 Mann–Whitney ＋ 自写重采样）复算 p57 的四个 AUROC 与配对差 |
| 55 | `indep_verify_p57_engine_check.py` | **独立**复算引擎的一致性检查 |
| 56 | `indep_verify_p57_mcnoise_check.py` | **独立**度量蒙特卡洛散布（3 seed × B=10k）＋ B=200k 收敛参考 ＋ 列恒等式 |
| 57 | `p59_survivor_scale_audit.py` | 幸存者代系审查（A4.2 工件）：四个千分位数过**代系 × 聚合健全性**两轴 |
| 58 | `p60_leak_scan.py` | 投稿包**泄漏闸门**（递归扫文件名/文本/docx-xlsx 内部 XML/文档属性/修订痕迹/隐藏表）＋ **全称声明专扫**（只列不判，供通读逐条核对覆盖范围） |
| 59 | `p61_transfer_level_analysis.py` | 转移级 vs 配对层面：**取向注册表** ＋ 已知答案测试 ＋ 带区间的分解 ＋ 不一致对 ＋ 覆盖/可用性 |
| 60 | `p62_seed_stability.py` | 端点**种子稳定性**（5 个非相邻种子）；「不得空口归因于蒙特卡洛」规则的第一个实例 |
| 61 | `p63_seed_fleet.py` | 种子稳定性**总表** —— ★ 状态：**第二遍作废**（INC-037：覆盖了 p6d_01 的种子，连带改了患者级划分）。第一遍的分诊两条有效。**保留在包内是为了留下这道事故的仪器。** |
| 62 | `p21e_render_tableS2.py` | 渲染 Supplementary Table S1 / S2 |
| 63 | `p21f_render_supplementary_methods.py` | 渲染 Supplementary Methods（估计器设置） |
| 64 | `p40_build_submission.py` | 投稿包构建（Markdown → Word） |
| 65 | `p41_verify_submission.py` | 投稿包验收（内容零丢失 / 数字守恒 / 稿件形态） |
| 66 | `p43_package_source_data.py` | Source Data 打包 |
| 67 | `p56_build_submission_folder.py` | 投稿上传文件夹装配（读 p40 产物，含新鲜度与缺件守卫） |

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
