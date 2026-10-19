# -*- coding: utf-8 -*-
r"""p21b — TRIPOD+AI 报告用数据画像（participant_flow / site_profile / task_profile）

================================ 这是什么，不是什么 ================================
**是**：TRIPOD+AI 报告完整性所需的 **descriptive reporting metadata**。
**不是**：重新分析。本脚本
  · 只读**已冻结**产物；
  · 不改变任何 episode / eligibility / feature / model / signal / outcome；
  · 不进入任何模型管线；
  · 不重算任何已冻结的结果。

================================ 先定 denominator，再算数 ================================
GPT6 裁定（2026-10-02）：**不要先冲着一个数字去算，先把口径写死。** 本脚本对每个量都输出
`denominator_definition` 与 `deduplication_key` 两列，**不允许只给数不给口径**。

三个层级的口径**互不相同，不得互相替代**：
  · participant_flow —— 报告 **unique patients** 与 **microbiological episodes**，两者分开，不与 task 计数混
  · site_profile     —— 分母是 **该站点在 38 个 organism–drug 宇宙内的去重 microbiological episodes**
                        （**在 task 展开之前**去重；同一 episode 不应因测了多药而被重复计权）
  · task_profile     —— 分母是 **organism–drug × site 的 episode 观测**
                        （此处重复是**设计的一部分**，不是错误）

================================ 输出未编号（避免与现有 Supplement 撞号）================================
figures/data/reporting_profiles/
    participant_flow.csv  ·  site_profile.csv  ·  task_profile.csv
用法：python p21b_reporting_profiles.py
"""
import os, sys
import numpy as np, pandas as pd

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); BASE = os.path.dirname(ROOT)
OUT = os.path.join(ROOT, "分析结果"); FIG = os.path.join(ROOT, "figures", "data")
DEST = os.path.join(FIG, "reporting_profiles"); os.makedirs(DEST, exist_ok=True)

PENDING = "PENDING — reporting-only derivation required; not present in frozen outputs"
rep_lines = []
def p(*a):
    line = " ".join(str(x) for x in a); print(line, flush=True); rep_lines.append(line)

B = pd.read_csv(os.path.join(OUT, "p6b_label_free_signals_blinded.csv"))
D = pd.read_csv(os.path.join(OUT, "p6d_outcome.csv"))
B["group"] = B.organism + "|" + B.drug

SRC_B, FLD = "p6b_label_free_signals_blinded.csv", "n_source / n_target / n_src_{train,val,test}"
FILT = ("episodes with a susceptible-or-resistant result for that organism–drug "
        "(p6b prep(): np.isin(X[drug], ['S','R']))")

# ---------------------------------------------------------------- 一致性断言
src = (B[["group", "organism", "drug", "source_site", "n_source", "n_src_train", "n_src_val", "n_src_test"]]
         .rename(columns={"source_site": "site"}))
tgt = (B[["group", "target_site", "n_target"]].rename(columns={"target_site": "site"}))
S = src.groupby(["group", "site"]).n_source.agg(["nunique", "first"]).reset_index()
T = tgt.groupby(["group", "site"]).n_target.agg(["nunique", "first"]).reset_index()
M = S.merge(T, on=["group", "site"], suffixes=("_s", "_t"))
assert (M.nunique_s == 1).all() and (M.nunique_t == 1).all(), "同一 (group,site) 的计数不唯一"
assert ((M.first_s - M.first_t).abs() == 0).all(), "n_source 与 n_target 在同一 (group,site) 上不一致"
p(f"[OK] 一致性断言：{len(M)} 个 (organism–drug, site) 组合，"
  f"n_source（作源径）与 n_target（作靶径）逐项相等，最大差 0")

# ---------------------------------------------------------------- task_profile
tp = (B.groupby(["organism", "drug", "group", "source_site"])
        .agg(n_as_source=("task_id", "size"), n_episodes_sr=("n_source", "first"),
             n_src_train=("n_src_train", "first"), n_src_val=("n_src_val", "first"),
             n_src_test=("n_src_test", "first")).reset_index()
        .rename(columns={"source_site": "site"}))
tp_t = (B.groupby(["group", "target_site"])
          .agg(n_as_target=("task_id", "size"), n_episodes_from_target=("n_target", "first"))
          .reset_index().rename(columns={"target_site": "site"}))
tp = tp.merge(tp_t, on=["group", "site"], how="outer")

# 【修正】202 个 eligible transfer 中有 13 个 (group, site) **只作靶、不作源**，
# 故 n_source 取不到值；但这些行的 n_target 存在，且已断言二者必然相等。
# 首版留下 13 行空白 —— 空白会被误读为"该站该组无数据"，而事实是**有数据、只是未从靶侧取**。
_both = tp.n_episodes_sr.notna() & tp.n_episodes_from_target.notna()
assert (tp.loc[_both, "n_episodes_sr"] == tp.loc[_both, "n_episodes_from_target"]).all(), \
    "双侧都有值时 n_source 与 n_target 必须相等"
_fill = tp.n_episodes_sr.isna()
tp.loc[_fill, "n_episodes_sr"] = tp.loc[_fill, "n_episodes_from_target"]
assert tp.n_episodes_sr.notna().all(), "补齐后不应再有缺失"
p(f"[OK] 从靶侧补齐 {int(_fill.sum())} 个仅在 eligible 宇宙中作靶的 (group, site) 计数")
tp["episode_count_source_role"] = np.where(_fill, "target-role only", "source-role")
tp["analysis_level"] = "organism–drug × site"
tp["denominator_definition"] = ("S/R-eligible microbiological episodes at this site for this "
                                "organism–drug; task expansion duplicates an episode across "
                                "its directed transfers BY DESIGN")
tp["deduplication_key"] = "one row per (organism, drug, site); episodes deduplicated within site×drug"
tp["source_file"] = SRC_B; tp["source_field"] = FLD; tp["filter_definition"] = FILT
tp["n_as_source"] = tp.n_as_source.fillna(0).astype(int)
tp["n_as_target"] = tp.n_as_target.fillna(0).astype(int)
tp = tp[["organism", "drug", "group", "site", "analysis_level", "n_episodes_sr",
         "episode_count_source_role", "n_src_train", "n_src_val", "n_src_test",
         "n_as_source", "n_as_target",
         "denominator_definition", "deduplication_key", "source_file", "source_field",
         "filter_definition"]].sort_values(["group", "site"])
tp.to_csv(os.path.join(DEST, "task_profile.csv"), index=False)
_tot = int(tp.n_episodes_sr.sum())
p(f"[OK] task_profile.csv — {len(tp)} 行（= 38 organism–drug × 3 sites，完整网格）；"
  f"其中 {int((tp.episode_count_source_role=='target-role only').sum())} 行计数取自靶侧")
p(f"     n_episodes_sr 合计 = {_tot:,} —— 口径是 **task-expanded episode observations**"
  f"（organism–drug × site 展开后）；同一 episode 因多药重复计入，"
  f"**不是** unique microbiological episodes，**不得**简称为 episodes")

# ---------------------------------------------------------------- site_profile
rows = []
for site in ["MGB", "SU", "UT"]:
    s = tp[tp.site == site]
    rows.append(dict(
        site=site,
        n_organism_drug_groups=int(s.group.nunique()),
        n_task_expanded_episode_observations=int(s.n_episodes_sr.sum()),
        n_unique_microbiological_episodes=PENDING,
        n_unique_patients=PENDING,
        median_episodes_per_group=float(s.n_episodes_sr.median()),
        min_episodes_per_group=int(s.n_episodes_sr.min()),
        max_episodes_per_group=int(s.n_episodes_sr.max()),
        n_eligible_transfers_as_source=int(s.n_as_source.sum()),
        n_eligible_transfers_as_target=int(s.n_as_target.sum()),
        age_distribution=PENDING, care_setting_distribution=PENDING,
        specimen_distribution=PENDING, outcome_S_R_counts=PENDING,
        missingness_by_field=PENDING,
        analysis_unit="microbiological episode (for the PENDING unique count)",
        deduplication_level="episode ID, BEFORE task expansion",
        denominator_definition=("n_unique_microbiological_episodes 的分母 = 该站点在 38 个 "
                                "organism–drug 宇宙内的**去重** microbiological episodes，"
                                "在 task 展开之前按 episode ID 去重；"
                                "**不可**用 n_task_expanded_episode_observations 代替"),
        caveat=("**n_task_expanded_episode_observations 不是『该站点的微生物学 episode 数』。** "
                "它是 organism–drug × site 展开后的 **task-expanded episode observations**："
                "同一 episode 若对多个药物有结果，会进入多个药物任务而被重复计数。"
                "**不得**简称为 episodes，**不得**用于 site_profile 的 unique episodes，"
                "**更不得**当作 patient/episode flow 的人数。"),
        deduplication_key=("去重键应为 (patient, specimen, episode anchor date)；"
                           "**当前冻结产物中无 episode 级表，故无法计算**"),
        source_file=SRC_B, source_field=FLD, filter_definition=FILT))
sp = pd.DataFrame(rows)
# ★ 2026-10-08（A3.2）：**暂存写**。真名 `site_profile.csv` 只由它的**真生产者**（p21a，读原始数据）
#   写；本步的产物是**降级版**（episodes/patients/age/… 均为 PENDING），
#   若直接写真名，会把 p21a 的完整产物**覆盖掉**，并让后续 p21c 报 KeyError。
#   病史：2026-10-07 实测踩中 —— 单跑 p21b 后 site_profile.csv 被换成 23 列的降级版，
#   p21c 随即崩（缺 17 列），我从存档恢复并补跑 p21a 才复原（INC-002 的同族）。
#   改为暂存名后，**p21a 与 p21b 的执行顺序不再重要**（约束是被制造出来的，此处解除）。
_pend = os.path.join(DEST, "site_profile.pending.csv")
genstamp.backup(_pend, "p21b")
sp.to_csv(_pend, index=False)
import genstamp  # noqa: E402
_genstamp_files = [_pend]
p(f"[OK] site_profile.**pending**.csv — {len(sp)} 行；分子可得、"
  f"**unique episodes / patients / 特征分布 / missingness 均为 {PENDING.split(' — ')[0]}**")
p("     ⚠ 本文件是**降级版暂存**；真名 `site_profile.csv` 由 p21a 写（需原始数据）。"
  " 二者**不得互换**，p21c 只读真名。")

# ---------------------------------------------------------------- participant_flow
pf = pd.DataFrame([
 dict(step="1. Microbiology records available at site", analysis_unit="record",
      value=PENDING, source_file="—", note="冻结产物中无原始记录计数"),
 dict(step="2. Records collapsed to microbiological episodes", analysis_unit="episode",
      value=PENDING, source_file="—",
      note="episode_definition.py 可复算，但**其输出未作为冻结产物保存**"),
 dict(step="3. Episodes with a susceptible-or-resistant result", analysis_unit="episode",
      value=PENDING, source_file="—", note="同上"),
 dict(step="4. Episodes within the 38 organism–drug benchmark universe",
      analysis_unit="episode", value=PENDING, source_file="—",
      note="需按 (organism, drug, site) 展开后去重"),
 dict(step="5. Unique patients contributing ≥1 such episode", analysis_unit="patient",
      value=PENDING, source_file="—", note="冻结产物中无患者标识"),
 dict(step="6. Eligible directed transfers (task level)", analysis_unit="directed transfer",
      value=202, source_file="fig1b_universe.csv (eligible_tasks)",
      note="**这是 task 计数，不是 participant 计数；不得与上列混用**"),
 dict(step="7. Organism–drug groups", analysis_unit="group", value=38,
      source_file="fig1b_universe.csv (n_organism_drug_groups)", note=""),
 dict(step="8. Candidate transfers before eligibility", analysis_unit="directed transfer",
      value=228, source_file="fig1b_universe.csv (candidate_tasks)", note=""),
])
pf["denominator_definition"] = ("participant_flow 必须**分开**报告 unique patients 与 "
                                "microbiological episodes；不得把 episodes 称作 participants")
pf["deduplication_key"] = "patients 去重键 = patient id；episodes 去重键 = (patient, specimen, anchor date)"
pf.to_csv(os.path.join(DEST, "participant_flow.csv"), index=False)
p(f"[OK] participant_flow.csv — 8 步中 **5 步为 PENDING**；仅第 6–8 步（task 计数）可由冻结产物支撑")

# ---------------------------------------------------------------- 报告
p("")
p("=== 口径声明（三层互不替代）===")
p(f"  site_profile     分母 = 该站点在 38 组宇宙内的**去重** microbiological episodes（task 展开前去重）")
p(f"  task_profile     分母 = organism–drug × site 的 episode 观测（重复是设计的一部分）")
p(f"  participant_flow 分开报告 unique patients 与 episodes，不与 task 计数混用")
p("")
p("=== 冻结产物中【不存在】的量（故标 PENDING，**不填 NA**）===")
for k in ["patient identifier", "episode-level table", "raw record counts",
          "age / care setting / specimen distributions", "R/S counts per task",
          "missingness", "age missingness rate (`age_miss` computed in p6b but never persisted)"]:
    p(f"  · {k}")
p("")
p("=== 下一步（p21a_reporting_profile_missingness.py 的输入需求）===")
p("  先固定 missingness 要报告哪些字段，再定 cohort。候选字段与其**设计上的**缺失可能性：")
for f, note in [("age", "**唯一**可缺失；以 source-train 中位数填补"),
                ("care setting", "由 ward flag 派生，设计上不缺失"),
                ("specimen", "由 culture_description 归一化派生，设计上不缺失"),
                ("history-derived predictors", "设计上编码为 none/全零，**非缺失**"),
                ("AST outcome", "进入 S/R cohort 后**理论上 0 缺失**（限定所致）"),
                ("other categorical", "需确认是否已在上游 harmonization 中剔除/归类")]:
    p(f"    - {f:<28} {note}")
p("  若最终确认只需派生 age，再运行 p21a；**不预先单独算 age_miss**。")

open(os.path.join(DEST, "reporting_profiles_report.txt"), "w", encoding="utf-8").write("\n".join(rep_lines))
print("\ndone ->", DEST)

try:
    for _p in _genstamp_files:
        genstamp.stamp(_p, __file__.split(chr(92))[-1].split('/')[-1].replace('.py',''))
    print('[genstamp] 代次头已写:', len(_genstamp_files), '个')
except Exception as _e:
    print('[genstamp] 跳过:', _e)
