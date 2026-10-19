# -*- coding: utf-8 -*-
r"""p21c — reporting profiles 定稿器

================================ 为什么需要单独一步 ================================
p21a 跑出数字后，`site_profile.csv` 与 `participant_flow.csv` 被 p21a 覆盖，
p21b 版本的 **task-expanded 列丢失**。本定稿器把两条口径**并列装进同一张表**，
并按 GPT6 要求补齐全部元数据列。

**关键原则（GPT6）**：数据已经算出来 ≠ 手稿已经报告出来。
本步只负责把数字落成**正式 reporting artifact**；是否"已报告"由后续 p21 / TRIPOD audit 判。

================================ 两条永久规则（写入每个产物的 metadata）================================
R-1  Episode partition equivalence must be assessed by **partition membership**,
     not raw block identifiers.
R-2  **Task-expanded episode observations must never be substituted for deduplicated
     microbiological episodes or for unique patients.**

用法：python p21c_reporting_profiles_finalize.py
"""
import os, sys, hashlib, io
import numpy as np, pandas as pd

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); FIG = os.path.join(ROOT, "figures", "data")
D = os.path.join(FIG, "reporting_profiles")
VER = "p21c v1.0 (2026-10-02)"

def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for c in iter(lambda: f.read(1 << 22), b""): h.update(c)
    return h.hexdigest()

rep = open(os.path.join(D, "p21c_finalize_report.txt"), "w", encoding="utf-8")
def p(*a):
    line = " ".join(str(x) for x in a); print(line, flush=True); rep.write(line + "\n"); rep.flush()

R = lambda f: pd.read_csv(os.path.join(D, f))
man  = R("p21a_source_manifest.csv")
sp_a = R("site_profile.csv")                    # p21a 版（unique 口径）
flowc= R("participant_flow_counts.csv")
xsd  = R("site_x_organism_drug_counts.csv")
tp   = R("task_profile.csv")                    # p21b 版（task-expanded 口径，未被覆盖）

MAN_SHA = sha(os.path.join(D, "p21a_source_manifest.csv"))
META = dict(
    status="REPORTING-ONLY DERIVATION",
    derivation_script="分析结果/p21a_reporting_profile_derivation.py + p21c_reporting_profiles_finalize.py",
    derivation_version=VER,
    source_manifest_file="reporting_profiles/p21a_source_manifest.csv",
    source_manifest_sha256=MAN_SHA,
    permanent_rule_R1=("Episode partition equivalence must be assessed by partition membership, "
                       "not raw block identifiers."),
    permanent_rule_R2=("Task-expanded episode observations must never be substituted for "
                       "deduplicated microbiological episodes or for unique patients."),
)

# ---------------------------------------------------------------- site_profile 定稿
te = tp.groupby("site").agg(n_task_expanded_episode_observations=("n_episodes_sr", "sum"),
                            n_eligible_transfers_as_source=("n_as_source", "sum"),
                            n_eligible_transfers_as_target=("n_as_target", "sum")).reset_index()
rs = xsd.groupby("site").agg(n_S_task_expanded=("n_S", "sum"), n_R_task_expanded=("n_R", "sum")).reset_index()
fin = sp_a.merge(te, on="site").merge(rs, on="site")
fin["n_organism_drug_groups"] = 38
fin = fin.rename(columns={"n_age_observations": "n_age_observed"})
order = ["site", "n_organism_drug_groups",
         "n_unique_microbiological_episodes", "n_unique_patients",
         "n_task_expanded_episode_observations",
         "n_eligible_transfers_as_source", "n_eligible_transfers_as_target",
         "n_age_observed", "n_age_missing", "pct_age_missing", "age_median", "age_p25", "age_p75",
         "setting_ER", "setting_Inpatient", "setting_Outpatient", "setting_Other", "setting_missing",
         "specimen_URINE", "specimen_BLOOD", "specimen_other",
         "n_S_task_expanded", "n_R_task_expanded",
         "analysis_unit", "deduplication_level", "denominator_definition",
         "source_file", "filter_definition"]
fin = fin[order]
for k, v in META.items(): fin[k] = v
fin.loc[fin.site == fin.site, "unit_of_unique_columns"] = (
    "n_unique_microbiological_episodes / n_unique_patients / age / setting / specimen "
    "为**去重 episode 口径**（episode ID，task 展开之前）；"
    "n_task_expanded_episode_observations / n_S_task_expanded / n_R_task_expanded 为"
    "**未去重、按测药数加权**的 task-expanded 口径。**两组不得混用。**")
genstamp.backup(os.path.join(D, "site_profile.csv"), "p21c")
fin.to_csv(os.path.join(D, "site_profile.csv"), index=False)
import genstamp  # noqa: E402
_genstamp_files = [os.path.join(D, "site_profile.csv")]
p(f"[OK] site_profile.csv 定稿 —— {len(fin)} 行 × {len(fin.columns)} 列")
p(fin[["site", "n_unique_microbiological_episodes", "n_unique_patients",
       "n_task_expanded_episode_observations", "n_age_missing", "pct_age_missing",
       "age_median", "n_S_task_expanded", "n_R_task_expanded"]].to_string(index=False))

# ---------------------------------------------------------------- participant_flow 定稿
rows = []
for _, f in flowc.iterrows():
    s = f["site"]
    raw = int(man[(man.site == s) & (man.role == "microbiology records")].raw_row_count.iat[0])
    rows += [
     dict(site=s, order=1, step="Records in microbiology source file (all organisms)",
          unit="record", value=raw,
          note="**文件级**计数，含非目标菌；MGB 的菌种过滤在 read 时同步完成"),
     dict(site=s, order=2, step="Records for the three focal organisms",
          unit="record", value=int(f.n_rows_focal_organisms),
          note="EC / KP / PA"),
     dict(site=s, order=3, step="Distinct (patient, date, specimen, antibiotic, result) records",
          unit="record", value=int(f.n_rows_dedup_abcd_s), note="p6b drop_duplicates 口径"),
     dict(site=s, order=4, step="Records with an interpretable AST result (S/I/R)",
          unit="record", value=int(f.n_rows_SRI), note=""),
     dict(site=s, order=5, step="Distinct (patient, specimen, date) units entering episode construction",
          unit="record", value=int(f.n_units), note="episode 划分的输入单元"),
     dict(site=s, order=6, step="Microbiological episodes constructed (3 focal organisms)",
          unit="microbiological episode", value=int(f.n_episodes_in_store),
          note="anchor 链，14 天；**task 展开之前**"),
     dict(site=s, order=7, step="Episodes with S/R for ≥1 of the 38 benchmark organism–drug pairs",
          unit="microbiological episode",
          value=int(sp_a[sp_a.site == s].n_unique_microbiological_episodes.iat[0]),
          note="**Tier-1 去重口径**；这是 site_profile 的分母"),
     dict(site=s, order=8, step="Unique patients contributing ≥1 such episode", unit="patient",
          value=int(sp_a[sp_a.site == s].n_unique_patients.iat[0]),
          note="**patient 口径，与 episode 数不同**"),
     dict(site=s, order=9, step="Eligible directed transfers involving this site (as source)",
          unit="directed transfer",
          value=int(tp[tp.site == s].n_as_source.sum()), note="**task 口径**"),
     dict(site=s, order=10, step="Eligible directed transfers involving this site (as target)",
          unit="directed transfer", value=int(tp[tp.site == s].n_as_target.sum()), note="task 口径"),
    ]
pf = pd.DataFrame(rows)
pf["analysis_unit"] = pf.unit
pf["deduplication_level"] = np.select(
    [pf.step.str.contains("Records in"), pf.step.str.contains("focal organisms"),
     pf.step.str.contains("Distinct \\(patient, date"), pf.step.str.contains("interpretable"),
     pf.step.str.contains("Distinct \\(patient, specimen"), pf.step.str.contains("constructed"),
     pf.step.str.contains("≥1 of the 38"), pf.step.str.contains("Unique patients"),
     pf.step.str.contains("directed transfers")],
    ["none（文件级）", "none", "(anon_id, cdate, spec, abx, s)", "none",
     "(anon_id, cdate, spec)", "(anon_id, spec, anchor block)",
     "(anon_id, spec, anchor block)", "patient id", "transfer"],
    default="—")
pf["denominator_definition"] = ("逐节点不同：见 analysis_unit 与 deduplication_level 两列；"
                                "**不得把 episode 数写成 participants**")
pf["source_file"] = "p21a_source_manifest.csv + p6b STORE（p21a 复用重建）"
for k, v in META.items(): pf[k] = v
pf.to_csv(os.path.join(D, "participant_flow.csv"), index=False)
p(f"\n[OK] participant_flow.csv 定稿 —— {len(pf)} 行（{pf.site.nunique()} 站 × {pf['order'].nunique()} 节点）")

# ---------------------------------------------------------------- task_profile：只加 provenance linkage
tp["analysis_unit"] = "organism–drug episode observation (TASK-EXPANDED, 未去重)"
tp["deduplication_level"] = "site × drug（各药独立计数；同一 episode 可进入多个药物任务）"
tp["denominator_definition"] = ("task-expanded episode observations —— "
                                "**与 site_profile 的 unique episode 口径不同，不可互换**")
tp["provenance_link"] = "site_profile.csv / participant_flow.csv / site_x_organism_drug_counts.csv"
tp["source_manifest_sha256"] = MAN_SHA
tp["derivation_script"] = META["derivation_script"]
tp["derivation_version"] = VER
tp["status"] = META["status"]
tp.to_csv(os.path.join(D, "task_profile.csv"), index=False)
p(f"[OK] task_profile.csv —— {len(tp)} 行；**数值未改动**，仅补 provenance linkage 与口径列")

for k, v in META.items():
    xsd[k] = v
xsd["analysis_unit"] = "organism–drug episode observation"
xsd["deduplication_level"] = "site × drug"
xsd.to_csv(os.path.join(D, "site_x_organism_drug_counts.csv"), index=False)

# ---------------------------------------------------------------- manifest
outs = ["site_profile.csv", "participant_flow.csv", "task_profile.csv",
        "site_x_organism_drug_counts.csv", "participant_flow_counts.csv", "p21a_source_manifest.csv"]
mf = pd.DataFrame([dict(artifact=f, sha256=sha(os.path.join(D, f)),
                        rows=len(pd.read_csv(os.path.join(D, f))),
                        status=META["status"], derivation_version=VER) for f in outs])
mf.to_csv(os.path.join(D, "profile_manifest.csv"), index=False)
p("\n[OK] profile_manifest.csv")
p(mf.to_string(index=False))

p("")
p("=== 两条永久规则（已写入每个产物的 metadata 列）===")
p(f"  R-1  {META['permanent_rule_R1']}")
p(f"  R-2  {META['permanent_rule_R2']}")
p("")
p("*** p21c 完成：数字已落成正式 reporting artifact；**尚未判定为『已报告』*** ")
p("    『是否已报告』由 p21 / TRIPOD+AI audit 依据**最终呈现位置**判定。")
rep.close(); print("\ndone ->", D)

try:
    for _p in _genstamp_files:
        genstamp.stamp(_p, __file__.split(chr(92))[-1].split('/')[-1].replace('.py',''))
    print('[genstamp] 代次头已写:', len(_genstamp_files), '个')
except Exception as _e:
    print('[genstamp] 跳过:', _e)
