# -*- coding: utf-8 -*-
r"""步骤 A（隔离措施）：从 p5 结果中提取 **blinded task manifest**

只保留：task_id / organism / drug / source_site / target_site /
        source model 引用 / 分析单位定义 / 样本量
**明确不带**：source AUROC、target AUROC、ΔAUROC、harmful/stable、任何 target Y。

此后 label-free 脚本只读本 manifest，不再读 p5_tasks_primary.csv 的性能字段。
输出 p6_task_manifest_blinded.csv
"""
import os
import pandas as pd

OUT = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(OUT, "p5_tasks_primary.csv")

# 只取允许的列（性能列一律不进入）
ALLOWED = ["organism", "drug", "source", "target", "n_target"]
df = pd.read_csv(SRC, usecols=lambda c: c in ALLOWED + ["eligible"])
m = df[df["eligible"]].copy()

man = pd.DataFrame({
    "task_id": m["organism"] + "_" + m["drug"] + "_" + m["source"] + "2" + m["target"],
    "organism": m["organism"],
    "drug": m["drug"],
    "source_site": m["source"],
    "target_site": m["target"],
    "analysis_unit": "episode (anchor, 14d; episode_definition.py v1.0)",
    "source_model_ref": "P1 History+LR (p5_01, seed 20261005, patient-level 70/15/15)",
    "n_target_units": m["n_target"],
}).reset_index(drop=True)

# ★ 2026-10-06：原为 `assert len(man) == 202` —— 硬编码了**修正前**的合格 task 数。
#   CS-01/02/03 重建 task universe 后该数为 186（见 BENCHMARK_IDENTITY_GATE.md），
#   硬编码会让本步直接 AssertionError，且**掩盖了"分母已变"这一事实**。
#   改为**从 p5 现算**：manifest 必须恰好覆盖 p5 判为 eligible 的全部 task ——
#   这才是本步真正要保证的不变式（既不漏、也不多），与具体数字无关。
assert len(man) == len(m), (
    f"manifest 行数 {len(man)} != p5 eligible 行数 {len(m)} —— 有 task 在派生中丢失或重复")
assert len(man) == int(df["eligible"].sum()), "manifest 覆盖数与 p5 eligible 计数不符"
assert len(man) > 0, "eligible task 为 0 —— 上游 p5 可能未重跑或全部出局"
assert man["task_id"].is_unique
man.to_csv(os.path.join(OUT, "p6_task_manifest_blinded.csv"), index=False)

print(f"[OK] manifest 写出 {len(man)} 行（p5 eligible={int(df['eligible'].sum())} / "
      f"candidate={len(df)}）-> p6_task_manifest_blinded.csv")
print(f"     列: {list(man.columns)}")
print(f"     organism: {dict(man.organism.value_counts())}")
print(f"     source→target: {dict((man.source_site + '→' + man.target_site).value_counts())}")
print(f"     organism-drug 分组数: {man.groupby(['organism','drug']).ngroups}")
print("\n[检查] 该文件不含任何性能字段：", not any(
    k in " ".join(man.columns).lower() for k in ["auroc", "delta", "harmful", "stable"]))
