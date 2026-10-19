# -*- coding: utf-8 -*-
r"""p21f — 渲染 Supplementary Methods：四个预设配置的估计器设置

================================ 为什么需要 ================================
Methods §5 末句：*"Exact estimator settings for every configuration are listed in
**Supplementary Methods**."* —— 而**该产物不存在**。

这是 p35 双向清点抓到的**悬空引用**（正文引用了不存在的产物），比「孤儿表」更硬：
审稿人按图索骥会扑空。

================================ 值从哪来 ================================
**不手工录入**。全部从权威实现 `p12_configs_corrected.py` 的源码常量**程序化抽取**，
每一条都 fail-closed（抽不到即抛错，不回退到默认值）。

配置名与分层另与**冻结产物** `figures/data/fig1e_configs.csv` 交叉核对 ——
防止脚本常量与已冻结结果脱节。

输出：figures/data/reporting_profiles/Supplementary_Methods.md
用法：python p21f_render_supplementary_methods.py
"""
import io, os, re, sys
import pandas as pd

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIG  = os.path.join(ROOT, "figures", "data")
D    = os.path.join(FIG, "reporting_profiles")
SRC  = os.path.join(ROOT, "分析结果", "p12_configs_corrected.py")
VER  = "p21f v1.0 (2026-10-07)"
os.makedirs(D, exist_ok=True)

SRC_TXT = io.open(SRC, encoding="utf-8").read()
ROWS = []
def rec(cid, item, ok, note=""):
    ROWS.append(dict(cid=cid, item=item, status="PASS" if ok else "FAIL", note=note))
    print("  [%-4s] %-60s %s" % ("PASS" if ok else "FAIL", item[:60], note), flush=True)
    return ok

def grab(pat, what, cast=str):
    """fail-closed 抽取：抽不到直接抛，不回退默认值。"""
    m = re.search(pat, SRC_TXT)
    if not m:
        raise RuntimeError("抽取失败（%s）：pattern=%r —— 权威实现可能已改动，须人工核对后更新 pattern" % (what, pat))
    return cast(m.group(1))

print("=" * 100); print("p21f — Supplementary Methods 渲染"); print("=" * 100)
print(f"权威实现：{SRC}")
print("")

# ★ 不加 ^ 锚点：p12 源码把这两个常量与 t0 写在同一行
#   （`t0 = time.time(); SEED = 20261005; AUROC_MIN = 0.65`），首版用 ^SEED 直接抽不到。
SEED      = grab(r"\bSEED\s*=\s*(\d+)", "SEED", int)
AUROC_MIN = grab(r"\bAUROC_MIN\s*=\s*([\d.]+)", "AUROC_MIN", float)
CONFIGS   = re.findall(r'\("([^"]+)",\s*"(hist|core)",\s*"(lr|lgbm)"\)', SRC_TXT)
LR_MAXIT  = grab(r"LogisticRegression\(max_iter=(\d+),\s*C=([\d.]+),\s*random_state=SEED\)", "LR max_iter", int)
LR_C      = grab(r"LogisticRegression\(max_iter=\d+,\s*C=([\d.]+),\s*random_state=SEED\)", "LR C", float)
GB_N      = grab(r"LGBMClassifier\(n_estimators=(\d+)", "LGBM n_estimators", int)
GB_L      = grab(r"num_leaves=(\d+)", "LGBM num_leaves", int)
GB_LR     = grab(r"learning_rate=([\d.]+)", "LGBM learning_rate", float)
GB_J      = grab(r"n_jobs=(\d+)", "LGBM n_jobs", int)
CAL_MAXIT = grab(r"cal = LogisticRegression\(max_iter=(\d+)\)", "Platt max_iter", int)
SPLIT     = re.findall(r"ntr, nva = int\(([\d.]+) \* len\(upat\)\), int\(([\d.]+) \* len\(upat\)\)", SRC_TXT)
MIN_EP    = grab(r"if len\(Xs\) < (\d+): continue", "min episodes", int)
MIN_SPLIT = grab(r"if min\(len\(i_tr\), len\(i_va\), len\(i_te\)\) < (\d+): continue", "min per split", int)

print(f"  抽取：SEED={SEED}  AUROC_MIN={AUROC_MIN}  CONFIGS={len(CONFIGS)}  "
      f"LR(C={LR_C}, max_iter={LR_MAXIT})  LGBM(n={GB_N}, leaves={GB_L}, lr={GB_LR}, jobs={GB_J})  "
      f"Platt(max_iter={CAL_MAXIT})  阈值 {MIN_EP}/{MIN_SPLIT}")
print("")

rec("E1", "SEED 已抽取", isinstance(SEED, int))
rec("E2", "四个配置已抽取（tier × algo 全组合）", len(CONFIGS) == 4,
    "｜".join("%s(%s,%s)" % c for c in CONFIGS))
rec("E3", "tier 取值集合 == {hist, core}", {c[1] for c in CONFIGS} == {"hist", "core"})
rec("E4", "algo 取值集合 == {lr, lgbm}", {c[2] for c in CONFIGS} == {"lr", "lgbm"})
rec("E5", "划分比例已抽取且合计 0.85（余下 0.15 为 internal test）",
    len(SPLIT) == 1 and abs(float(SPLIT[0][0]) + float(SPLIT[0][1]) - 0.85) < 1e-9,
    f"{SPLIT[0][0]}/{SPLIT[0][1]}" if SPLIT else "抽不到")

# 与冻结产物交叉核对
F1E = pd.read_csv(os.path.join(FIG, "fig1e_configs.csv"))
frozen_cfgs = set(F1E.config)
src_cfgs = {c[0] for c in CONFIGS}
rec("E6", "配置名集合 == 冻结产物 fig1e_configs.csv 的 config 列", src_cfgs == frozen_cfgs,
    f"脚本 {sorted(src_cfgs)} ｜ 冻结 {sorted(frozen_cfgs)}")
tier_map = dict(zip(F1E.config, F1E.tier)); algo_map = dict(zip(F1E.config, F1E.algo))
mism = [c[0] for c in CONFIGS if tier_map.get(c[0]) != c[1] or algo_map.get(c[0]) != c[2]]
rec("E7", "每个配置的 tier/algo 与冻结产物一致", not mism, f"不符={mism or '无'}")

# ---------------------------------------------------------------- 渲染
TIER = {"hist": "History-augmented (prior resistance history included)",
        "core": "Core (no prior-history features)"}
ALGO = {"lr": "L2-penalised logistic regression",
        "lgbm": "Gradient-boosted trees (LightGBM)"}
md = [
 "# Supplementary Methods — Estimator settings for the four prespecified configurations",
 "",
 "> All four configurations share one source-model development procedure, described in Methods "
 "§5. The tables below give the exact setting of every parameter for each configuration. "
 "Estimator settings are stated in the notation of the implementations used "
 "(scikit-learn for logistic regression, LightGBM for gradient-boosted trees).",
 "",
 "## Shared across all four configurations",
 "",
 "| Item | Setting |",
 "|---|---|",
 f"| Split unit | Patient (no patient spans two splits) |",
 f"| Split proportions | {float(SPLIT[0][0]):.0%} train / {float(SPLIT[0][1]):.0%} validation / "
 f"{1 - float(SPLIT[0][0]) - float(SPLIT[0][1]):.0%} internal test |",
 f"| Minimum eligible episodes (source) | {MIN_EP} |",
 f"| Minimum per split | {MIN_SPLIT} |",
 f"| Missing age | Imputed with the **source training-set median** |",
 f"| Continuous features | `StandardScaler` fitted on the **source training set only** |",
 f"| Random seed | `{SEED}` (shared by splitting, model fitting and every downstream step) |",
 f"| Platt calibration | Logistic regression (`max_iter={CAL_MAXIT}`) fitted on the "
 f"**source validation** split only; the map is frozen and applied unchanged |",
 f"| Source-performance screen | source internal-test AUROC ≥ `{AUROC_MIN}` |",
 "",
 "## Per-configuration settings",
 "",
 "| Configuration | Feature tier | Classifier | Classifier hyperparameters |",
 "|---|---|---|---|"]
for name, tier, algo in sorted(CONFIGS, key=lambda c: c[0]):
    if algo == "lr":
        hp = f"`max_iter={LR_MAXIT}`, `C={LR_C}` (L2), `random_state={SEED}`"
    else:
        hp = (f"`n_estimators={GB_N}`, `num_leaves={GB_L}`, `learning_rate={GB_LR}`, "
              f"`random_state={SEED}`, `n_jobs={GB_J}`")
    md.append(f"| `{name}` | {TIER[tier]} | {ALGO[algo]} | {hp} |")
md += ["",
       "> The two feature tiers differ **only** in whether prior resistance history is included: the "
       "history-augmented tier adds `prior_any`, `prior_cnt`, four prior-history time-bin indicators "
       "and three prior-result indicators for the same organism–drug pair. Every other feature is "
       "identical between tiers.",
       "",
       "> These settings are the values **actually used**; they were fixed in advance and not tuned "
       "against target performance (Methods §5).",
       ""]

# ★★ 2026-10-14（外部核验决定）：**把 Methods §15 的实现细节搬到这里**。
#   主文 §15 只留「两个缺陷」的披露（它们是对冻结方案的**偏离**，必须在主文可见）；
#   其余（九张表逐字节复现、SHA-256 侧车、自动一致性脚本、出处错误的修复史）
#   **原文照搬，不改写、不删** —— 它们是真披露，只是不属于主文粒度。
#   ★ 搬进来之后，主文那句指向 "Reproducibility and audit safeguards" 的**指针必须落在这个标题上**
#     —— 本仓的指针落空已出过两次（INC-032 族）。
md += [
    "",
    "## Reproducibility and audit safeguards",
    "",
    "**Further safeguards.** A dated analysis audit log records every analysis decision, defect and "
    "correction. Further figure-data provenance errors — including use of a superseded configuration "
    "source, an incorrectly implemented leave-one-organism-out panel, and superseded calibration "
    "binning — were identified by provenance auditing, traced to their code paths, and rebuilt. Each "
    "was a recomputation of quantities already reported and changed no signal definition, threshold or "
    "endpoint. A further pass audited the provenance of the derived summary tables themselves. Nine "
    "such tables were found to have been produced interactively rather than by a saved script; each was "
    "subsequently regenerated from its own frozen upstream artifacts by a dedicated reconstruction "
    "script, and all nine reproduced their predecessors **byte for byte**.",
    "",
    "Every data file underlying a figure — including the derived summary tables — is traced to a script "
    "that regenerates it from frozen upstream artifacts, and each regenerated file was verified to be "
    "byte-identical to its corresponding pre-existing frozen sidecar.",
    "",
    "The manuscript's reported numbers are checked against the frozen outputs by an automated "
    "consistency script that reads the frozen files without recomputing any statistic, covering every "
    "reported quantity together with qualitative assertions about analysis identity. Several of those "
    "assertions are made directly against the frozen outputs — for example, that the grouped evaluation "
    "in Figure 3A is identical, method by method, to the logistic-learned orientation and not to the "
    "native one — so that a later change in analysis semantics cannot pass unnoticed.",
    "",
    "The code-level blinding gate described in the main text makes the pre-unblinding status of the "
    "quality-control analysis verifiable rather than declarative.",
    "",
    "**Behaviour on incomplete or insufficient deployment input.** The pipeline has two defined "
    "responses, both of which were exercised rather than assumed. A target episode with a missing "
    "covariate is **not** dropped: a source-derived fill value is applied and the prediction proceeds, "
    "so a signal is still produced. A transfer whose target cohort is too small to support a stable "
    "estimate — fewer than 300 eligible target episodes — does **not** silently yield a signal; the "
    "transfer is recorded as failed with its failure reason and excluded from reporting rather than "
    "being carried forward as a value. In the benchmark as executed, all 186 eligible transfers "
    "returned a valid signal and none failed on this criterion.",
    "",
]
io.open(os.path.join(D, "Supplementary_Methods.md"), "w", encoding="utf-8").write("\n".join(md))

# 渲染件回读：表内出现的每个配置名必须来自抽取结果
rnd = io.open(os.path.join(D, "Supplementary_Methods.md"), encoding="utf-8").read()
rec("R1", "渲染件包含全部四个配置行", all(f"`{c[0]}`" in rnd for c in CONFIGS))
rec("R2", "渲染件含共享设置表与逐配置表", "## Shared across all four configurations" in rnd
    and "## Per-configuration settings" in rnd)
rec("R3", "渲染件中的 SEED 与抽取值一致", f"`{SEED}`" in rnd)

print("")
Dd = pd.DataFrame(ROWS); nbad = int((Dd.status != "PASS").sum())
print("=" * 100)
print(f"p21f 断言合计 {len(Dd)} 项：PASS {int((Dd.status=='PASS').sum())}｜**非 PASS {nbad}**")
print("*** p21f 通过 ***" if not nbad else "*** p21f 未通过 ***")
Dd.to_csv(os.path.join(D, "supplementary_methods_render_audit.csv"), index=False, encoding="utf-8-sig")
sys.exit(1 if nbad else 0)
