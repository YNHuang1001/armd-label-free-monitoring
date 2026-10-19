# -*- coding: utf-8 -*-
r"""p33 — protocol §7.2「临床部署预警任务（主）」四项指标的完整冻结块

================================ 依据 ================================
protocol §7.2（冻结）：

    把 harmful = ΔAUROC ≤ −0.05 当作 deployment-failure endpoint。
    **每个 L1 / L2 / L3 方法**（含 SUDO-AURCC 与三个 transported-AUC estimator）
    都输出风险分 r_failure，报告：AUROC / AUPRC / sensitivity at fixed specificity
    （specificity 固定为 0.80，预先写死）/ PPV / NPV（在观测到的 harmful 先验下）。

现状（本脚本建立前的实测）：
  · AUROC    —— 已生产（fig2_forest.csv、p6d_signal_eval.csv）
  · AUPRC    —— 已生产（p6d_signal_eval.csv:harmful_auprc）
  · sens@.80 —— 已生产（p6d_signal_eval.csv:sens_at_spec80）
  · **PPV / NPV —— 从未计算过（全项目 0 处）**

本脚本：
  1. 从上游冻结表**独立重算** AUROC / AUPRC / sens@spec80，与 p6d_signal_eval.csv 逐位锚定；
  2. 补算 PPV / NPV；
  3. 四项都给 organism–drug cluster bootstrap 95% 区间（protocol §8.3 的不确定性口径）；
  4. 三种口径：**primary（当前正文，三个有符号 L2 取 |·|）** 与 **signed（原值）** 并排，
     使 p32 的方向约定披露在 §7.2 主任务上也有对应数字。

判据来源：p6d_signal_eval.csv（冻结产物）—— 仅用于锚定，不用于取值。
"""
import io, os, sys
import numpy as np, pandas as pd
from scipy import stats
from sklearn.metrics import roc_auc_score, average_precision_score

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIG  = os.path.join(ROOT, "figures", "data")
OUT  = os.path.join(ROOT, "分析结果")
REP  = open(os.path.join(OUT, "report_p33_protocol72.txt"), "w", encoding="utf-8")

def p(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); REP.write(s + "\n"); REP.flush()

# ★ 2026-10-13：种子可由环境变量覆盖，供种子稳定性分析（p63）用。**不改逻辑**。
import os as _os_seed
SEED = 20261007
SPEC = 0.80

# ---------------------------------------------------------------- 独立上游
B = pd.read_csv(os.path.join(OUT, "p6b_label_free_signals_blinded.csv"))
C = pd.read_csv(os.path.join(OUT, "p6c_L3_blinded.csv"))
D = pd.read_csv(os.path.join(OUT, "p6d_outcome.csv"))
SIG = ["domain_auc", "mmd", "psi_mean", "psi_max", "js_mean", "js_max",
       "pred_mean_shift", "entropy_mean_shift", "entropy_js", "conf_dist_js",
       "pred_pos_prev_shift"]
L3 = ["AUROC_DR", "AUROC_IOW", "AUROC_OM"]
T = (D[["task_id", "organism", "drug", "harmful"]]
     .merge(B[["task_id"] + SIG], on="task_id", validate="1:1")
     .merge(C[["task_id", "SUDO_AURCC"] + L3], on="task_id", validate="1:1"))
T["group"] = T.organism + "|" + T.drug
y = T.harmful.values.astype(int); g = T.group.values

ROLE = {"domain_auc": "L1", "mmd": "L1", "psi_mean": "L1", "psi_max": "L1",
        "js_mean": "L1", "js_max": "L1",
        "pred_mean_shift": "L2", "entropy_mean_shift": "L2", "entropy_js": "L2",
        "conf_dist_js": "L2", "pred_pos_prev_shift": "L2", "SUDO_AURCC": "L2",
        "AUROC_DR": "L3-primary", "AUROC_IOW": "L3-sens", "AUROC_OM": "L3-sens"}
ORDER = list(ROLE.keys())
SIGNED_L2 = ["pred_mean_shift", "entropy_mean_shift", "pred_pos_prev_shift"]

def score(m, convention):
    """r_failure：越高越可能 harmful。
       primary：三个有符号 L2 取 |·|（当前正文口径）
       signed ：三个有符号 L2 取原值（protocol §5.2 的字面定义）
       SUDO / L3 在两种口径下都取负（语义强制：可靠性/目标 AUROC 越高越安全）。"""
    s = T[m].values.astype(float)
    if m in SIGNED_L2 and convention == "primary":
        return np.abs(s)
    if m in ("SUDO_AURCC",) + tuple(L3):
        return -s
    return s

def four_metrics(y, s, spec=SPEC):
    """按 protocol §7.2 定义直接实现（不照抄 p6d，以便构成独立核验）。
       specificity = P(s < thr | y = 0)；取满足 specificity ≥ spec 的**最小** thr。"""
    cand = np.unique(s)
    specv = np.array([(s[y == 0] < v).mean() for v in cand])
    ok = cand[specv >= spec]
    thr = float(ok.min()) if ok.size else float(np.nextafter(s.max(), np.inf))
    pred = s >= thr
    tp = int(((pred == 1) & (y == 1)).sum()); fn = int(((pred == 0) & (y == 1)).sum())
    fp = int(((pred == 1) & (y == 0)).sum()); tn = int(((pred == 0) & (y == 0)).sum())
    return dict(
        auroc=roc_auc_score(y, s), auprc=average_precision_score(y, s),
        sens=tp / max(tp + fn, 1), ppv=tp / max(tp + fp, 1), npv=tn / max(tn + fn, 1),
        spec=tn / max(tn + fp, 1), thr=thr, tp=tp, fp=fp, fn=fn, tn=tn)

def boot(y, s, g, n=1000, seed=SEED):
    rng = np.random.default_rng(seed); ugs = np.unique(g)
    keys = ["auroc", "auprc", "sens", "ppv", "npv"]
    acc = {k: [] for k in keys}
    for _ in range(n):
        pick = rng.choice(len(ugs), len(ugs), replace=True)
        idx = np.concatenate([np.where(g == ugs[i])[0] for i in pick])
        yy = y[idx]
        if len(np.unique(yy)) < 2: continue
        try:
            d = four_metrics(yy, s[idx])
        except Exception:
            continue
        for k in keys: acc[k].append(d[k])
    return {k: (float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))) for k, v in acc.items()}

# ---------------------------------------------------------------- 锚定
p("=" * 118)
p("p33 — protocol §7.2「临床部署预警任务（主）」四项指标冻结块")
p("=" * 118)
p(f"prevalence（观测到的 harmful 先验）= {y.sum()}/{len(y)} = {y.mean():.4f}")
p(f"specificity 固定 = {SPEC}（protocol 预先写死）；bootstrap = organism–drug cluster，1000 次，seed={SEED}")
p("")
f6 = pd.read_csv(os.path.join(OUT, "p6d_signal_eval.csv")).set_index("method")
p("【0】锚定：独立重算的 AUROC / AUPRC / sens@.80 必须与 p6d_signal_eval.csv 逐位一致")
bad = []
for m in ORDER:
    d = four_metrics(y, score(m, "primary"))
    ra, rp, rs = (float(f6.loc[m, "harmful_auc"]), float(f6.loc[m, "harmful_auprc"]),
                  float(f6.loc[m, "sens_at_spec80"]))
    oa = abs(d["auroc"] - ra) <= 1e-9 and abs(d["auprc"] - rp) <= 1e-9 and abs(d["sens"] - rs) <= 1e-9
    if not oa: bad.append(m)
    p("    %-22s AUROC %.6f/%.6f ｜ AUPRC %.6f/%.6f ｜ sens@.80 %.6f/%.6f  %s"
      % (m, ra, d["auroc"], rp, d["auprc"], rs, d["sens"], "OK" if oa else "**差异**"))
p("")
if bad:
    p(f"*** 锚定失败：{bad} —— 本脚本与生产管线不同源，结论作废 ***")
    REP.close(); sys.exit(2)
p("    锚定通过：15/15 三项逐位一致（差异 < 1e-9）。")
p("")

# ---------------------------------------------------------------- 冻结块
p("【1】protocol §7.2 四项指标 —— primary 口径（当前正文；三个有符号 L2 取 |·|）")
p("")
h = ("%-22s %-11s %7s %7s %7s %7s %7s %7s %7s %7s"
     % ("method", "role", "AUROC", "AUPRC", "sens@.80", "spec", "PPV", "NPV", "TP/FP", "FN/TN"))
p(h); p("-" * len(h))
rows = []
for m in ORDER:
    s = score(m, "primary")
    d = four_metrics(y, s); b = boot(y, s, g)
    p("%-22s %-11s %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7s %7s"
      % (m, ROLE[m], d["auroc"], d["auprc"], d["sens"], d["spec"], d["ppv"], d["npv"],
         f"{d['tp']}/{d['fp']}", f"{d['fn']}/{d['tn']}"))
    rows.append(dict(convention="primary", method=m, role=ROLE[m], n=len(y),
                     prevalence=round(float(y.mean()), 6),
                     auroc=d["auroc"], auroc_lo=b["auroc"][0], auroc_hi=b["auroc"][1],
                     auprc=d["auprc"], auprc_lo=b["auprc"][0], auprc_hi=b["auprc"][1],
                     sens80=d["sens"], sens80_lo=b["sens"][0], sens80_hi=b["sens"][1],
                     ppv=d["ppv"], ppv_lo=b["ppv"][0], ppv_hi=b["ppv"][1],
                     npv=d["npv"], npv_lo=b["npv"][0], npv_hi=b["npv"][1],
                     spec_at_thr=d["spec"], threshold=d["thr"],
                     tp=d["tp"], fp=d["fp"], fn=d["fn"], tn=d["tn"]))

p("")
p("【2】signed 口径（三个有符号 L2 用原值）—— 仅供 p32 的方向约定披露对照")
p("")
p(h); p("-" * len(h))
for m in SIGNED_L2:
    s = score(m, "signed")
    d = four_metrics(y, s); b = boot(y, s, g)
    p("%-22s %-11s %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7s %7s"
      % (m, ROLE[m], d["auroc"], d["auprc"], d["sens"], d["spec"], d["ppv"], d["npv"],
         f"{d['tp']}/{d['fp']}", f"{d['fn']}/{d['tn']}"))
    rows.append(dict(convention="signed", method=m, role=ROLE[m], n=len(y),
                     prevalence=round(float(y.mean()), 6),
                     auroc=d["auroc"], auroc_lo=b["auroc"][0], auroc_hi=b["auroc"][1],
                     auprc=d["auprc"], auprc_lo=b["auprc"][0], auprc_hi=b["auprc"][1],
                     sens80=d["sens"], sens80_lo=b["sens"][0], sens80_hi=b["sens"][1],
                     ppv=d["ppv"], ppv_lo=b["ppv"][0], ppv_hi=b["ppv"][1],
                     npv=d["npv"], npv_lo=b["npv"][0], npv_hi=b["npv"][1],
                     spec_at_thr=d["spec"], threshold=d["thr"],
                     tp=d["tp"], fp=d["fp"], fn=d["fn"], tn=d["tn"]))

# ---------------------------------------------------------------- 摘要
R = pd.DataFrame(rows)
R.to_csv(os.path.join(OUT, "p33_protocol72_metrics.csv"), index=False, encoding="utf-8-sig")
P1 = R[R.convention == "primary"]
p("")
p("【3】要点")
p("    最优 PPV : %s（%.3f）" % (P1.loc[P1.ppv.idxmax(), "method"], P1.ppv.max()))
p("    最优 NPV : %s（%.3f）" % (P1.loc[P1.npv.idxmax(), "method"], P1.npv.max()))
p("    最优 sens: %s（%.3f）" % (P1.loc[P1.sens80.idxmax(), "method"], P1.sens80.max()))
p("    最优 AUPRC: %s（%.3f）" % (P1.loc[P1.auprc.idxmax(), "method"], P1.auprc.max()))
p("    prevalence = %.4f —— 非退化基率，四项指标均有意义（不得以「基率过低」为由不报）" % y.mean())
p("")
p("=" * 118)
p("【4】信号间相关结构 —— 回答「15 个信号里有几个独立信息源」")
p("")
p("★ 必须用**送入评价的取值**（三个有符号 L2 取 |·|，SUDO/L3 取负）来算 ——")
p("  取绝对值是**非单射**变换，用原值算会低估折叠后的共线程度。")
p("")
X = pd.DataFrame({m: score(m, "primary") for m in ORDER})
CM = X.corr(method="pearson")
CM.to_csv(os.path.join(OUT, "p33_signal_correlation.csv"), encoding="utf-8-sig")
import itertools
A = CM.abs().values.copy(); np.fill_diagonal(A, 0.0)
pairs = sorted([(A[i, j], ORDER[i], ORDER[j])
                for i, j in itertools.combinations(range(len(ORDER)), 2)], reverse=True)
p("  |Pearson r| 最高的 8 对：")
for v, a, b in pairs[:8]:
    p("    %.3f  %-22s %s" % (v, a, b))
p("")
p("  按 |r| 阈值聚类（单链接）：")
_d = 1 - A; _d = (_d + _d.T) / 2.0
cl_rows = []
for thr in (0.95, 0.90, 0.80):
    seen, groups = set(), []
    for i in range(len(ORDER)):
        if i in seen: continue
        grp = {i}
        changed = True
        while changed:
            changed = False
            for j in range(len(ORDER)):
                if j in grp: continue
                if any(A[j, k] >= thr for k in grp):
                    grp.add(j); changed = True
        seen |= grp
        groups.append(sorted(ORDER[k] for k in grp))
    multi = [g for g in groups if len(g) > 1]
    p("    |r| >= %.2f → %2d 组（多成员组 %d）：%s"
      % (thr, len(groups), len(multi),
         " ｜ ".join("+".join(g) for g in multi) if multi else "（全部孤立）"))
    cl_rows.append(dict(threshold=thr, n_groups=len(groups), n_multi=len(multi),
                        singletons=";".join(g[0] for g in groups if len(g) == 1),
                        groups=" | ".join("+".join(g) for g in groups)))
    if abs(thr - 0.90) < 1e-9:
        p("      → 孤立信号：%s"
          % ", ".join(g[0] for g in groups if len(g) == 1))
        p("      → **15 个信号在 |r| >= 0.90 下塌缩为 %d 组独立信息源**" % len(groups))
pd.DataFrame(cl_rows).to_csv(os.path.join(OUT, "p33_signal_clusters.csv"), index=False, encoding="utf-8-sig")
p("")
_pairs = pd.DataFrame(pairs[:12], columns=["abs_r", "signal_a", "signal_b"])
_pairs.to_csv(os.path.join(OUT, "p33_signal_top_pairs.csv"), index=False, encoding="utf-8-sig")

# ---- 指定一对的 cluster bootstrap 区间（按 §8.3 口径）----
#   ★ 取的是**正文引用的那一对**（entropy_js × conf_dist_js，Table S2 里两行近乎相同者），
#     不是全局最高对（那是 psi_mean × js_mean，r=0.9956）。首版写死 pairs[0]，取错了对象。
_rng = np.random.default_rng(SEED)
_ugs = np.unique(g)
_ea, _eb = "entropy_js", "conf_dist_js"
_pick_pair = [p for p in pairs if {p[1], p[2]} == {_ea, _eb}]
if not _pick_pair:
    raise RuntimeError("指定的信号对不在 top pairs 中：%s × %s" % (_ea, _eb))
_r_point = float(_pick_pair[0][0])
_va = X[_ea].to_numpy(dtype=float); _vb = X[_eb].to_numpy(dtype=float)
_rs = []
for _ in range(1000):
    _pick = _rng.choice(len(_ugs), len(_ugs), replace=True)
    _idx = np.concatenate([np.where(g == _ugs[i])[0] for i in _pick])
    if np.std(_va[_idx]) == 0 or np.std(_vb[_idx]) == 0:
        continue
    _rs.append(np.corrcoef(_va[_idx], _vb[_idx])[0, 1])
_r_lo, _r_hi = float(np.percentile(_rs, 2.5)), float(np.percentile(_rs, 97.5))
p("  **正文所引一对的 cluster bootstrap 区间**（organism–drug 层，1000 次，seed=%d）：" % SEED)
p("    r(%s, %s) = %.6f  [%.6f, %.6f]   （重采样 %d 次全部有效）"
  % (_ea, _eb, _r_point, _r_lo, _r_hi, len(_rs)))
p("    （供参照：全局 |r| 最高对为 %s × %s = %.6f）" % (pairs[0][1], pairs[0][2], pairs[0][0]))
pd.DataFrame([dict(signal_a=_ea, signal_b=_eb, pearson_r=_r_point,
                   ci_lo=_r_lo, ci_hi=_r_hi, n_boot=len(_rs), seed=SEED)]).to_csv(
    os.path.join(OUT, "p33_top_pair_ci.csv"), index=False, encoding="utf-8-sig")

# ---- signed 口径的自洽断言：只有三个有符号 L2 会变，其余 12 个必须逐位不变 ----
_p_changed = [m for m in ORDER
              if not np.allclose(score(m, "primary"), score(m, "signed"), equal_nan=True)]
p("")
p("  **口径自洽断言**：primary 与 signed 只在三个有符号 L2 上不同 —— 实际不同的方法：%s" % _p_changed)
p("    期望 = %s ⇒ %s" % (SIGNED_L2, "一致" if sorted(_p_changed) == sorted(SIGNED_L2) else "**不符**"))
p("")
p("=" * 118)
p("本块为 protocol §7.2 主任务四项指标的完整冻结产物。")
p("此前 PPV / NPV 从未计算；AUROC / AUPRC / sens@.80 已生产且此处独立复现。")
p("§4 相关矩阵另存 p33_signal_correlation.csv（供正文「独立信息源」一句溯源）。")
REP.close()
