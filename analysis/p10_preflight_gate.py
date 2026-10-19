# -*- coding: utf-8 -*-
r"""§10 预检闸门 —— 在启动 Phase B 之前验证缓存正确

**为什么需要这道闸门**：上一轮 §10 失败，是因为 p10a 自行重写了 build_hist，
202/202 task 的 source model 与 p6b 不一致。第一道闸门（p10a 内）拦的是模型；
这道拦的是 **signal 通路**。

**可验证的等式**：§10 的 N = full 档 target 取全量，且 PSI/JS/prediction-space 的
bin 边界与 source 参考直方图是用 source-train ref + full target 固定下来的，
因此这些 signal 在 N=full 时**必须逐 task 精确复现 p6b**。

**不适用**（§10.1 有意改变了 source reference，不应相等）：
  domain_auc（source 固定 5,000）、mmd（source 固定 1,000）、
  SUDO / AUROC_IOW / AUROC_OM / AUROC_DR（含随机抽样）

容差 1e-9。任一 feature 不通过则 raise SystemExit(1)。
"""
import os, pickle
import numpy as np, pandas as pd

BASE = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))); OUT = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(OUT, "p10_cache")
FEATS = ["age", "set_IP", "set_OP", "set_ER", "sp_URINE", "sp_BLOOD",
         "prior_any", "prior_cnt", "t_15-30", "t_31-90", "t_91-180", "t_181-365",
         "past_S", "past_I", "past_R"]

def psijs(r, a):
    r = r.astype(float); a = a.astype(float)
    rp = r / max(r.sum(), 1); ap = a / max(a.sum(), 1); eps = 1e-6
    rpc = np.clip(rp, eps, 1); apc = np.clip(ap, eps, 1)
    psi = float(np.sum((apc - rpc) * np.log(apc / rpc)))
    M = 0.5 * (rp + ap)
    def kl(x, y):
        m = x > 0; return float(np.sum(x[m] * np.log2(x[m] / np.maximum(y[m], 1e-12))))
    return psi, 0.5 * kl(rp, M) + 0.5 * kl(ap, M)

def ent(pv):
    pv = np.clip(pv, 1e-12, 1 - 1e-12); return -(pv * np.log(pv) + (1 - pv) * np.log(1 - pv))

B = pd.read_csv(os.path.join(OUT, "p6b_label_free_signals_blinded.csv")).set_index("task_id")
r = open(os.path.join(OUT, "report_p10_preflight.txt"), "w", encoding="utf-8")
def p(*a):
    line = " ".join(str(x) for x in a); print(line, flush=True); r.write(line + "\n"); r.flush()

p("=" * 104)
p("§10 预检闸门 —— N=full 档 signal 通路 vs p6b（容差 1e-9）")
p("=" * 104)

tids = sorted(f[:-4] for f in os.listdir(CACHE) if f.endswith(".pkl"))
diffs = {k: [] for k in ["psi_mean", "psi_max", "js_mean", "js_max", "pred_mean_shift",
                         "entropy_mean_shift", "entropy_js", "conf_dist_js", "pred_pos_prev_shift"]}
per_feat = {f: [] for f in FEATS}
bad = []
for tid in tids:
    if tid not in B.index: bad.append((tid, "not in p6b")); continue
    A = pickle.load(open(os.path.join(CACHE, f"{tid}.pkl"), "rb"))
    n = A["n_target"]; idx = np.arange(n)          # N=full
    pt_s = A["pt_full"]; Dt_std_s = A["Dt_std"]
    psis, jss = [], []
    for j, f in enumerate(FEATS):
        a = np.bincount(A["tgt_binidx"][j], minlength=len(A["psi_refhist"][j]))
        pv, jv = psijs(A["psi_refhist"][j], a)
        psis.append(pv); jss.append(jv)
        per_feat[f].append(abs(pv - B.loc[tid, "psi_" + f]))
    pa, ja = np.array(psis, float), np.array(jss, float)
    got = {
        "psi_mean": float(np.nanmean(pa)), "psi_max": float(np.nanmax(pa)),
        "js_mean": float(np.nanmean(ja)), "js_max": float(np.nanmax(ja)),
        "pred_mean_shift": float(pt_s.mean() - A["ps"].mean()),
        "entropy_mean_shift": float(ent(pt_s).mean() - A["ps_ent"].mean()),
        "entropy_js": psijs(A["jsent"]["refhist"], np.bincount(A["jsent"]["tgtidx"], minlength=len(A["jsent"]["refhist"])))[1],
        "conf_dist_js": psijs(A["jsconf"]["refhist"], np.bincount(A["jsconf"]["tgtidx"], minlength=len(A["jsconf"]["refhist"])))[1],
        "pred_pos_prev_shift": float((pt_s > 0.5).mean() - (A["ps"] > 0.5).mean()),
    }
    for k, v in got.items():
        diffs[k].append(abs(v - float(B.loc[tid, k])))

p(f"{'signal':<24}{'max|diff|':>14}{'mean|diff|':>14}{'不一致 task 数':>16}")
allok = True
for k, v in diffs.items():
    v = np.array(v); nb = int((v > 1e-9).sum())
    if nb: allok = False
    p(f"{k:<24}{v.max():>14.3e}{v.mean():>14.3e}{nb:>16}")

p("")
p(f"{'per-feature PSI':<24}{'max|diff|':>14}")
for f in FEATS:
    v = np.array(per_feat[f]); nb = int((v > 1e-9).sum())
    if nb: allok = False
    p(f"  {f:<22}{v.max():>14.3e}   (不一致 {nb})")

p("")
if bad: p(f"  额外问题：{bad[:5]}")
p(f"  覆盖 task = {len(diffs['psi_mean'])}")
p("")
if allok:
    p("  *** 预检通过 —— N=full 档 L1c/L2 signal 与 p6b 逐 task 精确一致。可以启动 Phase B。 ***")
else:
    p("  *** 预检失败 —— 缓存 signal 通路与 p6b 不一致，禁止启动 Phase B。 ***")
r.close()
print("done" if allok else "FAILED")
raise SystemExit(0 if allok else 1)
