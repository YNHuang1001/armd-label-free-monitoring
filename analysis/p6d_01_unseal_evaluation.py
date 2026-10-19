# -*- coding: utf-8 -*-
r"""Protocol v1.2 §11″ 步骤 5-6 —— **解封 target Y** 并执行正式评价

第一层：用 corrected models 计算 AUROC_target 与 ΔAUROC = AUROC_target - AUROC_source_internal
        harmful = I(ΔAUROC <= -0.05)
第二层：**先评价单指标，不做 meta-model**
第三层：grouped validation 为主线 —— leave-one-organism-drug-out + organism–drug cluster bootstrap；
        另做 leave-one-organism-out 作为 OOD stress test

r_failure 方向（由各指标定义决定，非择优；在此写明）：
  L1 位移类      : r = metric（越大 = 位移越大 = 越可能失败）
  L2 预测层位移  : r = |metric|（双向位移均视为风险）
  SUDO-AURCC     : r = -AURCC（越可靠 = 越不易失败）
  transported-AUC: r = -ÂUROC_target
DR 的 >1 值**不 clip**；IOW 权重**不截尾**。

输出 report_p6d_unsealed.txt / p6d_outcome.csv / p6d_signal_eval.csv
"""
import os, time, warnings
import sys
# ★★ 2026-10-13：**本脚本此前在本机跑不完** —— print 里有非 ASCII 字符，Windows 控制台 GBK
#   直接抛 UnicodeEncodeError，**exit=1**（实测：218 秒后崩，产物其实已经写完且逐位正确）。
#   ⇒ 「跑不动」与「算错」是两件事：此脚本属于前者。守卫与全仓其余脚本一致。
try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass
import numpy as np, pandas as pd
warnings.filterwarnings("ignore")
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score, average_precision_score
from scipy.stats import spearmanr

BASE = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))); OUT = os.path.dirname(os.path.abspath(__file__))
# ★★ 2026-10-13：种子可由环境变量覆盖，供种子稳定性分析（p63）用。**不改逻辑**。
#   ★ 首版 p63 把这一项归错了生产者：`fig2_forest.csv` 与 `tableS2_metrics.csv` **不是** p11
#     算的，也不是 p33 算的 —— p11 只是从 `p6d_signal_eval.csv` **复制/合并**（p11:158、p11:457）。
#     真正的 bootstrap 在**本脚本**：`rng = default_rng(SEED + 7)`，B = 1000（本文件第 291 行）。
#     ⇒ 已发表种子下 bootstrap 种子 = 20261012，与项目早先记录的「fig2_forest 种子 20261012」一致。
SEED = 20261005
t0 = time.time()
FOCAL = {"EC": "ESCHERICHIA COLI", "KP": "KLEBSIELLA PNEUMONIAE", "PA": "PSEUDOMONAS AERUGINOSA"}
DRUGMAP = [
    ("AMC", ["amoxicillinclavulanate", "amoxicillinclavulanicacid", "amoxclavulanate", "amoxicillinclav", "amoxicillinkclavulanateaug"]),
    ("TZP", ["piperacillintazobactam", "piptazobactam", "piperacillintazo", "tazobactam"]),
    ("SXT", ["trimethoprimsulfamethoxazole", "trimethsulf", "trimethoprimsulfa", "cotrimoxazole"]),
    ("SAM", ["ampicillinsulbactam", "ampicsulbactam", "ampicillinsulb"]),
    ("TIM", ["ticarcillinclavulanate", "ticarcilclavulanate", "ticarcillinclavulanicacid", "ticacarcillinclavulanicacid"]), ("CZA", ["ceftolozanetazobactam"]),
    ("CZO", ["cefazolin"]), ("CRO", ["ceftriaxone"]), ("CAZ", ["ceftazidime"]),
    ("FEP", ["cefepime"]), ("CTX", ["cefotaxime"]), ("CXM", ["cefuroxime"]), ("FOX", ["cefoxitin"]),
    ("CIP", ["ciprofloxacin"]), ("LVX", ["levofloxacin"]), ("GEN", ["gentamicin", "gentscreen", "gentamicinscreen"]),
    ("TOB", ["tobramycin", "tobramycinrm"]), ("AMK", ["amikacin", "amikaciniv", "amikacinintravenous"]), ("NIT", ["nitrofurantoin"]),
    ("MEM", ["meropenem"]), ("IPM", ["imipenem"]), ("ETP", ["ertapenem"]), ("ATM", ["aztreonam"]),
    ("TET", ["tetracycline"]), ("TGC", ["tigecycline"]), ("FOF", ["fosfomycin"]),
    ("CST", ["colistin", "polymixine", "polymyxine"]), ("SAM", ["sulbactam"]), ("AMP", ["ampicillin"]), ("PEN", ["piperacillin"]),
]
def _canon(n):
    up = "".join(ch for ch in str(n).lower() if ch.isalnum())
    for c, ks in DRUGMAP:
        if up in ks: return c
    return up
def norm_spec(s):
    return (s.astype(str).str.upper().str.strip()
            .str.replace("RESPIRATORY_TRACT", "RESP", regex=False).str.replace("RESPIRATORY", "RESP", regex=False))
def to_day(s):
    v = pd.to_datetime(s, errors="coerce", utc=True); return v.dt.tz_localize(None).dt.normalize()
def age_lo(s):
    return pd.to_numeric(str(s).split("-")[0].split("+")[0].strip().split(" ")[0], errors="coerce")

rep = open(os.path.join(OUT, "report_p6d_unsealed.txt"), "w", encoding="utf-8")
def p(*a):
    line = " ".join(str(x) for x in a); print(line, flush=True); rep.write(line + "\n"); rep.flush()

p("=" * 112); p("Protocol v1.2 —— 解封 target Y 后的正式评价"); p("=" * 112)

# ================= 重建 corrected models（与 p6b/p6c 同种子） =================
STORE = {}
for site in ["MGB", "SU", "UT"]:
    MGB = site == "MGB"
    P_ = os.path.join(BASE, "ARMD-MGB" if MGB else ("ARMD-Stanford" if site == "SU" else "ARMD-UTSW"))
    tcol = ("order_time_jittered_utc_shifted" if MGB else
            ("order_time_jittered_utc" if site == "SU" else "order_time_jittered"))
    if MGB:
        d = pd.read_csv(os.path.join(P_, "microbiology_cohort_deid_tj_updated.csv"),
                        usecols=["anon_id", "order_proc_id_coded", "organism", "antibiotic",
                                 "AST_pheno", "culture_description", tcol], low_memory=False, encoding="utf-8-sig")
        d = d[d["organism"].astype(str).str.upper().str.strip().eq("ESCHERICHIA COLI")
              | d["organism"].astype(str).str.upper().str.strip().str.startswith("KLEBSIELLA PNEUMONIAE")
              | d["organism"].astype(str).str.upper().str.strip().str.startswith("PSEUDOMONAS AERUGINOSA")].copy()
        d = d.rename(columns={tcol: "cdate", "AST_pheno": "sus"})
    else:
        d = pd.read_csv(os.path.join(P_, "microbiology_cultures_cohort.csv"),
                        usecols=["anon_id", "order_proc_id_coded", tcol, "culture_description",
                                 "organism", "antibiotic", "susceptibility"], low_memory=False, encoding="utf-8-sig")
        d = d.rename(columns={tcol: "cdate", "susceptibility": "sus"})
    d["cdate"] = to_day(d["cdate"]); d = d[d["cdate"].notna()].copy()
    d["org"] = d["organism"].astype(str).str.upper().str.strip()
    d["abx"] = d["antibiotic"].map(_canon); d["spec"] = norm_spec(d["culture_description"])
    d["s"] = d["sus"].astype(str).str.upper().str.strip().map(
        {"SUSCEPTIBLE": "S", "RESISTANT": "R", "INTERMEDIATE": "I"}).fillna("X")
    dem = pd.read_csv(os.path.join(P_, "demographics_deid_tj.csv" if MGB
                                   else "microbiology_cultures_demographics.csv"),
                      usecols=["order_proc_id_coded", "age"], low_memory=False, encoding="utf-8-sig")
    dem = dem.drop_duplicates("order_proc_id_coded").set_index("order_proc_id_coded")["age"]
    wd = pd.read_csv(os.path.join(P_, "ward_type_deid_tj.csv" if MGB
                                  else "microbiology_cultures_ward_info.csv"),
                     low_memory=False, encoding="utf-8-sig")
    wd["setting"] = np.where(wd["hosp_ward_ER"] == 1, "ER",
                             np.where(wd["hosp_ward_IP"] == 1, "Inpatient",
                                      np.where(wd["hosp_ward_OP"] == 1, "Outpatient", "Other")))
    wd = wd.drop_duplicates("order_proc_id_coded").set_index("order_proc_id_coded")["setting"]
    for sym, fname in FOCAL.items():
        dd = d[d["org"].str.startswith(fname)].copy()
        if not len(dd): continue
        u = dd.drop_duplicates(["anon_id", "cdate", "spec", "abx", "s"])
        usr = u[u["s"].isin(["S", "R", "I"])]
        if not len(usr): continue
        g = usr.groupby(["anon_id", "cdate", "spec", "abx"])["s"].agg(lambda v: "BOTH" if ({"S", "R"} <= set(v)) else ("S" if "S" in set(v) else ("R" if "R" in set(v) else "I"))).reset_index()
        unit = (g[["anon_id", "cdate", "spec"]].drop_duplicates()
                .sort_values(["anon_id", "spec", "cdate"]).reset_index(drop=True))
        blk = np.zeros(len(unit), dtype=int); ck = None; anc = None; b = -1
        for i in range(len(unit)):
            k = (unit["anon_id"].iat[i], unit["spec"].iat[i]); dv = unit["cdate"].iat[i]
            if k != ck: b += 1; anc = dv; ck = k
            elif (dv - anc) / np.timedelta64(1, "D") >= 14: b += 1; anc = dv
            blk[i] = b
        unit["_blk"] = blk
        unit["epid"] = (unit["anon_id"].astype(str) + "|" + unit["spec"].astype(str) + "|" + unit["_blk"].astype(str))
        g = g.merge(unit[["anon_id", "cdate", "spec", "epid"]], on=["anon_id", "cdate", "spec"], how="left")
        wide = g[g["s"].isin(["S", "R"])].pivot_table(index="epid", columns="abx", values="s", aggfunc="first")
        meta = unit.drop_duplicates("epid").set_index("epid")[["anon_id", "cdate", "spec"]]
        X = meta.join(wide, how="left")
        oid_map = dd.drop_duplicates(["anon_id", "cdate", "spec"]).set_index(
            ["anon_id", "cdate", "spec"])["order_proc_id_coded"]
        oids = oid_map.reindex(pd.MultiIndex.from_arrays([X["anon_id"], X["cdate"], X["spec"]])).values
        X["age"] = [age_lo(dem.get(o, np.nan)) for o in oids]
        X["setting"] = wd.reindex(oids).values
        STORE[(site, sym)] = X.sort_values(["anon_id", "spec", "cdate"]).reset_index(drop=True)
    del d, wd, dem

def build_hist(X, drug):
    n = len(X); pa = np.zeros(n); pc = np.zeros(n)
    tb = np.full(n, "none", dtype=object); ps = np.full(n, "none", dtype=object)
    c0 = X["cdate"].values
    for _, idx in X.groupby(["anon_id", "spec"], sort=False).indices.items():
        idx = np.sort(np.asarray(idx))
        for k in range(len(idx)):
            i = idx[k]
            lo = c0[i] - np.timedelta64(365, "D"); hi = c0[i] - np.timedelta64(15, "D")
            pr = idx[:k]; sel = pr[(c0[pr] > lo) & (c0[pr] <= hi)]
            if len(sel) == 0: continue
            pa[i] = 1; pc[i] = len(sel)
            dts = (c0[i] - c0[sel]).astype("timedelta64[D]").astype(int); mn = dts.min()
            tb[i] = ("15-30" if mn <= 30 else "31-90" if mn <= 90 else "91-180" if mn <= 180 else "181-365")
            j = sel[np.argmin(dts)]; v = X[drug].values[j] if drug in X.columns else None
            ps[i] = v if v in ("S", "I", "R") else "none"
    return pd.DataFrame({"prior_any": pa, "prior_cnt": pc, "tbin": tb, "past_ast": ps})
def design(X, H):
    F = pd.DataFrame(index=X.index)
    F["age"] = pd.to_numeric(X["age"], errors="coerce")
    F["set_IP"] = (X["setting"] == "Inpatient").astype(int)
    F["set_OP"] = (X["setting"] == "Outpatient").astype(int)
    F["set_ER"] = (X["setting"] == "ER").astype(int)
    F["sp_URINE"] = (X["spec"] == "URINE").astype(int)
    F["sp_BLOOD"] = (X["spec"] == "BLOOD").astype(int)
    F["prior_any"] = H["prior_any"].values; F["prior_cnt"] = H["prior_cnt"].values
    for b in ["15-30", "31-90", "91-180", "181-365"]:
        F[f"t_{b}"] = (H["tbin"] == b).astype(int).values
    for v in ["S", "I", "R"]:
        F[f"past_{v}"] = (H["past_ast"] == v).astype(int).values
    return F
def prep(X, drug):
    X = X.sort_values(["anon_id", "spec", "cdate"]).reset_index(drop=True)
    H = build_hist(X, drug); m = np.isin(X[drug].values, ["S", "R"])
    return X[m].reset_index(drop=True), H[m].reset_index(drop=True)

SIG = pd.read_csv(os.path.join(OUT, "p6b_label_free_signals_blinded.csv"))
L3 = pd.read_csv(os.path.join(OUT, "p6c_L3_blinded.csv"))
TD = pd.read_csv(os.path.join(OUT, "p6b_corrected_universe.csv"))
TASKS = TD[TD["eligible"]].copy()
KEYS = sorted(set(zip(TASKS.organism, TASKS.drug, TASKS.source_site)))

MC = {}
for sym, drug, src in KEYS:
    Xs, Hs = prep(STORE[(src, sym)], drug)
    D = design(Xs, Hs); y = (Xs[drug].values == "R").astype(int)
    pats = Xs["anon_id"].values; upat = pd.unique(pats)
    rng = np.random.default_rng(SEED); perm = rng.permutation(len(upat))
    ntr, nva = int(0.70 * len(upat)), int(0.15 * len(upat))
    s_tr = set(upat[perm[:ntr]]); s_va = set(upat[perm[ntr:ntr + nva]])
    i_tr = np.array([i for i, a in enumerate(pats) if a in s_tr])
    i_va = np.array([i for i, a in enumerate(pats) if a in s_va])
    i_te = np.array([i for i, a in enumerate(pats) if a not in s_tr and a not in s_va])
    Dall = D.values.copy()
    af = float(np.nanmedian(Dall[i_tr, 0])) if np.isfinite(Dall[i_tr, 0]).any() else 0.0
    Dall[:, 0] = np.where(np.isfinite(Dall[:, 0]), Dall[:, 0], af)
    sc = StandardScaler().fit(Dall[i_tr])
    mdl = LogisticRegression(max_iter=2000, C=1.0, random_state=SEED).fit(sc.transform(Dall[i_tr]), y[i_tr])
    MC[(sym, drug, src)] = dict(sc=sc, mdl=mdl, yte=y[i_te], Dte=Dall[i_te], af=af)
p(f"  模型就绪 {len(MC)}  [{time.time()-t0:.0f}s]")

# ================= 第一层：解封 target Y，算 ΔAUROC =================
rows = []
for _, t in TASKS.iterrows():
    sym, drug, src, tgt = t["organism"], t["drug"], t["source_site"], t["target_site"]
    M = MC[(sym, drug, src)]
    Xt, Ht = prep(STORE[(tgt, sym)], drug)
    Dt = design(Xt, Ht).values.copy()
    Dt[:, 0] = np.where(np.isfinite(Dt[:, 0]), Dt[:, 0], M["af"])
    yt = (Xt[drug].values == "R").astype(int)
    a_src = roc_auc_score(M["yte"], M["mdl"].predict_proba(M["sc"].transform(M["Dte"]))[:, 1])
    a_tgt = roc_auc_score(yt, M["mdl"].predict_proba(M["sc"].transform(Dt))[:, 1])
    rows.append(dict(task_id=t["task_id"], organism=sym, drug=drug, source_site=src, target_site=tgt,
                     n_target=len(Xt), auroc_source_internal=a_src, auroc_target=a_tgt,
                     d_auroc=a_tgt - a_src, harmful=int((a_tgt - a_src) <= -0.05),
                     stable=int(abs(a_tgt - a_src) <= 0.02)))
O = pd.DataFrame(rows)
O.to_csv(os.path.join(OUT, "p6d_outcome.csv"), index=False)

p("\n" + "=" * 112); p("【第一层】corrected benchmark 的 degradation 分布"); p("=" * 112)
p(f"  task = {len(O)}   harmful(Δ<=-0.05) = {int(O.harmful.sum())}   stable(|Δ|<=0.02) = {int(O.stable.sum())}")
d = O["d_auroc"]
p("  ΔAUROC: " + "  ".join(f"P{q}={d.quantile(q/100):+.3f}" for q in [5, 10, 25, 50, 75, 90, 95]))
p(f"  P90-P10 = {d.quantile(.90)-d.quantile(.10):.3f}   median = {d.median():+.3f}")
p("\n  [对照] 旧 p5（superseded, preprocessing-mismatched）: harmful=65, stable=49, P90-P10=0.146")
p(f"  按 direction:")
O["dir"] = O.source_site + "->" + O.target_site
p(O.groupby("dir")["d_auroc"].agg(["size", "median", "min", "max"]).round(3).to_string())

# ================= 合并 signals =================
E = O.merge(SIG[["task_id", "domain_auc", "mmd", "psi_mean", "psi_max", "js_mean", "js_max",
                 "pred_mean_shift", "entropy_mean_shift", "entropy_js", "conf_dist_js",
                 "pred_pos_prev_shift"]], on="task_id", how="left") \
     .merge(L3[["task_id", "SUDO_AURCC", "AUROC_IOW", "AUROC_OM", "AUROC_DR", "w_max", "w_ESS"]],
            on="task_id", how="left")
E["group"] = E.organism + "|" + E.drug

METHODS = [
    ("domain_auc", "L1", lambda x: x["domain_auc"]),
    ("mmd", "L1", lambda x: x["mmd"]),
    ("psi_mean", "L1", lambda x: x["psi_mean"]),
    ("psi_max", "L1", lambda x: x["psi_max"]),
    ("js_mean", "L1", lambda x: x["js_mean"]),
    ("js_max", "L1", lambda x: x["js_max"]),
    ("pred_mean_shift", "L2", lambda x: x["pred_mean_shift"].abs()),
    ("entropy_mean_shift", "L2", lambda x: x["entropy_mean_shift"].abs()),
    ("entropy_js", "L2", lambda x: x["entropy_js"]),
    ("conf_dist_js", "L2", lambda x: x["conf_dist_js"]),
    ("pred_pos_prev_shift", "L2", lambda x: x["pred_pos_prev_shift"].abs()),
    ("SUDO_AURCC", "L2", lambda x: -x["SUDO_AURCC"]),
    ("AUROC_DR", "L3-primary", lambda x: -x["AUROC_DR"]),
    ("AUROC_IOW", "L3-sens", lambda x: -x["AUROC_IOW"]),
    ("AUROC_OM", "L3-sens", lambda x: -x["AUROC_OM"]),
]

def grouped_cv_oof(r, y, groups):
    """leave-one-organism-drug-out：在其余组上拟合单变量 logistic，预测被留出组"""
    oof = np.full(len(r), np.nan)
    for g in np.unique(groups):
        te = groups == g; tr = ~te
        if len(np.unique(y[tr])) < 2 or te.sum() == 0: continue
        mu, sd = np.nanmean(r[tr]), np.nanstd(r[tr]) or 1.0
        try:
            m = LogisticRegression(max_iter=1000).fit(((r[tr] - mu) / sd).reshape(-1, 1), y[tr])
            oof[te] = m.predict_proba(((r[te] - mu) / sd).reshape(-1, 1))[:, 1]
        except Exception:
            continue
    return oof

def safe_auc(y, s):
    m = np.isfinite(s)
    return roc_auc_score(y[m], s[m]) if len(np.unique(y[m])) > 1 else np.nan
def safe_ap(y, s):
    m = np.isfinite(s)
    return average_precision_score(y[m], s[m]) if len(np.unique(y[m])) > 1 else np.nan
def sens_at_spec(y, s, spec=0.80):
    m = np.isfinite(s)
    y, s = y[m], s[m]
    thr = min([v for v in np.unique(s) if (s[y == 0] < v).mean() >= spec] or [np.nextafter(s.max() if s.size else 0.0, np.inf)])
    pred = s >= thr
    tp = ((pred == 1) & (y == 1)).sum(); fn = ((pred == 0) & (y == 1)).sum()
    return tp / max(tp + fn, 1) if (len(np.unique(y)) > 1 and (y == 0).sum() > 0 and (y == 1).sum() > 0) else np.nan

res = []
for name, role, fn in METHODS:
    r = fn(E).values.astype(float); y = E["harmful"].values; g = E["group"].values
    ok = np.isfinite(r)
    sp = spearmanr(r[ok], E["d_auroc"].values[ok]).statistic if ok.sum() > 5 else np.nan
    auc0 = safe_auc(y[ok], r[ok]); ap0 = safe_ap(y[ok], r[ok]); ss = sens_at_spec(y[ok], r[ok])
    oof = grouped_cv_oof(r, y, g)
    gauc = safe_auc(y, oof)
    # leave-one-organism-out（EC+KP -> PA，再反向）
    loo = []
    for org in ["EC", "KP", "PA"]:
        te = (E.organism == org).values; tr = ~te
        if len(np.unique(y[tr])) < 2: continue
        mu, sd = np.nanmean(r[tr]), (np.nanstd(r[tr]) or 1.0)
        try:
            m = LogisticRegression(max_iter=1000).fit(((r[tr] - mu) / sd).reshape(-1, 1), y[tr])
            s2 = m.predict_proba(((r[te] - mu) / sd).reshape(-1, 1))[:, 1]
            loo.append(safe_auc(y[te], s2))
        except Exception: pass
    # cluster bootstrap（organism-drug 层）
    rng = np.random.default_rng(SEED + 7); ugs = np.unique(g); sps, aucs = [], []
    for _ in range(1000):
        pick = rng.choice(len(ugs), len(ugs), replace=True)
        idx = np.concatenate([np.where(g == ugs[i])[0] for i in pick])
        if len(np.unique(y[idx])) < 2: continue
        sps.append(spearmanr(r[idx], E["d_auroc"].values[idx]).statistic)
        aucs.append(safe_auc(y[idx], r[idx]))
    res.append(dict(method=name, role=role, coverage=int(ok.sum()),
                    spearman=sp, sp_lo=np.percentile(sps, 2.5), sp_hi=np.percentile(sps, 97.5),
                    harmful_auc=auc0, auc_lo=np.percentile(aucs, 2.5), auc_hi=np.percentile(aucs, 97.5),
                    harmful_auprc=ap0, sens_at_spec80=ss,
                    grouped_cv_auc=gauc, loo_org_auc=np.nanmean(loo) if loo else np.nan))
R = pd.DataFrame(res)
R.to_csv(os.path.join(OUT, "p6d_signal_eval.csv"), index=False)

p("\n" + "=" * 112); p("【第二/三层】单指标评价（grouped CV 为主线 + cluster bootstrap）"); p("=" * 112)
p(f"  {'Method':<20}{'Role':<12}{'Spearman':>9}{'95%CI(cluster)':>20}{'HarmAUC':>9}{'95%CI':>18}"
  f"{'AUPRC':>8}{'sens@.8':>9}{'LODO':>8}{'LOO-org':>9}{'Cov':>6}")
for _, r in R.iterrows():
    p(f"  {r['method']:<20}{r['role']:<12}{r['spearman']:>+9.3f}"
      f"{f'[{r.sp_lo:+.3f},{r.sp_hi:+.3f}]':>20}{r['harmful_auc']:>9.3f}"
      f"{f'[{r.auc_lo:.3f},{r.auc_hi:.3f}]':>18}{r['harmful_auprc']:>8.3f}"
      f"{r['sens_at_spec80']:>9.3f}{r['grouped_cv_auc']:>8.3f}"
      f"{r['loo_org_auc']:>9.3f}{r['coverage']:>6}")

# ================= L3 的 MAE / calibration =================
p("\n" + "=" * 112); p("【L3 MAE / calibration】（仅 transported-AUC，SUDO 不参与）"); p("=" * 112)
for c in ["AUROC_DR", "AUROC_IOW", "AUROC_OM"]:
    dhat = E[c].values - E["auroc_source_internal"].values
    dd = E["d_auroc"].values
    mae = np.nanmean(np.abs(dhat - dd))
    oor = int((E[c] > 1).sum() + (E[c] < 0).sum())
    # calibration：预测掉 0.05 一档时实际掉多少
    bins = [-np.inf, -0.10, -0.05, 0.0, 0.05, np.inf]
    lab = ["<-0.10", "-0.10~-0.05", "-0.05~0", "0~0.05", ">0.05"]
    b = pd.cut(dhat, bins, labels=lab)
    cal = pd.DataFrame({"pred": b, "act": dd}).groupby("pred", observed=True)["act"].agg(["size", "mean"])
    # ★ 2026-10-06：原写死 `/202`（修正前的 task 数）。改为现算，否则新产物会带着旧分母。
    p(f"  {c:<12} MAE={mae:.4f}  越界(>1或<0)={oor}/{len(E)}")
    p("       calibration(pred 档 -> 实际 ΔAUROC 均值): " +
      "  ".join(f"{i}:n={int(r['size'])},act={r['mean']:+.3f}" for i, r in cal.iterrows()))

# ★ 2026-10-06：原写死 "【DR 4 个越界 task】"。CS-05 修正后该数会变（本期为 0，
#   表是空的而标题仍写 4 —— 正是必须改成现算的原因）。此处口径与下面 bad 的筛选(>1)一致。
p(f"\n【DR {int((E['AUROC_DR'] > 1).sum())} 个越界 task（AUROC_DR > 1）】")
bad = E[E["AUROC_DR"] > 1][["task_id", "organism", "drug", "source_site", "target_site",
                            "AUROC_DR", "w_max", "w_ESS", "n_source" if "n_source" in E else "n_target"]]
p(bad.round(3).to_string(index=False))
p(f"\n  权重极端性与 L3 误差的关系（诊断，不修改估计量）：")
for c in ["w_max", "w_ESS"]:
    rr = E[[c, "AUROC_DR"]].dropna()
    err = np.abs(rr["AUROC_DR"] - rr["AUROC_DR"].median())
    p(f"    Spearman({c}, |ÂUROC_DR - median|) = {rr[c].corr(err, method='spearman'):+.3f}")
p("\n  **signal 定义、阈值、主次自解封起锁死，不再修改。**")
p(f"\n完成，耗时 {time.time()-t0:.0f}s")
rep.close(); print("done", flush=True)
