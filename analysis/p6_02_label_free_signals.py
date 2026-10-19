# -*- coding: utf-8 -*-
r"""Protocol v1.0 —— 执行 §5.1 + §5.2（blinded）

只计算：
  §5.1 X-space : domain_auc (+CI) / mmd / psi_* / js_*
  §5.2 pred-space: pred_mean_shift / entropy_mean_shift / entropy_js / conf_dist_js / pred_pos_prev_shift

不做：SUDO / meta-model / 与 ΔAUROC 的任何相关性 / 指标择优。
数据来源：**只读 p6_task_manifest_blinded.csv**，不读 p5 的任何性能字段。

冻结实现
  domain classifier : LogisticRegression(class_weight='balanced'), 5-fold StratifiedKFold OOF AUROC
                      （禁止训练集 AUROC；全 202 tasks 同一实现，不逐 task 调）
  MMD               : RBF；bandwidth = pooled subsample 的 median heuristic（单一，不扫描）
                      每侧 subsample m=1000，重复 K=5，取均值（固定实现常数，未见结果前写死）
  PSI / JS          : source **train** 的分位数切 10 等频箱（连续）；类别变量用 source train 取值集，
                      target 未见类别单列 UNSEEN；边界锁死
  预处理            : 一律 source-train 定义（StandardScaler 只 fit 在 source train）

输出 p6_label_free_signals_blinded.csv（**正好 202 行**）
     report_p6_label_free_blinded.txt
"""
import os, sys, time, warnings
import numpy as np, pandas as pd
warnings.filterwarnings("ignore")
from scipy.spatial.distance import pdist, cdist
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score

BASE = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))); OUT = os.path.dirname(os.path.abspath(__file__))
t0 = time.time()
SEED = 20261005
MMD_M, MMD_K = 1000, 5          # 冻结实现常数（未见任何结果前写死）
NBIN = 10

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

rep = open(os.path.join(OUT, "report_p6_label_free_blinded.txt"), "w", encoding="utf-8")
def p(*a):
    line = " ".join(str(x) for x in a); print(line, flush=True); rep.write(line + "\n"); rep.flush()

p("=" * 112); p("Protocol v1.0 §5.1 + §5.2 —— blinded label-free signals"); p("=" * 112)
p(f"冻结实现常数：MMD subsample m={MMD_M}×{MMD_K} repeats；PSI/JS 箱数={NBIN}（source-train 分位数）")

MAN = pd.read_csv(os.path.join(OUT, "p6_task_manifest_blinded.csv"))
p(f"manifest = {len(MAN)} tasks  （来源：p6_task_manifest_blinded.csv；不含任何性能字段）")
assert len(MAN) == 202

# ============ 数据构建（与 p5 完全一致，唯一差异 = episode 用 anchor）============
STORE = {}
for site in ["MGB", "SU", "UT"]:
    MGB = site == "MGB"
    P_ = os.path.join(BASE, "ARMD-MGB" if MGB else ("ARMD-Stanford" if site == "SU" else "ARMD-UTSW"))
    tcol = ("order_time_jittered_utc_shifted" if MGB else
            ("order_time_jittered_utc" if site == "SU" else "order_time_jittered"))
    if MGB:
        d = pd.read_csv(os.path.join(P_, "microbiology_cohort_deid_tj_updated.csv"),
                        usecols=["anon_id", "order_proc_id_coded", "organism", "antibiotic",
                                 "AST_pheno", "culture_description", tcol],
                        low_memory=False, encoding="utf-8-sig")
        d = d[d["organism"].astype(str).str.upper().str.strip().eq("ESCHERICHIA COLI")
              | d["organism"].astype(str).str.upper().str.strip().str.startswith("KLEBSIELLA PNEUMONIAE")
              | d["organism"].astype(str).str.upper().str.strip().str.startswith("PSEUDOMONAS AERUGINOSA")].copy()
        d = d.rename(columns={tcol: "cdate", "AST_pheno": "sus"})
    else:
        d = pd.read_csv(os.path.join(P_, "microbiology_cultures_cohort.csv"),
                        usecols=["anon_id", "order_proc_id_coded", tcol, "culture_description",
                                 "organism", "antibiotic", "susceptibility"],
                        low_memory=False, encoding="utf-8-sig")
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
    p(f"  [{site}] STORE 就绪  [{time.time()-t0:.0f}s]")

# ============ 特征构造 ============
def build_hist(X, drug):
    n = len(X)
    prior_any = np.zeros(n); prior_cnt = np.zeros(n)
    tbin = np.full(n, "none", dtype=object); past = np.full(n, "none", dtype=object)
    c0 = X["cdate"].values
    for _, idx in X.groupby(["anon_id", "spec"], sort=False).indices.items():
        idx = np.sort(np.asarray(idx))
        for k in range(len(idx)):
            i = idx[k]
            lo = c0[i] - np.timedelta64(365, "D"); hi = c0[i] - np.timedelta64(15, "D")
            pr = idx[:k]; sel = pr[(c0[pr] > lo) & (c0[pr] <= hi)]
            if len(sel) == 0: continue
            prior_any[i] = 1; prior_cnt[i] = len(sel)
            dts = (c0[i] - c0[sel]).astype("timedelta64[D]").astype(int); mn = dts.min()
            tbin[i] = ("15-30" if mn <= 30 else "31-90" if mn <= 90 else "91-180" if mn <= 180 else "181-365")
            j = sel[np.argmin(dts)]
            v = X[drug].values[j] if drug in X.columns else None
            past[i] = v if v in ("S", "I", "R") else "none"
    return pd.DataFrame({"prior_any": prior_any, "prior_cnt": prior_cnt, "tbin": tbin, "past_ast": past})

FEATS = ["age", "set_IP", "set_OP", "set_ER", "sp_URINE", "sp_BLOOD",
         "prior_any", "prior_cnt", "t_15-30", "t_31-90", "t_91-180", "t_181-365",
         "past_S", "past_I", "past_R"]
CONT = {"age", "prior_cnt"}

def design(X, H):
    F = pd.DataFrame(index=X.index)
    F["age"] = pd.to_numeric(X["age"], errors="coerce")
    F["set_IP"] = (X["setting"] == "Inpatient").astype(int)
    F["set_OP"] = (X["setting"] == "Outpatient").astype(int)
    F["set_ER"] = (X["setting"] == "ER").astype(int)
    F["sp_URINE"] = (X["spec"] == "URINE").astype(int)
    F["sp_BLOOD"] = (X["spec"] == "BLOOD").astype(int)
    F["prior_any"] = H["prior_any"].values
    F["prior_cnt"] = H["prior_cnt"].values
    for b in ["15-30", "31-90", "91-180", "181-365"]:
        F[f"t_{b}"] = (H["tbin"] == b).astype(int).values
    for v in ["S", "I", "R"]:
        F[f"past_{v}"] = (H["past_ast"] == v).astype(int).values
    return F

def prep(X, drug):
    X = X.sort_values(["anon_id", "spec", "cdate"]).reset_index(drop=True)
    H = build_hist(X, drug)
    m = np.isin(X[drug].values, ["S", "R"])
    return X[m].reset_index(drop=True), H[m].reset_index(drop=True)

# ============ 冻结的度量实现 ============
def psijs_from_bins(ref_r, act_r):
    r = ref_r.astype(float); a = act_r.astype(float)
    rp = r / max(r.sum(), 1); ap = a / max(a.sum(), 1)
    eps = 1e-6
    rpc = np.clip(rp, eps, 1.0); apc = np.clip(ap, eps, 1.0)
    psi = float(np.sum((apc - rpc) * np.log(apc / rpc)))
    M = 0.5 * (rp + ap)
    def kl(x, y):
        msk = x > 0
        return float(np.sum(x[msk] * np.log2(x[msk] / np.maximum(y[msk], 1e-12))))
    js = 0.5 * kl(rp, M) + 0.5 * kl(ap, M)
    return psi, js

def feat_psijs(ref, act, cont):
    ref = np.asarray(ref); act = np.asarray(act)
    if cont:
        qs = np.unique(np.quantile(ref, np.linspace(0, 1, NBIN + 1)))
        if len(qs) < 3:
            return np.nan, np.nan
        edges = qs.astype(float).copy(); edges[0] = -np.inf; edges[-1] = np.inf
        r = np.histogram(ref, bins=edges)[0]; a = np.histogram(act, bins=edges)[0]
    else:
        cats = np.unique(np.r_[np.unique(ref), np.unique(act)])
        r = np.bincount(np.searchsorted(cats, ref), minlength=len(cats))
        a = np.bincount(np.searchsorted(cats, act), minlength=len(cats))
    if r.sum() == 0 or a.sum() == 0: return np.nan, np.nan
    return psijs_from_bins(r, a)

def mmd_rbf(A, B, rng):
    vals = []
    for _ in range(MMD_K):
        ia = rng.choice(len(A), min(MMD_M, len(A)), replace=False)
        ib = rng.choice(len(B), min(MMD_M, len(B)), replace=False)
        a, b = A[ia], B[ib]
        P = np.vstack([a, b])
        s2 = float(np.median(pdist(P, "sqeuclidean")))
        if not np.isfinite(s2) or s2 <= 0: s2 = 1.0
        Kaa = np.exp(-pdist(a, "sqeuclidean") / s2)
        Kbb = np.exp(-pdist(b, "sqeuclidean") / s2)
        Kab = np.exp(-cdist(a, b, "sqeuclidean") / s2)
        n, m = len(a), len(b)
        vals.append(2 * Kaa.sum() / (n * (n - 1)) + 2 * Kbb.sum() / (m * (m - 1)) - 2 * Kab.mean())
    return float(np.mean(vals)), float(np.std(vals))

def js_cont(ref, act):
    qs = np.unique(np.quantile(ref, np.linspace(0, 1, NBIN + 1)))
    if len(qs) < 3: return np.nan
    edges = qs.astype(float).copy(); edges[0] = -np.inf; edges[-1] = np.inf
    _, j = psijs_from_bins(np.histogram(ref, bins=edges)[0], np.histogram(act, bins=edges)[0])
    return j

def ent(p):
    p = np.clip(p, 1e-12, 1 - 1e-12)
    return -(p * np.log(p) + (1 - p) * np.log(1 - p))

# ============ 主循环 ============
rows, fail_log = [], []
MODEL_CACHE = {}
for _, task in MAN.iterrows():
    sym, drug, src, tgt = task["organism"], task["drug"], task["source_site"], task["target_site"]
    rec = {c: np.nan for c in
           ["n_source", "n_target", "n_src_train", "n_src_val", "n_src_test",
            "domain_auc", "domain_auc_lo", "domain_auc_hi", "mmd", "mmd_sd",
            "psi_mean", "psi_max", "js_mean", "js_max",
            "pred_mean_shift", "entropy_mean_shift", "entropy_js", "conf_dist_js", "pred_pos_prev_shift"]}
    rec.update({"task_id": task["task_id"], "organism": sym, "drug": drug,
                "source_site": src, "target_site": tgt})
    reasons = []
    key = (src, sym, drug)
    try:
        if key not in MODEL_CACHE:
            Xs, Hs = prep(STORE[(src, sym)], drug)
            Ds = design(Xs, Hs)
            y = (Xs[drug].values == "R").astype(int)
            pats = Xs["anon_id"].values; upat = pd.unique(pats)
            rng = np.random.default_rng(SEED); perm = rng.permutation(len(upat))
            ntr, nva = int(0.70 * len(upat)), int(0.15 * len(upat))
            s_tr = set(upat[perm[:ntr]]); s_va = set(upat[perm[ntr:ntr + nva]])
            i_tr = np.array([i for i, a in enumerate(pats) if a in s_tr])
            i_va = np.array([i for i, a in enumerate(pats) if a in s_va])
            i_te = np.array([i for i, a in enumerate(pats) if a not in s_tr and a not in s_va])
            # 缺失处理必须 source-defined（protocol §4.2）：age 用 **source train 中位数** 填补，
            # 同一填补值应用到 source val/test、target、以及 domain classifier 的全部输入。
            # 注意：必须在**整块** Ds 上填补——domain classifier 用的是完整的 Ds，不是 Dtr。
            Ds_all = Ds.values.copy()
            age_miss = float(np.mean(~np.isfinite(Ds_all[:, 0])) * 100)
            age_fill = float(np.nanmedian(Ds_all[i_tr, 0])) if np.isfinite(Ds_all[i_tr, 0]).any() else 0.0
            Ds_all[:, 0] = np.where(np.isfinite(Ds_all[:, 0]), Ds_all[:, 0], age_fill)
            Dtr, Dva, Dte = Ds_all[i_tr], Ds_all[i_va], Ds_all[i_te]
            sc = StandardScaler().fit(Dtr)
            mdl = LogisticRegression(max_iter=2000, C=1.0, random_state=SEED)
            mdl.fit(sc.transform(Dtr), y[i_tr])
            ps = mdl.predict_proba(sc.transform(Dte))[:, 1]
            MODEL_CACHE[key] = dict(Ds=Ds_all, Dtr=Dtr, Dte=Dte, sc=sc, mdl=mdl, ps=ps,
                                    age_fill=age_fill, age_miss=age_miss,
                                    n_src=len(Xs), n_tr=len(i_tr), n_va=len(i_va), n_te=len(i_te))
        M = MODEL_CACHE[key]
        if drug not in STORE[(tgt, sym)].columns:
            raise ValueError(f"target 缺该药列: {tgt}/{sym}/{drug}")
        Xt, Ht = prep(STORE[(tgt, sym)], drug)
        Dt = design(Xt, Ht).values.copy()
        dt_age_miss = float(np.mean(~np.isfinite(Dt[:, 0])) * 100)
        Dt[:, 0] = np.where(np.isfinite(Dt[:, 0]), Dt[:, 0], M["age_fill"])   # source-train 填补值
        rec.update({"n_source": M["n_src"], "n_target": len(Xt),
                    "n_src_train": M["n_tr"], "n_src_val": M["n_va"], "n_src_test": M["n_te"],
                    "age_missing_src_pct": M["age_miss"], "age_missing_tgt_pct": dt_age_miss})

        # ---- §5.1a domain classifier（OOF，class_weight balanced）----
        Xd = np.vstack([M["Ds"], Dt])
        yd = np.r_[np.zeros(len(M["Ds"])), np.ones(len(Dt))]
        oof = np.zeros(len(yd))
        skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=SEED)
        for tr, te in skf.split(Xd, yd):
            dm = LogisticRegression(max_iter=2000, class_weight="balanced", random_state=SEED)
            dm.fit(Xd[tr], yd[tr]); oof[te] = dm.predict_proba(Xd[te])[:, 1]
        rec["domain_auc"] = float(roc_auc_score(yd, oof))
        rb = np.random.default_rng(SEED + 1); bs = []
        for _ in range(200):
            s = rb.integers(0, len(yd), len(yd))
            if len(np.unique(yd[s])) < 2: continue
            bs.append(roc_auc_score(yd[s], oof[s]))
        rec["domain_auc_lo"], rec["domain_auc_hi"] = np.percentile(bs, [2.5, 97.5])

        # ---- §5.1b MMD ----
        rm = np.random.default_rng(SEED + 2)
        mval, msd = mmd_rbf(M["Ds"], Dt, rm)
        rec["mmd"], rec["mmd_sd"] = mval, msd

        # ---- §5.1c PSI / JS（source TRAIN 为参考）----
        ref_train = M["Dtr"]; psis, jss = {}, {}
        for j, f in enumerate(FEATS):
            pv, jv = feat_psijs(ref_train[:, j], Dt[:, j], f in CONT)
            psis[f] = pv; jss[f] = jv
            rec[f"psi_{f}"] = pv; rec[f"js_{f}"] = jv
        pv_arr = np.array(list(psis.values()), float); jv_arr = np.array(list(jss.values()), float)
        rec["psi_mean"] = float(np.nanmean(pv_arr)); rec["psi_max"] = float(np.nanmax(pv_arr))
        rec["js_mean"] = float(np.nanmean(jv_arr)); rec["js_max"] = float(np.nanmax(jv_arr))

        # ---- §5.2 prediction-space ----
        pt = M["mdl"].predict_proba(M["sc"].transform(Dt))[:, 1]
        rec["pred_mean_shift"] = float(pt.mean() - M["ps"].mean())
        rec["entropy_mean_shift"] = float(ent(pt).mean() - ent(M["ps"]).mean())
        rec["entropy_js"] = js_cont(ent(M["ps"]), ent(pt))
        rec["conf_dist_js"] = js_cont(np.maximum(M["ps"], 1 - M["ps"]), np.maximum(pt, 1 - pt))
        rec["pred_pos_prev_shift"] = float((pt > 0.5).mean() - (M["ps"] > 0.5).mean())
    except Exception as e:
        reasons.append(f"{type(e).__name__}: {e}")
    rec["status"] = "OK" if not reasons else "FAILED"
    rec["failure_reason"] = "" if not reasons else " | ".join(reasons)
    if reasons: fail_log.append((task["task_id"], rec["failure_reason"]))
    rows.append(rec)
    if len(rows) % 25 == 0: p(f"    ... {len(rows)}/202  [{time.time()-t0:.0f}s]")

DF = pd.DataFrame(rows)
FRONT = ["task_id", "organism", "drug", "source_site", "target_site",
         "n_source", "n_target", "n_src_train", "n_src_val", "n_src_test",
         "domain_auc", "domain_auc_lo", "domain_auc_hi", "mmd", "mmd_sd",
         "psi_mean", "psi_max", "js_mean", "js_max", "pred_mean_shift",
         "entropy_mean_shift", "entropy_js", "conf_dist_js", "pred_pos_prev_shift"]
p("\n[失败原因汇总]")
if fail_log:
    import collections
    for msg, k in collections.Counter(b.split(": ", 1)[-1] for _, b in fail_log).most_common(10):
        p(f"   x{k:<4} {msg[:100]}")
    p(f"   示例: {fail_log[0][0]} -> {fail_log[0][1][:160]}")
else:
    p("   无")
COLS = [c for c in FRONT if c in DF.columns] + [c for c in DF.columns if c not in FRONT]
DF = DF[COLS]
DF.to_csv(os.path.join(OUT, "p6_label_free_signals_blinded.csv"), index=False)
p(f"\n[写出] p6_label_free_signals_blinded.csv  行数 = {len(DF)}  （必须 = 202）")

# ============ QC（不揭盲）============
p("\n" + "=" * 112); p("【QC-1】覆盖率：202 为分母，逐指标"); p("=" * 112)
p(f"  status == OK : {int((DF['status']=='OK').sum())} / {len(DF)}")
for c in FRONT[10:] + [c for c in DF.columns if c.startswith(("psi_", "js_")) and c not in FRONT]:
    n = int(DF[c].notna().sum())
    if n < len(DF): p(f"    {c:<20} 可算 {n}/202")
p(f"  完全无缺失的指标数 = {sum(1 for c in DF.columns if DF[c].notna().all())}/{len(DF.columns)}")
if fail_log:
    p("\n  失败 task："); [p(f"    {a}: {b}") for a, b in fail_log[:20]]
else:
    p("  无失败 task。")

p("\n" + "=" * 112); p("【QC-2】各 signal 自身分布（**不含与 ΔAUROC 的任何关系**）"); p("=" * 112)
SHOW = ["domain_auc", "mmd", "psi_mean", "psi_max", "js_mean", "js_max",
        "pred_mean_shift", "entropy_mean_shift", "entropy_js", "conf_dist_js", "pred_pos_prev_shift"]
p(DF[SHOW].describe(percentiles=[.05, .25, .5, .75, .95]).round(4).to_string())

p("\n" + "=" * 112); p("【QC-3】异常检查"); p("=" * 112)
for c in SHOW:
    v = DF[c].dropna()
    if not len(v): continue
    tag = []
    if v.nunique() == 1: tag.append("**常数**")
    if not np.isfinite(v).all(): tag.append("**非有限值**")
    if v.abs().max() > 1e6: tag.append("**爆值**")
    if tag: p(f"  {c}: {tag}")
p(f"  domain_auc ≈ 1.0 的 task 数（>0.99）= {int((DF['domain_auc']>0.99).sum())}")
p(f"  domain_auc ≈ 0.5 的 task 数（0.49–0.51）= {int(DF['domain_auc'].between(0.49,0.51).sum())}")
p(f"  domain_auc 的 CI 宽度中位 = {(DF['domain_auc_hi']-DF['domain_auc_lo']).median():.4f}")

p("\n" + "=" * 112); p("【QC-4】是否受 target N 控制（无标签层面）"); p("=" * 112)
for c in ["mmd", "psi_mean", "js_mean", "domain_auc"]:
    r = DF[["n_target", c]].dropna()
    p(f"  corr(n_target, {c:<18}) = {r['n_target'].corr(r[c]):+.3f}   Spearman = "
      f"{r['n_target'].corr(r[c], method='spearman'):+.3f}")

p("\n" + "=" * 112); p("【QC-5】source→target 方向的无标签结构差异"); p("=" * 112)
DF["dir"] = DF["source_site"] + "→" + DF["target_site"]
p(DF.groupby("dir")[["domain_auc", "mmd", "psi_mean", "js_max", "entropy_js"]].median().round(4).to_string())

p("\n" + "=" * 112); p("【QC-6】field QC"); p("=" * 112)
p(f"  psi_ 列数 = {sum(1 for c in DF.columns if c.startswith('psi_'))}；"
  f"js_ 列数 = {sum(1 for c in DF.columns if c.startswith('js_'))}（应为 15 特征 + mean/max = 17）")
p(f"  行数 = {len(DF)}；task_id 唯一 = {DF['task_id'].is_unique}")
p("\n  **本轮未计算 SUDO / meta-model，未计算任何与 ΔAUROC 的相关性，未做指标择优。**")
p(f"\n完成，耗时 {time.time()-t0:.0f}s")
rep.close(); print("done", flush=True)
