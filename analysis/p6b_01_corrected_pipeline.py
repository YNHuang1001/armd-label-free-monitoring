# -*- coding: utf-8 -*-
r"""Protocol v1.1 §11′ 步骤 1-4 —— corrected blinded pipeline

1. 从 228 candidate tasks 重建 **protocol-compliant** source models
   （age 中位数 / 编码 / 标准化 全部只来自 source train）
2. **只用 source Y** 计算 source internal AUROC
3. 按**原 0.65 阈值**决定 corrected eligible universe（不强行保留 202）
4. 在 corrected universe 上重算 §5.1 / §5.2（MMD 已按 v1.1 修正）

**全程不读取 target Y。**

MMD v1.1 修正：输入 = source-train scaler 标准化；bandwidth = source train 两两距离中位数；
estimator = unbiased U-statistic（负值原样报告，不 clip）。
domain classifier（§5.1）同样使用 source-train 标准化，与 MMD 输入保持一致。

输出 p6b_candidate_manifest_blinded.csv / p6b_corrected_universe.csv
     p6b_label_free_signals_blinded.csv / report_p6b_corrected.txt
"""
import os, time, warnings
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
MMD_M, MMD_K, NBIN = 1000, 5, 10
AUROC_MIN = 0.65                      # 阈值不变

# ---------- 步骤 1 前置：228 candidate manifest（只取 key，无性能字段） ----------
cand = pd.read_csv(os.path.join(OUT, "p5_tasks_primary.csv"),
                   usecols=lambda c: c in ("organism", "drug", "source", "target"))
assert len(cand) == 228, f"期望 228 candidate，实得 {len(cand)}"
CM = pd.DataFrame({
    "task_id": cand["organism"] + "_" + cand["drug"] + "_" + cand["source"] + "2" + cand["target"],
    "organism": cand["organism"], "drug": cand["drug"],
    "source_site": cand["source"], "target_site": cand["target"]})
assert CM["task_id"].is_unique
CM.to_csv(os.path.join(OUT, "p6b_candidate_manifest_blinded.csv"), index=False)

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

rep = open(os.path.join(OUT, "report_p6b_corrected.txt"), "w", encoding="utf-8")
def p(*a):
    line = " ".join(str(x) for x in a); print(line, flush=True); rep.write(line + "\n"); rep.flush()

p("=" * 112); p("Protocol v1.1 §11′ 步骤 1-4 —— corrected blinded pipeline（不读 target Y）"); p("=" * 112)
p(f"candidate = 228；eligibility 阈值 = source AUROC >= {AUROC_MIN}（不变）")

# ---------- 数据（与 p5/p6 同源，episode 用 anchor） ----------
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

def build_hist(X, drug):
    n = len(X); pr_any = np.zeros(n); pr_cnt = np.zeros(n)
    tbin = np.full(n, "none", dtype=object); past = np.full(n, "none", dtype=object)
    c0 = X["cdate"].values
    for _, idx in X.groupby(["anon_id", "spec"], sort=False).indices.items():
        idx = np.sort(np.asarray(idx))
        for k in range(len(idx)):
            i = idx[k]
            lo = c0[i] - np.timedelta64(365, "D"); hi = c0[i] - np.timedelta64(15, "D")
            pr = idx[:k]; sel = pr[(c0[pr] > lo) & (c0[pr] <= hi)]
            if len(sel) == 0: continue
            pr_any[i] = 1; pr_cnt[i] = len(sel)
            dts = (c0[i] - c0[sel]).astype("timedelta64[D]").astype(int); mn = dts.min()
            tbin[i] = ("15-30" if mn <= 30 else "31-90" if mn <= 90 else "91-180" if mn <= 180 else "181-365")
            j = sel[np.argmin(dts)]; v = X[drug].values[j] if drug in X.columns else None
            past[i] = v if v in ("S", "I", "R") else "none"
    return pd.DataFrame({"prior_any": pr_any, "prior_cnt": pr_cnt, "tbin": tbin, "past_ast": past})

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

def psijs(rr, aa):
    r = rr.astype(float); a = aa.astype(float)
    rp = r / max(r.sum(), 1); ap = a / max(a.sum(), 1)
    eps = 1e-6
    rpc = np.clip(rp, eps, 1); apc = np.clip(ap, eps, 1)
    psi = float(np.sum((apc - rpc) * np.log(apc / rpc)))
    M = 0.5 * (rp + ap)
    def kl(x, y):
        m = x > 0; return float(np.sum(x[m] * np.log2(x[m] / np.maximum(y[m], 1e-12))))
    return psi, 0.5 * kl(rp, M) + 0.5 * kl(ap, M)

def feat_psijs(ref, act, cont):
    ref = np.asarray(ref); act = np.asarray(act)
    if cont:
        qs = np.unique(np.quantile(ref, np.linspace(0, 1, NBIN + 1)))
        if len(qs) < 3: return np.nan, np.nan
        e = qs.astype(float).copy(); e[0] = -np.inf; e[-1] = np.inf
        r = np.histogram(ref, bins=e)[0]; a = np.histogram(act, bins=e)[0]
    else:
        cats = np.unique(np.r_[np.unique(ref), np.unique(act)])
        r = np.bincount(np.searchsorted(cats, ref), minlength=len(cats))
        a = np.bincount(np.searchsorted(cats, act), minlength=len(cats))
    if r.sum() == 0 or a.sum() == 0: return np.nan, np.nan
    return psijs(r, a)

def mmd_rbf(A, B, bw2, rng):
    """unbiased U-statistic MMD²；负值允许，不 clip。"""
    vals = []
    for _ in range(MMD_K):
        ia = rng.choice(len(A), min(MMD_M, len(A)), replace=False)
        ib = rng.choice(len(B), min(MMD_M, len(B)), replace=False)
        a, b = A[ia], B[ib]
        Kaa = np.exp(-pdist(a, "sqeuclidean") / bw2)
        Kbb = np.exp(-pdist(b, "sqeuclidean") / bw2)
        Kab = np.exp(-cdist(a, b, "sqeuclidean") / bw2)
        n, m = len(a), len(b)
        vals.append(2 * Kaa.sum() / (n * (n - 1)) + 2 * Kbb.sum() / (m * (m - 1)) - 2 * Kab.mean())
    return float(np.mean(vals)), float(np.std(vals))

def js_cont(ref, act):
    qs = np.unique(np.quantile(ref, np.linspace(0, 1, NBIN + 1)))
    if len(qs) < 3: return np.nan
    e = qs.astype(float).copy(); e[0] = -np.inf; e[-1] = np.inf
    return psijs(np.histogram(ref, bins=e)[0], np.histogram(act, bins=e)[0])[1]

def ent(pv):
    pv = np.clip(pv, 1e-12, 1 - 1e-12); return -(pv * np.log(pv) + (1 - pv) * np.log(1 - pv))

# ---------- 步骤 1-2：per (source, organism, drug) 建 protocol-compliant model ----------
MC, UNIV = {}, []
for _, r in CM.iterrows():
    sym, drug, src, tgt = r["organism"], r["drug"], r["source_site"], r["target_site"]
    key = (src, sym, drug)
    if key in MC: continue
    try:
        Xs, Hs = prep(STORE[(src, sym)], drug)
        if len(Xs) < 500: continue
        D = design(Xs, Hs); y = (Xs[drug].values == "R").astype(int)
        pats = Xs["anon_id"].values; upat = pd.unique(pats)
        rng = np.random.default_rng(SEED); perm = rng.permutation(len(upat))
        ntr, nva = int(0.70 * len(upat)), int(0.15 * len(upat))
        s_tr = set(upat[perm[:ntr]]); s_va = set(upat[perm[ntr:ntr + nva]])
        i_tr = np.array([i for i, a in enumerate(pats) if a in s_tr])
        i_va = np.array([i for i, a in enumerate(pats) if a in s_va])
        i_te = np.array([i for i, a in enumerate(pats) if a not in s_tr and a not in s_va])
        if min(len(i_tr), len(i_va), len(i_te)) < 50: continue
        Dall = D.values.copy()
        age_miss = float(np.mean(~np.isfinite(Dall[:, 0])) * 100)
        age_fill = float(np.nanmedian(Dall[i_tr, 0])) if np.isfinite(Dall[i_tr, 0]).any() else 0.0
        Dall[:, 0] = np.where(np.isfinite(Dall[:, 0]), Dall[:, 0], age_fill)   # 只用 source train 中位数
        Dtr, Dva, Dte = Dall[i_tr], Dall[i_va], Dall[i_te]
        sc = StandardScaler().fit(Dtr)                    # 只 fit 在 source train
        mdl = LogisticRegression(max_iter=2000, C=1.0, random_state=SEED)
        mdl.fit(sc.transform(Dtr), y[i_tr])
        yte = y[i_te]; ps = mdl.predict_proba(sc.transform(Dte))[:, 1]
        src_auroc = float(roc_auc_score(yte, ps)) if len(np.unique(yte)) > 1 else np.nan
        # bandwidth：只从 source train 估计（固定子样本）
        rbw = np.random.default_rng(SEED + 3)
        ib = rbw.choice(len(Dtr), min(MMD_M, len(Dtr)), replace=False)
        Dtr_std = sc.transform(Dtr)
        bw2 = float(np.median(pdist(Dtr_std[ib], "sqeuclidean")))
        if not np.isfinite(bw2) or bw2 <= 0: bw2 = 1.0
        MC[key] = dict(sc=sc, mdl=mdl, Ds=Dall, Dtr=Dtr, Dte=Dte, ps=ps, yte=yte,
                       src_auroc=src_auroc, bw2=bw2, age_miss=age_miss, age_fill=age_fill,
                       n_src=len(Xs), n_tr=len(i_tr), n_va=len(i_va), n_te=len(i_te))
        UNIV.append(dict(organism=sym, drug=drug, source_site=src, src_auroc=src_auroc))
    except Exception as e:
        p(f"  [warn] {key}: {type(e).__name__}: {e}")
p(f"\n  模型建成 {len(MC)} 个 (source, organism, drug)  [{time.time()-t0:.0f}s]")

U = pd.DataFrame(UNIV)
U["eligible"] = U["src_auroc"] >= AUROC_MIN
p("\n" + "=" * 112); p("【步骤 3】corrected eligible universe（阈值 0.65 不变）"); p("=" * 112)
p(f"  candidate transfer task = {len(CM)}")
elig_keys = set(zip(U.loc[U.eligible, "organism"], U.loc[U.eligible, "drug"], U.loc[U.eligible, "source_site"]))
CM["eligible"] = [ (a, b, c) in elig_keys for a, b, c in zip(CM.organism, CM.drug, CM.source_site) ]
p(f"  corrected eligible task = {int(CM['eligible'].sum())} / {len(CM)}")
p(f"  被排除 = {int((~CM['eligible']).sum())}")
if (U.src_auroc < AUROC_MIN).any():
    p(f"  被排除模型的 source AUROC 范围 = {U.loc[U.src_auroc<AUROC_MIN,'src_auroc'].min():.3f}"
      f" ~ {U.loc[U.src_auroc<AUROC_MIN,'src_auroc'].max():.3f}")
p(f"  覆盖 organism–drug 分组 = {U[U.eligible].groupby(['organism','drug']).ngroups}")
CM.to_csv(os.path.join(OUT, "p6b_corrected_universe.csv"), index=False)

# ---------- 步骤 4：corrected universe 上重算 §5.1/§5.2 ----------
TASKS = CM[CM["eligible"]]
p("\n" + "=" * 112); p("【步骤 4】corrected universe 上重算 §5.1 / §5.2"); p("=" * 112)
rows, fails = [], []
for _, t in TASKS.iterrows():
    sym, drug, src, tgt = t["organism"], t["drug"], t["source_site"], t["target_site"]
    M = MC.get((src, sym, drug))
    rec = {"task_id": t["task_id"], "organism": sym, "drug": drug,
           "source_site": src, "target_site": tgt, "source_auroc_internal": M["src_auroc"]}
    try:
        Xt, Ht = prep(STORE[(tgt, sym)], drug)
        if len(Xt) < 300: raise ValueError("target < 300")
        Dt = design(Xt, Ht).values.copy()
        Dt[:, 0] = np.where(np.isfinite(Dt[:, 0]), Dt[:, 0], M["age_fill"])   # source-train 填补值
        Ds_std = M["sc"].transform(M["Ds"]); Dt_std = M["sc"].transform(Dt)
        rec.update({"n_source": M["n_src"], "n_target": len(Xt),
                    "n_src_train": M["n_tr"], "n_src_val": M["n_va"], "n_src_test": M["n_te"]})
        # §5.1a domain classifier（source-train 标准化输入；OOF；class_weight balanced）
        Xd = np.vstack([Ds_std, Dt_std]); yd = np.r_[np.zeros(len(Ds_std)), np.ones(len(Dt_std))]
        oof = np.zeros(len(yd))
        for tr, te in StratifiedKFold(5, shuffle=True, random_state=SEED).split(Xd, yd):
            dm = LogisticRegression(max_iter=2000, class_weight="balanced", random_state=SEED)
            dm.fit(Xd[tr], yd[tr]); oof[te] = dm.predict_proba(Xd[te])[:, 1]
        rec["domain_auc"] = float(roc_auc_score(yd, oof))
        rb = np.random.default_rng(SEED + 1); bs = []
        for _ in range(200):
            s = rb.integers(0, len(yd), len(yd))
            if len(np.unique(yd[s])) < 2: continue
            bs.append(roc_auc_score(yd[s], oof[s]))
        rec["domain_auc_lo"], rec["domain_auc_hi"] = np.percentile(bs, [2.5, 97.5])
        # §5.1b MMD（v1.1 修正：标准化输入 + source-train bandwidth）
        mv, ms = mmd_rbf(Ds_std, Dt_std, M["bw2"], np.random.default_rng(SEED + 2))
        rec["mmd"], rec["mmd_sd"] = mv, ms
        # §5.1c PSI / JS（source TRAIN 参考）
        ref = M["Dtr"]; psis = {}; jss = {}
        for j, f in enumerate(FEATS):
            pv, jv = feat_psijs(ref[:, j], Dt[:, j], f in CONT)
            psis[f] = pv; jss[f] = jv; rec[f"psi_{f}"] = pv; rec[f"js_{f}"] = jv
        pa = np.array(list(psis.values()), float); ja = np.array(list(jss.values()), float)
        rec["psi_mean"] = float(np.nanmean(pa)); rec["psi_max"] = float(np.nanmax(pa))
        rec["js_mean"] = float(np.nanmean(ja)); rec["js_max"] = float(np.nanmax(ja))
        # §5.2 prediction-space（source internal TEST 为参考）
        pt = M["mdl"].predict_proba(M["sc"].transform(Dt))[:, 1]
        rec["pred_mean_shift"] = float(pt.mean() - M["ps"].mean())
        rec["entropy_mean_shift"] = float(ent(pt).mean() - ent(M["ps"]).mean())
        rec["entropy_js"] = js_cont(ent(M["ps"]), ent(pt))
        rec["conf_dist_js"] = js_cont(np.maximum(M["ps"], 1 - M["ps"]), np.maximum(pt, 1 - pt))
        rec["pred_pos_prev_shift"] = float((pt > 0.5).mean() - (M["ps"] > 0.5).mean())
        rec["status"] = "OK"
    except Exception as e:
        rec["status"] = "FAILED"; rec["failure_reason"] = f"{type(e).__name__}: {e}"
        fails.append((t["task_id"], rec["failure_reason"]))
    rows.append(rec)

DF = pd.DataFrame(rows)
DF.to_csv(os.path.join(OUT, "p6b_label_free_signals_blinded.csv"), index=False)
p(f"\n[写出] p6b_label_free_signals_blinded.csv  行数 = {len(DF)}  (corrected universe)")
p(f"  status OK = {int((DF['status']=='OK').sum())} / {len(DF)}")
if fails: p(f"  失败 {len(fails)}：{fails[:5]}")

SHOW = ["domain_auc", "mmd", "psi_mean", "psi_max", "js_mean", "js_max",
        "pred_mean_shift", "entropy_mean_shift", "entropy_js", "conf_dist_js", "pred_pos_prev_shift"]
p("\n【QC】corrected signal 分布")
p(DF[SHOW].describe(percentiles=[.05, .5, .95]).round(4).to_string())
p(f"\n  domain_auc > 0.99 的 task = {int((DF['domain_auc']>0.99).sum())}")
p(f"  mmd 全为负? {bool((DF.mmd<0).all())}   范围 {DF.mmd.min():.4f} ~ {DF.mmd.max():.4f}")
p("\n【QC】n_target 与各 metric（cross-task association，**不作机制结论**）")
for c in ["mmd", "psi_mean", "js_mean", "domain_auc"]:
    rr = DF[["n_target", c]].dropna()
    p(f"  rho(n_target, {c:<12}) = {rr['n_target'].corr(rr[c], method='spearman'):+.3f}")
p("\n  **未读取 target Y；未计算任何 signal 与 ΔAUROC 的关系。**")
p(f"\n完成，耗时 {time.time()-t0:.0f}s")
rep.close(); print("done", flush=True)
