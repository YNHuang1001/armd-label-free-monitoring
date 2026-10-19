# -*- coding: utf-8 -*-
r"""Protocol v1.2 §11″ 步骤 3-4 —— L3 计算（blinded，不读 target Y）

L3-A  SUDO-AURCC  —— reliability / failure-ranking proxy（不产出 Ahat_AUROC_target）
L3-B  IOW / outcome-model / DR —— formal transported-AUC estimators（产出 Ahat_AUROC_target）
      DR = primary；IOW 与 outcome-model = prespecified sensitivity

================================================================================
冻结公式（写死于脚本头部；此块在运行前固定）。★ CS-05 已更正 [OM]/[DR] 分母，见 CS05_amendment.txt
--------------------------------------------------------------------------------
记号：source train (X_s, y_s)；source internal test (X_e, y_e)；target covariates X_t
      source model 的 target 预测 pt = f(X_t)

[IOW]  e(x) = P(target | x)，由 domain classifier（logistic，source+target X）给出
       w(x) = e/(1−e)，对 source **internal test** 加权
       Ahat_IOW = Σ_{i:y=1} Σ_{j:y=0} w_i w_j 1[s_i > s_j] / [ (Σ_{y=1} w)(Σ_{y=0} w) ]

[OM]   m(x) = P(Y=1|x)，outcome model 在 source **train** 上拟合（logistic）
       Ahat_OM  = Σ_{i≠j∈T} m_i (1−m_j) 1[s_i > s_j] / Σ_{i≠j∈T} m_i (1−m_j)    [Li 2023 eq.(11)]

[DR]   Li et al. 2023 Biometrics eq.(14)+(15) 的 published estimator（分子与分母各自双重稳健）
       N_DR = Σ_{i≠j∈S_test} w_i w_j [ y_i(1−y_j) − m_i(1−m_j) ] 1[s_i>s_j]
              + Σ_{i≠j∈T} m_i (1−m_j) 1[s_i>s_j]
       D_DR = Σ_{i≠j∈S_test} w_i w_j [ y_i(1−y_j) − m_i(1−m_j) ] + Σ_{i≠j∈T} m_i (1−m_j)
       （分子分母同为 Σ_{i≠j}，故 eq.(16) 的 1/(n(n−1)) 相消；混合对贡献为 0）
       Ahat_DR  = N_DR / D_DR
       一致性：outcome model 或 density ratio 任一正确即可（双重稳健）。

[SUDO-AURCC]  §5.3.1 冻结实现
       target pt 按 10 个 score-quantile bins 分层；bin pseudo-label = (bin 内 pt 均值 ≥ 0.5)
       每 bin 每次抽该 bin 内 50% 的 target（上限 200），配等量、真实类别相反 (1−pseudo) 的 source-train
       auxiliary classifier = logistic；在 **source validation** 上评价（正类 = 真阳性，负类 = 真阴性）
       signed SUDO = 2·(AUC_aux − 0.5)；每 pseudo-label × bin 重复 5 次
       AURCC：bin 按 |bin均值 − 0.5| 降序累积，completeness = 累积 target 占比，
              reliability = 累积 |SUDO| 均值，AURCC = 该曲线下面积（梯形法）
================================================================================
输出 p6c_L3_blinded.csv / report_p6c_L3.txt
"""
import os, time, warnings
import numpy as np, pandas as pd
warnings.filterwarnings("ignore")
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score

BASE = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))); OUT = os.path.dirname(os.path.abspath(__file__))
t0 = time.time()
SEED = 20261005
NBIN, SUDO_REP, SUDO_CAP, SUDO_FRAC = 10, 5, 200, 0.5
NO_CLIP = True                       # 不截尾、不投影；异常值原样 flag

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

rep = open(os.path.join(OUT, "report_p6c_L3.txt"), "w", encoding="utf-8")
def p(*a):
    line = " ".join(str(x) for x in a); print(line, flush=True); rep.write(line + "\n"); rep.flush()

def wconcordance(score, a, b):
    """Σ_{i≠j} a_i b_j [1[score_i>score_j] + 0.5·1[score_i=score_j]]，按唯一 score 分组，O(n log n)。★ 并列记半分（Hanley–McNeil 惯例；Li 2023 eq.(15) 同）。"""
    o = np.argsort(score, kind="mergesort")
    s, aa, bb = score[o], a[o], b[o]
    uniq, inv = np.unique(s, return_inverse=True)
    A = np.bincount(inv, weights=aa); B = np.bincount(inv, weights=bb)
    cumB = np.concatenate([[0.0], np.cumsum(B)[:-1]])     # 严格小于当前 score 的 B 累积
    return float(np.sum(A * cumB) + 0.5 * (np.sum(A * B) - np.sum(aa * bb)))

p("=" * 112); p("Protocol v1.2 §11″ 步骤 3-4 —— L3 计算（blinded，不读 target Y）"); p("=" * 112)

# ---------------- 重建 STORE 与 corrected models（与 p6b 完全一致） ----------------
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

FEATS = ["age", "set_IP", "set_OP", "set_ER", "sp_URINE", "sp_BLOOD",
         "prior_any", "prior_cnt", "t_15-30", "t_31-90", "t_91-180", "t_181-365",
         "past_S", "past_I", "past_R"]
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

CM = pd.read_csv(os.path.join(OUT, "p6b_corrected_universe.csv"))
TASKS = CM[CM["eligible"]].copy()
KEYS = sorted(set(zip(TASKS.organism, TASKS.drug, TASKS.source_site)))

MC = {}
for sym, drug, src in KEYS:
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
        af = float(np.nanmedian(Dall[i_tr, 0])) if np.isfinite(Dall[i_tr, 0]).any() else 0.0
        Dall[:, 0] = np.where(np.isfinite(Dall[:, 0]), Dall[:, 0], af)
        Dtr, Dva, Dte = Dall[i_tr], Dall[i_va], Dall[i_te]
        sc = StandardScaler().fit(Dtr)
        mdl = LogisticRegression(max_iter=2000, C=1.0, random_state=SEED).fit(sc.transform(Dtr), y[i_tr])
        MC[(sym, drug, src)] = dict(Dtr=Dtr, Dva=Dva, Dte=Dte, sc=sc, mdl=mdl,
                                    ytr=y[i_tr], yva=y[i_va], yte=y[i_te], af=af)
    except Exception:
        continue
p(f"\n  模型就绪 {len(MC)} 个  [{time.time()-t0:.0f}s]")

# ---------------- 主循环 ----------------
rows = []
for _, t in TASKS.iterrows():
    sym, drug, src, tgt = t["organism"], t["drug"], t["source_site"], t["target_site"]
    rec = {"task_id": t["task_id"], "organism": sym, "drug": drug,
           "source_site": src, "target_site": tgt}
    reasons = []
    try:
        M = MC[(sym, drug, src)]
        Xt, Ht = prep(STORE[(tgt, sym)], drug)
        if len(Xt) < 300: raise ValueError("target<300")
        Dt = design(Xt, Ht).values.copy()
        Dt[:, 0] = np.where(np.isfinite(Dt[:, 0]), Dt[:, 0], M["af"])
        Ds_std = M["sc"].transform(np.vstack([M["Dtr"], M["Dva"], M["Dte"]]))
        Ds_tr_std = Ds_std[:len(M["ytr"])]          # source TRAIN 部分
        Dt_std = M["sc"].transform(Dt)
        rng = np.random.default_rng(SEED + 11)
        i_tr = np.arange(len(M["ytr"])); i_va = np.arange(len(M["yva"])) + len(i_tr)
        i_te = np.arange(len(M["yte"])) + len(i_tr) + len(i_va)
        rec.update({"n_source": len(i_tr) + len(i_va) + len(i_te), "n_target": len(Xt),
                    "n_src_train": len(i_tr), "n_src_val": len(i_va), "n_src_test": len(i_te)})
        # 目标模型对 target 的预测
        pt = M["mdl"].predict_proba(M["sc"].transform(Dt))[:, 1]
        sc_t = M["mdl"].decision_function(M["sc"].transform(Dt))

        # ============ L3-B: transported AUC ============
        # domain model e(x)（用于 density ratio）
        Xd = np.vstack([Ds_std, Dt_std]); yd = np.r_[np.zeros(len(Ds_std)), np.ones(len(Dt_std))]
        dm = LogisticRegression(max_iter=2000, class_weight="balanced", random_state=SEED)
        dm.fit(Xd, yd)
        e_s = np.clip(dm.predict_proba(Ds_std)[:, 1], 1e-6, 1 - 1e-6)
        w_all = e_s / (1 - e_s)
        w_te = w_all[i_te]
        rec["w_P50"] = float(np.percentile(w_all, 50)); rec["w_P95"] = float(np.percentile(w_all, 95))
        rec["w_P99"] = float(np.percentile(w_all, 99)); rec["w_max"] = float(w_all.max())
        rec["w_ESS"] = float((w_te.sum() ** 2) / (w_te ** 2).sum())
        # IOW：source internal test 上加权
        s_e = M["mdl"].decision_function(M["sc"].transform(M["Dte"]))
        y_e = M["yte"]
        num = wconcordance(s_e, w_te * (y_e == 1), w_te * (y_e == 0))
        den = (w_te * (y_e == 1)).sum() * (w_te * (y_e == 0)).sum()
        rec["AUROC_IOW"] = float(num / den) if den > 0 else np.nan
        # outcome model m(x)（source train 上拟合）
        om = LogisticRegression(max_iter=2000, random_state=SEED).fit(M["sc"].transform(M["Dtr"]), M["ytr"])
        m_t = np.clip(om.predict_proba(M["sc"].transform(Dt))[:, 1], 1e-9, 1 - 1e-9)
        m_e = np.clip(om.predict_proba(M["sc"].transform(M["Dte"]))[:, 1], 1e-9, 1 - 1e-9)
        # OM：target 侧
        numo = wconcordance(sc_t, m_t, 1 - m_t); deno = m_t.sum() * (1 - m_t).sum() - (m_t * (1 - m_t)).sum()
        rec["AUROC_OM"] = float(numo / deno) if deno > 0 else np.nan
        # DR：augmented IPW
        yt1 = (y_e == 1).astype(float); yt0 = (y_e == 0).astype(float); w1 = w_te * yt1; w0 = w_te * yt0; p1 = w_te * m_e; p0 = w_te * (1 - m_e)
        term_s = wconcordance(s_e, w1, w0) - wconcordance(s_e, p1, p0)
        term_t = wconcordance(sc_t, m_t, 1 - m_t)
        n_s, n_t = len(y_e), len(Dt)
        N_DR = term_s + term_t
        d_s = (w1.sum() * w0.sum() - (w1 * w0).sum()) - (p1.sum() * p0.sum() - (p1 * p0).sum())
        d_t = m_t.sum() * (1 - m_t).sum() - (m_t * (1 - m_t)).sum()
        rec["AUROC_DR"] = float(N_DR / (d_s + d_t)) if (d_s + d_t) > 0 else np.nan

        # ============ L3-A: SUDO-AURCC ============
        qs = np.quantile(pt, np.linspace(0, 1, NBIN + 1)); qs[0], qs[-1] = -np.inf, np.inf
        bid = np.digitize(pt, qs[1:-1])
        src_by_class = {0: np.where(M["ytr"] == 0)[0], 1: np.where(M["ytr"] == 1)[0]}
        sudo_bins, frac, conf = [], [], []
        for b in range(NBIN):
            idx = np.where(bid == b)[0]
            if len(idx) < 10: continue
            yb = 1 if pt[idx].mean() >= 0.5 else 0
            pool = src_by_class[1 - yb]
            if len(pool) < 10: continue
            aucs = []
            for r in range(SUDO_REP):
                nb = min(int(SUDO_FRAC * len(idx)), SUDO_CAP, len(idx), len(pool))
                ib = rng.choice(idx, nb, replace=False)
                sb = rng.choice(pool, nb, replace=False)
                # 正类 = **target** 观测（该 bin，伪标签 yb）；负类 = source-train 中真实类别相反者
                Xa = np.vstack([Dt_std[ib], Ds_tr_std[sb]])
                ya = np.r_[np.ones(nb), np.zeros(nb)]
                clf = LogisticRegression(max_iter=2000, random_state=SEED + r).fit(Xa, ya)
                # 在 source validation 上评价：正类 = 真实标签 == 伪标签 yb 者
                sv = M["sc"].transform(M["Dva"])
                yv = (M["yva"] == yb).astype(int)
                if len(np.unique(yv)) < 2: continue
                aucs.append(roc_auc_score(yv, clf.decision_function(sv)))
            if len(aucs) < 3: continue
            sudo_bins.append(2 * (float(np.mean(aucs)) - 0.5))
            frac.append(len(idx) / len(pt)); conf.append(abs(pt[idx].mean() - 0.5))
        if sudo_bins:
            order = np.argsort(-np.array(conf))
            fr = np.array(frac)[order]; rb = np.abs(np.array(sudo_bins))[order]
            comp = np.cumsum(fr) / np.sum(fr)
            rel = np.cumsum(rb * fr) / np.cumsum(fr)
            rec["SUDO_AURCC"] = float(np.trapezoid(rel, comp))
            rec["SUDO_nbin"] = len(sudo_bins)
            rec["SUDO_signed_min"] = float(np.min(sudo_bins))
            rec["SUDO_signed_max"] = float(np.max(sudo_bins))
            rec["SUDO_absmean"] = float(np.mean(np.abs(sudo_bins)))
        else:
            raise ValueError("SUDO: 无有效 bin")
        rec["status"] = "OK"
    except Exception as e:
        rec["status"] = "FAILED"; rec["failure_reason"] = f"{type(e).__name__}: {e}"
        reasons.append(str(e))
    rows.append(rec)
    # ★ 2026-10-06：分母原写死 202（修正前的 task 数），改为 len(TASKS) 现算。
    if len(rows) % 40 == 0: p(f"    ... {len(rows)}/{len(TASKS)}  [{time.time()-t0:.0f}s]")

DF = pd.DataFrame(rows)
DF.to_csv(os.path.join(OUT, "p6c_L3_blinded.csv"), index=False)
# ★ 2026-10-06：原写 "(corrected universe, 分母固定 202)" —— "分母固定 202" 是**错的说法**：
#   分母是 corrected universe 的规模，它随修正变化（本期 186）。改为如实报告来源。
p(f"\n[写出] p6c_L3_blinded.csv  行数 = {len(DF)}  "
  f"(corrected universe = p6b_corrected_universe.csv 中 eligible 的 task 数；本期 {len(TASKS)})")
p(f"  status OK = {int((DF['status']=='OK').sum())} / {len(DF)}")
if (DF["status"] != "OK").any():
    p("  失败：" + str(DF.loc[DF.status != "OK", ["task_id", "failure_reason"]].values.tolist()[:10]))

p("\n" + "=" * 112); p("【L3-A QC】SUDO-AURCC"); p("=" * 112)
s = DF["SUDO_AURCC"].dropna()
p(f"  覆盖率 = {len(s)}/{len(DF)}")
p(s.describe(percentiles=[.05, .25, .5, .75, .95]).round(4).to_string())
p(f"  有效 bin 数分布: {dict(DF['SUDO_nbin'].value_counts().sort_index())}")
p(f"  signed SUDO 极值: min={DF['SUDO_signed_min'].min():.3f}  max={DF['SUDO_signed_max'].max():.3f}")
p(f"  |signed SUDO|>1 的 task 数 = {int((DF['SUDO_signed_min'].abs()>1).sum() + (DF['SUDO_signed_max'].abs()>1).sum())}")

p("\n" + "=" * 112); p("【L3-B QC】IOW / OM / DR"); p("=" * 112)
for c in ["AUROC_IOW", "AUROC_OM", "AUROC_DR"]:
    v = DF[c].dropna()
    fin = np.isfinite(v)
    p(f"  {c:<12} 覆盖率={len(v)}/{len(DF)}  有限={int(fin.sum())}  "
      f"range=[{v.min():.4f}, {v.max():.4f}]  <0: {int((v<0).sum())}  >1: {int((v>1).sum())}")
p(f"\n  权重: P50={DF['w_P50'].median():.3f}  P95={DF['w_P95'].median():.3f}  "
  f"P99={DF['w_P99'].median():.3f}  max={DF['w_max'].max():.1f}")
p(f"  权重 ESS: 中位={DF['w_ESS'].median():.1f}  最小={DF['w_ESS'].min():.1f}  "
  f"占 n_src_test 比例中位={(DF['w_ESS']/DF['n_src_test']).median()*100:.1f}%")
p(f"  未做任何截尾（NO_CLIP={NO_CLIP}）；未对 A-hat(UROC) 做 [0,1] 投影")

p("\n" + "=" * 112); p("【无标签 QC】与 n_target / source→target 方向的结构（不碰真实性能）"); p("=" * 112)
for c in ["SUDO_AURCC", "AUROC_IOW", "AUROC_OM", "AUROC_DR"]:
    r = DF[["n_target", c]].dropna()
    p(f"  rho(n_target, {c:<12}) = {r['n_target'].corr(r[c], method='spearman'):+.3f}")
DF["dir"] = DF["source_site"] + "→" + DF["target_site"]
p("\n" + DF.groupby("dir")[["SUDO_AURCC", "AUROC_IOW", "AUROC_OM", "AUROC_DR"]].median().round(4).to_string())
p("\n  **SUDO 与 transported-AUC 角色不同，本轮不做互相择优；主次已在揭盲前冻结（DR = primary）。**")
p("\n  **未读取 target Y；未计算任何 L3 指标与 ΔAUROC / harmful 的关系。**")
p(f"\n完成，耗时 {time.time()-t0:.0f}s")
rep.close(); print("done", flush=True)
