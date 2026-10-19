# -*- coding: utf-8 -*-
r"""Protocol v1.3 §10.1 —— Phase A: 为每个 task 预计算并缓存 §10 所需的**冻结量**

§10.1 约束（v1.3 冻结）：
  - domain classifier / membership model : 每 task 固定 5,000 source units，5 档 N x 200 repeats 共用
  - MMD                                    : source reference 固定 m=1000（5 次 draw 固定），Kaa 预计算
  - PSI/JS、prediction-space、outcome-model: full source reference（fixed）
  - SUDO                                   : source train/val 原样，只变 target subsample
  - 只反复抽 target X；target Y 不参与（本阶段完全不写任何 outcome 字段）

预计算内容（都在 task 内固定，与 N / repeat 无关）：
  * Ds_ref_std  : 5,000 行固定 source reference（标准化后），seed = f(task_id)
  * mmd_src     : 5 次 draw 各 1000 行（标准化后，固定）；Kaa = 每次 draw 的 off-diag 均值
  * PSI/JS      : 每个 feature 的 bin 边界与 source-train 参考直方图 + target 每行的 bin 索引
                  （bin 用 source-train ref 与 **full target** 共同确定，与 p6b 完全一致）
  * 其余        : source train/val/test 标准化矩阵与标签、source test 预测 ps、
                  outcome model 在 source test / 全 target 上的预测、target 标准化矩阵与预测 pt

输出：p10_cache/<task_id>.pkl （每个 task 一个文件）
"""
import os, time, pickle, hashlib, warnings
import numpy as np, pandas as pd
warnings.filterwarnings("ignore")
from scipy.spatial.distance import pdist, cdist
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score

BASE = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))); OUT = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(OUT, "p10_cache")
os.makedirs(CACHE, exist_ok=True)
SEED = 20261005
MMD_M, MMD_K, NBIN = 1000, 5, 10
SRC_REF_N = 5000            # §10.1 固定 source reference 大小


def _prov_signature():
    """缓存**产出代码签名** = sha256(p10a 源码 + p10b 源码)。

    ★ 2026-10-06 新增。背景：2026-10-06 发现 `p10_cache/*.pkl` 里**没有任何来源标记**，
      而 p10b 直接 `pickle.load` —— 于是补丁 A3 改了 `Kaa` 的含义、CS-01 v2 改了药物身份
      之后，旧缓存仍会被静默当成新缓存使用。
    ★ 刻意用**源码整体**做签名：任何一方被编辑（哪怕只是注释）签名都会变。
      这是**保守方向** —— 宁可多算一次，也不可复用来源不明的缓存。
      p10b 侧有同名同实现的副本；两处必须一致。"""
    h = hashlib.sha256()
    here = os.path.dirname(os.path.abspath(__file__))
    for fn in ("p10a_prepare_cache.py", "p10b_subsample_curves.py"):
        with open(os.path.join(here, fn), "rb") as f:
            h.update(f.read())
    return h.hexdigest()


t0 = time.time()
rep = open(os.path.join(OUT, "report_p10a_prepare.txt"), "w", encoding="utf-8")
def p(*a):
    line = " ".join(str(x) for x in a); print(line, flush=True); rep.write(line + "\n"); rep.flush()

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
def task_seed(tid):
    """由 task_id 派生，跨进程/跨运行稳定（不用 Python 内置 hash）。"""
    return int(hashlib.md5(("p10|" + tid).encode()).hexdigest()[:8], 16)

# ---------------- STORE（与 p6b/p6c 完全一致） ----------------
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
    p(f"  [{site}] STORE ready  [{time.time()-t0:.0f}s]")

# ============================================================================
# 以下三个函数【逐字复制自 p6b_01_corrected_pipeline.py】（已核对 p6c 为同一实现）
# v1 曾自行重写 build_hist/design，导致 202/202 task 的 source model 与 p6b 不一致
#   （|d source AUROC| 均值 0.024，最大 0.095），该轮缓存与 §10 结果作废。
# 禁止再自行"等价改写"这段代码。
# ============================================================================
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
def ent(pv):
    pv = np.clip(pv, 1e-12, 1 - 1e-12); return -(pv * np.log(pv) + (1 - pv) * np.log(1 - pv))
def _bin_edges(ref, act):
    """连续特征：分位 bin（与 p6b feat_psijs 完全一致，含 act 影响）"""
    qs = np.unique(np.quantile(ref, np.linspace(0, 1, NBIN + 1)))
    if len(qs) < 3: return None
    e = qs.astype(float).copy(); e[0] = -np.inf; e[-1] = np.inf
    return e
def _cats(ref, act):
    return np.unique(np.r_[np.unique(ref), np.unique(act)])

CM = pd.read_csv(os.path.join(OUT, "p6b_corrected_universe.csv"))
TASKS = CM[CM["eligible"]].copy()
KEYS = sorted(set(zip(TASKS.organism, TASKS.drug, TASKS.source_site)))
p(f"\n  tasks = {len(TASKS)}   models needed = {len(KEYS)}")

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
        ps = mdl.predict_proba(sc.transform(Dte))[:, 1]
        rbw = np.random.default_rng(SEED + 3)
        ib = rbw.choice(len(Dtr), min(MMD_M, len(Dtr)), replace=False)
        Dtr_std = sc.transform(Dtr)
        bw2 = float(np.median(pdist(Dtr_std[ib], "sqeuclidean")))
        if not np.isfinite(bw2) or bw2 <= 0: bw2 = 1.0
        MC[(sym, drug, src)] = dict(sc=sc, mdl=mdl, Dtr=Dtr, Dva=Dva, Dte=Dte, Dall=Dall,
                                    ytr=y[i_tr], yva=y[i_va], yte=y[i_te],
                                    ps=ps, af=af, bw2=bw2, i_tr=i_tr, i_va=i_va, i_te=i_te,
                                    Dtr_std=Dtr_std)
    except Exception as e:
        p(f"  [warn] {(sym,drug,src)}: {type(e).__name__}: {e}")
p(f"  模型建成 {len(MC)} 个  [{time.time()-t0:.0f}s]")

# ============================================================================
# 【强制校验闸门】source internal AUROC 必须逐 task 复现 p6b（容差 1e-9）
# v1 的 build_hist 重写错误就是被这道闸门拦下的；不得删除。
# ============================================================================
_ref = pd.read_csv(os.path.join(OUT, "p6b_label_free_signals_blinded.csv")).set_index("task_id")
_dmax, _nbad, _miss = 0.0, 0, 0
for _, t in TASKS.iterrows():
    tid = t["task_id"]; key = (t["organism"], t["drug"], t["source_site"])
    if key not in MC or tid not in _ref.index: _miss += 1; continue
    M = MC[key]
    a_now = float(roc_auc_score(M["yte"], M["ps"]))
    d = abs(a_now - float(_ref.loc[tid, "source_auroc_internal"]))
    _dmax = max(_dmax, d)
    if d > 1e-9: _nbad += 1
p(f"\n  [校验闸门] source internal AUROC vs p6b：max|diff| = {_dmax:.3e}   不一致 task = {_nbad}   跳过 = {_miss}")
if _nbad > 0:
    p("  *** 校验失败：模型与 p6b 不一致，缓存不予写出。终止。 ***")
    rep.close(); raise SystemExit(1)
p("  [校验闸门] 通过 —— 与 p6b 逐 task 一致")

# ---------------- 逐 task 缓存 ----------------
ok, bad = 0, []
for n_done, (_, t) in enumerate(TASKS.iterrows(), 1):
    tid = t["task_id"]
    try:
        M = MC[(t["organism"], t["drug"], t["source_site"])]
        sc, mdl = M["sc"], M["mdl"]
        Ds_std = sc.transform(M["Dall"])          # 顺序 = [Dtr; Dva; Dte]
        Dtr_std, Dva_std, Dte_std = sc.transform(M["Dtr"]), sc.transform(M["Dva"]), sc.transform(M["Dte"])

        # --- 固定 5,000 source reference（domain classifier / membership model 共用）---
        r_ref = np.random.default_rng(task_seed(tid))
        n_ref = min(SRC_REF_N, len(Ds_std))
        ref_idx = r_ref.choice(len(Ds_std), n_ref, replace=False)
        Ds_ref_std = Ds_std[np.sort(ref_idx)]

        # --- 固定 MMD source reference：5 次 draw 各 1000，Kaa 预计算 ---
        mmd_src, Kaa = [], []
        for k in range(MMD_K):
            ia = r_ref.choice(len(Ds_std), min(MMD_M, len(Ds_std)), replace=False)
            a = Ds_std[ia]
            Ka = np.exp(-pdist(a, "sqeuclidean") / M["bw2"])
            Kaa.append(float(2 * Ka.sum() / (len(a) * (len(a) - 1))))
            mmd_src.append(a)

        # --- target ---
        Xt, Ht = prep(STORE[(t["target_site"], t["organism"])], t["drug"])
        if len(Xt) < 300: raise ValueError("target<300")
        Dt_raw = design(Xt, Ht).values.copy()
        Dt_raw[:, 0] = np.where(np.isfinite(Dt_raw[:, 0]), Dt_raw[:, 0], M["af"])

        # PSI/JS：bin 用 source-train ref 与 full target 共同确定（与 p6b 一致），之后固定
        psi_bins, psi_refhist, tgt_binidx = [], [], []
        for j, f in enumerate(FEATS):
            ref = M["Dtr"][:, j]; act = Dt_raw[:, j]
            if f in CONT:
                e = _bin_edges(ref, act)
                if e is None: raise ValueError(f"bin fail {f}")
                psi_bins.append(e)
                psi_refhist.append(np.histogram(ref, bins=e)[0])
                tgt_binidx.append(np.digitize(act, e[1:-1]).astype(np.int16))
            else:
                cats = _cats(ref, act)
                psi_bins.append(cats)
                psi_refhist.append(np.bincount(np.searchsorted(cats, ref), minlength=len(cats)))
                tgt_binidx.append(np.searchsorted(cats, act).astype(np.int16))

        # prediction-space：source internal TEST 为参考，bin 固定
        Dt_std_full = sc.transform(Dt_raw)
        pt_full = mdl.predict_proba(Dt_std_full)[:, 1]
        st_full = mdl.decision_function(Dt_std_full)      # 与 p6c sc_t 完全一致
        s_e = mdl.decision_function(Dte_std)              # source internal test 的 decision function
        ps_ent, pt_ent = ent(M["ps"]), ent(pt_full)
        ps_conf, pt_conf = np.maximum(M["ps"], 1 - M["ps"]), np.maximum(pt_full, 1 - pt_full)
        e1 = _bin_edges(ps_ent, pt_ent); e2 = _bin_edges(ps_conf, pt_conf)
        if e1 is None or e2 is None: raise ValueError("predspace bin fail")
        jsent = dict(edges=e1, refhist=np.histogram(ps_ent, bins=e1)[0],
                     tgtidx=np.digitize(pt_ent, e1[1:-1]).astype(np.int16))
        jsconf = dict(edges=e2, refhist=np.histogram(ps_conf, bins=e2)[0],
                      tgtidx=np.digitize(pt_conf, e2[1:-1]).astype(np.int16))

        # outcome model（source train 上拟合，冻结；不重训）
        om = LogisticRegression(max_iter=2000, random_state=SEED).fit(Dtr_std, M["ytr"])
        m_e = np.clip(om.predict_proba(Dte_std)[:, 1], 1e-9, 1 - 1e-9)
        m_t_full = np.clip(om.predict_proba(sc.transform(Dt_raw))[:, 1], 1e-9, 1 - 1e-9)

        art = dict(task_id=tid, organism=t["organism"], drug=t["drug"],
                   source_site=t["source_site"], target_site=t["target_site"],
                   n_target=len(Xt), n_source=len(Ds_std),
                   Ds_ref_std=Ds_ref_std, mmd_src=mmd_src, Kaa=Kaa, bw2=M["bw2"],
                   Dtr_std=Dtr_std, Dva_std=Dva_std, Dte_std=Dte_std,
                   ytr=M["ytr"], yva=M["yva"], yte=M["yte"],
                   ps=M["ps"], m_e=m_e, s_e=s_e,
                   ps_ent=ps_ent, ps_conf=ps_conf,
                   Dt_std=Dt_std_full, pt_full=pt_full, st_full=st_full, m_t_full=m_t_full,
                   psi_bins=psi_bins, psi_refhist=psi_refhist, tgt_binidx=tgt_binidx,
                   jsent=jsent, jsconf=jsconf)
        art["_prov"] = _prov_signature()     # ★ 2026-10-06：写入产出代码签名，供 p10b 校验
        with open(os.path.join(CACHE, f"{tid}.pkl"), "wb") as f:
            pickle.dump(art, f, protocol=4)
        ok += 1
    except Exception as e:
        bad.append((tid, f"{type(e).__name__}: {e}"))
    if n_done % 40 == 0: p(f"    ... {n_done}/{len(TASKS)}  [{time.time()-t0:.0f}s]")

p(f"\n[写出] 缓存成功 {ok} / {len(TASKS)}   失败 {len(bad)}")
if bad: p("  失败：" + str(bad[:8]))
sz = sum(os.path.getsize(os.path.join(CACHE, f)) for f in os.listdir(CACHE)) / 1e9
p(f"  缓存目录大小 = {sz:.2f} GB   [{time.time()-t0:.0f}s]")
# ★ 2026-10-06：**完成标记**。只有跑完且无失败时才写；p10b 用它判断"缓存是否整套来自当前代码"。
#   为什么不只看每个 pkl 里的 _prov：p10b 需要一次廉价检查"整套缓存是否可信"，
#   逐个 load 202 个 pkl 代价过高。该文件在缓存目录内，不是任何被报告/被引用的产物。
if ok == len(TASKS) and not bad:
    with open(os.path.join(CACHE, "_prov.txt"), "w", encoding="utf-8") as f:
        f.write(_prov_signature() + "\n")
    p(f"  [prov] 完成标记已写出：_prov.txt = {_prov_signature()[:16]}…")
else:
    # 失败时**删除**标记，避免"半套缓存 + 完成标记"这种最危险的组合
    _mk = os.path.join(CACHE, "_prov.txt")
    if os.path.exists(_mk):
        os.remove(_mk)
    p(f"  [prov] 有 {len(bad)} 个失败 ⇒ **不写**完成标记（并清除旧标记）")
rep.close()
print("done")
