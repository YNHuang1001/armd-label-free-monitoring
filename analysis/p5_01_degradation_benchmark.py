# -*- coding: utf-8 -*-
r"""第五方向 — 跨院模型失效 benchmark：真实 degradation 生成（不计算任何 label-free 指标）

冻结规格
  Feature tiers
    Core              : age + care setting + specimen
    History-augmented : Core + 4 类历史（窗口 = index episode 前 15-365 天）
                        prior same-organism episode 有/无；count；
                        time since most recent（分箱）；most recent prior AST for 同一目标药（None/S/I/R）
  Model 2x2（primary 冻结为 History + Logistic）
    P1  History + Logistic   <- PRIMARY
    S1  Core + Logistic
    S2  History + LightGBM
    S3  Core + LightGBM
  划分：patient-level；source site 固定 70/15/15；所有预处理/调参/校准只用 source train+val
  校准：source validation 上做一次 Platt，锁死后再应用于 target；target labels 绝不参与校准
  ΔAUROC = AUROC_target − AUROC_source_internal_test   （<0 表示下降）
  task 资格：AUROC_source >= 0.65
  闸门（仅看 PRIMARY）：N_harmful(Δ<=-0.05) >= 30 且 N_stable(|Δ|<=0.02) >= 30
                        且 P90(Δ) − P10(Δ) >= 0.10

输出 report_p5_degradation.txt / p5_tasks_primary.csv / p5_all_configs.csv
     p5_gate.csv / p5_task_characteristics.csv
"""
import os, time, warnings
import numpy as np
import pandas as pd
warnings.filterwarnings("ignore")
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, average_precision_score, log_loss, brier_score_loss
from sklearn.preprocessing import StandardScaler
import lightgbm as lgb

BASE = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
OUT = os.path.dirname(os.path.abspath(__file__))
t0 = time.time()
SEED = 20261005
HIST_LO, HIST_HI = 15, 365
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

rep = open(os.path.join(OUT, "report_p5_degradation.txt"), "w", encoding="utf-8")
def p(*a):
    line = " ".join(str(x) for x in a)
    print(line, flush=True); rep.write(line + "\n"); rep.flush()

p("=" * 112); p("第五方向 — 跨院模型失效 benchmark：真实 degradation（不计算 label-free 指标）"); p("=" * 112)

# ---------- 1. 逐站构建 episode + 历史特征 ----------
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
    if MGB:
        wd["setting"] = np.where(wd["hosp_ward_ER"] == 1, "ER",
                                 np.where(wd["hosp_ward_IP"] == 1, "Inpatient",
                                          np.where(wd["hosp_ward_OP"] == 1, "Outpatient", "Other")))
    else:
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
        # ---- 14 天 episode 必须建立在**唯一 (患者, 日期, 部位) 单元**上。
        #      若直接对 (日期 x 药物) 行做 cumsum，同一张培养单的不同药物行会落到不同 block，
        #      同一培养单被拆成多个 episode（实测会把 MGB/EC 从 87,590 虚增到 522,602）。
        unit = (g[["anon_id", "cdate", "spec"]].drop_duplicates()
                .sort_values(["anon_id", "spec", "cdate"]).reset_index(drop=True))
        # 确定性 anchor 算法：新 episode 当且仅当 (培养日 - 当前 episode 锚点日) >= 14。
        # 原迭代式启发式在 (0, 5, 20) 这类序列上会在两个解之间振荡，4 轮上限后停在错误答案。
        blk = np.zeros(len(unit), dtype=int)
        cur_key = None; anchor = None; b = -1
        for i in range(len(unit)):
            k = (unit["anon_id"].iat[i], unit["spec"].iat[i]); dv = unit["cdate"].iat[i]
            if k != cur_key:
                b += 1; anchor = dv; cur_key = k
            elif (dv - anchor) / np.timedelta64(1, "D") >= 14:
                b += 1; anchor = dv
            blk[i] = b
        unit["_blk"] = blk
        unit["epid"] = (unit["anon_id"].astype(str) + "|" + unit["spec"].astype(str)
                        + "|" + unit["_blk"].astype(str))
        g = g.merge(unit[["anon_id", "cdate", "spec", "epid"]],
                    on=["anon_id", "cdate", "spec"], how="left")
        wide = g[g["s"].isin(["S", "R"])].pivot_table(index="epid", columns="abx", values="s", aggfunc="first")
        meta = unit.drop_duplicates("epid").set_index("epid")[["anon_id", "cdate", "spec"]]
        X = meta.join(wide, how="left")
        oid_map = dd.drop_duplicates(["anon_id", "cdate", "spec"]).set_index(
            ["anon_id", "cdate", "spec"])["order_proc_id_coded"]
        oids = oid_map.reindex(pd.MultiIndex.from_arrays(
            [X["anon_id"], X["cdate"], X["spec"]])).values
        X["age"] = [age_lo(dem.get(o, np.nan)) for o in oids]
        X["setting"] = wd.reindex(oids).values
        X = X.sort_values(["anon_id", "spec", "cdate"]).reset_index(drop=True)
        STORE[(site, sym)] = X
        p(f"  [{site}/{sym}] unit(患者,日,部位) = {len(unit):,}  episode = {unit['epid'].nunique():,}  折叠 = {len(unit)-unit['epid'].nunique():,}  患者 = {X['anon_id'].nunique():,}  "
          f"药品种类 = {len([c for c in X.columns if c not in ('anon_id','cdate','spec','age','setting')])}")
    del d, wd, dem
    p(f"  [{time.time()-t0:.0f}s]")

# ---------- 2. task 列表（来自 p4 冻结结果） ----------
T4 = pd.read_csv(os.path.join(OUT, "p4_task_drug_site.csv"))
cnt = T4.pivot_table(index=["focal", "drug"], columns="site", values="qualify", aggfunc="first")
for s in ["MGB", "SU", "UT"]:
    if s not in cnt.columns: cnt[s] = False
COMBO = cnt[(cnt["MGB"] == True) & (cnt["SU"] == True) & (cnt["UT"] == True)].index.tolist()
p(f"\n  来自 p4 的冻结 task 组合 = {len(COMBO)}  (organism, drug)")
COMBOS = [(a, b) for a, b in COMBO if (("MGB", a) in STORE or True)]

def prep(X, drug):
    """返回行序一致的 (X_sorted, H)，并按该药 S/R 过滤。调用方必须用返回值。"""
    Xs = X.sort_values(["anon_id", "spec", "cdate"]).reset_index(drop=True)
    H = build_hist(Xs, drug)          # 假定 Xs 已排序
    m = np.isin(Xs[drug].values, ["S", "R"])
    return Xs[m].reset_index(drop=True), H[m].reset_index(drop=True)

def build_hist(X, drug):
    """要求 X 已按 (anon_id, spec, cdate) 排序；返回与 X 行序一致的 H。窗口 (index-365, index-15]"""
    n = len(X)
    prior_any = np.zeros(n); prior_cnt = np.zeros(n)
    tbin = np.full(n, "none", dtype=object)
    past = np.full(n, "none", dtype=object)
    c0 = X["cdate"].values
    for _, idx in X.groupby(["anon_id", "spec"], sort=False).indices.items():
        idx = np.sort(np.asarray(idx))
        for k in range(len(idx)):
            i = idx[k]
            lo = c0[i] - np.timedelta64(HIST_HI, "D")
            hi = c0[i] - np.timedelta64(HIST_LO, "D")
            pr = idx[:k]
            sel = pr[(c0[pr] > lo) & (c0[pr] <= hi)]
            if len(sel) == 0: continue
            prior_any[i] = 1; prior_cnt[i] = len(sel)
            dts = (c0[i] - c0[sel]).astype("timedelta64[D]").astype(int)
            mn = dts.min()
            tbin[i] = ("15-30" if mn <= 30 else "31-90" if mn <= 90
                       else "91-180" if mn <= 180 else "181-365")
            j = sel[np.argmin(dts)]
            v = X[drug].values[j] if drug in X.columns else None
            past[i] = v if v in ("S", "I", "R") else "none"
    H = pd.DataFrame({"prior_any": prior_any, "prior_cnt": prior_cnt,
                      "tbin": tbin, "past_ast": past})
    return H

def design(X, H, tier):
    F = pd.DataFrame(index=X.index)
    F["age"] = pd.to_numeric(X["age"], errors="coerce")
    F["age"] = F["age"].fillna(F["age"].median())
    F["set_IP"] = (X["setting"] == "Inpatient").astype(int)
    F["set_OP"] = (X["setting"] == "Outpatient").astype(int)
    F["set_ER"] = (X["setting"] == "ER").astype(int)
    F["sp_URINE"] = (X["spec"] == "URINE").astype(int)
    F["sp_BLOOD"] = (X["spec"] == "BLOOD").astype(int)
    if tier == "hist":
        F["prior_any"] = H["prior_any"].values
        F["prior_cnt"] = H["prior_cnt"].values
        for b in ["15-30", "31-90", "91-180", "181-365"]:
            F[f"t_{b}"] = (H["tbin"] == b).astype(int).values
        for v in ["S", "I", "R"]:
            F[f"past_{v}"] = (H["past_ast"] == v).astype(int).values
    return F.values.astype(float)

def fit_eval(Xtr, ytr, Xva, yva, Xte, yte, Xtg, ytg, algo):
    if algo == "lr":
        sc = StandardScaler().fit(Xtr)
        m = LogisticRegression(max_iter=2000, C=1.0, random_state=SEED)
        m.fit(sc.transform(Xtr), ytr)
        pv = m.predict_proba(sc.transform(Xva))[:, 1]
        f = lambda Z: m.predict_proba(sc.transform(Z))[:, 1]
    else:
        m = lgb.LGBMClassifier(n_estimators=300, num_leaves=31, learning_rate=0.05,
                               random_state=SEED, verbose=-1, n_jobs=4)
        m.fit(Xtr, ytr)
        pv = m.predict_proba(Xva)[:, 1]
        f = lambda Z: m.predict_proba(Z)[:, 1]
    # Platt 校准（仅用 source validation），锁死后应用
    eps = 1e-6
    lv = np.log(np.clip(pv, eps, 1 - eps) / (1 - np.clip(pv, eps, 1 - eps))).reshape(-1, 1)
    cal = LogisticRegression(max_iter=1000).fit(lv, yva)
    def calib(pr):
        z = np.log(np.clip(pr, eps, 1 - eps) / (1 - np.clip(pr, eps, 1 - eps))).reshape(-1, 1)
        return cal.predict_proba(z)[:, 1]
    out = {}
    for tag, X_, y_ in [("src_test", Xte, yte), ("target", Xtg, ytg)]:
        pr = calib(f(X_))
        out[tag] = dict(n=len(y_), prev=float(y_.mean()), auroc=float(roc_auc_score(y_, pr)),
                        auprc=float(average_precision_score(y_, pr)),
                        logloss=float(log_loss(y_, np.clip(pr, 1e-6, 1 - 1e-6))),
                        brier=float(brier_score_loss(y_, pr)))
    return out

psd = 20261005
rows = []
for (sym, drug) in COMBOS:
    if (("MGB", sym) not in STORE) or (("SU", sym) not in STORE) or (("UT", sym) not in STORE):
        continue
    for src in ["MGB", "SU", "UT"]:
        Xsrc = STORE[(src, sym)]
        if drug not in Xsrc.columns: continue
        Xs2, Hs2 = prep(Xsrc, drug)
        if len(Xs2) < 500: continue
        y = (Xs2[drug].values == "R").astype(int)
        pats = Xs2["anon_id"].values
        upat = pd.unique(pats)
        rng = np.random.default_rng(psd)
        perm = rng.permutation(len(upat))
        ntr, nva = int(0.70 * len(upat)), int(0.15 * len(upat))
        s_tr = set(upat[perm[:ntr]]); s_va = set(upat[perm[ntr:ntr + nva]])
        i_tr = np.array([i for i, a in enumerate(pats) if a in s_tr])
        i_va = np.array([i for i, a in enumerate(pats) if a in s_va])
        i_te = np.array([i for i, a in enumerate(pats) if a not in s_tr and a not in s_va])
        if min(len(i_tr), len(i_va), len(i_te)) < 50: continue
        for tgt in ["MGB", "SU", "UT"]:
            if tgt == src: continue
            Xtg_ = STORE[(tgt, sym)]
            if drug not in Xtg_.columns: continue
            Xt2, Ht2 = prep(Xtg_, drug)
            if len(Xt2) < 300: continue
            ytg = (Xt2[drug].values == "R").astype(int)
            for tier in ["core", "hist"]:
                D = design(Xs2, Hs2, tier); Dt = design(Xt2, Ht2, tier)
                for algo in ["lr", "lgb"]:
                    cfg = ("P1 History+LR" if (tier == "hist" and algo == "lr") else
                           "S1 Core+LR" if (tier == "core" and algo == "lr") else
                           "S2 History+LGBM" if (tier == "hist") else "S3 Core+LGBM")
                    try:
                        r = fit_eval(D[i_tr], y[i_tr], D[i_va], y[i_va],
                                     D[i_te], y[i_te], Dt, ytg, algo)
                    except Exception as e:
                        continue
                    d_auc = r["target"]["auroc"] - r["src_test"]["auroc"]
                    rows.append(dict(organism=sym, drug=drug, source=src, target=tgt, config=cfg,
                                     n_src_test=r["src_test"]["n"], src_auroc=r["src_test"]["auroc"],
                                     src_prev=r["src_test"]["prev"], n_target=r["target"]["n"],
                                     tgt_prev=r["target"]["prev"], tgt_auroc=r["target"]["auroc"],
                                     d_auroc=d_auc,
                                     d_auprc=r["target"]["auprc"] - r["src_test"]["auprc"],
                                     d_logloss=r["target"]["logloss"] - r["src_test"]["logloss"],
                                     d_brier=r["target"]["brier"] - r["src_test"]["brier"]))
        p(f"  [{sym}/{drug}] src={src} 完成  [{time.time()-t0:.0f}s]")

R = pd.DataFrame(rows)
R.to_csv(os.path.join(OUT, "p5_all_configs.csv"), index=False)
PR = R[R["config"] == "P1 History+LR"].copy()
PR["eligible"] = PR["src_auroc"] >= 0.65
PR.to_csv(os.path.join(OUT, "p5_tasks_primary.csv"), index=False)
PR[["organism", "drug", "source", "target", "n_target", "tgt_prev", "src_auroc",
    "tgt_auroc", "d_auroc", "eligible"]].to_csv(
    os.path.join(OUT, "p5_task_characteristics.csv"), index=False)

p("\n" + "=" * 112); p("【A】PRIMARY（History + Logistic）— 按 source/target 汇总"); p("=" * 112)
if len(PR):
    p(f"  总 task（有向）= {len(PR)}；source AUROC>=0.65 的 eligible = {int(PR['eligible'].sum())}")
    E = PR[PR["eligible"]]
    if len(E):
        p(f"  eligible 的 ΔAUROC: median={E['d_auroc'].median():+.3f}  "
          f"P10={E['d_auroc'].quantile(.10):+.3f}  P90={E['d_auroc'].quantile(.90):+.3f}")
        p(E.groupby(["source", "target"])["d_auroc"].agg(["size", "median", "min", "max"]).round(3).to_string())

p("\n" + "=" * 112); p("【B】闸门（仅 PRIMARY，eligible tasks）"); p("=" * 112)
gate = {}
if len(PR) and PR["eligible"].any():
    E = PR[PR["eligible"]]
    n_h = int((E["d_auroc"] <= -0.05).sum()); n_s = int((E["d_auroc"].abs() <= 0.02).sum())
    sp = float(E["d_auroc"].quantile(.90) - E["d_auroc"].quantile(.10))
    gate = dict(n_eligible=len(E), n_harmful=n_h, n_stable=n_s, spread_p90_p10=sp)
    p(f"  eligible tasks          = {len(E):<5} ")
    p(f"  N_harmful (Δ<=-0.05)    = {n_h:<5} (门槛 >=30)   {'PASS' if n_h>=30 else 'FAIL'}")
    p(f"  N_stable  (|Δ|<=0.02)   = {n_s:<5} (门槛 >=30)   {'PASS' if n_s>=30 else 'FAIL'}")
    p(f"  P90-P10(ΔAUROC)         = {sp:.3f} (门槛 >=0.10) {'PASS' if sp>=0.10 else 'FAIL'}")
    p(f"\n  总体判定 = {'PASS' if (n_h>=30 and n_s>=30 and sp>=0.10) else 'FAIL'}")
pd.DataFrame([gate]).to_csv(os.path.join(OUT, "p5_gate.csv"), index=False)

p("\n" + "=" * 112); p("【C】四配置对照（eligible tasks 的 ΔAUROC 中位数）"); p("=" * 112)
if len(R):
    R2 = R.merge(PR[["organism", "drug", "source", "target", "eligible"]],
                 on=["organism", "drug", "source", "target"], how="left")
    p(R2[R2["eligible"] == True].groupby("config")["d_auroc"]
      .agg(["size", "median", lambda x: (x <= -0.05).sum()]).round(3).to_string())

p(f"\n完成，耗时 {time.time()-t0:.0f}s | 输出目录 {OUT}")
rep.close(); print("done", flush=True)
