# -*- coding: utf-8 -*-
r"""Figure 1E —— 四个 model/feature 配置在**修正后管线**下的 degradation 对照

**为什么必须重跑**：
GPT6 指定的 Figure 1 Panel E 需要"History+LR 与 History+LGBM 的 degradation distribution 对照，
证明 benchmark failure 不是 classifier-specific"。现有 `p5_all_configs.csv` 属于
**protocol v1.1 明确标记为 superseded / preprocessing-mismatched 的 benchmark v0**，
不得进入主图（其 P1 History+LR 与 p6d 的 d_auroc 最大差 0.0084，仅 61/202 完全一致）。

**本脚本**：沿用 p6b_01_corrected_pipeline.py 的**逐字**数据加载、build_hist、prep；
仅把 p5 的 2x2 配置矩阵（tier x algo）重新执行一遍，并修正 preprocessing：
  * age 用 **source-train 中位数** 填补（p5 用的是含 target 的 F["age"].median()，违反 §4.2）
  * StandardScaler **只 fit 在 source train**
  * Platt 校准只用 source validation，锁死后套到 target
  * 其余（split seed、LGBM 超参、Platt 形式）与 p5 完全一致

**不改变任何设计**：只为 Figure 1E 产出 protocol-compliant 的数值。

输出：p12_configs_corrected.csv / report_p12_configs.txt
"""
import os, time, warnings
import numpy as np, pandas as pd
warnings.filterwarnings("ignore")
import lightgbm as lgb
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score, average_precision_score

BASE = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))); OUT = os.path.dirname(os.path.abspath(__file__))
t0 = time.time(); SEED = 20261005; AUROC_MIN = 0.65
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

rep = open(os.path.join(OUT, "report_p12_configs.txt"), "w", encoding="utf-8")
def p(*a):
    line = " ".join(str(x) for x in a); print(line, flush=True); rep.write(line + "\n"); rep.flush()

# ===== 数据加载（逐字复制自 p6b）=====
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

# ===== build_hist / prep（逐字复制自 p6b）=====
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

def design(X, H, tier):
    """与 p5 同构，但 age **不再在此处填补**（改为 source-train 中位数，见下）"""
    F = pd.DataFrame(index=X.index)
    F["age"] = pd.to_numeric(X["age"], errors="coerce")
    F["set_IP"] = (X["setting"] == "Inpatient").astype(int)
    F["set_OP"] = (X["setting"] == "Outpatient").astype(int)
    F["set_ER"] = (X["setting"] == "ER").astype(int)
    F["sp_URINE"] = (X["spec"] == "URINE").astype(int)
    F["sp_BLOOD"] = (X["spec"] == "BLOOD").astype(int)
    if tier == "hist":
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

CONFIGS = [("P1 History+LR", "hist", "lr"), ("S1 Core+LR", "core", "lr"),
           ("S2 History+LGBM", "hist", "lgbm"), ("S3 Core+LGBM", "core", "lgbm")]

# ===== 逐配置建模型 =====
CM = pd.read_csv(os.path.join(OUT, "p6b_corrected_universe.csv"))     # 228 candidate tasks
KEYS = sorted(set(zip(CM.organism, CM.drug, CM.source_site)))
p(f"\n  candidate tasks = {len(CM)}   模型键 = {len(KEYS)}   配置 = {len(CONFIGS)}")

rows = []
for cfgname, tier, algo in CONFIGS:
    MC = {}
    for sym, drug, src in KEYS:
        try:
            Xs, Hs = prep(STORE[(src, sym)], drug)
            if len(Xs) < 500: continue
            D = design(Xs, Hs, tier); y = (Xs[drug].values == "R").astype(int)
            pats = Xs["anon_id"].values; upat = pd.unique(pats)
            rng = np.random.default_rng(SEED); perm = rng.permutation(len(upat))
            ntr, nva = int(0.70 * len(upat)), int(0.15 * len(upat))
            s_tr = set(upat[perm[:ntr]]); s_va = set(upat[perm[ntr:ntr + nva]])
            i_tr = np.array([i for i, a in enumerate(pats) if a in s_tr])
            i_va = np.array([i for i, a in enumerate(pats) if a in s_va])
            i_te = np.array([i for i, a in enumerate(pats) if a not in s_tr and a not in s_va])
            if min(len(i_tr), len(i_va), len(i_te)) < 50: continue
            Dall = D.values.astype(float).copy()
            # === protocol §4.2：age 只用 source-train 中位数 ===
            af = float(np.nanmedian(Dall[i_tr, 0])) if np.isfinite(Dall[i_tr, 0]).any() else 0.0
            Dall[:, 0] = np.where(np.isfinite(Dall[:, 0]), Dall[:, 0], af)
            sc = StandardScaler().fit(Dall[i_tr])     # 只 fit 在 source train
            if algo == "lr":
                m = LogisticRegression(max_iter=2000, C=1.0, random_state=SEED).fit(sc.transform(Dall[i_tr]), y[i_tr])
                f = lambda Z, m=m, sc=sc: m.predict_proba(sc.transform(Z))[:, 1]
            else:
                m = lgb.LGBMClassifier(n_estimators=300, num_leaves=31, learning_rate=0.05,
                                       random_state=SEED, verbose=-1, n_jobs=4).fit(Dall[i_tr], y[i_tr])
                f = lambda Z, m=m: m.predict_proba(Z)[:, 1]
            eps = 1e-6
            pv = f(Dall[i_va])
            lv = np.log(np.clip(pv, eps, 1 - eps) / (1 - np.clip(pv, eps, 1 - eps))).reshape(-1, 1)
            cal = LogisticRegression(max_iter=1000).fit(lv, y[i_va])
            def calib(pr, cal=cal):
                z = np.log(np.clip(pr, eps, 1 - eps) / (1 - np.clip(pr, eps, 1 - eps))).reshape(-1, 1)
                return cal.predict_proba(z)[:, 1]
            src_auroc = float(roc_auc_score(y[i_te], calib(f(Dall[i_te]))))
            MC[(sym, drug, src)] = dict(Dall=Dall, f=f, calib=calib, af=af, src_auroc=src_auroc,
                                        yte=y[i_te], i_te=i_te)
        except Exception as e:
            p(f"  [warn] {cfgname} {(sym,drug,src)}: {type(e).__name__}: {e}")
    p(f"  [{cfgname}] 模型建成 {len(MC)}，eligible(source AUROC>={AUROC_MIN}) = "
      f"{sum(1 for v in MC.values() if v['src_auroc'] >= AUROC_MIN)}  [{time.time()-t0:.0f}s]")
    for _, t in CM.iterrows():
        sym, drug, src, tgt = t["organism"], t["drug"], t["source_site"], t["target_site"]
        M = MC.get((sym, drug, src))
        if M is None: continue
        try:
            Xt, Ht = prep(STORE[(tgt, sym)], drug)
            if len(Xt) < 300: continue
            Dt = design(Xt, Ht, tier).values.astype(float)
            Dt[:, 0] = np.where(np.isfinite(Dt[:, 0]), Dt[:, 0], M["af"])
            yt = (Xt[drug].values == "R").astype(int)
            tgt_auroc = float(roc_auc_score(yt, M["calib"](M["f"](Dt))))
            rows.append(dict(organism=sym, drug=drug, source=src, target=tgt,
                             task_id=f"{sym}_{drug}_{src}2{tgt}", config=cfgname,
                             tier=tier, algo=algo, n_target=len(Xt),
                             src_auroc=M["src_auroc"], tgt_auroc=tgt_auroc,
                             d_auroc=tgt_auroc - M["src_auroc"],
                             eligible=M["src_auroc"] >= AUROC_MIN))
        except Exception:
            continue

R = pd.DataFrame(rows)
R.to_csv(os.path.join(OUT, "p12_configs_corrected.csv"), index=False)
p(f"\n[写出] p12_configs_corrected.csv  行数 = {len(R)}")
p("\n" + "=" * 104)
p("【Figure 1E】四配置 degradation 对照（修正后管线）")
p("=" * 104)
p(f"{'config':<18}{'n_task':>8}{'n_elig':>8}{'median':>9}{'P10':>8}{'P90':>8}{'P90-P10':>9}"
  f"{'harmful':>9}{'harm%':>8}")
for cfgname, _, _ in CONFIGS:
    s = R[R.config == cfgname]
    e = s[s.eligible]
    p(f"{cfgname:<18}{len(s):>8}{len(e):>8}{e.d_auroc.median():>+9.3f}"
      f"{e.d_auroc.quantile(.1):>+8.3f}{e.d_auroc.quantile(.9):>+8.3f}"
      f"{e.d_auroc.quantile(.9)-e.d_auroc.quantile(.1):>9.3f}"
      f"{int((e.d_auroc<=-0.05).sum()):>9}{100*(e.d_auroc<=-0.05).mean():>7.1f}%")

p("\n  在 P1 的 202 eligible universe 上对齐（同一批 task，便于直接比较）:")
P1 = R[(R.config == "P1 History+LR") & R.eligible].task_id
p(f"{'config':<18}{'n':>6}{'median':>9}{'P90-P10':>9}{'harmful':>9}")
for cfgname, _, _ in CONFIGS:
    s = R[(R.config == cfgname) & (R.task_id.isin(P1))]
    p(f"{cfgname:<18}{len(s):>6}{s.d_auroc.median():>+9.3f}"
      f"{s.d_auroc.quantile(.9)-s.d_auroc.quantile(.1):>9.3f}"
      f"{int((s.d_auroc<=-0.05).sum()):>9}")

# 与 superseded p5 的对照
p("\n  [对照] superseded p5_all_configs.csv（不得进主图，仅作差异记录）:")
try:
    C5 = pd.read_csv(os.path.join(OUT, "p5_all_configs.csv"))
    C5 = C5[C5.organism.notna()]
    for cfgname, _, _ in CONFIGS:
        s5 = C5[C5.config == cfgname]
        s12 = R[(R.config == cfgname) & (R.task_id.isin(P1))]
        m = s5.merge(s12, on=["organism", "drug", "source", "target"], suffixes=("_p5", "_new"))
        if len(m):
            d = (m.d_auroc_p5 - m.d_auroc_new).abs()
            p(f"    {cfgname:<18} 配对 {len(m)}  |dAUC| max={d.max():.4f} mean={d.mean():.5f}")
except Exception as e:
    p(f"    （无法对照：{e}）")
rep.close()
print("done")
