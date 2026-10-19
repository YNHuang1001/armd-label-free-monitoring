# -*- coding: utf-8 -*-
r"""锚点口径复算（约束集 v1.2 的依据）
    —— 用冻结的 episode_definition.assign_episodes 重跑三项依赖 episode 划分的分析。

R1  p2_10 闸门 1：raw vs episode 的 S->R before/after
R2  p2_11 配对层面窗内暴露率（约束集 D-3 的 14-180 天口径）
R3  tt_01_02 landmark cohort + positivity（只确认 STOP 结论，不做 PS/不估效应）

输出 report_p6_anchor_rerun.txt / p6_R1_gate1.csv / p6_R2_pair_exposure.csv / p6_R3_positivity.csv
"""
import os, sys, time
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from episode_definition import assign_episodes

BASE = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))); OUT = os.path.dirname(os.path.abspath(__file__))
t0 = time.time()
SUS = {"SUSCEPTIBLE": "S", "RESISTANT": "R", "INTERMEDIATE": "I"}
DRUGS = {"CIP": "ciprofloxacin", "LVX": "levofloxacin", "CRO": "ceftriaxone",
         "FEP": "cefepime", "TZP": "piperacillin_tazobactam",
         "SXT": "trimethoprim_sulfamethoxazole", "GEN": "gentamicin"}
BINS = [(14, 30), (31, 90), (91, 180)]; BIN_LAB = [f"{a}-{b}" for a, b in BINS]
CANON = ["FQ", "Ceph3", "Ceph4", "Carbapenem", "BLBLI", "SXT", "Aminoglycoside"]
FQ_MAIN = {"CIPROFLOXACIN", "LEVOFLOXACIN"}
WASHOUT, D1, D2 = 180, 31, 180
MIN_EXPOSED, MIN_EVENTS = 200, 30

MGB_CODE = {"CIP": "FQ", "LVX": "FQ", "MXF": "FQ", "DLX": "FQ", "CRO": "Ceph3", "CAZ": "Ceph3",
            "CTX": "Ceph3", "CPD": "Ceph3", "CDR": "Ceph3", "FDC": "Ceph3", "FEP": "Ceph4",
            "CFZ": "Ceph1", "LEX": "Ceph1", "CFR": "Ceph1", "CXM": "Ceph2", "FOX": "Ceph2",
            "CTT": "Ceph2", "MEM": "Carbapenem", "IPM": "Carbapenem", "ETP": "Carbapenem",
            "TZP": "BLBLI", "SAM": "BLBLI", "AMC": "BLBLI", "CZA": "BLBLI", "C/T": "BLBLI",
            "MVB": "BLBLI", "I-R": "BLBLI", "SXT": "SXT", "TMP": "SXT", "GEN": "Aminoglycoside",
            "TOB": "Aminoglycoside", "AMK": "Aminoglycoside", "NIT": "Nitrofuran", "FOF": "Fosfomycin",
            "ATM": "Monobactam", "AMP": "Penicillin", "AMX": "Penicillin", "PEN": "Penicillin",
            "AZM": "Macrolide", "CLR": "Macrolide", "ERY": "Macrolide", "CLI": "Lincosamide",
            "DOX": "Tetracycline", "MIN": "Tetracycline", "TET": "Tetracycline", "OMC": "Tetracycline",
            "TGC": "Tetracycline", "MTZ": "Nitroimidazole", "VAN": "Glycopeptide", "DAP": "Lipopeptide",
            "LZD": "Oxazolidinone", "TZD": "Oxazolidinone", "CST": "Polymyxin", "FLC": "Azole",
            "VRC": "Azole", "POS": "Azole", "ITC": "Azole", "AMB": "Polyene", "MFG": "Echinocandin",
            "CAS": "Echinocandin", "NAF": "Antistaph", "OXA": "Antistaph", "DCX": "Antistaph", "DAL": "Antistaph"}
SUBTYPE = {"Fluoroquinolone": "FQ", "Cephalosporin Gen1": "Ceph1", "Cephalosporin Gen2": "Ceph2",
           "Cephalosporin Gen3": "Ceph3", "Cephalosporin Gen4": "Ceph4", "Carbapenem": "Carbapenem",
           "Beta Lactam Combo": "BLBLI", "Sulfonamide Combo": "SXT", "Sulfonamide": "SXT",
           "Folate Synthesis Inhibitor": "SXT", "Aminoglycoside": "Aminoglycoside",
           "Nitrofuran": "Nitrofuran", "Urinary Antiseptic": "Nitrofuran", "Fosfomycin": "Fosfomycin",
           "Monobactam": "Monobactam", "Penicillin": "Penicillin", "Macrolide": "Macrolide",
           "Lincosamide": "Lincosamide", "Tetracycline": "Tetracycline",
           "Nitroimidazole": "Nitroimidazole", "Glycopeptide": "Glycopeptide",
           "Oxazolidinone": "Oxazolidinone", "Polymyxin": "Polymyxin"}
NONSYS = ["SILVER SULFADIAZINE", "SULFADIAZINE SILVER", "SILVER", "CREAM", "OINTMENT", "TOPICAL",
          "LOTION", "PASTE", "POWDER", "MEDICATED", "BURN", "DERMAL", "CUTANEOUS", "BARRIER",
          "OPHTHALMIC", "OPTHALMIC", "OPTH", "OCULAR", "INTRAVITREAL", "-DEXAMETHASONE",
          "DEXAMETHASONE-", " OTIC", "EAR DROP", "EYE DROP", "IRRIGATION", "IRRIG", "RINSE",
          "WASH", "BLADDER", "INTRAVESICAL", "PERITONEAL DIALYSIS", "CATHETER", "IMPREGNATED",
          "BEADS", "CEMENT", "INHAL", "NEBULIZ", "AEROSOL", "TOBI", "VAGINAL", "VAGINA"]
AMBIG = {"OFLOXACIN", "GATIFLOXACIN", "BESIFLOXACIN", "NORFLOXACIN",
         "GENTAMICIN-SODIUM CITRATE", "NEOMYCIN", "POLYMYXIN B"}
def route_ok(n):
    up = str(n).upper().strip()
    return not (any(x in up for x in NONSYS) or up in AMBIG)
def to_day(s):
    v = pd.to_datetime(s, errors="coerce", utc=True); return v.dt.tz_localize(None).dt.normalize()
def norm_abx(s):
    return (s.astype(str).str.lower().str.strip()
            .str.replace(r"\s+", "_", regex=True).str.replace(r"\+", "_", regex=False))
def norm_spec(s):
    return (s.astype(str).str.upper().str.strip()
            .str.replace("RESPIRATORY_TRACT", "RESP", regex=False).str.replace("RESPIRATORY", "RESP", regex=False))

rep = open(os.path.join(OUT, "report_p6_anchor_rerun.txt"), "w", encoding="utf-8")
def p(*a):
    line = " ".join(str(x) for x in a); print(line, flush=True); rep.write(line + "\n"); rep.flush()

p("=" * 112); p("锚点口径复算 —— 约束集 v1.2 依据（算法 = episode_definition.assign_episodes）"); p("=" * 112)

R1, R2, R3 = [], [], []
for site in ["MGB", "SU", "UT"]:
    MGB = site == "MGB"
    P_ = os.path.join(BASE, "ARMD-MGB" if MGB else ("ARMD-Stanford" if site == "SU" else "ARMD-UTSW"))
    tcol = ("order_time_jittered_utc_shifted" if MGB else
            ("order_time_jittered_utc" if site == "SU" else "order_time_jittered"))
    if MGB:
        d = pd.read_csv(os.path.join(P_, "microbiology_cohort_deid_tj_updated.csv"),
                        usecols=["anon_id", "organism", "antibiotic", "AST_pheno",
                                 "culture_description", tcol], low_memory=False, encoding="utf-8-sig")
        d = d[d["organism"].astype(str).str.upper().str.strip().eq("ESCHERICHIA COLI")].rename(columns={"AST_pheno": "sus"})
    else:
        d = pd.read_csv(os.path.join(P_, "microbiology_cultures_cohort.csv"),
                        usecols=["anon_id", tcol, "culture_description", "organism", "antibiotic",
                                 "susceptibility"], low_memory=False, encoding="utf-8-sig")
        o = d["organism"].astype(str).str.upper().str.strip()
        d = d[o.str.startswith("ESCHERICHIA COLI") | o.eq("ESBL ESCHERICHIA COLI")].rename(columns={"susceptibility": "sus"})
    d["cdate"] = to_day(d[tcol]); d = d[d["cdate"].notna()].copy()
    d["abx"] = norm_abx(d["antibiotic"]); d["spec"] = norm_spec(d["culture_description"])
    d["s"] = d["sus"].astype(str).str.upper().str.strip().map(SUS).fillna("X")
    d["organism"] = "EC"

    key = ["anon_id", "cdate", "spec", "abx"]
    u = d.drop_duplicates(key + ["s"])
    us = u[u["s"].isin(["S", "R"])]
    g = pd.DataFrame({"n_all": u.groupby(key)["s"].nunique(),
                      "n_sr": us.groupby(key)["s"].nunique()}).fillna(0)
    ph = us.groupby(key)["s"].first().reindex(g.index)
    g["ph"] = ph.where((g["n_all"] == 1) & (g["n_sr"] == 1)).fillna("X")
    g = g.reset_index()
    wide0 = g.pivot_table(index=["anon_id", "cdate", "spec"], columns="abx", values="ph", aggfunc="first")
    X0 = wide0.reset_index(); X0["organism"] = "EC"
    n_unit = len(X0)
    X0["_blk"] = assign_episodes(X0)
    X0["epid"] = (X0["anon_id"].astype(str) + "|" + X0["spec"].astype(str) + "|" + X0["_blk"].astype(str))
    n_epi = X0["epid"].nunique()
    p(f"\n### {site}: unit = {n_unit:,}  ->  anchor episode = {n_epi:,}  "
      f"（折叠 {n_unit-n_epi:,} = {(n_unit-n_epi)/n_unit*100:.2f}%）")
    X1 = X0.sort_values("cdate").groupby("epid").head(1).drop(columns=["_blk"]).reset_index(drop=True)

    def pairs(X):
        x = X.reset_index(drop=True)
        o = np.lexsort((x["cdate"].values, x["spec"].astype(str).values, x["anon_id"].astype(str).values))
        x = x.iloc[o].reset_index(drop=True)
        same = (x["anon_id"].values[1:] == x["anon_id"].values[:-1]) & (x["spec"].values[1:] == x["spec"].values[:-1])
        gp = (x["cdate"].values[1:] - x["cdate"].values[:-1]).astype("timedelta64[D]").astype(float)
        k = np.where(same & (gp >= 14) & (gp <= 180))[0]
        return x.iloc[k].reset_index(drop=True), x.iloc[k + 1].reset_index(drop=True), gp[k]

    def trans(A, B, gap, tag):
        out = []
        for dk, col in DRUGS.items():
            if col not in A.columns: continue
            a = A[col].astype(str).values; b = B[col].astype(str).values
            m = ((a == "S") | (a == "R")) & ((b == "S") | (b == "R"))
            if m.sum() == 0: continue
            t = np.where(m, a + "->" + b, None)
            for (lo, hi), bl in zip(BINS, BIN_LAB):
                sel = m & (gap >= lo) & (gap <= hi)
                if sel.sum() == 0: continue
                tt = t[sel]; aa = a[sel]
                nS = int((aa == "S").sum()); nR = int((aa == "R").sum())
                out.append(dict(site=site, tag=tag, drug=dk, bin=bl, n=int(sel.sum()),
                                SR=int((tt == "S->R").sum()), RS=int((tt == "R->S").sum()),
                                n_S0=nS, n_R0=nR,
                                SR_pct=(int((tt == "S->R").sum()) / nS * 100) if nS else np.nan))
        return out

    A0, B0, gp0 = pairs(X0); A1, B1, gp1 = pairs(X1)
    R1 += trans(A0, B0, gp0, "raw") + trans(A1, B1, gp1, "anchor_episode")
    p(f"    配对：raw = {len(gp0):,}  anchor_episode = {len(gp1):,}")

    # ---- 暴露 ----
    if MGB:
        ex = pd.read_csv(os.path.join(P_, "prior_abx_deid_tj.csv"),
                         usecols=["anon_id", "drug_code", "medication_name", "last_dose_DT_shifted"],
                         low_memory=False, encoding="utf-8-sig")
        ex["canon"] = ex["drug_code"].astype(str).str.upper().str.strip().map(MGB_CODE)
        ex["expo_date"] = to_day(ex["last_dose_DT_shifted"])
    else:
        off = "medication_time_to_cultureTime" if site == "SU" else "medication_time_to_culturetime"
        ex = pd.read_csv(os.path.join(P_, "microbiology_cultures_antibiotic_subtype_exposure.csv"),
                         usecols=["anon_id", "antibiotic_subtype", "medication_name", tcol, off],
                         low_memory=False, encoding="utf-8-sig")
        ex["canon"] = ex["antibiotic_subtype"].astype(str).str.strip().map(SUBTYPE)
        ex["expo_date"] = to_day(ex[tcol]) - pd.to_timedelta(pd.to_numeric(ex[off], errors="coerce"), unit="D")
    ex["med"] = ex["medication_name"].astype(str).str.upper().str.strip()
    ex = ex[ex["expo_date"].notna() & ex["canon"].notna()]
    ex = ex[ex["med"].map(route_ok)].copy()
    ex["is_fq"] = ex["med"].isin(FQ_MAIN)
    exf = ex[["anon_id", "expo_date", "canon", "is_fq"]].drop_duplicates()
    pe = {k: (v["expo_date"].values, v["canon"].values, v["is_fq"].values)
          for k, v in exf.groupby("anon_id", sort=False)}

    # ---- R2：配对层面 (C0,C1] 暴露率（用 anchor episode 队列）----
    pr = pd.DataFrame({"anon_id": A1["anon_id"].values, "c0": A1["cdate"].values, "c1": B1["cdate"].values})
    hit = {c: np.zeros(len(pr), dtype=bool) for c in CANON}
    for i in range(len(pr)):
        r = pe.get(pr["anon_id"].iat[i])
        if r is None: continue
        dd, cc, _ = r
        sel = (dd > pr["c0"].iat[i]) & (dd <= pr["c1"].iat[i])
        if not sel.any(): continue
        for cv in pd.unique(cc[sel]):
            if cv in hit: hit[cv][i] = True
    row = dict(site=site, n_pairs=len(pr))
    for c in CANON: row[c] = hit[c].mean() * 100
    R2.append(row)
    p("    R2 配对层面暴露率: " + "  ".join(f"{c}={row[c]:.1f}%" for c in CANON))

    # ---- R3：landmark ----
    ep = X1.sort_values(["anon_id", "spec", "cdate"]).reset_index(drop=True)
    ep["prev_days"] = ep.groupby(["anon_id", "spec"])["cdate"].diff().dt.days
    nxt = np.full(len(ep), -1, dtype=int)
    for _, idx in ep.groupby(["anon_id", "spec"], sort=False).indices.items():
        idx = np.sort(np.asarray(idx)); dts = ep["cdate"].values[idx]
        for k in range(len(idx)):
            j = np.searchsorted(dts, dts[k] + np.timedelta64(D1, "D"), side="left")
            if j < len(idx) and int((dts[j] - dts[k]) / np.timedelta64(1, "D")) <= D2:
                nxt[idx[k]] = int(idx[j])
    ep["next_idx"] = nxt
    for od, ocol in {"CIP": "ciprofloxacin", "CRO": "ceftriaxone"}.items():
        if ocol not in ep.columns: continue
        b2 = ep[ep[ocol].isin(["S", "R"])]
        b2 = b2[b2[ocol] == "S"]
        wash = b2[(b2["prev_days"].isna()) | (b2["prev_days"] > WASHOUT)]
        arm = np.array(["U"] * len(wash), dtype=object)
        c0v = wash["cdate"].values
        for i, aid in enumerate(wash["anon_id"].values):
            r = pe.get(aid)
            if r is None: continue
            dd, cc, ff = r
            s = (dd > c0v[i]) & (dd <= c0v[i] + np.timedelta64(30, "D"))
            if s.any(): arm[i] = "T" if ff[s].any() else "C"
        wash = wash.assign(arm=arm)
        tc = wash[wash["arm"].isin(["T", "C"])].copy()
        okv = np.zeros(len(tc), dtype=bool); ov = np.array([None] * len(tc), dtype=object)
        epcol = ep[ocol].values
        for i, ni in enumerate(tc["next_idx"].values):
            if ni >= 0 and epcol[ni] in ("S", "R"): okv[i] = True; ov[i] = epcol[ni]
        tc = tc.assign(out_ok=okv, out=ov)
        for a in ["T", "C"]:
            x = tc[tc["arm"] == a]
            asc = int(x["out_ok"].sum()); ev = int((x.loc[x["out_ok"], "out"] == "R").sum())
            R3.append(dict(site=site, outcome=od, arm=a, n=len(x), n_asc=asc, events=ev,
                           rate=ev / max(asc, 1) * 100,
                           gate=("PASS" if (a == "T" and len(x) >= MIN_EXPOSED and ev >= MIN_EVENTS)
                                 else ("FAIL" if a == "T" else ""))))
        p(f"    R3 {od}: T n={int((tc['arm']=='T').sum()):,} 可判定={int(tc.loc[tc['arm']=='T','out_ok'].sum()):,} "
          f"events={int((tc.loc[(tc['arm']=='T')&tc['out_ok'],'out']=='R').sum()):,} | "
          f"C n={int((tc['arm']=='C').sum()):,} events={int((tc.loc[(tc['arm']=='C')&tc['out_ok'],'out']=='R').sum()):,}")
    del d, u, us, g, wide0, X0, X1, ex, exf, pe, ep
    p(f"  [{time.time()-t0:.0f}s]")

pd.DataFrame(R1).to_csv(os.path.join(OUT, "p6_R1_gate1.csv"), index=False)
pd.DataFrame(R2).to_csv(os.path.join(OUT, "p6_R2_pair_exposure.csv"), index=False)
pd.DataFrame(R3).to_csv(os.path.join(OUT, "p6_R3_positivity.csv"), index=False)

# ================= 对照展示 =================
p("\n" + "=" * 112); p("【R1】闸门 1：raw vs anchor_episode（S->R 事件数）"); p("=" * 112)
T1 = pd.DataFrame(R1)
p(f"  {'站点':<5}{'药':<5}{'窗':>8}{'raw S->R':>10}{'anchor S->R':>12}{'变化':>8}{'raw率':>9}{'anchor率':>10}")
for s in ["MGB", "SU", "UT"]:
    for dk in ["CIP", "LVX", "CRO", "SXT", "GEN"]:
        for bl in BIN_LAB:
            r = T1[(T1.site == s) & (T1.tag == "raw") & (T1.drug == dk) & (T1["bin"] == bl)]
            a = T1[(T1.site == s) & (T1.tag == "anchor_episode") & (T1.drug == dk) & (T1["bin"] == bl)]
            if not len(r) or not len(a): continue
            r, a = r.iloc[0], a.iloc[0]
            p(f"  {s:<5}{dk:<5}{bl:>8}{r['SR']:>10,}{a['SR']:>12,}{a['SR']-r['SR']:>+8,}"
              f"{r['SR_pct']:>8.1f}%{a['SR_pct']:>9.1f}%")

p("\n" + "=" * 112); p("【R2】配对层面窗内暴露率 —— 对照约束集 D-3 的旧值"); p("=" * 112)
OLD2 = {"MGB": {"FQ": 16.0, "Ceph3": 23.6, "Ceph4": 6.2, "Carbapenem": 5.5, "BLBLI": 11.3, "SXT": 16.4, "Aminoglycoside": 1.0},
        "SU":  {"FQ": 16.9, "Ceph3": 6.4, "Ceph4": 0.1, "Carbapenem": 1.3, "BLBLI": 11.3, "SXT": 13.0, "Aminoglycoside": 1.7},
        "UT":  {"FQ": 29.2, "Ceph3": 8.4, "Ceph4": 0.3, "Carbapenem": 0.7, "BLBLI": 5.9, "SXT": 1.1, "Aminoglycoside": 1.6}}
T2 = pd.DataFrame(R2).set_index("site")
p(f"  {'类':<18}{'MGB旧':>8}{'MGB新':>8}{'SU旧':>8}{'SU新':>8}{'UT旧':>8}{'UT新':>8}")
for c in CANON:
    p(f"  {c:<18}" + "".join(f"{OLD2[s][c]:>8.1f}{T2.loc[s, c]:>8.1f}" for s in ["MGB", "SU", "UT"]))
p(f"  {'配对 n':<18}" + "".join(f"{'':>8}{int(T2.loc[s,'n_pairs']):>8,}" for s in ["MGB", "SU", "UT"]))

p("\n" + "=" * 112); p("【R3】landmark positivity —— 对照 tt_01_02 的旧值"); p("=" * 112)
T3 = pd.DataFrame(R3)
OLD3 = {("CIP", "MGB"): (6172, 113), ("CIP", "SU"): (796, 20), ("CIP", "UT"): (3416, 72),
        ("CRO", "MGB"): (6486, 57), ("CRO", "SU"): (611, 6), ("CRO", "UT"): (3745, 18)}
p(f"  {'outcome':<9}{'站点':<5}{'n_T旧':>8}{'n_T新':>8}{'ev旧':>7}{'ev新':>7}   判定")
for od in ["CIP", "CRO"]:
    for s in ["MGB", "SU", "UT"]:
        r = T3[(T3.outcome == od) & (T3.site == s) & (T3.arm == "T")]
        if not len(r): continue
        r = r.iloc[0]; o = OLD3[(od, s)]
        ok = (r["n"] >= MIN_EXPOSED) and (r["events"] >= MIN_EVENTS)
        p(f"  {od:<9}{s:<5}{o[0]:>8,}{int(r['n']):>8,}{o[1]:>7,}{int(r['events']):>7,}   "
          f"{'PASS' if ok else 'FAIL'}  (旧判定 {'PASS' if o[1]>=MIN_EVENTS else 'FAIL'})")

p(f"\n完成，耗时 {time.time()-t0:.0f}s")
rep.close(); print("done", flush=True)
