# -*- coding: utf-8 -*-
r"""第四方向 feasibility scan — label-free 跨院模型失效预测：task-availability matrix

只出表，不训练模型。

表 organism x drug x site（focal = E. coli / K. pneumoniae / P. aeruginosa）
  direct S/R 总数、R 数、S 数、患者数、标本构成、共同 predictor 覆盖率

筛选：N >= 1000, N_R >= 200, N_S >= 200
再看能构成多少个完整的 source -> target transfer task（两端都过筛）
硬闸：>=15 个有效 task，覆盖 >=3 organisms、>=5 drugs、全部三个站点

单位 = isolate（order_proc_id_coded，该药 direct AST 为 S 或 R）
输出 report_p4_taskmatrix.txt / p4_task_drug_site.csv / p4_task_qualify.csv / p4_tasks.csv
"""
import os, time
import numpy as np
import pandas as pd

BASE = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
OUT = os.path.dirname(os.path.abspath(__file__))
t0 = time.time()
N_MIN, R_MIN, S_MIN = 1000, 200, 200
FOCAL = {"EC": "ESCHERICHIA COLI", "KP": "KLEBSIELLA PNEUMONIAE", "PA": "PSEUDOMONAS AERUGINOSA"}
SUSMAP = {"SUSCEPTIBLE": "S", "RESISTANT": "R", "INTERMEDIATE": "I"}

# 三站药名拼写不同（MGB: amoxicillin_clavulanate / TZP: piperacillin_tazobactam /
# SXT: trimethoprim_sulfamethoxazole，而 SU/UT 用斜杠或缩写形式）-> 必须显式映射。
# 顺序敏感：复方制剂必须排在其单组分之前。
DRUGMAP = [
    ("AMC", ["amoxicillinclavulanate", "amoxicillinclavulanicacid", "amoxclavulanate", "amoxicillinclav", "amoxicillinkclavulanateaug"]),
    ("TZP", ["piperacillintazobactam", "piptazobactam", "piperacillintazo", "tazobactam"]),
    ("SXT", ["trimethoprimsulfamethoxazole", "trimethsulf", "trimethoprimsulfa", "cotrimoxazole"]),
    ("SAM", ["ampicillinsulbactam", "ampicsulbactam", "ampicillinsulb"]),
    ("TIM", ["ticarcillinclavulanate", "ticarcilclavulanate", "ticarcillinclavulanicacid", "ticacarcillinclavulanicacid"]),
    ("CZA", ["ceftolozanetazobactam"]),
    ("CZO", ["cefazolin"]), ("CRO", ["ceftriaxone"]), ("CAZ", ["ceftazidime"]),
    ("FEP", ["cefepime"]), ("CTX", ["cefotaxime"]), ("CXM", ["cefuroxime"]),
    ("FOX", ["cefoxitin"]), ("LEX", ["cephalexin"]), ("CEP", ["cephalothin"]),
    ("CIP", ["ciprofloxacin"]), ("LVX", ["levofloxacin"]),
    ("GEN", ["gentamicin", "gentscreen", "gentamicinscreen"]), ("TOB", ["tobramycin", "tobramycinrm"]), ("AMK", ["amikacin", "amikaciniv", "amikacinintravenous"]),
    ("NIT", ["nitrofurantoin"]), ("MEM", ["meropenem"]), ("IPM", ["imipenem"]),
    ("ETP", ["ertapenem"]), ("ATM", ["aztreonam"]), ("TET", ["tetracycline"]),
    ("TGC", ["tigecycline"]), ("FOF", ["fosfomycin"]), ("CHL", ["chloramphenicol"]),
    ("CST", ["colistin", "polymixine", "polymyxine"]),
    ("SAM", ["sulbactam"]), ("AMP", ["ampicillin"]), ("PEN", ["piperacillin"]),
    ("SUL", ["sulfamethoxazole"]), ("TMP", ["trimethoprim"]),
]
def _canon(name):
    up = "".join(ch for ch in str(name).lower() if ch.isalnum())
    for code, keys in DRUGMAP:
        if up in keys:
            return code
    return up
def norm_abx(s):
    return s.map(_canon)
def norm_spec(s):
    return (s.astype(str).str.upper().str.strip()
            .str.replace("RESPIRATORY_TRACT", "RESP", regex=False).str.replace("RESPIRATORY", "RESP", regex=False))

rep = open(os.path.join(OUT, "report_p4_taskmatrix.txt"), "w", encoding="utf-8")
def p(*a):
    line = " ".join(str(x) for x in a)
    print(line, flush=True); rep.write(line + "\n"); rep.flush()

p("=" * 112)
p("第四方向 feasibility — task-availability matrix（label-free 跨院模型失效预测）")
p("=" * 112)
p(f"单位 = isolate（order_proc_id_coded，该药 direct AST ∈ {{S,R}}）。筛选：N>={N_MIN}, R>={R_MIN}, S>={S_MIN}")

rows, cov_rows = [], []
for site in ["MGB", "SU", "UT"]:
    p("\n" + "#" * 112); p(f"### {site}"); p("#" * 112)
    MGB = site == "MGB"
    P_ = os.path.join(BASE, "ARMD-MGB" if MGB else ("ARMD-Stanford" if site == "SU" else "ARMD-UTSW"))
    if MGB:
        d = pd.read_csv(os.path.join(P_, "microbiology_cohort_deid_tj_updated.csv"),
                        usecols=["anon_id", "order_proc_id_coded", "organism", "antibiotic",
                                 "AST_pheno", "culture_description"],
                        low_memory=False, encoding="utf-8-sig")
        d = d.rename(columns={"AST_pheno": "sus"})
    else:
        d = pd.read_csv(os.path.join(P_, "microbiology_cultures_cohort.csv"),
                        usecols=["anon_id", "order_proc_id_coded", "organism", "antibiotic",
                                 "susceptibility", "culture_description"],
                        low_memory=False, encoding="utf-8-sig")
        d = d.rename(columns={"susceptibility": "sus"})
    d["org"] = d["organism"].astype(str).str.upper().str.strip()
    d["abx"] = norm_abx(d["antibiotic"])
    d["spec"] = norm_spec(d["culture_description"])
    d["s"] = d["sus"].astype(str).str.upper().str.strip().map(SUSMAP).fillna("X")

    for sym, fname in FOCAL.items():
        hit = d[d["org"].str.startswith(fname)]
        if not len(hit):
            hit = d[d["org"].str.contains(fname.split()[0], na=False)]
        if not len(hit): continue
        lab = hit["org"].value_counts()
        p(f"\n  --- focal {sym} ({fname}) 标签变体 TOP6: {dict(lab.head(6))}")
        use = hit[hit["org"].str.startswith(fname)]        # 只用精确种名
        n_all = hit["order_proc_id_coded"].nunique()
        n_use = use["order_proc_id_coded"].nunique()
        if n_use < n_all:
            p(f"      精确种名匹配 {n_use:,} / 全部变体 {n_all:,} "
              f"（差 {n_all-n_use:,}，为形态学/耐药标签变体，本扫描排除）")
        sub = use[use["s"].isin(["S", "R"])]
        if not len(sub): continue
        iso = sub.drop_duplicates(["order_proc_id_coded", "abx"])
        g = iso.groupby("abx")["s"].agg(n="size", R=lambda x: int((x == "R").sum()),
                                        S=lambda x: int((x == "S").sum()))
        g["n_pat"] = sub.groupby("abx")["anon_id"].nunique()
        g["top_spec"] = sub.groupby("abx")["spec"].agg(lambda x: x.value_counts().index[0])
        g["top_spec_pct"] = sub.groupby("abx")["spec"].agg(lambda x: x.value_counts(normalize=True).iloc[0] * 100)
        g = g[(g["n"] >= 300)].sort_values("n", ascending=False)
        if not len(g): continue
        p(f"      {sym} 药敏可判 isolate 总数 = {iso['order_proc_id_coded'].nunique():,}；"
          f"药品种类 = {len(g)}")
        p(f"      {'drug':<34}{'N':>8}{'R':>7}{'S':>7}{'R%':>7}{'患者':>8}  {'主标本':<8}{'占比':>6}")
        for abx, r in g.iterrows():
            ok = (r["n"] >= N_MIN) and (r["R"] >= R_MIN) and (r["S"] >= S_MIN)
            p(f"      {abx[:32]:<34}{int(r['n']):>8,}{int(r['R']):>7,}{int(r['S']):>7,}"
              f"{r['R']/r['n']*100:>6.1f}%{int(r['n_pat']):>8,}  {str(r['top_spec'])[:7]:<8}{r['top_spec_pct']:>5.1f}%"
              f"{'  [过筛]' if ok else ''}")
            rows.append(dict(site=site, focal=sym, drug=abx, N=int(r["n"]), R=int(r["R"]), S=int(r["S"]),
                             R_pct=float(r["R"] / r["n"] * 100), n_pat=int(r["n_pat"]),
                             top_spec=r["top_spec"], top_spec_pct=float(r["top_spec_pct"]),
                             qualify=bool(ok)))

    # predictor 覆盖率（最小共同集：age + gender + ward + specimen）
    try:
        dem = pd.read_csv(os.path.join(P_, "demographics_deid_tj.csv" if MGB
                                       else "microbiology_cultures_demographics.csv"),
                          usecols=["order_proc_id_coded", "age", "gender"],
                          low_memory=False, encoding="utf-8-sig")
        dem = dem.drop_duplicates("order_proc_id_coded").set_index("order_proc_id_coded")
        wd = pd.read_csv(os.path.join(P_, "ward_type_deid_tj.csv" if MGB
                                      else "microbiology_cultures_ward_info.csv"),
                         usecols=["order_proc_id_coded"], low_memory=False, encoding="utf-8-sig")
        wset = set(wd["order_proc_id_coded"].unique())
        ids = d["order_proc_id_coded"].drop_duplicates()
        has_age = ids.isin(dem.index[dem["age"].notna()]).mean() * 100
        has_gen = ids.isin(dem.index[dem["gender"].astype(str).str.upper()
                                     .isin(["MALE", "FEMALE", "1", "0", "M", "F"])]).mean() * 100
        has_wd = ids.isin(wset).mean() * 100
        has_sp = d.drop_duplicates("order_proc_id_coded")["spec"].notna().mean() * 100
        cov_rows.append(dict(site=site, has_age=has_age, has_gender=has_gen,
                             has_ward=has_wd, has_specimen=has_sp))
        p(f"\n  [predictor 覆盖率] age={has_age:.1f}%  gender={has_gen:.1f}%  "
          f"ward={has_wd:.1f}%  specimen={has_sp:.1f}%")
    except Exception as e:
        p(f"  [warn] predictor 覆盖 {type(e).__name__}: {e}")
    del d
    p(f"  [{time.time()-t0:.0f}s]")

T = pd.DataFrame(rows)
T.to_csv(os.path.join(OUT, "p4_task_drug_site.csv"), index=False)
pd.DataFrame(cov_rows).to_csv(os.path.join(OUT, "p4_predictor_coverage.csv"), index=False)

p("\n" + "=" * 112); p("【A】过筛的 organism x drug x site cells"); p("=" * 112)
if len(T):
    Q = T[T["qualify"]]
    p(f"  过筛 cell 数 = {len(Q)} / {len(T)}")
    p(Q.pivot_table(index=["focal", "drug"], columns="site", values="N").fillna(0).astype(int).to_string())

p("\n" + "=" * 112); p("【B】能构成 transfer task 的 organism x drug（三站全部过筛）"); p("=" * 112)
if len(T):
    cnt = T.pivot_table(index=["focal", "drug"], columns="site", values="qualify", aggfunc="first")
    for s in ["MGB", "SU", "UT"]:
        if s not in cnt.columns: cnt[s] = False
    full = cnt[(cnt["MGB"] == True) & (cnt["SU"] == True) & (cnt["UT"] == True)]
    p(f"  三站全部过筛的 (organism, drug) 组合 = {len(full)}")
    if len(full):
        det = T.set_index(["focal", "drug"]).loc[full.index][["site", "N", "R", "S", "R_pct"]]
        p(det.round(1).to_string())
    n_task = len(full) * 6
    p(f"\n  -> 有向 transfer task 数（3 站 x 2 方向）= {len(full)} x 6 = {n_task}")
    if len(full):
        p(f"  -> 覆盖 organisms = {sorted(full.index.get_level_values(0).unique())}")
        p(f"  -> 覆盖 drugs     = {sorted(full.index.get_level_values(1).unique())}")

    p("\n" + "=" * 112); p("【C】硬闸判定"); p("=" * 112)
    n_org = full.index.get_level_values(0).nunique() if len(full) else 0
    n_drug = full.index.get_level_values(1).nunique() if len(full) else 0
    ok = (n_task >= 15) and (n_org >= 3) and (n_drug >= 5)
    p(f"  task 数 = {n_task}  (门槛 >=15)          {'PASS' if n_task >= 15 else 'FAIL'}")
    p(f"  organisms = {n_org}  (门槛 >=3)          {'PASS' if n_org >= 3 else 'FAIL'}")
    p(f"  drugs     = {n_drug}  (门槛 >=5)          {'PASS' if n_drug >= 5 else 'FAIL'}")
    p(f"  sites = 3 (MGB/SU/UT)                     PASS")
    p(f"\n  总体判定 = {'PASS' if ok else 'FAIL'}")

p(f"\n完成，耗时 {time.time()-t0:.0f}s | 输出目录 {OUT}")
rep.close(); print("done", flush=True)
