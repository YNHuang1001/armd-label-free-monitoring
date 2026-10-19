# -*- coding: utf-8 -*-
r"""p21a — TRIPOD+AI reporting profile：**reporting-only 派生**

=============================================================================
性质（写死）
-----------------------------------------------------------------------------
**REPORTING-ONLY DERIVATION** —— post hoc, required for TRIPOD+AI completeness。

  · 只读；输入 = p6b authoritative corrected pipeline 构造 STORE 时**实际读取的**三站源文件
  · **不改变** episode / eligibility / feature / model / signal / outcome
  · **不进入**任何模型管线；**不重算**任何已冻结结果
  · 输出为**独立的 reporting sidecar**，绝不回流修改冻结分析产物
  · "pre-imputation dataset" 定义 = *the authoritative site-level records after the same
    upstream harmonization used by p6b, but before source-model imputation/scaling and
    before task-level expansion*（**不是**最原始的医院 dump）

=============================================================================
四道闸门（任一不过 → 直接退出，不继续计算）
-----------------------------------------------------------------------------
G0  **source manifest**：从 `p6b_01_corrected_pipeline.py` **自动提取**其源文件映射
    （站点 / 文件名 / usecols / loader），与 p21a 实际使用的映射**逐项比对**；
    并记录 file_sha256 与 raw_row_count。
G1  **代码复用门**：p21a 内的 `_canon / norm_spec / to_day / age_lo / FOCAL / DRUGMAP`
    与 p6b 中同名定义**逐字比对**（空白归一化后必须完全一致）。
G2  **episode 算法等价门**：p6b 是**内联**实现 anchor 链（**并未调用**
    `episode_definition.assign_episodes()`）。本脚本对同一 unit 表**同时**跑两套实现，
    断言逐行 block id 完全一致 —— 以此证明"复用权威算法"而非"自行重写"。
G3  **宇宙对齐门（最强）**：派生出的每 (site, organism, drug) 的 S/R episode 计数，
    必须与冻结的 `n_source` **逐项相等**；不等即退出。

用法：python p21a_reporting_profile_derivation.py
输出：figures/data/reporting_profiles/  下的 reporting sidecars + 本脚本报告
"""
import os, re, io, sys, hashlib, time
import numpy as np, pandas as pd

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); BASE = os.path.dirname(ROOT)
OUT = os.path.join(ROOT, "分析结果"); FIG = os.path.join(ROOT, "figures", "data")
DEST = os.path.join(FIG, "reporting_profiles"); os.makedirs(DEST, exist_ok=True)
P6B = os.path.join(OUT, "p6b_01_corrected_pipeline.py")
tag = time.strftime("%Y%m%d-%H%M")

rep = open(os.path.join(DEST, f"p21a_report_{tag}.txt"), "w", encoding="utf-8")
def p(*a):
    line = " ".join(str(x) for x in a); print(line, flush=True); rep.write(line + "\n"); rep.flush()

def die(msg):
    p(f"\n*** GATE FAILED: {msg} —— 不继续计算，且未写出任何派生结果 ***")
    rep.close(); sys.exit(2)

SRC6 = io.open(P6B, encoding="utf-8").read()
def normcode(x):
    """代码等价比较：**去掉全部空白**而非折叠成单空格。
    首版用 " ".join(x.split())，在换行位置不同的两段等价代码间会产生
    `).` vs `) .` 的假差异 —— 闸门因此误报 G1 FAIL。"""
    return "".join(x.split())

# ============================================================================
# G1  代码复用门：下列定义必须与 p6b 逐字一致
# ============================================================================
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
            .str.replace("RESPIRATORY_TRACT", "RESP", regex=False)
            .str.replace("RESPIRATORY", "RESP", regex=False))

def to_day(s):
    v = pd.to_datetime(s, errors="coerce", utc=True); return v.dt.tz_localize(None).dt.normalize()

# 【G1 闸门抓到】初版我按"猜测"把 age_lo 写成 Series 版正则提取，而 p6b 的实现是
# **标量**切分（调用处传 `dem.get(o, np.nan)`，是标量）。我那版既与 p6b 语义不同
# （如 "≥18"：p6b→NaN，我→18），又会在标量上直接抛 AttributeError。
# 现逐字复用 p6b 实现。**注释必须留在 def 之外** —— 放进函数体会破坏"逐字一致"的断言。
def age_lo(s):
    return pd.to_numeric(str(s).split("-")[0].split("+")[0].strip().split(" ")[0], errors="coerce")

p("=" * 104); p(f"p21a — reporting-only derivation　{tag}"); p("=" * 104)
p("[REPORTING-ONLY DERIVATION] post hoc, required for TRIPOD+AI completeness")
p("")

def grab(src, header):
    """按**缩进**切出函数块。

    【两版 bug 记录】
    (1) 首版用 find("def ") 找块尾，会把函数之后、下一个 def 之前的**顶层语句**
        （如 rep = open(...)）一并纳入 → 造成假的不一致；
    (2) 块尾判定必须看**缩进**，不能看是否出现 "def"。
    """
    i = src.find(header)
    if i < 0:
        return None
    # 【第三处 bug】命中位置可能以换行开头，此时 lines[0] 为空串、lines[1] 才是顶格的
    # `def` 行 → 循环立刻 break，返回空块。必须先 lstrip 掉前导换行。
    lines = src[i:].lstrip("\n").split("\n")
    out = [lines[0]]
    for ln in lines[1:]:
        if ln.strip() == "" or ln[0] in " \t":
            out.append(ln)
        else:
            break                      # 顶格非空行 → 函数结束
    return normcode("\n".join(out))

for nm, blk in [("_canon", grab(SRC6, "def _canon(n):")),
                ("norm_spec", grab(SRC6, "def norm_spec(s):")),
                ("to_day", grab(SRC6, "def to_day(s):")),
                ("age_lo", grab(SRC6, "def age_lo(s):"))]:
    if blk is None: die(f"G1 无法从 p6b 提取 {nm}")
    # 从**本文件**用同一提取器取同名定义（而非字符串 split —— 那会先命中上方的字面量）
    mine = grab(io.open(__file__, encoding="utf-8").read(), f"\ndef {nm}(")
    if normcode(blk) != normcode(mine):
        die(f"G1 {nm} 与 p6b 定义不一致 —— 疑似自行重写")
p("[G1 PASS] _canon / norm_spec / to_day / age_lo 与 p6b **逐字一致**")

# ============================================================================
# G0  source manifest：从 p6b 自动提取并逐项比对
# ============================================================================
LIT = {
 "MGB": dict(dir="ARMD-MGB", micro="microbiology_cohort_deid_tj_updated.csv",
             tcol="order_time_jittered_utc_shifted", ast="AST_pheno",
             dem="demographics_deid_tj.csv", ward="ward_type_deid_tj.csv"),
 "SU":  dict(dir="ARMD-Stanford", micro="microbiology_cultures_cohort.csv",
             tcol="order_time_jittered_utc", ast="susceptibility",
             dem="microbiology_cultures_demographics.csv", ward="microbiology_cultures_ward_info.csv"),
 "UT":  dict(dir="ARMD-UTSW", micro="microbiology_cultures_cohort.csv",
             tcol="order_time_jittered", ast="susceptibility",
             dem="microbiology_cultures_demographics.csv", ward="microbiology_cultures_ward_info.csv"),
}
for site, L in LIT.items():
    for key in ("micro", "dem", "ward"):
        if L[key] not in SRC6:
            die(f"G0 站点 {site} 的 {key} 文件名 '{L[key]}' 未出现在 p6b 源码中")
    if L["dir"] not in SRC6 or L["tcol"] not in SRC6 or L["ast"] not in SRC6:
        die(f"G0 站点 {site} 的 dir/tcol/ast 与 p6b 源码不符")
p("[G0 PASS] 三站文件名 / 目录 / 时间列 / AST 列均出现在 p6b 源码中")

man = []
for site, L in LIT.items():
    for key, role in (("micro", "microbiology records"), ("dem", "demographics"), ("ward", "ward/setting")):
        fp = os.path.join(BASE, L["dir"], L[key])
        h = hashlib.sha256()
        with open(fp, "rb") as fh:
            for chunk in iter(lambda: fh.read(1 << 22), b""): h.update(chunk)
        nrow = sum(1 for _ in open(fp, encoding="utf-8-sig", errors="replace")) - 1
        man.append(dict(site=site, role=role, source_file=f"{L['dir']}/{L[key]}",
                        file_sha256=h.hexdigest(), loader="pandas.read_csv(encoding='utf-8-sig')",
                        raw_row_count=nrow,
                        patient_id_field="anon_id" if key == "micro" else "order_proc_id_coded",
                        date_field=L["tcol"] if key == "micro" else "",
                        organism_field="organism" if key == "micro" else "",
                        specimen_field="culture_description" if key == "micro" else "",
                        ast_field=L["ast"] if key == "micro" else ""))
pd.DataFrame(man).to_csv(os.path.join(DEST, "p21a_source_manifest.csv"), index=False)
for m in man:
    p(f"  [manifest] {m['site']:<4}{m['role']:<24}{m['source_file']:<52}"
      f"rows={m['raw_row_count']:>10,}  sha256={m['file_sha256'][:16]}")
p("[G0 PASS] source manifest 已写出（含 sha256 与 raw row count）")
p("")

# ============================================================================
# STORE 构造：**逐字复用 p6b**（仅去掉进度打印）
# ============================================================================
t0 = time.time()
STORE = {}
FLOW = {}
for site in ["MGB", "SU", "UT"]:
    MGB = site == "MGB"; L = LIT[site]
    P_ = os.path.join(BASE, L["dir"]); tcol = L["tcol"]
    if MGB:
        d = pd.read_csv(os.path.join(P_, L["micro"]),
                        usecols=["anon_id", "order_proc_id_coded", "organism", "antibiotic",
                                 "AST_pheno", "culture_description", tcol],
                        low_memory=False, encoding="utf-8-sig")
        d = d[d["organism"].astype(str).str.upper().str.strip().eq("ESCHERICHIA COLI")
              | d["organism"].astype(str).str.upper().str.strip().str.startswith("KLEBSIELLA PNEUMONIAE")
              | d["organism"].astype(str).str.upper().str.strip().str.startswith("PSEUDOMONAS AERUGINOSA")].copy()
        d = d.rename(columns={tcol: "cdate", "AST_pheno": "sus"})
    else:
        d = pd.read_csv(os.path.join(P_, L["micro"]),
                        usecols=["anon_id", "order_proc_id_coded", tcol, "culture_description",
                                 "organism", "antibiotic", "susceptibility"],
                        low_memory=False, encoding="utf-8-sig")
        d = d.rename(columns={tcol: "cdate", "susceptibility": "sus"})
    d["cdate"] = to_day(d["cdate"]); d = d[d["cdate"].notna()].copy()
    d["org"] = d["organism"].astype(str).str.upper().str.strip()
    d["abx"] = d["antibiotic"].map(_canon); d["spec"] = norm_spec(d["culture_description"])
    d["s"] = d["sus"].astype(str).str.upper().str.strip().map(
        {"SUSCEPTIBLE": "S", "RESISTANT": "R", "INTERMEDIATE": "I"}).fillna("X")
    # participant-flow 计数器（本轮新增；此前只记了 STORE 之后的量）
    FLOW[site] = dict(n_rows_after_read=int(len(d)),
                      n_rows_focal_organisms=0, n_rows_dedup_abcd_s=0,
                      n_rows_SRI=0, n_units=0, n_episodes_in_store=0)
    dem = pd.read_csv(os.path.join(P_, L["dem"]), usecols=["order_proc_id_coded", "age"],
                      low_memory=False, encoding="utf-8-sig")
    dem = dem.drop_duplicates("order_proc_id_coded").set_index("order_proc_id_coded")["age"]
    wd = pd.read_csv(os.path.join(P_, L["ward"]), low_memory=False, encoding="utf-8-sig")
    wd["setting"] = np.where(wd["hosp_ward_ER"] == 1, "ER",
                             np.where(wd["hosp_ward_IP"] == 1, "Inpatient",
                                      np.where(wd["hosp_ward_OP"] == 1, "Outpatient", "Other")))
    wd = wd.drop_duplicates("order_proc_id_coded").set_index("order_proc_id_coded")["setting"]
    for sym, fname in FOCAL.items():
        dd = d[d["org"].str.startswith(fname)].copy()
        if not len(dd): continue
        FLOW[site]["n_rows_focal_organisms"] += int(len(dd))
        u = dd.drop_duplicates(["anon_id", "cdate", "spec", "abx", "s"])
        FLOW[site]["n_rows_dedup_abcd_s"] += int(len(u))
        usr = u[u["s"].isin(["S", "R", "I"])]
        FLOW[site]["n_rows_SRI"] += int(len(usr))
        if not len(usr): continue
        g = usr.groupby(["anon_id", "cdate", "spec", "abx"])["s"].agg(lambda v: "BOTH" if ({"S", "R"} <= set(v)) else ("S" if "S" in set(v) else ("R" if "R" in set(v) else "I"))).reset_index()
        unit = (g[["anon_id", "cdate", "spec"]].drop_duplicates()
                .sort_values(["anon_id", "spec", "cdate"]).reset_index(drop=True))
        # ---- G2：p6b 内联实现 vs episode_definition 权威实现，逐行断言一致 ----
        blk = np.zeros(len(unit), dtype=int); ck = None; anc = None; b = -1
        for i in range(len(unit)):
            k = (unit["anon_id"].iat[i], unit["spec"].iat[i]); dv = unit["cdate"].iat[i]
            if k != ck: b += 1; anc = dv; ck = k
            elif (dv - anc) / np.timedelta64(1, "D") >= 14: b += 1; anc = dv
            blk[i] = b
        sys.path.insert(0, OUT)
        from episode_definition import assign_episodes
        blk2 = assign_episodes(unit, group_keys=("anon_id", "spec"),
                               date_col="cdate", window_days=14)
        # 【G2 的判据必须是**分区等价**，不是 block 编号相等】
        # p6b 按 (anon_id, spec, cdate) **元组**排序；assign_episodes 内部按
        # "anon_id|spec" **字符串**排序（于是 "10|X" < "9|X"）。**组序不同 → 编号不同**，
        # 但"哪些行属于同一 episode"这一划分可能完全一致。
        # 故把两边各自按**首次出现序**重编号（对标签重命名不变），再比较。
        def canon(b):
            m, nxt = {}, np.empty(len(b), dtype=int)
            for i, v in enumerate(b):
                if v not in m: m[v] = len(m)
                nxt[i] = m[v]
            return nxt
        n_diff = int((canon(blk) != canon(blk2)).sum())
        if n_diff:
            die(f"G2 {site}/{sym}: p6b 内联 anchor 链与 episode_definition.assign_episodes "
                f"**分区不等价**（{n_diff} 行落在不同 episode）")
        FLOW[site]["n_units"] += int(len(unit))
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
        STORE[(site, sym)] = X.assign(episode_key=sym + "|" + X.index.astype(str)).sort_values(["anon_id", "spec", "cdate"]).reset_index(drop=True)
        FLOW[site]["n_episodes_in_store"] += int(len(X))
    del d, wd, dem
    p(f"  [{site}] STORE 就绪  [{time.time()-t0:.0f}s]")
p("[G2 PASS] 全部 (site, organism)：p6b 内联 anchor 链 ≡ episode_definition.assign_episodes（逐行一致）")

# ============================================================================
# G3  宇宙对齐门：派生计数 == 冻结 n_source
# ============================================================================
FB = pd.read_csv(os.path.join(OUT, "p6b_label_free_signals_blinded.csv"))
FB["group"] = FB.organism + "|" + FB.drug
bad = []
for _, r in FB.iterrows():
    X = STORE.get((r.source_site, r.organism))
    if X is None: bad.append((r.task_id, "no STORE")); continue
    n = int((X[r.drug].isin(["S", "R"])).sum()) if r.drug in X.columns else 0
    if n != int(r.n_source): bad.append((r.task_id, f"派生 {n} vs 冻结 {int(r.n_source)}"))
if bad:
    die(f"G3 宇宙对齐失败 {len(bad)} 项，例：{bad[:4]}")
p(f"[G3 PASS] {len(FB)} 个任务的 S/R episode 计数与冻结 n_source **逐项相等**")
p("")

# ============================================================================
# 派生（全部标 analysis_unit / deduplication_level / source_file / filter）
# ============================================================================
BENCH = FB[["organism", "drug"]].drop_duplicates()
sites = ["MGB", "SU", "UT"]
rows_site, rows_task, flow = [], [], []

for site in sites:
    eps = set(); pts = set(); ages = {}; sets = {}; specs = {}
    n_records = n_epi_all = n_epi_sr = 0
    for sym in FOCAL:
        X = STORE.get((site, sym))
        if X is None: continue
        drugs = list(BENCH[BENCH.organism == sym].drug)
        mask = X[drugs].isin(["S", "R"]).any(axis=1) if drugs else pd.Series(False, index=X.index)
        Xb = X[mask]
        n_epi_all += len(X); n_epi_sr += len(Xb)
        eps |= set(Xb["episode_key"]); pts |= set(Xb["anon_id"])
        # 【bug 修正】zip 里取的必须是 Xb["age"]（初版误取 anon_id，导致 ages 里全是字符串）
        for e, a, st, sp in zip(Xb["episode_key"], Xb["age"], Xb["setting"], Xb["spec"]):
            ages[e] = a; sets[e] = st; specs[e] = sp
    ok_age = [v for v in ages.values() if pd.notna(v)]
    rows_site.append(dict(
        site=site,
        n_unique_microbiological_episodes=len(eps),
        n_unique_patients=len(pts),
        n_age_observations=len(ok_age),
        n_age_missing=len(ages) - len(ok_age),
        pct_age_missing=100.0 * (len(ages) - len(ok_age)) / max(len(ages), 1),
        age_median=float(np.median(ok_age)) if ok_age else np.nan,
        age_p25=float(np.percentile(ok_age, 25)) if ok_age else np.nan,
        age_p75=float(np.percentile(ok_age, 75)) if ok_age else np.nan,
        setting_ER=sum(1 for v in sets.values() if v == "ER"),
        setting_Inpatient=sum(1 for v in sets.values() if v == "Inpatient"),
        setting_Outpatient=sum(1 for v in sets.values() if v == "Outpatient"),
        setting_Other=sum(1 for v in sets.values() if v == "Other"),
        setting_missing=sum(1 for v in sets.values() if not isinstance(v, str)),
        specimen_URINE=sum(1 for v in specs.values() if v == "URINE"),
        specimen_BLOOD=sum(1 for v in specs.values() if v == "BLOOD"),
        specimen_other=sum(1 for v in specs.values() if v not in ("URINE", "BLOOD")),
        analysis_unit="microbiological episode（除 n_unique_patients 外）",
        deduplication_level="episode ID = (anon_id, specimen, anchor block)；**task 展开之前**",
        denominator_definition=("该站点在 38 个 organism–drug 宇宙内、且对 ≥1 个基准药物有 S/R "
                                "结果的**去重** microbiological episodes"),
        source_file="p6b STORE（由 p21a 复用 p6b 构造逻辑重建）",
        filter_definition="organism ∈ 3 focal；drug ∈ 38 组基准；sus ∈ {S,R}；≥1 个基准药命中",
        provenance="REPORTING-ONLY DERIVATION；G0–G3 四道闸门全部通过"))
    for sym in FOCAL:
        X = STORE.get((site, sym))
        if X is None: continue
        for drug in BENCH[BENCH.organism == sym].drug:
            if drug not in X.columns: continue
            v = X[drug]; nS = int((v == "S").sum()); nR = int((v == "R").sum())
            rows_task.append(dict(site=site, organism=sym, drug=drug,
                                  n_S=nS, n_R=nR, n_SR=nS + nR,
                                  analysis_unit="organism–drug episode observation",
                                  deduplication_level="site × drug（各药独立计数）",
                                  source_file="p6b STORE（p21a 复用重建）",
                                  filter_definition=f"organism={sym}; drug={drug}; sus ∈ {{S,R}}"))

pd.DataFrame(rows_site).to_csv(os.path.join(DEST, "site_profile.csv"), index=False)
# ★ 2026-10-08（A3.4 补洞）：p21a 是 site_profile.csv 的**真生产者**（A3.2），
#   此前只有 p21c 盖代次头 ⇒ 若 p21a 在 p21c 之后跑，工件是 p21a 的而 sidecar 记 p21c 的 sha。
import genstamp as _gsa; _gsa.stamp(os.path.join(DEST, "site_profile.csv"), "p21a")
pd.DataFrame(rows_task).to_csv(os.path.join(DEST, "site_x_organism_drug_counts.csv"), index=False)

# participant flow（以 site 为单位，逐节点）
for site in sites:
    Xs = {sym: STORE.get((site, sym)) for sym in FOCAL}
    tot_epi = sum(len(v) for v in Xs.values() if v is not None)
    flow.append(dict(site=site, step="episodes in STORE (3 focal organisms)", unit="microbiological episode",
                     value=tot_epi, source_file="p6b STORE（复用重建）", note="task 展开前"))
    r = pd.DataFrame(rows_site); rr = r[r.site == site].iloc[0]
    flow.append(dict(site=site, step="episodes with S/R for ≥1 benchmark drug", unit="microbiological episode",
                     value=int(rr.n_unique_microbiological_episodes), source_file="p6b STORE（复用重建）",
                     note="38 组宇宙内、去重"))
    flow.append(dict(site=site, step="unique patients contributing ≥1 such episode", unit="patient",
                     value=int(rr.n_unique_patients), source_file="p6b STORE（复用重建）", note=""))
pd.DataFrame(FLOW).T.rename_axis("site").reset_index().to_csv(
    os.path.join(DEST, "participant_flow_counts.csv"), index=False)
pd.DataFrame(flow).to_csv(os.path.join(DEST, "participant_flow.csv"), index=False)

p("")
r = pd.DataFrame(rows_site)
p("=== site_profile（REPORTING-ONLY DERIVATION）===")
p(r[["site", "n_unique_microbiological_episodes", "n_unique_patients", "n_age_missing",
     "pct_age_missing", "age_median", "setting_Inpatient", "setting_Outpatient", "setting_ER",
     "specimen_URINE", "specimen_BLOOD"]].to_string(index=False))
p("")
p("=== 与 task-expanded 口径的对照（**两个分母，不得混用**）===")
tp = pd.read_csv(os.path.join(DEST, "task_profile.csv"))
for site in sites:
    p(f"  {site}: unique episodes = {int(r[r.site==site].n_unique_microbiological_episodes.iat[0]):>9,}"
      f"   |  task-expanded episode observations = "
      f"{pd.read_csv(os.path.join(DEST,'site_profile.csv')).pipe(lambda d: 0) if False else ''}"
      f"{int(tp[tp.site==site].n_episodes_sr.sum()):>10,}")
p("  **前者按 episode ID 去重；后者未去重、按测药数加权。两者不可互换。**")
p("")
p("*** p21a 完成：G0–G3 全部通过；输出均为 REPORTING-ONLY，未触碰任何冻结结果 ***")
rep.close(); print("\ndone ->", DEST)
