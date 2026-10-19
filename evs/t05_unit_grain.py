# -*- coding: utf-8 -*-
r"""t05_unit_grain.py — 族代码 UNIT：分析单元（analytic unit）与派生报告计数的验证

================================================================================
本模块测什么（生产代码由 evs_load 逐字装载，绝不复制）
================================================================================
  * p6b_01_corrected_pipeline.py 第 83-117 行（站点装载：filter/rename/规格化），
    第 118-144 行（per-organism 建模表构造；其中 BLOCK_STORE = 第 121-144 行）——
    用 L.exec_block 逐字 exec（121-144 逐字包含在 118-144 内，因 BLOCK_STORE 自身第
    123 行有一条裸 `continue`，无法在顶层 exec；下面有断言校验逐字包含）。
  * p21a_reporting_profile_derivation.py 第 305-318 行（per-site episode 并集，
    含第 315 行 `eps |= set(Xb.index)`）—— 同样逐字 exec。

全流程只读；不写任何生产目录外的文件。

================================================================================
一个关键事实（本模块实测发现，决定了测什么）
================================================================================
BLOCK_STORE 里 `X` 与 `STORE[(site,sym)]` **不是同一个对象**：
    p6b:138   X  = meta.join(wide)          -> index = epid 字符串（**有重复行**）
    p6b:144   STORE[(site,sym)] = X.sort_values(...).reset_index(drop=True)
                                             -> index = 行号 0..n-1
下游（p6b 的 prep/design、p21a:276 的 STORE 重建）拿到的是**行号版本**。
本模块因此同时保留两者：
    T[sym]["X"]  = p6b:138 对象（epid index，用于 UNIT-01）
    T[sym]["Xr"] = p6b:144 对象（行号 index，用于 UNIT-04/05/06/07，即生产下游口径）
UNIT-06R 用 (anon_id,cdate,spec) 唯一键把生产 epid 贴回 Xr，并断言与独立参考划分等价。

================================================================================
brute-force 参考实现用了哪条定义（**不抄生产代码**）
================================================================================
  R1 分析单元（episode 划分）
       来源：Methods_v1.md:40 —— "The analytic unit is a microbiological episode,
             defined within each combination of patient, organism and normalized
             specimen site"；以及 episode_definition.py 的 docstring（2026-09-29
             冻结的 single source of truth）：
             组内按日期升序；第一条记录建立 episode、其日期为 anchor；
             后续记录当且仅当 (date - 当前 anchor) < 14 天时并入当前 episode，
             否则开新 episode 并以该记录日期为新 anchor；anchor 在 episode 内不移动。
       实现：本文件 ref_episode_labels()，独立重写（未看生产实现）；另按文件原样装载
       冻结 SSoT episode_definition.assign_episodes，两条独立划分都与生产划分对拍。
       本模块在**单一 organism 表**上运行，故协议 (patient, organism, specimen)
       分组退化为 (patient, specimen)。生产实现是 p6b:127-134 的内联 anchor 链；
       两者只比**划分**（canon 重编号后逐行相等），不比 block 编号（编号依赖组序，
       见 p21a G2 注释）。

  R2 报告计数
       定义（site_profile.csv 的 denominator_definition）：在有 >=1 个基准药 S/R
       结果的 episode 集合上，去重计数 (site, organism, episode) 三元组。
       基准药 = p6b_label_free_signals_blinded.csv 冻结的 (organism, drug) 组（只读）。

  R3 S/R 冲突诊断（reviewer 要求）
       对每个 (episode, drug)：取该 episode 内该药**非缺失**结果构成的集合 {S,I,R}；
       冲突 = 集合恰为 {S,R}、{S,I} 或 {R,I}。
       数据源：BLOCK_STORE 逐字 exec 产出的 `g` 表（per (anon_id,cdate,spec,abx) 的
       单条结果）；冲突判据由本模块独立给出。

================================================================================
输出契约
================================================================================
    RESULT|<test_id>|<PASS|FAIL|WARN>|<expected>|<observed>|<note>
字段内不得出现 "|"（分隔符）。

用法：python t05_unit_grain.py
"""
import importlib.util
import os
import sys
import time
import traceback

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import evs_load as L

BASE = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
T0 = time.time()
COVERED = []      # 实际完成重建的 site
FAILED = []       # 重建失败的 site（绝不虚报覆盖）

# ============================================================================
# 0. 逐字装载生产代码的公共件（顶层 def 与顶层常量）—— 属于"投喂 env"，不是参考实现
# ============================================================================
_canon = L.load_def(L.P6B, "_canon")
norm_spec = L.load_def(L.P6B, "norm_spec")
to_day = L.load_def(L.P6B, "to_day")
age_lo = L.load_def(L.P6B, "age_lo")

_ns = {"np": np, "pd": pd}
exec(compile(L.grab(L.P6B, 45, 59), "<verbatim:p6b:45-59>", "exec"), _ns)
FOCAL, DRUGMAP = _ns["FOCAL"], _ns["DRUGMAP"]
for _f in (_canon, norm_spec, to_day, age_lo):
    _f.__globals__.update({"pd": pd, "np": np, "DRUGMAP": DRUGMAP})

# 独立权威（冻结 SSoT，2026-09-29）：episode 定义的唯一实现来源。
# 按路径原样装载 episode_definition.py 本身（不抄写、不改写其规则），
# UNIT-06R 用它交叉核对参考划分，使 correct 目标不再只依赖自写的 paraphrase。
# 若补丁副本未携带该文件，则回退到真实生产目录中的同一冻结文件。
_EPI_PATH = L.EPI if os.path.exists(L.EPI) else os.path.join(L.OUT, "episode_definition.py")
_epi_spec = importlib.util.spec_from_file_location("evs_episode_definition_ssot", _EPI_PATH)
_epi_mod = importlib.util.module_from_spec(_epi_spec)
_epi_spec.loader.exec_module(_epi_mod)
assign_episodes = _epi_mod.assign_episodes

# 只读的冻结输入/冻结产物
FB_CSV = os.path.join(L.OUT, "p6b_label_free_signals_blinded.csv")
SITE_PROFILE_CSV = os.path.join(L.PROJ, "figures", "data", "reporting_profiles", "site_profile.csv")

SRC_UNIT_DEF = ("definition R1: Methods_v1.md:40 (analytic unit = microbiological episode) "
                "+ episode_definition.py docstring (anchor chain, 14-day window)")
SRC_R2 = ("definition R2: site_profile.csv denominator_definition = distinct "
          "(site, organism, episode) among episodes with S/R for >=1 benchmark drug")
SRC_R3 = ("definition R3 (reviewer diagnostic): set of non-missing AST results within "
          "episode x drug; conflict iff the set is exactly {S,R}, {S,I} or {R,I}")

# ============================================================================
# 1. brute-force 参考实现（从定义出发，独立重写）
# ============================================================================
WINDOW_DAYS = 14


def ref_episode_labels(df):
    """R1 的独立实现：anchor-chain episode 划分（对齐 df 行序的整数标签）。"""
    anon = df["anon_id"].astype(str).values
    spec = df["spec"].astype(str).values
    dts = pd.to_datetime(df["cdate"]).values
    groups = {}
    for i in range(len(df)):
        groups.setdefault((anon[i], spec[i]), []).append(i)
    lab = np.full(len(df), -1, dtype=int)
    b = -1
    for k in groups:
        idx = sorted(groups[k], key=lambda i: dts[i])
        anchor = None
        for n, i in enumerate(idx):
            if n == 0:
                b += 1
                anchor = dts[i]
            elif (dts[i] - anchor) / np.timedelta64(1, "D") >= WINDOW_DAYS:
                b += 1
                anchor = dts[i]
            lab[i] = b
    return lab


def canon(labels):
    """把任意标签重编号为首次出现序（对标签重命名不变），用于比较两份划分。"""
    m = {}
    out = np.empty(len(labels), dtype=int)
    for i, v in enumerate(labels):
        j = m.get(v)
        if j is None:
            j = len(m)
            m[v] = j
        out[i] = j
    return out


# ============================================================================
# 2. 生产代码重建（逐字 exec）
# ============================================================================
def build_site(site):
    """逐字 exec p6b:83-117，返回该站点的 (d, dem, wd, src)。"""
    env = dict(np=np, pd=pd, os=os, BASE=BASE, site=site, FOCAL=FOCAL,
               _canon=_canon, norm_spec=norm_spec, to_day=to_day)
    out, src = L.exec_block((L.P6B, 83, 117), env, names=["d", "dem", "wd"])
    return out["d"], out["dem"], out["wd"], src


def build_org(d, dem, wd, site, sym):
    """逐字 exec BLOCK_STORE（p6b:121-144，含在 118-144 区间内）。

    BLOCK_STORE 自身含裸 `continue`（第 123 行），不能在顶层 exec；
    故按生产原文 exec 118-144（`for sym, fname in FOCAL.items():` 单次迭代），
    并断言 121-144 逐字包含其中。
    """
    env = dict(np=np, pd=pd, d=d, dem=dem, wd=wd, site=site, age_lo=age_lo,
               FOCAL={sym: FOCAL[sym]}, STORE={})
    out, src = L.exec_block((L.P6B, 118, 144), env,
                            names=["X", "unit", "g", "wide", "meta", "STORE"])
    block = L.grab(*L.BLOCK_STORE)
    assert all(ln in src for ln in block.split("\n") if ln.strip()), \
        "BLOCK_STORE 文本未逐字包含在 exec 区间内"
    out["_src"] = src
    out["Xr"] = out["STORE"].get((site, sym)) if out.get("STORE") else None
    return out


def prod_site_union(site, STORE, BENCH):
    """逐字 exec p21a:305-318（第 315 行 eps |= set(Xb.index)）。"""
    env = dict(np=np, pd=pd, FOCAL=FOCAL, BENCH=BENCH, STORE=STORE, sites=[site])
    out, src = L.exec_block((L.P21A, 305, 318), env,
                            names=["eps", "pts", "n_epi_all", "n_epi_sr"])
    return out, src


def attach_epid(Xr, unit):
    """把生产 epid 贴回 Xr 行序（unit 与 Xr 一一对应；键 (anon_id,cdate,spec) 唯一）。"""
    umap = unit.set_index(["anon_id", "cdate", "spec"])["epid"]
    key = pd.MultiIndex.from_arrays([Xr["anon_id"].values, Xr["cdate"].values, Xr["spec"].values])
    return umap.reindex(key).values


# ============================================================================
# 3. 结果辅助
# ============================================================================
def warn(tid, expected, observed, note=""):
    L.emit(tid, "WARN", expected, observed, note)


def frozen_anchor(tid, ok, expected, observed, note=""):
    """冻结产物锚点（current-state 参照，非缺陷测试）：期望来自**未修补**管线派生的冻结产物。
    真实生产树上仍按硬门处理（重建必须复现冻结产物）；
    在补丁副本上（EVS_PROD_DIR）不硬判 —— 正确修补会合法地改变这些计数，故降为
    characterisation WARN，避免在正确实现上产生假 FAIL（套件不得变成状态钉死器）。"""
    if ok or not L.IS_PATCHED:
        return L.check(tid, ok, expected, observed, note)
    return warn(tid, "patched tree / frozen artifact predates the fix: " + expected,
                observed, note)


# ============================================================================
# UNIT-00：重建保真度（与冻结 n_source 逐项对账）
# ============================================================================
def fidelity(site, T, fb):
    rows = fb[fb.source_site == site]
    n_ok = n_bad = 0
    worst = 0
    bad_ex = []
    for _, r in rows.iterrows():
        t = T.get(r.organism)
        if t is None or r.drug not in t["Xr"].columns:
            n_bad += 1
            if len(bad_ex) < 3:
                bad_ex.append(f"{r.organism}/{r.drug}:no-table-or-column")
            continue
        n = int(t["Xr"][r.drug].isin(["S", "R"]).sum())
        if n == int(r.n_source):
            n_ok += 1
        else:
            n_bad += 1
            worst = max(worst, abs(n - int(r.n_source)))
            if len(bad_ex) < 3:
                bad_ex.append(f"{r.organism}/{r.drug}:rebuild={n} frozen={int(r.n_source)}")
    ok = (n_bad == 0)
    return ok, (f"{site}: matched {n_ok}/{len(rows)} frozen n_source values; mismatches={n_bad} "
                f"max_abs_diff={worst}; examples={bad_ex}")


# ============================================================================
# Part 1 — modeling-table grain（MGB x {EC,KP,PA}）
# ============================================================================
def part1_site_inflation(site, T):
    """UNIT-01 / UNIT-02：episode 粒度与膨胀。"""
    for sym, t in T.items():
        X, unit = t["X"], t["unit"]
        n_rows = len(X)
        n_ep = int(X.index.nunique())
        L.check("UNIT-01", n_rows == n_ep, SRC_UNIT_DEF,
                f"{site}/{sym}: len(X)={n_rows} nunique(epid)=X.index.nunique()={n_ep} "
                f"inflation={n_rows/max(n_ep,1):.6f}",
                "one row per episode required; p6b:138 joins an epid-indexed meta that "
                "carries duplicate epid rows for multi-date episodes (cross-check: "
                f"unit.epid.nunique()={int(unit['epid'].nunique())})")
        sz = unit.groupby("epid").size()
        warn("UNIT-02",
             "definition R1 table: culture events = distinct (anon_id,cdate,spec) records "
             "with >=1 S/I/R AST; episodes = anchor-chain episodes",
             f"{site}/{sym}: events={n_rows} episodes={n_ep} "
             f"pct_episodes_with_exactly_1_date={100.0*(sz == 1).mean():.3f} "
             f"max_dates_per_episode={int(sz.max())}",
             "max_dates>1 => modeling table is NOT at episode grain"
             if int(sz.max()) > 1 else
             "all episodes single-date at this site => rows == episodes here (no inflation)")


def part1_index_hazard(site, T):
    """UNIT-04（current-state characterisation）：直接观测两个真实对象的 index 载体 ——
    p6b:138 的 X（epid 字符串索引，多日期 episode 造成重复行）与 p6b:144 的 STORE
    （reset_index(drop=True) 后为行号索引）。"""
    for sym, t in T.items():
        X, Xr, unit = t["X"], t["Xr"], t["unit"]
        epid_idx_is_positions = isinstance(X.index, pd.RangeIndex)
        store_idx_is_positions = isinstance(Xr.index, pd.RangeIndex)
        n_ep = int(unit["epid"].nunique())
        L.check("UNIT-04",
                (not epid_idx_is_positions) and store_idx_is_positions
                and int(X.index.nunique()) == n_ep,
                "characterisation (current-state): p6b:144 / p21a:276 apply "
                "reset_index(drop=True) -> index = row POSITIONS; definition R1/R2: "
                "the index consumed downstream must identify the episode",
                f"{site}/{sym}: X.index type={type(X.index).__name__} "
                f"nunique={int(X.index.nunique())} len={len(X)} (dup-epid payload); "
                f"Xr.index type={type(Xr.index).__name__} len={len(Xr)}; nunique(epid)={n_ep}",
                "mechanism of the reporting bug: set(Xb.index) yields positions, not epids; "
                "this line characterises the current plumbing only — the verdict is carried "
                "by UNIT-01/UNIT-06",
                warn_only=True)


def part1_reference_check(site, T):
    """UNIT-06R：两条独立 episode 划分（本模块 paraphrase + 冻结 SSoT）必须都与生产
    划分逐行等价（否则 UNIT-06 的 correct 目标无效）。"""
    for sym, t in T.items():
        Xr, epid_x, ref = t["Xr"], t["epid_x"], t["ref"]
        if pd.isna(pd.Series(epid_x)).any():
            L.emit("UNIT-06R", "FAIL", "definition R1 partition equivalence",
                   f"{site}/{sym}: epid could not be attached to all Xr rows", "")
            continue
        # 冻结 SSoT 必须显式传 group_keys：其默认 ("anon_id","organism","spec") 需要
        # organism 列，而 Xr 无该列（本模块在单 organism 表上运行）。
        sso = assign_episodes(Xr, group_keys=("anon_id", "spec"))
        same_ref = bool(np.array_equal(canon(ref), canon(epid_x)))
        same_sso = bool(np.array_equal(canon(sso), canon(epid_x)))
        same = same_ref and same_sso
        L.check("UNIT-06R", same,
                "definition R1: independent anchor-chain implementations (module paraphrase + "
                "frozen SSoT episode_definition.assign_episodes) must induce the same PARTITION "
                "as the production inline chain (relabelled canonical compare)",
                f"{site}/{sym}: n_episodes_ref={len(set(ref))} "
                f"n_episodes_prod={len(set(epid_x))} n_episodes_sso={len(set(sso))} "
                f"rowwise_partition_equal={same} (ref_vs_prod={same_ref}, sso_vs_prod={same_sso})",
                "if False, the reference is invalid and UNIT-06 must not be believed")


def part1_key_safety(T):
    """UNIT-03：epid 字符串跨 organism 的键安全性。"""
    pairs = [("EC", "KP"), ("EC", "PA"), ("KP", "PA")]
    obs, notes = [], []
    for a, b in pairs:
        if a not in T or b not in T:
            continue
        u, v = T[a]["unit"], T[b]["unit"]
        ea, eb = set(u["epid"].astype(str)), set(v["epid"].astype(str))
        m1 = dict(zip(zip(u["anon_id"].astype(str), u["spec"]), u["_blk"]))
        m2 = dict(zip(zip(v["anon_id"].astype(str), v["spec"]), v["_blk"]))
        shared = set(m1) & set(m2)
        blk_eq = sum(1 for k in shared if m1[k] == m2[k])
        obs.append(f"{a}-{b}: n_epid={len(ea)}/{len(eb)} inter={len(ea & eb)} "
                   f"shared_pat_spec={len(shared)} blk_equal={blk_eq}")
        notes.append(f"{a}-{b} ex_blk(pair)=" + str([(k[1], m1[k], m2[k]) for k in sorted(shared)[:2]]))
    for sym in ("EC", "KP", "PA"):
        if sym in T:
            u = T[sym]["unit"]
            obs.append(f"{sym}: composite(anon_id,spec,_blk) unique-per-episode="
                       f"{u[['anon_id','spec','_blk']].drop_duplicates().shape[0] == int(u['epid'].nunique())}")
    warn("UNIT-03",
         "definition R1: episode identity = (patient, organism, specimen, block); the epid "
         "string is anon_id+spec+blk with blk restarting at 0 per organism (p6b:134)",
         "MGB " + " ; ".join(obs),
         "PREDICTION OF NON-ZERO INTERSECTION **NOT** CONFIRMED on MGB: inter=0 for all pairs. "
         "Reason: blk is a per-organism GLOBAL episode counter (0..n_episodes-1), not a "
         "per-group local counter, so a cross-organism collision requires an exact "
         "(anon_id, spec, blk) coincidence. The structural hazard is nevertheless real: the "
         "string carries no organism tag, so the correct cross-organism key MUST be "
         "(organism, anon_id, spec, block). " + " ".join(notes))


def part1_conflicts(site, T):
    """UNIT-05：episode x drug 的 S/R、S/I、R/I 冲突计数（MGB）。"""
    agg = {}
    for sym, t in T.items():
        Xr, g = t["Xr"], t["g"]
        drugs = [d for d in t["bench_drugs"] if d in Xr.columns]
        miss = [d for d in t["bench_drugs"] if d not in Xr.columns]
        if not len(g):
            continue
        gg = g.assign(isS=(g["s"] == "S"), isI=(g["s"] == "I"), isR=(g["s"] == "R"))
        tab = gg.groupby(["epid", "abx"], sort=False)[["isS", "isI", "isR"]].max()
        tab["n"] = tab[["isS", "isI", "isR"]].sum(axis=1)
        tab = tab.reset_index()
        tab["SR"] = tab.isS & tab.isR & (~tab.isI)
        tab["SI"] = tab.isS & tab.isI & (~tab.isR)
        tab["RI"] = tab.isR & tab.isI & (~tab.isS)
        epx = pd.Series(t["epid_x"])
        any_conf = tab[tab.SR | tab.SI | tab.RI]
        tot = dict(n_ep=0, n_md=0, n_SR=0, n_SI=0, n_RI=0,
                   n_any_conf=int(any_conf["epid"].nunique()) if len(any_conf) else 0)
        by_drug = {d: sub for d, sub in tab.groupby("abx", sort=False)}
        n_mismatch = 0
        n_gap = 0
        for d in drugs:
            sub = by_drug.get(d)
            if sub is None or not len(sub):
                warn("UNIT-05", SRC_R3,
                     f"{site}/{sym}/{d}: n_episode=0 (drug column present but no AST row in g)", "")
                continue
            n_ep = int(len(sub))
            n_md = int((sub["n"] >= 2).sum())
            n_sr, n_si, n_ri = int(sub["SR"].sum()), int(sub["SI"].sum()), int(sub["RI"].sum())
            warn("UNIT-05", SRC_R3,
                 f"{site}/{sym}/{d}: n_episode={n_ep} n_multidate={n_md} n_SR_conflict={n_sr} "
                 f"n_SI_conflict={n_si} n_RI_conflict={n_ri}",
                 "resolved silently by pivot_table(aggfunc='first') = earliest date wins"
                 if (n_sr + n_si + n_ri) else "")
            tot["n_ep"] += n_ep
            tot["n_md"] += n_md
            tot["n_SR"] += n_sr
            tot["n_SI"] += n_si
            tot["n_RI"] += n_ri
            n_ep_x = int(pd.Series(epx[Xr[d].notna().values]).nunique())
            # ★ 2026-10-05（CS-03 amendment）：不变式由 `n_ep_x == n_ep` 改为 `0 <= n_ep_x <= n_ep`。
            #   CS-03 裁决「取 episode 内最早携带**直接 S/R** 的日期；I-only 的日期跳过」
            #   ⇒ X 只在存在直接 S/R 时非空 ⇒ **仅 I 或同日 S+R 被排除的 episode 必然从 X 侧消失**。
            #   这不是缺陷，是 CS-03 的**定义**；旧的等式编码的是 CS-03 之前的语义。
            #   差 = 该药「无直接 S/R」的 episode 数（结构必然 >= 0）。
            n_gap += (n_ep - n_ep_x)
            if not (0 <= n_ep_x <= n_ep):
                n_mismatch += 1
                warn("UNIT-05Q", "definition R3 consistency",
                     f"{site}/{sym}/{d}: g-side n_episode={n_ep} vs X-side distinct epid={n_ep_x}",
                     "X-side must not EXCEED g-side: X can only lose episodes (I-only / ambiguous "
                     "results are excluded by CS-03), never gain them")
        if miss:
            warn("UNIT-05", SRC_R3, f"{site}/{sym}: benchmark drugs absent from X.columns = {miss}", "")
        L.check("UNIT-05Q", n_mismatch == 0, "definition R3 consistency（CS-03 后）",
                f"{site}/{sym}: X-side **超过** g-side 的药数 = {n_mismatch}；"
                f"X 侧相对 g 侧共少的 episode-drug 对 = {n_gap}"
                f"（= 该药「无直接 S/R」的 episode 数：仅 I 或同日 S+R 被排除）",
                "0 expected for the EXCESS direction: CS-03 后 X 只可能**少**，不可能多。"
                "旧判据 `n_ep_x == n_ep` 编码的是 CS-03 之前的语义，已按冻结裁决更新")
        agg[sym] = tot
    for sym, tot in agg.items():
        warn("UNIT-05A", SRC_R3,
             f"{site}/{sym} TOTALS (summed over benchmark drugs): n_episode_obs={tot['n_ep']} "
             f"n_multidate={tot['n_md']} n_SR={tot['n_SR']} n_SI={tot['n_SI']} n_RI={tot['n_RI']}",
             f"episodes with >=1 conflict of any type (dedup by epid) = {tot['n_any_conf']}; "
             "these counts decide whether an outcome-aggregation rule must be frozen")


# ============================================================================
# Part 2 — reporting counts（p21a:305-318）
# ============================================================================
def part2_counts(site, T, BENCH, frozen_val):
    # ★ v1.4：STORE 必须喂给 p21a:305-318 **与生产同形**的输入。
    #   修正后 p21a:276 给 X 挂了 canonical `episode_key = organism|epid`（CS-12），
    #   而 p6b:144 的 Xr 已经把 epid 丢成行号；若仍直接喂 Xr，块内 `Xb["episode_key"]`
    #   会 KeyError，且即使用行号补也会**复现被修正的那个 bug**、使 UNIT-06 失真。
    #   本 harness 因此按 p21a:276 的同构方式重建该列，键取自 t["epid_x"]
    #   （= attach_epid 把生产 epid 贴回 Xr 行序的结果）。
    STORE = {(site, sym): t["Xr"].assign(
                 episode_key=sym + "|" + pd.Series(list(t["epid_x"]),
                                                   index=t["Xr"].index).astype(str))
             for sym, t in T.items()}
    out, src = prod_site_union(site, STORE, BENCH)
    prod = len(out["eps"])
    per, unions, corr_keys = {}, set(), set()
    for sym, t in T.items():
        Xr, ref = t["Xr"], t["ref"]
        drugs = list(BENCH[BENCH.organism == sym].drug)
        mask = Xr[drugs].isin(["S", "R"]).any(axis=1) if drugs else pd.Series(False, index=Xr.index)
        mv = mask.values
        per[sym] = int(mv.sum())
        unions |= set(Xr.index[mv])
        corr_keys |= {(sym, int(l)) for l in ref[mv]}
    correct = len(corr_keys)
    max_len_x = max(len(t["Xr"]) for t in T.values())
    top = max(per, key=per.get) if per else None
    L.check("UNIT-06", prod == correct,
            f"{SRC_R2}; frozen sidecar {os.path.basename(SITE_PROFILE_CSV)} {site}={frozen_val}",
            f"{site}: production_value={prod} correct_union={correct} ratio={prod/max(correct,1):.6f} "
            f"(frozen_csv={frozen_val})",
            "production eps OR= set(Xb.index) unions row POSITIONS (p21a:315 with X from "
            "p21a:276 reset_index); the definition unions (organism, episode) keys")
    frozen_anchor("UNIT-06F", prod == frozen_val,
                  f"characterisation (anchor): frozen artifact {os.path.basename(SITE_PROFILE_CSV)} "
                  f"n_unique_microbiological_episodes (derived by the unpatched pipeline)",
                  f"{site}: rebuild_production_value={prod} frozen={frozen_val} "
                  f"equal={prod == frozen_val}",
                  "reproducibility anchor: on the real tree the frozen site_profile number "
                  "reproduces from this rebuild; under EVS_PROD_DIR the corrected union "
                  "legitimately differs (characterisation, not a defect test)")
    L.check("UNIT-07", prod < correct,
            "characterisation/diagnosis (current buggy position-union; p21a:306-315): eps = "
            "set-union of the POSITION sets of Xb, so the value is bounded by max over "
            "organisms of len(X) and can never reach a genuine cross-organism episode union",
            f"{site}: production_value={prod}; n_union_positions={len(unions)}; len(Xb) per "
            f"organism={per}; len(X) per organism="
            f"{ {s: len(t['Xr']) for s, t in T.items()} }; max_len_X={max_len_x}; "
            f"closest_organism={top} (prod - len(Xb) = {prod - per.get(top, 0)}); "
            f"correct_union={correct}",
            "diagnosis: eps==position-union confirmed by construction (prod and len(unions) "
            "printed above); the prediction prod == max_sym len(Xb) holds only approximately "
            "because Xb positions are masked subsets of 0..len(X)-1",
            warn_only=True)
    return prod, per, correct


# ============================================================================
# main
# ============================================================================
def main():
    L.banner("UNIT — analytic unit and derived reporting counts")
    print(f"[loads] p6b:118-144 verbatim (contains BLOCK_STORE {L.block_source(L.BLOCK_STORE)}); "
          f"p21a:305-318 verbatim; setup = p6b:83-117 verbatim", flush=True)

    fb = pd.read_csv(FB_CSV)
    BENCH = fb[["organism", "drug"]].drop_duplicates()
    sp = pd.read_csv(SITE_PROFILE_CSV).set_index("site")
    frozen = {s: int(sp.loc[s, "n_unique_microbiological_episodes"]) for s in sp.index}
    print(f"[inputs] benchmark groups={len(BENCH)}; frozen site counts={frozen}", flush=True)

    summary = []
    for site in ("MGB", "SU", "UT"):
        try:
            t_site = time.time()
            d, dem, wd, src_setup = build_site(site)
            T = {}
            for sym, fname in FOCAL.items():
                out = build_org(d, dem, wd, site, sym)
                if out["X"] is None or out["Xr"] is None or not len(out["Xr"]):
                    print(f"[skip] {site}/{sym}: no rows", flush=True)
                    continue
                out["epid_x"] = attach_epid(out["Xr"], out["unit"])
                out["ref"] = ref_episode_labels(out["Xr"])
                out["bench_drugs"] = list(BENCH[BENCH.organism == sym].drug)
                T[sym] = out
            del d, dem, wd
            COVERED.append(site)
            print(f"[rebuild] {site}: organisms={list(T)} "
                  f"len(X)/nunique(epid)/len(Xb)="
                  f"{ {k: (len(v['Xr']), int(v['X'].index.nunique())) for k, v in T.items()} } "
                  f"({time.time()-t_site:.1f}s)", flush=True)

            _exp00 = (f"characterisation (anchor): frozen artifact "
                      f"p6b_label_free_signals_blinded.csv n_source (verbatim "
                      f"p6b:229 n_src=len(Xs), prep() mask drug in S/R) -- real-tree gate for {site}")
            try:
                _ok00, _obs00 = fidelity(site, T, fb)
                frozen_anchor("UNIT-00", _ok00, _exp00, _obs00,
                              "setup-replication fidelity gate: on the real tree a mismatch "
                              "invalidates this site's numbers; under EVS_PROD_DIR the frozen "
                              "universe legitimately changes (characterisation, not a defect test)")
            except Exception as e:
                L.emit("UNIT-00", "FAIL", _exp00, f"EXCEPTION {type(e).__name__}: {e}",
                       traceback.format_exc().splitlines()[-1])

            for fn in (part1_site_inflation, part1_index_hazard, part1_reference_check):
                try:
                    fn(site, T)
                except Exception as e:
                    L.emit("UNIT-EXC", "FAIL", f"test block {fn.__name__} on {site}",
                           f"EXCEPTION {type(e).__name__}: {e}",
                           traceback.format_exc().splitlines()[-1])

            if site == "MGB":
                for fn in (part1_key_safety, lambda T_: part1_conflicts(site, T_)):
                    try:
                        fn(T)
                    except Exception as e:
                        L.emit("UNIT-EXC", "FAIL", f"test block {fn} on {site}",
                               f"EXCEPTION {type(e).__name__}: {e}",
                               traceback.format_exc().splitlines()[-1])

            try:
                prod, per, correct = part2_counts(site, T, BENCH, frozen.get(site, -1))
                summary.append(dict(site=site, prod=prod, per_organism_Xb=per,
                                    correct_union=correct, frozen=frozen.get(site, -1)))
            except Exception as e:
                L.emit("UNIT-EXC", "FAIL", "part2 for " + site,
                       f"EXCEPTION {type(e).__name__}: {e}",
                       traceback.format_exc().splitlines()[-1])
            del T
        except Exception as e:
            FAILED.append(site)
            L.emit("UNIT-EXC", "FAIL", f"site rebuild {site} (p6b:83-144 verbatim)",
                   f"EXCEPTION {type(e).__name__}: {e}",
                   traceback.format_exc().splitlines()[-1])

    print("", flush=True)
    warn("UNIT-COV", "coverage report (never claim more than achieved)",
         f"rebuilt sites={COVERED}; failed sites={FAILED}; summary={summary}",
         "all-site: UNIT-00/01/02/04/06/06F/06R/07; "
         "MGB-only: UNIT-03/UNIT-05/UNIT-05A/UNIT-05Q")
    print(f"[done] {time.time()-T0:.1f}s", flush=True)


if __name__ == "__main__":
    main()
