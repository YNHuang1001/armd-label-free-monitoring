# -*- coding: utf-8 -*-
r"""Protocol v1.3 §10 —— Phase B: 样本量敏感性曲线（全部 L1 / L2 / L3）

§10.1（v1.3 冻结）实现：
  - 5 个 target N 档：100 / 250 / 500 / 1000 / full，每档 B = 200 次重复
  - source reference 全部固定（见 p10a 缓存），**只反复抽 target X**
  - B=200 不减；不删任何方法；所有方法同一 B
  - 不读取 target Y 参与任何 signal 计算；Y 仅用于最终 §7.2 评价

每个 repeat = 一个"部署场景"：202 个 task 各抽一份大小为 N 的 target 子样本，
在该场景内计算全部 signal，再对 202 task 算 harmful 检测指标（§7.2）。
跨 200 个 repeat 汇总 mean / sd / 95% interval。

输出：report_p10b_subsample.txt / p10b_subsample_curves.csv / p10b_raw.csv
"""
import os, sys, time, glob, pickle, warnings
import numpy as np, pandas as pd
warnings.filterwarnings("ignore")
from multiprocessing import Pool
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score, average_precision_score
from scipy.spatial.distance import pdist, cdist

BASE = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))); OUT = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(OUT, "p10_cache")
SEED = 20261005
MMD_M, MMD_K, NBIN = 1000, 5, 10
SUDO_REP, SUDO_CAP, SUDO_FRAC = 5, 200, 0.5
NLEVELS = [100, 250, 500, 1000, 0]      # 0 = full
BREP = 200


def _prov_signature():
    """缓存**产出代码签名** = sha256(p10a 源码 + p10b 源码)。与 p10a 侧**同实现**，必须一致。

    ★ 2026-10-06 新增：原先 p10b 直接 `pickle.load`，对缓存的来源**零校验**。"""
    import hashlib as _h
    h = _h.sha256()
    here = os.path.dirname(os.path.abspath(__file__))
    for fn in ("p10a_prepare_cache.py", "p10b_subsample_curves.py"):
        with open(os.path.join(here, fn), "rb") as f:
            h.update(f.read())
    return h.hexdigest()
FEATS = ["age", "set_IP", "set_OP", "set_ER", "sp_URINE", "sp_BLOOD",
         "prior_any", "prior_cnt", "t_15-30", "t_31-90", "t_91-180", "t_181-365",
         "past_S", "past_I", "past_R"]
CONT = {"age", "prior_cnt"}

def ent(pv):
    pv = np.clip(pv, 1e-12, 1 - 1e-12); return -(pv * np.log(pv) + (1 - pv) * np.log(1 - pv))

def wconcordance(score, a, b):
    o = np.argsort(score, kind="mergesort")
    s, aa, bb = score[o], a[o], b[o]
    uniq, inv = np.unique(s, return_inverse=True)
    A = np.bincount(inv, weights=aa); B = np.bincount(inv, weights=bb)
    cumB = np.concatenate([[0.0], np.cumsum(B)[:-1]])
    return float(np.sum(A * cumB) + 0.5 * (np.sum(A * B) - np.sum(aa * bb)))

def psijs_from_hist(r, a):
    r = r.astype(float); a = a.astype(float)
    if r.sum() == 0 or a.sum() == 0: return np.nan, np.nan
    rp = r / max(r.sum(), 1); ap = a / max(a.sum(), 1)
    eps = 1e-6
    rpc = np.clip(rp, eps, 1); apc = np.clip(ap, eps, 1)
    psi = float(np.sum((apc - rpc) * np.log(apc / rpc)))
    M = 0.5 * (rp + ap)
    def kl(x, y):
        m = x > 0; return float(np.sum(x[m] * np.log2(x[m] / np.maximum(y[m], 1e-12))))
    return psi, 0.5 * kl(rp, M) + 0.5 * kl(ap, M)

# ---------------- worker ----------------
def run_task(tid):
    A = pickle.load(open(os.path.join(CACHE, f"{tid}.pkl"), "rb"))
    # ★ 2026-10-06：缓存来源校验。签名不符 = 该缓存由**不同版本**的 p10a/p10b 产出，
    #   直接**硬失败**而不是用旧缓存继续算（原先无任何校验，旧缓存会被静默当成新的用）。
    _want = _prov_signature()
    if A.get("_prov") != _want:
        raise SystemExit(
            f"[cache] {tid}.pkl 的产出代码签名不符（缓存 {str(A.get('_prov'))[:12]}… "
            f"!= 当前 {_want[:12]}…）—— 缓存由不同版本的 p10a/p10b 产生。\n"
            f"        必须重跑 p10a_prepare_cache.py 重建缓存。**拒绝使用来源不明的缓存。**")
    rows = []
    n_t = A["n_target"]
    for N in NLEVELS:
        n = n_t if N == 0 else min(N, n_t)
        for r in range(BREP):
            rng = np.random.default_rng((hash_seed(tid) + 7919 * r + 1) % (2**63))
            idx = rng.choice(n_t, n, replace=False)
            idx.sort()
            Dt_std_s = A["Dt_std"][idx]; pt_s = A["pt_full"][idx]; m_t_s = A["m_t_full"][idx]
            # level = 档位标签（"full"/"100"/...）；N = 实际抽样规模（full 档 = 该 task 真实 target 大小）
            rec = {"task_id": tid, "level": ("full" if N == 0 else str(N)),
                   "N": (n_t if N == 0 else N), "n_draw": n, "rep": r}

            # ===== L1a domain classifier（固定 5,000 source reference + target 子样本）=====
            Xd = np.vstack([A["Ds_ref_std"], Dt_std_s])
            yd = np.r_[np.zeros(len(A["Ds_ref_std"])), np.ones(n)]
            oof = np.zeros(len(yd))
            for tr, te in StratifiedKFold(5, shuffle=True, random_state=SEED).split(Xd, yd):
                dm = LogisticRegression(max_iter=2000, class_weight="balanced", random_state=SEED)
                dm.fit(Xd[tr], yd[tr]); oof[te] = dm.predict_proba(Xd[te])[:, 1]
            rec["domain_auc"] = float(roc_auc_score(yd, oof)) if len(np.unique(yd)) > 1 else np.nan

            # ===== L1b MMD（source reference 固定 m=1000，Kaa 预计算）=====
            vals = []
            for k in range(MMD_K):
                b = Dt_std_s[rng.choice(n, min(MMD_M, n), replace=False)]
                m_ = len(b)
                Kbb = np.exp(-pdist(b, "sqeuclidean") / A["bw2"])
                Kab = np.exp(-cdist(A["mmd_src"][k], b, "sqeuclidean") / A["bw2"])
                vals.append(A["Kaa"][k] + 2 * Kbb.sum() / (m_ * (m_ - 1)) - 2 * Kab.mean())
            rec["mmd"] = float(np.mean(vals))

            # ===== L1c PSI / JS（bin 与 source 参考直方图已固定）=====
            psis, jss = [], []
            for j, f in enumerate(FEATS):
                # bin 与 source 参考直方图在 p10a 已用 source-train ref + full target 固定；
                # 子样本只需 bincount 其预先算好的 bin 索引（与 p6b 的 full-target 情形逐位一致）
                a = np.bincount(A["tgt_binidx"][j][idx], minlength=len(A["psi_refhist"][j]))
                pv, jv = psijs_from_hist(A["psi_refhist"][j], a)
                psis.append(pv); jss.append(jv)
            pa, ja = np.array(psis, float), np.array(jss, float)
            rec["psi_mean"] = float(np.nanmean(pa)); rec["psi_max"] = float(np.nanmax(pa))
            rec["js_mean"] = float(np.nanmean(ja)); rec["js_max"] = float(np.nanmax(ja))

            # ===== L2 prediction-space =====
            rec["pred_mean_shift"] = float(pt_s.mean() - A["ps"].mean())
            rec["entropy_mean_shift"] = float(ent(pt_s).mean() - A["ps_ent"].mean())
            a_ent = np.bincount(A["jsent"]["tgtidx"][idx], minlength=len(A["jsent"]["refhist"]))
            rec["entropy_js"] = psijs_from_hist(A["jsent"]["refhist"], a_ent)[1]
            a_conf = np.bincount(A["jsconf"]["tgtidx"][idx], minlength=len(A["jsconf"]["refhist"]))
            rec["conf_dist_js"] = psijs_from_hist(A["jsconf"]["refhist"], a_conf)[1]
            rec["pred_pos_prev_shift"] = float((pt_s > 0.5).mean() - (A["ps"] > 0.5).mean())

            # ===== L2 SUDO-AURCC（source train/val 原样；只变 target 子样本）=====
            qs = np.quantile(pt_s, np.linspace(0, 1, NBIN + 1)); qs[0], qs[-1] = -np.inf, np.inf
            bid = np.digitize(pt_s, qs[1:-1])
            src_by_class = {0: np.where(A["ytr"] == 0)[0], 1: np.where(A["ytr"] == 1)[0]}
            sbins, frac, conf = [], [], []
            for b in range(NBIN):
                ib_ = np.where(bid == b)[0]
                if len(ib_) < 10: continue
                yb = 1 if pt_s[ib_].mean() >= 0.5 else 0
                pool = src_by_class[1 - yb]
                if len(pool) < 10: continue
                aucs = []
                for rr in range(SUDO_REP):
                    nb = min(int(SUDO_FRAC * len(ib_)), SUDO_CAP, len(ib_), len(pool))
                    xb = rng.choice(ib_, nb, replace=False); sb_ = rng.choice(pool, nb, replace=False)
                    Xa = np.vstack([Dt_std_s[xb], A["Dtr_std"][sb_]])
                    ya = np.r_[np.ones(nb), np.zeros(nb)]
                    clf = LogisticRegression(max_iter=2000, random_state=SEED + rr).fit(Xa, ya)
                    yv = (A["yva"] == yb).astype(int)
                    if len(np.unique(yv)) < 2: continue
                    aucs.append(roc_auc_score(yv, clf.decision_function(A["Dva_std"])))
                if len(aucs) < 3: continue
                sbins.append(2 * (float(np.mean(aucs)) - 0.5))
                frac.append(len(ib_) / n); conf.append(abs(pt_s[ib_].mean() - 0.5))
            if sbins:
                order = np.argsort(-np.array(conf))
                fr = np.array(frac)[order]; rb = np.abs(np.array(sbins))[order]
                comp = np.cumsum(fr) / np.sum(fr)
                rel = np.cumsum(rb * fr) / np.cumsum(fr)
                rec["SUDO_AURCC"] = float(np.trapezoid(rel, comp))
                rec["SUDO_nbin"] = len(sbins)
            else:
                rec["SUDO_AURCC"] = np.nan; rec["SUDO_nbin"] = 0

            # ===== L3 =====
            # membership model：固定 source reference + target 子样本
            dm2 = LogisticRegression(max_iter=2000, class_weight="balanced", random_state=SEED)
            dm2.fit(np.vstack([A["Ds_ref_std"], Dt_std_s]), np.r_[np.zeros(len(A["Ds_ref_std"])), np.ones(n)])
            w_te = np.clip(dm2.predict_proba(A["Dte_std"])[:, 1], 1e-6, 1 - 1e-6)
            w_te = w_te / (1 - w_te)
            rec["w_P95"] = float(np.percentile(w_te, 95)); rec["w_max"] = float(w_te.max())
            rec["w_ESS"] = float((w_te.sum() ** 2) / (w_te ** 2).sum())
            s_e, y_e = A["s_e"], A["yte"]
            num = wconcordance(s_e, w_te * (y_e == 1), w_te * (y_e == 0))
            den = (w_te * (y_e == 1)).sum() * (w_te * (y_e == 0)).sum()
            rec["AUROC_IOW"] = float(num / den) if den > 0 else np.nan
            m_e = A["m_e"]
            sc_t = A["st_full"][idx]            # 与 p6c 的 decision_function 完全一致
            numo = wconcordance(sc_t, m_t_s, 1 - m_t_s); deno = m_t_s.sum() * (1 - m_t_s).sum() - (m_t_s * (1 - m_t_s)).sum()
            rec["AUROC_OM"] = float(numo / deno) if deno > 0 else np.nan
            yt1 = (y_e == 1).astype(float); yt0 = (y_e == 0).astype(float); w1 = w_te * yt1; w0 = w_te * yt0; p1 = w_te * m_e; p0 = w_te * (1 - m_e)
            term_s = wconcordance(s_e, w1, w0) - wconcordance(s_e, p1, p0)
            term_t = wconcordance(sc_t, m_t_s, 1 - m_t_s)
            n_s = len(y_e)
            N_DR = term_s + term_t
            d_s = (w1.sum() * w0.sum() - (w1 * w0).sum()) - (p1.sum() * p0.sum() - (p1 * p0).sum())
            d_t = m_t_s.sum() * (1 - m_t_s).sum() - (m_t_s * (1 - m_t_s)).sum()
            rec["AUROC_DR"] = float(N_DR / (d_s + d_t)) if (d_s + d_t) > 0 else np.nan
            rec["status"] = "OK"
            rows.append(rec)
    return rows

def hash_seed(tid):
    import hashlib
    return int(hashlib.md5(("p10b|" + tid).encode()).hexdigest()[:12], 16)

def init():
    pass

if __name__ == "__main__":
    t0 = time.time()
    rep = open(os.path.join(OUT, "report_p10b_subsample.txt"), "w", encoding="utf-8")
    def p(*a):
        line = " ".join(str(x) for x in a); print(line, flush=True); rep.write(line + "\n"); rep.flush()

    tids = sorted(f[:-4] for f in os.listdir(CACHE) if f.endswith(".pkl"))
    p("=" * 116)
    p("[Protocol v1.3 §10.1] target 样本量敏感性 —— 全部 L1 / L2 / L3")
    p("=" * 116)
    p(f"  task = {len(tids)}   N 档 = {[('full' if v==0 else v) for v in NLEVELS]}   B = {BREP}")
    p(f"  source reference 全部固定（domain/membership 5,000；MMD m=1,000；PSI/JS/pred-space/OM full）")
    p(f"  只反复抽 target X。B=200 不减，不删任何方法。")

    NW = int(os.environ.get("P10_WORKERS", "24"))
    p(f"  workers = {NW}")
    RAWP = os.path.join(OUT, "p10b_raw.csv")
    # ---- resume 保护：若已有完整 raw 结果则直接复用，绝不重跑 2 小时计算 ----
    # ★★ 2026-10-06 修正：原先只比对**行数与 task 数** —— 形状对得上的**旧结果**会被静默复用。
    #   实测风险：补丁 A3 改了 `Kaa` 的含义、CS-01 v2 改了药物身份之后，一份形状正确的
    #   修正前 `p10b_raw.csv` 会被当成有效结果跳过全部计算。
    #   现增加**依赖检查**（真实依赖，不是挂钟启发式）：raw 由 p10_cache 派生，
    #   任何缓存文件比 raw 新 ⇒ raw 已过期 ⇒ 必须重算。
    #   ★ 判据有两条，缺一不可（只查 mtime 不够：实测 p10b_raw.csv(09-30 10:54) 比
    #     缓存(09-30 08:20) **更新**，所以"没有缓存比它新" ⇒ 看起来依赖没变，实则两者都是修正前的）：
    #       ① 整套缓存必须带**当前**产出代码签名（由 p10a 的完成标记 _prov.txt 回答）；
    #       ② 不能有缓存文件比 raw 新。
    _need = len(tids) * len(NLEVELS) * BREP
    _stale = None
    _mk = os.path.join(CACHE, "_prov.txt")
    if not os.path.exists(_mk):
        _stale = f"缓存完成标记 {os.path.basename(_mk)} 不存在 ⇒ 缓存来源不明（或 p10a 未跑完）"
    else:
        _have = open(_mk, encoding="utf-8").read().strip()
        if _have != _prov_signature():
            _stale = (f"缓存完成标记 {_have[:12]}… != 当前代码签名 {_prov_signature()[:12]}… "
                      f"⇒ 缓存由**不同版本**的 p10a/p10b 产出")
    if os.path.exists(RAWP) and _stale is None:
        _raw_m = os.path.getmtime(RAWP)
        _newer = [f for f in os.listdir(CACHE)
                  if f.endswith(".pkl") and os.path.getmtime(os.path.join(CACHE, f)) > _raw_m]
        if _newer:
            _stale = (f"{len(_newer)}/{len(tids)} 个缓存文件比 {os.path.basename(RAWP)} 新"
                      f"（例：{_newer[0]}）⇒ 该 raw 由**更早的缓存**派生，已过期")
    if os.path.exists(RAWP) and _stale is None:
        _r = pd.read_csv(RAWP, usecols=["task_id", "level", "rep"])
        if len(_r) == _need and _r.task_id.nunique() == len(tids):
            p(f"  [resume] 复用已有 {RAWP}（{len(_r)} 行 = {len(tids)} x {len(NLEVELS)} x {BREP}），"
              f"跳过计算；依赖检查通过（无缓存比它新）")
            RAW = pd.read_csv(RAWP)
        else:
            p(f"  [resume] 已有文件不完整（{len(_r)} != {_need}），重算")
            RAW = None
    else:
        if _stale:
            p(f"  [resume] **拒绝复用**：{_stale} ⇒ 重算")
        RAW = None
    if RAW is None:
        allrows = []
        with Pool(NW) as pool:
            for i, rows in enumerate(pool.imap_unordered(run_task, tids, chunksize=1), 1):
                allrows.extend(rows)
                if i % 10 == 0: p(f"    ... {i}/{len(tids)} tasks  [{time.time()-t0:.0f}s]")
        RAW = pd.DataFrame(allrows)
        RAW.to_csv(RAWP, index=False)
        p(f"\n[写出] p10b_raw.csv  行数 = {len(RAW)}  [{time.time()-t0:.0f}s]")

    # ---------------- 聚合 ----------------
    O = pd.read_csv(os.path.join(OUT, "p6d_outcome.csv"))[
        ["task_id", "organism", "drug", "d_auroc", "harmful", "stable"]]
    RAW = RAW.merge(O, on="task_id", how="left")
    # level 列已在 raw 行内直接写出（v1 bug：用 N==0 筛 full 档 → 空集 → percentile 崩）
    # 注意：full 档的标签是字符串 "full"，不是 str(0)；曾因写成 "0" 导致整列变 NaN
    RAW["level"] = pd.Categorical(RAW["level"], categories=["100", "250", "500", "1000", "full"], ordered=True)
    assert RAW["level"].isna().sum() == 0, "level 标签缺失"
    _c = RAW.groupby(["task_id", "level"], observed=True).size()
    assert (_c == BREP).all(), f"有 task 档位行数 != {BREP}"
    assert len(RAW) == RAW.task_id.nunique() * len(NLEVELS) * BREP, "行数不等于 task x 档位 x B"
    p(f"  档位自洽性检查：通过（每 task 每档 {BREP} 行）")
    ycol = "harmful"

    METHODS = [
        ("domain_auc", "L1", lambda x: x["domain_auc"]),
        ("mmd", "L1", lambda x: x["mmd"]),
        ("psi_mean", "L1", lambda x: x["psi_mean"]),
        ("psi_max", "L1", lambda x: x["psi_max"]),
        ("js_mean", "L1", lambda x: x["js_mean"]),
        ("js_max", "L1", lambda x: x["js_max"]),
        ("pred_mean_shift", "L2", lambda x: x["pred_mean_shift"].abs()),
        ("entropy_mean_shift", "L2", lambda x: x["entropy_mean_shift"].abs()),
        ("entropy_js", "L2", lambda x: x["entropy_js"]),
        ("conf_dist_js", "L2", lambda x: x["conf_dist_js"]),
        ("pred_pos_prev_shift", "L2", lambda x: x["pred_pos_prev_shift"].abs()),
        ("SUDO_AURCC", "L2", lambda x: -x["SUDO_AURCC"]),
        ("AUROC_DR", "L3-primary", lambda x: -x["AUROC_DR"]),
        ("AUROC_IOW", "L3-sens", lambda x: -x["AUROC_IOW"]),
        ("AUROC_OM", "L3-sens", lambda x: -x["AUROC_OM"]),
    ]

    def sens_at_spec(y, s, spec=0.80):
        m = np.isfinite(s)
        y, s = y[m], s[m]
        if len(np.unique(y)) < 2: return np.nan
        thr = min([v for v in np.unique(s) if (s[y == 0] < v).mean() >= spec] or [np.nextafter(s.max() if s.size else 0.0, np.inf)])
        pred = s >= thr
        tp = ((pred == 1) & (y == 1)).sum(); fn = ((pred == 0) & (y == 1)).sum()
        return tp / max(tp + fn, 1) if (len(np.unique(y)) > 1 and (y == 0).sum() > 0 and (y == 1).sum() > 0) else np.nan

    out = []
    for N in NLEVELS:
        lv = "full" if N == 0 else str(N)
        sub = RAW[RAW["level"] == lv]
        p(f"\n  ---- N = {lv}  (覆盖 task {sub.task_id.nunique()}, repeat {sub.rep.nunique()}) ----")
        for name, role, fn in METHODS:
            aucs, aps, sens, cov = [], [], [], []
            for r, g in sub.groupby("rep"):
                s = fn(g).values.astype(float); y = g[ycol].values.astype(int)
                m = np.isfinite(s)
                cov.append(m.sum())
                if len(np.unique(y[m])) < 2: continue
                aucs.append(roc_auc_score(y[m], s[m]))
                aps.append(average_precision_score(y[m], s[m]))
                sens.append(sens_at_spec(y, s))
            aucs = np.array(aucs, float)
            out.append(dict(method=name, role=role, N=lv,
                            n_repeat=len(aucs), coverage_mean=float(np.mean(cov)),
                            auroc_mean=float(np.mean(aucs)), auroc_sd=float(np.std(aucs)),
                            auroc_lo=float(np.percentile(aucs, 2.5)), auroc_hi=float(np.percentile(aucs, 97.5)),
                            auprc_mean=float(np.mean(aps)), sens80_mean=float(np.nanmean(sens))))
    CUR = pd.DataFrame(out)
    CUR.to_csv(os.path.join(OUT, "p10b_subsample_curves.csv"), index=False)

    p("\n" + "=" * 116)
    p("【主表】harmful 检测 AUROC 随可获得的 target 样本量变化  (mean [2.5%, 97.5%] over 200 repeats)")
    p("=" * 116)
    cols = ["full", "1000", "500", "250", "100"]
    hdr = f"{'method':<22}{'role':<12}" + "".join(f"{c:>18}" for c in cols) + f"{'Drift':>9}"
    p(hdr)
    for name, role, _ in METHODS:
        cells = []
        vals = []
        for c in cols:
            row = CUR[(CUR.method == name) & (CUR.N == c)]
            if len(row) == 0: cells.append(" " * 18); vals.append(np.nan); continue
            v = row.auroc_mean.iat[0]; vals.append(v)
            cells.append(f"{v:.3f}[{row.auroc_lo.iat[0]:.3f},{row.auroc_hi.iat[0]:.3f}]")
        drift = (np.nanmax(vals) - np.nanmin(vals)) if np.isfinite(vals).any() else np.nan
        p(f"{name:<22}{role:<12}" + "".join(f"{c:>18}" for c in cells) + f"{drift:>9.3f}")

    p("\n  AUPRC / sens@spec0.80 随 N 变化")
    for name, role, _ in METHODS:
        a = [CUR[(CUR.method == name) & (CUR.N == c)].auprc_mean.iat[0] if len(CUR[(CUR.method == name) & (CUR.N == c)]) else np.nan for c in cols]
        s = [CUR[(CUR.method == name) & (CUR.N == c)].sens80_mean.iat[0] if len(CUR[(CUR.method == name) & (CUR.N == c)]) else np.nan for c in cols]
        p(f"    {name:<22} AUPRC " + " ".join(f"{v:.3f}" for v in a) + "   sens80 " + " ".join(f"{v:.3f}" for v in s))

    # ★ 2026-10-06：原写死 "满额 202"。改为不写数字的表述 —— 该满额随 eligible 宇宙变化，
    #   写死则重跑后必错。此处不取数是因为本函数作用域内没有权威的 task 总数变量，
    #   与其猜一个不如不给数字（给错数字比不给更糟）。
    p("\n   Coverage（每 repeat 有限 signal 的 task 数；满额 = 本期 eligible task 数）")
    for name, role, _ in METHODS:
        c = [CUR[(CUR.method == name) & (CUR.N == cc)].coverage_mean.iat[0] if len(CUR[(CUR.method == name) & (CUR.N == cc)]) else np.nan for cc in cols]
        p(f"    {name:<22} " + " ".join(f"{v:7.1f}" for v in c))
    rep.close()
    print("done")
