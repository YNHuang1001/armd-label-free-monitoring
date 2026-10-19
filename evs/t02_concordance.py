# -*- coding: utf-8 -*-
r"""t02_concordance.py — family CONC — `wconcordance` 与 AUROC 语义的独立验证

生产被测对象（逐字抽取 exec，不复制、不修改）
--------------------------------------------
    WC_P6C  = L.load_def(L.P6C,  "wconcordance")   # p6c_01_L3_blinded.py:84-92
    WC_P10B = L.load_def(L.P10B, "wconcordance")   # p10b_subsample_curves.py:40-47
两者都声称同一泛函（生产 docstring p6c:85 逐字）：
    C(s, a, b) = Σ_{i≠j} a_i b_j 1[s_i > s_j]        （严格大于；并列 pair 记 0）

参考实现写自哪条定义 —— 不是抄生产代码
---------------------------------------
[R1] 一致性泛函（并列记 0.5）：显式二重循环逐对 (i, j) 判定 i != j；s_i > s_j 记 1，
     s_i = s_j 记 0.5，s_i < s_j 记 0。
     来源：与 [R2] 同一并列半分约定（Mann–Whitney 1947；Hanley & McNeil 1982），
     并由 CONC-06/06C 的 sklearn 独立佐证；生产把它作 AUROC 分子用（[PIPE]），
     故并列 pair 记 0.5——冻结 docstring p6c:85 的「严格 >」写法正是缺陷 (b)。
     生产用 argsort + bincount 向量化，本参考用最朴素的逐对枚举，实现路径完全不同。
[R2] AUROC（并列记 0.5 分，加权推广）：
        AUROC = [ Σ_{i:y=1} Σ_{j:y=0} w_i w_j ( 1[s_i>s_j] + 0.5·1[s_i=s_j] ) ]
                / [ (Σ_{i:y=1} w_i) · (Σ_{j:y=0} w_j) ]
     来源：Mann–Whitney U 统计量的归一化形式（Mann & Whitney 1947, Ann. Math. Stat.
     18:50–60）；Hanley & McNeil 1982（Radiology 143:29–36）关于并列半个记分的表述；
     scikit-learn 文档明示 roc_auc_score 在样本并列时给 0.5。
     R2 另由本次的 sklearn.metrics.roc_auc_score 独立交叉验证（CONC-06/06C）。
[PIPE] "ratio 形式"（生产用途）：num/den，a = w·1[y=1]，b = w·1[y=0]，
     den = (Σ w·1[y=1])(Σ w·1[y=0])，即 p6c:271-273 的 IOW 形式（对照 p6c 头部冻结
     公式块 [IOW]：Ahat_IOW = ΣΣ w_i w_j 1[s_i>s_j] / [(Σ_{y=1}w)(Σ_{y=0}w)]）。
     冻结公式块本身也写死了严格 `>`；本族由此同时暴露"公式块 headline 的严格 > 解释"
     与"该 ratio 被当作 AUROC 使用（要求并列 0.5 记分）"之间的语义缺口。

目标缺陷 (b)：生产用严格 `>`，并列 pair 记 0 分；标准 AUROC 记 0.5 分。
预期亮红：CONC-01 / CONC-04 / CONC-04B / CONC-05B / CONC-06B。
（并列语义断言全部按「并列记 0.5」书写；生产修复后它们应转绿，不得再按严格 > 判为回归。）

只读约束：p10_cache 只读，绝不写入。
"""
import os, sys, gc, glob, pickle, traceback
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import evs_load as L

RNG = np.random.default_rng(20261004)
CACHE_DIR = os.path.join(L.OUT, "p10_cache")      # <L.OUT>/p10_cache （只读）
TOL = 1e-9

# ---------------------------------------------------------------- 生产函数装载
LOAD_ERR = ""
WC_P6C = WC_P10B = None
try:
    WC_P6C = L.load_def(L.P6C, "wconcordance")
    WC_P10B = L.load_def(L.P10B, "wconcordance")
except Exception as e:                            # 装载失败不能整族挂掉
    LOAD_ERR = f"{type(e).__name__}: {e}"


# ---------------------------------------------------------------- 参考实现（R1/R2）
def ref_concordance_loop(score, a, b):
    r"""[R1] Σ_{i≠j} a_i b_j (1[s_i > s_j] + 0.5·1[s_i = s_j]) —— 显式二重循环逐对判定。

    来源：Mann–Whitney U 统计量的并列半分（mid-rank）约定（Mann & Whitney 1947；
    Hanley & McNeil 1982；与 [R2] 同一并列规则，由 sklearn 于 CONC-06/06C 独立佐证）。
    该泛函在管线中充当 AUROC 分子（[PIPE]），故并列 pair 记 0.5 而非 0。未参考生产实现。
    """
    s = np.asarray(score, dtype=float)
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    n = s.shape[0]
    tot = 0.0
    for i in range(n):
        for j in range(n):
            if i != j:
                if s[i] > s[j]:
                    tot += a[i] * b[j]
                elif s[i] == s[j]:
                    tot += 0.5 * a[i] * b[j]
    return float(tot)


def ref_auroc_ties(score, y, w=None):
    r"""[R2] 带并列 0.5 记分的加权 AUROC
           = Σ_{i:y=1}Σ_{j:y=0} w_i w_j (1[s_i>s_j] + 0.5·1[s_i=s_j]) / (W1·W0)。

    来源：Mann–Whitney U / Hanley–McNeil 的 AUROC 定义（并列记半分），
    由 sklearn.metrics.roc_auc_score 独立交叉验证。
    """
    s = np.asarray(score, dtype=float)
    y = np.asarray(y)
    if w is None:
        w = np.ones(s.shape[0], dtype=float)
    w = np.asarray(w, dtype=float)
    ip = np.where(y == 1)[0]
    ineg = np.where(y == 0)[0]
    W1 = float(w[ip].sum())
    W0 = float(w[ineg].sum())
    if W1 <= 0 or W0 <= 0:
        return float("nan"), 0.0
    num = 0.0
    for i in ip:
        si = s[i]
        wi = w[i]
        for j in ineg:
            sj = s[j]
            if si > sj:
                num += wi * w[j]
            elif si == sj:
                num += 0.5 * wi * w[j]
    return float(num / (W1 * W0)), W1 * W0


def prod_ratio(fn, score, w, y):
    """[PIPE] 生产用途的 ratio：num/den，den 与 IOW 形式（p6c:271-273）一致。"""
    s = np.asarray(score, dtype=float)
    y = np.asarray(y)
    w = np.asarray(w, dtype=float)
    a = w * (y == 1)
    b = w * (y == 0)
    num = fn(s, a, b)
    den = float(a.sum() * b.sum())
    return (float(num / den) if den > 0 else float("nan")), float(num), den


def _clean(x):
    return str(x).replace("|", "/")


def _guard_warn(tid, fn, expected, note=""):
    """L.guard 的 warn_only 变体（诊断项用）：异常也记 WARN，内部仍走 L.check。"""
    try:
        r = fn()
        ok, observed = r if isinstance(r, tuple) else (bool(r), str(r))
        return L.check(tid, ok, expected, _clean(observed), _clean(note), warn_only=True)
    except Exception as e:
        return L.emit(tid, "WARN", expected, _clean(f"EXCEPTION {type(e).__name__}: {e}"),
                      _clean(note + " | " + traceback.format_exc().splitlines()[-1]))


# ================================================================ CONC-01
def t_conc01():
    """Brute force agreement（功能层）：生产两副本 vs R1 显式二重循环（并列记 0.5），含大量精确并列。"""
    if WC_P6C is None or WC_P10B is None:
        raise RuntimeError(f"生产函数装载失败: {LOAD_ERR}")
    worst_c = worst_b = worst_cb = 0.0
    shapes = []
    for n, grid in ((50, 7), (91, 11), (137, 5), (200, 13)):
        s = RNG.integers(0, grid, size=n).astype(float)
        s[RNG.random(n) < 0.35] = 0.0                    # 人为制造更多精确并列
        a = RNG.random(n)
        b = RNG.random(n)
        r = ref_concordance_loop(s, a, b)
        vc = WC_P6C(s, a, b)
        vb = WC_P10B(s, a, b)
        sc = max(1.0, abs(r))
        worst_c = max(worst_c, abs(vc - r) / sc)
        worst_b = max(worst_b, abs(vb - r) / sc)
        worst_cb = max(worst_cb, abs(vc - vb) / sc)
        shapes.append(f"n={n}/uniq={len(np.unique(s))}")
    ok = (worst_c <= TOL) and (worst_b <= TOL) and (worst_cb <= TOL)
    observed = ("max_rel_err(p6c_vs_R1)=%.3e max_rel_err(p10b_vs_R1)=%.3e p6c_vs_p10b=%.3e ; %s"
                % (worst_c, worst_b, worst_cb, "; ".join(shapes)))
    return ok, observed


# ================================================================ CONC-02 / 03
def _sep_case(perfect):
    n1, n0 = 40, 60
    if perfect:
        s = np.r_[np.arange(n0) + 0.0, np.arange(n1) + 100.0]   # 正类全部严格高于负类
    else:
        s = np.r_[np.arange(n0) + 100.0, np.arange(n1) + 0.0]   # 完全反序
    y = np.r_[np.zeros(n0, dtype=int), np.ones(n1, dtype=int)]
    w = RNG.uniform(0.2, 5.0, size=n1 + n0)
    rw, _, _ = prod_ratio(WC_P6C, s, w, y)
    ru, _, _ = prod_ratio(WC_P6C, s, np.ones(n1 + n0), y)
    return rw, ru


def t_conc02():
    if WC_P6C is None:
        raise RuntimeError(f"生产函数装载失败: {LOAD_ERR}")
    rw, ru = _sep_case(True)
    ok = abs(rw - 1.0) <= 1e-12 and abs(ru - 1.0) <= 1e-12
    return ok, f"ratio(IOW加权,w~U[0.2,5])={rw:.15f} ratio(w=1)={ru:.15f} (n_pos=40,n_neg=60)"


def t_conc03():
    if WC_P6C is None:
        raise RuntimeError(f"生产函数装载失败: {LOAD_ERR}")
    rw, ru = _sep_case(False)
    ok = abs(rw) <= 1e-12 and abs(ru) <= 1e-12
    return ok, f"ratio(IOW加权,w~U[0.2,5])={rw:.15f} ratio(w=1)={ru:.15f} (n_pos=40,n_neg=60)"


# ================================================================ CONC-04
def t_conc04(fn, tag):
    """全并列 => 标准 AUROC 必须 0.5；生产严格 `>` 给 0。"""
    if fn is None:
        raise RuntimeError(f"生产函数装载失败: {LOAD_ERR}")
    n1, n0 = 60, 140
    s = np.full(n1 + n0, 2.5)
    y = np.r_[np.zeros(n0, dtype=int), np.ones(n1, dtype=int)]
    w = RNG.uniform(0.2, 5.0, size=n1 + n0)
    prod, num, den = prod_ratio(fn, s, w, y)
    ref, _ = ref_auroc_ties(s, y, w)
    ok = abs(prod - ref) <= 1e-12
    return ok, (f"{tag} production={prod:.6f} (num={num:.1f}/den={den:.1f}) "
                f"vs correct(R2,全并列=0.5)={ref:.6f}")


# ================================================================ CONC-05
def _load_sample(k=15):
    files = sorted(glob.glob(os.path.join(CACHE_DIR, "*.pkl")))
    if not files:
        raise FileNotFoundError(f"{CACHE_DIR} 下无 .pkl")
    step = max(1, len(files) // k)
    return files[::step][:k], len(files)


def t_conc05():
    """真实数据里的并列质量（诊断，不断言恒等式）。只读 p10_cache。"""
    sel, ntot = _load_sample(15)
    rows = []
    for f in sel:
        with open(f, "rb") as fh:
            d = pickle.load(fh)
        for key in ("s_e", "st_full"):
            v = np.asarray(d[key], dtype=float)
            u, c = np.unique(v, return_counts=True)
            rows.append((key, v.shape[0], u.shape[0], float(c[c > 1].sum()) / float(v.shape[0])))
        del d
        gc.collect()
    byk = {k: np.array([r[3] for r in rows if r[0] == k]) for k in ("s_e", "st_full")}
    n_rep = {k: np.array([r[1] for r in rows if r[0] == k]) for k in ("s_e", "st_full")}
    n_unq = {k: np.array([r[2] for r in rows if r[0] == k]) for k in ("s_e", "st_full")}
    allfrac = np.array([r[3] for r in rows])
    ok = float(allfrac.max()) < 0.01                  # 诊断：并列若无足轻重则 PASS
    observed = ("cache files %d/%d, tasks=%d | s_e n=%s uniq=%s tied-row-frac min/mean/max="
                "%.4f/%.4f/%.4f | st_full n=%s uniq=%s min/mean/max=%.4f/%.4f/%.4f"
                % (len(sel), ntot, len(sel),
                   f"{int(n_rep['s_e'].min())}-{int(n_rep['s_e'].max())}",
                   f"{int(n_unq['s_e'].min())}-{int(n_unq['s_e'].max())}",
                   byk["s_e"].min(), byk["s_e"].mean(), byk["s_e"].max(),
                   f"{int(n_rep['st_full'].min())}-{int(n_rep['st_full'].max())}",
                   f"{int(n_unq['st_full'].min())}-{int(n_unq['st_full'].max())}",
                   byk["st_full"].min(), byk["st_full"].mean(), byk["st_full"].max()))
    return ok, observed


_SCAN = {}


def _scan_real_data(k=8):
    """只读扫描 p10 缓存（s_e/ytes 与 st_full/m_t_full），量化并列对估计量的实际影响。"""
    if _SCAN:
        return _SCAN
    if WC_P6C is None:
        raise RuntimeError(f"生产函数装载失败: {LOAD_ERR}")
    sel, _ = _load_sample(k)
    src, tgt = [], []
    for f in sel:
        with open(f, "rb") as fh:
            d = pickle.load(fh)
        # --- 源侧：AUROC(s_e, yte)，无权（= IOW 在 w≡1 时的退化形式）---
        v = np.asarray(d["s_e"], dtype=float)
        y = np.asarray(d["yte"])
        ones = np.ones(v.shape[0])
        p, _, _ = prod_ratio(WC_P6C, v, ones, y)
        r, _ = ref_auroc_ties(v, y, ones)
        src.append((p, r))
        # --- 目标侧：OM 形式 Σ_{i≠j} m_i(1-m_j)1[s_i>s_j] / [(Σm)(Σ(1-m))] ---
        s = np.asarray(d["st_full"], dtype=float)
        m = np.asarray(d["m_t_full"], dtype=float)
        deno = float(m.sum() * (1.0 - m).sum())
        po = WC_P6C(s, m, 1.0 - m) / deno if deno > 0 else float("nan")
        # 定义式（并列记 0.5）：Σ_{i≠j} m_i(1-m_j)(1[s_i>s_j]+0.5·1[s_i=s_j])，按唯一分值分组精确求和；
        # 不依赖生产输出（否则生产修复后会重复计入并列修正）
        u, inv = np.unique(s, return_inverse=True)
        Mg = np.bincount(inv, weights=m)
        Ng = np.bincount(inv, weights=(1.0 - m))
        Qg = np.bincount(inv, weights=m * (1.0 - m))          # Σ_{i∈g} m_i(1-m_i)
        cumN = np.concatenate([[0.0], np.cumsum(Ng)[:-1]])    # 严格小于当前组的权重累积
        strict_pairs = float(np.sum(Mg * cumN))               # Σ_{i≠j: s_i>s_j} m_i(1-m_j)
        tie_pairs = float(np.sum(Mg * Ng - Qg))               # Σ_{i≠j: s_i=s_j} m_i(1-m_j)
        corr = (strict_pairs + 0.5 * tie_pairs) / deno if deno > 0 else float("nan")
        tgt.append((po, corr, tie_pairs / deno if deno > 0 else float("nan")))
        del d
        gc.collect()
    _SCAN["src"] = src
    _SCAN["tgt"] = tgt
    return _SCAN


def t_conc05b():
    """真实数据上缺陷 (b) 的幅度（源侧）：生产 ratio vs R2（s_e / yte）。"""
    sc = _scan_real_data(8)
    deltas = [r - p for p, r in sc["src"]]
    dmax, dmin = float(np.max(deltas)), float(np.min(deltas))
    ok = max(abs(dmin), abs(dmax)) <= 1e-9   # 并列 0.5 记分 => delta 必须双向为 0
    p0, r0 = sc["src"][0]
    return ok, ("tasks=%d delta=R2-production min/mean/max=%.4f/%.4f/%.4f ; e.g. "
                "production=%.4f vs correct(0.5-credit)=%.4f"
                % (len(deltas), dmin, float(np.mean(deltas)), dmax, p0, r0))


def t_conc05c():
    """目标侧诊断：OM 形式（st_full, m_t_full）上定义式（并列 0.5）与生产 ratio 的双向偏差。"""
    sc = _scan_real_data(8)
    dl = [c - p for p, c, _ in sc["tgt"]]
    tp = [x for _, _, x in sc["tgt"]]
    return (max(abs(min(dl)), abs(max(dl))) <= 1e-9,
            "tasks=%d OM: delta(def-0.5credit - production)= min/mean/max %.4f/%.4f/%.4f | "
            "tied ordered-pair mass/(Σm)(Σ(1-m)) min/mean/max %.4f/%.4f/%.4f"
            % (len(dl), float(np.min(dl)), float(np.mean(dl)), float(np.max(dl)),
               float(np.min(tp)), float(np.mean(tp)), float(np.max(tp))))


# ================================================================ CONC-06
def t_conc06():
    from sklearn.metrics import roc_auc_score
    rng = np.random.default_rng(777)
    n = 300
    s = rng.integers(0, 8, size=n).astype(float) * 0.5     # 只有 8 个取值 -> 大量并列
    y = (rng.random(n) < 0.4).astype(int)
    if y.sum() == 0 or y.sum() == n:
        raise RuntimeError("随机标签退化")
    skl = float(roc_auc_score(y, s))
    ref, _ = ref_auroc_ties(s, y, None)
    ok = abs(skl - ref) <= 1e-12
    return ok, f"sklearn={skl:.12f} reference_R2={ref:.12f} (n={n}, 取值数={len(np.unique(s))})"


def t_conc06b():
    """同一数据（与 CONC-06 同分布不同 seed）：sklearn（并列 0.5）必须等于生产 ratio；
    严格 > 生产漏记并列，diff 显著非零。"""
    from sklearn.metrics import roc_auc_score
    rng = np.random.default_rng(20261004)
    n = 300
    s = rng.integers(0, 8, size=n).astype(float) * 0.5
    y = (rng.random(n) < 0.4).astype(int)
    skl = float(roc_auc_score(y, s))
    prod, _, _ = prod_ratio(WC_P6C, s, np.ones(n), y)
    ok = abs(skl - prod) <= 1e-9
    return ok, f"sklearn={skl:.6f} production={prod:.6f} diff={skl - prod:+.6f} (n={n})"


def t_conc06c():
    from sklearn.metrics import roc_auc_score
    rng = np.random.default_rng(31415)
    n = 260
    s = rng.integers(0, 6, size=n).astype(float)
    y = (rng.random(n) < 0.45).astype(int)
    w = rng.uniform(0.2, 5.0, size=n)
    sklw = float(roc_auc_score(y, s, sample_weight=w))
    refw, _ = ref_auroc_ties(s, y, w)
    ok = abs(sklw - refw) <= 1e-9
    return ok, f"sklearn(sample_weight)={sklw:.12f} reference_R2(加权)={refw:.12f} (n={n})"


# ================================================================ 主流程
def main_body():
    L.banner("CONC — wconcordance / AUROC 并列语义")

    if LOAD_ERR:
        L.emit("CONC-00", "FAIL", "L.load_def(P6C/P10B,'wconcordance') 必须可装载",
               _clean(f"EXCEPTION {LOAD_ERR}"), "生产函数装载失败，后续测试无意义")
    else:
        L.emit("CONC-00", "PASS", "生产原文逐字抽取（load_def）成功，两副本非空",
               _clean("p6c=%s | p10b=%s" % (WC_P6C.__evs_origin__, WC_P10B.__evs_origin__)),
               "p6c:84-92 与 p10b:40-47 原文 exec（def 加尾随空行），非复制")

    L.guard("CONC-01", t_conc01,
            "AUROC 分子定义（Mann–Whitney 1947；Hanley&McNeil 1982）："
            "C(s,a,b)=Σ_{i≠j} a_i b_j (1[s_i>s_j]+0.5·1[s_i=s_j])：生产值必须等于 R1 显式"
            "二重循环（相对误差<=1e-9）；严格 >（冻结 docstring p6c:85 的写法）违反该定义",
            "功能层对拍：n=50/91/137/200，分值网格 5-13 档并叠加 35% 强制并列；"
            "R1 为逐对二重循环，与生产的 argsort+bincount 路径无关。"
            "本断言与 R1/expected 已按并列 0.5 语义同步（不再钉死冻结 docstring 的严格 >），"
            "生产按并列 0.5 修复后与 CONC-04/04B/05B/06B 一同转绿")

    L.guard("CONC-02", t_conc02,
            "AUROC 定义（Mann-Whitney 1947）：完全分离 => 1；管线形式 p6c:271-273 "
            "den=(Σw·1[y=1])(Σw·1[y=0])",
            "n_pos=40 n_neg=60，正类 score 全部严格高于负类（无跨类并列）")

    L.guard("CONC-03", t_conc03,
            "AUROC 定义（Mann-Whitney 1947）：完全反序 => 0",
            "n_pos=40 n_neg=60，正类 score 全部严格低于负类")

    L.guard("CONC-04", lambda: t_conc04(WC_P6C, "p6c副本"),
            "AUROC 并列语义（Hanley&McNeil 1982；sklearn 文档）：全并列 => 0.5。"
            "生产严格 > 使 num=0 -> 违反恒等式",
            "缺陷 (b) 的最小反例：n_pos=60 n_neg=140，所有 score 并列")

    L.guard("CONC-04B", lambda: t_conc04(WC_P10B, "p10b副本"),
            "同 CONC-04（p10b:40-47 是 p6c:84-92 的同源副本，应同罪）",
            "两副本必须同罪；p10b_subsample_curves.py 也用同一严格 > 实现")

    _guard_warn("CONC-05", t_conc05,
                "诊断项（data characterisation，非恒等式，WARN）：frozen p10 缓存 s_e / st_full "
                "的精确并列发生率；并列行占比 <1% 视为可忽略",
                "只读扫描 p10_cache（202 个任务中按文件顺序等距抽 15 个）；"
                "key 定义见 p10a_prepare_cache.py:307-309")

    if WC_P6C is not None:
        L.guard("CONC-05B", t_conc05b,
                "AUROC 恒等式（Mann–Whitney 1947；Hanley&McNeil 1982）：生产 ratio 必须等于"
                "并列记 0.5 的 AUROC（R2），双向 |delta|<=1e-9；严格 > 缺陷使 delta 恒正，"
                "过度记分（如并列记 1.0）变体使 delta 恒负",
                "真实 source internal test（s_e, yte，无权 = 加权 IOW 在 w≡1 时的退化形式）；"
                "严格 > 生产上实测低估 0.4-1.3 个百分点。加权 AUROC_IOW 的 w_te 不在 "
                "p10_cache 键中，该幅度不外推到加权形式；只读缓存")

        _guard_warn("CONC-05C", t_conc05c,
                    "诊断项（data characterisation，非恒等式）：目标侧 OM 形式（st_full, m_t_full）"
                    "上定义式（并列 0.5）与生产 ratio 的双向偏差（|delta|<=1e-9 方为一致）",
                    "OM 分子用同一 wconcordance（p6c:274-275），故并列漏记同样命中 OM/DR 的 target 项")

        L.guard("CONC-06", t_conc06,
                "sklearn.metrics.roc_auc_score（并列记 0.5）== R2 带并列 AUROC（sklearn 文档）",
                "外部独立对照：sklearn 是第三方实现，与 R2 无关")

        L.guard("CONC-06B", t_conc06b,
                "AUROC 并列语义（Hanley&McNeil 1982；sklearn 文档）：同一数据上 "
                "sklearn.roc_auc_score 必须等于生产 ratio（|diff|<=1e-9，并列统一记 0.5）",
                "缺陷 (b) 的外部对照：未修复生产（严格 >，并列记 0）diff 为 +0.06 量级；"
                "并列 0.5 修复后 diff=0 转绿")

        L.guard("CONC-06C", t_conc06c,
                "sklearn 的 sample_weight 加权 AUROC == R2 加权形式（sklearn 文档）",
                "加权形式的交叉验证（管线 IOW 用的是加权 ratio）")

    npass = sum(1 for _, s in L._SEEN if s == "PASS")
    nfail = sum(1 for _, s in L._SEEN if s == "FAIL")
    nwarn = sum(1 for _, s in L._SEEN if s == "WARN")
    print(f"SUMMARY|CONC|tests={len(L._SEEN)} PASS={npass} FAIL={nfail} WARN={nwarn}", flush=True)


if __name__ == "__main__":
    try:
        main_body()
    except Exception:
        traceback.print_exc()
    sys.exit(0)
