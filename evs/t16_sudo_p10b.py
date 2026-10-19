# -*- coding: utf-8 -*-
r"""t16_sudo_p10b.py — family SUDO-P10B：p10b 那份 SUDO 副本的 **exact-identity gate**

================================ 为什么单独一条 ================================
manifest 的 `SUDO_P10B` 在 2026-10-05 是 **UNKNOWN / 0 条 gate**：
唯一相关的 `FORK-03` 是**对原论文双分支 SUDO 的符合性测试**，而 CS-04 已裁决**保留**单分支实现
⇒ FORK-03 永不可能 PASS，故它不是 gate。p6c 侧尚有 `SUDO-04X`，p10b 侧没有对应物。

================================ GPT6 的约束（本模块的设计依据） ================================
    "独立 reference 必须直接从 fixture 原始输入按冻结的 single-branch + 自定义 AURCC 公式计算，
     不能调用生产 helper，也不能拿生产中间量再算面积。"

因此本模块**不写"把生产块抄一遍"的 reference** —— 那种做法只是同源自证。改为三条**不同性质**的判据：

  SUDO-P10B-01  块在受控 fixture 上**逐字执行**（BLOCK_P10B_SUDO）→ 值有限、nbin 一致
  SUDO-P10B-02  **独立实现 §5.3.1 冻结的聚合**（|bin 均值−0.5| 降序 → completeness/reliability
                → trapezoid），从块**自己**的 per-bin 归因 (sbins, frac, conf) 重算 →
                与块返回的 SUDO_AURCC **逐位相等**。这一条验的是"面积怎么算"。
  SUDO-P10B-03  **用 instrumentation 验冻结的单分支规格**（这一条才真正独立于生产代码）：
                  · 辅助分类器的 fit 次数必须 == SUDO_REP × n_valid_bins
                    （**单分支**每 bin 每 repeat 恰好 1 次；原论文双分支会是 2 倍）
                  · 每次抽样的 target 索引必须落在 **reference 独立算出的** bin 内
                  · 每次抽样的 source 索引必须落在 **reference 独立算出的** `1−yb` 类池内
                  reference 自己从 fixture 原始输入算 qs / bid / yb / pool，**不读生产中间量**。
  SUDO-P10B-04  **非空转对照**：把冻结聚合的排序方向篡改（conf 升序代替降序）后，
                SUDO-P10B-02 的判据必须判红。没有这一条，"逐位相等"只是声明。

★ 与本模块的关系：`t04` 的 SUDO-01/03/05 是**对原论文的偏离披露**（CS-04 要求写明两处偏离），
  它们**故意**在单分支实现上判红；本模块才是 p10b 副本的**正确性 gate**。

输出契约同 evs_load：`RESULT|tid|PASS|FAIL|WARN|expected|observed|note`
"""
import os, sys, io
import numpy as np

import evs_load as L

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score

P10B = L.P10B
BLOCK = L.BLOCK_P10B_SUDO                      # (P10B, 116, 148)

# 冻结常量：从**生产文件**逐字读取（不硬编码，避免与文件漂移）
SEED = None
NBIN = SUDO_REP = SUDO_CAP = SUDO_FRAC = None
_src = io.open(P10B, encoding="utf-8", errors="replace").read().split("\n")
for _ln in _src:
    _s = _ln.strip()
    if _s.startswith("SEED ="):
        SEED = int(_s.split("=", 1)[1].strip())
    if _s.startswith("MMD_M, MMD_K, NBIN"):
        NBIN = int(_s.split("=", 1)[1].split(",")[2].strip())
    if _s.startswith("SUDO_REP, SUDO_CAP, SUDO_FRAC"):
        _v = [x.strip() for x in _s.split("=", 1)[1].split(",")]
        SUDO_REP, SUDO_CAP, SUDO_FRAC = int(_v[0]), int(_v[1]), float(_v[2])

print("=" * 100, flush=True)
print("t16 — family SUDO-P10B：p10b SUDO 副本的 exact-identity gate", flush=True)
print("=" * 100, flush=True)
print(f"  生产树 = {L.PROD_DIR}   模式 = {'PATCHED' if L.IS_PATCHED else 'REAL'}", flush=True)
print(f"  块 = {L.block_source(BLOCK)}   冻结常量 SEED={SEED} NBIN={NBIN} "
      f"SUDO_REP={SUDO_REP} SUDO_CAP={SUDO_CAP} SUDO_FRAC={SUDO_FRAC}", flush=True)
print("", flush=True)


# ------------------------------------------------------------------ 受控 fixture
def make_fixture(seed=424242, n_tr=600, n_va=300, n_t=900, n_feat=4):
    """合成一个**两个类都有、bin 内样本充足、pool 足够**的 fixture。

    刻意让 target 分数与 source 类分布有可分辨的差异，使每个 bin 都获得有效 AUC。
    """
    rng = np.random.default_rng(seed)
    ytr = (rng.random(n_tr) < 0.5).astype(int)
    yva = (rng.random(n_va) < 0.5).astype(int)
    mu = lambda y: (np.where(y == 1, 0.9, -0.9))[:, None]          # 类间可分
    Dtr = mu(ytr) + rng.normal(0, 1.0, (n_tr, n_feat))
    Dva = mu(yva) + rng.normal(0, 1.0, (n_va, n_feat))
    z = rng.normal(0, 1.0, n_t)
    pt_full = 1.0 / (1.0 + np.exp(-(2.2 * z + rng.normal(0, 0.6, n_t))))   # 0..1 分数
    Dt = rng.normal(0, 1.6, (n_t, n_feat))
    A = dict(ytr=ytr, Dtr_std=Dtr, yva=yva, Dva_std=Dva, Dt_std=Dt, pt_full=pt_full)
    return A


A_FIX = make_fixture()
PT = A_FIX["pt_full"]
DT = A_FIX["Dt_std"]
N_T = len(PT)


# ------------------------------------------------------------------ 独立 reference：只从 fixture 原始输入算
def ref_bins(pt_s, nbin=NBIN):
    """**独立**重算分箱与每 bin 的成员 / 伪标签 / 类池。不读生产中间量。

    用 searchsorted 走一条与生产 `np.digitize` **不同**的代码路径（同样的等频分位数定义）。
    """
    qs = np.quantile(pt_s, np.linspace(0, 1, nbin + 1))
    inner = qs[1:-1]
    bid = np.searchsorted(inner, pt_s, side="right")
    out = []
    for b in range(nbin):
        ib = np.where(bid == b)[0]
        if len(ib) < 10:
            continue
        yb = 1 if float(pt_s[ib].mean()) >= 0.5 else 0
        out.append(dict(b=b, idx=ib, yb=yb, mean=float(pt_s[ib].mean())))
    return out


def ref_aurcc_from_parts(sbins, frac, conf):
    """**独立实现** §5.3.1 冻结的聚合（与生产代码不同的写法）：
       |bin 均值−0.5| 降序 → completeness = cumsum(frac)/Σfrac
                            → reliability  = cumsum(|S|·frac)/cumsum(frac)
                            → 梯形积分。
    这里用显式循环写梯形（生产用 np.trapezoid），以证明"逐位相等"不是同源自证。
    """
    order = sorted(range(len(conf)), key=lambda i: -conf[i])       # 降序，显式
    f = [frac[i] for i in order]
    r = [abs(sbins[i]) for i in order]
    tot = sum(f)
    comp, rel, cf, cr = [], [], 0.0, 0.0
    for fi, ri in zip(f, r):
        cf += fi
        comp.append(cf / tot)
        cr += ri * fi
        rel.append(cr / cf)
    area = 0.0
    for i in range(1, len(comp)):
        area += 0.5 * (rel[i] + rel[i - 1]) * (comp[i] - comp[i - 1])
    return float(area)


# ------------------------------------------------------------------ 带 instrumentation 的块执行
class _CountingLR(LogisticRegression):
    """计数的辅助分类器（t04 同款技术）：只记录 fit 次数与训练集大小，不改变计算。"""
    N_FIT = 0

    def fit(self, X, y, **kw):
        _CountingLR.N_FIT += 1
        return super().fit(X, y, **kw)


class _RecordingRng:
    """记录每次 choice 的 (population, size) —— 只记录，不改变抽签结果。"""
    def __init__(self, rng):
        self._r = rng
        self.calls = []

    def choice(self, a, size=None, replace=True, p=None):
        out = self._r.choice(a, size=size, replace=replace, p=p)
        self.calls.append((np.asarray(a).copy() if not np.isscalar(a) else a, size))
        return out

    def __getattr__(self, k):
        return getattr(self._r, k)


def run_block(rng_seed=777, block=BLOCK, extra_env=None):
    """逐字执行 p10b 的 SUDO 块。返回 (rec, sbins, frac, conf, recorder, n_fit)。"""
    _CountingLR.N_FIT = 0
    rec_rng = _RecordingRng(np.random.default_rng(rng_seed))
    env = dict(np=np, LogisticRegression=_CountingLR, roc_auc_score=roc_auc_score,
               pt_s=PT, Dt_std_s=DT, n=N_T, A=A_FIX, rng=rec_rng,
               NBIN=NBIN, SUDO_REP=SUDO_REP, SUDO_CAP=SUDO_CAP, SUDO_FRAC=SUDO_FRAC,
               SEED=SEED, rec={})
    if extra_env:
        env.update(extra_env)
    out, _src_txt = L.exec_block(block, env, names=["rec", "sbins", "frac", "conf"])
    return (out["rec"], out["sbins"] or [], out["frac"] or [], out["conf"] or [],
            rec_rng, _CountingLR.N_FIT)


REC, SB, FR, CF, RECRNG, NFIT = run_block()


# ------------------------------------------------------------------ SUDO-P10B-01 逐字执行
def t01():
    if REC is None or "SUDO_AURCC" not in REC:
        return False, "块未产出 rec['SUDO_AURCC']"
    v = REC["SUDO_AURCC"]; nb = REC.get("SUDO_nbin")
    ok = np.isfinite(v) and 0.0 <= v <= 1.0 and nb == len(SB) and nb >= 1
    return ok, (f"{L.block_source(BLOCK)} 逐字执行于受控 fixture（pt_s {N_T} 行、两 source 类、"
                f"4 特征）：SUDO_AURCC={v:.12g}，SUDO_nbin={nb}，per-bin 归因 len(sbins)={len(SB)}；"
                f"值域 [0,1] 与 nbin 一致性 {'成立' if ok else '**不成立**'}")


L.guard("SUDO-P10B-01", t01,
        "p10b 的 SUDO 副本必须在受控 fixture 上**逐字执行**并产出有限、在 [0,1] 内的 AURCC，"
        "且 rec['SUDO_nbin'] 与块自身的 per-bin 归因长度一致（副本可运行且自洽）",
        "块 = BLOCK_P10B_SUDO = p10b_subsample_curves.py:116-148（evs_load 逐字抽取）")


# ------------------------------------------------------------------ SUDO-P10B-02 聚合 identity
def t02():
    if not SB:
        return False, "无 per-bin 归因（fixture 未产生有效 bin）"
    ref = ref_aurcc_from_parts(SB, FR, CF)
    got = float(REC["SUDO_AURCC"])
    d = abs(ref - got)
    return d <= 1e-12, (
        f"独立实现 §5.3.1 聚合（显式降序排序 + 显式梯形循环；生产用 np.argsort/np.trapezoid）"
        f"从块自己的 per-bin 归因重算：ref={ref:.15g} vs block={got:.15g}，absdiff=={d:.3e}"
        f"；per-bin: n={len(SB)}，Σfrac={float(np.sum(FR)):.6f}，"
        f"conf 降序前 3 = {sorted(CF, reverse=True)[:3]}")


L.guard("SUDO-P10B-02", t02,
        "冻结聚合 identity（§5.3.1）：AURCC = 按 |bin 均值−0.5| 降序排列后，"
        "completeness=cumsum(frac)/Σfrac、reliability=cumsum(|S|·frac)/cumsum(frac) 的梯形积分；"
        "由**独立实现**从块的 per-bin 归因重算，必须与块返回的 SUDO_AURCC 逐位相等（≤1e-12）",
        "「独立」= 不同的排序写法与不同的梯形写法；这一条验的是**面积怎么算**")


# ------------------------------------------------------------------ SUDO-P10B-03 冻结单分支规格
def t03():
    rb = ref_bins(PT)
    valid = [d for d in rb if len(np.where(A_FIX["ytr"] == 1 - d["yb"])[0]) >= 10]
    expect_fits = SUDO_REP * len(valid)
    src0 = np.where(A_FIX["ytr"] == 0)[0]
    src1 = np.where(A_FIX["ytr"] == 1)[0]

    problems = []
    if NFIT != expect_fits:
        problems.append(f"辅助分类器 fit 次数 {NFIT} != SUDO_REP×n_valid_bins = {expect_fits}"
                        f"（**单分支**每 bin 每 repeat 恰好 1 次；双分支会是 2 倍）")

    # 逐次抽样核对：target 必须落在 reference 独立算出的 bin 内；source 必须在 1-yb 类池内
    bid = np.searchsorted(np.quantile(PT, np.linspace(0, 1, NBIN + 1))[1:-1], PT, side="right")
    # ★ 分类依据：**抽样次序**，不是索引范围。target 抽样的总体是 `ib_`（PT 的行号，≤899），
    #   source 抽样的总体是 `pool`（Dtr 的行号，≤599）——**两者范围重叠**，用"上界"区分会误判
    #   （我第一版就是这么错的：把 source 抽样当成 target，`bid[pop]` 自然散落 10 个 bin）。
    #   块内每轮 rr 固定两次抽样，先 target 后 source，故按 calls 的**奇偶**判定，并断言偶数次。
    calls = [c for c in RECRNG.calls if c[1] is not None]
    if len(calls) % 2 != 0:
        problems.append(f"抽样次数 {len(calls)} 为奇数（块内每轮应为 target+source 两次）")
    n_tgt = n_src = 0
    for i, (pop, size) in enumerate(calls):
        pop = np.asarray(pop)
        if pop.ndim == 0:
            continue
        if i % 2 == 0:                      # target 抽样
            bs = set(bid[pop].tolist())
            if len(bs) != 1:
                problems.append(f"一次 target 抽样跨越多个 bin：{sorted(bs)}")
            else:
                b = bs.pop()
                if not any(x["b"] == b for x in rb):
                    problems.append(f"target 抽样落在 reference 判定为**无效**的 bin {b}")
            n_tgt += 1
        else:                               # source 抽样
            hit0 = np.isin(pop, src0).all(); hit1 = np.isin(pop, src1).all()
            if not (hit0 or hit1):
                problems.append("一次 source 抽样**跨越两个真实类别**（伪标签池应按类分池）")
            n_src += 1
    n_chk = n_tgt + n_src
    return (not problems), (
        f"reference 独立算出的有效 bin = {len(valid)}（共 {len(rb)} 个 ≥10 样本的 bin）；"
        f"期望 fit 次数 = {SUDO_REP}×{len(valid)} = {expect_fits}，实测 {NFIT}"
        + ("" if not problems else f"；违规 {len(problems)} :: " + " | ".join(problems[:3]))
        + f"；逐次抽样核对 {n_chk} 次（target {n_tgt} / source {n_src}，按**次序**奇偶判定）")


L.guard("SUDO-P10B-03", t03,
        "冻结单分支规格的**独立**符合性（instrumentation）："
        "(a) 辅助分类器 fit 次数必须恰为 SUDO_REP × n_valid_bins（单分支=1/bin/repeat）；"
        "(b) 每次 target 抽样必须整体落在同一个 bin 内；"
        "(c) 每次 source 抽样必须整体落在**同一真实类别**的池内（伪标签池按 1−yb 分池）。"
        "bins / 类别池由 reference 自己从 fixture 原始输入算出，**不读生产中间量**",
        "★ 这一条才真正独立：即使有人把生产块的赋值重排，count 与池成员关系仍会被独立算出并比对")


# ------------------------------------------------------------------ SUDO-P10B-04 非空转对照
def _mutate_aggregation():
    """把块源码里的降序排序改为**升序**（已知错误方式），返回变异后的 block 文本。"""
    src = L.grab(*BLOCK)
    old = "order = np.argsort(-np.array(conf))"
    assert src.count(old) == 1, "找不到聚合排序行"
    return src.replace(old, "order = np.argsort(np.array(conf))", 1)


def t04():
    mut = _mutate_aggregation()
    env = dict(np=np, LogisticRegression=LogisticRegression, roc_auc_score=roc_auc_score,
               pt_s=PT, Dt_std_s=DT, n=N_T, A=A_FIX, rng=np.random.default_rng(777),
               NBIN=NBIN, SUDO_REP=SUDO_REP, SUDO_CAP=SUDO_CAP, SUDO_FRAC=SUDO_FRAC,
               SEED=SEED, rec={})
    ns = dict(env)
    exec(compile(mut, "<mutant:p10b-sudo-aggregation-order>", "exec"), ns)
    got = float(ns["rec"]["SUDO_AURCC"])
    ref = ref_aurcc_from_parts(ns["sbins"], ns["frac"], ns["conf"])
    d = abs(ref - got)
    # 对照：同一 reference 在**未篡改**块上必须逐位相等（已在 02 证明），在篡改块上必须不等
    ok = d > 1e-9
    return ok, (f"把冻结聚合的排序方向篡改（conf **升序**代替降序）后：块返回 {got:.15g}，"
                f"独立 reference（仍按冻结的降序定义）算出 {ref:.15g}，absdiff=={d:.3e} —— "
                f"{'被检出（判据非空转）' if ok else '**未被检出，判据空转**'}；"
                f"对照：同一 reference 在未篡改块上 absdiff=≤1e-12（见 SUDO-P10B-02）")


L.guard("SUDO-P10B-04", t04,
        "非空转对照：把 §5.3.1 冻结聚合的排序方向篡改（|bin 均值−0.5| 升序代替降序）后，"
        "SUDO-P10B-02 的独立 reference 必须与块返回的 AURCC **不相等**（absdiff=>1e-9）",
        "★ 没有这一条，「逐位相等」只是一句声明；有这一条，才证明 reference 真的在算面积")

print("\n[t16_sudo_p10b done]", flush=True)
