# -*- coding: utf-8 -*-
r"""evs_selftest — 验证**装载器本身**能把生产代码逐字取出并跑通

================================ 为什么需要它 ================================
"审计器也要被审计"。如果 evs_load 取错了行、或取出的块跑不起来，
后面六个测试族会集体给出假的 PASS/FAIL。本文件先证明：

  1. 顶层 def 抽取正确（与源文件对应行逐字一致）
  2. 内联块抽取正确（行号区间与原文件一致）
  3. 用合成数据能把 L3-B / SUDO 两个真实块跑通并产出有限值

用法：python evs_selftest.py
"""
import sys, os, hashlib
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import evs_load as L
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler

L.banner("装载器自测（infrastructure）")

# ---------- 1. 顶层 def 抽取 ----------
def t1():
    fn = L.load_def(L.P6B, "mmd_rbf")
    src = fn.__evs_src__
    ok_head = src.startswith("def mmd_rbf(")
    ok_body = "unbiased U-statistic" in src
    return (ok_head and ok_body), f"origin={fn.__evs_origin__}, len={len(src)}, head_ok={ok_head}, doc_ok={ok_body}"
L.guard("LOAD-01", t1, "P6B 顶层 def mmd_rbf 被逐字取出（首行 def、含 docstring）")

def t2():
    fn = L.load_def(L.P6C, "wconcordance")
    ok = fn.__evs_src__.startswith("def wconcordance(") and "严格小于" in fn.__evs_src__
    return ok, f"origin={fn.__evs_origin__}, len={len(fn.__evs_src__)}"
L.guard("LOAD-02", t2, "P6C 顶层 def wconcordance 被逐字取出")

# ---------- 2. 内联块抽取：与源文件同区间逐字一致 ----------
def t3():
    for blk, name in [(L.BLOCK_L3B, "L3-B"), (L.BLOCK_SUDO, "SUDO"), (L.BLOCK_STORE, "STORE")]:
        path, a, b = blk
        raw = "\n".join(L.read_lines(path)[a - 1:b])
        ded = L.grab(path, a, b)
        # dedent 只应去掉**公共前导缩进**，不得改动任何非空白内容。
        # 判据：两侧去空白后逐行相同，且去掉的缩进量在所有非空行上一致。
        ra, da = raw.split("\n"), ded.split("\n")
        if len(ra) != len(da):
            return False, f"{name} 区间 {a}-{b} 行数改变 {len(ra)}->{len(da)}"
        for i, (x, y) in enumerate(zip(ra, da)):
            if x.strip() != y.strip():
                return False, f"{name} 第 {a+i} 行内容被改动：{x.strip()[:60]!r} -> {y.strip()[:60]!r}"
        ind = {len(x) - len(x.lstrip()) for x in da if x.strip()}
        if ind and min(ind) != 0:
            return False, f"{name} 块去缩进后仍有公共前导空格 min={min(ind)}（exec 会失败）"
    return True, "L3-B 257-289 / SUDO 291-330 / STORE 121-144：去空白后逐行一致，且去缩进后无公共缩进"
L.guard("LOAD-03", t3, "三个内联块的行号区间与源文件逐行一致")

# ---------- 3. 用合成数据跑通真实块 ----------
rng = np.random.default_rng(0)
NF = 5
Dtr = rng.normal(size=(400, NF)); Dva = rng.normal(size=(100, NF)); Dte = rng.normal(size=(100, NF))
Dt = rng.normal(size=(300, NF))
w = rng.normal(size=NF); b0 = -0.3
ytr = (rng.random(400) < 1 / (1 + np.exp(-(Dtr @ w + b0)))).astype(int)
yva = (rng.random(100) < 0.5).astype(int)
yte = (rng.random(100) < 0.5).astype(int)
sc = StandardScaler().fit(Dtr)
mdl = LogisticRegression(max_iter=2000, random_state=1).fit(sc.transform(Dtr), ytr)
Ds_std = sc.transform(np.vstack([Dtr, Dva, Dte]))
Dt_std = sc.transform(Dt)
i_tr = np.arange(len(ytr)); i_va = np.arange(len(yva)) + len(i_tr); i_te = np.arange(len(yte)) + len(i_tr) + len(i_va)
M = dict(sc=sc, mdl=mdl, Dtr=Dtr, Dva=Dva, Dte=Dte, ytr=ytr, yva=yva, yte=yte)
wconcordance = L.load_def(L.P6C, "wconcordance")
SEED = 20261005
NBIN, SUDO_REP, SUDO_CAP, SUDO_FRAC = 10, 5, 200, 0.5
pt = mdl.predict_proba(Dt_std)[:, 1]
sc_t = mdl.decision_function(Dt_std)
rec = {}
env = dict(np=np, LogisticRegression=LogisticRegression, roc_auc_score=__import__(
    "sklearn.metrics", fromlist=["roc_auc_score"]).roc_auc_score,
    M=M, Ds_std=Ds_std, Dt=Dt, Dt_std=Dt_std, i_te=i_te, SEED=SEED, wconcordance=wconcordance,
    pt=pt, sc_t=sc_t, rec=rec, NBIN=NBIN, SUDO_REP=SUDO_REP, SUDO_CAP=SUDO_CAP,
    SUDO_FRAC=SUDO_FRAC, rng=np.random.default_rng(SEED + 11),
    Ds_tr_std=Ds_std[:len(ytr)])

def t4():
    out, _ = L.exec_block(L.BLOCK_L3B, env, names=["rec"])
    vals = {k: out["rec"].get(k) for k in ("AUROC_IOW", "AUROC_OM", "AUROC_DR")}
    ok = all(v is not None and np.isfinite(v) for v in vals.values())
    return ok, ", ".join(f"{k}={v:.4f}" for k, v in vals.items())
L.guard("LOAD-04", t4, "BLOCK_L3B（p6c 257-289）在合成数据上跑通并产出 3 个有限 AUROC")

def t5():
    env2 = dict(env); env2["rec"] = {}
    out, _ = L.exec_block(L.BLOCK_SUDO, env2, names=["rec"])
    v = out["rec"].get("SUDO_AURCC")
    return (v is not None and np.isfinite(v)), f"SUDO_AURCC={v}"
L.guard("LOAD-05", t5, "BLOCK_SUDO（p6c 291-330）在合成数据上跑通并产出有限 AURCC")

def t6():
    # 装载器必须能识别出块所需变量，且不静默吞掉 NameError
    try:
        L.exec_block(L.BLOCK_L3B, {"np": np}, names=["rec"])
        return False, "缺变量时**没有**抛异常（装载器会静默给假结果）"
    except NameError as e:
        return True, f"缺变量正确抛出 NameError: {e}"
    except Exception as e:
        return True, f"缺变量抛出 {type(e).__name__}: {e}"
L.guard("LOAD-06", t6, "装载器在缺变量时抛异常，而不是静默产出假结果")


# ---------- 3b. 块**首行指纹**：抓"行号错位指到别的代码" ----------
def t7():
    """对 evs_load 里**每一个** BLOCK_*，断言按当前行号抽出的块，其首行指纹 == 登记值。

    ★ 为什么必须有这条（2026-10-06 事故）：
      t14 的 MC-02 只验 **manifest 的 anchors** 能按文本解析，**不验** evs 代码里
      exec_block 的行号元组 —— 块被指到同一文件里的**另一段代码**时，MC-02 照样 PASS
      （只要那段也能 exec）。
      实测：给 p10a/p10b 插入代码后，`BLOCK_P10B_L3`(150,174) 变成了抽出 SUDO 块、
      `BLOCK_P10A_MMD`(273,280) 变成了抽出"[校验闸门]通过"，
      **没有任何判据报警**，是人工按内容对拍才发现的。
    ★ 覆盖面声明（写在此以免被当成更强的保证）：本判据看的是**块的首行**。
      若两个块的首行恰好相同，它会漏；若有人**故意**挪走行号并同步更新指纹，它也拦不住。
      它拦的是"改了行数却忘了同步块区间"这一**实际发生过**的失误，
      并强制把改块变成一次**有意识的行为**。
    ★ 非空转：指纹表为空、或某个 BLOCK_* 未登记指纹，都判红。"""
    fps = dict(getattr(L, "BLOCK_FINGERPRINTS", {}) or {})
    # ★ 只取**元组**形态的 BLOCK_*（块是 (path, a, b)）。
    #   不能只用 startswith("BLOCK_") —— `BLOCK_FINGERPRINTS` 本身也以 BLOCK_ 开头，
    #   第一版就是这么写的，于是判据把自己当成"未登记的块"报了红。
    names = sorted(n for n in dir(L)
                   if n.startswith("BLOCK_") and isinstance(getattr(L, n), tuple))
    if not fps:
        return False, "BLOCK_FINGERPRINTS 为空 —— 判据无靶子"
    missing = [n for n in names if n not in fps]
    if missing:
        return False, f"这些 BLOCK_* 没有登记指纹：{missing}"
    bad, okn = [], 0
    for n in names:
        blk = getattr(L, n)
        try:
            first = L.grab(*blk).split("\n")[0]
        except Exception as e:                                   # noqa: BLE001
            bad.append(f"{n}: 抽取失败 {type(e).__name__}: {e}")
            continue
        h = hashlib.sha256(first.encode("utf-8")).hexdigest()
        if h == fps[n]:
            okn += 1
        else:
            bad.append(f"{n}({L.block_source(blk)}): 首行={first.strip()[:44]!r} "
                       f"指纹 {h[:12]}… != 登记 {fps[n][:12]}…")
    return (not bad), (f"{okn}/{len(names)} 个块的首行指纹与登记一致"
                       + (f"；**错位**：{bad}" if bad else "；无错位"))
L.guard("LOAD-07", t7,
        "块首行指纹完整性：evs_load 里每个 BLOCK_* 按**当前行号**抽出的首行，"
        "其 sha256 必须等于 BLOCK_FINGERPRINTS 的登记值。"
        "这条堵的是「改了行数/行号却忘记同步块区间」——MC-02 看不见这种错位",
        "★ 2026-10-06 实测事故：给 p10a/p10b 插代码后 BLOCK_P10B_L3 变成了抽出 SUDO 块，无判据报警")

print("\n[selftest done]")
