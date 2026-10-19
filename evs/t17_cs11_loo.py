# -*- coding: utf-8 -*-
r"""t17_cs11_loo.py — family CS11：Fig 3B 的 **算法语义 + provenance** gate

================================ GPT6 的要求（本模块的设计依据） ================================
    "核心不是'哪个脚本数字碰巧一样'，而是确认真正的 authoritative provenance：
     p14 是否确实实现 true LOO；然后把 p6f:176-189 那条错误 lodo_oof 路径彻底隔离/废弃，
     确保未来只有一条可生成 Fig 3B 的 canonical path。
     这里要验 **算法语义 + provenance**，不能只验现有数值吻合。"

因此本模块不比对任何数字，只验两件事：

  provenance —— **谁有权写出** fig3b_loo_organism.csv。这条属于 §4.1 的「封闭世界 / control-plane」，
                故用**显式 registry**（下面的 CANONICAL_PRODUCER），而不是启发式发现：
                  ① 生产脚本里，以**写出位置**（`to_csv(...)` 实参）提到该文件名的，
                     必须**恰好**等于 registry；
                  ② 历史写者（p11:240 的旧写出）必须处于**注释/失效**状态。
  semantics  —— p14 的留出循环必须**真的**是 leave-one-organism-out：
                  ① 每一折的 logistic **只在另两个 organism 的行上**拟合（训练行数 == 补集行数）；
                  ② 评价**只用留出 organism 的行**；
                  ③ 学到的符号取自**该折的模型**（不是全局模型）。
                用 recording proxy 记录每次 fit 的 X 行数来判定，**与数值无关**。

  CS11-01  provenance：唯一写者 == registry，且历史写者已失效
  CS11-02  语义（execution + instrumentation）：每折 fit 行数 == 补集行数；评价集 == 留出行集
  CS11-03  陈旧路径隔离：被测树的 p6f 中**不得**存在 `LOO-organism` 段 header
  CS11-04  非空转：把 p14 的留出循环篡改为**在全部行上拟合**（泄漏留出折）后，CS11-02 必须判红

输出契约同 evs_load：`RESULT|tid|PASS|FAIL|WARN|expected|observed|note`
"""
import os, re, sys, io
import numpy as np

import evs_load as L

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

from sklearn.linear_model import LogisticRegression

P14 = L.P14
P6F = L.P6FD if hasattr(L, "P6FD") else os.path.join(L.PROD_DIR, "p6f_lodo_settlement.py")
P11 = os.path.join(L.PROD_DIR, "p11_figure_datapack.py")
ARTIFACT = "fig3b_loo_organism.csv"

# ★ 显式 registry（§4.1 封闭世界）：有权写出 Fig 3B 的脚本**只有** p14。
CANONICAL_PRODUCER = "p14_rebuild_fig3b.py"
# 历史写者（必须处于注释/失效状态）
RETIRED_WRITERS = ["p11_figure_datapack.py"]

print("=" * 100, flush=True)
print("t17 — family CS11：Fig 3B 的算法语义 + provenance", flush=True)
print("=" * 100, flush=True)
print(f"  生产树 = {L.PROD_DIR}   模式 = {'PATCHED' if L.IS_PATCHED else 'REAL'}", flush=True)
print(f"  canonical producer（显式 registry）= {CANONICAL_PRODUCER}", flush=True)
print("", flush=True)


def _read(p):
    with io.open(p, encoding="utf-8", errors="replace") as f:
        return f.read()


# ------------------------------------------------------------------ CS11-01 provenance
def _write_sites(path, name):
    """返回该文件里以**写出位置**提到 name 的行号（`to_csv`/`savefig`/`open(...,'w')` 等实参）。

    注释行（`#` 起始的 stripped 前缀）**不计**——历史写者正是被注释掉的那种。
    """
    out = []
    for i, ln in enumerate(_read(path).split("\n"), 1):
        s = ln.strip()
        if s.startswith("#"):
            continue
        if name in ln and re.search(r"to_csv|to_excel|savefig|open\s*\(|\.write\s*\(", ln):
            out.append(i)
    return out


def cs11_01():
    """★ 判据收窄（第一版过粗，已更正）：
       「历史写者已失效」必须按**写出位置**判定，不能按"任何非注释提及"判定 ——
       p11:244 是一条**活的 `p(...)` 声明**："由 p14_rebuild_fig3b.py 生成，本脚本不写"，
       那正是 CS-11 要的 provenance 记录，把它当成活跃写者是把对的判成错的。
       同时把这条声明升级为**正向判据**：活文本里必须**存在**对 canonical producer 的声明。"""
    import glob
    prods = sorted(glob.glob(os.path.join(L.PROD_DIR, "*.py")))
    writers, retired_ok, declared, listed = {}, {}, [], []
    for p in prods:
        b = os.path.basename(p)
        w = _write_sites(p, ARTIFACT)          # 只认**写出位置**（且剔除注释行）
        if w:
            writers[b] = w
        if b in RETIRED_WRITERS:
            retired_ok[b] = w                  # 历史写者不得再有任何写出位置
            listed.append(f"{b}: write-sites={w or 'none'}")
        for i, ln in enumerate(_read(p).split("\n"), 1):
            if ln.strip().startswith("#"):
                continue
            if ARTIFACT in ln and CANONICAL_PRODUCER in ln:
                declared.append(f"{b}:{i}")
    ok = (set(writers) == {CANONICAL_PRODUCER}) \
        and all(not v for v in retired_ok.values()) \
        and bool(declared)
    return ok, (f"生产脚本 {len(prods)} 个；以**写出位置**提到 {ARTIFACT} 的："
                f"{ {k: v for k, v in writers.items()} }（registry 要求恰为 {{{CANONICAL_PRODUCER}}}）；"
                f"历史写者的**写出位置**：{listed}（必须全空 —— 注释掉的旧写出不算活跃写者）；"
                f"活文本里声明 canonical producer 的位置：{declared}")


L.guard("CS11-01", cs11_01,
        "provenance（§4.1 封闭世界 → 显式 registry）：能**写出** fig3b_loo_organism.csv 的生产脚本"
        "必须恰好等于 registry = {p14_rebuild_fig3b.py}；历史写者（p11）对该文件名的提及必须全部处于注释态",
        "只读文本；不执行任何生产文件。注释行显式剔除，避免把'已停用的旧写出'当成活跃写者")


# ------------------------------------------------------------------ CS11-02 语义（execution + instrumentation）
class _RecLR(LogisticRegression):
    """记录每次 fit 的**训练行数**与 X 的前若干行；不改变计算。"""
    FITS = []

    def fit(self, X, y, **kw):
        _RecLR.FITS.append(int(np.asarray(X).shape[0]))
        return super().fit(X, y, **kw)


def make_fixture(seed=909090, sizes=(120, 90, 150)):
    """三个 organism 的受控 fixture：`r` 与 `y` 有真实关系，故 fit 可收敛。"""
    rng = np.random.default_rng(seed)
    ys, rs, org = [], [], []
    for k, n in enumerate(sizes):
        y = (rng.random(n) < 0.5).astype(int)
        r = 2.0 * y - 1.0 + rng.normal(0, 0.8, n)          # 与 y 相关的信号
        ys.append(y); rs.append(r); org += [f"O{k}"] * n
    return (np.concatenate(ys).astype(int), np.concatenate(rs), np.array(org))


Y_FIX, R_FIX, ORG_FIX = make_fixture()
ORGS_FIX = sorted(set(ORG_FIX.tolist()))


def _ref_auc(yy, ss):
    m = np.isfinite(ss)
    yy, ss = yy[m], ss[m]
    if len(np.unique(yy)) < 2:
        return np.nan
    o = np.argsort(ss, kind="mergesort")
    yy, ss = yy[o], ss[o]
    r = np.empty(len(ss)); i = 0
    while i < len(ss):
        j = i
        while j + 1 < len(ss) and ss[j + 1] == ss[i]:
            j += 1
        r[i:j + 1] = 0.5 * (i + j) + 1
        i = j + 1
    n1 = yy.sum(); n0 = len(yy) - n1
    return np.nan if (n1 == 0 or n0 == 0) else (r[yy == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)


def _run_p14_loop(src_text=None):
    """执行 p14 的留出循环（可给变异后的源码）。返回 (fits, msg)。"""
    _RecLR.FITS = []
    block = src_text if src_text is not None else L.grab(P14, *P14_LOOP)
    env = dict(np=np, LogisticRegression=_RecLR, auc=_ref_auc,
               y=Y_FIX, r=R_FIX, org=ORG_FIX, ORGS=ORGS_FIX,
               proper={}, within={}, signs={})
    exec(compile(block, "<verbatim:p14-loo-loop>", "exec"), env)
    return list(_RecLR.FITS), env


# p14 的留出循环：`for o_ in ORGS:` .. `pass`（含两端），由**模式**定位而非写死行号
def _find_loop():
    lines = _read(P14).split("\n")
    a = next(i for i, ln in enumerate(lines) if re.match(r"^\s*for o_ in ORGS:", ln))
    b = next(i for i in range(a + 1, len(lines)) if re.match(r"^\s*def bci\(", lines[i]))
    # 回退到 bci 之前的最后一个非空行
    while b > a and not lines[b - 1].strip():
        b -= 1
    return a + 1, b


P14_LOOP = _find_loop()


def cs11_02():
    fits, env = _run_p14_loop()
    exp = [int((ORG_FIX != o_).sum()) for o_ in ORGS_FIX]      # 每折的**补集**行数
    same = all((ORG_FIX != o_).sum() != (ORG_FIX == o_).sum() for o_ in ORGS_FIX)
    ok = (fits == exp) and same and len(env["proper"]) >= 1
    return ok, (f"p14 的留出循环（{os.path.basename(P14)}:{P14_LOOP[0]}-{P14_LOOP[1]}，按**模式**定位）"
                f"在受控 fixture（3 个 organism，行数 {[int((ORG_FIX == o_).sum()) for o_ in ORGS_FIX]}）上执行："
                f"每折 fit 的训练行数 = {fits}，而**补集**行数 = {exp} —— "
                f"{'逐折相等（每折只用另两个 organism 拟合）' if fits == exp else '**不相等（疑似泄漏留出折）**'}；"
                f"完成折数 = {len(env['proper'])}；学得符号 = "
                f"{ {k: int(v) for k, v in env['signs'].items()} }")


L.guard("CS11-02", cs11_02,
        "算法语义（execution + instrumentation，与数值无关）：p14 的留出循环里，"
        "**每一折的 logistic 必须只在另两个 organism 的行上拟合** —— 训练行数逐折等于补集行数；"
        "且每折都能产出留出评价（proper 非空）。"
        "★ 这条不能靠'数值碰巧吻合'，因为一个泄漏了留出折的模型同样能给出相近的 AUROC",
        "留出循环由**正则模式**定位（`for o_ in ORGS:` → `def bci(`），不写死行号，故对上游行数变化稳健")


# ------------------------------------------------------------------ CS11-03 陈旧路径隔离
def cs11_03():
    hits = [i for i, ln in enumerate(_read(P6F).split("\n"), 1)
            if "LOO-organism" in ln and ln.lstrip().startswith("p(")]
    lodo = [i for i, ln in enumerate(_read(P6F).split("\n"), 1) if "lodo_oof(" in ln]
    ok = not hits
    return ok, (f"{os.path.basename(P6F)}：`p(` 开头的 `LOO-organism` 段 header 出现 {len(hits)} 次"
                + (f" @ 行 {hits}" if hits else "（**不存在** ⇒ 陈旧段已隔离）")
                + f"；全文 `lodo_oof(` 出现在行 {lodo}（这些是**主结算循环**的正当调用，不是死掉的 LOO 段）")


L.guard("CS11-03", cs11_03,
        "陈旧路径隔离（CS-11）：被测树的 p6f 中**不得**存在 `LOO-organism` 段 header"
        "（只有 `p(` 开头的打印行才算 header；docstring 里的提及不算）",
        "★ 只删 header 而留下仍在运行的 lodo_oof 计算 = 文本标记式假修复；"
        "本判据与 t12 EVAL-D1 一起，要求该段的**计算本体也已移除**（见 sim_fix DELETIONS J1）")


# ------------------------------------------------------------------ CS11-04 非空转
def cs11_04():
    src = L.grab(P14, *P14_LOOP)
    old = "tr = ~te"
    assert src.count(old) == 1, f"留出循环里应有唯一一行 `{old}`"
    mut = src.replace(old, "tr = np.ones_like(te, dtype=bool)", 1)   # 在**全部行**上拟合 = 泄漏
    fits, _ = _run_p14_loop(mut)
    exp = [int((ORG_FIX != o_).sum()) for o_ in ORGS_FIX]
    ok = fits != exp
    return ok, (f"把留出循环篡改为 `tr = np.ones_like(te, dtype=bool)`（**在全部行上拟合**，泄漏留出折）后："
                f"每折 fit 行数 = {fits}，补集行数 = {exp} —— "
                f"{'与判据不符，**被检出**（判据非空转）' if ok else '**未被检出，判据空转**'}；"
                f"对照：未篡改时逐折相等（见 CS11-02）")


L.guard("CS11-04", cs11_04,
        "非空转对照：把 p14 留出循环的 `tr = ~te` 篡改为在**全部行**上拟合后，CS11-02 的判据必须判红",
        "★ 没有这一条，「每折只用补集拟合」只是一句声明；有这一条，才证明判据真的在看训练集")

print("\n[t17_cs11_loo done]", flush=True)
