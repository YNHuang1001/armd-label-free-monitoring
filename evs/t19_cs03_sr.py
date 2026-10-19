# -*- coding: utf-8 -*-
r"""t19_cs03_sr.py — family CS03：episode × drug 的 **S/R 冲突裁决** gate

================================ 性质声明（GPT6 2026-10-05 裁定） ================================
**这不是新的 audit family，而是 correction phase 内、为已冻结的 correction item
（CORRECTION_SPEC §1 CS-03）补齐其必要验收条件。**

    "不要写『新增 CS-03 审计模块』，而应该写『CS-03 correction item validation amendment』。
     原因：CS-03 在冻结前已经存在于 correction plan，现在只是发现它的验收覆盖不足。
     这和 census 重构 / manifest 重做性质不同 —— 那是验证架构错误，CS-03 是**验证缺口**。"

CS-01/02/03 共同构成 **BENCHMARK IDENTITY**：纳入哪些 episode、哪个 episode-drug 对进入分析、
哪个 S/R 结果作为最终标签。CS-01 有 gate、CS-02 有 gate，**CS-03 没有** ⇒ 冻结原则要求补上，
否则 benchmark identity 会变成"95% 审计过，5% 信我"。

================================ 被检验的语义（CORRECTION_SPEC §1 CS-03 裁决） ================================
1. **同日冲突消解**：`S+I` → **S**；`R+I` → **R**；**仅当同日同时存在 `S+R` 才置 ambiguous 并排除**
2. **跨日**：取 episode 内「**最早一个携带直接 S/R 的日期**」——I-only 的日期**跳过**
3. **后续日期的结果不得回填**（Day1=S、Day10=R ⇒ outcome 必须是 **S**）
4. **对文件行序不变**（可复现）

================================ 判据 ================================
  CS03-01  同日 S+R -> **排除**（既不得为 R，也不得为 S）
  CS03-02  同日 S+I -> S
  CS03-03  同日 R+I -> R
  CS03-04  后续日期的 R **不得**定义 index outcome；I-only 日期必须被跳过
  CS03-05  行序不变性：打乱输入行序 -> 结果逐位相同
  CS03-06  非空转：把冲突判据篡改（使 S+R 不再排除）后，CS03-01 必须判红

★ 被测的是**生产代码本身**：两行（同日消解 / episode 内取最早 S-R）由 `L.grab` 逐字抽出执行。
输出契约同 evs_load：`RESULT|tid|PASS|FAIL|WARN|expected|observed|note`
"""
import io, os, sys, re
import numpy as np
import pandas as pd

import evs_load as L

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

P6B = L.P6B
RAW = io.open(P6B, encoding="utf-8", errors="replace").read().split("\n")


def _find(pat):
    hits = [i + 1 for i, ln in enumerate(RAW) if re.search(pat, ln)]
    assert len(hits) == 1, f"p6b 中模式 {pat!r} 命中 {len(hits)} 次"
    return hits[0]


L_RESOLVE = _find(r'g = usr\.groupby\(\["anon_id", "cdate", "spec", "abx"\]\)\["s"\]')
# ★ 定位必须**形态无关**：CS-03 恰恰要改这一行（前面加 `g[g["s"].isin(["S","R"])]`），
#   第一版写死 `wide = g.pivot_table(`，改完就找不到 —— 与之前"锚点钉住被修正的旧文本"同一类错。
L_PIVOT = _find(r'^\s*wide = .*pivot_table\(index="epid", columns="abx", values="s"')

print("=" * 100, flush=True)
print("t19 — family CS03：episode × drug 的 S/R 冲突裁决", flush=True)
print("=" * 100, flush=True)
print(f"  生产树 = {L.PROD_DIR}   模式 = {'PATCHED' if L.IS_PATCHED else 'REAL'}", flush=True)
print(f"  被测行：同日消解 p6b:{L_RESOLVE} ｜ episode 内经 pivot 取最早 p6b:{L_PIVOT}", flush=True)
print("", flush=True)

EPID = "E1"
DRUG = "GEN"


def _mk(rows):
    return pd.DataFrame(rows, columns=["anon_id", "cdate", "spec", "abx", "s"])


def run_SR(rows, mutate_resolve=None):
    """逐字执行两行生产代码；返回 wide 表。"""
    usr = _mk(rows)
    src_a = L.grab(P6B, L_RESOLVE, L_RESOLVE)
    if mutate_resolve:
        src_a = mutate_resolve(src_a)
    ns = {"usr": usr, "pd": pd, "np": np}
    exec(compile(src_a, "<verbatim:p6b-resolve>", "exec"), ns)
    g = ns["g"]
    g = g.assign(epid=EPID)
    ns2 = {"g": g, "pd": pd, "np": np}
    exec(compile(L.grab(P6B, L_PIVOT, L_PIVOT), "<verbatim:p6b-wide>", "exec"), ns2)
    return ns2["wide"]


def _cell(wide, drug=DRUG):
    try:
        v = wide.loc[EPID, drug]
        return None if (v is None or (isinstance(v, float) and np.isnan(v))) else str(v)
    except (KeyError, IndexError):
        return None


def _r(date, s, drug=DRUG):
    return ("a1", pd.Timestamp(date), "U", drug, s)


# ------------------------------------------------------------------ CS03-01 同日 S+R
def cs03_01():
    w = run_SR([_r("2020-01-01", "S"), _r("2020-01-01", "R")])
    v = _cell(w)
    ok = v is None
    return ok, (f"同 episode 同日 同药 同时有 S 与 R -> 该 episode-drug 的结果 = {v!r}；"
                f"应然 = 排除（None），**不得为 'R'、也不得为 'S'**（ambigance 必须显式排除而非由行序裁决）")


L.guard("CS03-01", cs03_01,
        "CORRECTION_SPEC §1 CS-03：**仅当同一日期同时存在 S+R 才置 ambiguous 并排除**；"
        "不得用文件行序或任何优先级规则把它判成 S 或 R",
        "★ 这是本 gate 最关键的一条：它约束的是**输入空间**——哪个 episode-drug 对进入分析")


# ------------------------------------------------------------------ CS03-02 / 03
def cs03_02():
    w = run_SR([_r("2020-01-01", "S"), _r("2020-01-01", "I")])
    v = _cell(w)
    return (v == "S"), f"同日 S+I -> {v!r}；应然 = 'S'（I 的存在不得改变结论）"


L.guard("CS03-02", cs03_02, "CS-03：同日 `S+I` -> 取 **S**", "")


def cs03_03():
    w = run_SR([_r("2020-01-01", "R"), _r("2020-01-01", "I")])
    v = _cell(w)
    return (v == "R"), f"同日 R+I -> {v!r}；应然 = 'R'"


L.guard("CS03-03", cs03_03, "CS-03：同日 `R+I` -> 取 **R**", "")


# ------------------------------------------------------------------ CS03-04 后续日期不得回填 / 跳过 I-only
def cs03_04():
    # ④a 同一 episode 内（相隔 9 天 < 14 天阈值）：Day1=S、Day10=R -> 必须取 **S**
    w1 = run_SR([_r("2020-01-01", "S"), _r("2020-01-10", "R")])
    v1 = _cell(w1)
    # ④b I-only 的日期必须被**跳过**：Day1=I、Day10=R -> 必须取 **R**
    w2 = run_SR([_r("2020-01-01", "I"), _r("2020-01-10", "R")])
    v2 = _cell(w2)
    ok = (v1 == "S") and (v2 == "R")
    return ok, (f"同一 episode（相隔 9 天 < 14 天阈值）："
                f"Day1=S / Day10=R -> {v1!r}（应然 **S**：后续日期的 R 不得回填定义 index outcome）；"
                f"Day1=I / Day10=R -> {v2!r}（应然 **R**：I-only 的日期必须被跳过）")


L.guard("CS03-04", cs03_04,
        "CS-03：取 episode 内「**最早一个携带直接 S/R 的日期**」——"
        "(a) 后续日期的结果**不得**定义 index outcome；(b) I-only 的日期必须**跳过**",
        "★ 这条直接决定 episode identity：若后续 R 能回填，Day1 已给出 S 的 episode 会被改写成 R")


# ------------------------------------------------------------------ CS03-05 行序不变性
def cs03_05():
    base = [_r("2020-01-01", "S"), _r("2020-01-01", "I"), _r("2020-01-10", "R"),
            _r("2020-01-10", "S"), _r("2020-01-20", "R"), _r("2020-01-20", "I")]
    outs = []
    for perm in ([0, 1, 2, 3, 4, 5], [5, 4, 3, 2, 1, 0], [2, 0, 4, 1, 5, 3], [3, 5, 1, 4, 0, 2]):
        w = run_SR([base[i] for i in perm])
        outs.append(_cell(w))
    ok = (len(set(map(repr, outs))) == 1)
    return ok, (f"同一批记录按 {len(outs)} 种行序输入 -> 结果 = {outs}；"
                f"{'逐位相同（可复现）' if ok else '**随行序变化**（不可复现）'}")


L.guard("CS03-05", cs03_05,
        "CS-03 验收⑤：结果对**文件行序不变**（这是 reproducibility，不是统计问题）",
        "★ 现状的 `aggfunc=\"first\"` 对同日冲突按原始文件行序裁决 —— 不可复现")


# ------------------------------------------------------------------ CS03-06 非空转
def cs03_06():
    def mut(src):
        # ★ 变异必须**形态无关**：不依赖生产行的具体写法，而是把它整行换成一个
        #   "R 优先" 的实现 —— 正是 GPT6 点名要防的那种退化（`if "R" in outcomes: return R`）。
        #   第一版把变异写成「把冲突判据置 False」，在**旧**代码上找不到该片段，
        #   于是 gate 在 RED 阶段报"无法构造变异"，那是判据自身的缺陷。
        # ★ 必须与 `L.grab` 的输出同形：grab 会 dedent 成顶层语句，故这里也不带缩进。
        return ('g = usr.assign(_p=usr["s"].map({"R": 0, "S": 1, "I": 2}))'
                '.sort_values(["anon_id", "cdate", "spec", "abx", "_p"], kind="mergesort")'
                '.groupby(["anon_id", "cdate", "spec", "abx"])["s"]'
                '.agg(lambda v: v.iloc[0]).reset_index()')

    w = run_SR([_r("2020-01-01", "S"), _r("2020-01-01", "R")], mutate_resolve=mut)
    v_mut = _cell(w)
    w0 = run_SR([_r("2020-01-01", "S"), _r("2020-01-01", "R")])
    v0 = _cell(w0)
    ok = (v0 is None) and (v_mut is not None)
    return ok, (f"把同日冲突判据篡改为 False（S+R 不再排除）后：结果 = {v_mut!r}"
                f"（未篡改时 = {v0!r}）—— "
                f"{'被检出（判据非空转）' if ok else '**未被检出，判据空转**'}")


L.guard("CS03-06", cs03_06,
        "非空转对照：把同日消解整行替换成 R 优先 的实现后，CS03-01 的判据必须判红"
        "（即 S+R 会退化成某个确定值）",
        "★ GPT6：『否则这个 gate 还是容易退化』——例如有人日后写成 `if \"R\" in outcomes: return R`，必须红")

print("\n[t19_cs03_sr done]", flush=True)
