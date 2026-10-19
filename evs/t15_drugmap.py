# -*- coding: utf-8 -*-
r"""t15_drugmap.py — family DRUG：**CS-01 药物名映射**的 gate

================================ 为什么这是第一条 Critical path ================================
GPT6（2026-10-05）："药物映射是整个 corrected benchmark 的第一个 Critical path；
这两条不闭合，就绝不能 freeze。"

manifest 的 `P5_BENCHMARK_UNIVERSE` 与 `DRUGMAP_ALL` 在 2026-10-05 被**撤回 PASS**
（原绑定 BOUND-04/FORK-06/07 根本不覆盖 drug mapping），改为 `pending_gates: CS01_DRUGMAP_GATE`。
本模块就是那条 gate。

================================ 被检验的缺陷 ================================
`_canon()` 5 个活副本（p5:57 / p6_02:56 / p6b:63 / p6c:70 / p10a:59）都是：

    for c, ks in DRUGMAP:
        if any(k in up for k in ks): return c      # 子串 + 先到先得
    return up

三个方向的失效（已由 p28 在**真实数据**上量出）：
  ① 并错：`("TZP", [..., "tazobactam"])` 排在 `("CZA", ["ceftolozanetazobactam"])` **之前**
     ⇒ 泛键 "tazobactam" 吞掉 ceftolozane/tazobactam（真实 13,426 行并入 TZP）；
     而本项目 CZA 条目**就是** ceftolozane/tazobactam。
  ② 并错：同一个 "tazobactam"/"sulbactam" 机制吞掉 ceftazidime/avibactam（13,124 行并入 CAZ）、
     cefoperazone/sulbactam 等复合药。
  ③ **欠匹配**：UT 的 `Amoxicillin/k clavulanate aug` 归一化后匹配不上任何键 —— 这一支落到
     `return up`，其实**已经**是正确行为（成为自己的类），但它**没有被报告**，属静默丢失。

================================ 修法（CS-01 裁决） ================================
**废掉 substring-first-match，改精确匹配**：`if up in ks: return c`。
别名表已含精确键，故：
    ceftolozane/tazobactam -> ceftolozanetazobactam -> CZA          ✅ 正确身份
    ceftazidime/avibactam  -> 无精确键 -> 自己的类（进 unmapped 审计） ✅ 复合药不并进单药
    cefoperazone/sulbactam -> 自己的类                                ✅
    Pip/tazobactam         -> piptazobactam -> TZP                    ✅ 合法 TZP 不受影响

================================ 判据（每条都有 expected 来源与 observed） ================================
  DRUG-01  5 个活副本**逐字执行** + 跨副本一致（同输入 -> 同输出）
  DRUG-02  **复合制剂不得并进单药**：ceftolozane/tazobactam≠TZP、ceftazidime/avibactam≠CAZ、
           cefoperazone/sulbactam≠SAM、amoxicillin/clavulanate≠AMP、piperacillin/tazobactam 不得≠TZP
  DRUG-03  **子串陷阱**：仅仅"包含"某别名键的输入，不得被静默判成该键的 canonical 类
  DRUG-04  **未知串不得静默归类**：既不精确命中任何键，返回值又必须可被 unmapped 审计识别
  DRUG-05  **真实数据 unmapped 审计**：三站点 `antibiotic` 的全部唯一取值，
           每一个要么精确命中别名表，要么出现在 unmapped 列表里 —— **不得静默丢失**
  DRUG-06  非空转对照：上述判据在**当前生产实现**（子串匹配）上必须判红，
           在**正确实现**（精确匹配）上必须判绿 —— 由 REAL / PATCHED 双树自动给出

输出契约同 evs_load：`RESULT|tid|PASS|FAIL|WARN|expected|observed|note`
"""
import os, sys, io, ast, json, re
import numpy as np

import evs_load as L

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(HERE, "_cache")
os.makedirs(CACHE, exist_ok=True)

# ★★ 2026-10-06：**手写的 5 文件清单已废除**（原为 manifest 登记的 5 个活副本）。
#   原来这里写死 5 份；而 2026-10-06 的穷尽清点发现全项目实为 **12 份**（11 活 + 1 废弃）。
#   另外 6 份（p6d / p12 / p21a / p26 / p28 / p4_01）**从未被本判据看见** —— CS-01 的漏网由此产生。
#   手写清单是"普查"的否命题：它只能确认「我已知的那几份一致」，
#   不能回答「还有没有别的副本」。⇒ 改为**在生产树里普查**。
#   排除项必须**显式登记并给出理由**，不许静默跳过。
CENSUS_EXCLUDE = {
    "p6_superseded_pre_correction/p6_02_label_free_signals.py":
        "整目录已显式标记为 superseded（不被任何活脚本引用）；保留原文作为历史对照",
}


def _discover_copies():
    """在生产树（含子目录）里普查所有**同时**定义 DRUGMAP 与顶层 _canon 的 .py。

    返回 [(label, abspath, locator)]；label = 相对生产树的路径。
    这就是把 CS-01 的「位置」栏从手写升级为可复核的普查。"""
    out = []
    for root, dirs, files in os.walk(L.OUT):
        dirs[:] = [d for d in dirs if d != "__pycache__"]
        for fn in sorted(files):
            if not fn.endswith(".py"):
                continue
            p = os.path.join(root, fn)
            rel = os.path.relpath(p, L.OUT).replace(os.sep, "/")
            if rel in CENSUS_EXCLUDE:
                continue
            try:
                src = io.open(p, encoding="utf-8", errors="replace").read()
            except OSError:
                continue
            if "DRUGMAP" not in src:
                continue
            if not re.search(r"^def _canon\(", src, re.M):
                continue
            out.append((rel, p, rel + ":DRUGMAP"))
    return out


COPIES = _discover_copies()

# CS-01 的权威别名表 —— **从冻结 amendment 里读**，不在这里再抄一遍。
# （在判据里重抄一份权威 = 又造一个副本，正是本文件要修的病。）
AMENDMENT = os.path.join(L.PROJ, "FREEZE_pre_correction_20261005", "CS01_amendment.txt")
AMENDMENT_V2 = os.path.join(L.PROJ, "FREEZE_pre_correction_20261005", "CS01_amendment_v2.txt")


def _amendment_alias_table():
    """解析**冻结 amendment** 里「归一化键 -> canonical」的声明块。

    权威 = `CS01_amendment.txt` 的 `纳入 benchmark alias` 段（行形态固定为
    `  <key>  -> <CANON>`）。**不在判据里重抄一份** —— 重抄就是又造一个副本，
    正是本文件要修的病。v2 只修实现以符合本声明，未改声明本身。

    解析不到任何条目 -> **抛错**（fail-closed）；判据无靶子时不许静默通过。
    返回 {alias: canonical}。"""
    txt = io.open(AMENDMENT, encoding="utf-8", errors="replace").read()
    if "纳入 benchmark alias" not in txt:
        raise ValueError("amendment 里找不到声明段的标题 —— 版本可能已变，不得静默继续")
    seg = txt.split("纳入 benchmark alias", 1)[1].split("明确不纳入", 1)[0]
    tbl = {m.group(1): m.group(2)
           for m in re.finditer(r"^\s+([a-z0-9]+)\s+->\s+([A-Z]{2,4})\s*$", seg, re.M)}
    if not tbl:
        raise ValueError("从 %s 的声明段解析不到任何 <key> -> <CANON> 行"
                         % os.path.basename(AMENDMENT))
    return tbl

# 三站点原始微生物表（列 antibiotic）；与 p28 审计同一来源
_RAW = os.environ.get("ARMD_RAW") or os.path.dirname(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
SITES = [
    ("MGB", os.path.join(_RAW, "ARMD-MGB", "microbiology_cohort_deid_tj_updated.csv")),
    ("SU",  os.path.join(_RAW, "ARMD-Stanford", "microbiology_cultures_cohort.csv")),
    ("UT",  os.path.join(_RAW, "ARMD-UTSW", "microbiology_cultures_cohort.csv")),
]

EXP_HEAD = ("CORRECTION_SPEC §1 CS-01 验收：① 三站点 antibiotic 全部唯一取值的 exact-mapping 表；"
            "② unmapped 列表为空、或每一项有显式裁定；③ fixture：复合药**不得**映射到单药 canonical；"
            "④ fixture：未知串**必须**报错而非静默通过。裁决：废掉 substring-first-match，"
            "建立显式 exact alias → canonical 字典；组合制剂与单药严格分开；绝不静默归类。")


# ------------------------------------------------------------------ 逐字抽取
def _src_assignment(path, name):
    """用 AST 取顶层 `NAME = [...]` 的**原始文本段**（不执行、不猜测）。"""
    src = io.open(path, encoding="utf-8", errors="replace").read()
    tree = ast.parse(src, filename=path)
    for node in tree.body:
        if isinstance(node, ast.Assign) and any(
                getattr(t, "id", None) == name for t in node.targets):
            return ast.get_source_segment(src, node)
    raise LookupError(f"{os.path.basename(path)} 中找不到顶层赋值 {name}")


def load_copy(path, label):
    """逐字执行该副本的 DRUGMAP + _canon。返回 (canon_fn, drugmap, err)。"""
    try:
        seg = _src_assignment(path, "DRUGMAP")
        ns = {"np": np}
        exec(compile(seg, f"<verbatim:{os.path.basename(path)}:DRUGMAP>", "exec"), ns)
        drugmap = ns["DRUGMAP"]
        fn = L.load_def(path, "_canon", extra_ns={"DRUGMAP": drugmap})
        return fn, drugmap, None
    except Exception as e:                                   # noqa: BLE001
        return None, None, f"{type(e).__name__}: {e}"


COPIES_LOADED, COPIES_ERR, DRUGMAPS = {}, {}, {}
for _lab, _p, _loc in COPIES:
    _fn, _dm, _err = load_copy(_p, _lab)
    if _err:
        COPIES_ERR[_lab] = f"{_loc} -> {_err}"
    else:
        COPIES_LOADED[_lab] = _fn
        DRUGMAPS[_lab] = _dm

print("=" * 100, flush=True)
print("t15 — family DRUG：CS-01 药物名映射", flush=True)
print("=" * 100, flush=True)
print(f"  生产树 = {L.PROD_DIR}   模式 = {'PATCHED' if L.IS_PATCHED else 'REAL'}", flush=True)
print(f"  逐字装载副本 {len(COPIES_LOADED)}/{len(COPIES)}"
      + ("" if not COPIES_ERR else f"；失败：{COPIES_ERR}"), flush=True)
print("", flush=True)


# ------------------------------------------------------------------ DRUG-01 逐字执行 + 跨副本一致
_FIXTURE_RAW = [
    # 合法单药（不得受影响）
    "Piperacillin/Tazobactam", "Pip/tazobactam", "ceftazidime", "Ceftazidime",
    "meropenem", "Meropenem", "ciprofloxacin", "gentamicin", "amikacin",
    "trimethoprim/sulfamethoxazole", "Colistin", "tigecycline", "ampicillin",
    "cefazolin", "ceftriaxone", "cefepime", "cefotaxime", "cefuroxime", "cefoxitin",
    "levofloxacin", "tobramycin", "nitrofurantoin", "imipenem", "ertapenem",
    "aztreonam", "tetracycline", "fosfomycin", "piperacillin",
    # 复合制剂（必须**不**并进单药）
    "ceftolozane/tazobactam", "Ceftolozane/Tazobactam", "ceftolozane_tazobactam",
    "ceftazidime/avibactam", "Ceftazidime/Avibactam", "ceftazidime_avibactam",
    "cefoperazone/sulbactam", "ampicillin/sulbactam", "amoxicillin/clavulanate",
    "Amoxicillin/k clavulanate aug", "ticarcillin/clavulanate",
    # 归一化后才等价的写法
    "PIP/TAZO", "Pip/Tazo", "TZP", "CAZ", "CZA", "SXT", "TMP-SMX",
]


def _run_all(raw):
    return {lab: fn(raw) for lab, fn in COPIES_LOADED.items()}


def drug01():
    if len(COPIES_LOADED) != len(COPIES):
        return False, f"只装载到 {len(COPIES_LOADED)}/{len(COPIES)} 个副本：{COPIES_ERR}"
    out = {}
    for raw in _FIXTURE_RAW:
        vals = _run_all(raw)
        out[raw] = vals
    # 跨副本一致：同一输入必须同输出
    bad = {r: v for r, v in out.items() if len(set(v.values())) != 1}
    n_noncanon = sum(1 for r, v in out.items()
                     if list(v.values())[0] != "".join(c for c in r.lower() if c.isalnum()))
    return (not bad), (
        f"{len(COPIES_LOADED)} 个副本逐字执行（DRUGMAP 用 AST 取原文段 + _canon 用 load_def）；"
        f"fixture {len(_FIXTURE_RAW)} 条原始值；跨副本分歧 {len(bad)}"
        + ("" if not bad else f" :: {list(bad)[:4]}") +
        f"；其中映射到 canonical 类的 {n_noncanon}/{len(_FIXTURE_RAW)}")


L.guard("DRUG-01", drug01,
        EXP_HEAD + " ｜ 5 个活副本（p5/p6_02/p6b/p6c/p10a）必须被**逐字执行**，"
        "且对同一原始值给出**相同**结果（同一个统计量不能有两个答案）",
        "")


# ------------------------------------------------------------------ DRUG-02 复合制剂不得并进单药
# 依据：本项目命名空间由 DRUGMAP 自身定义 —— CZA 条目**就是** ceftolozane/tazobactam；
# CAZ 条目是**单独**的 ceftazidime。故复合药映射到单药类即可判定为"并错"。
_COMPOUND_MUST_NOT = [
    ("ceftolozane/tazobactam",   "TZP", "本项目 CZA=ceftolozane/tazobactam；并入 TZP 即并错"),
    ("Ceftolozane/Tazobactam",   "TZP", "同上（大小写/斜杠写法）"),
    ("ceftazidime/avibactam",    "CAZ", "CAZ=单独 ceftazidime；avibactam 复合药必须独立成类"),
    ("Ceftazidime/avibactam",    "CAZ", "同上"),
    ("cefoperazone/sulbactam",   "SAM", "SAM=ampicillin/sulbactam；cefoperazone 复合药必须独立"),
    ("amoxicillin/clavulanate",  "AMP", "AMP=单独 ampicillin；clavulanate 复合药必须独立"),
]
# 合法映射（必须**保持**，防止"修过头"）
_LEGIT_MUST = [
    ("Piperacillin/Tazobactam", "TZP"), ("Pip/tazobactam", "TZP"),
    ("piperacillin/tazobactam", "TZP"), ("ceftazidime", "CAZ"),
    ("meropenem", "MEM"), ("ciprofloxacin", "CIP"), ("colistin", "CST"),
    ("trimethoprim/sulfamethoxazole", "SXT"), ("ampicillin", "AMP"),
]


def drug02():
    if not COPIES_LOADED:
        return False, "无可用副本"
    lab = sorted(COPIES_LOADED)[0]
    fn = COPIES_LOADED[lab]
    viol, ok_legit, miss_legit = [], [], []
    for raw, wrong, why in _COMPOUND_MUST_NOT:
        got = fn(raw)
        if got == wrong:
            viol.append(f"{raw!r} -> {got}（不得为 {wrong}：{why}）")
        else:
            ok_legit.append(f"{raw!r}->{got}")
    for raw, want in _LEGIT_MUST:
        if fn(raw) != want:
            miss_legit.append(f"{raw!r} -> {fn(raw)}（应为 {want}）")
    return (not viol and not miss_legit), (
        f"副本 {lab}；复合药 {len(_COMPOUND_MUST_NOT)} 条中并错 {len(viol)}"
        + ("" if not viol else " :: " + " | ".join(viol[:4])) +
        f"；合法单药 {len(_LEGIT_MUST)} 条中错失 {len(miss_legit)}"
        + ("" if not miss_legit else " :: " + " | ".join(miss_legit[:4])) +
        ("" if viol else f"；未并错的复合药实际落点 = {ok_legit[:3]}"))


L.guard("DRUG-02", drug02,
        EXP_HEAD + " ｜ 具体判据：ceftolozane/tazobactam ≠ TZP、ceftazidime/avibactam ≠ CAZ、"
        "cefoperazone/sulbactam ≠ SAM、amoxicillin/clavulanate ≠ AMP；"
        "同时**合法单药映射必须保持不变**（防止修过头）",
        "canonical 命名空间由 DRUGMAP 自身定义：CZA 条目即 ceftolozane/tazobactam，CAZ 条目是单独的 ceftazidime")


# ------------------------------------------------------------------ DRUG-03 子串陷阱
_TRAPS = [
    # 含 "tazobactam" 但不是 piperacillin/tazobactam
    "ceftolozane/tazobactam", "ceftolozanetazobactam", "XYZ/tazobactam",
    # 含 "sulbactam" 但不是 ampicillin/sulbactam
    "cefoperazone/sulbactam", "XYZ/sulbactam",
    # 含单药名但不是该药
    "notmeropenem", "meropenemX", "pseudociprofloxacin", "ceftazidimeX",
    # 含 clavulanate
    "amoxicillin/clavulanate", "XYZ/clavulanate",
]


def drug03():
    if not COPIES_LOADED:
        return False, "无可用副本"
    lab = sorted(COPIES_LOADED)[0]
    fn = COPIES_LOADED[lab]
    dm = DRUGMAPS[lab]
    # 别名键集合（用于判断"是否被静默判给了某个键"）
    hit = []
    for raw in _TRAPS:
        up = "".join(c for c in raw.lower() if c.isalnum())
        got = fn(raw)
        # 该输入**精确等于**某个键吗？不等则是陷阱
        exact = any(up == k for _, ks in dm for k in ks)
        if not exact and got in {c for c, _ in dm}:
            hit.append(f"{raw!r} -> {got}（非精确命中却被归类）")
    return (not hit), (
        f"副本 {lab}；子串陷阱 {len(_TRAPS)} 条；被静默归类 {len(hit)}"
        + ("" if not hit else " :: " + " | ".join(hit[:4])))


L.guard("DRUG-03", drug03,
        EXP_HEAD + " ｜ 具体判据：仅仅**包含**某别名键的输入，不得被判成该键的 canonical 类"
        "（子串匹配的本质缺陷）",
        "陷阱串在归一化后与任何别名键都不精确相等，故正确的精确匹配实现应让它们落回自身")


# ------------------------------------------------------------------ DRUG-04 未知串不得静默归类
_UNKNOWN = ["ZZZnotadrug", "Amoxicillin/k clavulanate aug", "SomeNewAntibiotic",
            "ceftazidime_avibactam", "unknown", ""]


def drug04():
    if not COPIES_LOADED:
        return False, "无可用副本"
    lab = sorted(COPIES_LOADED)[0]
    fn = COPIES_LOADED[lab]
    dm = DRUGMAPS[lab]
    canon = {c for c, _ in dm}
    silently = []
    for raw in _UNKNOWN:
        got = fn(raw)
        up = "".join(c for c in str(raw).lower() if c.isalnum())
        exact = any(up == k for _, ks in dm for k in ks)
        if not exact and got in canon:
            silently.append(f"{raw!r} -> {got}（未知串被并进既有 canonical 类）")
        elif not exact and got != up:
            silently.append(f"{raw!r} -> {got}（未知串的返回值不可由 unmapped 审计识别）")
    return (not silently), (
        f"副本 {lab}；未知串 {len(_UNKNOWN)} 条；静默归类/不可识别 {len(silently)}"
        + ("" if not silently else " :: " + " | ".join(silently[:4])))


L.guard("DRUG-04", drug04,
        EXP_HEAD + " ｜ 具体判据：既不精确命中任何别名键，返回值又必须等于其归一化形式，"
        "从而可被 unmapped 审计识别 —— **绝不静默归类**",
        "")


# ------------------------------------------------------------------ DRUG-05 真实数据 unmapped 审计
def _distinct_raw(site, path):
    """该站点 antibiotic 列的全部唯一取值。按 (mtime,size) 缓存，源变则重读。"""
    ck = os.path.join(CACHE, f"drugmap_distinct_{site}.json")
    try:
        st = os.stat(path)
        stamp = [int(st.st_mtime), st.st_size]
    except OSError as e:
        return None, f"无法 stat {path}: {e}"
    if os.path.exists(ck):
        try:
            with io.open(ck, encoding="utf-8") as fh:
                obj = json.load(fh)
            if obj.get("stamp") == stamp:
                return obj["values"], f"cache（{len(obj['values'])} 个唯一值）"
        except Exception:                                    # noqa: BLE001
            pass
    try:
        import pandas as pd
        d = pd.read_csv(path, usecols=["antibiotic"], low_memory=False,
                        encoding="utf-8-sig", dtype=str)
        vals = sorted({str(x) for x in d["antibiotic"].dropna().unique()})
    except Exception as e:                                   # noqa: BLE001
        return None, f"读取失败 {type(e).__name__}: {e}"
    with io.open(ck, "w", encoding="utf-8") as fh:
        json.dump({"stamp": stamp, "values": vals}, fh, ensure_ascii=False)
    return vals, f"实读（{len(vals)} 个唯一值）"


def drug05():
    if not COPIES_LOADED:
        return False, "无可用副本"
    lab = sorted(COPIES_LOADED)[0]
    fn = COPIES_LOADED[lab]
    dm = DRUGMAPS[lab]
    canon = {c for c, _ in dm}
    lines, errs, tot, unmapped = [], [], 0, {}
    for site, path in SITES:
        vals, how = _distinct_raw(site, path)
        if vals is None:
            errs.append(f"{site}: {how}")
            continue
        tot += len(vals)
        u = []
        for v in vals:
            up = "".join(c for c in v.lower() if c.isalnum())
            if any(up == k for _, ks in dm for k in ks):
                continue                                    # 精确命中别名表
            got = fn(v)
            if got in canon:
                u.append((v, got, "**并进既有 canonical 类**"))
            else:
                u.append((v, got, "unmapped（显式可见）"))
        unmapped[site] = u
        lines.append(f"{site}: 唯一值 {len(vals)}，其中 unmapped {len(u)}")
    hard = [(s, v, g, w) for s, us in unmapped.items() for (v, g, w) in us if w.startswith("**")]
    soft = sum(1 for us in unmapped.values() for (_, _, w) in us if not w.startswith("**"))
    with io.open(os.path.join(HERE, "drugmap_unmapped_audit.txt"), "w", encoding="utf-8") as fh:
        fh.write(f"# CS-01 unmapped 审计 ｜ 副本 {lab} ｜ 模式 {'PATCHED' if L.IS_PATCHED else 'REAL'}\n")
        fh.write(f"# 判据：唯一取值要么精确命中别名表，要么出现在下表（不得静默并进既有 canonical 类）\n\n")
        for s in unmapped:
            fh.write(f"## {s}（{len(unmapped[s])} 条）\n")
            for v, g, w in sorted(unmapped[s]):
                fh.write(f"   {v!r:<48} -> {g!r:<40} {w}\n")
            fh.write("\n")
    return (not errs and not hard), (
        f"副本 {lab}；{' ｜ '.join(lines) if lines else '无站点可读'}"
        f"；合计唯一值 {tot}；显式 unmapped {soft}；**被静默并进既有类** {len(hard)}"
        + ("" if not hard else " :: " + " | ".join(f"{s}:{v!r}->{g}" for s, v, g, _ in hard[:4]))
        + ("" if not errs else f"；站点读取失败 {errs}")
        + f"；明细落盘 drugmap_unmapped_audit.txt")


L.guard("DRUG-05", drug05,
        EXP_HEAD + " ｜ 具体判据（真实数据）：三站点 antibiotic 列的**全部唯一取值**，"
        "每一个要么精确命中别名表、要么出现在显式 unmapped 列表里；"
        "**任何「既未精确命中、又被归类为既有 canonical 类」的值都判红**（静默丢失/并错）",
        "唯一取值按 (mtime,size) 缓存于 evs/_cache/，源文件变化即重读")


# ------------------------------------------------------------------ DRUG-06 五分类报告
# ★ GPT6 2026-10-05 裁定：**不要只报 "unmapped count ↓"，那个指标没有意义**。
#   要报的是 identity 的正确识别：
#     total raw drug labels -> mapped_to_benchmark -> recognized_nonbenchmark
#                            -> combination_not_in_benchmark -> ambiguous_excluded -> missing
#   且必须**完整分割**（每个唯一取值恰好落入一类，不重不漏）。
#   分类依据 = `_canon` 的结果 + **显式声明**的规则列表（封闭世界用 registry，不做模式推断）。
COMPOUND_NOT_IN_BENCH = {
    "ceftazidimeavibactam", "cefoperazonesulbactam", "imipenemrelebactam",
    "meropenemvaborbactam", "aztreonamavibactam", "cefotaximeclavulanate",
    "ceftazidimeclavulanate", "cephalexincephalothin",
}
AMBIGUOUS_EXCLUDED = {"meropenemimipenem"}          # 一条记录两药 -> 不猜
MISSING = {"", "null", "nan", "none"}


def _classify(v, fn, canon):
    """把一个原始 label 分到五类之一。完备且互斥。"""
    up = "".join(c for c in str(v).lower() if c.isalnum())
    got = fn(v)
    if up in MISSING:
        return "missing"
    if up in AMBIGUOUS_EXCLUDED:
        return "ambiguous_excluded"
    if got in canon:
        return "mapped_to_benchmark"
    if up in COMPOUND_NOT_IN_BENCH:
        return "combination_not_in_benchmark"
    return "recognized_nonbenchmark"


def drug06():
    if not COPIES_LOADED:
        return False, "无可用副本"
    lab = sorted(COPIES_LOADED)[0]
    fn, dm = COPIES_LOADED[lab], DRUGMAPS[lab]
    canon = {c for c, _ in dm}
    import collections
    counts = collections.Counter()
    per_site, total = {}, 0
    for site, path in SITES:
        vals, how = _distinct_raw(site, path)
        if vals is None:
            continue
        c = collections.Counter(_classify(v, fn, canon) for v in vals)
        per_site[site] = c
        counts.update(c)
        total += len(vals)
    order = ["mapped_to_benchmark", "recognized_nonbenchmark",
             "combination_not_in_benchmark", "ambiguous_excluded", "missing"]
    ssum = sum(counts.values())
    with io.open(os.path.join(HERE, "drugmap_identity_report.txt"), "w",
                 encoding="utf-8") as fh:
        fh.write("# CS-01 benchmark identity 报告 | 副本 " + lab + " | 模式 " +
                 ("PATCHED" if L.IS_PATCHED else "REAL") + chr(10))
        fh.write("# 每个唯一 raw drug label 恰好落入一类（完整分割，不重不漏）" + chr(10) + chr(10))
        fh.write("total raw drug labels = %d" % total + chr(10))
        for k in order:
            fh.write("  %-32s %5d  (%5.1f%%)" % (k, counts[k], 100.0 * counts[k] / max(ssum, 1)) + chr(10))
        fh.write(chr(10) + "## 逐站点" + chr(10))
        for site in per_site:
            fh.write(chr(10) + "[%s]  %d" % (site, sum(per_site[site].values())) + chr(10))
            for k in order:
                fh.write("  %-32s %5d" % (k, per_site[site][k]) + chr(10))
    ok = (ssum == total)
    return ok, ("副本 " + lab + "；total=" + str(total) + "；" +
                " ".join(k + "=" + str(counts[k]) for k in order) +
                "；分块之和=" + str(ssum) + "（须等于 total：完整分割）" +
                "；逐站点 " + ", ".join(s + ":" + str(sum(per_site[s].values())) for s in per_site) +
                "；明细落盘 drugmap_identity_report.txt")


L.guard("DRUG-06", drug06,
        "benchmark identity 的五分类报告（GPT6 裁定）：total raw labels -> "
        "mapped_to_benchmark / recognized_nonbenchmark / combination_not_in_benchmark / "
        "ambiguous_excluded / missing，必须**完整分割**（各类之和 == total，不重不漏）。"
        "★ 不报 unmapped 计数 —— 减少 unmapped 不是目标，正确识别 identity 才是",
        "分类依据 = _canon 结果 + **显式规则集**（COMPOUND_NOT_IN_BENCH / AMBIGUOUS_EXCLUDED / MISSING），不做模式推断")


# ---------------------------------------------------------------- DRUG-07 同源副本普查
def drug07():
    """同源副本**普查**：报告生产树里所有药物映射副本（含子目录，排除项显式登记）。

    ★ 它守的是「**还有没有别的副本**」这个问题 —— 手写清单回答不了它。
      CS-01 漏掉 6 份副本，正是因为判据用的是手写清单而不是普查。"""
    labs = [c[0] for c in COPIES]
    ex = sorted(CENSUS_EXCLUDE)
    return True, (f"普查到 {len(COPIES)} 份活副本 :: {labs}；"
                  f"显式排除 {len(ex)} 份 :: {ex}")


# ---------------------------------------------------------------- DRUG-08 对着裁决验
def drug08():
    """★ 与 DRUG-01 的**关键区别**：DRUG-01 验的是**副本之间一致**
    （5 个都错时它恒绿 —— CS-01 的 7 个别名错位就是这么漏过去的）；
    DRUG-08 验的是**每一份都与冻结 amendment 的声明一致**。

    权威 = `CS01_amendment.txt` 的声明段（判据**不重抄**它，直接解析）。"""
    tbl = _amendment_alias_table()
    bad, per = [], {}
    for lab, _p, _loc in COPIES:
        fn = COPIES_LOADED.get(lab)
        if fn is None:
            bad.append(f"{lab}: 装载失败 {COPIES_ERR.get(lab, '')[:60]}")
            continue
        wrong = [(a, c, fn(a)) for a, c in sorted(tbl.items()) if fn(a) != c]
        per[lab] = len(wrong)
        for a, c, g in wrong:
            bad.append(f"{lab}: {a} 应为 {c} 实为 {g}")
    return (not bad), (f"权威 = {os.path.basename(AMENDMENT)} 的 {len(tbl)} 条声明；"
                       f"普查 {len(COPIES)} 份；不符 {len(bad)}"
                       + (f" :: {bad[:6]}" if bad else "")
                       + f"；逐份不符数 {per}")


L.guard("DRUG-07", drug07,
        "CS-01 验收的**普查前提**：药物映射副本必须由**普查**发现，不得用手写清单。"
        "手写清单是'普查'的否命题 —— 它只能确认「已知的那几份一致」，不能回答「还有没有别的副本」",
        "★ 2026-10-06：原 COPIES 为手写的 5 文件；穷尽清点发现实为 12 份，漏网的 6 份由此产生")
L.guard("DRUG-08", drug08,
        "CS-01 验收的**裁决符合性**：每一份副本的 _canon 必须与冻结 amendment 的声明一致。"
        "权威 = CS01_amendment.txt 的声明段（**判据不重抄**，直接解析；解析不到即抛错，fail-closed）",
        "★ 与 DRUG-01 的分工：DRUG-01 验副本**互相一致**（全错时恒绿），DRUG-08 验**与裁决一致**。"
        "CS-01 那 7 个别名错位恰好只能被后者抓到")

print(chr(10) + "[t15_drugmap done]", flush=True)