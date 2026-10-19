# -*- coding: utf-8 -*-
r"""t14_manifest_closure — manifest 驱动的 registry/census（**只读文本，绝不执行生产文件**）

================================ 为什么整族重写（v1.4，2026-10-04） ================================
v1.2/v1.3 的 census 是 **heuristic discovery**：扫 *.py、按惯用式签名认副本、把认出来的
东西 exec 进隔离 namespace 再打分。三轮对抗测试证明这条路**两个方向同时失效**：

  过宽 → 假失败：无关的 clip+log 工具被当成 ent 副本；docstring 里出现四个 token 就被判 fork
  过窄 → 假通过：改名 fork 完全不被发现；SyntaxError 文件不参与判据
  自指 → 假通过：未登记模块的 claim 被计入覆盖；reachability 静默豁免未登记项
  ★ 最严重：**exec 被发现的副本会执行其模块级代码** —— 实测能让 `np.save` 建文件、
    `os.remove` 删文件、`setattr(np, ...)` 污染打分 oracle

所以本模块**换机制，不打补丁**：

  **登记表 = 手写的 `ACTIVE_PATH_MANIFEST.yaml`（生产脚本不得生成/修改）**
  **检查 = 存在性 + 锚点文本可解析 + 集合闭合（闭合在 run_evs.py 的 closure gate 做）**

本模块**不做**：import / exec / runpy / importlib / eval 任何生产文件；
不按签名认副本；不猜"这是不是一个 fork"。**未登记即未登记**，由 manifest 治理。

================================ 判据 ================================
  MC-00  manifest schema 合法（8 个必备字段 + 枚举 + Critical 必须有 required_tests 或 deferred_reason）
  MC-01  每条的 file 在 PROD_DIR / 项目根下存在（glob 至少命中 1）
  MC-02  每条的 anchors 可**按文本**解析（def / lines / regex / glob）
  MC-03  ★ 本模块源码中**不存在**动态执行原语，且**不调用装载器**（结构断言）
  MC-04a ★ 正对照：canary fixture 的模块级副作用**确实可被检出**（证明判据可满足）
  MC-04b ★ 负对照：registry 扫描后副作用**没有发生**（永久 negative test）
  MC-04c 扫描非空转：canary 文件确实被扫到了
  MC-05  manifest 身份（sha256 记录；且**不在 PROD_DIR 内**）
  MC-06  Critical 声明完整性
  MC-07  required_tests 的族前缀均为已知族（抓 test_id 拼写错误）
=============================================================================
"""
import os, sys, io, ast, re, glob, json, hashlib, shutil, tempfile, subprocess

import numpy as np
import yaml

import evs_load as L

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

HERE = os.path.dirname(os.path.abspath(__file__))
PROD = os.path.abspath(L.PROD_DIR)                      # 分析结果/ 或 _patched/
PROJ = os.path.dirname(HERE)           # 第二篇 新方案/
MANIFEST = os.path.join(PROJ, "ACTIVE_PATH_MANIFEST.yaml")
SELF = os.path.abspath(__file__)

KNOWN_FAMILIES = {
    "LOAD", "MMD", "CONC", "L3", "SUDO", "UNIT", "BOUND", "FORK",
    "PSI", "BW", "SENS", "EVAL", "EC", "MC", "DRUG", "CS11", "P12", "FIG5", "CS03",
    "ESS",     # ★ CS-09（2026-10-06）：ESS_P6C 此前 required_tests 里**一条 ESS 判据都没有**，
               #   族名当然也不在表里。加族名是补登记；真正的判据是 ESS-01/02/03。
    "T23",     # ★ 2026-10-07：渲染件的独立验证（t23_render_parity.py）。
               #   起因：p21d/p21e/p21f **不在 closure 集合内**，其产物只有渲染器自己的
               #   内部对账兜底 —— 那是**自证**，渲染器与被渲染物同错则同绿。
               #   T23 不 import 渲染器，独立解析 .md 并从冻结 CSV 重建每格。
    "T45",     # ★ 2026-10-07：**事后**分析 p45（分诊策略模拟）的独立验证（t45_triage.py）。
               #   与 T23 同构：不 import / 不 exec p45，自己从 p6b/p6c/p6d_outcome
               #   重建排名与命中数，并机械校验**预注册备忘的哈希未变**。
               #   其产物（Supplementary Figure S10 / Table S3）要进稿件 ⇒ 必须有独立判据。
    "T46",     # ★ 2026-10-07：补充表 S3 的独立验证（t46_tableS3.py）。
               #   p46 自带 8 条闭合断言，但那是**自证**（生成器与生成物同错则同绿）。
               #   t46 不 import p46：自己从 p6b/p6c/p6d 重建全部绩效格，
               #   并独立解析渲染出的 .md 逐格比对。
}
REQUIRED_FIELDS = ["path_id", "file", "symbol_or_region", "outputs_affected",
                   "risk_class", "validation_mode", "required_tests", "status"]

print("=" * 100, flush=True)
print("t14 — manifest 驱动的 registry（只读文本；不执行任何生产文件）", flush=True)
print("=" * 100, flush=True)
print(f"  PROD_DIR = {PROD}", flush=True)
print(f"  MANIFEST = {MANIFEST}", flush=True)
print(f"  模式     = {'PATCHED' if L.IS_PATCHED else 'REAL'}", flush=True)
print("", flush=True)


# ------------------------------------------------------------------ 基础工具（全部只读文本）
def read_text(path):
    """只读文本。用 errors='replace' —— 非 UTF-8 源文件**不得**让 census 崩溃
    （v1.3 已实测：GBK cookie 文件会让 encoding='utf-8' 抛 UnicodeDecodeError）。"""
    with open(path, "rb") as f:
        return f.read().decode("utf-8", errors="replace")


def root_of(anchor, tree_root):
    """锚点的解析根：`root: proj` -> 项目根；否则 -> **被扫的树**（PATCHED 时不是真实 PROD）。"""
    return PROJ if str(anchor.get("root", "prod")).lower() == "proj" else tree_root


def resolve_anchor(anchor, tree_root):
    """按文本解析单个锚点。返回 (ok, observed)。**不执行任何东西。**"""
    f = anchor.get("file")
    kind = anchor.get("kind")
    val = str(anchor.get("value", ""))
    contains = anchor.get("contains")
    base = root_of(anchor, tree_root)

    if kind == "glob":
        hits = glob.glob(os.path.join(base, f, val))
        return (len(hits) > 0), f"glob {f}/{val} -> {len(hits)} 命中"
    p = os.path.join(base, f)
    if not os.path.exists(p):
        return False, f"文件不存在 {f}"
    if kind == "regex":
        txt = read_text(p)
        n = len(re.findall(val, txt))
        return (n > 0), f"regex {val!r} -> {n} 命中"
    lines = read_text(p).split("\n")
    if kind == "def":
        ok = any(re.match(r"^def\s+" + re.escape(val) + r"\s*\(", ln) for ln in lines)
        return ok, f"顶层 def {val} -> {'存在' if ok else '缺失'}"
    if kind == "lines":
        try:
            if "-" in val:
                a, b = (int(x) for x in val.split("-", 1))
            else:
                a = b = int(val)
        except ValueError:
            return False, f"lines 值非法 {val!r}"
        inb = 1 <= a <= b <= len(lines)
        if not inb:
            return False, f"行号越界 {val}（文件 {len(lines)} 行）"
        if contains:
            seg = "\n".join(lines[a - 1:b])
            hit = re.search(contains, seg) is not None
            return hit, f"行段 {val} 在界内；contains={hit}"
        return True, f"行段 {val} 在界内（{len(lines)} 行）"
    return False, f"未知 anchor kind {kind!r}"


def scan(tree_root, entries):
    """对给定树做**存在性 + 锚点**扫描。返回 [(path_id, ok, observed)]。

    ★ 本函数**只调用 read_text / glob**，没有任何执行路径。
    下面的 MC-04 会用一个带模块级副作用的 canary 树跑**同一个函数**，
    证明「扫描不会执行被扫文件」。
    """
    out = []
    for e in entries:
        pid = e["path_id"]
        fspec = str(e.get("file", ""))
        if fspec.startswith("("):                       # (pending)
            out.append((pid, True, "deferred: file=(pending)"))
            continue
        files = [s.strip() for s in fspec.split(";") if s.strip()]
        miss = []
        for f in files:
            base = PROJ if str(e.get("root", "prod")).lower() == "proj" else tree_root
            pat = os.path.join(base, f)
            if any(c in f for c in "*?["):          # glob 形式的 file 字段
                if not glob.glob(pat):
                    miss.append(f"(glob 无命中) {f}")
            elif not os.path.exists(pat):
                miss.append(f)
        if miss:
            out.append((pid, False, f"缺文件 {miss}"))
            continue
        bad = []
        for a in (e.get("anchors") or []):
            ok, obs = resolve_anchor(a, tree_root)
            if not ok:
                bad.append(f"{a.get('file')}:{a.get('kind')}:{a.get('value')} -> {obs}")
        if bad:
            out.append((pid, False, "; ".join(bad[:3])))
        else:
            out.append((pid, True, f"{len(files)} 文件 + {len(e.get('anchors') or [])} 锚点 全部命中"))
    return out


# ------------------------------------------------------------------ 载入 manifest
_mtxt = read_text(MANIFEST) if os.path.exists(MANIFEST) else ""
_MSHA = hashlib.sha256(_mtxt.encode("utf-8")).hexdigest() if _mtxt else ""
_doc = None
_load_err = ""
try:
    _doc = L.load_manifest(MANIFEST)   # ★ 共享受保护 loader：重复键硬失败
except Exception as exc:                                 # noqa: BLE001
    _load_err = f"{type(exc).__name__}: {exc}"

if not isinstance(_doc, dict) or not isinstance(_doc.get("paths"), list):
    L.emit("MC-00", "FAIL", "manifest 必须可解析且含 paths 列表",
           f"err={_load_err or 'no paths'}", "t14 无法继续，后续判据全部不可引用")
    sys.exit(0)

PATHS = _doc["paths"]
META = _doc.get("meta") or {}
ENUM_RISK = set(META.get("risk_classes") or [])
ENUM_MODE = set((META.get("validation_modes") or {}).keys())
ENUM_STAT = set(META.get("statuses") or [])


# ------------------------------------------------------------------ MC-00 schema
def _mc00():
    bad, ids = [], []
    for i, e in enumerate(PATHS):
        tag = e.get("path_id") or f"<#{i}>"
        ids.append(tag)
        for k in REQUIRED_FIELDS:
            if k not in e:
                bad.append(f"{tag}: 缺字段 {k}")
        if e.get("risk_class") not in ENUM_RISK:
            bad.append(f"{tag}: risk_class={e.get('risk_class')!r} 不在枚举")
        if e.get("validation_mode") not in ENUM_MODE:
            bad.append(f"{tag}: validation_mode={e.get('validation_mode')!r} 不在枚举")
        if e.get("status") not in ENUM_STAT:
            bad.append(f"{tag}: status={e.get('status')!r} 不在枚举")
        rt = e.get("required_tests")
        if not isinstance(rt, list):
            bad.append(f"{tag}: required_tests 必须是列表")
            rt = []
        # ★ 2026-10-04 事故（审计器自身的错，closure 首轮抓到）：
        #   我把 BW-07 / EC-14 / UNIT-02,03,05 / L3-03b,07 / CONC-05,05C / SUDO-04
        #   **凭印象**填进了 required_tests，而它们在代码里自报 warn_only / characterisation。
        #   fail-closed 规则下，一个只发 WARN 的 required_test 会让该 path 永判 UNKNOWN。
        #   下面两条断言把这类错误变成**结构性不可能**：两者必须互斥、且都必须是字符串列表。
        ch = e.get("characterisations") or []
        if not isinstance(ch, list):
            bad.append(f"{tag}: characterisations 必须是列表")
            ch = []
        if set(rt) & set(ch):
            bad.append(f"{tag}: required_tests 与 characterisations 重叠 {sorted(set(rt) & set(ch))}"
                       "（一个测试不能既是 gate 又是 characterisation）")
        if e.get("risk_class") == "Critical":
            if not rt and not e.get("deferred_reason"):
                bad.append(f"{tag}: Critical 必须有 required_tests 或 deferred_reason")
            fi = e.get("fault_injection")
            if not isinstance(fi, dict) or not all(k in fi for k in ("id", "kind", "mutation", "expect")):
                bad.append(f"{tag}: Critical 必须有完整 fault_injection(id/kind/mutation/expect)")
            if fi and str(fi.get("expect", "")).upper() != "FAIL":
                bad.append(f"{tag}: fault_injection.expect 必须是 FAIL")
    # ★ 三阶段 closure 的 schema（2026-10-05）
    _PH = set(META.get("validation_phases") or [])
    _ET = set(META.get("evidence_type_enum") or [])
    _declared = _doc.get("evidence_types") or {}
    for e in PATHS:
        tag = e.get("path_id")
        if e.get("phase") is not None and e["phase"] not in _PH:
            bad.append(f"{tag}: phase={e['phase']!r} 不在 {sorted(_PH)}")
        for g in (e.get("deferred_gates") or []):
            if not isinstance(g, dict) or not all(k in g for k in ("test", "phase", "reason", "reentry")):
                bad.append(f"{tag}: deferred_gate 缺 test/phase/reason/reentry :: {g}")
            elif g["phase"] not in _PH:
                bad.append(f"{tag}: deferred_gate phase={g['phase']!r} 不在枚举")
            elif g["phase"] == (e.get("phase") or "pre_correction_freeze"):
                bad.append(f"{tag}: deferred_gate {g['test']} 的 phase 与本 path 相同 —— "
                           f"同阶段的东西不能叫 deferred（用 required_tests 或 pending_gates）")
        for g in (e.get("pending_gates") or []):
            if not isinstance(g, dict) or not all(k in g for k in ("test", "phase", "reason")):
                bad.append(f"{tag}: pending_gate 缺 test/phase/reason :: {g}")
        # 每条 required_test 必须有 evidence_type（族前缀或精确 test_id）
        for x in (e.get("required_tests") or []):
            if _evidence_of(x, _declared) is None:
                bad.append(f"{tag}: required_test {x} 在 evidence_types 中无声明")
    for k, v in _declared.items():
        for x in (v or []):
            if x not in _ET:
                bad.append(f"evidence_types[{k}]: {x!r} 不在枚举")
    dup = {x for x in ids if ids.count(x) > 1}
    if dup:
        bad.append(f"path_id 重复: {sorted(dup)}")
    return bad, ids


def _evidence_of(test_id, declared):
    """解析一个 test_id 的证据类型：精确 test_id 优先，其次最长匹配的族前缀。"""
    if test_id in declared:
        return set(declared[test_id])
    best = None
    for k in declared:
        if k.endswith("-") and test_id.startswith(k):
            if best is None or len(k) > len(best):
                best = k
    return set(declared[best]) if best else None


# validation_mode -> 必须具备的证据类型（Rendering 不强制）
_MODE_NEEDS = {
    "production-exec":    {"production_exec"},
    "independent-reimpl": {"independent_reference"},
    "both":               {"production_exec", "independent_reference"},
    "output-level-MNV":   {"artifact_validation"},
    "hash-render-check":  set(),
}


_bad00, PIDS = _mc00()
L.check("MC-00", not _bad00,
        "manifest schema 合法（8 必备字段 + 枚举 + Critical 的 required_tests/deferred_reason + fault_injection）",
        f"{len(PATHS)} 条 path；违规 {len(_bad00)} 项" + ("" if not _bad00 else f" :: {' | '.join(_bad00[:6])}"),
        "")
if _bad00:
    for b in _bad00:
        print(f"      !! {b}", flush=True)


# ------------------------------------------------------------------ MC-01 / MC-02
_res = scan(PROD, PATHS)
_bad01 = [(p, o) for p, ok, o in _res if not ok and "缺文件" in o]
_bad02 = [(p, o) for p, ok, o in _res if not ok and "缺文件" not in o]

L.check("MC-01", not _bad01,
        "每条 path 的 file 在登记根目录下存在",
        f"{len(PATHS) - len(_bad01)}/{len(PATHS)} 条文件齐全"
        + ("" if not _bad01 else f" :: {' | '.join(f'{p}({o})' for p, o in _bad01[:6])}"), "")

L.check("MC-02", not _bad02,
        "每条 path 的 anchors 可按**纯文本**解析（def/lines/regex/glob；不执行）",
        f"{len(PATHS) - len(_bad02)}/{len(PATHS)} 条锚点全中"
        + ("" if not _bad02 else f" :: {' | '.join(f'{p}({o})' for p, o in _bad02[:6])}"), "")
for p, o in (_bad01 + _bad02)[:12]:
    print(f"      !! {p} :: {o}", flush=True)


# ------------------------------------------------------------------ MC-03 结构断言
def _mc03():
    """本模块源码中不得存在动态执行原语，且不得调用装载器。

    ★ 为什么用 AST 而不是 grep：注释/docstring 里出现 "exec(" 不应误判。
    """
    src = read_text(SELF)
    tree = ast.parse(src, filename=SELF)
    viol = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Name):
            if node.func.id in ("exec", "eval", "compile", "__import__"):
                viol.append(f"{node.func.id}() @L{node.lineno}")
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute):
            if node.func.attr in ("load_def", "exec_block", "grab"):
                viol.append(f"L.{node.func.attr}() @L{node.lineno}")
        if isinstance(node, ast.Import):
            for a in node.names:
                if a.name.split(".")[0] in ("importlib", "runpy"):
                    viol.append(f"import {a.name} @L{node.lineno}")
        if isinstance(node, ast.ImportFrom):
            if (node.module or "").split(".")[0] in ("importlib", "runpy"):
                viol.append(f"from {node.module} import ... @L{node.lineno}")
        if isinstance(node, ast.Attribute) and node.attr in ("load_def", "exec_block"):
            viol.append(f"属性访问 .{node.attr} @L{node.lineno}")
    return viol


_v3 = _mc03()
L.check("MC-03", not _v3,
        "registry 模块源码中不存在 exec/eval/compile/__import__/importlib/runpy，且不调用装载器",
        "AST 扫描 0 违规" if not _v3 else f"{len(_v3)} 违规 :: {' | '.join(_v3[:6])}",
        "结构断言：即使将来有人把 exec 加回来，也会立刻在此判红")


# ------------------------------------------------------------------ MC-04 canary（永久 negative test）
scratch = tempfile.mkdtemp(prefix="evs_registry_canary_")
_tree = os.path.join(scratch, "tree")
os.makedirs(_tree, exist_ok=True)
_canary_out = os.path.join(scratch, "CANARY_CREATED.txt")
_np_out = os.path.join(scratch, "np_canary.npy")
_victim = os.path.join(scratch, "victim.txt")
with io.open(_victim, "w", encoding="utf-8") as _f:
    _f.write("keep me")

# canary fixture：**模块级**副作用（写入 / 删除 / np.save / 污染共享 np 对象）
# 同时长得像个"族内副本"（有这个家族的名字），保证扫描会碰到它。
_canary_py = os.path.join(_tree, "canary_probe.py")
with io.open(_canary_py, "w", encoding="utf-8") as _f:
    _f.write(
        "# -*- coding: utf-8 -*-\n"
        "import os, numpy as np\n"
        f"_C = open(r'{_canary_out}', 'w'); _C.write('EXECUTED'); _C.close()\n"
        f"_R = os.remove(r'{_victim}')\n"
        f"_N = np.save(r'{_np_out}', np.zeros(3))\n"
        "setattr(np, 'EVS_POISON_FLAG', 12345)\n"
        "def psijs(x):\n"
        "    return x\n"
    )

_canary_entry = [{
    "path_id": "CANARY", "file": "canary_probe.py",
    "anchors": [{"file": "canary_probe.py", "kind": "def", "value": "psijs"}],
}]


def _side_effects_present():
    return (os.path.exists(_canary_out), not os.path.exists(_victim),
            os.path.exists(_np_out), hasattr(np, "EVS_POISON_FLAG"))


def _reset_fixture():
    for p in (_canary_out, _np_out):
        if os.path.exists(p):
            os.remove(p)
    if not os.path.exists(_victim):
        with io.open(_victim, "w", encoding="utf-8") as f:
            f.write("keep me")
    if hasattr(np, "EVS_POISON_FLAG"):
        delattr(np, "EVS_POISON_FLAG")


def _mc04a():
    """正对照：**在子进程里直接运行 canary 脚本**（不是 exec，也不碰生产代码），
    证明这组副作用确实可被检出 —— 否则负对照是空转。"""
    _reset_fixture()
    r = subprocess.run([sys.executable, _canary_py], capture_output=True, text=True,
                       encoding="utf-8", errors="replace", cwd=scratch,
                       env={**os.environ, "PYTHONPATH": ""})
    created, removed, npsaved, poisoned = _side_effects_present()
    # 子进程污染不到父进程的 np；父进程 np 必须仍然干净
    return (created and removed and npsaved and not poisoned), \
        f"rc={r.returncode} created={created} victim_removed={removed} np_saved={npsaved} parent_np_clean={not poisoned}"


_ok4a, _obs4a = _mc04a()
L.check("MC-04a", _ok4a,
        "正对照：canary 的模块级副作用确实可被检出（判据非空转）",
        _obs4a, "子进程直接运行 canary 脚本；t14 自身不含 exec")
_reset_fixture()

# ---- 负对照：用**同一个 scan()** 扫 canary 树
_scanres = scan(_tree, _canary_entry)
_found = bool(_scanres) and _scanres[0][1]
_created, _removed, _npsaved, _poisoned = _side_effects_present()
_clean = (not _created) and (not _removed) and (not _npsaved) and (not _poisoned)

L.check("MC-04b", _clean,
        "registry 扫描 canary 树后 **副作用没有发生**（不执行被扫文件的模块级代码）",
        f"created={_created} victim_removed={_removed} np_saved={_npsaved} np_poisoned={_poisoned}",
        "★ 永久 negative test：若将来有人把 exec/import 加回 registry，此条立刻判红")
L.check("MC-04c", _found,
        "扫描非空转：canary 文件确实被扫到并解析了锚点",
        f"scan -> {_scanres}", "防止把 registry 写成什么都不做的空壳")
if not (_clean and _found):
    print(f"      !! canary 树保留在 {scratch} 以便复查（未清理）", flush=True)
else:
    shutil.rmtree(scratch, ignore_errors=True)


# ------------------------------------------------------------------ MC-05 manifest 身份
def _mc05():
    problems = []
    if not _mtxt:
        problems.append("manifest 不存在或为空")
    try:
        if os.path.abspath(MANIFEST).startswith(PROD + os.sep):
            problems.append("manifest 位于 PROD_DIR 内（违反「registry 不由被测代码生成」）")
    except Exception:                                    # noqa: BLE001
        pass
    au = str(META.get("authored_by", ""))
    if "hand" not in au.lower():
        problems.append(f"meta.authored_by 未声明手工撰写: {au!r}")
    if str(META.get("prod_dir", "")).strip() and \
            os.path.abspath(str(META["prod_dir"])) != PROD and not L.IS_PATCHED:
        problems.append(f"meta.prod_dir({META['prod_dir']}) 与实际 PROD_DIR({PROD}) 不一致")
    return problems


_p5 = _mc05()
L.check("MC-05", not _p5,
        "manifest 身份：可读、位于 PROD_DIR **之外**、声明手工撰写",
        f"sha256={_MSHA[:16]}…  size={len(_mtxt)}B"
        + ("" if not _p5 else f" :: {' | '.join(_p5)}"), "")
if _MSHA:
    # ★ 2026-10-07 修正：原先写 len(_mtxt)（**字符**数）却标 "bytes" —— manifest 含中文注释，
    #   两者相差 15039（50655 字符 vs 65694 字节），是个会误导人的数字。改为真字节数。
    _nbytes = os.path.getsize(MANIFEST) if os.path.exists(MANIFEST) else 0
    with io.open(os.path.join(HERE, "manifest_sha256.txt"), "w", encoding="utf-8") as _f:
        _f.write(f"{_MSHA}  ACTIVE_PATH_MANIFEST.yaml  ({_nbytes} bytes)\n"
                 f"# 本文件由 t14 写出，用于记录 closure 时点的 manifest 身份\n")


# ------------------------------------------------------------------ MC-06 Critical 声明完整性
_crit = [e for e in PATHS if e.get("risk_class") == "Critical"]
_deferred = [e["path_id"] for e in _crit if not e.get("required_tests")]
_sec = [e for e in PATHS if e.get("risk_class") == "Secondary"]
_ren = [e for e in PATHS if e.get("risk_class") == "Rendering"]
_n_tests = sum(len(e.get("required_tests") or []) for e in _crit)
L.check("MC-06", True,
        "Critical / Secondary / Rendering 分类计数与 deferred 项显式列出",
        f"Critical {len(_crit)} / Secondary {len(_sec)} / Rendering {len(_ren)}；"
        f"Critical 引用测试 {_n_tests} 条；deferred(将计为 UNKNOWN) = {_deferred or '无'}",
        "deferred 项**必须**在 closure 中呈现为 UNKNOWN，不得静默略过")


# ------------------------------------------------------------------ MC-07 test_id 族前缀
def _mc07():
    bad = []
    for e in PATHS:
        for t in (e.get("required_tests") or []):
            fam = re.split(r"[-:]", str(t), maxsplit=1)[0].strip()
            if fam not in KNOWN_FAMILIES:
                bad.append(f"{e['path_id']} -> {t}")
    return bad


_b7 = _mc07()
L.check("MC-07", not _b7,
        "required_tests 的族前缀均为已知族（抓 test_id 拼写错误）",
        f"已知族 {sorted(KNOWN_FAMILIES)}；未识别 {len(_b7)} 条"
        + ("" if not _b7 else f" :: {' | '.join(_b7[:6])}"), "")


# ------------------------------------------------------------------ MC-08 绑定依据
# ★ 2026-10-04 事故：`required_tests` 曾**凭印象**填写，把 10 条自报 warn_only / characterisation
#   的测试当成 gate 收进登记表。fail-closed 规则下，一个只发 WARN 的 required_test 会让该 path
#   永判 UNKNOWN —— 这是登记错误，不是覆盖漏洞。
#
# ★ 本判据把「依据」变成机械可验的：**一个 test_id 是 gate 还是 characterisation，由其所在模块
#   自己的发射语句决定**，不采信任何人的记忆。只看**实际发射点**（emit/warn 调用的实参位置），
#   不看 docstring、不看 expected 文本里顺带提到的 tid。
#
#   模块源码里对 tid 的发射形态：
#     L.emit("TID", "WARN", ...) / emit("TID","WARN",...)      -> characterisation
#     warn("TID", ...) / _run_warn("TID", ...)                 -> characterisation
#     <call>(..., "TID", ..., warn_only=True)                  -> characterisation
#     L.check("TID", ...) / L.guard("TID", ...)                -> gate
#   规则：manifest.required_tests 里的每一条都**不得**是 characterisation；
#         manifest.characterisations 里的每一条都**必须**是 characterisation。
def _emission_index():
    """扫验证套件自身的源码，返回 {tid: [(module, line, mode)]}，mode ∈ {gate, char}。"""
    idx = {}
    mods = [f for f in sorted(os.listdir(HERE)) if re.match(r"t\d\d_.*\.py$", f)]
    for m in mods:
        try:
            lines = read_text(os.path.join(HERE, m)).split("\n")
        except OSError:
            continue
        for ln, text in enumerate(lines, 1):
            for mm in re.finditer(r'"(LOAD|MMD|CONC|L3|SUDO|UNIT|BOUND|FORK|PSI|BW|SENS|EVAL|EC)'
                                  r'-[0-9A-Za-z_.:-]+"', text):
                tid = mm.group(0)[1:-1]
                before = text[:mm.start()]
                after = text[mm.end():]
                # 必须处在一次调用的**实参位置**：前缀里有未闭合的 `(`，且该 `(` 后紧跟的
                # callee 是 emit/warn/check/guard/run 一族
                call = re.search(r"([A-Za-z_][A-Za-z0-9_.]*)\s*\(([^()]*)$", before)
                if not call:
                    continue
                callee, first_arg = call.group(1), call.group(2)
                if not re.fullmatch(r"\s*", first_arg):     # tid 必须是**第一个**实参
                    continue
                # ★ 判据只建立在**能确证**的证据上，不做启发式推断：
                #   · 状态实参是**字面量 "WARN"**        -> 结构性 WARN（永不可能 PASS）=> char
                #   · callee 名字里含 warn               -> 同上 => char
                #   · 其余（含 PASS/FAIL 条件式、动态 tid、混合 helper）-> **不判定**，
                #     只登记发射点，交由人/GPT6 复核。**不猜。**
                lit_warn = bool(re.match(r'\s*,\s*"WARN"', after))
                is_char = lit_warn or ("warn" in callee.lower())
                idx.setdefault(tid, []).append(
                    (m, ln, f"{callee}(...)", "char" if is_char else "unjudged"))
    return idx


_EI = _emission_index()


def _is_char(tid):
    recs = _EI.get(tid) or []
    if not recs:
        return None, []
    chars = [r for r in recs if r[3] == "char"]
    others = [r for r in recs if r[3] == "unjudged"]
    # 只要**存在**一个「字面 WARN / warn-helper」发射点且**没有任何**非 WARN 发射点，
    # 才能确证它是 characterisation（=它永远不可能 PASS）。
    if chars and not others:
        return True, recs
    return False, recs


_b8 = []          # 确证违规 -> FAIL
_b8_ev = []       # 证据缺失 -> 呈现，不判红
for e in PATHS:
    for t in (e.get("required_tests") or []):
        is_char, recs = _is_char(t)
        if is_char is None:
            _b8_ev.append(f"{e['path_id']}::{t} 无字面发射点（可能为动态拼 tid，"
                          f"如 t06 的 \"BOUND-04::%s\" % fname）")
        elif is_char:
            loc = "; ".join(f"{m}:{ln} {c}" for m, ln, c, _ in recs)
            _b8.append(f"{e['path_id']}::required {t} 是**结构性 WARN**，永远不可能 PASS @ {loc}")

_n_req = sum(len(e.get("required_tests") or []) for e in PATHS)
_n_char = len([t for e in PATHS for t in (e.get("characterisations") or [])])
L.check("MC-08", not _b8,
        "登记为 required_test 的测试，其**发射语句不得是结构性 WARN**"
        "（字面 \"WARN\" 状态或 warn-helper）—— 依据是发射点本身，不采信任何人的记忆",
        f"发射点索引 {len(_EI)} 个 test_id；required {_n_req} 条 / characterisation {_n_char} 条；"
        f"确证违规 {len(_b8)}；无法从字面确证 {len(_b8_ev)}"
        + ("" if not _b8 else f" :: {' | '.join(_b8[:6])}"), "")
for _x in _b8[:12]:
    print(f"      !! {_x}", flush=True)
if _b8_ev:
    print(f"      （无法从字面确证 {len(_b8_ev)} 条，逐条列于 coverage_table 的 required_tests 列，"
          f"须人工/外部复核，**不自动判定**）", flush=True)
    with io.open(os.path.join(HERE, "mc08_evidence_missing.txt"), "w", encoding="utf-8") as _f:
        _f.write("\n".join(_b8_ev) + "\n")



# ------------------------------------------------------------------ MC-09 validation_mode 的名实相符
# ★ 2026-10-05：`validation_mode` 是**人工声明**。CONFIG_P12 声明 independent-reimpl，绑定却只有
#   两条 BOUND-04 冻结产物扫描 —— 声明 ≠ 证据。本判据把声明变成**可证实的证据类型**：
#   每条 required_test 有 evidence_type（族级声明，精确 test_id 优先），
#   而 path 的 validation_mode 所要求的最小证据类型集合必须是其 required_tests 证据类型的子集。
#   ★ 判据只在 required_tests 非空时生效；空集合的 path 已由 deferred_reason/pending_gates 计
#     UNKNOWN，不在此重复扣分。
_DECL = _doc.get("evidence_types") or {}
_b9 = []
for _e in PATHS:
    _req = list(_e.get("required_tests") or [])
    if not _req:
        continue
    _need = _MODE_NEEDS.get(_e.get("validation_mode"))
    if _need is None:
        _b9.append(f"{_e['path_id']}: 未知 validation_mode {_e.get('validation_mode')!r}")
        continue
    _have = set()
    for _x in _req:
        _ev = _evidence_of(_x, _DECL)
        if _ev:
            _have |= _ev
    _miss = _need - _have
    if _miss:
        _b9.append(f"{_e['path_id']}: mode={_e.get('validation_mode')} 需要 {sorted(_miss)}，"
                   f"但其 required_tests 只提供了 {sorted(_have)}")

L.check("MC-09", not _b9,
        "validation_mode 名实相符：声明的证据类型必须由其 required_tests 的 evidence_type 实际提供"
        "（依据 = evidence_types 声明，不采信 validation_mode 的字面）",
        f"证据族 {len(_DECL)} 个；检查 {len([e for e in PATHS if e.get('required_tests')])} 条有绑定 path；"
        f"违规 {len(_b9)}" + ("" if not _b9 else f" :: {' | '.join(_b9[:6])}"), "")
for _x in _b9[:12]:
    print(f"      !! {_x}", flush=True)


# ------------------------------------------------------------------ MC-10 每个 FAIL 都要有归属
# ★ 2026-10-05：closure 四条条件只对 Critical 的 required_tests 判定，因此**一个红着的测试
#   可以完全不被任何判据看见**。PATCHED 首次全绿时，t04/t06/t07/t12/t13 仍有 17 条 FAIL ——
#   逐条核对后发现 10 条 BOUND-01::* 与 EVAL-B1 **在登记表里没有归属**（我临时写脚本才看出来）。
#   本判据把这件事机械化：被测树上**任何 FAIL 的 test_id** 必须能归入下述四类之一，否则判红。
_INFRA = ("LOAD-", "MC-")      # 与 run_evs.py 同一口径：基础设施自测不属于任何生产 path


def _mc10():
    import csv as _csv
    rows = []
    for mode in ("real", "patched"):
        q = os.path.join(HERE, f"evs_results_{mode}.csv")
        if not os.path.exists(q):
            continue
        with io.open(q, encoding="utf-8") as fh:
            for r in _csv.DictReader(fh):
                if r.get("status") == "FAIL":
                    rows.append((mode, r["test_id"]))
    req, ch, dg, ub = set(), set(), set(), set()
    for e in PATHS:
        req |= set(e.get("required_tests") or [])
        ch |= set(e.get("characterisations") or [])
        dg |= {g["test"] for g in (e.get("deferred_gates") or [])}
    for g in (_doc.get("unbound_deferred") or []):
        ub.add(g["test"])

    def _known(x):
        # ★ 2026-10-06：基础设施自测（装载器 LOAD-*、registry MC-*）不属于任何生产 path。
        #   run_evs.py 早就这么算了（`_INFRA = ("LOAD-", "MC-")`），但本判据漏了这一步，
        #   于是 MC-10 一旦判红，它读到的 CSV 里就带着 `MC-10 FAIL`，而 MC-10 又不在登记表里
        #   ⇒ 下一轮它把自己算成「未归类」，**永久自锁**（2026-10-06 的闭包重跑实测到这一点）。
        #   豁免后与 run_evs 口径一致，且消除自指。
        if x.startswith(_INFRA):
            return True
        if x in req or x in ch or x in dg or x in ub:
            return True
        for pat in list(dg) + list(ub):        # 前缀登记（如 "BOUND-01::*"）
            if pat.endswith("*") and x.startswith(pat[:-1]):
                return True
        return False

    bad = sorted({f"{m}:{x}" for m, x in rows if not _known(x)})
    return bad, sorted({x for _, x in rows})


_b10, _allfail = _mc10()
L.check("MC-10", not _b10,
        "被测树上每一条 FAIL 的 test_id 都必须有归属：required_tests / characterisations / "
        "deferred_gates / unbound_deferred（含前缀登记）—— 不得有来历不明的红",
        f"两个模式下 FAIL 的 test_id 共 {len(_allfail)} 条；未归类 {len(_b10)}"
        + ("" if not _b10 else f" :: {_b10[:8]}"), "")
for _x in _b10[:12]:
    print(f"      !! 未归类：{_x}", flush=True)

print("", flush=True)
print(f"t14 done  (sha256 {_MSHA[:16]}…)", flush=True)