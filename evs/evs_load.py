# -*- coding: utf-8 -*-
r"""evs_load — estimator validation suite 的**逐字生产代码装载器**

================================ 为什么必须有它 ================================
GPT6 的硬性要求：

    "优化实现与论文公式的 brute-force 实现逐项对拍，**而不是代码副本对拍**"

所以要测的必须是**生产代码本身**，不是我把函数复制一份来测。本模块按**名称**
（顶层 def）或**行号区间**（内联块）把生产代码的**原始文本**抽出来并 exec，
保证测试跑的是真代码。

================================ 输出契约（强制） ================================
每个测试模块 print 形如（`|` 分隔，6 段）：

    RESULT|<test_id>|<PASS|FAIL|WARN>|<expected>|<observed>|<note>

- `test_id`  : 稳定标识，如 MMD-01
- `expected` : **期望来自哪里**（数学恒等式 / 定义 / 论文公式 / 代码规格），不是"应该通过"
- `observed` : 实测值
- `note`     : 补充

**禁止只打印 PASS 而不打印 expected/observed** —— 那是"审计器不被审计"。
`WARN` 用于 diagnostic bound（如 DR 越界），`FAIL` 用于 hard bound 或恒等式违反。

用法（测试模块内）：
    import evs_load as L
    fn = L.load_def(L.P6B, "mmd_rbf")
    L.emit("MMD-01", "PASS" if ok else "FAIL", "identical samples -> 0", f"{v:.6f}", "")
"""
import os, re, sys, textwrap, traceback
import numpy as np
import yaml

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

HERE = os.path.dirname(os.path.abspath(__file__))          # 第二篇 新方案\evs
OUT = os.path.join(os.path.dirname(HERE), "分析结果")        # 生产树
PROJ = os.path.dirname(HERE)                               # 第二篇 新方案
MANIFEST = os.path.join(PROJ, "ACTIVE_PATH_MANIFEST.yaml")


# ---------------------------------------------------------------- 受保护的 YAML 读取
# ★ 2026-10-05 事故：manifest 里 EPISODE_UNIT 出现了**两次** `characterisations:`（行 122/128）。
#   PyYAML 静默取最后一个，UNIT-07 被无声丢弃；而 MC-00 的互斥检查跑在**存活列表**上，通过了。
#   silent overwrite 会**直接篡改 closure 输入**，因此这里硬失败。
#   本函数是 manifest 的**唯一**入口：t14 与 run_evs 都必须用它，不得各自 yaml.safe_load。
class _NoDupLoader(yaml.SafeLoader):
    """SafeLoader + 重复键硬失败。"""


def _no_dup_mapping(loader, node, deep=False):
    mapping, seen = {}, {}
    for k_node, v_node in node.value:
        k = loader.construct_object(k_node, deep=deep)
        if k in seen:
            raise ValueError(
                f"duplicate key {k!r} —— 首次出现 line {seen[k] + 1}，"
                f"再次出现 line {k_node.start_mark.line + 1}。"
                f"重复键会静默覆盖前值，必须合并后重跑。")
        seen[k] = k_node.start_mark.line
        mapping[k] = loader.construct_object(v_node, deep=deep)
    return mapping


_NoDupLoader.add_constructor(
    yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG, _no_dup_mapping)


def load_manifest(path=None):
    """读取 manifest；重复键 -> 抛错（不静默覆盖）。"""
    with open(path or MANIFEST, encoding="utf-8") as f:
        return yaml.load(f, Loader=_NoDupLoader)

# ---------------------------------------------------------------- 生产代码根目录
# 默认 = 真实生产目录。设 EVS_PROD_DIR 可指向**打过补丁的副本**，
# 从而在**同一套测试**上验证"修好之后能否全绿"。
# 这是判定本套件是否为"下一代自证系统"的唯一手段：
#   红灯基线（真实生产） + 绿灯验证（补丁副本） 必须同时成立。
PROD_DIR = os.environ.get("EVS_PROD_DIR") or OUT
IS_PATCHED = os.path.abspath(PROD_DIR) != os.path.abspath(OUT)

P6B = os.path.join(PROD_DIR, "p6b_01_corrected_pipeline.py")
P6C = os.path.join(PROD_DIR, "p6c_01_L3_blinded.py")
P10A = os.path.join(PROD_DIR, "p10a_prepare_cache.py")
P10B = os.path.join(PROD_DIR, "p10b_subsample_curves.py")
P5 = os.path.join(PROD_DIR, "p5_01_degradation_benchmark.py")
P6D = os.path.join(PROD_DIR, "p6d_01_unseal_evaluation.py")
P21A = os.path.join(PROD_DIR, "p21a_reporting_profile_derivation.py")
EPI = os.path.join(PROD_DIR, "episode_definition.py")
P6_02 = os.path.join(PROD_DIR, "p6_02_label_free_signals.py")
P14 = os.path.join(PROD_DIR, "p14_rebuild_fig3b.py")
P6FD = os.path.join(PROD_DIR, "p6f_lodo_settlement.py")

# 生产代码里内联块的精确行号（1-indexed，含两端）
BLOCK_L3B = (P6C, 257, 289)      # L3-B: transported AUC（IOW / OM / DR）
BLOCK_SUDO = (P6C, 291, 330)     # L3-A: SUDO-AURCC
BLOCK_STORE = (P6B, 121, 144)    # unit / epid / wide / meta / X 构造

# ---- fork-parity：活的同源副本（同样的估计量在别的文件里被重写了一遍）----
BLOCK_P10B_SUDO = (P10B, 137, 169)   # p10b 的 SUDO-AURCC 内联块
BLOCK_P10B_L3 = (P10B, 171, 195)     # p10b 的 transported-AUC 内联块（IOW/OM/DR）
BLOCK_P6_02_MMD = (os.path.join(PROD_DIR, "p6_02_label_free_signals.py"), 223, 237)
BLOCK_P10A_MMD = (P10A, 291, 298)    # p10a 的 MMD source 内联块

# ---- v1.1：其余活估计量核（critic 扫出、此前无任何测试执行）----
BLOCK_BW_P6B = (P6B, 267, 272)       # 带宽中位数启发式（决定所有 MMD 的尺度）
BLOCK_BW_P10A = (P10A, 242, 246)     # 同上的 p10a 副本
BLOCK_P10B_L2 = (P10B, 128, 135)     # p10b 的 prediction-space L2 统计
BLOCK_P10B_SENS = (P10B, 305, 312)   # p10b 的 sens_at_spec（**嵌套** def，非顶层）
P10C = os.path.join(PROD_DIR, "p10c_aggregate.py")

# ★★ 2026-10-06：**块首行指纹**（由当时**已人工核验**的状态生成），用于 t14 之外的
#    `evs_selftest` LOAD-05 错位检测。
#   为什么必须有它：t14 的 MC-02 只验 **manifest 的 anchors** 能按文本解析，
#   **不验** evs 代码里 exec_block 的行号元组 —— 块被指到同一文件里的**另一段代码**时，
#   MC-02 照样 PASS（只要那段也能 exec）。实测事故：给 p10a/p10b 插入代码后，
#   `BLOCK_P10B_L3` 变成了抽出 SUDO 块、`BLOCK_P10A_MMD` 变成了抽出"[校验闸门]通过"——
#   是**人工按内容对拍**才发现的，没有任何判据报警。
#   指纹 = sha256(`L.grab(*block)` 的首行去缩进文本)。改行号时若不同步改指纹，LOAD-05 立刻判红。
BLOCK_FINGERPRINTS = {
    "BLOCK_BW_P10A": "1db19b4b78aa9b3181311e1311b31fe9fabba986f589d8b8bfe99599a10aaff0",
    "BLOCK_BW_P6B": "2ddefc74316b78f1c05366ac48bcb47de1ecacdcfed1b9ed066e4aef3ca3ac91",
    "BLOCK_L3B": "b7b01a964af5299372065bf9b5566fbef4d79c1a7c1a35460f16aedc4189f7b1",
    "BLOCK_P10A_MMD": "bd8d65280998ba9da739a3571b582dc8f55d79f8f46d48c85d1c59e7aa51c349",
    "BLOCK_P10B_L2": "1fd431eac7d2275947b917407c4c1d28d40ee74460ece1b2426fd34cf0da50e8",
    "BLOCK_P10B_L3": "848ef7ea8f53c1d23843d0444fe07c1c24d6e815d2d6c3f5ee211b84683ba758",
    "BLOCK_P10B_SENS": "64727384df93da8ffbd5c6035c6865bb5dbbdd8b1824f29bb4b2abd9e04d96cd",
    "BLOCK_P10B_SUDO": "9c556096246daa990239fe2ac4cdc9a8c7744f58fc8bc48c1a6dad6d3ccf2b0e",
    "BLOCK_P6_02_MMD": "286d320efdaafac19e9dcf5e619c251773f1c99437e918b2fdea5b10a608f34f",
    "BLOCK_STORE": "1d5a612fc3980ba9c0ff584832b91c1e0eea60000f4c282c0b46ba26fc61da68",
    "BLOCK_SUDO": "326fd238b0dcf0b99880b249b7106c7de9622fb1f4cef043f9f80a2d0a034523",
}

# v1.1 目标：这些**顶层 def** 需逐个用 load_def 抽出并逐字执行
# 注意：p10b 的 sens_at_spec 是嵌套 def，用 BLOCK_P10B_SENS 而非 load_def。
V11_TOPLEVEL_DEFS = {
    P6B:  ["psijs", "feat_psijs", "js_cont", "ent"],
    P10B: ["psijs_from_hist"],
    P10C: ["sens_at_spec"],
    P10A: [], P6C: [], P6_02: [],
}

# 需要做「同源副本清点」的估计量核：名字 -> 已知可能存在的 (文件, 抽取方式)
FORK_CENSUS = {
    "mmd_rbf": [
        (P6B, "def:216"), (P10A, "block:273-280"), (P10B, "inline:85-93"),
        (os.path.join(PROD_DIR, "p6_02_label_free_signals.py"), "def:223"),
        (os.path.join(PROD_DIR, "p6_superseded_pre_correction", "p6_02_label_free_signals.py"), "def:223"),
        (os.path.join(PROD_DIR, "p27_dr_formula_diagnostic.py"), None),   # 诊断脚本，可能无核
    ],
    "wconcordance": [
        (P6C, "def:84"), (P10B, "def:40"),
        (os.path.join(PROD_DIR, "p27_dr_formula_diagnostic.py"), "def"),
    ],
    # ★ CS-09（2026-10-06）：w_ESS 是**内联表达式**，不是顶层 def，所以只能按 block 抽取。
    #   ★ 但必须写明：**本 dict 目前没有任何代码程序化消费它**（全项目只有 t07:1206 在
    #     字符串里提到 "evs_load.FORK_CENSUS checked, not trusted"）。所以往这里加一行是
    #     **登记**，不是**强制** —— 若只加这一行而不加判据，等于什么都没做。
    #   真正的强制在两处：
    #     · FORK-07 缺陷表 (g) 行（枚举"哪些副本必须一致"）
    #     · t03 的 ESS-01（p6c 侧）+ t07 的 ESS-02（p10b 侧）（行为对拍）
    "w_ESS": [
        (P6C, "block:257-289"), (P10B, "block:150-174"),
    ],
}


def read_lines(path):
    with open(path, encoding="utf-8") as f:
        return f.read().split("\n")


def grab(path, a, b):
    """取 path 的第 a..b 行（1-indexed，含两端），并去掉公共缩进。"""
    lines = read_lines(path)
    if a < 1 or b > len(lines) or a > b:
        raise ValueError(f"行号越界 {path}: {a}-{b}（共 {len(lines)} 行）")
    return textwrap.dedent("\n".join(lines[a - 1:b]))


def load_def(path, name, extra_ns=None):
    """按名称抽取**顶层 def** 的原始文本并 exec，返回函数对象。

    抽取规则：从 `^def name(` 起，直到第一个**顶格且非空**的行之前。
    """
    lines = read_lines(path)
    pat = re.compile(r"^def\s+" + re.escape(name) + r"\s*\(")
    start = None
    for i, ln in enumerate(lines):
        if pat.match(ln):
            start = i
            break
    if start is None:
        raise LookupError(f"{path} 中找不到顶层 def {name}")
    end = start + 1
    while end < len(lines):
        ln = lines[end]
        if ln.strip() and not ln[:1].isspace():
            break
        end += 1
    src = "\n".join(lines[start:end])
    ns = {"np": np, "os": os, "re": re}
    if extra_ns:
        ns.update(extra_ns)
    exec(compile(src, f"<verbatim:{os.path.basename(path)}:{name}>", "exec"), ns)
    if name not in ns:
        raise RuntimeError(f"exec 后未得到 {name}")
    fn = ns[name]
    fn.__evs_src__ = src
    fn.__evs_origin__ = f"{os.path.basename(path)}:{start+1}-{end}"
    return fn


def exec_block(block, env, names=()):
    """逐字 exec 一个内联块，返回 names 指定的变量。

    block : (path, a, b)
    env   : 预先绑定的命名空间（np / sklearn / 生产代码里的常量和输入变量）
    names : 需要取回的变量名
    """
    path, a, b = block
    src = grab(path, a, b)
    ns = dict(env)
    ns.setdefault("np", np)
    exec(compile(src, "<verbatim-block:%s:%d-%d>" % (os.path.basename(path), a, b), "exec"), ns)
    return {k: ns.get(k, None) for k in names}, src


def block_source(block):
    path, a, b = block
    return f"{os.path.basename(path)}:{a}-{b}"


# ---------------- 结果输出 ----------------
_SEEN = []


def emit(tid, status, expected, observed, note=""):
    status = status.upper()
    assert status in ("PASS", "FAIL", "WARN"), status
    line = f"RESULT|{tid}|{status}|{expected}|{observed}|{note}"
    print(line, flush=True)
    _SEEN.append((tid, status))
    return status


def check(tid, ok, expected, observed, note="", warn_only=False):
    """ok=True -> PASS；ok=False -> FAIL（warn_only=True 时降为 WARN）。"""
    if ok:
        return emit(tid, "PASS", expected, observed, note)
    return emit(tid, "WARN" if warn_only else "FAIL", expected, observed, note)


def guard(tid, fn, expected, note=""):
    """把可能抛异常的测试包起来：抛异常 = FAIL（而不是让整族挂掉）。"""
    try:
        r = fn()
        if isinstance(r, tuple):
            ok, observed = r
        else:
            ok, observed = bool(r), str(r)
        return check(tid, ok, expected, observed, note)
    except Exception as e:
        return emit(tid, "FAIL", expected, f"EXCEPTION {type(e).__name__}: {e}",
                    (note + " | " + traceback.format_exc().splitlines()[-1]).strip(" |"))


def banner(family):
    print("=" * 100, flush=True)
    print(f"evs — {family}", flush=True)
    print("=" * 100, flush=True)
