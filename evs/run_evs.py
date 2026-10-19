# -*- coding: utf-8 -*-
r"""run_evs — estimator validation suite 聚合器 + **closure gate**

================================ 为什么必须看 RESULT 行，不能看退出码 ================================
模块按契约**恒以 exit 0 退出**（这样单个族的崩溃不会吞掉其他族）。
因此**任何按退出码判定的下游都会把整棵有 bug 的树读成全绿**。
本聚合器只解析 `RESULT|...` 行。

================================ 四态 ================================
  FAIL : hard bound 或恒等式被违反
  WARN : diagnostic（如 DR 此类未约束估计量越界）—— 只报告，**且不计为 PASS**
  PASS : 通过
  以及每族的 "无 RESULT 行 / 进程异常" 视为 BROKEN（最严重，必须显式暴露）

================================ v1.4 closure gate（CORRECTION_SPEC §6.3(3)） ================================
登记表 **ACTIVE_PATH_MANIFEST.yaml**（手写，t14 只读校验）× 本轮**实测** RESULT 行 → 覆盖表。
四条**机械**条件，全满足才允许 freeze：

  1. manifest_critical_paths  ==  coverage_table_critical_paths
  2. critical_paths_with_PASS ==  manifest_critical_paths
  3. Critical FAIL    == 0
  4. Critical UNKNOWN == 0

且每个 PASS 必须携带显式 validation_mode（来自 manifest，MC-00 已按枚举校验）。
Secondary / Rendering **不参与**上述集合等式，按 §6.1 的 residual 规则处理。

★ WARN 不计为 PASS（fail-closed）：一个 diagnostic 不能证明一条 Critical 路径。
★ 判读在哪棵树上做：REAL 上 Critical 本就应当 FAIL —— 那是**红灯基线**，不是 freeze 依据。
  故 freeze 依据是 **PATCHED（正确实现树）** 的 closure 全绿 + 每条 Critical 的 fault_injection
  注入后必须 FAIL（见 sim_fix.py 的 mutation acceptance）。

用法：python run_evs.py
      EVS_PROD_DIR=<patch tree> python run_evs.py
"""
import os, sys, io, re, csv as _csv
import concurrent.futures as _cf

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import evs_load as L          # 只为读出本次跑的是哪棵生产树

PROJ = os.path.dirname(HERE)
MANIFEST = os.path.join(PROJ, "ACTIVE_PATH_MANIFEST.yaml")

MODE = "PATCHED" if L.IS_PATCHED else "REAL"
MODULES = [
    ("evs_selftest.py",         "装载器自测（infrastructure）"),
    ("t01_mmd.py",              "MMD 估计量"),
    ("t02_concordance.py",      "AUC / concordance 并列语义"),
    ("t03_l3_estimators.py",    "IOW / OM / DR"),
    ("t04_sudo.py",             "SUDO-AURCC"),
    ("t05_unit_grain.py",       "分析单元 / 表格粒度 / 报告计数"),
    ("t06_bounds_scan.py",      "全冻结 CSV 值域扫描"),
    ("t07_fork_parity.py",      "活路径同源副本（fork parity）"),
    ("t08_psijs_family.py",     "PSI / JS / entropy"),
    ("t09_bandwidth.py",        "RBF 带宽中位数启发式"),
    ("t10_sens_at_spec.py",     "sens_at_spec80"),
    ("t12_eval_layer.py",       "评价/聚合层（bootstrap / Platt / LOO）"),
    ("t13_eval_core.py",        "评价核心（p11 / p6d / p6e）"),
    ("t14_manifest_closure.py", "manifest registry（**不执行生产文件**）"),
    ("t15_drugmap.py",          "CS-01 药物名映射"),
    ("t16_sudo_p10b.py",        "SUDO-P10B exact-identity"),
    ("t17_cs11_loo.py",         "CS-11 Fig3B 语义+provenance"),
    ("t18_reimpl.py",           "P12/FIG5 独立重实现"),
    ("t19_cs03_sr.py",          "CS-03 S/R 冲突裁决"),
    ("t22_cs04_naming.py",      "CS-04 SUDO 协议/稿件文本一致性（**不加载生产代码**）"),
    # ★ 2026-10-07 新增：补充材料**渲染件**的独立验证。
    #   动机：p21d/p21e/p21f 此前不在 closure 集合内，其产物只有渲染器的内部对账兜底 ——
    #   那是**自证**（渲染器与被渲染物同错则同绿）。t23 不 import 渲染器，
    #   自己解析 .md 并从冻结 CSV 重建每一格。
    ("t23_render_parity.py",    "补充材料渲染件独立验证（S1 / S2 / Supplementary Methods）"),
    # ★ 2026-10-07 新增（**事后**分析）：分诊策略模拟 p45 的独立验证。
    #   动机同 t23：p45 的产物要进 Supplementary Figure S10 / Table S3，
    #   不能只有 p45 自己的内部对账（自证）。
    #   t45 不 import / 不 exec p45 —— 自己从 p6b/p6c/p6d_outcome 重建排名与命中数，
    #   并校验**预注册备忘的哈希未变**（防「看完结果再改备忘」）。
    ("t45_triage.py",           "分诊策略模拟独立验证（预注册哈希 / 排名 / 曲线 / 主对比）"),
    # ★ 2026-10-07 新增：Supplementary Table S3（事后分诊）的独立验证。
    #   p46 自带 8 条闭合断言，但那是**自证**。t46 不 import p46，
    #   自己从 p6b/p6c/p6d 重建 45 个绩效格，并独立解析渲染出的 .md 逐格比对。
    ("t46_tableS3.py",          "Supplementary Table S3 独立验证（逐格重算 / 区间边界 / 舍入 / 自洽）"),
]
# 注（已 ls 核实）：t20 已被 t20_identity_gate.py 占用（该文件不在 MODULES 内，属手工运行的 gate 脚本）；
# t21 未使用；t22 为 CS-04；t23 取渲染件验证；t45 取事后分诊验证，均无冲突。
# t11_kernel_census.py 已退役 -> retired/（见 retired/README.md）

# 已知缺陷 → 应抓红它的测试（**红灯基线**判据；只在 REAL 模式下有意义）
EXPECTED_RED = {
    "MMD 因子 2（两个内项各半）": ["MMD-01", "MMD-02", "MMD-03", "MMD-04", "MMD-07", "MMD-08", "MMD-08b"],
    "wconcordance 并列记 0（应为 0.5）": ["CONC-04", "CONC-04B", "CONC-05B"],
    "DR 分子写成两残差之积": ["L3-03"],
    "SUDO 单分支（原文双分支）": ["SUDO-01", "SUDO-03", "SUDO-05"],
    "episode 粒度 / 报告层单元混淆": ["UNIT-01", "UNIT-06"],
    "MMD 越出理论下界（值域扫描）": ["BOUND-01-AGG"],
    "sens_at_spec80 分位数取反": ["SENS-03", "SENS-04", "SENS-06", "SENS-08"],
    "grouped bootstrap set() 去重": ["EVAL-A1"],
    "Platt 声明 vs 代码分歧": ["EVAL-B1"],
    "p6f LOO 段误调 lodo_oof": ["EVAL-D1"],
    # ---- fork-parity：三条**活路径**的同源副本 ----
    "DR 分子@p10b:167（活副本）": ["FORK-02", "FORK-02c"],
    "SUDO 单分支@p10b:116-148（活副本）": ["FORK-03"],
    "MMD 因子 2@p6_02:236（活副本）": ["FORK-04a", "FORK-05-MMD-p6_02-live"],
    "MMD 因子 2@其余活副本（清点）": ["FORK-05-MMD-p6b", "FORK-05-MMD-p10a", "FORK-05-MMD-p10b"],
    "wconcordance 并列@p6c/p10b（清点）": ["FORK-05-WC-p6c", "FORK-05-WC-p10b"],
    # ---- CS-05 新命名的两处偏离（此前被当成"约定差异"披露为 WARN，从未设闸）----
    "OM 分母含 i=j 对角线": ["L3-02"],
    "DR 分母为 D1*D0（边际均值之积）": ["L3-03"],
    "OM 分母含 i=j 对角线@p10b:164（活副本）": ["FORK-01"],
    # ---- CS-09：(g) 类缺陷 —— ESS 挂在错误样本域（公式对、域错）----
    #   注意 ESS-02 是 p10b 侧的**回归钉**（PASS->PASS），真正的判别项是 ESS-01。
    "ESS 样本域错（用 all-source 而非 source-test）": ["ESS-01"],
    # ---- CS-04：协议与稿件文本自称"忠于原论文"，而实现有两处已申报的偏离 ----
    #   注意这两条不是**代码**缺陷 -> 它们扫的是文档；列在这里是为了让红灯基线也覆盖文本层。
    "协议/稿件主张实现忠于原论文（实为两处偏离）": ["CS04-01", "CS04-02"],
}

rep = io.open(os.path.join(HERE, f"report_evs_{MODE.lower()}.txt"), "w", encoding="utf-8")
def p(*a):
    line = " ".join(str(x) for x in a); print(line, flush=True); rep.write(line + "\n"); rep.flush()

# 【2026-10-07 更正】头部此前写死「在**未修改的**生产代码上运行」/「真实未修改的生产代码」。
#   但 CS-05..CS-12 已于 2026-10-06 批量进生产，此后 REAL 树**就是修正后的树**，
#   该说法不实；且与本文件尾注（「本树是**修正后的生产树**」，见下方 _POSTMERGE 块）
#   自相矛盾 —— 审阅者已将此记为报告缺陷。现按同一个 _POSTMERGE 判据自动标注。
_POSTMERGE_HDR = os.path.isfile(os.path.join(
    PROJ, "FREEZE_pre_correction_20261005", "APPLIED_CORRECTIONS_20261006.txt"))
p("=" * 104)
if _POSTMERGE_HDR and MODE == "REAL":
    p("estimator validation suite —— **修正后生产树的复跑**（CS-05..CS-12 已进生产）")
else:
    p("estimator validation suite —— 红灯基线（在**未修改的**生产代码上运行）")
p("=" * 104)
p(f"生产代码树：{L.PROD_DIR}")
if MODE == "PATCHED":
    p("模式：PATCHED   ← 这是**打过补丁的副本**，用于验证「修好能否全绿」，**不得**称为红灯基线")
elif _POSTMERGE_HDR:
    p("模式：REAL      ← 生产树，**已含 CS-05..CS-12 修正**（非未修改基线；见尾注）")
else:
    p("模式：REAL      ← 真实未修改的生产代码")
p("说明：模块恒 exit 0；本表**只依据 RESULT 行**。BROKEN = 无 RESULT 行或进程异常。")
p("")
p("两种运行必须同时成立，本套件才成立：")
if _POSTMERGE_HDR and MODE == "REAL":
    p("  ① REAL    **修正前**在生产代码上 → 已知缺陷必须**全被抓红**"
      "（该次基线证据在修正期内留档；修正进生产后本节 MISSED 不携带信息，见尾注）")
else:
    p("  ① REAL    在未修改的生产代码上 → 已知缺陷必须**全被抓红**")
p("  ② PATCHED 在按正确方式修补的副本上 → 对应测试必须**转绿**，且 closure 四条全满足")
p("")

all_res = {}          # tid -> (family, worst_status)  ← 取**最差**状态，不取最后一行
_SEV = {"PASS": 0, "WARN": 1, "FAIL": 2}
fam_rows = []

JOBS = int(os.environ.get("EVS_JOBS") or max(1, min(8, (os.cpu_count() or 4) - 2)))
ONLY = [s.strip() for s in (os.environ.get("EVS_ONLY") or "").split(",") if s.strip()]
_run_mods = [m for m in MODULES if (not ONLY or m[0] in ONLY)]

def _run_one(fn):
    path = os.path.join(HERE, fn)
    if not os.path.exists(path):
        return fn, None, "", ""
    import subprocess
    r = subprocess.run([sys.executable, path], capture_output=True, text=True,
                       encoding="utf-8", errors="replace", cwd=HERE)
    return fn, r.returncode, (r.stdout or ""), (r.stderr or "")

p(f"  （并行执行 {len(_run_mods)} 个模块，workers={JOBS}"
  + (f"，EVS_ONLY={ONLY}" if ONLY else "") + "）")
with _cf.ThreadPoolExecutor(max_workers=JOBS) as _ex:
    _outs = {t[0]: t for t in _ex.map(_run_one, [m[0] for m in _run_mods])}

_struct = []          # 结构化行：module, test_id, status, expected, observed, note
for fn, label in _run_mods:
    rc, out, err = _outs[fn][1], _outs[fn][2], _outs[fn][3]
    lines = [l for l in out.split("\n") if l.startswith("RESULT|")]
    nf = sum(1 for l in lines if l.split("|")[2] == "FAIL")
    nw = sum(1 for l in lines if l.split("|")[2] == "WARN")
    np_ = sum(1 for l in lines if l.split("|")[2] == "PASS")
    for l in lines:
        parts = l.split("|")
        if len(parts) >= 6:
            tid, st = parts[1], parts[2]
            _struct.append([fn, tid, st, parts[3], parts[4], parts[5]])
            prev = all_res.get(tid)
            if prev is None or _SEV[st] > _SEV[prev[1]]:
                all_res[tid] = (fn, st)
    tail = err.strip().split("\n")[-1][:80] if rc not in (0, None) else ""
    state = "BROKEN" if (not lines or rc != 0) else ("RED" if nf else "green")
    fam_rows.append((fn, label, state, nf, nw, np_, tail))
for fn, label in MODULES:
    if not any(r[0] == fn for r in fam_rows):
        fam_rows.append((fn, label, "MISSING", 0, 0, 0, ""))

_struct_path = os.path.join(HERE, f"evs_results_{MODE.lower()}.csv")
with io.open(_struct_path, "w", encoding="utf-8", newline="") as _f:
    _w = _csv.writer(_f)
    _w.writerow(["module", "test_id", "status", "expected", "observed", "note"])
    _w.writerows(_struct)

# ---- 自检：本脚本自身的输出必须完整 ----
assert len(fam_rows) == len(MODULES), \
    f"模块表不完整：{len(fam_rows)} != {len(MODULES)}（本脚本自身有 bug，报告不可引用）"
assert all_res, "未收集到任何 RESULT 行（本脚本自身有 bug，报告不可引用）"
assert _struct, f"结构化结果为空：{_struct_path}（下游会读到空表，报告不可引用）"
_broken = [r[0] for r in fam_rows if r[2] in ("BROKEN", "MISSING")]
if _broken:
    p("!! 以下模块未产出可解析结果，整个报告不可引用：" + str(_broken))

p(f"  {'module':<26}{'family':<30}{'state':<9}{'FAIL':>6}{'WARN':>6}{'PASS':>6}")
for fn, label, state, nf, nw, np_, tail in fam_rows:
    mark = {"RED": "!!", "green": "  ", "BROKEN": "XX", "MISSING": "XX"}[state]
    p(f"  {mark}{fn:<24}{label:<30}{state:<9}{nf:>6}{nw:>6}{np_:>6}")
    if tail: p(f"      stderr: {tail}")

tot_f = sum(r[3] for r in fam_rows); tot_w = sum(r[4] for r in fam_rows); tot_p = sum(r[5] for r in fam_rows)
p("")
p(f"  合计：FAIL {tot_f} ｜ WARN {tot_w} ｜ PASS {tot_p}")

# ==================================================================================
# closure gate（§6.3(3)）
# ==================================================================================
p("")
p("=" * 104)
p("【closure gate】manifest（手写登记表）× 本轮实测 RESULT 行")
p("=" * 104)

# run_evs 自己读 manifest 只为生成覆盖表；**不执行任何生产文件**（t14 的结构断言另管 registry 侧）
import yaml
_cov_rows, _crit_ids, _cov_crit, _pass_crit = [], [], set(), set()
_crit_fail = _crit_unknown = 0
_unreg_tests = set()
_manifest_err = ""
try:
    with io.open(MANIFEST, encoding="utf-8") as _f:
        _doc = L.load_manifest(MANIFEST)   # ★ 共享受保护 loader：重复键硬失败
    _paths = _doc["paths"]
except Exception as _e:                                   # noqa: BLE001
    _manifest_err = f"{type(_e).__name__}: {_e}"
    _paths = []

_PHASE_DEFAULT = "pre_correction_freeze"
_SEV_ST = {"PASS": 0, "UNKNOWN": 1, "FAIL": 2}


def _status_of(tests):
    """一个测试集合的状态。WARN 与 ABSENT 都 fail-closed 记 UNKNOWN。"""
    sts = [all_res.get(x, (None, "ABSENT"))[1] for x in tests]
    c = {k: sum(1 for s in sts if s == k) for k in ("PASS", "FAIL", "WARN", "ABSENT")}
    if not tests or c["ABSENT"]:
        st = "UNKNOWN"
    elif c["FAIL"]:
        st = "FAIL"
    elif c["WARN"]:
        st = "UNKNOWN"
    else:
        st = "PASS"
    return st, c


if _paths:
    # ★ 2026-10-06 口径对齐（与 t14 MC-10 的 _known() 同口径）：
    #   原实现只把 required_tests 与 deferred_gates 计入 _referenced，
    #   **characterisations / pending_gates / unbound_deferred 全都漏了**，且**不支持前缀登记**
    #   （如 `BOUND-01::*`）。后果：120 条**已经登记**的 test_id 被报成"未登记"，
    #   真信号被噪声淹没 —— 与前面 FORK-01 未登记、MC-10 自锁是同一类"同一想法的两份实现不一致"。
    _referenced = set()
    for _e in _paths:
        _rt = list(_e.get("required_tests") or [])
        _dg = list(_e.get("deferred_gates") or [])
        _pg = list(_e.get("pending_gates") or [])
        _ch = list(_e.get("characterisations") or [])
        _referenced.update(_rt)
        _referenced.update(_ch)
        _referenced.update(_pg)
        _referenced.update(g["test"] for g in _dg)
        _phase = _e.get("phase") or _PHASE_DEFAULT
        # phase_status：**只看本阶段的 required_tests**（阶段落实到测试级，不整条 path 后移）
        _st, _c = _status_of(_rt)
        # deferred_status：后移阶段的门（跨阶段才允许存在，MC-00 已强制）
        _dst, _dc = (_status_of([g["test"] for g in _dg]) if _dg else (None, None))
        _overall = _st if _dst is None else max(_st, _dst, key=lambda s: _SEV_ST[s])
        _cov_rows.append([
            _e["path_id"], _e.get("risk_class"), _e.get("validation_mode"), _phase,
            _st, _overall,
            len(_rt), _c["PASS"], _c["FAIL"], _c["WARN"], _c["ABSENT"],
            (_dst or ""), ";".join(g["test"] for g in _dg),
            ";".join(g["test"] for g in _pg), ";".join(_rt),
            ";".join(f"{x}={all_res.get(x, (None, 'ABSENT'))[1]}" for x in _ch)])
        if _e.get("risk_class") == "Critical" and _phase == _PHASE_DEFAULT:
            _crit_ids.append(_e["path_id"])
            if _st == "PASS":
                _pass_crit.add(_e["path_id"])
            _crit_fail += int(_st == "FAIL")
            _crit_unknown += int(_st == "UNKNOWN")
    # ★ 2026-10-06：unbound_deferred 与**前缀登记**在循环之后补进 _referenced ——
    #   前者不属于任何 path，后者（如 `BOUND-01::*`）只需登记一次即覆盖整族。
    #   与 t14 MC-10 的 _known() 同口径。
    _referenced.update(g.get("test") for g in (_doc.get("unbound_deferred") or []))
    _referenced.discard(None)
    _ref_prefix = tuple(sorted(t[:-1] for t in _referenced
                               if isinstance(t, str) and t.endswith("*")))
    # coverage_table_critical_paths = **表里出现的全部 Critical ID**（不论阶段）。
    # 旧实现用 `n_obs > 0` 过滤，把「没有测试跑过」和「不在表里」混成同一件事，
    # 于是 SUDO_P10B / LOSO 从构造上永远缺席 —— 那不是 coverage hole，是**表的缺陷**。
    # 集合等式两侧都取**本阶段**集合（表里仍保留全部 Critical 行，靠 phase 列区分，
    # 于是 LOSO / SUDO_P10B 这类"未实现/他阶段"的项在表里**可见**却不污染本阶段判据）。
    _cov_crit = {r[0] for r in _cov_rows
                 if r[1] == "Critical" and r[3] == _PHASE_DEFAULT}
    _dup_ids = sorted({r[0] for r in _cov_rows if [x[0] for x in _cov_rows].count(r[0]) > 1})
    _crit_ids = sorted(set(_crit_ids))
    _man_crit = set(_crit_ids)
    # 基础设施自测（装载器 LOAD-*、registry MC-*）不属于任何生产 path，不计入「未登记」。
    _INFRA = ("LOAD-", "MC-")
    _unreg_tests = {t for t in all_res
                    if t not in _referenced
                    and not (_ref_prefix and t.startswith(_ref_prefix))
                    and not t.startswith(_INFRA)}

_cov_path = os.path.join(HERE, f"coverage_table_{MODE.lower()}.csv")
with io.open(_cov_path, "w", encoding="utf-8", newline="") as _f:
    _w = _csv.writer(_f)
    _w.writerow(["path_id", "risk_class", "validation_mode", "phase",
                 "phase_status", "overall_status",
                 "n_required", "n_pass", "n_fail", "n_warn", "n_absent",
                 "deferred_status", "deferred_tests", "pending_gates",
                 "required_tests", "characterisations"])
    _w.writerows(_cov_rows)

# §6.3(3)：**集合相等 + 无重复**，不是数个数
_man_crit, _cov_crit_s = set(_crit_ids), set(_cov_crit)
c1 = (_man_crit == _cov_crit_s) and not _dup_ids
if _dup_ids:
    p(f"  !! coverage table 有重复 path_id：{_dup_ids}")
c2 = (_pass_crit == _man_crit)
c3 = (_crit_fail == 0)
c4 = (_crit_unknown == 0)

def _cond(flag, text, detail):
    p(f"  [{'OK ' if flag else 'FAIL'}] {text}")
    p(f"        {detail}")
    return flag

p(f"  覆盖表：{_cov_path}（{len(_cov_rows)} 条 path）")
if _manifest_err:
    p(f"  !! manifest 不可读：{_manifest_err}")
p("")
p("  §6.3(3) 四条机械条件：")
_ok1 = _cond(c1, "1. manifest_critical_paths == coverage_table_critical_paths",
             f"manifest {len(_man_crit)} ｜ 覆盖表（有实测）{len(_cov_crit_s)}"
             + ("" if c1 else f" ｜ 只在 manifest: {sorted(_man_crit - _cov_crit_s)}"
                              f" ｜ 只在覆盖表: {sorted(_cov_crit_s - _man_crit)}"))
_ok2 = _cond(c2, "2. critical_paths_with_PASS == manifest_critical_paths",
             f"PASS {len(_pass_crit)}/{len(_man_crit)}"
             + ("" if c2 else f" ｜ 未 PASS: {sorted(_man_crit - _pass_crit)}"))
_ok3 = _cond(c3, "3. Critical FAIL == 0", f"FAIL = {_crit_fail}")
_ok4 = _cond(c4, "4. Critical UNKNOWN == 0", f"UNKNOWN = {_crit_unknown}")

p("")
p("  Critical 逐路径明细（phase / validation_mode 必须显式；本阶段集合 = Critical ∩ phase==pre_correction_freeze）：")
for _r in _cov_rows:
    if _r[1] != "Critical":
        continue
    p(f"    {_r[4]:<8} {_r[0]:<24} phase={_r[3]:<21} mode={_r[2]:<19} "
      f"tests {_r[7]}P/{_r[8]}F/{_r[9]}W/{_r[10]}A"
      + (f"  deferred[{_r[11]}]: {_r[12]}" if _r[12] else "")
      + (f"  PENDING: {_r[13]}" if _r[13] else ""))
_p_no_mode = [r[0] for r in _cov_rows if r[1] == "Critical" and not r[2]]
if _p_no_mode:
    p(f"  !! Critical 缺 validation_mode：{_p_no_mode}")
_p_warn_crit = [r[0] for r in _cov_rows if r[1] == "Critical" and r[10] == 0 and r[9] > 0]
if _p_warn_crit:
    p(f"  !! 仅靠 WARN 的 Critical 路径（fail-closed 计为 UNKNOWN）：{_p_warn_crit}")
p("")
p("  Secondary / Rendering（**不参与**集合等式，按 §6.1 residual 规则）：")
for _r in _cov_rows:
    if _r[1] != "Critical":
        p(f"    {_r[4]:<8} {_r[0]:<24} {_r[1]:<10} mode={_r[2]}")
if _unreg_tests:
    p("")
    p(f"  !! 有 RESULT 行但**未被任何 manifest 条目引用**的 test_id（{len(_unreg_tests)} 条）——")

    p("     登记表可能漏了它们所属的路径：" + ", ".join(sorted(_unreg_tests)[:20]))

p("")
p(f"  closure verdict [{MODE}]："
  + ("**四条全满足** → 允许 freeze" if (_ok1 and _ok2 and _ok3 and _ok4 and not _manifest_err)
     else "**未满足** → 不得 freeze；只对未 PASS/UNKNOWN 的 Critical path 做定点闭合"))

# ==================================================================================
p("")
p("=" * 104)
p("【红灯基线】每个已知缺陷是否被独立抓到")
p("=" * 104)
miss = []
for defect, tids in EXPECTED_RED.items():
    got = [(t, all_res.get(t, ("?", "ABSENT"))[1]) for t in tids]
    nfail = sum(1 for _, s in got if s == "FAIL")
    ok = nfail > 0
    p(f"  [{'CAUGHT' if ok else 'MISSED'}] {defect}")
    p(f"       抓到它的测试 = {', '.join(f'{t}({s})' for t, s in got if s == 'FAIL') or '（无）'}")
    if not ok: miss.append(defect)

p("")
p("=" * 104)
if miss:
    p(f"*** 有 {len(miss)} 个已知缺陷未被抓到 —— 套件不完整，必须修 ***")
    for m_ in miss: p(f"    MISSED: {m_}")
else:
    p("*** 全部已知缺陷均被独立抓到：红灯基线成立 ***")
# ★ 2026-10-06：CS-05..CS-12 已批量进生产，此后 REAL 树**就是修正后的树**。
#   于是"缺陷再也抓不到"成为**预期结果**，而不是失败。原判定只看 MODE=="PATCHED"，
#   合并后不再触发 ⇒ [MISSED] 会无解释地露出来，读者会误以为判据失效。
#   现按"是否存在合并存档"判定当前树状态，并给出对应说明。
_POSTMERGE = os.path.isfile(os.path.join(
    PROJ, "FREEZE_pre_correction_20261005", "APPLIED_CORRECTIONS_20261006.txt"))
if MODE == "PATCHED":
    p("")
    p("注：PATCHED 模式下 EXPECTED_RED 判据**不适用**（缺陷已被有意修好，本就该 MISSED）。")
    p("    本节的 MISSED 在此模式下**不构成失败**；判定见 sim_fix 的 REAL/PATCHED 三态表与 closure gate。")
elif _POSTMERGE:
    p("")
    p("注：本树是**修正后的生产树**（CS-05..CS-12 已于 2026-10-06 批量进生产；"
      "存档 FREEZE_.../APPLIED_CORRECTIONS_20261006.txt）。")
    p("    ⇒ 本节按设计**应当全部 MISSED**（缺陷已不存在，任何判据都不可能再抓到它们），"
      "**不构成失败**。")
    p("    红灯基线（判据能否抓到缺陷）的证据在修正期内已逐轮留档；进生产后本节的 MISSED "
      "不再携带信息。")
p("=" * 104)
p("")
p("注：WARN 不计入不通过，但**也不计为 PASS**（closure gate 里 fail-closed 计为 UNKNOWN）。")
p("    DR 此类未约束 augmented 估计量越界属 diagnostic，**不得**自动判为 implementation error。")
rep.close(); print("\ndone")

# 【2026-10-07 更正】退出码此前与 closure verdict **互相矛盾**：
#   修正进生产后，REAL 树上 EXPECTED_RED 那 18 项按设计**必然全部 MISSED**（尾注已说明
#   它们不构成失败），但第 428 行仍把 "miss 非空" 判为失败 → 脚本以 exit 1 结束，
#   而同一份报告的 closure gate 写的是「四条全满足 → 允许 freeze」。
#   任何按退出码判定的下游会把一个「允许 freeze」的状态读成失败 —— 与头部那处
#   （已同步更正）同属「同一状态三处信号不一致」。
#   现改为：**退出码由 closure verdict 决定**；MISSED 是否算失败，取决于当前树是否已进修正
#   （未进修正时 MISSED 仍是真失败，保持原语义）。
_broken_exit = bool(_broken)
if MODE == "PATCHED":
    sys.exit(1 if (_broken_exit or not (_ok1 and _ok2 and _ok3 and _ok4)) else 0)
if _POSTMERGE_HDR:
    # 已进修正：MISSED 为预期，退出码只看 closure 四条 + BROKEN
    sys.exit(1 if (_broken_exit or not (_ok1 and _ok2 and _ok3 and _ok4)) else 0)
sys.exit(0 if (not miss and not _broken_exit) else 1)
