# -*- coding: utf-8 -*-
r"""mutate_acceptance — **定点闭合的证明装置**（CORRECTION_SPEC §6.3(3) 的第二条腿）

================================ 要证什么 ================================
GPT6 的裁决（2026-10-04）：

    「t10/t12/t13 的修复目标不要写成"所有攻击都抓住"。只需要证明：
      **一个已登记 Critical implementation 被按已知错误方式修改时，对应 validation 必须 FAIL；
      正确实现必须 PASS。** 换句话说，只关 coverage hole，不要继续测试 census 自己能否
      识别千奇百怪的隐藏 fork。」

于是本驱动器对 **manifest 的每一条 Critical path** 做一件事：

    在正确树（_patched）上，按该 path 在 `fault_injection` 里**自己申报**的已知错误方式注入故障，
    然后重跑该 path 的 `required_tests`，要求覆盖表把这条 path 判成 **FAIL**。

  ★ 这不是"发现"机制：故障是 manifest **手写申报**的；驱动器只做执行与判读。
  ★ 也不再问"能不能找到第 34 个假通过"——那已不是本项目的 freeze criterion。

================================ 两种注入 ================================
  kind: revert-to-production   把该 path 涉及的**真实生产文本**覆盖到变异树上
                               （"退回未修正状态"本身就是已登记的已知错误）
  kind: synthetic              按 manifest 里的 (find, replace) 字面替换
                               （替换必须**恰好命中一次**，否则该条判 UNEXECUTABLE，绝不静默跳过）

================================ 判读 ================================
  CAUGHT        该 path 在变异树上被 required_tests 判为 FAIL  ← 目标
  ESCAPED       变异后仍 PASS/UNKNOWN                ← **coverage hole，必须呈现**
  UNEXECUTABLE  find 串未命中 / 文件不在树上 / 无 required_tests ← **必须呈现，不得静默**

用法：
  python mutate_acceptance.py                  # 全部 Critical
  python mutate_acceptance.py --only L1_MMD
  EVS_MUT_JOBS=4 python mutate_acceptance.py
输出：mutate_acceptance_report.txt / mutate_acceptance.csv
=============================================================================
"""
import os, sys, io, re, csv, json, shutil, subprocess, hashlib, argparse

import yaml

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(os.path.dirname(HERE), "分析结果")                       # 真实生产树
PROJ = os.path.dirname(HERE)
MANIFEST = os.path.join(PROJ, "ACTIVE_PATH_MANIFEST.yaml")
PATCHED = os.path.join(HERE, "_patched")
MUTROOT = os.path.join(HERE, "_mut")
PY = sys.executable

rep = io.open(os.path.join(HERE, "mutate_acceptance_report.txt"), "w", encoding="utf-8")


def p(*a):
    line = " ".join(str(x) for x in a)
    print(line, flush=True)
    rep.write(line + "\n")
    rep.flush()


def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 16), b""):
            h.update(b)
    return h.hexdigest()


ap = argparse.ArgumentParser()
ap.add_argument("--only", default="")
ap.add_argument("--keep", action="store_true", help="保留变异树以便复查")
args = ap.parse_args()

p("=" * 104)
p("mutate_acceptance — 对 manifest 每条 Critical path 注入其**自己申报**的已知错误")
p("=" * 104)
p(f"  正确树(_patched) = {PATCHED}")
p(f"  真实生产树        = {OUT}")
p(f"  manifest          = {MANIFEST}")
p("")

if not os.path.isdir(PATCHED):
    p("!! _patched 不存在 —— 请先跑 sim_fix.py 生成正确树（本驱动器只做注入，不生成正确树）")
    rep.close()
    sys.exit(2)

with io.open(MANIFEST, encoding="utf-8") as f:
    DOC = yaml.safe_load(f)
PATHS = DOC["paths"]

# test_id -> module（从真实跑的 RESULT 行取，缺失即视为该 test 不存在）
# ★ 只接受 **run_evs 的 schema**：`module` 必须是形如 `tNN_*.py` 的文件名。
#   事故（2026-10-05）：`sim_fix.py` 曾用**自己的 schema**（`module` = 族名，如 `MMD`）覆写
#   同名文件，使本驱动器把族名当模块名，`EVS_ONLY='FORK,MMD,...'` 匹配不到任何模块
#   ⇒ 覆盖表全 UNKNOWN ⇒ **7 条假的 ESCAPED**。现在 schema 不符即拒收并报错，不再静默。
def _run_evs_modules():
    """从 run_evs.py **源码**取 MODULES 清单（AST，不 import、不执行）。"""
    import ast as _ast
    src = io.open(os.path.join(HERE, "run_evs.py"), encoding="utf-8").read()
    for node in _ast.walk(_ast.parse(src)):
        if isinstance(node, _ast.Assign) and any(
                getattr(tg, "id", None) == "MODULES" for tg in node.targets):
            return {el.elts[0].value for el in node.value.elts}
    raise SystemExit("!! run_evs.py 中找不到 MODULES")


RUN_MODULES = _run_evs_modules()
TID2MOD = {}
for cand in ("evs_results_real.csv", "evs_results_patched.csv"):
    q = os.path.join(HERE, cand)
    if not os.path.exists(q):
        continue
    with io.open(q, encoding="utf-8") as f:
        rd = csv.DictReader(f)
        if "module" not in (rd.fieldnames or []):
            raise SystemExit(f"!! {cand} 缺 module 列（schema 不符，拒绝使用）")
        for r in rd:
            m = r["module"]
            # ★ 用 run_evs 的**显式 MODULES 清单**校验（AST 读源码，不 import），
            #   而不是名字正则——`evs_selftest.py` 是合法模块名却不符合 ^t\d\d_。
            #   这正是本项目自己的原则：封闭世界用显式 registry，不用模式。
            if m not in RUN_MODULES:
                raise SystemExit(f"!! {cand} 的 module={m!r} 不在 run_evs.MODULES 中"
                                 f"（疑似被别的写者覆写）。请先跑 run_evs 重建，再跑本驱动器。")
            TID2MOD.setdefault(r["test_id"], m)
p(f"  test_id -> module 映射：{len(TID2MOD)} 条（来自 evs_results_*.csv）")
p("")

CRIT = [e for e in PATHS if e.get("risk_class") == "Critical"]
if args.only:
    CRIT = [e for e in CRIT if e["path_id"] == args.only]
    p(f"  --only {args.only} -> {len(CRIT)} 条")
p(f"  待注入的 Critical path：{len(CRIT)}")
p("")


# ------------------------------------------------------------------ 注入实现
def files_of(entry):
    fs = [s.strip() for s in str(entry.get("file", "")).split(";") if s.strip()]
    return [f for f in fs if not f.startswith("(") and "*" not in f]


def apply_revert(tree, entry):
    """把该 path 涉及的**真实生产文本**覆盖到变异树上（退回未修正状态）。"""
    fs = files_of(entry)
    if not fs:
        return None, "该条 file 字段不是具体文件清单，revert 无法定位"
    done = []
    for f in fs:
        src, dst = os.path.join(OUT, f), os.path.join(tree, f)
        if not os.path.exists(src):
            return None, f"真实文件不存在：{f}"
        if not os.path.exists(dst):
            return None, f"变异树上没有 {f}（不在 _patched 里）"
        if sha(src) == sha(dst):
            continue                       # 该文件本就没有补丁，退回无意义但也不算错
        shutil.copy2(src, dst)
        done.append(f)
    if not done:
        return None, "该 path 涉及的文件在 _patched 上均**与生产一致**（没有可退回的补丁）"
    return done, ""


def apply_synthetic(tree, fi):
    """按 (find, replace) 做**恰好一次**的字面替换。"""
    f = fi.get("file")
    find, repl = fi.get("find"), fi.get("replace")
    if not f or find is None or repl is None:
        return None, "synthetic 注入缺少 file/find/replace 申报"
    dst = os.path.join(tree, f)
    if not os.path.exists(dst):
        return None, f"变异树上没有 {f}"
    with io.open(dst, encoding="utf-8", newline="") as fh:
        raw = fh.read()
    n = raw.count(find)
    if n != 1:
        return None, f"find 串在 {f} 中命中 {n} 次（必须恰好 1 次）"
    with io.open(dst, "w", encoding="utf-8", newline="") as fh:
        fh.write(raw.replace(find, repl, 1))
    return [f], ""


# ------------------------------------------------------------------ 逐条执行
rows = []
for e in CRIT:
    pid = e["path_id"]
    fi = e.get("fault_injection") or {}
    req = list(e.get("required_tests") or [])
    mods = sorted({TID2MOD[t] for t in req if t in TID2MOD})
    miss_tid = [t for t in req if t not in TID2MOD]

    p("-" * 104)
    p(f"[{pid}]  risk=Critical  mode={e.get('validation_mode')}  fault={fi.get('id')}")
    p(f"    required_tests = {len(req)} 条；可解析到模块 {len(mods)} 个；未在结果中出现 {len(miss_tid)} 条"
      + (f" :: {miss_tid[:6]}" if miss_tid else ""))

    rec = dict(path_id=pid, fault_id=fi.get("id", ""), kind=fi.get("kind", ""),
               required_tests=len(req), modules=";".join(mods), verdict="", detail="",
               cov_status="", cov_pass=0, cov_fail=0, cov_absent=0)

    if not req:
        rec.update(verdict="UNEXECUTABLE", detail="无 required_tests（deferred）—— 该 path 在 closure 中计 UNKNOWN")
        p(f"    -> UNEXECUTABLE :: {rec['detail']}")
        rows.append(rec)
        continue
    if not mods:
        rec.update(verdict="UNEXECUTABLE",
                   detail="required_tests 无一条出现在 evs_results_*.csv 中（测试根本不存在或未运行）")
        p(f"    -> UNEXECUTABLE :: {rec['detail']}")
        rows.append(rec)
        continue
    if not fi:
        rec.update(verdict="UNEXECUTABLE", detail="无 fault_injection 申报")
        p(f"    -> UNEXECUTABLE :: {rec['detail']}")
        rows.append(rec)
        continue

    tree = os.path.join(MUTROOT, pid)
    if os.path.isdir(tree):
        shutil.rmtree(tree)
    os.makedirs(MUTROOT, exist_ok=True)
    shutil.copytree(PATCHED, tree)

    kind = str(fi.get("kind", ""))
    if kind == "revert-to-production":
        done, err = apply_revert(tree, e)
    elif kind == "synthetic":
        done, err = apply_synthetic(tree, fi)
    else:
        done, err = None, f"未知 kind {kind!r}"

    if err:
        rec.update(verdict="UNEXECUTABLE", detail=err)
        p(f"    -> UNEXECUTABLE :: {err}")
        rows.append(rec)
        continue
    p(f"    注入完成：{done}")

    # 跑该 path 的 required_tests（+ t14 以便同时复核 registry 未被变异破坏）
    env = os.environ.copy()
    env["PYTHONIOENCODING"] = "utf-8"
    env["EVS_PROD_DIR"] = tree
    env["EVS_ONLY"] = ",".join(sorted(set(mods) | {"t14_manifest_closure.py"}))
    r = subprocess.run([PY, os.path.join(HERE, "run_evs.py")], cwd=HERE, env=env,
                       capture_output=True, text=True, encoding="utf-8", errors="replace")
    if r.returncode not in (0, 1):
        rec.update(verdict="UNEXECUTABLE", detail=f"run_evs rc={r.returncode}")
        p(f"    -> UNEXECUTABLE :: run_evs rc={r.returncode} :: {(r.stderr or '')[-300:]}")
        rows.append(rec)
        continue

    cov = os.path.join(HERE, "coverage_table_patched.csv")
    status = None
    if os.path.exists(cov):
        with io.open(cov, encoding="utf-8") as fh:
            for cr in csv.DictReader(fh):
                if cr["path_id"] == pid:
                    status = cr["status"]
                    rec.update(cov_pass=int(cr["n_pass"]), cov_fail=int(cr["n_fail"]),
                               cov_absent=int(cr["n_absent"]))
                    break
        shutil.copy2(cov, os.path.join(tree, "coverage_table.csv"))
    rec["cov_status"] = status or "NO-ROW"

    if status == "FAIL":
        rec.update(verdict="CAUGHT", detail="变异后该 path 被判 FAIL（required_tests 有 FAIL）")
    elif status is None:
        rec.update(verdict="UNEXECUTABLE", detail="覆盖表里没有这条 path")
    else:
        rec.update(verdict="ESCAPED",
                   detail=f"变异后该 path 仍为 {status}（pass={rec['cov_pass']} fail={rec['cov_fail']} "
                          f"absent={rec['cov_absent']}）—— **coverage hole**")
    p(f"    -> {rec['verdict']} :: 覆盖表 status={status} "
      f"(pass={rec['cov_pass']} fail={rec['cov_fail']} absent={rec['cov_absent']})")
    if rec["verdict"] == "ESCAPED":
        hits = [l for l in (r.stdout or "").split("\n") if f"RESULT|" in l
                and l.split("|")[2] in ("FAIL", "WARN")]
        for h in hits[:6]:
            p(f"       {h[:200]}")
    rows.append(rec)
    if not args.keep:
        shutil.rmtree(tree, ignore_errors=True)

# ------------------------------------------------------------------ 汇总
p("")
p("=" * 104)
p("【汇总】每条 Critical path 的注入判读")
p("=" * 104)
p(f"  {'path_id':<24}{'fault':<24}{'verdict':<14}{'cov_status':<10}detail")
for r in rows:
    p(f"  {r['path_id']:<24}{r['fault_id']:<24}{r['verdict']:<14}{r['cov_status']:<10}{r['detail'][:90]}")

n = {k: sum(1 for r in rows if r["verdict"] == k) for k in ("CAUGHT", "ESCAPED", "UNEXECUTABLE")}
p("")
p(f"  CAUGHT {n['CAUGHT']} ｜ ESCAPED {n['ESCAPED']} ｜ UNEXECUTABLE {n['UNEXECUTABLE']} ｜ 共 {len(rows)}")
p("")
if n["ESCAPED"] or n["UNEXECUTABLE"]:
    p("  *** 未闭合：ESCAPED = coverage hole（变异后仍然 PASS）；UNEXECUTABLE = 无法注入或测试不存在。")
    p("      两者都必须在 closure report 中呈现，**不得静默略过**。 ***")
else:
    p("  *** 全部 Critical path 均满足：注入已知错误 ⇒ validation 判 FAIL ***")

with io.open(os.path.join(HERE, "mutate_acceptance.csv"), "w", encoding="utf-8", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()) if rows else
                       ["path_id", "fault_id", "kind", "required_tests", "modules",
                        "verdict", "detail", "cov_status", "cov_pass", "cov_fail", "cov_absent"])
    w.writeheader()
    w.writerows(rows)
rep.close()
print("\ndone", flush=True)
