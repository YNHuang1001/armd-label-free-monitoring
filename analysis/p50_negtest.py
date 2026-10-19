# -*- coding: utf-8 -*-
r"""p50 的负向测试 —— **在副本上跑**（外部核验的三合一修法）

修法要点（一次解决三件事）：
  ① **在副本上跑** —— 不再改真 sidecar（此前两次直接改真身，第一次成了 INC-011）
  ② **tee 落盘** —— 发射留下耐久证据；「守卫咬过」不再只活在消息里
  ③ **引用时自带围栏** —— 日志本身就是可直接粘贴的工件

用法：python p50_negtest.py   （退出码 0 = 三个用例全部按预期判红）
产出：log_p50_negtest_<date>.txt
"""
import io
import json
import os
import shutil
import subprocess
import sys
import tempfile
from datetime import date

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REAL = os.path.join(ROOT, "figures", "data")
LOG = os.path.join(HERE, f"log_p50_negtest_{date.today().strftime('%Y%m%d')}.txt")

def _mutate(d, artifact, field, value):
    """在**副本**里改 sidecar 的某个字段。"""
    p = os.path.join(d, "reporting_profiles", artifact + ".gen.json")
    g = json.load(io.open(p, encoding="utf-8"))
    g[field] = value
    io.open(p, "w", encoding="utf-8").write(json.dumps(g, ensure_ascii=False))


CASES = []


def run(env):
    r = subprocess.run([sys.executable, os.path.join(HERE, "p50_dependency_edges.py")],
                       capture_output=True, text=True, encoding="utf-8", env=env)
    return r.returncode, r.stdout


lines = ["=" * 100,
         "p50 负向测试（**在副本上跑**；真产物未被触碰）",
         "=" * 100, ""]
ok_all = True

for name, mutate in [
    ("A. 改错 artifact_sha ⇒ 应判「盖章后被覆盖」",
     lambda d: _mutate(d, "site_profile.csv", "artifact_sha", "deadbeefdeadbeef")),
    ("B. 改错 universe_sha ⇒ 应判 STALE",
     lambda d: _mutate(d, "site_profile.csv", "universe_sha", "deadbeefdeadbeef")),
    ("C. 把产出者改成非声明收尾者 ⇒ 应判「错序」（Q3 检测）",
     lambda d: _mutate(d, "site_profile.csv", "producer", "p21a")),
]:
    tmp = tempfile.mkdtemp(prefix="p50neg_")
    try:
        shutil.copytree(os.path.join(REAL, "reporting_profiles"),
                        os.path.join(tmp, "reporting_profiles"))
        for f in ("tableS1_tasks.csv", "tableS3_triage.csv"):
            s = os.path.join(REAL, f)
            if os.path.exists(s):
                shutil.copy(s, os.path.join(tmp, f))
                for ext in (".gen.json",):
                    if os.path.exists(s + ext):
                        shutil.copy(s + ext, os.path.join(tmp, f + ext))
        mutate(tmp)
        env = dict(os.environ); env["P50_DATA_ROOT"] = tmp
        rc, out = run(env)
        hit = [l for l in out.splitlines() if "✘" in l]
        good = (rc != 0) and bool(hit)
        ok_all &= good
        lines += [f"【{name}】", f"  退出码 {rc} ｜ 判红行 {len(hit)} 条"]
        for l in hit[:1]:
            lines.append("   " + l.strip())
        lines += [f"  → {'✔ 按预期判红' if good else '✘ **未按预期**'}", ""]
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

lines += ["=" * 100,
          f"*** 负向测试：{'三项全部按预期判红' if ok_all else '**存在未按预期的用例**'} ***"]
io.open(LOG, "w", encoding="utf-8").write("\n".join(lines) + "\n")
print("\n".join(lines))
sys.exit(0 if ok_all else 1)
