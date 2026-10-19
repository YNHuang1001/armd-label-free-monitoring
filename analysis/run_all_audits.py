# -*- coding: utf-8 -*-
r"""run_all_audits — 统一运行全部审计器，并**验证报告是本次生成的**

================================ 为什么需要它 ================================
2026-10-02 事故：我给 p16 打补丁时，heredoc 把 `\n` 变成了真换行 → **SyntaxError**。
该错误发生在**解析期**，此时 `open(report, "w")` 尚未执行，
**磁盘上留着上一次运行的完整报告**。
我以 `>/dev/null` 运行、只看了 exit code（1），然后去读那份报告 ——
**读到的是陈旧内容，并据此以为问题仍然存在。**

exit code 只能说明"这次失败了"，**不能说明"你读到的是这次的结果"**。

================================ 本编排器做三件事 ================================
1. 逐脚本运行，记录 exit code；
2. **校验每个报告文件的 mtime 晚于本次运行开始时刻** —— 否则标 `STALE`；
3. 汇总 PASS / FAIL / STALE / MISSING 四态，任一非 PASS 即整体不通过。

用法：python run_all_audits.py
"""
import os, sys, time, subprocess

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)

# (脚本, 报告文件)  —— 报告一律在 分析结果/ 下
AUDITS = [
 ("p15_results_numeric_audit.py",              "report_p15_results_numeric_audit.txt"),
 ("p16_methods_audit.py",                      "report_p16_methods_audit.txt"),
 ("p17_sidecar_provenance_audit.py",           "report_p17_sidecar_provenance.txt"),
 ("p18_supplementary_citation_check.py",       "report_p18_supplementary_citation.txt"),
 ("p19_discussion_audit.py",                   "report_p19_discussion_audit.txt"),
 ("p20_style_normalize.py",                    "report_p20_style.txt"),
 ("p21_reporting_completeness_audit.py",       "report_p21_reporting_audit.txt"),
 ("p21f_layerB_supplement_reconcile.py",       "report_p21f_supplement_cells.txt"),
 ("p23_reference_merge_execute.py",            "report_p23_reference_merge.txt"),
 # ★ 2026-10-12：INC-035 的仪器 —— 与上面两条同属「引用」类，故并入本批。
 #   前两条查「号是否存在 / 条目是否合并」，**不查出现的先后**；本条补的正是那一格。
 ("p58_reference_order_audit.py",              "report_p58_reference_order.txt"),
 ("p59_survivor_scale_audit.py",               "report_p59_survivor_scale.txt"),
]

t0 = time.time()
rows = []
for script, report in AUDITS:
    sp = os.path.join(HERE, script)
    rp = os.path.join(HERE, report)
    t_start = time.time()
    r = subprocess.run([sys.executable, sp], capture_output=True, text=True,
                       encoding="utf-8", errors="replace", cwd=HERE)
    if not os.path.exists(rp):
        state = "MISSING"
    elif os.path.getmtime(rp) < t_start - 1:      # 留 1s 容差
        state = "STALE"                            # 报告不是本次写的 → 内容不可信
    elif r.returncode != 0:
        state = "FAIL"
    else:
        state = "PASS"
    rows.append((script, state, r.returncode, os.path.getsize(rp) if os.path.exists(rp) else 0,
                 (r.stdout or "")[-160:].replace("\n", " ")))

print("=" * 100)
print("run_all_audits —— 统一审计运行 + 报告新鲜度校验")
print("=" * 100)
for script, state, rc, size, tail in rows:
    mark = {"PASS": "  ", "FAIL": "!!", "STALE": "??", "MISSING": "XX"}[state]
    print(f"  {mark} {script:<44}{state:<8}rc={rc:<4}report={size:>7}B")

bad = [r for r in rows if r[1] != "PASS"]
print("")
print(f"  共 {len(rows)} 个审计器；非 PASS {len(bad)} 个")
for script, state, rc, _, tail in bad:
    print(f"    {state}: {script}")
    if tail: print(f"        …{tail[:150]}")
print("")
print("*** 全部通过且报告均为本次生成 ***" if not bad
      else "*** 存在非 PASS —— 注意 STALE 表示磁盘报告不是本次生成的，其内容不可引用 ***")
print(f"  （总耗时 {time.time()-t0:.0f}s）")
sys.exit(0 if not bad else 1)
