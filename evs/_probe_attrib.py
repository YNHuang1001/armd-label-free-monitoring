# -*- coding: utf-8 -*-
"""按 sim_fix.py 里**真实的**函数源码切片执行，验证新归因逻辑。不重写逻辑。"""
import io, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
t = io.open("sim_fix.py", encoding="utf-8").read()
i = t.index("def _manifest_registry():")
j = t.index("def family_summary(")
ns = {}
exec(compile(t[i:j], "<sim_fix slice>", "exec"), ns)
print("registry loaded: kind=%d  prefix=%d" % (len(ns["_MAN_KIND"]), len(ns["_MAN_PFX"])))
print()
cases = ["SUDO-01","SUDO-03","SUDO-05","FORK-03","BOUND-01-AGG",
         "BOUND-01::figS5_mmd_summary.csv::mmd_min::R07_mmd_floor",
         "EVAL-B1","EC-01","EC-03","MMD-02b","SENS-08b","UNIT-07",
         "L3-02","MC-10","NO-SUCH-TEST"]
print("%-22s %-24s %s" % ("test_id", "char_pin_reason?", "patched_fail_reason"))
print("-"*100)
for c in cases:
    cp = ns["char_pin_reason"](c)
    v, w = ns["patched_fail_reason"](c)
    print("%-22s %-24s %s" % (c, "yes" if cp else "no", v))
    if v == "expected":
        print("%-22s %s" % ("", w[:150]))
