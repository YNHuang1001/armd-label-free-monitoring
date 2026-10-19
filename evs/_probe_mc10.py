# -*- coding: utf-8 -*-
"""按 t14 的 MC-10 逐字逻辑，用**本轮**CSV 算出下轮 MC-10 的判定。只读。"""
import io, os, csv as _csv, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import evs_load as L

HERE = os.path.dirname(os.path.abspath(__file__))
_doc = L.load_manifest()
PATHS = _doc["paths"]

rows = []
for mode in ("real", "patched"):
    q = os.path.join(HERE, "evs_results_%s.csv" % mode)
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
    if x in req or x in ch or x in dg or x in ub:
        return True
    for pat in list(dg) + list(ub):
        if pat.endswith("*") and x.startswith(pat[:-1]):
            return True
    return False

bad = sorted({"%s:%s" % (m, x) for m, x in rows if not _known(x)})
allf = sorted({x for _, x in rows})
print("两个模式 FAIL test_id 共 %d 条；未归类 %d 条" % (len(allf), len(bad)))
for b in bad:
    print("   !! 未归类：", b)
print()
print("FAIL 清单：", allf)
