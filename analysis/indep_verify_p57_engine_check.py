# -*- coding: utf-8 -*-
"""Engine/stream sensitivity and boundary-stability diagnostics (not for reporting).

Reuses functions from indep_verify_p57_auroc_bootstrap.py (my own file).
"""
import importlib.util
import random

import numpy as np
import os

MAIN = os.path.join(os.path.dirname(os.path.abspath(__file__)), "indep_verify_p57_auroc_bootstrap.py")
spec = importlib.util.spec_from_file_location("indepv", MAIN)
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)

rows = mod.load_rows(mod.DATA)
for r in rows:
    r["harmful"] = 1 if r["d_auroc"] <= -0.05 else 0
    r["prespecified_DR"] = -r["AUROC_DR"]
    r["src_only"] = r["auroc_source_internal"]

print("closest |d_auroc + 0.05| among rows:",
      ["%.2e" % m for m in sorted(abs(r["d_auroc"] + 0.05) for r in rows)[:3]])

clusters = sorted(set((r["organism"], r["drug"]) for r in rows))
m = len(clusters)
cl_idx = {c: i for i, c in enumerate(clusters)}
for r in rows:
    r["cl"] = cl_idx[(r["organism"], r["drug"])]


def cp_for(sname):
    cp = [[] for _ in range(m)]
    for r in rows:
        cp[r["cl"]].append((r[sname], r["harmful"]))
    return cp


def pct(boot):
    boot.sort()
    return mod.percentile_linear(boot, 0.025), mod.percentile_linear(boot, 0.975)


def run_cp(cp, B, seed):
    rng = random.Random(seed)
    boot = []
    for _ in range(B):
        pool = []
        for _ in range(m):
            pool.extend(cp[rng.randrange(m)])
        a = mod.auroc_rank(pool)
        if a is not None:
            boot.append(a)
    return pct(boot)


cp = cp_for("prespecified_DR")

# stream check: default_rng batch vs per-replicate generation
r1 = np.random.default_rng(20261007)
d1 = r1.integers(0, m, size=(10000, m))
r2 = np.random.default_rng(20261007)
d2 = np.array([r2.integers(0, m, size=m) for _ in range(10000)])
print("default_rng batch == per-rep draws:", bool((d1 == d2).all()))

# legacy MT19937 batch draws
np.random.seed(20261007)
draws = np.random.randint(0, m, size=(10000, m))
boot = []
for draw in draws:
    pool = []
    for k in draw:
        pool.extend(cp[int(k)])
    a = mod.auroc_rank(pool)
    if a is not None:
        boot.append(a)
print("legacy RandomState(20261007) batch  prespecified_DR: [%.4f, %.4f]" % pct(boot))

# python random.choices stream
rng = random.Random(20261007)
boot = []
for _ in range(10000):
    pool = []
    for k in rng.choices(range(m), k=m):
        pool.extend(cp[k])
    a = mod.auroc_rank(pool)
    if a is not None:
        boot.append(a)
print("python random.choices(20261007)     prespecified_DR: [%.4f, %.4f]" % pct(boot))

# src_only lower-endpoint stability across 12 seeds at B=10k
cp2 = cp_for("src_only")
los, his = [], []
for s in range(20261007, 20261007 + 12):
    lo, hi = run_cp(cp2, 10000, s)
    los.append(lo)
    his.append(hi)
print("src_only 12 seeds B=10k: lo range [%.4f, %.4f]  hi range [%.4f, %.4f]"
      % (min(los), max(los), min(his), max(his)))
print("src_only lo values:", " ".join("%.4f" % v for v in los))
