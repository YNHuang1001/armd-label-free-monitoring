# -*- coding: utf-8 -*-
"""MC-noise / convergence checks for the independent AUROC cluster bootstrap.

Reuses functions from indep_verify_p57_auroc_bootstrap.py (my own file, written from
the brief; p57_posthoc_baselines.py never read). Not for reporting: diagnostics that
(a) quantify the run-to-run Monte-Carlo spread at B=10000 per score,
(b) give a high-precision converged interval (B=200000) per score as a reference,
(c) verify the exact algebraic identities between CSV columns and the brief's formulas.
"""
import importlib.util
import random
import time
import os

try:
    import numpy as np
    HAVE_NP = True
except Exception:
    HAVE_NP = False

MAIN = os.path.join(os.path.dirname(os.path.abspath(__file__)), "indep_verify_p57_auroc_bootstrap.py")
spec = importlib.util.spec_from_file_location("indepv", MAIN)
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)

rows = mod.load_rows(mod.DATA)
for r in rows:
    r["harmful"] = 1 if r["d_auroc"] <= -0.05 else 0

# exact-identity checks between columns and brief formulas
mx1 = max(abs(r["d_hat_DR"] - (r["AUROC_DR"] - r["auroc_source_internal"])) for r in rows)
mx2 = max(abs(r["d_auroc"] - (r["auroc_target"] - r["auroc_source_internal"])) for r in rows)
print("identity |d_hat_DR - (AUROC_DR - src_internal)| max = %.3e" % mx1)
print("identity |d_auroc - (auroc_target - src_internal)| max = %.3e" % mx2)

# score construction, exactly per brief
model_rows = {}
for r in rows:
    model_rows.setdefault((r["organism"], r["drug"], r["source_site"]), []).append(r)
for r in rows:
    grp = model_rows[(r["organism"], r["drug"], r["source_site"])]
    mean_dr = sum(x["AUROC_DR"] for x in grp) / len(grp)
    r["prespecified_DR"] = -r["AUROC_DR"]
    r["src_only"] = r["auroc_source_internal"]
    r["amean_control"] = r["auroc_source_internal"] - mean_dr
    r["posthoc_dhat_DR"] = -(r["AUROC_DR"] - r["auroc_source_internal"])

clusters = sorted(set((r["organism"], r["drug"]) for r in rows))
m = len(clusters)
cl_idx = {c: i for i, c in enumerate(clusters)}
for r in rows:
    r["cl"] = cl_idx[(r["organism"], r["drug"])]
sizes = sorted(len([1 for r in rows if r["cl"] == k]) for k in range(m))
print("m=%d  cluster sizes=%s  total rows=%d" % (m, sizes, sum(sizes)))


def cluster_pairs_for(sname):
    cp = [[] for _ in range(m)]
    for r in rows:
        cp[r["cl"]].append((r[sname], r["harmful"]))
    return cp


def run(cp, B, seed):
    rng = random.Random(seed)
    boot = []
    for _ in range(B):
        pool = []
        for _ in range(m):
            pool.extend(cp[rng.randrange(m)])
        a = mod.auroc_rank(pool)
        if a is not None:
            boot.append(a)
    boot.sort()
    return mod.percentile_linear(boot, 0.025), mod.percentile_linear(boot, 0.975), len(boot)


t0 = time.perf_counter()
for sname in mod.SCORES:
    cp = cluster_pairs_for(sname)
    val = mod.auroc_rank([(r[sname], r["harmful"]) for r in rows])
    lo0, hi0, _ = run(cp, 10000, 20261007)  # must equal the reported run (deterministic)
    los, his = [], []
    for seed in (20261007, 20261008, 20261009):
        lo, hi, kept = run(cp, 10000, seed)
        los.append(lo)
        his.append(hi)
    ref_lo, ref_hi, _ = run(cp, 200000, 7)  # high-precision bootstrap-percentile reference

    def sd_mean(v):
        mu = sum(v) / len(v)
        sd = (sum((x - mu) ** 2 for x in v) / (len(v) - 1)) ** 0.5
        return sd, mu

    slo, mlo = sd_mean(los)
    shi, mhi = sd_mean(his)
    print("%-18s point=%.6f  seed20261007 B=10k [%.4f, %.4f]" % (sname, val, lo0, hi0))
    print("   3 seeds B=10k: lo mean=%.4f sd=%.4f | hi mean=%.4f sd=%.4f | "
          "reference B=200k [%.4f, %.4f]  (%.0fs)"
          % (mlo, slo, mhi, shi, ref_lo, ref_hi, time.perf_counter() - t0), flush=True)

# diagnostic: same seed through numpy default_rng batch draws (engine comparison only)
if HAVE_NP:
    cp = cluster_pairs_for("prespecified_DR")
    rng = np.random.default_rng(20261007)
    draws = rng.integers(0, m, size=(10000, m))
    boot = []
    for draw in draws:
        pool = []
        for k in draw:
            pool.extend(cp[int(k)])
        a = mod.auroc_rank(pool)
        if a is not None:
            boot.append(a)
    boot.sort()
    print("numpy default_rng(20261007) B=10k, prespecified_DR: [%.4f, %.4f]"
          % (mod.percentile_linear(boot, 0.025), mod.percentile_linear(boot, 0.975)))
else:
    print("numpy not available; engine-comparison diagnostic skipped")
