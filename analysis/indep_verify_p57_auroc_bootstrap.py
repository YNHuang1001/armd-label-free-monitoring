# -*- coding: utf-8 -*-
"""Independent second implementation: harmful-detection AUROC + 95% cluster-bootstrap CI
for the four post-hoc risk scores (ARMD paper 2, p57 post-hoc baselines).

Independence pledge: this file was written from the task brief and the raw CSV only;
p57_posthoc_baselines.py was NOT read.

Route deliberately different from a default implementation:
  - no sklearn / no pandas: CSV parsed with the stdlib csv module;
  - AUROC computed by hand from the Mann-Whitney rank statistic with average ranks (ties);
  - a fully independent O(n_pos*n_neg) pairwise concordance routine cross-checks every
    point estimate (ties counted 0.5);
  - bootstrap loop written from scratch: resample the 37 organism-drug clusters with
    replacement via random.Random(20261007), B=10000;
  - percentile interval computed by hand with linear interpolation (numpy "linear" style).

Definitions used exactly as fixed in the brief:
  cluster unit for uncertainty = (organism, drug)                      -> 37 clusters
  source model                 = (organism, drug, source_site)         -> 93 models
  harmful                      = d_auroc <= -0.05
  prespecified_DR   = -AUROC_DR
  src_only          = +auroc_source_internal
  amean_control     = auroc_source_internal - mean(AUROC_DR) within the same source model
  posthoc_dhat_DR   = -(AUROC_DR - auroc_source_internal)
"""

import csv
import math
import random
import sys
import os

DATA = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "figures", "data", "fig4_l3_scatter.csv")
B = 10000
SEED = 20261007
SCORES = ["prespecified_DR", "src_only", "amean_control", "posthoc_dhat_DR"]


def load_rows(path):
    rows = []
    with open(path, newline="", encoding="utf-8-sig") as fh:
        rdr = csv.DictReader(fh)
        for rec in rdr:
            rows.append(
                dict(
                    organism=rec["organism"].strip(),
                    drug=rec["drug"].strip(),
                    source_site=rec["source_site"].strip(),
                    auroc_source_internal=float(rec["auroc_source_internal"]),
                    auroc_target=float(rec["auroc_target"]),
                    d_auroc=float(rec["d_auroc"]),
                    AUROC_DR=float(rec["AUROC_DR"]),
                    d_hat_DR=float(rec["d_hat_DR"]),
                )
            )
    return rows


def auroc_rank(pairs):
    """AUC (label 1 = harmful = positive) via Mann-Whitney rank statistic.

    Average ranks for ties; returns None if a resample lacks one class.
    """
    n = len(pairs)
    npos = 0
    for _, lab in pairs:
        npos += lab
    nneg = n - npos
    if npos == 0 or nneg == 0:
        return None
    order = sorted(range(n), key=lambda i: pairs[i][0])
    ranks = [0.0] * n
    i = 0
    while i < n:
        j = i
        v = pairs[order[i]][0]
        while j + 1 < n and pairs[order[j + 1]][0] == v:
            j += 1
        avg_rank = (i + j) / 2.0 + 1.0  # 1-based average rank
        for k in range(i, j + 1):
            ranks[order[k]] = avg_rank
        i = j + 1
    rsum_pos = 0.0
    for idx in range(n):
        if pairs[idx][1] == 1:
            rsum_pos += ranks[idx]
    return (rsum_pos - npos * (npos + 1) / 2.0) / (npos * nneg)


def auroc_pairwise(pairs):
    """Independent O(npos*nneg) concordance check; ties count 0.5. Point estimates only."""
    pos = [s for s, lab in pairs if lab == 1]
    neg = [s for s, lab in pairs if lab == 0]
    if not pos or not neg:
        return None
    c = 0.0
    for p in pos:
        for q in neg:
            if p > q:
                c += 1.0
            elif p == q:
                c += 0.5
    return c / (len(pos) * len(neg))


def percentile_linear(sorted_vals, q):
    """Percentile with linear interpolation (numpy default), q in [0,1]."""
    n = len(sorted_vals)
    if n == 1:
        return sorted_vals[0]
    pos = q * (n - 1)
    lo = int(math.floor(pos))
    hi = int(math.ceil(pos))
    if lo == hi:
        return sorted_vals[lo]
    frac = pos - lo
    return sorted_vals[lo] + frac * (sorted_vals[hi] - sorted_vals[lo])


def main():
    rows = load_rows(DATA)
    n = len(rows)

    for r in rows:
        r["harmful"] = 1 if r["d_auroc"] <= -0.05 else 0
    n_harm = sum(r["harmful"] for r in rows)

    # --- checks on dataset shape -------------------------------------------
    clusters = sorted(set((r["organism"], r["drug"]) for r in rows))
    cl_idx = {c: i for i, c in enumerate(clusters)}
    m = len(clusters)
    for r in rows:
        r["cl"] = cl_idx[(r["organism"], r["drug"])]

    model_rows = {}
    for r in rows:
        model_rows.setdefault((r["organism"], r["drug"], r["source_site"]), []).append(r)
    n_models = len(model_rows)
    model_sizes = sorted({len(v) for v in model_rows.values()})

    print("n_rows=%d  n_harmful=%d  n_clusters=%d  n_source_models=%d"
          % (n, n_harm, m, n_models))
    print("source-model row counts (distinct sizes): %s" % (model_sizes,))

    # --- score construction (exactly per brief) ----------------------------
    for r in rows:
        key = (r["organism"], r["drug"], r["source_site"])
        grp = model_rows[key]
        mean_dr = sum(x["AUROC_DR"] for x in grp) / len(grp)
        r["prespecified_DR"] = -r["AUROC_DR"]
        r["src_only"] = r["auroc_source_internal"]
        r["amean_control"] = r["auroc_source_internal"] - mean_dr
        r["posthoc_dhat_DR"] = -(r["AUROC_DR"] - r["auroc_source_internal"])

    # awareness check: derived posthoc score vs the d_hat_DR column
    md = max(abs(r["d_hat_DR"] - r["posthoc_dhat_DR"]) for r in rows)
    print("max|d_hat_DR - (auroc_source_internal - AUROC_DR)| = %.3e" % md)

    # --- point estimates + cluster bootstrap --------------------------------
    out = {}
    for sname in SCORES:
        all_pairs = [(r[sname], r["harmful"]) for r in rows]
        auc = auroc_rank(all_pairs)
        auc_check = auroc_pairwise(all_pairs)
        ok = (auc is not None and auc_check is not None
              and abs(auc - auc_check) < 1e-12)

        cluster_pairs = [[] for _ in range(m)]
        for r in rows:
            cluster_pairs[r["cl"]].append((r[sname], r["harmful"]))

        rng = random.Random(SEED)
        boot = []
        skipped = 0
        for _ in range(B):
            pool = []
            for _ in range(m):
                pool.extend(cluster_pairs[rng.randrange(m)])
            a = auroc_rank(pool)
            if a is None:
                skipped += 1
            else:
                boot.append(a)
        boot.sort()
        lo = percentile_linear(boot, 0.025)
        hi = percentile_linear(boot, 0.975)
        bmean = sum(boot) / len(boot)
        out[sname] = (auc, lo, hi)

        print("%-18s AUROC=%.6f  rank==pairwise:%s  bootmean=%.6f  "
              "95%%CI=[%.6f, %.6f]  kept=%d skipped=%d"
              % (sname, auc, ok, bmean, lo, hi, len(boot), skipped))

    print("---- machine readable (4 dp) ----")
    for sname in SCORES:
        auc, lo, hi = out[sname]
        print("RESULT auroc_%s=%.4f auroc_%s_lo=%.4f auroc_%s_hi=%.4f"
              % (sname, auc, sname, lo, sname, hi))


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    main()
