#!/usr/bin/env python3
"""Aggregate per-query benchmark output into per-run metrics.

  aggregate.py <results_dir> <out.csv> [workload_dir]

For every (dataset, method, mem, seed) run we compute, per workload family:
  nae      mean |est - truth| / N            (the metric of the OmniSketch paper)
  mre      mean |est - truth| / truth
  mre_d{j} mean relative error over queries whose rarest predicate has frequency in [10^j, 10^(j+1))
  mre_k{k} mean relative error over queries with k predicates
  cov95    (CoStrata only) fraction of queries with |est - truth| <= 1.96 * sd
  ci95     (CoStrata, CoStrata-H; runs with interval diagnostics) coverage of the small-count
           intervals: rule of three [0, min(n_c, ln(20) / pi_q)] when no sampled row matches, otherwise
           the Fay-Feuer gamma interval around the calibrated estimate, truncated to [m, n_c]
           (m = matching sampled rows, n_c = anchor count, both hard bounds on the answer)
"""
import glob
import os
import sys

import re

import numpy as np
import pandas as pd

try:
    from scipy.stats import gamma as _gamma
except ImportError:  # the new intervals need scipy; everything else does not
    _gamma = None

ALPHA = 0.05


def small_count_ci(est, ht, sd, m, wq, nb):
    """95% intervals for calibrated estimates (see the module docstring); returns lo, hi arrays."""
    est, ht, sd, m, wq, nb = (np.asarray(x, float) for x in (est, ht, sd, m, wq, nb))
    lo = np.zeros_like(est)
    hi = np.minimum(nb, np.log(1 / ALPHA) * wq)              # no matching row: P(m = 0) <= exp(-pi_q f)
    pos = m > 0
    y, v = est[pos], np.maximum(sd[pos] ** 2, 1e-12)
    scale = np.where(ht[pos] > 0, est[pos] / np.maximum(ht[pos], 1e-300), 1.0)
    wm = scale * wq[pos]                                     # calibrated weight of one more matching record
    with np.errstate(divide="ignore", invalid="ignore"):
        l_ = _gamma.ppf(ALPHA / 2, a=y * y / v, scale=v / y)
        h_ = _gamma.ppf(1 - ALPHA / 2, a=(y + wm) ** 2 / (v + wm * wm), scale=(v + wm * wm) / (y + wm))
    l_ = np.where(np.isfinite(l_), l_, y)
    lo[pos] = np.maximum(l_, m[pos])
    hi[pos] = np.minimum(np.where(np.isfinite(h_), h_, nb[pos]), nb[pos])
    return lo, np.maximum(hi, lo)

res_dir, out = sys.argv[1], sys.argv[2]
wl_dir = sys.argv[3] if len(sys.argv) > 3 else "work"

workloads = {}


def workload(name):
    if name not in workloads:
        w = pd.read_csv(os.path.join(wl_dir, name + ".wl"), sep=" ", header=None, usecols=[0, 1, 2, 3],
                        names=["kind", "npred", "truth", "minfreq"], engine="c")
        workloads[name] = w
    return workloads[name]


rows = []
for log in sorted(glob.glob(os.path.join(res_dir, "*.log"))):
    base = log[:-4]
    if not os.path.exists(base + ".csv"):
        continue
    runs = {}
    for line in open(log):
        if line.startswith("RUN,"):
            f = line.rstrip("\n").split(",", 11)
            runs[f[2]] = dict(p=int(f[5]), n=int(f[6]), ns_per_rec=float(f[8]), query_us=float(f[9]),
                              used_mb=float(f[10]), info=f[11])
    df = pd.read_csv(base + ".csv", header=None, names=["ds", "method", "mem", "seed", "qid", "est", "ht", "sd", "cm", "m", "wq"])
    wname = os.path.basename(base).split(".")[0]
    # p-scaling runs are named <dataset>_p<k>.<mem>.<seed>, merge runs <dataset>_m<parts>.<mem>.<seed>
    wname = re.sub(r"_m\d+$", "", wname)
    w = workload(wname)
    truth = w.truth.to_numpy(float)
    dec = np.floor(np.log10(np.maximum(w.minfreq.to_numpy(float), 1))).astype(int)
    for method, g in df.groupby("method", sort=False):
        if method not in runs:
            continue
        g = g.sort_values("qid")
        if len(g) != len(w):
            continue
        N = runs[method]["n"]
        est = g.est.to_numpy(float)
        err = np.abs(est - truth)
        rel = err / truth
        row = dict(dataset=g.ds.iloc[0], method=method, mem=float(g.mem.iloc[0]), seed=int(g.seed.iloc[0]), **runs[method])
        for kind in (0, 1, 2, 3, 4):
            sel = (w.kind == kind).to_numpy()
            if sel.sum() == 0:
                continue
            row[f"nae{kind}"] = err[sel].mean() / N
            row[f"mre{kind}"] = rel[sel].mean()
            row[f"medre{kind}"] = np.median(rel[sel])
            for j in range(1, 7):
                s2 = sel & (dec == j)
                row[f"mre{kind}_d{j}"] = rel[s2].mean() if s2.sum() >= 5 else np.nan
                row[f"cnt{kind}_d{j}"] = int(s2.sum())
            for k in range(1, 7):
                s2 = sel & (w.npred == k).to_numpy()
                row[f"mre{kind}_k{k}"] = rel[s2].mean() if s2.sum() >= 5 else np.nan
                row[f"nae{kind}_k{k}"] = err[s2].mean() / N if s2.sum() >= 5 else np.nan
        row["nae_all"] = err.mean() / N
        if method.startswith("cs") and g.cm.notna().all():
            # counters only: the Count-Min bound of the anchor cell, without the sample
            cme = np.abs(g.cm.to_numpy(float) - truth)
            for k in range(1, 7):
                s2 = ((w.kind == 0) & (w.npred == k)).to_numpy()
                row[f"cmnae0_k{k}"] = cme[s2].mean() / N if s2.sum() >= 5 else np.nan
        if g.ht.notna().all():
            ht = g.ht.to_numpy(float)
            sd = g.sd.to_numpy(float)
            for kind in (0, 1, 4):
                sel = (w.kind == kind).to_numpy()
                if sel.sum():
                    row[f"nae{kind}_ht"] = np.abs(ht - truth)[sel].mean() / N
                    row[f"mre{kind}_ht"] = (np.abs(ht - truth) / truth)[sel].mean()
            row["cov95"] = np.mean(err <= 1.96 * sd)
            big = truth >= 100
            row["cov95_big"] = np.mean(err[big] <= 1.96 * sd[big])
            row["bias_ht"] = ht.sum() / truth.sum() - 1.0
            if g.m.notna().all() and _gamma is not None:
                lo, hi = small_count_ci(est, ht, sd, g.m, g.wq, g.cm)
                inside = (lo <= truth) & (truth <= hi)
                zero = g.m.to_numpy(float) == 0
                row["ci95"] = inside.mean()
                row["ci95_big"] = inside[big].mean()
                row["ci95_small"] = inside[truth < 100].mean() if (truth < 100).any() else np.nan
                row["cov95_small"] = np.mean((err <= 1.96 * sd)[truth < 100]) if (truth < 100).any() else np.nan
                row["frac_zero"] = zero.mean()
                row["ci95_zero"] = inside[zero].mean() if zero.any() else np.nan
                row["ci_relwidth"] = np.median((hi - lo) / truth)
        rows.append(row)

pd.DataFrame(rows).to_csv(out, index=False)
print(f"wrote {out}: {len(rows)} runs")
