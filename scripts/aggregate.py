#!/usr/bin/env python3
"""Aggregate per-query benchmark output into per-run metrics.

  aggregate.py <results_dir> <out.csv> [workload_dir]

For every (dataset, method, mem, seed) run we compute, per workload family:
  nae      mean |est - truth| / N            (the metric of the OmniSketch paper)
  mre      mean |est - truth| / truth
  mre_d{j} mean relative error over queries whose rarest predicate has frequency in [10^j, 10^(j+1))
  mre_k{k} mean relative error over queries with k predicates
  cov95    (CoStrata only) fraction of queries with |est - truth| <= 1.96 * sd
"""
import glob
import os
import sys

import numpy as np
import pandas as pd

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
    df = pd.read_csv(base + ".csv", header=None, names=["ds", "method", "mem", "seed", "qid", "est", "ht", "sd", "cm"])
    wname = os.path.basename(base).split(".")[0]
    # p-scaling runs are named <dataset>_p<k>.<mem>.<seed>
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
        for kind in (0, 1, 2, 3):
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
        if g.ht.notna().all():
            ht = g.ht.to_numpy(float)
            sd = g.sd.to_numpy(float)
            row["nae0_ht"] = np.abs(ht - truth)[(w.kind == 0).to_numpy()].mean() / N
            row["nae1_ht"] = np.abs(ht - truth)[(w.kind == 1).to_numpy()].mean() / N
            row["cov95"] = np.mean(err <= 1.96 * sd)
            big = truth >= 100
            row["cov95_big"] = np.mean(err[big] <= 1.96 * sd[big])
            row["bias_ht"] = ht.sum() / truth.sum() - 1.0
        rows.append(row)

pd.DataFrame(rows).to_csv(out, index=False)
print(f"wrote {out}: {len(rows)} runs")
