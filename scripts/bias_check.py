#!/usr/bin/env python3
"""Analyse results/bias/*: are the Horvitz-Thompson-type estimators unbiased in practice?

For each method (cs = CoStrata, csh = CoStrata-H) and query, the estimates of S independent seeds
are averaged; z = (mean - truth) / (sd / sqrt(S)). Under unbiasedness, z is approximately standard
normal for queries whose estimator is not too skewed (we restrict to queries with truth >= 1000).
"ht" is the unbiased-by-design estimate (column ht), "cal" the calibrated ratio estimate (column est).
Writes results/bias.json.
"""
import glob
import json
from collections import defaultdict

import numpy as np
import pandas as pd

out = {}
runs = defaultdict(list)
for f in glob.glob("results/bias/*.csv"):
    name = f.split("/")[-1]
    ds = name.split(".")[0]
    runs[ds].append(f)
for ds, files in sorted(runs.items()):
    w = pd.read_csv(f"work/{ds}.wl", sep=" ", header=None, usecols=[0, 1, 2, 3], names=["kind", "npred", "truth", "minfreq"])
    truth = w.truth.to_numpy(float)
    acc = defaultdict(lambda: dict(cal=[], ht=[]))
    mem = None
    for f in files:
        d = pd.read_csv(f, header=None, names=["ds", "method", "mem", "seed", "qid", "est", "ht", "sd", "cm"])
        mem = float(d.mem.iloc[0])
        for method, g in d.groupby("method"):
            g = g.sort_values("qid")
            if len(g) != len(w):
                continue
            acc[method]["cal"].append(g.est.to_numpy(float))
            acc[method]["ht"].append(g.ht.to_numpy(float))
    big = truth >= 1000
    for method, a in sorted(acc.items()):
        res = dict(seeds=len(a["ht"]), queries=int(big.sum()), mem=mem)
        for label in ("ht", "cal"):
            E = np.array(a[label])
            S = len(E)
            mean, sd = E.mean(axis=0), E.std(axis=0, ddof=1)
            ok = big & (sd > 0)
            z = (mean[ok] - truth[ok]) / (sd[ok] / np.sqrt(S))
            tot = E[:, big].sum(axis=1)
            res[label] = dict(
                mean_z=float(z.mean()), frac_absz_lt_1_96=float(np.mean(np.abs(z) < 1.96)),
                rel_bias_total=float(tot.mean() / truth[big].sum() - 1),
                rel_bias_total_se=float(tot.std(ddof=1) / np.sqrt(S) / truth[big].sum()),
                median_abs_rel_bias=float(np.median(np.abs(mean[big] / truth[big] - 1))),
                max_abs_rel_bias=float(np.max(np.abs(mean[big] / truth[big] - 1))),
                mean_rel_bias=float(np.mean(mean[big] / truth[big] - 1)))
        out[f"{ds}/{method}"] = res
        print(ds, method, json.dumps(res, indent=1))
json.dump(out, open("results/bias.json", "w"), indent=1)
