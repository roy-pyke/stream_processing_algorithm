#!/usr/bin/env python3
"""Quick textual summary of a bench CSV: summarize.py <workload.txt> <results.csv> <N> [col]

col = 5 (default estimate) or 6 (Horvitz-Thompson estimate of CoStrata).
Relative error is reported per decade of the rarest-predicate frequency ("anchor").
"""
import sys
from collections import defaultdict

import numpy as np

wl, res, N = sys.argv[1], sys.argv[2], float(sys.argv[3])
col = int(sys.argv[4]) if len(sys.argv) > 4 else 5
kind, npred, truth, minf = [], [], [], []
for line in open(wl):
    f = line.split()
    kind.append(int(f[0])); npred.append(int(f[1])); truth.append(float(f[2])); minf.append(float(f[3]))
kind, npred, truth, minf = map(np.array, (kind, npred, truth, minf))
est = defaultdict(lambda: defaultdict(lambda: np.full(len(truth), np.nan)))
for line in open(res):
    f = line.rstrip("\n").split(",")
    if col != 5 and f[col] == "":
        continue
    est[(f[1], f[2])][f[3]][int(f[4])] = float(f[col])

dec = np.floor(np.log10(np.maximum(minf, 1))).astype(int)
print(f"{'method':18s} {'mem':>5s} | {'NAE(k0)':>9s} {'NAE(k1)':>9s} | mean rel.err by anchor-frequency decade 1..6, kind 1 | by #pred 1..6, kind 0")
for (m, mem), seeds in sorted(est.items(), key=lambda x: (float(x[0][1]), x[0][0])):
    E = np.stack(list(seeds.values()))
    abs_err = np.abs(E - truth)
    nae0 = abs_err[:, kind == 0].mean() / N
    nae1 = abs_err[:, kind == 1].mean() / N
    rel = abs_err / truth
    a = []
    for dd in range(1, 7):
        sel = (kind == 1) & (dec == dd)
        a.append(f"{rel[:, sel].mean():6.2f}" if sel.sum() >= 5 else "     -")
    b = []
    for k in range(1, 7):
        sel = (kind == 0) & (npred == k)
        b.append(f"{rel[:, sel].mean():6.2f}" if sel.sum() >= 5 else "     -")
    print(f"{m:18s} {mem:>5s} | {nae0:9.2e} {nae1:9.2e} | {' '.join(a)} | {' '.join(b)}")
