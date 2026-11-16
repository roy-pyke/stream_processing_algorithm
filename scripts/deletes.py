#!/usr/bin/env python3
"""Aggregate the bounded-deletion runs (scripts/run_deletes.sh, #25) into results/agg_deletes.csv.

  deletes.py [results/deletes] [work/delres.wl] [results/agg_deletes.csv]

One row per (perc, method, mem, seed); N = the live records at the end of the stream (the residue).
  nae0m    mean over k = 2..6 of the NAE (|est - truth| / N) of the record-sampled queries with k
           predicates (family 0; the primary metric of the main evaluation)
  nae1     NAE of family 1
  naej     NAE of the journal-style queries (k = 1..9 from sampled residue records)
  mrej     their mean relative error; mrej_lo / mrej_hi: queries with k <= 4 / k >= 6
  ns_upd   ingestion time per update (ns), mem_mb the memory the synopsis reports, held (ours) rows
           in the store at the end
Methods: s2, ash (the journal code), res = cs:0:1e9, rsc2 = cs:0:1e9:1:0.25, rsc3 = cs:0:1e9:1:0.15:0.15.
"""
import glob
import os
import re
import sys

import numpy as np
import pandas as pd

src = sys.argv[1] if len(sys.argv) > 1 else "results/deletes"
wlp = sys.argv[2] if len(sys.argv) > 2 else "work/delres.wl"
dst = sys.argv[3] if len(sys.argv) > 3 else "results/agg_deletes.csv"
SPEC = {"cs:0:1e9": "res", "cs:0:1e9:1:0.25": "rsc2", "cs:0:1e9:1:0.15:0.15": "rsc3"}

w = pd.DataFrame([l.split()[:4] for l in open(wlp)], columns=["kind", "npred", "truth", "minfreq"]).astype(np.int64)
truth = w.truth.to_numpy(float)
rows = []
for f in sorted(glob.glob(os.path.join(src, "*.csv"))):
    parts = os.path.basename(f)[:-4].split(".")  # <perc>.<mem>.<seed>.<method>, perc may contain a dot
    perc, mem, seed, which = ".".join(parts[:-3]), parts[-3], parts[-2], parts[-1]
    log = open(f[:-4] + ".log").read() if os.path.exists(f[:-4] + ".log") else ""
    if os.path.getsize(f) == 0:
        continue
    if which == "ours":
        d = pd.read_csv(f, header=None, usecols=[1, 4, 5], names=["spec", "q", "est"])
        runs = {SPEC.get(s, s): g.sort_values("q").est.to_numpy(float) for s, g in d.groupby("spec")}
        info = {SPEC.get(m.group(1), m.group(1)): m for m in re.finditer(r"^DRUN,[^,]*,([^,]*),[^,]*,[^,]*,[^,]*,(\d+),[^,]*,([\d.]+),[^,]*,([\d.]+),(.*)$", log, re.M)}
    else:
        d = pd.read_csv(f, header=None, names=["name", "mem", "rep", "q", "truth", "est"])
        assert np.array_equal(d.sort_values("q").truth.to_numpy(float), truth), f
        runs = {which: d.sort_values("q").est.to_numpy(float)}
        info = {which: m for m in re.finditer(r"^JDRUN,[^,]*,[^,]*,[^,]*,[^,]*,[^,]*,(\d+),(\d+),[^,]*,([\d.]+),[^,]*,([\d.]+),(.*)$", log, re.M)}
    for meth, est in runs.items():
        if len(est) != len(truth):
            print("incomplete:", f, meth)
            continue
        m = info.get(meth)
        N = float(m.group(2)) if m and which != "ours" else 4194304.0
        if m and which == "ours":
            N = float(re.search(r";live=(\d+)", m.group(5)).group(1))
        err = np.abs(est - truth)
        r = {"perc": float(perc), "method": meth, "mem": float(mem), "seed": int(seed), "N": N}
        f0 = (w.kind == 0).to_numpy()
        r["nae0m"] = np.mean([err[f0 & (w.npred == k).to_numpy()].mean() / N for k in range(2, 7) if (f0 & (w.npred == k).to_numpy()).sum() >= 5])
        f1 = (w.kind == 1).to_numpy()
        r["nae1"] = err[f1].mean() / N if f1.any() else np.nan
        fj = (w.kind == 3).to_numpy()
        rel = err / np.maximum(truth, 1)
        r["naej"] = err[fj].mean() / N
        r["mrej"] = rel[fj].mean()
        r["mrej_lo"] = rel[fj & (w.npred <= 4).to_numpy()].mean()
        r["mrej_hi"] = rel[fj & (w.npred >= 6).to_numpy()].mean()
        if m:
            r["ns_upd"] = float(m.group(3))
            r["mem_mb"] = float(m.group(4))
            h = re.search(r";held=(\d+)", m.group(5))
            if h:
                r["held"] = int(h.group(1))
                r["M"] = int(re.search(r";M=(\d+)", m.group(5)[m.group(5).find("held"):]).group(1))
        rows.append(r)
out = pd.DataFrame(rows).sort_values(["perc", "mem", "method", "seed"])
out.to_csv(dst, index=False)
print(f"{len(out)} runs -> {dst}")
print(out.groupby(["perc", "mem", "method"])[["nae0m", "naej", "mrej"]].mean().unstack("method").to_string(float_format=lambda x: f"{x:.3g}"))
