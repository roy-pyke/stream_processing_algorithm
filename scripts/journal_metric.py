#!/usr/bin/env python3
"""Journal-style workload on the journal-style stream: which metric decides the winner?

Reads results/journal/jz13p9-vldbj.<mem>.<seed>.csv and work/jz13p9-vldbj.wl and reports, by number of
predicates k, the mean absolute error (in records), the mean relative error, the fraction of queries
answered with exactly 0, and the mean true answer. Writes results/journal_by_k.csv and
paper/tables/journal_metric.tex (10 MB).
"""
import glob
import os

import numpy as np
import pandas as pd

WL = "jz13p9-vldbj"
w = pd.DataFrame([list(map(int, l.split()[:4])) for l in open(f"work/{WL}.wl")], columns=["kind", "k", "truth", "minf"])
w["q"] = range(len(w))
rows = []
for f in sorted(glob.glob(f"results/journal/{WL}.*.csv")):
    d = pd.read_csv(f, header=None, names=["ds", "method", "mem", "seed", "q", "est", "ht", "sd", "cm"])
    m = d.merge(w, on="q")
    m["ae"] = (m.est - m.truth).abs()
    m["re"] = m.ae / m.truth
    m["zero"] = (m.est == 0).astype(float)
    g = m.groupby(["method", "mem", "k"]).agg(ae=("ae", "mean"), re=("re", "mean"), zero=("zero", "mean"), truth=("truth", "mean"),
                                              medre=("re", "median")).reset_index()
    g["seed"] = int(d.seed.iloc[0])
    rows.append(g)
a = pd.concat(rows).groupby(["method", "mem", "k"], as_index=False).mean().drop(columns="seed")
a.to_csv("results/journal_by_k.csv", index=False)
pd.set_option("display.width", 220)
for mem in (10.0, 100.0):
    s = a[a.mem == mem]
    print(f"\n=== {mem:g} MB: mean relative error by k")
    print(s.pivot(index="method", columns="k", values="re").to_string(float_format=lambda v: "%.3g" % v))
    print(f"=== {mem:g} MB: mean absolute error (records) by k")
    print(s.pivot(index="method", columns="k", values="ae").to_string(float_format=lambda v: "%.3g" % v))
    print(f"=== {mem:g} MB: fraction of queries answered with 0")
    print(s.pivot(index="method", columns="k", values="zero").to_string(float_format=lambda v: "%.2f" % v))
print("\nmean truth by k:", a[a.method == "res"].groupby("k").truth.mean().round(1).to_dict())

import json
s10 = a[a.mem == 10.0].set_index(["method", "k"])
low = [1, 2, 3, 4]
high = [k for k in sorted(a.k.unique()) if k >= 6]
ratio = lambda m, col: [float(s10.loc[("omniv", k), col] / s10.loc[(m, k), col]) for k in low]
num = dict(
    ash_abs_lo=min(ratio("ash", "ae")), ash_abs_hi=max(ratio("ash", "ae")),
    res_abs_lo=min(ratio("res", "ae")), res_abs_hi=max(ratio("res", "ae")),
    omni_rel_high=float(np.mean([s10.loc[("omniv", k), "re"] for k in high])),
    ash_rel_high=float(np.mean([s10.loc[("ash", k), "re"] for k in high])),
    res_rel_high=float(np.mean([s10.loc[("res", k), "re"] for k in high])),
    omni_zero_high=float(np.mean([s10.loc[("omniv", k), "zero"] for k in high])),
    truth_k1=float(s10.loc[("res", 1), "truth"]), truth_high=float(np.mean([s10.loc[("res", k), "truth"] for k in high])),
    share_high=float((w.k >= 5).mean()))
json.dump(num, open("results/journal_numbers.json", "w"), indent=1)

lab = {"omniv": "OmniSketch (journal cfg.)", "ash": "aSH (journal acct.)", "resj": "RS (journal acct.)", "res": "RS",
       "cs:0:1e9:1:0.25": "RS+C2"}
os.makedirs("paper/tables", exist_ok=True)
s = a[a.mem == 10.0]
ks = [1, 2, 3, 4, 6, 9]
with open("paper/tables/journal_metric.tex", "w") as f:
    f.write("\\begin{tabular}{l" + "r" * (2 * len(ks)) + "}\n\\toprule\n")
    f.write(f" & \\multicolumn{{{len(ks)}}}{{c}}{{mean absolute error (records)}} & \\multicolumn{{{len(ks)}}}{{c}}{{mean relative error}} \\\\\n")
    f.write("$k$ & " + " & ".join(str(k) for k in ks) + " & " + " & ".join(str(k) for k in ks) + " \\\\\n\\midrule\n")
    t = s[s.method == "res"].set_index("k").truth
    num_ = lambda v: (f"{v:,.0f}".replace(",", "\\,") if v >= 100 else f"{v:.3g}")
    f.write("\\emph{mean true answer} & " + " & ".join(num_(t[k]) for k in ks) + " & " + " & ".join("" for _ in ks) + " \\\\\n\\midrule\n")
    for key, name in lab.items():
        r = s[s.method == key].set_index("k")
        if not len(r):
            continue
        f.write(name + " & " + " & ".join(num_(r.ae[k]) for k in ks) + " & " + " & ".join(f"{r.re[k]:.2f}" for k in ks) + " \\\\\n")
    f.write("\\bottomrule\n\\end{tabular}\n")
