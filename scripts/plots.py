#!/usr/bin/env python3
"""Produce all figures (PDF), LaTeX table fragments and quoted numbers of the paper.

  plots.py   reads results/agg_*.csv, results/time/*.log, results/bias.json and
             results/validation/numbers.json; writes paper/figs, paper/tables, paper/numbers.tex
"""
import json
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

os.makedirs("paper/figs", exist_ok=True)
os.makedirs("paper/tables", exist_ok=True)

RSC1, RSC2, RSC3 = "cs:0:1e9", "cs:0:1e9:1:0.25", "cs:0:1e9:1:0.15:0.15"
CS1, CS2, CSH = "cs", "cs:0:1:1:0.25", "csh:0:1:0.25"
USSC1, USSC2 = "ussc", "ussc:0:0.25"
# Fixed entity -> (label, colour, marker, linestyle). Colour follows the method in every figure.
STYLE = {
    "omnijs":   ("OmniSketch (authors' code)", "#e87ba4", "X", "--"),
    "omnip":    ("OmniSketch (ours, same $B$)", "#b05c7d", "x", "--"),
    "omniv":    ("OmniSketch (journal cfg.)",  "#9a6fb0", "+", "--"),
    "omnibest": ("OmniSketch (best)",          "#008300", "D", "--"),
    "ash":      ("aSH (journal acct.)",        "#8c8c84", "<", ":"),
    "ashl":     ("aSH (lean acct.)",           "#5c5c55", ">", ":"),
    "resj":     ("RS (journal acct.)",         "#c9a07a", "s", ":"),
    "res":      ("RS",                         "#eb6834", "s", "-"),
    "uss":      ("USS",                        "#1baf7a", "^", "-"),
    RSC1:       ("RS+C1",                      "#b0a4f0", "P", ":"),
    RSC2:       ("RS+C2",                      "#6250d6", "P", "-"),
    RSC3:       ("RS+C3",                      "#3a2c96", "*", "-"),
    USSC1:      ("USS+C1",                     "#f2c14e", "v", ":"),
    USSC2:      ("USS+C2",                     "#eda100", "v", "-"),
    CS1:        ("CoStrata (C1)",              "#7fb4ec", "o", ":"),
    CS2:        ("CoStrata",                   "#2a78d6", "o", "-"),
    CSH:        ("CoStrata-H",                 "#12a0b0", "h", "-"),
}
OMNI_GRID = ["omnil", "omni:28:3", "omni:64:2", "omni:256:2", "omni:1024:2", "omni:256:3", "omni:1024:1", "omni:4096:1", "omni:4096:2"]
OMNI_FAMILY = ["omnijs", "omnip", "omniv"] + OMNI_GRID
DS_LABEL = {"zipf10_5M": r"Zipf $\alpha{=}1.0$", "zipf13_5M": r"Zipf $\alpha{=}1.3$", "zipf15_5M": r"Zipf $\alpha{=}1.5$",
            "unif_5M": "Uniform", "taxi": "Taxi", "mawi": "MAWI", "dmv": "DMV"}
DATASETS = ["zipf10_5M", "zipf13_5M", "zipf15_5M", "unif_5M", "taxi", "mawi", "dmv"]
INK, MUTED, GRID = "#1a1a19", "#6b6a63", "#e4e2da"

plt.rcParams.update({
    "font.size": 8, "axes.titlesize": 8.5, "axes.labelsize": 8, "legend.fontsize": 7.5,
    "xtick.labelsize": 7.5, "ytick.labelsize": 7.5, "axes.edgecolor": MUTED, "axes.labelcolor": INK,
    "xtick.color": MUTED, "ytick.color": MUTED, "text.color": INK, "axes.linewidth": 0.6,
    "grid.color": GRID, "grid.linewidth": 0.6, "axes.spines.top": False, "axes.spines.right": False,
    "pdf.fonttype": 42, "font.family": "DejaVu Sans",
})

NUM = {}
NOTES = []  # free-form tables for results/summary.txt (not used by LaTeX)


def note(title, obj):
    NOTES.append(f"\n==== {title} ====\n{obj if isinstance(obj, str) else obj.to_string()}\n")


def mean_over_seeds(df, cols):
    cols = [c for c in cols if c in df.columns]
    return df.groupby(["dataset", "method", "mem"], as_index=False)[cols].mean()


def with_omnibest(m, cols, family=OMNI_FAMILY):
    """Add the pseudo-method 'omnibest': per (dataset, mem) and per metric, the best OmniSketch variant."""
    cols = [c for c in cols if c in m.columns]
    g = m[m.method.isin(family)]
    best = g.groupby(["dataset", "mem"], as_index=False)[cols].min()
    best["method"] = "omnibest"
    return pd.concat([m, best], ignore_index=True)


def line(ax, x, y, key, **kw):
    lab, col, mk, ls = STYLE[key]
    ax.plot(x, y, color=col, marker=mk, linestyle=ls, linewidth=1.6, markersize=4.5,
            markeredgecolor="white", markeredgewidth=0.6, label=lab, **kw)


def fmt_sci(x):
    if not np.isfinite(x):
        return "--"
    if x == 0:
        return "0"
    e = int(np.floor(np.log10(x)))
    return f"{x / 10 ** e:.1f}e{e}"


def fnum(x, d=1):
    """Ratio formatting: integers above 100, otherwise d decimals."""
    if not np.isfinite(x):
        return "??"
    return f"{x:.0f}" if abs(x) >= 100 else f"{x:.{d}f}"


gmean = lambda r: float(np.exp(np.log(r[np.isfinite(r) & (r > 0)]).mean()))


def put(tag, r, d=1):
    r = np.asarray(r, float)
    r = r[np.isfinite(r)]
    if len(r) == 0:
        NUM[tag + "Lo"] = NUM[tag + "Hi"] = NUM[tag + "Gm"] = "??"
        return
    NUM[tag + "Lo"], NUM[tag + "Hi"], NUM[tag + "Gm"] = fnum(r.min(), d), fnum(r.max(), d), fnum(gmean(r), d)


def tex_label(key):
    return STYLE[key][0].replace("(", "{\\footnotesize(").replace(")", ")}")


# ====================================================================== main accuracy matrix
acc = pd.read_csv("results/agg_acc.csv")
DATASETS = [d for d in DATASETS if d in set(acc.dataset)]
REAL = [d for d in DATASETS if d in ("taxi", "mawi", "dmv")]
# Primary metric: NAE on record-sampled queries with k >= 2 predicates (workload definition of the OmniSketch paper).
kmax = 6
acc["nae0m"] = acc[[f"nae0_k{k}" for k in range(2, kmax + 1)]].mean(axis=1)
dcols = [f"mre1_d{j}" for j in range(1, 7)]
kcols = [f"mre0_k{k}" for k in range(1, 7)]
nkcols = [f"nae0_k{k}" for k in range(1, 7)]
metrics = ["nae0", "nae0m", "nae1", "mre0", "mre1"] + dcols + kcols + nkcols
m = with_omnibest(mean_over_seeds(acc, metrics), metrics)
have = set(m.method.unique())
NUM["OmniFamilySize"] = str(len(OMNI_FAMILY))
NUM["RealStreams"] = {1: "one", 2: "two", 3: "three"}[len(REAL)]

# ---------------------------------------------------------------- Fig: NAE vs memory
order = ["omnijs", "omnibest", "res", RSC2, USSC2, CS2, CSH]
ncol_ = 4 if len(DATASETS) > 6 else 3
fig, axes = plt.subplots(2, ncol_, figsize=(7.0 if ncol_ == 3 else 8.6, 4.1), sharex=True)
for ax in axes.flat[len(DATASETS):]:
    ax.axis("off")
for ax, ds in zip(axes.flat, DATASETS):
    for key in order:
        s = m[(m.dataset == ds) & (m.method == key)].sort_values("mem")
        if len(s):
            line(ax, s.mem, s.nae0m, key)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_title(DS_LABEL[ds]); ax.grid(True, which="major")
    ax.set_xticks([0.5, 1, 2, 5, 10, 20, 50]); ax.set_xticklabels(["0.5", "1", "2", "5", "10", "20", "50"])
    ax.minorticks_off()
for ax in axes[-1]:
    ax.set_xlabel("memory (MB)")
for ax in axes[:, 0]:
    ax.set_ylabel("normalized abs. error")
h, l = axes[0, 0].get_legend_handles_labels()
fig.legend(h, l, loc="upper center", ncol=4, frameon=False, bbox_to_anchor=(0.5, 1.07))
fig.tight_layout()
fig.savefig("paper/figs/nae_vs_mem.pdf", bbox_inches="tight")
plt.close(fig)

# ---------------------------------------------------------------- Fig: relative error by anchor frequency
order2 = ["omnibest", "res", RSC2, USSC2, CS2, CSH]
panels = [("taxi", 1.0), ("taxi", 10.0), ("mawi", 10.0), ("zipf15_5M", 10.0)]
fig, axes = plt.subplots(1, 4, figsize=(7.0, 2.25), sharey=True)
xl = [r"$10^1$", r"$10^2$", r"$10^3$", r"$10^4$", r"$10^5$", r"$10^6$"]
for ax, (ds, mem) in zip(axes, panels):
    for key in order2:
        s = m[(m.dataset == ds) & (m.method == key) & (m.mem == mem)]
        if len(s):
            line(ax, np.arange(6), s[dcols].to_numpy(float)[0], key)
    ax.set_yscale("log"); ax.grid(True, which="major")
    ax.set_xticks(np.arange(6)); ax.set_xticklabels(xl)
    ax.set_title(f"{DS_LABEL[ds]}, {mem:g} MB"); ax.minorticks_off()
fig.supxlabel("frequency of the query's rarest predicate (decade)", fontsize=8, y=0.04)
axes[0].set_ylabel("mean relative error")
h, l = axes[0].get_legend_handles_labels()
fig.legend(h, l, loc="upper center", ncol=6, frameon=False, bbox_to_anchor=(0.5, 1.1))
fig.tight_layout()
fig.savefig("paper/figs/anchor.pdf", bbox_inches="tight")
plt.close(fig)

# ---------------------------------------------------------------- Table: main accuracy table at 1 and 10 MB
tab_methods = ["omnijs", "omnip", "omniv", "omnibest", "ash", "res", "uss", RSC1, RSC2, RSC3, USSC2, CS2, CSH]
with open("paper/tables/main_acc.tex", "w") as f:
    f.write("\\begin{tabular}{l" + "r" * len(DATASETS) + "}\n\\toprule\n")
    f.write("Method & " + " & ".join(DS_LABEL[d] for d in DATASETS) + " \\\\\n\\midrule\n")
    for mem in (1.0, 10.0):
        f.write(f"\\multicolumn{{{len(DATASETS) + 1}}}{{l}}{{\\emph{{{mem:g} MB}}}} \\\\\n")
        for key in tab_methods:
            cells = []
            for ds in DATASETS:
                s = m[(m.dataset == ds) & (m.mem == mem)]
                v = s[s.method == key].nae0m.to_numpy()
                best = s[s.method.isin(tab_methods)].nae0m.min()
                txt = fmt_sci(v[0]) if len(v) else "--"
                cells.append(f"\\textbf{{{txt}}}" if len(v) and v[0] <= best * 1.0001 else txt)
            f.write(tex_label(key) + " & " + " & ".join(cells) + " \\\\\n")
            if key in ("omnibest", "ash", "uss"):
                f.write("\\addlinespace[2pt]\n")
        if mem == 1.0:
            f.write("\\midrule\n")
    f.write("\\bottomrule\n\\end{tabular}\n")

# ---------------------------------------------------------------- Table: NAE by number of predicates (10 MB)
byk_methods = ["omnijs", "omnibest", "res", "uss", RSC1, RSC2, RSC3, USSC2, CS2, CSH]
with open("paper/tables/by_k.tex", "w") as f:
    f.write("\\begin{tabular}{l" + "r" * 12 + "}\n\\toprule\n")
    f.write(" & \\multicolumn{6}{c}{Taxi, number of predicates $k$} & \\multicolumn{6}{c}{MAWI, number of predicates $k$} \\\\\n")
    f.write("Method & " + " & ".join(str(k) for k in range(1, 7)) + " & " + " & ".join(str(k) for k in range(1, 7)) + " \\\\\n\\midrule\n")
    for key in byk_methods:
        cells = []
        for ds in ("taxi", "mawi"):
            s = m[(m.dataset == ds) & (m.mem == 10.0)]
            row = s[s.method == key]
            if not len(row):
                cells += ["--"] * 6
                continue
            v = row[nkcols].to_numpy(float)[0]
            best = s[s.method.isin(byk_methods)][nkcols].min().to_numpy(float)
            cells += [f"\\textbf{{{fmt_sci(x)}}}" if x <= b * 1.0001 else fmt_sci(x) for x, b in zip(v, best)]
        lab = STYLE[key][0].replace("OmniSketch", "Omni.").replace("(authors' code)", "(authors)")
        f.write(lab + " & " + " & ".join(cells) + " \\\\\n")
    f.write("\\bottomrule\n\\end{tabular}\n")

# ---------------------------------------------------------------- numbers quoted in the text
PA = m.pivot_table(index=["dataset", "mem"], columns="method", values="nae0")    # all k
P = m.pivot_table(index=["dataset", "mem"], columns="method", values="nae0m")    # k >= 2 (primary)
K = {k: m.pivot_table(index=["dataset", "mem"], columns="method", values=f"nae0_k{k}") for k in range(1, 7)}
nonunif = [i for i in P.index if i[0] != "unif_5M"]
col = lambda T, c: T[c] if c in T.columns else pd.Series(np.nan, index=T.index)

# OmniSketch variants versus the reservoir
put("JavaVsRes", col(P, "omnijs") / P["res"]);                      put("JavaVsResNonUnif", (col(P, "omnijs") / P["res"]).loc[nonunif])
put("OmniPaperVsRes", col(P, "omnip") / P["res"]);                  put("OmniJournalVsRes", (col(P, "omniv") / P["res"]).loc[nonunif])
put("OmniAnyVsRes", (P["omnibest"] / P["res"]).loc[nonunif]);       put("OmniAnyVsResUnif", (P["omnibest"] / P["res"]).loc["unif_5M"], 2)
put("JavaVsResAll", col(PA, "omnijs") / PA["res"], 0)
ten = [i for i in P.index if i[1] >= 10]            # the budgets of the original study
put("JavaVsResTen", (col(P, "omnijs") / P["res"]).loc[ten])
put("JavaVsResTenNonUnif", (col(P, "omnijs") / P["res"]).loc[[i for i in ten if i[0] != "unif_5M"]])
for ds, tag in [("mawi", "Mawi"), ("taxi", "Taxi"), ("dmv", "Dmv")]:
    if ds in DATASETS:
        put("JavaVsResTen" + tag, (col(P, "omnijs") / P["res"]).loc[[i for i in ten if i[0] == ds]])
put("OmniAnyVsCsh", (P["omnibest"] / col(P, CSH)).loc[nonunif])
put("OmniAnyVsCs", (P["omnibest"] / col(P, CS2)).loc[nonunif])
put("OmniAnyVsCalTwo", (P["omnibest"] / col(P, RSC2)).loc[nonunif])
put("JavaVsCalTwo", (col(P, "omnijs") / col(P, RSC2)).loc[nonunif])
for ds, tag in [("mawi", "Mawi"), ("taxi", "Taxi"), ("unif_5M", "Unif")]:
    put("JavaVsRes" + tag, (col(P, "omnijs") / P["res"]).loc[ds], 1)
    put("OmniAnyVsRes" + tag, (P["omnibest"] / P["res"]).loc[ds], 2 if ds == "unif_5M" else 1)
# how flat is the authors' configuration in memory?
for ds, tag in [("unif_5M", "Unif"), ("taxi", "Taxi"), ("zipf10_5M", "ZipfTen")]:
    if "omnijs" in P.columns:
        s = P["omnijs"].loc[ds]
        NUM[f"JavaFlat{tag}"] = fnum(s.loc[10.0] / s.loc[50.0], 2) if 10.0 in s.index and 50.0 in s.index else "??"
        r = P["res"].loc[ds]
        NUM[f"ResFlat{tag}"] = fnum(r.loc[10.0] / r.loc[50.0], 2) if 10.0 in r.index and 50.0 in r.index else "??"
# journal-version baselines on the main streams
put("AshVsRes", col(P, "ash") / P["res"], 2);       put("AshLeanVsRes", col(P, "ashl") / P["res"], 2)
put("ResjVsRes", col(P, "resj") / P["res"], 2)
put("OmniJournalVsAsh", (col(P, "omniv") / col(P, "ash")).loc[nonunif], 1)
put("OmniAnyVsAsh", (P["omnibest"] / col(P, "ash")).loc[nonunif], 1)
# calibration
put("CalOneVsRes", P["res"] / col(P, RSC1), 2);     put("CalOneVsResAll", PA["res"] / col(PA, RSC1), 2)
put("CalTwoVsRes", (P["res"] / col(P, RSC2)).loc[nonunif], 1);  put("CalTwoVsResUnif", (P["res"] / col(P, RSC2)).loc["unif_5M"], 2)
put("CalThreeVsRes", (P["res"] / col(P, RSC3)).loc[nonunif], 1)
put("CalThreeVsCalTwo", (col(P, RSC2) / col(P, RSC3)).loc[nonunif], 2)
put("CalOneKone", K[1]["res"] / col(K[1], RSC1), 0)
put("CalTwoKtwo", (K[2]["res"] / col(K[2], RSC2)).loc[nonunif], 0)
put("CalTwoKthree", (K[3]["res"] / col(K[3], RSC2)).loc[nonunif], 1)
put("CalThreeKthree", (K[3]["res"] / col(K[3], RSC3)).loc[nonunif], 1)
put("CalTwoKfour", (K[4]["res"] / col(K[4], RSC2)).loc[nonunif], 2)
put("CalThreeKfour", (K[4]["res"] / col(K[4], RSC3)).loc[nonunif], 2)
put("CalOneKtwo", (K[2]["res"] / col(K[2], RSC1)).loc[nonunif], 2)
put("UssCalTwoVsUss", (P["uss"] / col(P, USSC2)).loc[nonunif], 1)
for ds, tag in [("mawi", "Mawi"), ("taxi", "Taxi"), ("zipf10_5M", "ZipfTen"), ("zipf13_5M", "ZipfThirteen"), ("zipf15_5M", "ZipfFifteen")]:
    put("CalTwoVsRes" + tag, (P["res"] / col(P, RSC2)).loc[ds], 1)
    put("CalOneVsRes" + tag, (P["res"] / col(P, RSC1)).loc[ds], 2)
    put("UssVsRes" + tag, (P["res"] / P["uss"]).loc[ds], 1)
    put("UssCalTwoVsRes" + tag, (P["res"] / col(P, USSC2)).loc[ds], 1)
    put("CshVsRes" + tag, (P["res"] / col(P, CSH)).loc[ds], 1)
    put("CshVsUssCalTwo" + tag, (col(P, CSH) / col(P, USSC2)).loc[ds], 2)
    put("CsVsCalTwo" + tag, (col(P, CS2) / col(P, RSC2)).loc[ds], 2)
    put("CshVsCs" + tag, (col(P, CS2) / col(P, CSH)).loc[ds], 2)
put("CsVsCalTwo", col(P, CS2) / col(P, RSC2), 2)
put("CsVsRes", (P["res"] / col(P, CS2)).loc[nonunif], 1)
put("CshVsCalTwo", col(P, CSH) / col(P, RSC2), 2)
# best row-sample variant versus the reservoir
rowfam = [c for c in (RSC2, RSC3, USSC2, CS2, CSH) if c in P.columns]
if rowfam:
    put("BestRowVsRes", (P["res"] / P[rowfam].min(axis=1)).loc[nonunif], 1)
    put("OmniAnyVsBestRow", (P["omnibest"] / P[rowfam].min(axis=1)).loc[nonunif], 0)
    put("JavaVsBestRow", (col(P, "omnijs") / P[rowfam].min(axis=1)).loc[nonunif], 0)
    if CSH in P.columns:
        NUM["CshWorstFactor"] = fnum((P[CSH] / P[rowfam].min(axis=1)).loc[nonunif].max(), 1)
# counters only (the Count-Min bound of the anchor cell, no sample) versus the calibrated sample, k >= 3
cmc = [c for c in [f"cmnae0_k{k}" for k in range(3, 7)] if c in acc.columns]
if cmc:
    a2 = acc[acc.method == RSC2].copy()
    a2["cm3"] = a2[cmc].mean(axis=1)
    a2["nae3p"] = a2[[c.replace("cmnae0", "nae0") for c in cmc]].mean(axis=1)
    g2 = a2.groupby(["dataset", "mem"])[["cm3", "nae3p"]].mean()
    put("CmVsCalTwo", (g2.cm3 / g2.nae3p).loc[[i for i in g2.index if i[0] != "unif_5M"]], 1)

# rare anchors on Taxi
for mem, tag in [(5.0, "Five"), (10.0, "Ten"), (20.0, "Twenty"), (50.0, "Fifty")]:
    t = m[(m.dataset == "taxi") & (m.mem == mem)].set_index("method")
    if CS2 not in t.index:
        continue
    others = t.loc[[x for x in ["omnibest", "res", RSC1, RSC2, RSC3, "uss", USSC1, USSC2] if x in t.index]]
    for d, dn in [(1, "One"), (2, "Two")]:
        c = t.loc[CS2, f"mre1_d{d}"]
        o = others[f"mre1_d{d}"].min()
        NUM[f"TaxiD{dn}Cs{tag}"] = f"{c:.2g}"
        NUM[f"TaxiD{dn}Csh{tag}"] = f"{t.loc[CSH, f'mre1_d{d}']:.2g}" if CSH in t.index else "??"
        NUM[f"TaxiD{dn}Other{tag}"] = f"{o:.2g}"
        NUM[f"TaxiD{dn}Gain{tag}"] = fnum(o / c, 0) if o / c >= 10 else fnum(o / c, 1)
        NUM[f"TaxiD{dn}Res{tag}"] = f"{t.loc['res', f'mre1_d{d}']:.2g}"
tt = m[m.dataset == "taxi"].pivot_table(index="mem", columns="method", values=["mre1_d4", "mre1_d5"])
if CS2 in tt["mre1_d4"].columns:
    rr = pd.concat([tt["mre1_d4"][CS2] / tt["mre1_d4"][RSC2], tt["mre1_d5"][CS2] / tt["mre1_d5"][RSC2]])
    NUM["TaxiFreqPenaltyLo"], NUM["TaxiFreqPenaltyHi"] = fnum(rr.min(), 1), fnum(rr.max(), 1)
# anchors on MAWI / Zipf: CoStrata relative to RS+C2 over the two rarest decades
for ds, tag in [("mawi", "Mawi"), ("zipf15_5M", "ZipfFifteen")]:
    t = m[m.dataset == ds].pivot_table(index="mem", columns="method", values=["mre1_d1", "mre1_d2"])
    if CS2 in t["mre1_d1"].columns and RSC2 in t["mre1_d1"].columns:
        r = pd.concat([t["mre1_d1"][RSC2] / t["mre1_d1"][CS2], t["mre1_d2"][RSC2] / t["mre1_d2"][CS2]])
        NUM[f"AnchorNeutral{tag}Lo"], NUM[f"AnchorNeutral{tag}Hi"] = fnum(r.min(), 2), fnum(r.max(), 2)
# OmniSketch on rare anchors (best variant, rarest decade, 10 MB)
for ds, tag in [("taxi", "Taxi"), ("mawi", "Mawi")]:
    t = m[(m.dataset == ds) & (m.mem == 10.0)].set_index("method")
    NUM[f"OmniAnchorDOne{tag}"] = f"{t.loc['omnibest', 'mre1_d1']:.2g}" if "omnibest" in t.index else "??"

cst = acc[acc.method == CS2]
if len(cst):
    NUM["HtBiasMax"] = f"{100 * cst.groupby(['dataset', 'mem']).bias_ht.mean().abs().max():.1f}"

with open("paper/tables/ratios.txt", "w") as f:
    f.write("NAE (k>=2), mean over seeds\n" + P.to_string(float_format=lambda v: "%.3g" % v) + "\n")
    g = mean_over_seeds(acc[acc.method.isin(OMNI_FAMILY)], ["nae0m"])
    idx = g.groupby(["dataset", "mem"]).nae0m.idxmin()
    f.write("\nbest OmniSketch variant per (dataset, mem), NAE k>=2:\n" + g.loc[idx].to_string() + "\n")
    f.write("\nanchored queries, mean relative error by decade (10 MB):\n")
    f.write(m[m.mem == 10.0].set_index(["dataset", "method"])[dcols].to_string(float_format=lambda v: "%.3g" % v) + "\n")
    f.write("\nNAE by k (10 MB):\n")
    f.write(m[m.mem == 10.0].set_index(["dataset", "method"])[nkcols].to_string(float_format=lambda v: "%.3g" % v) + "\n")
    if len(cst):
        f.write("\nCoStrata 95% CI coverage (all queries / truth>=100), HT bias:\n")
        f.write(cst.groupby(["dataset", "mem"])[["cov95", "cov95_big", "bias_ht"]].mean().to_string() + "\n")
best_var = g.loc[idx]
NUM["BestVariantJavaCount"] = str(int((best_var.method == "omnijs").sum()))
NUM["BestVariantCells"] = str(len(best_var))

# ---------------------------------------------------------------- Table: CI coverage
if len(cst):
    cv = cst.groupby(["dataset", "mem"], as_index=False)[["cov95", "cov95_big"]].mean()
    with open("paper/tables/coverage.tex", "w") as f:
        mems = [1.0, 10.0, 50.0]
        f.write("\\begin{tabular}{l" + "rr" * len(mems) + "}\n\\toprule\n")
        f.write(" & " + " & ".join(f"\\multicolumn{{2}}{{c}}{{{x:g} MB}}" for x in mems) + " \\\\\n")
        f.write("Stream & " + " & ".join("all & $f{\\ge}100$" for _ in mems) + " \\\\\n\\midrule\n")
        for ds in DATASETS:
            cells = []
            for x in mems:
                s = cv[(cv.dataset == ds) & (cv.mem == x)]
                cells += [f"{100 * s.cov95.iloc[0]:.1f}", f"{100 * s.cov95_big.iloc[0]:.1f}"] if len(s) else ["--", "--"]
            f.write(DS_LABEL[ds] + " & " + " & ".join(cells) + " \\\\\n")
        f.write("\\bottomrule\n\\end{tabular}\n")

# ====================================================================== journal-version setting
if os.path.exists("results/agg_journal.csv"):
    jr = pd.read_csv("results/agg_journal.csv")
    jr["nae0m"] = jr[[c for c in [f"nae0_k{k}" for k in range(2, 7)] if c in jr.columns]].mean(axis=1)
    jm = jr.groupby(["dataset", "method", "mem"], as_index=False)[[c for c in ["nae3", "mre3", "nae0m", "nae0", "nae1", "mre1"] if c in jr.columns]].mean()
    jmeth = ["omniv", "omnijs", "omni:256:2", "omni:1024:2", "ash", "ashl", "resj", "res", RSC2]
    jlab = {"omniv": "OmniSketch (journal cfg.)", "omnijs": "OmniSketch (PVLDB code)", "omni:256:2": "OmniSketch ($256{\\times}2$, oracle)",
            "omni:1024:2": "OmniSketch ($1024{\\times}2$, oracle)", "ash": "aSH (journal acct.)", "ashl": "aSH (lean acct.)",
            "resj": "RS (journal acct.)", "res": "RS", RSC2: "RS+C2"}
    mems = [10.0, 25.0, 50.0, 100.0]
    with open("paper/tables/journal.tex", "w") as f:
        f.write("\\begin{tabular}{l" + "r" * 8 + "}\n\\toprule\n")
        f.write(" & \\multicolumn{4}{c}{journal-style queries ($k=1..9$)} & \\multicolumn{4}{c}{record-sampled queries ($k\\ge2$)} \\\\\n")
        f.write("Method & " + " & ".join(f"{x:g}" for x in mems) + " & " + " & ".join(f"{x:g}" for x in mems) + " MB \\\\\n\\midrule\n")
        for key in jmeth:
            cells = []
            for wl_, colname in (("jz13p9-vldbj", "nae3"), ("jz13p9_5M", "nae0m")):
                for x in mems:
                    s = jm[(jm.dataset == wl_) & (jm.method == key) & (jm.mem == x)]
                    cells.append(fmt_sci(s[colname].iloc[0]) if len(s) and colname in s else "--")
            f.write(jlab[key] + " & " + " & ".join(cells) + " \\\\\n")
            if key in ("omni:1024:2", "resj"):
                f.write("\\addlinespace[2pt]\n")
        f.write("\\bottomrule\n\\end{tabular}\n")
    R3 = jm[jm.dataset == "jz13p9-vldbj"].pivot_table(index="mem", columns="method", values="mre3")
    for key, tag in [("omniv", "Omni"), ("ash", "Ash"), ("res", "Res"), ("resj", "Resj"), (RSC2, "CalTwo")]:
        NUM[f"JSMre{tag}"] = f"{R3.loc[10.0, key]:.2f}" if key in R3.columns else "??"
    J3 = jm[jm.dataset == "jz13p9-vldbj"].pivot_table(index="mem", columns="method", values="nae3")
    J0 = jm[jm.dataset == "jz13p9_5M"].pivot_table(index="mem", columns="method", values="nae0m")
    note("journal setting, journal-style workload (NAE, all k)", J3)
    note("journal setting, record-sampled workload (NAE, k>=2)", J0)
    for T, tag in ((J3, "JS"), (J0, "JR")):
        c = lambda k: T[k] if k in T.columns else pd.Series(np.nan, index=T.index)
        put(f"{tag}OmniVsAsh", c("omniv") / c("ash"), 1)        # > 1: aSH (journal accounting) more accurate
        put(f"{tag}OmniVsRes", c("omniv") / c("res"), 1)
        put(f"{tag}AshVsRes", c("ash") / c("res"), 2)
        put(f"{tag}AshLeanVsRes", c("ashl") / c("res"), 2)
        put(f"{tag}ResjVsRes", c("resj") / c("res"), 2)
        put(f"{tag}OmniVsResj", c("omniv") / c("resj"), 1)
        put(f"{tag}PvldbVsRes", c("omnijs") / c("res"), 1)
        put(f"{tag}CalTwoVsRes", c("res") / c(RSC2), 1)

# ====================================================================== timing
if os.path.exists("results/time/time.log"):
    rows = []
    for ln in open("results/time/time.log"):
        if ln.startswith("TIME,"):
            f_ = ln.strip().split(",")
            rows.append(dict(dataset=f_[1], method=f_[2], mem=float(f_[3]), p=int(f_[4]), ns=float(f_[6]), mrps=float(f_[7])))
    t = pd.DataFrame(rows)
    qrows = []
    if os.path.exists("results/time/query.log"):
        for ln in open("results/time/query.log"):
            if ln.startswith("RUN,"):
                f_ = ln.strip().split(",")
                qrows.append(dict(dataset=f_[1], method=f_[2], mem=float(f_[3]), q_us=float(f_[9])))
    q = pd.DataFrame(qrows)
    tm = ["res", "uss", RSC1, RSC2, RSC3, CS1, CS2, CSH, USSC2, "ash", "omnij", "omnip", "omniv", "omnil", "omni:256:2", "omni:1024:2"]
    names = {"res": "RS", "uss": "USS", RSC1: "RS+C1", RSC2: "RS+C2", RSC3: "RS+C3", CS1: "CoStrata (C1)", CS2: "CoStrata",
             CSH: "CoStrata-H", USSC2: "USS+C2", "ash": "aSH", "omnij": "OmniSketch (port of authors' code)", "omnijs": "OmniSketch (port of authors' code)",
             "omnip": "OmniSketch (ours, same $B$)", "omniv": "OmniSketch (journal cfg.)",
             "omnil": "OmniSketch ($8{\\times}3$, lean)", "omni:256:2": "OmniSketch ($256{\\times}2$)",
             "omni:1024:2": "OmniSketch ($1024{\\times}2$)"}
    with open("paper/tables/throughput.tex", "w") as f:
        cols = [(ds, mem) for ds in ("taxi", "mawi") for mem in (1.0, 10.0, 50.0)]
        f.write("\\begin{tabular}{l" + "r" * len(cols) + "}\n\\toprule\n")
        f.write(" & \\multicolumn{3}{c}{Taxi ($p{=}12$)} & \\multicolumn{3}{c}{MAWI ($p{=}10$)} \\\\\n")
        f.write("Method & " + " & ".join(f"{mem:g} MB" for _, mem in cols) + " \\\\\n\\midrule\n")
        f.write(f"\\multicolumn{{{len(cols) + 1}}}{{l}}{{\\emph{{Update throughput (million records/s)}}}} \\\\\n")
        for key in tm:
            cells = []
            for ds, mem in cols:
                s = t[(t.dataset == ds) & (t.mem == mem) & (t.method == key)]
                cells.append((f"{s.mrps.iloc[0]:.1f}" if s.mrps.iloc[0] >= 1 else f"{s.mrps.iloc[0]:.2f}") if len(s) else "--")
            if any(c != "--" for c in cells):
                f.write(names[key] + " & " + " & ".join(cells) + " \\\\\n")
        if len(q):
            f.write("\\midrule\n")
            f.write(f"\\multicolumn{{{len(cols) + 1}}}{{l}}{{\\emph{{Query latency ($\\mu$s per query)}}}} \\\\\n")
            names[RSC2 + ":0:0:1"] = "RS+C2 with bucket index"
            for key in ["res", "uss", RSC2, RSC2 + ":0:0:1", CS2, CSH, "omnijs", "omnip", "omniv", "omni:256:2", "omni:1024:2"]:
                cells = []
                for ds, mem in cols:
                    s = q[(q.dataset == ds) & (q.mem == mem) & (q.method == key)]
                    cells.append((f"{s.q_us.iloc[0]:.1f}" if s.q_us.iloc[0] < 10 else f"{s.q_us.iloc[0]:.0f}") if len(s) else "--")
                if any(c != "--" for c in cells):
                    f.write(names[key] + " & " + " & ".join(cells) + " \\\\\n")
        f.write("\\bottomrule\n\\end{tabular}\n")
    tp = t[t.dataset.isin(["taxi", "mawi"])].pivot_table(index=["dataset", "mem"], columns="method", values="mrps")
    note("throughput (Mrec/s)", t[~t.dataset.str.startswith("census")].pivot_table(index=["dataset", "mem"], columns="method", values="mrps"))
    for meth, tag in [("res", "Res"), ("uss", "Uss"), (RSC1, "CalOne"), (RSC2, "CalTwo"), (RSC3, "CalThree"), (CS1, "CsOne"), (CS2, "Cs"),
                      (CSH, "Csh"), ("omnij", "Java"), ("omnip", "OmniPaper"), ("omnil", "OmniLean"), ("omniv", "OmniJournal"),
                      ("omni:256:2", "OmniMid"), ("ash", "Ash")]:
        if meth in tp.columns:
            lo, hi = tp[meth].min(), tp[meth].max()
            NUM[f"Tp{tag}Lo"], NUM[f"Tp{tag}Hi"] = (f"{lo:.1f}" if lo >= 1 else f"{lo:.2f}"), (f"{hi:.1f}" if hi >= 1 else f"{hi:.2f}")
    omni_cols = [c for c in ["omnij", "omnip", "omnil", "omniv", "omni:256:2", "omni:1024:2"] if c in tp.columns]
    if omni_cols:
        for meth, tag in [(RSC2, "CalTwo"), (CS2, "Cs"), (CS1, "CsOne"), (RSC3, "CalThree")]:
            if meth in tp.columns:
                r = tp[meth] / tp[omni_cols].max(axis=1)
                NUM[f"Tp{tag}VsOmniLo"], NUM[f"Tp{tag}VsOmniHi"] = fnum(r.min(), 1), fnum(r.max(), 1)
        if "omnij" in tp.columns and "omnip" in tp.columns:
            r = tp["omnip"] / tp["omnij"]
            NUM["TpOursVsJavaLo"], NUM["TpOursVsJavaHi"] = fnum(r.min(), 1), fnum(r.max(), 1)
    if RSC1 in tp.columns and RSC2 in tp.columns:
        r = tp[RSC1] / tp[RSC2]
        NUM["TpPairCostLo"], NUM["TpPairCostHi"] = fnum(r.min(), 1), fnum(r.max(), 1)
    if RSC2 in tp.columns and RSC3 in tp.columns:
        r = tp[RSC2] / tp[RSC3]
        NUM["TpTripleCostLo"], NUM["TpTripleCostHi"] = fnum(r.min(), 1), fnum(r.max(), 1)
    if len(q):
        ql = q.pivot_table(index=["dataset", "mem"], columns="method", values="q_us")
        note("query latency (us)", ql)
        oc = [c for c in ["omnijs", "omnip", "omniv", "omni:256:2", "omni:1024:2"] if c in ql.columns]
        for meth, tag in [(CS2, "Cs"), (RSC2, "CalTwo"), ("res", "Res")]:
            if meth in ql.columns and oc:
                r = ql[meth] / ql[oc].min(axis=1)
                NUM[f"Qlat{tag}VsOmniLo"], NUM[f"Qlat{tag}VsOmniHi"] = fnum(r.min(), 0), fnum(r.max(), 0)
                NUM[f"Qlat{tag}Lo"], NUM[f"Qlat{tag}Hi"] = fnum(ql[meth].min(), 0), fnum(ql[meth].max(), 0)
        if oc:
            NUM["QlatOmniLo"], NUM["QlatOmniHi"] = fnum(ql[oc].min().min(), 1), fnum(ql[oc].max().max(), 0)
        IDX = RSC2 + ":0:0:1"
        if IDX in ql.columns and oc:
            NUM["QlatIndexLo"], NUM["QlatIndexHi"] = fnum(ql[IDX].min(), 1), fnum(ql[IDX].max(), 0)
            r = ql[IDX] / ql[oc].min(axis=1)
            NUM["QlatIndexVsOmniLo"], NUM["QlatIndexVsOmniHi"] = fnum(r.min(), 1), fnum(r.max(), 1)
        if "omnijs" in ql.columns and CS2 in ql.columns:
            r = ql[CS2] / ql["omnijs"]
            NUM["QlatCsVsJavaLo"], NUM["QlatCsVsJavaHi"] = fnum(r.min(), 1), fnum(r.max(), 1)

    # attribute scaling on census
    if os.path.exists("results/agg_census.csv"):
        c = pd.read_csv("results/agg_census.csv")
        c["pp"] = c.dataset.str.extract(r"_p(\d+)").astype(int)
        c["nae0"] = c[[f"nae0_k{k}" for k in range(2, 5)]].mean(axis=1)   # k >= 2, as in the main results
        cm = c.groupby(["pp", "method"], as_index=False)[["nae0"]].mean()
        cg = cm[cm.method.isin(OMNI_FAMILY)].groupby("pp", as_index=False).nae0.min()
        cg["method"] = "omnibest"
        cm = pd.concat([cm, cg], ignore_index=True)
        tc = t[t.dataset.str.startswith("census_p")].copy()
        tc["pp"] = tc.dataset.str.extract(r"_p(\d+)").astype(int)
        fig, axes = plt.subplots(1, 2, figsize=(7.0, 2.3))
        for key in ["omnijs", "omnibest", "res", RSC2, USSC2, CS2]:
            s = cm[cm.method == key].sort_values("pp")
            if len(s):
                line(axes[0], s.pp, s.nae0, key)
        axes[0].set_xscale("log", base=2); axes[0].set_yscale("log")
        axes[0].set_ylabel("normalized abs. error"); axes[0].set_title("Accuracy at 10 MB")
        for key, tk in [("omnip", "omnip"), ("omnibest", "omni:256:2"), ("res", "res"), (RSC2, RSC2), (CS2, CS2), (RSC1, RSC1)]:
            s = tc[tc.method == tk].sort_values("pp")
            if len(s):
                line(axes[1], s.pp, s.mrps, key)
        if not len(tc):
            axes[1].plot([4, 68], [1, 1], alpha=0)
        axes[1].set_xscale("log", base=2); axes[1].set_yscale("log")
        axes[1].set_ylabel("million records / s"); axes[1].set_title("Update throughput at 10 MB")
        for ax in axes:
            ax.set_xticks([4, 8, 16, 32, 68]); ax.set_xticklabels(["4", "8", "16", "32", "68"])
            ax.set_xlabel("number of attributes $p$ (Census)"); ax.grid(True, which="major"); ax.minorticks_off()
        hs, ls_ = [], []
        for ax in axes:
            for h_, l_ in zip(*ax.get_legend_handles_labels()):
                if l_ not in ls_:
                    hs.append(h_); ls_.append(l_)
        fig.legend(hs, ls_, loc="upper center", ncol=4, frameon=False, bbox_to_anchor=(0.5, 1.2))
        fig.tight_layout()
        fig.savefig("paper/figs/census_scaling.pdf", bbox_inches="tight")
        plt.close(fig)
        import re as _re
        inf = c[(c.pp == 68) & (c.method == RSC2)]["info"]
        mw = _re.search(r"w2=(\d+)", inf.iloc[0]) if len(inf) else None
        NUM["CensusPairWidth"] = mw.group(1) if mw else "??"
        cp = cm.pivot(index="pp", columns="method", values="nae0")
        note("census scaling (NAE k>=2 by p)", cp)
        if len(tc):
            tpc = tc.pivot(index="pp", columns="method", values="mrps")
            note("census throughput (Mrec/s)", tpc)
            for meth, tag in [(CS2, "Cs"), (RSC2, "CalTwo"), (RSC1, "CalOne"), ("res", "Res"), ("omnip", "Omni"), ("uss", "Uss")]:
                if meth in tpc.columns:
                    for pp, pt in [(4, "Four"), (68, "SixtyEight")]:
                        v = tpc.loc[pp, meth]
                        NUM[f"CensusTp{tag}{pt}"] = f"{v:.1f}" if v >= 1 else f"{v:.2f}"
        for pp, tag in [(4, "Four"), (68, "SixtyEight")]:
            row = cp.loc[pp]
            g_ = lambda a, b: fnum(row[a] / row[b], 1) if a in row.index and b in row.index else "??"
            NUM[f"CensusOmniAnyVsRes{tag}"] = g_("omnibest", "res")
            NUM[f"CensusJavaVsRes{tag}"] = g_("omnijs", "res")
            NUM[f"CensusResVsCalTwo{tag}"] = g_("res", RSC2)
            NUM[f"CensusResVsCs{tag}"] = g_("res", CS2)
            NUM[f"CensusResVsCalOne{tag}"] = g_("res", RSC1)

# ====================================================================== ablations (10 MB)
if os.path.exists("results/agg_abl.csv"):
    a = pd.read_csv("results/agg_abl.csv")
    a["nae0m"] = a[[f"nae0_k{k}" for k in range(2, 7)]].mean(axis=1)
    am = a.groupby(["dataset", "method", "mem"], as_index=False)[["nae0", "nae0m", "nae0_ht", "nae0_k2", "nae0_k3", "mre1_d1", "mre1_d2"]].mean()
    rows = [("\\emph{pair-table share} $\\phi_2$ \\emph{(uniform sample)}", None),
            ("none (RS+C1)", "cs:0:1e9:1"), ("$\\phi_2=0.05$", "cs:0:1e9:1:0.05"), ("$\\phi_2=0.1$", "cs:0:1e9:1:0.1"),
            ("$\\phi_2=0.25$ (default)", "cs:0:1e9:1:0.25"), ("$\\phi_2=0.5$", "cs:0:1e9:1:0.5"),
            ("\\emph{pair and triple tables} $\\phi_2/\\phi_3$", None),
            ("$0.1/0.1$", "cs:0:1e9:1:0.1:0.1"), ("$0.15/0.15$ (RS+C3)", "cs:0:1e9:1:0.15:0.15"), ("$0.2/0.3$", "cs:0:1e9:1:0.2:0.3"),
            ("\\emph{all anchors (raking) instead of the lightest}", None),
            ("$\\phi_2=0.25$, all anchors", "cs:0:1e9:1:0.25:0:1"), ("$0.15/0.15$, all anchors", "cs:0:1e9:1:0.15:0.15:1"),
            ("\\emph{house weight} $\\lambda$ \\emph{(C1 only)}", None),
            ("$\\lambda=0$ (strata only)", "cs:0:0:1"), ("$\\lambda=1/4$", "cs:0:0.25:1"), ("$\\lambda=1$", "cs:0:1:1"),
            ("$\\lambda=4$", "cs:0:4:1"),
            ("\\emph{house weight with pair tables}", None),
            ("$\\lambda=1/4$, $\\phi_2=0.25$", "cs:0:0.25:1:0.25"), ("$\\lambda=1$, $\\phi_2=0.25$ (\\cs)", "cs:0:1:1:0.25"),
            ("\\emph{hold mode (}$\\phi_2=0.25$\\emph{)}", None),
            ("$\\lambda=1/4$", "csh:0:0.25:0.25"), ("$\\lambda=1$ (\\csh)", "csh:0:1:0.25"), ("$\\lambda=4$", "csh:0:4:0.25"),
            ("USS+C2 (no strata)", "ussc:0:0.25"),
            ("\\emph{allocation exponent} $\\gamma$ \\emph{(}$\\lambda=0$\\emph{, C1 only)}", None),
            ("$\\gamma=0.75$", "cs:0:0:0.75"), ("$\\gamma=0.5$", "cs:0:0:0.5"),
            ("\\emph{buckets per attribute} $w$ \\emph{(}$\\lambda=1$\\emph{, C1 only)}", None),
            ("$w=64$", "cs:64:1:1"), ("$w=256$", "cs:256:1:1"), ("$w=1024$", "cs:1024:1:1"), ("$w=4096$", "cs:4096:1:1"),
            ("$w=16384$", "cs:16384:1:1"),
            ("\\emph{reference}", None), ("RS", "res")]
    with open("paper/tables/ablation.tex", "w") as f:
        f.write("\\begin{tabular}{lrrrrrr}\n\\toprule\n")
        f.write(" & \\multicolumn{3}{c}{Taxi} & \\multicolumn{3}{c}{MAWI} \\\\\n")
        f.write("Variant & NAE & anchor $10^1$ & anchor $10^2$ & NAE & anchor $10^1$ & anchor $10^2$ \\\\\n\\midrule\n")
        for label, spec in rows:
            if spec is None:
                f.write(f"\\multicolumn{{7}}{{l}}{{{label}}} \\\\\n")
                continue
            cells = []
            for ds in ("taxi", "mawi"):
                r = am[(am.dataset == ds) & (am.mem == 10.0) & (am.method == spec)]
                cells += [fmt_sci(r.nae0m.iloc[0]), f"{r.mre1_d1.iloc[0]:.2g}", f"{r.mre1_d2.iloc[0]:.2g}"] if len(r) else ["--"] * 3
            f.write("\\quad " + label + " & " + " & ".join(cells) + " \\\\\n")
        f.write("\\bottomrule\n\\end{tabular}\n")
    for ds, tag in [("taxi", "Taxi"), ("mawi", "Mawi"), ("zipf15_5M", "Zipf")]:
        r = am[(am.dataset == ds) & (am.mem == 10.0) & (am.method == "cs:0:1:1:0.25")]
        if len(r):
            NUM["HtVsCal" + tag] = fnum(r.nae0_ht.iloc[0] / r.nae0.iloc[0], 1)
    note("ablation (mean over seeds)", am)

# ====================================================================== unbiasedness check
if os.path.exists("results/bias.json"):
    bj = json.load(open("results/bias.json"))
    for key, tag in [("taxi/cs", "Taxi"), ("zipf15_5M/cs", "Zipf"), ("mawi/cs", "Mawi"),
                     ("taxi/csh", "HoldTaxi"), ("zipf15_5M/csh", "HoldZipf"), ("mawi/csh", "HoldMawi")]:
        if key in bj:
            b = bj[key]
            NUM[f"BiasSeeds{tag}"] = str(b["seeds"])
            NUM[f"BiasQueries{tag}"] = str(b["queries"])
            NUM[f"BiasHtTotal{tag}"] = f"{100 * b['ht']['rel_bias_total']:+.2f}"
            NUM[f"BiasHtTotalSe{tag}"] = f"{100 * b['ht']['rel_bias_total_se']:.2f}"
            NUM[f"BiasHtZ{tag}"] = f"{100 * b['ht']['frac_absz_lt_1_96']:.0f}"
            NUM[f"BiasCalTotal{tag}"] = f"{100 * b['cal']['rel_bias_total']:+.2f}"
            NUM[f"BiasCalTotalSe{tag}"] = f"{100 * b['cal']['rel_bias_total_se']:.2f}"

if os.path.exists("results/bias_probe.txt"):
    import re as _re2
    blocks = open("results/bias_probe.txt").read().split("\n")
    cur, got = None, []
    for ln in blocks:
        mh = _re2.match(r"(\w+) N=(\d+) p=(\d+) M=(\d+) .*seeds=(\d+)", ln)
        if mh:
            cur = (mh.group(1), int(mh.group(3)), int(mh.group(5)))
        mt = _re2.match(r"\s+total N\s+mean rel. deviation ([+-][0-9.]+)%\s+\+- ([0-9.]+)%", ln)
        if mt and cur:
            got.append((cur, mt.group(1), mt.group(2)))
    for (variant, pp_, seeds_), mean_, se_ in got:
        tag = "Cs" if variant == "cs" else ("Csh" if pp_ == 3 else "CshRep")
        NUM[f"Probe{tag}"] = f"{float(mean_):+.3f}"
        NUM[f"Probe{tag}Se"] = f"{float(se_):.3f}"
        NUM["ProbeSeeds"] = f"{seeds_:,}".replace(",", "\\,")

# ====================================================================== validation against the authors' code
if os.path.exists("results/validation/numbers.json"):
    v = json.load(open("results/validation/numbers.json"))
    NUM["ValRuns"] = str(v["runs"])
    NUM["ValQueries"] = f"{v['queries']:,}".replace(",", "\\,")
    NUM["ValStreams"] = str(v["streams"])
    NUM["ValBoundFracLo"] = f"{100 * v['bound_frac_min']:.0f}"
    NUM["ValBoundFracHi"] = f"{100 * v['bound_frac_max']:.0f}"
    NUM["JavaThmNaeLo"] = f"{v['thm_nae_min']:.2f}"
    NUM["JavaThmNaeHi"] = f"{v['thm_nae_max']:.1f}"
    for ds, tag in [("jzipf10_5M", "ZipfTen"), ("jzipf13_5M", "ZipfThirteen"), ("jzipf15_5M", "ZipfFifteen"), ("taxi", "Taxi"),
                    ("mawi", "Mawi"), ("unif_5M", "Unif"), ("zipf10_5M", "OurZipfTen")]:
        NUM[f"RowOneShare{tag}"] = f"{100 * v['max_cell_share'][ds]:.0f}"
    NUM["JavaUsedTen"] = f"{v['java_used_mb_10']['jzipf10_5M']:.2f}"
    NUM["JavaNominalTen"] = f"{v['nominal_mb_10']['jzipf10_5M']:.2f}"

    def sci_tex(x):
        e = int(np.floor(np.log10(x)))
        return f"{x / 10 ** e:.1f}\\cdot10^{{{e}}}"
    for ds, tag in [("jzipf10_5M", "ZipfTen"), ("jzipf13_5M", "ZipfThirteen"), ("jzipf15_5M", "ZipfFifteen"), ("unif_5M", "Unif")]:
        for mem, mt in [(10, "Ten"), (50, "Fifty"), (100, "Hundred")]:
            NUM[f"FigSix{tag}{mt}"] = sci_tex(v["k3"][f"{ds}@{mem}"])
    jo = v["java_over_ours"]
    grp = lambda pre: np.array([x for k_, x in jo.items() if k_.split("@")[0] in pre and float(k_.split("@")[1]) >= 10])
    put("JavaOverOursZipfDag", grp(("jzipf10_5M", "jzipf13_5M", "jzipf15_5M")), 1)
    put("JavaOverOursZipf", grp(("zipf10_5M", "zipf13_5M", "zipf15_5M", "unif_5M")), 2)
    put("OursOverJavaTaxi", 1 / grp(("taxi",)), 1)
    put("JavaOverOursMawi", grp(("mawi",)), 2)

NUM.setdefault("CsThroughput", "??")
with open("paper/numbers.tex", "w") as f:
    f.write("% generated by scripts/plots.py -- do not edit\n")
    for k, v_ in sorted(NUM.items()):
        f.write(f"\\newcommand{{\\Num{k}}}{{{v_}}}\n")
with open("results/summary.txt", "w") as f:
    f.write("".join(NOTES))
print("figures and tables written;", len(NUM), "numbers")
