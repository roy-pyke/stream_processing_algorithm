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
USSC3, CSH3 = "ussc:0:0.15:0.15", "csh:0:1:0.15:0.15"
COCO, COCOC1, COCOC2 = "coco", "cococ", "cococ:0:0.25"
IX = "/ix=2"   # query index charged to the memory budget
# Fixed entity -> (label, colour, marker, linestyle). Colour follows the method in every figure.
STYLE = {
    "omnijs":   ("OmniSketch (authors' code)", "#e87ba4", "X", "--"),
    "omnip":    ("OmniSketch (ours, same $B$)", "#b05c7d", "x", "--"),
    "omniv":    ("OmniSketch (journal cfg.)",  "#9a6fb0", "+", "--"),
    "omnibest": ("OmniSketch (best)",          "#008300", "D", "--"),
    "hydra":    ("Hydra (authors' adaptation)", "#7a5c2e", "d", "-."),
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
    USSC3:      ("USS+C3",                     "#b87a00", "v", "--"),
    CSH3:       ("CoStrata-H (C3)",            "#0b6f7a", "h", "--"),
    COCO:       ("CocoSketch",                 "#d6457a", "p", "-"),
    COCOC2:     ("Coco+C2",                    "#8f1f4b", "p", "-"),
    CS2 + IX:   ("CoStrata, index charged",    "#2a78d6", "o", "--"),
    RSC2 + IX:  ("RS+C2, index charged",       "#6250d6", "P", "--"),
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
    if round(x / 10 ** e, 1) >= 10:
        e += 1
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
if os.path.exists("results/agg_acc3.csv"):   # triple-table variants of USS and CoStrata-H (run separately)
    acc = pd.concat([acc, pd.read_csv("results/agg_acc3.csv")], ignore_index=True)
if os.path.exists("results/agg_hydra.csv"):
    hy = pd.read_csv("results/agg_hydra.csv")
    acc = pd.concat([acc, hy[hy.dataset.isin(set(acc.dataset))]], ignore_index=True)
if os.path.exists("results/agg_acc2.csv"):   # CocoSketch and index-charged runs (scripts/run_presub.sh main)
    a2_ = pd.read_csv("results/agg_acc2.csv")
    acc = pd.concat([acc, a2_[a2_.method.isin([COCO, COCOC1, COCOC2]) | a2_.method.str.endswith(IX)]], ignore_index=True)
if os.path.exists("results/agg_dmv.csv"):    # full DMV matrix (scripts/run_presub.sh dmvfull) replaces the partial one
    acc = pd.concat([acc[acc.dataset != "dmv"], pd.read_csv("results/agg_dmv.csv")], ignore_index=True)
DATASETS =[d for d in DATASETS if d in set(acc.dataset)]
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
order = ["omnijs", "omnibest", "hydra", "res", COCO, RSC2, USSC3, CSH]
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
hs_, ls__ = [], []
for ax in axes.flat:                     # Hydra is missing in some panels: collect the legend from all
    for h_, l_ in zip(*ax.get_legend_handles_labels()):
        if l_ not in ls__:
            hs_.append(h_); ls__.append(l_)
fig.tight_layout()
fig.legend(hs_, ls__, loc="upper center", ncol=4, frameon=False, bbox_to_anchor=(0.5, 1.1))
fig.savefig("paper/figs/nae_vs_mem.pdf", bbox_inches="tight")
plt.close(fig)

# ---------------------------------------------------------------- Fig: relative error by anchor frequency
order2 = ["omnibest", "res", RSC2, USSC2, CS1, CS2, CSH]
panels = [("taxi", 10.0), ("taxi", 20.0), ("mawi", 10.0), ("zipf15_5M", 10.0)]
FLOOR = 2e-3   # smaller errors (essentially exact answers) are drawn at the floor
fig, axes = plt.subplots(1, 4, figsize=(7.0, 2.25), sharey=True)
xl = [r"$10^1$", r"$10^2$", r"$10^3$", r"$10^4$", r"$10^5$", r"$10^6$"]
for ax, (ds, mem) in zip(axes, panels):
    for key in order2:
        s = m[(m.dataset == ds) & (m.method == key) & (m.mem == mem)]
        if len(s):
            line(ax, np.arange(6), np.maximum(s[dcols].to_numpy(float)[0], FLOOR), key)
    ax.set_yscale("log"); ax.grid(True, which="major"); ax.set_ylim(FLOOR * 0.8, 4)
    ax.set_xticks(np.arange(6)); ax.set_xticklabels(xl)
    ax.set_title(f"{DS_LABEL[ds]}, {mem:g} MB"); ax.minorticks_off()
fig.supxlabel("frequency of the query's rarest predicate (decade)", fontsize=8, y=0.04)
axes[0].set_ylabel("mean relative error")
h, l = axes[0].get_legend_handles_labels()
fig.legend(h, l, loc="upper center", ncol=4, frameon=False, bbox_to_anchor=(0.5, 1.2))
fig.tight_layout()
fig.savefig("paper/figs/anchor.pdf", bbox_inches="tight")
plt.close(fig)

# ---------------------------------------------------------------- Table: main accuracy table at 1 and 10 MB
tab_methods = ["omnijs", "omnibest", "hydra", "ash", "res", "uss", COCO, RSC1, RSC2, RSC3, USSC2, USSC3, COCOC2, CS2, CSH, CSH3]
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
            if key in ("omnibest", "ash", COCO, RSC3, COCOC2):
                f.write("\\addlinespace[2pt]\n")
        if mem == 1.0:
            f.write("\\midrule\n")
    f.write("\\bottomrule\n\\end{tabular}\n")

# ---------------------------------------------------------------- Table: NAE by number of predicates (10 MB)
byk_methods = ["omnijs", "omnibest", "res", "uss", RSC1, RSC2, RSC3, USSC2, USSC3, CS2, CSH, CSH3]
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
put("JavaVsResTenSyn", (col(P, "omnijs") / P["res"]).loc[[i for i in ten if i[0].startswith(("zipf", "unif"))]], 0)
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
# Hydra: 5-attribute Zipf streams, and (scripts/run_hydra.sh with HYDRA_REAL=1) Taxi and MAWI
if "hydra" in P.columns:
    hyi = [i for i in P.index if np.isfinite(P.loc[i, "hydra"]) and i[0].startswith("zipf")]
    put("ResVsHydra", (P["res"] / P["hydra"]).loc[hyi], 1)            # > 1: Hydra more accurate than the reservoir
    put("OmniAnyVsHydra", (P["omnibest"] / P["hydra"]).loc[hyi], 1)
    put("HydraVsCalTwo", (P["hydra"] / col(P, RSC2)).loc[hyi], 2)     # > 1: RS+C2 more accurate than Hydra
    put("HydraVsCalThree", (P["hydra"] / col(P, RSC3)).loc[hyi], 2)
    hyr = [i for i in P.index if np.isfinite(P.loc[i, "hydra"]) and i[0] in ("taxi", "mawi")]
    if hyr:
        put("ResVsHydraReal", (P["res"] / P["hydra"]).loc[hyr], 1)
        put("OmniAnyVsHydraReal", (P["omnibest"] / P["hydra"]).loc[hyr], 1)
        put("HydraVsCalTwoReal", (P["hydra"] / col(P, RSC2)).loc[hyr], 1)
        put("HydraVsCsReal", (P["hydra"] / col(P, CS2)).loc[hyr], 1)
        if USSC3 in P.columns:
            put("HydraVsUssCalThreeReal", (P["hydra"] / P[USSC3]).loc[hyr], 1)
        hr_ = acc[(acc.method == "hydra") & acc.dataset.isin(["taxi", "mawi"])]
        for ds_, t_ in (("taxi", "Taxi"), ("mawi", "Mawi")):
            x_ = hr_[hr_.dataset == ds_]
            if len(x_):
                NUM[f"TpHydra{t_}"] = f"{1e3 / x_.ns_per_rec.median():.3f}"   # million records/s (accuracy runs, 4 in parallel)
    hrun = acc[acc.method == "hydra"]
    if len(hrun):
        NUM["TpHydraLo"], NUM["TpHydraHi"] = f"{1e3 / hrun.ns_per_rec.max():.2f}", f"{1e3 / hrun.ns_per_rec.min():.2f}"
if os.path.exists("results/time/hydra.log"):
    for ln in open("results/time/hydra.log"):
        if ln.startswith("HRUN,"):
            NUM["TpHydraIdle"] = f"{1e3 / float(ln.split(',')[7]):.2f}"
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
rowfam = [c for c in (RSC2, RSC3, USSC2, USSC3, CS2, CSH, CSH3) if c in P.columns]
for key_, tag_ in [(CSH3, "CshThree"), (USSC3, "UssCalThree")]:
    if key_ in P.columns:
        put("OmniAnyVs" + tag_, (P["omnibest"] / P[key_]).loc[nonunif])
        put("ResVs" + tag_, (P["res"] / P[key_]).loc[nonunif], 1)
        NUM[tag_ + "WorstFactor"] = fnum((P[key_] / P[rowfam].min(axis=1)).loc[nonunif].max(), 1)
        if "hydra" in P.columns:
            hyi_ = [i for i in P.index if np.isfinite(P.loc[i, "hydra"]) and i[0].startswith("zipf")]
            put("HydraVs" + tag_, (P["hydra"] / P[key_]).loc[hyi_], 2)
        for ds_, t2 in [("mawi", "Mawi"), ("taxi", "Taxi"), ("dmv", "Dmv"), ("zipf15_5M", "ZipfFifteen")]:
            if ds_ in DATASETS:
                put("ResVs" + tag_ + t2, (P["res"] / P[key_]).loc[ds_], 1)
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
BESTGAIN = []
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
        if CS1 in t.index:
            c1 = t.loc[CS1, f"mre1_d{d}"]
            NUM[f"TaxiD{dn}CsOne{tag}"] = f"{c1:.2g}"
            NUM[f"TaxiD{dn}GainOne{tag}"] = fnum(o / c1, 0) if o / c1 >= 10 else fnum(o / c1, 1)
            BESTGAIN.append(o / c1)
        BESTGAIN.append(o / c)
if BESTGAIN:
    NUM["TaxiAnchorBestGain"] = fnum(max(BESTGAIN), 0)
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
    if os.path.exists("results/journal_numbers.json"):
        jn = json.load(open("results/journal_numbers.json"))
        NUM["JKAshAbsLo"], NUM["JKAshAbsHi"] = fnum(jn["ash_abs_lo"], 1), fnum(jn["ash_abs_hi"], 1)
        NUM["JKResAbsLo"], NUM["JKResAbsHi"] = fnum(jn["res_abs_lo"], 1), fnum(jn["res_abs_hi"], 1)
        NUM["JKOmniRelHigh"], NUM["JKAshRelHigh"], NUM["JKResRelHigh"] = (f"{jn[k_]:.2f}" for k_ in ("omni_rel_high", "ash_rel_high", "res_rel_high"))
        NUM["JKOmniZeroHigh"] = f"{100 * jn['omni_zero_high']:.0f}"
        NUM["JKTruthHigh"] = f"{jn['truth_high']:.1f}"
        NUM["JKShareHigh"] = f"{100 * jn['share_high']:.0f}"
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
    t = pd.DataFrame(rows).drop_duplicates(["dataset", "method", "mem"])   # a few configurations are timed twice; keep the first
    qrows = []
    if os.path.exists("results/time/query.log"):
        for ln in open("results/time/query.log"):
            if ln.startswith("RUN,"):
                f_ = ln.strip().split(",")
                qrows.append(dict(dataset=f_[1], method=f_[2], mem=float(f_[3]), q_us=float(f_[9])))
    q = pd.DataFrame(qrows)
    tm = ["res", "uss", COCO, RSC1, RSC2, RSC3, CS1, CS2, CSH, CSH3, USSC2, USSC3, COCOC2, "ash", "omnij", "omnip", "omniv", "omnil", "omni:256:2", "omni:1024:2"]
    names = {"res": "RS", "uss": "USS", RSC1: "RS+C1", RSC2: "RS+C2", RSC3: "RS+C3", CS1: "CoStrata (C1)", CS2: "CoStrata",
             CSH: "CoStrata-H", CSH3: "CoStrata-H (C3)", USSC2: "USS+C2", USSC3: "USS+C3", "ash": "aSH",
             COCO: "CocoSketch", COCOC2: "Coco+C2", "omnij": "OmniSketch (port of authors' code)", "omnijs": "OmniSketch (port of authors' code)",
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
                      (CSH, "Csh"), (CSH3, "CshThree"), (USSC3, "UssCalThree"), (USSC2, "UssCalTwo"), ("omnij", "Java"), ("omnip", "OmniPaper"), ("omnil", "OmniLean"), ("omniv", "OmniJournal"),
                      ("omni:256:2", "OmniMid"), ("ash", "Ash")]:
        if meth in tp.columns:
            lo, hi = tp[meth].min(), tp[meth].max()
            NUM[f"Tp{tag}Lo"], NUM[f"Tp{tag}Hi"] = (f"{lo:.1f}" if lo >= 1 else f"{lo:.2f}"), (f"{hi:.1f}" if hi >= 1 else f"{hi:.2f}")
    tz = t[(t.dataset == "zipf13_5M")].pivot_table(index="mem", columns="method", values="mrps")
    for meth, tag in [(RSC2, "CalTwo"), (RSC3, "CalThree"), (RSC1, "CalOne"), (USSC3, "UssCalThree"), (CSH3, "CshThree")]:
        if meth in tz.columns and 10.0 in tz.index:
            NUM[f"TpZipf{tag}"] = f"{tz.loc[10.0, meth]:.1f}"
    omni_cols = [c for c in ["omnij", "omnip", "omnil", "omniv", "omni:256:2", "omni:1024:2"] if c in tp.columns]
    if omni_cols:
        for meth, tag in [(RSC2, "CalTwo"), (CS2, "Cs"), (CS1, "CsOne"), (RSC3, "CalThree")]:
            if meth in tp.columns:
                r = tp[meth] / tp[omni_cols].max(axis=1)
                NUM[f"Tp{tag}VsOmniLo"], NUM[f"Tp{tag}VsOmniHi"] = fnum(r.min(), 1), fnum(r.max(), 1)
        if RSC2 in tp.columns and CS2 in tp.columns:   # RS+C2 and CoStrata against the fastest OmniSketch, per budget
            r2 = pd.concat([tp[RSC2] / tp[omni_cols].max(axis=1), tp[CS2] / tp[omni_cols].max(axis=1)])
            for mem_, tag_ in ((1.0, "One"), (50.0, "Fifty")):
                rr_ = r2[r2.index.get_level_values("mem") == mem_]
                if len(rr_):
                    NUM[f"TpCalTwoVsOmni{tag_}Lo"], NUM[f"TpCalTwoVsOmni{tag_}Hi"] = fnum(rr_.min(), 1), fnum(rr_.max(), 1)
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
        # our implementation only: the port of the authors' tree-based code is listed separately
        oc = [c for c in ["omnip", "omniv", "omni:256:2", "omni:1024:2"] if c in ql.columns]
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
        fig, axes = plt.subplots(1, 2, figsize=(3.45, 2.0))   # one column of the two-column layout
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
            ax.set_xlabel("attributes $p$"); ax.grid(True, which="major"); ax.minorticks_off()
            ax.tick_params(labelsize=6.5); ax.title.set_fontsize(7.5); ax.yaxis.label.set_fontsize(7)
        hs, ls_ = [], []
        for ax in axes:
            for h_, l_ in zip(*ax.get_legend_handles_labels()):
                if l_ not in ls_:
                    hs.append(h_); ls_.append(l_)
        fig.tight_layout(w_pad=0.6)
        fig.legend(hs, ls_, loc="upper center", ncol=3, frameon=False, bbox_to_anchor=(0.5, 1.22), fontsize=6.3,
                   handlelength=1.8, columnspacing=0.8)
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

# sampled pair tables on the 68-attribute stream (negative result)
if os.path.exists("results/agg_census_cap.csv"):
    cc = pd.read_csv("results/agg_census_cap.csv").groupby("method").nae0_k2.mean()
    for key, tag in [(RSC1, "CalOne"), (RSC2, "Exact"), (RSC2 + ":0:0:0:256", "Sampled")]:
        if key in cc.index:
            e = int(np.floor(np.log10(cc[key])))
            NUM[f"CapKtwo{tag}"] = f"{cc[key] / 10 ** e:.1f}\\cdot10^{{{e}}}"

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

# ====================================================================== experiments added for the submission
# (scripts/run_presub.sh; aggregated into results/agg_{acc2,range,merge,pairs}.csv)
def geo(x):
    x = np.asarray(x, float)
    x = x[np.isfinite(x) & (x > 0)]
    return float(np.exp(np.log(x).mean())) if len(x) else np.nan


# ---------------------------------------------------------------- CocoSketch and the index charged to the budget
if COCO in P.columns:
    put("ResVsCoco", (P["res"] / P[COCO]).loc[nonunif], 2)              # > 1: CocoSketch more accurate than RS
    put("CocoVsUss", (P[COCO] / P["uss"]).loc[nonunif], 2)
    put("OmniAnyVsCoco", (P["omnibest"] / P[COCO]).loc[nonunif], 1)
    put("CocoVsBestRow", (P[COCO] / P[rowfam].min(axis=1)).loc[nonunif], 1)
    put("CocoVsCsh", (P[COCO] / col(P, CSH)).loc[nonunif], 1)
    if COCOC2 in P.columns:
        put("CocoCalTwoVsCoco", (P[COCO] / P[COCOC2]).loc[nonunif], 1)
        put("CocoCalTwoVsUssCalTwo", (P[COCOC2] / col(P, USSC2)).loc[nonunif], 2)
        put("CocoCalTwoVsBestRow", (P[COCOC2] / P[rowfam].min(axis=1)).loc[nonunif], 2)
        put("OmniAnyVsCocoCalTwo", (P["omnibest"] / P[COCOC2]).loc[nonunif], 0)
        for ds, tag in [("mawi", "Mawi"), ("taxi", "Taxi"), ("zipf15_5M", "ZipfFifteen"), ("zipf13_5M", "ZipfThirteen")]:
            put("CocoCalTwoVsCoco" + tag, (P[COCO] / P[COCOC2]).loc[ds], 1)
            put("CocoCalTwoVsUssCalTwo" + tag, (P[COCOC2] / col(P, USSC2)).loc[ds], 2)
    for ds, tag in [("mawi", "Mawi"), ("taxi", "Taxi")]:
        put("CocoVsUss" + tag, (P[COCO] / P["uss"]).loc[ds], 2)
for base, tag in [(CS2, "Cs"), (RSC2, "CalTwo"), (CS1, "CsOne")]:
    if base + IX in P.columns:
        put("IndexCost" + tag, (P[base + IX] / P[base]).loc[nonunif], 2)    # error with / without the index charged
        put("OmniAnyVsIndexed" + tag, (P["omnibest"] / P[base + IX]).loc[nonunif], 1)
        put("ResVsIndexed" + tag, (P["res"] / P[base + IX]).loc[nonunif], 1)
if CS2 + IX in P.columns:
    r = pd.concat([(P[CS2 + IX] / P[CS2]).loc[nonunif], (P[RSC2 + IX] / P[RSC2]).loc[nonunif]])
    put("IndexCost", r, 2)
    rr = pd.concat([(P["omnibest"] / P[CS2 + IX]).loc[nonunif], (P["omnibest"] / P[RSC2 + IX]).loc[nonunif]])
    put("OmniAnyVsIndexed", rr, 1)

# ---------------------------------------------------------------- small-count confidence intervals
if os.path.exists("results/agg_ci.csv"):   # scripts/run_presub.sh ci: runs with the anchor columns
    c2 = pd.read_csv("results/agg_ci.csv")
    c2 = c2[c2.method.isin([CS1, CS2, RSC1, RSC2, CSH]) & c2.ci95.notna()]
    cvm = c2.groupby(["method", "mem", "dataset"], as_index=False)[["cov95", "cov95_small", "cov95_big", "ci95", "ci95_small",
                                                                     "ci95_big", "ci95_zero", "ci95_f", "ci95_f_big", "frac_zero",
                                                                     "ci_relwidth", "ci_relwidth_f"]].mean()
    cva = cvm.groupby(["method", "mem"], as_index=False).mean(numeric_only=True)   # mean over the streams
    with open("paper/tables/coverage.tex", "w") as f:
        f.write("\\begin{tabular}{llrrrrrrr}\n\\toprule\n")
        f.write(" & & \\multicolumn{2}{c}{normal} & \\multicolumn{1}{c}{matches} & \\multicolumn{4}{c}{rarer side (Sec.~\\ref{sec:ci})} \\\\\n")
        f.write("\\cmidrule(lr){3-4}\\cmidrule(lr){5-5}\\cmidrule(lr){6-9}\n")
        f.write("Method & MB & $f{<}100$ & $f{\\ge}100$ & $f{\\ge}100$ & all & $f{<}100$ & $f{\\ge}100$ & $m{=}0$ \\\\\n\\midrule\n")
        for key, lab in [(CS2, "\\cs"), (RSC2, "RS+C2"), (CSH, "\\csh")]:
            for i, x in enumerate([1.0, 10.0, 50.0]):
                s = cva[(cva.method == key) & (cva.mem == x)]
                if not len(s):
                    continue
                s = s.iloc[0]
                pc = lambda v: f"{100 * v:.0f}" if np.isfinite(v) else "--"
                f.write((lab if i == 0 else "") + f" & {x:g} & " + " & ".join(pc(s[c]) for c in
                        ["cov95_small", "cov95_big", "ci95_f_big", "ci95", "ci95_small", "ci95_big", "ci95_zero"]) + " \\\\\n")
            if key != CSH:
                f.write("\\addlinespace[2pt]\n")
        f.write("\\bottomrule\n\\end{tabular}\n")
    note("interval coverage by stream", cvm)
    main3 = cvm[cvm.method.isin([CS2, RSC2, CSH])]
    NUM["CiNormSmallHi"] = f"{100 * main3.cov95_small.max():.0f}"
    NUM["CiNormSmallOneMb"] = f"{100 * cva[(cva.method == CS2) & (cva.mem == 1.0)].cov95_small.iloc[0]:.0f}"
    NUM["CiNormAllOneMb"] = f"{100 * cva[(cva.method == CS2) & (cva.mem == 1.0)].cov95.iloc[0]:.0f}"
    NUM["CiNewLo"], NUM["CiNewHi"] = f"{100 * main3.ci95.min():.0f}", f"{100 * main3.ci95.max():.0f}"
    NUM["CiNewSmallLo"] = f"{100 * main3.ci95_small.min():.0f}"
    NUM["CiNewBigLo"], NUM["CiNewBigHi"] = f"{100 * main3.ci95_big.min():.0f}", f"{100 * main3.ci95_big.max():.0f}"
    NUM["CiNewZeroLo"] = f"{100 * main3.ci95_zero.min():.0f}"
    NUM["CiFracZeroOneMb"] = f"{100 * cva[(cva.method == CS2) & (cva.mem == 1.0)].frac_zero.iloc[0]:.0f}"
    NUM["CiWidthTen"] = f"{cva[(cva.method == CS2) & (cva.mem == 10.0)].ci_relwidth.iloc[0]:.2f}"
    NUM["CiOneSidedBigLo"] = f"{100 * main3.ci95_f_big.min():.0f}"
    ws = main3.loc[main3.ci95_f_big.idxmin()]
    NUM["CiOneSidedBigWorst"] = DS_LABEL.get(ws.dataset, ws.dataset)
    NUM["CiWidthGain"] = f"{100 * ((main3.ci_relwidth / main3.ci_relwidth_f).median() - 1):.0f}"
    worst = main3.loc[main3.ci95.idxmin()]
    NUM["CiNewWorstStream"] = DS_LABEL.get(worst.dataset, worst.dataset)
    NUM["CiNewWorstMem"] = f"{worst.mem:g}"

# ---------------------------------------------------------------- range predicates
if os.path.exists("results/agg_range.csv"):
    rg_ = pd.read_csv("results/agg_range.csv")
    kc = [c for c in [f"nae4_k{k}" for k in range(1, 5)] if c in rg_.columns]
    rg_["naer"] = rg_[kc].mean(axis=1)
    kh = [c for c in [f"nae4_k{k}_ht" for k in range(1, 5)] if c in rg_.columns]
    rg_["naer_thr"] = rg_["nae4_ht"]     # OmniSketch: the threshold estimator is stored in the HT column
    rm = rg_.groupby(["dataset", "method", "mem"], as_index=False)[["naer", "naer_thr", "nae4", "ns_per_rec"]].mean()
    om = rm[rm.method.str.startswith("omnid")]
    ocode = om.groupby(["dataset", "mem"], as_index=False).naer.min().assign(method="omnid_code")
    othr = om.groupby(["dataset", "mem"], as_index=False).naer_thr.min().rename(columns={"naer_thr": "naer"}).assign(method="omnid_thr")
    rm = pd.concat([rm, ocode, othr], ignore_index=True)
    RP = rm.pivot_table(index=["dataset", "mem"], columns="method", values="naer")
    DY2, DY1, DYR = CS2 + "/dy=1", CS1 + "/dy=1", RSC2 + "/dy=1"
    rlab = [("res", "RS"), ("uss", "USS"), (COCO, "CocoSketch"), (USSC2, "USS+C2"), (COCOC2, "Coco+C2"), (CSH, "\\csh"),
            (CS1, "\\cs (C1)"), (CS2, "\\cs"), (RSC2, "RS+C2"),
            (DY1, "\\cs (C1) + dyadic"), (DY2, "\\cs + dyadic"), (DYR, "RS+C2 + dyadic"),
            ("omnid_code", "\\omni dyadic (code est.)"), ("omnid_thr", "\\omni dyadic (thresh.\\ est.)")]
    rcols = [(ds, x) for ds in ("taxi_range", "mawi_range") for x in (1.0, 10.0, 50.0)]
    with open("paper/tables/range.tex", "w") as f:
        f.write("\\begin{tabular}{l" + "r" * len(rcols) + "}\n\\toprule\n")
        f.write(" & \\multicolumn{3}{c}{Taxi} & \\multicolumn{3}{c}{MAWI} \\\\\n")
        f.write("Method & " + " & ".join(f"{x:g} MB" for _, x in rcols) + " \\\\\n\\midrule\n")
        for key, lab in rlab:
            if key not in RP.columns:
                continue
            cells = []
            for ds, x in rcols:
                v = RP[key].get((ds, x), np.nan)
                best = RP.loc[(ds, x)].min() if (ds, x) in RP.index else np.nan
                txt = fmt_sci(v)
                cells.append(f"\\textbf{{{txt}}}" if np.isfinite(v) and v <= best * 1.0001 else txt)
            f.write(lab + " & " + " & ".join(cells) + " \\\\\n")
            if key in (CSH, RSC2):
                f.write("\\addlinespace[2pt]\n")
        f.write("\\bottomrule\n\\end{tabular}\n")
    note("range predicates (NAE, equal weight per k)", RP)
    ours = RP[[c for c in (DY1, DY2, DYR) if c in RP.columns]].min(axis=1)
    obest = RP[["omnid_code", "omnid_thr"]].min(axis=1)
    put("RangeBestVsOmni", obest / ours, 0)
    put("RangeOmniThrVsCode", RP["omnid_code"] / RP["omnid_thr"], 1)
    put("RangeDyGain", RP[CS2] / RP[DY2], 1)
    put("RangeDyGainCalTwo", RP[RSC2] / RP[DYR], 1)
    put("RangeResVsDy", RP["res"] / ours, 1)
    put("RangeCocoCalVsDy", RP[COCOC2] / RP[DY2], 1) if COCOC2 in RP.columns else None
    put("RangeUssVsDy", RP["uss"] / RP[DY2], 1)
    for ds, tag in [("taxi_range", "Taxi"), ("mawi_range", "Mawi")]:
        put("RangeDyGain" + tag, (RP[CS2] / RP[DY2]).loc[ds], 1)
        put("RangeBestVsOmni" + tag, (obest / ours).loc[ds], 0)
        s = RP["omnid_code"].loc[ds]
        if 1.0 in s.index and 20.0 in s.index:
            NUM["RangeOmniCodeGrow" + tag] = fnum(s.loc[20.0] / s.loc[1.0], 2)
    # memory of the dyadic counters and the share of OmniSketch cells (from the info strings)
    import re as _re3
    inf = rg_[(rg_.method == DY2) & (rg_.mem == 10.0)]["info"]
    for ds, tag in [("taxi_range", "Taxi"), ("mawi_range", "Mawi")]:
        inf = rg_[(rg_.method == DY2) & (rg_.mem == 10.0) & (rg_.dataset == ds)]["info"]
        mm = _re3.search(r"dy_mb=([0-9.]+)", inf.iloc[0]) if len(inf) else None
        if mm:
            NUM["RangeDyKb" + tag] = f"{1000 * float(mm.group(1)):.0f}"

# ---------------------------------------------------------------- merging
if os.path.exists("results/agg_merge.csv"):
    mg = pd.read_csv("results/agg_merge.csv")
    mg["parts"] = mg.method.str.extract(r"@(\d+)$")[0].astype(float)    # bench merge names a spec "<spec>@<parts>"
    mg["stream"] = mg["dataset"]
    mg["method"] = mg.method.str.replace(r"@\d+$", "", regex=True)
    mg["nae0m"] = mg[[f"nae0_k{k}" for k in range(2, 7) if f"nae0_k{k}" in mg.columns]].mean(axis=1)
    mgm = mg.groupby(["stream", "method", "mem", "parts"], as_index=False).nae0m.mean()
    base = mgm[mgm.parts == 1].rename(columns={"nae0m": "b"}).drop(columns="parts")
    mgm = mgm.merge(base, on=["stream", "method", "mem"])
    mgm["ratio"] = mgm.nae0m / mgm.b
    # exactness of the OmniSketch merge: identical estimates per seed and query
    ks = [2, 4, 8, 16]
    mlab = [(CS1, "\\cs (C1)"), (CS2, "\\cs"), (RSC2, "RS+C2"), ("omni:256:2", "\\omni $256{\\times}2$"), ("omni:1024:2", "\\omni $1024{\\times}2$")]
    with open("paper/tables/merge.tex", "w") as f:
        f.write("\\begin{tabular}{l" + "r" * (2 * len(ks)) + "}\n\\toprule\n")
        f.write(" & \\multicolumn{4}{c}{1 MB, parts $k$} & \\multicolumn{4}{c}{10 MB, parts $k$} \\\\\n")
        f.write("Method & " + " & ".join(str(k) for k in ks) + " & " + " & ".join(str(k) for k in ks) + " \\\\\n\\midrule\n")
        for key, lab in mlab:
            cells = []
            for x in (1.0, 10.0):
                for k in ks:
                    s = mgm[(mgm.method == key) & (mgm.mem == x) & (mgm.parts == k)]
                    cells.append(f"{geo(s.ratio):.2f}" if len(s) else "--")
            f.write(lab + " & " + " & ".join(cells) + " \\\\\n")
        f.write("\\bottomrule\n\\end{tabular}\n")
    note("merging: NAE(k parts) / NAE(1 part)", mgm.pivot_table(index=["stream", "method", "mem"], columns="parts", values="ratio"))
    ourm = mgm[mgm.method.isin([CS1, CS2, RSC2]) & (mgm.parts > 1)]
    put("MergeCs", ourm.ratio, 2)
    omm = mgm[mgm.method.str.startswith("omni") & (mgm.parts > 1)]
    put("MergeOmni", omm.ratio, 3)
    NUM["MergeStreams"] = str(mgm.stream.nunique())

# ---------------------------------------------------------------- pair-table selection on Census
if os.path.exists("results/agg_pairs.csv"):
    pr = pd.read_csv("results/agg_pairs.csv")
    pr["pp"] = pr.dataset.str.extract(r"census_p(\d+)").astype(int)
    pr["wl"] = np.where(pr.dataset.str.endswith("_tmpl"), "tmpl", "rnd")
    pr["nae0m"] = pr[[f"nae0_k{k}" for k in range(2, 5) if f"nae0_k{k}" in pr.columns]].mean(axis=1)
    pm = pr.groupby(["pp", "wl", "method"], as_index=False)[["nae0m", "nae0_k2", "ns_per_rec"]].mean()
    sel_lab = {"wl": "workload", "mi": "mutual inf.", "sc": "anchor red.", "rnd": "random"}
    with open("paper/tables/pairs.tex", "w") as f:
        f.write("\\begin{tabular}{lrrrrrr}\n\\toprule\n")
        f.write(" & \\multicolumn{3}{c}{random workload} & \\multicolumn{3}{c}{template workload} \\\\\n")
        f.write("Pair tables & $p{=}16$ & $32$ & $68$ & $p{=}16$ & $32$ & $68$ \\\\\n\\midrule\n")
        def cell(pp_, wl_, meth):
            s = pm[(pm.pp == pp_) & (pm.wl == wl_) & (pm.method == meth)]
            return s.nae0m.iloc[0] if len(s) else np.nan
        rows_ = [("none (C1 only)", lambda pp_: CS1), ("all $\\binom p2$", lambda pp_: CS2)]
        for sel in ("wl", "mi", "sc", "rnd"):
            for mult, ml in ((1, "p"), (4, "4p")):
                rows_.append((f"{sel_lab[sel]}, $m={ml}$", lambda pp_, sel=sel, mult=mult: f"{CS2}/psel={sel}/np={pp_ * mult}"))
        for lab, fn in rows_:
            vals = [(pp_, wl_, cell(pp_, wl_, fn(pp_))) for wl_ in ("rnd", "tmpl") for pp_ in (16, 32, 68)]
            cells = []
            for pp_, wl_, v in vals:
                best = pm[(pm.pp == pp_) & (pm.wl == wl_)].nae0m.min()
                txt = fmt_sci(v)
                cells.append(f"\\textbf{{{txt}}}" if np.isfinite(v) and v <= best * 1.0001 else txt)
            f.write(lab + " & " + " & ".join(cells) + " \\\\\n")
            if lab.startswith("all") or lab.endswith("4p$"):
                f.write("\\addlinespace[2pt]\n")
        f.write("\\bottomrule\n\\end{tabular}\n")
    note("pair selection (NAE k>=2)", pm.pivot_table(index="method", columns=["wl", "pp"], values="nae0m"))
    for wl_, tag in (("tmpl", "Tmpl"), ("rnd", "Rnd")):
        for pp_, pt in ((68, "SixtyEight"), (32, "ThirtyTwo"), (16, "Sixteen")):
            g_ = pm[(pm.pp == pp_) & (pm.wl == wl_)].set_index("method").nae0m
            if CS2 not in g_.index:
                continue
            for sel in ("wl", "mi", "sc", "rnd"):
                for mult, ml in ((1, "P"), (4, "FourP")):
                    key = f"{CS2}/psel={sel}/np={pp_ * mult}"
                    if key in g_.index:
                        NUM[f"Pair{tag}{pt}{sel.capitalize()}{ml}"] = fnum(g_[CS2] / g_[key], 2)    # > 1: selection better than all pairs
            if CS1 in g_.index:
                NUM[f"Pair{tag}{pt}AllVsCOne"] = fnum(g_[CS1] / g_[CS2], 2)
    loss = []   # random workload: best selection rule with m = 4p versus all pairs (error ratio, > 1: selection worse)
    for pp_ in (16, 32, 68):
        g_ = pm[(pm.pp == pp_) & (pm.wl == "rnd")].set_index("method").nae0m
        sel4 = [g_.get(f"{CS2}/psel={sel}/np={4 * pp_}", np.nan) for sel in ("wl", "mi", "sc", "rnd")]
        if CS2 in g_.index:
            loss.append(np.nanmin(sel4) / g_[CS2])
    put("PairRndLoss", loss, 1)
    if os.path.exists("results/time/time.log"):
        tt_ = []
        for ln in open("results/time/time.log"):
            if ln.startswith("TIME,census"):
                f_ = ln.strip().split(",")
                tt_.append(dict(dataset=f_[1], method=f_[2], mrps=float(f_[7])))
        tt_ = pd.DataFrame(tt_).drop_duplicates(["dataset", "method"])
        if len(tt_):
            g_ = tt_[tt_.dataset == "census_p68"].set_index("method").mrps
            for key, tag in [(CS1, "COne"), (CS2, "All"), (f"{CS2}/psel=wl/np=68", "SelP"), (f"{CS2}/psel=wl/np=272", "SelFourP")]:
                if key in g_.index:
                    NUM[f"PairTp{tag}"] = f"{g_[key]:.1f}" if g_[key] >= 1 else f"{g_[key]:.2f}"

# ---------------------------------------------------------------- throughput of CocoSketch, dyadic ranges (scripts/run_timing.sh)
if os.path.exists("results/time/time.log"):
    tt2 = []
    for ln in open("results/time/time.log"):
        if ln.startswith("TIME,"):
            f_ = ln.strip().split(",")
            tt2.append(dict(dataset=f_[1], method=f_[2], mem=float(f_[3]), mrps=float(f_[7])))
    tt2 = pd.DataFrame(tt2).drop_duplicates(["dataset", "method", "mem"])
    note("throughput of CocoSketch / range configurations (Mrec/s)", tt2.pivot_table(index=["dataset", "mem"], columns="method", values="mrps"))
    rt = tt2[tt2.dataset.str.endswith("_range")].pivot_table(index=["dataset", "mem"], columns="method", values="mrps")
    if len(rt):
        oc_ = [c for c in rt.columns if c.startswith("omnid")]
        DY2 = CS2 + "/dy=1"
        if DY2 in rt.columns and oc_:
            put("RangeTpVsOmni", rt[DY2] / rt[oc_].max(axis=1), 1)        # > 1: CoStrata + dyadic faster than the fastest dyadic OmniSketch
            put("RangeTpDyCost", rt[CS2] / rt[DY2], 2)
            NUM["RangeTpDyLo"], NUM["RangeTpDyHi"] = f"{rt[DY2].min():.1f}", f"{rt[DY2].max():.1f}"
            NUM["RangeTpOmniLo"], NUM["RangeTpOmniHi"] = f"{rt[oc_].min().min():.1f}", f"{rt[oc_].max().max():.1f}"
    ct = tt2[tt2.dataset.isin(["taxi", "mawi"])].pivot_table(index=["dataset", "mem"], columns="method", values="mrps")
    for key, tag in [(COCO, "Coco"), (COCOC2, "CocoCalTwo"), ("uss", "Uss"), (USSC2, "UssCalTwo"), (CS2, "Cs")]:
        if key in ct.columns:
            NUM[f"TpNew{tag}Lo"], NUM[f"TpNew{tag}Hi"] = f"{ct[key].min():.1f}", f"{ct[key].max():.1f}"

NUM.setdefault("CsThroughput", "??")
with open("paper/numbers.tex", "w") as f:
    f.write("% generated by scripts/plots.py -- do not edit\n")
    for k, v_ in sorted(NUM.items()):
        f.write(f"\\newcommand{{\\Num{k}}}{{{v_}}}\n")
with open("results/summary.txt", "w") as f:
    f.write("".join(NOTES))
print("figures and tables written;", len(NUM), "numbers")
