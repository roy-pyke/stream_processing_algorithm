#!/usr/bin/env python3
"""Produce all figures (PDF) and LaTeX table fragments of the paper from aggregated results.

  plots.py   (reads results/agg_*.csv and results/time/*.log, writes paper/figs and paper/tables)
"""
import os
import re

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

os.makedirs("paper/figs", exist_ok=True)
os.makedirs("paper/tables", exist_ok=True)

# Fixed entity -> (label, colour, marker, linestyle). Colour follows the method in every figure.
STYLE = {
    "cs":       ("CoStrata",               "#2a78d6", "o", "-"),
    "res":      ("Reservoir (RS)",         "#eb6834", "s", "-"),
    "uss":      ("USS",                    "#1baf7a", "^", "-"),
    "ussc":     ("USS+C",                  "#eda100", "v", "-"),
    "omnip":    ("OmniSketch (paper)",     "#e87ba4", "X", "--"),
    "omnibest": ("OmniSketch (best cfg.)", "#008300", "D", "--"),
    "cs:0:1e9": ("RS+C",                   "#6250d6", "P", ":"),
}
OMNI_GRID = ["omnil", "omni:64:2", "omni:256:2", "omni:1024:2", "omni:256:3", "omni:1024:1", "omni:4096:1", "omni:4096:2"]
DS_LABEL = {"zipf10_5M": r"Zipf $\alpha{=}1.0$", "zipf13_5M": r"Zipf $\alpha{=}1.3$", "zipf15_5M": r"Zipf $\alpha{=}1.5$",
            "unif_5M": "Uniform", "taxi": "Taxi", "mawi": "MAWI"}
DATASETS = ["zipf10_5M", "zipf13_5M", "zipf15_5M", "unif_5M", "taxi", "mawi"]
INK, MUTED, GRID = "#1a1a19", "#6b6a63", "#e4e2da"

plt.rcParams.update({
    "font.size": 8, "axes.titlesize": 8.5, "axes.labelsize": 8, "legend.fontsize": 7.5,
    "xtick.labelsize": 7.5, "ytick.labelsize": 7.5, "axes.edgecolor": MUTED, "axes.labelcolor": INK,
    "xtick.color": MUTED, "ytick.color": MUTED, "text.color": INK, "axes.linewidth": 0.6,
    "grid.color": GRID, "grid.linewidth": 0.6, "axes.spines.top": False, "axes.spines.right": False,
    "pdf.fonttype": 42, "font.family": "DejaVu Sans",
})


def mean_over_seeds(df, cols):
    return df.groupby(["dataset", "method", "mem"], as_index=False)[cols].mean()


def with_omnibest(m, cols):
    """Add the pseudo-method 'omnibest': per (dataset, mem) and per metric, the best grid configuration."""
    g = m[m.method.isin(OMNI_GRID)]
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
    e = int(np.floor(np.log10(x)))
    return f"{x / 10 ** e:.1f}e{e}"


acc = pd.read_csv("results/agg_acc.csv")
# Primary metric: NAE on record-sampled queries with k >= 2 predicates (the workload definition of the OmniSketch paper).
acc["nae0m"] = acc[[f"nae0_k{k}" for k in range(2, 7)]].mean(axis=1)
dcols = [f"mre1_d{j}" for j in range(1, 7)]
kcols = [f"mre0_k{k}" for k in range(1, 7)]
nkcols = [f"nae0_k{k}" for k in range(1, 7)]
metrics = ["nae0", "nae0m", "nae1", "mre0", "mre1"] + dcols + kcols + nkcols
m = with_omnibest(mean_over_seeds(acc, metrics), metrics)

# ------------------------------------------------------------------ Fig: NAE vs memory
order = ["omnip", "omnibest", "res", "cs:0:1e9", "uss", "ussc", "cs"]
fig, axes = plt.subplots(2, 3, figsize=(7.0, 4.1), sharex=True)
for ax, ds in zip(axes.flat, DATASETS):
    for key in order:
        s = m[(m.dataset == ds) & (m.method == key)].sort_values("mem")
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

# ------------------------------------------------------------------ Fig: relative error by anchor frequency
order2 = ["omnibest", "res", "cs:0:1e9", "uss", "cs"]
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
fig.legend(h, l, loc="upper center", ncol=5, frameon=False, bbox_to_anchor=(0.5, 1.1))
fig.tight_layout()
fig.savefig("paper/figs/anchor.pdf", bbox_inches="tight")
plt.close(fig)

# ------------------------------------------------------------------ Table: main accuracy table at 10 MB
tab_methods = ["omnip", "omnibest", "res", "uss", "cs:0:1e9", "ussc", "cs"]
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
            f.write(STYLE[key][0].replace("(", "{\\footnotesize(").replace(")", ")}") + " & " + " & ".join(cells) + " \\\\\n")
        if mem == 1.0:
            f.write("\\midrule\n")
    f.write("\\bottomrule\n\\end{tabular}\n")

# ratios quoted in the text
with open("paper/tables/ratios.txt", "w") as f:
    for mem in sorted(m.mem.unique()):
        for ds in DATASETS:
            s = m[(m.dataset == ds) & (m.mem == mem)].set_index("method")
            if "res" not in s.index:
                continue
            r = lambda a, b, c="nae0": s.loc[a, c] / s.loc[b, c] if a in s.index and b in s.index else float("nan")
            f.write(f"{ds:10s} mem={mem:5g}  omnip/res={r('omnip', 'res'):7.1f}  omnibest/res={r('omnibest', 'res'):6.1f}  "
                    f"res/cs={r('res', 'cs'):5.2f}  res/rsc={r('res', 'cs:0:1e9'):5.2f}  uss/ussc={r('uss', 'ussc'):5.2f}  "
                    f"cs/ussc={r('cs', 'ussc'):5.2f}  cs/rsc={r('cs', 'cs:0:1e9'):5.2f} | k1: omnibest/res={r('omnibest', 'res', 'nae1'):6.1f} "
                    f"res/cs={r('res', 'cs', 'nae1'):5.2f} | d1 res/cs={r('res', 'cs', 'mre1_d1'):6.2f} d2 res/cs={r('res', 'cs', 'mre1_d2'):6.2f} "
                    f"d1 rsc/cs={r('cs:0:1e9', 'cs', 'mre1_d1'):6.2f} d1 omnibest/cs={r('omnibest', 'cs', 'mre1_d1'):6.2f}\n")
    # which OmniSketch configuration wins
    g = mean_over_seeds(acc[acc.method.isin(OMNI_GRID)], ["nae0"])
    idx = g.groupby(["dataset", "mem"]).nae0.idxmin()
    f.write("\nbest OmniSketch configuration per (dataset, mem):\n" + g.loc[idx].to_string() + "\n")
    cs = acc[acc.method == "cs"]
    f.write("\nCoStrata 95% CI coverage (all queries / truth>=100), HT bias:\n")
    f.write(cs.groupby(["dataset", "mem"])[["cov95", "cov95_big", "bias_ht"]].mean().to_string() + "\n")

# ------------------------------------------------------------------ Table: NAE by number of predicates (10 MB)
with open("paper/tables/by_k.tex", "w") as f:
    f.write("\\begin{tabular}{l" + "r" * 12 + "}\n\\toprule\n")
    f.write(" & \\multicolumn{6}{c}{Taxi, number of predicates $k$} & \\multicolumn{6}{c}{MAWI, number of predicates $k$} \\\\\n")
    f.write("Method & " + " & ".join(str(k) for k in range(1, 7)) + " & " + " & ".join(str(k) for k in range(1, 7)) + " \\\\\n\\midrule\n")
    for key in tab_methods:
        cells = []
        for ds in ("taxi", "mawi"):
            s = m[(m.dataset == ds) & (m.mem == 10.0)]
            v = s[s.method == key][nkcols].to_numpy(float)[0]
            best = s[s.method.isin(tab_methods)][nkcols].min().to_numpy(float)
            cells += [f"\\textbf{{{fmt_sci(x)}}}" if x <= b * 1.0001 else fmt_sci(x) for x, b in zip(v, best)]
        f.write(STYLE[key][0].replace("OmniSketch", "Omni.").replace("Reservoir (RS)", "RS") + " & " + " & ".join(cells) + " \\\\\n")
    f.write("\\bottomrule\n\\end{tabular}\n")

# ------------------------------------------------------------------ numbers quoted in the text (numbers.tex)
PA = m.pivot_table(index=["dataset", "mem"], columns="method", values="nae0")    # all k
P = m.pivot_table(index=["dataset", "mem"], columns="method", values="nae0m")    # k >= 2 (primary)
gmean = lambda r: float(np.exp(np.log(r).mean()))
def rng(r, d=1):
    f = lambda x: f"{x:.0f}" if x >= 100 else f"{x:.{d}f}"
    return f(r.min()), f(r.max()), f(gmean(r))
NUM = {}
def put(tag, r, d=1):
    NUM[tag + "Lo"], NUM[tag + "Hi"], NUM[tag + "Gm"] = rng(r, d)
nonunif = [i for i in P.index if i[0] != "unif_5M"]
put("OmniPaperVsRes", P["omnip"] / P["res"]);              put("OmniPaperVsResAll", PA["omnip"] / PA["res"], 0)
put("OmniBestVsRes", (P["omnibest"] / P["res"]).loc[nonunif]); put("OmniBestVsResAll", PA["omnibest"] / PA["res"])
put("OmniBestVsResUnif", (P["omnibest"] / P["res"]).loc["unif_5M"], 2)
put("OmniBestVsCs", (P["omnibest"] / P["cs"]).loc[nonunif])
put("CalibRs", P["res"] / P["cs:0:1e9"], 2);               put("CalibRsAll", PA["res"] / PA["cs:0:1e9"], 2)
put("CalibUss", P["uss"] / P["ussc"], 2);                  put("CalibUssAll", PA["uss"] / PA["ussc"], 2)
put("CsVsRsc", P["cs"] / P["cs:0:1e9"], 2);                put("CsVsRscAll", PA["cs"] / PA["cs:0:1e9"], 2)
put("CsVsRes", P["cs"] / P["res"], 2)
for ds, tag in [("mawi", "Mawi"), ("taxi", "Taxi")]:
    put("CalibRs" + tag, (P["res"] / P["cs:0:1e9"]).loc[ds], 2)
    put("UssVsRes" + tag, (P["res"] / P["uss"]).loc[ds], 1)
    put("OmniBestVsRes" + tag, (P["omnibest"] / P["res"]).loc[ds], 1)
syn = [i for i in P.index if i[0].startswith(("zipf", "unif"))]
put("CalibRsSyn", (P["res"] / P["cs:0:1e9"]).loc[syn], 2)
k1 = m.pivot_table(index=["dataset", "mem"], columns="method", values="nae0_k1")
put("CalibRsKone", k1["res"] / k1["cs:0:1e9"], 0)
put("OmniPaperFlat", P["omnip"].loc["unif_5M"] * 1e6, 2)
for ds, tag in [("zipf10_5M", "ZipfTen"), ("zipf13_5M", "ZipfThirteen"), ("zipf15_5M", "ZipfFifteen"), ("unif_5M", "Unif")]:
    v = m[(m.dataset == ds) & (m.mem == 10.0) & (m.method == "omnip")].nae0_k3.iloc[0]
    NUM["OmniPaperKthree" + tag] = fmt_sci(v)
for mem, tag in [(5.0, "Five"), (10.0, "Ten"), (20.0, "Twenty"), (50.0, "Fifty")]:
    t = m[(m.dataset == "taxi") & (m.mem == mem)].set_index("method")
    others = t.loc[["omnibest", "res", "cs:0:1e9", "uss", "ussc"]]
    for d, dn in [(1, "One"), (2, "Two")]:
        NUM[f"TaxiD{dn}Cs{tag}"] = f"{t.loc['cs', f'mre1_d{d}']:.2g}"
        NUM[f"TaxiD{dn}Other{tag}"] = f"{others[f'mre1_d{d}'].min():.2g}"
        NUM[f"TaxiD{dn}Gain{tag}"] = f"{others[f'mre1_d{d}'].min() / t.loc['cs', f'mre1_d{d}']:.0f}" if others[f'mre1_d{d}'].min() / t.loc['cs', f'mre1_d{d}'] >= 10 else f"{others[f'mre1_d{d}'].min() / t.loc['cs', f'mre1_d{d}']:.1f}"
tt = m[m.dataset == "taxi"].pivot_table(index="mem", columns="method", values=["mre1_d4", "mre1_d5"])
rr = pd.concat([tt["mre1_d4"]["cs"] / tt["mre1_d4"]["res"], tt["mre1_d5"]["cs"] / tt["mre1_d5"]["res"]])
NUM["TaxiFreqPenaltyLo"], NUM["TaxiFreqPenaltyHi"] = f"{rr.min():.1f}", f"{rr.max():.1f}"
cst = acc[acc.method == "cs"]
NUM["HtBiasMax"] = f"{100 * cst.groupby(['dataset', 'mem']).bias_ht.mean().abs().max():.1f}"
NUM["HtBiasMaxTen"] = f"{100 * cst[cst.mem >= 10].groupby(['dataset', 'mem']).bias_ht.mean().abs().max():.2f}"

# ------------------------------------------------------------------ Table: CI coverage
cs = acc[acc.method == "cs"].groupby(["dataset", "mem"], as_index=False)[["cov95", "cov95_big"]].mean()
with open("paper/tables/coverage.tex", "w") as f:
    mems = [1.0, 10.0, 50.0]
    f.write("\\begin{tabular}{l" + "rr" * len(mems) + "}\n\\toprule\n")
    f.write(" & " + " & ".join(f"\\multicolumn{{2}}{{c}}{{{x:g} MB}}" for x in mems) + " \\\\\n")
    f.write("Dataset & " + " & ".join("all & $f{\\ge}100$" for _ in mems) + " \\\\\n\\midrule\n")
    for ds in DATASETS:
        cells = []
        for x in mems:
            s = cs[(cs.dataset == ds) & (cs.mem == x)]
            cells += [f"{100 * s.cov95.iloc[0]:.1f}", f"{100 * s.cov95_big.iloc[0]:.1f}"] if len(s) else ["--", "--"]
        f.write(DS_LABEL[ds] + " & " + " & ".join(cells) + " \\\\\n")
    f.write("\\bottomrule\n\\end{tabular}\n")

# ------------------------------------------------------------------ Timing table + figure
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
    tm = ["res", "cs", "uss", "ussc", "omnip", "omnil", "omni:256:2", "omni:1024:2"]
    names = {"res": "Reservoir (RS)", "cs": "CoStrata", "uss": "USS", "ussc": "USS+C", "omnip": "OmniSketch (paper)",
             "omnil": "OmniSketch (paper $w,d$; lean)", "omni:256:2": "OmniSketch ($w{=}256,d{=}2$)",
             "omni:1024:2": "OmniSketch ($w{=}1024,d{=}2$)"}
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
                cells.append(f"{s.mrps.iloc[0]:.1f}" if len(s) else "--")
            f.write(names[key] + " & " + " & ".join(cells) + " \\\\\n")
        if len(q):
            f.write("\\midrule\n")
            f.write(f"\\multicolumn{{{len(cols) + 1}}}{{l}}{{\\emph{{Query latency ($\\mu$s per query)}}}} \\\\\n")
            for key in tm:
                cells = []
                for ds, mem in cols:
                    s = q[(q.dataset == ds) & (q.mem == mem) & (q.method == key)]
                    cells.append((f"{s.q_us.iloc[0]:.1f}" if s.q_us.iloc[0] < 10 else f"{s.q_us.iloc[0]:.0f}") if len(s) else "--")
                if any(c != "--" for c in cells):
                    f.write(names[key] + " & " + " & ".join(cells) + " \\\\\n")
        f.write("\\bottomrule\n\\end{tabular}\n")

    # attribute scaling on census
    if os.path.exists("results/agg_census.csv"):
        c = pd.read_csv("results/agg_census.csv")
        c["pp"] = c.dataset.str.extract(r"_p(\d+)").astype(int)
        c["nae0"] = c[[f"nae0_k{k}" for k in range(2, 5)]].mean(axis=1)   # k >= 2, as in the main results
        cm = c.groupby(["pp", "method"], as_index=False)[["nae0"]].mean()
        cg = cm[cm.method.isin(OMNI_GRID)].groupby("pp", as_index=False).nae0.min()
        cg["method"] = "omnibest"
        cm = pd.concat([cm, cg], ignore_index=True)
        tc = t[t.dataset.str.startswith("census_p")].copy()
        tc["pp"] = tc.dataset.str.extract(r"_p(\d+)").astype(int)
        fig, axes = plt.subplots(1, 2, figsize=(7.0, 2.3))
        for key in ["omnip", "omnibest", "res", "uss", "cs"]:
            s = cm[cm.method == key].sort_values("pp")
            line(axes[0], s.pp, s.nae0, key)
        axes[0].set_xscale("log", base=2); axes[0].set_yscale("log")
        axes[0].set_ylabel("normalized abs. error"); axes[0].set_title("Accuracy at 10 MB")
        for key, tk in [("omnip", "omnip"), ("omnibest", "omni:256:2"), ("res", "res"), ("uss", "uss"), ("cs", "cs")]:
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
        h, l = axes[0].get_legend_handles_labels()
        fig.legend(h, l, loc="upper center", ncol=5, frameon=False, bbox_to_anchor=(0.5, 1.1))
        fig.tight_layout()
        fig.savefig("paper/figs/census_scaling.pdf", bbox_inches="tight")
        plt.close(fig)
        with open("paper/tables/ratios.txt", "a") as f:
            f.write("\ncensus scaling (nae0 by p):\n" + cm.pivot(index="pp", columns="method", values="nae0").to_string() + "\n")
            if len(tc):
                f.write("\ncensus throughput (Mrec/s):\n" + tc.pivot(index="pp", columns="method", values="mrps").to_string() + "\n")
                tp = tc.pivot(index="pp", columns="method", values="mrps")
                for meth, tag in [("cs", "Cs"), ("res", "Res"), ("omnip", "Omni"), ("uss", "Uss")]:
                    NUM[f"CensusTp{tag}Four"] = f"{tp.loc[4, meth]:.1f}"
                    NUM[f"CensusTp{tag}SixtyEight"] = f"{tp.loc[68, meth]:.2f}"
        for pp, tag in [(4, "Four"), (68, "SixtyEight")]:
            row = cm[cm.pp == pp].set_index("method").nae0
            NUM[f"CensusOmniBestVsRes{tag}"] = f"{row['omnibest'] / row['res']:.1f}"
            NUM[f"CensusResVsCs{tag}"] = f"{row['res'] / row['cs']:.1f}"
    tp = t[t.dataset.isin(["taxi", "mawi"])].pivot_table(index=["dataset", "mem"], columns="method", values="mrps")  # as in the table
    if "cs" in tp.columns and len(tp.dropna()) > 0:
        for meth, tag in [("cs", "Cs"), ("res", "Res"), ("uss", "Uss"), ("omnip", "OmniPaper"), ("omnil", "OmniLean"), ("omni:256:2", "OmniMid")]:
            NUM[f"Tp{tag}Lo"], NUM[f"Tp{tag}Hi"] = f"{tp[meth].min():.1f}", f"{tp[meth].max():.1f}"
        NUM["CsThroughput"] = f"{tp['cs'].min():.0f}--{tp['cs'].max():.0f}"
        r = tp["cs"] / tp[["omnip", "omnil", "omni:256:2", "omni:1024:2"]].max(axis=1)
        NUM["TpCsVsOmniLo"], NUM["TpCsVsOmniHi"] = f"{r.min():.1f}", f"{r.max():.0f}"
        r = tp["res"] / tp["cs"]
        NUM["TpResVsCsLo"], NUM["TpResVsCsHi"] = f"{r.min():.1f}", f"{r.max():.0f}"
    if len(q):
        ql = q.pivot_table(index=["dataset", "mem"], columns="method", values="q_us")
        r = ql["cs"] / ql[["omnip", "omni:256:2", "omni:1024:2"]].min(axis=1)
        NUM["QlatCsVsOmniLo"], NUM["QlatCsVsOmniHi"] = f"{r.min():.0f}", f"{r.max():.0f}"
        r = ql["cs"] / ql["omnip"]
        NUM["QlatCsVsOmniPaperLo"], NUM["QlatCsVsOmniPaperHi"] = f"{r.min():.1f}", f"{r.max():.1f}"
        NUM["QlatCsLo"], NUM["QlatCsHi"] = f"{ql['cs'].min():.0f}", f"{ql['cs'].max():.0f}"
        NUM["QlatOmniLo"], NUM["QlatOmniHi"] = f"{ql[['omnip', 'omni:256:2', 'omni:1024:2']].min().min():.1f}", f"{ql[['omnip', 'omni:256:2', 'omni:1024:2']].max().max():.0f}"
        r = ql["res"] / ql["cs"]
        NUM["QlatResVsCsLo"], NUM["QlatResVsCsHi"] = f"{r.min():.1f}", f"{r.max():.1f}"
    with open("paper/tables/ratios.txt", "a") as f:
        f.write("\nthroughput (Mrec/s):\n" + t[~t.dataset.str.startswith("census")].pivot_table(index=["dataset", "mem"], columns="method", values="mrps").to_string() + "\n")
        if len(q):
            f.write("\nquery latency (us):\n" + q.pivot_table(index=["dataset", "mem"], columns="method", values="q_us").to_string() + "\n")

# ------------------------------------------------------------------ Ablation table (10 MB)
if os.path.exists("results/agg_abl.csv"):
    a = pd.read_csv("results/agg_abl.csv")
    a["nae0m"] = a[[f"nae0_k{k}" for k in range(2, 7)]].mean(axis=1)
    am = a.groupby(["dataset", "method", "mem"], as_index=False)[["nae0", "nae0m", "nae0_ht", "mre1_d1", "mre1_d2"]].mean()
    rows = [("\\emph{house weight} $\\lambda$ \\emph{(automatic} $w$\\emph{)}", None),
            ("$\\lambda=0$ (strata only)", "cs:0:0:1"), ("$\\lambda=1/4$", "cs:0:0.25:1"), ("$\\lambda=1$ (default)", "cs:0:1:1"),
            ("$\\lambda=4$", "cs:0:4:1"), ("$\\lambda=\\infty$ (RS+C)", "cs:0:1e9:1"),
            ("\\emph{allocation exponent} $\\gamma$ \\emph{(}$\\lambda=0$\\emph{)}", None),
            ("$\\gamma=0.75$", "cs:0:0:0.75"), ("$\\gamma=0.5$", "cs:0:0:0.5"),
            ("\\emph{buckets per attribute} $w$ \\emph{(}$\\lambda=1$\\emph{)}", None),
            ("$w=64$", "cs:64:1:1"), ("$w=256$", "cs:256:1:1"), ("$w=1024$", "cs:1024:1:1"), ("$w=4096$", "cs:4096:1:1"),
            ("$w=16384$", "cs:16384:1:1"),
            ("\\emph{reference}", None), ("Reservoir (RS)", "res")]
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
                cells += [fmt_sci(r.nae0m.iloc[0]), f"{r.mre1_d1.iloc[0]:.2g}", f"{r.mre1_d2.iloc[0]:.2g}"]
            f.write("\\quad " + label + " & " + " & ".join(cells) + " \\\\\n")
        f.write("\\bottomrule\n\\end{tabular}\n")
    for ds, tag in [("taxi", "Taxi"), ("mawi", "Mawi"), ("zipf15_5M", "Zipf")]:
        r = am[(am.dataset == ds) & (am.mem == 10.0) & (am.method == "cs:0:1:1")]
        NUM["HtVsCal" + tag] = f"{r.nae0_ht.iloc[0] / r.nae0.iloc[0]:.1f}"
    with open("paper/tables/ratios.txt", "a") as f:
        f.write("\nablation (mean over seeds):\n" + am.to_string() + "\n")

if os.path.exists("results/bias.json"):
    import json
    bj = json.load(open("results/bias.json"))
    for ds, tag in [("taxi", "Taxi"), ("zipf15_5M", "Zipf")]:
        if ds in bj:
            b = bj[ds]
            NUM[f"BiasSeeds{tag}"] = str(b["seeds"])
            NUM[f"BiasQueries{tag}"] = str(b["queries"])
            NUM[f"BiasHtTotal{tag}"] = f"{100 * b['ht']['rel_bias_total']:+.2f}"
            NUM[f"BiasHtTotalSe{tag}"] = f"{100 * b['ht']['rel_bias_total_se']:.2f}"
            NUM[f"BiasHtZ{tag}"] = f"{100 * b['ht']['frac_absz_lt_1_96']:.0f}"
            NUM[f"BiasCalTotal{tag}"] = f"{100 * b['cal']['rel_bias_total']:+.2f}"
            NUM[f"BiasCalTotalSe{tag}"] = f"{100 * b['cal']['rel_bias_total_se']:.2f}"
            NUM[f"BiasCalMean{tag}"] = f"{100 * b['cal']['mean_rel_bias']:+.2f}"
            NUM[f"BiasHtMean{tag}"] = f"{100 * b['ht']['mean_rel_bias']:+.2f}"

NUM.setdefault("CsThroughput", "??")
with open("paper/numbers.tex", "w") as f:
    f.write("% generated by scripts/plots.py -- do not edit\n")
    for k, v in sorted(NUM.items()):
        f.write(f"\\newcommand{{\\Num{k}}}{{{v}}}\n")
print("figures and tables written")
