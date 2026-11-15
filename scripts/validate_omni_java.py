#!/usr/bin/env python3
"""Compare the authors' Java OmniSketch with our C++ code, query by query.

Inputs (written by scripts/validate_omni_java.sh):
  results/validation/java/<ds>.<rep>.csv  name,mem,rep,q,truth,estimate,case2Estimate,intersectSize,thrm33Case2
  results/validation/cpp/<ds>.<rep>.csv   name,spec,mem,seed,q,est,ht,sd,cm
     omnij : est = Theorem 3.3 estimate, ht = always-scaled, sd = |S_cap|, cm = n_max
     omnip / omnil : est = always-scaled, ht = Theorem 3.3 estimate
     res   : est
Outputs:
  results/validation/summary.csv      one row per (dataset, mem, rep)
  results/validation/summary_agg.csv  mean over repetitions
  paper/tables/validation.tex
"""
import glob
import os
import sys

import numpy as np
import pandas as pd

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
V = os.path.join(ROOT, "results", "validation")

DS_LABEL = {
    "jzipf10_5M": r"Zipf 1.0$^\dagger$", "jzipf13_5M": r"Zipf 1.3$^\dagger$", "jzipf15_5M": r"Zipf 1.5$^\dagger$",
    "zipf10_5M": "Zipf 1.0", "zipf13_5M": "Zipf 1.3", "zipf15_5M": "Zipf 1.5", "unif_5M": "Uniform",
    "taxi": "Taxi", "mawi": "MAWI",
}


def load_workload(ds):
    rows = []
    with open(os.path.join(ROOT, "work", ds + ".wl")) as f:
        for i, line in enumerate(f):
            t = line.split()
            rows.append((i, int(t[0]), int(t[1]), int(t[2])))
    return pd.DataFrame(rows, columns=["q", "kind", "k", "truth"])


def main():
    out = []
    for jf in sorted(glob.glob(os.path.join(V, "java", "*.csv"))):
        base = os.path.basename(jf)[:-4]
        ds, rep = base.rsplit(".", 1)
        cf = os.path.join(V, "cpp", base + ".csv")
        if not os.path.exists(cf):
            continue
        j = pd.read_csv(jf, header=None, names=["name", "mem", "rep", "q", "truth", "est", "case2", "inter", "c2flag"])
        c = pd.read_csv(cf, header=None, names=["name", "spec", "mem", "seed", "q", "est", "ht", "sd", "cm"])
        wl = load_workload(ds)
        n = {"taxi": 7696617, "mawi": 11331914}.get(ds, 5000000)
        nline = [l for l in open(os.path.join(V, "cpp", base + ".log")) if l.startswith("RUN")]
        if nline:
            n = int(nline[0].split(",")[6])
        for mem, jm in j.groupby("mem"):
            cm = c[c.mem == mem]
            cj = cm[cm.spec == "omnij"].set_index("q").sort_index()
            jm = jm.set_index("q").sort_index()
            if len(cj) != len(jm):
                print("size mismatch", ds, rep, mem, file=sys.stderr)
                continue
            w = wl.set_index("q").loc[jm.index]
            assert (w.truth.values == jm.truth.values).all()
            same_thm = int((cj.est.values == jm.est.values).sum())
            same_c2 = int((cj.ht.values == jm.case2.values).sum())
            same_int = int((cj.sd.values == jm.inter.values).sum())
            # the workload of the original paper: record-sampled queries with >= 2 predicates
            sel = ((w.kind == 0) & (w.k >= 2)).values
            truth = jm.truth.values.astype(float)

            def nae(est):
                return float(np.mean(np.abs(est[sel] - truth[sel])) / n)

            row = dict(dataset=ds, mem=mem, rep=int(rep), queries=len(jm), identical_thm33=same_thm,
                       identical_scaled=same_c2, identical_intersection=same_int,
                       max_abs_diff_thm33=float(np.max(np.abs(cj.est.values - jm.est.values))),
                       max_abs_diff_scaled=float(np.max(np.abs(cj.ht.values - jm.case2.values))),
                       frac_bound_returned=float(np.mean(jm.c2flag.values[sel] == 0)),
                       nae_java_thm33=nae(jm.est.values.astype(float)),
                       nae_java_scaled=nae(jm.case2.values.astype(float)))
            for spec, col in (("omnip", "est"), ("omnil", "est"), ("res", "est")):
                x = cm[cm.spec == spec].set_index("q").sort_index()
                row["nae_" + spec] = nae(x[col].values.astype(float)) if len(x) == len(jm) else float("nan")
            x = cm[cm.spec == "omnip"].set_index("q").sort_index()
            row["nae_omnip_thm33"] = nae(x["ht"].values.astype(float)) if len(x) == len(jm) else float("nan")
            # three-predicate queries (the setting of Fig. 6 of the OmniSketch paper), scaled estimate
            s3 = ((w.kind == 0) & (w.k == 3)).values
            row["nae_java_scaled_k3"] = float(np.mean(np.abs(jm.case2.values[s3] - truth[s3])) / n)
            # share of the stream in the heaviest cell touched by any query (n_max / N)
            row["max_cell_share"] = float(cj.cm.values.max() / n)
            jr = [l for l in open(os.path.join(V, "java", base + ".log")) if l.startswith("JRUN")]
            used = {float(l.split(",")[2]): float(l.split(",")[9]) for l in jr}
            row["java_used_mb"] = used.get(float(mem), float("nan"))
            nom = [l for l in nline if l.split(",")[2] == "omnip" and float(l.split(",")[3]) == float(mem)]
            row["nominal_mb"] = float(nom[0].split(",")[10]) if nom else float("nan")
            out.append(row)
    df = pd.DataFrame(out).sort_values(["dataset", "mem", "rep"])
    df.to_csv(os.path.join(V, "summary.csv"), index=False)
    tot = df[["queries", "identical_thm33", "identical_scaled", "identical_intersection"]].sum()
    print("runs:", len(df), " queries:", int(tot.queries))
    print("identical Theorem-3.3 estimates :", int(tot.identical_thm33))
    print("identical scaled estimates      :", int(tot.identical_scaled))
    print("identical intersection sizes    :", int(tot.identical_intersection))
    print("max |diff| thm33 / scaled       :", df.max_abs_diff_thm33.max(), df.max_abs_diff_scaled.max())
    agg = df.groupby(["dataset", "mem"]).agg(
        frac_bound=("frac_bound_returned", "mean"), java_thm33=("nae_java_thm33", "mean"),
        java_scaled=("nae_java_scaled", "mean"), ours_paper_cfg=("nae_omnip", "mean"),
        ours_paper_cfg_thm33=("nae_omnip_thm33", "mean"), ours_lean=("nae_omnil", "mean"),
        reservoir=("nae_res", "mean")).reset_index()
    agg.to_csv(os.path.join(V, "summary_agg.csv"), index=False)
    # numbers quoted in the paper (read by scripts/plots.py)
    import json
    k3 = df.groupby(["dataset", "mem"]).nae_java_scaled_k3.mean()
    share = df.groupby("dataset").max_cell_share.max()
    a = agg.set_index(["dataset", "mem"])
    num = dict(
        runs=int(len(df)), queries=int(tot.queries), streams=int(df.dataset.nunique()),
        all_identical=bool(tot.queries == tot.identical_thm33 == tot.identical_scaled == tot.identical_intersection),
        bound_frac_min=float(df.frac_bound_returned.min()), bound_frac_max=float(df.frac_bound_returned.max()),
        thm_nae_min=float(agg.java_thm33.min()), thm_nae_max=float(agg.java_thm33.max()),
        k3={f"{d}@{m:g}": float(v) for (d, m), v in k3.items()},
        max_cell_share={d: float(v) for d, v in share.items()},
        java_used_mb_10={d: float(v) for d, v in df[df.mem == 10].groupby("dataset").java_used_mb.mean().items()},
        nominal_mb_10={d: float(v) for d, v in df[df.mem == 10].groupby("dataset").nominal_mb.mean().items()},
        # ratio > 1: our implementation (same B) is more accurate than the authors' code
        java_over_ours={f"{d}@{m:g}": float(a.loc[(d, m), "java_scaled"] / a.loc[(d, m), "ours_paper_cfg"]) for d, m in a.index},
        java_over_res={f"{d}@{m:g}": float(a.loc[(d, m), "java_scaled"] / a.loc[(d, m), "reservoir"])
                       for d, m in a.index if a.loc[(d, m), "reservoir"] > 0})
    json.dump(num, open(os.path.join(V, "numbers.json"), "w"), indent=1)
    pd.set_option("display.width", 200)
    print(agg.to_string(index=False, float_format=lambda v: "%.3g" % v))
    with open(os.path.join(V, "identity.txt"), "w") as f:
        f.write("runs %d queries %d identical_thm33 %d identical_scaled %d identical_intersection %d\n" % (
            len(df), tot.queries, tot.identical_thm33, tot.identical_scaled, tot.identical_intersection))

    # LaTeX table (10 MB and 100 MB)
    def fmt(v):
        if not np.isfinite(v):
            return "--"
        if v == 0:
            return "exact"
        e = int(np.floor(np.log10(v)))
        return r"$%.1f\!\cdot\!10^{%d}$" % (v / 10 ** e, e)

    order = [d for d in DS_LABEL if d in set(agg.dataset)]
    lines = [r"\begin{tabular}{llrrrrr}", r"\toprule",
             r" & & \multicolumn{2}{c}{authors' Java code} & \multicolumn{2}{c}{our \omni} & \\",
             r"\cmidrule(lr){3-4}\cmidrule(lr){5-6}",
             r"Stream & MB & Thm.~3.3 & scaled & same $B$ & lean & RS \\", r"\midrule"]
    for ds in order:
        for mem in (10, 50):
            r = agg[(agg.dataset == ds) & (agg.mem == mem)]
            if r.empty:
                continue
            r = r.iloc[0]
            lines.append("%s & %g & %s & %s & %s & %s & %s \\\\" % (
                DS_LABEL[ds], mem, fmt(r.java_thm33), fmt(r.java_scaled), fmt(r.ours_paper_cfg), fmt(r.ours_lean),
                fmt(r.reservoir)))
    lines += [r"\bottomrule", r"\end{tabular}"]
    os.makedirs(os.path.join(ROOT, "paper", "tables"), exist_ok=True)
    with open(os.path.join(ROOT, "paper", "tables", "validation.tex"), "w") as f:
        f.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
