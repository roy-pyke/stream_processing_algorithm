#!/usr/bin/env python3
"""Aggregate per-query benchmark output into per-run metrics.

  aggregate.py <results_dir> <out.csv> [workload_dir]

For every (dataset, method, mem, seed) run we compute, per workload family:
  nae      mean |est - truth| / N            (the metric of the OmniSketch paper)
  mre      mean |est - truth| / truth
  mre_d{j} mean relative error over queries whose rarest predicate has frequency in [10^j, 10^(j+1))
  mre_k{k} mean relative error over queries with k predicates
  cov95    (CoStrata only) fraction of queries with |est - truth| <= 1.96 * sd
  ci95     (CoStrata, CoStrata-H; runs with interval diagnostics) coverage of the intervals of
           Section 3.4. The calibrated estimate splits the anchor cell exactly, f_hat + g_hat = n_c,
           where g counts its non-matching records, and the interval is built on the rarer part:
           if f_hat <= n_c / 2 (or no sampled row matches), the rule of three
           [0, min(n_c, ln(20) / pi_q)] when m = 0 and otherwise the Fay-Feuer gamma interval
           around f_hat; if f_hat > n_c / 2, the same construction for g (rule of three with the
           inclusion-probability floor pi_a of the anchor's records when every sampled anchor row
           matches), mapped back by f = n_c - g. Truncated to [m, n_c] (m = matching sampled rows,
           n_c = anchor count, both hard bounds). Runs without the anchor columns (wa, wanc) get
           the one-sided version (the matches only).
  ci95_f   the one-sided version, for comparison
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


def _ff(y, v, wm):
    """Fay-Feuer gamma interval of a non-negative weighted total y with variance v; wm is the
    weight of one more contributing record (upper bound only)."""
    y, v, wm = np.asarray(y, float), np.maximum(np.asarray(v, float), 1e-12), np.asarray(wm, float)
    with np.errstate(divide="ignore", invalid="ignore"):
        lo = _gamma.ppf(ALPHA / 2, a=y * y / v, scale=v / np.maximum(y, 1e-300))
        hi = _gamma.ppf(1 - ALPHA / 2, a=(y + wm) ** 2 / (v + wm * wm), scale=(v + wm * wm) / (y + wm))
    lo = np.where(y > 0, np.where(np.isfinite(lo), lo, y), 0.0)
    return lo, np.where(np.isfinite(hi), hi, np.inf)


def _one_side(y, sd, k, w, scale):
    """Interval for a part of the anchor cell estimated by y from k >= 0 sampled records: rule of
    three when k = 0 (each record of the part is sampled with probability >= 1 / w), else Fay-Feuer."""
    lo, hi = _ff(y, sd ** 2, scale * w)
    z = k == 0
    return np.where(z, 0.0, lo), np.where(z, np.log(1 / ALPHA) * w, hi)


def small_count_ci(est, ht, sd, m, wq, nb, wa=None, wanc=None):
    """95% intervals for calibrated estimates (see the module docstring); returns lo, hi arrays."""
    est, ht, sd, m, wq, nb = (np.asarray(x, float) for x in (est, ht, sd, m, wq, nb))
    scale = np.where(ht > 0, est / np.maximum(ht, 1e-300), 1.0)   # n_c / W_a
    lo, hi = _one_side(est, sd, m, wq, scale)
    if wa is not None and wanc is not None and not np.isnan(np.asarray(wanc, float)).any():
        wa, wanc = np.asarray(wa, float), np.asarray(wanc, float)
        g = np.maximum(nb - est, 0.0)
        allmatch = g <= 1e-9 * np.maximum(nb, 1.0)                # no non-matching anchor row sampled
        sg = np.where(wa > 0, nb / np.maximum(wa, 1e-300), 1.0)
        glo, ghi = _one_side(g, sd, np.where(allmatch, 0, 1), wanc, sg)
        use_g = (m > 0) & (est > nb / 2)
        lo, hi = np.where(use_g, nb - ghi, lo), np.where(use_g, nb - glo, hi)
    lo = np.maximum(lo, m)
    hi = np.minimum(hi, nb)
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
    df = pd.read_csv(base + ".csv", header=None,
                     names=["ds", "method", "mem", "seed", "qid", "est", "ht", "sd", "cm", "m", "wq", "wa", "va", "wanc"])
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
                lo, hi = small_count_ci(est, ht, sd, g.m, g.wq, g.cm, g.wa, g.wanc)
                inside = (lo <= truth) & (truth <= hi)
                lo1, hi1 = small_count_ci(est, ht, sd, g.m, g.wq, g.cm)
                in1 = (lo1 <= truth) & (truth <= hi1)
                row["ci95_f"] = in1.mean()
                row["ci95_f_big"] = in1[big].mean()
                row["ci95_f_small"] = in1[truth < 100].mean() if (truth < 100).any() else np.nan
                row["ci_relwidth_f"] = np.median((hi1 - lo1) / truth)
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
