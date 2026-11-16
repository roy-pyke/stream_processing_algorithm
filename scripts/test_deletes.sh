#!/usr/bin/env bash
# Checks of the deletion support of the uniform designs (CoStrata::remove, bench del), on a small
# stream from data/zipf13_5M.bin: residue 200,000 records, 3 x 200,000 noise records inserted and
# deleted again (peak live 2.04 x residue).
#   1. the counters (C1 and the pair tables) after the update stream equal those of an insert-only
#      run over the residue: the Count-Min bound "cm" of every query is identical
#   2. over 20 seeds, the mean estimate of every query with true answer >= 200 is within sampling
#      error of the truth (z-scores with mean ~ 0 and sd ~ 1)
set -eu
cd "$(dirname "$0")/.."
T=work/test_deletes
mkdir -p $T
[ -s $T/s.ops ] || bin/mkdelstream data/zipf13_5M.bin 200000 3 7 $T/s.ops $T/res.bin
[ -s $T/res.wl ] || bin/mkworkload $T/res.bin $T/res.wl 100 4 11 2> /dev/null
: > $T/del.csv
: > $T/acc.csv
: > $T/many.csv
bin/bench del data/zipf13_5M.bin $T/s.ops $T/res.wl $T/del.csv t 1 1 cs:0:1e9 cs:0:1e9:1:0.25 cs:0:1e9:1:0.15:0.15 > /dev/null
bin/bench acc $T/res.bin $T/res.wl $T/acc.csv t 1 1 cs:0:1e9 cs:0:1e9:1:0.25 cs:0:1e9:1:0.15:0.15 > /dev/null
for s in $(seq 1 20); do bin/bench del data/zipf13_5M.bin $T/s.ops $T/res.wl $T/many.csv t 1 $s cs:0:1e9:1:0.25 > /dev/null; done
.venv/bin/python - "$T" <<'PY'
import sys
import numpy as np, pandas as pd
T = sys.argv[1]
cols = "name,spec,mem,seed,q,est,ht,sd,cm,m,wmax,wa,va,wanc".split(",")
truth = np.array([int(l.split()[2]) for l in open(f"{T}/res.wl")])
d, a = pd.read_csv(f"{T}/del.csv", names=cols), pd.read_csv(f"{T}/acc.csv", names=cols)
ok = True
for sp in d.spec.unique():
    same = np.array_equal(d[d.spec == sp].sort_values("q").cm.values, a[a.spec == sp].sort_values("q").cm.values)
    print(f"{sp}: counters identical to the residue-only run: {same}")
    ok &= same
piv = pd.read_csv(f"{T}/many.csv", names=cols).pivot(index="q", columns="seed", values="est").values
se = piv.std(1, ddof=1) / np.sqrt(piv.shape[1])
big = (truth >= 200) & (se > 0)
z = (piv.mean(1) - truth)[big] / se[big]
print(f"unbiasedness, {big.sum()} queries with truth >= 200: mean z {z.mean():.3f}, sd z {z.std():.3f}")
ok &= abs(z.mean()) < 0.25 and 0.7 < z.std() < 1.4
print("PASS" if ok else "FAIL")
sys.exit(0 if ok else 1)
PY
