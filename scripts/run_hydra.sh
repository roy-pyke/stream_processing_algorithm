#!/usr/bin/env bash
# Hydra baseline: the authors' adaptation of Hydra shipped with the OmniSketch code
# (omniVLDB.HYDRA_VLDB.HYDRA_VLDB.ImpHydraStruct, L1 only), sized by their RamToPar, run
# unmodified through repro/omni_java/Harness.java. Every record updates 2^p - 1 subpopulations,
# so this is run on the 5-attribute synthetic streams only (and once on the 9-attribute
# journal-style stream with HYDRA_P9=1). The sketch seeds are fixed in that code: one repetition.
# Output is converted to the format of `bench acc` (method "hydra") under results/hydra/conv
# and aggregated into results/agg_hydra.csv.
# Usage: scripts/run_hydra.sh [parallel_jobs]
set -u
cd "$(dirname "$0")/.."
J=${1:-2}
JH=$(cat ext/java_home.txt)
CP=$(cat ext/build/cp.txt)
OUT=results/hydra
MEMS=${MEMS:-"1 5 10 20 50"}
mkdir -p $OUT/conv
run_h() {
  ds=$1; wl=$2; shift 2
  f=$OUT/$wl
  if [ ! -s "$f.log" ] || [ "$(grep -c '^HRUN' "$f.log")" -lt "$#" ]; then
    rm -f "$f.csv"
    "$JH/bin/java" -Xmx6g -cp "$CP" omniVLDB.Harness hydra data/$ds.bin work/$wl.wl "$f.csv" "$wl" 1 "$@" > "$f.log" 2> "$f.err"
  fi
}
export -f run_h
export OUT JH CP
{ for ds in zipf10_5M zipf13_5M zipf15_5M unif_5M; do echo "$ds $ds $MEMS"; done
  if [ -n "${HYDRA_P9:-}" ]; then echo "jz13p9_5M jz13p9-vldbj 10"; echo "jz13p9_5M jz13p9_5M 10"; fi; } |
  xargs -P "$J" -L 1 bash -c 'run_h "$@"' _
# convert to the bench format: name,spec,mem,seed,query,est,,,  +  RUN line
for f in $OUT/*.csv; do
  wl=$(basename "$f" .csv)
  for mem in $(cut -d, -f2 "$f" | sort -u); do
    awk -F, -v m="$mem" '$2==m {printf "%s,hydra,%s,1,%s,%s,,,\n", $1, $2, $4, $6}' "$f" > "$OUT/conv/$wl.$mem.1.csv"
    grep "^HRUN,$wl,$mem," "$OUT/$wl.log" | awk -F, '{printf "RUN,%s,hydra,%s,1,%s,%s,%s,%s,%s,%s,%s\n", $2, $3, $5, $6, $7, $8, $9, $10, $11}' > "$OUT/conv/$wl.$mem.1.log"
  done
done
.venv/bin/python scripts/aggregate.py $OUT/conv results/agg_hydra.csv work
