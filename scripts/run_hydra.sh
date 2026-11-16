#!/usr/bin/env bash
# Hydra baseline: the authors' adaptation of Hydra shipped with the OmniSketch code
# (omniVLDB.HYDRA_VLDB.HYDRA_VLDB.ImpHydraStruct, L1 only), sized by their RamToPar, run
# unmodified through repro/omni_java/Harness.java. Every record updates 2^p - 1 subpopulations,
# so this is run on the 5-attribute synthetic streams (and once on the 9-attribute journal-style
# stream with HYDRA_P9=1); HYDRA_REAL=1 adds MAWI (10 attributes, 1,023 updates per record) and
# Taxi (12 attributes, 4,095), one JVM per budget (on an M4 Pro, 4 at a time: about 1 h per budget
# on MAWI and 2-2.6 h on Taxi; HYDRA_XMX=2g suffices). The sketch seeds are fixed in that
# code: one repetition. HYDRA_XMX sets the JVM heap (default 6g).
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
run_h() { # dataset workload tag mem...: output files $OUT/<tag>.{csv,log,err}
  ds=$1; wl=$2; tag=$3; shift 3
  f=$OUT/$tag
  if [ ! -s "$f.log" ] || [ "$(grep -c '^HRUN' "$f.log")" -lt "$#" ]; then
    : > "$f.csv"
    "$JH/bin/java" -Xmx${HYDRA_XMX:-6g} -cp "$CP" omniVLDB.Harness hydra data/$ds.bin work/$wl.wl "$f.csv" "$wl" 1 "$@" > "$f.log" 2> "$f.err"
  fi
}
export -f run_h
export OUT JH CP
{ for ds in zipf10_5M zipf13_5M zipf15_5M unif_5M; do echo "$ds $ds $ds $MEMS"; done
  if [ -n "${HYDRA_P9:-}" ]; then echo "jz13p9_5M jz13p9-vldbj jz13p9-vldbj 10"; echo "jz13p9_5M jz13p9_5M jz13p9_5M 10"; fi
  # real streams: one JVM per budget (tag <stream>@<mem>)
  if [ -n "${HYDRA_REAL:-}" ]; then for ds in taxi mawi; do for mem in $MEMS; do echo "$ds $ds $ds@$mem $mem"; done; done; fi; } |
  xargs -P "$J" -L 1 bash -c 'run_h "$@"' _
# convert to the bench format: name,spec,mem,seed,query,est,,,  +  RUN line
for f in $OUT/*.csv; do
  [ -s "$f" ] || continue
  log=${f%.csv}.log
  name=$(head -1 "$f" | cut -d, -f1)
  for mem in $(cut -d, -f2 "$f" | sort -u); do
    awk -F, -v m="$mem" '$2==m {printf "%s,hydra,%s,1,%s,%s,,,\n", $1, $2, $4, $6}' "$f" > "$OUT/conv/$name.$mem.1.csv"
    grep "^HRUN,$name,$mem," "$log" | awk -F, '{printf "RUN,%s,hydra,%s,1,%s,%s,%s,%s,%s,%s,%s\n", $2, $3, $5, $6, $7, $8, $9, $10, $11}' > "$OUT/conv/$name.$mem.1.log"
  done
done
.venv/bin/python scripts/aggregate.py $OUT/conv results/agg_hydra.csv work
