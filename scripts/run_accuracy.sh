#!/usr/bin/env bash
# Main accuracy matrix: every method at equal memory, 5 seeds.
# Usage: scripts/run_accuracy.sh [parallel_jobs]
#
# Method specs (see src/bench.cpp):
#   OmniSketch   omnijs               bit-exact port of the authors' Java code (their B, hashes, accounting)
#                omnip                our implementation, the authors' (w,d,B), 63-bit hashes
#                omniv                our implementation, configuration of the journal version's code
#                omnil omni:w:d       configuration search under lean, oracle accounting
#   journal      ash ashl resj        the journal version's sampling baseline and accountings
#   row samples  res uss              reservoir, Unbiased Space-Saving
#                cs:0:1e9[...]        calibrated uniform sample: C1 / C2 (pairs) / C3 (pairs + triples)
#                cs  cs:0:1:1:0.25    CoStrata with C1 / C2;  cs:0:0  without house stratum
#                ussc ussc:0:0.25     USS with C1 / C2
#                csh:0:1:0.25         CoStrata-H (hold variant) with C2
set -u
cd "$(dirname "$0")/.."
J=${1:-8}
OUT=${OUT:-results/acc}
mkdir -p "$OUT"
METHODS=${METHODS:-"res resj uss ash ashl cs:0:1e9 cs:0:1e9:1:0.25 cs:0:1e9:1:0.15:0.15 cs cs:0:1:1:0.25 cs:0:0 ussc ussc:0:0.25 csh:0:1:0.25 omnijs omnip omniv omnil omni:28:3 omni:64:2 omni:256:2 omni:1024:2 omni:256:3 omni:1024:1 omni:4096:1 omni:4096:2"}
DATASETS=${DATASETS:-"zipf10_5M zipf13_5M zipf15_5M unif_5M taxi mawi"}
MEMS=${MEMS:-"0.5 1 2 5 10 20 50"}
SEEDS=${SEEDS:-"1 2 3 4 5"}
NM=$(echo $METHODS | wc -w)
run_one() {
  ds=$1; mem=$2; seed=$3; shift 3
  f="$OUT/$ds.$mem.$seed"
  [ -s "$f.log" ] && [ "$(grep -c '^RUN' "$f.log")" -ge "$NM" ] && return 0
  rm -f "$f.csv"
  ./bin/bench acc data/$ds.bin work/$ds.wl "$f.csv" "$ds" "$mem" "$seed" "$@" > "$f.log" 2> "$f.err"
}
export -f run_one
export OUT NM
# small budgets first, so that partial results are useful early
for mem in $MEMS; do for ds in $DATASETS; do for seed in $SEEDS; do echo "$ds $mem $seed"; done; done; done |
  xargs -P "$J" -L 1 bash -c 'run_one $0 $1 $2 '"$METHODS"
echo "accuracy matrix done: $(ls $OUT/*.log | wc -l) jobs"
