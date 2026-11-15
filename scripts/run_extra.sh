#!/usr/bin/env bash
# Ablations and attribute scaling (census). Run after run_accuracy.sh; timing is separate
# (scripts/run_timing.sh, on an otherwise idle machine).
# Usage: scripts/run_extra.sh [parallel_jobs]
set -u
cd "$(dirname "$0")/.."
J=${1:-8}
mkdir -p results/abl results/census

# ---------------------------------------------------------------- ablations
#   house weight, allocation exponent, buckets per attribute (CoStrata with C1, as designed)
#   budget fraction of the pair tables (uniform sample) and of pair + triple tables
ABL="cs:64:1:1 cs:256:1:1 cs:1024:1:1 cs:4096:1:1 cs:16384:1:1 cs:0:0:1 cs:0:0.25:1 cs:0:1:1 cs:0:4:1 cs:0:1e9:1 cs:0:0:0.5 cs:0:0:0.75 cs:0:1e9:1:0.05 cs:0:1e9:1:0.1 cs:0:1e9:1:0.25 cs:0:1e9:1:0.5 cs:0:1e9:1:0.1:0.1 cs:0:1e9:1:0.15:0.15 cs:0:1e9:1:0.2:0.3 cs:0:1:1:0.25 cs:0:0.25:1:0.25 cs:0:1e9:1:0.25:0:1 cs:0:1e9:1:0.15:0.15:1 csh:0:0.25:0.25 csh:0:1:0.25 csh:0:4:0.25 csh uss ussc ussc:0:0.25 res"
NA=$(echo $ABL | wc -w)
run_abl() {
  ds=$1; mem=$2; seed=$3; shift 3
  f="results/abl/$ds.$mem.$seed"
  [ -s "$f.log" ] && [ "$(grep -c '^RUN' "$f.log")" -ge "$NA" ] && return 0
  rm -f "$f.csv"
  ./bin/bench acc data/$ds.bin work/$ds.wl "$f.csv" "$ds" "$mem" "$seed" "$@" > "$f.log" 2> "$f.err"
}
export -f run_abl
export NA
for ds in taxi mawi zipf15_5M; do for mem in 1 10; do for seed in 1 2 3 4 5; do echo "$ds $mem $seed"; done; done; done |
  xargs -P "$J" -L 1 bash -c 'run_abl $0 $1 $2 '"$ABL"

# ---------------------------------------------------------------- attribute scaling on census (p = 4..68)
CEN="res uss omnijs omnip omni:64:2 omni:256:2 omni:1024:2 cs:0:1e9 cs:0:1e9:1:0.25 cs cs:0:1:1:0.25 ussc:0:0.25 csh:0:1:0.25"
NC=$(echo $CEN | wc -w)
run_cen() {
  p=$1; seed=$2; shift 2
  f="results/census/census_p$p.10.$seed"
  [ -s "$f.log" ] && [ "$(grep -c '^RUN' "$f.log")" -ge "$NC" ] && return 0
  rm -f "$f.csv"
  ./bin/bench acc data/census.bin work/census_p$p.wl "$f.csv" "census_p$p" 10 "$seed" -p "$p" "$@" > "$f.log" 2> "$f.err"
}
export -f run_cen
export NC
for p in 4 8 16 32 68; do
  [ -s work/census_p$p.wl ] || ./bin/mkworkload data/census.bin work/census_p$p.wl 200 4 7 $p 2> work/census_p$p.wl.log
done
for p in 4 8 16 32 68; do for seed in 1 2 3; do echo "$p $seed"; done; done |
  xargs -P 4 -L 1 bash -c 'run_cen $0 $1 '"$CEN"
echo "extra experiments done"
