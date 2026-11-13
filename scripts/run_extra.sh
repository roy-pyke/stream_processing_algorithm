#!/usr/bin/env bash
# Ablations, attribute scaling (census), and timing. Run after run_accuracy.sh.
# Usage: scripts/run_extra.sh [parallel_jobs]
set -u
cd "$(dirname "$0")/.."
J=${1:-8}
mkdir -p results/abl results/census results/time

# ---------------------------------------------------------------- ablations (CoStrata design space)
ABL="cs:64:1:1 cs:256:1:1 cs:1024:1:1 cs:4096:1:1 cs:16384:1:1 cs:0:0:1 cs:0:0.25:1 cs:0:1:1 cs:0:4:1 cs:0:1e9:1 cs:0:0:0.5 cs:0:0:0.75 res"
run_abl() {
  ds=$1; mem=$2; seed=$3; shift 3
  f="results/abl/$ds.$mem.$seed"
  rm -f "$f.csv"
  ./bin/bench acc data/$ds.bin work/$ds.wl "$f.csv" "$ds" "$mem" "$seed" "$@" > "$f.log" 2> "$f.err"
}
export -f run_abl
for ds in taxi mawi zipf15_5M; do for mem in 1 10; do for seed in 1 2 3 4 5; do echo "$ds $mem $seed"; done; done; done |
  xargs -P "$J" -L 1 bash -c 'run_abl $0 $1 $2 '"$ABL"

# ---------------------------------------------------------------- attribute scaling on census (p = 4..68)
CEN="res uss omnip omni:64:2 omni:256:2 omni:1024:2 cs cs:0:1e9"
run_cen() {
  p=$1; seed=$2; shift 2
  f="results/census/census_p$p.10.$seed"
  rm -f "$f.csv"
  ./bin/bench acc data/census.bin work/census_p$p.wl "$f.csv" "census_p$p" 10 "$seed" -p "$p" "$@" > "$f.log" 2> "$f.err"
}
export -f run_cen
for p in 4 8 16 32 68; do
  [ -s work/census_p$p.wl ] || ./bin/mkworkload data/census.bin work/census_p$p.wl 200 4 7 $p 2> work/census_p$p.wl.log
done
for p in 4 8 16 32 68; do for seed in 1 2 3; do echo "$p $seed"; done; done |
  xargs -P 4 -L 1 bash -c 'run_cen $0 $1 '"$CEN"

if [ -n "${SKIP_TIMING:-}" ]; then echo "extra experiments done (timing skipped)"; exit 0; fi
bash scripts/run_timing.sh
echo "extra experiments done"
