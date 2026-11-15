#!/usr/bin/env bash
# Empirical unbiasedness check: mean of the Horvitz-Thompson-type estimate over many independent
# seeds, for CoStrata (Theorem 1 applies) and for the hold variant CoStrata-H (no proof given).
#   taxi 1 MB, seeds 101-700;  zipf15_5M 0.5 MB, seeds 101-400;  mawi 1 MB, seeds 101-400
# Usage: scripts/run_bias.sh [parallel_jobs]
set -u
cd "$(dirname "$0")/.."
J=${1:-8}
mkdir -p results/bias
run_b() {
  ds=$1; mem=$2; seed=$3
  f="results/bias/$ds.$mem.$seed"
  [ -s "$f.log" ] && [ "$(grep -c '^RUN' "$f.log")" -ge 2 ] && return 0
  rm -f "$f.csv"
  ./bin/bench acc data/$ds.bin work/$ds.wl "$f.csv" "$ds" "$mem" "$seed" cs csh > "$f.log" 2> "$f.err"
}
export -f run_b
{ for seed in $(seq 101 700); do echo "taxi 1 $seed"; done
  for seed in $(seq 101 400); do echo "zipf15_5M 0.5 $seed"; done
  for seed in $(seq 101 400); do echo "mawi 1 $seed"; done; } | xargs -P "$J" -L 1 bash -c 'run_b $0 $1 $2'
echo "bias runs done"
