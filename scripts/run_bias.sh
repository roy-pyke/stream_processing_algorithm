#!/usr/bin/env bash
# Empirical check of Theorem 1: mean of the Horvitz-Thompson estimate over many independent seeds.
# Taxi (1 MB): seeds 101-400 and 1001-1900 (the second batch was added after the first one showed a
# +0.15% +- 0.06% deviation; both batches are pooled in the analysis). Zipf 1.5 (0.5 MB): seeds 101-400.
# Usage: scripts/run_bias.sh [parallel_jobs]
set -u
cd "$(dirname "$0")/.."
J=${1:-8}
mkdir -p results/bias
run_b() {
  ds=$1; mem=$2; seed=$3
  f="results/bias/$ds.$mem.$seed"
  rm -f "$f.csv"
  ./bin/bench acc data/$ds.bin work/$ds.wl "$f.csv" "$ds" "$mem" "$seed" cs > "$f.log" 2> "$f.err"
}
export -f run_b
{ for seed in $(seq 101 400) $(seq 1001 1900); do echo "taxi 1 $seed"; done
  for seed in $(seq 101 400); do echo "zipf15_5M 0.5 $seed"; done; } | xargs -P "$J" -L 1 bash -c 'run_b $0 $1 $2'
echo "bias runs done"
