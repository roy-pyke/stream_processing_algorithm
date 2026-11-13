#!/usr/bin/env bash
# Update throughput and query latency on one pinned core. Run on an otherwise idle machine.
set -u
cd "$(dirname "$0")/.."
mkdir -p results/time
# ---------------------------------------------------------------- timing (single pinned core, machine otherwise idle)
# Timings are only meaningful on a quiet, dedicated machine: set SKIP_TIMING=1 on shared VMs.
T=results/time/time.log
: > "$T"
for ds in taxi mawi zipf13_5M; do
  for mem in 1 10 50; do
    taskset -c 0 ./bin/bench time data/$ds.bin $ds $mem 1 3 res uss ussc omnip omnil omni:256:2 omni:1024:2 cs >> "$T"
  done
done
for p in 4 8 16 32 68; do
  taskset -c 0 ./bin/bench time data/census.bin census_p$p 10 1 3 -p $p res uss omnip omni:256:2 cs >> "$T"
done
# query latency (pinned, sequential)
Q=results/time/query.log
: > "$Q"
for ds in taxi mawi; do
  for mem in 1 10 50; do
    taskset -c 0 ./bin/bench acc data/$ds.bin work/$ds.wl /dev/null $ds $mem 1 res uss omnip omni:256:2 omni:1024:2 cs >> "$Q"
  done
done
echo "timing done"
