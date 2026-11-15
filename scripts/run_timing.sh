#!/usr/bin/env bash
# Update throughput and query latency on one pinned core. Run on an otherwise idle machine.
set -u
cd "$(dirname "$0")/.."
mkdir -p results/time
T=results/time/time.log
: > "$T"
TM="res uss ussc ussc:0:0.25 omnij omnip omnil omniv omni:256:2 omni:1024:2 ash cs:0:1e9 cs:0:1e9:1:0.25 cs:0:1e9:1:0.15:0.15 cs cs:0:1:1:0.25 csh:0:1:0.25"
for ds in taxi mawi zipf13_5M; do
  for mem in 1 10 50; do
    taskset -c 0 ./bin/bench time data/$ds.bin $ds $mem 1 3 $TM >> "$T"
  done
done
for p in 4 8 16 32 68; do
  taskset -c 0 ./bin/bench time data/census.bin census_p$p 10 1 3 -p $p res uss omnip omni:256:2 cs:0:1e9 cs:0:1e9:1:0.25 cs cs:0:1:1:0.25 >> "$T"
done
# query latency (pinned, sequential)
Q=results/time/query.log
: > "$Q"
for ds in taxi mawi; do
  for mem in 1 10 50; do
    taskset -c 0 ./bin/bench acc data/$ds.bin work/$ds.wl /dev/null $ds $mem 1 res uss omnijs omnip omniv omni:256:2 omni:1024:2 cs:0:1e9:1:0.25 cs cs:0:1:1:0.25 csh:0:1:0.25 >> "$Q"
  done
done
echo "timing done"
