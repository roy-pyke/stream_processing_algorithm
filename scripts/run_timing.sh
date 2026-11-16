#!/usr/bin/env bash
# Update throughput and query latency on one pinned core. Run on an otherwise idle machine.
set -u
cd "$(dirname "$0")/.."
mkdir -p results/time
# pin to one core where taskset exists (Linux); on macOS the scheduler decides
PIN=""
command -v taskset >/dev/null 2>&1 && PIN="taskset -c 0"
T=results/time/time.log
: > "$T"
TM="res uss coco ussc ussc:0:0.25 cococ:0:0.25 omnij omnip omnil omniv omni:256:2 omni:1024:2 ash cs:0:1e9 cs:0:1e9:1:0.25 cs:0:1e9:1:0.15:0.15 cs cs:0:1:1:0.25 csh:0:1:0.25 ussc:0:0.15:0.15 csh:0:1:0.15:0.15"
for ds in taxi mawi zipf13_5M; do
  for mem in 1 10 50; do
    $PIN ./bin/bench time data/$ds.bin $ds $mem 1 3 $TM >> "$T"
  done
done
for p in 4 8 16 32 68; do
  $PIN ./bin/bench time data/census.bin census_p$p 10 1 3 -p $p res uss omnip omni:256:2 cs:0:1e9 cs:0:1e9:1:0.25 cs cs:0:1:1:0.25 >> "$T"
done
# range predicates (dyadic anchors versus OmniSketch with dyadic ranges) and pair-table selection;
# both need a workload (-w) for the range attributes and for psel=wl
for ds in taxi mawi; do
  for mem in 1 10 50; do
    $PIN ./bin/bench time data/$ds.bin ${ds}_range $mem 1 3 -w work/${ds}_range.wl cs:0:1:1:0.25 cs/dy=1 cs:0:1:1:0.25/dy=1 \
      omnid:64:2 omnid:256:2 omnid:1024:2 omnid:256:3 omnid:1024:1 omnid:4096:1 >> "$T"
  done
done
$PIN ./bin/bench time data/census.bin census_p68 10 1 3 -p 68 -w work/census_p68_tmpl.wl cs cs:0:1:1:0.25 \
  cs:0:1:1:0.25/psel=wl/np=68 cs:0:1:1:0.25/psel=wl/np=272 cs:0:1:1:0.25/psel=sc/np=272 >> "$T"
# query latency (pinned, sequential)
Q=results/time/query.log
: > "$Q"
for ds in taxi mawi; do
  for mem in 1 10 50; do
    $PIN ./bin/bench acc data/$ds.bin work/$ds.wl /dev/null $ds $mem 1 res uss omnijs omnip omniv omni:256:2 omni:1024:2 cs:0:1e9:1:0.25 cs:0:1e9:1:0.25:0:0:1 cs cs:0:1:1:0.25 csh:0:1:0.25 coco cococ:0:0.25 >> "$Q"
  done
done
echo "timing done"
