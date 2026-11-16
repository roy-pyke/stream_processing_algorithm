#!/usr/bin/env bash
# Experiments added before the SIGMOD submission (issues #18-#23). Run after run_accuracy.sh;
# each part writes its own result directory and can be run alone: run_presub.sh <part> [jobs].
#   main   (#18, #21, #23) CocoSketch (plain and calibrated), CoStrata with the query index charged
#          to the budget, and the CoStrata family again with interval diagnostics (the estimates
#          are identical to results/acc; the extra columns feed the small-count intervals)
#   range  (#19) range predicates on Taxi and MAWI, against OmniSketch with dyadic ranges
#   merge  (#22) synopses built on 2..16 parts of the stream and merged
#   pairs  (#20) pair-table selection on Census, p = 16/32/68, random and template workloads
#   dmv    the configurations of "main" that are new (CocoSketch, index charged) on DMV: 3 budgets, 3 seeds
#   ci     (#18) the methods of the coverage table again, with the anchor columns that the
#          two-sided intervals need (estimates identical to results/acc)
#   timing update throughput of the new configurations (sequential, one job; run on an idle machine)
set -u
cd "$(dirname "$0")/.."
PART=${1:-all}
J=${2:-8}

run_jobs() { # dir; job lines on stdin: <name> <data> <wl> <mem> <seed> [-p k] <spec>...
  local dir=$1
  mkdir -p "$dir"
  export dir
  xargs -P "$J" -L 1 bash -c '
    name=$0; data=$1; wl=$2; mem=$3; seed=$4; shift 4
    nm=$#; [ "${1:-}" = -p ] && nm=$((nm - 2))
    f="$dir/$name.$mem.$seed"
    [ -s "$f.log" ] && [ "$(grep -c "^RUN" "$f.log")" -ge "$nm" ] && exit 0
    : > "$f.csv"
    ./bin/bench acc "data/$data.bin" "$wl" "$f.csv" "$name" "$mem" "$seed" "$@" > "$f.log" 2> "$f.err"'
}

DATASETS="zipf10_5M zipf13_5M zipf15_5M unif_5M taxi mawi"
MEMS="0.5 1 2 5 10 20 50"

if [ "$PART" = main ] || [ "$PART" = all ]; then
  M="coco cococ cococ:0:0.25 cs cs:0:1:1:0.25 cs:0:1e9 cs:0:1e9:1:0.25 csh:0:1:0.25 cs/ix=2 cs:0:1:1:0.25/ix=2 cs:0:1e9:1:0.25/ix=2"
  NM=$(echo $M | wc -w)
  for mem in $MEMS; do for ds in $DATASETS; do for seed in 1 2 3 4 5; do
    echo "$ds $ds work/$ds.wl $mem $seed $M"; done; done; done | run_jobs results/acc2
fi

if [ "$PART" = range ] || [ "$PART" = all ]; then
  [ -s work/taxi_range.wl ] || ./bin/mkworkload data/taxi.bin work/taxi_range.wl range 100 41 6,8,9,11
  [ -s work/mawi_range.wl ] || ./bin/mkworkload data/mawi.bin work/mawi_range.wl range 100 42 2,3,5
  M="res uss coco cococ:0:0.25 ussc:0:0.25 cs cs:0:1:1:0.25 cs/dy=1 cs:0:1:1:0.25/dy=1 cs:0:1e9:1:0.25 cs:0:1e9:1:0.25/dy=1 csh:0:1:0.25 omnid:64:2 omnid:256:2 omnid:1024:2 omnid:256:3 omnid:1024:1 omnid:4096:1"
  NM=$(echo $M | wc -w)
  for mem in $MEMS; do for ds in taxi mawi; do for seed in 1 2 3; do
    echo "${ds}_range $ds work/${ds}_range.wl $mem $seed $M"; done; done; done | run_jobs results/range
fi

if [ "$PART" = merge ] || [ "$PART" = all ]; then
  M="cs cs:0:1:1:0.25 cs:0:1e9:1:0.25 omni:256:2 omni:1024:2"
  NM=$(echo $M | wc -w)
  mkdir -p results/merge
  export NM
  for ds in taxi mawi zipf13_5M; do for mem in 1 10; do for k in 1 2 4 8 16; do for seed in 1 2 3 4 5; do
    echo "$ds $mem $k $seed"; done; done; done; done |
    xargs -P "$J" -L 1 bash -c '
      ds=$0; mem=$1; k=$2; seed=$3
      f="results/merge/${ds}_m$k.$mem.$seed"
      [ -s "$f.log" ] && [ "$(grep -c "^RUN" "$f.log")" -ge "$NM" ] && exit 0
      : > "$f.csv"
      ./bin/bench merge "data/$ds.bin" "work/$ds.wl" "$f.csv" "$ds" "$mem" "$seed" "$k" '"$M"' > "$f.log" 2> "$f.err"'
fi

if [ "$PART" = pairs ] || [ "$PART" = all ]; then
  for p in 16 32 68; do
    [ -s work/census_p$p.wl ] || ./bin/mkworkload data/census.bin work/census_p$p.wl 200 4 7 $p 2> work/census_p$p.wl.log
    [ -s work/census_p${p}_tmpl.wl ] || MKW_TEMPLATES=20 ./bin/mkworkload data/census.bin work/census_p${p}_tmpl.wl 200 4 7 $p 2> work/census_p${p}_tmpl.wl.log
  done
  for p in 16 32 68; do
    case $p in 16) NP="16 32 64" ;; 32) NP="32 64 128 256" ;; 68) NP="68 136 272 544" ;; esac
    M="cs cs:0:1:1:0.25"
    for sel in wl mi sc rnd; do for np in $NP; do M="$M cs:0:1:1:0.25/psel=$sel/np=$np"; done; done
    NM=$(echo $M | wc -w)
    for wl in census_p$p census_p${p}_tmpl; do for seed in 1 2 3; do
      echo "$wl census work/$wl.wl 10 $seed -p $p $M"; done; done
  done | run_jobs results/pairs
fi
if [ "$PART" = dmv ] || [ "$PART" = all ]; then
  M="coco cococ cococ:0:0.25 cs/ix=2 cs:0:1:1:0.25/ix=2 cs:0:1e9:1:0.25/ix=2"
  for mem in 1 10 50; do for seed in 1 2 3; do
    echo "dmv dmv work/dmv.wl $mem $seed $M"; done; done | run_jobs results/acc2
fi

if [ "$PART" = ci ] || [ "$PART" = all ]; then
  M="cs cs:0:1:1:0.25 cs:0:1e9:1:0.25 csh:0:1:0.25"
  for mem in $MEMS; do for ds in $DATASETS; do for seed in 1 2 3 4 5; do
    echo "$ds $ds work/$ds.wl $mem $seed $M"; done; done; done | run_jobs results/ci
fi

if [ "$PART" = timing ]; then   # not part of "all": needs an otherwise idle machine
  T=results/time/presub_time.log
  mkdir -p results/time
  : > "$T"
  for ds in taxi mawi; do for mem in 1 10 50; do
    ./bin/bench time data/$ds.bin $ds $mem 1 3 res uss coco ussc:0:0.25 cococ:0:0.25 cs cs:0:1:1:0.25 >> "$T"
    ./bin/bench time data/$ds.bin ${ds}_range $mem 1 3 -w work/${ds}_range.wl cs:0:1:1:0.25 cs/dy=1 cs:0:1:1:0.25/dy=1 \
      omnid:64:2 omnid:256:2 omnid:1024:2 omnid:256:3 omnid:1024:1 omnid:4096:1 >> "$T"
  done; done
  for p in 32 68; do
    case $p in 32) NP="32 128" ;; 68) NP="68 272" ;; esac
    S="cs cs:0:1:1:0.25"
    for np in $NP; do S="$S cs:0:1:1:0.25/psel=wl/np=$np cs:0:1:1:0.25/psel=sc/np=$np"; done
    ./bin/bench time data/census.bin census_p$p 10 1 3 -p $p -w work/census_p${p}_tmpl.wl $S >> "$T"
  done
fi
echo "presubmission experiments ($PART) done"
