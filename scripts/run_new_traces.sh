#!/usr/bin/env bash
# The two real streams of the OmniSketch paper that we do not have yet (#12): SNMP (#11) and CAIDA (#10).
# Put the raw data here, then run this script; every step is skipped if its output exists.
#   data/raw/snmp/    fall03.tar.gz extracted (day directories 031101 ... with *.snmp.gz)
#   data/raw/caida/   CAIDA anonymised trace files (*.pcap or *.pcap.gz), converted in name order;
#                     the OmniSketch paper uses 109 M packets (about 4-5 one-minute files of a 2019 trace)
# Steps: convert (scripts/convert.py snmp, bin/pcap2bin --caida), workloads (bin/mkworkload, 200 queries
# per family and k, k <= 6), the accuracy matrix of run_accuracy.sh plus the triple-table variants,
# CocoSketch and the index-charged runs, and the aggregation into results/agg_new.csv, which
# scripts/plots.py merges into the figures and tables.
# CAIDA holds 109 M x 6 x 4 bytes = 2.6 GB per process: CAIDA_JOBS (default 2) limits its parallelism.
# Usage: scripts/run_new_traces.sh [parallel_jobs]
set -u
cd "$(dirname "$0")/.."
J=${1:-6}
CJ=${CAIDA_JOBS:-2}
OUT=results/new
mkdir -p "$OUT" work
DS=""
if [ -d data/raw/snmp ]; then
  [ -s data/snmp.bin ] || .venv/bin/python scripts/convert.py snmp data/raw/snmp data/snmp.bin || exit 1
  DS="$DS snmp"
fi
if [ -d data/raw/caida ]; then
  if [ ! -s data/caida.bin ]; then
    files=$(ls data/raw/caida/*.pcap data/raw/caida/*.pcap.gz 2>/dev/null | sort)
    [ -n "$files" ] || { echo "no pcap files in data/raw/caida"; exit 1; }
    bin/pcap2bin --caida data/caida.bin $files || exit 1
  fi
  DS="$DS caida"
fi
[ -n "$DS" ] || { echo "neither data/raw/snmp nor data/raw/caida exists (see docs/datasets.md)"; exit 1; }
for ds in $DS; do
  [ -s work/$ds.wl ] || bin/mkworkload data/$ds.bin work/$ds.wl 200 6 4242 || exit 1
done
M="res resj uss ash ashl cs:0:1e9 cs:0:1e9:1:0.25 cs:0:1e9:1:0.15:0.15 cs cs:0:1:1:0.25 cs:0:0 ussc ussc:0:0.25 csh:0:1:0.25"
M="$M omnijs omnip omniv omnil omni:28:3 omni:64:2 omni:256:2 omni:1024:2 omni:256:3 omni:1024:1 omni:4096:1 omni:4096:2"
M="$M ussc:0:0.15:0.15 csh:0:1:0.15:0.15 coco cococ cococ:0:0.25 cs/ix=2 cs:0:1:1:0.25/ix=2 cs:0:1e9:1:0.25/ix=2"
for ds in $DS; do
  jobs=$J; [ "$ds" = caida ] && jobs=$CJ
  OUT=$OUT METHODS="$M" DATASETS="$ds" scripts/run_accuracy.sh "$jobs"
done
.venv/bin/python scripts/aggregate.py $OUT results/agg_new.csv work
echo "new traces done:$DS"
