#!/usr/bin/env bash
# Bounded deletions in the setting of the journal version's code (#25; omni-deletes,
# Configurations/StrictTurnstile/SynthBounded*.json, and their data pipeline dataset-pipeline):
#   stream    9 attributes, Zipf(1.3) over {1..10000} from their generator (DeletesHarness gen); a residue
#             of 2^22 records (their 2^23-record stream, of which splitDeletesSynth.py keeps half) and
#             perc x 2^22 noise records, each inserted and deleted again at random later positions
#             (bin/mkdelstream; their spreadOutDeletes), perc in {0, 0.33, 1, 3, 9} (their "percs")
#   queries   on the residue: ours (200 per family and k, k <= 6) and the journal-style procedure
#             (k = 1..9 from 300 sampled residue records), in one file (work/delres.wl)
#   methods   s2   their OmniSketch with batched deletions (w = 28, d = 3, b = 31), their code
#             ash  their adaptive Sample-and-Hold, their code
#             cs:0:1e9[...]  our uniform designs RS, RS+C2, RS+C3 with deletions (bench del)
#   budgets   10, 25, 50 MB (their ramVals), seeds 1-3 (= their repetition)
# Results: results/deletes/<perc>.<mem>.<seed>.<method>.{csv,log} and results/deletes_streams.txt (updates, peak live);
# scripts/deletes.py aggregates them.
# Usage: scripts/run_deletes.sh [parallel_jobs]    (scripts/setup_omni_deletes.sh first)
set -u
cd "$(dirname "$0")/.."
J=${1:-3}
OUT=results/deletes
mkdir -p $OUT work
JH=$(cat ext/java_home.txt)
CP=$(cat ext/build/deletes_cp.txt)
NRES=4194304
PERCS=${PERCS:-"0 0.33 1 3 9"}
MEMS=${MEMS:-"10 25 50"}
SEEDS=${SEEDS:-"1 2 3"}
[ -s data/del_base.bin ] || "$JH/bin/java" -Xmx2g -cp "$CP" omni.DeletesHarness gen data/del_base.bin $((NRES * 10)) 9 1.3 || exit 1
for perc in $PERCS; do
  if [ ! -s work/del_$perc.ops ]; then
    res=""; [ -s data/del_res.bin ] || res=data/del_res.bin
    bin/mkdelstream data/del_base.bin $NRES $perc 25 work/del_$perc.ops $res > work/del_$perc.txt || exit 1
    grep -v "^$perc " results/deletes_streams.txt 2>/dev/null > work/del_streams.tmp; sed "s/^/$perc /" work/del_$perc.txt >> work/del_streams.tmp
    sort -g work/del_streams.tmp > results/deletes_streams.txt
  fi
done
if [ ! -s work/delres.wl ]; then
  bin/mkworkload data/del_res.bin work/delres-own.wl 200 6 4242 2>/dev/null || exit 1
  bin/mkworkload data/del_res.bin work/delres-vldbj.wl vldbj 300 777 2>/dev/null || exit 1
  cat work/delres-own.wl work/delres-vldbj.wl > work/delres.wl
fi
OURS="cs:0:1e9 cs:0:1e9:1:0.25 cs:0:1e9:1:0.15:0.15"
run_one() { # perc mem seed method
  local perc=$1 mem=$2 seed=$3 m=$4 f="$OUT/$1.$2.$3.$4"
  [ -s "$f.log" ] && grep -Eq "^J?DRUN" "$f.log" && return 0
  : > "$f.csv"
  if [ "$m" = ours ]; then
    bin/bench del data/del_base.bin work/del_$perc.ops work/delres.wl "$f.csv" del$perc $mem $seed $OURS > "$f.log" 2>&1
  else
    "$JH/bin/java" -Xmx5g -cp "$CP" omni.DeletesHarness run $m data/del_base.bin work/del_$perc.ops work/delres.wl \
      "$f.csv" del$perc $seed $mem 2>&1 | grep -v "^Found\|^OmniSketch Conf\|^K: \|^Running" > "$f.log"
  fi
}
export -f run_one
export OUT JH CP OURS
for perc in $PERCS; do for mem in $MEMS; do for seed in $SEEDS; do for m in ours s2 ash; do
  echo "$perc $mem $seed $m"
done; done; done; done | xargs -P "$J" -L 1 bash -c 'run_one "$@"' _
n=$(grep -El "^J?DRUN" $OUT/*.log 2>/dev/null | wc -l)
echo "deletion runs done: $n logs"
