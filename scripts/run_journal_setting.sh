#!/usr/bin/env bash
# Insert-only comparison in the setting of the OmniSketch journal version's code
# (repository omni-deletes, Configurations/CashRegister/SynthCashRegister.json):
# 9 Zipf(1.3) attributes over {1..10000} from the authors' generator, queries with
# k = 1..9 predicates taken from sampled records, budgets 10-100 MB, and
#   omniv  OmniSketch with their configuration (w=28, d=3, 31 bits per stored hash)
#   ash    their sampling baseline (adaptive Sample-and-Hold on full tuples), their accounting
#   ashl   the same algorithm with lean accounting
#   resj   row reservoir with their accounting (64p+32 bits per row)
#   res    row reservoir, 32 bits per attribute value
# plus the PVLDB configuration (omnip), USS and CoStrata for reference.
# Both the journal-style workload and our own workload are run on this stream; the other
# streams are covered by scripts/run_accuracy.sh.
# Usage: scripts/run_journal_setting.sh [parallel_jobs]
set -u
cd "$(dirname "$0")/.."
J=${1:-3}
OUT=results/journal
mkdir -p $OUT work
JH=$(cat ext/java_home.txt)
CP=$(cat ext/build/cp.txt)
[ -s data/jz13p9_5M.bin ] || "$JH/bin/java" -Xmx6g -cp "$CP" omniVLDB.Harness gen data/jz13p9_5M.bin 5000000 9 1.3
[ -s work/jz13p9-vldbj.wl ] || bin/mkworkload data/jz13p9_5M.bin work/jz13p9-vldbj.wl vldbj 300 777 2>/dev/null
[ -s work/jz13p9_5M.wl ] || bin/mkworkload data/jz13p9_5M.bin work/jz13p9_5M.wl 200 6 4242 2>/dev/null
METHODS="res resj ash ashl omniv omnijs omnip omni:256:2 omni:1024:2 uss cs:0:1e9:1:0.25 cs:0:1:1:0.25"
NM=$(echo $METHODS | wc -w)
run_one() {
  ds=$1; wl=$2; mem=$3; seed=$4; shift 4
  f="$OUT/$wl.$mem.$seed"
  [ -s "$f.log" ] && [ "$(grep -c '^RUN' "$f.log")" -ge "$NM" ] && return 0
  rm -f "$f.csv"
  bin/bench acc data/$ds.bin work/$wl.wl "$f.csv" "$wl" "$mem" "$seed" "$@" > "$f.log" 2> "$f.err"
}
export -f run_one
export OUT NM
for pair in "jz13p9_5M jz13p9-vldbj" "jz13p9_5M jz13p9_5M"; do
  for mem in 10 25 50 100; do for seed in 1 2 3; do echo "$pair $mem $seed"; done; done
done | xargs -P "$J" -L 1 bash -c 'run_one $0 $1 $2 $3 '"$METHODS"
echo "journal-setting runs done: $(ls $OUT/*.log | wc -l)"
