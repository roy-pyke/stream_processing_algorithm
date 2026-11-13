#!/usr/bin/env bash
# Main accuracy matrix: every method at equal memory, 5 seeds.
# Usage: scripts/run_accuracy.sh [parallel_jobs]
set -u
cd "$(dirname "$0")/.."
J=${1:-8}
SKIP_LIST=${SKIP_LIST:-}
OUT=results/acc
mkdir -p "$OUT"
METHODS="res uss ussc omnip omnil omni:64:2 omni:256:2 omni:1024:2 omni:256:3 omni:1024:1 omni:4096:1 omni:4096:2 cs cs:0:1e9 cs:0:0"
jobs=()
for ds in zipf10_5M zipf13_5M zipf15_5M unif_5M taxi mawi; do
  for mem in 0.5 1 2 5 10 20 50; do
    for seed in 1 2 3 4 5; do
      jobs+=("$ds $mem $seed")
    done
  done
done
run_one() {
  ds=$1; mem=$2; seed=$3; shift 3
  f="$OUT/$ds.$mem.$seed"
  # jobs already finished on another machine (results are seed-deterministic) can be listed in $SKIP_LIST
  [ -n "${SKIP_LIST:-}" ] && grep -qx "$ds.$mem.$seed" "$SKIP_LIST" && return 0
  [ -s "$f.log" ] && [ "$(wc -l < "$f.log")" -ge 15 ] && return 0
  rm -f "$f.csv"
  ./bin/bench acc data/$ds.bin work/$ds.wl "$f.csv" "$ds" "$mem" "$seed" "$@" > "$f.log" 2> "$f.err"
}
export -f run_one
export OUT SKIP_LIST
printf '%s\n' "${jobs[@]}" | xargs -P "$J" -L 1 bash -c 'run_one $0 $1 $2 '"$METHODS"
echo "accuracy matrix done: $(ls $OUT/*.log | wc -l) jobs"
