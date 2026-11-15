#!/usr/bin/env bash
# Cross-validation against the authors' Java implementation.
# For every dataset, memory budget and repetition this runs, on the same records and queries,
#   (1) the unmodified upstream Java OmniSketch (through repro/omni_java/Harness.java),
#   (2) `bench ... omnij`, our bit-exact C++ port of it,
#   (3) `bench ... omnip omnil res`, our optimised OmniSketch and the row reservoir.
# scripts/validate_omni_java.py then compares (1) and (2) query by query and tabulates errors.
# Usage: scripts/validate_omni_java.sh [parallel_jobs]
set -u
cd "$(dirname "$0")/.."
J=${1:-3}
JH=$(cat ext/java_home.txt)
CP=$(cat ext/build/cp.txt)
OUT=results/validation
MEMS=${MEMS:-"1 10 50 100"}
DATASETS=${DATASETS:-"jzipf10_5M jzipf13_5M jzipf15_5M zipf10_5M zipf13_5M zipf15_5M unif_5M taxi mawi"}
REPS=${REPS:-"1 2 3"}
mkdir -p $OUT/java $OUT/cpp

# hash test vectors: JVM vs. C++ port
ok=1
for cfg in "1 25 8 3" "3 34 8 3" "7 31 8 3" "2 20 5 4"; do
  if ! diff -q <("$JH/bin/java" -cp "$CP" omniVLDB.Harness hashes $cfg) <(bin/jhash_test $cfg) >/dev/null; then ok=0; fi
done
echo "hash test vectors identical: $ok" | tee $OUT/hash_vectors.txt

run_one() {
  ds=$1; rep=$2
  jf=$OUT/java/$ds.$rep; cf=$OUT/cpp/$ds.$rep
  if [ ! -s "$jf.log" ] || [ "$(grep -c '^JRUN' "$jf.log")" -lt "$(echo $MEMS | wc -w)" ]; then
    rm -f "$jf.csv"
    "$JH/bin/java" -Xmx4g -cp "$CP" omniVLDB.Harness run data/$ds.bin work/$ds.wl "$jf.csv" "$ds" "$rep" $MEMS > "$jf.log" 2> "$jf.err"
  fi
  if [ ! -s "$cf.log" ] || [ "$(grep -c '^RUN' "$cf.log")" -lt "$((4 * $(echo $MEMS | wc -w)))" ]; then
    rm -f "$cf.csv"; : > "$cf.log"
    for mem in $MEMS; do
      bin/bench acc data/$ds.bin work/$ds.wl "$cf.csv" "$ds" "$mem" "$rep" omnij omnip omnil res >> "$cf.log" 2>> "$cf.err"
    done
  fi
}
export -f run_one
export OUT MEMS JH CP
for ds in $DATASETS; do for rep in $REPS; do echo "$ds $rep"; done; done | xargs -P "$J" -L 1 bash -c 'run_one $0 $1'
echo "validation runs done: $(ls $OUT/java/*.log | wc -l) java, $(ls $OUT/cpp/*.log | wc -l) cpp"
