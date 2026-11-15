#!/usr/bin/env bash
# High-precision unbiasedness probes on small synthetic streams (many compactions, many seeds).
#   1. CoStrata, 20,000 records, 3 attributes, 500 stored rows
#   2. CoStrata-H, same stream (almost no repeated tuples)
#   3. CoStrata-H, 2 attributes over 50 values: at most 2,500 distinct tuples, 200 entries,
#      so that held counts and re-anchoring dominate
# Writes results/bias_probe.txt. Usage: scripts/run_bias_probe.sh [seeds]
set -u
cd "$(dirname "$0")/.."
S=${1:-400000}
OUT=results/bias_probe.txt
{ ./bin/bias_probe cs 20000 3 500 64 1 "$S" > "$OUT.1" & ./bin/bias_probe csh 20000 3 500 64 1 "$S" > "$OUT.2" &
  ./bin/bias_probe csh 20000 2 200 64 1 "$S" 50 > "$OUT.3" & wait; }
cat "$OUT.1" "$OUT.2" "$OUT.3" > "$OUT" && rm -f "$OUT.1" "$OUT.2" "$OUT.3"
cat "$OUT"
