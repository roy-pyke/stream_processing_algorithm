#!/usr/bin/env bash
# Synthetic streams produced by the authors' generator (omniVLDB.ZipfGenerator), 5M records,
# 5 attributes, domain {1..10000}, plus query workloads from our generator.
# The upstream generator is not seedable; the md5 sums of the files used for the results in
# this repository are recorded in results/validation/md5_jzipf.txt.
set -euo pipefail
cd "$(dirname "$0")/.."
JH=$(cat ext/java_home.txt)
CP=$(cat ext/build/cp.txt)
mkdir -p work results/validation
for a in 1.0 1.3 1.5; do
  n=jzipf$(echo $a | tr -d .)_5M
  [ -s data/$n.bin ] || "$JH/bin/java" -Xmx4g -cp "$CP" omniVLDB.Harness gen data/$n.bin 5000000 5 $a
  bin/mkworkload data/$n.bin work/$n.wl 200 5 4242 2>/dev/null
done
md5sum data/jzipf*.bin > results/validation/md5_jzipf.txt
