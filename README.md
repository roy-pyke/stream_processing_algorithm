# Calibrated row samples vs. sketch-in-sketch for ad-hoc predicate counting on streams

Code, aggregated results and paper draft for a study of COUNT queries with conjunctive predicates
chosen at query time on multi-attribute streams:

- an anatomy of **OmniSketch** (PVLDB 17(3); VLDB Journal 35(1)) and a reproduction with the
  **authors' Java implementation**, which is run unmodified and ported bit-exactly to C++;
- **calibrated row samples**: a row sample (uniform, USS, or the coordinated hash-stratified sampler
  **CoStrata**, optionally in sample-and-hold mode) whose estimates are calibrated to exact counters of
  hashed one-, two- and three-dimensional marginals;
- an experimental comparison at equal memory, including Hydra and the journal version's sampling
  baseline (adaptive Sample-and-Hold).

## Layout

| Path | Contents |
|---|---|
| `src/synopses.hpp` | reservoir, OmniSketch (ours), USS, aSH (journal baseline), marginal tables (`MargCal`), CoStrata, CoStrata-H |
| `src/omni_java_compat.hpp` | bit-exact C++ port of the authors' Java OmniSketch (hashes, parameter search, both estimators, accounting) |
| `src/bench.cpp` | benchmark driver (`acc`, `time`, `merge`); method specs are documented at the top of the file |
| `src/mkworkload.cpp`, `gen_synth.cpp`, `pcap2bin.cpp` | workload generator (incl. the query procedures of both OmniSketch papers), synthetic streams, pcap converter |
| `src/bias_probe.cpp`, `jhash_test.cpp` | unbiasedness probe; test vectors for the ported Java hash functions |
| `repro/omni_java/Harness.java` | driver for the authors' unmodified Java code (OmniSketch and their Hydra adaptation) |
| `scripts/` | setup of the upstream code, dataset conversion, experiment drivers, aggregation, plotting, paper checks |
| `results/` | aggregated results (`agg_*.csv`, `validation/`, `bias.json`, timing logs); raw per-run logs are not tracked |
| `paper/` | LaTeX draft (`paper.tex`, `results.tex`), generated `figs/`, `tables/`, `numbers.tex` |
| `docs/` | [`datasets.md`](docs/datasets.md): every stream (type, source, terms, checksums); [`references.md`](docs/references.md): where to get the papers and code we worked from |

## Build

```bash
make                                   # bin/bench bin/gen_synth bin/mkworkload bin/pcap2bin bin/jhash_test
g++ -O3 -march=native -std=c++17 -o bin/bias_probe src/bias_probe.cpp
python -m venv .venv && .venv/bin/pip install -r requirements.txt
scripts/setup_upstream.sh              # clones github.com/wiegerpunter/omnisketch @ 8dc67b5 into ext/, builds it
                                       # (uses the jars shipped in that repository; downloads Temurin JDK 21 if no javac)
```

## Data (not tracked)

See [`docs/datasets.md`](docs/datasets.md) for a description of every stream, its terms of use and the
checksums of the files behind the published results. Synthetic streams (including those from the authors'
generator, which cannot be regenerated bit for bit), query workloads and raw per-run logs are attached to the
release [`artifacts-2026-11-15`](https://github.com/roy-pyke/stream_processing_algorithm/releases/tag/artifacts-2026-11-15).

Real streams (place downloads in `data/raw/`):

- NYC Yellow Taxi 2019-01: https://d37ci6vzurychx.cloudfront.net/trip-data/yellow_tripdata_2019-01.parquet
- US Census 1990 (UCI): https://archive.ics.uci.edu/static/public/116/us+census+data+1990.zip
- MAWI trace: http://mawi.wide.ad.jp/mawi/samplepoint-F/2010/201001011400.dump.gz
- NY DMV registrations: https://data.ny.gov/api/views/w4pv-hbkt/rows.csv?accessType=DOWNLOAD

```bash
python scripts/convert.py taxi   data/raw/yellow_2019-01.parquet   data/taxi.bin
python scripts/convert.py census data/raw/USCensus1990.data.txt    data/census.bin
python scripts/convert.py dmv    data/raw/dmv/dmv.csv              data/dmv.bin
bin/pcap2bin <decompressed MAWI dump> data/mawi.bin

bin/gen_synth data/unif_5M.bin   5000000 5 0   10000 4
bin/gen_synth data/zipf10_5M.bin 5000000 5 1.0 10000 1
bin/gen_synth data/zipf13_5M.bin 5000000 5 1.3 10000 2
bin/gen_synth data/zipf15_5M.bin 5000000 5 1.5 10000 3
scripts/gen_upstream_synth.sh          # Zipf streams from the authors' generator (data/jzipf*_5M.bin)

for ds in zipf10_5M zipf13_5M zipf15_5M unif_5M; do bin/mkworkload data/$ds.bin work/$ds.wl 200 5 4242; done
for ds in taxi mawi dmv; do bin/mkworkload data/$ds.bin work/$ds.wl 200 6 4242; done
# range workloads (kind 4: one or two ranges plus equalities, k = 1..4; attributes by index)
bin/mkworkload data/taxi.bin work/taxi_range.wl range 100 41 6,8,9,11    # hour, distance, fare, duration
bin/mkworkload data/mawi.bin work/mawi_range.wl range 100 42 2,3,5       # ports, packet length
# Census template workloads (20 attribute sets; used for pair-table selection)
for p in 16 32 68; do MKW_TEMPLATES=20 bin/mkworkload data/census.bin work/census_p${p}_tmpl.wl 200 4 7 $p; done
```
(`scripts/run_presub.sh` creates the range and template workloads itself if they are missing.)

The SNMP (CRAWDAD Dartmouth campus, IEEE DataPort, free registration) and CAIDA (application required)
traces of the original evaluation are not included.

## Reproduce

```bash
scripts/validate_omni_java.sh 4 && python scripts/validate_omni_java.py   # Java code vs. C++ port, query by query
scripts/run_accuracy.sh 8        # accuracy matrix (5 seeds)
DATASETS=dmv MEMS="1 10 50" SEEDS="1 2 3" scripts/run_accuracy.sh 2          # DMV (large): 3 budgets, 3 seeds
OUT=results/acc3 METHODS="ussc:0:0.15:0.15 csh:0:1:0.15:0.15" scripts/run_accuracy.sh 8   # triple-table variants
OUT=results/acc3 METHODS="ussc:0:0.15:0.15 csh:0:1:0.15:0.15" DATASETS=dmv MEMS="1 10 50" SEEDS="1 2 3" scripts/run_accuracy.sh 4
scripts/run_journal_setting.sh 3 # journal-version setting (9 attributes, their baselines and accounting)
scripts/run_hydra.sh 2           # the authors' Hydra adaptation on the 5-attribute streams
scripts/run_extra.sh 4           # ablations, attribute scaling
scripts/run_bias.sh 8 && python scripts/bias_check.py ; scripts/run_bias_probe.sh
scripts/run_timing.sh            # update throughput / query latency (idle machine)
scripts/run_presub.sh all 8      # CocoSketch, index charged, interval diagnostics, ranges, merging, pair selection, DMV
scripts/run_presub.sh timing     # throughput of those configurations (sequential; idle machine)

python scripts/aggregate.py results/acc     results/agg_acc.csv     work
python scripts/aggregate.py results/acc3    results/agg_acc3.csv    work
python scripts/aggregate.py results/journal results/agg_journal.csv work
python scripts/aggregate.py results/abl     results/agg_abl.csv     work
python scripts/aggregate.py results/census  results/agg_census.csv  work
python scripts/aggregate.py results/acc2    results/agg_acc2.csv    work   # needs scipy (small-count intervals)
python scripts/aggregate.py results/ci      results/agg_ci.csv      work   # needs scipy (coverage table)
python scripts/aggregate.py results/range   results/agg_range.csv   work
python scripts/aggregate.py results/merge   results/agg_merge.csv   work
python scripts/aggregate.py results/pairs   results/agg_pairs.csv   work
python scripts/journal_metric.py
python scripts/plots.py          # -> paper/figs, paper/tables, paper/numbers.tex, results/summary.txt
scripts/build_paper.sh           # paper/paper.pdf and paper/appendix.pdf (acmart; needs a TeX installation)
python scripts/check_paper.py    # static checks of the LaTeX sources (macros, references, environments, page limit)
```

`ext/` (upstream checkouts, JDK, class files), `data/` and `work/` are not tracked.

The scripts run on Linux and macOS (Apple silicon included). Accuracy results depend only on the data,
the workloads and the seeds, not on the hardware; across platforms they agree statistically but need not be
bit-identical, because a few decisions go through the platform's `log`/`pow`. Timing results
(`scripts/run_timing.sh`) must be taken on an otherwise idle machine.
