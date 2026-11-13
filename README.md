# CoStrata: COUNT queries with ad-hoc conjunctive predicates on multi-attribute streams

Code, aggregated results and paper draft comparing **OmniSketch** (PVLDB 17(3)) against
row sampling (reservoir, USS) and **CoStrata**, a coordinated hash-stratified row sample
with counter calibration.

## Layout

| Path | Contents |
|---|---|
| `src/` | C++17 benchmark (`bench`), synthetic generator, workload generator, pcap converter; synopses in `synopses.hpp` |
| `scripts/` | dataset conversion, experiment drivers (`run_*.sh`), aggregation and plotting |
| `colab/costrata_experiments.ipynb` | end-to-end notebook: download data, build, run the accuracy matrix |
| `results/` | aggregated results (`agg_*.csv`, `bias.json`, timing logs); raw per-run logs are not tracked |
| `paper/` | LaTeX draft (`paper.tex`, `results.tex`), generated `figs/`, `tables/`, `numbers.tex` |

## Build

```bash
make                      # builds bin/bench bin/gen_synth bin/mkworkload bin/pcap2bin
python -m venv .venv && .venv/bin/pip install -r requirements.txt
```

## Data (not tracked, ~2.7 GB)

Real datasets (place downloads in `data/raw/`):

- NYC Yellow Taxi 2019-01: https://d37ci6vzurychx.cloudfront.net/trip-data/yellow_tripdata_2019-01.parquet
- US Census 1990 (UCI): https://archive.ics.uci.edu/static/public/116/us+census+data+1990.zip
- MAWI trace: http://mawi.wide.ad.jp/mawi/samplepoint-F/2010/201001011400.dump.gz

```bash
python scripts/convert.py taxi   data/raw/yellow_2019-01.parquet   data/taxi.bin
python scripts/convert.py census data/raw/USCensus1990.data.txt    data/census.bin
bin/pcap2bin <decompressed MAWI dump> data/mawi.bin

bin/gen_synth data/unif_5M.bin   5000000 5 0   10000 4
bin/gen_synth data/zipf10_5M.bin 5000000 5 1.0 10000 1
bin/gen_synth data/zipf13_5M.bin 5000000 5 1.3 10000 2
bin/gen_synth data/zipf15_5M.bin 5000000 5 1.5 10000 3
```

Query workloads (`work/*.wl`) are generated with `bin/mkworkload`; the exact commands are in
the Colab notebook.

## Reproduce

```bash
scripts/run_accuracy.sh 8      # accuracy matrix (5 seeds, results are seed-deterministic)
scripts/run_timing.sh          # update/query throughput
scripts/run_bias.sh            # unbiasedness check
python scripts/aggregate.py    # -> results/agg_*.csv
python scripts/plots.py        # -> paper/figs, paper/tables, paper/numbers.tex
```
