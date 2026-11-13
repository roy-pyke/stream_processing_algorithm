#!/usr/bin/env python3
"""Build a self-contained Colab notebook that reproduces the accuracy experiments.

  make_colab_notebook.py   (writes colab/costrata_experiments.ipynb)

The notebook writes all sources with %%writefile, downloads and converts the public datasets,
runs scripts/run_accuracy.sh and scripts/run_extra.sh (timing skipped: shared VMs are too noisy),
mirrors results to Google Drive and aggregates them.
"""
import json
import os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FILES = ["Makefile", "src/common.hpp", "src/synopses.hpp", "src/bench.cpp", "src/gen_synth.cpp", "src/mkworkload.cpp",
         "src/pcap2bin.cpp", "scripts/run_accuracy.sh", "scripts/run_extra.sh", "scripts/aggregate.py",
         "scripts/convert.py", "colab/done_local.txt"]
MD5 = {}
for line in open(os.path.join(ROOT, "colab/md5_local.txt")):
    h, path = line.split()
    MD5[path] = h

cells = []


def md(text):
    cells.append({"cell_type": "markdown", "metadata": {}, "source": text})


def code(text):
    cells.append({"cell_type": "code", "metadata": {}, "execution_count": None, "outputs": [], "source": text})


md("""# CoStrata experiments (accuracy matrix, ablations, attribute scaling)

**Before running:** `Runtime -> Change runtime type -> CPU`, and enable **High-RAM** (gives 8 vCPUs; GPUs do not help, the
workload is single-threaded C++ run as parallel jobs). Then `Runtime -> Run all`.

* Results are mirrored to `MyDrive/costrata_results/` every two minutes, and the run resumes from there if the
  runtime is recycled (just *Run all* again).
* Jobs already finished on the laptop are skipped if the regenerated data is byte-identical (checked by MD5).
* Update-throughput timings are **not** measured here (shared VMs are too noisy); they are run locally.""")

code("""import os, subprocess
from google.colab import drive
drive.mount('/content/drive')
ROOT = '/content/costrata'
DRIVE = '/content/drive/MyDrive/costrata_results'
for d in ('src', 'scripts', 'colab', 'data/raw', 'work', 'results'):
    os.makedirs(f'{ROOT}/{d}', exist_ok=True)
os.makedirs(f'{DRIVE}/results', exist_ok=True)
os.chdir(ROOT)
print('vCPUs:', os.cpu_count())
print(subprocess.run("grep -m1 'model name' /proc/cpuinfo; free -g | sed -n 2p", shell=True, capture_output=True, text=True).stdout)""")

for f in FILES:
    code(f"%%writefile {f}\n" + open(os.path.join(ROOT, f)).read())

code("""%%bash
set -e
cd /content/costrata
chmod +x scripts/*.sh
make 2>&1 | tail -5
ls -la bin""")

code("""%%bash
set -e
cd /content/costrata
# synthetic streams (generator and seeds identical to the local run)
[ -s data/zipf10_5M.bin ] || ./bin/gen_synth data/zipf10_5M.bin 5000000 5 1.0 10000 1
[ -s data/zipf13_5M.bin ] || ./bin/gen_synth data/zipf13_5M.bin 5000000 5 1.3 10000 2
[ -s data/zipf15_5M.bin ] || ./bin/gen_synth data/zipf15_5M.bin 5000000 5 1.5 10000 3
[ -s data/unif_5M.bin ]   || ./bin/gen_synth data/unif_5M.bin   5000000 5 0   10000 4
# NYC yellow taxi, January 2019
if [ ! -s data/taxi.bin ]; then
  curl -sL -o data/raw/yellow_2019-01.parquet https://d37ci6vzurychx.cloudfront.net/trip-data/yellow_tripdata_2019-01.parquet
  python scripts/convert.py taxi data/raw/yellow_2019-01.parquet data/taxi.bin
fi
# MAWI backbone trace 2010-01-01 14:00 (samplepoint F)
if [ ! -s data/mawi.bin ]; then
  curl -sL -o data/raw/mawi_201001011400.dump.gz http://mawi.wide.ad.jp/mawi/samplepoint-F/2010/201001011400.dump.gz
  ./bin/pcap2bin data/raw/mawi_201001011400.dump.gz data/mawi.bin
fi
# US Census 1990 (UCI)
if [ ! -s data/census.bin ]; then
  curl -sL -o data/raw/census1990.zip "https://archive.ics.uci.edu/static/public/116/us+census+data+1990.zip"
  unzip -o -q data/raw/census1990.zip USCensus1990.data.txt -d data/raw
  python scripts/convert.py census data/raw/USCensus1990.data.txt data/census.bin | tail -1
fi
for ds in zipf10_5M zipf13_5M zipf15_5M unif_5M; do [ -s work/$ds.wl ] || ./bin/mkworkload data/$ds.bin work/$ds.wl 200 5 7 2>&1 | tail -1; done
for ds in taxi mawi; do [ -s work/$ds.wl ] || ./bin/mkworkload data/$ds.bin work/$ds.wl 200 6 7 2>&1 | tail -1; done
ls -la data/*.bin work/*.wl""")

code("""import hashlib
LOCAL_MD5 = """ + json.dumps(MD5, indent=1) + """
def md5(path):
    h = hashlib.md5()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 22), b''):
            h.update(chunk)
    return h.hexdigest()
same = True
for path, want in LOCAL_MD5.items():
    got = md5(path)
    ok = got == want
    if not ok and 'census' not in path:
        same = False
    print(('OK      ' if ok else 'DIFFERS ') + path)
SKIP_LIST = f'{ROOT}/colab/done_local.txt' if same else ''
print('\\nData identical to the local run -> skipping the 61 jobs already finished there.' if same
      else '\\nData differs from the local run -> running the full matrix here (local partial results must then be discarded).')
open(f'{DRIVE}/data_identical.txt', 'w').write(str(same))""")

code("""# restore earlier progress from Drive (if the runtime was recycled), then launch everything in the background
subprocess.run(f'rsync -a {DRIVE}/results/ {ROOT}/results/', shell=True)
env = dict(os.environ, SKIP_LIST=SKIP_LIST, SKIP_TIMING='1')
J = os.cpu_count()
log = open(f'{ROOT}/run.out', 'w')
proc = subprocess.Popen(['bash', '-c', f'bash scripts/run_accuracy.sh {J} && bash scripts/run_extra.sh {J}'],
                        cwd=ROOT, env=env, stdout=log, stderr=subprocess.STDOUT)
print('started, pid', proc.pid, 'parallel jobs:', J)""")

code("""# progress monitor: mirrors results to Drive every 2 minutes; keep this cell running until it prints DONE
import glob, time
def complete(d, need):
    n = 0
    for f in glob.glob(f'{ROOT}/results/{d}/*.log'):
        with open(f) as fh:
            n += sum(1 for _ in fh) >= need
    return n
t0 = time.time()
while True:
    done = proc.poll() is not None
    subprocess.run(f'rsync -a {ROOT}/results/ {DRIVE}/results/', shell=True)
    print(f'[{(time.time() - t0) / 60:6.1f} min] accuracy {complete("acc", 15)}/210 (61 may be skipped) | '
          f'ablation {complete("abl", 13)}/30 | census {complete("census", 8)}/15', flush=True)
    if done:
        break
    time.sleep(120)
print('DONE, exit code', proc.returncode)
print(open(f'{ROOT}/run.out').read()[-2000:])""")

code("""%%bash
cd /content/costrata
python scripts/aggregate.py results/acc results/agg_acc.csv
python scripts/aggregate.py results/abl results/agg_abl.csv
python scripts/aggregate.py results/census results/agg_census.csv
cat results/*/*.err 2>/dev/null | head -20
D=/content/drive/MyDrive/costrata_results
cp results/agg_*.csv $D/
cat results/acc/*.log results/abl/*.log results/census/*.log > $D/run_lines.csv
md5sum data/*.bin work/*.wl > $D/md5_colab.txt
ls -la $D""")

md("""When the last cell has finished, `MyDrive/costrata_results/` contains `agg_acc.csv`, `agg_abl.csv`, `agg_census.csv`
(a few hundred KB in total) plus the raw per-query output under `results/`. Only the three `agg_*.csv` files,
`run_lines.csv` and `data_identical.txt` are needed to finish the paper.""")

nb = {"cells": [dict(c, source=c["source"].splitlines(keepends=True)) for c in cells],
      "metadata": {"colab": {"name": "costrata_experiments.ipynb", "provenance": []},
                   "kernelspec": {"name": "python3", "display_name": "Python 3"}, "language_info": {"name": "python"}},
      "nbformat": 4, "nbformat_minor": 0}
out = os.path.join(ROOT, "colab/costrata_experiments.ipynb")
json.dump(nb, open(out, "w"), indent=1)
print("wrote", out, os.path.getsize(out), "bytes,", len(cells), "cells")
