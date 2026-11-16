# Datasets

The streams are not stored in this repository: the real ones are third-party data with their own terms of
use, and together they take about 3 GB after conversion. This page lists every stream used in the paper:
what it is, where to get it, how it is converted, and checksums of the files behind the published results.
Synthetic streams, query workloads and raw experiment logs are attached to the GitHub release
[`artifacts-2026-11-15`](https://github.com/roy-pyke/stream_processing_algorithm/releases/tag/artifacts-2026-11-15)
(see [Release assets](#release-assets)).

All converted streams use one binary format (`scripts/convert.py`, `src/common.hpp`):
`"STRM" | u32 p | u64 n | n*p u32`, row-major, every attribute mapped to a non-negative 32-bit integer
(equality semantics preserved). Put downloads in `data/raw/` and converted files in `data/`.

## Real streams

| Stream | Type | Records | Attributes | Source (public) | Terms |
|---|---|---:|---:|---|---|
| Taxi | relational table (trip records), stream order = pick-up time | 7,696,617 | 12 | NYC TLC Yellow Taxi trips, January 2019: <https://d37ci6vzurychx.cloudfront.net/trip-data/yellow_tripdata_2019-01.parquet> (index: <https://www.nyc.gov/site/tlc/about/tlc-trip-record-data.page>) | NYC Open Data terms of use |
| MAWI | network packet trace (IPv4 headers), 15 min, trans-Pacific backbone | 11,331,916 | 10 | MAWI archive, samplepoint-F, 2010-01-01 14:00: <http://mawi.wide.ad.jp/mawi/samplepoint-F/2010/201001011400.dump.gz> (page: <http://mawi.wide.ad.jp/mawi/samplepoint-F/2010/201001011400.html>) | MAWI / WIDE archive terms; addresses are anonymised by the archive |
| DMV | relational table (registrations), stream order = file order | 12,540,479 | 15 | New York State vehicle, snowmobile and boat registrations: <https://data.ny.gov/Transportation/Vehicle-Snowmobile-and-Boat-Registrations/w4pv-hbkt>, CSV export <https://data.ny.gov/api/views/w4pv-hbkt/rows.csv?accessType=DOWNLOAD> | NY Open Data terms of use |
| Census | relational table (1990 US census extract), used for attribute scaling | 2,458,285 | 68 | UCI ML Repository, US Census Data (1990): <https://archive.ics.uci.edu/dataset/116/us+census+data+1990>, file <https://archive.ics.uci.edu/static/public/116/us+census+data+1990.zip> | CC BY 4.0 |

Attributes:

- **Taxi**: VendorID, passenger_count, RatecodeID, PULocationID, DOLocationID, payment_type, pick-up hour,
  pick-up weekday, trip distance (0.1-mile bins), fare ($0.50 bins), tip ($1 bins), duration (minutes, capped at 600).
- **MAWI**: source and destination IPv4 address, source and destination port, protocol, packet length,
  source and destination /16 prefix, TTL, TCP flags (`src/pcap2bin.cpp`).
- **DMV**: record type, registration class, city, state, ZIP, county, model year, make, body type, fuel type,
  registration valid date, colour, scofflaw / suspension / revocation indicators. The VIN and the sparse
  weight and passenger columns are dropped; values are dictionary-coded in order of first appearance.
- **Census**: all 68 categorical attributes of `USCensus1990.data.txt` (the `caseid` column is dropped);
  the attribute-scaling experiment uses the first p of them.

Conversion:

```bash
python scripts/convert.py taxi   data/raw/yellow_2019-01.parquet               data/taxi.bin
python scripts/convert.py census data/raw/census_x/USCensus1990.data.txt     data/census.bin   # unzip census1990.zip first
python scripts/convert.py dmv    data/raw/dmv/dmv.csv                          data/dmv.bin
gunzip -k data/raw/mawi_201001011400.dump.gz && bin/pcap2bin data/raw/mawi_201001011400.dump data/mawi.bin
```

Checksums (MD5) of the files used for the published results:

| File | MD5 |
|---|---|
| `data/raw/yellow_2019-01.parquet` (110 MB) | `5d72c723f22a934780ab4afe81903497` |
| `data/raw/mawi_201001011400.dump.gz` (309 MB) | `344b5e9cb8413dddbc50ec2bd505c195` |
| `data/raw/dmv/dmv.csv` (1.38 GB) | `6e57f90976de043f34cf5341d7295516` |
| `data/raw/census1990.zip` (169 MB) | `1c75a465195592ce335705e904653e47` |
| `data/taxi.bin` | `2c5a25834c3573172bf934d0928e357a` |
| `data/mawi.bin` | `7491f634d22289f00992a4cd3aec9570` |
| `data/dmv.bin` | `f319b8803082fa36f840de9d99af65e3` |
| `data/census.bin` | `e69071ef0ba0b78716727c035aa62217` |

The DMV table on data.ny.gov is updated monthly, so a download made later will not match the checksum above
(it was retrieved on 2026-11-15). The exact file is kept as release
[`data-dmv-2026-11`](https://github.com/roy-pyke/stream_processing_algorithm/releases/tag/data-dmv-2026-11)
(`dmv-2026-11.csv.xz`, 148 MB, SHA-256 `de1e0631…2e2513d82be`; the OPEN-NY Terms of Use allow sharing):

```bash
xz -dc dmv-2026-11.csv.xz > data/raw/dmv/dmv.csv     # MD5 6e57f90976de043f34cf5341d7295516
```

The Taxi file has also been republished by the TLC in the past; it is not mirrored here, because the NYC
terms of use do not grant a license to redistribute it. If a checksum differs, the results will differ in
detail but not in kind.

## Streams of the original evaluation (not available to us)

The OmniSketch paper evaluates on two further real streams. Neither is used in our results yet: SNMP
needs a free IEEE DataPort account (#11), CAIDA an approved application (#10).

| Stream | Type | Size in the OmniSketch paper | Access |
|---|---|---|---|
| SNMP | SNMP polling records of the Dartmouth College wireless network, fall 2003 | 8.2 M records, 11 attributes | CRAWDAD `dartmouth/campus` on IEEE DataPort (free account required), file `fall03.tar.gz` of `dartmouth/campus/snmp`: <https://ieee-dataport.org/open-access/crawdad-dartmouthcampus-v-2009-09-09>, doi:10.15783/C7F59T |
| CAIDA | anonymised backbone packet traces (CAIDA UCSD Anonymized Internet Traces) | 109 M records, 6 attributes | by application (academic and US government researchers), usually processed in a few business days: <https://www.caida.org/catalog/datasets/request_user_info_forms/passive_dataset_request/>; catalogue: <https://www.caida.org/catalog/datasets/passive_dataset/>. Subject to the CAIDA Acceptable Use Agreement. The paper does not state which year's trace it uses. |

The pipeline for both is ready (#12); `scripts/run_new_traces.sh` converts, generates the workloads, runs
the accuracy matrix and aggregates as soon as the raw data is in place:

```bash
mkdir -p data/raw/snmp && tar -xzf fall03.tar.gz -C data/raw/snmp     # day directories 031101 ... *.snmp.gz
mkdir -p data/raw/caida && mv <trace files>.pcap.gz data/raw/caida/  # converted in name order
scripts/run_new_traces.sh 6
```

- **SNMP** (`scripts/convert.py snmp`) follows the authors' loader (`omniVLDB.RecordSNMP`,
  `ProcessedStreamLoaderGeneric`, `Dataset.SNMPDataset`): the days 031101-031119 (their file condition
  `03110_OR_03111`); a record is a `sys` line with the fields of the last well-formed `if` line before the
  next `sys` line; records missing a non-optional field are dropped; stream order = timestamp. The 11
  attributes are those of their experiments (`CompareBaselines`): ifOutOctets, ifOutErrors, ifInDiscards,
  AP, ifDescr, ifOutUcastPkts, ifInUcastPkts, ifInOctets, ifInErrors, awcFtClientSTASelf, ifOutDiscards.
- **CAIDA** (`bin/pcap2bin --caida out.bin trace...`): the 6 attributes of their experiments, dstIP,
  protocol, dst/16, srcIP, src&0xffff, src/16. The authors encode an address by deleting its dots, which
  merges some addresses (1.23.4.5 and 12.3.4.5); we keep exact addresses.

## Synthetic streams

All have 5 million records. They are attached to the release because the authors' generator cannot be seeded.

| File | Generator | Attributes | Distribution | MD5 |
|---|---|---:|---|---|
| `zipf10_5M.bin` | ours, `bin/gen_synth data/zipf10_5M.bin 5000000 5 1.0 10000 1` | 5 | Zipf α=1.0 over 10,000 values, ranks randomly permuted | `82fdc8845ed2df2d518c705c594a07f0` |
| `zipf13_5M.bin` | ours, `... 5 1.3 10000 2` | 5 | Zipf α=1.3 | `5b312733dd299abf6b16ffca707abb95` |
| `zipf15_5M.bin` | ours, `... 5 1.5 10000 3` | 5 | Zipf α=1.5 | `13d59487418153094e2910e58f48e530` |
| `unif_5M.bin` | ours, `... 5 0 10000 4` | 5 | uniform over 10,000 values | `a493edffacb10298843df8bfb41a63e0` |
| `jzipf10_5M.bin` | authors' `ZipfGenerator` (`scripts/gen_upstream_synth.sh`) | 5 | Zipf α=1.0, value k has rank k | `c88e2cd0f414bce716b617653fbf273a` |
| `jzipf13_5M.bin` | authors' generator | 5 | Zipf α=1.3 | `f00a7268daf30d51ee70571e372d2263` |
| `jzipf15_5M.bin` | authors' generator | 5 | Zipf α=1.5 | `b682b9eb92fde4f23a765cc9cc99322e` |
| `jz13p9_5M.bin` | authors' generator (`scripts/run_journal_setting.sh`) | 9 | Zipf α=1.3, journal-version setting | `87679f1c16fe1df2ccc17abcaac16883` |

Our generator is deterministic given its seed (last argument); the authors' generator is not, so the
`jzipf*` and `jz13p9` files can only be obtained from the release.

## Query workloads

`bin/mkworkload` draws the queries and computes exact answers by a full scan; it is deterministic given the
stream and the seed (see the README for the commands). The workload files used for the published results
(`work/*.wl`) are in the release asset `workloads.tar.zst`.

## Release assets

Release [`artifacts-2026-11-15`](https://github.com/roy-pyke/stream_processing_algorithm/releases/tag/artifacts-2026-11-15)
(`SHA256SUMS` included):

| Asset | Contents |
|---|---|
| `workloads.tar.zst` | `work/*.wl` and generator logs for all streams |
| `synthetic-authors-generator.tar.zst` | `data/jzipf{10,13,15}_5M.bin`, `data/jz13p9_5M.bin` |
| `synthetic-ours.tar.zst` | `data/zipf{10,13,15}_5M.bin`, `data/unif_5M.bin` |
| `raw-acc.tar.zst` | per-run logs of the main accuracy matrix (`results/acc/`) |
| `raw-acc3.tar.zst` | triple-table variants (`results/acc3/`) |
| `raw-abl.tar.zst` | ablations (`results/abl/`) |
| `raw-census.tar.zst`, `raw-census_cap.tar.zst` | attribute scaling, sampled table updates |
| `raw-journal.tar.zst` | journal-version setting (`results/journal/`) |
| `raw-hydra.tar.zst` | the authors' Hydra adaptation (`results/hydra/`) |
| `raw-validation.tar.zst` | authors' Java code vs. C++ port, query by query (`results/validation/{java,cpp}/`) |
| `raw-bias.tar.zst` | bias checks, including the first batch with identical seeds that is mentioned in the paper |
| `colab-raw-results-2026-11-13.zip` | raw logs of the first (Colab) run of the initial version |

Unpack in the repository root, e.g. `zstd -dc raw-acc.tar.zst | tar -x`; `scripts/aggregate.py` then
rebuilds the `results/agg_*.csv` files.
