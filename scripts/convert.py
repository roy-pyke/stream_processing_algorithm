#!/usr/bin/env python3
"""Convert public datasets to the binary stream format used by the benchmarks.

  convert.py taxi   <yellow_tripdata_2019-01.parquet> <out.bin>
  convert.py census <USCensus1990.data.txt>           <out.bin>

Output: "STRM" | u32 p | u64 n | n*p u32 (row-major). All attributes are mapped
to non-negative 32-bit integers (equality semantics are preserved).
"""
import struct
import sys

import numpy as np


def write(path, cols, names):
    n, p = len(cols[0]), len(cols)
    arr = np.empty((n, p), dtype=np.uint32)
    for j, c in enumerate(cols):
        arr[:, j] = c
        print(f"  attr {j:2d} {names[j]:24s} distinct={len(np.unique(c))}")
    with open(path, "wb") as f:
        f.write(b"STRM")
        f.write(struct.pack("<IQ", p, n))
        arr.tofile(f)
    print(f"wrote {path}: n={n} p={p}")


def taxi(src, dst):
    import pyarrow.parquet as pq

    t = pq.read_table(src)
    g = lambda name: t.column(name).to_numpy(zero_copy_only=False)
    pick = g("tpep_pickup_datetime").astype("datetime64[s]").astype(np.int64)
    drop = g("tpep_dropoff_datetime").astype("datetime64[s]").astype(np.int64)
    dur_min = np.clip((drop - pick) // 60, 0, 600)
    hour = (pick // 3600) % 24
    dow = ((pick // 86400) + 4) % 7
    nz = lambda x: np.nan_to_num(x, nan=0.0)
    cols = [
        g("VendorID"),
        np.clip(nz(g("passenger_count")), 0, 15),
        np.clip(nz(g("RatecodeID")), 0, 120),
        g("PULocationID"),
        g("DOLocationID"),
        g("payment_type"),
        hour,
        dow,
        np.clip(np.round(nz(g("trip_distance")) * 10), 0, 5000),   # 0.1 mile bins
        np.clip(np.round(nz(g("fare_amount")) * 2), 0, 2000),      # $0.50 bins
        np.clip(np.round(nz(g("tip_amount"))), 0, 200),            # $1 bins
        dur_min,
    ]
    names = ["VendorID", "passenger_count", "RatecodeID", "PULocationID", "DOLocationID", "payment_type",
             "pickup_hour", "pickup_dow", "trip_distance_x10", "fare_x2", "tip_dollars", "duration_min"]
    # stream order = pickup time
    order = np.argsort(pick, kind="stable")
    write(dst, [np.asarray(c, dtype=np.int64)[order].astype(np.uint32) for c in cols], names)


def census(src, dst):
    with open(src) as f:
        names = f.readline().strip().split(",")[1:]
        data = np.loadtxt(f, delimiter=",", dtype=np.int64)
    data = data[:, 1:]  # drop caseid
    data = data - data.min(axis=0)
    write(dst, [data[:, j].astype(np.uint32) for j in range(data.shape[1])], names)


if __name__ == "__main__":
    {"taxi": taxi, "census": census}[sys.argv[1]](sys.argv[2], sys.argv[3])
