#!/usr/bin/env python3
"""Convert public datasets to the binary stream format used by the benchmarks.

  convert.py taxi   <yellow_tripdata_2019-01.parquet> <out.bin>
  convert.py census <USCensus1990.data.txt>           <out.bin>
  convert.py dmv    <Vehicle__Snowmobile__and_Boat_Registrations.csv> <out.bin>
  convert.py snmp   <directory of the extracted fall03.tar.gz>       <out.bin>

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


def dmv(src, dst):
    """New York State vehicle, snowmobile and boat registrations (data.ny.gov, dataset w4pv-hbkt),
    the "DMV" table of the selectivity-estimation literature. 15 categorical attributes (the VIN and
    the sparse weight/passenger columns are dropped); values are dictionary-coded in order of first
    appearance; the stream order is the file order."""
    import pandas as pd

    names = ["Record Type", "Registration Class", "City", "State", "Zip", "County", "Model Year", "Make", "Body Type",
             "Fuel Type", "Reg Valid Date", "Color", "Scofflaw Indicator", "Suspension Indicator", "Revocation Indicator"]
    codes = [dict() for _ in names]
    parts = [[] for _ in names]
    for chunk in pd.read_csv(src, usecols=names, dtype=str, keep_default_na=False, chunksize=500000):
        for j, name in enumerate(names):
            c, u = pd.factorize(chunk[name], sort=False)
            d = codes[j]
            m = np.empty(len(u), dtype=np.uint32)
            for i, val in enumerate(u):
                x = d.get(val)
                if x is None:
                    x = d[val] = len(d)
                m[i] = x
            parts[j].append(m[c])
    write(dst, [np.concatenate(x) for x in parts], names)


# SNMP fields of a record, in the order of the authors' Parser.attrMap (indices 0..23)
SNMP_FIELDS = ["timestamp", "AP", "sysUpTime", "sysDescr", "ifIndex", "ifDescr", "ifType", "ifSpeed", "ifInOctets",
               "ifInUcastPkts", "ifInErrors", "ifInDiscards", "ifOutOctets", "ifOutUcastPkts", "ifOutErrors",
               "ifOutDiscards", "awcDot11AssociatedStationCount", "awcDot11ReassociatedStationCount",
               "awcDot11RoamedStationCount", "awcDot11DeauthenticateCount", "awcDot11DisassociateCount",
               "awcFtClientSTASelf", "awcFtBridgeSelf", "awcFtRepeaterSelf"]
SNMP_ATTRS = [12, 14, 11, 1, 5, 13, 9, 8, 10, 21, 15]      # the 11 attributes of the OmniSketch paper (CompareBaselines)
SNMP_OPTIONAL = {3, 16, 17, 18, 19, 20}                     # may be missing; every other field must be present
SNMP_DAYS = ("03110", "03111")                              # fileStartCondition "03110_OR_03111": 1-19 Nov 2003


def _jsplit(line):
    """String.split(",") of Java: trailing empty strings are removed."""
    a = line.rstrip("\r\n").split(",")
    while a and a[-1] == "":
        a.pop()
    return a


def snmp_records(path):
    """Records of one .snmp.gz file as the authors' loader builds them (omniVLDB.RecordSNMP and
    ProcessedStreamLoaderGeneric): the first 4 lines and comment lines are skipped; a record starts at
    a "sys" line (4 or 5 fields: sys, timestamp, AP, sysUpTime[, sysDescr]) and takes the interface fields
    of the last well-formed "if" line (23 fields) before the next "sys" line. Strings are kept as strings
    (the authors use String.hashCode; we dictionary-code later)."""
    import gzip

    rec = None
    with gzip.open(path, "rt", errors="replace") as f:
        for k, line in enumerate(f):
            if k < 4 or not line.strip() or line[0] == "#":
                continue
            a = _jsplit(line)
            if not a:
                continue
            if a[0][:1] == "s":
                if rec is not None:
                    yield rec
                rec = None
                if a[0] == "sys" and len(a) in (4, 5):
                    rec = [None] * 24
                    try:
                        rec[0], rec[1], rec[2] = int(a[1]), a[2], int(a[3])
                    except ValueError:
                        rec = None
                        continue
                    if len(a) == 5:
                        rec[3] = a[4]
            elif a[0] == "if" and rec is not None and len(a) == 23:
                try:
                    vals = [int(a[3]), a[4]] + [int(x) for x in a[5:15]] + [int(x) if x != "" else None for x in a[15:22]] + [int(a[22])]
                except ValueError:
                    continue
                rec[4:24] = vals
    if rec is not None:
        yield rec


def snmp(src, dst):
    """Dartmouth campus SNMP trace (CRAWDAD dartmouth/campus, fall03.tar.gz), the "SNMP" stream of the
    OmniSketch paper: the days 031101-031119, records with every non-optional field present, ordered by
    timestamp (stable), the 11 attributes the authors use, dictionary-coded in order of first appearance."""
    import os

    days = []
    for root, dirs, _ in os.walk(src):
        days += [os.path.join(root, d) for d in dirs if d.startswith(SNMP_DAYS)]
    days.sort(key=os.path.basename)
    if not days:
        sys.exit(f"no day directories 03110*/03111* below {src}")
    recs, files, dropped = [], 0, 0
    for day in days:
        for fn in sorted(os.listdir(day)):
            if not fn.endswith(".snmp.gz"):
                continue
            files += 1
            for r in snmp_records(os.path.join(day, fn)):
                if any(r[i] is None for i in range(24) if i not in SNMP_OPTIONAL):
                    dropped += 1
                else:
                    recs.append(r)
    print(f"{len(days)} days, {files} files, {len(recs)} records ({dropped} dropped for missing fields)")
    recs.sort(key=lambda r: r[0])
    cols = []
    for i in SNMP_ATTRS:
        d = {}
        cols.append(np.fromiter((d.setdefault(r[i], len(d)) for r in recs), dtype=np.uint32, count=len(recs)))
    write(dst, cols, [SNMP_FIELDS[i] for i in SNMP_ATTRS])


if __name__ == "__main__":
    {"taxi": taxi, "census": census, "dmv": dmv, "snmp": snmp}[sys.argv[1]](sys.argv[2], sys.argv[3])
