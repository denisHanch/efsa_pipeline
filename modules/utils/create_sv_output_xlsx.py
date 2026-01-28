#!/usr/bin/env python3
"""
Same description as before, but outputs CSV files instead of a single XLSX.

Outputs:
  <outdir>/Insertions.csv
  <outdir>/Deletions.csv
  <outdir>/Replacements.csv
  <outdir>/Inversions.csv
  <outdir>/Translocations.csv
  <outdir>/Other.csv (if needed)
"""
from __future__ import annotations

import argparse
import os
import re
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd


# ----------------------------
# SV type standardization
# ----------------------------

TAB_BY_TYPE = {
    "INS": "Insertions",
    "DEL": "Deletions",
    "RPL": "Replacements",
    "INV": "Inversions",
    "TRA": "Translocations",
}

INFO_MAP = {
    "INS": "INS",
    "DEL": "DEL",
    "INV": "INV",
    "TRA": "TRA",
    "TRANS": "TRA",
    "BND": "TRA",
    "CTX": "TRA",
    "DUP": "RPL",
    "DUP:TANDEM": "RPL",
    "DUP:INT": "RPL",
    "RPL": "RPL",
    "REPL": "RPL",
    "SUB": "RPL",
    "SNV": "RPL",
}

ASM_PREFIX_MAP = [
    (r"^DEL", "DEL"),
    (r"^INS", "INS"),
    (r"^INV", "INV"),
    (r"^TRANS", "TRA"),
    (r"^TRA", "TRA"),
    (r"^BND", "TRA"),
    (r"^DUP", "RPL"),
    (r"^CPG", "INS"),
    (r"^CPL", "DEL"),
    (r"^SYN", "RPL"),
]


def _to_int(x):
    try:
        if pd.isna(x):
            return None
        return int(float(x))
    except Exception:
        return None


def _to_float(x):
    try:
        if pd.isna(x):
            return None
        return float(x)
    except Exception:
        return None


def standardize_type(raw_svtype, info_svtype, source):
    if info_svtype is not None:
        key = str(info_svtype).strip().upper()
        if key in INFO_MAP:
            return INFO_MAP[key]

    raw_u = str(raw_svtype).upper()

    for token, mapped in INFO_MAP.items():
        if token in raw_u:
            return mapped

    if source == "asm":
        for pat, mapped in ASM_PREFIX_MAP:
            if re.search(pat, raw_u):
                return mapped

    return "RPL"


# ----------------------------
# Data models
# ----------------------------

@dataclass
class Record:
    source: str
    chrom: str
    start: int
    end: int
    std_type: str
    raw_svtype: str
    info_svtype: Optional[str] = None
    supporting_reads: Optional[int] = None
    score: Optional[float] = None
    copy_number: Optional[int] = None


@dataclass
class EventCluster:
    chrom: str
    std_type: str
    start: int
    end: int
    members: List[Record]


# ----------------------------
# Clustering
# ----------------------------

def intervals_overlap(a_start, a_end, b_start, b_end, tol):
    overlap = (b_start <= a_end + tol) and (b_end >= a_start - tol)
    if not overlap:
        return False
    return abs(a_start - b_start) <= tol or abs(a_end - b_end) <= tol


def cluster_records(records, tol):
    clusters = []
    key_to_recs = {}

    for r in records:
        key_to_recs.setdefault((r.chrom, r.std_type), []).append(r)

    for (chrom, std_type), recs in key_to_recs.items():
        recs = sorted(recs, key=lambda r: (r.start, r.end))
        cur = [recs[0]]

        for r in recs[1:]:
            last = cur[-1]
            if intervals_overlap(last.start, last.end, r.start, r.end, tol):
                cur.append(r)
            else:
                clusters.append(
                    EventCluster(
                        chrom,
                        std_type,
                        min(x.start for x in cur),
                        max(x.end for x in cur),
                        cur,
                    )
                )
                cur = [r]

        clusters.append(
            EventCluster(
                chrom,
                std_type,
                min(x.start for x in cur),
                max(x.end for x in cur),
                cur,
            )
        )

    return sorted(clusters, key=lambda c: (c.std_type, c.chrom, c.start))


def choose_best_record(recs):
    if not recs:
        return None

    def key(r):
        return (
            r.supporting_reads or -1,
            r.score or -1,
            -(r.end - r.start + 1),
        )

    return max(recs, key=key)


# ----------------------------
# IO
# ----------------------------

def read_tsv(path):
    return pd.read_csv(path, sep="\t", dtype=str, comment="#")


def load_records(path, source):
    if not path:
        return []

    df = read_tsv(path)
    out = []

    for _, row in df.iterrows():
        start = _to_int(row["start"])
        end = _to_int(row["end"])
        if start is None or end is None:
            continue
        if start > end:
            start, end = end, start

        std = standardize_type(row["svtype"], row.get("info_svtype"), source)

        copy_number = None
        if source == "short":
            copy_number = _to_int(row.get("RDCN"))

        out.append(
            Record(
                source,
                row["chrom"],
                start,
                end,
                std,
                row["svtype"],
                row.get("info_svtype"),
                _to_int(row.get("supporting_reads")),
                _to_float(row.get("score")),
                _to_int(row.get("RDCN")),
                
            )
        )
    return out


def build_output_table(clusters):
    rows = []
    counters = {k: 0 for k in TAB_BY_TYPE}

    for c in clusters:
        counters[c.std_type] += 1
        eid = f"{c.std_type}_{counters[c.std_type]}"

        by_src = {"asm": [], "long": [], "short": []}
        for m in c.members:
            by_src[m.source].append(m)

        asm = choose_best_record(by_src["asm"])
        lng = choose_best_record(by_src["long"])
        sht = choose_best_record(by_src["short"])

        rows.append({
            "event_id": eid,
            "chrom": c.chrom,
            "std_svtype": c.std_type,
            "event_start": c.start,
            "event_end": c.end,
            "event_length_bp": c.end - c.start + 1,

            "asm_start": asm.start if asm else np.nan,
            "asm_end": asm.end if asm else np.nan,
            "asm_svtype_raw": asm.raw_svtype if asm else "",
            "asm_score" : 1,

            "long_start": lng.start if lng else np.nan,
            "long_end": lng.end if lng else np.nan,
            "long_svtype_raw": lng.raw_svtype if lng else "",
            "long_info_svtype": lng.info_svtype if lng else "",
            "long_score": lng.score if lng else np.nan,
            "long_supporting_reads": lng.supporting_reads if lng else np.nan,

            "short_start": sht.start if sht else np.nan,
            "short_end": sht.end if sht else np.nan,
            "short_svtype_raw": sht.raw_svtype if sht else "",
            "short_info_svtype": sht.info_svtype if sht else "",
            "short_score": sht.score if sht else np.nan,
            "short_supporting_reads": sht.supporting_reads if sht else np.nan,
            "short_reads_copy_number_estimate": (sht.copy_number if sht else np.nan)
        })

    return pd.DataFrame(rows)


def write_csv_tables(df, outdir):
    os.makedirs(outdir, exist_ok=True)

    for std_type, name in TAB_BY_TYPE.items():
        sub = df[df["std_svtype"] == std_type]
        if not sub.empty:
            path = os.path.join(outdir, f"{name}.csv")
            sub.sort_values(["chrom", "event_start", "event_end"]).to_csv(path, index=False)

    other = df[~df["std_svtype"].isin(TAB_BY_TYPE)]
    if not other.empty:
        other.to_csv(os.path.join(outdir, "Other.csv"), index=False)


# ----------------------------
# Main
# ----------------------------

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--asm")
    p.add_argument("--long", dest="long_reads")
    p.add_argument("--short", dest="short_reads")
    p.add_argument("--out", required=True, help="Output directory")
    p.add_argument("--tol", type=int, default=10)
    args = p.parse_args()

    records = []
    records += load_records(args.asm, "asm")
    records += load_records(args.long_reads, "long")
    records += load_records(args.short_reads, "short")

    clusters = cluster_records(records, args.tol)
    df = build_output_table(clusters)
    mask = df["long_info_svtype"].notna() & (df["long_info_svtype"] != "")
    df.loc[mask, "long_score"] = df.loc[mask, "long_score"].fillna(0)
    write_csv_tables(df, args.out)

    print(f"Wrote CSV tables to: {args.out}")


if __name__ == "__main__":
    main()
