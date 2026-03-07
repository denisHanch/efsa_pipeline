#!/usr/bin/env python3
"""
Merge SV summaries from:
  1) assembly/syri summary (required columns: chrom, start, end, svtype)
  2) long-read SV summary (required columns: chrom, start, end, svtype; optional: info_svtype, supporting_reads, score) both pacbio and ont
  3) short-read SV summary (required columns: chrom, start, end, svtype; optional: info_svtype, supporting_reads, score)

Outputs a folder with single CSV for each type of variations:
  Insertions, Deletions, Replacements, Inversions, Translocations

Events are clustered by (chrom, standardized_svtype) using interval overlap with an optional tolerance (bp),
PLUS breakpoint proximity (start or end must be within tol bp). This avoids merging very large calls with
many unrelated smaller calls.

Within each event, at most one record per source is chosen (best by supporting_reads then score).

Usage:
    python create_sv_output.py --asm assembly.tsv --long_ont long_ont.tsv --long_pacbio long_pb.tsv --short short.tsv --out outdir --tol 10

Any of the inputs may be omitted; the script will process whichever of the assembly, long_ont,
long_pacbio or short tables are provided. Long-read inputs (ONT and PacBio) are treated
the same for clustering and merged into the "long" source.
"""
from __future__ import annotations

import argparse
import os
import re
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import pandas as pd
import numpy as np


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


# Data models


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
    supporting_methods: Optional[str] = None



@dataclass
class EventCluster:
    chrom: str
    std_type: str
    start: int
    end: int
    members: List[Record]
    # sequence of overlap percentages (floats, 0-100) collected during daisy-chain merges
    percentage_overlaps: List[float] = field(default_factory=list)


# Clustering

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
        cur_percs: List[float] = []

        for r in recs[1:]:
            last = cur[-1]
            if intervals_overlap(last.start, last.end, r.start, r.end, tol):
                overlap_len = max(0, min(last.end, r.end) - max(last.start, r.start) + 1)
                len_last = last.end - last.start + 1
                len_r = r.end - r.start + 1
                denom = max(len_last, len_r) if max(len_last, len_r) > 0 else 1
                pct = (overlap_len / denom) * 100 if denom else 0.0
                cur.append(r)
                cur_percs.append(pct)
            else:
                clusters.append(
                    EventCluster(
                        chrom,
                        std_type,
                        min(x.start for x in cur),
                        max(x.end for x in cur),
                        cur,
                        percentage_overlaps=cur_percs,
                    )
                )
                cur = [r]
                cur_percs = []

        clusters.append(
            EventCluster(
                chrom,
                std_type,
                min(x.start for x in cur),
                max(x.end for x in cur),
                cur,
                percentage_overlaps=cur_percs,
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



def read_tsv(path):
    return pd.read_csv(path, sep="\t", dtype=str, comment="#")


def annotate_event_relationships(clusters):
    """
    Robust overlap classification using connected components.

    For each (chrom, std_type):
        - Build overlap graph
        - Extract connected components
        - Within each component:
            - Longest interval = anchor
            - Others classified relative to anchor
            - If anchor has no other variants, label as 'isolated'
    """

    grouped = {}
    for c in clusters:
        grouped.setdefault((c.chrom, c.std_type), []).append(c)

    relationship = {}

    for key, group in grouped.items():

        group = sorted(group, key=lambda c: (c.start, c.end))

        adj = {id(c): [] for c in group}

        for i in range(len(group)):
            for j in range(i + 1, len(group)):
                c1 = group[i]
                c2 = group[j]

                if c1.start <= c2.end and c1.end >= c2.start:
                    adj[id(c1)].append(c2)
                    adj[id(c2)].append(c1)

        visited = set()

        for c in group:
            if id(c) in visited:
                continue

            stack = [c]
            component = []

            while stack:
                node = stack.pop()
                if id(node) in visited:
                    continue
                visited.add(id(node))
                component.append(node)
                stack.extend(adj[id(node)])

            if len(component) == 1:
                relationship[id(component[0])] = "isolated"
                continue

            anchor = max(component, key=lambda x: (x.end - x.start))

            for member in component:
                if member is anchor:
                    if len(component) == 1:
                        relationship[id(member)] = "isolated"
                    else:
                        relationship[id(member)] = f"contains_{len(component)-1}_variants"
                elif member.start >= anchor.start and member.end <= anchor.end:
                    relationship[id(member)] = (
                        f"nested_in_{anchor.std_type}_{anchor.start}_{anchor.end}"
                    )
                else:
                    relationship[id(member)] = (
                        f"overlapping_with_{anchor.std_type}_{anchor.start}_{anchor.end}"
                    )

    return relationship


def load_records(path, source):
    """Load records from a TSV file into Record objects.

    - path: path to tsv (str or None). If falsy, returns empty list.
    - source: one of 'asm', 'long_ont', 'long_pacbio', 'long' or 'short'.
      Sources starting with 'long' are treated as long-read sources and stored
      separately (e.g. 'long_ont' vs 'long_pacbio').
    """
    if not path:
        return []

    try:
        df = read_tsv(path)
    except Exception:
        return []

    out = []

    # defensive: accept files that may be missing optional columns
    for _, row in df.iterrows():
        chrom = row.get("chrom") or row.get("#chrom")
        raw_svtype = row.get("svtype")
        if chrom is None or raw_svtype is None:
            # malformed row -> skip
            continue

        start = _to_int(row.get("start"))
        end = _to_int(row.get("end"))
        if start is None or end is None:
            continue
        if start > end:
            start, end = end, start

        std = standardize_type(raw_svtype, row.get("info_svtype"), source)

        # treat any long* sources as long for supporting_methods
        supporting_methods = row.get("supporting_methods") if str(source).startswith("long") else None
        copy_number = _to_int(row.get("RDCN")) if source == "short" else None

        out.append(
            Record(
                source,
                chrom,
                start,
                end,
                std,
                raw_svtype,
                row.get("info_svtype"),
                _to_int(row.get("supporting_reads")),
                _to_float(row.get("score")),
                copy_number,
                supporting_methods,
            )
        )

    return out


def build_output_table(clusters, relationship_map):
    rows = []
    counters = {k: 0 for k in TAB_BY_TYPE}
    
    for c in clusters:
        counters[c.std_type] += 1
        eid = f"{c.std_type}_{counters[c.std_type]}"

        by_src = {"asm": [], "long_ont": [], "long_pacbio": [], "short": []}
        for m in c.members:
            src = m.source
            if src not in by_src:
                by_src.setdefault(src, []).append(m)
            else:
                by_src[src].append(m)

        asm = choose_best_record(by_src.get("asm", []))
        long_ont = choose_best_record(by_src.get("long_ont", []))
        long_pacbio = choose_best_record(by_src.get("long_pacbio", []))
        sht = choose_best_record(by_src.get("short", []))
        pipelines_confirmed = sum(x is not None for x in (asm, long_ont, long_pacbio, sht))

        # format percentage overlaps collected during clustering into a human-friendly string
        percs = getattr(c, "percentage_overlaps", []) or []
        if percs:
            def fmt(p):
                try:
                    pf = float(p)
                except Exception:
                    return str(p)
                if pf.is_integer():
                    return f"{int(pf)}%"
                return f"{pf:.2f}%"

            pct_str = ", ".join(fmt(p) for p in percs)
        else:
            pct_str = ""

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

            "long_ont_start": long_ont.start if long_ont else np.nan,
            "long_ont_end": long_ont.end if long_ont else np.nan,
            "long_ont_svtype_raw": long_ont.raw_svtype if long_ont else "",
            "long_ont_info_svtype": long_ont.info_svtype if long_ont else "",
            "long_ont_score": long_ont.score if long_ont else np.nan,
            "long_ont_supporting_reads": long_ont.supporting_reads if long_ont else np.nan,
            "long_ont_supporting_methods": long_ont.supporting_methods if long_ont else np.nan,

            "long_pacbio_start": long_pacbio.start if long_pacbio else np.nan,
            "long_pacbio_end": long_pacbio.end if long_pacbio else np.nan,
            "long_pacbio_svtype_raw": long_pacbio.raw_svtype if long_pacbio else "",
            "long_pacbio_info_svtype": long_pacbio.info_svtype if long_pacbio else "",
            "long_pacbio_score": long_pacbio.score if long_pacbio else np.nan,
            "long_pacbio_supporting_reads": long_pacbio.supporting_reads if long_pacbio else np.nan,
            "long_pacbio_supporting_methods": long_pacbio.supporting_methods if long_pacbio else np.nan,

            "short_start": sht.start if sht else np.nan,
            "short_end": sht.end if sht else np.nan,
            "short_svtype_raw": sht.raw_svtype if sht else "",
            "short_info_svtype": sht.info_svtype if sht else "",
            "short_score": sht.score if sht else np.nan,
            "short_supporting_reads": sht.supporting_reads if sht else np.nan,
            "short_reads_copy_number_estimate": (sht.copy_number if sht else np.nan),
            "percentage_overlap": pct_str,

            "support_score": pipelines_confirmed,
            "region_relationship": relationship_map.get(id(c), "isolated"),
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


# Main

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--asm")
    p.add_argument("--long_ont", dest="long_ont")
    p.add_argument("--long_pacbio", dest="long_pacbio")
    p.add_argument("--short", dest="short_reads")
    p.add_argument("--out", required=True, help="Output directory")
    p.add_argument("--tol", type=int, default=10)
    args = p.parse_args()
    records = []
    records += load_records(args.asm, "asm")

    # long inputs: accept both ONT and PacBio; also support legacy --long
    records += load_records(getattr(args, "long_reads", None), "long")
    records += load_records(getattr(args, "long_ont", None), "long_ont")
    records += load_records(getattr(args, "long_pacbio", None), "long_pacbio")

    records += load_records(args.short_reads, "short")

    if not records:
        os.makedirs(args.out, exist_ok=True)
        print("No valid input records found; created output directory and exiting.")
        return

    clusters = cluster_records(records, args.tol)

    if not clusters:
        df = pd.DataFrame()
    else:
        relationship_map = annotate_event_relationships(clusters)
        df = build_output_table(clusters, relationship_map)

    # fix scores for long_ont and long_pacbio when info_svtype present
    if "long_ont_info_svtype" in df.columns:
        mask = df["long_ont_info_svtype"].notna() & (df["long_ont_info_svtype"] != "")
        df.loc[mask, "long_ont_score"] = df.loc[mask, "long_ont_score"].fillna(0)
    if "long_pacbio_info_svtype" in df.columns:
        mask = df["long_pacbio_info_svtype"].notna() & (df["long_pacbio_info_svtype"] != "")
        df.loc[mask, "long_pacbio_score"] = df.loc[mask, "long_pacbio_score"].fillna(0)

    write_csv_tables(df, args.out)

    print(f"Wrote CSV tables to: {args.out}")


if __name__ == "__main__":
    main()
