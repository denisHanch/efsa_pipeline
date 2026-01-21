#!/usr/bin/env python3
"""
Merge SV summaries from:
  1) assembly/syri summary (required columns: chrom, start, end, svtype)
  2) long-read SV summary (required columns: chrom, start, end, svtype; optional: info_svtype, supporting_reads, score)
  3) short-read SV summary (required columns: chrom, start, end, svtype; optional: info_svtype, supporting_reads, score)

Outputs a single XLSX with 5 tabs:
  Insertions, Deletions, Replacements, Inversions, Translocations

Events are clustered by (chrom, standardized_svtype) using interval overlap with an optional tolerance (bp),
PLUS breakpoint proximity (start or end must be within tol bp). This avoids merging very large calls with
many unrelated smaller calls.

Within each event, at most one record per source is chosen (best by supporting_reads then score).

Usage:
  python merge_sv_results_v2.py --asm assembly.tsv --long long.tsv --short short.tsv --out merged.xlsx --tol 10
Any of --asm/--long/--short may be omitted; columns for missing sources are still present but empty.
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

# Common caller types -> standardized
INFO_MAP = {
    "INS": "INS",
    "DEL": "DEL",
    "INV": "INV",
    "TRA": "TRA",
    "TRANS": "TRA",
    "BND": "TRA",
    "CTX": "TRA",
    # IMPORTANT: user requested DUP events to be treated as REPLACEMENTS
    "DUP": "RPL",
    "DUP:TANDEM": "RPL",
    "DUP:INT": "RPL",
    "RPL": "RPL",
    "REPL": "RPL",
    "SUB": "RPL",
    "SNV": "RPL",
}

# Assembly/syri-ish svtype prefixes -> standardized (best-effort; adjust if needed)
ASM_PREFIX_MAP = [
    (r"^DEL", "DEL"),
    (r"^INS", "INS"),
    (r"^INV", "INV"),
    (r"^TRANS", "TRA"),
    (r"^TRA", "TRA"),
    (r"^BND", "TRA"),
    (r"^DUP", "RPL"),
    # example types seen in some assembly/syri summaries:
    (r"^CPG", "INS"),  # best-effort: "copy gain" -> insertion-like
    (r"^CPL", "DEL"),  # best-effort: "copy loss" -> deletion-like
    (r"^SYN", "RPL"),  # best-effort: treat synteny-alignment blocks as replacements
]


def _to_int(x) -> Optional[int]:
    if x is None or (isinstance(x, float) and np.isnan(x)):
        return None
    s = str(x).strip()
    if s == "" or s == "." or s.lower() == "nan":
        return None
    try:
        return int(float(s))
    except Exception:
        return None


def _to_float(x) -> Optional[float]:
    if x is None or (isinstance(x, float) and np.isnan(x)):
        return None
    s = str(x).strip()
    if s == "" or s == "." or s.lower() == "nan":
        return None
    try:
        return float(s)
    except Exception:
        return None


def standardize_type(raw_svtype: str, info_svtype: Optional[str], source: str) -> str:
    """
    Standardize to one of: INS, DEL, RPL, INV, TRA
    """
    # Prefer info_svtype when present (long/short summaries)
    if info_svtype is not None:
        key = str(info_svtype).strip().upper()
        if key in INFO_MAP:
            return INFO_MAP[key]

    raw = ("" if raw_svtype is None else str(raw_svtype)).strip()
    raw_u = raw.upper()

    # Generic parsing: if it contains something like ".DEL." etc.
    for token, mapped in INFO_MAP.items():
        if token in raw_u:
            return mapped

    if source == "asm":
        for pat, mapped in ASM_PREFIX_MAP:
            if re.search(pat, raw_u):
                return mapped

    # Fallback
    return "RPL"


# ----------------------------
# Event clustering
# ----------------------------

@dataclass
class Record:
    source: str  # asm, long, short
    chrom: str
    start: int
    end: int
    std_type: str
    raw_svtype: str
    info_svtype: Optional[str] = None
    supporting_reads: Optional[int] = None
    score: Optional[float] = None


@dataclass
class EventCluster:
    chrom: str
    std_type: str
    start: int
    end: int
    members: List[Record]


def intervals_overlap(a_start: int, a_end: int, b_start: int, b_end: int, tol: int) -> bool:
    """True if intervals likely describe the same SV event.

    Require:
      1) interval overlap (with gap tolerance), AND
      2) at least one breakpoint close (start or end within tol)
    """
    overlap = (b_start <= a_end + tol) and (b_end >= a_start - tol)
    if not overlap:
        return False
    return (abs(a_start - b_start) <= tol) or (abs(a_end - b_end) <= tol)


def cluster_records(records: List[Record], tol: int) -> List[EventCluster]:
    clusters: List[EventCluster] = []
    if not records:
        return clusters

    key_to_recs: Dict[Tuple[str, str], List[Record]] = {}
    for r in records:
        key_to_recs.setdefault((r.chrom, r.std_type), []).append(r)

    for (chrom, std_type), recs in key_to_recs.items():
        recs_sorted = sorted(recs, key=lambda x: (x.start, x.end))
        cur_start = recs_sorted[0].start
        cur_end = recs_sorted[0].end
        cur_members = [recs_sorted[0]]

        for r in recs_sorted[1:]:
            if intervals_overlap(cur_start, cur_end, r.start, r.end, tol):
                cur_end = max(cur_end, r.end)
                cur_start = min(cur_start, r.start)
                cur_members.append(r)
            else:
                clusters.append(
                    EventCluster(chrom=chrom, std_type=std_type, start=cur_start, end=cur_end, members=cur_members)
                )
                cur_start, cur_end, cur_members = r.start, r.end, [r]

        clusters.append(EventCluster(chrom=chrom, std_type=std_type, start=cur_start, end=cur_end, members=cur_members))

    clusters.sort(key=lambda c: (c.std_type, c.chrom, c.start, c.end))
    return clusters


def choose_best_record(recs: List[Record]) -> Optional[Record]:
    if not recs:
        return None

    def key(r: Record):
        sup = r.supporting_reads if r.supporting_reads is not None else -1
        sc = r.score if r.score is not None else -1.0
        length = (r.end - r.start + 1)
        return (sup, sc, -length)

    return sorted(recs, key=key, reverse=True)[0]


# ----------------------------
# IO helpers
# ----------------------------

def read_tsv(path: str) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype=str, comment="#")


def load_records(path: Optional[str], source: str) -> List[Record]:
    if not path:
        return []
    if not os.path.exists(path):
        raise FileNotFoundError(f"{source} file not found: {path}")

    df = read_tsv(path)

    required = {"chrom", "start", "end", "svtype"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"{source} file is missing required columns: {sorted(missing)}")

    out: List[Record] = []
    for _, row in df.iterrows():
        chrom = str(row["chrom"])
        start = _to_int(row["start"])
        end = _to_int(row["end"])
        if start is None or end is None:
            continue
        if start > end:
            start, end = end, start

        raw_svtype = str(row["svtype"]) if row["svtype"] is not None else ""
        info_svtype = row.get("info_svtype", None)
        std_type = standardize_type(raw_svtype, info_svtype, source)

        supporting_reads = _to_int(row.get("supporting_reads", None))
        score = _to_float(row.get("score", None))

        out.append(
            Record(
                source=source,
                chrom=chrom,
                start=start,
                end=end,
                std_type=std_type,
                raw_svtype=raw_svtype,
                info_svtype=(None if info_svtype is None else str(info_svtype)),
                supporting_reads=supporting_reads,
                score=score,
            )
        )
    return out


def build_output_table(clusters: List[EventCluster]) -> pd.DataFrame:
    rows: List[Dict[str, object]] = []

    counters: Dict[str, int] = {t: 0 for t in TAB_BY_TYPE.keys()}

    for c in clusters:
        counters.setdefault(c.std_type, 0)
        counters[c.std_type] += 1
        # IMPORTANT: no zero padding (INV_1 not INV_000001)
        event_id = f"{c.std_type}_{counters[c.std_type]}"

        by_source: Dict[str, List[Record]] = {"asm": [], "long": [], "short": []}
        for m in c.members:
            by_source.setdefault(m.source, []).append(m)

        asm = choose_best_record(by_source.get("asm", []))
        lng = choose_best_record(by_source.get("long", []))
        sht = choose_best_record(by_source.get("short", []))

        row = {
            "event_id": event_id,
            "chrom": c.chrom,
            "std_svtype": c.std_type,
            "event_start": c.start,
            "event_end": c.end,
            "event_length_bp": (c.end - c.start + 1),

            # assembly/syri
            "asm_start": (asm.start if asm else np.nan),
            "asm_end": (asm.end if asm else np.nan),
            "asm_svtype_raw": (asm.raw_svtype if asm else ""),
            "asm_length_bp": ((asm.end - asm.start + 1) if asm else np.nan),

            # long reads
            "long_start": (lng.start if lng else np.nan),
            "long_end": (lng.end if lng else np.nan),
            "long_svtype_raw": (lng.raw_svtype if lng else ""),
            "long_info_svtype": (lng.info_svtype if lng else ""),
            "long_length_bp": ((lng.end - lng.start + 1) if lng else np.nan),
            "long_score": (lng.score if lng else np.nan),
            "long_supporting_reads": (lng.supporting_reads if lng else np.nan),

            # short reads
            "short_start": (sht.start if sht else np.nan),
            "short_end": (sht.end if sht else np.nan),
            "short_svtype_raw": (sht.raw_svtype if sht else ""),
            "short_info_svtype": (sht.info_svtype if sht else ""),
            "short_length_bp": ((sht.end - sht.start + 1) if sht else np.nan),
            "short_score": (sht.score if sht else np.nan),
            "short_supporting_reads": (sht.supporting_reads if sht else np.nan),
        }
        rows.append(row)

    df = pd.DataFrame(rows)

    col_order = [
        "event_id", "chrom", "std_svtype",
        "event_start", "event_end", "event_length_bp",
        "asm_start", "asm_end", "asm_svtype_raw", "asm_length_bp",
        "long_start", "long_end", "long_svtype_raw", "long_info_svtype",
        "long_length_bp", "long_score", "long_supporting_reads",
        "short_start", "short_end", "short_svtype_raw", "short_info_svtype",
        "short_length_bp", "short_score", "short_supporting_reads",
    ]
    for c in col_order:
        if c not in df.columns:
            df[c] = np.nan
    df = df[col_order]
    return df


def write_xlsx(df: pd.DataFrame, out_path: str):
    with pd.ExcelWriter(out_path, engine="openpyxl") as writer:
        for std_type, sheet in TAB_BY_TYPE.items():
            sub = df[df["std_svtype"] == std_type].copy()
            sub.sort_values(["chrom", "event_start", "event_end"], inplace=True, kind="mergesort")
            sub.to_excel(writer, sheet_name=sheet, index=False)

        unknown = df[~df["std_svtype"].isin(list(TAB_BY_TYPE.keys()))]
        if len(unknown) > 0:
            sub = unknown.copy()
            sub.sort_values(["chrom", "event_start", "event_end"], inplace=True, kind="mergesort")
            sub.to_excel(writer, sheet_name="Other", index=False)


def main():
    p = argparse.ArgumentParser(description="Merge SV summaries into a single XLSX with tabs by SV type.")
    p.add_argument("--asm", default=None, help="Assembly/syri SV summary TSV")
    p.add_argument("--long", dest="long_reads", default=None, help="Long-read SV summary TSV")
    p.add_argument("--short", dest="short_reads", default=None, help="Short-read SV summary TSV")
    p.add_argument("--out", required=True, help="Output XLSX path")
    p.add_argument("--tol", type=int, default=10, help="Clustering tolerance in bp (default: 10)")
    args = p.parse_args()

    records: List[Record] = []
    records.extend(load_records(args.asm, "asm"))
    records.extend(load_records(args.long_reads, "long"))
    records.extend(load_records(args.short_reads, "short"))

    clusters = cluster_records(records, tol=args.tol)
    df = build_output_table(clusters)

    write_xlsx(df, args.out)

    print(f"Wrote: {args.out}")
    print(f"Total events: {len(df)}")
    for t, sheet in TAB_BY_TYPE.items():
        print(f"  {sheet}: {int((df['std_svtype']==t).sum())}")


if __name__ == "__main__":
    main()
