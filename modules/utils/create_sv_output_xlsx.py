#!/usr/bin/env python3

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd


SHEETS_IN_ORDER = ["Insertions", "Deletions", "Replacements", "Inversions", "Translocations"]

# Curated column headers (row 1 in Excel)
DISPLAY_HEADERS: List[str] = [
    "SV unique ID",
    "which type of event is taking place (possibly dismissable)",
    "which fasta sequence the given SV is occuring",
    "assembly-based SV start",
    "assembly-based SV end",
    "assembly-based SV length",
    "long-reads-based SV start (median)",
    "long-reads-based SV start (standard deviation)",
    "long-reads-based SV end (median)",
    "long-reads-based SV end (standard deviation)",
    "long-reads-based SV length (median)",
    "short-reads-based SV start (median)",
    "short-reads-based SV start (standard deviation)",
    "short-reads-based SV end (median)",
    "short-reads-based SV end (standard deviation)",
    "short-reads-based SV length (median)",
    "long-reads-based support (how many long-reads supporting the SV)",
    "long-reads-based confidence score (%, p-value, e-value, etc)",
    "short-reads-based support (how many short-reads supporting the SV)",
    "short-reads-based confidence score (%, p-value, e-value, etc)",
    "assembly-based estimate number of copies",
    "long-reads-based estimate number of copies",
    "short-reads-based estimate number of copies",
]

# Internal field names row (row 2 in Excel)
FIELD_ROW: List[str] = [
    "event_id",
    "event_type",
    "sequence_id",
    "asm_start",
    "asm_end",
    "asm_length_bp",
    "long_reads_start_median",
    "long_reads_start_std",
    "long_reads_end_median",
    "long_reads_end_std",
    "long_reads_length_median",
    "short_reads_start_median",
    "short_reads_start_std",
    "short_reads_end_median",
    "short_reads_end_std",
    "short_reads_length_median",
    "long_reads_support",
    "long_reads_score",
    "short_reads_support",
    "short_reads_score",
    "asm_copy_number_estimate",
    "long_reads_copy_number_estimate",
    "short_reads_copy_number_estimate",
]

# Mapping from standardized SVTYPE (info_svtype) to sheet + event_id prefix + event_type label
SVTYPE_MAP: Dict[str, Tuple[str, str, str]] = {
    "INS": ("Insertions", "INS", "insertion"),
    "DEL": ("Deletions", "DEL", "deletion"),
    "DUP": ("Replacements", "REP", "replacement"),  # best-fit mapping for this workflow
    "INV": ("Inversions", "INV", "inversion"),
    "TRA": ("Translocations", "TRA", "translocation"),
    "BND": ("Translocations", "TRA", "translocation"),  # common breakend label
}


def _to_numeric(series: pd.Series) -> pd.Series:
    """Convert to numeric where possible; non-numeric become NaN."""
    return pd.to_numeric(series, errors="coerce")


def _median(series: pd.Series) -> Optional[float]:
    s = _to_numeric(series).dropna()
    if s.empty:
        return None
    return float(np.median(s.to_numpy()))


def _std(series: pd.Series) -> Optional[float]:
    s = _to_numeric(series).dropna()
    if s.empty:
        return None
    # population std (ddof=0) so n=1 yields 0.0 rather than NaN
    return float(np.std(s.to_numpy(), ddof=0))


def _lengths(starts: pd.Series, ends: pd.Series) -> pd.Series:
    return _to_numeric(ends) - _to_numeric(starts) + 1


def _safe_int_sum(series: pd.Series) -> Optional[int]:
    s = _to_numeric(series).dropna()
    if s.empty:
        return None
    return int(s.sum())


def _score_aggregate(series: pd.Series) -> Optional[float]:
    """
    Aggregate score column.
    - If any numeric values exist -> median of numeric values
    - Else -> None
    """
    s = _to_numeric(series).dropna()
    if s.empty:
        return None
    return float(np.median(s.to_numpy()))


def _overlap_mask(asm_start: int, asm_end: int, starts: pd.Series, ends: pd.Series) -> pd.Series:
    s = _to_numeric(starts)
    e = _to_numeric(ends)
    return (s <= asm_end) & (e >= asm_start)


def read_inputs(assembly_path: Path, long_path: Path, short_path: Path) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    asm = pd.read_csv(assembly_path, sep="\t", dtype={"chrom": str})
    long_df = pd.read_csv(long_path, sep="\t", dtype={"chrom": str})
    short_df = pd.read_csv(short_path, sep="\t", dtype={"chrom": str})

    required_asm = {"chrom", "start", "end", "svtype"}
    required_reads = {"chrom", "start", "end", "info_svtype", "supporting_reads", "score"}

    missing_asm = required_asm - set(asm.columns)
    if missing_asm:
        raise ValueError(f"Assembly file is missing required columns: {sorted(missing_asm)}")

    missing_long = required_reads - set(long_df.columns)
    if missing_long:
        raise ValueError(f"Long-read file is missing required columns: {sorted(missing_long)}")

    missing_short = required_reads - set(short_df.columns)
    if missing_short:
        raise ValueError(f"Short-read file is missing required columns: {sorted(missing_short)}")

    # Normalize svtype strings
    long_df["info_svtype"] = long_df["info_svtype"].astype(str).str.upper()
    short_df["info_svtype"] = short_df["info_svtype"].astype(str).str.upper()

    return asm, long_df, short_df


def build_event_rows_for_interval(
    asm_row: pd.Series,
    long_df: pd.DataFrame,
    short_df: pd.DataFrame,
    svtype: str,
    event_id: str,
    event_type_label: str,
) -> Dict[str, object]:
    asm_start = int(asm_row["start"])
    asm_end = int(asm_row["end"])
    chrom = str(asm_row["chrom"])

    long_sub = long_df[(long_df["chrom"].astype(str) == chrom) & (long_df["info_svtype"] == svtype)]
    short_sub = short_df[(short_df["chrom"].astype(str) == chrom) & (short_df["info_svtype"] == svtype)]

    long_ov = long_sub[_overlap_mask(asm_start, asm_end, long_sub["start"], long_sub["end"])]
    short_ov = short_sub[_overlap_mask(asm_start, asm_end, short_sub["start"], short_sub["end"])]

    # Stats for long reads
    long_start_med = _median(long_ov["start"]) if not long_ov.empty else None
    long_start_std = _std(long_ov["start"]) if not long_ov.empty else None
    long_end_med = _median(long_ov["end"]) if not long_ov.empty else None
    long_end_std = _std(long_ov["end"]) if not long_ov.empty else None
    long_len_med = _median(_lengths(long_ov["start"], long_ov["end"])) if not long_ov.empty else None
    long_support = _safe_int_sum(long_ov["supporting_reads"]) if not long_ov.empty else None
    long_score = _score_aggregate(long_ov["score"]) if not long_ov.empty else None

    # Stats for short reads
    short_start_med = _median(short_ov["start"]) if not short_ov.empty else None
    short_start_std = _std(short_ov["start"]) if not short_ov.empty else None
    short_end_med = _median(short_ov["end"]) if not short_ov.empty else None
    short_end_std = _std(short_ov["end"]) if not short_ov.empty else None
    short_len_med = _median(_lengths(short_ov["start"], short_ov["end"])) if not short_ov.empty else None
    short_support = _safe_int_sum(short_ov["supporting_reads"]) if not short_ov.empty else None
    short_score = _score_aggregate(short_ov["score"]) if not short_ov.empty else None

    asm_length_bp = asm_end - asm_start + 1

    return {
        "event_id": event_id,
        "event_type": event_type_label,
        "sequence_id": chrom,
        "asm_start": asm_start,
        "asm_end": asm_end,
        "asm_length_bp": asm_length_bp,
        "long_reads_start_median": long_start_med,
        "long_reads_start_std": long_start_std,
        "long_reads_end_median": long_end_med,
        "long_reads_end_std": long_end_std,
        "long_reads_length_median": long_len_med,
        "short_reads_start_median": short_start_med,
        "short_reads_start_std": short_start_std,
        "short_reads_end_median": short_end_med,
        "short_reads_end_std": short_end_std,
        "short_reads_length_median": short_len_med,
        "long_reads_support": long_support,
        "long_reads_score": long_score,
        "short_reads_support": short_support,
        "short_reads_score": short_score,
        "asm_copy_number_estimate": None,
        "long_reads_copy_number_estimate": None,
        "short_reads_copy_number_estimate": None,
    }


def merge_to_workbook(asm: pd.DataFrame, long_df: pd.DataFrame, short_df: pd.DataFrame) -> Dict[str, pd.DataFrame]:
    """
    Returns dict: sheet_name -> dataframe ready to write,
    including the internal field-name row as first record (after headers).
    """
    # Prepare empty containers per sheet
    sheet_rows: Dict[str, List[Dict[str, object]]] = {s: [] for s in SHEETS_IN_ORDER}

    # Counters for event IDs per SVTYPE_MAP entry/prefix
    counters: Dict[str, int] = {}

    # Determine which svtypes appear in reads (restrict to recognized)
    observed_svtypes = sorted(set(long_df["info_svtype"]).union(set(short_df["info_svtype"])))
    observed_svtypes = [t for t in observed_svtypes if t in SVTYPE_MAP]

    for _, asm_row in asm.iterrows():
        asm_start = int(asm_row["start"])
        asm_end = int(asm_row["end"])
        chrom = str(asm_row["chrom"])

        # For each SV type, check if there is any overlapping call in either long or short
        for svt in observed_svtypes:
            # Filter quickly by chrom+svt
            long_sub = long_df[(long_df["chrom"].astype(str) == chrom) & (long_df["info_svtype"] == svt)]
            short_sub = short_df[(short_df["chrom"].astype(str) == chrom) & (short_df["info_svtype"] == svt)]

            has_long = False
            has_short = False
            if not long_sub.empty:
                has_long = bool(_overlap_mask(asm_start, asm_end, long_sub["start"], long_sub["end"]).any())
            if not short_sub.empty:
                has_short = bool(_overlap_mask(asm_start, asm_end, short_sub["start"], short_sub["end"]).any())

            if not (has_long or has_short):
                continue

            sheet_name, prefix, event_type_label = SVTYPE_MAP[svt]
            counters.setdefault(prefix, 0)
            counters[prefix] += 1
            event_id = f"{prefix}_{counters[prefix]}"

            row = build_event_rows_for_interval(
                asm_row=asm_row,
                long_df=long_df,
                short_df=short_df,
                svtype=svt,
                event_id=event_id,
                event_type_label=event_type_label,
            )
            sheet_rows[sheet_name].append(row)

    # Convert to DataFrames and add internal field-name row at the top
    out: Dict[str, pd.DataFrame] = {}
    for sheet_name in SHEETS_IN_ORDER:
        df = pd.DataFrame(sheet_rows[sheet_name], columns=FIELD_ROW)
        # Prepend the internal field-name row (matches the example GMO_output format)
        field_df = pd.DataFrame([dict(zip(FIELD_ROW, FIELD_ROW))], columns=FIELD_ROW)
        df2 = pd.concat([field_df, df], ignore_index=True)
        out[sheet_name] = df2

    return out


def write_excel(sheets: Dict[str, pd.DataFrame], output_path: Path) -> None:
    with pd.ExcelWriter(output_path, engine="openpyxl") as writer:
        for sheet_name in SHEETS_IN_ORDER:
            df = sheets[sheet_name]
            # Write with display headers as the header row
            df.to_excel(writer, sheet_name=sheet_name, index=False, header=DISPLAY_HEADERS)


def main() -> None:
    p = argparse.ArgumentParser(description="Merge assembly + long/short SV summaries into GMO-style Excel output.")
    p.add_argument("--assembly", required=True, type=Path, help="assembly_ref_x_modsyri_sv_summary.tsv")
    p.add_argument("--long", required=True, type=Path, help="long_read_sv_summary.tsv")
    p.add_argument("--short", required=True, type=Path, help="short_read_sv_summary.tsv")
    p.add_argument("--output", required=True, type=Path, help="output .xlsx path")
    args = p.parse_args()

    asm, long_df, short_df = read_inputs(args.assembly, args.long, args.short)
    sheets = merge_to_workbook(asm, long_df, short_df)
    write_excel(sheets, args.output)
    print(f"Wrote: {args.output}")


if __name__ == "__main__":
    main()
