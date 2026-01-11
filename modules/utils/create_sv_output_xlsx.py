#!/usr/bin/env python3
"""Creates SV output file

Inputs:
  - Assembly (SyRI) TSV: chrom, start, end, svtype
  - Short-read TSV: chrom, start, end, svtype, [info_svtype], [debreak_type], supporting_reads
  - Long-read TSV: chrom, start, end, svtype, [info_svtype], [debreak_type], supporting_reads

Main functionality implemented:
  - Separate sheets: Insertions, Deletions, Replacements, Inversions, Translocations
  - Unique IDs per sheet (e.g., INS1, INS2 ...)
  - SV length computed as end - start
  - Validate INS/DEL: keep assembly events that intersect BOTH short-read and long-read events of same type
  - SV type derived from (in order): info_svtype, svtype, debreak_type; DeBreak stores type in ALT (debreak_type)
  - For long-read supporting_reads with 4 comma-separated numbers, the 3rd number is used
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd
import openpyxl


TEMPLATE_HEADERS: Dict[str, List[str]] = {
  "Insertions": [
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
    "short-reads-based estimate number of copies"
  ],
  "Deletions": [
    "SV unique ID",
    "which type of event is taking place (possibly dismissable)",
    "which fasta sequence the given SV is occuring",
    "assembly-based SV start",
    "assembly-based SV end",
    "assembly-based SV length",
    "long reads-based SV start (median)",
    "long reads-based SV start (standard deviation)",
    "long reads-based SV end (median)",
    "long reads-based SV end (standard deviation)",
    "long reads-based SV length (median)",
    "short reads-based SV start (median)",
    "short reads-based SV start (standard deviation)",
    "short reads-based SV end (median)",
    "short reads-based SV end (standard deviation)",
    "short reads-based SV length (median)",
    "assembly-based confidence score (%, e-value, p-value)",
    "long readss-based support (how many long readss supporting the SV)",
    "long readss-based confidence score (%, p-value, e-value, etc)",
    "short reads-based support (how many short reads supporting the SV)",
    "short reads-based confidence score (%, p-value, e-value, etc)",
    "assembly-based estimate number of copies",
    "long readss-based estimate number of copies",
    "short reads-based estimate number of copies"
  ],
  "Replacements": [
    "SV unique ID",
    "which type of event is taking place (possibly dismissable)",
    "which fasta sequence the given SV is occuring",
    "assembly-based start of deleted region",
    "assembly-based end of deleted region",
    "assembly-based length of deleted region",
    "assembly-based start of inserted region",
    "assembly-based end of inserted region",
    "assembly-based length of inserted region",
    "long reads-based start of deletion (median)",
    "long reads-based start of deletion (standard deviation)",
    "long reads-based deletion end (median)",
    "long reads-based deletion end (standard deviation)",
    "long reads-based length of deletion (standard deviation)",
    "long reads-based start of insertion (median)",
    "long reads-based start of insertion (standard deviation)",
    "long reads-based end of insertion (median)",
    "long reads-based end of insertion (standard deviation)",
    "long reads-based length of insertion (median)",
    "short reads-based start of deletion (median)",
    "short reads-based start of deletion (standard deviation)",
    "short reads-based deletion end (median)",
    "short reads-based deletion end (standard deviation)",
    "short reads-based length of deletion (standard deviation)",
    "short reads-based start of insertion (median)",
    "short reads-based start of insertion (standard deviation)",
    "short reads-based end of insertion (median)",
    "short reads-based end of insertion (standard deviation)",
    "short reads-based length of insertion (median)",
    "assembly-based SV support",
    "assembly-based confidence score (%, e-value, p-value)",
    "long reads-based support (how many long reads supporting the SV)",
    "long reads-based confidence score (%, p-value, e-value, etc)",
    "short reads-based support (how many short reads supporting the SV)",
    "short reads-based confidence score (%, p-value, e-value, etc)",
    "assembly-based estimate number of copies",
    "long reads-based estimate number of copies",
    "short reads-based estimate number of copies"
  ],
  "Inversions": [
    "SV unique ID",
    "which type of event is taking place (possibly dismissable)",
    "which fasta sequence the given SV is occuring",
    "assembly-based SV start",
    "assembly-based SV end",
    "assembly-based SV length",
    "long reads-based SV start (median)",
    "long reads-based SV start (standard deviation)",
    "long reads-based SV end (median)",
    "long reads-based SV end (standard deviation)",
    "long reads-based SV length (median)",
    "short reads-based SV start (median)",
    "short reads-based SV start (standard deviation)",
    "short reads-based SV end (median)",
    "short reads-based SV end (standard deviation)",
    "short reads-based SV length (median)",
    "assembly-based SV support",
    "assembly-based confidence score (%, e-value, p-value)",
    "long reads-based support (how many long reads supporting the SV)",
    "long reads-based confidence score (%, p-value, e-value, etc)",
    "short reads-based support (how many short reads supporting the SV)",
    "long reads-based confidence score (%, p-value, e-value, etc)",
    "assembly-based estimate number of copies",
    "long reads-based estimate number of copies",
    "short long reads-based estimate number of copies"
  ],
  "Translocations": [
    "SV unique ID",
    "which type of event is taking place (possibly dismissable)",
    "sequence from which the DNA fragment originates",
    "sequence to which the DNA fragment is translocated",
    "assembly-based start coordinate of origin breakpoint",
    "assembly-based end coordinate of origin breakpoint",
    "assembly-based start coordinate of destination breakpoint",
    "assembly-based end coordinate of destination breakpoint",
    "assembly-based length of translocated fragment",
    "long reads-based start coordinate of origin breakpoint (median)",
    "long reads-based start coordinate of origin breakpoint (standard deviation)",
    "long reads-based end coordinate of origin breakpoint (median)",
    "long reads-based end coordinate of origin breakpoint (standard deviation)",
    "long reads-based start coordinate of destination breakpoint (median)",
    "long reads-based start coordinate of destination breakpoint (standard deviation)",
    "long reads-based end coordinate of destination breakpoint (median)",
    "long reads-based end coordinate of destination breakpoint (standard deviation)",
    "long reads-based length of translocated fragment",
    "short reads-based start coordinate of origin breakpoint (median)",
    "short reads-based start coordinate of origin breakpoint (standard deviation)",
    "short reads-based end coordinate of origin breakpoint (median)",
    "short reads-based end coordinate of origin breakpoint (standard deviation)",
    "short reads-based start coordinate of destination breakpoint (median)",
    "short reads-based start coordinate of destination breakpoint (standard deviation)",
    "short reads-based end coordinate of destination breakpoint (median)",
    "short reads-based end coordinate of destination breakpoint (standard deviation)",
    "short reads-based length of translocated fragment",
    "assembly-based SV support",
    "assembly-based confidence score (%, e-value, p-value)",
    "long reads-based support (how many long reads supporting the SV)",
    "long reads-based confidence score (%, p-value, e-value, etc)",
    "short reads-based support (how many short reads supporting the SV)",
    "short reads-based confidence score (%, p-value, e-value, etc)",
    "assembly-based estimate number of copies",
    "long reads-based estimate number of copies",
    "short reads-based estimate number of copies"
  ]
}


def _read_tsv(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    for c in ["start", "end"]:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    return df


def _clean_svtype(raw: str) -> str:
    s = (raw or "").strip()
    if not s:
        return ""
    s = s.strip("<>")
    s = s.split(",")[0].strip()
    s = s.upper()
    if s == "BND":
        return "TRA"
    return s


def _svtype_from_row(row: pd.Series) -> str:
    for key in ["info_svtype", "svtype", "debreak_type"]:
        if key in row and str(row[key]).strip():
            return _clean_svtype(str(row[key]))
    return ""


def _parse_supporting_reads(value: str, *, long_read_special: bool) -> Optional[int]:
    s = (value or "").strip()
    if not s:
        return None

    if long_read_special and "," in s:
        parts = [p.strip() for p in s.split(",")]
        if len(parts) >= 3:
            try:
                return int(float(parts[2]))
            except Exception:
                pass

    try:
        return int(float(s))
    except Exception:
        return None


def _intervals_intersect(a_start: float, a_end: float, b_start: float, b_end: float) -> bool:
    if any(
        x is None or (isinstance(x, float) and math.isnan(x))
        for x in [a_start, a_end, b_start, b_end]
    ):
        return False
    return max(a_start, b_start) <= min(a_end, b_end)


def _overlapping_rows(df: pd.DataFrame, chrom: str, start: float, end: float, svtype: str) -> pd.DataFrame:
    if df.empty:
        return df
    sub = df[df["chrom"] == chrom]
    if "svtype_norm" in sub.columns:
        sub = sub[sub["svtype_norm"] == svtype]
    mask = sub.apply(lambda r: _intervals_intersect(start, end, r["start"], r["end"]), axis=1)
    return sub[mask]


def _median(values: List[float]) -> Optional[float]:
    vals = [v for v in values if v is not None and not (isinstance(v, float) and math.isnan(v))]
    if not vals:
        return None
    vals_sorted = sorted(vals)
    n = len(vals_sorted)
    mid = n // 2
    if n % 2 == 1:
        return float(vals_sorted[mid])
    return float((vals_sorted[mid - 1] + vals_sorted[mid]) / 2)


def build_workbook(assembly_tsv: Path, short_tsv: Path, long_tsv: Path, out_xlsx: Path) -> None:
    asm = _read_tsv(assembly_tsv)
    sr = _read_tsv(short_tsv)
    lr = _read_tsv(long_tsv)

    # Standardize column names for assembly: expect chrom/start/end/svtype
    asm_cols = {str(c).lower(): c for c in asm.columns}
    asm = asm.rename(
        columns={
            asm_cols.get("chrom", "chrom"): "chrom",
            asm_cols.get("start", "start"): "start",
            asm_cols.get("end", "end"): "end",
            asm_cols.get("svtype", "svtype"): "svtype",
        }
    )
    asm["svtype_norm"] = asm["svtype"].apply(_clean_svtype)

    def prep_reads(df: pd.DataFrame, long_read_special: bool) -> pd.DataFrame:
        out = df.copy()
        if "info_svtype" not in out.columns:
            out["info_svtype"] = ""
        if "debreak_type" not in out.columns:
            out["debreak_type"] = ""
        if "svtype" not in out.columns:
            out["svtype"] = ""
        out["svtype_norm"] = out.apply(_svtype_from_row, axis=1)

        if "supporting_reads" not in out.columns:
            out["supporting_reads"] = ""
        out["supporting_reads_norm"] = out["supporting_reads"].apply(
            lambda x: _parse_supporting_reads(str(x), long_read_special=long_read_special)
        )
        return out

    sr = prep_reads(sr, long_read_special=False)
    lr = prep_reads(lr, long_read_special=True)

    wb = openpyxl.Workbook()
    wb.remove(wb.active)

    def add_sheet(name: str):
        ws = wb.create_sheet(title=name)
        for col_idx, h in enumerate(TEMPLATE_HEADERS[name], start=1):
            ws.cell(row=1, column=col_idx, value=h)
        return ws

    sheets = {name: add_sheet(name) for name in TEMPLATE_HEADERS.keys()}

    bucket_map = {
        "Insertions": {"INS"},
        "Deletions": {"DEL"},
        "Replacements": {"REP", "REPLACEMENT"},
        "Inversions": {"INV"},
        "Translocations": {"TRA", "TRANSLOCATION"},
    }

    prefix_map = {
        "Insertions": "INS",
        "Deletions": "DEL",
        "Replacements": "REP",
        "Inversions": "INV",
        "Translocations": "TRA",
    }

    counters = {name: 0 for name in TEMPLATE_HEADERS.keys()}

    def write_row(sheet_name: str, row_dict: Dict[str, object]) -> None:
        ws = sheets[sheet_name]
        counters[sheet_name] += 1
        event_id = f"{prefix_map[sheet_name]}{counters[sheet_name]}"

        row_dict = dict(row_dict)
        row_dict["SV unique ID"] = event_id

        r = ws.max_row + 1
        for c, h in enumerate(TEMPLATE_HEADERS[sheet_name], start=1):
            ws.cell(row=r, column=c, value=row_dict.get(h, ""))

    def summarize(df: pd.DataFrame) -> Tuple[Optional[float], Optional[float], Optional[float], Optional[int]]:
        if df.empty:
            return None, None, None, None
        s_med = _median(df["start"].tolist())
        e_med = _median(df["end"].tolist())
        l_med = None
        if s_med is not None and e_med is not None:
            l_med = e_med - s_med

        supp_vals = [v for v in df.get("supporting_reads_norm", pd.Series([], dtype=object)).tolist() if v is not None]
        supp = int(max(supp_vals)) if supp_vals else None
        return s_med, e_med, l_med, supp

    for _, a in asm.iterrows():
        chrom = str(a.get("chrom", "")).strip()
        start = a.get("start", float("nan"))
        end = a.get("end", float("nan"))
        svt = str(a.get("svtype_norm", "")).strip().upper()

        sheet_name = None
        for sname, types in bucket_map.items():
            if svt in types:
                sheet_name = sname
                break
        if sheet_name is None:
            continue

        sr_ov = _overlapping_rows(sr, chrom, start, end, svt)
        lr_ov = _overlapping_rows(lr, chrom, start, end, svt)

        if svt in {"INS", "DEL"}:
            if sr_ov.empty or lr_ov.empty:
                continue

        lr_s, lr_e, lr_len, lr_supp = summarize(lr_ov)
        sr_s, sr_e, sr_len, sr_supp = summarize(sr_ov)

        asm_len = None
        if pd.notna(start) and pd.notna(end):
            asm_len = float(end) - float(start)

        base = {
            "which type of event is taking place (possibly dismissable)": svt,
            "which fasta sequence the given SV is occuring": chrom,
            "assembly-based SV start": None if pd.isna(start) else float(start),
            "assembly-based SV end": None if pd.isna(end) else float(end),
            "assembly-based SV length": asm_len,
            "long-reads-based SV start (median)": lr_s,
            "long-reads-based SV end (median)": lr_e,
            "long-reads-based SV length (median)": lr_len,
            "short-reads-based SV start (median)": sr_s,
            "short-reads-based SV end (median)": sr_e,
            "short-reads-based SV length (median)": sr_len,
            "long reads-based support (how many long reads supporting the SV)": lr_supp,
            "short reads-based support (how many short reads supporting the SV)": sr_supp,
        }

        if sheet_name == "Deletions":
            base["assembly-based deleted segment length"] = asm_len

        write_row(sheet_name, base)

    wb.save(out_xlsx)


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--assembly", required=True, type=Path, help="SyRI assembly SV summary TSV")
    p.add_argument("--short", required=True, type=Path, help="Short-read SV summary TSV")
    p.add_argument("--long", required=True, type=Path, help="Long-read SV summary TSV")
    p.add_argument("--out", required=True, type=Path, help="Output XLSX path")
    args = p.parse_args()
    build_workbook(args.assembly, args.short, args.long, args.out)


if __name__ == "__main__":
    main()
