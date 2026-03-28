#!/usr/bin/env python3
"""
Merge SV summaries from:
  1) assembly/syri summary (required columns: chrom, start, end, svtype)
  2) long-read SV summary (required columns: chrom, start, end, svtype; optional: info_svtype, supporting_reads, score) both pacbio and ont
  3) short-read SV summary (required columns: chrom, start, end, svtype; optional: info_svtype, supporting_reads, score)

Outputs a folder with single CSV for each type of variations:
  Insertions, Deletions, Duplications, Replacements, Inversions, Translocations

SHORT-READ VARIANT TYPE CONVENTIONS (from delly):
  - DEL (deletion): real interval (start < end)
  - DUP (duplication): real interval (start < end)
  - INV (inversion): real interval (start < end)
  - INS (insertion): single position (start == end, represents insertion point)
  - TRA (translocation): breakpoint only (start == end, represents breakpoint position)
  
  For INS and TRA, the variant length is captured in svlen, not from end-start calculation.
  The code correctly handles all these conventions by:
  - Preserving start and end as reported
  - Using svlen for event_length_bp calculation (not end-start)
  - Clustering by position proximity with configurable tolerance

Events are clustered by (chrom, standardized_svtype) using interval overlap with an optional tolerance (bp),
PLUS breakpoint proximity (start or end must be within tol bp). This avoids merging very large calls with
many unrelated smaller calls.

Within each event, at most one record per source is chosen (best by supporting_reads then score).

A final pass scans all final SV rows and adds `linked_event` entries for any overlapping events
on the same chromosome. This captures exact coordinate matches, partial overlaps, and nested
events for both same-type and cross-type SV rows. Set --cross_type_tol to a small positive
integer to also link near-identical final coordinates.

Usage:
    python create_sv_output.py --asm assembly.tsv --long_ont long_ont.tsv --long_pacbio long_pb.tsv --short short.tsv --out outdir --tol 10 --cross_type_tol 0

Any of the inputs may be omitted; the script will process whichever of the assembly, long_ont,
long_pacbio or short tables are provided. Long-read inputs (ONT and PacBio) are treated
the same for clustering, but preserved as separate `long_ont_*` and `long_pacbio_*` columns in the final output tables.
"""
from __future__ import annotations

import atexit
import argparse
import logging
import glob
import os
import re
from collections import Counter
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import pandas as pd
import structlog

# ----------------------------
# SV type standardization
# ----------------------------

TAB_BY_TYPE = {
    "INS": "Insertions",
    "DEL": "Deletions",
    "DUP": "Duplications",
    "RPL": "Replacements",
    "INV": "Inversions",
    "TRA": "Translocations",
}

INFO_MAP = {
    "INS": "INS",
    "DEL": "DEL",
    "DUP": "DUP",
    "DUP:TANDEM": "DUP",
    "DUP:INT": "DUP",
    "INV": "INV",
    "TRA": "TRA",
    "TRANS": "TRA",
    "BND": "TRA",
    "CTX": "TRA",
    "RPL": "RPL",
    "REPL": "RPL",
    "SUB": "RPL",
    "SNV": "RPL",
}

ASM_PREFIX_MAP = [
    (r"^DEL", "DEL"),
    (r"^INS", "INS"),
    (r"^INV", "INV"),
    (r"^DUP", "DUP"),
    (r"^TRANS", "TRA"),
    (r"^TRA", "TRA"),
    (r"^BND", "TRA"),
    (r"^CPG", "INS"),
    (r"^CPL", "DEL"),
    (r"^SYN", "RPL"),
]

DEFAULT_LOG_DIR = Path("data") / "outputs" / "logs"
ROW_LOG_FIELDS = (
    "chrom",
    "#chrom",
    "start",
    "end",
    "svtype",
    "info_svtype",
    "svlen",
    "supporting_reads",
    "score",
    "supporting_methods",
    "RDCN",
    "chr2",
    "pos2",
    "start_mod",
    "end_mod",
)


def _resolve_log_dir() -> Path:
    """Resolve a writable log directory under data/outputs/logs/."""
    env_log_dir = os.environ.get("SV_OUTPUT_LOG_DIR")
    candidates: List[Path] = []

    if env_log_dir:
        candidates.append(Path(env_log_dir))

    try:
        candidates.append(Path(__file__).resolve().parents[2] / DEFAULT_LOG_DIR)
    except Exception:
        pass

    candidates.append(Path.cwd() / DEFAULT_LOG_DIR)

    last_error = None
    seen = set()
    for candidate in candidates:
        if candidate in seen:
            continue
        seen.add(candidate)
        try:
            candidate.mkdir(parents=True, exist_ok=True)
            return candidate
        except OSError as exc:
            last_error = exc

    raise OSError(
        f"Unable to create a writable log directory under {DEFAULT_LOG_DIR}"
    ) from last_error


def _render_log_entry(_, __, event_dict):
    """Render one log line as: [timestamp] {json-without-redundant-fields}."""
    event_dict = dict(event_dict)
    timestamp = event_dict.pop("timestamp", None)
    prefix = f"[{timestamp}] " if timestamp else ""
    return prefix + structlog.processors.JSONRenderer(sort_keys=True)(None, None, event_dict)


def setup_logging():
    """Configure a file-backed structlog logger for create_sv_output.py."""
    log_dir = _resolve_log_dir()
    timestamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    log_path = log_dir / f"create_sv_output_{timestamp}_{os.getpid()}.log"

    log_file = log_path.open("w", encoding="utf-8")
    atexit.register(log_file.close)

    structlog.configure(
        processors=[
            structlog.processors.TimeStamper(fmt="iso", utc=True),
            structlog.processors.add_log_level,
            structlog.processors.StackInfoRenderer(),
            structlog.processors.format_exc_info,
            _render_log_entry,
        ],
        wrapper_class=structlog.make_filtering_bound_logger(logging.INFO),
        context_class=dict,
        logger_factory=structlog.PrintLoggerFactory(file=log_file),
        cache_logger_on_first_use=True,
    )

    logger = structlog.get_logger("create_sv_output")
    logger.info("logging_initialized")
    return logger, log_path


def _json_safe(value):
    if value is None:
        return None

    if isinstance(value, Path):
        return str(value)

    if isinstance(value, np.generic):
        return value.item()

    try:
        if pd.isna(value):
            return None
    except Exception:
        pass

    return value


def _row_snapshot(row: Dict[str, object]) -> Dict[str, object]:
    snapshot: Dict[str, object] = {}
    for key in ROW_LOG_FIELDS:
        if key not in row:
            continue
        value = _json_safe(row.get(key))
        if value in (None, ""):
            continue
        snapshot[key] = value
    return snapshot


def _clean_dict(values: Dict[str, object]) -> Dict[str, object]:
    cleaned: Dict[str, object] = {}
    for key, value in values.items():
        safe_value = _json_safe(value)
        if safe_value in (None, "", [], {}):
            continue
        cleaned[key] = safe_value
    return cleaned


def _is_missing(value) -> bool:
    if value is None:
        return True
    if isinstance(value, str) and value.strip() == "":
        return True
    try:
        return bool(pd.isna(value))
    except Exception:
        return False

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

def standardize_type_details(raw_svtype, info_svtype, source):
    if info_svtype is not None:
        key = str(info_svtype).strip().upper()
        if key in INFO_MAP:
            return INFO_MAP[key], "info_svtype", key

    raw_u = str(raw_svtype).upper()

    for token, mapped in INFO_MAP.items():
        if token in raw_u:
            return mapped, "raw_svtype", token

    if source == "asm":
        for pat, mapped in ASM_PREFIX_MAP:
            if re.search(pat, raw_u):
                return mapped, "asm_prefix", pat

    return "RPL", "fallback_default", None


def standardize_type(raw_svtype, info_svtype, source):
    return standardize_type_details(raw_svtype, info_svtype, source)[0]

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
    chr2: Optional[str] = None
    pos2: Optional[int] = None
    supporting_methods: Optional[str] = None
    svlen: Optional[int] = None
    start_mod: Optional[int] = None
    end_mod: Optional[int] = None

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
    """Check if two intervals overlap or are close enough to cluster together.
    
    Works correctly for all variant types:
    - Interval variants (DEL, DUP, INV): real intervals with start < end
    - Point variants (INS): start == end at insertion point
    - Breakpoint variants (TRA): start == end at breakpoint
    
    For point variants, the condition simplifies correctly since start == end.
    """
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
            -(r.svlen if r.svlen is not None else -1),
        )

    return max(recs, key=key)

def read_tsv(path):
    return pd.read_csv(path, sep="\t", dtype=str, comment="#")


# ----------------------------
# Supporting-read resolution
# ----------------------------

def _load_supp_reads_files(paths):
    """Load all per-caller supporting-reads TSVs into a single DataFrame.

    Expected columns per file: chrom, start, end, svtype, supporting_reads.
    Returns a DataFrame with columns: chrom, start, svtype, std_type, reads.
    """
    frames = []
    for path in paths:
        try:
            df = read_tsv(path)
        except Exception:
            continue
        if df.empty:
            continue
        df = df.rename(columns={"supporting_reads": "reads"})
        df["start"] = pd.to_numeric(df["start"], errors="coerce")
        df["reads"] = pd.to_numeric(df["reads"], errors="coerce").fillna(0).astype(int)
        df = df.dropna(subset=["chrom", "start"])
        df["start"] = df["start"].astype(int)
        # Standardise caller svtype so it matches the record's std_type
        df["std_type"] = df["svtype"].apply(
            lambda s: standardize_type(str(s).strip().upper(), None, "long")
        )
        frames.append(df[["chrom", "start", "std_type", "reads"]])

    if not frames:
        return pd.DataFrame(columns=["chrom", "start", "std_type", "reads"])
    return pd.concat(frames, ignore_index=True)


def resolve_supporting_reads(records, supp_files, tol=1000):
    """Replace missing supporting_reads on long-read Records using pandas merge_asof.

    For each long-read record, find matching caller entries by chrom + svtype
    within *tol* bp of the start position, and sum the reads across all callers.
    """
    if not supp_files:
        return records

    caller_df = _load_supp_reads_files(supp_files)
    if caller_df.empty:
        return records

    # Build a DataFrame of long-read records that need resolution
    long_indices = [
        i for i, r in enumerate(records)
        if r.source.startswith("long")
    ]
    if not long_indices:
        return records

    long_rows = pd.DataFrame([
        {"_idx": i, "chrom": records[i].chrom,
         "start": records[i].start, "std_type": records[i].std_type}
        for i in long_indices
    ])

    # Sum caller reads per (chrom, std_type, start) so each position
    # carries the total across sniffles + cuteSV + debreak
    caller_agg = (
        caller_df
        .groupby(["chrom", "std_type", "start"], as_index=False)["reads"]
        .sum()
        .sort_values("start")
    )

    # merge_asof: nearest match within tolerance, grouped by chrom + svtype
    long_rows = long_rows.sort_values("start")
    merged = pd.merge_asof(
        long_rows,
        caller_agg,
        on="start",
        by=["chrom", "std_type"],
        tolerance=tol,
        direction="nearest",
    )

    # Write resolved values back into the Record objects
    for _, row in merged.iterrows():
        reads = row.get("reads")
        if pd.notna(reads) and int(reads) > 0:
            records[int(row["_idx"])].supporting_reads = int(reads)

    return records

def load_records(path, source, logger=None):
    """Load records from a TSV file into Record objects.

    - path: path to tsv (str or None). If falsy, returns empty list.
    - source: one of 'asm', 'long_ont', 'long_pacbio', 'long' or 'short'.
      Sources starting with 'long' are treated as long-read sources and stored
      separately (e.g. 'long_ont' vs 'long_pacbio').
    
    COORDINATE CONVENTIONS:
      For short reads from delly:
      - Interval variants (DEL, DUP, INV): start and end define the real interval
      - Point variants (INS): start == end at insertion point; length is in svlen
      - Breakpoint variants (TRA): start == end at breakpoint; length is in svlen
      
      No swapping of start/end occurs for point variants since start == end.
    """
    if not path:
        if logger is not None:
            logger.info("load_records_input_missing", source=source, reason="path_not_provided")
        return []

    path_str = str(path)
    source_logger = logger.bind(source=source, input_path=path_str) if logger is not None else None

    try:
        df = read_tsv(path)
    except Exception as exc:
        if source_logger is not None:
            source_logger.error(
                "load_records_read_failed",
                error=str(exc),
            )
        return []

    out = []
    rows_seen = 0
    skipped_records = 0
    modified_records = 0
    skip_reasons: Counter = Counter()
    modification_reasons: Counter = Counter()

    if source_logger is not None:
        source_logger.info(
            "load_records_started",
            total_rows=int(len(df.index)),
        )

    # defensive: accept files that may be missing optional columns
    for row_number, row in enumerate(df.to_dict(orient="records"), start=2):
        rows_seen += 1
        row_before = _row_snapshot(row)

        chrom = row.get("chrom")
        if _is_missing(chrom):
            chrom = row.get("#chrom")
        raw_svtype = row.get("svtype")
        if _is_missing(chrom) or _is_missing(raw_svtype):
            # malformed row -> skip
            skipped_records += 1
            skip_reasons["missing_required_fields"] += 1
            if source_logger is not None:
                source_logger.warning(
                    "load_records_row_skipped",
                    row_number=row_number,
                    reason="missing_required_fields",
                    row=row_before,
                )
            continue

        chrom = str(chrom).strip()
        raw_svtype = str(raw_svtype).strip()

        start = _to_int(row.get("start"))
        end = _to_int(row.get("end"))
        if start is None or end is None:
            skipped_records += 1
            skip_reasons["invalid_coordinates"] += 1
            if source_logger is not None:
                source_logger.warning(
                    "load_records_row_skipped",
                    row_number=row_number,
                    reason="invalid_coordinates",
                    row=row_before,
                )
            continue

        row_modifications = []
        modification_details: Dict[str, object] = {}

        # Only swap if start > end; for INS/TRA with start==end, preserve as-is
        if start > end:
            original_start, original_end = start, end
            start, end = end, start
            row_modifications.append("coordinates_swapped")
            modification_reasons["coordinates_swapped"] += 1
            modification_details["coordinates_swapped"] = {
                "from": {"start": original_start, "end": original_end},
                "to": {"start": start, "end": end},
            }

        std, std_resolution, std_match = standardize_type_details(
            raw_svtype,
            row.get("info_svtype"),
            source,
        )
        if std_resolution == "fallback_default":
            row_modifications.append("svtype_defaulted_to_rpl")
            modification_reasons["svtype_defaulted_to_rpl"] += 1
            modification_details["svtype_defaulted_to_rpl"] = {
                "raw_svtype": _json_safe(raw_svtype),
                "info_svtype": _json_safe(row.get("info_svtype")),
                "std_type": std,
            }

        # treat any long* sources as long for supporting_methods
        supporting_methods = row.get("supporting_methods") if str(source).startswith("long") else None
        copy_number = _to_int(row.get("RDCN")) if source == "short" else None
        chr2 = row.get("chr2") if source == "short" else None
        pos2 = _to_int(row.get("pos2")) if source == "short" else None
        start_mod = _to_int(row.get("start_mod")) if source == "asm" else None
        end_mod = _to_int(row.get("end_mod")) if source == "asm" else None
        svlen = _to_int(row.get("svlen"))

        # Derive missing svlen from coordinates.
        if svlen is None and start is not None and end is not None:
            if std in ("DEL", "DUP", "INV", "RPL"):
                svlen = (end - start + 1) if end >= start else None
            elif std == "INS":
                svlen = None
            elif std == "TRA":
                svlen = 0
            else:
                svlen = (end - start + 1) if end >= start else None

            if svlen is not None:
                row_modifications.append("svlen_derived")
                modification_reasons["svlen_derived"] += 1
                modification_details["svlen_derived"] = {
                    "std_type": std,
                    "derived_svlen": svlen,
                    "start": start,
                    "end": end,
                }

        if row_modifications:
            modified_records += 1
            if source_logger is not None:
                row_after = dict(row_before)
                row_after.update(
                    _clean_dict(
                        {
                            "chrom": chrom,
                            "start": start,
                            "end": end,
                            "svtype": raw_svtype,
                            "info_svtype": row.get("info_svtype"),
                            "svlen": svlen,
                            "std_type": std,
                        }
                    )
                )
                source_logger.info(
                    "load_records_row_modified",
                    row_number=row_number,
                    modifications=row_modifications,
                    details=_clean_dict(modification_details),
                    row_before=row_before,
                    row_after=row_after,
                    std_type=std,
                    std_type_resolution=std_resolution,
                    std_type_match=std_match,
                )

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
                chr2,
                pos2,
                supporting_methods,
                svlen,
                start_mod,
                end_mod,
            )
        )

    if source_logger is not None:
        source_logger.info(
            "load_records_completed",
            total_rows=rows_seen,
            loaded_records=len(out),
            skipped_records=skipped_records,
            modified_records=modified_records,
            skip_reasons=dict(skip_reasons),
            modification_reasons=dict(modification_reasons),
        )

    return out

def build_output_table(clusters):
    rows = []
    counters = {k: 0 for k in TAB_BY_TYPE}

    for c in clusters:
        counters[c.std_type] += 1
        eid = f"{c.std_type}_{counters[c.std_type]}"

        by_src = {"asm": [], "long_ont": [], "long_pacbio": [], "short": []}
        for m in c.members:
            by_src.setdefault(m.source, []).append(m)

        asm = choose_best_record(by_src.get("asm", []))
        long_ont = choose_best_record(by_src.get("long_ont", []))
        long_pacbio = choose_best_record(by_src.get("long_pacbio", []))
        sht = choose_best_record(by_src.get("short", []))
        pipelines_confirmed = sum(x is not None for x in (asm, long_ont, long_pacbio, sht))

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

        # Calculate lengths for each source.
        # Restore last week's insertion handling: INS lengths come only from reported svlen values.
        asm_length = asm.svlen if (asm and pd.notna(asm.svlen)) else np.nan
        long_ont_length = long_ont.svlen if (long_ont and pd.notna(long_ont.svlen)) else np.nan
        long_pacbio_length = long_pacbio.svlen if (long_pacbio and pd.notna(long_pacbio.svlen)) else np.nan
        sht_length = sht.svlen if (sht and pd.notna(sht.svlen)) else np.nan

        # Calculate event coordinates.
        # Restore last week's insertion handling: INS uses cluster bounds rather than a single consensus point.
        member_starts = [m.start for m in c.members]
        member_ends = [m.end for m in c.members]
        if c.std_type == "INS":
            overlap_start = c.start
            overlap_end = c.end
        elif c.std_type == "TRA":
            consensus_pos = int(round(float(np.median(member_starts)))) if member_starts else c.start
            overlap_start = consensus_pos
            overlap_end = consensus_pos
        else:
            overlap_start = max(member_starts) if member_starts else c.start
            overlap_end = min(member_ends) if member_ends else c.end
        
        lengths = [asm_length, long_ont_length, long_pacbio_length, sht_length]
        valid_lengths = [l for l in lengths if pd.notna(l)]
        if c.std_type == "INS":
            event_length_bp = min(valid_lengths) if valid_lengths else np.nan
        else:
            # event_length_bp uses svlen (from variant callers) rather than coordinate difference
            # This correctly handles TRA where start==end but svlen still reports or defaults the size
            event_length_bp = min(valid_lengths) if valid_lengths else np.nan
        row = {
            "event_id": eid,
            "chrom": c.chrom,
            "std_svtype": c.std_type,
            "event_start": overlap_start,
            "event_end": overlap_end,
            "event_length_bp": event_length_bp,

            "asm_start": asm.start if asm else np.nan,
            "asm_end": asm.end if asm else np.nan,
            "asm_start_mod": asm.start_mod if asm else np.nan,
            "asm_end_mod": asm.end_mod if asm else np.nan,
            "asm_length": asm_length,
            "asm_svtype_raw": asm.raw_svtype if asm else "",
            "asm_score": 1,

            "long_ont_start": long_ont.start if long_ont else np.nan,
            "long_ont_end": long_ont.end if long_ont else np.nan,
            "long_ont_length": long_ont_length,
            "long_ont_info_svtype": long_ont.info_svtype if long_ont else "",
            "long_ont_score": long_ont.score if long_ont else np.nan,
            "long_ont_supporting_reads": long_ont.supporting_reads if long_ont else np.nan,
            "long_ont_supporting_methods": long_ont.supporting_methods if long_ont else np.nan,

            "long_pacbio_start": long_pacbio.start if long_pacbio else np.nan,
            "long_pacbio_end": long_pacbio.end if long_pacbio else np.nan,
            "long_pacbio_length": long_pacbio_length,
            "long_pacbio_info_svtype": long_pacbio.info_svtype if long_pacbio else "",
            "long_pacbio_score": long_pacbio.score if long_pacbio else np.nan,
            "long_pacbio_supporting_reads": long_pacbio.supporting_reads if long_pacbio else np.nan,
            "long_pacbio_supporting_methods": long_pacbio.supporting_methods if long_pacbio else np.nan,

            "short_start": sht.start if sht else np.nan,
            "short_end": sht.end if sht else np.nan,
            "short_length": sht_length,
            "short_svtype_raw": sht.raw_svtype if sht else "",
            "short_info_svtype": sht.info_svtype if sht else "",
            "short_chr2": sht.chr2 if sht else np.nan,
            "short_pos2": sht.pos2 if sht else np.nan,
            "short_score": sht.score if sht else np.nan,
            "short_supporting_reads": sht.supporting_reads if sht else np.nan,
            "short_reads_copy_number_estimate": (sht.copy_number if sht else np.nan),

            "percentage_overlap": pct_str,
            "support_score": pipelines_confirmed,
        }

        rows.append(row)

    return pd.DataFrame(rows)

def _format_linked_event(event_id, std_svtype, chrom, start, end, relation):
    return f"{event_id} ({std_svtype}, {chrom}:{start}-{end}, {relation})"

def _interval_relation(left_start, left_end, right_start, right_end):
    if left_start == right_start and left_end == right_end:
        return "exact_coordinates"
    if left_start >= right_start and left_end <= right_end:
        return "nested_in"
    if left_start <= right_start and left_end >= right_end:
        return "contains"
    return "overlap"

def _linked_event_relation(left_start, left_end, right_start, right_end, coord_tol=0):
    if left_start == right_start and left_end == right_end:
        return "exact_coordinates"

    overlaps = (right_start <= left_end) and (right_end >= left_start)
    if overlaps:
        return _interval_relation(left_start, left_end, right_start, right_end)

    coord_tol = max(int(coord_tol or 0), 0)
    if coord_tol > 0 and abs(left_start - right_start) <= coord_tol and abs(left_end - right_end) <= coord_tol:
        return f"same_coordinates_within_{coord_tol}bp"

    return None

def annotate_linked_events(df, coord_tol=0):
    """
    Final post-processing step that links overlapping final SV rows using only
    the event coordinates written to the per-type tables.

    Rules:
      - same chromosome
      - intervals overlap in any way, including exact matches and nesting
      - links are added for both same-type and cross-type event pairs
      - when coord_tol > 0, near-identical coordinates can also be linked

    The resulting linked_event column contains a semicolon-separated list of
    linked event IDs with their type, coordinates, and interval relation from
    the current row's point of view.
    """
    df = df.copy()

    if df.empty:
        df["linked_event"] = pd.Series(dtype=str)
        return df

    coord_tol = max(int(coord_tol or 0), 0)
    sort_cols = ["chrom", "event_start", "event_end", "std_svtype", "event_id"]
    work = df.sort_values(sort_cols).copy()
    link_map = {idx: [] for idx in work.index}

    for chrom, group in work.groupby("chrom", sort=False):
        rows = list(group.itertuples())
        active = []

        for row in rows:
            row_start = int(row.event_start)
            row_end = int(row.event_end)

            active = [other for other in active if int(other.event_end) >= row_start - coord_tol]

            for other in active:
                other_start = int(other.event_start)
                other_end = int(other.event_end)

                if other_start > row_end + coord_tol or row_start > other_end + coord_tol:
                    continue

                row_relation = _linked_event_relation(
                    row_start,
                    row_end,
                    other_start,
                    other_end,
                    coord_tol,
                )
                if row_relation is None:
                    continue

                other_relation = _linked_event_relation(
                    other_start,
                    other_end,
                    row_start,
                    row_end,
                    coord_tol,
                )

                link_map[row.Index].append(
                    _format_linked_event(
                        other.event_id,
                        other.std_svtype,
                        chrom,
                        other_start,
                        other_end,
                        row_relation,
                    )
                )
                link_map[other.Index].append(
                    _format_linked_event(
                        row.event_id,
                        row.std_svtype,
                        chrom,
                        row_start,
                        row_end,
                        other_relation,
                    )
                )

            active.append(row)

    def _dedupe_keep_order(items):
        seen = set()
        out = []
        for item in items:
            if item not in seen:
                seen.add(item)
                out.append(item)
        return out

    df["linked_event"] = [
        "; ".join(_dedupe_keep_order(link_map[idx]))
        for idx in df.index
    ]
    return df

def write_csv_tables(df, outdir):
    os.makedirs(outdir, exist_ok=True)

    for std_type, name in TAB_BY_TYPE.items():
        sub = df[df["std_svtype"] == std_type].copy()
        if not sub.empty:
            # Drop type-specific columns
            cols_to_drop = ["std_svtype"]
            if std_type != "TRA":
                cols_to_drop.extend(["asm_start_mod", "asm_end_mod"])
            sub = sub.drop(columns=cols_to_drop, errors="ignore")
            path = os.path.join(outdir, f"{name}.csv")
            sub.sort_values(["chrom", "event_start", "event_end"]).to_csv(path, index=False)

    other = df[~df["std_svtype"].isin(TAB_BY_TYPE)].copy()
    if not other.empty:
        other = other.drop(columns=["std_svtype", "asm_start_mod", "asm_end_mod"], errors="ignore")
        other.to_csv(os.path.join(outdir, "Other.csv"), index=False)

# Main

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--asm")
    p.add_argument("--long_ont", dest="long_ont")
    p.add_argument("--long_pacbio", dest="long_pacbio")
    p.add_argument("--short", dest="short_reads")
    p.add_argument("--out", required=True, help="Output directory")
    p.add_argument("--tol", type=int, default=10, help="Within-type clustering tolerance in bp")
    p.add_argument(
        "--cross_type_tol",
        type=int,
        default=0,
        help="Tolerance in bp for linking near-identical final event coordinates in linked_event. Default 0 keeps overlap-only linking.",
    )
    args = p.parse_args()

    logger, log_path = setup_logging()
    logger.info(
        "create_sv_output_started",
        output_dir=str(args.out),
        tol=args.tol,
        cross_type_tol=args.cross_type_tol,
        inputs=_clean_dict(
            {
                "asm": args.asm,
                "long_reads": getattr(args, "long_reads", None),
                "long_ont": getattr(args, "long_ont", None),
                "long_pacbio": getattr(args, "long_pacbio", None),
                "short": args.short_reads,
            }
        ),
    )

    records = []
    records += load_records(args.asm, "asm", logger=logger)

    # long inputs: accept both ONT and PacBio; also support legacy --long
    records += load_records(getattr(args, "long_reads", None), "long", logger=logger)
    records += load_records(getattr(args, "long_ont", None), "long_ont", logger=logger)
    records += load_records(getattr(args, "long_pacbio", None), "long_pacbio", logger=logger)

    records += load_records(args.short_reads, "short", logger=logger)

    # Resolve long-read supporting reads from per-caller TSVs staged in the work directory
    supp_files = sorted(glob.glob("*_supporting_reads.tsv"))
    if supp_files:
        logger.info(
            "resolve_supporting_reads_started",
            supporting_reads_files=supp_files,
            tolerance=args.tol,
        )
        records = resolve_supporting_reads(records, supp_files, tol=args.tol)
        logger.info("resolve_supporting_reads_completed", total_records=len(records))

    if not records:
        os.makedirs(args.out, exist_ok=True)
        logger.warning(
            "create_sv_output_no_valid_records",
            output_dir=str(args.out),
        )
        print("No valid input records found; created output directory and exiting.")
        print(f"Load-record audit log: {log_path}")
        return

    clusters = cluster_records(records, args.tol)
    logger.info(
        "cluster_records_completed",
        total_records=len(records),
        total_clusters=len(clusters),
    )

    if not clusters:
        df = pd.DataFrame()
    else:
        df = build_output_table(clusters)

    if "long_ont_info_svtype" in df.columns:
        mask = df["long_ont_info_svtype"].notna() & (df["long_ont_info_svtype"] != "")
        df.loc[mask, "long_ont_score"] = df.loc[mask, "long_ont_score"].fillna(0)
    if "long_pacbio_info_svtype" in df.columns:
        mask = df["long_pacbio_info_svtype"].notna() & (df["long_pacbio_info_svtype"] != "")
        df.loc[mask, "long_pacbio_score"] = df.loc[mask, "long_pacbio_score"].fillna(0)

    write_csv_tables(df, args.out)

    df = annotate_linked_events(df, args.cross_type_tol)
    write_csv_tables(df, args.out)

    logger.info(
        "create_sv_output_completed",
        output_dir=str(args.out),
        total_records=len(records),
        total_clusters=len(clusters),
        total_output_rows=int(len(df.index)),
    )

    print(f"Wrote CSV tables to: {args.out}")
    print(f"Load-record audit log: {log_path}")

if __name__ == "__main__":
    main()