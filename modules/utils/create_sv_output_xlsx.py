#!/usr/bin/env python3
"""
Create SV output workbook split across tabs by SV type, combining long-read and short-read calls.

High-level logic (per business requirements):
- Read long.tsv and short.tsv produced by modules/sv_calling.nf (vcf_to_table).
- Split calls into tabs: Insertions, Deletions, Replacements, Inversions, Translocations.
- Assign unique event_id per tab based on LONG-read events ordered by (chrom, start).
- If a short-read event overlaps ("is within") a long-read event of the same SV type on the same chromosome,
  it is treated as evidence for the same event (same row / same event_id).
- Optionally, use syri.tsv (assembly) to fill asm_* coordinates by best-overlap with the long-read event.
- Write results into an Excel workbook. If --template is provided, keep the first two header rows intact and
  overwrite rows 3+ only.

Notes on upstream TSV formats:
- TSV columns: chrom, start, end, svtype, debreak_type, supporting_reads
- SV type normalization:
  - prefers a recognized SV token parsed from svtype; if not recognized, falls back to debreak_type (e.g. <INS>)
  - DUP is mapped to INS (template has no DUP tab)
  - BND is mapped to TRA
  - SyRI codes: TRANS* -> TRA, CPL* -> REP
- supporting_reads in long-read summary may be "a,b,c,d"; the 3rd number (index 2) is used.

Dependencies: openpyxl only (no pandas required).
"""
from __future__ import annotations

import argparse
import csv
import math
import re
import statistics
from collections import defaultdict
from typing import Any, Dict, List, Optional, Tuple

from openpyxl import load_workbook, Workbook


KNOWN = {"INS", "DEL", "INV", "TRA", "REP"}


def to_int(x: Any) -> Optional[int]:
    try:
        if x is None:
            return None
        s = str(x).strip()
        if s == "" or s.lower() == "nan":
            return None
        return int(float(s))
    except Exception:
        return None


def parse_support(x: Any) -> Optional[int]:
    if x is None:
        return None
    s = str(x).strip()
    if s == "" or s.lower() == "nan":
        return None
    if "," in s:
        parts = [p for p in s.split(",") if p != ""]
        if len(parts) >= 3:
            return to_int(parts[2])
        return to_int(parts[-1]) if parts else None
    return to_int(s)


def normalize_type(svtype: Any, debreak_type: Any) -> Optional[str]:
    def norm_one(v: Any) -> Optional[str]:
        if v is None:
            return None
        t = str(v).strip()
        if t == "" or t.lower() == "nan":
            return None

        # <INS>
        m = re.match(r"^<(\w+)>$", t)
        if m:
            tok = m.group(1).upper()
            if tok == "BND":
                tok = "TRA"
            if tok == "DUP":
                tok = "INS"
            return tok

        # IDs like INV0000001, DEL123, TRA...
        m = re.match(r"^(INS|DEL|INV|TRA|REP|DUP|BND)", t.upper())
        if m:
            tok = m.group(1)
            if tok == "BND":
                tok = "TRA"
            if tok == "DUP":
                tok = "INS"
            return tok

        # caller like cuteSV.INS.0
        m = re.search(r"\b(INS|DEL|INV|TRA|REP|DUP|BND)\b", t.upper())
        if m:
            tok = m.group(1)
            if tok == "BND":
                tok = "TRA"
            if tok == "DUP":
                tok = "INS"
            return tok

        # SyRI codes
        up = t.upper()
        if up.startswith("CPL"):
            return "REP"
        if up.startswith("TRANS"):
            return "TRA"
        if up.startswith("INV"):
            return "INV"
        if up.startswith("DEL"):
            return "DEL"
        if up.startswith("INS"):
            return "INS"

        return None

    a = norm_one(svtype)
    b = norm_one(debreak_type)

    if a in KNOWN:
        return a
    if b in KNOWN:
        return b
    return a or b


def calc_length(start: Any, end: Any, svtype: Any, debreak_type: Any) -> Optional[int]:
    s = to_int(start)
    e = to_int(end)
    if s is None or e is None:
        return None

    t = normalize_type(svtype, debreak_type)
    if t == "INS":
        # If debreak_type contains the inserted sequence (not <INS>), use its length
        if debreak_type and not re.match(r"^<\w+>$", str(debreak_type).strip()):
            seq = str(debreak_type).strip()
            if re.fullmatch(r"[ACGTNacgtn]+", seq):
                return len(seq)
        # Otherwise coordinates don't encode insertion length
        return 0

    return e - s


def interval_overlap(a_start: int, a_end: int, b_start: int, b_end: int) -> int:
    return max(0, min(a_end, b_end) - max(a_start, b_start) + 1)


def overlaps(short_row: Dict[str, Any], long_row: Dict[str, Any]) -> bool:
    if short_row.get("chrom") != long_row.get("chrom"):
        return False
    ss, se = short_row.get("start"), short_row.get("end")
    ls, le = long_row.get("start"), long_row.get("end")
    if None in (ss, se, ls, le):
        return False
    return interval_overlap(ss, se, ls, le) > 0


def stats(vals: List[Optional[int]]) -> Tuple[Optional[float], Optional[float]]:
    vv = [v for v in vals if v is not None]
    if not vv:
        return None, None
    if len(vv) == 1:
        return float(vv[0]), 0.0
    med = float(statistics.median(vv))
    mean = sum(vv) / len(vv)
    var = sum((v - mean) ** 2 for v in vv) / (len(vv) - 1)
    return med, math.sqrt(var)


def read_tsv(path: str) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    with open(path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for r in reader:
            rows.append(dict(r))
    return rows


def enrich(rows: List[Dict[str, Any]], source: str) -> List[Dict[str, Any]]:
    out: List[Dict[str, Any]] = []
    for r in rows:
        rr = dict(r)
        rr["start"] = to_int(rr.get("start"))
        rr["end"] = to_int(rr.get("end"))
        rr["support"] = parse_support(rr.get("supporting_reads"))
        rr["type"] = normalize_type(rr.get("svtype"), rr.get("debreak_type"))
        rr["length"] = calc_length(rr.get("start"), rr.get("end"), rr.get("svtype"), rr.get("debreak_type"))
        rr["source"] = source
        out.append(rr)
    return out


def enrich_syri(rows: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    out: List[Dict[str, Any]] = []
    for r in rows:
        rr = dict(r)
        rr["start"] = to_int(rr.get("start"))
        rr["end"] = to_int(rr.get("end"))
        rr["type"] = normalize_type(rr.get("svtype"), None)
        if rr["start"] is not None and rr["end"] is not None:
            rr["length"] = rr["end"] - rr["start"]
        else:
            rr["length"] = None
        out.append(rr)
    return out


def best_asm_match(
    longr: Dict[str, Any],
    svtype: str,
    syri_by_type: Dict[str, List[Dict[str, Any]]],
) -> Optional[Dict[str, Any]]:
    candidates = syri_by_type.get(svtype, [])
    ls, le = longr.get("start"), longr.get("end")
    if ls is None or le is None:
        return None
    best = None
    bestov = 0
    for c in candidates:
        if c.get("chrom") != longr.get("chrom"):
            continue
        cs, ce = c.get("start"), c.get("end")
        if cs is None or ce is None:
            continue
        ov = interval_overlap(ls, le, cs, ce)
        if ov > bestov:
            bestov = ov
            best = c
    return best


def load_template_or_create(template_path: Optional[str]) -> Workbook:
    if template_path:
        return load_workbook(template_path)

    # No template: create workbook with expected sheets + 2 header rows (minimal)
    wb = Workbook()
    wb.remove(wb.active)

    # Minimal headers (row 2 keys; row 1 descriptions are left blank)
    headers = {
        "Insertions": [
            "event_id","event_type","sequence_id","asm_start","asm_end","asm_length_bp",
            "long_reads_start_median","long_reads_start_std","long_reads_end_median","long_reads_end_std","long_reads_length_median",
            "short_reads_start_median","short_reads_start_std","short_reads_end_median","short_reads_end_std","short_reads_length_median",
            "long_reads_support","long_reads_score","short_reads_support","short_reads_score",
            "asm_copy_number_estimate","long_reads_copy_number_estimate","short_reads_copy_number_estimate"
        ],
        "Deletions": [
            "event_id","event_type","sequence_id","asm_start","asm_end","asm_length_bp",
            "long_reads_start_median","long_reads_start_std","long_reads_end_median","long_reads_end_std","long_reads_length_median",
            "short_reads_start_median","short_reads_start_std","short_reads_end_median","short_reads_end_std","short_reads_length_median",
            "assembly_support","assembly_score","long_reads_support","long_reads_score","short_reads_support","short_reads_score",
            "asm_copy_number_estimate","long_reads_copy_number_estimate","short_reads_copy_number_estimate"
        ],
        "Inversions": [
            "event_id","event_type","sequence_id","asm_start","asm_end","asm_length_bp",
            "long_reads_start_median","long_reads_start_std","long_reads_end_median","long_reads_end_std","long_reads_length_median",
            "short_reads_start_median","short_reads_start_std","short_reads_end_median","short_reads_end_std","short_reads_length_median",
            "assembly_support","assembly_score","long_reads_support","long_reads_score","short_reads_support","short_reads_score",
            "asm_copy_number_estimate","long_reads_copy_number_estimate","short_reads_copy_number_estimate"
        ],
        "Replacements": [
            "event_id","event_type","sequence_id",
            "asm_del_start","asm_del_end","asm_del_length","asm_ins_start","asm_ins_end","asm_ins_length",
            "long_reads_del_start_median","long_reads_del_start_std","long_reads_del_end_median","long_reads_del_end_std","long_reads_del_length_median",
            "long_reads_ins_start_median","long_reads_ins_start_std","long_reads_ins_end_median","long_reads_ins_end_std","long_reads_ins_length_median",
            "short_reads_del_start_median","short_reads_del_start_std","short_reads_del_end_median","short_reads_del_end_std","short_reads_del_length_median",
            "short_reads_ins_start_median","short_reads_ins_start_std","short_reads_ins_end_median","short_reads_ins_end_std","short_reads_ins_length_median",
            "assembly_support","assembly_score","long_reads_support","long_reads_score","short_reads_support","short_reads_score",
            "asm_copy_number_estimate","long_reads_copy_number_estimate","short_reads_copy_number_estimate"
        ],
        "Translocations": [
            "event_id","event_type","origin_sequence_id","destination_sequence_id",
            "asm_origin_start","asm_origin_end","asm_dest_start","asm_dest_end","asm_length",
            "long_reads_origin_start_median","long_reads_origin_start_std","long_reads_origin_end_median","long_reads_origin_end_std",
            "long_reads_dest_start_median","long_reads_dest_start_std","long_reads_dest_end_median","long_reads_dest_end_std","long_reads_length_median",
            "short_reads_origin_start_median","short_reads_origin_start_std","short_reads_origin_end_median","short_reads_origin_end_std",
            "short_reads_dest_start_median","short_reads_dest_start_std","short_reads_dest_end_median","short_reads_dest_end_std","short_reads_length_median",
            "assembly_support","assembly_score","long_reads_support","long_reads_score","short_reads_support","short_reads_score",
        ],
    }

    for name, keyrow in headers.items():
        ws = wb.create_sheet(title=name)
        ws.append(["" for _ in keyrow])  # row 1: descriptions blank
        ws.append(keyrow)                # row 2: keys

    return wb


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--assembly", required=True, help="SyRI SV summary TSV (syri.tsv)")
    ap.add_argument("--short", required=True, help="Short-read SV summary TSV (short.tsv)")
    ap.add_argument("--long", required=True, help="Long-read SV summary TSV (long.tsv)")
    ap.add_argument("--out", required=True, help="Output XLSX path")
    ap.add_argument("--template", required=False, default=None, help="Optional template XLSX (output.xlsx)")
    args = ap.parse_args()

    short_rows = enrich(read_tsv(args.short), "short")
    long_rows = enrich(read_tsv(args.long), "long")
    syri_rows = enrich_syri(read_tsv(args.assembly))

    short_by_type: Dict[str, List[Dict[str, Any]]] = defaultdict(list)
    for r in short_rows:
        if r.get("type") in KNOWN:
            short_by_type[r["type"]].append(r)

    long_by_type: Dict[str, List[Dict[str, Any]]] = defaultdict(list)
    for r in long_rows:
        if r.get("type") in KNOWN:
            long_by_type[r["type"]].append(r)
    for t in long_by_type:
        long_by_type[t].sort(key=lambda x: (x.get("chrom") or "", x.get("start") or -1))

    syri_by_type: Dict[str, List[Dict[str, Any]]] = defaultdict(list)
    for r in syri_rows:
        if r.get("type") in KNOWN:
            syri_by_type[r["type"]].append(r)

    sheet_for_type = {"INS": "Insertions", "DEL": "Deletions", "INV": "Inversions", "TRA": "Translocations", "REP": "Replacements"}
    prefix_for_type = {"INS": "INS", "DEL": "DEL", "INV": "INV", "TRA": "TRL", "REP": "REP"}  # matches example style

    events_by_sheet: Dict[str, List[Dict[str, Any]]] = defaultdict(list)

    for svt, long_list in long_by_type.items():
        sheet = sheet_for_type[svt]
        prefix = prefix_for_type[svt]

        for idx, longr in enumerate(long_list, 1):
            eid = f"{prefix}_{idx}"
            shorts = [sr for sr in short_by_type.get(svt, []) if overlaps(sr, longr)]

            asm = best_asm_match(longr, svt, syri_by_type)
            if asm:
                asm_start, asm_end, asm_len = asm.get("start"), asm.get("end"), asm.get("length")
            else:
                asm_start, asm_end, asm_len = longr.get("start"), longr.get("end"), longr.get("length")

            row: Dict[str, Any] = {}

            if sheet in ("Insertions", "Deletions", "Inversions"):
                row.update({
                    "event_id": eid,
                    "event_type": svt,
                    "sequence_id": longr.get("chrom"),
                    "asm_start": asm_start,
                    "asm_end": asm_end,
                    "asm_length_bp": asm_len,
                    "long_reads_start_median": longr.get("start"),
                    "long_reads_start_std": 0.0,
                    "long_reads_end_median": longr.get("end"),
                    "long_reads_end_std": 0.0,
                    "long_reads_length_median": longr.get("length"),
                    "long_reads_support": longr.get("support"),
                    "long_reads_score": None,
                })

                smed, sstd = stats([s.get("start") for s in shorts])
                emed, estd = stats([s.get("end") for s in shorts])
                lmed, _ = stats([s.get("length") for s in shorts])
                supmed, _ = stats([s.get("support") for s in shorts])

                row.update({
                    "short_reads_start_median": smed,
                    "short_reads_start_std": sstd,
                    "short_reads_end_median": emed,
                    "short_reads_end_std": estd,
                    "short_reads_length_median": lmed,
                    "short_reads_support": supmed,
                    "short_reads_score": None,
                })

            elif sheet == "Replacements":
                row.update({
                    "event_id": eid,
                    "event_type": svt,
                    "sequence_id": longr.get("chrom"),
                    "asm_del_start": asm_start,
                    "asm_del_end": asm_end,
                    "asm_del_length": asm_len,
                    "asm_ins_start": None,
                    "asm_ins_end": None,
                    "asm_ins_length": None,
                    "long_reads_del_start_median": longr.get("start"),
                    "long_reads_del_start_std": 0.0,
                    "long_reads_del_end_median": longr.get("end"),
                    "long_reads_del_end_std": 0.0,
                    "long_reads_del_length_median": longr.get("length"),
                    "long_reads_ins_start_median": None,
                    "long_reads_ins_start_std": None,
                    "long_reads_ins_end_median": None,
                    "long_reads_ins_end_std": None,
                    "long_reads_ins_length_median": None,
                    "long_reads_support": longr.get("support"),
                    "long_reads_score": None,
                })

                smed, sstd = stats([s.get("start") for s in shorts])
                emed, estd = stats([s.get("end") for s in shorts])
                lmed, _ = stats([s.get("length") for s in shorts])
                supmed, _ = stats([s.get("support") for s in shorts])

                row.update({
                    "short_reads_del_start_median": smed,
                    "short_reads_del_start_std": sstd,
                    "short_reads_del_end_median": emed,
                    "short_reads_del_end_std": estd,
                    "short_reads_del_length_median": lmed,
                    "short_reads_ins_start_median": None,
                    "short_reads_ins_start_std": None,
                    "short_reads_ins_end_median": None,
                    "short_reads_ins_end_std": None,
                    "short_reads_ins_length_median": None,
                    "short_reads_support": supmed,
                    "short_reads_score": None,
                })

            elif sheet == "Translocations":
                row.update({
                    "event_id": eid,
                    "event_type": svt,
                    "origin_sequence_id": longr.get("chrom"),
                    "destination_sequence_id": None,
                    "asm_origin_start": asm_start,
                    "asm_origin_end": asm_end,
                    "asm_dest_start": None,
                    "asm_dest_end": None,
                    "asm_length": asm_len,
                    "long_reads_origin_start_median": longr.get("start"),
                    "long_reads_origin_start_std": 0.0,
                    "long_reads_origin_end_median": longr.get("end"),
                    "long_reads_origin_end_std": 0.0,
                    "long_reads_dest_start_median": None,
                    "long_reads_dest_start_std": None,
                    "long_reads_dest_end_median": None,
                    "long_reads_dest_end_std": None,
                    "long_reads_length_median": longr.get("length"),
                    "long_reads_support": longr.get("support"),
                    "long_reads_score": None,
                })

                smed, sstd = stats([s.get("start") for s in shorts])
                emed, estd = stats([s.get("end") for s in shorts])
                lmed, _ = stats([s.get("length") for s in shorts])
                supmed, _ = stats([s.get("support") for s in shorts])

                row.update({
                    "short_reads_origin_start_median": smed,
                    "short_reads_origin_start_std": sstd,
                    "short_reads_origin_end_median": emed,
                    "short_reads_origin_end_std": estd,
                    "short_reads_dest_start_median": None,
                    "short_reads_dest_start_std": None,
                    "short_reads_dest_end_median": None,
                    "short_reads_dest_end_std": None,
                    "short_reads_length_median": lmed,
                    "short_reads_support": supmed,
                    "short_reads_score": None,
                })

            events_by_sheet[sheet].append(row)

    wb = load_template_or_create(args.template)

    for sheet_name in wb.sheetnames:
        ws = wb[sheet_name]

        # Determine key row from template row 2 (preferred)
        max_col = ws.max_column
        key_row = [ws.cell(row=2, column=c).value for c in range(1, max_col + 1)]
        key_row = [k for k in key_row if k is not None]

        # Clear rows 3+ only (keep header rows 1-2)
        if ws.max_row > 2:
            ws.delete_rows(3, ws.max_row - 2)

        for r in events_by_sheet.get(sheet_name, []):
            ws.append([r.get(k) for k in key_row])

    wb.save(args.out)


if __name__ == "__main__":
    main()
