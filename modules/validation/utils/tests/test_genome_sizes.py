"""
Tests for ref_genome_size_bp / mod_genome_size_bp calculation.

Covers two layers:
  1. build_params() correctly extracts genome sizes from validator metadata
     (strict mode → total_genome_size; trust mode → sum(sequence_lengths);
      minimal mode → None).
  2. The mod_genome_size correction block in main.py produces chromosome-only
     sizes by using contig_orientations from GXG results.

The correction logic is replicated here as _apply_mod_size_correction() so
the tests are self-contained and don't need to invoke main().
"""

import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

_VALIDATION_ROOT = Path(__file__).parent.parent.parent  # modules/validation/
if str(_VALIDATION_ROOT) not in sys.path:
    sys.path.insert(0, str(_VALIDATION_ROOT))

from nextflow_params_handler import build_params, NextflowParams


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _genome_meta(
    total_genome_size=None,
    sequence_lengths=None,
    fragmented=False,
    output_file="/genome.fasta",
):
    return SimpleNamespace(
        output_file=output_file,
        total_genome_size=total_genome_size,
        sequence_lengths=sequence_lengths,
        fragmented=fragmented,
    )


def _gxg(contig_orientations=None, contig_files=None, passed=True, plasmid_file=None):
    return {
        "passed": passed,
        "metadata": {
            "contig_files": contig_files or [],
            "contig_orientations": contig_orientations or {},
            "plasmid_file": plasmid_file,
        },
    }


def _results(ref=None, mod=None, gxg=None):
    return {
        "ref_genome": ref or _genome_meta(output_file="/ref/genome.fasta"),
        "mod_genome": mod,
        "ref_plasmid": None,
        "mod_plasmid": None,
        "genomexgenome": gxg,
        "reads": [],
        "ref_feature": None,
    }


def _apply_mod_size_correction(genomexgenome_res, mod_genome_res):
    """Replicate the correction block from main.py (must stay in sync)."""
    if genomexgenome_res is None or mod_genome_res is None:
        return
    gxg_meta = genomexgenome_res.get("metadata") or {}
    chr_ids = set((gxg_meta.get("contig_orientations") or {}).keys())
    seq_lengths = getattr(mod_genome_res, "sequence_lengths", None) or {}
    if chr_ids and seq_lengths:
        mod_chr_size = sum(seq_lengths[sid] for sid in chr_ids if sid in seq_lengths)
        if mod_chr_size > 0:
            mod_genome_res.total_genome_size = mod_chr_size


# ---------------------------------------------------------------------------
# build_params: genome size extraction
# ---------------------------------------------------------------------------

class TestGenomeSizeExtraction:

    def test_strict_mode_uses_total_genome_size(self):
        ref = _genome_meta(total_genome_size=4_500_000, output_file="/ref/genome.fasta")
        p = build_params(_results(ref=ref))
        assert p.ref_genome_size_bp == 4_500_000

    def test_trust_mode_sums_sequence_lengths(self):
        ref = _genome_meta(
            total_genome_size=None,
            sequence_lengths={"chr": 3_000_000, "chr1": 1_500_000},
            output_file="/ref/genome.fasta",
        )
        p = build_params(_results(ref=ref))
        assert p.ref_genome_size_bp == 4_500_000

    def test_minimal_mode_returns_none(self):
        ref = _genome_meta(total_genome_size=None, sequence_lengths=None, output_file="/ref/genome.fasta")
        p = build_params(_results(ref=ref))
        assert p.ref_genome_size_bp is None

    def test_size_absent_from_json_when_none(self):
        ref = _genome_meta(total_genome_size=None, sequence_lengths=None, output_file="/ref/genome.fasta")
        p = build_params(_results(ref=ref))
        d = p.to_dict()
        assert "ref_genome_size_bp" not in d
        assert "mod_genome_size_bp" not in d

    def test_both_ref_and_mod_sizes_present(self):
        ref = _genome_meta(total_genome_size=4_000_000, output_file="/ref/genome.fasta")
        mod = _genome_meta(total_genome_size=3_900_000, output_file="/mod/genome.fasta")
        p = build_params(_results(ref=ref, mod=mod, gxg=_gxg(passed=True)))
        assert p.ref_genome_size_bp == 4_000_000
        assert p.mod_genome_size_bp == 3_900_000

    def test_fragmented_mod_size_is_total_assembly(self):
        """Scenarios 3/4: GXG did not run, size = total of all parsed contigs."""
        mod = _genome_meta(
            sequence_lengths={"chr": 2_000_000, "chr1": 1_500_000, "chr2": 300_000},
            fragmented=True,
            output_file="/mod/genome.fasta",
        )
        p = build_params(_results(mod=mod))
        assert p.mod_genome_size_bp == 3_800_000

    def test_no_mod_genome_gives_none(self):
        p = build_params(_results(mod=None))
        assert p.mod_genome_size_bp is None


# ---------------------------------------------------------------------------
# mod_genome_size correction (main.py logic)
# ---------------------------------------------------------------------------

class TestModGenomeSizeCorrection:

    def test_scenario1_single_chromosome_plus_plasmid(self):
        """Scenario 1: 1 chromosomal contig + 1 plasmid → chromosome only."""
        mod = _genome_meta(
            total_genome_size=5_500_000,
            sequence_lengths={"chr": 5_000_000, "chr1": 500_000},
        )
        gxg = _gxg(contig_orientations={"chr": "+"})
        _apply_mod_size_correction(gxg, mod)
        assert mod.total_genome_size == 5_000_000

    def test_scenario2_multiple_chromosome_contigs(self):
        """Scenario 2: 2 chromosomal contigs + 1 plasmid → sum of chromosome contigs."""
        mod = _genome_meta(
            total_genome_size=8_500_000,
            sequence_lengths={"chr": 5_000_000, "chr1": 3_000_000, "chr2": 500_000},
        )
        gxg = _gxg(contig_orientations={"chr": "+", "chr1": "-"})
        _apply_mod_size_correction(gxg, mod)
        assert mod.total_genome_size == 8_000_000

    def test_reverse_complement_contig_included(self):
        """Reverse-complemented contigs (strand='-') count toward chromosome size."""
        mod = _genome_meta(
            sequence_lengths={"chr": 4_000_000, "chr1": 600_000},
        )
        gxg = _gxg(contig_orientations={"chr": "+", "chr1": "-"})
        _apply_mod_size_correction(gxg, mod)
        assert mod.total_genome_size == 4_600_000

    def test_trust_mode_sets_total_from_none(self):
        """Trust mode: total_genome_size starts None; correction sets it."""
        mod = _genome_meta(
            total_genome_size=None,
            sequence_lengths={"chr": 5_000_000, "chr1": 500_000},
        )
        gxg = _gxg(contig_orientations={"chr": "+"})
        _apply_mod_size_correction(gxg, mod)
        assert mod.total_genome_size == 5_000_000

    def test_no_correction_when_gxg_not_run(self):
        """Scenarios 3/4: GXG=None, size unchanged (total assembly)."""
        mod = _genome_meta(
            total_genome_size=10_000_000,
            sequence_lengths={"chr": 6_000_000, "chr1": 4_000_000},
            fragmented=True,
        )
        _apply_mod_size_correction(None, mod)
        assert mod.total_genome_size == 10_000_000

    def test_no_correction_when_all_sequences_unmapped(self):
        """All sequences are plasmids (no contig_orientations entries): unchanged."""
        mod = _genome_meta(
            total_genome_size=500_000,
            sequence_lengths={"chr": 500_000},
        )
        gxg = _gxg(contig_orientations={})
        _apply_mod_size_correction(gxg, mod)
        assert mod.total_genome_size == 500_000

    def test_no_correction_when_mod_genome_res_none(self):
        gxg = _gxg(contig_orientations={"chr": "+"})
        _apply_mod_size_correction(gxg, None)  # must not raise

    def test_corrected_value_flows_into_build_params(self):
        """End-to-end: correction + build_params → correct mod_genome_size_bp in JSON."""
        mod = _genome_meta(
            total_genome_size=8_500_000,
            sequence_lengths={"chr": 5_000_000, "chr1": 3_000_000, "chr2": 500_000},
            output_file="/mod/genome.fasta",
        )
        gxg = _gxg(contig_orientations={"chr": "+", "chr1": "-"}, passed=True)
        _apply_mod_size_correction(gxg, mod)

        ref = _genome_meta(total_genome_size=4_000_000, output_file="/ref/genome.fasta")
        p = build_params(_results(ref=ref, mod=mod, gxg=gxg))
        assert p.mod_genome_size_bp == 8_000_000
        assert p.to_dict()["mod_genome_size_bp"] == 8_000_000

    def test_ref_size_unchanged_by_correction(self):
        """Correction must never touch ref_genome_size_bp."""
        ref = _genome_meta(total_genome_size=4_000_000, output_file="/ref/genome.fasta")
        mod = _genome_meta(
            total_genome_size=5_500_000,
            sequence_lengths={"chr": 5_000_000, "chr1": 500_000},
            output_file="/mod/genome.fasta",
        )
        gxg = _gxg(contig_orientations={"chr": "+"}, passed=True)
        _apply_mod_size_correction(gxg, mod)
        p = build_params(_results(ref=ref, mod=mod, gxg=gxg))
        assert p.ref_genome_size_bp == 4_000_000
        assert p.mod_genome_size_bp == 5_000_000
