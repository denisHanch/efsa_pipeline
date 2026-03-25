"""
Tests for nextflow_params.py

Tests cover:
- run_ref_x_mod conditions
- Read type detection (illumina, ont, pacbio)
- GFF / run_vcf_annotation
- Conditional keys (ref_fasta_validated, mod_fasta_validated, pacbio_fastq, gff)
- write_params serialises valid JSON
"""

import json
import tempfile
from pathlib import Path
from types import SimpleNamespace

import pytest

from nextflow_params_handler import build_params, write_params


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _meta(path: str, ngs_type: str = None) -> SimpleNamespace:
    """Minimal metadata object with output_file (and optional ngs_type)."""
    m = SimpleNamespace(output_file=path)
    if ngs_type is not None:
        m.ngs_type = ngs_type
    return m


def _gxg(contig_files: list, passed: bool = True) -> dict:
    return {"metadata": {"contig_files": contig_files}, "passed": passed}


def _base(contig_files=None, passed=True):
    """Minimal validation_results with ref + mod genome and given contigs."""
    return {
        "ref_genome":    _meta("/ref/genome.fasta"),
        "mod_genome":    _meta("/mod/genome.fasta"),
        "genomexgenome": _gxg(contig_files or [], passed=passed),
        "reads":         [],
        "ref_feature":   None,
    }


# ---------------------------------------------------------------------------
# run_ref_x_mod
# ---------------------------------------------------------------------------

class TestRunRefXMod:

    def test_both_present_and_passed(self):
        p = build_params(_base(contig_files=["c1"], passed=True))
        assert p.run_ref_x_mod is True

    def test_gxg_not_passed(self):
        p = build_params(_base(contig_files=["c1"], passed=False))
        assert p.run_ref_x_mod is False

    def test_no_mod_genome(self):
        r = _base()
        r["mod_genome"] = None
        p = build_params(r)
        assert p.run_ref_x_mod is False

    def test_no_ref_genome(self):
        r = _base()
        r["ref_genome"] = None
        p = build_params(r)
        assert p.run_ref_x_mod is False

    def test_no_genomexgenome(self):
        r = _base()
        r["genomexgenome"] = None
        p = build_params(r)
        assert p.run_ref_x_mod is False


# ---------------------------------------------------------------------------
# FASTA conditional keys
# ---------------------------------------------------------------------------

class TestFastaAvailability:

    def test_ref_and_mod_available(self):
        p = build_params(_base())
        assert p.ref_fasta_validated == "/ref/genome.fasta"
        assert p.mod_fasta_validated == "/mod/genome.fasta"

    def test_ref_only(self):
        r = _base()
        r["mod_genome"] = None
        p = build_params(r)
        assert p.ref_fasta_validated is not None
        assert p.mod_fasta_validated is None

    def test_no_genomes(self):
        r = _base()
        r["ref_genome"] = None
        r["mod_genome"] = None
        p = build_params(r)
        assert p.ref_fasta_validated is None
        assert p.mod_fasta_validated is None


# ---------------------------------------------------------------------------
# Reads
# ---------------------------------------------------------------------------

class TestReads:

    def test_illumina_read(self):
        r = _base()
        r["reads"] = [_meta("/reads/R1.fastq", ngs_type="illumina")]
        p = build_params(r)
        assert p.run_illumina is True
        assert p.run_nanopore is False
        assert p.run_pacbio is False

    def test_ont_read(self):
        r = _base()
        r["reads"] = [_meta("/reads/nano.fastq", ngs_type="ont")]
        p = build_params(r)
        assert p.run_nanopore is True
        assert p.nanopore_fastq == "/reads/nano.fastq"
        assert p.run_illumina is False

    def test_pacbio_read(self):
        r = _base()
        r["reads"] = [_meta("/reads/pb.fastq", ngs_type="pacbio")]
        p = build_params(r)
        assert p.run_pacbio is True
        assert p.pacbio_fastq == "/reads/pb.fastq"

    def test_no_pacbio_key_when_absent(self):
        r = _base()
        r["reads"] = [_meta("/reads/R1.fastq", ngs_type="illumina")]
        p = build_params(r)
        assert p.pacbio_fastq is None

    def test_multiple_reads_same_type_keeps_first(self):
        r = _base()
        r["reads"] = [
            _meta("/reads/ont1.fastq", ngs_type="ont"),
            _meta("/reads/ont2.fastq", ngs_type="ont"),
        ]
        p = build_params(r)
        assert p.nanopore_fastq == "/reads/ont1.fastq"

    def test_mixed_read_types(self):
        r = _base()
        r["reads"] = [
            _meta("/reads/R1.fastq", ngs_type="illumina"),
            _meta("/reads/nano.fastq", ngs_type="ont"),
            _meta("/reads/pb.fastq", ngs_type="pacbio"),
        ]
        p = build_params(r)
        assert p.run_illumina is True
        assert p.run_nanopore is True
        assert p.run_pacbio is True

    def test_no_reads(self):
        p = build_params(_base())
        assert p.run_illumina is False
        assert p.run_nanopore is False
        assert p.run_pacbio is False
        assert p.nanopore_fastq is None


# ---------------------------------------------------------------------------
# GFF / feature annotation
# ---------------------------------------------------------------------------

class TestFeatureAnnotation:

    def test_gff_present(self):
        r = _base()
        r["ref_feature"] = _meta("/features/ref.gff")
        p = build_params(r)
        assert p.run_vcf_annotation is True
        assert p.gff == "/features/ref.gff"

    def test_no_gff(self):
        p = build_params(_base())
        assert p.run_vcf_annotation is False
        assert p.gff is None


# ---------------------------------------------------------------------------
# run_truvari always False
# ---------------------------------------------------------------------------

class TestRunTruvari:

    def test_always_false(self):
        p = build_params(_base())
        assert p.run_truvari is False


# ---------------------------------------------------------------------------
# write_params
# ---------------------------------------------------------------------------

class TestWriteParams:

    def test_writes_valid_json(self):
        params = build_params(_base(contig_files=["c1"]))
        with tempfile.TemporaryDirectory() as tmpdir:
            out = Path(tmpdir) / "params.json"
            write_params(params, out)
            data = json.loads(out.read_text())
        assert data["run_ref_x_mod"] == params.run_ref_x_mod

    def test_output_is_readable_dict(self):
        params = build_params(_base(contig_files=["c1"]))
        with tempfile.TemporaryDirectory() as tmpdir:
            out = Path(tmpdir) / "params.json"
            write_params(params, out)
            loaded = json.loads(out.read_text())
        assert isinstance(loaded, dict)
        assert "run_ref_x_mod" in loaded
        assert "run_truvari" in loaded
