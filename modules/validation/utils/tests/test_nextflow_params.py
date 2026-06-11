"""
Tests for nextflow_params.py

Tests cover:
- run_ref_x_mod conditions
- Read type detection (illumina, ont, pacbio)
- Conditional keys (ref_fasta_validated, mod_fasta_validated, pacbio_fastqs)
- write_params serialises valid JSON
"""

import json
import tempfile
from pathlib import Path
from types import SimpleNamespace

import sys
from pathlib import Path

_VALIDATION_ROOT = Path(__file__).parent.parent.parent  # modules/validation/
if str(_VALIDATION_ROOT) not in sys.path:
    sys.path.insert(0, str(_VALIDATION_ROOT))

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


def _gxg(contig_files: list, passed: bool = True, plasmid_file: str = None) -> dict:
    return {"metadata": {"contig_files": contig_files, "plasmid_file": plasmid_file}, "passed": passed}


def _base(contig_files=None, passed=True, plasmid_file=None):
    """Minimal validation_results with ref + mod genome and given contigs."""
    return {
        "ref_genome":    _meta("/ref/genome.fasta"),
        "mod_genome":    _meta("/mod/genome.fasta"),
        "genomexgenome": _gxg(contig_files or [], passed=passed, plasmid_file=plasmid_file),
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
        assert "/reads/nano.fastq" in p.ont_fastqs
        assert p.run_illumina is False

    def test_pacbio_read(self):
        r = _base()
        r["reads"] = [_meta("/reads/pb.fastq", ngs_type="pacbio")]
        p = build_params(r)
        assert p.run_pacbio is True
        assert "/reads/pb.fastq" in p.pacbio_fastqs

    def test_no_pacbio_key_when_absent(self):
        r = _base()
        r["reads"] = [_meta("/reads/R1.fastq", ngs_type="illumina")]
        p = build_params(r)
        assert p.pacbio_fastqs == []

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
        assert p.ont_fastqs == []


# ---------------------------------------------------------------------------
# Plasmid paths
# ---------------------------------------------------------------------------

class TestPlasmidPaths:

    def test_explicit_plasmid_config(self):
        r = _base()
        r["ref_plasmid"] = _meta("/valid/ref_plasmid.fasta")
        r["mod_plasmid"] = _meta("/valid/mod_plasmid.fasta")
        p = build_params(r)
        assert p.ref_plasmid_fasta == "/valid/ref_plasmid.fasta"
        assert p.mod_plasmid_fasta == "/valid/mod_plasmid.fasta"

    def test_ref_plasmid_from_genome_validator_split(self):
        """ref_plasmid_fasta falls back to plasmid_filenames on ref_genome when no explicit config."""
        r = _base()
        r["ref_genome"] = SimpleNamespace(
            output_file="/valid/ref.fasta",
            plasmid_filenames=["/valid/ref_plasmid.fasta"],
        )
        p = build_params(r)
        assert p.ref_plasmid_fasta == "/valid/ref_plasmid.fasta"

    def test_mod_plasmid_from_gxg_characterisation(self):
        """mod_plasmid_fasta falls back to plasmid_file in GXG metadata when no explicit config."""
        r = _base(plasmid_file="/valid/mod_plasmid.fasta")
        p = build_params(r)
        assert p.mod_plasmid_fasta == "/valid/mod_plasmid.fasta"

    def test_explicit_config_takes_priority_over_fallback(self):
        r = _base(plasmid_file="/gxg/mod_plasmid.fasta")
        r["mod_plasmid"] = _meta("/explicit/mod_plasmid.fasta")
        p = build_params(r)
        assert p.mod_plasmid_fasta == "/explicit/mod_plasmid.fasta"

    def test_no_plasmids(self):
        p = build_params(_base())
        assert p.ref_plasmid_fasta is None
        assert p.mod_plasmid_fasta is None

    # -- regression tests: output_file set by is_plasmid validator fix --------

    def test_mod_plasmid_from_output_file_not_fallback(self):
        """mod_plasmid_fasta is read from mod_plasmid.output_file, not from the
        gxg plasmid_file fallback — this is the scenario fixed by setting
        output_file correctly in is_plasmid validator mode."""
        gxg_plasmid = "/gxg/should_not_be_used.fasta"
        r = _base(plasmid_file=gxg_plasmid)
        r["mod_plasmid"] = _meta("/valid/mod_plasmid_mod_plasmid.fasta")
        p = build_params(r)
        assert p.mod_plasmid_fasta == "/valid/mod_plasmid_mod_plasmid.fasta"

    def test_ref_plasmid_from_output_file_not_genome_fallback(self):
        """ref_plasmid_fasta is read from ref_plasmid.output_file, not from
        ref_genome.plasmid_filenames — ensures explicit config takes precedence
        over the genome-extraction fallback."""
        r = _base()
        r["ref_genome"] = SimpleNamespace(
            output_file="/valid/ref.fasta",
            plasmid_filenames=["/valid/genome_extracted_plasmid.fasta"],
        )
        r["ref_plasmid"] = _meta("/valid/ref_plasmid_ref_plasmid.fasta")
        p = build_params(r)
        assert p.ref_plasmid_fasta == "/valid/ref_plasmid_ref_plasmid.fasta"

    def test_mod_plasmid_absent_from_json_when_output_file_none_and_no_fallback(self):
        """When mod_plasmid.output_file is None and no gxg plasmid_file exists,
        mod_plasmid_fasta must be absent from the serialised JSON."""
        r = _base()  # gxg plasmid_file is None, no mod_plasmid key
        r["mod_plasmid"] = SimpleNamespace(output_file=None)
        p = build_params(r)
        d = p.to_dict()
        assert "mod_plasmid_fasta" not in d

    def test_mod_plasmid_present_in_json_when_output_file_set(self):
        """When mod_plasmid.output_file is set, mod_plasmid_fasta must appear in
        the serialised JSON — end-to-end check of to_dict()."""
        r = _base()
        r["mod_plasmid"] = _meta("/valid/mod_plasmid.fasta")
        p = build_params(r)
        d = p.to_dict()
        assert "mod_plasmid_fasta" in d
        assert d["mod_plasmid_fasta"] == "/valid/mod_plasmid.fasta"


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
