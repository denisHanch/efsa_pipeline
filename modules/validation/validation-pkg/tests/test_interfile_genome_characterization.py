"""
Tests for genome characterization (contigs vs plasmids via minimap2),
which is now integrated into genomexgenome_validation.

Tests cover:
- characterize=True/False on GenomeXGenomeSettings
- Successful run with mixed contigs and plasmids
- Edge cases: all mapped, none mapped
- minimap2 unavailability (demoted to warning, validation still passes)
- Correct sequence content in output files
"""

import pytest
from pathlib import Path
from unittest.mock import patch, MagicMock
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from validation_pkg.validators.interfile_genome import (
    GenomeXGenomeSettings,
    genomexgenome_validation,
)
from validation_pkg.validators.genome_validator import GenomeOutputMetadata
from validation_pkg.exceptions import GenomeValidationError


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_records(names):
    """Create simple SeqRecord objects."""
    return [
        SeqRecord(Seq("ACGT" * 25), id=name, description="")
        for name in names
    ]


def _write_fasta(path: Path, records):
    with open(path, "w") as fh:
        SeqIO.write(records, fh, "fasta")


def _read_ids(path: Path):
    return [r.id for r in SeqIO.parse(str(path), "fasta")]


def _make_paf_line(query_id: str) -> str:
    return (
        f"{query_id}\t100\t0\t100\t+\t"
        "ref_chr1\t1000000\t0\t100\t100\t100\t255"
    )


# ---------------------------------------------------------------------------
# Settings tests
# ---------------------------------------------------------------------------

class TestGenomeXGenomeSettingsCharacterize:
    def test_characterize_default_true(self):
        assert GenomeXGenomeSettings().characterize is True

    def test_characterize_disabled(self):
        assert GenomeXGenomeSettings(characterize=False).characterize is False

    def test_update_preserves_original(self):
        s = GenomeXGenomeSettings()
        updated = s.update(characterize=False)
        assert updated.characterize is False
        assert s.characterize is True  # original unchanged


# ---------------------------------------------------------------------------
# Core characterization tests (via genomexgenome_validation)
# ---------------------------------------------------------------------------

class TestCharacterizationInGenomeXGenome:

    @pytest.fixture
    def fasta_setup(self, tmp_path):
        """Create ref and mod FASTA files; disable sequence-count check."""
        ref_seqs = _make_records(["chr1", "chr2"])
        mod_seqs = _make_records(["chr1", "chr2", "plasmid1"])

        ref_path = tmp_path / "ref.fasta"
        mod_path = tmp_path / "mod.fasta"
        _write_fasta(ref_path, ref_seqs)
        _write_fasta(mod_path, mod_seqs)

        ref_result = GenomeOutputMetadata(
            output_file=str(ref_path),
            output_filename="ref.fasta",
            num_sequences=2,
        )
        mod_result = GenomeOutputMetadata(
            output_file=str(mod_path),
            output_filename="mod.fasta",
            num_sequences=3,
        )
        # Disable count check so mismatched counts don't raise
        settings = GenomeXGenomeSettings(same_number_of_sequences=False, characterize=True)
        return tmp_path, ref_result, mod_result, settings

    # ------------------------------------------------------------------
    # Successful run: 2 contigs + 1 plasmid
    # ------------------------------------------------------------------

    def test_successful_run(self, fasta_setup):
        tmp_path, ref_result, mod_result, settings = fasta_setup
        paf = "\n".join([_make_paf_line("chr1"), _make_paf_line("chr2")])

        with patch("validation_pkg.validators.interfile_genome.check_tool_available", return_value=True), \
             patch("validation_pkg.validators.interfile_genome.subprocess.run") as mock_run:
            mock_run.return_value = MagicMock(stdout=paf, returncode=0)
            result = genomexgenome_validation(ref_result, mod_result, settings)

        assert result["passed"] is True
        assert result["metadata"]["contigs_found"] == 2
        assert result["metadata"]["plasmids_found"] == 1
        assert len(result["metadata"]["contig_files"]) == 2
        assert result["metadata"]["plasmid_file"] is not None

    # ------------------------------------------------------------------
    # All sequences map → 0 plasmids
    # ------------------------------------------------------------------

    def test_all_mapped_no_plasmids(self, fasta_setup):
        tmp_path, ref_result, mod_result, settings = fasta_setup
        paf = "\n".join([
            _make_paf_line("chr1"),
            _make_paf_line("chr2"),
            _make_paf_line("plasmid1"),
        ])

        with patch("validation_pkg.validators.interfile_genome.check_tool_available", return_value=True), \
             patch("validation_pkg.validators.interfile_genome.subprocess.run") as mock_run:
            mock_run.return_value = MagicMock(stdout=paf, returncode=0)
            result = genomexgenome_validation(ref_result, mod_result, settings)

        assert result["metadata"]["contigs_found"] == 3
        assert result["metadata"]["plasmids_found"] == 0
        assert result["metadata"]["plasmid_file"] is None

    # ------------------------------------------------------------------
    # No sequences map → all plasmids
    # ------------------------------------------------------------------

    def test_no_mapped_all_plasmids(self, fasta_setup):
        tmp_path, ref_result, mod_result, settings = fasta_setup

        with patch("validation_pkg.validators.interfile_genome.check_tool_available", return_value=True), \
             patch("validation_pkg.validators.interfile_genome.subprocess.run") as mock_run:
            mock_run.return_value = MagicMock(stdout="", returncode=0)
            result = genomexgenome_validation(ref_result, mod_result, settings)

        assert result["metadata"]["contigs_found"] == 0
        assert result["metadata"]["plasmids_found"] == 3
        assert result["metadata"]["contig_files"] == []
        assert result["metadata"]["plasmid_file"] is not None

    # ------------------------------------------------------------------
    # minimap2 not available → demoted to warning, validation still passes
    # ------------------------------------------------------------------

    def test_minimap2_not_available_becomes_warning(self, fasta_setup):
        tmp_path, ref_result, mod_result, settings = fasta_setup

        with patch("validation_pkg.validators.interfile_genome.check_tool_available", return_value=False):
            result = genomexgenome_validation(ref_result, mod_result, settings)

        assert result["passed"] is True
        assert any("minimap2" in w for w in result["warnings"])
        assert result["metadata"]["contigs_found"] is None

    # ------------------------------------------------------------------
    # minimap2 subprocess failure → demoted to warning
    # ------------------------------------------------------------------

    def test_minimap2_subprocess_failure_becomes_warning(self, fasta_setup):
        import subprocess as _subprocess
        tmp_path, ref_result, mod_result, settings = fasta_setup

        with patch("validation_pkg.validators.interfile_genome.check_tool_available", return_value=True), \
             patch("validation_pkg.validators.interfile_genome.subprocess.run",
                   side_effect=_subprocess.CalledProcessError(1, "minimap2", stderr="error")):
            result = genomexgenome_validation(ref_result, mod_result, settings)

        assert result["passed"] is True
        assert any("minimap2 failed" in w for w in result["warnings"])

    # ------------------------------------------------------------------
    # characterize=False → minimap2 never called
    # ------------------------------------------------------------------

    def test_characterize_disabled(self, fasta_setup):
        tmp_path, ref_result, mod_result, _ = fasta_setup
        settings = GenomeXGenomeSettings(same_number_of_sequences=False, characterize=False)

        with patch("validation_pkg.validators.interfile_genome.subprocess.run") as mock_run:
            result = genomexgenome_validation(ref_result, mod_result, settings)

        mock_run.assert_not_called()
        assert "contigs_found" not in result["metadata"]

    # ------------------------------------------------------------------
    # Output files contain correct sequences
    # ------------------------------------------------------------------

    def test_output_files_correct_sequences(self, fasta_setup):
        tmp_path, ref_result, mod_result, settings = fasta_setup
        paf = "\n".join([_make_paf_line("chr1"), _make_paf_line("chr2")])

        with patch("validation_pkg.validators.interfile_genome.check_tool_available", return_value=True), \
             patch("validation_pkg.validators.interfile_genome.subprocess.run") as mock_run:
            mock_run.return_value = MagicMock(stdout=paf, returncode=0)
            result = genomexgenome_validation(ref_result, mod_result, settings)

        for contig_file in result["metadata"]["contig_files"]:
            ids = _read_ids(Path(contig_file))
            assert len(ids) == 1
            assert ids[0] in {"chr1", "chr2"}

        plasmid_ids = _read_ids(Path(result["metadata"]["plasmid_file"]))
        assert plasmid_ids == ["plasmid1"]
