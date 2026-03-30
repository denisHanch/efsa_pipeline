"""Tests for genome characterization (contigs vs plasmids via minimap2)."""

import pytest
from pathlib import Path
from unittest.mock import patch, MagicMock
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from validation_pkg.validators.interfile_genome import (
    GenomeXGenomeSettings,
    genomexgenome_validation,
    _parse_paf_best_hits,
)
from validation_pkg.validators.genome_validator import GenomeOutputMetadata
from validation_pkg.exceptions import InterFileValidationError


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


def _make_paf_line(
    query_id: str,
    strand: str = "+",
    ref_name: str = "ref_chr1",
    alignment_len: int = 100,
    q_start: int = 0,
    r_start: int = 0,
) -> str:
    """Build a minimal valid PAF line (12 tab-separated columns)."""
    return (
        f"{query_id}\t{alignment_len}\t{q_start}\t{q_start + alignment_len}\t{strand}\t"
        f"{ref_name}\t1000000\t{r_start}\t{r_start + alignment_len}\t{alignment_len}\t{alignment_len}\t255"
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
    # minimap2 not available → added to errors, raises InterFileValidationError
    # ------------------------------------------------------------------

    def test_minimap2_not_available_raises(self, fasta_setup):
        tmp_path, ref_result, mod_result, settings = fasta_setup

        with patch("validation_pkg.validators.interfile_genome.check_tool_available", return_value=False):
            with pytest.raises(InterFileValidationError) as exc_info:
                genomexgenome_validation(ref_result, mod_result, settings)

        assert "minimap2" in str(exc_info.value)

    # ------------------------------------------------------------------
    # minimap2 subprocess failure → added to errors, raises InterFileValidationError
    # ------------------------------------------------------------------

    def test_minimap2_subprocess_failure_raises(self, fasta_setup):
        import subprocess as _subprocess
        tmp_path, ref_result, mod_result, settings = fasta_setup

        with patch("validation_pkg.validators.interfile_genome.check_tool_available", return_value=True), \
             patch("validation_pkg.validators.interfile_genome.subprocess.run",
                   side_effect=_subprocess.CalledProcessError(1, "minimap2", stderr="error")):
            with pytest.raises(InterFileValidationError) as exc_info:
                genomexgenome_validation(ref_result, mod_result, settings)

        assert "minimap2 failed" in str(exc_info.value)

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
            assert ids[0] == "chr"

        plasmid_ids = _read_ids(Path(result["metadata"]["plasmid_file"]))
        assert plasmid_ids == ["plasmid1"]


# ---------------------------------------------------------------------------
# Unit tests for _parse_paf_best_hits
# ---------------------------------------------------------------------------

class TestParsePafBestHits:

    def test_single_hit_forward(self):
        paf = _make_paf_line("chr1", strand="+", ref_name="ref1", alignment_len=500)
        hits = _parse_paf_best_hits(paf)
        assert hits["chr1"] == {"ref_name": "ref1", "strand": "+", "alignment_len": 500}

    def test_single_hit_reverse(self):
        paf = _make_paf_line("chr1", strand="-", ref_name="ref1", alignment_len=500)
        hits = _parse_paf_best_hits(paf)
        assert hits["chr1"]["strand"] == "-"

    def test_best_hit_selected_by_alignment_len(self):
        """When a query has multiple hits, the one with the largest alignment_len wins."""
        lines = [
            _make_paf_line("chr1", strand="+", ref_name="ref_small", alignment_len=200),
            _make_paf_line("chr1", strand="-", ref_name="ref_large", alignment_len=5000),
            _make_paf_line("chr1", strand="+", ref_name="ref_mid",   alignment_len=1000),
        ]
        hits = _parse_paf_best_hits("\n".join(lines))
        assert hits["chr1"]["ref_name"] == "ref_large"
        assert hits["chr1"]["strand"] == "-"
        assert hits["chr1"]["alignment_len"] == 5000

    def test_ribosomal_rna_many_small_hits(self):
        """Multiple ~6k ribosomal hits: only the largest alignment survives."""
        rrna_hits = [
            _make_paf_line("rRNA_contig", strand="+", ref_name=f"rrna_copy_{i}", alignment_len=5900 + i)
            for i in range(5)
        ]
        hits = _parse_paf_best_hits("\n".join(rrna_hits))
        assert hits["rRNA_contig"]["ref_name"] == "rrna_copy_4"
        assert hits["rRNA_contig"]["alignment_len"] == 5904

    def test_multiple_queries_independent(self):
        lines = [
            _make_paf_line("chr1", strand="+", ref_name="ref1", alignment_len=1000),
            _make_paf_line("chr2", strand="-", ref_name="ref2", alignment_len=2000),
        ]
        hits = _parse_paf_best_hits("\n".join(lines))
        assert hits["chr1"]["strand"] == "+"
        assert hits["chr2"]["strand"] == "-"

    def test_empty_output(self):
        assert _parse_paf_best_hits("") == {}

    def test_short_lines_skipped(self):
        """Lines with fewer than 12 columns are ignored."""
        assert _parse_paf_best_hits("chr1\t100\t0\t100\t+\tref1") == {}

    def test_duplicate_pair_same_starts_ok(self):
        """Same (query, ref) pair with identical q_start and r_start should not raise."""
        lines = [
            _make_paf_line("chr1", ref_name="ref1", alignment_len=100, q_start=0, r_start=0),
            _make_paf_line("chr1", ref_name="ref1", alignment_len=200, q_start=0, r_start=0),
        ]
        hits = _parse_paf_best_hits("\n".join(lines))
        assert hits["chr1"]["alignment_len"] == 200

    # def test_conflicting_q_start_raises(self):
    #     """Same (query, ref) pair with different q_start values must raise."""
    #     lines = [
    #         _make_paf_line("chr1", ref_name="ref1", q_start=0, r_start=0),
    #         _make_paf_line("chr1", ref_name="ref1", q_start=500, r_start=0),
    #     ]
    #     with pytest.raises(InterFileValidationError, match="q_start"):
    #         _parse_paf_best_hits("\n".join(lines))

    # def test_conflicting_r_start_raises(self):
    #     """Same (query, ref) pair with different r_start values must raise."""
    #     lines = [
    #         _make_paf_line("chr1", ref_name="ref1", q_start=0, r_start=0),
    #         _make_paf_line("chr1", ref_name="ref1", q_start=0, r_start=1000),
    #     ]
    #     with pytest.raises(InterFileValidationError, match="r_start"):
    #         _parse_paf_best_hits("\n".join(lines))

    def test_different_ref_names_not_conflicting(self):
        """Same query mapping to two different refs with different starts should not raise."""
        lines = [
            _make_paf_line("chr1", ref_name="ref1", q_start=0,   r_start=0),
            _make_paf_line("chr1", ref_name="ref2", q_start=500, r_start=999),
        ]
        hits = _parse_paf_best_hits("\n".join(lines))
        # Best hit is the one with larger alignment_len (both default to 100, so first wins)
        assert hits["chr1"] is not None


# ---------------------------------------------------------------------------
# Orientation: metadata and reverse-complement in output files
# ---------------------------------------------------------------------------

class TestOrientationHandling:

    @pytest.fixture
    def fasta_setup(self, tmp_path):
        # chr1 maps forward, chr2 maps reverse, plasmid1 is unmapped
        ref_seqs = _make_records(["chr1", "chr2"])
        mod_seqs = _make_records(["chr1", "chr2", "plasmid1"])
        ref_path = tmp_path / "ref.fasta"
        mod_path = tmp_path / "mod.fasta"
        _write_fasta(ref_path, ref_seqs)
        _write_fasta(mod_path, mod_seqs)
        ref_result = GenomeOutputMetadata(
            output_file=str(ref_path), output_filename="ref.fasta", num_sequences=2
        )
        mod_result = GenomeOutputMetadata(
            output_file=str(mod_path), output_filename="mod.fasta", num_sequences=3
        )
        settings = GenomeXGenomeSettings(same_number_of_sequences=False, characterize=True)
        return tmp_path, ref_result, mod_result, settings

    def test_orientation_stored_in_metadata(self, fasta_setup):
        tmp_path, ref_result, mod_result, settings = fasta_setup
        paf = "\n".join([
            _make_paf_line("chr1", strand="+"),
            _make_paf_line("chr2", strand="-"),
        ])
        with patch("validation_pkg.validators.interfile_genome.check_tool_available", return_value=True), \
             patch("validation_pkg.validators.interfile_genome.subprocess.run") as mock_run:
            mock_run.return_value = MagicMock(stdout=paf, returncode=0)
            result = genomexgenome_validation(ref_result, mod_result, settings)

        orientations = result["metadata"]["contig_orientations"]
        assert orientations["chr1"] == "+"
        assert orientations["chr2"] == "-"

    def test_ref_names_stored_in_metadata(self, fasta_setup):
        tmp_path, ref_result, mod_result, settings = fasta_setup
        paf = "\n".join([
            _make_paf_line("chr1", ref_name="ref_chr1"),
            _make_paf_line("chr2", ref_name="ref_chr2"),
        ])
        with patch("validation_pkg.validators.interfile_genome.check_tool_available", return_value=True), \
             patch("validation_pkg.validators.interfile_genome.subprocess.run") as mock_run:
            mock_run.return_value = MagicMock(stdout=paf, returncode=0)
            result = genomexgenome_validation(ref_result, mod_result, settings)

        ref_names = result["metadata"]["contig_ref_names"]
        assert ref_names["chr1"] == "ref_chr1"
        assert ref_names["chr2"] == "ref_chr2"

    def test_reverse_strand_sequence_is_reverse_complemented(self, fasta_setup):
        """A contig with '-' strand should be written as its reverse complement."""
        tmp_path, ref_result, mod_result, settings = fasta_setup

        # Use a non-palindromic sequence so RC is clearly different
        known_seq = "AAAACCCCGGGGTTTT"  # RC = AAAACCCCGGGGTTTT... let's pick better
        known_seq = "AAACCCGGG"         # RC = CCCGGGTTT
        expected_rc = "CCCGGGTTT"

        # Overwrite mod.fasta with a sequence we control for chr2
        mod_path = Path(tmp_path / "mod.fasta")
        records = [
            SeqRecord(Seq("ACGTACGT"), id="chr1", description=""),
            SeqRecord(Seq(known_seq),  id="chr2", description=""),
            SeqRecord(Seq("TTTTAAAA"), id="plasmid1", description=""),
        ]
        _write_fasta(mod_path, records)

        paf = "\n".join([
            _make_paf_line("chr1", strand="+"),
            _make_paf_line("chr2", strand="-"),
        ])
        with patch("validation_pkg.validators.interfile_genome.check_tool_available", return_value=True), \
             patch("validation_pkg.validators.interfile_genome.subprocess.run") as mock_run:
            mock_run.return_value = MagicMock(stdout=paf, returncode=0)
            result = genomexgenome_validation(ref_result, mod_result, settings)

        # Find the contig file for chr2 (id is renamed to 'chr' after digit-stripping)
        chr2_file = next(
            f for f in result["metadata"]["contig_files"]
            if "_contig_1" in f
        )
        written_seq = str(SeqIO.read(chr2_file, "fasta").seq)
        assert written_seq == expected_rc

    def test_forward_strand_sequence_unchanged(self, fasta_setup):
        """A contig with '+' strand should be written as-is."""
        tmp_path, ref_result, mod_result, settings = fasta_setup
        known_seq = "AAACCCGGG"
        mod_path = Path(tmp_path / "mod.fasta")
        records = [
            SeqRecord(Seq(known_seq), id="chr1", description=""),
            SeqRecord(Seq("ACGT"),    id="chr2", description=""),
            SeqRecord(Seq("TTTT"),    id="plasmid1", description=""),
        ]
        _write_fasta(mod_path, records)

        paf = "\n".join([
            _make_paf_line("chr1", strand="+"),
            _make_paf_line("chr2", strand="+"),
        ])
        with patch("validation_pkg.validators.interfile_genome.check_tool_available", return_value=True), \
             patch("validation_pkg.validators.interfile_genome.subprocess.run") as mock_run:
            mock_run.return_value = MagicMock(stdout=paf, returncode=0)
            result = genomexgenome_validation(ref_result, mod_result, settings)

        chr1_file = next(
            f for f in result["metadata"]["contig_files"]
            if "_contig_0" in f
        )
        written_seq = str(SeqIO.read(chr1_file, "fasta").seq)
        assert written_seq == known_seq
