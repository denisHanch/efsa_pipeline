"""
Tests for genome characterization (contigs vs plasmids via minimap2),
which is now integrated into genomexgenome_validation.

Tests cover:
- characterize=True/False on GenomeXGenomeSettings
- Successful run with mixed contigs and plasmids
- Edge cases: all mapped, none mapped
- minimap2 unavailability (demoted to warning, validation still passes)
- Correct sequence content in output files
- PAF best-hit selection (multi-hit queries, ribosomal RNA scenario)
- Orientation tracking ('+' and '-' strand)
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
    _parse_paf_best_hits,
    _filter_and_group_hits,
)
from validation_pkg.validators.genome_validator import GenomeOutputMetadata
from validation_pkg.exceptions import InterFileValidationError


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _mock_minimap2(paf_content: str):
    """Return a side_effect for subprocess.run that writes PAF content to the stdout file."""
    def side_effect(cmd, stdout=None, **kwargs):
        if stdout is not None and hasattr(stdout, 'write'):
            stdout.write(paf_content)
        return MagicMock(returncode=0)
    return side_effect


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
    q_len: int = None,
    r_end: int = None,
    residue_matches: int = None,
    mapq: int = 255,
) -> str:
    """Build a minimal valid PAF line (12 tab-separated columns)."""
    q_len = q_len if q_len is not None else alignment_len
    r_end = r_end if r_end is not None else r_start + alignment_len
    residue_matches = residue_matches if residue_matches is not None else alignment_len
    return (
        f"{query_id}\t{q_len}\t{q_start}\t{q_start + alignment_len}\t{strand}\t"
        f"{ref_name}\t1000000\t{r_start}\t{r_end}\t{residue_matches}\t{alignment_len}\t{mapq}"
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
            mock_run.side_effect = _mock_minimap2(paf)
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
            mock_run.side_effect = _mock_minimap2(paf)
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
            mock_run.side_effect = _mock_minimap2("")
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
            mock_run.side_effect = _mock_minimap2(paf)
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
        assert "ref1" in hits["chr1"]
        assert hits["chr1"]["ref1"]["strand"] == "+"
        assert hits["chr1"]["ref1"]["alignment_len"] == 500

    def test_single_hit_reverse(self):
        paf = _make_paf_line("chr1", strand="-", ref_name="ref1", alignment_len=500)
        hits = _parse_paf_best_hits(paf)
        assert hits["chr1"]["ref1"]["strand"] == "-"

    def test_all_ref_hits_kept_per_query(self):
        """When a query maps to multiple refs, all (query, ref) pairs are kept."""
        lines = [
            _make_paf_line("chr1", strand="+", ref_name="ref_small", alignment_len=200),
            _make_paf_line("chr1", strand="-", ref_name="ref_large", alignment_len=5000),
            _make_paf_line("chr1", strand="+", ref_name="ref_mid",   alignment_len=1000),
        ]
        hits = _parse_paf_best_hits("\n".join(lines))
        # All three refs are retained as separate entries
        assert set(hits["chr1"].keys()) == {"ref_small", "ref_large", "ref_mid"}
        assert hits["chr1"]["ref_small"]["alignment_len"] == 200
        assert hits["chr1"]["ref_large"]["alignment_len"] == 5000
        assert hits["chr1"]["ref_large"]["strand"] == "-"
        assert hits["chr1"]["ref_mid"]["alignment_len"] == 1000

    def test_ribosomal_rna_all_ref_hits_kept(self):
        """Multiple ribosomal hits: all unique ref names produce separate entries."""
        rrna_hits = [
            _make_paf_line("rRNA_contig", strand="+", ref_name=f"rrna_copy_{i}", alignment_len=5900 + i)
            for i in range(5)
        ]
        hits = _parse_paf_best_hits("\n".join(rrna_hits))
        assert set(hits["rRNA_contig"].keys()) == {f"rrna_copy_{i}" for i in range(5)}
        assert hits["rRNA_contig"]["rrna_copy_4"]["alignment_len"] == 5904

    def test_multiple_queries_independent(self):
        lines = [
            _make_paf_line("chr1", strand="+", ref_name="ref1", alignment_len=1000),
            _make_paf_line("chr2", strand="-", ref_name="ref2", alignment_len=2000),
        ]
        hits = _parse_paf_best_hits("\n".join(lines))
        assert hits["chr1"]["ref1"]["strand"] == "+"
        assert hits["chr2"]["ref2"]["strand"] == "-"

    def test_empty_output(self):
        assert _parse_paf_best_hits("") == {}

    def test_short_lines_skipped(self):
        """Lines with fewer than 12 columns are ignored."""
        assert _parse_paf_best_hits("chr1\t100\t0\t100\t+\tref1") == {}

    def test_duplicate_pair_best_alignment_kept(self):
        """Same (query, ref) pair: only the hit with the larger alignment_len is kept."""
        lines = [
            _make_paf_line("chr1", ref_name="ref1", alignment_len=100, q_start=0, r_start=0),
            _make_paf_line("chr1", ref_name="ref1", alignment_len=200, q_start=0, r_start=0),
        ]
        hits = _parse_paf_best_hits("\n".join(lines))
        assert hits["chr1"]["ref1"]["alignment_len"] == 200

    def test_different_ref_names_both_kept(self):
        """Same query mapping to two different refs: both entries are retained."""
        lines = [
            _make_paf_line("chr1", ref_name="ref1", q_start=0,   r_start=0),
            _make_paf_line("chr1", ref_name="ref2", q_start=500, r_start=999),
        ]
        hits = _parse_paf_best_hits("\n".join(lines))
        assert "ref1" in hits["chr1"]
        assert "ref2" in hits["chr1"]


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
            _make_paf_line("chr1", strand="+", ref_name="ref_chr1"),
            _make_paf_line("chr2", strand="-", ref_name="ref_chr2"),
        ])
        with patch("validation_pkg.validators.interfile_genome.check_tool_available", return_value=True), \
             patch("validation_pkg.validators.interfile_genome.subprocess.run") as mock_run:
            mock_run.side_effect = _mock_minimap2(paf)
            result = genomexgenome_validation(ref_result, mod_result, settings)

        orientations = result["metadata"]["contig_orientations"]
        # Individual contig entries are keyed by (query_name, ref_name) tuple
        assert orientations[("chr1", "ref_chr1")] == "+"
        assert orientations[("chr2", "ref_chr2")] == "-"

    def test_ref_names_stored_in_metadata(self, fasta_setup):
        tmp_path, ref_result, mod_result, settings = fasta_setup
        paf = "\n".join([
            _make_paf_line("chr1", ref_name="ref_chr1"),
            _make_paf_line("chr2", ref_name="ref_chr2"),
        ])
        with patch("validation_pkg.validators.interfile_genome.check_tool_available", return_value=True), \
             patch("validation_pkg.validators.interfile_genome.subprocess.run") as mock_run:
            mock_run.side_effect = _mock_minimap2(paf)
            result = genomexgenome_validation(ref_result, mod_result, settings)

        ref_names = result["metadata"]["contig_ref_names"]
        # Individual contig entries are keyed by (query_name, ref_name) tuple
        assert ref_names[("chr1", "ref_chr1")] == "ref_chr1"
        assert ref_names[("chr2", "ref_chr2")] == "ref_chr2"

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
            mock_run.side_effect = _mock_minimap2(paf)
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
            mock_run.side_effect = _mock_minimap2(paf)
            result = genomexgenome_validation(ref_result, mod_result, settings)

        chr1_file = next(
            f for f in result["metadata"]["contig_files"]
            if "_contig_0" in f
        )
        written_seq = str(SeqIO.read(chr1_file, "fasta").seq)
        assert written_seq == known_seq



# ---------------------------------------------------------------------------
# Ref plasmid files included in minimap2 command
# ---------------------------------------------------------------------------

class TestRefPlasmidFilesInMinimap2:
    """Ref plasmid files (from plasmids_to_one=True processing) must be
    passed to minimap2 so mod sequences homologous to non-main ref
    chromosomes are not mis-classified as plasmids."""

    @pytest.fixture
    def setup(self, tmp_path):
        ref_path = tmp_path / "ref.fasta"
        plasmid_path = tmp_path / "ref_plasmid.fasta"
        mod_path = tmp_path / "mod.fasta"

        _write_fasta(ref_path, _make_records(["ref_chr"]))
        _write_fasta(plasmid_path, _make_records(["ref_plasmid"]))
        _write_fasta(mod_path, _make_records(["chr", "chr1", "chr2"]))

        ref_result = GenomeOutputMetadata(
            output_file=str(ref_path),
            output_filename="ref.fasta",
            num_sequences=1,
            plasmid_filenames=["ref_plasmid.fasta"],
        )
        mod_result = GenomeOutputMetadata(
            output_file=str(mod_path),
            output_filename="mod.fasta",
            num_sequences=3,
        )
        settings = GenomeXGenomeSettings(same_number_of_sequences=False, characterize=True)
        return tmp_path, ref_result, mod_result, settings

    def test_plasmid_file_appended_to_minimap2_command(self, setup):
        """minimap2 receives both the main ref file and the plasmid file."""
        tmp_path, ref_result, mod_result, settings = setup
        paf = "\n".join([_make_paf_line("chr"), _make_paf_line("chr1"), _make_paf_line("chr2")])

        with patch("validation_pkg.validators.interfile_genome.check_tool_available", return_value=True), \
             patch("validation_pkg.validators.interfile_genome.subprocess.run") as mock_run:
            mock_run.side_effect = _mock_minimap2(paf)
            genomexgenome_validation(ref_result, mod_result, settings)

        cmd = mock_run.call_args[0][0]
        assert str(tmp_path / "ref.fasta") in cmd
        assert str(tmp_path / "ref_plasmid.fasta") in cmd

    def test_missing_plasmid_file_skipped_gracefully(self, setup):
        """A plasmid filename that doesn't exist on disk is silently skipped."""
        tmp_path, ref_result, mod_result, settings = setup
        ref_result.plasmid_filenames = ["nonexistent_plasmid.fasta"]
        paf = _make_paf_line("chr")

        with patch("validation_pkg.validators.interfile_genome.check_tool_available", return_value=True), \
             patch("validation_pkg.validators.interfile_genome.subprocess.run") as mock_run:
            mock_run.side_effect = _mock_minimap2(paf)
            result = genomexgenome_validation(ref_result, mod_result, settings)

        cmd = mock_run.call_args[0][0]
        assert "nonexistent_plasmid.fasta" not in " ".join(cmd)
        assert result["passed"] is True

    def test_no_plasmid_filenames_uses_only_main_ref(self, setup):
        """When ref has no plasmid files, only the main ref is passed to minimap2."""
        tmp_path, ref_result, mod_result, settings = setup
        ref_result.plasmid_filenames = []
        paf = _make_paf_line("chr")

        with patch("validation_pkg.validators.interfile_genome.check_tool_available", return_value=True), \
             patch("validation_pkg.validators.interfile_genome.subprocess.run") as mock_run:
            mock_run.side_effect = _mock_minimap2(paf)
            genomexgenome_validation(ref_result, mod_result, settings)

        cmd = mock_run.call_args[0][0]
        ref_files_in_cmd = [c for c in cmd if c.endswith(".fasta") and "mod" not in c]
        assert ref_files_in_cmd == [str(tmp_path / "ref.fasta")]


# ---------------------------------------------------------------------------
# Settings tests for concatenation
# ---------------------------------------------------------------------------

class TestGenomeXGenomeSettingsConcatenate:

    def test_concatenate_default_false(self):
        assert GenomeXGenomeSettings().concatenate is False

    def test_concatenate_enabled(self):
        assert GenomeXGenomeSettings(concatenate=True).concatenate is True

    def test_min_query_coverage_invalid_raises(self):
        with pytest.raises(Exception, match="min_query_coverage"):
            GenomeXGenomeSettings(min_query_coverage=1.5)

    def test_min_query_coverage_negative_raises(self):
        with pytest.raises(Exception, match="min_query_coverage"):
            GenomeXGenomeSettings(min_query_coverage=-0.1)

    def test_min_identity_invalid_raises(self):
        with pytest.raises(Exception, match="min_identity"):
            GenomeXGenomeSettings(min_identity=2.0)

    def test_min_mapping_quality_negative_raises(self):
        with pytest.raises(Exception, match="min_mapping_quality"):
            GenomeXGenomeSettings(min_mapping_quality=-1)

    def test_max_ref_overlap_negative_raises(self):
        with pytest.raises(Exception, match="max_ref_overlap"):
            GenomeXGenomeSettings(max_ref_overlap=-1)

    def test_settings_roundtrip(self):
        s = GenomeXGenomeSettings(concatenate=True, min_query_coverage=0.8, min_identity=0.95)
        d = s.to_dict()
        s2 = GenomeXGenomeSettings.from_dict(d)
        assert s2.concatenate is True
        assert s2.min_query_coverage == 0.8
        assert s2.min_identity == 0.95


# ---------------------------------------------------------------------------
# Unit tests for _filter_and_group_hits
# ---------------------------------------------------------------------------

class TestFilterAndGroupHits:

    def _make_logger(self):
        import logging
        return logging.getLogger("test")

    def _hit(self, strand="+", alignment_len=1000,
              q_len=1000, r_start=0, r_end=1000, residue_matches=1000, mapq=60):
        """Return a single hit dict (ref_name is now the key, not a field)."""
        return {
            'strand': strand,
            'alignment_len': alignment_len, 'q_len': q_len,
            'q_end': alignment_len, 'r_start': r_start, 'r_end': r_end,
            'residue_matches': residue_matches, 'mapq': mapq,
        }

    def test_single_hit_passes(self):
        hits = {"ctg1": {"chrX": self._hit()}}
        settings = GenomeXGenomeSettings(concatenate=True)
        groups, individual_ids = _filter_and_group_hits(hits, settings, self._make_logger())
        assert "chrX" in groups
        assert groups["chrX"] == ["ctg1"]
        assert individual_ids == set()

    def test_mapq_below_threshold_goes_individual(self):
        hits = {"ctg1": {"chrX": self._hit(mapq=5)}}
        settings = GenomeXGenomeSettings(concatenate=True, min_mapping_quality=10)
        groups, individual_ids = _filter_and_group_hits(hits, settings, self._make_logger())
        assert groups == {}
        assert ("ctg1", "chrX") in individual_ids

    def test_query_coverage_below_threshold_goes_individual(self):
        # alignment_len=400, q_len=1000 → coverage=0.4 < 0.5
        hits = {"ctg1": {"chrX": self._hit(alignment_len=400, q_len=1000)}}
        settings = GenomeXGenomeSettings(concatenate=True, min_query_coverage=0.5)
        groups, individual_ids = _filter_and_group_hits(hits, settings, self._make_logger())
        assert groups == {}
        assert ("ctg1", "chrX") in individual_ids

    def test_identity_below_threshold_goes_individual(self):
        # residue_matches=800, alignment_len=1000 → identity=0.8 < 0.9
        hits = {"ctg1": {"chrX": self._hit(residue_matches=800, alignment_len=1000)}}
        settings = GenomeXGenomeSettings(concatenate=True, min_identity=0.9)
        groups, individual_ids = _filter_and_group_hits(hits, settings, self._make_logger())
        assert groups == {}
        assert ("ctg1", "chrX") in individual_ids

    def test_multiple_queries_same_ref_grouped_sorted_by_r_start(self):
        hits = {
            "ctg_b": {"chr1": self._hit(r_start=5000, r_end=9000)},
            "ctg_a": {"chr1": self._hit(r_start=0,    r_end=4000)},
            "ctg_c": {"chr1": self._hit(r_start=9500, r_end=13000)},
        }
        settings = GenomeXGenomeSettings(concatenate=True)
        groups, individual_ids = _filter_and_group_hits(hits, settings, self._make_logger())
        assert groups["chr1"] == ["ctg_a", "ctg_b", "ctg_c"]
        assert individual_ids == set()

    def test_different_refs_form_separate_groups(self):
        hits = {
            "ctg1": {"chr1": self._hit()},
            "ctg2": {"chr2": self._hit()},
        }
        settings = GenomeXGenomeSettings(concatenate=True)
        groups, individual_ids = _filter_and_group_hits(hits, settings, self._make_logger())
        assert set(groups.keys()) == {"chr1", "chr2"}
        assert individual_ids == set()

    def test_overlap_exceeds_max_demotes_entire_group(self):
        # ctg_a ends at 1000, ctg_b starts at 800 → overlap=200 > max_ref_overlap=0
        hits = {
            "ctg_a": {"chr1": self._hit(r_start=0,   r_end=1000)},
            "ctg_b": {"chr1": self._hit(r_start=800, r_end=2000)},
        }
        settings = GenomeXGenomeSettings(concatenate=True, max_ref_overlap=0)
        groups, individual_ids = _filter_and_group_hits(hits, settings, self._make_logger())
        assert groups == {}
        assert {("ctg_a", "chr1"), ("ctg_b", "chr1")} == individual_ids

    def test_overlap_within_max_passes(self):
        # overlap=100, max_ref_overlap=200 → passes
        hits = {
            "ctg_a": {"chr1": self._hit(r_start=0,   r_end=1000)},
            "ctg_b": {"chr1": self._hit(r_start=900, r_end=2000)},
        }
        settings = GenomeXGenomeSettings(concatenate=True, max_ref_overlap=200)
        groups, individual_ids = _filter_and_group_hits(hits, settings, self._make_logger())
        assert "chr1" in groups
        assert individual_ids == set()

    def test_overlap_in_one_group_does_not_affect_other(self):
        hits = {
            "bad_a": {"chr1": self._hit(r_start=0,   r_end=1000)},
            "bad_b": {"chr1": self._hit(r_start=500, r_end=2000)},
            "good":  {"chr2": self._hit(r_start=0,   r_end=1000)},
        }
        settings = GenomeXGenomeSettings(concatenate=True, max_ref_overlap=0)
        groups, individual_ids = _filter_and_group_hits(hits, settings, self._make_logger())
        assert "chr1" not in groups
        assert "chr2" in groups
        assert {("bad_a", "chr1"), ("bad_b", "chr1")} == individual_ids


# ---------------------------------------------------------------------------
# Integration tests for concatenation mode
# ---------------------------------------------------------------------------

class TestConcatenateMode:

    @pytest.fixture
    def setup(self, tmp_path):
        # 4 mod sequences: ctg_a + ctg_b → chr1, ctg_c → chr2, plasmid stays unmapped
        ref_seqs = _make_records(["chr1", "chr2"])
        mod_seqs = [
            SeqRecord(Seq("AAAA" * 25), id="ctg_a", description=""),
            SeqRecord(Seq("CCCC" * 25), id="ctg_b", description=""),
            SeqRecord(Seq("GGGG" * 25), id="ctg_c", description=""),
            SeqRecord(Seq("TTTT" * 25), id="plasmid1", description=""),
        ]
        ref_path = tmp_path / "ref.fasta"
        mod_path = tmp_path / "mod.fasta"
        _write_fasta(ref_path, ref_seqs)
        _write_fasta(mod_path, mod_seqs)

        ref_result = GenomeOutputMetadata(
            output_file=str(ref_path), output_filename="ref.fasta", num_sequences=2,
        )
        mod_result = GenomeOutputMetadata(
            output_file=str(mod_path), output_filename="mod.fasta", num_sequences=4,
        )
        return tmp_path, ref_result, mod_result

    def _run(self, ref_result, mod_result, paf, **setting_kwargs):
        settings = GenomeXGenomeSettings(
            same_number_of_sequences=False, characterize=True,
            concatenate=True, **setting_kwargs
        )
        with patch("validation_pkg.validators.interfile_genome.check_tool_available", return_value=True), \
             patch("validation_pkg.validators.interfile_genome.subprocess.run") as mock_run:
            mock_run.side_effect = _mock_minimap2(paf)
            return genomexgenome_validation(ref_result, mod_result, settings)

    def test_two_queries_same_ref_produce_one_contig_file(self, setup):
        tmp_path, ref_result, mod_result = setup
        paf = "\n".join([
            _make_paf_line("ctg_a", ref_name="chr1", r_start=0,    r_end=100, alignment_len=100),
            _make_paf_line("ctg_b", ref_name="chr1", r_start=100,  r_end=200, alignment_len=100),
            _make_paf_line("ctg_c", ref_name="chr2", r_start=0,    r_end=100, alignment_len=100),
        ])
        result = self._run(ref_result, mod_result, paf)
        meta = result["metadata"]
        # chr1 group → 1 file, chr2 individual → 1 file = 2 contig files total
        assert len(meta["contig_files"]) == 2
        assert meta["plasmids_found"] == 1
        assert meta["concat_groups"] is not None
        assert "chr1" in meta["concat_groups"]
        assert meta["concat_groups"]["chr1"]["query_names"] == ["ctg_a", "ctg_b"]

    def test_concat_sequence_is_in_r_start_order(self, setup):
        tmp_path, ref_result, mod_result = setup
        # ctg_b maps earlier on ref (r_start=0), ctg_a maps later (r_start=100)
        paf = "\n".join([
            _make_paf_line("ctg_a", ref_name="chr1", r_start=100, r_end=200, alignment_len=100),
            _make_paf_line("ctg_b", ref_name="chr1", r_start=0,   r_end=100, alignment_len=100),
            _make_paf_line("ctg_c", ref_name="chr2", r_start=0,   r_end=100, alignment_len=100),
        ])
        result = self._run(ref_result, mod_result, paf)
        meta = result["metadata"]
        # ctg_b (r_start=0) should come before ctg_a (r_start=100)
        assert meta["concat_groups"]["chr1"]["query_names"] == ["ctg_b", "ctg_a"]

    def test_concat_file_contains_merged_sequence(self, setup):
        tmp_path, ref_result, mod_result = setup
        paf = "\n".join([
            _make_paf_line("ctg_a", ref_name="chr1", r_start=0,   r_end=100, alignment_len=100),
            _make_paf_line("ctg_b", ref_name="chr1", r_start=100, r_end=200, alignment_len=100),
            _make_paf_line("ctg_c", ref_name="chr2", r_start=0,   r_end=100, alignment_len=100),
        ])
        result = self._run(ref_result, mod_result, paf)
        chr1_file = result["metadata"]["concat_groups"]["chr1"]["contig_file"]
        records = list(SeqIO.parse(chr1_file, "fasta"))
        assert len(records) == 1
        # ctg_a = AAAA*25 (100bp), ctg_b = CCCC*25 (100bp) → 200bp total
        assert len(records[0].seq) == 200

    def test_query_failing_threshold_written_individually(self, setup):
        tmp_path, ref_result, mod_result = setup
        # ctg_a has low MAPQ → falls back to individual file; ctg_b alone on chr1
        paf = "\n".join([
            _make_paf_line("ctg_a", ref_name="chr1", r_start=0,   r_end=100,
                           alignment_len=100, mapq=2),
            _make_paf_line("ctg_b", ref_name="chr1", r_start=100, r_end=200,
                           alignment_len=100, mapq=60),
            _make_paf_line("ctg_c", ref_name="chr2", r_start=0,   r_end=100,
                           alignment_len=100, mapq=60),
        ])
        result = self._run(ref_result, mod_result, paf, min_mapping_quality=10)
        meta = result["metadata"]
        # chr1 group has only ctg_b (ctg_a failed), chr2 has ctg_c → 2 concat groups
        # ctg_a written individually → 1 extra file
        # total = 2 + 1 = 3 contig files
        assert len(meta["contig_files"]) == 3

    def test_overlapping_group_written_individually(self, setup):
        tmp_path, ref_result, mod_result = setup
        # ctg_a and ctg_b overlap on chr1 → both written individually
        paf = "\n".join([
            _make_paf_line("ctg_a", ref_name="chr1", r_start=0,   r_end=150, alignment_len=100),
            _make_paf_line("ctg_b", ref_name="chr1", r_start=100, r_end=200, alignment_len=100),
            _make_paf_line("ctg_c", ref_name="chr2", r_start=0,   r_end=100, alignment_len=100),
        ])
        result = self._run(ref_result, mod_result, paf, max_ref_overlap=0)
        meta = result["metadata"]
        assert "chr1" not in (meta["concat_groups"] or {})
        # ctg_a + ctg_b individual + ctg_c individual = 3 files
        assert len(meta["contig_files"]) == 3

    def test_contig_orientations_nested_for_concat_groups(self, setup):
        tmp_path, ref_result, mod_result = setup
        paf = "\n".join([
            _make_paf_line("ctg_a", strand="+", ref_name="chr1", r_start=0,   r_end=100, alignment_len=100),
            _make_paf_line("ctg_b", strand="-", ref_name="chr1", r_start=100, r_end=200, alignment_len=100),
            _make_paf_line("ctg_c", ref_name="chr2", r_start=0,  r_end=100, alignment_len=100),
        ])
        result = self._run(ref_result, mod_result, paf)
        orientations = result["metadata"]["contig_orientations"]
        # chr1 is a concat group → nested dict
        assert isinstance(orientations["chr1"], dict)
        assert orientations["chr1"]["ctg_a"] == "+"
        assert orientations["chr1"]["ctg_b"] == "-"

    def test_concatenate_false_behaves_as_before(self, setup):
        tmp_path, ref_result, mod_result = setup
        paf = "\n".join([
            _make_paf_line("ctg_a", ref_name="chr1"),
            _make_paf_line("ctg_b", ref_name="chr1"),
            _make_paf_line("ctg_c", ref_name="chr2"),
        ])
        settings = GenomeXGenomeSettings(same_number_of_sequences=False, characterize=True)
        with patch("validation_pkg.validators.interfile_genome.check_tool_available", return_value=True), \
             patch("validation_pkg.validators.interfile_genome.subprocess.run") as mock_run:
            mock_run.side_effect = _mock_minimap2(paf)
            result = genomexgenome_validation(ref_result, mod_result, settings)
        meta = result["metadata"]
        assert len(meta["contig_files"]) == 3
        assert meta["concat_groups"] is None
        # In non-concat mode all entries are individual: keys are (query, ref) tuples, values are strand strings
        assert all(isinstance(k, tuple) for k in meta["contig_orientations"].keys())
        assert all(isinstance(v, str) for v in meta["contig_orientations"].values())

    def test_query_mapping_two_refs_produces_two_contig_files(self, setup):
        """A single query that maps to two distinct refs gets a unique contig file per ref."""
        tmp_path, ref_result, mod_result = setup
        paf = "\n".join([
            # ctg_a maps to BOTH chr1 and chr2
            _make_paf_line("ctg_a", ref_name="chr1", r_start=0,   r_end=100, alignment_len=100),
            _make_paf_line("ctg_a", ref_name="chr2", r_start=0,   r_end=100, alignment_len=80),
            # ctg_b maps only to chr1
            _make_paf_line("ctg_b", ref_name="chr1", r_start=100, r_end=200, alignment_len=100),
        ])
        settings = GenomeXGenomeSettings(same_number_of_sequences=False, characterize=True)
        with patch("validation_pkg.validators.interfile_genome.check_tool_available", return_value=True), \
             patch("validation_pkg.validators.interfile_genome.subprocess.run") as mock_run:
            mock_run.side_effect = _mock_minimap2(paf)
            result = genomexgenome_validation(ref_result, mod_result, settings)

        meta = result["metadata"]
        # ctg_a → chr1, ctg_a → chr2, ctg_b → chr1 = 3 individual contig files
        assert len(meta["contig_files"]) == 3
        # Both (ctg_a, chr1) and (ctg_a, chr2) are in orientations
        assert ("ctg_a", "chr1") in meta["contig_orientations"]
        assert ("ctg_a", "chr2") in meta["contig_orientations"]
        assert ("ctg_b", "chr1") in meta["contig_orientations"]
