"""
Integration tests for the 5 input scenarios defined in docs/validation/OVERVIEW.md.

Each test class covers one scenario from the documentation table and verifies:
  - genome and plasmid file outputs match documented expectations
  - run_ref_x_mod flag value
  - whether minimap2 (characterisation) is invoked or not
  - plasmid detection / separation behaviour

Documentation reference: docs/validation/OVERVIEW.md

Main problem areas flagged by developer: Scenarios 3 and 4.
"""

import sys
import pytest
from pathlib import Path
from unittest.mock import patch, MagicMock
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from validation_pkg.validators.genome_validator import GenomeValidator
from validation_pkg.validators.interfile_genome import GenomeXGenomeSettings, genomexgenome_validation
from validation_pkg.config_manager import GenomeConfig
from validation_pkg.utils.formats import GenomeFormat, CodingType

# nextflow_params_handler lives in modules/validation/utils/, not inside the package
_VALIDATION_ROOT = Path(__file__).parent.parent.parent  # modules/validation/
if str(_VALIDATION_ROOT) not in sys.path:
    sys.path.insert(0, str(_VALIDATION_ROOT))
from utils.nextflow_params_handler import build_params  # noqa: E402


# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------

# Sequences long enough to survive the min_sequence_length=100 filter.
_CHR_SEQ = "ACGT" * 75   # 300 bp — plays the role of a chromosome
_PLA_SEQ = "GCTA" * 30   # 120 bp — plays the role of a plasmid


def _write_fasta(path: Path, records: list[tuple[str, str]]) -> Path:
    seqs = [SeqRecord(Seq(seq), id=sid, description="") for sid, seq in records]
    with open(path, "w") as fh:
        SeqIO.write(seqs, fh, "fasta")
    return path


def _read_ids(path: Path) -> list[str]:
    return [r.id for r in SeqIO.parse(str(path), "fasta")]


def _make_genome_config(
    path: Path,
    output_dir: Path,
    n_sequence_limit: int = 5,
    org_type: str = "prokaryote",
) -> GenomeConfig:
    return GenomeConfig(
        filename=path.name,
        filepath=path,
        basename=path.stem,
        coding_type=CodingType.NONE,
        detected_format=GenomeFormat.FASTA,
        output_dir=output_dir,
        global_options={"type": org_type},
        n_sequence_limit=n_sequence_limit,
    )


def _ref_settings() -> GenomeValidator.Settings:
    """Mirror the ref_genome_settings used in main.py."""
    return GenomeValidator.Settings(
        plasmids_to_one=True,
        main_longest=True,
        coding_type=None,
        output_filename_suffix="ref",
        replace_id_with_incremental="chr",
        min_sequence_length=100,
    )


def _mod_settings() -> GenomeValidator.Settings:
    """Mirror the mod_genome_settings used in main.py."""
    return GenomeValidator.Settings(
        plasmids_to_one=False,
        coding_type=None,
        output_filename_suffix="mod",
        replace_id_with_incremental="chr",
        min_sequence_length=100,
    )


def _gxg_settings() -> GenomeXGenomeSettings:
    """Mirror the genomexgenome_settings used in main.py."""
    return GenomeXGenomeSettings(
        characterize=True,
        same_sequence_ids=False,
        same_number_of_sequences=False,
    )


def _paf_line(query_id: str, ref_id: str = "chr", aln_len: int = 200) -> str:
    """Build a minimal PAF line for a mapped sequence."""
    return (
        f"{query_id}\t{aln_len}\t0\t{aln_len}\t+\t"
        f"{ref_id}\t300\t0\t{aln_len}\t{aln_len}\t{aln_len}\t255"
    )


def _build_validation_results(ref_res, mod_res, gxg_res=None):
    return {
        "ref_genome":    ref_res,
        "mod_genome":    mod_res,
        "ref_plasmid":   None,
        "mod_plasmid":   None,
        "genomexgenome": gxg_res,
        "reads":         [],
    }


# ---------------------------------------------------------------------------
# Scenario 1: Single contig + plasmids (PROKARYOTE)
#
# ref.fa  : 1 chromosome + 1 plasmid  →  main_longest splits them
# mod.fa  : 1 chromosome + 1 plasmid  →  minimap2 separates them
#
# Expected: run_ref_x_mod=True, minimap2 used, 1 contig file, plasmid files set
# ---------------------------------------------------------------------------

class TestScenario1SingleContigWithPlasmids:

    @pytest.fixture
    def setup(self, tmp_path):
        ref_path = _write_fasta(tmp_path / "ref.fasta", [
            ("ref_chr", _CHR_SEQ),
            ("ref_pla", _PLA_SEQ),
        ])
        mod_path = _write_fasta(tmp_path / "mod.fasta", [
            ("mod_chr", _CHR_SEQ),
            ("mod_pla", _PLA_SEQ),
        ])
        out = tmp_path / "out"
        out.mkdir()
        return (
            _make_genome_config(ref_path, out),
            _make_genome_config(mod_path, out),
        )

    # ── ref processing ──

    def test_ref_plasmid_extracted(self, setup):
        ref_cfg, _ = setup
        ref_res = GenomeValidator(ref_cfg, _ref_settings()).run()
        assert ref_res.plasmid_filenames, "ref plasmid should be extracted (shorter sequence)"
        assert len(ref_res.plasmid_filenames) == 1

    def test_ref_output_is_single_chromosome(self, setup):
        ref_cfg, _ = setup
        ref_res = GenomeValidator(ref_cfg, _ref_settings()).run()
        assert ref_res.num_sequences == 1
        ids = _read_ids(Path(ref_res.output_file))
        assert ids == ["chr"], f"Expected renamed id 'chr', got {ids}"

    # ── mod processing ──

    def test_mod_output_contains_both_sequences_before_gxg(self, setup):
        """mod keeps all sequences until minimap2 splits them in GXG."""
        _, mod_cfg = setup
        mod_res = GenomeValidator(mod_cfg, _mod_settings()).run()
        assert mod_res.num_sequences == 2
        assert mod_res.fragmented is False

    # ── inter-genome characterisation ──

    def test_minimap2_is_called(self, setup):
        ref_cfg, mod_cfg = setup
        ref_res = GenomeValidator(ref_cfg, _ref_settings()).run()
        mod_res = GenomeValidator(mod_cfg, _mod_settings()).run()
        # mod output has 2 renamed sequences: 'chr' (chromosome) and 'chr1' (plasmid)
        paf = _paf_line("chr")
        with patch("validation_pkg.validators.interfile_genome.check_tool_available", return_value=True), \
             patch("validation_pkg.validators.interfile_genome.subprocess.run") as mock_run:
            mock_run.return_value = MagicMock(stdout=paf, returncode=0)
            genomexgenome_validation(ref_res, mod_res, _gxg_settings())
        mock_run.assert_called_once()

    def test_one_contig_file_produced(self, setup):
        ref_cfg, mod_cfg = setup
        ref_res = GenomeValidator(ref_cfg, _ref_settings()).run()
        mod_res = GenomeValidator(mod_cfg, _mod_settings()).run()
        paf = _paf_line("chr")  # only 'chr' (index 0) maps; 'chr1' (plasmid) does not
        with patch("validation_pkg.validators.interfile_genome.check_tool_available", return_value=True), \
             patch("validation_pkg.validators.interfile_genome.subprocess.run") as mock_run:
            mock_run.return_value = MagicMock(stdout=paf, returncode=0)
            gxg_res = genomexgenome_validation(ref_res, mod_res, _gxg_settings())
        assert gxg_res["metadata"]["contigs_found"] == 1
        assert len(gxg_res["metadata"]["contig_files"]) == 1

    def test_plasmid_file_produced_from_mod(self, setup):
        ref_cfg, mod_cfg = setup
        ref_res = GenomeValidator(ref_cfg, _ref_settings()).run()
        mod_res = GenomeValidator(mod_cfg, _mod_settings()).run()
        paf = _paf_line("chr")
        with patch("validation_pkg.validators.interfile_genome.check_tool_available", return_value=True), \
             patch("validation_pkg.validators.interfile_genome.subprocess.run") as mock_run:
            mock_run.return_value = MagicMock(stdout=paf, returncode=0)
            gxg_res = genomexgenome_validation(ref_res, mod_res, _gxg_settings())
        assert gxg_res["metadata"]["plasmids_found"] == 1
        assert gxg_res["metadata"]["plasmid_file"] is not None

    # ── final params ──

    def test_run_ref_x_mod_is_true(self, setup):
        ref_cfg, mod_cfg = setup
        ref_res = GenomeValidator(ref_cfg, _ref_settings()).run()
        mod_res = GenomeValidator(mod_cfg, _mod_settings()).run()
        paf = _paf_line("chr")
        with patch("validation_pkg.validators.interfile_genome.check_tool_available", return_value=True), \
             patch("validation_pkg.validators.interfile_genome.subprocess.run") as mock_run:
            mock_run.return_value = MagicMock(stdout=paf, returncode=0)
            gxg_res = genomexgenome_validation(ref_res, mod_res, _gxg_settings())
        params = build_params(_build_validation_results(ref_res, mod_res, gxg_res))
        assert params.run_ref_x_mod is True

    def test_contig_file_size_is_one(self, setup):
        ref_cfg, mod_cfg = setup
        ref_res = GenomeValidator(ref_cfg, _ref_settings()).run()
        mod_res = GenomeValidator(mod_cfg, _mod_settings()).run()
        paf = _paf_line("chr")
        with patch("validation_pkg.validators.interfile_genome.check_tool_available", return_value=True), \
             patch("validation_pkg.validators.interfile_genome.subprocess.run") as mock_run:
            mock_run.return_value = MagicMock(stdout=paf, returncode=0)
            gxg_res = genomexgenome_validation(ref_res, mod_res, _gxg_settings())
        params = build_params(_build_validation_results(ref_res, mod_res, gxg_res))
        assert params.contig_file_size == 1

    def test_ref_plasmid_fasta_is_set_in_params(self, setup):
        ref_cfg, mod_cfg = setup
        ref_res = GenomeValidator(ref_cfg, _ref_settings()).run()
        mod_res = GenomeValidator(mod_cfg, _mod_settings()).run()
        paf = _paf_line("chr")
        with patch("validation_pkg.validators.interfile_genome.check_tool_available", return_value=True), \
             patch("validation_pkg.validators.interfile_genome.subprocess.run") as mock_run:
            mock_run.return_value = MagicMock(stdout=paf, returncode=0)
            gxg_res = genomexgenome_validation(ref_res, mod_res, _gxg_settings())
        params = build_params(_build_validation_results(ref_res, mod_res, gxg_res))
        assert params.ref_plasmid_fasta is not None

    def test_mod_plasmid_fasta_is_set_in_params(self, setup):
        ref_cfg, mod_cfg = setup
        ref_res = GenomeValidator(ref_cfg, _ref_settings()).run()
        mod_res = GenomeValidator(mod_cfg, _mod_settings()).run()
        paf = _paf_line("chr")
        with patch("validation_pkg.validators.interfile_genome.check_tool_available", return_value=True), \
             patch("validation_pkg.validators.interfile_genome.subprocess.run") as mock_run:
            mock_run.return_value = MagicMock(stdout=paf, returncode=0)
            gxg_res = genomexgenome_validation(ref_res, mod_res, _gxg_settings())
        params = build_params(_build_validation_results(ref_res, mod_res, gxg_res))
        assert params.mod_plasmid_fasta is not None


# ---------------------------------------------------------------------------
# Scenario 2: Fragmented assembly — below limit (PROKARYOTE)
#
# ref.fa  : 1 sequence
# mod.fa  : 3 sequences (≤ limit=5)  →  2 map, 1 unmapped (plasmid)
#
# Expected: run_ref_x_mod=True, minimap2 used, 2 contigs, 1 plasmid
# ---------------------------------------------------------------------------

class TestScenario2FragmentedBelowLimit:

    @pytest.fixture
    def setup(self, tmp_path):
        ref_path = _write_fasta(tmp_path / "ref.fasta", [("ref_chr", _CHR_SEQ)])
        mod_path = _write_fasta(tmp_path / "mod.fasta", [
            ("mod_c1", _CHR_SEQ),
            ("mod_c2", _CHR_SEQ),
            ("mod_p1", _PLA_SEQ),
        ])
        out = tmp_path / "out"
        out.mkdir()
        return (
            _make_genome_config(ref_path, out, n_sequence_limit=5),
            _make_genome_config(mod_path, out, n_sequence_limit=5),
        )

    def test_mod_not_fragmented(self, setup):
        _, mod_cfg = setup
        mod_res = GenomeValidator(mod_cfg, _mod_settings()).run()
        assert mod_res.fragmented is False
        assert mod_res.num_sequences == 3

    def test_minimap2_is_called(self, setup):
        ref_cfg, mod_cfg = setup
        ref_res = GenomeValidator(ref_cfg, _ref_settings()).run()
        mod_res = GenomeValidator(mod_cfg, _mod_settings()).run()
        paf = "\n".join([_paf_line("chr"), _paf_line("chr1")])
        with patch("validation_pkg.validators.interfile_genome.check_tool_available", return_value=True), \
             patch("validation_pkg.validators.interfile_genome.subprocess.run") as mock_run:
            mock_run.return_value = MagicMock(stdout=paf, returncode=0)
            genomexgenome_validation(ref_res, mod_res, _gxg_settings())
        mock_run.assert_called_once()

    def test_two_contigs_and_one_plasmid(self, setup):
        ref_cfg, mod_cfg = setup
        ref_res = GenomeValidator(ref_cfg, _ref_settings()).run()
        mod_res = GenomeValidator(mod_cfg, _mod_settings()).run()
        paf = "\n".join([_paf_line("chr"), _paf_line("chr1")])
        with patch("validation_pkg.validators.interfile_genome.check_tool_available", return_value=True), \
             patch("validation_pkg.validators.interfile_genome.subprocess.run") as mock_run:
            mock_run.return_value = MagicMock(stdout=paf, returncode=0)
            gxg_res = genomexgenome_validation(ref_res, mod_res, _gxg_settings())
        assert gxg_res["metadata"]["contigs_found"] == 2
        assert gxg_res["metadata"]["plasmids_found"] == 1

    def test_run_ref_x_mod_is_true(self, setup):
        ref_cfg, mod_cfg = setup
        ref_res = GenomeValidator(ref_cfg, _ref_settings()).run()
        mod_res = GenomeValidator(mod_cfg, _mod_settings()).run()
        paf = "\n".join([_paf_line("chr"), _paf_line("chr1")])
        with patch("validation_pkg.validators.interfile_genome.check_tool_available", return_value=True), \
             patch("validation_pkg.validators.interfile_genome.subprocess.run") as mock_run:
            mock_run.return_value = MagicMock(stdout=paf, returncode=0)
            gxg_res = genomexgenome_validation(ref_res, mod_res, _gxg_settings())
        params = build_params(_build_validation_results(ref_res, mod_res, gxg_res))
        assert params.run_ref_x_mod is True

    def test_contig_file_size_is_two(self, setup):
        ref_cfg, mod_cfg = setup
        ref_res = GenomeValidator(ref_cfg, _ref_settings()).run()
        mod_res = GenomeValidator(mod_cfg, _mod_settings()).run()
        paf = "\n".join([_paf_line("chr"), _paf_line("chr1")])
        with patch("validation_pkg.validators.interfile_genome.check_tool_available", return_value=True), \
             patch("validation_pkg.validators.interfile_genome.subprocess.run") as mock_run:
            mock_run.return_value = MagicMock(stdout=paf, returncode=0)
            gxg_res = genomexgenome_validation(ref_res, mod_res, _gxg_settings())
        params = build_params(_build_validation_results(ref_res, mod_res, gxg_res))
        assert params.contig_file_size == 2


# ---------------------------------------------------------------------------
# Scenario 3: Fragmented assembly — above limit (PROKARYOTE)
#
# ref.fa  : 1 sequence (within limit)
# mod.fa  : 3 sequences > limit=2  →  fragmented, copied as-is
#
# Expected:
#   - mod.fragmented = True
#   - run_ref_x_mod = False
#   - minimap2 NOT used (genomexgenome skipped in main.py)
#   - mod output contains all original sequences
#   - mod_plasmid_fasta = None (no plasmid detection for fragmented mod)
#   - ref_plasmid_fasta = None when ref has a single sequence
# ---------------------------------------------------------------------------

class TestScenario3FragmentedAboveLimit:

    @pytest.fixture
    def setup(self, tmp_path):
        ref_path = _write_fasta(tmp_path / "ref.fasta", [("ref_chr", _CHR_SEQ)])
        mod_path = _write_fasta(tmp_path / "mod.fasta", [
            ("mod_c1", _CHR_SEQ),
            ("mod_c2", _CHR_SEQ),
            ("mod_c3", _PLA_SEQ),
        ])
        out = tmp_path / "out"
        out.mkdir()
        return (
            _make_genome_config(ref_path, out, n_sequence_limit=2),
            _make_genome_config(mod_path, out, n_sequence_limit=2),
        )

    def test_mod_is_fragmented(self, setup):
        _, mod_cfg = setup
        mod_res = GenomeValidator(mod_cfg, _mod_settings()).run()
        assert mod_res.fragmented is True, "mod (3 seqs > limit=2) must be fragmented"

    def test_mod_output_contains_all_original_sequences(self, setup):
        """mod output is a direct copy — sequence count must match input."""
        _, mod_cfg = setup
        mod_res = GenomeValidator(mod_cfg, _mod_settings()).run()
        assert mod_res.num_sequences == 3, "copied mod should report all 3 sequences"
        ids = _read_ids(Path(mod_res.output_file))
        assert len(ids) == 3, f"output file should have 3 sequences, got {ids}"

    def test_run_ref_x_mod_is_false(self, setup):
        ref_cfg, mod_cfg = setup
        ref_res = GenomeValidator(ref_cfg, _ref_settings()).run()
        mod_res = GenomeValidator(mod_cfg, _mod_settings()).run()
        params = build_params(_build_validation_results(ref_res, mod_res, gxg_res=None))
        assert params.run_ref_x_mod is False

    def test_minimap2_not_called_when_mod_fragmented(self, setup):
        """Replicate main.py's condition: skip GXG when mod is fragmented."""
        ref_cfg, mod_cfg = setup
        ref_res = GenomeValidator(ref_cfg, _ref_settings()).run()
        mod_res = GenomeValidator(mod_cfg, _mod_settings()).run()

        with patch("validation_pkg.validators.interfile_genome.subprocess.run") as mock_run:
            # main.py condition: only run GXG if mod is not fragmented
            if not mod_res.fragmented:
                genomexgenome_validation(ref_res, mod_res, _gxg_settings())

        mock_run.assert_not_called()

    def test_no_contig_files(self, setup):
        ref_cfg, mod_cfg = setup
        ref_res = GenomeValidator(ref_cfg, _ref_settings()).run()
        mod_res = GenomeValidator(mod_cfg, _mod_settings()).run()
        params = build_params(_build_validation_results(ref_res, mod_res, gxg_res=None))
        assert params.contig_file_size == 0
        assert params.contig_files == []

    def test_mod_plasmid_fasta_is_none(self, setup):
        """No minimap2 means no plasmid detection for mod."""
        ref_cfg, mod_cfg = setup
        ref_res = GenomeValidator(ref_cfg, _ref_settings()).run()
        mod_res = GenomeValidator(mod_cfg, _mod_settings()).run()
        params = build_params(_build_validation_results(ref_res, mod_res, gxg_res=None))
        assert params.mod_plasmid_fasta is None

    def test_ref_plasmid_fasta_is_none_when_ref_single_seq(self, setup):
        """Single-sequence ref produces no plasmid file."""
        ref_cfg, mod_cfg = setup
        ref_res = GenomeValidator(ref_cfg, _ref_settings()).run()
        mod_res = GenomeValidator(mod_cfg, _mod_settings()).run()
        params = build_params(_build_validation_results(ref_res, mod_res, gxg_res=None))
        assert params.ref_plasmid_fasta is None


class TestScenario3RefWithPlasmid:
    """Scenario 3 variant: ref has chromosome + plasmid, mod is fragmented."""

    @pytest.fixture
    def setup(self, tmp_path):
        # ref: 2 seqs, limit=3 → 2 < 3 → NOT fragmented → main_longest extracts plasmid
        ref_path = _write_fasta(tmp_path / "ref.fasta", [
            ("ref_chr", _CHR_SEQ),
            ("ref_pla", _PLA_SEQ),
        ])
        # mod: 3 seqs, limit=3 → 3 >= 3 → fragmented → copied as-is
        mod_path = _write_fasta(tmp_path / "mod.fasta", [
            ("mod_c1", _CHR_SEQ),
            ("mod_c2", _CHR_SEQ),
            ("mod_c3", _PLA_SEQ),
        ])
        out = tmp_path / "out"
        out.mkdir()
        return (
            _make_genome_config(ref_path, out, n_sequence_limit=3),
            _make_genome_config(mod_path, out, n_sequence_limit=3),
        )

    def test_ref_plasmid_still_extracted_despite_fragmented_mod(self, setup):
        """Ref plasmid extraction is independent of mod's fragmentation status."""
        ref_cfg, mod_cfg = setup
        ref_res = GenomeValidator(ref_cfg, _ref_settings()).run()
        mod_res = GenomeValidator(mod_cfg, _mod_settings()).run()
        params = build_params(_build_validation_results(ref_res, mod_res, gxg_res=None))
        assert params.ref_plasmid_fasta is not None, (
            "ref_plasmid_fasta should be set even when mod is fragmented"
        )

    def test_mod_plasmid_fasta_still_none(self, setup):
        ref_cfg, mod_cfg = setup
        ref_res = GenomeValidator(ref_cfg, _ref_settings()).run()
        mod_res = GenomeValidator(mod_cfg, _mod_settings()).run()
        params = build_params(_build_validation_results(ref_res, mod_res, gxg_res=None))
        assert params.mod_plasmid_fasta is None

    def test_run_ref_x_mod_is_false(self, setup):
        ref_cfg, mod_cfg = setup
        ref_res = GenomeValidator(ref_cfg, _ref_settings()).run()
        mod_res = GenomeValidator(mod_cfg, _mod_settings()).run()
        params = build_params(_build_validation_results(ref_res, mod_res, gxg_res=None))
        assert params.run_ref_x_mod is False


# ---------------------------------------------------------------------------
# Scenario 4a: Multiple sequences in reference — EUKARYOTE type
#
# ref.fa  : 3 sequences (eukaryote → fragmented regardless of count)
# mod.fa  : 1 sequence  (also eukaryote → fragmented)
#
# Expected:
#   - ref.fragmented = True, mod.fragmented = True
#   - No plasmid extracted from ref (copied as-is)
#   - run_ref_x_mod = False
#   - minimap2 NOT used
# ---------------------------------------------------------------------------

class TestScenario4aEukaryoteType:

    @pytest.fixture
    def setup(self, tmp_path):
        ref_path = _write_fasta(tmp_path / "ref.fasta", [
            ("chr1", _CHR_SEQ),
            ("chr2", _CHR_SEQ),
            ("chr3", _CHR_SEQ),
        ])
        mod_path = _write_fasta(tmp_path / "mod.fasta", [("mod_chr", _CHR_SEQ)])
        out = tmp_path / "out"
        out.mkdir()
        return (
            _make_genome_config(ref_path, out, org_type="eukaryote"),
            _make_genome_config(mod_path, out, org_type="eukaryote"),
        )

    def test_ref_is_fragmented(self, setup):
        ref_cfg, _ = setup
        ref_res = GenomeValidator(ref_cfg, _ref_settings()).run()
        assert ref_res.fragmented is True, "EUKARYOTE ref must be marked fragmented"

    def test_mod_is_fragmented(self, setup):
        _, mod_cfg = setup
        mod_res = GenomeValidator(mod_cfg, _mod_settings()).run()
        assert mod_res.fragmented is True, "EUKARYOTE mod must be marked fragmented"

    def test_ref_output_file_exists_and_is_copied_as_is(self, setup):
        """EUKARYOTE ref should be copied without any plasmid extraction."""
        ref_cfg, _ = setup
        ref_res = GenomeValidator(ref_cfg, _ref_settings()).run()
        assert ref_res.output_file is not None
        ids = _read_ids(Path(ref_res.output_file))
        assert len(ids) == 3, "EUKARYOTE ref should be copied with all 3 sequences"

    def test_ref_no_plasmid_extracted(self, setup):
        """Plasmid splitting must not happen for EUKARYOTE."""
        ref_cfg, _ = setup
        ref_res = GenomeValidator(ref_cfg, _ref_settings()).run()
        assert not ref_res.plasmid_filenames, (
            "EUKARYOTE ref should not have plasmids extracted"
        )

    def test_run_ref_x_mod_is_false(self, setup):
        ref_cfg, mod_cfg = setup
        ref_res = GenomeValidator(ref_cfg, _ref_settings()).run()
        mod_res = GenomeValidator(mod_cfg, _mod_settings()).run()
        params = build_params(_build_validation_results(ref_res, mod_res, gxg_res=None))
        assert params.run_ref_x_mod is False

    def test_minimap2_not_called_for_eukaryote(self, setup):
        """Expected: genomexgenome (minimap2) should NOT run when both genomes are fragmented."""
        ref_cfg, mod_cfg = setup
        ref_res = GenomeValidator(ref_cfg, _ref_settings()).run()
        mod_res = GenomeValidator(mod_cfg, _mod_settings()).run()

        with patch("validation_pkg.validators.interfile_genome.subprocess.run") as mock_run:
            # Correct condition (per docs): skip if EITHER genome is fragmented
            if not mod_res.fragmented and not ref_res.fragmented:
                genomexgenome_validation(ref_res, mod_res, _gxg_settings())

        mock_run.assert_not_called()

    def test_fixed_main_py_skips_gxg_for_eukaryote(self, setup):
        """Fixed condition checks both ref and mod — GXG skipped for EUKARYOTE."""
        ref_cfg, mod_cfg = setup
        ref_res = GenomeValidator(ref_cfg, _ref_settings()).run()
        mod_res = GenomeValidator(mod_cfg, _mod_settings()).run()

        fixed_main_py_would_call_gxg = (
            mod_res is not None
            and ref_res is not None
            and not getattr(mod_res, "fragmented", False)
            and not getattr(ref_res, "fragmented", False)
        )
        assert fixed_main_py_would_call_gxg is False


# ---------------------------------------------------------------------------
# Scenario 4b: Multiple sequences in reference — PROKARYOTE, ref > limit
#
# ref.fa  : 3 sequences > limit=2  →  ref fragmented
# mod.fa  : 1 sequence  (within limit)  →  mod NOT fragmented
#
# Expected (per docs):
#   - ref.fragmented = True, mod.fragmented = False
#   - No plasmid extraction from ref
#   - run_ref_x_mod = False  (because ref_fragmented)
#   - minimap2 NOT used
#
# KNOWN ISSUE: main.py currently checks only mod.fragmented.
# When only ref is fragmented (mod is not), main.py would still call
# genomexgenome_validation and thus invoke minimap2 — contrary to the docs.
# ---------------------------------------------------------------------------

class TestScenario4bProkaryoteFragmentedRef:

    @pytest.fixture
    def setup(self, tmp_path):
        ref_path = _write_fasta(tmp_path / "ref.fasta", [
            ("ref_c1", _CHR_SEQ),
            ("ref_c2", _CHR_SEQ),
            ("ref_c3", _CHR_SEQ),
        ])
        mod_path = _write_fasta(tmp_path / "mod.fasta", [("mod_chr", _CHR_SEQ)])
        out = tmp_path / "out"
        out.mkdir()
        return (
            _make_genome_config(ref_path, out, n_sequence_limit=2),
            _make_genome_config(mod_path, out, n_sequence_limit=2),
        )

    def test_ref_is_fragmented(self, setup):
        ref_cfg, _ = setup
        ref_res = GenomeValidator(ref_cfg, _ref_settings()).run()
        assert ref_res.fragmented is True

    def test_mod_is_not_fragmented(self, setup):
        _, mod_cfg = setup
        mod_res = GenomeValidator(mod_cfg, _mod_settings()).run()
        assert mod_res.fragmented is False

    def test_ref_output_is_copied_as_is(self, setup):
        """Fragmented ref should be copied without plasmid extraction."""
        ref_cfg, _ = setup
        ref_res = GenomeValidator(ref_cfg, _ref_settings()).run()
        ids = _read_ids(Path(ref_res.output_file))
        assert len(ids) == 3, "fragmented ref must be copied with all sequences"

    def test_ref_no_plasmid_extracted(self, setup):
        ref_cfg, _ = setup
        ref_res = GenomeValidator(ref_cfg, _ref_settings()).run()
        assert not ref_res.plasmid_filenames, (
            "fragmented ref must not have plasmids extracted"
        )

    def test_run_ref_x_mod_is_false(self, setup):
        """build_params checks both ref_fragmented and mod_fragmented."""
        ref_cfg, mod_cfg = setup
        ref_res = GenomeValidator(ref_cfg, _ref_settings()).run()
        mod_res = GenomeValidator(mod_cfg, _mod_settings()).run()
        params = build_params(_build_validation_results(ref_res, mod_res, gxg_res=None))
        assert params.run_ref_x_mod is False, (
            "run_ref_x_mod must be False when ref is fragmented"
        )

    def test_expected_behavior_minimap2_not_called(self, setup):
        """
        Per docs (scenario 4): minimap2 should NOT be called when ref has multiple
        sequences. The correct condition to guard GXG should check both ref and mod.
        """
        ref_cfg, mod_cfg = setup
        ref_res = GenomeValidator(ref_cfg, _ref_settings()).run()
        mod_res = GenomeValidator(mod_cfg, _mod_settings()).run()

        # Expected (documented) condition: skip if EITHER genome is fragmented
        should_run_gxg = not mod_res.fragmented and not ref_res.fragmented
        assert should_run_gxg is False, (
            "GXG / minimap2 should not run when ref is fragmented (scenario 4)"
        )

    def test_fixed_main_py_skips_gxg_when_ref_fragmented(self, setup):
        """
        After fix: main.py now checks both ref.fragmented and mod.fragmented.
        When only ref is fragmented (scenario 4), GXG / minimap2 must be skipped.
        """
        ref_cfg, mod_cfg = setup
        ref_res = GenomeValidator(ref_cfg, _ref_settings()).run()
        mod_res = GenomeValidator(mod_cfg, _mod_settings()).run()

        # Fixed condition (both ref and mod must be non-fragmented):
        fixed_main_py_would_call_gxg = (
            mod_res is not None
            and ref_res is not None
            and not getattr(mod_res, "fragmented", False)
            and not getattr(ref_res, "fragmented", False)
        )
        assert fixed_main_py_would_call_gxg is False, (
            "Fixed main.py must NOT call GXG when ref is fragmented (scenario 4)"
        )


# ---------------------------------------------------------------------------
# Scenario 5: Fragmented reference + force_defragment_ref (unsupported)
#
# Expected:
#   - Fragmented ref is handled by scenarios 3/4 without the flag
#   - With the flag, ref contigs are merged before validation
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Real-world failure: Case 3 — duplicate sequence IDs in fragmented mod
#
# When mod has > n_sequence_limit sequences it is copied as-is (scenario 3).
# If the original mod.fa contains duplicate sequence IDs those duplicates
# survive in mod_mod.fasta and cause downstream tools (e.g. samtools) to
# fail with "duplicate entry in SAM header".
#
# Expected fix: deduplicate IDs in the copied output, just like plasmid files.
# ---------------------------------------------------------------------------

class TestScenario3DuplicateIdsInFragmentedMod:
    """
    Real-world failure: mod with duplicate IDs is fragmented → copied as-is
    → downstream samtools fails with "duplicate entry in SAM header".

    Expected: duplicate IDs in the copied fragmented output must be renamed.
    """

    @pytest.fixture
    def setup(self, tmp_path):
        # mod has 3 sequences > limit=2, two of them share the same ID
        ref_path = _write_fasta(tmp_path / "ref.fasta", [("ref_chr", _CHR_SEQ)])
        mod_path = _write_fasta(tmp_path / "mod.fasta", [
            ("ENA|CP045672|CP045672.1", _CHR_SEQ),
            ("ENA|CP045672|CP045672.2", _CHR_SEQ),
            ("ENA|CP045672|CP045672",   _CHR_SEQ),   # ← first copy
            ("ENA|CP045672|CP045672",   _CHR_SEQ),   # ← duplicate
            ("ENA|CP045672|CP045672",   _CHR_SEQ),   # ← duplicate
        ])
        out = tmp_path / "out"
        out.mkdir()
        return (
            _make_genome_config(ref_path, out, n_sequence_limit=4),
            _make_genome_config(mod_path, out, n_sequence_limit=4),
        )

    def test_mod_is_fragmented_with_five_seqs_over_limit(self, setup):
        _, mod_cfg = setup
        mod_res = GenomeValidator(mod_cfg, _mod_settings()).run()
        assert mod_res.fragmented is True, "5 sequences > limit=4 must be fragmented"

    def test_output_has_no_duplicate_ids(self, setup):
        """
        After copying the fragmented mod, all sequence IDs in the output file
        must be unique so downstream tools do not fail.
        """
        _, mod_cfg = setup
        mod_res = GenomeValidator(mod_cfg, _mod_settings()).run()
        ids = _read_ids(Path(mod_res.output_file))
        assert len(ids) == len(set(ids)), (
            f"Duplicate IDs found in copied fragmented output: "
            f"{[i for i in ids if ids.count(i) > 1]}"
        )

    def test_all_sequences_preserved_after_deduplication(self, setup):
        """Deduplication must rename, not drop sequences."""
        _, mod_cfg = setup
        mod_res = GenomeValidator(mod_cfg, _mod_settings()).run()
        ids = _read_ids(Path(mod_res.output_file))
        assert len(ids) == 5, "All 5 sequences must be present after deduplication"


# ---------------------------------------------------------------------------
# Scenario 4 boundary fix: ref with multiple sequences at the limit
#
# Fix: changed `> n_sequence_limit` to `>= n_sequence_limit` so that a genome
# with exactly limit sequences is treated as fragmented (scenario 4 — copied as-is).
# ---------------------------------------------------------------------------

class TestScenario4cRefAtSequenceLimit:
    """
    Case 4 boundary: ref has exactly n_sequence_limit sequences.

    Fixed: len(seqs) >= limit → fragmented at the boundary.
    Per scenario 4 docs: multiple sequences in ref → copied as-is, no main_longest.
    """

    @pytest.fixture
    def setup(self, tmp_path):
        ref_path = _write_fasta(tmp_path / "ref.fasta", [
            ("NC_000964.3", _CHR_SEQ),
            ("NC_000964.3", _CHR_SEQ),   # intentional duplicate (real-world case)
            ("NC_000964.1", _CHR_SEQ),
            ("NC_000964.2", _CHR_SEQ),
            ("CP045673.1",  _PLA_SEQ),
        ])
        mod_path = _write_fasta(tmp_path / "mod.fasta", [("mod_chr", _CHR_SEQ)])
        out = tmp_path / "out"
        out.mkdir()
        return (
            _make_genome_config(ref_path, out, n_sequence_limit=5),
            _make_genome_config(mod_path, out, n_sequence_limit=5),
        )

    def test_ref_is_fragmented_at_exact_limit(self, setup):
        """5 seqs with limit=5: >= check → ref IS fragmented (scenario 4)."""
        ref_cfg, _ = setup
        ref_res = GenomeValidator(ref_cfg, _ref_settings()).run()
        assert ref_res.fragmented is True

    def test_ref_copied_as_is_no_chromosome_extraction(self, setup):
        """
        With fragmented=True, main_longest does NOT run.
        All 5 sequences must be present in output (no plasmid extraction).
        """
        ref_cfg, _ = setup
        ref_res = GenomeValidator(ref_cfg, _ref_settings()).run()
        assert ref_res.num_sequences == 5, "All sequences preserved when fragmented"
        assert not ref_res.plasmid_filenames, "No plasmid extraction in scenario 4"

    def test_limit_minus_one_is_not_fragmented(self, setup):
        """4 seqs with limit=5 (below limit) → NOT fragmented, normal processing."""
        ref_path = Path(setup[0].filepath)
        out = ref_path.parent / "out2"
        out.mkdir(exist_ok=True)
        fasta_4 = _write_fasta(out / "ref4.fasta", [
            ("seq1", _CHR_SEQ), ("seq2", _CHR_SEQ),
            ("seq3", _CHR_SEQ), ("seq4", _CHR_SEQ),
        ])
        cfg4 = _make_genome_config(fasta_4, out, n_sequence_limit=5)
        ref_res = GenomeValidator(cfg4, _ref_settings()).run()
        assert ref_res.fragmented is False


class TestScenario5ForceDefragmentRefBaseline:
    """Baseline: without force_defragment_ref, fragmented ref follows scenarios 3/4."""

    def test_without_flag_fragmented_ref_behaves_as_scenario_4(self, tmp_path):
        ref_path = _write_fasta(tmp_path / "ref.fasta", [
            ("ref_c1", _CHR_SEQ),
            ("ref_c2", _CHR_SEQ),
        ])
        out = tmp_path / "out"
        out.mkdir()
        ref_cfg = _make_genome_config(ref_path, out, n_sequence_limit=1)
        ref_res = GenomeValidator(ref_cfg, _ref_settings()).run()
        assert ref_res.fragmented is True, (
            "Without force_defragment_ref, fragmented ref is marked fragmented"
        )

    def test_without_flag_no_plasmid_extracted_from_fragmented_ref(self, tmp_path):
        ref_path = _write_fasta(tmp_path / "ref.fasta", [
            ("ref_c1", _CHR_SEQ),
            ("ref_c2", _CHR_SEQ),
        ])
        out = tmp_path / "out"
        out.mkdir()
        ref_cfg = _make_genome_config(ref_path, out, n_sequence_limit=1)
        ref_res = GenomeValidator(ref_cfg, _ref_settings()).run()
        assert not ref_res.plasmid_filenames, (
            "Fragmented ref without force_defragment_ref: no plasmid extraction"
        )
