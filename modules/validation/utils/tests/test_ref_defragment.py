"""Unit tests for ref_defragment.py.

Tests cover:
- FASTA input: merging, sequence concatenation, TSV content, output paths
- GenBank input: merging from ORIGIN sections
- Edge cases: empty file, single contig
- TSV columns and offset tracking
- Output file locations (FASTA next to input, TSV in outputs/tables/)
"""

import sys
import gzip
import pytest
from pathlib import Path
from unittest.mock import MagicMock
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Ensure validation_utils and validation_pkg are importable when running
# pytest from this directory or from modules/validation/.
_VALIDATION_ROOT = Path(__file__).parent.parent.parent  # modules/validation/
if str(_VALIDATION_ROOT) not in sys.path:
    sys.path.insert(0, str(_VALIDATION_ROOT))

from validation_pkg.config_manager import GenomeConfig
from validation_pkg.utils.formats import GenomeFormat, CodingType
from ref_defragment import defragment_reference


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def logger():
    return MagicMock()


@pytest.fixture
def tmp_inputs(tmp_path):
    """Create the expected directory layout: <root>/inputs/ and <root>/outputs/tables/."""
    inputs = tmp_path / "inputs"
    inputs.mkdir()
    (tmp_path / "outputs" / "tables").mkdir(parents=True)
    return inputs


def _make_genome_config(filepath: Path, fmt: GenomeFormat) -> GenomeConfig:
    return GenomeConfig(
        filename=filepath.name,
        filepath=filepath,
        coding_type=CodingType.NONE,
        detected_format=fmt,
        output_dir=filepath.parent,
        global_options={},
    )


def _write_fasta(path: Path, records: list[tuple[str, str]]) -> None:
    with open(path, "w") as fh:
        for seq_id, seq in records:
            fh.write(f">{seq_id}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i : i + 60] + "\n")


def _write_genbank(path: Path, records: list[tuple[str, str]]) -> None:
    """Write a minimal multi-record GenBank file."""
    biopython_records = [
        SeqRecord(
            Seq(seq),
            id=seq_id,
            name=seq_id,
            description="test",
            annotations={"molecule_type": "DNA"},
        )
        for seq_id, seq in records
    ]
    with open(path, "w") as fh:
        SeqIO.write(biopython_records, fh, "genbank")


# ---------------------------------------------------------------------------
# FASTA tests
# ---------------------------------------------------------------------------

class TestDefragmentFasta:

    CONTIGS = [
        ("contig_1", "ATCGATCGATCG"),
        ("contig_2", "GGGGCCCCTTTT"),
        ("contig_3", "AAAATTTTCCCC"),
    ]
    EXPECTED_SEQ = "ATCGATCGATCGGGGGCCCCTTTTAAAATTTTCCCC"

    @pytest.fixture
    def fasta_config(self, tmp_inputs):
        fasta_path = tmp_inputs / "ref.fasta"
        _write_fasta(fasta_path, self.CONTIGS)
        return _make_genome_config(fasta_path, GenomeFormat.FASTA)

    def test_returns_two_paths(self, fasta_config, logger):
        merged, tsv = defragment_reference(fasta_config, logger)
        assert isinstance(merged, Path)
        assert isinstance(tsv, Path)

    def test_merged_fasta_exists(self, fasta_config, logger):
        merged, _ = defragment_reference(fasta_config, logger)
        assert merged.exists()

    def test_tsv_exists(self, fasta_config, logger):
        _, tsv = defragment_reference(fasta_config, logger)
        assert tsv.exists()

    def test_merged_fasta_has_one_sequence(self, fasta_config, logger):
        merged, _ = defragment_reference(fasta_config, logger)
        records = list(SeqIO.parse(merged, "fasta"))
        assert len(records) == 1

    def test_merged_sequence_is_concatenation(self, fasta_config, logger):
        merged, _ = defragment_reference(fasta_config, logger)
        records = list(SeqIO.parse(merged, "fasta"))
        assert str(records[0].seq) == self.EXPECTED_SEQ

    def test_merged_fasta_id_contains_basename(self, fasta_config, logger):
        merged, _ = defragment_reference(fasta_config, logger)
        records = list(SeqIO.parse(merged, "fasta"))
        assert "ref" in records[0].id

    def test_fasta_written_next_to_input(self, fasta_config, logger):
        merged, _ = defragment_reference(fasta_config, logger)
        assert merged.parent == fasta_config.filepath.parent

    def test_tsv_written_to_outputs_tables(self, fasta_config, logger):
        _, tsv = defragment_reference(fasta_config, logger)
        assert "outputs" in tsv.parts
        assert "tables" in tsv.parts

    def test_tsv_header(self, fasta_config, logger):
        _, tsv = defragment_reference(fasta_config, logger)
        lines = tsv.read_text().splitlines()
        assert lines[0] == "seq_id\tlength\tstart"

    def test_tsv_row_count(self, fasta_config, logger):
        _, tsv = defragment_reference(fasta_config, logger)
        lines = tsv.read_text().splitlines()
        assert len(lines) == 1 + len(self.CONTIGS)  # header + one row per contig

    def test_tsv_seq_ids(self, fasta_config, logger):
        _, tsv = defragment_reference(fasta_config, logger)
        rows = tsv.read_text().splitlines()[1:]
        ids = [r.split("\t")[0] for r in rows]
        assert ids == [c[0] for c in self.CONTIGS]

    def test_tsv_lengths(self, fasta_config, logger):
        _, tsv = defragment_reference(fasta_config, logger)
        rows = tsv.read_text().splitlines()[1:]
        lengths = [int(r.split("\t")[1]) for r in rows]
        assert lengths == [len(seq) for _, seq in self.CONTIGS]

    def test_tsv_start_offsets(self, fasta_config, logger):
        _, tsv = defragment_reference(fasta_config, logger)
        rows = tsv.read_text().splitlines()[1:]
        starts = [int(r.split("\t")[2]) for r in rows]
        expected = [0, 12, 24]
        assert starts == expected

    def test_logger_warns_for_each_contig(self, fasta_config, logger):
        defragment_reference(fasta_config, logger)
        # At least one warning call per contig
        warning_calls = logger.warning.call_count
        assert warning_calls >= len(self.CONTIGS)


# ---------------------------------------------------------------------------
# GenBank tests
# ---------------------------------------------------------------------------

class TestDefragmentGenbank:

    CONTIGS = [
        ("rec1", "ATCGATCGATCGATCGATCG"),
        ("rec2", "TTTTAAAACCCCGGGG"),
    ]
    EXPECTED_SEQ = "ATCGATCGATCGATCGATCGTTTTAAAACCCCGGGG"

    @pytest.fixture
    def gbk_config(self, tmp_inputs):
        gbk_path = tmp_inputs / "ref.gbk"
        _write_genbank(gbk_path, self.CONTIGS)
        return _make_genome_config(gbk_path, GenomeFormat.GENBANK)

    def test_merged_fasta_has_one_sequence(self, gbk_config, logger):
        merged, _ = defragment_reference(gbk_config, logger)
        records = list(SeqIO.parse(merged, "fasta"))
        assert len(records) == 1

    def test_merged_sequence_is_concatenation(self, gbk_config, logger):
        merged, _ = defragment_reference(gbk_config, logger)
        records = list(SeqIO.parse(merged, "fasta"))
        assert str(records[0].seq) == self.EXPECTED_SEQ

    def test_tsv_row_count(self, gbk_config, logger):
        _, tsv = defragment_reference(gbk_config, logger)
        lines = tsv.read_text().splitlines()
        assert len(lines) == 1 + len(self.CONTIGS)

    def test_tsv_start_offsets(self, gbk_config, logger):
        _, tsv = defragment_reference(gbk_config, logger)
        rows = tsv.read_text().splitlines()[1:]
        starts = [int(r.split("\t")[2]) for r in rows]
        assert starts == [0, 20]


# ---------------------------------------------------------------------------
# Compressed input tests
# ---------------------------------------------------------------------------

class TestDefragmentCompressed:

    CONTIGS = [("c1", "ATCG" * 10), ("c2", "GCTA" * 10)]

    @pytest.fixture
    def gz_config(self, tmp_inputs):
        gz_path = tmp_inputs / "ref.fasta.gz"
        fasta_lines = []
        for seq_id, seq in self.CONTIGS:
            fasta_lines.append(f">{seq_id}\n{seq}\n")
        with gzip.open(gz_path, "wt") as fh:
            fh.write("".join(fasta_lines))
        config = _make_genome_config(gz_path, GenomeFormat.FASTA)
        config.coding_type = CodingType.GZIP
        return config

    def test_merged_sequence_correct(self, gz_config, logger):
        merged, _ = defragment_reference(gz_config, logger)
        records = list(SeqIO.parse(merged, "fasta"))
        expected = "".join(seq for _, seq in self.CONTIGS)
        assert str(records[0].seq) == expected


# ---------------------------------------------------------------------------
# Edge case tests
# ---------------------------------------------------------------------------

class TestDefragmentEdgeCases:

    def test_empty_fasta_raises(self, tmp_inputs, logger):
        empty = tmp_inputs / "empty.fasta"
        empty.write_text("")
        config = _make_genome_config(empty, GenomeFormat.FASTA)
        with pytest.raises(ValueError, match="No sequences found"):
            defragment_reference(config, logger)

    def test_single_contig_works(self, tmp_inputs, logger):
        fasta_path = tmp_inputs / "single.fasta"
        _write_fasta(fasta_path, [("only", "ATCGATCG")])
        config = _make_genome_config(fasta_path, GenomeFormat.FASTA)
        merged, tsv = defragment_reference(config, logger)
        records = list(SeqIO.parse(merged, "fasta"))
        assert len(records) == 1
        assert str(records[0].seq) == "ATCGATCG"
        rows = tsv.read_text().splitlines()
        assert len(rows) == 2  # header + 1 row

    def test_empty_fasta_cleans_up_partial_files(self, tmp_inputs, logger):
        empty = tmp_inputs / "empty.fasta"
        empty.write_text("")
        config = _make_genome_config(empty, GenomeFormat.FASTA)
        with pytest.raises(ValueError):
            defragment_reference(config, logger)
        merged_path = tmp_inputs / "empty_defragmented.fasta"
        assert not merged_path.exists()
