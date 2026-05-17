"""Tests for inter-file genome validation."""

import tempfile
from pathlib import Path
from types import SimpleNamespace

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from validation_pkg.validators.interfile_genome import (
    GenomeXGenomeSettings,
    genomexgenome_validation,
    _deduplicate_plasmid_ids,
)
from validation_pkg.validators.genome_validator import GenomeOutputMetadata
from validation_pkg.exceptions import InterFileValidationError


class TestGenomeXGenomeSettings:
    """Test GenomeXGenomeSettings configuration."""

    def test_default_settings(self):
        """Test default settings values."""
        settings = GenomeXGenomeSettings()
        assert settings.same_number_of_sequences is True
        assert settings.same_sequence_ids is False
        assert settings.same_sequence_lengths is False
        assert settings.characterize is True

    def test_update_settings(self):
        """Test immutable update pattern."""
        settings = GenomeXGenomeSettings()
        updated = settings.update(same_sequence_ids=True)

        assert updated.same_sequence_ids is True
        assert settings.same_sequence_ids is False  # Original unchanged

    # def test_length_check_requires_id_check(self):
    #     """Test that same_sequence_lengths requires same_sequence_ids=True."""
    #     with pytest.raises(ValueError) as exc_info:
    #         GenomeXGenomeSettings(
    #             same_sequence_lengths=True,
    #             same_sequence_ids=False
    #         )

    #     assert 'requires same_sequence_ids=True' in str(exc_info.value)


class TestSequenceCountValidation:
    """Test sequence count matching validation."""

    def test_same_count_passes(self):
        """Test that same sequence counts pass."""
        ref_result = GenomeOutputMetadata(
            output_file='ref_genome.fasta',
            num_sequences=3,
            sequence_ids=['chr1', 'chr2', 'chr3'],
            sequence_lengths={'chr1': 5000000, 'chr2': 3000000, 'chr3': 2000000}
        )
        mod_result = GenomeOutputMetadata(
            output_file='mod_genome.fasta',
            num_sequences=3,
            sequence_ids=['chr1', 'chr2', 'chr3'],
            sequence_lengths={'chr1': 5000000, 'chr2': 3000000, 'chr3': 2000000}
        )

        result = genomexgenome_validation(ref_result, mod_result, GenomeXGenomeSettings(characterize=False))

        assert result['passed'] is True
        assert len(result['errors']) == 0

    def test_different_count_fails(self):
        """Test that different sequence counts raise error."""
        ref_result = GenomeOutputMetadata(
            output_file='ref_genome.fasta',
            num_sequences=3,
            sequence_ids=['chr1', 'chr2', 'chr3'],
            sequence_lengths={}
        )
        mod_result = GenomeOutputMetadata(
            output_file='mod_genome.fasta',
            num_sequences=2,
            sequence_ids=['chr1', 'chr2'],
            sequence_lengths={}
        )

        with pytest.raises(InterFileValidationError) as exc_info:
            genomexgenome_validation(ref_result, mod_result)

        error_msg = str(exc_info.value)
        assert 'sequence count mismatch' in error_msg.lower()
        assert '3' in error_msg
        assert '2' in error_msg


class TestSequenceIDValidation:
    """Test sequence ID matching validation."""

    def test_same_ids_passes(self):
        """Test that matching sequence IDs pass."""
        ref_result = GenomeOutputMetadata(
            output_file='ref_genome.fasta',
            num_sequences=2,
            sequence_ids=['chr1', 'chr2'],
            sequence_lengths={'chr1': 5000000, 'chr2': 3000000}
        )
        mod_result = GenomeOutputMetadata(
            output_file='mod_genome.fasta',
            num_sequences=2,
            sequence_ids=['chr2', 'chr1'],
            sequence_lengths={'chr1': 5000000, 'chr2': 3000000}
        )

        settings = GenomeXGenomeSettings(same_sequence_ids=True, characterize=False)
        result = genomexgenome_validation(ref_result, mod_result, settings)

        assert result['passed'] is True
        assert len(result['errors']) == 0
        assert set(result['metadata']['common_sequence_ids']) == {'chr1', 'chr2'}

    def test_different_ids_fails(self):
        """Test that mismatched sequence IDs raise error."""
        ref_result = GenomeOutputMetadata(
            output_file='ref_genome.fasta',
            num_sequences=2,
            sequence_ids=['chr1', 'chr2'],
            sequence_lengths={}
        )
        mod_result = GenomeOutputMetadata(
            output_file='mod_genome.fasta',
            num_sequences=2,
            sequence_ids=['chr1', 'chr3'],
            sequence_lengths={}
        )

        settings = GenomeXGenomeSettings(same_sequence_ids=True)

        with pytest.raises(InterFileValidationError) as exc_info:
            genomexgenome_validation(ref_result, mod_result, settings)

        error_msg = str(exc_info.value)
        assert 'sequence id mismatch' in error_msg.lower()  # All lowercase since we call .lower()
        assert 'chr2' in error_msg  # Reference-only
        assert 'chr3' in error_msg  # Modified-only

    def test_ref_only_ids(self):
        """Test detection of reference-only sequence IDs."""
        ref_result = GenomeOutputMetadata(
            output_file='ref_genome.fasta',
            num_sequences=3,
            sequence_ids=['chr1', 'chr2', 'plasmid1'],
            sequence_lengths={}
        )
        mod_result = GenomeOutputMetadata(
            output_file='mod_genome.fasta',
            num_sequences=2,
            sequence_ids=['chr1', 'chr2'],
            sequence_lengths={}
        )

        settings = GenomeXGenomeSettings(same_sequence_ids=True)

        with pytest.raises(InterFileValidationError) as exc_info:
            genomexgenome_validation(ref_result, mod_result, settings)

        error_msg = str(exc_info.value)
        assert 'plasmid1' in error_msg
        assert 'reference-only' in error_msg.lower()

    def test_mod_only_ids(self):
        """Test detection of modified-only sequence IDs."""
        ref_result = GenomeOutputMetadata(
            output_file='ref_genome.fasta',
            num_sequences=2,
            sequence_ids=['chr1', 'chr2'],
            sequence_lengths={}
        )
        mod_result = GenomeOutputMetadata(
            output_file='mod_genome.fasta',
            num_sequences=3,
            sequence_ids=['chr1', 'chr2', 'insert1'],
            sequence_lengths={}
        )

        settings = GenomeXGenomeSettings(same_sequence_ids=True)

        with pytest.raises(InterFileValidationError) as exc_info:
            genomexgenome_validation(ref_result, mod_result, settings)

        error_msg = str(exc_info.value)
        assert 'insert1' in error_msg
        assert 'modified-only' in error_msg.lower()


class TestSequenceLengthValidation:
    """Test sequence length matching validation."""

    def test_same_lengths_passes(self):
        """Test that matching sequence lengths pass."""
        ref_result = GenomeOutputMetadata(
            output_file='ref_genome.fasta',
            num_sequences=2,
            sequence_ids=['chr1', 'chr2'],
            sequence_lengths={'chr1': 5000000, 'chr2': 3000000}
        )
        mod_result = GenomeOutputMetadata(
            output_file='mod_genome.fasta',
            num_sequences=2,
            sequence_ids=['chr1', 'chr2'],
            sequence_lengths={'chr1': 5000000, 'chr2': 3000000}
        )

        settings = GenomeXGenomeSettings(
            same_sequence_ids=True,
            same_sequence_lengths=True,
            characterize=False
        )
        result = genomexgenome_validation(ref_result, mod_result, settings)

        assert result['passed'] is True
        assert len(result['errors']) == 0
        assert len(result['metadata']['length_mismatches']) == 0

    def test_different_lengths_fails(self):
        """Test that different sequence lengths raise error."""
        ref_result = GenomeOutputMetadata(
            output_file='ref_genome.fasta',
            num_sequences=2,
            sequence_ids=['chr1', 'chr2'],
            sequence_lengths={'chr1': 5000000, 'chr2': 3000000}
        )
        mod_result = GenomeOutputMetadata(
            output_file='mod_genome.fasta',
            num_sequences=2,
            sequence_ids=['chr1', 'chr2'],
            sequence_lengths={'chr1': 5001000, 'chr2': 3000000}
        )

        settings = GenomeXGenomeSettings(
            same_sequence_ids=True,
            same_sequence_lengths=True
        )

        with pytest.raises(InterFileValidationError) as exc_info:
            genomexgenome_validation(ref_result, mod_result, settings)

        error_msg = str(exc_info.value)
        assert 'length mismatch' in error_msg.lower()
        assert 'chr1' in error_msg
        assert '5000000' in error_msg
        assert '5001000' in error_msg

    def test_length_difference_reported(self):
        """Test that length differences are calculated correctly."""
        ref_result = GenomeOutputMetadata(
            output_file='ref_genome.fasta',
            num_sequences=1,
            sequence_ids=['chr1'],
            sequence_lengths={'chr1': 5000000}
        )
        mod_result = GenomeOutputMetadata(
            output_file='mod_genome.fasta',
            num_sequences=1,
            sequence_ids=['chr1'],
            sequence_lengths={'chr1': 5001500}
        )

        settings = GenomeXGenomeSettings(
            same_sequence_ids=True,
            same_sequence_lengths=True
        )

        with pytest.raises(InterFileValidationError) as exc_info:
            genomexgenome_validation(ref_result, mod_result, settings)

        error_msg = str(exc_info.value)
        assert '+1500' in error_msg  # Positive difference

    def test_missing_length_info_warns(self):
        """Test that missing length information produces warning."""
        ref_result = GenomeOutputMetadata(
            output_file='ref_genome.fasta',
            num_sequences=1,
            sequence_ids=['chr1'],
            sequence_lengths={}
        )
        mod_result = GenomeOutputMetadata(
            output_file='mod_genome.fasta',
            num_sequences=1,
            sequence_ids=['chr1'],
            sequence_lengths={'chr1': 5000000}
        )

        settings = GenomeXGenomeSettings(
            same_sequence_ids=True,
            same_sequence_lengths=True,
            characterize=False
        )

        result = genomexgenome_validation(ref_result, mod_result, settings)

        assert len(result['warnings']) > 0
        assert 'Missing length information' in result['warnings'][0]


# ---------------------------------------------------------------------------
# Helpers shared by duplicate-ID tests
# ---------------------------------------------------------------------------

def _make_plasmid_fasta(tmp_path: Path, records: list) -> Path:
    """Write a FASTA file from (id, seq) tuples and return its Path."""
    p = tmp_path / "plasmid.fasta"
    seqs = [SeqRecord(Seq(seq), id=sid, description='') for sid, seq in records]
    with open(p, 'w') as fh:
        SeqIO.write(seqs, fh, 'fasta')
    return p


def _read_ids(path: Path) -> list:
    return [r.id for r in SeqIO.parse(str(path), 'fasta')]


def _base_genomes():
    ref = GenomeOutputMetadata(
        output_file='ref.fasta',
        num_sequences=1,
        sequence_ids=['chr'],
        sequence_lengths={'chr': 1000},
    )
    mod = GenomeOutputMetadata(
        output_file='mod.fasta',
        num_sequences=1,
        sequence_ids=['chr'],
        sequence_lengths={'chr': 1000},
    )
    return ref, mod


# ---------------------------------------------------------------------------
# _deduplicate_plasmid_ids — unit tests
# ---------------------------------------------------------------------------

class TestDeduplicatePlasmidIds:

    def test_no_duplicates_returns_empty(self, tmp_path):
        path = _make_plasmid_fasta(tmp_path, [('pA', 'ATCG'), ('pB', 'GCTA')])
        renamed = _deduplicate_plasmid_ids(str(path), None)
        assert renamed == []
        assert _read_ids(path) == ['pA', 'pB']

    def test_single_duplicate_renamed(self, tmp_path):
        path = _make_plasmid_fasta(tmp_path, [('pA', 'ATCG'), ('pA', 'GCTA'), ('pB', 'TTTT')])
        renamed = _deduplicate_plasmid_ids(str(path), None)
        assert len(renamed) == 1
        assert renamed[0] == ('pA', 'pA_1')
        assert _read_ids(path) == ['pA', 'pA_1', 'pB']

    def test_multiple_duplicates_renamed_incrementally(self, tmp_path):
        path = _make_plasmid_fasta(tmp_path, [
            ('pA', 'ATCG'), ('pA', 'GCTA'), ('pA', 'TTTT'),
        ])
        renamed = _deduplicate_plasmid_ids(str(path), None)
        assert [new for _, new in renamed] == ['pA_1', 'pA_2']
        assert _read_ids(path) == ['pA', 'pA_1', 'pA_2']

    def test_two_distinct_duplicate_ids(self, tmp_path):
        path = _make_plasmid_fasta(tmp_path, [
            ('pA', 'ATCG'), ('pB', 'GCTA'), ('pA', 'TTTT'), ('pB', 'CCCC'),
        ])
        renamed = _deduplicate_plasmid_ids(str(path), None)
        assert len(renamed) == 2
        ids_after = _read_ids(path)
        assert 'pA' in ids_after and 'pA_1' in ids_after
        assert 'pB' in ids_after and 'pB_1' in ids_after

    def test_file_not_rewritten_when_no_duplicates(self, tmp_path):
        path = _make_plasmid_fasta(tmp_path, [('pA', 'ATCG')])
        mtime_before = path.stat().st_mtime_ns
        _deduplicate_plasmid_ids(str(path), None)
        assert path.stat().st_mtime_ns == mtime_before


# ---------------------------------------------------------------------------
# genomexgenome_validation — duplicate-ID integration tests
# ---------------------------------------------------------------------------

class TestPlasmidDuplicateIdCheck:

    def test_no_warning_when_no_duplicates(self, tmp_path):
        path = _make_plasmid_fasta(tmp_path, [('pA', 'ATCG'), ('pB', 'GCTA')])
        ref, mod = _base_genomes()
        ref_plasmid = SimpleNamespace(output_file=str(path))
        result = genomexgenome_validation(
            ref, mod,
            GenomeXGenomeSettings(characterize=False, same_number_of_sequences=False),
            ref_plasmid_result=ref_plasmid,
        )
        assert not any('duplicate' in w.lower() for w in result['warnings'])

    def test_warning_emitted_for_ref_plasmid_duplicates(self, tmp_path):
        path = _make_plasmid_fasta(tmp_path, [('pA', 'ATCG'), ('pA', 'GCTA')])
        ref, mod = _base_genomes()
        ref_plasmid = SimpleNamespace(output_file=str(path))
        result = genomexgenome_validation(
            ref, mod,
            GenomeXGenomeSettings(characterize=False, same_number_of_sequences=False),
            ref_plasmid_result=ref_plasmid,
        )
        assert any('ref' in w and 'duplicate' in w.lower() for w in result['warnings'])

    def test_warning_emitted_for_mod_plasmid_duplicates(self, tmp_path):
        path = _make_plasmid_fasta(tmp_path, [('pA', 'ATCG'), ('pA', 'GCTA')])
        ref, mod = _base_genomes()
        mod_plasmid = SimpleNamespace(output_file=str(path))
        result = genomexgenome_validation(
            ref, mod,
            GenomeXGenomeSettings(characterize=False, same_number_of_sequences=False),
            mod_plasmid_result=mod_plasmid,
        )
        assert any('mod' in w and 'duplicate' in w.lower() for w in result['warnings'])

    def test_duplicates_renamed_in_place(self, tmp_path):
        path = _make_plasmid_fasta(tmp_path, [('pA', 'ATCG'), ('pA', 'GCTA')])
        ref, mod = _base_genomes()
        ref_plasmid = SimpleNamespace(output_file=str(path))
        genomexgenome_validation(
            ref, mod,
            GenomeXGenomeSettings(characterize=False, same_number_of_sequences=False),
            ref_plasmid_result=ref_plasmid,
        )
        assert _read_ids(path) == ['pA', 'pA_1']

    def test_validation_still_passes_despite_duplicates(self, tmp_path):
        path = _make_plasmid_fasta(tmp_path, [('pA', 'ATCG'), ('pA', 'GCTA')])
        ref, mod = _base_genomes()
        ref_plasmid = SimpleNamespace(output_file=str(path))
        result = genomexgenome_validation(
            ref, mod,
            GenomeXGenomeSettings(characterize=False, same_number_of_sequences=False),
            ref_plasmid_result=ref_plasmid,
        )
        assert result['passed'] is True


class TestMetadataErrors:
    """Test error handling for missing metadata."""

    def test_missing_ref_metadata_fails(self):
        """Test that missing reference metadata raises error."""
        ref_result = GenomeOutputMetadata(
            output_file='mod_genome.fasta',
            num_sequences=1,
            sequence_ids=[],
            sequence_lengths={}
        )
