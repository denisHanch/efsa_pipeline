"""
Tests for inter-file genome validation.

Tests cover:
- Genome-to-genome consistency checks
- Sequence ID matching between reference and modified genomes
- GenomeXGenomeSettings configuration
- genomexgenome_validation function
- Error handling for mismatched genomes
"""

import pytest
from validation_pkg.validators.interfile_genome import GenomeXGenomeSettings, genomexgenome_validation
from validation_pkg.validators.genome_validator import GenomeOutputMetadata
from validation_pkg.exceptions import GenomeValidationError


class TestGenomeXGenomeSettings:
    """Test GenomeXGenomeSettings configuration."""

    def test_default_settings(self):
        """Test default settings values."""
        settings = GenomeXGenomeSettings()
        assert settings.same_number_of_sequences is True
        assert settings.same_sequence_ids is False
        assert settings.same_sequence_lengths is False

    def test_update_settings(self):
        """Test immutable update pattern."""
        settings = GenomeXGenomeSettings()
        updated = settings.update(same_sequence_ids=True)

        assert updated.same_sequence_ids is True
        assert settings.same_sequence_ids is False  # Original unchanged

    def test_length_check_requires_id_check(self):
        """Test that same_sequence_lengths requires same_sequence_ids=True."""
        with pytest.raises(ValueError) as exc_info:
            GenomeXGenomeSettings(
                same_sequence_lengths=True,
                same_sequence_ids=False
            )

        assert 'requires same_sequence_ids=True' in str(exc_info.value)


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

        result = genomexgenome_validation(ref_result, mod_result)

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

        with pytest.raises(GenomeValidationError) as exc_info:
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

        settings = GenomeXGenomeSettings(same_sequence_ids=True)
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

        with pytest.raises(GenomeValidationError) as exc_info:
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

        with pytest.raises(GenomeValidationError) as exc_info:
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

        with pytest.raises(GenomeValidationError) as exc_info:
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
            same_sequence_lengths=True
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

        with pytest.raises(GenomeValidationError) as exc_info:
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

        with pytest.raises(GenomeValidationError) as exc_info:
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
            same_sequence_lengths=True
        )

        result = genomexgenome_validation(ref_result, mod_result, settings)

        assert len(result['warnings']) > 0
        assert 'Missing length information' in result['warnings'][0]


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
