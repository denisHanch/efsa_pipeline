"""
Tests for inter-file read validation.

Tests cover:
- Paired-end read completeness validation
- Read ID matching between R1 and R2 files
- ReadXReadSettings configuration
- readxread_validation function
- Error handling for incomplete pairs
"""

import pytest
from validation_pkg.validators.interfile_read import ReadXReadSettings, readxread_validation
from validation_pkg.validators.read_validator import ReadOutputMetadata
from validation_pkg.exceptions import ReadValidationError


class TestReadXReadSettings:
    """Test ReadXReadSettings configuration."""

    def test_default_settings(self):
        """Test default settings values."""
        settings = ReadXReadSettings()
        assert settings.pair_end_basename is True
        assert settings.allow_missing_r1 is False

    def test_update_settings(self):
        """Test immutable update pattern."""
        settings = ReadXReadSettings()
        updated = settings.update(allow_missing_r1=True)

        assert updated.allow_missing_r1 is True
        assert settings.allow_missing_r1 is False  # Original unchanged


class TestPairedEndValidation:
    """Test paired-end completeness validation."""

    def test_complete_pairs_pass(self):
        """Test that complete R1+R2 pairs pass validation."""
        reads_results = [
            ReadOutputMetadata(
                output_file='sample_R1.fastq.gz',
                base_name='sample',
                read_number=1,
                illumina_pairing_detected='illumina',
                num_reads=1000
            ),
            ReadOutputMetadata(
                output_file='sample_R2.fastq.gz',
                base_name='sample',
                read_number=2,
                illumina_pairing_detected='illumina',
                num_reads=1000
            )
        ]

        result = readxread_validation(reads_results)

        assert result['passed'] is True
        assert len(result['errors']) == 0
        assert result['metadata']['pairs_checked'] == 1
        assert 'sample' in result['metadata']['complete_pairs']

    def test_missing_r1_fails(self):
        """Test that R2 without R1 raises error."""
        reads_results = [
            ReadOutputMetadata(
                output_file='sample_R2.fastq.gz',
                base_name='sample',
                read_number=2,
                illumina_pairing_detected='illumina',
                num_reads=1000
            )
        ]

        with pytest.raises(ReadValidationError) as exc_info:
            readxread_validation(reads_results)

        assert 'without matching r1' in str(exc_info.value).lower()  # All lowercase since we call .lower()
        assert 'sample' in str(exc_info.value)

    def test_r1_only_passes(self):
        """Test that R1-only files pass (no R2 required)."""
        reads_results = [
            ReadOutputMetadata(
                output_file='sample_R1.fastq.gz',
                base_name='sample',
                read_number=1,
                illumina_pairing_detected='illumina',
                num_reads=1000
            )
        ]

        result = readxread_validation(reads_results)

        assert result['passed'] is True
        assert len(result['errors']) == 0

    def test_multiple_missing_r1(self):
        """Test multiple R2 files without R1."""
        reads_results = [
            ReadOutputMetadata(
                output_file='sample1_R2.fastq.gz',
                base_name='sample1',
                read_number=2,
                illumina_pairing_detected='illumina',
                num_reads=1000
            ),
            ReadOutputMetadata(
                output_file='sample2_R2.fastq.gz',
                base_name='sample2',
                read_number=2,
                illumina_pairing_detected='illumina',
                num_reads=1000
            )
        ]

        with pytest.raises(ReadValidationError) as exc_info:
            readxread_validation(reads_results)

        error_msg = str(exc_info.value)
        assert 'sample1' in error_msg
        assert 'sample2' in error_msg

    def test_mixed_complete_incomplete_fails(self):
        """Test mix of complete and incomplete pairs."""
        reads_results = [
            # Complete pair
            ReadOutputMetadata(
                output_file='sample1_R1.fastq.gz',
                base_name='sample1',
                read_number=1,
                illumina_pairing_detected='illumina',
                num_reads=1000
            ),
            ReadOutputMetadata(
                output_file='sample1_R2.fastq.gz',
                base_name='sample1',
                read_number=2,
                illumina_pairing_detected='illumina',
                num_reads=1000
            ),
            # Incomplete pair (R2 only)
            ReadOutputMetadata(
                output_file='sample2_R2.fastq.gz',
                base_name='sample2',
                read_number=2,
                illumina_pairing_detected='illumina',
                num_reads=1000
            )
        ]

        with pytest.raises(ReadValidationError) as exc_info:
            readxread_validation(reads_results)

        assert 'sample2' in str(exc_info.value)

    def test_no_patterns_detected_passes(self):
        """Test that files without paired-end patterns pass."""
        reads_results = [
            ReadOutputMetadata(
                output_file='sample.fastq.gz',
                illumina_pairing_detected='ont',
                num_reads=1000
            )
        ]

        result = readxread_validation(reads_results)

        assert result['passed'] is True
        assert result['metadata']['pairs_checked'] == 0

    def test_allow_missing_r1_passes(self):
        """Test that allow_missing_r1=True allows orphan R2."""
        reads_results = [
            ReadOutputMetadata(
                output_file='sample_R2.fastq.gz',
                base_name='sample',
                read_number=2,
                illumina_pairing_detected='illumina',
                num_reads=1000
            )
        ]

        settings = ReadXReadSettings(allow_missing_r1=True)
        result = readxread_validation(reads_results, settings)

        assert result['passed'] is True
        assert len(result['warnings']) > 0  # Should warn but not error
        assert 'sample' in result['warnings'][0]

    def test_duplicate_r1_warns(self):
        """Test that multiple R1 files with same base name warn."""
        reads_results = [
            ReadOutputMetadata(
                output_file='sample_R1_copy1.fastq.gz',
                base_name='sample',
                read_number=1,
                illumina_pairing_detected='illumina',
                num_reads=1000
            ),
            ReadOutputMetadata(
                output_file='sample_R1_copy2.fastq.gz',
                base_name='sample',
                read_number=1,
                illumina_pairing_detected='illumina',
                num_reads=1000
            )
        ]

        result = readxread_validation(reads_results)

        assert result['passed'] is True
        assert len(result['warnings']) > 0
        assert 'Multiple R1' in result['warnings'][0]

    def test_duplicate_r2_warns(self):
        """Test that multiple R2 files with same base name warn."""
        reads_results = [
            ReadOutputMetadata(
                output_file='sample_R1.fastq.gz',
                base_name='sample',
                read_number=1,
                illumina_pairing_detected='illumina',
                num_reads=1000
            ),
            ReadOutputMetadata(
                output_file='sample_R2_copy1.fastq.gz',
                base_name='sample',
                read_number=2,
                illumina_pairing_detected='illumina',
                num_reads=1000
            ),
            ReadOutputMetadata(
                output_file='sample_R2_copy2.fastq.gz',
                base_name='sample',
                read_number=2,
                illumina_pairing_detected='illumina',
                num_reads=1000
            )
        ]

        result = readxread_validation(reads_results)

        assert result['passed'] is True
        assert len(result['warnings']) > 0
        assert 'Multiple R2' in result['warnings'][0]

    def test_empty_results_passes(self):
        """Test that empty results list passes."""
        result = readxread_validation([])

        assert result['passed'] is True
        assert result['metadata']['pairs_checked'] == 0

    def test_check_disabled_passes(self):
        """Test that validation passes when check is disabled."""
        reads_results = [
            ReadOutputMetadata(
                output_file='sample_R2.fastq.gz',
                base_name='sample',
                read_number=2,
                illumina_pairing_detected='illumina',
                num_reads=1000
            )
        ]

        settings = ReadXReadSettings(pair_end_basename=False)
        result = readxread_validation(reads_results, settings)

        assert result['passed'] is True
        assert len(result['errors']) == 0
