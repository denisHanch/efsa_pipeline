"""
Unit tests for the custom exceptions module.

Tests the exception hierarchy and ensures all custom exceptions:
1. Can be raised and caught correctly
2. Inherit from the proper base classes
3. Store and display error messages correctly
4. Follow the documented hierarchy
"""

import pytest
from validation_pkg.exceptions import (
    ValidationError,
    ConfigurationError,
    FileNotFoundError,
    FileFormatError,
    FastaFormatError,
    GenBankFormatError,
    BedFormatError,
    GffFormatError,
    FastqFormatError,
    BamFormatError,
    CompressionError,
    GenomeValidationError,
    FeatureValidationError,
    ReadValidationError,
    InterFileValidationError,
)


class TestValidationErrorBase:
    """Tests for the base ValidationError exception."""

    def test_can_raise_validation_error(self):
        """Test that ValidationError can be raised."""
        with pytest.raises(ValidationError):
            raise ValidationError("Test error")

    def test_validation_error_inherits_from_exception(self):
        """Test that ValidationError inherits from Exception."""
        assert issubclass(ValidationError, Exception)

    def test_validation_error_stores_message(self):
        """Test that ValidationError stores and returns error message."""
        message = "This is a test error message"
        try:
            raise ValidationError(message)
        except ValidationError as e:
            assert str(e) == message

    def test_validation_error_with_no_message(self):
        """Test that ValidationError can be raised without a message."""
        with pytest.raises(ValidationError):
            raise ValidationError()


class TestConfigurationError:
    """Tests for ConfigurationError exception."""

    def test_can_raise_configuration_error(self):
        """Test that ConfigurationError can be raised."""
        with pytest.raises(ConfigurationError):
            raise ConfigurationError("Invalid config")

    def test_inherits_from_validation_error(self):
        """Test that ConfigurationError inherits from ValidationError."""
        assert issubclass(ConfigurationError, ValidationError)

    def test_can_be_caught_as_validation_error(self):
        """Test that ConfigurationError can be caught as ValidationError."""
        with pytest.raises(ValidationError):
            raise ConfigurationError("Missing required field")

    def test_stores_message(self):
        """Test that ConfigurationError stores error message."""
        message = "Missing required field: ref_genome_filename"
        try:
            raise ConfigurationError(message)
        except ConfigurationError as e:
            assert str(e) == message


class TestFileNotFoundError:
    """Tests for FileNotFoundError exception."""

    def test_can_raise_file_not_found_error(self):
        """Test that FileNotFoundError can be raised."""
        with pytest.raises(FileNotFoundError):
            raise FileNotFoundError("File not found")

    def test_inherits_from_validation_error(self):
        """Test that FileNotFoundError inherits from ValidationError."""
        assert issubclass(FileNotFoundError, ValidationError)

    def test_can_be_caught_as_validation_error(self):
        """Test that FileNotFoundError can be caught as ValidationError."""
        with pytest.raises(ValidationError):
            raise FileNotFoundError("Genome file doesn't exist")

    def test_stores_message(self):
        """Test that FileNotFoundError stores error message."""
        message = "Genome file doesn't exist: /path/to/genome.fasta"
        try:
            raise FileNotFoundError(message)
        except FileNotFoundError as e:
            assert str(e) == message


class TestFileFormatError:
    """Tests for FileFormatError base exception."""

    def test_can_raise_file_format_error(self):
        """Test that FileFormatError can be raised."""
        with pytest.raises(FileFormatError):
            raise FileFormatError("Invalid format")

    def test_inherits_from_validation_error(self):
        """Test that FileFormatError inherits from ValidationError."""
        assert issubclass(FileFormatError, ValidationError)

    def test_can_be_caught_as_validation_error(self):
        """Test that FileFormatError can be caught as ValidationError."""
        with pytest.raises(ValidationError):
            raise FileFormatError("File format invalid")

    def test_stores_message(self):
        """Test that FileFormatError stores error message."""
        message = "File has invalid format"
        try:
            raise FileFormatError(message)
        except FileFormatError as e:
            assert str(e) == message


class TestFastaFormatError:
    """Tests for FastaFormatError exception."""

    def test_can_raise_fasta_format_error(self):
        """Test that FastaFormatError can be raised."""
        with pytest.raises(FastaFormatError):
            raise FastaFormatError("Invalid FASTA")

    def test_inherits_from_file_format_error(self):
        """Test that FastaFormatError inherits from FileFormatError."""
        assert issubclass(FastaFormatError, FileFormatError)

    def test_inherits_from_validation_error(self):
        """Test that FastaFormatError inherits from ValidationError."""
        assert issubclass(FastaFormatError, ValidationError)

    def test_can_be_caught_as_file_format_error(self):
        """Test that FastaFormatError can be caught as FileFormatError."""
        with pytest.raises(FileFormatError):
            raise FastaFormatError("Missing header line")

    def test_can_be_caught_as_validation_error(self):
        """Test that FastaFormatError can be caught as ValidationError."""
        with pytest.raises(ValidationError):
            raise FastaFormatError("Invalid characters in sequence")

    def test_stores_message(self):
        """Test that FastaFormatError stores error message."""
        message = "Invalid FASTA format at line 42: expected '>' but found 'A'"
        try:
            raise FastaFormatError(message)
        except FastaFormatError as e:
            assert str(e) == message


class TestGenBankFormatError:
    """Tests for GenBankFormatError exception."""

    def test_can_raise_genbank_format_error(self):
        """Test that GenBankFormatError can be raised."""
        with pytest.raises(GenBankFormatError):
            raise GenBankFormatError("Invalid GenBank")

    def test_inherits_from_file_format_error(self):
        """Test that GenBankFormatError inherits from FileFormatError."""
        assert issubclass(GenBankFormatError, FileFormatError)

    def test_can_be_caught_as_file_format_error(self):
        """Test that GenBankFormatError can be caught as FileFormatError."""
        with pytest.raises(FileFormatError):
            raise GenBankFormatError("Missing LOCUS line")

    def test_stores_message(self):
        """Test that GenBankFormatError stores error message."""
        message = "Missing LOCUS line in GenBank file"
        try:
            raise GenBankFormatError(message)
        except GenBankFormatError as e:
            assert str(e) == message


class TestBedFormatError:
    """Tests for BedFormatError exception."""

    def test_can_raise_bed_format_error(self):
        """Test that BedFormatError can be raised."""
        with pytest.raises(BedFormatError):
            raise BedFormatError("Invalid BED")

    def test_inherits_from_file_format_error(self):
        """Test that BedFormatError inherits from FileFormatError."""
        assert issubclass(BedFormatError, FileFormatError)

    def test_can_be_caught_as_file_format_error(self):
        """Test that BedFormatError can be caught as FileFormatError."""
        with pytest.raises(FileFormatError):
            raise BedFormatError("Wrong number of columns")

    def test_stores_message(self):
        """Test that BedFormatError stores error message."""
        message = "Invalid coordinates: start > end at line 15"
        try:
            raise BedFormatError(message)
        except BedFormatError as e:
            assert str(e) == message


class TestGffFormatError:
    """Tests for GffFormatError exception."""

    def test_can_raise_gff_format_error(self):
        """Test that GffFormatError can be raised."""
        with pytest.raises(GffFormatError):
            raise GffFormatError("Invalid GFF")

    def test_inherits_from_file_format_error(self):
        """Test that GffFormatError inherits from FileFormatError."""
        assert issubclass(GffFormatError, FileFormatError)

    def test_can_be_caught_as_file_format_error(self):
        """Test that GffFormatError can be caught as FileFormatError."""
        with pytest.raises(FileFormatError):
            raise GffFormatError("Missing required attributes")

    def test_stores_message(self):
        """Test that GffFormatError stores error message."""
        message = "Wrong number of columns in GFF line 8"
        try:
            raise GffFormatError(message)
        except GffFormatError as e:
            assert str(e) == message


class TestFastqFormatError:
    """Tests for FastqFormatError exception."""

    def test_can_raise_fastq_format_error(self):
        """Test that FastqFormatError can be raised."""
        with pytest.raises(FastqFormatError):
            raise FastqFormatError("Invalid FASTQ")

    def test_inherits_from_file_format_error(self):
        """Test that FastqFormatError inherits from FileFormatError."""
        assert issubclass(FastqFormatError, FileFormatError)

    def test_can_be_caught_as_file_format_error(self):
        """Test that FastqFormatError can be caught as FileFormatError."""
        with pytest.raises(FileFormatError):
            raise FastqFormatError("Missing quality scores")

    def test_stores_message(self):
        """Test that FastqFormatError stores error message."""
        message = "Length mismatch between sequence and quality at read 42"
        try:
            raise FastqFormatError(message)
        except FastqFormatError as e:
            assert str(e) == message


class TestBamFormatError:
    """Tests for BamFormatError exception."""

    def test_can_raise_bam_format_error(self):
        """Test that BamFormatError can be raised."""
        with pytest.raises(BamFormatError):
            raise BamFormatError("Invalid BAM")

    def test_inherits_from_file_format_error(self):
        """Test that BamFormatError inherits from FileFormatError."""
        assert issubclass(BamFormatError, FileFormatError)

    def test_can_be_caught_as_file_format_error(self):
        """Test that BamFormatError can be caught as FileFormatError."""
        with pytest.raises(FileFormatError):
            raise BamFormatError("Corrupted BAM header")

    def test_stores_message(self):
        """Test that BamFormatError stores error message."""
        message = "Invalid BAM magic number"
        try:
            raise BamFormatError(message)
        except BamFormatError as e:
            assert str(e) == message


class TestCompressionError:
    """Tests for CompressionError exception."""

    def test_can_raise_compression_error(self):
        """Test that CompressionError can be raised."""
        with pytest.raises(CompressionError):
            raise CompressionError("Decompression failed")

    def test_inherits_from_validation_error(self):
        """Test that CompressionError inherits from ValidationError."""
        assert issubclass(CompressionError, ValidationError)

    def test_can_be_caught_as_validation_error(self):
        """Test that CompressionError can be caught as ValidationError."""
        with pytest.raises(ValidationError):
            raise CompressionError("Corrupted gzip file")

    def test_stores_message(self):
        """Test that CompressionError stores error message."""
        message = "Corrupted gzip file: /path/to/file.gz"
        try:
            raise CompressionError(message)
        except CompressionError as e:
            assert str(e) == message


class TestGenomeValidationError:
    """Tests for GenomeValidationError exception."""

    def test_can_raise_genome_validation_error(self):
        """Test that GenomeValidationError can be raised."""
        with pytest.raises(GenomeValidationError):
            raise GenomeValidationError("Genome validation failed")

    def test_inherits_from_validation_error(self):
        """Test that GenomeValidationError inherits from ValidationError."""
        assert issubclass(GenomeValidationError, ValidationError)

    def test_can_be_caught_as_validation_error(self):
        """Test that GenomeValidationError can be caught as ValidationError."""
        with pytest.raises(ValidationError):
            raise GenomeValidationError("Duplicate sequence IDs")

    def test_stores_message(self):
        """Test that GenomeValidationError stores error message."""
        message = "Duplicate sequence ID: 'chr1' found at positions 1 and 5"
        try:
            raise GenomeValidationError(message)
        except GenomeValidationError as e:
            assert str(e) == message


class TestFeatureValidationError:
    """Tests for FeatureValidationError exception."""

    def test_can_raise_feature_validation_error(self):
        """Test that FeatureValidationError can be raised."""
        with pytest.raises(FeatureValidationError):
            raise FeatureValidationError("Feature validation failed")

    def test_inherits_from_validation_error(self):
        """Test that FeatureValidationError inherits from ValidationError."""
        assert issubclass(FeatureValidationError, ValidationError)

    def test_can_be_caught_as_validation_error(self):
        """Test that FeatureValidationError can be caught as ValidationError."""
        with pytest.raises(ValidationError):
            raise FeatureValidationError("Overlapping features")

    def test_stores_message(self):
        """Test that FeatureValidationError stores error message."""
        message = "Overlapping features not allowed: 'gene1' overlaps with 'gene2'"
        try:
            raise FeatureValidationError(message)
        except FeatureValidationError as e:
            assert str(e) == message


class TestReadValidationError:
    """Tests for ReadValidationError exception."""

    def test_can_raise_read_validation_error(self):
        """Test that ReadValidationError can be raised."""
        with pytest.raises(ReadValidationError):
            raise ReadValidationError("Read validation failed")

    def test_inherits_from_validation_error(self):
        """Test that ReadValidationError inherits from ValidationError."""
        assert issubclass(ReadValidationError, ValidationError)

    def test_can_be_caught_as_validation_error(self):
        """Test that ReadValidationError can be caught as ValidationError."""
        with pytest.raises(ValidationError):
            raise ReadValidationError("Empty read file")

    def test_stores_message(self):
        """Test that ReadValidationError stores error message."""
        message = "Inconsistent paired-end reads: R1 has 100 reads but R2 has 95"
        try:
            raise ReadValidationError(message)
        except ReadValidationError as e:
            assert str(e) == message


class TestInterFileValidationError:
    """Tests for InterFileValidationError exception."""

    def test_can_raise_inter_file_validation_error(self):
        """Test that InterFileValidationError can be raised."""
        with pytest.raises(InterFileValidationError):
            raise InterFileValidationError("Inter-file validation failed")

    def test_inherits_from_validation_error(self):
        """Test that InterFileValidationError inherits from ValidationError."""
        assert issubclass(InterFileValidationError, ValidationError)

    def test_can_be_caught_as_validation_error(self):
        """Test that InterFileValidationError can be caught as ValidationError."""
        with pytest.raises(ValidationError):
            raise InterFileValidationError("Feature coordinates outside genome bounds")

    def test_stores_message(self):
        """Test that InterFileValidationError stores error message."""
        message = "Feature 'gene1' at position 5000 exceeds genome length of 4500"
        try:
            raise InterFileValidationError(message)
        except InterFileValidationError as e:
            assert str(e) == message


class TestExceptionHierarchy:
    """Tests for the overall exception hierarchy."""

    def test_all_exceptions_inherit_from_validation_error(self):
        """Test that all custom exceptions inherit from ValidationError."""
        exceptions = [
            ConfigurationError,
            FileNotFoundError,
            FileFormatError,
            FastaFormatError,
            GenBankFormatError,
            BedFormatError,
            GffFormatError,
            FastqFormatError,
            BamFormatError,
            CompressionError,
            GenomeValidationError,
            FeatureValidationError,
            ReadValidationError,
            InterFileValidationError,
        ]

        for exc in exceptions:
            assert issubclass(exc, ValidationError), \
                f"{exc.__name__} does not inherit from ValidationError"

    def test_file_format_exceptions_inherit_from_file_format_error(self):
        """Test that format-specific exceptions inherit from FileFormatError."""
        format_exceptions = [
            FastaFormatError,
            GenBankFormatError,
            BedFormatError,
            GffFormatError,
            FastqFormatError,
            BamFormatError,
        ]

        for exc in format_exceptions:
            assert issubclass(exc, FileFormatError), \
                f"{exc.__name__} does not inherit from FileFormatError"

    def test_catch_all_file_format_errors(self):
        """Test that all file format errors can be caught with FileFormatError."""
        format_exceptions = [
            FastaFormatError("test"),
            GenBankFormatError("test"),
            BedFormatError("test"),
            GffFormatError("test"),
            FastqFormatError("test"),
            BamFormatError("test"),
        ]

        for exc in format_exceptions:
            with pytest.raises(FileFormatError):
                raise exc

    def test_catch_all_exceptions_with_validation_error(self):
        """Test that all exceptions can be caught with ValidationError."""
        all_exceptions = [
            ConfigurationError("test"),
            FileNotFoundError("test"),
            FileFormatError("test"),
            FastaFormatError("test"),
            GenBankFormatError("test"),
            BedFormatError("test"),
            GffFormatError("test"),
            FastqFormatError("test"),
            BamFormatError("test"),
            CompressionError("test"),
            GenomeValidationError("test"),
            FeatureValidationError("test"),
            ReadValidationError("test"),
            InterFileValidationError("test"),
        ]

        for exc in all_exceptions:
            with pytest.raises(ValidationError):
                raise exc


class TestExceptionUsagePatterns:
    """Tests for common usage patterns with the exceptions."""

    def test_try_except_specific_then_general(self):
        """Test catching specific exception first, then general."""
        def raise_fasta_error():
            raise FastaFormatError("Invalid FASTA")

        caught_specific = False
        try:
            raise_fasta_error()
        except FastaFormatError:
            caught_specific = True
        except FileFormatError:
            pass
        except ValidationError:
            pass

        assert caught_specific

    def test_try_except_with_multiple_exception_types(self):
        """Test catching multiple specific exception types."""
        def raise_config_error():
            raise ConfigurationError("Invalid config")

        caught = False
        try:
            raise_config_error()
        except (ConfigurationError, FileNotFoundError):
            caught = True

        assert caught

    def test_exception_with_formatted_message(self):
        """Test exceptions with formatted messages."""
        filename = "/path/to/file.fasta"
        line_num = 42
        message = f"Invalid FASTA format at line {line_num} in {filename}"

        try:
            raise FastaFormatError(message)
        except FastaFormatError as e:
            assert str(e) == message
            assert str(line_num) in str(e)
            assert filename in str(e)

    def test_re_raising_exception(self):
        """Test re-raising an exception after catching it."""
        with pytest.raises(ValidationError):
            try:
                raise ConfigurationError("Test error")
            except ValidationError:
                # Log or process the error
                raise  # Re-raise the same exception

    def test_exception_chaining(self):
        """Test exception chaining with 'from' clause."""
        try:
            try:
                raise ValueError("Original error")
            except ValueError as e:
                raise ConfigurationError("Configuration failed") from e
        except ConfigurationError as e:
            assert e.__cause__ is not None
            assert isinstance(e.__cause__, ValueError)
            assert str(e.__cause__) == "Original error"
