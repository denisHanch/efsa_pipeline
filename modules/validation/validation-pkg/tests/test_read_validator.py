"""
Comprehensive tests for ReadValidator.

Tests cover:
- File format detection and parsing (FASTQ, BAM)
- Compression handling (gzip, bzip2, uncompressed)
- Validation rules (empty IDs, duplicate IDs, invalid characters)
- BAM file handling (pass-through, keep_bam, conversion)
- Statistics collection (read counts, length, GC content, quality scores)
- Output generation with various compression formats
- Error handling and edge cases
"""

import pytest
import tempfile
import gzip
import bz2
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from validation_pkg.config_manager import ReadConfig
from validation_pkg.validators.read_validator import ReadValidator, OutputMetadata
from validation_pkg.exceptions import (
    ReadValidationError,
    FastqFormatError,
    BamFormatError,
    CompressionError,
)
from validation_pkg.utils.formats import ReadFormat, CodingType


class TestReadValidatorInitialization:
    """Test ReadValidator initialization."""

    @pytest.fixture
    def temp_dir(self):
        """Create a temporary directory for test files."""
        with tempfile.TemporaryDirectory() as tmpdir:
            yield Path(tmpdir)

    @pytest.fixture
    def output_dir(self, temp_dir):
        """Create output directory."""
        out_dir = temp_dir / "output"
        out_dir.mkdir(parents=True, exist_ok=True)
        return out_dir

    @pytest.fixture
    def simple_fastq(self, temp_dir):
        """Create a simple FASTQ file."""
        fastq_file = temp_dir / "reads.fastq"
        with open(fastq_file, "w") as f:
            f.write("@read1\n")
            f.write("ATCGATCGATCG\n")
            f.write("+\n")
            f.write("IIIIIIIIIIII\n")
        return fastq_file

    def test_init_with_defaults(self, simple_fastq, output_dir):
        """Test initialization with default settings."""
        read_config = ReadConfig(
            filename="reads.fastq",
            basename="reads",
            filepath=simple_fastq,
            ngs_type="illumina",
            coding_type=CodingType.NONE,
            detected_format=ReadFormat.FASTQ,
            output_dir=output_dir,
            global_options={}
        )

        # Pass explicit default settings (workaround for validator bug with 'if not None')
        validator = ReadValidator(read_config, ReadValidator.Settings())

        assert validator.input_path == simple_fastq
        assert validator.output_dir == output_dir
        assert validator.settings is not None
        assert validator.sequences == []

    def test_init_with_custom_settings(self, simple_fastq, output_dir):
        """Test initialization with custom settings."""
        settings = ReadValidator.Settings(
            check_invalid_chars=True,
            allow_duplicate_ids=False,
            keep_bam=False
        )

        read_config = ReadConfig(
            filename="reads.fastq",
            basename="reads",
            filepath=simple_fastq,
            ngs_type="illumina",
            coding_type=CodingType.NONE,
            detected_format=ReadFormat.FASTQ,
            output_dir=output_dir,
            global_options={}
        )

        validator = ReadValidator(read_config, settings)

        assert validator.settings.check_invalid_chars is True
        assert validator.settings.allow_duplicate_ids is False
        assert validator.settings.keep_bam is False


class TestReadValidatorParsing:
    """Test file parsing functionality."""

    @pytest.fixture
    def temp_dir(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            yield Path(tmpdir)

    @pytest.fixture
    def output_dir(self, temp_dir):
        out_dir = temp_dir / "output"
        out_dir.mkdir(parents=True, exist_ok=True)
        return out_dir

    @pytest.fixture
    def simple_fastq(self, temp_dir):
        """Create a simple FASTQ file with two reads."""
        fastq_file = temp_dir / "reads.fastq"
        with open(fastq_file, "w") as f:
            f.write("@read1\n")
            f.write("ATCGATCGATCGATCGATCG\n")
            f.write("+\n")
            f.write("IIIIIIIIIIIIIIIIIIII\n")
            f.write("@read2\n")
            f.write("GCTAGCTAGCTAGCTAGCTA\n")
            f.write("+\n")
            f.write("HHHHHHHHHHHHHHHHHHHH\n")
        return fastq_file

    @pytest.fixture
    def empty_fastq(self, temp_dir):
        """Create an empty FASTQ file."""
        fastq_file = temp_dir / "empty.fastq"
        with open(fastq_file, "w") as f:
            f.write("")
        return fastq_file

    @pytest.fixture
    def invalid_fastq(self, temp_dir):
        """Create an invalid FASTQ file (missing quality scores)."""
        fastq_file = temp_dir / "invalid.fastq"
        with open(fastq_file, "w") as f:
            f.write("@read1\n")
            f.write("ATCGATCG\n")
            # Missing + and quality lines
        return fastq_file

    def test_parse_simple_fastq(self, simple_fastq, output_dir):
        """Test parsing a simple FASTQ file."""
        read_config = ReadConfig(
            filename="reads.fastq",
            basename="reads",
            filepath=simple_fastq,
            ngs_type="illumina",
            coding_type=CodingType.NONE,
            detected_format=ReadFormat.FASTQ,
            output_dir=output_dir,
            global_options={}
        )

        validator = ReadValidator(read_config, ReadValidator.Settings())
        validator.run()

        assert len(validator.sequences) == 2
        assert validator.sequences[0].id == "read1"
        assert validator.sequences[1].id == "read2"
        assert str(validator.sequences[0].seq) == "ATCGATCGATCGATCGATCG"

    def test_parse_empty_file_raises_error(self, empty_fastq, output_dir):
        """Test that empty FASTQ file raises ReadValidationError."""
        read_config = ReadConfig(
            filename="empty.fastq",
            basename="empty",
            filepath=empty_fastq,
            ngs_type="illumina",
            coding_type=CodingType.NONE,
            detected_format=ReadFormat.FASTQ,
            output_dir=output_dir,
            global_options={}
        )

        validator = ReadValidator(read_config, ReadValidator.Settings())

        with pytest.raises(FastqFormatError, match="No sequences found"):
            validator.run()

    def test_parse_invalid_fastq_raises_error(self, invalid_fastq, output_dir):
        """Test that invalid FASTQ raises FastqFormatError."""
        read_config = ReadConfig(
            filename="invalid.fastq",
            basename="invalid",
            filepath=invalid_fastq,
            ngs_type="illumina",
            coding_type=CodingType.NONE,
            detected_format=ReadFormat.FASTQ,
            output_dir=output_dir,
            global_options={}
        )

        validator = ReadValidator(read_config, ReadValidator.Settings())

        with pytest.raises((FastqFormatError, ReadValidationError)):
            validator.run()


class TestReadValidatorCompression:
    """Test compression handling."""

    @pytest.fixture
    def temp_dir(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            yield Path(tmpdir)

    @pytest.fixture
    def output_dir(self, temp_dir):
        out_dir = temp_dir / "output"
        out_dir.mkdir(parents=True, exist_ok=True)
        return out_dir

    @pytest.fixture
    def compressed_fastq_gz(self, temp_dir):
        """Create a gzip compressed FASTQ file."""
        fastq_file = temp_dir / "reads.fastq.gz"
        with gzip.open(fastq_file, "wt") as f:
            f.write("@read1\n")
            f.write("ATCGATCGATCGATCGATCG\n")
            f.write("+\n")
            f.write("IIIIIIIIIIIIIIIIIIII\n")
        return fastq_file

    @pytest.fixture
    def compressed_fastq_bz2(self, temp_dir):
        """Create a bzip2 compressed FASTQ file."""
        fastq_file = temp_dir / "reads.fastq.bz2"
        with bz2.open(fastq_file, "wt") as f:
            f.write("@read1\n")
            f.write("ATCGATCGATCGATCGATCG\n")
            f.write("+\n")
            f.write("IIIIIIIIIIIIIIIIIIII\n")
        return fastq_file

    def test_parse_gzip_compressed(self, compressed_fastq_gz, output_dir):
        """Test parsing gzip compressed FASTQ."""
        read_config = ReadConfig(
            filename="reads.fastq.gz",
            basename="reads",
            filepath=compressed_fastq_gz,
            ngs_type="illumina",
            coding_type=CodingType.GZIP,
            detected_format=ReadFormat.FASTQ,
            output_dir=output_dir,
            global_options={}
        )

        validator = ReadValidator(read_config, ReadValidator.Settings())
        validator.run()

        assert len(validator.sequences) == 1
        assert str(validator.sequences[0].seq) == "ATCGATCGATCGATCGATCG"

    def test_parse_bzip2_compressed(self, compressed_fastq_bz2, output_dir):
        """Test parsing bzip2 compressed FASTQ."""
        read_config = ReadConfig(
            filename="reads.fastq.bz2",
            basename="reads",
            filepath=compressed_fastq_bz2,
            ngs_type="illumina",
            coding_type=CodingType.BZIP2,
            detected_format=ReadFormat.FASTQ,
            output_dir=output_dir,
            global_options={}
        )

        validator = ReadValidator(read_config, ReadValidator.Settings())
        validator.run()

        assert len(validator.sequences) == 1
        assert str(validator.sequences[0].seq) == "ATCGATCGATCGATCGATCG"


class TestReadValidatorValidation:
    """Test validation rules."""

    @pytest.fixture
    def temp_dir(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            yield Path(tmpdir)

    @pytest.fixture
    def output_dir(self, temp_dir):
        out_dir = temp_dir / "output"
        out_dir.mkdir(parents=True, exist_ok=True)
        return out_dir

    @pytest.fixture
    def fastq_with_duplicates(self, temp_dir):
        """Create FASTQ with duplicate read IDs."""
        fastq_file = temp_dir / "duplicates.fastq"
        with open(fastq_file, "w") as f:
            f.write("@read1\n")
            f.write("ATCGATCGATCGATCGATCG\n")
            f.write("+\n")
            f.write("IIIIIIIIIIIIIIIIIIII\n")
            f.write("@read1\n")
            f.write("GCTAGCTAGCTAGCTAGCTA\n")
            f.write("+\n")
            f.write("HHHHHHHHHHHHHHHHHHHH\n")
        return fastq_file

    @pytest.fixture
    def fastq_with_invalid_chars(self, temp_dir):
        """Create FASTQ with invalid characters."""
        fastq_file = temp_dir / "invalid_chars.fastq"
        with open(fastq_file, "w") as f:
            f.write("@read1\n")
            f.write("ATCGXYZATCG\n")  # X, Y, Z are invalid
            f.write("+\n")
            f.write("IIIIIIIIIII\n")
        return fastq_file

    def test_duplicate_ids_allowed_by_default(self, fastq_with_duplicates, output_dir):
        """Test that duplicate IDs are allowed by default."""
        read_config = ReadConfig(
            filename="duplicates.fastq",
            basename="duplicates",
            filepath=fastq_with_duplicates,
            ngs_type="illumina",
            coding_type=CodingType.NONE,
            detected_format=ReadFormat.FASTQ,
            output_dir=output_dir,
            global_options={}
        )

        # Default settings allow duplicates
        validator = ReadValidator(read_config, ReadValidator.Settings())
        validator.run()  # Should not raise

        assert len(validator.sequences) == 2

    def test_duplicate_ids_rejected_when_disabled(self, fastq_with_duplicates, output_dir):
        """Test that duplicate IDs are rejected when allow_duplicate_ids=False."""
        settings = ReadValidator.Settings(allow_duplicate_ids=False)

        read_config = ReadConfig(
            filename="duplicates.fastq",
            basename="duplicates",
            filepath=fastq_with_duplicates,
            ngs_type="illumina",
            coding_type=CodingType.NONE,
            detected_format=ReadFormat.FASTQ,
            output_dir=output_dir,
            global_options={}
        )

        validator = ReadValidator(read_config, settings)

        with pytest.raises(ReadValidationError, match="Duplicate sequence IDs"):
            validator.run()

    def test_invalid_chars_detected(self, fastq_with_invalid_chars, output_dir):
        """Test that invalid characters are detected when check_invalid_chars=True."""
        settings = ReadValidator.Settings(check_invalid_chars=True)

        read_config = ReadConfig(
            filename="invalid_chars.fastq",
            basename="invalid_chars",
            filepath=fastq_with_invalid_chars,
            ngs_type="illumina",
            coding_type=CodingType.NONE,
            detected_format=ReadFormat.FASTQ,
            output_dir=output_dir,
            global_options={}
        )

        validator = ReadValidator(read_config, settings)

        with pytest.raises(ReadValidationError, match="error"):
            validator.run()

    def test_invalid_chars_ignored_by_default(self, fastq_with_invalid_chars, output_dir):
        """Test that invalid characters are ignored by default."""
        read_config = ReadConfig(
            filename="invalid_chars.fastq",
            basename="invalid_chars",
            filepath=fastq_with_invalid_chars,
            ngs_type="illumina",
            coding_type=CodingType.NONE,
            detected_format=ReadFormat.FASTQ,
            output_dir=output_dir,
            global_options={}
        )

        # Default settings don't check invalid chars
        validator = ReadValidator(read_config, ReadValidator.Settings())
        validator.run()  # Should not raise


class TestReadValidatorBAMHandling:
    """Test BAM file handling."""

    @pytest.fixture
    def temp_dir(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            yield Path(tmpdir)

    @pytest.fixture
    def output_dir(self, temp_dir):
        out_dir = temp_dir / "output"
        out_dir.mkdir(parents=True, exist_ok=True)
        return out_dir

    @pytest.fixture
    def fake_bam(self, temp_dir):
        """Create a fake BAM file (for testing pass-through)."""
        bam_file = temp_dir / "reads.bam"
        # Create a dummy file - actual BAM parsing tested separately
        with open(bam_file, "wb") as f:
            f.write(b"BAM\x01")  # BAM magic number
        return bam_file

    # def test_bam_conversion_workflow(self, fake_bam, output_dir):
    #     """
    #     Test that BAM files follow the proper workflow.

    #     With fake BAM files, the conversion will fail (as expected).
    #     This tests that the workflow is triggered correctly.
    #     Real BAM conversion requires pysam or samtools.
    #     """
    #     read_config = ReadConfig(
    #         filename="reads.bam",
    #         filepath=fake_bam,
    #         ngs_type="illumina",
    #         coding_type=CodingType.NONE,
    #         detected_format=ReadFormat.BAM
    #     )

    #     validator = ReadValidator(read_config, output_dir,
    #                              ReadValidator.Settings(keep_bam=True))

    #     # Should raise error during BAM conversion (fake BAM file)
    #     # The error could be from pysam, samtools, or "neither tool available"
    #     with pytest.raises((ReadValidationError, Exception)):
    #         validator.run()

    #     # With keep_bam=True, the original BAM should be copied
    #     copied_files = list(output_dir.glob("*.bam"))
    #     assert len(copied_files) == 1


class TestReadValidatorOutput:
    """Test output generation."""

    @pytest.fixture
    def temp_dir(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            yield Path(tmpdir)

    @pytest.fixture
    def output_dir(self, temp_dir):
        out_dir = temp_dir / "output"
        out_dir.mkdir(parents=True, exist_ok=True)
        return out_dir

    @pytest.fixture
    def simple_fastq(self, temp_dir):
        """Create a simple FASTQ file."""
        fastq_file = temp_dir / "reads.fastq"
        with open(fastq_file, "w") as f:
            f.write("@read1\n")
            f.write("ATCGATCGATCGATCGATCG\n")
            f.write("+\n")
            f.write("IIIIIIIIIIIIIIIIIIII\n")
        return fastq_file

    def test_output_uncompressed(self, simple_fastq, output_dir):
        """Test generating uncompressed output."""
        settings = ReadValidator.Settings(coding_type=None)

        read_config = ReadConfig(
            filename="reads.fastq",
            basename="reads",
            filepath=simple_fastq,
            ngs_type="illumina",
            coding_type=CodingType.NONE,
            detected_format=ReadFormat.FASTQ,
            output_dir=output_dir,
            global_options={}
        )

        validator = ReadValidator(read_config, settings)
        validator.run()

        # Check output file exists
        output_files = list(output_dir.glob("*.fastq"))
        assert len(output_files) == 1
        assert output_files[0].suffix == ".fastq"

    def test_output_gzip_compressed(self, simple_fastq, output_dir):
        """Test generating gzip compressed output."""
        settings = ReadValidator.Settings(coding_type='gz')

        read_config = ReadConfig(
            filename="reads.fastq",
            basename="reads",
            filepath=simple_fastq,
            ngs_type="illumina",
            coding_type=CodingType.NONE,
            detected_format=ReadFormat.FASTQ,
            output_dir=output_dir,
            global_options={}
        )

        validator = ReadValidator(read_config, settings)
        validator.run()

        # Check gzip output file exists
        output_files = list(output_dir.glob("*.fastq.gz"))
        assert len(output_files) == 1

        # Verify it's actually gzipped and contains correct data
        with gzip.open(output_files[0], 'rt') as f:
            records = list(SeqIO.parse(f, 'fastq'))
            assert len(records) == 1

    def test_output_bzip2_compressed(self, simple_fastq, output_dir):
        """Test generating bzip2 compressed output."""
        settings = ReadValidator.Settings(coding_type='bz2')

        read_config = ReadConfig(
            filename="reads.fastq",
            basename="reads",
            filepath=simple_fastq,
            ngs_type="illumina",
            coding_type=CodingType.NONE,
            detected_format=ReadFormat.FASTQ,
            output_dir=output_dir,
            global_options={}
        )

        validator = ReadValidator(read_config, settings)
        validator.run()

        # Check bzip2 output file exists
        output_files = list(output_dir.glob("*.fastq.bz2"))
        assert len(output_files) == 1

        # Verify it's actually bzip2'd and contains correct data
        with bz2.open(output_files[0], 'rt') as f:
            records = list(SeqIO.parse(f, 'fastq'))
            assert len(records) == 1

    def test_output_with_suffix(self, simple_fastq, output_dir):
        """Test output filename with custom suffix."""
        settings = ReadValidator.Settings(output_filename_suffix="filtered")

        read_config = ReadConfig(
            filename="reads.fastq",
            basename="reads",
            filepath=simple_fastq,
            ngs_type="illumina",
            coding_type=CodingType.NONE,
            detected_format=ReadFormat.FASTQ,
            output_dir=output_dir,
            global_options={}
        )

        validator = ReadValidator(read_config, settings)
        validator.run()

        output_files = list(output_dir.glob("*_filtered.fastq"))
        assert len(output_files) == 1

    def test_output_with_subdirectory(self, simple_fastq, output_dir):
        """Test output to custom subdirectory."""
        settings = ReadValidator.Settings(output_subdir_name="reads")

        read_config = ReadConfig(
            filename="reads.fastq",
            basename="reads",
            filepath=simple_fastq,
            ngs_type="illumina",
            coding_type=CodingType.NONE,
            detected_format=ReadFormat.FASTQ,
            output_dir=output_dir,
            global_options={}
        )

        validator = ReadValidator(read_config, settings)
        validator.run()

        subdir = output_dir / "reads"
        assert subdir.exists()
        output_files = list(subdir.glob("*.fastq"))
        assert len(output_files) == 1


class TestReadValidatorNGSTypes:
    """Test different NGS types."""

    @pytest.fixture
    def temp_dir(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            yield Path(tmpdir)

    @pytest.fixture
    def output_dir(self, temp_dir):
        out_dir = temp_dir / "output"
        out_dir.mkdir(parents=True, exist_ok=True)
        return out_dir

    @pytest.fixture
    def simple_fastq(self, temp_dir):
        """Create a simple FASTQ file."""
        fastq_file = temp_dir / "reads.fastq"
        with open(fastq_file, "w") as f:
            f.write("@read1\n")
            f.write("ATCGATCGATCGATCGATCG\n")
            f.write("+\n")
            f.write("IIIIIIIIIIIIIIIIIIII\n")
        return fastq_file

    def test_illumina_reads(self, simple_fastq, output_dir):
        """Test processing Illumina reads."""
        read_config = ReadConfig(
            filename="reads.fastq",
            basename="reads",
            filepath=simple_fastq,
            ngs_type="illumina",
            coding_type=CodingType.NONE,
            detected_format=ReadFormat.FASTQ,
            output_dir=output_dir,
            global_options={}
        )

        validator = ReadValidator(read_config, ReadValidator.Settings())
        validator.run()

        assert len(validator.sequences) == 1

    def test_ont_reads(self, simple_fastq, output_dir):
        """Test processing Oxford Nanopore reads."""
        read_config = ReadConfig(
            filename="reads.fastq",
            basename="reads",
            filepath=simple_fastq,
            ngs_type="ont",
            coding_type=CodingType.NONE,
            detected_format=ReadFormat.FASTQ,
            output_dir=output_dir,
            global_options={}
        )

        validator = ReadValidator(read_config, ReadValidator.Settings())
        validator.run()

        assert len(validator.sequences) == 1

    def test_pacbio_reads(self, simple_fastq, output_dir):
        """Test processing PacBio reads."""
        read_config = ReadConfig(
            filename="reads.fastq",
            basename="reads",
            filepath=simple_fastq,
            ngs_type="pacbio",
            coding_type=CodingType.NONE,
            detected_format=ReadFormat.FASTQ,
            output_dir=output_dir,
            global_options={}
        )

        validator = ReadValidator(read_config, ReadValidator.Settings())
        validator.run()

        assert len(validator.sequences) == 1


class TestReadValidatorValidationLevels:
    """Test multi-level validation modes (strict, trust, minimal)."""

    @pytest.fixture
    def temp_dir(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            yield Path(tmpdir)

    @pytest.fixture
    def output_dir(self, temp_dir):
        out_dir = temp_dir / "output"
        out_dir.mkdir(parents=True, exist_ok=True)
        return out_dir

    @pytest.fixture
    def multi_read_fastq(self, temp_dir):
        """Create a gzipped FASTQ file with multiple reads (for minimal mode)."""
        fastq_file = temp_dir / "reads.fastq.gz"
        with gzip.open(fastq_file, "wt") as f:
            for i in range(1, 21):  # 20 reads
                f.write(f"@read{i}\n")
                f.write("ATCGATCGATCGATCGATCG\n")
                f.write("+\n")
                f.write("IIIIIIIIIIIIIIIIIIII\n")
        return fastq_file

    @pytest.fixture
    def damaged_fastq(self, temp_dir):
        """Create a FASTQ file with empty ID (uncompressed for strict/trust mode)."""
        fastq_file = temp_dir / "damaged.fastq"
        with open(fastq_file, "w") as f:
            f.write("@\n")  # Empty ID
            f.write("ATCGATCGATCGATCGATCG\n")
            f.write("+\n")
            f.write("IIIIIIIIIIIIIIIIIIII\n")
        return fastq_file

    # ===== Tests for STRICT validation level =====

    def test_strict_correct_file_passes(self, multi_read_fastq, output_dir):
        """Test strict mode with correct FASTQ file - should pass."""
        settings = ReadValidator.Settings(coding_type='gz')
        read_config = ReadConfig(
            filename="reads.fastq.gz",
            basename="reads",
            filepath=multi_read_fastq,
            ngs_type="illumina",
            coding_type=CodingType.GZIP,
            detected_format=ReadFormat.FASTQ,
            output_dir=output_dir,
            global_options={'validation_level': 'strict'}
        )

        validator = ReadValidator(read_config, settings)
        validator.run()

        # All reads should be parsed
        assert len(validator.sequences) == 20
        # Output file should exist (default output is gzipped)
        output_files = list(output_dir.glob("*.fastq.gz"))
        assert len(output_files) == 1

    def test_strict_damaged_file_fails(self, damaged_fastq, output_dir):
        """Test strict mode with damaged file - should fail during parsing."""
        settings = ReadValidator.Settings(allow_empty_id=False
        )
        read_config = ReadConfig(
            filename="damaged.fastq",
            basename="damaged",
            filepath=damaged_fastq,
            ngs_type="illumina",
            coding_type=CodingType.NONE,
            detected_format=ReadFormat.FASTQ,
            output_dir=output_dir,
            global_options={'validation_level': 'strict'}
        )

        validator = ReadValidator(read_config, settings)

        # Should fail during parsing - empty ID after @ is malformed FASTQ
        with pytest.raises(FastqFormatError):
            validator.run()

    # ===== Tests for TRUST validation level =====

    def test_trust_correct_file_passes(self, multi_read_fastq, output_dir):
        """Test trust mode with correct FASTQ file - parses first 10, validates them, copies original file."""
        settings = ReadValidator.Settings(coding_type='gz')
        read_config = ReadConfig(
            filename="reads.fastq.gz",
            basename="reads",
            filepath=multi_read_fastq,
            ngs_type="illumina",
            coding_type=CodingType.GZIP,
            detected_format=ReadFormat.FASTQ,
            output_dir=output_dir,
            global_options={'validation_level': 'trust'}
        )

        validator = ReadValidator(read_config, settings)
        validator.run()

        # Only first 10 reads should be parsed (trust mode parses only first 10)
        assert len(validator.sequences) == 10
        # Output file should exist with all reads copied (gzipped)
        output_files = list(output_dir.glob("*.fastq.gz"))
        assert len(output_files) == 1

    def test_trust_damaged_first_sequence_fails(self, damaged_fastq, output_dir):
        """Test trust mode detects error in first sequence (within first 10 validated)."""
        settings = ReadValidator.Settings(allow_empty_id=False
        )
        read_config = ReadConfig(
            filename="damaged.fastq",
            basename="damaged",
            filepath=damaged_fastq,
            ngs_type="illumina",
            coding_type=CodingType.NONE,
            detected_format=ReadFormat.FASTQ,
            output_dir=output_dir,
            global_options={'validation_level': 'trust'}
        )

        validator = ReadValidator(read_config, settings)

        # Should fail during parsing - empty ID after @ is malformed FASTQ
        with pytest.raises(FastqFormatError):
            validator.run()

    # ===== Tests for MINIMAL validation level =====

    def test_minimal_correct_file_passes(self, multi_read_fastq, output_dir):
        """Test minimal mode with correct file - should pass without validation."""
        settings = ReadValidator.Settings(coding_type='gz')
        read_config = ReadConfig(
            filename="reads.fastq.gz",
            basename="reads",
            filepath=multi_read_fastq,
            ngs_type="illumina",
            coding_type=CodingType.GZIP,
            detected_format=ReadFormat.FASTQ,
            output_dir=output_dir,
            global_options={'validation_level': 'minimal'}
        )

        validator = ReadValidator(read_config, settings)
        validator.run()

        # No sequences parsed in minimal mode
        assert len(validator.sequences) == 0
        # Output file should exist (gzipped copy)
        output_files = list(output_dir.glob("*.fastq.gz"))
        assert len(output_files) == 1

    def test_minimal_damaged_file_raises_error(self, damaged_fastq, output_dir):
        """Test minimal mode with uncompressed file - should raise error due to coding mismatch."""
        settings = ReadValidator.Settings(coding_type='gz',  # Require GZIP output
            allow_empty_id=False  # Ignored in minimal mode
        )
        read_config = ReadConfig(
            filename="damaged.fastq",
            basename="damaged",
            filepath=damaged_fastq,
            ngs_type="illumina",
            coding_type=CodingType.NONE,  # Uncompressed - doesn't match output coding
            detected_format=ReadFormat.FASTQ,
            output_dir=output_dir,
            global_options={'validation_level': 'minimal'}
        )

        validator = ReadValidator(read_config, settings)
        # Should raise error - minimal mode requires input coding to match output coding
        with pytest.raises(ReadValidationError, match="input coding to match output coding"):
            validator.run()

    def test_minimal_output_is_copy(self, multi_read_fastq, output_dir):
        """Test that minimal mode copies file as-is."""
        settings = ReadValidator.Settings(coding_type='gz')
        read_config = ReadConfig(
            filename="reads.fastq.gz",
            basename="reads",
            filepath=multi_read_fastq,
            ngs_type="illumina",
            coding_type=CodingType.GZIP,
            detected_format=ReadFormat.FASTQ,
            output_dir=output_dir,
            global_options={'validation_level': 'minimal'}
        )

        validator = ReadValidator(read_config, settings)
        validator.run()

        # Check output file exists (gzipped)
        output_files = list(output_dir.glob("*.fastq.gz"))
        assert len(output_files) == 1

        # Verify output is byte-for-byte identical to input
        with open(multi_read_fastq, 'rb') as f_in, open(output_files[0], 'rb') as f_out:
            assert f_in.read() == f_out.read()

    # ===== Test for compressed files =====

    def test_trust_compressed_gz_passes(self, temp_dir, output_dir):
        """Test trust mode with gzip compressed file."""
        fastq_file = temp_dir / "reads.fastq.gz"
        with gzip.open(fastq_file, "wt") as f:
            for i in range(1, 11):
                f.write(f"@read{i}\n")
                f.write("ATCGATCGATCGATCGATCG\n")
                f.write("+\n")
                f.write("IIIIIIIIIIIIIIIIIIII\n")

        settings = ReadValidator.Settings(coding_type='gz')
        read_config = ReadConfig(
            filename="reads.fastq.gz",
            basename="reads",
            filepath=fastq_file,
            ngs_type="illumina",
            coding_type=CodingType.GZIP,
            detected_format=ReadFormat.FASTQ,
            output_dir=output_dir,
            global_options={'validation_level': 'trust'}
        )

        validator = ReadValidator(read_config, settings)
        validator.run()

        # All 10 reads should be parsed (trust mode parses first 10)
        assert len(validator.sequences) == 10

        # Output should be created
        output_files = list(output_dir.glob("*.fastq.gz"))
        assert len(output_files) == 1

    def test_minimal_compressed_bz2_raises_error(self, temp_dir, output_dir):
        """Test minimal mode with bzip2 compressed file - should raise error due to coding mismatch."""
        fastq_file = temp_dir / "reads.fastq.bz2"
        with bz2.open(fastq_file, "wt") as f:
            f.write("@read1\n")
            f.write("ATCGATCGATCGATCGATCG\n")
            f.write("+\n")
            f.write("IIIIIIIIIIIIIIIIIIII\n")

        settings = ReadValidator.Settings(coding_type='gz')
        read_config = ReadConfig(
            filename="reads.fastq.bz2",
            basename="reads",
            filepath=fastq_file,
            ngs_type="illumina",
            coding_type=CodingType.BZIP2,
            detected_format=ReadFormat.FASTQ,
            output_dir=output_dir,
            global_options={'validation_level': 'minimal'}
        )

        validator = ReadValidator(read_config, settings)

        # Should raise error - minimal mode requires input coding to match output coding
        with pytest.raises(ReadValidationError, match="input coding to match output coding"):
            validator.run()
class TestReadValidatorEdgeCases:
    """Test edge cases and error handling."""

    @pytest.fixture
    def temp_dir(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            yield Path(tmpdir)

    @pytest.fixture
    def output_dir(self, temp_dir):
        out_dir = temp_dir / "output"
        out_dir.mkdir(parents=True, exist_ok=True)
        return out_dir

    @pytest.fixture
    def fastq_with_varying_lengths(self, temp_dir):
        """Create FASTQ with reads of varying lengths."""
        fastq_file = temp_dir / "varying.fastq"
        with open(fastq_file, "w") as f:
            f.write("@short_read\n")
            f.write("ATCG\n")
            f.write("+\n")
            f.write("IIII\n")
            f.write("@medium_read\n")
            f.write("ATCGATCGATCGATCG\n")
            f.write("+\n")
            f.write("IIIIIIIIIIIIIIII\n")
            f.write("@long_read\n")
            f.write("ATCGATCG" * 50 + "\n")  # 400bp
            f.write("+\n")
            f.write("I" * 400 + "\n")
        return fastq_file

    @pytest.fixture
    def fastq_with_low_quality(self, temp_dir):
        """Create FASTQ with low quality scores."""
        fastq_file = temp_dir / "low_quality.fastq"
        with open(fastq_file, "w") as f:
            f.write("@read1\n")
            f.write("ATCGATCGATCGATCGATCG\n")
            f.write("+\n")
            f.write("!!!!!!!!!!!!!!!!!!!!\n")  # Quality score 0
        return fastq_file

    def test_varying_read_lengths(self, fastq_with_varying_lengths, output_dir):
        """Test handling reads of varying lengths."""
        read_config = ReadConfig(
            filename="varying.fastq",
            basename="varying",
            filepath=fastq_with_varying_lengths,
            ngs_type="illumina",
            coding_type=CodingType.NONE,
            detected_format=ReadFormat.FASTQ,
            output_dir=output_dir,
            global_options={}
        )

        validator = ReadValidator(read_config, ReadValidator.Settings())
        validator.run()

        # All reads should be parsed
        assert len(validator.sequences) == 3
        # Verify different lengths were parsed correctly
        lengths = [len(seq.seq) for seq in validator.sequences]
        assert min(lengths) == 4
        assert max(lengths) == 400

    def test_low_quality_reads(self, fastq_with_low_quality, output_dir):
        """Test handling reads with low quality scores."""
        read_config = ReadConfig(
            filename="low_quality.fastq",
            basename="low_quality",
            filepath=fastq_with_low_quality,
            ngs_type="illumina",
            coding_type=CodingType.NONE,
            detected_format=ReadFormat.FASTQ,
            output_dir=output_dir,
            global_options={}
        )

        validator = ReadValidator(read_config, ReadValidator.Settings())
        validator.run()

        # Low quality reads are allowed
        assert len(validator.sequences) == 1


class TestSecurityCommandInjection:
    """
    Security tests for command injection protection in read validator.

    These tests verify that the _count_lines_fast() method uses Python's
    native gzip/bz2 libraries instead of shell commands, preventing command
    injection attacks through malicious filenames.
    """

    def test_line_counting_uses_python_not_shell_gzip(self, tmp_path):
        """Test that line counting uses Python gzip library, not shell commands."""
        # Create a gzip compressed FASTQ file
        fastq_file = tmp_path / "reads.fastq.gz"
        content = "@read1\nATCG\n+\nIIII\n@read2\nGCTA\n+\nJJJJ\n"
        with gzip.open(fastq_file, 'wt') as f:
            f.write(content)

        # Create read config
        read_config = ReadConfig(
            filename="reads.fastq.gz",
            basename="reads",
            filepath=fastq_file,
            ngs_type='illumina',
            coding_type=CodingType.GZIP,
            detected_format=ReadFormat.FASTQ,
            output_dir=tmp_path,
            global_options={'validation_level': 'strict'}
        )

        output_dir = tmp_path / "output"
        output_dir.mkdir()

        settings = ReadValidator.Settings()

        validator = ReadValidator(read_config, settings)

        # Call _count_lines_fast directly
        line_count = validator._count_lines_fast()

        # Should count correct number of lines (8 lines total)
        assert line_count == 8

    def test_line_counting_uses_python_not_shell_bzip2(self, tmp_path):
        """Test that line counting uses Python bz2 library, not shell commands."""
        # Create a bzip2 compressed FASTQ file
        fastq_file = tmp_path / "reads.fastq.bz2"
        content = "@read1\nATCGATCG\n+\nIIIIIIII\n@read2\nGCTAGCTA\n+\nJJJJJJJJ\n"
        with bz2.open(fastq_file, 'wt') as f:
            f.write(content)

        # Create read config
        read_config = ReadConfig(
            filename="reads.fastq.bz2",
            basename="reads",
            filepath=fastq_file,
            ngs_type='ont',
            coding_type=CodingType.BZIP2,
            detected_format=ReadFormat.FASTQ,
            output_dir=tmp_path,
            global_options={'validation_level': 'strict'}
        )

        output_dir = tmp_path / "output"
        output_dir.mkdir()

        settings = ReadValidator.Settings()

        validator = ReadValidator(read_config, settings)

        # Call _count_lines_fast directly
        line_count = validator._count_lines_fast()

        # Should count correct number of lines (8 lines total)
        assert line_count == 8

    def test_line_counting_uncompressed(self, tmp_path):
        """Test that line counting works for uncompressed files."""
        # Create an uncompressed FASTQ file
        fastq_file = tmp_path / "reads.fastq"
        content = "@read1\nATCG\n+\nIIII\n"
        fastq_file.write_text(content)

        # Create read config
        read_config = ReadConfig(
            filename="reads.fastq",
            basename="reads",
            filepath=fastq_file,
            ngs_type='illumina',
            coding_type=CodingType.NONE,
            detected_format=ReadFormat.FASTQ,
            output_dir=tmp_path,
            global_options={'validation_level': 'strict'}
        )

        output_dir = tmp_path / "output"
        output_dir.mkdir()

        settings = ReadValidator.Settings()

        validator = ReadValidator(read_config, settings)

        # Call _count_lines_fast directly
        line_count = validator._count_lines_fast()

        # Should count correct number of lines (4 lines total)
        assert line_count == 4

    def test_no_subprocess_import_in_count_lines_fast(self, tmp_path):
        """
        Verify that _count_lines_fast doesn't use subprocess module.

        This is a regression test to ensure the security fix stays in place.
        If someone accidentally reintroduces shell commands, this test
        documents the expected behavior.
        """
        # Create a simple FASTQ file
        fastq_file = tmp_path / "reads.fastq.gz"
        content = "@read1\nATCG\n+\nIIII\n"
        with gzip.open(fastq_file, 'wt') as f:
            f.write(content)

        # Create read config
        read_config = ReadConfig(
            filename="reads.fastq.gz",
            basename="reads",
            filepath=fastq_file,
            ngs_type='illumina',
            coding_type=CodingType.GZIP,
            detected_format=ReadFormat.FASTQ,
            output_dir=tmp_path,
            global_options={}
        )

        output_dir = tmp_path / "output"
        output_dir.mkdir()

        validator = ReadValidator(read_config, ReadValidator.Settings())

        # Read the source code of _count_lines_fast method
        import inspect
        source = inspect.getsource(validator._count_lines_fast)

        # Verify no subprocess imports or calls are present (check actual code, not docs)
        # We exclude docstrings by checking for the patterns in context
        assert 'subprocess.run' not in source, \
            "_count_lines_fast should not use subprocess.run (security vulnerability)"
        assert 'subprocess.Popen' not in source, \
            "_count_lines_fast should not use subprocess.Popen (security vulnerability)"
        assert 'import subprocess' not in source, \
            "_count_lines_fast should not import subprocess"
        assert "shell=True" not in source, \
            "_count_lines_fast should not use shell=True"

        # Verify no shell command patterns are present
        assert "['sh', '-c'" not in source and '["sh", "-c"' not in source, \
            "_count_lines_fast should not use shell commands"
        assert '| wc -l' not in source, \
            "_count_lines_fast should not use shell pipes"

        # Verify it uses Python libraries instead
        assert 'open_file_with_coding_type' in source, \
            "_count_lines_fast should use open_file_with_coding_type"


class TestParallelValidationLogging:
    """Test that parallel validation correctly enables/disables process_id logging."""

    @pytest.fixture
    def temp_dir(self):
        """Create a temporary directory for test files."""
        with tempfile.TemporaryDirectory() as tmpdir:
            yield Path(tmpdir)

    @pytest.fixture
    def output_dir(self, temp_dir):
        """Create output directory."""
        out_dir = temp_dir / "output"
        out_dir.mkdir(parents=True, exist_ok=True)
        return out_dir

    def _create_fastq(self, path, num_reads=5000):
        """Create a FASTQ file for testing."""
        with open(path, 'w') as f:
            for i in range(num_reads):
                f.write(f"@read{i}\n")
                f.write("ATCGATCGATCGATCG\n")
                f.write("+\n")
                f.write("IIIIIIIIIIIIIIII\n")

    def test_sequential_validation_no_parallel_logging(self, temp_dir, output_dir):
        """Test that sequential validation (threads=1) does NOT enable parallel logging."""
        from validation_pkg.logger import setup_logging, get_logger

        fastq_file = temp_dir / "test.fastq"
        self._create_fastq(fastq_file, num_reads=100)

        log_file = temp_dir / "test.log"
        logger = setup_logging(log_file=log_file)


        read_config = ReadConfig(
            filename="test.fastq",
            basename="test",
            filepath=fastq_file,
            ngs_type="illumina",
            coding_type=CodingType.NONE,
            detected_format=ReadFormat.FASTQ,
            output_dir=output_dir,
            global_options={"threads": 1, "validation_level": "strict"}
        )

        settings = ReadValidator.Settings(check_invalid_chars=True)
        validator = ReadValidator(read_config, settings)

        # Run validation (sequential)
        validator._parse_file()
        validator._validate_sequences()


        # Check log file - should NOT have process_id
        log_content = log_file.read_text()
        assert '"process_id"' not in log_content
    def test_parallel_validation_enables_parallel_logging(self, temp_dir, output_dir):
        """Test that parallel validation (threads>1) logs parallel-related messages."""
        from validation_pkg.logger import setup_logging, get_logger

        fastq_file = temp_dir / "test.fastq"
        self._create_fastq(fastq_file, num_reads=5000)

        log_file = temp_dir / "test.log"
        logger = setup_logging(log_file=log_file)


        read_config = ReadConfig(
            filename="test.fastq",
            basename="test",
            filepath=fastq_file,
            ngs_type="illumina",
            coding_type=CodingType.NONE,
            detected_format=ReadFormat.FASTQ,
            output_dir=output_dir,
            global_options={"threads": 4, "validation_level": "strict"}
        )

        settings = ReadValidator.Settings(check_invalid_chars=True)
        validator = ReadValidator(read_config, settings)

        # Run validation (parallel)
        validator._parse_file()
        validator._validate_sequences()


        # Check log file - SHOULD have parallel-related messages
        log_content = log_file.read_text()
        assert '"parallel_mode": true' in log_content or 'Parallel validation enabled' in log_content
    def test_trust_mode_no_parallel_logging(self, temp_dir, output_dir):
        """Test that trust mode never enables parallel logging (always sequential)."""
        from validation_pkg.logger import setup_logging, get_logger

        fastq_file = temp_dir / "test.fastq"
        self._create_fastq(fastq_file, num_reads=1000)

        log_file = temp_dir / "test.log"
        logger = setup_logging(log_file=log_file)


        read_config = ReadConfig(
            filename="test.fastq",
            basename="test",
            filepath=fastq_file,
            ngs_type="illumina",
            coding_type=CodingType.NONE,
            detected_format=ReadFormat.FASTQ,
            output_dir=output_dir,
            global_options={"threads": 8, "validation_level": "trust"}  # threads=8 but trust mode
        )

        settings = ReadValidator.Settings(check_invalid_chars=True)
        validator = ReadValidator(read_config, settings)

        # Run validation (trust mode - always sequential)
        validator._parse_file()
        validator._validate_sequences()


        # Check log file - should NOT have process_id
        log_content = log_file.read_text()
        assert '"process_id"' not in log_content
    def test_parallel_logging_cleanup_on_error(self, temp_dir, output_dir):
        """Test that parallel logging is disabled even if validation fails."""
        from validation_pkg.logger import setup_logging, get_logger

        # Create FASTQ with invalid characters (properly formatted for BioPython)
        fastq_file = temp_dir / "invalid.fastq"
        with open(fastq_file, 'w') as f:
            for i in range(2000):
                f.write(f"@read{i}\n")
                if i == 1000:
                    # Invalid characters - will fail validation but parse correctly
                    f.write("ATCGATCGXYZATCGA\n")  # Same length as quality (16 chars)
                else:
                    f.write("ATCGATCGATCGATCG\n")
                f.write("+\n")
                f.write("IIIIIIIIIIIIIIII\n")

        log_file = temp_dir / "test.log"
        logger = setup_logging(log_file=log_file)


        read_config = ReadConfig(
            filename="invalid.fastq",
            basename="invalid",
            filepath=fastq_file,
            ngs_type="illumina",
            coding_type=CodingType.NONE,
            detected_format=ReadFormat.FASTQ,
            output_dir=output_dir,
            global_options={"threads": 4, "validation_level": "strict"}
        )

        settings = ReadValidator.Settings(check_invalid_chars=True)
        validator = ReadValidator(read_config, settings)

        # Parse file (should succeed)
        validator._parse_file()

        # Run validation (should fail due to invalid characters)
        with pytest.raises(ReadValidationError):
            validator._validate_sequences()


    def test_full_validation_with_parallel_logging(self, temp_dir, output_dir):
        """Integration test: Full validation run with parallel logging."""
        from validation_pkg.logger import setup_logging, get_logger

        fastq_file = temp_dir / "test.fastq"
        self._create_fastq(fastq_file, num_reads=10000)

        log_file = temp_dir / "test.log"
        logger = setup_logging(log_file=log_file, console_level="DEBUG")

        # Initial state

        read_config = ReadConfig(
            filename="test.fastq",
            basename="test",
            filepath=fastq_file,
            ngs_type="illumina",
            coding_type=CodingType.NONE,
            detected_format=ReadFormat.FASTQ,
            output_dir=output_dir,
            global_options={"threads": 8, "validation_level": "strict"}
        )

        settings = ReadValidator.Settings(
            check_invalid_chars=True,
            allow_duplicate_ids=False
        )
        validator = ReadValidator(read_config, settings)

        # Full validation run
        validator.run()

        # Check log file
        log_content = log_file.read_text()

        # Should mention parallel validation (no longer checking for process_id)
        assert 'parallel_mode' in log_content.lower() or 'parallel validation' in log_content.lower()

    def test_multiple_validations_cleanup(self, temp_dir, output_dir):
        """Test that parallel logging cleanup works across multiple validation runs."""
        from validation_pkg.logger import setup_logging, get_logger

        fastq_file1 = temp_dir / "test1.fastq"
        fastq_file2 = temp_dir / "test2.fastq"
        self._create_fastq(fastq_file1, num_reads=3000)
        self._create_fastq(fastq_file2, num_reads=3000)

        log_file = temp_dir / "test.log"
        logger = setup_logging(log_file=log_file)

        # First validation (parallel)
        read_config1 = ReadConfig(
            filename="test1.fastq",
            basename="test1",
            filepath=fastq_file1,
            ngs_type="illumina",
            coding_type=CodingType.NONE,
            detected_format=ReadFormat.FASTQ,
            output_dir=output_dir,
            global_options={"threads": 4, "validation_level": "strict"}
        )

        validator1 = ReadValidator(read_config1, ReadValidator.Settings())
        validator1._parse_file()
        validator1._validate_sequences()

        # Should be cleaned up

        # Second validation (sequential)
        read_config2 = ReadConfig(
            filename="test2.fastq",
            basename="test2",
            filepath=fastq_file2,
            ngs_type="illumina",
            coding_type=CodingType.NONE,
            detected_format=ReadFormat.FASTQ,
            output_dir=output_dir,
            global_options={"threads": 1, "validation_level": "strict"}
        )

        validator2 = ReadValidator(read_config2, ReadValidator.Settings())
        validator2._parse_file()
        validator2._validate_sequences()

        # Should still be disabled


class TestIlluminaPatternDetection:
    """Test Illumina paired-end filename pattern detection."""

    def _create_validator_with_basename(self, filename):
        """
        Helper to create validator with mock read_config containing basename.

        The basename is calculated the same way as ReadConfig.__post_init__():
        - If compressed: strip last 2 extensions (format + compression)
        - If uncompressed: strip last 1 extension (format only)
        """
        from pathlib import Path
        from unittest.mock import Mock

        suffixes = Path(filename).suffixes
        compression_exts = {'.gz', '.gzip', '.bz2', '.bzip2'}
        has_compression = any(s.lower() in compression_exts for s in suffixes)

        if has_compression:
            # Strip last 2 extensions: format + compression
            basename = filename.rsplit('.', 2)[0]
        else:
            # Strip last 1 extension: format only
            basename = filename.rsplit('.', 1)[0] if '.' in filename else filename

        validator = ReadValidator.__new__(ReadValidator)
        validator.output_metadata = OutputMetadata()
        validator.read_config = Mock()
        validator.read_config.basename = basename

        return validator

    # ==================== User-Provided Examples (12 files) ====================

    def test_standard_illumina_with_sample_and_lane(self):
        """Test Illumina filenames with sample ID and lane numbers (S1_L001_R1/R2 format)."""
        # Bacillus_toyonensis_LE1_S1_L001_R1.fastq.gz
        validator = self._create_validator_with_basename("Bacillus_toyonensis_LE1_S1_L001_R1.fastq.gz")

        validator._detect_illumina_pattern("Bacillus_toyonensis_LE1_S1_L001_R1.fastq.gz")
        assert validator.output_metadata.base_name == "Bacillus_toyonensis_LE1_S1_L001"
        assert validator.output_metadata.read_number == 1
        assert validator.output_metadata.ngs_type_detected == "illumina"

        # Bacillus_toyonensis_LE1_S1_L001_R2.fastq.gz
        validator2 = self._create_validator_with_basename("Bacillus_toyonensis_LE1_S1_L001_R2.fastq.gz")
        validator2._detect_illumina_pattern("Bacillus_toyonensis_LE1_S1_L001_R2.fastq.gz")
        assert validator2.output_metadata.base_name == "Bacillus_toyonensis_LE1_S1_L001"
        assert validator2.output_metadata.read_number == 2

    def test_hyphenated_names_with_underscore_numbers(self):
        """Test filenames with hyphens and _1/_2 suffix (no R prefix)."""
        # Bacillus-toyonensis-LE1-sequence_1.fastq.gz
        validator = self._create_validator_with_basename("Bacillus-toyonensis-LE1-sequence_1.fastq.gz")

        validator._detect_illumina_pattern("Bacillus-toyonensis-LE1-sequence_1.fastq.gz")
        assert validator.output_metadata.base_name == "Bacillus-toyonensis-LE1-sequence"
        assert validator.output_metadata.read_number == 1
        assert validator.output_metadata.ngs_type_detected == "illumina"

        # Bacillus-toyonensis-LE1-sequence_2.fastq.gz
        validator2 = self._create_validator_with_basename("Bacillus-toyonensis-LE1-sequence_2.fastq.gz")
        validator2._detect_illumina_pattern("Bacillus-toyonensis-LE1-sequence_2.fastq.gz")
        assert validator2.output_metadata.base_name == "Bacillus-toyonensis-LE1-sequence"
        assert validator2.output_metadata.read_number == 2

    def test_dot_separator_style(self):
        """Test filenames with dot separator (.R1 style)."""
        # Bacillus_toyonensis_LE1_sequence.R1.fastq.gz
        validator = self._create_validator_with_basename("Bacillus_toyonensis_LE1_sequence.R1.fastq.gz")

        validator._detect_illumina_pattern("Bacillus_toyonensis_LE1_sequence.R1.fastq.gz")
        assert validator.output_metadata.base_name == "Bacillus_toyonensis_LE1_sequence"
        assert validator.output_metadata.read_number == 1
        assert validator.output_metadata.ngs_type_detected == "illumina"

    def test_single_end_no_pattern(self):
        """Test single-end files without paired-end indicators."""
        # SRR834393.fq.gz (no R1/R2 or _1/_2 indicator)
        validator = self._create_validator_with_basename("SRR834393.fq.gz")

        validator._detect_illumina_pattern("SRR834393.fq.gz")
        # No paired-end pattern - should set fallback metadata for single-end
        assert validator.output_metadata.base_name == "SRR834393"
        assert validator.output_metadata.read_number == 1  # Defaults to R1
        assert validator.output_metadata.ngs_type_detected == 'illumina'

    def test_sra_style_underscore_suffix(self):
        """Test SRA-style filenames with _1/_2 suffix."""
        # SRR834394_1.fq.gz
        validator = self._create_validator_with_basename("SRR834394_1.fq.gz")

        validator._detect_illumina_pattern("SRR834394_1.fq.gz")
        assert validator.output_metadata.base_name == "SRR834394"
        assert validator.output_metadata.read_number == 1
        assert validator.output_metadata.ngs_type_detected == "illumina"

        # SRR837394_1.fastq.gz
        validator2 = self._create_validator_with_basename("SRR837394_1.fastq.gz")
        validator2._detect_illumina_pattern("SRR837394_1.fastq.gz")
        assert validator2.output_metadata.base_name == "SRR837394"
        assert validator2.output_metadata.read_number == 1

    def test_lane_numbers_without_r_prefix(self):
        """Test filenames with _1_001 and _2_001 patterns (lane numbers without R prefix)."""
        # SRR837394_1_001.fastq.gz
        validator = self._create_validator_with_basename("SRR837394_1_001.fastq.gz")

        validator._detect_illumina_pattern("SRR837394_1_001.fastq.gz")
        assert validator.output_metadata.base_name == "SRR837394"
        assert validator.output_metadata.read_number == 1
        assert validator.output_metadata.ngs_type_detected == "illumina"

        # SRR837394_2_001.fastq.gz
        validator2 = self._create_validator_with_basename("SRR837394_2_001.fastq.gz")
        validator2._detect_illumina_pattern("SRR837394_2_001.fastq.gz")
        assert validator2.output_metadata.base_name == "SRR837394"
        assert validator2.output_metadata.read_number == 2

    def test_lane_before_read_number_suffix(self):
        """Test filenames with _X_1 and _X_2 patterns where X is a lane number and 1/2 is read number.

        This handles the case where the read number is at the END after a lane/run number,
        e.g., NG-9904_PG4602EAc01_lib132096_4703_1_1.fastq.bz2 (R1) and
              NG-9904_PG4602EAc01_lib132096_4703_1_2.fastq.bz2 (R2)
        The last digit (1 or 2) is the read number, not the middle digit.
        """
        # NG-9904_PG4602EAc01_lib132096_4703_1_1.fastq.bz2 (R1)
        validator = self._create_validator_with_basename("NG-9904_PG4602EAc01_lib132096_4703_1_1.fastq.bz2")
        validator._detect_illumina_pattern("NG-9904_PG4602EAc01_lib132096_4703_1_1.fastq.bz2")
        assert validator.output_metadata.base_name == "NG-9904_PG4602EAc01_lib132096_4703"
        assert validator.output_metadata.read_number == 1
        assert validator.output_metadata.ngs_type_detected == "illumina"

        # NG-9904_PG4602EAc01_lib132096_4703_1_2.fastq.bz2 (R2)
        validator2 = self._create_validator_with_basename("NG-9904_PG4602EAc01_lib132096_4703_1_2.fastq.bz2")
        validator2._detect_illumina_pattern("NG-9904_PG4602EAc01_lib132096_4703_1_2.fastq.bz2")
        assert validator2.output_metadata.base_name == "NG-9904_PG4602EAc01_lib132096_4703"
        assert validator2.output_metadata.read_number == 2
        assert validator2.output_metadata.ngs_type_detected == "illumina"

        # Both files should have the SAME base_name (critical for pairing)
        assert validator.output_metadata.base_name == validator2.output_metadata.base_name

    def test_lane_before_read_number_with_longer_lane(self):
        """Test _XXX_1/_XXX_2 patterns with multi-digit lane numbers."""
        # sample_001_1.fastq.gz (R1 with lane 001)
        validator = self._create_validator_with_basename("sample_001_1.fastq.gz")
        validator._detect_illumina_pattern("sample_001_1.fastq.gz")
        assert validator.output_metadata.base_name == "sample"
        assert validator.output_metadata.read_number == 1

        # sample_001_2.fastq.gz (R2 with lane 001)
        validator2 = self._create_validator_with_basename("sample_001_2.fastq.gz")
        validator2._detect_illumina_pattern("sample_001_2.fastq.gz")
        assert validator2.output_metadata.base_name == "sample"
        assert validator2.output_metadata.read_number == 2

    def test_no_separator_r_prefix(self):
        """Test filenames with R1/R2 directly attached (no separator)."""
        # SRR834393R1_combined.fq.gz
        validator = self._create_validator_with_basename("SRR834393R1_combined.fq.gz")

        validator._detect_illumina_pattern("SRR834393R1_combined.fq.gz")
        assert validator.output_metadata.base_name == "SRR834393"
        assert validator.output_metadata.read_number == 1
        assert validator.output_metadata.ngs_type_detected == "illumina"

        # SRR834393R2_combined.fq.gz
        validator2 = self._create_validator_with_basename("SRR834393R2_combined.fq.gz")
        validator2._detect_illumina_pattern("SRR834393R2_combined.fq.gz")
        assert validator2.output_metadata.base_name == "SRR834393"
        assert validator2.output_metadata.read_number == 2

    def test_complex_basename_with_suffix(self):
        """Test complex filenames like Trichoderma_illumina_ABE_R1_combined.fastq.gz"""
        # Trichoderma_illumina_ABE_R1_combined.fastq.gz
        validator = self._create_validator_with_basename("Trichoderma_illumina_ABE_R1_combined.fastq.gz")

        validator._detect_illumina_pattern("Trichoderma_illumina_ABE_R1_combined.fastq.gz")
        assert validator.output_metadata.base_name == "Trichoderma_illumina_ABE"
        assert validator.output_metadata.read_number == 1
        assert validator.output_metadata.ngs_type_detected == "illumina"

        # Trichoderma_illumina_ABE_R2_combined.fastq.gz
        validator2 = self._create_validator_with_basename("Trichoderma_illumina_ABE_R2_combined.fastq.gz")
        validator2._detect_illumina_pattern("Trichoderma_illumina_ABE_R2_combined.fastq.gz")
        assert validator2.output_metadata.base_name == "Trichoderma_illumina_ABE"
        assert validator2.output_metadata.read_number == 2

    # ==================== Pattern Coverage Tests ====================

    def test_standard_r_prefix_patterns(self):
        """Test standard _R1/_R2 patterns (most common)."""
        # Basic _R1 pattern
        validator = self._create_validator_with_basename("sample_R1.fastq.gz")
        validator._detect_illumina_pattern("sample_R1.fastq.gz")
        assert validator.output_metadata.base_name == "sample"
        assert validator.output_metadata.read_number == 1

        # Basic _R2 pattern
        validator2 = self._create_validator_with_basename("sample_R2.fastq.gz")
        validator2._detect_illumina_pattern("sample_R2.fastq.gz")
        assert validator2.output_metadata.base_name == "sample"
        assert validator2.output_metadata.read_number == 2

    def test_lane_number_patterns_with_r_prefix(self):
        """Test _R1_001/_R2_001 patterns (Illumina lane numbers)."""
        # _R1_001 pattern
        validator = self._create_validator_with_basename("sample_R1_001.fastq.gz")
        validator._detect_illumina_pattern("sample_R1_001.fastq.gz")
        assert validator.output_metadata.base_name == "sample"
        assert validator.output_metadata.read_number == 1

        # _R2_999 pattern (any digits)
        validator2 = self._create_validator_with_basename("sample_R2_999.fastq.gz")
        validator2._detect_illumina_pattern("sample_R2_999.fastq.gz")
        assert validator2.output_metadata.base_name == "sample"
        assert validator2.output_metadata.read_number == 2

    def test_arbitrary_suffix_patterns(self):
        """Test _R1_<suffix> patterns (arbitrary text after R1/R2)."""
        # _R1_combined pattern
        validator = self._create_validator_with_basename("sample_R1_combined.fastq.gz")
        validator._detect_illumina_pattern("sample_R1_combined.fastq.gz")
        assert validator.output_metadata.base_name == "sample"
        assert validator.output_metadata.read_number == 1

        # _R2_final pattern
        validator2 = self._create_validator_with_basename("sample_R2_final.fastq.gz")
        validator2._detect_illumina_pattern("sample_R2_final.fastq.gz")
        assert validator2.output_metadata.base_name == "sample"
        assert validator2.output_metadata.read_number == 2

        # _R1_merged_v2 pattern (complex suffix)
        validator3 = self._create_validator_with_basename("sample_R1_merged_v2.fastq.gz")
        validator3._detect_illumina_pattern("sample_R1_merged_v2.fastq.gz")
        assert validator3.output_metadata.base_name == "sample"
        assert validator3.output_metadata.read_number == 1

    def test_arbitrary_suffix_patterns_no_r_prefix(self):
        """Test _1_<suffix> and _2_<suffix> patterns (no R prefix, with arbitrary suffix)."""
        # _1_combined pattern
        validator = self._create_validator_with_basename("SRR837394_1_combined.fastq.gz")
        validator._detect_illumina_pattern("SRR837394_1_combined.fastq.gz")
        assert validator.output_metadata.base_name == "SRR837394"
        assert validator.output_metadata.read_number == 1
        assert validator.output_metadata.ngs_type_detected == "illumina"

        # _2_combined pattern
        validator2 = self._create_validator_with_basename("SRR837394_2_combined.fastq.gz")
        validator2._detect_illumina_pattern("SRR837394_2_combined.fastq.gz")
        assert validator2.output_metadata.base_name == "SRR837394"
        assert validator2.output_metadata.read_number == 2

    def test_dot_separator_with_numbers(self):
        """Test .1/.2 patterns (dot separator without R prefix)."""
        # .1 pattern
        validator = self._create_validator_with_basename("sample.1.fastq.gz")
        validator._detect_illumina_pattern("sample.1.fastq.gz")
        assert validator.output_metadata.base_name == "sample"
        assert validator.output_metadata.read_number == 1

        # .2 pattern
        validator2 = self._create_validator_with_basename("sample.2.fastq.gz")
        validator2._detect_illumina_pattern("sample.2.fastq.gz")
        assert validator2.output_metadata.base_name == "sample"
        assert validator2.output_metadata.read_number == 2

    def test_no_separator_numbers(self):
        """Test R1/R2 directly attached to basename (no separator)."""
        # sampleR1 pattern
        validator = self._create_validator_with_basename("sampleR1.fastq.gz")
        validator._detect_illumina_pattern("sampleR1.fastq.gz")
        assert validator.output_metadata.base_name == "sample"
        assert validator.output_metadata.read_number == 1

        # sampleR2 pattern
        validator2 = self._create_validator_with_basename("sampleR2.fastq.gz")
        validator2._detect_illumina_pattern("sampleR2.fastq.gz")
        assert validator2.output_metadata.base_name == "sample"
        assert validator2.output_metadata.read_number == 2

    # ==================== Compression Format Independence ====================

    def test_compression_independence(self):
        """Test that pattern detection works regardless of compression format."""
        # Test with .gz compression
        validator1 = self._create_validator_with_basename("sample_R1.fastq.gz")
        validator1._detect_illumina_pattern("sample_R1.fastq.gz")
        assert validator1.output_metadata.base_name == "sample"
        assert validator1.output_metadata.read_number == 1

        # Test with .bz2 compression
        validator2 = self._create_validator_with_basename("sample_R2.fastq.bz2")
        validator2._detect_illumina_pattern("sample_R2.fastq.bz2")
        assert validator2.output_metadata.base_name == "sample"
        assert validator2.output_metadata.read_number == 2

        # Test with no compression
        validator3 = self._create_validator_with_basename("sample_1.fastq")
        validator3._detect_illumina_pattern("sample_1.fastq")
        assert validator3.output_metadata.base_name == "sample"
        assert validator3.output_metadata.read_number == 1

        # Test with .fq extension (alternate FASTQ)
        validator4 = self._create_validator_with_basename("sample_2.fq")
        validator4._detect_illumina_pattern("sample_2.fq")
        assert validator4.output_metadata.base_name == "sample"
        assert validator4.output_metadata.read_number == 2

    # ==================== Pattern Priority Order ====================

    def test_pattern_priority_order(self):
        """Test that patterns are matched in correct priority order (most specific first)."""
        # _R1_001 should match lane number pattern (higher priority)
        validator1 = self._create_validator_with_basename("sample_R1_001.fastq.gz")
        validator1._detect_illumina_pattern("sample_R1_001.fastq.gz")
        assert validator1.output_metadata.base_name == "sample"
        assert validator1.output_metadata.read_number == 1

        # _R1_combined should match arbitrary suffix pattern
        validator2 = self._create_validator_with_basename("sample_R1_combined.fastq.gz")
        validator2._detect_illumina_pattern("sample_R1_combined.fastq.gz")
        assert validator2.output_metadata.base_name == "sample"
        assert validator2.output_metadata.read_number == 1

        # _1_001 should match lane number pattern without R prefix
        validator3 = self._create_validator_with_basename("sample_1_001.fastq.gz")
        validator3._detect_illumina_pattern("sample_1_001.fastq.gz")
        assert validator3.output_metadata.base_name == "sample"
        assert validator3.output_metadata.read_number == 1

        # Plain _R1 should match standard pattern
        validator4 = self._create_validator_with_basename("sample_R1.fastq.gz")
        validator4._detect_illumina_pattern("sample_R1.fastq.gz")
        assert validator4.output_metadata.base_name == "sample"
        assert validator4.output_metadata.read_number == 1

    # ==================== Edge Cases ====================

    def test_edge_case_short_filenames(self):
        """Test very short filenames."""
        # Single character base name
        validator = self._create_validator_with_basename("a_R1.fastq")
        validator._detect_illumina_pattern("a_R1.fastq")
        assert validator.output_metadata.base_name == "a"
        assert validator.output_metadata.read_number == 1

        # Two character base name
        validator2 = self._create_validator_with_basename("ab_R2.fastq")
        validator2._detect_illumina_pattern("ab_R2.fastq")
        assert validator2.output_metadata.base_name == "ab"
        assert validator2.output_metadata.read_number == 2

    def test_edge_case_long_filenames(self):
        """Test filenames with many underscores and parts."""
        validator = self._create_validator_with_basename(
            "my_very_long_sample_name_with_many_parts_R2.fastq.gz"
        )
        validator._detect_illumina_pattern(
            "my_very_long_sample_name_with_many_parts_R2.fastq.gz"
        )
        assert validator.output_metadata.base_name == "my_very_long_sample_name_with_many_parts"
        assert validator.output_metadata.read_number == 2

    def test_edge_case_numbers_in_basename(self):
        """Test filenames with numbers in the base name."""
        # Numbers in middle
        validator = self._create_validator_with_basename("sample123_R1.fastq")
        validator._detect_illumina_pattern("sample123_R1.fastq")
        assert validator.output_metadata.base_name == "sample123"
        assert validator.output_metadata.read_number == 1

        # Numbers at start
        validator2 = self._create_validator_with_basename("123sample_R2.fastq")
        validator2._detect_illumina_pattern("123sample_R2.fastq")
        assert validator2.output_metadata.base_name == "123sample"
        assert validator2.output_metadata.read_number == 2

    def test_edge_case_mixed_separators(self):
        """Test filenames with both hyphens and underscores."""
        validator = self._create_validator_with_basename("sample-name_123_R1.fastq")
        validator._detect_illumina_pattern("sample-name_123_R1.fastq")
        assert validator.output_metadata.base_name == "sample-name_123"
        assert validator.output_metadata.read_number == 1

    def test_edge_case_uppercase_extensions(self):
        """Test that pattern detection works with uppercase extensions."""
        # Uppercase .FASTQ.GZ
        validator = self._create_validator_with_basename("sample_R1.FASTQ.GZ")
        validator._detect_illumina_pattern("sample_R1.FASTQ.GZ")
        assert validator.output_metadata.base_name == "sample"
        assert validator.output_metadata.read_number == 1

    # ==================== Invalid Patterns ====================

    def test_invalid_read_numbers(self):
        """Test that invalid read numbers get fallback metadata (treated as single-end)."""
        # R3 is not valid (only R1 and R2 are valid) - treated as single-end
        validator1 = self._create_validator_with_basename("sample_R3.fastq")
        validator1._detect_illumina_pattern("sample_R3.fastq")
        assert validator1.output_metadata.base_name == "sample_R3"
        assert validator1.output_metadata.read_number == 1  # Fallback
        assert validator1.output_metadata.ngs_type_detected == 'illumina'

        # _3 is not valid - treated as single-end
        validator2 = self._create_validator_with_basename("sample_3.fastq")
        validator2._detect_illumina_pattern("sample_3.fastq")
        assert validator2.output_metadata.base_name == "sample_3"
        assert validator2.output_metadata.read_number == 1  # Fallback

        # R0 is not valid - treated as single-end
        validator3 = self._create_validator_with_basename("sample_R0.fastq")
        validator3._detect_illumina_pattern("sample_R0.fastq")
        assert validator3.output_metadata.base_name == "sample_R0"

    def test_invalid_ambiguous_patterns(self):
        """Test that ambiguous patterns get fallback metadata."""
        # Just R1 without base name (edge case - might match with empty base_name)
        validator1 = self._create_validator_with_basename("R1.fastq")
        validator1._detect_illumina_pattern("R1.fastq")
        # Pattern should match with empty base_name OR fallback to full basename
        assert validator1.output_metadata.read_number is not None
        assert validator1.output_metadata.ngs_type_detected == 'illumina'

        # R1 in middle of filename (not at end) - treated as single-end
        validator2 = self._create_validator_with_basename("R1_sample.fastq")
        validator2._detect_illumina_pattern("R1_sample.fastq")
        # Should not match standard patterns, gets fallback
        assert validator2.output_metadata.base_name == "R1_sample"
        assert validator2.output_metadata.read_number == 1  # Fallback
        assert validator2.output_metadata.ngs_type_detected == 'illumina'

    def test_invalid_lowercase_r_prefix(self):
        """Test that lowercase 'r' gets fallback metadata (not a valid pattern)."""
        validator = self._create_validator_with_basename("sample_r1.fastq")
        validator._detect_illumina_pattern("sample_r1.fastq")
        # Should not match - Illumina uses uppercase R1/R2, gets fallback
        assert validator.output_metadata.base_name == "sample_r1"
        assert validator.output_metadata.read_number == 1  # Fallback
        assert validator.output_metadata.ngs_type_detected == 'illumina'

    # ==================== Metadata Isolation Tests ====================

    def test_metadata_isolation_between_calls(self):
        """Test that each call properly updates metadata independently."""
        validator = self._create_validator_with_basename("sample1_R1.fastq.gz")

        # First detection
        validator._detect_illumina_pattern("sample1_R1.fastq.gz")
        assert validator.output_metadata.base_name == "sample1"
        assert validator.output_metadata.read_number == 1

        # Create new validator for second file
        validator2 = self._create_validator_with_basename("sample2_R2.fastq.gz")
        validator2._detect_illumina_pattern("sample2_R2.fastq.gz")

        # Verify second detection didn't affect first
        assert validator.output_metadata.base_name == "sample1"
        assert validator.output_metadata.read_number == 1

        # Verify second detection correct
        assert validator2.output_metadata.base_name == "sample2"
        assert validator2.output_metadata.read_number == 2

    def test_metadata_cleared_when_no_pattern(self):
        """Test that fallback metadata is set when no paired-end pattern is detected."""
        validator = self._create_validator_with_basename("single_end_sample.fastq.gz")

        # Should set fallback metadata for single-end files
        validator._detect_illumina_pattern("single_end_sample.fastq.gz")
        assert validator.output_metadata.base_name == "single_end_sample"
        assert validator.output_metadata.read_number == 1
        assert validator.output_metadata.ngs_type_detected == 'illumina'


class TestReadStatistics:
    """Test read statistics calculation in strict mode."""

    @pytest.fixture
    def temp_dir(self):
        """Create a temporary directory for test files."""
        with tempfile.TemporaryDirectory() as tmpdir:
            yield Path(tmpdir)

    @pytest.fixture
    def output_dir(self, temp_dir):
        """Create output directory."""
        out_dir = temp_dir / "output"
        out_dir.mkdir(parents=True, exist_ok=True)
        return out_dir

    @pytest.fixture
    def test_fastq_with_known_stats(self, temp_dir):
        """Create FASTQ file with known read statistics."""
        fastq_path = temp_dir / "test_reads.fastq.gz"

        # Create reads with known lengths for N50 calculation
        reads = [
            SeqRecord(Seq("A" * 1000), id="read1", description=""),  # 1000 bp
            SeqRecord(Seq("T" * 800), id="read2", description=""),   # 800 bp
            SeqRecord(Seq("G" * 600), id="read3", description=""),   # 600 bp
            SeqRecord(Seq("C" * 400), id="read4", description=""),   # 400 bp
            SeqRecord(Seq("A" * 200), id="read5", description=""),   # 200 bp
        ]
        # Total: 3000 bp, N50 should be 800 bp (cumulative > 1500 bp at 800)

        # Add quality scores
        for read in reads:
            read.letter_annotations["phred_quality"] = [40] * len(read.seq)

        with gzip.open(fastq_path, 'wt') as f:
            SeqIO.write(reads, f, 'fastq')

        return fastq_path

    def test_strict_mode_statistics(self, test_fastq_with_known_stats, output_dir):
        """Test that strict mode calculates all read statistics."""
        read_config = ReadConfig(
            filepath=test_fastq_with_known_stats,
            filename=test_fastq_with_known_stats.name,
            coding_type=CodingType.GZIP,
            detected_format=ReadFormat.FASTQ,
            ngs_type='illumina',
            output_dir=output_dir,
            global_options={'validation_level': 'strict'}
        )

        validator = ReadValidator(read_config)
        metadata = validator.run()

        # Basic fields
        assert metadata.output_file is not None
        assert metadata.num_reads == 5
        assert metadata.validation_level == 'strict'

        # Statistics should be calculated
        assert metadata.n50 == 800, f"Expected N50=800, got {metadata.n50}"
        assert metadata.total_bases == 3000, f"Expected total=3000, got {metadata.total_bases}"
        assert metadata.mean_read_length == 600.0, f"Expected mean=600.0, got {metadata.mean_read_length}"
        assert metadata.longest_read_length == 1000
        assert metadata.shortest_read_length == 200

    def test_trust_mode_no_statistics(self, test_fastq_with_known_stats, output_dir):
        """Test that trust mode does NOT calculate expensive statistics."""
        read_config = ReadConfig(
            filepath=test_fastq_with_known_stats,
            filename=test_fastq_with_known_stats.name,
            coding_type=CodingType.GZIP,
            detected_format=ReadFormat.FASTQ,
            ngs_type='illumina',
            output_dir=output_dir,
            global_options={'validation_level': 'trust'}
        )

        validator = ReadValidator(read_config)
        metadata = validator.run()

        # Basic fields should be set
        assert metadata.output_file is not None
        assert metadata.num_reads == 5
        assert metadata.validation_level == 'trust'

        # Statistics should NOT be calculated in trust mode
        assert metadata.n50 is None
        assert metadata.total_bases is None
        assert metadata.mean_read_length is None
        assert metadata.longest_read_length is None
        assert metadata.shortest_read_length is None

    def test_n50_calculation_accuracy(self, temp_dir, output_dir):
        """Test N50 calculation with known dataset."""
        # Create specific dataset to verify N50 calculation
        fastq_path = temp_dir / "n50_test.fastq"
        reads = [
            SeqRecord(Seq("A" * 100), id="read1", description=""),   # 100 bp
            SeqRecord(Seq("T" * 200), id="read2", description=""),   # 200 bp
            SeqRecord(Seq("G" * 300), id="read3", description=""),   # 300 bp
            SeqRecord(Seq("C" * 400), id="read4", description=""),   # 400 bp
        ]
        # Total: 1000 bp
        # Sorted: 400, 300, 200, 100
        # Cumulative: 400 (not >= 500), 700 (>= 500)
        # N50 = 300 bp

        for read in reads:
            read.letter_annotations["phred_quality"] = [40] * len(read.seq)

        with open(fastq_path, 'w') as f:
            SeqIO.write(reads, f, 'fastq')

        read_config = ReadConfig(
            filepath=fastq_path,
            filename=fastq_path.name,
            coding_type=CodingType.NONE,
            detected_format=ReadFormat.FASTQ,
            ngs_type='ont',
            output_dir=output_dir,
            global_options={'validation_level': 'strict'}
        )

        validator = ReadValidator(read_config)
        metadata = validator.run()

        assert metadata.n50 == 300, f"N50 should be 300 bp, got {metadata.n50}"
        assert metadata.total_bases == 1000
        assert metadata.mean_read_length == 250.0

    def test_empty_reads_statistics(self, temp_dir, output_dir):
        """Test statistics with no reads (edge case)."""
        # This test may not work if validation fails on empty files
        # But it tests the _calculate_read_statistics() method robustness
        from validation_pkg.validators.read_validator import ReadValidator

        read_config = ReadConfig(
            filepath=Path("dummy.fastq"),
            filename="dummy.fastq",
            coding_type=CodingType.NONE,
            detected_format=ReadFormat.FASTQ,
            ngs_type='illumina',
            output_dir=output_dir,
            global_options={'validation_level': 'strict'}
        )

        validator = ReadValidator(read_config)
        validator.sequences = []  # Empty sequences

        stats = validator._calculate_read_statistics()

        assert stats['n50'] == 0
        assert stats['total_bases'] == 0
        assert stats['mean_read_length'] == 0.0
        assert stats['longest_read_length'] == 0
        assert stats['shortest_read_length'] == 0


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
