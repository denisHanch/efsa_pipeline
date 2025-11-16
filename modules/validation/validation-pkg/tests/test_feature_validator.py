"""
Comprehensive tests for FeatureValidator.

Tests cover:
- File format detection and parsing (GFF, GTF, BED)
- Compression handling (gzip, bzip2, uncompressed)
- Validation rules (coordinate validation, zero-length features)
- Editing specifications (sorting by position)
- Statistics collection (counts by type, strand, length)
- Format conversion (BED to GFF coordinate conversion)
- Output generation with various compression formats
- Error handling and edge cases
"""

import pytest
import tempfile
import gzip
import bz2
from pathlib import Path

from validation_pkg.config_manager import FeatureConfig
from validation_pkg.validators.feature_validator import FeatureValidator, Feature
from validation_pkg.exceptions import (
    FeatureValidationError,
    GffFormatError,
    BedFormatError,
    CompressionError,
)
from validation_pkg.utils.formats import FeatureFormat, CodingType


class TestFeatureValidatorInitialization:
    """Test FeatureValidator initialization."""

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
    def simple_gff(self, temp_dir):
        """Create a simple GFF file."""
        gff_file = temp_dir / "features.gff"
        with open(gff_file, "w") as f:
            f.write("##gff-version 3\n")
            f.write("chr1\ttest\tgene\t100\t200\t.\t+\t.\tID=gene1\n")
        return gff_file

    def test_init_with_defaults(self, simple_gff, output_dir):
        """Test initialization with default settings."""
        feature_config = FeatureConfig(
            filename="features.gff",
            filepath=simple_gff,
            coding_type=CodingType.NONE,
            detected_format=FeatureFormat.GFF,
            output_dir=output_dir,
            global_options={}
        )

        validator = FeatureValidator(feature_config, FeatureValidator.Settings())

        assert validator.input_path == simple_gff
        assert validator.output_dir == output_dir
        assert validator.settings is not None
        assert validator.features == []

    def test_init_with_custom_settings(self, simple_gff, output_dir):
        """Test initialization with custom settings."""
        settings = FeatureValidator.Settings(
            sort_by_position=True,
            check_coordinates=True,
        )

        feature_config = FeatureConfig(
            filename="features.gff",
            filepath=simple_gff,
            coding_type=CodingType.NONE,
            detected_format=FeatureFormat.GFF,
            output_dir=output_dir,
            global_options={}
        )

        validator = FeatureValidator(feature_config, settings)

        assert validator.settings.sort_by_position is True
        assert validator.settings.check_coordinates is True


class TestFeatureValidatorParsing:
    """Test feature file parsing."""

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

    def test_parse_gff_file(self, temp_dir, output_dir):
        """Test parsing a valid GFF file."""
        gff_file = temp_dir / "test.gff"
        with open(gff_file, "w") as f:
            f.write("##gff-version 3\n")
            f.write("chr1\ttest\tgene\t100\t200\t.\t+\t.\tID=gene1;Name=TestGene\n")
            f.write("chr1\ttest\tCDS\t120\t180\t.\t+\t0\tID=cds1;Parent=gene1\n")

        feature_config = FeatureConfig(
            filename="test.gff",
            filepath=gff_file,
            coding_type=CodingType.NONE,
            detected_format=FeatureFormat.GFF,
            output_dir=output_dir,
            global_options={}
        )

        validator = FeatureValidator(feature_config, FeatureValidator.Settings())
        validator.run()

        assert len(validator.features) == 2
        assert validator.features[0].seqname == "chr1"
        assert validator.features[0].feature_type == "gene"
        assert validator.features[0].start == 100
        assert validator.features[0].end == 200
        assert validator.features[1].feature_type == "CDS"

    def test_parse_bed_file(self, temp_dir, output_dir):
        """Test parsing a valid BED file."""
        bed_file = temp_dir / "test.bed"
        with open(bed_file, "w") as f:
            f.write("chr1\t99\t200\tfeature1\t100\t+\n")  # BED is 0-based
            f.write("chr2\t500\t600\tfeature2\t200\t-\n")

        feature_config = FeatureConfig(
            filename="test.bed",
            filepath=bed_file,
            coding_type=CodingType.NONE,
            detected_format=FeatureFormat.BED,
            output_dir=output_dir,
            global_options={}
        )

        validator = FeatureValidator(feature_config, FeatureValidator.Settings())
        validator.run()

        assert len(validator.features) == 2
        assert validator.features[0].seqname == "chr1"
        assert validator.features[0].start == 100  # Converted to 1-based
        assert validator.features[0].end == 200
        assert validator.features[0].strand == "+"
        assert validator.features[1].seqname == "chr2"
        assert validator.features[1].strand == "-"

    def test_parse_empty_file(self, temp_dir, output_dir):
        """Test parsing an empty GFF file."""
        gff_file = temp_dir / "empty.gff"
        with open(gff_file, "w") as f:
            f.write("##gff-version 3\n")
            # No features

        feature_config = FeatureConfig(
            filename="empty.gff",
            filepath=gff_file,
            coding_type=CodingType.NONE,
            detected_format=FeatureFormat.GFF,
            output_dir=output_dir,
            global_options={}
        )

        validator = FeatureValidator(feature_config, FeatureValidator.Settings())
        validator.run()  # Should succeed with warning

        assert len(validator.features) == 0

    def test_parse_malformed_gff(self, temp_dir, output_dir):
        """Test parsing malformed GFF lines."""
        gff_file = temp_dir / "malformed.gff"
        with open(gff_file, "w") as f:
            f.write("##gff-version 3\n")
            f.write("chr1\ttest\tgene\t100\t200\t.\t+\t.\tID=gene1\n")  # Valid
            f.write("chr2\ttest\tgene\n")  # Invalid - too few columns
            f.write("chr3\ttest\tgene\tinvalid\t200\t.\t+\t.\tID=gene3\n")  # Invalid start

        feature_config = FeatureConfig(
            filename="malformed.gff",
            filepath=gff_file,
            coding_type=CodingType.NONE,
            detected_format=FeatureFormat.GFF,
            output_dir=output_dir,
            global_options={}
        )

        validator = FeatureValidator(feature_config, FeatureValidator.Settings())
        validator.run()

        # Should parse valid line, skip invalid ones
        assert len(validator.features) == 1
        assert validator.features[0].seqname == "chr1"


class TestFeatureValidatorCompression:
    """Test compressed file handling."""

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

    def test_gzip_compressed_gff(self, temp_dir, output_dir):
        """Test reading gzip-compressed GFF file."""
        gff_content = b"##gff-version 3\nchr1\ttest\tgene\t100\t200\t.\t+\t.\tID=gene1\n"
        gff_file = temp_dir / "test.gff.gz"

        with gzip.open(gff_file, "wb") as f:
            f.write(gff_content)

        feature_config = FeatureConfig(
            filename="test.gff.gz",
            filepath=gff_file,
            coding_type=CodingType.GZIP,
            detected_format=FeatureFormat.GFF,
            output_dir=output_dir,
            global_options={}
        )

        validator = FeatureValidator(feature_config, FeatureValidator.Settings())
        validator.run()

        assert len(validator.features) == 1
        assert validator.features[0].seqname == "chr1"

    def test_bzip2_compressed_bed(self, temp_dir, output_dir):
        """Test reading bzip2-compressed BED file."""
        bed_content = b"chr1\t99\t200\tfeature1\t100\t+\n"
        bed_file = temp_dir / "test.bed.bz2"

        with bz2.open(bed_file, "wb") as f:
            f.write(bed_content)

        feature_config = FeatureConfig(
            filename="test.bed.bz2",
            filepath=bed_file,
            coding_type=CodingType.BZIP2,
            detected_format=FeatureFormat.BED,
            output_dir=output_dir,
            global_options={}
        )

        validator = FeatureValidator(feature_config, FeatureValidator.Settings())
        validator.run()

        assert len(validator.features) == 1
        assert validator.features[0].seqname == "chr1"


class TestFeatureValidatorValidation:
    """Test feature validation rules."""

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

    def test_invalid_coordinates_start_greater_than_end(self, temp_dir, output_dir):
        """Test validation fails when start > end."""
        gff_file = temp_dir / "invalid.gff"
        with open(gff_file, "w") as f:
            f.write("##gff-version 3\n")
            f.write("chr1\ttest\tgene\t200\t100\t.\t+\t.\tID=gene1\n")  # start > end

        feature_config = FeatureConfig(
            filename="invalid.gff",
            filepath=gff_file,
            coding_type=CodingType.NONE,
            detected_format=FeatureFormat.GFF,
            output_dir=output_dir,
            global_options={}
        )

        settings = FeatureValidator.Settings(check_coordinates=True)
        validator = FeatureValidator(feature_config, settings)

        with pytest.raises(FeatureValidationError, match="start > end"):
            validator.run()

    def test_zero_length_feature_not_allowed(self, temp_dir, output_dir):
        """Test validation fails for zero-length features (hardcoded check)."""
        gff_file = temp_dir / "zero_length.gff"
        with open(gff_file, "w") as f:
            f.write("##gff-version 3\n")
            f.write("chr1\ttest\tgene\t100\t100\t.\t+\t.\tID=gene1\n")  # start == end

        feature_config = FeatureConfig(
            filename="zero_length.gff",
            filepath=gff_file,
            coding_type=CodingType.NONE,
            detected_format=FeatureFormat.GFF,
            output_dir=output_dir,
            global_options={}
        )

        settings = FeatureValidator.Settings(check_coordinates=True)
        validator = FeatureValidator(feature_config, settings)

        with pytest.raises(FeatureValidationError, match="Zero-length feature"):
            validator.run()

    def test_negative_coordinates(self, temp_dir, output_dir):
        """Test validation fails for negative coordinates (hardcoded check)."""
        gff_file = temp_dir / "negative.gff"
        with open(gff_file, "w") as f:
            f.write("##gff-version 3\n")
            f.write("chr1\ttest\tgene\t-10\t100\t.\t+\t.\tID=gene1\n")

        feature_config = FeatureConfig(
            filename="negative.gff",
            filepath=gff_file,
            coding_type=CodingType.NONE,
            detected_format=FeatureFormat.GFF,
            output_dir=output_dir,
            global_options={}
        )

        settings = FeatureValidator.Settings(check_coordinates=True)
        validator = FeatureValidator(feature_config, settings)

        with pytest.raises(FeatureValidationError, match="Negative coordinates"):
            validator.run()


class TestFeatureValidatorEditing:
    """Test editing operations."""

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

    def test_sort_by_position(self, temp_dir, output_dir):
        """Test sorting features by position."""
        gff_file = temp_dir / "unsorted.gff"
        with open(gff_file, "w") as f:
            f.write("##gff-version 3\n")
            f.write("chr2\ttest\tgene\t100\t200\t.\t+\t.\tID=gene2\n")
            f.write("chr1\ttest\tgene\t500\t600\t.\t+\t.\tID=gene3\n")
            f.write("chr1\ttest\tgene\t100\t200\t.\t+\t.\tID=gene1\n")

        feature_config = FeatureConfig(
            filename="unsorted.gff",
            filepath=gff_file,
            coding_type=CodingType.NONE,
            detected_format=FeatureFormat.GFF,
            output_dir=output_dir,
            global_options={}
        )

        settings = FeatureValidator.Settings(sort_by_position=True)
        validator = FeatureValidator(feature_config, settings)
        validator.run()

        # Check order: chr1:100, chr1:500, chr2:100
        assert validator.features[0].seqname == "chr1"
        assert validator.features[0].start == 100
        assert validator.features[1].seqname == "chr1"
        assert validator.features[1].start == 500
        assert validator.features[2].seqname == "chr2"

    def test_no_sorting(self, temp_dir, output_dir):
        """Test features maintain original order when sorting disabled."""
        gff_file = temp_dir / "unsorted.gff"
        with open(gff_file, "w") as f:
            f.write("##gff-version 3\n")
            f.write("chr2\ttest\tgene\t100\t200\t.\t+\t.\tID=gene2\n")
            f.write("chr1\ttest\tgene\t100\t200\t.\t+\t.\tID=gene1\n")

        feature_config = FeatureConfig(
            filename="unsorted.gff",
            filepath=gff_file,
            coding_type=CodingType.NONE,
            detected_format=FeatureFormat.GFF,
            output_dir=output_dir,
            global_options={}
        )

        settings = FeatureValidator.Settings(sort_by_position=False)
        validator = FeatureValidator(feature_config, settings)
        validator.run()

        # Check order maintained
        assert validator.features[0].seqname == "chr2"
        assert validator.features[1].seqname == "chr1"


class TestFeatureValidatorFormatConversion:
    """Test format conversion between BED and GFF."""

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

    def test_bed_to_gff_conversion(self, temp_dir, output_dir):
        """Test BED to GFF coordinate conversion."""
        bed_file = temp_dir / "test.bed"
        with open(bed_file, "w") as f:
            f.write("chr1\t0\t100\tfeature1\t100\t+\n")  # BED: 0-based [0, 100)
            f.write("chr1\t150\t200\tfeature2\t100\t-\n")  # BED: [150, 200)

        feature_config = FeatureConfig(
            filename="test.bed",
            filepath=bed_file,
            coding_type=CodingType.NONE,
            detected_format=FeatureFormat.BED,
            output_dir=output_dir,
            global_options={}
        )

        settings = FeatureValidator.Settings()
        validator = FeatureValidator(feature_config, settings)
        validator.run()

        # Check conversion to GFF (1-based)
        assert validator.features[0].start == 1  # 0 -> 1
        assert validator.features[0].end == 100
        assert validator.features[1].start == 151  # 150 -> 151
        assert validator.features[1].end == 200

        # Verify output is GFF format
        output_file = output_dir / "test.gff"
        assert output_file.exists()

        with open(output_file, "r") as f:
            content = f.read()
            assert "##gff-version 3" in content
            assert "\t1\t100\t" in content  # GFF coordinates


class TestFeatureValidatorOutput:
    """Test output generation."""

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

    def test_output_uncompressed(self, temp_dir, output_dir):
        """Test uncompressed output."""
        gff_file = temp_dir / "test.gff"
        with open(gff_file, "w") as f:
            f.write("##gff-version 3\n")
            f.write("chr1\ttest\tgene\t100\t200\t.\t+\t.\tID=gene1\n")

        feature_config = FeatureConfig(
            filename="test.gff",
            filepath=gff_file,
            coding_type=CodingType.NONE,
            detected_format=FeatureFormat.GFF,
            output_dir=output_dir,
            global_options={}
        )

        settings = FeatureValidator.Settings(coding_type=None)
        validator = FeatureValidator(feature_config, settings)
        validator.run()

        output_file = output_dir / "test.gff"
        assert output_file.exists()
        assert output_file.suffix == ".gff"

    def test_output_gzip_compressed(self, temp_dir, output_dir):
        """Test gzip-compressed output."""
        gff_file = temp_dir / "test.gff"
        with open(gff_file, "w") as f:
            f.write("##gff-version 3\n")
            f.write("chr1\ttest\tgene\t100\t200\t.\t+\t.\tID=gene1\n")

        feature_config = FeatureConfig(
            filename="test.gff",
            filepath=gff_file,
            coding_type=CodingType.NONE,
            detected_format=FeatureFormat.GFF,
            output_dir=output_dir,
            global_options={}
        )

        settings = FeatureValidator.Settings(coding_type='gz')
        validator = FeatureValidator(feature_config, settings)
        validator.run()

        output_file = output_dir / "test.gff.gz"
        assert output_file.exists()

        # Verify content is compressed
        with gzip.open(output_file, "rt") as f:
            content = f.read()
            assert "##gff-version 3" in content

    def test_output_bzip2_compressed(self, temp_dir, output_dir):
        """Test bzip2-compressed output."""
        gff_file = temp_dir / "test.gff"
        with open(gff_file, "w") as f:
            f.write("##gff-version 3\n")
            f.write("chr1\ttest\tgene\t100\t200\t.\t+\t.\tID=gene1\n")

        feature_config = FeatureConfig(
            filename="test.gff",
            filepath=gff_file,
            coding_type=CodingType.NONE,
            detected_format=FeatureFormat.GFF,
            output_dir=output_dir,
            global_options={}
        )

        settings = FeatureValidator.Settings(coding_type='bz2')
        validator = FeatureValidator(feature_config, settings)
        validator.run()

        output_file = output_dir / "test.gff.bz2"
        assert output_file.exists()

        # Verify content is compressed
        with bz2.open(output_file, "rt") as f:
            content = f.read()
            assert "##gff-version 3" in content

    def test_output_with_suffix(self, temp_dir, output_dir):
        """Test output filename with custom suffix."""
        gff_file = temp_dir / "test.gff"
        with open(gff_file, "w") as f:
            f.write("##gff-version 3\n")
            f.write("chr1\ttest\tgene\t100\t200\t.\t+\t.\tID=gene1\n")

        feature_config = FeatureConfig(
            filename="test.gff",
            filepath=gff_file,
            coding_type=CodingType.NONE,
            detected_format=FeatureFormat.GFF,
            output_dir=output_dir,
            global_options={}
        )

        settings = FeatureValidator.Settings(output_filename_suffix="validated")
        validator = FeatureValidator(feature_config, settings)
        validator.run()

        output_file = output_dir / "test_validated.gff"
        assert output_file.exists()

    def test_output_in_subdirectory(self, temp_dir, output_dir):
        """Test output in custom subdirectory."""
        gff_file = temp_dir / "test.gff"
        with open(gff_file, "w") as f:
            f.write("##gff-version 3\n")
            f.write("chr1\ttest\tgene\t100\t200\t.\t+\t.\tID=gene1\n")

        feature_config = FeatureConfig(
            filename="test.gff",
            filepath=gff_file,
            coding_type=CodingType.NONE,
            detected_format=FeatureFormat.GFF,
            output_dir=output_dir,
            global_options={}
        )

        settings = FeatureValidator.Settings(output_subdir_name="features")
        validator = FeatureValidator(feature_config, settings)
        validator.run()

        output_file = output_dir / "features" / "test.gff"
        assert output_file.exists()


class TestFeatureValidatorValidationLevels:
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
    def multi_feature_gff(self, temp_dir):
        """Create a GFF file with multiple features."""
        gff_file = temp_dir / "features.gff"
        with open(gff_file, "w") as f:
            f.write("##gff-version 3\n")
            for i in range(1, 21):  # 20 features
                f.write(f"chr1\ttest\tgene\t{i*100}\t{i*100+50}\t.\t+\t.\tID=gene{i}\n")
        return gff_file

    @pytest.fixture
    def damaged_gff(self, temp_dir):
        """Create a GFF file with invalid coordinates in first feature."""
        gff_file = temp_dir / "damaged.gff"
        with open(gff_file, "w") as f:
            f.write("##gff-version 3\n")
            f.write("chr1\ttest\tgene\t200\t100\t.\t+\t.\tID=gene1\n")  # start > end (invalid)
            f.write("chr1\ttest\tgene\t300\t400\t.\t+\t.\tID=gene2\n")  # valid
        return gff_file

    # ===== Tests for STRICT validation level =====

    def test_strict_correct_file_passes(self, multi_feature_gff, output_dir):
        """Test strict mode with correct GFF file - should pass."""
        settings = FeatureValidator.Settings()
        feature_config = FeatureConfig(
            filename="features.gff",
            filepath=multi_feature_gff,
            coding_type=CodingType.NONE,
            detected_format=FeatureFormat.GFF,
            output_dir=output_dir,
            global_options={'validation_level': 'strict'}
        )

        validator = FeatureValidator(feature_config, settings)
        validator.run()

        # All features should be parsed and validated
        assert len(validator.features) == 20
        # Output file should exist
        output_files = list(output_dir.glob("*.gff"))
        assert len(output_files) == 1

    def test_strict_damaged_file_fails(self, damaged_gff, output_dir):
        """Test strict mode with damaged file - should fail."""
        settings = FeatureValidator.Settings(check_coordinates=True
        )
        feature_config = FeatureConfig(
            filename="damaged.gff",
            filepath=damaged_gff,
            coding_type=CodingType.NONE,
            detected_format=FeatureFormat.GFF,
            output_dir=output_dir,
            global_options={'validation_level': 'strict'}
        )

        validator = FeatureValidator(feature_config, settings)

        with pytest.raises(FeatureValidationError, match="start > end"):
            validator.run()

    def test_strict_applies_edits(self, multi_feature_gff, output_dir):
        """Test strict mode applies all edits."""
        settings = FeatureValidator.Settings(sort_by_position=True
        )
        feature_config = FeatureConfig(
            filename="features.gff",
            filepath=multi_feature_gff,
            coding_type=CodingType.NONE,
            detected_format=FeatureFormat.GFF,
            output_dir=output_dir,
            global_options={'validation_level': 'strict'}
        )

        validator = FeatureValidator(feature_config, settings)
        validator.run()

        # Check that sorting was applied
        for i in range(len(validator.features) - 1):
            assert validator.features[i].start <= validator.features[i+1].start

    # ===== Tests for TRUST validation level =====

    def test_trust_correct_file_passes(self, multi_feature_gff, output_dir):
        """Test trust mode with correct GFF file - should pass."""
        settings = FeatureValidator.Settings()
        feature_config = FeatureConfig(
            filename="features.gff",
            filepath=multi_feature_gff,
            coding_type=CodingType.NONE,
            detected_format=FeatureFormat.GFF,
            output_dir=output_dir,
            global_options={'validation_level': 'trust'}
        )

        validator = FeatureValidator(feature_config, settings)
        validator.run()

        # All features should be parsed (20), only first 10 validated
        assert len(validator.features) == 20
        # Output file should exist with all features
        output_files = list(output_dir.glob("*.gff"))
        assert len(output_files) == 1

    def test_trust_damaged_first_feature_fails(self, damaged_gff, output_dir):
        """Test trust mode detects error in first feature."""
        settings = FeatureValidator.Settings(check_coordinates=True
        )
        feature_config = FeatureConfig(
            filename="damaged.gff",
            filepath=damaged_gff,
            coding_type=CodingType.NONE,
            detected_format=FeatureFormat.GFF,
            output_dir=output_dir,
            global_options={'validation_level': 'trust'}
        )

        validator = FeatureValidator(feature_config, settings)

        # Should fail because first feature has invalid coordinates
        with pytest.raises(FeatureValidationError, match="start > end"):
            validator.run()

    def test_trust_applies_edits_to_all(self, multi_feature_gff, output_dir):
        """Test trust mode applies edits to all features."""
        settings = FeatureValidator.Settings(sort_by_position=True
        )
        feature_config = FeatureConfig(
            filename="features.gff",
            filepath=multi_feature_gff,
            coding_type=CodingType.NONE,
            detected_format=FeatureFormat.GFF,
            output_dir=output_dir,
            global_options={'validation_level': 'trust'}
        )

        validator = FeatureValidator(feature_config, settings)
        validator.run()

        # All 20 features should have edits applied
        assert len(validator.features) == 20
        # Check that sorting was applied to ALL features
        for i in range(len(validator.features) - 1):
            assert validator.features[i].start <= validator.features[i+1].start

    # ===== Tests for MINIMAL validation level =====

    def test_minimal_correct_file_passes(self, multi_feature_gff, output_dir):
        """Test minimal mode with correct file - should pass without validation."""
        settings = FeatureValidator.Settings()
        feature_config = FeatureConfig(
            filename="features.gff",
            filepath=multi_feature_gff,
            coding_type=CodingType.NONE,
            detected_format=FeatureFormat.GFF,
            output_dir=output_dir,
            global_options={'validation_level': 'minimal'}
        )

        validator = FeatureValidator(feature_config, settings)
        validator.run()

        # No features parsed in minimal mode
        assert len(validator.features) == 0
        # Output file should exist (copy)
        output_files = list(output_dir.glob("*.gff"))
        assert len(output_files) == 1

    def test_minimal_damaged_file_passes(self, damaged_gff, output_dir):
        """Test minimal mode with damaged file - should pass (no validation)."""
        settings = FeatureValidator.Settings(check_coordinates=True  # Ignored in minimal mode
        )
        feature_config = FeatureConfig(
            filename="damaged.gff",
            filepath=damaged_gff,
            coding_type=CodingType.NONE,
            detected_format=FeatureFormat.GFF,
            output_dir=output_dir,
            global_options={'validation_level': 'minimal'}
        )

        validator = FeatureValidator(feature_config, settings)
        # Should NOT raise - minimal mode doesn't validate
        validator.run()

        # Output should be created
        output_files = list(output_dir.glob("*.gff"))
        assert len(output_files) == 1

    def test_minimal_output_is_copy(self, multi_feature_gff, output_dir):
        """Test that minimal mode copies file as-is."""
        settings = FeatureValidator.Settings()
        feature_config = FeatureConfig(
            filename="features.gff",
            filepath=multi_feature_gff,
            coding_type=CodingType.NONE,
            detected_format=FeatureFormat.GFF,
            output_dir=output_dir,
            global_options={'validation_level': 'minimal'}
        )

        validator = FeatureValidator(feature_config, settings)
        validator.run()

        # Check output file exists
        output_files = list(output_dir.glob("*.gff"))
        assert len(output_files) == 1

        # Verify output is byte-for-byte identical to input
        with open(multi_feature_gff, 'rb') as f_in, open(output_files[0], 'rb') as f_out:
            assert f_in.read() == f_out.read()

    def test_minimal_does_not_apply_edits(self, multi_feature_gff, output_dir):
        """Test minimal mode does NOT apply edits."""
        settings = FeatureValidator.Settings(sort_by_position=True  # Should be ignored
        )
        feature_config = FeatureConfig(
            filename="features.gff",
            filepath=multi_feature_gff,
            coding_type=CodingType.NONE,
            detected_format=FeatureFormat.GFF,
            output_dir=output_dir,
            global_options={'validation_level': 'minimal'}
        )

        validator = FeatureValidator(feature_config, settings)
        validator.run()

        # Output should be identical to input (no edits applied)
        output_file = list(output_dir.glob("*.gff"))[0]
        with open(multi_feature_gff, 'r') as f_in, open(output_file, 'r') as f_out:
            assert f_in.read() == f_out.read()

    # ===== Test for compressed files =====

    def test_trust_compressed_gz_passes(self, temp_dir, output_dir):
        """Test trust mode with gzip compressed file."""
        gff_file = temp_dir / "features.gff.gz"
        with gzip.open(gff_file, "wt") as f:
            f.write("##gff-version 3\n")
            for i in range(1, 11):
                f.write(f"chr1\ttest\tgene\t{i*100}\t{i*100+50}\t.\t+\t.\tID=gene{i}\n")

        settings = FeatureValidator.Settings()
        feature_config = FeatureConfig(
            filename="features.gff.gz",
            filepath=gff_file,
            coding_type=CodingType.GZIP,
            detected_format=FeatureFormat.GFF,
            output_dir=output_dir,
            global_options={'validation_level': 'trust'}
        )

        validator = FeatureValidator(feature_config, settings)
        validator.run()

        # All features should be parsed
        assert len(validator.features) == 10

    def test_minimal_compressed_bz2_raises_error(self, temp_dir, output_dir):
        """Test minimal mode with bzip2 compressed file - should raise error."""
        gff_file = temp_dir / "features.gff.bz2"
        with bz2.open(gff_file, "wt") as f:
            f.write("##gff-version 3\n")
            f.write("chr1\ttest\tgene\t100\t200\t.\t+\t.\tID=gene1\n")

        settings = FeatureValidator.Settings()
        feature_config = FeatureConfig(
            filename="features.gff.bz2",
            filepath=gff_file,
            coding_type=CodingType.BZIP2,
            detected_format=FeatureFormat.GFF,
            output_dir=output_dir,
            global_options={'validation_level': 'minimal'}
        )

        validator = FeatureValidator(feature_config, settings)

        # Should raise error - minimal mode doesn't support compressed files
        with pytest.raises(FeatureValidationError):
            validator.run()
class TestFeatureValidatorEdgeCases:
    """Test edge cases and error handling."""

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

    def test_bed_with_minimal_columns(self, temp_dir, output_dir):
        """Test BED file with only required 3 columns."""
        bed_file = temp_dir / "minimal.bed"
        with open(bed_file, "w") as f:
            f.write("chr1\t0\t100\n")  # Only 3 columns

        feature_config = FeatureConfig(
            filename="minimal.bed",
            filepath=bed_file,
            coding_type=CodingType.NONE,
            detected_format=FeatureFormat.BED,
            output_dir=output_dir,
            global_options={}
        )

        validator = FeatureValidator(feature_config, FeatureValidator.Settings())
        validator.run()

        assert len(validator.features) == 1
        assert validator.features[0].seqname == "chr1"
        assert validator.features[0].attributes == ""
        assert validator.features[0].score == "."
        assert validator.features[0].strand == "."

    def test_gff_with_comments(self, temp_dir, output_dir):
        """Test GFF file with comments and empty lines."""
        gff_file = temp_dir / "comments.gff"
        with open(gff_file, "w") as f:
            f.write("##gff-version 3\n")
            f.write("# This is a comment\n")
            f.write("\n")
            f.write("chr1\ttest\tgene\t100\t200\t.\t+\t.\tID=gene1\n")
            f.write("\n")
            f.write("# Another comment\n")
            f.write("chr1\ttest\tCDS\t120\t180\t.\t+\t0\tID=cds1\n")

        feature_config = FeatureConfig(
            filename="comments.gff",
            filepath=gff_file,
            coding_type=CodingType.NONE,
            detected_format=FeatureFormat.GFF,
            output_dir=output_dir,
            global_options={}
        )

        validator = FeatureValidator(feature_config, FeatureValidator.Settings())
        validator.run()

        # Should skip comments and empty lines
        assert len(validator.features) == 2

    def test_strand_validation(self, temp_dir, output_dir):
        """Test handling of invalid strand values."""
        gff_file = temp_dir / "strand.gff"
        with open(gff_file, "w") as f:
            f.write("##gff-version 3\n")
            f.write("chr1\ttest\tgene\t100\t200\t.\t+\t.\tID=gene1\n")
            f.write("chr1\ttest\tgene\t300\t400\t.\t-\t.\tID=gene2\n")
            f.write("chr1\ttest\tgene\t500\t600\t.\t.\t.\tID=gene3\n")
            f.write("chr1\ttest\tgene\t700\t800\t.\tX\t.\tID=gene4\n")  # Invalid strand

        feature_config = FeatureConfig(
            filename="strand.gff",
            filepath=gff_file,
            coding_type=CodingType.NONE,
            detected_format=FeatureFormat.GFF,
            output_dir=output_dir,
            global_options={}
        )

        validator = FeatureValidator(feature_config, FeatureValidator.Settings())
        validator.run()

        assert len(validator.features) == 4
        assert validator.features[0].strand == "+"
        assert validator.features[1].strand == "-"
        assert validator.features[2].strand == "."
        assert validator.features[3].strand == "."  # Invalid converted to "."
