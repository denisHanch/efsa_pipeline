"""Tests for validation helper utilities."""

import pytest
from pathlib import Path
import tempfile

from validation_pkg.utils.validation_helpers import (
    validate_minimal_mode_requirements,
    copy_file_minimal_mode
)
from validation_pkg.utils.formats import GenomeFormat, CodingType
from validation_pkg.exceptions import ValidationError
from validation_pkg.logger import get_logger


class TestValidateMinimalModeRequirements:
    """Test minimal mode validation."""

    def test_valid_format_and_coding(self):
        """Test that matching format and coding passes."""
        logger = get_logger()
        # Should not raise
        validate_minimal_mode_requirements(
            GenomeFormat.FASTA,
            GenomeFormat.FASTA,
            CodingType.GZIP,
            CodingType.GZIP,
            "genome.fasta.gz",
            logger,
            "genome"
        )

    def test_format_mismatch_raises(self):
        """Test that format mismatch raises error."""
        logger = get_logger()
        with pytest.raises(ValidationError, match="requires.*FASTA.*format"):
            validate_minimal_mode_requirements(
                GenomeFormat.GENBANK,
                GenomeFormat.FASTA,
                CodingType.GZIP,
                CodingType.GZIP,
                "genome.gb.gz",
                logger,
                "genome"
            )

    def test_coding_mismatch_raises(self):
        """Test that coding mismatch raises error."""
        logger = get_logger()
        with pytest.raises(ValidationError, match="coding to match output coding"):
            validate_minimal_mode_requirements(
                GenomeFormat.FASTA,
                GenomeFormat.FASTA,
                CodingType.NONE,
                CodingType.GZIP,
                "genome.fasta",
                logger,
                "genome"
            )


class TestCopyFileMinimalMode:
    """Test file copying for minimal mode."""

    def test_copy_file(self):
        """Test basic file copying."""
        logger = get_logger()

        with tempfile.TemporaryDirectory() as tmpdir:
            # Create source file
            src = Path(tmpdir) / "input.txt"
            src.write_text("test content")

            # Copy to destination
            dst = Path(tmpdir) / "output.txt"
            result = copy_file_minimal_mode(src, dst, logger)

            assert result == dst
            assert dst.exists()
            assert dst.read_text() == "test content"
