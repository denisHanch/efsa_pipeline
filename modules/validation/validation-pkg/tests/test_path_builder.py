"""Tests for path builder utilities."""

import pytest
from pathlib import Path
import tempfile

from validation_pkg.utils.file_handler import (
    strip_all_extensions,
    build_output_filename,
    build_output_path
)
from validation_pkg.utils.formats import CodingType
from validation_pkg.exceptions import ConfigurationError


class TestStripAllExtensions:
    """Test extension stripping."""

    def test_single_extension(self):
        result = strip_all_extensions("genome.fasta")
        assert result == "genome"

    def test_double_extension(self):
        result = strip_all_extensions("genome.fasta.gz")
        assert result == "genome"

    def test_no_extension(self):
        result = strip_all_extensions("genome")
        assert result == "genome"

    def test_with_path_object(self):
        path = Path("genome.fasta.gz")
        result = strip_all_extensions("genome.fasta.gz", path)
        assert result == "genome"

    def test_complex_name(self):
        result = strip_all_extensions("my.genome.v2.fasta.bz2")
        assert result == "my"


class TestBuildOutputFilename:
    """Test output filename building."""

    def test_basic_filename(self):
        result = build_output_filename("genome.fa", "fasta")
        assert result == "genome.fasta"

    def test_with_compression(self):
        result = build_output_filename("genome.fa", "fasta", CodingType.GZIP)
        assert result == "genome.fasta.gz"

    def test_with_suffix(self):
        result = build_output_filename("genome.fa", "fasta", suffix="validated")
        assert result == "genome_validated.fasta"

    def test_with_suffix_and_compression(self):
        result = build_output_filename(
            "genome.fa.gz",
            "fasta",
            CodingType.GZIP,
            "filtered"
        )
        assert result == "genome_filtered.fasta.gz"


class TestBuildOutputPath:
    """Test complete output path building."""

    @pytest.fixture
    def temp_base(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            yield Path(tmpdir)

    def test_basic_path(self, temp_base):
        output_path = build_output_path(
            temp_base,
            "genome.fasta",
            "fasta"
        )
        assert output_path == temp_base / "genome.fasta"

    def test_with_subdirectory(self, temp_base):
        output_path = build_output_path(
            temp_base,
            "genome.fasta",
            "fasta",
            subdir_name="validated"
        )
        assert output_path == temp_base / "validated" / "genome.fasta"
        assert output_path.parent.exists()

    def test_blocks_path_traversal(self, temp_base):
        with pytest.raises(ValueError):
            build_output_path(
                temp_base,
                "genome.fasta",
                "fasta",
                subdir_name="../escape"
            )

    def test_full_features(self, temp_base):
        output_path = build_output_path(
            temp_base,
            "genome.fa.gz",
            "fasta",
            CodingType.GZIP,
            subdir_name="results",
            filename_suffix="filtered"
        )
        expected = temp_base / "results" / "genome_filtered.fasta.gz"
        assert output_path == expected
