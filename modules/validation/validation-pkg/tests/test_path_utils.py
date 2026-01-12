"""
Tests for path utility functions.

This module tests path resolution and security utilities in the path_utils module,
including path traversal protection.
"""

import pytest
import tempfile
from pathlib import Path

from validation_pkg.utils.path_utils import *
from validation_pkg.exceptions import ConfigurationError


class TestResolveFilepath:
    """Test path resolution with security protections."""

    @pytest.fixture
    def temp_dir(self):
        """Create a temporary directory for test files."""
        with tempfile.TemporaryDirectory() as tmpdir:
            yield Path(tmpdir)

    def test_resolve_simple_filename(self, temp_dir):
        """Test resolving a simple filename in base directory."""
        filepath = resolve_filepath(temp_dir, "genome.fasta")
        assert filepath == temp_dir / "genome.fasta"

    def test_resolve_relative_path(self, temp_dir):
        """Test resolving a relative path within base directory."""
        filepath = resolve_filepath(temp_dir, "subdir/genome.fasta")
        assert filepath == temp_dir / "subdir" / "genome.fasta"

    def test_resolve_with_dot_notation(self, temp_dir):
        """Test resolving path with ./ notation."""
        filepath = resolve_filepath(temp_dir, "./genome.fasta")
        assert filepath == temp_dir / "genome.fasta"

    def test_path_traversal_with_dotdot_blocked(self, temp_dir):
        """Test that ../ path traversal is blocked."""
        with pytest.raises(ConfigurationError, match="Path traversal detected"):
            resolve_filepath(temp_dir, "../../../etc/passwd")

    def test_path_traversal_complex_blocked(self, temp_dir):
        """Test that complex path traversal is blocked."""
        with pytest.raises(ConfigurationError, match="Path traversal detected"):
            resolve_filepath(temp_dir, "subdir/../../etc/passwd")

    def test_path_traversal_absolute_outside_blocked(self, temp_dir):
        """Test that absolute paths outside base directory are blocked."""
        with pytest.raises(ConfigurationError, match="Path traversal detected"):
            resolve_filepath(temp_dir, "/etc/passwd")

    def test_symlink_outside_base_blocked(self, temp_dir):
        """Test that symlinks pointing outside base directory are blocked."""
        # Create a symlink pointing outside the base directory
        outside_file = Path("/etc/passwd")
        symlink_path = temp_dir / "malicious_symlink"

        # Only run this test on systems where /etc/passwd exists
        if outside_file.exists():
            symlink_path.symlink_to(outside_file)

            with pytest.raises(ConfigurationError, match="Path traversal detected"):
                resolve_filepath(temp_dir, "malicious_symlink")

    def test_valid_symlink_within_base_allowed(self, temp_dir):
        """Test that symlinks within base directory are allowed."""
        # Create a real file
        real_file = temp_dir / "genome.fasta"
        real_file.write_text(">seq1\nATCG\n")

        # Create a symlink to it within the same directory
        symlink_path = temp_dir / "genome_link.fasta"
        symlink_path.symlink_to(real_file)

        # Should resolve successfully (symlinks are resolved to their targets)
        filepath = resolve_filepath(temp_dir, "genome_link.fasta")
        # Symlink resolves to the real file path
        assert filepath == temp_dir / "genome.fasta"
        assert filepath.exists()

    def test_path_with_multiple_slashes(self, temp_dir):
        """Test that paths with multiple slashes are normalized."""
        filepath = resolve_filepath(temp_dir, "subdir//genome.fasta")
        assert filepath == temp_dir / "subdir" / "genome.fasta"

    def test_error_message_includes_details(self, temp_dir):
        """Test that error message includes helpful details."""
        try:
            resolve_filepath(temp_dir, "../../../etc/passwd")
            pytest.fail("Should have raised ConfigurationError")
        except ConfigurationError as e:
            error_msg = str(e)
            assert "Path traversal detected" in error_msg
            assert "resolves outside config directory" in error_msg
            assert str(temp_dir) in error_msg


class TestResolveFilepathEdgeCases:
    """Test edge cases for path resolution."""

    @pytest.fixture
    def temp_dir(self):
        """Create a temporary directory for test files."""
        with tempfile.TemporaryDirectory() as tmpdir:
            yield Path(tmpdir)

    def test_empty_filename(self, temp_dir):
        """Test resolving empty filename."""
        filepath = resolve_filepath(temp_dir, "")
        assert filepath == temp_dir

    def test_dot_as_filename(self, temp_dir):
        """Test resolving . as filename."""
        filepath = resolve_filepath(temp_dir, ".")
        assert filepath == temp_dir

    def test_unicode_filename(self, temp_dir):
        """Test resolving filename with unicode characters."""
        filepath = resolve_filepath(temp_dir, "genome_\u03B1\u03B2.fasta")
        assert filepath == temp_dir / "genome_\u03B1\u03B2.fasta"

    def test_filename_with_spaces(self, temp_dir):
        """Test resolving filename with spaces."""
        filepath = resolve_filepath(temp_dir, "my genome file.fasta")
        assert filepath == temp_dir / "my genome file.fasta"

    def test_deeply_nested_path(self, temp_dir):
        """Test resolving deeply nested path within base directory."""
        deep_path = "a/b/c/d/e/f/genome.fasta"
        filepath = resolve_filepath(temp_dir, deep_path)
        assert filepath == temp_dir / "a" / "b" / "c" / "d" / "e" / "f" / "genome.fasta"


class TestPathIncrement:
    """Test suite for get_incremented_path() function."""

    def test_increment_nonexistent_file(self):
        """Test that original path is returned if file doesn't exist."""
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "report.txt"
            result = get_incremented_path(path)
            assert result == path
            assert str(result) == str(path)

    def test_increment_existing_file(self):
        """Test that _001 suffix is added if file exists."""
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "report.txt"
            path.write_text("test")

            result = get_incremented_path(path)
            assert result.name == "report_001.txt"
            assert result.parent == path.parent
            assert not result.exists()

    def test_increment_multiple_times(self):
        """Test correct incrementing: 001 → 002 → 003."""
        with tempfile.TemporaryDirectory() as tmpdir:
            base_path = Path(tmpdir) / "report.txt"
            base_path.write_text("test")

            # First increment
            path1 = get_incremented_path(base_path)
            assert path1.name == "report_001.txt"
            path1.write_text("test1")

            # Second increment from base
            path2 = get_incremented_path(base_path)
            assert path2.name == "report_002.txt"
            path2.write_text("test2")

            # Third increment from base
            path3 = get_incremented_path(base_path)
            assert path3.name == "report_003.txt"

    def test_increment_preserves_extension(self):
        """Test that file extension is preserved correctly."""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Test .txt
            txt_path = Path(tmpdir) / "report.txt"
            txt_path.write_text("test")
            result = get_incremented_path(txt_path)
            assert result.suffix == ".txt"
            assert result.name == "report_001.txt"

            # Test .log
            log_path = Path(tmpdir) / "validation.log"
            log_path.write_text("test")
            result = get_incremented_path(log_path)
            assert result.suffix == ".log"
            assert result.name == "validation_001.log"

    def test_increment_with_existing_number(self):
        """Test incrementing a file that already has a number."""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Create report_001.txt
            path1 = Path(tmpdir) / "report_001.txt"
            path1.write_text("test")

            # Incrementing it should give report_002.txt
            result = get_incremented_path(path1)
            assert result.name == "report_002.txt"

    def test_increment_with_custom_separator(self):
        """Test using a custom separator."""
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "report.txt"
            path.write_text("test")

            result = get_incremented_path(path, separator="-")
            assert result.name == "report-001.txt"

    def test_increment_no_extension(self):
        """Test incrementing a file without extension."""
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "report"
            path.write_text("test")

            result = get_incremented_path(path)
            assert result.name == "report_001"

    def test_increment_multiple_dots(self):
        """Test with files that have multiple dots in name."""
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "my.report.txt"
            path.write_text("test")

            result = get_incremented_path(path)
            # Should preserve stem and only last extension
            assert result.name == "my.report_001.txt"

    def test_increment_safety_limit(self):
        """Test that function raises error after reaching counter limit."""
        import unittest.mock as mock

        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "report.txt"
            path.write_text("test")

            # Mock Path.exists to always return True, simulating all files exist
            # This will cause the counter to keep incrementing until it hits 10000
            original_exists = Path.exists
            def mock_exists(self):
                # Allow the original path to exist, but all numbered paths also exist
                return True

            with mock.patch.object(Path, 'exists', mock_exists):
                with pytest.raises(RuntimeError, match="Too many incremented files"):
                    get_incremented_path(path)
