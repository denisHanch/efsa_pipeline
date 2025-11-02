"""
Tests for path utility functions.

This module tests path resolution and security utilities in the path_utils module,
including path traversal protection.
"""

import pytest
import tempfile
from pathlib import Path

from validation_pkg.utils.path_utils import resolve_filepath
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
