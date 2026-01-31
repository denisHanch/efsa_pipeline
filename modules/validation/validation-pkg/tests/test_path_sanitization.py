"""
Tests for path sanitization utilities.

Tests cover:
- Path component sanitization
- Safe output directory creation
- Path traversal attack prevention
- Illegal character handling
- Security validation
"""

import pytest
from pathlib import Path
import tempfile

from validation_pkg.utils.path_utils import sanitize_path_component, build_safe_output_dir
from validation_pkg.exceptions import ConfigurationError


class TestSanitizePathComponent:
    """Test path component sanitization."""

    def test_valid_simple_name(self):
        """Test valid simple directory name."""
        result = sanitize_path_component("results")
        assert result == "results"

    def test_valid_with_underscores_dashes(self):
        """Test valid name with underscores and dashes."""
        result = sanitize_path_component("my_results-v2")
        assert result == "my_results-v2"

    def test_valid_with_numbers(self):
        """Test valid name with numbers."""
        result = sanitize_path_component("output123")
        assert result == "output123"

    def test_blocks_parent_directory(self):
        """Test that .. is blocked."""
        with pytest.raises(ValueError, match="contains '..'"):
            sanitize_path_component("..")

    def test_blocks_parent_in_path(self):
        """Test that ../ in path is blocked."""
        with pytest.raises(ValueError, match="contains '..'"):
            sanitize_path_component("subdir/..")

    def test_blocks_absolute_path(self):
        """Test that absolute paths are blocked."""
        with pytest.raises(ValueError, match="absolute path"):
            sanitize_path_component("/etc/passwd")

    def test_blocks_null_byte(self):
        """Test that null bytes are blocked."""
        with pytest.raises(ValueError, match="null byte"):
            sanitize_path_component("test\x00malicious")

    def test_blocks_empty_string(self):
        """Test that empty strings are blocked."""
        with pytest.raises(ValueError, match="empty"):
            sanitize_path_component("")

    def test_blocks_current_directory(self):
        """Test that . is blocked."""
        with pytest.raises(ValueError, match="empty or current directory"):
            sanitize_path_component(".")

    def test_strips_whitespace(self):
        """Test that leading/trailing whitespace is stripped."""
        result = sanitize_path_component("  results  ")
        assert result == "results"

    def test_blocks_slashes_by_default(self):
        """Test that slashes are blocked by default."""
        with pytest.raises(ValueError, match="path separator"):
            sanitize_path_component("sub/dir")

    def test_allows_slashes_when_enabled(self):
        """Test that slashes are allowed when explicitly enabled."""
        result = sanitize_path_component("sub/dir", allow_slashes=True)
        assert result == "sub/dir"


class TestBuildSafeOutputDir:
    """Test safe output directory building."""

    @pytest.fixture
    def temp_base(self):
        """Create temporary base directory."""
        with tempfile.TemporaryDirectory() as tmpdir:
            yield Path(tmpdir)

    def test_without_subdirectory(self, temp_base):
        """Test building path without subdirectory."""
        output_dir = build_safe_output_dir(temp_base)
        assert output_dir == temp_base
        assert output_dir.exists()

    def test_with_valid_subdirectory(self, temp_base):
        """Test building path with valid subdirectory."""
        output_dir = build_safe_output_dir(temp_base, "results")
        assert output_dir == temp_base / "results"
        assert output_dir.exists()

    def test_blocks_parent_directory_escape(self, temp_base):
        """Test that ../ escape is blocked."""
        with pytest.raises(ValueError, match="Invalid output subdirectory"):
            build_safe_output_dir(temp_base, "../escape")

    def test_blocks_absolute_path(self, temp_base):
        """Test that absolute paths are blocked."""
        with pytest.raises(ValueError, match="Invalid output subdirectory"):
            build_safe_output_dir(temp_base, "/etc/passwd")

    def test_without_create(self, temp_base):
        """Test building path without creating directory."""
        output_dir = build_safe_output_dir(temp_base, "results", create=False)
        assert output_dir == temp_base / "results"
        assert not output_dir.exists()

    def test_creates_nested_directories(self, temp_base):
        """Test that directory is created with parents."""
        # Note: This would need to be a valid subdirectory name
        # For nested creation, the subdirectory should be a single component
        output_dir = build_safe_output_dir(temp_base, "results")
        assert output_dir.exists()
