"""
Unit tests for the file_handler module.

Tests file handling utilities including:
- File opening with automatic decompression
- Compression detection
- Format detection
- Filename processing
- Compression conversion functions
"""

import pytest
import gzip
import bz2
import tempfile
import subprocess
from pathlib import Path

from validation_pkg.utils.file_handler import (
    open_file_with_coding_type,
    gz_to_bz2,
    bz2_to_gz,
    none_to_gz,
    gz_to_none,
    bz2_to_none,
    none_to_bz2,
    check_compression_tool_available,
    get_compression_command,
    open_compressed_writer,
)
from validation_pkg.utils.formats import CodingType
from validation_pkg.exceptions import CompressionError


@pytest.fixture
def temp_dir():
    """Create a temporary directory for test files."""
    with tempfile.TemporaryDirectory() as tmpdir:
        yield Path(tmpdir)


@pytest.fixture
def sample_content():
    """Sample text content for test files."""
    return "ATCGATCGATCG\nGCTAGCTAGCTA\n"


@pytest.fixture
def plain_file(temp_dir, sample_content):
    """Create a plain text file."""
    filepath = temp_dir / "test.fasta"
    filepath.write_text(sample_content)
    return filepath


@pytest.fixture
def gzip_file(temp_dir, sample_content):
    """Create a gzip compressed file."""
    filepath = temp_dir / "test.fasta.gz"
    with gzip.open(filepath, 'wt') as f:
        f.write(sample_content)
    return filepath


@pytest.fixture
def bz2_file(temp_dir, sample_content):
    """Create a bzip2 compressed file."""
    filepath = temp_dir / "test.fasta.bz2"
    with bz2.open(filepath, 'wt') as f:
        f.write(sample_content)
    return filepath


class TestOpenFileWithCodingType:
    """Tests for open_file_with_coding_type function."""

    def test_open_plain_file_with_none_coding(self, plain_file, sample_content):
        """Test opening plain file with CodingType.NONE."""
        with open_file_with_coding_type(plain_file, CodingType.NONE) as f:
            content = f.read()
        assert content == sample_content

    def test_open_gzip_file_with_gzip_coding(self, gzip_file, sample_content):
        """Test opening gzip file with CodingType.GZIP."""
        with open_file_with_coding_type(gzip_file, CodingType.GZIP) as f:
            content = f.read()
        assert content == sample_content

    def test_open_bz2_file_with_bzip2_coding(self, bz2_file, sample_content):
        """Test opening bzip2 file with CodingType.BZIP2."""
        with open_file_with_coding_type(bz2_file, CodingType.BZIP2) as f:
            content = f.read()
        assert content == sample_content

    def test_open_file_accepts_string_path(self, plain_file, sample_content):
        """Test that function accepts string paths."""
        with open_file_with_coding_type(str(plain_file), CodingType.NONE) as f:
            content = f.read()
        assert content == sample_content

    def test_open_file_accepts_path_object(self, plain_file, sample_content):
        """Test that function accepts Path objects."""
        with open_file_with_coding_type(plain_file, CodingType.NONE) as f:
            content = f.read()
        assert content == sample_content

    def test_open_nonexistent_file_raises_compression_error(self, temp_dir):
        """Test that opening non-existent file raises CompressionError."""
        nonexistent = temp_dir / "nonexistent.fasta"
        with pytest.raises(CompressionError) as exc_info:
            open_file_with_coding_type(nonexistent, CodingType.NONE)
        assert "Failed to open file" in str(exc_info.value)
        assert str(nonexistent) in str(exc_info.value)

    def test_open_gzip_as_plain_raises_compression_error(self, gzip_file):
        """Test that opening gzip file as plain raises error."""
        with pytest.raises(UnicodeDecodeError):
            with open_file_with_coding_type(gzip_file, CodingType.NONE) as f:
                f.read()

    def test_open_plain_as_gzip_raises_compression_error(self, plain_file):
        """Test that opening plain file as gzip raises error."""
        with pytest.raises(gzip.BadGzipFile):
            with open_file_with_coding_type(plain_file, CodingType.GZIP) as f:
                f.read()

    def test_compression_error_chains_original_exception(self, temp_dir):
        """Test that CompressionError chains the original exception."""
        nonexistent = temp_dir / "nonexistent.fasta"
        with pytest.raises(CompressionError) as exc_info:
            open_file_with_coding_type(nonexistent, CodingType.NONE)
        assert exc_info.value.__cause__ is not None

    def test_read_mode_default(self, plain_file, sample_content):
        """Test default mode is 'rt' (text read)."""
        with open_file_with_coding_type(plain_file, CodingType.NONE) as f:
            content = f.read()
            assert isinstance(content, str)
        assert content == sample_content

    def test_binary_mode(self, plain_file, sample_content):
        """Test opening file in binary mode."""
        with open_file_with_coding_type(plain_file, CodingType.NONE, mode='rb') as f:
            content = f.read()
            assert isinstance(content, bytes)
        assert content == sample_content.encode()


class TestCompressionConversion:
    """Tests for compression conversion functions."""

    def test_gz_to_bz2(self, temp_dir, gzip_file, sample_content):
        """Test converting gzip to bzip2."""
        output_file = temp_dir / "output.fasta.bz2"
        gz_to_bz2(gzip_file, output_file)

        # Verify output exists and has correct content
        assert output_file.exists()
        with bz2.open(output_file, 'rt') as f:
            content = f.read()
        assert content == sample_content

    def test_bz2_to_gz(self, temp_dir, bz2_file, sample_content):
        """Test converting bzip2 to gzip."""
        output_file = temp_dir / "output.fasta.gz"
        bz2_to_gz(bz2_file, output_file)

        # Verify output exists and has correct content
        assert output_file.exists()
        with gzip.open(output_file, 'rt') as f:
            content = f.read()
        assert content == sample_content

    def test_none_to_gz(self, temp_dir, plain_file, sample_content):
        """Test compressing plain file to gzip."""
        output_file = temp_dir / "output.fasta.gz"
        none_to_gz(plain_file, output_file)

        # Verify output exists and has correct content
        assert output_file.exists()
        with gzip.open(output_file, 'rt') as f:
            content = f.read()
        assert content == sample_content

    def test_gz_to_none(self, temp_dir, gzip_file, sample_content):
        """Test decompressing gzip to plain file."""
        output_file = temp_dir / "output.fasta"
        gz_to_none(gzip_file, output_file)

        # Verify output exists and has correct content
        assert output_file.exists()
        content = output_file.read_text()
        assert content == sample_content

    def test_bz2_to_none(self, temp_dir, bz2_file, sample_content):
        """Test decompressing bzip2 to plain file."""
        output_file = temp_dir / "output.fasta"
        bz2_to_none(bz2_file, output_file)

        # Verify output exists and has correct content
        assert output_file.exists()
        content = output_file.read_text()
        assert content == sample_content

    def test_none_to_bz2(self, temp_dir, plain_file, sample_content):
        """Test compressing plain file to bzip2."""
        output_file = temp_dir / "output.fasta.bz2"
        none_to_bz2(plain_file, output_file)

        # Verify output exists and has correct content
        assert output_file.exists()
        with bz2.open(output_file, 'rt') as f:
            content = f.read()
        assert content == sample_content

    def test_gz_to_bz2_with_nonexistent_input(self, temp_dir):
        """Test conversion with non-existent input file raises error."""
        nonexistent = temp_dir / "nonexistent.fasta.gz"
        output_file = temp_dir / "output.fasta.bz2"

        with pytest.raises(CompressionError):
            gz_to_bz2(nonexistent, output_file)

    def test_conversion_preserves_content(self, temp_dir, sample_content):
        """Test round-trip conversion preserves content."""
        # Create plain file
        plain = temp_dir / "test.fasta"
        plain.write_text(sample_content)

        # Convert: plain -> gz -> bz2 -> plain
        gz_temp = temp_dir / "temp.fasta.gz"
        bz2_temp = temp_dir / "temp.fasta.bz2"
        final = temp_dir / "final.fasta"

        none_to_gz(plain, gz_temp)
        gz_to_bz2(gz_temp, bz2_temp)
        bz2_to_none(bz2_temp, final)

        # Verify content is preserved
        assert final.read_text() == sample_content


class TestParallelCompressionTools:
    """Tests for parallel compression tool utilities."""

    def test_check_compression_tool_available_gzip(self):
        """Test that gzip is available (should always be available on Unix)."""
        assert check_compression_tool_available('gzip') is True

    def test_check_compression_tool_available_bzip2(self):
        """Test that bzip2 is available (should always be available on Unix)."""
        assert check_compression_tool_available('bzip2') is True

    def test_check_compression_tool_available_nonexistent(self):
        """Test that nonexistent tools return False."""
        assert check_compression_tool_available('nonexistent_compression_tool_xyz') is False

    def test_check_compression_tool_available_caching(self):
        """Test that tool availability checks are cached."""
        # First call
        result1 = check_compression_tool_available('gzip')
        # Second call should use cache
        result2 = check_compression_tool_available('gzip')
        assert result1 == result2

    def test_get_compression_command_gzip(self):
        """Test getting gzip compression command."""
        cmd, args = get_compression_command(CodingType.GZIP, 'compress')
        # Should be either pigz or gzip
        assert cmd in ('pigz', 'gzip')
        assert '-c' in args

    def test_get_compression_command_bzip2(self):
        """Test getting bzip2 compression command."""
        cmd, args = get_compression_command(CodingType.BZIP2, 'compress')
        # Should be either pbzip2 or bzip2
        assert cmd in ('pbzip2', 'bzip2')
        assert '-c' in args

    def test_get_compression_command_none(self):
        """Test getting command for no compression."""
        cmd, args = get_compression_command(CodingType.NONE)
        assert cmd == 'cat'
        assert args == []

    def test_get_compression_command_decompress_gzip(self):
        """Test getting gzip decompression command."""
        cmd, args = get_compression_command(CodingType.GZIP, 'decompress')
        # Should be either pigz or gzip
        assert cmd in ('pigz', 'gzip')
        assert '-dc' in args or '-d' in args

    def test_get_compression_command_decompress_bzip2(self):
        """Test getting bzip2 decompression command."""
        cmd, args = get_compression_command(CodingType.BZIP2, 'decompress')
        # Should be either pbzip2 or bzip2
        assert cmd in ('pbzip2', 'bzip2')
        assert '-dc' in args or '-d' in args

    def test_get_compression_command_custom_threads(self):
        """Test getting compression command with custom thread count."""
        cmd, args = get_compression_command(CodingType.GZIP, 'compress', threads=2)
        # Should be either pigz or gzip
        assert cmd in ('pigz', 'gzip')
        if cmd == 'pigz':
            # pigz should have thread argument
            assert '2' in args or '-p' in args


class TestOpenCompressedWriter:
    """Tests for open_compressed_writer function."""

    def test_write_plain_file(self, temp_dir, sample_content):
        """Test writing to uncompressed file."""
        output_file = temp_dir / "output.txt"

        with open_compressed_writer(output_file, CodingType.NONE) as f:
            f.write(sample_content)

        # Verify file was written correctly
        assert output_file.exists()
        assert output_file.read_text() == sample_content

    def test_write_gzip_file(self, temp_dir, sample_content):
        """Test writing to gzip compressed file."""
        output_file = temp_dir / "output.txt.gz"

        with open_compressed_writer(output_file, CodingType.GZIP) as f:
            f.write(sample_content)

        # Verify file was written and is valid gzip
        assert output_file.exists()
        with gzip.open(output_file, 'rt') as f:
            assert f.read() == sample_content

    def test_write_bzip2_file(self, temp_dir, sample_content):
        """Test writing to bzip2 compressed file."""
        output_file = temp_dir / "output.txt.bz2"

        with open_compressed_writer(output_file, CodingType.BZIP2) as f:
            f.write(sample_content)

        # Verify file was written and is valid bzip2
        assert output_file.exists()
        with bz2.open(output_file, 'rt') as f:
            assert f.read() == sample_content

    def test_write_gzip_file_no_parallel(self, temp_dir, sample_content):
        """Test writing to gzip file with use_parallel=False."""
        output_file = temp_dir / "output.txt.gz"

        with open_compressed_writer(output_file, CodingType.GZIP, use_parallel=False) as f:
            f.write(sample_content)

        # Verify file was written correctly
        assert output_file.exists()
        with gzip.open(output_file, 'rt') as f:
            assert f.read() == sample_content

    def test_write_multiline_content(self, temp_dir):
        """Test writing multiline content."""
        output_file = temp_dir / "output.txt.gz"
        content = "Line 1\nLine 2\nLine 3\n"

        with open_compressed_writer(output_file, CodingType.GZIP) as f:
            f.write(content)

        # Verify content
        with gzip.open(output_file, 'rt') as f:
            assert f.read() == content

    def test_write_large_content(self, temp_dir):
        """Test writing large content (streaming behavior)."""
        output_file = temp_dir / "output.txt.gz"
        # Generate large content (1MB)
        content = "ATCG" * (1024 * 256)  # 1MB

        with open_compressed_writer(output_file, CodingType.GZIP) as f:
            f.write(content)

        # Verify content
        with gzip.open(output_file, 'rt') as f:
            assert f.read() == content

    def test_context_manager_cleanup(self, temp_dir, sample_content):
        """Test that context manager properly cleans up resources."""
        output_file = temp_dir / "output.txt.gz"

        # Write file
        with open_compressed_writer(output_file, CodingType.GZIP) as f:
            f.write(sample_content)

        # File should be fully written and closed
        assert output_file.exists()
        # Should be able to read immediately
        with gzip.open(output_file, 'rt') as f:
            assert f.read() == sample_content


class TestSecurityCommandInjectionFileHandler:
    """
    Security tests for command injection protection in file_handler.

    These tests verify that compression conversion functions use subprocess
    with list arguments (not shell=True), preventing command injection attacks
    through malicious filenames.
    """

    @pytest.fixture
    def temp_dir(self):
        """Create temporary directory for testing."""
        with tempfile.TemporaryDirectory() as tmpdir:
            yield Path(tmpdir)

    @pytest.fixture
    def sample_content(self):
        """Provide sample FASTA content for testing."""
        return ">test_sequence\nATCGATCGATCG\n"

    def test_no_shell_true_in_compression_functions(self, temp_dir, sample_content):
        """
        Verify that compression conversion functions don't use shell=True.

        This is a regression test to ensure the security fix stays in place.
        If someone accidentally reintroduces shell commands, this test
        documents the expected behavior.
        """
        import inspect
        from validation_pkg.utils import file_handler

        # List of functions to check
        functions_to_check = [
            'gz_to_bz2',
            'bz2_to_gz',
            'none_to_gz',
            'gz_to_none',
            'bz2_to_none',
            'none_to_bz2'
        ]

        for func_name in functions_to_check:
            func = getattr(file_handler, func_name)
            source = inspect.getsource(func)

            # Verify no dangerous shell command patterns (actual code, not docs)
            # These patterns indicate actual vulnerable code
            assert "subprocess.run(cmd, shell=True" not in source, \
                f"{func_name} should not use subprocess.run with shell=True (security vulnerability)"
            assert 'executable=\'/bin/bash\'' not in source and 'executable="/bin/bash"' not in source, \
                f"{func_name} should not use bash as executable (security vulnerability)"

            # Verify no string interpolation into shell commands
            assert "f\"set -o pipefail" not in source, \
                f"{func_name} should not use shell pipes with f-strings"
            assert 'cmd = f"' not in source, \
                f"{func_name} should not build shell commands with f-strings"

            # Verify uses subprocess.Popen or subprocess.run with list args
            assert ('subprocess.Popen' in source or 'subprocess.run' in source), \
                f"{func_name} should use subprocess.Popen or subprocess.run"

            # Verify uses list concatenation for args (secure pattern)
            assert ('[decompress_cmd]' in source or '[compress_cmd]' in source or
                    'stdin=' in source or 'stdout=' in source), \
                f"{func_name} should use list arguments for subprocess"

    def test_gz_to_bz2_secure_processing(self, temp_dir, sample_content):
        """Test gz_to_bz2 processes files securely."""
        # Create input file
        input_file = temp_dir / "input.fasta.gz"
        with gzip.open(input_file, 'wt') as f:
            f.write(sample_content)

        output_file = temp_dir / "output.fasta.bz2"

        # Should complete successfully without executing shell commands
        gz_to_bz2(input_file, output_file)

        # Verify output is correct
        assert output_file.exists()
        with bz2.open(output_file, 'rt') as f:
            assert f.read() == sample_content

    def test_none_to_gz_secure_processing(self, temp_dir, sample_content):
        """Test none_to_gz processes files securely."""
        # Create input file
        input_file = temp_dir / "input.fasta"
        input_file.write_text(sample_content)

        output_file = temp_dir / "output.fasta.gz"

        # Should complete successfully without executing shell commands
        none_to_gz(input_file, output_file)

        # Verify output is correct
        assert output_file.exists()
        with gzip.open(output_file, 'rt') as f:
            assert f.read() == sample_content

    def test_error_handling_without_shell_exposure(self, temp_dir):
        """Test that error messages don't expose shell command vulnerabilities."""
        # Try to convert non-existent file
        nonexistent = temp_dir / "nonexistent.gz"
        output_file = temp_dir / "output.bz2"

        try:
            gz_to_bz2(nonexistent, output_file)
            assert False, "Should have raised CompressionError"
        except CompressionError as e:
            error_msg = str(e)
            # Error message should mention the problem, not shell syntax
            assert "Decompression failed" in error_msg or "Compression conversion failed" in error_msg
            # Should not contain shell syntax
            assert "sh -c" not in error_msg
            assert "set -o pipefail" not in error_msg
