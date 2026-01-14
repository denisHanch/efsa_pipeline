"""
Tests for the file_handler module.

Tests cover:
- File opening with automatic decompression
- Compression detection and type handling
- Format detection for genome, read, and feature files
- Filename processing and path utilities
- Compression conversion functions (gzip, bzip2)
- Parallel compression tool selection and usage
"""

import pytest
import gzip
import bz2
import tempfile
import subprocess
from pathlib import Path

from validation_pkg.utils.file_handler import *
from validation_pkg.utils.formats import CodingType,ReadFormat,FeatureFormat,GenomeFormat
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

    def test_check_tool_available_gzip(self):
        """Test that gzip is available (should always be available on Unix)."""
        assert check_tool_available('gzip') is True

    def test_check_tool_available_bzip2(self):
        """Test that bzip2 is available (should always be available on Unix)."""
        assert check_tool_available('bzip2') is True

    def test_check_tool_available_nonexistent(self):
        """Test that nonexistent tools return False."""
        assert check_tool_available('nonexistent_compression_tool_xyz') is False

    def test_check_tool_available_caching(self):
        """Test that tool availability checks are cached."""
        # First call
        result1 = check_tool_available('gzip')
        # Second call should use cache
        result2 = check_tool_available('gzip')
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


class TestDetectCompressionType:
    """Test compression type detection from file paths."""

    def test_detect_compression_type_none(self):
        """Test detecting no compression."""
        from validation_pkg.utils.file_handler import detect_compression_type
        filepath = Path("genome.fasta")
        coding_type = detect_compression_type(filepath)
        assert coding_type == CodingType.NONE

    def test_detect_compression_type_gzip(self):
        """Test detecting gzip compression."""
        from validation_pkg.utils.file_handler import detect_compression_type
        filepath = Path("genome.fasta.gz")
        coding_type = detect_compression_type(filepath)
        assert coding_type == CodingType.GZIP

    def test_detect_compression_type_gzip_alt(self):
        """Test detecting gzip with .gzip extension."""
        from validation_pkg.utils.file_handler import detect_compression_type
        filepath = Path("genome.fasta.gzip")
        coding_type = detect_compression_type(filepath)
        assert coding_type == CodingType.GZIP

    def test_detect_compression_type_bzip2(self):
        """Test detecting bzip2 compression."""
        from validation_pkg.utils.file_handler import detect_compression_type
        filepath = Path("genome.fasta.bz2")
        coding_type = detect_compression_type(filepath)
        assert coding_type == CodingType.BZIP2

    def test_detect_compression_type_bzip2_alt(self):
        """Test detecting bzip2 with .bzip2 extension."""
        from validation_pkg.utils.file_handler import detect_compression_type
        filepath = Path("genome.fasta.bzip2")
        coding_type = detect_compression_type(filepath)
        assert coding_type == CodingType.BZIP2

    def test_detect_compression_case_insensitive(self):
        """Test that compression detection is case insensitive."""
        from validation_pkg.utils.file_handler import detect_compression_type
        filepath = Path("genome.fasta.GZ")
        coding_type = detect_compression_type(filepath)
        assert coding_type == CodingType.GZIP


class TestPairedEndFormatDetection:
    """Test format detection with paired-end indicators."""
    
    @pytest.fixture
    def temp_dir(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            yield Path(tmpdir)
    
    def test_standard_naming(self, temp_dir):
        """Test format detection for standard filenames."""
        # Standard 2-extension files
        test_cases = [
            ('sample.fastq.gz', CodingType.GZIP, ReadFormat.FASTQ),
            ('sample.fq.bz2', CodingType.BZIP2, ReadFormat.FASTQ),
            ('sample.bam', CodingType.NONE, ReadFormat.BAM),
            ('genome.fasta.gz', CodingType.GZIP, GenomeFormat.FASTA),
            ('features.gff.gz', CodingType.GZIP, FeatureFormat.GFF),
        ]

        for filename, expected_coding, expected_format in test_cases:
            filepath = temp_dir / filename
            filepath.touch()

            coding = detect_compression_type(filepath)
            # Determine correct format enum based on expected format type
            if expected_format in [ReadFormat.FASTQ, ReadFormat.BAM]:
                format_enum = ReadFormat
            elif expected_format in [GenomeFormat.FASTA, GenomeFormat.GENBANK]:
                format_enum = GenomeFormat
            else:  # FeatureFormat
                format_enum = FeatureFormat
            fmt = detect_file_format(filepath, format_enum)

            assert coding == expected_coding, f"Wrong compression for {filename}"
            assert fmt == expected_format, f"Wrong format for {filename}"
    
    def test_paired_end_dot_separator(self, temp_dir):
        """Test format detection for dot-separated paired-end files (.R1, .R2, .1, .2)."""
        test_cases = [
            ('sample.R1.fastq.gz', CodingType.GZIP, ReadFormat.FASTQ),
            ('sample.R2.fastq.gz', CodingType.GZIP, ReadFormat.FASTQ),
            ('sample.1.fq.bz2', CodingType.BZIP2, ReadFormat.FASTQ),
            ('sample.2.fq.bz2', CodingType.BZIP2, ReadFormat.FASTQ),
            ('sample.R1.fastq', CodingType.NONE, ReadFormat.FASTQ),
            ('sample.R2.fq', CodingType.NONE, ReadFormat.FASTQ),
        ]
        
        for filename, expected_coding, expected_format in test_cases:
            filepath = temp_dir / filename
            filepath.touch()
            
            coding = detect_compression_type(filepath)
            fmt = detect_file_format(filepath, ReadFormat)
            
            assert coding == expected_coding, f"Wrong compression for {filename}"
            assert fmt == expected_format, f"Wrong format for {filename}"
    
    def test_paired_end_underscore_separator(self, temp_dir):
        """Test format detection for underscore-separated paired-end files (_R1, _R2, _1, _2)."""
        # These should work without any changes (underscore is not a dot separator)
        test_cases = [
            ('sample_R1.fastq.gz', CodingType.GZIP, ReadFormat.FASTQ),
            ('sample_R2.fastq.gz', CodingType.GZIP, ReadFormat.FASTQ),
            ('sample_1.fq.bz2', CodingType.BZIP2, ReadFormat.FASTQ),
            ('sample_2.fq.bz2', CodingType.BZIP2, ReadFormat.FASTQ),
        ]
        
        for filename, expected_coding, expected_format in test_cases:
            filepath = temp_dir / filename
            filepath.touch()
            
            coding = detect_compression_type(filepath)
            fmt = detect_file_format(filepath, ReadFormat)
            
            assert coding == expected_coding, f"Wrong compression for {filename}"
            assert fmt == expected_format, f"Wrong format for {filename}"
    
    def test_arbitrary_prefixes(self, temp_dir):
        """Test that arbitrary prefix extensions are ignored."""
        test_cases = [
            ('sample.processed.fastq.gz', CodingType.GZIP, ReadFormat.FASTQ),
            ('sample.v2.final.fastq.gz', CodingType.GZIP, ReadFormat.FASTQ),
            ('sample.filtered.R1.fastq.gz', CodingType.GZIP, ReadFormat.FASTQ),
        ]
        
        for filename, expected_coding, expected_format in test_cases:
            filepath = temp_dir / filename
            filepath.touch()
            
            coding = detect_compression_type(filepath)
            fmt = detect_file_format(filepath, ReadFormat)
            
            assert coding == expected_coding, f"Wrong compression for {filename}"
            assert fmt == expected_format, f"Wrong format for {filename}"
    
    def test_user_provided_examples(self, temp_dir):
        """Test all user-provided filename examples."""
        test_cases = [
            'sample_S1_L001_R1_001.fastq.gz'
            'Bacillus_toyonensis_LE1_S1_L001_R1.fastq.gz',
            'Bacillus_toyonensis_LE1_S1_L001_R2.fastq.gz',
            'Bacillus-toyonensis-LE1-sequence_1.fastq.gz',
            'Bacillus-toyonensis-LE1-sequence_2.fastq.gz',
            'Bacillus_toyonensis_LE1_sequence.R1.fastq.gz',  # Dot separator
            'SRR834393.fq.gz',
            'SRR834394_1.fq.gz',
            'SRR837394_1.fastq.gz',
            'SRR837394_1_001.fastq.gz',
            'SRR837394_2_001.fastq.gz',
            'SRR834393R1_combined.fq.gz',
            'SRR834393R2_combined.fq.gz',
        ]
        
        for filename in test_cases:
            filepath = temp_dir / filename
            filepath.touch()
            
            # Should not raise any exceptions
            coding = detect_compression_type(filepath)
            fmt = detect_file_format(filepath, ReadFormat)
            
            # All should be FASTQ with GZIP compression
            assert coding == CodingType.GZIP, f"Wrong compression for {filename}"
            assert fmt == ReadFormat.FASTQ, f"Wrong format for {filename}"


class TestDetectFileFormat:
    """Test file format detection from file paths."""

    def test_detect_file_format_fasta(self):
        """Test detecting FASTA format."""
        from validation_pkg.utils.file_handler import detect_file_format
        from validation_pkg.utils.formats import GenomeFormat
        filepath = Path("genome.fasta")
        file_format = detect_file_format(filepath, GenomeFormat)
        assert file_format == GenomeFormat.FASTA

    def test_detect_file_format_fasta_compressed(self):
        """Test detecting FASTA format from compressed file."""
        from validation_pkg.utils.file_handler import detect_file_format
        from validation_pkg.utils.formats import GenomeFormat
        filepath = Path("genome.fasta.gz")
        file_format = detect_file_format(filepath, GenomeFormat)
        assert file_format == GenomeFormat.FASTA

    def test_detect_file_format_genbank(self):
        """Test detecting GenBank format."""
        from validation_pkg.utils.file_handler import detect_file_format
        from validation_pkg.utils.formats import GenomeFormat
        filepath = Path("genome.gbk")
        file_format = detect_file_format(filepath, GenomeFormat)
        assert file_format == GenomeFormat.GENBANK

    def test_detect_file_format_fastq(self):
        """Test detecting FASTQ format."""
        from validation_pkg.utils.file_handler import detect_file_format
        from validation_pkg.utils.formats import ReadFormat
        filepath = Path("reads.fastq.gz")
        file_format = detect_file_format(filepath, ReadFormat)
        assert file_format == ReadFormat.FASTQ

    def test_detect_file_format_gff(self):
        """Test detecting GFF format."""
        from validation_pkg.utils.file_handler import detect_file_format
        from validation_pkg.utils.formats import FeatureFormat
        filepath = Path("features.gff")
        file_format = detect_file_format(filepath, FeatureFormat)
        assert file_format == FeatureFormat.GFF

    def test_detect_file_format_no_extension(self):
        """Test error when no extension found."""
        from validation_pkg.utils.file_handler import detect_file_format
        from validation_pkg.utils.formats import GenomeFormat
        filepath = Path("genome")
        with pytest.raises(ValueError, match="no extension found"):
            detect_file_format(filepath, GenomeFormat)


class TestParseConfigFileValue:
    """Test parsing of configuration file values."""

    def test_parse_config_file_value_dict(self):
        """Test parsing config value as dict."""
        from validation_pkg.utils.file_handler import parse_config_file_value
        value = {"filename": "genome.fasta", "extra_key": "value"}
        filename, extra = parse_config_file_value(value, "ref_genome")
        assert filename == "genome.fasta"
        assert extra == {"extra_key": "value"}

    def test_parse_config_file_value_string(self):
        """Test parsing config value as string."""
        from validation_pkg.utils.file_handler import parse_config_file_value
        value = "genome.fasta"
        filename, extra = parse_config_file_value(value, "ref_genome")
        assert filename == "genome.fasta"
        assert extra == {}

    def test_parse_config_file_value_missing_filename(self):
        """Test error when dict missing filename field."""
        from validation_pkg.utils.file_handler import parse_config_file_value
        value = {"other_key": "value"}
        with pytest.raises(ValueError, match="must contain 'filename' field"):
            parse_config_file_value(value, "ref_genome")

    def test_parse_config_file_value_invalid_type(self):
        """Test error with invalid value type."""
        from validation_pkg.utils.file_handler import parse_config_file_value
        value = 123  # Not a dict or string
        with pytest.raises(ValueError, match="must be a dict or string"):
            parse_config_file_value(value, "ref_genome")

    def test_parse_config_file_value_dict_excludes_filename(self):
        """Test that filename is excluded from extra dict."""
        from validation_pkg.utils.file_handler import parse_config_file_value
        value = {"filename": "genome.fasta", "key1": "val1", "key2": "val2"}
        filename, extra = parse_config_file_value(value, "ref_genome")
        assert filename == "genome.fasta"
        assert "filename" not in extra
        assert extra == {"key1": "val1", "key2": "val2"}

