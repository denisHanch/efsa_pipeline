"""
Utility functions for file handling with compression support.

Provides unified utilities for:
- File opening with automatic decompression
- Compression detection
- Format detection
- File validation
- ConfigManager parsing helpers
"""

import gzip
import bz2
from pathlib import Path
from typing import Union, TextIO, Type, Tuple, Any, Dict

from validation_pkg.utils.formats import CodingType, GenomeFormat, ReadFormat, FeatureFormat
from validation_pkg.exceptions import CompressionError

import subprocess
import shutil


# Cache for tool availability checks to avoid repeated subprocess calls
_TOOL_CACHE = {}
_LOGGER_INITIALIZED = False

# Compression tool configuration mapping
# Format: coding_type -> (parallel_tool, standard_tool, install_message)
_COMPRESSION_TOOLS = {
    CodingType.GZIP: {
        'parallel': 'pigz',
        'standard': 'gzip',
        'install_msg': 'sudo apt-get install pigz'
    },
    CodingType.BZIP2: {
        'parallel': 'pbzip2',
        'standard': 'bzip2',
        'install_msg': 'sudo apt-get install pbzip2'
    }
}

# Command arguments for different tools and modes
# Format: (tool_name, mode) -> args_generator_function
_COMMAND_ARGS = {
    ('pigz', 'compress'): lambda threads: ['-c', '-p', str(threads)],
    ('pigz', 'decompress'): lambda threads: ['-dc', '-p', str(threads)],
    ('gzip', 'compress'): lambda threads: ['-c'],
    ('gzip', 'decompress'): lambda threads: ['-dc'],
    ('pbzip2', 'compress'): lambda threads: ['-c', '-p' + str(threads)],
    ('pbzip2', 'decompress'): lambda threads: ['-dc', '-p' + str(threads)],
    ('bzip2', 'compress'): lambda threads: ['-c'],
    ('bzip2', 'decompress'): lambda threads: ['-dc'],
    ('cat', 'compress'): lambda threads: [],
    ('cat', 'decompress'): lambda threads: []
}


def check_compression_tool_available(tool_name: str) -> bool:
    """
    Check if a compression tool is available on the system.

    Results are cached to avoid repeated subprocess calls.

    Args:
        tool_name: Name of the tool to check ('pigz', 'pbzip2', 'gzip', 'bzip2')

    Returns:
        True if tool is available, False otherwise

    Example:
        >>> if check_compression_tool_available('pigz'):
        ...     print("pigz is available for faster compression")
    """
    if tool_name in _TOOL_CACHE:
        return _TOOL_CACHE[tool_name]

    # Check if tool exists using shutil.which
    available = shutil.which(tool_name) is not None
    _TOOL_CACHE[tool_name] = available

    return available


def _log_compression_tool(tool_name: str, threads: int, is_parallel: bool, install_msg: str = None):
    """
    Log compression tool usage (one-time initialization message).

    Args:
        tool_name: Name of the compression tool
        threads: Number of threads to use
        is_parallel: Whether this is a parallel tool
        install_msg: Installation instruction for parallel tool (if not available)
    """
    global _LOGGER_INITIALIZED

    if _LOGGER_INITIALIZED:
        return

    try:
        from validation_pkg.logger import get_logger
        logger = get_logger()

        if is_parallel:
            logger.info(f"Using {tool_name} for compression ({threads} threads)")
        else:
            logger.info(f"Using standard {tool_name} (install parallel tool for better performance: {install_msg})")

        _LOGGER_INITIALIZED = True
    except:
        pass


def _select_compression_tool(coding_type: CodingType) -> tuple:
    """
    Select the best available compression tool for the given coding type.

    Args:
        coding_type: CodingType enum (GZIP, BZIP2, or NONE)

    Returns:
        Tuple of (tool_name, is_parallel, install_msg)
    """
    # Handle no compression case
    if coding_type == CodingType.NONE:
        return ('cat', False, None)

    # Get tool configuration
    tool_config = _COMPRESSION_TOOLS.get(coding_type)
    if not tool_config:
        return ('cat', False, None)

    # Try parallel tool first
    parallel_tool = tool_config['parallel']
    if check_compression_tool_available(parallel_tool):
        return (parallel_tool, True, None)

    # Fallback to standard tool
    standard_tool = tool_config['standard']
    install_msg = tool_config['install_msg']
    return (standard_tool, False, install_msg)


def get_compression_command(coding_type: CodingType, mode: str = 'compress', threads: int = None) -> tuple:
    """
    Get the best available compression command for the given coding type.

    Prefers parallel tools (pigz, pbzip2) over standard tools (gzip, bzip2).
    Falls back gracefully if parallel tools are not available.

    Args:
        coding_type: CodingType enum (GZIP, BZIP2, or NONE)
        mode: 'compress' or 'decompress'
        threads: Number of threads to use (None = use default of 1)

    Returns:
        Tuple of (command_name, [args]) for subprocess

    Example:
        >>> cmd, args = get_compression_command(CodingType.GZIP, 'compress', threads=4)
        >>> # Returns ('pigz', ['-c', '-p', '4']) if pigz available
        >>> # Otherwise ('gzip', ['-c'])
    """
    # Default to 1 thread if not specified
    if threads is None:
        threads = 1

    # Select the best available tool
    tool_name, is_parallel, install_msg = _select_compression_tool(coding_type)

    # Log tool selection (one-time initialization)
    _log_compression_tool(tool_name, threads, is_parallel, install_msg)

    # Get command arguments from lookup table
    args_generator = _COMMAND_ARGS.get((tool_name, mode))
    if args_generator is None:
        # Fallback for unknown tool/mode combinations
        return ('cat', [])

    args = args_generator(threads)
    return (tool_name, args)


def open_file_with_coding_type(
    filepath: Union[str, Path],
    coding_type: CodingType,
    mode: str = 'rt'
) -> TextIO:
    """
    Open a file with automatic decompression based on CodingType enum.

    This function provides centralized file opening logic for all validators,
    eliminating code duplication and ensuring consistent compression handling.

    Args:
        filepath: Path to file
        coding_type: CodingType enum indicating compression
        mode: Opening mode (default: 'rt' for text read)

    Returns:
        File handle with automatic decompression

    Raises:
        CompressionError: If file cannot be opened or decompressed

    Example:
        >>> from validation_pkg.utils.formats import CodingType
        >>> from pathlib import Path
        >>> filepath = Path('genome.fasta.gz')
        >>> coding = CodingType.GZIP
        >>> with open_file_with_coding_type(filepath, coding) as f:
        ...     content = f.read()
    """
    filepath = Path(filepath)

    try:
        if coding_type == CodingType.GZIP:
            return gzip.open(filepath, mode)
        elif coding_type == CodingType.BZIP2:
            return bz2.open(filepath, mode)
        else:
            # No compression or unknown
            return open(filepath, mode)
    except Exception as e:
        raise CompressionError(f"Failed to open file {filepath}: {e}") from e


def detect_compression_type(filepath: Path) -> CodingType:
    """
    Detect compression type from file path and return CodingType enum.

    Checks the last extension to determine compression:
    - .gz or .gzip → CodingType.GZIP
    - .bz2 or .bzip2 → CodingType.BZIP2
    - other → CodingType.NONE

    Args:
        filepath: Path to file

    Returns:
        CodingType enum indicating compression

    Examples:
        >>> detect_compression_type(Path('genome.fasta'))
        CodingType.NONE
        >>> detect_compression_type(Path('genome.fasta.gz'))
        CodingType.GZIP
        >>> detect_compression_type(Path('genome.fasta.gzip'))
        CodingType.GZIP
        >>> detect_compression_type(Path('reads.fastq.bz2'))
        CodingType.BZIP2

    Note:
        The function will detect .gz from file.tar.gz, which will
        cause issues since validators do not handle TAR extraction.
        Users must extract TAR archives before processing.
    """
    suffixes = Path(filepath).suffixes

    if not suffixes:
        return CodingType.NONE

    # Check last extension for compression
    last_ext = suffixes[-1].lower()

    if last_ext in ['.gz', '.gzip']:
        return CodingType.GZIP
    elif last_ext in ['.bz2', '.bzip2']:
        return CodingType.BZIP2
    else:
        return CodingType.NONE


def detect_file_format(filepath: Path, format_enum: Type[Union[GenomeFormat, ReadFormat, FeatureFormat]]) -> Union[GenomeFormat, ReadFormat, FeatureFormat]:
    """
    Detect file format from file path and return appropriate enum.

    Checks the first extension (before compression) to determine format.

    Args:
        filepath: Path to file
        format_enum: Format enum class (GenomeFormat, ReadFormat, or FeatureFormat)

    Returns:
        Format enum member

    Raises:
        ValueError: If format cannot be determined

    Examples:
        >>> detect_file_format(Path('genome.fasta'), GenomeFormat)
        GenomeFormat.FASTA
        >>> detect_file_format(Path('genome.fasta.gz'), GenomeFormat)
        GenomeFormat.FASTA
        >>> detect_file_format(Path('reads.fastq.bz2'), ReadFormat)
        ReadFormat.FASTQ
        >>> detect_file_format(Path('features.gff'), FeatureFormat)
        FeatureFormat.GFF
    """
    suffixes = Path(filepath).suffixes

    if not suffixes:
        raise ValueError(f"Cannot determine format: no extension found in {filepath.name}")

    # Get first extension (the format, before compression)
    format_ext = suffixes[0]

    # Let the format enum's _missing_ method handle the conversion
    return format_enum(format_ext)


def parse_config_file_value(
    value: Any,
    field_name: str
) -> Tuple[str, Dict[str, Any]]:
    """
    Parse file configuration value from JSON config.

    Handles both dict and string formats:
    - Dict: {"filename": "genome.fasta", "extra_key": "value"}
    - String: "genome.fasta" (backwards compatibility)

    Args:
        value: Configuration value (dict or string)
        field_name: Name of configuration field (for error messages)

    Returns:
        Tuple of (filename, extra_fields_dict)

    Raises:
        ValueError: If value is invalid format

    Examples:
        >>> parse_config_file_value({"filename": "genome.fasta"}, "ref_genome")
        ("genome.fasta", {})
        >>> parse_config_file_value({"filename": "reads.fastq", "ngs_type": "illumina"}, "reads")
        ("reads.fastq", {"ngs_type": "illumina"})
        >>> parse_config_file_value("genome.fasta", "ref_genome")
        ("genome.fasta", {})
    """
    filename = None
    extra = {}

    if isinstance(value, dict):
        if 'filename' not in value:
            raise ValueError(f"{field_name} must contain 'filename' field")
        filename = value['filename']

        # Extract extra keys (excluding 'filename')
        extra = {k: v for k, v in value.items() if k != 'filename'}

    elif isinstance(value, str):
        # Accept plain string for backwards compatibility
        filename = value
    else:
        raise ValueError(f"{field_name} must be a dict or string, got {type(value).__name__}")

    return filename, extra


def gz_to_bz2(gz_file: Path, bz2_file: Path, threads: int = None):
    """
    Convert gzip compressed file to bzip2.

    Uses parallel tools (pigz/pbzip2) if available, falls back to gzip/bzip2.
    Securely pipes decompression to compression using subprocess.Popen without shell.

    Args:
        gz_file: Path to input .gz file
        bz2_file: Path to output .bz2 file
        threads: Number of threads to use (None = auto-detect)

    Raises:
        CompressionError: If conversion fails

    Security:
        Uses subprocess.Popen with list arguments (no shell=True) to prevent
        command injection attacks through malicious filenames.
    """
    # Get best available compression commands
    decompress_cmd, decompress_args = get_compression_command(CodingType.GZIP, 'decompress', threads)
    compress_cmd, compress_args = get_compression_command(CodingType.BZIP2, 'compress', threads)

    try:
        # Create decompression process (reads input file, writes to stdout)
        decompress_proc = subprocess.Popen(
            [decompress_cmd] + decompress_args + [str(gz_file)],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        )

        # Create compression process (reads from stdin, writes to output file)
        with open(bz2_file, 'wb') as output_file:
            compress_proc = subprocess.Popen(
                [compress_cmd] + compress_args,
                stdin=decompress_proc.stdout,
                stdout=output_file,
                stderr=subprocess.PIPE
            )

            # Close decompress stdout in parent to allow SIGPIPE propagation
            decompress_proc.stdout.close()

            # Wait for both processes to complete
            compress_stdout, compress_stderr = compress_proc.communicate()
            decompress_returncode = decompress_proc.wait()

            # Check for errors
            if decompress_returncode != 0:
                _, decompress_stderr = decompress_proc.communicate()
                raise CompressionError(
                    f"Decompression failed (gzip): {decompress_stderr.decode('utf-8', errors='replace')}"
                )

            if compress_proc.returncode != 0:
                raise CompressionError(
                    f"Compression failed (bzip2): {compress_stderr.decode('utf-8', errors='replace')}"
                )

    except (OSError, FileNotFoundError) as e:
        raise CompressionError(f"Compression tool not found: {e}") from e
    except Exception as e:
        raise CompressionError(f"Compression conversion failed: {e}") from e


def bz2_to_gz(bz2_file: Path, gz_file: Path, threads: int = None):
    """
    Convert bzip2 compressed file to gzip.

    Uses parallel tools (pbzip2/pigz) if available, falls back to bzip2/gzip.
    Securely pipes decompression to compression using subprocess.Popen without shell.

    Args:
        bz2_file: Path to input .bz2 file
        gz_file: Path to output .gz file
        threads: Number of threads to use (None = auto-detect)

    Raises:
        CompressionError: If conversion fails

    Security:
        Uses subprocess.Popen with list arguments (no shell=True) to prevent
        command injection attacks through malicious filenames.
    """
    # Get best available compression commands
    decompress_cmd, decompress_args = get_compression_command(CodingType.BZIP2, 'decompress', threads)
    compress_cmd, compress_args = get_compression_command(CodingType.GZIP, 'compress', threads)

    try:
        # Create decompression process (reads input file, writes to stdout)
        decompress_proc = subprocess.Popen(
            [decompress_cmd] + decompress_args + [str(bz2_file)],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        )

        # Create compression process (reads from stdin, writes to output file)
        with open(gz_file, 'wb') as output_file:
            compress_proc = subprocess.Popen(
                [compress_cmd] + compress_args,
                stdin=decompress_proc.stdout,
                stdout=output_file,
                stderr=subprocess.PIPE
            )

            # Close decompress stdout in parent to allow SIGPIPE propagation
            decompress_proc.stdout.close()

            # Wait for both processes to complete
            compress_stdout, compress_stderr = compress_proc.communicate()
            decompress_returncode = decompress_proc.wait()

            # Check for errors
            if decompress_returncode != 0:
                _, decompress_stderr = decompress_proc.communicate()
                raise CompressionError(
                    f"Decompression failed (bzip2): {decompress_stderr.decode('utf-8', errors='replace')}"
                )

            if compress_proc.returncode != 0:
                raise CompressionError(
                    f"Compression failed (gzip): {compress_stderr.decode('utf-8', errors='replace')}"
                )

    except (OSError, FileNotFoundError) as e:
        raise CompressionError(f"Compression tool not found: {e}") from e
    except Exception as e:
        raise CompressionError(f"Compression conversion failed: {e}") from e


def none_to_gz(none_file: Path, gz_file: Path, threads: int = None):
    """
    Compress uncompressed file to gzip.

    Uses pigz if available for parallel compression, falls back to gzip.

    Args:
        none_file: Path to input uncompressed file
        gz_file: Path to output .gz file
        threads: Number of threads to use (None = auto-detect)

    Raises:
        CompressionError: If compression fails

    Security:
        Uses subprocess with list arguments (no shell=True) to prevent
        command injection attacks through malicious filenames.
    """
    compress_cmd, compress_args = get_compression_command(CodingType.GZIP, 'compress', threads)

    try:
        with open(none_file, 'rb') as input_file, open(gz_file, 'wb') as output_file:
            subprocess.run(
                [compress_cmd] + compress_args,
                stdin=input_file,
                stdout=output_file,
                stderr=subprocess.PIPE,
                check=True
            )
    except subprocess.CalledProcessError as e:
        raise CompressionError(
            f"Compression failed (gzip): {e.stderr.decode('utf-8', errors='replace')}"
        ) from e
    except (OSError, FileNotFoundError) as e:
        raise CompressionError(f"Compression tool not found: {e}") from e


def gz_to_none(gz_file: Path, none_file: Path, threads: int = None):
    """
    Decompress gzip file to uncompressed file.

    Uses pigz if available for parallel decompression, falls back to gzip.

    Args:
        gz_file: Path to input .gz file
        none_file: Path to output uncompressed file
        threads: Number of threads to use (None = auto-detect)

    Raises:
        CompressionError: If decompression fails

    Security:
        Uses subprocess with list arguments (no shell=True) to prevent
        command injection attacks through malicious filenames.
    """
    decompress_cmd, decompress_args = get_compression_command(CodingType.GZIP, 'decompress', threads)

    try:
        with open(none_file, 'wb') as output_file:
            subprocess.run(
                [decompress_cmd] + decompress_args + [str(gz_file)],
                stdout=output_file,
                stderr=subprocess.PIPE,
                check=True
            )
    except subprocess.CalledProcessError as e:
        raise CompressionError(
            f"Decompression failed (gzip): {e.stderr.decode('utf-8', errors='replace')}"
        ) from e
    except (OSError, FileNotFoundError) as e:
        raise CompressionError(f"Compression tool not found: {e}") from e


def bz2_to_none(bz2_file: Path, none_file: Path, threads: int = None):
    """
    Decompress bzip2 file to uncompressed file.

    Uses pbzip2 if available for parallel decompression, falls back to bzip2.

    Args:
        bz2_file: Path to input .bz2 file
        none_file: Path to output uncompressed file
        threads: Number of threads to use (None = auto-detect)

    Raises:
        CompressionError: If decompression fails

    Security:
        Uses subprocess with list arguments (no shell=True) to prevent
        command injection attacks through malicious filenames.
    """
    decompress_cmd, decompress_args = get_compression_command(CodingType.BZIP2, 'decompress', threads)

    try:
        with open(none_file, 'wb') as output_file:
            subprocess.run(
                [decompress_cmd] + decompress_args + [str(bz2_file)],
                stdout=output_file,
                stderr=subprocess.PIPE,
                check=True
            )
    except subprocess.CalledProcessError as e:
        raise CompressionError(
            f"Decompression failed (bzip2): {e.stderr.decode('utf-8', errors='replace')}"
        ) from e
    except (OSError, FileNotFoundError) as e:
        raise CompressionError(f"Compression tool not found: {e}") from e


def none_to_bz2(none_file: Path, bz2_file: Path, threads: int = None):
    """
    Compress uncompressed file to bzip2.

    Uses pbzip2 if available for parallel compression, falls back to bzip2.

    Args:
        none_file: Path to input uncompressed file
        bz2_file: Path to output .bz2 file
        threads: Number of threads to use (None = auto-detect)

    Raises:
        CompressionError: If compression fails

    Security:
        Uses subprocess with list arguments (no shell=True) to prevent
        command injection attacks through malicious filenames.
    """
    compress_cmd, compress_args = get_compression_command(CodingType.BZIP2, 'compress', threads)

    try:
        with open(none_file, 'rb') as input_file, open(bz2_file, 'wb') as output_file:
            subprocess.run(
                [compress_cmd] + compress_args,
                stdin=input_file,
                stdout=output_file,
                stderr=subprocess.PIPE,
                check=True
            )
    except subprocess.CalledProcessError as e:
        raise CompressionError(
            f"Compression failed (bzip2): {e.stderr.decode('utf-8', errors='replace')}"
        ) from e
    except (OSError, FileNotFoundError) as e:
        raise CompressionError(f"Compression tool not found: {e}") from e


def open_compressed_writer(filepath: Union[str, Path], coding_type: CodingType, use_parallel: bool = True, threads: int = None):
    """
    Open a file handle for writing compressed data efficiently.

    Prefers parallel compression tools (pigz/pbzip2) when use_parallel=True.
    Falls back to Python gzip/bz2 libraries for compatibility.

    This function returns a context manager that should be used with 'with' statement.

    Args:
        filepath: Path to output file
        coding_type: CodingType enum (GZIP, BZIP2, or NONE)
        use_parallel: If True, use parallel tools (pigz/pbzip2) when available
        threads: Number of threads to use (None = auto-detect)

    Returns:
        File handle for writing (text mode)

    Example:
        >>> from validation_pkg.utils.formats import CodingType
        >>> from pathlib import Path
        >>> with open_compressed_writer('output.gz', CodingType.GZIP, threads=4) as f:
        ...     f.write("Hello, world!\\n")

    Note:
        For parallel tools to work, pigz or pbzip2 must be installed.
        If not available or use_parallel=False, falls back to Python libraries.
    """
    filepath = Path(filepath)

    # If parallel compression is requested and available, use subprocess approach
    if use_parallel and coding_type in (CodingType.GZIP, CodingType.BZIP2):
        # Check if parallel tools are available
        use_subprocess = False

        if coding_type == CodingType.GZIP:
            use_subprocess = check_compression_tool_available('pigz')
        elif coding_type == CodingType.BZIP2:
            use_subprocess = check_compression_tool_available('pbzip2')

        if use_subprocess:
            # Use subprocess pipe for parallel compression
            compress_cmd, compress_args = get_compression_command(coding_type, 'compress', threads)

            # Create subprocess with pipe to stdin
            process = subprocess.Popen(
                [compress_cmd] + compress_args,
                stdin=subprocess.PIPE,
                stdout=open(filepath, 'wb'),
                stderr=subprocess.PIPE,
                text=False  # Binary mode for pipe
            )

            # Wrap in a context manager that handles text encoding and cleanup
            class SubprocessWriter:
                def __init__(self, proc):
                    self.proc = proc
                    self.stdin = proc.stdin

                def write(self, text: str):
                    """Write text to subprocess stdin."""
                    self.stdin.write(text.encode('utf-8'))

                def __enter__(self):
                    return self

                def __exit__(self, exc_type, exc_val, exc_tb):
                    self.stdin.close()
                    self.proc.wait()
                    if self.proc.returncode != 0:
                        stderr_output = self.proc.stderr.read().decode('utf-8')
                        raise CompressionError(
                            f"Compression failed with return code {self.proc.returncode}: {stderr_output}"
                        )
                    return False

            return SubprocessWriter(process)

    # Fallback to Python compression libraries
    if coding_type == CodingType.GZIP:
        return gzip.open(filepath, 'wt')
    elif coding_type == CodingType.BZIP2:
        return bz2.open(filepath, 'wt')
    else:
        # No compression
        return open(filepath, 'w')