"""Utility functions for file handling with compression support."""

import gzip
import bz2
import re
import threading
from pathlib import Path
from typing import Union, TextIO, Type, Tuple, Any, Dict, Optional

from validation_pkg.utils.formats import CodingType, GenomeFormat, ReadFormat, FeatureFormat
from validation_pkg.exceptions import CompressionError

import subprocess
import shutil

__all__ = [
    # Tool availability and compression commands
    'check_tool_available',
    'get_compression_command',

    # File operations
    'open_file_with_coding_type',
    'open_compressed_writer',

    # Detection functions
    'detect_compression_type',
    'detect_file_format',

    # ConfigManager helpers
    'parse_config_file_value',

    # Compression conversion functions
    'gz_to_bz2',
    'bz2_to_gz',
    'none_to_gz',
    'gz_to_none',
    'bz2_to_none',
    'none_to_bz2',

    # Path utilities
    'get_incremented_path',
    'strip_all_extensions',
    'build_output_filename',
    'build_output_path',
]


# Cache for tool availability checks to avoid repeated subprocess calls
# Protected by _TOOL_CACHE_LOCK for thread-safe access
_TOOL_CACHE = {}
_TOOL_CACHE_LOCK = threading.Lock()
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


def check_tool_available(tool_name: str) -> bool:
    """Check if a tool is available on the system."""
    # Fast path: check cache without lock (safe for reads)
    if tool_name in _TOOL_CACHE:
        return _TOOL_CACHE[tool_name]

    # Slow path: check tool and update cache with lock
    with _TOOL_CACHE_LOCK:
        # Double-check inside lock (another thread might have updated it)
        if tool_name in _TOOL_CACHE:
            return _TOOL_CACHE[tool_name]

        # Check if tool exists using shutil.which
        available = shutil.which(tool_name) is not None
        _TOOL_CACHE[tool_name] = available

        return available


def _log_compression_tool(tool_name: str, threads: int, is_parallel: bool, install_msg: str = None):
    """Log compression tool usage (one-time initialization message)."""
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
    except Exception:
        # If logger initialization fails, continue without logging
        # (This allows the module to work even if logger is unavailable)
        pass


def _select_compression_tool(coding_type: CodingType) -> tuple:
    """Select the best available compression tool for the given coding type."""
    # Handle no compression case
    if coding_type == CodingType.NONE:
        return ('cat', False, None)

    # Get tool configuration
    tool_config = _COMPRESSION_TOOLS.get(coding_type)
    if not tool_config:
        return ('cat', False, None)

    # Try parallel tool first
    parallel_tool = tool_config['parallel']
    if check_tool_available(parallel_tool):
        return (parallel_tool, True, None)

    # Fallback to standard tool
    standard_tool = tool_config['standard']
    install_msg = tool_config['install_msg']
    return (standard_tool, False, install_msg)


def get_compression_command(coding_type: CodingType, mode: str = 'compress', threads: int = None) -> tuple:
    """Get the best available compression command for the given coding type."""
    # Default to 1 thread if not specified
    if threads is None:
        threads = 1

    # For single thread, always use standard tools (better performance)
    # Parallel tools have overhead that makes them slower with threads=1
    if threads == 1:
        if coding_type == CodingType.GZIP:
            tool_name = 'gzip'
            is_parallel = False
        elif coding_type == CodingType.BZIP2:
            tool_name = 'bzip2'
            is_parallel = False
        else:
            tool_name = 'cat'
            is_parallel = False
        install_msg = None
    else:
        # Select the best available tool for multi-threading
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
    """Open a file with automatic decompression based on CodingType enum."""
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
    """Detect compression type from file path and return CodingType enum."""
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
    """Detect file format from extension."""
    suffixes = Path(filepath).suffixes

    if not suffixes:
        raise ValueError(f"Cannot determine format: no extension found in {filepath.name}")

    # Check if last extension is a known compression type
    compression_extensions = {'.gz', '.gzip', '.bz2', '.bzip2'}
    last_ext_lower = suffixes[-1].lower()
    has_compression = last_ext_lower in compression_extensions

    # Determine format extension based on compression status
    if len(suffixes) == 1:
        # Only one extension: must be the format (no compression)
        # Example: sample.fastq
        format_ext = suffixes[0]
    elif has_compression:
        # Last extension is compression: format is second-to-last
        # Example: sample.fastq.gz → suffixes[-2] = .fastq
        # Example: sample.R1.fastq.gz → suffixes[-2] = .fastq (3+ extensions)
        format_ext = suffixes[-2]
    else:
        # Multiple extensions but last is NOT compression: format is last
        # Example: sample.R1.fastq → suffixes[-1] = .fastq
        # Example: sample.processed.fasta → suffixes[-1] = .fasta
        format_ext = suffixes[-1]

    # Let the format enum's _missing_ method handle the conversion
    try:
        return format_enum(format_ext)
    except ValueError as e:
        raise ValueError(
            f"Cannot determine format for {filepath.name}: {e}\n"
            f"Supported formats: {', '.join([fmt.name for fmt in format_enum])}"
        )


def parse_config_file_value(
    value: Any,
    field_name: str
) -> Tuple[str, Dict[str, Any]]:
    """Parse file configuration value from JSON config."""
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
    """Convert gzip compressed file to bzip2."""
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
    """Convert bzip2 compressed file to gzip."""
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
    """Compress uncompressed file to gzip."""
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
    """Decompress gzip file to uncompressed file."""
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
    """Decompress bzip2 file to uncompressed file."""
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
    """Compress uncompressed file to bzip2."""
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
    """Open a file handle for writing compressed data efficiently."""
    filepath = Path(filepath)

    # If parallel compression is requested and available, use subprocess approach
    if use_parallel and coding_type in (CodingType.GZIP, CodingType.BZIP2):
        # Check if parallel tools are available
        use_subprocess = False

        if coding_type == CodingType.GZIP:
            use_subprocess = check_tool_available('pigz')
        elif coding_type == CodingType.BZIP2:
            use_subprocess = check_tool_available('pbzip2')

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


def get_incremented_path(path: Path, separator: str = "_") -> Path:
    """Get next available filename by auto-incrementing if file exists."""
    # Ensure path is a Path object
    path = Path(path)

    # If file doesn't exist, return original path
    if not path.exists():
        return path

    stem = path.stem
    suffix = path.suffix
    parent = path.parent

    # Check if stem already has increment pattern (e.g., report_001)
    match = re.match(r'^(.+)_(\d+)$', stem)
    if match:
        base_stem = match.group(1)
        start_counter = int(match.group(2)) + 1
    else:
        base_stem = stem
        start_counter = 1

    # Find next available number
    counter = start_counter
    while True:
        new_name = f"{base_stem}{separator}{counter:03d}{suffix}"
        new_path = parent / new_name
        if not new_path.exists():
            return new_path
        counter += 1

        # Safety check to avoid infinite loop
        if counter > 9999:
            raise RuntimeError(f"Too many incremented files for {path}. Maximum is 9999.")


def strip_all_extensions(filename: str, path: Optional[Path] = None) -> str:
    """
    Strip all extensions from a filename.

    Args:
        filename: The filename to process
        path: Optional Path object to get suffixes from

    Returns:
        Filename with all extensions removed

    Examples:
        >>> strip_all_extensions("genome.fasta.gz")
        'genome'
        >>> strip_all_extensions("reads_R1.fastq")
        'reads_R1'
    """
    if not filename:
        return filename

    # If path provided, use its suffixes list (more reliable)
    if path:
        base_name = filename
        for suffix in path.suffixes:
            base_name = base_name.replace(suffix, '', 1)
        return base_name

    # Otherwise, iteratively remove extensions
    result = filename
    while True:
        name_without_ext = Path(result).stem
        if name_without_ext == result:
            break
        result = name_without_ext

    return result


def build_output_filename(
    input_filename: str,
    output_format: str,
    coding_type: Optional[CodingType] = None,
    suffix: Optional[str] = None,
    input_path: Optional[Path] = None
) -> str:
    """
    Build output filename with consistent naming convention.

    Args:
        input_filename: Original input filename
        output_format: Output format extension (e.g., 'fasta', 'fastq', 'gff')
        coding_type: Optional compression type (adds .gz, .bz2, etc.)
        suffix: Optional custom suffix to add before extension
        input_path: Optional input Path object for accurate extension stripping

    Returns:
        Output filename with proper extensions

    Examples:
        >>> build_output_filename("genome.fa.gz", "fasta", CodingType.GZIP, "validated")
        'genome_validated.fasta.gz'
        >>> build_output_filename("reads.fq", "fastq", CodingType.NONE)
        'reads.fastq'
    """
    # Strip all extensions from input filename
    base_name = strip_all_extensions(input_filename, input_path)

    # Add custom suffix if provided
    if suffix:
        filename = f"{base_name}_{suffix}.{output_format}"
    else:
        filename = f"{base_name}.{output_format}"

    # Add compression extension if specified
    if coding_type:
        filename += coding_type.to_extension()

    return filename


def build_output_path(
    base_dir: Path,
    input_filename: str,
    output_format: str,
    coding_type: Optional[CodingType] = None,
    subdir_name: Optional[str] = None,
    filename_suffix: Optional[str] = None,
    input_path: Optional[Path] = None,
    create_dirs: bool = True
) -> Path:
    """
    Build complete output path with consistent naming and directory structure.

    This is the main entry point for building output paths. It handles:
    - Safe subdirectory creation (with path traversal protection)
    - Extension stripping
    - Custom suffix addition
    - Compression extension handling

    Args:
        base_dir: Base output directory
        input_filename: Original input filename
        output_format: Output format extension (e.g., 'fasta', 'fastq', 'gff')
        coding_type: Optional compression type
        subdir_name: Optional subdirectory name (will be sanitized)
        filename_suffix: Optional suffix to add to filename
        input_path: Optional input Path for accurate extension handling
        create_dirs: Whether to create directories (default: True)

    Returns:
        Complete output path

    Raises:
        ValueError: If subdir_name contains illegal characters
        ConfigurationError: If path would escape base_dir

    Examples:
        >>> build_output_path(
        ...     Path("/output"),
        ...     "genome.fasta.gz",
        ...     "fasta",
        ...     CodingType.GZIP,
        ...     subdir_name="validated",
        ...     filename_suffix="filtered"
        ... )
        Path("/output/validated/genome_filtered.fasta.gz")
    """
    from validation_pkg.utils.path_utils import build_safe_output_dir

    # Build safe output directory (with path traversal protection)
    output_dir = build_safe_output_dir(base_dir, subdir_name, create=create_dirs)

    # Build output filename
    output_filename = build_output_filename(
        input_filename,
        output_format,
        coding_type,
        filename_suffix,
        input_path
    )

    return output_dir / output_filename
