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
    'copy_file',

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
    'convert_file_compression',
]


# ===== Cache for Tool Availability =====
# Protected by _TOOL_CACHE_LOCK for thread-safe access
_TOOL_CACHE = {}
_TOOL_CACHE_LOCK = threading.Lock()
_LOGGER_INITIALIZED = False

# ===== Compression Tool Configuration =====
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

# ===== Command Arguments for Tools =====
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


def convert_file_compression(
    input_path: Union[str, Path],
    output_path: Union[str, Path],
    input_coding: CodingType,
    output_coding: CodingType,
    threads: int = None
) -> None:
    """Convert file from one compression type to another."""
    input_path = Path(input_path)
    output_path = Path(output_path)

    # If coding types match, just copy the file
    if input_coding == output_coding:
        import shutil
        shutil.copy2(input_path, output_path)
        return

    # Map (input, output) pairs to conversion functions
    conversion_map = {
        (CodingType.BZIP2, CodingType.GZIP): bz2_to_gz,
        (CodingType.NONE, CodingType.GZIP): none_to_gz,
        (CodingType.GZIP, CodingType.NONE): gz_to_none,
        (CodingType.BZIP2, CodingType.NONE): bz2_to_none,
        (CodingType.NONE, CodingType.BZIP2): none_to_bz2,
        (CodingType.GZIP, CodingType.BZIP2): gz_to_bz2,
    }

    conversion_func = conversion_map.get((input_coding, output_coding))
    if conversion_func:
        conversion_func(input_path, output_path, threads=threads)
    else:
        raise CompressionError(
            f"Unsupported compression conversion: {input_coding} -> {output_coding}"
        )


def copy_file(
    input_path: Union[str, Path],
    output_path: Union[str, Path],
    logger=None
) -> Path:
    """Copy file with optional logging."""
    input_path = Path(input_path)
    output_path = Path(output_path)

    if logger:
        logger.debug(f"Copying {input_path} to {output_path}")

    shutil.copy2(input_path, output_path)

    if logger:
        logger.info(f"File copied: {output_path}")

    return output_path


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

