"""
Path resolution and security utilities for file operations.

This module provides utilities for safe path resolution with security protections
against path traversal attacks.
"""

from pathlib import Path
from validation_pkg.exceptions import ConfigurationError

__all__ = [
    'resolve_filepath',
]


def resolve_filepath(base_dir: Path, filename: str) -> Path:
    """
    Resolve filepath relative to base directory with path traversal protection.

    This function prevents path traversal attacks by ensuring the resolved path
    stays within the base directory. It handles relative paths with '..' and
    symlinks by resolving to absolute canonical paths, then validating the
    resolved path is still within the allowed base directory.

    Args:
        base_dir: Base directory (e.g., config_dir) - files must be under this
        filename: Filename or relative path from config file

    Returns:
        Resolved absolute Path object

    Raises:
        ConfigurationError: If path traverses outside base directory

    Examples:
        >>> resolve_filepath(Path("/home/user/project"), "genome.fasta")
        Path("/home/user/project/genome.fasta")

        >>> resolve_filepath(Path("/home/user/project"), "../../../etc/passwd")
        ConfigurationError: Path traversal detected

    Security:
        - Prevents directory traversal attacks (../..)
        - Resolves symlinks to their targets
        - Validates resolved path is within base_dir
    """
    # Resolve to absolute canonical paths
    filepath = (base_dir / filename).resolve()
    base_dir_resolved = base_dir.resolve()

    # Check if filepath is within base_dir
    # Use try/except for Python 3.7/3.8 compatibility (is_relative_to added in 3.9)
    try:
        filepath.relative_to(base_dir_resolved)
    except ValueError:
        # Path is outside base directory - security violation
        raise ConfigurationError(
            f"Path traversal detected: '{filename}' resolves outside config directory.\n"
            f"Resolved path: {filepath}\n"
            f"Config directory: {base_dir_resolved}\n"
            f"Only files within the config directory are allowed."
        )

    return filepath
