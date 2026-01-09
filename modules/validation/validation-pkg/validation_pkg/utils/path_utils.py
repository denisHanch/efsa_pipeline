"""Path resolution and security utilities for file operations."""

from pathlib import Path
from validation_pkg.exceptions import ConfigurationError

__all__ = [
    'resolve_filepath',
]


def resolve_filepath(base_dir: Path, filename: str) -> Path:
    """Resolve filepath relative to base directory with path traversal protection."""
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
