"""Path resolution and security utilities for file operations."""

from pathlib import Path
from typing import Optional
from validation_pkg.exceptions import ConfigurationError

__all__ = [
    'resolve_filepath',
    'sanitize_path_component',
    'build_safe_output_dir',
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


def sanitize_path_component(component: str, allow_slashes: bool = False) -> str:
    """
    Sanitize a path component to prevent directory traversal attacks.

    Args:
        component: The path component to sanitize
        allow_slashes: Whether to allow forward slashes (for subdirectories)

    Returns:
        Sanitized path component

    Raises:
        ValueError: If component contains illegal characters or patterns

    Examples:
        >>> sanitize_path_component("valid_subdir")
        'valid_subdir'
        >>> sanitize_path_component("../etc/passwd")
        ValueError: Invalid path component: contains '..'
        >>> sanitize_path_component("/etc/passwd")
        ValueError: Invalid path component: contains absolute path
    """
    # Remove leading/trailing whitespace
    component = component.strip()

    # Check for empty after stripping
    if not component or component == '.':
        raise ValueError("Invalid path component: empty or current directory")

    # Block parent directory references
    if '..' in component:
        raise ValueError(f"Invalid path component: contains '..' (potential path traversal)")

    # Block absolute paths
    if component.startswith('/') or (len(component) > 1 and component[1] == ':'):
        raise ValueError(f"Invalid path component: absolute path not allowed")

    # Block null bytes
    if '\0' in component:
        raise ValueError("Invalid path component: contains null byte")

    # Platform-specific checks
    import platform
    if platform.system() == 'Windows':
        # Block Windows reserved characters
        invalid_chars = ['<', '>', '"', '|', '?', '*']
        for char in invalid_chars:
            if char in component:
                raise ValueError(f"Invalid path component: contains illegal character '{char}'")

        # Block Windows reserved names
        reserved_names = ['CON', 'PRN', 'AUX', 'NUL', 'COM1', 'COM2', 'COM3', 'COM4',
                         'COM5', 'COM6', 'COM7', 'COM8', 'COM9', 'LPT1', 'LPT2',
                         'LPT3', 'LPT4', 'LPT5', 'LPT6', 'LPT7', 'LPT8', 'LPT9']
        if component.upper() in reserved_names:
            raise ValueError(f"Invalid path component: Windows reserved name '{component}'")

    # If slashes not allowed, check for path separators
    if not allow_slashes:
        if '/' in component or '\\' in component:
            raise ValueError(f"Invalid path component: contains path separator")

    return component


def build_safe_output_dir(
    base_dir: Path,
    subdir_name: Optional[str] = None,
    create: bool = True
) -> Path:
    """
    Build output directory with path traversal protection.

    Args:
        base_dir: Base output directory (must exist)
        subdir_name: Optional subdirectory name (will be sanitized)
        create: Whether to create the directory if it doesn't exist

    Returns:
        Safe output directory path

    Raises:
        ValueError: If subdir_name contains illegal characters
        ConfigurationError: If path would escape base_dir

    Examples:
        >>> build_safe_output_dir(Path("/output"), "results")
        Path("/output/results")
        >>> build_safe_output_dir(Path("/output"), "../etc")
        ValueError: Invalid path component: contains '..'
    """
    # Start with base directory
    output_dir = base_dir

    # Add sanitized subdirectory if provided
    if subdir_name:
        # Sanitize the subdirectory name
        try:
            safe_subdir = sanitize_path_component(subdir_name, allow_slashes=False)
        except ValueError as e:
            raise ValueError(f"Invalid output subdirectory name: {e}")

        output_dir = output_dir / safe_subdir

        # Verify the result is still within base_dir using resolve_filepath
        # This provides defense-in-depth against any sanitization bypasses
        try:
            # resolve_filepath already validates the path stays within base
            output_dir = resolve_filepath(base_dir, safe_subdir)
        except ConfigurationError:
            raise ConfigurationError(
                f"Output subdirectory '{subdir_name}' would escape base directory '{base_dir}'"
            )

    # Create directory if requested
    if create:
        output_dir.mkdir(parents=True, exist_ok=True)

    return output_dir
