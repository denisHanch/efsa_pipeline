"""Shared formatting helpers for validation metadata output."""


def format_metadata_value(value) -> str:
    """Format a metadata value for human-readable report output."""
    if isinstance(value, float):
        return f"{value:.2f}"
    if isinstance(value, int) and value > 999:
        return f"{value:,}"
    return str(value)
