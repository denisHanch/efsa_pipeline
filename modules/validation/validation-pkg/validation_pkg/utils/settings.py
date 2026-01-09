"""Base settings infrastructure for validators."""

from dataclasses import asdict, fields
from copy import deepcopy
from typing import Dict, Any
from abc import ABC

__all__ = [
    'BaseSettings',
]

# ===== Base Settings Class =====
class BaseSettings(ABC):
    """Base class for all validator settings with common functionality."""

    def copy(self):
        """Return a deep copy of settings."""
        return deepcopy(self)

    def update(self, **kwargs):
        """Update settings and return new instance (immutable pattern)."""
        new_settings = self.copy()

        # Get list of valid field names
        valid_fields = {f.name for f in fields(new_settings)}

        # Check for unknown fields
        unknown = set(kwargs.keys()) - valid_fields
        if unknown:
            allowed = ', '.join(sorted(valid_fields))
            unknown_str = ', '.join(f"'{k}'" for k in sorted(unknown))
            raise ValueError(
                f"Unknown setting(s) {unknown_str} for {self.__class__.__name__}. "
                f"Allowed settings: {allowed}"
            )

        # Update fields
        for key, value in kwargs.items():
            setattr(new_settings, key, value)

        return new_settings

    def to_dict(self) -> Dict[str, Any]:
        """Convert settings to dictionary for serialization."""
        return asdict(self)

    @classmethod
    def from_dict(cls, data: Dict[str, Any]):
        """Create settings instance from dictionary."""
        # Get valid field names
        valid_fields = {f.name for f in fields(cls)}

        # Check for unknown fields
        unknown = set(data.keys()) - valid_fields
        if unknown:
            allowed = ', '.join(sorted(valid_fields))
            unknown_str = ', '.join(f"'{k}'" for k in sorted(unknown))
            raise ValueError(
                f"Unknown setting(s) {unknown_str} for {cls.__name__}. "
                f"Allowed settings: {allowed}"
            )

        return cls(**data)

    def __str__(self) -> str:
        """Pretty print settings for inspection."""
        lines = [f"{self.__class__.__name__}:"]
        for key, value in self.to_dict().items():
            lines.append(f"  {key}: {value}")
        return '\n'.join(lines)

    def __repr__(self) -> str:
        """Return repr string for debugging."""
        params = ', '.join(f"{k}={v!r}" for k, v in self.to_dict().items())
        return f"{self.__class__.__name__}({params})"
    
