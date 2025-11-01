"""
Base settings infrastructure for validators.

Provides BaseSettings abstract base class with immutable update pattern,
deep copy support, and dictionary serialization for all validator Settings classes.

Features:
    - Immutable update pattern (prevents shared mutable state bugs)
    - Deep copy support
    - Dictionary serialization (to_dict/from_dict)
    - Field name validation (prevents typos)
    - Pretty printing for inspection
"""

from dataclasses import asdict, fields
from copy import deepcopy
from typing import Dict, Any
from abc import ABC

# ===== Base Settings Class =====
class BaseSettings(ABC):
    """
    Base class for all validator settings with common functionality.

    Features:
    - Deep copy support
    - Immutable update pattern (returns new instance)
    - Dictionary conversion for serialization
    - Pretty printing
    - from_dict() class method for loading from dictionaries
    - Validation of setting names (prevents typos)

    Usage:
        # Create settings
        settings = MyValidator.Settings()

        # Update (returns new instance)
        new_settings = settings.update(param1=value1, param2=value2)

        # Print to inspect
        print(settings)

        # Convert to dict
        settings_dict = settings.to_dict()

        # Load from dict
        settings = MyValidator.Settings.from_dict({'param1': value1})
    """

    def copy(self):
        """
        Return a deep copy of settings.

        Returns:
            New settings instance with copied values
        """
        return deepcopy(self)

    def update(self, **kwargs):
        """
        Update settings and return new instance (immutable pattern).

        This method prevents mutation of the original settings object,
        which avoids bugs from shared mutable state.

        Args:
            **kwargs: Settings to update

        Returns:
            New settings instance with updated values

        Raises:
            ValueError: If an unknown setting name is provided

        Example:
            >>> settings = GenomeValidator.Settings()
            >>> new_settings = settings.update(sequence_prefix="chr1")
            >>> settings.sequence_prefix  # Original unchanged
            None
            >>> new_settings.sequence_prefix
            'chr1'
        """
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

    def to_dict(self, include_unset: bool = True) -> Dict[str, Any]:
        """
        Convert settings to dictionary for serialization.

        Args:
            include_unset: If False, exclude fields with UNSET values

        Returns:
            Dictionary with all settings (or only set values if include_unset=False)

        Example:
            >>> settings.to_dict()
            {'plasmid_split': True, 'sequence_prefix': None, ...}
            >>> settings.to_dict(include_unset=False)
            {'plasmid_split': True}  # Excludes UNSET fields
        """
        result = asdict(self)
        if not include_unset:
            # Filter out UNSET values
            result = {k: v for k, v in result.items() if not isinstance(v, _UnsetType)}
        return result

    @classmethod
    def from_dict(cls, data: Dict[str, Any]):
        """
        Create settings instance from dictionary.

        Args:
            data: Dictionary with settings

        Returns:
            New settings instance

        Raises:
            ValueError: If dictionary contains unknown settings

        Example:
            >>> data = {'sequence_prefix': 'chr1', 'plasmid_split': False}
            >>> settings = GenomeValidator.Settings.from_dict(data)
        """
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
        """
        Pretty print settings for inspection.

        Returns:
            Formatted string representation

        Example:
            >>> print(settings)
            GenomeValidator.Settings:
              plasmid_split: True
              sequence_prefix: None
              min_sequence_length: 100
              ...
        """
        lines = [f"{self.__class__.__name__}:"]
        for key, value in self.to_dict().items():
            lines.append(f"  {key}: {value}")
        return '\n'.join(lines)

    def __repr__(self) -> str:
        """
        Return repr string for debugging.

        Returns:
            String representation suitable for debugging
        """
        params = ', '.join(f"{k}={v!r}" for k, v in self.to_dict().items())
        return f"{self.__class__.__name__}({params})"