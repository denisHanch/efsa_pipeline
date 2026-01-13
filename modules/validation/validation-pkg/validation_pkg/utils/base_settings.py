"""Base settings infrastructure for validators."""

from dataclasses import asdict, fields, dataclass
from copy import deepcopy
from typing import Dict, Any, Optional
from abc import ABC
from validation_pkg.utils.formats import CodingType

__all__ = [
    'BaseSettings',
    'BaseOutputMetadata',
    'BaseValidatorSettings',
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

        # Update fields with automatic type normalization
        for key, value in kwargs.items():
            # Get the field type annotation
            field_info = next((f for f in fields(new_settings) if f.name == key), None)
            if field_info:
                # Check if this field is a CodingType enum and normalize it
                field_type = field_info.type
                # Handle Optional[CodingType] or CodingType annotations
                if 'CodingType' in str(field_type):
                    # Import here to avoid circular imports
                    from validation_pkg.utils.formats import CodingType
                    value = CodingType.normalize(value)

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

        # Normalize CodingType fields before creating instance
        normalized_data = {}
        for key, value in data.items():
            field_info = next((f for f in fields(cls) if f.name == key), None)
            if field_info and 'CodingType' in str(field_info.type):
                # Import here to avoid circular imports
                from validation_pkg.utils.formats import CodingType
                value = CodingType.normalize(value)
            normalized_data[key] = value

        return cls(**normalized_data)

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


# ===== Base Output Metadata =====
@dataclass
class BaseOutputMetadata(BaseSettings):
    """Base metadata class with common fields for all validator output metadata."""
    input_file: str = None
    output_file: str = None
    output_filename: str = None
    validation_level: str = None
    elapsed_time: float = None

    def format_common_fields(self, indent: str = "  ") -> list[str]:
        """
        Format common fields that are present in all validators.
        Returns a list of formatted strings ready for report output.

        Args:
            indent: Indentation string (default: 2 spaces)

        Returns:
            List of formatted strings (input file, output file, elapsed time)
        """
        lines = []

        if self.input_file is not None:
            lines.append(f"{indent}Input:  {self.input_file}")

        if self.output_file is not None:
            lines.append(f"{indent}Output: {self.output_file}")

        if self.elapsed_time is not None:
            lines.append(f"{indent}Time:   {self.elapsed_time:.2f}s")

        return lines


# ===== Base Validator Settings =====
@dataclass
class BaseValidatorSettings(BaseSettings):
    """Base settings class with common fields for all validators."""
    coding_type: Optional[CodingType] = None 
    output_filename_suffix: Optional[str] = None
    output_subdir_name: Optional[str] = None

    def _normalize_coding_type(self):
        """Normalize coding_type to handle strings and None."""
        self.coding_type = CodingType.normalize(self.coding_type)

