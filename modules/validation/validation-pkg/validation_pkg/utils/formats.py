"""File format and compression type enumerations."""

from enum import Enum
from pathlib import Path

__all__ = [
    'CodingType',
    'GenomeFormat',
    'ReadFormat',
    'OrganismType',
    'ValidationLevel',
    'LoggingLevel',
    'NgsType',
]


class CodingType(Enum):
    """Supported compression types for genomic data files."""
    GZIP = "gzip"
    BZIP2 = "bzip2"
    NONE = "none"

    def to_extension(self) -> str:
        """Return file extension for this compression type."""
        extension_map = {
            CodingType.GZIP: '.gz',
            CodingType.BZIP2: '.bz2',
            CodingType.NONE: ''
        }
        return extension_map[self]

    @classmethod
    def normalize(cls, value):
        """Convert various input formats to CodingType enum."""
        # If already a CodingType, return as-is
        if isinstance(value, cls):
            return value

        # Handle None
        if value is None:
            return cls.NONE

        # Normalize the input
        value_lower = str(value).lower().strip()

        # Remove leading dot if present
        if value_lower.startswith('.'):
            value_lower = value_lower[1:]

        # Extension mapping
        extension_map = {
            'gz': cls.GZIP,
            'gzip': cls.GZIP,
            'bz2': cls.BZIP2,
            'bzip2': cls.BZIP2,
            'none': cls.NONE,
            '': cls.NONE
        }

        # Try to find matching enum
        if value_lower in extension_map:
            return extension_map[value_lower]

        # Check if it's a filename with extension
        if '.' in value_lower:
            path = Path(value)
            ext = path.suffix.lower()
            if ext:
                # Remove dot and try again
                return cls.normalize(ext)

        # Default to NONE if nothing matches
        return cls.NONE

    @classmethod
    def _missing_(cls, value):
        """Called when enum lookup fails - allows flexible input formats."""
        return cls.normalize(value)


class GenomeFormat(Enum):
    FASTA = "fasta"
    GENBANK = "genbank"

    def to_biopython(self) -> str:
        """Return format string for BioPython SeqIO."""
        return self.value
    
    def to_extension(self) -> str:
        """Return file extension for this genome format."""
        return f'.{self.value}'
    
    @classmethod
    def _missing_(cls, value):
        """Handle flexible input formats for GenomeFormat."""
        value_lower = str(value).lower().strip()
        
        # Remove leading dot
        if value_lower.startswith('.'):
            value_lower = value_lower[1:]
        
        # Extension mapping
        extension_map = {
            'fa': cls.FASTA,
            'fasta': cls.FASTA,
            'fna': cls.FASTA,
            'genbank': cls.GENBANK,
            'gb': cls.GENBANK,
            'gbk': cls.GENBANK,
        }
        
        # Direct match
        if value_lower in extension_map:
            return extension_map[value_lower]
        
        # If it looks like a filename, extract extension
        if '.' in value_lower:
            ext = Path(value).suffix.lower()[1:]  # Remove dot
            if ext in extension_map:
                return extension_map[ext]
        
        raise ValueError(f"'{value}' is not a valid {cls.__name__}")


class ReadFormat(Enum):
    FASTQ = "fastq"
    BAM = "bam"

    def to_biopython(self) -> str:
        """Return format string for BioPython SeqIO."""
        return self.value

    def to_extension(self) -> str:
        """Return file extension for this read format."""
        return f'.{self.value}'
    
    @classmethod
    def _missing_(cls, value):
        """Handle flexible input formats for ReadFormat."""
        value_lower = str(value).lower().strip()
        
        # Remove leading dot
        if value_lower.startswith('.'):
            value_lower = value_lower[1:]
        
        # Extension mapping
        extension_map = {
            'fq': cls.FASTQ,
            'fastq': cls.FASTQ,
            'bam': cls.BAM,
        }
        
        # Direct match
        if value_lower in extension_map:
            return extension_map[value_lower]
        
        # If it looks like a filename, extract extension
        if '.' in value_lower:
            ext = Path(value).suffix.lower()[1:]  # Remove dot
            if ext in extension_map:
                return extension_map[ext]
        
        raise ValueError(f"'{value}' is not a valid {cls.__name__}")



class OrganismType(Enum):
    """Supported organism types for genomic validation."""
    PROKARYOTE = "prokaryote"
    EUKARYOTE = "eukaryote"

    @classmethod
    def normalize(cls, value):
        """Convert various input formats to OrganismType enum."""
        if isinstance(value, cls):
            return value
        if value is None:
            return cls.PROKARYOTE
        value_lower = str(value).lower().strip()
        mapping = {
            'prokaryote': cls.PROKARYOTE,
            'eukaryote': cls.EUKARYOTE,
        }
        if value_lower in mapping:
            return mapping[value_lower]
        raise ValueError(
            f"'{value}' is not a valid OrganismType. Must be one of: prokaryote, eukaryote"
        )

    @classmethod
    def _missing_(cls, value):
        return cls.normalize(value)


class ValidationLevel(Enum):
    """Supported validation depth levels."""
    STRICT = "strict"
    TRUST = "trust"
    MINIMAL = "minimal"

    @classmethod
    def normalize(cls, value):
        """Convert various input formats to ValidationLevel enum."""
        if isinstance(value, cls):
            return value
        if value is None:
            return cls.TRUST
        value_lower = str(value).lower().strip()
        mapping = {
            'strict': cls.STRICT,
            'trust': cls.TRUST,
            'minimal': cls.MINIMAL,
        }
        if value_lower in mapping:
            return mapping[value_lower]
        raise ValueError(
            f"'{value}' is not a valid ValidationLevel. Must be one of: strict, trust, minimal"
        )

    @classmethod
    def _missing_(cls, value):
        return cls.normalize(value)


class LoggingLevel(Enum):
    """Supported logging verbosity levels."""
    DEBUG = "DEBUG"
    INFO = "INFO"
    WARNING = "WARNING"
    ERROR = "ERROR"

    @classmethod
    def normalize(cls, value):
        """Convert various input formats to LoggingLevel enum."""
        if isinstance(value, cls):
            return value
        if value is None:
            return cls.INFO
        value_upper = str(value).upper().strip()
        mapping = {
            'DEBUG': cls.DEBUG,
            'INFO': cls.INFO,
            'WARNING': cls.WARNING,
            'ERROR': cls.ERROR,
        }
        if value_upper in mapping:
            return mapping[value_upper]
        raise ValueError(
            f"'{value}' is not a valid LoggingLevel. Must be one of: DEBUG, INFO, WARNING, ERROR"
        )

    @classmethod
    def _missing_(cls, value):
        return cls.normalize(value)


class NgsType(Enum):
    """Supported next-generation sequencing technology types."""
    ILLUMINA = "illumina"
    ONT = "ont"
    PACBIO_HIFI = "pacbio-hifi"
    PACBIO_CLR = "pacbio-clr"

    @classmethod
    def normalize(cls, value):
        """Convert various input formats to NgsType enum."""
        if isinstance(value, cls):
            return value
        if value is None:
            return None
        value_lower = str(value).lower().strip()
        mapping = {
            'illumina': cls.ILLUMINA,
            'ont': cls.ONT,
            'pacbio-hifi': cls.PACBIO_HIFI,
            'pacbio-clr': cls.PACBIO_CLR,
        }
        if value_lower in mapping:
            return mapping[value_lower]
        raise ValueError(
            f"'{value}' is not a valid NgsType. Must be one of: illumina, ont, pacbio-hifi, pacbio-clr"
        )

    @classmethod
    def _missing_(cls, value):
        return cls.normalize(value)
