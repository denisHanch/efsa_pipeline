"""File format and compression type enumerations."""

from enum import Enum
from pathlib import Path

__all__ = [
    'CodingType',
    'GenomeFormat',
    'ReadFormat',
    'FeatureFormat',
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
    def _missing_(cls, value):
        """Called when enum lookup fails - allows flexible input formats."""
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
                return cls._missing_(ext)

        # Default to NONE if nothing matches
        return cls.NONE


class GenomeFormat(Enum):
    FASTA = "fasta"
    GENBANK = "genbank"

    def to_biopython(self) -> str:
        """Return format string for BioPython SeqIO."""
        return self.value
    
    def to_extension(self) -> str:
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
            'faa': cls.FASTA,
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


class FeatureFormat(Enum):
    GFF = "gff"
    GTF = "gtf"
    BED = "bed"

    def to_biopython(self) -> str:
        """Return format string for file parsing."""
        return self.value
    
    def to_extension(self) -> str:
        return f'.{self.value}'

    @classmethod
    def _missing_(cls, value):
        """Handle flexible input formats for FeatureFormat."""
        value_lower = str(value).lower().strip()
        
        # Remove leading dot
        if value_lower.startswith('.'):
            value_lower = value_lower[1:]
        
        # Extension mapping
        extension_map = {
            'gff': cls.GFF,
            'gff3': cls.GFF,
            'gtf': cls.GTF,
            'gff2': cls.GTF,
            'bed': cls.BED,
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