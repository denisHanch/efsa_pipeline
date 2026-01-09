"""Custom exceptions for the bioinformatics validation package."""


class ValidationError(Exception):
    """Base exception for all validation errors."""
    pass


class ConfigurationError(ValidationError):
    """Raised when there are errors in the configuration file."""
    pass


class FileNotFoundError(ValidationError):
    """Raised when required files are missing."""
    pass


class FileFormatError(ValidationError):
    """Base exception for file format errors."""
    pass


class FastaFormatError(FileFormatError):
    """Raised when FASTA file has invalid format."""
    pass


class GenBankFormatError(FileFormatError):
    """Raised when GenBank file has invalid format."""
    pass


class BedFormatError(FileFormatError):
    """Raised when BED file has invalid format."""
    pass


class GffFormatError(FileFormatError):
    """Raised when GFF/GTF file has invalid format."""
    pass


class FastqFormatError(FileFormatError):
    """Raised when FASTQ file has invalid format."""
    pass


class BamFormatError(FileFormatError):
    """Raised when BAM file has invalid format."""
    pass


class CompressionError(ValidationError):
    """Raised when there are errors decompressing files."""
    pass


class GenomeValidationError(ValidationError):
    """Raised when genome validation fails."""
    pass


class FeatureValidationError(ValidationError):
    """Raised when feature validation fails."""
    pass


class ReadValidationError(ValidationError):
    """Raised when read validation fails."""
    pass


class InterFileValidationError(ValidationError):
    """Raised when inter-file consistency checks fail."""
    pass
