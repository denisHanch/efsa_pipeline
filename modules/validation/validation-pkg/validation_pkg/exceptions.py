"""
Custom exceptions for the bioinformatics validation package.

Provides a hierarchical exception system for type-specific error handling:

Exception hierarchy:
    ValidationError (base)
    ├── ConfigurationError (config file issues)
    ├── FileNotFoundError (missing files)
    ├── FileFormatError (invalid file formats)
    │   ├── FastaFormatError
    │   ├── GenBankFormatError
    │   ├── BedFormatError
    │   ├── GffFormatError
    │   ├── FastqFormatError
    │   └── BamFormatError
    ├── CompressionError (decompression failures)
    ├── GenomeValidationError (genome validation failures)
    ├── FeatureValidationError (feature validation failures)
    ├── ReadValidationError (read validation failures)
    └── InterFileValidationError (inter-file consistency errors)

Usage:
    Catch ValidationError to handle all package-specific errors:

    try:
        validate_genome(config)
    except ValidationError as e:
        print(f"Validation failed: {e}")
"""


class ValidationError(Exception):
    """
    Base exception for all validation errors.
    
    All custom exceptions in this package inherit from this.
    """
    pass


class ConfigurationError(ValidationError):
    """
    Raised when there are errors in the configuration file.
    
    Examples:
    - Missing required fields
    - Invalid values (output_coding, ngs_type)
    - Malformed JSON
    """
    pass


class FileNotFoundError(ValidationError):
    """
    Raised when required files are missing.
    
    Note: Inherits from ValidationError, not built-in FileNotFoundError
    to maintain our exception hierarchy.
    
    Examples:
    - Genome file doesn't exist
    - Read directory not found
    """
    pass


class FileFormatError(ValidationError):
    """
    Base exception for file format errors.
    
    Raised when a file exists but has invalid format.
    """
    pass


class FastaFormatError(FileFormatError):
    """
    Raised when FASTA file has invalid format.
    
    Examples:
    - Missing header line (>)
    - Invalid characters in sequence
    - Empty sequences
    """
    pass


class GenBankFormatError(FileFormatError):
    """
    Raised when GenBank file has invalid format.
    
    Examples:
    - Missing LOCUS line
    - Invalid feature table
    - Missing ORIGIN section
    """
    pass


class BedFormatError(FileFormatError):
    """
    Raised when BED file has invalid format.
    
    Examples:
    - Wrong number of columns
    - Invalid coordinates (start > end)
    - Non-numeric score
    """
    pass


class GffFormatError(FileFormatError):
    """
    Raised when GFF/GTF file has invalid format.
    
    Examples:
    - Wrong number of columns
    - Invalid coordinates
    - Missing required attributes
    """
    pass


class FastqFormatError(FileFormatError):
    """
    Raised when FASTQ file has invalid format.
    
    Examples:
    - Missing quality scores
    - Length mismatch between sequence and quality
    - Invalid quality encoding
    """
    pass


class BamFormatError(FileFormatError):
    """
    Raised when BAM file has invalid format.
    
    Examples:
    - Missing quality scores
    - Length mismatch between sequence and quality
    - Invalid quality encoding
    """
    pass


class CompressionError(ValidationError):
    """
    Raised when there are errors decompressing files.
    
    Examples:
    - Corrupted gzip file
    - Unsupported compression format
    """
    pass


class GenomeValidationError(ValidationError):
    """
    Raised when genome validation fails.
    
    Examples:
    - Duplicate sequence IDs
    - Invalid sequence composition
    - Empty genome file
    """
    pass


class FeatureValidationError(ValidationError):
    """
    Raised when feature validation fails.
    
    Examples:
    - Overlapping features (when not allowed)
    - Invalid feature type
    - Missing required attributes
    """
    pass


class ReadValidationError(ValidationError):
    """
    Raised when read validation fails.
    
    Examples:
    - Empty read file
    - Inconsistent paired-end reads
    - Invalid read names
    """
    pass


class InterFileValidationError(ValidationError):
    """
    Raised when inter-file consistency checks fail.
    
    Examples:
    - Feature coordinates outside genome bounds
    - Feature references non-existent sequence
    - Reads map to sequences not in genome
    """
    pass
