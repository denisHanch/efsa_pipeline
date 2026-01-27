"""
Bioinformatics Validation Package
==================================

A comprehensive validation package for genomic data files used in bioinformatics
pipelines. Validates and standardizes genome files, sequencing reads, and
feature annotations with detailed error reporting and quality checks.

Supported File Types
-------------------
- **Genome files:** FASTA (.fasta, .fa, .fna), GenBank (.gb, .gbk, .genbank)
- **Read files:** FASTQ (.fastq, .fq), BAM (.bam)
- **Feature files:** GFF (.gff, .gff3, .gtf), BED (.bed)
- **Compression:** gzip (.gz), bzip2 (.bz2)

Features
--------
- Format detection and validation
- Automatic compression/decompression
- Sequence quality checks (duplicate IDs, invalid characters, empty sequences)
- Coordinate validation for features
- Format conversion (e.g., GenBank to FASTA, BED to GFF)
- Plasmid splitting for bacterial genomes
- BAM to FASTQ conversion
- Detailed validation reports with structured error logging
- Flexible configuration via JSON

Quick Start
-----------

**Functional API with Settings:**

>>> from validation_pkg import (
...     ConfigManager, ValidationReport, GenomeValidator,
...     validate_genome, validate_reads, genomexgenome_validation
... )
>>> config = ConfigManager.load("config.json")
>>> report = ValidationReport("logs/report.txt")
>>>
>>> # Configure and validate genomes
>>> ref_settings = GenomeValidator.Settings(plasmids_to_one=True, replace_id_with='chr')
>>> ref_result = validate_genome(config.ref_genome, ref_settings)
>>> report.write(ref_result, file_type="genome")
>>>
>>> # Validate reads (using defaults)
>>> reads_results = validate_reads(config.reads)
>>> report.write(reads_results, file_type="read")
>>> report.flush()

Configuration Example
--------------------

Create a config.json file:

.. code-block:: json

    {
      "output_dir": "./output",
      "ref_genome_filename": {"filename": "reference.fasta"},
      "mod_genome_filename": {"filename": "modified.fasta"},
      "reads": [
        {"filename": "sample1.fastq.gz", "ngs_type": "illumina"},
        {"filename": "sample2.fastq.gz", "ngs_type": "illumina"}
      ],
      "ref_feature_filename": {"filename": "annotations.gff"}
    }

See docs/config_guide.md for detailed configuration options.

Package Structure
----------------
- validation_pkg.config_manager: Configuration loading and parsing
- validation_pkg.validators: Individual file validators
    - genome_validator: Genome file validation
    - read_validator: Read file validation
    - feature_validator: Feature annotation validation
- validation_pkg.utils: Utility functions and formats
- validation_pkg.logger: Structured logging system

Workflow
--------
1. **Load configuration:** Parse JSON config and detect file formats
2. **Validate genomes:** Check reference and modified genome files
3. **Validate features:** Validate feature annotations (optional)
4. **Validate reads:** Process all sequencing read files

Error Handling
-------------
The package uses a hierarchical exception system:

- ValidationError: Base exception for all validation errors
    - GenomeValidationError: Genome-specific errors
        - FastaFormatError
        - GenBankFormatError
    - ReadValidationError: Read-specific errors
        - FastqFormatError
        - BamFormatError
    - FeatureValidationError: Feature-specific errors
    - ConfigurationError: Config file errors
    - FileNotFoundError: Missing file errors

All errors are logged with structured details for debugging.

Version Information
------------------
"""

__version__ = "0.1.0"
__author__ = "Dominika Bohuslavova"
__license__ = "EUPL-1.2 license"

# Public API exports
from validation_pkg.config_manager import ConfigManager, Config
from validation_pkg.validators.genome_validator import GenomeValidator
from validation_pkg.validators.genome_validator import GenomeOutputMetadata
from validation_pkg.validators.read_validator import ReadValidator
from validation_pkg.validators.read_validator import ReadOutputMetadata
from validation_pkg.validators.feature_validator import FeatureValidator
from validation_pkg.validators.feature_validator import FeatureOutputMetadata
from validation_pkg.validators.interfile_read import ReadXReadSettings, readxread_validation
from validation_pkg.validators.interfile_genome import GenomeXGenomeSettings, genomexgenome_validation
from validation_pkg.logger import setup_logging, get_logger
from validation_pkg.report import ValidationReport

# Functional API imports
from typing import Optional, List

# ============================================================================
# Functional API - Simplified wrapper functions
# ============================================================================

def validate_genome(
    genome_config,
    settings: Optional[GenomeValidator.Settings] = None
) -> GenomeOutputMetadata:
    """
    Validate a genome file with optional custom settings.

    This is a simplified wrapper around GenomeValidator for easier usage.
    """
    validator = GenomeValidator(genome_config, settings)
    return validator.run()

def validate_read(
    read_config,
    settings: Optional[ReadValidator.Settings] = None
) -> ReadOutputMetadata:
    """
    Validate a single read file with optional custom settings.

    This is a simplified wrapper around ReadValidator for easier usage.
    """
    validator = ReadValidator(read_config, settings)
    output_metadata = validator.run()
    return output_metadata

def validate_reads(
    read_configs: List,
    settings: Optional[ReadValidator.Settings] = None
) -> List[ReadOutputMetadata]:
    """Validate multiple read files with optional custom settings."""
    results = []
    for read_config in read_configs:
        validator = ReadValidator(read_config, settings)
        result = validator.run()
        results.append(result)
    return results

def validate_genomes(
    genome_configs: List,
    settings: Optional[GenomeValidator.Settings] = None
) -> List[GenomeOutputMetadata]:
    """Validate multiple genome files with optional custom settings."""
    results = []
    for genome_config in genome_configs:
        validator = GenomeValidator(genome_config, settings)
        result = validator.run()
        results.append(result)
    return results

def validate_feature(
    feature_config,
    settings: Optional[FeatureValidator.Settings] = None
) -> FeatureOutputMetadata:
    """
    Validate a feature annotation file with optional custom settings.

    This is a simplified wrapper around FeatureValidator for easier usage.
    """
    validator = FeatureValidator(feature_config, settings)
    return validator.run()

def validate_features(
    feature_configs: List,
    settings: Optional[FeatureValidator.Settings] = None
) -> List[FeatureOutputMetadata]:
    """Validate multiple feature annotation files with optional custom settings."""
    results = []
    for feature_config in feature_configs:
        validator = FeatureValidator(feature_config, settings)
        result = validator.run()
        results.append(result)
    return results

__all__ = [
    # Configuration
    'ConfigManager',
    'Config',

    # Validators
    'GenomeValidator',
    'ReadValidator',
    'FeatureValidator',

    # Functional API (Primary Interface)
    'validate_genome',
    'validate_genomes',
    'validate_read',
    'validate_reads',
    'validate_feature',
    'validate_features',

    # Inter-file Validation
    'ReadXReadSettings',
    'readxread_validation',
    'GenomeXGenomeSettings',
    'genomexgenome_validation',

    # Logging
    'setup_logging',
    'get_logger',
    'ValidationReport'

    # Version info
    '__version__',
    '__author__',
    '__license__',
]
