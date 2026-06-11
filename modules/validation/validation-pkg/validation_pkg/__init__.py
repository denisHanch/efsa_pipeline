"""
Validation package for genomic data files (genomes, reads, feature annotations).

Supports FASTA/GenBank, FASTQ/BAM, GFF/GTF/BED with gzip/bzip2 compression.
See docs/CONFIG_GUIDE.md for configuration options.
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
from validation_pkg.utils.logger import setup_logging, get_logger
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
    """Validate a genome file with optional custom settings."""
    validator = GenomeValidator(genome_config, settings)
    return validator.run()

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

def validate_feature(
    feature_config,
    settings: Optional[FeatureValidator.Settings] = None
) -> FeatureOutputMetadata:
    """Validate a feature annotation file with optional custom settings."""
    validator = FeatureValidator(feature_config, settings)
    return validator.run()

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
    'validate_reads',
    'validate_feature',

    # Inter-file Validation
    'ReadXReadSettings',
    'readxread_validation',
    'GenomeXGenomeSettings',
    'genomexgenome_validation',

    # Logging
    'setup_logging',
    'get_logger',
    'ValidationReport',

    # Version info
    '__version__',
    '__author__',
    '__license__',
]
