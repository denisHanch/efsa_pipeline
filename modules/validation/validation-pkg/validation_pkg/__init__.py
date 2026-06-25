"""
Validation package for genomic data files (genomes, reads).

Supports FASTA/GenBank, FASTQ/BAM with gzip/bzip2 compression.
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
from validation_pkg.validators.interfile_read import ReadXReadSettings, readxread_validation
from validation_pkg.validators.interfile_genome import GenomeXGenomeSettings, genomexgenome_validation
from validation_pkg.utils.logger import setup_logging, get_logger

__all__ = [
    # Configuration
    'ConfigManager',
    'Config',

    # Validators
    'GenomeValidator',
    'ReadValidator',

    # Inter-file Validation
    'ReadXReadSettings',
    'readxread_validation',
    'GenomeXGenomeSettings',
    'genomexgenome_validation',

    # Logging
    'setup_logging',
    'get_logger',

    # Version info
    '__version__',
    '__author__',
    '__license__',
]
