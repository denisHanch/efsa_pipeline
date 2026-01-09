"""Configuration management for bioinformatics file validation."""

import json
import os
from pathlib import Path
from typing import List, Optional, Dict, Any
from dataclasses import dataclass

from validation_pkg.utils.formats import CodingType, GenomeFormat, ReadFormat, FeatureFormat
from validation_pkg.utils import file_handler, path_utils
from validation_pkg.logger import get_logger
from validation_pkg.exceptions import (
    ConfigurationError,
    FileNotFoundError as ValidationFileNotFoundError
)

# ===== Global Configuration Constants =====
# Only these fields can be specified in config.json "options" section
# These are the only settings that make sense to apply globally to all files
ALLOWED_GLOBAL_OPTIONS = {'threads', 'validation_level', 'logging_level'}
MAX_RECOMMENDED_THREADS = 16

# Default thread count when not specified in config
DEFAULT_THREADS = 8

@dataclass
class GenomeConfig:
    """Configuration for genome or plasmid files."""
    filename: str
    filepath: Path
    basename: str = None
    coding_type: CodingType = None
    detected_format: GenomeFormat = None
    output_dir: Path = None
    global_options: Dict[str, Any] = None

    def __post_init__(self):
        if self.global_options is None:
            self.global_options = {}

        # Extract basename (filename without extension)
        # Handle files with or without extensions safely
        if self.coding_type:
            # Compressed file: remove both compression and format extensions (e.g., .fasta.gz)
            parts = self.filename.rsplit('.', 2)
            self.basename = parts[0] if len(parts) >= 2 else self.filename
        else:
            # Uncompressed file: remove only format extension (e.g., .fasta)
            parts = self.filename.rsplit('.', 1)
            self.basename = parts[0] if len(parts) >= 2 else self.filename


@dataclass
class ReadConfig:
    """Configuration for sequencing read files."""
    filename: str
    filepath: Path
    basename: str = None
    ngs_type: str = None
    coding_type: CodingType = None
    detected_format: ReadFormat = None
    output_dir: Path = None
    global_options: Dict[str, Any] = None

    def __post_init__(self):
        if self.ngs_type not in ["illumina", "ont", "pacbio"]:
            raise ValueError(f"Invalid ngs_type: {self.ngs_type}")
        if self.global_options is None:
            self.global_options = {}

        # Extract basename (filename without extension)
        # Handle files with or without extensions safely
        if self.coding_type:
            # Compressed file: remove both compression and format extensions (e.g., .fastq.gz)
            parts = self.filename.rsplit('.', 2)
            self.basename = parts[0] if len(parts) >= 2 else self.filename
        else:
            # Uncompressed file: remove only format extension (e.g., .fastq)
            parts = self.filename.rsplit('.', 1)
            self.basename = parts[0] if len(parts) >= 2 else self.filename


@dataclass
class FeatureConfig:
    """Configuration for feature annotation files."""
    filename: str
    filepath: Path
    basename: str = None
    coding_type: CodingType = None
    detected_format: FeatureFormat = None
    output_dir: Path = None
    global_options: Dict[str, Any] = None

    def __post_init__(self):
        if self.global_options is None:
            self.global_options = {}

        # Extract basename (filename without extension)
        # Handle files with or without extensions safely
        if self.coding_type:
            # Compressed file: remove both compression and format extensions (e.g., .gff.gz)
            parts = self.filename.rsplit('.', 2)
            self.basename = parts[0] if len(parts) >= 2 else self.filename
        else:
            # Uncompressed file: remove only format extension (e.g., .gff)
            parts = self.filename.rsplit('.', 1)
            self.basename = parts[0] if len(parts) >= 2 else self.filename


class Config:
    """Main configuration container for validation workflow."""

    def __init__(self):
        # Required fields
        self.ref_genome: Optional[GenomeConfig] = None
        self.mod_genome: Optional[GenomeConfig] = None
        self.reads: List[ReadConfig] = []

        # Optional fields
        self.ref_plasmid: Optional[GenomeConfig] = None
        self.mod_plasmid: Optional[GenomeConfig] = None
        self.ref_feature: Optional[FeatureConfig] = None
        self.mod_feature: Optional[FeatureConfig] = None
        self.options: Dict[str, Any] = {}

        # Internal fields
        self.config_dir: Optional[Path] = None
        self.output_dir: Optional[Path] = None

    def get_threads(self) -> Optional[int]:
        """Get thread count from options, or None for auto-detection."""
        return self.options.get('threads')

    def get_logging_level(self) -> str:
        """Get logging level from options, or default to 'INFO'."""
        return self.options.get('logging_level', 'INFO')

    def __repr__(self):
        return (
            f"Config(\n"
            f"  ref_genome={self.ref_genome},\n"
            f"  mod_genome={self.mod_genome},\n"
            f"  reads={len(self.reads)} file(s),\n"
            f"  ref_plasmid={self.ref_plasmid},\n"
            f"  mod_plasmid={self.mod_plasmid},\n"
            f"  ref_feature={self.ref_feature},\n"
            f"  mod_feature={self.mod_feature}\n"
            f")"
        )


class ConfigManager:
    """Configuration file parser and validator."""
    
    @staticmethod
    def load(config_path: str) -> Config:
        """Load and validate configuration from JSON file."""
        logger = get_logger()
        logger.info(f"Loading configuration from: {config_path}")
        
        config_file = Path(config_path)
        
        if not config_file.exists():
            error_msg = f"Configuration file not found: {config_path}"
            logger.error(error_msg)
            raise ValidationFileNotFoundError(error_msg)

        # Load JSON
        try:
            logger.debug(f"Reading configuration file...")
            with open(config_file, 'r', encoding='utf-8') as f:
                data = json.load(f)
            logger.debug(f"Configuration file parsed successfully")
        except json.JSONDecodeError as e:
            error_msg = f"Invalid JSON in configuration file: {e}"
            logger.error(error_msg)
            logger.add_validation_issue(
                level='ERROR',
                category='configuration',
                message='Malformed JSON in config file',
                details={'file': str(config_path), 'error': str(e)}
            )
            raise ConfigurationError(error_msg) from e
        
        # Create Config object
        config = Config()
        config.config_dir = config_file.parent
        
        # Validate and parse
        try:
            logger.debug("Validating configuration structure...")
            ConfigManager._validate_required_fields(data)
            logger.debug("Parsing options...")
            ConfigManager._parse_options(data, config)

            # Reconfigure logging level based on config options (if specified)
            if config.options and 'logging_level' in config.options:
                logging_level = config.options['logging_level']
                logger.info(f"Applying logging level from config: {logging_level}")
                logger.reconfigure_level(console_level=logging_level)

            logger.debug("Setup base for outputs")
            ConfigManager._setup_output_directory(config)
            logger.debug("Parsing genome configurations...")
            ConfigManager._parse_genome_configs(data, config)
            logger.debug("Parsing reads configurations...")
            ConfigManager._parse_reads_configs(data, config)
            logger.debug("Parsing feature configurations...")
            ConfigManager._parse_feature_configs(data, config)

            logger.info("✓ Configuration loaded and validated successfully")
            return config
            
        except (ValueError, KeyError) as e:
            error_msg = f"Configuration validation failed: {e}"
            logger.error(error_msg)
            raise ConfigurationError(error_msg) from e

    @staticmethod
    def _validate_required_fields(data: dict):
        """Validate that required top-level fields exist (internal method)."""
        required = ['ref_genome_filename', 'mod_genome_filename', 'reads']
        
        for field in required:
            if field not in data:
                raise ValueError(f"Missing required field: {field}")
        
        # Reads must be non-empty
        if not data['reads']:
            raise ValueError("'reads' must be a non-empty list")

    @staticmethod
    def _parse_genome_configs(data: dict, config: Config):
        """Parse genome and plasmid configurations (internal method)."""
        # Required genomes
        config.ref_genome = ConfigManager._parse_genome_config(
            data['ref_genome_filename'], 'ref_genome_filename', config.config_dir, config.output_dir, config.options
        )
        config.mod_genome = ConfigManager._parse_genome_config(
            data['mod_genome_filename'], 'mod_genome_filename', config.config_dir, config.output_dir, config.options
        )

        # Optional plasmids
        if 'ref_plasmid_filename' in data and data['ref_plasmid_filename']:
            config.ref_plasmid = ConfigManager._parse_genome_config(
                data['ref_plasmid_filename'], 'ref_plasmid_filename', config.config_dir, config.output_dir, config.options
            )

        if 'mod_plasmid_filename' in data and data['mod_plasmid_filename']:
            config.mod_plasmid = ConfigManager._parse_genome_config(
                data['mod_plasmid_filename'], 'mod_plasmid_filename', config.config_dir, config.output_dir, config.options
            )

    @staticmethod
    def _parse_genome_config(value: Any, field_name: str, config_dir: Path, output_dir: Path, global_options: Dict[str, Any] = None) -> GenomeConfig:
        """Parse a genome configuration entry."""
        logger = get_logger()

        # Parse filename and extra fields using unified utility
        filename, extra = file_handler.parse_config_file_value(value, field_name)

        # Extract file-level global options (threads, validation_level only)
        filelvl_options = ConfigManager._merge_options(field_name, global_options, extra)

        # Resolve absolute path with security validation
        filepath = path_utils.resolve_filepath(config_dir, filename)

        if not filepath.exists():
            raise ValidationFileNotFoundError(
                f"The following file was not found: {filepath}\n")

        # Detect compression and format using unified utilities
        coding_type = file_handler.detect_compression_type(filepath)
        detected_format = file_handler.detect_file_format(filepath, GenomeFormat)

        return GenomeConfig(
            filename=filepath.name,
            filepath=filepath,
            coding_type=coding_type,
            detected_format=detected_format,
            output_dir=output_dir,
            global_options=filelvl_options,
        )

    @staticmethod
    def _parse_reads_configs(data: dict, config: Config):
        """Parse reads configuration (internal method)."""
        logger = get_logger()
        reads_data = data['reads']

        if not isinstance(reads_data, list):
            raise ValueError("'reads' must be a list")

        for idx, read_entry in enumerate(reads_data):
            if not isinstance(read_entry, dict):
                raise ValueError(f"reads[{idx}] must be a dict")
            
            try:
                # Get filename or directory
                directory = read_entry.get('directory')
                filename = read_entry.get('filename')

                if not (filename or directory):
                    raise ValueError("'filename' or 'directory' field is required")

                if filename and directory:
                    raise ValueError("Cannot specify both 'filename' and 'directory' - use one or the other")

                if filename:
                    config.reads.append(ConfigManager._parse_read_config(read_entry,'filename',config.config_dir, config.output_dir, config.options))

                if directory:
                    #   Resolve absolute path directory
                    dirpath = config.config_dir / directory
                    if not dirpath.exists():
                        raise ValidationFileNotFoundError(
                            f"Directory not found: {dirpath}"
                        )

                    if not dirpath.is_dir():
                        raise ValueError(f"Path is not a directory: {dirpath}")

                    files = [Path(f) for f in dirpath.iterdir()]
                    if not files:
                        logger.warning(
                            f"No valid read files found in directory: {dirpath}. "
                        )
                        continue

                    logger.info(f"Found {len(files)} read file(s) in directory: {dirpath}")

                    for file in files:
                        # Create a dict with filename and ngs_type from directory entry
                        file_entry = {
                            'filename': str(Path(directory) / file.name),
                            'ngs_type': read_entry.get('ngs_type')
                        }
                        # Include any other extra fields from the directory entry
                        for key, value in read_entry.items():
                            if key not in ['directory', 'filename', 'ngs_type']:
                                file_entry[key] = value

                        config.reads.append(ConfigManager._parse_read_config(file_entry, "", config.config_dir, config.output_dir, config.options))

            except ValueError as e:
                raise ValueError(f"Invalid reads[{idx}]: {e}")

    @staticmethod
    def _parse_read_config(value: Any, field_name: str, config_dir: Path, output_dir: Path, global_options: Dict[str, Any] = None) -> ReadConfig:
        """Parse a read configuration entry."""
        logger = get_logger()

        # Parse filename and extra fields using unified utility
        filename, extra = file_handler.parse_config_file_value(value, field_name)

        # Extract ngs_type first (it's a required field, not a global option)
        ngs_type = extra.pop('ngs_type', 'illumina')

        if not ngs_type:
            raise ValueError("Missing required 'ngs_type'")

        # Extract file-level global options (threads, validation_level only)
        filelvl_options = ConfigManager._merge_options(field_name,global_options,extra)

        # Resolve absolute path with security validation
        filepath = path_utils.resolve_filepath(config_dir, filename)

        if not filepath.exists():
            raise ValidationFileNotFoundError(
                f"The following file was not found: {filepath}\n")

        # Detect compression and format using unified utilities
        coding_type = file_handler.detect_compression_type(filepath)
        detected_format = file_handler.detect_file_format(filepath, ReadFormat)

        return ReadConfig(
            filename=filepath.name,
            filepath=filepath,
            ngs_type=ngs_type,
            coding_type=coding_type,
            detected_format=detected_format,
            output_dir=output_dir,
            global_options=filelvl_options
        )

    @staticmethod
    def _parse_feature_configs(data: dict, config: Config):
        """Parse feature configurations (internal method)."""
        if 'ref_feature_filename' in data and data['ref_feature_filename']:
            config.ref_feature = ConfigManager._parse_feature_config(
                data['ref_feature_filename'], 'ref_feature_filename', config.config_dir, config.output_dir, config.options
            )

        if 'mod_feature_filename' in data and data['mod_feature_filename']:
            config.mod_feature = ConfigManager._parse_feature_config(
                data['mod_feature_filename'], 'mod_feature_filename', config.config_dir, config.output_dir, config.options
            )

    @staticmethod
    def _parse_feature_config(value: Any, field_name: str, config_dir: Path, output_dir : Path, global_options: Dict[str, Any] = None) -> FeatureConfig:
        """Parse a single feature config entry and resolve paths."""
        logger = get_logger()

        # Parse filename and extra fields using unified utility
        filename, extra = file_handler.parse_config_file_value(value, field_name)

        # Extract file-level global options (threads, validation_level only)
        filelvl_options = ConfigManager._merge_options(field_name,global_options,extra)

        # Resolve absolute path with security validation
        filepath = path_utils.resolve_filepath(config_dir, filename)
        if not filepath.exists():
            raise ValidationFileNotFoundError(
                f"The following file was not found: {filepath}\n")

        # Detect compression and format using unified utilities
        coding_type = file_handler.detect_compression_type(filepath)
        detected_format = file_handler.detect_file_format(filepath, FeatureFormat)

        return FeatureConfig(
            filename=filepath.name,
            filepath=filepath,
            coding_type=coding_type,
            detected_format=detected_format,
            output_dir=output_dir,
            global_options=filelvl_options
        )

    @staticmethod
    def _parse_options(data: dict, config: Config):
        """Parse and validate global configuration options."""
        logger = get_logger()

        # Get options dict from config, or empty dict if not present
        if 'options' not in data:
            logger.debug("No options specified in config - using defaults")
            options = {}
        else:
            options = data['options']

            if not isinstance(options, dict):
                raise ConfigurationError("'options' must be a dictionary")

            # Validate that ONLY allowed fields are present
            invalid_keys = set(options.keys()) - ALLOWED_GLOBAL_OPTIONS

            if invalid_keys:
                allowed_list = ', '.join(f"'{k}'" for k in sorted(ALLOWED_GLOBAL_OPTIONS))
                raise ConfigurationError(
                    f"Invalid global options: {invalid_keys}. "
                    f"Only {allowed_list} are allowed in global options. "
                    f"Other settings should be specified per-file in the config."
                )

        # Parse threads option (with default)
        if 'threads' in options:
            threads = options['threads']

            # Allow null/None to mean auto-detect
            if threads is None:
                logger.info("Thread count: auto-detect (null specified)")
                config.options['threads'] = None
            elif isinstance(threads, int):
                # Validate positive integer
                if threads <= 0:
                    raise ConfigurationError(f"'threads' must be a positive integer, got {threads}")

                # Detect system CPU cores
                system_cores = os.cpu_count()
                if system_cores is None:
                    system_cores_msg = "unknown"
                else:
                    system_cores_msg = f"{system_cores} cores"

                # Warn if threads exceed system CPU cores
                if system_cores and threads > system_cores:
                    logger.warning(
                        f"Requested {threads} threads but system only has {system_cores} CPU cores. "
                        f"Performance may degrade due to context switching overhead. "
                        f"Consider using threads ≤ {system_cores} for optimal performance."
                    )

                # Warn if excessive (diminishing returns beyond MAX_RECOMMENDED_THREADS)
                if threads > MAX_RECOMMENDED_THREADS:
                    logger.warning(
                        f"Thread count {threads} is high - diminishing returns beyond {MAX_RECOMMENDED_THREADS} threads. "
                        f"Consider using 4-8 threads for optimal performance."
                    )

                logger.info(f"Global option: threads={threads} (system has {system_cores_msg} available)")
                config.options['threads'] = threads
            else:
                raise ConfigurationError(
                    f"'threads' must be an integer or null, got {type(threads).__name__}: {threads}"
                )
        else:
            logger.debug("Thread count: auto-detect (not specified)")

        # Parse validation_level option (with default)
        if 'validation_level' in options:
            validation_level = options['validation_level']

            VALID_LEVELS = {'strict', 'trust', 'minimal'}

            if not isinstance(validation_level, str):
                raise ConfigurationError(
                    f"'validation_level' must be a string, got {type(validation_level).__name__}: {validation_level}"
                )

            if validation_level not in VALID_LEVELS:
                raise ConfigurationError(
                    f"Invalid validation_level '{validation_level}'. "
                    f"Must be one of: {', '.join(sorted(VALID_LEVELS))}"
                )

            logger.info(f"Global option: validation_level={validation_level}")
            config.options['validation_level'] = validation_level
        else:
            logger.debug("validation_level not specified in global options")

        # Parse logging_level option
        if 'logging_level' in options:
            logging_level = options['logging_level']

            VALID_LOGGING_LEVELS = {'DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'}

            if not isinstance(logging_level, str):
                raise ConfigurationError(
                    f"'logging_level' must be a string, got {type(logging_level).__name__}: {logging_level}"
                )

            # Normalize to uppercase
            logging_level = logging_level.upper()

            if logging_level not in VALID_LOGGING_LEVELS:
                raise ConfigurationError(
                    f"Invalid logging_level '{logging_level}'. "
                    f"Must be one of: {', '.join(sorted(VALID_LOGGING_LEVELS))}"
                )

            logger.info(f"Global option: logging_level={logging_level}")
            config.options['logging_level'] = logging_level
        else:
            logger.debug("logging_level not specified in global options (default: INFO)")

    @staticmethod
    def _merge_options(field_name: str, global_options: Dict[str, Any], extra: Dict[str, Any]) -> Dict[str, Any]:        
        logger = get_logger()

        # Extract file-level global options (threads, validation_level only)
        filelvl_options = {}
        remaining_extra = {}

        filelvl_options.update(global_options)
        for key, value_item in extra.items():
            if key in ALLOWED_GLOBAL_OPTIONS:
                # Only log override if the key existed in global_options
                if key in filelvl_options:
                    logger.warning(
                        f"{field_name}: File-level option '{key}={value_item}' overrides "
                        f"global option '{key}={filelvl_options[key]}'"
                    )
                else:
                    logger.debug(
                        f"{field_name}: File-level option '{key}={value_item}' specified"
                    )
                filelvl_options[key] = value_item
            else:
                remaining_extra[key] = value_item

        # Log file-level options
        if filelvl_options:
            logger.debug(f"{field_name}: Found file-level options: {list(filelvl_options.keys())}")

        # Warn about remaining extra fields
        if remaining_extra:
            logger.warning(
                f"{field_name}: Found extra fields that will be ignored: {list(remaining_extra.keys())}"
            )

        return filelvl_options

    @staticmethod
    def _setup_output_directory(config: Config):
        """Set up output directory for validation results."""
        logger = get_logger()

        # Create output directory path
        output_dir = config.config_dir.parent / "valid"

        # Create directory if it doesn't exist
        try:
            output_dir.mkdir(parents=True, exist_ok=True)
            config.output_dir = output_dir
            logger.info(f"Output directory: {output_dir}")
            logger.debug(f"Output directory created/verified: {output_dir}")
        except OSError as e:
            error_msg = f"Failed to create output directory: {e}"
            logger.error(error_msg)
            raise ConfigurationError(error_msg) from e