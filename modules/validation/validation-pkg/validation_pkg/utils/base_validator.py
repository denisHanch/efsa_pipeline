"""Common validation helper functions."""
"""Abstract base class for all validators."""

from abc import ABC, abstractmethod
from pathlib import Path
from typing import IO, Any, Type, Optional

from validation_pkg.logger import get_logger
from validation_pkg.utils.base_settings import BaseSettings
from validation_pkg.exceptions import ValidationError, CompressionError
from validation_pkg.utils.file_handler import open_file_with_coding_type
from validation_pkg.utils.path_utils import strip_all_extensions


class BaseValidator(ABC):
    """
    Abstract base class providing common infrastructure for all validators.

    Subclasses must implement:
    - _validator_type property (returns 'genome', 'read', or 'feature')
    - OutputMetadata property (returns the metadata class)
    - _parse_file() method
    - _validate_sequences() method
    - _apply_edits() method
    - _save_output() method
    - _fill_output_metadata() method
    """
    # Subclasses must define their own Settings class
    Settings: Type[BaseSettings]

    def __init__(self, config: Any, settings: Optional[BaseSettings] = None):
        """
        Initialize common validator infrastructure.

        Args:
            config: Configuration object (genome_config, read_config, or feature_config)
            settings: Validator-specific settings (optional, uses defaults if not provided)
        """
        # Common infrastructure
        self.logger = get_logger()
        self.config = config

        # Extract common attributes from config
        self.output_dir = config.output_dir
        self.validation_level = config.global_options.get("validation_level")
        self.threads = config.global_options.get("threads")
        self.input_path = config.filepath
        
        # Initialize settings (subclass-specific)
        self.settings = settings if settings is not None else self.Settings()
        self.output_path = None

        # Initialize output metadata (subclass-specific)
        self.output_metadata = self.OutputMetadata()

        # Apply defaults
        self._initialize_defaults()

    def _initialize_defaults(self) -> None:
        """Set default values for validation_level and threads."""
        if not self.validation_level:
            self.validation_level = 'strict'

        if not self.threads:
            from validation_pkg.config_manager import DEFAULT_THREADS
            self.threads = DEFAULT_THREADS

    def run(self) -> Any:
        """
        Execute validation workflow with timing and error handling.

        Template method that orchestrates:
        1. Build output path (early)
        2. Start timer and log start
        3. Handle minimal mode OR call _run_validation()
        4. Stop timer and log completion
        5. Record file timing
        6. Set elapsed time in metadata
        7. Return metadata

        Returns:
            OutputMetadata object with validation results
        """
        # Build output path once at the beginning
        self.output_path = self._build_output_path()

        timer_name = f"{self._validator_type}_validation"
        self.logger.start_timer(timer_name)
        self.logger.info(f"Processing {self._validator_type} file: {self.config.filename}")
        self.logger.debug(
            f"Format: {self.config.detected_format}, "
            f"Compression: {self.config.coding_type}"
        )

        try:
            # Handle minimal mode separately
            if self.validation_level == 'minimal':
                output_path = self._handle_minimal_mode()
            else:
                # Trust/Strict modes - call subclass implementation
                output_path = self._run_validation()

            # Stop timer and log success
            elapsed = self.logger.stop_timer(timer_name)
            self.logger.info(
                f"✓ {self._validator_type.capitalize()} validation completed "
                f"in {elapsed:.2f}s"
            )

            # Record timing for report
            self.logger.add_file_timing(
                self.config.filename,
                self._validator_type,
                elapsed
            )

            # Fill metadata
            self._fill_output_metadata(output_path)
            self.output_metadata.elapsed_time = elapsed

            return self.output_metadata

        except Exception as e:
            self.logger.error(
                f"{self._validator_type.capitalize()} validation failed: {e}"
            )
            raise

    def _open_file(self, mode: str = 'rt') -> IO:
        """
        Open file with automatic decompression based on config.

        Args:
            mode: File open mode (default: 'rt' for text reading)

        Returns:
            File handle with automatic decompression

        Raises:
            CompressionError: If decompression fails
        """
        try:
            return open_file_with_coding_type(
                self.input_path,
                self.config.coding_type,
                mode
            )
        except CompressionError as e:
            self.logger.add_validation_issue(
                level='ERROR',
                category=self._validator_type,
                message=str(e),
                details={'file': str(self.input_path)}
            )
            raise

    def _fill_base_metadata(self, output_path: Path) -> None:
        """
        Fill common output metadata fields.

        Subclasses should call this first, then add their specific fields.

        Args:
            output_path: Path to output file
        """
        self.output_metadata.input_file = self.config.filename
        self.output_metadata.output_file = str(output_path) if output_path else None
        self.output_metadata.output_filename = output_path.name if output_path else None
        self.output_metadata.validation_level = self.validation_level

    def _build_output_filename(
        self,
        input_filename: str,
        output_format: str,
        coding_type: Optional[Any] = None,
        suffix: Optional[str] = None,
        input_path: Optional[Path] = None
    ) -> str:
        """
        Build output filename with consistent naming convention.

        Args:
            input_filename: Original input filename
            output_format: Output format extension (e.g., 'fasta', 'fastq', 'gff')
            coding_type: Optional compression type (adds .gz, .bz2, etc.)
            suffix: Optional custom suffix to add before extension
            input_path: Optional input Path object for accurate extension stripping

        Returns:
            Output filename with proper extensions

        Examples:
            >>> build_output_filename("genome.fa.gz", "fasta", CodingType.GZIP, "validated")
            'genome_validated.fasta.gz'
            >>> build_output_filename("reads.fq", "fastq", CodingType.NONE)
            'reads.fastq'
        """
        # Strip all extensions from input filename
        base_name = strip_all_extensions(input_filename, input_path)

        # Add custom suffix if provided
        if suffix:
            filename = f"{base_name}_{suffix}.{output_format}"
        else:
            filename = f"{base_name}.{output_format}"

        # Add compression extension if specified
        if coding_type:
            filename += coding_type.to_extension()

        return filename

    def _build_output_path(self) -> Path:
        """
        Build complete output path with consistent naming and directory structure.

        This is the main entry point for building output paths. It handles:
        - Safe subdirectory creation (with path traversal protection)
        - Extension stripping
        - Custom suffix addition
        - Compression extension handling

        Args:
            base_dir: Base output directory
            input_filename: Original input filename
            output_format: Output format extension (e.g., 'fasta', 'fastq', 'gff')
            coding_type: Optional compression type
            subdir_name: Optional subdirectory name (will be sanitized)
            filename_suffix: Optional suffix to add to filename
            input_path: Optional input Path for accurate extension handling
            create_dirs: Whether to create directories (default: True)

        Returns:
            Complete output path

        Raises:
            ValueError: If subdir_name contains illegal characters
            ConfigurationError: If path would escape base_dir

        Examples:
            >>> build_output_path(
            ...     Path("/output"),
            ...     "genome.fasta.gz",
            ...     "fasta",
            ...     CodingType.GZIP,
            ...     subdir_name="validated",
            ...     filename_suffix="filtered"
            ... )
            Path("/output/validated/genome_filtered.fasta.gz")
        """
        from validation_pkg.utils.path_utils import build_safe_output_dir

        # Build safe output directory (with path traversal protection)
        output_dir = build_safe_output_dir(self.output_dir, self.settings.output_subdir_name)

        # Build output filename
        output_filename = self._build_output_filename(
            self.config.filename,
            self._output_format,
            self.settings.coding_type,
            self.settings.output_filename_suffix,
            self.input_path
        )

        return output_dir / output_filename

    def _handle_minimal_mode(self) -> Path:
        """
        Handle minimal validation mode.

        Validates that input format and compression match requirements,
        then copies file without parsing.

        Returns:
            Path to output file

        Raises:
            ValidationError: If format or compression doesn't match requirements
        """
        self.logger.debug("Minimal mode - validating format and coding requirements")

        try:
            """
            Validate that minimal mode requirements are met.

            Minimal mode requires:
            1. Input format matches expected output format (no conversion)
            2. Input compression matches required output compression (no recompression)

            Args:
                detected_format: Detected format from input file
                expected_format: Expected format enum value
                input_coding: Detected coding type from input
                required_coding: Required coding type for output
                filename: Filename for error messages
                logger: Logger instance for validation issues
                category: Category for logging (genome/read/feature)

            Raises:
                ValidationError: If requirements not met
            """
            # Check format matches
            if self.config.detected_format != self._expected_format:
                error_msg = (
                    f'Minimal mode requires {self._expected_format} format, '
                    f'got {self.config.detected_format}. '
                    f'Use validation_level "trust" or "strict" to convert.'
                )
                self.logger.add_validation_issue(
                    level='ERROR',
                    category=self._validator_type,
                    message=error_msg,
                    details={
                        'file': self.config.filename,
                        'self.config.detected_format': str(self.config.detected_format),
                        'self._expected_format': str(self._expected_format)
                    }
                )
                raise ValidationError(error_msg)

            # Check coding matches
            if self.config.coding_type != self.settings.coding_type:
                error_msg = (
                    f'Minimal mode requires input coding to match output coding. '
                    f'Input: {self.config.coding_type}, Required: {self.settings.coding_type}. '
                    f'Use validation_level "trust" or "strict" to change compression.'
                )
                self.logger.add_validation_issue(
                    level='ERROR',
                    category=self._validator_type,
                    message=error_msg,
                    details={
                        'file': self.config.filename,
                        'input coding_type': str(self.config.coding_type),
                        'required coding_type': str(self.settings.coding_type)
                    }
                )
                raise ValidationError(error_msg)

        except ValidationError as e:
            # Re-raise as validator-specific exception
            raise self._get_validator_exception()(str(e)) from e

        # Use file_handler copy utility
        from validation_pkg.utils.file_handler import copy_file
        return copy_file(self.input_path, self.output_path, self.logger)

    # Abstract properties and methods that subclasses must implement

    @property
    @abstractmethod
    def _validator_type(self) -> str:
        """
        Return validator type string.

        Returns:
            'genome', 'read', or 'feature'
        """
        pass

    @property
    @abstractmethod
    def OutputMetadata(self) -> Type:
        """
        Return OutputMetadata class for this validator.

        Returns:
            OutputMetadata class (not instance)
        """
        pass

    @property
    @abstractmethod
    def _output_format(self) -> str:
        """
        Return output format string for build_output_path.

        Returns:
            Format string: 'fasta', 'fastq', or 'gff'
        """
        pass

    @property
    @abstractmethod
    def _expected_format(self) -> Any:
        """
        Return expected format for minimal mode validation.

        Returns:
            Format enum value (e.g., GenomeFormat.FASTA, ReadFormat.FASTQ, FeatureFormat.GFF)
        """
        pass

    @abstractmethod
    def _get_validator_exception(self) -> Type[Exception]:
        """
        Return the validator-specific exception class.

        Returns:
            Exception class: GenomeValidationError, ReadValidationError, or FeatureValidationError
        """
        pass

    @abstractmethod
    def _run_validation(self) -> Path:
        """
        Execute the validation workflow for trust/strict modes.

        This is where the main validation logic goes:
        - Parse file
        - Validate sequences
        - Apply edits
        - Save output

        Returns:
            Path to output file
        """
        pass

    @abstractmethod
    def _write_output(self) -> Path:
        """
        Write processed data to output file.

        Subclass-specific implementation for writing output.
        Uses self.output_path which is already built.

        Returns:
            Path to output file (should return self.output_path)
        """
        pass

    @abstractmethod
    def _fill_output_metadata(self, output_path: Path) -> None:
        """
        Populate output metadata with validation results.

        Should call self._fill_base_metadata(output_path) first,
        then add validator-specific fields.

        Args:
            output_path: Path to output file
        """
        pass
