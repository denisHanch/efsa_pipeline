"""
Feature file validator and processor.

This module provides comprehensive validation and processing for genomic feature annotation files
in GFF, GTF, and BED formats. It handles parsing, coordinate validation, format conversion, and
various editing operations for genomic annotations.

Key Features:
    - Multi-format support: GFF, GTF, and BED
    - Compression handling: gzip, bzip2, and uncompressed files
    - Three validation levels: strict (full validation), trust (fast validation), minimal (copy only)
    - Coordinate validation: start < end, no negative coordinates
    - BED to GFF conversion (0-based to 1-based coordinate system)
    - Feature sorting by genomic position
    - Sequence ID replacement with original tracking
    - Parallel compression support: automatic detection of pigz/pbzip2

Classes:
    Feature: Data container for genomic feature information
    FeatureValidator: Main validator class for feature annotation files
    FeatureValidator.Settings: Configuration dataclass for validation behavior
"""

from pathlib import Path
from typing import Optional, Dict, List, IO, Union
from dataclasses import dataclass
import shutil

from validation_pkg.logger import get_logger
from validation_pkg.utils.settings import BaseSettings
from validation_pkg.exceptions import (
    FeatureValidationError,
    FileFormatError,
    BedFormatError,
    GffFormatError,
    CompressionError,
)
from validation_pkg.utils.formats import CodingType as CT
from validation_pkg.utils.formats import FeatureFormat
from validation_pkg.utils.file_handler import open_compressed_writer


@dataclass
class OutputMetadata(BaseSettings):
    """
    Metadata returned from feature validation.

    Attributes:
        input_file: Full path to input file
        output_file: Full path to main output file
        output_filename: Name of output file
        num_features: Total number of features in output
        feature_types: List of unique feature types (gene, CDS, exon, etc.)
        sequence_ids: List of sequence IDs referenced by features (for inter-file validation)
        validation_level: Validation level used ('strict'/'trust'/'minimal')
        elapsed_time: Time taken for validation in seconds
    """
    input_file: str = None
    output_file: str = None
    output_filename: str = None
    num_features: int = None
    feature_types: List[str] = None
    sequence_ids: List[str] = None
    validation_level: str = None
    elapsed_time: float = None

    def __str__(self):
        parts = [f"Validation Level: {self.validation_level or 'N/A'}"]
        parts.append(f"Output File: {self.output_file or 'N/A'}")
        parts.append(f"Output Filename: {self.output_filename or 'N/A'}")

        if self.num_features is not None:
            parts.append(f"Number of Features: {self.num_features}")

        if self.feature_types:
            parts.append(f"Feature Types: {', '.join(self.feature_types)}")

        if self.sequence_ids:
            parts.append(f"Sequence IDs: {len(self.sequence_ids)} sequences referenced")

        if self.elapsed_time is not None:
            parts.append(f"Elapsed Time: {self.elapsed_time:.2f} seconds")

        return "\n".join(parts)
    
    
@dataclass
class Feature:
    """
    Represents a genomic feature (generic container for GFF/BED data).

    This class stores feature annotation data in a format-agnostic way,
    allowing seamless conversion between GFF and BED formats.

    Attributes:
        seqname: Sequence/chromosome name (e.g., 'chr1', 'contig_1')
        start: Start position (coordinate system depends on source format)
        end: End position (coordinate system depends on source format)
        feature_type: Feature type (e.g., 'gene', 'exon', 'CDS')
        score: Feature score (often '.' for no score)
        strand: Strand orientation ('+', '-', or '.')
        source: Source/program that generated the feature
        frame: Reading frame for CDS features (0, 1, 2, or '.')
        attributes: Semicolon-separated attribute string (GFF3 format)

    Properties:
        length: Calculated feature length (end - start)
    """
    seqname: str
    start: int
    end: int
    feature_type: str = ""
    score: str = "."
    strand: str = "+"
    source: str = "."
    frame: str = "."
    attributes: str = ""

    @property
    def length(self) -> int:
        """
        Calculate feature length in base pairs.

        Returns:
            Feature length (end - start)
        """
        return self.end - self.start


class FeatureValidator:
    """
    Validates and processes genomic feature annotation files in GFF, GTF, and BED formats.

    This validator provides three validation levels:
        - 'strict': Parse and validate all features, apply all edits (most thorough)
        - 'trust': Parse all features, validate only first 10, apply all edits (faster)
        - 'minimal': No parsing/validation, direct file copy (fastest, requires GFF format)

    Processing Workflow:
        1. Parse features from file (format auto-detected by ConfigManager)
        2. Validate coordinates: start < end, no negative values, no zero-length
        3. Apply edits: sort by position, replace sequence IDs
        4. Convert BED to GFF format (0-based → 1-based coordinates)
        5. Save as GFF3 with optional compression

    Key Features:
        - Coordinate system conversion: BED (0-based, half-open) → GFF (1-based, closed)
        - Feature sorting by genomic position (seqname, start, end)
        - Sequence ID replacement with original tracking in attributes
        - GFF3 attribute parsing and reconstruction
        - Graceful handling of malformed lines (logs warnings, continues parsing)

    Attributes:
        feature_config: FeatureConfig object with file path and format info
        output_dir: Directory for output files
        settings: Settings object controlling validation and processing behavior
        features: List of parsed Feature objects (populated during validation)

    Example:
        >>> from validation_pkg import ConfigManager, FeatureValidator
        >>> config = ConfigManager.load("config.json")
        >>> settings = FeatureValidator.Settings(sort_by_position=True, replace_id_with='chr1')
        >>> validator = FeatureValidator(config.ref_feature, config.output_dir, settings)
        >>> validator.validate()
    """

    @dataclass
    class Settings(BaseSettings):
        """
        Settings for feature validation and processing.

        Attributes:
            sort_by_position: Sort features by position (default: True)
            check_coordinates: Validate that start < end (default: True)
            allow_zero_length: Allow zero-length features (start == end) (default: False)
            allow_negative_coords: Allow negative coordinates (default: False)
            replace_id_with: Replace the seqname/chromosome field (column 1) for all features
                           with this value. Original seqname stored in Note attribute.
                           (e.g., 'chr1' replaces all seqnames with 'chr1') (default: None)
            coding_type: Output compression type: 'gz', 'bz2', or None (default: None)
            output_filename_suffix: Suffix to add to output filename (default: None)
            output_subdir_name: Subdirectory name for output files (default: None)
        """
        # Validation options
        sort_by_position: bool = True
        check_coordinates: bool = True

        # Editing options
        replace_id_with: Optional[str] = None

        # Output format
        coding_type: Optional[str] = None
        output_filename_suffix: Optional[str] = None
        output_subdir_name: Optional[str] = None


    def __init__(self, feature_config, settings: Optional[Settings] = None) -> None:
        """
        Initialize feature validator.

        Args:
            feature_config: FeatureConfig object from ConfigManager with file info
            settings: Settings object with validation parameters.
                If None, uses options from feature_config.options (threads, validation_level),
                otherwise uses defaults.
        """
        self.logger = get_logger()
       
        # From genome global configuration
        self.feature_config = feature_config
        self.output_dir = feature_config.output_dir
        self.validation_level = feature_config.global_options.get("validation_level")   
        self.threads = feature_config.global_options.get("threads") 
        self.input_path = feature_config.filepath

        # From settings
        self.settings = settings if settings is not None else self.Settings() 
        self.output_metadata = OutputMetadata()

        # Parsed data
        self.features: List[Feature] = []

        if not self.validation_level:
            self.validation_level = 'strict'    #   default global value

        if not self.threads:
            from validation_pkg.config_manager import DEFAULT_THREADS
            self.threads = DEFAULT_THREADS    #   default global value

    def _fill_output_metadata(self, output_path: Path) -> None:
        """
        Create OutputMetadata object based on validation level.

        Args:
            output_path: Path to the main output file

        """
        self.output_metadata.input_file = self.feature_config.filename
        self.output_metadata.output_file = str(output_path) if output_path else None
        self.output_metadata.output_filename = output_path.name if output_path else None
        self.output_metadata.num_features = len(self.features)
        self.output_metadata.feature_types = list(set(f.feature_type for f in self.features))
        self.output_metadata.sequence_ids = list(set(f.seqname for f in self.features))
        self.output_metadata.validation_level = self.validation_level

    def run(self) -> OutputMetadata:
        """
        Main validation and processing workflow.

        Uses feature_config data (format, compression) provided by ConfigManager.

        Returns:
            OutputMetadata: Metadata object containing:
                - output_file: Full path to output file
                - output_filename: Name of output file
                - num_features: Total number of features
                - feature_types: List of unique feature types
                - sequence_ids: List of sequence IDs (for inter-file validation)
                - validation_level: Validation mode used

        Raises:
            FeatureValidationError: If validation fails
        """
        self.logger.start_timer("feature_validation")
        self.logger.info(f"Processing feature file: {self.feature_config.filename}")
        self.logger.debug(f"Format: {self.feature_config.detected_format}, Compression: {self.feature_config.coding_type}")

        try:
            self._parse_file()

            self._validate_features()

            self._convert_to_gff()

            self._apply_edits()

            output_path = self._save_output()

            elapsed = self.logger.stop_timer("feature_validation")
            self.logger.info(f"✓ Feature validation completed in {elapsed:.2f}s")

            # Record file timing for report
            self.logger.add_file_timing(
                self.feature_config.filename,
                "feature",
                elapsed
            )

            # Create and return OutputMetadata
            self._fill_output_metadata(output_path)
            self.output_metadata.elapsed_time = elapsed
            return self.output_metadata

        except Exception as e:
            self.logger.error(f"Feature validation failed: {e}")
            raise

    def _open_file(self, mode: str = 'rt') -> IO:
        """
        Open file with automatic decompression based on feature_config.

        Uses centralized file opening from file_handler.py to eliminate code duplication.

        Args:
            mode: File opening mode (default: 'rt' for text read)

        Returns:
            File handle

        Raises:
            CompressionError: If file cannot be opened or decompressed
        """
        try:
            from validation_pkg.utils.file_handler import open_file_with_coding_type
            return open_file_with_coding_type(
                self.input_path,
                self.feature_config.coding_type,
                mode
            )
        except CompressionError as e:
            self.logger.add_validation_issue(
                level='ERROR',
                category='feature',
                message=str(e),
                details={'file': str(self.input_path)}
            )
            raise

    def _parse_file(self) -> None:
        """
        Parse feature file based on detected format from feature_config.

        Behavior depends on validation_level:
        - 'strict': Full parsing and validation of all features
        - 'trust': Parse all features, but validate only first 10
        - 'minimal': Skip parsing entirely
        """
        format_str = str(self.feature_config.detected_format)
        self.logger.debug(f"Parsing {format_str} file (validation_level={self.validation_level})...")

        # Minimal mode - skip parsing
        if self.validation_level == 'minimal':
            self.logger.info("Minimal validation mode - skipping file parsing")
            self.features = []
            return

        try:
            with self._open_file() as handle:
                if 'gff' in format_str.lower():
                    self.features = self._parse_gff(handle)
                elif 'bed' in format_str.lower():
                    self.features = self._parse_bed(handle)
                else:
                    raise FileFormatError(f"Unsupported feature format: {format_str}")

            # Validate we got features
            if not self.features:
                error_msg = f"No features found in {format_str} file"
                self.logger.add_validation_issue(
                    level='WARNING',
                    category='feature',
                    message=error_msg,
                    details={'file': self.feature_config.filename, 'format': format_str}
                )
                # Empty feature files are allowed, just warn

            if self.validation_level == 'trust':
                self.logger.info(f"Trust mode - parsed {len(self.features)} feature(s), will validate first 10 only")
            else:
                self.logger.debug(f"Parsed {len(self.features)} feature(s)")

        except (FileFormatError, FeatureValidationError):
            raise
        except Exception as e:
            error_msg = f"Failed to parse {format_str} file: {e}"

            if 'gff' in format_str.lower():
                exception_class = GffFormatError
            elif 'bed' in format_str.lower():
                exception_class = BedFormatError
            else:
                exception_class = FileFormatError

            self.logger.add_validation_issue(
                level='ERROR',
                category='feature',
                message=error_msg,
                details={
                    'file': self.feature_config.filename,
                    'format': format_str,
                    'error': str(e)
                }
            )
            raise exception_class(error_msg) from e

    def _normalize_strand(self, strand: str) -> str:
        """
        Normalize strand value to valid GFF3 format.

        Args:
            strand: Raw strand value from input file

        Returns:
            Normalized strand ('+', '-', or '.' for unknown)
        """
        return strand if strand in ['+', '-', '.'] else '.'

    def _skip_line(self, line: str) -> bool:
        """
        Check if line should be skipped during parsing.

        Args:
            line: Input line to check

        Returns:
            True if line is empty or a comment (starts with '#')
        """
        return not line or line.startswith('#')

    def _parse_gff(self, handle) -> List[Feature]:
        """
        Parse GFF/GTF format file.

        GFF format (tab-separated):
        seqname source feature start end score strand frame attributes

        Parses all features regardless of validation_level.

        Returns:
            List of Feature objects
        """
        features = []

        for line_num, line in enumerate(handle, start=1):
            line = line.strip()
            if self._skip_line(line):
                continue

            try:
                parts = line.split('\t')
                if len(parts) < 8:
                    self._log_parse_warning(
                        f"Invalid GFF line (expected 8-9 columns, got {len(parts)})",
                        line_num, line
                    )
                    continue

                feature = Feature(
                    seqname=parts[0],
                    source=parts[1],
                    feature_type=parts[2],
                    start=int(parts[3]),
                    end=int(parts[4]),
                    score=parts[5],
                    strand=self._normalize_strand(parts[6]),
                    frame=parts[7],
                    attributes=parts[8] if len(parts) > 8 else "",
                )
                features.append(feature)

            except (ValueError, IndexError) as e:
                self._log_parse_warning(f"Failed to parse GFF line: {e}", line_num)
                continue

        return features

    def _parse_bed(self, handle) -> List[Feature]:
        """
        Parse BED format file.

        BED format (tab-separated, minimum 3 columns):
        chrom start end [name] [score] [strand] [...]

        Parses all features regardless of validation_level.

        Returns:
            List of Feature objects
        """
        features = []

        for line_num, line in enumerate(handle, start=1):
            line = line.strip()
            if self._skip_line(line):
                continue

            try:
                parts = line.split('\t')
                if len(parts) < 3:
                    self._log_parse_warning(
                        f"Invalid BED line (expected at least 3 columns, got {len(parts)})",
                        line_num, line
                    )
                    continue

                feature = Feature(
                    seqname=parts[0],
                    start=int(parts[1]),
                    end=int(parts[2]),
                    attributes=parts[3] if len(parts) > 3 else "",
                    score=parts[4] if len(parts) > 4 else ".",
                    strand=self._normalize_strand(parts[5]) if len(parts) > 5 else '.',
                    source=".",
                    feature_type="region",
                    frame=".",
                )
                features.append(feature)

            except (ValueError, IndexError) as e:
                self._log_parse_warning(f"Failed to parse BED line: {e}", line_num)
                continue

        return features

    def _log_parse_warning(self, message: str, line_num: int, content: str = "") -> None:
        """Log a parsing warning."""
        details = {'line': line_num}
        if content:
            details['content'] = content[:100]
        self.logger.add_validation_issue(
            level='WARNING',
            category='feature',
            message=message,
            details=details
        )

    def _validate_features(self) -> None:
        """
        Validate parsed features.

        Behavior depends on validation_level:
        - 'strict': Validate all features
        - 'trust': Validate only first 10 features
        - 'minimal': Skip validation (no parsing done)

        Note: Coordinate validation is trivial math - parallelization overhead exceeds gains.
        """
        self.logger.debug("Validating features...")

        # Minimal mode - no features to validate
        if self.validation_level == 'minimal':
            self.logger.debug("Minimal mode - skipping feature validation")
            return

        # Determine how many features to validate
        if self.validation_level == 'trust':
            # Trust mode - validate only first 10 features
            features_to_validate = min(10, len(self.features))
            self.logger.debug(f"Trust mode - validating first {features_to_validate} of {len(self.features)} features")
        else:
            # Strict mode - validate all features
            features_to_validate = len(self.features)

        # Validate features
        if self.settings.check_coordinates:
            for idx in range(features_to_validate):
                self._validate_coordinates(self.features[idx], idx)

        self.logger.debug("✓ Feature validation passed")

    def _validate_coordinates(self, feature: Feature, idx: int) -> None:
        """
        Validate feature genomic coordinates.

        Performs three checks:
            1. Start must be less than or equal to end
            2. No negative coordinates allowed
            3. No zero-length features (start == end)

        Args:
            feature: Feature object to validate
            idx: Feature index (for error reporting)

        Raises:
            FeatureValidationError: If any validation check fails
        """
        if feature.start > feature.end:
            error_message = f"Feature #{idx}: start > end ({feature.seqname}:{feature.start}-{feature.end})"
            self.logger.add_validation_issue(
                level='ERROR',
                category='feature',
                message=error_message,
                details={'index': idx}
            )
            raise FeatureValidationError(error_message)

        if feature.start < 0 or feature.end < 0:
            error_message = f"Feature #{idx}: negative coordinates ({feature.seqname}:{feature.start}-{feature.end})"
            self.logger.add_validation_issue(
                level='ERROR',
                category='feature',
                message=error_message,
                details={'index': idx}
            )
            raise FeatureValidationError(error_message)

        if feature.start == feature.end:
            error_message = f"Feature #{idx}: zero-length ({feature.seqname}:{feature.start})"
            self.logger.add_validation_issue(
                level='ERROR',
                category='feature',
                message=error_message,
                details={'index': idx}
            )
            raise FeatureValidationError(error_message)

    def _parse_gff_attributes(self, attr_string: str) -> Dict[str, str]:
        """
        Parse GFF attributes string into dictionary.

        GFF attributes are semicolon-separated key=value pairs:
        "ID=gene1;Name=testGene;Note=description"

        Args:
            attr_string: GFF attributes string

        Returns:
            Dictionary of attribute key-value pairs

        Example:
            >>> attrs = self._parse_gff_attributes("ID=gene1;Name=test")
            >>> attrs
            {'ID': 'gene1', 'Name': 'test'}
        """
        if not attr_string or attr_string == '.':
            return {}

        attrs = {}
        # Split by semicolon
        for pair in attr_string.split(';'):
            pair = pair.strip()
            if not pair:
                continue

            # Split by first '=' only (values might contain '=')
            if '=' in pair:
                key, value = pair.split('=', 1)
                attrs[key.strip()] = value.strip()
            else:
                # Malformed attribute, store as-is
                self.logger.debug(f"Malformed GFF attribute (no '='): {pair}")

        return attrs

    def _build_gff_attributes(self, attr_dict: Dict[str, str]) -> str:
        """
        Build GFF attributes string from dictionary.

        Creates semicolon-separated key=value pairs, with ID first if present.

        Args:
            attr_dict: Dictionary of attribute key-value pairs

        Returns:
            GFF attributes string

        Example:
            >>> attr_string = self._build_gff_attributes({'ID': 'gene1', 'Name': 'test'})
            >>> attr_string
            'ID=gene1;Name=test'
        """
        if not attr_dict:
            return '.'

        # Build list of key=value pairs, with ID first
        parts = []

        # ID should come first (GFF3 convention)
        if 'ID' in attr_dict:
            parts.append(f"ID={attr_dict['ID']}")

        # Add all other attributes in sorted order (for consistency)
        for key in sorted(attr_dict.keys()):
            if key != 'ID':  # Skip ID, already added
                parts.append(f"{key}={attr_dict[key]}")

        return ';'.join(parts)

    def _apply_edits(self) -> None:
        """
        Apply editing specifications to features based on settings.

        Applies the following edits (if enabled in settings):
        - Sort features by position (sort_by_position)
        - Replace feature seqname/chromosome field (replace_id_with)

        Behavior depends on validation_level:
        - 'strict': Apply all edits to all features
        - 'trust': Apply all edits to parsed features (first 10)
        - 'minimal': Skip edits (file will be copied as-is)

        Note: Consider parallelizing ID replacement if profiling shows it's a bottleneck for large GFF files (>100K features).
        """
        self.logger.debug("Applying editing specifications from settings...")

        # Minimal mode - skip edits, file will be copied as-is
        if self.validation_level == 'minimal':
            self.logger.debug("Minimal mode - skipping edits")
            return

        # Sort by position (seqname, then start, then end)
        if self.settings.sort_by_position:
            original_count = len(self.features)
            self.features.sort(key=lambda f: (f.seqname, f.start, f.end))
            self.logger.debug(f"Sorted {original_count} features by position")

        # Replace feature seqnames (chromosome/sequence identifiers)
        if self.settings.replace_id_with:
            new_seqname = self.settings.replace_id_with
            self.logger.debug(f"Replacing feature seqnames with '{new_seqname}'...")

            for feature in self.features:
                # Store original seqname in attributes
                old_seqname = feature.seqname

                # Parse GFF attributes
                attrs = self._parse_gff_attributes(feature.attributes)

                # Store original seqname in Note field (GFF3 standard)
                note_text = f"Original_seqname:{old_seqname}"
                if 'Note' in attrs and attrs['Note']:
                    # Append to existing Note
                    attrs['Note'] = f"{attrs['Note']};{note_text}"
                else:
                    # Create new Note
                    attrs['Note'] = note_text

                # Rebuild attributes string
                feature.attributes = self._build_gff_attributes(attrs)

                # Replace seqname for both GFF and BED
                feature.seqname = new_seqname

            self.logger.debug(f"Replaced {len(self.features)} feature seqnames with '{new_seqname}'")

        self.logger.debug(f"✓ Edits applied, {len(self.features)} feature(s) remaining")

    def _convert_to_gff(self) -> None:
        """
        Convert features to GFF format if needed.

        BED and GFF use different coordinate systems:
            - BED: 0-based, half-open [start, end)
            - GFF: 1-based, closed [start, end]

        Conversion: BED [start, end) → GFF [start+1, end]

        Example:
            BED: chr1 100 200 (covers positions 100-199 in 0-based)
            GFF: chr1 . . 101 200 (covers positions 101-200 in 1-based)
        """
        if self.feature_config.detected_format == FeatureFormat.GFF:
            self.logger.debug("Keeping GFF format")
            return

        self.logger.debug(f"Converting {len(self.features)} BED features to GFF format...")

        # BED is 0-based half-open, GFF is 1-based closed
        # Convert: BED [start, end) -> GFF [start+1, end]
        for feature in self.features:
            feature.start += 1  # Convert from 0-based to 1-based

        self.logger.debug("✓ BED to GFF conversion complete")

    def _save_output(self) -> Path:
        """
        Save processed features to output directory using settings.

        Behavior depends on validation_level:
        - 'strict': Write features using parser (may convert format)
        - 'trust': Write parsed features (first 10 only) - NOT RECOMMENDED, use strict instead
        - 'minimal': Copy file as-is (preserve original format)

        Returns:
            Path to output file
        """
        self.logger.debug("Saving output file...")

        # Determine output directory (with optional subdirectory)
        output_dir = self.output_dir
        if self.settings.output_subdir_name:
            output_dir = output_dir / self.settings.output_subdir_name

        # Create output directory
        output_dir.mkdir(parents=True, exist_ok=True)

        # Generate output filename from original filename (without compression extension)
        # Get base name without any extensions
        base_name = self.feature_config.filename

        # Remove all suffixes (.gff.gz -> remove .gz then .gff)
        for suffix in self.input_path.suffixes:
            base_name = base_name.replace(suffix, '')
        # Add suffix if specified
        if self.settings.output_filename_suffix:
            output_filename = f"{base_name}_{self.settings.output_filename_suffix}.gff"
        else:
            output_filename = f"{base_name}.gff"

        # Add compression extension if requested (from settings, not input)
        # Support both string ('gz', 'bz2') and enum (CodingType.GZIP, CodingType.BZIP2)
        coding = self.settings.coding_type

        if coding in ('gz', 'gzip', CT.GZIP):
            output_filename += '.gz'
        elif coding in ('bz2', 'bzip2', CT.BZIP2):
            output_filename += '.bz2'

        # Minimal mode - copy file as-is without parsing
        # Required: GFF format + coding must match settings.coding_type
        if self.validation_level == 'minimal':
            self.logger.debug("Minimal mode - validating format and coding requirements")

            # Check format - must be GFF (reject BED files)
            if self.feature_config.detected_format != FeatureFormat.GFF:
                error_msg = f'Minimal mode requires GFF format, got {self.feature_config.detected_format}. Use validation_level "trust" or "strict" to convert.'
                self.logger.add_validation_issue(
                    level='ERROR',
                    category='feature',
                    message=error_msg,
                    details={'file': self.feature_config.filename, 'detected_format': str(self.feature_config.detected_format)}
                )
                raise FeatureValidationError(error_msg)

            # Check coding - must match settings.coding_type
            # Normalize both for comparison
            input_coding = self.feature_config.coding_type
            required_coding = coding  # This is already normalized from settings.coding_type

            # Map None to CT.NONE for comparison
            if input_coding is None:
                input_coding = CT.NONE
            if required_coding is None:
                required_coding = CT.NONE

            # Convert string to CT if needed
            if isinstance(required_coding, str):
                if required_coding in ('gz', 'gzip'):
                    required_coding = CT.GZIP
                elif required_coding in ('bz2', 'bzip2'):
                    required_coding = CT.BZIP2

            if input_coding != required_coding:
                error_msg = f'Minimal mode requires input coding to match output coding. Input: {input_coding}, Required: {required_coding}. Use validation_level "trust" or "strict" to change compression.'
                self.logger.add_validation_issue(
                    level='ERROR',
                    category='feature',
                    message=error_msg,
                    details={'file': self.feature_config.filename, 'input_coding': str(input_coding), 'required_coding': str(required_coding)}
                )
                raise FeatureValidationError(error_msg)

            # Use the output_filename already constructed (includes correct extension)
            output_path_minimal = output_dir / output_filename

            # Copy file
            self.logger.debug(f"Copying {self.input_path} to {output_path_minimal}")
            shutil.copy2(self.input_path, output_path_minimal)

            self.logger.info(f"Output saved: {output_path_minimal}")
            return output_path_minimal

        # Strict and Trust modes - write edited features
        output_path = output_dir / output_filename
        self.logger.debug(f"Writing output to: {output_path}")

        # Normalize coding type to CodingType enum
        coding_enum = CT.NONE
        if self.settings.coding_type in ('gz', 'gzip'):
            coding_enum = CT.GZIP
        elif self.settings.coding_type in ('bz2', 'bzip2'):
            coding_enum = CT.BZIP2

        # Use optimized compression writer
        with open_compressed_writer(output_path, coding_enum, threads=self.threads) as handle:
            # Write GFF header
            handle.write("##gff-version 3\n")

            # Write features
            for feature in self.features:
                line = '\t'.join([
                    feature.seqname,
                    feature.source,
                    feature.feature_type,
                    str(feature.start),
                    str(feature.end),
                    feature.score,
                    feature.strand,
                    feature.frame,
                    feature.attributes
                ])
                handle.write(line + '\n')

        self.logger.info(f"Output saved: {output_path}")
        return output_path