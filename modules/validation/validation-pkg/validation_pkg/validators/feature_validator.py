"""Feature file validator and processor for GFF, GTF, and BED formats."""

from pathlib import Path
from typing import Optional, List, Type, Any
from dataclasses import dataclass
import tempfile
import subprocess

from validation_pkg.utils.base_settings import BaseOutputMetadata, BaseValidatorSettings
from validation_pkg.exceptions import FeatureValidationError
from validation_pkg.utils.formats import CodingType as CT, FeatureFormat
from validation_pkg.utils.file_handler import open_compressed_writer, check_tool_available
from validation_pkg.utils.base_validator import BaseValidator

@dataclass
class OutputMetadata(BaseOutputMetadata):
    """Metadata returned from feature validation."""
    # Feature-specific fields
    num_features: int = None
    feature_types: List[str] = None
    sequence_ids: List[str] = None

    def format_statistics(self, indent: str = "    ", input_settings: dict = None) -> list[str]:
        """
        Format feature-specific statistics for report output.

        Args:
            indent: Indentation string (default: 4 spaces)
            input_settings: Optional settings dict (not currently used but kept for consistency)

        Returns:
            List of formatted strings with feature statistics
        """
        lines = []

        # Helper to format values
        def format_value(value):
            if isinstance(value, float):
                return f"{value:.2f}"
            elif isinstance(value, int) and value > 999:
                return f"{value:,}"
            return str(value)

        # Iterate through all fields, skipping common ones
        skip_fields = {'input_file', 'output_file', 'output_filename', 'validation_level', 'elapsed_time'}
        special_fields = {'feature_types', 'sequence_ids'}

        data = self.to_dict()

        for key, value in data.items():
            if key in skip_fields or value is None:
                continue

            # Special handling for specific fields
            if key == 'feature_types' and isinstance(value, list):
                if len(value) <= 5:
                    lines.append(f"{indent}feature_types: {', '.join(value)}")
                else:
                    lines.append(f"{indent}feature_types: {', '.join(list(value)[:5])}, ... (+{len(value)-5} more)")

            elif key == 'sequence_ids' and isinstance(value, list):
                if len(value) <= 3:
                    lines.append(f"{indent}sequence_ids: {', '.join(value)}")
                else:
                    lines.append(f"{indent}sequence_ids: {value[0]}, {value[1]}, ... (+{len(value)-2} more)")

            elif key not in special_fields:
                # Generic field formatting
                formatted_value = format_value(value)
                lines.append(f"{indent}{key}: {formatted_value}")

        return lines

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
    """Genomic feature container for GFF/BED data."""
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
        """Feature length in base pairs."""
        return self.end - self.start


class FeatureValidator(BaseValidator):
    """Validates and processes GFF, GTF, and BED feature annotation files."""

    @dataclass
    class Settings(BaseValidatorSettings):
        """Settings for feature validation and processing."""
        # Feature-specific settings
        sort_by_position: bool = True
        check_coordinates: bool = True
        replace_id_with: Optional[str] = None

        def __post_init__(self):
            """Normalize settings after initialization."""
            # Normalize coding_type from base class
            self._normalize_coding_type()

    def __init__(self, feature_config, settings: Optional[Settings] = None) -> None:
        # Call base class initialization
        super().__init__(feature_config, settings)

        # Keep feature_config for type safety and specific access
        self.feature_config = feature_config

        # Validator-specific data
        self.features: List[Feature] = []

    # Required abstract properties and methods from BaseValidator

    @property
    def _validator_type(self) -> str:
        """Return validator type string."""
        return 'feature'

    @property
    def OutputMetadata(self) -> Type:
        """Return OutputMetadata class for this validator."""
        return OutputMetadata

    @property
    def _output_format(self) -> str:
        """Return output format string for build_output_path."""
        return 'gff'

    @property
    def _expected_format(self) -> Any:
        """Return expected format for minimal mode validation."""
        return FeatureFormat.GFF

    def _get_validator_exception(self) -> Type[Exception]:
        """Return the validator-specific exception class."""
        return FeatureValidationError

    def _fill_output_metadata(self, output_path: Path) -> None:
        """Populate output metadata with validation results."""
        # Fill common fields from base class
        self._fill_base_metadata(output_path)

        # Add feature-specific fields
        self.output_metadata.num_features = len(self.features) if self.features else None
        self.output_metadata.feature_types = list(set(f.feature_type for f in self.features)) if self.features else None
        self.output_metadata.sequence_ids = list(set(f.seqname for f in self.features)) if self.features else None

    def _run_validation(self) -> Path:
        """Execute validation and processing workflow for trust/strict modes."""
        # Trust/Strict modes - full validation and processing
        self._parse_input()
        self._edit_features()
        output_path = self._write_output()
        return output_path

    def _parse_gff(self, handle) -> List[Feature]:
        """Parse GFF3 format file."""
        features = []
        progress_interval = 50000

        for line in handle:
            line = line.strip()
            if not line or line.startswith('#'):
                continue

            parts = line.split('\t')
            if len(parts) != 9:
                self.logger.debug(f"Skipping malformed line: {line[:50]}...")
                continue

            try:
                feature = Feature(
                    seqname=parts[0],
                    source=parts[1],
                    feature_type=parts[2],
                    start=int(parts[3]),
                    end=int(parts[4]),
                    score=parts[5],
                    strand=parts[6],
                    frame=parts[7],
                    attributes=parts[8],
                )
                features.append(feature)

                if len(features) % progress_interval == 0:
                    self.logger.info(f"Parsing progress: {len(features):,} features parsed...")

            except (ValueError, IndexError) as e:
                self.logger.debug(f"Failed to parse GFF line: {e}")
                continue

        return features

    def _parse_input(self) -> None:
        """Parse and convert input file to GFF3 using gffread."""

        self.logger.debug("Parsing input using gffread...")

        temp_files = []
        try:
            temp_input = self._prepare_input_file()
            if temp_input != self.input_path:
                temp_files.append(temp_input)

            temp_gff3_file = tempfile.NamedTemporaryFile(mode='w', suffix='.gff3', delete=False)
            temp_gff3_file.close()
            temp_gff3 = Path(temp_gff3_file.name)
            temp_files.append(temp_gff3)

            self._run_gffread(temp_input, temp_gff3)

            with open(temp_gff3, 'r') as f:
                self.features = self._parse_gff(f)

            self.logger.info(f"Parsed {len(self.features)} features")

        finally:
            for temp_file in temp_files:
                if temp_file and isinstance(temp_file, Path):
                    temp_file.unlink(missing_ok=True)

    def _edit_features(self) -> None:
        """Apply feature edits (sorting, ID replacement) in strict mode."""
        if self.validation_level != 'strict':
            self.logger.debug(f"{self.validation_level.capitalize()} mode - skipping feature editing")
            return

        if self.settings.sort_by_position:
            self.logger.debug("Sorting features by position...")
            self.features.sort(key=lambda f: (f.seqname, f.start, f.end))
            self.logger.info("✓ Features sorted by position")

        if self.settings.replace_id_with:
            new_id = self.settings.replace_id_with
            self.logger.debug(f"Replacing sequence IDs with '{new_id}'...")

            for feature in self.features:
                feature.seqname = new_id

            self.logger.info(f"✓ Replaced sequence IDs with '{new_id}'")

    def _run_gffread(
        self,
        input_file: Path,
        output_file: Path
    ) -> Path:
        """Run gffread to parse and validate input file. """

        if not check_tool_available('gffread'):
            raise FeatureValidationError(
                "gffread tool required."
            )

        cmd = [
            'gffread',
            str(input_file),
            '-O',
            '-v',
            '-E',
            '-o', str(output_file)
        ]

        self.logger.debug(f"Running gffread: {' '.join(cmd)}")

        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=True
            )

            if result.stderr:
                for line in result.stderr.strip().split('\n'):
                    if line:
                        self.logger.debug(f"gffread: {line}")

            self.logger.info("✓ gffread processing complete")
            return output_file

        except subprocess.CalledProcessError as e:
            error_msg = f"gffread failed: {e.stderr}"
            self.logger.error(error_msg)
            raise FeatureValidationError(error_msg) from e

    def _prepare_input_file(self) -> Path:
        """Decompress input file if needed for gffread processing."""

        if self.feature_config.coding_type == CT.NONE:
            return self.input_path

        suffix = self.feature_config.detected_format.to_extension()
        temp_file = tempfile.NamedTemporaryFile(mode='w', suffix=suffix, delete=False)

        self.logger.debug(f"Decompressing {self.input_path} to {temp_file.name}")

        with self._open_file('rt') as input_handle:
            with temp_file as output_handle:
                for line in input_handle:
                    output_handle.write(line)

        return Path(temp_file.name)

    def _write_output(self) -> Path:
        """Write processed features to output file."""
        self.logger.debug(f"Writing output to: {self.output_path}")

        with open_compressed_writer(self.output_path, self.settings.coding_type, threads=self.threads) as handle:
            handle.write("##gff-version 3\n")

            feature_count = len(self.features)
            if feature_count > 10000:
                self.logger.info(f"Writing {feature_count} features to output file...")

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

        self.logger.info(f"Output saved: {self.output_path}")
        return self.output_path