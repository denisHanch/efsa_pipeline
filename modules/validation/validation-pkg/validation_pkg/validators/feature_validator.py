"""Feature file validator and processor for GFF, GTF, and BED formats."""

from pathlib import Path
from typing import Optional, List, IO
from dataclasses import dataclass
import shutil

from validation_pkg.logger import get_logger
from validation_pkg.utils.settings import BaseSettings
from validation_pkg.exceptions import (
    FeatureValidationError,
    CompressionError,
)
from validation_pkg.utils.formats import CodingType as CT
from validation_pkg.utils.formats import FeatureFormat
from validation_pkg.utils.file_handler import open_compressed_writer, check_tool_available, open_file_with_coding_type


@dataclass
class OutputMetadata(BaseSettings):
    """Metadata returned from feature validation."""
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


class FeatureValidator:
    """Validates and processes GFF, GTF, and BED feature annotation files."""

    @dataclass
    class Settings(BaseSettings):
        """Settings for feature validation and processing."""
        sort_by_position: bool = True
        check_coordinates: bool = True
        replace_id_with: Optional[str] = None
        coding_type: Optional[str] = None
        output_filename_suffix: Optional[str] = None
        output_subdir_name: Optional[str] = None

    def __init__(self, feature_config, settings: Optional[Settings] = None) -> None:
        self.logger = get_logger()
        self.feature_config = feature_config
        self.output_dir = feature_config.output_dir
        self.validation_level = feature_config.global_options.get("validation_level")
        self.threads = feature_config.global_options.get("threads")
        self.input_path = feature_config.filepath
        self.settings = settings if settings is not None else self.Settings()
        self.output_metadata = OutputMetadata()
        self.features: List[Feature] = []

        if not self.validation_level:
            self.validation_level = 'strict'

    def _fill_output_metadata(self, output_path: Path) -> None:
        """Populate output metadata with validation results."""
        self.output_metadata.input_file = self.feature_config.filename
        self.output_metadata.output_file = str(output_path) if output_path else None
        self.output_metadata.output_filename = output_path.name if output_path else None
        self.output_metadata.num_features = len(self.features) if self.features else None
        self.output_metadata.feature_types = list(set(f.feature_type for f in self.features)) if self.features else None
        self.output_metadata.sequence_ids = list(set(f.seqname for f in self.features)) if self.features else None
        self.output_metadata.validation_level = self.validation_level

    def run(self) -> OutputMetadata:
        """Execute validation and processing workflow."""
        self.logger.start_timer("feature_validation")
        self.logger.info(f"Processing feature file: {self.feature_config.filename}")
        self.logger.debug(f"Format: {self.feature_config.detected_format}, Compression: {self.feature_config.coding_type}")

        try:
            self._parse_input()
            self._edit_features()
            output_path = self._save_output()

            elapsed = self.logger.stop_timer("feature_validation")
            self.logger.info(f"✓ Feature validation completed in {elapsed:.2f}s")

            self.logger.add_file_timing(
                self.feature_config.filename,
                "feature",
                elapsed
            )

            self._fill_output_metadata(output_path)
            self.output_metadata.elapsed_time = elapsed
            return self.output_metadata

        except Exception as e:
            self.logger.error(f"Feature validation failed: {e}")
            raise

    def _open_file(self, mode: str = 'rt') -> IO:
        """Open file with automatic decompression based on feature_config."""
        try:
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
        if self.validation_level == 'minimal':
            return

        self.logger.debug("Parsing input using gffread...")

        temp_files = []
        try:
            temp_input = self._prepare_input_file()
            if temp_input != self.input_path:
                temp_files.append(temp_input)

            import tempfile
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
        """
        Run gffread to parse and validate input file.

        Args:
            input_file: Input GFF/GTF/BED file (uncompressed)
            output_file: Output GFF3 file

        Returns:
            Path to output GFF3 file

        Raises:
            FeatureValidationError: If gffread fails
        """
        import subprocess

        if not check_tool_available('gffread'):
            raise FeatureValidationError(
                "gffread tool required. Install via: conda install -c bioconda gffread"
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
        import tempfile

        if self.feature_config.coding_type == CT.NONE:
            return self.input_path

        format_suffixes = {
            FeatureFormat.GFF: '.gff',
            FeatureFormat.GTF: '.gtf',
            FeatureFormat.BED: '.bed',
        }
        suffix = format_suffixes.get(self.feature_config.detected_format, '.gff')

        temp_file = tempfile.NamedTemporaryFile(mode='w', suffix=suffix, delete=False)

        self.logger.debug(f"Decompressing {self.input_path} to {temp_file.name}")

        with self._open_file() as input_handle:
            with temp_file as output_handle:
                for line in input_handle:
                    output_handle.write(line)

        return Path(temp_file.name)

    def _cleanup_temp_files(self, temp_files: list) -> None:
        """Remove temporary files."""
        for temp_file in temp_files:
            if temp_file and temp_file != self.input_path and isinstance(temp_file, Path):
                temp_file.unlink(missing_ok=True)

    def _save_output(self) -> Path:
        """Save processed features to output directory."""
        self.logger.debug("Saving output file...")

        output_dir = self.output_dir
        if self.settings.output_subdir_name:
            output_dir = output_dir / self.settings.output_subdir_name

        output_dir.mkdir(parents=True, exist_ok=True)

        base_name = self.feature_config.filename
        for suffix in self.input_path.suffixes:
            base_name = base_name.replace(suffix, '')

        if self.settings.output_filename_suffix:
            output_filename = f"{base_name}_{self.settings.output_filename_suffix}.gff"
        else:
            output_filename = f"{base_name}.gff"

        coding = self.settings.coding_type

        if coding in ('gz', 'gzip', CT.GZIP):
            output_filename += '.gz'
        elif coding in ('bz2', 'bzip2', CT.BZIP2):
            output_filename += '.bz2'

        if self.validation_level == 'minimal':
            self.logger.debug("Minimal mode - validating format and coding requirements")

            if self.feature_config.detected_format != FeatureFormat.GFF:
                error_msg = f'Minimal mode requires GFF format, got {self.feature_config.detected_format}. Use validation_level "trust" or "strict" to convert.'
                self.logger.add_validation_issue(
                    level='ERROR',
                    category='feature',
                    message=error_msg,
                    details={'file': self.feature_config.filename, 'detected_format': str(self.feature_config.detected_format)}
                )
                raise FeatureValidationError(error_msg)

            input_coding = self.feature_config.coding_type
            required_coding = coding

            if input_coding is None:
                input_coding = CT.NONE
            if required_coding is None:
                required_coding = CT.NONE

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

            output_path_minimal = output_dir / output_filename
            self.logger.debug(f"Copying {self.input_path} to {output_path_minimal}")
            shutil.copy2(self.input_path, output_path_minimal)
            self.logger.info(f"Output saved: {output_path_minimal}")
            return output_path_minimal

        output_path = output_dir / output_filename
        self.logger.debug(f"Writing output to: {output_path}")

        coding_enum = CT.NONE
        if self.settings.coding_type in ('gz', 'gzip'):
            coding_enum = CT.GZIP
        elif self.settings.coding_type in ('bz2', 'bzip2'):
            coding_enum = CT.BZIP2

        with open_compressed_writer(output_path, coding_enum, threads=self.threads) as handle:
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

        self.logger.info(f"Output saved: {output_path}")
        return output_path