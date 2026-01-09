"""Genome file validator and processor for FASTA and GenBank formats."""

from pathlib import Path
from typing import Optional, List, IO, Union
from dataclasses import dataclass
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import gc_fraction
import shutil

from validation_pkg.logger import get_logger
from validation_pkg.utils.settings import BaseSettings
from validation_pkg.utils.formats import CodingType as CT, GenomeFormat
from validation_pkg.exceptions import (
    ValidationError,
    GenomeValidationError,
    FileFormatError,
    FastaFormatError,
    GenBankFormatError,
    CompressionError
)
from validation_pkg.utils.file_handler import open_compressed_writer


@dataclass
class OutputMetadata(BaseSettings):
    """Metadata returned from genome validation."""
    input_file: str = None
    output_file: str = None
    output_filename: str = None
    num_sequences: int = None
    total_genome_size: int = None  # strict only
    longest_sequence_length: int = None
    longest_sequence_id: str = None
    gc_content: float = None  # strict only
    n50: int = None  # strict only
    plasmid_count: int = None
    plasmid_filenames: List[str] = None
    num_sequences_filtered: int = None
    validation_level: str = None
    elapsed_time: float = None

    # Inter-file validation fields
    sequence_ids: List[str] = None
    sequence_lengths: dict = None

class GenomeValidator:
    """Validates and processes genome files in FASTA and GenBank formats."""

    @dataclass
    class Settings(BaseSettings):
        """Settings for genome validation and processing."""
        # Validation thresholds
        allow_empty_sequences: bool = False
        allow_empty_id: bool = False
        warn_n_sequences: int = 2

        # Editing specifications
        is_plasmid: bool = False
        plasmid_split: bool = False
        plasmids_to_one: bool = False
        main_longest: bool = True
        main_first: bool = False

        replace_id_with: Optional[str] = None
        min_sequence_length: int = 100

        # Output format
        coding_type: Optional[CT] = CT.NONE
        output_filename_suffix: Optional[str] = None
        output_subdir_name: Optional[str] = None

        def __post_init__(self):
            """Validate and normalize settings after initialization."""
            # Normalize coding_type if needed (handles direct instantiation)
            self.coding_type = CT.normalize(self.coding_type)

            if self.plasmid_split and self.plasmids_to_one:
                raise ValueError(
                    "plasmid_split and plasmids_to_one cannot both be True. "
                    "Choose one: plasmid_split creates separate files for each plasmid, "
                    "plasmids_to_one merges all plasmids into a single file."
                )

            if self.main_longest and self.main_first:
                raise ValueError(
                    "main_longest and main_first cannot both be True. "
                    "Choose one: main_longest selects the longest sequence as main, "
                    "main_first selects the first sequence as main."
                )

    def __init__(self, genome_config, settings: Optional[Settings] = None) -> None:
        self.logger = get_logger()

        # From genome global configuration
        self.genome_config = genome_config
        self.output_dir = genome_config.output_dir
        self.validation_level = genome_config.global_options.get("validation_level")   
        self.threads = genome_config.global_options.get("threads") 
        self.input_path = genome_config.filepath

        # From settings
        self.settings = settings if settings is not None else self.Settings() 
        self.output_metadata = OutputMetadata()

        # Parsed data
        self.sequences = []  # List of SeqRecord objects

        if not self.validation_level:
            self.validation_level = 'strict'

        if not self.threads:
            from validation_pkg.config_manager import DEFAULT_THREADS
            self.threads = DEFAULT_THREADS

        # Tracking for metadata
        self.num_sequences_filtered = 0
        self.plasmid_filenames = []

    def _calculate_gc_content(self, sequences: List[SeqRecord]) -> float:
        """Calculate GC content percentage for all sequences using BioPython."""
        if not sequences:
            return 0.0

        total_gc = 0.0
        total_bases = 0

        for seq in sequences:
            seq_length = len(seq.seq)
            if seq_length > 0:
                # gc_fraction returns value between 0 and 1
                total_gc += gc_fraction(seq.seq) * seq_length
                total_bases += seq_length

        if total_bases == 0:
            return 0.0

        # Return as percentage (0-100)
        return (total_gc / total_bases) * 100

    def _calculate_n50(self, sequences: List[SeqRecord]) -> int:
        """Calculate N50 assembly quality metric."""
        if not sequences:
            return 0

        # Get lengths and sort in descending order
        lengths = sorted([len(seq.seq) for seq in sequences], reverse=True)
        total_length = sum(lengths)

        # Find N50
        cumulative_length = 0
        for length in lengths:
            cumulative_length += length
            if cumulative_length >= total_length / 2:
                return length

        return 0

    def _fill_output_metadata(self, output_path: Path) -> None:
        """Populate output metadata with validation results."""
        self.output_metadata.input_file = self.genome_config.filename
        self.output_metadata.validation_level = self.validation_level
        self.output_metadata.output_file = str(output_path) if output_path else None
        self.output_metadata.output_filename = output_path.name if output_path else None

        # Minimal mode - only basic info available
        if self.validation_level == 'minimal':
            return self.output_metadata

        # Trust and Strict modes - sequences were parsed
        if self.sequences:
            self.output_metadata.num_sequences = len(self.sequences)

            # Find longest sequence
            if self.sequences:
                longest_seq = max(self.sequences, key=lambda x: len(x.seq))
                self.output_metadata.longest_sequence_length = len(longest_seq.seq)
                self.output_metadata.longest_sequence_id = str(longest_seq.id)

            # Inter-file validation fields (trust and strict modes)
            self.output_metadata.sequence_ids = [str(seq.id) for seq in self.sequences]
            self.output_metadata.sequence_lengths = {str(seq.id): len(seq.seq) for seq in self.sequences}

        # Plasmid information (available in both trust and strict)
        if self.plasmid_filenames:
            self.output_metadata.plasmid_count = len(self.plasmid_filenames)
            self.output_metadata.plasmid_filenames = self.plasmid_filenames
        elif self.num_sequences_filtered > 0:
            # No plasmids split, but sequences were filtered
            pass

        self.output_metadata.num_sequences_filtered = self.num_sequences_filtered

        # Strict mode only - compute expensive statistics
        if self.validation_level == 'strict' and self.sequences:
            # Total genome size
            self.output_metadata.total_genome_size = sum(len(seq.seq) for seq in self.sequences)

            # GC content
            self.output_metadata.gc_content = self._calculate_gc_content(self.sequences)

            # N50
            self.output_metadata.n50 = self._calculate_n50(self.sequences)

    def run(self) -> OutputMetadata:
        """Execute validation and processing workflow."""
        self.logger.start_timer("genome_validation")
        self.logger.info(f"Processing genome file: {self.genome_config.filename}")
        self.logger.debug(f"Format: {self.genome_config.detected_format}, Compression: {self.genome_config.coding_type}")

        try:
            self._parse_file()

            self._validate_sequences()

            self._apply_edits() # include plasmid handle

            output_path = self._save_output()

            elapsed = self.logger.stop_timer("genome_validation")
            self.logger.info(f"✓ Genome validation completed in {elapsed:.2f}s")

            # Record file timing for report
            self.logger.add_file_timing(
                self.genome_config.filename,
                "genome",
                elapsed
            )

            # Create and return OutputMetadata
            self._fill_output_metadata(output_path)
            self.output_metadata.elapsed_time = elapsed
            return self.output_metadata

        except Exception as e:
            self.logger.error(f"Genome validation failed: {e}")
            raise

    def _open_file(self, mode: str = 'rt') -> IO:
        """Open file with automatic decompression based on genome_config."""
        try:
            from validation_pkg.utils.file_handler import open_file_with_coding_type
            return open_file_with_coding_type(
                self.input_path,
                self.genome_config.coding_type,
                mode
            )
        except CompressionError as e:
            self.logger.add_validation_issue(
                level='ERROR',
                category='genome',
                message=str(e),
                details={'file': str(self.input_path)}
            )
            raise
    
    def _parse_file(self) -> None:
        """Parse file using BioPython and validate format from genome_config."""
        self.logger.debug(f"Parsing {self.genome_config.detected_format} file (validation_level={self.validation_level})...")

        # Minimal mode - skip parsing
        if self.validation_level == 'minimal':
            self.logger.info("Minimal validation mode - skipping file parsing")
            self.sequences = []
            return

        # Trust and Strict modes - parse all sequences
        try:
            with self._open_file() as handle:
                # Parse using BioPython with clean enum conversion
                biopython_format = self.genome_config.detected_format.to_biopython()
                self.sequences = list(SeqIO.parse(handle, biopython_format))

            # Validate sequences
            if not self.sequences:
                error_msg = f"No sequences found in {self.genome_config.detected_format} file"
                self.logger.add_validation_issue(
                    level='ERROR',
                    category='genome',
                    message=error_msg,
                    details={'file': self.genome_config.filename, 'format': str(self.genome_config.detected_format)}
                )
                raise GenomeValidationError(error_msg)

            if self.validation_level == 'trust':
                self.logger.info(f"Trust mode - parsed {len(self.sequences)} sequence(s), will validate first only")
            else:
                self.logger.debug(f"Parsed {len(self.sequences)} sequence(s)")

        except FileFormatError:
            raise
        except Exception as e:
            error_msg = f"Failed to parse {self.genome_config.detected_format} file: {e}"

            # Determine exception type based on format using enum
            from validation_pkg.utils.formats import GenomeFormat
            if self.genome_config.detected_format == GenomeFormat.FASTA:
                exception_class = FastaFormatError
            elif self.genome_config.detected_format == GenomeFormat.GENBANK:
                exception_class = GenBankFormatError
            else:
                exception_class = FileFormatError

            self.logger.add_validation_issue(
                level='ERROR',
                category='genome',
                message=error_msg,
                details={
                    'file': self.genome_config.filename,
                    'format': str(self.genome_config.detected_format),
                    'error': str(e)
                }
            )
            raise exception_class(error_msg) from e
    
    def _validate_sequences(self) -> None:
        """Validate parsed sequences."""
        self.logger.debug("Validating sequences...")

        # Minimal mode - no sequences to validate
        if self.validation_level == 'minimal':
            self.logger.debug("Minimal mode - skipping sequence validation")
            return

        # Trust mode - validate only first sequence
        if self.validation_level == 'trust':
            self.logger.debug("Trust mode - validating first sequence only")
            if len(self.sequences) > 0:
                record = self.sequences[0]
                # Check sequence ID
                if not record.id and not self.settings.allow_empty_id:
                    error_msg = f"First sequence has no ID"
                    self.logger.add_validation_issue(
                        level='ERROR',
                        category='genome',
                        message=error_msg,
                        details={'sequence_index': 0, 'sequence_id': str(record.id)}
                    )
                    raise GenomeValidationError(error_msg)

                # Check sequence length
                if len(record.seq) == 0 and not self.settings.allow_empty_sequences:
                    error_msg = f"First sequence '{record.id}' has zero length"
                    self.logger.add_validation_issue(
                        level='ERROR',
                        category='genome',
                        message=error_msg,
                        details={'sequence_id': record.id, 'index': 0}
                    )
                    raise GenomeValidationError(error_msg)

                self.logger.debug(f"✓ First sequence validated: {record.id} ({len(record.seq)} bp)")
            return

        # Strict mode - full validation
        # Warn about number of sequences
        if len(self.sequences) >= self.settings.warn_n_sequences:
            self.logger.add_validation_issue(
                level='WARNING',
                category='genome',
                message=f"High number of sequences: {len(self.sequences)}",
                details={
                    'num_sequences': len(self.sequences),
                    'threshold': self.settings.warn_n_sequences
                }
            )

        for idx, record in enumerate(self.sequences):
            # Check sequence ID
            if not record.id and not self.settings.allow_empty_id:
                error_msg = f"Sequence at index {idx} has no ID"
                self.logger.add_validation_issue(
                    level='ERROR',
                    category='genome',
                    message=error_msg,
                    details={'sequence_index': idx, 'sequence_id': str(record.id)}
                )
                raise GenomeValidationError(error_msg)

            # Check sequence length
            if len(record.seq) == 0 and not self.settings.allow_empty_sequences:
                error_msg = f"Sequence '{record.id}' has zero length"
                self.logger.add_validation_issue(
                    level='ERROR',
                    category='genome',
                    message=error_msg,
                    details={'sequence_id': record.id, 'index': idx}
                )
                raise GenomeValidationError(error_msg)

        self.logger.debug("✓ Sequence validation passed")
    
    def _apply_edits(self) -> None:
        """Apply editing specifications to sequences based on settings."""
        self.logger.debug("Applying editing specifications from settings...")

        # Minimal mode - skip edits, file will be copied as-is
        if self.validation_level == 'minimal':
            self.logger.debug("Minimal mode - skipping edits")
            return

        # 1. Remove short sequences
        if self.settings.min_sequence_length > 0:
            min_length = self.settings.min_sequence_length
            original_count = len(self.sequences)
            self.sequences = [seq for seq in self.sequences if len(seq.seq) >= min_length]

            if len(self.sequences) < original_count:
                removed = original_count - len(self.sequences)
                self.num_sequences_filtered = removed  # Track for metadata
                self.logger.add_validation_issue(
                    level='WARNING',
                    category='genome',
                    message=f'Removed {removed} sequence(s) shorter than {min_length}bp',
                    details={'min_length': min_length, 'removed_count': removed}
                )

            # If all sequences were filtered out, return early
            if len(self.sequences) == 0:
                self.logger.debug("All sequences filtered out")
                return

        # 2. Replace sequence IDs with auto-increment
        if self.settings.replace_id_with:
            prefix = self.settings.replace_id_with
            for idx, record in enumerate(self.sequences):
                # Store original ID in description
                record.description = f"{record.id}"
                # First sequence: just prefix, subsequent: prefix + number (no separator)
                if idx == 0:
                    record.id = f"{prefix}"
                else:
                    record.id = f"{prefix}{idx}"
            self.logger.debug(f"Replaced sequence IDs with '{prefix}' (auto-incremented for multiple sequences)")

        # 3. Handle plasmid sequences
        if self.settings.is_plasmid:
            # Treat all sequences as plasmids (no main chromosome)
            self._handle_plasmids(self.sequences)
            self.sequences = []
        elif self.settings.plasmid_split or self.settings.plasmids_to_one:
            # Only select main sequence if we're actually doing plasmid handling
            self.main_sequence, plasmid_sequences = self._select_main_sequence(self.sequences)
            self._handle_plasmids(plasmid_sequences)
            self.sequences = [self.main_sequence]
        # else: keep all sequences in main output file (default behavior)

        self.logger.debug(f"✓ Edits applied, {len(self.sequences)} sequence(s) remaining")

    def _select_main_sequence(self, sequences: List[SeqRecord]) -> tuple[SeqRecord, List[SeqRecord]]:
        """Select main chromosome from sequences based on settings."""
        if self.settings.main_longest:
            # Sort by length (longest first)
            sorted_sequences = sorted(sequences, key=lambda x: len(x.seq), reverse=True)
            main_sequence = sorted_sequences[0]
            plasmid_sequences = sorted_sequences[1:]
            self.logger.debug(f"Selected longest sequence as main: {main_sequence.id} ({len(main_sequence.seq)} bp)")
        elif self.settings.main_first:
            # Keep original order, first is main
            main_sequence = sequences[0]
            plasmid_sequences = sequences[1:]
            self.logger.debug(f"Selected first sequence as main: {main_sequence.id} ({len(main_sequence.seq)} bp)")
        else:
            # This should not happen due to __post_init__ validation
            raise ValueError("Either main_longest or main_first must be True")

        return main_sequence, plasmid_sequences

    def _handle_plasmids(self, plasmid_sequences: List[SeqRecord]) -> None:
        """Handle plasmid sequences according to settings."""
        if len(plasmid_sequences) == 0:
            return

        self.logger.debug(f"Handling plasmid file with {len(plasmid_sequences)} sequence(s)...")

        if self.settings.plasmid_split:
            # Save each plasmid to separate file
            self.logger.info(
                f"Processing plasmid file: saving {len(plasmid_sequences)} plasmid(s) "
                f"to separate files"
            )
            for i, plasmid in enumerate(plasmid_sequences):
                # Save each plasmid individually (pass as list to reuse existing method)
                self._save_plasmid_file([plasmid], i)

        elif self.settings.plasmids_to_one:
            # Save all plasmids to one merged file
            self.logger.info(
                f"Processing plasmid file: merging {len(plasmid_sequences)} plasmid(s) "
                f"into one file"
            )
            self._save_plasmid_file(plasmid_sequences, "")

        else:
            # Keep all sequences in main output file
            self.logger.info(
                f"Processing plasmid file: keeping all {len(plasmid_sequences)} "
                f"sequence(s) in main output file"
            )

    def _save_plasmid_file(self, plasmid_sequences: List[SeqRecord], index: str) -> None:
        """Save plasmid sequences to a separate file."""
        from validation_pkg.utils.path_utils import build_safe_output_dir
        from validation_pkg.utils.file_handler import strip_all_extensions

        # Build safe output directory (with path traversal protection)
        # Only use subdir if it's not "plasmid" itself
        subdir = None
        if self.settings.output_subdir_name and self.settings.output_subdir_name != "plasmid":
            subdir = self.settings.output_subdir_name
        output_dir = build_safe_output_dir(self.genome_config.output_dir, subdir, create=True)

        # Generate plasmid filename using utility
        base_name = strip_all_extensions(self.genome_config.filename, self.input_path)

        # Add plasmid suffix and optional user suffix
        if self.settings.output_filename_suffix:
            plasmid_filename = f"{base_name}_{self.settings.output_filename_suffix}_plasmid{index}.fasta"
        else:
            plasmid_filename = f"{base_name}_plasmid{index}.fasta"

        plasmid_filename += self.settings.coding_type.to_extension()

        plasmid_path = output_dir / plasmid_filename

        # Write plasmid sequences with appropriate compression
        self.logger.debug(f"Writing plasmid sequences to: {plasmid_path}")

        # Use optimized compression writer
        with open_compressed_writer(plasmid_path, self.settings.coding_type, threads=self.threads) as handle:
            SeqIO.write(plasmid_sequences, handle, 'fasta')

        self.logger.info(f"Plasmid sequences saved: {plasmid_path}")

        # Track plasmid filename for metadata
        self.plasmid_filenames.append(plasmid_filename)

        # Log details about each plasmid
        for seq in plasmid_sequences:
            self.logger.debug(f"  Plasmid: {seq.id} ({len(seq.seq)} bp)")
    
    def _save_output(self) -> Path:
        """Save processed genome to output directory using settings."""
        from validation_pkg.utils.file_handler import build_output_path
        from validation_pkg.utils.validation_helpers import (
            validate_minimal_mode_requirements,
            copy_file_minimal_mode
        )

        self.logger.debug("Saving output file...")

        # Build output path using utility (with path traversal protection)
        output_path = build_output_path(
            base_dir=self.genome_config.output_dir,
            input_filename=self.genome_config.filename,
            output_format="fasta",
            coding_type=self.settings.coding_type,
            subdir_name=self.settings.output_subdir_name,
            filename_suffix=self.settings.output_filename_suffix,
            input_path=self.input_path
        )

        # Minimal mode - copy file as-is without parsing
        # Required: FASTA format + coding must match settings.coding_type
        if self.validation_level == 'minimal':
            self.logger.debug("Minimal mode - validating format and coding requirements")

            try:
                validate_minimal_mode_requirements(
                    self.genome_config.detected_format,
                    GenomeFormat.FASTA,
                    self.genome_config.coding_type,
                    self.settings.coding_type,
                    self.genome_config.filename,
                    self.logger,
                    'genome'
                )
            except ValidationError as e:
                # Re-raise as GenomeValidationError for consistency
                raise GenomeValidationError(str(e)) from e

            return copy_file_minimal_mode(self.input_path, output_path, self.logger)

        # Strict and Trust modes - write sequences using BioPython
        # Check if we have sequences to write
        if self.sequences == []:
            self.logger.warning("No sequences to write")
            return None

        # Write output with appropriate compression
        self.logger.debug(f"Writing output to: {output_path}")

        # Use optimized compression writer
        with open_compressed_writer(output_path, self.settings.coding_type, threads=self.threads) as handle:
            SeqIO.write(self.sequences, handle, 'fasta')

        self.logger.info(f"Output saved: {output_path}")

        return output_path
