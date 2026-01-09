"""Read file validator and processor for FASTQ and BAM formats."""

from pathlib import Path
from typing import Optional, List, IO, Union, Dict, Any
from dataclasses import dataclass, asdict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import shutil
from multiprocessing import Pool
from functools import partial
import re


from validation_pkg.logger import get_logger
from validation_pkg.utils.settings import BaseSettings
from validation_pkg.exceptions import (
    ValidationError,
    ReadValidationError,
    FileFormatError,
    FastqFormatError,
    BamFormatError,
    CompressionError
)
from validation_pkg.utils.formats import ReadFormat
from validation_pkg.utils.formats import CodingType as CT

# Constants
TRUST_MODE_SAMPLE_SIZE = 10  # Number of records to validate in trust mode
PARALLEL_CHUNK_MULTIPLIER = 4  # Multiplier for calculating chunk size (sequences // (threads * 4))
MIN_CHUNK_SIZE = 1000  # Minimum chunk size for parallel processing

# Illumina paired-end filename patterns (compiled regex for performance)
# Pattern format: (regex, read_number_group_index)
ILLUMINA_PAIRED_END_PATTERNS = [
    # With lane numbers: _R1_001, _R2_001
    (re.compile(r'(.+)_(R[12])_\d+$'), 2),
    # With lane numbers: _1_001, _2_001
    (re.compile(r'(.+)_([12])_\d+$'), 2),
    # With arbitrary suffix: _R1_combined, _R2_final
    (re.compile(r'(.+)_(R[12])_[A-Za-z0-9_]+$'), 2),
    # With arbitrary suffix (no R prefix): _1_combined, _2_final
    (re.compile(r'(.+)_([12])_[A-Za-z0-9_]+$'), 2),
    # No separator with suffix: sampleR1_combined
    (re.compile(r'(.+)(R[12])_[A-Za-z0-9_]+$'), 2),
    # Standard Illumina: _R1, _R2
    (re.compile(r'(.+)_(R[12])$'), 2),
    # SRA style: _1, _2
    (re.compile(r'(.+)_([12])$'), 2),
    # Dot separator: .R1, .R2
    (re.compile(r'(.+)\.(R[12])$'), 2),
    # Dot with numbers: .1, .2
    (re.compile(r'(.+)\.([12])$'), 2),
    # No separator (must be at end): sampleR1, sampleR2
    (re.compile(r'(.+)(R[12])$'), 2),
]

from validation_pkg.utils.file_handler import (
    bz2_to_gz,
    none_to_gz,
    gz_to_bz2,
    gz_to_none,
    bz2_to_none,
    none_to_bz2,
    open_compressed_writer
)


# Standalone function for parallel processing (must be picklable - defined at module level)
def _validate_single_read(record, check_invalid_chars: bool, allow_empty_id: bool):
    """Validate a single read record (parallelizable)."""
    # Check empty ID
    if not record.id and not allow_empty_id:
        return {
            'error': f"Sequence has no ID",
            'record_id': None,
            'details': {'sequence_index': None}  # Index will be added by caller
        }

    # Check invalid characters
    if check_invalid_chars:
        valid_chars = set('ATCGNatcgn')
        seq_str = str(record.seq)
        invalid_chars = set(seq_str) - valid_chars

        if invalid_chars:
            char_count = sum(1 for c in seq_str if c in invalid_chars)
            return {
                'error': f"Contains {char_count} invalid character(s): {', '.join(sorted(invalid_chars))}",
                'record_id': record.id,
                'details': {
                    'invalid_chars': list(invalid_chars),
                    'count': char_count
                }
            }

    return {'success': True, 'record_id': record.id}


@dataclass
class OutputMetadata(BaseSettings):
    """Metadata returned from read validation."""
    input_file: str = None
    output_file: str = None
    output_filename: str = None
    base_name: str = None
    read_number: int = None
    ngs_type_detected: str = None
    num_reads: int = None
    validation_level: str = None
    elapsed_time: float = None

    # Strict mode statistics
    n50: int = None
    total_bases: int = None
    mean_read_length: float = None
    longest_read_length: int = None
    shortest_read_length: int = None

    def __str__(self):
        parts = [f"Validation Level: {self.validation_level or 'N/A'}"]
        parts.append(f"Output File: {self.output_file or 'N/A'}")
        parts.append(f"Output Filename: {self.output_filename or 'N/A'}")

        # Add Illumina pattern info if available
        if self.base_name or self.read_number or self.ngs_type_detected or self.num_reads is not None:
            parts.append("\n--- Pattern & Read Info ---")
            parts.append(f"Base Name: {self.base_name or 'N/A'}")
            parts.append(f"Read Number: {self.read_number if self.read_number is not None else 'N/A'}")
            parts.append(f"NGS Type Detected: {self.ngs_type_detected or 'N/A'}")
            parts.append(f"Number of Reads: {self.num_reads if self.num_reads is not None else 'N/A'}")

        # Add strict statistics if present
        if any(v is not None for v in [self.n50, self.total_bases, self.mean_read_length,
                                       self.longest_read_length, self.shortest_read_length]):
            parts.append("\n--- Strict Statistics ---")
            parts.append(f"N50: {self.n50 if self.n50 is not None else 'N/A'} bp")
            parts.append(f"Total Bases: {self.total_bases if self.total_bases is not None else 'N/A'} bp")
            parts.append(f"Mean Read Length: {self.mean_read_length if self.mean_read_length is not None else 'N/A'} bp")
            parts.append(f"Longest Read Length: {self.longest_read_length if self.longest_read_length is not None else 'N/A'} bp")
            parts.append(f"Shortest Read Length: {self.shortest_read_length if self.shortest_read_length is not None else 'N/A'} bp")

        return "\n".join(parts)


class ReadValidator:
    """Validates and processes sequencing read files in FASTQ and BAM formats."""

    @dataclass
    class Settings(BaseSettings):
        """Settings for read validation and processing."""
        # Validation options
        check_invalid_chars: bool = False
        allow_empty_id: bool = False
        allow_duplicate_ids: bool = True

        # Editing specifications
        keep_bam: bool = True
        ignore_bam: bool = True

        # Output format
        coding_type: Optional[CT] = CT.NONE
        output_filename_suffix: Optional[str] = None
        output_subdir_name: Optional[str] = None
        outdir_by_ngs_type: bool = False

        def __post_init__(self):
            """Normalize settings after initialization."""
            # Normalize coding_type to handle strings and None
            self.coding_type = CT.normalize(self.coding_type)

    def __init__(self, read_config, settings: Optional[Settings] = None) -> None:
        self.logger = get_logger()

        # From genome global configuration
        self.read_config = read_config
        self.output_dir = read_config.output_dir
        self.validation_level = read_config.global_options.get("validation_level")   
        self.threads = read_config.global_options.get("threads") 
        self.input_path = read_config.filepath

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

        # Apply outdir_by_ngs_type if enabled
        if self.settings.outdir_by_ngs_type:
            # Override output_subdir_name with ngs_type
            self.settings = self.settings.update(output_subdir_name=self.read_config.ngs_type)
            self.logger.debug(f"Applied outdir_by_ngs_type: output_subdir_name set to '{self.read_config.ngs_type}'")

    def _calculate_read_statistics(self) -> dict:
        """Calculate read statistics for strict mode."""
        if not self.sequences:
            return {
                'n50': 0,
                'total_bases': 0,
                'mean_read_length': 0.0,
                'longest_read_length': 0,
                'shortest_read_length': 0
            }

        # Get read lengths
        lengths = [len(seq.seq) for seq in self.sequences]

        # Calculate N50
        sorted_lengths = sorted(lengths, reverse=True)
        total_length = sum(sorted_lengths)
        cumulative_length = 0
        n50 = 0
        for length in sorted_lengths:
            cumulative_length += length
            if cumulative_length >= total_length / 2:
                n50 = length
                break

        return {
            'n50': n50,
            'total_bases': total_length,
            'mean_read_length': total_length / len(lengths) if lengths else 0.0,
            'longest_read_length': max(lengths) if lengths else 0,
            'shortest_read_length': min(lengths) if lengths else 0
        }

    def _fill_output_metadata(self, output_path: Path) -> None:
        """Populate output metadata with validation results."""
        # Populate and return OutputMetadata
        self.output_metadata.input_file = self.read_config.filename
        self.output_metadata.output_file = str(output_path) if output_path else None
        self.output_metadata.output_filename = output_path.name if output_path else None
        self.output_metadata.num_reads = len(self.sequences)
        self.output_metadata.validation_level = self.validation_level

        # Calculate read statistics in strict mode only
        if self.validation_level == 'strict' and self.sequences:
            stats = self._calculate_read_statistics()
            self.output_metadata.n50 = stats['n50']
            self.output_metadata.total_bases = stats['total_bases']
            self.output_metadata.mean_read_length = stats['mean_read_length']
            self.output_metadata.longest_read_length = stats['longest_read_length']
            self.output_metadata.shortest_read_length = stats['shortest_read_length']

    def run(self) -> OutputMetadata:
        """Execute validation and processing workflow."""
        self.logger.start_timer("read_validation")
        self.logger.info(f"Processing read file: {self.read_config.filename}")
        self.logger.debug(f"Format: {self.read_config.detected_format}, Compression: {self.read_config.coding_type}")

        try:
            # Special workflow for BAM files
            if self.read_config.detected_format == ReadFormat.BAM:
                if self.settings.ignore_bam:
                    self.logger.warning(f"Cannot process BAM files, the file will be ignored")
                    # Return minimal metadata indicating the file was skipped
                    return OutputMetadata(
                        input_file=str(self.input_path),
                        output_file=None,
                        output_filename=None,
                        num_reads=0,
                        validation_level=self.settings.validation_level,
                        elapsed_time=0.0
                    )

                # Step 1: Copy original BAM to output (if keep_bam is enabled)
                if self.settings.keep_bam:
                    self._copy_bam_to_output()

                # Step 2: Convert BAM to FASTQ (populates self.sequences)
                self._convert_bam_to_fastq()

                # Step 3: Validate sequences
                self._validate_sequences()

                # Step 4: Apply editing specifications
                self._apply_edits()

                # Step 6: Save FASTQ output
                output_path = self._save_output()

            else:
                # Standard FASTQ workflow
                # Step 1: Parse and validate
                self._parse_file()

                # Validate sequences
                self._validate_sequences()

                # Step 2: Apply editing specifications
                self._apply_edits()

                # Step 3: If ilumina, detect file read number and basename
                if self.read_config.ngs_type == 'illumina':
                    self._detect_illumina_pattern(self.read_config.filename)

                # Step 4: Save to output directory (already FASTQ)
                output_path = self._save_output()

            elapsed = self.logger.stop_timer("read_validation")
            self.logger.info(f"✓ Read validation completed in {elapsed:.2f}s")

            # Record file timing for report
            self.logger.add_file_timing(
                self.read_config.filename,
                "read",
                elapsed
            )

            # Create and return OutputMetadata
            self._fill_output_metadata(output_path)
            self.output_metadata.elapsed_time = elapsed
            return self.output_metadata

        except Exception as e:
            self.logger.error(f"Read validation failed: {e}")
            raise

    def _open_file(self, mode: str = 'rt') -> IO:
        """Open file with automatic decompression based on read_config."""
        try:
            from validation_pkg.utils.file_handler import open_file_with_coding_type
            return open_file_with_coding_type(
                self.input_path,
                self.read_config.coding_type,
                mode
            )
        except CompressionError as e:
            self.logger.add_validation_issue(
                level='ERROR',
                category='read',
                message=str(e),
                details={'file': str(self.input_path)}
            )
            raise

    def _count_lines_fast(self) -> int:
        """Fast line counting using Python file iteration."""
        from validation_pkg.utils.file_handler import open_file_with_coding_type

        try:
            with open_file_with_coding_type(self.input_path, self.read_config.coding_type, mode='rt') as f:
                return sum(1 for _ in f)
        except Exception as e:
            raise ReadValidationError(f"Line counting failed: {e}")

    def _parse_file(self) -> None:
        """Parse FASTQ file using BioPython and validate format from read_config."""
        self.logger.debug(f"Parsing {self.read_config.detected_format} file (validation_level={self.validation_level})...")

        # Minimal mode - skip all parsing and validation except CodingType
        if self.validation_level == 'minimal':
            self.logger.info("Minimal validation mode - skipping file parsing")
            # Set empty sequences list - file will be copied as-is in _save_output
            self.sequences = []
            return

        # Trust mode: Parse only first TRUST_MODE_SAMPLE_SIZE sequences for validation
        # Strict mode: Parse all sequences
        if self.validation_level == 'trust':
            self.logger.info(f"Trust mode - parsing first {TRUST_MODE_SAMPLE_SIZE} sequences only for validation")
            parse_limit = TRUST_MODE_SAMPLE_SIZE
        else:
            # Strict mode - parse all
            parse_limit = None

        try:
            # Get total line count for progress reporting (strict mode only)
            if self.validation_level == 'strict':
                try:
                    line_count = self._count_lines_fast()
                    estimated_sequences = line_count // 4  # FASTQ has 4 lines per sequence
                    self.logger.info(f"Processing {estimated_sequences:,} reads from FASTQ file...")
                    progress_interval = max(1, estimated_sequences // 20)  # Report every 5%
                except Exception as e:
                    # If counting fails, just proceed without progress reporting
                    self.logger.debug(f"Could not count lines for progress reporting: {e}")
                    estimated_sequences = None
                    progress_interval = 100000  # Report every 100k reads
            else:
                estimated_sequences = None
                progress_interval = None

            # Open file with automatic decompression
            handle = self._open_file()
            try:
                # Parse using BioPython with clean enum conversion
                biopython_format = self.read_config.detected_format.to_biopython()

                # Parse sequences
                self.sequences = []
                processed = 0
                for record in SeqIO.parse(handle, biopython_format):
                    self.sequences.append(record)
                    processed += 1

                    # Trust mode - stop after parsing limit
                    if parse_limit and processed >= parse_limit:
                        break

                    # Progress reporting (strict mode only)
                    if progress_interval and processed % progress_interval == 0:
                        if estimated_sequences:
                            percent = (processed / estimated_sequences) * 100
                            self.logger.info(f"Progress: {processed:,}/{estimated_sequences:,} reads ({percent:.1f}%)")
                        else:
                            self.logger.info(f"Progress: {processed:,} reads processed...")
            finally:
                handle.close()

            # Validate we got sequences
            if not self.sequences:
                error_msg = f"No sequences found in {self.read_config.detected_format} file"
                self.logger.add_validation_issue(
                    level='ERROR',
                    category='read',
                    message=error_msg,
                    details={'file': self.read_config.filename, 'format': self.read_config.detected_format}
                )
                raise ReadValidationError(error_msg)

            self.logger.info(f"✓ Parsed {len(self.sequences):,} sequence(s)")

        except FileFormatError:
            raise
        except Exception as e:
            error_msg = f"Failed to parse {self.read_config.detected_format} file: {e}"

            if self.read_config.detected_format == ReadFormat.FASTQ:
                exception_class = FastqFormatError
            elif self.read_config.detected_format == ReadFormat.BAM :
                exception_class = BamFormatError
            else:
                exception_class = FileFormatError

            self.logger.add_validation_issue(
                level='ERROR',
                category='read',
                message=error_msg,
                details={
                    'file': self.read_config.filename,
                    'format': self.read_config.detected_format,
                    'error': str(e)
                }
            )
            raise exception_class(error_msg) from e
    
    def _validate_sequences(self) -> None:
        """Validate parsed sequences with optional parallelization."""
        # Determine how many sequences to validate
        if self.validation_level == 'trust':
            validate_count = min(10, len(self.sequences))
            self.logger.info(f"Trust mode - validating first {validate_count} of {len(self.sequences):,} sequences (sequential)")
            use_parallel = False
        else:
            validate_count = len(self.sequences)
            use_parallel = (self.threads and self.threads > 1)

            if use_parallel:
                self.logger.info(
                    f"Parallel validation enabled: {validate_count:,} sequences across {self.threads} workers",
                    parallel_mode=True,
                    workers=self.threads,
                    total_sequences=validate_count
                )
            else:
                self.logger.info(f"Sequential validation: {validate_count:,} sequences (threads=1)")

        # Parallel validation (strict mode with threads > 1)
        if use_parallel:
            # Determine chunk size: MIN_CHUNK_SIZE records per chunk, or split into PARALLEL_CHUNK_MULTIPLIER chunks per worker
            chunk_size = max(MIN_CHUNK_SIZE, validate_count // (self.threads * PARALLEL_CHUNK_MULTIPLIER))

            self.logger.debug(
                f"Parallel processing configuration: {self.threads} workers, chunk_size={chunk_size:,}",
                workers=self.threads,
                chunk_size=chunk_size,
                estimated_chunks=validate_count // chunk_size
            )

            # Create partial function with settings
            validator_func = partial(
                _validate_single_read,
                check_invalid_chars=self.settings.check_invalid_chars,
                allow_empty_id=self.settings.allow_empty_id
            )

            # Validate in parallel
            self.logger.debug(f"Starting parallel validation across {self.threads} worker processes...")
            with Pool(processes=self.threads) as pool:
                results = pool.map(validator_func, self.sequences[:validate_count], chunksize=chunk_size)

            self.logger.debug(
                f"Parallel validation completed: {validate_count:,} sequences processed",
                sequences_validated=validate_count
            )

            # Collect errors
            errors = []
            for idx, result in enumerate(results):
                if 'error' in result:
                    # Add index to details if missing
                    if result['details'].get('sequence_index') is None:
                        result['details']['sequence_index'] = idx

                    # Log validation issue
                    self.logger.add_validation_issue(
                        level='ERROR',
                        category='read',
                        message=f"Sequence '{result['record_id']}': {result['error']}",
                        details=result['details']
                    )
                    errors.append(result)

            # Raise if any errors found
            if errors:
                error_msg = f"{len(errors)} validation error(s) found in sequences (parallel validation)"
                self.logger.error(error_msg, error_count=len(errors), total_sequences=validate_count)
                raise ReadValidationError(error_msg)
            else:
                self.logger.info(f"✓ All {validate_count:,} sequences validated successfully (parallel mode)")

        # Sequential validation (trust mode or single-threaded)
        else:
            for idx in range(validate_count):
                record = self.sequences[idx]

                # Use the same validation function as parallel mode
                result = _validate_single_read(
                    record,
                    check_invalid_chars=self.settings.check_invalid_chars,
                    allow_empty_id=self.settings.allow_empty_id
                )

                if 'error' in result:
                    error_msg = f"Sequence '{result['record_id']}': {result['error']}"
                    self.logger.add_validation_issue(
                        level='ERROR',
                        category='read',
                        message=error_msg,
                        details={**result['details'], 'sequence_index': idx}
                    )
                    raise ReadValidationError(error_msg)

        # Check for duplicate IDs (always sequential - requires full list)
        if not self.settings.allow_duplicate_ids:
            from collections import Counter
            seq_ids = [record.id for record in self.sequences]
            if len(seq_ids) != len(set(seq_ids)):
                # Use Counter for O(n) instead of O(n²) duplicate detection
                id_counts = Counter(seq_ids)
                duplicates = [sid for sid, count in id_counts.items() if count > 1]
                error_msg = f"Duplicate sequence IDs not allowed: {duplicates}"
                self.logger.add_validation_issue(
                    level='ERROR',
                    category='read',
                    message=error_msg,
                    details={'duplicate_ids': duplicates}
                )
                raise ReadValidationError(error_msg)

        # Final success message (if not already logged in parallel mode)
        if not use_parallel:
            self.logger.debug("✓ Sequence validation passed")
    
    def _apply_edits(self) -> None:
        """Apply editing specifications to sequences based on settings."""
        self.logger.debug("Applying editing specifications from settings...")

        # Currently no edits are implemented
        # This is a placeholder for future functionality

        self.logger.debug(f"✓ Edits applied, {len(self.sequences)} sequence(s) remaining")
    
    def _convert_bam_to_fastq(self) -> None:
        """Convert BAM file to FASTQ format using pysam or samtools."""
        import subprocess

        # Minimal mode - skip conversion
        if self.validation_level == 'minimal':
            self.logger.info("Minimal validation mode - skipping BAM conversion")
            self.sequences = []
            return

        # Trust mode - validate header and first few reads only
        if self.validation_level == 'trust':
            self.logger.info("Trust mode - validating BAM header and first records only")
            try:
                import pysam  # type: ignore

                with pysam.AlignmentFile(str(self.input_path), "rb") as bam_file:
                    # Validate BAM header exists
                    if not bam_file.header:
                        raise BamFormatError("BAM file has no header")

                    self.logger.debug(f"✓ BAM header validated: {len(bam_file.header.get('SQ', []))} sequences in reference")

                    # Read and validate first few records (up to 10)
                    first_records = []
                    for i, read in enumerate(bam_file):
                        if i >= 10:
                            break

                        # Skip unmapped or secondary/supplementary alignments
                        if read.is_unmapped or read.is_secondary or read.is_supplementary:
                            continue

                        # Create SeqRecord from BAM read
                        seq = Seq(read.query_sequence if read.query_sequence else "")
                        qualities = read.query_qualities if read.query_qualities else []

                        record = SeqRecord(
                            seq,
                            id=read.query_name,
                            description="",
                            letter_annotations={"phred_quality": list(qualities)} if qualities else {}
                        )
                        first_records.append(record)

                    if not first_records:
                        raise BamFormatError("No valid reads found in BAM file")

                    self.logger.info(f"✓ Trust mode BAM validation passed - validated {len(first_records)} records")
                    self.sequences = first_records
                    return

            except ImportError:
                # If pysam not available, try samtools for header validation
                try:
                    result = subprocess.run(
                        ["samtools", "view", "-H", str(self.input_path)],
                        capture_output=True,
                        text=True,
                        timeout=30
                    )

                    if result.returncode != 0:
                        raise BamFormatError(f"Failed to read BAM header: {result.stderr}")

                    if not result.stdout.strip():
                        raise BamFormatError("BAM file has no header")

                    self.logger.info("✓ Trust mode BAM validation passed (header only, pysam not available)")
                    self.sequences = []
                    return

                except (FileNotFoundError, subprocess.TimeoutExpired) as e:
                    raise ReadValidationError(
                        "Trust mode BAM validation requires either 'pysam' or 'samtools'. "
                        "Please install one of them."
                    ) from e

        # Strict mode - full conversion (original behavior)
        # Try using pysam first (optional dependency)
        try:
            import pysam  # type: ignore

            self.logger.debug("Using pysam for BAM to FASTQ conversion")
            self.sequences = []

            with pysam.AlignmentFile(str(self.input_path), "rb") as bam_file:
                # Try to get total read count for progress reporting
                try:
                    total_reads = bam_file.count(until_eof=True)
                    bam_file.reset()  # Reset file pointer after counting
                    self.logger.info(f"Processing {total_reads:,} reads from BAM file...")
                    progress_interval = max(1, total_reads // 20)  # Report every 5%
                except Exception:
                    # If counting fails (e.g., truncated file, permission issues),
                    # just proceed without progress reporting
                    total_reads = None
                    progress_interval = 100000  # Report every 100k reads

                processed = 0
                for read in bam_file:
                    # Skip unmapped or secondary/supplementary alignments
                    if read.is_unmapped or read.is_secondary or read.is_supplementary:
                        continue

                    # Create SeqRecord from BAM read
                    seq = Seq(read.query_sequence if read.query_sequence else "")
                    qualities = read.query_qualities if read.query_qualities else []

                    record = SeqRecord(
                        seq,
                        id=read.query_name,
                        description="",
                        letter_annotations={"phred_quality": list(qualities)} if qualities else {}
                    )
                    self.sequences.append(record)

                    processed += 1
                    # Progress reporting
                    if processed % progress_interval == 0:
                        if total_reads:
                            percent = (processed / total_reads) * 100
                            self.logger.info(f"Progress: {processed:,}/{total_reads:,} reads ({percent:.1f}%)")
                        else:
                            self.logger.info(f"Progress: {processed:,} reads processed...")

            self.logger.info(f"✓ Converted {len(self.sequences):,} reads from BAM to FASTQ")
            return

        except ImportError:
            self.logger.debug("pysam not available, trying samtools...")

        # Fall back to samtools
        try:
            # Check if samtools is available
            result = subprocess.run(
                ["samtools", "--version"],
                capture_output=True,
                text=True,
                timeout=5
            )

            if result.returncode != 0:
                raise FileNotFoundError("samtools not found")

            self.logger.debug("Using samtools for BAM to FASTQ conversion")

            # Use samtools fastq to convert BAM to FASTQ
            # samtools fastq writes to stdout by default
            result = subprocess.run(
                ["samtools", "fastq", str(self.input_path)],
                capture_output=True,
                text=True,
                timeout=300  # 5 minutes timeout
            )

            if result.returncode != 0:
                error_msg = f"samtools conversion failed: {result.stderr}"
                self.logger.add_validation_issue(
                    level='ERROR',
                    category='read',
                    message=error_msg,
                    details={'file': self.read_config.filename, 'error': result.stderr}
                )
                raise ReadValidationError(error_msg)

            # Parse FASTQ output from samtools
            from io import StringIO
            fastq_data = StringIO(result.stdout)
            self.sequences = list(SeqIO.parse(fastq_data, "fastq"))

            self.logger.info(f"Converted {len(self.sequences)} reads from BAM to FASTQ using samtools")

        except (FileNotFoundError, subprocess.TimeoutExpired, subprocess.SubprocessError) as e:
            error_msg = (
                "BAM to FASTQ conversion requires either 'pysam' Python package or 'samtools' command-line tool. "
                "Please install one of them:\n"
                "  - pip install pysam\n"
                "  - or install samtools (https://www.htslib.org/)"
            )
            self.logger.add_validation_issue(
                level='ERROR',
                category='read',
                message=error_msg,
                details={'error': str(e)}
            )
            raise ReadValidationError(error_msg) from e
    
    def _copy_bam_to_output(self) -> Path:
        """Copy BAM file to output directory without processing."""
        from validation_pkg.utils.path_utils import build_safe_output_dir
        from validation_pkg.utils.file_handler import strip_all_extensions

        self.logger.debug("Copying BAM file to output directory...")

        # Build safe output directory (with path traversal protection)
        output_dir = build_safe_output_dir(
            self.output_dir,
            self.settings.output_subdir_name,
            create=True
        )

        # Generate output filename - keep original name and compression
        output_filename = self.read_config.filename
        if self.settings.output_filename_suffix:
            # Insert suffix before extensions
            base_name = strip_all_extensions(self.read_config.filename, self.input_path)
            # Reconstruct with suffix and original extensions
            extensions = ''.join(self.input_path.suffixes)
            output_filename = f"{base_name}_{self.settings.output_filename_suffix}{extensions}"

        output_path = output_dir / output_filename

        # Copy file
        self.logger.debug(f"Copying {self.input_path} to {output_path}")
        shutil.copy2(self.input_path, output_path)

        self.logger.info(f"BAM file copied: {output_path}")
        return output_path

    def _save_output(self) -> Path:
        """Save processed read to output directory using settings."""
        from validation_pkg.utils.path_utils import build_safe_output_dir
        from validation_pkg.utils.validation_helpers import (
            validate_minimal_mode_requirements,
            copy_file_minimal_mode
        )

        self.logger.debug("Saving output file...")

        # Build safe output directory (with path traversal protection)
        output_dir = build_safe_output_dir(
            self.output_dir,
            self.settings.output_subdir_name,
            create=True
        )

        # Generate output filename from original filename (without compression extension)
        # Get base name without any extensions
        if self.output_metadata.base_name and self.output_metadata.read_number:
            # Paired-end detected: use detected base name with R1/R2 suffix
            base_name = f"{self.output_metadata.base_name}_R{self.output_metadata.read_number}"
        else:
            # Single-end (no pattern detected): use basename with _R1 suffix
            # This ensures all Illumina reads have consistent _R1 or _R2 naming
            base_name = f"{self.read_config.basename}_R1"

        # Add suffix if specified
        if self.settings.output_filename_suffix:
            output_filename = f"{base_name}_{self.settings.output_filename_suffix}.fastq"
        else:
            output_filename = f"{base_name}.fastq"

        output_filename += self.settings.coding_type.to_extension()
        output_path = output_dir / output_filename

        # Store input and required coding for later use
        input_coding = self.read_config.coding_type
        required_coding = self.settings.coding_type

        # Minimal mode - copy file as-is without parsing
        # Required: FASTQ format + coding must match settings.coding_type
        if self.validation_level == 'minimal':
            self.logger.debug("Minimal mode - validating format and coding requirements")

            # Validate format and coding requirements
            try:
                validate_minimal_mode_requirements(
                    self.read_config.detected_format,
                    ReadFormat.FASTQ,
                    input_coding,
                    required_coding,
                    self.read_config.filename,
                    self.logger,
                    'read'
                )
            except ValidationError as e:
                # Re-raise as ReadValidationError for consistency
                raise ReadValidationError(str(e)) from e

            # Additional FASTQ-specific validation: check line count divisible by 4
            try:
                line_count = self._count_lines_fast()
                if line_count % 4 != 0:
                    error_msg = f"Invalid FASTQ: {line_count:,} lines not divisible by 4"
                    self.logger.add_validation_issue(
                        level='ERROR',
                        category='read',
                        message=error_msg,
                        details={'file': self.read_config.filename, 'line_count': line_count}
                    )
                    raise FastqFormatError(error_msg)
                num_sequences = line_count // 4
                self.logger.info(f"Minimal mode - verified {num_sequences:,} sequences (line count check)")
            except FastqFormatError:
                raise
            except Exception as e:
                error_msg = f"Minimal mode validation failed: {e}"
                self.logger.add_validation_issue(
                    level='ERROR',
                    category='read',
                    message=error_msg,
                    details={'file': self.read_config.filename, 'error': str(e)}
                )
                raise ReadValidationError(error_msg) from e

            return copy_file_minimal_mode(self.input_path, output_path, self.logger)

        # Trust mode - copy original file with coding conversion
        # Trust mode parsed only first 10 sequences for validation, so we copy the original file
        if self.validation_level == 'trust':
            self.logger.debug("Trust mode - copying original file with coding conversion")

            # If input and output coding match, simple copy
            if input_coding == required_coding:
                self.logger.debug(f"Copying {self.input_path} to {output_path} (same coding)")
                shutil.copy2(self.input_path, output_path)
            else:
                # Convert coding using file_handler utilities
                self.logger.debug(f"Converting {input_coding} -> {required_coding}")

                # Map (input, output) pairs to conversion functions
                conversion_map = {
                    (CT.BZIP2, CT.GZIP): bz2_to_gz,
                    (CT.NONE, CT.GZIP): none_to_gz,
                    (CT.GZIP, CT.NONE): gz_to_none,
                    (CT.BZIP2, CT.NONE): bz2_to_none,
                    (CT.NONE, CT.BZIP2): none_to_bz2,
                    (CT.GZIP, CT.BZIP2): gz_to_bz2,
                }

                conversion_func = conversion_map.get((input_coding, required_coding))
                if conversion_func:
                    conversion_func(self.input_path, output_path, threads=self.threads)
                else:
                    error_msg = f"Unsupported coding conversion: {input_coding} -> {required_coding}"
                    raise ReadValidationError(error_msg)

            self.logger.info(f"Output saved: {output_path}")
            return output_path

        # Strict mode - optimize based on whether edits were applied
        # If input is FASTQ, coding matches, and no edits were made, use simple copy
        # Otherwise, write sequences using BioPython

        input_format = self.read_config.detected_format

        # Check if we can use optimized copy (no edits, same format and coding)
        # Note: Currently _apply_edits() does nothing, so no edits are ever made
        # In the future, check if edits were actually applied (e.g., sequences removed/modified)
        can_optimize = (
            input_format == ReadFormat.FASTQ and  # Input is FASTQ
            input_coding == required_coding  # Compression matches
            # Future: add check for "no edits were actually applied"
        )

        if can_optimize:
            # Optimized path: simple copy (no need to parse/write)
            self.logger.debug(f"Optimized copy: {self.input_path} to {output_path} (no edits, same format/coding)")
            shutil.copy2(self.input_path, output_path)
        else:
            # Standard path: write sequences using BioPython
            # Used when: BAM conversion, compression change, or future edits
            self.logger.debug(f"Writing output to: {output_path}")

            # Use optimized compression writer
            with open_compressed_writer(output_path, self.settings.coding_type, threads=self.threads) as handle:
                SeqIO.write(self.sequences, handle, 'fastq')

        self.logger.info(f"Output saved: {output_path}")

        return output_path
    
    def _detect_illumina_pattern(self, filename: str) -> None:
        """Detect Illumina paired-end naming patterns in filename."""

        # Strip only the last 2 extensions (format + compression)
        base = self.read_config.basename

        # Use pre-compiled patterns from module-level constant
        for pattern, read_group_idx in ILLUMINA_PAIRED_END_PATTERNS:
            match = pattern.match(base)
            if match:
                base_name = match.group(1)
                read_indicator = match.group(read_group_idx)

                # Extract read number (1 or 2)
                if read_indicator in ('R1', '1'):
                    read_number = 1
                elif read_indicator in ('R2', '2'):
                    read_number = 2
                else:
                    continue  # Should not happen with our patterns

                self.output_metadata.base_name = base_name
                self.output_metadata.read_number = read_number
                self.output_metadata.ngs_type_detected = 'illumina'

                if hasattr(self, 'logger') and self.logger:
                    self.logger.debug(
                        f"Detected Illumina paired-end pattern: {base_name}_R{read_number}"
                    )
                return  # Early return after successful match

        # No paired-end pattern matched - fallback for single-end or non-standard naming
        # Always set ngs_type_detected for Illumina files, even without pairing info
        self.output_metadata.base_name = base
        self.output_metadata.read_number = 1  # Default to R1 for single-end
        self.output_metadata.ngs_type_detected = 'illumina'

        if hasattr(self, 'logger') and self.logger:
            self.logger.debug(
                f"No paired-end pattern detected in '{base}' - treating as single-end R1"
            )
