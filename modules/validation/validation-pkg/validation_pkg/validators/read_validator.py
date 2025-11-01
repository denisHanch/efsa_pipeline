"""
Read file validator and processor.

This module provides comprehensive validation and processing for sequencing read files in FASTQ
and BAM formats. It handles decompression, parsing, validation, and various quality control
operations for next-generation sequencing data.

Key Features:
    - Multi-format support: FASTQ and BAM (with optional conversion)
    - Compression handling: gzip, bzip2, and uncompressed files
    - Three validation levels: strict (full validation), trust (fast validation), minimal (copy only)
    - BAM to FASTQ conversion using pysam or samtools
    - Fast line counting for trust mode validation
    - NGS type-based output directory organization
    - Parallel compression support: automatic detection of pigz/pbzip2

Classes:
    ReadValidator: Main validator class for sequencing read files
    ReadValidator.Settings: Configuration dataclass for validation behavior
"""

from pathlib import Path
from typing import Optional, List, IO, Union
from dataclasses import dataclass
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import shutil
from multiprocessing import Pool
from functools import partial

from validation_pkg.logger import get_logger
from validation_pkg.utils.settings import BaseSettings
from validation_pkg.exceptions import (
    ReadValidationError,
    FileFormatError,
    FastqFormatError,
    BamFormatError,
    CompressionError
)
from validation_pkg.utils.formats import ReadFormat
from validation_pkg.utils.formats import CodingType as CT
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
    """
    Validate a single read record (parallelizable).

    This function is standalone (not a method) to be picklable for multiprocessing.

    Args:
        record: BioPython SeqRecord object
        check_invalid_chars: If True, check for non-ATCGN characters
        allow_empty_id: If True, allow empty sequence IDs

    Returns:
        dict: {'success': True} if valid, or {'error': str, 'record_id': str, 'details': dict} if invalid
    """
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

class ReadValidator:
    """
    Validates and processes sequencing read files in FASTQ and BAM formats.

    This validator provides three validation levels:
        - 'strict': Full parsing and validation of all reads (slowest, most thorough)
        - 'trust': Fast validation - check line count (FASTQ) or header (BAM) and first records (~10-15x faster)
        - 'minimal': Only verify file format and copy as-is (fastest, FASTQ only)

    Processing Workflow:
        FASTQ: Parse → Validate → Apply edits → Save with optional compression
        BAM: Optionally copy original → Convert to FASTQ → Parse → Validate → Save

    Key Features:
        - FASTQ line count validation (must be divisible by 4)
        - Invalid character detection (non-ATCGN nucleotides)
        - Duplicate ID detection
        - BAM to FASTQ conversion (requires pysam or samtools)
        - NGS type-based output directory organization
        - Efficient compression conversion without full decompression

    Attributes:
        read_config: ReadConfig object with file path, format, and NGS type
        output_dir: Directory for output files
        settings: Settings object controlling validation and processing behavior
        sequences: List of parsed SeqRecord objects (populated during validation)

    Example:
        >>> from validation_pkg import ConfigManager, ReadValidator
        >>> config = ConfigManager.load("config.json")
        >>> settings = ReadValidator.Settings(validation_level='trust', outdir_by_ngs_type=True)
        >>> validator = ReadValidator(config.reads[0], config.output_dir, settings)
        >>> validator.validate()
    """

    @dataclass
    class Settings(BaseSettings):
        """
        Settings for read validation and processing.

        Attributes:
            check_invalid_chars: If True, validates that sequences contain only valid nucleotides (ATCGN)
            allow_empty_id: If True, allows sequences without IDs
            allow_duplicate_ids: If True, allows duplicate sequence IDs across the file
            keep_bam: If True, keeps a copy of the original BAM file when converting to FASTQ
            coding_type: Output compression type ('gz', 'bz2', or None for uncompressed)
            output_filename_suffix: Optional suffix to add to output filename (e.g., 'filtered')
            output_subdir_name: Optional subdirectory name within output_dir for saving files
            outdir_by_ngs_type: If True, automatically set output_subdir_name to ngs_type (default: False)

        Note:
            When outdir_by_ngs_type=True, output_subdir_name will be automatically set to the
            read's ngs_type (illumina/ont/pacbio), overriding any manually set value.

        """
        # Validation options
        check_invalid_chars: bool = False 
        allow_empty_id: bool = False
        allow_duplicate_ids: bool = True

        # Editing specifications
        keep_bam: bool = True 
        ignore_bam: bool = True

        # Output format
        coding_type: Optional[str] = None
        output_filename_suffix: Optional[str] = None
        output_subdir_name: Optional[str] = None
        outdir_by_ngs_type: bool = False

    def __init__(self, read_config, settings: Optional[Settings] = None) -> None:
        """
        Initialize read validator.

        Args:
            read_config: ReadConfig object from ConfigManager with file info
            settings: Settings object with validation parameters.
                If None, uses options from read_config.options (threads, validation_level),
                otherwise uses defaults.
        """
        self.logger = get_logger()

        # From genome global configuration
        self.read_config = read_config
        self.output_dir = read_config.output_dir
        self.validation_level = read_config.global_options.get("validation_level")   
        self.threads = read_config.global_options.get("threads") 
        self.input_path = read_config.filepath

        # From settings
        self.settings = settings if not None else self.Settings() 

        # Parsed data
        self.sequences = []  # List of SeqRecord objects

        print("post init")
        if not self.validation_level:
            self.validation_level = 'strict'    #   default global value

        if not self.threads:
            self.threads = 1    #   default global value

        # Apply outdir_by_ngs_type if enabled
        if self.settings.outdir_by_ngs_type:
            # Override output_subdir_name with ngs_type
            self.settings = self.settings.update(output_subdir_name=self.read_config.ngs_type)
            self.logger.debug(f"Applied outdir_by_ngs_type: output_subdir_name set to '{self.read_config.ngs_type}'")

    def run(self) -> None:
        """
        Main validation and processing workflow.

        Uses read_config data (format, compression) provided by ConfigManager.

        Workflow differs based on input format:
        - FASTQ: Parse → Validate → Edit → Save
        - BAM: Copy original → Convert to FASTQ → Parse → Validate → Edit → Save

        Raises:
            ReadValidationError: If validation fails
        """
        self.logger.info(f"Processing read file: {self.read_config.filename}")
        self.logger.debug(f"Format: {self.read_config.detected_format}, Compression: {self.read_config.coding_type}")

        try:
            # Special workflow for BAM files
            if self.read_config.detected_format == ReadFormat.BAM:
                if self.settings.ignore_bam:
                    self.logger.warning(f"Cannot proccess BAM files, the file will be ignored")
                    return

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

                # Step 4: Save to output directory (already FASTQ)
                self._save_output()

            self.logger.info(f"✓ Read validation completed")

        except Exception as e:
            self.logger.error(f"Read validation failed: {e}")
            raise

    def _open_file(self, mode: str = 'rt') -> IO:
        """
        Open file with automatic decompression based on read_config.

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
        """
        Fast line counting using Python file iteration.

        Uses Python's built-in gzip/bz2 libraries for secure line counting:
        - gzip (.gz): gzip.open() with text mode
        - bzip2 (.bz2): bz2.open() with text mode
        - uncompressed: open() with text mode

        This method is secure against command injection attacks (no shell commands)
        and works efficiently by iterating over lines without loading entire file
        into memory.

        Returns:
            Number of lines in the file

        Raises:
            ReadValidationError: If line counting fails

        Security:
            This implementation uses Python libraries exclusively (no subprocess or
            shell commands), preventing command injection vulnerabilities.
        """
        from validation_pkg.utils.file_handler import open_file_with_coding_type

        try:
            # Count lines using Python file iteration (no shell commands)
            # This is secure against command injection and works for all compression types
            with open_file_with_coding_type(self.input_path, self.read_config.coding_type, mode='rt') as f:
                line_count = sum(1 for _ in f)

            return line_count

        except TimeoutError:
            raise ReadValidationError("Line counting timed out - file too large")
        except (ValueError, IndexError) as e:
            raise ReadValidationError(f"Failed to parse line count: {e}")
        except Exception as e:
            raise ReadValidationError(f"Line counting failed: {e}")

    def _validate_first_fastq_record(self) -> Optional[SeqRecord]:
        """
        Validate only the first FASTQ record for format correctness.

        This method performs fast validation by checking only the first record
        without parsing the entire file. Used in trust mode for performance.

        Validation Checks:
            1. Line 1 starts with '@' (sequence ID)
            2. Line 2 contains sequence
            3. Line 3 starts with '+' (separator)
            4. Line 4 contains quality scores (same length as sequence)

        Returns:
            First SeqRecord if valid, None if validation fails

        Raises:
            FastqFormatError: If first record is not in valid FASTQ format
        """
        self.logger.debug("Validating first FASTQ record...")

        try:
            handle = self._open_file()
            try:
                # Read first 4 lines
                lines = []
                for i in range(4):
                    line = handle.readline()
                    if not line:
                        raise FastqFormatError(
                            f"File has fewer than 4 lines - not a valid FASTQ file"
                        )
                    lines.append(line.rstrip('\n\r'))

                # Validate FASTQ format
                if not lines[0].startswith('@'):
                    raise FastqFormatError(
                        f"First line must start with '@', found: {lines[0][:50]}"
                    )

                if not lines[2].startswith('+'):
                    raise FastqFormatError(
                        f"Third line must start with '+', found: {lines[2][:50]}"
                    )

                # Check that sequence and quality have same length
                seq_len = len(lines[1])
                qual_len = len(lines[3])
                if seq_len != qual_len:
                    raise FastqFormatError(
                        f"Sequence length ({seq_len}) doesn't match quality length ({qual_len})"
                    )

                # Parse first record using BioPython to ensure it's valid
                from io import StringIO
                first_record_text = '\n'.join(lines)
                record = next(SeqIO.parse(StringIO(first_record_text), 'fastq'))

                self.logger.debug(f"✓ First record validated: {record.id}, length={len(record.seq)}")
                return record

            finally:
                handle.close()

        except FastqFormatError:
            raise
        except Exception as e:
            error_msg = f"Failed to validate first FASTQ record: {e}"
            self.logger.add_validation_issue(
                level='ERROR',
                category='read',
                message=error_msg,
                details={'file': self.read_config.filename, 'error': str(e)}
            )
            raise FastqFormatError(error_msg) from e
    
    def _parse_file(self) -> None:
        """
        Parse FASTQ file using BioPython and validate format from read_config.

        Behavior depends on validation_level:
        - 'strict': Full parsing and validation of all sequences
        - 'trust': Only validate first record and check line count
        - 'minimal': No parsing or validation

        Note: BAM files are handled separately in the validate() method.
        This method only processes FASTQ files.
        """
        self.logger.debug(f"Parsing {self.read_config.detected_format} file (validation_level={self.validation_level})...")

        # Minimal mode - skip all parsing and validation except CodingType
        if self.validation_level == 'minimal':
            self.logger.info("Minimal validation mode - skipping file parsing")
            # Set empty sequences list - file will be copied as-is in _save_output
            self.sequences = []
            return

        # Trust mode: Parse only first 10 sequences for validation
        # Strict mode: Parse all sequences
        if self.validation_level == 'trust':
            self.logger.info("Trust mode - parsing first 10 sequences only for validation")
            parse_limit = 10
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
        """
        Validate parsed sequences with optional parallelization.

        Behavior depends on validation_level:
        - 'strict': Validate all sequences (optionally in parallel if threads > 1)
        - 'trust': Validate only first 10 sequences (sequential)
        - 'minimal': Skip (no sequences parsed)

        Parallelization:
        - Only enabled in strict mode when self.threads > 1
        - Uses multiprocessing.Pool for true parallelism (bypasses GIL)
        - Chunk-based processing to amortize overhead
        - All validation errors are collected before raising
        """
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
            # Enable parallel logging mode (adds process_id and thread_id to logs)
            self.logger.enable_parallel_logging()

            try:
                # Determine chunk size: min 1000 records per chunk, or split into 4 chunks per worker
                chunk_size = max(1000, validate_count // (self.threads * 4))

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

            finally:
                # Always disable parallel logging mode when done
                self.logger.disable_parallel_logging()

        # Sequential validation (trust mode or single-threaded)
        else:
            for idx in range(validate_count):
                record = self.sequences[idx]

                # Check sequence ID
                if not record.id and not self.settings.allow_empty_id:
                    error_msg = f"Sequence at index {idx} has no ID"
                    self.logger.add_validation_issue(
                        level='ERROR',
                        category='read',
                        message=error_msg,
                        details={'sequence_index': idx, 'sequence_id': str(record.id)}
                    )
                    raise ReadValidationError(error_msg)

                # Check for valid nucleotides
                if self.settings.check_invalid_chars:
                    valid_chars = set('ATCGNatcgn')
                    seq_str = str(record.seq)
                    invalid_chars = set(seq_str) - valid_chars

                    if invalid_chars:
                        char_count = sum(1 for c in seq_str if c in invalid_chars)
                        error_msg = f"Sequence '{record.id}' contains {char_count} invalid character(s): {', '.join(sorted(invalid_chars))}"
                        self.logger.add_validation_issue(
                            level='ERROR',
                            category='read',
                            message=error_msg,
                            details={
                                'sequence_id': record.id,
                                'invalid_chars': list(invalid_chars),
                                'count': char_count
                            }
                        )
                        raise ReadValidationError(error_msg)

        # Check for duplicate IDs (always sequential - requires full list)
        if not self.settings.allow_duplicate_ids:
            seq_ids = [record.id for record in self.sequences]
            if len(seq_ids) != len(set(seq_ids)):
                duplicates = [sid for sid in set(seq_ids) if seq_ids.count(sid) > 1]
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
        """
        Apply editing specifications to sequences based on settings.

        Future edits may include:
        - Remove short sequences (min_read_length)
        - Quality filtering
        - Adapter trimming
        - Read ID prefix/suffix modification

        Note: Parallelization should be added here when edits are implemented.
        """
        self.logger.debug("Applying editing specifications from settings...")

        # Currently no edits are implemented
        # This is a placeholder for future functionality

        self.logger.debug(f"✓ Edits applied, {len(self.sequences)} sequence(s) remaining")
    
    def _convert_bam_to_fastq(self) -> None:
        """
        Convert BAM file to FASTQ format using pysam or samtools.

        This method attempts to use pysam (Python library) first for better
        performance and control. If pysam is not available, it falls back to
        the samtools command-line tool.

        Behavior by validation level:
            - 'strict': Full conversion of all reads (skips unmapped/secondary/supplementary)
            - 'trust': Only validate BAM header and first 10 reads for format correctness
            - 'minimal': Skip conversion entirely

        Progress Reporting:
            Reports progress every 5% or 100k reads (whichever is appropriate)

        Raises:
            ReadValidationError: If conversion fails or neither pysam nor samtools is available
            BamFormatError: If BAM file is malformed or missing header
        """
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
                except:
                    # If counting fails, just proceed without progress reporting
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
        """
        Copy BAM file to output directory without processing.

        Returns:
            Path to output file
        """
        self.logger.debug("Copying BAM file to output directory...")

        # Determine output directory (with optional subdirectory)
        output_dir = self.output_dir
        if self.settings.output_subdir_name:
            output_dir = output_dir / self.settings.output_subdir_name

        # Create output directory
        output_dir.mkdir(parents=True, exist_ok=True)

        # Generate output filename - keep original name and compression
        output_filename = self.read_config.filename
        if self.settings.output_filename_suffix:
            # Insert suffix before extensions
            base_name = self.read_config.filename
            for suffix in self.input_path.suffixes:
                base_name = base_name.replace(suffix, '')
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
        """
        Save processed read to output directory using settings.

        Behavior depends on validation_level:
        - 'strict': Write sequences using BioPython
        - 'trust' or 'minimal': Copy file as-is (preserve original format)

        Note: Compression handling will be moved to utils/file_handler.py in future.

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
        base_name = self.read_config.filename

        # Remove all suffixes (.fastq.gz -> remove .gz then .fastq)
        for suffix in self.input_path.suffixes:
            base_name = base_name.replace(suffix, '')
        # Add suffix if specified
        if self.settings.output_filename_suffix:
            output_filename = f"{base_name}_{self.settings.output_filename_suffix}.fastq"
        else:
            output_filename = f"{base_name}.fastq"


        # Add compression extension if requested (from settings)
        # Normalize settings.coding_type to CodingType enum
        coding = CT(self.settings.coding_type) if self.settings.coding_type else CT.NONE

        if coding == CT.GZIP:
            output_filename += '.gz'
        elif coding == CT.BZIP2:
            output_filename += '.bz2'

        output_path = output_dir / output_filename

        # Minimal mode - copy file as-is without parsing
        # Required: FASTQ format + coding must match settings.coding_type
        if self.validation_level == 'minimal':
            self.logger.debug("Minimal mode - validating format and coding requirements")

            # Check format - must be FASTQ
            if self.read_config.detected_format != ReadFormat.FASTQ:
                error_msg = f'Minimal mode requires FASTQ format, got {self.read_config.detected_format}. Use validation_level "trust" or "strict" to convert.'
                self.logger.add_validation_issue(
                    level='ERROR',
                    category='read',
                    message=error_msg,
                    details={'file': self.read_config.filename, 'detected_format': str(self.read_config.detected_format)}
                )
                raise ReadValidationError(error_msg)

            # Check coding - must match settings.coding_type
            # read_config.coding_type is CodingType enum, coding is normalized above
            input_coding = self.read_config.coding_type
            required_coding = coding

            if input_coding != required_coding:
                error_msg = f'Minimal mode requires input coding to match output coding. Input: {input_coding}, Required: {required_coding}. Use validation_level "trust" or "strict" to change compression.'
                self.logger.add_validation_issue(
                    level='ERROR',
                    category='read',
                    message=error_msg,
                    details={'file': self.read_config.filename, 'input_coding': str(input_coding), 'required_coding': str(required_coding)}
                )
                raise ReadValidationError(error_msg)

            # Check line count divisible by 4 (fast validation)
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

            # Use the output_filename already constructed (includes correct extension)
            output_path_minimal = output_dir / output_filename

            # Copy file
            self.logger.debug(f"Copying {self.input_path} to {output_path_minimal}")
            shutil.copy2(self.input_path, output_path_minimal)

            self.logger.info(f"Output saved: {output_path_minimal}")
            return output_path_minimal

        # Trust mode - copy original file with coding conversion
        # Trust mode parsed only first 10 sequences for validation, so we copy the original file
        if self.validation_level == 'trust':
            self.logger.debug("Trust mode - copying original file with coding conversion")

            # Both are CodingType enum: read_config.coding_type and coding (normalized above)
            input_coding = self.read_config.coding_type
            output_coding = coding

            # If input and output coding match, simple copy
            if input_coding == output_coding:
                self.logger.debug(f"Copying {self.input_path} to {output_path} (same coding)")
                shutil.copy2(self.input_path, output_path)
            else:
                # Convert coding using file_handler utilities
                self.logger.debug(f"Converting {input_coding} -> {output_coding}")

                # Map (input, output) pairs to conversion functions
                conversion_map = {
                    (CT.BZIP2, CT.GZIP): bz2_to_gz,
                    (CT.NONE, CT.GZIP): none_to_gz,
                    (CT.GZIP, CT.NONE): gz_to_none,
                    (CT.BZIP2, CT.NONE): bz2_to_none,
                    (CT.NONE, CT.BZIP2): none_to_bz2,
                    (CT.GZIP, CT.BZIP2): gz_to_bz2,
                }

                conversion_func = conversion_map.get((input_coding, output_coding))
                if conversion_func:
                    conversion_func(self.input_path, output_path, threads=self.threads)
                else:
                    error_msg = f"Unsupported coding conversion: {input_coding} -> {output_coding}"
                    raise ReadValidationError(error_msg)

            self.logger.info(f"Output saved: {output_path}")
            return output_path

        # Strict mode - write sequences using BioPython (original behavior)
        # Write output with appropriate compression
        self.logger.debug(f"Writing output to: {output_path}")

        # Use optimized compression writer
        with open_compressed_writer(output_path, coding, threads=self.threads) as handle:
            SeqIO.write(self.sequences, handle, 'fastq')

        self.logger.info(f"Output saved: {output_path}")

        return output_path
    