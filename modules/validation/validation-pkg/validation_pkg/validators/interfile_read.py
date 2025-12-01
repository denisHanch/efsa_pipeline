"""
Inter-file validation for read files.

Provides validation functions for checking consistency across multiple
read files, particularly for paired-end read completeness (R1 ↔ R2 matching).

Usage:
    from validation_pkg import validate_reads, readxread_validation, ReadXReadSettings

    # Validate individual files
    reads_results = validate_reads(config.reads, settings)

    # Check inter-file consistency
    readxread_settings = ReadXReadSettings(pair_end_basename=True)
    result = readxread_validation(reads_results, readxread_settings)
"""

from dataclasses import dataclass
from typing import List, Dict, Optional, Any
from ..utils.settings import BaseSettings
from ..exceptions import ReadValidationError
from ..logger import get_logger


@dataclass
class ReadXReadSettings(BaseSettings):
    """
    Settings for read-to-read inter-file validation.

    Attributes:
        pair_end_basename: Check that R2 files have matching R1 files
        allow_missing_r1: Allow R2 files without matching R1 (for incomplete data)
    """
    pair_end_basename: bool = True
    allow_missing_r1: bool = False

    def __post_init__(self):
        """Validate settings."""
        # Nothing to validate currently, but follows pattern
        pass


def readxread_validation(
    reads_results: List[Dict[str, Any]],
    settings: Optional[ReadXReadSettings] = None
) -> Dict[str, Any]:
    """
    Validate consistency across multiple read files.

    Currently supports:
    - Paired-end completeness: R2 files must have matching R1 files

    Args:
        reads_results: List of result dicts from validate_reads()
                      Each dict should contain 'output_metadata' key with
                      base_name, read_number, ngs_type_detected fields
        settings: Validation settings (uses defaults if not provided)

    Returns:
        Dict with validation results and metadata

    Raises:
        ReadValidationError: If critical validation fails (missing R1 for R2)
    """
    settings = settings or ReadXReadSettings()
    logger = get_logger()

    warnings = []
    errors = []
    metadata = {
        'pairs_checked': 0,
        'complete_pairs': [],
        'missing_r1': [],
        'duplicate_r1': [],
        'duplicate_r2': []
    }

    logger.info("Running inter-file validation: read-to-read consistency")

    if not settings.pair_end_basename:
        logger.debug("Paired-end basename check disabled, skipping")
        return {
            'passed': True,
            'warnings': warnings,
            'errors': errors,
            'metadata': metadata
        }

    # Extract metadata from results (OutputMetadata objects)
    read_metadata = []
    for result in reads_results:
        # Extract fields from OutputMetadata object
        output_file = result.output_file or 'unknown'
        base_name = result.base_name
        read_number = result.read_number
        ngs_type_detected = result.ngs_type_detected

        # Skip if no pattern detected (not paired-end or not Illumina)
        if read_number is None:
            logger.debug(f"No paired-end pattern detected for {output_file}, skipping")
            continue

        read_metadata.append({
            'output_file': output_file,
            'base_name': base_name,
            'read_number': read_number,
            'ngs_type_detected': ngs_type_detected
        })

    if not read_metadata:
        logger.info("No paired-end patterns detected in results, validation passed (nothing to check)")
        return {
            'passed': True,
            'warnings': warnings,
            'errors': errors,
            'metadata': metadata
        }

    logger.debug(f"Found {len(read_metadata)} files with paired-end patterns")

    # Build mapping: {base_name: {1: [files], 2: [files]}}
    pairs: Dict[str, Dict[int, List[str]]] = {}

    for meta in read_metadata:
        base_name = meta['base_name']
        read_num = meta['read_number']
        output_file = meta['output_file']

        if base_name not in pairs:
            pairs[base_name] = {1: [], 2: []}

        pairs[base_name][read_num].append(output_file)

    metadata['pairs_checked'] = len(pairs)

    # Validate completeness: every R2 must have matching R1
    for base_name, reads in pairs.items():
        r1_files = reads[1]
        r2_files = reads[2]

        # Check for R2 without R1
        if r2_files and not r1_files:
            if not settings.allow_missing_r1:
                error_msg = f"Found R2 file(s) without matching R1 for base name '{base_name}': {r2_files}"
                errors.append(error_msg)
                metadata['missing_r1'].extend(r2_files)
                logger.error(error_msg)
            else:
                warning_msg = f"Found R2 file(s) without matching R1 for base name '{base_name}': {r2_files} (allowed by settings)"
                warnings.append(warning_msg)
                logger.warning(warning_msg)

        # Check for complete pairs
        elif r1_files and r2_files:
            metadata['complete_pairs'].append(base_name)
            logger.debug(f"Complete pair found for '{base_name}': R1={r1_files}, R2={r2_files}")

        # R1 only (no R2) - this is valid, just log
        elif r1_files and not r2_files:
            logger.debug(f"R1-only file for '{base_name}': {r1_files} (valid, no R2 required)")

        # Warn about multiple files with same read number
        if len(r1_files) > 1:
            warning_msg = f"Multiple R1 files found for base name '{base_name}': {r1_files}"
            warnings.append(warning_msg)
            metadata['duplicate_r1'].extend(r1_files)
            logger.warning(warning_msg)

        if len(r2_files) > 1:
            warning_msg = f"Multiple R2 files found for base name '{base_name}': {r2_files}"
            warnings.append(warning_msg)
            metadata['duplicate_r2'].extend(r2_files)
            logger.warning(warning_msg)

    # Determine if validation passed (no ERROR-level issues)
    passed = len(errors) == 0

    if passed:
        logger.info(f"✓ Read inter-file validation passed: {metadata['pairs_checked']} base name(s) checked, {len(metadata['complete_pairs'])} complete pair(s)")
    else:
        logger.error(f"✗ Read inter-file validation failed: {len(errors)} error(s), {len(warnings)} warning(s)")

    result = {
        'passed': passed,
        'warnings': warnings,
        'errors': errors,
        'metadata': metadata
    }

    # Raise exception if failed
    if not passed:
        error_summary = '\n  - '.join([''] + errors)
        raise ReadValidationError(f"Read inter-file validation failed:{error_summary}")

    return result
