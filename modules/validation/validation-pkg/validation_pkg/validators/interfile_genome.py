"""
Inter-file validation for genome files.

Provides validation functions for checking consistency between genome
files, such as comparing reference vs modified genomes.
"""

from dataclasses import dataclass
from typing import Dict, Optional, Any
from ..utils.settings import BaseSettings
from ..exceptions import GenomeValidationError
from ..logger import get_logger


def _get_metadata_field(metadata_obj, field_name, default=None):
    """Helper to get field from either OutputMetadata object or dict."""
    if hasattr(metadata_obj, field_name):
        # OutputMetadata object - use attribute access
        return getattr(metadata_obj, field_name, default)
    else:
        # Dictionary - use dict access
        return metadata_obj.get(field_name, default)


@dataclass
class GenomeXGenomeSettings(BaseSettings):
    """
    Settings for genome-to-genome inter-file validation.

    Attributes:
        same_number_of_sequences: Require same number of sequences (chromosomes/contigs)
        same_sequence_ids: Require matching sequence IDs (order-independent)
        same_sequence_lengths: Require matching sequence lengths for common IDs
    """
    same_number_of_sequences: bool = True
    same_sequence_ids: bool = False
    same_sequence_lengths: bool = False

    def __post_init__(self):
        """Validate settings."""
        # If checking lengths, must also check IDs
        if self.same_sequence_lengths and not self.same_sequence_ids:
            raise ValueError(
                "same_sequence_lengths requires same_sequence_ids=True "
                "(cannot compare lengths without matching IDs)"
            )


def genomexgenome_validation(
    ref_genome_result,  # OutputMetadata or Dict[str, Any]
    mod_genome_result,  # OutputMetadata or Dict[str, Any]
    settings: Optional[GenomeXGenomeSettings] = None
) -> Dict[str, Any]:
    """
    Validate consistency between two genome files (e.g., reference vs modified).

    Currently supports:
    - Sequence count matching
    - Sequence ID matching (order-independent)
    - Sequence length matching for common IDs

    Args:
        ref_genome_result: Result dict from validate_genome() for reference genome
        mod_genome_result: Result dict from validate_genome() for modified genome
        settings: Validation settings (uses defaults if not provided)

    Returns:
        Dict with validation results and metadata

    Raises:
        GenomeValidationError: If critical validation fails
    """
    settings = settings or GenomeXGenomeSettings()
    logger = get_logger()

    warnings = []
    errors = []

    logger.info("Running inter-file validation: genome-to-genome consistency")

    # Extract metadata from OutputMetadata objects
    ref_meta = ref_genome_result
    mod_meta = mod_genome_result

    # Initialize result metadata
    metadata = {
        'ref_num_sequences': _get_metadata_field(ref_meta, 'num_sequences', 0),
        'mod_num_sequences': _get_metadata_field(mod_meta, 'num_sequences', 0),
        'common_sequence_ids': [],
        'ref_only_ids': [],
        'mod_only_ids': [],
        'length_mismatches': {}
    }

    # Check 1: Same number of sequences
    if settings.same_number_of_sequences:
        ref_count = _get_metadata_field(ref_meta, 'num_sequences', 0)
        mod_count = _get_metadata_field(mod_meta, 'num_sequences', 0)

        if ref_count != mod_count:
            error_msg = (
                f"Genome sequence count mismatch: "
                f"reference has {ref_count} sequence(s), "
                f"modified has {mod_count} sequence(s)"
            )
            errors.append(error_msg)
            logger.error(error_msg)
        else:
            logger.debug(f"Sequence count match: {ref_count} sequence(s)")

    # Check 2: Same sequence IDs
    if settings.same_sequence_ids:
        ref_ids = set(_get_metadata_field(ref_meta, 'sequence_ids', []))
        mod_ids = set(_get_metadata_field(mod_meta, 'sequence_ids', []))

        common_ids = ref_ids & mod_ids
        ref_only = ref_ids - mod_ids
        mod_only = mod_ids - ref_ids

        metadata['common_sequence_ids'] = sorted(common_ids)
        metadata['ref_only_ids'] = sorted(ref_only)
        metadata['mod_only_ids'] = sorted(mod_only)

        if ref_only or mod_only:
            error_parts = []
            if ref_only:
                error_parts.append(f"reference-only: {sorted(ref_only)}")
            if mod_only:
                error_parts.append(f"modified-only: {sorted(mod_only)}")

            error_msg = f"Genome sequence ID mismatch: {', '.join(error_parts)}"
            errors.append(error_msg)
            logger.error(error_msg)
        else:
            logger.debug(f"Sequence IDs match: {len(common_ids)} common ID(s)")

        # Check 3: Same sequence lengths for common IDs
        if settings.same_sequence_lengths and common_ids:
            ref_lengths = _get_metadata_field(ref_meta, 'sequence_lengths', {})
            mod_lengths = _get_metadata_field(mod_meta, 'sequence_lengths', {})

            for seq_id in common_ids:
                ref_len = ref_lengths.get(seq_id)
                mod_len = mod_lengths.get(seq_id)

                if ref_len is None or mod_len is None:
                    warning_msg = f"Missing length information for sequence '{seq_id}'"
                    warnings.append(warning_msg)
                    logger.warning(warning_msg)
                    continue

                if ref_len != mod_len:
                    metadata['length_mismatches'][seq_id] = {
                        'ref_length': ref_len,
                        'mod_length': mod_len,
                        'difference': mod_len - ref_len
                    }
                    error_msg = (
                        f"Sequence length mismatch for '{seq_id}': "
                        f"reference={ref_len} bp, modified={mod_len} bp "
                        f"(diff={mod_len - ref_len:+d} bp)"
                    )
                    errors.append(error_msg)
                    logger.error(error_msg)

            if not metadata['length_mismatches']:
                logger.debug(f"Sequence lengths match for {len(common_ids)} sequence(s)")

    # Determine if validation passed (no ERROR-level issues)
    passed = len(errors) == 0

    if passed:
        logger.info(f"✓ Genome inter-file validation passed")
    else:
        logger.error(f"✗ Genome inter-file validation failed: {len(errors)} error(s), {len(warnings)} warning(s)")

    result = {
        'passed': passed,
        'warnings': warnings,
        'errors': errors,
        'metadata': metadata
    }

    # Raise exception if failed
    if not passed:
        error_summary = '\n  - '.join([''] + errors)
        raise GenomeValidationError(f"Genome inter-file validation failed:{error_summary}")

    return result
