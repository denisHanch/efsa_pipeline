"""Common validation helper functions."""

from pathlib import Path
from typing import Any
import shutil

from validation_pkg.exceptions import ValidationError


def validate_minimal_mode_requirements(
    detected_format: Any,
    expected_format: Any,
    input_coding: Any,
    required_coding: Any,
    filename: str,
    logger,
    category: str
) -> None:
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
    if detected_format != expected_format:
        error_msg = (
            f'Minimal mode requires {expected_format} format, '
            f'got {detected_format}. '
            f'Use validation_level "trust" or "strict" to convert.'
        )
        logger.add_validation_issue(
            level='ERROR',
            category=category,
            message=error_msg,
            details={
                'file': filename,
                'detected_format': str(detected_format),
                'expected_format': str(expected_format)
            }
        )
        raise ValidationError(error_msg)

    # Check coding matches
    if input_coding != required_coding:
        error_msg = (
            f'Minimal mode requires input coding to match output coding. '
            f'Input: {input_coding}, Required: {required_coding}. '
            f'Use validation_level "trust" or "strict" to change compression.'
        )
        logger.add_validation_issue(
            level='ERROR',
            category=category,
            message=error_msg,
            details={
                'file': filename,
                'input_coding': str(input_coding),
                'required_coding': str(required_coding)
            }
        )
        raise ValidationError(error_msg)


def copy_file_minimal_mode(
    input_path: Path,
    output_path: Path,
    logger
) -> Path:
    """
    Copy file for minimal validation mode.

    Args:
        input_path: Source file path
        output_path: Destination file path
        logger: Logger instance

    Returns:
        Output path
    """
    logger.debug(f"Minimal mode - copying {input_path} to {output_path}")
    shutil.copy2(input_path, output_path)
    logger.info(f"Output saved: {output_path}")
    return output_path
