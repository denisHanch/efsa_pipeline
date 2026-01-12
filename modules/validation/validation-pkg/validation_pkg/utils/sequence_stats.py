"""Utility functions for sequence statistics calculations."""

from typing import List


def calculate_n50(lengths: List[int]) -> int:
    """Calculate N50 assembly/read quality metric.

    N50 is defined as the sequence length of the shortest contig/read at which
    at least 50% of the total assembly/read set is contained in contigs/reads
    of that length or longer.

    Args:
        lengths: List of sequence lengths (unsorted)

    Returns:
        N50 value as integer. Returns 0 if lengths list is empty.

    Examples:
        >>> calculate_n50([100, 200, 300, 400, 500])
        400
        >>> calculate_n50([10, 20, 30])
        20
        >>> calculate_n50([])
        0
    """
    if not lengths:
        return 0

    # Sort lengths in descending order
    sorted_lengths = sorted(lengths, reverse=True)
    total_length = sum(sorted_lengths)

    # Find N50: the length at which cumulative length >= 50% of total
    cumulative_length = 0
    for length in sorted_lengths:
        cumulative_length += length
        if cumulative_length >= total_length / 2:
            return length

    return 0
