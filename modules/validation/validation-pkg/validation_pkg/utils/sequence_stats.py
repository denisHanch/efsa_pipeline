"""Utility functions for sequence statistics calculations."""

from typing import List


def calculate_n50(lengths: List[int]) -> int:
    """Calculate N50 assembly/read quality metric."""
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
