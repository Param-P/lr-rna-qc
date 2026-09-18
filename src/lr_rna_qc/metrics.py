"""Transcriptome-specific QC metrics: N50, base yields, and summary dataclasses."""

from dataclasses import dataclass
from typing import List, Sequence


@dataclass
class RunSummary:
    """Summary metrics of a Nanopore RNA sequencing run."""
    total_reads: int
    total_bases: int
    n50: int
    mean_q: float
    median_q: float
    min_length: int
    max_length: int
    pass_reads: int
    fail_reads: int


def calculate_n50(lengths: Sequence[int]) -> int:
    """Calculate the N50 metric (length-weighted median).

    N50 is the read length L such that sum of reads with length >= L
    accounts for at least 50% of the total sequenced bases.
    """
    if not lengths:
        return 0
    sorted_lens = sorted(lengths, reverse=True)
    half_total = sum(sorted_lens) / 2.0
    cum_sum = 0
    for length in sorted_lens:
        cum_sum += length
        if cum_sum >= half_total:
            return length
    return sorted_lens[-1]
