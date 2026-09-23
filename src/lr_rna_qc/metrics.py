from dataclasses import dataclass, asdict
import math
import statistics
from typing import Sequence, Dict, Any


@dataclass
class RunSummary:
    total_reads: int
    total_bases: int
    n50: int
    mean_q: float
    median_q: float
    min_length: int
    max_length: int
    pass_reads: int
    fail_reads: int

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


def calculate_n50(lengths: Sequence[int]) -> int:
    """Calculate the N50 metric.

    Length at which cumulative yield from descending reads reaches >= 50% of total bases.
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


def compute_run_summary(
    lengths: Sequence[int],
    mean_qs: Sequence[float],
    pass_reads: int = 0,
    fail_reads: int = 0,
) -> RunSummary:
    """Compute aggregate QC metrics for a sequenced run."""
    if not lengths:
        return RunSummary(
            total_reads=0,
            total_bases=0,
            n50=0,
            mean_q=0.0,
            median_q=0.0,
            min_length=0,
            max_length=0,
            pass_reads=pass_reads,
            fail_reads=fail_reads,
        )

    total_reads = len(lengths)
    total_bases = sum(lengths)
    n50 = calculate_n50(lengths)

    # Average Phred score in probability space to avoid logarithmic bias
    mean_prob = sum(math.pow(10.0, -q / 10.0) for q in mean_qs) / total_reads
    mean_q = -10.0 * math.log10(mean_prob) if mean_prob > 0.0 else 93.0
    median_q = float(statistics.median(mean_qs))

    return RunSummary(
        total_reads=total_reads,
        total_bases=total_bases,
        n50=n50,
        mean_q=round(mean_q, 2),
        median_q=round(median_q, 2),
        min_length=min(lengths),
        max_length=max(lengths),
        pass_reads=pass_reads,
        fail_reads=fail_reads,
    )
