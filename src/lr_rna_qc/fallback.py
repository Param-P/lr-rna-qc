"""Pure-Python fallback engine for Phred quality and sequence metric calculations."""

import math
from typing import NamedTuple, List

class ReadQC(NamedTuple):
    mean_q: float
    length: int


# Precalculate 256-element Look-Up Table (ASCII -> error probability)
# P = 10^(-Q / 10), where Q = ASCII - 33
ASCII_TO_PROB = [1.0] * 256
for i in range(33, 127):
    q = float(i - 33)
    ASCII_TO_PROB[i] = 10.0 ** (-q / 10.0)


def process_quality_string(qual: str) -> ReadQC:
    """Calculate exact mean quality score using probability Look-Up Table.

    Q_mean = -10 * log10(mean(P_error))
    """
    length = len(qual)
    if length == 0:
        return ReadQC(mean_q=0.0, length=0)

    # Lookup error probability for each character byte
    sum_p = sum(ASCII_TO_PROB[ord(c)] for c in qual)
    mean_p = sum_p / length

    if mean_p <= 0.0:
        mean_q = 93.0
    else:
        mean_q = -10.0 * math.log10(mean_p)

    return ReadQC(mean_q=mean_q, length=length)


def process_fastq_chunk(
    headers: List[str],
    seqs: List[str],
    quals: List[str],
) -> List[ReadQC]:
    """Batch-process FASTQ chunk in pure Python."""
    return [process_quality_string(q) for q in quals]
