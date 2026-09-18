"""Unit tests for transcriptome QC metrics (N50, Phred probabilities)."""

import math
import pytest
from lr_rna_qc.metrics import calculate_n50
from lr_rna_qc.fallback import process_quality_string


def test_n50_analytical():
    """Verify N50 calculation against known hand-calculated array."""
    # Lengths: [100, 200, 300, 400], total = 1000, half = 500
    # Descending: 400 (cum 400 < 500), 300 (cum 700 >= 500) -> N50 = 300
    lengths = [100, 200, 300, 400]
    assert calculate_n50(lengths) == 300


def test_n50_single_read():
    assert calculate_n50([1500]) == 1500


def test_n50_empty():
    assert calculate_n50([]) == 0


def test_phred_exact_vs_arithmetic():
    """Verify exact logarithmic error probability conversion avoids arithmetic bias.

    If a read has one Q=10 base (P=0.1) and one Q=30 base (P=0.001):
    Naive arithmetic mean Q = (10 + 30) / 2 = 20.0
    True mean error probability = (0.1 + 0.001) / 2 = 0.0505
    True mean Q = -10 * log10(0.0505) ~= 12.967
    """
    # ASCII 33 + 10 = 43 ('+'), ASCII 33 + 30 = 63 ('?')
    qual = "+?"
    qc = process_quality_string(qual)

    naive_mean = 20.0
    expected_exact = -10.0 * math.log10((0.1 + 0.001) / 2.0)

    assert abs(qc.mean_q - expected_exact) < 1e-4
    assert abs(qc.mean_q - naive_mean) > 6.0, "Arithmetic mean falsely overestimates quality!"
