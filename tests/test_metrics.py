import math
from lr_rna_qc.metrics import calculate_n50, compute_run_summary
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
    qual = "+?"  # ASCII 33 + 10 = 43 ('+'), ASCII 33 + 30 = 63 ('?')
    qc = process_quality_string(qual)

    naive_mean = 20.0
    expected_exact = -10.0 * math.log10((0.1 + 0.001) / 2.0)

    assert abs(qc.mean_q - expected_exact) < 1e-4
    assert abs(qc.mean_q - naive_mean) > 6.0


def test_compute_run_summary():
    lengths = [500, 1000, 1500]
    mean_qs = [10.0, 12.0, 14.0]
    summary = compute_run_summary(lengths, mean_qs, pass_reads=2, fail_reads=1)

    assert summary.total_reads == 3
    assert summary.total_bases == 3000
    assert summary.min_length == 500
    assert summary.max_length == 1500
    assert summary.pass_reads == 2
    assert summary.fail_reads == 1
    assert summary.n50 == 1500

    d = summary.to_dict()
    assert d["total_reads"] == 3
    assert d["n50"] == 1500


def test_compute_run_summary_empty():
    summary = compute_run_summary([], [])
    assert summary.total_reads == 0
    assert summary.total_bases == 0
    assert summary.n50 == 0
