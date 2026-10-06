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


def test_n50_identical_lengths():
    assert calculate_n50([800, 800, 800, 800]) == 800


def test_n50_skewed_transcriptome_yield():
    # One long transcript (9000 bp) and 10 short fragments (100 bp each)
    # Total bases = 10,000 bp. 50% threshold = 5,000 bp.
    # The 9000 bp read alone exceeds 5000 bp -> N50 must be 9000
    lengths = [9000] + [100] * 10
    assert calculate_n50(lengths) == 9000


def test_phred_exact_vs_arithmetic():
    qual = "+?"  # ASCII 33 + 10 = 43 ('+'), ASCII 33 + 30 = 63 ('?')
    qc = process_quality_string(qual)

    naive_mean = 20.0
    expected_exact = -10.0 * math.log10((0.1 + 0.001) / 2.0)

    assert abs(qc.mean_q - expected_exact) < 1e-4
    assert abs(qc.mean_q - naive_mean) > 6.0


def test_phred_severe_low_quality_burst():
    # 5 bases of Q=2 ('#', P=0.630957) and 5 bases of Q=30 ('?', P=0.001)
    # Naive mean = (5*2 + 5*30)/10 = 16.0
    # True mean P = (5 * 0.63095734 + 5 * 0.001) / 10 = 0.31597867
    # True Q_mean = -10 * log10(0.31597867) ~= 5.003
    qual = "#####?????"
    qc = process_quality_string(qual)

    p_low = 10.0 ** (-2.0 / 10.0)
    p_high = 10.0 ** (-30.0 / 10.0)
    expected_exact = -10.0 * math.log10((5.0 * p_low + 5.0 * p_high) / 10.0)

    assert abs(qc.mean_q - expected_exact) < 1e-4
    assert qc.mean_q < 6.0  # Naive 16.0 is dangerously misleading


def test_phred_uniform_quality_no_drift():
    # If all bases are Q=25 (ASCII 33+25 = 58 ':'), mean Q must be exactly 25.0
    qual = ":" * 50
    qc = process_quality_string(qual)
    assert abs(qc.mean_q - 25.0) < 1e-5


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
    assert summary.median_q == 12.0

    d = summary.to_dict()
    assert d["total_reads"] == 3
    assert d["n50"] == 1500


def test_compute_run_summary_empty():
    summary = compute_run_summary([], [])
    assert summary.total_reads == 0
    assert summary.total_bases == 0
    assert summary.n50 == 0
    assert summary.mean_q == 0.0
