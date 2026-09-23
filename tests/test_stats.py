from pathlib import Path
from lr_rna_qc.stats import compute_correlation, compare_distributions
from lr_rna_qc.plot import generate_qc_report


def test_correlation_positive_and_negative():
    lengths = [100, 200, 300, 400, 500]
    # Positively correlated
    qs_pos = [10.0, 11.0, 12.0, 13.0, 14.0]
    res_pos = compute_correlation(lengths, qs_pos)
    assert res_pos["pearson_r"] > 0.99
    assert res_pos["spearman_rho"] > 0.99
    assert res_pos["pearson_p"] < 0.01

    # Negatively correlated (simulating sensor decay)
    qs_neg = [15.0, 14.0, 13.0, 12.0, 11.0]
    res_neg = compute_correlation(lengths, qs_neg)
    assert res_neg["pearson_r"] < -0.99
    assert res_neg["spearman_rho"] < -0.99


def test_correlation_edge_cases():
    # Empty
    assert compute_correlation([], [])["pearson_r"] == 0.0
    # Single element
    assert compute_correlation([100], [12.0])["pearson_r"] == 0.0
    # Zero variance in lengths or qualities
    res_flat = compute_correlation([100, 100, 100], [12.0, 12.0, 12.0])
    assert res_flat["pearson_r"] == 0.0
    assert res_flat["pearson_p"] == 1.0


def test_distribution_comparison():
    # Identical distributions
    lens = [500, 1000, 1500, 2000, 2500]
    qs = [10.0, 12.0, 11.0, 13.0, 12.5]
    comp_same = compare_distributions(lens, qs, lens, qs)
    assert comp_same["length_comparison"]["significant"] is False
    assert comp_same["quality_comparison"]["significant"] is False

    # Distinct distributions (e.g. wild-type vs severe truncation mutant)
    lens_short = [50, 70, 90, 110, 130]
    qs_low = [4.0, 5.0, 4.5, 5.5, 4.0]
    comp_diff = compare_distributions(lens, qs, lens_short, qs_low)
    assert comp_diff["length_comparison"]["significant"] is True
    assert comp_diff["quality_comparison"]["significant"] is True
    assert comp_diff["length_comparison"]["mann_whitney_p"] < 0.05


def test_distribution_comparison_edge_cases():
    comp = compare_distributions([], [], [100], [12.0])
    assert comp["length_comparison"]["significant"] is False
    assert comp["quality_comparison"]["significant"] is False


def test_generate_qc_report(tmp_path: Path):
    out_png = tmp_path / "qc_report.png"
    lengths = [200, 500, 1200, 1800, 3000, 8000]
    qs = [12.5, 11.8, 10.9, 10.2, 9.5, 8.8]
    n50 = 3000

    report_path = generate_qc_report(lengths, qs, n50, out_png, sample_name="Test Flowcell")
    assert report_path.exists()
    assert report_path.stat().st_size > 1000  # Non-trivial PNG file size


def test_generate_qc_report_empty(tmp_path: Path):
    out_png = tmp_path / "empty_report.png"
    report_path = generate_qc_report([], [], 0, out_png)
    assert report_path.exists()
    assert report_path.stat().st_size > 500
