from pathlib import Path
from click.testing import CliRunner
from lr_rna_qc.cli import main
from tests.make_synthetic_fastq import generate_synthetic_fastq


def test_cli_help():
    runner = CliRunner()
    result = runner.invoke(main, ["--help"])
    assert result.exit_code == 0
    assert "scan" in result.output
    assert "compare" in result.output
    assert "benchmark" in result.output


def test_cli_scan(tmp_path: Path):
    runner = CliRunner()
    fq = generate_synthetic_fastq(tmp_path / "sample.fastq.gz", num_reads=30)
    plot_path = tmp_path / "report.png"
    json_path = tmp_path / "summary.json"
    filt_path = tmp_path / "filtered.fastq.gz"

    result = runner.invoke(
        main,
        [
            "scan",
            str(fq),
            "--min-q", "10",
            "--min-len", "200",
            "--out-filtered", str(filt_path),
            "--plot", str(plot_path),
            "--json", str(json_path),
            "--backend", "auto",
        ],
    )
    assert result.exit_code == 0
    assert "TRANSCRIPTOME QC METRICS" in result.output
    assert "FILTERING SUMMARY" in result.output
    assert "SENSOR DEGRADATION CORRELATION" in result.output
    assert plot_path.exists()
    assert json_path.exists()
    assert filt_path.exists()


def test_cli_compare(tmp_path: Path):
    runner = CliRunner()
    s1 = generate_synthetic_fastq(tmp_path / "sample1.fastq.gz", num_reads=25, seed=1)
    s2 = generate_synthetic_fastq(tmp_path / "sample2.fastq.gz", num_reads=25, seed=2)
    comp_json = tmp_path / "compare.json"

    result = runner.invoke(
        main,
        [
            "compare",
            str(s1),
            str(s2),
            "--out-json", str(comp_json),
        ],
    )
    assert result.exit_code == 0
    assert "Comparative Transcriptome Hypothesis Testing" in result.output
    assert "READ LENGTH DISTRIBUTION TEST" in result.output
    assert "PHRED QUALITY DISTRIBUTION TEST" in result.output
    assert comp_json.exists()


def test_cli_benchmark(tmp_path: Path):
    runner = CliRunner()
    fq = generate_synthetic_fastq(tmp_path / "bench.fastq.gz", num_reads=50)

    result = runner.invoke(main, ["benchmark", str(fq)])
    assert result.exit_code == 0
    assert "Backend Engine Performance Benchmark" in result.output
    assert "Python" in result.output
