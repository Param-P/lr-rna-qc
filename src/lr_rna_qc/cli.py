import json
import time
from pathlib import Path
import click
from lr_rna_qc.core import get_backend, scan_fastq
from lr_rna_qc.stats import compute_correlation, compare_distributions
from lr_rna_qc.plot import generate_qc_report


@click.group()
@click.version_option(package_name="lr-rna-qc")
def main():
    """lr-rna-qc: Oxford Nanopore RNA-seq Quality Control & Splicing Analysis Tool."""
    pass


@main.command()
@click.argument("fastq_file", type=click.Path(exists=True, dir_okay=False))
@click.option("--min-q", default=10.0, show_default=True, type=float, help="Minimum mean Phred quality threshold.")
@click.option("--min-len", default=200, show_default=True, type=int, help="Minimum read length threshold.")
@click.option("--out-filtered", type=click.Path(dir_okay=False), help="Output path for filtered FASTQ.")
@click.option("--plot", type=click.Path(dir_okay=False), help="Path to save multi-panel QC plot.")
@click.option("--json", "json_out", type=click.Path(dir_okay=False), help="Path to save JSON metrics.")
@click.option("--backend", type=click.Choice(["auto", "cpp", "python"]), default="auto", show_default=True)
def scan(fastq_file, min_q, min_len, out_filtered, plot, json_out, backend):
    """Scan and compute QC metrics for a FASTQ/FASTQ.GZ file."""
    path = Path(fastq_file)
    backend_override = None if backend == "auto" else backend

    click.echo("=" * 80)
    click.echo("           lr-rna-qc: Oxford Nanopore RNA Quality Control Toolkit")
    click.echo("=" * 80)
    click.echo(f"Input Target  : {path.name}")
    active_backend = backend if backend != "auto" else get_backend()
    click.echo(f"Engine Backend: {active_backend}")
    click.echo(f"Filter Bounds : Min Quality = Q{min_q:.1f} | Min Read Length = {min_len:,} bp")
    click.echo("-" * 80)

    start_time = time.perf_counter()
    summary, lengths, mean_qs = scan_fastq(
        path,
        min_q=min_q,
        min_len=min_len,
        out_filtered=out_filtered,
        backend_override=backend_override,
    )
    elapsed = time.perf_counter() - start_time

    corr = compute_correlation(lengths, mean_qs)

    total_gb = summary.total_bases / 1e9
    pass_pct = (summary.pass_reads / summary.total_reads * 100.0) if summary.total_reads > 0 else 0.0
    fail_pct = (summary.fail_reads / summary.total_reads * 100.0) if summary.total_reads > 0 else 0.0
    mean_len = int(summary.total_bases / summary.total_reads) if summary.total_reads > 0 else 0

    click.echo("--------------------------- TRANSCRIPTOME QC METRICS ---------------------------")
    click.echo(f"Total Reads Scanned      : {summary.total_reads:,}")
    click.echo(f"Total Bases Sequenced    : {summary.total_bases:,} bp ({total_gb:.2f} Gbp)")
    click.echo(f"Mean Read Length         : {mean_len:,} bp")
    click.echo(f"Transcript N50           : {summary.n50:,} bp  (Length-weighted median)")
    click.echo(f"Maximum Read Length      : {summary.max_length:,} bp")
    click.echo(f"Minimum Read Length      : {summary.min_length:,} bp")
    click.echo(f"Mean Phred Quality (Q)   : {summary.mean_q:.2f}  (Exact probability calculation)")
    click.echo(f"Median Phred Quality (Q) : {summary.median_q:.2f}")

    click.echo("---------------------------- FILTERING SUMMARY ---------------------------------")
    pass_msg = f"{summary.pass_reads:,} reads ({pass_pct:.1f}%)"
    if out_filtered:
        pass_msg += f" -> saved to {out_filtered}"
    click.echo(f"Passed Quality & Length  : {pass_msg}")
    click.echo(f"Failed (Filtered Out)    : {summary.fail_reads:,} reads ({fail_pct:.1f}%)")

    click.echo("----------------------- SENSOR DEGRADATION CORRELATION -------------------------")
    click.echo(f"Pearson r  (Length vs Q) : {corr['pearson_r']:.3f} (p = {corr['pearson_p']:.2e})")
    click.echo(f"Spearman rho (Length vs Q): {corr['spearman_rho']:.3f} (p = {corr['spearman_p']:.2e})")

    if plot:
        plot_path = generate_qc_report(lengths, mean_qs, summary.n50, plot, sample_name=path.name)
        click.echo(f"[SUCCESS] Saved 3-panel publication figure -> {plot_path}")

    if json_out:
        data = {
            "sample_name": path.name,
            "backend": active_backend,
            "filters": {"min_q": min_q, "min_len": min_len},
            "metrics": summary.to_dict(),
            "correlation": corr,
            "runtime_seconds": round(elapsed, 4),
        }
        json_path = Path(json_out)
        json_path.parent.mkdir(parents=True, exist_ok=True)
        with open(json_path, "w", encoding="utf-8") as f:
            json.dump(data, f, indent=2)
        click.echo(f"[SUCCESS] Saved machine-readable metrics   -> {json_path}")

    throughput = summary.total_reads / elapsed if elapsed > 0 else 0
    click.echo(f"Completed in {elapsed:.2f}s ({throughput:,.0f} reads/s).")


@main.command()
@click.argument("sample1", type=click.Path(exists=True, dir_okay=False))
@click.argument("sample2", type=click.Path(exists=True, dir_okay=False))
@click.option("--out-json", type=click.Path(dir_okay=False), help="Path to write statistical test results.")
@click.option("--backend", type=click.Choice(["auto", "cpp", "python"]), default="auto", show_default=True)
def compare(sample1, sample2, out_json, backend):
    """Statistically compare read distributions between two flowcells or conditions."""
    p1 = Path(sample1)
    p2 = Path(sample2)
    backend_override = None if backend == "auto" else backend

    click.echo("=" * 80)
    click.echo("           lr-rna-qc: Comparative Transcriptome Hypothesis Testing")
    click.echo("=" * 80)
    click.echo(f"Sample 1: {p1.name}")
    click.echo(f"Sample 2: {p2.name}")
    click.echo("-" * 80)

    sum1, lens1, qs1 = scan_fastq(p1, backend_override=backend_override)
    sum2, lens2, qs2 = scan_fastq(p2, backend_override=backend_override)

    results = compare_distributions(lens1, qs1, lens2, qs2)
    len_comp = results["length_comparison"]
    q_comp = results["quality_comparison"]

    click.echo("------------------------- SAMPLE OVERVIEW COMPARISON ---------------------------")
    click.echo(f"Sample 1 ({p1.name}): {sum1.total_reads:,} reads | N50: {sum1.n50:,} bp | Median Q: {sum1.median_q}")
    click.echo(f"Sample 2 ({p2.name}): {sum2.total_reads:,} reads | N50: {sum2.n50:,} bp | Median Q: {sum2.median_q}")

    click.echo("----------------------- READ LENGTH DISTRIBUTION TEST --------------------------")
    len_sig = "SIGNIFICANT SHIFT" if len_comp["significant"] else "NOT SIGNIFICANT"
    click.echo(f"Mann-Whitney U Statistic : {len_comp['mann_whitney_u']:,.1f} | p-value: {len_comp['mann_whitney_p']:.2e}")
    click.echo(f"Kolmogorov-Smirnov D     : {len_comp['ks_statistic']:.4f} | p-value: {len_comp['ks_p']:.2e}")
    click.echo(f"=> Result: {len_sig} in length distribution (alpha = 0.05)")

    click.echo("---------------------- PHRED QUALITY DISTRIBUTION TEST -------------------------")
    q_sig = "SIGNIFICANT SHIFT" if q_comp["significant"] else "NOT SIGNIFICANT"
    click.echo(f"Mann-Whitney U Statistic : {q_comp['mann_whitney_u']:,.1f} | p-value: {q_comp['mann_whitney_p']:.2e}")
    click.echo(f"Kolmogorov-Smirnov D     : {q_comp['ks_statistic']:.4f} | p-value: {q_comp['ks_p']:.2e}")
    click.echo(f"=> Result: {q_sig} in quality score distribution (alpha = 0.05)")

    if out_json:
        out_data = {
            "sample1": {"name": p1.name, "summary": sum1.to_dict()},
            "sample2": {"name": p2.name, "summary": sum2.to_dict()},
            "comparison": results,
        }
        out_p = Path(out_json)
        out_p.parent.mkdir(parents=True, exist_ok=True)
        with open(out_p, "w", encoding="utf-8") as f:
            json.dump(out_data, f, indent=2)
        click.echo(f"[SUCCESS] Saved comparison statistics -> {out_p}")


@main.command()
@click.argument("fastq_file", type=click.Path(exists=True, dir_okay=False))
def benchmark(fastq_file):
    """Benchmark C++ vs pure Python backend throughput on a FASTQ file."""
    path = Path(fastq_file)
    file_size_mb = path.stat().st_size / (1024 * 1024)

    click.echo("=" * 80)
    click.echo("             lr-rna-qc: Backend Engine Performance Benchmark")
    click.echo("=" * 80)
    click.echo(f"Target File: {path.name} ({file_size_mb:.2f} MB)")
    click.echo("-" * 80)

    # 1. Benchmark Pure-Python Backend
    click.echo("Running Pure-Python Backend...")
    t0_py = time.perf_counter()
    sum_py, _, _ = scan_fastq(path, backend_override="python")
    t_py = time.perf_counter() - t0_py
    py_tput_reads = sum_py.total_reads / t_py if t_py > 0 else 0
    py_tput_mb = file_size_mb / t_py if t_py > 0 else 0

    # 2. Benchmark C++ Accelerated Backend
    has_cpp = (get_backend() == "C++")
    if not has_cpp:
        click.echo("[WARNING] C++ backend (_fast_qc) is not available on this platform.")
        click.echo(f"Python Throughput: {py_tput_reads:,.0f} reads/s ({py_tput_mb:.2f} MB/s)")
        return

    click.echo("Running C++ Accelerated Backend...")
    t0_cpp = time.perf_counter()
    sum_cpp, _, _ = scan_fastq(path, backend_override="cpp")
    t_cpp = time.perf_counter() - t0_cpp
    cpp_tput_reads = sum_cpp.total_reads / t_cpp if t_cpp > 0 else 0
    cpp_tput_mb = file_size_mb / t_cpp if t_cpp > 0 else 0

    speedup = t_py / t_cpp if t_cpp > 0 else 1.0

    click.echo("\n--------------------------- BENCHMARK RESULTS ----------------------------------")
    click.echo(f"{'Backend':<12} | {'Time (s)':<10} | {'Throughput (reads/s)':<22} | {'Throughput (MB/s)':<18}")
    click.echo("-" * 72)
    click.echo(f"{'Python':<12} | {t_py:<10.3f} | {py_tput_reads:<22,.0f} | {py_tput_mb:<18.2f}")
    click.echo(f"{'C++':<12} | {t_cpp:<10.3f} | {cpp_tput_reads:<22,.0f} | {cpp_tput_mb:<18.2f}")
    click.echo("-" * 72)
    click.echo(f"Acceleration Speedup: {speedup:.2f}x faster execution with C++ backend")

    # Verify numerical agreement
    q_diff = abs(sum_cpp.mean_q - sum_py.mean_q)
    n50_diff = abs(sum_cpp.n50 - sum_py.n50)
    parity_ok = (q_diff < 0.05 and n50_diff == 0)
    click.echo(f"Numerical Parity    : {'PASSED (Concordant Q and N50)' if parity_ok else 'FAILED'}")


if __name__ == "__main__":
    main()
