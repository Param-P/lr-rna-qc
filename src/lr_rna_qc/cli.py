import click

@click.group()
@click.version_option(package_name="lr-rna-qc")
def main():
    """lr-rna-qc: Oxford Nanopore RNA-seq Quality Control & Splicing Analysis Tool."""
    pass


@main.command()
@click.argument("fastq_file", type=click.Path(exists=True))
@click.option("--min-q", default=10.0, help="Minimum mean Phred quality threshold.")
@click.option("--min-len", default=200, help="Minimum read length threshold.")
@click.option("--out-filtered", type=click.Path(), help="Output path for filtered FASTQ.")
@click.option("--plot", type=click.Path(), help="Path to save multi-panel QC plot.")
@click.option("--json", "json_out", type=click.Path(), help="Path to save JSON metrics.")
@click.option("--backend", type=click.Choice(["auto", "cpp", "python"]), default="auto")
def scan(fastq_file, min_q, min_len, out_filtered, plot, json_out, backend):
    """Scan and compute QC metrics for a FASTQ/FASTQ.GZ file."""
    click.echo(f"Scanning {fastq_file} with backend={backend}...")


@main.command()
@click.argument("sample1", type=click.Path(exists=True))
@click.argument("sample2", type=click.Path(exists=True))
@click.option("--out-json", type=click.Path(), help="Path to write statistical test results.")
def compare(sample1, sample2, out_json):
    """Statistically compare read distributions between two flowcells or conditions."""
    click.echo(f"Comparing {sample1} and {sample2}...")


@main.command()
@click.argument("fastq_file", type=click.Path(exists=True))
def benchmark(fastq_file):
    """Benchmark C++ vs pure Python backend throughput on a FASTQ file."""
    click.echo(f"Benchmarking C++ vs Python on {fastq_file}...")


if __name__ == "__main__":
    main()
