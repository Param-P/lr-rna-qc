# lr-rna-qc

A production-grade, hybrid C++/Python CLI tool designed to inspect, benchmark, and quality-filter raw Oxford Nanopore FASTQ files for alternative RNA splicing analysis (tailored for the Brooks Lab, UC Santa Cruz).

## Core Architectural Goals

1. **High-Performance C++ Core**: A fast C++ extension (bound via `pybind11`) with a precomputed 256-element Look-Up Table (LUT) that converts ASCII-33 Phred quality characters to error probabilities ($P = 10^{-Q/10}$), computes the mathematically exact mean quality score ($Q_{\text{mean}} = -10 \log_{10} \bar{P}$), and extracts read lengths without heap thrashing.
2. **Pure-Python Fallback**: A seamless pure-Python fallback implementation so the package installs and runs without error even if a C++ compilation toolchain is missing.
3. **Transcriptome-Specific Metrics**: Exact computation of N50 (length-weighted median), total gigabases, read length distributions, and pass/fail filter thresholds.
4. **Statistical Analysis Layer**:
   - Correlation: Spearman/Pearson correlation between read length and mean quality score (testing sensor degradation over long molecules).
   - Comparative Hypothesis Testing: `compare` subcommand executing a Mann-Whitney U test and two-sample Kolmogorov-Smirnov test to evaluate distributional shift between two flowcells or experimental conditions.
5. **Automated Reporting**: A publication-ready multi-panel visualization (matplotlib/seaborn) displaying length distributions, cumulative volume curves indicating the N50 cutoff, and quality heatmaps.
6. **Zero-Config Packaging**: Modern `pyproject.toml` build system using `setuptools` and `pybind11`.

## Project Structure

```text
lr-rna-qc/
├── setup.py
├── pyproject.toml
├── README.md
├── src/
│   ├── cpp/
│   │   ├── fast_qc.hpp
│   │   ├── fast_qc.cpp
│   │   └── bindings.cpp
│   └── lr_rna_qc/
│       ├── __init__.py
│       ├── core.py         # Handles C++ loading with fallback
│       ├── fallback.py     # Pure Python implementation
│       ├── metrics.py      # N50 and throughput logic
│       ├── stats.py        # Statistical tests (scipy.stats)
│       ├── plot.py         # Visual summary generator
│       └── cli.py          # Click CLI entry point
└── tests/
    ├── test_core.py
    ├── test_metrics.py
    └── make_synthetic_fastq.py
```

## Quickstart

```bash
# Install with uv
uv pip install -e .

# Run test suite
uv run pytest -v

# CLI usage
uv run lr-rna-qc --help
```
