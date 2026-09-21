from typing import Sequence
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np


def generate_qc_report(
    lengths: Sequence[int],
    mean_qs: Sequence[float],
    n50: int,
    output_path: str | Path,
    sample_name: str = "Nanopore RNA Sample",
) -> None:
    """Generate a 3-panel publication-ready QC summary figure."""
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    fig.suptitle(f"Transcriptome QC Report: {sample_name}", fontsize=14, fontweight="bold")

    # Panel 1: Read length distribution (log scale)
    ax1 = axes[0]
    ax1.hist(lengths, bins=50, color="#1f77b4", alpha=0.75, edgecolor="black")
    ax1.axvline(n50, color="red", linestyle="--", linewidth=2, label=f"N50: {n50:,} bp")
    ax1.set_xlabel("Read Length (bp)")
    ax1.set_ylabel("Read Count")
    ax1.set_title("Read Length Distribution")
    ax1.legend()

    # Panel 2: Cumulative Base Yield vs Read Length
    ax2 = axes[1]
    sorted_lens = np.sort(lengths)[::-1]
    cum_bases = np.cumsum(sorted_lens) / 1e9
    ax2.plot(range(1, len(sorted_lens) + 1), cum_bases, color="#2ca02c", linewidth=2)
    ax2.set_xlabel("Cumulative Read Count (Sorted Descending)")
    ax2.set_ylabel("Cumulative Yield (Gigabases)")
    ax2.set_title("Yield Accumulation")
    ax2.grid(True, linestyle=":", alpha=0.6)

    # Panel 3: Quality vs Length Scatter / Hexbin
    ax3 = axes[2]
    if len(lengths) > 2000:
        hb = ax3.hexbin(lengths, mean_qs, gridsize=40, cmap="viridis", mincnt=1)
        fig.colorbar(hb, ax=ax3, label="Count")
    else:
        ax3.scatter(lengths, mean_qs, alpha=0.4, color="#9467bd", edgecolors="none")
    ax3.set_xlabel("Read Length (bp)")
    ax3.set_ylabel("Mean Phred Quality (Q)")
    ax3.set_title("Quality vs. Read Length")
    ax3.grid(True, linestyle=":", alpha=0.6)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close(fig)
