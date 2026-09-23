from pathlib import Path
from typing import Sequence
import matplotlib.pyplot as plt
import numpy as np
from lr_rna_qc.stats import compute_correlation


def generate_qc_report(
    lengths: Sequence[int],
    mean_qs: Sequence[float],
    n50: int,
    output_path: str | Path,
    sample_name: str = "Nanopore RNA Sample",
) -> Path:
    """Generate a publication-ready 3-panel QC summary figure."""
    out = Path(output_path)
    out.parent.mkdir(parents=True, exist_ok=True)

    if not lengths:
        fig, ax = plt.subplots(figsize=(8, 4))
        ax.text(0.5, 0.5, "No reads available to plot", ha="center", va="center")
        plt.tight_layout()
        plt.savefig(out, dpi=300)
        plt.close(fig)
        return out

    lens_arr = np.array(lengths, dtype=np.int64)
    qs_arr = np.array(mean_qs, dtype=np.float64)

    fig, axes = plt.subplots(1, 3, figsize=(18, 5.2), dpi=300)
    fig.suptitle(f"Transcriptome QC Report: {sample_name}", fontsize=14, fontweight="bold", y=0.98)


    # Panel 1: Read length distribution and N50
    ax1 = axes[0]
    min_l = max(10, int(np.min(lens_arr)))
    max_l = int(np.max(lens_arr))
    if max_l > min_l:
        bins = np.logspace(np.log10(min_l), np.log10(max_l), num=45)
    else:
        bins = 30

    ax1.hist(lens_arr, bins=bins, color="#2b5c8f", edgecolor="#1a365d", alpha=0.85)
    ax1.set_xscale("log")
    ax1.axvline(n50, color="#d9381e", linestyle="--", linewidth=2, label=f"N50: {n50:,} bp")
    ax1.set_xlabel("Read Length (bp, log scale)", fontsize=11)
    ax1.set_ylabel("Read Count", fontsize=11)
    ax1.set_title("Read Length Distribution", fontsize=12, fontweight="semibold")
    ax1.legend(loc="upper right", frameon=True)
    ax1.grid(True, which="both", linestyle=":", alpha=0.5)


    # Panel 2: Cumulative Base Yield vs. Read Length
    ax2 = axes[1]
    sorted_idx = np.argsort(lens_arr)[::-1]
    sorted_lens = lens_arr[sorted_idx]
    cum_yield_gb = np.cumsum(sorted_lens) / 1e9
    total_gb = cum_yield_gb[-1]
    cum_pct = (cum_yield_gb / total_gb) * 100.0

    ax2.plot(sorted_lens, cum_pct, color="#2e7d32", linewidth=2.2, label="Cumulative Yield (%)")
    ax2.set_xscale("log")
    ax2.invert_xaxis()  # From longest read down to shortest

    # Mark N50 threshold (50% yield)
    ax2.axhline(50.0, color="#d9381e", linestyle=":", linewidth=1.5)
    ax2.axvline(n50, color="#d9381e", linestyle="--", linewidth=1.5, label=f"N50 (50% Yield: {n50:,} bp)")
    ax2.plot(n50, 50.0, marker="o", markersize=7, color="#d9381e")

    ax2.set_xlabel("Read Length Threshold (bp, descending)", fontsize=11)
    ax2.set_ylabel("Cumulative Base Yield (%)", fontsize=11)
    ax2.set_title(f"Yield Accumulation (Total: {total_gb:.2f} Gbp)", fontsize=12, fontweight="semibold")
    ax2.legend(loc="lower left", frameon=True)
    ax2.grid(True, which="both", linestyle=":", alpha=0.5)

    # Panel 3: Quality vs. Read Length with Trend Line
    ax3 = axes[2]
    if len(lens_arr) > 1500:
        hb = ax3.hexbin(
            lens_arr, qs_arr,
            gridsize=40,
            xscale="log",
            cmap="viridis",
            mincnt=1,
            edgecolors="none"
        )
        cbar = fig.colorbar(hb, ax=ax3, fraction=0.046, pad=0.04)
        cbar.set_label("Read Count", fontsize=10)
    else:
        ax3.scatter(lens_arr, qs_arr, alpha=0.5, color="#6b46c1", edgecolors="none", s=25)
        ax3.set_xscale("log")

    # Fit and display linear regression trend line to evaluate sensor decay
    corr = compute_correlation(lens_arr, qs_arr)
    if len(lens_arr) >= 2 and max_l > min_l:
        # Fit in log10 space for length
        log_lens = np.log10(lens_arr)
        slope, intercept = np.polyfit(log_lens, qs_arr, 1)
        x_fit = np.logspace(np.log10(min_l), np.log10(max_l), 100)
        y_fit = slope * np.log10(x_fit) + intercept
        ax3.plot(x_fit, y_fit, color="#e53e3e", linewidth=2.0, linestyle="-", label="Trend Line")

    # Inset stats box
    stats_text = (
        f"Pearson r: {corr['pearson_r']:.3f} (p={corr['pearson_p']:.2e})\n"
        f"Spearman \u03c1: {corr['spearman_rho']:.3f} (p={corr['spearman_p']:.2e})"
    )
    ax3.text(
        0.05, 0.05, stats_text,
        transform=ax3.transAxes,
        fontsize=9,
        verticalalignment="bottom",
        bbox=dict(boxstyle="round,pad=0.4", facecolor="white", alpha=0.85, edgecolor="#cbd5e0")
    )

    ax3.set_xlabel("Read Length (bp, log scale)", fontsize=11)
    ax3.set_ylabel("Mean Phred Quality (Q)", fontsize=11)
    ax3.set_title("Quality vs. Length (Sensor Decay)", fontsize=12, fontweight="semibold")
    ax3.grid(True, which="both", linestyle=":", alpha=0.5)

    plt.tight_layout()
    plt.savefig(out, dpi=300)
    plt.close(fig)
    return out
