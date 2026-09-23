from typing import Dict, Sequence, Any
from scipy import stats


def compute_correlation(lengths: Sequence[int], mean_qs: Sequence[float]) -> Dict[str, float]:
    """Compute Pearson and Spearman correlation between read length and mean quality."""
    if len(lengths) < 2 or len(set(lengths)) <= 1 or len(set(mean_qs)) <= 1:
        return {
            "pearson_r": 0.0,
            "pearson_p": 1.0,
            "spearman_rho": 0.0,
            "spearman_p": 1.0,
        }

    p_r, p_p = stats.pearsonr(lengths, mean_qs)
    s_r, s_p = stats.spearmanr(lengths, mean_qs)

    return {
        "pearson_r": round(float(p_r), 4),
        "pearson_p": float(p_p),
        "spearman_rho": round(float(s_r), 4),
        "spearman_p": float(s_p),
    }


def compare_distributions(
    lens1: Sequence[int],
    qs1: Sequence[float],
    lens2: Sequence[int],
    qs2: Sequence[float],
) -> Dict[str, Any]:
    """Perform Mann-Whitney U and Kolmogorov-Smirnov tests comparing two flowcells or conditions."""
    if len(lens1) < 2 or len(lens2) < 2:
        return {
            "length_comparison": {
                "mann_whitney_u": 0.0,
                "mann_whitney_p": 1.0,
                "ks_statistic": 0.0,
                "ks_p": 1.0,
                "significant": False,
            },
            "quality_comparison": {
                "mann_whitney_u": 0.0,
                "mann_whitney_p": 1.0,
                "ks_statistic": 0.0,
                "ks_p": 1.0,
                "significant": False,
            },
        }

    mw_len = stats.mannwhitneyu(lens1, lens2, alternative="two-sided")
    ks_len = stats.ks_2samp(lens1, lens2)

    mw_q = stats.mannwhitneyu(qs1, qs2, alternative="two-sided")
    ks_q = stats.ks_2samp(qs1, qs2)

    return {
        "length_comparison": {
            "mann_whitney_u": float(mw_len.statistic),
            "mann_whitney_p": float(mw_len.pvalue),
            "ks_statistic": round(float(ks_len.statistic), 4),
            "ks_p": float(ks_len.pvalue),
            "significant": bool(mw_len.pvalue < 0.05 or ks_len.pvalue < 0.05),
        },
        "quality_comparison": {
            "mann_whitney_u": float(mw_q.statistic),
            "mann_whitney_p": float(mw_q.pvalue),
            "ks_statistic": round(float(ks_q.statistic), 4),
            "ks_p": float(ks_q.pvalue),
            "significant": bool(mw_q.pvalue < 0.05 or ks_q.pvalue < 0.05),
        },
    }
