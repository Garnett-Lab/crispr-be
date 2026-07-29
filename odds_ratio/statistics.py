"""
Statistical computations for odds ratio analysis.

Provides Fisher's exact test, confidence intervals, and multiple testing correction.
"""

from dataclasses import dataclass
from typing import List, Tuple

import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.stats.multitest as smm

from .config import MAX_ODDS_RATIO


@dataclass
class OddsRatioResult:
    """Container for odds ratio computation results."""
    odds_ratio: float
    p_value: float
    ci_lower: float
    ci_upper: float
    contingency_table: pd.DataFrame
    ks_location: float


def compute_confidence_interval_95(
    odds_ratio: float,
    contingency_table: pd.DataFrame,
) -> Tuple[float, float]:
    """Compute 95% confidence interval for odds ratio.
    
    Args:
        odds_ratio: Computed odds ratio
        contingency_table: 2x2 contingency table with counts
        
    Returns:
        Tuple of (lower, upper) confidence bounds
    """
    # Standard error using log-odds method
    se = np.sqrt(
        (1 / contingency_table["high_impact"]["site of interest"])
        + (1 / contingency_table["high_impact"]["rest"])
        + (1 / contingency_table["low_impact"]["site of interest"])
        + (1 / contingency_table["low_impact"]["rest"])
    )
    
    # 95% CI bounds
    log_or = np.log(odds_ratio)
    lower = np.exp(log_or - 1.96 * se)
    upper = np.exp(log_or + 1.96 * se)
    
    return lower, upper


def compute_odds_ratio(
    df: pd.DataFrame,
    pocket_rnos: List[int],
    center_rno: int,
    uniprot_id: str,
    max_or: float = MAX_ODDS_RATIO,
) -> OddsRatioResult:
    """Compute odds ratio for hit enrichment in a structural patch.
    
    Compares the proportion of significant hits in a patch vs. rest of structure.
    
    Args:
        df: DataFrame with LFC data, must have 'rno' and 'lt_cutoff' columns
        pocket_rnos: List of residue numbers in the patch
        center_rno: Center residue number (for logging)
        uniprot_id: UniProt ID (for logging)
        max_or: Maximum odds ratio (used when cells are zero)
        
    Returns:
        OddsRatioResult with odds ratio, p-value, CI, and contingency table
    """
    # Split data into patch (site of interest) and rest
    in_pocket = df["rno"].isin(pocket_rnos)
    site_df = df[in_pocket]
    rest_df = df[~in_pocket]
    
    # Count hits (lt_cutoff=True) and non-hits by unique residue
    selected_columns = ["rno"]
    
    high_impact_site = site_df[site_df["lt_cutoff"]][selected_columns].drop_duplicates().shape[0]
    high_impact_rest = rest_df[rest_df["lt_cutoff"]][selected_columns].drop_duplicates().shape[0]
    low_impact_site = site_df[~site_df["lt_cutoff"]][selected_columns].drop_duplicates().shape[0]
    low_impact_rest = rest_df[~rest_df["lt_cutoff"]][selected_columns].drop_duplicates().shape[0]
    
    # Build contingency table
    contingency = [[high_impact_site, high_impact_rest], 
                   [low_impact_site, low_impact_rest]]
    
    contingency_df = pd.DataFrame(
        dict(zip(["high_impact", "low_impact"], contingency))
    ).rename(index={0: "site of interest", 1: "rest"})
    
    # Fisher's exact test (one-sided: greater)
    odds_ratio, p_value = stats.fisher_exact(contingency, alternative="greater")
    
    # Kolmogorov-Smirnov test for LFC distributions
    site_lfcs = site_df.drop_duplicates(selected_columns)["LFC"].values
    other_lfcs = rest_df.drop_duplicates(selected_columns)["LFC"].values
    
    if len(site_lfcs) > 0 and len(other_lfcs) > 0:
        ks_result = stats.ks_2samp(site_lfcs, other_lfcs)
        # Find location of max KS statistic
        all_values = np.sort(np.concatenate([site_lfcs, other_lfcs]))
        cdf1 = np.searchsorted(np.sort(site_lfcs), all_values, side='right') / len(site_lfcs)
        cdf2 = np.searchsorted(np.sort(other_lfcs), all_values, side='right') / len(other_lfcs)
        ks_location = all_values[np.argmax(np.abs(cdf1 - cdf2))]
    else:
        ks_location = 0.0
    
    # Handle zero cells in contingency table
    has_zero = any(v == 0 for v in [high_impact_site, high_impact_rest, 
                                      low_impact_site, low_impact_rest])
    
    if has_zero:
        ci_lower = 0.0
        ci_upper = 0.0
        odds_ratio = max_or
    else:
        ci_lower, ci_upper = compute_confidence_interval_95(odds_ratio, contingency_df)
    
    return OddsRatioResult(
        odds_ratio=odds_ratio,
        p_value=p_value,
        ci_lower=ci_lower,
        ci_upper=ci_upper,
        contingency_table=contingency_df,
        ks_location=ks_location,
    )


def apply_multiple_testing_correction(
    df: pd.DataFrame,
    method: str = "fdr_bh",
    alpha: float = 0.05,
) -> pd.DataFrame:
    """Apply multiple testing correction to p-values.
    
    Args:
        df: DataFrame with 'pval', 'central_residue', 'radius' columns
        method: Correction method (default: Benjamini-Hochberg)
        alpha: Significance level
        
    Returns:
        DataFrame with added 'bonf_pval' column
    """
    # Get unique combinations for correction
    subset_df = df.drop_duplicates(['pval', 'central_residue', 'radius'])
    
    # Apply correction
    _, corrected_pvals, _, _ = smm.multipletests(
        subset_df['pval'], 
        method=method, 
        alpha=alpha
    )
    
    # Create correction lookup
    correction_df = pd.DataFrame({
        'radius': subset_df['radius'].values,
        'central_residue': subset_df['central_residue'].values,
        'bonf_pval': corrected_pvals,
    })
    
    # Merge back to original data
    result = pd.merge(df, correction_df, on=['central_residue', 'radius'], how='inner')
    
    return result
