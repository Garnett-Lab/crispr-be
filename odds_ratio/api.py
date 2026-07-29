"""
Public API for odds ratio analysis.

Provides the main entry point for computing odds ratios via the web API.
"""

import logging
from typing import Any, Dict, List

import pandas as pd

from .fetch_lfc_values import add_n_hits, load_data, preprocess_df
from .config import CORRECTED_PVALUE_THRESHOLD, DEFAULT_PATCH_CENTER, DEFAULT_RADIUS
from .patch_analysis import (
    analyze_all_patches,
    postprocess_results,
    save_intermediate_results,
)

log = logging.getLogger(__name__)


def calculate_odds_ratios(
    datasets_desc: List[str],
    mutation_type: str,
    pdb_code: str,
    is_sig: bool = False,
    remove_ot: bool = True,
    radius: float = DEFAULT_RADIUS,
    patch_center: str = DEFAULT_PATCH_CENTER,
) -> Dict[str, Any]:
    """Calculate odds ratios for all patches in a structure.
    
    Args:
        datasets_desc: List of dataset descriptions (cell_line_editor format)
        mutation_type: Type of mutation to analyze (e.g., 'missense_variant')
        pdb_code: PDB code or UniProt ID
        is_sig: Filter to significant hits only
        remove_ot: Remove off-target hits
        radius: Patch radius in Angstroms
        patch_center: Center type ('residue', 'fpocket', 'ligand')
        
    Returns:
        Dictionary mapping chains to PatchResult objects
    """
    # Load and prepare data
    all_data, guide_level_data = load_data(pdb_code)
    
    # Parse dataset descriptions to filter by cell line and editor
    description_rows = []
    for description in datasets_desc:
        parts = description.split("_")
        cell_line = parts[2]
        editor = parts[4][:3]
        description_rows.append({"cell_line": cell_line, "editor": editor})
    
    description_df = pd.DataFrame(description_rows)
    df = all_data.merge(description_df, on=["cell_line", "editor"])
    guide_level_data = guide_level_data.merge(description_df, on=["cell_line", "editor"])
    
    # Preprocess data
    df = preprocess_df(df, is_sig, remove_ot)
    
    # Filter by mutation type
    mutation_guides = guide_level_data[guide_level_data["main_conseq"] == mutation_type]
    df = df[
        (df["sgRNA_ID"].isin(mutation_guides["sgRNA_ID"].values))
        & (df["Consequence"] == mutation_type)
    ]
    
    # Add n_hits column
    mut_guide_level_data = guide_level_data.merge(
        mutation_guides[["sgRNA_ID", "editor"]].drop_duplicates(),
        on=["sgRNA_ID", "editor"],
    )
    df = add_n_hits(df, mutation_guides, mut_guide_level_data, mutation_type)
    
    # Maximum possible hits is 1.0 (normalized)
    max_n_hits = 1.0
    
    # Analyze all patches
    results = analyze_all_patches(df, pdb_code, radius, max_n_hits)
    
    # Apply PDB-specific postprocessing
    results = postprocess_results(results, pdb_code)
    
    return results


def get_or_formatted(
    datasets_desc: List[str],
    mutation_type: str,
    pdb_code: str,
    is_sig: bool,
    remove_ot: bool,
    radius: float,
    patch_center: str,
    center_only: bool,
) -> Dict[str, Dict[str, Any]]:
    """Get formatted odds ratio results for API response.
    
    This is the main entry point called by the Flask API.
    
    Args:
        datasets_desc: List of dataset descriptions
        mutation_type: Type of mutation to analyze
        pdb_code: PDB code or UniProt ID
        is_sig: Filter to significant hits only
        remove_ot: Remove off-target hits
        radius: Patch radius in Angstroms
        patch_center: Center type ('residue', 'fpocket', 'ligand')
        center_only: If True, only return center residue OR values
        
    Returns:
        Dictionary mapping chains to {values: {rno: OR}, uniprot_id: str}
    """
    # Calculate odds ratios
    patch_results = calculate_odds_ratios(
        datasets_desc=datasets_desc,
        mutation_type=mutation_type,
        pdb_code=pdb_code,
        is_sig=is_sig,
        remove_ot=remove_ot,
        radius=radius,
        patch_center=patch_center,
    )
    
    # Format results for API response
    final_results = {}
    
    for chain, result in patch_results.items():
        df = result.df
        uniprot_id = result.uniprot_id
        
        # Filter by center_only
        if center_only:
            column_pocket = "central_residue"
            df = df[df["central_residue"] == df["pocket_residues"]]
        else:
            column_pocket = "pocket_residues"
        
        # Apply significance filter
        stat_significant = df["bonf_pval"] < CORRECTED_PVALUE_THRESHOLD
        higher_odds_ratio = df["OR"] > 1
        
        significant_df = df[stat_significant & higher_odds_ratio]
        
        # Get maximum OR per residue
        values = significant_df.groupby(column_pocket)["OR"].max().to_dict()
        
        final_results[chain] = {
            "values": values,
            "uniprot_id": uniprot_id,
        }
    
    # Save intermediate results for later analysis
    try:
        save_intermediate_results(patch_results, pdb_code, datasets_desc, mutation_type)
    except Exception as e:
        log.warning(f"Failed to save intermediate results: {e}")
    
    return final_results
