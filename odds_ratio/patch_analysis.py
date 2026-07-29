"""
Patch analysis for odds ratio computation.

Handles the analysis of structural patches/pockets for hit enrichment.
"""

import logging
from dataclasses import dataclass
from typing import Dict, List, Optional
import pandas as pd

from .config import (
    MIN_PATCH_RESIDUES,
    MAX_PATCH_RESIDUES,
    get_processed_files_dirpath,
)
from .structure import StructureData, process_structure, find_residues_in_radius
from .statistics import compute_odds_ratio, apply_multiple_testing_correction

log = logging.getLogger(__name__)


@dataclass
class PatchResult:
    """Container for patch analysis results."""
    df: pd.DataFrame
    uniprot_id: str


def get_common_hit_flag(count: float, max_n_hits: float) -> Optional[bool]:
    """Determine if a residue is a common hit across datasets.
    
    Args:
        count: Number of hits for this residue
        max_n_hits: Maximum possible hits (number of datasets)
        
    Returns:
        True if common hit, False if no hit, None otherwise
    """
    if count >= max_n_hits:
        return True
    elif count == 0:
        return False
    else:
        return None  # Hit in some but not all cell lines


def analyze_patches_for_chain(
    df: pd.DataFrame,
    structure_data: StructureData,
    chain: str,
    uniprot_id: str,
    radius: float,
    max_n_hits: float,
    min_patch_residues: int = MIN_PATCH_RESIDUES,
    max_patch_residues: int = MAX_PATCH_RESIDUES,
) -> pd.DataFrame:
    """Analyze all residue-centered patches for a single chain.
    
    Args:
        df: LFC data for this chain
        structure_data: Processed structure data
        chain: Chain identifier
        uniprot_id: UniProt ID
        radius: Patch radius in Angstroms
        max_n_hits: Maximum possible hits
        min_patch_residues: Minimum residues required in patch
        max_patch_residues: Maximum residues allowed in patch
        
    Returns:
        DataFrame with odds ratios for all patches
    """
    # Apply common hit classification
    df = df.copy()
    df["common_hit"] = df["n_hits"].apply(lambda x: get_common_hit_flag(x, max_n_hits))
    df = df[df["common_hit"].isin([True, False])]
    
    # Keep only the highest common_hit value per position
    df = df.sort_values("common_hit", ascending=False).drop_duplicates(["Gene", "rno"])
    
    protein_rnos = structure_data.protein_rnos
    
    # Skip if too few residues
    if len(protein_rnos) <= 30:
        log.warning(f"Chain {chain} has too few residues ({len(protein_rnos)}), skipping")
        return pd.DataFrame()
    
    all_rows = []
    
    for resi, rno in enumerate(protein_rnos):
        if rno % 100 == 0:
            log.debug(f"Processing residue {rno}")
        
        # Find residues in patch
        pocket_rnos = find_residues_in_radius(structure_data, resi, radius)
        
        if len(pocket_rnos) > max_patch_residues:
            continue
        
        # Check minimum guides in patch
        in_pocket = df["rno"].isin(pocket_rnos)
        site = df[in_pocket]
        n_guides = site["rno"].nunique()
        
        # Require minimum guides and at least 2 high-impact residues
        high_impact_count = site[site["lt_cutoff"]]["rno"].nunique()
        
        if n_guides < min_patch_residues or high_impact_count < 2:
            continue
        
        # Compute odds ratio
        or_result = compute_odds_ratio(df, pocket_rnos, rno, uniprot_id)
        
        # Build result rows (one per residue in patch)
        for pocket_rno in pocket_rnos:
            row = {
                "radius": radius,
                "central_residue": rno,
                "pocket_residues": pocket_rno,
                "OR": or_result.odds_ratio,
                "pval": or_result.p_value,
                "low": or_result.ci_lower,
                "up": or_result.ci_upper,
                "ks_location": or_result.ks_location,
                "site_of_interest_high_impact": or_result.contingency_table["high_impact"]["site of interest"],
                "site_of_interest_low_impact": or_result.contingency_table["low_impact"]["site of interest"],
                "rest_high_impact": or_result.contingency_table["high_impact"]["rest"],
                "rest_low_impact": or_result.contingency_table["low_impact"]["rest"],
            }
            all_rows.append(row)
    
    if not all_rows:
        return pd.DataFrame()
    
    return pd.DataFrame(all_rows)


def analyze_all_patches(
    df: pd.DataFrame,
    pdb_code: str,
    radius: float,
    max_n_hits: float,
) -> Dict[str, PatchResult]:
    """Analyze patches for all chains in a structure.
    
    Args:
        df: LFC data with chain and residue information
        pdb_code: PDB code or UniProt ID
        radius: Patch radius in Angstroms
        max_n_hits: Maximum possible hits
        
    Returns:
        Dictionary mapping chain IDs to PatchResult objects
    """
    chain_uniprot_df = df[["chain", "uniprot_id"]].drop_duplicates()
    chain_uniprot_df = chain_uniprot_df.sort_values("chain")
    
    all_results = {}
    
    for _, row in chain_uniprot_df.iterrows():
        chain = row["chain"]
        uniprot_id = row["uniprot_id"]
        
        log.info(f"Processing chain {chain} (UniProt: {uniprot_id})")
        
        try:
            structure_data = process_structure(pdb_code, chain)
        except Exception as e:
            log.error(f"Failed to process structure for chain {chain}: {e}")
            continue
        
        chain_df = df[df["chain"] == chain]
        
        patch_df = analyze_patches_for_chain(
            chain_df,
            structure_data,
            chain,
            uniprot_id,
            radius,
            max_n_hits,
        )
        
        if len(patch_df) == 0:
            continue
        
        # Apply multiple testing correction
        corrected_df = apply_multiple_testing_correction(patch_df)
        
        all_results[chain] = PatchResult(df=corrected_df, uniprot_id=uniprot_id)
    
    return all_results


def postprocess_results(
    results: Dict[str, PatchResult],
    pdb_code: str,
) -> Dict[str, PatchResult]:
    """Apply PDB-specific postprocessing to results.
    
    Args:
        results: Dictionary of chain results
        pdb_code: PDB code
        
    Returns:
        Filtered results dictionary
    """
    # Handle special cases (e.g., 6bcu has specific chains of interest)
    if pdb_code == "6bcu":
        valid_chains = {"A", "D", "W", "S"}
        results = {
            chain: result 
            for chain, result in results.items() 
            if chain in valid_chains
        }
    
    return results


def save_intermediate_results(
    results: Dict[str, PatchResult],
    pdb_code: str,
    datasets_desc: List[str],
    mutation_type: str,
) -> None:
    """Save intermediate results to processed_data_files.
    
    Args:
        results: Dictionary of chain results
        pdb_code: PDB code
        datasets_desc: List of dataset descriptions
        mutation_type: Mutation type analyzed
    """
    output_dir = get_processed_files_dirpath() / "odds_ratio_results"
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Combine all chain results
    all_dfs = []
    for chain, result in results.items():
        chain_df = result.df.copy()
        chain_df["chain"] = chain
        chain_df["uniprot_id"] = result.uniprot_id
        all_dfs.append(chain_df)
    
    if not all_dfs:
        return
    
    combined_df = pd.concat(all_dfs, axis=0)
    combined_df = combined_df.sort_values(["chain", "central_residue", "pocket_residues"])
    combined_df = combined_df.round(5)
    
    # Generate filename from parameters
    dataset_hash = "_".join(sorted(datasets_desc))[:50]  # Truncate long names
    filename = f"odds_ratio_{pdb_code}_{mutation_type}.csv"
    
    output_path = output_dir / filename
    combined_df.to_csv(output_path, index=False)
    log.info(f"Saved odds ratio results to {output_path}")
