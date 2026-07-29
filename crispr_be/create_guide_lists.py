#!/usr/bin/env python
"""Create guide lists by merging used guides with all possible edit combinations."""

from pathlib import Path

import hydra
import pandas as pd
from omegaconf import DictConfig

# Get absolute path to config directory
CONFIG_DIR = Path(__file__).parent.parent / "config"


def select_gene(gene: str, input_filepath: str | Path) -> pd.DataFrame:
    """Select guides for a specific gene from the input file.

    Args:
        gene: Gene symbol to filter for.
        input_filepath: Path to the CSV file containing guide data.

    Returns:
        DataFrame containing guides for the specified gene.
    """
    str_columns = ["guide_in_CDS"]
    used_guides = pd.read_csv(input_filepath, dtype={col: str for col in str_columns})

    # Standardize column names
    used_guides.rename(
        columns={
            "CRISPR+PAM Sequence": "CRISPR_PAM_Sequence",
            "gRNA Target Sequence": "gRNA_Target_Sequence",
        },
        inplace=True,
    )

    # Fix gene name inconsistencies
    gene_replacements = {"PRR5l": "PRR5L", "RHEB2": "RHEB"}
    for wrong_gene, correct_gene in gene_replacements.items():
        used_guides["Hugo_Symbol"] = used_guides["Hugo_Symbol"].str.replace(
            wrong_gene, correct_gene, regex=False
        )
        used_guides["gRNA_ID"] = used_guides["gRNA_ID"].str.replace(
            wrong_gene, correct_gene, regex=False
        )

    gene_selection = used_guides[used_guides["Hugo_Symbol"] == gene]
    return gene_selection


def merge_guides(config: DictConfig) -> None:
    """Merge used guides with all possible edit combinations.

    For each gene and editor combination, merges the guide data with
    BEstimate edit predictions and saves the result.

    Args:
        config: Hydra configuration containing paths and parameters.
    """
    target_genes = config["target_genes"]
    editors = config["editors"]
    guide_data_filepath = config["guide_data_filepath"]
    bestimate_data_dirpath = Path(config["bestimate_data_dirpath"])
    guide_with_edits_dirpath = Path(config["guide_with_edits_dirpath"])

    # Ensure output directory exists
    guide_with_edits_dirpath.mkdir(parents=True, exist_ok=True)

    for gene in target_genes:
        for editor in editors:
            print(f"Processing {gene} with {editor}...")

            # Select guides for this gene
            used_guides = select_gene(gene, guide_data_filepath)

            # Load BEstimate edit predictions
            bestimate_filepath = bestimate_data_dirpath / f"{gene}_{editor}_edit_df.csv"
            if not bestimate_filepath.exists():
                print(f"  Warning: {bestimate_filepath} not found, skipping.")
                continue

            all_edits_df = pd.read_csv(bestimate_filepath)

            # Standardize column names for merging
            used_guides.rename(
                columns={
                    "CRISPR+PAM Sequence": "CRISPR_PAM_Sequence",
                    "gRNA Target Sequence": "gRNA_Target_Sequence",
                    "gRNA_ID": "sgRNA_ID",
                    "Transcript_ID": "Canonical_transcript",
                },
                inplace=True,
            )

            # Merge on common columns
            merged = pd.merge(
                all_edits_df,
                used_guides,
                on=[
                    "CRISPR_PAM_Sequence",
                    "Direction",
                    "Gene_ID",
                    "gRNA_Target_Sequence",
                    "Location",
                    "Hugo_Symbol",
                ],
                how="inner",
            )

            # Sanity check
            n_used = len(used_guides["gRNA_Target_Sequence"].unique())
            n_edits = len(all_edits_df["gRNA_Target_Sequence"].unique())
            n_merged = len(merged["gRNA_Target_Sequence"].unique())
            print(f"  Used guides: {n_used}, All edits: {n_edits}, Merged: {n_merged}")

            if n_merged == 0:
                print(f"  Warning: No guides merged for {gene} {editor}!")
                continue

            # Save merged data
            output_filepath = (
                guide_with_edits_dirpath / f"{gene}_{editor}_guides_w_edits.csv"
            )
            merged.to_csv(output_filepath, index=False)
            print(f"  Saved to {output_filepath}")


@hydra.main(version_base=None, config_path=str(CONFIG_DIR), config_name="config")
def main(cfg: DictConfig) -> None:
    """Main entry point."""
    merge_guides(cfg)


if __name__ == "__main__":
    main()
