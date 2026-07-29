"""Compute gene to ENSG mapping from guide data."""

import json
import logging
from pathlib import Path

import hydra
import pandas as pd
from omegaconf import DictConfig

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Get absolute path to config directory
CONFIG_DIR = Path(__file__).parent.parent / "config"


def correct_gene_names(df: pd.DataFrame, corrections: dict[str, str]) -> pd.DataFrame:
    """Apply gene name corrections using exact matching."""
    df = df.copy()
    df["Hugo_Symbol"] = df["Hugo_Symbol"].replace(corrections)
    df["gRNA_ID"] = df["gRNA_ID"].replace(corrections, regex=True)
    return df


def compute_gene2ensg(guides_df: pd.DataFrame) -> dict[str, str]:
    """Compute mapping from gene symbol to ENSG identifier."""
    gene2ensg = (
        guides_df[["Hugo_Symbol", "Gene_ID"]]
        .drop_duplicates()
        .set_index("Hugo_Symbol")["Gene_ID"]
        .to_dict()
    )
    return gene2ensg


def compute_gene2enst(guides_df: pd.DataFrame) -> dict[str, str]:
    """Compute mapping from gene symbol to ENST identifier."""
    gene2enst = (
        guides_df[["Hugo_Symbol", "Transcript_ID"]]
        .drop_duplicates()
        .set_index("Hugo_Symbol")["Transcript_ID"]
        .to_dict()
    )
    return gene2enst


@hydra.main(version_base=None, config_path=str(CONFIG_DIR), config_name="config")
def main(cfg: DictConfig) -> None:
    """Main entry point.

    Args:
        cfg: Hydra configuration.
    """
    guides_df_path = Path(cfg.guide_data_filepath)
    gene2ensg_path = Path(cfg.gene2ensg_filepath)
    gene2enst_path = Path(cfg.gene2enst_filepath)
    target_genes = frozenset(cfg.target_genes)
    gene_name_corrections = cfg.gene_name_correction

    logger.info(f"Loading guides from {guides_df_path}")
    str_columns = ["guide_in_CDS"]
    guides_df = pd.read_csv(guides_df_path, dtype={col: str for col in str_columns})

    guides_df = correct_gene_names(guides_df, gene_name_corrections)
    guides_df = guides_df[guides_df["Hugo_Symbol"].isin(target_genes)]

    # Validate all target genes found
    found_genes = set(guides_df["Hugo_Symbol"].unique())
    missing_genes = target_genes - found_genes
    if missing_genes:
        logger.warning(f"Missing genes in data: {missing_genes}")

    gene2ensg = compute_gene2ensg(guides_df)

    gene2ensg_path.parent.mkdir(parents=True, exist_ok=True)
    with open(gene2ensg_path, "w") as f:
        json.dump(gene2ensg, f, indent=2)

    logger.info(f"Saved gene2ensg mapping ({len(gene2ensg)} genes) to {gene2ensg_path}")
    logger.info(f"Unique Gene IDs: {guides_df['Gene_ID'].nunique()}")

    gene2enst = compute_gene2enst(guides_df)

    gene2enst_path.parent.mkdir(parents=True, exist_ok=True)
    with open(gene2enst_path, "w") as f:
        json.dump(gene2enst, f, indent=2)

    logger.info(f"Saved gene2enst mapping ({len(gene2enst)} genes) to {gene2enst_path}")
    logger.info(f"Unique Transcript IDs: {guides_df['Transcript_ID'].nunique()}")


if __name__ == "__main__":
    main()
