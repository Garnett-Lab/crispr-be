# import os
from ast import literal_eval

# from concurrent.futures import ProcessPoolExecutor, as_completed
# from functools import partial
from pathlib import Path

import hydra
import pandas as pd
from omegaconf import DictConfig
from tqdm import tqdm

"""
Select guides with 4x 0MM or 1MM or less, now counting with the essential genes.
"""

CONFIG_DIR = Path(__file__).parent.parent.parent / "config"


def get_essential_genes(
    gene_essentials_path: Path, cell_lines: list[str]
) -> dict[str, set[str]]:
    """Returns a dict mapping cell_line -> set of essential genes (set for O(1) lookup)."""
    essential_genes = {}
    for cell_line in cell_lines:
        ess = pd.read_csv(gene_essentials_path / f"{cell_line}_essential.csv")
        ess.rename(columns={"Unnamed: 0": "Gene"}, inplace=True)
        essential_genes[cell_line] = set(ess["Gene"].values)  # Use set for fast lookup
    return essential_genes


def is_essential_vectorized(
    conseq_series: pd.Series, current_gene: str, gene_set: set
) -> pd.Series:
    """Vectorized version of is_essential check."""

    def check_single(strg):
        if not isinstance(strg, str):
            return False
        try:
            genes_in_conseq = literal_eval(strg)  # safer than eval()
            # Check if any essential gene is in conseq, but current_gene is not
            if current_gene in genes_in_conseq:
                return False
            return bool(gene_set & set(genes_in_conseq))  # set intersection
        except (ValueError, SyntaxError):
            return False

    return conseq_series.apply(check_single)


def select_guides_4(
    gene: str, cell_line: str, essential_genes_set: set, off_targets_dirpath: Path
):
    """Process a single gene-cell_line combination."""
    path = off_targets_dirpath / gene / "consequences.csv"

    # Read only needed columns if possible (reduces memory & I/O)
    conseq_df = pd.read_csv(path, low_memory=False)
    conseq_df = conseq_df.rename(columns={"id": "sgRNA_ID"})

    # Select sgrna with zero or one nucleotide mismatches
    ot_df = conseq_df[conseq_df["mismatches"] <= 1].copy()

    count_df = ot_df.groupby("sgRNA_ID", as_index=False)["mismatches"].count()
    count_df["is_ot_count"] = count_df["mismatches"] > 4

    ot_df = ot_df.drop("mismatches", axis=1)
    with_count_df = ot_df.merge(count_df, on="sgRNA_ID", how="inner")

    # Vectorized essential check
    with_count_df["has_essential"] = is_essential_vectorized(
        with_count_df["conseq"], gene, essential_genes_set
    )

    with_count_df["is_ot"] = (
        with_count_df["has_essential"] | with_count_df["is_ot_count"]
    )

    allg = with_count_df[
        ["sgRNA_ID", "mismatches", "is_ot_count", "has_essential", "is_ot"]
    ].copy()
    allg[["is_ot_count", "has_essential", "is_ot"]] = allg[
        ["is_ot_count", "has_essential", "is_ot"]
    ].astype(int)

    allg = allg.sort_values(["is_ot", "has_essential", "mismatches"], ascending=False)
    allg = allg.drop_duplicates("sgRNA_ID")

    allg.to_csv(off_targets_dirpath / gene / f"{gene}_{cell_line}.csv", index=False)
    return gene, cell_line


def process_gene_cell_line(args, essential_genes: dict, off_targets_dirpath: Path):
    """Wrapper for multiprocessing."""
    gene, cell_line = args
    select_guides_4(gene, cell_line, essential_genes[cell_line], off_targets_dirpath)
    return gene, cell_line


@hydra.main(version_base=None, config_path=str(CONFIG_DIR), config_name="config")
def main(cfg: DictConfig) -> None:

    target_genes = list(cfg.target_genes)
    off_targets_dirpath = Path(cfg.off_targets_dirpath)
    gene_essentials_path = Path(cfg.essential_genes_off_targets_dirpath)
    cell_lines = list(cfg.cell_lines)

    essential_genes = get_essential_genes(gene_essentials_path, cell_lines)

    # Single-threaded processing
    for gene in tqdm(target_genes, desc="Processing genes"):
        for cell_line in cell_lines:
            select_guides_4(
                gene, cell_line, essential_genes[cell_line], off_targets_dirpath
            )

    # # Parallel processing (commented out)
    # # Build list of all (gene, cell_line) combinations
    # tasks = [(gene, cell_line) for gene in target_genes for cell_line in cell_lines]
    #
    # # Parallel processing with progress bar
    # n_workers = min(os.cpu_count() or 4, len(tasks))
    #
    # with ProcessPoolExecutor(max_workers=n_workers) as executor:
    #     func = partial(process_gene_cell_line,
    #                    essential_genes=essential_genes,
    #                    off_targets_dirpath=off_targets_dirpath)
    #     futures = {executor.submit(func, task): task for task in tasks}
    #
    #     for future in tqdm(as_completed(futures), total=len(tasks), desc="Processing"):
    #         future.result()  # Raise any exceptions


if __name__ == "__main__":
    main()
