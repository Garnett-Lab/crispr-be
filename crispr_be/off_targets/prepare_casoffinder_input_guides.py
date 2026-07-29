from pathlib import Path

import hydra
import pandas as pd
from omegaconf import DictConfig

CONFIG_DIR = Path(__file__).parent.parent.parent / "config"


def generate_guide_file_for_gene(
    gene_guides_df: pd.DataFrame, output_path: str, grch38_dirpath: str
):
    gene_guides_df["c"] = gene_guides_df["CRISPR_PAM_Sequence"].str[:-3] + "NNN"
    gene_guides_df["b"] = "4"

    t = gene_guides_df[["c", "b", "gRNA_ID"]].drop_duplicates()

    with open(output_path, "w") as f:
        f.write(f"{grch38_dirpath}\n")
        f.write("NNNNNNNNNNNNNNNNNNNNNGN\n")
        t.to_csv(f, index=False, sep=" ", header=False)


@hydra.main(version_base=None, config_path=str(CONFIG_DIR), config_name="config")
def main(cfg: DictConfig) -> None:
    str_columns = ["guide_in_CDS"]
    guides_df = pd.read_csv(
        cfg.guide_data_filepath, dtype={col: str for col in str_columns}
    )
    off_targets_dirpath = cfg.off_targets_dirpath
    grch38_dirpath = cfg.grch38_dirpath

    gene_name_correction = cfg.gene_name_correction
    for wrong_gene, correct_gene in gene_name_correction.items():
        guides_df["Hugo_Symbol"] = guides_df["Hugo_Symbol"].str.replace(
            wrong_gene, correct_gene, regex=False
        )
        guides_df["gRNA_ID"] = guides_df["gRNA_ID"].str.replace(
            wrong_gene, correct_gene, regex=False
        )

    for gene in cfg.target_genes:
        gene = gene_name_correction.get(gene, gene)
        gene_guides_df = guides_df[guides_df["Hugo_Symbol"] == gene]
        output_path = Path(off_targets_dirpath) / gene / "guides.csv"
        output_path.parent.mkdir(parents=True, exist_ok=True)
        generate_guide_file_for_gene(gene_guides_df, output_path, grch38_dirpath)


if __name__ == "__main__":
    main()
