from pathlib import Path

import hydra
import pandas as pd
from omegaconf import DictConfig
from tqdm import tqdm

"""
Get gene consequences of targeted positions using ensembl annotations
Changed to add the most relevant off targets
"""

CONFIG_DIR = Path(__file__).parent.parent.parent / "config"


def process_alignment(alignment_path: str):
    columns_names = ["guide", "desc", "ini", "target", "direction", "mismatches", "id"]
    alg = pd.read_csv(alignment_path, sep="\t", header=None, names=columns_names)
    alg["chr"] = alg["desc"].apply(
        lambda x: x.split(":")[3] if "chromosome" in x else ""
    )
    return alg


@hydra.main(version_base=None, config_path=str(CONFIG_DIR), config_name="config")
def main(cfg: DictConfig) -> None:

    off_targets_dirpath = Path(cfg.off_targets_dirpath)
    chr_gff3_dirpath = Path(cfg.annotated_grch38_dirpath)

    target_genes = cfg.target_genes
    chrs = cfg.chromosomes

    # Keeping only these types, i.e. removing lnc_RNA, transcript with unknown CDS
    # CHR_TYPES = ['exon', 'three_prime_UTR', 'five_prime_UTR']
    CHR_TYPES = ["mRNA"]

    CHR_COLUMN_NAMES = [
        "seqid",
        "source",
        "type",
        "start",
        "end",
        "score",
        "strand",
        "phase",
        "attributes",
    ]

    print("Read gffs")
    chr_gffs: dict[str, pd.DataFrame] = {}
    for chr in tqdm(chrs):
        chr_file = chr_gff3_dirpath / f"Homo_sapiens.GRCh38.115.chromosome.{chr}.gff3"
        gff = pd.read_csv(
            chr_file, sep="\t", skiprows=9, header=None, names=CHR_COLUMN_NAMES
        )
        gff.dropna(axis=0, inplace=True)
        gff = gff[gff["type"].isin(CHR_TYPES)]
        gff["id"] = gff["attributes"].apply(lambda x: x.split(";")[0].split(":")[1])
        gff["info2"] = gff["attributes"].apply(
            lambda x: (
                x.split("Name")[1].split(";")[0].strip("=") if "Name" in x else ""
            )
        )
        gff.drop_duplicates(inplace=True)
        gff["start"] = gff["start"].astype(int)
        gff["end"] = gff["end"].astype(int)

        chr_gffs[chr] = gff

    for gene in tqdm(target_genes[22:]):

        alignment_path = off_targets_dirpath / gene / "alignment_guides.txt"
        alignment = process_alignment(alignment_path)

        alignment["targeted_dna"] = alignment["target"].str.upper()
        alignment["ini"] = alignment["ini"].astype(int)
        chr_alignments = []
        for chr, chr_alignment in alignment.groupby("chr"):
            if len(chr) > 0:
                gff = chr_gffs[chr]

                # Create IntervalIndex for fast range lookups
                intervals = pd.IntervalIndex.from_arrays(
                    gff["start"], gff["end"], closed="neither"
                )
                gff_indexed = gff.set_index(intervals)

                def get_conseq(position):
                    try:
                        matches = gff_indexed.loc[position]
                        if isinstance(matches, pd.Series):
                            return [matches["info2"]]
                        return matches["info2"].unique().tolist()
                    except KeyError:
                        return None

                chr_alignment["conseq"] = chr_alignment["ini"].apply(get_conseq)
            else:
                chr_alignment["conseq"] = None

            chr_alignments.append(chr_alignment)

        alignment = pd.concat(chr_alignments)
        alignment_path = off_targets_dirpath / gene / "consequences.csv"
        alignment.to_csv(alignment_path, index=False)


if __name__ == "__main__":
    main()
