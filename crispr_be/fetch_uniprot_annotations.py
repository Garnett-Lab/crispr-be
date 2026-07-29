"""Fetch and process UniProt annotations for target genes."""

import json
import logging
from pathlib import Path

import hydra
import pandas as pd
import requests
from omegaconf import DictConfig

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)

# Get absolute path to config directory
CONFIG_DIR = Path(__file__).parent.parent / "config"


def make_request(upid: str) -> dict:
    """Fetch UniProt annotation data for a given UniProt ID.

    Args:
        upid: UniProt ID.

    Returns:
        JSON response from UniProt API.
    """
    url = f"http://www.uniprot.org/uniprot/{upid}.json"
    response = requests.get(url)
    json_response = json.loads(response.text)
    return json_response


def create_dataframe(json_data: dict, output_path: Path, gene: str) -> None:
    """Create and save annotation dataframe from UniProt JSON data.

    Args:
        json_data: JSON response from UniProt API.
        output_path: Path to save the CSV file.
        gene: Gene symbol.
    """
    dflist = []
    for elem in json_data["features"]:
        if elem["type"] == "Disulfide bond":
            entry = [
                elem["type"],
                [
                    int(elem["location"]["start"]["value"]),
                    int(elem["location"]["end"]["value"]),
                ],
                elem["description"],
            ]
        else:
            entry = [
                elem["type"],
                list(
                    range(
                        int(elem["location"]["start"]["value"]),
                        int(elem["location"]["end"]["value"]) + 1,
                    )
                ),
                elem["description"],
            ]

        if elem["type"] == "Binding site":
            entry.append(elem["ligand"]["name"])
        else:
            entry.append("")

        dflist.append(entry)

    df = pd.DataFrame(
        dflist,
        columns=["annotation_type", "uniprot_rno", "annotation_note", "ligand"],
    )
    dataframe_clean = df.explode("uniprot_rno").reset_index(drop=True)

    sequence = json_data["sequence"]["value"]
    rno2aa = {
        i + 1: aa for i, aa in enumerate(sequence)
    }  # residue counting starts at 1
    dataframe_clean["AA"] = dataframe_clean["uniprot_rno"].apply(lambda x: rno2aa[x])

    output_file = output_path / f"{gene}_annotations.csv"
    dataframe_clean.to_csv(output_file, index=False)
    logger.info(f"Saved annotations for {gene} to {output_file}")


def fetch_annotations(upid: str, output_path: Path, gene: str) -> None:
    """Fetch and process UniProt annotations for a gene.

    Args:
        upid: UniProt ID.
        output_path: Output directory path.
        gene: Gene symbol.
    """
    logger.info(f"Fetching annotations for {gene} ({upid})")
    req = make_request(upid)
    create_dataframe(req, output_path, gene)


@hydra.main(version_base=None, config_path=str(CONFIG_DIR), config_name="config")
def main(cfg: DictConfig) -> None:
    """Main entry point.

    Args:
        cfg: Hydra configuration.
    """
    # Get gene to UniProt ID mapping from config
    gene_name_to_uniprot_id = dict(cfg.gene_name_to_uniprot_id)

    # Get output directory from config
    output_dir = Path(cfg.uniprot_annotations_dirpath)
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info(f"Fetching annotations for {len(gene_name_to_uniprot_id)} genes")
    logger.info(f"Output directory: {output_dir}")

    # Fetch annotations for each gene
    for gene_name, uniprot_id in sorted(gene_name_to_uniprot_id.items()):
        try:
            fetch_annotations(uniprot_id, output_dir, gene_name)
        except Exception as e:
            logger.error(f"Error processing {gene_name} ({uniprot_id}): {e}")
            continue

    logger.info("Finished fetching all annotations")


if __name__ == "__main__":
    main()
