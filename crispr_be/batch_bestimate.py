"""Batch processing script for BEstimate analysis across multiple genes."""

import logging
import subprocess
import sys
from pathlib import Path
from typing import Literal

import hydra
from omegaconf import DictConfig

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)

# Get absolute path to config directory
CONFIG_DIR = Path(__file__).parent.parent / "config"

EditorType = Literal["ABE", "CBE"]


def run_bestimate(
    gene: str,
    editor: EditorType,
    editor_config: DictConfig,
    gene2ensg_path: Path,
    output_dir: Path,
    assembly: str,
    pam_seq: str,
    pam_window: str,
    activity_window: str,
    protospacer_length: int,
) -> bool:
    """
    Run BEstimate for a single gene and editor type.

    Args:
        gene: Hugo symbol of the gene
        editor: Base editor type ('ABE' or 'CBE')
        gene2ensg_path: Path to gene2ensg JSON mapping file
        output_dir: Output directory for results
        assembly: Genome assembly (default: GRCh38)
        pam_seq: PAM sequence pattern
        pam_window: PAM window range
        activity_window: Activity window range
        protospacer_length: Length of protospacer

    Returns:
        True if successful, False otherwise
    """
    editor_cfg = editor_config[editor]
    output_file = f"{gene}_{editor_cfg['suffix']}"

    cmd = [
        "python",
        "bestimate/beestimate_changed.py",
        "-gene",
        gene,
        "-assembly",
        assembly,
        "-pamseq",
        pam_seq,
        "-pamwin",
        pam_window,
        "-actwin",
        activity_window,
        "-protolen",
        str(protospacer_length),
        "-edit",
        editor_cfg["edit"],
        "-edit_to",
        editor_cfg["edit_to"],
        "-ofile",
        output_file,
        "-o",
        str(output_dir),
        "-gene2ensg",
        str(gene2ensg_path),
    ]

    logger.info(f"Processing {gene} with {editor}...")
    logger.debug(f"Command: {' '.join(cmd)}")

    try:
        _ = subprocess.run(
            cmd, check=True, capture_output=True, text=True, cwd=Path(__file__).parent
        )
        logger.info(f"Successfully processed {gene} with {editor}")
        return True
    except subprocess.CalledProcessError as e:
        logger.error(f"Failed to process {gene} with {editor}")
        logger.error(f"Error: {e.stderr}")
        return False
    except Exception as e:
        logger.error(f"Unexpected error processing {gene} with {editor}: {e}")
        return False


@hydra.main(version_base=None, config_path=str(CONFIG_DIR), config_name="config")
def main(cfg: DictConfig) -> int:
    """
    Main entry point for batch BEstimate processing.

    Args:
        cfg: Hydra configuration

    Returns:
        Exit code (0 for success, 1 for any failures)
    """
    # Get parameters from config
    genes = cfg.target_genes
    editors = cfg.editors
    gene2ensg_path = Path(cfg.gene2ensg_filepath)
    output_dir = Path(cfg.bestimate_data_dirpath)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Validate gene2ensg file exists
    if not gene2ensg_path.exists():
        logger.error(f"gene2ensg file not found: {gene2ensg_path}")
        logger.error("Please run compute_gene2ensg.py first")
        return 1

    logger.info(f"Starting batch processing for {len(genes)} genes with {editors}")
    logger.info(f"Output directory: {output_dir}")

    # Track results
    total = len(genes) * len(editors)
    successful = 0
    failed = []

    # Get editor configuration
    editor_config = cfg.editor_config

    # Get BEstimate parameters
    bestimate_params = cfg.bestimate_parameters

    for gene in genes:
        for editor in editors:
            success = run_bestimate(
                gene=gene,
                editor=editor,
                editor_config=editor_config,
                gene2ensg_path=gene2ensg_path,
                output_dir=output_dir,
                assembly=bestimate_params["assembly"],
                pam_seq=bestimate_params["pam_seq"],
                pam_window=bestimate_params["pam_window"],
                activity_window=bestimate_params["activity_window"],
                protospacer_length=bestimate_params["protospacer_length"],
            )

            if success:
                successful += 1
            else:
                failed.append((gene, editor))

    # Summary
    logger.info("\n" + "=" * 60)
    logger.info(f"Batch processing complete: {successful}/{total} successful")

    if failed:
        logger.warning(f"\n{len(failed)} job(s) failed:")
        for gene, editor in failed:
            logger.warning(f"  - {gene} ({editor})")
        return 1
    else:
        logger.info("All jobs completed successfully!")
        return 0


if __name__ == "__main__":
    sys.exit(main())
