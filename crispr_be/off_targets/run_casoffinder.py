import subprocess
from pathlib import Path

import hydra
from omegaconf import DictConfig
from tqdm import tqdm

CONFIG_DIR = Path(__file__).parent.parent.parent / "config"


def run_casoffinder(
    guides_path: Path,
    output_path: Path,
    casoffinder_path: str,
    device: str = "G0",
) -> None:
    """Run cas-offinder on a guides file."""
    if output_path.exists():
        print(f"Skipping {output_path.parent.name} (output already exists)")
        return

    if not guides_path.exists():
        print(f"Skipping {output_path.parent.name} (guides file not found)")
        return

    cmd = [casoffinder_path, str(guides_path), device, str(output_path)]
    print(f"Running cas-offinder for {output_path.parent.name}...")

    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        print(f"Error running cas-offinder for {output_path.parent.name}:")
        print(result.stderr)
    else:
        print(f"Completed {output_path.parent.name}")


@hydra.main(version_base=None, config_path=str(CONFIG_DIR), config_name="config")
def main(cfg: DictConfig) -> None:
    off_targets_dirpath = Path(cfg.off_targets_dirpath)
    casoffinder_path = cfg.casoffinder_path
    target_genes = cfg.target_genes

    for gene in tqdm(target_genes, desc="Running cas-offinder"):
        gene_dir = off_targets_dirpath / gene
        guides_path = gene_dir / "guides.csv"
        output_path = gene_dir / "alignment_guides.txt"

        run_casoffinder(guides_path, output_path, casoffinder_path)

    print("\nCas-OFFinder analysis complete")


if __name__ == "__main__":
    main()
