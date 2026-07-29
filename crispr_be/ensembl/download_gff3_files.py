import gzip
import shutil
import urllib.request
from pathlib import Path

import hydra
from omegaconf import DictConfig
from tqdm import tqdm

CONFIG_DIR = Path(__file__).parent.parent.parent / "config"


def download_and_extract(url: str, output_path: Path) -> None:
    """Download a gzipped file and extract it."""
    extracted_path = output_path.with_suffix("")  # Remove .gz extension

    # Skip if already extracted
    if extracted_path.exists():
        print(f"Skipping {extracted_path.name} (already exists)")
        return

    # Download if .gz doesn't exist
    if not output_path.exists():
        print(f"Downloading {output_path.name}...")
        urllib.request.urlretrieve(url, output_path)
        print(f"Downloaded {output_path.name}")

    # Extract
    print(f"Extracting {output_path.name}...")
    with gzip.open(output_path, "rb") as f_in:
        with open(extracted_path, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)

    # Remove .gz file after extraction
    output_path.unlink()
    print(f"Extracted to {extracted_path.name}")


@hydra.main(version_base=None, config_path=str(CONFIG_DIR), config_name="config")
def main(cfg: DictConfig) -> None:
    base_url = cfg.chr_gff3_download_url
    chromosomes = cfg.chromosomes
    output_dir = Path(cfg.annotated_grch38_dirpath)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Download chromosome files
    for chrom in tqdm(chromosomes, desc="Downloading chromosomes"):
        filename = f"Homo_sapiens.GRCh38.115.chromosome.{chrom}.gff3.gz"
        url = f"{base_url}{filename}"
        output_path = output_dir / filename
        download_and_extract(url, output_path)

    print(f"\nAll files downloaded to {output_dir}")


if __name__ == "__main__":
    main()
