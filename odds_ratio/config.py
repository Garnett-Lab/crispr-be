"""
Configuration utilities for odds ratio analysis.

Uses hydra configuration from the crispr_be package.
"""

from pathlib import Path
from typing import Dict

from hydra import compose, initialize_config_dir
from omegaconf import DictConfig

# THIS NEEDS TO BE CHANGED TO THE CORRECT PATH FOR YOUR CONFIGURATION
CONFIG_DIR = "/my/config/path"

# Paths for structural data
PDB_DATABASE_PATH = Path("my/data/pdb/all")
CIF_DATABASE_PATH = Path("my/data/pdb/entries")
ALPHAFOLD_DATABASE_PATH = Path("my/data/pdb/alphafold2/")

# Module-level config cache
_cfg: DictConfig = None


def get_config() -> DictConfig:
    """Load hydra configuration from crispr_be config directory."""
    with initialize_config_dir(config_dir=CONFIG_DIR, version_base=None):
        cfg = compose(config_name="config")
    return cfg


def get_cached_config() -> DictConfig:
    """Get cached configuration, loading if not already cached."""
    global _cfg
    if _cfg is None:
        _cfg = get_config()
    return _cfg


def get_processed_files_dirpath() -> Path:
    """Get the processed files directory path from config."""
    cfg = get_cached_config()
    return Path(cfg.processed_files_dirpath)


def get_gene_to_uniprot() -> Dict[str, str]:
    """Get gene name to uniprot_id mapping from config."""
    cfg = get_cached_config()
    return dict(cfg.gene_name_to_uniprot_id)


def get_uniprot_to_gene() -> Dict[str, str]:
    """Get uniprot_id to gene name mapping from config."""
    gene_to_uniprot = get_gene_to_uniprot()
    return {v: k for k, v in gene_to_uniprot.items()}


def get_amino_acid_mapping() -> Dict[str, str]:
    """Get 3-letter to 1-letter amino acid mapping from config."""
    cfg = get_cached_config()
    return dict(cfg.amino_acid_3l_to_1l)


# Standard protein residue names (3-letter codes)
PROTEIN_RESIDUES = frozenset([
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY", "HIS", "HIE", "HID",
    "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
])

# Default analysis parameters
DEFAULT_RADIUS = 6.0
DEFAULT_PATCH_CENTER = "residue"
MIN_PATCH_RESIDUES = 10
MAX_PATCH_RESIDUES = 100
CORRECTED_PVALUE_THRESHOLD = 0.05
MAX_ODDS_RATIO = 1000.0
