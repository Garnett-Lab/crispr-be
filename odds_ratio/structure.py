"""
Protein structure handling utilities.

Provides functions for loading and processing protein structures from PDB/CIF files.
"""

import gzip
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import List, Tuple

import numpy as np
from Bio.PDB import MMCIFParser, PDBParser, Selection
from scipy.spatial.distance import pdist, squareform

from .config import (
    ALPHAFOLD_DATABASE_PATH,
    CIF_DATABASE_PATH,
    PDB_DATABASE_PATH,
    PROTEIN_RESIDUES,
)


@dataclass
class StructureData:
    """Container for processed structure data."""
    atoms: List
    atom_to_residue_idx: List[int]
    residue_to_atom_idxs: List[List[int]]
    protein_residues: List
    protein_rnos: np.ndarray
    distance_matrix: np.ndarray


def get_cif_file_path(pdb_code: str) -> Path:
    """Get path to CIF file for a PDB code.
    
    Args:
        pdb_code: 4-letter PDB code
        
    Returns:
        Path to CIF file if it exists, None otherwise
    """
    entry_subdir = pdb_code[1:3]
    cif_path = CIF_DATABASE_PATH / entry_subdir / f"{pdb_code}.cif.gz"
    return cif_path if cif_path.exists() else None


def load_structure(pdb_code: str, chain_id: str):
    """Load a protein structure from PDB or CIF file.
    
    Args:
        pdb_code: PDB code (4-letter) or UniProt ID (for AlphaFold)
        chain_id: Chain identifier
        
    Returns:
        Bio.PDB Structure object
    """
    if len(pdb_code) == 4:
        # Standard PDB code
        pdb_path = PDB_DATABASE_PATH / f"{pdb_code}.pdb"
        
        if pdb_path.exists():
            parser = PDBParser(QUIET=True)
            structure = parser.get_structure(pdb_code, str(pdb_path))
        else:
            # Try CIF format
            cif_path = get_cif_file_path(pdb_code)
            if cif_path is None:
                raise FileNotFoundError(f"No structure file found for {pdb_code}")
            
            # Decompress CIF file to temp location
            temp_cif = Path("temp.cif")
            with gzip.open(cif_path, "rb") as f_in:
                with open(temp_cif, "wb") as f_out:
                    shutil.copyfileobj(f_in, f_out)
            
            parser = MMCIFParser(QUIET=True)
            structure = parser.get_structure(pdb_code, str(temp_cif))
    else:
        # AlphaFold model (UniProt ID)
        uniprot_id = pdb_code
        af_path = ALPHAFOLD_DATABASE_PATH / f"AF-{uniprot_id}-F1-model_v1.pdb"
        
        if not af_path.exists():
            raise FileNotFoundError(f"AlphaFold model not found for {uniprot_id}")
        
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure(uniprot_id, str(af_path))
    
    return structure


def extract_protein_residues(chain) -> List:
    """Extract protein residues from a chain.
    
    Args:
        chain: Bio.PDB Chain object
        
    Returns:
        List of protein residue objects
    """
    all_residues = Selection.unfold_entities(chain, "R")
    protein_residues = [
        res for res in all_residues 
        if res.get_resname() in PROTEIN_RESIDUES
    ]
    return protein_residues


def process_structure(pdb_code: str, chain_id: str) -> StructureData:
    """Load and process a protein structure for patch analysis.
    
    Args:
        pdb_code: PDB code or UniProt ID
        chain_id: Chain identifier
        
    Returns:
        StructureData with atoms, residue mappings, and distance matrix
    """
    structure = load_structure(pdb_code, chain_id)
    chain = structure[0][chain_id]
    protein_residues = extract_protein_residues(chain)
    
    # Build atom-to-residue and residue-to-atom mappings
    atom_to_residue_idx = []
    atoms = []
    
    for resi, residue in enumerate(protein_residues):
        for atom in Selection.unfold_entities(residue, "A"):
            atom_to_residue_idx.append(resi)
            atoms.append(atom)
    
    # Build reverse mapping (residue to atom indices)
    residue_to_atom_idxs = []
    current_atomi = 0
    for resi in range(len(protein_residues)):
        resi_atomidxs = []
        while current_atomi < len(atom_to_residue_idx) and atom_to_residue_idx[current_atomi] == resi:
            resi_atomidxs.append(current_atomi)
            current_atomi += 1
        residue_to_atom_idxs.append(resi_atomidxs)
    
    # Extract residue numbers
    protein_rnos = np.array([res.id[1] for res in protein_residues])
    
    # Compute all-atom distance matrix
    all_coords = [atom.coord for atom in atoms]
    distance_matrix = squareform(pdist(all_coords))
    
    return StructureData(
        atoms=atoms,
        atom_to_residue_idx=atom_to_residue_idx,
        residue_to_atom_idxs=residue_to_atom_idxs,
        protein_residues=protein_residues,
        protein_rnos=protein_rnos,
        distance_matrix=distance_matrix,
    )


def find_residues_in_radius(
    structure_data: StructureData,
    center_resi: int,
    radius: float,
) -> np.ndarray:
    """Find all residues within a radius of a center residue.
    
    Args:
        structure_data: Processed structure data
        center_resi: Index of center residue (0-based)
        radius: Search radius in Angstroms
        
    Returns:
        Array of residue numbers within the radius
    """
    atom_idxs = structure_data.residue_to_atom_idxs[center_resi]
    
    # Get minimum distance from center residue atoms to all other atoms
    distance_to_center = structure_data.distance_matrix[atom_idxs].min(axis=0)
    
    # Find atoms within radius
    atoms_in_radius = np.argwhere(distance_to_center < radius).flatten()
    
    # Map atoms to residue indices
    residue_idxs_in_radius = set(
        structure_data.atom_to_residue_idx[atomi] 
        for atomi in atoms_in_radius
    )
    
    # Get residue numbers
    pocket_rnos = structure_data.protein_rnos[list(residue_idxs_in_radius)]
    
    return pocket_rnos
