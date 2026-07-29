from pathlib import Path
from typing import Any, Dict, List

import pandas as pd
from hydra import compose, initialize_config_dir
from omegaconf import DictConfig
from .config import CONFIG_DIR

N_HITS_COLNAME = "n_hits"
AF_CHAIN = "A"


def get_config() -> DictConfig:
    """Load hydra configuration."""
    with initialize_config_dir(config_dir=CONFIG_DIR, version_base=None):
        cfg = compose(config_name="config")
    return cfg


# Load configuration at module level
_cfg = None


def get_cached_config() -> DictConfig:
    global _cfg
    if _cfg is None:
        _cfg = get_config()
    return _cfg


def get_gene_to_uniprot() -> Dict[str, str]:
    """Get gene name to uniprot_id mapping from config."""
    cfg = get_cached_config()
    return dict(cfg.gene_name_to_uniprot_id)


def get_uniprot_to_gene() -> Dict[str, str]:
    """Get uniprot_id to gene name mapping from config."""
    gene_to_uniprot = get_gene_to_uniprot()
    return {v: k for k, v in gene_to_uniprot.items()}


def match_pdb_2_uniprot(pdb_code: str) -> pd.DataFrame:
    """Get residue mapping from PDB code to uniprot using UniprotMapper."""
    
    df = pd.DataFrame(columns=["code", "chain", "uniprot_id", "rno", "uniprot_rno"])
    
    # THIS NEEDS TO BE FILLED WITH CODE TO GET THE MAPPING FROM PDB TO UNIPROT
    # FROM UNIPROT RNO (RESIDUE NUMBER) TO PDB RNO
    
    return df


def get_processed_files_dirpath() -> Path:
    """Get the processed files directory path from config."""
    cfg = get_cached_config()
    return Path(cfg.processed_files_dirpath)


def load_data(pdb_code):
    """Load LFC data for genes associated with a PDB code or UniProt ID."""
    
    gene_to_uniprot = get_gene_to_uniprot()
    uniprot_to_gene = get_uniprot_to_gene()
    processed_files_dirpath = get_processed_files_dirpath()

    # Load the uniprot_id/gene/residues information
    if len(pdb_code) == 4:  # this is a real pdb code
        residue_match = match_pdb_2_uniprot(pdb_code)
        print(residue_match.columns)
        print(residue_match.head())
        residue_match["gene"] = residue_match["uniprot_id"].map(uniprot_to_gene)
        uniprot_genes = residue_match[["uniprot_id", "gene"]].drop_duplicates()

        uniprot_genes = [
            (row["uniprot_id"], row["gene"]) for row_i, row in uniprot_genes.iterrows()
        ]
    else:  # this is an AF: we have the uniprot_id
        uniprot_id = pdb_code
        gene = uniprot_to_gene[uniprot_id]
        uniprot_genes = [(uniprot_id, gene)]
        residue_match = pd.DataFrame(
            [(uniprot_id, gene, AF_CHAIN)], columns=["uniprot_id", "gene", "chain"]
        )

    all_data_dfs = []
    guide_data_dfs = []
    for uniprot_id, gene in uniprot_genes:
        if gene in gene_to_uniprot:
            # Use processed_files_dirpath/{gene}/{gene}_all_results.csv
            gene_dirpath = processed_files_dirpath / gene
            all_results_filepath = gene_dirpath / f"{gene}_all_results.csv"
            guide_results_filepath = gene_dirpath / f"{gene}_results_no_dup_guide.csv"
            
            data_df = pd.read_csv(all_results_filepath)
            data_df["uniprot_id"] = uniprot_id
            all_data_dfs.append(data_df)
            
            guide_df = pd.read_csv(guide_results_filepath)
            guide_df["uniprot_id"] = uniprot_id
            guide_data_dfs.append(guide_df)

    all_data = pd.concat(all_data_dfs)
    all_data["uniprot_rno"] = all_data["Protein_position"]

    guide_level_data = pd.concat(guide_data_dfs)
    guide_level_data["uniprot_rno"] = guide_level_data["Protein_position"]

    if len(pdb_code) == 4:  # this is a real pdb code
        all_data = all_data.merge(residue_match, on=["uniprot_rno", "uniprot_id"])
        guide_level_data = guide_level_data.merge(
            residue_match, on=["uniprot_rno", "uniprot_id"]
        )
    else:
        all_data = all_data.merge(residue_match, on=["uniprot_id"])
        all_data["rno"] = all_data["uniprot_rno"]
        guide_level_data = guide_level_data.merge(residue_match, on=["uniprot_id"])
        guide_level_data["rno"] = guide_level_data["uniprot_rno"]

    return all_data, guide_level_data


def preprocess_df(df_all, is_sig, remove_ot):
    if is_sig:
        df_all = df_all[df_all["logP"] >= 2]

    if remove_ot:
        df_all = df_all[df_all["is_ot"] != 1]

    df_all = df_all.dropna(subset=["uniprot_rno"])
    df_all["uniprot_rno"] = df_all["uniprot_rno"].astype(int)

    return df_all


def filter_by_main_conseq(df_all, guide_level_data: pd.DataFrame, mutation_type):
    mutation_guides = guide_level_data[guide_level_data["main_conseq"] == mutation_type]
    df_all = df_all[
        (df_all["sgRNA_ID"].isin(mutation_guides["sgRNA_ID"].values))
        & (df_all["Consequence"] == mutation_type)
    ]
    return df_all


def postprocess_data(pdb_code, results_d):
    if pdb_code == "6bcu":
        chains = ["A", "D", "W", "S"]
        results_d = {chain: rs for chain, rs in results_d.items() if chain in chains}
    return results_d


def add_n_hits(
    df_all, mutation_guides: pd.DataFrame, mut_guide_level_data, mutation_type
):
    mut_guide_level_data = mut_guide_level_data.drop_duplicates(
        subset=["sgRNA_ID", "cell_line", "editor"]
    )
    # n_hit_per_guide = mut_guide_level_data.groupby(['sgRNA_ID', 'editor'])['lt_cutoff'].sum()
    n_hit_per_guide = mut_guide_level_data.groupby(["sgRNA_ID", "editor"])[
        "lt_cutoff"
    ].mean()
    assert n_hit_per_guide.max() <= 1

    n_hit_df = n_hit_per_guide.reset_index().rename(
        {"lt_cutoff": N_HITS_COLNAME}, axis=1
    )

    df_all = df_all[
        (df_all["sgRNA_ID"].isin(mutation_guides["sgRNA_ID"].values))
        & (df_all["Consequence"] == mutation_type)
    ]
    df_all = df_all.merge(n_hit_df, on=["sgRNA_ID", "editor"])
    return df_all

def get_all_values(
    dataset_names: List[str],
    mutation_type: str,
    aggregation: str,
    pdb_code: str,
    is_sig: bool = False,
    remove_ot: bool = True,
    value_col: str = "LFC",
) -> Dict[str, Dict[str, Any]]:

    # Get data for corresponding genes in pdb file chains
    all_data, guide_level_data = load_data(pdb_code)

    # Match to selected cell_line - editor combos
    description_rows = []
    for description in dataset_names:
        d_split = description.split("_")
        cell_line = d_split[2]
        editor = d_split[4][:3]
        description_rows.append({"cell_line": cell_line, "editor": editor})
    description_df = pd.DataFrame(description_rows)
    df_all = all_data.merge(description_df, on=["cell_line", "editor"])
    guide_level_data = guide_level_data.merge(
        description_df, on=["cell_line", "editor"]
    )

    df_all = preprocess_df(df_all, is_sig, remove_ot)
    mutation_guides = guide_level_data[guide_level_data["main_conseq"] == mutation_type]
    df_all = df_all[
        (df_all["sgRNA_ID"].isin(mutation_guides["sgRNA_ID"].values))
        & (df_all["Consequence"] == mutation_type)
    ]

    chain_uniprot_df: pd.DataFrame = df_all[["chain", "uniprot_id"]].drop_duplicates()
    chain_uniprot_df = chain_uniprot_df.sort_values("chain")

    results_d = {}
    for uniprot_id in chain_uniprot_df["uniprot_id"].unique():
        subset_chain_df = chain_uniprot_df[chain_uniprot_df["uniprot_id"] == uniprot_id]
        chain = subset_chain_df["chain"].values[0]
        df_gene: pd.DataFrame = df_all[df_all["chain"] == chain]
        if aggregation == "minimum":
            values = df_gene.groupby("rno")[value_col].min().to_dict()
            # table2 = []
        elif aggregation == "median":
            values = df_gene.groupby("rno")[value_col].median().to_dict()
            # table2 = []
        elif aggregation == "average":
            values = df_gene.groupby("rno")[value_col].mean().to_dict()
            # table2 = []
        elif aggregation == "maximum":
            values = df_gene.groupby("rno")[value_col].max().to_dict()
            # table2 = []

        for chain in subset_chain_df["chain"].unique():
            results_d[chain] = {"values": values, "uniprot_id": uniprot_id}

    results_d = postprocess_data(pdb_code, results_d)

    print(values)

    return results_d


def get_lfcs(dataset_name, mutation_type, aggregation, pdb_code, is_sig, remove_ot):
    return get_all_values(
        [dataset_name],  # dataset_names is a list
        mutation_type,
        aggregation,
        pdb_code,
        is_sig,
        remove_ot,
        value_col="LFC",
    )


def get_zscores(dataset_names, mutation_type, aggregation, pdb_code, is_sig, remove_ot):
    return get_all_values(
        dataset_names,
        mutation_type,
        aggregation,
        pdb_code,
        is_sig,
        remove_ot,
        value_col="zscore",
    )


def get_n_hits(
    dataset_names, mutation_type, aggregation, pdb_code, remove_ot, is_sig: bool = False
):

    all_data, guide_level_data = load_data(pdb_code)

    description_rows = []
    for description in dataset_names:
        d_split = description.split("_")
        cell_line = d_split[2]
        editor = d_split[4][:3]
        description_rows.append({"cell_line": cell_line, "editor": editor})

    description_df = pd.DataFrame(description_rows)
    df_all = all_data.merge(description_df, on=["cell_line", "editor"])
    guide_level_data = guide_level_data.merge(
        description_df, on=["cell_line", "editor"]
    )

    df_all = preprocess_df(df_all, is_sig, remove_ot)
    print(df_all.shape)

    mutation_guides = guide_level_data[guide_level_data["main_conseq"] == mutation_type]
    df_all = df_all[
        (df_all["sgRNA_ID"].isin(mutation_guides["sgRNA_ID"].values))
        & (df_all["Consequence"] == mutation_type)
    ]

    # Keep guides that have the main selected mutation
    mut_guide_level_data = guide_level_data.merge(
        mutation_guides[["sgRNA_ID", "editor"]].drop_duplicates(), on=["sgRNA_ID", "editor"]
    )

    df_all = add_n_hits(df_all, mutation_guides, mut_guide_level_data, mutation_type)

    # We are only interested in n_hits per residue, we remove chain-uniprot_rno duplicates
    # This removes the cell_line and editor dimensions
    df_all = df_all.sort_values(
        N_HITS_COLNAME,
        ascending=False,
    ).drop_duplicates(["chain", "Gene", "uniprot_rno"])

    chain_uniprot_df: pd.DataFrame = df_all[["chain", "uniprot_id"]].drop_duplicates()
    chain_uniprot_df = chain_uniprot_df.sort_values("chain")

    # we want the maximum n_hits per residue
    aggregation = "maximum"

    results_d = {}
    for uniprot_id in chain_uniprot_df["uniprot_id"].unique():
        df_gene: pd.DataFrame = df_all[df_all["uniprot_id"] == uniprot_id]

        for chain, df_gene_chain in df_gene.groupby("chain"):
            if aggregation == "minimum":
                values = df_gene_chain.groupby("rno")[N_HITS_COLNAME].min().to_dict()
            elif aggregation == "median":
                values = df_gene_chain.groupby("rno")[N_HITS_COLNAME].median().to_dict()
            elif aggregation == "average":
                values = df_gene_chain.groupby("rno")[N_HITS_COLNAME].mean().to_dict()
            elif aggregation == "maximum":
                values = df_gene_chain.groupby("rno")[N_HITS_COLNAME].max().to_dict()

            if chain in results_d:  # adds the results to the chain
                results_d[chain]["values"] = {**results_d[chain]["values"], **values}
                results_d[chain]["uniprot_id"] = uniprot_id
                # there could be multiple uniprot_id for a single chain because
                # of chimeric constructions
            else:  # this should be the default in 99% case, 1% is chimeric construct
                results_d[chain] = {"values": values, "uniprot_id": uniprot_id}

    results_d = postprocess_data(pdb_code, results_d)

    print(results_d.keys())

    return results_d
