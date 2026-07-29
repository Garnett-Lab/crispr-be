#!/usr/bin/env python
"""Prepare variant files for local VEP (Variant Effect Predictor) analysis.

This script generates single and multiple edit variant files from guide data
for use with Ensembl's VEP tool.
"""

import itertools
from pathlib import Path

import hydra
import pandas as pd
from omegaconf import DictConfig

# Get absolute path to config directory
CONFIG_DIR = Path(__file__).parent.parent.parent / "config"


def create_multiple_combinations(list_positions: list) -> list:
    """Create all combinations of positions with more than one element.

    Args:
        list_positions: List of positions to combine.

    Returns:
        List of tuples containing all combinations with length > 1.
    """
    all_combinations = itertools.chain(
        *map(
            lambda x: itertools.combinations(list_positions, x),
            range(len(list_positions) + 1),
        )
    )
    return [i for i in all_combinations if len(i) > 1]


def get_real_edit_position_in_guide(
    guide: str, ini: int, direction: str
) -> list[tuple]:
    """Get the mapping between guide positions and genomic coordinates.

    Args:
        guide: The guide sequence.
        ini: Initial genomic position.
        direction: Direction of the guide ('left' or 'right').

    Returns:
        List of tuples mapping guide positions to genomic coordinates.
    """
    if direction == "left":
        match = list(
            zip(
                [[i + 1, k] for i, k in enumerate(guide)],
                list(range(ini + 22, ini + 2, -1)),
            )
        )
    elif direction == "right":
        match = list(
            zip(
                [[i + 1, k] for i, k in enumerate(guide)],
                list(range(ini, ini + 20)),
            )
        )
    return match


def map_combination_to_grna(
    list_positions: list,
    ini: int,
    guide: str,
    direction: str,
    wt_base: str,
    edited_base: str,
) -> list[str]:
    """Map edit combinations to edited guide sequences.

    Args:
        list_positions: List of edit positions.
        ini: Initial genomic position.
        guide: Original guide sequence.
        direction: Direction of the guide.
        wt_base: Wild-type base.
        edited_base: Edited base.

    Returns:
        List of edited guide sequences.
    """
    match = get_real_edit_position_in_guide(guide, ini, direction)
    combos = create_multiple_combinations(list_positions)
    all_combo_guides = []

    for combo in combos:
        current_combo = sum(
            list(
                map(
                    lambda elem: [
                        match[k][0][0] for k, x in enumerate(match) if elem in x
                    ],
                    combo,
                )
            ),
            [],
        )
        indices_combo = list(map(lambda x: x - 1, current_combo))

        edited_guide = []
        for h, nuc in enumerate(guide):
            for elem in indices_combo:
                if elem == h and nuc == wt_base:
                    nuc = edited_base
            edited_guide.append(nuc)

        all_combo_guides.append("".join(edited_guide))

    return all_combo_guides


def output_single_edits(
    input_filepath: str,
    wt_base: str,
    edited_base: str,
    output_path: str | Path,
    suffix: str,
) -> pd.DataFrame:
    """Generate single edit variant file for VEP.

    Args:
        input_filepath: Path to input CSV file with guide edits.
        wt_base: Wild-type base (A or C).
        edited_base: Edited base (G or T).
        output_path: Directory path for output files.
        suffix: Suffix for output filename (e.g., 'AKT1_ABE').

    Returns:
        DataFrame with single edit variants.
    """
    df = pd.read_csv(input_filepath)
    subset_df = df[["Location", "Edit_Location", "Direction"]].copy()
    subset_df.drop_duplicates(inplace=True)
    subset_df["chr"] = subset_df["Location"].apply(lambda x: x.split(":")[0])
    subset_df["sense"] = subset_df["Direction"].apply(
        lambda x: "+" if x == "right" else "-"
    )
    string_edition = f"{wt_base}/{edited_base}"
    subset_df["nuc"] = string_edition
    subset_df.sort_values(by=["Edit_Location"], ascending=True, inplace=True)

    final_df = subset_df[["chr", "Edit_Location", "Edit_Location", "nuc", "sense"]]
    output_filepath = Path(output_path) / f"single_edits_{suffix}.csv"
    print(f"Writing single edits to: {output_filepath}")

    final_df.drop_duplicates().to_csv(
        output_filepath, header=False, index=False, sep=" "
    )
    return final_df


def generate_multiple_edits(
    df: pd.DataFrame,
    wt_base: str,
    edited_base: str,
    output_path: str | Path,
) -> pd.DataFrame:
    """Generate multiple edit combinations from guide data.

    Args:
        df: DataFrame with guide edit data.
        wt_base: Wild-type base.
        edited_base: Edited base.
        output_path: Directory path for output files.

    Returns:
        DataFrame with all multiple edit combinations.
    """
    subset = df[
        [
            "gRNA_Target_Sequence",
            "CRISPR_PAM_Sequence",
            "Direction",
            "Location",
            "Edit_Location",
            "Strand",
        ]
    ]
    grouped = (
        subset.drop_duplicates().groupby("CRISPR_PAM_Sequence")["Edit_Location"].count()
    )

    all_multiples = pd.DataFrame()

    for seq in grouped[grouped > 1].index:
        seq_df = subset[subset["CRISPR_PAM_Sequence"] == seq].drop_duplicates()

        list_positions = seq_df["Edit_Location"].to_list()
        guide = seq_df["CRISPR_PAM_Sequence"].unique()[0]
        direction = seq_df["Direction"].unique()[0]
        ini = (
            seq_df["Location"]
            .apply(lambda x: x.split("-")[0].split(":")[1])
            .astype(int)
            .unique()[0]
        )

        combos = create_multiple_combinations(list_positions)
        edited_guides = map_combination_to_grna(
            list_positions, ini, guide, direction, wt_base, edited_base
        )

        edited_all = pd.DataFrame.from_dict(
            {
                "CRISPR_PAM_Sequence": [guide] * len(combos),
                "Combination": combos,
                "Edited_guide": edited_guides,
            }
        )
        all_multiples = pd.concat([all_multiples, edited_all], axis=0)

    # output_filepath = Path(output_path) / "multiple_edits_list_only.csv"
    # all_multiples.drop_duplicates().to_csv(output_filepath, index=False)
    return all_multiples


def _restrict_mutation(row: pd.Series) -> pd.Series:
    """Restrict mutation to only the changed nucleotides.

    Args:
        row: DataFrame row with mutation data.

    Returns:
        Series with restricted mutation coordinates and sequences.
    """
    original_sequence = row["CRISPR_PAM_Sequence"]
    modified_sequence = row["Edited_guide"]
    direction = row["Direction"]
    start_location = row["Edit_Location"]
    end_location = row["Edit_end"]
    locations = range(start_location, end_location + 1)

    if direction != "right":
        locations = reversed(locations)

    new_start_location = None
    new_end_location = None
    current_modif = ["", ""]
    buffer = ""
    z = zip(locations, original_sequence, modified_sequence)

    for location, original_n, modified_n in z:
        if original_n != modified_n:
            if new_start_location is None:
                new_start_location = location
                current_modif[0] = original_n
                current_modif[1] = modified_n
            else:
                current_modif[0] = current_modif[0] + buffer + original_n
                current_modif[1] = current_modif[1] + buffer + modified_n
                buffer = ""
            new_end_location = location
        else:
            if new_start_location is not None:
                buffer += modified_n

    assert new_end_location is not None

    if direction != "right":
        new_start_location, new_end_location = new_end_location, new_start_location

    return pd.Series(
        {
            "current_modif": current_modif,
            "new_start_location": new_start_location,
            "new_end_location": new_end_location,
        }
    )


def _reverse_guide(guide: str) -> str:
    """Reverse complement a guide sequence.

    Args:
        guide: Guide sequence to reverse complement.

    Returns:
        Reverse complemented sequence.
    """
    complement = {"A": "T", "T": "A", "C": "G", "G": "C"}
    return "".join([complement[n] for n in guide[::-1]])


def _reverse_replacement(replacement: str) -> str:
    """Reverse complement a replacement string (format: 'REF/ALT').

    Args:
        replacement: Replacement string in 'REF/ALT' format.

    Returns:
        Reverse complemented replacement string.
    """
    parts = replacement.split("/")
    return _reverse_guide(parts[0]) + "/" + _reverse_guide(parts[1])


def output_multiple_edits(
    input_filepath: str,
    wt_base: str,
    edited_base: str,
    output_path: str | Path,
    suffix: str,
) -> None:
    """Generate multiple edit variant file for VEP.

    Args:
        input_filepath: Path to input CSV file with guide edits.
        wt_base: Wild-type base (A or C).
        edited_base: Edited base (G or T).
        output_path: Directory path for output files.
        suffix: Suffix for output filename (e.g., 'AKT1_ABE').
    """
    input_df = pd.read_csv(input_filepath)
    multiple = generate_multiple_edits(input_df, wt_base, edited_base, output_path)

    df = pd.merge(input_df, multiple, on="CRISPR_PAM_Sequence", how="inner")

    subset_df = df[
        [
            "Location",
            "Combination",
            "Direction",
            "Edited_guide",
            "CRISPR_PAM_Sequence",
            "gRNA_Target_Sequence",
            "sgRNA_ID",
        ]
    ].copy()
    subset_df.drop_duplicates(inplace=True)

    subset_df["chr"] = subset_df["Location"].apply(lambda x: x.split(":")[0])
    subset_df["Edit_Location"] = subset_df["Location"].apply(
        lambda x: int(x.split(":")[1].split("-")[0])
    )
    subset_df["Edit_end"] = subset_df["Location"].apply(
        lambda x: int(x.split(":")[1].split("-")[1])
    )
    subset_df["sense"] = subset_df["Direction"].apply(
        lambda x: "+" if x == "right" else "-"
    )

    # Apply mutation restriction
    new_df = subset_df.apply(_restrict_mutation, axis=1)
    subset_df = pd.concat([subset_df, new_df], axis=1)

    subset_df["replacement"] = subset_df["current_modif"].apply(
        lambda x: x[0] + "/" + x[1]
    )
    subset_df["vep_replacement"] = subset_df[["replacement", "sense"]].apply(
        lambda s: (
            _reverse_replacement(s["replacement"])
            if s["sense"] == "-"
            else s["replacement"]
        ),
        axis=1,
    )
    subset_df["vep_sense"] = "+"

    # Output files
    output_path = Path(output_path)
    final_df = subset_df[
        [
            "chr",
            "new_start_location",
            "new_end_location",
            "vep_replacement",
            "vep_sense",
        ]
    ]
    match_df = subset_df[
        [
            "chr",
            "new_start_location",
            "new_end_location",
            "vep_replacement",
            "vep_sense",
            "CRISPR_PAM_Sequence",
            "sgRNA_ID",
            "replacement",
            "sense",
        ]
    ]

    # TODO: remove duplicates in final_df

    vep_filepath = output_path / f"multiple_edits_{suffix}.csv"
    match_filepath = output_path / f"multiple_edits_{suffix}_match.csv"

    match_df.to_csv(match_filepath, index=False)
    final_df.to_csv(vep_filepath, header=False, index=False, sep=" ")
    print(f"Writing multiple edits to: {vep_filepath}")


@hydra.main(version_base=None, config_path=str(CONFIG_DIR), config_name="config")
def main(cfg: DictConfig) -> None:
    """Main entry point for VEP variant file preparation."""
    target_genes = cfg.target_genes
    editors = cfg.editors
    editor_config = cfg.editor_config
    guide_with_edits_dirpath = Path(cfg.guide_with_edits_dirpath)
    vep_data_dirpath = Path(cfg.vep_data_dirpath)

    # Ensure output directory exists
    vep_data_dirpath.mkdir(parents=True, exist_ok=True)

    for gene in target_genes:
        for editor in editors:
            wt_base = editor_config[editor]["edit"]
            edited_base = editor_config[editor]["edit_to"]

            input_filepath = (
                guide_with_edits_dirpath / f"{gene}_{editor}_guides_w_edits.csv"
            )

            if not input_filepath.exists():
                print(f"Warning: {input_filepath} not found, skipping.")
                continue

            print(f"Processing {gene} with {editor}...")
            suffix = f"{gene}_{editor}"

            output_single_edits(
                input_filepath=str(input_filepath),
                wt_base=wt_base,
                edited_base=edited_base,
                output_path=vep_data_dirpath,
                suffix=suffix,
            )
            output_multiple_edits(
                input_filepath=str(input_filepath),
                wt_base=wt_base,
                edited_base=edited_base,
                output_path=vep_data_dirpath,
                suffix=suffix,
            )


if __name__ == "__main__":
    main()
