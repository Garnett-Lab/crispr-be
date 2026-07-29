# Analysis pipeline for the CRISPR-BE screening of the PI3K pathway

Authors: Barbara Abreu and Benoit Baillif

This is the code to reproduce the analysis done in the manuscript
"Structural analysis of PI3K signalling cascade CRISPR base editing screens to inform cancer drug discovery"

**TODO: Add link to the manuscript**

## Required input

**TODO: Put the data on figshare, and add the link**

- BElib2_full_library.csv: list of screened sgRNAs
- Controls: list of controls sgRNAs
- Validated_controls.xlsx: list of validated controls
- MAGECK_processed: files with LFC
- data_on_off-targets: contains gene essentiality for each cell line computed with CoRe
- gene_effects_updated_all.csv: DEPMAP data

## Installation

### Conda environment

We have an environment.yml that lists all dependencies for reproducibility.  
This is installing the environment in the venv folder in your current directory
`mamba env create -f environment.yml -p ./venv`


### VEP Installation

This will be installed to a separate venv

```bash
# Go in the dir where you want on install the vep conda environment
cd /my/home/directory
mamba create -p ./venv_vep ensembl-vep -c bioconda
conda activate ./venv_vep

# Go in a dir where you want to install big data (26GB)
cd /my/data/directory
vep_install -a cf -s homo_sapiens -y GRCh38 -c ./vep_install/ --CONVERT
```

### Cas-OFFinder
Also requires [Cas-OFFinder](https://github.com/snugel/cas-offinder)  
You need the additional files (different from the ones for VEP)  
You need to be on a gpu server to run it (needs the OpenCL libs + uses 1 GPU)

## How to use this repo ?

The configuration required for all scripts is placed in the config.yaml
You need to setup the input paths in this file, and verify that all
intermediate files are up to date. This is also defining the constants
like the list of target genes, uniprot_ids, etc...

In the main scripts, the config.yaml is loaded with the hydra package.
Use the scripts in crispr-be in the following order:
- compute_gene_to_ensembl_ids.py
- batch_bestimate.py (uses a custom [Bestimate](https://github.com/Garnett-Lab/BEstimate) version in bestimate/bestimated.py)
- create_guide_lists.py
- vep/prepare_variant_files.py
- vep/run_vep.py
- ensembl/download_*
- off_targets/prepare_casoffinder_input_guides.py
- off_targets/run_casoffinder.py (on gpu server)
- off_targets/map_location.py
- fetch_uniprot_annotations.py

The current repo does not include the code to retrieve pdb ids from
uniprot_ids, you might need this when you need to link Uniprot residues numbers to
PDB residue indices 

Then the main analysis are in the Jupyter notebooks in clean_workflow:
- control_and_merge_data.ipynb
- data_analysis.ipynb

The odds_ratio patch code is isolated. It was embedded in another repo,
this is a shareable version. You need to fill some code in 
config.py and fetch_lfc_values.py
In this code:
Dataset_name needs to be in the format `"EXP_0_HGC27_0_ABE8e_0"` as long as it matches the split
mutation_type can be synonymous, missense, low_imp_splice, high_imp_splice, stop
aggregation can be minimum, median, average, maximum