# Phylogenetic coherence in microbiome composition across environmental gradients

"Phylogenetic coherence in microbiome composition across environmental gradients"
Milena S. Chakraverti-Wuerthwein, Alissa Domenig, Seppe Kuehn,
bioRxiv 2026.06.07.730742; doi: https://doi.org/10.64898/2026.06.07.730742.

This repository contains code and analyses for studying how microbial communities respond to environmental gradients across soil and marine biomes. Using canonical correlation analysis (CCA) on large-scale microbiome surveys, we show that species' responses along the dominant environmental gradient are phylogenetically conserved, closely related taxa behave similarly, while secondary environmental axes show little to no such structure. These findings support a two-scale model of microbial community assembly, where deeply conserved evolutionary traits drive major compositional shifts, and more evolutionarily labile traits govern finer-scale variation.

---

## Overview

This repository contains the full analysis pipeline for the manuscript. The workflow proceeds in six stages:

0. **Environmental data preparation** (`code/read_data/`) — Loads and standardizes raw environmental metadata (soil Excel file; TARA ocean XLSX tables) into `envdata.csv`, used by CCA. Run once per dataset/config before Step 3.
1. **DADA2 denoising** (`code/DADA/`) — Amplicon sequence variant (ASV) inference from raw 16S reads, run per flowcell on a SLURM cluster. Soil dataset only; the ocean dataset uses already inferred miTAG sequences.
2. **Greengenes2 mapping** (`code/GG2/`) — ASV/miTAG sequences are aligned to the Greengenes2 backbone at a range of sequence identity thresholds using vsearch, then placed onto the GG2 reference phylogeny via QIIME2.
3. **Regularized CCA** (`code/CCA/`) — Canonical correlation analysis between microbial composition and environmental variables, run across taxonomic levels and cross-validation folds. Uses Elena Tuzhilina's [RCCA](https://github.com/ElenaTuzhilina/RCCA) package.
4. **Phylogenetic coherence** (`code/CCA/scripts/06_consenTRAIT.R`, `07_pagelsLambda.R`) — Adjusted consenTRAIT and Pagel's lambda applied to CCA abundance loadings.
5. **Reference distances** (`code/reference_distances/`) — Computes inter- and intra-group cophenetic distances at each taxonomic level using the CCA-filtered feature set. Establishes the phylogenetic distance scale that calibrates interpretation of consenTRAIT and Pagel's lambda results. Run on the cluster; output is **not committed** (too large for GitHub). Regenerate before figures if not present (~15 min on SLURM).
6. **Figures** (`code/manuscript_plotting/`) — Generates all manuscript and supplementary figures from committed result files.

Re-running the full pipeline from raw data is only needed if modifying the upstream steps.

---

## Setup

Two conda environments are required:

| Environment | Used for |
|---|---|
| `microbiome-phylo-coherence` | DADA2, vsearch/GG2 R steps, CCA, figures |
| `qiime2-amplicon-2025.7` | GG2 QIIME2 steps (01a, 03b, 04, 05) |

### 1. Project environment

```bash
conda env create -f environment.yml
conda activate <path-to-envs>/microbiome-phylo-coherence
Rscript setup_r_packages.R   # installs RCCA from GitHub and phloylm through CRAN; run once
```

### 2. QIIME2 environment

Follow the official installation instructions at https://docs.qiime2.org and install the `qiime2-amplicon-2025.7` distribution for your OS and install the GreenGenes2 plugin. The conda install command will look like:

```bash
# Create/activate a clean QIIME 2 env (example: 2025.7)
conda env create \
  --name qiime2-amplicon-2025.7 \
  --file https://raw.githubusercontent.com/qiime2/distributions/refs/heads/dev/2025.7/amplicon/released/qiime2-amplicon-ubuntu-latest-conda.yml

conda activate qiime2-amplicon-2025.7
qiime info

# Install the plugin without pulling extra build deps
pip install q2-greengenes2 --no-deps --no-build-isolation
# Add some supporting deps needed for this analysis
pip install redbiom msgpack nltk regex --no-build-isolation

# Refresh QIIME 2 plugin cache
qiime dev refresh-cache
qiime greengenes2 --help
```

Verify the exact URL and filename at the QIIME2 documentation for your platform.

---

## Data

### Greengenes2 backbone

The GG2 reference files are not included in this repository.

Download the following reference files from the official GG2 FTP and place them under `greengenes/` at the project root before running GG2 steps:

```bash
mkdir -p greengenes
wget -P greengenes/ \
  https://ftp.microbio.me/greengenes_release/current/2024.09.backbone.full-length.fna.qza \
  https://ftp.microbio.me/greengenes_release/current/2024.09.phylogeny.asv.nwk.qza \
  https://ftp.microbio.me/greengenes_release/current/2024.09.taxonomy.id.nwk.qza
```

### Raw sequencing and environmental data

**Soil dataset:** "A global atlas of the dominant bacteria found in soil".
Science 359, eaap9516 (2018). DOI:10.1126/science.aap9516

Raw 16S FASTQ files and environmental data are available from [figshare](https://figshare.com/s/82a2d3f5d38ace925492). 

16S rRNA amplicon sequencing data should be downloaded and then placed into `soil/data/downloaded_data/16S`. You can then run
```bash
python code/partition_by_flowcell.py \
    --input-dir soil/data/downloaded_data/16S \
    --output-dir soil/data/downloaded_data/16S \
    --flowcell-regex '[^-]+-(.+)'
```
to sort the sequencing files by flow cell, as reported in the sequence headers. Files whose headers are missing or unreadable are placed into an `UNKNOWN` subdirectory. The script prints a summary of flowcells and file counts, and writes a `flowcells.tsv` to the output directory.

Environmental data (`Dataset_figshare.xlsx`) should be downloaded and then placed into `soil/data/downloaded_data` and renamed: `env_metadata.xlsx`.

**Ocean dataset:** TARA Oceans Expedition

16S miTAG reference sequences (`16S.OTU.SILVA.reference.sequences.fna`) and OTU table (`miTAG.taxonomic.profiles.release.tsv`) are available from the TARA Oceans project via their [companion website](https://ocean-microbiome.embl.de/companion.html). Place these two files under `ocean/data/downloaded_data/16S`.

The environmental data can be found through two different sources:
- [Companion website](https://ocean-microbiome.embl.de/companion.html) contains associated metadata tables W1-W8 (`OM.CompanionTables.xlsx`) which includes sampling depth and nutrient levels.
- [Pangea database](https://doi.pangaea.de/10.1594/PANGAEA.875567) contains carbonate chemistry data for many of the samples (`TARA_SAMPLES_CONTEXT_ENV-DEPTH-CARB_20170515.xlsx`). Note: this data was not used in the final analysis.

Place these files under `ocean/data/downloaded_data/environmental`.

---

## Reproducing figures

All result files needed for figures are committed to this repository, **except the reference distance matrices** required for Fig. S13. Run Step 5 first if those outputs are not present. With the project environment active and working directory at the repo root:

```bash
source activate <path-to-envs>/microbiome-phylo-coherence
bash code/manuscript_plotting/run_all_plots.sh
```

For individual figures and exact commands, see [`code/manuscript_plotting/README.md`](code/manuscript_plotting/README.md).

---

## Running the full pipeline from scratch

Run all steps from the repo root with the appropriate environment active. DADA2 and GG2 steps are designed for SLURM; the CCA can be run locally or on a cluster.

### Step 0 — Environmental data preparation

Processes raw environmental metadata into `envdata.csv` required by CCA. Run once for each dataset before Step 3.

```bash
bash code/read_data/read_all.sh filtered   # or: full
```

See [`code/read_data/README.md`](code/read_data/README.md) for config structure and variable selection details.

### Step 1 — DADA2 denoising (soil only, SLURM)

```bash
cd code/DADA/sbatch && sbatch soil.sbatch
# After all array jobs complete:
Rscript code/DADA/09_mergeFlowcells.R --config code/DADA/config/soil.R
```

See [`code/DADA/README.md`](code/DADA/README.md) for config options, per-step details, and SLURM setup.

### Step 2 — Greengenes2 mapping (soil and ocean, SLURM)

```bash
bash code/GG2/sbatch/run_gg2_vsearch_sweep.sh --config code/GG2/config/soil.R
bash code/GG2/sbatch/run_gg2_vsearch_sweep.sh --config code/GG2/config/ocean.R
```

See [`code/GG2/README.md`](code/GG2/README.md) for config options, identity threshold sweep, and diagnostic outputs.

### Step 3 — CCA (SLURM, one job per taxonomic level)

```bash
bash code/CCA/sbatch/run_cca_tax_sweep.sh --config code/CCA/config/soil.R --perc-identity 0.90
bash code/CCA/sbatch/run_cca_tax_sweep.sh --config code/CCA/config/ocean.R --perc-identity 0.90
```

See [`code/CCA/README.md`](code/CCA/README.md) for the full step-by-step breakdown, config structure, and coherence scripts.

### Step 4 — Phylogenetic coherence (run locally after CCA)

See [`code/CCA/README.md`](code/CCA/README.md) for consenTRAIT (`06_consenTRAIT.R`) and Pagel's lambda (`07_pagelsLambda.R`) commands — they require specific `--direction` and `--mode` arguments detailed there.

### Step 5 — Reference distances (SLURM)

Computes inter- and intra-group cophenetic distances using the CCA-filtered feature set. Run once per dataset; output is not committed to the repository. Re-run whenever upstream data changes or after a fresh clone.

```bash
cd code/reference_distances
DATASET=soil  sbatch --job-name=ref_dist_soil  --output=./output_files/soil.out  --error=./output_files/soil.err  compute_reference_dists.sbatch
DATASET=ocean sbatch --job-name=ref_dist_ocean --output=./output_files/ocean.out --error=./output_files/ocean.err compute_reference_dists.sbatch
cd -
```

### Step 6 — Figures

```bash
bash code/manuscript_plotting/run_all_plots.sh
```

---

## Tests

Tests cover utility functions, CCA alignment functions, and consenTRAIT. No external data files are required.

```bash
REPO_ROOT=$(pwd) Rscript tests/run_tests.R
```

Expected: 0 failures across 23 tests.

---

## Repository structure

```
microbiome-phylo-coherence/
├── environment.yml          # lean conda env (R + vsearch)
├── setup_r_packages.R       # installs RCCA from GitHub; run once
├── code/
│   ├── setup.R / setup.sh   # repo-root and R_HOME resolution
│   ├── env_config.sh        # user-specific conda env path (edit before running sbatch)
│   ├── utility_functions.R
│   ├── partition_by_flowcell.py  # sorts soil FASTQs by flowcell header
│   ├── DADA/                # DADA2 denoising (soil only)
│   │   ├── sbatch/          # soil.sbatch (SLURM array job)
│   │   ├── config/soil.R
│   │   └── 01_*.R ... 09_*.R
│   ├── GG2/                 # Greengenes2 mapping (soil + ocean)
│   │   ├── sbatch/          # gg2_vsearch.sbatch + sweep wrapper
│   │   ├── scripts/
│   │   ├── functions/
│   │   └── config/
│   ├── CCA/                 # Regularized CCA + coherence
│   │   ├── sbatch/          # run_cca.sbatch + sweep wrapper
│   │   ├── scripts/
│   │   ├── functions/       # CCA_functions.R, consentrait_signed.R
│   │   └── config/
│   ├── reference_distances/ # inter/intra-group cophenetic distances (reference baseline)
│   ├── read_data/           # environmental data loading + diagnostics
│   └── manuscript_plotting/ # figure scripts (see README inside)
├── manuscript/              # committed figure PDFs produced by manuscript_plotting/
│   ├── soil/                # main soil figures
│   ├── ocean/               # main ocean figures
│   └── SI/                  # supplementary figures
├── soil/
│   ├── data/
│   │   ├── downloaded_data/ # gitignored; see Data section above
│   │   └── processed_data/
│   │       ├── 16S/         # gitignored; DADA2 + GG2 intermediate outputs
│   │       └── environmental/
│   └── results/
│       ├── CCA/             # CCA outputs (committed)
│       └── reference_distances/ # inter/intra-group cophenetic distances (not committed; run Step 5)
├── ocean/                   # same structure; no DADA2 outputs
├── greengenes/              # gitignored; download separately
└── tests/
    ├── run_tests.R
    └── testthat/
```
