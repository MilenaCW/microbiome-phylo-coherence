# GreenGenes2 (GG2) Backbone Pipeline

## Overview
This folder contains a generalized GG2 backbone workflow that can be reused across
datasets (e.g., ocean miTAG, soil DADA2) via simple config files.

The pipeline:
- Starts from an existing sequence table + representative sequences.
- Maps features to the GG2 backbone with vsearch.
- Builds a GG2-based feature table and pruned tree.
- Produces cleaned `seqtab.csv`, `taxonomy.csv`, and `dna-sequences.fasta`
  suitable for downstream analyses (e.g., CCA).

See the sections below for script layout, configuration, environments, and outputs.

---

## 1. Setup

Edit these files before running the pipeline:

| File | What to edit |
|---|---|
| `code/setup.sh` / `code/setup.R` | Repo root path for your machine |
| `code/env_config.sh` | `CONDA_ENVS_PATH` — path to your conda environments directory |
| `code/GG2/sbatch/gg2_vsearch.sbatch` | SLURM `--account`, `--partition`, `--mail-user` |

Paths in config files (`input$table`, `output$directory`, etc.) are relative to the
project root set in `setup.sh` / `setup.R`.

---

## 2. Environments

See the [repo-level README](../../README.md) for full environment setup instructions
(both `microbiome-phylo-coherence` and `qiime2-amplicon-2025.7`).

This pipeline uses both environments:

| Environment | Used for |
|---|---|
| `microbiome-phylo-coherence` | R steps (01b, 03a, 06, 07, 08) and vsearch (02) |
| `qiime2-amplicon-2025.7` | QIIME2 steps (01a, 03b, 04, 05) |

Conda environment activation and deactivation happens inside `sbatch/gg2_vsearch.sbatch`.

---

## 3. GG2 Reference Files

Download the following reference files from the official GG2 FTP and place them under
`greengenes/` at the project root before running GG2 steps:

```bash
mkdir -p greengenes
wget -P greengenes/ \
  https://ftp.microbio.me/greengenes_release/current/2024.09.backbone.full-length.fna.qza \
  https://ftp.microbio.me/greengenes_release/current/2024.09.phylogeny.asv.nwk.qza \
  https://ftp.microbio.me/greengenes_release/current/2024.09.taxonomy.id.nwk.qza
```

---

## 4. Configuration

Configs live in `config/` and each evaluates to an R `list`. The config files assume
file paths originating from the repo base. The repo root is set in `code/setup.sh` and `code/setup.R` (edit those for your environment).

Shell scripts read config via `functions/access_config.R` (dot-path). From repo root:
`Rscript ./code/GG2/functions/access_config.R <config> <variable>`, e.g.
`Rscript ./code/GG2/functions/access_config.R config/ocean.R output.directory`.

Two examples:
- `config/ocean.R` (miTAG-style table):
  - `input` block sets format and where to find the sequence table + sequence file
    - `input$format = "mitag"`
    - `input$table = "./ocean/data/downloaded_data/16S/miTAG.taxonomic.profiles.release.tsv"` (miTAG taxonomic profiles TSV)
    - `input$sequences = "./ocean/data/downloaded_data/16S/16S.OTU.SILVA.reference.sequences.fna"` (SILVA reference sequences FASTA)
  - `output` block sets where the output files will be saved
    - `output$directory = "./ocean/data/processed_data/16S/GG2"`
  - `greengenes` block sets where to find the backbone, tree, taxonomy QZA files
    - `greengenes$backbone_qza = "./greengenes/2024.09.backbone.full-length.fna.qza"`
    - `greengenes$tree_qza = "./greengenes/2024.09.phylogeny.asv.nwk.qza"`
    - `greengenes$taxonomy_qza = "./greengenes/2024.09.taxonomy.id.nwk.qza"`
  - `vsearch$default_perc_identity = 0.99`
  - `sample_mapping` block specifies the sample label mapping (necessary for matching environment sampling to the sequence sampling)
    - `sample_mapping$lookup_file = "./ocean/data/processed_data/environmental/filtered/sample_id_lookup.csv"`
    - `sample_mapping$sample_label_column = "sample_label"` (sets the column to look up the study specified sample label)
    - `sample_mapping$analysis_id_column = "analysis_id"` (sets the column to look up the sample numbering used for this analysis)

- `config/soil.R` (DADA2 ASV table):
  - `input`
    - `input$format = "dada2"`
    - `input$table = "./soil/data/processed_data/16S/DADA/merged/merged_seqtab_asv.csv"` (merged ASV table; columns: `sample_id, ASV_1, ...`)
    - `input$sequences = "./soil/data/processed_data/16S/DADA/merged/merged_asv_sequences.fna"` (ASV sequences FASTA)
  - `output` block sets where the output files will be saved
    - `output$directory` = `./soil/data/processed_data/16S/GG2`
  - `greengenes` block sets where to find the backbone, tree, taxonomy QZA files
    - Same as ocean
  - `vsearch$default_perc_identity = 0.99`
  - `sample_mapping = NULL` block specifies the sample label mapping (not necessary for soil)

Each dataset uses the same pipeline; only the config changes.

---

## 5. Script Layout

All scripts live under `scripts/` and are numbered to reflect their order:

- `01a_prepare_backbone.sh` (shell, QIIME2 env via sbatch)
  - Exports GG2 backbone sequences from a QIIME2 artifact to
    `greengenes/backbone_full-length.fasta`.

- `01b_prepare_data.R` (R, microbiome-phylo-coherence env via sbatch)
  - Reads config and calls `functions/reformat_input.R` to harmonize input formats
    (miTAG vs DADA2) into a common `Feature ID × sample_N` table. Writes standardized
    inputs to `<output$directory>/input/`.

- `02_vsearch_mapping.sh` (shell, microbiome-phylo-coherence env via sbatch)
  - Runs vsearch closed-reference mapping against the GG2 backbone FASTA. Writes a
    UC file to `PI/intermediate/exports/`.

- `03a_build_gg2_seqtab.R` (R, microbiome-phylo-coherence env via sbatch)
  - Parses the UC file and aggregates counts from the input table to GG2 features.
    Writes the GG2 feature table, representative sequences, and ASV-to-GG2 mapping to
    `PI/intermediate/exports/`.

- `03b_import_to_qiime.sh` (shell, QIIME2 env via sbatch)
  - Converts the feature table to BIOM format and imports it as a QIIME2 artifact to
    `PI/intermediate/`.

- `04_tree_and_taxonomy.sh` (shell, QIIME2 env via sbatch)
  - Prunes the GG2 tree to mapped features (`phylogeny filter-tree`) and assigns
    taxonomy using the GG2 reference (`greengenes2 taxonomy-from-table`). Writes QZA
    artifacts to `PI/intermediate/`.

- `05_export_results.sh` (shell, QIIME2 env via sbatch)
  - Exports the pruned tree and taxonomy from QIIME2 artifacts to plain files in
    `PI/intermediate/exports/`.

- `06_clean_outputs.R` (R, microbiome-phylo-coherence env via sbatch)
  - Reads intermediate outputs and produces cleaned, CCA-ready files in `PI/final/`.

- `07_gg2_diagnostic.R` (R, microbiome-phylo-coherence env via sbatch)
  - Reads the cleaned outputs and produces occupancy, read distribution, and
    phylogenetic tree diagnostic plots in `PI/diagnostics/`.

- `08_gg2_summary_diagnostic.R` (R, microbiome-phylo-coherence env via sbatch)
  - Run after step 07 for multiple perc-identity values (e.g. after a sweep).
    Aggregates diagnostics across identity thresholds and writes summary plots and a
    combined TSV to `<output>/summary_diagnostics/`.

Shared helpers live in:
- `functions/reformat_input.R`
- `functions/diagnostic_functions.R`

---

## 6. SLURM Entry Points

### Main pipeline

`sbatch/gg2_vsearch.sbatch` orchestrates all steps. Examples:

```bash
cd code/GG2/sbatch

# Run full pipeline (01–07) for ocean config at default perc-identity
sbatch --output=path/to/output --error=/path/to/error gg2_vsearch.sbatch --config ./code/GG2/config/ocean.R

# Run mapping - diagnostics (steps 2-7) at specific identity (0.98)
sbatch --output=path/to/output --error=/path/to/error gg2_vsearch.sbatch \
  --config ./code/GG2/config/soil.R \
  --step vsearch_mapping --step build_seqtab --step tree_and_taxonomy \
  --step export_results --step clean_outputs --step run_diagnostics \
  --perc-identity 0.98
```

All conda environment activation/deactivation happens **inside this sbatch
script**, not in the R or shell helper scripts.

### Parameter sweep

`sbatch/run_gg2_vsearch_sweep.sh` submits a sweep across perc-identity values:

```bash
cd code/GG2/sbatch

bash run_gg2_vsearch_sweep.sh \
  --config ./code/GG2/config/soil.R \
  --start-id 0.80 --end-id 0.99 --step-size 0.01 \
  --steps all \
  --output-dir ./sweep_logs
```

After the sweep jobs complete, run step 08 once to aggregate diagnostics across
all perc-identity values (outputs in `<output>/summary_diagnostics/`).

```bash
# After a perc-identity sweep, run step 08 to aggregate diagnostics
sbatch --output=./sweep_logs/summary.out --error=./sweep_logs/summary.err gg2_vsearch.sbatch \
  --config ./code/GG2/config/ocean.R \
  --step summary_diagnostics
```

---

## 7. Output Structure

For each dataset (config) and identity threshold `PI` (e.g. `0.99`):

```text
<output$directory>/
  input/
    table_biom_ready.tsv          # BIOM-ready input table (Feature ID × sample_N)
    sequences_filtered.fna        # filtered input sequences
    sample_id_mapping.csv         # sample_N ↔ original_sample_id (+ sample_label for ocean)

  PI/
    intermediate/
      gg2_mapped_table.biom       # BIOM (FeatureTable[Frequency])
      gg2_mapped_table.qza        # QIIME2 artifact
      gg2_pruned_tree.qza         # pruned GG2 tree
      gg2_taxonomy.qza            # GG2 taxonomy for mapped features
      exports/
        otu_to_gg2.uc             # vsearch mappings (input feature → GG2 feature)
        feature-table.tsv         # GG2 feature table; first col 'Feature ID', cols sample_N
        dna-sequences.fasta       # representative GG2 sequences
        asv_to_gg2_map.tsv        # mapping from input Feature ID to GG2 Feature ID
        tree.nwk                  # pruned GG2 tree (Newick)
        taxonomy.tsv              # raw taxonomy; first col 'Feature ID'

    final/
      seqtab.csv                  # sample_id × feature table; feature columns are true IDs
      taxonomy.csv                # cleaned taxonomy; first col Feature_ID, matches seqtab columns
      dna-sequences.fasta         # copy of representative sequences
      tree.nwk                    # copy of pruned GG2 tree (Newick)

    diagnostics/
      gg2_occupancy_distribution.jpg
      gg2_occupancy_data.tsv
      read_distribution_counts.jpg
      read_distribution_fractions.jpg
      read_distribution_scatter.jpg
      gg2_tree_phyla.jpg

  summary_diagnostics/          # (after step 08; when multiple PI exist)
    gg2_mapping_success_vs_perc_identity.jpg
    gg2_occupancy_vs_perc_identity.jpg
    gg2_sweep_summary.tsv
```

---

## 8. Interpretation of Diagnostics

- **Occupancy distribution**:
  - How many original features (ASVs/OTUs) map to each GG2 feature.
  - High occupancy features represent strong clustering in GG2 space.

- **Read distribution (counts/fractions)**:
  - For each sample, fraction and count of reads retained after GG2 mapping.
  - Helps assess how stringent the chosen identity threshold is.

- **Read distribution scatter**:
  - Compares total reads per sample before vs. after GG2 mapping.
  - Includes an additional "GG2 Filtered" curve with low-abundance/noise
    features removed.

- **Tree colored by phylum**:
  - Circular plot showing GG2 tips colored by phylum.
  - Uses cleaned `taxonomy.csv` (top 20 phyla + "Other").
