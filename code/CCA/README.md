# CCA Pipeline (code/CCA)

## Overview

This folder contains a regularized CCA (RCCA) pipeline for linking compositional (microbiome) and environmental data. It takes `envdata.csv` and GG2 outputs as inputs and produces cross-validated canonical loadings, a shuffle-based null distribution, and diagnostic plots.

**High-level workflow**
1. **Read data** (`00_read_data.R`) -- load and scale X (composition) and Y (environment)
2. **Hyperparameter search** (`01_hyperparam_search.R`) -- k-fold CV over lambda grid
3. **Loadings** (`02_loadings.R`) -- extract, normalize, and align across folds
4. **Null distribution** (`03_null_distribution.R`) -- shuffle-based null for CD1
5. **Model performance plots** (`04_model_performance_plots.R`) -- true vs null, env + tree loadings
6. **Cross-taxonomic compare** (`05_crosstaxonomic_compare.R`) -- compare results across tax levels (optional)
7. **consenTRAIT** (`06_consenTRAIT.R`) -- phylogenetic trait conservation (OTU only)
8. **Pagel's lambda** (`07_pagelsLambda.R`) -- phylogenetic signal in loadings (OTU only)
9. **PCA comparison** (`08_PCA_compare.R`) -- angle and EVR vs PCA (optional)

See `CCA_Notes.tex` for methodological background on RCCA and the cross-validation procedure.

---

## 1. Setup

Edit these files before running the pipeline:

| File | What to edit |
|---|---|
| `code/setup.sh` / `code/setup.R` | Repo root path for your machine |
| `code/env_config.sh` | `CONDA_ENVS_PATH` -- path to your conda environments directory |
| `code/CCA/sbatch/run_cca.sbatch` | SLURM `--account`, `--partition`, `--mail-user` |
| `code/CCA/sbatch/08_PCA_compare.sbatch` | SLURM `--account`, `--partition`, `--mail-user` |

Paths in config files (`data_path`, `results_path`, etc.) are relative to the project root set in `setup.sh` / `setup.R`.

---

## 2. Environment

This pipeline runs in the `microbiome-phylo-coherence` conda environment.

See the [repo-level README](../../README.md) for full environment setup instructions.

---

## 3. Input Data

Both input datasets must exist before running step 0.

### Environmental data

`envdata.csv` is produced by the `read_data` pipeline and must exist at `data_path/environmental/filtered/envdata.csv`. Override the location with `env_file` in the config if it lives elsewhere. Step 0 will error if the file is not found.

### Compositional data

GG2 pipeline outputs must exist under `data_path/16S/GG2/<perc_identity>/final/`:

- `seqtab.csv` -- sample x feature table (rows: samples, columns: GG2 feature IDs)
- `taxonomy.csv` -- taxonomy per feature (first column: `Feature_ID`)
- `tree.nwk` -- pruned GG2 phylogenetic tree

See `code/GG2/README.md` for instructions on running the GG2 pipeline.

---

## 4. Configuration

Configs live in `config/` and each evaluates to an R `list`. Paths are relative to the repo root.

| Field | Description |
|---|---|
| `dataset` | Label used in output file names |
| `data_path` | Root for input data (env + GG2 composition) |
| `results_path` | Root for CCA outputs |
| `env_file` | Override path to `envdata.csv` (default: `data_path/environmental/filtered/envdata.csv`) |
| `hyperparam$k` | Number of CV folds |
| `hyperparam$lambda1_range` | lambda1 grid (L1 penalty on composition loadings) |
| `hyperparam$lambda2_range` | lambda2 grid (`0` = no penalty on env loadings) |
| `hyperparam$seed` | RNG seed for fold assignment |
| `null$n_shuffles` | Number of shuffle iterations for the null distribution |
| `null$seeds` | Per-shuffle RNG seeds (length must equal `n_shuffles`) |
| `plotting_params$env_var_labels` | Named vector mapping column names to display labels for env loadings plots |
| `plotting_params$tax_palette` | Color palette for taxonomic levels in cross-taxonomic plots |

Two examples are provided:
- `config/soil.R`: DADA2 ASV input, soil variables (pH, C, N, P, C:N, clay-silt)
- `config/ocean.R`: miTAG input, ocean variables (depth, temperature, salinity, oxygen, nutrients)

---

## 5. Script Layout

All scripts live under `scripts/` and are numbered in order. All take `--config` (required), `--tax_level`, `--perc_identity`, and `--verbose`.

- `00_read_data.R` (R)
  - Loads `envdata.csv` and GG2 composition, coarse-grains by taxonomic level, scales, and saves X/Y matrices to `step0_data/`.

- `01_hyperparam_search.R` (R)
  - k-fold CV over the lambda1 (and optionally lambda2) grid; saves best hyperparameters and CV plot to `step1_hyperparam/`. Accepts `--n_cores` to parallelize.

- `02_loadings.R` (R)
  - Extracts loadings at best hyperparameters, normalizes to unit Euclidean norm, and aligns across folds to resolve sign ambiguity; saves loadings CSVs and alignment assessment plot to `step2_loadings/`. Accepts `--debug` to save intermediate fold-level files.

- `03_null_distribution.R` (R)
  - Per-seed row shuffle of X, followed by a full hyperparameter search and loadings extraction; saves null correlations and null loadings to `step3_null/`. Accepts `--n_cores` to parallelize.

- `04_model_performance_plots.R` (R)
  - Plots true vs null test correlation, env loadings scatter, and tree loadings (one file per canonical direction); writes to `step4_plots/`. Also accepts `--loadings_style continuous|binary`.

- `05_crosstaxonomic_compare.R` (R)
  - Discovers all tax levels run under a given `perc_identity` and compares canonical correlations and env loadings across levels; writes to `crosstax_compare/`. Takes `--config`, `--perc_identity`, `--verbose` (no `--tax_level`).

- `06_consenTRAIT.R` (R, OTU only)
  - Runs `consentrait_signed()` for each (canonical direction x fold) combination; saves aggregated tauD CSVs and diagnostic plots to `step6_consentrait/`. Takes `--config`, `--perc_identity`, `--verbose` (no `--tax_level`; always operates on OTU).

- `07_pagelsLambda.R` (R, OTU only)
  - Pagel's lambda via `phylolm` per fold for canonical directions 1..N. Takes `--n_cds N` to process all N directions in a single run; writes per-direction CSVs and then produces the aggregate `pagel_lambda.jpeg` automatically. Optional null via shuffling (`--n_shuffles`, `--seed`).

- `08_PCA_compare.R` (R)
  - For each available tax level, compares CCA canonical directions to PCA: computes the angle between each CCA canonical direction and PCA PC1 (with a shuffle null distribution), and cumulative explained variance ratios (EVR) for the CCA vs PCA subspaces on train and test splits. Processes all available tax levels in a single run. Takes `--dataset` (required; e.g. `soil`), `--perc_identity`, `--n_cds` (number of significant CDs; default 4), `--n_shuffles` (default 100), `--n_cores`, `--blas_threads`. Does **not** take `--config` or `--tax_level`.

Shared helpers:
- `functions/CCA_functions.R` -- alignment, CV, loadings extraction, `shuffle_sample`, CV plot helper, alignment assessment plot
- `functions/consentrait_signed.R` -- consenTRAIT implementation used by step 6

---

## 6. Running the Pipeline

### Core steps (0-4, per tax level)

Run from the repo root:

```bash
Rscript code/CCA/scripts/00_read_data.R             --config code/CCA/config/soil.R --tax_level OTU --perc_identity 0.90 --verbose
Rscript code/CCA/scripts/01_hyperparam_search.R      --config code/CCA/config/soil.R --tax_level OTU --perc_identity 0.90 --verbose
Rscript code/CCA/scripts/02_loadings.R               --config code/CCA/config/soil.R --tax_level OTU --perc_identity 0.90 --verbose
Rscript code/CCA/scripts/03_null_distribution.R      --config code/CCA/config/soil.R --tax_level OTU --perc_identity 0.90 --verbose
Rscript code/CCA/scripts/04_model_performance_plots.R --config code/CCA/config/soil.R --tax_level OTU --perc_identity 0.90 --verbose --loadings_style continuous
```

Step 3 can be slow (many shuffles x CV). Steps 1 and 3 accept `--n_cores` to parallelize across available cores.

### Cross-taxonomic compare (step 5)

After running steps 0-4 for each tax level to compare:

```bash
Rscript code/CCA/scripts/05_crosstaxonomic_compare.R --config code/CCA/config/soil.R --perc_identity 0.90 --verbose
```

### consenTRAIT (step 6, OTU only)

```bash
Rscript code/CCA/scripts/06_consenTRAIT.R --config code/CCA/config/soil.R --perc_identity 0.90 --verbose
```

### Pagel's lambda (step 7, OTU only)

```bash
Rscript code/CCA/scripts/07_pagelsLambda.R --config code/CCA/config/soil.R --perc_identity 0.90 --n_cds 4 --verbose
```

### PCA comparison (step 8)

After steps 0-4 have been run for at least one tax level:

```bash
Rscript code/CCA/scripts/08_PCA_compare.R \
  --dataset soil --perc_identity 0.90 --n_cds 4 --n_shuffles 50
```

Processes all available tax levels in one run. `--n_cds` should match the number of interpretable canonical directions (those above the null envelope in step 4). Pass `--n_cores` to parallelise the fold loop.

### SLURM batch jobs

`sbatch/run_cca.sbatch` runs steps 0-4 for a single tax level:

```bash
cd code/CCA/sbatch
sbatch run_cca.sbatch --config code/CCA/config/soil.R --tax-level OTU --perc-identity 0.90
```

`sbatch/run_cca_tax_sweep.sh` submits one job per tax level:

```bash
cd code/CCA/sbatch
bash run_cca_tax_sweep.sh --config code/CCA/config/soil.R --perc-identity 0.90 --levels all
# or a subset:
bash run_cca_tax_sweep.sh --config code/CCA/config/soil.R --perc-identity 0.90 --levels OTU,Genus,Phylum
```

`sbatch/08_PCA_compare.sbatch` runs step 8 for both datasets sequentially (soil then ocean), processing all tax levels in each:

```bash
cd code/CCA/sbatch
sbatch 08_PCA_compare.sbatch
# optional overrides:
sbatch 08_PCA_compare.sbatch --perc-identity 0.97 --n-cds 3
```

Uses 16 cores with BLAS threads fixed at 1, so all cores go to fold-level parallelism via `mclapply`.

---

## 7. Output Structure

Results live under `results_path/<perc_identity>/<tax_level>/`:

```text
<results_path>/
  <perc_identity>/
    <tax_level>/
      step0_data/
        X_matrix.csv                     # scaled composition matrix (samples x features)
        Y_matrix.csv                     # scaled env matrix (samples x variables)
        sample_id.csv                    # sample identifiers

      step1_hyperparam/
        cv_results.csv                   # test correlation by fold and lambda
        best_hyperparam.csv              # selected lambda1, lambda2
        cv_plot.jpg                      # CV correlation vs lambda grid

      step2_loadings/
        correlations_per_fold.csv        # per-fold canonical correlations
        env_loadings.csv                 # env variable loadings (long format)
        abundance_loadings.csv           # microbial abundance loadings (long format)
        alignment_combined.jpg           # raw vs aligned loadings by fold

      step3_null/
        null_correlations_per_fold.csv   # null canonical correlations (with seed)
        null_env_loadings.csv            # null env loadings (long format, with seed)

      step4_plots/
        test_correlation_vs_null.jpg     # true test correlation vs null distribution
        env_loadings.jpg                 # env variable loadings scatter
        abundance_loadings_tree_cd1.jpg  # microbial loadings on GG2 tree, CD1
        abundance_loadings_tree_cd2.jpg  # ... one file per canonical direction

      step6_consentrait/                 # (OTU only)
        tau_D.csv                        # observed tauD per canonical direction and fold
        null_tau_D.csv                   # null tauD per canonical direction, fold, and shuffle
        cd1/fold1/
          clades.csv                     # coherent clades for this (cd, fold)
          null_clades.csv                # null clades for this (cd, fold)
        ...                              # one cd{N}/fold{F}/ subdirectory per combination
        cluster_sizes.jpeg               # observed clade sizes vs canonical direction
        tauD.jpeg                        # observed tauD vs null distribution

      step7_pagelsLambda/                # (OTU only)
        stats/
          true_lambda_cd1.csv            # Pagel's lambda per fold, direction 1
          null_lambda_cd1.csv            # null lambda (if --n_shuffles > 0)
        pagel_lambda.jpeg                # aggregate plot across directions

      step8_PCA_compare/                 # (step 8; per tax level)
        angle_results.csv                # fold, cd, angle (radians vs PCA PC1)
        null_results.csv                 # fold, cd, null_mean, null_sd
        evr_results.csv                  # fold, j, method (CCA|PCA), split (train|test), evr

    crosstax_compare/                    # (step 5; one directory per perc_identity)
      correlation_comparative.jpg
      env_loadings_comparative.jpg
```

All outputs are CSV/TSV (no RDS) for portability. Use `readr::read_csv` (not base `read.csv`) when reading `X_matrix.csv`, `Y_matrix.csv`, or any file where column names are feature IDs -- base R mangles names containing hyphens by default.

---

## 8. Intermediate Diagnostics

Several plots are written during the pipeline as quality checks, before the final step-4 output plots.

**`cv_plot.jpg`** (step 1)
CV test correlation as a function of lambda. Check that the best-lambda region is stable (not at the boundary of `lambda1_range`) and that mean test correlation is clearly above zero. If the optimum sits at the grid edge, expand the range and rerun.

**`alignment_combined.jpg`** (step 2)
Raw loadings by fold (left panel) vs aligned loadings (right panel). All folds should show consistent sign and relative magnitude after alignment. Visible sign flips or large fold-to-fold variance indicate instability in the CCA solution at the chosen hyperparameters.

**`test_correlation_vs_null.jpg`** (step 4)
True test correlation (mean +/- SD across folds) per canonical direction, overlaid with the null distribution (mean +/- SE across shuffles) for CD1. Canonical directions where the true correlation does not exceed the null envelope are not interpretable.

**`env_loadings.jpg`** (step 4)
Environmental variable loadings per canonical direction and fold. High-loading variables are the dominant environmental drivers of that direction. Tight clustering of per-fold points indicates a robust result.

**`abundance_loadings_tree_cd*.jpg`** (step 4)
GG2 phylogenetic tree with tips colored by abundance loading magnitude (continuous) or sign (binary). Clustering of high-loading tips within clades suggests phylogenetic signal. One file per canonical direction.

**`correlation_comparative.jpg`** (step 5)
True test correlations across taxonomic levels for each canonical direction. Coarser levels (Phylum) typically show higher correlations due to aggregation; check whether signal is preserved at finer resolutions.

**`env_loadings_comparative.jpg`** (step 5)
Environmental loadings side-by-side across taxonomic levels. Consistent sign and rank across levels supports robustness; reversals suggest the signal is level-dependent.
