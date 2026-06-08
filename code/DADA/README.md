# DADA2 pipeline (code/DADA)

This folder contains a DADA2-based amplicon pipeline used to process paired-end 16S reads by flow cell. Each step is a standalone R script that reads a config file and writes results to a structured output directory. A SLURM script is provided for the soil dataset.

## High-level workflow
1. **Plot quality profiles** (`01_plotQualityProfile.R`)
2. **Filter + trim reads** (`02_filterAndTrim.R`)
3. **Learn error model** (`03_learnError.R`)
4. **Denoise / infer ASVs** (`04_sampleInference.R`)
5. **Merge paired reads** (`05_merge.R`)
6. **Build sequence table** (`06_buildTable.R`)
7. **Remove chimeras** (`07_removeChimera.R`)
8. **Diagnostics (pipeline)** (`08a_dada_diagnostics.R`)
9. **Diagnostics (seqtab/ASVs)** (`08b_seqtab_diagnostics.R`)
10. **Merge flow cells** (`09_mergeFlowcells.R`)

## Setup

Edit these files before running the pipeline:

| File | What to edit |
|---|---|
| `code/setup.sh` / `code/setup.R` | Repo root path for your machine |
| `code/env_config.sh` | `CONDA_ENVS_PATH` — path to your conda environments directory |
| `code/DADA/sbatch/soil.sbatch` | SLURM `--account`, `--partition`, `--mail-user` |

Paths in config files (`input$directory`, `output$directory`, etc.) are relative to the project root set in `setup.sh` / `setup.R`.

## Environment

This pipeline runs entirely in the `microbiome-phylo-coherence` conda environment.

See the [repo-level README](../../README.md) for environment setup instructions.

## Sequencing input format

Raw FASTQ files must be sorted into per-flowcell subdirectories before running the pipeline. For the soil dataset, use `code/partition_by_flowcell.py` to sort by the flowcell ID encoded in the sequence headers. For example, for the soil dataset:

```bash
python code/partition_by_flowcell.py \
    --input-dir soil/data/downloaded_data/16S \
    --output-dir soil/data/downloaded_data/16S \
    --flowcell-regex '[^-]+-(.+)'
```

The script prints a summary of flowcells and file counts, writes a `flowcells.tsv` to the output directory, and places files with missing or unreadable headers into an `UNKNOWN` subdirectory.

The pipeline expects the following directory structure:

```
<input_directory>/
  <FLOWCELL_1>/
    <sample>_R1.fastq.gz
    <sample>_R2.fastq.gz
  <FLOWCELL_2>/
    ...
```

The filename suffixes for forward/reverse reads are configurable via `patterns$forward` and `patterns$reverse` in the config.

## Output layout (per flow cell)
Each script creates a numbered subfolder under `<output_directory>/<FLOWCELL>/`:

```text
<output_directory>/<FLOWCELL>/
  01_qualityProfile/
    forward/*.jpeg                             # forward read quality profiles
    reverse/*.jpeg                             # reverse read quality profiles
  02_filterAndTrim/
    *_F_filt.fastq.gz                          # filtered forward reads
    *_R_filt.fastq.gz                          # filtered reverse reads
    <FLOWCELL>_filter_trim.csv                 # per-sample read counts before/after filtering
  03_learnError/
    <FLOWCELL>_errF.rds                        # forward error model
    <FLOWCELL>_errR.rds                        # reverse error model
    errF.jpeg                                  # forward error model plot
    errR.jpeg                                  # reverse error model plot
  04_sampleInference/
    <FLOWCELL>_dadaFs.rds                      # denoised forward reads
    <FLOWCELL>_dadaRs.rds                      # denoised reverse reads
  05_merge/
    <FLOWCELL>_mergers.rds                     # merged paired-end reads
  06_buildTable/
    <FLOWCELL>_seqtab.rds                      # sequence table (all lengths)
    <FLOWCELL>_seqtab_filtered.rds             # sequence table (length-filtered)
  07_removeChimera/
    <FLOWCELL>_seqtab_noChimeras.rds           # chimera-filtered sequence table
  08a_dada_diagnostics/
    <FLOWCELL>_track_reads.csv                 # per-sample read counts through pipeline
    <FLOWCELL>_track_reads.jpeg                # read tracking plot
    <FLOWCELL>_track_reads_fraction.jpeg       # read fraction tracking plot
    <FLOWCELL>_sequence_lengths.jpeg           # sequence length distribution
  08b_seqtab_diagnostics/
    <FLOWCELL>_num_samples_per_ASV.jpg         # number of samples per ASV
    <FLOWCELL>_n_ASVs_per_sample.jpg           # number of ASVs per sample
    <FLOWCELL>_frac_ASVs_per_sample.jpg        # fraction of ASVs per sample
    <FLOWCELL>_ASV_n_reads_per_sample.jpg      # ASV read counts per sample
    <FLOWCELL>_ASV_frac_reads_per_sample.jpg   # ASV read fractions per sample
    <FLOWCELL>_ASV_mean_relative_abundance.jpg # mean relative abundance per ASV
```

The final merge step writes a combined table under `<output_directory>/merged/`:

```text
<output_directory>/merged/
  merged_seqtab_asv.csv       # combined ASV table (all flowcells; rows: samples, cols: ASVs)
  merged_asv_sequences.fna    # FASTA of all ASV sequences
```

## Configuration
Config files live in `config/` and must evaluate to an R `list`. The provided example is `config/soil.R`.

Key fields:
- `input$directory`: base input path with per-flowcell subfolders
- `output$directory`: base output path
- `patterns$forward`, `patterns$reverse`: file suffixes (regex)
- `filtering`: `truncQ`, `trimLeft` (optional), `truncLen`, `maxN`, `maxEE`
- `sequence_table$length_range`: acceptable amplicon lengths
- `sequence_table$expected_length`: used in diagnostics plots
- `error_model$method`: `"standard"` or `"loess"`

## Running the pipeline
Each script can run a single flow cell or all flow cells in the input/output directory.

Single flow cell:
```
Rscript 01_plotQualityProfile.R --config code/DADA/config/soil.R --flowcell AA37J
Rscript 02_filterAndTrim.R     --config code/DADA/config/soil.R --flowcell AA37J
Rscript 03_learnError.R        --config code/DADA/config/soil.R --flowcell AA37J
Rscript 04_sampleInference.R   --config code/DADA/config/soil.R --flowcell AA37J
Rscript 05_merge.R             --config code/DADA/config/soil.R --flowcell AA37J
Rscript 06_buildTable.R        --config code/DADA/config/soil.R --flowcell AA37J
Rscript 07_removeChimera.R     --config code/DADA/config/soil.R --flowcell AA37J
Rscript 08a_dada_diagnostics.R --config code/DADA/config/soil.R --flowcell AA37J
Rscript 08b_seqtab_diagnostics.R --config code/DADA/config/soil.R --flowcell AA37J
```

All flow cells (omit `--flowcell`):
```
Rscript 02_filterAndTrim.R --config code/DADA/config/soil.R
```

Merge flow cells into a single ASV table (after all flow cells finish):
```
Rscript 09_mergeFlowcells.R --config code/DADA/config/soil.R
```

## SLURM batch job (soil)
`sbatch/soil.sbatch` runs steps 01–08 for one flow cell per array task. It loads the `microbiome-phylo-coherence` conda env and calls each R script with `--flowcell` set from the array index.
