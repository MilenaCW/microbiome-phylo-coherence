# Manuscript Figures

All commands run from the repository root. Activate the project environment first:

```bash
source activate <path-to-envs>/microbiome-phylo-coherence
```

## Main figures

### Figure 2a — Sampling site map
```bash
Rscript code/manuscript_plotting/plot_site_map.R --soil --ocean
```

### Figure 2c — Out-of-sample CCA correlation (soil and ocean)
```bash
Rscript code/manuscript_plotting/plot_corr_performance.R \
  --dataset soil  --perc_identity 0.90 --tax_level OTU
Rscript code/manuscript_plotting/plot_corr_performance.R \
  --dataset ocean --perc_identity 0.90 --tax_level OTU
```
Also writes `manuscript/<dataset>/<dataset>_corr_pvalues_<perc>_<tax>.csv`: per CD, the fraction of folds whose test correlation is <= the CD1 null mean. An asterisk marks CDs where no fold falls at or below the null (p < 1/n_folds).

### Figure 3b — Cross-taxonomic environmental loadings (soil)
```bash
Rscript code/manuscript_plotting/plot_crosstax_env_loadings.R \
  --dataset soil --perc_identity 0.90
```

### Figure 3c — Binarized abundance loadings on phylogenetic tree, CD1 and CD3 (soil)
```bash
Rscript code/manuscript_plotting/plot_abu_loadings.R \
  --dataset soil --perc_identity 0.99 --direction 1
Rscript code/manuscript_plotting/plot_abu_loadings.R \
  --dataset soil --perc_identity 0.99 --direction 3
```

### Figure 3d — Adjusted consenTRAIT (soil)
```bash
Rscript code/manuscript_plotting/plot_coherence.R \
  --dataset soil --method consentrait
```

### Figure 4a — Cross-taxonomic environmental loadings (ocean)
```bash
Rscript code/manuscript_plotting/plot_crosstax_env_loadings.R \
  --dataset ocean --perc_identity 0.90
```

### Figure 4b — Binarized abundance loadings, CD1 and CD3 (ocean)
```bash
Rscript code/manuscript_plotting/plot_abu_loadings.R \
  --dataset ocean --perc_identity 0.99 --direction 1
Rscript code/manuscript_plotting/plot_abu_loadings.R \
  --dataset ocean --perc_identity 0.99 --direction 3
```

### Figure 4c — Adjusted consenTRAIT (ocean)
```bash
Rscript code/manuscript_plotting/plot_coherence.R \
  --dataset ocean --method consentrait
```

## Supplementary figures

### Fig. S1 — Correlogram of filtered soil environmental variables
```bash
Rscript code/manuscript_plotting/SI_env_diagnostics.R --figure S1
```

### Fig. S2 — Missingness structure of full ocean environmental variable set
```bash
Rscript code/manuscript_plotting/SI_env_diagnostics.R --figure S2
```

### Fig. S3 — Correlogram of core ocean environmental variables
```bash
Rscript code/manuscript_plotting/SI_env_diagnostics.R --figure S3
```

### Fig. S4 — GreenGenes2 mapping performance vs. sequence identity threshold
```bash
Rscript code/manuscript_plotting/SI_GG2_performance.R --dataset soil
Rscript code/manuscript_plotting/SI_GG2_performance.R --dataset ocean
```

### Fig. S5 — Percent OTUs assigned per taxonomic level
```bash
Rscript code/manuscript_plotting/SI_taxonomic_assignment.R
```

### Fig. S6 — Cross-validation and hyperparameter selection (example)
```bash
Rscript code/manuscript_plotting/SI_example_hyperparameter_search.R \
  --dataset soil --tax_level OTU
```

### Fig. S7 — Sign-flip alignment of loading vectors across folds (example)
```bash
Rscript code/manuscript_plotting/SI_example_alignment.R \
  --dataset soil --tax_level OTU
```

### Fig. S9 — Out-of-sample CCA correlation across taxonomic levels
```bash
Rscript code/manuscript_plotting/SI_crosstax_corr.R
```

### Fig. S10 — CCA vs. PCA comparison (soil)
```bash
Rscript code/manuscript_plotting/SI_PCA_comparison.R \
  --dataset soil --n_cds 4
```

### Fig. S11 — CCA vs. PCA comparison (ocean)
```bash
Rscript code/manuscript_plotting/SI_PCA_comparison.R \
  --dataset ocean --n_cds 3
```

### Fig. S12 and S13 — Cross-taxonomic stability of canonical directions
```bash
Rscript code/manuscript_plotting/SI_crosstax_stability.R
```
*(S12 and S13 are both produced in a single run.)*

### Fig. S14 — Cophenetic distance distributions across taxonomic levels

> **Prerequisite:** requires Step 5 outputs (`{dataset}/results/reference_distances/`). Run Step 5 if not present.

```bash
Rscript code/manuscript_plotting/SI_interintra_dist.R
```

### Fig. S15 — Pagel's λ across canonical directions
```bash
Rscript code/manuscript_plotting/plot_coherence.R \
  --dataset soil  --method pagel
Rscript code/manuscript_plotting/plot_coherence.R \
  --dataset ocean --method pagel
```

## Running all figures at once

```bash
bash code/manuscript_plotting/run_all_plots.sh
```
