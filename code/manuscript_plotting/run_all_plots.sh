#!/usr/bin/env bash
set -euo pipefail

PERC=0.90
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

run() { echo ">>> $*"; Rscript "$@"; }

# --- Main figures ---
run "$SCRIPT_DIR/plot_site_map.R" --soil --ocean

for DS in soil ocean; do
  if [ "$DS" = "soil" ]; then
    N_CDS=4
  else
    N_CDS=3
  fi
  run "$SCRIPT_DIR/plot_coherence.R"             --dataset "$DS" --n_cds "$N_CDS" --method consenTRAIT --perc_identity "$PERC"
  run "$SCRIPT_DIR/plot_crosstax_env_loadings.R" --dataset "$DS" --n_cds "$N_CDS" --perc_identity "$PERC"
  run "$SCRIPT_DIR/plot_corr_performance.R"      --dataset "$DS" --perc_identity "$PERC"
  run "$SCRIPT_DIR/SI_PCA_comparison.R"          --dataset "$DS" --perc_identity "$PERC"
done

for DIR in 1 3; do
  run "$SCRIPT_DIR/plot_abu_loadings.R" --dataset soil  --perc_identity "$PERC" --direction "$DIR"
  run "$SCRIPT_DIR/plot_abu_loadings.R" --dataset ocean --perc_identity "$PERC" --direction "$DIR"
done

# --- SI figures ---
run "$SCRIPT_DIR/SI_env_diagnostics.R" --figure S1
run "$SCRIPT_DIR/SI_env_diagnostics.R" --figure S2
run "$SCRIPT_DIR/SI_env_diagnostics.R" --figure S3
run "$SCRIPT_DIR/SI_GG2_performance.R"
run "$SCRIPT_DIR/SI_coherence.R"                    --method pagel --perc_identity "$PERC"
run "$SCRIPT_DIR/SI_crosstax_corr.R"                --perc_identity "$PERC"
run "$SCRIPT_DIR/SI_crosstax_stability.R"           --perc_identity "$PERC"
run "$SCRIPT_DIR/SI_example_alignment.R"
run "$SCRIPT_DIR/SI_example_hyperparameter_search.R" --dataset soil --perc_identity "$PERC"
echo "Running final SI figure (takes a while)..."
run "$SCRIPT_DIR/SI_interintra_dist.R"              --perc_identity "$PERC"
