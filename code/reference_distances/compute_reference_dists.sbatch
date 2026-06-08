#!/bin/bash
#SBATCH --job-name=reference_dists
#SBATCH --output=./output_files/reference_dists.out
#SBATCH --error=./output_files/reference_dists.err
#SBATCH --account=<account-name>
#SBATCH --partition=<partition-name>
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=<your-email@domain>

# Reference cophenetic distance computation.
# Set DATASET before submitting (soil or ocean).
# sbatch options on the command line override the #SBATCH directives above.
#
# Usage (from code/reference_distances/):
#   DATASET=soil  sbatch --job-name=ref_dist_soil  --output=./output_files/soil.out  --error=./output_files/soil.err  compute_reference_dists.sbatch
#   DATASET=ocean sbatch --job-name=ref_dist_ocean --output=./output_files/ocean.out --error=./output_files/ocean.err compute_reference_dists.sbatch
#
# Optional env var overrides (set before sbatch, same as DATASET):
#   N_SAMPLE=100000   -- pairs to sample per level per type (-1 = all; default)
#   PERC_IDENTITY=0.90
#
# Output: {dataset}/results/reference_distances/{perc_identity}/dist_{level}.rds

set -euo pipefail

if [[ -z "${DATASET:-}" ]]; then
  echo "ERROR: DATASET must be set to 'soil' or 'ocean' before submitting." >&2
  echo "  Example: DATASET=soil sbatch compute_reference_dists.sbatch" >&2
  exit 1
fi

if [[ "$DATASET" != "soil" && "$DATASET" != "ocean" ]]; then
  echo "ERROR: DATASET must be 'soil' or 'ocean'. Got: ${DATASET}" >&2
  exit 1
fi

# Resolve setup.sh relative to the submit directory.
# Convention: submit from code/reference_distances/ so ../setup.sh -> code/setup.sh.
if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
  SETUP_SH="${SLURM_SUBMIT_DIR}/../setup.sh"
else
  SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
  SETUP_SH="${SCRIPT_DIR}/../setup.sh"
fi
source "$SETUP_SH"
source "$(dirname "$SETUP_SH")/env_config.sh"

N_SAMPLE="${N_SAMPLE:--1}"
PERC_IDENTITY="${PERC_IDENTITY:-0.90}"

echo "Reference distance computation"
echo "    dataset        = ${DATASET}"
echo "    perc_identity  = ${PERC_IDENTITY}"
echo "    n_sample       = ${N_SAMPLE}"

mkdir -p "code/reference_distances/output_files"

module load python # adjust as needed for your cluster
source activate "${PROJECT_ENV}"

Rscript code/reference_distances/compute_reference_dists.R \
  --dataset "${DATASET}" \
  --perc_identity "${PERC_IDENTITY}" \
  --n_sample "${N_SAMPLE}" \
  --verbose

conda deactivate

echo "Done."
