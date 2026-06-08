#!/bin/bash
#SBATCH --job-name=pca_compare
#SBATCH --output=./output_files/pca_compare.out
#SBATCH --error=./output_files/pca_compare.err
#SBATCH --account=<account-name>
#SBATCH --partition=<partition-name>
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=4G
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=<your-email@domain>

# Step 8: PCA comparison. Runs soil then ocean sequentially; folds parallelised within each dataset.
#
# Usage:
#   sbatch 08_PCA_compare.sbatch [OPTIONS]
#
# Options:
#   --perc-identity ID   Percent identity threshold (default: 0.90)
#   --n-cds N            Number of significant canonical directions (default: 4)
#   --n-shuffles N       Shuffles per fold for angle null distribution (default: 100)
#   --help, -h           Show this help message
#
# Examples:
#   sbatch 08_PCA_compare.sbatch
#   sbatch 08_PCA_compare.sbatch --perc-identity 0.97 --n-cds 3

set -euo pipefail

usage() {
  sed -n '15,28p' "$0" | sed 's/^#//'
  exit 1
}

# Source shared setup (strict mode, R_HOME, repo root).
# Under sbatch the script runs from a spool copy so $0 is wrong; use SLURM_SUBMIT_DIR
# (directory from which sbatch was invoked). When not under SLURM, use script location.
# Convention: submit from code/CCA/sbatch so ../../setup.sh resolves to code/setup.sh.
if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
  SETUP_SH="${SLURM_SUBMIT_DIR}/../../setup.sh"
else
  SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
  SETUP_SH="${SCRIPT_DIR}/../../setup.sh"
fi
source "$SETUP_SH"
source "$(dirname "$SETUP_SH")/env_config.sh"

# Default values
PERC_IDENTITY="0.90"
N_CDS=""  # default varies by dataset: soil=4, ocean=3; override with --n-cds
N_SHUFFLES=100

while [[ $# -gt 0 ]]; do
  case "$1" in
    --perc-identity)
      PERC_IDENTITY="$2"
      shift 2
      ;;
    --n-cds)
      N_CDS="$2"
      shift 2
      ;;
    --n-shuffles)
      N_SHUFFLES="$2"
      shift 2
      ;;
    --help|-h)
      usage
      ;;
    *)
      echo "Unknown option: $1" >&2
      usage
      ;;
  esac
done

# BLAS threads fixed at 1; all cores go to fold-level parallelism via mclapply
BLAS_T=1
if [[ -n "${SLURM_CPUS_PER_TASK:-}" ]]; then
  N_CORES=$(( SLURM_CPUS_PER_TASK / BLAS_T ))
  [[ $N_CORES -lt 1 ]] && N_CORES=1
else
  N_CORES=1
fi

export OMP_NUM_THREADS="${BLAS_T}"
export MKL_NUM_THREADS="${BLAS_T}"
export OPENBLAS_NUM_THREADS="${BLAS_T}"
export VECLIB_MAXIMUM_THREADS="${BLAS_T}"
export NUMEXPR_NUM_THREADS="${BLAS_T}"

echo "08_PCA_compare running with parameters:"
echo "    perc_identity=${PERC_IDENTITY}"
echo "    n_cds=${N_CDS:-<dataset default: soil=4, ocean=3>}"
echo "    n_shuffles=${N_SHUFFLES}"
echo "    n_cores=${N_CORES}  (blas_threads=${BLAS_T})"
echo "    SLURM_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK:-}"
echo
env | grep -E 'OMP_NUM_THREADS|MKL_NUM_THREADS|OPENBLAS_NUM_THREADS|VECLIB_MAXIMUM_THREADS|NUMEXPR_NUM_THREADS' || true
echo

module load python # adjust as needed for your specific cluster
source activate "${PROJECT_ENV}"

DATASETS=(soil ocean)

for DATASET in "${DATASETS[@]}"; do
  # Resolve n_cds: CLI override takes precedence; otherwise use dataset-specific default
  if [[ -n "${N_CDS}" ]]; then
    n_cds="${N_CDS}"
  else
    case "${DATASET}" in
      soil)  n_cds=4 ;;
      ocean) n_cds=3 ;;
      *)     n_cds=4 ;;
    esac
  fi

  echo
  echo "=== Dataset: ${DATASET} (n_cds=${n_cds}) ==="
  echo

  start_time=$(date +%s)

  Rscript code/CCA/scripts/08_PCA_compare.R \
    --dataset       "${DATASET}"       \
    --perc_identity "${PERC_IDENTITY}" \
    --n_cds         "${n_cds}"         \
    --n_shuffles    "${N_SHUFFLES}"    \
    --n_cores       "${N_CORES}"       \
    --blas_threads  "${BLAS_T}"

  end_time=$(date +%s)
  elapsed=$(( end_time - start_time ))
  if [[ $elapsed -lt 60 ]]; then
    echo "=== ${DATASET} completed in ${elapsed} seconds ==="
  elif [[ $elapsed -lt 3600 ]]; then
    min=$((elapsed / 60))
    sec=$((elapsed % 60))
    echo "=== ${DATASET} completed in ${elapsed} seconds (${min} min ${sec} sec) ==="
  else
    hr=$((elapsed / 3600))
    min=$(((elapsed % 3600) / 60))
    sec=$((elapsed % 60))
    echo "=== ${DATASET} completed in ${elapsed} seconds (${hr} hr ${min} min ${sec} sec) ==="
  fi
done

conda deactivate

echo "Done."
