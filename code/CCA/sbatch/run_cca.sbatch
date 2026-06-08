#!/bin/bash
#SBATCH --job-name=cca
#SBATCH --output=./output_files/cca.out
#SBATCH --error=./output_files/cca.err
#SBATCH --account=<account-name>
#SBATCH --partition=<partition-name>
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem-per-cpu=4G
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=<your-email@domain>

# CCA pipeline. Pass flags after the script name; sbatch forwards them to the script.
#
# Usage:
#   sbatch run_cca.sbatch --config code/CCA/config/soil.R [OPTIONS]
#
# Options:
#   --config CFG         Path to CCA config (R list; default: code/CCA/config/soil.R)
#   --tax-level LEVEL    Taxonomy level, e.g. OTU, Genus (default: OTU)
#   --perc-identity ID   Percent identity threshold (default: 0.90)
#   --step STEP          Run specific step(s). May be repeated. (default: all)
#                        Steps: 00_read_data, 01_hyperparam_search, 02_loadings,
#                               03_null_distribution, 04_model_performance_plots, all
#   --help, -h           Show this help message
#
# Examples:
#   sbatch run_cca.sbatch --config code/CCA/config/soil.R
#   sbatch run_cca.sbatch --config code/CCA/config/soil.R --perc-identity 0.97 --tax-level Genus
#   sbatch run_cca.sbatch --config code/CCA/config/soil.R --step 00_read_data --step 01_hyperparam_search

set -euo pipefail

usage() {
  sed -n '13,36p' "$0" | sed 's/^#//'
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
CONFIG="code/CCA/config/soil.R"
TAX_LEVEL="OTU"
PERC_IDENTITY="0.90"
STEPS=()

while [[ $# -gt 0 ]]; do
  case "$1" in
    --config)
      CONFIG="$2"
      shift 2
      ;;
    --tax-level)
      TAX_LEVEL="$2"
      shift 2
      ;;
    --perc-identity)
      PERC_IDENTITY="$2"
      shift 2
      ;;
    --step)
      STEPS+=("$2")
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

if [[ ${#STEPS[@]} -eq 0 ]]; then
  STEPS=("00_read_data" "01_hyperparam_search" "02_loadings" "03_null_distribution" "04_model_performance_plots")
else
  # Expand "all" to every step
  if printf '%s\n' "${STEPS[@]}" | grep -qx "all"; then
    STEPS=("00_read_data" "01_hyperparam_search" "02_loadings" "03_null_distribution" "04_model_performance_plots")
  fi
fi

if [[ ! -f "${CONFIG}" ]]; then
  echo "ERROR: Config file not found: ${CONFIG}" >&2
  exit 1
fi

echo "CCA pipeline running with parameters:"
echo "    config=$CONFIG"
echo "    tax_level=$TAX_LEVEL"
echo "    perc_identity=$PERC_IDENTITY"
echo "    steps=${STEPS[*]}"
echo "    SLURM_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK:-}"

# BLAS/OpenMP thread counts: parallel steps (01, 03) vs others (defaults; override by exporting before sbatch)
# these are both capped at the SLURM_CPUS_PER_TASK
normal_BLAS_T="${normal_BLAS_T:-${SLURM_CPUS_PER_TASK:-1}}"
parallel_BLAS_T="${parallel_BLAS_T:-1}"
if [[ -n "${SLURM_CPUS_PER_TASK:-}" ]]; then
  [[ "$normal_BLAS_T" -gt "$SLURM_CPUS_PER_TASK" ]] && normal_BLAS_T="$SLURM_CPUS_PER_TASK"
  [[ "$parallel_BLAS_T" -gt "$SLURM_CPUS_PER_TASK" ]] && parallel_BLAS_T="$SLURM_CPUS_PER_TASK"
fi

run_step() {
  local step="$1"
  echo
  echo "=== Step $step ==="
  echo

  if [[ "$step" == "01_hyperparam_search" || "$step" == "03_null_distribution" ]]; then
    export OMP_NUM_THREADS="${parallel_BLAS_T}"
    export MKL_NUM_THREADS="${parallel_BLAS_T}"
    export OPENBLAS_NUM_THREADS="${parallel_BLAS_T}"
    export VECLIB_MAXIMUM_THREADS="${parallel_BLAS_T}"
    export NUMEXPR_NUM_THREADS="${parallel_BLAS_T}"
  else
    export OMP_NUM_THREADS="${normal_BLAS_T}"
    export MKL_NUM_THREADS="${normal_BLAS_T}"
    export OPENBLAS_NUM_THREADS="${normal_BLAS_T}"
    export VECLIB_MAXIMUM_THREADS="${normal_BLAS_T}"
    export NUMEXPR_NUM_THREADS="${normal_BLAS_T}"
  fi

  echo "Thread env (SLURM export at step $step):"
  env | egrep 'OMP_NUM_THREADS|MKL_NUM_THREADS|OPENBLAS_NUM_THREADS|VECLIB_MAXIMUM_THREADS|NUMEXPR_NUM_THREADS' || true
  echo

  local start_time end_time elapsed
  start_time=$(date +%s)

  local rscript_args=(--config "$CONFIG" --tax_level "$TAX_LEVEL" --perc_identity "$PERC_IDENTITY" --verbose)

  # number of cores for parallel processing * the number of BLAS threads for internal multithreading cannot be more than the SLURM_CPUS_PER_TASK
  if [[ ("$step" == "01_hyperparam_search" || "$step" == "03_null_distribution") && -n "${SLURM_CPUS_PER_TASK:-}" ]]; then
    local n_cores=$(( SLURM_CPUS_PER_TASK / parallel_BLAS_T ))
    [[ $n_cores -lt 1 ]] && n_cores=1
    rscript_args+=(--n_cores "$n_cores")
    echo "  n_cores for step $step: $n_cores"
  fi

  Rscript "code/CCA/scripts/${step}.R" "${rscript_args[@]}"

  end_time=$(date +%s)
  elapsed=$(( end_time - start_time ))
  if [[ $elapsed -lt 60 ]]; then
    echo "=== Step $step completed in ${elapsed} seconds ==="
  elif [[ $elapsed -lt 3600 ]]; then
    min=$((elapsed / 60))
    sec=$((elapsed % 60))
    echo "=== Step $step completed in ${elapsed} seconds (${min} min ${sec} sec) ==="
  else
    hr=$((elapsed / 3600))
    min=$(((elapsed % 3600) / 60))
    sec=$((elapsed % 60))
    echo "=== Step $step completed in ${elapsed} seconds (${hr} hr ${min} min ${sec} sec) ==="
  fi
}

module load python # adjust as needed for your specific cluster
source activate "${PROJECT_ENV}"

for step in "${STEPS[@]}"; do
  run_step "$step"
done

conda deactivate

echo "Done."
