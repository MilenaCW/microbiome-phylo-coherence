#!/bin/bash
#SBATCH --job-name=gg2_vsearch
#SBATCH --output=../output_files/gg2_vsearch.out
#SBATCH --error=../output_files/gg2_vsearch.err
#SBATCH --account=<account-name>
#SBATCH --partition=<partition-name>
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=4G
#SBATCH --time=3:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=<your-email@domain>

# Source shared setup (strict mode, R_HOME, repo root).
# Under sbatch the script runs from a spool copy so $0 is wrong; use SLURM_SUBMIT_DIR
# (directory from which sbatch was invoked). When not under SLURM, use script location.
# Convention: submit from code/GG2/sbatch so ../../setup.sh resolves to code/setup.sh.
if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
  SETUP_SH="${SLURM_SUBMIT_DIR}/../../setup.sh"
else
  SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
  SETUP_SH="${SCRIPT_DIR}/../../setup.sh"
fi
source "$SETUP_SH"
source "$(dirname "$SETUP_SH")/env_config.sh"

usage() {
  echo "Usage: $0 [OPTIONS]"
  echo ""
  echo "Options:"
  echo "  --config CFG        Path to GG2 config (R list; required)"
  echo "  --step STEP         Run specific step(s). Options:"
  echo "                        - prepare_data       (01a + 01b)"
  echo "                        - vsearch_mapping    (02)"
  echo "                        - build_seqtab       (03a + 03b)"
  echo "                        - tree_and_taxonomy  (04)"
  echo "                        - export_results     (05)"
  echo "                        - clean_outputs      (06)"
  echo "                        - run_diagnostics    (07)"
  echo "                        - summary_diagnostics (08; aggregates 07 across perc-identity)"
  echo "                        - all                (default; runs 01–08)"
  echo "  --perc-identity ID  Percentage identity threshold (default: 0.99)."
  echo "                      Will be normalized to 0.XX (e.g., .95 -> 0.95, 95 -> 0.95)."
  echo "  --help, -h          Show this help message."
  echo ""
  echo "Examples:"
  echo "  sbatch gg2_vsearch.sbatch --config ./code/GG2/config/ocean.R"
  echo "  sbatch gg2_vsearch.sbatch --config ./code/GG2/config/soil.R --step prepare_data"
  echo "  sbatch gg2_vsearch.sbatch --config ./code/GG2/config/soil.R --step all --perc-identity 0.98"
  exit 1
}

normalize_perc_identity() {
  local input="$1"
  input=$(echo "$input" | xargs)
  if [[ "$input" =~ ^\. ]]; then
    input="0$input"
  fi
  if [[ ! "$input" =~ ^0\. ]]; then
    input="0.$input"
  fi
  if [[ "$input" =~ ^0\.[0-9]$ ]]; then
    input="${input}0"
  fi
  echo "$input"
}

STEPS=()
PERC_IDENTITY=0.99
CONFIG_PATH=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --config)
      CONFIG_PATH="$2"
      shift 2
      ;;
    --step)
      STEPS+=("$2")
      shift 2
      ;;
    --perc-identity)
      PERC_IDENTITY=$(normalize_perc_identity "$2")
      shift 2
      ;;
    --help|-h)
      usage
      ;;
    *)
      echo "Unknown option: $1"
      usage
      ;;
  esac
done

if [[ -z "${CONFIG_PATH}" ]]; then
  echo "ERROR: --config <path> is required." >&2
  usage
fi

if [[ ! -f "${CONFIG_PATH}" ]]; then
  echo "ERROR: Config file not found: ${CONFIG_PATH}" >&2
  exit 1
fi

if [ ${#STEPS[@]} -eq 0 ]; then
  STEPS=("all")
fi

print_step_elapsed() {
  local step_name="$1"
  local start="$2"
  local end=$(date +%s)
  local elapsed=$((end - start))
  if [[ $elapsed -ge 3600 ]]; then
    echo "[TIMER] ${step_name} completed in ${elapsed}s ($((elapsed/3600))h $(((elapsed%3600)/60))m $((elapsed%60))s)"
  elif [[ $elapsed -ge 60 ]]; then
    echo "[TIMER] ${step_name} completed in ${elapsed}s ($((elapsed/60))m $((elapsed%60))s)"
  else
    echo "[TIMER] ${step_name} completed in ${elapsed}s"
  fi
  echo
}

should_run_step() {
  local target_step="$1"
  for step in "${STEPS[@]}"; do
    if [[ "$step" == "all" || "$step" == "$target_step" ]]; then
      return 0
    fi
  done
  return 1
}

# adjust as needed for your specific cluster
module load python # required for conda activation
module load R      # provides clean Rscript fallback; prevents MRO from being used on deactivate
N_THREADS=${SLURM_CPUS_PER_TASK:-1}

echo "Using ${N_THREADS} threads and 4G memory/thread"
echo
echo "$(date)"
echo

echo "***********************"
echo "Running GG2 vsearch pipeline"
echo "Config file           : ${CONFIG_PATH}"
echo "Selected steps        : ${STEPS[*]}"
echo "Percentage identity   : ${PERC_IDENTITY}"
echo "***********************"
echo

# Paths from config via functions/access_config.R (dot-path)
GG2_SCRIPTS="./code/GG2/scripts"
ACCESS_CONFIG_R="./code/GG2/functions/access_config.R"
source activate "${PROJECT_ENV}"
OUT_BASE=$(Rscript "${ACCESS_CONFIG_R}" "${CONFIG_PATH}" output.directory | tail -n1 | xargs)
conda deactivate

########################################
## Step 01a + 01b: prepare_data      ##
########################################

if should_run_step "prepare_data"; then
  STEP_START=$(date +%s)
  echo "---- Step 01a/01b: Prepare data and backbone ----"
  echo

  # 01a: backbone export under QIIME2 env
  source activate "${QIIME_ENV}"
  BACKBONE_QZA=$(Rscript "${ACCESS_CONFIG_R}" "${CONFIG_PATH}" greengenes.backbone_qza | tail -n1 | xargs)
  BACKBONE_FASTA="$(dirname "${BACKBONE_QZA}")/backbone_full-length.fasta"
  bash ${GG2_SCRIPTS}/01a_prepare_backbone.sh "${BACKBONE_QZA}" "${BACKBONE_FASTA}"
  conda deactivate

  # 01b: table + sequences reformat under project env
  source activate "${PROJECT_ENV}"
  Rscript ${GG2_SCRIPTS}/01b_prepare_data.R \
    --config "${CONFIG_PATH}" \
    --verbose
  conda deactivate
  print_step_elapsed "01a/01b Prepare data and backbone" "${STEP_START}"
fi

########################################
## Step 02: vsearch_mapping          ##
########################################

if should_run_step "vsearch_mapping"; then
  STEP_START=$(date +%s)
  echo "---- Step 02: vsearch mapping ----"
  echo
  source activate "${PROJECT_ENV}"

  bash ${GG2_SCRIPTS}/02_vsearch_mapping.sh "${CONFIG_PATH}" "${PERC_IDENTITY}"

  conda deactivate
  print_step_elapsed "02 vsearch mapping" "${STEP_START}"
fi

########################################
## Step 03: build GG2 seqtab         ##
########################################

if should_run_step "build_seqtab"; then
  STEP_START=$(date +%s)
  echo "---- Step 03a/03b: Build GG2 seqtab ----"
  echo

  PI_DIR="${OUT_BASE}/${PERC_IDENTITY}"
  INTER_DIR="${PI_DIR}/intermediate"
  EXPORT_DIR="${INTER_DIR}/exports"
  mkdir -p "${EXPORT_DIR}"

  UC_FILE="${EXPORT_DIR}/otu_to_gg2.uc"
  TABLE_IN="${OUT_BASE}/input/table_biom_ready.tsv"
  BACKBONE_QZA=$(Rscript "${ACCESS_CONFIG_R}" "${CONFIG_PATH}" greengenes.backbone_qza | tail -n1 | xargs)
  BACKBONE_FASTA="$(dirname "${BACKBONE_QZA}")/backbone_full-length.fasta"

  # 03a: aggregation under project env
  source activate "${PROJECT_ENV}"
  Rscript ${GG2_SCRIPTS}/03a_build_gg2_seqtab.R \
    --uc-file "${UC_FILE}" \
    --table "${TABLE_IN}" \
    --gg2-backbone "${BACKBONE_FASTA}" \
    --output-dir "${PI_DIR}" \
    --verbose
  conda deactivate

  # 03b: BIOM + QIIME import under QIIME2 env
  source activate "${QIIME_ENV}"
  bash ${GG2_SCRIPTS}/03b_import_to_qiime.sh "${PI_DIR}"
  conda deactivate
  print_step_elapsed "03a/03b Build GG2 seqtab" "${STEP_START}"
fi

########################################
## Step 04: tree_and_taxonomy        ##
########################################

if should_run_step "tree_and_taxonomy"; then
  STEP_START=$(date +%s)
  echo "---- Step 04: Tree pruning and taxonomy ----"
  echo
  source activate "${QIIME_ENV}"
  bash ${GG2_SCRIPTS}/04_tree_and_taxonomy.sh \
    "${CONFIG_PATH}" "${PERC_IDENTITY}"
  conda deactivate
  print_step_elapsed "04 Tree pruning and taxonomy" "${STEP_START}"
fi

########################################
## Step 05: export_results            ##
########################################

if should_run_step "export_results"; then
  STEP_START=$(date +%s)
  echo "---- Step 05: Export tree + taxonomy ----"
  echo
  source activate "${QIIME_ENV}"
  bash ${GG2_SCRIPTS}/05_export_results.sh \
    "${CONFIG_PATH}" "${PERC_IDENTITY}"
  conda deactivate
  print_step_elapsed "05 Export tree + taxonomy" "${STEP_START}"
fi

########################################
## Step 06: clean_outputs             ##
########################################

if should_run_step "clean_outputs"; then
  STEP_START=$(date +%s)
  echo "---- Step 06: Create cleaned seqtab.csv + taxonomy.csv ----"
  echo
  source activate "${PROJECT_ENV}"
  Rscript ${GG2_SCRIPTS}/06_clean_outputs.R \
    --config "${CONFIG_PATH}" \
    --perc-identity "${PERC_IDENTITY}" \
    --verbose
  conda deactivate
  print_step_elapsed "06 Clean outputs" "${STEP_START}"
fi

########################################
## Step 07: run_diagnostics           ##
########################################

if should_run_step "run_diagnostics"; then
  STEP_START=$(date +%s)
  echo "---- Step 07: Run diagnostics ----"
  echo
  source activate "${PROJECT_ENV}"
  # Diagnostics output base mirrors existing pattern: ../../results/gg2_diagnostic/
  Rscript ${GG2_SCRIPTS}/07_gg2_diagnostic.R \
    --config "${CONFIG_PATH}" \
    --perc-identity "${PERC_IDENTITY}" \
    --verbose
  conda deactivate
  print_step_elapsed "07 Run diagnostics" "${STEP_START}"
fi

########################################
## Step 08: summary_diagnostics       ##
########################################

if should_run_step "summary_diagnostics"; then
  STEP_START=$(date +%s)
  echo "---- Step 08: Summary diagnostics (aggregate across perc-identity) ----"
  echo
  source activate "${PROJECT_ENV}"
  Rscript ${GG2_SCRIPTS}/08_gg2_summary_diagnostic.R \
    --config "${CONFIG_PATH}" \
    --verbose
  conda deactivate
  print_step_elapsed "08 Summary diagnostics" "${STEP_START}"
fi

echo "***********************"
echo "GG2 vsearch pipeline completed."
echo
echo "$(date)"

