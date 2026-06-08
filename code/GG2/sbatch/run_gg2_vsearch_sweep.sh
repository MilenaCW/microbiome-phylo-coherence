# Wrapper script to run the generalized GG2 vsearch pipeline with different
# perc-identity values. This script submits multiple SLURM jobs with
# different identity thresholds and custom output/error files for each job.
#
# Usage:
#   bash run_gg2_vsearch_sweep.sh --config config/soil.R [OPTIONS]
#
# Options:
#   --config CFG        Path to GG2 config (R list; required)
#   --start-id ID       Starting perc-identity value (default: 0.90)
#   --end-id ID         Ending perc-identity value (default: 0.99)
#   --step-size STEP    Step size between identity values (default: 0.01)
#   --steps STEPS       Comma-separated list of specific steps to run (default: all)
#                       Available steps: prepare_data, vsearch_mapping, build_seqtab,
#                                        tree_and_taxonomy, export_results,
#                                        clean_outputs, run_diagnostics,
#                                        summary_diagnostics, all
#   --output-dir DIR    Directory for output files (default: ./output_files)
#   --dry-run           Show what would be submitted without actually submitting
#   --help, -h          Show this help message

# Exit immediately if anything fails (prevents silent failure with job still showing as "completed")
set -euo pipefail

usage() {
  sed -n '1,40p' "$0" | sed 's/^#//'
  exit 1
}

CONFIG_PATH=""
START_ID=0.90
END_ID=0.99
STEP_SIZE=0.01
STEPS="all"
OUTPUT_DIR="./output_files"
DRY_RUN=false

while [[ $# -gt 0 ]]; do
  case "$1" in
    --config)
      CONFIG_PATH="$2"
      shift 2
      ;;
    --start-id)
      START_ID="$2"
      shift 2
      ;;
    --end-id)
      END_ID="$2"
      shift 2
      ;;
    --step-size)
      STEP_SIZE="$2"
      shift 2
      ;;
    --steps)
      STEPS="$2"
      shift 2
      ;;
    --output-dir)
      OUTPUT_DIR="$2"
      shift 2
      ;;
    --dry-run)
      DRY_RUN=true
      shift
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

# Config existence is checked only in gg2_vsearch.sbatch (after setup.sh sets CWD to repo root).

if (( $(echo "$START_ID >= $END_ID" | bc -l) )); then
  echo "Error: start-id ($START_ID) must be less than end-id ($END_ID)" >&2
  exit 1
fi

if (( $(echo "$STEP_SIZE <= 0" | bc -l) )); then
  echo "Error: step-size ($STEP_SIZE) must be positive" >&2
  exit 1
fi

# Resolve sbatch script the same way gg2_vsearch.sbatch resolves setup.sh: by script location.
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
SBATCH_SCRIPT="${SCRIPT_DIR}/gg2_vsearch.sbatch"
if [[ ! -f "$SBATCH_SCRIPT" ]]; then
  echo "Error: gg2_vsearch.sbatch not found at ${SBATCH_SCRIPT}" >&2
  exit 1
fi

# Resolve OUTPUT_DIR to absolute path when relative, so it is independent of the cd below.
if [[ "$OUTPUT_DIR" != /* ]]; then
  OUTPUT_DIR="$(pwd)/${OUTPUT_DIR#./}"
fi

# --- Pre-submission: run prepare_data once in the shell ---
# prepare_data (01a + 01b) writes shared files used by all sweep jobs. Running it
# here serially before submitting avoids races on those shared outputs.
if [[ "$STEPS" == "all" ]] || [[ ",$STEPS," == *",prepare_data,"* ]]; then
  ENV_CONFIG="${SCRIPT_DIR}/../../env_config.sh"
  if [[ ! -f "${ENV_CONFIG}" ]]; then
    echo "ERROR: env_config.sh not found at ${ENV_CONFIG}" >&2
    exit 1
  fi
  source "${ENV_CONFIG}"

  # Save SCRIPT_DIR before sourcing setup.sh: setup.sh defines its own SCRIPT_DIR
  # (pointing to code/) which would overwrite ours (pointing to code/GG2/sbatch/).
  SBATCH_DIR="$SCRIPT_DIR"

  # setup.sh sets strict mode and changes CWD to repo root; scripts expect to run from there.
  source "${SBATCH_DIR}/../../setup.sh"

  ACCESS_CONFIG_R="./code/GG2/functions/access_config.R"
  GG2_SCRIPTS="./code/GG2/scripts"

  echo "=========================================="
  echo "Pre-submission: prepare_data"
  echo "=========================================="

  # 01a: backbone export under QIIME2 env
  module load python
  source activate "${QIIME_ENV}"
  BACKBONE_QZA=$(Rscript "${ACCESS_CONFIG_R}" "${CONFIG_PATH}" greengenes.backbone_qza | tail -n1 | xargs)
  BACKBONE_FASTA="$(dirname "${BACKBONE_QZA}")/backbone_full-length.fasta"
  bash "${GG2_SCRIPTS}/01a_prepare_backbone.sh" "${BACKBONE_QZA}" "${BACKBONE_FASTA}"
  conda deactivate

  # 01b: data reformat under project env
  source activate "${PROJECT_ENV}"
  Rscript "${GG2_SCRIPTS}/01b_prepare_data.R" \
    --config "${CONFIG_PATH}" \
    --verbose
  conda deactivate
  module unload python

  echo "Pre-submission prepare_data complete."
  echo "=========================================="
  echo

  # Remove prepare_data from the steps passed to submitted jobs.
  if [[ "$STEPS" == "all" ]]; then
    STEPS="vsearch_mapping,build_seqtab,tree_and_taxonomy,export_results,clean_outputs,run_diagnostics"
  else
    STEPS=$(echo "$STEPS" | tr ',' '\n' | grep -v "^prepare_data$" | tr '\n' ',' | sed 's/,$//')
  fi

  # If prepare_data was the only step, nothing left to submit.
  if [[ -z "$STEPS" ]]; then
    echo "prepare_data was the only requested step; no jobs to submit."
    exit 0
  fi
fi

# Run sbatch from the sbatch directory so SLURM_SUBMIT_DIR is correct (sbatch expects submit from code/GG2/sbatch).
# Use SBATCH_DIR (saved before setup.sh overwrote SCRIPT_DIR) if the pre-submission block ran;
# otherwise SCRIPT_DIR is still the correct value.
cd "${SBATCH_DIR:-$SCRIPT_DIR}"
mkdir -p "$OUTPUT_DIR"

echo "=========================================="
echo "GG2 vsearch Parameter Sweep"
echo "=========================================="
echo "Config file      : $CONFIG_PATH"
echo "Start ID         : $START_ID"
echo "End ID           : $END_ID"
echo "Step size        : $STEP_SIZE"
echo "Steps to run     : $STEPS"
echo "Output directory : $OUTPUT_DIR"
echo "Dry run          : $DRY_RUN"
echo "=========================================="
echo

config_name_formatted=$(basename "${CONFIG_PATH}" .R)

current_id=$START_ID
job_count=0
submitted_jobs=()

while (( $(echo "$current_id <= $END_ID" | bc -l) )); do
  perc_id_formatted=$(printf "%.2f" "$current_id")

  output_file="$OUTPUT_DIR/gg2_${perc_id_formatted}.out"
  error_file="$OUTPUT_DIR/gg2_${perc_id_formatted}.err"
  job_name="gg2_vsearch_${perc_id_formatted}_${config_name_formatted}"

  sbatch_cmd="sbatch"
  sbatch_cmd="$sbatch_cmd --output=$output_file"
  sbatch_cmd="$sbatch_cmd --error=$error_file"
  sbatch_cmd="$sbatch_cmd --job-name=$job_name"
  sbatch_cmd="$sbatch_cmd $SBATCH_SCRIPT"
  sbatch_cmd="$sbatch_cmd --config $CONFIG_PATH"
  sbatch_cmd="$sbatch_cmd --perc-identity $current_id"

  if [ "$STEPS" != "all" ]; then
    IFS=',' read -ra STEP_ARRAY <<< "$STEPS"
    for step in "${STEP_ARRAY[@]}"; do
      sbatch_cmd="$sbatch_cmd --step $step"
    done
  fi

  echo "Job $((++job_count)): perc-identity=$current_id"
  echo "  Output : $output_file"
  echo "  Error  : $error_file"
  echo "  Command: $sbatch_cmd"
  echo

  if [ "$DRY_RUN" = false ]; then
    job_id=$(eval "$sbatch_cmd")
    submitted_jobs+=("$job_id")
    echo "  Submitted as job: $job_id"
  else
    echo "  [DRY RUN] Would submit: $sbatch_cmd"
  fi
  echo

  current_id=$(echo "$current_id + $STEP_SIZE" | bc -l)
done

echo "=========================================="
if [ "$DRY_RUN" = false ]; then
  echo "Submitted $job_count jobs:"
  for job_id in "${submitted_jobs[@]}"; do
    echo "  $job_id"
  done
  echo ""
  echo "Monitor jobs with:"
  echo "  squeue -u \$USER"
  echo ""
  echo "Check job status with:"
  echo "  sacct"
else
  echo "DRY RUN: Would have submitted $job_count jobs"
fi
echo "=========================================="

