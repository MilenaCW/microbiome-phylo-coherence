# Wrapper script to run the CCA pipeline with different tax_level values.
# Submits multiple SLURM jobs (one per tax level) with custom output/error
# files for each job.
#
# Usage:
#   bash run_cca_tax_sweep.sh --config ./code/CCA/config/soil.R [OPTIONS]
#
# Options:
#   --config CFG        Path to CCA config (from repo root; required)
#   --perc-identity ID  Percent identity threshold (default: 0.90)
#   --levels LEVELS     Comma-separated tax levels to run (default: all)
#                       Use "all" for all levels: Phylum, Class, Order,
#                       Family, Genus, Species, OTU
#   --steps STEPS       Comma-separated list of specific steps to run (default: all)
#                       Steps: 00_read_data, 01_hyperparam_search, 02_loadings,
#                              03_null_distribution, 04_model_performance_plots, all
#   --output-dir DIR    Directory for output files (default: ./output_files)
#   --dry-run           Show what would be submitted without actually submitting
#   --help, -h          Show this help message

set -euo pipefail

usage() {
  sed -n '1,35p' "$0" | sed 's/^#//'
  exit 1
}

# All allowed tax levels (must match CCA scripts).
ALL_LEVELS="Phylum,Class,Order,Family,Genus,Species,OTU"

CONFIG_PATH=""
PERC_IDENTITY="0.90"
LEVELS="all"
STEPS="all"
OUTPUT_DIR="./output_files"
DRY_RUN=false

while [[ $# -gt 0 ]]; do
  case "$1" in
    --config)
      CONFIG_PATH="$2"
      shift 2
      ;;
    --perc-identity)
      PERC_IDENTITY="$2"
      shift 2
      ;;
    --levels)
      LEVELS="$2"
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

# Build list of levels to run.
if [[ "$LEVELS" == "all" ]]; then
  IFS=',' read -ra LEVEL_ARRAY <<< "$ALL_LEVELS"
else
  IFS=',' read -ra LEVEL_ARRAY <<< "$LEVELS"
fi

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
SBATCH_SCRIPT="${SCRIPT_DIR}/run_cca.sbatch"
if [[ ! -f "$SBATCH_SCRIPT" ]]; then
  echo "Error: run_cca.sbatch not found at ${SBATCH_SCRIPT}" >&2
  exit 1
fi

if [[ "$OUTPUT_DIR" != /* ]]; then
  OUTPUT_DIR="$(pwd)/${OUTPUT_DIR#./}"
fi

cd "$SCRIPT_DIR"
mkdir -p "$OUTPUT_DIR"

echo "=========================================="
echo "CCA Tax Level Sweep"
echo "=========================================="
echo "Config file      : $CONFIG_PATH"
echo "Perc identity    : $PERC_IDENTITY"
echo "Levels           : ${LEVEL_ARRAY[*]}"
echo "Steps to run     : $STEPS"
echo "Output directory : $OUTPUT_DIR"
echo "Dry run          : $DRY_RUN"
echo "=========================================="
echo

job_count=0
submitted_jobs=()

trim() { echo "$1" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//'; }

for level in "${LEVEL_ARRAY[@]}"; do
  level_trimmed=$(trim "$level")
  [[ -z "$level_trimmed" ]] && continue

  output_file="$OUTPUT_DIR/cca_${level_trimmed}.out"
  error_file="$OUTPUT_DIR/cca_${level_trimmed}.err"
  job_name="cca_${level_trimmed}"

  sbatch_cmd="sbatch"
  sbatch_cmd="$sbatch_cmd --output=$output_file"
  sbatch_cmd="$sbatch_cmd --error=$error_file"
  sbatch_cmd="$sbatch_cmd --job-name=$job_name"
  sbatch_cmd="$sbatch_cmd $SBATCH_SCRIPT"
  sbatch_cmd="$sbatch_cmd --config $CONFIG_PATH"
  sbatch_cmd="$sbatch_cmd --perc-identity $PERC_IDENTITY"
  sbatch_cmd="$sbatch_cmd --tax-level $level_trimmed"

  if [[ "$STEPS" != "all" ]]; then
    IFS=',' read -ra STEP_ARRAY <<< "$STEPS"
    for step in "${STEP_ARRAY[@]}"; do
      step_trimmed=$(trim "$step")
      [[ -n "$step_trimmed" ]] && sbatch_cmd="$sbatch_cmd --step $step_trimmed"
    done
  fi

  echo "Job $((++job_count)): tax-level=$level_trimmed"
  echo "  Output : $output_file"
  echo "  Error  : $error_file"
  echo "  Command: $sbatch_cmd"
  echo

  if [[ "$DRY_RUN" == false ]]; then
    job_id=$(eval "$sbatch_cmd")
    submitted_jobs+=("$job_id")
    echo "  Submitted as job: $job_id"
  else
    echo "  [DRY RUN] Would submit: $sbatch_cmd"
  fi
  echo
done

echo "=========================================="
if [[ "$DRY_RUN" == false ]]; then
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
