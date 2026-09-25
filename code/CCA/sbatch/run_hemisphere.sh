#!/bin/bash
# Submit the soil CCA pipeline (OTU level) separately for each hemisphere.
# Results go to soil/results/CCA_hemisphere/<hemisphere>/ (existing soil/results/CCA is untouched).
#
# Usage (from code/CCA/sbatch):
#   bash run_hemisphere.sh [--hemispheres north,south] [--perc-identity 0.90] [--dry-run]
#
# Runs north and south by default. The full-data run is not part of this script: use the
# regular pipeline (run_cca.sbatch --config code/CCA/config/soil.R) with its usual output location.
# Steps 01 and 03 dominate runtime; raise --time in run_cca.sbatch if needed.

set -euo pipefail

HEMISPHERES="north,south"
PERC_IDENTITY="0.90"
DRY_RUN=false

while [[ $# -gt 0 ]]; do
  case "$1" in
    --hemispheres) HEMISPHERES="$2"; shift 2 ;;
    --perc-identity) PERC_IDENTITY="$2"; shift 2 ;;
    --dry-run) DRY_RUN=true; shift ;;
    --help|-h) sed -n '2,11p' "$0" | sed 's/^# \{0,1\}//'; exit 1 ;;
    *) echo "Unknown option: $1" >&2; exit 1 ;;
  esac
done

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"
mkdir -p output_files

IFS=',' read -ra HEMI_ARRAY <<< "$HEMISPHERES"
for h in "${HEMI_ARRAY[@]}"; do
  cmd=(sbatch --output="./output_files/cca_${h}.out" --error="./output_files/cca_${h}.err" --job-name="cca_${h}"
       run_cca.sbatch --config "code/CCA/config/soil_${h}.R" --tax-level OTU --perc-identity "$PERC_IDENTITY")
  echo "${cmd[*]}"
  if [[ "$DRY_RUN" == false ]]; then "${cmd[@]}"; fi
done
