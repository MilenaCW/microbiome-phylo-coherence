#!/bin/bash
#SBATCH --job-name=soil_dada
#SBATCH --output=./output_files/soil_dada_%a.out
#SBATCH --error=./output_files/soil_dada_%a.err
#SBATCH --account=<account-name>
#SBATCH --partition=<partition-name>
#SBATCH --array=1-6
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=4G
#SBATCH --time=36:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=<your-email@domain>

# =============================================================================
# DADA2 denoising pipeline — soil dataset only.
# Ocean dataset uses miTAG sequences (pre-denoised); DADA2 is not run for ocean.
# =============================================================================
# Array Jobs: 6 tasks (1-6), one per flowcell
# Resources per task:
#   - CPUs: 32 cores
#   - Memory: 128 GB total (4 GB per CPU)
#   - Time limit: 36 hours
# =============================================================================

# Load required modules
module load python # adjust as needed for your specific cluster

# Source environment config (sets PROJECT_ENV; edit code/env_config.sh once per user).
# Convention: submit from code/DADA/sbatch so ../../env_config.sh resolves to code/env_config.sh.
if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
  ENV_CONFIG_SH="${SLURM_SUBMIT_DIR}/../../env_config.sh"
else
  ENV_CONFIG_SH="$(cd "$(dirname "$0")" && pwd)/../../env_config.sh"
fi
source "$ENV_CONFIG_SH"
source activate "${PROJECT_ENV}"

# Create output directory if it doesn't exist
mkdir -p ./output_files

# Function to print timestamped messages
print_timestamp() {
    echo "=========================================="
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"
    echo "=========================================="
}

# Function to run R script with error checking
run_r_script() {
    local script_name=$1
    local flow_cell=$2
    local config_path=$3
    local description=$4

    print_timestamp "$description"
    if Rscript "$script_name" --config "$config_path" --flowcell "$flow_cell"; then
        echo "✓ Successfully completed: $description"
    else
        echo "✗ Failed: $description"
        exit 1
    fi
}

# Flowcell IDs are specific to this dataset. To adapt for a different dataset,
# sort your FASTQ files by flowcell and update the array size (#SBATCH --array)
# and case entries below accordingly.
case $SLURM_ARRAY_TASK_ID in
    1)
        print_timestamp "Processing flow cell: AA37J"
        flow_cell="AA37J"
        ;;
    2)
        print_timestamp "Processing flow cell: AA3RF"
        flow_cell="AA3RF"
        ;;
    3)
        print_timestamp "Processing flow cell: AK2PM"
        flow_cell="AK2PM"
        ;;
    4)
        print_timestamp "Processing flow cell: AN72A"
        flow_cell="AN72A"
        ;;
    5)
        print_timestamp "Processing flow cell: ATY9B"
        flow_cell="ATY9B"
        ;;
    6)
        print_timestamp "Processing flow cell: AWR5C"
        flow_cell="AWR5C"
        ;;
    *)
        echo "Error: Invalid array task ID: $SLURM_ARRAY_TASK_ID"
        exit 1
        ;;
esac

config_path="code/DADA/config/soil.R"

print_timestamp "Starting DADA2 pipeline for flow cell: $flow_cell"

# Run DADA2 pipeline steps
run_r_script "../01_plotQualityProfile.R" "$flow_cell" "$config_path" "Plotting quality profiles"
run_r_script "../02_filterAndTrim.R" "$flow_cell" "$config_path" "Filtering and trimming reads"
run_r_script "../03_learnError.R" "$flow_cell" "$config_path" "Learning error model"
run_r_script "../04_sampleInference.R" "$flow_cell" "$config_path" "Denoising and inferring sequences"
run_r_script "../05_merge.R" "$flow_cell" "$config_path" "Merging forward and reverse reads"
run_r_script "../06_buildTable.R" "$flow_cell" "$config_path" "Building sequence table"
run_r_script "../07_removeChimera.R" "$flow_cell" "$config_path" "Removing chimeras"
run_r_script "../08a_dada_diagnostics.R" "$flow_cell" "$config_path" "Creating dada diagnostic plots"
run_r_script "../08b_seqtab_diagnostics.R" "$flow_cell" "$config_path" "Creating sequence table diagnostic plots"

print_timestamp "DADA2 pipeline completed successfully for flow cell: $flow_cell"
