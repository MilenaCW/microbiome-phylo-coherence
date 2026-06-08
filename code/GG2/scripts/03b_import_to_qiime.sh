# 03b_import_to_qiime.sh
# ----------------------
# Convert GG2 feature-table.tsv to BIOM and import into QIIME2.
# Environment activation is handled by the calling sbatch script; this script
# assumes the correct QIIME2 environment is already active.
#
# Usage:
#   bash 03b_import_to_qiime.sh <pi_dir>
#
# Arguments:
#   pi_dir  PI-specific directory (e.g., <output$directory>/0.99)

# Source shared setup (strict mode, R_HOME, repo root). Path is relative to this script.
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
source "${SCRIPT_DIR}/../../setup.sh"

if [[ "$#" -ne 1 ]]; then
  echo "Usage: $0 <pi_dir>" >&2
  exit 1
fi

PI_DIR="$1"
EXPORT_DIR="${PI_DIR}/intermediate/exports"

FEATURE_TSV="${EXPORT_DIR}/feature-table.tsv"
BIOM_OUT="${PI_DIR}/intermediate/gg2_mapped_table.biom"
QZA_OUT="${PI_DIR}/intermediate/gg2_mapped_table.qza"

if [[ ! -f "${FEATURE_TSV}" ]]; then
  echo "[03b_import_to_qiime] ERROR: feature-table.tsv not found at ${FEATURE_TSV}" >&2
  exit 1
fi

mkdir -p "$(dirname "${BIOM_OUT}")"

echo "[03b_import_to_qiime] Converting feature-table.tsv to BIOM..."
biom convert \
  -i "${FEATURE_TSV}" \
  -o "${BIOM_OUT}" \
  --to-hdf5

echo "[03b_import_to_qiime] Importing BIOM table into QIIME2..."
qiime tools import \
  --input-path "${BIOM_OUT}" \
  --output-path "${QZA_OUT}" \
  --type 'FeatureTable[Frequency]'

echo "[03b_import_to_qiime] Completed BIOM conversion and QIIME2 import."

