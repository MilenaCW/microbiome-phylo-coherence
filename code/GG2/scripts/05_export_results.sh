# 05_export_results.sh
# ---------------------
# Export the pruned GG2 tree and taxonomy from QIIME2 artifacts into
# plain-text formats under the PI/intermediate/exports/ directory.
# Environment activation is handled by the calling sbatch script; this
# script assumes the correct QIIME2 environment is already active.
#
# Usage:
#   bash 05_export_results.sh <config_R> <perc_identity>
#
# Arguments:
#   config_R      Path to GG2 config (R list)
#   perc_identity Percentage identity (e.g., 0.99)

# Source shared setup (strict mode, R_HOME, repo root). Path is relative to this script.
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
source "${SCRIPT_DIR}/../../setup.sh"

if [[ "$#" -ne 2 ]]; then
  echo "Usage: $0 <config_R> <perc_identity>" >&2
  exit 1
fi

CONFIG_R="$1"
PERC_ID="$2"

if [[ ! -f "${CONFIG_R}" ]]; then
  echo "[05_export_results] ERROR: Config file not found: ${CONFIG_R}" >&2
  exit 1
fi

# Read path from config via functions/access_config.R (dot-path)
ACCESS_CONFIG_R="./code/GG2/functions/access_config.R"
OUT_BASE=$(Rscript "${ACCESS_CONFIG_R}" "${CONFIG_R}" output.directory | tail -n1 | xargs)
PI_DIR="${OUT_BASE}/${PERC_ID}"
INTER_DIR="${PI_DIR}/intermediate"
EXPORT_DIR="${INTER_DIR}/exports"

PRUNED_TREE_QZA="${INTER_DIR}/gg2_pruned_tree.qza"
TAX_QZA="${INTER_DIR}/gg2_taxonomy.qza"

if [[ ! -f "${PRUNED_TREE_QZA}" ]]; then
  echo "[05_export_results] ERROR: gg2_pruned_tree.qza not found at ${PRUNED_TREE_QZA}" >&2
  exit 1
fi
if [[ ! -f "${TAX_QZA}" ]]; then
  echo "[05_export_results] ERROR: gg2_taxonomy.qza not found at ${TAX_QZA}" >&2
  exit 1
fi

mkdir -p "${EXPORT_DIR}"

echo "[05_export_results] Exporting pruned tree to Newick..."
qiime tools export \
  --input-path "${PRUNED_TREE_QZA}" \
  --output-path "${EXPORT_DIR}"

echo "[05_export_results] Exporting taxonomy to TSV..."
qiime tools export \
  --input-path "${TAX_QZA}" \
  --output-path "${EXPORT_DIR}"

echo "[05_export_results] Completed export of tree and taxonomy."

