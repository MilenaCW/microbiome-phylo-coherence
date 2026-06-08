# 04_tree_and_taxonomy.sh
# ------------------------
# Prune the GG2 tree to mapped features and assign taxonomy using the GG2
# reference taxonomy. Environment activation is handled by the calling
# sbatch script; this script assumes the correct QIIME2 environment is
# already active.
#
# Usage:
#   bash 04_tree_and_taxonomy.sh <config_R> <perc_identity>
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
  echo "[04_tree_and_taxonomy] ERROR: Config file not found: ${CONFIG_R}" >&2
  exit 1
fi

# Read paths from config via functions/access_config.R (dot-path)
ACCESS_CONFIG_R="./code/GG2/functions/access_config.R"
OUT_BASE=$(Rscript "${ACCESS_CONFIG_R}" "${CONFIG_R}" output.directory | tail -n1 | xargs)
TREE_QZA=$(Rscript "${ACCESS_CONFIG_R}" "${CONFIG_R}" greengenes.tree_qza | tail -n1 | xargs)
TAX_QZA=$(Rscript "${ACCESS_CONFIG_R}" "${CONFIG_R}" greengenes.taxonomy_qza | tail -n1 | xargs)

PI_DIR="${OUT_BASE}/${PERC_ID}"
INTER_DIR="${PI_DIR}/intermediate"

TABLE_QZA="${INTER_DIR}/gg2_mapped_table.qza"
PRUNED_TREE_QZA="${INTER_DIR}/gg2_pruned_tree.qza"
OUT_TAX_QZA="${INTER_DIR}/gg2_taxonomy.qza"

if [[ ! -f "${TREE_QZA}" ]]; then
  echo "[04_tree_and_taxonomy] ERROR: GG2 tree artifact not found: ${TREE_QZA}" >&2
  exit 1
fi
if [[ ! -f "${TAX_QZA}" ]]; then
  echo "[04_tree_and_taxonomy] ERROR: GG2 taxonomy artifact not found: ${TAX_QZA}" >&2
  exit 1
fi
if [[ ! -f "${TABLE_QZA}" ]]; then
  echo "[04_tree_and_taxonomy] ERROR: gg2_mapped_table.qza not found at ${TABLE_QZA}" >&2
  exit 1
fi

mkdir -p "${INTER_DIR}"

echo "[04_tree_and_taxonomy] Pruning GG2 tree..."
qiime phylogeny filter-tree \
  --i-tree "${TREE_QZA}" \
  --i-table "${TABLE_QZA}" \
  --o-filtered-tree "${PRUNED_TREE_QZA}"

echo "[04_tree_and_taxonomy] Assigning taxonomy from GG2 reference..."
qiime greengenes2 taxonomy-from-table \
  --i-table "${TABLE_QZA}" \
  --i-reference-taxonomy "${TAX_QZA}" \
  --o-classification "${OUT_TAX_QZA}"

echo "[04_tree_and_taxonomy] Completed tree pruning and taxonomy assignment."

