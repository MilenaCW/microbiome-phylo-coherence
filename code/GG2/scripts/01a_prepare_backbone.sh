# 01a_prepare_backbone.sh
# ------------------------
# Export GG2 backbone sequences from a QIIME2 artifact if they do not
# already exist as a plain FASTA file. Environment activation is handled
# by the calling sbatch script; this script assumes the correct QIIME2
# environment is already active.
#
# Usage:
#   bash 01a_prepare_backbone.sh <backbone_qza> <output_fasta>
#
# Arguments:
#   backbone_qza   Path to GG2 backbone QZA (e.g., 2024.09.backbone.full-length.fna.qza)
#   output_fasta   Path to write backbone_full-length.fasta

# Source shared setup (strict mode, R_HOME, repo root). Path is relative to this script.
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
source "${SCRIPT_DIR}/../../setup.sh"

# Confirm that there are two inputs given (if not throw error and print proper usage)
if [[ "$#" -ne 2 ]]; then
  echo "Usage: $0 <backbone_qza> <output_fasta>" >&2
  exit 1 # note: exit 1 indicates an error has occurred
fi

BACKBONE_QZA="$1"
OUT_FASTA="$2"

# Check if the backbone QZA exists, if not throw an error
if [[ ! -f "${BACKBONE_QZA}" ]]; then
  echo "[01a_prepare_backbone] ERROR: backbone QZA not found: ${BACKBONE_QZA}" >&2
  exit 1
fi

# pull out the directory name for the output FASTA path and make that directory if it doesn't already exist
OUT_DIR="$(dirname "${OUT_FASTA}")"
mkdir -p "${OUT_DIR}"

LOCK_FILE="${OUT_DIR}/.backbone_export.lock"

(
  flock -x 200

  if [[ -f "${OUT_FASTA}" ]]; then
    echo "[01a_prepare_backbone] Backbone FASTA already exists at ${OUT_FASTA}; skipping export."
    exit 0
  fi

  echo "[01a_prepare_backbone] Exporting GG2 backbone from QIIME2 artifact..."
  qiime tools export \
    --input-path "${BACKBONE_QZA}" \
    --output-path "${OUT_DIR}"

  if [[ ! -f "${OUT_DIR}/dna-sequences.fasta" ]]; then
    echo "[01a_prepare_backbone] ERROR: dna-sequences.fasta not found in export directory ${OUT_DIR}" >&2
    exit 1
  fi

  mv "${OUT_DIR}/dna-sequences.fasta" "${OUT_FASTA}"
  echo "[01a_prepare_backbone] Wrote backbone FASTA to ${OUT_FASTA}"

) 200>"${LOCK_FILE}"