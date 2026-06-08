# 02_vsearch_mapping.sh
# ----------------------
# Run vsearch closed-reference mapping of input features/ASVs to the GG2 backbone.
# Environment activation is handled by the calling sbatch script; this script
# assumes the correct vsearch-capable environment is already active.
#
# Usage:
#   bash 02_vsearch_mapping.sh <config_R> <perc_identity>
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
  echo "[02_vsearch_mapping] ERROR: Config file not found: ${CONFIG_R}" >&2
  exit 1
fi

# Read paths from config via functions/access_config.R (dot-path)
ACCESS_CONFIG_R="./code/GG2/functions/access_config.R"
# Only the last line of stdout is the config value; setwd messages go to stderr (access_config.R)
OUT_BASE=$(Rscript "${ACCESS_CONFIG_R}" "${CONFIG_R}" output.directory | tail -n1 | xargs)
BACKBONE_QZA=$(Rscript "${ACCESS_CONFIG_R}" "${CONFIG_R}" greengenes.backbone_qza | tail -n1 | xargs)
INPUT_SEQS="${OUT_BASE}/input/sequences_filtered.fna"

if [[ ! -f "${INPUT_SEQS}" ]]; then
  echo "[02_vsearch_mapping] ERROR: Input sequences not found: ${INPUT_SEQS}" >&2
  exit 1
fi

# 01a has already exported backbone_full-length.fasta
BACKBONE_FASTA_DIR="$(dirname "${BACKBONE_QZA}")"
BACKBONE_FASTA="${BACKBONE_FASTA_DIR}/backbone_full-length.fasta"

if [[ ! -f "${BACKBONE_FASTA}" ]]; then
  echo "[02_vsearch_mapping] ERROR: Expected backbone FASTA at: ${BACKBONE_FASTA}" >&2
  echo "Did you run 01a_prepare_backbone.sh?" >&2
  exit 1
fi

PI_DIR="${OUT_BASE}/${PERC_ID}"
EXPORT_DIR="${PI_DIR}/intermediate/exports"
mkdir -p "${EXPORT_DIR}"

OUT_UC="${EXPORT_DIR}/otu_to_gg2.uc"

echo "[02_vsearch_mapping] Mapping ${INPUT_SEQS} to GG2 backbone at ${BACKBONE_FASTA}"
echo "[02_vsearch_mapping] Output UC: ${OUT_UC}"

vsearch \
  --usearch_global "${INPUT_SEQS}" \
  --db "${BACKBONE_FASTA}" \
  --id "${PERC_ID}" \
  --strand plus \
  --top_hits_only \
  --uc "${OUT_UC}" \
  --threads "${SLURM_CPUS_PER_TASK:-1}"

echo "[02_vsearch_mapping] Completed vsearch mapping."

