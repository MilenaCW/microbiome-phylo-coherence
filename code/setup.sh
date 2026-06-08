# setup.sh
# ---------
# Shared setup for scripts under code/. Source this file using a path
# relative to the calling script (e.g. from GG2/scripts: source ../../setup.sh).
#
# - Sets strict mode (set -euo pipefail).
# - Ensures rpy2 (used by q2_composition) finds conda's R when applicable.
# - Sets working directory to the repo root (edit the paths below for your machine).
# - Exports REPO_ROOT for child processes.

set -euo pipefail

# Ensure rpy2 (used by q2_composition) finds conda's R and its base libraries.
# Without this, "qiime tools export" can fail with: shared object 'methods.dylib' not found
if [[ -n "${CONDA_PREFIX:-}" && -d "${CONDA_PREFIX}/lib/R" ]]; then
  export R_HOME="${CONDA_PREFIX}/lib/R"
fi

# Repo root: derive from this script's own location (code/ is one level below root).
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
REPO_ROOT="$(dirname "$SCRIPT_DIR")"
cd "$REPO_ROOT"
export REPO_ROOT
echo "[setup.sh] Working directory set to: ${REPO_ROOT}" >&2
