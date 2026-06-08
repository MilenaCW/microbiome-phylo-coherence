# env_config.sh
# --------------
# User-specific conda environment paths. Set CONDA_ENVS_PATH to the directory
# containing your conda environments (e.g. /home/user/miniconda3/envs).
# Convention: submit sbatch jobs from code/<pipeline>/sbatch/ so that
# ../../env_config.sh resolves to this file.

export CONDA_ENVS_PATH="<path-to-envs>"
export QIIME_ENV="${CONDA_ENVS_PATH}/qiime2-amplicon-2025.7"
export PROJECT_ENV="${CONDA_ENVS_PATH}/microbiome-phylo-coherence"
