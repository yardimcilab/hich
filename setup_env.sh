CONDA_BASE=$(conda info --base)
SNAKEMAKE_BASE=$(pwd)
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate "${SNAKEMAKE_BASE}/resources/snakemake"
