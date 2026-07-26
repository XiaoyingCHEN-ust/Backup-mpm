#!/bin/bash
#SBATCH --job-name=mpm_rebuild
#SBATCH --partition=granularmech
#SBATCH --account=comgranmech
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --gres=gpu:1
#SBATCH --time=02:00:00
#SBATCH --mem=16G
#SBATCH --output=rebuild-%j.out
#SBATCH --error=rebuild-%j.err

set -eo pipefail

module load miniconda3/24.3.0-quc3pyu
eval "$("$(command -v conda)" shell.bash hook)"
conda activate cbgeo_tbb
set -u

case_dir="${SLURM_SUBMIT_DIR:?SLURM_SUBMIT_DIR is not set}"
installer="${case_dir}/program_patch/install_hpc4.sh"
if [[ ! -f "${installer}" ]]; then
  echo "Installer is missing: ${installer}" >&2
  echo "Submit rebuild_hpc4.sh from the Siemens validation case directory" >&2
  exit 2
fi
bash "${installer}"
