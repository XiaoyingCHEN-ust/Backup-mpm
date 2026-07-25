#!/bin/bash
#SBATCH --job-name=siemens3p
#SBATCH --partition=granularmech
#SBATCH --account=comgranmech
#SBATCH --array=0-1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --gres=gpu:1
#SBATCH --time=24:00:00
#SBATCH --mem=16G
#SBATCH --output=slurm-%A_%a.out
#SBATCH --error=slurm-%A_%a.err

set -eo pipefail

case_dir="${SLURM_SUBMIT_DIR:?SLURM_SUBMIT_DIR is not set}"
cd "${case_dir}"

module load miniconda3/24.3.0-quc3pyu
eval "$("$(command -v conda)" shell.bash hook)"
conda activate cbgeo_tbb

mpm_bin="${MPM_BIN:-/home/xchenjm/chen/MPM/mpm/build/mpm}"
cases=(open closed)
case_name="${cases[${SLURM_ARRAY_TASK_ID:-0}]}"
threads="${SLURM_CPUS_PER_TASK:-32}"

if [[ ! -x "${mpm_bin}" ]]; then
  echo "MPM executable is missing or not executable: ${mpm_bin}" >&2
  exit 126
fi

echo "Running Siemens et al. (2013) ${case_name} column test"
echo "Node: $(hostname)"
echo "Executable: ${mpm_bin}"
echo "Case directory: ${case_dir}"
echo "Threads: ${threads}"
nvidia-smi || true

"${mpm_bin}" \
  -f "${case_dir}/" \
  -i "mpm_${case_name}.json" \
  -p "${threads}"
