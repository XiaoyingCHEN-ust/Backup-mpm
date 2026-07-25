#!/bin/bash
#SBATCH --job-name=siemens3p
#SBATCH --array=0-1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --time=24:00:00
#SBATCH --mem=16G
#SBATCH --output=slurm-%A_%a.out

set -euo pipefail

if [[ -z "${MPM_BIN:-}" ]]; then
  echo "Set MPM_BIN to the absolute path of the compiled mpm executable."
  exit 2
fi

case_dir="${SLURM_SUBMIT_DIR:-.}"
cases=(open closed)
case_name="${cases[${SLURM_ARRAY_TASK_ID:-0}]}"
threads="${SLURM_CPUS_PER_TASK:-24}"

echo "Running Siemens et al. (2013) ${case_name} column test"
echo "Executable: ${MPM_BIN}"
echo "Case directory: ${case_dir}"

srun "${MPM_BIN}" \
  -f "${case_dir}/" \
  -i "mpm_${case_name}.json" \
  -p "${threads}"
