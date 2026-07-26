#!/bin/bash
#SBATCH --job-name=siemens_dt
#SBATCH --partition=granularmech
#SBATCH --account=comgranmech
#SBATCH --array=0-7
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --gres=gpu:1
#SBATCH --time=04:00:00
#SBATCH --mem=16G
#SBATCH --output=stability-%A_%a.out
#SBATCH --error=stability-%A_%a.err

set -euo pipefail

case_dir="${SLURM_SUBMIT_DIR:?SLURM_SUBMIT_DIR is not set}"
cd "${case_dir}"

module load miniconda3/24.3.0-quc3pyu
eval "$("$(command -v conda)" shell.bash hook)"
conda activate cbgeo_tbb

mpm_source="${MPM_SOURCE:-/home/xchenjm/chen/MPM/mpm}"
mpm_bin="${MPM_BIN:-${mpm_source}/build/mpm}"
threads="${SLURM_CPUS_PER_TASK:-32}"
labels=(
  open_5e-04 open_5e-05 open_5e-06 open_1e-06
  closed_5e-04 closed_5e-05 closed_5e-06 closed_1e-06
)
label="${labels[${SLURM_ARRAY_TASK_ID:-0}]}"
input="stability_inputs/mpm_${label}.json"
uuid="siemens2013-stability-${label}"

if [[ ! -x "${mpm_bin}" ]]; then
  echo "MPM executable is missing or not executable: ${mpm_bin}" >&2
  exit 126
fi
for relative_path in \
  include/particles/particle_threephase_new.tcc \
  include/solvers/particle_threephase_new.tcc; do
  source_file="${mpm_source}/${relative_path}"
  if ! grep -q "initial_suction_from_swrc" "${source_file}"; then
    echo "Validation patch is absent from ${source_file}" >&2
    exit 3
  fi
  if [[ "${source_file}" -nt "${mpm_bin}" ]]; then
    echo "Executable is older than ${source_file}; rebuild before submitting" >&2
    exit 3
  fi
done
if [[ ! -f "${input}" ]]; then
  echo "Missing ${input}; run python prepare_stability_checks.py first" >&2
  exit 4
fi

echo "Running ${label} on $(hostname) with ${threads} threads"
"${mpm_bin}" -f "${case_dir}/" -i "${input}" -p "${threads}" 2>&1 | \
  awk '/uuid : .*Step:/ {step_count++; if (step_count % 10000 != 0) next} {print}'
python check_vtp_ranges.py "stability_results/${uuid}"
