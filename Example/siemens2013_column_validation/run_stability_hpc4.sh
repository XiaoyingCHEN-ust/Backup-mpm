#!/bin/bash
#SBATCH --job-name=siemens_dt
#SBATCH --partition=granularmech
#SBATCH --account=comgranmech
#SBATCH --array=0-5
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --gres=gpu:1
#SBATCH --time=01:00:00
#SBATCH --mem=16G
#SBATCH --output=stability-%A_%a.out
#SBATCH --error=stability-%A_%a.err

set -eo pipefail

case_dir="${SLURM_SUBMIT_DIR:?SLURM_SUBMIT_DIR is not set}"
cd "${case_dir}"

module load miniconda3/24.3.0-quc3pyu
eval "$("$(command -v conda)" shell.bash hook)"
conda activate cbgeo_tbb
set -u

mpm_source="${MPM_SOURCE:-/home/xchenjm/chen/MPM/mpm}"
mpm_bin="${MPM_BIN:-${mpm_source}/build/mpm}"
threads="${SLURM_CPUS_PER_TASK:-32}"
task_index="${SLURM_ARRAY_TASK_ID:-0}"
manifest="stability_inputs/manifest.csv"
if [[ ! -f "${manifest}" ]]; then
  echo "Missing ${manifest}; run python prepare_stability_checks.py first" >&2
  exit 4
fi
manifest_row="$(awk -F, -v target="$((task_index + 2))" \
  'NR == target {print; exit}' "${manifest}")"
if [[ -z "${manifest_row}" ]]; then
  echo "No manifest entry for array index ${task_index}" >&2
  exit 4
fi
IFS=',' read -r manifest_index label case_name dt duration nsteps \
  output_steps uuid input <<< "${manifest_row}"
input="${input%$'\r'}"
if [[ "${manifest_index}" != "${task_index}" ]]; then
  echo "Manifest index mismatch: expected ${task_index}, got ${manifest_index}" >&2
  exit 4
fi

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
echo "dt=${dt} s, duration=${duration} s, nsteps=${nsteps}"
"${mpm_bin}" -f "${case_dir}/" -i "${input}" -p "${threads}" 2>&1 | \
  awk '/uuid : .*Step:/ {step_count++; if (step_count % 10000 != 0) next} {print}'
python check_vtp_ranges.py "stability_results/${uuid}"
