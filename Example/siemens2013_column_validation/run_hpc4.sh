#!/bin/bash
#SBATCH --job-name=siemens3p
#SBATCH --partition=granularmech
#SBATCH --account=comgranmech
#SBATCH --array=0-1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=24:00:00
#SBATCH --mem=16G
#SBATCH --output=slurm-%A_%a.out
#SBATCH --error=slurm-%A_%a.err

set -eo pipefail

if [[ "${ALLOW_UNSEGMENTED_FULL_RUN:-0}" != "1" ]]; then
  echo "The direct full-duration run is disabled because it exceeds the 24-hour limit." >&2
  echo "Use: python submit_production_hpc4.py" >&2
  exit 4
fi

case_dir="${SLURM_SUBMIT_DIR:?SLURM_SUBMIT_DIR is not set}"
cd "${case_dir}"

module load miniconda3/24.3.0-quc3pyu
eval "$("$(command -v conda)" shell.bash hook)"
conda activate cbgeo_tbb

mpm_source="${MPM_SOURCE:-/home/xchenjm/chen/MPM/mpm}"
mpm_bin="${MPM_BIN:-${mpm_source}/build/mpm}"
cases=(open closed)
case_name="${cases[${SLURM_ARRAY_TASK_ID:-0}]}"
threads="${SLURM_CPUS_PER_TASK:-32}"

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
  if grep -Fq \
      "PIC_liquid_pressure_ - (-this->liquid_density_ * 9.81" \
      "${source_file}"; then
    echo "Obsolete double-gravity liquid-force expression remains in ${source_file}" >&2
    exit 3
  fi
  if [[ "${source_file}" -nt "${mpm_bin}" ]]; then
    echo "Executable is older than ${source_file}; rebuild before submitting" >&2
    exit 3
  fi
done

echo "Running Siemens et al. (2013) ${case_name} column test"
echo "Node: $(hostname)"
echo "Executable: ${mpm_bin}"
echo "Case directory: ${case_dir}"
echo "Threads: ${threads}"
"${mpm_bin}" \
  -f "${case_dir}/" \
  -i "mpm_${case_name}.json" \
  -p "${threads}" 2>&1 | \
  awk '/uuid : .*Step:/ {step_count++; if (step_count % 10000 != 0) next} {print}'

python check_vtp_ranges.py "results/siemens2013-${case_name}"
