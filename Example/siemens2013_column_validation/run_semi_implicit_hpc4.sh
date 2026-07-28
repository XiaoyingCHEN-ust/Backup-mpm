#!/bin/bash
#SBATCH --job-name=siemens_semi
#SBATCH --partition=granularmech
#SBATCH --account=comgranmech
#SBATCH --array=0-3
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --gres=gpu:1
#SBATCH --time=00:30:00
#SBATCH --mem=16G
#SBATCH --output=semi-%A_%a.out
#SBATCH --error=semi-%A_%a.err

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
manifest="semi_implicit_inputs/manifest.csv"

if [[ ! -f "${manifest}" ]]; then
  echo "Missing ${manifest}; run python prepare_semi_implicit_checks.py first" >&2
  exit 4
fi
manifest_row="$(awk -F, -v target="$((task_index + 2))" \
  'NR == target {print; exit}' "${manifest}")"
if [[ -z "${manifest_row}" ]]; then
  echo "No manifest entry for array index ${task_index}" >&2
  exit 4
fi
IFS=',' read -r manifest_index label case_name method dt duration nsteps \
  output_steps uuid input boundary_penalty <<< "${manifest_row}"
boundary_penalty="${boundary_penalty%$'\r'}"
if [[ "${manifest_index}" != "${task_index}" ]]; then
  echo "Manifest index mismatch: expected ${task_index}, got ${manifest_index}" >&2
  exit 4
fi

if [[ ! -x "${mpm_bin}" ]]; then
  echo "MPM executable is missing or not executable: ${mpm_bin}" >&2
  exit 126
fi
source_files=(
  "include/particles/threephase_pressure_state.h"
  "include/particles/particle_threephase_new.h"
  "include/particles/particle_threephase_new.tcc"
  "include/solvers/particle_threephase_new.h"
  "include/solvers/particle_threephase_new.tcc"
  "include/solvers/thm_mpm_explicit_threephase_new.h"
  "include/solvers/thm_mpm_explicit_threephase_new.tcc"
)
for relative_path in "${source_files[@]}"; do
  source_file="${mpm_source}/${relative_path}"
  if [[ ! -f "${source_file}" ]]; then
    echo "Required source file is missing: ${source_file}" >&2
    exit 3
  fi
  if [[ "${source_file}" -nt "${mpm_bin}" ]]; then
    echo "Executable is older than ${source_file}; rebuild before submitting" >&2
    exit 3
  fi
done
grep -Fq "ThreePhasePressureState" \
  "${mpm_source}/include/particles/particle_threephase_new.h"
grep -Fq "solve_semi_implicit_pressure" \
  "${mpm_source}/include/solvers/thm_mpm_explicit_threephase_new.tcc"
grep -Fq "compact_liquid_pressure_increment" \
  "${mpm_source}/include/solvers/thm_mpm_explicit_threephase_new.tcc"

if [[ ! -f "${input}" ]]; then
  echo "Missing ${input}; run python prepare_semi_implicit_checks.py first" >&2
  exit 4
fi

echo "Running ${label} on $(hostname) with ${threads} threads"
echo "method=${method}, dt=${dt} s, duration=${duration} s, nsteps=${nsteps}"
echo "boundary penalty=${boundary_penalty}"
"${mpm_bin}" -f "${case_dir}/" -i "${input}" -p "${threads}" 2>&1 | \
  awk '/uuid : .*Step:/ {step_count++; if (step_count % 1000 != 0) next} {print}'

result_dir="stability_results/${uuid}"
if [[ "${method}" == "semi_implicit" ]]; then
  # particle00000.vtp is written after the first, potentially much larger,
  # semi-implicit step; row 0 already verifies the common initial suction.
  python check_vtp_ranges.py "${result_dir}" --skip-initial-suction-check
else
  python check_vtp_ranges.py "${result_dir}"
fi
cat > "${result_dir}/SEMI_IMPLICIT_RUN_COMPLETED.txt" <<EOF
case=${case_name}
method=${method}
dt_s=${dt}
duration_s=${duration}
nsteps=${nsteps}
boundary_penalty=${boundary_penalty}
EOF
