#!/bin/bash
#SBATCH --job-name=siemens_semi
#SBATCH --partition=granularmech
#SBATCH --account=comgranmech
#SBATCH --array=0-3
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
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
  output_steps uuid input boundary_penalty gravity_scale mesh_variant \
  bounded_pressure_transfer reconstruct_pressure_gradient \
  reconstruct_pressure_force reconstruct_darcy_velocity \
  pressure_projection_rate <<< "${manifest_row}"
bounded_pressure_transfer="${bounded_pressure_transfer%$'\r'}"
reconstruct_pressure_gradient="${reconstruct_pressure_gradient%$'\r'}"
reconstruct_pressure_force="${reconstruct_pressure_force%$'\r'}"
reconstruct_darcy_velocity="${reconstruct_darcy_velocity%$'\r'}"
pressure_projection_rate="${pressure_projection_rate%$'\r'}"
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
  "include/solvers/mpm_base.tcc"
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
grep -Fq "direct_gradient_force" \
  "${mpm_source}/include/solvers/thm_mpm_explicit_threephase_new.tcc"
grep -Fq 'pressure_options["reconstruct_darcy_velocity"]' \
  "${mpm_source}/include/solvers/thm_mpm_explicit_threephase_new.tcc"
grep -Fq 'pressure_options["projection_rate"]' \
  "${mpm_source}/include/solvers/thm_mpm_explicit_threephase_new.tcc"
grep -Fq "time-step-invariant, pressure-only PIC correction" \
  "${mpm_source}/include/particles/particle_threephase_new.tcc"
grep -Fq "Bound the FLIP pressure transfer" \
  "${mpm_source}/include/particles/particle_threephase_new.tcc"
grep -Fq 'liquid_vtk_allowed.emplace_back("force_liquid_pressures");' \
  "${mpm_source}/include/solvers/mpm_base.tcc"
grep -Fq 'liquid_vtk_allowed.emplace_back("force_gas_pressures");' \
  "${mpm_source}/include/solvers/mpm_base.tcc"
grep -Fq '"force_liquid_pressures",' \
  "${mpm_source}/include/solvers/mpm_base.tcc"
grep -Fq '"force_gas_pressures",' \
  "${mpm_source}/include/solvers/mpm_base.tcc"
for source_file in \
  "${mpm_source}/include/particles/particle_threephase_new.tcc" \
  "${mpm_source}/include/solvers/particle_threephase_new.tcc"; do
  grep -Fq "liquid_permeability_ / std::max(liquid_viscosity_" \
    "${source_file}"
  if grep -Fq "0.981 * liquid_viscosity_" "${source_file}"; then
    echo "Obsolete 0.981 mobility factor remains in ${source_file}" >&2
    exit 3
  fi
  if ! grep -Fq \
      "liquid_force = -shapefn_[i] * liquid_pressure_gradient_;" \
      "${source_file}"; then
    echo "Direct pressure-gradient force is absent from ${source_file}" >&2
    exit 3
  fi
  if ! grep -Fq \
      "liquid_density_ * pgravity_ - liquid_pressure_gradient_" \
      "${source_file}"; then
    echo "Darcy-consistent phase velocity is absent from ${source_file}" >&2
    exit 3
  fi
  if ! grep -Fq \
      "time-step-invariant, pressure-only PIC correction" \
      "${source_file}"; then
    echo "Pressure projection-rate correction is absent from ${source_file}" >&2
    exit 3
  fi
done

if [[ ! -f "${input}" ]]; then
  echo "Missing ${input}; run python prepare_semi_implicit_checks.py first" >&2
  exit 4
fi

echo "Running ${label} on $(hostname) with ${threads} threads"
echo "method=${method}, dt=${dt} s, duration=${duration} s, nsteps=${nsteps}"
echo "boundary penalty=${boundary_penalty}"
echo "gravity scale=${gravity_scale:-1}"
echo "mesh variant=${mesh_variant:-coarse}"
echo "bounded pressure transfer=${bounded_pressure_transfer:-False}"
echo "reconstruct pressure gradient=${reconstruct_pressure_gradient:-False}"
echo "reconstruct pressure force=${reconstruct_pressure_force:-False}"
echo "reconstruct Darcy velocity=${reconstruct_darcy_velocity:-False}"
echo "pressure projection rate=${pressure_projection_rate:-0} 1/s"
"${mpm_bin}" -f "${case_dir}/" -i "${input}" -p "${threads}" 2>&1 | \
  awk '/uuid : .*Step:/ {step_count++; if (step_count % 1000 != 0) next} {print}'

result_dir="stability_results/${uuid}"
range_check_args=()
if [[ "${method}" == "semi_implicit" ]]; then
  # particle00000.vtp is written after the first, potentially much larger,
  # semi-implicit step; row 0 already verifies the common initial suction.
  range_check_args+=(--skip-initial-suction-check)
fi
if [[ "${reconstruct_pressure_force}" == "True" ]]; then
  range_check_args+=(--require-force-pressure-fields)
fi
if [[ "${reconstruct_pressure_gradient}" == "True" ]]; then
  range_check_args+=(--require-pressure-gradient-fields)
fi
if [[ "${reconstruct_pressure_gradient}" == "True" && \
      "${reconstruct_pressure_force}" == "True" ]]; then
  if awk -v rate="${pressure_projection_rate:-0}" \
      'BEGIN {exit !(rate > 0)}'; then
    range_check_args+=(
      --allow-transient-dry-gas-velocity-sign-alternation
    )
  else
    range_check_args+=(--reject-dry-gas-velocity-sign-alternation)
  fi
fi
if [[ "${reconstruct_darcy_velocity}" == "True" ]]; then
  range_check_args+=(--reject-dry-liquid-velocity-sign-alternation)
fi
python check_vtp_ranges.py "${result_dir}" "${range_check_args[@]}"
cat > "${result_dir}/SEMI_IMPLICIT_RUN_COMPLETED.txt" <<EOF
case=${case_name}
method=${method}
dt_s=${dt}
duration_s=${duration}
nsteps=${nsteps}
boundary_penalty=${boundary_penalty}
gravity_scale=${gravity_scale:-1}
mesh_variant=${mesh_variant:-coarse}
bounded_pressure_transfer=${bounded_pressure_transfer:-False}
reconstruct_pressure_gradient=${reconstruct_pressure_gradient:-False}
reconstruct_pressure_force=${reconstruct_pressure_force:-False}
reconstruct_darcy_velocity=${reconstruct_darcy_velocity:-False}
pressure_projection_rate=${pressure_projection_rate:-0}
EOF
