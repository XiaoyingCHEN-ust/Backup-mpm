#!/bin/bash
#SBATCH --job-name=siemens_prod
#SBATCH --partition=granularmech
#SBATCH --account=comgranmech
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --gres=gpu:1
#SBATCH --time=24:00:00
#SBATCH --mem=16G
#SBATCH --output=production-%j.out
#SBATCH --error=production-%j.err

set -eo pipefail

case_dir="${SLURM_SUBMIT_DIR:?SLURM_SUBMIT_DIR is not set}"
cd "${case_dir}"

module load miniconda3/24.3.0-quc3pyu
eval "$("$(command -v conda)" shell.bash hook)"
conda activate cbgeo_tbb
set -u

row_target="${PRODUCTION_ROW:?PRODUCTION_ROW is not set}"
manifest="production_inputs/manifest.csv"
if [[ ! -f "${manifest}" ]]; then
  echo "Missing ${manifest}; run python prepare_production_segments.py first" >&2
  exit 4
fi
manifest_row="$(awk -F, -v target="$((row_target + 2))" \
  'NR == target {print; exit}' "${manifest}")"
if [[ -z "${manifest_row}" ]]; then
  echo "No production manifest entry for row ${row_target}" >&2
  exit 4
fi
IFS=',' read -r row_index case_name segment_index segment_count start_step \
  end_step total_steps dt nominal_start nominal_end resume_time uuid input \
  result_dir checkpoint final_particle final_checkpoint <<< "${manifest_row}"
final_checkpoint="${final_checkpoint%$'\r'}"
if [[ "${row_index}" != "${row_target}" ]]; then
  echo "Manifest row mismatch: expected ${row_target}, got ${row_index}" >&2
  exit 4
fi
if [[ ! -f "${input}" ]]; then
  echo "Production input is missing: ${input}" >&2
  exit 4
fi
if [[ -e "${result_dir}" ]]; then
  echo "Refusing to mix output with existing directory: ${result_dir}" >&2
  exit 4
fi
if [[ "${checkpoint}" != "none" ]]; then
  if [[ ! -f "${checkpoint}" ]]; then
    echo "Previous checkpoint is missing: ${checkpoint}" >&2
    exit 4
  fi
  previous_result_dir="$(dirname "${checkpoint}")"
  if [[ ! -f "${previous_result_dir}/SEGMENT_COMPLETED.txt" ]]; then
    echo "Previous segment lacks its completion marker: ${previous_result_dir}" >&2
    exit 4
  fi
fi

mpm_source="${MPM_SOURCE:-/home/xchenjm/chen/MPM/mpm}"
mpm_bin="${MPM_BIN:-${mpm_source}/build/mpm}"
threads="${SLURM_CPUS_PER_TASK:-32}"
if [[ ! -x "${mpm_bin}" ]]; then
  echo "MPM executable is missing or not executable: ${mpm_bin}" >&2
  exit 126
fi

particle_source="${mpm_source}/include/particles/particle_threephase_new.tcc"
solver_particle_source="${mpm_source}/include/solvers/particle_threephase_new.tcc"
mesh_source="${mpm_source}/include/mesh/mesh.tcc"
solver_source="${mpm_source}/include/solvers/thm_mpm_explicit_threephase_new.tcc"
for source_file in "${particle_source}" "${solver_particle_source}" \
  "${mesh_source}" "${solver_source}"; do
  if [[ ! -f "${source_file}" ]]; then
    echo "Required MPM source is missing: ${source_file}" >&2
    exit 3
  fi
  if [[ "${source_file}" -nt "${mpm_bin}" ]]; then
    echo "Executable is older than ${source_file}; rebuild before submitting" >&2
    exit 3
  fi
done
for source_file in "${particle_source}" "${solver_particle_source}"; do
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
done
grep -Fq "std::vector<HDF5Particle> dst_buf(nparticles);" "${mesh_source}"
grep -Fq "for (step_ = start_step; step_ <= nsteps_; ++step_)" "${solver_source}"

echo "Running ${case_name} production segment $((segment_index + 1))/${segment_count}"
echo "Global steps ${start_step}-${end_step}; nominal time ${nominal_start}-${nominal_end} s"
echo "Node $(hostname); threads ${threads}; executable ${mpm_bin}"
"${mpm_bin}" -f "${case_dir}/" -i "${input}" -p "${threads}" 2>&1 | \
  awk '/uuid : .*Step:/ {step_count++; if (step_count % 10000 != 0) next} {print}'

if [[ ! -f "${final_particle}" ]]; then
  echo "Final particle output is missing: ${final_particle}" >&2
  exit 5
fi
if [[ ! -f "${final_checkpoint}" ]]; then
  echo "Final HDF5 checkpoint is missing: ${final_checkpoint}" >&2
  exit 5
fi
if [[ "${segment_index}" == "0" ]]; then
  python check_vtp_ranges.py "${result_dir}"
else
  python check_vtp_ranges.py --skip-initial-suction-check "${result_dir}"
fi

{
  echo "case=${case_name}"
  echo "segment_index=${segment_index}"
  echo "segment_count=${segment_count}"
  echo "start_step=${start_step}"
  echo "end_step=${end_step}"
  echo "nominal_end_time_s=${nominal_end}"
  echo "uuid=${uuid}"
} > "${result_dir}/SEGMENT_COMPLETED.txt"
echo "Completed ${uuid}"
