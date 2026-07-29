#!/bin/bash
#SBATCH --job-name=siemens_r20
#SBATCH --partition=granularmech
#SBATCH --account=comgranmech
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=02:00:00
#SBATCH --mem=16G
#SBATCH --output=checkpoint-%j.out
#SBATCH --error=checkpoint-%j.err

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
segment="${R20_SEGMENT:?R20_SEGMENT must be a or b}"

case "${segment}" in
  a)
    input="stability_inputs/mpm_r20_checkpoint_open_segment_a.json"
    uuid="siemens2013-r20-checkpoint-open-segment-a"
    final_vtp="particle600000.vtp"
    ;;
  b)
    input="stability_inputs/mpm_r20_checkpoint_open_segment_b.json"
    uuid="siemens2013-r20-checkpoint-open-segment-b"
    final_vtp="particle1200000.vtp"
    checkpoint="stability_results/siemens2013-r20-checkpoint-open-segment-a/particles600000.h5"
    if [[ ! -f "${checkpoint}" ]]; then
      echo "Required segment-A checkpoint is missing: ${checkpoint}" >&2
      exit 4
    fi
    ;;
  *)
    echo "Unknown R20_SEGMENT=${segment}; expected a or b" >&2
    exit 4
    ;;
esac

if [[ ! -f "${input}" ]]; then
  echo "Missing ${input}; run python prepare_checkpoint_validation.py first" >&2
  exit 4
fi
if [[ ! -x "${mpm_bin}" ]]; then
  echo "MPM executable is missing or not executable: ${mpm_bin}" >&2
  exit 126
fi

particle_source="${mpm_source}/include/particles/particle_threephase_new.tcc"
mesh_source="${mpm_source}/include/mesh/mesh.tcc"
solver_source="${mpm_source}/include/solvers/thm_mpm_explicit_threephase_new.tcc"
for source_file in "${particle_source}" "${mesh_source}" "${solver_source}"; do
  if [[ ! -f "${source_file}" ]]; then
    echo "Required MPM source is missing: ${source_file}" >&2
    exit 3
  fi
  if [[ "${source_file}" -nt "${mpm_bin}" ]]; then
    echo "Executable is older than ${source_file}; rebuild before submitting" >&2
    exit 3
  fi
done
grep -q "initial_suction_from_swrc" "${particle_source}"
grep -Fq "std::vector<HDF5Particle> dst_buf(nparticles);" "${mesh_source}"
grep -Fq "for (step_ = start_step; step_ <= nsteps_; ++step_)" "${solver_source}"

echo "Running r20 checkpoint segment ${segment} on $(hostname) with ${threads} threads"
"${mpm_bin}" -f "${case_dir}/" -i "${input}" -p "${threads}" 2>&1 | \
  awk '/uuid : .*Step:/ {step_count++; if (step_count % 10000 != 0) next} {print}'

result_dir="stability_results/${uuid}"
expected_vtp="${result_dir}/${final_vtp}"
if [[ ! -f "${expected_vtp}" ]]; then
  echo "Final segment output is missing: ${expected_vtp}" >&2
  exit 5
fi
if [[ "${segment}" == "a" ]]; then
  python check_vtp_ranges.py "${result_dir}"
else
  python check_vtp_ranges.py --skip-initial-suction-check "${result_dir}"
fi

if [[ "${segment}" == "a" ]]; then
  expected_hdf5="${result_dir}/particles600000.h5"
  if [[ ! -f "${expected_hdf5}" ]]; then
    echo "Segment-A HDF5 checkpoint is missing: ${expected_hdf5}" >&2
    exit 5
  fi
else
  continuous_dir="stability_results/siemens2013-stability-r19b-open-3s-open_2p5e-06"
  python compare_checkpoint_restart.py "${continuous_dir}" "${result_dir}"
fi
