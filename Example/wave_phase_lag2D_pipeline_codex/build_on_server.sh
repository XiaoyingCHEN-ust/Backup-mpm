#!/usr/bin/env bash
# Build the patched source snapshot and run the material/checkpoint preflight.
set -eo pipefail

module load miniconda3/24.3.0-quc3pyu
eval "$(conda shell.bash hook)"
conda activate cbgeo_tbb
set -u

case_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source_dir=${MPM_SOURCE_DIR:-"${case_dir}/../../mpm_hpc_source"}
build_dir=${MPM_BUILD_DIR:-"${source_dir}/build-pipeline"}
jobs=${BUILD_JOBS:-${SLURM_CPUS_PER_TASK:-16}}

cmake -S "${source_dir}" -B "${build_dir}" \
  -DCMAKE_BUILD_TYPE=Release \
  -DMPM_BUILD_TESTING=ON \
  -DMPM_BUILD_LIB=OFF
cmake --build "${build_dir}" --parallel "${jobs}"

registered_tests=$(ctest --test-dir "${build_dir}" -N)
for required_test in \
  sanisand_material_test \
  mohr_coulomb_handoff_test \
  hdf5_particle_test \
  prescribed_pressure_validation_test \
  rigid_pipeline2d_test; do
  if ! grep -q "${required_test}" <<< "${registered_tests}"; then
    printf 'Required CTest target was not registered: %s\n' \
      "${required_test}" >&2
    exit 3
  fi
done

ctest --test-dir "${build_dir}" --output-on-failure \
  -R '^(sanisand_material_test|mohr_coulomb_handoff_test|hdf5_particle_test|prescribed_pressure_validation_test|rigid_pipeline2d_test)$'

"${build_dir}/sanisand_low_pressure_driver" \
  > "${case_dir}/material_preflight_3kPa.csv"

printf 'Validated executable: %s\n' "${build_dir}/mpm"
printf 'Inspect low-confinement diagnostic: %s\n' \
  "${case_dir}/material_preflight_3kPa.csv"
printf 'Export before the run if using another location:\n'
printf '  export MPM_BIN=%s\n' "${build_dir}/mpm"
