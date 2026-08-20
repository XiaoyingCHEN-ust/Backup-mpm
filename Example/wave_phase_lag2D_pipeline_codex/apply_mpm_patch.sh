#!/bin/bash
set -euo pipefail

case_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source_dir=$(cd "${case_dir}/../../mpm" && pwd)
patch_file="${case_dir}/mpm_pipeline_no_flux.patch"
particle_source="${source_dir}/include/particles/particle_threephase_lag.tcc"
solver_source="${source_dir}/include/solvers/thm_mpm_explicit_threephase_lag.tcc"

if [[ ! -s "${patch_file}" ]]; then
  echo "Patch file is missing or empty: ${patch_file}" >&2
  exit 2
fi
if ! grep -Fq "apply_rigid_circle_phase_no_flux" "${patch_file}" ||
   ! grep -Fq "apply_rigid_pipeline_phase_no_flux" "${patch_file}"; then
  echo "Patch file does not contain the required no-flux changes." >&2
  exit 2
fi

if grep -Fq "apply_rigid_circle_phase_no_flux" "${particle_source}" &&
   grep -Fq "apply_rigid_pipeline_phase_no_flux" "${solver_source}"; then
  echo "Moving pipeline no-flux patch is already applied."
elif git -C "${source_dir}" apply --check "${patch_file}"; then
  git -C "${source_dir}" apply "${patch_file}"
  echo "Applied moving pipeline no-flux patch to ${source_dir}."
else
  echo "Patch cannot be applied cleanly to ${source_dir}." >&2
  echo "Do not run the case with an unverified solver source tree." >&2
  exit 2
fi

grep -Fq "apply_rigid_circle_phase_no_flux" "${particle_source}"
grep -Fq "apply_rigid_pipeline_phase_no_flux" "${solver_source}"
