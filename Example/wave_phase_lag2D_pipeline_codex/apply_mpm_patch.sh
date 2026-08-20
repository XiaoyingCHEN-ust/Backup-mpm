#!/bin/bash
set -euo pipefail

case_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source_dir=$(cd "${case_dir}/../../mpm" && pwd)
patch_file="${case_dir}/mpm_pipeline_no_flux.patch"
pressure_patch_file="${case_dir}/mpm_initial_hydrostatic_pressure.patch"
surface_patch_file="${case_dir}/mpm_submerged_surface_traction.patch"
particle_header="${source_dir}/include/particles/particle_threephase_lag.h"
particle_source="${source_dir}/include/particles/particle_threephase_lag.tcc"
solver_source="${source_dir}/include/solvers/thm_mpm_explicit_threephase_lag.tcc"

# `mpm` can either be its own Git worktree or a subdirectory of a larger
# worktree.  Run git-apply from the actual worktree root and prepend the source
# directory when needed.  Running `git -C <subdir> apply` can silently skip all
# patch paths while still returning success.  The selected upstream files also
# contain mixed LF/CRLF regions, so context matching must ignore line-ending
# whitespace while the post-apply marker checks enforce the intended result.
repo_root=$(git -C "${source_dir}" rev-parse --show-toplevel)
repo_root=$(cd "${repo_root}" && pwd)
if [[ "${source_dir}" == "${repo_root}" ]]; then
  apply_directory_args=()
elif [[ "${source_dir}" == "${repo_root}/"* ]]; then
  source_prefix=${source_dir#"${repo_root}/"}
  apply_directory_args=(--directory="${source_prefix}")
else
  echo "MPM source is outside its Git worktree: ${source_dir}" >&2
  exit 2
fi

if [[ ! -s "${patch_file}" ]]; then
  echo "Patch file is missing or empty: ${patch_file}" >&2
  exit 2
fi
if [[ ! -s "${pressure_patch_file}" ]]; then
  echo "Patch file is missing or empty: ${pressure_patch_file}" >&2
  exit 2
fi
if [[ ! -s "${surface_patch_file}" ]]; then
  echo "Patch file is missing or empty: ${surface_patch_file}" >&2
  exit 2
fi
if ! grep -Fq "apply_rigid_circle_phase_no_flux" "${patch_file}" ||
   ! grep -Fq "apply_rigid_pipeline_phase_no_flux" "${patch_file}"; then
  echo "Patch file does not contain the required no-flux changes." >&2
  exit 2
fi
if ! grep -Fq "has_input_initial_liquid_pressure_" "${pressure_patch_file}" ||
   ! grep -Fq "registered capillary suction and target saturation" "${pressure_patch_file}"; then
  echo "Patch file does not contain the required phase-pressure fixes." >&2
  exit 2
fi
if ! grep -Fq "apply_submerged_surface_pressure_traction" "${surface_patch_file}" ||
   ! grep -Fq "submerged_surface_pressure_traction_" "${surface_patch_file}"; then
  echo "Patch file does not contain the required surface-traction fixes." >&2
  exit 2
fi

if grep -Fq "apply_rigid_circle_phase_no_flux" "${particle_source}" &&
   grep -Fq "apply_rigid_pipeline_phase_no_flux" "${solver_source}"; then
  echo "Moving pipeline no-flux patch is already applied."
elif git -C "${repo_root}" apply "${apply_directory_args[@]}" \
    --ignore-space-change --ignore-whitespace \
    --check --verbose "${patch_file}"; then
  git -C "${repo_root}" apply "${apply_directory_args[@]}" \
    --ignore-space-change --ignore-whitespace \
    --verbose "${patch_file}"
else
  echo "Patch cannot be applied cleanly to ${source_dir}." >&2
  echo "Do not run the case with an unverified solver source tree." >&2
  exit 2
fi

if ! grep -Fq "apply_rigid_circle_phase_no_flux" "${particle_source}" ||
   ! grep -Fq "apply_rigid_pipeline_phase_no_flux" "${solver_source}"; then
  echo "Patch command completed without installing the required markers." >&2
  exit 3
fi

if grep -Fq "has_input_initial_liquid_pressure_" "${particle_header}" &&
   grep -Fq "input_initial_liquid_pressure_" "${particle_source}" &&
   grep -Fq "registered capillary suction and target saturation" "${particle_source}"; then
  echo "Hydrostatic initial-pressure patch is already applied."
elif git -C "${repo_root}" apply "${apply_directory_args[@]}" \
    --unidiff-zero \
    --ignore-space-change --ignore-whitespace \
    --check --verbose "${pressure_patch_file}"; then
  git -C "${repo_root}" apply "${apply_directory_args[@]}" \
    --unidiff-zero \
    --ignore-space-change --ignore-whitespace \
    --verbose "${pressure_patch_file}"
else
  echo "Hydrostatic pressure patch cannot be applied cleanly to ${source_dir}." >&2
  echo "Do not run the case with an unverified initial phase-pressure state." >&2
  exit 2
fi

if ! grep -Fq "has_input_initial_liquid_pressure_" "${particle_header}" ||
   ! grep -Fq "input_initial_liquid_pressure_" "${particle_source}" ||
   ! grep -Fq "registered capillary suction and target saturation" "${particle_source}"; then
  echo "Patch command completed without installing the pressure markers." >&2
  exit 3
fi

if grep -Fq "apply_submerged_surface_pressure_traction" "${particle_source}" &&
   grep -Fq "submerged_surface_pressure_traction_" "${solver_source}"; then
  echo "Matched submerged-surface traction patch is already applied."
elif git -C "${repo_root}" apply "${apply_directory_args[@]}" \
    --unidiff-zero \
    --ignore-space-change --ignore-whitespace \
    --check --verbose "${surface_patch_file}"; then
  git -C "${repo_root}" apply "${apply_directory_args[@]}" \
    --unidiff-zero \
    --ignore-space-change --ignore-whitespace \
    --verbose "${surface_patch_file}"
else
  echo "Surface-traction patch cannot be applied cleanly to ${source_dir}." >&2
  echo "Do not run the case with an unbalanced submerged boundary." >&2
  exit 2
fi

if ! grep -Fq "apply_submerged_surface_pressure_traction" "${particle_source}" ||
   ! grep -Fq "submerged_surface_pressure_traction_" "${solver_source}"; then
  echo "Patch command completed without installing the surface-traction markers." >&2
  exit 3
fi

echo "Verified moving no-flux, hydrostatic-pressure, and matched surface-traction patches in ${source_dir}."
