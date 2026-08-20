#!/bin/bash
set -euo pipefail

case_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source_dir=$(cd "${case_dir}/../../mpm" && pwd)
patch_file="${case_dir}/mpm_pipeline_no_flux.patch"
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
if ! grep -Fq "apply_rigid_circle_phase_no_flux" "${patch_file}" ||
   ! grep -Fq "apply_rigid_pipeline_phase_no_flux" "${patch_file}"; then
  echo "Patch file does not contain the required no-flux changes." >&2
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

echo "Verified moving pipeline no-flux patch in ${source_dir}."
