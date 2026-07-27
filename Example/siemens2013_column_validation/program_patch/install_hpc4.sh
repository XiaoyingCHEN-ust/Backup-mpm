#!/bin/bash

# Install the validated three-phase particle replacement in the source tree
# used by the HPC4 executable, then rebuild that executable.

set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "${script_dir}/../../.." && pwd)"
mpm_source="${MPM_SOURCE:-/home/xchenjm/chen/MPM/mpm}"
threads="${SLURM_CPUS_PER_TASK:-32}"

replacement_root="${repo_root}/server_overrides/mpm/include"
source_files=(
  "particles/particle_threephase_new.tcc"
  "solvers/particle_threephase_new.tcc"
)

if [[ ! -d "${mpm_source}/build" ]]; then
  echo "MPM build directory is missing: ${mpm_source}/build" >&2
  exit 2
fi

for relative_path in "${source_files[@]}"; do
  replacement="${replacement_root}/${relative_path}"
  target="${mpm_source}/include/${relative_path}"
  backup="${target}.pre-siemens-validation"

  if [[ ! -f "${replacement}" ]]; then
    echo "Missing replacement for ${relative_path}: ${replacement}" >&2
    exit 2
  fi
  if [[ -f "${target}" && ! -f "${backup}" ]]; then
    cp --preserve=mode,timestamps "${target}" "${backup}"
  elif [[ ! -f "${target}" ]]; then
    echo "Restoring missing source file: ${target}"
  fi
  install -D -m 0644 "${replacement}" "${target}"
done

grep -q "initial_suction_from_swrc" \
  "${mpm_source}/include/particles/particle_threephase_new.tcc"
grep -q "initial_suction_from_swrc" \
  "${mpm_source}/include/solvers/particle_threephase_new.tcc"
for relative_path in "${source_files[@]}"; do
  source_file="${mpm_source}/include/${relative_path}"
  if grep -Fq \
      "PIC_liquid_pressure_ - (-this->liquid_density_ * 9.81" \
      "${source_file}"; then
    echo "Obsolete double-gravity liquid-force expression remains in ${source_file}" >&2
    exit 3
  fi
done

cmake --build "${mpm_source}/build" --parallel "${threads}"

mpm_bin="${mpm_source}/build/mpm"
if [[ ! -x "${mpm_bin}" ]]; then
  echo "Build completed without an executable at ${mpm_bin}" >&2
  exit 2
fi

echo "Installed Siemens validation patch and rebuilt ${mpm_bin}"
