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
  "particles/threephase_pressure_state.h"
  "particles/particle_threephase_new.h"
  "particles/particle_threephase_new.tcc"
  "solvers/particle_threephase_new.h"
  "solvers/particle_threephase_new.tcc"
  "solvers/thm_mpm_explicit_threephase_new.h"
  "solvers/thm_mpm_explicit_threephase_new.tcc"
)

if [[ ! -d "${mpm_source}/build" ]]; then
  echo "MPM build directory is missing: ${mpm_source}/build" >&2
  exit 2
fi

for relative_path in "${source_files[@]}"; do
  replacement="${replacement_root}/${relative_path}"
  target="${mpm_source}/include/${relative_path}"

  if [[ ! -f "${replacement}" ]]; then
    echo "Missing replacement for ${relative_path}: ${replacement}" >&2
    exit 2
  fi
  if [[ ! -f "${target}" ]]; then
    echo "Restoring missing source file: ${target}"
  fi
  install -D -m 0644 "${replacement}" "${target}"
done

checkpoint_patcher="${script_dir}/patch_checkpoint_resume.py"
if [[ ! -f "${checkpoint_patcher}" ]]; then
  echo "Checkpoint patcher is missing: ${checkpoint_patcher}" >&2
  exit 2
fi
python "${checkpoint_patcher}" --mpm-source "${mpm_source}"

grep -q "initial_suction_from_swrc" \
  "${mpm_source}/include/particles/particle_threephase_new.tcc"
grep -q "initial_suction_from_swrc" \
  "${mpm_source}/include/solvers/particle_threephase_new.tcc"
grep -q "ThreePhasePressureState" \
  "${mpm_source}/include/particles/threephase_pressure_state.h"
grep -q "ThreePhasePressureState" \
  "${mpm_source}/include/particles/particle_threephase_new.h"
grep -q "ThreePhasePressureState" \
  "${mpm_source}/include/solvers/particle_threephase_new.h"
grep -q "solve_semi_implicit_pressure" \
  "${mpm_source}/include/solvers/thm_mpm_explicit_threephase_new.h"
grep -q 'pressure_integration == "semi_implicit"' \
  "${mpm_source}/include/solvers/thm_mpm_explicit_threephase_new.tcc"
grep -Fq "No supported pressure degrees of freedom" \
  "${mpm_source}/include/solvers/thm_mpm_explicit_threephase_new.tcc"
grep -Fq "compact_liquid_pressure_increment" \
  "${mpm_source}/include/solvers/thm_mpm_explicit_threephase_new.tcc"
grep -Fq "Row-sum lump the pressure-storage matrix" \
  "${mpm_source}/include/solvers/thm_mpm_explicit_threephase_new.tcc"
particle_implementations=(
  "particles/particle_threephase_new.tcc"
  "solvers/particle_threephase_new.tcc"
)
for relative_path in "${particle_implementations[@]}"; do
  source_file="${mpm_source}/include/${relative_path}"
  if ! grep -Fq \
      "Fixed-gas formulation: solve only the liquid mass balance." \
      "${source_file}"; then
    echo "Corrected fixed-gas formulation is absent from ${source_file}" >&2
    exit 3
  fi
  if grep -Fq \
      "PIC_liquid_pressure_ - (-this->liquid_density_ * 9.81" \
      "${source_file}"; then
    echo "Obsolete double-gravity liquid-force expression remains in ${source_file}" >&2
    exit 3
  fi
  if ! grep -Fq \
      "liquid_permeability_ / std::max(liquid_viscosity_" \
      "${source_file}"; then
    echo "Correct Darcy mobility is absent from ${source_file}" >&2
    exit 3
  fi
  if grep -Fq "0.981 * liquid_viscosity_" "${source_file}"; then
    echo "Obsolete 0.981 mobility factor remains in ${source_file}" >&2
    exit 3
  fi
done
grep -Fq "std::vector<HDF5Particle> dst_buf(nparticles);" \
  "${mpm_source}/include/mesh/mesh.tcc"
grep -Fq "for (step_ = start_step; step_ <= nsteps_; ++step_)" \
  "${mpm_source}/include/solvers/thm_mpm_explicit_threephase_new.tcc"

cmake --build "${mpm_source}/build" --parallel "${threads}"

mpm_bin="${mpm_source}/build/mpm"
if [[ ! -x "${mpm_bin}" ]]; then
  echo "Build completed without an executable at ${mpm_bin}" >&2
  exit 2
fi

# The replacements are versioned in Backup-mpm, so historical side-by-side
# copies only create ambiguity about which source was compiled.  Remove the
# two legacy backup naming schemes after a successful rebuild.
for relative_path in \
  "particles/particle_threephase_new.h" \
  "particles/particle_threephase_new.tcc" \
  "solvers/particle_threephase_new.h" \
  "solvers/particle_threephase_new.tcc"; do
  target="${mpm_source}/include/${relative_path}"
  for backup in "${target}.pre-siemens-validation" "${target}".before-*; do
    if [[ -f "${backup}" ]]; then
      rm -- "${backup}"
      echo "Removed obsolete source backup: ${backup}"
    fi
  done
done

echo "Installed Siemens validation patch and rebuilt ${mpm_bin}"
