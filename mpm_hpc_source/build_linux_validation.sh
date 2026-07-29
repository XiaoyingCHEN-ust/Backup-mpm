#!/usr/bin/env bash

# Configure and build the MPM executable used by the Siemens column validation.
# Dependency locations may be supplied through the normal CMake environment or
# as additional arguments to this script.

set -euo pipefail

source_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
build_dir="${MPM_BUILD_DIR:-${source_dir}/build}"
build_jobs="${MPM_BUILD_JOBS:-$(nproc)}"

cmake -S "${source_dir}" -B "${build_dir}" \
  -DCMAKE_BUILD_TYPE=Release \
  -DMPM_BUILD_TESTING=OFF \
  -DMPM_BUILD_LIB=OFF \
  "$@"

cmake --build "${build_dir}" --parallel "${build_jobs}"

if [[ ! -x "${build_dir}/mpm" ]]; then
  echo "Build completed without an executable at ${build_dir}/mpm" >&2
  exit 2
fi

echo "Built Siemens-validation executable: ${build_dir}/mpm"
