#!/usr/bin/env bash
# Execute one pipeline case on a local Linux workstation with the same
# validation and completion-audit contract used by the HPC wrapper.
set -euo pipefail

case_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)
cd "${case_dir}"

validation_args=()
if [[ "${1:-}" == "--smoke" ]]; then
  validation_args=(--local-smoke)
  shift
fi
if [[ $# -ne 1 ]]; then
  printf 'Usage: run_case_local.sh [--smoke] configs/group/case.json\n' >&2
  exit 2
fi
config=$1
config_abs=$(realpath -e -- "${config}")
case "${config_abs}" in
  "${case_dir}"/*) config_relative=${config_abs#"${case_dir}"/} ;;
  *)
    printf 'Configuration must be inside the pipeline case directory: %s\n' \
      "${config_abs}" >&2
    exit 2
    ;;
esac

threads=${MPM_THREADS:-8}
python_bin=${PYTHON_BIN:-python3}
mpm_bin=${MPM_BIN:-"${case_dir}/../../mpm_hpc_source/build-pipeline/mpm"}

if ! [[ "${threads}" =~ ^[1-9][0-9]*$ ]]; then
  printf 'MPM_THREADS must be a positive integer, got: %s\n' "${threads}" >&2
  exit 2
fi
if [[ ! -x "${mpm_bin}" ]]; then
  printf 'MPM executable is missing or not executable: %s\n' "${mpm_bin}" >&2
  exit 2
fi
if ! command -v "${python_bin}" >/dev/null 2>&1; then
  printf 'Python executable is unavailable: %s\n' "${python_bin}" >&2
  exit 2
fi
if ! command -v flock >/dev/null 2>&1; then
  printf 'flock is required to prevent concurrent writes to one result target\n' >&2
  exit 2
fi
if ! command -v sha256sum >/dev/null 2>&1; then
  printf 'sha256sum is required to derive the result-directory lock\n' >&2
  exit 2
fi

lock_root="${case_dir}/results/.local_run_locks"
mkdir -p -- "${lock_root}"
result_target=$("${python_bin}" validate_case.py "${config_abs}" \
  --print-result-directory "${validation_args[@]}")
lock_id=$(printf '%s' "${result_target}" | sha256sum | cut -d' ' -f1)
exec 9>"${lock_root}/${lock_id}.lock"
if ! flock -n 9; then
  printf 'This local case is already running: %s\n' "${config_relative}" >&2
  exit 3
fi

"${python_bin}" validate_case.py "${config_abs}" \
  --runtime --prepare-output --clear-completion "${validation_args[@]}"
export OMP_NUM_THREADS="${threads}"
export TBB_NUM_THREADS="${threads}"
"${mpm_bin}" -p "${threads}" -f ./ -i "${config_relative}"
"${python_bin}" validate_case.py "${config_abs}" \
  --write-completion "${validation_args[@]}"
