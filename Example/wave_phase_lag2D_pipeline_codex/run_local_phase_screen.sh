#!/usr/bin/env bash
# Run only the registered baseline screen EQ/SANISAND/phase-replay DAG locally.
set -euo pipefail

case_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)
cd "${case_dir}"

if [[ $# -ne 0 ]]; then
  printf 'Usage: run_local_phase_screen.sh\n' >&2
  exit 2
fi
if ! command -v flock >/dev/null 2>&1; then
  printf 'flock is required to serialize the local phase-screen controller\n' >&2
  exit 2
fi

python_bin=${PYTHON_BIN:-python3}
lock_root="${case_dir}/results/.local_run_locks"
mkdir -p -- "${lock_root}"
exec 9>"${lock_root}/formal-phase-screen.lock"
if ! flock -n 9; then
  printf 'The formal local phase-screen controller is already running\n' >&2
  exit 3
fi

exec "${python_bin}" "${case_dir}/local_phase_screen.py"
