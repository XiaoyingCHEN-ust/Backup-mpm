#!/usr/bin/env bash
# Generic HPC4 batch entry point. Resource details are supplied by the submitter.
#SBATCH --partition=granularmech
#SBATCH --account=comgranmech
#SBATCH --nodes=1
#SBATCH --ntasks=1

set -eo pipefail

case_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
cd "${case_dir}"

if [[ "$#" -eq 0 ]]; then
  printf 'Usage: sbatch run_task.sbatch command [args ...]\n' >&2
  exit 2
fi

exec "$@"
