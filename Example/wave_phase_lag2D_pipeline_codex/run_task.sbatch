#!/usr/bin/env bash
# Generic HPC4 batch entry point. Resource details are supplied by the submitter.
#SBATCH --partition=granularmech
#SBATCH --account=comgranmech
#SBATCH --nodes=1
#SBATCH --ntasks=1

set -eo pipefail

# Slurm copies the submitted script into its spool before execution, so
# BASH_SOURCE points at the spool copy rather than this case directory.  sbatch
# exports the directory from which it was invoked as SLURM_SUBMIT_DIR; the
# Python submitter always invokes sbatch with this case directory as cwd.
if [[ -z "${SLURM_SUBMIT_DIR:-}" ]]; then
  printf 'SLURM_SUBMIT_DIR is unavailable; submit this wrapper through sbatch.\n' >&2
  exit 2
fi
case_dir=$(cd "${SLURM_SUBMIT_DIR}" && pwd -P)
cd "${case_dir}"

if [[ "$#" -eq 0 ]]; then
  printf 'Usage: sbatch run_task.sbatch command [args ...]\n' >&2
  exit 2
fi

exec "$@"
