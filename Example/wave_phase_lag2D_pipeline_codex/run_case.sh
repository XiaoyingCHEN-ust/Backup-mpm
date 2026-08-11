#!/usr/bin/env bash
# Execute one case inside an HPC4 Slurm batch allocation.
set -eo pipefail

case_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)
cd "${case_dir}"

config=${1:?Usage: run_case.sh configs/group/case.json}
config_abs=$(realpath -e -- "${config}")
case "${config_abs}" in
  "${case_dir}"/*) config_relative=${config_abs#"${case_dir}"/} ;;
  *)
    printf 'Configuration must be inside the pipeline case directory: %s\n' \
      "${config_abs}" >&2
    exit 2
    ;;
esac
threads=${SLURM_CPUS_PER_TASK:-${MPM_THREADS:-16}}
mpm_bin=${MPM_BIN:-"${case_dir}/../../mpm_hpc_source/build-pipeline/mpm"}

# Conda activation scripts inspect optional variables, so nounset remains off.
module load miniconda3/24.3.0-quc3pyu
eval "$(conda shell.bash hook)"
conda activate cbgeo_tbb

if [[ ! -x "${mpm_bin}" ]]; then
  printf 'MPM executable is missing or not executable: %s\n' "${mpm_bin}" >&2
  exit 2
fi

python validate_case.py "${config_abs}" --runtime --clear-completion
export OMP_NUM_THREADS="${threads}"
export TBB_NUM_THREADS="${threads}"
"${mpm_bin}" -p "${threads}" -f ./ -i "${config_relative}"
python validate_case.py "${config_abs}" --write-completion
