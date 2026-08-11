#!/usr/bin/env bash
# Reduce one completed study group inside an HPC4 Slurm batch allocation.
set -eo pipefail

case_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
cd "${case_dir}"

group=${1:?Usage: analyze_results.sh group physical-or-full}
stage=${2:?Usage: analyze_results.sh group physical-or-full}
if [[ "${stage}" != "physical" && "${stage}" != "full" ]]; then
  printf 'Analysis stage must be physical or full, got: %s\n' "${stage}" >&2
  exit 2
fi

module load miniconda3/24.3.0-quc3pyu
eval "$(conda shell.bash hook)"
conda activate cbgeo_tbb

configs=(
  "configs/${group}/02_LS.json"
  "configs/${group}/02_HS.json"
  "configs/${group}/02_HM.json"
)
if [[ "${stage}" == "full" ]]; then
  configs+=(
    "configs/${group}/04_RL.json"
    "configs/${group}/04_RM.json"
    "configs/${group}/04_RE.json"
  )
fi

python analyze_study.py "${configs[@]}" --output "analysis/${group}"
