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
  "configs/${group}/03_HM.json"
)
if [[ "${stage}" == "full" ]]; then
  configs+=(
    "configs/${group}/04_RL.json"
    "configs/${group}/04_RM.json"
    "configs/${group}/04_RE.json"
  )
fi

python analyze_study.py "${configs[@]}" --output "analysis/${group}"

if [[ "${stage}" == "full" ]]; then
  python synthesize_manuscript_evidence.py \
    --analysis-dir "analysis/${group}" \
    --phase-metadata \
      "pressure_databases/${group}/phase_erased/phase_erased_metadata.json" \
    --output-json "analysis/${group}/manuscript_evidence.json" \
    --output-markdown "analysis/${group}/manuscript_evidence.md"
  if python -c 'import matplotlib' >/dev/null 2>&1; then
    python plot_manuscript_figures.py \
      --analysis-dir "analysis/${group}" \
      --output-dir "analysis/${group}/figures"
  else
    printf '%s\n' \
      'Matplotlib unavailable; evidence tables were written, figure drafts skipped.' >&2
  fi
fi
