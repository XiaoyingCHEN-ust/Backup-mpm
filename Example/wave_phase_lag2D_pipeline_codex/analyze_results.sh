#!/usr/bin/env bash
# Reduce one completed study group inside an HPC4 Slurm batch allocation.
set -eo pipefail

case_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
cd "${case_dir}"

group=${1:?Usage: analyze_results.sh group physical-or-full}
stage=${2:?Usage: analyze_results.sh group physical-or-phase-or-full}
if [[ "${stage}" != "physical" && "${stage}" != "phase" && "${stage}" != "full" ]]; then
  printf 'Analysis stage must be physical, phase or full, got: %s\n' "${stage}" >&2
  exit 2
fi

if [[ "${PIPELINE_LOCAL:-0}" != "1" ]]; then
  module load miniconda3/24.3.0-quc3pyu
  eval "$(conda shell.bash hook)"
  conda activate cbgeo_tbb
fi
python_bin=${PYTHON_BIN:-python3}

configs=(
  "configs/${group}/02_LS.json"
  "configs/${group}/02_HS.json"
)
if [[ "${stage}" == "phase" ]]; then
  configs+=(
    "configs/${group}/04_RL.json"
    "configs/${group}/04_RE.json"
  )
elif [[ "${stage}" == "full" ]]; then
  configs+=(
    "configs/${group}/03_HM.json"
    "configs/${group}/04_RL.json"
    "configs/${group}/04_RM.json"
    "configs/${group}/04_RE.json"
  )
fi

"${python_bin}" analyze_study.py "${configs[@]}" --output "analysis/${group}"

if [[ "${stage}" == "phase" || "${stage}" == "full" ]]; then
  synthesis_options=()
  plotting_options=()
  if [[ "${stage}" == "phase" ]]; then
    synthesis_options+=(--phase-only)
    plotting_options+=(--phase-only)
  fi
  "${python_bin}" synthesize_manuscript_evidence.py \
    --analysis-dir "analysis/${group}" \
    --phase-metadata \
      "pressure_databases/${group}/phase_erased/phase_erased_metadata.json" \
    --output-json "analysis/${group}/manuscript_evidence.json" \
    --output-markdown "analysis/${group}/manuscript_evidence.md" \
    "${synthesis_options[@]}"
  if "${python_bin}" -c 'import matplotlib' >/dev/null 2>&1; then
    "${python_bin}" plot_manuscript_figures.py \
      --analysis-dir "analysis/${group}" \
      --output-dir "analysis/${group}/figures" \
      "${plotting_options[@]}"
  else
    printf '%s\n' \
      'Matplotlib unavailable; evidence tables were written, figure drafts skipped.' >&2
  fi
fi
