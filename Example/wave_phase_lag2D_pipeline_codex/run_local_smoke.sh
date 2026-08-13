#!/usr/bin/env bash
# Generate and run the four required local smoke stages. MC_EQ is available
# only as an explicit fail-fast diagnostic; no MC dynamic smoke is defined.
set -euo pipefail

case_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)
cd "${case_dir}"

include_mc_diagnostic=0
if [[ ${1:-} == "--mc-diagnostic" ]]; then
  include_mc_diagnostic=1
  shift
fi
if [[ $# -ne 1 ]]; then
  printf 'Usage: run_local_smoke.sh [--mc-diagnostic] unique_lowercase_label\n' >&2
  exit 2
fi
label=$1
if ! [[ "${label}" =~ ^[a-z0-9][a-z0-9_-]{0,31}$ ]]; then
  printf 'Invalid smoke label: %s\n' "${label}" >&2
  exit 2
fi

python_bin=${PYTHON_BIN:-python3}
prepare_args=(--label "${label}")
if [[ "${include_mc_diagnostic}" -eq 1 ]]; then
  prepare_args+=(--include-mc-diagnostic)
fi
"${python_bin}" prepare_local_smoke.py "${prepare_args[@]}"
manifest="configs/local_smoke/${label}/manifest.json"
"${python_bin}" prepare_local_smoke.py --validate-manifest "${manifest}"
"${python_bin}" validate_case.py --local-smoke --manifest "${manifest}"

for filename in \
  01_EQ_LS.json \
  01_EQ_HS.json \
  02_SANI_LS.json \
  02_SANI_HS.json
do
  bash run_case_local.sh --smoke "configs/local_smoke/${label}/${filename}"
done

if [[ "${include_mc_diagnostic}" -eq 1 ]]; then
  printf '%s\n' \
    'Running opt-in MC_EQ diagnostic; it is excluded from manuscript evidence.'
  # run_case_local publishes no completion and exits nonzero if the unchanged
  # max|v| <= 1e-3 m/s gate or any other structured QA check fails.
  bash run_case_local.sh --smoke \
    "configs/local_smoke/${label}/02_MC_EQ.json"
fi
