#!/usr/bin/env bash
# Generate and run the exact six-stage local smoke DAG. Any failed stage stops
# the sequence, so no dependent checkpoint consumer can start after failure.
set -euo pipefail

case_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)
cd "${case_dir}"

if [[ $# -ne 1 ]]; then
  printf 'Usage: run_local_smoke.sh unique_lowercase_label\n' >&2
  exit 2
fi
label=$1
if ! [[ "${label}" =~ ^[a-z0-9][a-z0-9_-]{0,31}$ ]]; then
  printf 'Invalid smoke label: %s\n' "${label}" >&2
  exit 2
fi

python_bin=${PYTHON_BIN:-python3}
"${python_bin}" prepare_local_smoke.py --label "${label}"
manifest="configs/local_smoke/${label}/manifest.json"
"${python_bin}" validate_case.py --local-smoke --manifest "${manifest}"

for filename in \
  01_EQ_LS.json \
  01_EQ_HS.json \
  02_SANI_LS.json \
  02_SANI_HS.json \
  02_MC_EQ.json \
  03_MC_DYNAMIC.json
do
  bash run_case_local.sh --smoke "configs/local_smoke/${label}/${filename}"
done
