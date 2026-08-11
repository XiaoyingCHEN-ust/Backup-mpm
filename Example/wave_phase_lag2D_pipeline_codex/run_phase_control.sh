#!/usr/bin/env bash
# Transform and validate the phase-erased database inside a Slurm batch allocation.
set -eo pipefail

case_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
cd "${case_dir}"

group=${1:?Usage: run_phase_control.sh screen_or_production_group fit_end_s}
fit_end=${2:?Usage: run_phase_control.sh group fit_end_s}

module load miniconda3/24.3.0-quc3pyu
eval "$(conda shell.bash hook)"
conda activate cbgeo_tbb

python phase_controls.py \
  --source-dir "pressure_databases/${group}/lagged" \
  --source-prefix pressure \
  --output-dir "pressure_databases/${group}/phase_erased" \
  --output-prefix phase_erased \
  --period 1.3 \
  --fit-start 2.6 \
  --fit-end "${fit_end}" \
  --ramp-time 1.3

python validate_case.py "configs/${group}/04_RE.json" --runtime
