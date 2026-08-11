#!/usr/bin/env bash
# Validate one generated configuration group inside an HPC4 Slurm batch allocation.
set -eo pipefail
case_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
cd "${case_dir}"
group=${1:?Usage: validate_inputs.sh config_group}

module load miniconda3/24.3.0-quc3pyu
eval "$(conda shell.bash hook)"
conda activate cbgeo_tbb
python validate_case.py --manifest "configs/${group}/manifest.json"
