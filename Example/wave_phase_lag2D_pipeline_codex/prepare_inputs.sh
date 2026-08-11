#!/usr/bin/env bash
# Generate geometry/configuration and validate them inside an HPC4 allocation.
set -eo pipefail
case_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
cd "${case_dir}"

tier=${1:?Usage: prepare_inputs.sh tier label mesh_dir cell_size particle_spacing wave_height}
label=${2:?Missing label}
mesh_dir=${3:?Missing mesh directory}
cell_size=${4:?Missing cell size}
particle_spacing=${5:?Missing particle spacing}
wave_height=${6:?Missing wave height}

module load miniconda3/24.3.0-quc3pyu
eval "$(conda shell.bash hook)"
conda activate cbgeo_tbb

python mesh.py \
  --mesh-size "${cell_size}" \
  --particle-spacing "${particle_spacing}" \
  --output-dir "${mesh_dir}"
python prepare_study.py \
  --tier "${tier}" \
  --mesh-dir "${mesh_dir}" \
  --cell-size "${cell_size}" \
  --particle-spacing "${particle_spacing}" \
  --wave-height "${wave_height}" \
  --label "${label}"

suffix=""
if [[ "${label}" != "baseline" ]]; then
  suffix="_${label}"
fi
python validate_case.py --manifest "configs/${tier}${suffix}/manifest.json"
