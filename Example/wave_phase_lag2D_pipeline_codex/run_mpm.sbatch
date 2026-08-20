#!/bin/bash
#SBATCH --job-name=pipe_094_1e11_noflux
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=24:00:00
#SBATCH --partition=granularmech
#SBATCH --account=comgranmech
#SBATCH --output=/home/xchenjm/chen/MPM/Example/wave_phase_lag2D_pipeline_codex/mpm_sanisand_%j.out
#SBATCH --error=/home/xchenjm/chen/MPM/Example/wave_phase_lag2D_pipeline_codex/mpm_sanisand_%j.err

# Conda's activate.d scripts inspect optional variables such as ADDR2LINE, so
# nounset (-u) must remain disabled during environment activation.
set -eo pipefail
cd "$SLURM_SUBMIT_DIR"

module load miniconda3/24.3.0-quc3pyu
eval "$(conda shell.bash hook)"
conda activate cbgeo_tbb

python check_pipeline_wall.py

echo "job=${SLURM_JOB_ID} host=$(hostname) started=$(date --iso-8601=seconds)"
echo "workdir=$(pwd)"
runtime_log="${SLURM_SUBMIT_DIR}/mpm_sanisand_runtime_${SLURM_JOB_ID}.log"
mpm_executable=/home/xchenjm/chen/MPM/mpm/build/mpm
solver_source=/home/xchenjm/chen/MPM/mpm/include/solvers/thm_mpm_explicit_threephase_lag.tcc
if [[ ! -x "$mpm_executable" || "$solver_source" -nt "$mpm_executable" ]]; then
  echo "MPM binary is missing or older than the no-flux solver source." >&2
  echo "Run: cd /home/xchenjm/chen/MPM/mpm/build && make -j8" >&2
  exit 2
fi

set +e
stdbuf -oL -eL "$mpm_executable" \
  -p 32 -f ./ -i mpm-3p.json 2>&1 | tee "$runtime_log"
status=${PIPESTATUS[0]}
set -e

echo "MPM_EXIT_CODE=${status} finished=$(date --iso-8601=seconds)"
exit "$status"
