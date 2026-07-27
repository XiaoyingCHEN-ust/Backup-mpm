#!/bin/bash

set -euo pipefail

case_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${case_dir}"

reference="stability_results/siemens2013-stability-r19b-open-3s-open_2p5e-06/particle1200000.vtp"
if [[ ! -f "${reference}" ]]; then
  echo "Continuous r19b reference is missing: ${reference}" >&2
  exit 4
fi

for uuid in \
  siemens2013-r20-checkpoint-open-segment-a \
  siemens2013-r20-checkpoint-open-segment-b; do
  if [[ -e "stability_results/${uuid}" ]]; then
    echo "Refusing to mix a rerun with existing results: stability_results/${uuid}" >&2
    exit 4
  fi
done

python prepare_checkpoint_validation.py

job_a_raw="$(sbatch --parsable --export=ALL,R20_SEGMENT=a run_checkpoint_segment_hpc4.sh)"
job_a="${job_a_raw%%;*}"
job_b_raw="$(sbatch --parsable --dependency="afterok:${job_a}" \
  --export=ALL,R20_SEGMENT=b run_checkpoint_segment_hpc4.sh)"
job_b="${job_b_raw%%;*}"

echo "Submitted r20 segment A: ${job_a}"
echo "Submitted r20 segment B: ${job_b} (afterok:${job_a})"
echo "Segment B will run the final continuous-versus-resumed comparison."
