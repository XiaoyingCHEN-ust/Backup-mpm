# Linux/RTX 5080 continuation handoff

This file is the authoritative continuation note for branch
`codex/pipeline-phase-lag-study`. It deliberately distinguishes verified code
from simulation evidence that has not yet been generated.

## State at handoff (2026-08-12)

- Linux build completed with CMake in `mpm_hpc_source/build-pipeline`.
- C++: all 7 registered CTests passed.
- Python: all 80 pipeline tests passed.
- All 20 registered `screen`/`production` JSON files pass the strict validator.
- No post-fix smoke run is complete. A local `EQ_LS` smoke was interrupted at
  the user's request before a completion sentinel was written. It is ignored
  and is not manuscript evidence.
- No post-fix registered `screen` result is complete. Do not analyse older or
  partial result directories.
- The local planning PDF `Manuscript-Discussion-0811.pdf` is intentionally not
  committed. The scientific requirements needed to continue are recorded
  below.

The RTX 5080 is not directly used by this CPU/TBB/OpenMP MPM executable. A
faster CPU and adequate RAM are what improve local solve time; set
`MPM_THREADS` to the number of physical CPU cores allocated to a case.

## Non-negotiable scientific contract

1. Do not modify `Example/slope_infiltration` or other unrelated examples.
2. Run a shortened, isolated smoke DAG before any registered long calculation.
3. Equilibrium uses `LinearElastic2D`, pure PIC (`PIC=1`, `APIC=false`), fixed
   pipeline, `dt=1e-4 s`, 40,000 registered steps and Cundall damping `5 s^-1`.
4. `EQ_LS` uses `Sw=0.993`; `EQ_HS` uses `Sw=0.94`.
5. Dynamic cases switch to SANISAND or Mohr-Coulomb and use APIC, zero PIC
   blending and zero damping. MC dynamics must first pass through `MC_EQ`.
6. Never relax the equilibrium stability gate `max|v| <= 1e-3 m/s`.
7. A result is usable only when its hash-bound `pipeline_completion.json`,
   final VTP/HDF5 and dependency checks all pass. Smoke results are always
   excluded from manuscript evidence.
8. The paper's hypothesis is that hydraulic phase lag can facilitate
   liquefaction and that SANISAND exposes cyclic mechanisms missed by a simple
   Mohr-Coulomb ablation. Report this only if the pre-registered gates support
   it; do not force the conclusion.
9. The final result package must include registered-time, common-scale 2-D
   comparisons. In particular, use excess pressure relative to the checkpoint,
   mask `Rsigma` where initial confinement is below 100 Pa, plot the actual
   moving pipe position, and label RL/RE differences as ID-matched Lagrangian
   differences.

## 1. Clone and verify

```bash
git clone --branch codex/pipeline-phase-lag-study --single-branch \
  https://github.com/XiaoyingCHEN-ust/Backup-mpm.git Backup-mpm
cd Backup-mpm
git status --short --branch
git rev-parse --short HEAD
uname -a
```

If system dependencies are missing on Ubuntu, install the equivalent of C++17,
CMake, Boost filesystem/system, Eigen3, HDF5 C/C++, TBB and VTK. Python needs
NumPy; final figures also need Matplotlib.

## 2. Build and rerun both test suites

Use a conservative build parallelism first; this source tree can consume
substantial RAM while compiling.

```bash
cmake -S mpm_hpc_source -B mpm_hpc_source/build-pipeline \
  -DCMAKE_BUILD_TYPE=Release -DMPM_BUILD_TESTING=ON -DMPM_BUILD_LIB=OFF
cmake --build mpm_hpc_source/build-pipeline --parallel 2
ctest --test-dir mpm_hpc_source/build-pipeline --output-on-failure

cd Example/wave_phase_lag2D_pipeline_codex
python3 -m unittest discover -s tests -p 'test_*.py'
python3 validate_case.py --manifest configs/screen/manifest.json
```

Expected: CTest `7/7`, Python `80/80`, and 10 validated screen configs.

## 3. Mandatory isolated smoke DAG

Choose a new lowercase label; never reuse the interrupted label
`audit-20260812a`.

```bash
MPM_THREADS=8 bash run_local_smoke.sh rtx5080-20260812
find results/local_smoke/rtx5080-20260812 \
  -name pipeline_completion.json -print | sort
```

The command runs six fail-fast 5,000-step stages in dependency order:
`EQ_LS`, `EQ_HS`, `SANI_LS`, `SANI_HS`, `MC_EQ`, `MC_DYNAMIC`. Exactly six
completion files must be printed. The completion audit checks particle count,
mesh bounds, finite key fields, final VTP/HDF5 and the unchanged equilibrium
velocity/displacement/porosity gate. Stop and diagnose any failure; do not
start the registered screen after a failed smoke.

## 4. Registered screen calculation

Run inside `tmux` or another persistent terminal. The examples below assume a
16-physical-core CPU. Do not oversubscribe: reduce per-case threads if several
cases run together.

First run both registered equilibria:

```bash
mkdir -p logs/local
MPM_THREADS=8 bash run_case_local.sh configs/screen/01_EQ_LS.json \
  >logs/local/EQ_LS.log 2>&1 & p1=$!
MPM_THREADS=8 bash run_case_local.sh configs/screen/01_EQ_HS.json \
  >logs/local/EQ_HS.log 2>&1 & p2=$!
wait "$p1" "$p2"
```

Only after both stable sentinels exist, run physical SANISAND cases, the fixed
driver, and MC handoff relaxation. Four concurrent cases below use four CPU
threads each:

```bash
MPM_THREADS=4 bash run_case_local.sh configs/screen/02_LS.json \
  >logs/local/LS.log 2>&1 & p1=$!
MPM_THREADS=4 bash run_case_local.sh configs/screen/02_HS.json \
  >logs/local/HS.log 2>&1 & p2=$!
MPM_THREADS=4 bash run_case_local.sh configs/screen/02_HD.json \
  >logs/local/HD.log 2>&1 & p3=$!
MPM_THREADS=4 bash run_case_local.sh configs/screen/02_MC_EQ.json \
  >logs/local/MC_EQ.log 2>&1 & p4=$!
wait "$p1" "$p2" "$p3" "$p4"
```

Create and validate the phase-erased replay only after HD completes:

```bash
PIPELINE_LOCAL=1 PYTHON_BIN=python3 bash run_phase_control.sh screen 10.4
```

Then run the constitutive and replay controls:

```bash
MPM_THREADS=4 bash run_case_local.sh configs/screen/03_HM.json \
  >logs/local/HM.log 2>&1 & p1=$!
MPM_THREADS=4 bash run_case_local.sh configs/screen/04_RL.json \
  >logs/local/RL.log 2>&1 & p2=$!
MPM_THREADS=4 bash run_case_local.sh configs/screen/04_RM.json \
  >logs/local/RM.log 2>&1 & p3=$!
MPM_THREADS=4 bash run_case_local.sh configs/screen/04_RE.json \
  >logs/local/RE.log 2>&1 & p4=$!
wait "$p1" "$p2" "$p3" "$p4"
```

Do not use `--force` or delete results to recover from failure. Read the failing
log and validator message first. The per-result lock prevents two writers from
using the same UUID target.

## 5. Evidence synthesis and mandatory 2-D figures

After all ten registered screen sentinels exist:

```bash
PIPELINE_LOCAL=1 PYTHON_BIN=python3 bash analyze_results.sh screen full
```

Review, at minimum:

- `analysis/screen/manuscript_evidence.json`
- `analysis/screen/manuscript_evidence.md`
- `analysis/screen/figures/figure_manifest.json`
- all PNG/PDF files under `analysis/screen/figures/`
- `pressure_databases/screen/phase_erased/phase_erased_metadata.json`

The phase-only claim requires pressure mean/amplitude agreement within 2%, at
least 5 degrees of lag reduction, a hydraulic-trigger difference beyond one
particle-area by one saved-frame resolution, and a realised same-particle/time
`IF >= 1` plus `Rsigma <= 0.05` chain. The SANISAND discussion must use cyclic
`q-p'`, plastic accumulation and contact evolution; internal Alpha/Z fabric is
not exported and must not be described as directly observed.

Only proceed to the registered production tier after the complete screen and
all evidence gates have been reviewed. Populate the remaining manuscript
placeholders from audited JSON/CSV outputs, not from partial VTP files.

## Known implementation limitation

The current HDF5 schema does not preserve the initial liquid/gas/pore-pressure
reference fields or initial vertical effective stress for arbitrary same-stage
continuation. Registered EQ-to-dynamic transitions intentionally rebase those
references and are supported. Do not claim that an arbitrary mid-stage restart
is bitwise/physically identical without first extending the schema and adding
an explicit `rebase` versus `preserve` policy.

## HPC4 fallback (not needed for the requested local run)

The fixed allocation contract remains:

```bash
python3 submit_study_sbatch.py screen --prepare --build --analyze \
  --cpus 16 --gres gpu:1 --nodelist gpu30
```

The submitter does not expose partition/account overrides: it always emits and
hard-checks `partition=granularmech` and `account=comgranmech`.
