# Linux/RTX 5080 continuation handoff

This file is the authoritative continuation note for branch
`codex/pipeline-phase-lag-study`. It deliberately distinguishes verified code
from simulation evidence that has not yet been generated.

## State at handoff (updated 2026-08-13)

- Linux build completed with CMake in `mpm_hpc_source/build-pipeline`.
- C++: all 18 registered CTests passed after adding executable thread-limit,
  facet-traction, free-surface, HDF5 and material-point regressions.
- Python: all 102 pipeline tests passed.
- All 20 registered `screen`/`production` JSON files pass the strict validator.
- The earlier combined smoke stopped correctly at `MC_EQ`: before the facet
  correction it reached `max|v|=7.519e-3 m/s`; after correcting the facet load
  mapping it still reached `max|v|=1.8469e-3 m/s > 1e-3 m/s`. The latter is an
  MC handoff diagnostic only, not evidence and not a reason to weaken the
  gate. The smoke workflow no longer creates or runs an MC dynamic stage.
- The current binary passed a fresh four-stage required smoke under label
  `finalfs3-20260813`, with exactly four completion-v3 sentinels. EQ_LS/HS
  reached `max|v|=1.2685e-5/5.9047e-6 m/s`, respectively; both had 7,376
  in-grid particles, zero tensile particles in each audited quiet top row, and
  matching HDF5/VTP state. SANI_LS/HS also retained all 7,376 particles in the
  mesh and ended at `max|v|=4.0161e-5/5.2854e-5 m/s`. These ignored local
  outputs are smoke evidence only, not manuscript evidence. A different
  machine should still rerun the four stages with its own fresh label before
  starting a registered calculation.
- No post-fix registered `screen` result is complete. Do not analyse older or
  partial result directories.
- The local planning PDF `Manuscript-Discussion-0811.pdf` is intentionally not
  committed. The scientific requirements needed to continue are recorded
  below.

The RTX 5080 is not directly used by this CPU/TBB/OpenMP MPM executable. A
faster CPU and adequate RAM are what improve local solve time; set
`MPM_THREADS` to the number of physical CPU cores allocated to a case. The
executable now applies that value to both OpenMP and TBB; the regression test
checks the complete `-p -> IO -> TBB` path.

## Non-negotiable scientific contract

1. Do not modify `Example/slope_infiltration` or other unrelated examples.
2. Run the shortened, isolated four-stage required smoke DAG before any
   registered long calculation. `MC_EQ` is an optional local diagnostic, not a
   prerequisite for declaring the required smoke successful.
3. Equilibrium uses `LinearElastic2D`, pure PIC (`PIC=1`, `APIC=false`), fixed
   pipeline, `dt=1e-4 s`, 40,000 registered steps and Cundall damping `5 s^-1`.
4. `EQ_LS` uses `Sw=0.993`; `EQ_HS` uses `Sw=0.94`.
5. Dynamic cases switch to SANISAND or Mohr-Coulomb and use APIC, zero PIC
   blending and zero damping. MC dynamics must first pass through `MC_EQ`.
6. Never relax the equilibrium stability gate `max|v| <= 1e-3 m/s`.
7. A result is usable only when its hash-bound `pipeline_completion.json`,
   complete VTP grid, unique history CSV, HDF5/VTP cross-check and dependency
   checks all pass. Smoke results are always
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

Expected: CTest `18/18`, Python `102/102`, and 10 validated screen configs.

## 3. Mandatory isolated required smoke DAG

Choose a new lowercase label; never reuse any earlier diagnostic label or the
verified local label `finalfs3-20260813`.

```bash
MPM_THREADS=4 bash run_local_smoke.sh rtx5080-YYYYMMDD
find results/local_smoke/rtx5080-YYYYMMDD \
  -name pipeline_completion.json -print | sort
```

The default command runs four fail-fast 5,000-step stages in dependency order:
`EQ_LS`, `EQ_HS`, `SANI_LS`, `SANI_HS`. Exactly four completion files must be
printed. The completion audit checks particle count, mesh bounds, finite key
fields, the complete VTP grid, unique history CSV, HDF5/VTP cross-check and the
unchanged structured stability QA. Stop and
diagnose any required-stage failure; do not start the registered screen after
a failed required smoke.

To rerun the Mohr-Coulomb handoff as an explicit diagnostic, use a different
fresh label:

```bash
MPM_THREADS=4 bash run_local_smoke.sh --mc-diagnostic mcdiag-YYYYMMDD
```

This first runs the same four required stages and then `MC_EQ`. Any
`max|v| > 1e-3 m/s` or other QA failure exits nonzero and writes no MC_EQ
completion. It never generates or runs an MC dynamic smoke. Its outputs are
always excluded from manuscript evidence, whether it passes or fails. A local
diagnostic failure does not retroactively invalidate the four required smoke
stages, but the formal HM/RM branch still cannot proceed unless its registered
`MC_EQ` stage passes the unchanged gate.

For the current phase-lag paper, proceed with the fail-closed seven-case
SANISAND/phase-only controller after the four required smoke stages pass:

```bash
MPM_THREADS=4 PHASE_SCREEN_MAX_TOTAL_THREADS=8 \
  bash run_local_phase_screen.sh
```

This runs `EQ_LS/EQ_HS -> LS/HS -> HD -> phase transform -> RL/RE`, never
selects MC_EQ/HM/RM, waits each PID separately, verifies the exact VTP grid and
completion after every case, and refuses any non-empty unaudited target.

## 4. Registered phase-only screen calculation

Run inside `tmux` or another persistent terminal. This workstation reached
98--99 degrees C with one 12-thread case, so use the audited controller with at
most two four-thread cases:

```bash
MPM_THREADS=4 PHASE_SCREEN_MAX_TOTAL_THREADS=8 \
  bash run_local_phase_screen.sh
```

The controller alone owns the dependency-safe schedule
`EQ_LS/EQ_HS -> LS/HS -> HD -> phase transform -> RL/RE`. It waits for both
PIDs in every pair, revalidates the exact output grid and hashes, skips only
fully audited completions, and stops on a non-empty unaudited target. It never
selects `MC_EQ`, `HM` or `RM`.

Do not replace this command with the older hand-written full-DAG sequence: the
post-fix 0.5 s `MC_EQ` diagnostic still fails the unchanged velocity gate. Do
not use `--force` or delete results to recover from failure; inspect the log and
validator output first.

## 5. Evidence synthesis and mandatory 2-D figures

After all seven registered phase-screen sentinels exist:

```bash
PIPELINE_LOCAL=1 PYTHON_BIN=python3 bash analyze_results.sh screen phase
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
python3 submit_study_sbatch.py screen --stage phase --prepare --build --analyze \
  --cpus 16 --gres gpu:1 --nodelist gpu30
```

The submitter does not expose partition/account overrides: it always emits and
hard-checks `partition=granularmech` and `account=comgranmech`.
