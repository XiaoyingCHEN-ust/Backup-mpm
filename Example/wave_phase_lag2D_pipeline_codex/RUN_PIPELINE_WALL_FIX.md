# Moving impermeable pipeline rerun

## What was wrong

The downloaded v2 initialization completed, but it was not a valid equilibrium.
Although `initial_liquid_pressures.txt` contains the intended hydrostatic field
(about 9.76 kPa at the base), the three-phase particle's empty
`initial_pore_pressure()` override discarded it. The first VTP therefore had
zero liquid pressure through most of the bed. During four seconds of numerical
relaxation, the mean vertical skeleton stress changed from about -2.03 kPa
(compression) to +2.32 kPa (tension), and the compressive-particle fraction
fell from 1.0 to 0.00027.

The following SANISAND run inherited that tensile checkpoint, stayed on its
low-confinement elastic fallback (`eps_p_q=0` everywhere), and aborted at
0.7266 s with `SANISAND elastic fallback became invalid`. This was before pipe
release and was not physical large deformation.

The older fixed-wall model had an additional discretisation artifact: 22 fixed
Cartesian nodes represented the circle with radii from 0.050 to 0.0671 m for a
real 0.060 m pipe. Its wall-saturation spread was larger than the moving
particle-level no-flux result, so reverting to the old Euler/node wall would
not fix the hydrostatic error.

## What changed

- The fixed node-set-4/Euler-angle wall has been removed.
- The solver now projects liquid and gas particle velocity onto the tangent of
  the actual circular wall. The condition is
  `(v_phase - v_pipe) dot n = 0`, and follows pipe translation/rotation.
- The solver also supports projecting liquid/gas affine maps if APIC is enabled
  in a later sensitivity run. The present diagnostic keeps APIC disabled.
- Initialization remains `LinearElastic2D`, pure `PIC=1`,
  `APIC=false`, fixed pipe and Cundall damping. Its permeability is
  `9.79e-12 m2`, identical to the subsequent dynamic stage, so checkpoint
  handoff does not introduce a permeability jump. The wave stage resumes the
  fresh checkpoint, returns to `SANISAND2D`, and uses
  `PIC=0, APIC=false` with the same permeability.
- Fresh UUIDs prevent old partial results from being reused.
- A second source patch preserves the hydrostatic liquid pressure until the
  liquid/gas properties are available, then derives gas pressure by adding the
  saturation-dependent capillary suction.
- The submerged seabed pressure boundary now adds water/wave pressure equally
  to liquid and gas phases instead of forcing suction to zero. This prevents
  the top two particle rows next to the shallow pipe crown from jumping from
  `Sr=0.94` to about 0.999 in the first time step.
- `check_initial_state.py` gates the SANISAND stage. It rejects a missing
  hydrostatic field, a mostly tensile skeleton checkpoint, non-decayed solid
  velocity, a final pressure RMSE above 500 Pa, pipe-wall saturation more than
  0.02 from the prescribed 0.94, or nonzero pipe-wall phase flux.

GIMP is deliberately not enabled in this diagnostic. The downloaded result
never approached one background-cell displacement, so changing interpolation
and mesh at the same time would hide whether the wall correction worked. Use
GIMP only after a corrected run actually develops mesh-crossing soil motion.

## Run on HPC4

Upload this complete directory to:

`/home/xchenjm/chen/MPM/Example/wave_phase_lag2D_pipeline_codex`

Apply the solver patch and compile interactively:

```bash
cd /home/xchenjm/chen/MPM/Example/wave_phase_lag2D_pipeline_codex
bash apply_mpm_patch.sh

cd /home/xchenjm/chen/MPM/mpm/build
make -j8
./mpm --help

cd /home/xchenjm/chen/MPM/Example/wave_phase_lag2D_pipeline_codex
python mesh.py
python check_pipeline_wall.py
```

Then submit only initialization and wave calculation as a dependency chain:

```bash
init_job=$(sbatch --parsable run_mpm_ini.sbatch)
init_job=${init_job%%;*}

wave_job=$(sbatch --parsable --dependency=afterok:${init_job} run_mpm.sbatch)
wave_job=${wave_job%%;*}

printf 'initial=%s wave=%s\n' "$init_job" "$wave_job"
```

Patch application is idempotent: the moving no-flux and hydrostatic-pressure
patches are verified/applied independently. Empty or truncated patches are
rejected. The two calculation jobs request partition
`granularmech`, account `comgranmech` and 32 CPUs.

Monitor without hiding failed batch steps:

```bash
squeue -j "${init_job},${wave_job}" \
  -o "%.10i %.22j %.2t %.10M %.28R"
sacct -j "${init_job},${wave_job}" \
  --format=JobID,JobName%24,State,ExitCode,Elapsed,NodeList -n -P
```

Do not submit the wave job manually if initialization fails; inspect
`mpm_initial_JOBID.err` and `mpm_initial_runtime_JOBID.log`.

## Outputs and acceptance check

Initialization:
`results/Initial_094_1E11_pipeline_hydrostatic_noflux_v3`

Dynamic:
`results/Wave2D_SN_094_1E11_pipeline_hydrostatic_noflux_v3`

The initialization job writes `initial_state_qa.json` and returns nonzero if
the checkpoint is unacceptable. It must contain `"passed": true` before the
dependent dynamic job can start.

After the dynamic job:

```bash
python analyze_pipeline_wall.py
```

The analyzer selects only particles inside
`R + sqrt(2) * particle_radius`, uses the moving pipe velocity, and writes
`pipeline_wall_diagnostics.csv`. For saved frames after step 0, the liquid
and gas relative normal velocities should be near floating-point noise. The
pipe-region saturation standard deviation should lose the old cross/grid
imprint. Soil displacement divided by the 0.02 m cell width remains the gate
for deciding whether a later GIMP rerun is justified.

The legacy `liquid_seepage_velocities` field is not used as this acceptance
metric because its current implementation omits the gravity term and is
nonzero even under a hydrostatic pressure gradient.

This v3 rerun uses the target permeability in both stages. If a later
equilibration sensitivity is needed, change only `mpm-initial.json` and use a
fresh UUID; do not silently change the scientific target in `mpm-3p.json`.

## Result retention and Git synchronization

Raw results, Slurm logs and Python caches are intentionally ignored by Git.
For each new run, keep the following on HPC4 until the result has been
downloaded and checked:

- all `particle*.vtp` and matching `pipeline*.vtp` frames from both stages;
- the initialization final `particles40000.h5` checkpoint;
- `pipeline-history*.csv`;
- Slurm `.out`, `.err`, and `*_runtime_*.log` files.

The complete particle time series is needed to detect when wall leakage,
saturation oscillation or displacement first develops. Do not prune it before
`analyze_pipeline_wall.py` has run. After a successful audit, duplicate
`node*.vtp`, `mesh*.vtp`, and `geometry*.vtp` files may be deleted; keep
the diagnostic CSV, initial/release/peak/final particle and pipeline frames,
the final initialization checkpoint, and the logs. These retained results stay
local/HPC-only and are not committed to GitHub.

The v1 and v2 results are failure diagnostics only and must not be used for
physical interpretation. The v3 run is acceptable only when source/build
checks pass, `initial_state_qa.json` passes, skeleton stress stays predominantly
compressive, and post-step-0 wall-normal relative phase velocities remain near
floating-point noise.
