# Moving impermeable pipeline rerun

## What was wrong

The failed diagnostic run did not show physical large deformation before the
abort. At step 26000 the maximum soil displacement was only
`4.29e-5 m = 0.00215` background cells. The initialization then aborted at
step 27742 because a particle left the mesh.

The visible pipe-region pattern was mainly a boundary-discretisation artifact:
22 fixed Cartesian background nodes (old node set 4) represented a circular
pipe. Their radii ranged from 0.050 to 0.0671 m for a real 0.060 m pipe, so the
stair-step wall imposed a grid-shaped liquid/gas constraint and could not
follow a moving pipe. The earlier analyzer also used an overly wide wall ring
and mixed 132 particles into its statistics.

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

The patch application is idempotent: rerunning it reports that the patch is
already applied only after both required source markers are found. Empty or
truncated patch files are rejected. The two calculation jobs request partition
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
`results/Initial_094_1E11_pipeline_moving_noflux_v2`

Dynamic:
`results/Wave2D_SN_094_1E11_pipeline_moving_noflux_v2`

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

This v2 rerun uses the target permeability in both stages. If a later
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

The incomplete v1 initialization downloaded on 2026-08-20 was generated by a
stale executable: its saved wall-normal liquid and gas velocities were not
zero even though the JSON enabled `phase_no_flux`. It is retained only as a
failure diagnostic and must not be used for physical interpretation. The v2
run is acceptable only if the build/source checks pass and the post-step-0
wall-normal relative velocities are near floating-point noise.
