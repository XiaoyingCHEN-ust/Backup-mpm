# Phase-lag-controlled pipeline instability study

This directory contains the final case-study workflow for the manuscript.  The
case is a shallow-buried, empty PE100 SDR17 pipeline that becomes a sensitive
indicator of support loss when the surrounding seabed liquefies.  The proposed
contribution is:

> Phase-resolved three-phase MPM linking hydraulic phase lag to post-liquefaction
> pipeline migration, with a phase-erased numerical counterfactual and a cyclic
> constitutive-model ablation.

The code does **not** assume that SANISAND must predict more displacement than
Mohr-Coulomb.  The comparison tests whether cyclic memory, state dependence,
fabric evolution and cyclic mobility change the trigger, accumulation rate and
post-trigger path.  A result is useful even when final displacement ordering is
not monotonic.

## Study matrix

| Code | Phase input | Soil model | Status and purpose |
|---|---|---|---|
| LS | Fully coupled, near-saturated `Sw=0.999` | SANISAND | Physical low-lag reference; report the measured lag, never call it exact zero |
| HS | Fully coupled, `Sw=0.94` | SANISAND | Physical high-lag main case |
| HM | Fully coupled, `Sw=0.94` | Mohr-Coulomb | Simple constitutive ablation against HS |
| HD | Fully coupled, fixed pipe, `Sw=0.94` | SANISAND | Generates the pressure database |
| RL | Original HD pressure history, free pipe | SANISAND | One-way replay control |
| RM | Original HD pressure history, free pipe | Mohr-Coulomb | Matched-pressure constitutive ablation against RL |
| RE | Same mean/amplitude as RL, fundamental lag erased | SANISAND | One-way numerical counterfactual; **not** a physical fully coupled solution |

LS–HS is the physical comparison, but saturation also changes attenuation and
drainage.  RL–RE is therefore the strict phase-only diagnostic.  RL–RM applies
the same pressure database to both constitutive models and is the primary
constitutive ablation; HS–HM is retained only as fully coupled context.

All principal cases use `kappa=9.79e-12 m2` (approximately
`Kh=9.60e-5 m/s` for water), porosity `n=0.485`, wave `H=0.12 m`,
`T=1.3 s`, pipe `D=0.12 m`, and `cover/D=0.25`.  The pipe invert is 0.15 m
below the bed, deliberately intersecting the deepest phase-lag-sensitive zone
identified in the preceding manuscript section.  MC uses `c=0`, `psi=0` and
`phi=31.1474 degrees`, obtained from the SANISAND compression critical-state
slope `M=1.25`.  To reduce elastic-stiffness confounding, its constant tangent
is matched to the implemented SANISAND tangent at the pipe-relevant reference
state `p'=3 kPa`, `n=0.485`: `G=4.58693 MPa`, `K=2.98605 MPa`, giving
`E=9.10081 MPa` and `nu=-0.0079625`.  This is a local reference-state match;
SANISAND remains pressure dependent.  The generated JSON records the formulas
and calibration values.

## Server workflow

Run from this directory after checking out the branch.  Windows does not need
to build the C++ code.  The top-level Python runner is lightweight; all input
generation, validation, compilation and simulation work is launched inside an
HPC4 allocation.

```bash
# Submit the complete eight-cycle screen dependency chain and return immediately.
python submit_study_sbatch.py screen --prepare --build --analyze \
  --cpus 16 --gres gpu:1 --nodelist gpu30
```

The two equilibrium jobs deliberately use `LinearElastic2D`, `PIC=1.0`, a
fixed pipe and positive Cundall damping.  `mesh.py` also generates a
hydrostatic liquid-pressure profile for every particle plus separate LS/HS
effective self-weight stress fields.  Validation checks the single-column MPM
pressure format and the stress/pressure values before a job is submitted.  The
dynamic jobs then restore the successful equilibrium HDF5 and switch to
SANISAND or Mohr-Coulomb as registered in the study matrix.

The submitter uses `sbatch --parsable` and records the dependency graph and job
IDs under `submissions/`.  Every job is forced onto
`--partition=granularmech --account=comgranmech`; the defaults also request
`gpu30`, 16 CPU and `gres:gpu:1`.  It submits both equilibria, physical cases,
the fixed-pipe driver, phase transform, replays and analysis with `afterok`
dependencies, then returns so the terminal may be disconnected or reused for
Code Tunnel.  Completed final VTP/HDF5 outputs are skipped unless `--force` is
specified.  A file is not considered complete merely because it exists:
`run_case.sh` publishes `pipeline_completion.json` atomically only after the
solver exits successfully and the final VTP, required HDF5, configuration hash
and (for HD) complete pressure-database frames have been audited.  Resubmission
rechecks those sizes and SHA-256 hashes.  Dynamic-case sentinels also fingerprint
the exact equilibrium HDF5/QA VTP and any replay pressure database (including
phase metadata for RE).  If an equilibrium, HD driver or phase transform is
resubmitted, the DAG marks all of its descendants dirty and reruns them even if
their old output sentinel still exists.  A scheduled rerun clears its prior
sentinel before the solver starts, so a failed forced run cannot leave a stale
success marker.  The phase-erased control is reused only
when its metadata hashes still match both the transformed output and the current
lagged source database.  The cases use the locally built
`../../mpm_hpc_source/build-pipeline/mpm`; set `MPM_BIN` before the run to
use another executable.  Every dynamic job verifies that the requested HDF5
checkpoint exists and that its final equilibrium VTP has finite fields,
`max|v| <= 1e-3 m/s`, displacement below one coarse cell, and valid porosity;
failure stops the dependent job.  The server build also rejects a legacy
159-field SANISAND checkpoint that declares more history variables than it
actually stores; the registered workflow always creates a fresh linear-elastic
equilibrium checkpoint before the SANISAND/MC handoff.
Each batch job writes `logs/<job-name>-<job-id>.out/.err`.
The submit command returns after queueing; monitor later with `squeue -u "$USER"`
or inspect the saved `submissions/<group>-<UTC timestamp>.json` while using Code
Tunnel for other work.

To resubmit only the reduction after jobs finish:

```bash
sbatch --partition=granularmech --account=comgranmech \
  --nodes=1 --ntasks=1 --cpus-per-task=16 --gres=gpu:1 \
  --nodelist=gpu30 --time=02:00:00 --job-name=plp_screen_analyze \
  run_task.sbatch bash analyze_results.sh screen full
```

Please return `analysis/screen/`, `material_preflight_3kPa.csv`, the newest
submission JSON, the
`phase_erased_metadata.json` file, and the SLURM logs.  The large VTP/result
directories and binary pressure databases do not need to be transferred or
committed.  The low-confinement CSV is a single-material-point diagnostic at
the approximately 3 kPa effective stress relevant around the pipe; inspect its
stress reversals and plastic/fabric accumulation before interpreting the field
comparison.

For only the fully coupled physical screen (without the replay diagnostic), use
`--stage physical`.  Only proceed to the 18-cycle production tier after
inspecting the complete screen:

```bash
python submit_study_sbatch.py production --analyze \
  --cpus 16 --gres gpu:1 --nodelist gpu30
```

## Pre-registered decision gates

The pipeline is fixed for three wave periods and then released.  The screen is
five additional periods; production is fifteen additional periods.

- A strong MPM large-deformation result requires either maximum
  `|uy|/D >= 0.5`, or pipe breakout plus a documented loss/reformation of the
  soil-contact network.  Also report post-release soil `max|u|/h` and the
  fraction of material points crossing at least one background cell, both
  globally and for the pre-registered support-zone cohort.  Cross-cell travel
  documents mesh-distortion burden but is not itself a strain measure.
  Otherwise describe the case as small deformation and do not claim that
  conventional FEM is incapable.
- If baseline HS does not cross the gate, retain and report it.  A registered
  stronger-wave matrix may then be generated with the command below.  It is a
  sensitivity case, not a replacement for an inconvenient baseline.

```bash
python submit_study_sbatch.py screen --label h015 --prepare --analyze \
  --wave-height 0.15 --nodelist gpu30
```

- Once a production baseline is complete, use one refinement point:

```bash
python submit_study_sbatch.py production --label fine --prepare --analyze \
  --mesh-dir mesh_fine --cell-size 0.01 --particle-spacing 0.005 \
  --nodelist gpu30
```

## Output definitions

`analyze_study.py` writes a JSON summary, probe history, per-cycle pipeline
history and an aggregate CSV.  The manuscript comparison should use only:

1. measured mixture-pore-pressure amplitude ratio and phase lag at the
   crown/shoulder/invert;
2. the upward excess-seepage-force index `IF`, with `IF>=1` as the hydraulic
   critical condition;
3. the actual vertical skeleton effective-stress remaining ratio
   `Rsigma=sigma'_v(t)/sigma'_v0`, with `Rsigma<=0.05` as near-total stress
   loss;
4. representative `q-p'` path and `eps_p_q` (SANISAND) or `pdstrain` (MC);
5. per-cycle pipe `uy/D`, rotation, contacts, breakout time and criterion area,
   area-time and resolved-duration histories;
6. post-release soil `max|u|/h` and `fraction(|u|>=h)` as method-choice
   evidence, reported separately from the deformation/strain claim.

The phase metric uses `PIC_pore_pressures` (the saturation-weighted mixture
pressure used in the parametric manuscript section); liquid and gas pressure
histories remain available separately.  The force index removes the local
hydrostatic component before normalisation:

`IF = max(0, (f_seepage + rho_l*g) dot (-g/|g|) / gamma_sub)`.

The stress criterion is evaluated only where the initial confinement satisfies
`|sigma'_v0|>=100 Pa`, avoiding a near-zero surface denominator.  Both primary
criteria are volume weighted and are reported over the full seabed and the
pre-registered pipeline-support region
`|x-xc|<=1.5D, bed_y-2D<=y<=bed_y`.  `ru` remains in probe histories only
as an auxiliary pore-pressure response; it is not used to label liquefaction,
compute its area/duration, or support the conclusion.  The legacy solver field
`liquefaction_potentials` is an effective-stress-loss quantity, not a
seepage-force index, and is retained only for QA.

## Recommended result presentation

Use three composite figures to make the causal chain visible without another
large parameter sweep:

1. **Trigger chain:** phase/amplitude histories, `IF` and `Rsigma` fields at the
   same critical phase, followed by pipeline contact and `uy/D`.  This shows
   phase lag -> excess gradient force -> realised effective-stress loss ->
   support loss.
2. **Physical and phase controls:** LS/HS and RL/RE curves of support-region
   criterion area, area-time, contact count and per-cycle uplift.  RL/RE must
   include the fitted mean/amplitude QA beside the response comparison.
3. **Constitutive and large-deformation response:** RL/RM matched-pressure
   `q-p'` paths, plastic-history variables, pipeline motion and soil `|u|/h`,
   plus particle snapshots immediately before stress loss, at breakout and
   after migration.

Do not use a single final-displacement bar as the SANISAND argument.  Show its
cyclic path, accumulation and contact evolution; if the large-deformation gate
is not crossed, state that result directly and keep the MPM claim modest.

## Scope and limitations

The rigid pipe is loaded by self-weight, full submerged Archimedes buoyancy and
soil-particle contact.  Direct oscillatory wave pressure/drag on an exposed pipe
is not implemented.  The calculation therefore isolates seabed-mediated
support loss and buoyant breakout; post-breakout migration must be described as
an idealised response, not complete wave–pipe FSI.

The phase-erased replay rotates each particle's fitted fundamental liquid and
gas pressure components to the phase of its local surface point.  It retains
the point mean, fundamental amplitude, residual/higher-harmonic signal and the
progressive phase along the wave direction.  SHA-256 hashes and transformation
parameters are written beside the database.  It intentionally breaks full
two-way conservation and is only a diagnostic counterfactual.

The engineering motivation follows observed liquefaction-driven pipeline
flotation and displacement ([Sumer et al., 1999](https://doi.org/10.1016/S0378-3839(99)00024-1),
[Sumer et al., 2006](https://doi.org/10.1061/(ASCE)0733-950X(2006)132:4(266))).
The constitutive ablation is grounded in the original SANISAND framework
([Dafalias and Manzari, 2004](https://doi.org/10.1061/(ASCE)0733-9399(2004)130:6(622))).
For positioning, MPM has already been used for large-displacement pipe uplift
([Wang et al., 2022](https://doi.org/10.1016/j.tust.2021.104203)); the novelty
claimed here is the phase-resolved three-phase trigger and matched phase
counterfactual, not merely using MPM for a pipe.
