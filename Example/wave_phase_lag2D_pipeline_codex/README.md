# Phase-lag-controlled pipeline instability study

For an exact Linux workstation continuation, including the verified baseline,
mandatory smoke DAG, registered run order and local analysis commands, read
[`LINUX_5080_HANDOFF.md`](LINUX_5080_HANDOFF.md) first.

This directory contains the audited case-study workflow for the manuscript.
The case is a shallow-buried, empty PE100 SDR17 pipeline used as an indicator
of support loss. The registered screen tests, rather than assumes, the proposed
phase-lag mechanism:

> Phase-resolved three-phase MPM testing whether hydraulic phase lag changes
> near-pipeline pressure-gradient forcing, skeleton-stress loss and subsequent
> pipeline migration, with a phase-erased numerical counterfactual and a cyclic
> constitutive-model ablation.

The earlier registered seven-case phase-only screen is retained as a
sensitivity analysis, not as the primary Section 7.5 result. The raw
phase-erased database preserved point means and fundamental amplitudes and
reduced the mean probe lag by 45.37 degrees, but the response-amplitude control
failed. An unsmoothed raw-database derivative gave RL larger local positive-part
activity, whereas the VTK-smoothed derivative gave RE the larger `IF>=1`
area--time; pipe motion was materially indistinguishable. This derivative
sensitivity and the absence of a net-resultant test prevent a lag-to-
liquefaction or uplift claim from that screen. SANISAND's additional cyclic
mechanisms are demonstrated only by the material-point diagnostic;
field-scale superiority is not claimed because MC_EQ failed the unchanged
stability gate.

The primary current-runner permeability screen uses four fixed points at
`Sw=0.94`, disables both solver pressure-smoothing switches, and evaluates raw
pressure-database gradients over two complete post-ramp cycles. All four points
fail the net-uplift engineering gate. The fresh `k=1e-13 m2` point resolves a
`+49.90 deg` crown phase difference with a crown/surface amplitude ratio of
`0.132`, yet its net-uplift impulse changes by `-0.0239%`; the other three net
changes are `-0.0131%`, `-0.1399%` and `-0.9038%`. The result supports the
presence of phase lag, but not a resolved amplification of the pipeline-zone
net-uplift driver in this matrix.

A separate phase-conditioned diagnostic now answers the narrower engineering
question that the complete-cycle integral cancels: short upward drag while the
same-column seabed pressure is negative and recovering from its trough. At
`k=3e-12 m2`, the two repeated recovery branches increase raw local
upward-force activity by `6.77%` and `6.67%`; their signed support-force
impulses increase by `1.27%` and `1.16%`. The other permeability points do not
pass the fixed short-time mechanism gate. Thus the data support a transient
hydraulic-risk window, not a monotonic low-permeability rule and not realised
liquefaction. Reproduce the audited CSV and common-time two-dimensional fields
with:

```bash
python3 analyze_phase_conditioned_uplift.py
```

Numeric derivatives and metrics remain unsmoothed; `sigma=1.25 pixels` is used
only for the displayed two-dimensional raster. The resulting HTML technical
report and PNG/PDF figures are written under
`analysis/phase_conditioned_uplift/`.

The code does **not** assume that SANISAND must predict more displacement than
Mohr-Coulomb.  The comparison tests whether cyclic memory, state dependence,
fabric evolution and cyclic mobility change the trigger, accumulation rate and
post-trigger path.  A result is useful even when final displacement ordering is
not monotonic.

## Study matrix

| Code | Phase input | Soil model | Status and purpose |
|---|---|---|---|
| MC_EQ | No wave; resumes stable EQ_HS | Mohr-Coulomb | Pure-PIC, damped, fixed-pipe handoff relaxation; numerical staging only |
| LS | Fully coupled, near-saturated `Sw=0.993` | SANISAND | Physical low-lag reference; report the measured lag, never call it exact zero |
| HS | Fully coupled, `Sw=0.94` | SANISAND | Physical high-lag main case |
| HM | Fully coupled, `Sw=0.94` | Mohr-Coulomb | Simple constitutive ablation against HS |
| HD | Fully coupled, fixed pipe, `Sw=0.94` | SANISAND | Generates the pressure database |
| RL | Original HD pressure history, free pipe | SANISAND | One-way replay control |
| RM | Original HD pressure history, free pipe | Mohr-Coulomb | Matched-pressure constitutive ablation against RL |
| RE | Same mean/amplitude as RL, fundamental lag erased | SANISAND | One-way numerical counterfactual; **not** a physical fully coupled solution |

LS–HS is the physical comparison, but saturation also changes attenuation and
drainage.  RL–RE is therefore a controlled one-way phase-component diagnostic.
RL–RM applies
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
and calibration values.  The separate `critical_timestep_modulus=23.8 MPa` is
only a conservative wave-speed bound for the explicit time-step audit; it is
not used by the SANISAND stress update.

## Local Linux workflow

The same fail-fast validation and completion audit can be used without Slurm,
environment modules or Conda. Build `mpm_hpc_source/build-pipeline/mpm`, then
run one registered case from this directory:

```bash
MPM_THREADS=8 bash run_case_local.sh configs/screen/01_EQ_LS.json
```

`run_case_local.sh` validates runtime dependencies, removes only that case's
old completion sentinel, executes the solver, and publishes a new sentinel only
after the complete configured particle-VTP grid, unique pipeline-history CSV,
same-step HDF5/VTP cross-check and any unchanged equilibrium stability gate pass.
Run the four required shortened stages before the registered 40,000-step
equilibria or 104,000-step screen:

```bash
MPM_THREADS=4 bash run_local_smoke.sh new-lowercase-label
```

The default fail-fast DAG contains only `EQ_LS`, `EQ_HS`, `SANI_LS` and
`SANI_HS`.  `MC_EQ` is a separate, explicit diagnostic:

```bash
MPM_THREADS=4 bash run_local_smoke.sh --mc-diagnostic new-diagnostic-label
```

The diagnostic uses the same unchanged `max|v| <= 1e-3 m/s` gate and exits
nonzero if it fails.  It never launches an MC dynamic smoke case.  Neither the
required smoke nor the optional MC diagnostic is manuscript evidence; formal
analysis rejects their validation profile even when every completion audit
passes.  The v2 smoke manifest labels every case `required` or `diagnostic`
and hash-binds the canonical payload, complete config file and manifest; the
wrapper revalidates those hashes before launching any stage.

The post-facet-correction MC_EQ diagnostic reached
`max|v|=1.8469e-3 m/s`, so it remains a failed diagnostic rather than evidence
or a reason to relax the gate.  After separating the assigned top-boundary
phase kinematics from geometric cavity free surfaces, the current binary
passed all four required stages under the fresh local label
`finalfs3-20260813`.  EQ_LS/HS remained far below the unchanged velocity gate,
the audited quiet top rows stayed in compression, and SANI_LS/HS retained all
7,376 particles inside the mesh.  These ignored results certify the local
smoke only; a new machine must use a new label.  The optional local MC result
does not certify the registered MC branch: formal HM/RM work still depends on
a separately gated registered `MC_EQ` checkpoint.

After a fresh four-stage smoke passes, the guarded local formal controller runs
only the seven registered phase-study cases (no MC branch), with at most two
four-thread solvers and exact completion checks after every case:

```bash
MPM_THREADS=4 PHASE_SCREEN_MAX_TOTAL_THREADS=8 \
  bash run_local_phase_screen.sh
```

It follows `EQ_LS/EQ_HS -> LS/HS -> HD -> phase transform -> RL/RE`. A
non-empty result directory without a valid completion is a hard stop for manual
audit; the controller never deletes or overwrites partial formal output.

## Server workflow

Run from this directory after checking out the branch.  Windows does not need
to build the C++ code.  The top-level Python runner is lightweight; all input
generation, validation, compilation and simulation work is launched inside an
HPC4 allocation.

```bash
# Submit the registered phase-only screen dependency chain and return immediately.
python submit_study_sbatch.py screen --stage phase --prepare --build --analyze \
  --cpus 16 --gres gpu:1 --nodelist gpu30
```

The two primary equilibrium jobs deliberately use `LinearElastic2D`, pure
`PIC=1.0` with `APIC=false`, a fixed pipe and Cundall damping of `5 s^-1`
for 4 s.  The particle update keeps this damping active in the pure-PIC branch;
dynamic cases use zero damping.
The equilibrium elastic constants use the same `E=9.10081 MPa`,
`nu=-0.0079625` reference tangent as the registered Mohr-Coulomb ablation;
the plane-strain initial state sets both horizontal effective stresses to
`K0 sigma_yy`.  The separate `23.8 MPa` timestep modulus remains a conservative
numerical bound.
The low-lag case uses `Sw=0.993`, avoiding the singular near-saturated endpoint
while retaining the previously demonstrated low-lag regime.  `mesh.py` also generates a
hydrostatic liquid-pressure profile for every particle plus separate LS/HS
effective self-weight stress fields.  The submerged bed surface carries the
matching downward `4.905 kPa` water-column traction. Both the prescribed phase
pressure and the matching total-pressure traction are assigned only to the
single physical top particle row; the two-row near-surface set is diagnostic
only. Wave-pressure increments receive the same total normal traction on
that physical top row; this is a seabed boundary load, not direct wave drag on
the pipe. The phase momentum equations consistently use pressure increments
relative to the equilibrium checkpoint for both liquid and gas, so the
initial gas suction/pressure cannot create an artificial unbalanced force.
Validation checks the single-column MPM
pressure format and the stress/pressure values before a job is submitted.
SANISAND dynamics restore the successful LinearElastic2D HDF5 directly.  The
cross-model handoff derives the SANISAND void ratio from restored porosity and
centres its finite `m_iso` yield surface on the restored effective stress.
Every stage also suppresses velocity updates only on emerging nodes carrying
less than 1% of a full-cell skeleton support; the physical equilibrium gate
remains the unchanged `max|v| <= 1e-3 m/s`.
zero-cohesion MC comparison instead passes through the explicit `MC_EQ` stage:
it restores EQ_HS, switches to Mohr-Coulomb, retains pure PIC, damping
`5 s^-1`, the fixed pipe and zero wave loading for another 4 s, and writes a
new checkpoint.  This stage is numerical handoff relaxation, not a scientific
comparison result.  HM/RM restore that MC checkpoint and only then return to
`APIC=true`, `PIC=0`, zero damping and their registered dynamic loading.

The submitter uses `sbatch --parsable` and records the dependency graph and job
IDs under `submissions/`.  Every job is forced onto
`--partition=granularmech --account=comgranmech`; the defaults also request
`gpu30`, 16 CPU and `gres:gpu:1`. With `--stage phase`, it submits the two
primary equilibria, LS/HS, the fixed-pipe driver, phase transform, RL/RE and
phase analysis with `afterok` dependencies, then returns so the terminal may
be disconnected or reused for Code Tunnel. The separate `--stage full` mode
also includes MC_EQ/HM/RM, but must not be used unless MC_EQ first passes the
unchanged registered gate. Completed outputs are skipped only when the complete particle-VTP
grid, pipeline-history CSV, HDF5/VTP cross-check and dependencies remain
hash-valid, unless `--force` is
specified.  A file is not considered complete merely because it exists:
`run_case.sh` publishes `pipeline_completion.json` atomically only after the
solver exits successfully and the complete VTP grid, unique history CSV,
required HDF5/VTP cross-check, configuration hash
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
checkpoint exists and that its final QA VTP has finite fields,
`max|v| <= 1e-3 m/s`, displacement below one coarse cell, and valid porosity.
Both LinearElastic2D equilibrium jobs and MC_EQ apply this unchanged gate while
publishing their own completion sentinel; failure makes the job fail and stops
all dependent `afterok` jobs.  The server build also rejects a legacy
159-field SANISAND checkpoint that declares more history variables than it
actually stores; the registered workflow always creates a fresh linear-elastic
equilibrium checkpoint before either constitutive handoff.
Each batch job writes `logs/<job-name>-<job-id>.out/.err`.
The submit command returns after queueing; monitor later with `squeue -u "$USER"`
or inspect the saved `submissions/<group>-<UTC timestamp>.json` while using Code
Tunnel for other work.

To resubmit only the reduction after jobs finish:

```bash
sbatch --partition=granularmech --account=comgranmech \
  --nodes=1 --ntasks=1 --cpus-per-task=16 --gres=gpu:1 \
  --nodelist=gpu30 --time=02:00:00 --job-name=plp_screen_analyze \
  run_task.sbatch bash analyze_results.sh screen phase
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
history and an aggregate CSV.  For a full analysis, `analyze_results.sh` also
runs `synthesize_manuscript_evidence.py` and `plot_manuscript_figures.py`.  The
former writes `manuscript_evidence.json/.md` with pre-registered claim gates;
the latter writes PNG/PDF curve composites and common-scale two-dimensional
fields (when Matplotlib is available), plus JSON manifests selecting common
registered RL/RE times.
The Section 7.5 draft and the
source-PDF correction audit are in `MANUSCRIPT_SECTION_7_5_DRAFT.md` and
`MANUSCRIPT_PDF_REVISION_NOTES.md`.  The manuscript comparison should use only:

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

For an explicit two-dimensional liquid-phase pressure and phase-lag diagnostic
of the completed physical HS screen, run:

```bash
python3 plot_liquid_pressure_phase_lag.py
```

This subtracts the audited HS_EQ `PIC_liquid_pressures` checkpoint by particle
ID, draws four common-scale quarter-cycle fields, fits the fundamental over the
registered 1.3--3.9 s pre-release window, and writes same-column-surface phase
differences plus a time--depth section under
`analysis/screen/liquid_pressure_phase_lag/`. The same-column reference isolates
subsurface lag; a comparison to the central surface also contains the
horizontal travelling-wave phase.

### Exploratory low-permeability, no-smoothing replay

The current four-point, fail-closed permeability screen and its exact gates are
recorded in `PHASE_LAG_PERMEABILITY_SCREEN.md`. Use
`run_phase_lag_exploration.py` for new calculations; the longer direct-solver
sequence below documents the original `1e-13 m2` legacy run and must not be
used to overwrite or resume a partial exploration namespace.

The registered screen/production configurations remain immutable.  A separate
exploratory namespace can test a more strongly diffusion-limited pressure
response using `Sw=0.94`, `intrinsic_permeability=1e-13 m^2`, and both pressure-
smoothing switches disabled.  Here `1e-13` is the solver's intrinsic-
permeability input in square metres; it is not hydraulic conductivity in m/s.
The equilibrium remains LinearElastic2D, pure PIC, fixed-pipe, no-wave and
`5 s^-1` damping, with the unchanged `max|v|<=1e-3 m/s` gate.  The dynamic
driver and replays use SANISAND, APIC, zero damping and a fixed pipe.

```bash
python3 prepare_phase_lag_exploration.py \
  --label k1e-13_sw094_nosmooth \
  --permeability 1e-13 --saturation 0.94 \
  --no-pressure-smoothing --paired-replay

../../mpm_hpc_source/build-pipeline/mpm -p 3 -f ./ \
  -i configs/phase_lag_exploratory/k1e-13_sw094_nosmooth/01_EQ.json
../../mpm_hpc_source/build-pipeline/mpm -p 3 -f ./ \
  -i configs/phase_lag_exploratory/k1e-13_sw094_nosmooth/03_HD.json

python3 phase_controls.py \
  --source-dir pressure_databases/phase_lag_exploratory/k1e-13_sw094_nosmooth/lagged \
  --source-prefix pressure \
  --output-dir pressure_databases/phase_lag_exploratory/k1e-13_sw094_nosmooth/phase_erased \
  --output-prefix phase_erased --period 1.3 \
  --fit-start 2.6 --fit-end 3.9 --ramp-time 1.3

../../mpm_hpc_source/build-pipeline/mpm -p 3 -f ./ \
  -i configs/phase_lag_exploratory/k1e-13_sw094_nosmooth/04_RL.json
../../mpm_hpc_source/build-pipeline/mpm -p 3 -f ./ \
  -i configs/phase_lag_exploratory/k1e-13_sw094_nosmooth/05_RE.json
python3 analyze_phase_lag_replay_exploration.py \
  --label k1e-13_sw094_nosmooth
```

`03_HD` writes the fully coupled pressure history; `04_RL` and `05_RE` read
the original and locally phase-erased versions of that exact history.  The
transform retains each particle's mean, amplitude, residual/higher harmonics
and along-wave progressive phase.  Consequently no-smoothing spatial texture
is common to both replays and cannot by itself be counted as a phase-lag
effect.  The analysis requires the exact 60-frame grid, finite/unique particle
IDs, no particles outside the mesh, valid porosity and the primary `IF`,
`R_sigma` and same-particle joint criteria.  These outputs are explicitly
exploratory and do not replace registered manuscript evidence.

The legacy WAVE-v2 `k=1e-13 m^2`, `Sw=0.94`, no-smoothing run resolves a
credible local pressure lag at the pipe crown: `+50.20 deg` (the later HD
phase-v3 reconstruction gives `+50.02 deg`), amplitude `37.48 Pa`
(`0.132` of the same-column surface amplitude), and harmonic `R^2=0.9836`.
The shoulder and far-field same-column lags are `-10.20 deg` and `+7.23 deg`;
the invert estimate is rejected by the predeclared harmonic-quality gate
(`R^2=0.5526`).  Removing phase while retaining fitted mean and amplitude gives
a maximum framewise particle `q95` pressure difference of `29.37 Pa` at
`t=3.705 s`.  Both 39000-step replays retain all 7376 particles in the mesh;
the lagged run has `max|v|=1.370e-4 m/s` and porosity
`0.484404--0.485223`.

This stronger local lag does **not support the claim** that retained lag makes
liquefaction easier in the short, fixed-pipe exploratory replay. Relative to
phase-erased, lagged `IF>=1`
area-time is `0.35%` smaller; the sparse `R_sigma<=0.05` area-time is `8.11%`
smaller, and the same-particle/same-time joint criterion is zero in both cases.
The result therefore supports the existence of a local hydraulic phase lag,
but not the manuscript's proposed lag-to-liquefaction causal claim. Its
historical equilibrium-config identity was not frozen by the current-runner
schema, so it is excluded from the primary four-point trend. The large
raw no-smoothing snapshot roughness (`q95=137 Pa` between adjacent saved
snapshots) is reported separately and is not interpreted as phase lag.

The evidence synthesis fixes its thresholds before reading the field result:
RL/RE fitted pressure means and amplitudes must agree within 2%, removal must
reduce the mean crown/shoulder/invert absolute phase lag by at least 5 degrees,
and a response difference must exceed 5% to be called material.  RL/RM also
must pass a 5% first-frame `p'-q` initial-state audit.  This last check matters
because RM passes through MC_EQ while RL resumes the elastic checkpoint
directly; a failed audit downgrades the result to a comparison of complete
handoff/model chains rather than a clean constitutive isolation.

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

Use the completed figures to test the proposed chain without implying that the
chain passed:

1. **Trigger-chain audit:** phase/amplitude histories, `IF`, `Rsigma` and the
   same-particle joint indicator at common registered times, followed by
   pipeline contact and `uy/D`. In the earlier registered seven-case screen,
   the separate metrics do not establish a lag-driven causal chain.
2. **Physical and phase controls:** LS/HS and RL/RE curves of support-region
   criterion area, area-time, contact count and per-cycle uplift.  RL/RE must
   include the fitted mean/amplitude QA beside the response comparison.
3. **Constitutive mechanism and handoff audit:** identical material-point
   cyclic paths for SANISAND and Mohr--Coulomb, together with the failed MC_EQ
   field handoff. Do not present an RL/RM field response because no admissible
   RM checkpoint exists.

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
