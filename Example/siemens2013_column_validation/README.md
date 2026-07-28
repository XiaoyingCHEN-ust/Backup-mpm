# Siemens et al. (2013) pore-air infiltration validation

This case validates `ThermoMPMExplicitThreePhaseNew2D` against the paired
open/closed coarse-sand column tests reported by:

G. A. Siemens, S. B. Peters, and W. A. Take (2013), “Comparison of confined
and unconfined infiltration in transparent porous media,” *Water Resources
Research*, 49, 851–863. <https://doi.org/10.1002/wrcr.20101>

## Why this case was selected

The experiment supplies physical measurements of the variables needed to
test the gas phase directly: six pore-pressure transducers, saturation
profiles at 5 s intervals, and wetting-front trajectories. It also provides
two boundary-condition regimes in the same apparatus:

- `mpm_open.json`: pore air is maintained at atmospheric pressure and the
  cropped, initially unsaturated interval has a no-flow liquid boundary at
  its lower wet-interface datum.
- `mpm_closed.json`: the sides and base are impermeable; compressed air can
  escape only by counterflow through the wetted zone and the surface pond.

The closed test measured an approximately uniform 1.2–1.6 kPa pore-air
pressure during counterflow, a transmission-zone saturation of about 0.78,
and a wetting front roughly three times slower than the open test. These are
independent validation targets, not calibration inputs.

## Model idealisation

The computational domain represents the approximately 1.075 m initially
unsaturated interval from the surface of the coarse sand to the lower wet
boundary. Because the experiment and validation observables are one-dimensional,
the numerical domain is a representative 11.4 mm-wide plane-strain slice with
one cell across its width and 94 cells along its height. Standard 2 x 2 particle
integration gives two particles per horizontal layer, 188 layers, and a 5.7 mm
vertical particle spacing. The modeled height is 1.0716 m. Side constraints
enforce one-dimensional flow, and the skeleton is fixed to isolate hydraulic
liquid–gas coupling. The reduced width changes total modeled volume but not
pressure, saturation, or wetting-front quantities defined per unit area.

Published inputs used directly are:

- porosity: 0.50 (reported range 0.48–0.51)
- oil density: 845 kg/m³
- oil viscosity at 25 °C: 0.0101 Pa s
- oil saturated hydraulic conductivity: 1.0e-3 m/s
- derived intrinsic permeability: 1.2188e-9 m²
- van Genuchten parameters: `p0 = 170 Pa`, `m = 0.84`
- residual liquid and gas saturations: 0.03 and 0.02
- surface ponding head: 0.12 m (`994.734 Pa`)

The initial saturation is set just above the reported residual value
(`Sw = 0.0301`) to keep the retention relation finite. Oil compressibility,
air properties, and the low elastic modulus of the fully constrained
skeleton are numerical/supporting properties not reported by the experiment;
they are not fitted to the response curves.

A relative-permeability floor of `1e-6` is applied only when a phase reaches
its residual endpoint. Without this regularisation the model evaluates
viscosity divided by zero permeability at the saturated surface (gas phase)
and in the residual-water zone (liquid phase), which makes the nodal drag
matrix non-finite. The floor keeps the residual phase effectively locked and
is not a fitted infiltration parameter.

## Files

- `preprocess.py` generates the mesh, particles, volumes, temperatures, and
  entity sets.
- `mpm_open.json` and `mpm_closed.json` are the two solver inputs.
- `reference_data/` contains values read from Figures 7, 8, and 11 and the
  pressure range stated in the paper. The wetting-front points are approximate
  digitisation targets and should not be presented as exact tabulated data.
- `run_stability_hpc4.sh` is a manifest-driven, short-duration time-step check.
- `prepare_checkpoint_validation.py`, `run_checkpoint_segment_hpc4.sh`, and
  `submit_checkpoint_validation_hpc4.sh` generate and submit the dependent r20
  restart-continuity test.
- `compare_checkpoint_restart.py` compares every reported hydraulic field in
  the resumed final state with the continuous 3 s reference.
- `compare_permeability_time_scaling.py` tests a proposed permeability/time
  acceleration after normalising phase permeability and velocity fields.
- `prepare_production_segments.py`, `run_production_segment_hpc4.sh`, and
  `submit_production_hpc4.py` generate, check, and submit the restartable
  full-duration dependency chains.
- `run_hpc4.sh` is retained only as a fail-safe: the direct full-duration
  array is disabled because neither case fits the 24-hour job limit.
- `postprocess_validation.py` extracts wetting-front depth and pore-air
  pressure at the six experimental PPT elevations from the particle VTP files.
  A sensor curve terminates when the wetting front arrives and the instrument
  would begin measuring pore-oil pressure.
- `program_patch/` installs the backward-compatible source replacement,
  rebuilds the executable, and prevents a stale executable being used again.
- `validation_subsection_draft.md` is intentionally left with placeholders
  until the HPC results are available.

## Generate and check the case

From this directory:

```bash
python preprocess.py
```

Expected counts are 94 cells, 376 particles, four top-boundary particles (the
complete top cell), and two particles marking the lower datum. The latter are
retained in the generated entity sets for diagnostics but are not assigned a
liquid-pressure constraint in the cropped-domain input.

## Rebuild and run the mandatory short check on HPC4

The downloaded first run is not valid for comparison. Its initial dry-zone
suction is 137340 Pa instead of the input-consistent 972.99 Pa, proving that
the executable was not rebuilt with the validation patch. The pressures then
grow far beyond the experimental scale. Install the two replacement files in
the source tree used by the executable and rebuild:

```bash
sbatch program_patch/rebuild_hpc4.sh
```

Wait until this one rebuild job completes successfully. It requests the same
`granularmech`/`comgranmech` resources and 32 CPUs as the calculation jobs.

The default first stage is deliberately small: six 0.01 s checks covering the
open and closed cases at `5e-5`, `5e-6`, and `1e-6` s. They require only 200,
2000, and 10000 steps per case:

```bash
python prepare_stability_checks.py
sbatch run_stability_hpc4.sh
```

The first-stage UUIDs contain `r3-smoke`. The completed r3 runs remained finite
but failed the time-step consistency check. At the common time 0.001 s, reducing
the step caused larger rather than converging pressure and saturation changes.
The source input also smooths liquid pressure every 100 computational steps, so
the three runs apply this operation 2, 20, and 100 times over 0.01 s. This is a
time-step-dependent algorithm and cannot be used to establish convergence.

Run the next, still smaller, diagnostic with pressure smoothing disabled and
PIC phase-velocity transfer. These settings leave all physical properties and
boundary conditions unchanged. The duration is only 0.001 s, corresponding to
20, 200, and 1000 steps per case:

```bash
python prepare_stability_checks.py --duration 0.001 \
  --revision r4-pic-micro --pic 1 --pressure-smoothing false
sbatch run_stability_hpc4.sh
```

The r4 results for the original four-cell-wide mesh are finite and time-step
convergent: between `5e-6` and `1e-6 s`, the maximum open-case liquid-pressure
and saturation differences are 2.75 Pa and `1.43e-6`; the corresponding
closed-case differences are 0.97 Pa and `5.37e-7`, with a maximum gas-pressure
difference of 0.79 Pa. However, r4 is only a diagnostic because pure PIC
transfer introduces numerical dissipation. It is not adopted for the
experimental validation.

Isolate the effect of step-count-based pressure smoothing with a second 0.001 s
micro test. This restores both transfer ratios to zero and changes only
`pressure_smoothing` relative to the original numerical settings:

```bash
python prepare_stability_checks.py --duration 0.001 \
  --revision r5-flip-micro-1d --pic 0 --pic-t 0 --pressure-smoothing false
sbatch run_stability_hpc4.sh
```

The r5 reduced-mesh test confirms FLIP convergence when step-count-based
pressure smoothing is disabled. At 0.001 s, reducing the step from `5e-5` to
`5e-6 s` and then `1e-6 s` reduces the closed-case gas-pressure RMS field
difference from 3.06 Pa to 0.177 Pa. The medium-step gas-pressure field differs
from the fine reference by at most 2.41 Pa. The `5e-5 s` step is not selected
at this stage.

The imposed surface boundary now covers all four particles in the complete top
cell. Regenerate the geometry and repeat the 0.001 s FLIP check under a new UUID
before extending physical time:

```bash
python preprocess.py
python prepare_stability_checks.py --duration 0.001 \
  --revision r6-flip-micro-1d-topcell --pic 0 --pic-t 0 \
  --pressure-smoothing false
sbatch run_stability_hpc4.sh
```

The r6 gas-pressure field converges strongly. In the closed case, the `5e-6 s`
field differs from the `1e-6 s` reference by at most 0.056 Pa and by 0.0207 Pa
in RMS; the coarse-to-medium versus medium-to-fine RMS difference decreases by
a factor of 24. The `5e-5 s` solution is again too coarse. Imposing pressure on
the complete top cell also produces a convergent but still-growing local liquid
velocity: at 0.001 s it is 0.535 m/s in the top row and 0.346 m/s immediately
below the boundary cell. Confirm that this startup transient decays or remains
bounded over 0.01 s before selecting a production step:

```bash
python prepare_stability_checks.py --duration 0.01 \
  --revision r7-flip-short-1d-topcell --time-steps 5e-6 1e-6 \
  --pic 0 --pic-t 0 --pressure-smoothing false
sbatch --array=0-3 run_stability_hpc4.sh
```

The r7 transient is bounded and time-step consistent. The maximum liquid
velocity occurs near 0.002 s at 0.596 m/s and then decreases; at 0.01 s the
top-row and first-interior-layer values are 0.483 and 0.309 m/s. The two steps'
final interior liquid-velocity components differ by at most 0.00126 m/s. In the
closed case, the `5e-6 s` final gas-pressure field differs from the `1e-6 s`
reference by at most 0.863 Pa and by 0.368 Pa in RMS.

Although convergent, `5e-6 s` would require approximately 175 million steps for
an 875 s closed test. Use the existing r7 fine run as the 0.01 s reference and
screen practical larger steps with only 20, 100, and 200 steps per case:

```bash
python prepare_stability_checks.py --duration 0.01 \
  --revision r8-large-dt-micro --time-steps 5e-4 1e-4 5e-5 \
  --pic 0 --pic-t 0 --pressure-smoothing false
sbatch run_stability_hpc4.sh
```

The r8 test rejects both `5e-4` and `1e-4 s`. They exceed the 1 MPa safety limit
within approximately 0.0025 and 0.002 s in the open case and diverge even more
rapidly in the closed case, reaching final pressures of order `1e23--1e30 Pa`.
The `5e-5 s` runs remain finite. Their liquid velocity peaks near 0.002 s at
about 0.600 m/s and decreases to 0.496 m/s by 0.01 s. The closed-case interior
gas-pressure range at 0.01 s is 3.2--10.1 Pa, compared with approximately
-0.7--4.8 Pa in the previous fine reference. This absolute difference is small,
but the coarse solution has an early liquid-pressure overshoot, so extend the
comparison by one order of magnitude before selecting it:

```bash
python prepare_stability_checks.py --duration 0.1 \
  --revision r9-flip-0p1s --time-steps 5e-5 5e-6 \
  --pic 0 --pic-t 0 --pressure-smoothing false
sbatch --array=0-3 run_stability_hpc4.sh
```

The r9 test rejects `5e-5 s` as a common production step. The open case diverges
between 0.05 and 0.07 s and reaches `7.3e16 Pa`. The closed case remains finite,
but at 0.1 s its interior gas pressure is 39.8--44.7 Pa, versus 9.7--18.1 Pa
with `5e-6 s`; the pointwise gas-pressure field differs by 31.8 Pa maximum and
27.9 Pa RMS. Both `5e-6 s` cases remain stable. Test intermediate steps in one
batch so that every comparison field is retained together:

```bash
python prepare_stability_checks.py --duration 0.1 \
  --revision r10-intermediate-dt --time-steps 2e-5 1e-5 5e-6 \
  --pic 0 --pic-t 0 --pressure-smoothing false
sbatch run_stability_hpc4.sh
```

All six r10 runs remain finite. Relative to the `5e-6 s` field at 0.1 s, the
closed-case gas-pressure maximum/RMS differences are 3.59/3.14 Pa for `1e-5 s`
and 10.73/9.40 Pa for `2e-5 s`. The open case is more time-step sensitive: its
saturation and liquid-pressure RMS differences are 0.00345 and 39.7 Pa for
`1e-5 s`, versus 0.0335 and 111.9 Pa for `2e-5 s`. The closed gas-pressure
errors approximately halve as the step halves, consistent with first-order
convergence. Treat `1e-5 s` as the accuracy candidate and retain `2e-5 s` only
as a performance candidate. Extend all three steps to 3 s before choosing:

```bash
python prepare_stability_checks.py --duration 3 \
  --revision r11-convergence-3s --time-steps 2e-5 1e-5 5e-6 \
  --pic 0 --pic-t 0 --pressure-smoothing false
sbatch run_stability_hpc4.sh
```

This creates six tasks with 150000, 300000, or 600000 steps. Do not proceed to
the experimental pressure-rise interval until these results have been reviewed.

The r11 range checks initially passed, but a direct comparison of the two
particles in every horizontal layer exposed a nonphysical lateral mode in the
open case. The maximum within-layer saturation difference reaches 0.95. It is
already visible at the first 0.15 s output for `dt=2e-5 s`, at 0.30 s for
`dt=1e-5 s`, and between 1.95 and 2.10 s for `dt=5e-6 s`. Consequently, none
of the open r11 runs is accepted and the production step is not yet selected.
The closed runs retain layer symmetry to roundoff, remain bounded, and show a
continued decrease in the gas-pressure difference as the step is reduced.

Siemens et al. represent their open test using constant heads at the top and
base plus manometer ports open along the column. In this reduced model,
`fixed_gas_pressure=true` already homogenizes that distributed venting as zero
gauge gas pressure. The side-wall gas velocity is therefore now constrained in
the horizontal direction, as is required by the one-dimensional equivalent
model, to avoid representing lateral venting twice. Test this correction only
with the previously most stable step before repeating a convergence study:

```bash
python prepare_stability_checks.py --duration 3 \
  --revision r12-open-gas-x --time-steps 5e-6 \
  --pic 0 --pic-t 0 --pressure-smoothing false
sbatch --array=0 run_stability_hpc4.sh
```

This creates open and closed manifest rows, but `--array=0` submits only the
open diagnostic (600000 steps). The job script reads the generated manifest,
so its array index must match the desired row.

The r12 open diagnostic passes the strengthened check. Across all 21 outputs,
the largest saturation difference between the two particles in a horizontal
layer is `2.71e-14`, compared with 0.95 in r11. Horizontal liquid and gas
velocities remain zero. The maximum absolute phase pressure is 1.260 kPa, and
the liquid-velocity magnitude peaks at 0.126 m/s at 0.15 s before decreasing
to 0.0526 m/s at 3 s. The central dry interval retains its initial saturation;
the saturation increase below it is the separate base-connected wet region
created by the open experiment's lower constant-head boundary. The top-
connected front has not moved beyond the prescribed top cell by 3 s.

Complete the corrected open-case time-step comparison without rerunning the
existing r12 `5e-6 s` reference:

```bash
python prepare_stability_checks.py --duration 3 \
  --revision r13-open-convergence --time-steps 2e-5 1e-5 \
  --pic 0 --pic-t 0 --pressure-smoothing false
sbatch --array=0-1 run_stability_hpc4.sh
```

This submits only the two open rows, with 150000 and 300000 steps. Compare
their 21 common output times with the r12 reference before choosing a longer
physical-time run.

Both r13 runs are finite and retain horizontal symmetry. At 3 s, the lower
connected wet-region depths are 0.2280, 0.2508, and 0.2622 m for time steps
`2e-5`, `1e-5`, and `5e-6 s`, respectively. Relative to `5e-6 s`, the
`1e-5 s` mean saturation differs by 0.00367 and its lower-front position by
one cell (0.0114 m); the interior liquid-pressure RMS difference is 76.6 Pa.
The `2e-5 s` errors are appreciably larger and that step is rejected.

The time-step trend also exposed a boundary-model error that must be corrected
before extending the run. Figure 5 of Siemens et al. shows only a thin wet
layer at the lower constant-head screen at both 0 and 45 s, whereas every r12/
r13 solution develops a 0.23--0.26 m base-connected wet region by only 3 s.
The validation domain represents the initially unsaturated interval above the
lower fluid surface and the simulation terminates when the top-connected front
first reaches that surface. Before contact, impose zero vertical liquid
velocity at the lower boundary rather than forcing the bottom particles to
zero suction. Distributed open-system air venting remains represented by
`fixed_gas_pressure=true`.

First test this lower-boundary correction over only 0.6 s, which is long enough
to cover the onset of the spurious lower wet region seen in r12/r13:

```bash
python prepare_stability_checks.py --duration 0.6 \
  --revision r14-open-dry-base --time-steps 1e-5 5e-6 \
  --pic 0 --pic-t 0 --pressure-smoothing false
sbatch --array=0-1 run_stability_hpc4.sh
```

The two open tasks contain 60000 and 120000 steps. Do not extend to the
experimental time scale until the lower saturation profile has been checked.

The r14 dry-base runs are finite, retain layer symmetry to roundoff, and do
not develop a lower connected wet region. At 0.6 s, the `1e-5` versus `5e-6 s`
maximum/RMS interior saturation differences are `3.85e-4` and `2.86e-5`;
the corresponding liquid-pressure differences are 194 and 21.9 Pa. However,
the entire upper half of each r14 field is identical to r12 at the common
times to approximately `1e-12 Pa` and `1e-15 m/s`. Thus r14 used the previous
liquid-force implementation and tests only the lower-boundary change.

Code review then identified a separate two-dimensional liquid momentum error.
Liquid gravity is already mapped as an external body force, while the internal
force additionally subtracted a hard-coded hydrostatic pressure. Because
`PIC_liquid_pressure_` stores actual gauge pressure rather than excess pressure,
this introduced an extra gravity gradient. The two-dimensional internal force
now uses `PIC_liquid_pressure_` directly, consistent with the three-dimensional
implementation. Restore the experimental lower constant-head condition and
repeat the same 0.6 s matrix after installing and rebuilding the updated patch:

```bash
sbatch program_patch/rebuild_hpc4.sh

python prepare_stability_checks.py --duration 0.6 \
  --revision r15-force-corrected --time-steps 1e-5 5e-6 \
  --pic 0 --pic-t 0 --pressure-smoothing false
sbatch --array=0-1 run_stability_hpc4.sh
```

Do not reuse r14 as the fine reference: r15 changes the governing liquid
momentum force. The run scripts now reject source trees that still contain the
obsolete hard-coded hydrostatic subtraction.

Both r15 jobs produced all 21 requested snapshots. The `5e-6 s` run remained
bounded and laterally symmetric, with maximum absolute pressure `1.253 kPa`,
maximum liquid-velocity component `0.06343 m/s`, and within-layer saturation
spread `1.15e-14`. The `1e-5 s` run lost lateral symmetry at 0.6 s: its maximum
within-layer saturation spread reached `0.00543`, and its mean saturation
dropped nonphysically from `0.1140` at 0.54 s to `0.0894` at 0.6 s.

The fine r15 field also exposed a boundary-initialisation incompatibility that
the previous range checks did not detect. At 0.6 s, its lowest layer is
saturated, the next layer remains at `Sw=0.0354`, and a separate nearly
saturated band occupies approximately 0.014--0.100 m. This disconnected band
first appears by 0.06 s. Figure 5 of Siemens et al. already contains a lower
wet zone at the experimental initial time, whereas the reduced model contains
only the initially unsaturated interval above that datum. Suddenly imposing
zero suction on its lowest particle row therefore creates an artificial
pressure/saturation transient.

Retain the corrected liquid force but return the cropped lower boundary to
zero vertical liquid velocity. Repeat the 0.6 s comparison before any longer
run:

```bash
python prepare_stability_checks.py --duration 0.6 \
  --revision r16-force-corrected-dry-base --time-steps 1e-5 5e-6 \
  --pic 0 --pic-t 0 --pressure-smoothing false
sbatch --array=0-1 run_stability_hpc4.sh
```

The range checker now also rejects any layer with mean saturation at or above
`0.40` that is not connected continuously to either the upper or lower model
boundary. This prevents a bounded but disconnected wet band such as r15 from
being accepted. A strict experimental constant-head base would instead require
initialising the measured lower wet profile; applying the pressure to more dry
particles would not supply that missing initial state.

Both r16 open runs pass the symmetry, range, and wet-region connectivity
checks. Their top-connected wet region remains one cell deep and the lower
boundary stays dry. Relative to `5e-6 s`, the `1e-5 s` saturation maximum/RMS
differences at 0.6 s are `3.85e-4` and `4.00e-5`; the liquid-pressure maximum
difference is 193 Pa and its RMS difference is 23.1 Pa. The maximum
liquid-velocity component differs by `0.00142 m/s`.

The r16 closed runs are likewise finite and one-dimensional, but a subsequent
boundary audit invalidated them for physical comparison. At 0.6 s the dry-zone
mean gas pressures are 40.4 and 32.8 Pa for `1e-5` and `5e-6 s`, respectively,
while the four prescribed top particles have `p_g=994.734 Pa`. The previous
`surface_gas_pressure_ratio=1` applied the liquid ponding pressure directly to
the residual gas phase. Siemens et al. state that confined air is released
upward through the transmission zone and that some bubbles reach the surface
and dissipate into the atmosphere. Thus the closed test has no lateral or base
gas venting, but its surface gas outlet is atmospheric.

Keep `fixed_gas_pressure=false` so that interior gas pressure evolves, and set
`surface_gas_pressure_ratio=0` so that air can leave only through the upper
wet zone. Repeat only the two closed 0.6 s diagnostics:

```bash
python prepare_stability_checks.py --duration 0.6 \
  --revision r17-atmospheric-gas-outlet --time-steps 1e-5 5e-6 \
  --pic 0 --pic-t 0 --pressure-smoothing false
sbatch --array=2-3 run_stability_hpc4.sh
```

The checker now requires all four top particles to retain zero gauge gas
pressure. This retroactively rejects the r16 closed outputs and prevents the
experiment's measured pressure plateau from being imposed as a numerical
boundary condition.

Both r17 closed runs pass the range, symmetry, wet-region connectivity, and
atmospheric-outlet checks. At 0.6 s the dry-zone mean gas pressures are 21.03
and 13.48 Pa for `1e-5` and `5e-6 s`, respectively; the maximum and RMS gas
pressure differences are 7.72 and 7.55 Pa. The saturation maximum/RMS
differences are `6.24e-4` and `6.52e-5`, and the maximum liquid-pressure
difference is 76.3 Pa. Because the coarser-step mean gas pressure is about 56%
above the finer result, reject `1e-5 s` for the gas-pressure validation and
test whether `5e-6 s` is converged with one geometric-refinement reference:

```bash
python prepare_stability_checks.py --duration 0.6 \
  --revision r18-gas-dt-refinement --time-steps 2.5e-6 \
  --pic 0 --pic-t 0 --pressure-smoothing false
sbatch --array=1 run_stability_hpc4.sh
```

With one requested time step, manifest row 0 is the open case and row 1 is the
closed case, so only array task 1 is submitted. It advances 240,000 steps and
does not repeat the already-screened open column.

The r18 closed reference passes every range and boundary check. At 0.6 s its
dry-zone mean gas pressure is 9.67 Pa, compared with 13.48 Pa at `5e-6 s`.
The coarse-to-medium and medium-to-fine gas-pressure RMS differences are 7.51
and 3.79 Pa, giving an observed temporal order of 0.99. The corresponding
fine-step error estimate is 3.87 Pa and the safety-factored grid-convergence
index (GCI) is 4.84 Pa, or 0.40% of the experiment's lower 1.2 kPa pressure
plateau. The saturation RMS differences likewise decrease from `6.49e-5` to
`3.23e-5`, with a fine-step GCI of `4.02e-5`. The liquid-pressure RMS
difference between the last two steps is 9.31 Pa; its GCI is 30.6 Pa, or 3.08%
of the 994.734 Pa applied ponding pressure. The largest local liquid-pressure
difference, 68.4 Pa, occurs in the dry layer immediately below the wetted top
cell. These results reject `5e-6 s` but support `2.5e-6 s` as the production
candidate on the experimental pressure and saturation scales.

Before launching the 400 and 900 s production runs, extend both final boundary
configurations to 3 s at the selected candidate step:

```bash
python prepare_stability_checks.py --duration 3 \
  --revision r19-production-step-3s --time-steps 2.5e-6 \
  --pic 0 --pic-t 0 --pressure-smoothing false
sbatch --array=0-1 run_stability_hpc4.sh
```

The open and closed tasks each advance 1,200,000 steps and write 21 particle
outputs. If both pass the range, layer-symmetry, wet-connectivity, and surface
gas-pressure checks, update the full-duration inputs to the selected step.

The r19 closed task completes all 1,200,000 steps and passes every check. Its
maximum absolute pressure and velocity component are 994.734 Pa and
`0.06278 m/s`; layer symmetry is retained to `1.39e-17`, the wettest
unconnected layer remains at saturation 0.0320, and the four top particles
remain at zero gauge gas pressure. The open task remains finite and symmetric
through step 1,140,000 (2.85 s), with maximum absolute pressure and velocity
of 1.412 kPa and `0.06279 m/s`, but reaches the original one-hour Slurm limit
before writing the required 3 s output and pass marker. Repeat only the open
case under a fresh result identifier; the stability script now requests two
hours:

```bash
python prepare_stability_checks.py --duration 3 \
  --revision r19b-open-3s --time-steps 2.5e-6 \
  --pic 0 --pic-t 0 --pressure-smoothing false
sbatch --array=0 run_stability_hpc4.sh
```

Do not rerun manifest row 1: the original r19 closed directory is already the
accepted 3 s result. The new revision also prevents the incomplete 20-file
open directory from being mistaken for a complete rerun.

The r19b open rerun completes all 1,200,000 steps and passes every check. Its
maximum pressure and velocity component are 1.412 kPa and `0.06279 m/s`, the
maximum within-layer saturation difference is `6.94e-18`, no disconnected wet
region develops, and the surface gas pressure remains zero. Its first 20
outputs reproduce the incomplete r19 run to roundoff (zero gas-pressure
difference and at most about `1e-12` in the other pressure arrays). Together
with the accepted r19 closed result, this selects `2.5e-6 s` for production.

At the observed throughput, however, the 400 and 900 s runs cannot fit within
one 24-hour job. Validate HDF5 restart continuity before constructing the
production chain. The r20 test runs the open column from 0 to 1.5 s, writes a
checkpoint, resumes with global steps 600001--1200000, and compares its final
state against the continuous r19b 3 s result. The installation patch first
sizes the HDF5 read buffer and prevents the three-phase solver from resetting
the resumed step counter. Pull and rebuild once, then submit both dependent
segments:

```bash
sbatch program_patch/rebuild_hpc4.sh
# Wait for the rebuild to complete successfully, then:
bash submit_checkpoint_validation_hpc4.sh
```

The submission helper prepares both inputs and submits segment B with an
`afterok` dependency on segment A. Segment B checks all field ranges and writes
`CHECKPOINT_COMPARISON_PASSED.txt` under
`stability_results/siemens2013-r20-checkpoint-open-segment-b/` only if its
coordinates, saturations, phase pressures, permeabilities, and velocities
match the continuous reference within the documented tolerances.

The r20 test passes. Both segments pass all range checks, and the resumed 3 s
state is roundoff-equivalent to the continuous r19b result: coordinates, gas
pressure, and both permeability fields are identical; the largest liquid
pressure/suction difference is `2.84e-12 Pa`; the largest saturation
difference is `3.47e-18`; and the largest phase-velocity difference is below
`8.48e-16 m/s`. HDF5 restart is therefore accepted for production.

After r20 was completed, the fixed-gas formulation used by the open column
was corrected to solve the reduced liquid balance `f_w / K_ww` while holding
the gas-pressure rate at zero. This change does not alter the HDF5 restart
layout or invalidate r20's restart-continuity result, but it does change the
open-column evolution. Rebuild and repeat only the selected 3 s open check
under a fresh UUID before production:

```bash
cd /home/xchenjm/chen/MPM/Backup-mpm
git pull origin agent/siemens-column-validation
cd Example/siemens2013_column_validation
sbatch program_patch/rebuild_hpc4.sh
# Wait for a successful rebuild, then:
python prepare_stability_checks.py --duration 3 \
  --revision r21-fixed-gas-3s --time-steps 2.5e-6 \
  --pic 0 --pic-t 0 --pressure-smoothing false
sbatch --array=0 run_stability_hpc4.sh
```

Do not submit the production chains until the r21 directory contains
`RANGE_CHECK_PASSED.txt`. The closed case uses evolving gas pressure and does
not need to be repeated for this source change.

Each array task requests one GPU and 32 CPUs from `granularmech` under
`comgranmech`. It exits nonzero when the executable is stale, the initial
suction is not 972.99 Pa, saturation leaves [0, 1], or any phase pressure
exceeds the deliberately loose 1 MPa safety bound. It also rejects non-finite
arrays, nonpositive phase permeability, inconsistent phase saturations, and
velocity components above the loose 10 m/s bound. The checker now also rejects
a saturation difference greater than `1e-4` between the two particles in any
horizontal layer, a disconnected wet band, or non-atmospheric gas pressure in
the four top particles. Passing directories contain `RANGE_CHECK_PASSED.txt`
under `stability_results/`.

The checker decodes the compressed VTP arrays themselves rather than trusting
the XML `RangeMin`/`RangeMax` metadata, because VTK range metadata silently
ignores NaN values.

## Restartable production run on HPC4

The committed production inputs use `dt=2.5e-6 s`, `PIC=PIC_T=0`, no
step-count-based pressure smoothing, and HDF5 output. The 400 s open case has
160,000,000 updates and the 900 s closed case has 360,000,000 updates. Based
on the measured 3 s throughput, a 50 s segment should take about 17.5 hours,
leaving useful margin under the 24-hour limit. The generator creates eight
open segments and eighteen closed segments. The two cases run independently,
while each case's segments are sequential `afterok` jobs.

After r21 passes, first inspect the exact production commands without
submitting:

```bash
cd /home/xchenjm/chen/MPM/Backup-mpm/Example/siemens2013_column_validation
python submit_production_hpc4.py --dry-run
```

If the dry run reports 8 open and 18 closed jobs, submit both chains once:

```bash
python submit_production_hpc4.py
```

Every job requests one GPU and 32 CPUs from `granularmech` under
`comgranmech`. It checks both installed source copies, rejects a stale
executable, verifies the preceding HDF5 checkpoint and completion marker,
runs the range checks, and writes `SEGMENT_COMPLETED.txt` only after the final
VTP and HDF5 files exist. Inputs and the submitted job IDs are recorded under
`production_inputs/`; outputs use unique directories such as
`results/siemens2013-production-open-seg000/`.

If a dependency chain stops, inspect the failed segment first. Do not mix an
incomplete directory with a rerun. After moving the incomplete result and its
log out of the case directory, resume the completed prefix and submit the
remaining jobs with:

```bash
python submit_production_hpc4.py --resume
```

Set and export `MPM_SOURCE` or `MPM_BIN` only if the executable locations
move. Do not submit `run_hpc4.sh`; it intentionally exits before starting the
unsegmented calculation.

### Production segment 0 audit and permeability/time acceleration test

The first 50 s open production segment is numerically complete and passes all
range checks, but it fails the experimental wetting-front comparison. It has
41 particle outputs through global step 20,000,000, zero gas pressure, a
maximum absolute pressure of 994.734 Pa, and no lateral or disconnected-band
instability. Nevertheless, only the prescribed top cell exceeds `Sw=0.40` at
50 s. The first interior cell reaches `Sw=0.0767`, and material about 100 mm
below the boundary remains at the initial `Sw=0.0301`. The simulated wetting
depth is therefore 0 mm at the standard threshold (about 11 mm even at a
threshold of 0.05), whereas the digitised experiment gives approximately
212 mm at 50 s. A numerical range pass is not an experimental validation.

Do not resume the production chains until this mismatch is resolved. A full
three-phase semi-implicit conversion would require two coupled pressure
unknowns and is not the existing single-pore-pressure two-phase solver. First
test the cheaper alternative: multiply the common intrinsic permeability by
10 and 100 while reducing computed time by the same factors. This is accepted
only if it reproduces the unscaled state after mapping back to the 50 s
reference time; it is not a calibration of the experimental permeability.

Generate the two open tests (computed durations 5 and 0.5 s) and submit only
manifest rows 0 and 1:

```bash
python prepare_stability_checks.py --duration 50 \
  --revision r22-k-time-scaling --time-steps 2.5e-6 \
  --permeability-scales 10 100 \
  --pic 0 --pic-t 0 --pressure-smoothing false
sbatch --array=0-1 --time=04:00:00 run_stability_hpc4.sh
```

The generator records the original 50 s duration, actual computed duration,
and permeability multiplier in both JSON and the manifest. The comparison
script divides the scaled permeability and velocity arrays by the multiplier
before comparing every hydraulic field. The `k x 100` route would reduce the
event durations to 4 s open and 9 s closed only if this overlap test passes.

Both r22 runs are stable but fail that equivalence requirement. At the mapped
50 s state, the first interior layer has `Sw=0.07670` in the unscaled run,
`0.03835` for `k x 10`, and `0.03584` for `k x 100`. The maximum liquid-
pressure differences are 244 and 415 Pa. The discrepancy begins at the first
2.5 s comparison output and grows monotonically, even though the input files
differ only in permeability, duration, and output controls. Permeability/time
scaling is therefore rejected and the corresponding closed tests must not be
submitted.

### Semi-implicit three-phase pressure smoke test

The validation patch now provides an opt-in global backward-Euler pressure
solve with liquid and gas pressure as separate unknowns. Solid mechanics and
the final three-phase momentum update remain explicit. The pressure system
uses the existing SWRC tangent, liquid compressibility, ideal-gas storage,
phase permeabilities, gravity, and the same `viscosity / permeability` drag
convention as the explicit momentum equation. The open case reduces to
one liquid-pressure block when `fixed_gas_pressure=true`; the closed case
retains the full coupled two-pressure block. The historical solver remains
unchanged unless the input contains:

```json
"pressure_integration": "semi_implicit"
```

The first r23 semi-implicit launch stopped before step 0 because a few GIMP
support nodes were mechanically active but had zero pressure projection
weight. The r24 correction uses a compact pressure-specific DOF map containing
only nodes with nonzero shape-function or gradient support. Gradient-only
nodes are initialised with the volume-weighted particle pressure mean; this
removes the invalid zero-weight division without changing the pressure
boundary conditions or material parameters.

All r24 open jobs then produced eleven outputs, but they did not pass the
comparison: the `2.5e-5 s` run had larger saturation error than the `5e-4 s`
run and advanced the wetting front much farther. This inverse time-step trend
identified repeated absolute particle-node-particle pressure projection as
artificial per-step diffusion. The r25 correction solves the same global
backward-Euler system but transfers only the nodal pressure increment to the
particles. With zero physical pressure change, the transferred increment is
exactly zero and the particle pressure is not smoothed.

The r25 open results reject `5e-4 s`: pressure grew to `8.42 MPa`, velocity to
about `86 m/s`, and the disconnected-wet-layer check failed. The `2.5e-5` and
`1e-4 s` runs both remained bounded and were nearly identical. At `0.05 s`,
only four particles immediately below the prescribed top cell had saturation
errors above `0.001`; all remaining particles were below `1e-4`. Two of those
four particles were fully wetted in the semi-implicit solution but remained
dry in the explicit solution. This local error points to the `1e6` weak
boundary penalty, not continuing bulk projection diffusion.

Before extending the duration or running the closed column, keep the r25
executable and perform a `0.005 s`, `dt=1e-4 s` penalty sweep. No rebuild is
needed after pulling this input-script update:

```bash
python prepare_semi_implicit_checks.py \
  --duration 0.005 \
  --semi-implicit-dts 1e-4 \
  --boundary-penalties 1 1e2 1e4 1e6 \
  --revision r26-boundary-penalty-sweep
sbatch --array=0-4 run_semi_implicit_hpc4.sh
```

Compare the five open jobs with:

```bash
python compare_semi_implicit.py --case open
```

The r26 sweep rejects penalty `1`, which reached `301 MPa` and failed the
range checks. Penalties `100`, `1e4`, and `1e6` were all bounded and nearly
identical; `100` gave the smallest pressure and velocity errors and is selected
as the least stiff stable boundary penalty. The existing passed r21 fixed-gas
open run already supplies a 3 s explicit baseline at `dt=2.5e-6 s`, so the
next stage should not repeat its 1.2 million explicit steps.

Generate the r27 manifest and submit only row 1, the 3 s semi-implicit open
case (`30,000` steps):

```bash
python prepare_semi_implicit_checks.py \
  --duration 3 \
  --semi-implicit-dts 1e-4 \
  --boundary-penalty 100 \
  --revision r27-open-3s
sbatch --array=1 run_semi_implicit_hpc4.sh
```

If the old r21 result remains on the server, compare r27 at the eleven shared
physical times without rerunning the explicit reference:

```bash
python compare_semi_implicit.py --case open \
  --open-reference-dir \
  stability_results/siemens2013-stability-r21-fixed-gas-3s-open_2p5e-06 \
  --open-reference-dt 2.5e-6
```

The r27 output did not pass the range checker. At `0.3 s`, it already contained
a fully saturated layer near `y=0.989 m` disconnected from the prescribed top
water layer. The left/right saturation spread reached `0.949` by `0.6 s`.
These alternating wet and dry layers are a non-monotone consistent-mass
oscillation at the sharp infiltration front, rather than physical one-
dimensional infiltration. The next solver revision row-sum lumps the liquid/
gas storage and weak pressure-boundary matrices but leaves the spatial Darcy
diffusion matrix unchanged.

After rebuilding, first run only a `0.3 s` r29 open smoke test. This is long
enough to reach the first r27 disconnected band but costs only 3,000
semi-implicit steps:

```bash
python prepare_semi_implicit_checks.py \
  --duration 0.3 \
  --semi-implicit-dts 1e-4 \
  --boundary-penalty 100 \
  --revision r29-lumped-pressure-0p3s
sbatch --array=1 run_semi_implicit_hpc4.sh
```

The r29 run produced all eleven outputs and passed every range check. Pressure
remained below `1.343 kPa`, the largest phase-velocity component was
`0.0802 m/s`, and the saturation difference between the two particles in each
horizontal layer was below `1.4e-13`. The apparent maximum non-connected
saturation of `0.394` occurred in the single transition layer immediately
below the top-connected wet region at `0.27 s`; at `0.30 s` it crossed the
`0.40` front threshold and joined that region. The profile is therefore a
continuous advancing front rather than a detached band. The front depth at
`0.30 s` was `17.1 mm`.

During this audit, the semi-implicit Darcy mobility was found to contain an
obsolete factor of `0.981` that is absent from the current three-phase drag
coefficient. The r30 correction uses exactly `permeability / viscosity` for
both phases. Rebuild it, then repeat the short open test at two time steps to
check that the lumped formulation is time-step convergent before extending
the duration:

```bash
sbatch program_patch/rebuild_hpc4.sh
# After the rebuild succeeds:
python prepare_semi_implicit_checks.py \
  --duration 0.3 \
  --semi-implicit-dts 2.5e-5 1e-4 \
  --boundary-penalty 100 \
  --revision r30-correct-darcy-mobility-0p3s
sbatch --array=1-2 run_semi_implicit_hpc4.sh
```

After both rows pass, compare `dt=1e-4 s` directly with the finer
semi-implicit result:

```bash
python compare_semi_implicit.py --case open \
  --semi-implicit-reference-dt 2.5e-5
```

Both r30 jobs produced eleven outputs and passed every range check. The
`dt=2.5e-5 s` and `dt=1e-4 s` solutions gave the same top-connected wetting-
front depth at every output time, including `17.1 mm` at `0.30 s`. The maximum
saturation difference (`0.0479`) was confined to the single 5.7 mm transition
layer. The nominal maximum pressure difference (`505 Pa` at `0.18 s`) occurred
in that same layer when the coarse-step solution had just reached full
saturation and its capillary suction became zero; it was not a bulk dry-zone
pressure discrepancy. The RMS saturation and liquid-pressure differences
over all particles and times were only `0.00268` and `15.0 Pa`, respectively.
The `1e-4 s` step is therefore accepted for the longer smoke test.

No rebuild is needed after r30. Generate and submit only the 3 s open
semi-implicit row (`30,000` steps):

```bash
python prepare_semi_implicit_checks.py \
  --duration 3 \
  --semi-implicit-dts 1e-4 \
  --boundary-penalty 100 \
  --revision r31-correct-darcy-mobility-3s
sbatch --array=1 run_semi_implicit_hpc4.sh
```

The rebuild installs the tracked particle headers and implementations plus
the pressure-solver header and implementation. After a successful build it
also removes obsolete `.pre-siemens-validation` and `.before-*` copies of the
two duplicated particle sources; the current sources remain recoverable from
the tracked replacements. Each smoke-test task performs the usual decoded-VTP
range checks. The r23--r27 sequence above supersedes the original four-job
smoke-test command. Do not start the closed or production calculation until
the r31 three-second open run passes.

## Post-process after both runs

The script needs Python packages `vtk`, `numpy`, and `matplotlib`:

```bash
python postprocess_validation.py
```

It requires all expected global output steps (321 open and 361 closed), merges
the segment directories by numeric global step, and rejects missing,
duplicate, incomplete, or out-of-order segments. It then writes two summary
CSV files and `siemens2013_validation.png` under `validation_results/`. The
wetting front is the deepest layer in the wet
region connected continuously to the top. This excludes any lower wet region
from the front measurement and, together with the range checker, prevents an
isolated numerical band from being treated as infiltration. The default
mean-saturation threshold is 0.40 and can be changed with
`--saturation-threshold`. The script also rejects outputs when the SWRC
initialisation is absent or pressures have blown up, so an invalid run cannot
silently produce a manuscript figure.
