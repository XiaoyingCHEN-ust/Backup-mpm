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
  bottom boundary has constant pressure.
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
- `run_hpc4.sh` is the two-task full-duration Slurm array script.
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
complete top cell), and two bottom-boundary particles.

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

Each array task requests one GPU and 32 CPUs from `granularmech` under
`comgranmech`. It exits nonzero when the executable is stale, the initial
suction is not 972.99 Pa, saturation leaves [0, 1], or any phase pressure
exceeds the deliberately loose 1 MPa safety bound. It also rejects non-finite
arrays, nonpositive phase permeability, inconsistent phase saturations, and
velocity components above the loose 10 m/s bound. The checker now also rejects
a saturation difference greater than `1e-4` between the two particles in any
horizontal layer. Passing directories contain `RANGE_CHECK_PASSED.txt` under
`stability_results/`.

The checker decodes the compressed VTP arrays themselves rather than trusting
the XML `RangeMin`/`RangeMax` metadata, because VTK range metadata silently
ignores NaN values.

Download `stability_results/` after this short array. The stable step will be
selected before the full inputs are changed; do not resubmit `run_hpc4.sh`
with the current `dt=5e-4 s` inputs.

## Full run on HPC4 (after the short check)

The full script loads `miniconda3/24.3.0-quc3pyu`, activates `cbgeo_tbb`, and
defaults to `/home/xchenjm/chen/MPM/mpm/build/mpm`. It also verifies that both
source copies contain the patch and are not newer than the executable. Once
the selected time step has been committed, submit with:

```bash
sbatch run_hpc4.sh
```

Set and export `MPM_SOURCE` or `MPM_BIN` only if those locations move. Array
index 0 runs the open test and index 1 the closed test. Full results are
written to `results/siemens2013-open/` and
`results/siemens2013-closed/`.

## Post-process after both runs

The script needs Python packages `vtk`, `numpy`, and `matplotlib`:

```bash
python postprocess_validation.py
```

It writes two summary CSV files and `siemens2013_validation.png` under
`validation_results/`. The wetting front is the deepest layer in the wet
region connected continuously to the top. This excludes the separate
saturated constant-head layer at the base of the open column. The default
mean-saturation threshold is 0.40 and can be changed with
`--saturation-threshold`. The script also rejects outputs when the SWRC
initialisation is absent or pressures have blown up, so an invalid run cannot
silently produce a manuscript figure.
