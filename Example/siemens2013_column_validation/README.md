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
boundary. It is a 45.6 mm-wide plane-strain rectangle. Four 11.4 mm cells are
used across the width and 94 along the height, giving a modeled height of
1.0716 m. The skeleton is fixed to isolate hydraulic liquid–gas coupling.

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
- `run_stability_hpc4.sh` is an eight-task, short-duration time-step check.
- `run_hpc4.sh` is the two-task full-duration Slurm array script.
- `postprocess_validation.py` extracts wetting-front depth and dry-zone
  pore-air pressure from the particle VTP files.
- `program_patch/` installs the backward-compatible source replacement,
  rebuilds the executable, and prevents a stale executable being used again.
- `validation_subsection_draft.md` is intentionally left with placeholders
  until the HPC results are available.

## Generate and check the case

From this directory:

```bash
python preprocess.py
```

Expected counts are 376 cells, 1504 particles, eight top-boundary particles,
and eight bottom-boundary particles.

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

Next generate and submit eight 3 s checks (open/closed, each at
`5e-4`, `5e-5`, `5e-6`, and `1e-6` s):

```bash
python prepare_stability_checks.py
sbatch run_stability_hpc4.sh
```

The 3 s window covers the time by which the first run had already diverged;
the job scripts retain only every 10000th routine step message to keep the
Slurm logs small.

Each array task requests one GPU and 32 CPUs from `granularmech` under
`comgranmech`. It exits nonzero when the executable is stale, the initial
suction is not 972.99 Pa, saturation leaves [0, 1], or any phase pressure
exceeds the deliberately loose 1 MPa safety bound. It also rejects non-finite
arrays, nonpositive phase permeability, inconsistent phase saturations, and
velocity components above the loose 10 m/s bound. Passing directories contain
`RANGE_CHECK_PASSED.txt` under `stability_results/`.

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
