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

## Files

- `preprocess.py` generates the mesh, particles, volumes, temperatures, and
  entity sets.
- `mpm_open.json` and `mpm_closed.json` are the two solver inputs.
- `reference_data/` contains values read from Figures 7, 8, and 11 and the
  pressure range stated in the paper. The wetting-front points are approximate
  digitisation targets and should not be presented as exact tabulated data.
- `run_hpc4.sh` is a two-task Slurm array script.
- `postprocess_validation.py` extracts wetting-front depth and dry-zone
  pore-air pressure from the particle VTP files.
- `program_patch/` contains the small, backward-compatible source patch
  required to initialise suction consistently with the published saturation.
- `validation_subsection_draft.md` is intentionally left with placeholders
  until the HPC results are available.

## Generate and check the case

From this directory:

```bash
python preprocess.py
```

Expected counts are 376 cells, 1504 particles, eight top-boundary particles,
and eight bottom-boundary particles.

## Run on hpc4

Apply the patch in `program_patch/` and compile the program. The current HPC4
script loads `miniconda3/24.3.0-quc3pyu`, activates `cbgeo_tbb`, requests one
GPU and 32 CPUs, and defaults to the executable at
`/home/xchenjm/chen/MPM/mpm/build/mpm`. Submit both cases with:

```bash
sbatch run_hpc4.sh
```

Set and export `MPM_BIN` only if the executable is stored elsewhere.

The array index 0 runs the open test and index 1 runs the closed test. Results
are written to `results/siemens2013-open/` and
`results/siemens2013-closed/`.

## Post-process after both runs

The script needs Python packages `vtk`, `numpy`, and `matplotlib`:

```bash
python postprocess_validation.py
```

It writes two summary CSV files and `siemens2013_validation.png` under
`validation_results/`. The default wetting-front definition is the deepest
horizontal layer with mean liquid saturation of at least 0.40; the threshold
can be changed with `--saturation-threshold` for a bounded interpretation
check.
