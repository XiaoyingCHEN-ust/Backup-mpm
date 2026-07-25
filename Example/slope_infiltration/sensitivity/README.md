# Reviewer sensitivity analysis

## Purpose

The 0.5 m value in the current manuscript is a configuration-specific point at
which the two onset-time curves begin to separate. It is not a general failure
criterion. The smallest useful revision is:

1. a one-at-a-time permeability sensitivity with
   \(k/k_0=0.1,1,10\), where \(k_0=1.5\times10^{-12}\ {\rm m^2}\);
2. an optional gas-boundary sensitivity with
   \(P_{g,b}/P_{l,b}=0,0.5,1\) at \(k=k_0\);
3. an explicit limitation covering initial saturation, SWRC parameters,
   geometry/slope angle, strength and other gas boundary conditions.

The same six surface heads used in the manuscript are retained:
0, 0.5, 1.0, 1.5, 2.0 and 2.5 m. Every variable-gas case is paired with a
fixed-gas reference at the same permeability and water head.

## New JSON controls

The calculation JSON now controls the former hard-coded C++ choices:

- `fixed_gas_pressure`: `false` evolves gas pressure; `true` freezes the
  initial equilibrated gas-pressure field.
- `surface_liquid_pressure`: imposed surface water pressure in Pa.
- `surface_gas_pressure_ratio`: \(P_{g,b}/P_{l,b}\) for the variable-gas
  model. `1.0` reproduces the co-pressurized manuscript boundary and `0.0`
  represents an atmospheric/vented gas boundary.

The defaults reproduce the current variable-gas 25 kPa case.

## Server workflow

From `Example/slope_infiltration`:

```bash
python3 sensitivity/generate_cases.py
```

Rebuild the MPM executable once after replacing the two
`particle_threephase_new.tcc` files:

```bash
cd /path/to/mpm
cmake --build build -j 12
```

Stage 1 runs 36 permeability cases (three permeabilities, six heads and two gas
models):

```bash
cd /path/to/Backup-mpm/Example/slope_infiltration
python3 sensitivity/run_cases.py \
  --mpm-bin /path/to/mpm/build/mpm \
  --family permeability \
  --jobs 4
```

Set `--jobs 1` if memory is limited. The optional Stage 2 adds 12 new
gas-boundary cases; the 12 shared reference cases from Stage 1 are skipped
automatically because their result directories already exist:

```bash
python3 sensitivity/run_cases.py \
  --mpm-bin /path/to/mpm/build/mpm \
  --family gas_boundary \
  --jobs 4
```

Do not use `--rerun` unless an existing result directory is intentionally to
be reused by the solver. `--skip-legacy-reference` is provided only for quick
regression checks against the old manuscript results. It is not recommended
for the final sensitivity figure because the current restart and source files
should be applied consistently to every curve.

## Reduce and return results

The extraction step needs Python packages `vtk` and `numpy`:

```bash
python3 -m pip install vtk numpy
python3 sensitivity/extract_metrics.py --family permeability
python3 sensitivity/estimate_onset.py
```

After Stage 2, extract all cases instead:

```bash
python3 sensitivity/extract_metrics.py
python3 sensitivity/estimate_onset.py
python3 sensitivity/plot_sensitivity.py
```

Only these small files need to be returned for analysis:

- `sensitivity/manifest.csv`
- `sensitivity/metrics.csv`
- `sensitivity/onset_times.csv`
- `sensitivity/sensitivity_summary.csv` if plotting was run

The default reproducible onset measure is the linearly interpolated time when
the 99th-percentile particle displacement has increased by 0.20 m from the
equilibrated restart value. The plotted light bands repeat the calculation at
0.15 and 0.25 m so the conclusion is not tied to one displacement cutoff.

## Expected manuscript figure

`plot_sensitivity.py` creates one two-panel figure:

- panel (a): \(\Delta t=t_{\rm fixed}-t_{\rm variable}\) versus surface head
  for the three permeabilities;
- panel (b): the same quantity for the three surface gas-pressure ratios.

If only Stage 1 is required, panel (a) alone is sufficient together with a
clear limitation paragraph.
