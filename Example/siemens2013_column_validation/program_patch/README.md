# Program patch for `mpm_hpc`

`initial_suction_from_swrc.patch` changes only
`include/particles/particle_threephase_new.tcc`. It removes the mandatory
20 psi initial suction for this particle type while preserving that value as
the default for all existing inputs.

The patch adds two optional material properties:

- `initial_suction`: an explicitly prescribed initial suction in Pa;
- `initial_suction_from_swrc`: when `true`, initialise suction from the input
  liquid saturation and van Genuchten parameters (`para_p0`, `para_m`).

The Siemens inputs use the second option because the published initial state
is close to residual saturation and is incompatible with a hard-coded 20 psi
value.

## Apply on hpc4

From the root of the `mpm_hpc` source checkout, run:

```bash
git apply --check --ignore-space-change --ignore-whitespace \
  /absolute/path/to/Backup-mpm/Example/siemens2013_column_validation/program_patch/initial_suction_from_swrc.patch
git apply --ignore-space-change --ignore-whitespace \
  /absolute/path/to/Backup-mpm/Example/siemens2013_column_validation/program_patch/initial_suction_from_swrc.patch
```

Then rebuild `mpm_hpc`. If the patch has already been applied,
`git apply --check` will stop without changing the source.
