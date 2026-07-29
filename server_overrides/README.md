# MPM source replacements

The particle header and implementation are identical copies kept in both
source locations used by the current checkout. The three-phase solver header
and implementation additionally contain the opt-in coupled semi-implicit
pressure path. For HPC4, the Siemens validation installer performs the
replacement and rebuild automatically:

```bash
cd /home/xchenjm/chen/MPM/Backup-mpm/Example/siemens2013_column_validation
sbatch program_patch/rebuild_hpc4.sh
```

The tracked replacements are the recovery copies. After a successful rebuild,
the installer removes obsolete `.pre-siemens-validation` and `.before-*`
particle-source backups so a stale file cannot be mistaken for the compiled
source.

The replacement adds three optional liquid-material properties:

- `fixed_gas_pressure`
- `surface_liquid_pressure`
- `surface_gas_pressure_ratio`

When they are omitted, the defaults preserve the current variable-gas 25 kPa
surface boundary.

For `fixed_gas_pressure: true`, the gas-pressure rate is set to zero and the
liquid-pressure rate is evaluated from the reduced balance `f_w / K_ww`, which
matches the fixed-pressure formulation used for the manuscript comparison.

For the laboratory validation it additionally adds `initial_suction` and
`initial_suction_from_swrc`, initialises the retention-curve derivative, and
updates effective saturation consistently so that phase permeabilities can
evolve during infiltration. The optional `relative_permeability_floor`
regularises residual phase endpoints so the drag calculation never divides
by zero permeability.

Set `analysis.pressure_integration` to `"semi_implicit"` to assemble the
liquid/gas backward-Euler pressure system. The default remains `"explicit"`.
Its Darcy operator uses the same `permeability / viscosity` mobility implied
by the current three-phase drag coefficient. The new integration path remains
experimental. The r30 open time-step check, r31 three-second open test, and
coupled-pressure r32 short closed time-step check passed. Longer closed runs
revealed a particle-scale alternating mode near the wetting front. Boundary-
penalty, diffusion-matrix, gravity-removal, and time-step diagnostics did not
eliminate that mode. The r39 local pressure bound removed the oscillation but
approximately doubled the early wetting-front depth, so it is retained only
as an opt-in diagnostic (`semi_implicit_pressure.bounded_transfer`, default
`false`). The r40 check keeps pure incremental pressure transfer and instead
uses one particle per background cell while preserving two particles per
horizontal layer and four surface-loaded particles.
