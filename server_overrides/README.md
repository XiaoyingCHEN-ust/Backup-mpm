# MPM source replacements

These two files are identical copies kept in both source locations used by the
current checkout. For HPC4, the Siemens validation installer performs the
replacement, backup, and rebuild automatically:

```bash
cd /home/xchenjm/chen/MPM/Backup-mpm/Example/siemens2013_column_validation
sbatch program_patch/rebuild_hpc4.sh
```

The `.pre-siemens-validation` files are local backups and should not be
committed.

The replacement adds three optional liquid-material properties:

- `fixed_gas_pressure`
- `surface_liquid_pressure`
- `surface_gas_pressure_ratio`

When they are omitted, the defaults preserve the current variable-gas 25 kPa
surface boundary.

For the laboratory validation it additionally adds `initial_suction` and
`initial_suction_from_swrc`, initialises the retention-curve derivative, and
updates effective saturation consistently so that phase permeabilities can
evolve during infiltration. The optional `relative_permeability_floor`
regularises residual phase endpoints so the drag calculation never divides
by zero permeability.
