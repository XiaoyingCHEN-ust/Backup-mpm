# MPM source replacements

These two files are identical copies kept in both source locations used by the
current checkout. On the Linux server, replace the existing files and rebuild:

```bash
MPM_SOURCE=/path/to/mpm

cp "$MPM_SOURCE/include/particles/particle_threephase_new.tcc" \
   "$MPM_SOURCE/include/particles/particle_threephase_new.tcc.pre-sensitivity"
cp "$MPM_SOURCE/include/solvers/particle_threephase_new.tcc" \
   "$MPM_SOURCE/include/solvers/particle_threephase_new.tcc.pre-sensitivity"

cp server_overrides/mpm/include/particles/particle_threephase_new.tcc \
   "$MPM_SOURCE/include/particles/particle_threephase_new.tcc"
cp server_overrides/mpm/include/solvers/particle_threephase_new.tcc \
   "$MPM_SOURCE/include/solvers/particle_threephase_new.tcc"

cmake --build "$MPM_SOURCE/build" -j 12
```

The `.pre-sensitivity` files are local backups and should not be committed.

The replacement adds three optional liquid-material properties:

- `fixed_gas_pressure`
- `surface_liquid_pressure`
- `surface_gas_pressure_ratio`

When they are omitted, the defaults preserve the current variable-gas 25 kPa
surface boundary.
