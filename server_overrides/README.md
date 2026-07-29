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

The installer additionally patches the shared `mpm_base.tcc` VTK whitelist
in place. This permits the opt-in `force_liquid_pressures` and
`force_gas_pressures` diagnostics in both the input-validation whitelist and
the scalar writer while leaving the default output fields for other particle
types unchanged. The guarded patcher is tracked with the Siemens case and
refuses an unrecognised or partially patched source layout.

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
`false`). The r40 and r41 mesh checks retained an unstable incrementally
accumulated gas-pressure gradient irrespective of whether a horizontal layer
occupied one or two cells. The independent opt-in setting
`semi_implicit_pressure.reconstruct_gradient` reconstructs only pressure
gradients from the current nodal solution; particle pressure values remain
pure incremental FLIP and are not bounded or projected.

The r42 gradient-only result was unchanged because the current 2-D phase
momentum internal force is assembled from particle pressure values rather
than the stored pressure-gradient vectors. The independent r43 option
`semi_implicit_pressure.reconstruct_force` supplies nodally reconstructed
pressure only to that internal force. It does not feed the reconstructed
value back into pressure storage, retention, saturation, or gas density.

For the r46 diagnostic, enabling both reconstruction switches maps the phase
pressure force from particle-quadrature `p grad(N)` to the directly
reconstructed strong form `-N grad(p)`. The separate one-switch behaviours
remain unchanged. This is an opt-in numerical diagnostic rather than a new
default.

The r46 short test preserved physical pressures and saturations to round-off
and eliminated the dry-zone gas-velocity sign alternation. It remains an
experimental opt-in path pending the phase-velocity trend in the subsequent
`0.1 s` test.

The r47 extension confirmed that pure particle FLIP can retain a phase-
velocity null mode even with direct pressure-gradient forces. The independent
`semi_implicit_pressure.reconstruct_darcy_velocity` option reconstructs fluid
velocities from the Darcy relation already assumed by the semi-implicit
pressure operator. It leaves solid `PIC=0` and defaults to `false`.

The r48 `0.05 s` closed test passed with zero adjacent dry-layer velocity sign
changes in both phases. After its first output, both phase maxima decayed, and
physical pressure and saturation matched the pre-velocity-reconstruction
reference to round-off. Continue with the same executable for only `0.1 s`;
no source reinstall or rebuild is needed.

The r49 `0.1 s` extension passed with zero gas and liquid sign changes. Both
phase maxima continued to decay through the interval where the earlier FLIP
mode had reappeared, and the saturation profile remained connected and
monotone. The next accepted step is an unchanged `0.3 s` closed-column run;
it also requires no rebuild.

The r50 `0.3 s` run passed the range, profile, and sign-change checks. Dry-zone
gas velocity continued to decay, while upward gas velocity localized in the
connected wetting-front layer rose to `0.0436 m/s`. The checker now separates
post-initial dry- and wet-zone gas maxima. Continue with the unchanged binary
for only `0.5 s`; no rebuild is needed.

The r51 `0.5 s` result confirmed a finite wet-front gas pulse, peaking at
`0.0832 m/s` before collapsing when that layer approached full saturation.
Dry-zone velocities remained bounded and nonalternating. Continue to one
second with `0.05 s` VTK output spacing; no rebuild is needed.

The r52 one-second run also passed. The next layer joined the connected front
at `0.6 s`, dry-zone gas pressure remained smooth, and both phase sign-change
counts stayed zero. Because the second wet-layer gas pulse was still rising
at `1.0 s`, extend only to `1.2 s` with `0.05 s` outputs. No rebuild is needed.

The r53 `1.2 s` result is rejected. Its second wet-front gas pulse collapsed,
but a dry-column pressure mode appeared after `1.1 s`, reaching 35 adjacent
gas-velocity sign changes at `1.2 s`. Run the unchanged binary at `dt=8e-5 s`
and `5e-5 s` over the same interval to distinguish a time-step error from a
pressure-space mode. No rebuild is needed.

The r54 time-step check confirmed a spatial pressure mode: every time step
first failed at `1.15 s`, and refinement increased the sign-change count. The
new opt-in `semi_implicit_pressure.projection_rate` supplies a time-step-
invariant pressure-only relaxation toward the resolved nodal solution. The
default is zero, so existing inputs and mechanical `PIC=0` are unchanged.
Rebuild before the r55 projection-rate sweep.

All r55 rates removed the pressure mode, but `0.5--10 1/s` altered the front-
layer saturation history and the stronger rates delayed front advance. None
is accepted. Reuse the binary for a weaker `0.01--0.2 1/s` sweep; no rebuild
is required.

The r56 sweep identified `0.1 1/s` as the weakest passing rate; `0.01` and
`0.05 1/s` retained the late mode. The `0.1 1/s` connected-front history was
unchanged and its dry-zone mean gas pressure remained close to r52. Verify it
at `dt=5e-5 s` before extending the duration.

The r57 half-step check passed with the same front history, zero sign changes,
and small pressure/saturation differences. Accept `0.1 1/s` provisionally and
return to `dt=1e-4 s` for a three-second repeated-front diagnostic. No rebuild
is required.
