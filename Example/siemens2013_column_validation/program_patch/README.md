# Program patch used by the Siemens validation

The first HPC4 outputs proved that the executable still used the historical
hard-coded 20 psi suction: the initial VTP value was 137340 Pa rather than the
972.99 Pa required by the input saturation and SWRC. The complete replacement
is therefore stored in both include locations under `server_overrides/`, and
`install_hpc4.sh` installs both copies and rebuilds the exact executable used
by the job script. If either target template has been removed from the MPM
source tree, the installer recreates its parent path and restores the tracked
replacement before rebuilding.

The installer also applies two guarded checkpoint-resume corrections needed
for full-duration jobs longer than the 24-hour queue limit. It sizes the HDF5
read buffer before loading particle records, and makes the three-phase solver
start its loop at the resumed global step instead of silently restarting the
step counter at zero. `patch_checkpoint_resume.py` is idempotent and refuses
to edit an unrecognised source layout.

It also applies a guarded VTK whitelist correction in `mpm_base.tcc`. The
r44 job originally stopped at its first output because the new opt-in
`force_liquid_pressures` and `force_gas_pressures` particle fields were
registered by the three-phase particle but were absent from the global liquid
VTK whitelist. The first whitelist fix allowed r44 to run, but inspection of
the completed VTP files exposed a second hard-coded scalar-writer list that
silently omitted both fields. `patch_vtk_force_fields.py` now updates both
lists without adding the diagnostics to the default output used by other
particle types. The run script verifies both changes and refuses to use an
older executable.

The replacement adds two optional material properties:

- `initial_suction`: an explicitly prescribed initial suction in Pa;
- `initial_suction_from_swrc`: when `true`, initialise suction from the input
  liquid saturation and van Genuchten parameters (`para_p0`, `para_m`).

It also initialises the SWRC tangent before the first pressure update and
updates effective saturation after every saturation update. The latter is
required because the liquid and gas relative permeabilities otherwise remain
frozen at their initial values.

The optional `relative_permeability_floor` prevents division by zero in the
phase-drag coefficient when effective saturation is exactly zero or one. Its
backward-compatible default is `1e-12`; the validation uses `1e-6` to keep the
residual phase effectively immobile without producing an ill-conditioned
drag matrix.

The replacement also corrects the two-dimensional liquid internal force. The
solver already maps liquid gravity through the external body force, and
`PIC_liquid_pressure_` is the actual gauge pressure. The previous expression
subtracted a hard-coded hydrostatic pressure inside the internal force and
therefore introduced an additional gravity gradient. Both two-dimensional
source copies now use `PIC_liquid_pressure_` directly, consistent with the
three-dimensional implementation. Installation and run scripts explicitly
reject the obsolete expression.

For the open column's fixed-gas formulation, the replacement now solves the
reduced liquid balance `f_w / K_ww` and sets the gas-pressure rate to zero.
This avoids solving the coupled 2 x 2 pressure system and then discarding its
gas-pressure increment. Installation and run scripts require this corrected
branch explicitly.

The Siemens inputs use `initial_suction_from_swrc` because the published
initial state is close to residual saturation and is incompatible with a
hard-coded 20 psi value. Existing inputs retain 20 psi by default.

The patch also contains an opt-in semi-implicit pressure integration path for
the three-phase solver. It assembles a backward-Euler finite-element system
for liquid and gas gauge pressures while retaining explicit solid mechanics.
The open validation case automatically reduces to the liquid block when its
gas pressure is fixed; the closed case solves the coupled two-pressure block.
Activate it only with `analysis.pressure_integration = "semi_implicit"`.
Omitting the key, or setting it to `"explicit"`, preserves the historical
update. The r23 run exposed mechanically active GIMP support nodes with zero
pressure projection weight. The r24 path therefore built a compact pressure
DOF map from nodes with nonzero shape-function or gradient support. Although
r24 completed, its smaller time steps moved the wetting front farther because
absolute particle pressure was reprojected through the grid on every step.
The r25 path transfers only the solved nodal pressure increment back to each
particle, so a zero physical increment produces no projection smoothing. The
r27 three-second run subsequently exposed non-monotone saturated bands
disconnected from the top boundary. The pressure storage and weak boundary
terms now use row-sum mass lumping, while the Darcy diffusion blocks remain
consistent. This is intended to preserve a sharp dry/wet front without the
overshoot produced by the consistent pressure mass matrix. The r29 short open
test passed; r30 removes an obsolete `0.981` factor so the implicit Darcy
mobility exactly matches the current three-phase `viscosity / permeability`
drag coefficient. The r30 two-time-step open check passed with an identical
wetting-front position at every output time, and the r31 three-second open
test remained one-dimensional and free of detached wet bands. The r32 closed
time-step check also passed; its dry-zone gas pressure and wetting-front
history converged at `dt=1e-4 s`. The r33 three-second closed test stayed
bounded, but a later adjacent-layer audit detected the onset of the same
nonmonotone mode that invalidated r34. All four r35 boundary penalties failed
at the same layer. The r36 diffusion correction was inactive, zero gravity
only reduced the oscillation, and r38 still failed after halving the step. The
The r39 diagnostic bounded pressure-FLIP extrema outside each particle's
local nodal liquid, gas, and capillary pressure ranges. Although it removed
the alternating mode, it approximately doubled the early wetting-front depth
and is therefore not accepted. The bound is retained only behind
`semi_implicit_pressure.bounded_transfer`, whose default is `false`. The r40
and r41 meshes showed that the short-time gas-velocity instability is
independent of whether each horizontal layer occupies one or two cells. The
r42 option `semi_implicit_pressure.reconstruct_gradient` rebuilt only particle
pressure gradients from the current nodal solution, but code-path inspection
showed that those vectors do not enter the two-dimensional phase momentum
internal force and r42 consequently matched r41. The r43 option
`semi_implicit_pressure.reconstruct_force` uses nodally reconstructed pressure
only for that internal force. Storage, suction, saturation, density, and
reported particle pressures remain pure incremental FLIP. Mechanical `PIC`
and `PIC_T` remain zero.

The r45 diagnostics showed that reconstructed force pressures were smooth,
but mapping them through the particle-quadrature weak form `p grad(N)` still
generated alternating gas velocities. When both `reconstruct_gradient` and
`reconstruct_force` are true, the r46 diagnostic instead maps the directly
reconstructed strong-form pressure force `-N grad(p)`. Either option alone
retains its earlier behaviour. This path is opt-in and must pass the short
hydraulic and velocity comparisons before it can be considered for a longer
validation run.

r46 passed the short hydraulic comparison: it preserved the r43 physical
pressure, saturation, and `17.1 mm` connected wet depth to round-off while
removing all adjacent dry-layer gas-velocity sign changes. Its gas velocity
remained near `0.038 m/s`; the next accepted step is only a `0.1 s` extension
because the liquid-phase velocity had not yet plateaued.

That r47 extension showed the particle FLIP velocity mode returning after
`0.06 s` even though the reconstructed pressure gradients remained bounded.
The independent r48 option `reconstruct_darcy_velocity` reconstructs phase
velocities from the same permeability/viscosity Darcy relation used by the
semi-implicit pressure equation. It does not change solid `PIC=0`, physical
particle pressures, or the default solver path.

## Apply on HPC4

From the Siemens case directory in the `Backup-mpm` clone, submit the rebuild
through Slurm so that it receives the requested 32 CPUs:

```bash
sbatch program_patch/rebuild_hpc4.sh
```

The default source and executable are
`/home/xchenjm/chen/MPM/mpm` and `/home/xchenjm/chen/MPM/mpm/build/mpm`.
Set and export `MPM_SOURCE` before submission only if that source checkout
moves. The replacements are versioned in `Backup-mpm`, so the installer no
longer creates side-by-side source backups. After a successful build it
removes legacy `.pre-siemens-validation` and `.before-*` copies of the two
duplicated particle sources. It builds with 32 parallel jobs by default.

Do not run the full 400/900 s calculations immediately after rebuilding.
First run the short time-step matrix described in the case README; every
short job checks the initial suction and pressure ranges automatically.
