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
history converged at `dt=1e-4 s`. The r33 three-second closed test must pass
before the first experimental-time checkpoint.

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
