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

The Siemens inputs use `initial_suction_from_swrc` because the published
initial state is close to residual saturation and is incompatible with a
hard-coded 20 psi value. Existing inputs retain 20 psi by default.

## Apply on HPC4

From the Siemens case directory in the `Backup-mpm` clone, submit the rebuild
through Slurm so that it receives the requested 32 CPUs:

```bash
sbatch program_patch/rebuild_hpc4.sh
```

The default source and executable are
`/home/xchenjm/chen/MPM/mpm` and `/home/xchenjm/chen/MPM/mpm/build/mpm`.
Set and export `MPM_SOURCE` before submission only if that source checkout
moves. The installer creates one `.pre-siemens-validation` backup beside each
original source file before replacement and builds with 32 parallel jobs by
default.

Do not run the full 400/900 s calculations immediately after rebuilding.
First run the short time-step matrix described in the case README; every
short job checks the initial suction and pressure ranges automatically.
