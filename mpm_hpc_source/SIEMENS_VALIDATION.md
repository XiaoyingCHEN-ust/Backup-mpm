# Siemens column validation source snapshot

This directory is a clean, buildable snapshot of the MPM source used for the
Siemens et al. (2013) three-phase column validation. It contains the accepted
semi-implicit pressure implementation, pressure projection, checkpoint-resume
fix, and force-pressure VTK output support.

The machine-specific `build/` directory, previous Git object database,
temporary inputs, and numerical results are intentionally excluded. Configure
and rebuild them on the target Linux machine:

```bash
cd mpm_hpc
MPM_BUILD_JOBS=32 bash build_linux_validation.sh
```

The validation case, preparation scripts, reference data, and result checks
are stored in `../Example/siemens2013_column_validation/`. The independently
tracked recovery copies of the modified source files remain in
`../server_overrides/`.
