# Backup-mpm: slope infiltration manuscript case

This repository contains the compact, reproducible subset needed to revise the
paper *The Infiltration Paradox: Pore-gas Pressure Controls on Early Slope
Failure during Rapid Hydraulic Loading*.

Included:

- the stable-slope input `mpm-3p-initial.json`;
- the resumed loading input `mpm-3p.json`;
- the mesh, particle and entity-set files required by those inputs;
- a sensitivity workflow for permeability and gas boundary conditions;
- two replacement C++ files that remove the former hard-coded pressure
  choices.

Excluded:

- approximately 16 GB of VTP/HDF5 result files;
- local build products and toolchains;
- the separate full `mpm` source checkout.

Start with
[`Example/slope_infiltration/sensitivity/README.md`](Example/slope_infiltration/sensitivity/README.md).
The baseline adopted for the new calculations is
\(k=1.5\times10^{-12}\ {\rm m^2}\) and
\(\Delta t=1.0\times10^{-4}\ {\rm s}\).
