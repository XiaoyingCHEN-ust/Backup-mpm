#!/usr/bin/env python3
"""Allow the opt-in three-phase force-pressure fields in liquid VTK output."""

from __future__ import annotations

import argparse
from pathlib import Path


ANCHOR = '  std::vector<std::string> liquid_vtk_allowed = liquid_vtk;\n'
INSERTION = ANCHOR + (
    "  // Force-only pressures are opt-in three-phase diagnostics.  Keep them out\n"
    "  // of the default list because other particle types do not register them.\n"
    '  liquid_vtk_allowed.emplace_back("force_liquid_pressures");\n'
    '  liquid_vtk_allowed.emplace_back("force_gas_pressures");\n'
)
LIQUID_MARKER = (
    'liquid_vtk_allowed.emplace_back("force_liquid_pressures");'
)
GAS_MARKER = 'liquid_vtk_allowed.emplace_back("force_gas_pressures");'


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--mpm-source",
        type=Path,
        required=True,
        help="Root of the MPM source tree",
    )
    args = parser.parse_args()

    target = args.mpm_source / "include" / "solvers" / "mpm_base.tcc"
    if not target.is_file():
        raise FileNotFoundError(f"MPM base implementation is missing: {target}")

    source = target.read_text(encoding="utf-8")
    has_liquid = LIQUID_MARKER in source
    has_gas = GAS_MARKER in source
    if has_liquid and has_gas:
        print(f"Already patched: allow force-pressure VTK fields ({target})")
        return
    if has_liquid != has_gas:
        raise RuntimeError(
            f"Partial force-pressure VTK patch detected in {target}; refusing to edit"
        )
    if source.count(ANCHOR) != 1:
        raise RuntimeError(
            f"Unrecognised liquid VTK whitelist layout in {target}; refusing to edit"
        )

    target.write_text(source.replace(ANCHOR, INSERTION), encoding="utf-8")
    print(f"Patched: allow force-pressure VTK fields ({target})")


if __name__ == "__main__":
    main()
