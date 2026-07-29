#!/usr/bin/env python3
"""Allow the opt-in three-phase force-pressure fields in liquid VTK output."""

from __future__ import annotations

import argparse
from pathlib import Path


ALLOWED_ANCHOR = '  std::vector<std::string> liquid_vtk_allowed = liquid_vtk;\n'
ALLOWED_INSERTION = ALLOWED_ANCHOR + (
    "  // Force-only pressures are opt-in three-phase diagnostics.  Keep them out\n"
    "  // of the default list because other particle types do not register them.\n"
    '  liquid_vtk_allowed.emplace_back("force_liquid_pressures");\n'
    '  liquid_vtk_allowed.emplace_back("force_gas_pressures");\n'
)
ALLOWED_LIQUID_MARKER = (
    'liquid_vtk_allowed.emplace_back("force_liquid_pressures");'
)
ALLOWED_GAS_MARKER = (
    'liquid_vtk_allowed.emplace_back("force_gas_pressures");'
)
SCALAR_LIQUID_ANCHOR = '                                               "PIC_liquid_pressures",\n'
SCALAR_LIQUID_INSERTION = SCALAR_LIQUID_ANCHOR + (
    '                                               "force_liquid_pressures",\n'
)
SCALAR_GAS_ANCHOR = '                                              "PIC_gas_pressures",\n'
SCALAR_GAS_INSERTION = SCALAR_GAS_ANCHOR + (
    '                                              "force_gas_pressures",\n'
)
SCALAR_LIQUID_MARKER = '"force_liquid_pressures",'
SCALAR_GAS_MARKER = '"force_gas_pressures",'


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
    changed = False

    has_allowed_liquid = ALLOWED_LIQUID_MARKER in source
    has_allowed_gas = ALLOWED_GAS_MARKER in source
    if has_allowed_liquid != has_allowed_gas:
        raise RuntimeError(
            f"Partial force-pressure VTK whitelist patch in {target}; refusing to edit"
        )
    if not has_allowed_liquid:
        if source.count(ALLOWED_ANCHOR) != 1:
            raise RuntimeError(
                f"Unrecognised liquid VTK whitelist layout in {target}; refusing to edit"
            )
        source = source.replace(ALLOWED_ANCHOR, ALLOWED_INSERTION)
        changed = True

    has_scalar_liquid = SCALAR_LIQUID_MARKER in source
    has_scalar_gas = SCALAR_GAS_MARKER in source
    if has_scalar_liquid != has_scalar_gas:
        raise RuntimeError(
            f"Partial force-pressure VTK writer patch in {target}; refusing to edit"
        )
    if not has_scalar_liquid:
        if source.count(SCALAR_LIQUID_ANCHOR) != 1:
            raise RuntimeError(
                f"Unrecognised liquid scalar-writer layout in {target}; refusing to edit"
            )
        if source.count(SCALAR_GAS_ANCHOR) != 1:
            raise RuntimeError(
                f"Unrecognised gas scalar-writer layout in {target}; refusing to edit"
            )
        source = source.replace(SCALAR_LIQUID_ANCHOR, SCALAR_LIQUID_INSERTION)
        source = source.replace(SCALAR_GAS_ANCHOR, SCALAR_GAS_INSERTION)
        changed = True

    if changed:
        target.write_text(source, encoding="utf-8")
        print(f"Patched: validate and write force-pressure VTK fields ({target})")
    else:
        print(f"Already patched: validate and write force-pressure VTK fields ({target})")


if __name__ == "__main__":
    main()
