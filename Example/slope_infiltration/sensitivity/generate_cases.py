#!/usr/bin/env python3
"""Generate the paired sensitivity cases used for the reviewer response."""

from __future__ import annotations

import argparse
import copy
import csv
import json
from pathlib import Path


REFERENCE_PERMEABILITY = 1.5e-12
REFERENCE_DT = 1.0e-4
HEAD_PRESSURES_PA = (0, 5000, 10000, 15000, 20000, 25000)
PERMEABILITY_FACTORS = (0.1, 1.0, 10.0)
SURFACE_GAS_PRESSURE_RATIOS = (0.0, 0.5, 1.0)


def number_token(value: float) -> str:
    return f"{value:g}".replace(".", "p").replace("-", "m")


def find_material(config: dict, material_id: int) -> dict:
    for material in config["materials"]:
        if material["id"] == material_id:
            return material
    raise ValueError(f"Material id {material_id} was not found")


def register_case(
    cases: dict,
    *,
    family: str,
    permeability_factor: float,
    gas_model: str,
    surface_gas_pressure_ratio: float,
    head_pa: int,
) -> None:
    if gas_model not in {"variable", "fixed"}:
        raise ValueError(f"Unknown gas model: {gas_model}")

    ratio_token = (
        f"g{number_token(surface_gas_pressure_ratio)}"
        if gas_model == "variable"
        else "gfixed"
    )
    case_id = (
        f"k{number_token(permeability_factor)}-{ratio_token}-"
        f"{gas_model}-h{head_pa:05d}"
    )
    if case_id not in cases:
        cases[case_id] = {
            "case_id": case_id,
            "families": set(),
            "permeability_factor": permeability_factor,
            "gas_model": gas_model,
            "surface_gas_pressure_ratio": surface_gas_pressure_ratio,
            "head_pa": head_pa,
        }
    cases[case_id]["families"].add(family)


def build_case_matrix() -> list[dict]:
    cases: dict[str, dict] = {}

    # Figure panel (a): one-at-a-time intrinsic-permeability sensitivity.
    for factor in PERMEABILITY_FACTORS:
        for head_pa in HEAD_PRESSURES_PA:
            register_case(
                cases,
                family="permeability",
                permeability_factor=factor,
                gas_model="variable",
                surface_gas_pressure_ratio=1.0,
                head_pa=head_pa,
            )
            register_case(
                cases,
                family="permeability",
                permeability_factor=factor,
                gas_model="fixed",
                surface_gas_pressure_ratio=0.0,
                head_pa=head_pa,
            )

    # Figure panel (b): gas-boundary sensitivity at the reference permeability.
    for ratio in SURFACE_GAS_PRESSURE_RATIOS:
        for head_pa in HEAD_PRESSURES_PA:
            register_case(
                cases,
                family="gas_boundary",
                permeability_factor=1.0,
                gas_model="variable",
                surface_gas_pressure_ratio=ratio,
                head_pa=head_pa,
            )
            register_case(
                cases,
                family="gas_boundary",
                permeability_factor=1.0,
                gas_model="fixed",
                surface_gas_pressure_ratio=0.0,
                head_pa=head_pa,
            )

    return [cases[key] for key in sorted(cases)]


def write_case(base: dict, case: dict, case_dir: Path) -> dict:
    config = copy.deepcopy(base)
    soil = find_material(config, 0)
    fluid = find_material(config, 1)

    permeability = REFERENCE_PERMEABILITY * case["permeability_factor"]
    soil["intrinsic_permeability"] = permeability
    fluid["fixed_gas_pressure"] = case["gas_model"] == "fixed"
    fluid["surface_liquid_pressure"] = float(case["head_pa"])
    fluid["surface_gas_pressure_ratio"] = float(
        case["surface_gas_pressure_ratio"]
    )

    analysis = config["analysis"]
    analysis["uuid"] = f"sens-{case['case_id']}"
    analysis["dt"] = REFERENCE_DT
    analysis["nsteps"] = 385000

    # Retain only the fields needed to determine onset and audit the gas response.
    output = config["post_processing"]
    output["write_hdf5"] = False
    output["write_vtk"] = True
    output["output_steps"] = 5000
    output["vtk"] = [
        "porosities",
        "displacements",
        "velocities",
        "strains",
        "free_surfaces",
    ]
    output["liquid_vtk"] = [
        "liquid_pressures",
        "liquid_saturations",
        "gas_pressures",
    ]

    case_path = case_dir / f"{case['case_id']}.json"
    case_path.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")

    return {
        "case_id": case["case_id"],
        "case_file": case_path.as_posix(),
        "uuid": analysis["uuid"],
        "families": "+".join(sorted(case["families"])),
        "permeability_m2": f"{permeability:.8e}",
        "permeability_factor": f"{case['permeability_factor']:g}",
        "gas_model": case["gas_model"],
        "surface_gas_pressure_ratio": f"{case['surface_gas_pressure_ratio']:g}",
        "head_pa": case["head_pa"],
        "head_m": f"{case['head_pa'] / 10000.0:g}",
        "dt_s": f"{REFERENCE_DT:.8e}",
        "nsteps": analysis["nsteps"],
    }


def main() -> None:
    script_dir = Path(__file__).resolve().parent
    input_dir = script_dir.parent

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--base",
        type=Path,
        default=input_dir / "mpm-3p.json",
        help="Baseline calculation JSON",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=script_dir / "cases",
        help="Directory for generated JSON files",
    )
    args = parser.parse_args()

    base = json.loads(args.base.read_text(encoding="utf-8"))
    soil = find_material(base, 0)
    if float(soil["intrinsic_permeability"]) != REFERENCE_PERMEABILITY:
        raise ValueError(
            "Baseline intrinsic_permeability must be 1.5e-12 m^2; "
            f"found {soil['intrinsic_permeability']!r}"
        )
    if float(base["analysis"]["dt"]) != REFERENCE_DT:
        raise ValueError(
            "Baseline dt must be 1.0e-4 s; "
            f"found {base['analysis']['dt']!r}"
        )

    args.output_dir.mkdir(parents=True, exist_ok=True)
    rows = [write_case(base, case, args.output_dir) for case in build_case_matrix()]

    manifest_path = script_dir / "manifest.csv"
    with manifest_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    print(f"Generated {len(rows)} unique cases in {args.output_dir}")
    print(f"Manifest: {manifest_path}")


if __name__ == "__main__":
    main()
