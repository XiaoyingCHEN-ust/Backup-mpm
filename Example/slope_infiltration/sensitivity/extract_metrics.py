#!/usr/bin/env python3
"""Reduce particle VTP files to a small, transferable metrics CSV."""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy


STEP_PATTERN = re.compile(r"(\d+)\.vtp$")


def read_array(point_data, name: str) -> np.ndarray:
    array = point_data.GetArray(name)
    if array is None:
        raise KeyError(f"VTK array {name!r} was not found")
    return vtk_to_numpy(array)


def read_metrics(path: Path) -> dict[str, float]:
    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(str(path))
    reader.Update()
    data = reader.GetOutput()
    if data.GetNumberOfPoints() == 0:
        raise ValueError(f"No particle points were read from {path}")

    point_data = data.GetPointData()
    displacement = read_array(point_data, "displacements")[:, :2]
    velocity = read_array(point_data, "velocities")[:, :2]
    strain = read_array(point_data, "strains")
    gas_pressure = read_array(point_data, "gas_pressures")
    liquid_pressure = read_array(point_data, "liquid_pressures")

    displacement_magnitude = np.linalg.norm(displacement, axis=1)
    velocity_magnitude = np.linalg.norm(velocity, axis=1)
    strain_magnitude = np.linalg.norm(strain, axis=1)

    return {
        "p95_displacement_m": np.percentile(displacement_magnitude, 95),
        "p99_displacement_m": np.percentile(displacement_magnitude, 99),
        "max_displacement_m": np.max(displacement_magnitude),
        "p95_velocity_m_s": np.percentile(velocity_magnitude, 95),
        "p99_velocity_m_s": np.percentile(velocity_magnitude, 99),
        "max_velocity_m_s": np.max(velocity_magnitude),
        "p99_strain_magnitude": np.percentile(strain_magnitude, 99),
        "max_strain_magnitude": np.max(strain_magnitude),
        "p01_gas_pressure_pa": np.percentile(gas_pressure, 1),
        "median_gas_pressure_pa": np.median(gas_pressure),
        "p99_gas_pressure_pa": np.percentile(gas_pressure, 99),
        "p99_liquid_pressure_pa": np.percentile(liquid_pressure, 99),
    }


def legacy_result_dir(input_dir: Path, row: dict[str, str]) -> Path | None:
    """Map reference cases to the existing manuscript result directories."""
    if float(row["permeability_factor"]) != 1.0:
        return None
    head_pa = int(row["head_pa"])
    if row["gas_model"] == "fixed":
        return input_dir / "results" / f"stability3p-pl={head_pa}-new"
    if (
        row["gas_model"] == "variable"
        and float(row["surface_gas_pressure_ratio"]) == 1.0
    ):
        return input_dir / "results" / f"stability3p-gas-pl={head_pa}-new"
    return None


def main() -> None:
    script_dir = Path(__file__).resolve().parent
    input_dir = script_dir.parent

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--reuse-legacy-reference",
        action="store_true",
        help="Use existing k/k0=1 manuscript result directories when needed",
    )
    parser.add_argument(
        "--family",
        choices=("all", "permeability", "gas_boundary"),
        default="all",
    )
    args = parser.parse_args()

    with (script_dir / "manifest.csv").open(
        newline="", encoding="utf-8"
    ) as stream:
        manifest = list(csv.DictReader(stream))
    if args.family != "all":
        manifest = [
            row for row in manifest if args.family in row["families"].split("+")
        ]

    rows: list[dict] = []
    missing: list[str] = []
    for case in manifest:
        result_dir = input_dir / "results" / case["uuid"]
        if not result_dir.is_dir() and args.reuse_legacy_reference:
            fallback = legacy_result_dir(input_dir, case)
            if fallback is not None:
                result_dir = fallback
        if not result_dir.is_dir():
            missing.append(case["case_id"])
            continue

        files = sorted(result_dir.glob("particle*.vtp"))
        if not files:
            missing.append(case["case_id"])
            continue

        print(f"{case['case_id']}: {len(files)} VTP files", flush=True)
        for path in files:
            match = STEP_PATTERN.search(path.name)
            if match is None:
                continue
            step = int(match.group(1))
            metrics = read_metrics(path)
            rows.append(
                {
                    "case_id": case["case_id"],
                    "step": step,
                    "time_s": step * float(case["dt_s"]),
                    "source_result_dir": result_dir.name,
                    **metrics,
                }
            )

    output_path = script_dir / "metrics.csv"
    if not rows:
        raise RuntimeError("No metrics were extracted")
    with output_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote {len(rows)} rows to {output_path}")
    if missing:
        print(f"Missing {len(missing)} case(s):")
        for case_id in missing:
            print(f"  {case_id}")


if __name__ == "__main__":
    main()
