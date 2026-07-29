"""Compare semi-implicit histories with explicit or finer-step baselines."""

from __future__ import annotations

import argparse
import csv
import math
import re
from pathlib import Path

from check_vtp_ranges import point_arrays, point_coordinates


STEP_PATTERN = re.compile(r"particle(\d+)\.vtp$")
FIELDS = (
    "liquid_saturations",
    "gas_saturations",
    "gas_pressures",
    "liquid_pressures",
    "suction_pressures",
    "liquid_permeabilities",
    "gas_permeabilities",
    "liquid_velocities",
    "gas_velocities",
)


def output_files(directory: Path, dt: float) -> dict[float, Path]:
    outputs: dict[float, Path] = {}
    for path in directory.glob("particle*.vtp"):
        match = STEP_PATTERN.fullmatch(path.name)
        if match:
            outputs[round(int(match.group(1)) * dt, 12)] = path
    if not outputs:
        raise FileNotFoundError(f"No particle VTP files in {directory}")
    return outputs


def error_metrics(reference, candidate) -> tuple[float, float]:
    if len(reference) != len(candidate):
        raise RuntimeError("Compared arrays have different sizes")
    differences = [abs(a - b) for a, b in zip(reference, candidate)]
    return max(differences, default=0.0), math.sqrt(
        sum(value * value for value in differences) / len(differences)
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--manifest", type=Path, default=Path("semi_implicit_inputs/manifest.csv")
    )
    parser.add_argument(
        "--results", type=Path, default=Path("stability_results")
    )
    parser.add_argument(
        "--case",
        choices=("open", "closed", "both"),
        default="both",
        help="Compare only the selected boundary condition (default: both)",
    )
    parser.add_argument(
        "--open-reference-dir",
        type=Path,
        default=None,
        help="reuse an existing open explicit result directory",
    )
    parser.add_argument(
        "--open-reference-dt",
        type=float,
        default=2.5e-6,
        help="time step used by --open-reference-dir",
    )
    parser.add_argument(
        "--semi-implicit-reference-dt",
        type=float,
        default=None,
        help=(
            "compare the other semi-implicit rows with the generated "
            "semi-implicit row at this time step"
        ),
    )
    args = parser.parse_args()

    with args.manifest.open(encoding="utf-8", newline="") as stream:
        rows = list(csv.DictReader(stream))
    case_names = ("open", "closed") if args.case == "both" else (args.case,)
    for case_name in case_names:
        case_rows = [row for row in rows if row["case"] == case_name]
        semi_implicit_rows = [
            row for row in case_rows if row["method"] == "semi_implicit"
        ]
        if args.semi_implicit_reference_dt is not None:
            matching_rows = [
                row
                for row in semi_implicit_rows
                if math.isclose(
                    float(row["dt_s"]),
                    args.semi_implicit_reference_dt,
                    rel_tol=1.0e-12,
                    abs_tol=0.0,
                )
            ]
            if len(matching_rows) != 1:
                raise RuntimeError(
                    f"Expected one {case_name} semi-implicit reference at "
                    f"dt={args.semi_implicit_reference_dt:g} s; "
                    f"found {len(matching_rows)}"
                )
            reference_row = matching_rows[0]
            reference_files = output_files(
                args.results / reference_row["uuid"],
                float(reference_row["dt_s"]),
            )
            rows_to_compare = [
                row for row in semi_implicit_rows if row is not reference_row
            ]
        else:
            explicit_rows = [
                row for row in case_rows if row["method"] == "explicit"
            ]
            if len(explicit_rows) != 1:
                raise RuntimeError(f"Expected one {case_name} explicit baseline")
            reference_row = explicit_rows[0]
            if case_name == "open" and args.open_reference_dir is not None:
                reference_files = output_files(
                    args.open_reference_dir, args.open_reference_dt
                )
            else:
                reference_files = output_files(
                    args.results / reference_row["uuid"],
                    float(reference_row["dt_s"]),
                )
            rows_to_compare = semi_implicit_rows

        for row in rows_to_compare:
            candidate_files = output_files(
                args.results / row["uuid"], float(row["dt_s"])
            )
            if not set(candidate_files).issubset(reference_files):
                raise RuntimeError(
                    f"Reference lacks candidate output times for {row['uuid']}: "
                    f"missing={sorted(set(candidate_files) - set(reference_files))}"
                )

            totals = {
                name: {"maximum": 0.0, "squares": 0.0, "count": 0}
                for name in FIELDS
            }
            coordinate_maximum = 0.0
            for time in sorted(candidate_files):
                reference_path = reference_files[time]
                candidate_path = candidate_files[time]
                reference = point_arrays(reference_path, FIELDS)
                candidate = point_arrays(candidate_path, FIELDS)
                reference_coordinates = tuple(
                    value
                    for xyz in point_coordinates(reference_path)
                    for value in xyz
                )
                candidate_coordinates = tuple(
                    value
                    for xyz in point_coordinates(candidate_path)
                    for value in xyz
                )
                coordinate_maximum = max(
                    coordinate_maximum,
                    error_metrics(reference_coordinates, candidate_coordinates)[0],
                )
                for name in FIELDS:
                    maximum, _ = error_metrics(reference[name], candidate[name])
                    differences = [
                        a - b for a, b in zip(reference[name], candidate[name])
                    ]
                    totals[name]["maximum"] = max(
                        totals[name]["maximum"], maximum
                    )
                    totals[name]["squares"] += sum(
                        value * value for value in differences
                    )
                    totals[name]["count"] += len(differences)

            print(
                f"{case_name} dt={float(row['dt_s']):g} s, "
                f"penalty={float(row['boundary_penalty']):g}, "
                f"reference={reference_row['method']} "
                f"dt={float(reference_row['dt_s']):g} s, "
                f"outputs={len(candidate_files)}, "
                f"coordinate_max={coordinate_maximum:.6g} m"
            )
            for name in FIELDS:
                total = totals[name]
                rms = math.sqrt(total["squares"] / total["count"])
                print(
                    f"  {name}: max={total['maximum']:.6g}, rms={rms:.6g}"
                )


if __name__ == "__main__":
    main()
