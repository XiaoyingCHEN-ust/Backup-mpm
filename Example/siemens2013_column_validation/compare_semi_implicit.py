"""Compare semi-implicit smoke-test histories with their explicit baselines."""

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
    args = parser.parse_args()

    with args.manifest.open(encoding="utf-8", newline="") as stream:
        rows = list(csv.DictReader(stream))
    for case_name in ("open", "closed"):
        case_rows = [row for row in rows if row["case"] == case_name]
        explicit_rows = [row for row in case_rows if row["method"] == "explicit"]
        if len(explicit_rows) != 1:
            raise RuntimeError(f"Expected one {case_name} explicit baseline")
        reference_row = explicit_rows[0]
        reference_files = output_files(
            args.results / reference_row["uuid"], float(reference_row["dt_s"])
        )

        for row in case_rows:
            if row["method"] != "semi_implicit":
                continue
            candidate_files = output_files(
                args.results / row["uuid"], float(row["dt_s"])
            )
            if set(candidate_files) != set(reference_files):
                raise RuntimeError(
                    f"Physical output times differ for {row['uuid']}: "
                    f"reference={sorted(reference_files)}, "
                    f"candidate={sorted(candidate_files)}"
                )

            totals = {
                name: {"maximum": 0.0, "squares": 0.0, "count": 0}
                for name in FIELDS
            }
            coordinate_maximum = 0.0
            for time in sorted(reference_files):
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
