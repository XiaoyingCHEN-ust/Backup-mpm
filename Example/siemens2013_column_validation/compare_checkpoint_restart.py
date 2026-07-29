"""Compare the final VTP state of a resumed run with a continuous reference."""

from __future__ import annotations

import argparse
import math
from pathlib import Path

from check_vtp_ranges import point_arrays, point_coordinates


ARRAY_TOLERANCES = {
    "liquid_saturations": (1.0e-7, 1.0e-7),
    "gas_saturations": (1.0e-7, 1.0e-7),
    "gas_pressures": (5.0e-2, 1.0e-6),
    "liquid_pressures": (5.0e-2, 1.0e-6),
    "suction_pressures": (5.0e-2, 1.0e-6),
    "liquid_permeabilities": (1.0e-18, 1.0e-5),
    "gas_permeabilities": (1.0e-18, 1.0e-5),
    "liquid_velocities": (1.0e-7, 1.0e-5),
    "gas_velocities": (1.0e-7, 1.0e-5),
}


def comparison_metrics(reference, resumed):
    if len(reference) != len(resumed):
        raise RuntimeError(
            f"Array sizes differ: reference={len(reference)}, resumed={len(resumed)}"
        )
    differences = [abs(a - b) for a, b in zip(reference, resumed)]
    maximum = max(differences, default=0.0)
    rms = math.sqrt(sum(value * value for value in differences) / len(differences))
    scale = max((abs(value) for value in reference), default=0.0)
    return maximum, rms, scale


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("continuous_dir", type=Path)
    parser.add_argument("resumed_dir", type=Path)
    parser.add_argument("--step", type=int, default=1_200_000)
    args = parser.parse_args()

    filename = f"particle{args.step:07d}.vtp"
    continuous_path = args.continuous_dir / filename
    resumed_path = args.resumed_dir / filename
    if not continuous_path.is_file():
        parser.error(f"continuous reference is missing: {continuous_path}")
    if not resumed_path.is_file():
        parser.error(f"resumed result is missing: {resumed_path}")

    marker = args.resumed_dir / "CHECKPOINT_COMPARISON_PASSED.txt"
    marker.unlink(missing_ok=True)
    names = tuple(ARRAY_TOLERANCES)
    continuous = point_arrays(continuous_path, names)
    resumed = point_arrays(resumed_path, names)

    coordinate_metrics = comparison_metrics(
        tuple(value for xyz in point_coordinates(continuous_path) for value in xyz),
        tuple(value for xyz in point_coordinates(resumed_path) for value in xyz),
    )
    if coordinate_metrics[0] > 1.0e-10:
        raise RuntimeError(
            "Checkpoint coordinate mismatch: "
            f"max difference={coordinate_metrics[0]:.6g} m"
        )

    lines = [
        f"continuous={continuous_path}",
        f"resumed={resumed_path}",
        f"coordinate_max_abs_difference_m={coordinate_metrics[0]:.12g}",
    ]
    failures = []
    for name, (absolute_tolerance, relative_tolerance) in ARRAY_TOLERANCES.items():
        maximum, rms, scale = comparison_metrics(continuous[name], resumed[name])
        limit = absolute_tolerance + relative_tolerance * scale
        lines.append(f"{name}_max_abs_difference={maximum:.12g}")
        lines.append(f"{name}_rms_difference={rms:.12g}")
        lines.append(f"{name}_limit={limit:.12g}")
        if not math.isfinite(maximum) or maximum > limit:
            failures.append(f"{name}: max={maximum:.6g}, limit={limit:.6g}")

    if failures:
        raise RuntimeError("Checkpoint comparison failed: " + "; ".join(failures))

    marker.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"PASS: checkpoint result matches continuous reference ({marker})")


if __name__ == "__main__":
    main()
