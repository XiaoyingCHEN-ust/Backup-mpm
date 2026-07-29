"""Test whether permeability/time scaling reproduces an unscaled MPM state."""

from __future__ import annotations

import argparse
import math
import re
from pathlib import Path

from check_vtp_ranges import point_arrays, point_coordinates


STEP_PATTERN = re.compile(r"particle(\d+)\.vtp$")
ARRAY_RULES = {
    "liquid_saturations": (1.0e-3, 0.0, 1.0),
    "gas_saturations": (1.0e-3, 0.0, 1.0),
    "gas_pressures": (20.0, 1.0e-2, 1.0),
    "liquid_pressures": (20.0, 1.0e-2, 1.0),
    "suction_pressures": (20.0, 1.0e-2, 1.0),
    "liquid_permeabilities": (1.0e-18, 2.0e-2, "permeability"),
    "gas_permeabilities": (1.0e-18, 2.0e-2, "permeability"),
    "liquid_velocities": (5.0e-4, 5.0e-2, "velocity"),
    "gas_velocities": (5.0e-4, 5.0e-2, "velocity"),
}


def particle_file(result_dir: Path, requested_step: int | None) -> tuple[Path, int]:
    candidates = []
    for path in result_dir.glob("particle*.vtp"):
        match = STEP_PATTERN.fullmatch(path.name)
        if match is not None:
            candidates.append((int(match.group(1)), path))
    if not candidates:
        raise FileNotFoundError(f"No particle VTP files in {result_dir}")
    if requested_step is None:
        return max(candidates, key=lambda item: item[0])[1], max(
            step for step, _ in candidates
        )
    matches = [path for step, path in candidates if step == requested_step]
    if len(matches) != 1:
        raise FileNotFoundError(
            f"Expected one particle output at step {requested_step} in "
            f"{result_dir}, found {len(matches)}"
        )
    return matches[0], requested_step


def metrics(reference, scaled):
    if len(reference) != len(scaled):
        raise RuntimeError(
            f"Array sizes differ: reference={len(reference)}, scaled={len(scaled)}"
        )
    differences = [abs(a - b) for a, b in zip(reference, scaled)]
    maximum = max(differences, default=0.0)
    rms = math.sqrt(sum(value * value for value in differences) / len(differences))
    scale = max((abs(value) for value in reference), default=0.0)
    return maximum, rms, scale


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("reference_dir", type=Path)
    parser.add_argument("scaled_dir", type=Path)
    parser.add_argument("--permeability-scale", type=float, required=True)
    parser.add_argument("--reference-step", type=int)
    parser.add_argument("--scaled-step", type=int)
    args = parser.parse_args()
    if args.permeability_scale <= 0.0:
        parser.error("--permeability-scale must be positive")

    reference_path, reference_step = particle_file(
        args.reference_dir, args.reference_step
    )
    scaled_path, scaled_step = particle_file(args.scaled_dir, args.scaled_step)
    names = tuple(ARRAY_RULES)
    reference = point_arrays(reference_path, names)
    scaled = point_arrays(scaled_path, names)

    coordinate_metrics = metrics(
        tuple(value for xyz in point_coordinates(reference_path) for value in xyz),
        tuple(value for xyz in point_coordinates(scaled_path) for value in xyz),
    )
    failures = []
    if coordinate_metrics[0] > 1.0e-10:
        failures.append(
            f"coordinates: max={coordinate_metrics[0]:.6g} m, limit=1e-10 m"
        )

    lines = [
        f"reference={reference_path}",
        f"scaled={scaled_path}",
        f"reference_step={reference_step}",
        f"scaled_step={scaled_step}",
        f"permeability_scale={args.permeability_scale:.12g}",
        f"coordinate_max_abs_difference_m={coordinate_metrics[0]:.12g}",
    ]
    for name, (absolute_tolerance, relative_tolerance, scaling) in ARRAY_RULES.items():
        scaled_values = scaled[name]
        if scaling in ("permeability", "velocity"):
            scaled_values = tuple(
                value / args.permeability_scale for value in scaled_values
            )
        maximum, rms, field_scale = metrics(reference[name], scaled_values)
        limit = absolute_tolerance + relative_tolerance * field_scale
        lines.extend(
            (
                f"{name}_max_abs_difference={maximum:.12g}",
                f"{name}_rms_difference={rms:.12g}",
                f"{name}_limit={limit:.12g}",
            )
        )
        if not math.isfinite(maximum) or maximum > limit:
            failures.append(f"{name}: max={maximum:.6g}, limit={limit:.6g}")

    lines.insert(0, f"status={'FAIL' if failures else 'PASS'}")
    report = args.scaled_dir / "TIME_SCALING_COMPARISON.txt"
    report.write_text("\n".join(lines) + "\n", encoding="utf-8")
    if failures:
        raise RuntimeError(
            "Permeability/time scaling comparison failed: " + "; ".join(failures)
        )
    print(f"PASS: permeability/time scaling is equivalent ({report})")


if __name__ == "__main__":
    main()
