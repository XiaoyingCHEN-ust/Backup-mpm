"""Check VTP metadata without requiring the VTK Python package."""

from __future__ import annotations

import argparse
import math
import xml.etree.ElementTree as ET
from pathlib import Path


def point_ranges(path: Path):
    root = ET.parse(path).getroot()
    point_data = root.find("./PolyData/Piece/PointData")
    if point_data is None:
        raise ValueError(f"PointData is missing in {path}")
    ranges = {}
    for array in point_data.findall("DataArray"):
        name = array.get("Name")
        if name:
            ranges[name] = (float(array.get("RangeMin")), float(array.get("RangeMax")))
    return ranges


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("result_dir", type=Path)
    parser.add_argument("--expected-suction-pa", type=float, default=972.9920291627446)
    parser.add_argument("--pressure-limit-pa", type=float, default=1.0e6)
    args = parser.parse_args()

    files = sorted(args.result_dir.glob("particle*.vtp"))
    if not files:
        raise FileNotFoundError(f"No particle VTP files in {args.result_dir}")

    initial = point_ranges(files[0])
    observed_suction = initial["suction_pressures"][1]
    if not math.isclose(
        observed_suction, args.expected_suction_pa, rel_tol=0.02, abs_tol=1.0
    ):
        raise RuntimeError(
            "Executable does not contain the validation patch: initial maximum "
            f"suction is {observed_suction:.6g} Pa, expected "
            f"{args.expected_suction_pa:.6g} Pa"
        )

    maximum_pressure = 0.0
    saturation_min = math.inf
    saturation_max = -math.inf
    for path in files:
        ranges = point_ranges(path)
        for name in ("gas_pressures", "liquid_pressures"):
            low, high = ranges[name]
            maximum_pressure = max(maximum_pressure, abs(low), abs(high))
        low, high = ranges["liquid_saturations"]
        saturation_min = min(saturation_min, low)
        saturation_max = max(saturation_max, high)

    if not (0.0 <= saturation_min <= saturation_max <= 1.0):
        raise RuntimeError(
            f"Nonphysical saturation range: [{saturation_min}, {saturation_max}]"
        )
    if not math.isfinite(maximum_pressure) or maximum_pressure >= args.pressure_limit_pa:
        raise RuntimeError(
            f"Pressure stability check failed: max |p|={maximum_pressure:.6g} Pa"
        )

    marker = args.result_dir / "RANGE_CHECK_PASSED.txt"
    marker.write_text(
        f"files={len(files)}\n"
        f"initial_max_suction_pa={observed_suction:.12g}\n"
        f"max_abs_pressure_pa={maximum_pressure:.12g}\n"
        f"saturation_range={saturation_min:.12g},{saturation_max:.12g}\n",
        encoding="utf-8",
    )
    print(f"PASS: {args.result_dir} (max |p|={maximum_pressure:.6g} Pa)")


if __name__ == "__main__":
    main()
