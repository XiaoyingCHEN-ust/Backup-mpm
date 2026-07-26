"""Check actual VTP arrays without requiring the VTK Python package."""

from __future__ import annotations

import argparse
import base64
import math
import struct
import xml.etree.ElementTree as ET
import zlib
from pathlib import Path


VTK_STRUCT_TYPES = {"Float32": "f", "Float64": "d"}


def decode_binary_array(element, byte_order: str, header_type: str, compressor: str):
    payload = "".join(element.itertext()).strip()
    endian = "<" if byte_order == "LittleEndian" else ">"
    header_code = "I" if header_type == "UInt32" else "Q"
    header_size = struct.calcsize(header_code)

    if compressor != "vtkZLibDataCompressor":
        raise ValueError(f"Unsupported VTP compressor: {compressor}")

    probe_chars = 4 * math.ceil(header_size / 3)
    probe = base64.b64decode(payload[:probe_chars])
    number_of_blocks = struct.unpack(endian + header_code, probe[:header_size])[0]
    header_bytes = header_size * (3 + number_of_blocks)
    header_chars = 4 * math.ceil(header_bytes / 3)
    header = base64.b64decode(payload[:header_chars])
    fields = struct.unpack(
        endian + header_code * (3 + number_of_blocks), header[:header_bytes]
    )
    block_size, last_block_size = fields[1], fields[2]
    compressed_sizes = fields[3:]
    compressed = base64.b64decode(payload[header_chars:])

    raw_blocks = []
    offset = 0
    for index, compressed_size in enumerate(compressed_sizes):
        block = zlib.decompress(compressed[offset : offset + compressed_size])
        expected = last_block_size if index == number_of_blocks - 1 else block_size
        if len(block) != expected:
            raise ValueError(
                f"Decoded VTP block has {len(block)} bytes; expected {expected}"
            )
        raw_blocks.append(block)
        offset += compressed_size
    return b"".join(raw_blocks), endian


def point_arrays(path: Path, requested_names):
    root = ET.parse(path).getroot()
    point_data = root.find("./PolyData/Piece/PointData")
    if point_data is None:
        raise ValueError(f"PointData is missing in {path}")
    arrays = {}
    for element in point_data.findall("DataArray"):
        name = element.get("Name")
        if name not in requested_names:
            continue
        if element.get("format") != "binary":
            raise ValueError(f"Unsupported VTP array format for {name}")
        vtk_type = element.get("type")
        if vtk_type not in VTK_STRUCT_TYPES:
            raise ValueError(f"Unsupported VTP data type for {name}: {vtk_type}")
        raw, endian = decode_binary_array(
            element,
            root.get("byte_order", "LittleEndian"),
            root.get("header_type", "UInt32"),
            root.get("compressor", ""),
        )
        code = VTK_STRUCT_TYPES[vtk_type]
        value_size = struct.calcsize(code)
        if len(raw) % value_size:
            raise ValueError(f"Decoded byte count is invalid for {name} in {path}")
        arrays[name] = tuple(
            value[0] for value in struct.iter_unpack(endian + code, raw)
        )
    missing = set(requested_names) - set(arrays)
    if missing:
        raise KeyError(f"Missing VTP arrays in {path}: {sorted(missing)}")
    return arrays


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("result_dir", type=Path)
    parser.add_argument("--expected-suction-pa", type=float, default=972.9920291627446)
    parser.add_argument("--pressure-limit-pa", type=float, default=1.0e6)
    args = parser.parse_args()

    files = sorted(args.result_dir.glob("particle*.vtp"))
    if not files:
        raise FileNotFoundError(f"No particle VTP files in {args.result_dir}")

    marker = args.result_dir / "RANGE_CHECK_PASSED.txt"
    marker.unlink(missing_ok=True)
    names = (
        "liquid_saturations",
        "gas_pressures",
        "liquid_pressures",
        "suction_pressures",
    )
    initial = point_arrays(files[0], names)
    if any(not math.isfinite(value) for value in initial["suction_pressures"]):
        raise RuntimeError("Initial suction array contains non-finite values")
    observed_suction = max(initial["suction_pressures"])
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
        arrays = point_arrays(path, names)
        for name, values in arrays.items():
            nonfinite = sum(not math.isfinite(value) for value in values)
            if nonfinite:
                raise RuntimeError(
                    f"Non-finite {name} values in {path.name}: {nonfinite}"
                )
        for name in ("gas_pressures", "liquid_pressures"):
            maximum_pressure = max(
                maximum_pressure, max(abs(value) for value in arrays[name])
            )
        low = min(arrays["liquid_saturations"])
        high = max(arrays["liquid_saturations"])
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
