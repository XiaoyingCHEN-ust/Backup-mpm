#!/usr/bin/env python3
"""Measure pipe-wall phase leakage and local saturation oscillation."""

import argparse
import base64
import csv
import json
import math
import re
import struct
import xml.etree.ElementTree as ET
import zlib
from pathlib import Path

import numpy as np

PATTERN = re.compile(r"particle(\d+)\.vtp$")
VTK_DTYPES = {
    "Float32": np.dtype("f4"),
    "Float64": np.dtype("f8"),
    "Int32": np.dtype("i4"),
    "UInt32": np.dtype("u4"),
    "Int64": np.dtype("i8"),
    "UInt64": np.dtype("u8"),
}
WANTED_ARRAYS = {
    "liquid_saturations",
    "liquid_velocities",
    "gas_velocities",
    "velocities",
    "displacements",
}


def _decode_compressed_binary(text, byte_order, header_type):
    encoded = "".join(text.split())
    endian = "<" if byte_order == "LittleEndian" else ">"
    code = {"UInt32": "I", "UInt64": "Q"}.get(header_type)
    if code is None:
        raise ValueError(f"Unsupported VTK header type {header_type}")
    header_size = struct.calcsize(code)
    first_chars = 4 * math.ceil(header_size / 3)
    first = base64.b64decode(encoded[:first_chars])
    (blocks,) = struct.unpack(endian + code, first[:header_size])
    header_bytes = (3 + blocks) * header_size
    header_chars = 4 * math.ceil(header_bytes / 3)
    header = base64.b64decode(encoded[:header_chars])[:header_bytes]
    values = struct.unpack(endian + code * (3 + blocks), header)
    payload = base64.b64decode(encoded[header_chars:])
    output = bytearray()
    offset = 0
    for compressed_size in values[3:]:
        end = offset + compressed_size
        if end > len(payload):
            raise ValueError("Truncated compressed VTK payload")
        output.extend(zlib.decompress(payload[offset:end]))
        offset = end
    return bytes(output)


def _decode_array(element, byte_order, header_type, compressor, point_count):
    vtk_type = element.attrib["type"]
    if vtk_type not in VTK_DTYPES:
        raise ValueError(f"Unsupported VTK scalar type {vtk_type}")
    components = int(element.attrib.get("NumberOfComponents", "1"))
    dtype = VTK_DTYPES[vtk_type].newbyteorder(
        "<" if byte_order == "LittleEndian" else ">"
    )
    data_format = element.attrib.get("format", "ascii")
    text = element.text or ""
    if data_format == "ascii":
        values = np.fromstring(text, sep=" ", dtype=dtype)
    elif data_format == "binary" and compressor == "vtkZLibDataCompressor":
        values = np.frombuffer(
            _decode_compressed_binary(text, byte_order, header_type), dtype=dtype
        )
    else:
        raise ValueError(
            f"Unsupported VTK encoding format={data_format}, compressor={compressor}"
        )
    expected = point_count * components
    if values.size != expected:
        name = element.attrib.get("Name", "<points>")
        raise ValueError(f"VTK array {name} has {values.size} values; expected {expected}")
    if components > 1:
        values = values.reshape(point_count, components)
    if not np.all(np.isfinite(values)):
        raise ValueError(f"VTK array {element.attrib.get('Name')} has NaN/Inf")
    return values


def read_vtp(path):
    root = ET.parse(path).getroot()
    byte_order = root.attrib.get("byte_order", "LittleEndian")
    header_type = root.attrib.get("header_type", "UInt32")
    compressor = root.attrib.get("compressor", "")
    piece = root.find(".//Piece")
    if piece is None:
        raise ValueError(f"No PolyData Piece in {path}")
    count = int(piece.attrib.get("NumberOfPoints", "0"))
    if count <= 0:
        raise ValueError(f"No particle points in {path}")
    points_element = piece.find("./Points/DataArray")
    if points_element is None:
        raise ValueError(f"No Points array in {path}")
    points = _decode_array(
        points_element, byte_order, header_type, compressor, count
    )
    arrays = {}
    for element in piece.findall("./PointData/DataArray"):
        name = element.attrib.get("Name", "")
        if name in WANTED_ARRAYS:
            arrays[name] = _decode_array(
                element, byte_order, header_type, compressor, count
            )
    return points, arrays


def pipeline_state(path):
    piece = ET.parse(path).getroot().find(".//Piece")
    positions = np.fromstring(
        piece.find("./Points/DataArray").text, sep=" "
    ).reshape(-1, 3)[:, :2]
    velocities = np.fromstring(
        piece.find("./PointData/DataArray").text, sep=" "
    ).reshape(-1, 3)[:, :2]
    return positions.mean(axis=0), velocities.mean(axis=0)


def load_pipe_config(path):
    with path.open(encoding="utf-8") as stream:
        config = json.load(stream)
    pipe = config["analysis"]["rigid_pipeline"]
    if not pipe.get("phase_no_flux", False):
        raise RuntimeError(f"{path}: phase_no_flux is not enabled")
    return {
        "radius": float(pipe["outer_radius"]),
        "particle_radius": float(pipe["particle_radius"]),
        "band_factor": float(pipe["phase_no_flux_band_factor"]),
        "cell_size": float(config["mesh"]["cellsize_min"]),
        "uuid": config["analysis"]["uuid"],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "result",
        type=Path,
        nargs="?",
        default=Path("results/Wave2D_SN_094_1E11_pipeline_hydrostatic_noflux_v3"),
    )
    parser.add_argument("--dt", type=float, default=1.0e-4)
    parser.add_argument("--config", type=Path, default=Path("mpm-3p.json"))
    args = parser.parse_args()

    pipe_config = load_pipe_config(args.config)
    influence_radius = (
        pipe_config["radius"]
        + pipe_config["band_factor"] * pipe_config["particle_radius"]
    )
    rows = []
    reference_centre = None
    for particle_path in sorted(args.result.glob("particle*.vtp")):
        match = PATTERN.match(particle_path.name)
        if match is None:
            continue
        step = int(match.group(1))
        pipeline_path = args.result / f"pipeline{step:0{len(match.group(1))}d}.vtp"
        if not pipeline_path.is_file():
            raise RuntimeError(f"Missing matching pipeline output: {pipeline_path}")

        points, arrays = read_vtp(particle_path)
        required = {
            "liquid_saturations",
            "liquid_velocities",
            "gas_velocities",
            "velocities",
            "displacements",
        }
        missing = required - arrays.keys()
        if missing:
            raise RuntimeError(f"{particle_path.name} is missing {sorted(missing)}")

        centre, pipe_velocity = pipeline_state(pipeline_path)
        if reference_centre is None:
            reference_centre = centre.copy()
        radial = points[:, :2] - centre
        radius = np.linalg.norm(radial, axis=1)
        wall = radius <= influence_radius + 1.0e-12
        if not np.any(wall):
            raise RuntimeError(f"{particle_path.name} has no particles near the pipe")
        normal = radial[wall] / radius[wall, None]

        liquid_normal = np.sum(
            (arrays["liquid_velocities"][wall, :2] - pipe_velocity) * normal,
            axis=1,
        )
        gas_normal = np.sum(
            (arrays["gas_velocities"][wall, :2] - pipe_velocity) * normal,
            axis=1,
        )
        saturation = arrays["liquid_saturations"]
        far = (
            (points[:, 0] >= 0.95)
            & (points[:, 0] <= 1.15)
            & (points[:, 1] >= centre[1] - 0.09)
            & (points[:, 1] <= centre[1] + 0.09)
        )
        displacement = np.linalg.norm(arrays["displacements"][:, :2], axis=1)
        rows.append(
            {
                "step": step,
                "time_s": step * args.dt,
                "pipe_displacement_m": float(
                    np.linalg.norm(centre - reference_centre)
                ),
                "wall_particle_count": int(np.sum(wall)),
                "wall_particle_radius_min_m": float(np.min(radius[wall])),
                "wall_particle_radius_max_m": float(np.max(radius[wall])),
                "wall_mean_abs_liquid_normal_velocity_m_s": float(
                    np.mean(np.abs(liquid_normal))
                ),
                "wall_max_abs_liquid_normal_velocity_m_s": float(
                    np.max(np.abs(liquid_normal))
                ),
                "wall_mean_abs_gas_normal_velocity_m_s": float(
                    np.mean(np.abs(gas_normal))
                ),
                "wall_max_abs_gas_normal_velocity_m_s": float(
                    np.max(np.abs(gas_normal))
                ),
                "wall_saturation_mean": float(np.mean(saturation[wall])),
                "wall_saturation_std": float(np.std(saturation[wall])),
                "wall_saturation_min": float(np.min(saturation[wall])),
                "wall_saturation_max": float(np.max(saturation[wall])),
                "far_saturation_std": float(np.std(saturation[far])),
                "max_soil_displacement_m": float(np.max(displacement)),
                "max_soil_displacement_over_cell": float(
                    np.max(displacement) / pipe_config["cell_size"]
                ),
            }
        )

    if not rows:
        raise RuntimeError(f"No particle VTP files found in {args.result}")
    output = args.result / "pipeline_wall_diagnostics.csv"
    with output.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)

    print(f"wrote {output}")
    post_update = [row for row in rows if row["step"] > 0]
    if not post_update:
        raise RuntimeError("No post-update particle frames were found")
    print(
        "max wall |v_liquid,n-v_pipe,n| = "
        f"{max(r['wall_max_abs_liquid_normal_velocity_m_s'] for r in post_update):.6g} m/s"
    )
    print(
        "max wall |v_gas,n-v_pipe,n| = "
        f"{max(r['wall_max_abs_gas_normal_velocity_m_s'] for r in post_update):.6g} m/s"
    )
    print(
        "max wall saturation std = "
        f"{max(r['wall_saturation_std'] for r in rows):.6g}"
    )
    print(
        "max pipe displacement = "
        f"{max(r['pipe_displacement_m'] for r in rows):.6g} m"
    )
    print(
        "max soil displacement/cell = "
        f"{max(r['max_soil_displacement_over_cell'] for r in rows):.6g}"
    )


if __name__ == "__main__":
    main()
