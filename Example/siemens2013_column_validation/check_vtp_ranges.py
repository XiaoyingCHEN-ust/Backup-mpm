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
    # DataArray can contain nested InformationKey values; only element.text
    # is the encoded payload. Including descendant text silently appends the
    # range values to the base64 stream.
    payload = "".join((element.text or "").split())
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
    encoded = payload[header_chars:]
    encoded += "=" * (-len(encoded) % 4)
    compressed = base64.b64decode(encoded)

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


def point_coordinates(path: Path):
    root = ET.parse(path).getroot()
    element = root.find("./PolyData/Piece/Points/DataArray")
    if element is None:
        raise ValueError(f"Point coordinates are missing in {path}")
    vtk_type = element.get("type")
    if vtk_type not in VTK_STRUCT_TYPES:
        raise ValueError(f"Unsupported point-coordinate type: {vtk_type}")
    raw, endian = decode_binary_array(
        element,
        root.get("byte_order", "LittleEndian"),
        root.get("header_type", "UInt32"),
        root.get("compressor", ""),
    )
    code = VTK_STRUCT_TYPES[vtk_type]
    values = tuple(value[0] for value in struct.iter_unpack(endian + code, raw))
    if len(values) % 3:
        raise ValueError(f"Invalid point-coordinate count in {path}")
    return tuple(zip(values[0::3], values[1::3], values[2::3]))


def maximum_layer_spread(path: Path, values):
    coordinates = point_coordinates(path)
    if len(coordinates) != len(values):
        raise ValueError(f"Coordinate/value count mismatch in {path}")

    layers = {}
    for coordinates_i, value in zip(coordinates, values):
        layers.setdefault(round(coordinates_i[1], 8), []).append(value)
    if not layers or any(
        len(layer_values) != 2 for layer_values in layers.values()
    ):
        raise RuntimeError(
            "The reduced column checker expects exactly two particles per "
            f"horizontal layer in {path.name}"
        )
    return max(
        max(layer_values) - min(layer_values)
        for layer_values in layers.values()
    )


def maximum_unconnected_layer_saturation(path: Path, values, wet_threshold):
    """Return the wettest layer not connected to either column boundary."""
    coordinates = point_coordinates(path)
    if len(coordinates) != len(values):
        raise ValueError(f"Coordinate/value count mismatch in {path}")

    layers = {}
    for coordinates_i, value in zip(coordinates, values):
        layers.setdefault(round(coordinates_i[1], 8), []).append(value)
    ordered_layers = sorted(layers)
    if not ordered_layers or any(len(layers[y]) != 2 for y in ordered_layers):
        raise RuntimeError(
            "The reduced column checker expects exactly two particles per "
            f"horizontal layer in {path.name}"
        )

    layer_saturations = [sum(layers[y]) / len(layers[y]) for y in ordered_layers]
    lower_connected = 0
    while (
        lower_connected < len(layer_saturations)
        and layer_saturations[lower_connected] >= wet_threshold
    ):
        lower_connected += 1

    upper_connected = len(layer_saturations)
    while (
        upper_connected > lower_connected
        and layer_saturations[upper_connected - 1] >= wet_threshold
    ):
        upper_connected -= 1

    if lower_connected == upper_connected:
        return 0.0, None
    interior = layer_saturations[lower_connected:upper_connected]
    local_index = max(range(len(interior)), key=interior.__getitem__)
    layer_index = lower_connected + local_index
    return interior[local_index], ordered_layers[layer_index]


def maximum_downward_layer_saturation_increase(path: Path, values):
    """Return the largest saturation increase from one layer to the one below."""
    coordinates = point_coordinates(path)
    if len(coordinates) != len(values):
        raise ValueError(f"Coordinate/value count mismatch in {path}")

    layers = {}
    for coordinates_i, value in zip(coordinates, values):
        layers.setdefault(round(coordinates_i[1], 8), []).append(value)
    ordered_layers = sorted(layers)
    if len(ordered_layers) < 2 or any(
        len(layers[y]) != 2 for y in ordered_layers
    ):
        raise RuntimeError(
            "The reduced column checker expects exactly two particles per "
            f"horizontal layer in {path.name}"
        )

    layer_saturations = [
        sum(layers[y]) / len(layers[y]) for y in ordered_layers
    ]
    increase, lower_index = max(
        (layer_saturations[index] - layer_saturations[index + 1], index)
        for index in range(len(layer_saturations) - 1)
    )
    return max(increase, 0.0), (
        ordered_layers[lower_index],
        ordered_layers[lower_index + 1],
    )


def maximum_top_cell_abs_value(path: Path, values):
    """Return the maximum absolute scalar value in the four top particles."""
    coordinates = point_coordinates(path)
    if len(coordinates) != len(values):
        raise ValueError(f"Coordinate/value count mismatch in {path}")

    layers = {}
    for index, coordinates_i in enumerate(coordinates):
        layers.setdefault(round(coordinates_i[1], 8), []).append(index)
    ordered_layers = sorted(layers)
    if len(ordered_layers) < 2 or any(
        len(layers[y]) != 2 for y in ordered_layers
    ):
        raise RuntimeError(
            "The reduced column checker expects exactly two particles per "
            f"horizontal layer in {path.name}"
        )
    top_cell_ids = [
        index for y in ordered_layers[-2:] for index in layers[y]
    ]
    return max(abs(values[index]) for index in top_cell_ids)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("result_dir", type=Path)
    parser.add_argument("--expected-suction-pa", type=float, default=972.9920291627446)
    parser.add_argument(
        "--skip-initial-suction-check",
        action="store_true",
        help=(
            "skip the zero-time suction-patch check for a resumed directory "
            "whose first output occurs after the checkpoint"
        ),
    )
    parser.add_argument(
        "--require-force-pressure-fields",
        action="store_true",
        help=(
            "require and validate the force_liquid_pressures and "
            "force_gas_pressures diagnostic arrays"
        ),
    )
    parser.add_argument(
        "--require-pressure-gradient-fields",
        action="store_true",
        help=(
            "require and validate the liquid_pressure_gradients and "
            "gas_pressure_gradients diagnostic arrays"
        ),
    )
    parser.add_argument("--pressure-limit-pa", type=float, default=1.0e6)
    parser.add_argument("--velocity-limit-mps", type=float, default=10.0)
    parser.add_argument(
        "--layer-saturation-spread-limit",
        type=float,
        default=1.0e-4,
        help="maximum allowed saturation difference within a horizontal layer",
    )
    parser.add_argument(
        "--connected-wet-saturation-threshold",
        type=float,
        default=0.40,
        help=(
            "layer-mean saturation used to reject wet bands that are not "
            "connected to either column boundary"
        ),
    )
    parser.add_argument(
        "--downward-saturation-increase-limit",
        type=float,
        default=0.05,
        help=(
            "maximum allowed increase in layer-mean saturation from an upper "
            "layer to the adjacent lower layer"
        ),
    )
    parser.add_argument(
        "--surface-gas-pressure-limit-pa",
        type=float,
        default=1.0e-6,
        help="maximum absolute gas pressure allowed in the four top particles",
    )
    args = parser.parse_args()

    files = sorted(args.result_dir.glob("particle*.vtp"))
    if not files:
        raise FileNotFoundError(f"No particle VTP files in {args.result_dir}")

    marker = args.result_dir / "RANGE_CHECK_PASSED.txt"
    marker.unlink(missing_ok=True)
    names = (
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
    if args.require_force_pressure_fields:
        names += ("force_liquid_pressures", "force_gas_pressures")
    if args.require_pressure_gradient_fields:
        names += ("liquid_pressure_gradients", "gas_pressure_gradients")
    initial = point_arrays(files[0], names)
    if any(not math.isfinite(value) for value in initial["suction_pressures"]):
        raise RuntimeError("Initial suction array contains non-finite values")
    observed_suction = max(initial["suction_pressures"])
    if not args.skip_initial_suction_check and not math.isclose(
        observed_suction, args.expected_suction_pa, rel_tol=0.02, abs_tol=1.0
    ):
        raise RuntimeError(
            "Executable does not contain the validation patch: initial maximum "
            f"suction is {observed_suction:.6g} Pa, expected "
            f"{args.expected_suction_pa:.6g} Pa"
        )

    maximum_pressure = 0.0
    maximum_velocity = 0.0
    saturation_min = math.inf
    saturation_max = -math.inf
    gas_saturation_min = math.inf
    gas_saturation_max = -math.inf
    minimum_permeability = math.inf
    layer_saturation_spread = 0.0
    maximum_unconnected_saturation = 0.0
    maximum_unconnected_location = (None, None)
    maximum_downward_saturation_increase = 0.0
    maximum_downward_saturation_location = (None, None, None)
    maximum_surface_gas_pressure = 0.0
    maximum_surface_gas_pressure_file = None
    for path in files:
        arrays = point_arrays(path, names)
        for name, values in arrays.items():
            nonfinite = sum(not math.isfinite(value) for value in values)
            if nonfinite:
                raise RuntimeError(
                    f"Non-finite {name} values in {path.name}: {nonfinite}"
                )
        pressure_names = ["gas_pressures", "liquid_pressures"]
        if args.require_force_pressure_fields:
            pressure_names.extend(
                ("force_gas_pressures", "force_liquid_pressures")
            )
        for name in pressure_names:
            maximum_pressure = max(
                maximum_pressure, max(abs(value) for value in arrays[name])
            )
        for name in ("liquid_velocities", "gas_velocities"):
            maximum_velocity = max(
                maximum_velocity, max(abs(value) for value in arrays[name])
            )
        for name in ("liquid_permeabilities", "gas_permeabilities"):
            minimum_permeability = min(minimum_permeability, min(arrays[name]))
        low = min(arrays["liquid_saturations"])
        high = max(arrays["liquid_saturations"])
        saturation_min = min(saturation_min, low)
        saturation_max = max(saturation_max, high)
        gas_saturation_min = min(gas_saturation_min, min(arrays["gas_saturations"]))
        gas_saturation_max = max(gas_saturation_max, max(arrays["gas_saturations"]))
        maximum_saturation_sum_error = max(
            abs(liquid + gas - 1.0)
            for liquid, gas in zip(
                arrays["liquid_saturations"], arrays["gas_saturations"]
            )
        )
        if maximum_saturation_sum_error > 1.0e-10:
            raise RuntimeError(
                "Liquid and gas saturations do not sum to one in "
                f"{path.name}: error={maximum_saturation_sum_error:.6g}"
            )
        layer_saturation_spread = max(
            layer_saturation_spread,
            maximum_layer_spread(path, arrays["liquid_saturations"]),
        )
        unconnected_saturation, unconnected_y = (
            maximum_unconnected_layer_saturation(
                path,
                arrays["liquid_saturations"],
                args.connected_wet_saturation_threshold,
            )
        )
        if unconnected_saturation > maximum_unconnected_saturation:
            maximum_unconnected_saturation = unconnected_saturation
            maximum_unconnected_location = (path.name, unconnected_y)
        downward_increase, adjacent_y = (
            maximum_downward_layer_saturation_increase(
                path, arrays["liquid_saturations"]
            )
        )
        if downward_increase > maximum_downward_saturation_increase:
            maximum_downward_saturation_increase = downward_increase
            maximum_downward_saturation_location = (
                path.name,
                adjacent_y[0],
                adjacent_y[1],
            )
        surface_gas_pressure = maximum_top_cell_abs_value(
            path, arrays["gas_pressures"]
        )
        if surface_gas_pressure > maximum_surface_gas_pressure:
            maximum_surface_gas_pressure = surface_gas_pressure
            maximum_surface_gas_pressure_file = path.name

    if not (0.0 <= saturation_min <= saturation_max <= 1.0):
        raise RuntimeError(
            f"Nonphysical saturation range: [{saturation_min}, {saturation_max}]"
        )
    if not (0.0 <= gas_saturation_min <= gas_saturation_max <= 1.0):
        raise RuntimeError(
            "Nonphysical gas saturation range: "
            f"[{gas_saturation_min}, {gas_saturation_max}]"
        )
    if minimum_permeability <= 0.0:
        raise RuntimeError(
            f"Nonpositive phase permeability: {minimum_permeability:.6g} m^2"
        )
    if layer_saturation_spread > args.layer_saturation_spread_limit:
        raise RuntimeError(
            "One-dimensional layer symmetry check failed: maximum liquid-"
            f"saturation spread={layer_saturation_spread:.6g}, limit="
            f"{args.layer_saturation_spread_limit:.6g}"
        )
    if (
        maximum_unconnected_saturation
        >= args.connected_wet_saturation_threshold
    ):
        filename, layer_y = maximum_unconnected_location
        raise RuntimeError(
            "Disconnected wet-layer check failed: maximum layer-mean "
            f"saturation={maximum_unconnected_saturation:.6g} at y="
            f"{layer_y:.6g} m in {filename}; threshold="
            f"{args.connected_wet_saturation_threshold:.6g}"
        )
    if (
        maximum_downward_saturation_increase
        > args.downward_saturation_increase_limit
    ):
        filename, lower_y, upper_y = maximum_downward_saturation_location
        raise RuntimeError(
            "Nonmonotone wetting-profile check failed: layer-mean saturation "
            f"increased downward by {maximum_downward_saturation_increase:.6g} "
            f"between upper y={upper_y:.6g} m and lower y={lower_y:.6g} m "
            f"in {filename}; limit="
            f"{args.downward_saturation_increase_limit:.6g}"
        )
    if maximum_surface_gas_pressure > args.surface_gas_pressure_limit_pa:
        raise RuntimeError(
            "Atmospheric surface gas-pressure check failed: maximum |p_g|="
            f"{maximum_surface_gas_pressure:.6g} Pa in "
            f"{maximum_surface_gas_pressure_file}; limit="
            f"{args.surface_gas_pressure_limit_pa:.6g} Pa"
        )
    if not math.isfinite(maximum_pressure) or maximum_pressure >= args.pressure_limit_pa:
        raise RuntimeError(
            f"Pressure stability check failed: max |p|={maximum_pressure:.6g} Pa"
        )
    if maximum_velocity >= args.velocity_limit_mps:
        raise RuntimeError(
            "Velocity stability check failed: max component="
            f"{maximum_velocity:.6g} m/s"
        )

    marker.write_text(
        f"files={len(files)}\n"
        f"force_pressure_fields="
        f"{'required' if args.require_force_pressure_fields else 'not_required'}\n"
        f"pressure_gradient_fields="
        f"{'required' if args.require_pressure_gradient_fields else 'not_required'}\n"
        f"initial_suction_check="
        f"{'skipped' if args.skip_initial_suction_check else 'passed'}\n"
        f"initial_max_suction_pa={observed_suction:.12g}\n"
        f"max_abs_pressure_pa={maximum_pressure:.12g}\n"
        f"max_abs_velocity_component_mps={maximum_velocity:.12g}\n"
        f"min_phase_permeability_m2={minimum_permeability:.12g}\n"
        f"saturation_range={saturation_min:.12g},{saturation_max:.12g}\n"
        f"max_layer_saturation_spread={layer_saturation_spread:.12g}\n"
        f"max_unconnected_layer_saturation="
        f"{maximum_unconnected_saturation:.12g}\n"
        f"max_downward_layer_saturation_increase="
        f"{maximum_downward_saturation_increase:.12g}\n"
        f"max_abs_surface_gas_pressure_pa="
        f"{maximum_surface_gas_pressure:.12g}\n",
        encoding="utf-8",
    )
    print(f"PASS: {args.result_dir} (max |p|={maximum_pressure:.6g} Pa)")


if __name__ == "__main__":
    main()
