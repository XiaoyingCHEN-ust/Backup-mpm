#!/usr/bin/env python3
"""Reduce pipeline-study output to the manuscript comparison metrics.

Liquefaction is assessed from actual skeleton-stress loss and upward excess
seepage force.  Pore-pressure ratio (ru) and legacy solver flags are retained
only as diagnostics.  Pipeline displacement/contact topology and soil travel
relative to the background-cell width are reported in parallel as evidence for
the MPM choice; mesh-crossing travel is not itself a strain measure and does not
imply that an FEM treatment is impossible.
"""

from __future__ import annotations

import argparse
import base64
import csv
import hashlib
import json
import math
import re
import struct
import xml.etree.ElementTree as ET
import zlib
from copy import deepcopy
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

import numpy as np


PARTICLE_PATTERN = re.compile(r"particle(\d+)\.vtp$")
COMPLETION_SCHEMA = "pipeline-case-completion-v2"
COMPLETION_FILENAME = "pipeline_completion.json"
REGISTERED_STUDY_VALIDATION_PROFILE = "registered-study-v1"
STRESS_LOSS_RATIO_THRESHOLD = 0.05
MINIMUM_REFERENCE_VERTICAL_STRESS_PA = 100.0
THRESHOLD_RATIO_TOLERANCE = 1.0e-12
UPWARD_SEEPAGE_FORCE_THRESHOLD = 1.0
VTK_DTYPES = {
    "Float32": np.dtype("f4"),
    "Float64": np.dtype("f8"),
    "Int32": np.dtype("i4"),
    "UInt32": np.dtype("u4"),
    "Int64": np.dtype("i8"),
    "UInt64": np.dtype("u8"),
}
WANTED_ARRAYS = {
    "ids",
    "volumes",
    "PIC_pore_pressures",
    "PIC_pore_pressure_excess",
    "PIC_liquid_pressures",
    "PIC_gas_pressures",
    "PIC_ru",
    "initial_vertical_effective_stresses",
    "vertical_effective_stress_remaining_ratios",
    "liquefaction_potentials",
    "momentary_liquefied",
    "liquid_seepage_forces",
    "liquid_densities",
    "gamma_sub",
    "stresses",
    "displacements",
    "velocities",
    "porosities",
    "eps_p_q",
    "void_ratio",
    "pdstrain",
    "phi",
    "psi",
    "cohesion",
}


@dataclass(frozen=True)
class Probe:
    name: str
    index: int
    requested_x: float
    requested_y: float
    actual_x: float
    actual_y: float


def _decode_compressed_binary(text: str, byte_order: str, header_type: str) -> bytes:
    encoded = "".join(text.split())
    endian = "<" if byte_order == "LittleEndian" else ">"
    if header_type == "UInt32":
        code = "I"
    elif header_type == "UInt64":
        code = "Q"
    else:
        raise ValueError(f"Unsupported VTK header type {header_type}")
    header_size = struct.calcsize(code)
    first_chars = 4 * math.ceil(header_size / 3)
    first = base64.b64decode(encoded[:first_chars])
    (blocks,) = struct.unpack(endian + code, first[:header_size])
    header_bytes = (3 + blocks) * header_size
    header_chars = 4 * math.ceil(header_bytes / 3)
    header = base64.b64decode(encoded[:header_chars])[:header_bytes]
    values = struct.unpack(endian + code * (3 + blocks), header)
    output = bytearray()
    offset = 0
    payload = base64.b64decode(encoded[header_chars:])
    for compressed_size in values[3:]:
        end = offset + compressed_size
        if end > len(payload):
            raise ValueError("Truncated compressed VTK payload")
        output.extend(zlib.decompress(payload[offset:end]))
        offset = end
    return bytes(output)


def _decode_array(
    element: ET.Element,
    *,
    byte_order: str,
    header_type: str,
    compressor: str,
    point_count: int,
) -> np.ndarray:
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
        raise ValueError(
            f"VTK array {element.attrib.get('Name', '<points>')} has "
            f"{values.size} values; expected {expected}"
        )
    if components > 1:
        values = values.reshape(point_count, components)
    if not np.all(np.isfinite(values)):
        raise ValueError(f"VTK array {element.attrib.get('Name')} has NaN/Inf")
    return values


def read_vtp(path: Path) -> tuple[np.ndarray, dict[str, np.ndarray]]:
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
        points_element,
        byte_order=byte_order,
        header_type=header_type,
        compressor=compressor,
        point_count=count,
    )
    arrays: dict[str, np.ndarray] = {}
    for element in piece.findall("./PointData/DataArray"):
        name = element.attrib.get("Name", "")
        if name not in WANTED_ARRAYS:
            continue
        arrays[name] = _decode_array(
            element,
            byte_order=byte_order,
            header_type=header_type,
            compressor=compressor,
            point_count=count,
        )
    return points, arrays


def particle_step(path: Path) -> int:
    match = PARTICLE_PATTERN.match(path.name)
    if match is None:
        raise ValueError(f"Cannot parse output step from {path.name}")
    return int(match.group(1))


def file_sha256(path: Path) -> str:
    """Return the SHA-256 fingerprint used by completion sentinels."""

    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def expected_particle_steps(config: dict[str, Any]) -> list[int]:
    """Return the complete configured particle-output step grid."""

    nsteps = int(config["analysis"]["nsteps"])
    output_steps = int(config["post_processing"]["output_steps"])
    if nsteps <= 0 or output_steps <= 0:
        raise ValueError("Configured nsteps and output_steps must be positive")
    steps = list(range(output_steps, nsteps + 1, output_steps))
    if not steps or steps[-1] != nsteps:
        steps.append(nsteps)
    return steps


def _resolve_recorded_path(value: str | Path, case_root: Path) -> Path:
    path = Path(value)
    resolved = (path if path.is_absolute() else case_root / path).resolve()
    try:
        resolved.relative_to(case_root.resolve())
    except ValueError as error:
        raise ValueError(f"Recorded artifact escapes the case root: {path}") from error
    return resolved


def _verify_nested_artifact_audits(
    value: Any, case_root: Path, audited: list[dict[str, Any]]
) -> None:
    """Recursively verify every {path,size_bytes,sha256} sentinel record."""

    if isinstance(value, dict):
        if {"path", "size_bytes", "sha256"}.issubset(value):
            path = _resolve_recorded_path(value["path"], case_root)
            if not path.is_file():
                raise FileNotFoundError(f"Audited runtime artifact is missing: {path}")
            size = path.stat().st_size
            digest = file_sha256(path)
            if size <= 0 or int(value["size_bytes"]) != size:
                raise ValueError(f"Audited runtime artifact size is stale: {path}")
            if str(value["sha256"]) != digest:
                raise ValueError(f"Audited runtime artifact hash is stale: {path}")
            audited.append(
                {"path": str(path), "size_bytes": size, "sha256": digest}
            )
            return
        for child in value.values():
            _verify_nested_artifact_audits(child, case_root, audited)
    elif isinstance(value, list):
        for child in value:
            _verify_nested_artifact_audits(child, case_root, audited)


def validate_completed_result(
    config_path: Path,
    case_root: Path,
    config: dict[str, Any],
    result_directory: Path,
    vtp_files: Iterable[Path],
) -> dict[str, Any]:
    """Reject unaudited, stale, final-frame-only, or partial result sets.

    A valid manuscript input must have the exact configured particle-frame
    grid and a completion sentinel that fingerprints both the current config
    and the final VTP.  Merely finding a final-looking filename is not enough.
    """

    config_path = config_path.resolve()
    case_root = case_root.resolve()
    result_directory = result_directory.resolve()
    files = sorted((Path(path).resolve() for path in vtp_files), key=particle_step)
    actual_steps = [particle_step(path) for path in files]
    expected_steps = expected_particle_steps(config)
    if actual_steps != expected_steps:
        missing = sorted(set(expected_steps) - set(actual_steps))
        unexpected = sorted(set(actual_steps) - set(expected_steps))
        raise ValueError(
            "Particle VTP time grid is incomplete or unexpected: "
            f"expected {len(expected_steps)} frames through step {expected_steps[-1]}, "
            f"found {len(actual_steps)}; missing={missing[:8]}, "
            f"unexpected={unexpected[:8]}"
        )

    nsteps = int(config["analysis"]["nsteps"])
    final_vtp = result_directory / f"particle{nsteps:0{len(str(nsteps))}d}.vtp"
    if files[-1] != final_vtp.resolve() or not final_vtp.is_file():
        raise FileNotFoundError(f"Configured final particle VTP is missing: {final_vtp}")

    sentinel_path = result_directory / COMPLETION_FILENAME
    if not sentinel_path.is_file():
        raise FileNotFoundError(
            f"Audited completion sentinel is missing: {sentinel_path}"
        )
    sentinel = json.loads(sentinel_path.read_text(encoding="utf-8"))
    if sentinel.get("schema") != COMPLETION_SCHEMA:
        raise ValueError("Completion sentinel schema is invalid")
    record = sentinel.get("config")
    if not isinstance(record, dict):
        raise ValueError("Completion sentinel has no config audit")
    if _resolve_recorded_path(record["path"], case_root) != config_path:
        raise ValueError("Completion sentinel references another config")
    config_hash = file_sha256(config_path)
    if str(record.get("sha256")) != config_hash:
        raise ValueError("Completion sentinel config hash is stale")
    if str(record.get("uuid")) != str(config["analysis"]["uuid"]):
        raise ValueError("Completion sentinel UUID is stale")
    if int(record.get("nsteps", -1)) != nsteps:
        raise ValueError("Completion sentinel nsteps is stale")
    if record.get("validation_profile") != REGISTERED_STUDY_VALIDATION_PROFILE:
        raise ValueError(
            "Completion sentinel is not registered-study manuscript evidence: "
            f"validation_profile={record.get('validation_profile')!r}"
        )

    artifact = sentinel.get("artifacts", {}).get("final_vtp")
    if not isinstance(artifact, dict):
        raise ValueError("Completion sentinel has no final-VTP audit")
    if _resolve_recorded_path(artifact["path"], case_root) != final_vtp.resolve():
        raise ValueError("Completion sentinel final-VTP path is stale")
    size = final_vtp.stat().st_size
    final_hash = file_sha256(final_vtp)
    if size <= 0 or int(artifact.get("size_bytes", -1)) != size:
        raise ValueError("Completion sentinel final-VTP size is stale")
    if str(artifact.get("sha256")) != final_hash:
        raise ValueError("Completion sentinel final-VTP hash is stale")

    audited_dependencies: list[dict[str, Any]] = []
    _verify_nested_artifact_audits(
        sentinel.get("runtime_dependencies", {}),
        case_root,
        audited_dependencies,
    )

    dt = float(config["analysis"]["dt"])
    return {
        "schema": COMPLETION_SCHEMA,
        "sentinel_path": str(sentinel_path.resolve()),
        "sentinel_sha256": file_sha256(sentinel_path),
        "config_sha256": config_hash,
        "final_vtp_sha256": final_hash,
        "expected_particle_steps": expected_steps,
        "expected_particle_times_s": [step * dt for step in expected_steps],
        "runtime_dependencies": sentinel.get("runtime_dependencies", {}),
        "verified_runtime_artifacts": audited_dependencies,
    }


def read_initial_coordinates(path: Path) -> np.ndarray:
    with path.open(encoding="utf-8") as stream:
        first = stream.readline().strip()
        expected = int(first)
        coordinates = np.loadtxt(stream, dtype=np.float64)
    if coordinates.ndim == 1:
        coordinates = coordinates.reshape(1, -1)
    if len(coordinates) != expected:
        raise ValueError(
            f"{path} declares {expected} particles but contains {len(coordinates)}"
        )
    return coordinates[:, :2]


def resolve_particle_input_path(config: dict[str, Any], case_root: Path) -> Path:
    """Resolve the configured particle generator location against the case root."""

    location = Path(config["particles"][0]["generator"]["location"])
    return location if location.is_absolute() else case_root / location


def checkpoint_pressure_reference(
    completion_audit: dict[str, Any], case_root: Path, expected_count: int
) -> tuple[np.ndarray | None, dict[str, Any]]:
    """Load the audited equilibrium pressure reference in particle-ID order."""

    runtime = completion_audit.get("runtime_dependencies", {})
    resume = runtime.get("resume_equilibrium")
    if not isinstance(resume, dict) or not isinstance(resume.get("qa_vtp"), dict):
        return None, {"available": False, "reason": "case has no resume checkpoint"}
    record = resume["qa_vtp"]
    path = _resolve_recorded_path(record["path"], case_root)
    _, arrays = read_vtp(path)
    required = ("ids", "PIC_pore_pressures")
    missing = [name for name in required if name not in arrays]
    if missing:
        return None, {
            "available": False,
            "reason": f"audited checkpoint VTP lacks {missing}",
            "path": str(path),
            "sha256": record.get("sha256"),
        }
    ids = np.asarray(arrays["ids"], dtype=np.float64).reshape(-1)
    pressure = np.asarray(arrays["PIC_pore_pressures"], dtype=np.float64).reshape(-1)
    if ids.size != expected_count or pressure.size != expected_count:
        raise ValueError("Checkpoint pressure-reference particle count is invalid")
    integer_ids = np.rint(ids).astype(np.int64)
    order = np.argsort(integer_ids)
    if not np.array_equal(integer_ids[order], np.arange(expected_count)):
        raise ValueError("Checkpoint pressure-reference particle IDs are invalid")
    if not np.all(np.isfinite(pressure)):
        raise ValueError("Checkpoint pressure reference contains NaN/Inf")
    return pressure[order], {
        "available": True,
        "method": "audited_resume_checkpoint_PIC_pore_pressures_by_particle_id",
        "path": str(path),
        "sha256": record.get("sha256"),
    }


def nearest_probe(
    coordinates: np.ndarray, name: str, x: float, y: float
) -> Probe:
    distances = np.sum((coordinates - np.asarray([x, y])) ** 2, axis=1)
    index = int(np.argmin(distances))
    actual = coordinates[index]
    return Probe(name, index, x, y, float(actual[0]), float(actual[1]))


def build_probes(config: dict[str, Any], coordinates: np.ndarray) -> list[Probe]:
    pipeline = config["analysis"]["rigid_pipeline"]
    center_x, center_y = map(float, pipeline["center"])
    radius = float(pipeline["outer_radius"])
    soil_top = float(np.max(coordinates[:, 1]))
    offset = max(0.005, float(config["mesh"]["cellsize_min"]) / 4.0)
    requests = (
        ("surface", center_x, soil_top),
        ("crown", center_x, center_y + radius + offset),
        ("shoulder", center_x + radius + offset, center_y),
        ("invert", center_x, center_y - radius - offset),
        ("farfield", center_x + 3.0 * 2.0 * radius, center_y),
    )
    return [nearest_probe(coordinates, *request) for request in requests]


def read_pipeline_history(result_directory: Path) -> list[dict[str, float]]:
    paths = sorted(result_directory.glob("pipeline-history*.csv"))
    if not paths:
        raise FileNotFoundError(f"No pipeline history CSV in {result_directory}")
    by_step: dict[int, dict[str, float]] = {}
    for path in paths:
        with path.open(newline="", encoding="utf-8") as stream:
            for raw in csv.DictReader(stream):
                row = {key: float(value) for key, value in raw.items()}
                by_step[int(row["step"])] = row
    return [by_step[step] for step in sorted(by_step)]


def nearest_history_row(
    history: list[dict[str, float]], target_time: float
) -> dict[str, float]:
    return min(history, key=lambda row: abs(row["time"] - target_time))


def validate_pipeline_history_grid(
    history: list[dict[str, float]], config: dict[str, Any]
) -> dict[str, Any]:
    """Require every configured rigid-pipeline history sample through nsteps."""

    if not history:
        raise ValueError("Rigid-pipeline history is empty")
    analysis = config["analysis"]
    interval = int(analysis["rigid_pipeline"]["history_interval"])
    nsteps = int(analysis["nsteps"])
    dt = float(analysis["dt"])
    if interval <= 0:
        raise ValueError("Rigid-pipeline history interval must be positive")
    expected_steps = list(range(interval, nsteps + 1, interval))
    actual_steps = [int(row["step"]) for row in history]
    if actual_steps != expected_steps:
        missing = sorted(set(expected_steps) - set(actual_steps))
        unexpected = sorted(set(actual_steps) - set(expected_steps))
        raise ValueError(
            "Rigid-pipeline history grid is incomplete or unexpected: "
            f"expected {len(expected_steps)} rows, found {len(actual_steps)}; "
            f"missing={missing[:8]}, unexpected={unexpected[:8]}"
        )
    for row, step in zip(history, expected_steps):
        expected_time = step * dt
        if not math.isclose(
            float(row["time"]), expected_time, rel_tol=1.0e-10, abs_tol=1.0e-12
        ):
            raise ValueError(
                f"Rigid-pipeline history time at step {step} is "
                f"{row['time']}; expected {expected_time}"
            )
    return {
        "expected_pipeline_history_steps": expected_steps,
        "expected_pipeline_history_times_s": [step * dt for step in expected_steps],
    }


def pipeline_metrics(
    history: list[dict[str, float]],
    diameter: float,
    release: float,
    period: float,
    cover: float,
) -> tuple[dict[str, Any], list[dict[str, float]]]:
    after_release = [row for row in history if row["time"] >= release]
    if not after_release:
        raise ValueError("Pipeline history does not reach its release time")
    max_uplift = max(row["displacement_y"] for row in after_release)
    max_abs_vertical = max(abs(row["displacement_y"]) for row in after_release)
    max_abs_horizontal = max(abs(row["displacement_x"]) for row in after_release)
    max_abs_rotation = max(abs(row["angle"]) for row in after_release)
    no_contact = sum(row["contacts"] < 0.5 for row in after_release)
    breakout_rows = [row for row in after_release if row["displacement_y"] >= cover]
    cycle_rows: list[dict[str, float]] = []
    cycles = int(math.floor((history[-1]["time"] - release) / period))
    for cycle in range(cycles + 1):
        row = nearest_history_row(history, release + cycle * period)
        cycle_rows.append(
            {
                "cycle_after_release": float(cycle),
                "time_s": row["time"],
                "uplift_over_D": row["displacement_y"] / diameter,
                "horizontal_over_D": row["displacement_x"] / diameter,
                "rotation_rad": row["angle"],
                "contacts": row["contacts"],
            }
        )
    normalized = max_abs_vertical / diameter
    no_contact_fraction = no_contact / len(after_release)
    breakout_reached = bool(breakout_rows)
    topology_change = breakout_reached and no_contact_fraction >= 0.05
    return (
        {
            "final_uplift_m": after_release[-1]["displacement_y"],
            "final_uplift_over_D": after_release[-1]["displacement_y"] / diameter,
            "maximum_uplift_m": max_uplift,
            "maximum_abs_vertical_displacement_over_D": normalized,
            "maximum_abs_horizontal_displacement_over_D": max_abs_horizontal
            / diameter,
            "maximum_abs_rotation_rad": max_abs_rotation,
            "post_release_no_contact_fraction": no_contact_fraction,
            "maximum_contact_force_N_per_m": max(
                math.hypot(row["contact_force_x"], row["contact_force_y"])
                for row in after_release
            ),
            "large_deformation_threshold_over_D": 0.5,
            "breakout_uplift_over_D": cover / diameter,
            "breakout_reached": breakout_reached,
            "first_breakout_time_s": breakout_rows[0]["time"] if breakout_rows else None,
            "contact_topology_change_demonstrated": topology_change,
            "large_deformation_demonstrated": normalized >= 0.5 or topology_change,
            "model_scope_note": (
                "pipe loading includes gravity, full submerged buoyancy and soil contact; "
                "direct oscillatory hydrodynamic traction is not included, so post-breakout "
                "motion is an idealised support-loss response"
            ),
        },
        cycle_rows,
    )


def harmonic_fit(times: np.ndarray, values: np.ndarray, period: float) -> tuple[float, float, float]:
    omega = 2.0 * math.pi / period
    design = np.column_stack(
        (np.ones_like(times), np.cos(omega * times), np.sin(omega * times))
    )
    coefficients, _, _, _ = np.linalg.lstsq(design, values, rcond=None)
    mean, cosine, sine = map(float, coefficients)
    amplitude = math.hypot(cosine, sine)
    delay_phase = math.atan2(sine, cosine)
    return mean, amplitude, delay_phase


def wrap_degrees(angle_rad: float) -> float:
    return math.degrees((angle_rad + math.pi) % (2.0 * math.pi) - math.pi)


def stress_invariants(
    stress: np.ndarray, compression_sign: float = -1.0
) -> tuple[float, float]:
    xx, yy, zz, xy, yz, xz = map(float, stress[:6])
    if compression_sign not in (-1.0, 1.0):
        raise ValueError("Compression sign must be -1 or +1")
    mean_effective = compression_sign * (xx + yy + zz) / 3.0
    j2 = (
        ((xx - yy) ** 2 + (yy - zz) ** 2 + (zz - xx) ** 2) / 6.0
        + xy**2
        + yz**2
        + xz**2
    )
    return mean_effective, math.sqrt(max(0.0, 3.0 * j2))


def vertical_stress_loss_state(
    vertical_stress: np.ndarray,
    reference_vertical_stress: np.ndarray,
    *,
    minimum_reference_pa: float = MINIMUM_REFERENCE_VERTICAL_STRESS_PA,
    threshold: float = STRESS_LOSS_RATIO_THRESHOLD,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return R_sigma, eligible points and the 95% skeleton-stress-loss mask.

    R_sigma = sigma'_yy(t) / sigma'_yy(0), so both positive- and
    negative-compression conventions are supported without an absolute-value
    transformation.  A sign reversal is consequently treated as full loss.
    """

    vertical_stress = np.asarray(vertical_stress, dtype=np.float64)
    reference_vertical_stress = np.asarray(
        reference_vertical_stress, dtype=np.float64
    )
    if vertical_stress.ndim != 1 or vertical_stress.shape != reference_vertical_stress.shape:
        raise ValueError("Current and reference vertical stresses must be equal 1D arrays")
    if not np.all(np.isfinite(vertical_stress)) or not np.all(
        np.isfinite(reference_vertical_stress)
    ):
        raise ValueError("Vertical skeleton stresses contain NaN/Inf")
    if not math.isfinite(minimum_reference_pa) or minimum_reference_pa <= 0.0:
        raise ValueError("Minimum reference vertical stress must be positive")
    if not math.isfinite(threshold) or not 0.0 < threshold < 1.0:
        raise ValueError("Stress-loss ratio threshold must lie between zero and one")

    eligible = np.abs(reference_vertical_stress) >= minimum_reference_pa
    ratios = np.full(vertical_stress.shape, np.nan, dtype=np.float64)
    ratios[eligible] = (
        vertical_stress[eligible] / reference_vertical_stress[eligible]
    )
    lost = eligible & (ratios <= threshold + THRESHOLD_RATIO_TOLERANCE)
    return ratios, eligible, lost


def upward_excess_seepage_force_index(
    liquid_seepage_forces: np.ndarray,
    liquid_densities: np.ndarray,
    gamma_sub: np.ndarray,
    gravity: np.ndarray,
) -> np.ndarray:
    """Compute wave-induced upward excess seepage force / submerged weight.

    The saved liquid seepage force is ``-grad(p_l)`` and therefore includes the
    hydrostatic gradient.  Adding ``rho_l * g`` removes that static component;
    the residual is projected along the direction opposite gravity before it
    is normalised by ``gamma_sub``.
    """

    forces = np.asarray(liquid_seepage_forces, dtype=np.float64)
    densities = np.asarray(liquid_densities, dtype=np.float64)
    submerged_weights = np.asarray(gamma_sub, dtype=np.float64)
    gravity = np.asarray(gravity, dtype=np.float64)
    if forces.ndim != 2 or forces.shape[1] < 2:
        raise ValueError("Liquid seepage forces must be an N-by-dimension array")
    dimension = gravity.size
    if dimension < 2 or forces.shape[1] < dimension:
        raise ValueError("Gravity dimension is incompatible with seepage forces")
    if densities.ndim != 1 or submerged_weights.ndim != 1:
        raise ValueError("Liquid densities and gamma_sub must be 1D arrays")
    if densities.size != forces.shape[0] or submerged_weights.size != forces.shape[0]:
        raise ValueError("Seepage-force, density and gamma_sub counts must match")
    if not np.all(np.isfinite(forces[:, :dimension])):
        raise ValueError("Liquid seepage forces contain NaN/Inf")
    if not np.all(np.isfinite(densities)) or np.any(densities <= 0.0):
        raise ValueError("Liquid densities must be finite and positive")
    if not np.all(np.isfinite(submerged_weights)) or np.any(submerged_weights <= 0.0):
        raise ValueError("gamma_sub must be finite and positive")
    if not np.all(np.isfinite(gravity)):
        raise ValueError("Gravity contains NaN/Inf")
    gravity_magnitude = float(np.linalg.norm(gravity))
    if gravity_magnitude <= 0.0:
        raise ValueError("Gravity magnitude must be positive")

    excess_force = forces[:, :dimension] + densities[:, np.newaxis] * gravity
    upward = -gravity / gravity_magnitude
    signed_index = (excess_force @ upward) / submerged_weights
    return np.maximum(0.0, signed_index)


def pipeline_support_roi(
    config: dict[str, Any], coordinates: np.ndarray
) -> np.ndarray | None:
    """Return the pre-registered support-zone mask, when pipe data exist."""

    pipeline = config.get("analysis", {}).get("rigid_pipeline")
    if not pipeline:
        return None
    liquid = config.get("materials", [{}, {}])[1]
    if "sea_level" not in liquid or "depth_left" not in liquid:
        return None
    center_x = float(pipeline["center"][0])
    diameter = 2.0 * float(pipeline["outer_radius"])
    bed_y = float(liquid["sea_level"]) - float(liquid["depth_left"])
    coordinates = np.asarray(coordinates, dtype=np.float64)
    return (
        (np.abs(coordinates[:, 0] - center_x) <= 1.5 * diameter)
        & (coordinates[:, 1] >= bed_y - 2.0 * diameter)
        & (coordinates[:, 1] <= bed_y)
    )


def nominal_reference_particle_area(coordinates: np.ndarray) -> float:
    """Return the represented area of one point on the initial regular grid.

    The pipeline study uses a uniform Cartesian material-point lattice.  Missing
    points in the pipe cavity do not change the underlying grid spacing.  This
    is deliberately a *reference* area: without particle volumes in the saved
    VTP it must not be presented as an exact current-configuration area after
    large deformation.
    """

    coordinates = np.asarray(coordinates, dtype=np.float64)
    if coordinates.ndim != 2 or coordinates.shape[1] < 2:
        raise ValueError("Particle coordinates must be an N-by-2 array")

    spacings: list[float] = []
    for axis in range(2):
        unique = np.unique(np.round(coordinates[:, axis], decimals=12))
        differences = np.diff(unique)
        differences = differences[differences > 1.0e-12]
        if differences.size == 0:
            raise ValueError("Cannot infer a two-dimensional particle spacing")
        spacing = float(np.min(differences))
        multiples = differences / spacing
        if not np.allclose(multiples, np.rint(multiples), rtol=1.0e-7, atol=1.0e-7):
            raise ValueError(
                "Initial particles do not lie on a regular Cartesian lattice; "
                "a count-based reference area would be misleading"
            )
        spacings.append(spacing)
    area = spacings[0] * spacings[1]
    if not math.isfinite(area) or area <= 0.0:
        raise ValueError("Inferred particle reference area is invalid")
    return area


def threshold_spatiotemporal_metrics(
    times: np.ndarray,
    masks: np.ndarray,
    particle_area: float,
    current_volumes: np.ndarray | None = None,
) -> dict[str, Any]:
    """Summarise a saved-frame threshold mask without inventing crossings.

    Integrated quantities use trapezoidal interpolation between saved frames.
    Resolved durations are conservative lower bounds: an interval is counted
    only when the threshold is met at both of its saved endpoints.
    """

    times = np.asarray(times, dtype=np.float64)
    masks = np.asarray(masks, dtype=bool)
    if times.ndim != 1 or masks.ndim != 2 or masks.shape[0] != times.size:
        raise ValueError("Threshold masks must have shape (number of times, particles)")
    if times.size == 0 or masks.shape[1] == 0:
        raise ValueError("Threshold history is empty")
    if not np.all(np.isfinite(times)):
        raise ValueError("Threshold-history times contain NaN/Inf")
    intervals = np.diff(times)
    if np.any(intervals <= 0.0):
        raise ValueError("Threshold-history times must be strictly increasing")
    if not math.isfinite(particle_area) or particle_area <= 0.0:
        raise ValueError("Particle reference area must be positive")

    counts = np.count_nonzero(masks, axis=1)
    fractions = counts.astype(np.float64) / masks.shape[1]
    areas = counts.astype(np.float64) * particle_area
    active = counts > 0

    if times.size > 1:
        integrated_fraction = float(
            np.sum(0.5 * (fractions[:-1] + fractions[1:]) * intervals)
        )
        integrated_area = float(
            np.sum(0.5 * (areas[:-1] + areas[1:]) * intervals)
        )
        resolved_intervals = active[:-1] & active[1:]
        resolved_duration = float(np.sum(intervals[resolved_intervals]))

        longest_resolved_duration = 0.0
        current_duration = 0.0
        for duration, resolved in zip(intervals, resolved_intervals):
            if resolved:
                current_duration += float(duration)
                longest_resolved_duration = max(
                    longest_resolved_duration, current_duration
                )
            else:
                current_duration = 0.0

    else:
        integrated_fraction = 0.0
        integrated_area = 0.0
        resolved_duration = 0.0
        longest_resolved_duration = 0.0

    active_indices = np.flatnonzero(active)
    result: dict[str, Any] = {
        "max_fraction": float(np.max(fractions)),
        "max_nominal_reference_area_m2": float(np.max(areas)),
        "time_integrated_fraction_s": integrated_fraction,
        "time_integrated_nominal_reference_area_m2_s": integrated_area,
        "resolved_duration_anywhere_lower_bound_s": resolved_duration,
        "longest_resolved_duration_anywhere_lower_bound_s": (
            longest_resolved_duration
        ),
        "active_saved_frames": int(np.count_nonzero(active)),
        "first_active_saved_time_s": (
            float(times[active_indices[0]]) if active_indices.size else None
        ),
        "last_active_saved_time_s": (
            float(times[active_indices[-1]]) if active_indices.size else None
        ),
    }
    if current_volumes is not None:
        current_volumes = np.asarray(current_volumes, dtype=np.float64)
        if current_volumes.shape != masks.shape:
            raise ValueError(
                "Current particle volumes must have the same frame/particle "
                "shape as the threshold masks"
            )
        if not np.all(np.isfinite(current_volumes)):
            raise ValueError("Current particle volumes contain NaN/Inf")
        if np.any(current_volumes <= 0.0):
            raise ValueError("Current particle volumes must all be positive")
        total_current_areas = np.sum(current_volumes, axis=1)
        current_areas = np.sum(np.where(masks, current_volumes, 0.0), axis=1)
        current_area_fractions = current_areas / total_current_areas
        if times.size > 1:
            integrated_current_area = float(
                np.sum(
                    0.5
                    * (current_areas[:-1] + current_areas[1:])
                    * intervals
                )
            )
            integrated_current_fraction = float(
                np.sum(
                    0.5
                    * (
                        current_area_fractions[:-1]
                        + current_area_fractions[1:]
                    )
                    * intervals
                )
            )
        else:
            integrated_current_area = 0.0
            integrated_current_fraction = 0.0
        result.update(
            {
                "area_basis": "current_particle_volumes",
                "max_area_m2": float(np.max(current_areas)),
                "max_area_fraction": float(np.max(current_area_fractions)),
                "time_integrated_area_m2_s": integrated_current_area,
                "time_integrated_area_fraction_s": integrated_current_fraction,
                "max_current_area_m2": float(np.max(current_areas)),
                "max_current_area_fraction": float(
                    np.max(current_area_fractions)
                ),
                "time_integrated_current_area_m2_s": integrated_current_area,
                "time_integrated_current_area_fraction_s": (
                    integrated_current_fraction
                ),
            }
        )
    else:
        result.update(
            {
                "area_basis": "nominal_reference_particle_area_fallback",
                "max_area_m2": float(np.max(areas)),
                "max_area_fraction": float(np.max(fractions)),
                "time_integrated_area_m2_s": integrated_area,
                "time_integrated_area_fraction_s": integrated_fraction,
            }
        )
    return result


def _linear_integral_between(
    times: np.ndarray, values: np.ndarray, start: float, end: float
) -> float:
    """Integrate a linearly interpolated saved-frame series on [start, end]."""

    lower = max(float(start), float(times[0]))
    upper = min(float(end), float(times[-1]))
    if upper <= lower:
        return 0.0
    interior = times[(times > lower) & (times < upper)]
    sample_times = np.concatenate(([lower], interior, [upper]))
    sample_values = np.interp(sample_times, times, values)
    intervals = np.diff(sample_times)
    return float(
        np.sum(0.5 * (sample_values[:-1] + sample_values[1:]) * intervals)
    )


def threshold_release_segments(
    times: np.ndarray,
    masks: np.ndarray,
    particle_area: float,
    release_time: float,
    current_volumes: np.ndarray | None = None,
) -> dict[str, Any]:
    """Split a threshold area history at the pre-registered pipe release."""

    times = np.asarray(times, dtype=np.float64)
    masks = np.asarray(masks, dtype=bool)
    if masks.ndim != 2 or masks.shape[0] != times.size:
        raise ValueError("Threshold masks and times have inconsistent shapes")
    if not math.isfinite(release_time):
        raise ValueError("Pipeline release time must be finite")
    if current_volumes is None:
        areas = np.count_nonzero(masks, axis=1).astype(float) * particle_area
        fractions = np.mean(masks, axis=1)
        basis = "nominal_reference_particle_area_fallback"
    else:
        volumes = np.asarray(current_volumes, dtype=np.float64)
        if volumes.shape != masks.shape:
            raise ValueError("Current volumes and threshold masks must match")
        areas = np.sum(np.where(masks, volumes, 0.0), axis=1)
        totals = np.sum(volumes, axis=1)
        fractions = areas / totals
        basis = "current_particle_volumes"

    output: dict[str, Any] = {"release_segment_area_basis": basis}
    spans = {
        "pre_release": (float(times[0]), min(release_time, float(times[-1]))),
        "post_release": (max(release_time, float(times[0])), float(times[-1])),
    }
    for name, (start, end) in spans.items():
        saved = (times >= start - THRESHOLD_RATIO_TOLERANCE) & (
            times <= end + THRESHOLD_RATIO_TOLERANCE
        )
        active_saved = saved & (np.count_nonzero(masks, axis=1) > 0)
        output[f"{name}_saved_frames"] = int(np.count_nonzero(saved))
        output[f"{name}_active_saved_frames"] = int(np.count_nonzero(active_saved))
        output[f"{name}_first_active_saved_time_s"] = (
            float(times[np.flatnonzero(active_saved)[0]])
            if np.any(active_saved)
            else None
        )
        output[f"{name}_max_area_m2"] = (
            float(np.max(areas[saved])) if np.any(saved) else 0.0
        )
        output[f"{name}_time_integrated_area_m2_s"] = _linear_integral_between(
            times, areas, start, end
        )
        output[f"{name}_time_integrated_area_fraction_s"] = (
            _linear_integral_between(times, fractions, start, end)
        )
    return output


def trigger_to_stress_loss_timing(
    times: np.ndarray,
    seepage_masks: np.ndarray,
    stress_loss_masks: np.ndarray,
    wave_period: float,
) -> dict[str, Any]:
    """Quantify particle-level IF ordering relative to first stress loss.

    The preceding-cycle test is strict in time: a simultaneous saved-frame
    exceedance is captured by the joint criterion but is not retroactively
    counted as a preceding trigger.
    """

    times = np.asarray(times, dtype=np.float64)
    seepage_masks = np.asarray(seepage_masks, dtype=bool)
    stress_loss_masks = np.asarray(stress_loss_masks, dtype=bool)
    if (
        times.ndim != 1
        or seepage_masks.shape != stress_loss_masks.shape
        or seepage_masks.ndim != 2
        or seepage_masks.shape[0] != times.size
    ):
        raise ValueError("IF/stress-loss histories have inconsistent shapes")
    if not math.isfinite(wave_period) or wave_period <= 0.0:
        raise ValueError("Wave period must be finite and positive")

    particle_count = seepage_masks.shape[1]
    first_if = np.full(particle_count, np.nan)
    first_stress = np.full(particle_count, np.nan)
    for index in range(particle_count):
        if_indices = np.flatnonzero(seepage_masks[:, index])
        stress_indices = np.flatnonzero(stress_loss_masks[:, index])
        if if_indices.size:
            first_if[index] = times[if_indices[0]]
        if stress_indices.size:
            first_stress[index] = times[stress_indices[0]]

    lost = np.isfinite(first_stress)
    lost_count = int(np.count_nonzero(lost))
    if_before = lost & np.isfinite(first_if) & (
        first_if < first_stress - THRESHOLD_RATIO_TOLERANCE
    )
    simultaneous = lost & np.isfinite(first_if) & np.isclose(
        first_if,
        first_stress,
        rtol=0.0,
        atol=THRESHOLD_RATIO_TOLERANCE,
    )
    stress_before = lost & np.isfinite(first_if) & (
        first_stress < first_if - THRESHOLD_RATIO_TOLERANCE
    )
    no_if = lost & ~np.isfinite(first_if)
    observable = lost & (
        first_stress - times[0] >= wave_period - THRESHOLD_RATIO_TOLERANCE
    )
    preceding = np.zeros(particle_count, dtype=bool)
    for index in np.flatnonzero(lost):
        start = first_stress[index] - wave_period
        preceding[index] = bool(
            np.any(
                seepage_masks[:, index]
                & (times >= start - THRESHOLD_RATIO_TOLERANCE)
                & (times < first_stress[index] - THRESHOLD_RATIO_TOLERANCE)
            )
        )

    def fraction(mask: np.ndarray, denominator: int) -> float | None:
        return float(np.count_nonzero(mask) / denominator) if denominator else None

    observable_count = int(np.count_nonzero(observable))
    return {
        "first_IF_active_saved_time_s": (
            float(np.nanmin(first_if)) if np.any(np.isfinite(first_if)) else None
        ),
        "first_stress_loss_active_saved_time_s": (
            float(np.nanmin(first_stress)) if lost_count else None
        ),
        "first_stress_loss_minus_first_IF_active_saved_time_s": (
            float(np.nanmin(first_stress) - np.nanmin(first_if))
            if lost_count and np.any(np.isfinite(first_if))
            else None
        ),
        "stress_loss_particles": lost_count,
        "stress_loss_particles_IF_first": int(np.count_nonzero(if_before)),
        "stress_loss_particles_simultaneous_first": int(
            np.count_nonzero(simultaneous)
        ),
        "stress_loss_particles_stress_first": int(np.count_nonzero(stress_before)),
        "stress_loss_particles_without_IF": int(np.count_nonzero(no_if)),
        "fraction_stress_loss_particles_IF_first": fraction(if_before, lost_count),
        "stress_loss_particles_with_IF_in_preceding_cycle": int(
            np.count_nonzero(preceding)
        ),
        "fraction_stress_loss_particles_with_IF_in_preceding_cycle": fraction(
            preceding, lost_count
        ),
        "preceding_cycle_fully_observable_stress_loss_particles": observable_count,
        "fraction_fully_observable_stress_loss_particles_with_IF_in_preceding_cycle": (
            fraction(preceding & observable, observable_count)
        ),
        "preceding_cycle_window_s": wave_period,
        "preceding_cycle_rule": (
            "same particle has IF>=1 at a saved time in "
            "[first R_sigma<=0.05 - T, first R_sigma<=0.05); "
            "simultaneous onset is reported separately"
        ),
    }


def add_frame_threshold_metrics(
    row: dict[str, Any],
    name: str,
    mask: np.ndarray,
    particle_area: float,
    current_volumes: np.ndarray | None,
    domain_mask: np.ndarray | None = None,
) -> None:
    """Add count/reference-area and preferred current-area frame metrics."""

    mask = np.asarray(mask, dtype=bool)
    domain = (
        np.ones(mask.shape, dtype=bool)
        if domain_mask is None
        else np.asarray(domain_mask, dtype=bool)
    )
    if mask.ndim != 1 or domain.shape != mask.shape:
        raise ValueError("Threshold and reporting-domain masks must match")
    domain_count = int(np.count_nonzero(domain))
    if domain_count == 0:
        raise ValueError(f"Reporting domain for {name} contains no particles")
    affected = mask & domain
    affected_count = int(np.count_nonzero(affected))
    count_fraction = affected_count / domain_count
    nominal_area = affected_count * particle_area
    row[f"fraction_{name}"] = count_fraction
    row[f"nominal_reference_area_{name}_m2"] = nominal_area
    if current_volumes is not None:
        domain_area = float(np.sum(current_volumes[domain]))
        current_area = float(np.sum(current_volumes[affected]))
        current_fraction = current_area / domain_area
        row[f"current_area_{name}_m2"] = current_area
        row[f"current_area_fraction_{name}"] = current_fraction
        row[f"area_{name}_m2"] = current_area
        row[f"area_fraction_{name}"] = current_fraction
    else:
        row[f"area_{name}_m2"] = nominal_area
        row[f"area_fraction_{name}"] = count_fraction


def displacement_mesh_crossing_metrics(
    displacements: np.ndarray,
    particle_count: int,
    cell_size: float,
    domain_mask: np.ndarray | None = None,
) -> dict[str, Any]:
    """Measure material-point travel relative to one background-cell width."""

    displacements = np.asarray(displacements, dtype=np.float64)
    if (
        displacements.ndim != 2
        or displacements.shape[0] != particle_count
        or displacements.shape[1] not in (2, 3)
    ):
        raise ValueError(
            "Particle displacements must have shape (particle_count, 2 or 3)"
        )
    if not np.all(np.isfinite(displacements)):
        raise ValueError("Particle displacements contain NaN/Inf")
    if not math.isfinite(cell_size) or cell_size <= 0.0:
        raise ValueError("Background cell size must be finite and positive")
    domain = (
        np.ones(particle_count, dtype=bool)
        if domain_mask is None
        else np.asarray(domain_mask, dtype=bool)
    )
    if domain.ndim != 1 or domain.size != particle_count:
        raise ValueError("Displacement reporting-domain mask has an invalid shape")
    if not np.any(domain):
        raise ValueError("Displacement reporting domain contains no particles")

    magnitudes = np.linalg.norm(displacements[:, :2], axis=1)
    selected = magnitudes[domain]
    maximum = float(np.max(selected))
    return {
        "max_abs_displacement_m": maximum,
        "max_abs_displacement_over_h": maximum / cell_size,
        "fraction_abs_displacement_ge_h": float(np.mean(selected >= cell_size)),
        "crossed_one_background_cell": bool(maximum >= cell_size),
    }


def particle_histories(
    files: Iterable[Path],
    config: dict[str, Any],
    coordinates: np.ndarray,
    probes: list[Probe],
    checkpoint_pressure: np.ndarray | None = None,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    dt = float(config["analysis"]["dt"])
    rows: list[dict[str, Any]] = []
    expected_count = len(coordinates)
    if checkpoint_pressure is not None:
        checkpoint_pressure = np.asarray(
            checkpoint_pressure, dtype=np.float64
        ).reshape(-1)
        if checkpoint_pressure.size != expected_count or not np.all(
            np.isfinite(checkpoint_pressure)
        ):
            raise ValueError(
                "Checkpoint pore-pressure reference does not match particle input"
            )
    particle_area = nominal_reference_particle_area(coordinates)
    stress_threshold_name = "stress_loss_Rsigma_le_0p05"
    seepage_threshold_name = "upward_seepage_IF_ge_1"
    joint_threshold_name = "joint_IF_ge_1_and_Rsigma_le_0p05"
    threshold_masks: dict[str, list[np.ndarray | None]] = {
        stress_threshold_name: [],
        seepage_threshold_name: [],
        joint_threshold_name: [],
    }
    volume_frames: list[np.ndarray | None] = []
    support_roi = pipeline_support_roi(config, coordinates)
    if support_roi is not None and not np.any(support_roi):
        raise ValueError("The pre-registered pipeline support ROI is empty")
    gravity = np.asarray(
        config.get("external_loading_conditions", {}).get(
            "gravity", [0.0, -9.81]
        ),
        dtype=np.float64,
    )
    reference_vertical_stress: np.ndarray | None = None
    reference_eligible: np.ndarray | None = None
    reference_time: float | None = None
    reference_basis: str | None = None
    solver_initial_vertical_present: list[bool] = []
    solver_ratio_max_abs_differences: list[float] = []
    solver_ratio_present: list[bool] = []
    particle_ids_present: list[bool] = []
    displacement_frames_present: list[bool] = []
    cell_size_value = config.get("mesh", {}).get("cellsize_min")
    cell_size = float(cell_size_value) if cell_size_value is not None else None
    release_time = float(
        config.get("analysis", {}).get("rigid_pipeline", {}).get(
            "release_time", 0.0
        )
    )
    for path in sorted(files, key=particle_step):
        step = particle_step(path)
        points, arrays = read_vtp(path)
        if len(points) != expected_count:
            raise ValueError(f"Particle count changed in {path}")
        ids_available = "ids" in arrays
        particle_ids_present.append(ids_available)
        if ids_available:
            raw_ids = np.asarray(arrays["ids"], dtype=np.float64).reshape(-1)
            if raw_ids.size != expected_count or not np.all(np.isfinite(raw_ids)):
                raise ValueError(f"Particle ids in {path} are invalid")
            integer_ids = np.rint(raw_ids).astype(np.int64)
            if not np.array_equal(raw_ids, integer_ids.astype(np.float64)):
                raise ValueError(f"Particle ids in {path} are not integral")
            if np.unique(integer_ids).size != expected_count:
                raise ValueError(f"Particle ids in {path} are not unique")
            order = np.argsort(integer_ids)
            if not np.array_equal(integer_ids[order], np.arange(expected_count)):
                raise ValueError(
                    f"Particle ids in {path} do not match initial ids 0..N-1"
                )
            points = points[order]
            arrays = {name: values[order] for name, values in arrays.items()}
        row: dict[str, Any] = {"step": float(step), "time_s": step * dt}
        current_volumes: np.ndarray | None = None
        if "volumes" in arrays:
            current_volumes = np.asarray(arrays["volumes"], dtype=np.float64)
            if current_volumes.ndim != 1 or current_volumes.size != expected_count:
                raise ValueError(
                    f"Particle volume count in {path} does not match its "
                    f"{expected_count} points"
                )
            if not np.all(np.isfinite(current_volumes)):
                raise ValueError(f"Particle volumes in {path} contain NaN/Inf")
            if np.any(current_volumes <= 0.0):
                raise ValueError(f"Particle volumes in {path} must all be positive")
            row["current_total_area_m2"] = float(np.sum(current_volumes))
            volume_frames.append(current_volumes)
        else:
            volume_frames.append(None)

        displacements_present = "displacements" in arrays
        displacement_frames_present.append(displacements_present)
        if displacements_present:
            if cell_size is None:
                raise ValueError(
                    f"Particle displacements are present in {path}, but "
                    "config mesh.cellsize_min is missing"
                )
            soil_displacement = displacement_mesh_crossing_metrics(
                arrays["displacements"], expected_count, cell_size
            )
            if row["time_s"] >= release_time - THRESHOLD_RATIO_TOLERANCE:
                row["soil_max_abs_displacement_m"] = soil_displacement[
                    "max_abs_displacement_m"
                ]
                row["soil_max_abs_displacement_over_h"] = soil_displacement[
                    "max_abs_displacement_over_h"
                ]
                row["soil_fraction_abs_displacement_ge_h"] = soil_displacement[
                    "fraction_abs_displacement_ge_h"
                ]
                if support_roi is not None:
                    support_displacement = displacement_mesh_crossing_metrics(
                        arrays["displacements"],
                        expected_count,
                        cell_size,
                        support_roi,
                    )
                    row["support_ROI_soil_max_abs_displacement_m"] = (
                        support_displacement["max_abs_displacement_m"]
                    )
                    row["support_ROI_soil_max_abs_displacement_over_h"] = (
                        support_displacement["max_abs_displacement_over_h"]
                    )
                    row["support_ROI_soil_fraction_abs_displacement_ge_h"] = (
                        support_displacement[
                            "fraction_abs_displacement_ge_h"
                        ]
                    )

        if "PIC_ru" in arrays:
            ru = np.asarray(arrays["PIC_ru"], dtype=np.float64).reshape(-1)
            if ru.size != expected_count or not np.all(np.isfinite(ru)):
                raise ValueError(f"PIC_ru in {path} is invalid")
            row["max_ru_diagnostic"] = float(np.max(ru))

        initial_field_name = "initial_vertical_effective_stresses"
        initial_field_present = initial_field_name in arrays
        solver_initial_vertical_present.append(initial_field_present)
        initial_vertical: np.ndarray | None = None
        if initial_field_present:
            initial_vertical = np.asarray(
                arrays[initial_field_name], dtype=np.float64
            ).reshape(-1)
            if initial_vertical.size != expected_count or not np.all(
                np.isfinite(initial_vertical)
            ):
                raise ValueError(
                    f"Initial vertical skeleton stress field in {path} is invalid"
                )
            if reference_vertical_stress is None:
                reference_vertical_stress = initial_vertical.copy()
                reference_time = 0.0
                reference_basis = (
                    "checkpoint_initial_vertical_effective_stresses"
                )
            elif reference_basis == "checkpoint_initial_vertical_effective_stresses":
                if not np.allclose(
                    initial_vertical,
                    reference_vertical_stress,
                    rtol=1.0e-12,
                    atol=1.0e-8,
                ):
                    raise ValueError(
                        "Initial vertical skeleton stress field changed between "
                        "saved frames"
                    )
            else:
                raise ValueError(
                    "Initial vertical stress field appears only after a first-frame "
                    "reference fallback"
                )
        elif reference_basis == "checkpoint_initial_vertical_effective_stresses":
            raise ValueError(
                "Initial vertical skeleton stress field is missing from a saved frame"
            )

        stress_ratios: np.ndarray | None = None
        stress_loss_mask: np.ndarray | None = None
        solver_ratio_name = "vertical_effective_stress_remaining_ratios"
        solver_ratio_present.append(solver_ratio_name in arrays)
        if "stresses" in arrays:
            stresses = np.asarray(arrays["stresses"], dtype=np.float64)
            if (
                stresses.ndim != 2
                or stresses.shape[0] != expected_count
                or stresses.shape[1] < 3
                or not np.all(np.isfinite(stresses))
            ):
                raise ValueError(f"Skeleton stresses in {path} are invalid")
            vertical_stress = stresses[:, 1]
            if reference_vertical_stress is None:
                reference_vertical_stress = vertical_stress.copy()
                reference_time = row["time_s"]
                reference_basis = "first_saved_dynamic_stress_fallback"
            stress_ratios, reference_eligible, stress_loss_mask = (
                vertical_stress_loss_state(
                    vertical_stress, reference_vertical_stress
                )
            )
            threshold_masks[stress_threshold_name].append(stress_loss_mask)
            row["stress_reference_eligible_fraction"] = float(
                np.mean(reference_eligible)
            )
            if np.any(reference_eligible):
                row["minimum_vertical_stress_remaining_ratio"] = float(
                    np.min(stress_ratios[reference_eligible])
                )
            if solver_ratio_name in arrays:
                solver_ratios = np.asarray(
                    arrays[solver_ratio_name], dtype=np.float64
                ).reshape(-1)
                if solver_ratios.size != expected_count or not np.all(
                    np.isfinite(solver_ratios)
                ):
                    raise ValueError(f"Solver R_sigma field in {path} is invalid")
                difference = np.abs(
                    solver_ratios[reference_eligible]
                    - stress_ratios[reference_eligible]
                )
                max_difference = float(np.max(difference))
                solver_ratio_max_abs_differences.append(max_difference)
                if not np.allclose(
                    solver_ratios[reference_eligible],
                    stress_ratios[reference_eligible],
                    rtol=1.0e-10,
                    atol=1.0e-12,
                ):
                    raise ValueError(
                        f"Solver and Python R_sigma disagree in {path}; "
                        f"max |difference|={max_difference:.6g}"
                    )
            add_frame_threshold_metrics(
                row,
                stress_threshold_name,
                stress_loss_mask,
                particle_area,
                current_volumes,
            )
            if support_roi is not None:
                add_frame_threshold_metrics(
                    row,
                    f"support_ROI_{stress_threshold_name}",
                    stress_loss_mask,
                    particle_area,
                    current_volumes,
                    support_roi,
                )
        else:
            threshold_masks[stress_threshold_name].append(None)

        seepage_required = (
            "liquid_seepage_forces",
            "liquid_densities",
            "gamma_sub",
        )
        seepage_present = [name in arrays for name in seepage_required]
        seepage_index: np.ndarray | None = None
        seepage_mask: np.ndarray | None = None
        if any(seepage_present) and not all(seepage_present):
            missing = [
                name
                for name, present in zip(seepage_required, seepage_present)
                if not present
            ]
            raise ValueError(
                f"Cannot compute upward seepage IF in {path}; missing {missing}"
            )
        if all(seepage_present):
            seepage_index = upward_excess_seepage_force_index(
                arrays["liquid_seepage_forces"],
                arrays["liquid_densities"],
                arrays["gamma_sub"],
                gravity,
            )
            if seepage_index.size != expected_count:
                raise ValueError(f"Upward seepage IF count in {path} is invalid")
            seepage_mask = (
                seepage_index
                >= UPWARD_SEEPAGE_FORCE_THRESHOLD - THRESHOLD_RATIO_TOLERANCE
            )
            threshold_masks[seepage_threshold_name].append(seepage_mask)
            row["maximum_upward_seepage_IF"] = float(np.max(seepage_index))
            add_frame_threshold_metrics(
                row,
                seepage_threshold_name,
                seepage_mask,
                particle_area,
                current_volumes,
            )
            if support_roi is not None:
                add_frame_threshold_metrics(
                    row,
                    f"support_ROI_{seepage_threshold_name}",
                    seepage_mask,
                    particle_area,
                    current_volumes,
                    support_roi,
                )
        else:
            threshold_masks[seepage_threshold_name].append(None)

        joint_mask: np.ndarray | None = None
        if stress_loss_mask is not None and seepage_mask is not None:
            joint_mask = stress_loss_mask & seepage_mask
            threshold_masks[joint_threshold_name].append(joint_mask)
            add_frame_threshold_metrics(
                row,
                joint_threshold_name,
                joint_mask,
                particle_area,
                current_volumes,
            )
            if support_roi is not None:
                add_frame_threshold_metrics(
                    row,
                    f"support_ROI_{joint_threshold_name}",
                    joint_mask,
                    particle_area,
                    current_volumes,
                    support_roi,
                )
        else:
            threshold_masks[joint_threshold_name].append(None)

        if "liquefaction_potentials" in arrays:
            legacy_potential = np.asarray(
                arrays["liquefaction_potentials"], dtype=np.float64
            ).reshape(-1)
            if legacy_potential.size != expected_count or not np.all(
                np.isfinite(legacy_potential)
            ):
                raise ValueError(f"Legacy stress-loss potential in {path} is invalid")
            row["max_legacy_stress_loss_potential_QA"] = float(
                np.max(legacy_potential)
            )
            row["fraction_legacy_stress_loss_potential_ge_1_QA"] = float(
                np.mean(legacy_potential >= 1.0)
            )
        if "momentary_liquefied" in arrays:
            momentary_mask = np.asarray(
                arrays["momentary_liquefied"] >= 0.5
            ).reshape(-1)
            if momentary_mask.size != expected_count:
                raise ValueError(f"Legacy momentary flag in {path} is invalid")
            row["fraction_legacy_momentary_flag_QA"] = float(
                np.mean(momentary_mask)
            )
        pressure_name = (
            "PIC_pore_pressures"
            if "PIC_pore_pressures" in arrays
            else "PIC_liquid_pressures"
            if "PIC_liquid_pressures" in arrays
            else None
        )
        pressure_excess: np.ndarray | None = None
        if "PIC_pore_pressure_excess" in arrays:
            pressure_excess = np.asarray(
                arrays["PIC_pore_pressure_excess"], dtype=np.float64
            ).reshape(-1)
            if pressure_excess.size != expected_count or not np.all(
                np.isfinite(pressure_excess)
            ):
                raise ValueError(f"PIC pore-pressure excess in {path} is invalid")
        elif pressure_name is not None and checkpoint_pressure is not None:
            pressure_excess = (
                np.asarray(arrays[pressure_name], dtype=np.float64).reshape(-1)
                - checkpoint_pressure
            )
        for probe in probes:
            index = probe.index
            if pressure_name:
                row[f"{probe.name}_pressure_pa"] = float(arrays[pressure_name][index])
            if pressure_excess is not None:
                row[f"{probe.name}_pressure_excess_pa"] = float(
                    pressure_excess[index]
                )
            if "PIC_ru" in arrays:
                row[f"{probe.name}_ru"] = float(arrays["PIC_ru"][index])
            if "liquefaction_potentials" in arrays:
                row[f"{probe.name}_legacy_stress_loss_potential_QA"] = float(
                    arrays["liquefaction_potentials"][index]
                )
            if "stresses" in arrays:
                compression_sign = (
                    1.0
                    if reference_vertical_stress is not None
                    and reference_vertical_stress[index] > 0.0
                    else -1.0
                )
                mean, deviator = stress_invariants(
                    arrays["stresses"][index], compression_sign
                )
                row[f"{probe.name}_p_effective_pa"] = mean
                row[f"{probe.name}_q_pa"] = deviator
                row[f"{probe.name}_vertical_skeleton_stress_pa"] = float(
                    arrays["stresses"][index][1]
                )
                if stress_ratios is not None and reference_eligible is not None:
                    if reference_eligible[index]:
                        row[f"{probe.name}_vertical_stress_remaining_ratio"] = (
                            float(stress_ratios[index])
                        )
            if seepage_index is not None:
                row[f"{probe.name}_upward_seepage_IF"] = float(
                    seepage_index[index]
                )
            for state_name in ("eps_p_q", "void_ratio", "pdstrain"):
                if state_name in arrays:
                    row[f"{probe.name}_{state_name}"] = float(
                        arrays[state_name][index]
                    )
        rows.append(row)
    if not rows:
        raise FileNotFoundError("No particle VTP output was found")

    threshold_summaries: dict[str, dict[str, Any]] = {}
    frame_times = np.asarray([row["time_s"] for row in rows], dtype=np.float64)
    volumes_present = [volumes is not None for volumes in volume_frames]
    if any(volumes_present) and not all(volumes_present):
        raise ValueError("Particle volumes are missing from some saved particle frames")
    stacked_volumes = (
        np.stack([volumes for volumes in volume_frames if volumes is not None])
        if all(volumes_present)
        else None
    )
    if any(particle_ids_present) and not all(particle_ids_present):
        raise ValueError("Particle ids are missing from some saved particle frames")
    if any(displacement_frames_present) and not all(displacement_frames_present):
        raise ValueError(
            "Particle displacements are missing from some saved particle frames"
        )
    if any(solver_ratio_present) and not all(solver_ratio_present):
        raise ValueError(
            "Solver vertical-stress remaining ratios are missing from some frames"
        )
    stacked_threshold_masks: dict[str, np.ndarray] = {}
    for name, saved_masks in threshold_masks.items():
        present = [mask is not None for mask in saved_masks]
        if any(present) and not all(present):
            raise ValueError(
                f"Threshold field {name} is missing from some saved particle frames"
            )
        if all(present):
            all_masks = np.stack(
                [mask for mask in saved_masks if mask is not None]
            )
            stacked_threshold_masks[name] = all_masks
            metrics = threshold_spatiotemporal_metrics(
                frame_times,
                all_masks,
                particle_area,
                stacked_volumes,
            )
            metrics.update(
                threshold_release_segments(
                    frame_times,
                    all_masks,
                    particle_area,
                    release_time,
                    stacked_volumes,
                )
            )
            threshold_summaries[name] = metrics
            if support_roi is not None:
                support_volumes = (
                    stacked_volumes[:, support_roi]
                    if stacked_volumes is not None
                    else None
                )
                support_metrics = threshold_spatiotemporal_metrics(
                    frame_times,
                    all_masks[:, support_roi],
                    particle_area,
                    support_volumes,
                )
                support_metrics.update(
                    threshold_release_segments(
                        frame_times,
                        all_masks[:, support_roi],
                        particle_area,
                        release_time,
                        support_volumes,
                    )
                )
                threshold_summaries[f"support_ROI_{name}"] = support_metrics

    materials = config.get("materials", [])
    wave_period = (
        float(materials[1].get("wave_period", 0.0))
        if len(materials) > 1
        else 0.0
    )
    trigger_timing: dict[str, Any] = {}
    if (
        stress_threshold_name in stacked_threshold_masks
        and seepage_threshold_name in stacked_threshold_masks
        and wave_period > 0.0
    ):
        trigger_timing = trigger_to_stress_loss_timing(
            frame_times,
            stacked_threshold_masks[seepage_threshold_name],
            stacked_threshold_masks[stress_threshold_name],
            wave_period,
        )
        if support_roi is not None:
            support_timing = trigger_to_stress_loss_timing(
                frame_times,
                stacked_threshold_masks[seepage_threshold_name][:, support_roi],
                stacked_threshold_masks[stress_threshold_name][:, support_roi],
                wave_period,
            )
            trigger_timing.update(
                {f"support_ROI_{key}": value for key, value in support_timing.items()}
            )

    material = config["materials"][0]["type"]
    criterion_area_basis = (
        "current_particle_volumes"
        if stacked_volumes is not None
        else "nominal_reference_particle_area_fallback"
    )
    summary: dict[str, Any] = {
        "material": material,
        "particle_count": expected_count,
        "frames": len(rows),
        "last_particle_time_s": rows[-1]["time_s"],
        "nominal_reference_area_per_particle_m2": particle_area,
        "total_nominal_reference_area_m2": particle_area * expected_count,
        "current_particle_volumes_available": stacked_volumes is not None,
        "criterion_area_basis": criterion_area_basis,
        "liquefaction_area_basis": criterion_area_basis,
        "liquefaction_area_basis_compatibility_note": (
            "legacy alias of criterion_area_basis for the two primary "
            "stress-loss and seepage-force criteria; ru has no area metric"
        ),
        "primary_liquefaction_criteria": [
            stress_threshold_name,
            seepage_threshold_name,
        ],
        "joint_liquefaction_criterion": joint_threshold_name,
        "joint_liquefaction_definition": (
            "same particle at the same saved physical time satisfies both "
            "IF>=1 and eligible R_sigma<=0.05"
        ),
        "stress_reference_basis": reference_basis,
        "stress_reference_time_s": reference_time,
        "stress_reference_fallback_used": (
            reference_basis == "first_saved_dynamic_stress_fallback"
        ),
        "minimum_abs_reference_vertical_stress_pa": (
            MINIMUM_REFERENCE_VERTICAL_STRESS_PA
        ),
        "stress_loss_remaining_ratio_threshold": STRESS_LOSS_RATIO_THRESHOLD,
        "stress_loss_ratio_numerical_tolerance": THRESHOLD_RATIO_TOLERANCE,
        "solver_initial_vertical_effective_stress_field_used": (
            reference_basis == "checkpoint_initial_vertical_effective_stresses"
        ),
        "solver_initial_vertical_effective_stress_field_available_all_frames": (
            bool(solver_initial_vertical_present)
            and all(solver_initial_vertical_present)
        ),
        "particle_id_basis": (
            "validated_unique_ids_reordered_to_initial_0_to_N_minus_1"
            if particle_ids_present and all(particle_ids_present)
            else "VTK_order_fallback_no_ids"
        ),
        "solver_Rsigma_QA_available_all_frames": (
            bool(solver_ratio_present) and all(solver_ratio_present)
        ),
        "mesh_crossing_displacement_fields_available_all_frames": (
            bool(displacement_frames_present)
            and all(displacement_frames_present)
        ),
        "mesh_crossing_displacement_interpretation": (
            "Material-point displacement relative to the background cell width "
            "documents mesh-crossing / mesh-distortion burden. It is not a "
            "strain measure and does not imply that FEM is impossible."
        ),
    }
    displacement_rows = [
        row for row in rows if "soil_max_abs_displacement_over_h" in row
    ]
    if displacement_rows:
        summary.update(
            {
                "mesh_crossing_displacement_available": True,
                "mesh_crossing_cell_size_m": cell_size,
                "mesh_crossing_cell_size_source": "config.mesh.cellsize_min",
                "mesh_crossing_release_time_s": release_time,
                "post_release_displacement_saved_frames": len(
                    displacement_rows
                ),
                "post_release_max_soil_abs_displacement_m": max(
                    row["soil_max_abs_displacement_m"]
                    for row in displacement_rows
                ),
                "post_release_max_soil_abs_displacement_over_h": max(
                    row["soil_max_abs_displacement_over_h"]
                    for row in displacement_rows
                ),
                "post_release_max_soil_fraction_abs_displacement_ge_h": max(
                    row["soil_fraction_abs_displacement_ge_h"]
                    for row in displacement_rows
                ),
            }
        )
        summary["post_release_soil_crossed_one_background_cell"] = bool(
            summary["post_release_max_soil_abs_displacement_over_h"] >= 1.0
        )
        support_displacement_rows = [
            row
            for row in displacement_rows
            if "support_ROI_soil_max_abs_displacement_over_h" in row
        ]
        if support_displacement_rows:
            summary[
                "post_release_max_support_ROI_soil_abs_displacement_m"
            ] = max(
                row["support_ROI_soil_max_abs_displacement_m"]
                for row in support_displacement_rows
            )
            summary[
                "post_release_max_support_ROI_soil_abs_displacement_over_h"
            ] = max(
                row["support_ROI_soil_max_abs_displacement_over_h"]
                for row in support_displacement_rows
            )
            summary[
                "post_release_max_support_ROI_soil_fraction_abs_displacement_ge_h"
            ] = max(
                row["support_ROI_soil_fraction_abs_displacement_ge_h"]
                for row in support_displacement_rows
            )
            summary[
                "post_release_support_ROI_soil_crossed_one_background_cell"
            ] = bool(
                summary[
                    "post_release_max_support_ROI_soil_abs_displacement_over_h"
                ]
                >= 1.0
            )
    else:
        summary["mesh_crossing_displacement_available"] = False
        summary["mesh_crossing_displacement_unavailable_reason"] = (
            "saved particle frames do not reach pipeline release time"
            if displacement_frames_present and all(displacement_frames_present)
            else "displacements were not written to particle VTP"
        )
    if solver_ratio_max_abs_differences:
        summary["solver_Rsigma_QA_max_abs_difference"] = max(
            solver_ratio_max_abs_differences
        )
    diagnostic_maxima = {
        "max_ru_diagnostic": "max_ru_diagnostic",
        "max_upward_seepage_IF": "maximum_upward_seepage_IF",
        "max_legacy_stress_loss_potential_QA": (
            "max_legacy_stress_loss_potential_QA"
        ),
        "max_fraction_legacy_momentary_flag_QA": (
            "fraction_legacy_momentary_flag_QA"
        ),
    }
    for output_name, row_name in diagnostic_maxima.items():
        values = [row[row_name] for row in rows if row_name in row]
        if values:
            summary[output_name] = max(values)
    ratio_values = [
        row["minimum_vertical_stress_remaining_ratio"]
        for row in rows
        if "minimum_vertical_stress_remaining_ratio" in row
    ]
    if ratio_values:
        summary["minimum_vertical_stress_remaining_ratio"] = min(ratio_values)
    if reference_eligible is not None:
        summary["stress_reference_eligible_fraction"] = float(
            np.mean(reference_eligible)
        )
        summary["stress_reference_eligible_particles"] = int(
            np.count_nonzero(reference_eligible)
        )
        if not np.any(reference_eligible):
            raise ValueError(
                "No particle has an eligible initial vertical skeleton stress"
            )
    if support_roi is not None:
        summary["support_ROI_particles"] = int(np.count_nonzero(support_roi))
        summary["support_ROI_nominal_reference_area_m2"] = float(
            np.count_nonzero(support_roi) * particle_area
        )
    if stacked_volumes is not None:
        total_current_areas = np.sum(stacked_volumes, axis=1)
        summary["minimum_total_current_area_m2"] = float(
            np.min(total_current_areas)
        )
        summary["maximum_total_current_area_m2"] = float(
            np.max(total_current_areas)
        )
        if support_roi is not None:
            support_current_areas = np.sum(
                stacked_volumes[:, support_roi], axis=1
            )
            summary["support_ROI_minimum_total_current_area_m2"] = float(
                np.min(support_current_areas)
            )
            summary["support_ROI_maximum_total_current_area_m2"] = float(
                np.max(support_current_areas)
            )
    for threshold_name, metrics in threshold_summaries.items():
        for metric_name, value in metrics.items():
            summary[f"{threshold_name}_{metric_name}"] = value
    summary.update(trigger_timing)
    return rows, summary


def phase_metrics(
    rows: list[dict[str, float]], probes: list[Probe], config: dict[str, Any]
) -> dict[str, Any]:
    liquid = config["materials"][1]
    period = float(liquid["wave_period"])
    ramp = float(liquid.get("wave_ramp_time", period))
    release = float(config["analysis"]["rigid_pipeline"].get("release_time", 3 * period))
    fit_end = min(release, rows[-1]["time_s"])
    fit_start = max(ramp, fit_end - 2.0 * period)
    selected = [row for row in rows if fit_start <= row["time_s"] <= fit_end]
    if len(selected) < 6:
        return {
            "available": False,
            "reason": "fewer than six frames before pipeline release",
            "fit_start_s": fit_start,
            "fit_end_s": fit_end,
        }
    times = np.asarray([row["time_s"] for row in selected])
    surface_key = "surface_pressure_excess_pa"
    if surface_key not in selected[0]:
        return {
            "available": False,
            "reason": "checkpoint-relative pressure was not available in VTP analysis",
        }
    _, surface_amplitude, surface_phase = harmonic_fit(
        times, np.asarray([row[surface_key] for row in selected]), period
    )
    result: dict[str, Any] = {
        "available": True,
        "fit_start_s": fit_start,
        "fit_end_s": fit_end,
        "surface_amplitude_pa": surface_amplitude,
        "probes": {},
    }
    for probe in probes:
        key = f"{probe.name}_pressure_excess_pa"
        if key not in selected[0]:
            continue
        mean, amplitude, phase = harmonic_fit(
            times, np.asarray([row[key] for row in selected]), period
        )
        result["probes"][probe.name] = {
            "mean_pressure_pa": mean,
            "amplitude_pa": amplitude,
            "amplitude_ratio_to_surface": (
                amplitude / surface_amplitude if surface_amplitude > 0.0 else None
            ),
            "phase_lag_deg": wrap_degrees(phase - surface_phase),
        }
    return result


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        return
    keys: list[str] = []
    for row in rows:
        for key in row:
            if key not in keys:
                keys.append(key)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)


def flatten_summary(summary: dict[str, Any]) -> dict[str, Any]:
    output: dict[str, Any] = {"case": summary["case"], "uuid": summary["uuid"]}
    for group in ("pipeline", "particles"):
        for key, value in summary[group].items():
            if isinstance(value, (str, int, float, bool)):
                output[f"{group}_{key}"] = value
    phase = summary["phase"]
    if phase.get("available"):
        for probe, values in phase["probes"].items():
            output[f"phase_{probe}_lag_deg"] = values["phase_lag_deg"]
            output[f"phase_{probe}_amplitude_ratio"] = values[
                "amplitude_ratio_to_surface"
            ]
    return output


def _resume_dependency_fingerprints(summary: dict[str, Any]) -> tuple[str, str]:
    runtime = summary["completion_audit"].get("runtime_dependencies", {})
    resume = runtime.get("resume_equilibrium")
    if not isinstance(resume, dict):
        raise ValueError(
            f"{summary['case']} completion audit has no resume-equilibrium fingerprint"
        )
    return (
        str(resume["checkpoint_hdf5"]["sha256"]),
        str(resume["qa_vtp"]["sha256"]),
    )


def _pressure_replay_grid(summary: dict[str, Any]) -> tuple[Any, ...]:
    runtime = summary["completion_audit"].get("runtime_dependencies", {})
    pressure = runtime.get("read_pressure_database")
    if not isinstance(pressure, dict) or not isinstance(pressure.get("header"), dict):
        raise ValueError(
            f"{summary['case']} completion audit has no pressure-replay header"
        )
    header = pressure["header"]
    return tuple(
        header.get(name)
        for name in (
            "format_version",
            "particle_count",
            "step_interval",
            "max_step",
            "source_dt_s",
            "frame_count",
        )
    )


def _phase_comparison_config(config_path: Path) -> dict[str, Any]:
    """Strip only the pre-registered RL/RE identity and forcing-source fields."""

    config = deepcopy(json.loads(config_path.read_text(encoding="utf-8")))
    config.pop("title", None)
    analysis = config["analysis"]
    analysis.pop("uuid", None)
    pressure = analysis["prescribed_phase_pressures"]
    pressure.pop("path", None)
    pressure.pop("file_prefix", None)
    return config


def validate_RL_RE_comparison(summaries: Iterable[dict[str, Any]]) -> dict[str, Any] | None:
    """Require RL and RE to share the exact dynamic comparison basis."""

    by_code: dict[str, dict[str, Any]] = {}
    for summary in summaries:
        stem = Path(summary["config"]).stem
        if stem.endswith("_RL"):
            by_code["RL"] = summary
        elif stem.endswith("_RE"):
            by_code["RE"] = summary
    if not by_code:
        return None
    if set(by_code) != {"RL", "RE"}:
        raise ValueError("RL/RE phase comparison requires both complete cases")
    rl, re_case = by_code["RL"], by_code["RE"]
    rl_audit = rl["completion_audit"]
    re_audit = re_case["completion_audit"]
    checks = {
        "dt_s": (
            float(rl["analysis_dt_s"]),
            float(re_case["analysis_dt_s"]),
        ),
        "particle_steps": (
            rl_audit["expected_particle_steps"],
            re_audit["expected_particle_steps"],
        ),
        "particle_times_s": (
            rl_audit["expected_particle_times_s"],
            re_audit["expected_particle_times_s"],
        ),
        "particle_count": (
            int(rl["particles"]["particle_count"]),
            int(re_case["particles"]["particle_count"]),
        ),
        "resume_checkpoint_fingerprints": (
            _resume_dependency_fingerprints(rl),
            _resume_dependency_fingerprints(re_case),
        ),
        "pressure_replay_grid": (
            _pressure_replay_grid(rl),
            _pressure_replay_grid(re_case),
        ),
        "config_except_registered_forcing_source": (
            _phase_comparison_config(Path(rl["config"])),
            _phase_comparison_config(Path(re_case["config"])),
        ),
    }
    mismatches = [name for name, (left, right) in checks.items() if left != right]
    if mismatches:
        raise ValueError(
            "RL/RE comparison basis does not match for: " + ", ".join(mismatches)
        )
    return {
        "schema": "RL-RE-comparison-grid-v1",
        "RL_config_sha256": rl_audit["config_sha256"],
        "RE_config_sha256": re_audit["config_sha256"],
        "dt_s": checks["dt_s"][0],
        "particle_steps": checks["particle_steps"][0],
        "particle_times_s": checks["particle_times_s"][0],
        "particle_count": checks["particle_count"][0],
        "resume_checkpoint_fingerprints": checks[
            "resume_checkpoint_fingerprints"
        ][0],
        "non_forcing_config_identical": True,
    }


def analyse_case(config_path: Path, output_directory: Path) -> dict[str, Any]:
    config_path = config_path.resolve()
    case_root = config_path.parent
    while not (case_root / "prepare_study.py").exists() and case_root != case_root.parent:
        case_root = case_root.parent
    if not (case_root / "prepare_study.py").exists():
        raise FileNotFoundError(f"Cannot find case root above {config_path}")
    config = json.loads(config_path.read_text(encoding="utf-8"))
    uuid = str(config["analysis"]["uuid"])
    result_base = Path(config["post_processing"].get("path", "results"))
    if not result_base.is_absolute():
        result_base = case_root / result_base
    result_directory = result_base / uuid
    if not result_directory.exists():
        raise FileNotFoundError(f"Result directory does not exist: {result_directory}")

    particle_location = resolve_particle_input_path(config, case_root)
    coordinates = read_initial_coordinates(particle_location)
    probes = build_probes(config, coordinates)
    vtp_files = sorted(result_directory.glob("particle*.vtp"), key=particle_step)
    completion_audit = validate_completed_result(
        config_path,
        case_root,
        config,
        result_directory,
        vtp_files,
    )
    pressure_reference, pressure_reference_audit = checkpoint_pressure_reference(
        completion_audit, case_root, len(coordinates)
    )
    particle_rows, particle_summary = particle_histories(
        vtp_files,
        config,
        coordinates,
        probes,
        checkpoint_pressure=pressure_reference,
    )
    particle_summary["pore_pressure_excess_reference"] = pressure_reference_audit
    history = read_pipeline_history(result_directory)
    completion_audit.update(validate_pipeline_history_grid(history, config))
    pipeline = config["analysis"]["rigid_pipeline"]
    diameter = 2.0 * float(pipeline["outer_radius"])
    bed_elevation = float(config["materials"][1]["sea_level"]) - float(
        config["materials"][1]["depth_left"]
    )
    cover = bed_elevation - (
        float(pipeline["center"][1]) + float(pipeline["outer_radius"])
    )
    release = float(pipeline.get("release_time", 0.0))
    period = float(config["materials"][1]["wave_period"])
    pipeline_summary, cycle_rows = pipeline_metrics(
        history, diameter, release, period, cover
    )
    summary = {
        "case": config_path.stem,
        "uuid": uuid,
        "config": str(config_path),
        "result_directory": str(result_directory),
        "analysis_dt_s": float(config["analysis"]["dt"]),
        "completion_audit": completion_audit,
        "definitions": {
            "ru": (
                "excess pore-pressure ratio retained only as a diagnostic "
                "probe/history maximum; it is not a liquefaction criterion"
            ),
            "vertical_stress_loss": (
                "R_sigma = sigma'_yy(t) / sigma'_yy(0); the checkpoint-resumed "
                "initial_vertical_effective_stresses field is preferred and "
                "represents dynamic t=0 and must remain constant, with first-"
                "frame stress (at that saved frame's time) an explicit legacy "
                "fallback; particles require |sigma'_yy(0)| >= 100 Pa and the "
                "primary threshold is R_sigma <= 0.05 + 1e-12"
            ),
            "upward_seepage_IF": (
                "IF = dot(f_liquid_seepage + rho_l*g, -g/|g|) / gamma_sub; "
                "the hydrostatic component is removed, negative/downward "
                "values are clipped to zero, and the primary threshold is IF >= 1"
            ),
            "wave_excess_pore_pressure": (
                "PIC pore pressure relative to the checkpoint-resumed initial "
                "pore pressure for the same particle ID; solver-saved excess is "
                "preferred and the audited equilibrium VTP is the fallback"
            ),
            "joint_liquefaction": (
                "same particle at the same saved physical time satisfies "
                "IF>=1 and eligible R_sigma<=0.05; first-event ordering and "
                "the preceding-wave-cycle association are reported separately"
            ),
            "legacy_solver_fields": (
                "liquefaction_potentials and momentary_liquefied are retained "
                "only as legacy stress-loss QA; neither is labelled seepage LI"
            ),
            "liquefied_area": (
                "sum of per-frame particle volumes (2D unit-thickness current "
                "area) when VTP volumes are available; otherwise threshold-"
                "particle count times the initial Cartesian spacing area is "
                "reported as an explicit nominal reference-area fallback"
            ),
            "time_integrated_area": (
                "trapezoidal integral between saved VTP frames, in m2 s"
            ),
            "resolved_duration": (
                "conservative lower bound: an interval counts only when the "
                "threshold is met at both saved endpoints"
            ),
            "large_deformation": "maximum |vertical pipeline displacement| / D >= 0.5",
            "mesh_crossing_displacement": (
                "for saved frames at or after pipeline release, soil |u| is "
                "normalised by config mesh.cellsize_min = h; max |u|/h and "
                "the fraction with |u| >= h are reported for the full bed and "
                "the pre-registered initial support-ROI cohort. This is evidence "
                "of material-point mesh crossing / mesh-distortion burden, not "
                "a strain measure or a claim that FEM is impossible"
            ),
        },
        "probes": [probe.__dict__ for probe in probes],
        "pipeline": pipeline_summary,
        "particles": particle_summary,
        "phase": phase_metrics(particle_rows, probes, config),
        "method_evidence": {
            "pipeline_maximum_abs_vertical_displacement_over_D": (
                pipeline_summary["maximum_abs_vertical_displacement_over_D"]
            ),
            "pipeline_contact_topology_change_demonstrated": (
                pipeline_summary["contact_topology_change_demonstrated"]
            ),
            "pipeline_large_deformation_or_topology_gate": (
                pipeline_summary["large_deformation_demonstrated"]
            ),
            "soil_mesh_crossing_displacement_available": particle_summary[
                "mesh_crossing_displacement_available"
            ],
            "soil_full_domain_crossed_one_background_cell_after_release": (
                particle_summary.get(
                    "post_release_soil_crossed_one_background_cell"
                )
            ),
            "soil_support_ROI_crossed_one_background_cell_after_release": (
                particle_summary.get(
                    "post_release_support_ROI_soil_crossed_one_background_cell"
                )
            ),
            "interpretation": (
                "Soil mesh-crossing displacement is assessed in parallel with "
                "pipeline |uy|/D and the contact-topology gate when motivating "
                "MPM. It measures mesh-distortion burden rather than large "
                "strain, and it does not assert that an FEM treatment is "
                "impossible."
            ),
        },
    }
    output_directory.mkdir(parents=True, exist_ok=True)
    stem = config_path.stem
    (output_directory / f"{stem}_summary.json").write_text(
        json.dumps(summary, indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )
    write_csv(output_directory / f"{stem}_probe_history.csv", particle_rows)
    write_csv(output_directory / f"{stem}_pipeline_cycles.csv", cycle_rows)
    return summary


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("configs", nargs="+", type=Path)
    parser.add_argument("--output", type=Path, default=Path("analysis"))
    return parser


def main() -> int:
    args = build_parser().parse_args()
    summaries = [analyse_case(path, args.output) for path in args.configs]
    comparison_audit = validate_RL_RE_comparison(summaries)
    if comparison_audit is not None:
        (args.output / "RL_RE_comparison_grid.json").write_text(
            json.dumps(comparison_audit, indent=2, allow_nan=False) + "\n",
            encoding="utf-8",
        )
    write_csv(args.output / "study_summary.csv", [flatten_summary(s) for s in summaries])
    for summary in summaries:
        displacement = summary["pipeline"][
            "maximum_abs_vertical_displacement_over_D"
        ]
        status = "large deformation" if displacement >= 0.5 else "small deformation"
        print(
            f"{summary['case']}: max |uy|/D={displacement:.4g} ({status}); "
            f"min R_sigma="
            f"{summary['particles']['minimum_vertical_stress_remaining_ratio']:.4g}; "
            f"max upward seepage IF="
            f"{summary['particles']['max_upward_seepage_IF']:.4g}; "
            f"max ru (diagnostic)="
            f"{summary['particles'].get('max_ru_diagnostic', float('nan')):.4g}"
        )
    print(f"Wrote analysis tables to {args.output.resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
