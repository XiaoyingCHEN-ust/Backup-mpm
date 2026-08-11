#!/usr/bin/env python3
"""Fail-fast validation for generated pipeline study inputs and dependencies."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np

from analyze_study import read_vtp
import prepare_study as study
from phase_controls import (
    database_paths,
    iter_frames,
    read_header,
    read_points,
    validate_database,
)


CASE_ROOT = Path(__file__).resolve().parent
COMPLETION_SCHEMA = "pipeline-case-completion-v2"
COMPLETION_FILENAME = "pipeline_completion.json"
HDF5_SIGNATURE = b"\x89HDF\r\n\x1a\n"


def resolve_case_path(value: str | Path) -> Path:
    path = Path(value)
    return path if path.is_absolute() else CASE_ROOT / path


def first_count(path: Path) -> int:
    with path.open(encoding="utf-8") as stream:
        return int(stream.readline().strip())


def read_ascii_table(path: Path, columns: int) -> np.ndarray:
    """Read a counted MPM ASCII table and reject silently ignored columns."""
    rows: list[list[float]] = []
    declared_count: int | None = None
    with path.open(encoding="utf-8") as stream:
        for line in stream:
            stripped = line.strip()
            if not stripped or stripped.startswith(("#", "!")):
                continue
            if declared_count is None:
                declared_count = int(stripped)
                continue
            values = stripped.split()
            if len(values) != columns:
                raise ValueError(
                    f"{path} row has {len(values)} columns; expected {columns}: "
                    f"{stripped}"
                )
            rows.append([float(value) for value in values])
    if declared_count is None:
        raise ValueError(f"{path} has no declared row count")
    values = np.asarray(rows, dtype=float)
    if values.shape != (declared_count, columns):
        raise ValueError(
            f"{path} declares {declared_count} rows but contains {values.shape[0]}"
        )
    if not np.all(np.isfinite(values)):
        raise ValueError(f"{path} contains non-finite values")
    return values


def mesh_counts(path: Path) -> tuple[int, int]:
    with path.open(encoding="utf-8") as stream:
        for line in stream:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            values = stripped.split()
            if len(values) != 2:
                raise ValueError(f"Invalid mesh count row in {path}: {stripped}")
            return int(values[0]), int(values[1])
    raise ValueError(f"Mesh file has no node/cell count row: {path}")


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def checked_case_path(path: Path) -> Path:
    """Resolve a path and reject completion artifacts outside this case tree."""

    resolved = path.resolve()
    root = CASE_ROOT.resolve()
    try:
        resolved.relative_to(root)
    except ValueError as error:
        raise ValueError(f"Completion artifact escapes the case directory: {path}") from error
    return resolved


def case_relative(path: Path) -> str:
    return checked_case_path(path).relative_to(CASE_ROOT.resolve()).as_posix()


def artifact_audit(path: Path) -> dict[str, Any]:
    path = checked_case_path(path)
    if not path.is_file() or path.stat().st_size <= 0:
        raise FileNotFoundError(f"Completion artifact is missing or empty: {path}")
    return {
        "path": case_relative(path),
        "size_bytes": path.stat().st_size,
        "sha256": file_sha256(path),
    }


def result_directory(config: dict[str, Any]) -> Path:
    base = resolve_case_path(config["post_processing"].get("path", "results/"))
    return checked_case_path(base / str(config["analysis"]["uuid"]))


def completion_path(config: dict[str, Any]) -> Path:
    return result_directory(config) / COMPLETION_FILENAME


def prepare_output_base(config_path: Path) -> Path:
    """Create the configured output base recursively inside the case tree."""

    config_path = checked_case_path(config_path)
    config = json.loads(config_path.read_text(encoding="utf-8"))
    output_base = checked_case_path(
        resolve_case_path(config["post_processing"].get("path", "results/"))
    )
    if output_base.exists() and not output_base.is_dir():
        raise ValueError(f"Configured output path is not a directory: {output_base}")
    output_base.mkdir(parents=True, exist_ok=True)
    return output_base


def clear_case_completion(config_path: Path) -> Path:
    """Remove only this config's audited sentinel before a new solver run."""

    config_path = checked_case_path(config_path)
    config = json.loads(config_path.read_text(encoding="utf-8"))
    sentinel = completion_path(config)
    if sentinel.exists() or sentinel.is_symlink():
        sentinel = checked_case_path(sentinel)
        if sentinel.name != COMPLETION_FILENAME or not sentinel.is_file():
            raise ValueError(f"Refusing to clear unexpected completion path: {sentinel}")
        sentinel.unlink()
    return sentinel


def final_artifact_paths(config: dict[str, Any]) -> tuple[Path, Path | None]:
    nsteps = int(config["analysis"]["nsteps"])
    digits = len(str(nsteps))
    directory = result_directory(config)
    vtp = directory / f"particle{nsteps:0{digits}d}.vtp"
    hdf5 = None
    if config["post_processing"].get("write_hdf5"):
        hdf5 = directory / f"particles{nsteps:0{digits}d}.h5"
    return vtp, hdf5


def validate_final_vtp(path: Path) -> None:
    points, _ = read_vtp(checked_case_path(path))
    if len(points) == 0 or not np.all(np.isfinite(points)):
        raise ValueError(f"Final VTP has no finite particle coordinates: {path}")


def validate_final_hdf5(path: Path) -> None:
    path = checked_case_path(path)
    if not path.is_file() or path.stat().st_size <= len(HDF5_SIGNATURE):
        raise ValueError(f"Final HDF5 is missing or too short: {path}")
    with path.open("rb") as stream:
        if stream.read(len(HDF5_SIGNATURE)) != HDF5_SIGNATURE:
            raise ValueError(f"Final HDF5 has an invalid signature: {path}")


def written_pressure_database_audit(config: dict[str, Any]) -> dict[str, Any] | None:
    pressure = config["analysis"].get("prescribed_phase_pressures")
    if pressure is None or not bool(pressure.get("write")):
        return None
    directory = checked_case_path(resolve_case_path(pressure["path"]))
    prefix = str(pressure["file_prefix"])
    points_path, values_path = database_paths(directory, prefix)
    points_path = checked_case_path(points_path)
    values_path = checked_case_path(values_path)
    with values_path.open("rb") as stream:
        header = read_header(stream)
        for _ in iter_frames(stream, header):
            pass
    validate_database(values_path, header)
    points = read_points(points_path, header.dimension)
    if len(points.particle_ids) != header.particle_count:
        raise ValueError("Written pressure point and value counts differ")
    if header.step_interval != int(pressure["step_interval"]):
        raise ValueError("Written pressure database step interval differs from config")
    if header.max_step != int(pressure["max_step"]):
        raise ValueError("Written pressure database max step differs from config")
    if not math.isclose(
        header.source_dt, float(pressure["source_dt"]), rel_tol=0.0, abs_tol=1.0e-15
    ):
        raise ValueError("Written pressure database source dt differs from config")
    return {
        "points": artifact_audit(points_path),
        "values": artifact_audit(values_path),
        "header": {
            "format_version": header.format_version,
            "dimension": header.dimension,
            "particle_count": header.particle_count,
            "step_interval": header.step_interval,
            "max_step": header.max_step,
            "source_dt_s": header.source_dt,
            "frame_count": header.frame_count,
        },
    }


def read_pressure_database_audit(config: dict[str, Any]) -> dict[str, Any] | None:
    pressure = config["analysis"].get("prescribed_phase_pressures")
    if pressure is None or not bool(pressure.get("enable")):
        return None
    directory = checked_case_path(resolve_case_path(pressure["path"]))
    prefix = str(pressure["file_prefix"])
    points_path, values_path = database_paths(directory, prefix)
    points_path = checked_case_path(points_path)
    values_path = checked_case_path(values_path)
    with values_path.open("rb") as stream:
        header = read_header(stream)
        for _ in iter_frames(stream, header):
            pass
    validate_database(values_path, header)
    points = read_points(points_path, header.dimension)
    if len(points.particle_ids) != header.particle_count:
        raise ValueError("Read pressure point and value counts differ")
    if header.step_interval != int(pressure["step_interval"]):
        raise ValueError("Read pressure database step interval differs from config")
    if header.max_step < int(pressure["max_step"]):
        raise ValueError("Read pressure database ends before the configured replay")
    if not math.isclose(
        header.source_dt, float(pressure["source_dt"]), rel_tol=0.0, abs_tol=1.0e-15
    ):
        raise ValueError("Read pressure database source dt differs from config")
    audit: dict[str, Any] = {
        "points": artifact_audit(points_path),
        "values": artifact_audit(values_path),
        "header": {
            "format_version": header.format_version,
            "dimension": header.dimension,
            "particle_count": header.particle_count,
            "step_interval": header.step_interval,
            "max_step": header.max_step,
            "source_dt_s": header.source_dt,
            "frame_count": header.frame_count,
        },
    }
    if prefix == "phase_erased":
        metadata_path = checked_case_path(directory / f"{prefix}_metadata.json")
        audit["metadata"] = artifact_audit(metadata_path)
    return audit


def runtime_dependency_audit(config: dict[str, Any]) -> dict[str, Any]:
    """Fingerprint every external state consumed by this solver invocation."""

    dependencies: dict[str, Any] = {}
    resume = config["analysis"].get("resume", {})
    if bool(resume.get("resume")):
        checkpoint = checked_case_path(checkpoint_path(config))
        validate_final_hdf5(checkpoint)
        equilibrium_vtp = checkpoint.with_name(
            checkpoint.name.replace("particles", "particle")
        ).with_suffix(".vtp")
        validate_equilibrium_vtp(checkpoint)
        dependencies["resume_equilibrium"] = {
            "checkpoint_hdf5": artifact_audit(checkpoint),
            "qa_vtp": artifact_audit(equilibrium_vtp),
        }
    pressure = read_pressure_database_audit(config)
    if pressure is not None:
        dependencies["read_pressure_database"] = pressure
    return dependencies


def write_case_completion(config_path: Path) -> Path:
    """Validate final artifacts and atomically publish a restart sentinel."""

    config_path = checked_case_path(config_path)
    config = json.loads(config_path.read_text(encoding="utf-8"))
    final_vtp, final_hdf5 = final_artifact_paths(config)
    validate_final_vtp(final_vtp)
    artifacts: dict[str, Any] = {"final_vtp": artifact_audit(final_vtp)}
    if final_hdf5 is not None:
        validate_final_hdf5(final_hdf5)
        artifacts["final_hdf5"] = artifact_audit(final_hdf5)
    pressure_database = written_pressure_database_audit(config)
    if pressure_database is not None:
        artifacts["written_pressure_database"] = pressure_database
    runtime_dependencies = runtime_dependency_audit(config)

    payload = {
        "schema": COMPLETION_SCHEMA,
        "completed_at_utc": datetime.now(timezone.utc).isoformat(),
        "config": {
            "path": case_relative(config_path),
            "sha256": file_sha256(config_path),
            "uuid": str(config["analysis"]["uuid"]),
            "nsteps": int(config["analysis"]["nsteps"]),
        },
        "artifacts": artifacts,
        "runtime_dependencies": runtime_dependencies,
    }
    sentinel = completion_path(config)
    sentinel.parent.mkdir(parents=True, exist_ok=True)
    temporary = sentinel.with_name(f".{sentinel.name}.{os.getpid()}.tmp")
    try:
        with temporary.open("w", encoding="utf-8", newline="\n") as stream:
            json.dump(payload, stream, indent=2, sort_keys=True)
            stream.write("\n")
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, sentinel)
    finally:
        if temporary.exists():
            temporary.unlink()
    return sentinel


def validate_geometry(config: dict[str, Any]) -> list[str]:
    messages: list[str] = []
    mesh = config["mesh"]
    generator = config["particles"][0]["generator"]
    required = {
        "mesh": resolve_case_path(mesh["mesh"]),
        "entity sets": resolve_case_path(mesh["entity_sets"]),
        "particles": resolve_case_path(generator["location"]),
        "temperatures": resolve_case_path(mesh["particles_temperatures"]),
    }
    for description, path in required.items():
        if not path.is_file():
            raise FileNotFoundError(
                f"Missing {description} input {path}. Run mesh.py with matching "
                "resolution/mesh-dir arguments first."
            )
    node_count, cell_count = mesh_counts(required["mesh"])
    particle_count = first_count(required["particles"])
    temperature_count = first_count(required["temperatures"])
    if node_count <= 0 or particle_count <= 0:
        raise ValueError("Mesh/particle files declare no entities")
    if temperature_count != particle_count:
        raise ValueError(
            f"Temperature count {temperature_count} != particle count {particle_count}"
        )
    temperatures = read_ascii_table(required["temperatures"], 1)[:, 0]
    if not np.allclose(temperatures, 0.0, rtol=0.0, atol=1.0e-12):
        raise ValueError("Initial particle temperatures must be zero")

    initial_stress_name = mesh.get("particles_stresses")
    initial_pressure_spec = mesh.get("particles_pore_pressures")
    if bool(initial_stress_name) != bool(initial_pressure_spec):
        raise ValueError("Initial stress and liquid-pressure fields must be paired")
    if initial_stress_name:
        stress_path = resolve_case_path(initial_stress_name)
        pressure_path = resolve_case_path(initial_pressure_spec.get("file", ""))
        for description, path in (
            ("initial effective stresses", stress_path),
            ("initial liquid pressures", pressure_path),
        ):
            if not path.is_file():
                raise FileNotFoundError(f"Missing {description} input {path}")
        coordinates = read_ascii_table(required["particles"], 2)
        stresses = read_ascii_table(stress_path, 6)
        liquid_pressures = read_ascii_table(pressure_path, 1)[:, 0]
        if stresses.shape[0] != particle_count or liquid_pressures.size != particle_count:
            raise ValueError("Initial-field count differs from particle count")

        soil, fluid = config["materials"]
        gravity = 9.81
        liquid_density = float(fluid["density"])
        expected_pressure = liquid_density * gravity * np.maximum(
            float(fluid["sea_level"]) - coordinates[:, 1], 0.0
        )
        if not np.allclose(
            liquid_pressures, expected_pressure, rtol=0.0, atol=1.0e-6
        ):
            raise ValueError("Initial liquid-pressure file is not hydrostatic")

        porosity = float(soil["porosity"])
        liquid_saturation = float(fluid["liquid_saturation"])
        gas_saturation = float(fluid["gas_saturation"])
        gas_density = (
            float(fluid["gas_molar_mass"])
            * float(soil["p_ref"])
            / float(fluid["gas_constant"])
            / (273.15 + float(soil["initial_temperature"]))
        )
        mixture_density = (
            (1.0 - porosity) * float(soil["density"])
            + porosity
            * (
                liquid_saturation * liquid_density
                + gas_saturation * gas_density
            )
        )
        effective_unit_weight = (mixture_density - liquid_density) * gravity
        seabed = float(fluid["sea_level"]) - float(fluid["depth_left"])
        burial_depth = np.maximum(seabed - coordinates[:, 1], 0.0)
        expected_yy = -effective_unit_weight * burial_depth
        expected_xx = study.INITIAL_EFFECTIVE_K0 * expected_yy
        if not np.allclose(stresses[:, 1], expected_yy, rtol=0.0, atol=1.0e-6):
            raise ValueError("Initial vertical stress is not effective overburden")
        if not np.allclose(stresses[:, 0], expected_xx, rtol=0.0, atol=1.0e-6):
            raise ValueError("Initial horizontal effective stress has stale K0")
        if not np.allclose(stresses[:, 2:], 0.0, rtol=0.0, atol=1.0e-12):
            raise ValueError("Initial out-of-plane/shear stresses must be zero")
        messages.append(
            "initial state: hydrostatic liquid pressure + effective self-weight stress"
        )
    entity_sets = json.loads(required["entity sets"].read_text(encoding="utf-8"))
    node_ids = [int(value) for item in entity_sets["node_sets"] for value in item["set"]]
    particle_ids = [
        int(value) for item in entity_sets["particle_sets"] for value in item["set"]
    ]
    if not node_ids or min(node_ids) < 0 or max(node_ids) >= node_count:
        raise ValueError("Node entity set contains invalid IDs")
    if not particle_ids or min(particle_ids) < 0 or max(particle_ids) >= particle_count:
        raise ValueError("Particle entity set contains invalid IDs")

    geometry_path = required["particles"].parent / "pipeline_geometry.json"
    if geometry_path.is_file():
        geometry = json.loads(geometry_path.read_text(encoding="utf-8"))
        pipeline = config["analysis"]["rigid_pipeline"]
        for key in ("outer_radius", "wall_thickness", "density"):
            if not math.isclose(
                float(geometry[key]), float(pipeline[key]), rel_tol=0.0, abs_tol=1.0e-12
            ):
                raise ValueError(f"Pipeline {key} differs between mesh and config")
        if any(
            not math.isclose(float(a), float(b), rel_tol=0.0, abs_tol=1.0e-12)
            for a, b in zip(geometry["center"], pipeline["center"])
        ):
            raise ValueError("Pipeline center differs between mesh and config")
        if not math.isclose(
            float(geometry["cover_over_diameter"]),
            study.PIPE_COVER_RATIO,
            abs_tol=1.0e-12,
        ):
            raise ValueError("Geometry does not use the registered cover/D")
    messages.append(
        f"geometry: {node_count} nodes, {cell_count} cells, {particle_count} particles"
    )
    return messages


def checkpoint_path(config: dict[str, Any]) -> Path:
    resume = config["analysis"]["resume"]
    post_path = resolve_case_path(config["post_processing"].get("path", "results/"))
    digits = int(math.log10(int(resume["nsteps"]))) + 1
    filename = f"particles{int(resume['step']):0{digits}d}.h5"
    return post_path / str(resume["uuid"]) / filename


def validate_equilibrium_vtp(checkpoint: Path) -> str:
    vtp = checkpoint.with_name(checkpoint.name.replace("particles", "particle")).with_suffix(
        ".vtp"
    )
    if not vtp.is_file():
        raise FileNotFoundError(f"Equilibrium QA VTP is missing: {vtp}")
    _, arrays = read_vtp(vtp)
    for required in ("velocities", "displacements", "porosities"):
        if required not in arrays:
            raise ValueError(f"Equilibrium QA field {required} is absent in {vtp}")
    maximum_velocity = float(np.max(np.linalg.norm(arrays["velocities"][:, :2], axis=1)))
    maximum_displacement = float(
        np.max(np.linalg.norm(arrays["displacements"][:, :2], axis=1))
    )
    porosity_min = float(np.min(arrays["porosities"]))
    porosity_max = float(np.max(arrays["porosities"]))
    if maximum_velocity > 1.0e-3:
        raise ValueError(
            f"Equilibrium max velocity {maximum_velocity:.3e} m/s exceeds 1e-3 m/s"
        )
    if maximum_displacement > 2.0e-2:
        raise ValueError(
            f"Equilibrium max displacement {maximum_displacement:.3e} m exceeds one cell"
        )
    if porosity_min <= 0.0 or porosity_max >= 1.0:
        raise ValueError("Equilibrium porosity lies outside (0, 1)")
    return (
        f"equilibrium QA: max|v|={maximum_velocity:.3e} m/s, "
        f"max|u|={maximum_displacement:.3e} m, "
        f"n=[{porosity_min:.5f}, {porosity_max:.5f}]"
    )


def validate_pressure_dependency(config: dict[str, Any], runtime: bool) -> list[str]:
    messages: list[str] = []
    pressure = config["analysis"].get("prescribed_phase_pressures")
    if pressure is None:
        return messages
    read_mode = bool(pressure["enable"])
    write_mode = bool(pressure["write"])
    if read_mode == write_mode:
        raise ValueError("Pressure database must be in exactly one of read/write modes")
    directory = resolve_case_path(pressure["path"])
    prefix = str(pressure["file_prefix"])
    points_path, values_path = database_paths(directory, prefix)
    if write_mode:
        messages.append(f"pressure writer: {values_path}")
        return messages
    if not runtime:
        messages.append(f"pressure reader dependency (deferred): {values_path}")
        return messages
    if not points_path.is_file() or not values_path.is_file():
        raise FileNotFoundError(
            f"Pressure replay input is missing: {points_path} / {values_path}"
        )
    with values_path.open("rb") as stream:
        header = read_header(stream)
    validate_database(values_path, header)
    points = read_points(points_path, header.dimension)
    if len(points.particle_ids) != header.particle_count:
        raise ValueError("Pressure point and value counts differ")
    expected = config["analysis"]["prescribed_phase_pressures"]
    if header.step_interval != int(expected["step_interval"]):
        raise ValueError("Pressure database step interval differs from config")
    if header.max_step < int(expected["max_step"]):
        raise ValueError("Pressure database ends before the requested replay")
    if not math.isclose(header.source_dt, float(expected["source_dt"]), abs_tol=1e-15):
        raise ValueError("Pressure database source dt differs from config")
    if prefix == "phase_erased":
        metadata_path = directory / f"{prefix}_metadata.json"
        if not metadata_path.is_file():
            raise FileNotFoundError("Phase-erased pressure database lacks audit metadata")
        metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
        if metadata.get("schema") != "phase-erased-pressure-control-v1":
            raise ValueError("Phase-erased database has unexpected metadata schema")
        if metadata.get("classification") != "one-way numerical counterfactual; not fully coupled":
            raise ValueError("Phase-erased database has unexpected classification")
        recorded_output_points = checked_case_path(
            resolve_case_path(metadata["output"]["points"])
        )
        recorded_output_values = checked_case_path(
            resolve_case_path(metadata["output"]["values"])
        )
        if recorded_output_points != points_path.resolve() or recorded_output_values != values_path.resolve():
            raise ValueError("Phase-erased metadata points to unexpected output files")
        if metadata["output"]["points_sha256"] != file_sha256(points_path):
            raise ValueError("Phase-erased point-file hash differs from metadata")
        if metadata["output"]["values_sha256"] != file_sha256(values_path):
            raise ValueError("Phase-erased value-file hash differs from metadata")
        source_points = checked_case_path(
            resolve_case_path(metadata["source"]["points"])
        )
        source_values = checked_case_path(
            resolve_case_path(metadata["source"]["values"])
        )
        if metadata["source"]["points_sha256"] != file_sha256(source_points):
            raise ValueError("Current lagged point-file hash differs from phase metadata")
        if metadata["source"]["values_sha256"] != file_sha256(source_values):
            raise ValueError("Current lagged value-file hash differs from phase metadata")
        with source_values.open("rb") as stream:
            source_header = read_header(stream)
        validate_database(source_values, source_header)
        if source_header != header:
            raise ValueError("Lagged and phase-erased pressure headers differ")
        recorded_header = metadata["header"]
        expected_header = {
            "format_version": header.format_version,
            "dimension": header.dimension,
            "particle_count": header.particle_count,
            "step_interval": header.step_interval,
            "max_step": header.max_step,
        }
        for key, expected_value in expected_header.items():
            if recorded_header.get(key) != expected_value:
                raise ValueError(f"Phase metadata header field {key} is stale")
        if not math.isclose(
            float(recorded_header["source_dt_s"]),
            header.source_dt,
            rel_tol=0.0,
            abs_tol=1.0e-15,
        ):
            raise ValueError("Phase metadata source dt is stale")
    messages.append(
        f"pressure reader: {header.frame_count} frames, {header.particle_count} particles"
    )
    return messages


def validate_runtime_dependencies(config: dict[str, Any], runtime: bool) -> list[str]:
    messages: list[str] = []
    resume = config["analysis"].get("resume", {})
    if runtime and resume.get("resume"):
        checkpoint = checkpoint_path(config)
        if not checkpoint.is_file():
            raise FileNotFoundError(
                f"Required equilibrium checkpoint does not exist: {checkpoint}"
            )
        messages.append(f"checkpoint: {checkpoint}")
        messages.append(validate_equilibrium_vtp(checkpoint))
    messages.extend(validate_pressure_dependency(config, runtime))
    return messages


def validate_one(config_path: Path, runtime: bool) -> list[str]:
    config_path = config_path.resolve()
    config = json.loads(config_path.read_text(encoding="utf-8"))
    study.validate_config(config)
    messages = [f"config: {config_path}", f"uuid: {config['analysis']['uuid']}"]
    messages.extend(validate_geometry(config))
    messages.extend(validate_runtime_dependencies(config, runtime))
    material = config["materials"][0]
    fluid = config["materials"][1]
    messages.append(
        "physics: "
        f"{material['type']}, Sw={fluid['liquid_saturation']}, "
        f"kappa={material['intrinsic_permeability']:.3e} m^2"
    )
    return messages


def config_paths(arguments: argparse.Namespace) -> list[Path]:
    paths = [resolve_case_path(path) for path in arguments.configs]
    if arguments.manifest:
        manifest_path = resolve_case_path(arguments.manifest)
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        paths.extend(manifest_path.parent / case["config"] for case in manifest["cases"])
    if not paths:
        paths = sorted((CASE_ROOT / "configs" / "screen").glob("*.json"))
        paths = [path for path in paths if path.name != "manifest.json"]
    unique: list[Path] = []
    for path in paths:
        if path not in unique:
            unique.append(path)
    return unique


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("configs", nargs="*", type=Path)
    parser.add_argument("--manifest", type=Path)
    parser.add_argument(
        "--runtime",
        action="store_true",
        help="also require checkpoints and pressure databases needed at launch",
    )
    parser.add_argument(
        "--write-completion",
        action="store_true",
        help="after runtime validation, audit final outputs and atomically publish completion",
    )
    parser.add_argument(
        "--clear-completion",
        action="store_true",
        help="after validation, remove only this case's old completion before rerunning",
    )
    parser.add_argument(
        "--prepare-output",
        action="store_true",
        help="recursively create this case's configured output base before launch",
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()
    paths = config_paths(args)
    if not paths:
        raise FileNotFoundError("No study configurations were selected")
    if (
        args.write_completion
        or args.clear_completion
        or args.prepare_output
    ) and len(paths) != 1:
        raise ValueError("Per-case operations require exactly one configuration")
    if args.write_completion and args.clear_completion:
        raise ValueError("Cannot write and clear a completion in the same invocation")
    runtime = args.runtime or args.write_completion
    for path in paths:
        for message in validate_one(path, runtime):
            print(f"[{path.stem}] {message}")
    if args.write_completion:
        sentinel = write_case_completion(paths[0])
        print(f"Published completion sentinel: {sentinel}")
    if args.prepare_output:
        output_base = prepare_output_base(paths[0])
        print(f"Prepared output base: {output_base}")
    if args.clear_completion:
        sentinel = clear_case_completion(paths[0])
        print(f"Cleared prior completion sentinel if present: {sentinel}")
    print(f"Validated {len(paths)} configuration(s); runtime={runtime}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
