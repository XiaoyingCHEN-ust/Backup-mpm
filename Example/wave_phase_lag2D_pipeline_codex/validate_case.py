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
import prepare_local_smoke as local_smoke
import qa_contract as qa
import hdf5_checkpoint_crosscheck as hdf5_crosscheck
from phase_controls import (
    database_paths,
    iter_frames,
    read_header,
    read_points,
    validate_database,
)


CASE_ROOT = Path(__file__).resolve().parent
COMPLETION_SCHEMA = qa.COMPLETION_SCHEMA
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


def mesh_node_bounds(path: Path) -> tuple[np.ndarray, np.ndarray]:
    """Read the node block of an MPM mesh and return its finite 2-D bounds."""

    node_count: int | None = None
    coordinates: list[list[float]] = []
    with path.open(encoding="utf-8") as stream:
        for line in stream:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            values = stripped.split()
            if node_count is None:
                if len(values) != 2:
                    raise ValueError(f"Invalid mesh count row in {path}: {stripped}")
                node_count = int(values[0])
                if node_count <= 0:
                    raise ValueError(f"Mesh contains no nodes: {path}")
                continue
            if len(coordinates) >= node_count:
                break
            if len(values) < 2:
                raise ValueError(f"Invalid mesh node row in {path}: {stripped}")
            coordinates.append([float(values[0]), float(values[1])])
    if node_count is None or len(coordinates) != node_count:
        raise ValueError(f"Mesh node block is truncated: {path}")
    points = np.asarray(coordinates, dtype=float)
    if not np.all(np.isfinite(points)):
        raise ValueError(f"Mesh node block contains non-finite values: {path}")
    return np.min(points, axis=0), np.max(points, axis=0)


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


def expected_particle_steps(config: dict[str, Any]) -> list[int]:
    """Return the exact configured particle-output grid through ``nsteps``."""

    nsteps = int(config["analysis"]["nsteps"])
    output_steps = int(config["post_processing"]["output_steps"])
    if nsteps <= 0 or output_steps <= 0:
        raise ValueError("Configured nsteps and output_steps must be positive")
    steps = list(range(output_steps, nsteps + 1, output_steps))
    if not steps or steps[-1] != nsteps:
        steps.append(nsteps)
    return steps


def expected_particle_vtp_paths(config: dict[str, Any]) -> list[Path]:
    directory = result_directory(config)
    digits = len(str(int(config["analysis"]["nsteps"])))
    return [
        directory / f"particle{step:0{digits}d}.vtp"
        for step in expected_particle_steps(config)
    ]


def exact_particle_vtp_paths(config: dict[str, Any]) -> list[Path]:
    """Reject missing, stale, or extra particle frames in a result directory."""

    expected = expected_particle_vtp_paths(config)
    actual = sorted(result_directory(config).glob("particle*.vtp"))
    if [path.name for path in actual] != [path.name for path in expected]:
        raise ValueError(
            "Particle VTP grid is incomplete or unexpected: "
            f"expected={[path.name for path in expected]}, "
            f"found={[path.name for path in actual]}"
        )
    return expected


def expected_pipeline_history_path(config: dict[str, Any]) -> Path:
    digits = len(str(int(config["analysis"]["nsteps"])))
    return result_directory(config) / f"pipeline-history{0:0{digits}d}.csv"


def exact_pipeline_history_path(config: dict[str, Any]) -> Path:
    """Require the solver's one exact, configured pipeline-history CSV."""

    expected = expected_pipeline_history_path(config)
    actual = sorted(result_directory(config).glob("pipeline-history*.csv"))
    if [path.name for path in actual] != [expected.name]:
        raise ValueError(
            "Pipeline-history CSV inventory is incomplete or unexpected: "
            f"expected={[expected.name]}, found={[path.name for path in actual]}"
        )
    return expected


def validate_final_vtp(path: Path, config: dict[str, Any] | None = None) -> None:
    points, arrays = read_vtp(checked_case_path(path))
    if len(points) == 0 or not np.all(np.isfinite(points)):
        raise ValueError(f"Final VTP has no finite particle coordinates: {path}")
    if config is not None:
        if "mesh" not in config or "particles" not in config:
            return
        mesh_path = resolve_case_path(config["mesh"]["mesh"])
        particle_path = resolve_case_path(
            config["particles"][0]["generator"]["location"]
        )
        # Direct completion-unit tests may intentionally omit the full runtime
        # geometry. The CLI always calls validate_one first, so real runs reach
        # this branch with both inputs present.
        if not mesh_path.is_file() or not particle_path.is_file():
            return
        lower, upper = mesh_node_bounds(mesh_path)
        tolerance = 1.0e-9
        if np.any(points[:, :2] < lower - tolerance) or np.any(
            points[:, :2] > upper + tolerance
        ):
            raise ValueError(f"Final VTP contains particles outside the mesh: {path}")
        expected_particles = first_count(particle_path)
        if len(points) != expected_particles:
            raise ValueError(
                f"Final VTP particle count {len(points)} differs from input "
                f"count {expected_particles}: {path}"
            )
        essential = {
            "ids",
            "porosities",
            "volumes",
            "displacements",
            "velocities",
            "stresses",
        }
        missing = essential - set(arrays)
        if missing:
            raise ValueError(f"Final VTP lacks essential fields {sorted(missing)}: {path}")


def validate_final_hdf5(path: Path) -> None:
    path = checked_case_path(path)
    if not path.is_file() or path.stat().st_size <= len(HDF5_SIGNATURE):
        raise ValueError(f"Final HDF5 is missing or too short: {path}")
    with path.open("rb") as stream:
        if stream.read(len(HDF5_SIGNATURE)) != HDF5_SIGNATURE:
            raise ValueError(f"Final HDF5 has an invalid signature: {path}")


def hdf5_auditor_path() -> Path:
    """Resolve the required structural HDF5 auditor without a fallback pass."""

    configured = os.environ.get("MPM_HDF5_AUDITOR")
    path = (
        Path(configured).expanduser()
        if configured
        else hdf5_crosscheck.DEFAULT_AUDITOR
    ).resolve()
    if not path.is_file():
        source = "MPM_HDF5_AUDITOR" if configured else "default build-pipeline path"
        raise FileNotFoundError(f"HDF5 checkpoint auditor from {source} is missing: {path}")
    if not os.access(path, os.X_OK):
        raise PermissionError(f"HDF5 checkpoint auditor is not executable: {path}")
    return path


def hdf5_vtp_crosscheck_audit(hdf5: Path, vtp: Path) -> dict[str, Any]:
    """Structurally inspect HDF5 and compare every core field to same-step VTP."""

    hdf5 = checked_case_path(hdf5)
    vtp = checked_case_path(vtp)
    record = hdf5_crosscheck.crosscheck(hdf5, vtp, hdf5_auditor_path())
    return qa.validate_hdf5_vtp_audit(
        record, expected_hdf5=str(hdf5), expected_vtp=str(vtp)
    )


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
        stability_qa = validate_equilibrium_vtp(
            checkpoint, config, mode=resume_stability_mode(config)
        )
        dependencies["resume_equilibrium"] = {
            "checkpoint_hdf5": artifact_audit(checkpoint),
            "qa_vtp": artifact_audit(equilibrium_vtp),
            "hdf5_vtp_crosscheck": hdf5_vtp_crosscheck_audit(
                checkpoint, equilibrium_vtp
            ),
            "stability_qa": stability_qa,
        }
    pressure = read_pressure_database_audit(config)
    if pressure is not None:
        dependencies["read_pressure_database"] = pressure
    return dependencies


def write_case_completion(
    config_path: Path, *, validation_profile: str = "registered-study-v1"
) -> Path:
    """Validate final artifacts and atomically publish a restart sentinel."""

    config_path = checked_case_path(config_path)
    config = json.loads(config_path.read_text(encoding="utf-8"))
    final_vtp, final_hdf5 = final_artifact_paths(config)
    particle_vtps = exact_particle_vtp_paths(config)
    history_csv = exact_pipeline_history_path(config)
    validate_final_vtp(final_vtp, config)
    artifacts: dict[str, Any] = {
        "particle_vtp_grid": {
            "steps": expected_particle_steps(config),
            "files": [artifact_audit(path) for path in particle_vtps],
        },
        "pipeline_history_csv": artifact_audit(history_csv),
    }
    if final_hdf5 is not None:
        validate_final_hdf5(final_hdf5)
        artifacts["final_hdf5"] = artifact_audit(final_hdf5)
        artifacts["hdf5_vtp_crosscheck"] = hdf5_vtp_crosscheck_audit(
            final_hdf5, final_vtp
        )
    stability_qa = None
    if config["analysis"].get("stability_gate"):
        if final_hdf5 is None:
            raise ValueError("A stability-gated stage must write an HDF5 checkpoint")
        stability_qa = validate_equilibrium_vtp(
            final_hdf5, config, mode=own_stability_mode(config)
        )
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
            "validation_profile": validation_profile,
        },
        "artifacts": artifacts,
        "runtime_dependencies": runtime_dependencies,
    }
    if stability_qa is not None:
        payload["stability_qa"] = stability_qa
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
        if not np.allclose(stresses[:, 2], expected_xx, rtol=0.0, atol=1.0e-6):
            raise ValueError("Initial out-of-plane effective stress has stale K0")
        if not np.allclose(stresses[:, 3:], 0.0, rtol=0.0, atol=1.0e-12):
            raise ValueError("Initial shear stresses must be zero")
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


def own_stability_mode(config: dict[str, Any]) -> str:
    """Return the registered QA mode for this stability-gated stage."""

    analysis = config["analysis"]
    material_type = str(config["materials"][0]["type"])
    uuid = str(analysis["uuid"])
    resume_enabled = bool(analysis.get("resume", {}).get("resume"))
    handoff = bool(analysis.get("handoff_relaxation"))
    if (
        not resume_enabled
        and not handoff
        and material_type == "LinearElastic2D"
        and uuid.endswith("_EQ")
    ):
        mode = qa.LINEAR_EQUILIBRIUM_MODE
    elif (
        resume_enabled
        and handoff
        and material_type == "MohrCoulomb2D"
        and uuid.endswith("_MC_RELAX")
    ):
        mode = qa.MC_HANDOFF_MODE
    else:
        raise ValueError(f"{uuid}: stability-gated stage has no registered QA mode")
    if analysis.get("stability_qa_contract") != qa.stability_contract(mode):
        raise ValueError(f"{uuid}: own-stage stability QA contract is stale")
    return mode


def resume_stability_mode(config: dict[str, Any]) -> str:
    """Return the registered QA mode required of a resume checkpoint."""

    analysis = config["analysis"]
    resume = analysis.get("resume", {})
    if not bool(resume.get("resume")):
        raise ValueError(f"{analysis['uuid']}: case has no enabled resume checkpoint")
    source_uuid = str(resume.get("uuid", ""))
    if source_uuid.endswith("_MC_RELAX"):
        mode = qa.MC_HANDOFF_MODE
    elif source_uuid.endswith("_EQ"):
        mode = qa.LINEAR_EQUILIBRIUM_MODE
    else:
        raise ValueError(
            f"{analysis['uuid']}: resume UUID has no registered stability QA mode"
        )
    if analysis.get("resume_stability_qa_contract") != qa.stability_contract(mode):
        raise ValueError(f"{analysis['uuid']}: resume stability QA contract is stale")
    return mode


def _effective_unit_weight(config: dict[str, Any]) -> float:
    soil, fluid = config["materials"]
    porosity = float(soil["porosity"])
    liquid_density = float(fluid["density"])
    gas_density = (
        float(fluid["gas_molar_mass"])
        * float(soil["p_ref"])
        / float(fluid["gas_constant"])
        / (273.15 + float(soil["initial_temperature"]))
    )
    mixture_density = (1.0 - porosity) * float(soil["density"]) + porosity * (
        float(fluid["liquid_saturation"]) * liquid_density
        + float(fluid["gas_saturation"]) * gas_density
    )
    unit_weight = (mixture_density - liquid_density) * study.GRAVITY
    if not math.isfinite(unit_weight) or unit_weight <= 0.0:
        raise ValueError("Equilibrium effective unit weight must be positive")
    return unit_weight


def _linear_initial_stress_path(config: dict[str, Any]) -> Path:
    configured = config["mesh"].get("particles_stresses")
    if configured:
        return resolve_case_path(configured)
    saturation = float(config["materials"][1]["liquid_saturation"])
    if math.isclose(
        saturation, study.LOW_LAG_SATURATION, rel_tol=0.0, abs_tol=1.0e-12
    ):
        state = "LS"
    elif math.isclose(
        saturation, study.HIGH_LAG_SATURATION, rel_tol=0.0, abs_tol=1.0e-12
    ):
        state = "HS"
    else:
        raise ValueError("Resume checkpoint saturation has no registered stress state")
    particles_path = resolve_case_path(
        config["particles"][0]["generator"]["location"]
    )
    return particles_path.parent / f"initial_effective_stresses_{state}.txt"


def _particle_ids(arrays: dict[str, np.ndarray], count: int, vtp: Path) -> np.ndarray:
    if "ids" not in arrays:
        raise ValueError(f"Equilibrium QA field ids is absent in {vtp}")
    raw = np.asarray(arrays["ids"]).reshape(-1)
    if raw.size != count or not np.all(np.isfinite(raw)):
        raise ValueError("Equilibrium VTP particle IDs are incomplete or non-finite")
    identifiers = raw.astype(np.int64)
    if not np.array_equal(raw, identifiers):
        raise ValueError("Equilibrium VTP particle IDs are not integral")
    if not np.array_equal(np.sort(identifiers), np.arange(count, dtype=np.int64)):
        raise ValueError("Equilibrium VTP particle IDs are missing or duplicated")
    return identifiers


def _particle_set(entity_sets: dict[str, Any], set_id: int) -> np.ndarray:
    matches = [
        item for item in entity_sets.get("particle_sets", []) if int(item["id"]) == set_id
    ]
    if len(matches) != 1:
        raise ValueError(f"Entity sets must contain exactly one particle set {set_id}")
    identifiers = np.asarray(matches[0].get("set", []), dtype=np.int64).reshape(-1)
    if identifiers.size == 0 or np.unique(identifiers).size != identifiers.size:
        raise ValueError(f"Particle set {set_id} is empty or contains duplicates")
    return identifiers


def _normal_stress_qa(
    config: dict[str, Any], arrays: dict[str, np.ndarray], identifiers: np.ndarray
) -> dict[str, Any]:
    mesh_path = resolve_case_path(config["mesh"]["mesh"])
    particle_path = resolve_case_path(
        config["particles"][0]["generator"]["location"]
    )
    stress_path = _linear_initial_stress_path(config)
    entity_sets_path = resolve_case_path(config["mesh"]["entity_sets"])
    for description, path in (
        ("mesh", mesh_path),
        ("initial particles", particle_path),
        ("initial stresses", stress_path),
        ("entity sets", entity_sets_path),
    ):
        if not path.is_file():
            raise FileNotFoundError(f"Equilibrium QA {description} is missing: {path}")

    coordinates = read_ascii_table(particle_path, 2)
    initial_stresses = read_ascii_table(stress_path, 6)
    count = len(identifiers)
    if coordinates.shape != (count, 2) or initial_stresses.shape != (count, 6):
        raise ValueError("Equilibrium QA input count differs from VTP particle count")
    stresses = np.asarray(arrays.get("stresses"))
    if stresses.shape != (count, 6):
        raise ValueError(f"Equilibrium QA field stresses is absent or malformed")

    entity_sets = json.loads(entity_sets_path.read_text(encoding="utf-8"))
    near_surface_ids = _particle_set(entity_sets, 0)
    top_row_ids = _particle_set(entity_sets, 4)
    if (
        np.any(near_surface_ids < 0)
        or np.any(near_surface_ids >= count)
        or np.any(top_row_ids < 0)
        or np.any(top_row_ids >= count)
    ):
        raise ValueError("Equilibrium QA surface particle set contains invalid IDs")
    elevations = np.unique(coordinates[:, 1])
    if elevations.size < 2:
        raise ValueError("Equilibrium QA requires at least two initial particle rows")
    highest_rows = np.sort(elevations[-2:])
    coordinate_tolerance = max(1.0e-12, float(config["mesh"]["cellsize_min"]) * 1.0e-9)
    expected_near_surface = np.flatnonzero(
        np.isclose(coordinates[:, 1, None], highest_rows[None, :], rtol=0.0,
                   atol=coordinate_tolerance).any(axis=1)
    )
    expected_top = np.flatnonzero(
        np.isclose(
            coordinates[:, 1], highest_rows[-1], rtol=0.0, atol=coordinate_tolerance
        )
    )
    if not np.array_equal(np.sort(near_surface_ids), expected_near_surface):
        raise ValueError("Particle set 0 is not exactly the highest two initial rows")
    if not np.array_equal(np.sort(top_row_ids), expected_top):
        raise ValueError("Particle set 4 is not exactly the highest initial row")

    cell_size = float(config["mesh"]["cellsize_min"])
    if not math.isfinite(cell_size) or cell_size <= 0.0:
        raise ValueError("Equilibrium QA background cell size is invalid")
    lower, upper = mesh_node_bounds(mesh_path)
    pipeline = config["analysis"]["rigid_pipeline"]
    pipe_x = float(pipeline["center"][0])
    pipe_radius = float(pipeline["outer_radius"])
    side_buffer = qa.NORMAL_STRESS_QUIET_SIDE_BUFFER_CELLS * cell_size
    pipe_buffer = qa.NORMAL_STRESS_QUIET_PIPE_BUFFER_CELLS * cell_size
    x = coordinates[:, 0]
    quiet_initial_ids = (
        np.isin(np.arange(count), near_surface_ids)
        & (x >= lower[0] + side_buffer - coordinate_tolerance)
        & (x <= upper[0] - side_buffer + coordinate_tolerance)
        & (np.abs(x - pipe_x) >= pipe_radius + pipe_buffer - coordinate_tolerance)
    )

    unit_weight = _effective_unit_weight(config)
    seabed = float(config["materials"][1]["sea_level"]) - float(
        config["materials"][1]["depth_left"]
    )
    expected_initial_yy = -unit_weight * np.maximum(
        seabed - coordinates[:, 1], 0.0
    )
    if not np.allclose(
        initial_stresses[:, 1], expected_initial_yy, rtol=0.0, atol=1.0e-6
    ):
        raise ValueError("Equilibrium QA initial vertical stresses are not geostatic")

    final_yy_by_initial_id = np.empty(count, dtype=float)
    final_yy_by_initial_id[identifiers] = stresses[:, 1]
    tolerance = (
        qa.NORMAL_STRESS_MEDIAN_DRIFT_CELL_WEIGHTS * unit_weight * cell_size
    )
    rows: list[dict[str, Any]] = []
    for elevation in highest_rows:
        row_ids = quiet_initial_ids & np.isclose(
            coordinates[:, 1], elevation, rtol=0.0, atol=coordinate_tolerance
        )
        row_count = int(np.count_nonzero(row_ids))
        if row_count < qa.NORMAL_STRESS_MINIMUM_PARTICLES_PER_ROW:
            raise ValueError(
                f"Equilibrium normal-stress quiet row {elevation:g} m has only "
                f"{row_count} particles"
            )
        drift = final_yy_by_initial_id[row_ids] - initial_stresses[row_ids, 1]
        median_drift = float(np.median(drift))
        tensile_fraction = float(np.mean(final_yy_by_initial_id[row_ids] > 0.0))
        if abs(median_drift) > tolerance:
            raise ValueError(
                f"Equilibrium row y={elevation:g} m median vertical-stress drift "
                f"{median_drift:.3e} Pa exceeds gamma_eff*h={tolerance:.3e} Pa"
            )
        if tensile_fraction > qa.NORMAL_STRESS_MAXIMUM_TENSILE_FRACTION:
            raise ValueError(
                f"Equilibrium row y={elevation:g} m tensile fraction "
                f"{tensile_fraction:.3%} exceeds 5%"
            )
        rows.append(
            {
                "initial_y_m": float(elevation),
                "quiet_particle_count": row_count,
                "initial_sigma_yy_median_pa": float(
                    np.median(initial_stresses[row_ids, 1])
                ),
                "final_sigma_yy_median_pa": float(
                    np.median(final_yy_by_initial_id[row_ids])
                ),
                "median_drift_pa": median_drift,
                "tensile_fraction": tensile_fraction,
            }
        )
    return {
        "inputs": {
            "mesh": artifact_audit(mesh_path),
            "particles": artifact_audit(particle_path),
            "initial_stresses": artifact_audit(stress_path),
            "entity_sets": artifact_audit(entity_sets_path),
        },
        "particle_id_basis": "initial-particle-file",
        "background_cell_size_m": cell_size,
        "effective_unit_weight_n_m3": unit_weight,
        "median_drift_tolerance_pa": tolerance,
        "selection": {
            "initial_surface_pset_id": 0,
            "initial_top_row_pset_id": 4,
            "side_buffer_m": side_buffer,
            "pipe_buffer_m": pipe_buffer,
            "domain_x_m": [float(lower[0]), float(upper[0])],
            "pipe_center_x_m": pipe_x,
            "pipe_outer_radius_m": pipe_radius,
        },
        "rows": rows,
    }


def _mc_feasibility_qa(
    config: dict[str, Any], arrays: dict[str, np.ndarray], count: int
) -> dict[str, Any]:
    required_shapes = {
        "stresses": (count, 6),
        "phi": (count,),
        "cohesion": (count,),
    }
    values: dict[str, np.ndarray] = {}
    for name, shape in required_shapes.items():
        value = np.asarray(arrays.get(name)).reshape(shape) if name in arrays else None
        if value is None or value.shape != shape or not np.all(np.isfinite(value)):
            raise ValueError(f"Mohr-Coulomb QA field {name} is absent or malformed")
        values[name] = value.astype(float, copy=False)
    stresses = values["stresses"]
    phi = values["phi"]
    cohesion = values["cohesion"]
    mean_stress = np.mean(stresses[:, :3], axis=1)
    deviator = stresses.copy()
    deviator[:, :3] -= mean_stress[:, None]
    j2 = (
        (stresses[:, 0] - stresses[:, 1]) ** 2
        + (stresses[:, 1] - stresses[:, 2]) ** 2
        + (stresses[:, 0] - stresses[:, 2]) ** 2
    ) / 6.0 + np.sum(stresses[:, 3:] ** 2, axis=1)
    j3 = (
        deviator[:, 0] * deviator[:, 1] * deviator[:, 2]
        - deviator[:, 2] * stresses[:, 3] ** 2
        + 2.0 * stresses[:, 3] * stresses[:, 4] * stresses[:, 5]
        - deviator[:, 0] * stresses[:, 4] ** 2
        - deviator[:, 1] * stresses[:, 5] ** 2
    )
    theta_argument = np.zeros(count, dtype=float)
    nonzero_j2 = np.abs(j2) > 0.0
    theta_argument[nonzero_j2] = (
        1.5 * math.sqrt(3.0) * j3[nonzero_j2] / np.power(j2[nonzero_j2], 1.5)
    )
    theta = np.arccos(np.clip(theta_argument, -1.0, 1.0)) / 3.0
    rho = np.sqrt(2.0 * j2)
    epsilon = np.sum(stresses[:, :3], axis=1) / math.sqrt(3.0)
    configured_tension = float(config["materials"][0]["tension_cutoff"])
    tan_phi = np.tan(phi)
    strength_limit = np.full(count, math.inf, dtype=float)
    nonzero_tan = np.abs(tan_phi) > np.finfo(float).tiny
    strength_limit[nonzero_tan] = cohesion[nonzero_tan] / tan_phi[nonzero_tan]
    tension_cutoff = np.minimum(configured_tension, strength_limit)
    tension = (
        math.sqrt(2.0 / 3.0) * np.cos(theta) * rho
        + epsilon / math.sqrt(3.0)
        - tension_cutoff
    )
    shear = (
        math.sqrt(1.5)
        * rho
        * (
            np.sin(theta + math.pi / 3.0) / (math.sqrt(3.0) * np.cos(phi))
            + np.cos(theta + math.pi / 3.0) * np.tan(phi) / 3.0
        )
        + epsilon / math.sqrt(3.0) * np.tan(phi)
        - cohesion
    )
    if not np.all(np.isfinite(tension)) or not np.all(np.isfinite(shear)):
        raise ValueError("Mohr-Coulomb QA yield residual is non-finite")
    positive = np.maximum.reduce((np.zeros(count), tension, shear))
    maximum_positive = float(np.max(positive))
    violations = int(
        np.count_nonzero(positive > qa.MC_MAXIMUM_POSITIVE_YIELD_RESIDUAL_PA)
    )
    if violations:
        raise ValueError(
            f"Mohr-Coulomb positive yield residual {maximum_positive:.3e} Pa "
            f"exceeds {qa.MC_MAXIMUM_POSITIVE_YIELD_RESIDUAL_PA:.0e} Pa for "
            f"{violations} particles"
        )
    return {
        "particle_count": count,
        "maximum_tension_residual_pa": float(np.max(tension)),
        "maximum_shear_residual_pa": float(np.max(shear)),
        "maximum_positive_residual_pa": maximum_positive,
        "violating_particle_count": violations,
    }


def validate_equilibrium_vtp(
    checkpoint: Path, config: dict[str, Any], *, mode: str
) -> dict[str, Any]:
    """Validate a static checkpoint and return its versioned QA audit."""

    if mode not in qa.STABILITY_MODES:
        raise ValueError(f"Unknown equilibrium QA mode: {mode}")
    vtp = checkpoint.with_name(checkpoint.name.replace("particles", "particle")).with_suffix(
        ".vtp"
    )
    if not vtp.is_file():
        raise FileNotFoundError(f"Equilibrium QA VTP is missing: {vtp}")
    points, arrays = read_vtp(vtp)
    count = len(points)
    for required in ("velocities", "displacements", "porosities"):
        if required not in arrays:
            raise ValueError(f"Equilibrium QA field {required} is absent in {vtp}")
    velocities = np.asarray(arrays["velocities"])
    displacements = np.asarray(arrays["displacements"])
    porosities = np.asarray(arrays["porosities"]).reshape(-1)
    if (
        velocities.shape != (count, 3)
        or displacements.shape != (count, 3)
        or porosities.shape != (count,)
    ):
        raise ValueError("Equilibrium QA kinematic field shape is malformed")
    identifiers = _particle_ids(arrays, count, vtp)
    maximum_velocity = float(np.max(np.linalg.norm(velocities[:, :2], axis=1)))
    maximum_displacement = float(
        np.max(np.linalg.norm(displacements[:, :2], axis=1))
    )
    porosity_min = float(np.min(porosities))
    porosity_max = float(np.max(porosities))
    if maximum_velocity > qa.MAXIMUM_VELOCITY_M_S:
        raise ValueError(
            f"Equilibrium max velocity {maximum_velocity:.3e} m/s exceeds 1e-3 m/s"
        )
    if maximum_displacement > qa.MAXIMUM_DISPLACEMENT_M:
        raise ValueError(
            f"Equilibrium max displacement {maximum_displacement:.3e} m exceeds "
            f"{qa.MAXIMUM_DISPLACEMENT_M:g} m"
        )
    if porosity_min <= qa.POROSITY_MIN_EXCLUSIVE or porosity_max >= qa.POROSITY_MAX_EXCLUSIVE:
        raise ValueError("Equilibrium porosity lies outside (0, 1)")
    observed = {
        "maximum_velocity_m_s": maximum_velocity,
        "maximum_displacement_m": maximum_displacement,
        "porosity_min": porosity_min,
        "porosity_max": porosity_max,
        "particle_count": count,
    }
    record: dict[str, Any] = {
        "schema": qa.STABILITY_QA_SCHEMA,
        "mode": mode,
        "limits": qa.stability_contract(mode)["limits"],
        "observed": observed,
    }
    if mode == qa.LINEAR_EQUILIBRIUM_MODE:
        detail = _normal_stress_qa(config, arrays, identifiers)
        record["normal_stress"] = detail
        row_text = ", ".join(
            f"y={row['initial_y_m']:.6g}:N={row['quiet_particle_count']},"
            f"dsyy50={row['median_drift_pa']:.3e} Pa,"
            f"ft={row['tensile_fraction']:.3%}"
            for row in detail["rows"]
        )
        mode_text = f"normal stress [{row_text}]"
    else:
        detail = _mc_feasibility_qa(config, arrays, count)
        record["mc_feasibility"] = detail
        mode_text = (
            "MC max positive yield residual="
            f"{detail['maximum_positive_residual_pa']:.3e} Pa"
        )
    record["summary"] = (
        f"equilibrium QA ({mode}): max|v|={maximum_velocity:.3e} m/s, "
        f"max|u|={maximum_displacement:.3e} m, "
        f"n=[{porosity_min:.5f}, {porosity_max:.5f}], {mode_text}"
    )
    return qa.validate_stability_qa(record, expected_mode=mode)


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
        stability_qa = validate_equilibrium_vtp(
            checkpoint, config, mode=resume_stability_mode(config)
        )
        messages.append(stability_qa["summary"])
    messages.extend(validate_pressure_dependency(config, runtime))
    return messages


def validate_local_smoke_config(config: dict[str, Any]) -> None:
    """Accept only an exact hash-bound config generated by the smoke profile."""

    marker = config.get("local_smoke")
    if not isinstance(marker, dict):
        raise ValueError("Local smoke marker is absent or malformed")
    if set(marker) != {
        "schema",
        "label",
        "role",
        "baseline_sha256",
        "result_eligibility",
    }:
        raise ValueError("Local smoke marker fields are missing or unexpected")
    if marker.get("schema") != local_smoke.SMOKE_SCHEMA:
        raise ValueError("Local smoke schema is not registered")
    label = local_smoke.validate_label(str(marker.get("label", "")))
    role = str(marker.get("role", ""))
    if role not in local_smoke.ROLE_FILENAMES:
        raise ValueError(f"Local smoke role is not allowed: {role}")
    if marker.get("result_eligibility") != (
        "local-smoke-only; excluded from manuscript evidence"
    ):
        raise ValueError("Local smoke evidence exclusion is absent")

    payload = dict(config)
    payload.pop("local_smoke")
    observed_hash = local_smoke.payload_sha256(payload)
    expected = local_smoke.canonical_config(label, role)
    expected_hash = expected["local_smoke"]["baseline_sha256"]
    if marker.get("baseline_sha256") != observed_hash:
        raise ValueError("Local smoke payload differs from its recorded baseline hash")
    if observed_hash != expected_hash or config != expected:
        raise ValueError(
            f"Local smoke {label}/{role} differs from the exact canonical profile"
        )


def validate_one(
    config_path: Path, runtime: bool, *, use_local_smoke_profile: bool = False
) -> list[str]:
    config_path = config_path.resolve()
    config = json.loads(config_path.read_text(encoding="utf-8"))
    if use_local_smoke_profile:
        validate_local_smoke_config(config)
    else:
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
    parser.add_argument(
        "--local-smoke",
        action="store_true",
        help=(
            "validate only an exact hash-bound pipeline-local-smoke-v1 config; "
            "never relax the registered study validator"
        ),
    )
    parser.add_argument(
        "--print-result-directory",
        action="store_true",
        help="structurally validate one config and print only its resolved result directory",
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
        or args.print_result_directory
    ) and len(paths) != 1:
        raise ValueError("Per-case operations require exactly one configuration")
    if args.write_completion and args.clear_completion:
        raise ValueError("Cannot write and clear a completion in the same invocation")
    if args.print_result_directory:
        config = json.loads(paths[0].read_text(encoding="utf-8"))
        if args.local_smoke:
            validate_local_smoke_config(config)
        else:
            study.validate_config(config)
        print(result_directory(config))
        return 0
    runtime = args.runtime or args.write_completion
    for path in paths:
        for message in validate_one(
            path, runtime, use_local_smoke_profile=args.local_smoke
        ):
            print(f"[{path.stem}] {message}")
    if args.write_completion:
        profile = (
            local_smoke.SMOKE_SCHEMA if args.local_smoke else "registered-study-v1"
        )
        sentinel = write_case_completion(paths[0], validation_profile=profile)
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
