#!/usr/bin/env python3
"""Audit and plot an exploratory lagged/phase-erased replay pair."""

from __future__ import annotations

import argparse
import copy
import csv
import hashlib
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np

import phase_controls as phase_control_module
from analyze_study import (
    THRESHOLD_RATIO_TOLERANCE,
    expected_particle_steps,
    particle_step,
    read_initial_coordinates,
    read_pipeline_history,
    read_vtp,
    upward_excess_seepage_force_index,
    vertical_stress_loss_state,
    validate_pipeline_history_grid,
)
from plot_liquid_pressure_phase_lag import draw_field
from phase_controls import (
    particle_ids_sha256,
    read_counted_particle_ids,
    validate_phase_alignment_qa,
)


SCHEMA = "pipeline-phase-lag-exploratory-replay-analysis-v2"
EXPECTED_TRANSFORM_PARAMETERS = {
    "period_s": 1.3,
    "fit_start_s": 2.6,
    "fit_end_s": 3.9,
    "ramp_time_s": 1.3,
    "surface_band_m": 0.015,
    "minimum_reference_amplitude_pa": 1.0,
}
TRANSFORM_INVARIANCE_TOLERANCE_PA = 1.0e-10
REQUIRED_ARRAYS = (
    "ids",
    "volumes",
    "velocities",
    "displacements",
    "porosities",
    "stresses",
    "initial_vertical_effective_stresses",
    "vertical_effective_stress_remaining_ratios",
    "PIC_pore_pressure_excess",
    "liquid_seepage_forces",
    "liquid_densities",
    "gamma_sub",
)


@dataclass(frozen=True)
class Replay:
    role: str
    config_name: str


@dataclass(frozen=True)
class ReplayParameters:
    intrinsic_permeability_m2: float
    liquid_saturation: float
    gas_saturation: float
    pressure_smoothing: bool
    pressure_smoothing_in_loop: bool
    pipeline_fixed: bool

    def as_dict(self) -> dict[str, float | bool]:
        return {
            "intrinsic_permeability_m2": self.intrinsic_permeability_m2,
            "liquid_saturation": self.liquid_saturation,
            "gas_saturation": self.gas_saturation,
            "pressure_smoothing": self.pressure_smoothing,
            "pressure_smoothing_in_loop": self.pressure_smoothing_in_loop,
            "pipeline_fixed": self.pipeline_fixed,
        }


REPLAYS = (Replay("lagged", "04_RL.json"), Replay("phase_erased", "05_RE.json"))


def replay_parameters(config: dict[str, Any], role: str) -> ReplayParameters:
    materials = config.get("materials")
    if not isinstance(materials, list) or len(materials) != 2:
        raise ValueError(f"{role}: expected exactly two material definitions")
    soil, fluid = materials
    try:
        permeability = float(soil["intrinsic_permeability"])
        saturation = float(fluid["liquid_saturation"])
        gas_saturation = float(fluid["gas_saturation"])
        analysis = config["analysis"]
        pressure_smoothing = analysis["pressure_smoothing"]
        pressure_smoothing_in_loop = analysis["pressure_smoothing_in_loop"]
        pipeline_fixed = analysis["rigid_pipeline"]["fixed"]
    except (KeyError, TypeError, ValueError) as error:
        raise ValueError(f"{role}: incomplete hydraulic replay parameters") from error
    if not math.isfinite(permeability) or permeability <= 0.0:
        raise ValueError(f"{role}: intrinsic permeability must be finite and positive")
    if not math.isfinite(saturation) or not 0.0 < saturation < 1.0:
        raise ValueError(f"{role}: liquid saturation must lie strictly inside (0,1)")
    if not math.isfinite(gas_saturation) or not math.isclose(
        saturation + gas_saturation, 1.0, rel_tol=0.0, abs_tol=1.0e-12
    ):
        raise ValueError(f"{role}: liquid and gas saturations must sum to one")
    # JSON booleans are required here.  Numeric zero or a missing key must not
    # silently pass as "false" because smoothing changes the evidence pathway.
    if pressure_smoothing is not False or pressure_smoothing_in_loop is not False:
        raise ValueError(f"{role}: pressure smoothing must be explicitly disabled")
    if pipeline_fixed is not True:
        raise ValueError(f"{role}: pipeline must be explicitly fixed")
    return ReplayParameters(
        intrinsic_permeability_m2=permeability,
        liquid_saturation=saturation,
        gas_saturation=gas_saturation,
        pressure_smoothing=pressure_smoothing,
        pressure_smoothing_in_loop=pressure_smoothing_in_loop,
        pipeline_fixed=pipeline_fixed,
    )


def validate_replay_pair_configs(
    lagged_config: dict[str, Any], erased_config: dict[str, Any]
) -> ReplayParameters:
    lagged_parameters = replay_parameters(lagged_config, "lagged")
    erased_parameters = replay_parameters(erased_config, "phase_erased")
    if lagged_parameters != erased_parameters:
        raise ValueError(
            "Replay configs disagree on permeability, saturation, smoothing, "
            "or pipeline constraint"
        )

    lagged_physics = copy.deepcopy(lagged_config)
    erased_physics = copy.deepcopy(erased_config)
    for config in (lagged_physics, erased_physics):
        config["title"] = "paired-replay"
        config["analysis"]["uuid"] = "paired-replay"
        pressure = config["analysis"]["prescribed_phase_pressures"]
        pressure["path"] = "paired-replay"
        pressure["file_prefix"] = "paired-replay"
    if lagged_physics != erased_physics:
        raise ValueError(
            "Replay configs differ in fields other than pressure source identity"
        )
    return lagged_parameters


def parameter_label(parameters: ReplayParameters) -> str:
    return (
        rf"$k={parameters.intrinsic_permeability_m2:.3e}\ \mathrm{{m}}^2$, "
        rf"$S_w={parameters.liquid_saturation:.3f}$, no smoothing; fixed pipeline"
    )


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def configured_pressure_artifact_paths(
    root: Path, label: str, *, role: str
) -> dict[str, Path]:
    """Resolve the exact pressure database selected by one replay config."""

    config_name = {"lagged": "04_RL.json", "phase_erased": "05_RE.json"}.get(
        role
    )
    if config_name is None:
        raise ValueError(f"Unsupported replay pressure role: {role}")
    config_path = root / "configs/phase_lag_exploratory" / label / config_name
    try:
        config = json.loads(config_path.read_text(encoding="utf-8"))
        pressure = config["analysis"]["prescribed_phase_pressures"]
        enabled = pressure["enable"]
        write = pressure["write"]
        directory_value = pressure["path"]
        prefix = pressure["file_prefix"]
    except (OSError, json.JSONDecodeError, KeyError, TypeError) as error:
        raise ValueError(f"{role}: replay pressure config is incomplete") from error
    if enabled is not True or write is not False:
        raise ValueError(f"{role}: replay pressure config is not read-only")
    if not isinstance(directory_value, str) or not directory_value:
        raise ValueError(f"{role}: replay pressure path is invalid")
    if not isinstance(prefix, str) or not prefix:
        raise ValueError(f"{role}: replay pressure prefix is invalid")
    directory = Path(directory_value)
    if not directory.is_absolute():
        directory = root / directory
    return {
        "points": (directory / f"{prefix}_points.txt").resolve(),
        "values": (directory / f"{prefix}_values.bin").resolve(),
    }


def _finite_metadata_number(value: Any, context: str) -> float:
    if isinstance(value, bool):
        raise ValueError(f"{context} must be finite")
    try:
        number = float(value)
    except (TypeError, ValueError) as error:
        raise ValueError(f"{context} must be finite") from error
    if not math.isfinite(number):
        raise ValueError(f"{context} must be finite")
    return number


def validate_phase_metadata(root: Path, label: str) -> dict[str, Any]:
    root = root.resolve()
    path = (
        root
        / "pressure_databases/phase_lag_exploratory"
        / label
        / "phase_erased/phase_erased_metadata.json"
    )
    metadata = json.loads(path.read_text(encoding="utf-8"))
    if metadata.get("schema") != "phase-erased-pressure-control-v1":
        raise ValueError("Phase-erased metadata schema is invalid")
    producer = metadata.get("producer")
    producer_path = Path(phase_control_module.__file__).resolve()
    if not isinstance(producer, dict) or set(producer) != {"path", "sha256"}:
        raise ValueError("Phase-erased metadata producer binding is incomplete")
    recorded_producer_path = Path(str(producer.get("path", "")))
    if (
        not recorded_producer_path.is_absolute()
        or recorded_producer_path.resolve() != producer_path
        or not producer_path.is_file()
        or producer.get("sha256") != sha256(producer_path)
    ):
        raise ValueError("Phase-erased metadata producer binding is stale")
    header = metadata.get("header", {})
    expected_header = {
        "format_version": "V2",
        "dimension": 2,
        "particle_count": 7376,
        "step_interval": 130,
        "max_step": 39000,
        "source_dt_s": 1.0e-4,
        "sample_dt_s": 0.013000000000000001,
    }
    if header != expected_header:
        raise ValueError(f"Phase-erased database header changed: {header}")
    transform = metadata.get("transform", {})
    if not isinstance(transform, dict):
        raise ValueError("Phase transform parameters are missing")
    for name, expected in EXPECTED_TRANSFORM_PARAMETERS.items():
        observed = _finite_metadata_number(
            transform.get(name), f"Phase transform {name}"
        )
        if not math.isclose(
            observed, expected, rel_tol=0.0, abs_tol=1.0e-12
        ):
            raise ValueError(f"Phase transform {name} changed")
    reference = transform.get("surface_reference", {})
    if (
        reference.get("method") != "registered_top_surface_particle_set"
        or int(reference.get("reference_particle_count", 0)) != 150
    ):
        raise ValueError("Phase transform does not use the 150-point top row")
    registered_path = (root / "top_surface_traction_particle_id.txt").resolve()
    if not registered_path.is_file() or registered_path.stat().st_size <= 0:
        raise ValueError(
            f"Registered top-surface particle file is missing or empty: {registered_path}"
        )
    source = reference.get("source")
    if not isinstance(source, dict):
        raise ValueError("Phase transform surface-reference source is missing")
    source_value = source.get("path")
    if not isinstance(source_value, str) or not source_value:
        raise ValueError("Phase transform surface-reference source path is missing")
    source_path = Path(source_value)
    if not source_path.is_absolute():
        source_path = root / source_path
    source_path = source_path.resolve()
    registered_size = registered_path.stat().st_size
    registered_sha256 = sha256(registered_path)
    if (
        source_path != registered_path
        or not source_path.is_file()
        or source_path.stat().st_size != registered_size
        or sha256(source_path) != registered_sha256
        or source.get("sha256") != registered_sha256
    ):
        raise ValueError(
            "Phase transform surface-reference source is not the registered "
            "top-surface particle file"
        )
    recorded_size = source.get("size_bytes")
    if recorded_size is not None and (
        isinstance(recorded_size, bool)
        or not isinstance(recorded_size, int)
        or recorded_size != registered_size
    ):
        raise ValueError("Phase transform surface-reference source size mismatch")
    registered_ids = read_counted_particle_ids(registered_path)
    reference_ids = reference.get("reference_particle_ids")
    if (
        registered_ids.size != 150
        or not isinstance(reference_ids, list)
        or any(isinstance(value, bool) or not isinstance(value, int) for value in reference_ids)
        or reference_ids != [int(value) for value in registered_ids]
        or reference.get("reference_particle_ids_sha256")
        != particle_ids_sha256(registered_ids)
    ):
        raise ValueError(
            "Phase transform particle IDs do not match the registered top-surface set"
        )
    configured_artifacts = {
        "source": configured_pressure_artifact_paths(
            root, label, role="lagged"
        ),
        "output": configured_pressure_artifact_paths(
            root, label, role="phase_erased"
        ),
    }
    for group in ("source", "output"):
        record = metadata.get(group, {})
        if not isinstance(record, dict):
            raise ValueError(f"Phase metadata {group} record is missing")
        for kind in ("points", "values"):
            recorded_path = record.get(kind)
            if not isinstance(recorded_path, str) or not recorded_path:
                raise ValueError(f"Phase metadata {group} {kind} path is missing")
            artifact = Path(recorded_path)
            if not artifact.is_absolute():
                artifact = root / artifact
            artifact = artifact.resolve()
            expected_artifact = configured_artifacts[group][kind]
            if artifact != expected_artifact:
                raise ValueError(
                    f"Phase metadata {group} {kind} is not the database selected "
                    "by the replay config"
                )
            if (
                not artifact.is_file()
                or sha256(artifact) != record.get(f"{kind}_sha256")
            ):
                raise ValueError(f"Phase metadata {group} {kind} hash mismatch")
    qa = metadata.get("qa", {})
    if not isinstance(qa, dict):
        raise ValueError("Phase transform QA is missing")
    mean_change = _finite_metadata_number(
        qa.get("max_abs_mean_change_pa"), "Phase transform mean-change QA"
    )
    amplitude_change = _finite_metadata_number(
        qa.get("max_abs_fundamental_amplitude_change_pa"),
        "Phase transform amplitude-change QA",
    )
    pressure_change = _finite_metadata_number(
        qa.get("max_instantaneous_pressure_change_pa"),
        "Phase transform instantaneous-change QA",
    )
    if (
        mean_change < 0.0
        or amplitude_change < 0.0
        or mean_change > TRANSFORM_INVARIANCE_TOLERANCE_PA
        or amplitude_change > TRANSFORM_INVARIANCE_TOLERANCE_PA
    ):
        raise ValueError("Phase transform changed fitted mean or amplitude")
    if pressure_change <= 0.0:
        raise ValueError("Phase transform records no positive pressure change")
    validate_phase_alignment_qa(metadata)
    return {
        "path": str(path.resolve()),
        "sha256": sha256(path),
        "surface_reference_source": {
            "path": str(registered_path),
            "size_bytes": registered_size,
            "sha256": registered_sha256,
        },
        "metadata": metadata,
    }


def ordered_frame(path: Path, count: int) -> tuple[np.ndarray, dict[str, np.ndarray]]:
    points, arrays = read_vtp(path)
    missing = [name for name in REQUIRED_ARRAYS if name not in arrays]
    if missing:
        raise ValueError(f"{path}: missing required arrays {missing}")
    raw_ids = np.asarray(arrays["ids"], dtype=np.float64).reshape(-1)
    integer_ids = np.rint(raw_ids).astype(np.int64)
    order = np.argsort(integer_ids)
    if raw_ids.size != count or not np.array_equal(
        integer_ids[order], np.arange(count)
    ):
        raise ValueError(f"{path}: particle IDs are not unique 0..N-1")
    ordered = {name: np.asarray(values)[order] for name, values in arrays.items()}
    points = np.asarray(points, dtype=np.float64)[order, :2]
    if not np.all(np.isfinite(points)):
        raise ValueError(f"{path}: particle coordinates contain NaN/Inf")
    for name in REQUIRED_ARRAYS[1:]:
        if not np.all(np.isfinite(ordered[name])):
            raise ValueError(f"{path}: {name} contains NaN/Inf")
    return points, ordered


def integrate(times: np.ndarray, values: np.ndarray) -> float:
    return float(np.trapz(values, times))


def mesh_bounds(path: Path) -> tuple[np.ndarray, np.ndarray]:
    with path.open(encoding="utf-8") as stream:
        line = next(stream)
        while line.lstrip().startswith("#"):
            line = next(stream)
        node_count = int(line.split()[0])
        nodes = np.asarray(
            [list(map(float, next(stream).split()[:2])) for _ in range(node_count)],
            dtype=np.float64,
        )
    if nodes.shape != (node_count, 2) or not np.all(np.isfinite(nodes)):
        raise ValueError(f"{path}: invalid mesh-node table")
    return np.min(nodes, axis=0), np.max(nodes, axis=0)


def validate_replay_frame_health(
    path: Path,
    points: np.ndarray,
    arrays: dict[str, np.ndarray],
    mesh_min: np.ndarray,
    mesh_max: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, int]:
    """Reject geometry or represented-volume defects before area integration."""

    coordinates = np.asarray(points, dtype=np.float64)
    volumes = np.asarray(arrays["volumes"], dtype=np.float64).reshape(-1)
    porosity = np.asarray(arrays["porosities"], dtype=np.float64).reshape(-1)
    if (
        coordinates.ndim != 2
        or coordinates.shape[1] < 2
        or coordinates.shape[0] != volumes.size
        or porosity.size != volumes.size
        or not np.all(np.isfinite(coordinates[:, :2]))
    ):
        raise ValueError(f"{path}: particle geometry is invalid")
    if not np.all(np.isfinite(volumes)) or np.any(volumes <= 0.0):
        raise ValueError(f"{path}: particle volume must be finite and positive")
    if (
        not np.all(np.isfinite(porosity))
        or np.any(porosity <= 0.0)
        or np.any(porosity >= 1.0)
    ):
        raise ValueError(f"{path}: porosity lies outside (0,1)")
    outside = np.any(
        (coordinates[:, :2] < np.asarray(mesh_min, dtype=np.float64) - 1.0e-8)
        | (coordinates[:, :2] > np.asarray(mesh_max, dtype=np.float64) + 1.0e-8),
        axis=1,
    )
    outside_count = int(np.count_nonzero(outside))
    if outside_count:
        raise ValueError(f"{path}: {outside_count} particles are outside the mesh")
    return volumes, porosity, outside_count


def read_replay(root: Path, label: str, replay: Replay) -> dict[str, Any]:
    config_path = root / "configs/phase_lag_exploratory" / label / replay.config_name
    config = json.loads(config_path.read_text(encoding="utf-8"))
    analysis = config["analysis"]
    pressure = analysis.get("prescribed_phase_pressures", {})
    parameters = replay_parameters(config, replay.role)
    if bool(config["materials"][1]["wave_pressure"]):
        raise ValueError(f"{replay.role}: replay must disable physical wave forcing")
    if pressure.get("enable") is not True or pressure.get("write") is not False:
        raise ValueError(f"{replay.role}: invalid pressure read/write mode")
    pressure_artifact_paths = configured_pressure_artifact_paths(
        root, label, role=replay.role
    )
    metadata_path = pressure_artifact_paths["values"].with_name(
        f"{pressure['file_prefix']}_metadata.json"
    )
    if replay.role == "phase_erased":
        pressure_artifact_paths["metadata"] = metadata_path
    missing_pressure_artifacts = [
        str(path) for path in pressure_artifact_paths.values() if not path.is_file()
    ]
    if missing_pressure_artifacts:
        raise FileNotFoundError(
            f"{replay.role}: missing pressure artifacts {missing_pressure_artifacts}"
        )

    coordinates = read_initial_coordinates(root / config["particles"][0]["generator"]["location"])
    count = len(coordinates)
    result_directory = root / config["post_processing"]["path"] / analysis["uuid"]
    files = sorted(result_directory.glob("particle*.vtp"), key=particle_step)
    expected_steps = expected_particle_steps(config)
    actual_steps = [particle_step(path) for path in files]
    if actual_steps != expected_steps:
        raise ValueError(
            f"{replay.role}: expected frames {expected_steps}, got {actual_steps}"
        )

    histories = sorted(result_directory.glob("pipeline-history*.csv"))
    if len(histories) != 1:
        raise ValueError(
            f"{replay.role}: expected one pipeline history, found {histories}"
        )
    history_audit = validate_pipeline_history_grid(
        read_pipeline_history(result_directory), config
    )

    mesh_min, mesh_max = mesh_bounds(root / config["mesh"]["mesh"])
    gravity = np.asarray(config["external_loading_conditions"]["gravity"], dtype=float)
    times: list[float] = []
    points_rows: list[np.ndarray] = []
    pressure_rows: list[np.ndarray] = []
    seepage_rows: list[np.ndarray] = []
    ratio_rows: list[np.ndarray] = []
    volume_rows: list[np.ndarray] = []
    velocity_maxima: list[float] = []
    displacement_maxima: list[float] = []
    porosity_minima: list[float] = []
    porosity_maxima: list[float] = []
    outside_counts: list[int] = []
    eligible_reference: np.ndarray | None = None

    for path in files:
        points, arrays = ordered_frame(path, count)
        volumes, porosity, outside_count = validate_replay_frame_health(
            path, points, arrays, mesh_min, mesh_max
        )
        stress_ratio, eligible, lost = vertical_stress_loss_state(
            np.asarray(arrays["stresses"], dtype=float)[:, 1],
            np.asarray(arrays["initial_vertical_effective_stresses"], dtype=float),
        )
        solver_ratio = np.asarray(
            arrays["vertical_effective_stress_remaining_ratios"], dtype=float
        ).reshape(-1)
        if not np.allclose(
            stress_ratio[eligible],
            solver_ratio[eligible],
            rtol=1.0e-10,
            atol=1.0e-12,
        ):
            raise ValueError(f"{path}: solver/Python R_sigma disagreement")
        if eligible_reference is None:
            eligible_reference = eligible
        elif not np.array_equal(eligible_reference, eligible):
            raise ValueError(f"{path}: R_sigma eligibility changed")
        seepage = upward_excess_seepage_force_index(
            arrays["liquid_seepage_forces"],
            arrays["liquid_densities"],
            arrays["gamma_sub"],
            gravity,
        )
        velocities = np.asarray(arrays["velocities"], dtype=float)[:, :2]
        displacements = np.asarray(arrays["displacements"], dtype=float)[:, :2]
        times.append(particle_step(path) * float(analysis["dt"]))
        points_rows.append(points)
        pressure_rows.append(
            np.asarray(arrays["PIC_pore_pressure_excess"], dtype=float).reshape(-1)
        )
        seepage_rows.append(seepage)
        ratio_rows.append(stress_ratio)
        volume_rows.append(volumes)
        velocity_maxima.append(float(np.max(np.linalg.norm(velocities, axis=1))))
        displacement_maxima.append(
            float(np.max(np.linalg.norm(displacements, axis=1)))
        )
        porosity_minima.append(float(np.min(porosity)))
        porosity_maxima.append(float(np.max(porosity)))
        outside_counts.append(outside_count)

    times_array = np.asarray(times)
    points_array = np.asarray(points_rows)
    pressure_array = np.asarray(pressure_rows)
    seepage_array = np.asarray(seepage_rows)
    ratio_array = np.asarray(ratio_rows)
    volume_array = np.asarray(volume_rows)
    assert eligible_reference is not None
    if_mask = seepage_array >= 1.0 - THRESHOLD_RATIO_TOLERANCE
    stress_mask = eligible_reference[np.newaxis, :] & (
        ratio_array <= 0.05 + THRESHOLD_RATIO_TOLERANCE
    )
    joint_mask = if_mask & stress_mask
    masks = {"IF_ge_1": if_mask, "R_sigma_le_0p05": stress_mask, "joint": joint_mask}
    criteria: dict[str, Any] = {}
    for name, mask in masks.items():
        fractions = np.mean(mask, axis=1)
        areas = np.sum(np.where(mask, volume_array, 0.0), axis=1)
        criteria[name] = {
            "maximum_fraction": float(np.max(fractions)),
            "maximum_current_area_m2": float(np.max(areas)),
            "time_integrated_fraction_s": integrate(times_array, fractions),
            "time_integrated_current_area_m2_s": integrate(times_array, areas),
            "active_saved_frames": int(np.count_nonzero(np.any(mask, axis=1))),
        }

    return {
        "role": replay.role,
        "config": config,
        "parameters": parameters,
        "config_path": config_path,
        "result_directory": result_directory,
        "files": files,
        "history_path": histories[0],
        "history_audit": history_audit,
        "pressure_artifact_paths": pressure_artifact_paths,
        "times": times_array,
        "points": points_array,
        "pressure": pressure_array,
        "seepage": seepage_array,
        "stress_ratio": ratio_array,
        "masks": masks,
        "criteria": criteria,
        "health": {
            "particle_count": count,
            "saved_frames": len(files),
            "outside_mesh_maximum_count": max(outside_counts),
            "maximum_velocity_m_s": max(velocity_maxima),
            "maximum_displacement_m": max(displacement_maxima),
            "minimum_porosity": min(porosity_minima),
            "maximum_porosity": max(porosity_maxima),
        },
    }


def robust_symmetric_limit(*fields: np.ndarray) -> float:
    values = np.concatenate([np.abs(np.asarray(field).reshape(-1)) for field in fields])
    values = values[np.isfinite(values)]
    if values.size == 0:
        raise ValueError("Cannot scale a field without finite values")
    return max(float(np.percentile(values, 99.5)), np.finfo(float).eps)


def write_history(path: Path, lagged: dict[str, Any], erased: dict[str, Any]) -> None:
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(
            [
                "time_s",
                "lagged_IF_fraction",
                "erased_IF_fraction",
                "lagged_Rsigma_fraction",
                "erased_Rsigma_fraction",
                "lagged_joint_fraction",
                "erased_joint_fraction",
            ]
        )
        for index, time in enumerate(lagged["times"]):
            writer.writerow(
                [
                    float(time),
                    float(np.mean(lagged["masks"]["IF_ge_1"][index])),
                    float(np.mean(erased["masks"]["IF_ge_1"][index])),
                    float(np.mean(lagged["masks"]["R_sigma_le_0p05"][index])),
                    float(np.mean(erased["masks"]["R_sigma_le_0p05"][index])),
                    float(np.mean(lagged["masks"]["joint"][index])),
                    float(np.mean(erased["masks"]["joint"][index])),
                ]
            )


def plot_history(
    output: Path,
    lagged: dict[str, Any],
    erased: dict[str, Any],
    parameters: ReplayParameters,
) -> None:
    figure, axes = plt.subplots(3, 1, figsize=(8.2, 8.0), sharex=True, constrained_layout=True)
    rows = (
        ("IF_ge_1", "Area fraction with upward seepage IF ≥ 1"),
        ("R_sigma_le_0p05", r"Area fraction with $R_\sigma\leq0.05$"),
        ("joint", r"Joint area fraction: IF ≥ 1 and $R_\sigma\leq0.05$"),
    )
    for axis, (name, ylabel) in zip(axes, rows):
        for case, label, colour in (
            (lagged, "Original phase lag", "#b2182b"),
            (erased, "Phase-erased", "#2166ac"),
        ):
            axis.plot(
                case["times"],
                np.mean(case["masks"][name], axis=1),
                label=label,
                color=colour,
                linewidth=1.8,
            )
        axis.set_ylabel(ylabel)
        axis.grid(alpha=0.25)
    axes[0].legend(loc="upper right")
    axes[-1].set_xlabel("Time (s)")
    figure.suptitle(f"Exploratory replay: {parameter_label(parameters)}")
    for suffix in ("png", "pdf"):
        figure.savefig(output / f"lagged_vs_phase_erased_liquefaction_history.{suffix}", dpi=220)
    plt.close(figure)


def plot_fields(
    output: Path,
    lagged: dict[str, Any],
    erased: dict[str, Any],
    index: int,
    parameters: ReplayParameters,
) -> None:
    figure, axes = plt.subplots(3, 3, figsize=(15.5, 9.6), constrained_layout=True)
    config = lagged["config"]
    center = np.asarray(config["analysis"]["rigid_pipeline"]["center"], dtype=float)
    coordinates = read_initial_coordinates(
        Path(__file__).resolve().parent
        / config["particles"][0]["generator"]["location"]
    )
    time = float(lagged["times"][index])
    fields = (
        (
            lagged["pressure"][index] / 1000.0,
            erased["pressure"][index] / 1000.0,
            "Liquid excess pressure (kPa)",
            "RdBu_r",
        ),
        (
            lagged["seepage"][index],
            erased["seepage"][index],
            "Upward seepage IF (-)",
            "magma",
        ),
        (
            lagged["stress_ratio"][index],
            erased["stress_ratio"][index],
            r"Vertical stress ratio $R_\sigma$ (-)",
            "viridis",
        ),
    )
    column_titles = ("Original phase lag", "Phase-erased", "ID-matched difference (lagged - erased)")
    for row, (left, right, label, cmap) in enumerate(fields):
        if row == 0:
            common_limit = robust_symmetric_limit(left, right)
            common_limits = (-common_limit, common_limit)
        else:
            combined = np.concatenate((left, right))
            combined = combined[np.isfinite(combined)]
            common_limits = (
                float(np.percentile(combined, 0.5)),
                float(np.percentile(combined, 99.5)),
            )
            if math.isclose(*common_limits):
                common_limits = (common_limits[0] - 1.0, common_limits[1] + 1.0)
        difference = left - right
        difference_limit = robust_symmetric_limit(difference)
        artists = []
        for column, values in enumerate((left, right)):
            artists.append(
                draw_field(
                    axes[row, column],
                    coordinates,
                    values,
                    config,
                    center,
                    title=column_titles[column] if row == 0 else "",
                    cmap=cmap,
                    limits=common_limits,
                )
            )
        difference_artist = draw_field(
            axes[row, 2],
            coordinates,
            difference,
            config,
            center,
            title=column_titles[2] if row == 0 else "",
            cmap="RdBu_r",
            limits=(-difference_limit, difference_limit),
        )
        axes[row, 0].set_ylabel(f"{label}\ny (m)")
        figure.colorbar(artists[0], ax=axes[row, :2], label=label, shrink=0.82)
        figure.colorbar(
            difference_artist,
            ax=axes[row, 2],
            label=f"Δ {label}",
            shrink=0.82,
        )
    figure.suptitle(
        f"Exploratory lagged / phase-erased 2-D comparison at t={time:.3f} s\n"
        f"{parameter_label(parameters)}; fields shown at initial material-point coordinates"
    )
    for suffix in ("png", "pdf"):
        figure.savefig(output / f"lagged_vs_phase_erased_2d.{suffix}", dpi=220)
    plt.close(figure)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--label", required=True)
    args = parser.parse_args()
    root = Path(__file__).resolve().parent
    output = root / "analysis/phase_lag_exploratory" / args.label / "replay_comparison"
    output.mkdir(parents=True, exist_ok=True)
    lagged = read_replay(root, args.label, REPLAYS[0])
    erased = read_replay(root, args.label, REPLAYS[1])
    parameters = validate_replay_pair_configs(lagged["config"], erased["config"])
    phase_metadata_audit = validate_phase_metadata(root, args.label)
    if not np.array_equal(lagged["times"], erased["times"]):
        raise ValueError("Lagged and phase-erased saved-time grids differ")
    pressure_difference_q95 = np.percentile(
        np.abs(lagged["pressure"] - erased["pressure"]), 95, axis=1
    )
    selected_index = int(np.argmax(pressure_difference_q95))
    write_history(output / "lagged_vs_phase_erased_history.csv", lagged, erased)
    plot_history(output, lagged, erased, parameters)
    plot_fields(output, lagged, erased, selected_index, parameters)

    comparison: dict[str, Any] = {}
    for criterion in ("IF_ge_1", "R_sigma_le_0p05", "joint"):
        comparison[criterion] = {
            "lagged": lagged["criteria"][criterion],
            "phase_erased": erased["criteria"][criterion],
            "lagged_minus_phase_erased_time_integrated_fraction_s": (
                lagged["criteria"][criterion]["time_integrated_fraction_s"]
                - erased["criteria"][criterion]["time_integrated_fraction_s"]
            ),
            "lagged_minus_phase_erased_time_integrated_area_m2_s": (
                lagged["criteria"][criterion]["time_integrated_current_area_m2_s"]
                - erased["criteria"][criterion]["time_integrated_current_area_m2_s"]
            ),
        }
    audit = {
        "schema": SCHEMA,
        "status": "exploratory-only-not-registered-manuscript-evidence",
        "solver_binary": {
            "path": str(
                (root.parents[1] / "mpm_hpc_source/build-pipeline/mpm").resolve()
            ),
            "sha256": sha256(
                root.parents[1] / "mpm_hpc_source/build-pipeline/mpm"
            ),
        },
        "parameters": parameters.as_dict(),
        "parameter_label": parameter_label(parameters),
        "frame_selection": {
            "rule": "maximum across frames of particlewise q95 absolute lagged-minus-erased excess-pressure difference",
            "time_s": float(lagged["times"][selected_index]),
            "q95_absolute_pressure_difference_pa": float(
                pressure_difference_q95[selected_index]
            ),
        },
        "phase_controls_producer": phase_metadata_audit["metadata"]["producer"],
        "phase_transform": phase_metadata_audit,
        "cases": {
            case["role"]: {
                "config": str(case["config_path"].resolve()),
                "config_sha256": sha256(case["config_path"]),
                "result_directory": str(case["result_directory"].resolve()),
                "health": case["health"],
                "criteria": case["criteria"],
                "frame_sha256": [sha256(path) for path in case["files"]],
                "pipeline_history": {
                    "path": str(case["history_path"].resolve()),
                    "sha256": sha256(case["history_path"]),
                    "audit": case["history_audit"],
                },
                "pressure_artifacts": {
                    name: {
                        "path": str(path.resolve()),
                        "sha256": sha256(path),
                        "size_bytes": path.stat().st_size,
                    }
                    for name, path in case["pressure_artifact_paths"].items()
                },
            }
            for case in (lagged, erased)
        },
        "comparison": comparison,
        "interpretation_guardrail": (
            "The replay pair is a fixed-pipeline numerical counterfactual. A "
            "liquefaction-effect statement requires the hydraulic IF gate and "
            "the same-particle joint IF/R_sigma evidence; raw spatial roughness "
            "or a selected 2-D frame alone is insufficient."
        ),
    }
    (output / "lagged_vs_phase_erased.audit.json").write_text(
        json.dumps(audit, indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )
    print(json.dumps({"frame_selection": audit["frame_selection"], "comparison": comparison}, indent=2))
    print(f"Wrote exploratory replay comparison to {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
