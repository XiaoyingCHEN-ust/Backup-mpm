#!/usr/bin/env python3
"""Plot audited two-dimensional liquid-pressure phase-lag diagnostics.

The script uses the completed physical HS screen, subtracts the audited HS_EQ
``PIC_liquid_pressures`` checkpoint by particle ID, and fits the fundamental
wave component over the registered pre-release window.  It produces common-
scale quarter-cycle fields, a local-surface-referenced phase-difference field,
and a time--depth section at the pre-registered far-field x location.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Mapping

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle

from analyze_study import build_probes, expected_particle_steps, read_vtp
from plot_manuscript_figures import _masked_triangulation


SCHEMA = "pipeline-liquid-pressure-phase-lag-audit-v3"
ARTIFACT_SCHEMA = "pipeline-liquid-pressure-phase-lag-artifacts-v2"
INPUT_SCHEMA = "pipeline-liquid-pressure-phase-lag-inputs-v1"
EXPECTED_PHASE_FRAMES = 60
CASE = "HS"
WAVE_PERIOD_S = 1.3
FIT_START_S = 1.3
FIT_END_S = 3.9
SNAPSHOT_TIMES_S = (2.6, 2.925, 3.25, 3.575)
MINIMUM_LOCAL_AMPLITUDE_RATIO = 0.05
MINIMUM_HARMONIC_R_SQUARED = 0.80
FARFIELD_X_M = 1.06
PLOT_SCRIPT_PATH = Path(__file__).resolve()
ANALYZER_PATH = PLOT_SCRIPT_PATH.with_name("analyze_study.py")


class PhaseAuditError(RuntimeError):
    """Raised when phase-field artifact provenance is incomplete or stale."""


@dataclass(frozen=True)
class PlotParameters:
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


def plot_parameters(
    config: dict[str, Any], *, require_unsmoothed: bool = False
) -> PlotParameters:
    materials = config.get("materials")
    if not isinstance(materials, list) or len(materials) != 2:
        raise ValueError("Phase-lag config must define exactly two materials")
    soil, fluid = materials
    try:
        permeability = float(soil["intrinsic_permeability"])
        liquid_saturation = float(fluid["liquid_saturation"])
        gas_saturation = float(fluid["gas_saturation"])
        analysis = config["analysis"]
        pressure_smoothing = analysis["pressure_smoothing"]
        pressure_smoothing_in_loop = analysis["pressure_smoothing_in_loop"]
        pipeline_fixed = analysis["rigid_pipeline"]["fixed"]
    except (KeyError, TypeError, ValueError) as error:
        raise ValueError("Phase-lag config has incomplete hydraulic parameters") from error
    if not math.isfinite(permeability) or permeability <= 0.0:
        raise ValueError("Intrinsic permeability must be finite and positive")
    if not math.isfinite(liquid_saturation) or not 0.0 < liquid_saturation < 1.0:
        raise ValueError("Liquid saturation must lie strictly inside (0,1)")
    if not math.isfinite(gas_saturation) or not math.isclose(
        liquid_saturation + gas_saturation,
        1.0,
        rel_tol=0.0,
        abs_tol=1.0e-12,
    ):
        raise ValueError("Liquid and gas saturations must sum to one")
    for name, value in (
        ("pressure_smoothing", pressure_smoothing),
        ("pressure_smoothing_in_loop", pressure_smoothing_in_loop),
        ("rigid_pipeline.fixed", pipeline_fixed),
    ):
        if type(value) is not bool:
            raise ValueError(f"{name} must be an explicit JSON boolean")
    if require_unsmoothed and (
        pressure_smoothing is not False
        or pressure_smoothing_in_loop is not False
    ):
        raise ValueError(
            "Exploratory phase-lag plotting requires both pressure-smoothing "
            "switches to be false"
        )
    if require_unsmoothed and pipeline_fixed is not True:
        raise ValueError(
            "Exploratory phase-lag plotting requires an explicitly fixed pipeline"
        )
    return PlotParameters(
        intrinsic_permeability_m2=permeability,
        liquid_saturation=liquid_saturation,
        gas_saturation=gas_saturation,
        pressure_smoothing=pressure_smoothing,
        pressure_smoothing_in_loop=pressure_smoothing_in_loop,
        pipeline_fixed=pipeline_fixed,
    )


def parameter_subtitle(parameters: PlotParameters) -> str:
    return (
        rf"$k={parameters.intrinsic_permeability_m2:.3e}\ \mathrm{{m}}^2$, "
        rf"$S_w={parameters.liquid_saturation:.3f}$, "
        rf"$S_g={parameters.gas_saturation:.3f}$; "
        f"smoothing(output={str(parameters.pressure_smoothing).lower()}, "
        f"in-loop={str(parameters.pressure_smoothing_in_loop).lower()}); "
        f"pipeline fixed={str(parameters.pipeline_fixed).lower()}"
    )


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def artifact_record(path: Path) -> dict[str, str | int]:
    """Describe one non-empty regular file using canonical provenance."""

    resolved = path.resolve()
    if not resolved.is_file():
        raise PhaseAuditError(f"Required artifact is missing: {resolved}")
    size = resolved.stat().st_size
    if size <= 0:
        raise PhaseAuditError(f"Required artifact is empty: {resolved}")
    return {
        "path": str(resolved),
        "size_bytes": size,
        "sha256": sha256(resolved),
    }


def _validated_artifact_record(
    value: Any, *, name: str, expected_path: Path | None = None
) -> Path:
    if not isinstance(value, dict) or set(value) != {
        "path",
        "size_bytes",
        "sha256",
    }:
        raise PhaseAuditError(
            f"{name} must contain exactly path, size_bytes, and sha256"
        )
    raw_path = value["path"]
    if not isinstance(raw_path, str) or not raw_path:
        raise PhaseAuditError(f"{name} path must be a non-empty string")
    path = Path(raw_path)
    if not path.is_absolute() or raw_path != str(path.resolve()):
        raise PhaseAuditError(f"{name} path is not canonical and absolute")
    if expected_path is not None and path != expected_path.resolve():
        raise PhaseAuditError(
            f"{name} path mismatch: expected {expected_path.resolve()}, got {path}"
        )
    if not path.is_file():
        raise PhaseAuditError(f"{name} is missing: {path}")
    recorded_size = value["size_bytes"]
    if type(recorded_size) is not int or recorded_size <= 0:
        raise PhaseAuditError(f"{name} size_bytes must be a positive integer")
    observed_size = path.stat().st_size
    if observed_size != recorded_size:
        raise PhaseAuditError(
            f"{name} size mismatch: recorded {recorded_size}, observed {observed_size}"
        )
    recorded_digest = value["sha256"]
    if not isinstance(recorded_digest, str) or re.fullmatch(
        r"[0-9a-f]{64}", recorded_digest
    ) is None:
        raise PhaseAuditError(f"{name} sha256 must be 64 lowercase hex characters")
    observed_digest = sha256(path)
    if observed_digest != recorded_digest:
        raise PhaseAuditError(f"{name} sha256 mismatch: {path}")
    if path.stat().st_size != observed_size:
        raise PhaseAuditError(f"{name} changed while it was being validated: {path}")
    return path


def checkpoint_config_path(
    analysis_config_path: Path, config: dict[str, Any]
) -> Path:
    """Resolve the unique config that produced the supplied resume checkpoint."""

    try:
        resume = config["analysis"]["resume"]
        checkpoint_uuid = str(resume["uuid"])
        checkpoint_step = int(resume["step"])
        enabled = resume["resume"]
    except (KeyError, TypeError, ValueError) as error:
        raise PhaseAuditError(
            "Analysis config does not identify its equilibrium checkpoint"
        ) from error
    if enabled is not True or not checkpoint_uuid or checkpoint_step <= 0:
        raise PhaseAuditError(
            "Analysis config must explicitly resume a positive-step checkpoint"
        )

    matches: list[Path] = []
    for candidate in sorted(analysis_config_path.resolve().parent.glob("*.json")):
        if candidate.resolve() == analysis_config_path.resolve():
            continue
        try:
            payload = read_json(candidate)
            analysis = payload["analysis"]
            uuid = str(analysis["uuid"])
            nsteps = int(analysis["nsteps"])
        except (OSError, json.JSONDecodeError, KeyError, TypeError, ValueError):
            continue
        if uuid == checkpoint_uuid and nsteps == checkpoint_step:
            matches.append(candidate.resolve())
    if len(matches) != 1:
        raise PhaseAuditError(
            "Expected exactly one config for checkpoint "
            f"{checkpoint_uuid}@{checkpoint_step}, found {len(matches)}"
        )
    return matches[0]


def _expected_input_path(
    expected: Mapping[str, Any] | None, name: str
) -> Path | None:
    if expected is None or name not in expected:
        return None
    value = expected[name]
    if not isinstance(value, Path):
        raise PhaseAuditError(f"Expected input path {name} must be a Path")
    return value.resolve()


def _validated_phase_inputs(
    audit: dict[str, Any],
    inputs: Any,
    *,
    expected_input_paths: Mapping[str, Any] | None,
) -> dict[str, Any]:
    """Validate every byte-bearing input used to construct a phase field."""

    required = {
        "schema",
        "analysis_config",
        "checkpoint_config",
        "checkpoint",
        "particles",
        "pipeline_history",
        "result_directory",
        "result_particle_vtp",
    }
    if not isinstance(inputs, dict) or set(inputs) != required:
        raise PhaseAuditError(
            "phase inputs must contain exactly schema, both configs, checkpoint, "
            "particles, pipeline history, result directory, and result VTPs"
        )
    if inputs["schema"] != INPUT_SCHEMA:
        raise PhaseAuditError(f"Expected phase input schema {INPUT_SCHEMA}")

    records = {
        name: _validated_artifact_record(
            inputs[name],
            name=f"inputs.{name}",
            expected_path=_expected_input_path(expected_input_paths, name),
        )
        for name in (
            "analysis_config",
            "checkpoint_config",
            "checkpoint",
            "particles",
            "pipeline_history",
        )
    }
    if records["analysis_config"] == records["checkpoint_config"]:
        raise PhaseAuditError("Analysis and checkpoint configs must be distinct")

    raw_result_directory = inputs["result_directory"]
    if not isinstance(raw_result_directory, str) or not raw_result_directory:
        raise PhaseAuditError("inputs.result_directory must be an absolute path")
    result_directory = Path(raw_result_directory)
    if (
        not result_directory.is_absolute()
        or raw_result_directory != str(result_directory.resolve())
        or not result_directory.is_dir()
    ):
        raise PhaseAuditError(
            "inputs.result_directory must be a canonical existing directory"
        )
    expected_result_directory = _expected_input_path(
        expected_input_paths, "result_directory"
    )
    if (
        expected_result_directory is not None
        and result_directory != expected_result_directory
    ):
        raise PhaseAuditError("inputs.result_directory path mismatch")

    raw_frames = inputs["result_particle_vtp"]
    if not isinstance(raw_frames, list) or len(raw_frames) != EXPECTED_PHASE_FRAMES:
        raise PhaseAuditError(
            f"inputs.result_particle_vtp must contain exactly {EXPECTED_PHASE_FRAMES} records"
        )
    expected_frames: tuple[Path, ...] | None = None
    if expected_input_paths is not None and "result_particle_vtp" in expected_input_paths:
        raw_expected = expected_input_paths["result_particle_vtp"]
        if not isinstance(raw_expected, (list, tuple)):
            raise PhaseAuditError("Expected result_particle_vtp paths must be a sequence")
        expected_frames = tuple(Path(path).resolve() for path in raw_expected)
        if len(expected_frames) != EXPECTED_PHASE_FRAMES:
            raise PhaseAuditError(
                f"Expected result_particle_vtp must contain {EXPECTED_PHASE_FRAMES} paths"
            )
    frames = tuple(
        _validated_artifact_record(
            value,
            name=f"inputs.result_particle_vtp[{index}]",
            expected_path=(None if expected_frames is None else expected_frames[index]),
        )
        for index, value in enumerate(raw_frames)
    )
    if len(set(frames)) != EXPECTED_PHASE_FRAMES:
        raise PhaseAuditError("Result particle VTP paths must be unique")
    if any(path.parent != result_directory for path in frames):
        raise PhaseAuditError("Every result particle VTP must be in result_directory")
    try:
        ordered_frames = tuple(
            sorted(result_directory.glob("particle*.vtp"), key=particle_step)
        )
    except ValueError as error:
        raise PhaseAuditError("Result directory contains an invalid particle VTP name") from error
    if ordered_frames != frames:
        raise PhaseAuditError(
            "Recorded result particle VTPs are not the exact ordered directory set"
        )
    if int(audit.get("frame_count", -1)) != EXPECTED_PHASE_FRAMES:
        raise PhaseAuditError(
            f"Phase audit frame_count must equal {EXPECTED_PHASE_FRAMES}"
        )
    if len({*records.values(), *frames}) != len(records) + len(frames):
        raise PhaseAuditError("Phase input artifact paths must be unique")
    return {
        **records,
        "result_directory": result_directory,
        "result_particle_vtp": frames,
    }


def validate_phase_audit_artifacts(
    audit: dict[str, Any],
    *,
    expected_figure_paths: Iterable[Path],
    analyzer_path: Path = ANALYZER_PATH,
    plot_script_path: Path = PLOT_SCRIPT_PATH,
    expected_input_paths: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Fail closed unless phase-v3 binds every input and exact four plot files."""

    if not isinstance(audit, dict):
        raise PhaseAuditError("Phase audit must be a JSON object")
    if audit.get("schema") != SCHEMA:
        raise PhaseAuditError(f"Expected phase audit schema {SCHEMA}")
    expected_figures = tuple(path.resolve() for path in expected_figure_paths)
    if (
        len(expected_figures) != 4
        or len(set(expected_figures)) != 4
        or sorted(path.suffix.lower() for path in expected_figures)
        != [".pdf", ".pdf", ".png", ".png"]
    ):
        raise PhaseAuditError(
            "Expected exactly four unique figure paths (two PNG and two PDF)"
        )
    artifacts = audit.get("artifacts")
    if not isinstance(artifacts, dict) or set(artifacts) != {
        "schema",
        "analyzer",
        "plot_script",
        "figures",
        "inputs",
    }:
        raise PhaseAuditError(
            "artifacts must contain exactly schema, analyzer, plot_script, figures, and inputs"
        )
    if artifacts["schema"] != ARTIFACT_SCHEMA:
        raise PhaseAuditError(f"Expected artifact schema {ARTIFACT_SCHEMA}")
    analyzer = _validated_artifact_record(
        artifacts["analyzer"], name="analyzer", expected_path=analyzer_path
    )
    plot_script = _validated_artifact_record(
        artifacts["plot_script"],
        name="plot_script",
        expected_path=plot_script_path,
    )
    figures = artifacts["figures"]
    if not isinstance(figures, list) or len(figures) != 4:
        raise PhaseAuditError("figures must contain exactly four artifact records")
    observed_figures = [
        _validated_artifact_record(value, name=f"figures[{index}]")
        for index, value in enumerate(figures)
    ]
    if len(set(observed_figures)) != len(observed_figures):
        raise PhaseAuditError("Figure artifact paths must be unique")
    if set(observed_figures) != set(expected_figures):
        raise PhaseAuditError("Figure artifact paths do not match the expected set")
    if len({analyzer, plot_script, *observed_figures}) != 6:
        raise PhaseAuditError("All analyzer, plot-script, and figure paths must be unique")
    inputs = _validated_phase_inputs(
        audit,
        artifacts["inputs"],
        expected_input_paths=expected_input_paths,
    )
    if audit.get("config") != artifacts["inputs"]["analysis_config"]:
        raise PhaseAuditError(
            "Top-level config record must exactly alias inputs.analysis_config"
        )
    if audit.get("checkpoint") != artifacts["inputs"]["checkpoint"]:
        raise PhaseAuditError(
            "Top-level checkpoint record must exactly alias inputs.checkpoint"
        )
    if set((analyzer, plot_script, *observed_figures)) & {
        *(
            path
            for name, path in inputs.items()
            if name != "result_particle_vtp"
        ),
        *inputs["result_particle_vtp"],
    }:
        raise PhaseAuditError("Phase input and output artifact paths must be disjoint")
    return inputs


def atomic_write_json(path: Path, value: dict[str, Any]) -> None:
    """Atomically replace a JSON audit without exposing a partial document."""

    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    try:
        with temporary.open("x", encoding="utf-8", newline="\n") as stream:
            json.dump(value, stream, indent=2, allow_nan=False)
            stream.write("\n")
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
    finally:
        if temporary.exists():
            temporary.unlink()


def figure_output_paths(output: Path) -> tuple[Path, Path]:
    """Return the exact PNG/PDF paths written for one figure base path."""

    return output.with_suffix(".png"), output.with_suffix(".pdf")


def particle_step(path: Path) -> int:
    match = re.fullmatch(r"particle(\d+)\.vtp", path.name)
    if match is None:
        raise ValueError(f"Unexpected particle filename: {path}")
    return int(match.group(1))


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def read_initial_coordinates(path: Path) -> np.ndarray:
    coordinates = np.loadtxt(path, skiprows=1, dtype=np.float64)
    if coordinates.ndim != 2 or coordinates.shape[1] != 2:
        raise ValueError(f"Invalid particle coordinate table: {path}")
    return coordinates


def align_liquid_pressure(
    path: Path, expected_ids: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    points, arrays = read_vtp(path)
    if "ids" not in arrays or "PIC_liquid_pressures" not in arrays:
        raise KeyError(f"{path} lacks ids or PIC_liquid_pressures")
    raw_ids = np.asarray(arrays["ids"], dtype=np.float64).reshape(-1)
    ids = np.rint(raw_ids).astype(np.int64)
    if not np.all(np.isfinite(raw_ids)) or not np.array_equal(
        raw_ids, ids.astype(np.float64)
    ):
        raise ValueError(f"{path} contains non-finite or non-integral IDs")
    order = np.argsort(ids)
    if not np.array_equal(ids[order], expected_ids):
        raise ValueError(f"Particle IDs in {path} do not match the checkpoint")
    pressure = np.asarray(arrays["PIC_liquid_pressures"], dtype=np.float64)[order]
    aligned_points = np.asarray(points, dtype=np.float64)[order, :2]
    if not np.all(np.isfinite(pressure)) or not np.all(np.isfinite(aligned_points)):
        raise ValueError(f"{path} contains non-finite pressure or coordinates")
    return aligned_points, pressure


def harmonic_fields(
    times: np.ndarray, values: np.ndarray, period: float
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    omega = 2.0 * math.pi / period
    design = np.column_stack(
        (np.ones_like(times), np.cos(omega * times), np.sin(omega * times))
    )
    coefficients = np.linalg.pinv(design) @ values
    mean = coefficients[0]
    amplitude = np.hypot(coefficients[1], coefficients[2])
    phase = np.arctan2(coefficients[2], coefficients[1])
    return mean, amplitude, phase


def harmonic_r_squared(
    times: np.ndarray, values: np.ndarray, period: float
) -> np.ndarray:
    """Return per-particle coherence with the fitted fundamental harmonic."""

    omega = 2.0 * math.pi / period
    design = np.column_stack(
        (np.ones_like(times), np.cos(omega * times), np.sin(omega * times))
    )
    fitted = design @ (np.linalg.pinv(design) @ values)
    residual = np.sum(np.square(values - fitted), axis=0)
    total = np.sum(np.square(values - np.mean(values, axis=0)), axis=0)
    return 1.0 - residual / np.maximum(total, np.finfo(float).eps)


def wrap_degrees(values: np.ndarray) -> np.ndarray:
    return np.degrees((values + math.pi) % (2.0 * math.pi) - math.pi)


def same_column_surface_values(
    coordinates: np.ndarray, values: np.ndarray
) -> np.ndarray:
    output = np.empty_like(values, dtype=np.float64)
    for x in np.unique(coordinates[:, 0]):
        indices = np.flatnonzero(np.isclose(coordinates[:, 0], x, atol=1.0e-10))
        surface = indices[np.argmax(coordinates[indices, 1])]
        output[indices] = values[surface]
    return output


def pipeline_pose(
    history: list[dict[str, float]], time_s: float, initial_center: np.ndarray
) -> np.ndarray:
    row = min(history, key=lambda item: abs(item["time"] - time_s))
    return initial_center + np.asarray([row["displacement_x"], row["displacement_y"]])


def read_pipeline_history(path: Path) -> list[dict[str, float]]:
    import csv

    with path.open(newline="", encoding="utf-8") as stream:
        rows = [{key: float(value) for key, value in row.items()} for row in csv.DictReader(stream)]
    if not rows:
        raise ValueError(f"No pipeline history rows in {path}")
    return rows


def draw_field(
    axis: plt.Axes,
    points: np.ndarray,
    values: np.ndarray,
    config: dict[str, Any],
    center: np.ndarray,
    *,
    title: str,
    cmap: str,
    limits: tuple[float, float],
) -> Any:
    triangulation = _masked_triangulation(points, config, center)
    finite = np.isfinite(values)
    if not np.all(finite):
        invalid = np.any(~finite[triangulation.triangles], axis=1)
        existing = triangulation.mask
        triangulation.set_mask(invalid if existing is None else (existing | invalid))
    artist = axis.tripcolor(
        triangulation,
        np.where(finite, values, 0.0),
        shading="gouraud",
        cmap=cmap,
        vmin=limits[0],
        vmax=limits[1],
        rasterized=True,
    )
    radius = float(config["analysis"]["rigid_pipeline"]["outer_radius"])
    axis.add_patch(Circle(center, radius, facecolor="white", edgecolor="black", lw=0.9))
    axis.set_xlim(0.0, 1.5)
    axis.set_ylim(0.0, 0.52)
    axis.set_aspect("equal")
    axis.set_xlabel("x (m)")
    axis.set_ylabel("y (m)")
    axis.set_title(title, fontsize=10)
    return artist


def plot_fields(
    output: Path,
    config: dict[str, Any],
    snapshots: list[dict[str, Any]],
    initial_coordinates: np.ndarray,
    amplitude_ratio: np.ndarray,
    phase_difference_deg: np.ndarray,
    case_label: str,
    parameters: PlotParameters,
) -> None:
    figure = plt.figure(figsize=(16.0, 8.0), constrained_layout=True)
    grid = figure.add_gridspec(2, 4, height_ratios=(1.0, 1.12))
    snapshot_axes = [figure.add_subplot(grid[0, column]) for column in range(4)]
    amplitude_axis = figure.add_subplot(grid[1, 0:2])
    phase_axis = figure.add_subplot(grid[1, 2:4])
    pressure_limit = max(
        float(np.max(np.abs(snapshot["pressure_excess_pa"]))) for snapshot in snapshots
    )
    pressure_artist = None
    for column, snapshot in enumerate(snapshots):
        pressure_artist = draw_field(
            snapshot_axes[column],
            snapshot["points"],
            snapshot["pressure_excess_pa"] / 1000.0,
            config,
            snapshot["pipeline_center"],
            title=(
                f"t={snapshot['time_s']:.3f} s "
                f"(phase={snapshot['phase_fraction']:.2f}T)"
            ),
            cmap="RdBu_r",
            limits=(-pressure_limit / 1000.0, pressure_limit / 1000.0),
        )
    assert pressure_artist is not None
    figure.colorbar(
        pressure_artist,
        ax=snapshot_axes,
        label="Liquid excess pressure, $p_l-p_{l,EQ}$ (kPa)",
        shrink=0.86,
    )

    center = np.asarray(config["analysis"]["rigid_pipeline"]["center"], dtype=float)
    amplitude_max = max(1.0, float(np.nanmax(amplitude_ratio)))
    amplitude_artist = draw_field(
        amplitude_axis,
        initial_coordinates,
        amplitude_ratio,
        config,
        center,
        title="Fundamental amplitude relative to same-column surface",
        cmap="viridis",
        limits=(0.0, amplitude_max),
    )
    figure.colorbar(
        amplitude_artist, ax=amplitude_axis, label="Amplitude ratio", shrink=0.85
    )
    phase_artist = draw_field(
        phase_axis,
        initial_coordinates,
        phase_difference_deg,
        config,
        center,
        title="Fundamental phase difference from same-column surface",
        cmap="twilight_shifted",
        limits=(-180.0, 180.0),
    )
    figure.colorbar(
        phase_artist,
        ax=phase_axis,
        label="Phase difference (degrees; positive = delayed)",
        shrink=0.85,
    )
    probes = build_probes(config, initial_coordinates)
    annotation_offsets = {
        "crown": (-8, -20),
        "shoulder": (12, 10),
        "invert": (10, -16),
        "farfield": (8, 8),
    }
    for probe in probes:
        if probe.name == "surface":
            continue
        value = float(phase_difference_deg[probe.index])
        phase_axis.scatter(
            [probe.actual_x],
            [probe.actual_y],
            marker="x",
            s=22,
            linewidths=1.0,
            color="black",
            zorder=5,
        )
        phase_axis.annotate(
            (
                f"{probe.name} {value:+.1f}°"
                if math.isfinite(value)
                else f"{probe.name} masked (amplitude/$R^2$ gate)"
            ),
            (probe.actual_x, probe.actual_y),
            xytext=annotation_offsets[probe.name],
            textcoords="offset points",
            fontsize=7.2,
            color="black",
            ha="center" if probe.name == "crown" else "left",
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.78, "pad": 1.2},
            zorder=6,
        )
    figure.suptitle(
        f"{case_label} two-dimensional PIC liquid-pressure response\n"
        f"{parameter_subtitle(parameters)}\n"
        "common snapshot scale; harmonic fit 1.3–3.9 s; phase masked where "
        "amplitude <5% of surface or harmonic $R^2<0.8$",
        fontsize=13,
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    png_path, pdf_path = figure_output_paths(output)
    figure.savefig(png_path, dpi=300, bbox_inches="tight")
    figure.savefig(pdf_path, bbox_inches="tight")
    plt.close(figure)


def plot_time_depth(
    output: Path,
    times: np.ndarray,
    pressure_excess: np.ndarray,
    initial_coordinates: np.ndarray,
    amplitude: np.ndarray,
    phase: np.ndarray,
    harmonic_r2: np.ndarray,
    case_label: str,
    parameters: PlotParameters,
) -> dict[str, Any]:
    x_values = np.unique(initial_coordinates[:, 0])
    chosen_x = float(x_values[np.argmin(np.abs(x_values - FARFIELD_X_M))])
    column = np.isclose(initial_coordinates[:, 0], chosen_x, atol=1.0e-10)
    indices = np.flatnonzero(column)
    indices = indices[np.argsort(initial_coordinates[indices, 1])]
    depth = 0.5 - initial_coordinates[indices, 1]
    selected_time = (times >= FIT_START_S - 1.0e-12) & (times <= FIT_END_S + 1.0e-12)
    section = pressure_excess[selected_time][:, indices].T
    limit = float(np.max(np.abs(section)))
    surface = indices[np.argmin(depth)]
    phase_difference = wrap_degrees(phase[indices] - phase[surface])
    amplitude_ratio = amplitude[indices] / max(amplitude[surface], np.finfo(float).eps)
    phase_difference[
        (amplitude_ratio < MINIMUM_LOCAL_AMPLITUDE_RATIO)
        | (harmonic_r2[indices] < MINIMUM_HARMONIC_R_SQUARED)
    ] = np.nan

    figure, axes = plt.subplots(1, 2, figsize=(11.5, 5.4), constrained_layout=True)
    artist = axes[0].pcolormesh(
        times[selected_time],
        depth,
        section,
        shading="auto",
        cmap="RdBu_r",
        vmin=-limit,
        vmax=limit,
        rasterized=True,
    )
    axes[0].invert_yaxis()
    axes[0].set_xlabel("Time (s)")
    axes[0].set_ylabel("Depth below seabed (m)")
    axes[0].set_title(f"Liquid excess pressure at x={chosen_x:.3f} m")
    figure.colorbar(artist, ax=axes[0], label="$p_l-p_{l,EQ}$ (Pa)")

    axes[1].plot(phase_difference, depth, color="#0072B2", lw=1.8, label="Phase difference")
    axes[1].axvline(0.0, color="0.4", lw=0.8, ls="--")
    axes[1].invert_yaxis()
    axes[1].set_xlabel("Phase difference from surface (degrees)")
    axes[1].set_ylabel("Depth below seabed (m)")
    axes[1].set_title("Fundamental phase profile")
    amplitude_axis = axes[1].twiny()
    amplitude_axis.plot(
        amplitude_ratio,
        depth,
        color="#E69F00",
        lw=1.3,
        ls="--",
        label="Amplitude ratio",
    )
    amplitude_axis.set_xlabel("Amplitude / surface amplitude", color="#A35F00")
    amplitude_axis.tick_params(axis="x", colors="#A35F00")
    figure.suptitle(
        f"{case_label} liquid-pressure time–depth propagation "
        "(harmonic-fit window)\n"
        f"{parameter_subtitle(parameters)}",
        fontsize=13,
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    png_path, pdf_path = figure_output_paths(output)
    figure.savefig(png_path, dpi=300, bbox_inches="tight")
    figure.savefig(pdf_path, bbox_inches="tight")
    plt.close(figure)
    return {
        "requested_x_m": FARFIELD_X_M,
        "selected_x_m": chosen_x,
        "phase_difference_min_deg": float(np.nanmin(phase_difference)),
        "phase_difference_max_deg": float(np.nanmax(phase_difference)),
        "maximum_amplitude_ratio": float(np.nanmax(amplitude_ratio)),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("analysis/screen/liquid_pressure_phase_lag"),
    )
    parser.add_argument(
        "--config", type=Path, default=Path("configs/screen/02_HS.json")
    )
    parser.add_argument(
        "--result-dir",
        type=Path,
        default=Path("results/screen/PLP_SCREEN_HS_SANI"),
    )
    parser.add_argument(
        "--checkpoint",
        type=Path,
        default=Path("results/screen/PLP_SCREEN_HS_EQ/particle40000.vtp"),
    )
    parser.add_argument("--particles", type=Path, default=Path("particles.txt"))
    parser.add_argument("--case-label", default=CASE)
    parser.add_argument(
        "--require-unsmoothed",
        action="store_true",
        help=(
            "Require both pressure-smoothing switches to be explicit false and "
            "the rigid pipeline to be explicitly fixed (exploratory runs)."
        ),
    )
    args = parser.parse_args()

    root = Path(__file__).resolve().parent
    resolve = lambda path: path if path.is_absolute() else root / path
    config_path = resolve(args.config)
    config = read_json(config_path)
    checkpoint_config = checkpoint_config_path(config_path, config)
    parameters = plot_parameters(
        config, require_unsmoothed=args.require_unsmoothed
    )
    result = resolve(args.result_dir)
    checkpoint = resolve(args.checkpoint)
    particle_file = resolve(args.particles)
    history_paths = sorted(result.glob("pipeline-history*.csv"))
    if len(history_paths) != 1:
        raise ValueError(
            f"Expected exactly one HS pipeline history, found {len(history_paths)}"
        )
    history_path = history_paths[0]
    output = resolve(args.output_dir).resolve()

    initial_coordinates = read_initial_coordinates(particle_file)
    _, checkpoint_arrays = read_vtp(checkpoint)
    raw_checkpoint_ids = np.asarray(
        checkpoint_arrays["ids"], dtype=np.float64
    ).reshape(-1)
    checkpoint_ids = np.rint(raw_checkpoint_ids).astype(np.int64)
    if not np.all(np.isfinite(raw_checkpoint_ids)) or not np.array_equal(
        raw_checkpoint_ids, checkpoint_ids.astype(np.float64)
    ):
        raise ValueError("Checkpoint particle IDs are non-finite or non-integral")
    checkpoint_order = np.argsort(checkpoint_ids)
    expected_ids = checkpoint_ids[checkpoint_order]
    if not np.array_equal(expected_ids, np.arange(initial_coordinates.shape[0])):
        raise ValueError("Checkpoint particle IDs are not contiguous 0..N-1")
    checkpoint_pressure = np.asarray(
        checkpoint_arrays["PIC_liquid_pressures"], dtype=np.float64
    )[checkpoint_order]
    if checkpoint_pressure.size != initial_coordinates.shape[0] or not np.all(
        np.isfinite(checkpoint_pressure)
    ):
        raise ValueError("Checkpoint liquid pressure is invalid")

    files = sorted(result.glob("particle*.vtp"), key=particle_step)
    expected_steps = expected_particle_steps(config)
    actual_steps = [particle_step(path) for path in files]
    if actual_steps != expected_steps:
        raise ValueError(
            f"Expected {args.case_label} steps {expected_steps}, got {actual_steps}"
        )
    dt = float(config["analysis"]["dt"])
    times = np.asarray([particle_step(path) * dt for path in files])
    pressure_rows: list[np.ndarray] = []
    points_by_step: dict[int, np.ndarray] = {}
    for path in files:
        points, pressure = align_liquid_pressure(path, expected_ids)
        pressure_rows.append(pressure - checkpoint_pressure)
        step = particle_step(path)
        if any(abs(step * dt - target) < 0.5 * dt for target in SNAPSHOT_TIMES_S):
            points_by_step[step] = points
    pressure_excess = np.asarray(pressure_rows)
    fit = (times >= FIT_START_S - 1.0e-12) & (times <= FIT_END_S + 1.0e-12)
    _, amplitude, phase = harmonic_fields(times[fit], pressure_excess[fit], WAVE_PERIOD_S)
    harmonic_r2 = harmonic_r_squared(
        times[fit], pressure_excess[fit], WAVE_PERIOD_S
    )
    surface_amplitude = same_column_surface_values(initial_coordinates, amplitude)
    surface_phase = same_column_surface_values(initial_coordinates, phase)
    amplitude_ratio = amplitude / np.maximum(surface_amplitude, np.finfo(float).eps)
    phase_difference = wrap_degrees(phase - surface_phase)
    phase_difference[
        (amplitude_ratio < MINIMUM_LOCAL_AMPLITUDE_RATIO)
        | (harmonic_r2 < MINIMUM_HARMONIC_R_SQUARED)
    ] = np.nan

    history = read_pipeline_history(history_path)
    initial_center = np.asarray(
        config["analysis"]["rigid_pipeline"]["center"], dtype=np.float64
    )
    snapshots: list[dict[str, Any]] = []
    for time_s in SNAPSHOT_TIMES_S:
        index = int(np.argmin(np.abs(times - time_s)))
        if abs(times[index] - time_s) > 0.5 * dt:
            raise ValueError(f"No exact saved frame for t={time_s}")
        step = particle_step(files[index])
        snapshots.append(
            {
                "time_s": float(times[index]),
                "step": step,
                "phase_fraction": ((times[index] - SNAPSHOT_TIMES_S[0]) / WAVE_PERIOD_S) % 1.0,
                "points": points_by_step[step],
                "pressure_excess_pa": pressure_excess[index],
                "pipeline_center": pipeline_pose(
                    history, float(times[index]), initial_center
                ),
                "path": str(files[index].resolve()),
            }
        )

    plot_fields(
        output / f"{args.case_label}_liquid_pressure_2d_phase_lag",
        config,
        snapshots,
        initial_coordinates,
        amplitude_ratio,
        phase_difference,
        args.case_label,
        parameters,
    )
    time_depth = plot_time_depth(
        output / f"{args.case_label}_liquid_pressure_time_depth",
        times,
        pressure_excess,
        initial_coordinates,
        amplitude,
        phase,
        harmonic_r2,
        args.case_label,
        parameters,
    )
    figure_paths = tuple(
        path.resolve()
        for base in (
            output / f"{args.case_label}_liquid_pressure_2d_phase_lag",
            output / f"{args.case_label}_liquid_pressure_time_depth",
        )
        for path in figure_output_paths(base)
    )
    phase_inputs = {
        "schema": INPUT_SCHEMA,
        "analysis_config": artifact_record(config_path),
        "checkpoint_config": artifact_record(checkpoint_config),
        "checkpoint": artifact_record(checkpoint),
        "particles": artifact_record(particle_file),
        "pipeline_history": artifact_record(history_path),
        "result_directory": str(result.resolve()),
        "result_particle_vtp": [artifact_record(path) for path in files],
    }
    artifacts = {
        "schema": ARTIFACT_SCHEMA,
        "analyzer": artifact_record(ANALYZER_PATH),
        "plot_script": artifact_record(PLOT_SCRIPT_PATH),
        "figures": [artifact_record(path) for path in figure_paths],
        "inputs": phase_inputs,
    }

    probes = build_probes(config, initial_coordinates)
    surface_index = next(probe.index for probe in probes if probe.name == "surface")
    probe_audit: dict[str, Any] = {}
    for probe in probes:
        difference = float(wrap_degrees(np.asarray([phase[probe.index] - phase[surface_index]]))[0])
        same_column_difference = float(phase_difference[probe.index])
        probe_audit[probe.name] = {
            "particle_id": probe.index,
            "x_m": probe.actual_x,
            "y_m": probe.actual_y,
            "fundamental_amplitude_pa": float(amplitude[probe.index]),
            "harmonic_r_squared": float(harmonic_r2[probe.index]),
            "phase_difference_from_registered_surface_deg": difference,
            "phase_difference_from_same_column_surface_deg": (
                same_column_difference
                if math.isfinite(same_column_difference)
                else None
            ),
        }

    audit = {
        "schema": SCHEMA,
        "case": args.case_label,
        "parameters": parameters.as_dict(),
        "parameter_subtitle": parameter_subtitle(parameters),
        "require_unsmoothed": args.require_unsmoothed,
        "artifacts": artifacts,
        "field": "PIC_liquid_pressures",
        "excess_definition": (
            "current PIC_liquid_pressures minus the supplied equilibrium "
            "checkpoint value for the same particle ID"
        ),
        "particle_count": int(initial_coordinates.shape[0]),
        "frame_count": len(files),
        "wave_period_s": WAVE_PERIOD_S,
        "fit_window_s": [FIT_START_S, FIT_END_S],
        "minimum_local_amplitude_ratio_for_phase": MINIMUM_LOCAL_AMPLITUDE_RATIO,
        "minimum_harmonic_r_squared_for_phase": MINIMUM_HARMONIC_R_SQUARED,
        "checkpoint": phase_inputs["checkpoint"],
        "config": phase_inputs["analysis_config"],
        "snapshots": [
            {
                "time_s": item["time_s"],
                "step": item["step"],
                "phase_fraction": item["phase_fraction"],
                "path": item["path"],
                "sha256": sha256(Path(item["path"])),
            }
            for item in snapshots
        ],
        "probe_results": probe_audit,
        "time_depth": time_depth,
        "phase_present": bool(
            any(
                item["phase_difference_from_same_column_surface_deg"] is not None
                and abs(item["phase_difference_from_same_column_surface_deg"]) >= 5.0
                for item in probe_audit.values()
            )
        ),
        "interpretation": (
            f"Phase differences are directly resolved in the {args.case_label} "
            "liquid-pressure "
            "fundamental. The same-column surface difference isolates vertical/"
            "subsurface lag; the registered central-surface difference also contains "
            "the horizontal travelling-wave phase. Sign follows the harmonic delay "
            "convention; use absolute magnitude when describing the presence of lag."
        ),
    }
    audit_path = output / f"{args.case_label}_liquid_pressure_phase_lag.audit.json"
    validate_phase_audit_artifacts(audit, expected_figure_paths=figure_paths)
    atomic_write_json(audit_path, audit)
    validate_phase_audit_artifacts(
        read_json(audit_path), expected_figure_paths=figure_paths
    )
    print(json.dumps(audit["probe_results"], indent=2))
    print(f"Wrote liquid-pressure phase-lag figures to {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
