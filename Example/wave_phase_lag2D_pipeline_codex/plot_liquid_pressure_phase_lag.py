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
import re
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle

from analyze_study import build_probes, read_vtp
from plot_manuscript_figures import _masked_triangulation


SCHEMA = "pipeline-liquid-pressure-phase-lag-audit-v1"
CASE = "HS"
WAVE_PERIOD_S = 1.3
FIT_START_S = 1.3
FIT_END_S = 3.9
SNAPSHOT_TIMES_S = (2.6, 2.925, 3.25, 3.575)
MINIMUM_LOCAL_AMPLITUDE_RATIO = 0.05
FARFIELD_X_M = 1.06


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


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
    ids = np.rint(np.asarray(arrays["ids"], dtype=np.float64)).astype(np.int64)
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


def pipeline_pose(history: list[dict[str, float]], time_s: float) -> np.ndarray:
    row = min(history, key=lambda item: abs(item["time"] - time_s))
    return np.asarray([0.7 + row["displacement_x"], 0.41 + row["displacement_y"]])


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
            f"{probe.name} {value:+.1f}°",
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
        "HS two-dimensional PIC liquid-pressure response\n"
        "common snapshot scale; harmonic fit 1.3–3.9 s; phase masked where local amplitude <5% of surface",
        fontsize=13,
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output.with_suffix(".png"), dpi=300, bbox_inches="tight")
    figure.savefig(output.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(figure)


def plot_time_depth(
    output: Path,
    times: np.ndarray,
    pressure_excess: np.ndarray,
    initial_coordinates: np.ndarray,
    amplitude: np.ndarray,
    phase: np.ndarray,
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
    phase_difference[amplitude_ratio < MINIMUM_LOCAL_AMPLITUDE_RATIO] = np.nan

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
        "HS liquid-pressure time–depth propagation (registered pre-release window)",
        fontsize=13,
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output.with_suffix(".png"), dpi=300, bbox_inches="tight")
    figure.savefig(output.with_suffix(".pdf"), bbox_inches="tight")
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
    args = parser.parse_args()

    root = Path(__file__).resolve().parent
    config_path = root / "configs/screen/02_HS.json"
    config = read_json(config_path)
    result = root / "results/screen/PLP_SCREEN_HS_SANI"
    checkpoint = root / "results/screen/PLP_SCREEN_HS_EQ/particle40000.vtp"
    particle_file = root / "particles.txt"
    history_paths = sorted(result.glob("pipeline-history*.csv"))
    if len(history_paths) != 1:
        raise ValueError(
            f"Expected exactly one HS pipeline history, found {len(history_paths)}"
        )
    history_path = history_paths[0]
    output = (root / args.output_dir).resolve() if not args.output_dir.is_absolute() else args.output_dir

    initial_coordinates = read_initial_coordinates(particle_file)
    _, checkpoint_arrays = read_vtp(checkpoint)
    checkpoint_ids = np.rint(
        np.asarray(checkpoint_arrays["ids"], dtype=np.float64)
    ).astype(np.int64)
    checkpoint_order = np.argsort(checkpoint_ids)
    expected_ids = checkpoint_ids[checkpoint_order]
    if not np.array_equal(expected_ids, np.arange(initial_coordinates.shape[0])):
        raise ValueError("Checkpoint particle IDs are not contiguous 0..N-1")
    checkpoint_pressure = np.asarray(
        checkpoint_arrays["PIC_liquid_pressures"], dtype=np.float64
    )[checkpoint_order]

    files = sorted(result.glob("particle*.vtp"), key=particle_step)
    if len(files) != 160:
        raise ValueError(f"Expected 160 HS frames, found {len(files)}")
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
    surface_amplitude = same_column_surface_values(initial_coordinates, amplitude)
    surface_phase = same_column_surface_values(initial_coordinates, phase)
    amplitude_ratio = amplitude / np.maximum(surface_amplitude, np.finfo(float).eps)
    phase_difference = wrap_degrees(phase - surface_phase)
    phase_difference[amplitude_ratio < MINIMUM_LOCAL_AMPLITUDE_RATIO] = np.nan

    history = read_pipeline_history(history_path)
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
                "pipeline_center": pipeline_pose(history, float(times[index])),
                "path": str(files[index].resolve()),
            }
        )

    plot_fields(
        output / "HS_liquid_pressure_2d_phase_lag",
        config,
        snapshots,
        initial_coordinates,
        amplitude_ratio,
        phase_difference,
    )
    time_depth = plot_time_depth(
        output / "HS_liquid_pressure_time_depth",
        times,
        pressure_excess,
        initial_coordinates,
        amplitude,
        phase,
    )

    probes = build_probes(config, initial_coordinates)
    surface_index = next(probe.index for probe in probes if probe.name == "surface")
    probe_audit: dict[str, Any] = {}
    for probe in probes:
        difference = float(wrap_degrees(np.asarray([phase[probe.index] - phase[surface_index]]))[0])
        probe_audit[probe.name] = {
            "particle_id": probe.index,
            "x_m": probe.actual_x,
            "y_m": probe.actual_y,
            "fundamental_amplitude_pa": float(amplitude[probe.index]),
            "phase_difference_from_registered_surface_deg": difference,
            "phase_difference_from_same_column_surface_deg": float(
                phase_difference[probe.index]
            ),
        }

    audit = {
        "schema": SCHEMA,
        "case": CASE,
        "field": "PIC_liquid_pressures",
        "excess_definition": (
            "current PIC_liquid_pressures minus audited HS_EQ checkpoint value "
            "for the same particle ID"
        ),
        "particle_count": int(initial_coordinates.shape[0]),
        "frame_count": len(files),
        "wave_period_s": WAVE_PERIOD_S,
        "fit_window_s": [FIT_START_S, FIT_END_S],
        "minimum_local_amplitude_ratio_for_phase": MINIMUM_LOCAL_AMPLITUDE_RATIO,
        "checkpoint": {
            "path": str(checkpoint.resolve()),
            "sha256": sha256(checkpoint),
        },
        "config": {"path": str(config_path.resolve()), "sha256": sha256(config_path)},
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
            max(
                abs(item["phase_difference_from_same_column_surface_deg"])
                for item in probe_audit.values()
            )
            >= 5.0
        ),
        "interpretation": (
            "Phase differences are directly resolved in the HS liquid-pressure "
            "fundamental. The same-column surface difference isolates vertical/"
            "subsurface lag; the registered central-surface difference also contains "
            "the horizontal travelling-wave phase. Sign follows the harmonic delay "
            "convention; use absolute magnitude when describing the presence of lag."
        ),
    }
    (output / "HS_liquid_pressure_phase_lag.audit.json").write_text(
        json.dumps(audit, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(audit["probe_results"], indent=2))
    print(f"Wrote liquid-pressure phase-lag figures to {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
