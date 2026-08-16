#!/usr/bin/env python3
"""Plot manuscript-ready curve panels and two-dimensional field comparisons.

The script consumes only ``analyze_study.py`` outputs.  It creates three
curve-based composite figures corresponding to the Section 7.5 draft, renders
LS/HS and RL/RE particle fields at registered common physical times, and writes
a JSON manifest recording every selected VTP frame and colour limit.  It never
chooses separate best-looking times per case.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.tri as mtri  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.patches import Circle, Rectangle  # noqa: E402

from analyze_study import (
    MINIMUM_REFERENCE_VERTICAL_STRESS_PA,
    read_pipeline_history,
    read_vtp,
    upward_excess_seepage_force_index,
)


IF_AREA = "area_support_ROI_upward_seepage_IF_ge_1_m2"
STRESS_AREA = "area_support_ROI_stress_loss_Rsigma_le_0p05_m2"
JOINT_AREA = "area_support_ROI_joint_IF_ge_1_and_Rsigma_le_0p05_m2"
PARTICLE_PATTERN = re.compile(r"particle(\d+)\.vtp$")
CASE_PREFIX = {
    "LS": "02_LS",
    "HS": "02_HS",
    "HM": "03_HM",
    "RL": "04_RL",
    "RM": "04_RM",
    "RE": "04_RE",
}
COLORS = {
    "LS": "#0072B2",
    "HS": "#D55E00",
    "HM": "#CC79A7",
    "RL": "#009E73",
    "RM": "#E69F00",
    "RE": "#56B4E9",
}


def read_csv(path: Path) -> list[dict[str, float]]:
    with path.open(newline="", encoding="utf-8") as stream:
        rows = [
            {key: float(value) for key, value in row.items() if value != ""}
            for row in csv.DictReader(stream)
        ]
    if not rows:
        raise ValueError(f"No data rows in {path}")
    return rows


def read_summary(directory: Path, code: str) -> dict[str, Any]:
    path = directory / f"{CASE_PREFIX[code]}_summary.json"
    if not path.is_file():
        raise FileNotFoundError(path)
    return json.loads(path.read_text(encoding="utf-8"))


def probe_rows(directory: Path, code: str) -> list[dict[str, float]]:
    return read_csv(directory / f"{CASE_PREFIX[code]}_probe_history.csv")


def cycle_rows(directory: Path, code: str) -> list[dict[str, float]]:
    return read_csv(directory / f"{CASE_PREFIX[code]}_pipeline_cycles.csv")


def series(rows: list[dict[str, float]], key: str) -> tuple[list[float], list[float]]:
    selected = [(row["time_s"], row[key]) for row in rows if key in row]
    if not selected:
        raise KeyError(f"No {key} values in history")
    return [item[0] for item in selected], [item[1] for item in selected]


def style_axis(axis: plt.Axes, label: str) -> None:
    axis.text(
        0.0,
        1.04,
        label,
        transform=axis.transAxes,
        fontsize=10,
        fontweight="bold",
        va="bottom",
    )
    axis.grid(True, alpha=0.25, linewidth=0.6)


def save_figure(figure: plt.Figure, output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.tight_layout()
    figure.savefig(output.with_suffix(".png"), dpi=300, bbox_inches="tight")
    figure.savefig(output.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(figure)


def plot_hydraulic_bridge(directory: Path, output: Path) -> None:
    histories = {code: probe_rows(directory, code) for code in ("LS", "HS")}
    summaries = {code: read_summary(directory, code) for code in ("LS", "HS")}
    figure, axes = plt.subplots(2, 2, figsize=(10.2, 6.8))

    for code, rows in histories.items():
        for probe, linestyle in (("surface", "--"), ("crown", "-"), ("invert", ":")):
            time, pressure = series(rows, f"{probe}_pressure_excess_pa")
            axes[0, 0].plot(
                time,
                pressure,
                color=COLORS[code],
                linestyle=linestyle,
                linewidth=1.1,
                label=f"{code} {probe}",
            )
    axes[0, 0].set_xlabel("Time (s)")
    axes[0, 0].set_ylabel("Wave-induced excess mixture pressure (Pa)")
    axes[0, 0].legend(ncol=2, fontsize=7)
    style_axis(axes[0, 0], "(a)")

    probes = ("crown", "shoulder", "invert")
    positions = list(range(len(probes)))
    width = 0.36
    for offset, code in ((-width / 2, "LS"), (width / 2, "HS")):
        values = [
            abs(summaries[code]["phase"]["probes"][probe]["phase_lag_deg"])
            for probe in probes
        ]
        axes[0, 1].bar(
            [position + offset for position in positions],
            values,
            width,
            color=COLORS[code],
            label=code,
        )
    axes[0, 1].set_xticks(positions, probes)
    axes[0, 1].set_ylabel("Absolute phase lag (degrees)")
    axes[0, 1].legend()
    style_axis(axes[0, 1], "(b)")

    for offset, code in ((-width / 2, "LS"), (width / 2, "HS")):
        values = [
            summaries[code]["phase"]["probes"][probe][
                "amplitude_ratio_to_surface"
            ]
            for probe in probes
        ]
        axes[1, 0].bar(
            [position + offset for position in positions],
            values,
            width,
            color=COLORS[code],
            label=code,
        )
    axes[1, 0].set_xticks(positions, probes)
    axes[1, 0].set_ylabel("Amplitude / surface amplitude")
    axes[1, 0].set_ylim(bottom=0.0)
    axes[1, 0].legend()
    style_axis(axes[1, 0], "(c)")

    for code in ("LS", "HS"):
        rows = cycle_rows(directory, code)
        axes[1, 1].plot(
            [row["cycle_after_release"] for row in rows],
            [row["uplift_over_D"] for row in rows],
            marker="o",
            markersize=3,
            color=COLORS[code],
            label=code,
        )
    axes[1, 1].axhline(0.5, color="0.3", linestyle="--", linewidth=0.9)
    axes[1, 1].set_xlabel("Cycle after release")
    axes[1, 1].set_ylabel("Pipeline uplift, $u_y/D$")
    axes[1, 1].legend()
    style_axis(axes[1, 1], "(d)")
    save_figure(figure, output)


def plot_phase_control(directory: Path, output: Path) -> None:
    histories = {code: probe_rows(directory, code) for code in ("RL", "RE")}
    figure, axes = plt.subplots(2, 2, figsize=(10.2, 6.8))
    for code, rows in histories.items():
        time, pressure = series(rows, "invert_pressure_excess_pa")
        axes[0, 0].plot(time, pressure, color=COLORS[code], label=code)
        time, area = series(rows, IF_AREA)
        axes[0, 1].plot(time, area, color=COLORS[code], label=code)
        time, joint_area = series(rows, JOINT_AREA)
        axes[0, 1].plot(
            time,
            joint_area,
            color=COLORS[code],
            linestyle="--",
            linewidth=1.0,
            label=f"{code} joint IF + $R_\\sigma$",
        )
        time, area = series(rows, STRESS_AREA)
        axes[1, 0].plot(time, area, color=COLORS[code], label=code)
    axes[0, 0].set_xlabel("Time (s)")
    axes[0, 0].set_ylabel("Invert wave-induced excess mixture pressure (Pa)")
    axes[0, 0].legend()
    style_axis(axes[0, 0], "(a)")
    axes[0, 1].set_xlabel("Time (s)")
    axes[0, 1].set_ylabel("Support-zone IF / joint area (m$^2$)")
    axes[0, 1].legend()
    style_axis(axes[0, 1], "(b)")
    axes[1, 0].set_xlabel("Time (s)")
    axes[1, 0].set_ylabel("Support-zone area with $R_\\sigma\\leq0.05$ (m$^2$)")
    axes[1, 0].legend()
    style_axis(axes[1, 0], "(c)")
    for code in ("RL", "RE"):
        rows = cycle_rows(directory, code)
        axes[1, 1].plot(
            [row["cycle_after_release"] for row in rows],
            [row["uplift_over_D"] for row in rows],
            marker="o",
            markersize=3,
            color=COLORS[code],
            label=code,
        )
    axes[1, 1].set_xlabel("Cycle after release")
    axes[1, 1].set_ylabel("Pipeline uplift, $u_y/D$")
    axes[1, 1].legend()
    style_axis(axes[1, 1], "(d)")
    save_figure(figure, output)


def plot_constitutive(directory: Path, output: Path) -> None:
    histories = {code: probe_rows(directory, code) for code in ("RL", "RM")}
    figure, axes = plt.subplots(2, 2, figsize=(10.2, 6.8))
    for code, rows in histories.items():
        p = [row["invert_p_effective_pa"] / 1000.0 for row in rows]
        q = [row["invert_q_pa"] / 1000.0 for row in rows]
        axes[0, 0].plot(p, q, color=COLORS[code], label=code, linewidth=1.0)
        state_key = "invert_eps_p_q" if code == "RL" else "invert_pdstrain"
        time, state = series(rows, state_key)
        axes[0, 1].plot(time, state, color=COLORS[code], label=code)
        time, area = series(rows, STRESS_AREA)
        axes[1, 0].plot(time, area, color=COLORS[code], label=code)
    axes[0, 0].set_xlabel("Mean effective stress, $p'$ (kPa)")
    axes[0, 0].set_ylabel("Deviatoric stress, $q$ (kPa)")
    axes[0, 0].legend()
    style_axis(axes[0, 0], "(a)")
    axes[0, 1].set_xlabel("Time (s)")
    axes[0, 1].set_ylabel("Plastic-history variable")
    axes[0, 1].legend(title="RL: $\\epsilon^p_q$; RM: pdstrain", fontsize=8)
    style_axis(axes[0, 1], "(b)")
    axes[1, 0].set_xlabel("Time (s)")
    axes[1, 0].set_ylabel("Support-zone area with $R_\\sigma\\leq0.05$ (m$^2$)")
    axes[1, 0].legend()
    style_axis(axes[1, 0], "(c)")
    for code in ("RL", "RM"):
        rows = cycle_rows(directory, code)
        axes[1, 1].plot(
            [row["cycle_after_release"] for row in rows],
            [row["uplift_over_D"] for row in rows],
            marker="o",
            markersize=3,
            color=COLORS[code],
            label=code,
        )
    axes[1, 1].set_xlabel("Cycle after release")
    axes[1, 1].set_ylabel("Pipeline uplift, $u_y/D$")
    axes[1, 1].legend()
    style_axis(axes[1, 1], "(d)")
    save_figure(figure, output)


def particle_step(path: Path) -> int:
    match = PARTICLE_PATTERN.match(path.name)
    if match is None:
        raise ValueError(f"Cannot parse step from {path}")
    return int(match.group(1))


def select_common_field_frames(
    directory: Path,
    *,
    codes: tuple[str, str],
    driver: str,
    metric: str,
) -> dict[str, Any]:
    driver_rows = probe_rows(directory, driver)
    if metric not in driver_rows[0] and not any(metric in row for row in driver_rows):
        raise KeyError(f"{driver} history has no {metric}")
    critical = max(
        (row for row in driver_rows if metric in row),
        key=lambda row: row[metric],
    )
    target_time = critical["time_s"]
    output: dict[str, Any] = {
        "selection_rule": (
            f"common physical time at the maximum {driver} {metric}; "
            f"the same target time is used for {' and '.join(codes)}"
        ),
        "target_time_s": target_time,
        "driver": driver,
        "selection_metric": metric,
        "driver_metric_value": critical[metric],
        "cases": {},
    }
    for code in codes:
        summary = read_summary(directory, code)
        dt = float(summary.get("analysis_dt_s", 0.0))
        result_directory = Path(summary["result_directory"])
        paths = sorted(result_directory.glob("particle*.vtp"), key=particle_step)
        if not paths:
            raise FileNotFoundError(f"No VTP frames in {result_directory}")
        if dt <= 0.0:
            config = Path(summary["config"])
            dt = float(json.loads(config.read_text(encoding="utf-8"))["analysis"]["dt"])
        target_step = int(round(target_time / dt))
        chosen = min(paths, key=lambda path: abs(particle_step(path) - target_step))
        chosen_step = particle_step(chosen)
        chosen_time = chosen_step * dt
        if abs(chosen_time - target_time) > 0.5 * dt + 1.0e-12:
            raise ValueError(
                f"{code} has no particle frame on the registered target time "
                f"{target_time}; nearest is {chosen_time}"
            )
        output["cases"][code] = {
            "path": str(chosen.resolve()),
            "step": chosen_step,
            "time_s": chosen_step * dt,
            "time_offset_from_target_s": chosen_step * dt - target_time,
        }
    selected_times = [item["time_s"] for item in output["cases"].values()]
    if not all(
        math.isclose(value, selected_times[0], rel_tol=0.0, abs_tol=1.0e-12)
        for value in selected_times[1:]
    ):
        raise ValueError(
            f"{' and '.join(codes)} do not share a common saved physical time"
        )
    output["selected_common_time_s"] = selected_times[0]
    return output


def common_field_frames(directory: Path) -> dict[str, Any]:
    """Return the pre-registered common RL/RE field-frame selection."""

    output = select_common_field_frames(
        directory,
        codes=("RL", "RE"),
        driver="RL",
        metric=IF_AREA,
    )
    # Preserve the established manifest key for downstream consumers.
    output["RL_support_IF_area_m2"] = output["driver_metric_value"]
    return output


def registered_phase_field_frames(directory: Path) -> dict[str, dict[str, Any]]:
    """Register the three a-priori causal field times on RL and reuse for RE."""

    return {
        "t_IF": select_common_field_frames(
            directory,
            codes=("RL", "RE"),
            driver="RL",
            metric=IF_AREA,
        ),
        "t_joint": select_common_field_frames(
            directory,
            codes=("RL", "RE"),
            driver="RL",
            metric=JOINT_AREA,
        ),
        "t_Rsigma": select_common_field_frames(
            directory,
            codes=("RL", "RE"),
            driver="RL",
            metric=STRESS_AREA,
        ),
    }


def _aligned_difference(
    left_ids: np.ndarray,
    left_values: np.ndarray,
    right_ids: np.ndarray,
    right_values: np.ndarray,
) -> np.ndarray:
    """Return left-minus-right values in left particle-ID order."""

    left_ids = np.asarray(left_ids, dtype=np.int64).reshape(-1)
    right_ids = np.asarray(right_ids, dtype=np.int64).reshape(-1)
    left_values = np.asarray(left_values, dtype=np.float64).reshape(-1)
    right_values = np.asarray(right_values, dtype=np.float64).reshape(-1)
    if left_ids.size != left_values.size or right_ids.size != right_values.size:
        raise ValueError("Particle IDs and field values have inconsistent counts")
    if len(np.unique(left_ids)) != left_ids.size or len(np.unique(right_ids)) != right_ids.size:
        raise ValueError("Particle IDs must be unique for a field difference")
    order = np.argsort(right_ids)
    sorted_ids = right_ids[order]
    locations = np.searchsorted(sorted_ids, left_ids)
    if np.any(locations >= sorted_ids.size) or not np.array_equal(
        sorted_ids[locations], left_ids
    ):
        raise ValueError("Compared VTP frames do not contain the same particle IDs")
    return left_values - right_values[order[locations]]


def _masked_triangulation(
    points: np.ndarray,
    config: dict[str, Any],
    pipeline_center: np.ndarray | None = None,
) -> mtri.Triangulation:
    x = np.asarray(points[:, 0], dtype=np.float64)
    y = np.asarray(points[:, 1], dtype=np.float64)
    triangulation = mtri.Triangulation(x, y)
    triangles = triangulation.triangles
    triangle_points = points[triangles, :2]
    edges = np.stack(
        (
            triangle_points[:, 1] - triangle_points[:, 0],
            triangle_points[:, 2] - triangle_points[:, 1],
            triangle_points[:, 0] - triangle_points[:, 2],
        ),
        axis=1,
    )
    maximum_edge = np.max(np.linalg.norm(edges, axis=2), axis=1)
    cell_size = float(config["mesh"]["cellsize_min"])
    pipeline = config["analysis"]["rigid_pipeline"]
    center = (
        np.asarray(pipeline["center"], dtype=np.float64)
        if pipeline_center is None
        else np.asarray(pipeline_center, dtype=np.float64)
    )
    radius = float(pipeline["outer_radius"])
    centroids = np.mean(triangle_points, axis=1)
    inside_pipe = np.linalg.norm(centroids - center, axis=1) < radius
    triangulation.set_mask((maximum_edge > 1.5 * cell_size) | inside_pipe)
    return triangulation


def _pipeline_pose(
    summary: dict[str, Any], config: dict[str, Any], target_time: float
) -> dict[str, Any]:
    """Return the saved rigid-pipeline pose nearest the selected VTP time."""

    history = read_pipeline_history(Path(summary["result_directory"]))
    row = min(history, key=lambda item: abs(item["time"] - target_time))
    initial = np.asarray(
        config["analysis"]["rigid_pipeline"]["center"], dtype=np.float64
    )
    center = initial + np.asarray(
        [row["displacement_x"], row["displacement_y"]], dtype=np.float64
    )
    return {
        "center": center,
        "angle_rad": float(row["angle"]),
        "history_time_s": float(row["time"]),
        "history_time_offset_s": float(row["time"] - target_time),
        "source": "nearest saved rigid-pipeline history row",
    }


def _eligible_stress_ratio(
    ratios: np.ndarray,
    initial_vertical_stress: np.ndarray,
    minimum_reference_pa: float = MINIMUM_REFERENCE_VERTICAL_STRESS_PA,
) -> np.ndarray:
    """Mask undefined R_sigma values before plotting or differencing."""

    ratios = np.asarray(ratios, dtype=np.float64).reshape(-1).copy()
    reference = np.asarray(initial_vertical_stress, dtype=np.float64).reshape(-1)
    if ratios.shape != reference.shape:
        raise ValueError("R_sigma and initial-stress arrays must have equal counts")
    ratios[np.abs(reference) < minimum_reference_pa] = np.nan
    return ratios


def _case_root(config_path: Path) -> Path:
    root = config_path.resolve().parent
    while not (root / "prepare_study.py").exists() and root != root.parent:
        root = root.parent
    if not (root / "prepare_study.py").exists():
        raise FileNotFoundError(f"Cannot locate case root above {config_path}")
    return root


def _checkpoint_pressure_excess(
    summary: dict[str, Any],
    current_ids: np.ndarray,
    current_pressure: np.ndarray,
) -> tuple[np.ndarray, dict[str, Any]]:
    """Subtract the audited resume-checkpoint pressure in particle-ID space."""

    runtime = summary.get("completion_audit", {}).get("runtime_dependencies", {})
    resume = runtime.get("resume_equilibrium")
    if not isinstance(resume, dict) or not isinstance(resume.get("qa_vtp"), dict):
        raise ValueError(
            f"{summary['case']} has neither PIC_pore_pressure_excess nor an "
            "audited resume-equilibrium QA VTP"
        )
    record = resume["qa_vtp"]
    reference_path = Path(record["path"])
    if not reference_path.is_absolute():
        reference_path = _case_root(Path(summary["config"])) / reference_path
    reference_path = reference_path.resolve()
    _, reference_arrays = read_vtp(reference_path)
    missing = [
        name for name in ("ids", "PIC_pore_pressures") if name not in reference_arrays
    ]
    if missing:
        raise KeyError(f"Resume-checkpoint VTP is missing {missing}: {reference_path}")
    excess = _aligned_difference(
        current_ids,
        current_pressure,
        reference_arrays["ids"],
        reference_arrays["PIC_pore_pressures"],
    )
    return excess, {
        "method": "current_minus_audited_resume_checkpoint_by_particle_id",
        "checkpoint_vtp": str(reference_path),
        "checkpoint_vtp_sha256": record.get("sha256"),
    }


def _field_frame(
    directory: Path, code: str, frame: dict[str, Any]
) -> tuple[
    dict[str, Any],
    np.ndarray,
    dict[str, np.ndarray],
    dict[str, np.ndarray],
    dict[str, Any],
]:
    summary = read_summary(directory, code)
    config = json.loads(Path(summary["config"]).read_text(encoding="utf-8"))
    points, arrays = read_vtp(Path(frame["path"]))
    required = (
        "ids",
        "liquid_seepage_forces",
        "liquid_densities",
        "gamma_sub",
        "initial_vertical_effective_stresses",
        "vertical_effective_stress_remaining_ratios",
    )
    missing = [name for name in required if name not in arrays]
    if missing:
        raise KeyError(f"{code} field frame is missing {missing}")
    if "PIC_pore_pressure_excess" in arrays:
        pressure_excess = np.asarray(
            arrays["PIC_pore_pressure_excess"], dtype=float
        )
        pressure_audit = {
            "method": "solver_saved_current_minus_checkpoint_initial_pressure"
        }
    elif "PIC_pore_pressures" in arrays:
        pressure_excess, pressure_audit = _checkpoint_pressure_excess(
            summary,
            arrays["ids"],
            arrays["PIC_pore_pressures"],
        )
    else:
        raise KeyError(
            f"{code} field frame contains neither PIC_pore_pressure_excess nor "
            "PIC_pore_pressures"
        )
    ratios = _eligible_stress_ratio(
        arrays["vertical_effective_stress_remaining_ratios"],
        arrays["initial_vertical_effective_stresses"],
    )
    seepage = upward_excess_seepage_force_index(
        arrays["liquid_seepage_forces"],
        arrays["liquid_densities"],
        arrays["gamma_sub"],
        np.asarray(config["external_loading_conditions"]["gravity"], dtype=float),
    )
    fields = {
        "wave_excess_mixture_pressure_kpa": np.asarray(
            pressure_excess, dtype=float
        ) / 1000.0,
        "upward_seepage_IF": seepage,
        "vertical_stress_remaining_ratio": ratios,
        "joint_liquefaction_indicator": (
            (seepage >= 1.0) & np.isfinite(ratios) & (ratios <= 0.05)
        ).astype(float),
    }
    pose = _pipeline_pose(summary, config, float(frame["time_s"]))
    pose["pressure_reference_audit"] = pressure_audit
    return config, points, arrays, fields, pose


def _robust_limits(
    values: list[np.ndarray], *, include: float | None = None
) -> tuple[float, float]:
    finite = np.concatenate([np.asarray(value)[np.isfinite(value)] for value in values])
    if finite.size == 0:
        raise ValueError("Cannot plot a field without finite values")
    lower, upper = np.percentile(finite, [1.0, 99.0])
    if include is not None:
        lower = min(lower, include)
        upper = max(upper, include)
    if math.isclose(float(lower), float(upper), rel_tol=0.0, abs_tol=1.0e-15):
        padding = max(abs(float(lower)) * 0.05, 1.0e-12)
        lower -= padding
        upper += padding
    return float(lower), float(upper)


def _draw_field(
    axis: plt.Axes,
    points: np.ndarray,
    values: np.ndarray,
    config: dict[str, Any],
    *,
    title: str,
    cmap: str,
    limits: tuple[float, float],
    threshold: float | None = None,
    pipeline_pose: dict[str, Any] | None = None,
    spatial_limits: tuple[float, float, float, float] | None = None,
    draw_pipeline: bool = True,
) -> Any:
    pipeline_center = (
        np.asarray(pipeline_pose["center"], dtype=float)
        if pipeline_pose is not None
        else np.asarray(config["analysis"]["rigid_pipeline"]["center"], dtype=float)
    )
    triangulation = _masked_triangulation(points, config, pipeline_center)
    finite = np.isfinite(values)
    if not np.all(finite):
        invalid_triangles = np.any(~finite[triangulation.triangles], axis=1)
        existing = triangulation.mask
        triangulation.set_mask(
            invalid_triangles if existing is None else (existing | invalid_triangles)
        )
    plot_values = np.where(finite, values, 0.0)
    artist = axis.tripcolor(
        triangulation,
        plot_values,
        shading="gouraud",
        cmap=cmap,
        vmin=limits[0],
        vmax=limits[1],
        rasterized=True,
    )
    finite_values = np.asarray(values)[finite]
    if (
        threshold is not None
        and finite_values.size
        and float(np.min(finite_values)) < threshold
        and float(np.max(finite_values)) > threshold
    ):
        axis.tricontour(
            triangulation,
            plot_values,
            levels=[threshold],
            colors="black",
            linewidths=0.8,
        )
    pipeline = config["analysis"]["rigid_pipeline"]
    center_x, center_y = map(float, pipeline_center)
    radius = float(pipeline["outer_radius"])
    diameter = 2.0 * radius
    fluid = config["materials"][1]
    bed_y = float(fluid["sea_level"]) - float(fluid["depth_left"])
    if draw_pipeline:
        axis.add_patch(
            Circle((center_x, center_y), radius, fc="white", ec="black", lw=0.9)
        )
    initial_center_x = float(pipeline["center"][0])
    axis.add_patch(
        Rectangle(
            (initial_center_x - 1.5 * diameter, bed_y - 2.0 * diameter),
            3.0 * diameter,
            2.0 * diameter,
            fill=False,
            edgecolor="#E69F00",
            linestyle="--",
            linewidth=0.9,
        )
    )
    if spatial_limits is None:
        spatial_limits = (
            float(np.min(points[:, 0])),
            float(np.max(points[:, 0])),
            float(np.min(points[:, 1])),
            float(np.max(points[:, 1])) + 0.005,
        )
    axis.set_xlim(spatial_limits[0], spatial_limits[1])
    axis.set_ylim(spatial_limits[2], spatial_limits[3])
    axis.set_aspect("equal")
    axis.set_title(title, fontsize=9)
    axis.set_xlabel("x (m)")
    axis.set_ylabel("y (m)")
    return artist


def _plot_field_triptych(
    axes: np.ndarray,
    *,
    left_code: str,
    right_code: str,
    left: tuple[
        dict[str, Any],
        np.ndarray,
        dict[str, np.ndarray],
        dict[str, np.ndarray],
        dict[str, Any],
    ],
    right: tuple[
        dict[str, Any],
        np.ndarray,
        dict[str, np.ndarray],
        dict[str, np.ndarray],
        dict[str, Any],
    ],
    field: str,
    label: str,
    cmap: str,
    difference_label: str,
    difference_sign: float = 1.0,
    threshold: float | None = None,
) -> dict[str, Any]:
    left_config, left_points, left_arrays, left_fields, left_pose = left
    right_config, right_points, right_arrays, right_fields, right_pose = right
    left_values = left_fields[field]
    right_values = right_fields[field]
    limits = (
        (0.0, 1.0)
        if field == "joint_liquefaction_indicator"
        else _robust_limits([left_values, right_values], include=threshold)
    )
    difference = difference_sign * _aligned_difference(
        left_arrays["ids"], left_values, right_arrays["ids"], right_values
    )
    finite_difference = difference[np.isfinite(difference)]
    difference_limit = float(np.percentile(np.abs(finite_difference), 99.0))
    if field == "joint_liquefaction_indicator":
        difference_limit = max(difference_limit, 1.0)
    else:
        difference_limit = max(difference_limit, 1.0e-12)
    difference_title = difference_label
    if field == "vertical_stress_remaining_ratio":
        difference_title += f" (positive = more loss in {left_code})"
    spatial_limits = (
        float(min(np.min(left_points[:, 0]), np.min(right_points[:, 0]))),
        float(max(np.max(left_points[:, 0]), np.max(right_points[:, 0]))),
        float(min(np.min(left_points[:, 1]), np.min(right_points[:, 1]))),
        float(max(np.max(left_points[:, 1]), np.max(right_points[:, 1])) + 0.005),
    )
    for axis, code, payload, values in (
        (axes[0], left_code, left, left_values),
        (axes[1], right_code, right, right_values),
    ):
        config, points, _, _, pose = payload
        artist = _draw_field(
            axis,
            points,
            values,
            config,
            title=code,
            cmap=cmap,
            limits=limits,
            threshold=threshold,
            pipeline_pose=pose,
            spatial_limits=spatial_limits,
        )
    plt.colorbar(artist, ax=list(axes[:2]), fraction=0.025, pad=0.015, label=label)
    difference_artist = _draw_field(
        axes[2],
        left_points,
        difference,
        left_config,
        title=f"Lagrangian ID-matched: {difference_title}\n(left current coordinates)",
        cmap="RdBu_r",
        limits=(-difference_limit, difference_limit),
        pipeline_pose=left_pose,
        spatial_limits=spatial_limits,
        draw_pipeline=False,
    )
    if finite_difference.size and np.allclose(finite_difference, 0.0):
        axes[2].text(
            0.5,
            0.08,
            "No ID-matched difference at this registered frame",
            transform=axes[2].transAxes,
            ha="center",
            va="bottom",
            fontsize=7.5,
            color="#333333",
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.85},
        )
    plt.colorbar(
        difference_artist,
        ax=axes[2],
        fraction=0.05,
        pad=0.02,
        label=difference_label,
    )
    return {
        "field": field,
        "shared_limits": list(limits),
        "threshold": threshold,
        "difference_convention": difference_label,
        "difference_basis": (
            "Lagrangian particle-ID-matched scalar difference, rendered at "
            f"the {left_code} current particle coordinates; no pipeline is drawn "
            "in the difference panel"
        ),
        "difference_limits": [-difference_limit, difference_limit],
        "difference_is_identically_zero": bool(
            finite_difference.size and np.allclose(finite_difference, 0.0)
        ),
        "shared_spatial_limits_m": list(spatial_limits),
        "pipeline_pose_sources": {
            left_code: {
                "center_m": np.asarray(left_pose["center"]).tolist(),
                "history_time_s": left_pose["history_time_s"],
                "history_time_offset_s": left_pose["history_time_offset_s"],
                "pressure_reference_audit": left_pose[
                    "pressure_reference_audit"
                ],
            },
            right_code: {
                "center_m": np.asarray(right_pose["center"]).tolist(),
                "history_time_s": right_pose["history_time_s"],
                "history_time_offset_s": right_pose["history_time_offset_s"],
                "pressure_reference_audit": right_pose[
                    "pressure_reference_audit"
                ],
            },
        },
        "finite_particles": {
            left_code: int(np.count_nonzero(np.isfinite(left_values))),
            right_code: int(np.count_nonzero(np.isfinite(right_values))),
            "ID_matched_difference": int(
                np.count_nonzero(np.isfinite(difference))
            ),
        },
    }


def plot_2d_field_comparisons(directory: Path, output: Path) -> dict[str, Any]:
    """Render physical and phase-only two-dimensional field comparisons."""

    physical_selection = select_common_field_frames(
        directory, codes=("HS", "LS"), driver="HS", metric=IF_AREA
    )
    phase_selections = registered_phase_field_frames(directory)
    manifest: dict[str, Any] = {
        "schema": "pipeline-2d-field-comparisons-v1",
        "pressure_reference": (
            "PIC_pore_pressure_excess = current mixture pressure minus the "
            "checkpoint-resumed initial pore pressure, evaluated per particle"
        ),
        "R_sigma_eligibility": (
            f"plot only particles with |initial sigma'_yy| >= "
            f"{MINIMUM_REFERENCE_VERTICAL_STRESS_PA:g} Pa; other values are masked"
        ),
        "support_zone_outline": (
            "orange dashed rectangle is the pre-registered initial support cohort: "
            "|x-x_pipe_initial|<=1.5D and bed_y-2D<=y<=bed_y"
        ),
        "pipeline_pose": (
            "case panels use the nearest saved rigid-pipeline history pose; "
            "difference panels omit the pipeline"
        ),
        "comparisons": {},
    }
    rows = (
        (
            "wave_excess_mixture_pressure_kpa",
            "Excess mixture pressure (kPa)",
            "viridis",
            "$\\Delta p$ = {left} - {right} (kPa)",
            1.0,
            None,
        ),
        (
            "upward_seepage_IF",
            "Upward seepage index, IF",
            "magma",
            "$\\Delta IF$ = {left} - {right}",
            1.0,
            1.0,
        ),
        (
            "vertical_stress_remaining_ratio",
            "Stress remaining ratio, $R_\\sigma$",
            "viridis",
            "$\\Delta R_\\sigma$ = {right} - {left}",
            -1.0,
            0.05,
        ),
        (
            "joint_liquefaction_indicator",
            "Joint indicator, IF$\\geq$1 and $R_\\sigma\\leq$0.05",
            "cividis",
            "$\\Delta$joint = {left} - {right}",
            1.0,
            0.5,
        ),
    )
    definitions: list[
        tuple[str, str, str, str, dict[str, Any], str]
    ] = [
        (
            "physical_LS_HS_t_IF",
            "HS",
            "LS",
            "physical t_IF",
            physical_selection,
            "figure7_2d_physical_fields_t_IF",
        )
    ]
    for registered_name, selection in phase_selections.items():
        definitions.append(
            (
                f"phase_only_RL_RE_{registered_name}",
                "RL",
                "RE",
                registered_name,
                selection,
                f"figure8_2d_phase_fields_{registered_name}",
            )
        )

    for key, left_code, right_code, registered_name, selection, filename in definitions:
        left = _field_frame(directory, left_code, selection["cases"][left_code])
        right = _field_frame(directory, right_code, selection["cases"][right_code])
        figure, axes = plt.subplots(4, 3, figsize=(15.2, 10.2), constrained_layout=True)
        figure.suptitle(
            f"{left_code} vs {right_code}: {registered_name}, "
            f"target t={selection['target_time_s']:.4g} s",
            fontsize=11,
        )
        field_rows = []
        for row_index, (field, label, cmap, difference_label, sign, threshold) in enumerate(rows):
            field_rows.append(
                _plot_field_triptych(
                    axes[row_index],
                    left_code=left_code,
                    right_code=right_code,
                    left=left,
                    right=right,
                    field=field,
                    label=label,
                    cmap=cmap,
                    difference_label=difference_label.format(
                        left=left_code, right=right_code
                    ),
                    difference_sign=sign,
                    threshold=threshold,
                )
            )
        output.mkdir(parents=True, exist_ok=True)
        figure.savefig(output / f"{filename}.png", dpi=300, bbox_inches="tight")
        figure.savefig(output / f"{filename}.pdf", bbox_inches="tight")
        plt.close(figure)
        manifest["comparisons"][key] = {
            "registered_time": registered_name,
            "selection": selection,
            "fields": field_rows,
            "figure_png": str((output / f"{filename}.png").resolve()),
            "figure_pdf": str((output / f"{filename}.pdf").resolve()),
        }
    (output / "two_dimensional_field_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return manifest


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--analysis-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--phase-only",
        action="store_true",
        help="skip the unavailable RL/RM constitutive curve figure",
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    plot_hydraulic_bridge(args.analysis_dir, args.output_dir / "figure7_hydraulic_bridge")
    plot_phase_control(args.analysis_dir, args.output_dir / "figure8_phase_control")
    if not args.phase_only:
        plot_constitutive(args.analysis_dir, args.output_dir / "figure9_constitutive")
    frames = common_field_frames(args.analysis_dir)
    (args.output_dir / "figure8_common_field_frames.json").write_text(
        json.dumps(frames, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    plot_2d_field_comparisons(args.analysis_dir, args.output_dir)
    print(f"Wrote manuscript figure drafts to {args.output_dir.resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
