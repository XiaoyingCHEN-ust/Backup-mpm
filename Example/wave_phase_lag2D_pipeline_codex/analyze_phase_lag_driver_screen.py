#!/usr/bin/env python3
"""Audit a fixed-state, pressure-only phase-lag gradient-force screen.

The lagged and locally phase-erased pressure databases are differentiated at
the same saved HD material-point coordinates. Particle state, density, volume,
stress and support membership are therefore identical in each pair; only the
prescribed pressure field changes. This is an inexpensive counterfactual
diagnostic, not a substitute for independently integrated RL/RE replay cases
or a released-pipeline engineering result.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.patches import Circle  # noqa: E402

import analyze_phase_lag_replay_exploration as replay_analysis  # noqa: E402
import plot_liquid_pressure_phase_lag as phase_plot  # noqa: E402
import plot_pipeline_gradient_force as gradient_force_module  # noqa: E402
import run_phase_lag_exploration as runner  # noqa: E402
from analyze_study import (  # noqa: E402
    expected_particle_steps,
    nominal_reference_particle_area,
    particle_step,
    pipeline_support_roi,
    read_initial_coordinates,
    vertical_stress_loss_state,
)
from plot_pipeline_gradient_force import (  # noqa: E402
    display_smoothed_grid,
    ordered_frame,
    raw_finite_difference_gradient,
    raw_force_metrics,
    read_raw_pressure_database,
    upward_excess_force_density,
)


SCHEMA = "pipeline-phase-lag-driver-screen-v2"
FIT_START_S = 1.3
FIT_END_S = 3.9
PERIOD_S = 1.3
EXPECTED_FRAMES = 60
MINIMUM_FORCE_IMPULSE_FRACTION = 0.05
MINIMUM_JOINT_AREA_TIME_FRACTION = 0.05
AREA_COMPARISON_ABS_TOL_M2 = 1.0e-12
MINIMUM_CROWN_PHASE_LAG_DEG = 18.0
MINIMUM_CROWN_AMPLITUDE_RATIO = 0.05
MINIMUM_CROWN_HARMONIC_R_SQUARED = 0.8
PHASE_RUNNER_BINDING_SCHEMA = "pipeline-phase-field-runner-input-binding-v1"


class ScreenError(RuntimeError):
    """Raised when a driver screen is incomplete or cannot be trusted."""


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _relative_fraction(lagged: float, erased: float) -> float | None:
    if not math.isfinite(lagged) or not math.isfinite(erased):
        raise ScreenError("Integrated screen metric is non-finite")
    if erased <= 0.0:
        return None
    return lagged / erased - 1.0


def integrate_window(
    rows: list[dict[str, float]], start_s: float, end_s: float
) -> dict[str, Any]:
    """Integrate paired metrics on one exact, inclusive saved-time window."""

    if not rows or not (math.isfinite(start_s) and math.isfinite(end_s)):
        raise ScreenError("Integration window is empty or non-finite")
    if end_s <= start_s:
        raise ScreenError("Integration window must have positive duration")
    selected = [
        row
        for row in rows
        if start_s - 1.0e-12 <= row["time_s"] <= end_s + 1.0e-12
    ]
    if not selected:
        raise ScreenError("No saved rows lie in the requested window")
    times = np.asarray([row["time_s"] for row in selected], dtype=np.float64)
    if (
        not np.all(np.isfinite(times))
        or np.any(np.diff(times) <= 0.0)
        or not math.isclose(float(times[0]), start_s, abs_tol=1.0e-12)
        or not math.isclose(float(times[-1]), end_s, abs_tol=1.0e-12)
    ):
        raise ScreenError("Saved-time window is incomplete or unordered")

    output: dict[str, Any] = {
        "start_s": start_s,
        "end_s": end_s,
        "saved_frames": len(selected),
    }
    for metric, output_name, units in (
        ("signed_force_n_per_m", "signed_force_impulse", "N s/m"),
        ("net_uplift_force_n_per_m", "net_uplift_force_impulse", "N s/m"),
        (
            "positive_force_n_per_m",
            "local_positive_force_activity_impulse",
            "N s/m",
        ),
        ("critical_area_m2", "IF_ge_1_area_time", "m2 s"),
        ("joint_area_m2", "shared_HD_joint_area_time", "m2 s"),
    ):
        lagged = np.asarray(
            [row[f"lagged_{metric}"] for row in selected], dtype=np.float64
        )
        erased = np.asarray(
            [row[f"phase_erased_{metric}"] for row in selected], dtype=np.float64
        )
        if not np.all(np.isfinite(lagged)) or not np.all(np.isfinite(erased)):
            raise ScreenError(f"{metric} contains NaN/Inf")
        lagged_integral = float(np.trapz(lagged, times))
        erased_integral = float(np.trapz(erased, times))
        output[output_name] = {
            "units": units,
            "lagged": lagged_integral,
            "phase_erased": erased_integral,
            "lagged_minus_phase_erased": lagged_integral - erased_integral,
            "lagged_minus_phase_erased_fraction": _relative_fraction(
                lagged_integral, erased_integral
            ),
        }
    return output


def select_post_ramp_maximum(
    rows: list[dict[str, float]], metric: str
) -> dict[str, float]:
    """Select a registered post-ramp frame, never an early ramp transient."""

    eligible = [
        row
        for row in rows
        if FIT_START_S - 1.0e-12 <= row["time_s"] <= FIT_END_S + 1.0e-12
    ]
    if not eligible:
        raise ScreenError("No post-ramp row is available for selection")
    if any(metric not in row or not math.isfinite(row[metric]) for row in eligible):
        raise ScreenError(f"Post-ramp selection metric {metric} is missing/non-finite")
    return max(eligible, key=lambda row: row[metric])


def phase_field_qualification(audit: dict[str, Any]) -> dict[str, Any]:
    """Apply the locked crown phase-resolution qualification."""

    try:
        crown = audit["probe_results"]["crown"]
        surface = audit["probe_results"]["surface"]
        crown_amplitude = float(crown["fundamental_amplitude_pa"])
        surface_amplitude = float(surface["fundamental_amplitude_pa"])
        harmonic_r_squared = float(crown["harmonic_r_squared"])
    except (KeyError, TypeError, ValueError) as error:
        raise ScreenError("Phase-field audit lacks finite crown/surface probes") from error
    phase_value = crown.get("phase_difference_from_same_column_surface_deg")
    if phase_value is None:
        phase: float | None = None
    else:
        try:
            phase = float(phase_value)
        except (TypeError, ValueError) as error:
            raise ScreenError("Phase-field crown phase is invalid") from error
        if not math.isfinite(phase):
            raise ScreenError("Phase-field crown phase is non-finite")
    values = (crown_amplitude, surface_amplitude, harmonic_r_squared)
    if (
        not all(math.isfinite(value) for value in values)
        or crown_amplitude < 0.0
        or surface_amplitude <= 0.0
    ):
        raise ScreenError("Phase-field crown/surface probes are non-finite")
    amplitude_ratio = crown_amplitude / surface_amplitude
    checks = {
        "crown_abs_same_column_phase_lag_ge_18deg": (
            phase is not None and abs(phase) >= MINIMUM_CROWN_PHASE_LAG_DEG
        ),
        "crown_amplitude_ratio_ge_0p05": (
            amplitude_ratio >= MINIMUM_CROWN_AMPLITUDE_RATIO
        ),
        "crown_harmonic_r_squared_ge_0p8": (
            harmonic_r_squared >= MINIMUM_CROWN_HARMONIC_R_SQUARED
        ),
    }
    return {
        "passed": all(checks.values()),
        "checks": checks,
        "observed": {
            "crown_same_column_phase_lag_deg": phase,
            "crown_amplitude_ratio": amplitude_ratio,
            "crown_harmonic_r_squared": harmonic_r_squared,
        },
    }


def validate_phase_field_audit(
    root: Path,
    label: str,
    hd_path: Path,
    parameters: replay_analysis.ReplayParameters,
    runner_state: dict[str, Any],
) -> dict[str, Any]:
    """Bind the driver screen to one exact unsmoothed phase-field audit."""

    directory = root / "analysis/phase_lag_exploratory" / label / "phase_v3"
    candidates = sorted(directory.glob("*liquid_pressure_phase_lag.audit.json"))
    if len(candidates) != 1:
        raise ScreenError("Expected exactly one phase-v3 field audit")
    path = candidates[0]
    audit = read_json(path)
    case = audit.get("case")
    if not isinstance(case, str) or not case or Path(case).name != case:
        raise ScreenError("Phase-field audit case identity is invalid")
    expected_figure_paths = tuple(
        path.parent / f"{case}_{stem}{suffix}"
        for stem in (
            "liquid_pressure_2d_phase_lag",
            "liquid_pressure_time_depth",
        )
        for suffix in (".png", ".pdf")
    )
    try:
        stages = runner_state["stages"]
        eq_stage = stages["eq"]
        hd_stage = stages["hd"]
        eq_audit = eq_stage["audit"]
        hd_audit = hd_stage["audit"]
        particle_input = eq_audit["stability_qa"]["normal_stress"]["inputs"][
            "particles"
        ]
        raw_particle_path = Path(str(particle_input["path"]))
        particle_path = (
            raw_particle_path
            if raw_particle_path.is_absolute()
            else root / raw_particle_path
        ).resolve()
        runner_particles = {
            "path": str(particle_path),
            "size_bytes": int(particle_input["size_bytes"]),
            "sha256": str(particle_input["sha256"]),
        }
        runner_binding = {
            "schema": PHASE_RUNNER_BINDING_SCHEMA,
            "eq": {
                "checkpoint_config": eq_stage["config"],
                "checkpoint": eq_audit["final_vtp"],
                "particles": runner_particles,
            },
            "hd": {
                "analysis_config": hd_stage["config"],
                "result_directory": str(Path(hd_audit["result_directory"]).resolve()),
                "result_particle_vtp": hd_audit["particle_vtp"],
                "pipeline_history": hd_audit["pipeline_history"],
            },
        }
    except (KeyError, TypeError, ValueError) as error:
        raise ScreenError("Runner EQ/HD input provenance is incomplete") from error
    expected_input_paths = {
        "analysis_config": Path(runner_binding["hd"]["analysis_config"]["path"]),
        "checkpoint_config": Path(
            runner_binding["eq"]["checkpoint_config"]["path"]
        ),
        "checkpoint": Path(runner_binding["eq"]["checkpoint"]["path"]),
        "particles": Path(runner_binding["eq"]["particles"]["path"]),
        "pipeline_history": Path(
            runner_binding["hd"]["pipeline_history"]["path"]
        ),
        "result_directory": Path(runner_binding["hd"]["result_directory"]),
        "result_particle_vtp": [
            Path(record["path"])
            for record in runner_binding["hd"]["result_particle_vtp"]
        ],
    }
    try:
        phase_plot.validate_phase_audit_artifacts(
            audit,
            expected_figure_paths=expected_figure_paths,
            expected_input_paths=expected_input_paths,
        )
    except (phase_plot.PhaseAuditError, OSError) as error:
        raise ScreenError(
            "Phase-field audit artifacts are incomplete or stale"
        ) from error
    if (
        audit.get("schema") != "pipeline-liquid-pressure-phase-lag-audit-v3"
        or audit.get("require_unsmoothed") is not True
        or int(audit.get("particle_count", 0)) != 7376
        or int(audit.get("frame_count", 0)) != EXPECTED_FRAMES
        or audit.get("parameters") != parameters.as_dict()
        or audit.get("fit_window_s") != [FIT_START_S, FIT_END_S]
        or not math.isclose(
            float(audit.get("wave_period_s", math.nan)), PERIOD_S, abs_tol=1.0e-12
        )
    ):
        raise ScreenError("Phase-field audit parameters/schema are incompatible")
    phase_inputs = audit["artifacts"]["inputs"]
    exact_bindings = {
        "analysis_config": runner_binding["hd"]["analysis_config"],
        "checkpoint_config": runner_binding["eq"]["checkpoint_config"],
        "checkpoint": runner_binding["eq"]["checkpoint"],
        "particles": runner_binding["eq"]["particles"],
        "pipeline_history": runner_binding["hd"]["pipeline_history"],
    }
    for name, runner_record in exact_bindings.items():
        if phase_inputs.get(name) != runner_record:
            raise ScreenError(
                f"Phase-field {name} does not exactly bind the runner artifact"
            )
    if (
        phase_inputs.get("result_directory")
        != runner_binding["hd"]["result_directory"]
        or phase_inputs.get("result_particle_vtp")
        != runner_binding["hd"]["result_particle_vtp"]
    ):
        raise ScreenError(
            "Phase-field result VTP inputs do not exactly bind all runner HD artifacts"
        )
    if (
        Path(str(exact_bindings["analysis_config"]["path"])).resolve()
        != hd_path.resolve()
    ):
        raise ScreenError("Phase-field audit is not bound to the HD config")
    return {
        "path": str(path.resolve()),
        "sha256": sha256(path),
        "qualification": phase_field_qualification(audit),
        "runner_binding": runner_binding,
    }


def diagnostic_gate(
    rows: list[dict[str, float]],
    windows: list[dict[str, Any]],
    *,
    minimum_resolved_area_m2: float,
    phase_qualification: dict[str, Any],
) -> dict[str, Any]:
    """Apply the predeclared pressure-only screening guardrails."""

    if len(windows) != 3:
        raise ScreenError("Diagnostic gate requires two cycles and their union")
    if not math.isfinite(minimum_resolved_area_m2) or minimum_resolved_area_m2 <= 0.0:
        raise ScreenError("Minimum resolved area must be finite and positive")

    combined = windows[-1]
    net_uplift_fraction = combined["net_uplift_force_impulse"][
        "lagged_minus_phase_erased_fraction"
    ]
    joint_fraction = combined["shared_HD_joint_area_time"][
        "lagged_minus_phase_erased_fraction"
    ]
    net_uplift_each_cycle_positive = all(
        window["net_uplift_force_impulse"]["lagged_minus_phase_erased"] > 0.0
        for window in windows[:-1]
    )
    signed_force_each_cycle_positive = all(
        window["signed_force_impulse"]["lagged_minus_phase_erased"] > 0.0
        for window in windows[:-1]
    )
    if_area_direction_positive = (
        combined["IF_ge_1_area_time"]["lagged_minus_phase_erased"] > 0.0
    )
    selected = [
        row
        for row in rows
        if FIT_START_S - 1.0e-12 <= row["time_s"] <= FIT_END_S + 1.0e-12
    ]
    area_tolerance = max(
        AREA_COMPARISON_ABS_TOL_M2,
        64.0 * np.finfo(np.float64).eps * max(1.0, minimum_resolved_area_m2),
    )
    resolved_threshold = minimum_resolved_area_m2 - area_tolerance
    joint_advantage_resolved = np.asarray(
        [
            row["lagged_joint_area_m2"] >= resolved_threshold
            and row["lagged_joint_area_m2"]
            - row["phase_erased_joint_area_m2"]
            >= resolved_threshold
            for row in selected
        ],
        dtype=bool,
    )
    consecutive_joint_advantage = bool(
        joint_advantage_resolved.size >= 2
        and np.any(joint_advantage_resolved[:-1] & joint_advantage_resolved[1:])
    )
    checks = {
        "crown_phase_lag_is_resolved": bool(phase_qualification.get("passed")),
        "both_post_ramp_cycles_have_positive_signed_force_delta": (
            signed_force_each_cycle_positive
        ),
        "combined_signed_force_impulse_delta_positive": (
            combined["signed_force_impulse"]["lagged_minus_phase_erased"] > 0.0
        ),
        "both_post_ramp_cycles_have_positive_net_uplift_delta": (
            net_uplift_each_cycle_positive
        ),
        "combined_net_uplift_impulse_increase_ge_5pct": (
            net_uplift_fraction is not None
            and net_uplift_fraction >= MINIMUM_FORCE_IMPULSE_FRACTION
        ),
        "combined_IF_ge_1_area_time_direction_positive": if_area_direction_positive,
        "combined_shared_HD_joint_area_time_nonzero": (
            combined["shared_HD_joint_area_time"]["lagged"] > 0.0
        ),
        "combined_shared_HD_joint_area_time_increase_ge_5pct": (
            joint_fraction is not None
            and joint_fraction >= MINIMUM_JOINT_AREA_TIME_FRACTION
        ),
        "resolved_joint_advantage_for_two_consecutive_saved_frames": (
            consecutive_joint_advantage
        ),
    }
    return {
        "passed": all(checks.values()),
        "checks": checks,
        "classification": (
            "pressure-only-screen-passed"
            if all(checks.values())
            else "pressure-only-screen-not-passed"
        ),
        "minimum_resolved_joint_advantage_area_m2": minimum_resolved_area_m2,
        "area_comparison_absolute_tolerance_m2": area_tolerance,
        "phase_qualification": phase_qualification,
        "guardrail": (
            "Passing permits an independently integrated RL/RE diagnostic; it "
            "does not establish released-pipeline response, a causal liquefaction "
            "difference, or manuscript evidence. The joint metric is only spatial "
            "co-location with the one shared HD stress-loss field."
        ),
    }


def validate_phase_database_bindings(
    root: Path,
    phase_audit: dict[str, Any],
    lagged_database: dict[str, Any],
    erased_database: dict[str, Any],
) -> None:
    """Bind transform metadata to the exact databases used in this analysis."""

    metadata = phase_audit.get("metadata", {})
    for group, observed in (
        ("source", lagged_database),
        ("output", erased_database),
    ):
        recorded = metadata.get(group, {})
        for kind in ("points", "values"):
            path = Path(str(recorded.get(kind, "")))
            if not path.is_absolute():
                path = root / path
            if (
                path.resolve() != Path(str(observed[kind])).resolve()
                or recorded.get(f"{kind}_sha256")
                != observed[f"{kind}_sha256"]
            ):
                raise ScreenError(
                    f"Phase metadata {group} {kind} is not the compared database"
                )

    for field in ("step_interval", "max_step"):
        if lagged_database[field] != erased_database[field]:
            raise ScreenError(f"Pressure databases disagree on {field}")
    if not math.isclose(
        float(lagged_database["source_dt_s"]),
        float(erased_database["source_dt_s"]),
        rel_tol=0.0,
        abs_tol=1.0e-15,
    ):
        raise ScreenError("Pressure databases disagree on source_dt_s")

    alignment = metadata.get("qa", {}).get("phase_alignment", {}).get("liquid", {})
    phase_defined = int(alignment.get("phase_defined_particle_count", 0))
    changed = int(alignment.get("changed_particle_count", 0))
    source_difference = float(
        alignment.get("max_abs_source_to_surface_phase_difference_deg", math.nan)
    )
    residual = float(
        alignment.get("max_abs_residual_to_surface_phase_difference_deg", math.nan)
    )
    if (
        phase_defined <= 0
        or changed <= 0
        or not math.isfinite(source_difference)
        or source_difference <= 0.0
        or not math.isfinite(residual)
        or residual > 1.0e-8
    ):
        raise ScreenError(
            "Phase metadata does not prove a non-trivial liquid-phase alignment"
        )


def plot_diagnostic(
    output_dir: Path,
    rows: list[dict[str, float]],
    selected_field: dict[str, Any],
    config: dict[str, Any],
    support: np.ndarray,
    display_smoothing_sigma: float,
    diagnostic_gate_result: dict[str, Any],
    combined_net_uplift_fraction: float | None,
) -> dict[str, Any]:
    """Plot raw metrics and display-only rasters for the selected frame."""

    if not math.isfinite(display_smoothing_sigma) or display_smoothing_sigma < 0.0:
        raise ScreenError("Display smoothing sigma must be finite and non-negative")
    center = np.asarray(
        config["analysis"]["rigid_pipeline"]["center"], dtype=np.float64
    )
    radius = float(config["analysis"]["rigid_pipeline"]["outer_radius"])
    diameter = 2.0 * radius
    fluid = config["materials"][1]
    bed_y = float(fluid["sea_level"]) - float(fluid["depth_left"])
    bounds = (
        center[0] - 1.5 * diameter,
        center[0] + 1.5 * diameter,
        bed_y - 2.0 * diameter,
        bed_y,
    )

    lagged_if = np.asarray(selected_field["lagged_upward_if"], dtype=np.float64)
    erased_if = np.asarray(
        selected_field["phase_erased_upward_if"], dtype=np.float64
    )
    difference = lagged_if - erased_if
    common_max = max(
        float(np.max(lagged_if[support])),
        float(np.max(erased_if[support])),
        1.0,
    )
    difference_max = max(float(np.max(np.abs(difference[support]))), 1.0e-12)

    figure = plt.figure(figsize=(13.2, 7.8), constrained_layout=True)
    grid = figure.add_gridspec(2, 3, height_ratios=(1.0, 0.9))
    field_axes = [figure.add_subplot(grid[0, index]) for index in range(3)]
    history_axes = [figure.add_subplot(grid[1, index]) for index in range(3)]
    field_specs = (
        (lagged_if, "Lagged raw upward-force index", "magma", 0.0, common_max),
        (
            erased_if,
            "Phase-erased raw upward-force index",
            "magma",
            0.0,
            common_max,
        ),
        (
            difference,
            "Instantaneous raw IF difference (lagged - erased)",
            "RdBu_r",
            -difference_max,
            difference_max,
        ),
    )
    artists = []
    for axis, (values, title, cmap, lower, upper) in zip(field_axes, field_specs):
        _, _, raster = display_smoothed_grid(
            selected_field["points"],
            values,
            config,
            center,
            bounds,
            display_smoothing_sigma,
        )
        artist = axis.imshow(
            raster,
            origin="lower",
            extent=bounds,
            cmap=cmap,
            vmin=lower,
            vmax=upper,
            interpolation="nearest",
            aspect="equal",
        )
        artists.append(artist)
        axis.add_patch(Circle(center, radius, facecolor="white", edgecolor="0.1"))
        axis.set_title(title, fontsize=9.5)
        axis.set_xlabel("x (m)")
        axis.set_ylabel("y (m)")
    figure.colorbar(
        artists[0], ax=field_axes[:2], label="Upward seepage-force index, IF"
    )
    figure.colorbar(artists[2], ax=field_axes[2], label="Difference in IF")

    times = np.asarray([row["time_s"] for row in rows], dtype=np.float64)
    colors = {"lagged": "#B13A3A", "phase_erased": "#2468A2"}
    labels = {"lagged": "Lagged", "phase_erased": "Phase-erased"}
    joint_axis = history_axes[2].twinx()
    for role in ("lagged", "phase_erased"):
        history_axes[0].plot(
            times,
            [row[f"{role}_signed_force_n_per_m"] for row in rows],
            color=colors[role],
            label=labels[role],
        )
        history_axes[1].plot(
            times,
            [row[f"{role}_positive_force_n_per_m"] for row in rows],
            color=colors[role],
            label=labels[role],
        )
        history_axes[2].plot(
            times,
            [row[f"{role}_critical_area_m2"] for row in rows],
            color=colors[role],
            label=f"{labels[role]} IF>=1",
        )
        joint_axis.plot(
            times,
            [row[f"{role}_joint_area_m2"] for row in rows],
            color=colors[role],
            linestyle="--",
            label=f"{labels[role]} co-location",
        )
    history_axes[0].set_title("Signed support-zone force")
    history_axes[0].set_ylabel("Signed force (N/m)")
    history_axes[1].set_title("Local upward-positive force activity")
    history_axes[1].set_ylabel("Local positive-part activity (N/m)")
    history_axes[2].set_title("IF threshold and shared-HD co-location")
    history_axes[2].set_ylabel("IF>=1 area (m$^2$)")
    joint_axis.set_ylabel("Shared-HD co-location area (m$^2$)")
    for axis in history_axes:
        for time_s in (FIT_START_S, FIT_START_S + PERIOD_S, FIT_END_S):
            axis.axvline(time_s, color="0.55", linewidth=0.7, linestyle=":")
        axis.axvline(
            selected_field["time_s"], color="0.2", linewidth=0.8, linestyle="--"
        )
        axis.set_xlabel("Time (s)")
        axis.grid(color="0.88", linewidth=0.6)
        if axis is not history_axes[2]:
            axis.legend(fontsize=7.2)
    joint_axis.axvline(
        selected_field["time_s"], color="0.2", linewidth=0.8, linestyle="--"
    )
    handles, legend_labels = history_axes[2].get_legend_handles_labels()
    joint_handles, joint_labels = joint_axis.get_legend_handles_labels()
    history_axes[2].legend(
        handles + joint_handles, legend_labels + joint_labels, fontsize=6.8
    )

    parameters = replay_analysis.replay_parameters(config, "driver")
    passed = bool(diagnostic_gate_result["passed"])
    if combined_net_uplift_fraction is None:
        combined_delta = "undefined"
    else:
        combined_delta = f"{100.0 * combined_net_uplift_fraction:+.3f}%"
    figure.suptitle(
        "Fixed-state pressure-lag gradient-force diagnostic\n"
        f"k={parameters.intrinsic_permeability_m2:.3e} m$^2$, "
        f"$S_w$={parameters.liquid_saturation:.3f}; raw metrics; "
        f"raster display smoothing only, sigma={display_smoothing_sigma:g} pixels\n"
        f"SCREEN {'PASSED' if passed else 'FAILED'}: combined post-ramp "
        f"net-uplift impulse delta={combined_delta}; selected t="
        f"{selected_field['time_s']:.3f} s is the instantaneous maximum only",
        fontsize=12.5,
        color="0.12" if passed else "#8B1A1A",
    )
    base = output_dir / "pressure_only_gradient_force"
    png = base.with_suffix(".png")
    pdf = base.with_suffix(".pdf")
    figure.savefig(png, dpi=300, bbox_inches="tight")
    figure.savefig(pdf, bbox_inches="tight")
    plt.close(figure)
    return {
        "display_smoothing_sigma_pixels": display_smoothing_sigma,
        "display_smoothing_used_in_metrics": False,
        "png": {"path": str(png.resolve()), "sha256": sha256(png)},
        "pdf": {"path": str(pdf.resolve()), "sha256": sha256(pdf)},
    }


def verify_runner_state(root: Path, label: str) -> tuple[dict[str, Any], Path]:
    audit_path = runner.state_path(root, label)
    state = runner.load_state(root, label)
    if state.get("stage") != "hd_complete":
        raise ScreenError(f"HD runner stage is not complete: {state.get('stage')}")
    runner.verify_provenance(root, label, state)
    hd_audit = state.get("stages", {}).get("hd", {}).get("audit", {})
    if (
        hd_audit.get("health", {}).get("saved_frames") != EXPECTED_FRAMES
        or hd_audit.get("health", {}).get("outside_mesh_maximum_count") != 0
    ):
        raise ScreenError("Recorded HD health audit is incomplete")
    particle_records = hd_audit.get("particle_vtp", [])
    if not isinstance(particle_records, list) or len(particle_records) != EXPECTED_FRAMES:
        raise ScreenError("Runner audit does not bind exactly 60 HD particle VTPs")
    for record in particle_records:
        runner.verify_artifact(record, "HD particle VTP")
    runner.verify_artifact(hd_audit["pipeline_history"], "HD pipeline history")
    runner.verify_artifact(
        hd_audit["pressure_database"]["points"], "lagged pressure points"
    )
    runner.verify_artifact(
        hd_audit["pressure_database"]["values"], "lagged pressure values"
    )
    return state, audit_path


def analyze(
    root: Path, label: str, *, display_smoothing_sigma: float = 0.75
) -> dict[str, Any]:
    root = root.resolve()
    state, runner_audit_path = verify_runner_state(root, label)
    config_dir = root / "configs/phase_lag_exploratory" / label
    hd_path = config_dir / "03_HD.json"
    rl_path = config_dir / "04_RL.json"
    re_path = config_dir / "05_RE.json"
    hd = read_json(hd_path)
    rl = read_json(rl_path)
    re_case = read_json(re_path)
    parameters = replay_analysis.validate_replay_pair_configs(rl, re_case)
    if hd["analysis"]["pressure_smoothing"] is not False or hd["analysis"][
        "pressure_smoothing_in_loop"
    ] is not False:
        raise ScreenError("HD pressure smoothing must be explicitly disabled")
    if hd["analysis"]["rigid_pipeline"]["fixed"] is not True:
        raise ScreenError("HD pipeline must be explicitly fixed")
    phase_audit = replay_analysis.validate_phase_metadata(root, label)
    phase_field_audit = validate_phase_field_audit(
        root, label, hd_path, parameters, state
    )

    result = root / hd["post_processing"]["path"] / hd["analysis"]["uuid"]
    files = sorted(result.glob("particle*.vtp"), key=particle_step)
    expected = expected_particle_steps(hd)
    if (
        len(files) != EXPECTED_FRAMES
        or [particle_step(path) for path in files] != expected
    ):
        raise ScreenError("HD particle frame grid is not exactly 60 saved frames")
    recorded_files = [
        Path(str(record["path"])).resolve()
        for record in state["stages"]["hd"]["audit"]["particle_vtp"]
    ]
    if recorded_files != [path.resolve() for path in files]:
        raise ScreenError(
            "Runner VTP artifact list differs from the analyzed frame grid"
        )

    reference = read_initial_coordinates(
        root / hd["particles"][0]["generator"]["location"]
    )
    support = pipeline_support_roi(hd, reference)
    if support is None or not np.any(support):
        raise ScreenError("Pipeline support cohort is empty")
    lagged_pressure, lagged_database = read_raw_pressure_database(
        root, rl, reference
    )
    erased_pressure, erased_database = read_raw_pressure_database(
        root, re_case, reference
    )
    validate_phase_database_bindings(
        root, phase_audit, lagged_database, erased_database
    )

    gravity = np.asarray(hd["external_loading_conditions"]["gravity"], dtype=float)
    minimum_resolved_joint_area = 2.0 * nominal_reference_particle_area(reference)
    rows: list[dict[str, float]] = []
    selected_field: dict[str, Any] | None = None
    best_net_uplift_delta = -math.inf
    for path in files:
        step = particle_step(path)
        if step not in lagged_pressure or step not in erased_pressure:
            raise ScreenError(f"Pressure databases omit saved HD step {step}")
        points, arrays = ordered_frame(path, len(reference))
        if not np.all(np.isfinite(points)):
            raise ScreenError(f"{path} has non-finite coordinates")
        volumes = np.asarray(arrays["volumes"], dtype=np.float64).reshape(-1)
        initial_vertical = np.asarray(
            arrays["initial_vertical_effective_stresses"], dtype=np.float64
        ).reshape(-1)
        solver_stress_ratio = np.asarray(
            arrays["vertical_effective_stress_remaining_ratios"], dtype=np.float64
        ).reshape(-1)
        if not (
            np.all(np.isfinite(volumes))
            and np.all(np.isfinite(initial_vertical))
            and np.all(np.isfinite(solver_stress_ratio))
            and np.all(volumes > 0.0)
        ):
            raise ScreenError(f"{path} has invalid mechanical support fields")
        if "stresses" not in arrays:
            raise ScreenError(f"{path} lacks stresses for R_sigma cross-check")
        stress_ratio, eligible, stress_loss = vertical_stress_loss_state(
            np.asarray(arrays["stresses"], dtype=np.float64)[:, 1],
            initial_vertical,
        )
        if not np.allclose(
            stress_ratio[eligible],
            solver_stress_ratio[eligible],
            rtol=1.0e-10,
            atol=1.0e-12,
        ):
            raise ScreenError(f"{path}: solver/Python R_sigma disagreement")

        row: dict[str, float] = {
            "time_s": step * float(hd["analysis"]["dt"])
        }
        upward_fields: dict[str, np.ndarray] = {}
        for role, pressure in (
            ("lagged", lagged_pressure[step]),
            ("phase_erased", erased_pressure[step]),
        ):
            gradient = raw_finite_difference_gradient(reference, points, pressure)
            excess = upward_excess_force_density(
                gradient, arrays["liquid_densities"], gravity
            )
            metrics, upward_if = raw_force_metrics(
                excess, arrays["gamma_sub"], volumes, support
            )
            upward_fields[role] = upward_if
            joint = support & stress_loss & (upward_if >= 1.0)
            row.update(
                {
                    f"{role}_signed_force_n_per_m": metrics.signed_force_n_per_m,
                    f"{role}_net_uplift_force_n_per_m": max(
                        metrics.signed_force_n_per_m, 0.0
                    ),
                    f"{role}_positive_force_n_per_m": (
                        metrics.positive_force_n_per_m
                    ),
                    f"{role}_critical_area_m2": metrics.critical_area_m2,
                    f"{role}_maximum_upward_if": metrics.maximum_upward_if,
                    f"{role}_joint_area_m2": float(np.sum(volumes[joint])),
                }
            )
        row["delta_positive_force_n_per_m"] = (
            row["lagged_positive_force_n_per_m"]
            - row["phase_erased_positive_force_n_per_m"]
        )
        row["delta_signed_force_n_per_m"] = (
            row["lagged_signed_force_n_per_m"]
            - row["phase_erased_signed_force_n_per_m"]
        )
        row["delta_net_uplift_force_n_per_m"] = (
            row["lagged_net_uplift_force_n_per_m"]
            - row["phase_erased_net_uplift_force_n_per_m"]
        )
        if (
            FIT_START_S - 1.0e-12 <= row["time_s"] <= FIT_END_S + 1.0e-12
            and row["delta_net_uplift_force_n_per_m"] > best_net_uplift_delta
        ):
            best_net_uplift_delta = row["delta_net_uplift_force_n_per_m"]
            selected_field = {
                "time_s": row["time_s"],
                "step": step,
                "points": points.copy(),
                "lagged_upward_if": upward_fields["lagged"].copy(),
                "phase_erased_upward_if": upward_fields["phase_erased"].copy(),
            }
        rows.append(row)

    windows = [
        integrate_window(rows, FIT_START_S, FIT_START_S + PERIOD_S),
        integrate_window(rows, FIT_START_S + PERIOD_S, FIT_END_S),
        integrate_window(rows, FIT_START_S, FIT_END_S),
    ]
    if selected_field is None:
        raise ScreenError("No post-ramp field was selected for the diagnostic figure")
    selected = select_post_ramp_maximum(
        rows, "delta_net_uplift_force_n_per_m"
    )

    output_dir = root / "analysis/phase_lag_exploratory" / label / "driver_screen"
    output_dir.mkdir(parents=True, exist_ok=True)
    csv_path = output_dir / "pressure_only_gradient_force_history.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    gate = diagnostic_gate(
        rows,
        windows,
        minimum_resolved_area_m2=minimum_resolved_joint_area,
        phase_qualification=phase_field_audit["qualification"],
    )
    combined_net_fraction = windows[-1]["net_uplift_force_impulse"][
        "lagged_minus_phase_erased_fraction"
    ]
    figure = plot_diagnostic(
        output_dir,
        rows,
        selected_field,
        hd,
        support,
        display_smoothing_sigma,
        gate,
        combined_net_fraction,
    )

    audit = {
        "schema": SCHEMA,
        "classification": "exploratory-pressure-only-fixed-state-diagnostic",
        "label": label,
        "analyzer": {
            "path": str(Path(__file__).resolve()),
            "sha256": sha256(Path(__file__).resolve()),
        },
        "implementation_dependencies": {
            "gradient_force_module": {
                "path": str(Path(gradient_force_module.__file__).resolve()),
                "sha256": sha256(Path(gradient_force_module.__file__).resolve()),
            },
            "replay_analyzer": {
                "path": str(Path(replay_analysis.__file__).resolve()),
                "sha256": sha256(Path(replay_analysis.__file__).resolve()),
            },
            "phase_plot": {
                "path": str(Path(phase_plot.__file__).resolve()),
                "sha256": sha256(Path(phase_plot.__file__).resolve()),
            },
            "runner": {
                "path": str(Path(runner.__file__).resolve()),
                "sha256": sha256(Path(runner.__file__).resolve()),
            },
        },
        "parameters": parameters.as_dict(),
        "definition": (
            "Lagged and phase-erased raw pressure databases are differentiated "
            "on the deformed material-point lattice by solving the complete 2x2 "
            "directional-difference system at each particle, then integrated at "
            "identical HD coordinates with the same particle volumes, densities, "
            "stresses and support cohort."
        ),
        "metric_definitions": {
            "signed_force_n_per_m": (
                "spatial integral of signed upward excess pressure-gradient force "
                "density over the registered support cohort"
            ),
            "net_uplift_force_n_per_m": (
                "positive part of the spatially integrated signed support force"
            ),
            "positive_force_n_per_m": (
                "spatial integral of each particle's local positive force part; "
                "not the net support resultant"
            ),
            "shared_HD_joint_area_m2": (
                "IF>=1 co-located with R_sigma<=0.05 in the one shared HD state"
            ),
        },
        "limitations": (
            "This isolates instantaneous pressure-phase forcing but omits the "
            "separate mechanical histories of RL/RE and all released-pipe feedback. "
            "The joint metric is co-location against one shared HD stress-loss "
            "field and cannot show that retained lag caused additional liquefaction."
        ),
        "screen_window": {
            "start_s": FIT_START_S,
            "end_s": FIT_END_S,
            "period_s": PERIOD_S,
            "classification": "two post-ramp cycles",
        },
        "frame_count": len(rows),
        "particle_count": len(reference),
        "support_particle_count": int(np.count_nonzero(support)),
        "runner_audit": {
            "path": str(runner_audit_path.resolve()),
            "sha256": sha256(runner_audit_path),
        },
        "configs": {
            "HD": {"path": str(hd_path.resolve()), "sha256": sha256(hd_path)},
            "RL": {"path": str(rl_path.resolve()), "sha256": sha256(rl_path)},
            "RE": {"path": str(re_path.resolve()), "sha256": sha256(re_path)},
        },
        "phase_transform": phase_audit,
        "phase_field": phase_field_audit,
        "pressure_databases": {
            "lagged": lagged_database,
            "phase_erased": erased_database,
        },
        "integration_windows": windows,
        "maximum_post_ramp_net_uplift_force_delta_frame": selected,
        "diagnostic_gate": gate,
        "figure": figure,
        "history_csv": str(csv_path.resolve()),
        "history_csv_sha256": sha256(csv_path),
    }
    audit_path = output_dir / "pressure_only_gradient_force.audit.json"
    audit["audit_path"] = str(audit_path.resolve())
    audit_path.write_text(
        json.dumps(audit, indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )
    return audit


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--label", required=True)
    parser.add_argument("--display-smoothing-sigma", type=float, default=0.75)
    args = parser.parse_args()
    root = Path(__file__).resolve().parent
    audit = analyze(
        root, args.label, display_smoothing_sigma=args.display_smoothing_sigma
    )
    print(
        json.dumps(
            {
                "diagnostic_gate": audit["diagnostic_gate"],
                "integration_windows": audit["integration_windows"],
                "audit_path": audit["audit_path"],
            },
            indent=2,
            allow_nan=False,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
