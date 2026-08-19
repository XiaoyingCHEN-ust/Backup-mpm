#!/usr/bin/env python3
"""Compare strong- and weak-phase-lag engineering parameter combinations.

This is deliberately separate from the same-parameter phase-erased
counterfactual.  One case has lower permeability and saturation; the physical
reference has higher permeability and saturation.  Their raw physical pressure
fields are compared directly at the same wave phase.

Because permeability and saturation both change, the result is an engineering
parameter-combination contrast, not an isolated causal estimate of phase lag.
All numerical metrics use unsmoothed pressure data.  A one-particle-layer
local affine reconstruction and Gaussian smoothing are used only for the
two-dimensional display.
"""

from __future__ import annotations

import argparse
import copy
import csv
import hashlib
import json
import math
import os
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.patches import Circle  # noqa: E402

import analyze_phase_conditioned_uplift as conditioned  # noqa: E402
from plot_pipeline_gradient_force import display_smoothed_grid  # noqa: E402


SCHEMA = "pipeline-phase-lag-parameter-contrast-v1"
# The prior same-parameter counterfactual selected this permeability before the
# present cross-parameter comparison.  It therefore avoids choosing the strong
# case after looking at the direct engineering contrast.  The lower-transmission
# 1e-13 case remains available through --strong-label as a sensitivity.
DEFAULT_LOWER_LABEL = "k7e-12_sw094_nosmooth"
DEFAULT_HIGHER_LABEL = "k9p79e-12_sw0993_nosmooth"
DISPLAY_SIGMA_DEFAULT = 3.0
SURFACE_RESPONSE_RELATIVE_TOLERANCE = 1.0e-3


class ParameterContrastError(RuntimeError):
    """Raised when the cross-parameter comparison is incomplete or invalid."""


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def file_record(path: Path) -> dict[str, Any]:
    path = path.resolve()
    if not path.is_file() or path.stat().st_size <= 0:
        raise ParameterContrastError(f"Missing or empty artifact: {path}")
    return {
        "path": str(path),
        "size_bytes": path.stat().st_size,
        "sha256": sha256(path),
    }


def finite(value: Any, name: str) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError) as error:
        raise ParameterContrastError(f"{name} is not numeric") from error
    if not math.isfinite(result):
        raise ParameterContrastError(f"{name} is not finite")
    return result


def relative_fraction(strong: float, weak: float) -> float | None:
    """Return strong/weak - 1, or null for a non-positive denominator."""

    strong = finite(strong, "strong metric")
    weak = finite(weak, "weak metric")
    if weak <= 0.0:
        return None
    result = strong / weak - 1.0
    if not math.isfinite(result):
        raise ParameterContrastError("Relative parameter contrast is non-finite")
    return result


def metric_contrast(strong: float, weak: float, units: str) -> dict[str, Any]:
    strong = finite(strong, "strong-lag metric")
    weak = finite(weak, "weak-lag metric")
    return {
        "units": units,
        "lower_k_lower_sw": strong,
        "higher_k_higher_sw": weak,
        "lower_minus_higher": strong - weak,
        "lower_minus_higher_fraction": relative_fraction(strong, weak),
    }


def _config_invariants(case: dict[str, Any]) -> dict[str, Any]:
    hd = case["hd"]
    water = hd["materials"][1]
    pipe = hd["analysis"]["rigid_pipeline"]
    return {
        "analysis_type": hd["analysis"]["type"],
        "dt_s": hd["analysis"]["dt"],
        "nsteps": hd["analysis"]["nsteps"],
        "APIC": hd["analysis"]["APIC"],
        "PIC": hd["analysis"]["PIC"],
        "pressure_smoothing": hd["analysis"]["pressure_smoothing"],
        "pressure_smoothing_in_loop": hd["analysis"][
            "pressure_smoothing_in_loop"
        ],
        "pipeline_fixed": pipe["fixed"],
        "pipeline_center": pipe["center"],
        "pipeline_outer_radius_m": pipe["outer_radius"],
        "wave_pressure": water["wave_pressure"],
        "wave_height_m": water["wave_height_ini_"],
        "wave_period_s": water["wave_period"],
        "wave_ramp_time_s": water["wave_ramp_time"],
        "domain_length_x_m": water["domain_length_x"],
        "depth_left_m": water["depth_left"],
        "depth_right_m": water["depth_right"],
        "gravity_m_s2": hd["external_loading_conditions"]["gravity"],
        "output_interval_steps": hd["post_processing"]["output_steps"],
    }


def _normalized_hd_config(case: dict[str, Any]) -> dict[str, Any]:
    """Remove only the two intended parameters and case-identity fields."""

    normalized = copy.deepcopy(case["hd"])
    normalized["title"] = "<case-title>"
    normalized["materials"][0]["intrinsic_permeability"] = "<varied-k>"
    normalized["materials"][1]["liquid_saturation"] = "<varied-Sw>"
    normalized["materials"][1]["gas_saturation"] = "<derived-from-Sw>"
    analysis = normalized["analysis"]
    analysis["uuid"] = "<case-uuid>"
    analysis["resume"]["uuid"] = "<equilibrium-uuid>"
    analysis["prescribed_phase_pressures"]["path"] = "<case-pressure-path>"
    normalized["post_processing"]["path"] = "<case-result-path>"
    return normalized


def validate_parameter_pair(
    strong: dict[str, Any], weak: dict[str, Any]
) -> dict[str, Any]:
    """Validate a declared lower-k/lower-Sw versus higher-k/higher-Sw pair."""

    strong_parameters = strong["validated_case"]
    weak_parameters = weak["validated_case"]
    strong_k = finite(
        strong_parameters["intrinsic_permeability_m2"], "strong permeability"
    )
    weak_k = finite(
        weak_parameters["intrinsic_permeability_m2"], "weak permeability"
    )
    strong_sw = finite(strong_parameters["liquid_saturation"], "strong saturation")
    weak_sw = finite(weak_parameters["liquid_saturation"], "weak saturation")
    if not strong_k < weak_k:
        raise ParameterContrastError(
            "Lower-k/lower-Sw case must have lower intrinsic permeability"
        )
    if not strong_sw < weak_sw:
        raise ParameterContrastError(
            "Lower-k/lower-Sw case must have lower liquid saturation"
        )
    for name, case in (("strong", strong), ("weak", weak)):
        parameters = case["validated_case"]
        if parameters["pressure_smoothing"] is not False:
            raise ParameterContrastError(f"{name} case uses pressure smoothing")
        if parameters["pipeline_fixed"] is not True:
            raise ParameterContrastError(f"{name} case does not keep the pipe fixed")

    strong_invariants = _config_invariants(strong)
    weak_invariants = _config_invariants(weak)
    if strong_invariants != weak_invariants:
        raise ParameterContrastError("Wave, geometry, solver, or output invariants differ")
    if _normalized_hd_config(strong) != _normalized_hd_config(weak):
        raise ParameterContrastError(
            "HD configurations differ beyond permeability, saturation, and case identity"
        )
    strong_reference = np.asarray(strong["reference"], dtype=np.float64)
    weak_reference = np.asarray(weak["reference"], dtype=np.float64)
    if strong_reference.shape != weak_reference.shape or not np.array_equal(
        strong_reference, weak_reference
    ):
        raise ParameterContrastError("Initial particle coordinates differ")
    strong_times = np.asarray([row["time_s"] for row in strong["rows"]])
    weak_times = np.asarray([row["time_s"] for row in weak["rows"]])
    if not np.array_equal(strong_times, weak_times):
        raise ParameterContrastError("Saved time grids differ")
    strong_surface = np.asarray(
        [row["surface_excess_pressure_pa"] for row in strong["rows"]]
    )
    weak_surface = np.asarray(
        [row["surface_excess_pressure_pa"] for row in weak["rows"]]
    )
    maximum_surface_difference = float(np.max(np.abs(strong_surface - weak_surface)))
    surface_response_scale = max(float(np.max(np.abs(weak_surface))), 1.0e-12)
    maximum_surface_difference_fraction = (
        maximum_surface_difference / surface_response_scale
    )
    if maximum_surface_difference_fraction > SURFACE_RESPONSE_RELATIVE_TOLERANCE:
        raise ParameterContrastError(
            "The compared cases have materially different surface-pressure responses"
        )
    for strong_branch, weak_branch in zip(
        strong["branch_rows"], weak["branch_rows"]
    ):
        if [row["time_s"] for row in strong_branch] != [
            row["time_s"] for row in weak_branch
        ]:
            raise ParameterContrastError("Negative-pressure recovery windows differ")
    return {
        "lower_k_lower_sw": {
            "label": strong["label"],
            "intrinsic_permeability_m2": strong_k,
            "liquid_saturation": strong_sw,
        },
        "higher_k_higher_sw": {
            "label": weak["label"],
            "intrinsic_permeability_m2": weak_k,
            "liquid_saturation": weak_sw,
        },
        "shared_invariants": strong_invariants,
        "maximum_surface_response_difference_pa": maximum_surface_difference,
        "maximum_surface_response_difference_fraction": (
            maximum_surface_difference_fraction
        ),
        "surface_response_relative_tolerance": (
            SURFACE_RESPONSE_RELATIVE_TOLERANCE
        ),
    }


def _branch_value(
    case: dict[str, Any],
    branch_index: int,
    metric: str,
    *,
    local_linear: bool = False,
) -> float:
    if metric in {
        "horizontal_absolute_force_activity_impulse",
        "vertical_absolute_force_activity_impulse",
    }:
        source = (
            case["secondary_local_linear_branches"]
            if local_linear
            else case["raw_directional_branches"]
        )
    else:
        source = case["branches"]
    return finite(source[branch_index][metric]["lagged"], metric)


def compare_recoveries(
    strong: dict[str, Any], weak: dict[str, Any]
) -> list[dict[str, Any]]:
    """Return two raw, matched-time engineering recovery contrasts."""

    output: list[dict[str, Any]] = []
    for branch_index in range(2):
        strong_branch = strong["branches"][branch_index]
        weak_branch = weak["branches"][branch_index]
        record: dict[str, Any] = {
            "recovery_index": branch_index + 1,
            "start_s": strong_branch["start_s"],
            "end_s": strong_branch["end_s"],
            "saved_frames": strong_branch["saved_frames"],
        }
        if (
            strong_branch["start_s"] != weak_branch["start_s"]
            or strong_branch["end_s"] != weak_branch["end_s"]
            or strong_branch["saved_frames"] != weak_branch["saved_frames"]
        ):
            raise ParameterContrastError("Recovery branch metadata differs")
        for metric, output_name in (
            ("local_upward_activity_impulse", "local_upward_activity"),
            ("signed_support_force_impulse", "signed_support_force"),
            (
                "horizontal_absolute_force_activity_impulse",
                "raw_horizontal_absolute_force_activity",
            ),
            (
                "vertical_absolute_force_activity_impulse",
                "raw_vertical_absolute_force_activity",
            ),
        ):
            record[output_name] = metric_contrast(
                _branch_value(strong, branch_index, metric),
                _branch_value(weak, branch_index, metric),
                "N s/m",
            )
        for metric, output_name in (
            (
                "horizontal_absolute_force_activity_impulse",
                "one_layer_horizontal_absolute_force_activity",
            ),
            (
                "vertical_absolute_force_activity_impulse",
                "one_layer_vertical_absolute_force_activity",
            ),
        ):
            record[output_name] = metric_contrast(
                _branch_value(
                    strong, branch_index, metric, local_linear=True
                ),
                _branch_value(weak, branch_index, metric, local_linear=True),
                "N s/m",
            )
        for prefix in ("raw", "one_layer"):
            horizontal = record[f"{prefix}_horizontal_absolute_force_activity"]
            vertical = record[f"{prefix}_vertical_absolute_force_activity"]
            strong_ratio = (
                vertical["lower_k_lower_sw"]
                / horizontal["lower_k_lower_sw"]
            )
            weak_ratio = (
                vertical["higher_k_higher_sw"]
                / horizontal["higher_k_higher_sw"]
            )
            record[f"{prefix}_vertical_to_horizontal_ratio"] = metric_contrast(
                strong_ratio, weak_ratio, "dimensionless"
            )
        output.append(record)
    return output


def _probe_record(case: dict[str, Any], name: str) -> dict[str, Any]:
    probe = case["phase_audit"]["probe_results"][name]
    surface = case["phase_audit"]["probe_results"]["surface"]
    phase = probe.get("phase_difference_from_same_column_surface_deg")
    amplitude = finite(probe["fundamental_amplitude_pa"], f"{name} amplitude")
    surface_amplitude = finite(
        surface["fundamental_amplitude_pa"], "surface amplitude"
    )
    return {
        "phase_difference_from_same_column_surface_deg": (
            finite(phase, f"{name} phase") if phase is not None else None
        ),
        "fundamental_amplitude_pa": amplitude,
        "surface_amplitude_ratio": amplitude / surface_amplitude,
        "harmonic_r_squared": finite(probe["harmonic_r_squared"], f"{name} R2"),
    }


def phase_contrast(strong: dict[str, Any], weak: dict[str, Any]) -> dict[str, Any]:
    probes: dict[str, Any] = {}
    for name in ("crown", "shoulder", "invert"):
        strong_probe = _probe_record(strong, name)
        weak_probe = _probe_record(weak, name)
        strong_phase = strong_probe["phase_difference_from_same_column_surface_deg"]
        weak_phase = weak_probe["phase_difference_from_same_column_surface_deg"]
        probes[name] = {
            "lower_k_lower_sw": strong_probe,
            "higher_k_higher_sw": weak_probe,
            "absolute_phase_difference_contrast_deg": (
                abs(strong_phase) - abs(weak_phase)
                if strong_phase is not None and weak_phase is not None
                else None
            ),
        }
    crown = probes["crown"]
    strong_crown = crown["lower_k_lower_sw"][
        "phase_difference_from_same_column_surface_deg"
    ]
    weak_crown = crown["higher_k_higher_sw"][
        "phase_difference_from_same_column_surface_deg"
    ]
    return {
        "probes": probes,
        "lower_k_lower_sw_has_larger_crown_phase_magnitude": (
            strong_crown is not None
            and weak_crown is not None
            and abs(strong_crown) > abs(weak_crown)
        ),
        "higher_k_higher_sw_is_exactly_zero_phase": (
            weak_crown is not None and abs(weak_crown) <= 1.0e-12
        ),
    }


def _case_label(case: dict[str, Any]) -> str:
    values = case["validated_case"]
    return (
        f"$k_s={values['intrinsic_permeability_m2']:.3g}$ m$^2$, "
        f"$S_w={values['liquid_saturation']:.3f}$"
    )


def atomic_savefig(figure: plt.Figure, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.stem}.{os.getpid()}.tmp{path.suffix}")
    try:
        figure.savefig(temporary, dpi=220, bbox_inches="tight", format=path.suffix[1:])
        os.replace(temporary, path)
    finally:
        if temporary.exists():
            temporary.unlink()


def plot_summary(
    strong: dict[str, Any], weak: dict[str, Any], recoveries: list[dict[str, Any]], output: Path
) -> tuple[Path, Path]:
    figure = plt.figure(figsize=(14.6, 7.4), constrained_layout=True)
    grid = figure.add_gridspec(2, 3, height_ratios=(1.2, 1.0))
    history_axis = figure.add_subplot(grid[0, :])
    local_axis = figure.add_subplot(grid[1, 0])
    signed_axis = figure.add_subplot(grid[1, 1])
    orientation_axis = figure.add_subplot(grid[1, 2])
    colours = ("#2F5597", "#D97941")
    for case, colour, name in (
        (strong, colours[0], "lower k / lower Sw"),
        (weak, colours[1], "higher k / higher Sw"),
    ):
        rows = [
            row
            for row in case["rows"]
            if conditioned.PRIMARY_CYCLE[0] - 1.0e-12
            <= row["time_s"]
            <= conditioned.PRIMARY_CYCLE[1] + 1.0e-12
        ]
        history_axis.plot(
            [row["time_s"] for row in rows],
            [row["lagged_positive_force_n_per_m"] for row in rows],
            color=colour,
            linewidth=2.2,
            label=f"{name}: {_case_label(case)}",
        )
    branch = strong["branch_rows"][-1]
    history_axis.axvspan(
        branch[0]["time_s"], branch[-1]["time_s"], color="#A8DADC", alpha=0.32
    )
    history_axis.set_xlim(*conditioned.PRIMARY_CYCLE)
    history_axis.set_xlabel("Time (s)")
    history_axis.set_ylabel("Raw local upward-force activity (N/m)")
    history_axis.set_title("(a) Physical pressure fields during the second post-ramp cycle")
    history_axis.grid(alpha=0.22)
    history_axis.legend(frameon=False, fontsize=8.5, loc="upper left")

    x = np.arange(2, dtype=float)
    width = 0.34
    recovery_labels = ("recovery 1", "recovery 2")
    for axis, key, title in (
        (local_axis, "local_upward_activity", "(b) Local upward-activity impulse"),
        (signed_axis, "signed_support_force", "(c) Signed support-force impulse"),
    ):
        axis.bar(
            x - width / 2,
            [row[key]["lower_k_lower_sw"] for row in recoveries],
            width,
            color=colours[0],
            label="lower k / lower Sw",
        )
        axis.bar(
            x + width / 2,
            [row[key]["higher_k_higher_sw"] for row in recoveries],
            width,
            color=colours[1],
            label="higher k / higher Sw",
        )
        axis.set_xticks(x, recovery_labels)
        axis.set_ylabel("Impulse (N s/m)")
        axis.set_title(title)
        axis.grid(axis="y", alpha=0.22)
        axis.legend(frameon=False, fontsize=8.0)

    primary = recoveries[-1]
    orientation_keys = (
        "raw_vertical_to_horizontal_ratio",
        "one_layer_vertical_to_horizontal_ratio",
    )
    orientation_axis.bar(
        x - width / 2,
        [primary[key]["lower_k_lower_sw"] for key in orientation_keys],
        width,
        color=colours[0],
        label="lower k / lower Sw",
    )
    orientation_axis.bar(
        x + width / 2,
        [primary[key]["higher_k_higher_sw"] for key in orientation_keys],
        width,
        color=colours[1],
        label="higher k / higher Sw",
    )
    orientation_axis.set_xticks(x, ("raw derivative", "one-layer affine"))
    orientation_axis.set_ylabel(r"Vertical/horizontal activity ratio $|f_y|/|f_x|$")
    orientation_axis.set_title("(d) Recovery-2 force orientation")
    orientation_axis.grid(axis="y", alpha=0.22)
    orientation_axis.legend(frameon=False, fontsize=8.0)
    figure.suptitle(
        "Direct engineering contrast: lower-k/lower-Sw versus higher-k/higher-Sw",
        fontsize=14,
    )
    figure.text(
        0.5,
        -0.012,
        "Raw unsmoothed metrics; changing k and Sw means this is not a "
        "phase-lag-only causal estimate.",
        ha="center",
        fontsize=9,
    )
    png = output / "phase_lag_parameter_contrast_summary.png"
    pdf = output / "phase_lag_parameter_contrast_summary.pdf"
    atomic_savefig(figure, png)
    atomic_savefig(figure, pdf)
    plt.close(figure)
    return png, pdf


def plot_2d(
    strong: dict[str, Any], weak: dict[str, Any], output: Path, display_sigma: float
) -> tuple[Path, Path, dict[str, Any]]:
    if not math.isfinite(display_sigma) or display_sigma < 0.0:
        raise ParameterContrastError("Display sigma must be finite and non-negative")
    fields = [
        conditioned._field_for_case(case, conditioned.COMMON_FIELD_TIME_S, local_linear=True)
        for case in (strong, weak)
    ]
    signed_fields: list[np.ndarray] = []
    for case, field in zip((strong, weak), fields):
        gradient = conditioned.local_linear_pressure_gradient(
            case["reference"],
            field["points"],
            case["lagged_pressure"][field["step"]],
            lattice_radius=1,
        )
        force = conditioned.upward_excess_force_density(
            gradient,
            field["arrays"]["liquid_densities"],
            np.asarray(case["hd"]["external_loading_conditions"]["gravity"]),
        )
        gamma = np.asarray(
            field["arrays"]["gamma_sub"], dtype=np.float64
        ).reshape(-1)
        signed_if = force / gamma
        if not np.all(np.isfinite(signed_if)):
            raise ParameterContrastError("Signed upward-force index is non-finite")
        signed_fields.append(signed_if)
    hd = strong["hd"]
    pipe = hd["analysis"]["rigid_pipeline"]
    center = np.asarray(pipe["center"], dtype=np.float64)
    diameter = 2.0 * float(pipe["outer_radius"])
    water = hd["materials"][1]
    bed_y = float(water["sea_level"]) - float(water["depth_left"])
    bounds = (
        center[0] - 1.5 * diameter,
        center[0] + 1.5 * diameter,
        bed_y - 2.0 * diameter,
        bed_y,
    )
    rasters = []
    for case, field, signed_if in zip((strong, weak), fields, signed_fields):
        rasters.append(
            display_smoothed_grid(
                field["points"],
                signed_if,
                case["hd"],
                center,
                bounds,
                display_sigma,
                nx=330,
                ny=230,
            )
        )
    for left, right in zip(rasters[0][:2], rasters[1][:2]):
        if not np.array_equal(left, right):
            raise ParameterContrastError("Two-dimensional raster grids differ")
    difference = rasters[0][2] - rasters[1][2]
    signed_values = np.concatenate(
        [grid[2][np.isfinite(grid[2])] for grid in rasters]
    )
    difference_values = difference[np.isfinite(difference)]
    signed_limit = max(
        float(np.quantile(np.abs(signed_values), 0.995)), 1.0e-12
    )
    difference_limit = max(
        float(np.quantile(np.abs(difference_values), 0.995)), 1.0e-12
    )
    figure, axes = plt.subplots(1, 3, figsize=(11.6, 3.6), constrained_layout=True)
    signed_mappable = None
    difference_mappable = None
    for axis, values, title in zip(
        axes,
        (rasters[0][2], rasters[1][2], difference),
        (
            "Lower k / lower Sw",
            "Higher k / higher Sw",
            "Lower-k/lower-Sw minus reference",
        ),
    ):
        if axis is axes[2]:
            mappable = axis.pcolormesh(
                rasters[0][0],
                rasters[0][1],
                values,
                cmap="RdBu_r",
                shading="auto",
                vmin=-difference_limit,
                vmax=difference_limit,
            )
            difference_mappable = mappable
        else:
            mappable = axis.pcolormesh(
                rasters[0][0],
                rasters[0][1],
                values,
                cmap="RdBu_r",
                shading="auto",
                vmin=-signed_limit,
                vmax=signed_limit,
            )
            signed_mappable = mappable
        axis.add_patch(
            Circle(
                tuple(center),
                float(pipe["outer_radius"]),
                facecolor="white",
                edgecolor="black",
                linewidth=1.0,
            )
        )
        axis.set_aspect("equal")
        axis.set_xlim(bounds[0], bounds[1])
        axis.set_ylim(bounds[2], bounds[3])
        axis.set_title(title)
        axis.set_xlabel("x (m)")
    axes[0].set_ylabel("y (m)")
    if signed_mappable is None or difference_mappable is None:
        raise ParameterContrastError("No parameter-contrast field was rendered")
    figure.colorbar(
        signed_mappable,
        ax=axes[:2],
        shrink=0.82,
        label="Signed upward hydraulic force index $f_y/\\gamma'$",
    )
    figure.colorbar(
        difference_mappable,
        ax=axes[2],
        shrink=0.82,
        label="Direct parameter-combination difference",
    )
    figure.suptitle(
        "Matched negative-pressure recovery at "
        f"t={conditioned.COMMON_FIELD_TIME_S:.3f} s",
        fontsize=13,
    )
    figure.text(
        0.5,
        -0.03,
        "One-layer affine reconstruction + display-only Gaussian "
        f"sigma={display_sigma:g} px; raw metrics are unsmoothed.",
        ha="center",
        fontsize=8.5,
    )
    png = output / "phase_lag_parameter_contrast_2d.png"
    pdf = output / "phase_lag_parameter_contrast_2d.pdf"
    atomic_savefig(figure, png)
    atomic_savefig(figure, pdf)
    plt.close(figure)
    return png, pdf, {
        "time_s": conditioned.COMMON_FIELD_TIME_S,
        "display_smoothing_sigma_pixels": display_sigma,
        "numeric_smoothing_applied": False,
        "spatial_gradient_reconstruction": {
            "method": "inverse-distance-weighted local affine pressure plane",
            "reference_lattice_radius": 1,
            "used_in_reported_metrics": False,
        },
        "signed_colour_limit": signed_limit,
        "difference_colour_limit": difference_limit,
        "source_vtp": [file_record(field["path"]) for field in fields],
    }


def atomic_write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    try:
        with temporary.open("x", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
    finally:
        if temporary.exists():
            temporary.unlink()


def atomic_write_json(path: Path, value: dict[str, Any]) -> None:
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


def summary_rows(recoveries: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for recovery in recoveries:
        row: dict[str, Any] = {
            "recovery_index": recovery["recovery_index"],
            "start_s": recovery["start_s"],
            "end_s": recovery["end_s"],
            "saved_frames": recovery["saved_frames"],
        }
        for key in (
            "local_upward_activity",
            "signed_support_force",
            "raw_horizontal_absolute_force_activity",
            "raw_vertical_absolute_force_activity",
            "raw_vertical_to_horizontal_ratio",
            "one_layer_horizontal_absolute_force_activity",
            "one_layer_vertical_absolute_force_activity",
            "one_layer_vertical_to_horizontal_ratio",
        ):
            record = recovery[key]
            row[f"{key}_lower_k_lower_sw"] = record["lower_k_lower_sw"]
            row[f"{key}_higher_k_higher_sw"] = record["higher_k_higher_sw"]
            row[f"{key}_lower_minus_higher"] = record["lower_minus_higher"]
            row[f"{key}_lower_minus_higher_percent"] = (
                None
                if record["lower_minus_higher_fraction"] is None
                else 100.0 * record["lower_minus_higher_fraction"]
            )
        rows.append(row)
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--lower-label",
        "--strong-label",
        dest="lower_label",
        default=DEFAULT_LOWER_LABEL,
    )
    parser.add_argument(
        "--higher-label",
        "--weak-label",
        dest="higher_label",
        default=DEFAULT_HIGHER_LABEL,
    )
    parser.add_argument(
        "--output-dir", default="analysis/phase_lag_parameter_contrast"
    )
    parser.add_argument(
        "--display-smoothing-sigma", type=float, default=DISPLAY_SIGMA_DEFAULT
    )
    args = parser.parse_args()
    root = Path(__file__).resolve().parent
    if args.lower_label == args.higher_label:
        raise ParameterContrastError("Lower and higher parameter labels must differ")
    strong = conditioned.load_case(
        root, args.lower_label, expected_saturation=None
    )
    weak = conditioned.load_case(root, args.higher_label, expected_saturation=None)
    pair = validate_parameter_pair(strong, weak)
    recoveries = compare_recoveries(strong, weak)
    phases = phase_contrast(strong, weak)
    output = (root / args.output_dir).resolve()
    output.mkdir(parents=True, exist_ok=True)
    csv_path = output / "phase_lag_parameter_contrast.csv"
    atomic_write_csv(csv_path, summary_rows(recoveries))
    summary_png, summary_pdf = plot_summary(strong, weak, recoveries, output)
    field_png, field_pdf, field = plot_2d(
        strong, weak, output, args.display_smoothing_sigma
    )
    audit_path = output / "phase_lag_parameter_contrast.audit.json"
    audit = {
        "schema": SCHEMA,
        "classification": "cross-parameter-engineering-contrast",
        "technical_question": (
            "How do raw phase, vertical hydraulic forcing, and force orientation "
            "differ between a lower-k/lower-Sw combination and a higher-k/higher-Sw reference?"
        ),
        "claim_boundary": (
            "Permeability and saturation both change. The contrast is useful for "
            "engineering bracketing but cannot attribute the full difference to phase lag alone."
        ),
        "parameter_pair": pair,
        "phase_contrast": phases,
        "negative_pressure_recoveries": recoveries,
        "field_snapshot": field,
        "direct_metric_pressure_role": "physical lagged HD pressure database",
        "phase_erased_values_used_in_direct_metrics": False,
        "numeric_pressure_smoothing": False,
        "display_smoothing_only": True,
        "sources": {
            "lower_k_lower_sw_runner": file_record(strong["runner_path"]),
            "lower_k_lower_sw_phase": file_record(strong["phase_path"]),
            "lower_k_lower_sw_driver": file_record(
                Path(strong["validated_case"]["sources"]["driver_screen"]["path"])
            ),
            "higher_k_higher_sw_runner": file_record(weak["runner_path"]),
            "higher_k_higher_sw_phase": file_record(weak["phase_path"]),
            "higher_k_higher_sw_driver": file_record(
                Path(weak["validated_case"]["sources"]["driver_screen"]["path"])
            ),
        },
        "implementation": {
            "analyzer": file_record(Path(__file__)),
            "conditioned_analyzer": file_record(Path(conditioned.__file__)),
        },
        "artifacts": {
            "csv": file_record(csv_path),
            "summary_png": file_record(summary_png),
            "summary_pdf": file_record(summary_pdf),
            "field_png": file_record(field_png),
            "field_pdf": file_record(field_pdf),
        },
    }
    atomic_write_json(audit_path, audit)
    print(audit_path)


if __name__ == "__main__":
    main()
