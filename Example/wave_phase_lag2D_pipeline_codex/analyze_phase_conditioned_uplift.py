#!/usr/bin/env python3
"""Audit short upward hydraulic forcing during surface-pressure recovery.

This analysis answers a deliberately narrower question than the existing
two-cycle driver screen: after the wave ramp, does retained subsurface phase
structure increase the *short-duration* upward hydraulic forcing while the
same-column seabed pressure is still negative and recovering from its trough?

The paired pressure fields are evaluated on the same fixed-pipeline HD state.
All derivatives and metrics use raw pressure-database values.  Gaussian
smoothing is permitted only while rasterising the two-dimensional figure.
The result is a fixed-state hydraulic-risk diagnostic, not a realised
liquefaction or released-pipeline response.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import sqlite3
from pathlib import Path
from typing import Any, Iterable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.patches import Circle  # noqa: E402

import plot_pipeline_gradient_force as gradient_force_module  # noqa: E402
import summarize_phase_lag_permeability_screen as permeability_summary  # noqa: E402
from analyze_study import (  # noqa: E402
    particle_step,
    pipeline_support_roi,
    read_initial_coordinates,
)
from plot_pipeline_gradient_force import (  # noqa: E402
    display_smoothed_grid,
    excess_force_density,
    local_linear_pressure_gradient,
    ordered_frame,
    raw_finite_difference_gradient,
    raw_force_metrics,
    read_raw_pressure_database,
    upward_excess_force_density,
)


SCHEMA = "pipeline-phase-conditioned-uplift-v1"
DEFAULT_LABELS = tuple(permeability_summary.DEFAULT_LABELS)
UPPER_BOUND_SENSITIVITY_LABELS = (
    "k5e-12_sw094_nosmooth",
    "k7e-12_sw094_nosmooth",
)
EXTENDED_LABELS = DEFAULT_LABELS + UPPER_BOUND_SENSITIVITY_LABELS
PERIOD_S = 1.3
RECOVERY_CYCLES = ((1.3, 2.6), (2.6, 3.9))
PRIMARY_CYCLE = RECOVERY_CYCLES[-1]
SURFACE_NEGATIVE_CORE_FRACTION = 0.10
MINIMUM_MATERIAL_EFFECT_FRACTION = 0.05
MINIMUM_RELEVANT_PHASE_DEG = 15.0
MINIMUM_RELEVANT_AMPLITUDE_RATIO = 0.05
MINIMUM_RELEVANT_R_SQUARED = 0.8
COMMON_FIELD_TIME_S = 3.835
EXPECTED_FRAME_COUNT = 60
EXPECTED_FRAME_DT_S = 0.065
DISPLAY_SIGMA_DEFAULT = 2.0


class ConditionedUpliftError(RuntimeError):
    """Raised when the short-time evidence is incomplete or inconsistent."""


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def file_record(path: Path) -> dict[str, Any]:
    path = path.resolve()
    if not path.is_file():
        raise ConditionedUpliftError(f"Missing artifact: {path}")
    return {
        "path": str(path),
        "size_bytes": path.stat().st_size,
        "sha256": sha256(path),
    }


def read_json(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ConditionedUpliftError(f"Cannot read JSON {path}") from error
    if not isinstance(value, dict):
        raise ConditionedUpliftError(f"JSON root must be an object: {path}")
    return value


def _finite(value: Any, name: str) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError) as error:
        raise ConditionedUpliftError(f"{name} is not numeric") from error
    if not math.isfinite(result):
        raise ConditionedUpliftError(f"{name} is not finite")
    return result


def relative_fraction(lagged: float, erased: float) -> float | None:
    lagged = _finite(lagged, "lagged metric")
    erased = _finite(erased, "phase-erased metric")
    if erased <= 0.0:
        return None
    result = lagged / erased - 1.0
    if not math.isfinite(result):
        raise ConditionedUpliftError("Relative effect is non-finite")
    return result


HISTORY_FIELDS = (
    "time_s",
    "lagged_signed_force_n_per_m",
    "lagged_net_uplift_force_n_per_m",
    "lagged_positive_force_n_per_m",
    "lagged_critical_area_m2",
    "lagged_maximum_upward_if",
    "lagged_joint_area_m2",
    "phase_erased_signed_force_n_per_m",
    "phase_erased_net_uplift_force_n_per_m",
    "phase_erased_positive_force_n_per_m",
    "phase_erased_critical_area_m2",
    "phase_erased_maximum_upward_if",
    "phase_erased_joint_area_m2",
    "delta_positive_force_n_per_m",
    "delta_signed_force_n_per_m",
    "delta_net_uplift_force_n_per_m",
)


def read_history(path: Path) -> list[dict[str, float]]:
    try:
        with path.open(newline="", encoding="utf-8") as stream:
            reader = csv.DictReader(stream)
            if tuple(reader.fieldnames or ()) != HISTORY_FIELDS:
                raise ConditionedUpliftError(
                    f"Unexpected driver-history columns in {path}"
                )
            rows = [
                {name: _finite(value, f"{path.name}:{name}") for name, value in row.items()}
                for row in reader
            ]
    except OSError as error:
        raise ConditionedUpliftError(f"Cannot read driver history {path}") from error
    if len(rows) != EXPECTED_FRAME_COUNT:
        raise ConditionedUpliftError(
            f"Driver history must contain {EXPECTED_FRAME_COUNT} frames"
        )
    times = np.asarray([row["time_s"] for row in rows], dtype=np.float64)
    if (
        not np.all(np.isfinite(times))
        or np.any(np.diff(times) <= 0.0)
        or not np.allclose(
            np.diff(times), EXPECTED_FRAME_DT_S, rtol=0.0, atol=1.0e-12
        )
    ):
        raise ConditionedUpliftError("Driver history time grid is incomplete")
    return rows


def attach_surface_pressure(
    rows: list[dict[str, float]],
    pressure_frames: dict[int, np.ndarray],
    *,
    surface_particle_id: int,
    dt_s: float,
) -> list[dict[str, float]]:
    """Attach same-column surface excess pressure to each saved driver row."""

    if 0 not in pressure_frames:
        raise ConditionedUpliftError("Pressure database has no step-zero baseline")
    baseline = _finite(
        pressure_frames[0][surface_particle_id], "surface step-zero pressure"
    )
    output: list[dict[str, float]] = []
    for row in rows:
        step = int(round(row["time_s"] / dt_s))
        if step not in pressure_frames:
            raise ConditionedUpliftError(f"Pressure database omits saved step {step}")
        pressure = _finite(
            pressure_frames[step][surface_particle_id],
            f"surface pressure at step {step}",
        )
        output.append(
            {
                **row,
                "surface_excess_pressure_pa": pressure - baseline,
                "step": float(step),
            }
        )
    return output


def recovery_branch(
    rows: list[dict[str, float]],
    *,
    cycle_start_s: float,
    cycle_end_s: float,
    surface_amplitude_pa: float,
    negative_core_fraction: float = SURFACE_NEGATIVE_CORE_FRACTION,
) -> list[dict[str, float]]:
    """Return the resolved trough-to-recovery branch inside one fixed cycle.

    The branch begins at the minimum saved same-column surface pressure and
    continues while that pressure remains at least ten percent of the fitted
    amplitude below equilibrium.  This excludes a numerical sign decision at
    the zero crossing and gives identical physical-phase samples in both
    post-ramp cycles.
    """

    values = (
        cycle_start_s,
        cycle_end_s,
        surface_amplitude_pa,
        negative_core_fraction,
    )
    if not all(math.isfinite(value) for value in values):
        raise ConditionedUpliftError("Recovery-branch parameters are non-finite")
    if cycle_end_s <= cycle_start_s or surface_amplitude_pa <= 0.0:
        raise ConditionedUpliftError("Recovery-branch parameters are invalid")
    if not 0.0 < negative_core_fraction < 1.0:
        raise ConditionedUpliftError("Negative-core fraction must lie in (0,1)")
    cycle = [
        row
        for row in rows
        if cycle_start_s - 1.0e-12
        <= row["time_s"]
        <= cycle_end_s + 1.0e-12
    ]
    if len(cycle) < 3:
        raise ConditionedUpliftError("Cycle has too few saved frames")
    pressure = np.asarray(
        [row["surface_excess_pressure_pa"] for row in cycle], dtype=np.float64
    )
    minimum_index = int(np.argmin(pressure))
    threshold = -negative_core_fraction * surface_amplitude_pa
    branch = [
        row
        for row in cycle[minimum_index:]
        if row["surface_excess_pressure_pa"] <= threshold + 1.0e-12
    ]
    if len(branch) < 3 or branch[0] is not cycle[minimum_index]:
        raise ConditionedUpliftError("Negative-pressure recovery is unresolved")
    branch_times = np.asarray([row["time_s"] for row in branch], dtype=np.float64)
    branch_pressure = np.asarray(
        [row["surface_excess_pressure_pa"] for row in branch], dtype=np.float64
    )
    if (
        not np.allclose(
            np.diff(branch_times), EXPECTED_FRAME_DT_S, rtol=0.0, atol=1.0e-12
        )
        or np.any(np.diff(branch_pressure) <= 0.0)
        or not math.isclose(branch_times[-1], cycle_end_s, abs_tol=1.0e-12)
    ):
        raise ConditionedUpliftError(
            "Negative-pressure recovery must be contiguous, rising, and reach the cycle end"
        )
    return branch


def integrate_branch(rows: list[dict[str, float]]) -> dict[str, Any]:
    """Integrate paired short-time metrics on one recovery branch."""

    if len(rows) < 3:
        raise ConditionedUpliftError("Recovery branch has too few frames")
    times = np.asarray([row["time_s"] for row in rows], dtype=np.float64)
    output: dict[str, Any] = {
        "start_s": float(times[0]),
        "end_s": float(times[-1]),
        "saved_frames": len(rows),
        "duration_s": float(times[-1] - times[0]),
        "surface_excess_pressure_range_pa": [
            float(min(row["surface_excess_pressure_pa"] for row in rows)),
            float(max(row["surface_excess_pressure_pa"] for row in rows)),
        ],
    }
    metrics = (
        ("positive_force_n_per_m", "local_upward_activity_impulse", "N s/m"),
        ("signed_force_n_per_m", "signed_support_force_impulse", "N s/m"),
        ("net_uplift_force_n_per_m", "net_uplift_impulse", "N s/m"),
        ("critical_area_m2", "IF_ge_1_area_time", "m2 s"),
        ("maximum_upward_if", "maximum_upward_IF_time_integral", "s"),
    )
    if all(
        f"{role}_{source}" in row
        for row in rows
        for role in ("lagged", "phase_erased")
        for source in (
            "horizontal_absolute_force_n_per_m",
            "vertical_absolute_force_n_per_m",
        )
    ):
        metrics += (
            (
                "horizontal_absolute_force_n_per_m",
                "horizontal_absolute_force_activity_impulse",
                "N s/m",
            ),
            (
                "vertical_absolute_force_n_per_m",
                "vertical_absolute_force_activity_impulse",
                "N s/m",
            ),
        )
    for source, name, units in metrics:
        lagged = np.asarray(
            [row[f"lagged_{source}"] for row in rows], dtype=np.float64
        )
        erased = np.asarray(
            [row[f"phase_erased_{source}"] for row in rows], dtype=np.float64
        )
        lagged_integral = float(np.trapz(lagged, times))
        erased_integral = float(np.trapz(erased, times))
        paired_delta = lagged - erased
        output[name] = {
            "units": units,
            "lagged": lagged_integral,
            "phase_erased": erased_integral,
            "lagged_minus_phase_erased": lagged_integral - erased_integral,
            "lagged_minus_phase_erased_fraction": relative_fraction(
                lagged_integral, erased_integral
            ),
            "lagged_peak": float(np.max(lagged)),
            "phase_erased_peak": float(np.max(erased)),
            "peak_to_peak_fraction": relative_fraction(
                float(np.max(lagged)), float(np.max(erased))
            ),
            "maximum_paired_delta": float(np.max(paired_delta)),
            "maximum_paired_delta_fraction": float(
                np.max(
                    np.divide(
                        paired_delta,
                        erased,
                        out=np.full_like(paired_delta, -math.inf),
                        where=erased > 0.0,
                    )
                )
            ),
        }
    if "horizontal_absolute_force_activity_impulse" in output:
        horizontal = output["horizontal_absolute_force_activity_impulse"]
        vertical = output["vertical_absolute_force_activity_impulse"]
        if horizontal["lagged"] <= 0.0 or horizontal["phase_erased"] <= 0.0:
            raise ConditionedUpliftError(
                "Horizontal force-activity impulse must be positive"
            )
        lagged_ratio = vertical["lagged"] / horizontal["lagged"]
        erased_ratio = vertical["phase_erased"] / horizontal["phase_erased"]
        output["vertical_to_horizontal_force_activity_ratio"] = {
            "lagged": lagged_ratio,
            "phase_erased": erased_ratio,
            "lagged_minus_phase_erased": lagged_ratio - erased_ratio,
            "lagged_minus_phase_erased_fraction": relative_fraction(
                lagged_ratio, erased_ratio
            ),
        }
    return output


def relevant_pipe_phase(phase_audit: dict[str, Any]) -> dict[str, Any]:
    """Return eligible phase probes around the pipe, not just its crown."""

    try:
        surface_amplitude = _finite(
            phase_audit["probe_results"]["surface"]["fundamental_amplitude_pa"],
            "surface fundamental amplitude",
        )
        probes = phase_audit["probe_results"]
    except (KeyError, TypeError) as error:
        raise ConditionedUpliftError("Phase audit lacks probe results") from error
    eligible: list[dict[str, Any]] = []
    for name in ("crown", "shoulder", "invert"):
        probe = probes.get(name)
        if not isinstance(probe, dict):
            raise ConditionedUpliftError(f"Phase audit lacks {name} probe")
        phase_value = probe.get("phase_difference_from_same_column_surface_deg")
        if phase_value is None:
            continue
        phase = _finite(phase_value, f"{name} phase")
        amplitude = _finite(probe.get("fundamental_amplitude_pa"), f"{name} amplitude")
        r_squared = _finite(probe.get("harmonic_r_squared"), f"{name} R2")
        amplitude_ratio = amplitude / surface_amplitude
        if (
            r_squared >= MINIMUM_RELEVANT_R_SQUARED
            and amplitude_ratio >= MINIMUM_RELEVANT_AMPLITUDE_RATIO
            and abs(phase) >= MINIMUM_RELEVANT_PHASE_DEG
        ):
            eligible.append(
                {
                    "probe": name,
                    "same_column_phase_difference_deg": phase,
                    "absolute_phase_difference_deg": abs(phase),
                    "amplitude_ratio": amplitude_ratio,
                    "harmonic_r_squared": r_squared,
                }
            )
    eligible.sort(key=lambda item: item["absolute_phase_difference_deg"], reverse=True)
    return {
        "passed": bool(eligible),
        "criteria": {
            "minimum_absolute_phase_difference_deg": MINIMUM_RELEVANT_PHASE_DEG,
            "minimum_amplitude_ratio": MINIMUM_RELEVANT_AMPLITUDE_RATIO,
            "minimum_harmonic_r_squared": MINIMUM_RELEVANT_R_SQUARED,
        },
        "eligible_probes": eligible,
        "most_delayed_eligible_probe": eligible[0] if eligible else None,
    }


def short_time_gate(
    branches: list[dict[str, Any]], phase_qualification: dict[str, Any]
) -> dict[str, Any]:
    """Evaluate a mechanism gate without claiming realised liquefaction."""

    if len(branches) != 2:
        raise ConditionedUpliftError("Exactly two post-ramp recoveries are required")
    local = [branch["local_upward_activity_impulse"] for branch in branches]
    signed = [branch["signed_support_force_impulse"] for branch in branches]
    checks = {
        "relevant_pipe_phase_is_resolved": bool(phase_qualification["passed"]),
        "both_recoveries_have_positive_local_upward_activity_delta": all(
            item["lagged_minus_phase_erased"] > 0.0 for item in local
        ),
        "both_recoveries_increase_local_upward_activity_by_at_least_5pct": all(
            item["lagged_minus_phase_erased_fraction"] is not None
            and item["lagged_minus_phase_erased_fraction"]
            >= MINIMUM_MATERIAL_EFFECT_FRACTION
            for item in local
        ),
        "both_recoveries_have_positive_signed_support_force_delta": all(
            item["lagged_minus_phase_erased"] > 0.0 for item in signed
        ),
        "primary_recovery_has_at_least_5pct_paired_local_peak_advantage": (
            local[-1]["maximum_paired_delta_fraction"]
            >= MINIMUM_MATERIAL_EFFECT_FRACTION
        ),
    }
    return {
        "passed": all(checks.values()),
        "classification": (
            "short-time-hydraulic-mechanism-supported"
            if all(checks.values())
            else "short-time-hydraulic-mechanism-not-resolved"
        ),
        "checks": checks,
        "minimum_material_effect_fraction": MINIMUM_MATERIAL_EFFECT_FRACTION,
        "claim_boundary": (
            "A PASS supports a fixed-state short-time upward hydraulic-risk "
            "mechanism only; it is not realised liquefaction or pipeline uplift."
        ),
    }


def _history_path(case: dict[str, Any]) -> Path:
    try:
        return Path(
            case["sources"]["driver_provenance"]["history_csv"]["path"]
        ).resolve()
    except (KeyError, TypeError) as error:
        raise ConditionedUpliftError("Validated case lacks driver history") from error


def load_case(
    root: Path,
    label: str,
    *,
    expected_saturation: float | None = permeability_summary.EXPECTED_SATURATION,
) -> dict[str, Any]:
    """Load one fully audited case and calculate its recovery metrics."""

    try:
        validated = permeability_summary.load_case(
            root, label, expected_saturation=expected_saturation
        )
    except Exception as error:  # preserve the detailed underlying cause
        raise ConditionedUpliftError(f"{label}: source audits are invalid") from error
    runner_path = Path(validated["runner_audit_path"]).resolve()
    runner_audit = read_json(runner_path)
    phase_path = Path(validated["sources"]["phase_v3"]["path"]).resolve()
    phase_audit = read_json(phase_path)
    try:
        hd_path = Path(runner_audit["configs"]["03_HD.json"]["path"]).resolve()
        rl_path = Path(runner_audit["configs"]["04_RL.json"]["path"]).resolve()
        re_path = Path(runner_audit["configs"]["05_RE.json"]["path"]).resolve()
    except (KeyError, TypeError) as error:
        raise ConditionedUpliftError(f"{label}: runner configs are incomplete") from error
    hd = read_json(hd_path)
    rl = read_json(rl_path)
    re_case = read_json(re_path)
    reference = read_initial_coordinates(
        root / hd["particles"][0]["generator"]["location"]
    )
    lagged_pressure, lagged_database = read_raw_pressure_database(root, rl, reference)
    erased_pressure, erased_database = read_raw_pressure_database(root, re_case, reference)
    try:
        surface_particle_id = int(
            phase_audit["probe_results"]["surface"]["particle_id"]
        )
        surface_amplitude = _finite(
            phase_audit["probe_results"]["surface"]["fundamental_amplitude_pa"],
            "surface amplitude",
        )
    except (KeyError, TypeError, ValueError) as error:
        raise ConditionedUpliftError(f"{label}: surface probe is invalid") from error
    if not 0 <= surface_particle_id < len(reference):
        raise ConditionedUpliftError(f"{label}: surface particle ID is out of range")
    common_steps = sorted(set(lagged_pressure) & set(erased_pressure))
    if common_steps != sorted(lagged_pressure) or common_steps != sorted(erased_pressure):
        raise ConditionedUpliftError(f"{label}: pressure time grids differ")
    surface_difference = max(
        abs(
            float(lagged_pressure[step][surface_particle_id])
            - float(erased_pressure[step][surface_particle_id])
        )
        for step in common_steps
    )
    if surface_difference > 1.0e-10:
        raise ConditionedUpliftError(
            f"{label}: phase erasure changed the registered surface driver"
        )
    history_path = _history_path(validated)
    rows = attach_surface_pressure(
        read_history(history_path),
        lagged_pressure,
        surface_particle_id=surface_particle_id,
        dt_s=float(hd["analysis"]["dt"]),
    )
    branch_rows = [
        recovery_branch(
            rows,
            cycle_start_s=start,
            cycle_end_s=end,
            surface_amplitude_pa=surface_amplitude,
        )
        for start, end in RECOVERY_CYCLES
    ]
    branches = [integrate_branch(branch) for branch in branch_rows]
    phase_qualification = relevant_pipe_phase(phase_audit)
    gate = short_time_gate(branches, phase_qualification)
    case = {
        "label": label,
        "validated_case": validated,
        "runner_audit": runner_audit,
        "runner_path": runner_path,
        "phase_audit": phase_audit,
        "phase_path": phase_path,
        "history_path": history_path,
        "hd_path": hd_path,
        "rl_path": rl_path,
        "re_path": re_path,
        "hd": hd,
        "rl": rl,
        "re": re_case,
        "reference": reference,
        "lagged_pressure": lagged_pressure,
        "erased_pressure": erased_pressure,
        "lagged_database": lagged_database,
        "erased_database": erased_database,
        "rows": rows,
        "branch_rows": branch_rows,
        "branches": branches,
        "surface_particle_id": surface_particle_id,
        "surface_amplitude_pa": surface_amplitude,
        "surface_driver_maximum_pair_difference_pa": surface_difference,
        "phase_qualification": phase_qualification,
        "gate": gate,
    }
    case["raw_directional_branches"] = spatial_recovery_diagnostics(
        case, local_linear=False
    )
    case["secondary_local_linear_branches"] = spatial_recovery_diagnostics(
        case, local_linear=True
    )
    return case


def spatial_recovery_diagnostics(
    case: dict[str, Any], *, local_linear: bool
) -> list[dict[str, Any]]:
    """Re-evaluate recovery metrics and pressure-force orientation.

    The raw route independently reproduces the registered upward metrics and
    adds horizontal/vertical force activity.  The local-linear route is a
    deliberately secondary one-layer reconstruction.  It changes neither the
    raw history, the registered recovery windows nor the short-time gate.
    """

    hd = case["hd"]
    support = pipeline_support_roi(hd, case["reference"])
    if support is None or not np.any(support):
        raise ConditionedUpliftError("Local-linear support cohort is empty")
    records = case["runner_audit"]["stages"]["hd"]["audit"]["particle_vtp"]
    by_step = {
        particle_step(Path(record["path"])): Path(record["path"])
        for record in records
    }
    gravity = np.asarray(hd["external_loading_conditions"]["gravity"], dtype=float)
    output: list[dict[str, Any]] = []
    for branch in case["branch_rows"]:
        reconstructed_rows: list[dict[str, float]] = []
        for raw_row in branch:
            step = int(round(raw_row["time_s"] / float(hd["analysis"]["dt"])))
            path = by_step.get(step)
            if path is None:
                raise ConditionedUpliftError(
                    f"{case['label']}: missing VTP for local-linear step {step}"
                )
            points, arrays = ordered_frame(path, len(case["reference"]))
            row: dict[str, float] = {
                "time_s": raw_row["time_s"],
                "surface_excess_pressure_pa": raw_row[
                    "surface_excess_pressure_pa"
                ],
            }
            for role, pressure in (
                ("lagged", case["lagged_pressure"][step]),
                ("phase_erased", case["erased_pressure"][step]),
            ):
                if local_linear:
                    gradient = local_linear_pressure_gradient(
                        case["reference"], points, pressure, lattice_radius=1
                    )
                else:
                    gradient = raw_finite_difference_gradient(
                        case["reference"], points, pressure
                    )
                vector_force = excess_force_density(
                    gradient, arrays["liquid_densities"], gravity
                )
                force = upward_excess_force_density(
                    gradient, arrays["liquid_densities"], gravity
                )
                metrics, _ = raw_force_metrics(
                    force,
                    arrays["gamma_sub"],
                    arrays["volumes"],
                    support,
                )
                row[f"{role}_positive_force_n_per_m"] = (
                    metrics.positive_force_n_per_m
                )
                row[f"{role}_signed_force_n_per_m"] = metrics.signed_force_n_per_m
                row[f"{role}_net_uplift_force_n_per_m"] = max(
                    metrics.signed_force_n_per_m, 0.0
                )
                row[f"{role}_critical_area_m2"] = metrics.critical_area_m2
                row[f"{role}_maximum_upward_if"] = metrics.maximum_upward_if
                volumes = np.asarray(
                    arrays["volumes"], dtype=np.float64
                ).reshape(-1)
                row[f"{role}_horizontal_absolute_force_n_per_m"] = float(
                    np.sum(np.abs(vector_force[support, 0]) * volumes[support])
                )
                row[f"{role}_vertical_absolute_force_n_per_m"] = float(
                    np.sum(np.abs(force[support]) * volumes[support])
                )
                if not local_linear:
                    for observed, expected, name in (
                        (
                            metrics.positive_force_n_per_m,
                            raw_row[f"{role}_positive_force_n_per_m"],
                            "positive force",
                        ),
                        (
                            metrics.signed_force_n_per_m,
                            raw_row[f"{role}_signed_force_n_per_m"],
                            "signed force",
                        ),
                    ):
                        if not math.isclose(
                            observed, expected, rel_tol=1.0e-10, abs_tol=1.0e-10
                        ):
                            raise ConditionedUpliftError(
                                f"{case['label']}: raw {name} reconstruction drift"
                            )
            reconstructed_rows.append(row)
        output.append(integrate_branch(reconstructed_rows))
    return output


def _case_colour(index: int) -> str:
    return ("#2F5597", "#008C95", "#E07B39", "#A23B72", "#6B8E23", "#6A5ACD")[
        index
    ]


def _label_k(case: dict[str, Any]) -> str:
    permeability = case["validated_case"]["intrinsic_permeability_m2"]
    return f"$k_s={permeability:.0e}$ m$^2$"


def atomic_savefig(figure: plt.Figure, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.stem}.{os.getpid()}.tmp{path.suffix}")
    try:
        figure.savefig(temporary, dpi=220, bbox_inches="tight", format=path.suffix[1:])
        os.replace(temporary, path)
    finally:
        if temporary.exists():
            temporary.unlink()


def plot_summary(cases: list[dict[str, Any]], output: Path) -> tuple[Path, Path]:
    figure = plt.figure(figsize=(14.6, 7.4), constrained_layout=True)
    grid = figure.add_gridspec(2, 3, height_ratios=(1.2, 1.0))
    history_axis = figure.add_subplot(grid[0, :])
    local_axis = figure.add_subplot(grid[1, 0])
    signed_axis = figure.add_subplot(grid[1, 1])
    direction_axis = figure.add_subplot(grid[1, 2])

    for index, case in enumerate(cases):
        rows = [
            row
            for row in case["rows"]
            if PRIMARY_CYCLE[0] - 1.0e-12
            <= row["time_s"]
            <= PRIMARY_CYCLE[1] + 1.0e-12
        ]
        history_axis.plot(
            [row["time_s"] for row in rows],
            [row["delta_positive_force_n_per_m"] for row in rows],
            color=_case_colour(index),
            linewidth=2.0,
            label=_label_k(case),
        )
    primary = cases[0]["branch_rows"][-1]
    history_axis.axvspan(
        primary[0]["time_s"],
        primary[-1]["time_s"],
        color="#A8DADC",
        alpha=0.32,
        label="negative-pressure recovery",
    )
    history_axis.axvline(COMMON_FIELD_TIME_S, color="0.25", linestyle="--", linewidth=1.2)
    history_axis.axhline(0.0, color="0.35", linewidth=0.9)
    history_axis.set_xlim(*PRIMARY_CYCLE)
    history_axis.set_ylabel("Lagged − erased local upward activity (N/m)")
    history_axis.set_xlabel("Time (s)")
    history_axis.set_title("(a) Phase-conditioned short-time upward forcing")
    history_axis.grid(alpha=0.22)
    history_axis.legend(ncol=3, fontsize=8.5, frameon=False, loc="upper left")

    x = np.arange(len(cases), dtype=float)
    width = 0.34
    labels = [f"{case['validated_case']['intrinsic_permeability_m2']:.0e}" for case in cases]
    first_local = [
        100.0
        * case["branches"][0]["local_upward_activity_impulse"][
            "lagged_minus_phase_erased_fraction"
        ]
        for case in cases
    ]
    second_local = [
        100.0
        * case["branches"][1]["local_upward_activity_impulse"][
            "lagged_minus_phase_erased_fraction"
        ]
        for case in cases
    ]
    first_signed = [
        100.0
        * case["branches"][0]["signed_support_force_impulse"][
            "lagged_minus_phase_erased_fraction"
        ]
        for case in cases
    ]
    second_signed = [
        100.0
        * case["branches"][1]["signed_support_force_impulse"][
            "lagged_minus_phase_erased_fraction"
        ]
        for case in cases
    ]
    for axis, first, second, title, ylabel in (
        (
            local_axis,
            first_local,
            second_local,
            "(b) Local upward-activity impulse",
            "Lagged − erased (%)",
        ),
        (
            signed_axis,
            first_signed,
            second_signed,
            "(c) Signed support-force impulse",
            "Lagged − erased (%)",
        ),
    ):
        axis.bar(x - width / 2, first, width, color="#8FB9D0", label="recovery 1")
        axis.bar(x + width / 2, second, width, color="#D95F59", label="recovery 2")
        axis.axhline(0.0, color="0.25", linewidth=0.9)
        if axis is local_axis:
            axis.axhline(5.0, color="#4A7C59", linestyle="--", linewidth=1.1)
        axis.set_xticks(x, labels)
        axis.set_xlabel("Intrinsic permeability $k_s$ (m$^2$)")
        axis.set_ylabel(ylabel)
        axis.set_title(title)
        axis.grid(axis="y", alpha=0.22)
        axis.legend(frameon=False, fontsize=8.5)
    raw_horizontal = [
        100.0
        * case["raw_directional_branches"][-1][
            "horizontal_absolute_force_activity_impulse"
        ]["lagged_minus_phase_erased_fraction"]
        for case in cases
    ]
    raw_vertical = [
        100.0
        * case["raw_directional_branches"][-1][
            "vertical_absolute_force_activity_impulse"
        ]["lagged_minus_phase_erased_fraction"]
        for case in cases
    ]
    local_horizontal = [
        100.0
        * case["secondary_local_linear_branches"][-1][
            "horizontal_absolute_force_activity_impulse"
        ]["lagged_minus_phase_erased_fraction"]
        for case in cases
    ]
    local_vertical = [
        100.0
        * case["secondary_local_linear_branches"][-1][
            "vertical_absolute_force_activity_impulse"
        ]["lagged_minus_phase_erased_fraction"]
        for case in cases
    ]
    direction_axis.bar(
        x - width / 2,
        raw_horizontal,
        width,
        color="#8093A8",
        label=r"raw $|f_x|$",
    )
    direction_axis.bar(
        x + width / 2,
        raw_vertical,
        width,
        color="#D95F59",
        label=r"raw $|f_y|$",
    )
    direction_axis.scatter(
        x - width / 2,
        local_horizontal,
        marker="D",
        s=28,
        facecolor="white",
        edgecolor="#26384A",
        zorder=4,
        label="one-layer check",
    )
    direction_axis.scatter(
        x + width / 2,
        local_vertical,
        marker="D",
        s=28,
        facecolor="white",
        edgecolor="#7B1E1E",
        zorder=4,
    )
    direction_axis.axhline(0.0, color="0.25", linewidth=0.9)
    direction_axis.set_xticks(x, labels)
    direction_axis.set_xlabel("Intrinsic permeability $k_s$ (m$^2$)")
    direction_axis.set_ylabel("Lagged − erased (%)")
    direction_axis.set_title("(d) Recovery-2 force orientation")
    direction_axis.grid(axis="y", alpha=0.22)
    direction_axis.legend(frameon=False, fontsize=8.0)
    figure.suptitle(
        "Negative-boundary-pressure recovery: raw fixed-state hydraulic comparison",
        fontsize=14,
    )
    figure.text(
        0.5,
        -0.012,
        "Metrics use unsmoothed pressure data; the shaded branch starts at the surface-pressure trough.",
        ha="center",
        fontsize=9,
    )
    png = output / "phase_conditioned_uplift_summary.png"
    pdf = output / "phase_conditioned_uplift_summary.pdf"
    atomic_savefig(figure, png)
    atomic_savefig(figure, pdf)
    plt.close(figure)
    return png, pdf


def _field_for_case(
    case: dict[str, Any], time_s: float, *, local_linear: bool = False
) -> dict[str, Any]:
    hd = case["hd"]
    dt = float(hd["analysis"]["dt"])
    step = int(round(time_s / dt))
    records = case["runner_audit"]["stages"]["hd"]["audit"]["particle_vtp"]
    matches = [Path(record["path"]) for record in records if particle_step(Path(record["path"])) == step]
    if len(matches) != 1:
        raise ConditionedUpliftError(
            f"{case['label']}: expected one particle VTP for step {step}"
        )
    points, arrays = ordered_frame(matches[0], len(case["reference"]))
    gravity = np.asarray(hd["external_loading_conditions"]["gravity"], dtype=float)
    fields: dict[str, np.ndarray] = {}
    for role, pressure in (
        ("lagged", case["lagged_pressure"][step]),
        ("phase_erased", case["erased_pressure"][step]),
    ):
        if local_linear:
            gradient = local_linear_pressure_gradient(
                case["reference"], points, pressure, lattice_radius=1
            )
        else:
            gradient = raw_finite_difference_gradient(
                case["reference"], points, pressure
            )
        force = upward_excess_force_density(
            gradient, arrays["liquid_densities"], gravity
        )
        gamma = np.asarray(arrays["gamma_sub"], dtype=np.float64).reshape(-1)
        signed_if = force / gamma
        fields[role] = np.maximum(signed_if, 0.0)
    fields["difference"] = fields["lagged"] - fields["phase_erased"]
    return {
        "step": step,
        "time_s": time_s,
        "points": points,
        "arrays": arrays,
        "fields": fields,
        "path": matches[0],
    }


def plot_2d(
    cases: list[dict[str, Any]], output: Path, display_sigma: float
) -> tuple[Path, Path, dict[str, Any]]:
    if not math.isfinite(display_sigma) or display_sigma < 0.0:
        raise ConditionedUpliftError("Display sigma must be finite and non-negative")
    reconstructed_fields = [
        _field_for_case(case, COMMON_FIELD_TIME_S, local_linear=True)
        for case in cases
    ]
    rasters: list[dict[str, Any]] = []
    all_positive: list[np.ndarray] = []
    all_difference: list[np.ndarray] = []
    for case, field in zip(cases, reconstructed_fields):
        pipeline = case["hd"]["analysis"]["rigid_pipeline"]
        center = np.asarray(pipeline["center"], dtype=float)
        diameter = 2.0 * float(pipeline["outer_radius"])
        bed_y = float(case["hd"]["materials"][1]["sea_level"]) - float(
            case["hd"]["materials"][1]["depth_left"]
        )
        bounds = (
            center[0] - 1.5 * diameter,
            center[0] + 1.5 * diameter,
            bed_y - 2.0 * diameter,
            bed_y,
        )
        grids: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]] = {}
        for role in ("lagged", "phase_erased", "difference"):
            grids[role] = display_smoothed_grid(
                field["points"],
                field["fields"][role],
                case["hd"],
                center,
                bounds,
                display_sigma,
                nx=330,
                ny=230,
            )
            finite = grids[role][2][np.isfinite(grids[role][2])]
            if role == "difference":
                all_difference.append(finite)
            else:
                all_positive.append(finite)
        rasters.append({"bounds": bounds, "center": center, "grids": grids})
    positive_limit = float(np.quantile(np.concatenate(all_positive), 0.995))
    difference_limit = float(
        np.quantile(np.abs(np.concatenate(all_difference)), 0.995)
    )
    positive_limit = max(positive_limit, 1.0e-12)
    difference_limit = max(difference_limit, 1.0e-12)

    figure, axes = plt.subplots(
        len(cases),
        3,
        figsize=(11.3, 2.55 * len(cases) + 1.4),
        sharex=True,
        sharey=True,
        constrained_layout=True,
    )
    positive_mappable = None
    difference_mappable = None
    for row, (case, raster) in enumerate(zip(cases, rasters)):
        for column, role in enumerate(("lagged", "phase_erased", "difference")):
            axis = axes[row, column]
            xx, yy, values = raster["grids"][role]
            if role == "difference":
                mappable = axis.pcolormesh(
                    xx,
                    yy,
                    values,
                    shading="auto",
                    cmap="RdBu_r",
                    vmin=-difference_limit,
                    vmax=difference_limit,
                )
                difference_mappable = mappable
            else:
                mappable = axis.pcolormesh(
                    xx,
                    yy,
                    values,
                    shading="auto",
                    cmap="viridis",
                    vmin=0.0,
                    vmax=positive_limit,
                )
                positive_mappable = mappable
            pipeline = case["hd"]["analysis"]["rigid_pipeline"]
            axis.add_patch(
                Circle(
                    tuple(raster["center"]),
                    float(pipeline["outer_radius"]),
                    facecolor="white",
                    edgecolor="black",
                    linewidth=1.0,
                    zorder=5,
                )
            )
            axis.set_aspect("equal")
            axis.set_xlim(raster["bounds"][0], raster["bounds"][1])
            axis.set_ylim(raster["bounds"][2], raster["bounds"][3])
            if row == 0:
                axis.set_title(
                    ("Lagged $I_F^+$", "Phase-erased $I_F^+$", "$\Delta I_F^+$")[column]
                )
            if column == 0:
                effect = 100.0 * case["branches"][-1]["local_upward_activity_impulse"][
                    "lagged_minus_phase_erased_fraction"
                ]
                axis.set_ylabel(f"{_label_k(case)}\ny (m)\nrecovery $\Delta J_+={effect:+.2f}\%$")
            if row == len(cases) - 1:
                axis.set_xlabel("x (m)")
    if positive_mappable is None or difference_mappable is None:
        raise ConditionedUpliftError("No two-dimensional field was rendered")
    figure.colorbar(
        positive_mappable,
        ax=axes[:, :2],
        shrink=0.72,
        label="Upward hydraulic force index $I_F^+$",
    )
    figure.colorbar(
        difference_mappable,
        ax=axes[:, 2],
        shrink=0.72,
        label="Lagged − erased $I_F^+$",
    )
    pressure = cases[0]["lagged_pressure"][int(round(COMMON_FIELD_TIME_S / 1.0e-4))][
        cases[0]["surface_particle_id"]
    ] - cases[0]["lagged_pressure"][0][cases[0]["surface_particle_id"]]
    figure.suptitle(
        f"Common negative-pressure recovery phase at t={COMMON_FIELD_TIME_S:.3f} s "
        f"(surface excess pressure {pressure:+.0f} Pa)",
        fontsize=14,
    )
    figure.text(
        0.5,
        -0.008,
        f"Field uses a one-lattice-layer local affine gradient plus display-only Gaussian sigma={display_sigma:g} pixels; all registered metrics use raw nearest-neighbour gradients.",
        ha="center",
        fontsize=9,
    )
    png = output / "phase_conditioned_uplift_2d.png"
    pdf = output / "phase_conditioned_uplift_2d.pdf"
    atomic_savefig(figure, png)
    atomic_savefig(figure, pdf)
    plt.close(figure)
    return png, pdf, {
        "time_s": COMMON_FIELD_TIME_S,
        "surface_excess_pressure_pa": pressure,
        "display_smoothing_sigma_pixels": display_sigma,
        "numeric_smoothing_applied": False,
        "spatial_gradient_reconstruction": {
            "method": "inverse-distance-weighted local affine pressure plane",
            "reference_lattice_radius": 1,
            "used_in_registered_metrics_or_gate": False,
        },
        "positive_colour_limit": positive_limit,
        "difference_colour_limit": difference_limit,
        "source_vtp": [
            file_record(field["path"]) for field in reconstructed_fields
        ],
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


def _summary_rows(cases: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for case in cases:
        first, second = case["branches"]
        phase = case["phase_qualification"]["most_delayed_eligible_probe"]
        raw_direction = case["raw_directional_branches"][-1]
        local_direction = case["secondary_local_linear_branches"][-1]
        rows.append(
            {
                "label": case["label"],
                "intrinsic_permeability_m2": case["validated_case"][
                    "intrinsic_permeability_m2"
                ],
                "liquid_saturation": case["validated_case"]["liquid_saturation"],
                "recovery_1_local_upward_activity_effect_percent": 100.0
                * first["local_upward_activity_impulse"][
                    "lagged_minus_phase_erased_fraction"
                ],
                "recovery_2_local_upward_activity_effect_percent": 100.0
                * second["local_upward_activity_impulse"][
                    "lagged_minus_phase_erased_fraction"
                ],
                "recovery_1_signed_support_effect_percent": 100.0
                * first["signed_support_force_impulse"][
                    "lagged_minus_phase_erased_fraction"
                ],
                "recovery_2_signed_support_effect_percent": 100.0
                * second["signed_support_force_impulse"][
                    "lagged_minus_phase_erased_fraction"
                ],
                "recovery_2_maximum_paired_local_effect_percent": 100.0
                * second["local_upward_activity_impulse"][
                    "maximum_paired_delta_fraction"
                ],
                "recovery_2_maximum_IF_effect_percent": 100.0
                * second["maximum_upward_IF_time_integral"][
                    "lagged_minus_phase_erased_fraction"
                ],
                "recovery_2_raw_horizontal_absolute_effect_percent": 100.0
                * raw_direction["horizontal_absolute_force_activity_impulse"][
                    "lagged_minus_phase_erased_fraction"
                ],
                "recovery_2_raw_vertical_absolute_effect_percent": 100.0
                * raw_direction["vertical_absolute_force_activity_impulse"][
                    "lagged_minus_phase_erased_fraction"
                ],
                "recovery_2_one_layer_horizontal_absolute_effect_percent": 100.0
                * local_direction["horizontal_absolute_force_activity_impulse"][
                    "lagged_minus_phase_erased_fraction"
                ],
                "recovery_2_one_layer_vertical_absolute_effect_percent": 100.0
                * local_direction["vertical_absolute_force_activity_impulse"][
                    "lagged_minus_phase_erased_fraction"
                ],
                "most_delayed_eligible_probe": phase["probe"] if phase else "",
                "most_delayed_eligible_phase_deg": (
                    phase["same_column_phase_difference_deg"] if phase else ""
                ),
                "short_time_mechanism_status": (
                    "PASS" if case["gate"]["passed"] else "FAIL"
                ),
            }
        )
    return rows


def report_artifact(
    cases: list[dict[str, Any]], generated_at: str, display_sigma: float
) -> dict[str, Any]:
    """Return the canonical bounded data-analytics technical-report payload."""

    summary_rows = _summary_rows(cases)
    recovery_rows: list[dict[str, Any]] = []
    history_rows: list[dict[str, Any]] = []
    for case, summary_row in zip(cases, summary_rows):
        permeability_label = f"{summary_row['intrinsic_permeability_m2']:.0e}"
        for index in (1, 2):
            recovery_rows.append(
                {
                    "label": case["label"],
                    "permeability_label": permeability_label,
                    "intrinsic_permeability_m2": summary_row[
                        "intrinsic_permeability_m2"
                    ],
                    "recovery": f"recovery {index}",
                    "local_effect_percent": summary_row[
                        f"recovery_{index}_local_upward_activity_effect_percent"
                    ],
                    "signed_effect_percent": summary_row[
                        f"recovery_{index}_signed_support_effect_percent"
                    ],
                    "liquid_saturation": summary_row["liquid_saturation"],
                    "mechanism_status": summary_row[
                        "short_time_mechanism_status"
                    ],
                }
            )
        for row in case["rows"]:
            if PRIMARY_CYCLE[0] <= row["time_s"] <= PRIMARY_CYCLE[1]:
                history_rows.append(
                    {
                        "label": case["label"],
                        "permeability_label": permeability_label,
                        "time_s": row["time_s"],
                        "surface_excess_pressure_pa": row[
                            "surface_excess_pressure_pa"
                        ],
                        "local_force_difference_n_per_m": row[
                            "delta_positive_force_n_per_m"
                        ],
                        "signed_force_difference_n_per_m": row[
                            "delta_signed_force_n_per_m"
                        ],
                        "is_primary_negative_recovery": bool(
                            case["branch_rows"][-1][0]["time_s"]
                            <= row["time_s"]
                            <= case["branch_rows"][-1][-1]["time_s"]
                        ),
                    }
                )
    passed_rows = [
        row
        for row in summary_rows
        if row["short_time_mechanism_status"] == "PASS"
    ]
    selected = max(
        passed_rows or summary_rows,
        key=lambda row: min(
            row["recovery_1_local_upward_activity_effect_percent"],
            row["recovery_2_local_upward_activity_effect_percent"],
        ),
    )
    selected_k = f"{selected['intrinsic_permeability_m2']:.0e}"
    supported_k = ", ".join(
        f"`{row['intrinsic_permeability_m2']:.0e} m2`" for row in passed_rows
    )
    if not supported_k:
        supported_k = "none of the audited permeability points"
    def sql_literal(value: Any) -> str:
        if value is None:
            return "NULL"
        if isinstance(value, bool):
            return "1" if value else "0"
        if isinstance(value, (int, float)):
            number = float(value)
            if not math.isfinite(number):
                raise ConditionedUpliftError("Report SQL cannot contain NaN/Inf")
            return repr(value)
        return "'" + str(value).replace("'", "''") + "'"

    def literal_query(name: str, rows: list[dict[str, Any]]) -> str:
        if not rows:
            raise ConditionedUpliftError(f"Report dataset {name} is empty")
        columns = list(rows[0])
        if any(list(row) != columns for row in rows):
            raise ConditionedUpliftError(f"Report dataset {name} columns differ")
        values = ",\n    ".join(
            "(" + ", ".join(sql_literal(row[column]) for column in columns) + ")"
            for row in rows
        )
        return (
            f"WITH {name} ({', '.join(columns)}) AS (\n  VALUES\n    {values}\n)\n"
            f"SELECT * FROM {name};"
        )

    headline_rows = [
        {
            "permeability_m2": selected["intrinsic_permeability_m2"],
            "local_effect_fraction": selected[
                "recovery_2_local_upward_activity_effect_percent"
            ]
            / 100.0,
            "repeat_local_effect_fraction": selected[
                "recovery_1_local_upward_activity_effect_percent"
            ]
            / 100.0,
            "signed_effect_fraction": selected[
                "recovery_2_signed_support_effect_percent"
            ]
            / 100.0,
            "paired_peak_fraction": selected[
                "recovery_2_maximum_paired_local_effect_percent"
            ]
            / 100.0,
        }
    ]
    query_by_dataset = {
        "headline": literal_query("headline", headline_rows),
        "case_summary": literal_query("case_summary", summary_rows),
        "recovery_effects": literal_query("recovery_effects", recovery_rows),
    }
    # The portable report validator requires an actual query for every native
    # card/chart/table. Execute the bounded literal queries here so the saved
    # snapshot is exactly query-produced; the physics calculation remains the
    # separately hash-bound Python audit rather than being misrepresented as SQL.
    for name, query in query_by_dataset.items():
        with sqlite3.connect(":memory:") as connection:
            connection.row_factory = sqlite3.Row
            observed = [dict(row) for row in connection.execute(query)]
        expected = {
            "headline": headline_rows,
            "case_summary": summary_rows,
            "recovery_effects": recovery_rows,
        }[name]
        if len(observed) != len(expected) or list(observed[0]) != list(expected[0]):
            raise ConditionedUpliftError(f"Report SQL did not reproduce {name}")

    common_query_metadata = {
        "engine": "SQLite",
        "language": "SQL",
        "description": (
            "Bounded report snapshot of hash-audited Python results. Physics "
            "derivatives and integrals remain documented in the companion audit."
        ),
        "filters": [
            "Sw=0.94",
            "both solver pressure-smoothing switches false",
            "fixed pipeline",
            "post-ramp cycles 1.3-2.6 s and 2.6-3.9 s",
            "surface pressure <= -10% of fitted amplitude after its trough",
        ],
        "metric_definitions": [
            "Local upward activity is the support-zone spatial integral of max(f_y,0).",
            "Signed support force is the spatial integral of signed upward hydraulic force density.",
            "Relative effect is lagged/phase-erased minus one on the same recovery branch.",
        ],
    }
    sources = [
        {
            "id": f"{name}_source",
            "label": f"Audited {name.replace('_', ' ')} snapshot",
            "query": {**common_query_metadata, "sql": query},
        }
        for name, query in query_by_dataset.items()
    ]
    title = "Short-time pipeline uplift risk during negative seabed pressure"
    return {
        "surface": "report",
        "manifest": {
            "version": 1,
            "surface": "report",
            "title": title,
            "description": (
                "Technical evidence report for the unsmoothed Sw=0.94 "
                "phase-conditioned fixed-pipeline comparison."
            ),
            "generatedAt": generated_at,
            "sources": sources,
            "cards": [
                {
                    "id": "local_effect_card",
                    "description": (
                        "Second negative-pressure recovery at the only point that "
                        "passes the short-time mechanism gate."
                    ),
                    "dataset": "headline",
                    "sourceId": "headline_source",
                    "metrics": [
                        {
                            "label": "Local upward-activity increase",
                            "field": "local_effect_fraction",
                            "format": "percent",
                            "signed": True,
                        },
                        {
                            "label": "Repeat recovery",
                            "field": "repeat_local_effect_fraction",
                            "format": "percent",
                            "signed": True,
                        },
                    ],
                },
                {
                    "id": "signed_effect_card",
                    "description": "Paired signed support-force impulse during recovery.",
                    "dataset": "headline",
                    "sourceId": "headline_source",
                    "metrics": [
                        {
                            "label": "Signed support-force increase",
                            "field": "signed_effect_fraction",
                            "format": "percent",
                            "signed": True,
                        },
                        {
                            "label": "Peak paired local advantage",
                            "field": "paired_peak_fraction",
                            "format": "percent",
                            "signed": True,
                        },
                    ],
                },
            ],
            "charts": [
                {
                    "id": "recovery_effect_chart",
                    "title": "Recovery-conditioned upward-force effects",
                    "subtitle": (
                        "Two repeated post-ramp recoveries; percentage relative to "
                        "the phase-erased pressure field"
                    ),
                    "intent": "comparison",
                    "question": (
                        "At which permeability does retained phase structure "
                        "increase short-time local upward forcing?"
                    ),
                    "rationale": (
                        "Grouped bars expose repeatability and the non-monotonic "
                        "permeability response without selecting one favourable frame."
                    ),
                    "type": "bar",
                    "dataset": "recovery_effects",
                    "sourceId": "recovery_effects_source",
                    "encodings": {
                        "x": {
                            "field": "permeability_label",
                            "type": "ordinal",
                            "label": "Intrinsic permeability k_s (m2)",
                        },
                        "y": {
                            "field": "local_effect_percent",
                            "type": "quantitative",
                            "label": "Lagged - erased (%)",
                            "unit": "%",
                        },
                        "color": {
                            "field": "recovery",
                            "type": "nominal",
                            "label": "Recovery branch",
                        },
                        "tooltip": [
                            {"field": "signed_effect_percent", "label": "Signed effect", "unit": "%"},
                            {"field": "mechanism_status", "label": "Mechanism gate"},
                        ],
                    },
                    "layout": "full",
                }
            ],
            "tables": [
                {
                    "id": "case_table",
                    "title": f"Audited {len(cases)}-point recovery comparison",
                    "subtitle": "Unsmoothed pressure data, Sw=0.94, fixed pipeline",
                    "dataset": "case_summary",
                    "defaultSort": {
                        "field": "intrinsic_permeability_m2",
                        "direction": "asc",
                    },
                    "density": "spacious",
                    "sourceId": "case_summary_source",
                    "layout": "full",
                    "columns": [
                        {"field": "intrinsic_permeability_m2", "label": "k_s (m2)", "format": "number"},
                        {"field": "recovery_1_local_upward_activity_effect_percent", "label": "Recovery 1 local (%)", "format": "number", "movement": True},
                        {"field": "recovery_2_local_upward_activity_effect_percent", "label": "Recovery 2 local (%)", "format": "number", "movement": True},
                        {"field": "recovery_2_signed_support_effect_percent", "label": "Recovery 2 signed (%)", "format": "number", "movement": True},
                        {"field": "most_delayed_eligible_probe", "label": "Resolved probe", "type": "text"},
                        {"field": "most_delayed_eligible_phase_deg", "label": "Probe phase (deg)", "format": "number"},
                        {"field": "short_time_mechanism_status", "label": "Short-time gate", "type": "text"},
                    ],
                }
            ],
            "blocks": [
                {"id": "title", "type": "markdown", "body": f"# {title}"},
                {
                    "id": "technical_summary",
                    "type": "markdown",
                    "sourceId": "case_summary_source",
                    "body": (
                        f"## Strongest repeatable short-time response: {selected_k} m2\n\n"
                        "At this audited point, retained phase structure changes "
                        "local upward-force activity by "
                        f"**{selected['recovery_1_local_upward_activity_effect_percent']:+.2f}%** "
                        "and "
                        f"**{selected['recovery_2_local_upward_activity_effect_percent']:+.2f}%** "
                        "in the two negative-pressure recoveries. The corresponding "
                        "signed support-force changes are "
                        f"**{selected['recovery_1_signed_support_effect_percent']:+.2f}%** "
                        "and "
                        f"**{selected['recovery_2_signed_support_effect_percent']:+.2f}%**."
                    ),
                },
                {"id": "metrics", "type": "metric-strip", "cardIds": ["local_effect_card", "signed_effect_card"]},
                {
                    "id": "finding",
                    "type": "markdown",
                    "sourceId": "case_summary_source",
                    "body": (
                        "## Recovery-phase forcing has a permeability window\n\n"
                        "The chart compares both repeated trough-to-recovery branches. "
                        f"The complete short-time mechanism gate is supported at {supported_k}. "
                        "This is a transient redistribution result: a complete-cycle "
                        "net-uplift failure remains a guardrail because positive and "
                        "negative parts can cancel over a full period."
                    ),
                },
                {"id": "chart", "type": "chart", "chartId": "recovery_effect_chart", "layout": "full"},
                {
                    "id": "definitions",
                    "type": "markdown",
                    "body": (
                        "## Scope and metric definitions\n\n"
                        "The comparison uses two post-ramp cycles (1.3-2.6 s and 2.6-3.9 s). Each recovery begins at the saved same-column surface-pressure trough and continues while the raw surface excess pressure remains at least 10% of its fitted amplitude below equilibrium. `Local upward activity` integrates `max(f_y,0)` over the registered pipeline-support zone; `signed support force` integrates signed `f_y` over the same cohort. They are paired fixed-state counterfactuals on identical SANISAND particle states."
                    ),
                },
                {"id": "table", "type": "table", "tableId": "case_table", "layout": "full"},
                {
                    "id": "method",
                    "type": "markdown",
                    "body": (
                        "## Raw derivatives and repeatability make the comparison auditable\n\n"
                        "Liquid-pressure gradients are reconstructed directly from ID-aligned raw pressure-database values at the saved material-point coordinates. Both solver smoothing switches are false. Gaussian sigma="
                        f"{display_sigma:g} pixels is used only in the companion two-dimensional raster; it never changes a derivative, integral, frame, threshold, or gate. The two recovery branches contain the same five physical-phase samples and provide an internal repeatability check."
                    ),
                },
                {
                    "id": "limits",
                    "type": "markdown",
                    "body": (
                        "## The result indicates potential risk, not realised liquefaction\n\n"
                        "This pressure-only fixed-state comparison does not evolve "
                        "separate lagged and phase-erased skeleton stresses, release "
                        "the pipe, or establish additional liquefaction. The supported "
                        "claim is a short-duration hydraulic-demand mechanism. Independent "
                        "paired replay and released-pipeline calculations are required "
                        "for a causal liquefaction or displacement claim."
                    ),
                },
                {
                    "id": "next_steps",
                    "type": "markdown",
                    "body": (
                        "## Recommended next step\n\n"
                        f"Use `{selected_k} m2`, `Sw=0.94`, and unsmoothed raw pressure "
                        "for the first independent lagged/phase-erased SANISAND replay. "
                        "Preserve the negative-pressure recovery metrics as registered "
                        "secondary outcomes and keep the complete-cycle net result as a "
                        "guardrail. Stop the stronger manuscript claim if separate "
                        "skeleton states do not reproduce the short-time advantage."
                    ),
                },
                {
                    "id": "questions",
                    "type": "markdown",
                    "body": (
                        "## Further questions\n\n"
                        "Does the recovery-phase local forcing persist after two-way stress coupling? Does it enlarge a spatially resolved `IF>=1` and low-effective-stress overlap? Does a released pipeline respond measurably within the same sub-cycle window?"
                    ),
                },
            ],
        },
        "snapshot": {
            "version": 1,
            "generatedAt": generated_at,
            "status": "ready",
            "datasets": {
                "headline": headline_rows,
                "case_summary": summary_rows,
                "recovery_effects": recovery_rows,
                "history_effects": history_rows,
            },
        },
        "sources": sources,
    }


def run(root: Path, labels: Iterable[str], output: Path, display_sigma: float) -> Path:
    root = root.resolve()
    labels = tuple(labels)
    if set(labels) == set(DEFAULT_LABELS) and len(labels) == len(DEFAULT_LABELS):
        registered_order = DEFAULT_LABELS
    elif set(labels) == set(EXTENDED_LABELS) and len(labels) == len(EXTENDED_LABELS):
        registered_order = EXTENDED_LABELS
    else:
        raise ConditionedUpliftError(
            "The conditioned screen requires either the exact registered four "
            "labels or the exact six-label upper-bound sensitivity extension"
        )
    cases_by_label = {label: load_case(root, label) for label in labels}
    cases = [cases_by_label[label] for label in registered_order]
    summary_png, summary_pdf = plot_summary(cases, output)
    field_png, field_pdf, field_audit = plot_2d(cases, output, display_sigma)
    csv_path = output / "phase_conditioned_uplift_summary.csv"
    atomic_write_csv(csv_path, _summary_rows(cases))

    case_audits: list[dict[str, Any]] = []
    for case in cases:
        case_audits.append(
            {
                "label": case["label"],
                "intrinsic_permeability_m2": case["validated_case"][
                    "intrinsic_permeability_m2"
                ],
                "liquid_saturation": case["validated_case"]["liquid_saturation"],
                "pressure_smoothing": False,
                "pipeline_fixed": True,
                "surface_particle_id": case["surface_particle_id"],
                "surface_fundamental_amplitude_pa": case["surface_amplitude_pa"],
                "surface_driver_maximum_pair_difference_pa": case[
                    "surface_driver_maximum_pair_difference_pa"
                ],
                "phase_qualification": case["phase_qualification"],
                "negative_pressure_recoveries": case["branches"],
                "raw_directional_recoveries": case["raw_directional_branches"],
                "secondary_one_layer_local_linear_recoveries": case[
                    "secondary_local_linear_branches"
                ],
                "short_time_mechanism_gate": case["gate"],
                "sources": case["validated_case"]["sources"],
            }
        )
    audit_path = output / "phase_conditioned_uplift.audit.json"
    audit = {
        "schema": SCHEMA,
        "classification": "fixed-state-negative-pressure-short-time-diagnostic",
        "technical_question": (
            "During the post-ramp recovery from a negative same-column seabed "
            "pressure trough, does retained phase structure increase short-time "
            "upward hydraulic forcing near the fixed pipeline?"
        ),
        "analyzer": file_record(Path(__file__)),
        "implementation_dependencies": {
            "permeability_summary": file_record(
                Path(permeability_summary.__file__).resolve()
            ),
            "gradient_force_module": file_record(
                Path(gradient_force_module.__file__).resolve()
            ),
        },
        "metric_contract": {
            "post_ramp_cycles_s": [list(value) for value in RECOVERY_CYCLES],
            "primary_cycle_s": list(PRIMARY_CYCLE),
            "surface_pressure_definition": (
                "raw liquid pressure at registered same-column top-surface particle "
                "minus its step-zero value"
            ),
            "recovery_definition": (
                "contiguous saved frames from the surface-pressure trough to the "
                "cycle end while p_excess <= -0.10 times fitted surface amplitude"
            ),
            "negative_core_fraction": SURFACE_NEGATIVE_CORE_FRACTION,
            "local_upward_activity": (
                "integral over the registered support cohort of max(f_y,0); "
                "not the spatial net resultant"
            ),
            "signed_support_force": (
                "spatial integral over the same cohort of signed upward hydraulic "
                "force density"
            ),
            "force_orientation": (
                "separate spatial integrals of absolute horizontal and absolute "
                "vertical hydraulic force density over the same support cohort; "
                "reported for the raw derivative and the one-layer sensitivity"
            ),
            "numeric_pressure_smoothing": False,
            "secondary_spatial_reconstruction": {
                "method": "inverse-distance-weighted local affine pressure plane",
                "reference_lattice_radius": 1,
                "purpose": "lattice-layering sensitivity and two-dimensional display",
                "used_in_registered_metrics_or_gate": False,
                "two_layer_reconstruction_permitted_for_inference": False,
            },
            "minimum_material_effect_fraction": MINIMUM_MATERIAL_EFFECT_FRACTION,
        },
        "claim_boundary": (
            "This paired fixed-state counterfactual can diagnose transient hydraulic "
            "risk.  It cannot establish realised liquefaction, additional skeleton "
            "stress loss, or released-pipeline motion."
        ),
        "cases": case_audits,
        "field_snapshot": field_audit,
        "artifacts": {
            "summary_csv": file_record(csv_path),
            "summary_png": file_record(summary_png),
            "summary_pdf": file_record(summary_pdf),
            "field_png": file_record(field_png),
            "field_pdf": file_record(field_pdf),
        },
        "audit_path": str(audit_path.resolve()),
    }
    atomic_write_json(audit_path, audit)
    from datetime import datetime, timezone

    generated_at = datetime.now(timezone.utc).isoformat()
    atomic_write_json(
        output / "phase_conditioned_uplift_report_artifact.json",
        report_artifact(cases, generated_at, display_sigma),
    )
    return audit_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--root", type=Path, default=Path(__file__).resolve().parent
    )
    parser.add_argument("--labels", nargs="+", default=list(DEFAULT_LABELS))
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("analysis/phase_conditioned_uplift"),
    )
    parser.add_argument(
        "--display-smoothing-sigma", type=float, default=DISPLAY_SIGMA_DEFAULT
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    root = args.root.resolve()
    output = args.output_dir
    if not output.is_absolute():
        output = root / output
    audit = run(root, args.labels, output.resolve(), args.display_smoothing_sigma)
    print(audit)


if __name__ == "__main__":
    main()
