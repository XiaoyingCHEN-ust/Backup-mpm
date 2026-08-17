#!/usr/bin/env python3
"""Plot the audited phase-lag-to-pipeline gradient-force bridge.

All reported metrics are computed from the raw particle fields.  Optional
Gaussian smoothing is applied only to a raster used to display the three field
panels; it never feeds the force, impulse, threshold-area, or pipeline-response
calculations.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.tri as mtri  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.patches import Circle  # noqa: E402
from scipy.ndimage import gaussian_filter  # noqa: E402

from analyze_study import (  # noqa: E402
    expected_particle_steps,
    nominal_reference_particle_area,
    particle_step,
    pipeline_support_roi,
    read_initial_coordinates,
    read_pipeline_history,
    read_vtp,
    resolve_particle_input_path,
    validate_completed_result,
    validate_pipeline_history_grid,
)
from plot_manuscript_figures import _masked_triangulation  # noqa: E402


SCHEMA = "pipeline-gradient-force-bridge-v1"
CASES = {
    "RL": "configs/screen/04_RL.json",
    "RE": "configs/screen/04_RE.json",
}
COLORS = {"RL": "#B13A3A", "RE": "#2468A2"}


@dataclass(frozen=True)
class RawForceMetrics:
    signed_force_n_per_m: float
    positive_force_n_per_m: float
    positive_force_index_area_m2: float
    maximum_upward_if: float
    critical_area_m2: float


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def upward_excess_force_density(
    liquid_seepage_forces: np.ndarray,
    liquid_densities: np.ndarray,
    gravity: np.ndarray,
) -> np.ndarray:
    """Return signed upward wave-induced pressure-gradient force (N/m3)."""

    forces = np.asarray(liquid_seepage_forces, dtype=np.float64)
    densities = np.asarray(liquid_densities, dtype=np.float64).reshape(-1)
    gravity = np.asarray(gravity, dtype=np.float64).reshape(-1)
    if forces.ndim != 2 or forces.shape[0] != densities.size:
        raise ValueError("Seepage-force and density arrays are incompatible")
    if gravity.size < 2 or forces.shape[1] < gravity.size:
        raise ValueError("Gravity and seepage-force dimensions are incompatible")
    if not np.all(np.isfinite(forces[:, : gravity.size])):
        raise ValueError("Seepage force contains NaN/Inf")
    if not np.all(np.isfinite(densities)) or np.any(densities <= 0.0):
        raise ValueError("Liquid density must be finite and positive")
    magnitude = float(np.linalg.norm(gravity))
    if not math.isfinite(magnitude) or magnitude <= 0.0:
        raise ValueError("Gravity magnitude must be finite and positive")
    upward = -gravity / magnitude
    return (forces[:, : gravity.size] + densities[:, None] * gravity) @ upward


def raw_force_metrics(
    excess_upward_force_density: np.ndarray,
    gamma_sub: np.ndarray,
    represented_volumes: np.ndarray,
    support_mask: np.ndarray,
) -> tuple[RawForceMetrics, np.ndarray]:
    """Integrate raw force density over the initial-ID support cohort."""

    force = np.asarray(excess_upward_force_density, dtype=np.float64).reshape(-1)
    weights = np.asarray(gamma_sub, dtype=np.float64).reshape(-1)
    volumes = np.asarray(represented_volumes, dtype=np.float64).reshape(-1)
    mask = np.asarray(support_mask, dtype=bool).reshape(-1)
    if not (force.size == weights.size == volumes.size == mask.size):
        raise ValueError("Force, weight, volume, and mask counts must match")
    if not np.any(mask):
        raise ValueError("Pipeline support cohort is empty")
    if not np.all(np.isfinite(force)):
        raise ValueError("Excess force density contains NaN/Inf")
    if not np.all(np.isfinite(weights)) or np.any(weights <= 0.0):
        raise ValueError("Submerged unit weight must be finite and positive")
    if not np.all(np.isfinite(volumes)) or np.any(volumes <= 0.0):
        raise ValueError("Particle volumes must be finite and positive")
    signed_if = force / weights
    upward_if = np.maximum(signed_if, 0.0)
    selected = mask
    metrics = RawForceMetrics(
        signed_force_n_per_m=float(np.sum(force[selected] * volumes[selected])),
        positive_force_n_per_m=float(
            np.sum(np.maximum(force[selected], 0.0) * volumes[selected])
        ),
        positive_force_index_area_m2=float(
            np.sum(upward_if[selected] * volumes[selected])
        ),
        maximum_upward_if=float(np.max(upward_if[selected])),
        critical_area_m2=float(
            np.sum(volumes[selected & (upward_if >= 1.0)])
        ),
    )
    return metrics, upward_if


def display_smoothed_grid(
    points: np.ndarray,
    values: np.ndarray,
    config: dict[str, Any],
    pipeline_center: np.ndarray,
    bounds: tuple[float, float, float, float],
    sigma: float,
    *,
    nx: int = 420,
    ny: int = 280,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return a display-only interpolated raster without modifying raw values."""

    if not math.isfinite(sigma) or sigma < 0.0:
        raise ValueError("Display smoothing sigma must be finite and non-negative")
    raw = np.asarray(values, dtype=np.float64)
    if not np.all(np.isfinite(raw)):
        raise ValueError("Display field contains NaN/Inf")
    triangulation = _masked_triangulation(points, config, pipeline_center)
    x0, x1, y0, y1 = bounds
    x = np.linspace(x0, x1, nx)
    y = np.linspace(y0, y1, ny)
    xx, yy = np.meshgrid(x, y)
    interpolated = mtri.LinearTriInterpolator(triangulation, raw)(xx, yy)
    data = np.asarray(interpolated.filled(np.nan), dtype=np.float64)
    valid = np.isfinite(data)
    if sigma > 0.0:
        numerator = gaussian_filter(np.where(valid, data, 0.0), sigma=sigma)
        denominator = gaussian_filter(valid.astype(np.float64), sigma=sigma)
        smoothed = np.divide(
            numerator,
            denominator,
            out=np.full_like(numerator, np.nan),
            where=denominator > 1.0e-8,
        )
        data = np.where(valid, smoothed, np.nan)
    return xx, yy, data


def ordered_frame(path: Path, particle_count: int) -> tuple[np.ndarray, dict[str, np.ndarray]]:
    points, arrays = read_vtp(path)
    required = (
        "ids",
        "liquid_seepage_forces",
        "liquid_densities",
        "gamma_sub",
        "volumes",
    )
    missing = [name for name in required if name not in arrays]
    if missing:
        raise KeyError(f"{path} is missing {missing}")
    raw_ids = np.asarray(arrays["ids"], dtype=np.float64).reshape(-1)
    ids = np.rint(raw_ids).astype(np.int64)
    order = np.argsort(ids)
    if raw_ids.size != particle_count or not np.array_equal(
        raw_ids, ids.astype(np.float64)
    ) or not np.array_equal(ids[order], np.arange(particle_count)):
        raise ValueError(f"{path} particle IDs are not exactly 0..N-1")
    ordered_points = np.asarray(points, dtype=np.float64)[order, :2]
    ordered_arrays = {
        name: np.asarray(values)[order] for name, values in arrays.items()
    }
    return ordered_points, ordered_arrays


def read_case(root: Path, code: str, config_relative: str) -> dict[str, Any]:
    config_path = (root / config_relative).resolve()
    config = json.loads(config_path.read_text(encoding="utf-8"))
    result_base = Path(config["post_processing"]["path"])
    if not result_base.is_absolute():
        result_base = root / result_base
    result_directory = (result_base / config["analysis"]["uuid"]).resolve()
    files = sorted(result_directory.glob("particle*.vtp"), key=particle_step)
    completion = validate_completed_result(
        config_path, root, config, result_directory, files
    )
    expected_steps = expected_particle_steps(config)
    if [particle_step(path) for path in files] != expected_steps:
        raise ValueError(f"{code} saved-frame grid is incomplete")
    coordinates = read_initial_coordinates(resolve_particle_input_path(config, root))
    reference_particle_area = nominal_reference_particle_area(coordinates)
    support = pipeline_support_roi(config, coordinates)
    if support is None:
        raise ValueError(f"{code} has no registered pipeline support cohort")
    gravity = np.asarray(config["external_loading_conditions"]["gravity"], dtype=float)
    rows: list[dict[str, float]] = []
    selected_fields: dict[float, dict[str, Any]] = {}
    dt = float(config["analysis"]["dt"])
    for path in files:
        points, arrays = ordered_frame(path, len(coordinates))
        excess = upward_excess_force_density(
            arrays["liquid_seepage_forces"], arrays["liquid_densities"], gravity
        )
        metrics, upward_if = raw_force_metrics(
            excess, arrays["gamma_sub"], arrays["volumes"], support
        )
        time_s = particle_step(path) * dt
        rows.append(
            {
                "time_s": time_s,
                "signed_force_n_per_m": metrics.signed_force_n_per_m,
                "positive_force_n_per_m": metrics.positive_force_n_per_m,
                "positive_force_index_area_m2": metrics.positive_force_index_area_m2,
                "maximum_upward_if": metrics.maximum_upward_if,
                "critical_area_m2": metrics.critical_area_m2,
            }
        )
        selected_fields[time_s] = {
            "path": path,
            "points": points,
            "upward_if": upward_if,
            "ids": np.arange(len(coordinates), dtype=np.int64),
        }
    history = read_pipeline_history(result_directory)
    history_audit = validate_pipeline_history_grid(history, config)
    return {
        "code": code,
        "config": config,
        "config_path": config_path,
        "result_directory": result_directory,
        "completion": completion,
        "coordinates": coordinates,
        "reference_particle_area": reference_particle_area,
        "support": support,
        "rows": rows,
        "fields": selected_fields,
        "history": history,
        "history_audit": history_audit,
    }


def nearest_history_row(history: list[dict[str, float]], target_time: float) -> dict[str, float]:
    return min(history, key=lambda row: abs(row["time"] - target_time))


def write_csv(path: Path, cases: dict[str, dict[str, Any]]) -> None:
    rl_rows = cases["RL"]["rows"]
    re_rows = cases["RE"]["rows"]
    if len(rl_rows) != len(re_rows):
        raise ValueError("RL and RE histories have different lengths")
    metric_names = tuple(name for name in rl_rows[0] if name != "time_s")
    fieldnames = ["time_s"]
    for name in metric_names:
        fieldnames.extend((f"RL_{name}", f"RE_{name}", f"delta_{name}"))
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        for left, right in zip(rl_rows, re_rows):
            if not math.isclose(left["time_s"], right["time_s"], abs_tol=1.0e-12):
                raise ValueError("RL and RE histories are not time aligned")
            row: dict[str, float] = {"time_s": left["time_s"]}
            for name in metric_names:
                row[f"RL_{name}"] = left[name]
                row[f"RE_{name}"] = right[name]
                row[f"delta_{name}"] = left[name] - right[name]
            writer.writerow(row)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir", type=Path, default=Path("analysis/screen/figures")
    )
    parser.add_argument(
        "--display-smoothing-sigma",
        type=float,
        default=0.75,
        help="Gaussian sigma in raster pixels; display only, never used in metrics",
    )
    args = parser.parse_args()

    root = Path(__file__).resolve().parent
    output = args.output_dir if args.output_dir.is_absolute() else root / args.output_dir
    output.mkdir(parents=True, exist_ok=True)
    cases = {
        code: read_case(root, code, relative) for code, relative in CASES.items()
    }
    if not np.array_equal(cases["RL"]["support"], cases["RE"]["support"]):
        raise ValueError("RL and RE pipeline-support cohorts differ")
    if not math.isclose(
        cases["RL"]["reference_particle_area"],
        cases["RE"]["reference_particle_area"],
        rel_tol=0.0,
        abs_tol=1.0e-15,
    ):
        raise ValueError("RL and RE reference particle areas differ")
    times = np.asarray([row["time_s"] for row in cases["RL"]["rows"]])
    if not np.array_equal(
        times, np.asarray([row["time_s"] for row in cases["RE"]["rows"]])
    ):
        raise ValueError("RL and RE particle times differ")

    registered_path = output / "figure8_common_field_frames.json"
    registered = json.loads(registered_path.read_text(encoding="utf-8"))
    target_time = float(registered["selected_common_time_s"])
    target_index = int(np.argmin(np.abs(times - target_time)))
    if not math.isclose(times[target_index], target_time, abs_tol=1.0e-12):
        raise ValueError("Registered t_IF is absent from the raw histories")

    target = {
        code: case["fields"][target_time] for code, case in cases.items()
    }
    rl_if = target["RL"]["upward_if"]
    re_if = target["RE"]["upward_if"]
    difference_if = rl_if - re_if
    config = cases["RL"]["config"]
    pipeline = config["analysis"]["rigid_pipeline"]
    initial_center = np.asarray(pipeline["center"], dtype=float)
    radius = float(pipeline["outer_radius"])
    diameter = 2.0 * radius
    liquid = config["materials"][1]
    bed_y = float(liquid["sea_level"]) - float(liquid["depth_left"])
    bounds = (
        initial_center[0] - 1.5 * diameter,
        initial_center[0] + 1.5 * diameter,
        bed_y - 2.0 * diameter,
        bed_y,
    )
    poses: dict[str, np.ndarray] = {}
    for code, case in cases.items():
        history_row = nearest_history_row(case["history"], target_time)
        poses[code] = initial_center + np.asarray(
            [history_row["displacement_x"], history_row["displacement_y"]]
        )

    figure = plt.figure(figsize=(13.4, 8.2), constrained_layout=True)
    grid = figure.add_gridspec(2, 3, height_ratios=(1.0, 0.92))
    top_axes = [figure.add_subplot(grid[0, column]) for column in range(3)]
    bottom_axes = [figure.add_subplot(grid[1, column]) for column in range(3)]
    common_max = float(
        max(
            np.max(rl_if[cases["RL"]["support"]]),
            np.max(re_if[cases["RE"]["support"]]),
            1.0,
        )
    )
    difference_max = float(
        max(np.max(np.abs(difference_if[cases["RL"]["support"]])), 1.0e-12)
    )
    field_specs = (
        ("RL", rl_if, "(a) Retained phase lag (RL)", "magma", 0.0, common_max),
        ("RE", re_if, "(b) Phase-erased (RE)", "magma", 0.0, common_max),
        (
            "RL",
            difference_if,
            "(c) Raw-ID difference, RL - RE",
            "RdBu_r",
            -difference_max,
            difference_max,
        ),
    )
    artists = []
    for axis, (code, values, title, cmap, lower, upper) in zip(top_axes, field_specs):
        _, _, raster = display_smoothed_grid(
            target[code]["points"],
            values,
            config,
            poses[code],
            bounds,
            args.display_smoothing_sigma,
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
        axis.add_patch(
            Circle(poses[code], radius, facecolor="white", edgecolor="0.1", lw=1.2)
        )
        if code in ("RL", "RE") and title[1:2] in ("a", "b"):
            critical = cases[code]["support"] & (values >= 1.0)
            axis.scatter(
                target[code]["points"][critical, 0],
                target[code]["points"][critical, 1],
                s=12,
                facecolors="none",
                edgecolors="white",
                linewidths=0.7,
                label="raw IF >= 1",
            )
            axis.legend(loc="lower left", fontsize=7, frameon=True)
        axis.set_title(title, fontsize=10)
        axis.set_xlabel("x (m)")
        axis.set_ylabel("y (m)")
    figure.colorbar(
        artists[0], ax=top_axes[:2], label="Upward excess seepage-force index, IF"
    )
    figure.colorbar(artists[2], ax=top_axes[2], label="Delta IF (RL - RE)")

    raw = {
        code: {
            name: np.asarray([row[name] for row in case["rows"]], dtype=float)
            for name in case["rows"][0]
            if name != "time_s"
        }
        for code, case in cases.items()
    }
    release_time = float(pipeline["release_time"])
    for code in ("RL", "RE"):
        bottom_axes[0].plot(
            times,
            raw[code]["positive_force_n_per_m"],
            color=COLORS[code],
            linewidth=1.35,
            label=code,
        )
    bottom_axes[0].axvline(release_time, color="0.35", ls="--", lw=0.9)
    bottom_axes[0].axvline(target_time, color="0.35", ls=":", lw=0.9)
    bottom_axes[0].set_title("(d) Raw support-zone positive upward force", fontsize=10)
    bottom_axes[0].set_xlabel("Time (s)")
    bottom_axes[0].set_ylabel("$F^+_{up}$ (N m$^{-1}$)")
    bottom_axes[0].legend()

    delta_force = (
        raw["RL"]["positive_force_n_per_m"]
        - raw["RE"]["positive_force_n_per_m"]
    )
    bottom_axes[1].axhline(0.0, color="0.2", lw=0.8)
    bottom_axes[1].plot(times, delta_force, color="0.25", lw=1.0)
    bottom_axes[1].fill_between(
        times, 0.0, delta_force, where=delta_force >= 0.0, color="#D58B30", alpha=0.75,
        label="RL > RE",
    )
    bottom_axes[1].fill_between(
        times, 0.0, delta_force, where=delta_force < 0.0, color="#5B8DB8", alpha=0.75,
        label="RE > RL",
    )
    bottom_axes[1].axvline(release_time, color="0.35", ls="--", lw=0.9)
    bottom_axes[1].axvline(target_time, color="0.35", ls=":", lw=0.9)
    bottom_axes[1].set_title("(e) Phase-specific force redistribution", fontsize=10)
    bottom_axes[1].set_xlabel("Time (s)")
    bottom_axes[1].set_ylabel("Delta $F^+_{up}$, RL - RE (N m$^{-1}$)")
    bottom_axes[1].legend(fontsize=8)

    for code, case in cases.items():
        history_times = np.asarray([row["time"] for row in case["history"]])
        uplift = np.asarray([row["displacement_y"] for row in case["history"]]) / diameter
        bottom_axes[2].plot(
            history_times,
            uplift,
            color=COLORS[code],
            linewidth=1.35,
            label=code,
        )
    bottom_axes[2].axvline(release_time, color="0.35", ls="--", lw=0.9)
    bottom_axes[2].axvline(target_time, color="0.35", ls=":", lw=0.9)
    bottom_axes[2].set_title("(f) Pipeline engineering response", fontsize=10)
    bottom_axes[2].set_xlabel("Time (s)")
    bottom_axes[2].set_ylabel("Pipeline uplift, $u_y/D$")
    bottom_axes[2].legend()
    for axis in bottom_axes:
        axis.grid(True, color="0.88", linewidth=0.6)

    figure.suptitle(
        "Phase-lag control of pipeline-zone upward pressure-gradient force\n"
        f"registered $t_{{IF}}={target_time:.3f}$ s; all metrics use raw particles; "
        f"field raster only: display smoothing sigma={args.display_smoothing_sigma:g} pixels",
        fontsize=13,
    )
    base = output / "figure8_pipeline_gradient_force_bridge"
    figure.savefig(base.with_suffix(".png"), dpi=300, bbox_inches="tight")
    figure.savefig(base.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(figure)

    csv_path = output / "figure8_pipeline_gradient_force_bridge.csv"
    write_csv(csv_path, cases)
    target_metrics = {
        code: case["rows"][target_index] for code, case in cases.items()
    }
    positive_impulse = {
        code: float(
            np.trapz(raw[code]["positive_force_n_per_m"], times)
        )
        for code in ("RL", "RE")
    }
    peak_positive_force = {
        code: float(np.max(raw[code]["positive_force_n_per_m"]))
        for code in ("RL", "RE")
    }
    pipeline_max = {}
    for code, case in cases.items():
        uplift = np.abs(
            np.asarray([row["displacement_y"] for row in case["history"]])
        ) / diameter
        pipeline_max[code] = float(np.max(uplift))
    audit = {
        "schema": SCHEMA,
        "status": "audited-registered-screen-derived-figure",
        "definitions": {
            "upward_excess_force_density": (
                "dot(-grad(p_l) + rho_l*g, -g/|g|) in N/m3"
            ),
            "positive_support_force": (
                "integral over the registered initial-ID pipeline support cohort "
                "of max(upward excess force density, 0), using raw current "
                "particle volumes; N/m for unit out-of-plane thickness"
            ),
            "positive_force_impulse": "trapezoidal time integral of raw positive support force",
            "display_smoothing": (
                "linear triangulation plus Gaussian raster smoothing for display only; "
                "not used by any metric, threshold, colour limit, or frame selection"
            ),
        },
        "display_smoothing_sigma_pixels": args.display_smoothing_sigma,
        "registered_frame": {
            "time_s": target_time,
            "selection_manifest": str(registered_path.resolve()),
            "selection_manifest_sha256": sha256(registered_path),
            "selection_rule": registered["selection_rule"],
        },
        "support_cohort": {
            "particle_count": int(np.count_nonzero(cases["RL"]["support"])),
            "initial_reference_area_m2": float(
                np.count_nonzero(cases["RL"]["support"])
                * cases["RL"]["reference_particle_area"]
            ),
            "reference_particle_area_m2": float(
                cases["RL"]["reference_particle_area"]
            ),
        },
        "cases": {
            code: {
                "config": str(case["config_path"]),
                "config_sha256": sha256(case["config_path"]),
                "completion_sentinel": str(
                    (case["result_directory"] / "pipeline_completion.json").resolve()
                ),
                "completion_sentinel_sha256": sha256(
                    case["result_directory"] / "pipeline_completion.json"
                ),
                "selected_vtp": str(target[code]["path"].resolve()),
                "selected_vtp_sha256": sha256(target[code]["path"]),
                "target_raw_metrics": target_metrics[code],
                "peak_positive_force_n_per_m": peak_positive_force[code],
                "positive_force_impulse_n_s_per_m": positive_impulse[code],
                "maximum_abs_pipeline_uplift_over_D": pipeline_max[code],
            }
            for code, case in cases.items()
        },
        "lagged_minus_phase_erased": {
            "target_positive_force_fraction": (
                target_metrics["RL"]["positive_force_n_per_m"]
                / target_metrics["RE"]["positive_force_n_per_m"]
                - 1.0
            ),
            "target_maximum_IF_fraction": (
                target_metrics["RL"]["maximum_upward_if"]
                / target_metrics["RE"]["maximum_upward_if"]
                - 1.0
            ),
            "target_critical_area_ratio": (
                target_metrics["RL"]["critical_area_m2"]
                / target_metrics["RE"]["critical_area_m2"]
            ),
            "positive_force_impulse_fraction": (
                positive_impulse["RL"] / positive_impulse["RE"] - 1.0
            ),
            "peak_positive_force_fraction": (
                peak_positive_force["RL"] / peak_positive_force["RE"] - 1.0
            ),
            "maximum_abs_pipeline_uplift_over_D_fraction": (
                pipeline_max["RL"] / pipeline_max["RE"] - 1.0
            ),
        },
        "interpretation": (
            "The retained lag focuses the upward gradient more strongly at the "
            "registered t_IF phase, but does not increase the full-record positive "
            "force impulse, peak positive force, or engineering-scale pipe uplift. "
            "The defensible mechanism is phase-specific redistribution, not a "
            "universal increase in liquefaction severity."
        ),
        "csv": str(csv_path.resolve()),
        "csv_sha256": sha256(csv_path),
    }
    audit_path = output / "figure8_pipeline_gradient_force_bridge.audit.json"
    audit_path.write_text(json.dumps(audit, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(audit["lagged_minus_phase_erased"], indent=2))
    print(f"Wrote pipeline gradient-force bridge to {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
