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


SCHEMA = "pipeline-phase-lag-exploratory-replay-analysis-v1"
DEFAULT_LABEL = "k1e-13_sw094_nosmooth"
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


REPLAYS = (Replay("lagged", "04_RL.json"), Replay("phase_erased", "05_RE.json"))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def validate_phase_metadata(root: Path, label: str) -> dict[str, Any]:
    path = (
        root
        / "pressure_databases/phase_lag_exploratory"
        / label
        / "phase_erased/phase_erased_metadata.json"
    )
    metadata = json.loads(path.read_text(encoding="utf-8"))
    if metadata.get("schema") != "phase-erased-pressure-control-v1":
        raise ValueError("Phase-erased metadata schema is invalid")
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
    for name, expected in (
        ("period_s", 1.3),
        ("fit_start_s", 2.6),
        ("fit_end_s", 3.9),
        ("ramp_time_s", 1.3),
    ):
        if not math.isclose(float(transform.get(name, math.nan)), expected):
            raise ValueError(f"Phase transform {name} changed")
    reference = transform.get("surface_reference", {})
    if (
        reference.get("method") != "registered_top_surface_particle_set"
        or int(reference.get("reference_particle_count", 0)) != 150
    ):
        raise ValueError("Phase transform does not use the 150-point top row")
    for group in ("source", "output"):
        record = metadata.get(group, {})
        for kind in ("points", "values"):
            artifact = Path(str(record.get(kind, "")))
            if not artifact.is_absolute():
                artifact = root / artifact
            if not artifact.is_file() or sha256(artifact) != record.get(
                f"{kind}_sha256"
            ):
                raise ValueError(f"Phase metadata {group} {kind} hash mismatch")
    qa = metadata.get("qa", {})
    if (
        float(qa.get("max_abs_mean_change_pa", math.inf)) > 1.0e-10
        or float(qa.get("max_abs_fundamental_amplitude_change_pa", math.inf))
        > 1.0e-10
    ):
        raise ValueError("Phase transform changed fitted mean or amplitude")
    return {"path": str(path.resolve()), "sha256": sha256(path), "metadata": metadata}


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


def read_replay(root: Path, label: str, replay: Replay) -> dict[str, Any]:
    config_path = root / "configs/phase_lag_exploratory" / label / replay.config_name
    config = json.loads(config_path.read_text(encoding="utf-8"))
    analysis = config["analysis"]
    pressure = analysis.get("prescribed_phase_pressures", {})
    if bool(analysis.get("pressure_smoothing")) or bool(
        analysis.get("pressure_smoothing_in_loop")
    ):
        raise ValueError(f"{replay.role}: smoothing must remain disabled")
    if not bool(analysis["rigid_pipeline"]["fixed"]):
        raise ValueError(f"{replay.role}: pipeline must remain fixed")
    if bool(config["materials"][1]["wave_pressure"]):
        raise ValueError(f"{replay.role}: replay must disable physical wave forcing")
    if not bool(pressure.get("enable")) or bool(pressure.get("write")):
        raise ValueError(f"{replay.role}: invalid pressure read/write mode")
    pressure_directory = root / pressure["path"]
    pressure_prefix = str(pressure["file_prefix"])
    pressure_artifact_paths = {
        "points": pressure_directory / f"{pressure_prefix}_points.txt",
        "values": pressure_directory / f"{pressure_prefix}_values.bin",
    }
    metadata_path = pressure_directory / f"{pressure_prefix}_metadata.json"
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
        porosity = np.asarray(arrays["porosities"], dtype=float).reshape(-1)
        outside = np.any(
            (points < mesh_min - 1.0e-8) | (points > mesh_max + 1.0e-8), axis=1
        )
        times.append(particle_step(path) * float(analysis["dt"]))
        points_rows.append(points)
        pressure_rows.append(
            np.asarray(arrays["PIC_pore_pressure_excess"], dtype=float).reshape(-1)
        )
        seepage_rows.append(seepage)
        ratio_rows.append(stress_ratio)
        volume_rows.append(np.asarray(arrays["volumes"], dtype=float).reshape(-1))
        velocity_maxima.append(float(np.max(np.linalg.norm(velocities, axis=1))))
        displacement_maxima.append(
            float(np.max(np.linalg.norm(displacements, axis=1)))
        )
        porosity_minima.append(float(np.min(porosity)))
        porosity_maxima.append(float(np.max(porosity)))
        outside_counts.append(int(np.count_nonzero(outside)))
        if np.any(porosity <= 0.0) or np.any(porosity >= 1.0):
            raise ValueError(f"{path}: porosity lies outside (0,1)")

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


def plot_history(output: Path, lagged: dict[str, Any], erased: dict[str, Any]) -> None:
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
    figure.suptitle(
        r"Exploratory replay: $k=10^{-13}$ m$^2$, $S_w=0.94$, no smoothing; fixed pipeline"
    )
    for suffix in ("png", "pdf"):
        figure.savefig(output / f"lagged_vs_phase_erased_liquefaction_history.{suffix}", dpi=220)
    plt.close(figure)


def plot_fields(
    output: Path,
    lagged: dict[str, Any],
    erased: dict[str, Any],
    index: int,
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
        r"$k=10^{-13}$ m$^2$, $S_w=0.94$, no smoothing; fields shown at initial material-point coordinates"
    )
    for suffix in ("png", "pdf"):
        figure.savefig(output / f"lagged_vs_phase_erased_2d.{suffix}", dpi=220)
    plt.close(figure)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--label", default=DEFAULT_LABEL)
    args = parser.parse_args()
    root = Path(__file__).resolve().parent
    output = root / "analysis/phase_lag_exploratory" / args.label / "replay_comparison"
    output.mkdir(parents=True, exist_ok=True)
    phase_metadata_audit = validate_phase_metadata(root, args.label)
    lagged = read_replay(root, args.label, REPLAYS[0])
    erased = read_replay(root, args.label, REPLAYS[1])
    lagged_physics = copy.deepcopy(lagged["config"])
    erased_physics = copy.deepcopy(erased["config"])
    for config in (lagged_physics, erased_physics):
        config["title"] = "paired-replay"
        config["analysis"]["uuid"] = "paired-replay"
        pressure = config["analysis"]["prescribed_phase_pressures"]
        pressure["path"] = "paired-replay"
        pressure["file_prefix"] = "paired-replay"
    if lagged_physics != erased_physics:
        raise ValueError("Replay configs differ in fields other than pressure source identity")
    if not np.array_equal(lagged["times"], erased["times"]):
        raise ValueError("Lagged and phase-erased saved-time grids differ")
    pressure_difference_q95 = np.percentile(
        np.abs(lagged["pressure"] - erased["pressure"]), 95, axis=1
    )
    selected_index = int(np.argmax(pressure_difference_q95))
    write_history(output / "lagged_vs_phase_erased_history.csv", lagged, erased)
    plot_history(output, lagged, erased)
    plot_fields(output, lagged, erased, selected_index)

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
        "parameters": {
            "intrinsic_permeability_m2": 1.0e-13,
            "liquid_saturation": 0.94,
            "pressure_smoothing": False,
            "pipeline_fixed": True,
        },
        "frame_selection": {
            "rule": "maximum across frames of particlewise q95 absolute lagged-minus-erased excess-pressure difference",
            "time_s": float(lagged["times"][selected_index]),
            "q95_absolute_pressure_difference_pa": float(
                pressure_difference_q95[selected_index]
            ),
        },
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
