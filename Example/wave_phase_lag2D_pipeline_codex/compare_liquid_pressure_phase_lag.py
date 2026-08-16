#!/usr/bin/env python3
"""Build a common-scale 2-D comparison of exploratory phase-lag probes."""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np

from analyze_study import build_probes, expected_particle_steps, read_vtp
from plot_liquid_pressure_phase_lag import (
    FIT_END_S,
    FIT_START_S,
    MINIMUM_HARMONIC_R_SQUARED,
    MINIMUM_LOCAL_AMPLITUDE_RATIO,
    WAVE_PERIOD_S,
    align_liquid_pressure,
    draw_field,
    harmonic_fields,
    harmonic_r_squared,
    particle_step,
    read_initial_coordinates,
    same_column_surface_values,
    wrap_degrees,
)


SCHEMA = "pipeline-liquid-pressure-phase-lag-comparison-v1"
SNAPSHOT_TIME_S = 2.925


@dataclass(frozen=True)
class Case:
    key: str
    title: str
    config: str
    result: str
    checkpoint: str


CASES = (
    Case(
        "baseline",
        "Baseline\n$k=9.79\\times10^{-12}$, $S_w=0.94$, smooth",
        "configs/screen/02_HS.json",
        "results/screen/PLP_SCREEN_HS_SANI",
        "results/screen/PLP_SCREEN_HS_EQ/particle40000.vtp",
    ),
    Case(
        "k3_sw094",
        "$k=3\\times10^{-12}$, $S_w=0.94$\nsmooth",
        "configs/phase_lag_exploratory/k3e-12_sw094/02_WAVE.json",
        "results/phase_lag_exploratory/k3e-12_sw094/PLP_EXP_K3E_12_SW094_WAVE",
        "results/phase_lag_exploratory/k3e-12_sw094/PLP_EXP_K3E_12_SW094_EQ/particle5000.vtp",
    ),
    Case(
        "k3_sw090",
        "$k=3\\times10^{-12}$, $S_w=0.90$\nsmooth",
        "configs/phase_lag_exploratory/k3e-12_sw090/02_WAVE.json",
        "results/phase_lag_exploratory/k3e-12_sw090/PLP_EXP_K3E_12_SW090_WAVE",
        "results/phase_lag_exploratory/k3e-12_sw090/PLP_EXP_K3E_12_SW090_EQ/particle5000.vtp",
    ),
    Case(
        "k1e13_sw094_raw",
        "$k=10^{-13}$, $S_w=0.94$\nno smoothing",
        "configs/phase_lag_exploratory/k1e-13_sw094_nosmooth/02_WAVE.json",
        "results/phase_lag_exploratory/k1e-13_sw094_nosmooth/PLP_EXP_K1E_13_SW094_NOSMOOTH_WAVE",
        "results/phase_lag_exploratory/k1e-13_sw094_nosmooth/PLP_EXP_K1E_13_SW094_NOSMOOTH_EQ/particle5000.vtp",
    ),
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_case(root: Path, case: Case, coordinates: np.ndarray) -> dict[str, Any]:
    config_path = root / case.config
    result_path = root / case.result
    checkpoint_path = root / case.checkpoint
    config = json.loads(config_path.read_text(encoding="utf-8"))
    _, checkpoint_arrays = read_vtp(checkpoint_path)
    raw_identifiers = np.asarray(
        checkpoint_arrays["ids"], dtype=np.float64
    ).reshape(-1)
    identifiers = np.rint(raw_identifiers).astype(np.int64)
    if not np.all(np.isfinite(raw_identifiers)) or not np.array_equal(
        raw_identifiers, identifiers.astype(np.float64)
    ):
        raise ValueError(f"{case.key}: checkpoint IDs are invalid")
    order = np.argsort(identifiers)
    expected_ids = identifiers[order]
    if not np.array_equal(expected_ids, np.arange(coordinates.shape[0])):
        raise ValueError(f"{case.key}: checkpoint IDs are not 0..N-1")
    checkpoint_pressure = np.asarray(
        checkpoint_arrays["PIC_liquid_pressures"], dtype=np.float64
    )[order]
    if checkpoint_pressure.size != coordinates.shape[0] or not np.all(
        np.isfinite(checkpoint_pressure)
    ):
        raise ValueError(f"{case.key}: checkpoint pressure is invalid")
    files = sorted(result_path.glob("particle*.vtp"), key=particle_step)
    expected_steps = expected_particle_steps(config)
    actual_steps = [particle_step(path) for path in files]
    if actual_steps != expected_steps:
        raise ValueError(
            f"{case.key}: expected steps {expected_steps}, got {actual_steps}"
        )
    dt = float(config["analysis"]["dt"])
    times = np.asarray([particle_step(path) * dt for path in files])
    fit_indices = np.flatnonzero(
        (times >= FIT_START_S - 1.0e-12) & (times <= FIT_END_S + 1.0e-12)
    )
    pressure_rows = []
    for index in fit_indices:
        _, pressure = align_liquid_pressure(files[index], expected_ids)
        pressure_rows.append(pressure - checkpoint_pressure)
    fitted_values = np.asarray(pressure_rows)
    _, amplitude, phase = harmonic_fields(
        times[fit_indices], fitted_values, WAVE_PERIOD_S
    )
    r_squared = harmonic_r_squared(
        times[fit_indices], fitted_values, WAVE_PERIOD_S
    )
    surface_amplitude = same_column_surface_values(coordinates, amplitude)
    surface_phase = same_column_surface_values(coordinates, phase)
    amplitude_ratio = amplitude / np.maximum(
        surface_amplitude, np.finfo(float).eps
    )
    raw_phase = wrap_degrees(phase - surface_phase)
    eligible = (
        (amplitude_ratio >= MINIMUM_LOCAL_AMPLITUDE_RATIO)
        & (r_squared >= MINIMUM_HARMONIC_R_SQUARED)
    )
    phase_difference = raw_phase.copy()
    phase_difference[~eligible] = np.nan

    snapshot_index = int(np.argmin(np.abs(times - SNAPSHOT_TIME_S)))
    if abs(times[snapshot_index] - SNAPSHOT_TIME_S) > 0.5 * dt:
        raise ValueError(f"{case.key}: no exact {SNAPSHOT_TIME_S:g} s frame")
    snapshot_points, snapshot_pressure = align_liquid_pressure(
        files[snapshot_index], expected_ids
    )
    snapshot_excess = snapshot_pressure - checkpoint_pressure

    adjacent_differences: list[float] = []
    for elevation in np.unique(np.round(snapshot_points[:, 1], 5)):
        row = np.flatnonzero(
            np.isclose(snapshot_points[:, 1], elevation, atol=2.0e-4)
        )
        row = row[np.argsort(snapshot_points[row, 0])]
        if row.size > 1:
            dx = np.diff(snapshot_points[row, 0])
            dp = np.abs(np.diff(snapshot_pressure[row]))
            adjacent_differences.extend(dp[dx < 0.011].tolist())
    adjacent = np.asarray(adjacent_differences)
    absolute_phase = np.abs(raw_phase[eligible])
    probes: dict[str, Any] = {}
    for probe in build_probes(config, coordinates):
        probes[probe.name] = {
            "amplitude_pa": float(amplitude[probe.index]),
            "amplitude_ratio": float(amplitude_ratio[probe.index]),
            "harmonic_r_squared": float(r_squared[probe.index]),
            "same_column_phase_deg": (
                float(phase_difference[probe.index])
                if np.isfinite(phase_difference[probe.index])
                else None
            ),
        }
    return {
        "case": case,
        "config": config,
        "config_path": config_path,
        "result_path": result_path,
        "checkpoint_path": checkpoint_path,
        "snapshot_path": files[snapshot_index],
        "snapshot_points": snapshot_points,
        "snapshot_excess": snapshot_excess,
        "amplitude_ratio": amplitude_ratio,
        "phase_difference": phase_difference,
        "eligible": eligible,
        "metrics": {
            "eligible_particle_count": int(np.count_nonzero(eligible)),
            "eligible_fraction": float(np.mean(eligible)),
            "absolute_phase_quantiles_deg": {
                key: float(value)
                for key, value in zip(
                    ("q50", "q75", "q90", "q95", "max"),
                    np.percentile(absolute_phase, (50, 75, 90, 95, 100)),
                )
            },
            "fraction_eligible_abs_phase_ge_20deg": float(
                np.mean(absolute_phase >= 20.0)
            ),
            "fraction_eligible_abs_phase_ge_30deg": float(
                np.mean(absolute_phase >= 30.0)
            ),
            "adjacent_snapshot_abs_dp_q95_pa": float(
                np.percentile(adjacent, 95)
            ),
            "adjacent_snapshot_abs_dp_q99_pa": float(
                np.percentile(adjacent, 99)
            ),
            "probes": probes,
        },
    }


def main() -> int:
    root = Path(__file__).resolve().parent
    output = root / "analysis/phase_lag_exploratory/comparison"
    output.mkdir(parents=True, exist_ok=True)
    coordinates = read_initial_coordinates(root / "particles.txt")
    loaded = [load_case(root, case, coordinates) for case in CASES]
    pressure_limit = max(
        float(np.max(np.abs(item["snapshot_excess"]))) for item in loaded
    )
    amplitude_limit = max(
        1.0, max(float(np.nanpercentile(item["amplitude_ratio"], 99.5)) for item in loaded)
    )
    figure, axes = plt.subplots(
        3, len(loaded), figsize=(18.0, 9.1), constrained_layout=True
    )
    pressure_artist = phase_artist = amplitude_artist = None
    for column, item in enumerate(loaded):
        center = np.asarray(
            item["config"]["analysis"]["rigid_pipeline"]["center"], dtype=float
        )
        pressure_artist = draw_field(
            axes[0, column],
            item["snapshot_points"],
            item["snapshot_excess"] / 1000.0,
            item["config"],
            center,
            title=item["case"].title,
            cmap="RdBu_r",
            limits=(-pressure_limit / 1000.0, pressure_limit / 1000.0),
        )
        phase_artist = draw_field(
            axes[1, column],
            coordinates,
            item["phase_difference"],
            item["config"],
            center,
            title=(
                "Coherent same-column phase\n"
                f"|phase|≥20°: "
                f"{100*item['metrics']['fraction_eligible_abs_phase_ge_20deg']:.1f}%"
            ),
            cmap="twilight_shifted",
            limits=(-180.0, 180.0),
        )
        amplitude_artist = draw_field(
            axes[2, column],
            coordinates,
            item["amplitude_ratio"],
            item["config"],
            center,
            title=(
                "Amplitude / same-column surface\n"
                f"snapshot roughness q95="
                f"{item['metrics']['adjacent_snapshot_abs_dp_q95_pa']:.0f} Pa"
            ),
            cmap="viridis",
            limits=(0.0, amplitude_limit),
        )
    assert pressure_artist is not None and phase_artist is not None
    assert amplitude_artist is not None
    figure.colorbar(
        pressure_artist,
        ax=axes[0, :],
        label=f"Liquid excess pressure at t={SNAPSHOT_TIME_S:g} s (kPa)",
        shrink=0.88,
    )
    figure.colorbar(
        phase_artist,
        ax=axes[1, :],
        label="Same-column phase difference (degrees)",
        shrink=0.88,
    )
    figure.colorbar(
        amplitude_artist,
        ax=axes[2, :],
        label="Fundamental amplitude ratio",
        shrink=0.88,
    )
    figure.suptitle(
        "Two-dimensional liquid-pressure phase-lag parameter comparison\n"
        "fit 1.3–3.9 s; phase requires amplitude ≥5% and harmonic $R^2≥0.8$; "
        "all panels use common row-wise colour scales",
        fontsize=14,
    )
    figure_path = output / "liquid_pressure_phase_lag_parameter_comparison"
    figure.savefig(figure_path.with_suffix(".png"), dpi=300, bbox_inches="tight")
    figure.savefig(figure_path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(figure)

    audit = {
        "schema": SCHEMA,
        "snapshot_time_s": SNAPSHOT_TIME_S,
        "fit_window_s": [FIT_START_S, FIT_END_S],
        "minimum_amplitude_ratio": MINIMUM_LOCAL_AMPLITUDE_RATIO,
        "minimum_harmonic_r_squared": MINIMUM_HARMONIC_R_SQUARED,
        "cases": {
            item["case"].key: {
                "title": item["case"].title.replace("\n", " "),
                "config": {
                    "path": str(item["config_path"].resolve()),
                    "sha256": sha256(item["config_path"]),
                },
                "checkpoint": {
                    "path": str(item["checkpoint_path"].resolve()),
                    "sha256": sha256(item["checkpoint_path"]),
                },
                "snapshot": {
                    "path": str(item["snapshot_path"].resolve()),
                    "sha256": sha256(item["snapshot_path"]),
                },
                "metrics": item["metrics"],
            }
            for item in loaded
        },
    }
    (output / "liquid_pressure_phase_lag_parameter_comparison.audit.json").write_text(
        json.dumps(audit, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps({key: value["metrics"] for key, value in audit["cases"].items()}, indent=2))
    print(f"Wrote common-scale comparison to {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
