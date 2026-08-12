#!/usr/bin/env python3
"""Plot manuscript-ready curve panels and register common field frames.

The script consumes only ``analyze_study.py`` outputs.  It creates three
curve-based composite figures corresponding to the Section 7.5 draft and a
JSON manifest identifying RL/RE VTP files at one common critical time for the
final field rendering.  It never chooses separate best-looking times per case.
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


IF_AREA = "area_support_ROI_upward_seepage_IF_ge_1_m2"
STRESS_AREA = "area_support_ROI_stress_loss_Rsigma_le_0p05_m2"
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
        -0.12,
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
            time, pressure = series(rows, f"{probe}_pressure_pa")
            axes[0, 0].plot(
                time,
                pressure,
                color=COLORS[code],
                linestyle=linestyle,
                linewidth=1.1,
                label=f"{code} {probe}",
            )
    axes[0, 0].set_xlabel("Time (s)")
    axes[0, 0].set_ylabel("Mixture pore pressure (Pa)")
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
        time, pressure = series(rows, "invert_pressure_pa")
        axes[0, 0].plot(time, pressure, color=COLORS[code], label=code)
        time, area = series(rows, IF_AREA)
        axes[0, 1].plot(time, area, color=COLORS[code], label=code)
        time, area = series(rows, STRESS_AREA)
        axes[1, 0].plot(time, area, color=COLORS[code], label=code)
    axes[0, 0].set_xlabel("Time (s)")
    axes[0, 0].set_ylabel("Invert mixture pressure (Pa)")
    axes[0, 0].legend()
    style_axis(axes[0, 0], "(a)")
    axes[0, 1].set_xlabel("Time (s)")
    axes[0, 1].set_ylabel("Support-zone area with $IF\\geq1$ (m$^2$)")
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


def common_field_frames(directory: Path) -> dict[str, Any]:
    rl_rows = probe_rows(directory, "RL")
    if IF_AREA not in rl_rows[0] and not any(IF_AREA in row for row in rl_rows):
        raise KeyError(f"RL history has no {IF_AREA}")
    critical = max(
        (row for row in rl_rows if IF_AREA in row),
        key=lambda row: row[IF_AREA],
    )
    target_time = critical["time_s"]
    output: dict[str, Any] = {
        "selection_rule": (
            "common physical time at the maximum RL support-zone area with IF>=1; "
            "the same time is used for RL and RE"
        ),
        "target_time_s": target_time,
        "RL_support_IF_area_m2": critical[IF_AREA],
        "cases": {},
    }
    for code in ("RL", "RE"):
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
        output["cases"][code] = {
            "path": str(chosen.resolve()),
            "step": chosen_step,
            "time_s": chosen_step * dt,
            "time_offset_from_target_s": chosen_step * dt - target_time,
        }
    return output


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--analysis-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    plot_hydraulic_bridge(args.analysis_dir, args.output_dir / "figure7_hydraulic_bridge")
    plot_phase_control(args.analysis_dir, args.output_dir / "figure8_phase_control")
    plot_constitutive(args.analysis_dir, args.output_dir / "figure9_constitutive")
    frames = common_field_frames(args.analysis_dir)
    (args.output_dir / "figure8_common_field_frames.json").write_text(
        json.dumps(frames, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(f"Wrote manuscript figure drafts to {args.output_dir.resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
