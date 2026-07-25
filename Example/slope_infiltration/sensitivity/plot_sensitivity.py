#!/usr/bin/env python3
"""Create the two-panel sensitivity figure for the manuscript."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ONSET_COLUMNS = (
    "onset_at_0p15m_s",
    "onset_s",
    "onset_at_0p25m_s",
)


def load_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def select_one(
    rows: list[dict[str, str]],
    *,
    permeability_factor: float,
    gas_model: str,
    head_pa: int,
    surface_gas_pressure_ratio: float | None = None,
) -> dict[str, str]:
    candidates = [
        row
        for row in rows
        if float(row["permeability_factor"]) == permeability_factor
        and row["gas_model"] == gas_model
        and int(row["head_pa"]) == head_pa
        and (
            surface_gas_pressure_ratio is None
            or float(row["surface_gas_pressure_ratio"])
            == surface_gas_pressure_ratio
        )
    ]
    if len(candidates) != 1:
        raise ValueError(
            f"Expected one matching row, found {len(candidates)}: "
            f"k={permeability_factor}, model={gas_model}, head={head_pa}, "
            f"ratio={surface_gas_pressure_ratio}"
        )
    return candidates[0]


def delta_triplet(fixed: dict[str, str], variable: dict[str, str]) -> np.ndarray:
    if any(not fixed[column] or not variable[column] for column in ONSET_COLUMNS):
        return np.full(3, np.nan)
    return np.array(
        [
            float(fixed[column]) - float(variable[column])
            for column in ONSET_COLUMNS
        ]
    )


def plot_curve(
    axis,
    heads_m: np.ndarray,
    values: np.ndarray,
    *,
    label: str,
    marker: str,
    color,
) -> None:
    central = values[:, 1]
    lower = np.nanmin(values, axis=1)
    upper = np.nanmax(values, axis=1)
    axis.plot(
        heads_m,
        central,
        marker=marker,
        linewidth=1.8,
        markersize=5,
        label=label,
        color=color,
    )
    axis.fill_between(heads_m, lower, upper, color=color, alpha=0.12)


def main() -> None:
    script_dir = Path(__file__).resolve().parent
    rows = load_rows(script_dir / "onset_times.csv")
    heads_pa = np.array(sorted({int(row["head_pa"]) for row in rows}))
    heads_m = heads_pa / 10000.0

    figure, axes = plt.subplots(1, 2, figsize=(10.2, 4.0), sharey=True)
    colors = plt.get_cmap("viridis")(np.linspace(0.15, 0.85, 3))

    summary_rows: list[dict[str, float | str]] = []
    for color, factor, marker in zip(colors, (0.1, 1.0, 10.0), ("o", "s", "^")):
        values = []
        for head_pa in heads_pa:
            fixed = select_one(
                rows,
                permeability_factor=factor,
                gas_model="fixed",
                head_pa=int(head_pa),
            )
            variable = select_one(
                rows,
                permeability_factor=factor,
                gas_model="variable",
                surface_gas_pressure_ratio=1.0,
                head_pa=int(head_pa),
            )
            delta = delta_triplet(fixed, variable)
            values.append(delta)
            summary_rows.append(
                {
                    "panel": "permeability",
                    "series": f"k/k0={factor:g}",
                    "head_m": head_pa / 10000.0,
                    "delta_t_s": delta[1],
                    "delta_t_0p15m_s": delta[0],
                    "delta_t_0p25m_s": delta[2],
                }
            )
        plot_curve(
            axes[0],
            heads_m,
            np.asarray(values),
            label=f"$k/k_0={factor:g}$",
            marker=marker,
            color=color,
        )

    for color, ratio, marker in zip(colors, (0.0, 0.5, 1.0), ("o", "s", "^")):
        values = []
        for head_pa in heads_pa:
            fixed = select_one(
                rows,
                permeability_factor=1.0,
                gas_model="fixed",
                head_pa=int(head_pa),
            )
            variable = select_one(
                rows,
                permeability_factor=1.0,
                gas_model="variable",
                surface_gas_pressure_ratio=ratio,
                head_pa=int(head_pa),
            )
            delta = delta_triplet(fixed, variable)
            values.append(delta)
            summary_rows.append(
                {
                    "panel": "gas_boundary",
                    "series": f"Pg/Pl={ratio:g}",
                    "head_m": head_pa / 10000.0,
                    "delta_t_s": delta[1],
                    "delta_t_0p15m_s": delta[0],
                    "delta_t_0p25m_s": delta[2],
                }
            )
        plot_curve(
            axes[1],
            heads_m,
            np.asarray(values),
            label=f"$P_{{g,b}}/P_{{l,b}}={ratio:g}$",
            marker=marker,
            color=color,
        )

    for axis in axes:
        axis.axhline(0.0, color="0.45", linewidth=1.0, linestyle="--")
        axis.set_xlabel("Surface water head (m)")
        axis.set_xlim(heads_m.min(), heads_m.max())
        axis.grid(alpha=0.18)
        axis.legend(frameon=False, fontsize=9)
    axes[0].set_ylabel(
        r"Onset-time difference, $\Delta t=t_{\rm fixed}-t_{\rm variable}$ (s)"
    )
    axes[0].set_title("(a) Intrinsic permeability")
    axes[1].set_title("(b) Surface gas boundary")
    figure.tight_layout()

    figure_dir = script_dir / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)
    figure.savefig(figure_dir / "sensitivity_delta_t.png", dpi=400)
    figure.savefig(figure_dir / "sensitivity_delta_t.pdf")

    summary_path = script_dir / "sensitivity_summary.csv"
    with summary_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(summary_rows[0]))
        writer.writeheader()
        writer.writerows(summary_rows)
    print(f"Wrote figure files to {figure_dir}")
    print(f"Wrote paired differences to {summary_path}")


if __name__ == "__main__":
    main()
