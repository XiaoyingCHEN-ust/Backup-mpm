"""Extract wetting-front and pore-air-pressure histories from MPM VTP output."""

from __future__ import annotations

import argparse
import csv
import json
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy


CASE_DIR = Path(__file__).resolve().parent
STEP_PATTERN = re.compile(r"particle(\d+)\.vtp$")


def read_point_array(polydata, name: str):
    array = polydata.GetPointData().GetArray(name)
    if array is None:
        raise KeyError(f"VTP array '{name}' was not found")
    return vtk_to_numpy(array)


def read_vtp(path: Path):
    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(str(path))
    reader.Update()
    return reader.GetOutput()


def layer_average(y: np.ndarray, values: np.ndarray):
    rounded_y = np.round(y, decimals=8)
    layers = np.unique(rounded_y)
    averages = np.array([values[rounded_y == layer].mean() for layer in layers])
    return layers, averages


def extract_case(case_name: str, saturation_threshold: float):
    config_path = CASE_DIR / f"mpm_{case_name}.json"
    with config_path.open(encoding="utf-8") as stream:
        config = json.load(stream)

    dt = float(config["analysis"]["dt"])
    uuid = config["analysis"]["uuid"]
    result_dir = CASE_DIR / config["post_processing"]["path"] / uuid
    files = sorted(result_dir.glob("particle*.vtp"))
    if not files:
        raise FileNotFoundError(f"No particle VTP files found in {result_dir}")

    height = 94 * float(config["mesh"]["cellsize_min"])
    cell_size = float(config["mesh"]["cellsize_min"])
    rows = []

    for path in files:
        match = STEP_PATTERN.search(path.name)
        if match is None:
            continue
        step = int(match.group(1))
        polydata = read_vtp(path)
        points = vtk_to_numpy(polydata.GetPoints().GetData())
        y = points[:, 1]
        saturation = read_point_array(polydata, "liquid_saturations").ravel()
        gas_pressure = read_point_array(polydata, "gas_pressures").ravel()

        layers, layer_saturation = layer_average(y, saturation)
        wetted_layers = layers[layer_saturation >= saturation_threshold]
        if wetted_layers.size:
            front_y = float(np.min(wetted_layers))
            wetting_depth = float(np.clip(height - front_y, 0.0, height))
        else:
            front_y = height
            wetting_depth = 0.0

        dry_mask = (saturation < 0.10) & (y < front_y - 0.5 * cell_size)
        if np.any(dry_mask):
            dry_air = gas_pressure[dry_mask] / 1000.0
            air_mean = float(np.mean(dry_air))
            air_p10 = float(np.percentile(dry_air, 10.0))
            air_p90 = float(np.percentile(dry_air, 90.0))
        else:
            air_mean = air_p10 = air_p90 = float("nan")

        transmission_mask = (
            (saturation >= saturation_threshold)
            & (y < height - cell_size)
            & (y > front_y + cell_size)
        )
        transmission_saturation = (
            float(np.mean(saturation[transmission_mask]))
            if np.any(transmission_mask)
            else float("nan")
        )

        rows.append(
            {
                "step": step,
                "time_s": step * dt,
                "wetting_depth_m": wetting_depth,
                "front_y_m": front_y,
                "dry_air_pressure_mean_kpa": air_mean,
                "dry_air_pressure_p10_kpa": air_p10,
                "dry_air_pressure_p90_kpa": air_p90,
                "transmission_saturation": transmission_saturation,
            }
        )

    output_dir = CASE_DIR / "validation_results"
    output_dir.mkdir(exist_ok=True)
    output_csv = output_dir / f"simulation_summary_{case_name}.csv"
    with output_csv.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return rows, output_csv


def load_reference(path: Path):
    return np.genfromtxt(path, delimiter=",", names=True)


def plot_comparison(summaries):
    reference_dir = CASE_DIR / "reference_data"
    output_dir = CASE_DIR / "validation_results"
    fig, axes = plt.subplots(1, 2, figsize=(10.0, 4.0), constrained_layout=True)

    styles = {"open": ("tab:blue", "o"), "closed": ("tab:red", "s")}
    for case_name, rows in summaries.items():
        color, marker = styles[case_name]
        reference = load_reference(reference_dir / f"wetting_front_{case_name}.csv")
        axes[0].plot(
            reference["time_s"],
            reference["wetting_depth_from_surface_mm"] / 1000.0,
            marker=marker,
            linestyle="none",
            markerfacecolor="none",
            color=color,
            label=f"Experiment ({case_name})",
        )
        axes[0].plot(
            [row["time_s"] for row in rows],
            [row["wetting_depth_m"] for row in rows],
            color=color,
            label=f"MPM ({case_name})",
        )

    axes[0].set_xlabel("Time (s)")
    axes[0].set_ylabel("Wetting-front depth (m)")
    axes[0].set_xlim(left=0.0)
    axes[0].set_ylim(bottom=0.0)
    axes[0].legend(frameon=False, fontsize=8)

    envelope = load_reference(reference_dir / "closed_air_pressure_envelope.csv")
    axes[1].fill_between(
        envelope["time_s"],
        envelope["lower_kpa"],
        envelope["upper_kpa"],
        color="0.85",
        label="Measured range (closed)",
    )
    for case_name, rows in summaries.items():
        color, _ = styles[case_name]
        axes[1].plot(
            [row["time_s"] for row in rows],
            [row["dry_air_pressure_mean_kpa"] for row in rows],
            color=color,
            label=f"MPM ({case_name})",
        )
    axes[1].set_xlabel("Time (s)")
    axes[1].set_ylabel("Mean pore-air pressure ahead of front (kPa)")
    axes[1].set_xlim(left=0.0)
    axes[1].set_ylim(bottom=0.0)
    axes[1].legend(frameon=False, fontsize=8)

    output = output_dir / "siemens2013_validation.png"
    fig.savefig(output, dpi=300)
    plt.close(fig)
    return output


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--cases", nargs="+", choices=("open", "closed"), default=("open", "closed")
    )
    parser.add_argument("--saturation-threshold", type=float, default=0.40)
    args = parser.parse_args()

    summaries = {}
    for case_name in args.cases:
        rows, output_csv = extract_case(case_name, args.saturation_threshold)
        summaries[case_name] = rows
        print(f"Wrote {output_csv}")
    print(f"Wrote {plot_comparison(summaries)}")


if __name__ == "__main__":
    main()
