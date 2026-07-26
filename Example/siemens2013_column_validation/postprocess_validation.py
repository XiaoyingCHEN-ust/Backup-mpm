"""Extract wetting-front and pore-air-pressure histories from MPM VTP output."""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy


CASE_DIR = Path(__file__).resolve().parent
STEP_PATTERN = re.compile(r"particle(\d+)\.vtp$")
EXPERIMENT_BOTTOM_ELEVATION_MM = 200.0


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


def top_connected_front(
    layers: np.ndarray,
    layer_saturation: np.ndarray,
    saturation_threshold: float,
    infiltration_origin_y: float,
):
    """Return the deepest wet layer connected continuously to the top.

    The open experiment also has a saturated constant-head layer at the base.
    Treating the deepest wet layer anywhere in the column as the front would
    therefore report full infiltration at the initial output.
    """
    order = np.argsort(layers)[::-1]
    layers_from_top = layers[order]
    wet_from_top = layer_saturation[order] >= saturation_threshold
    if not wet_from_top[0]:
        return infiltration_origin_y, 0.0

    first_dry = np.flatnonzero(~wet_from_top)
    if first_dry.size:
        connected_count = int(first_dry[0])
        if connected_count == 0:
            return infiltration_origin_y, 0.0
        front_y = float(
            0.5
            * (
                layers_from_top[connected_count - 1]
                + layers_from_top[connected_count]
            )
        )
    else:
        layer_spacing = float(np.median(np.abs(np.diff(layers_from_top))))
        front_y = max(0.0, float(layers_from_top[-1] - 0.5 * layer_spacing))
    return front_y, float(
        np.clip(infiltration_origin_y - front_y, 0.0, infiltration_origin_y)
    )


def load_ppt_locations(height: float):
    path = CASE_DIR / "reference_data" / "ppt_locations.csv"
    locations = []
    with path.open(encoding="utf-8", newline="") as stream:
        for row in csv.DictReader(stream):
            sensor = row["sensor"]
            model_y = (
                float(row["elevation_mm"]) - EXPERIMENT_BOTTOM_ELEVATION_MM
            ) / 1000.0
            locations.append((sensor, float(np.clip(model_y, 0.0, height))))
    return locations


def expected_initial_suction(config: dict):
    liquid = next(
        material for material in config["materials"] if "liquid_saturation" in material
    )
    if not liquid.get("initial_suction_from_swrc", False):
        return float(liquid.get("initial_suction", 137340.0))

    sw = float(liquid["liquid_saturation"])
    swr = float(liquid["liquid_saturation_res"])
    sgr = float(liquid["gas_saturation_res"])
    p0 = float(liquid["para_p0"])
    m = float(liquid["para_m"])
    effective = np.clip((sw - swr) / (1.0 - swr - sgr), 1.0e-12, 1.0 - 1.0e-12)
    return float(p0 * (effective ** (-1.0 / m) - 1.0) ** (1.0 - m))


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

    cell_size = float(config["mesh"]["cellsize_min"])
    height = float("nan")
    ppt_locations = []
    with (CASE_DIR / config["mesh"]["entity_sets"]).open(encoding="utf-8") as stream:
        entity_sets = json.load(stream)
    top_boundary_ids = np.asarray(
        next(
            item["set"]
            for item in entity_sets["particle_sets"]
            if item["id"] == 1
        ),
        dtype=int,
    )
    rows = []
    initial_dry_suction_median_pa = float("nan")
    infiltration_origin_y = float("nan")

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
        liquid_pressure = read_point_array(polydata, "liquid_pressures").ravel()
        suction_pressure = read_point_array(polydata, "suction_pressures").ravel()

        layers, layer_saturation = layer_average(y, saturation)
        _, layer_gas_pressure = layer_average(y, gas_pressure)
        if not rows:
            layer_spacing = float(np.min(np.diff(layers)))
            height = float(np.max(layers) + 0.5 * layer_spacing)
            ppt_locations = load_ppt_locations(height)
            infiltration_origin_y = float(
                np.min(y[top_boundary_ids]) - 0.5 * layer_spacing
            )
        front_y, wetting_depth = top_connected_front(
            layers,
            layer_saturation,
            saturation_threshold,
            infiltration_origin_y,
        )

        if not rows:
            initial_dry = saturation < 0.10
            if np.any(initial_dry):
                initial_dry_suction_median_pa = float(
                    np.median(suction_pressure[initial_dry])
                )

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

        row = {
            "step": step,
            "time_s": step * dt,
            "wetting_depth_m": wetting_depth,
            "front_y_m": front_y,
            "dry_air_pressure_mean_kpa": air_mean,
            "dry_air_pressure_p10_kpa": air_p10,
            "dry_air_pressure_p90_kpa": air_p90,
            "transmission_saturation": transmission_saturation,
            "saturation_min": float(np.min(saturation)),
            "saturation_max": float(np.max(saturation)),
            "max_abs_gas_pressure_kpa": float(
                np.max(np.abs(gas_pressure)) / 1000.0
            ),
            "max_abs_liquid_pressure_kpa": float(
                np.max(np.abs(liquid_pressure)) / 1000.0
            ),
        }
        for sensor, sensor_y in ppt_locations:
            layer_index = int(np.argmin(np.abs(layers - sensor_y)))
            sensor_saturation = float(layer_saturation[layer_index])
            key = sensor.lower()
            row[f"{key}_saturation"] = sensor_saturation
            row[f"{key}_gas_pressure_kpa"] = (
                float(layer_gas_pressure[layer_index] / 1000.0)
                if sensor_saturation < saturation_threshold
                else float("nan")
            )
        rows.append(row)

    output_dir = CASE_DIR / "validation_results"
    output_dir.mkdir(exist_ok=True)
    output_csv = output_dir / f"simulation_summary_{case_name}.csv"
    with output_csv.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    expected_suction = expected_initial_suction(config)
    checks = {
        "case": case_name,
        "output_files": len(rows),
        "expected_initial_suction_pa": expected_suction,
        "observed_initial_dry_suction_median_pa": initial_dry_suction_median_pa,
        "initialization_patch_active": math.isclose(
            initial_dry_suction_median_pa, expected_suction, rel_tol=0.02, abs_tol=1.0
        ),
        "finite_values": bool(
            all(
                np.isfinite(row[key])
                for row in rows
                for key in (
                    "saturation_min",
                    "saturation_max",
                    "max_abs_gas_pressure_kpa",
                    "max_abs_liquid_pressure_kpa",
                )
            )
        ),
        "saturation_in_unit_interval": bool(
            min(row["saturation_min"] for row in rows) >= 0.0
            and max(row["saturation_max"] for row in rows) <= 1.0
        ),
        "gas_sensor_samples": {
            sensor: sum(
                math.isfinite(row[f"{sensor.lower()}_gas_pressure_kpa"])
                for row in rows
            )
            for sensor, _ in ppt_locations
        },
        "gas_sensor_model_y_m": {
            sensor: sensor_y for sensor, sensor_y in ppt_locations
        },
        # A deliberately loose guard compared with the measured 1.2--1.6 kPa.
        # It catches numerical blow-up without being an agreement criterion.
        "pressures_below_safety_limit": bool(
            max(
                max(row["max_abs_gas_pressure_kpa"] for row in rows),
                max(row["max_abs_liquid_pressure_kpa"] for row in rows),
            )
            < 1000.0
        ),
    }
    checks["valid_for_comparison"] = bool(
        checks["initialization_patch_active"]
        and checks["finite_values"]
        and checks["saturation_in_unit_interval"]
        and checks["pressures_below_safety_limit"]
    )
    return rows, output_csv, checks


def load_reference(path: Path):
    return np.genfromtxt(path, delimiter=",", names=True)


def plot_comparison(summaries):
    reference_dir = CASE_DIR / "reference_data"
    output_dir = CASE_DIR / "validation_results"
    output_dir.mkdir(exist_ok=True)
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
    axes[0].set_title("(a) Wetting-front migration")
    axes[0].set_xlim(left=0.0)
    axes[0].set_ylim(bottom=0.0)
    axes[0].legend(frameon=False, fontsize=8)

    envelope = load_reference(reference_dir / "closed_air_pressure_envelope.csv")
    axes[1].fill_between(
        envelope["time_s"],
        envelope["lower_kpa"],
        envelope["upper_kpa"],
        color="0.85",
        label="Measured range (closed; six PPTs)",
    )
    if "closed" in summaries:
        closed_rows = summaries["closed"]
        sensors = [sensor for sensor, _ in load_ppt_locations(height=1.0716)]
        sensor_colors = plt.cm.Blues(np.linspace(0.35, 0.90, len(sensors)))
        sensor_styles = ("-", "--", "-.", ":", (0, (5, 2)), (0, (1, 1)))
        for sensor, color, line_style in zip(
            sensors, sensor_colors, sensor_styles
        ):
            pressure = np.asarray(
                [row[f"{sensor.lower()}_gas_pressure_kpa"] for row in closed_rows]
            )
            if not np.any(np.isfinite(pressure)):
                continue
            axes[1].plot(
                [row["time_s"] for row in closed_rows],
                pressure,
                color=color,
                linestyle=line_style,
                linewidth=1.0,
                label=f"MPM {sensor}",
            )
    if "open" in summaries:
        open_rows = summaries["open"]
        axes[1].plot(
            [row["time_s"] for row in open_rows],
            [row["dry_air_pressure_mean_kpa"] for row in open_rows],
            color="0.2",
            linestyle="--",
            linewidth=1.2,
            label="MPM open (dry-zone mean)",
        )
    axes[1].set_xlabel("Time (s)")
    axes[1].set_ylabel("Pore-air pressure at PPT elevations (kPa)")
    axes[1].set_title("(b) Pore-air pressure histories")
    axes[1].set_xlim(left=0.0)
    axes[1].set_ylim(bottom=0.0)
    axes[1].legend(frameon=False, fontsize=7, ncol=2)

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
    checks = {}
    for case_name in args.cases:
        rows, output_csv, case_checks = extract_case(
            case_name, args.saturation_threshold
        )
        summaries[case_name] = rows
        checks[case_name] = case_checks
        print(f"Wrote {output_csv}")

    status_path = CASE_DIR / "validation_results" / "validation_status.json"
    with status_path.open("w", encoding="utf-8") as stream:
        json.dump(checks, stream, indent=2)
    print(f"Wrote {status_path}")

    invalid = [name for name, status in checks.items() if not status["valid_for_comparison"]]
    if invalid:
        details = ", ".join(
            f"{name} (patch={checks[name]['initialization_patch_active']}, "
            f"pressure_safe={checks[name]['pressures_below_safety_limit']})"
            for name in invalid
        )
        raise RuntimeError(
            "Outputs are not valid for experimental comparison: " + details
        )

    print(f"Wrote {plot_comparison(summaries)}")


if __name__ == "__main__":
    main()
