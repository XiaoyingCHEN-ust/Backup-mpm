#!/usr/bin/env python3
"""Gate the SANISAND stage on a physically valid hydrostatic checkpoint."""

import argparse
import json
from pathlib import Path

import numpy as np

import analyze_pipeline_wall as wall_tools


def vector_norm(values):
    return np.linalg.norm(np.asarray(values)[:, :2], axis=1)


def validate_frame(path, points, arrays, required):
    """Reject incomplete or non-finite VTP data before computing QA metrics."""
    points = np.asarray(points)
    if points.ndim != 2 or points.shape[1] < 2 or not np.all(np.isfinite(points)):
        raise RuntimeError(f"{path.name} contains invalid particle coordinates")
    missing = required - arrays.keys()
    if missing:
        raise RuntimeError(f"{path.name} is missing {sorted(missing)}")
    scalar_fields = {"liquid_pressures", "liquid_saturations"}
    for field in required:
        values = np.asarray(arrays[field])
        if values.shape[0] != points.shape[0] or not np.all(np.isfinite(values)):
            raise RuntimeError(f"{path.name} contains invalid {field}")
        if field in scalar_fields and values.ndim != 1:
            raise RuntimeError(f"{path.name} has malformed scalar field {field}")
        if field not in scalar_fields and (values.ndim != 2 or values.shape[1] < 2):
            raise RuntimeError(f"{path.name} has malformed vector field {field}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=Path("mpm-initial.json"))
    args = parser.parse_args()

    config = json.loads(args.config.read_text(encoding="utf-8"))
    analysis = config["analysis"]
    result = Path(config["post_processing"]["path"]) / analysis["uuid"]
    final_step = int(analysis["nsteps"])
    initial_path = result / "particle00000.vtp"
    final_path = result / f"particle{final_step:05d}.vtp"
    pipeline_path = result / f"pipeline{final_step:05d}.vtp"
    for path in (initial_path, final_path, pipeline_path):
        if not path.is_file() or path.stat().st_size == 0:
            raise RuntimeError(f"Required initialization result is missing: {path}")

    wall_tools.WANTED_ARRAYS.update(
        {"stresses", "liquid_pressures", "liquid_saturations"}
    )
    points0, arrays0 = wall_tools.read_vtp(initial_path)
    points, arrays = wall_tools.read_vtp(final_path)
    required = {
        "stresses",
        "liquid_pressures",
        "liquid_saturations",
        "liquid_velocities",
        "gas_velocities",
        "velocities",
        "displacements",
    }
    validate_frame(initial_path, points0, arrays0, required)
    validate_frame(final_path, points, arrays, required)

    expected0 = 1000.0 * 9.81 * np.maximum(1.0 - points0[:, 1], 0.0)
    pressure_error0 = arrays0["liquid_pressures"] - expected0
    sigma_yy = arrays["stresses"][:, 1]
    compressive_fraction = float(np.mean(sigma_yy <= 1.0e-6))

    pipe = analysis["rigid_pipeline"]
    centre, pipe_velocity = wall_tools.pipeline_state(pipeline_path)
    radius = np.linalg.norm(points[:, :2] - centre, axis=1)
    influence_radius = float(pipe["outer_radius"]) + float(
        pipe["phase_no_flux_band_factor"]
    ) * float(pipe["particle_radius"])
    mask = radius <= influence_radius + 1.0e-12
    if not np.any(mask) or np.any(radius[mask] <= 1.0e-12):
        raise RuntimeError("No valid particle cohort was found at the pipeline wall")
    normal = (points[mask, :2] - centre) / radius[mask, None]
    liquid_normal = np.sum(
        (arrays["liquid_velocities"][mask, :2] - pipe_velocity) * normal,
        axis=1,
    )
    gas_normal = np.sum(
        (arrays["gas_velocities"][mask, :2] - pipe_velocity) * normal,
        axis=1,
    )

    expected_final = 1000.0 * 9.81 * np.maximum(1.0 - points[:, 1], 0.0)
    final_pressure_rmse = float(
        np.sqrt(np.mean((arrays["liquid_pressures"] - expected_final) ** 2))
    )
    target_saturation = float(config["materials"][1]["liquid_saturation"])
    wall_saturation = arrays["liquid_saturations"][mask]
    metrics = {
        "schema": "pipeline-initial-state-qa-v1",
        "result": str(result),
        "particle_count": int(points.shape[0]),
        "initial_liquid_pressure_max_error_pa": float(
            np.max(np.abs(pressure_error0))
        ),
        "initial_liquid_pressure_rmse_pa": float(
            np.sqrt(np.mean(pressure_error0**2))
        ),
        "final_liquid_pressure_rmse_from_hydrostatic_pa": final_pressure_rmse,
        "final_vertical_effective_stress_mean_pa": float(np.mean(sigma_yy)),
        "final_vertical_effective_stress_min_pa": float(np.min(sigma_yy)),
        "final_vertical_effective_stress_max_pa": float(np.max(sigma_yy)),
        "final_compressive_particle_fraction": compressive_fraction,
        "final_max_solid_speed_m_s": float(np.max(vector_norm(arrays["velocities"]))),
        "final_max_solid_displacement_m": float(
            np.max(vector_norm(arrays["displacements"]))
        ),
        "wall_particle_count": int(np.sum(mask)),
        "wall_max_abs_liquid_normal_velocity_m_s": float(
            np.max(np.abs(liquid_normal))
        ),
        "wall_max_abs_gas_normal_velocity_m_s": float(np.max(np.abs(gas_normal))),
        "wall_saturation_target": target_saturation,
        "wall_saturation_std": float(np.std(wall_saturation)),
        "wall_saturation_min": float(np.min(wall_saturation)),
        "wall_saturation_max": float(np.max(wall_saturation)),
        "wall_saturation_max_abs_target_error": float(
            np.max(np.abs(wall_saturation - target_saturation))
        ),
    }
    failures = []
    if metrics["initial_liquid_pressure_max_error_pa"] > 50.0:
        failures.append("the input hydrostatic liquid pressure was not installed")
    if compressive_fraction < 0.90:
        failures.append("fewer than 90% of particles retain compressive vertical stress")
    if metrics["final_vertical_effective_stress_mean_pa"] >= -100.0:
        failures.append("mean final vertical effective stress is not compressive")
    if metrics["final_max_solid_speed_m_s"] > 5.0e-3:
        failures.append("initialization has not reached a low solid velocity")
    if final_pressure_rmse > 500.0:
        failures.append("the final liquid-pressure field is not near hydrostatic")
    if metrics["wall_saturation_max_abs_target_error"] > 0.02:
        failures.append("pipeline-wall saturation drifted too far from its target")
    if max(
        metrics["wall_max_abs_liquid_normal_velocity_m_s"],
        metrics["wall_max_abs_gas_normal_velocity_m_s"],
    ) > 1.0e-6:
        failures.append("the fixed pipeline wall is not phase-impermeable")
    metrics["passed"] = not failures
    metrics["failures"] = failures

    output = result / "initial_state_qa.json"
    output.write_text(json.dumps(metrics, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(metrics, indent=2))
    if failures:
        raise RuntimeError("Initialization QA failed: " + "; ".join(failures))


if __name__ == "__main__":
    main()
