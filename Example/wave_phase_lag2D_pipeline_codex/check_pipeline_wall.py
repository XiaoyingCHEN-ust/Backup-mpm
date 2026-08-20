#!/usr/bin/env python3
"""Fail fast on stale geometry or an invalid moving pipe no-flux setup."""

import json
import math
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parent
CONFIGS = (ROOT / "mpm-initial.json", ROOT / "mpm-3p.json")
TARGET_SATURATION = 0.94
TARGET_DYNAMIC_PERMEABILITY = 9.79e-12
MAX_EQUILIBRIUM_PERMEABILITY_FACTOR = 5.0
TARGET_THREADS = 32


def load(path):
    with path.open(encoding="utf-8") as stream:
        return json.load(stream)


def read_particles(path):
    with path.open(encoding="utf-8") as stream:
        count = int(stream.readline())
        points = np.loadtxt(stream, ndmin=2)
    if points.shape != (count, 2) or not np.all(np.isfinite(points)):
        raise RuntimeError("particles.txt is malformed")
    return points


def read_scalars(path, count):
    with path.open(encoding="utf-8") as stream:
        declared = int(stream.readline())
        values = np.loadtxt(stream, ndmin=1)
    if declared != count or values.shape != (count,) or not np.all(np.isfinite(values)):
        raise RuntimeError(f"{path.name} is malformed")
    return values


def check_batch(path):
    text = path.read_text(encoding="utf-8")
    required = (
        "#SBATCH --partition=granularmech",
        "#SBATCH --account=comgranmech",
        f"#SBATCH --cpus-per-task={TARGET_THREADS}",
    )
    missing = [line for line in required if line not in text]
    if missing:
        raise RuntimeError(f"{path.name}: missing {missing}")


def main():
    source_root = ROOT.parent.parent / "mpm"
    source_checks = (
        (
            source_root / "include/particles/particle_threephase_lag.h",
            "has_input_initial_liquid_pressure_",
        ),
        (
            source_root / "include/particles/particle_threephase_lag.tcc",
            "input_initial_liquid_pressure_",
        ),
        (
            source_root / "include/particles/particle_threephase_lag.tcc",
            "registered capillary suction and target saturation",
        ),
        (
            source_root / "include/particles/particle_threephase_lag.tcc",
            "apply_rigid_circle_phase_no_flux",
        ),
        (
            source_root / "include/solvers/thm_mpm_explicit_threephase_lag.tcc",
            "apply_rigid_pipeline_phase_no_flux",
        ),
        (
            source_root / "include/particles/particle_threephase_lag.tcc",
            "apply_submerged_surface_pressure_traction",
        ),
        (
            source_root / "include/solvers/thm_mpm_explicit_threephase_lag.tcc",
            "submerged_surface_pressure_traction_",
        ),
    )
    for source, marker in source_checks:
        if not source.is_file():
            raise RuntimeError(f"Expected solver source file is missing: {source}")
        if marker not in source.read_text(encoding="utf-8", errors="replace"):
            raise RuntimeError(
                f"Required solver marker {marker!r} is absent from {source}; "
                "run bash apply_mpm_patch.sh and rebuild MPM"
            )

    geometry = load(ROOT / "pipeline_geometry.json")
    entity_sets = load(ROOT / "entity_sets.json")
    node_sets = {entry["id"]: entry["set"] for entry in entity_sets["node_sets"]}
    if 4 in node_sets:
        raise RuntimeError("Stale stair-step pipeline node set 4 is still present")
    if (ROOT / "pipeline_wall_euler_angles.txt").exists():
        raise RuntimeError("Stale pipeline_wall_euler_angles.txt is still present")

    particles = read_particles(ROOT / "particles.txt")
    initial_pressures = read_scalars(
        ROOT / "initial_liquid_pressures.txt", particles.shape[0]
    )
    expected_pressures = 1000.0 * 9.81 * np.maximum(1.0 - particles[:, 1], 0.0)
    if not np.allclose(initial_pressures, expected_pressures, atol=1.0e-8, rtol=0.0):
        raise RuntimeError("initial_liquid_pressures.txt is not hydrostatic")
    particle_sets = {
        entry["id"]: np.asarray(entry["set"], dtype=int)
        for entry in entity_sets["particle_sets"]
    }
    summary = load(ROOT / "generation_summary.json")
    particle_spacing = float(summary["particle_spacing"])
    surface_y = 0.5
    surface_traction_ids = np.flatnonzero(
        np.abs(
            (surface_y - particles[:, 1]) - 0.5 * particle_spacing
        )
        <= 0.05 * particle_spacing
    )
    free_surface_ids = particle_sets.get(0, np.empty(0, dtype=int))
    if surface_traction_ids.size == 0 or not set(surface_traction_ids).issubset(
        set(free_surface_ids)
    ):
        raise RuntimeError(
            "The uppermost matched-traction row is missing from particle set 0"
        )
    contact_ids = particle_sets.get(3, np.empty(0, dtype=int))
    if (
        contact_ids.size == 0
        or contact_ids.min() < 0
        or contact_ids.max() >= particles.shape[0]
    ):
        raise RuntimeError("Pipeline contact particle set 3 is invalid")

    configs = [load(path) for path in CONFIGS]
    phase_wall_ids = None
    for stage, path, config in zip(("equilibrium", "dynamic"), CONFIGS, configs):
        mesh = config["mesh"]
        bc = mesh["boundary_conditions"]
        if "nodal_euler_angles" in bc:
            raise RuntimeError(f"{path.name}: stale nodal Euler wall is enabled")
        if any(item["nset_id"] == 4 for item in bc["velocity_constraints"]):
            raise RuntimeError(f"{path.name}: stale node-set-4 wall is enabled")

        soil, fluid = config["materials"]
        permeability = float(soil["intrinsic_permeability"])
        if stage == "dynamic":
            if not np.isclose(
                permeability,
                TARGET_DYNAMIC_PERMEABILITY,
                rtol=1.0e-12,
                atol=0.0,
            ):
                raise RuntimeError(
                    f"{path.name}: target intrinsic permeability must be "
                    f"{TARGET_DYNAMIC_PERMEABILITY:.3e} m^2"
                )
        elif not (
            TARGET_DYNAMIC_PERMEABILITY
            <= permeability
            <= MAX_EQUILIBRIUM_PERMEABILITY_FACTOR
            * TARGET_DYNAMIC_PERMEABILITY
        ):
            raise RuntimeError(
                f"{path.name}: equilibrium permeability must remain between "
                "1x and 5x the dynamic value"
            )
        if not np.isclose(fluid["liquid_saturation"], TARGET_SATURATION):
            raise RuntimeError(f"{path.name}: liquid saturation drifted")
        if not np.isclose(
            fluid["liquid_saturation"] + fluid["gas_saturation"], 1.0
        ):
            raise RuntimeError(f"{path.name}: phase saturations do not sum to one")

        pipe = config["analysis"]["rigid_pipeline"]
        if not config["analysis"].get(
            "submerged_surface_pressure_traction", False
        ):
            raise RuntimeError(
                f"{path.name}: matching submerged surface traction is disabled"
            )
        if not pipe.get("phase_no_flux", False):
            raise RuntimeError(f"{path.name}: moving phase no-flux is disabled")
        factor = float(pipe.get("phase_no_flux_band_factor", 0.0))
        if not np.isclose(factor, math.sqrt(2.0), atol=1.0e-12):
            raise RuntimeError(f"{path.name}: unexpected no-flux band factor")
        if not np.allclose(pipe["center"], geometry["center"], atol=1.0e-12):
            raise RuntimeError(f"{path.name}: rigid-pipe centre disagrees with mesh")
        if not np.isclose(pipe["outer_radius"], geometry["outer_radius"]):
            raise RuntimeError(f"{path.name}: rigid-pipe radius disagrees with mesh")

        center = np.asarray(pipe["center"], dtype=float)
        radii = np.linalg.norm(particles - center, axis=1)
        influence_radius = (
            float(pipe["outer_radius"]) + factor * float(pipe["particle_radius"])
        )
        selected = np.flatnonzero(radii <= influence_radius + 1.0e-12)
        if selected.size < 12:
            raise RuntimeError(f"{path.name}: too few particles in no-flux band")
        if not set(selected).issubset(set(contact_ids)):
            raise RuntimeError(
                f"{path.name}: no-flux particles escape contact/free-surface set 3"
            )
        if phase_wall_ids is None:
            phase_wall_ids = selected
        elif not np.array_equal(phase_wall_ids, selected):
            raise RuntimeError("Initialization and dynamic no-flux cohorts differ")

    initial, dynamic = configs
    initial_pressure_file = initial["mesh"].get(
        "particles_pore_pressures", {}
    ).get("file")
    if initial_pressure_file != "initial_liquid_pressures.txt":
        raise RuntimeError("Initialization does not register the hydrostatic pressure file")
    if initial["materials"][0]["type"] != "LinearElastic2D":
        raise RuntimeError("Initialization must use LinearElastic2D")
    if initial["analysis"]["PIC"] != 1 or initial["analysis"].get("APIC", False):
        raise RuntimeError("Initialization must use pure PIC=1 with APIC disabled")
    if not initial["analysis"]["rigid_pipeline"]["fixed"]:
        raise RuntimeError("Initialization pipeline must be fixed")
    if dynamic["materials"][0]["type"] != "SANISAND2D":
        raise RuntimeError("Dynamic stage must use SANISAND2D")
    if dynamic["analysis"]["PIC"] != 0 or dynamic["analysis"].get("APIC", False):
        raise RuntimeError("Dynamic stage must use PIC=0 with APIC disabled")
    if "ids" in dynamic["post_processing"]["vtk"]:
        raise RuntimeError(
            "Dynamic VTK requests unsupported 'ids'; particle IDs are not "
            "required by analyze_pipeline_wall.py"
        )
    if dynamic["analysis"]["resume"]["uuid"] != initial["analysis"]["uuid"]:
        raise RuntimeError("Dynamic checkpoint UUID does not match initialization")
    if dynamic["analysis"]["resume"]["step"] != initial["analysis"]["nsteps"]:
        raise RuntimeError("Dynamic checkpoint step does not match initialization")

    for batch in (
        ROOT / "run_mpm_ini.sbatch",
        ROOT / "run_mpm.sbatch",
    ):
        check_batch(batch)

    print(
        "Validated moving impermeable pipeline case: "
        f"phase_wall_particles={phase_wall_ids.size}, "
        f"surface_traction_particles={surface_traction_ids.size}, "
        f"Sw={TARGET_SATURATION}, "
        f"kappa_eq={initial['materials'][0]['intrinsic_permeability']:.3e} m^2, "
        f"kappa_dynamic={TARGET_DYNAMIC_PERMEABILITY:.3e} m^2, "
        f"threads={TARGET_THREADS}"
    )


if __name__ == "__main__":
    main()
