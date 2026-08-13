#!/usr/bin/env python3
"""Generate the pre-registered pipeline phase-lag study configurations.

The primary comparisons are deliberately sparse:

* LS vs HS: physical low/high phase-lag contrast (saturation only);
* HS vs HM: SANISAND/Mohr-Coulomb constitutive ablation;
* RL vs RE: matched one-way replay with original/erased fundamental phase.

Run HD before RL/RE, then use ``phase_controls.py`` to create the RE pressure
database.  The replay pair is a diagnostic counterfactual, not a fully coupled
physical pair.
"""

from __future__ import annotations

import argparse
import copy
import json
import math
from pathlib import Path
from typing import Any


SCRIPT_DIR = Path(__file__).resolve().parent
PERIOD = 1.3
DT = 1.0e-4
STEPS_PER_CYCLE = int(round(PERIOD / DT))
PRESSURE_INTERVAL = STEPS_PER_CYCLE // 100
VTK_INTERVAL = STEPS_PER_CYCLE // 20
EQUILIBRIUM_STEPS = 40_000
EQUILIBRIUM_DAMPING_FACTOR = 5.0
MC_HANDOFF_STEPS = 40_000
STABILITY_MAX_VELOCITY = 1.0e-3
STABILITY_MAX_DISPLACEMENT = 2.0e-2
CRITICAL_TIMESTEP_MODULUS = 23_800_000.0
POROSITY = 0.485
PERMEABILITY = 9.79e-12
LOW_LAG_SATURATION = 0.993
HIGH_LAG_SATURATION = 0.94
WATER_DENSITY = 1000.0
GRAVITY = 9.81
SEA_LEVEL = 1.0
SEABED_ELEVATION = 0.5
SUBMERGED_SURFACE_TRACTION = -WATER_DENSITY * GRAVITY * (
    SEA_LEVEL - SEABED_ELEVATION
)
SURFACE_TRACTION_PSET_ID = 4
PIPE_RADIUS = 0.06
PIPE_DIAMETER = 2.0 * PIPE_RADIUS
PIPE_COVER_RATIO = 0.25
PIPE_CENTER = [
    0.7,
    SEABED_ELEVATION - PIPE_RADIUS - PIPE_COVER_RATIO * PIPE_DIAMETER,
]
RELEASE_TIME = 3.0 * PERIOD
SANISAND_MC = 1.25
MC_FRICTION_DEG = math.degrees(math.asin(3.0 * SANISAND_MC / (6.0 + SANISAND_MC)))
SANISAND_G0 = 125.0
SANISAND_K0 = 150.0
SANISAND_PATM = 100_000.0
MC_STIFFNESS_REFERENCE_PRESSURE = 3_000.0
INITIAL_EFFECTIVE_K0 = 0.72
MINIMUM_NODAL_SUPPORT_FRACTION = 0.01

REGISTERED_DYNAMIC_STAGES = {
    "_LS_SANI": "LS",
    "_HS_SANI": "HS",
    "_HS_MC": "HM",
    "_HD_DRIVER": "HD",
    "_RL_SANI": "RL",
    "_RM_MC": "RM",
    "_RE_SANI": "RE",
}


def sanisand_elastic_moduli(
    mean_effective_stress: float,
    porosity: float = POROSITY,
) -> tuple[float, float]:
    """Return the implemented SANISAND tangent shear and bulk moduli."""
    if mean_effective_stress <= 0.0:
        raise ValueError("mean effective stress must be positive")
    if not 0.0 < porosity < 1.0:
        raise ValueError("porosity must be in (0, 1)")
    void_ratio = porosity / (1.0 - porosity)
    specific_volume = 1.0 + void_ratio
    pressure_ratio = mean_effective_stress / SANISAND_PATM
    shear_modulus = (
        math.sqrt(pressure_ratio)
        * SANISAND_G0
        * SANISAND_PATM
        * (2.97 - void_ratio) ** 2
        / specific_volume
    )
    bulk_modulus = (
        pressure_ratio ** (2.0 / 3.0)
        * SANISAND_K0
        * SANISAND_PATM
        * specific_volume
        / void_ratio
    )
    return shear_modulus, bulk_modulus


def isotropic_parameters_from_moduli(
    shear_modulus: float,
    bulk_modulus: float,
) -> tuple[float, float]:
    """Convert positive isotropic G and K to Young's modulus and Poisson ratio."""
    if shear_modulus <= 0.0 or bulk_modulus <= 0.0:
        raise ValueError("elastic moduli must be positive")
    denominator = 3.0 * bulk_modulus + shear_modulus
    youngs_modulus = 9.0 * bulk_modulus * shear_modulus / denominator
    poisson_ratio = (3.0 * bulk_modulus - 2.0 * shear_modulus) / (
        2.0 * denominator
    )
    return youngs_modulus, poisson_ratio


MC_MATCHED_SHEAR_MODULUS, MC_MATCHED_BULK_MODULUS = sanisand_elastic_moduli(
    MC_STIFFNESS_REFERENCE_PRESSURE
)
MC_MATCHED_YOUNGS_MODULUS, MC_MATCHED_POISSON_RATIO = isotropic_parameters_from_moduli(
    MC_MATCHED_SHEAR_MODULUS, MC_MATCHED_BULK_MODULUS
)


def mc_stiffness_calibration() -> dict[str, Any]:
    """Machine-readable provenance for the MC/SANISAND stiffness match."""
    void_ratio = POROSITY / (1.0 - POROSITY)
    return {
        "method": "SANISAND initial elastic tangent at reference mean effective stress",
        "reference_mean_effective_stress_pa": MC_STIFFNESS_REFERENCE_PRESSURE,
        "reference_porosity": POROSITY,
        "reference_void_ratio": void_ratio,
        "reference_specific_volume": 1.0 + void_ratio,
        "sanisand_G0": SANISAND_G0,
        "sanisand_K0": SANISAND_K0,
        "sanisand_Patm_pa": SANISAND_PATM,
        "matched_shear_modulus_pa": MC_MATCHED_SHEAR_MODULUS,
        "matched_bulk_modulus_pa": MC_MATCHED_BULK_MODULUS,
        "mc_youngs_modulus_pa": MC_MATCHED_YOUNGS_MODULUS,
        "mc_poisson_ratio": MC_MATCHED_POISSON_RATIO,
        "shear_formula": "G=(p'/Patm)^0.5*G0*Patm*(2.97-e)^2/(1+e)",
        "bulk_formula": "K=(p'/Patm)^(2/3)*K0*Patm*(1+e)/e",
    }


TIERS = {
    "screen": {"cycles": 8, "walltime": "12:00:00"},
    # 3 fixed cycles followed by 15 released cycles.
    "production": {"cycles": 18, "walltime": "36:00:00"},
}


def prefixed(mesh_directory: str, filename: str) -> str:
    return str(Path(mesh_directory) / filename) if mesh_directory != "." else filename


def mesh_block(mesh_directory: str, cell_size: float, particle_spacing: float) -> dict[str, Any]:
    bottom_y_1 = 0.5 * particle_spacing
    bottom_y_2 = 1.5 * particle_spacing
    hydrostatic = lambda y: WATER_DENSITY * GRAVITY * (SEA_LEVEL - y)
    return {
        "mesh": prefixed(mesh_directory, "mesh.txt"),
        "isoparametric": False,
        "entity_sets": prefixed(mesh_directory, "entity_sets.json"),
        "cell_type": "ED2Q4",
        "cellsize_min": cell_size,
        "io_type": "Ascii2D",
        "node_type": "N2D3PHASE&LAG",
        "particles_temperatures": prefixed(mesh_directory, "initial_temperature.txt"),
        "boundary_conditions": {
            "velocity_constraints": [
                {"nset_id": 0, "dir": direction, "velocity": 0.0}
                for direction in range(6)
            ]
            + [
                {"nset_id": side, "dir": direction, "velocity": 0.0}
                for side in (1, 3)
                for direction in (0, 2, 4)
            ],
            "particles_at_free_surface": [{"pset_id": 0, "nonfree_pset_id": 3}],
            "particle_pore_pressure_constraints": [
                {"pset_id": 1, "pore_pressure": hydrostatic(bottom_y_1)},
                {"pset_id": 2, "pore_pressure": hydrostatic(bottom_y_2)},
            ],
            "temperature_constraints": [{"nset_id": 2, "temperature": 0.0}],
        },
    }


def particle_block(mesh_directory: str) -> list[dict[str, Any]]:
    return [
        {
            "generator": {
                "check_duplicates": True,
                "location": prefixed(mesh_directory, "particles.txt"),
                "io_type": "Ascii2D",
                "particle_type": "P2D3PHASE&LAG",
                "material_id": 0,
                "liquid_material_id": 1,
                "type": "file",
            },
            "set_id": 0,
        }
    ]


def common_soil(material_type: str) -> dict[str, Any]:
    material: dict[str, Any] = {
        "id": 0,
        "name": "soil",
        "type": material_type,
        "density": 2650.0,
        # LinearElastic2D equilibrium uses the same registered reference
        # tangent as the MC handoff. The separate field below remains the
        # conservative wave-speed bound used by the timestep audit.
        "youngs_modulus": MC_MATCHED_YOUNGS_MODULUS,
        "poisson_ratio": MC_MATCHED_POISSON_RATIO,
        # Numerical wave-speed bound used by the explicit time-step check.
        # This is not a constitutive Young's modulus for SANISAND.
        "critical_timestep_modulus": CRITICAL_TIMESTEP_MODULUS,
        "eta": 10.0,
        "eta_T": 0.2,
        "eta_p": 1.0e-5,
        "porosity": POROSITY,
        "intrinsic_permeability": PERMEABILITY,
        "theta": 0.0,
        "thermal_expansivity": 1.0e-5,
        "thermal_conductivity": 100_000.0,
        "specific_heat": 2080.0,
        "density_ratio": 0.8,
        "initial_temperature": 0.0,
        "p_ref": 100_000.0,
    }
    if material_type != "LinearElastic2D":
        material["resume_state_policy"] = "reinitialize_from_particle"
    return material


def sanisand_material() -> dict[str, Any]:
    material = common_soil("SANISAND2D")
    # SANISAND computes its elastic tangent from G0/K0, p' and void ratio.
    # Remove inactive linear-elastic fields so the generated JSON cannot imply
    # that E/nu are part of the constitutive comparison.
    material.pop("youngs_modulus")
    material.pop("poisson_ratio")
    material.update(
        {
            "G0": SANISAND_G0,
            "K0": SANISAND_K0,
            "Mc": SANISAND_MC,
            "Lambda": 0.37,
            "N_c": 240.902666166616,
            "alpha_c": 3_370_000.0,
            "n_b": 1.25,
            "ch": 0.968,
            "n_d": 2.3,
            "h0": 12.0,
            "A0": 0.4,
            "Me": 0.89,
            "cz": 600.0,
            "zmax": 4.0,
            "m_iso": 0.01,
            "Patm": SANISAND_PATM,
            "P_min": 100.0,
            "FTOL": 1.0e-5,
            "STOL": 1.0e-5,
            "LTOL": 1.0e-6,
        }
    )
    return material


def mohr_coulomb_material() -> dict[str, Any]:
    material = common_soil("MohrCoulomb2D")
    material.update(
        {
            "youngs_modulus": MC_MATCHED_YOUNGS_MODULUS,
            "poisson_ratio": MC_MATCHED_POISSON_RATIO,
            "friction": MC_FRICTION_DEG,
            "cohesion": 0.0,
            "residual_friction": MC_FRICTION_DEG,
            "residual_cohesion": 0.0,
            "dilation": 0.0,
            "residual_dilation": 0.0,
            "softening": False,
            "softening_type": "none",
            "peak_pdstrain": 0.01,
            "residual_pdstrain": 0.05,
            "tension_cutoff": 0.0,
            "stiffness_calibration": mc_stiffness_calibration(),
        }
    )
    return material


def fluid_material(saturation: float, wave: bool, wave_height: float, particle_spacing: float) -> dict[str, Any]:
    if not 0.0 < saturation < 1.0:
        raise ValueError("A small gas phase is retained, so saturation must be in (0, 1)")
    return {
        "id": 1,
        "name": "water",
        "type": "Newtonian2D",
        "density": WATER_DENSITY,
        "liquid_saturation": saturation,
        "liquid_viscosity": 1.0e-3,
        "bulk_modulus": 2.2e9,
        "mu": 1.0e-3,
        "stab_para": 0.0,
        "liquid_thermal_conductivity": 100_000.0,
        "liquid_specific_heat": 4190.0,
        "liquid_expansivity": 7.0e-5,
        "liquid_compressibility": 0.4545e-9,
        "liquid_molar_mass": 0.018,
        "liquid_saturation_res": 0.001,
        "constant_mw": 2.0,
        "suction": 10.0,
        "two_phase": False,
        "name_3": "gas",
        "gas_thermal_conductivity": 0.0335,
        "gas_specific_heat": 2100.0,
        "gas_constant": 8.314,
        "gas_molar_mass": 0.029,
        "gas_saturation_res": 0.001,
        "gas_saturation": 1.0 - saturation,
        "gas_viscosity": 1.0e-5,
        "constant_mg": 3.0,
        "para_p0": 1000.0,
        "para_m": 0.5,
        "wave_pressure": wave,
        "domain_length_x": 1.5,
        "Nx": 1.5 / particle_spacing + 1.0,
        "depth_left": SEA_LEVEL - SEABED_ELEVATION,
        "depth_right": SEA_LEVEL - SEABED_ELEVATION,
        "sea_level": SEA_LEVEL,
        "wave_height_ini_": wave_height,
        "wave_period": PERIOD,
        "wave_ramp_time": PERIOD,
    }


def rigid_pipeline(*, fixed: bool, release_time: float, particle_spacing: float) -> dict[str, Any]:
    return {
        "enable": True,
        "fixed": fixed,
        "release_time": release_time,
        "allow_rotation": True,
        "center": PIPE_CENTER,
        "outer_radius": PIPE_RADIUS,
        "wall_thickness": PIPE_DIAMETER / 17.0,
        "density": 959.0,
        "contents_mass_per_length": 0.0,
        "contents_inertia_per_length": 0.0,
        "surface_markers": 128,
        "particle_radius": 0.5 * particle_spacing,
        "friction_coefficient": 0.15,
        "normal_penalty": 238_000.0,
        "normal_damping": 75.0,
        "tangential_damping": 75.0,
        "fluid_density": WATER_DENSITY,
        "translational_damping": 0.0,
        "rotational_damping": 0.0,
        "history_interval": PRESSURE_INTERVAL,
    }


def external_loading() -> dict[str, Any]:
    return {
        "gravity": [0.0, -GRAVITY],
        "particle_surface_traction": [
            {
                "pset_id": SURFACE_TRACTION_PSET_ID,
                "dir": 1,
                "traction": SUBMERGED_SURFACE_TRACTION,
                "facet": 2,
            }
        ],
    }


def resume_block(
    enable: bool,
    uuid: str,
    *,
    step: int = EQUILIBRIUM_STEPS,
    nsteps: int = EQUILIBRIUM_STEPS,
) -> dict[str, Any]:
    return {
        "resume": enable,
        "uuid": uuid,
        "step": step,
        "nsteps": nsteps,
        "start_from_this_step": False,
        "this_step": 0,
        "current_time": 0.0,
    }


def post_processing(
    result_path: str,
    material_type: str,
    *,
    initial: bool,
    write_checkpoint: bool = False,
) -> dict[str, Any]:
    if initial:
        solid = [
            "ids",
            "porosities",
            "volumes",
            "displacements",
            "velocities",
            "stresses",
        ]
        liquid = [
            "PIC_pore_pressures",
            "PIC_liquid_pressures",
            "PIC_gas_pressures",
            "PIC_ru",
        ]
        output_steps = EQUILIBRIUM_STEPS
    else:
        solid = [
            "ids",
            "porosities",
            "volumes",
            "displacements",
            "velocities",
            "stresses",
            "vertical_effective_stress_remaining_ratios",
        ]
        if material_type == "SANISAND2D":
            solid.extend(["eps_p_q", "void_ratio"])
        elif material_type == "MohrCoulomb2D":
            solid.extend(["pdstrain", "phi", "psi", "cohesion"])
        liquid = [
            "PIC_pore_pressures",
            # Checkpoint-relative mixture pressure used by the registered
            # excess-pressure curves and 2-D fields.
            "PIC_pore_pressure_excess",
            "PIC_liquid_pressures",
            "PIC_gas_pressures",
            "PIC_ru",
            "liquid_densities",
            "initial_vertical_effective_stresses",
            "dynamic_vertical_effective_stresses",
            "liquefaction_potentials",
            "gamma_sub",
            "momentary_liquefied",
            "liquid_pressure_gradients",
            "liquid_seepage_forces",
        ]
        output_steps = MC_HANDOFF_STEPS if write_checkpoint else VTK_INTERVAL
    return {
        "path": result_path,
        "write_hdf5": initial or write_checkpoint,
        "write_vtk": True,
        "vtk": solid,
        "liquid_vtk": liquid,
        "log_output_steps": False,
        "output_steps": output_steps,
    }


def base_config(mesh_directory: str, cell_size: float, particle_spacing: float) -> dict[str, Any]:
    return {
        "mesh": mesh_block(mesh_directory, cell_size, particle_spacing),
        "particles": particle_block(mesh_directory),
        "external_loading_conditions": external_loading(),
    }


def equilibrium_config(
    code: str,
    saturation: float,
    result_path: str,
    mesh_directory: str,
    cell_size: float,
    particle_spacing: float,
    wave_height: float,
) -> dict[str, Any]:
    config = base_config(mesh_directory, cell_size, particle_spacing)
    if math.isclose(saturation, LOW_LAG_SATURATION, abs_tol=1.0e-12):
        stress_state = "LS"
    elif math.isclose(saturation, HIGH_LAG_SATURATION, abs_tol=1.0e-12):
        stress_state = "HS"
    else:
        raise ValueError(
            "Equilibrium saturation has no registered initial effective-stress field"
        )
    config["mesh"]["particles_stresses"] = prefixed(
        mesh_directory, f"initial_effective_stresses_{stress_state}.txt"
    )
    config["mesh"]["particles_pore_pressures"] = {
        "file": prefixed(mesh_directory, "initial_liquid_pressures.txt")
    }
    uuid = f"PLP_{code}_EQ"
    config.update(
        {
            "title": f"Pipeline equilibrium {code}",
            "materials": [
                common_soil("LinearElastic2D"),
                fluid_material(saturation, False, wave_height, particle_spacing),
            ],
            "analysis": {
                "type": "ThermoMPMExplicitThreePhaseLag2D",
                "stress_update": "usf",
                # Use true pure PIC for equilibrium. The affine velocity field
                # amplifies the penalty-contact discontinuity around the pipe.
                "APIC": False,
                "PIC": 1.0,
                "PIC_T": 0.0,
                "free_surface": {
                    "free_surface_particle": "assign",
                    "volume_tolerance": 0.25,
                },
                "damping": {
                    "type": "Cundall",
                    "damping_factor": EQUILIBRIUM_DAMPING_FACTOR,
                },
                "pressure_smoothing": True,
                "pressure_smoothing_iterations": 1,
                "pressure_smoothing_in_loop": False,
                "minimum_nodal_support_fraction": MINIMUM_NODAL_SUPPORT_FRACTION,
                "rigid_pipeline": rigid_pipeline(
                    fixed=True,
                    release_time=(EQUILIBRIUM_STEPS + 1) * DT,
                    particle_spacing=particle_spacing,
                ),
                "uuid": uuid,
                "dt": DT,
                "nsteps": EQUILIBRIUM_STEPS,
                "stability_gate": True,
                "resume": resume_block(False, uuid),
            },
            "post_processing": post_processing(result_path, "LinearElastic2D", initial=True),
        }
    )
    return config


def mc_handoff_config(
    *,
    code: str,
    saturation: float,
    equilibrium_uuid: str,
    result_path: str,
    mesh_directory: str,
    cell_size: float,
    particle_spacing: float,
    wave_height: float,
) -> dict[str, Any]:
    """Create an explicit, damped MC relaxation from the elastic checkpoint."""

    config = base_config(mesh_directory, cell_size, particle_spacing)
    uuid = f"PLP_{code}_MC_RELAX"
    config.update(
        {
            "title": "High-lag Mohr-Coulomb handoff relaxation",
            "materials": [
                mohr_coulomb_material(),
                fluid_material(saturation, False, wave_height, particle_spacing),
            ],
            "analysis": {
                "type": "ThermoMPMExplicitThreePhaseLag2D",
                "stress_update": "usf",
                "APIC": False,
                "PIC": 1.0,
                "PIC_T": 0.0,
                "free_surface": {
                    "free_surface_particle": "assign",
                    "volume_tolerance": 0.25,
                },
                "damping": {
                    "type": "Cundall",
                    "damping_factor": EQUILIBRIUM_DAMPING_FACTOR,
                },
                "pressure_smoothing": True,
                "pressure_smoothing_iterations": 1,
                "pressure_smoothing_in_loop": False,
                "minimum_nodal_support_fraction": MINIMUM_NODAL_SUPPORT_FRACTION,
                "rigid_pipeline": rigid_pipeline(
                    fixed=True,
                    release_time=(MC_HANDOFF_STEPS + 1) * DT,
                    particle_spacing=particle_spacing,
                ),
                "uuid": uuid,
                "dt": DT,
                "nsteps": MC_HANDOFF_STEPS,
                "handoff_relaxation": True,
                "stability_gate": True,
                "resume": resume_block(True, equilibrium_uuid),
            },
            "post_processing": post_processing(
                result_path,
                "MohrCoulomb2D",
                initial=False,
                write_checkpoint=True,
            ),
        }
    )
    return config


def dynamic_config(
    *,
    code: str,
    title: str,
    saturation: float,
    material_type: str,
    equilibrium_uuid: str,
    result_path: str,
    nsteps: int,
    mesh_directory: str,
    cell_size: float,
    particle_spacing: float,
    wave_height: float,
    physical_wave: bool,
    equilibrium_steps: int = EQUILIBRIUM_STEPS,
    fixed_pipeline: bool = False,
    pressure_mode: str | None = None,
    pressure_path: str | None = None,
    pressure_prefix: str = "pressure",
) -> dict[str, Any]:
    config = base_config(mesh_directory, cell_size, particle_spacing)
    if material_type == "SANISAND2D":
        soil = sanisand_material()
    elif material_type == "MohrCoulomb2D":
        soil = mohr_coulomb_material()
    else:
        raise ValueError(material_type)
    analysis: dict[str, Any] = {
        "type": "ThermoMPMExplicitThreePhaseLag2D",
        "stress_update": "usf",
        "APIC": True,
        "PIC": 0.0,
        "PIC_T": 0.0,
        "free_surface": {
            "free_surface_particle": "assign",
            "volume_tolerance": 0.25,
        },
        "damping": {"type": "Cundall", "damping_factor": 0.0},
        "pressure_smoothing": True,
        "pressure_smoothing_iterations": 1,
        "pressure_smoothing_in_loop": False,
        "minimum_nodal_support_fraction": MINIMUM_NODAL_SUPPORT_FRACTION,
        "rigid_pipeline": rigid_pipeline(
            fixed=fixed_pipeline,
            release_time=(nsteps + 1) * DT if fixed_pipeline else RELEASE_TIME,
            particle_spacing=particle_spacing,
        ),
        "uuid": f"PLP_{code}",
        "dt": DT,
        "nsteps": nsteps,
        "resume": resume_block(
            True,
            equilibrium_uuid,
            step=equilibrium_steps,
            nsteps=equilibrium_steps,
        ),
    }
    if pressure_mode:
        if pressure_path is None:
            raise ValueError("pressure_path is required for database mode")
        analysis["prescribed_phase_pressures"] = {
            # Solver semantics are mutually exclusive: enable means read;
            # write records the coupled solution while normal pressure update
            # remains active.
            "enable": pressure_mode == "read",
            "write": pressure_mode == "write",
            "path": pressure_path,
            "file_prefix": pressure_prefix,
            "source_dt": DT,
            "step_interval": PRESSURE_INTERVAL,
            "max_step": nsteps,
            "mapping": "id",
        }
    config.update(
        {
            "title": title,
            "materials": [
                soil,
                fluid_material(saturation, physical_wave, wave_height, particle_spacing),
            ],
            "analysis": analysis,
            "post_processing": post_processing(result_path, material_type, initial=False),
        }
    )
    return config


def registered_tier(uuid: str) -> str | None:
    """Return the registered study tier encoded in a generated UUID."""

    for tier in TIERS:
        if uuid.startswith(f"PLP_{tier.upper()}_"):
            return tier
    return None


def registered_dynamic_stage(uuid: str) -> str | None:
    """Return the registered dynamic case encoded in a generated UUID."""

    for suffix, stage in REGISTERED_DYNAMIC_STAGES.items():
        if uuid.endswith(suffix):
            return stage
    return None


def registered_dynamic_prefix(uuid: str) -> str | None:
    """Return the UUID prefix shared by a dynamic case and its prerequisites."""

    for suffix in REGISTERED_DYNAMIC_STAGES:
        if uuid.endswith(suffix):
            return uuid[: -len(suffix)]
    return None


def require_registered_saturation(
    analysis: dict[str, Any], fluid: dict[str, Any], expected: float
) -> None:
    saturation = float(fluid["liquid_saturation"])
    if not math.isclose(saturation, expected, rel_tol=0.0, abs_tol=1.0e-12):
        raise ValueError(
            f"{analysis['uuid']}: liquid saturation must be {expected:g}, "
            f"not {saturation:g}"
        )


def validate_pressure_registration(
    analysis: dict[str, Any], *, tier: str, stage: str
) -> None:
    """Validate mutually exclusive physical-wave and replay database roles."""

    pressure = analysis.get("prescribed_phase_pressures")
    if stage in {"LS", "HS", "HM"}:
        if pressure is not None:
            raise ValueError(
                f"{analysis['uuid']}: fully coupled {stage} cannot use a pressure database"
            )
        return

    if not isinstance(pressure, dict):
        raise ValueError(f"{analysis['uuid']}: {stage} pressure database is absent")
    expected_read = stage in {"RL", "RM", "RE"}
    if bool(pressure.get("enable")) != expected_read:
        raise ValueError(f"{analysis['uuid']}: {stage} pressure read mode is wrong")
    if bool(pressure.get("write")) != (stage == "HD"):
        raise ValueError(f"{analysis['uuid']}: {stage} pressure write mode is wrong")
    expected_prefix = "phase_erased" if stage == "RE" else "pressure"
    expected_leaf = "phase_erased" if stage == "RE" else "lagged"
    path = str(pressure.get("path", ""))
    if (
        pressure.get("file_prefix") != expected_prefix
        or not path.startswith(f"pressure_databases/{tier}")
        or not path.endswith(f"/{expected_leaf}")
        or not math.isclose(float(pressure.get("source_dt", -1.0)), DT, abs_tol=1.0e-15)
        or int(pressure.get("step_interval", -1)) != PRESSURE_INTERVAL
        or int(pressure.get("max_step", -1)) != int(analysis["nsteps"])
        or pressure.get("mapping") != "id"
    ):
        raise ValueError(
            f"{analysis['uuid']}: {stage} pressure database registration is stale"
        )


def validate_config(config: dict[str, Any]) -> None:
    soil, fluid = config["materials"]
    analysis = config["analysis"]
    pipeline = analysis["rigid_pipeline"]
    resume_enabled = bool(analysis["resume"].get("resume"))
    handoff_relaxation = bool(analysis.get("handoff_relaxation", False))
    uuid = str(analysis["uuid"])
    if analysis.get("type") != "ThermoMPMExplicitThreePhaseLag2D":
        raise ValueError(f"{uuid}: unexpected analysis type")
    if not math.isclose(float(analysis.get("dt", -1.0)), DT, abs_tol=1.0e-15):
        raise ValueError(f"{uuid}: dt must remain {DT:g} s")
    if pipeline.get("enable") is not True:
        raise ValueError(f"{uuid}: rigid pipeline must be enabled")
    support_fraction = float(analysis.get("minimum_nodal_support_fraction", -1.0))
    if not math.isclose(
        support_fraction, MINIMUM_NODAL_SUPPORT_FRACTION, abs_tol=1.0e-12
    ):
        raise ValueError(
            f"{analysis['uuid']}: minimum_nodal_support_fraction must be "
            f"{MINIMUM_NODAL_SUPPORT_FRACTION:g}"
        )
    if "/APIC" in analysis or not isinstance(analysis.get("APIC"), bool):
        raise ValueError(f"{analysis['uuid']}: APIC key is absent or malformed")
    if not resume_enabled:
        if handoff_relaxation:
            raise ValueError(
                f"{analysis['uuid']}: MC handoff relaxation must resume the elastic checkpoint"
            )
        if (
            soil["type"] != "LinearElastic2D"
            or analysis["APIC"]
            or float(analysis["PIC"]) != 1.0
            or float(analysis.get("PIC_T", -1.0)) != 0.0
            or not pipeline.get("fixed")
            or bool(fluid.get("wave_pressure"))
        ):
            raise ValueError(
                f"{analysis['uuid']}: equilibrium must use LinearElastic2D, "
                "APIC=false, PIC=1, a fixed pipe and no wave"
            )
        if int(analysis["nsteps"]) != EQUILIBRIUM_STEPS:
            raise ValueError(
                f"{analysis['uuid']}: equilibrium must run {EQUILIBRIUM_STEPS} steps"
            )
        damping_factor = float(analysis["damping"].get("damping_factor", 0.0))
        if not math.isclose(
            damping_factor, EQUILIBRIUM_DAMPING_FACTOR, abs_tol=1.0e-12
        ) or analysis["damping"].get("type") != "Cundall":
            raise ValueError(
                f"{analysis['uuid']}: equilibrium damping must be "
                f"{EQUILIBRIUM_DAMPING_FACTOR:g} 1/s"
            )
        initial_stress = config["mesh"].get("particles_stresses", "")
        initial_pressure = config["mesh"].get(
            "particles_pore_pressures", {}
        ).get("file", "")
        if not initial_stress or not initial_pressure:
            raise ValueError(
                f"{analysis['uuid']}: equilibrium initial stress/pressure field is absent"
            )
        if analysis.get("stability_gate") is not True:
            raise ValueError(f"{analysis['uuid']}: equilibrium stability gate is disabled")
        saturation = float(fluid["liquid_saturation"])
        if math.isclose(saturation, LOW_LAG_SATURATION, rel_tol=0.0, abs_tol=1.0e-12):
            expected_state = "LS"
        elif math.isclose(
            saturation, HIGH_LAG_SATURATION, rel_tol=0.0, abs_tol=1.0e-12
        ):
            expected_state = "HS"
        else:
            raise ValueError(
                f"{analysis['uuid']}: equilibrium saturation is not registered as LS/HS"
            )
        if registered_tier(uuid) is None or not uuid.endswith(f"_{expected_state}_EQ"):
            raise ValueError(f"{uuid}: equilibrium tier/case registration is absent")
        if not str(initial_stress).endswith(
            f"initial_effective_stresses_{expected_state}.txt"
        ):
            raise ValueError(
                f"{analysis['uuid']}: equilibrium saturation/stress field mismatch"
            )
        expected_release = (EQUILIBRIUM_STEPS + 1) * DT
        if not math.isclose(
            float(pipeline.get("release_time", -1.0)),
            expected_release,
            abs_tol=1.0e-12,
        ):
            raise ValueError(f"{analysis['uuid']}: equilibrium pipe release is stale")
        if not config["post_processing"].get("write_hdf5"):
            raise ValueError(f"{analysis['uuid']}: equilibrium checkpoint is disabled")
    elif handoff_relaxation:
        damping_factor = float(analysis["damping"].get("damping_factor", 0.0))
        if (
            soil["type"] != "MohrCoulomb2D"
            or analysis["APIC"]
            or float(analysis["PIC"]) != 1.0
            or float(analysis.get("PIC_T", -1.0)) != 0.0
            or not pipeline.get("fixed")
            or bool(fluid.get("wave_pressure"))
            or int(analysis["nsteps"]) != MC_HANDOFF_STEPS
            or not math.isclose(
                damping_factor,
                EQUILIBRIUM_DAMPING_FACTOR,
                abs_tol=1.0e-12,
            )
            or analysis["damping"].get("type") != "Cundall"
        ):
            raise ValueError(
                f"{analysis['uuid']}: MC handoff must use MohrCoulomb2D, "
                "APIC=false, PIC=1, damping=5 1/s, fixed pipe and no wave"
            )
        if analysis.get("stability_gate") is not True:
            raise ValueError(f"{analysis['uuid']}: MC handoff stability gate is disabled")
        if not config["post_processing"].get("write_hdf5"):
            raise ValueError(f"{analysis['uuid']}: MC handoff checkpoint is disabled")
        if analysis.get("prescribed_phase_pressures") is not None:
            raise ValueError(
                f"{analysis['uuid']}: MC handoff cannot read or write wave pressures"
            )
        require_registered_saturation(analysis, fluid, HIGH_LAG_SATURATION)
        resume = analysis["resume"]
        expected_equilibrium_uuid = uuid.removesuffix("_MC_RELAX") + "_EQ"
        if (
            registered_tier(uuid) is None
            or not uuid.endswith("_HS_MC_RELAX")
            or str(resume.get("uuid", "")) != expected_equilibrium_uuid
            or int(resume.get("step", -1)) != EQUILIBRIUM_STEPS
            or int(resume.get("nsteps", -1)) != EQUILIBRIUM_STEPS
            or not math.isclose(
                float(pipeline.get("release_time", -1.0)),
                (MC_HANDOFF_STEPS + 1) * DT,
                abs_tol=1.0e-12,
            )
        ):
            raise ValueError(f"{analysis['uuid']}: MC handoff registration is stale")
    else:
        tier = registered_tier(uuid)
        stage = registered_dynamic_stage(uuid)
        dynamic_prefix = registered_dynamic_prefix(uuid)
        if tier is None or stage is None or dynamic_prefix is None:
            raise ValueError(f"{uuid}: dynamic tier/case registration is absent")
        expected_nsteps = int(TIERS[tier]["cycles"] * STEPS_PER_CYCLE)
        damping_factor = float(analysis["damping"].get("damping_factor", -1.0))
        if (
            not analysis["APIC"]
            or float(analysis["PIC"]) != 0.0
            or float(analysis.get("PIC_T", -1.0)) != 0.0
            or damping_factor != 0.0
            or analysis["damping"].get("type") != "Cundall"
            or int(analysis["nsteps"]) != expected_nsteps
        ):
            raise ValueError(
                f"{uuid}: dynamic {tier}/{stage} must use APIC=true, PIC=0, "
                f"zero damping and {expected_nsteps} steps"
            )
        if analysis.get("stability_gate"):
            raise ValueError(f"{uuid}: dynamic phase cannot publish a static gate")

        expected_material = (
            "MohrCoulomb2D" if stage in {"HM", "RM"} else "SANISAND2D"
        )
        expected_saturation = (
            LOW_LAG_SATURATION if stage == "LS" else HIGH_LAG_SATURATION
        )
        if soil["type"] != expected_material:
            raise ValueError(f"{uuid}: {stage} material registration is wrong")
        require_registered_saturation(analysis, fluid, expected_saturation)

        physical_wave = stage in {"LS", "HS", "HM", "HD"}
        if bool(fluid.get("wave_pressure")) != physical_wave:
            raise ValueError(f"{uuid}: {stage} physical-wave registration is wrong")
        expected_fixed = stage == "HD"
        expected_release = (expected_nsteps + 1) * DT if expected_fixed else RELEASE_TIME
        if bool(pipeline.get("fixed")) != expected_fixed or not math.isclose(
            float(pipeline.get("release_time", -1.0)), expected_release, abs_tol=1.0e-12
        ):
            raise ValueError(f"{uuid}: {stage} pipe release registration is wrong")
        validate_pressure_registration(analysis, tier=tier, stage=stage)

        resume = analysis["resume"]
        expected_resume_uuid = (
            f"{dynamic_prefix}_HS_MC_RELAX"
            if stage in {"HM", "RM"}
            else f"{dynamic_prefix}_LS_EQ"
            if stage == "LS"
            else f"{dynamic_prefix}_HS_EQ"
        )
        if (
            not resume_enabled
            or str(resume.get("uuid", "")) != expected_resume_uuid
            or int(resume.get("step", -1)) != EQUILIBRIUM_STEPS
            or int(resume.get("nsteps", -1)) != EQUILIBRIUM_STEPS
        ):
            raise ValueError(f"{uuid}: {stage} resume registration is wrong")
    if abs(float(soil["intrinsic_permeability"]) - PERMEABILITY) > 1.0e-20:
        raise ValueError(f"{analysis['uuid']}: permeability changed")
    if not math.isclose(
        float(soil["critical_timestep_modulus"]),
        CRITICAL_TIMESTEP_MODULUS,
        abs_tol=1.0e-12,
    ):
        raise ValueError(f"{analysis['uuid']}: critical timestep modulus changed")
    saturation_sum = float(fluid["liquid_saturation"]) + float(fluid["gas_saturation"])
    if not math.isclose(saturation_sum, 1.0, abs_tol=1.0e-12):
        raise ValueError(f"{analysis['uuid']}: phase saturations do not sum to one")
    if not all(
        math.isclose(float(a), float(b), abs_tol=1.0e-12)
        for a, b in zip(pipeline["center"], PIPE_CENTER)
    ):
        raise ValueError(f"{analysis['uuid']}: pipe center is inconsistent")
    nsteps = int(analysis["nsteps"])
    if float(analysis["nsteps"]) != nsteps:
        raise ValueError(f"{analysis['uuid']}: nsteps is not integral")
    pressure = analysis.get("prescribed_phase_pressures")
    if pressure and (nsteps % int(pressure["step_interval"]) != 0):
        raise ValueError(f"{analysis['uuid']}: pressure database ends mid-frame")
    if pressure and pressure.get("enable") and pressure.get("write"):
        raise ValueError(f"{analysis['uuid']}: pressure read/write modes conflict")
    requested = set(config["post_processing"]["vtk"])
    if analysis["resume"].get("resume"):
        required_solid = {
            "ids",
            "volumes",
            "stresses",
            "vertical_effective_stress_remaining_ratios",
        }
        required_liquid = {
            "PIC_pore_pressure_excess",
            "liquid_densities",
            "gamma_sub",
            "liquid_seepage_forces",
        }
        missing_solid = required_solid - requested
        missing_liquid = required_liquid - set(
            config["post_processing"]["liquid_vtk"]
        )
        if missing_solid or missing_liquid:
            raise ValueError(
                f"{analysis['uuid']}: missing primary liquefaction outputs "
                f"solid={sorted(missing_solid)}, liquid={sorted(missing_liquid)}"
            )
    if soil["type"] == "SANISAND2D" and "pdstrain" in requested:
        raise ValueError("SANISAND config requests MC-only state")
    if soil["type"] == "SANISAND2D" and (
        "youngs_modulus" in soil or "poisson_ratio" in soil
    ):
        raise ValueError("SANISAND config contains inactive E/nu parameters")
    if soil["type"] == "MohrCoulomb2D" and "eps_p_q" in requested:
        raise ValueError("MC config requests SANISAND-only state")
    if soil["type"] == "MohrCoulomb2D":
        shear = float(soil["youngs_modulus"]) / (2.0 * (1.0 + float(soil["poisson_ratio"])))
        bulk = float(soil["youngs_modulus"]) / (
            3.0 * (1.0 - 2.0 * float(soil["poisson_ratio"]))
        )
        if not math.isclose(shear, MC_MATCHED_SHEAR_MODULUS, rel_tol=1.0e-12):
            raise ValueError("MC shear stiffness is not matched to SANISAND")
        if not math.isclose(bulk, MC_MATCHED_BULK_MODULUS, rel_tol=1.0e-12):
            raise ValueError("MC bulk stiffness is not matched to SANISAND")
        if soil.get("stiffness_calibration") != mc_stiffness_calibration():
            raise ValueError("MC stiffness calibration metadata is absent or stale")


def generate_tier(
    tier: str,
    output_root: Path,
    *,
    mesh_directory: str,
    cell_size: float,
    particle_spacing: float,
    wave_height: float,
    label: str,
) -> Path:
    settings = TIERS[tier]
    nsteps = int(settings["cycles"] * STEPS_PER_CYCLE)
    suffix = "" if label == "baseline" else f"_{label}"
    directory = output_root / f"{tier}{suffix}"
    directory.mkdir(parents=True, exist_ok=True)
    result_path = f"results/{tier}{suffix}/"
    pressure_root = f"pressure_databases/{tier}{suffix}"
    label_code = "" if label == "baseline" else f"_{label.upper()}"

    low_eq_code = f"{tier.upper()}{label_code}_LS"
    high_eq_code = f"{tier.upper()}{label_code}_HS"
    low_eq_uuid = f"PLP_{low_eq_code}_EQ"
    high_eq_uuid = f"PLP_{high_eq_code}_EQ"
    mc_handoff_code = f"{tier.upper()}{label_code}_HS"
    mc_handoff_uuid = f"PLP_{mc_handoff_code}_MC_RELAX"
    configs: dict[str, dict[str, Any]] = {
        "EQ_LS": equilibrium_config(
            low_eq_code,
            LOW_LAG_SATURATION,
            result_path,
            mesh_directory,
            cell_size,
            particle_spacing,
            wave_height,
        ),
        "EQ_HS": equilibrium_config(
            high_eq_code,
            HIGH_LAG_SATURATION,
            result_path,
            mesh_directory,
            cell_size,
            particle_spacing,
            wave_height,
        ),
        "MC_EQ": mc_handoff_config(
            code=mc_handoff_code,
            saturation=HIGH_LAG_SATURATION,
            equilibrium_uuid=high_eq_uuid,
            result_path=result_path,
            mesh_directory=mesh_directory,
            cell_size=cell_size,
            particle_spacing=particle_spacing,
            wave_height=wave_height,
        ),
    }
    dynamic_common = {
        "result_path": result_path,
        "nsteps": nsteps,
        "mesh_directory": mesh_directory,
        "cell_size": cell_size,
        "particle_spacing": particle_spacing,
        "wave_height": wave_height,
    }
    configs.update(
        {
            "LS": dynamic_config(
                code=f"{tier.upper()}{label_code}_LS_SANI",
                title="LS physical low-lag reference: near saturated SANISAND",
                saturation=LOW_LAG_SATURATION,
                material_type="SANISAND2D",
                equilibrium_uuid=low_eq_uuid,
                physical_wave=True,
                **dynamic_common,
            ),
            "HS": dynamic_config(
                code=f"{tier.upper()}{label_code}_HS_SANI",
                title="HS physical high-lag main case: S=0.94 SANISAND",
                saturation=HIGH_LAG_SATURATION,
                material_type="SANISAND2D",
                equilibrium_uuid=high_eq_uuid,
                physical_wave=True,
                **dynamic_common,
            ),
            "HM": dynamic_config(
                code=f"{tier.upper()}{label_code}_HS_MC",
                title="HM fully coupled high-lag Mohr-Coulomb context",
                saturation=HIGH_LAG_SATURATION,
                material_type="MohrCoulomb2D",
                equilibrium_uuid=mc_handoff_uuid,
                equilibrium_steps=MC_HANDOFF_STEPS,
                physical_wave=True,
                **dynamic_common,
            ),
            "HD": dynamic_config(
                code=f"{tier.upper()}{label_code}_HD_DRIVER",
                title="HD fixed-pipe high-lag pressure database driver",
                saturation=HIGH_LAG_SATURATION,
                material_type="SANISAND2D",
                equilibrium_uuid=high_eq_uuid,
                physical_wave=True,
                fixed_pipeline=True,
                pressure_mode="write",
                pressure_path=f"{pressure_root}/lagged",
                pressure_prefix="pressure",
                **dynamic_common,
            ),
            "RL": dynamic_config(
                code=f"{tier.upper()}{label_code}_RL_SANI",
                title="RL one-way replay control: original lag, SANISAND",
                saturation=HIGH_LAG_SATURATION,
                material_type="SANISAND2D",
                equilibrium_uuid=high_eq_uuid,
                physical_wave=False,
                pressure_mode="read",
                pressure_path=f"{pressure_root}/lagged",
                pressure_prefix="pressure",
                **dynamic_common,
            ),
            "RM": dynamic_config(
                code=f"{tier.upper()}{label_code}_RM_MC",
                title=(
                    "RM matched-pressure constitutive ablation: original lag, "
                    "Mohr-Coulomb"
                ),
                saturation=HIGH_LAG_SATURATION,
                material_type="MohrCoulomb2D",
                equilibrium_uuid=mc_handoff_uuid,
                equilibrium_steps=MC_HANDOFF_STEPS,
                physical_wave=False,
                pressure_mode="read",
                pressure_path=f"{pressure_root}/lagged",
                pressure_prefix="pressure",
                **dynamic_common,
            ),
            "RE": dynamic_config(
                code=f"{tier.upper()}{label_code}_RE_SANI",
                title="RE one-way counterfactual: fundamental phase erased",
                saturation=HIGH_LAG_SATURATION,
                material_type="SANISAND2D",
                equilibrium_uuid=high_eq_uuid,
                physical_wave=False,
                pressure_mode="read",
                pressure_path=f"{pressure_root}/phase_erased",
                pressure_prefix="phase_erased",
                **dynamic_common,
            ),
        }
    )

    rows = [
        ("EQ_LS", 1, "physical equilibrium", "none"),
        ("EQ_HS", 1, "physical equilibrium", "none"),
        ("MC_EQ", 2, "MC handoff relaxation", "EQ_HS"),
        ("LS", 2, "physical low-lag reference", "EQ_LS"),
        ("HS", 2, "physical high-lag main case", "EQ_HS"),
        ("HM", 3, "fully-coupled MC context", "MC_EQ"),
        ("HD", 2, "fixed-pipe replay driver", "EQ_HS"),
        ("RL", 4, "one-way lagged replay", "HD"),
        ("RM", 4, "matched-pressure constitutive ablation", "HD + MC_EQ"),
        ("RE", 4, "one-way phase-erased counterfactual", "HD + phase_controls.py"),
    ]
    manifest: dict[str, Any] = {
        "schema": "pipeline-phase-lag-study-v1",
        "tier": tier,
        "label": label,
        "scientific_scope": (
            "MC_EQ is numerical handoff relaxation and not a comparison case; "
            "LS/HS/HM are fully coupled physical cases; RL/RM/RE are one-way "
            "diagnostic replays and must not be described as physical solutions"
        ),
        "constants": {
            "wave_height_m": wave_height,
            "wave_period_s": PERIOD,
            "dt_s": DT,
            "cycles": settings["cycles"],
            "release_time_s": RELEASE_TIME,
            "post_release_cycles": settings["cycles"] - 3,
            "equilibrium_steps": EQUILIBRIUM_STEPS,
            "equilibrium_duration_s": EQUILIBRIUM_STEPS * DT,
            "equilibrium_pic": 1.0,
            "equilibrium_damping_factor_per_s": EQUILIBRIUM_DAMPING_FACTOR,
            "mc_handoff_steps": MC_HANDOFF_STEPS,
            "mc_handoff_duration_s": MC_HANDOFF_STEPS * DT,
            "stability_max_velocity_m_s": STABILITY_MAX_VELOCITY,
            "stability_max_displacement_m": STABILITY_MAX_DISPLACEMENT,
            "porosity": POROSITY,
            "intrinsic_permeability_m2": PERMEABILITY,
            "hydraulic_conductivity_approx_m_s": (
                PERMEABILITY * WATER_DENSITY * GRAVITY / 1.0e-3
            ),
            "pipe_diameter_m": PIPE_DIAMETER,
            "cover_over_diameter": PIPE_COVER_RATIO,
            "pipe_invert_depth_m": (
                SEABED_ELEVATION - (PIPE_CENTER[1] - PIPE_RADIUS)
            ),
            "submerged_surface_traction_pa": SUBMERGED_SURFACE_TRACTION,
            "mc_friction_deg_from_Mc": MC_FRICTION_DEG,
            "mc_stiffness_calibration": mc_stiffness_calibration(),
            "background_cell_size_m": cell_size,
            "particle_spacing_m": particle_spacing,
        },
        "cases": [],
    }
    for code, order, role, dependency in rows:
        config = configs[code]
        validate_config(config)
        filename = f"{order:02d}_{code}.json"
        (directory / filename).write_text(
            json.dumps(config, indent=2) + "\n", encoding="utf-8"
        )
        manifest["cases"].append(
            {
                "code": code,
                "run_order": order,
                "role": role,
                "dependency": dependency,
                "config": filename,
                "uuid": config["analysis"]["uuid"],
            }
        )
    manifest_path = directory / "manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    return manifest_path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tier", choices=("screen", "production", "both"), default="both")
    parser.add_argument("--output-root", type=Path, default=Path("configs"))
    parser.add_argument("--mesh-dir", default=".")
    parser.add_argument("--cell-size", type=float, default=0.02)
    parser.add_argument("--particle-spacing", type=float, default=0.01)
    parser.add_argument("--wave-height", type=float, default=0.12)
    parser.add_argument("--label", default="baseline")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    output_root = args.output_root
    if not output_root.is_absolute():
        output_root = SCRIPT_DIR / output_root
    tiers = TIERS if args.tier == "both" else (args.tier,)
    manifests = []
    for tier in tiers:
        manifests.append(
            generate_tier(
                tier,
                output_root,
                mesh_directory=args.mesh_dir,
                cell_size=args.cell_size,
                particle_spacing=args.particle_spacing,
                wave_height=args.wave_height,
                label=args.label,
            )
        )
    for manifest in manifests:
        print(f"Generated and validated {manifest}")
    print(
        "Run screen cases first. Do not claim MPM necessity unless the reported "
        "pipeline motion reaches |uy|/D >= 0.5 or shows documented flow/contact "
        "topology change."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
