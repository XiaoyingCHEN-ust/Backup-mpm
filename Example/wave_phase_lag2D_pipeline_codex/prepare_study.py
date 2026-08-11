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
EQUILIBRIUM_STEPS = 10_000
POROSITY = 0.485
PERMEABILITY = 9.79e-12
PIPE_RADIUS = 0.06
PIPE_DIAMETER = 2.0 * PIPE_RADIUS
PIPE_COVER_RATIO = 0.25
PIPE_CENTER = [0.7, 0.5 - PIPE_RADIUS - PIPE_COVER_RATIO * PIPE_DIAMETER]
RELEASE_TIME = 3.0 * PERIOD
SANISAND_MC = 1.25
MC_FRICTION_DEG = math.degrees(math.asin(3.0 * SANISAND_MC / (6.0 + SANISAND_MC)))
SANISAND_G0 = 125.0
SANISAND_K0 = 150.0
SANISAND_PATM = 100_000.0
MC_STIFFNESS_REFERENCE_PRESSURE = 3_000.0


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
    hydrostatic = lambda y: 1000.0 * 9.81 * (1.0 - y)
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
        "youngs_modulus": 23_800_000.0,
        "poisson_ratio": 0.4,
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
        "density": 1000.0,
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
        "depth_left": 0.5,
        "depth_right": 0.5,
        "sea_level": 1.0,
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
        "fluid_density": 1000.0,
        "translational_damping": 0.0,
        "rotational_damping": 0.0,
        "history_interval": PRESSURE_INTERVAL,
    }


def external_loading() -> dict[str, Any]:
    return {
        "gravity": [0.0, -9.81],
        "particle_surface_traction": [
            {"pset_id": 0, "dir": 1, "traction": 0.0, "facet": 2}
        ],
    }


def resume_block(enable: bool, uuid: str) -> dict[str, Any]:
    return {
        "resume": enable,
        "uuid": uuid,
        "step": EQUILIBRIUM_STEPS,
        "nsteps": EQUILIBRIUM_STEPS,
        "start_from_this_step": False,
        "this_step": 0,
        "current_time": 0.0,
    }


def post_processing(result_path: str, material_type: str, *, initial: bool) -> dict[str, Any]:
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
        output_steps = VTK_INTERVAL
    return {
        "path": result_path,
        "write_hdf5": initial,
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
                "APIC": True,
                "PIC": 1.0,
                "PIC_T": 0.0,
                "free_surface": {
                    "free_surface_particle": "assign",
                    "volume_tolerance": 0.25,
                },
                "damping": {"type": "Cundall", "damping_factor": 0.05},
                "pressure_smoothing": True,
                "pressure_smoothing_iterations": 1,
                "pressure_smoothing_in_loop": False,
                "rigid_pipeline": rigid_pipeline(
                    fixed=True,
                    release_time=(EQUILIBRIUM_STEPS + 1) * DT,
                    particle_spacing=particle_spacing,
                ),
                "uuid": uuid,
                "dt": DT,
                "nsteps": EQUILIBRIUM_STEPS,
                "resume": resume_block(False, uuid),
            },
            "post_processing": post_processing(result_path, "LinearElastic2D", initial=True),
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
        "rigid_pipeline": rigid_pipeline(
            fixed=fixed_pipeline,
            release_time=(nsteps + 1) * DT if fixed_pipeline else RELEASE_TIME,
            particle_spacing=particle_spacing,
        ),
        "uuid": f"PLP_{code}",
        "dt": DT,
        "nsteps": nsteps,
        "resume": resume_block(True, equilibrium_uuid),
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


def validate_config(config: dict[str, Any]) -> None:
    soil, fluid = config["materials"]
    analysis = config["analysis"]
    pipeline = analysis["rigid_pipeline"]
    if "/APIC" in analysis or analysis.get("APIC") is not True:
        raise ValueError(f"{analysis['uuid']}: APIC key is absent or malformed")
    if abs(float(soil["intrinsic_permeability"]) - PERMEABILITY) > 1.0e-20:
        raise ValueError(f"{analysis['uuid']}: permeability changed")
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
    configs: dict[str, dict[str, Any]] = {
        "EQ_LS": equilibrium_config(
            low_eq_code,
            0.999,
            result_path,
            mesh_directory,
            cell_size,
            particle_spacing,
            wave_height,
        ),
        "EQ_HS": equilibrium_config(
            high_eq_code,
            0.94,
            result_path,
            mesh_directory,
            cell_size,
            particle_spacing,
            wave_height,
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
                saturation=0.999,
                material_type="SANISAND2D",
                equilibrium_uuid=low_eq_uuid,
                physical_wave=True,
                **dynamic_common,
            ),
            "HS": dynamic_config(
                code=f"{tier.upper()}{label_code}_HS_SANI",
                title="HS physical high-lag main case: S=0.94 SANISAND",
                saturation=0.94,
                material_type="SANISAND2D",
                equilibrium_uuid=high_eq_uuid,
                physical_wave=True,
                **dynamic_common,
            ),
            "HM": dynamic_config(
                code=f"{tier.upper()}{label_code}_HS_MC",
                title="HM fully coupled high-lag Mohr-Coulomb context",
                saturation=0.94,
                material_type="MohrCoulomb2D",
                equilibrium_uuid=high_eq_uuid,
                physical_wave=True,
                **dynamic_common,
            ),
            "HD": dynamic_config(
                code=f"{tier.upper()}{label_code}_HD_DRIVER",
                title="HD fixed-pipe high-lag pressure database driver",
                saturation=0.94,
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
                saturation=0.94,
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
                saturation=0.94,
                material_type="MohrCoulomb2D",
                equilibrium_uuid=high_eq_uuid,
                physical_wave=False,
                pressure_mode="read",
                pressure_path=f"{pressure_root}/lagged",
                pressure_prefix="pressure",
                **dynamic_common,
            ),
            "RE": dynamic_config(
                code=f"{tier.upper()}{label_code}_RE_SANI",
                title="RE one-way counterfactual: fundamental phase erased",
                saturation=0.94,
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
        ("LS", 2, "physical low-lag reference", "EQ_LS"),
        ("HS", 2, "physical high-lag main case", "EQ_HS"),
        ("HM", 2, "fully-coupled MC context", "EQ_HS"),
        ("HD", 2, "fixed-pipe replay driver", "EQ_HS"),
        ("RL", 4, "one-way lagged replay", "HD"),
        ("RM", 4, "matched-pressure constitutive ablation", "HD"),
        ("RE", 4, "one-way phase-erased counterfactual", "HD + phase_controls.py"),
    ]
    manifest: dict[str, Any] = {
        "schema": "pipeline-phase-lag-study-v1",
        "tier": tier,
        "label": label,
        "scientific_scope": (
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
            "porosity": POROSITY,
            "intrinsic_permeability_m2": PERMEABILITY,
            "hydraulic_conductivity_approx_m_s": PERMEABILITY * 1000.0 * 9.81 / 1.0e-3,
            "pipe_diameter_m": PIPE_DIAMETER,
            "cover_over_diameter": PIPE_COVER_RATIO,
            "pipe_invert_depth_m": 0.5 - (PIPE_CENTER[1] - PIPE_RADIUS),
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
