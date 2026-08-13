from __future__ import annotations

import copy
import json
import math
import sys
import tempfile
import unittest
from pathlib import Path


CASE_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(CASE_DIR))

import prepare_study as study  # noqa: E402


class PrepareStudyTest(unittest.TestCase):
    def test_mc_friction_matches_sanisand_M(self):
        sine_phi = math.sin(math.radians(study.MC_FRICTION_DEG))
        recovered_M = 6.0 * sine_phi / (3.0 - sine_phi)
        self.assertAlmostEqual(recovered_M, study.SANISAND_MC, places=12)

    def test_mc_stiffness_matches_sanisand_tangent_at_three_kpa(self):
        shear, bulk = study.sanisand_elastic_moduli(3_000.0, study.POROSITY)
        self.assertAlmostEqual(shear, 4_586_927.52449314, places=8)
        self.assertAlmostEqual(bulk, 2_986_049.6431101956, places=8)

        material = study.mohr_coulomb_material()
        youngs = material["youngs_modulus"]
        poisson = material["poisson_ratio"]
        recovered_shear = youngs / (2.0 * (1.0 + poisson))
        recovered_bulk = youngs / (3.0 * (1.0 - 2.0 * poisson))
        self.assertAlmostEqual(recovered_shear, shear, places=8)
        self.assertAlmostEqual(recovered_bulk, bulk, places=8)
        self.assertGreater(poisson, -1.0)
        self.assertLess(poisson, 0.5)
        self.assertEqual(material["stiffness_calibration"], study.mc_stiffness_calibration())

        sanisand = study.sanisand_material()
        self.assertNotIn("youngs_modulus", sanisand)
        self.assertNotIn("poisson_ratio", sanisand)

    def test_generated_matrix_is_structurally_valid(self):
        with tempfile.TemporaryDirectory() as temporary:
            manifest = study.generate_tier(
                "screen",
                Path(temporary),
                mesh_directory=".",
                cell_size=0.02,
                particle_spacing=0.01,
                wave_height=0.12,
                label="baseline",
            )
            self.assertTrue(manifest.exists())
            self.assertEqual(len(list(manifest.parent.glob("*.json"))), 11)

            def load(name: str):
                return json.loads((manifest.parent / name).read_text(encoding="utf-8"))

            low = load("02_LS.json")
            high = load("02_HS.json")
            mc_handoff = load("02_MC_EQ.json")
            mc = load("03_HM.json")
            equilibrium = load("01_EQ_LS.json")
            equilibrium_high = load("01_EQ_HS.json")
            for config in (equilibrium, low, high, mc_handoff, mc):
                self.assertIn("volumes", config["post_processing"]["vtk"])
                self.assertIn("ids", config["post_processing"]["vtk"])
            self.assertEqual(equilibrium["materials"][0]["type"], "LinearElastic2D")
            self.assertAlmostEqual(
                equilibrium["materials"][0]["youngs_modulus"],
                study.MC_MATCHED_YOUNGS_MODULUS,
            )
            self.assertAlmostEqual(
                equilibrium["materials"][0]["poisson_ratio"],
                study.MC_MATCHED_POISSON_RATIO,
            )
            self.assertFalse(equilibrium["analysis"]["APIC"])
            self.assertEqual(equilibrium["analysis"]["PIC"], 1.0)
            self.assertEqual(
                equilibrium["analysis"]["minimum_nodal_support_fraction"],
                study.MINIMUM_NODAL_SUPPORT_FRACTION,
            )
            self.assertTrue(low["analysis"]["APIC"])
            self.assertEqual(
                low["analysis"]["minimum_nodal_support_fraction"],
                study.MINIMUM_NODAL_SUPPORT_FRACTION,
            )
            self.assertEqual(
                mc_handoff["materials"][0]["type"], "MohrCoulomb2D"
            )
            self.assertTrue(mc_handoff["analysis"]["handoff_relaxation"])
            self.assertTrue(mc_handoff["analysis"]["stability_gate"])
            self.assertFalse(mc_handoff["analysis"]["APIC"])
            self.assertEqual(mc_handoff["analysis"]["PIC"], 1.0)
            self.assertTrue(mc_handoff["analysis"]["rigid_pipeline"]["fixed"])
            self.assertFalse(mc_handoff["materials"][1]["wave_pressure"])
            self.assertEqual(
                mc_handoff["analysis"]["damping"]["damping_factor"],
                study.EQUILIBRIUM_DAMPING_FACTOR,
            )
            self.assertEqual(
                mc_handoff["analysis"]["nsteps"], study.MC_HANDOFF_STEPS
            )
            self.assertTrue(mc_handoff["post_processing"]["write_hdf5"])
            self.assertEqual(
                mc_handoff["analysis"]["resume"]["uuid"],
                equilibrium_high["analysis"]["uuid"],
            )
            self.assertTrue(mc["analysis"]["APIC"])
            self.assertEqual(mc["analysis"]["PIC"], 0.0)
            self.assertEqual(
                mc["analysis"]["resume"]["uuid"],
                mc_handoff["analysis"]["uuid"],
            )
            self.assertNotIn("youngs_modulus", low["materials"][0])
            self.assertEqual(
                low["materials"][0]["critical_timestep_modulus"],
                study.CRITICAL_TIMESTEP_MODULUS,
            )
            self.assertEqual(
                equilibrium["analysis"]["damping"]["damping_factor"],
                study.EQUILIBRIUM_DAMPING_FACTOR,
            )
            self.assertEqual(equilibrium["analysis"]["nsteps"], 40_000)
            self.assertEqual(
                low["materials"][1]["liquid_saturation"],
                study.LOW_LAG_SATURATION,
            )
            self.assertEqual(
                equilibrium["mesh"]["particles_stresses"],
                "initial_effective_stresses_LS.txt",
            )
            self.assertEqual(
                equilibrium_high["mesh"]["particles_stresses"],
                "initial_effective_stresses_HS.txt",
            )
            self.assertEqual(
                equilibrium["mesh"]["particles_pore_pressures"],
                {"file": "initial_liquid_pressures.txt"},
            )
            expected_surface_loading = [
                {
                    "pset_id": study.SURFACE_TRACTION_PSET_ID,
                    "dir": 1,
                    "traction": study.SUBMERGED_SURFACE_TRACTION,
                    "facet": 2,
                }
            ]
            for config in (
                equilibrium,
                equilibrium_high,
                low,
                high,
                mc_handoff,
                mc,
            ):
                self.assertEqual(
                    config["external_loading_conditions"][
                        "particle_surface_traction"
                    ],
                    expected_surface_loading,
                )
                self.assertEqual(
                    config["mesh"]["boundary_conditions"][
                        "particles_at_free_surface"
                    ],
                    [
                        {
                            "pset_id": study.PHYSICAL_SURFACE_PSET_ID,
                            "nonfree_pset_id": 3,
                        }
                    ],
                )
            self.assertEqual(study.SUBMERGED_SURFACE_TRACTION, -4_905.0)
            self.assertNotIn("particles_stresses", low["mesh"])
            self.assertNotIn("particles_pore_pressures", low["mesh"])
            for config in (low, high, mc):
                self.assertIn(
                    "vertical_effective_stress_remaining_ratios",
                    config["post_processing"]["vtk"],
                )
                self.assertIn(
                    "liquid_densities",
                    config["post_processing"]["liquid_vtk"],
                )
                self.assertIn(
                    "PIC_pore_pressure_excess",
                    config["post_processing"]["liquid_vtk"],
                )
            self.assertNotIn(
                "PIC_pore_pressure_excess",
                equilibrium["post_processing"]["liquid_vtk"],
            )
            low_fluid = copy.deepcopy(low["materials"][1])
            high_fluid = copy.deepcopy(high["materials"][1])
            for fluid in (low_fluid, high_fluid):
                fluid.pop("liquid_saturation")
                fluid.pop("gas_saturation")
            self.assertEqual(low_fluid, high_fluid)
            self.assertEqual(low["materials"][0], high["materials"][0])
            self.assertAlmostEqual(
                mc["materials"][0]["youngs_modulus"],
                study.MC_MATCHED_YOUNGS_MODULUS,
            )
            self.assertAlmostEqual(
                mc["materials"][0]["poisson_ratio"],
                study.MC_MATCHED_POISSON_RATIO,
            )
            manifest_data = json.loads(manifest.read_text(encoding="utf-8"))
            self.assertEqual(
                manifest_data["constants"]["mc_stiffness_calibration"],
                study.mc_stiffness_calibration(),
            )
            self.assertEqual(
                manifest_data["constants"]["equilibrium_steps"],
                study.EQUILIBRIUM_STEPS,
            )
            self.assertEqual(
                manifest_data["constants"]["equilibrium_damping_factor_per_s"],
                study.EQUILIBRIUM_DAMPING_FACTOR,
            )
            self.assertEqual(
                manifest_data["constants"]["mc_handoff_steps"],
                study.MC_HANDOFF_STEPS,
            )
            self.assertEqual(
                manifest_data["constants"]["stability_max_velocity_m_s"],
                study.STABILITY_MAX_VELOCITY,
            )

            lagged = load("04_RL.json")
            matched_mc = load("04_RM.json")
            erased = load("04_RE.json")
            for config in (lagged, erased):
                config.pop("title")
                config["analysis"].pop("uuid")
                pressure = config["analysis"]["prescribed_phase_pressures"]
                pressure.pop("path")
                pressure.pop("file_prefix")
            self.assertEqual(lagged, erased)

            self.assertEqual(matched_mc["materials"][0]["type"], "MohrCoulomb2D")
            self.assertEqual(matched_mc["materials"][1], load("04_RL.json")["materials"][1])
            self.assertEqual(
                matched_mc["analysis"]["resume"]["uuid"],
                mc_handoff["analysis"]["uuid"],
            )
            self.assertEqual(
                matched_mc["analysis"]["prescribed_phase_pressures"],
                load("04_RL.json")["analysis"]["prescribed_phase_pressures"],
            )
            manifest_roles = {
                case["code"]: (case["role"], case["dependency"])
                for case in manifest_data["cases"]
            }
            self.assertEqual(
                manifest_roles["RM"],
                ("matched-pressure constitutive ablation", "HD + MC_EQ"),
            )
            self.assertEqual(
                manifest_roles["MC_EQ"],
                ("MC handoff relaxation", "EQ_HS"),
            )

    def test_equilibrium_ids_are_accepted_and_written_by_base_solver(self):
        source_path = (
            CASE_DIR.parents[1]
            / "mpm_hpc_source"
            / "include"
            / "solvers"
            / "mpm_base.tcc"
        )
        source = source_path.read_text(encoding="utf-8")
        initialise = source[
            source.index("bool mpm::MPMBase<Tdim>::initialise_vtk()") :
            source.index("bool mpm::MPMBase<Tdim>::initialise_vtk_twophase()")
        ]
        writer = source[source.index("void mpm::MPMBase<Tdim>::write_vtk(") :]
        scalar_list = writer[
            writer.index("std::vector<std::string> vtk_scalar_data") :
            writer.index("//! VTK vector variable")
        ]
        self.assertIn('"ids"', initialise)
        self.assertIn('"ids"', scalar_list)

    def test_threephase_initial_pressure_survives_compute_mass_initialisation(self):
        header = (
            CASE_DIR.parents[1]
            / "mpm_hpc_source"
            / "include"
            / "particles"
            / "particle_threephase_lag.h"
        ).read_text(encoding="utf-8")
        source = (
            CASE_DIR.parents[1]
            / "mpm_hpc_source"
            / "include"
            / "particles"
            / "particle_threephase_lag.tcc"
        ).read_text(encoding="utf-8")
        self.assertIn("has_input_initial_liquid_pressure_ = true", header)
        self.assertIn("has_input_initial_liquid_pressure_ ?", source)
        self.assertIn("input_initial_liquid_pressure_", source)

    def test_primary_stress_ratio_is_registered_for_base_vtk_validation(self):
        source = (
            CASE_DIR.parents[1]
            / "mpm_hpc_source"
            / "include"
            / "solvers"
            / "mpm_base.tcc"
        ).read_text(encoding="utf-8")
        self.assertGreaterEqual(
            source.count('"vertical_effective_stress_remaining_ratios"'), 2
        )
        for field in ("phi", "psi", "cohesion", "pdstrain"):
            self.assertGreaterEqual(source.count(f'"{field}"'), 2)

    def test_cundall_damping_is_applied_to_pure_pic_velocities(self):
        particle_source = (
            CASE_DIR.parents[1]
            / "mpm_hpc_source"
            / "include"
            / "particles"
            / "particle.tcc"
        ).read_text(encoding="utf-8")
        threephase_source = (
            CASE_DIR.parents[1]
            / "mpm_hpc_source"
            / "include"
            / "particles"
            / "particle_threephase_lag.tcc"
        ).read_text(encoding="utf-8")
        self.assertIn(
            "pic_velocity -= damping_factor * this->velocity_ * dt",
            particle_source,
        )
        self.assertIn(
            "pic_liquid_velocity -= damping_factor * this->liquid_velocity_ * dt",
            threephase_source,
        )
        self.assertIn(
            "pic_gas_velocity -= damping_factor * this->gas_velocity_ * dt",
            threephase_source,
        )

    def test_validator_rejects_old_slash_apic_key(self):
        config = study.equilibrium_config(
            "TEST",
            0.94,
            "results/test/",
            ".",
            0.02,
            0.01,
            0.12,
        )
        config = copy.deepcopy(config)
        config["analysis"]["/APIC"] = config["analysis"].pop("APIC")
        with self.assertRaisesRegex(ValueError, "APIC"):
            study.validate_config(config)

    def test_validator_requires_pure_pic_equilibrium(self):
        config = study.equilibrium_config(
            "TEST", 0.94, "results/test/", ".", 0.02, 0.01, 0.12
        )
        config["analysis"]["APIC"] = True
        with self.assertRaisesRegex(ValueError, "APIC=false"):
            study.validate_config(config)

    def test_validator_requires_pure_pic_mc_handoff(self):
        config = study.mc_handoff_config(
            code="TEST",
            saturation=study.HIGH_LAG_SATURATION,
            equilibrium_uuid="EQ",
            result_path="results/test/",
            mesh_directory=".",
            cell_size=0.02,
            particle_spacing=0.01,
            wave_height=0.12,
        )
        config["analysis"]["APIC"] = True
        with self.assertRaisesRegex(ValueError, "MC handoff"):
            study.validate_config(config)

    def test_validator_requires_primary_liquefaction_outputs(self):
        config = study.dynamic_config(
            code="SCREEN_HS_SANI",
            title="test",
            saturation=0.94,
            material_type="SANISAND2D",
            equilibrium_uuid="PLP_SCREEN_HS_EQ",
            result_path="results/test/",
            nsteps=study.TIERS["screen"]["cycles"] * study.STEPS_PER_CYCLE,
            mesh_directory=".",
            cell_size=0.02,
            particle_spacing=0.01,
            wave_height=0.12,
            physical_wave=True,
        )
        config["post_processing"]["liquid_vtk"].remove(
            "PIC_pore_pressure_excess"
        )
        with self.assertRaisesRegex(ValueError, "primary liquefaction outputs"):
            study.validate_config(config)

    def test_validator_rejects_registered_stage_drift(self):
        with tempfile.TemporaryDirectory() as temporary:
            manifest = study.generate_tier(
                "screen",
                Path(temporary),
                mesh_directory=".",
                cell_size=0.02,
                particle_spacing=0.01,
                wave_height=0.12,
                label="baseline",
            )

            def load(name: str) -> dict:
                return json.loads((manifest.parent / name).read_text(encoding="utf-8"))

            mutations = [
                ("01_EQ_LS.json", lambda c: c["analysis"].__setitem__("dt", 2.0e-4)),
                (
                    "01_EQ_LS.json",
                    lambda c: c["analysis"]["rigid_pipeline"].__setitem__("fixed", False),
                ),
                (
                    "01_EQ_LS.json",
                    lambda c: c["materials"][1].__setitem__("wave_pressure", True),
                ),
                (
                    "01_EQ_LS.json",
                    lambda c: c["materials"][1].__setitem__("liquid_saturation", 0.98),
                ),
                (
                    "01_EQ_LS.json",
                    lambda c: c["analysis"]["stability_qa_contract"][
                        "limits"
                    ].__setitem__("maximum_velocity_m_s", 2.0e-3),
                ),
                (
                    "01_EQ_LS.json",
                    lambda c: c["mesh"]["boundary_conditions"].__setitem__(
                        "particles_at_free_surface",
                        [{"pset_id": 0, "nonfree_pset_id": 3}],
                    ),
                ),
                (
                    "01_EQ_LS.json",
                    lambda c: c["external_loading_conditions"][
                        "particle_surface_traction"
                    ][0].__setitem__("pset_id", 0),
                ),
                (
                    "02_MC_EQ.json",
                    lambda c: c["analysis"]["resume"].__setitem__("uuid", "WRONG_EQ"),
                ),
                ("02_HS.json", lambda c: c["analysis"].__setitem__("nsteps", 13_000)),
                ("02_HS.json", lambda c: c["analysis"].__setitem__("PIC", 0.1)),
                (
                    "02_HS.json",
                    lambda c: c["analysis"].pop(
                        "resume_stability_qa_contract"
                    ),
                ),
                (
                    "02_HS.json",
                    lambda c: c["analysis"]["damping"].__setitem__(
                        "damping_factor", 1.0
                    ),
                ),
                (
                    "02_HS.json",
                    lambda c: c["materials"][0].__setitem__("type", "MohrCoulomb2D"),
                ),
                (
                    "02_HS.json",
                    lambda c: c["analysis"]["rigid_pipeline"].__setitem__(
                        "release_time", 4.0
                    ),
                ),
                (
                    "02_HS.json",
                    lambda c: c["materials"][1].__setitem__("wave_pressure", False),
                ),
                (
                    "02_HD.json",
                    lambda c: c["analysis"]["rigid_pipeline"].__setitem__("fixed", False),
                ),
                (
                    "04_RL.json",
                    lambda c: c["analysis"].pop("prescribed_phase_pressures"),
                ),
                (
                    "04_RE.json",
                    lambda c: c["analysis"]["prescribed_phase_pressures"].__setitem__(
                        "file_prefix", "pressure"
                    ),
                ),
                (
                    "04_RM.json",
                    lambda c: c["materials"][0].__setitem__("type", "SANISAND2D"),
                ),
            ]
            for filename, mutate in mutations:
                with self.subTest(filename=filename, mutation=mutate):
                    config = load(filename)
                    mutate(config)
                    with self.assertRaises(ValueError):
                        study.validate_config(config)


if __name__ == "__main__":
    unittest.main()
