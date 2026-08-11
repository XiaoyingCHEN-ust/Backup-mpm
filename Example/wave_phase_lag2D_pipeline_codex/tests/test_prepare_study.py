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
            self.assertEqual(len(list(manifest.parent.glob("*.json"))), 10)

            def load(name: str):
                return json.loads((manifest.parent / name).read_text(encoding="utf-8"))

            low = load("02_LS.json")
            high = load("02_HS.json")
            mc = load("02_HM.json")
            equilibrium = load("01_EQ_LS.json")
            for config in (equilibrium, low, high, mc):
                self.assertIn("volumes", config["post_processing"]["vtk"])
                self.assertIn("ids", config["post_processing"]["vtk"])
            for config in (low, high, mc):
                self.assertIn(
                    "vertical_effective_stress_remaining_ratios",
                    config["post_processing"]["vtk"],
                )
                self.assertIn(
                    "liquid_densities",
                    config["post_processing"]["liquid_vtk"],
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
                matched_mc["analysis"]["resume"],
                load("04_RL.json")["analysis"]["resume"],
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
                ("matched-pressure constitutive ablation", "HD"),
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

    def test_validator_requires_primary_liquefaction_outputs(self):
        config = study.dynamic_config(
            code="TEST",
            title="test",
            saturation=0.94,
            material_type="SANISAND2D",
            equilibrium_uuid="EQ",
            result_path="results/test/",
            nsteps=study.STEPS_PER_CYCLE,
            mesh_directory=".",
            cell_size=0.02,
            particle_spacing=0.01,
            wave_height=0.12,
            physical_wave=True,
        )
        config["post_processing"]["liquid_vtk"].remove(
            "liquid_seepage_forces"
        )
        with self.assertRaisesRegex(ValueError, "primary liquefaction outputs"):
            study.validate_config(config)


if __name__ == "__main__":
    unittest.main()
