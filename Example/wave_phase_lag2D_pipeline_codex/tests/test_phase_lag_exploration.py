from __future__ import annotations

import hashlib
import json
import shutil
import tempfile
import unittest
from pathlib import Path

import numpy as np

import prepare_phase_lag_exploration as exploration
import analyze_phase_lag_replay_exploration as replay_analysis


class PhaseLagExplorationTests(unittest.TestCase):
    def setUp(self) -> None:
        self.source_root = Path(__file__).resolve().parents[1]
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        (self.root / "configs/screen").mkdir(parents=True)
        for relative in (
            "configs/screen/01_EQ_HS.json",
            "configs/screen/02_HS.json",
            "particles.txt",
        ):
            destination = self.root / relative
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(self.source_root / relative, destination)

    def tearDown(self) -> None:
        self.temporary.cleanup()

    @staticmethod
    def read(path: Path) -> dict:
        return json.loads(path.read_text(encoding="utf-8"))

    @staticmethod
    def digest(path: Path) -> str:
        return hashlib.sha256(path.read_bytes()).hexdigest()

    def test_generates_strict_unsmoothed_paired_replay(self) -> None:
        manifest = exploration.generate(
            self.root,
            label="unit_k1e13_sw094",
            permeability_m2=1.0e-13,
            saturation=0.94,
            equilibrium_steps=5000,
            wave_steps=39000,
            pressure_smoothing=False,
            paired_replay=True,
        )
        directory = (
            self.root / "configs/phase_lag_exploratory/unit_k1e13_sw094"
        )
        equilibrium = self.read(directory / "01_EQ.json")
        wave = self.read(directory / "02_WAVE.json")
        cases = {
            role: self.read(directory / filename)
            for role, filename in (
                ("HD", "03_HD.json"),
                ("RL", "04_RL.json"),
                ("RE", "05_RE.json"),
            )
        }

        self.assertEqual(equilibrium["materials"][0]["type"], "LinearElastic2D")
        self.assertEqual(equilibrium["analysis"]["PIC"], 1.0)
        self.assertFalse(equilibrium["analysis"]["APIC"])
        self.assertEqual(
            equilibrium["analysis"]["damping"]["damping_factor"], 5.0
        )
        self.assertTrue(equilibrium["analysis"]["rigid_pipeline"]["fixed"])
        self.assertFalse(equilibrium["materials"][1]["wave_pressure"])
        self.assertEqual(equilibrium["analysis"]["nsteps"], 5000)

        for config in (equilibrium, wave, *cases.values()):
            self.assertEqual(config["materials"][0]["intrinsic_permeability"], 1e-13)
            self.assertEqual(config["materials"][1]["liquid_saturation"], 0.94)
            self.assertEqual(
                config["materials"][1]["gas_saturation"], 0.06000000000000005
            )
            self.assertFalse(config["analysis"]["pressure_smoothing"])
            self.assertFalse(config["analysis"]["pressure_smoothing_in_loop"])

        self.assertEqual(wave["materials"][0]["type"], "SANISAND2D")
        self.assertTrue(wave["materials"][1]["wave_pressure"])
        for role, config in cases.items():
            pressure = config["analysis"]["prescribed_phase_pressures"]
            self.assertTrue(config["analysis"]["rigid_pipeline"]["fixed"])
            self.assertEqual(pressure["step_interval"], 130)
            self.assertEqual(pressure["max_step"], 39000)
            self.assertEqual(pressure["mapping"], "id")
            self.assertEqual(config["analysis"]["nsteps"], 39000)
            if role == "HD":
                self.assertTrue(config["materials"][1]["wave_pressure"])
                self.assertTrue(pressure["write"])
                self.assertFalse(pressure["enable"])
            else:
                self.assertFalse(config["materials"][1]["wave_pressure"])
                self.assertFalse(pressure["write"])
                self.assertTrue(pressure["enable"])

        generated = manifest["generated"]
        self.assertEqual(
            generated["hd_sha256"], self.digest(directory / "03_HD.json")
        )
        self.assertEqual(
            generated["rl_sha256"], self.digest(directory / "04_RL.json")
        )
        self.assertEqual(
            generated["re_sha256"], self.digest(directory / "05_RE.json")
        )
        self.assertTrue(manifest["parameters"]["paired_replay"])
        self.assertEqual(
            manifest["parameters"]["pressure_database_interval_steps"], 130
        )

    def test_replay_analysis_helpers_are_finite_and_nan_safe(self) -> None:
        mesh = self.root / "small_mesh.txt"
        mesh.write_text(
            "#! elementShape quadrilateral\n"
            "4 1\n"
            "0 0\n"
            "2 0\n"
            "0 1\n"
            "2 1\n"
            "0 1 3 2\n",
            encoding="utf-8",
        )
        lower, upper = replay_analysis.mesh_bounds(mesh)
        self.assertEqual(lower.tolist(), [0.0, 0.0])
        self.assertEqual(upper.tolist(), [2.0, 1.0])
        self.assertAlmostEqual(
            replay_analysis.integrate(
                np.asarray([0.0, 1.0, 2.0]), np.asarray([0.0, 1.0, 0.0])
            ),
            1.0,
        )
        self.assertAlmostEqual(
            replay_analysis.robust_symmetric_limit(
                np.asarray([np.nan, -2.0, 2.0])
            ),
            2.0,
        )


if __name__ == "__main__":
    unittest.main()
