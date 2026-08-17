from __future__ import annotations

import json
import shutil
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import numpy as np

import qa_contract as qa
import run_phase_lag_exploration as runner
from phase_controls import FRAME_STEP, PressureHeader, write_header


class PhaseLagExplorationRunnerTests(unittest.TestCase):
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
        self.binary = self.root / "fake-mpm"
        self.binary.write_bytes(b"fake solver binary\n")
        self.binary.chmod(0o755)

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def prepare(self, label: str = "unit_k3e13") -> dict:
        return runner.prepare_stage(
            self.root,
            label=label,
            permeability_m2=3.0e-13,
            saturation=0.94,
            equilibrium_steps=5000,
            mpm_binary=self.binary,
        )

    def test_prepare_is_unsmoothed_hash_bound_and_single_use(self) -> None:
        state = self.prepare()
        self.assertEqual(state["stage"], "prepared")
        self.assertIs(state["parameters"]["pressure_smoothing"], False)
        self.assertEqual(state["solver_binary"]["sha256"], runner.sha256(self.binary))
        directory = self.root / "configs/phase_lag_exploratory/unit_k3e13"
        for name in runner.CONFIG_NAMES:
            config = json.loads((directory / name).read_text(encoding="utf-8"))
            self.assertIs(config["analysis"]["pressure_smoothing"], False)
            self.assertIs(config["analysis"]["pressure_smoothing_in_loop"], False)
        with self.assertRaisesRegex(runner.SafetyError, "single-use"):
            self.prepare()

    def test_generated_family_rejects_nonboolean_false_smoothing(self) -> None:
        for index, (key, value) in enumerate(
            (
                ("pressure_smoothing", 0),
                ("pressure_smoothing", None),
                ("pressure_smoothing_in_loop", 0),
                ("pressure_smoothing_in_loop", None),
            )
        ):
            label = f"unit_bad_smoothing_{index}"
            self.prepare(label)
            directory = self.root / "configs/phase_lag_exploratory" / label
            config_path = directory / "01_EQ.json"
            config = json.loads(config_path.read_text(encoding="utf-8"))
            config["analysis"][key] = value
            config_path.write_text(json.dumps(config), encoding="utf-8")
            manifest_path = directory / "manifest.json"
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            manifest["generated"]["equilibrium_sha256"] = runner.sha256(config_path)
            manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
            with self.subTest(key=key, value=value), self.assertRaisesRegex(
                runner.SafetyError, "pressure smoothing must be explicitly disabled"
            ):
                runner.validate_generated_family(
                    self.root,
                    label,
                    permeability_m2=3.0e-13,
                    saturation=0.94,
                    equilibrium_steps=5000,
                )

    def test_prepare_refuses_preexisting_partial_namespace(self) -> None:
        partial = self.root / "results/phase_lag_exploratory/partial"
        partial.mkdir(parents=True)
        (partial / "particle00650.vtp").write_bytes(b"partial")
        with self.assertRaisesRegex(runner.SafetyError, "single-use"):
            self.prepare("partial")
        self.assertFalse(
            (self.root / "configs/phase_lag_exploratory/partial").exists()
        )

    def test_provenance_tamper_and_hd_before_eq_are_hard_stops(self) -> None:
        self.prepare()
        with mock.patch.object(runner, "run_solver") as launch:
            with self.assertRaisesRegex(runner.SafetyError, "requires runner stage"):
                runner.solver_stage(
                    self.root, label="unit_k3e13", stage="hd", threads=1
                )
            launch.assert_not_called()
        config = (
            self.root
            / "configs/phase_lag_exploratory/unit_k3e13/03_HD.json"
        )
        config.write_text(config.read_text(encoding="utf-8") + " ", encoding="utf-8")
        state = runner.load_state(self.root, "unit_k3e13")
        with self.assertRaisesRegex(runner.SafetyError, "config 03_HD.json changed"):
            runner.verify_provenance(self.root, "unit_k3e13", state)

    def test_equilibrium_audit_keeps_original_velocity_gate(self) -> None:
        config = {
            "analysis": {"uuid": "UNIT_EQ", "nsteps": 5000},
            "post_processing": {
                "path": "results/unit/",
                "output_steps": 5000,
            },
        }
        config_path = self.root / "eq.json"
        config_path.write_text(json.dumps(config), encoding="utf-8")
        result = self.root / "results/unit/UNIT_EQ"
        result.mkdir(parents=True)
        for name in ("particle5000.vtp", "particles5000.h5", "pipeline-history0000.csv"):
            (result / name).write_bytes(b"nonempty")
        stability = {
            "observed": {"maximum_velocity_m_s": qa.MAXIMUM_VELOCITY_M_S + 1.0e-9}
        }
        with (
            mock.patch.object(runner.validate_case, "validate_final_vtp"),
            mock.patch.object(runner.validate_case, "validate_final_hdf5"),
            mock.patch.object(
                runner.validate_case,
                "validate_equilibrium_vtp",
                return_value=stability,
            ),
            mock.patch.object(
                runner.hdf5_checkpoint_crosscheck,
                "crosscheck",
                return_value={"passed": True},
            ),
        ):
            with self.assertRaisesRegex(runner.SafetyError, "unchanged 1e-3"):
                runner.audit_equilibrium(self.root, config_path)

    @staticmethod
    def _hd_arrays(count: int) -> dict[str, np.ndarray]:
        return {
            "volumes": np.full(count, 1.0e-4),
            "velocities": np.zeros((count, 3)),
            "displacements": np.zeros((count, 3)),
            "porosities": np.full(count, 0.485),
        }

    def _synthetic_hd(self) -> tuple[Path, np.ndarray, dict[str, np.ndarray]]:
        particles = self.root / "tiny_particles.txt"
        particles.write_text("2\n0.25 0.25\n0.75 0.75\n", encoding="utf-8")
        mesh = self.root / "tiny_mesh.txt"
        mesh.write_text(
            "#! elementShape quadrilateral\n"
            "4 1\n0 0\n1 0\n0 1\n1 1\n0 1 3 2\n",
            encoding="utf-8",
        )
        config = {
            "analysis": {
                "uuid": "UNIT_HD",
                "nsteps": runner.WAVE_STEPS,
                "dt": runner.DT,
                "prescribed_phase_pressures": {
                    "enable": False,
                    "write": True,
                    "path": "pressure_databases/unit/lagged",
                    "file_prefix": "pressure",
                    "source_dt": runner.DT,
                    "step_interval": runner.PRESSURE_INTERVAL,
                    "max_step": runner.WAVE_STEPS,
                },
            },
            "post_processing": {
                "path": "results/unit/",
                "output_steps": runner.WAVE_OUTPUT_INTERVAL,
            },
            "particles": [{"generator": {"location": particles.name}}],
            "mesh": {"mesh": mesh.name},
        }
        config_path = self.root / "hd.json"
        config_path.write_text(json.dumps(config), encoding="utf-8")
        result = self.root / "results/unit/UNIT_HD"
        result.mkdir(parents=True)
        for step in range(runner.WAVE_OUTPUT_INTERVAL, runner.WAVE_STEPS + 1, runner.WAVE_OUTPUT_INTERVAL):
            (result / f"particle{step:05d}.vtp").write_bytes(b"frame")
        (result / "pipeline-history00000.csv").write_bytes(b"history")

        database = self.root / "pressure_databases/unit/lagged"
        database.mkdir(parents=True)
        points = database / "pressure_points.txt"
        points.write_text("0 0.25 0.25\n1 0.75 0.75\n", encoding="utf-8")
        header = PressureHeader(
            dimension=2,
            values_per_particle=2,
            particle_count=2,
            step_interval=runner.PRESSURE_INTERVAL,
            max_step=runner.WAVE_STEPS,
            source_dt=runner.DT,
            format_version="V2",
        )
        zeros = np.zeros(2, dtype="<f8").tobytes()
        with (database / "pressure_values.bin").open("wb") as stream:
            write_header(stream, header)
            for step in range(0, runner.WAVE_STEPS + 1, runner.PRESSURE_INTERVAL):
                stream.write(FRAME_STEP.pack(step))
                stream.write(zeros)
                stream.write(zeros)
        return config_path, np.asarray([[0.25, 0.25], [0.75, 0.75]]), self._hd_arrays(2)

    def test_hd_audit_requires_exact_grid_complete_database_and_inside_mesh(self) -> None:
        config_path, points, arrays = self._synthetic_hd()
        with (
            mock.patch.object(
                runner, "ordered_vtp_frame", return_value=(points, arrays)
            ),
            mock.patch.object(runner, "read_pipeline_history", return_value=[]),
            mock.patch.object(
                runner,
                "validate_pipeline_history_grid",
                return_value={"rows": 39001},
            ),
        ):
            audit = runner.audit_hd(self.root, config_path)
        self.assertEqual(audit["health"]["saved_frames"], 60)
        self.assertEqual(audit["health"]["outside_mesh_maximum_count"], 0)
        self.assertEqual(audit["pressure_database"]["header"]["frame_count"], 301)

        outside = points.copy()
        outside[0, 0] = 1.1
        with (
            mock.patch.object(
                runner, "ordered_vtp_frame", return_value=(outside, arrays)
            ),
            mock.patch.object(runner, "read_pipeline_history", return_value=[]),
            mock.patch.object(
                runner,
                "validate_pipeline_history_grid",
                return_value={"rows": 39001},
            ),
        ):
            with self.assertRaisesRegex(runner.SafetyError, "outside the mesh"):
                runner.audit_hd(self.root, config_path)


if __name__ == "__main__":
    unittest.main()
