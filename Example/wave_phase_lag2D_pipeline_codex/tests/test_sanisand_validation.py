from __future__ import annotations

import csv
import math
import subprocess
import tempfile
import unittest
from pathlib import Path

import run_sanisand_validation as validation


class SanisandValidationTests(unittest.TestCase):
    def test_protocol_and_reference_identity_are_frozen(self) -> None:
        self.assertEqual(
            validation.EXPECTED_REFERENCE_COMMIT,
            "205c13b0a8fe5ffdcc1404b6fe63a59a67e613b9",
        )
        self.assertEqual(validation.INITIAL_PRESSURE_PA, 294000.0)
        self.assertEqual(validation.Q_AMPLITUDE_PA, 114200.0)
        self.assertEqual(validation.TIGHT_TOLERANCE, 1.0e-7)

    @staticmethod
    def monotonic_rows() -> list[dict]:
        return [
            {
                "implementation": validation.IMPLEMENTATION_CURRENT,
                "path_id": validation.PATH_MONOTONIC,
                "tolerance": validation.TIGHT_TOLERANCE,
                "step": step,
                "status": "ok",
            }
            for step in range(0, 1001, 5)
        ]

    def test_registered_row_grid_is_exact_and_missing_frame_fails(self) -> None:
        rows = self.monotonic_rows()
        validation.validate_row_grids(rows)
        with self.assertRaisesRegex(validation.ValidationError, "unexpected step grid"):
            validation.validate_row_grids(rows[:-1])

    def test_normalized_error_uses_only_common_ok_steps(self) -> None:
        left = [
            {"step": 0, "status": "ok", "p": 10.0},
            {"step": 1, "status": "ok", "p": 12.0},
            {"step": 2, "status": "controller_limit", "p": 99.0},
        ]
        right = [
            {"step": 0, "status": "ok", "p": 10.0},
            {"step": 1, "status": "ok", "p": 10.0},
            {"step": 2, "status": "ok", "p": 10.0},
        ]
        metric = validation.normalized_error(left, right, "p", 10.0)
        self.assertEqual(metric["paired_samples"], 2)
        self.assertAlmostEqual(metric["maximum_absolute_normalized"], 0.2)

    def test_terminal_record_ignores_nan_initial_controller_residual(self) -> None:
        rows = [
            {
                "status": "ok",
                "controller_residual_pa": math.nan,
                "step": 0,
                "cycle_time": 0.0,
                "p_effective_pa": 294000.0,
                "ru": 0.0,
                "axial_strain": 0.0,
            },
            {
                "status": "ok",
                "controller_residual_pa": 12.0,
                "step": 1,
                "cycle_time": 0.1,
                "p_effective_pa": 290000.0,
                "ru": 0.01,
                "axial_strain": -1e-4,
            },
            {
                "status": "controller_limit",
                "controller_residual_pa": 900.0,
                "step": 2,
                "cycle_time": 0.2,
                "p_effective_pa": 1000.0,
                "ru": 0.99,
                "axial_strain": 0.01,
            },
        ]
        record = validation.terminal_record(rows)
        self.assertEqual(record["maximum_accepted_controller_residual_pa"], 12.0)
        self.assertEqual(record["step"], 2)

    def test_q_target_history_drift_is_rejected(self) -> None:
        rows = []
        for step, cycle_time, status, residual in (
            (0, 0.0, "ok", 0.0),
            (1, 0.0025, "controller_limit", 1000.0),
        ):
            target = validation.Q_AMPLITUDE_PA * validation.triangle(cycle_time)
            rows.append(
                {
                    "implementation": validation.IMPLEMENTATION_CURRENT,
                    "path_id": validation.PATH_LITERATURE,
                    "tolerance": validation.TIGHT_TOLERANCE,
                    "step": step,
                    "cycle_time": cycle_time,
                    "axial_strain": 0.0,
                    "p_effective_pa": 294000.0,
                    "q_signed_pa": 0.0,
                    "ru": 0.0,
                    "void_ratio": 0.808,
                    "plastic_strain": 0.0,
                    "alpha_norm": 0.0,
                    "alpha_in_norm": 0.0,
                    "alpha_memory_norm": math.nan,
                    "memory_radius": math.nan,
                    "fabric_norm": 0.0,
                    "controller_target_q_pa": target,
                    "controller_residual_pa": residual,
                    "status": status,
                }
            )
        rows[1]["controller_target_q_pa"] = 0.0
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "driver.csv"
            with path.open("w", newline="", encoding="utf-8") as stream:
                writer = csv.DictWriter(stream, fieldnames=validation.EXPECTED_FIELDS)
                writer.writeheader()
                writer.writerows(rows)
            with self.assertRaisesRegex(
                validation.ValidationError, "target history drift"
            ):
                validation.load_rows(path, validation.IMPLEMENTATION_CURRENT)

    def test_artifact_validator_detects_equal_size_hash_drift(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "figure.png"
            path.write_bytes(b"abc")
            record = validation.artifact(path)
            validation.validate_artifacts([record])
            path.write_bytes(b"abd")
            with self.assertRaisesRegex(validation.ValidationError, "artifact drift"):
                validation.validate_artifacts([record])

    def test_reference_root_must_be_the_pinned_commit(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary) / "reference"
            root.mkdir()
            subprocess.run(["git", "init", "-q", str(root)], check=True)
            subprocess.run(
                ["git", "-C", str(root), "config", "user.email", "test@example.com"],
                check=True,
            )
            subprocess.run(
                ["git", "-C", str(root), "config", "user.name", "Test"],
                check=True,
            )
            (root / "README").write_text("not the reference\n", encoding="utf-8")
            subprocess.run(["git", "-C", str(root), "add", "README"], check=True)
            subprocess.run(
                ["git", "-C", str(root), "commit", "-q", "-m", "fixture"],
                check=True,
            )
            adapter = Path(validation.__file__).resolve().parent / (
                "sanisand_validation/pisano_udsm_material_point_driver.f90"
            )
            with tempfile.TemporaryDirectory() as build:
                with self.assertRaisesRegex(
                    validation.ValidationError, "reference commit must be"
                ):
                    validation.compile_reference(root, adapter, Path(build))

    def test_manuscript_does_not_mislabel_current_model_as_ms(self) -> None:
        manuscript = (
            Path(validation.__file__).resolve().parent
            / "SANISAND_VALIDATION_SUBSECTION.md"
        ).read_text(encoding="utf-8")
        self.assertIn("MPM SANISAND04-style", manuscript)
        self.assertIn("must not attribute", manuscript)
        self.assertIn("2.505", manuscript)
        self.assertIn("3.500", manuscript)
        self.assertIn("no numerical or display smoothing", manuscript)


if __name__ == "__main__":
    unittest.main()
