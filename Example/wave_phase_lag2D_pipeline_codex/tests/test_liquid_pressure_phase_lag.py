from __future__ import annotations

import copy
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np

import plot_liquid_pressure_phase_lag as phase_plot


class LiquidPressurePhaseLagTests(unittest.TestCase):
    @staticmethod
    def config(
        *,
        permeability: float = 3.0e-13,
        liquid_saturation: float = 0.937,
        pressure_smoothing: bool = False,
        pressure_smoothing_in_loop: bool = False,
        pipeline_fixed: bool = True,
    ) -> dict:
        return {
            "materials": [
                {"intrinsic_permeability": permeability},
                {
                    "liquid_saturation": liquid_saturation,
                    "gas_saturation": 1.0 - liquid_saturation,
                },
            ],
            "analysis": {
                "pressure_smoothing": pressure_smoothing,
                "pressure_smoothing_in_loop": pressure_smoothing_in_loop,
                "rigid_pipeline": {"fixed": pipeline_fixed},
            },
        }

    @staticmethod
    def artifact_fixture(
        directory: Path,
    ) -> tuple[dict, tuple[Path, Path, Path, Path], Path, Path]:
        analyzer = directory / "analyze_study.py"
        plot_script = directory / "plot_liquid_pressure_phase_lag.py"
        analyzer.write_bytes(b"analyzer source\n")
        plot_script.write_bytes(b"plot source\n")
        figures = (
            directory / "case_liquid_pressure_2d_phase_lag.png",
            directory / "case_liquid_pressure_2d_phase_lag.pdf",
            directory / "case_liquid_pressure_time_depth.png",
            directory / "case_liquid_pressure_time_depth.pdf",
        )
        for index, path in enumerate(figures):
            path.write_bytes(f"figure-{index}\n".encode("ascii"))
        analysis_config = directory / "03_HD.json"
        checkpoint_config = directory / "01_EQ.json"
        checkpoint = directory / "equilibrium/particle5000.vtp"
        particles = directory / "particles.txt"
        history = directory / "result/pipeline-history0000.csv"
        result = history.parent
        checkpoint.parent.mkdir()
        result.mkdir()
        analysis_config.write_bytes(b"analysis config\n")
        checkpoint_config.write_bytes(b"checkpoint config\n")
        checkpoint.write_bytes(b"checkpoint\n")
        particles.write_bytes(b"particles\n")
        history.write_bytes(b"history\n")
        frames = []
        for index in range(1, phase_plot.EXPECTED_PHASE_FRAMES + 1):
            frame = result / f"particle{650 * index:05d}.vtp"
            frame.write_bytes(f"frame-{index}\n".encode("ascii"))
            frames.append(frame)
        audit = {
            "schema": phase_plot.SCHEMA,
            "frame_count": phase_plot.EXPECTED_PHASE_FRAMES,
            "config": phase_plot.artifact_record(analysis_config),
            "checkpoint": phase_plot.artifact_record(checkpoint),
            "artifacts": {
                "schema": phase_plot.ARTIFACT_SCHEMA,
                "analyzer": phase_plot.artifact_record(analyzer),
                "plot_script": phase_plot.artifact_record(plot_script),
                "figures": [
                    phase_plot.artifact_record(path) for path in figures
                ],
                "inputs": {
                    "schema": phase_plot.INPUT_SCHEMA,
                    "analysis_config": phase_plot.artifact_record(analysis_config),
                    "checkpoint_config": phase_plot.artifact_record(
                        checkpoint_config
                    ),
                    "checkpoint": phase_plot.artifact_record(checkpoint),
                    "particles": phase_plot.artifact_record(particles),
                    "pipeline_history": phase_plot.artifact_record(history),
                    "result_directory": str(result.resolve()),
                    "result_particle_vtp": [
                        phase_plot.artifact_record(path) for path in frames
                    ],
                },
            },
        }
        return audit, figures, analyzer, plot_script

    def test_harmonic_fit_and_coherence_separate_phase_from_noise(self) -> None:
        period = 1.3
        times = np.linspace(0.0, 2.0 * period, 81)
        omega = 2.0 * np.pi / period
        exact = 25.0 + 80.0 * np.cos(omega * times - 0.7)
        noisy = exact + 200.0 * np.where(np.arange(times.size) % 2, 1.0, -1.0)
        values = np.column_stack((exact, noisy))

        _, amplitude, fitted_phase = phase_plot.harmonic_fields(
            times, values, period
        )
        coherence = phase_plot.harmonic_r_squared(times, values, period)

        self.assertAlmostEqual(amplitude[0], 80.0, places=10)
        self.assertAlmostEqual(fitted_phase[0], 0.7, places=10)
        self.assertGreater(coherence[0], 1.0 - 1.0e-12)
        self.assertLess(coherence[1], phase_plot.MINIMUM_HARMONIC_R_SQUARED)

    def test_same_column_reference_uses_highest_particle(self) -> None:
        coordinates = np.asarray(
            [[0.0, 0.1], [0.0, 0.3], [1.0, 0.2], [1.0, 0.4]], dtype=float
        )
        values = np.asarray([1.0, 2.0, 3.0, 4.0])
        mapped = phase_plot.same_column_surface_values(coordinates, values)
        np.testing.assert_array_equal(mapped, np.asarray([2.0, 2.0, 4.0, 4.0]))

    def test_wrap_degrees_uses_signed_principal_interval(self) -> None:
        values = np.radians(np.asarray([-540.0, -190.0, 190.0, 540.0]))
        np.testing.assert_allclose(
            phase_plot.wrap_degrees(values),
            np.asarray([-180.0, 170.0, -170.0, -180.0]),
            atol=1.0e-12,
        )

    def test_exploratory_parameters_and_subtitle_come_from_config(self) -> None:
        parameters = phase_plot.plot_parameters(
            self.config(), require_unsmoothed=True
        )
        self.assertEqual(
            parameters.as_dict(),
            {
                "intrinsic_permeability_m2": 3.0e-13,
                "liquid_saturation": 0.937,
                "gas_saturation": 1.0 - 0.937,
                "pressure_smoothing": False,
                "pressure_smoothing_in_loop": False,
                "pipeline_fixed": True,
            },
        )
        subtitle = phase_plot.parameter_subtitle(parameters)
        for expected in (
            "3.000e-13",
            "0.937",
            "0.063",
            "output=false",
            "in-loop=false",
            "pipeline fixed=true",
        ):
            self.assertIn(expected, subtitle)

    def test_formal_smoothed_mobile_config_remains_supported_by_default(self) -> None:
        parameters = phase_plot.plot_parameters(
            self.config(
                permeability=9.79e-12,
                liquid_saturation=0.94,
                pressure_smoothing=True,
                pipeline_fixed=False,
            )
        )
        self.assertTrue(parameters.pressure_smoothing)
        self.assertFalse(parameters.pressure_smoothing_in_loop)
        self.assertFalse(parameters.pipeline_fixed)

    def test_exploratory_gate_rejects_smoothing_or_mobile_pipeline(self) -> None:
        changes = (
            ("pressure_smoothing", True),
            ("pressure_smoothing_in_loop", True),
            ("fixed", False),
        )
        for name, value in changes:
            config = self.config()
            if name == "fixed":
                config["analysis"]["rigid_pipeline"][name] = value
            else:
                config["analysis"][name] = value
            with self.subTest(name=name), self.assertRaises(ValueError):
                phase_plot.plot_parameters(config, require_unsmoothed=True)

    def test_parameter_validation_rejects_nonboolean_and_saturation_error(self) -> None:
        nonboolean = self.config()
        nonboolean["analysis"]["pressure_smoothing"] = 0
        with self.assertRaisesRegex(ValueError, "explicit JSON boolean"):
            phase_plot.plot_parameters(nonboolean)

        inconsistent = copy.deepcopy(self.config())
        inconsistent["materials"][1]["gas_saturation"] = 0.2
        with self.assertRaisesRegex(ValueError, "saturations must sum to one"):
            phase_plot.plot_parameters(inconsistent)

    def test_artifact_validator_accepts_exact_bound_outputs(self) -> None:
        with TemporaryDirectory() as temporary:
            directory = Path(temporary)
            audit, figures, analyzer, plot_script = self.artifact_fixture(directory)

            phase_plot.validate_phase_audit_artifacts(
                audit,
                expected_figure_paths=figures,
                analyzer_path=analyzer,
                plot_script_path=plot_script,
            )
            audit_path = directory / "phase.audit.json"
            phase_plot.atomic_write_json(audit_path, audit)
            phase_plot.validate_phase_audit_artifacts(
                phase_plot.read_json(audit_path),
                expected_figure_paths=figures,
                analyzer_path=analyzer,
                plot_script_path=plot_script,
            )
            self.assertFalse(
                list(directory.glob(f".{audit_path.name}.*.tmp"))
            )

    def test_artifact_validator_rejects_figure_hash_drift_or_missing_file(self) -> None:
        for change in ("hash drift", "missing"):
            with self.subTest(change=change), TemporaryDirectory() as temporary:
                directory = Path(temporary)
                audit, figures, analyzer, plot_script = self.artifact_fixture(
                    directory
                )
                target = figures[0]
                if change == "hash drift":
                    target.write_bytes(b"x" * target.stat().st_size)
                    message = "sha256 mismatch"
                else:
                    target.unlink()
                    message = "missing"
                with self.assertRaisesRegex(phase_plot.PhaseAuditError, message):
                    phase_plot.validate_phase_audit_artifacts(
                        audit,
                        expected_figure_paths=figures,
                        analyzer_path=analyzer,
                        plot_script_path=plot_script,
                    )

    def test_artifact_validator_rejects_script_hash_drift_or_missing_file(self) -> None:
        for role in ("analyzer", "plot_script"):
            for change in ("hash drift", "missing"):
                with (
                    self.subTest(role=role, change=change),
                    TemporaryDirectory() as temporary,
                ):
                    directory = Path(temporary)
                    audit, figures, analyzer, plot_script = self.artifact_fixture(
                        directory
                    )
                    target = analyzer if role == "analyzer" else plot_script
                    if change == "hash drift":
                        target.write_bytes(b"x" * target.stat().st_size)
                        message = "sha256 mismatch"
                    else:
                        target.unlink()
                        message = "missing"
                    with self.assertRaisesRegex(
                        phase_plot.PhaseAuditError, message
                    ):
                        phase_plot.validate_phase_audit_artifacts(
                            audit,
                            expected_figure_paths=figures,
                            analyzer_path=analyzer,
                            plot_script_path=plot_script,
                        )

    def test_artifact_validator_rejects_size_and_figure_set_mismatch(self) -> None:
        with TemporaryDirectory() as temporary:
            directory = Path(temporary)
            audit, figures, analyzer, plot_script = self.artifact_fixture(directory)
            wrong_size = copy.deepcopy(audit)
            wrong_size["artifacts"]["figures"][0]["size_bytes"] += 1
            with self.assertRaisesRegex(phase_plot.PhaseAuditError, "size mismatch"):
                phase_plot.validate_phase_audit_artifacts(
                    wrong_size,
                    expected_figure_paths=figures,
                    analyzer_path=analyzer,
                    plot_script_path=plot_script,
                )

            duplicate = copy.deepcopy(audit)
            duplicate["artifacts"]["figures"][-1] = copy.deepcopy(
                duplicate["artifacts"]["figures"][0]
            )
            with self.assertRaisesRegex(phase_plot.PhaseAuditError, "unique"):
                phase_plot.validate_phase_audit_artifacts(
                    duplicate,
                    expected_figure_paths=figures,
                    analyzer_path=analyzer,
                    plot_script_path=plot_script,
                )

            unexpected = directory / "unexpected.png"
            unexpected.write_bytes(b"unexpected figure\n")
            wrong_set = copy.deepcopy(audit)
            wrong_set["artifacts"]["figures"][-1] = phase_plot.artifact_record(
                unexpected
            )
            with self.assertRaisesRegex(
                phase_plot.PhaseAuditError, "expected set"
            ):
                phase_plot.validate_phase_audit_artifacts(
                    wrong_set,
                    expected_figure_paths=figures,
                    analyzer_path=analyzer,
                    plot_script_path=plot_script,
                )

    def test_artifact_validator_binds_every_phase_input_and_exact_frame_set(self) -> None:
        for mutation, message in (
            ("checkpoint_hash", "sha256 mismatch"),
            ("missing_frame_record", "exactly 60"),
            ("unrecorded_frame", "exact ordered directory set"),
        ):
            with (
                self.subTest(mutation=mutation),
                TemporaryDirectory() as temporary,
            ):
                directory = Path(temporary)
                audit, figures, analyzer, plot_script = self.artifact_fixture(
                    directory
                )
                inputs = audit["artifacts"]["inputs"]
                if mutation == "checkpoint_hash":
                    Path(inputs["checkpoint"]["path"]).write_bytes(
                        b"x" * inputs["checkpoint"]["size_bytes"]
                    )
                elif mutation == "missing_frame_record":
                    inputs["result_particle_vtp"].pop()
                else:
                    result = Path(inputs["result_directory"])
                    (result / "particle40000.vtp").write_bytes(b"extra frame\n")
                with self.assertRaisesRegex(phase_plot.PhaseAuditError, message):
                    phase_plot.validate_phase_audit_artifacts(
                        audit,
                        expected_figure_paths=figures,
                        analyzer_path=analyzer,
                        plot_script_path=plot_script,
                    )


if __name__ == "__main__":
    unittest.main()
