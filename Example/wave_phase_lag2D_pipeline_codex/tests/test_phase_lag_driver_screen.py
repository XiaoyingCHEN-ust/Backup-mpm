from __future__ import annotations

import copy
import json
import tempfile
import unittest
from pathlib import Path

import numpy as np

import analyze_phase_lag_driver_screen as screen
from phase_controls import rotated_coefficients


def row(
    time_s: float,
    *,
    lagged_force: float,
    erased_force: float,
    lagged_if_area: float,
    erased_if_area: float,
    lagged_joint: float,
    erased_joint: float,
) -> dict[str, float]:
    return {
        "time_s": time_s,
        "lagged_signed_force_n_per_m": lagged_force,
        "phase_erased_signed_force_n_per_m": erased_force,
        "lagged_net_uplift_force_n_per_m": max(lagged_force, 0.0),
        "phase_erased_net_uplift_force_n_per_m": max(erased_force, 0.0),
        "lagged_positive_force_n_per_m": lagged_force,
        "phase_erased_positive_force_n_per_m": erased_force,
        "lagged_critical_area_m2": lagged_if_area,
        "phase_erased_critical_area_m2": erased_if_area,
        "lagged_joint_area_m2": lagged_joint,
        "phase_erased_joint_area_m2": erased_joint,
    }


class PhaseLagDriverScreenTests(unittest.TestCase):
    @staticmethod
    def qualified_phase() -> dict[str, object]:
        return {
            "passed": True,
            "checks": {"synthetic": True},
            "observed": {},
        }

    def _write_phase_field_audit(
        self, root: Path
    ) -> tuple[
        Path,
        Path,
        screen.replay_analysis.ReplayParameters,
        dict,
    ]:
        label = "case"
        case = label.upper()
        phase_dir = (
            root / "analysis/phase_lag_exploratory" / label / "phase_v3"
        )
        phase_dir.mkdir(parents=True)
        hd_path = root / "configs/phase_lag_exploratory" / label / "03_HD.json"
        hd_path.parent.mkdir(parents=True)
        hd_path.write_text("{}\n", encoding="utf-8")
        eq_path = hd_path.with_name("01_EQ.json")
        eq_path.write_text("{}\n", encoding="utf-8")
        particles = root / "particles.txt"
        particles.write_text("particles\n", encoding="utf-8")
        eq_result = root / "results/case/EQ"
        hd_result = root / "results/case/HD"
        eq_result.mkdir(parents=True)
        hd_result.mkdir(parents=True)
        checkpoint = eq_result / "particle5000.vtp"
        checkpoint.write_text("checkpoint\n", encoding="utf-8")
        history = hd_result / "pipeline-history0000.csv"
        history.write_text("history\n", encoding="utf-8")
        frames = []
        for index in range(1, screen.EXPECTED_FRAMES + 1):
            frame = hd_result / f"particle{650 * index:05d}.vtp"
            frame.write_text(f"frame {index}\n", encoding="utf-8")
            frames.append(frame)
        parameters = screen.replay_analysis.ReplayParameters(
            intrinsic_permeability_m2=1.0e-13,
            liquid_saturation=0.94,
            gas_saturation=0.06,
            pressure_smoothing=False,
            pressure_smoothing_in_loop=False,
            pipeline_fixed=True,
        )
        figure_paths = tuple(
            phase_dir / f"{case}_{stem}{suffix}"
            for stem in (
                "liquid_pressure_2d_phase_lag",
                "liquid_pressure_time_depth",
            )
            for suffix in (".png", ".pdf")
        )
        for index, figure_path in enumerate(figure_paths):
            figure_path.write_bytes(f"synthetic figure {index}\n".encode())
        audit = {
            "schema": screen.phase_plot.SCHEMA,
            "case": case,
            "parameters": parameters.as_dict(),
            "require_unsmoothed": True,
            "particle_count": 7376,
            "frame_count": screen.EXPECTED_FRAMES,
            "fit_window_s": [screen.FIT_START_S, screen.FIT_END_S],
            "wave_period_s": screen.PERIOD_S,
            "config": screen.phase_plot.artifact_record(hd_path),
            "checkpoint": screen.phase_plot.artifact_record(checkpoint),
            "probe_results": {
                "surface": {"fundamental_amplitude_pa": 100.0},
                "crown": {
                    "fundamental_amplitude_pa": 10.0,
                    "harmonic_r_squared": 0.95,
                    "phase_difference_from_same_column_surface_deg": -20.0,
                },
            },
            "artifacts": {
                "schema": screen.phase_plot.ARTIFACT_SCHEMA,
                "analyzer": screen.phase_plot.artifact_record(
                    screen.phase_plot.ANALYZER_PATH
                ),
                "plot_script": screen.phase_plot.artifact_record(
                    screen.phase_plot.PLOT_SCRIPT_PATH
                ),
                "figures": [
                    screen.phase_plot.artifact_record(path)
                    for path in figure_paths
                ],
                "inputs": {
                    "schema": screen.phase_plot.INPUT_SCHEMA,
                    "analysis_config": screen.phase_plot.artifact_record(hd_path),
                    "checkpoint_config": screen.phase_plot.artifact_record(eq_path),
                    "checkpoint": screen.phase_plot.artifact_record(checkpoint),
                    "particles": screen.phase_plot.artifact_record(particles),
                    "pipeline_history": screen.phase_plot.artifact_record(history),
                    "result_directory": str(hd_result.resolve()),
                    "result_particle_vtp": [
                        screen.phase_plot.artifact_record(path) for path in frames
                    ],
                },
            },
        }
        runner_state = {
            "stages": {
                "eq": {
                    "config": screen.phase_plot.artifact_record(eq_path),
                    "audit": {
                        "final_vtp": screen.phase_plot.artifact_record(checkpoint),
                        "stability_qa": {
                            "normal_stress": {
                                "inputs": {
                                    "particles": {
                                        "path": "particles.txt",
                                        "size_bytes": particles.stat().st_size,
                                        "sha256": screen.sha256(particles),
                                    }
                                }
                            }
                        },
                    },
                },
                "hd": {
                    "config": screen.phase_plot.artifact_record(hd_path),
                    "audit": {
                        "result_directory": str(hd_result.resolve()),
                        "particle_vtp": [
                            screen.phase_plot.artifact_record(path) for path in frames
                        ],
                        "pipeline_history": screen.phase_plot.artifact_record(history),
                    },
                },
            }
        }
        audit_path = phase_dir / f"{case}_liquid_pressure_phase_lag.audit.json"
        audit_path.write_text(json.dumps(audit), encoding="utf-8")
        return audit_path, hd_path, parameters, runner_state

    def test_phase_field_audit_accepts_complete_artifact_set(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            _audit_path, hd_path, parameters, runner_state = (
                self._write_phase_field_audit(root)
            )
            result = screen.validate_phase_field_audit(
                root, "case", hd_path, parameters, runner_state
            )
            self.assertTrue(result["qualification"]["passed"])
            self.assertEqual(
                result["runner_binding"]["schema"],
                screen.PHASE_RUNNER_BINDING_SCHEMA,
            )

    def test_phase_field_audit_artifact_failures_are_screen_errors(self) -> None:
        for mutation in ("missing_artifacts", "missing_figure", "script_drift"):
            with (
                self.subTest(mutation=mutation),
                tempfile.TemporaryDirectory() as temporary,
            ):
                root = Path(temporary)
                audit_path, hd_path, parameters, runner_state = (
                    self._write_phase_field_audit(root)
                )
                audit = json.loads(audit_path.read_text(encoding="utf-8"))
                if mutation == "missing_artifacts":
                    audit.pop("artifacts")
                    audit_path.write_text(json.dumps(audit), encoding="utf-8")
                elif mutation == "missing_figure":
                    Path(audit["artifacts"]["figures"][0]["path"]).unlink()
                else:
                    audit["artifacts"]["plot_script"]["sha256"] = "0" * 64
                    audit_path.write_text(json.dumps(audit), encoding="utf-8")
                with self.assertRaisesRegex(
                    screen.ScreenError, "artifacts are incomplete or stale"
                ) as caught:
                    screen.validate_phase_field_audit(
                        root, "case", hd_path, parameters, runner_state
                    )
                self.assertIsInstance(
                    caught.exception.__cause__, screen.phase_plot.PhaseAuditError
                )

    def test_phase_field_inputs_must_match_runner_artifacts_item_by_item(self) -> None:
        for mutation, message in (
            ("checkpoint", "checkpoint does not exactly bind"),
            ("history", "pipeline_history does not exactly bind"),
            ("frame", "result VTP inputs do not exactly bind"),
        ):
            with (
                self.subTest(mutation=mutation),
                tempfile.TemporaryDirectory() as temporary,
            ):
                root = Path(temporary)
                _audit_path, hd_path, parameters, runner_state = (
                    self._write_phase_field_audit(root)
                )
                if mutation == "checkpoint":
                    runner_state["stages"]["eq"]["audit"]["final_vtp"][
                        "sha256"
                    ] = "0" * 64
                elif mutation == "history":
                    runner_state["stages"]["hd"]["audit"]["pipeline_history"][
                        "sha256"
                    ] = "0" * 64
                else:
                    runner_state["stages"]["hd"]["audit"]["particle_vtp"][17][
                        "sha256"
                    ] = "0" * 64
                with self.assertRaisesRegex(screen.ScreenError, message):
                    screen.validate_phase_field_audit(
                        root, "case", hd_path, parameters, runner_state
                    )

    def test_phase_rotation_records_nontrivial_alignment_and_zero_residual(self) -> None:
        coefficients = np.zeros((3, 3, 2), dtype=float)
        coefficients[0] = 100.0
        coefficients[1, :, 0] = [10.0, 0.0, -10.0]
        coefficients[2, :, 0] = [0.0, 10.0, 0.0]
        coefficients[1, :, 1] = [5.0, 0.0, -5.0]
        coefficients[2, :, 1] = [0.0, 5.0, 0.0]
        rotated, audit = rotated_coefficients(
            coefficients,
            np.asarray([0, 0, 0], dtype=np.int64),
            minimum_reference_amplitude=1.0,
        )
        self.assertTrue(
            np.allclose(
                np.hypot(rotated[1], rotated[2]),
                np.hypot(coefficients[1], coefficients[2]),
            )
        )
        liquid = audit["phase_alignment"]["liquid"]
        self.assertEqual(liquid["phase_defined_particle_count"], 3)
        self.assertEqual(liquid["changed_particle_count"], 2)
        self.assertGreater(
            liquid["max_abs_source_to_surface_phase_difference_deg"], 0.0
        )
        self.assertLessEqual(
            liquid["max_abs_residual_to_surface_phase_difference_deg"], 1.0e-12
        )

    def test_integrates_exact_inclusive_window(self) -> None:
        rows = [
            row(
                time,
                lagged_force=2.0,
                erased_force=1.0,
                lagged_if_area=0.2,
                erased_if_area=0.1,
                lagged_joint=0.02,
                erased_joint=0.01,
            )
            for time in (1.3, 1.95, 2.6)
        ]
        result = screen.integrate_window(rows, 1.3, 2.6)
        self.assertEqual(result["saved_frames"], 3)
        self.assertAlmostEqual(
            result["local_positive_force_activity_impulse"]["lagged"], 2.6
        )
        self.assertAlmostEqual(
            result["local_positive_force_activity_impulse"][
                "lagged_minus_phase_erased_fraction"
            ],
            1.0,
        )

    def test_incomplete_window_is_rejected(self) -> None:
        rows = [
            row(
                1.95,
                lagged_force=2.0,
                erased_force=1.0,
                lagged_if_area=0.2,
                erased_if_area=0.1,
                lagged_joint=0.02,
                erased_joint=0.01,
            )
        ]
        with self.assertRaisesRegex(screen.ScreenError, "incomplete"):
            screen.integrate_window(rows, 1.3, 2.6)

    def test_gate_requires_both_cycles_direction_and_joint_continuity(self) -> None:
        rows = []
        for time in (1.3, 1.95, 2.6, 3.25, 3.9):
            rows.append(
                row(
                    time,
                    lagged_force=1.1,
                    erased_force=1.0,
                    lagged_if_area=0.11,
                    erased_if_area=0.1,
                    lagged_joint=0.011,
                    erased_joint=0.01,
                )
            )
        windows = [
            screen.integrate_window(rows, 1.3, 2.6),
            screen.integrate_window(rows, 2.6, 3.9),
            screen.integrate_window(rows, 1.3, 3.9),
        ]
        gate = screen.diagnostic_gate(
            rows,
            windows,
            minimum_resolved_area_m2=0.001,
            phase_qualification=self.qualified_phase(),
        )
        self.assertTrue(gate["passed"])
        broken = [dict(item) for item in rows]
        for item in broken:
            if item["time_s"] > 2.6:
                item["lagged_signed_force_n_per_m"] = 0.9
                item["lagged_net_uplift_force_n_per_m"] = 0.9
        broken_windows = [
            screen.integrate_window(broken, 1.3, 2.6),
            screen.integrate_window(broken, 2.6, 3.9),
            screen.integrate_window(broken, 1.3, 3.9),
        ]
        self.assertFalse(
            screen.diagnostic_gate(
                broken,
                broken_windows,
                minimum_resolved_area_m2=0.001,
                phase_qualification=self.qualified_phase(),
            )["passed"]
        )

    def test_gate_rejects_positive_hotspots_when_net_uplift_decreases(self) -> None:
        rows = []
        for time in (1.3, 1.95, 2.6, 3.25, 3.9):
            item = row(
                time,
                lagged_force=1.1,
                erased_force=1.0,
                lagged_if_area=0.11,
                erased_if_area=0.1,
                lagged_joint=0.011,
                erased_joint=0.01,
            )
            item["lagged_signed_force_n_per_m"] = -0.2
            item["phase_erased_signed_force_n_per_m"] = -0.1
            item["lagged_net_uplift_force_n_per_m"] = 0.0
            item["phase_erased_net_uplift_force_n_per_m"] = 0.0
            rows.append(item)
        windows = [
            screen.integrate_window(rows, 1.3, 2.6),
            screen.integrate_window(rows, 2.6, 3.9),
            screen.integrate_window(rows, 1.3, 3.9),
        ]
        gate = screen.diagnostic_gate(
            rows,
            windows,
            minimum_resolved_area_m2=0.001,
            phase_qualification=self.qualified_phase(),
        )
        self.assertFalse(gate["passed"])
        self.assertFalse(
            gate["checks"][
                "both_post_ramp_cycles_have_positive_signed_force_delta"
            ]
        )

    def test_post_ramp_selection_ignores_larger_ramp_transient(self) -> None:
        rows = [
            {"time_s": 0.65, "delta": 100.0},
            {"time_s": 1.3, "delta": 2.0},
            {"time_s": 2.6, "delta": 3.0},
            {"time_s": 3.9, "delta": 1.0},
        ]
        selected = screen.select_post_ramp_maximum(rows, "delta")
        self.assertEqual(selected["time_s"], 2.6)

    def test_phase_field_qualification_enforces_locked_crown_thresholds(self) -> None:
        audit = {
            "probe_results": {
                "surface": {"fundamental_amplitude_pa": 100.0},
                "crown": {
                    "fundamental_amplitude_pa": 10.0,
                    "harmonic_r_squared": 0.95,
                    "phase_difference_from_same_column_surface_deg": -20.0,
                },
            }
        }
        self.assertTrue(screen.phase_field_qualification(audit)["passed"])
        for field, value in (
            ("phase_difference_from_same_column_surface_deg", 17.9),
            ("fundamental_amplitude_pa", 4.9),
            ("harmonic_r_squared", 0.79),
        ):
            changed = copy.deepcopy(audit)
            changed["probe_results"]["crown"][field] = value
            with self.subTest(field=field):
                self.assertFalse(
                    screen.phase_field_qualification(changed)["passed"]
                )

    def test_masked_crown_phase_is_retained_as_a_qualification_failure(self) -> None:
        audit = {
            "probe_results": {
                "surface": {"fundamental_amplitude_pa": 100.0},
                "crown": {
                    "fundamental_amplitude_pa": 2.0,
                    "harmonic_r_squared": 0.7,
                    "phase_difference_from_same_column_surface_deg": None,
                },
            }
        }
        qualification = screen.phase_field_qualification(audit)
        self.assertFalse(qualification["passed"])
        self.assertIsNone(
            qualification["observed"]["crown_same_column_phase_lag_deg"]
        )
        self.assertFalse(
            qualification["checks"]["crown_abs_same_column_phase_lag_ge_18deg"]
        )

    def test_phase_database_binding_rejects_substituted_source(self) -> None:
        metadata = {
            "source": {
                "points": "/expected/source_points.txt",
                "points_sha256": "a" * 64,
                "values": "/expected/source_values.bin",
                "values_sha256": "b" * 64,
            },
            "output": {
                "points": "/expected/output_points.txt",
                "points_sha256": "c" * 64,
                "values": "/expected/output_values.bin",
                "values_sha256": "d" * 64,
            },
            "qa": {
                "phase_alignment": {
                    "liquid": {
                        "phase_defined_particle_count": 10,
                        "changed_particle_count": 8,
                        "max_abs_source_to_surface_phase_difference_deg": 20.0,
                        "max_abs_residual_to_surface_phase_difference_deg": 0.0,
                    }
                }
            },
        }
        lagged = {
            "points": "/substituted/source_points.txt",
            "points_sha256": "a" * 64,
            "values": "/expected/source_values.bin",
            "values_sha256": "b" * 64,
            "step_interval": 130,
            "max_step": 39000,
            "source_dt_s": 1.0e-4,
        }
        erased = {
            "points": "/expected/output_points.txt",
            "points_sha256": "c" * 64,
            "values": "/expected/output_values.bin",
            "values_sha256": "d" * 64,
            "step_interval": 130,
            "max_step": 39000,
            "source_dt_s": 1.0e-4,
        }
        phase_audit = {"metadata": copy.deepcopy(metadata)}
        with self.assertRaisesRegex(screen.ScreenError, "not the compared database"):
            screen.validate_phase_database_bindings(
                Path("/"), phase_audit, lagged, erased
            )


if __name__ == "__main__":
    unittest.main()
