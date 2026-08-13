from __future__ import annotations

import math
import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np


CASE_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(CASE_DIR))

import analyze_study as analysis  # noqa: E402


class StudyAnalysisTest(unittest.TestCase):
    def test_harmonic_delay_convention(self):
        times = np.linspace(0.0, 4.0, 161)
        delay = math.radians(27.0)
        values = 12.0 + 7.0 * np.cos(2.0 * math.pi * times - delay)
        mean, amplitude, phase = analysis.harmonic_fit(times, values, 1.0)
        self.assertAlmostEqual(mean, 12.0, places=10)
        self.assertAlmostEqual(amplitude, 7.0, places=10)
        self.assertAlmostEqual(analysis.wrap_degrees(phase), 27.0, places=10)

    def test_stress_invariants_supports_both_compression_signs(self):
        negative = np.asarray([-3000.0, -3000.0, -3000.0, 0.0, 0.0, 0.0])
        positive = -negative
        mean, deviator = analysis.stress_invariants(negative)
        self.assertAlmostEqual(mean, 3000.0)
        self.assertAlmostEqual(deviator, 0.0)
        mean, deviator = analysis.stress_invariants(positive, compression_sign=1.0)
        self.assertAlmostEqual(mean, 3000.0)
        self.assertAlmostEqual(deviator, 0.0)

    def test_vertical_stress_loss_sign_reference_and_threshold(self):
        reference = np.asarray([-1000.0, 1000.0, 99.0, -1000.0])
        current = np.asarray([-50.0, 50.0, 0.0, 10.0])
        ratios, eligible, lost = analysis.vertical_stress_loss_state(
            current, reference
        )
        np.testing.assert_array_equal(eligible, [True, True, False, True])
        self.assertAlmostEqual(ratios[0], 0.05)
        self.assertAlmostEqual(ratios[1], 0.05)
        self.assertTrue(np.isnan(ratios[2]))
        self.assertAlmostEqual(ratios[3], -0.01)
        np.testing.assert_array_equal(lost, [True, True, False, True])

    def test_upward_excess_seepage_force_removes_hydrostatic_and_clips(self):
        forces = np.asarray(
            [[500.0, 10_000.0], [0.0, 11_000.0], [0.0, 9_000.0]]
        )
        index = analysis.upward_excess_seepage_force_index(
            forces,
            np.asarray([1000.0, 1000.0, 1000.0]),
            np.asarray([1000.0, 1000.0, 1000.0]),
            np.asarray([0.0, -10.0]),
        )
        np.testing.assert_allclose(index, [0.0, 1.0, 0.0])
        with self.assertRaisesRegex(ValueError, "densities"):
            analysis.upward_excess_seepage_force_index(
                forces,
                np.asarray([1000.0, 0.0, 1000.0]),
                np.asarray([1000.0, 1000.0, 1000.0]),
                np.asarray([0.0, -10.0]),
            )

    def test_pipeline_large_deformation_gate(self):
        history = [
            {
                "time": float(i),
                "displacement_y": 0.03 * i,
                "displacement_x": 0.0,
                "angle": 0.0,
                "contacts": 1.0,
                "contact_force_x": 0.0,
                "contact_force_y": 1.0,
            }
            for i in range(4)
        ]
        metrics, _ = analysis.pipeline_metrics(history, 0.12, 1.0, 1.0, 0.03)
        self.assertTrue(metrics["large_deformation_demonstrated"])
        self.assertTrue(metrics["breakout_reached"])
        self.assertAlmostEqual(
            metrics["maximum_abs_vertical_displacement_over_D"], 0.75
        )

    def test_nominal_reference_area_uses_initial_lattice(self):
        coordinates = np.asarray(
            [
                [0.05, 0.05],
                [0.05, 0.15],
                [0.15, 0.05],
                [0.25, 0.05],
                [0.25, 0.15],
            ]
        )
        self.assertAlmostEqual(
            analysis.nominal_reference_particle_area(coordinates), 0.01
        )

    def test_configured_fine_particle_location_is_resolved(self):
        config = {
            "particles": [
                {"generator": {"location": "mesh_fine/particles.txt"}}
            ]
        }
        case_root = Path("case_root")
        self.assertEqual(
            analysis.resolve_particle_input_path(config, case_root),
            case_root / "mesh_fine" / "particles.txt",
        )

    def test_pre_registered_pipeline_support_roi(self):
        config = {
            "analysis": {
                "rigid_pipeline": {"center": [0.7, 0.41], "outer_radius": 0.06}
            },
            "materials": [{}, {"sea_level": 1.0, "depth_left": 0.5}],
        }
        coordinates = np.asarray(
            [[0.7, 0.3], [0.51, 0.3], [0.7, 0.25], [0.7, 0.51]]
        )
        np.testing.assert_array_equal(
            analysis.pipeline_support_roi(config, coordinates),
            [True, False, False, False],
        )

    def test_displacement_mesh_crossing_validates_shape_and_finite_values(self):
        metrics = analysis.displacement_mesh_crossing_metrics(
            np.asarray([[0.06, 0.08], [0.2, 0.0]]), 2, 0.1
        )
        self.assertAlmostEqual(metrics["max_abs_displacement_over_h"], 2.0)
        self.assertAlmostEqual(metrics["fraction_abs_displacement_ge_h"], 1.0)
        self.assertTrue(metrics["crossed_one_background_cell"])
        with self.assertRaisesRegex(ValueError, "shape"):
            analysis.displacement_mesh_crossing_metrics(
                np.asarray([[0.1], [0.2]]), 2, 0.1
            )
        with self.assertRaisesRegex(ValueError, "NaN/Inf"):
            analysis.displacement_mesh_crossing_metrics(
                np.asarray([[np.nan, 0.0], [0.2, 0.0]]), 2, 0.1
            )

    def test_threshold_spatiotemporal_metrics(self):
        times = np.asarray([0.0, 1.0, 3.0, 4.0])
        masks = np.asarray(
            [
                [False, False, False, False],
                [True, False, False, False],
                [True, True, False, False],
                [False, False, False, False],
            ]
        )
        metrics = analysis.threshold_spatiotemporal_metrics(
            times, masks, particle_area=0.02
        )
        self.assertAlmostEqual(metrics["max_fraction"], 0.5)
        self.assertAlmostEqual(metrics["max_nominal_reference_area_m2"], 0.04)
        self.assertAlmostEqual(
            metrics["time_integrated_nominal_reference_area_m2_s"], 0.09
        )
        self.assertAlmostEqual(
            metrics["resolved_duration_anywhere_lower_bound_s"], 2.0
        )

    def test_release_segments_and_particle_level_trigger_ordering(self):
        times = np.asarray([0.0, 1.0, 2.0, 3.0])
        seepage = np.asarray(
            [
                [True, False, False, False],
                [False, True, False, False],
                [False, False, True, False],
                [False, False, False, False],
            ]
        )
        stress = np.asarray(
            [
                [False, False, False, False],
                [True, True, True, False],
                [True, True, True, True],
                [False, False, False, False],
            ]
        )
        timing = analysis.trigger_to_stress_loss_timing(
            times, seepage, stress, wave_period=1.0
        )
        self.assertEqual(timing["stress_loss_particles"], 4)
        self.assertEqual(timing["stress_loss_particles_IF_first"], 1)
        self.assertEqual(timing["stress_loss_particles_simultaneous_first"], 1)
        self.assertEqual(timing["stress_loss_particles_stress_first"], 1)
        self.assertEqual(timing["stress_loss_particles_without_IF"], 1)
        self.assertAlmostEqual(
            timing[
                "fraction_stress_loss_particles_with_IF_in_preceding_cycle"
            ],
            0.25,
        )

        segments = analysis.threshold_release_segments(
            times,
            seepage & stress,
            particle_area=0.02,
            release_time=1.5,
        )
        self.assertEqual(
            segments["release_segment_area_basis"],
            "nominal_reference_particle_area_fallback",
        )
        self.assertAlmostEqual(
            segments["pre_release_time_integrated_area_m2_s"], 0.02
        )
        self.assertAlmostEqual(
            segments["post_release_time_integrated_area_m2_s"], 0.02
        )

    def test_completed_result_requires_exact_grid_and_matching_hashes(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            config_path = root / "04_RL.json"
            config = {
                "analysis": {"uuid": "CASE", "nsteps": 20, "dt": 0.1},
                "post_processing": {"output_steps": 10},
            }
            config_path.write_text(json.dumps(config), encoding="utf-8")
            result = root / "results" / "CASE"
            result.mkdir(parents=True)
            files = [result / "particle10.vtp", result / "particle20.vtp"]
            for path in files:
                path.write_text(f"frame-{path.stem}", encoding="utf-8")
            final = files[-1]
            sentinel = {
                "schema": analysis.COMPLETION_SCHEMA,
                "config": {
                    "path": config_path.name,
                    "sha256": analysis.file_sha256(config_path),
                    "uuid": "CASE",
                    "nsteps": 20,
                    "validation_profile": "registered-study-v1",
                },
                "artifacts": {
                    "final_vtp": {
                        "path": str(final.relative_to(root)),
                        "size_bytes": final.stat().st_size,
                        "sha256": analysis.file_sha256(final),
                    }
                },
                "runtime_dependencies": {},
            }
            (result / analysis.COMPLETION_FILENAME).write_text(
                json.dumps(sentinel), encoding="utf-8"
            )
            audit = analysis.validate_completed_result(
                config_path, root, config, result, files
            )
            self.assertEqual(audit["expected_particle_steps"], [10, 20])
            self.assertEqual(audit["expected_particle_times_s"], [1.0, 2.0])
            sentinel["config"]["validation_profile"] = "pipeline-local-smoke-v1"
            (result / analysis.COMPLETION_FILENAME).write_text(
                json.dumps(sentinel), encoding="utf-8"
            )
            with self.assertRaisesRegex(ValueError, "not registered-study"):
                analysis.validate_completed_result(
                    config_path, root, config, result, files
                )
            sentinel["config"]["validation_profile"] = "registered-study-v1"
            (result / analysis.COMPLETION_FILENAME).write_text(
                json.dumps(sentinel), encoding="utf-8"
            )
            with self.assertRaisesRegex(ValueError, "time grid"):
                analysis.validate_completed_result(
                    config_path, root, config, result, files[-1:]
                )
            final.write_text("corrupt", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "size|hash"):
                analysis.validate_completed_result(
                    config_path, root, config, result, files
                )

    def test_RL_RE_comparison_rejects_different_time_grids(self):
        def summary(code: str) -> dict[str, object]:
            return {
                "case": code,
                "config": f"/case/04_{code}.json",
                "analysis_dt_s": 0.1,
                "particles": {"particle_count": 4},
                "completion_audit": {
                    "config_sha256": code,
                    "expected_particle_steps": [10, 20],
                    "expected_particle_times_s": [1.0, 2.0],
                    "runtime_dependencies": {
                        "resume_equilibrium": {
                            "checkpoint_hdf5": {"sha256": "h5"},
                            "qa_vtp": {"sha256": "vtp"},
                        },
                        "read_pressure_database": {
                            "header": {
                                "format_version": "V2",
                                "particle_count": 4,
                                "step_interval": 10,
                                "max_step": 20,
                                "source_dt_s": 0.1,
                                "frame_count": 3,
                            }
                        },
                    },
                },
            }

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            for code, prefix in (("RL", "pressure"), ("RE", "phase_erased")):
                config_path = root / f"04_{code}.json"
                config_path.write_text(
                    json.dumps(
                        {
                            "title": code,
                            "analysis": {
                                "uuid": code,
                                "dt": 0.1,
                                "prescribed_phase_pressures": {
                                    "path": f"database/{code}",
                                    "file_prefix": prefix,
                                    "step_interval": 10,
                                },
                            },
                        }
                    ),
                    encoding="utf-8",
                )
            rl = summary("RL")
            re_case = summary("RE")
            rl["config"] = str(root / "04_RL.json")
            re_case["config"] = str(root / "04_RE.json")
            audit = analysis.validate_RL_RE_comparison([rl, re_case])
            self.assertEqual(audit["particle_steps"], [10, 20])
            re_case["completion_audit"]["expected_particle_steps"] = [10]
            with self.assertRaisesRegex(ValueError, "particle_steps"):
                analysis.validate_RL_RE_comparison([rl, re_case])

    @staticmethod
    def _frame(
        coordinates: np.ndarray,
        order: np.ndarray,
        reference: np.ndarray,
        vertical: np.ndarray,
        seepage_index: np.ndarray,
        volumes: np.ndarray,
        ru: np.ndarray,
        displacements: np.ndarray,
    ) -> tuple[np.ndarray, dict[str, np.ndarray]]:
        stresses = np.zeros((len(reference), 6), dtype=np.float64)
        stresses[:, 0] = vertical
        stresses[:, 1] = vertical
        stresses[:, 2] = vertical
        densities = np.full(len(reference), 1000.0)
        gamma_sub = np.full(len(reference), 1000.0)
        forces = np.zeros((len(reference), 3), dtype=np.float64)
        forces[:, 1] = 10_000.0 + 1000.0 * seepage_index
        arrays = {
            "ids": np.arange(len(reference), dtype=np.float64),
            "volumes": volumes,
            "PIC_ru": ru,
            "PIC_pore_pressure_excess": 100.0 * ru,
            "stresses": stresses,
            "initial_vertical_effective_stresses": reference,
            "vertical_effective_stress_remaining_ratios": vertical / reference,
            "liquid_seepage_forces": forces,
            "liquid_densities": densities,
            "gamma_sub": gamma_sub,
            "liquefaction_potentials": np.maximum(0.0, 1.0 - vertical / reference),
            "momentary_liquefied": (vertical / reference <= 0.05).astype(float),
            "displacements": displacements,
        }
        return coordinates[order], {
            name: values[order] for name, values in arrays.items()
        }

    def test_particle_histories_uses_two_primary_criteria_and_ids(self):
        coordinates = np.asarray(
            [[0.05, 0.05], [0.05, 0.15], [0.15, 0.05], [0.15, 0.15]]
        )
        reference = np.asarray([-1000.0, 1000.0, -50.0, 2000.0])
        frames = [
            self._frame(
                coordinates,
                np.asarray([2, 0, 3, 1]),
                reference,
                reference,
                np.asarray([0.0, 0.5, 1.0, 0.0]),
                np.asarray([0.009, 0.011, 0.01, 0.01]),
                np.asarray([10.0, 0.0, 0.0, 0.0]),
                np.zeros((4, 2)),
            ),
            self._frame(
                coordinates,
                np.asarray([1, 3, 0, 2]),
                reference,
                np.asarray([-50.0, 50.0, -1.0, 101.0]),
                np.asarray([1.2, 0.0, 2.0, 0.0]),
                np.asarray([0.012, 0.01, 0.01, 0.01]),
                np.asarray([20.0, 20.0, 0.0, 0.0]),
                np.asarray([[0.05, 0.0], [0.1, 0.0], [0.2, 0.0], [0.0, 0.0]]),
            ),
            self._frame(
                coordinates,
                np.asarray([3, 2, 1, 0]),
                reference,
                np.asarray([-40.0, -10.0, -1.0, 80.0]),
                np.asarray([0.0, 1.0, 1.5, 2.0]),
                np.asarray([0.013, 0.009, 0.01, 0.01]),
                np.asarray([0.0, 0.0, 0.0, 0.0]),
                np.asarray([[0.3, 0.0], [0.05, 0.0], [0.1, 0.0], [0.2, 0.0]]),
            ),
        ]
        config = {
            "analysis": {
                "dt": 0.1,
                "rigid_pipeline": {
                    "release_time": 1.0,
                    "center": [0.05, 0.1],
                    "outer_radius": 0.025,
                },
            },
            "mesh": {"cellsize_min": 0.1},
            "materials": [
                {"type": "SANISAND2D"},
                {"sea_level": 0.15, "depth_left": 0.0, "wave_period": 1.0},
            ],
            "external_loading_conditions": {"gravity": [0.0, -10.0]},
        }
        paths = [Path(f"particle{step:07d}.vtp") for step in (0, 10, 20)]
        with patch.object(analysis, "read_vtp", side_effect=frames):
            rows, summary = analysis.particle_histories(
                paths,
                config,
                coordinates,
                probes=[analysis.Probe("surface", 0, 0.05, 0.05, 0.05, 0.05)],
            )

        stress_name = "stress_loss_Rsigma_le_0p05"
        seepage_name = "upward_seepage_IF_ge_1"
        joint_name = "joint_IF_ge_1_and_Rsigma_le_0p05"
        self.assertAlmostEqual(rows[1][f"current_area_{stress_name}_m2"], 0.022)
        self.assertAlmostEqual(rows[1][f"current_area_{seepage_name}_m2"], 0.022)
        self.assertAlmostEqual(rows[1][f"current_area_{joint_name}_m2"], 0.012)
        self.assertAlmostEqual(rows[2][f"current_area_{joint_name}_m2"], 0.019)
        self.assertAlmostEqual(
            summary[f"{stress_name}_time_integrated_current_area_m2_s"], 0.038
        )
        self.assertAlmostEqual(
            summary[f"{seepage_name}_time_integrated_current_area_m2_s"], 0.0415
        )
        self.assertAlmostEqual(
            summary[f"{stress_name}_resolved_duration_anywhere_lower_bound_s"],
            1.0,
        )
        self.assertAlmostEqual(
            summary[f"{seepage_name}_resolved_duration_anywhere_lower_bound_s"],
            2.0,
        )
        self.assertAlmostEqual(
            summary[f"{joint_name}_pre_release_time_integrated_area_m2_s"],
            0.006,
        )
        self.assertAlmostEqual(
            summary[f"{joint_name}_post_release_time_integrated_area_m2_s"],
            0.0155,
        )
        self.assertEqual(summary["stress_loss_particles"], 3)
        self.assertEqual(summary["stress_loss_particles_simultaneous_first"], 2)
        self.assertEqual(
            summary["stress_reference_basis"],
            "checkpoint_initial_vertical_effective_stresses",
        )
        self.assertEqual(summary["stress_reference_time_s"], 0.0)
        self.assertEqual(summary["criterion_area_basis"], "current_particle_volumes")
        self.assertIn(
            "ru has no area metric",
            summary["liquefaction_area_basis_compatibility_note"],
        )
        self.assertEqual(
            summary["particle_id_basis"],
            "validated_unique_ids_reordered_to_initial_0_to_N_minus_1",
        )
        self.assertAlmostEqual(summary["solver_Rsigma_QA_max_abs_difference"], 0.0)
        self.assertEqual(summary["max_ru_diagnostic"], 20.0)
        self.assertAlmostEqual(
            rows[1]["surface_pressure_excess_pa"], 2000.0
        )
        self.assertNotIn("soil_max_abs_displacement_over_h", rows[0])
        self.assertAlmostEqual(rows[1]["soil_max_abs_displacement_over_h"], 2.0)
        self.assertAlmostEqual(
            rows[1]["soil_fraction_abs_displacement_ge_h"], 0.5
        )
        self.assertAlmostEqual(
            rows[1]["support_ROI_soil_max_abs_displacement_over_h"], 1.0
        )
        self.assertAlmostEqual(
            rows[1]["support_ROI_soil_fraction_abs_displacement_ge_h"], 0.5
        )
        self.assertEqual(summary["post_release_displacement_saved_frames"], 2)
        self.assertAlmostEqual(
            summary["post_release_max_soil_abs_displacement_over_h"], 3.0
        )
        self.assertAlmostEqual(
            summary["post_release_max_soil_fraction_abs_displacement_ge_h"],
            0.75,
        )
        self.assertTrue(summary["post_release_soil_crossed_one_background_cell"])
        self.assertAlmostEqual(
            summary[
                "post_release_max_support_ROI_soil_abs_displacement_over_h"
            ],
            3.0,
        )
        self.assertAlmostEqual(
            summary[
                "post_release_max_support_ROI_soil_fraction_abs_displacement_ge_h"
            ],
            0.5,
        )
        self.assertTrue(
            summary[
                "post_release_support_ROI_soil_crossed_one_background_cell"
            ]
        )
        self.assertIn(
            "not a strain", summary["mesh_crossing_displacement_interpretation"]
        )
        forbidden = ("ru_ge", "area_ru", "momentary_liquefied_flag_time_")
        for key in [*summary, *(key for row in rows for key in row)]:
            self.assertFalse(any(token in key for token in forbidden), key)

    def test_first_frame_reference_and_order_are_explicit_legacy_fallbacks(self):
        coordinates = np.asarray(
            [[0.05, 0.05], [0.05, 0.15], [0.15, 0.05], [0.15, 0.15]]
        )
        stress = np.zeros((4, 6))
        stress[:, 1] = [-1000.0, -1000.0, -1000.0, -1000.0]
        config = {"analysis": {"dt": 0.1}, "materials": [{"type": "SANISAND2D"}]}
        with patch.object(
            analysis,
            "read_vtp",
            return_value=(coordinates.copy(), {"stresses": stress}),
        ):
            _, summary = analysis.particle_histories(
                [Path("particle0000065.vtp")], config, coordinates, probes=[]
            )
        self.assertEqual(
            summary["stress_reference_basis"],
            "first_saved_dynamic_stress_fallback",
        )
        self.assertTrue(summary["stress_reference_fallback_used"])
        self.assertAlmostEqual(summary["stress_reference_time_s"], 6.5)
        self.assertEqual(summary["particle_id_basis"], "VTK_order_fallback_no_ids")

    def test_checkpoint_reference_time_is_dynamic_zero_not_first_saved_time(self):
        coordinates = np.asarray(
            [[0.05, 0.05], [0.05, 0.15], [0.15, 0.05], [0.15, 0.15]]
        )
        stress = np.zeros((4, 6))
        stress[:, 1] = -1000.0
        config = {"analysis": {"dt": 0.001}, "materials": [{"type": "SANISAND2D"}]}
        with patch.object(
            analysis,
            "read_vtp",
            return_value=(
                coordinates.copy(),
                {
                    "stresses": stress,
                    "initial_vertical_effective_stresses": np.full(4, -1000.0),
                },
            ),
        ):
            rows, summary = analysis.particle_histories(
                [Path("particle0000065.vtp")], config, coordinates, probes=[]
            )
        self.assertAlmostEqual(rows[0]["time_s"], 0.065)
        self.assertEqual(summary["stress_reference_time_s"], 0.0)

    def test_initial_vertical_reference_must_remain_constant(self):
        coordinates = np.asarray(
            [[0.05, 0.05], [0.05, 0.15], [0.15, 0.05], [0.15, 0.15]]
        )
        stress = np.zeros((4, 6))
        stress[:, 1] = -1000.0
        first = np.full(4, -1000.0)
        changed = first.copy()
        changed[2] = -900.0
        frames = [
            (
                coordinates.copy(),
                {"stresses": stress, "initial_vertical_effective_stresses": first},
            ),
            (
                coordinates.copy(),
                {"stresses": stress, "initial_vertical_effective_stresses": changed},
            ),
        ]
        config = {"analysis": {"dt": 0.1}, "materials": [{"type": "SANISAND2D"}]}
        with patch.object(analysis, "read_vtp", side_effect=frames):
            with self.assertRaisesRegex(ValueError, "changed between"):
                analysis.particle_histories(
                    [Path("particle0000000.vtp"), Path("particle0000001.vtp")],
                    config,
                    coordinates,
                    probes=[],
                )

    def test_particle_ids_must_be_unique(self):
        coordinates = np.asarray(
            [[0.05, 0.05], [0.05, 0.15], [0.15, 0.05], [0.15, 0.15]]
        )
        config = {"analysis": {"dt": 0.1}, "materials": [{"type": "SANISAND2D"}]}
        with patch.object(
            analysis,
            "read_vtp",
            return_value=(coordinates.copy(), {"ids": np.asarray([0, 1, 1, 3])}),
        ):
            with self.assertRaisesRegex(ValueError, "not unique"):
                analysis.particle_histories(
                    [Path("particle0000000.vtp")], config, coordinates, probes=[]
                )

    def test_particle_histories_rejects_invalid_current_volumes(self):
        coordinates = np.asarray(
            [[0.05, 0.05], [0.05, 0.15], [0.15, 0.05], [0.15, 0.15]]
        )
        config = {"analysis": {"dt": 0.1}, "materials": [{"type": "SANISAND2D"}]}
        path = Path("particle0000000.vtp")
        invalid = (
            (np.asarray([0.01, 0.01, 0.01]), "volume count"),
            (np.asarray([0.01, np.nan, 0.01, 0.01]), "NaN/Inf"),
            (np.asarray([0.01, 0.0, 0.01, 0.01]), "positive"),
        )
        for volumes, message in invalid:
            with self.subTest(message=message), patch.object(
                analysis,
                "read_vtp",
                return_value=(coordinates.copy(), {"volumes": volumes}),
            ):
                with self.assertRaisesRegex(ValueError, message):
                    analysis.particle_histories(
                        [path], config, coordinates, probes=[]
                    )

    def test_particle_histories_rejects_partial_displacement_frames(self):
        coordinates = np.asarray(
            [[0.05, 0.05], [0.05, 0.15], [0.15, 0.05], [0.15, 0.15]]
        )
        frames = [
            (coordinates.copy(), {"displacements": np.zeros((4, 2))}),
            (coordinates.copy(), {}),
        ]
        config = {
            "analysis": {"dt": 0.1},
            "mesh": {"cellsize_min": 0.1},
            "materials": [{"type": "SANISAND2D"}],
        }
        with patch.object(analysis, "read_vtp", side_effect=frames):
            with self.assertRaisesRegex(ValueError, "missing from some"):
                analysis.particle_histories(
                    [Path("particle0000000.vtp"), Path("particle0000001.vtp")],
                    config,
                    coordinates,
                    probes=[],
                )


if __name__ == "__main__":
    unittest.main()
