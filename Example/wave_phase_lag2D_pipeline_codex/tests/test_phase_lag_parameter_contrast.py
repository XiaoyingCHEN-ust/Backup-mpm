from __future__ import annotations

import unittest

import numpy as np

import analyze_phase_lag_parameter_contrast as contrast
import summarize_phase_lag_permeability_screen as summary


def _metric(value: float) -> dict[str, float]:
    return {"lagged": value}


def _case(local: tuple[float, float], signed: tuple[float, float]) -> dict:
    branches = []
    raw = []
    local_linear = []
    for index in range(2):
        branches.append(
            {
                "start_s": 2.34 + 1.3 * index,
                "end_s": 2.60 + 1.3 * index,
                "saved_frames": 5,
                "local_upward_activity_impulse": _metric(local[index]),
                "signed_support_force_impulse": _metric(signed[index]),
            }
        )
        raw.append(
            {
                "horizontal_absolute_force_activity_impulse": _metric(100.0),
                "vertical_absolute_force_activity_impulse": _metric(120.0),
            }
        )
        local_linear.append(
            {
                "horizontal_absolute_force_activity_impulse": _metric(90.0),
                "vertical_absolute_force_activity_impulse": _metric(135.0),
            }
        )
    return {
        "branches": branches,
        "raw_directional_branches": raw,
        "secondary_local_linear_branches": local_linear,
    }


def _pair_case(
    permeability: float, saturation: float, *, surface_shift: float = 0.0
) -> dict:
    rows = [
        {
            "time_s": 0.065 * index,
            "surface_excess_pressure_pa": float(index) + surface_shift,
        }
        for index in range(4)
    ]
    return {
        "label": f"k{permeability}_sw{saturation}",
        "validated_case": {
            "intrinsic_permeability_m2": permeability,
            "liquid_saturation": saturation,
            "pressure_smoothing": False,
            "pipeline_fixed": True,
        },
        "reference": np.asarray([[0.0, 0.0], [1.0, 0.0]]),
        "rows": rows,
        "branch_rows": [rows[:3], rows[1:]],
        "hd": {
            "materials": [
                {"intrinsic_permeability": permeability},
                {
                    "liquid_saturation": saturation,
                    "gas_saturation": 1.0 - saturation,
                    "wave_pressure": True,
                    "wave_height_ini_": 0.12,
                    "wave_period": 1.3,
                    "wave_ramp_time": 1.3,
                    "domain_length_x": 1.5,
                    "depth_left": 0.5,
                    "depth_right": 0.5,
                },
            ],
            "analysis": {
                "type": "ThermoMPMExplicitThreePhaseLag2D",
                "dt": 1.0e-4,
                "nsteps": 39000,
                "APIC": True,
                "PIC": 0.0,
                "pressure_smoothing": False,
                "pressure_smoothing_in_loop": False,
                "rigid_pipeline": {
                    "fixed": True,
                    "center": [0.7, 0.41],
                    "outer_radius": 0.06,
                },
                "uuid": f"CASE_{permeability}_{saturation}_HD",
                "resume": {"uuid": f"CASE_{permeability}_{saturation}_EQ"},
                "prescribed_phase_pressures": {
                    "path": f"pressure/{permeability}_{saturation}"
                },
            },
            "title": f"case {permeability} {saturation}",
            "external_loading_conditions": {"gravity": [0.0, -9.81]},
            "post_processing": {
                "output_steps": 650,
                "path": f"results/{permeability}_{saturation}",
            },
        },
    }


class PhaseLagParameterContrastTest(unittest.TestCase):
    def test_cross_saturation_loader_does_not_relax_registered_default(self) -> None:
        audit = {
            "parameters": {
                "intrinsic_permeability_m2": 9.79e-12,
                "liquid_saturation": 0.993,
                "gas_saturation": 0.007,
                "pressure_smoothing": False,
                "pressure_smoothing_in_loop": False,
                "pipeline_fixed": True,
            }
        }
        with self.assertRaises(summary.SummaryError):
            summary._validated_parameters(audit, "registered")
        observed = summary._validated_parameters(
            audit, "cross-parameter", expected_saturation=None
        )
        self.assertEqual(observed["liquid_saturation"], 0.993)

    def test_parameter_pair_requires_lower_k_lower_sw_and_same_driver(self) -> None:
        strong = _pair_case(1.0e-13, 0.94)
        weak = _pair_case(9.79e-12, 0.993)
        result = contrast.validate_parameter_pair(strong, weak)
        self.assertEqual(result["maximum_surface_response_difference_pa"], 0.0)
        with self.assertRaisesRegex(
            contrast.ParameterContrastError, "lower intrinsic permeability"
        ):
            contrast.validate_parameter_pair(weak, strong)
        shifted = _pair_case(9.79e-12, 0.993, surface_shift=1.0e-2)
        with self.assertRaisesRegex(
            contrast.ParameterContrastError, "surface-pressure responses"
        ):
            contrast.validate_parameter_pair(strong, shifted)
        altered = _pair_case(9.79e-12, 0.993)
        altered["hd"]["materials"][0]["G0"] = 126.0
        with self.assertRaisesRegex(
            contrast.ParameterContrastError, "differ beyond permeability"
        ):
            contrast.validate_parameter_pair(strong, altered)

    def test_relative_fraction_handles_non_positive_reference(self) -> None:
        self.assertAlmostEqual(contrast.relative_fraction(120.0, 100.0), 0.2)
        self.assertIsNone(contrast.relative_fraction(1.0, 0.0))
        with self.assertRaises(contrast.ParameterContrastError):
            contrast.relative_fraction(float("nan"), 1.0)

    def test_display_vector_clipping_preserves_direction(self) -> None:
        vectors = np.asarray([[3.0, 4.0], [-6.0, 8.0], [0.5, 0.0]])
        clipped, mask = contrast.clip_vector_magnitudes(vectors, 5.0)
        np.testing.assert_allclose(np.linalg.norm(clipped, axis=1), [5.0, 5.0, 0.5])
        np.testing.assert_array_equal(mask, [False, True, False])
        self.assertGreater(float(np.dot(vectors[1], clipped[1])), 0.0)
        self.assertAlmostEqual(clipped[1, 0] / clipped[1, 1], -0.75)
        with self.assertRaises(contrast.ParameterContrastError):
            contrast.clip_vector_magnitudes(vectors, 0.0)

    def test_compare_recoveries_reports_direct_orientation(self) -> None:
        strong = _case((120.0, 130.0), (60.0, 70.0))
        weak = _case((100.0, 100.0), (50.0, 50.0))
        weak["raw_directional_branches"][1][
            "vertical_absolute_force_activity_impulse"
        ] = _metric(100.0)
        weak["secondary_local_linear_branches"][1][
            "vertical_absolute_force_activity_impulse"
        ] = _metric(100.0)
        rows = contrast.compare_recoveries(strong, weak)
        self.assertEqual(len(rows), 2)
        self.assertAlmostEqual(
            rows[1]["local_upward_activity"]["lower_minus_higher_fraction"], 0.3
        )
        self.assertAlmostEqual(
            rows[1]["raw_vertical_absolute_force_activity"][
                "lower_minus_higher_fraction"
            ],
            0.2,
        )
        self.assertGreater(
            rows[1]["one_layer_vertical_to_horizontal_ratio"][
                "lower_k_lower_sw"
            ],
            rows[1]["one_layer_vertical_to_horizontal_ratio"][
                "higher_k_higher_sw"
            ],
        )

    def test_compare_recoveries_rejects_different_windows(self) -> None:
        strong = _case((120.0, 130.0), (60.0, 70.0))
        weak = _case((100.0, 100.0), (50.0, 50.0))
        weak["branches"][0]["start_s"] += 0.01
        with self.assertRaisesRegex(
            contrast.ParameterContrastError, "branch metadata"
        ):
            contrast.compare_recoveries(strong, weak)


if __name__ == "__main__":
    unittest.main()
