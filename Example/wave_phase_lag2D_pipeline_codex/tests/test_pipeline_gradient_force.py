from __future__ import annotations

import json
import math
import unittest

import numpy as np

import plot_pipeline_gradient_force as bridge


class PipelineGradientForceTests(unittest.TestCase):
    def test_relative_difference_preserves_positive_denominator_results(self) -> None:
        self.assertAlmostEqual(bridge.relative_difference(12.0, 10.0), 0.2)
        self.assertAlmostEqual(
            bridge.relative_difference(12.0, 10.0, subtract_one=False), 1.2
        )

    def test_relative_difference_uses_null_for_nonpositive_denominator(self) -> None:
        self.assertIsNone(bridge.relative_difference(1.0, 0.0))
        self.assertIsNone(bridge.relative_difference(1.0, -1.0))
        self.assertEqual(
            json.dumps(
                {"relative_difference": bridge.relative_difference(1.0, 0.0)},
                allow_nan=False,
            ),
            '{"relative_difference": null}',
        )

    def test_relative_difference_rejects_nonfinite_values(self) -> None:
        for numerator, denominator in (
            (math.nan, 1.0),
            (math.inf, 1.0),
            (1.0, math.nan),
            (1.0, math.inf),
        ):
            with self.subTest(numerator=numerator, denominator=denominator):
                with self.assertRaisesRegex(ValueError, "must be finite"):
                    bridge.relative_difference(numerator, denominator)

        with self.assertRaisesRegex(ValueError, "must be finite"):
            bridge.relative_difference(np.finfo(float).max, np.finfo(float).tiny)

    def test_hydrostatic_force_is_removed_before_upward_projection(self) -> None:
        gravity = np.asarray([0.0, -10.0])
        density = np.asarray([1000.0, 1000.0])
        hydrostatic_gradient = np.asarray([[0.0, -10000.0], [0.0, -10000.0]])
        np.testing.assert_allclose(
            bridge.upward_excess_force_density(
                hydrostatic_gradient, density, gravity
            ),
            np.zeros(2),
        )

    def test_raw_finite_difference_gradient_recovers_affine_field(self) -> None:
        reference = np.asarray(
            [
                [0.0, 0.0],
                [1.0, 0.0],
                [2.0, 0.0],
                [0.0, 1.0],
                [1.0, 1.0],
                [2.0, 1.0],
            ]
        )
        current = reference + np.asarray([0.2, -0.1])
        pressures = 7.0 + 3.0 * current[:, 0] - 4.0 * current[:, 1]
        gradient = bridge.raw_finite_difference_gradient(
            reference, current, pressures
        )
        np.testing.assert_allclose(gradient, np.tile([3.0, -4.0], (6, 1)))

    def test_raw_gradient_recovers_affine_field_on_skewed_lattice(self) -> None:
        reference = np.asarray(
            [
                [0.0, 0.0],
                [1.0, 0.0],
                [2.0, 0.0],
                [0.0, 1.0],
                [1.0, 1.0],
                [2.0, 1.0],
                [0.0, 2.0],
                [1.0, 2.0],
                [2.0, 2.0],
            ]
        )
        deformation = np.asarray([[1.1, 0.35], [-0.2, 0.9]])
        current = reference @ deformation + np.asarray([0.2, -0.1])
        pressures = 7.0 + 3.0 * current[:, 0] - 4.0 * current[:, 1]
        gradient = bridge.raw_finite_difference_gradient(
            reference, current, pressures
        )
        np.testing.assert_allclose(
            gradient, np.tile([3.0, -4.0], (len(reference), 1)), atol=1.0e-12
        )

    def test_raw_gradient_rejects_collapsed_direction_system(self) -> None:
        reference = np.asarray(
            [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]]
        )
        current = np.column_stack(
            (reference[:, 0] + reference[:, 1], np.zeros(len(reference)))
        )
        pressures = np.arange(len(reference), dtype=float)
        with self.assertRaisesRegex(ValueError, "direction system is singular"):
            bridge.raw_finite_difference_gradient(reference, current, pressures)

    def test_raw_force_metrics_use_signed_and_positive_integrals(self) -> None:
        force = np.asarray([20.0, -10.0, 40.0])
        gamma = np.asarray([10.0, 10.0, 20.0])
        volume = np.asarray([0.1, 0.2, 0.3])
        support = np.asarray([True, True, False])
        metrics, upward_if = bridge.raw_force_metrics(
            force, gamma, volume, support
        )
        self.assertAlmostEqual(metrics.signed_force_n_per_m, 0.0)
        self.assertAlmostEqual(metrics.positive_force_n_per_m, 2.0)
        self.assertAlmostEqual(metrics.positive_force_index_area_m2, 0.2)
        self.assertAlmostEqual(metrics.maximum_upward_if, 2.0)
        self.assertAlmostEqual(metrics.critical_area_m2, 0.1)
        np.testing.assert_allclose(upward_if, np.asarray([2.0, 0.0, 2.0]))

    def test_display_smoothing_does_not_mutate_raw_values(self) -> None:
        points = np.asarray(
            [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]], dtype=float
        )
        values = np.asarray([0.0, 1.0, 2.0, 1.0])
        original = values.copy()
        config = {
            "mesh": {"cellsize_min": 1.0},
            "analysis": {
                "rigid_pipeline": {
                    "center": [10.0, 10.0],
                    "outer_radius": 0.1,
                }
            },
        }
        _, _, raster = bridge.display_smoothed_grid(
            points,
            values,
            config,
            np.asarray([10.0, 10.0]),
            (0.0, 1.0, 0.0, 1.0),
            0.75,
            nx=20,
            ny=20,
        )
        np.testing.assert_array_equal(values, original)
        self.assertTrue(np.any(np.isfinite(raster)))


if __name__ == "__main__":
    unittest.main()
