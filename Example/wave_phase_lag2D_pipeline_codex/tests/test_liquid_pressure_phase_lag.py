from __future__ import annotations

import unittest

import numpy as np

import plot_liquid_pressure_phase_lag as phase_plot


class LiquidPressurePhaseLagTests(unittest.TestCase):
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


if __name__ == "__main__":
    unittest.main()
