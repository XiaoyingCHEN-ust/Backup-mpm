from __future__ import annotations

import unittest

import analyze_phase_conditioned_uplift as conditioned


def _row(time_s: float, pressure: float, lagged: float, erased: float) -> dict[str, float]:
    row = {
        "time_s": time_s,
        "surface_excess_pressure_pa": pressure,
    }
    for metric in (
        "positive_force_n_per_m",
        "signed_force_n_per_m",
        "net_uplift_force_n_per_m",
        "critical_area_m2",
        "maximum_upward_if",
    ):
        row[f"lagged_{metric}"] = lagged
        row[f"phase_erased_{metric}"] = erased
    return row


class PhaseConditionedUpliftTest(unittest.TestCase):
    def test_recovery_branch_begins_at_trough_and_reaches_cycle_end(self) -> None:
        rows = [
            _row(2.600 + 0.065 * index, pressure, 106.0, 100.0)
            for index, pressure in enumerate(
                (-90.0, -1.0, 90.0, 200.0, 280.0, 200.0, 90.0, 1.0, -90.0,
                 -170.0, -230.0, -284.0, -270.0, -230.0, -168.0, -89.0,
                 -1.0, 90.0, 200.0, 280.0, 200.0)
            )
        ]
        # Use the first 21-point synthetic period and set the requested end to
        # the final core-negative frame, matching the registered saved grid.
        branch = conditioned.recovery_branch(
            rows,
            cycle_start_s=2.6,
            cycle_end_s=3.575,
            surface_amplitude_pa=284.0,
        )
        self.assertTrue(
            all(
                abs(observed - expected) <= 1.0e-12
                for observed, expected in zip(
                    [row["time_s"] for row in branch],
                    [3.315, 3.38, 3.445, 3.51, 3.575],
                )
            )
        )
        self.assertTrue(all(row["surface_excess_pressure_pa"] <= -28.4 for row in branch))

    def test_integrate_branch_reports_local_and_signed_effects(self) -> None:
        rows = [
            _row(3.64 + 0.065 * index, -284.0 + 50.0 * index, 106.0, 100.0)
            for index in range(5)
        ]
        result = conditioned.integrate_branch(rows)
        local = result["local_upward_activity_impulse"]
        signed = result["signed_support_force_impulse"]
        self.assertAlmostEqual(local["lagged_minus_phase_erased_fraction"], 0.06)
        self.assertAlmostEqual(signed["lagged_minus_phase_erased_fraction"], 0.06)
        self.assertAlmostEqual(local["maximum_paired_delta_fraction"], 0.06)

    def test_gate_requires_repeatability_signed_direction_and_phase(self) -> None:
        branch = conditioned.integrate_branch(
            [
                _row(3.64 + 0.065 * index, -284.0 + 50.0 * index, 106.0, 100.0)
                for index in range(5)
            ]
        )
        passed = conditioned.short_time_gate([branch, branch], {"passed": True})
        self.assertTrue(passed["passed"])
        failed = conditioned.short_time_gate([branch, branch], {"passed": False})
        self.assertFalse(failed["passed"])
        weaker = conditioned.integrate_branch(
            [
                _row(3.64 + 0.065 * index, -284.0 + 50.0 * index, 101.0, 100.0)
                for index in range(5)
            ]
        )
        self.assertFalse(
            conditioned.short_time_gate([weaker, weaker], {"passed": True})["passed"]
        )

    def test_relevant_phase_can_be_resolved_below_the_crown(self) -> None:
        audit = {
            "probe_results": {
                "surface": {"fundamental_amplitude_pa": 280.0},
                "crown": {
                    "fundamental_amplitude_pa": 45.0,
                    "harmonic_r_squared": 0.99,
                    "phase_difference_from_same_column_surface_deg": -5.0,
                },
                "shoulder": {
                    "fundamental_amplitude_pa": 10.0,
                    "harmonic_r_squared": 0.99,
                    "phase_difference_from_same_column_surface_deg": -30.0,
                },
                "invert": {
                    "fundamental_amplitude_pa": 28.0,
                    "harmonic_r_squared": 0.98,
                    "phase_difference_from_same_column_surface_deg": -54.0,
                },
            }
        }
        result = conditioned.relevant_pipe_phase(audit)
        self.assertTrue(result["passed"])
        self.assertEqual(result["most_delayed_eligible_probe"]["probe"], "invert")

    def test_relative_fraction_rejects_non_finite_and_handles_zero(self) -> None:
        self.assertIsNone(conditioned.relative_fraction(0.0, 0.0))
        with self.assertRaises(conditioned.ConditionedUpliftError):
            conditioned.relative_fraction(float("nan"), 1.0)


if __name__ == "__main__":
    unittest.main()
