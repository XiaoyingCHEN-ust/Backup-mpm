from __future__ import annotations

import unittest
from pathlib import Path

import run_sanisand_method_optimization as study
import run_sanisand_validation as base


def error(rmse: float, maximum: float) -> dict:
    return {
        "paired_samples": 10,
        "rmse_normalized": rmse,
        "maximum_absolute_normalized": maximum,
    }


class SanisandMethodOptimizationTests(unittest.TestCase):
    def test_variant_matrix_is_exact(self) -> None:
        self.assertEqual(
            [item[0] for item in study.VARIANTS],
            [
                "mapped_legacy",
                "direct_csl",
                "direct_elasticity",
                "direct_csl_and_elasticity",
            ],
        )
        self.assertEqual(study.HOLDOUT_REGRESSION_LIMIT, 1.10)

    def test_safe_ratio_rejects_nonpositive_denominator(self) -> None:
        with self.assertRaisesRegex(study.MethodStudyError, "denominator"):
            study.safe_ratio(1.0, 0.0)

    def test_maximum_error_uses_both_fields(self) -> None:
        value = study.maximum_error(
            {
                "pressure": error(0.01, 0.02),
                "deviator_stress": error(0.03, 0.09),
            }
        )
        self.assertEqual(value, 0.09)

    @staticmethod
    def selection_records() -> dict:
        records = {}
        for variant, _, _ in study.VARIANTS:
            records[variant] = {
                "monotonic_maximum_normalized": 0.10,
                "cyclic_holdout": {
                    "pressure": {"rmse_normalized": 0.20},
                    "deviator_stress": {"rmse_normalized": 0.30},
                },
                "q_control_terminal_cycle_gap": 1.0,
            }
        return records

    def test_calibration_gain_cannot_rescue_holdout_regression(self) -> None:
        records = self.selection_records()
        records["direct_elasticity"]["monotonic_maximum_normalized"] = 0.08
        records["direct_elasticity"]["cyclic_holdout"]["pressure"][
            "rmse_normalized"
        ] = 0.23
        selected = study.select_parameterization(records)
        self.assertEqual(selected, "mapped_legacy")
        self.assertFalse(
            records["direct_elasticity"]["selection"][
                "eligible_to_replace_baseline"
            ]
        )

    def test_candidate_must_improve_and_pass_both_holdouts(self) -> None:
        records = self.selection_records()
        records["direct_elasticity"]["monotonic_maximum_normalized"] = 0.08
        records["direct_elasticity"]["cyclic_holdout"]["pressure"][
            "rmse_normalized"
        ] = 0.21
        records["direct_elasticity"]["cyclic_holdout"]["deviator_stress"][
            "rmse_normalized"
        ] = 0.31
        records["direct_elasticity"]["q_control_terminal_cycle_gap"] = 1.05
        selected = study.select_parameterization(records)
        self.assertEqual(selected, "direct_elasticity")

    def test_current_manuscript_declares_method_selection_limitation(self) -> None:
        manuscript = (
            Path(base.__file__).resolve().parent
            / "SANISAND_VALIDATION_SUBSECTION.md"
        ).read_text(encoding="utf-8")
        self.assertIn("memory-surface", manuscript)
        self.assertIn("not numerically equivalent", manuscript)


if __name__ == "__main__":
    unittest.main()
