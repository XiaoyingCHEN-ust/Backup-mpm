from __future__ import annotations

import csv
import json
import sys
import tempfile
import unittest
from pathlib import Path


CASE_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(CASE_DIR))

import plot_manuscript_figures as plotting  # noqa: E402


class ManuscriptFigureTest(unittest.TestCase):
    def write_inputs(self, root: Path) -> Path:
        analysis = root / "analysis"
        analysis.mkdir()
        result_root = root / "results"
        result_root.mkdir()
        for code, prefix in plotting.CASE_PREFIX.items():
            result = result_root / code
            result.mkdir()
            for step in (0, 10, 20):
                (result / f"particle{step:07d}.vtp").write_text(
                    "field", encoding="utf-8"
                )
            config = root / f"{code}.json"
            config.write_text(
                json.dumps({"analysis": {"dt": 0.1}}), encoding="utf-8"
            )
            summary = {
                "case": prefix,
                "config": str(config),
                "result_directory": str(result),
                "phase": {
                    "surface_amplitude_pa": 600.0,
                    "probes": {
                        probe: {
                            "phase_lag_deg": 5.0 if code == "LS" else 25.0,
                            "amplitude_ratio_to_surface": 0.8,
                        }
                        for probe in ("crown", "shoulder", "invert")
                    },
                },
            }
            (analysis / f"{prefix}_summary.json").write_text(
                json.dumps(summary), encoding="utf-8"
            )

            history_fields = [
                "time_s",
                "surface_pressure_pa",
                "crown_pressure_pa",
                "invert_pressure_pa",
                plotting.IF_AREA,
                plotting.STRESS_AREA,
                "invert_p_effective_pa",
                "invert_q_pa",
                "invert_eps_p_q",
                "invert_pdstrain",
            ]
            with (analysis / f"{prefix}_probe_history.csv").open(
                "w", newline="", encoding="utf-8"
            ) as stream:
                writer = csv.DictWriter(stream, fieldnames=history_fields)
                writer.writeheader()
                for index in range(3):
                    writer.writerow(
                        {
                            "time_s": index,
                            "surface_pressure_pa": 500.0 * index,
                            "crown_pressure_pa": 400.0 * index,
                            "invert_pressure_pa": 300.0 * index,
                            plotting.IF_AREA: 0.2 if index == 1 else 0.1,
                            plotting.STRESS_AREA: 0.1 * index,
                            "invert_p_effective_pa": 3000.0 - 100.0 * index,
                            "invert_q_pa": 500.0 + 50.0 * index,
                            "invert_eps_p_q": 0.01 * index,
                            "invert_pdstrain": 0.02 * index,
                        }
                    )
            with (analysis / f"{prefix}_pipeline_cycles.csv").open(
                "w", newline="", encoding="utf-8"
            ) as stream:
                writer = csv.DictWriter(
                    stream,
                    fieldnames=(
                        "cycle_after_release",
                        "time_s",
                        "uplift_over_D",
                    ),
                )
                writer.writeheader()
                for index in range(3):
                    writer.writerow(
                        {
                            "cycle_after_release": index,
                            "time_s": index,
                            "uplift_over_D": 0.1 * index,
                        }
                    )
        return analysis

    def test_plots_and_common_time_manifest_are_written(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            analysis = self.write_inputs(root)
            output = root / "figures"
            plotting.plot_hydraulic_bridge(
                analysis, output / "figure7_hydraulic_bridge"
            )
            plotting.plot_phase_control(
                analysis, output / "figure8_phase_control"
            )
            plotting.plot_constitutive(
                analysis, output / "figure9_constitutive"
            )
            frames = plotting.common_field_frames(analysis)
            self.assertEqual(frames["target_time_s"], 1.0)
            self.assertEqual(frames["cases"]["RL"]["step"], 10)
            self.assertEqual(frames["cases"]["RE"]["step"], 10)
            for stem in (
                "figure7_hydraulic_bridge",
                "figure8_phase_control",
                "figure9_constitutive",
            ):
                self.assertGreater((output / f"{stem}.png").stat().st_size, 1000)
                self.assertGreater((output / f"{stem}.pdf").stat().st_size, 1000)


if __name__ == "__main__":
    unittest.main()
