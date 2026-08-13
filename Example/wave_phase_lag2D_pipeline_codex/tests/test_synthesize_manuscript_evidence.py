from __future__ import annotations

import csv
import json
import sys
import tempfile
import unittest
from pathlib import Path


CASE_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(CASE_DIR))

import synthesize_manuscript_evidence as synthesis  # noqa: E402


def threshold_particles(if_area: float, stress_area: float) -> dict[str, float]:
    return {
        "max_upward_seepage_IF": 1.5,
        "minimum_vertical_stress_remaining_ratio": 0.01,
        "nominal_reference_area_per_particle_m2": 0.01,
        "saved_particle_frame_times_s": [0.0, 0.1, 0.2],
        (
            "support_ROI_upward_seepage_IF_ge_1_"
            "time_integrated_area_m2_s"
        ): if_area,
        (
            "support_ROI_stress_loss_Rsigma_le_0p05_"
            "time_integrated_area_m2_s"
        ): stress_area,
    }


def phase(lag: float) -> dict[str, object]:
    probes = {
        name: {
            "mean_pressure_pa": 1000.0 + index,
            "amplitude_pa": 500.0 - 10.0 * index,
            "amplitude_ratio_to_surface": 0.8 - 0.05 * index,
            "phase_lag_deg": lag,
        }
        for index, name in enumerate(("crown", "shoulder", "invert"))
    }
    return {
        "available": True,
        "surface_amplitude_pa": 600.0,
        "probes": probes,
    }


def summary(code: str, lag: float, if_area: float, stress_area: float) -> dict:
    return {
        "case": f"04_{code}" if code in ("RL", "RM", "RE") else f"02_{code}",
        "phase": phase(lag),
        "particles": threshold_particles(if_area, stress_area),
        "pipeline": {
            "maximum_abs_vertical_displacement_over_D": 0.6,
            "post_release_no_contact_fraction": 0.2,
            "maximum_abs_rotation_rad": 0.1,
        },
    }


class ManuscriptEvidenceTest(unittest.TestCase):
    def write_study(self, root: Path, *, reverse_phase_effect: bool = False) -> Path:
        analysis = root / "analysis"
        analysis.mkdir()
        cases = {
            "LS": summary("LS", 5.0, 0.2, 0.1),
            "HS": summary("HS", 25.0, 0.5, 0.4),
            "HM": summary("HM", 25.0, 0.4, 0.3),
            "RL": summary(
                "RL", 25.0, 0.2 if reverse_phase_effect else 0.8, 0.7
            ),
            "RM": summary("RM", 25.0, 0.4, 0.3),
            "RE": summary(
                "RE", 1.0, 0.9 if reverse_phase_effect else 0.2, 0.2
            ),
        }
        cases["RM"]["pipeline"][
            "maximum_abs_vertical_displacement_over_D"
        ] = 0.2
        for code, data in cases.items():
            prefix = "04" if code in ("RL", "RM", "RE") else "02"
            (analysis / f"{prefix}_{code}_summary.json").write_text(
                json.dumps(data), encoding="utf-8"
            )

        fields = [
            "time_s",
            "crown_p_effective_pa",
            "crown_q_pa",
            "crown_eps_p_q",
            "crown_pdstrain",
            "shoulder_p_effective_pa",
            "shoulder_q_pa",
            "shoulder_eps_p_q",
            "shoulder_pdstrain",
            "invert_p_effective_pa",
            "invert_q_pa",
            "invert_eps_p_q",
            "invert_pdstrain",
        ]
        for code, state in (("RL", "eps_p_q"), ("RM", "pdstrain")):
            with (analysis / f"04_{code}_probe_history.csv").open(
                "w", newline="", encoding="utf-8"
            ) as stream:
                writer = csv.DictWriter(stream, fieldnames=fields)
                writer.writeheader()
                for index in range(4):
                    row = {"time_s": index * 0.1}
                    for probe_index, probe in enumerate(
                        ("crown", "shoulder", "invert")
                    ):
                        row[f"{probe}_p_effective_pa"] = 3000 - 200 * index
                        row[f"{probe}_q_pa"] = (
                            500 + 100 * ((-1) ** index) + probe_index
                        )
                        row[f"{probe}_{state}"] = index * 0.01
                    writer.writerow(row)

        metadata = {
            "schema": "phase-erased-pressure-control-v1",
            "classification": "one-way numerical counterfactual; not fully coupled",
            "source": {"points_sha256": "same"},
            "output": {"points_sha256": "same"},
            "qa": {
                "max_abs_mean_change_pa": 0.0,
                "max_abs_fundamental_amplitude_change_pa": 1.0e-12,
            },
        }
        metadata_path = root / "phase_erased_metadata.json"
        metadata_path.write_text(json.dumps(metadata), encoding="utf-8")
        return metadata_path

    def test_full_evidence_supports_phase_and_mechanistic_claims(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            metadata = self.write_study(root)
            evidence = synthesis.synthesize(root / "analysis", metadata)
        gate = evidence["claim_gate"]
        self.assertTrue(gate["phase_lag_hydraulic_trigger_supported"])
        self.assertTrue(gate["phase_lag_realised_liquefaction_supported"])
        self.assertTrue(gate["sanisand_mechanistic_advantage_supported"])
        self.assertEqual(
            evidence["constitutive_RL_RM"]["status"],
            "mechanistic_advantage_observed",
        )
        resolution = evidence["phase_only_RL_RE"]["hydraulic_trigger"][
            "absolute_resolution"
        ]
        self.assertTrue(resolution["available"])
        self.assertAlmostEqual(
            resolution["required_absolute_difference_m2_s"], 0.001
        )
        self.assertFalse(
            evidence["constitutive_RL_RM"]["fabric_state_outputs_available"]
        )
        report = synthesis.markdown_report(evidence)
        self.assertIn("one-way numerical counterfactual", report)
        self.assertIn("not universal superiority", report)

    def test_opposite_phase_effect_blocks_favourable_claim(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            metadata = self.write_study(root, reverse_phase_effect=True)
            evidence = synthesis.synthesize(root / "analysis", metadata)
        self.assertFalse(
            evidence["claim_gate"]["phase_lag_hydraulic_trigger_supported"]
        )
        self.assertFalse(
            evidence["claim_gate"]["phase_lag_realised_liquefaction_supported"]
        )
        self.assertEqual(
            evidence["phase_only_RL_RE"]["hydraulic_trigger"]["status"],
            "opposite",
        )

    def test_absolute_area_time_resolution_blocks_unresolved_five_percent_change(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            metadata = self.write_study(root)
            for code, if_area, stress_area in (
                ("RL", 0.0105, 0.0105),
                ("RE", 0.0100, 0.0100),
            ):
                path = root / "analysis" / f"04_{code}_summary.json"
                data = json.loads(path.read_text(encoding="utf-8"))
                data["particles"].update(threshold_particles(if_area, stress_area))
                path.write_text(json.dumps(data), encoding="utf-8")
            evidence = synthesis.synthesize(root / "analysis", metadata)
        phase = evidence["phase_only_RL_RE"]
        self.assertEqual(
            phase["hydraulic_trigger"]["status"], "not_materially_resolved"
        )
        self.assertEqual(
            phase["realised_skeleton_stress_loss"]["status"],
            "not_materially_resolved",
        )
        self.assertFalse(
            evidence["claim_gate"]["phase_lag_hydraulic_trigger_supported"]
        )
        self.assertFalse(
            evidence["claim_gate"]["phase_lag_realised_liquefaction_supported"]
        )

    def test_monotonic_plastic_accumulation_is_not_cyclic_state_evidence(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            metadata = self.write_study(root)
            path = root / "analysis" / "04_RL_probe_history.csv"
            with path.open(newline="", encoding="utf-8") as stream:
                rows = list(csv.DictReader(stream))
            for index, row in enumerate(rows):
                for probe_index, probe in enumerate(
                    ("crown", "shoulder", "invert")
                ):
                    row[f"{probe}_q_pa"] = str(
                        600.0 + probe_index + 100.0 * index
                    )
            with path.open("w", newline="", encoding="utf-8") as stream:
                writer = csv.DictWriter(stream, fieldnames=rows[0].keys())
                writer.writeheader()
                writer.writerows(rows)
            evidence = synthesis.synthesize(root / "analysis", metadata)
        constitutive = evidence["constitutive_RL_RM"]
        self.assertFalse(constitutive["sanisand_cyclic_state_evolution_active"])
        self.assertEqual(constitutive["cyclic_probe_count"], 0)
        self.assertEqual(constitutive["status"], "field_evidence_insufficient")
        self.assertFalse(
            evidence["claim_gate"]["sanisand_mechanistic_advantage_supported"]
        )
        report = synthesis.markdown_report(evidence)
        self.assertNotIn("resolves cyclic state evolution", report)

    def test_state_activity_without_response_separation_is_not_overclaimed(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            metadata = self.write_study(root)
            rm_path = root / "analysis" / "04_RM_summary.json"
            rl_path = root / "analysis" / "04_RL_summary.json"
            rm = json.loads(rm_path.read_text(encoding="utf-8"))
            rl = json.loads(rl_path.read_text(encoding="utf-8"))
            rm["particles"] = rl["particles"]
            rm["pipeline"] = rl["pipeline"]
            rm_path.write_text(json.dumps(rm), encoding="utf-8")
            evidence = synthesis.synthesize(root / "analysis", metadata)
        self.assertEqual(
            evidence["constitutive_RL_RM"]["status"],
            "cyclic_state_active_but_response_not_materially_separated",
        )
        self.assertFalse(
            evidence["claim_gate"]["sanisand_mechanistic_advantage_supported"]
        )

    def test_initial_state_mismatch_blocks_clean_constitutive_claim(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            metadata = self.write_study(root)
            path = root / "analysis" / "04_RM_probe_history.csv"
            with path.open(newline="", encoding="utf-8") as stream:
                rows = list(csv.DictReader(stream))
            rows[0]["crown_p_effective_pa"] = "1000.0"
            with path.open("w", newline="", encoding="utf-8") as stream:
                writer = csv.DictWriter(stream, fieldnames=rows[0].keys())
                writer.writeheader()
                writer.writerows(rows)
            evidence = synthesis.synthesize(root / "analysis", metadata)
        constitutive = evidence["constitutive_RL_RM"]
        self.assertFalse(constitutive["initial_state_QA"]["passed"])
        self.assertEqual(
            constitutive["status"],
            "initial_state_mismatch_confounds_constitutive_ablation",
        )
        self.assertFalse(
            evidence["claim_gate"]["sanisand_mechanistic_advantage_supported"]
        )

    def test_missing_full_replay_summary_is_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            metadata = self.write_study(root)
            (root / "analysis" / "04_RE_summary.json").unlink()
            with self.assertRaisesRegex(FileNotFoundError, "RE"):
                synthesis.synthesize(root / "analysis", metadata)

    def test_phase_only_evidence_does_not_require_mc_branch(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            metadata = self.write_study(root)
            for code, prefix in (("HM", "02"), ("RM", "04")):
                (root / "analysis" / f"{prefix}_{code}_summary.json").unlink()
            (root / "analysis" / "04_RM_probe_history.csv").unlink()
            evidence = synthesis.synthesize(
                root / "analysis", metadata, include_constitutive=False
            )
        self.assertTrue(
            evidence["claim_gate"]["phase_lag_hydraulic_trigger_supported"]
        )
        self.assertTrue(
            evidence["claim_gate"]["phase_lag_realised_liquefaction_supported"]
        )
        self.assertFalse(
            evidence["claim_gate"]["sanisand_mechanistic_advantage_supported"]
        )
        self.assertEqual(evidence["constitutive_RL_RM"]["status"], "unavailable")
        self.assertIn("No SANISAND-versus-Mohr-Coulomb", synthesis.markdown_report(evidence))


if __name__ == "__main__":
    unittest.main()
