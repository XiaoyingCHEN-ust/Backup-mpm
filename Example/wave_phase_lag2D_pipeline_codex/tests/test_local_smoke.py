from __future__ import annotations

import copy
import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch


CASE_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(CASE_DIR))

import prepare_local_smoke as smoke  # noqa: E402
import prepare_study as study  # noqa: E402
import submit_study_sbatch as submitter  # noqa: E402
import validate_case as validate  # noqa: E402


EQUILIBRIUM_VTP = """<?xml version="1.0"?>
<VTKFile type="PolyData" version="0.1" byte_order="LittleEndian">
  <PolyData>
    <Piece NumberOfPoints="1" NumberOfVerts="0" NumberOfLines="0" NumberOfStrips="0" NumberOfPolys="0">
      <PointData>
        <DataArray type="Float64" Name="velocities" NumberOfComponents="3" format="ascii">0 0 0</DataArray>
        <DataArray type="Float64" Name="displacements" NumberOfComponents="3" format="ascii">0 0 0</DataArray>
        <DataArray type="Float64" Name="porosities" format="ascii">0.485</DataArray>
      </PointData>
      <Points><DataArray type="Float64" NumberOfComponents="3" format="ascii">0 0 0</DataArray></Points>
    </Piece>
  </PolyData>
</VTKFile>
"""


class LocalSmokeProfileTest(unittest.TestCase):
    def test_all_roles_have_exact_registered_smoke_physics(self):
        configs = {
            role: smoke.canonical_config("audit", role)
            for role in smoke.ROLE_FILENAMES
        }
        for role, config in configs.items():
            validate.validate_local_smoke_config(config)
            analysis = config["analysis"]
            soil, fluid = config["materials"]
            self.assertEqual(analysis["dt"], 1.0e-4)
            self.assertEqual(analysis["nsteps"], 5_000)
            self.assertEqual(
                config["post_processing"]["path"],
                "results/local_smoke/audit/",
            )
            self.assertNotIn("prescribed_phase_pressures", analysis)
            self.assertEqual(config["local_smoke"]["role"], role)

            if role in {"EQ_LS", "EQ_HS", "MC_EQ"}:
                self.assertFalse(analysis["APIC"])
                self.assertEqual(analysis["PIC"], 1.0)
                self.assertEqual(analysis["PIC_T"], 0.0)
                self.assertEqual(analysis["damping"]["damping_factor"], 5.0)
                self.assertTrue(analysis["rigid_pipeline"]["fixed"])
                self.assertFalse(fluid["wave_pressure"])
                self.assertTrue(analysis["stability_gate"])
                self.assertTrue(config["post_processing"]["write_hdf5"])
                self.assertEqual(
                    config["post_processing"]["output_steps"], 5_000
                )
            else:
                self.assertTrue(analysis["APIC"])
                self.assertEqual(analysis["PIC"], 0.0)
                self.assertEqual(analysis["PIC_T"], 0.0)
                self.assertEqual(analysis["damping"]["damping_factor"], 0.0)
                self.assertFalse(analysis["rigid_pipeline"]["fixed"])
                self.assertTrue(fluid["wave_pressure"])
                self.assertNotIn("stability_gate", analysis)
                self.assertIn(
                    "PIC_pore_pressure_excess",
                    config["post_processing"]["liquid_vtk"],
                )

            if role.startswith("EQ_"):
                self.assertEqual(soil["type"], "LinearElastic2D")
                self.assertFalse(analysis["resume"]["resume"])
                # Even the disabled self-resume metadata is shortened, so a
                # future toggle cannot silently point at the 40k checkpoint.
                self.assertEqual(analysis["resume"]["step"], 5_000)
                self.assertEqual(analysis["resume"]["nsteps"], 5_000)

        prefix = smoke.uuid_prefix("audit")
        self.assertEqual(
            configs["SANI_LS"]["analysis"]["resume"]["uuid"],
            f"{prefix}_LS_EQ",
        )
        self.assertEqual(
            configs["SANI_HS"]["analysis"]["resume"]["uuid"],
            f"{prefix}_HS_EQ",
        )
        self.assertEqual(
            configs["MC_EQ"]["analysis"]["resume"]["uuid"],
            f"{prefix}_HS_EQ",
        )
        self.assertEqual(
            configs["MC_DYNAMIC"]["analysis"]["resume"]["uuid"],
            f"{prefix}_HS_MC_RELAX",
        )
        self.assertEqual(
            configs["MC_DYNAMIC"]["materials"][0]["type"],
            "MohrCoulomb2D",
        )

    def test_recomputed_marker_hash_cannot_authorise_parameter_drift(self):
        config = smoke.canonical_config("tamper", "EQ_LS")
        config["analysis"]["nsteps"] = 4_999
        payload = copy.deepcopy(config)
        payload.pop("local_smoke")
        config["local_smoke"]["baseline_sha256"] = smoke.payload_sha256(payload)
        with self.assertRaisesRegex(ValueError, "exact canonical profile"):
            validate.validate_local_smoke_config(config)

    def test_marker_cannot_convert_formal_or_arbitrary_json_to_smoke(self):
        formal = study.equilibrium_config(
            "SCREEN_LS",
            study.LOW_LAG_SATURATION,
            "results/screen/",
            ".",
            0.02,
            0.01,
            0.12,
        )
        formal["local_smoke"] = {
            "schema": smoke.SMOKE_SCHEMA,
            "label": "audit",
            "role": "EQ_LS",
            "baseline_sha256": "0" * 64,
            "result_eligibility": (
                "local-smoke-only; excluded from manuscript evidence"
            ),
        }
        with self.assertRaisesRegex(ValueError, "baseline hash"):
            validate.validate_local_smoke_config(formal)

    def test_formal_validator_still_rejects_short_smoke_uuid(self):
        config = smoke.canonical_config("formal-separation", "EQ_LS")
        with self.assertRaises(ValueError):
            study.validate_config(config)

    def test_generation_is_non_overwriting_and_manifest_is_auditable(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary) / "configs"
            manifest_path = smoke.generate("run1", root)
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            self.assertEqual(manifest["schema"], "pipeline-local-smoke-manifest-v1")
            self.assertEqual(manifest["steps_per_stage"], 5_000)
            self.assertEqual(manifest["dt_s"], 1.0e-4)
            self.assertEqual(manifest["result_path"], "results/local_smoke/run1/")
            self.assertFalse(manifest["manuscript_evidence_eligible"])
            self.assertEqual(len(manifest["cases"]), 6)
            with self.assertRaises(FileExistsError):
                smoke.generate("run1", root)

    def test_local_completion_audits_vtp_hdf5_but_hpc_skip_rejects_it(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            config = smoke.canonical_config("complete", "EQ_LS")
            config_path = root / "configs" / "local_smoke" / "complete" / "eq.json"
            config_path.parent.mkdir(parents=True)
            config_path.write_text(json.dumps(config), encoding="utf-8")
            result = (
                root
                / "results"
                / "local_smoke"
                / "complete"
                / config["analysis"]["uuid"]
            )
            result.mkdir(parents=True)
            (result / "particle5000.vtp").write_text(
                EQUILIBRIUM_VTP, encoding="utf-8"
            )
            (result / "particles5000.h5").write_bytes(
                validate.HDF5_SIGNATURE + b"local-smoke-checkpoint"
            )
            with patch.object(validate, "CASE_ROOT", root):
                sentinel = validate.write_case_completion(
                    config_path, validation_profile=smoke.SMOKE_SCHEMA
                )
            payload = json.loads(sentinel.read_text(encoding="utf-8"))
            self.assertEqual(
                payload["config"]["validation_profile"], smoke.SMOKE_SCHEMA
            )
            self.assertIn("final_vtp", payload["artifacts"])
            self.assertIn("final_hdf5", payload["artifacts"])
            self.assertIn("max|v|=0.000e+00", payload["stability_qa"])
            with patch.object(submitter, "CASE_DIR", root):
                self.assertFalse(submitter.config_complete(config_path))

    def test_wrappers_make_smoke_explicit_locked_and_fail_fast(self):
        local_wrapper = (CASE_DIR / "run_case_local.sh").read_text(encoding="utf-8")
        sequence = (CASE_DIR / "run_local_smoke.sh").read_text(encoding="utf-8")
        self.assertIn("set -euo pipefail", local_wrapper)
        self.assertIn("--local-smoke", local_wrapper)
        self.assertIn("flock -n", local_wrapper)
        self.assertLess(
            local_wrapper.index("flock -n"),
            local_wrapper.index('"${mpm_bin}" -p'),
        )
        self.assertLess(
            local_wrapper.index("--clear-completion"),
            local_wrapper.index('"${mpm_bin}" -p'),
        )
        self.assertGreater(
            local_wrapper.index("--write-completion"),
            local_wrapper.index('"${mpm_bin}" -p'),
        )
        self.assertIn("set -euo pipefail", sequence)
        self.assertIn("bash run_case_local.sh --smoke", sequence)


if __name__ == "__main__":
    unittest.main()
