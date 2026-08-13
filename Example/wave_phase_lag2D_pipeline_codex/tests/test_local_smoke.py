from __future__ import annotations

import copy
import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch


CASE_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(CASE_DIR))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import prepare_local_smoke as smoke  # noqa: E402
import prepare_study as study  # noqa: E402
import qa_contract as qa  # noqa: E402
import submit_study_sbatch as submitter  # noqa: E402
import validate_case as validate  # noqa: E402
from qa_test_helpers import (  # noqa: E402
    write_fake_hdf5_auditor,
    write_linear_qa_fixture,
)


class LocalSmokeProfileTest(unittest.TestCase):
    def test_required_and_diagnostic_roles_have_exact_smoke_physics(self):
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
                self.assertEqual(validate.expected_particle_steps(config), [5_000])
                expected_mode = (
                    qa.MC_HANDOFF_MODE
                    if role == "MC_EQ"
                    else qa.LINEAR_EQUILIBRIUM_MODE
                )
                self.assertEqual(
                    analysis["stability_qa_contract"],
                    qa.stability_contract(expected_mode),
                )
            else:
                self.assertTrue(analysis["APIC"])
                self.assertEqual(analysis["PIC"], 0.0)
                self.assertEqual(analysis["PIC_T"], 0.0)
                self.assertEqual(analysis["damping"]["damping_factor"], 0.0)
                self.assertEqual(
                    validate.expected_particle_steps(config),
                    list(range(500, 5_001, 500)),
                )
                self.assertFalse(analysis["rigid_pipeline"]["fixed"])
                self.assertTrue(fluid["wave_pressure"])
                self.assertNotIn("stability_gate", analysis)
                self.assertIn(
                    "PIC_pore_pressure_excess",
                    config["post_processing"]["liquid_vtk"],
                )

            if role not in {"EQ_LS", "EQ_HS"}:
                self.assertEqual(
                    analysis["resume_stability_qa_contract"],
                    qa.stability_contract(qa.LINEAR_EQUILIBRIUM_MODE),
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
            tuple(smoke.REQUIRED_ROLES),
            ("EQ_LS", "EQ_HS", "SANI_LS", "SANI_HS"),
        )
        self.assertEqual(tuple(smoke.DIAGNOSTIC_ROLES), ("MC_EQ",))
        self.assertNotIn("MC_DYNAMIC", smoke.ROLE_FILENAMES)
        with self.assertRaisesRegex(ValueError, "Unknown local smoke role"):
            smoke.canonical_config("audit", "MC_DYNAMIC")

    def test_recomputed_marker_hash_cannot_authorise_parameter_drift(self):
        config = smoke.canonical_config("tamper", "EQ_LS")
        config["analysis"]["nsteps"] = 4_999
        payload = copy.deepcopy(config)
        payload.pop("local_smoke")
        config["local_smoke"]["baseline_sha256"] = smoke.payload_sha256(payload)
        with self.assertRaisesRegex(ValueError, "exact canonical profile"):
            validate.validate_local_smoke_config(config)

    def test_recomputed_marker_hash_cannot_authorise_qa_contract_drift(self):
        config = smoke.canonical_config("qa-tamper", "EQ_LS")
        config["analysis"]["stability_qa_contract"]["limits"][
            "maximum_velocity_m_s"
        ] = 2.0e-3
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
            self.assertEqual(manifest["schema"], smoke.MANIFEST_SCHEMA)
            self.assertEqual(manifest["steps_per_stage"], 5_000)
            self.assertEqual(manifest["dt_s"], 1.0e-4)
            self.assertEqual(manifest["result_path"], "results/local_smoke/run1/")
            self.assertFalse(manifest["manuscript_evidence_eligible"])
            self.assertEqual(
                manifest["required_roles"], list(smoke.REQUIRED_ROLES)
            )
            self.assertEqual(manifest["diagnostic_roles"], [])
            self.assertEqual(len(manifest["cases"]), 4)
            self.assertTrue(
                all(
                    case["classification"] == "required"
                    for case in manifest["cases"]
                )
            )
            self.assertTrue(
                all(case["required_for_smoke_pass"] for case in manifest["cases"])
            )
            self.assertEqual(smoke.validate_manifest(manifest_path), manifest)
            self.assertFalse((root / "02_MC_EQ.json").exists())
            self.assertFalse((root / "03_MC_DYNAMIC.json").exists())
            with self.assertRaises(FileExistsError):
                smoke.generate("run1", root)

    def test_opt_in_manifest_adds_only_fail_fast_mc_eq_diagnostic(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary) / "configs"
            manifest_path = smoke.generate(
                "diagnostic", root, include_mc_diagnostic=True
            )
            manifest = smoke.validate_manifest(manifest_path)
            self.assertEqual(manifest["diagnostic_roles"], ["MC_EQ"])
            self.assertEqual(len(manifest["cases"]), 5)
            diagnostic = manifest["cases"][-1]
            self.assertEqual(diagnostic["role"], "MC_EQ")
            self.assertEqual(diagnostic["classification"], "diagnostic")
            self.assertFalse(diagnostic["required_for_smoke_pass"])
            self.assertEqual(diagnostic["dependencies"], ["EQ_HS"])
            config = json.loads(
                (root / diagnostic["config"]).read_text(encoding="utf-8")
            )
            self.assertEqual(
                config["local_smoke"]["result_eligibility"],
                "local-smoke-only; excluded from manuscript evidence",
            )
            self.assertFalse((root / "03_MC_DYNAMIC.json").exists())

    def test_manifest_and_config_hashes_reject_tampering(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary) / "configs"
            manifest_path = smoke.generate("hash-audit", root)
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            manifest["cases"][0]["classification"] = "diagnostic"
            manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "manifest hash"):
                smoke.validate_manifest(manifest_path)

            manifest_path = smoke.generate("hash-audit", root, force=True)
            first_config = root / smoke.ROLE_FILENAMES[smoke.REQUIRED_ROLES[0]]
            first_config.write_text(
                first_config.read_text(encoding="utf-8") + "\n",
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "manifest hash"):
                smoke.validate_manifest(manifest_path)

    def test_force_refuses_stale_removed_mc_dynamic_config(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary) / "configs"
            smoke.generate("stale", root)
            (root / "03_MC_DYNAMIC.json").write_text("{}", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "unexpected files"):
                smoke.generate("stale", root, force=True)

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
            write_linear_qa_fixture(root, config, result / "particles5000.h5")
            auditor = write_fake_hdf5_auditor(root)
            with patch.dict(os.environ, {"MPM_HDF5_AUDITOR": str(auditor)}), patch.object(
                validate, "CASE_ROOT", root
            ):
                sentinel = validate.write_case_completion(
                    config_path, validation_profile=smoke.SMOKE_SCHEMA
                )
            payload = json.loads(sentinel.read_text(encoding="utf-8"))
            self.assertEqual(
                payload["config"]["validation_profile"], smoke.SMOKE_SCHEMA
            )
            self.assertEqual(
                payload["artifacts"]["particle_vtp_grid"]["steps"], [5_000]
            )
            self.assertIn("pipeline_history_csv", payload["artifacts"])
            self.assertIn("final_hdf5", payload["artifacts"])
            self.assertTrue(
                payload["artifacts"]["hdf5_vtp_crosscheck"]["passed"]
            )
            self.assertIn(
                "max|v|=0.000e+00", payload["stability_qa"]["summary"]
            )
            with patch.object(submitter, "CASE_DIR", root):
                self.assertFalse(submitter.config_complete(config_path))

            sentinel.unlink()
            with patch.dict(
                os.environ,
                {"MPM_HDF5_AUDITOR": str(root / "missing-auditor")},
            ), patch.object(validate, "CASE_ROOT", root):
                with self.assertRaises(FileNotFoundError):
                    validate.write_case_completion(
                        config_path, validation_profile=smoke.SMOKE_SCHEMA
                    )
            self.assertFalse(sentinel.exists())

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
        self.assertIn("--mc-diagnostic", sequence)
        self.assertIn("--validate-manifest", sequence)
        self.assertIn("02_MC_EQ.json", sequence)
        self.assertNotIn("MC_DYNAMIC", sequence)
        self.assertNotIn("|| true", sequence)

    def test_default_runs_four_required_roles_and_mc_failure_is_nonzero(self):
        with tempfile.TemporaryDirectory() as temporary:
            fake_root = Path(temporary)
            log = fake_root / "calls.txt"
            fake_python = fake_root / "python"
            fake_python.write_text("#!/bin/sh\nexit 0\n", encoding="utf-8")
            fake_python.chmod(0o755)
            fake_bash = fake_root / "bash"
            fake_bash.write_text(
                "#!/bin/sh\n"
                'printf "%s\\n" "$*" >> "$SMOKE_CALL_LOG"\n'
                'case "$*" in *02_MC_EQ.json*) exit "${MC_STATUS:-0}";; esac\n'
                "exit 0\n",
                encoding="utf-8",
            )
            fake_bash.chmod(0o755)
            environment = os.environ.copy()
            environment.update(
                {
                    "PATH": f"{fake_root}:{environment['PATH']}",
                    "PYTHON_BIN": str(fake_python),
                    "SMOKE_CALL_LOG": str(log),
                }
            )
            wrapper = CASE_DIR / "run_local_smoke.sh"
            default = subprocess.run(
                ["/bin/bash", str(wrapper), "default-test"],
                cwd=CASE_DIR,
                env=environment,
                check=False,
                capture_output=True,
                text=True,
            )
            self.assertEqual(default.returncode, 0, default.stderr)
            default_calls = log.read_text(encoding="utf-8").splitlines()
            self.assertEqual(len(default_calls), 4)
            self.assertFalse(any("MC_EQ" in call for call in default_calls))

            log.write_text("", encoding="utf-8")
            environment["MC_STATUS"] = "17"
            diagnostic = subprocess.run(
                [
                    "/bin/bash",
                    str(wrapper),
                    "--mc-diagnostic",
                    "diagnostic-test",
                ],
                cwd=CASE_DIR,
                env=environment,
                check=False,
                capture_output=True,
                text=True,
            )
            self.assertEqual(diagnostic.returncode, 17)
            diagnostic_calls = log.read_text(encoding="utf-8").splitlines()
            self.assertEqual(len(diagnostic_calls), 5)
            self.assertIn("02_MC_EQ.json", diagnostic_calls[-1])


if __name__ == "__main__":
    unittest.main()
