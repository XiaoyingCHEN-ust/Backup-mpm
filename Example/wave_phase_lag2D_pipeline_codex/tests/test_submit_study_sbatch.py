from __future__ import annotations

import json
import struct
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch


CASE_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(CASE_DIR))

import submit_study_sbatch as runner  # noqa: E402


class SbatchWorkflowTest(unittest.TestCase):
    def test_hpc4_resource_flags_are_fixed_and_defaults_target_gpu30(self):
        args = runner.build_parser().parse_args(["screen"])
        command = runner.sbatch_prefix(args, "test", "01:00:00")
        self.assertEqual(command[:2], ["sbatch", "--parsable"])
        self.assertIn("--partition=granularmech", command)
        self.assertIn("--account=comgranmech", command)
        self.assertIn("--cpus-per-task=16", command)
        self.assertIn("--gres=gpu:1", command)
        self.assertIn("--nodelist=gpu30", command)
        self.assertIn("--output=logs/test-%j.out", command)
        self.assertIn("--error=logs/test-%j.err", command)

    def test_afterok_contains_all_dependency_job_ids(self):
        args = runner.build_parser().parse_args(["production"])
        command = runner.sbatch_command(
            args,
            "dependent",
            "36:00:00",
            ["bash", "run_case.sh", "configs/production/02_HS.json"],
            ["12345", "67890", "12345"],
        )
        self.assertIn("--dependency=afterok:12345:67890", command)
        self.assertEqual(command[-4:], [
            runner.RUN_TASK,
            "bash",
            "run_case.sh",
            "configs/production/02_HS.json",
        ])

    def test_full_dag_has_driver_phase_and_constitutive_dependencies(self):
        args = runner.build_parser().parse_args(
            ["screen", "--prepare", "--build", "--analyze", "--dry-run"]
        )
        queue = runner.submit_workflow(args)

        self.assertEqual(queue.jobs["build"].dependency_keys, ["prepare"])
        self.assertEqual(queue.jobs["LS"].dependency_keys, ["EQ_LS"])
        for key in ("HS", "HD"):
            self.assertEqual(queue.jobs[key].dependency_keys, ["EQ_HS"])
        self.assertEqual(queue.jobs["MC_EQ"].dependency_keys, ["EQ_HS"])
        self.assertEqual(queue.jobs["HM"].dependency_keys, ["MC_EQ"])
        for key in ("phase", "RL"):
            self.assertEqual(queue.jobs[key].dependency_keys, ["HD"])
        self.assertEqual(queue.jobs["RM"].dependency_keys, ["HD", "MC_EQ"])
        self.assertEqual(queue.jobs["RE"].dependency_keys, ["phase"])
        self.assertEqual(
            queue.jobs["analysis"].dependency_keys,
            [key for key in queue.jobs if key != "analysis"],
        )

    def test_phase_dag_excludes_failed_mohr_coulomb_branch(self):
        args = runner.build_parser().parse_args(
            [
                "screen",
                "--stage",
                "phase",
                "--prepare",
                "--build",
                "--analyze",
                "--dry-run",
            ]
        )
        queue = runner.submit_workflow(args)
        for code in ("LS", "HS", "HD", "phase", "RL", "RE", "analysis"):
            self.assertIn(code, queue.jobs)
        for code in ("MC_EQ", "HM", "RM"):
            self.assertNotIn(code, queue.jobs)
        self.assertEqual(queue.jobs["RL"].dependency_keys, ["HD"])
        self.assertEqual(queue.jobs["RE"].dependency_keys, ["phase"])

    def test_new_label_does_not_require_a_manifest_at_submit_time(self):
        args = runner.build_parser().parse_args(
            ["screen", "--label", "not_generated_yet", "--prepare", "--dry-run"]
        )
        queue = runner.submit_workflow(args)
        self.assertIn("RM", queue.jobs)
        self.assertEqual(
            queue.jobs["MC_EQ"].command[-1],
            "configs/screen_not_generated_yet/02_MC_EQ.json",
        )
        self.assertEqual(
            queue.jobs["RM"].command[-1],
            "configs/screen_not_generated_yet/04_RM.json",
        )
        self.assertEqual(
            queue.jobs["HM"].command[-1],
            "configs/screen_not_generated_yet/03_HM.json",
        )
        for code in (
            "EQ_LS",
            "EQ_HS",
            "MC_EQ",
            "LS",
            "HS",
            "HM",
            "HD",
            "RL",
            "RM",
            "RE",
        ):
            config_argument = Path(queue.jobs[code].command[-1])
            self.assertFalse(config_argument.is_absolute())
            self.assertTrue(config_argument.as_posix().startswith("configs/"))

    def test_dirty_low_equilibrium_forces_low_dynamic_rerun(self):
        args = runner.build_parser().parse_args(["screen", "--dry-run"])

        def complete(path: Path) -> bool:
            return path.name != "01_EQ_LS.json"

        with patch.object(runner, "config_complete", side_effect=complete), patch.object(
            runner, "phase_control_complete", return_value=True
        ):
            queue = runner.submit_workflow(args)
        self.assertIn("EQ_LS", queue.jobs)
        self.assertIn("LS", queue.jobs)
        self.assertEqual(queue.jobs["LS"].dependency_keys, ["EQ_LS"])
        for code in (
            "EQ_HS",
            "MC_EQ",
            "HS",
            "HM",
            "HD",
            "phase",
            "RL",
            "RM",
            "RE",
        ):
            self.assertNotIn(code, queue.jobs)

    def test_dirty_high_equilibrium_propagates_through_driver_and_replays(self):
        args = runner.build_parser().parse_args(["screen", "--dry-run"])

        def complete(path: Path) -> bool:
            return path.name != "01_EQ_HS.json"

        with patch.object(runner, "config_complete", side_effect=complete), patch.object(
            runner, "phase_control_complete", return_value=True
        ):
            queue = runner.submit_workflow(args)
        for code in (
            "EQ_HS",
            "MC_EQ",
            "HS",
            "HM",
            "HD",
            "phase",
            "RL",
            "RM",
            "RE",
        ):
            self.assertIn(code, queue.jobs)
        self.assertEqual(queue.jobs["HD"].dependency_keys, ["EQ_HS"])
        self.assertEqual(queue.jobs["MC_EQ"].dependency_keys, ["EQ_HS"])
        self.assertEqual(queue.jobs["HM"].dependency_keys, ["MC_EQ"])
        self.assertEqual(queue.jobs["RL"].dependency_keys, ["HD"])
        self.assertEqual(queue.jobs["RM"].dependency_keys, ["HD", "MC_EQ"])
        self.assertEqual(queue.jobs["RE"].dependency_keys, ["phase"])

    def test_dirty_mc_handoff_forces_both_mc_descendants(self):
        args = runner.build_parser().parse_args(["screen", "--dry-run"])

        def complete(path: Path) -> bool:
            return path.name != "02_MC_EQ.json"

        with patch.object(runner, "config_complete", side_effect=complete), patch.object(
            runner, "phase_control_complete", return_value=True
        ):
            queue = runner.submit_workflow(args)
        for code in ("MC_EQ", "HM", "RM"):
            self.assertIn(code, queue.jobs)
        self.assertEqual(queue.jobs["HM"].dependency_keys, ["MC_EQ"])
        self.assertEqual(queue.jobs["RM"].dependency_keys, ["validate", "MC_EQ"])
        for code in ("EQ_HS", "HS", "HD", "phase", "RL", "RE"):
            self.assertNotIn(code, queue.jobs)

    def test_dirty_phase_transform_forces_replay_counterfactual_rerun(self):
        args = runner.build_parser().parse_args(["screen", "--dry-run"])
        with patch.object(runner, "config_complete", return_value=True), patch.object(
            runner, "phase_control_complete", return_value=False
        ):
            queue = runner.submit_workflow(args)
        self.assertIn("phase", queue.jobs)
        self.assertIn("RE", queue.jobs)
        self.assertEqual(queue.jobs["RE"].dependency_keys, ["phase"])
        for code in ("HD", "RL", "RM"):
            self.assertNotIn(code, queue.jobs)

    def test_batch_wrapper_pins_partition_and_account(self):
        wrapper = (CASE_DIR / runner.RUN_TASK).read_text(encoding="utf-8")
        self.assertIn("#SBATCH --partition=granularmech", wrapper)
        self.assertIn("#SBATCH --account=comgranmech", wrapper)
        self.assertIn('case_dir=$(cd "${SLURM_SUBMIT_DIR}" && pwd -P)', wrapper)
        self.assertNotIn('dirname "${BASH_SOURCE[0]}"', wrapper)
        self.assertIn('exec "$@"', wrapper)

    def test_case_wrapper_publishes_completion_only_after_solver(self):
        wrapper = (CASE_DIR / "run_case.sh").read_text(encoding="utf-8")
        clear_position = wrapper.index("--clear-completion")
        solver_position = wrapper.index('"${mpm_bin}" -p')
        completion_position = wrapper.index("--write-completion")
        self.assertLess(clear_position, solver_position)
        self.assertGreater(completion_position, solver_position)

    def test_config_complete_requires_matching_sentinel_and_artifact_hash(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            config_path = root / "configs" / "screen" / "case.json"
            config_path.parent.mkdir(parents=True)
            config = {
                "analysis": {"uuid": "TEST", "nsteps": 10},
                "post_processing": {"path": "results/screen/", "write_hdf5": False},
            }
            config_path.write_text(json.dumps(config), encoding="utf-8")
            result = root / "results" / "screen" / "TEST"
            result.mkdir(parents=True)
            final_vtp = result / "particle10.vtp"
            final_vtp.write_text("complete-vtp", encoding="utf-8")
            sentinel = {
                "schema": runner.COMPLETION_SCHEMA,
                "config": {
                    "path": "configs/screen/case.json",
                    "sha256": runner.file_sha256(config_path),
                    "uuid": "TEST",
                    "nsteps": 10,
                    "validation_profile": "registered-study-v1",
                },
                "artifacts": {
                    "final_vtp": {
                        "path": "results/screen/TEST/particle10.vtp",
                        "size_bytes": final_vtp.stat().st_size,
                        "sha256": runner.file_sha256(final_vtp),
                    }
                },
                "runtime_dependencies": {},
            }
            (result / runner.COMPLETION_FILENAME).write_text(
                json.dumps(sentinel), encoding="utf-8"
            )
            with patch.object(runner, "CASE_DIR", root), patch.object(
                runner.study, "validate_config"
            ):
                self.assertTrue(runner.config_complete(config_path))
                sentinel["config"]["validation_profile"] = "pipeline-local-smoke-v1"
                (result / runner.COMPLETION_FILENAME).write_text(
                    json.dumps(sentinel), encoding="utf-8"
                )
                self.assertFalse(runner.config_complete(config_path))
                sentinel["config"]["validation_profile"] = "registered-study-v1"
                (result / runner.COMPLETION_FILENAME).write_text(
                    json.dumps(sentinel), encoding="utf-8"
                )
                final_vtp.write_text("corrupt-vtp", encoding="utf-8")
                self.assertFalse(runner.config_complete(config_path))

    def test_config_complete_rejects_checkpoint_and_read_pressure_hash_drift(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            config_path = root / "configs" / "screen" / "replay.json"
            config_path.parent.mkdir(parents=True)
            config = {
                "analysis": {
                    "uuid": "DYNAMIC",
                    "nsteps": 10,
                    "resume": {
                        "resume": True,
                        "uuid": "EQ",
                        "step": 10,
                        "nsteps": 10,
                    },
                    "prescribed_phase_pressures": {
                        "enable": True,
                        "write": False,
                        "path": "pressure_databases/screen/lagged",
                        "file_prefix": "pressure",
                        "source_dt": 0.1,
                        "step_interval": 10,
                        "max_step": 10,
                    },
                },
                "post_processing": {
                    "path": "results/screen/",
                    "write_hdf5": False,
                },
            }
            config_path.write_text(json.dumps(config), encoding="utf-8")
            dynamic = root / "results" / "screen" / "DYNAMIC"
            equilibrium = root / "results" / "screen" / "EQ"
            dynamic.mkdir(parents=True)
            equilibrium.mkdir(parents=True)
            final_vtp = dynamic / "particle10.vtp"
            checkpoint = equilibrium / "particles10.h5"
            qa_vtp = equilibrium / "particle10.vtp"
            final_vtp.write_text("dynamic-result", encoding="utf-8")
            checkpoint.write_bytes(b"checkpoint-state")
            qa_vtp.write_text("equilibrium-qa", encoding="utf-8")

            database = root / "pressure_databases" / "screen" / "lagged"
            database.mkdir(parents=True)
            points = database / "pressure_points.txt"
            values = database / "pressure_values.bin"
            points.write_text("0 0.0 0.0\n", encoding="utf-8")
            pressure_payload = bytearray(
                runner.PRESSURE_HEADER.pack(
                    b"MPM_PRESSURE_V2\0", 2, 2, 1, 10, 10, 0.1
                )
            )
            for step in (0, 10):
                pressure_payload.extend(
                    struct.pack("<Qdd", step, 1000.0 + step, 100.0 + step)
                )
            values.write_bytes(pressure_payload)

            def audit(path: Path, relative: str) -> dict[str, object]:
                return {
                    "path": relative,
                    "size_bytes": path.stat().st_size,
                    "sha256": runner.file_sha256(path),
                }

            with patch.object(runner, "CASE_DIR", root), patch.object(
                runner.study, "validate_config"
            ):
                pressure_header = runner.pressure_database_header(values)
                sentinel = {
                    "schema": runner.COMPLETION_SCHEMA,
                    "config": {
                        "path": "configs/screen/replay.json",
                        "sha256": runner.file_sha256(config_path),
                        "uuid": "DYNAMIC",
                        "nsteps": 10,
                        "validation_profile": "registered-study-v1",
                    },
                    "artifacts": {
                        "final_vtp": audit(
                            final_vtp, "results/screen/DYNAMIC/particle10.vtp"
                        )
                    },
                    "runtime_dependencies": {
                        "resume_equilibrium": {
                            "checkpoint_hdf5": audit(
                                checkpoint, "results/screen/EQ/particles10.h5"
                            ),
                            "qa_vtp": audit(
                                qa_vtp, "results/screen/EQ/particle10.vtp"
                            ),
                        },
                        "read_pressure_database": {
                            "points": audit(
                                points,
                                "pressure_databases/screen/lagged/pressure_points.txt",
                            ),
                            "values": audit(
                                values,
                                "pressure_databases/screen/lagged/pressure_values.bin",
                            ),
                            "header": pressure_header,
                        },
                    },
                }
                (dynamic / runner.COMPLETION_FILENAME).write_text(
                    json.dumps(sentinel), encoding="utf-8"
                )
                self.assertTrue(runner.config_complete(config_path))

                original_checkpoint = checkpoint.read_bytes()
                checkpoint.write_bytes(b"changed-checkpoint")
                self.assertFalse(runner.config_complete(config_path))
                checkpoint.write_bytes(original_checkpoint)
                self.assertTrue(runner.config_complete(config_path))

                changed_pressure = bytearray(values.read_bytes())
                changed_pressure[-1] ^= 1
                values.write_bytes(changed_pressure)
                self.assertFalse(runner.config_complete(config_path))

    def test_phase_completion_checks_current_source_provenance(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            group_root = root / "pressure_databases" / "screen"
            source = group_root / "lagged"
            output = group_root / "phase_erased"
            source.mkdir(parents=True)
            output.mkdir(parents=True)
            source_points = source / "pressure_points.txt"
            source_values = source / "pressure_values.bin"
            output_points = output / "phase_erased_points.txt"
            output_values = output / "phase_erased_values.bin"
            source_points.write_text("0 0.0 0.0\n", encoding="utf-8")
            output_points.write_text("0 0.0 0.0\n", encoding="utf-8")
            payload = bytearray(
                runner.PRESSURE_HEADER.pack(
                    b"MPM_PRESSURE_V2\0", 2, 2, 1, 10, 10, 0.1
                )
            )
            for step in (0, 10):
                payload.extend(struct.pack("<Qdd", step, 1000.0 + step, 100.0 + step))
            source_values.write_bytes(payload)
            output_values.write_bytes(payload)
            metadata = {
                "schema": "phase-erased-pressure-control-v1",
                "classification": "one-way numerical counterfactual; not fully coupled",
                "source": {
                    "points": "pressure_databases/screen/lagged/pressure_points.txt",
                    "values": "pressure_databases/screen/lagged/pressure_values.bin",
                    "points_sha256": runner.file_sha256(source_points),
                    "values_sha256": runner.file_sha256(source_values),
                },
                "output": {
                    "points": "pressure_databases/screen/phase_erased/phase_erased_points.txt",
                    "values": "pressure_databases/screen/phase_erased/phase_erased_values.bin",
                    "points_sha256": runner.file_sha256(output_points),
                    "values_sha256": runner.file_sha256(output_values),
                },
                "header": {
                    "format_version": "V2",
                    "dimension": 2,
                    "particle_count": 1,
                    "step_interval": 10,
                    "max_step": 10,
                    "source_dt_s": 0.1,
                },
                "transform": {
                    "period_s": 1.3,
                    "fit_start_s": 2.6,
                    "fit_end_s": 10.4,
                    "ramp_time_s": 1.3,
                    "surface_band_m": 0.015,
                    "minimum_reference_amplitude_pa": 1.0,
                },
            }
            (output / "phase_erased_metadata.json").write_text(
                json.dumps(metadata), encoding="utf-8"
            )
            with patch.object(runner, "CASE_DIR", root):
                self.assertTrue(runner.phase_control_complete("screen", 10.4))
                corrupted = bytearray(source_values.read_bytes())
                corrupted[-1] ^= 1
                source_values.write_bytes(corrupted)
                self.assertFalse(runner.phase_control_complete("screen", 10.4))


if __name__ == "__main__":
    unittest.main()
