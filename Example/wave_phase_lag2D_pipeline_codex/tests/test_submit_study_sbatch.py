from __future__ import annotations

import json
import hashlib
import struct
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch


CASE_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(CASE_DIR))

import submit_study_sbatch as runner  # noqa: E402


def artifact_record(path: Path, root: Path) -> dict[str, object]:
    return {
        "path": path.resolve().relative_to(root.resolve()).as_posix(),
        "size_bytes": path.stat().st_size,
        "sha256": runner.file_sha256(path),
    }


def passing_hdf5_vtp_audit(
    hdf5: Path, vtp: Path, particle_count: int = 1
) -> dict[str, object]:
    identifier_hash = hashlib.sha256(
        b"".join(index.to_bytes(8, "little") for index in range(particle_count))
    ).hexdigest()
    return {
        "schema": runner.qa.HDF5_VTP_AUDIT_SCHEMA,
        "passed": True,
        "hdf5": str(hdf5.resolve()),
        "vtp": str(vtp.resolve()),
        "particle_count_hdf5": particle_count,
        "particle_count_vtp": particle_count,
        "ids_sha256_hdf5": identifier_hash,
        "ids_sha256_vtp": identifier_hash,
        "compared_fields": list(runner.qa.HDF5_VTP_COMPARED_FIELDS),
        "maximum_absolute_differences": {
            field: 0.0 for field in runner.qa.HDF5_VTP_COMPARED_FIELDS
        },
        "tolerances": dict(runner.qa.HDF5_VTP_ABSOLUTE_TOLERANCES),
        "hdf5_structure": {
            "schema": runner.qa.HDF5_STRUCTURE_AUDIT_SCHEMA,
            "passed": True,
            "hdf5": str(hdf5.resolve()),
            "table": "table",
            "table_fields": runner.qa.HDF5_FORMAL_TABLE_FIELDS,
            "particle_count": particle_count,
            "minimum_id": 0,
            "maximum_id": particle_count - 1,
            "ids_contiguous_unique": True,
            "formal_schema": True,
            "legacy_schema": False,
            "compared_fields": list(runner.qa.HDF5_STRUCTURE_COMPARED_FIELDS),
        },
    }


def passing_mc_qa(particle_count: int = 1) -> dict[str, object]:
    return {
        "schema": runner.qa.STABILITY_QA_SCHEMA,
        "mode": runner.qa.MC_HANDOFF_MODE,
        "limits": runner.qa.stability_contract(runner.qa.MC_HANDOFF_MODE)[
            "limits"
        ],
        "observed": {
            "maximum_velocity_m_s": 0.0,
            "maximum_displacement_m": 0.0,
            "porosity_min": 0.485,
            "porosity_max": 0.485,
            "particle_count": particle_count,
        },
        "mc_feasibility": {
            "particle_count": particle_count,
            "maximum_tension_residual_pa": 0.0,
            "maximum_shear_residual_pa": 0.0,
            "maximum_positive_residual_pa": 0.0,
            "violating_particle_count": 0,
        },
        "summary": "passing synthetic MC handoff QA",
    }


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
                "analysis": {"uuid": "TEST", "nsteps": 20},
                "post_processing": {
                    "path": "results/screen/",
                    "write_hdf5": False,
                    "output_steps": 10,
                },
            }
            config_path.write_text(json.dumps(config), encoding="utf-8")
            result = root / "results" / "screen" / "TEST"
            result.mkdir(parents=True)
            middle_vtp = result / "particle10.vtp"
            final_vtp = result / "particle20.vtp"
            middle_vtp.write_text("middle-vtp", encoding="utf-8")
            final_vtp.write_text("complete-vtp", encoding="utf-8")
            history = result / "pipeline-history00.csv"
            history.write_text("step,time\n0,0\n", encoding="utf-8")
            sentinel = {
                "schema": runner.COMPLETION_SCHEMA,
                "config": {
                    "path": "configs/screen/case.json",
                    "sha256": runner.file_sha256(config_path),
                    "uuid": "TEST",
                    "nsteps": 20,
                    "validation_profile": "registered-study-v1",
                },
                "artifacts": {
                    "particle_vtp_grid": {
                        "steps": [10, 20],
                        "files": [
                            artifact_record(middle_vtp, root),
                            artifact_record(final_vtp, root),
                        ],
                    },
                    "pipeline_history_csv": artifact_record(history, root),
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
                middle_original = middle_vtp.read_text(encoding="utf-8")
                middle_vtp.write_text("tampered-middle", encoding="utf-8")
                self.assertFalse(runner.config_complete(config_path))
                middle_vtp.write_text(middle_original, encoding="utf-8")
                self.assertTrue(runner.config_complete(config_path))

                history_original = history.read_text(encoding="utf-8")
                history.write_text("tampered-history", encoding="utf-8")
                self.assertFalse(runner.config_complete(config_path))
                history.write_text(history_original, encoding="utf-8")
                self.assertTrue(runner.config_complete(config_path))

                sentinel["schema"] = "pipeline-case-completion-v2"
                (result / runner.COMPLETION_FILENAME).write_text(
                    json.dumps(sentinel), encoding="utf-8"
                )
                self.assertFalse(runner.config_complete(config_path))
                sentinel["schema"] = runner.COMPLETION_SCHEMA
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
                        "uuid": "EQ_MC_RELAX",
                        "step": 10,
                        "nsteps": 10,
                    },
                    "resume_stability_qa_contract": runner.qa.stability_contract(
                        runner.qa.MC_HANDOFF_MODE
                    ),
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
                    "output_steps": 10,
                },
            }
            config_path.write_text(json.dumps(config), encoding="utf-8")
            dynamic = root / "results" / "screen" / "DYNAMIC"
            equilibrium = root / "results" / "screen" / "EQ_MC_RELAX"
            dynamic.mkdir(parents=True)
            equilibrium.mkdir(parents=True)
            final_vtp = dynamic / "particle10.vtp"
            checkpoint = equilibrium / "particles10.h5"
            qa_vtp = equilibrium / "particle10.vtp"
            final_vtp.write_text("dynamic-result", encoding="utf-8")
            history = dynamic / "pipeline-history00.csv"
            history.write_text("step,time\n0,0\n", encoding="utf-8")
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
                        "particle_vtp_grid": {
                            "steps": [10],
                            "files": [artifact_record(final_vtp, root)],
                        },
                        "pipeline_history_csv": artifact_record(history, root),
                    },
                    "runtime_dependencies": {
                        "resume_equilibrium": {
                            "checkpoint_hdf5": artifact_record(checkpoint, root),
                            "qa_vtp": artifact_record(qa_vtp, root),
                            "hdf5_vtp_crosscheck": passing_hdf5_vtp_audit(
                                checkpoint, qa_vtp
                            ),
                            "stability_qa": passing_mc_qa(),
                        },
                        "read_pressure_database": {
                            "points": artifact_record(points, root),
                            "values": artifact_record(values, root),
                            "header": pressure_header,
                        },
                    },
                }
                (dynamic / runner.COMPLETION_FILENAME).write_text(
                    json.dumps(sentinel), encoding="utf-8"
                )
                self.assertTrue(runner.config_complete(config_path))

                substitute = equilibrium / "substitute.h5"
                substitute.write_bytes(checkpoint.read_bytes())
                replaced_resume = json.loads(json.dumps(sentinel))
                replaced_resume["runtime_dependencies"]["resume_equilibrium"][
                    "checkpoint_hdf5"
                ] = artifact_record(substitute, root)
                (dynamic / runner.COMPLETION_FILENAME).write_text(
                    json.dumps(replaced_resume), encoding="utf-8"
                )
                self.assertFalse(runner.config_complete(config_path))

                missing_pressure = json.loads(json.dumps(sentinel))
                del missing_pressure["runtime_dependencies"][
                    "read_pressure_database"
                ]["header"]
                (dynamic / runner.COMPLETION_FILENAME).write_text(
                    json.dumps(missing_pressure), encoding="utf-8"
                )
                self.assertFalse(runner.config_complete(config_path))

                alternate_points = database / "alternate_points.txt"
                alternate_points.write_bytes(points.read_bytes())
                replaced_pressure = json.loads(json.dumps(sentinel))
                replaced_pressure["runtime_dependencies"][
                    "read_pressure_database"
                ]["points"] = artifact_record(alternate_points, root)
                (dynamic / runner.COMPLETION_FILENAME).write_text(
                    json.dumps(replaced_pressure), encoding="utf-8"
                )
                self.assertFalse(runner.config_complete(config_path))

                (dynamic / runner.COMPLETION_FILENAME).write_text(
                    json.dumps(sentinel), encoding="utf-8"
                )

                missing_qa = json.loads(json.dumps(sentinel))
                del missing_qa["runtime_dependencies"]["resume_equilibrium"][
                    "stability_qa"
                ]
                (dynamic / runner.COMPLETION_FILENAME).write_text(
                    json.dumps(missing_qa), encoding="utf-8"
                )
                self.assertFalse(runner.config_complete(config_path))
                (dynamic / runner.COMPLETION_FILENAME).write_text(
                    json.dumps(sentinel), encoding="utf-8"
                )

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
