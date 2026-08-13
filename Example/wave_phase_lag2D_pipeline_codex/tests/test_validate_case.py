from __future__ import annotations

import json
import os
import struct
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch


CASE_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(CASE_DIR))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import prepare_study as study  # noqa: E402
import qa_contract as qa  # noqa: E402
import validate_case as validate  # noqa: E402
from phase_controls import FRAME_STEP, PressureHeader, write_header  # noqa: E402
from qa_test_helpers import (  # noqa: E402
    linear_equilibrium_config,
    write_fake_hdf5_auditor,
    write_linear_qa_fixture,
    write_mc_qa_fixture,
    write_pipeline_history_fixture,
)


MINIMAL_VTP = """<?xml version="1.0"?>
<VTKFile type="PolyData" version="0.1" byte_order="LittleEndian">
  <PolyData>
    <Piece NumberOfPoints="1" NumberOfVerts="0" NumberOfLines="0" NumberOfStrips="0" NumberOfPolys="0">
      <PointData />
      <Points><DataArray type="Float64" NumberOfComponents="3" format="ascii">0 0 0</DataArray></Points>
    </Piece>
  </PolyData>
</VTKFile>
"""

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

FULL_FINAL_VTP = """<?xml version="1.0"?>
<VTKFile type="PolyData" version="0.1" byte_order="LittleEndian">
  <PolyData>
    <Piece NumberOfPoints="1" NumberOfVerts="0" NumberOfLines="0" NumberOfStrips="0" NumberOfPolys="0">
      <PointData>
        <DataArray type="Int64" Name="ids" format="ascii">0</DataArray>
        <DataArray type="Float64" Name="porosities" format="ascii">0.485</DataArray>
        <DataArray type="Float64" Name="volumes" format="ascii">0.0001</DataArray>
        <DataArray type="Float64" Name="displacements" NumberOfComponents="3" format="ascii">0 0 0</DataArray>
        <DataArray type="Float64" Name="velocities" NumberOfComponents="3" format="ascii">0 0 0</DataArray>
        <DataArray type="Float64" Name="stresses" NumberOfComponents="6" format="ascii">0 0 0 0 0 0</DataArray>
      </PointData>
      <Points><DataArray type="Float64" NumberOfComponents="3" format="ascii">0.5 0.5 0</DataArray></Points>
    </Piece>
  </PolyData>
</VTKFile>
"""


class ValidateCaseTest(unittest.TestCase):
    def test_final_vtp_rejects_particles_outside_mesh(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            mesh = root / "mesh.txt"
            mesh.write_text(
                "4 1\n0 0\n1 0\n1 1\n0 1\n0 1 2 3\n", encoding="utf-8"
            )
            particles = root / "particles.txt"
            particles.write_text("1\n0.5 0.5\n", encoding="utf-8")
            vtp = root / "particle1.vtp"
            vtp.write_text(FULL_FINAL_VTP, encoding="utf-8")
            config = {
                "mesh": {"mesh": "mesh.txt"},
                "particles": [{"generator": {"location": "particles.txt"}}],
            }
            with patch.object(validate, "CASE_ROOT", root):
                validate.validate_final_vtp(vtp, config)
                vtp.write_text(
                    FULL_FINAL_VTP.replace(">0.5 0.5 0<", ">1.1 0.5 0<"),
                    encoding="utf-8",
                )
                with self.assertRaisesRegex(ValueError, "outside the mesh"):
                    validate.validate_final_vtp(vtp, config)

    def test_ascii_scalar_reader_rejects_legacy_id_value_rows(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "pressures.txt"
            path.write_text("2\n0 9810\n1 9712\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "2 columns; expected 1"):
                validate.read_ascii_table(path, 1)

    def test_output_base_is_created_recursively_and_cannot_escape_case(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            config_path = root / "configs" / "case.json"
            config_path.parent.mkdir(parents=True)
            config = {
                "analysis": {"uuid": "TEST", "nsteps": 1},
                "post_processing": {"path": "results/screen/"},
            }
            config_path.write_text(json.dumps(config), encoding="utf-8")

            with patch.object(validate, "CASE_ROOT", root):
                output_base = validate.prepare_output_base(config_path)
                self.assertEqual(output_base, (root / "results" / "screen").resolve())
                self.assertTrue(output_base.is_dir())

                config["post_processing"]["path"] = "../outside/"
                config_path.write_text(json.dumps(config), encoding="utf-8")
                with self.assertRaisesRegex(ValueError, "escapes the case directory"):
                    validate.prepare_output_base(config_path)

    def test_checkpoint_filename_uses_solver_padding(self):
        config = study.equilibrium_config(
            "TEST", 0.94, "results/test/", ".", 0.02, 0.01, 0.12
        )
        config["analysis"]["resume"]["resume"] = True
        config["analysis"]["resume"]["uuid"] = "EQ"
        path = validate.checkpoint_path(config)
        self.assertEqual(path.name, "particles40000.h5")
        self.assertEqual(path.parent.name, "EQ")

    def test_generated_manifest_selects_ten_cases(self):
        with tempfile.TemporaryDirectory() as temporary:
            manifest = study.generate_tier(
                "screen",
                Path(temporary),
                mesh_directory=".",
                cell_size=0.02,
                particle_spacing=0.01,
                wave_height=0.12,
                label="baseline",
            )
            data = json.loads(manifest.read_text(encoding="utf-8"))
            self.assertEqual(len(data["cases"]), 10)

    def test_completion_is_atomic_and_audits_vtp_hdf5_and_config(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            config_path = root / "configs" / "case.json"
            config_path.parent.mkdir(parents=True)
            config = linear_equilibrium_config(uuid="TEST_EQ")
            config_path.write_text(json.dumps(config), encoding="utf-8")
            result = root / "results" / "test" / "TEST_EQ"
            write_linear_qa_fixture(
                root, config, result / "particles10.h5"
            )
            auditor = write_fake_hdf5_auditor(root)

            with patch.dict(os.environ, {"MPM_HDF5_AUDITOR": str(auditor)}), patch.object(
                validate, "CASE_ROOT", root
            ):
                sentinel = validate.write_case_completion(config_path)
                payload = json.loads(sentinel.read_text(encoding="utf-8"))
                self.assertEqual(payload["schema"], validate.COMPLETION_SCHEMA)
                self.assertEqual(payload["config"]["uuid"], "TEST_EQ")
                self.assertEqual(
                    payload["artifacts"]["particle_vtp_grid"]["steps"], [10]
                )
                self.assertEqual(
                    len(payload["artifacts"]["particle_vtp_grid"]["files"]), 1
                )
                self.assertIn("pipeline_history_csv", payload["artifacts"])
                self.assertIn("final_hdf5", payload["artifacts"])
                self.assertTrue(
                    payload["artifacts"]["hdf5_vtp_crosscheck"]["passed"]
                )
                self.assertIn(
                    "max|v|=0.000e+00", payload["stability_qa"]["summary"]
                )
                self.assertEqual(
                    payload["stability_qa"]["mode"],
                    qa.LINEAR_EQUILIBRIUM_MODE,
                )
                self.assertFalse(list(result.glob(".*.tmp")))

                unrelated = result / "keep-me.txt"
                unrelated.write_text("unrelated", encoding="utf-8")
                cleared = validate.clear_case_completion(config_path)
                self.assertEqual(cleared, sentinel)
                self.assertFalse(sentinel.exists())
                self.assertTrue(unrelated.is_file())

    def test_completion_rejects_unstable_gated_checkpoint(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            config_path = root / "configs" / "case.json"
            config_path.parent.mkdir(parents=True)
            config = linear_equilibrium_config(uuid="UNSTABLE_EQ")
            config_path.write_text(json.dumps(config), encoding="utf-8")
            result = root / "results" / "test" / "UNSTABLE_EQ"
            write_linear_qa_fixture(
                root,
                config,
                result / "particles10.h5",
                maximum_velocity=0.002,
            )
            auditor = write_fake_hdf5_auditor(root)

            with patch.dict(os.environ, {"MPM_HDF5_AUDITOR": str(auditor)}), patch.object(
                validate, "CASE_ROOT", root
            ):
                with self.assertRaisesRegex(ValueError, "exceeds 1e-3"):
                    validate.write_case_completion(config_path)
            self.assertFalse((result / validate.COMPLETION_FILENAME).exists())

    def test_completion_rejects_nonfinite_gated_checkpoint(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            config_path = root / "configs" / "case.json"
            config_path.parent.mkdir(parents=True)
            config = linear_equilibrium_config(uuid="NONFINITE_EQ")
            config_path.write_text(json.dumps(config), encoding="utf-8")
            result = root / "results" / "test" / "NONFINITE_EQ"
            write_linear_qa_fixture(
                root,
                config,
                result / "particles10.h5",
                maximum_velocity=float("nan"),
            )
            auditor = write_fake_hdf5_auditor(root)

            with patch.dict(os.environ, {"MPM_HDF5_AUDITOR": str(auditor)}), patch.object(
                validate, "CASE_ROOT", root
            ):
                with self.assertRaisesRegex(ValueError, "non-finite|NaN/Inf"):
                    validate.write_case_completion(config_path)
            self.assertFalse((result / validate.COMPLETION_FILENAME).exists())

    def test_completion_rejects_incomplete_written_pressure_database(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            config_path = root / "configs" / "driver.json"
            config_path.parent.mkdir(parents=True)
            config = {
                "analysis": {
                    "uuid": "TEST_HD",
                    "nsteps": 10,
                    "prescribed_phase_pressures": {
                        "enable": False,
                        "write": True,
                        "path": "pressure_databases/test/lagged",
                        "file_prefix": "pressure",
                        "source_dt": 0.1,
                        "step_interval": 10,
                        "max_step": 10,
                    },
                },
                "post_processing": {
                    "path": "results/test/",
                    "write_hdf5": False,
                    "output_steps": 10,
                },
            }
            config_path.write_text(json.dumps(config), encoding="utf-8")
            result = root / "results" / "test" / "TEST_HD"
            result.mkdir(parents=True)
            (result / "particle10.vtp").write_text(MINIMAL_VTP, encoding="utf-8")
            write_pipeline_history_fixture(config, result)
            database = root / "pressure_databases" / "test" / "lagged"
            database.mkdir(parents=True)
            (database / "pressure_points.txt").write_text("0 0.0 0.0\n", encoding="utf-8")
            values = database / "pressure_values.bin"
            header = PressureHeader(2, 2, 1, 10, 10, 0.1)
            with values.open("wb") as stream:
                write_header(stream, header)
                for step in (0, 10):
                    stream.write(FRAME_STEP.pack(step))
                    stream.write(struct.pack("<dd", 1000.0 + step, 100.0 + step))

            with patch.object(validate, "CASE_ROOT", root):
                sentinel = validate.write_case_completion(config_path)
                payload = json.loads(sentinel.read_text(encoding="utf-8"))
                self.assertIn(
                    "written_pressure_database", payload["artifacts"]
                )
                values.write_bytes(values.read_bytes()[:-1])
                with self.assertRaisesRegex(ValueError, "truncated|incomplete"):
                    validate.write_case_completion(config_path)

    def test_completion_records_resume_and_read_pressure_dependencies(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            config_path = root / "configs" / "replay.json"
            config_path.parent.mkdir(parents=True)
            config = linear_equilibrium_config(uuid="REPLAY")
            config["analysis"].pop("stability_gate")
            config["analysis"].pop("stability_qa_contract")
            config["analysis"]["resume"] = {
                "resume": True,
                "uuid": "SOURCE_EQ",
                "step": 10,
                "nsteps": 10,
            }
            config["analysis"]["resume_stability_qa_contract"] = (
                qa.stability_contract(qa.LINEAR_EQUILIBRIUM_MODE)
            )
            config["analysis"]["prescribed_phase_pressures"] = {
                "enable": True,
                "write": False,
                "path": "pressure_databases/test/lagged",
                "file_prefix": "pressure",
                "source_dt": 0.1,
                "step_interval": 10,
                "max_step": 10,
            }
            config["post_processing"]["write_hdf5"] = False
            config_path.write_text(json.dumps(config), encoding="utf-8")
            equilibrium = root / "results" / "test" / "SOURCE_EQ"
            replay = root / "results" / "test" / "REPLAY"
            source = write_linear_qa_fixture(
                root, config, equilibrium / "particles10.h5"
            )
            replay.mkdir(parents=True)
            (replay / "particle10.vtp").write_text(
                source["vtp"].read_text(encoding="utf-8"), encoding="utf-8"
            )
            write_pipeline_history_fixture(config, replay)
            auditor = write_fake_hdf5_auditor(root)

            database = root / "pressure_databases" / "test" / "lagged"
            database.mkdir(parents=True)
            (database / "pressure_points.txt").write_text(
                "0 0.0 0.0\n", encoding="utf-8"
            )
            values = database / "pressure_values.bin"
            header = PressureHeader(2, 2, 1, 10, 10, 0.1)
            with values.open("wb") as stream:
                write_header(stream, header)
                for step in (0, 10):
                    stream.write(FRAME_STEP.pack(step))
                    stream.write(struct.pack("<dd", 1000.0 + step, 100.0 + step))

            with patch.dict(os.environ, {"MPM_HDF5_AUDITOR": str(auditor)}), patch.object(
                validate, "CASE_ROOT", root
            ):
                sentinel = validate.write_case_completion(config_path)
            payload = json.loads(sentinel.read_text(encoding="utf-8"))
            dependencies = payload["runtime_dependencies"]
            self.assertIn("resume_equilibrium", dependencies)
            self.assertIn("read_pressure_database", dependencies)
            self.assertIn(
                "checkpoint_hdf5", dependencies["resume_equilibrium"]
            )
            self.assertTrue(
                dependencies["resume_equilibrium"]["hdf5_vtp_crosscheck"][
                    "passed"
                ]
            )
            self.assertEqual(
                dependencies["resume_equilibrium"]["stability_qa"]["mode"],
                qa.LINEAR_EQUILIBRIUM_MODE,
            )
            self.assertEqual(
                dependencies["read_pressure_database"]["header"]["frame_count"],
                2,
            )

    def test_completion_fails_closed_without_hdf5_auditor(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            config_path = root / "configs" / "case.json"
            config_path.parent.mkdir(parents=True)
            config = linear_equilibrium_config(uuid="NO_AUDITOR_EQ")
            config_path.write_text(json.dumps(config), encoding="utf-8")
            result = root / "results" / "test" / "NO_AUDITOR_EQ"
            write_linear_qa_fixture(root, config, result / "particles10.h5")
            missing = root / "missing-hdf5-auditor"

            with patch.dict(
                os.environ, {"MPM_HDF5_AUDITOR": str(missing)}
            ), patch.object(validate, "CASE_ROOT", root):
                with self.assertRaisesRegex(FileNotFoundError, "auditor.*missing"):
                    validate.write_case_completion(config_path)
            self.assertFalse((result / validate.COMPLETION_FILENAME).exists())

    def test_normal_stress_rejects_old_positive_1p2_kpa_offset(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            config = linear_equilibrium_config(uuid="DRIFT_EQ")
            checkpoint = root / "results" / "test" / "DRIFT_EQ" / "particles10.h5"
            write_linear_qa_fixture(
                root, config, checkpoint, stress_offset_pa=1_200.0
            )
            with patch.object(validate, "CASE_ROOT", root):
                with self.assertRaisesRegex(ValueError, "gamma_eff\\*h"):
                    validate.validate_equilibrium_vtp(
                        checkpoint,
                        config,
                        mode=qa.LINEAR_EQUILIBRIUM_MODE,
                    )

    def test_normal_stress_rejects_more_than_five_percent_tension(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            config = linear_equilibrium_config(uuid="TENSION_EQ")
            checkpoint = root / "results" / "test" / "TENSION_EQ" / "particles10.h5"
            fixture = write_linear_qa_fixture(
                root, config, checkpoint, tensile_quiet_per_row=8
            )
            self.assertGreater(
                8 / len(fixture["quiet_ids_by_row"][0]),
                qa.NORMAL_STRESS_MAXIMUM_TENSILE_FRACTION,
            )
            with patch.object(validate, "CASE_ROOT", root):
                with self.assertRaisesRegex(ValueError, "tensile fraction"):
                    validate.validate_equilibrium_vtp(
                        checkpoint,
                        config,
                        mode=qa.LINEAR_EQUILIBRIUM_MODE,
                    )

    def test_normal_stress_aligns_shuffled_vtp_rows_by_particle_id(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            config = linear_equilibrium_config(uuid="SHUFFLED_EQ")
            checkpoint = root / "results" / "test" / "SHUFFLED_EQ" / "particles10.h5"
            write_linear_qa_fixture(
                root,
                config,
                checkpoint,
                row_spacing=0.05,
                shuffled_ids=True,
            )
            with patch.object(validate, "CASE_ROOT", root):
                record = validate.validate_equilibrium_vtp(
                    checkpoint,
                    config,
                    mode=qa.LINEAR_EQUILIBRIUM_MODE,
                )
            self.assertEqual(
                [row["median_drift_pa"] for row in record["normal_stress"]["rows"]],
                [0.0, 0.0],
            )

    def test_normal_stress_tolerance_uses_cellsize_not_particle_spacing(self):
        records = []
        counts = []
        for spacing in (0.01, 0.005):
            with tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                config = linear_equilibrium_config(
                    uuid="PPC_EQ", particle_spacing=spacing
                )
                checkpoint = root / "results" / "test" / "PPC_EQ" / "particles10.h5"
                fixture = write_linear_qa_fixture(
                    root,
                    config,
                    checkpoint,
                    row_spacing=spacing,
                    x_spacing=spacing,
                )
                with patch.object(validate, "CASE_ROOT", root):
                    records.append(
                        validate.validate_equilibrium_vtp(
                            checkpoint,
                            config,
                            mode=qa.LINEAR_EQUILIBRIUM_MODE,
                        )
                    )
                counts.append(fixture["particle_count"])
        tolerances = [
            record["normal_stress"]["median_drift_tolerance_pa"]
            for record in records
        ]
        self.assertNotEqual(counts[0], counts[1])
        self.assertAlmostEqual(tolerances[0], tolerances[1], places=12)
        self.assertAlmostEqual(
            tolerances[0],
            records[0]["normal_stress"]["effective_unit_weight_n_m3"] * 0.02,
            places=12,
        )

    def test_mohr_coulomb_handoff_residual_passes_and_fails_at_one_micro_pa(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            config = linear_equilibrium_config(uuid="UNIT_MC_RELAX")
            config["materials"][0]["type"] = "MohrCoulomb2D"
            config["materials"][0]["tension_cutoff"] = 0.0
            checkpoint = root / "results" / "test" / "UNIT_MC_RELAX" / "particles10.h5"
            write_mc_qa_fixture(root, config, checkpoint)
            with patch.object(validate, "CASE_ROOT", root):
                record = validate.validate_equilibrium_vtp(
                    checkpoint,
                    config,
                    mode=qa.MC_HANDOFF_MODE,
                )
            self.assertEqual(
                record["mc_feasibility"]["maximum_positive_residual_pa"], 0.0
            )

            write_mc_qa_fixture(
                root,
                config,
                checkpoint,
                tensile_stress_pa=2.0e-6,
            )
            with patch.object(validate, "CASE_ROOT", root):
                with self.assertRaisesRegex(ValueError, "positive yield residual"):
                    validate.validate_equilibrium_vtp(
                        checkpoint,
                        config,
                        mode=qa.MC_HANDOFF_MODE,
                    )

    def test_registered_velocity_limit_remains_one_millimetre_per_second(self):
        self.assertEqual(qa.MAXIMUM_VELOCITY_M_S, 1.0e-3)
        for mode in qa.STABILITY_MODES:
            self.assertEqual(
                qa.stability_contract(mode)["limits"]["maximum_velocity_m_s"],
                1.0e-3,
            )


if __name__ == "__main__":
    unittest.main()
