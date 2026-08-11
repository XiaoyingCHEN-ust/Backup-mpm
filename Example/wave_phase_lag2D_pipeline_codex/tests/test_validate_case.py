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

import prepare_study as study  # noqa: E402
import validate_case as validate  # noqa: E402
from phase_controls import FRAME_STEP, PressureHeader, write_header  # noqa: E402


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


class ValidateCaseTest(unittest.TestCase):
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
        self.assertEqual(path.name, "particles10000.h5")
        self.assertEqual(path.parent.name, "EQ")

    def test_generated_manifest_selects_nine_cases(self):
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
            self.assertEqual(len(data["cases"]), 9)

    def test_completion_is_atomic_and_audits_vtp_hdf5_and_config(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            config_path = root / "configs" / "case.json"
            config_path.parent.mkdir(parents=True)
            config = {
                "analysis": {"uuid": "TEST_EQ", "nsteps": 10},
                "post_processing": {"path": "results/test/", "write_hdf5": True},
            }
            config_path.write_text(json.dumps(config), encoding="utf-8")
            result = root / "results" / "test" / "TEST_EQ"
            result.mkdir(parents=True)
            (result / "particle10.vtp").write_text(MINIMAL_VTP, encoding="utf-8")
            (result / "particles10.h5").write_bytes(
                validate.HDF5_SIGNATURE + b"completion-test"
            )

            with patch.object(validate, "CASE_ROOT", root):
                sentinel = validate.write_case_completion(config_path)
                payload = json.loads(sentinel.read_text(encoding="utf-8"))
                self.assertEqual(payload["schema"], validate.COMPLETION_SCHEMA)
                self.assertEqual(payload["config"]["uuid"], "TEST_EQ")
                self.assertIn("final_vtp", payload["artifacts"])
                self.assertIn("final_hdf5", payload["artifacts"])
                self.assertFalse(list(result.glob(".*.tmp")))

                unrelated = result / "keep-me.txt"
                unrelated.write_text("unrelated", encoding="utf-8")
                cleared = validate.clear_case_completion(config_path)
                self.assertEqual(cleared, sentinel)
                self.assertFalse(sentinel.exists())
                self.assertTrue(unrelated.is_file())

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
                "post_processing": {"path": "results/test/", "write_hdf5": False},
            }
            config_path.write_text(json.dumps(config), encoding="utf-8")
            result = root / "results" / "test" / "TEST_HD"
            result.mkdir(parents=True)
            (result / "particle10.vtp").write_text(MINIMAL_VTP, encoding="utf-8")
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
            config = {
                "analysis": {
                    "uuid": "REPLAY",
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
                        "path": "pressure_databases/test/lagged",
                        "file_prefix": "pressure",
                        "source_dt": 0.1,
                        "step_interval": 10,
                        "max_step": 10,
                    },
                },
                "post_processing": {"path": "results/test/", "write_hdf5": False},
            }
            config_path.write_text(json.dumps(config), encoding="utf-8")
            equilibrium = root / "results" / "test" / "EQ"
            replay = root / "results" / "test" / "REPLAY"
            equilibrium.mkdir(parents=True)
            replay.mkdir(parents=True)
            (equilibrium / "particles10.h5").write_bytes(
                validate.HDF5_SIGNATURE + b"equilibrium-checkpoint"
            )
            (equilibrium / "particle10.vtp").write_text(
                EQUILIBRIUM_VTP, encoding="utf-8"
            )
            (replay / "particle10.vtp").write_text(MINIMAL_VTP, encoding="utf-8")

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

            with patch.object(validate, "CASE_ROOT", root):
                sentinel = validate.write_case_completion(config_path)
            payload = json.loads(sentinel.read_text(encoding="utf-8"))
            dependencies = payload["runtime_dependencies"]
            self.assertIn("resume_equilibrium", dependencies)
            self.assertIn("read_pressure_database", dependencies)
            self.assertIn(
                "checkpoint_hdf5", dependencies["resume_equilibrium"]
            )
            self.assertEqual(
                dependencies["read_pressure_database"]["header"]["frame_count"],
                2,
            )


if __name__ == "__main__":
    unittest.main()
