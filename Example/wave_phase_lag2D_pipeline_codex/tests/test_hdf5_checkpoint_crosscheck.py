from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np


CASE_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(CASE_DIR))

import hdf5_checkpoint_crosscheck as check  # noqa: E402


class Hdf5CheckpointCrosscheckTest(unittest.TestCase):
    def setUp(self) -> None:
        self.summary = {
            "schema": "mpm-hdf5-checkpoint-audit-v1",
            "passed": True,
            "hdf5": "",
            "table": "table",
            "particle_count": 2,
            "table_fields": 173,
            "minimum_id": 0,
            "maximum_id": 1,
            "ids_contiguous_unique": True,
            "formal_schema": True,
            "legacy_schema": False,
            "compared_fields": [
                "coordinates",
                "velocities",
                "displacements",
                "stresses",
                "porosities",
                "volumes",
                "status",
            ],
        }
        self.fields = {
            "ids": np.asarray([0, 1], dtype=np.int64),
            "coordinates": np.asarray([[0.0, 0.2, 0.0], [0.1, 0.3, 0.0]]),
            "velocities": np.asarray([[0.0, 0.0, 0.0], [3e-4, -4e-4, 0.0]]),
            "displacements": np.asarray([[0.0, 0.0, 0.0], [1e-3, -2e-3, 0.0]]),
            "stresses": np.asarray(
                [[-1000.0, -2000.0, -1500.0, 0.0, 0.0, 0.0],
                 [-1001.0, -2001.0, -1501.0, 25.0, 0.0, 0.0]]
            ),
            "porosities": np.asarray([[0.45], [0.46]]),
            "volumes": np.asarray([[1e-4], [1.01e-4]]),
        }

    def run_check(self, points: np.ndarray, arrays: dict[str, np.ndarray]):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            hdf5 = root / "particles10.h5"
            vtp = root / "particle10.vtp"
            auditor = root / "hdf5_checkpoint_audit"
            for path in (hdf5, vtp, auditor):
                path.write_bytes(b"fixture")
            summary = dict(self.summary)
            summary["hdf5"] = str(hdf5.resolve())
            with patch.object(
                check, "read_hdf5_jsonl", return_value=(summary, self.fields)
            ), patch.object(check, "read_vtp", return_value=(points, arrays)):
                return check.crosscheck(hdf5, vtp, auditor)

    def test_identical_fields_pass_even_when_vtp_rows_are_shuffled(self):
        order = np.asarray([1, 0])
        arrays = {"ids": self.fields["ids"][order]}
        for name in (
            "velocities",
            "displacements",
            "stresses",
            "porosities",
            "volumes",
        ):
            arrays[name] = self.fields[name][order]
        result = self.run_check(self.fields["coordinates"][order], arrays)
        self.assertTrue(result["passed"])
        self.assertEqual(result["particle_count_hdf5"], 2)
        self.assertEqual(result["ids_sha256_hdf5"], result["ids_sha256_vtp"])
        self.assertEqual(result["maximum_absolute_differences"]["stresses"], 0.0)

    def test_stale_vtp_field_fails(self):
        arrays = {"ids": self.fields["ids"]}
        for name in (
            "velocities",
            "displacements",
            "stresses",
            "porosities",
            "volumes",
        ):
            arrays[name] = self.fields[name].copy()
        arrays["stresses"][1, 0] += 1.0
        with self.assertRaisesRegex(ValueError, "stresses difference"):
            self.run_check(self.fields["coordinates"], arrays)

    def test_count_mismatch_fails(self):
        arrays = {
            "ids": np.asarray([0]),
            "velocities": self.fields["velocities"][:1],
            "displacements": self.fields["displacements"][:1],
            "stresses": self.fields["stresses"][:1],
            "porosities": self.fields["porosities"][:1],
            "volumes": self.fields["volumes"][:1],
        }
        with self.assertRaisesRegex(ValueError, "particle counts differ"):
            self.run_check(self.fields["coordinates"][:1], arrays)


if __name__ == "__main__":
    unittest.main()
