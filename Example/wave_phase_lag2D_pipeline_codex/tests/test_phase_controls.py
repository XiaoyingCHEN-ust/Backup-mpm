from __future__ import annotations

import math
import struct
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np


CASE_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(CASE_DIR))

import phase_controls as controls  # noqa: E402


class PhaseControlTest(unittest.TestCase):
    def make_database(self, directory: Path) -> tuple[Path, Path]:
        points, values = controls.database_paths(directory, "pressure")
        points.write_text(
            "# particle_count 4\n"
            "# id x0 x1\n"
            "0 0.0 0.5\n"
            "1 0.0 0.4\n"
            "2 1.0 0.5\n"
            "3 1.0 0.4\n",
            encoding="utf-8",
        )
        header = controls.PressureHeader(2, 2, 4, 1, 39, 0.1)
        omega = 2.0 * math.pi
        surface_phase = np.asarray([0.0, 0.0, 0.6, 0.6])
        lag = np.asarray([0.0, 0.4, 0.0, 0.4])
        amplitude = np.asarray([10.0, 7.0, 12.0, 8.0])
        with values.open("wb") as stream:
            controls.write_header(stream, header)
            for step in range(40):
                time = step * header.source_dt
                phase = surface_phase - lag
                liquid = 100.0 + amplitude * np.cos(omega * time + phase)
                gas = 20.0 + 0.5 * amplitude * np.cos(omega * time + phase)
                stream.write(struct.pack("<Q", step))
                stream.write(np.asarray(liquid, dtype="<f8").tobytes())
                stream.write(np.asarray(gas, dtype="<f8").tobytes())
        return points, values

    def test_transform_preserves_mean_and_amplitude_and_erases_depth_lag(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            _, source_values = self.make_database(directory)
            metadata = controls.transform_database(
                directory,
                "pressure",
                directory / "control",
                "phase_erased",
                period=1.0,
                fit_start=0.0,
                fit_end=3.9,
                ramp_time=0.0,
                surface_band=0.01,
                minimum_reference_amplitude=0.01,
            )
            self.assertTrue(metadata.exists())

            with source_values.open("rb") as stream:
                source_header = controls.read_header(stream)
            transformed_values = directory / "control" / "phase_erased_values.bin"
            original = controls.fit_harmonic(
                source_values, source_header, 1.0, 0.0, 3.9
            )
            transformed = controls.fit_harmonic(
                transformed_values, source_header, 1.0, 0.0, 3.9
            )
            np.testing.assert_allclose(transformed[0], original[0], atol=1.0e-11)
            np.testing.assert_allclose(
                np.hypot(transformed[1], transformed[2]),
                np.hypot(original[1], original[2]),
                atol=1.0e-11,
            )
            for phase in range(2):
                for surface, depth in ((0, 1), (2, 3)):
                    cross = (
                        transformed[1, surface, phase]
                        * transformed[2, depth, phase]
                        - transformed[2, surface, phase]
                        * transformed[1, depth, phase]
                    )
                    self.assertAlmostEqual(float(cross), 0.0, places=10)

    def test_rejects_incomplete_database(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            _, values = self.make_database(directory)
            values.write_bytes(values.read_bytes()[:-1])
            with values.open("rb") as stream:
                header = controls.read_header(stream)
            with self.assertRaisesRegex(ValueError, "incomplete"):
                controls.validate_database(values, header)

    def test_v1_frame_zero_has_legacy_one_step_time_offset(self):
        v1 = controls.PressureHeader(2, 2, 4, 5, 40, 0.02, "V1")
        v2 = controls.PressureHeader(2, 2, 4, 5, 40, 0.02, "V2")
        self.assertAlmostEqual(v1.physical_time(0), 0.02)
        self.assertAlmostEqual(v1.physical_time(5), 0.12)
        self.assertAlmostEqual(v2.physical_time(0), 0.0)
        self.assertAlmostEqual(v2.physical_time(5), 0.10)


if __name__ == "__main__":
    unittest.main()
