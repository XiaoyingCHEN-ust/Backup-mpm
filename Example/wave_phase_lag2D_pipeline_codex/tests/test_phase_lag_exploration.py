from __future__ import annotations

import copy
import hashlib
import json
import math
import shutil
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import numpy as np

import prepare_phase_lag_exploration as exploration
import analyze_phase_lag_replay_exploration as replay_analysis
import phase_controls


class PhaseLagExplorationTests(unittest.TestCase):
    def setUp(self) -> None:
        self.source_root = Path(__file__).resolve().parents[1]
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        (self.root / "configs/screen").mkdir(parents=True)
        for relative in (
            "configs/screen/01_EQ_HS.json",
            "configs/screen/02_HS.json",
            "particles.txt",
            "top_surface_traction_particle_id.txt",
        ):
            destination = self.root / relative
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(self.source_root / relative, destination)

    def tearDown(self) -> None:
        self.temporary.cleanup()

    @staticmethod
    def read(path: Path) -> dict:
        return json.loads(path.read_text(encoding="utf-8"))

    @staticmethod
    def digest(path: Path) -> str:
        return hashlib.sha256(path.read_bytes()).hexdigest()

    @staticmethod
    def valid_phase_alignment_metadata() -> dict:
        record = {
            "surface_reference_eligible_particle_count": 12,
            "phase_defined_particle_count": 10,
            "changed_particle_count": 7,
            "max_abs_source_to_surface_phase_difference_deg": 18.0,
            "max_abs_residual_to_surface_phase_difference_deg": 0.0,
        }
        return {
            "qa": {
                "phase_alignment": {
                    "liquid": copy.deepcopy(record),
                    "gas": copy.deepcopy(record),
                }
            }
        }

    def write_phase_metadata(
        self, label: str, *, surface_source: Path | None = None
    ) -> Path:
        self.generate_replay_configs(label=label)
        lagged = (
            self.root
            / "pressure_databases/phase_lag_exploratory"
            / label
            / "lagged"
        )
        erased = lagged.parent / "phase_erased"
        lagged.mkdir(parents=True)
        erased.mkdir(parents=True)
        files = {
            "source_points": lagged / "pressure_points.txt",
            "source_values": lagged / "pressure_values.bin",
            "output_points": erased / "phase_erased_points.txt",
            "output_values": erased / "phase_erased_values.bin",
        }
        for name, path in files.items():
            path.write_bytes(f"{name}\n".encode("ascii"))
        registered = self.root / "top_surface_traction_particle_id.txt"
        reference_ids = [
            int(value)
            for value in registered.read_text(encoding="utf-8").splitlines()[1:]
        ]
        source_path = surface_source or registered
        metadata = {
            "schema": "phase-erased-pressure-control-v1",
            "producer": {
                "path": str(Path(phase_controls.__file__).resolve()),
                "sha256": self.digest(Path(phase_controls.__file__).resolve()),
            },
            "header": {
                "format_version": "V2",
                "dimension": 2,
                "particle_count": 7376,
                "step_interval": 130,
                "max_step": 39000,
                "source_dt_s": 1.0e-4,
                "sample_dt_s": 0.013000000000000001,
            },
            "transform": {
                "period_s": 1.3,
                "fit_start_s": 2.6,
                "fit_end_s": 3.9,
                "ramp_time_s": 1.3,
                "surface_band_m": 0.015,
                "minimum_reference_amplitude_pa": 1.0,
                "surface_reference": {
                    "method": "registered_top_surface_particle_set",
                    "source": {
                        "path": str(source_path.resolve()),
                        "sha256": self.digest(source_path),
                    },
                    "reference_particle_count": 150,
                    "reference_particle_ids": reference_ids,
                    "reference_particle_ids_sha256": (
                        phase_controls.particle_ids_sha256(
                            np.asarray(reference_ids, dtype=np.uint64)
                        )
                    ),
                },
            },
            "source": {
                "points": str(files["source_points"].resolve()),
                "points_sha256": self.digest(files["source_points"]),
                "values": str(files["source_values"].resolve()),
                "values_sha256": self.digest(files["source_values"]),
            },
            "output": {
                "points": str(files["output_points"].resolve()),
                "points_sha256": self.digest(files["output_points"]),
                "values": str(files["output_values"].resolve()),
                "values_sha256": self.digest(files["output_values"]),
            },
            "qa": {
                "max_abs_mean_change_pa": 0.0,
                "max_abs_fundamental_amplitude_change_pa": 0.0,
                "max_instantaneous_pressure_change_pa": 10.0,
                **self.valid_phase_alignment_metadata()["qa"],
            },
        }
        path = erased / "phase_erased_metadata.json"
        path.write_text(json.dumps(metadata), encoding="utf-8")
        return path

    def generate_replay_configs(
        self,
        *,
        label: str = "unit_variable_parameters",
        permeability_m2: float = 3.0e-13,
        saturation: float = 0.937,
    ) -> tuple[dict, dict, dict]:
        manifest = exploration.generate(
            self.root,
            label=label,
            permeability_m2=permeability_m2,
            saturation=saturation,
            equilibrium_steps=5000,
            wave_steps=39000,
            pressure_smoothing=False,
            paired_replay=True,
        )
        directory = self.root / "configs/phase_lag_exploratory" / label
        return (
            manifest,
            self.read(directory / "04_RL.json"),
            self.read(directory / "05_RE.json"),
        )

    def test_generates_strict_unsmoothed_paired_replay(self) -> None:
        label = "unit_variable_parameters"
        manifest, lagged, erased = self.generate_replay_configs(
            label=label,
            permeability_m2=3.0e-13,
            saturation=0.937,
        )
        directory = self.root / "configs/phase_lag_exploratory" / label
        equilibrium = self.read(directory / "01_EQ.json")
        wave = self.read(directory / "02_WAVE.json")
        cases = {
            "HD": self.read(directory / "03_HD.json"),
            "RL": lagged,
            "RE": erased,
        }

        self.assertEqual(equilibrium["materials"][0]["type"], "LinearElastic2D")
        self.assertEqual(equilibrium["analysis"]["PIC"], 1.0)
        self.assertFalse(equilibrium["analysis"]["APIC"])
        self.assertEqual(
            equilibrium["analysis"]["damping"]["damping_factor"], 5.0
        )
        self.assertTrue(equilibrium["analysis"]["rigid_pipeline"]["fixed"])
        self.assertFalse(equilibrium["materials"][1]["wave_pressure"])
        self.assertEqual(equilibrium["analysis"]["nsteps"], 5000)

        for config in (equilibrium, wave, *cases.values()):
            self.assertEqual(config["materials"][0]["intrinsic_permeability"], 3e-13)
            self.assertEqual(config["materials"][1]["liquid_saturation"], 0.937)
            self.assertEqual(
                config["materials"][1]["gas_saturation"], 1.0 - 0.937
            )
            self.assertIs(config["analysis"]["pressure_smoothing"], False)
            self.assertIs(config["analysis"]["pressure_smoothing_in_loop"], False)

        self.assertEqual(wave["materials"][0]["type"], "SANISAND2D")
        self.assertTrue(wave["materials"][1]["wave_pressure"])
        for role, config in cases.items():
            pressure = config["analysis"]["prescribed_phase_pressures"]
            self.assertTrue(config["analysis"]["rigid_pipeline"]["fixed"])
            self.assertEqual(pressure["step_interval"], 130)
            self.assertEqual(pressure["max_step"], 39000)
            self.assertEqual(pressure["mapping"], "id")
            self.assertEqual(config["analysis"]["nsteps"], 39000)
            if role == "HD":
                self.assertTrue(config["materials"][1]["wave_pressure"])
                self.assertTrue(pressure["write"])
                self.assertFalse(pressure["enable"])
            else:
                self.assertFalse(config["materials"][1]["wave_pressure"])
                self.assertFalse(pressure["write"])
                self.assertTrue(pressure["enable"])

        generated = manifest["generated"]
        self.assertEqual(
            generated["hd_sha256"], self.digest(directory / "03_HD.json")
        )
        self.assertEqual(
            generated["rl_sha256"], self.digest(directory / "04_RL.json")
        )
        self.assertEqual(
            generated["re_sha256"], self.digest(directory / "05_RE.json")
        )
        self.assertTrue(manifest["parameters"]["paired_replay"])
        self.assertEqual(
            manifest["parameters"]["pressure_database_interval_steps"], 130
        )
        parameters = replay_analysis.validate_replay_pair_configs(lagged, erased)
        self.assertEqual(parameters.intrinsic_permeability_m2, 3.0e-13)
        self.assertEqual(parameters.liquid_saturation, 0.937)
        self.assertEqual(
            parameters.as_dict(),
            {
                "intrinsic_permeability_m2": 3.0e-13,
                "liquid_saturation": 0.937,
                "gas_saturation": 1.0 - 0.937,
                "pressure_smoothing": False,
                "pressure_smoothing_in_loop": False,
                "pipeline_fixed": True,
            },
        )
        display = replay_analysis.parameter_label(parameters)
        self.assertIn("3.000e-13", display)
        self.assertIn("0.937", display)

    def test_replay_pair_rejects_parameter_and_physics_mismatches(self) -> None:
        _, lagged, erased = self.generate_replay_configs()
        mismatches = (
            ("permeability", ("materials", 0, "intrinsic_permeability"), 4.0e-13),
            ("time step", ("analysis", "dt"), 2.0e-4),
        )
        for name, keys, value in mismatches:
            changed = copy.deepcopy(erased)
            target = changed
            for key in keys[:-1]:
                target = target[key]
            target[keys[-1]] = value
            with self.subTest(name=name), self.assertRaises(ValueError):
                replay_analysis.validate_replay_pair_configs(lagged, changed)

    def test_replay_parameters_require_explicit_no_smoothing(self) -> None:
        _, lagged, _ = self.generate_replay_configs()
        for key, value in (
            ("pressure_smoothing", True),
            ("pressure_smoothing", 0),
            ("pressure_smoothing", None),
            ("pressure_smoothing_in_loop", True),
            ("pressure_smoothing_in_loop", 0),
            ("pressure_smoothing_in_loop", None),
        ):
            changed = copy.deepcopy(lagged)
            changed["analysis"][key] = value
            with self.subTest(key=key, value=value), self.assertRaisesRegex(
                ValueError, "pressure smoothing must be explicitly disabled"
            ):
                replay_analysis.replay_parameters(changed, "lagged")

    def test_phase_metadata_binds_registered_surface_reference_artifact(self) -> None:
        label = "unit_phase_metadata"
        self.write_phase_metadata(label)
        audit = replay_analysis.validate_phase_metadata(self.root, label)
        registered = (self.root / "top_surface_traction_particle_id.txt").resolve()
        self.assertEqual(
            audit["surface_reference_source"],
            {
                "path": str(registered),
                "size_bytes": registered.stat().st_size,
                "sha256": self.digest(registered),
            },
        )

    def test_phase_metadata_rejects_arbitrary_150_id_surface_file(self) -> None:
        arbitrary = self.root / "arbitrary_150_particle_ids.txt"
        arbitrary.write_text(
            "150\n" + "\n".join(str(value) for value in range(150)) + "\n",
            encoding="utf-8",
        )
        label = "unit_arbitrary_surface"
        self.write_phase_metadata(label, surface_source=arbitrary)
        with self.assertRaisesRegex(
            ValueError, "not the registered top-surface particle file"
        ):
            replay_analysis.validate_phase_metadata(self.root, label)

    def test_phase_metadata_binds_source_and_output_to_replay_configs(self) -> None:
        label = "unit_database_binding"
        metadata_path = self.write_phase_metadata(label)
        metadata = self.read(metadata_path)
        substitute = self.root / "substitute_points.txt"
        substitute.write_text("substitute\n", encoding="utf-8")
        metadata["source"]["points"] = str(substitute.resolve())
        metadata["source"]["points_sha256"] = self.digest(substitute)
        metadata_path.write_text(json.dumps(metadata), encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "selected by the replay config"):
            replay_analysis.validate_phase_metadata(self.root, label)

    def test_phase_metadata_locks_every_numeric_transform_parameter(self) -> None:
        for name, replacement in (
            ("period_s", 1.4),
            ("fit_start_s", 2.5),
            ("fit_end_s", 4.0),
            ("ramp_time_s", 1.2),
            ("surface_band_m", 0.02),
            ("minimum_reference_amplitude_pa", 2.0),
        ):
            with self.subTest(name=name):
                label = f"unit_transform_{name}"
                metadata_path = self.write_phase_metadata(label)
                metadata = self.read(metadata_path)
                metadata["transform"][name] = replacement
                metadata_path.write_text(json.dumps(metadata), encoding="utf-8")
                with self.assertRaisesRegex(ValueError, rf"{name} changed"):
                    replay_analysis.validate_phase_metadata(self.root, label)

    def test_phase_metadata_rejects_nonfinite_or_zero_qa(self) -> None:
        cases = (
            ("max_abs_mean_change_pa", math.nan, "must be finite"),
            ("max_abs_fundamental_amplitude_change_pa", math.inf, "must be finite"),
            ("max_instantaneous_pressure_change_pa", math.nan, "must be finite"),
            ("max_instantaneous_pressure_change_pa", 0.0, "no positive pressure change"),
        )
        for index, (name, replacement, message) in enumerate(cases):
            with self.subTest(name=name, replacement=replacement):
                label = f"unit_qa_{index}"
                metadata_path = self.write_phase_metadata(label)
                metadata = self.read(metadata_path)
                metadata["qa"][name] = replacement
                metadata_path.write_text(json.dumps(metadata), encoding="utf-8")
                with self.assertRaisesRegex(ValueError, message):
                    replay_analysis.validate_phase_metadata(self.root, label)

    def test_phase_metadata_requires_current_phase_controls_producer(self) -> None:
        for mutation in ("missing", "hash"):
            with self.subTest(mutation=mutation):
                label = f"unit_producer_{mutation}"
                metadata_path = self.write_phase_metadata(label)
                metadata = self.read(metadata_path)
                if mutation == "missing":
                    metadata.pop("producer")
                    message = "producer binding is incomplete"
                else:
                    metadata["producer"]["sha256"] = "0" * 64
                    message = "producer binding is stale"
                metadata_path.write_text(json.dumps(metadata), encoding="utf-8")
                with self.assertRaisesRegex(ValueError, message):
                    replay_analysis.validate_phase_metadata(self.root, label)

    def test_atomic_json_rejects_nan_without_replacing_target(self) -> None:
        path = self.root / "atomic_metadata.json"
        phase_controls.atomic_write_json(path, {"valid": 1.0})
        original = path.read_bytes()
        with self.assertRaises(ValueError):
            phase_controls.atomic_write_json(path, {"invalid": math.nan})
        self.assertEqual(path.read_bytes(), original)
        self.assertEqual(list(path.parent.glob(f".{path.name}.*.tmp")), [])

    def test_transform_publication_replaces_data_then_metadata(self) -> None:
        output_points = self.root / "phase_erased_points.txt"
        output_values = self.root / "phase_erased_values.bin"
        metadata_path = self.root / "phase_erased_metadata.json"
        output_points.write_bytes(b"old points")
        output_values.write_bytes(b"old values")
        metadata_path.write_text('{"generation": "old"}\n', encoding="utf-8")

        candidate_points = self.root / "candidate_points.txt"
        candidate_values = self.root / "candidate_values.bin"
        candidate_points.write_bytes(b"new points")
        candidate_values.write_bytes(b"new values")
        staged_points = phase_controls._stage_file_copy(
            candidate_points, output_points
        )
        staged_values = phase_controls._stage_file_copy(
            candidate_values, output_values
        )

        real_replace = phase_controls.os.replace
        replacements: list[tuple[Path, Path]] = []

        def record_replace(source: Path, target: Path) -> None:
            replacements.append((Path(source), Path(target)))
            real_replace(source, target)

        with mock.patch.object(
            phase_controls.os, "replace", side_effect=record_replace
        ):
            phase_controls._publish_transform_files(
                staged_points=staged_points,
                staged_values=staged_values,
                output_points=output_points,
                output_values=output_values,
                metadata_path=metadata_path,
                metadata={"generation": "new"},
            )

        self.assertEqual(output_points.read_bytes(), b"new points")
        self.assertEqual(output_values.read_bytes(), b"new values")
        self.assertEqual(self.read(metadata_path), {"generation": "new"})
        self.assertEqual(
            [target for _, target in replacements],
            [output_values, output_points, metadata_path],
        )
        for source, target in replacements:
            self.assertEqual(source.parent.resolve(), target.parent.resolve())
            self.assertTrue(source.name.endswith(".tmp"))
        self.assertEqual(list(self.root.glob(".*.tmp")), [])

    def test_transform_publication_failure_preserves_old_trio(self) -> None:
        output_points = self.root / "phase_erased_points.txt"
        output_values = self.root / "phase_erased_values.bin"
        metadata_path = self.root / "phase_erased_metadata.json"
        old = {
            output_points: b"old points",
            output_values: b"old values",
            metadata_path: b'{"generation": "old"}\n',
        }
        old_modes: dict[Path, int] = {}
        for index, (path, content) in enumerate(old.items()):
            path.write_bytes(content)
            path.chmod(0o600 + index * 0o20)
            old_modes[path] = path.stat().st_mode & 0o7777

        candidate_points = self.root / "candidate_points.txt"
        candidate_values = self.root / "candidate_values.bin"
        candidate_points.write_bytes(b"new points")
        candidate_values.write_bytes(b"new values")
        staged_points = phase_controls._stage_file_copy(
            candidate_points, output_points
        )
        staged_values = phase_controls._stage_file_copy(
            candidate_values, output_values
        )

        real_replace = phase_controls.os.replace
        injected = False

        def fail_during_points_publish(source: Path, target: Path) -> None:
            nonlocal injected
            if Path(target) == output_points and not injected:
                injected = True
                raise OSError("simulated mid-publication failure")
            real_replace(source, target)

        with mock.patch.object(
            phase_controls.os,
            "replace",
            side_effect=fail_during_points_publish,
        ), self.assertRaisesRegex(OSError, "simulated mid-publication failure"):
            phase_controls._publish_transform_files(
                staged_points=staged_points,
                staged_values=staged_values,
                output_points=output_points,
                output_values=output_values,
                metadata_path=metadata_path,
                metadata={"generation": "new"},
            )

        for path, content in old.items():
            self.assertEqual(path.read_bytes(), content)
            self.assertEqual(path.stat().st_mode & 0o7777, old_modes[path])
        self.assertEqual(list(self.root.glob(".*.tmp")), [])

    def test_transform_rejects_source_output_path_aliases(self) -> None:
        parameters = {
            "period": 1.3,
            "fit_start": 0.0,
            "fit_end": 1.0,
            "ramp_time": 0.0,
            "surface_band": 0.1,
            "minimum_reference_amplitude": 1.0,
        }
        with self.assertRaisesRegex(ValueError, "source pressure database"):
            phase_controls.transform_database(
                self.root,
                "pressure",
                self.root,
                "pressure",
                **parameters,
            )

        source_points, source_values = phase_controls.database_paths(
            self.root, "pressure"
        )
        _, output_values = phase_controls.database_paths(self.root, "phase_erased")
        source_points.write_text("0 0 0\n", encoding="utf-8")
        source_values.write_bytes(b"invalid but alias-checked first")
        phase_controls.os.link(source_values, output_values)
        with self.assertRaisesRegex(ValueError, "source pressure database"):
            phase_controls.transform_database(
                self.root,
                "pressure",
                self.root,
                "phase_erased",
                **parameters,
            )

    def test_replay_frame_health_rejects_outside_or_nonpositive_volume(self) -> None:
        points = np.asarray([[0.25, 0.25], [0.75, 0.75]], dtype=float)
        arrays = {
            "volumes": np.asarray([0.1, 0.1]),
            "porosities": np.asarray([0.4, 0.4]),
        }
        lower = np.asarray([0.0, 0.0])
        upper = np.asarray([1.0, 1.0])
        volumes, porosity, outside = replay_analysis.validate_replay_frame_health(
            Path("frame.vtp"), points, arrays, lower, upper
        )
        np.testing.assert_array_equal(volumes, arrays["volumes"])
        np.testing.assert_array_equal(porosity, arrays["porosities"])
        self.assertEqual(outside, 0)

        outside_points = points.copy()
        outside_points[1, 0] = 1.01
        with self.assertRaisesRegex(ValueError, "outside the mesh"):
            replay_analysis.validate_replay_frame_health(
                Path("frame.vtp"), outside_points, arrays, lower, upper
            )

        invalid_arrays = copy.deepcopy(arrays)
        invalid_arrays["volumes"][1] = 0.0
        with self.assertRaisesRegex(ValueError, "finite and positive"):
            replay_analysis.validate_replay_frame_health(
                Path("frame.vtp"), points, invalid_arrays, lower, upper
            )

    def test_replay_analysis_helpers_are_finite_and_nan_safe(self) -> None:
        mesh = self.root / "small_mesh.txt"
        mesh.write_text(
            "#! elementShape quadrilateral\n"
            "4 1\n"
            "0 0\n"
            "2 0\n"
            "0 1\n"
            "2 1\n"
            "0 1 3 2\n",
            encoding="utf-8",
        )
        lower, upper = replay_analysis.mesh_bounds(mesh)
        self.assertEqual(lower.tolist(), [0.0, 0.0])
        self.assertEqual(upper.tolist(), [2.0, 1.0])
        self.assertAlmostEqual(
            replay_analysis.integrate(
                np.asarray([0.0, 1.0, 2.0]), np.asarray([0.0, 1.0, 0.0])
            ),
            1.0,
        )
        self.assertAlmostEqual(
            replay_analysis.robust_symmetric_limit(
                np.asarray([np.nan, -2.0, 2.0])
            ),
            2.0,
        )

    def test_phase_alignment_gate_rejects_missing_fields(self) -> None:
        valid = self.valid_phase_alignment_metadata()
        accepted = phase_controls.validate_phase_alignment_qa(valid)
        self.assertEqual(accepted["phase_defined_particle_count"], 10)
        missing_paths = (
            ("qa",),
            ("qa", "phase_alignment"),
            ("qa", "phase_alignment", "liquid"),
            (
                "qa",
                "phase_alignment",
                "liquid",
                "phase_defined_particle_count",
            ),
            (
                "qa",
                "phase_alignment",
                "liquid",
                "max_abs_source_to_surface_phase_difference_deg",
            ),
        )
        for keys in missing_paths:
            changed = copy.deepcopy(valid)
            target = changed
            for key in keys[:-1]:
                target = target[key]
            del target[keys[-1]]
            with self.subTest(keys=keys), self.assertRaisesRegex(
                ValueError, "missing"
            ):
                phase_controls.validate_phase_alignment_qa(changed)

    def test_phase_alignment_gate_requires_defined_changed_source_phase(self) -> None:
        cases = (
            ("phase_defined_particle_count", 0, "defined phase"),
            ("changed_particle_count", 0, "no changed particles"),
            (
                "max_abs_source_to_surface_phase_difference_deg",
                0.0,
                "not positive",
            ),
        )
        for key, value, message in cases:
            metadata = self.valid_phase_alignment_metadata()
            metadata["qa"]["phase_alignment"]["liquid"][key] = value
            with self.subTest(key=key), self.assertRaisesRegex(ValueError, message):
                phase_controls.validate_phase_alignment_qa(metadata)

    def test_phase_alignment_gate_rejects_nonzero_residual_above_tolerance(self) -> None:
        metadata = self.valid_phase_alignment_metadata()
        metadata["qa"]["phase_alignment"]["liquid"][
            "max_abs_residual_to_surface_phase_difference_deg"
        ] = 10.0 * phase_controls.PHASE_ALIGNMENT_RESIDUAL_TOLERANCE_DEG
        with self.assertRaisesRegex(ValueError, "residual exceeds"):
            phase_controls.validate_phase_alignment_qa(metadata)

    def test_phase_alignment_gate_validates_gas_qa_too(self) -> None:
        for field, value, message in (
            ("changed_particle_count", 0, "Gas.*no changed particles"),
            (
                "max_abs_source_to_surface_phase_difference_deg",
                math.nan,
                "Gas.*invalid",
            ),
        ):
            metadata = self.valid_phase_alignment_metadata()
            metadata["qa"]["phase_alignment"]["gas"][field] = value
            with self.subTest(field=field), self.assertRaisesRegex(
                ValueError, message
            ):
                phase_controls.validate_phase_alignment_qa(metadata)


if __name__ == "__main__":
    unittest.main()
