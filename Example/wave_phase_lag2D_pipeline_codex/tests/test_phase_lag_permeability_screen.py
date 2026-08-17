from __future__ import annotations

import copy
import csv
import json
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import summarize_phase_lag_permeability_screen as summary


class PhaseLagPermeabilityScreenTests(unittest.TestCase):
    def test_default_labels_use_the_four_current_runner_cases(self) -> None:
        self.assertEqual(
            summary.DEFAULT_LABELS,
            (
                "k1e-13_sw094_nosmooth_fresh",
                "k3e-13_sw094_nosmooth",
                "k1e-12_sw094_nosmooth",
                "k3e-12_sw094_nosmooth",
            ),
        )

    @staticmethod
    def _parameters(permeability: float) -> dict:
        return {
            "intrinsic_permeability_m2": permeability,
            "liquid_saturation": 0.94,
            "gas_saturation": 0.06,
            "pressure_smoothing": False,
            "pressure_smoothing_in_loop": False,
            "pipeline_fixed": True,
        }

    @staticmethod
    def _integral(lagged: float, erased: float, units: str) -> dict:
        return {
            "units": units,
            "lagged": lagged,
            "phase_erased": erased,
            "lagged_minus_phase_erased": lagged - erased,
            "lagged_minus_phase_erased_fraction": (
                None if erased <= 0.0 else lagged / erased - 1.0
            ),
        }

    @staticmethod
    def _artifact(path: Path, payload: bytes | None = None) -> dict:
        path.parent.mkdir(parents=True, exist_ok=True)
        if payload is not None:
            path.write_bytes(payload)
        return {"path": str(path.resolve()), "sha256": summary.sha256(path)}

    @staticmethod
    def _sized_artifact(path: Path, payload: bytes | None = None) -> dict:
        path.parent.mkdir(parents=True, exist_ok=True)
        if payload is not None:
            path.write_bytes(payload)
        return summary.phase_plot.artifact_record(path)

    def _write_case(
        self,
        root: Path,
        label: str,
        permeability: float,
        *,
        effect_fraction: float = 0.1,
        net_effect_fraction: float | None = None,
        crown_phase: float | None = None,
        crown_phase_undefined: bool = False,
        zero_denominator_metric: str | None = None,
    ) -> tuple[Path, Path]:
        base = root / "analysis/phase_lag_exploratory" / label
        phase_dir = base / "phase_v3"
        driver_dir = base / "driver_screen"
        phase_dir.mkdir(parents=True)
        driver_dir.mkdir(parents=True)
        parameters = self._parameters(permeability)
        case = label.upper()
        if net_effect_fraction is None:
            net_effect_fraction = effect_fraction
        if crown_phase_undefined:
            crown_phase = None
        elif crown_phase is None:
            crown_phase = -20.0 - permeability * 1.0e13
        figure_paths = tuple(
            phase_dir / f"{case}_{stem}{suffix}"
            for stem in (
                "liquid_pressure_2d_phase_lag",
                "liquid_pressure_time_depth",
            )
            for suffix in (".png", ".pdf")
        )
        for index, figure_path in enumerate(figure_paths):
            figure_path.write_bytes(f"synthetic figure {index}\n".encode())
        phase = {
            "schema": summary.PHASE_SCHEMA,
            "case": case,
            "parameters": copy.deepcopy(parameters),
            "require_unsmoothed": True,
            "frame_count": 60,
            "particle_count": 100,
            "wave_period_s": 1.3,
            "fit_window_s": [1.3, 3.9],
            "probe_results": {
                "surface": {"fundamental_amplitude_pa": 200.0},
                "crown": {
                    "fundamental_amplitude_pa": 50.0 + permeability * 1.0e13,
                    "harmonic_r_squared": 0.95,
                    "phase_difference_from_same_column_surface_deg": crown_phase,
                },
            },
            "phase_present": crown_phase is not None,
            "artifacts": {
                "schema": summary.phase_plot.ARTIFACT_SCHEMA,
                "analyzer": summary.phase_plot.artifact_record(
                    summary.phase_plot.ANALYZER_PATH
                ),
                "plot_script": summary.phase_plot.artifact_record(
                    summary.phase_plot.PLOT_SCRIPT_PATH
                ),
                "figures": [
                    summary.phase_plot.artifact_record(path)
                    for path in figure_paths
                ],
            },
        }
        particles_path = root / "particles.txt"
        if not particles_path.exists():
            particles_path.write_bytes(b"synthetic particles\n")
        eq_result = root / "results/phase_lag_exploratory" / label / "EQ"
        hd_result = root / "results/phase_lag_exploratory" / label / "HD"
        checkpoint = self._sized_artifact(
            eq_result / "particle5000.vtp", b"equilibrium checkpoint\n"
        )
        pipeline_history = self._sized_artifact(
            hd_result / "pipeline-history0000.csv", b"pipeline history\n"
        )
        phase_frames = [
            self._sized_artifact(
                hd_result / f"particle{650 * index:05d}.vtp",
                f"phase frame {index}\n".encode(),
            )
            for index in range(1, summary.EXPECTED_FRAME_COUNT + 1)
        ]

        base_rates = {
            "signed_force_impulse": 4.0,
            "net_uplift_force_impulse": 4.0,
            "local_positive_force_activity_impulse": 5.0,
            "IF_ge_1_area_time": 1.0,
            "shared_HD_joint_area_time": 1.0,
        }
        effect_fractions = {
            metric: (
                net_effect_fraction
                if metric in {"signed_force_impulse", "net_uplift_force_impulse"}
                else effect_fraction
            )
            for metric in base_rates
        }
        if zero_denominator_metric is not None:
            if zero_denominator_metric not in base_rates:
                raise AssertionError("unknown synthetic zero-denominator metric")
            base_rates[zero_denominator_metric] = 0.0

        raw_names = {
            output: raw for output, raw, _units in summary.HISTORY_METRICS
        }
        fieldnames = list(summary.HISTORY_REQUIRED_FIELDS)
        rows = []
        for index in range(1, 61):
            row = {"time_s": 0.065 * index}
            for output_name, raw_name, _units in summary.HISTORY_METRICS:
                erased = base_rates[output_name]
                if output_name == zero_denominator_metric:
                    lagged = 0.10
                else:
                    lagged = erased * (1.0 + effect_fractions[output_name])
                row[f"lagged_{raw_name}"] = lagged
                row[f"phase_erased_{raw_name}"] = erased
            row["lagged_maximum_upward_if"] = 2.0
            row["phase_erased_maximum_upward_if"] = 1.8
            row["delta_positive_force_n_per_m"] = (
                row["lagged_positive_force_n_per_m"]
                - row["phase_erased_positive_force_n_per_m"]
            )
            row["delta_signed_force_n_per_m"] = (
                row["lagged_signed_force_n_per_m"]
                - row["phase_erased_signed_force_n_per_m"]
            )
            row["delta_net_uplift_force_n_per_m"] = (
                row["lagged_net_uplift_force_n_per_m"]
                - row["phase_erased_net_uplift_force_n_per_m"]
            )
            rows.append(row)

        history_path = driver_dir / "pressure_only_gradient_force_history.csv"
        with history_path.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(rows)

        windows = []
        for window_index, (start, end) in enumerate(summary.EXPECTED_WINDOWS_S):
            duration = end - start
            window = {
                "start_s": start,
                "end_s": end,
                "saved_frames": 41 if window_index == 2 else 21,
            }
            for output_name, _raw_name, units in summary.HISTORY_METRICS:
                erased = base_rates[output_name] * duration
                lagged_rate = rows[0][f"lagged_{raw_names[output_name]}"]
                window[output_name] = self._integral(
                    lagged_rate * duration, erased, units
                )
            windows.append(window)

        config_dir = root / "configs/phase_lag_exploratory" / label
        configs = {}
        for role, filename in (
            ("HD", "03_HD.json"),
            ("RL", "04_RL.json"),
            ("RE", "05_RE.json"),
        ):
            configs[role] = self._artifact(
                config_dir / filename,
                json.dumps({"role": role, "parameters": parameters}).encode(),
            )
        eq_config = self._sized_artifact(
            config_dir / "01_EQ.json",
            json.dumps({"role": "EQ", "parameters": parameters}).encode(),
        )
        hd_config_sized = summary.phase_plot.artifact_record(
            Path(configs["HD"]["path"])
        )
        phase["artifacts"]["inputs"] = {
            "schema": summary.phase_plot.INPUT_SCHEMA,
            "analysis_config": hd_config_sized,
            "checkpoint_config": eq_config,
            "checkpoint": checkpoint,
            "particles": summary.phase_plot.artifact_record(particles_path),
            "pipeline_history": pipeline_history,
            "result_directory": str(hd_result.resolve()),
            "result_particle_vtp": phase_frames,
        }
        phase["config"] = hd_config_sized
        phase["checkpoint"] = checkpoint

        database_root = root / "pressure_databases/phase_lag_exploratory" / label
        database_paths = {
            "lagged": {
                "points": database_root / "lagged/pressure_points.txt",
                "values": database_root / "lagged/pressure_values.bin",
            },
            "phase_erased": {
                "points": database_root / "phase_erased/phase_erased_points.txt",
                "values": database_root / "phase_erased/phase_erased_values.bin",
            },
        }
        database_records = {}
        for role, paths in database_paths.items():
            points = self._artifact(paths["points"], f"{role} points\n".encode())
            values = self._artifact(paths["values"], f"{role} values\n".encode())
            database_records[role] = {
                "points": points["path"],
                "points_sha256": points["sha256"],
                "values": values["path"],
                "values_sha256": values["sha256"],
                "step_interval": 130,
                "max_step": 39000,
                "source_dt_s": 1.0e-4,
            }

        surface_reference = self._artifact(
            root / "top_surface_traction_particle_id.txt", b"1\n2\n"
        )
        metadata = {
            "schema": "phase-erased-pressure-control-v1",
            "source": {
                "points": database_records["lagged"]["points"],
                "points_sha256": database_records["lagged"]["points_sha256"],
                "values": database_records["lagged"]["values"],
                "values_sha256": database_records["lagged"]["values_sha256"],
            },
            "output": {
                "points": database_records["phase_erased"]["points"],
                "points_sha256": database_records["phase_erased"]["points_sha256"],
                "values": database_records["phase_erased"]["values"],
                "values_sha256": database_records["phase_erased"]["values_sha256"],
            },
            "transform": {"surface_reference": {"source": surface_reference}},
        }
        transform_path = database_root / "phase_erased/phase_erased_metadata.json"
        transform_path.write_text(json.dumps(metadata), encoding="utf-8")

        runner_path = base / summary.RUNNER_AUDIT_NAME
        runner = {
            "schema": summary.RUNNER_SCHEMA,
            "stage": summary.RUNNER_STAGE,
            "label": label,
            "parameters": {
                "intrinsic_permeability_m2": permeability,
                "liquid_saturation": 0.94,
                "pressure_smoothing": False,
            },
            "stages": {
                "eq": {
                    "config": eq_config,
                    "audit": {
                        "final_vtp": checkpoint,
                        "stability_qa": {
                            "normal_stress": {
                                "inputs": {
                                    "particles": {
                                        "path": "particles.txt",
                                        "size_bytes": particles_path.stat().st_size,
                                        "sha256": summary.sha256(particles_path),
                                    }
                                }
                            }
                        },
                    },
                },
                "hd": {
                    "config": hd_config_sized,
                    "audit": {
                        "result_directory": str(hd_result.resolve()),
                        "particle_vtp": phase_frames,
                        "pipeline_history": pipeline_history,
                    },
                },
            },
        }
        runner_path.write_text(json.dumps(runner), encoding="utf-8")
        runner_binding = {
            "schema": summary.PHASE_RUNNER_BINDING_SCHEMA,
            "eq": {
                "checkpoint_config": eq_config,
                "checkpoint": checkpoint,
                "particles": summary.phase_plot.artifact_record(particles_path),
            },
            "hd": {
                "analysis_config": hd_config_sized,
                "result_directory": str(hd_result.resolve()),
                "result_particle_vtp": phase_frames,
                "pipeline_history": pipeline_history,
            },
        }

        phase_checks = {
            "crown_abs_same_column_phase_lag_ge_18deg": (
                crown_phase is not None and abs(crown_phase) >= 18.0
            ),
            "crown_amplitude_ratio_ge_0p05": True,
            "crown_harmonic_r_squared_ge_0p8": True,
        }
        qualification = {
            "passed": all(phase_checks.values()),
            "checks": phase_checks,
            "observed": {
                "crown_same_column_phase_lag_deg": crown_phase,
                "crown_amplitude_ratio": (
                    phase["probe_results"]["crown"]["fundamental_amplitude_pa"]
                    / phase["probe_results"]["surface"]["fundamental_amplitude_pa"]
                ),
                "crown_harmonic_r_squared": 0.95,
            },
        }
        net_fraction = windows[-1]["net_uplift_force_impulse"][
            "lagged_minus_phase_erased_fraction"
        ]
        joint_fraction = windows[-1]["shared_HD_joint_area_time"][
            "lagged_minus_phase_erased_fraction"
        ]
        gate_checks = {
            "crown_phase_lag_is_resolved": qualification["passed"],
            "both_post_ramp_cycles_have_positive_signed_force_delta": (
                net_effect_fraction > 0.0
            ),
            "combined_signed_force_impulse_delta_positive": (
                net_effect_fraction > 0.0
            ),
            "both_post_ramp_cycles_have_positive_net_uplift_delta": (
                net_effect_fraction > 0.0
            ),
            "combined_net_uplift_impulse_increase_ge_5pct": (
                net_fraction is not None and net_fraction >= 0.05
            ),
            "combined_IF_ge_1_area_time_direction_positive": True,
            "combined_shared_HD_joint_area_time_nonzero": True,
            "combined_shared_HD_joint_area_time_increase_ge_5pct": (
                joint_fraction is not None and joint_fraction >= 0.05
            ),
            "resolved_joint_advantage_for_two_consecutive_saved_frames": True,
        }
        gate_passed = all(gate_checks.values())

        driver_png = self._artifact(
            driver_dir / "pressure_only_gradient_force.png", b"driver png\n"
        )
        driver_pdf = self._artifact(
            driver_dir / "pressure_only_gradient_force.pdf", b"driver pdf\n"
        )
        driver = {
            "schema": summary.DRIVER_SCHEMA,
            "classification": "exploratory-pressure-only-fixed-state-diagnostic",
            "label": label,
            "parameters": copy.deepcopy(parameters),
            "screen_window": {"start_s": 1.3, "end_s": 3.9, "period_s": 1.3},
            "frame_count": 60,
            "particle_count": 100,
            "support_particle_count": 20,
            "runner_audit": {
                "path": str(runner_path.resolve()),
                "sha256": summary.sha256(runner_path),
            },
            "analyzer": self._artifact(summary.DRIVER_ANALYZER_PATH),
            "implementation_dependencies": {
                name: self._artifact(path)
                for name, path in summary.IMPLEMENTATION_DEPENDENCY_PATHS.items()
            },
            "configs": configs,
            "phase_transform": {
                "path": str(transform_path.resolve()),
                "sha256": summary.sha256(transform_path),
                "surface_reference_source": {
                    **surface_reference,
                    "size_bytes": Path(surface_reference["path"]).stat().st_size,
                },
                "metadata": metadata,
            },
            "phase_field": {
                "path": str((phase_dir / f"{case}{summary.PHASE_AUDIT_SUFFIX}").resolve()),
                "sha256": "pending",
                "qualification": qualification,
                "runner_binding": runner_binding,
            },
            "pressure_databases": database_records,
            "integration_windows": windows,
            "diagnostic_gate": {
                "passed": gate_passed,
                "checks": gate_checks,
                "classification": (
                    "pressure-only-screen-passed"
                    if gate_passed
                    else "pressure-only-screen-not-passed"
                ),
                "minimum_resolved_joint_advantage_area_m2": 0.01,
                "area_comparison_absolute_tolerance_m2": 1.0e-12,
                "phase_qualification": qualification,
            },
            "figure": {
                "display_smoothing_sigma_pixels": 100.0,
                "display_smoothing_used_in_metrics": False,
                "png": driver_png,
                "pdf": driver_pdf,
            },
            "history_csv": str(history_path.resolve()),
            "history_csv_sha256": summary.sha256(history_path),
        }
        phase_path = phase_dir / f"{case}{summary.PHASE_AUDIT_SUFFIX}"
        driver_path = driver_dir / summary.DRIVER_AUDIT_NAME
        phase_path.write_text(json.dumps(phase), encoding="utf-8")
        driver["phase_field"]["sha256"] = summary.sha256(phase_path)
        driver_path.write_text(json.dumps(driver), encoding="utf-8")
        return phase_path, driver_path

    def _four_cases(self, root: Path, *, failed_index: int | None = None) -> list[str]:
        specifications = list(
            zip(summary.DEFAULT_LABELS, summary.EXPECTED_PERMEABILITIES_M2)
        )
        labels = []
        for index, (label, permeability) in enumerate(specifications):
            self._write_case(
                root,
                label,
                permeability,
                effect_fraction=0.10 + 0.01 * index,
                net_effect_fraction=(
                    0.01 if index == failed_index else 0.10 + 0.01 * index
                ),
            )
            labels.append(label)
        return labels

    def test_cases_are_sorted_by_permeability_and_raw_effect_is_recomputed(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            labels = self._four_cases(root)
            cases = summary.collect_cases(root, list(reversed(labels)))
            self.assertEqual([case["label"] for case in cases], labels)
            self.assertEqual(
                [case["intrinsic_permeability_m2"] for case in cases],
                [1.0e-13, 3.0e-13, 1.0e-12, 3.0e-12],
            )
            self.assertAlmostEqual(cases[2]["combined_net_uplift_effect_percent"], 12.0)
            self.assertEqual(
                cases[0]["evidence_role"],
                "primary/current-runner",
            )
            self.assertTrue(cases[0]["primary_screen_point"])
            self.assertTrue(cases[0]["current_runner_point"])
            self.assertTrue(cases[1]["primary_screen_point"])

    def test_main_matrix_rejects_repeat_substitution_duplicate_or_extra_labels(self) -> None:
        registered = list(summary.DEFAULT_LABELS)
        cases = (
            (
                [
                    "k1e-13_sw094_nosmooth_repeat",
                    *registered[1:],
                ],
                "exactly match the registered DEFAULT_LABELS",
            ),
            (
                [registered[0], registered[0], *registered[2:]],
                "labels must be unique",
            ),
            (
                [*registered, "unregistered-extra"],
                "requires exactly 4 labels",
            ),
            (registered[:-1], "requires exactly 4 labels"),
        )
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            for labels, message in cases:
                with self.subTest(labels=labels), self.assertRaisesRegex(
                    summary.SummaryError, message
                ):
                    summary.collect_cases(root, labels)

    def test_phase_artifacts_fail_closed_on_missing_or_drift(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            phase_path, _driver_path = self._write_case(root, "k1", 1.0e-13)
            self.assertEqual(summary.load_case(root, "k1")["label"], "k1")

        for mutation in ("missing_artifacts", "missing_figure", "script_drift"):
            with (
                self.subTest(mutation=mutation),
                tempfile.TemporaryDirectory() as temporary,
            ):
                root = Path(temporary)
                phase_path, _driver_path = self._write_case(root, "k1", 1.0e-13)
                phase = json.loads(phase_path.read_text(encoding="utf-8"))
                if mutation == "missing_artifacts":
                    phase.pop("artifacts")
                    phase_path.write_text(json.dumps(phase), encoding="utf-8")
                elif mutation == "missing_figure":
                    Path(phase["artifacts"]["figures"][0]["path"]).unlink()
                else:
                    phase["artifacts"]["plot_script"]["sha256"] = "0" * 64
                    phase_path.write_text(json.dumps(phase), encoding="utf-8")
                with self.assertRaisesRegex(
                    summary.SummaryError, "artifacts are incomplete or stale"
                ) as caught:
                    summary.load_case(root, "k1")
                self.assertIsInstance(
                    caught.exception.__cause__, summary.phase_plot.PhaseAuditError
                )

    def test_phase_inputs_and_driver_runner_binding_fail_closed(self) -> None:
        for mutation, message in (
            ("phase_history", "artifacts are incomplete or stale"),
            ("driver_frame", "binding differs from runner EQ/HD artifacts"),
            ("runner_frame", "runner HD VTP 9 artifact hash differs"),
        ):
            with (
                self.subTest(mutation=mutation),
                tempfile.TemporaryDirectory() as temporary,
            ):
                root = Path(temporary)
                phase_path, driver_path = self._write_case(
                    root, "k1", 1.0e-13
                )
                phase = json.loads(phase_path.read_text(encoding="utf-8"))
                driver = json.loads(driver_path.read_text(encoding="utf-8"))
                runner_path = (
                    root
                    / "analysis/phase_lag_exploratory/k1"
                    / summary.RUNNER_AUDIT_NAME
                )
                runner = json.loads(runner_path.read_text(encoding="utf-8"))
                if mutation == "phase_history":
                    phase["artifacts"]["inputs"]["pipeline_history"][
                        "sha256"
                    ] = "0" * 64
                    phase_path.write_text(json.dumps(phase), encoding="utf-8")
                    driver["phase_field"]["sha256"] = summary.sha256(phase_path)
                elif mutation == "driver_frame":
                    driver["phase_field"]["runner_binding"]["hd"][
                        "result_particle_vtp"
                    ][9]["sha256"] = "0" * 64
                else:
                    runner["stages"]["hd"]["audit"]["particle_vtp"][9][
                        "sha256"
                    ] = "0" * 64
                    runner_path.write_text(json.dumps(runner), encoding="utf-8")
                    driver["runner_audit"]["sha256"] = summary.sha256(runner_path)
                driver_path.write_text(json.dumps(driver), encoding="utf-8")
                with self.assertRaisesRegex(summary.SummaryError, message):
                    summary.load_case(root, "k1")

    def test_runner_binding_is_revalidated_and_exposed_as_a_source(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            _phase_path, _driver_path = self._write_case(root, "k1", 1.0e-13)
            runner_path = (
                root
                / "analysis/phase_lag_exploratory/k1"
                / summary.RUNNER_AUDIT_NAME
            )
            case = summary.load_case(root, "k1")
            self.assertTrue(case["runner_provenance_verified"])
            self.assertEqual(case["runner_audit_schema"], summary.RUNNER_SCHEMA)
            self.assertEqual(case["runner_audit_stage"], summary.RUNNER_STAGE)
            self.assertEqual(
                case["sources"]["runner_audit"],
                {
                    "path": str(runner_path.resolve()),
                    "sha256": summary.sha256(runner_path),
                    "schema": summary.RUNNER_SCHEMA,
                    "stage": summary.RUNNER_STAGE,
                },
            )

    def test_runner_binding_fails_closed_on_missing_or_drift(self) -> None:
        mutations = (
            ("missing_record", "runner_audit must be a JSON object"),
            ("missing_file", "Missing .* current runner audit"),
            ("missing_hash", "runner_audit.sha256"),
            ("hash_drift", "runner audit hash differs"),
            ("path_drift", "runner_audit.path does not bind"),
            ("schema_drift", "runner audit schema must be"),
            ("stage_drift", "runner audit stage must be"),
            ("label_drift", "runner audit label identity differs"),
            ("permeability_drift", "runner/driver parameter.*differs"),
            ("saturation_drift", "runner/driver parameter.*differs"),
            ("smoothing_drift", "runner/driver parameter.*differs"),
        )
        for mutation, message in mutations:
            with (
                self.subTest(mutation=mutation),
                tempfile.TemporaryDirectory() as temporary,
            ):
                root = Path(temporary)
                _phase_path, driver_path = self._write_case(root, "k1", 1.0e-13)
                runner_path = (
                    root
                    / "analysis/phase_lag_exploratory/k1"
                    / summary.RUNNER_AUDIT_NAME
                )
                driver = json.loads(driver_path.read_text(encoding="utf-8"))
                runner = json.loads(runner_path.read_text(encoding="utf-8"))
                rewrite_runner = False
                if mutation == "missing_record":
                    driver.pop("runner_audit")
                elif mutation == "missing_file":
                    runner_path.unlink()
                elif mutation == "missing_hash":
                    driver["runner_audit"].pop("sha256")
                elif mutation == "hash_drift":
                    driver["runner_audit"]["sha256"] = "0" * 64
                elif mutation == "path_drift":
                    other_path = runner_path.with_name("other_runner_audit.json")
                    other_path.write_text(json.dumps(runner), encoding="utf-8")
                    driver["runner_audit"] = {
                        "path": str(other_path.resolve()),
                        "sha256": summary.sha256(other_path),
                    }
                elif mutation == "schema_drift":
                    runner["schema"] = "pipeline-phase-lag-exploration-runner-v0"
                    rewrite_runner = True
                elif mutation == "stage_drift":
                    runner["stage"] = "eq_complete"
                    rewrite_runner = True
                elif mutation == "label_drift":
                    runner["label"] = "another-label"
                    rewrite_runner = True
                elif mutation == "permeability_drift":
                    runner["parameters"]["intrinsic_permeability_m2"] = 2.0e-13
                    rewrite_runner = True
                elif mutation == "saturation_drift":
                    runner["parameters"]["liquid_saturation"] = 0.93
                    rewrite_runner = True
                else:
                    runner["parameters"]["pressure_smoothing"] = True
                    rewrite_runner = True
                if rewrite_runner:
                    runner_path.write_text(json.dumps(runner), encoding="utf-8")
                    driver["runner_audit"]["sha256"] = summary.sha256(runner_path)
                driver_path.write_text(json.dumps(driver), encoding="utf-8")
                with self.assertRaisesRegex(summary.SummaryError, message):
                    summary.load_case(root, "k1")

    def test_phase_driver_parameter_drift_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            labels = self._four_cases(root)
            driver_path = (
                root
                / "analysis/phase_lag_exploratory"
                / labels[0]
                / "driver_screen"
                / summary.DRIVER_AUDIT_NAME
            )
            driver = json.loads(driver_path.read_text(encoding="utf-8"))
            driver["parameters"]["intrinsic_permeability_m2"] = 2.0e-13
            driver_path.write_text(json.dumps(driver), encoding="utf-8")
            with self.assertRaisesRegex(summary.SummaryError, "parameter.*differs"):
                summary.collect_cases(root, labels)

    def test_phase_field_is_strongly_bound_to_driver_path_and_hash(self) -> None:
        for mutation, message in (
            ("path", "phase_field.path does not bind"),
            ("hash", "phase_field artifact hash differs"),
        ):
            with (
                self.subTest(mutation=mutation),
                tempfile.TemporaryDirectory() as temporary,
            ):
                root = Path(temporary)
                phase_path, driver_path = self._write_case(root, "k1", 1.0e-13)
                driver = json.loads(driver_path.read_text(encoding="utf-8"))
                if mutation == "path":
                    substitute = phase_path.with_name("substitute.audit.json")
                    substitute.write_bytes(phase_path.read_bytes())
                    driver["phase_field"].update(
                        {
                            "path": str(substitute.resolve()),
                            "sha256": summary.sha256(substitute),
                        }
                    )
                else:
                    driver["phase_field"]["sha256"] = "0" * 64
                driver_path.write_text(json.dumps(driver), encoding="utf-8")
                with self.assertRaisesRegex(summary.SummaryError, message):
                    summary.load_case(root, "k1")

    def test_driver_file_provenance_fails_closed_for_every_required_family(self) -> None:
        mutations = (
            ("analyzer", "driver analyzer artifact hash differs"),
            ("dependency_missing", "must contain exactly"),
            ("dependency", "implementation dependency runner artifact hash differs"),
            ("config", "HD config artifact hash differs"),
            ("database", "lagged database values artifact hash differs"),
            ("history", "driver history artifact hash differs"),
            ("figure", "driver figure png artifact hash differs"),
            ("transform", "phase transform metadata artifact hash differs"),
        )
        for mutation, message in mutations:
            with (
                self.subTest(mutation=mutation),
                tempfile.TemporaryDirectory() as temporary,
            ):
                root = Path(temporary)
                _phase_path, driver_path = self._write_case(root, "k1", 1.0e-13)
                driver = json.loads(driver_path.read_text(encoding="utf-8"))
                if mutation == "analyzer":
                    driver["analyzer"]["sha256"] = "0" * 64
                elif mutation == "dependency_missing":
                    driver["implementation_dependencies"].pop("runner")
                elif mutation == "dependency":
                    driver["implementation_dependencies"]["runner"]["sha256"] = (
                        "0" * 64
                    )
                elif mutation == "config":
                    driver["configs"]["HD"]["sha256"] = "0" * 64
                elif mutation == "database":
                    driver["pressure_databases"]["lagged"]["values_sha256"] = (
                        "0" * 64
                    )
                elif mutation == "history":
                    driver["history_csv_sha256"] = "0" * 64
                elif mutation == "figure":
                    driver["figure"]["png"]["sha256"] = "0" * 64
                else:
                    driver["phase_transform"]["sha256"] = "0" * 64
                driver_path.write_text(json.dumps(driver), encoding="utf-8")
                with self.assertRaisesRegex(summary.SummaryError, message):
                    summary.load_case(root, "k1")

    def test_complete_gate_and_raw_window_recomputation_fail_closed(self) -> None:
        mutations = (
            ("missing_check", "complete preregistered set"),
            ("extra_check", "complete preregistered set"),
            ("check_drift", "differs from raw-history recomputation"),
            ("cycle_integral_drift", "raw-history recomputation"),
            ("history_drift", "raw-history recomputation"),
        )
        for mutation, message in mutations:
            with (
                self.subTest(mutation=mutation),
                tempfile.TemporaryDirectory() as temporary,
            ):
                root = Path(temporary)
                _phase_path, driver_path = self._write_case(root, "k1", 1.0e-13)
                driver = json.loads(driver_path.read_text(encoding="utf-8"))
                checks = driver["diagnostic_gate"]["checks"]
                if mutation == "missing_check":
                    checks.pop("combined_signed_force_impulse_delta_positive")
                elif mutation == "extra_check":
                    checks["unregistered_override"] = True
                elif mutation == "check_drift":
                    checks[
                        "both_post_ramp_cycles_have_positive_signed_force_delta"
                    ] = False
                    driver["diagnostic_gate"]["passed"] = False
                    driver["diagnostic_gate"]["classification"] = (
                        "pressure-only-screen-not-passed"
                    )
                elif mutation == "cycle_integral_drift":
                    driver["integration_windows"][0]["signed_force_impulse"][
                        "lagged"
                    ] += 1.0
                else:
                    history_path = Path(driver["history_csv"])
                    with history_path.open(encoding="utf-8") as stream:
                        rows = list(csv.DictReader(stream))
                    rows[20]["lagged_signed_force_n_per_m"] = "99.0"
                    rows[20]["lagged_net_uplift_force_n_per_m"] = "99.0"
                    rows[20]["lagged_positive_force_n_per_m"] = "100.0"
                    rows[20]["delta_signed_force_n_per_m"] = str(
                        99.0 - float(rows[20]["phase_erased_signed_force_n_per_m"])
                    )
                    rows[20]["delta_net_uplift_force_n_per_m"] = str(
                        99.0
                        - float(rows[20]["phase_erased_net_uplift_force_n_per_m"])
                    )
                    rows[20]["delta_positive_force_n_per_m"] = str(
                        100.0
                        - float(rows[20]["phase_erased_positive_force_n_per_m"])
                    )
                    with history_path.open("w", newline="", encoding="utf-8") as stream:
                        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
                        writer.writeheader()
                        writer.writerows(rows)
                    driver["history_csv_sha256"] = summary.sha256(history_path)
                driver_path.write_text(json.dumps(driver), encoding="utf-8")
                with self.assertRaisesRegex(summary.SummaryError, message):
                    summary.load_case(root, "k1")

    def test_undefined_crown_phase_is_a_legal_qualification_fail(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            self._write_case(
                root, "k1", 1.0e-13, crown_phase_undefined=True
            )
            case = summary.load_case(root, "k1")
            self.assertIsNone(case["crown_phase_difference_same_column_deg"])
            self.assertIsNone(case["crown_phase_lag_magnitude_deg"])
            self.assertFalse(case["phase_qualification"]["passed"])
            self.assertFalse(case["diagnostic_gate_passed"])
            self.assertFalse(
                case["diagnostic_gate_checks"]["crown_phase_lag_is_resolved"]
            )

    def test_zero_denominator_is_a_legal_na_and_cannot_pass_gate(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            self._write_case(
                root,
                "k1",
                1.0e-13,
                zero_denominator_metric="shared_HD_joint_area_time",
            )
            case = summary.load_case(root, "k1")
            self.assertIsNone(
                case["combined_shared_HD_joint_area_time_effect_percent"]
            )
            self.assertIsNone(
                case["combined_raw_integrals"]["shared_HD_joint_area_time"][
                    "lagged_minus_phase_erased_fraction"
                ]
            )
            self.assertFalse(
                case["diagnostic_gate_checks"][
                    "combined_shared_HD_joint_area_time_increase_ge_5pct"
                ]
            )
            self.assertFalse(case["diagnostic_gate_passed"])

    def test_missing_or_legacy_schema_is_never_silently_included(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            labels = self._four_cases(root)
            driver_path = (
                root
                / "analysis/phase_lag_exploratory"
                / labels[1]
                / "driver_screen"
                / summary.DRIVER_AUDIT_NAME
            )
            driver = json.loads(driver_path.read_text(encoding="utf-8"))
            driver.pop("schema")
            driver_path.write_text(json.dumps(driver), encoding="utf-8")
            with self.assertRaisesRegex(summary.SummaryError, "schema must be"):
                summary.collect_cases(root, labels)

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            labels = self._four_cases(root)
            phase_v3 = (
                root
                / "analysis/phase_lag_exploratory"
                / labels[2]
                / "phase_v3"
            )
            phase_path = next(phase_v3.glob(f"*{summary.PHASE_AUDIT_SUFFIX}"))
            legacy_path = phase_v3.parent / phase_path.name
            legacy_path.write_text(phase_path.read_text(encoding="utf-8"), encoding="utf-8")
            phase_path.unlink()
            with self.assertRaisesRegex(summary.SummaryError, "Legacy or substitute"):
                summary.collect_cases(root, labels)

    def test_failed_gate_is_reported_instead_of_filtered_out(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            labels = self._four_cases(root, failed_index=1)
            output = root / "summary"
            audit = summary.write_summary(root, labels, output, dpi=72)
            self.assertEqual(audit["case_count"], 4)
            self.assertEqual(audit["passed_case_count"], 3)
            self.assertTrue(audit["evidence_policy"]["all_four_points_primary"])
            self.assertEqual(
                audit["evidence_policy"]["primary_role"],
                "primary/current-runner",
            )
            self.assertFalse(
                audit["evidence_policy"][
                    "external_repeatability_references_included"
                ]
            )
            failed = [
                case for case in audit["cases"] if not case["diagnostic_gate_passed"]
            ]
            self.assertEqual(
                [case["label"] for case in failed],
                [summary.DEFAULT_LABELS[1]],
            )
            self.assertEqual(failed[0]["diagnostic_gate_status"], "FAIL")
            self.assertEqual(failed[0]["net_uplift_gate_status"], "FAIL")
            self.assertGreater(
                failed[0]["combined_local_positive_activity_effect_percent"], 0.0
            )
            self.assertFalse(
                failed[0]["local_positive_activity_can_override_net_gate"]
            )
            with Path(audit["outputs"]["csv"]["path"]).open(
                newline="", encoding="utf-8"
            ) as stream:
                rows = list(csv.DictReader(stream))
            self.assertEqual(len(rows), 4)
            self.assertEqual(rows[1]["diagnostic_gate_status"], "FAIL")
            self.assertEqual(
                {row["evidence_role"] for row in rows},
                {"primary/current-runner"},
            )
            self.assertEqual({row["primary_screen_point"] for row in rows}, {"True"})
            self.assertEqual({row["current_runner_point"] for row in rows}, {"True"})
            self.assertEqual(
                {row["runner_audit_stage"] for row in rows},
                {summary.RUNNER_STAGE},
            )

    def test_failed_render_leaves_complete_old_output_set_untouched(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            labels = self._four_cases(root)
            output = root / "summary"
            output.mkdir()
            targets = (
                output / "phase_lag_permeability_screen.csv",
                output / "phase_lag_permeability_screen.png",
                output / "phase_lag_permeability_screen.pdf",
                output / "phase_lag_permeability_screen.audit.json",
            )
            old_payloads = {}
            for index, path in enumerate(targets):
                payload = f"old artifact {index}\n".encode()
                path.write_bytes(payload)
                old_payloads[path] = payload
            with mock.patch.object(
                summary, "plot_summary", side_effect=RuntimeError("synthetic render fail")
            ):
                with self.assertRaisesRegex(RuntimeError, "synthetic render fail"):
                    summary.write_summary(root, labels, output, dpi=72)
            self.assertEqual(
                {path: path.read_bytes() for path in targets}, old_payloads
            )
            self.assertFalse(
                list(output.parent.glob(".phase-lag-summary-stage-*"))
            )

    def test_failed_commit_rolls_back_complete_old_output_set(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            labels = self._four_cases(root)
            output = root / "summary"
            output.mkdir()
            targets = (
                output / "phase_lag_permeability_screen.csv",
                output / "phase_lag_permeability_screen.png",
                output / "phase_lag_permeability_screen.pdf",
                output / "phase_lag_permeability_screen.audit.json",
            )
            old_payloads = {}
            for index, path in enumerate(targets):
                payload = f"old commit artifact {index}\n".encode()
                path.write_bytes(payload)
                old_payloads[path] = payload

            real_replace = summary.os.replace
            failed = False

            def fail_once_on_staged_png(source: Path, destination: Path) -> None:
                nonlocal failed
                source_path = Path(source)
                if (
                    not failed
                    and source_path.suffix == ".png"
                    and source_path.parent.name.startswith(
                        ".phase-lag-summary-stage-"
                    )
                ):
                    failed = True
                    raise OSError("synthetic commit fail")
                real_replace(source, destination)

            with mock.patch.object(
                summary.os, "replace", side_effect=fail_once_on_staged_png
            ):
                with self.assertRaisesRegex(
                    summary.SummaryError, "output transaction failed"
                ):
                    summary.write_summary(root, labels, output, dpi=72)
            self.assertTrue(failed)
            self.assertEqual(
                {path: path.read_bytes() for path in targets}, old_payloads
            )

    def test_display_settings_do_not_change_numerical_values(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            labels = self._four_cases(root)
            low_dpi = summary.write_summary(root, labels, root / "low", dpi=72)
            high_dpi = summary.write_summary(root, labels, root / "high", dpi=144)
            self.assertEqual(low_dpi["cases"], high_dpi["cases"])
            self.assertFalse(low_dpi["display_contract"]["smoothing"])
            self.assertFalse(low_dpi["display_contract"]["interpolation"])
            self.assertFalse(
                low_dpi["display_contract"]["display_changes_numerical_values"]
            )
            self.assertNotEqual(
                low_dpi["outputs"]["png"]["sha256"],
                high_dpi["outputs"]["png"]["sha256"],
            )


if __name__ == "__main__":
    unittest.main()
