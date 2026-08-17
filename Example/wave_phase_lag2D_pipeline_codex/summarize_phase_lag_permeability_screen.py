#!/usr/bin/env python3
"""Summarize an audited four-point phase-lag permeability screen.

Only current-runner cases with matching phase-v3 and pressure-only
driver-screen audit files are accepted.  The source simulations and their
result directories are never opened or modified; all plotted values are
revalidated from the raw integrals recorded in the two analysis audit
families, while every driver audit is rebound to its recorded runner audit.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import tempfile
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

import plot_liquid_pressure_phase_lag as phase_plot  # noqa: E402


SCHEMA = "pipeline-phase-lag-permeability-screen-summary-v1"
PHASE_SCHEMA = "pipeline-liquid-pressure-phase-lag-audit-v3"
DRIVER_SCHEMA = "pipeline-phase-lag-driver-screen-v2"
RUNNER_SCHEMA = "pipeline-phase-lag-exploration-runner-v1"
RUNNER_STAGE = "hd_complete"
PHASE_RUNNER_BINDING_SCHEMA = "pipeline-phase-field-runner-input-binding-v1"
EXPECTED_CASE_COUNT = 4
EXPECTED_SATURATION = 0.94
EXPECTED_FRAME_COUNT = 60
EXPECTED_PERIOD_S = 1.3
EXPECTED_COMBINED_WINDOW_S = (1.3, 3.9)
EXPECTED_PERMEABILITIES_M2 = (1.0e-13, 3.0e-13, 1.0e-12, 3.0e-12)
DEFAULT_LABELS = (
    "k1e-13_sw094_nosmooth_fresh",
    "k3e-13_sw094_nosmooth",
    "k1e-12_sw094_nosmooth",
    "k3e-12_sw094_nosmooth",
)
PHASE_AUDIT_SUFFIX = "_liquid_pressure_phase_lag.audit.json"
DRIVER_AUDIT_NAME = "pressure_only_gradient_force.audit.json"
RUNNER_AUDIT_NAME = "runner_audit.json"

EFFECT_METRICS = (
    (
        "net_uplift_force_impulse",
        "combined_net_uplift_effect_percent",
        "Net uplift",
        "N s/m",
    ),
    (
        "local_positive_force_activity_impulse",
        "combined_local_positive_activity_effect_percent",
        "Local + activity",
        "N s/m",
    ),
    (
        "IF_ge_1_area_time",
        "combined_IF_ge_1_area_time_effect_percent",
        r"$IF\geq1$ area-time",
        "m2 s",
    ),
    (
        "shared_HD_joint_area_time",
        "combined_shared_HD_joint_area_time_effect_percent",
        "Shared-HD joint",
        "m2 s",
    ),
)
NET_GATE_CHECKS = (
    "both_post_ramp_cycles_have_positive_net_uplift_delta",
    "combined_net_uplift_impulse_increase_ge_5pct",
)
GATE_CHECKS = (
    "crown_phase_lag_is_resolved",
    "both_post_ramp_cycles_have_positive_signed_force_delta",
    "combined_signed_force_impulse_delta_positive",
    "both_post_ramp_cycles_have_positive_net_uplift_delta",
    "combined_net_uplift_impulse_increase_ge_5pct",
    "combined_IF_ge_1_area_time_direction_positive",
    "combined_shared_HD_joint_area_time_nonzero",
    "combined_shared_HD_joint_area_time_increase_ge_5pct",
    "resolved_joint_advantage_for_two_consecutive_saved_frames",
)
PHASE_QUALIFICATION_CHECKS = (
    "crown_abs_same_column_phase_lag_ge_18deg",
    "crown_amplitude_ratio_ge_0p05",
    "crown_harmonic_r_squared_ge_0p8",
)
MINIMUM_CROWN_PHASE_LAG_DEG = 18.0
MINIMUM_CROWN_AMPLITUDE_RATIO = 0.05
MINIMUM_CROWN_HARMONIC_R_SQUARED = 0.8
MINIMUM_FORCE_IMPULSE_FRACTION = 0.05
MINIMUM_JOINT_AREA_TIME_FRACTION = 0.05
AREA_COMPARISON_ABS_TOL_M2 = 1.0e-12
EXPECTED_WINDOWS_S = ((1.3, 2.6), (2.6, 3.9), (1.3, 3.9))
DRIVER_ANALYZER_PATH = Path(__file__).resolve().with_name(
    "analyze_phase_lag_driver_screen.py"
)
IMPLEMENTATION_DEPENDENCY_PATHS = {
    "gradient_force_module": Path(__file__).resolve().with_name(
        "plot_pipeline_gradient_force.py"
    ),
    "replay_analyzer": Path(__file__).resolve().with_name(
        "analyze_phase_lag_replay_exploration.py"
    ),
    "phase_plot": Path(__file__).resolve().with_name(
        "plot_liquid_pressure_phase_lag.py"
    ),
    "runner": Path(__file__).resolve().with_name("run_phase_lag_exploration.py"),
}
HISTORY_METRICS = (
    (
        "signed_force_impulse",
        "signed_force_n_per_m",
        "N s/m",
    ),
    (
        "net_uplift_force_impulse",
        "net_uplift_force_n_per_m",
        "N s/m",
    ),
    (
        "local_positive_force_activity_impulse",
        "positive_force_n_per_m",
        "N s/m",
    ),
    ("IF_ge_1_area_time", "critical_area_m2", "m2 s"),
    ("shared_HD_joint_area_time", "joint_area_m2", "m2 s"),
)
HISTORY_REQUIRED_FIELDS = (
    "time_s",
    *(
        f"{role}_{metric}"
        for _output, metric, _units in HISTORY_METRICS
        for role in ("lagged", "phase_erased")
    ),
    "lagged_maximum_upward_if",
    "phase_erased_maximum_upward_if",
    "delta_positive_force_n_per_m",
    "delta_signed_force_n_per_m",
    "delta_net_uplift_force_n_per_m",
)


class SummaryError(RuntimeError):
    """Raised when any permeability-screen input is incomplete or inconsistent."""


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _mapping(value: Any, context: str) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise SummaryError(f"{context} must be a JSON object")
    return value


def _finite_float(value: Any, context: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise SummaryError(f"{context} must be a finite number")
    result = float(value)
    if not math.isfinite(result):
        raise SummaryError(f"{context} must be a finite number")
    return result


def _optional_finite_float(value: Any, context: str) -> float | None:
    if value is None:
        return None
    return _finite_float(value, context)


def _finite_csv_float(value: Any, context: str) -> float:
    if not isinstance(value, str) or not value.strip():
        raise SummaryError(f"{context} must be a finite CSV number")
    try:
        result = float(value)
    except ValueError as error:
        raise SummaryError(f"{context} must be a finite CSV number") from error
    if not math.isfinite(result):
        raise SummaryError(f"{context} must be a finite CSV number")
    return result


def _positive_int(value: Any, context: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
        raise SummaryError(f"{context} must be a positive integer")
    return value


def _explicit_bool(value: Any, context: str) -> bool:
    if type(value) is not bool:
        raise SummaryError(f"{context} must be an explicit JSON boolean")
    return value


def _sha256_text(value: Any, context: str) -> str:
    if (
        not isinstance(value, str)
        or len(value) != 64
        or any(character not in "0123456789abcdef" for character in value)
    ):
        raise SummaryError(f"{context} must be a lowercase SHA-256")
    return value


def _validated_file_record(
    value: Any,
    context: str,
    *,
    expected_path: Path | None = None,
) -> dict[str, str]:
    """Verify one absolute path/hash record against the bytes on disk."""

    record = _mapping(value, context)
    if set(record) != {"path", "sha256"}:
        raise SummaryError(f"{context} must contain exactly path and sha256")
    path_value = record.get("path")
    if not isinstance(path_value, str) or not path_value:
        raise SummaryError(f"{context}.path must be an absolute path")
    path = Path(path_value)
    if not path.is_absolute():
        raise SummaryError(f"{context}.path must be an absolute path")
    resolved = path.resolve()
    if expected_path is not None and resolved != expected_path.resolve():
        raise SummaryError(f"{context}.path does not bind the expected artifact")
    recorded_hash = _sha256_text(record.get("sha256"), f"{context}.sha256")
    if not resolved.is_file():
        raise SummaryError(f"Missing {context} artifact: {resolved}")
    if sha256(resolved) != recorded_hash:
        raise SummaryError(f"{context} artifact hash differs")
    return {"path": str(resolved), "sha256": recorded_hash}


def _validated_path_and_hash(
    path_value: Any,
    hash_value: Any,
    context: str,
    *,
    expected_path: Path,
) -> dict[str, str]:
    return _validated_file_record(
        {"path": path_value, "sha256": hash_value},
        context,
        expected_path=expected_path,
    )


def _validated_sized_file_record(
    value: Any,
    context: str,
    *,
    expected_path: Path,
) -> dict[str, str | int]:
    record = _mapping(value, context)
    if set(record) != {"path", "size_bytes", "sha256"}:
        raise SummaryError(
            f"{context} must contain exactly path, size_bytes, and sha256"
        )
    validated = _validated_path_and_hash(
        record.get("path"),
        record.get("sha256"),
        context,
        expected_path=expected_path,
    )
    size = _positive_int(record.get("size_bytes"), f"{context}.size_bytes")
    if Path(validated["path"]).stat().st_size != size:
        raise SummaryError(f"{context} artifact size differs")
    return {**validated, "size_bytes": size}


def _same_float(
    left: float,
    right: float,
    context: str,
    *,
    relative_tolerance: float = 1.0e-12,
    absolute_tolerance: float = 1.0e-12,
) -> None:
    if not math.isclose(
        left,
        right,
        rel_tol=relative_tolerance,
        abs_tol=absolute_tolerance,
    ):
        raise SummaryError(f"{context} differs: {left!r} versus {right!r}")


def _safe_label(label: str) -> str:
    if not label or Path(label).name != label or label in {".", ".."}:
        raise SummaryError(f"Unsafe permeability-screen label: {label!r}")
    return label


def _read_audit(path: Path, expected_schema: str, context: str) -> dict[str, Any]:
    if not path.is_file():
        raise SummaryError(
            f"Missing {context}: {path}. Legacy or substitute audits are not accepted."
        )
    try:
        audit = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise SummaryError(f"Cannot read {context}: {path}") from error
    audit = _mapping(audit, context)
    if audit.get("schema") != expected_schema:
        raise SummaryError(
            f"{context} schema must be {expected_schema!r}, got "
            f"{audit.get('schema')!r}"
        )
    return audit


def _phase_audit_path(root: Path, label: str) -> Path:
    phase_dir = root / "analysis/phase_lag_exploratory" / label / "phase_v3"
    matches = sorted(phase_dir.glob(f"*{PHASE_AUDIT_SUFFIX}"))
    if len(matches) != 1:
        raise SummaryError(
            f"Expected exactly one phase-v3 audit for {label}, found "
            f"{len(matches)}. Legacy or substitute audits are not accepted."
        )
    return matches[0]


def _driver_audit_path(root: Path, label: str) -> Path:
    return (
        root
        / "analysis/phase_lag_exploratory"
        / label
        / "driver_screen"
        / DRIVER_AUDIT_NAME
    )


def _runner_audit_path(root: Path, label: str) -> Path:
    return root / "analysis/phase_lag_exploratory" / label / RUNNER_AUDIT_NAME


def _validated_runner_source(
    root: Path,
    label: str,
    driver_audit: dict[str, Any],
    driver_parameters: dict[str, Any],
) -> dict[str, Any]:
    """Rebind a driver audit to the exact current-runner audit it recorded."""

    record = _mapping(
        driver_audit.get("runner_audit"), f"{label} driver runner_audit"
    )
    recorded_path_value = record.get("path")
    if not isinstance(recorded_path_value, str) or not recorded_path_value:
        raise SummaryError(f"{label} driver runner_audit.path must be absolute")
    recorded_path = Path(recorded_path_value)
    if not recorded_path.is_absolute():
        raise SummaryError(f"{label} driver runner_audit.path must be absolute")

    expected_path = _runner_audit_path(root, label)
    if not expected_path.is_file():
        raise SummaryError(f"Missing {label} current runner audit: {expected_path}")
    expected_path = expected_path.resolve()
    if recorded_path.resolve() != expected_path:
        raise SummaryError(
            f"{label} driver runner_audit.path does not bind the label's "
            f"current runner audit"
        )

    recorded_hash = record.get("sha256")
    if (
        not isinstance(recorded_hash, str)
        or len(recorded_hash) != 64
        or any(character not in "0123456789abcdef" for character in recorded_hash)
    ):
        raise SummaryError(
            f"{label} driver runner_audit.sha256 must be a lowercase SHA-256"
        )
    try:
        payload = expected_path.read_bytes()
    except OSError as error:
        raise SummaryError(
            f"Cannot read {label} current runner audit: {expected_path}"
        ) from error
    observed_hash = hashlib.sha256(payload).hexdigest()
    if observed_hash != recorded_hash:
        raise SummaryError(
            f"{label} current runner audit hash differs from the driver binding"
        )
    try:
        runner_audit = json.loads(payload.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise SummaryError(
            f"Cannot read {label} current runner audit: {expected_path}"
        ) from error
    runner_audit = _mapping(runner_audit, f"{label} current runner audit")
    if runner_audit.get("schema") != RUNNER_SCHEMA:
        raise SummaryError(
            f"{label} current runner audit schema must be {RUNNER_SCHEMA!r}, got "
            f"{runner_audit.get('schema')!r}"
        )
    if runner_audit.get("stage") != RUNNER_STAGE:
        raise SummaryError(
            f"{label} current runner audit stage must be {RUNNER_STAGE!r}, got "
            f"{runner_audit.get('stage')!r}"
        )
    if runner_audit.get("label") != label:
        raise SummaryError(f"{label} current runner audit label identity differs")

    runner_parameters = _mapping(
        runner_audit.get("parameters"), f"{label} current runner parameters"
    )
    runner_permeability = _finite_float(
        runner_parameters.get("intrinsic_permeability_m2"),
        f"{label} current runner intrinsic_permeability_m2",
    )
    runner_saturation = _finite_float(
        runner_parameters.get("liquid_saturation"),
        f"{label} current runner liquid_saturation",
    )
    runner_smoothing = _explicit_bool(
        runner_parameters.get("pressure_smoothing"),
        f"{label} current runner pressure_smoothing",
    )
    _same_float(
        runner_permeability,
        float(driver_parameters["intrinsic_permeability_m2"]),
        f"{label} runner/driver parameter intrinsic_permeability_m2",
        absolute_tolerance=0.0,
    )
    _same_float(
        runner_saturation,
        float(driver_parameters["liquid_saturation"]),
        f"{label} runner/driver parameter liquid_saturation",
        relative_tolerance=0.0,
    )
    if runner_smoothing is not driver_parameters["pressure_smoothing"]:
        raise SummaryError(
            f"{label} runner/driver parameter pressure_smoothing differs"
        )

    return {
        "path": str(expected_path),
        "sha256": observed_hash,
        "schema": RUNNER_SCHEMA,
        "stage": RUNNER_STAGE,
        "_audit": runner_audit,
    }


def _validate_phase_runner_input_binding(
    root: Path,
    label: str,
    phase_audit: dict[str, Any],
    driver_audit: dict[str, Any],
    runner_audit: dict[str, Any],
) -> None:
    """Recheck every phase input against both the driver and runner records."""

    try:
        stages = _mapping(runner_audit.get("stages"), f"{label} runner stages")
        if set(stages) != {"eq", "hd"}:
            raise SummaryError(f"{label} runner stages must contain exactly eq and hd")
        eq_stage = _mapping(stages["eq"], f"{label} runner EQ stage")
        hd_stage = _mapping(stages["hd"], f"{label} runner HD stage")
        eq_audit = _mapping(eq_stage.get("audit"), f"{label} runner EQ audit")
        hd_audit = _mapping(hd_stage.get("audit"), f"{label} runner HD audit")
        particle_input = _mapping(
            eq_audit["stability_qa"]["normal_stress"]["inputs"]["particles"],
            f"{label} runner EQ particle input",
        )
    except (KeyError, TypeError) as error:
        raise SummaryError(f"{label} runner EQ/HD input provenance is incomplete") from error

    config_dir = root / "configs/phase_lag_exploratory" / label
    eq_config = _validated_sized_file_record(
        eq_stage.get("config"),
        f"{label} runner EQ config",
        expected_path=config_dir / "01_EQ.json",
    )
    hd_config = _validated_sized_file_record(
        hd_stage.get("config"),
        f"{label} runner HD config",
        expected_path=config_dir / "03_HD.json",
    )
    checkpoint_record = _mapping(
        eq_audit.get("final_vtp"), f"{label} runner EQ checkpoint"
    )
    checkpoint_path = Path(str(checkpoint_record.get("path", "")))
    if not checkpoint_path.is_absolute():
        raise SummaryError(f"{label} runner EQ checkpoint path must be absolute")
    checkpoint = _validated_sized_file_record(
        checkpoint_record,
        f"{label} runner EQ checkpoint",
        expected_path=checkpoint_path,
    )

    raw_particle_path = Path(str(particle_input.get("path", "")))
    particle_path = (
        raw_particle_path if raw_particle_path.is_absolute() else root / raw_particle_path
    ).resolve()
    particles = _validated_sized_file_record(
        {
            "path": str(particle_path),
            "size_bytes": particle_input.get("size_bytes"),
            "sha256": particle_input.get("sha256"),
        },
        f"{label} runner EQ particles",
        expected_path=particle_path,
    )

    result_directory_value = hd_audit.get("result_directory")
    if not isinstance(result_directory_value, str) or not result_directory_value:
        raise SummaryError(f"{label} runner HD result_directory is invalid")
    result_directory = Path(result_directory_value)
    if (
        not result_directory.is_absolute()
        or result_directory_value != str(result_directory.resolve())
        or not result_directory.is_dir()
    ):
        raise SummaryError(
            f"{label} runner HD result_directory must be canonical and existing"
        )
    raw_frames = hd_audit.get("particle_vtp")
    if not isinstance(raw_frames, list) or len(raw_frames) != EXPECTED_FRAME_COUNT:
        raise SummaryError(
            f"{label} runner HD must bind exactly {EXPECTED_FRAME_COUNT} particle VTPs"
        )
    frames: list[dict[str, str | int]] = []
    for index, raw_record in enumerate(raw_frames):
        record = _mapping(raw_record, f"{label} runner HD VTP {index}")
        frame_path = Path(str(record.get("path", "")))
        if not frame_path.is_absolute():
            raise SummaryError(f"{label} runner HD VTP {index} path must be absolute")
        frames.append(
            _validated_sized_file_record(
                record,
                f"{label} runner HD VTP {index}",
                expected_path=frame_path,
            )
        )
    if len({record["path"] for record in frames}) != EXPECTED_FRAME_COUNT:
        raise SummaryError(f"{label} runner HD VTP paths must be unique")
    if any(Path(str(record["path"])).parent != result_directory for record in frames):
        raise SummaryError(f"{label} runner HD VTP escaped result_directory")

    history_record = _mapping(
        hd_audit.get("pipeline_history"), f"{label} runner HD pipeline history"
    )
    history_path = Path(str(history_record.get("path", "")))
    if not history_path.is_absolute():
        raise SummaryError(f"{label} runner HD pipeline history path must be absolute")
    pipeline_history = _validated_sized_file_record(
        history_record,
        f"{label} runner HD pipeline history",
        expected_path=history_path,
    )

    expected_binding = {
        "schema": PHASE_RUNNER_BINDING_SCHEMA,
        "eq": {
            "checkpoint_config": eq_config,
            "checkpoint": checkpoint,
            "particles": particles,
        },
        "hd": {
            "analysis_config": hd_config,
            "result_directory": str(result_directory),
            "result_particle_vtp": frames,
            "pipeline_history": pipeline_history,
        },
    }
    try:
        phase_inputs = phase_audit["artifacts"]["inputs"]
    except (KeyError, TypeError) as error:
        raise SummaryError(f"{label} phase-v3 input provenance is incomplete") from error
    expected_phase_inputs = {
        "schema": phase_plot.INPUT_SCHEMA,
        "analysis_config": hd_config,
        "checkpoint_config": eq_config,
        "checkpoint": checkpoint,
        "particles": particles,
        "pipeline_history": pipeline_history,
        "result_directory": str(result_directory),
        "result_particle_vtp": frames,
    }
    if phase_inputs != expected_phase_inputs:
        raise SummaryError(
            f"{label} phase-v3 inputs do not exactly bind runner EQ/HD artifacts"
        )
    phase_field = _mapping(
        driver_audit.get("phase_field"), f"{label} driver phase_field"
    )
    if phase_field.get("runner_binding") != expected_binding:
        raise SummaryError(
            f"{label} driver phase-field binding differs from runner EQ/HD artifacts"
        )


def _validated_driver_file_provenance(
    root: Path,
    label: str,
    driver_audit: dict[str, Any],
    phase_path: Path,
) -> dict[str, Any]:
    """Reverify every driver-recorded input/output needed by the summary."""

    analyzer = _validated_file_record(
        driver_audit.get("analyzer"),
        f"{label} driver analyzer",
        expected_path=DRIVER_ANALYZER_PATH,
    )

    dependency_records = _mapping(
        driver_audit.get("implementation_dependencies"),
        f"{label} implementation_dependencies",
    )
    if set(dependency_records) != set(IMPLEMENTATION_DEPENDENCY_PATHS):
        raise SummaryError(
            f"{label} implementation_dependencies must contain exactly "
            f"{sorted(IMPLEMENTATION_DEPENDENCY_PATHS)}"
        )
    dependencies = {
        name: _validated_file_record(
            dependency_records[name],
            f"{label} implementation dependency {name}",
            expected_path=expected_path,
        )
        for name, expected_path in IMPLEMENTATION_DEPENDENCY_PATHS.items()
    }

    config_records = _mapping(driver_audit.get("configs"), f"{label} configs")
    if set(config_records) != {"HD", "RL", "RE"}:
        raise SummaryError(f"{label} configs must contain exactly HD, RL, and RE")
    config_dir = root / "configs/phase_lag_exploratory" / label
    configs = {
        role: _validated_file_record(
            config_records[role],
            f"{label} {role} config",
            expected_path=config_dir / filename,
        )
        for role, filename in (
            ("HD", "03_HD.json"),
            ("RL", "04_RL.json"),
            ("RE", "05_RE.json"),
        )
    }

    database_records = _mapping(
        driver_audit.get("pressure_databases"), f"{label} pressure_databases"
    )
    if set(database_records) != {"lagged", "phase_erased"}:
        raise SummaryError(
            f"{label} pressure_databases must contain exactly lagged and phase_erased"
        )
    database_root = root / "pressure_databases/phase_lag_exploratory" / label
    expected_database_paths = {
        "lagged": {
            "points": database_root / "lagged/pressure_points.txt",
            "values": database_root / "lagged/pressure_values.bin",
        },
        "phase_erased": {
            "points": database_root / "phase_erased/phase_erased_points.txt",
            "values": database_root / "phase_erased/phase_erased_values.bin",
        },
    }
    databases: dict[str, Any] = {}
    for role in ("lagged", "phase_erased"):
        record = _mapping(database_records[role], f"{label} {role} database")
        for required in (
            "points",
            "points_sha256",
            "values",
            "values_sha256",
            "step_interval",
            "max_step",
            "source_dt_s",
        ):
            if required not in record:
                raise SummaryError(f"{label} {role} database omits {required}")
        databases[role] = {
            "points": _validated_path_and_hash(
                record.get("points"),
                record.get("points_sha256"),
                f"{label} {role} database points",
                expected_path=expected_database_paths[role]["points"],
            ),
            "values": _validated_path_and_hash(
                record.get("values"),
                record.get("values_sha256"),
                f"{label} {role} database values",
                expected_path=expected_database_paths[role]["values"],
            ),
            "step_interval": _positive_int(
                record.get("step_interval"), f"{label} {role} step_interval"
            ),
            "max_step": _positive_int(
                record.get("max_step"), f"{label} {role} max_step"
            ),
            "source_dt_s": _finite_float(
                record.get("source_dt_s"), f"{label} {role} source_dt_s"
            ),
        }
    for field in ("step_interval", "max_step", "source_dt_s"):
        _same_float(
            float(databases["lagged"][field]),
            float(databases["phase_erased"][field]),
            f"{label} pressure database {field}",
            relative_tolerance=0.0,
            absolute_tolerance=1.0e-15,
        )

    driver_dir = root / "analysis/phase_lag_exploratory" / label / "driver_screen"
    history = _validated_path_and_hash(
        driver_audit.get("history_csv"),
        driver_audit.get("history_csv_sha256"),
        f"{label} driver history",
        expected_path=driver_dir / "pressure_only_gradient_force_history.csv",
    )
    figure_record = _mapping(driver_audit.get("figure"), f"{label} driver figure")
    figures = {
        suffix: _validated_file_record(
            figure_record.get(suffix),
            f"{label} driver figure {suffix}",
            expected_path=driver_dir / f"pressure_only_gradient_force.{suffix}",
        )
        for suffix in ("png", "pdf")
    }

    phase_field = _mapping(
        driver_audit.get("phase_field"), f"{label} phase_field"
    )
    phase_field_record = _validated_path_and_hash(
        phase_field.get("path"),
        phase_field.get("sha256"),
        f"{label} driver phase_field",
        expected_path=phase_path,
    )

    transform_record = _mapping(
        driver_audit.get("phase_transform"), f"{label} phase_transform"
    )
    if set(transform_record) != {
        "path",
        "sha256",
        "surface_reference_source",
        "metadata",
    }:
        raise SummaryError(
            f"{label} phase_transform must contain exactly path, sha256, "
            "surface_reference_source, and metadata"
        )
    transform_path = database_root / "phase_erased/phase_erased_metadata.json"
    transform = _validated_path_and_hash(
        transform_record.get("path"),
        transform_record.get("sha256"),
        f"{label} phase transform metadata",
        expected_path=transform_path,
    )
    try:
        observed_metadata = json.loads(transform_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise SummaryError(f"Cannot read {label} phase transform metadata") from error
    recorded_metadata = _mapping(
        transform_record.get("metadata"), f"{label} phase transform metadata payload"
    )
    if observed_metadata != recorded_metadata:
        raise SummaryError(
            f"{label} phase transform embedded metadata differs from its file"
        )
    if recorded_metadata.get("schema") != "phase-erased-pressure-control-v1":
        raise SummaryError(f"{label} phase transform schema is not current")
    for group, role in (("source", "lagged"), ("output", "phase_erased")):
        group_record = _mapping(
            recorded_metadata.get(group), f"{label} phase transform {group}"
        )
        for kind in ("points", "values"):
            candidate = Path(str(group_record.get(kind, "")))
            if not candidate.is_absolute():
                candidate = root / candidate
            expected = Path(databases[role][kind]["path"])
            if candidate.resolve() != expected:
                raise SummaryError(
                    f"{label} phase transform {group} {kind} path differs"
                )
            if group_record.get(f"{kind}_sha256") != databases[role][kind]["sha256"]:
                raise SummaryError(
                    f"{label} phase transform {group} {kind} hash differs"
                )
    try:
        surface_source = recorded_metadata["transform"]["surface_reference"]["source"]
    except (KeyError, TypeError) as error:
        raise SummaryError(
            f"{label} phase transform surface-reference dependency is incomplete"
        ) from error
    surface_reference = _validated_sized_file_record(
        transform_record.get("surface_reference_source"),
        f"{label} recorded phase transform surface-reference source",
        expected_path=root / "top_surface_traction_particle_id.txt",
    )
    embedded_surface_reference = _validated_file_record(
        surface_source,
        f"{label} embedded phase transform surface-reference dependency",
        expected_path=root / "top_surface_traction_particle_id.txt",
    )
    if (
        embedded_surface_reference["path"] != surface_reference["path"]
        or embedded_surface_reference["sha256"] != surface_reference["sha256"]
    ):
        raise SummaryError(
            f"{label} recorded and embedded surface-reference bindings differ"
        )

    return {
        "analyzer": analyzer,
        "implementation_dependencies": dependencies,
        "configs": configs,
        "pressure_databases": databases,
        "history_csv": history,
        "figures": figures,
        "phase_field": phase_field_record,
        "phase_transform": transform,
        "surface_reference": surface_reference,
    }


def _validated_parameters(audit: dict[str, Any], context: str) -> dict[str, Any]:
    parameters = _mapping(audit.get("parameters"), f"{context}.parameters")
    permeability = _finite_float(
        parameters.get("intrinsic_permeability_m2"),
        f"{context}.parameters.intrinsic_permeability_m2",
    )
    saturation = _finite_float(
        parameters.get("liquid_saturation"),
        f"{context}.parameters.liquid_saturation",
    )
    gas_saturation = _finite_float(
        parameters.get("gas_saturation"),
        f"{context}.parameters.gas_saturation",
    )
    if permeability <= 0.0:
        raise SummaryError(f"{context} permeability must be positive")
    _same_float(
        saturation,
        EXPECTED_SATURATION,
        f"{context} liquid saturation",
        relative_tolerance=0.0,
    )
    _same_float(
        saturation + gas_saturation,
        1.0,
        f"{context} saturation closure",
        relative_tolerance=0.0,
    )
    output_smoothing = _explicit_bool(
        parameters.get("pressure_smoothing"),
        f"{context}.parameters.pressure_smoothing",
    )
    in_loop_smoothing = _explicit_bool(
        parameters.get("pressure_smoothing_in_loop"),
        f"{context}.parameters.pressure_smoothing_in_loop",
    )
    fixed = _explicit_bool(
        parameters.get("pipeline_fixed"),
        f"{context}.parameters.pipeline_fixed",
    )
    if output_smoothing is not False or in_loop_smoothing is not False:
        raise SummaryError(f"{context} must have both pressure-smoothing flags false")
    if fixed is not True:
        raise SummaryError(f"{context} must have an explicitly fixed pipeline")
    return {
        "intrinsic_permeability_m2": permeability,
        "liquid_saturation": saturation,
        "gas_saturation": gas_saturation,
        "pressure_smoothing": output_smoothing,
        "pressure_smoothing_in_loop": in_loop_smoothing,
        "pipeline_fixed": fixed,
    }


def _validate_parameter_pair(
    phase: dict[str, Any], driver: dict[str, Any], label: str
) -> dict[str, Any]:
    for key in (
        "intrinsic_permeability_m2",
        "liquid_saturation",
        "gas_saturation",
    ):
        _same_float(
            float(phase[key]),
            float(driver[key]),
            f"{label} phase/driver parameter {key}",
            absolute_tolerance=(
                0.0 if key == "intrinsic_permeability_m2" else 1.0e-12
            ),
        )
    for key in (
        "pressure_smoothing",
        "pressure_smoothing_in_loop",
        "pipeline_fixed",
    ):
        if phase[key] is not driver[key]:
            raise SummaryError(f"{label} phase/driver parameter {key} differs")
    return phase


def _validated_phase_metrics(
    phase_audit: dict[str, Any], label: str
) -> dict[str, Any]:
    if str(phase_audit.get("case", "")).upper() != label.upper():
        raise SummaryError(f"{label} phase-v3 case identity differs")
    if _explicit_bool(
        phase_audit.get("require_unsmoothed"),
        f"{label} phase-v3 require_unsmoothed",
    ) is not True:
        raise SummaryError(f"{label} phase-v3 was not run with --require-unsmoothed")
    if _positive_int(
        phase_audit.get("frame_count"), f"{label} phase-v3 frame_count"
    ) != EXPECTED_FRAME_COUNT:
        raise SummaryError(f"{label} phase-v3 must bind exactly 60 saved frames")
    particle_count = _positive_int(
        phase_audit.get("particle_count"), f"{label} phase-v3 particle_count"
    )
    period = _finite_float(
        phase_audit.get("wave_period_s"), f"{label} phase-v3 wave_period_s"
    )
    _same_float(period, EXPECTED_PERIOD_S, f"{label} phase-v3 wave period")
    fit_window = phase_audit.get("fit_window_s")
    if not isinstance(fit_window, list) or len(fit_window) != 2:
        raise SummaryError(f"{label} phase-v3 fit_window_s is incomplete")
    for observed, expected in zip(fit_window, EXPECTED_COMBINED_WINDOW_S):
        _same_float(
            _finite_float(observed, f"{label} phase-v3 fit window"),
            expected,
            f"{label} phase-v3 fit window",
        )

    probes = _mapping(phase_audit.get("probe_results"), f"{label} probe_results")
    surface = _mapping(probes.get("surface"), f"{label} surface probe")
    crown = _mapping(probes.get("crown"), f"{label} crown probe")
    surface_amplitude = _finite_float(
        surface.get("fundamental_amplitude_pa"), f"{label} surface amplitude"
    )
    crown_amplitude = _finite_float(
        crown.get("fundamental_amplitude_pa"), f"{label} crown amplitude"
    )
    if surface_amplitude <= 0.0 or crown_amplitude < 0.0:
        raise SummaryError(f"{label} phase amplitudes are outside their valid range")
    crown_phase = _optional_finite_float(
        crown.get("phase_difference_from_same_column_surface_deg"),
        f"{label} crown same-column phase",
    )
    crown_r_squared = _finite_float(
        crown.get("harmonic_r_squared"), f"{label} crown harmonic R-squared"
    )
    if crown_r_squared < 0.0 or crown_r_squared > 1.0 + 1.0e-12:
        raise SummaryError(f"{label} crown harmonic R-squared is outside [0,1]")
    phase_present = _explicit_bool(
        phase_audit.get("phase_present"), f"{label} phase_present"
    )
    amplitude_ratio = crown_amplitude / surface_amplitude
    qualification_checks = {
        "crown_abs_same_column_phase_lag_ge_18deg": (
            crown_phase is not None
            and abs(crown_phase) >= MINIMUM_CROWN_PHASE_LAG_DEG
        ),
        "crown_amplitude_ratio_ge_0p05": (
            amplitude_ratio >= MINIMUM_CROWN_AMPLITUDE_RATIO
        ),
        "crown_harmonic_r_squared_ge_0p8": (
            crown_r_squared >= MINIMUM_CROWN_HARMONIC_R_SQUARED
        ),
    }
    qualification = {
        "passed": all(qualification_checks.values()),
        "checks": qualification_checks,
        "observed": {
            "crown_same_column_phase_lag_deg": crown_phase,
            "crown_amplitude_ratio": amplitude_ratio,
            "crown_harmonic_r_squared": crown_r_squared,
        },
    }
    return {
        "particle_count": particle_count,
        "crown_phase_difference_same_column_deg": crown_phase,
        "crown_phase_lag_magnitude_deg": (
            None if crown_phase is None else abs(crown_phase)
        ),
        "crown_fundamental_amplitude_pa": crown_amplitude,
        "surface_fundamental_amplitude_pa": surface_amplitude,
        "crown_to_surface_amplitude_percent": 100.0 * amplitude_ratio,
        "crown_harmonic_r_squared": crown_r_squared,
        "phase_present": phase_present,
        "phase_qualification": qualification,
    }


def _relative_fraction(lagged: float, erased: float) -> float | None:
    return None if erased <= 0.0 else lagged / erased - 1.0


def _read_raw_history(path: Path, label: str) -> list[dict[str, float]]:
    """Read the hashed raw history and validate its physical bookkeeping."""

    try:
        with path.open(newline="", encoding="utf-8") as stream:
            reader = csv.DictReader(stream)
            if reader.fieldnames is None or len(reader.fieldnames) != len(
                set(reader.fieldnames)
            ):
                raise SummaryError(f"{label} driver history header is invalid")
            if set(reader.fieldnames) != set(HISTORY_REQUIRED_FIELDS):
                raise SummaryError(
                    f"{label} driver history fields differ from the raw contract"
                )
            raw_rows = list(reader)
    except (OSError, csv.Error) as error:
        raise SummaryError(f"Cannot read {label} driver history") from error
    if len(raw_rows) != EXPECTED_FRAME_COUNT:
        raise SummaryError(
            f"{label} driver history must contain exactly {EXPECTED_FRAME_COUNT} rows"
        )
    rows: list[dict[str, float]] = []
    for index, raw_row in enumerate(raw_rows):
        row = {
            field: _finite_csv_float(
                raw_row.get(field), f"{label} history row {index} {field}"
            )
            for field in HISTORY_REQUIRED_FIELDS
        }
        for role in ("lagged", "phase_erased"):
            signed = row[f"{role}_signed_force_n_per_m"]
            net = row[f"{role}_net_uplift_force_n_per_m"]
            positive = row[f"{role}_positive_force_n_per_m"]
            _same_float(
                net,
                max(signed, 0.0),
                f"{label} history row {index} {role} net/signed identity",
                relative_tolerance=1.0e-10,
                absolute_tolerance=1.0e-12,
            )
            if positive < net - 1.0e-12:
                raise SummaryError(
                    f"{label} history row {index} local positive activity is below net uplift"
                )
            for metric in ("critical_area_m2", "joint_area_m2"):
                if row[f"{role}_{metric}"] < 0.0:
                    raise SummaryError(
                        f"{label} history row {index} {role}_{metric} is negative"
                    )
        for delta_name, metric in (
            ("delta_positive_force_n_per_m", "positive_force_n_per_m"),
            ("delta_signed_force_n_per_m", "signed_force_n_per_m"),
            ("delta_net_uplift_force_n_per_m", "net_uplift_force_n_per_m"),
        ):
            _same_float(
                row[delta_name],
                row[f"lagged_{metric}"] - row[f"phase_erased_{metric}"],
                f"{label} history row {index} {delta_name}",
                relative_tolerance=1.0e-10,
                absolute_tolerance=1.0e-12,
            )
        rows.append(row)
    times = np.asarray([row["time_s"] for row in rows], dtype=np.float64)
    if not np.all(np.diff(times) > 0.0):
        raise SummaryError(f"{label} driver history times are not strictly increasing")
    return rows


def _trapezoid(values: np.ndarray, times: np.ndarray) -> float:
    return float(np.sum(0.5 * (values[:-1] + values[1:]) * np.diff(times)))


def _integrate_raw_window(
    rows: list[dict[str, float]], start_s: float, end_s: float, label: str
) -> dict[str, Any]:
    selected = [
        row
        for row in rows
        if start_s - 1.0e-12 <= row["time_s"] <= end_s + 1.0e-12
    ]
    if len(selected) < 2:
        raise SummaryError(f"{label} raw history does not resolve {start_s}--{end_s} s")
    times = np.asarray([row["time_s"] for row in selected], dtype=np.float64)
    if not math.isclose(float(times[0]), start_s, abs_tol=1.0e-12) or not math.isclose(
        float(times[-1]), end_s, abs_tol=1.0e-12
    ):
        raise SummaryError(
            f"{label} raw history omits a boundary of {start_s}--{end_s} s"
        )
    window: dict[str, Any] = {
        "start_s": start_s,
        "end_s": end_s,
        "saved_frames": len(selected),
    }
    for output_name, raw_metric, units in HISTORY_METRICS:
        lagged = _trapezoid(
            np.asarray(
                [row[f"lagged_{raw_metric}"] for row in selected],
                dtype=np.float64,
            ),
            times,
        )
        erased = _trapezoid(
            np.asarray(
                [row[f"phase_erased_{raw_metric}"] for row in selected],
                dtype=np.float64,
            ),
            times,
        )
        window[output_name] = {
            "units": units,
            "lagged": lagged,
            "phase_erased": erased,
            "lagged_minus_phase_erased": lagged - erased,
            "lagged_minus_phase_erased_fraction": _relative_fraction(
                lagged, erased
            ),
        }
    return window


def _validate_integral_record(
    recorded: Any,
    recomputed: dict[str, Any],
    context: str,
) -> None:
    record = _mapping(recorded, context)
    expected_keys = {
        "units",
        "lagged",
        "phase_erased",
        "lagged_minus_phase_erased",
        "lagged_minus_phase_erased_fraction",
    }
    if set(record) != expected_keys:
        raise SummaryError(f"{context} integral record is incomplete")
    if record.get("units") != recomputed["units"]:
        raise SummaryError(f"{context} units differ from the registered definition")
    for key in ("lagged", "phase_erased", "lagged_minus_phase_erased"):
        _same_float(
            _finite_float(record.get(key), f"{context}.{key}"),
            float(recomputed[key]),
            f"{context}.{key} raw-history recomputation",
            relative_tolerance=1.0e-10,
            absolute_tolerance=1.0e-12,
        )
    recorded_fraction = record.get("lagged_minus_phase_erased_fraction")
    recomputed_fraction = recomputed["lagged_minus_phase_erased_fraction"]
    if recomputed_fraction is None:
        if recorded_fraction is not None:
            raise SummaryError(
                f"{context} fraction must be null for a non-positive denominator"
            )
    else:
        _same_float(
            _finite_float(recorded_fraction, f"{context}.fraction"),
            float(recomputed_fraction),
            f"{context}.fraction raw-history recomputation",
            relative_tolerance=1.0e-10,
            absolute_tolerance=1.0e-12,
        )


def _validated_windows(
    driver_audit: dict[str, Any], rows: list[dict[str, float]], label: str
) -> list[dict[str, Any]]:
    recorded = driver_audit.get("integration_windows")
    if not isinstance(recorded, list) or len(recorded) != len(EXPECTED_WINDOWS_S):
        raise SummaryError(f"{label} driver screen must contain exactly three windows")
    recomputed_windows = [
        _integrate_raw_window(rows, start, end, label)
        for start, end in EXPECTED_WINDOWS_S
    ]
    for index, (recorded_value, recomputed) in enumerate(
        zip(recorded, recomputed_windows)
    ):
        record = _mapping(recorded_value, f"{label} integration_windows[{index}]")
        expected_keys = {
            "start_s",
            "end_s",
            "saved_frames",
            *(output for output, _raw, _units in HISTORY_METRICS),
        }
        if set(record) != expected_keys:
            raise SummaryError(
                f"{label} integration_windows[{index}] is incomplete or has extra fields"
            )
        for key in ("start_s", "end_s"):
            _same_float(
                _finite_float(record.get(key), f"{label} window {index} {key}"),
                float(recomputed[key]),
                f"{label} window {index} {key}",
            )
        if record.get("saved_frames") != recomputed["saved_frames"]:
            raise SummaryError(
                f"{label} window {index} saved-frame count differs from raw history"
            )
        for output_name, _raw_metric, _units in HISTORY_METRICS:
            _validate_integral_record(
                record.get(output_name),
                recomputed[output_name],
                f"{label} window {index} {output_name}",
            )
    return recomputed_windows


def _validate_phase_qualification_record(
    recorded: Any, recomputed: dict[str, Any], context: str
) -> None:
    record = _mapping(recorded, context)
    if set(record) != {"passed", "checks", "observed"}:
        raise SummaryError(f"{context} must contain passed, checks, and observed")
    if _explicit_bool(record.get("passed"), f"{context}.passed") is not recomputed[
        "passed"
    ]:
        raise SummaryError(f"{context}.passed differs from phase-v3 recomputation")
    checks = _mapping(record.get("checks"), f"{context}.checks")
    if set(checks) != set(PHASE_QUALIFICATION_CHECKS):
        raise SummaryError(f"{context}.checks is not the complete registered set")
    for name in PHASE_QUALIFICATION_CHECKS:
        if _explicit_bool(checks.get(name), f"{context}.{name}") is not recomputed[
            "checks"
        ][name]:
            raise SummaryError(f"{context} check {name} differs from recomputation")
    observed = _mapping(record.get("observed"), f"{context}.observed")
    if set(observed) != set(recomputed["observed"]):
        raise SummaryError(f"{context}.observed is incomplete")
    for name, expected in recomputed["observed"].items():
        value = _optional_finite_float(observed.get(name), f"{context}.observed.{name}")
        if expected is None:
            if value is not None:
                raise SummaryError(f"{context}.observed.{name} must be null")
        elif value is None:
            raise SummaryError(f"{context}.observed.{name} must be finite")
        else:
            _same_float(value, float(expected), f"{context}.observed.{name}")


def _validated_driver_metrics(
    driver_audit: dict[str, Any],
    label: str,
    phase_particle_count: int,
    phase_qualification: dict[str, Any],
    history_path: Path,
) -> dict[str, Any]:
    if driver_audit.get("classification") != (
        "exploratory-pressure-only-fixed-state-diagnostic"
    ):
        raise SummaryError(f"{label} driver-screen classification is not current")
    if driver_audit.get("label") != label:
        raise SummaryError(f"{label} driver-screen label identity differs")
    if _positive_int(
        driver_audit.get("frame_count"), f"{label} driver frame_count"
    ) != EXPECTED_FRAME_COUNT:
        raise SummaryError(f"{label} driver screen must bind exactly 60 saved frames")
    driver_particle_count = _positive_int(
        driver_audit.get("particle_count"), f"{label} driver particle_count"
    )
    if driver_particle_count != phase_particle_count:
        raise SummaryError(f"{label} phase/driver particle counts differ")
    support_count = _positive_int(
        driver_audit.get("support_particle_count"), f"{label} support count"
    )
    if support_count > driver_particle_count:
        raise SummaryError(f"{label} support count exceeds particle count")

    screen_window = _mapping(
        driver_audit.get("screen_window"), f"{label} screen_window"
    )
    for key, expected in (
        ("start_s", EXPECTED_COMBINED_WINDOW_S[0]),
        ("end_s", EXPECTED_COMBINED_WINDOW_S[1]),
        ("period_s", EXPECTED_PERIOD_S),
    ):
        _same_float(
            _finite_float(screen_window.get(key), f"{label} screen_window.{key}"),
            expected,
            f"{label} screen_window.{key}",
        )

    figure = _mapping(driver_audit.get("figure"), f"{label} driver figure")
    if _explicit_bool(
        figure.get("display_smoothing_used_in_metrics"),
        f"{label} display_smoothing_used_in_metrics",
    ) is not False:
        raise SummaryError(f"{label} display smoothing contaminated driver metrics")

    rows = _read_raw_history(history_path, label)
    if len(rows) != int(driver_audit["frame_count"]):
        raise SummaryError(f"{label} driver frame count differs from raw history")
    windows = _validated_windows(driver_audit, rows, label)
    first_cycle, second_cycle, combined = windows

    gate = _mapping(driver_audit.get("diagnostic_gate"), f"{label} gate")
    gate_passed = _explicit_bool(gate.get("passed"), f"{label} gate.passed")
    checks = _mapping(gate.get("checks"), f"{label} gate.checks")
    if set(checks) != set(GATE_CHECKS):
        raise SummaryError(
            f"{label} gate checks are not the complete preregistered set"
        )
    recorded_check_values = {
        str(name): _explicit_bool(value, f"{label} gate check {name}")
        for name, value in checks.items()
    }
    _validate_phase_qualification_record(
        _mapping(driver_audit.get("phase_field"), f"{label} phase_field").get(
            "qualification"
        ),
        phase_qualification,
        f"{label} phase_field.qualification",
    )
    _validate_phase_qualification_record(
        gate.get("phase_qualification"),
        phase_qualification,
        f"{label} gate.phase_qualification",
    )
    minimum_resolved_area = _finite_float(
        gate.get("minimum_resolved_joint_advantage_area_m2"),
        f"{label} minimum resolved joint advantage area",
    )
    if minimum_resolved_area <= 0.0:
        raise SummaryError(f"{label} minimum resolved joint advantage area is invalid")
    recomputed_area_tolerance = max(
        AREA_COMPARISON_ABS_TOL_M2,
        64.0
        * np.finfo(np.float64).eps
        * max(1.0, minimum_resolved_area),
    )
    _same_float(
        _finite_float(
            gate.get("area_comparison_absolute_tolerance_m2"),
            f"{label} area comparison tolerance",
        ),
        recomputed_area_tolerance,
        f"{label} area comparison tolerance",
        relative_tolerance=0.0,
        absolute_tolerance=1.0e-18,
    )
    resolved_threshold = minimum_resolved_area - recomputed_area_tolerance
    selected_rows = [
        row
        for row in rows
        if EXPECTED_COMBINED_WINDOW_S[0] - 1.0e-12
        <= row["time_s"]
        <= EXPECTED_COMBINED_WINDOW_S[1] + 1.0e-12
    ]
    resolved_joint_advantage = [
        row["lagged_joint_area_m2"] >= resolved_threshold
        and row["lagged_joint_area_m2"] - row["phase_erased_joint_area_m2"]
        >= resolved_threshold
        for row in selected_rows
    ]
    consecutive_joint_advantage = any(
        left and right
        for left, right in zip(
            resolved_joint_advantage[:-1], resolved_joint_advantage[1:]
        )
    )
    net_fraction = combined["net_uplift_force_impulse"][
        "lagged_minus_phase_erased_fraction"
    ]
    joint_fraction = combined["shared_HD_joint_area_time"][
        "lagged_minus_phase_erased_fraction"
    ]
    recomputed_checks = {
        "crown_phase_lag_is_resolved": bool(phase_qualification["passed"]),
        "both_post_ramp_cycles_have_positive_signed_force_delta": all(
            window["signed_force_impulse"]["lagged_minus_phase_erased"] > 0.0
            for window in (first_cycle, second_cycle)
        ),
        "combined_signed_force_impulse_delta_positive": (
            combined["signed_force_impulse"]["lagged_minus_phase_erased"] > 0.0
        ),
        "both_post_ramp_cycles_have_positive_net_uplift_delta": all(
            window["net_uplift_force_impulse"]["lagged_minus_phase_erased"]
            > 0.0
            for window in (first_cycle, second_cycle)
        ),
        "combined_net_uplift_impulse_increase_ge_5pct": (
            net_fraction is not None
            and net_fraction >= MINIMUM_FORCE_IMPULSE_FRACTION
        ),
        "combined_IF_ge_1_area_time_direction_positive": (
            combined["IF_ge_1_area_time"]["lagged_minus_phase_erased"] > 0.0
        ),
        "combined_shared_HD_joint_area_time_nonzero": (
            combined["shared_HD_joint_area_time"]["lagged"] > 0.0
        ),
        "combined_shared_HD_joint_area_time_increase_ge_5pct": (
            joint_fraction is not None
            and joint_fraction >= MINIMUM_JOINT_AREA_TIME_FRACTION
        ),
        "resolved_joint_advantage_for_two_consecutive_saved_frames": (
            consecutive_joint_advantage
        ),
    }
    for name in GATE_CHECKS:
        if recorded_check_values[name] is not recomputed_checks[name]:
            raise SummaryError(
                f"{label} gate check {name} differs from raw-history recomputation"
            )
    recomputed_gate_passed = all(recomputed_checks.values())
    if gate_passed is not recomputed_gate_passed:
        raise SummaryError(f"{label} gate.passed disagrees with recomputed checks")
    net_uplift_gate_passed = all(
        recomputed_checks[name] for name in NET_GATE_CHECKS
    )
    expected_classification = (
        "pressure-only-screen-passed"
        if gate_passed
        else "pressure-only-screen-not-passed"
    )
    if gate.get("classification") != expected_classification:
        raise SummaryError(f"{label} gate classification disagrees with gate.passed")

    effects: dict[str, float | None] = {}
    for metric, output_name, _display_name, units in EFFECT_METRICS:
        if combined[metric]["units"] != units:
            raise SummaryError(f"{label} recomputed units unexpectedly drifted")
        fraction = combined[metric]["lagged_minus_phase_erased_fraction"]
        effects[output_name] = (
            None if fraction is None else 100.0 * float(fraction)
        )
    return {
        **effects,
        "recomputed_integration_windows": windows,
        "combined_raw_integrals": {
            metric: combined[metric]
            for metric, _raw_metric, _units in HISTORY_METRICS
        },
        "diagnostic_gate_passed": gate_passed,
        "diagnostic_gate_status": "PASS" if gate_passed else "FAIL",
        "diagnostic_gate_classification": expected_classification,
        "diagnostic_gate_checks": recomputed_checks,
        "net_uplift_gate_passed": net_uplift_gate_passed,
        "net_uplift_gate_status": (
            "PASS" if net_uplift_gate_passed else "FAIL"
        ),
        "local_positive_activity_can_override_net_gate": False,
    }


def load_case(root: Path, label: str) -> dict[str, Any]:
    """Load one exact current-runner/phase-v3/driver-screen audit set."""

    root = root.resolve()
    label = _safe_label(label)
    phase_path = _phase_audit_path(root, label)
    driver_path = _driver_audit_path(root, label)
    phase_audit = _read_audit(phase_path, PHASE_SCHEMA, f"{label} phase-v3 audit")
    case = phase_audit.get("case")
    if not isinstance(case, str) or not case or Path(case).name != case:
        raise SummaryError(f"{label} phase-v3 case identity is invalid")
    expected_figure_paths = tuple(
        phase_path.parent / f"{case}_{stem}{suffix}"
        for stem in (
            "liquid_pressure_2d_phase_lag",
            "liquid_pressure_time_depth",
        )
        for suffix in (".png", ".pdf")
    )
    try:
        phase_plot.validate_phase_audit_artifacts(
            phase_audit,
            expected_figure_paths=expected_figure_paths,
        )
    except (phase_plot.PhaseAuditError, OSError) as error:
        raise SummaryError(
            f"{label} phase-v3 artifacts are incomplete or stale"
        ) from error
    driver_audit = _read_audit(
        driver_path, DRIVER_SCHEMA, f"{label} driver-screen audit"
    )
    phase_parameters = _validated_parameters(phase_audit, f"{label} phase-v3")
    driver_parameters = _validated_parameters(driver_audit, f"{label} driver")
    parameters = _validate_parameter_pair(
        phase_parameters, driver_parameters, label
    )
    runner_source = _validated_runner_source(
        root, label, driver_audit, driver_parameters
    )
    runner_audit = runner_source.pop("_audit")
    _validate_phase_runner_input_binding(
        root, label, phase_audit, driver_audit, runner_audit
    )
    driver_provenance = _validated_driver_file_provenance(
        root, label, driver_audit, phase_path
    )
    phase_metrics = _validated_phase_metrics(phase_audit, label)
    driver_metrics = _validated_driver_metrics(
        driver_audit,
        label,
        int(phase_metrics["particle_count"]),
        _mapping(
            phase_metrics["phase_qualification"],
            f"{label} recomputed phase qualification",
        ),
        Path(driver_provenance["history_csv"]["path"]),
    )
    return {
        "label": label,
        **parameters,
        **phase_metrics,
        **driver_metrics,
        "runner_provenance_verified": True,
        "runner_audit_schema": runner_source["schema"],
        "runner_audit_stage": runner_source["stage"],
        "runner_audit_path": runner_source["path"],
        "runner_audit_sha256": runner_source["sha256"],
        "sources": {
            "runner_audit": runner_source,
            "phase_v3": {
                "path": str(phase_path.resolve()),
                "sha256": sha256(phase_path),
                "schema": PHASE_SCHEMA,
            },
            "driver_screen": {
                "path": str(driver_path.resolve()),
                "sha256": sha256(driver_path),
                "schema": DRIVER_SCHEMA,
            },
            "driver_provenance": driver_provenance,
        },
    }


def collect_cases(root: Path, labels: list[str] | tuple[str, ...]) -> list[dict[str, Any]]:
    """Return exactly four validated cases sorted by physical permeability."""

    if len(labels) != EXPECTED_CASE_COUNT:
        raise SummaryError(
            f"The permeability screen requires exactly {EXPECTED_CASE_COUNT} labels"
        )
    if len(set(labels)) != len(labels):
        raise SummaryError("Permeability-screen labels must be unique")
    observed_labels = set(labels)
    required_labels = set(DEFAULT_LABELS)
    if observed_labels != required_labels:
        missing = sorted(required_labels - observed_labels)
        unexpected = sorted(observed_labels - required_labels)
        raise SummaryError(
            "Permeability-screen labels must exactly match the registered "
            f"DEFAULT_LABELS; missing={missing}, unexpected={unexpected}"
        )
    cases = [load_case(root, label) for label in labels]
    cases.sort(key=lambda item: float(item["intrinsic_permeability_m2"]))
    permeabilities = [float(item["intrinsic_permeability_m2"]) for item in cases]
    if any(
        not left < right
        for left, right in zip(permeabilities[:-1], permeabilities[1:])
    ):
        raise SummaryError("The four permeability values must be unique")
    for observed, expected in zip(permeabilities, EXPECTED_PERMEABILITIES_M2):
        _same_float(
            observed,
            expected,
            "Four-point permeability grid",
            absolute_tolerance=0.0,
        )
    for case in cases:
        case["evidence_role"] = "primary/current-runner"
        case["primary_screen_point"] = True
        case["current_runner_point"] = True
    particle_counts = {int(item["particle_count"]) for item in cases}
    if len(particle_counts) != 1:
        raise SummaryError("Particle count drifts across permeability cases")
    return cases


def _csv_rows(cases: list[dict[str, Any]]) -> list[dict[str, Any]]:
    fields = (
        "label",
        "intrinsic_permeability_m2",
        "liquid_saturation",
        "gas_saturation",
        "pressure_smoothing",
        "pressure_smoothing_in_loop",
        "pipeline_fixed",
        "runner_provenance_verified",
        "runner_audit_schema",
        "runner_audit_stage",
        "runner_audit_path",
        "runner_audit_sha256",
        "crown_phase_difference_same_column_deg",
        "crown_phase_lag_magnitude_deg",
        "crown_fundamental_amplitude_pa",
        "surface_fundamental_amplitude_pa",
        "crown_to_surface_amplitude_percent",
        "crown_harmonic_r_squared",
        "phase_present",
        "combined_net_uplift_effect_percent",
        "combined_local_positive_activity_effect_percent",
        "combined_IF_ge_1_area_time_effect_percent",
        "combined_shared_HD_joint_area_time_effect_percent",
        "diagnostic_gate_passed",
        "diagnostic_gate_status",
        "net_uplift_gate_passed",
        "net_uplift_gate_status",
        "local_positive_activity_can_override_net_gate",
        "evidence_role",
        "primary_screen_point",
        "current_runner_point",
    )
    return [{field: case[field] for field in fields} for case in cases]


def plot_summary(
    cases: list[dict[str, Any]], png_path: Path, pdf_path: Path, *, dpi: int
) -> None:
    """Render validated values without smoothing, interpolation, or mutation."""

    if isinstance(dpi, bool) or not isinstance(dpi, int) or not 72 <= dpi <= 600:
        raise SummaryError("Figure DPI must be an integer in [72, 600]")
    permeability = np.asarray(
        [case["intrinsic_permeability_m2"] for case in cases], dtype=float
    )
    phase = np.asarray(
        [
            np.nan
            if case["crown_phase_lag_magnitude_deg"] is None
            else case["crown_phase_lag_magnitude_deg"]
            for case in cases
        ],
        dtype=float,
    )
    amplitude = np.asarray(
        [case["crown_to_surface_amplitude_percent"] for case in cases], dtype=float
    )
    r_squared = np.asarray(
        [case["crown_harmonic_r_squared"] for case in cases], dtype=float
    )
    if not all(
        np.all(np.isfinite(values))
        for values in (permeability, amplitude, r_squared)
    ) or np.any(np.isinf(phase)):
        raise SummaryError("Validated chart arrays unexpectedly contain NaN/Inf")

    blue = "#2D5F8B"
    orange = "#C66A2B"
    gold = "#B58A28"
    olive = "#6F7B3A"
    charcoal = "#27313A"
    quiet = "#D8DEE3"
    effect_colors = (blue, orange, gold, olive)

    figure, (phase_axis, effect_axis) = plt.subplots(
        2,
        1,
        figsize=(10.8, 8.2),
        sharex=True,
        gridspec_kw={"height_ratios": (1.0, 1.25)},
    )
    figure.subplots_adjust(
        left=0.10, right=0.90, top=0.85, bottom=0.18, hspace=0.22
    )
    figure.suptitle(
        "Four-point phase-lag permeability screen", fontsize=15, y=0.975
    )
    figure.text(
        0.5,
        0.925,
        (
            r"$S_w=0.94$; pressure smoothing=false (output and in-loop); "
            "fixed pipe; all points primary/current-runner; raw audit metrics only"
        ),
        ha="center",
        va="top",
        fontsize=9,
        color=charcoal,
    )

    phase_axis.set_xscale("log")
    phase_line = phase_axis.plot(
        permeability,
        phase,
        color=blue,
        marker="o",
        markerfacecolor="white",
        markeredgewidth=1.5,
        linewidth=1.8,
        label="Crown phase-lag magnitude",
    )[0]
    phase_axis.set_ylabel("|Crown phase lag| (deg)", color=blue)
    phase_axis.tick_params(axis="y", colors=blue)
    phase_axis.margins(y=0.16)
    phase_axis.grid(axis="y", color=quiet, linewidth=0.7)
    phase_axis.set_title(
        "Crown fundamental: same-column phase and transmitted amplitude",
        fontsize=11,
        loc="left",
    )
    amplitude_axis = phase_axis.twinx()
    amplitude_line = amplitude_axis.plot(
        permeability,
        amplitude,
        color=orange,
        marker="s",
        markerfacecolor="white",
        markeredgewidth=1.5,
        linestyle="--",
        linewidth=1.8,
        label="Crown/surface amplitude",
    )[0]
    amplitude_axis.set_ylabel("Crown / surface amplitude (%)", color=orange)
    amplitude_axis.tick_params(axis="y", colors=orange)
    amplitude_axis.margins(y=0.16)
    phase_axis.legend(
        (phase_line, amplitude_line),
        (phase_line.get_label(), amplitude_line.get_label()),
        loc="upper left",
        frameon=False,
        fontsize=8,
    )
    if not np.any(np.isfinite(phase)):
        phase_axis.set_ylim(0.0, 1.0)
    for x_value, y_value, r2_value in zip(permeability, phase, r_squared):
        if math.isfinite(float(y_value)):
            phase_axis.annotate(
                rf"$R^2={r2_value:.3f}$",
                (x_value, y_value),
                xytext=(0, 8),
                textcoords="offset points",
                ha="center",
                fontsize=7,
                color=charcoal,
            )
        else:
            phase_axis.text(
                x_value,
                0.06,
                "phase N/A",
                transform=phase_axis.get_xaxis_transform(),
                ha="center",
                va="bottom",
                fontsize=7,
                color=charcoal,
            )

    for (_, output_name, display_name, _units), color in zip(
        EFFECT_METRICS, effect_colors
    ):
        values = np.asarray(
            [
                np.nan if case[output_name] is None else case[output_name]
                for case in cases
            ],
            dtype=float,
        )
        if np.any(np.isinf(values)):
            raise SummaryError(f"{display_name} contains Inf")
        effect_axis.plot(
            permeability,
            values,
            color=color,
            linewidth=1.5,
            marker="o",
            markersize=5,
            label=display_name,
        )
        for x_value, value in zip(permeability, values):
            if not math.isfinite(float(value)):
                effect_axis.text(
                    x_value,
                    0.04 + 0.04 * effect_colors.index(color),
                    "N/A",
                    transform=effect_axis.get_xaxis_transform(),
                    ha="center",
                    va="bottom",
                    fontsize=6.5,
                    color=color,
                )
    effect_axis.axhline(0.0, color=charcoal, linewidth=0.9)
    effect_axis.grid(axis="y", color=quiet, linewidth=0.7)
    effect_axis.set_ylabel("Lagged relative to phase-erased (%)")
    effect_axis.set_xlabel(r"Intrinsic permeability, $k$ ($\mathrm{m}^2$; log scale)")
    effect_axis.set_title(
        "Combined two-cycle raw-integral effects (1.3–3.9 s)",
        fontsize=11,
        loc="left",
    )
    effect_axis.legend(ncol=2, loc="best", frameon=False, fontsize=8)

    effect_axis.set_xticks(permeability)
    effect_axis.set_xticklabels(
        [
            (
                f"{value:.0e}\nPRIMARY · CURRENT RUNNER\n"
                f"NET {case['net_uplift_gate_status']} · "
                f"SCREEN {case['diagnostic_gate_status']}"
            )
            for value, case in zip(permeability, cases)
        ],
        fontsize=8,
    )
    for tick, case in zip(effect_axis.get_xticklabels(), cases):
        tick.set_fontweight("bold" if case["net_uplift_gate_passed"] else "normal")
        if not case["net_uplift_gate_passed"]:
            tick.set_bbox(
                {"facecolor": "white", "edgecolor": charcoal, "pad": 2.0, "linewidth": 0.7}
            )

    figure.text(
        0.5,
        0.025,
        (
            "Points are discrete screen cases; lines only guide the eye. Local + is "
            "local positive-force activity and cannot rescue NET FAIL. Shared-HD "
            "joint is co-location, not a causal liquefaction result.\n"
            "All four points are primary/current-runner cases with rebound runner "
            "provenance; external repeatability references are excluded. No display "
            "smoothing or interpolation enters any metric."
        ),
        ha="center",
        va="bottom",
        fontsize=7.3,
        color=charcoal,
        wrap=True,
    )
    figure.savefig(png_path, dpi=dpi, bbox_inches="tight", facecolor="white")
    figure.savefig(pdf_path, bbox_inches="tight", facecolor="white")
    plt.close(figure)


def _commit_output_transaction(
    staged_to_final: list[tuple[Path, Path]], output_dir: Path
) -> None:
    """Install a complete artifact set, restoring the old set on any error."""

    output_dir.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(
        prefix=".phase-lag-summary-backup-", dir=output_dir.parent
    ) as backup_name:
        backup_dir = Path(backup_name)
        backups: list[tuple[Path, Path]] = []
        installed: list[Path] = []
        try:
            for _staged, final in staged_to_final:
                if final.exists():
                    if not final.is_file():
                        raise SummaryError(
                            f"Summary target is not a regular file: {final}"
                        )
                    backup = backup_dir / final.name
                    os.replace(final, backup)
                    backups.append((backup, final))
            for staged, final in staged_to_final:
                os.replace(staged, final)
                installed.append(final)
        except Exception as error:
            rollback_errors: list[Exception] = []
            for final in reversed(installed):
                try:
                    if final.is_file():
                        final.unlink()
                except OSError as rollback_error:
                    rollback_errors.append(rollback_error)
            for backup, final in reversed(backups):
                try:
                    if backup.is_file():
                        os.replace(backup, final)
                except OSError as rollback_error:
                    rollback_errors.append(rollback_error)
            if rollback_errors:
                raise SummaryError(
                    "Summary transaction failed and rollback was incomplete"
                ) from error
            if isinstance(error, SummaryError):
                raise
            raise SummaryError("Summary output transaction failed") from error


def write_summary(
    root: Path,
    labels: list[str] | tuple[str, ...],
    output_dir: Path,
    *,
    dpi: int = 300,
) -> dict[str, Any]:
    """Validate all inputs first, then write the four summary artifacts."""

    root = root.resolve()
    cases = collect_cases(root, labels)
    output_dir = output_dir.resolve()
    results_dir = (root / "results").resolve()
    if output_dir == results_dir or results_dir in output_dir.parents:
        raise SummaryError("Summary output must not be written under results/")
    base = output_dir / "phase_lag_permeability_screen"
    csv_path = base.with_suffix(".csv")
    png_path = base.with_suffix(".png")
    pdf_path = base.with_suffix(".pdf")
    audit_path = output_dir / "phase_lag_permeability_screen.audit.json"

    output_dir.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(
        prefix=".phase-lag-summary-stage-", dir=output_dir.parent
    ) as staging_name:
        staging = Path(staging_name)
        staged_csv = staging / csv_path.name
        staged_png = staging / png_path.name
        staged_pdf = staging / pdf_path.name
        staged_audit = staging / audit_path.name
        csv_rows = _csv_rows(cases)
        with staged_csv.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(csv_rows[0]))
            writer.writeheader()
            writer.writerows(csv_rows)
        plot_summary(cases, staged_png, staged_pdf, dpi=dpi)

        audit = {
            "schema": SCHEMA,
            "classification": "exploratory-four-point-fixed-state-screen-summary",
            "labels_requested": list(labels),
            "labels_sorted_by_permeability": [case["label"] for case in cases],
            "case_count": len(cases),
            "passed_case_count": sum(
                bool(case["diagnostic_gate_passed"]) for case in cases
            ),
            "net_uplift_passed_case_count": sum(
                bool(case["net_uplift_gate_passed"]) for case in cases
            ),
            "required_parameters": {
                "liquid_saturation": EXPECTED_SATURATION,
                "pressure_smoothing": False,
                "pressure_smoothing_in_loop": False,
                "pipeline_fixed": True,
            },
            "combined_window_s": list(EXPECTED_COMBINED_WINDOW_S),
            "effect_definition": (
                "100 * (lagged raw integral / phase-erased raw integral - 1); "
                "null/N/A when the phase-erased denominator is non-positive; "
                "every integral and gate is recomputed from the hashed raw history"
            ),
            "evidence_policy": {
                "all_four_points_primary": True,
                "primary_role": "primary/current-runner",
                "current_runner_required": True,
                "required_runner_schema": RUNNER_SCHEMA,
                "required_runner_stage": RUNNER_STAGE,
                "external_repeatability_references_included": False,
            },
            "pending_policy": (
                "fail-closed: all four current runner, phase-v3, and driver-screen "
                "audits, their complete path/hash chains, raw histories, registered "
                "windows, and preregistered gates must validate before a transactional "
                "replacement; missing scores remain N/A and cannot pass a gate"
            ),
            "gate_policy": (
                "All nine preregistered checks are independently recomputed. "
                "Net-uplift gate status is reported independently; positive local "
                "force activity cannot override a net-uplift FAIL."
            ),
            "display_contract": {
                "renderer": "matplotlib-static",
                "dpi": dpi,
                "x_scale": "logarithmic permeability",
                "smoothing": False,
                "interpolation": False,
                "display_changes_numerical_values": False,
                "status_encoding": "PASS/FAIL/N/A text; FAIL has an outlined tick label",
            },
            "limitations": (
                "The driver values are pressure-only fixed-state diagnostics. Local "
                "positive-force activity is not net pipeline uplift, and shared-HD "
                "joint area-time is co-location rather than a causal liquefaction "
                "difference or released-pipeline response."
            ),
            "cases": cases,
            "outputs": {
                "csv": {"path": str(csv_path), "sha256": sha256(staged_csv)},
                "png": {"path": str(png_path), "sha256": sha256(staged_png)},
                "pdf": {"path": str(pdf_path), "sha256": sha256(staged_pdf)},
            },
            "audit_path": str(audit_path),
        }
        staged_audit.write_text(
            json.dumps(audit, indent=2, allow_nan=False) + "\n", encoding="utf-8"
        )
        _commit_output_transaction(
            [
                (staged_csv, csv_path),
                (staged_png, png_path),
                (staged_pdf, pdf_path),
                (staged_audit, audit_path),
            ],
            output_dir,
        )
    return audit


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--labels",
        nargs="+",
        default=list(DEFAULT_LABELS),
        help="Exactly four exploratory labels (default: the declared k screen)",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Summary directory (default: analysis/phase_lag_exploratory/permeability_screen)",
    )
    parser.add_argument("--dpi", type=int, default=300)
    args = parser.parse_args()
    root = Path(__file__).resolve().parent
    output_dir = args.output_dir or (
        root / "analysis/phase_lag_exploratory/permeability_screen"
    )
    audit = write_summary(root, args.labels, output_dir, dpi=args.dpi)
    print(
        json.dumps(
            {
                "labels_sorted_by_permeability": audit[
                    "labels_sorted_by_permeability"
                ],
                "passed_case_count": audit["passed_case_count"],
                "audit_path": audit["audit_path"],
                "outputs": audit["outputs"],
            },
            indent=2,
            allow_nan=False,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
