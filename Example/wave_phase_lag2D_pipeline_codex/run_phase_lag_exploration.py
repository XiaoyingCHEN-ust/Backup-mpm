#!/usr/bin/env python3
"""Fail-closed, stage-by-stage runner for exploratory phase-lag probes.

The registered study validator intentionally rejects exploratory UUIDs.  This
runner therefore supplies the missing local safety envelope without weakening
that validator: a label can be prepared only once, every solver stage requires
an empty target and an unchanged provenance record, and a failed/partial stage
is never restarted in place.  Only the equilibrium and fully coupled pressure
driver are runnable here; pressure replay cases remain explicit later steps.
"""

from __future__ import annotations

import argparse
import contextlib
import fcntl
import hashlib
import json
import math
import os
import re
import subprocess
from pathlib import Path
from typing import Any, Iterator

import numpy as np

import hdf5_checkpoint_crosscheck
import prepare_phase_lag_exploration as generator
import qa_contract as qa
import validate_case
from analyze_study import (
    expected_particle_steps,
    particle_step,
    read_initial_coordinates,
    read_pipeline_history,
    read_vtp,
    validate_pipeline_history_grid,
)
from phase_controls import (
    database_paths,
    iter_frames,
    read_header,
    read_points,
    validate_database,
)


SCHEMA = "pipeline-phase-lag-exploration-runner-v1"
WAVE_STEPS = 39_000
WAVE_OUTPUT_INTERVAL = 650
PRESSURE_INTERVAL = 130
DT = 1.0e-4
EXPECTED_WAVE_FRAMES = 60
CONFIG_NAMES = ("01_EQ.json", "02_WAVE.json", "03_HD.json", "04_RL.json", "05_RE.json")
HD_REQUIRED_ARRAYS = (
    "ids",
    "volumes",
    "velocities",
    "displacements",
    "porosities",
    "stresses",
    "initial_vertical_effective_stresses",
    "vertical_effective_stress_remaining_ratios",
    "PIC_pore_pressure_excess",
    "liquid_seepage_forces",
    "liquid_densities",
    "gamma_sub",
)


class SafetyError(RuntimeError):
    """Raised when continuing could trust or overwrite unaudited state."""


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def artifact(path: Path) -> dict[str, Any]:
    path = path.resolve()
    if not path.is_file() or path.stat().st_size <= 0:
        raise SafetyError(f"Required artifact is missing or empty: {path}")
    return {
        "path": str(path),
        "size_bytes": path.stat().st_size,
        "sha256": sha256(path),
    }


def label_token(label: str) -> str:
    if re.fullmatch(r"[a-z0-9][a-z0-9_-]*", label) is None:
        raise SafetyError("label must match [a-z0-9][a-z0-9_-]*")
    return label.upper().replace("-", "_")


def state_path(root: Path, label: str) -> Path:
    return root / "analysis/phase_lag_exploratory" / label / "runner_audit.json"


def atomic_write_json(path: Path, value: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    try:
        with temporary.open("x", encoding="utf-8", newline="\n") as stream:
            json.dump(value, stream, indent=2, sort_keys=True, allow_nan=False)
            stream.write("\n")
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
    finally:
        if temporary.exists():
            temporary.unlink()


@contextlib.contextmanager
def label_lock(root: Path, label: str) -> Iterator[None]:
    token = label_token(label)
    lock_root = root / "results/.phase_lag_exploration_locks"
    lock_root.mkdir(parents=True, exist_ok=True)
    lock_path = lock_root / f"{token}.lock"
    with lock_path.open("a+", encoding="utf-8") as stream:
        try:
            fcntl.flock(stream.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as error:
            raise SafetyError(f"Exploration label is already locked: {label}") from error
        try:
            yield
        finally:
            fcntl.flock(stream.fileno(), fcntl.LOCK_UN)


def namespace_paths(root: Path, label: str) -> tuple[Path, ...]:
    return (
        root / "configs/phase_lag_exploratory" / label,
        root / "results/phase_lag_exploratory" / label,
        root / "pressure_databases/phase_lag_exploratory" / label,
        root / "analysis/phase_lag_exploratory" / label,
    )


def require_new_label(root: Path, label: str) -> None:
    label_token(label)
    existing = [path for path in namespace_paths(root, label) if path.exists() or path.is_symlink()]
    log_root = root / "logs/phase_lag_exploratory"
    existing.extend(sorted(log_root.glob(f"{label}_*.log")) if log_root.is_dir() else [])
    if existing:
        raise SafetyError(
            "Exploration labels are single-use; refusing existing namespace(s): "
            + ", ".join(str(path) for path in existing)
        )


def require_empty_target(path: Path, description: str) -> None:
    if path.is_symlink():
        raise SafetyError(f"Refusing symlinked {description}: {path}")
    if path.exists() and not path.is_dir():
        raise SafetyError(f"Expected a directory for {description}: {path}")
    if path.is_dir() and next(path.iterdir(), None) is not None:
        raise SafetyError(f"Refusing non-empty {description}: {path}")


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _same_number(actual: Any, expected: float) -> bool:
    return math.isclose(float(actual), expected, rel_tol=0.0, abs_tol=1.0e-15)


def validate_generated_family(
    root: Path,
    label: str,
    permeability_m2: float,
    saturation: float,
    equilibrium_steps: int,
) -> dict[str, Any]:
    directory = root / "configs/phase_lag_exploratory" / label
    manifest_path = directory / "manifest.json"
    manifest = read_json(manifest_path)
    if manifest.get("schema") != generator.SCHEMA or manifest.get("label") != label:
        raise SafetyError("Generated exploration manifest identity is invalid")
    parameters = manifest.get("parameters", {})
    expected_parameters = {
        "intrinsic_permeability_m2": permeability_m2,
        "liquid_saturation": saturation,
        "gas_saturation": 1.0 - saturation,
        "equilibrium_steps": equilibrium_steps,
        "wave_steps": WAVE_STEPS,
        "dt_s": DT,
        "wave_output_interval_steps": WAVE_OUTPUT_INTERVAL,
        "pressure_smoothing": False,
        "paired_replay": True,
        "pressure_database_interval_steps": PRESSURE_INTERVAL,
    }
    if parameters != expected_parameters:
        raise SafetyError(f"Generated manifest parameters changed: {parameters}")

    baseline = manifest.get("baseline", {})
    for prefix in ("equilibrium", "wave"):
        source = root / str(baseline.get(f"{prefix}_config", ""))
        if not source.is_file() or sha256(source) != baseline.get(f"{prefix}_sha256"):
            raise SafetyError(f"Generated manifest baseline hash is stale: {prefix}")

    generated = manifest.get("generated", {})
    key_by_name = {
        "01_EQ.json": "equilibrium",
        "02_WAVE.json": "wave",
        "03_HD.json": "hd",
        "04_RL.json": "rl",
        "05_RE.json": "re",
    }
    configs: dict[str, dict[str, Any]] = {}
    for name, key in key_by_name.items():
        path = directory / name
        recorded = root / str(generated.get(f"{key}_config", ""))
        if path.resolve() != recorded.resolve() or sha256(path) != generated.get(f"{key}_sha256"):
            raise SafetyError(f"Generated config hash/path mismatch: {name}")
        configs[name] = read_json(path)
    stress_path = directory / "initial_effective_stresses.txt"
    if (
        stress_path.resolve()
        != (root / str(generated.get("initial_effective_stresses", ""))).resolve()
        or sha256(stress_path) != generated.get("initial_effective_stresses_sha256")
    ):
        raise SafetyError("Generated initial-stress hash/path mismatch")

    token = label_token(label)
    expected_result = f"results/phase_lag_exploratory/{label}/"
    expected_eq_uuid = f"PLP_EXP_{token}_EQ"
    for name, config in configs.items():
        analysis = config["analysis"]
        soil, fluid = config["materials"]
        if not _same_number(soil["intrinsic_permeability"], permeability_m2):
            raise SafetyError(f"{name}: permeability changed")
        if not _same_number(fluid["liquid_saturation"], saturation) or not _same_number(
            fluid["gas_saturation"], 1.0 - saturation
        ):
            raise SafetyError(f"{name}: saturation changed")
        if (
            analysis.get("pressure_smoothing") is not False
            or analysis.get("pressure_smoothing_in_loop") is not False
        ):
            raise SafetyError(
                f"{name}: pressure smoothing must be explicitly disabled"
            )
        if config["post_processing"]["path"] != expected_result:
            raise SafetyError(f"{name}: result namespace changed")
        if not bool(analysis["rigid_pipeline"]["fixed"]):
            raise SafetyError(f"{name}: pipeline must remain fixed")

    equilibrium = configs["01_EQ.json"]
    eq_analysis = equilibrium["analysis"]
    if (
        equilibrium["materials"][0]["type"] != "LinearElastic2D"
        or bool(eq_analysis["APIC"])
        or not _same_number(eq_analysis["PIC"], 1.0)
        or equilibrium["materials"][1]["wave_pressure"]
        or int(eq_analysis["nsteps"]) != equilibrium_steps
        or eq_analysis["uuid"] != expected_eq_uuid
        or bool(eq_analysis["resume"]["resume"])
        or not _same_number(eq_analysis["damping"]["damping_factor"], 5.0)
        or eq_analysis.get("stability_qa_contract")
        != qa.stability_contract(qa.LINEAR_EQUILIBRIUM_MODE)
    ):
        raise SafetyError("Generated equilibrium protocol changed")

    role_modes = {
        "03_HD.json": ("HD", True, False, True, "lagged", "pressure"),
        "04_RL.json": ("RL", False, True, False, "lagged", "pressure"),
        "05_RE.json": ("RE", False, True, False, "phase_erased", "phase_erased"),
    }
    for name, (role, wave_pressure, enable, write, leaf, prefix) in role_modes.items():
        config = configs[name]
        analysis = config["analysis"]
        pressure = analysis.get("prescribed_phase_pressures", {})
        expected_path = f"pressure_databases/phase_lag_exploratory/{label}/{leaf}"
        resume = analysis["resume"]
        if (
            config["materials"][0]["type"] != "SANISAND2D"
            or not bool(analysis["APIC"])
            or not _same_number(analysis["PIC"], 0.0)
            or not _same_number(analysis["damping"]["damping_factor"], 0.0)
            or int(analysis["nsteps"]) != WAVE_STEPS
            or analysis["uuid"] != f"PLP_EXP_{token}_{role}"
            or bool(config["materials"][1]["wave_pressure"]) != wave_pressure
            or not bool(resume["resume"])
            or resume["uuid"] != expected_eq_uuid
            or int(resume["step"]) != equilibrium_steps
            or int(resume["nsteps"]) != equilibrium_steps
            or pressure.get("enable") is not enable
            or pressure.get("write") is not write
            or pressure.get("path") != expected_path
            or pressure.get("file_prefix") != prefix
            or int(pressure.get("step_interval", -1)) != PRESSURE_INTERVAL
            or int(pressure.get("max_step", -1)) != WAVE_STEPS
            or not _same_number(pressure.get("source_dt", -1.0), DT)
            or pressure.get("mapping") != "id"
        ):
            raise SafetyError(f"Generated {role} protocol changed")
    return manifest


def prepare_stage(
    root: Path,
    *,
    label: str,
    permeability_m2: float,
    saturation: float,
    equilibrium_steps: int,
    mpm_binary: Path,
) -> dict[str, Any]:
    root = root.resolve()
    mpm_binary = mpm_binary.resolve()
    with label_lock(root, label):
        require_new_label(root, label)
        if not mpm_binary.is_file() or not os.access(mpm_binary, os.X_OK):
            raise SafetyError(f"MPM binary is missing or not executable: {mpm_binary}")
        generator.generate(
            root,
            label=label,
            permeability_m2=permeability_m2,
            saturation=saturation,
            equilibrium_steps=equilibrium_steps,
            wave_steps=WAVE_STEPS,
            pressure_smoothing=False,
            paired_replay=True,
        )
        validate_generated_family(
            root, label, permeability_m2, saturation, equilibrium_steps
        )
        config_dir = root / "configs/phase_lag_exploratory" / label
        state = {
            "schema": SCHEMA,
            "label": label,
            "stage": "prepared",
            "parameters": {
                "intrinsic_permeability_m2": permeability_m2,
                "liquid_saturation": saturation,
                "pressure_smoothing": False,
                "equilibrium_steps": equilibrium_steps,
                "wave_steps": WAVE_STEPS,
                "dt_s": DT,
            },
            "solver_binary": artifact(mpm_binary),
            "manifest": artifact(config_dir / "manifest.json"),
            "configs": {
                name: artifact(config_dir / name) for name in CONFIG_NAMES
            },
            "initial_effective_stresses": artifact(
                config_dir / "initial_effective_stresses.txt"
            ),
            "stages": {},
        }
        atomic_write_json(state_path(root, label), state)
        return state


def load_state(root: Path, label: str) -> dict[str, Any]:
    path = state_path(root, label)
    if not path.is_file():
        raise SafetyError(f"Runner audit is missing; prepare the label first: {path}")
    state = read_json(path)
    if state.get("schema") != SCHEMA or state.get("label") != label:
        raise SafetyError("Runner audit identity/schema changed")
    return state


def verify_artifact(record: dict[str, Any], description: str) -> Path:
    path = Path(str(record.get("path", "")))
    if (
        not path.is_absolute()
        or not path.is_file()
        or path.stat().st_size != int(record.get("size_bytes", -1))
        or sha256(path) != record.get("sha256")
    ):
        raise SafetyError(f"Recorded {description} changed: {path}")
    return path


def verify_provenance(root: Path, label: str, state: dict[str, Any]) -> None:
    verify_artifact(state["solver_binary"], "solver binary")
    verify_artifact(state["manifest"], "manifest")
    verify_artifact(state["initial_effective_stresses"], "initial stresses")
    for name, record in state["configs"].items():
        expected = root / "configs/phase_lag_exploratory" / label / name
        if verify_artifact(record, f"config {name}").resolve() != expected.resolve():
            raise SafetyError(f"Config escaped the prepared label: {name}")
    parameters = state["parameters"]
    validate_generated_family(
        root,
        label,
        float(parameters["intrinsic_permeability_m2"]),
        float(parameters["liquid_saturation"]),
        int(parameters["equilibrium_steps"]),
    )


def result_directory(root: Path, config: dict[str, Any]) -> Path:
    return root / config["post_processing"]["path"] / config["analysis"]["uuid"]


def exact_history(result: Path, config: dict[str, Any]) -> tuple[Path, dict[str, Any]]:
    histories = sorted(result.glob("pipeline-history*.csv"))
    if len(histories) != 1:
        raise SafetyError(f"Expected exactly one pipeline history in {result}")
    audit = validate_pipeline_history_grid(read_pipeline_history(result), config)
    return histories[0], audit


def audit_equilibrium(root: Path, config_path: Path) -> dict[str, Any]:
    config = read_json(config_path)
    nsteps = int(config["analysis"]["nsteps"])
    result = result_directory(root, config)
    steps = sorted((particle_step(path), path) for path in result.glob("particle*.vtp"))
    expected = expected_particle_steps(config)
    if [step for step, _ in steps] != expected:
        raise SafetyError(f"Equilibrium VTP grid differs: {[step for step, _ in steps]}")
    digits = len(str(nsteps))
    vtp = result / f"particle{nsteps:0{digits}d}.vtp"
    hdf5 = result / f"particles{nsteps:0{digits}d}.h5"
    validate_case.validate_final_vtp(vtp, config)
    validate_case.validate_final_hdf5(hdf5)
    stability = validate_case.validate_equilibrium_vtp(
        hdf5, config, mode=qa.LINEAR_EQUILIBRIUM_MODE
    )
    maximum_velocity = float(stability["observed"]["maximum_velocity_m_s"])
    if maximum_velocity > qa.MAXIMUM_VELOCITY_M_S:
        raise SafetyError(
            f"Equilibrium max|v|={maximum_velocity:.6g} exceeds unchanged 1e-3 m/s gate"
        )
    crosscheck = hdf5_checkpoint_crosscheck.crosscheck(hdf5, vtp)
    history, history_grid = exact_history(result, config)
    return {
        "result_directory": str(result.resolve()),
        "particle_steps": expected,
        "final_vtp": artifact(vtp),
        "final_hdf5": artifact(hdf5),
        "pipeline_history": artifact(history),
        "pipeline_history_grid": history_grid,
        "stability_qa": stability,
        "hdf5_vtp_crosscheck": crosscheck,
    }


def mesh_bounds(path: Path) -> tuple[np.ndarray, np.ndarray]:
    with path.open(encoding="utf-8") as stream:
        line = next(stream)
        while line.lstrip().startswith("#"):
            line = next(stream)
        node_count = int(line.split()[0])
        nodes = np.asarray(
            [list(map(float, next(stream).split()[:2])) for _ in range(node_count)],
            dtype=np.float64,
        )
    if nodes.shape != (node_count, 2) or not np.all(np.isfinite(nodes)):
        raise SafetyError(f"Invalid mesh-node table: {path}")
    return np.min(nodes, axis=0), np.max(nodes, axis=0)


def ordered_vtp_frame(path: Path, count: int) -> tuple[np.ndarray, dict[str, np.ndarray]]:
    points, arrays = read_vtp(path)
    missing = [name for name in HD_REQUIRED_ARRAYS if name not in arrays]
    if missing:
        raise SafetyError(f"{path}: missing HD arrays {missing}")
    raw_ids = np.asarray(arrays["ids"], dtype=np.float64).reshape(-1)
    ids = np.rint(raw_ids).astype(np.int64)
    if not np.all(np.isfinite(raw_ids)) or not np.array_equal(raw_ids, ids.astype(float)):
        raise SafetyError(f"{path}: particle IDs are non-finite/non-integral")
    order = np.argsort(ids)
    if raw_ids.size != count or not np.array_equal(ids[order], np.arange(count)):
        raise SafetyError(f"{path}: particle IDs are not unique 0..N-1")
    ordered = {name: np.asarray(values)[order] for name, values in arrays.items()}
    points = np.asarray(points, dtype=np.float64)[order, :2]
    if not np.all(np.isfinite(points)):
        raise SafetyError(f"{path}: particle coordinates contain NaN/Inf")
    for name in HD_REQUIRED_ARRAYS[1:]:
        if not np.all(np.isfinite(ordered[name])):
            raise SafetyError(f"{path}: {name} contains NaN/Inf")
    return points, ordered


def audit_hd(root: Path, config_path: Path) -> dict[str, Any]:
    config = read_json(config_path)
    analysis = config["analysis"]
    result = result_directory(root, config)
    files = sorted(result.glob("particle*.vtp"), key=particle_step)
    expected_steps = expected_particle_steps(config)
    actual_steps = [particle_step(path) for path in files]
    if actual_steps != expected_steps or len(files) != EXPECTED_WAVE_FRAMES:
        raise SafetyError(
            f"HD requires exact {EXPECTED_WAVE_FRAMES}-frame grid; got {actual_steps}"
        )
    coordinates = read_initial_coordinates(
        root / config["particles"][0]["generator"]["location"]
    )
    count = len(coordinates)
    lower, upper = mesh_bounds(root / config["mesh"]["mesh"])
    outside_max = 0
    velocity_max = 0.0
    displacement_max = 0.0
    porosity_min = math.inf
    porosity_max = -math.inf
    for path in files:
        points, arrays = ordered_vtp_frame(path, count)
        outside = np.any((points < lower - 1.0e-8) | (points > upper + 1.0e-8), axis=1)
        outside_count = int(np.count_nonzero(outside))
        outside_max = max(outside_max, outside_count)
        if outside_count:
            raise SafetyError(f"{path}: {outside_count} particles are outside the mesh")
        volumes = np.asarray(arrays["volumes"], dtype=float).reshape(-1)
        porosity = np.asarray(arrays["porosities"], dtype=float).reshape(-1)
        if np.any(volumes <= 0.0):
            raise SafetyError(f"{path}: non-positive particle volume")
        if np.any(porosity <= 0.0) or np.any(porosity >= 1.0):
            raise SafetyError(f"{path}: porosity lies outside (0,1)")
        velocity = np.asarray(arrays["velocities"], dtype=float)[:, :2]
        displacement = np.asarray(arrays["displacements"], dtype=float)[:, :2]
        velocity_max = max(velocity_max, float(np.max(np.linalg.norm(velocity, axis=1))))
        displacement_max = max(
            displacement_max, float(np.max(np.linalg.norm(displacement, axis=1)))
        )
        porosity_min = min(porosity_min, float(np.min(porosity)))
        porosity_max = max(porosity_max, float(np.max(porosity)))

    history, history_grid = exact_history(result, config)
    pressure = analysis.get("prescribed_phase_pressures", {})
    if bool(pressure.get("enable")) or not bool(pressure.get("write")):
        raise SafetyError("HD pressure database is not in write-only mode")
    directory = root / str(pressure["path"])
    points_path, values_path = database_paths(directory, str(pressure["file_prefix"]))
    with values_path.open("rb") as stream:
        header = read_header(stream)
        frame_count = sum(1 for _ in iter_frames(stream, header))
    validate_database(values_path, header)
    pressure_points = read_points(points_path, header.dimension)
    ordered_ids = np.sort(pressure_points.particle_ids.astype(np.int64))
    if (
        header.format_version != "V2"
        or header.dimension != 2
        or header.particle_count != count
        or header.step_interval != PRESSURE_INTERVAL
        or header.max_step != WAVE_STEPS
        or not _same_number(header.source_dt, DT)
        or frame_count != WAVE_STEPS // PRESSURE_INTERVAL + 1
        or not np.array_equal(ordered_ids, np.arange(count))
        or len(pressure_points.particle_ids) != count
        or not np.all(np.isfinite(pressure_points.coordinates))
    ):
        raise SafetyError("HD pressure database header/IDs/frame grid is invalid")

    return {
        "result_directory": str(result.resolve()),
        "particle_steps": actual_steps,
        "particle_vtp": [artifact(path) for path in files],
        "pipeline_history": artifact(history),
        "pipeline_history_grid": history_grid,
        "pressure_database": {
            "points": artifact(points_path),
            "values": artifact(values_path),
            "header": {
                "format_version": header.format_version,
                "dimension": header.dimension,
                "particle_count": header.particle_count,
                "step_interval": header.step_interval,
                "max_step": header.max_step,
                "source_dt_s": header.source_dt,
                "frame_count": frame_count,
            },
        },
        "health": {
            "particle_count": count,
            "saved_frames": len(files),
            "outside_mesh_maximum_count": outside_max,
            "maximum_velocity_m_s": velocity_max,
            "maximum_displacement_m": displacement_max,
            "minimum_porosity": porosity_min,
            "maximum_porosity": porosity_max,
        },
    }


def run_solver(
    root: Path,
    *,
    label: str,
    stage: str,
    config_path: Path,
    binary: Path,
    threads: int,
) -> tuple[list[str], Path]:
    if threads <= 0:
        raise SafetyError("threads must be a positive integer")
    config = read_json(config_path)
    target = result_directory(root, config)
    require_empty_target(target, f"{stage} result target")
    if stage == "hd":
        pressure = config["analysis"]["prescribed_phase_pressures"]
        require_empty_target(root / pressure["path"], "HD pressure-database target")
    log_root = root / "logs/phase_lag_exploratory"
    log_root.mkdir(parents=True, exist_ok=True)
    log_path = log_root / f"{label}_{stage.upper()}.log"
    if log_path.exists() or log_path.is_symlink():
        raise SafetyError(f"Refusing existing stage log: {log_path}")
    relative_config = config_path.resolve().relative_to(root.resolve())
    command = [str(binary), "-p", str(threads), "-f", "./", "-i", str(relative_config)]
    environment = os.environ.copy()
    environment["OMP_NUM_THREADS"] = str(threads)
    environment["TBB_NUM_THREADS"] = str(threads)
    with log_path.open("x", encoding="utf-8") as log:
        completed = subprocess.run(
            command,
            cwd=root,
            env=environment,
            stdout=log,
            stderr=subprocess.STDOUT,
            check=False,
        )
    if completed.returncode != 0:
        raise SafetyError(
            f"{stage.upper()} solver failed with code {completed.returncode}; "
            f"partial state is quarantined in place: {log_path}"
        )
    return command, log_path


def solver_stage(
    root: Path,
    *,
    label: str,
    stage: str,
    threads: int,
) -> dict[str, Any]:
    if stage not in {"eq", "hd"}:
        raise SafetyError(f"Unsupported solver stage: {stage}")
    root = root.resolve()
    with label_lock(root, label):
        state = load_state(root, label)
        expected_state = "prepared" if stage == "eq" else "eq_complete"
        if state.get("stage") != expected_state:
            raise SafetyError(
                f"{stage.upper()} requires runner stage {expected_state}; got {state.get('stage')}"
            )
        verify_provenance(root, label, state)
        if stage == "hd":
            for name in ("final_vtp", "final_hdf5", "pipeline_history"):
                verify_artifact(state["stages"]["eq"]["audit"][name], f"EQ {name}")
            maximum_velocity = float(
                state["stages"]["eq"]["audit"]["stability_qa"]["observed"][
                    "maximum_velocity_m_s"
                ]
            )
            if maximum_velocity > qa.MAXIMUM_VELOCITY_M_S:
                raise SafetyError("Recorded EQ stability gate is no longer valid")
        filename = "01_EQ.json" if stage == "eq" else "03_HD.json"
        config_path = root / "configs/phase_lag_exploratory" / label / filename
        binary = verify_artifact(state["solver_binary"], "solver binary")
        command, log_path = run_solver(
            root,
            label=label,
            stage=stage,
            config_path=config_path,
            binary=binary,
            threads=threads,
        )
        audit = (
            audit_equilibrium(root, config_path)
            if stage == "eq"
            else audit_hd(root, config_path)
        )
        state["stages"][stage] = {
            "command": command,
            "threads": threads,
            "log": artifact(log_path),
            "config": artifact(config_path),
            "solver_binary": artifact(binary),
            "audit": audit,
        }
        state["stage"] = f"{stage}_complete"
        atomic_write_json(state_path(root, label), state)
        return state


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="stage", required=True)
    prepare = subparsers.add_parser("prepare", help="create one new unsmoothed paired family")
    prepare.add_argument("--label", required=True)
    prepare.add_argument("--permeability", required=True, type=float)
    prepare.add_argument("--saturation", required=True, type=float)
    prepare.add_argument("--equilibrium-steps", type=int, default=generator.DEFAULT_EQ_STEPS)
    prepare.add_argument("--mpm-bin", type=Path)
    for name in ("eq", "hd"):
        stage_parser = subparsers.add_parser(name, help=f"run and audit only the {name.upper()} stage")
        stage_parser.add_argument("--label", required=True)
        stage_parser.add_argument("--threads", type=int, default=3)
    status = subparsers.add_parser("status", help="print the immutable stage audit")
    status.add_argument("--label", required=True)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    root = Path(__file__).resolve().parent
    if args.stage == "prepare":
        default_binary = root.parents[1] / "mpm_hpc_source/build-pipeline/mpm"
        state = prepare_stage(
            root,
            label=args.label,
            permeability_m2=args.permeability,
            saturation=args.saturation,
            equilibrium_steps=args.equilibrium_steps,
            mpm_binary=args.mpm_bin or default_binary,
        )
    elif args.stage in {"eq", "hd"}:
        state = solver_stage(
            root, label=args.label, stage=args.stage, threads=args.threads
        )
    else:
        state = load_state(root, args.label)
    print(json.dumps(state, indent=2, sort_keys=True, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
