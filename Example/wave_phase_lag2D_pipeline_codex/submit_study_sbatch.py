#!/usr/bin/env python3
"""Submit a restart-aware pipeline study DAG to HPC4 with ``sbatch``."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import re
import shlex
import struct
import subprocess
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable


CASE_DIR = Path(__file__).resolve().parent
RUN_TASK = "run_task.sbatch"
COMPLETION_SCHEMA = "pipeline-case-completion-v2"
COMPLETION_FILENAME = "pipeline_completion.json"
PRESSURE_HEADER = struct.Struct("<16sIIQQQd")
PRESSURE_FRAME_STEP_BYTES = struct.calcsize("<Q")
PRESSURE_MAGIC = {b"MPM_PRESSURE_V1\0", b"MPM_PRESSURE_V2\0"}

# These names are deliberately independent of manifest.json. A preparation job may
# create a new labelled group only after this submitter has already returned.
CASE_CONFIGS = {
    "EQ_LS": "01_EQ_LS.json",
    "EQ_HS": "01_EQ_HS.json",
    "LS": "02_LS.json",
    "HS": "02_HS.json",
    "HM": "02_HM.json",
    "HD": "02_HD.json",
    "RL": "04_RL.json",
    "RM": "04_RM.json",
    "RE": "04_RE.json",
}


@dataclass(frozen=True)
class SubmittedJob:
    """One queued Slurm job and the dependencies used at submission time."""

    key: str
    job_id: str
    job_name: str
    dependency_keys: list[str]
    dependency_ids: list[str]
    command: list[str]


def _unique(values: Iterable[str]) -> list[str]:
    return list(dict.fromkeys(value for value in values if value))


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def checked_case_path(path: Path) -> Path:
    resolved = path.resolve()
    try:
        resolved.relative_to(CASE_DIR.resolve())
    except ValueError as error:
        raise ValueError(f"Recorded artifact escapes the case directory: {path}") from error
    return resolved


def resolve_case_path(value: str | Path) -> Path:
    path = Path(value)
    return checked_case_path(path if path.is_absolute() else CASE_DIR / path)


def verify_artifact(record: Any, expected_path: Path) -> None:
    if not isinstance(record, dict):
        raise ValueError("Completion artifact record is not an object")
    expected_path = checked_case_path(expected_path)
    recorded_path = resolve_case_path(record["path"])
    if recorded_path != expected_path:
        raise ValueError("Completion artifact path differs from the configured path")
    if not expected_path.is_file():
        raise FileNotFoundError(expected_path)
    size = expected_path.stat().st_size
    if size <= 0 or int(record["size_bytes"]) != size:
        raise ValueError("Completion artifact size differs from the sentinel")
    if str(record["sha256"]) != file_sha256(expected_path):
        raise ValueError("Completion artifact hash differs from the sentinel")


def pressure_database_header(values_path: Path) -> dict[str, Any]:
    values_path = checked_case_path(values_path)
    with values_path.open("rb") as stream:
        raw = stream.read(PRESSURE_HEADER.size)
    if len(raw) != PRESSURE_HEADER.size:
        raise ValueError("Pressure database header is truncated")
    magic, dimension, values_per_particle, count, interval, max_step, source_dt = (
        PRESSURE_HEADER.unpack(raw)
    )
    if magic not in PRESSURE_MAGIC or dimension not in (2, 3):
        raise ValueError("Pressure database header is invalid")
    if values_per_particle != 2 or count == 0 or interval == 0 or source_dt <= 0.0:
        raise ValueError("Pressure database header has invalid dimensions")
    frame_count = max_step // interval + 1
    frame_bytes = PRESSURE_FRAME_STEP_BYTES + 2 * count * 8
    expected_size = PRESSURE_HEADER.size + frame_count * frame_bytes
    if values_path.stat().st_size != expected_size:
        raise ValueError("Pressure database has an incomplete frame payload")
    return {
        "format_version": "V1" if magic == b"MPM_PRESSURE_V1\0" else "V2",
        "dimension": dimension,
        "particle_count": count,
        "step_interval": interval,
        "max_step": max_step,
        "source_dt_s": source_dt,
        "frame_count": frame_count,
    }


def validate_pressure_points(points_path: Path, header: dict[str, Any]) -> None:
    points_path = checked_case_path(points_path)
    identifiers: set[int] = set()
    count = 0
    with points_path.open(encoding="utf-8") as stream:
        for line in stream:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            values = stripped.replace(",", " ").split()
            if len(values) < int(header["dimension"]) + 1:
                raise ValueError("Pressure point row has too few columns")
            raw_identifier = float(values[0])
            identifier = int(raw_identifier)
            if not math.isfinite(raw_identifier) or raw_identifier != identifier:
                raise ValueError("Pressure point identifier is not an integer")
            coordinates = [float(value) for value in values[1 : int(header["dimension"]) + 1]]
            if not all(math.isfinite(value) for value in coordinates):
                raise ValueError("Pressure point coordinate is not finite")
            if identifier in identifiers:
                raise ValueError("Pressure point identifiers are duplicated")
            identifiers.add(identifier)
            count += 1
    if count != int(header["particle_count"]):
        raise ValueError("Pressure point count differs from the binary header")


def verify_pressure_database(config: dict[str, Any], record: Any) -> None:
    if not isinstance(record, dict):
        raise ValueError("Pressure database audit is absent")
    pressure = config["analysis"]["prescribed_phase_pressures"]
    directory = resolve_case_path(pressure["path"])
    prefix = str(pressure["file_prefix"])
    points_path = directory / f"{prefix}_points.txt"
    values_path = directory / f"{prefix}_values.bin"
    verify_artifact(record["points"], points_path)
    verify_artifact(record["values"], values_path)
    observed = pressure_database_header(values_path)
    validate_pressure_points(points_path, observed)
    if record.get("header") != observed:
        raise ValueError("Pressure header differs from its completion audit")
    if observed["step_interval"] != int(pressure["step_interval"]):
        raise ValueError("Pressure step interval differs from config")
    if bool(pressure.get("write")):
        if observed["max_step"] != int(pressure["max_step"]):
            raise ValueError("Written pressure max step differs from config")
    elif observed["max_step"] < int(pressure["max_step"]):
        raise ValueError("Read pressure database ends before the configured replay")
    if not math.isclose(
        float(observed["source_dt_s"]),
        float(pressure["source_dt"]),
        rel_tol=0.0,
        abs_tol=1.0e-15,
    ):
        raise ValueError("Pressure source dt differs from config")
    if prefix == "phase_erased":
        metadata_path = directory / f"{prefix}_metadata.json"
        verify_artifact(record["metadata"], metadata_path)


def resume_checkpoint_paths(config: dict[str, Any]) -> tuple[Path, Path]:
    resume = config["analysis"]["resume"]
    result_base = resolve_case_path(
        config["post_processing"].get("path", "results/")
    )
    digits = len(str(int(resume["nsteps"])))
    checkpoint = (
        result_base
        / str(resume["uuid"])
        / f"particles{int(resume['step']):0{digits}d}.h5"
    )
    qa_vtp = checkpoint.with_name(
        checkpoint.name.replace("particles", "particle")
    ).with_suffix(".vtp")
    return checked_case_path(checkpoint), checked_case_path(qa_vtp)


def verify_runtime_dependencies(config: dict[str, Any], record: Any) -> None:
    if not isinstance(record, dict):
        raise ValueError("Runtime dependency audit is absent")
    expected_keys: set[str] = set()
    resume = config["analysis"].get("resume", {})
    if bool(resume.get("resume")):
        expected_keys.add("resume_equilibrium")
        checkpoint, qa_vtp = resume_checkpoint_paths(config)
        resume_record = record["resume_equilibrium"]
        verify_artifact(resume_record["checkpoint_hdf5"], checkpoint)
        verify_artifact(resume_record["qa_vtp"], qa_vtp)
    pressure = config["analysis"].get("prescribed_phase_pressures")
    if pressure is not None and bool(pressure.get("enable")):
        expected_keys.add("read_pressure_database")
        verify_pressure_database(config, record["read_pressure_database"])
    if set(record) != expected_keys:
        raise ValueError("Runtime dependency audit contains missing or unexpected inputs")


def config_complete(config_path: Path) -> bool:
    """Return true only for an audited, configuration-matched completion."""

    if not config_path.is_file():
        return False
    try:
        config_path = checked_case_path(config_path)
        config = json.loads(config_path.read_text(encoding="utf-8"))
        analysis = config["analysis"]
        nsteps = int(analysis["nsteps"])
        digits = len(str(nsteps))
        result_base = resolve_case_path(
            config["post_processing"].get("path", "results/")
        )
        result_directory = result_base / str(analysis["uuid"])
        final_vtp = result_directory / f"particle{nsteps:0{digits}d}.vtp"
        sentinel_path = result_directory / COMPLETION_FILENAME
        sentinel = json.loads(sentinel_path.read_text(encoding="utf-8"))
        if sentinel.get("schema") != COMPLETION_SCHEMA:
            raise ValueError("Completion sentinel schema is invalid")
        sentinel_config = sentinel["config"]
        if resolve_case_path(sentinel_config["path"]) != config_path:
            raise ValueError("Completion sentinel references another config")
        if str(sentinel_config["sha256"]) != file_sha256(config_path):
            raise ValueError("Completion sentinel config hash is stale")
        if str(sentinel_config["uuid"]) != str(analysis["uuid"]):
            raise ValueError("Completion sentinel UUID is stale")
        if int(sentinel_config["nsteps"]) != nsteps:
            raise ValueError("Completion sentinel nsteps is stale")
        artifacts = sentinel["artifacts"]
        verify_artifact(artifacts["final_vtp"], final_vtp)
        if config["post_processing"].get("write_hdf5"):
            final_hdf5 = result_directory / f"particles{nsteps:0{digits}d}.h5"
            verify_artifact(artifacts["final_hdf5"], final_hdf5)
        pressure = analysis.get("prescribed_phase_pressures")
        if pressure is not None and bool(pressure.get("write")):
            verify_pressure_database(config, artifacts["written_pressure_database"])
        verify_runtime_dependencies(config, sentinel["runtime_dependencies"])
    except (
        KeyError,
        TypeError,
        ValueError,
        json.JSONDecodeError,
        OSError,
        OverflowError,
    ):
        return False
    return True


def phase_control_complete(group: str, fit_end: float | str | None = None) -> bool:
    """Verify output, source provenance, binary shape and transform metadata."""

    try:
        root = resolve_case_path(Path("pressure_databases") / group)
        source_points = root / "lagged" / "pressure_points.txt"
        source_values = root / "lagged" / "pressure_values.bin"
        output_points = root / "phase_erased" / "phase_erased_points.txt"
        output_values = root / "phase_erased" / "phase_erased_values.bin"
        metadata_path = root / "phase_erased" / "phase_erased_metadata.json"
        metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
        if metadata.get("schema") != "phase-erased-pressure-control-v1":
            raise ValueError("Unexpected phase metadata schema")
        if metadata.get("classification") != (
            "one-way numerical counterfactual; not fully coupled"
        ):
            raise ValueError("Unexpected phase-control classification")
        expected_paths = {
            "source": {"points": source_points, "values": source_values},
            "output": {"points": output_points, "values": output_values},
        }
        for section, paths in expected_paths.items():
            for name, path in paths.items():
                if resolve_case_path(metadata[section][name]) != checked_case_path(path):
                    raise ValueError("Phase metadata records an unexpected path")
                if not path.is_file() or path.stat().st_size <= 0:
                    raise FileNotFoundError(path)
                if metadata[section][f"{name}_sha256"] != file_sha256(path):
                    raise ValueError("Phase database hash differs from metadata")

        source_header = pressure_database_header(source_values)
        output_header = pressure_database_header(output_values)
        validate_pressure_points(source_points, source_header)
        validate_pressure_points(output_points, output_header)
        if source_header != output_header:
            raise ValueError("Lagged and phase-erased database headers differ")
        recorded_header = metadata["header"]
        for key in (
            "format_version",
            "dimension",
            "particle_count",
            "step_interval",
            "max_step",
            "source_dt_s",
        ):
            if recorded_header.get(key) != source_header[key]:
                raise ValueError(f"Phase metadata header field {key} is stale")

        transform = metadata["transform"]
        expected_fit_end = (
            float(fit_end)
            if fit_end is not None
            else (10.4 if group.startswith("screen") else 23.4)
        )
        expected_transform = {
            "period_s": 1.3,
            "fit_start_s": 2.6,
            "fit_end_s": expected_fit_end,
            "ramp_time_s": 1.3,
            "surface_band_m": 0.015,
            "minimum_reference_amplitude_pa": 1.0,
        }
        for key, expected in expected_transform.items():
            if not math.isclose(
                float(transform[key]), expected, rel_tol=0.0, abs_tol=1.0e-12
            ):
                raise ValueError(f"Phase transform field {key} is stale")
    except (
        KeyError,
        TypeError,
        ValueError,
        json.JSONDecodeError,
        OSError,
        OverflowError,
    ):
        return False
    return True


def sbatch_prefix(
    args: argparse.Namespace, job_name: str, walltime: str
) -> list[str]:
    """Build the immutable HPC4 resource portion of every sbatch command."""

    command = [
        "sbatch",
        "--parsable",
        "--partition=granularmech",
        "--account=comgranmech",
        "--nodes=1",
        "--ntasks=1",
        f"--cpus-per-task={args.cpus}",
        f"--time={walltime}",
        f"--job-name={job_name}",
        f"--output=logs/{job_name}-%j.out",
        f"--error=logs/{job_name}-%j.err",
    ]
    if args.gres:
        command.append(f"--gres={args.gres}")
    if args.nodelist:
        command.append(f"--nodelist={args.nodelist}")
    return command


def sbatch_command(
    args: argparse.Namespace,
    job_name: str,
    walltime: str,
    payload: list[str],
    dependency_ids: Iterable[str] = (),
) -> list[str]:
    command = sbatch_prefix(args, job_name, walltime)
    dependencies = _unique(dependency_ids)
    if dependencies:
        command.append(f"--dependency=afterok:{':'.join(dependencies)}")
    command.extend((RUN_TASK, *payload))
    return command


def submit_sbatch(command: list[str]) -> str:
    """Submit one job and parse ``sbatch --parsable`` output."""

    completed = subprocess.run(
        command,
        cwd=CASE_DIR,
        check=True,
        capture_output=True,
        text=True,
        env=os.environ.copy(),
    )
    output = completed.stdout.strip()
    match = re.fullmatch(r"([0-9]+)(?:;[^\s]+)?", output)
    if match is None:
        raise RuntimeError(f"Unexpected sbatch --parsable output: {output!r}")
    return match.group(1)


class Submitter:
    """Queue jobs while retaining an auditable dependency graph."""

    def __init__(self, args: argparse.Namespace):
        self.args = args
        self.jobs: dict[str, SubmittedJob] = {}
        self.skipped: dict[str, str] = {}
        self._dry_job_id = 9000000

    def submit(
        self,
        key: str,
        job_name: str,
        walltime: str,
        payload: list[str],
        dependency_keys: Iterable[str] = (),
    ) -> str:
        if key in self.jobs:
            raise ValueError(f"Duplicate workflow job key: {key}")
        keys = _unique(dependency_keys)
        missing = [dependency for dependency in keys if dependency not in self.jobs]
        if missing:
            raise KeyError(f"Unknown dependency keys for {key}: {missing}")
        dependency_ids = [self.jobs[dependency].job_id for dependency in keys]
        command = sbatch_command(
            self.args, job_name, walltime, payload, dependency_ids
        )
        print("+", shlex.join(command), flush=True)
        if self.args.dry_run:
            self._dry_job_id += 1
            job_id = str(self._dry_job_id)
        else:
            job_id = submit_sbatch(command)
        record = SubmittedJob(
            key=key,
            job_id=job_id,
            job_name=job_name,
            dependency_keys=keys,
            dependency_ids=dependency_ids,
            command=command,
        )
        self.jobs[key] = record
        print(f"SUBMITTED {key}: {job_id}", flush=True)
        return key

    def skip(self, key: str, reason: str) -> None:
        self.skipped[key] = reason
        print(f"SKIP {key}: {reason}", flush=True)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("tier", choices=("screen", "production"))
    parser.add_argument("--label", default="baseline")
    parser.add_argument("--stage", choices=("physical", "full"), default="full")
    parser.add_argument("--cpus", type=int, default=16)
    parser.add_argument("--gres", default="gpu:1")
    parser.add_argument("--nodelist", default="gpu30")
    parser.add_argument("--prepare", action="store_true")
    parser.add_argument("--build", action="store_true")
    parser.add_argument("--analyze", action="store_true")
    parser.add_argument("--mesh-dir", default=".")
    parser.add_argument("--cell-size", type=float, default=0.02)
    parser.add_argument("--particle-spacing", type=float, default=0.01)
    parser.add_argument("--wave-height", type=float, default=0.12)
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    return parser


def _case_needs_run(args: argparse.Namespace, config_path: Path) -> bool:
    # Preparing a group may alter its geometry or loading while retaining UUIDs.
    return args.force or args.prepare or not config_complete(config_path)


def submit_workflow(
    args: argparse.Namespace, submitter: Submitter | None = None
) -> Submitter:
    """Submit the DAG and return all job IDs without waiting for computation."""

    if args.cpus <= 0:
        raise ValueError("--cpus must be positive")
    if not re.fullmatch(r"[A-Za-z0-9_-]+", args.label):
        raise ValueError("--label may contain only letters, digits, '_' and '-'")
    if len(args.label) > 32:
        raise ValueError("--label must contain no more than 32 characters")

    suffix = "" if args.label == "baseline" else f"_{args.label}"
    group = f"{args.tier}{suffix}"
    config_dir = CASE_DIR / "configs" / group
    walltime = "12:00:00" if args.tier == "screen" else "36:00:00"
    fit_end = "10.4" if args.tier == "screen" else "23.4"
    prefix = f"plp_{group}"
    queue = submitter or Submitter(args)

    if args.prepare:
        setup_key = queue.submit(
            "prepare",
            f"{prefix}_prepare",
            "02:00:00",
            [
                "bash",
                "prepare_inputs.sh",
                args.tier,
                args.label,
                args.mesh_dir,
                str(args.cell_size),
                str(args.particle_spacing),
                str(args.wave_height),
            ],
        )
    else:
        setup_key = queue.submit(
            "validate",
            f"{prefix}_validate",
            "01:00:00",
            ["bash", "validate_inputs.sh", group],
        )

    base_key = setup_key
    if args.build:
        base_key = queue.submit(
            "build",
            "plp_build",
            "04:00:00",
            ["bash", "build_on_server.sh"],
            [setup_key],
        )

    def submit_case(
        code: str,
        dependencies: Iterable[str],
        *,
        upstream_dirty: bool = False,
    ) -> str | None:
        config_path = config_dir / CASE_CONFIGS[code]
        if not upstream_dirty and not _case_needs_run(args, config_path):
            queue.skip(code, f"complete result for {config_path}")
            return None
        config_argument = config_path.relative_to(CASE_DIR).as_posix()
        return queue.submit(
            code,
            f"{prefix}_{code}",
            walltime,
            ["bash", "run_case.sh", config_argument],
            dependencies,
        )

    eq_ls = submit_case("EQ_LS", [base_key])
    eq_hs = submit_case("EQ_HS", [base_key])
    low_parent = eq_ls or base_key
    high_parent = eq_hs or base_key

    submit_case("LS", [low_parent], upstream_dirty=eq_ls is not None)
    submit_case("HS", [high_parent], upstream_dirty=eq_hs is not None)
    submit_case("HM", [high_parent], upstream_dirty=eq_hs is not None)

    if args.stage == "full":
        hd = submit_case("HD", [high_parent], upstream_dirty=eq_hs is not None)
        driver_parent = hd or high_parent

        phase_required = (
            args.force
            or args.prepare
            or hd is not None
            or not phase_control_complete(group, fit_end)
        )
        if phase_required:
            phase = queue.submit(
                "phase",
                f"{prefix}_phase",
                "04:00:00",
                ["bash", "run_phase_control.sh", group, fit_end],
                [driver_parent],
            )
        else:
            queue.skip("phase", f"complete phase-erased database for {group}")
            phase = None

        # RL and RM intentionally share the lagged HD pressure history. They can
        # run concurrently once the driver database is complete.
        submit_case("RL", [driver_parent], upstream_dirty=hd is not None)
        submit_case("RM", [driver_parent], upstream_dirty=hd is not None)
        submit_case(
            "RE",
            [phase or driver_parent],
            upstream_dirty=phase is not None,
        )

    if args.analyze:
        # Depend on every submitted job, including setup/build and intermediates.
        # This remains correct when restart-aware skips remove any case from the DAG.
        queue.submit(
            "analysis",
            f"{prefix}_analyze",
            "02:00:00",
            ["bash", "analyze_results.sh", group, args.stage],
            list(queue.jobs),
        )

    return queue


def save_submission(
    args: argparse.Namespace, queue: Submitter, output_path: Path
) -> None:
    payload: dict[str, Any] = {
        "schema": "pipeline-sbatch-submission-v1",
        "submitted_at_utc": datetime.now(timezone.utc).isoformat(),
        "tier": args.tier,
        "label": args.label,
        "stage": args.stage,
        "resources": {
            "partition": "granularmech",
            "account": "comgranmech",
            "cpus_per_task": args.cpus,
            "gres": args.gres,
            "nodelist": args.nodelist,
        },
        "jobs": [asdict(job) for job in queue.jobs.values()],
        "skipped": queue.skipped,
    }
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if not args.dry_run:
        (CASE_DIR / "logs").mkdir(exist_ok=True)
    queue = submit_workflow(args)
    if args.dry_run:
        print("Dry run only; no jobs were submitted and no submission file was written.")
        return 0

    timestamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    suffix = "" if args.label == "baseline" else f"_{args.label}"
    group = f"{args.tier}{suffix}"
    output_path = CASE_DIR / "submissions" / f"{group}-{timestamp}.json"
    save_submission(args, queue, output_path)
    print(f"Submission record: {output_path}")
    print("All jobs are queued; this command does not wait for computation.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
