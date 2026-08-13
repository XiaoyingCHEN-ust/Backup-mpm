#!/usr/bin/env python3
"""Run the formal, non-Mohr-Coulomb screen DAG on a local Linux host.

This controller is deliberately fixed to the registered baseline ``screen``
study.  It never prepares inputs, deletes outputs, clears sentinels, or selects
any Mohr-Coulomb case.  Incomplete non-empty outputs are therefore a hard stop
for manual audit rather than an invitation to overwrite them.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Mapping, Protocol

import analyze_study
import submit_study_sbatch
import validate_case


CASE_ROOT = Path(__file__).resolve().parent
GROUP = "screen"
FIT_END_S = "10.4"


class SafetyError(RuntimeError):
    """Raised when continuing could overwrite or trust unaudited state."""


class PairRunError(RuntimeError):
    """Raised after every process in one logical pair has been waited."""


class WaitableProcess(Protocol):
    pid: int

    def wait(self) -> int: ...


@dataclass(frozen=True)
class CaseSpec:
    key: str
    filename: str
    uuid: str
    material: str
    dependencies: tuple[str, ...] = ()


CASE_SPECS: tuple[CaseSpec, ...] = (
    CaseSpec("EQ_LS", "01_EQ_LS.json", "PLP_SCREEN_LS_EQ", "LinearElastic2D"),
    CaseSpec("EQ_HS", "01_EQ_HS.json", "PLP_SCREEN_HS_EQ", "LinearElastic2D"),
    CaseSpec("LS", "02_LS.json", "PLP_SCREEN_LS_SANI", "SANISAND2D", ("EQ_LS",)),
    CaseSpec("HS", "02_HS.json", "PLP_SCREEN_HS_SANI", "SANISAND2D", ("EQ_HS",)),
    CaseSpec(
        "HD", "02_HD.json", "PLP_SCREEN_HD_DRIVER", "SANISAND2D", ("EQ_HS",)
    ),
    CaseSpec("RL", "04_RL.json", "PLP_SCREEN_RL_SANI", "SANISAND2D", ("HD",)),
    CaseSpec("RE", "04_RE.json", "PLP_SCREEN_RE_SANI", "SANISAND2D", ("HD",)),
)
SPEC_BY_KEY = {spec.key: spec for spec in CASE_SPECS}
RUN_STAGES: tuple[tuple[str, ...], ...] = (
    ("EQ_LS", "EQ_HS"),
    ("LS", "HS"),
    ("HD",),
)
REPLAY_STAGE = ("RL", "RE")


def _positive_integer(value: str, name: str) -> int:
    if not value.isascii() or not value.isdigit() or value.startswith("0"):
        raise SafetyError(f"{name} must be a positive base-10 integer, got {value!r}")
    number = int(value)
    if number <= 0:
        raise SafetyError(f"{name} must be positive, got {value!r}")
    return number


def _checked_path(root: Path, value: str | Path) -> Path:
    path = Path(value)
    resolved = (path if path.is_absolute() else root / path).resolve()
    try:
        resolved.relative_to(root.resolve())
    except ValueError as error:
        raise SafetyError(f"Configured path escapes the pipeline directory: {path}") from error
    return resolved


def _directory_nonempty(path: Path) -> bool:
    if path.is_symlink():
        raise SafetyError(f"Refusing a symlinked output directory: {path}")
    if not path.exists():
        return False
    if not path.is_dir():
        raise SafetyError(f"Expected an output directory but found another file type: {path}")
    return next(path.iterdir(), None) is not None


def _default_formal_validate(config_path: Path) -> None:
    # No --local-smoke and no relaxed path: this is the registered validator.
    validate_case.validate_one(config_path, runtime=False)


def _default_exact_grid_verify(config_path: Path) -> None:
    config_path = config_path.resolve()
    config = json.loads(config_path.read_text(encoding="utf-8"))
    result = _result_directory(CASE_ROOT, config)
    files = sorted(result.glob("particle*.vtp"), key=analyze_study.particle_step)
    expected = analyze_study.expected_particle_steps(config)
    actual = [analyze_study.particle_step(path) for path in files]
    if actual != expected:
        raise SafetyError(
            f"{config_path.name}: analyzer particle grid differs from config; "
            f"expected={expected}, found={actual}"
        )
    analyze_study.validate_completed_result(
        config_path, CASE_ROOT, config, result, files
    )


def _result_directory(root: Path, config: dict[str, Any]) -> Path:
    base = _checked_path(root, config["post_processing"].get("path", "results/"))
    return _checked_path(root, base / str(config["analysis"]["uuid"]))


class LocalPhaseScreenController:
    """Fail-closed orchestration for the seven registered phase-only cases."""

    def __init__(
        self,
        *,
        case_root: Path = CASE_ROOT,
        environment: Mapping[str, str] | None = None,
        formal_validate: Callable[[Path], None] | None = None,
        config_complete: Callable[[Path], bool] | None = None,
        exact_grid_verify: Callable[[Path], None] | None = None,
        phase_complete: Callable[[], bool] | None = None,
        case_launcher: Callable[[CaseSpec, dict[str, str]], WaitableProcess]
        | None = None,
        phase_runner: Callable[[dict[str, str]], int] | None = None,
    ) -> None:
        self.case_root = case_root.resolve()
        self.environment = dict(os.environ if environment is None else environment)
        self.formal_validate = formal_validate or _default_formal_validate
        self.config_complete = config_complete or submit_study_sbatch.config_complete
        self.exact_grid_verify = exact_grid_verify or _default_exact_grid_verify
        self.phase_complete = phase_complete or (
            lambda: submit_study_sbatch.phase_control_complete(GROUP, FIT_END_S)
        )
        self.case_launcher = case_launcher or self._launch_case
        self.phase_runner = phase_runner or self._run_phase_transform
        self.completed: set[str] = set()
        self.configs: dict[str, dict[str, Any]] = {}

        self.threads = _positive_integer(
            self.environment.get("MPM_THREADS", "4"), "MPM_THREADS"
        )
        host_threads = os.cpu_count() or 1
        self.total_thread_limit = _positive_integer(
            self.environment.get("PHASE_SCREEN_MAX_TOTAL_THREADS", str(host_threads)),
            "PHASE_SCREEN_MAX_TOTAL_THREADS",
        )
        if self.threads > self.total_thread_limit:
            raise SafetyError(
                f"MPM_THREADS={self.threads} exceeds the safe total-thread limit "
                f"{self.total_thread_limit}"
            )
        self.parallel_limit = min(2, self.total_thread_limit // self.threads)
        if self.parallel_limit < 1:
            raise SafetyError("The total-thread limit cannot run even one case")

    def config_path(self, spec: CaseSpec) -> Path:
        return self.case_root / "configs" / GROUP / spec.filename

    def _load_and_check_manifest(self) -> None:
        manifest_path = self.case_root / "configs" / GROUP / "manifest.json"
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        if (
            manifest.get("schema") != "pipeline-phase-lag-study-v1"
            or manifest.get("tier") != GROUP
            or manifest.get("label") != "baseline"
        ):
            raise SafetyError("configs/screen/manifest.json is not the baseline screen manifest")
        cases = manifest.get("cases")
        if not isinstance(cases, list):
            raise SafetyError("The screen manifest case inventory is malformed")
        recorded: dict[str, list[str]] = {}
        for row in cases:
            if isinstance(row, dict) and isinstance(row.get("code"), str):
                recorded.setdefault(row["code"], []).append(str(row.get("config", "")))
        for spec in CASE_SPECS:
            if recorded.get(spec.key) != [spec.filename]:
                raise SafetyError(
                    f"The screen manifest does not bind {spec.key} to {spec.filename}"
                )

    def _load_and_check_config(self, spec: CaseSpec) -> dict[str, Any]:
        path = self.config_path(spec)
        config = json.loads(path.read_text(encoding="utf-8"))
        analysis = config["analysis"]
        material = config["materials"][0]
        if analysis.get("uuid") != spec.uuid:
            raise SafetyError(f"{spec.filename} has an unexpected registered UUID")
        if material.get("type") != spec.material:
            raise SafetyError(f"{spec.filename} has an unexpected material type")
        if "MohrCoulomb" in str(material.get("type", "")):
            raise SafetyError(f"Mohr-Coulomb is forbidden in this controller: {spec.filename}")
        result_base = _checked_path(
            self.case_root, config["post_processing"].get("path", "results/")
        )
        if result_base != (self.case_root / "results" / GROUP).resolve():
            raise SafetyError(f"{spec.filename} does not target results/screen")

        pressure = analysis.get("prescribed_phase_pressures")
        expected_pressure = {
            "HD": (False, True, "pressure_databases/screen/lagged", "pressure"),
            "RL": (True, False, "pressure_databases/screen/lagged", "pressure"),
            "RE": (
                True,
                False,
                "pressure_databases/screen/phase_erased",
                "phase_erased",
            ),
        }.get(spec.key)
        if expected_pressure is not None:
            if not isinstance(pressure, dict):
                raise SafetyError(f"{spec.filename} has no registered pressure control")
            observed = (
                bool(pressure.get("enable")),
                bool(pressure.get("write")),
                str(pressure.get("path")),
                str(pressure.get("file_prefix")),
            )
            if observed != expected_pressure:
                raise SafetyError(f"{spec.filename} pressure-control registration is stale")
        return config

    def _owned_write_directories(self, config: dict[str, Any]) -> list[Path]:
        pressure = config["analysis"].get("prescribed_phase_pressures")
        if isinstance(pressure, dict) and bool(pressure.get("write")):
            return [_checked_path(self.case_root, pressure["path"])]
        return []

    def _require_empty_pending_outputs(self, spec: CaseSpec) -> None:
        config = self.configs[spec.key]
        paths = [_result_directory(self.case_root, config)]
        paths.extend(self._owned_write_directories(config))
        for path in paths:
            if _directory_nonempty(path):
                raise SafetyError(
                    f"{spec.key} is not validly complete but output is non-empty: {path}. "
                    "Stop for manual audit; nothing was overwritten."
                )

    def preflight(self) -> None:
        """Validate every selected config and output before starting any process."""

        self._load_and_check_manifest()
        for spec in CASE_SPECS:
            config = self._load_and_check_config(spec)
            self.configs[spec.key] = config
            self.formal_validate(self.config_path(spec))

        for spec in CASE_SPECS:
            path = self.config_path(spec)
            if self.config_complete(path):
                self.exact_grid_verify(path)
                self.completed.add(spec.key)
                print(f"SKIP {spec.key}: valid completion and exact analyzer grid", flush=True)
            else:
                self._require_empty_pending_outputs(spec)

        phase_output = self.case_root / "pressure_databases" / GROUP / "phase_erased"
        if not self.phase_complete() and _directory_nonempty(phase_output):
            raise SafetyError(
                f"Phase transform is invalid but its output is non-empty: {phase_output}. "
                "Stop for manual audit; nothing was overwritten."
            )

    def _launch_case(self, spec: CaseSpec, environment: dict[str, str]) -> WaitableProcess:
        runner = Path(
            environment.get(
                "PHASE_SCREEN_CASE_RUNNER", str(self.case_root / "run_case_local.sh")
            )
        )
        relative_config = self.config_path(spec).relative_to(self.case_root)
        return subprocess.Popen(
            ["bash", str(runner), str(relative_config)],
            cwd=self.case_root,
            env=environment,
        )

    def _run_phase_transform(self, environment: dict[str, str]) -> int:
        runner = Path(
            environment.get(
                "PHASE_SCREEN_PHASE_RUNNER", str(self.case_root / "run_phase_control.sh")
            )
        )
        completed = subprocess.run(
            ["bash", str(runner), GROUP, FIT_END_S],
            cwd=self.case_root,
            env=environment,
            check=False,
        )
        return int(completed.returncode)

    def _case_environment(self) -> dict[str, str]:
        environment = dict(self.environment)
        environment["MPM_THREADS"] = str(self.threads)
        environment["OMP_NUM_THREADS"] = str(self.threads)
        environment["TBB_NUM_THREADS"] = str(self.threads)
        environment["PIPELINE_LOCAL"] = "1"
        return environment

    def _check_dependencies(self, spec: CaseSpec) -> None:
        missing = [key for key in spec.dependencies if key not in self.completed]
        if missing:
            raise SafetyError(f"{spec.key} dependencies are incomplete: {missing}")

    def run_pair(self, keys: tuple[str, ...]) -> None:
        """Run a logical pair, wait every PID, then fail if either did not pass."""

        pending = [SPEC_BY_KEY[key] for key in keys if key not in self.completed]
        for spec in pending:
            self._check_dependencies(spec)
            # Close the preflight-to-launch race without deleting anything.
            self._require_empty_pending_outputs(spec)
        if not pending:
            return

        failures: list[str] = []
        environment = self._case_environment()
        for offset in range(0, len(pending), self.parallel_limit):
            chunk = pending[offset : offset + self.parallel_limit]
            running: list[tuple[CaseSpec, WaitableProcess]] = []
            for spec in chunk:
                try:
                    process = self.case_launcher(spec, dict(environment))
                    running.append((spec, process))
                    print(f"START {spec.key}: pid={process.pid}", flush=True)
                except Exception as error:  # launch failure must not hide its pair mate
                    failures.append(f"{spec.key} launch failed: {error}")

            statuses: list[tuple[CaseSpec, int]] = []
            for spec, process in running:
                try:
                    status = int(process.wait())
                except Exception as error:
                    failures.append(f"{spec.key} wait failed for pid {process.pid}: {error}")
                    continue
                statuses.append((spec, status))
                print(
                    f"WAIT {spec.key}: pid={process.pid} status={status}", flush=True
                )

            for spec, status in statuses:
                if status != 0:
                    failures.append(f"{spec.key} pid exited with status {status}")
                    continue
                path = self.config_path(spec)
                try:
                    if not self.config_complete(path):
                        raise SafetyError("config_complete returned false")
                    self.exact_grid_verify(path)
                    self.completed.add(spec.key)
                    print(f"PASS {spec.key}: completion and exact grid verified", flush=True)
                except Exception as error:
                    failures.append(f"{spec.key} post-run audit failed: {error}")

        if failures:
            raise PairRunError("; ".join(failures))

    def ensure_phase_transform(self) -> None:
        if "HD" not in self.completed:
            raise SafetyError("HD must be complete before the phase transform")
        if self.phase_complete():
            print("SKIP phase: existing transform passed provenance audit", flush=True)
            return
        output = self.case_root / "pressure_databases" / GROUP / "phase_erased"
        if _directory_nonempty(output):
            raise SafetyError(
                f"Refusing to overwrite invalid phase-transform output: {output}"
            )
        status = self.phase_runner(self._case_environment())
        if status != 0:
            raise PairRunError(f"phase transform exited with status {status}")
        if not self.phase_complete():
            raise SafetyError("Phase transform exited zero but provenance audit failed")
        print("PASS phase: transform and provenance verified", flush=True)

    def run(self) -> None:
        self.preflight()
        print(
            f"LOCAL PHASE SCREEN: {self.threads} threads/case, "
            f"parallel={self.parallel_limit}, total<={self.total_thread_limit}",
            flush=True,
        )
        for stage in RUN_STAGES:
            self.run_pair(stage)
        self.ensure_phase_transform()
        self.run_pair(REPLAY_STAGE)
        if self.completed != set(SPEC_BY_KEY):
            raise SafetyError("The seven-case phase-only DAG did not fully complete")
        print("PASS: registered seven-case local phase-only screen DAG", flush=True)


def main() -> int:
    if len(sys.argv) != 1:
        print("Usage: run_local_phase_screen.sh", file=sys.stderr)
        return 2
    try:
        LocalPhaseScreenController().run()
    except (SafetyError, PairRunError, KeyError, TypeError, ValueError, OSError) as error:
        print(f"LOCAL PHASE SCREEN FAILED: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
