from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path
from typing import Any


CASE_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(CASE_DIR))

import local_phase_screen as local_screen  # noqa: E402


KEY_BY_FILENAME = {spec.filename: spec.key for spec in local_screen.CASE_SPECS}


def write_registered_fixture(root: Path) -> None:
    config_dir = root / "configs" / "screen"
    config_dir.mkdir(parents=True)
    manifest_cases: list[dict[str, str]] = []
    for spec in local_screen.CASE_SPECS:
        analysis: dict[str, Any] = {"uuid": spec.uuid}
        if spec.key == "HD":
            analysis["prescribed_phase_pressures"] = {
                "enable": False,
                "write": True,
                "path": "pressure_databases/screen/lagged",
                "file_prefix": "pressure",
            }
        elif spec.key == "RL":
            analysis["prescribed_phase_pressures"] = {
                "enable": True,
                "write": False,
                "path": "pressure_databases/screen/lagged",
                "file_prefix": "pressure",
            }
        elif spec.key == "RE":
            analysis["prescribed_phase_pressures"] = {
                "enable": True,
                "write": False,
                "path": "pressure_databases/screen/phase_erased",
                "file_prefix": "phase_erased",
            }
        config = {
            "analysis": analysis,
            "materials": [{"type": spec.material}],
            "post_processing": {"path": "results/screen/"},
        }
        (config_dir / spec.filename).write_text(
            json.dumps(config), encoding="utf-8"
        )
        manifest_cases.append({"code": spec.key, "config": spec.filename})

    # The real registered manifest contains these MC cases.  Their presence in
    # the inventory must never make the phase-only controller invoke them.
    manifest_cases.extend(
        [
            {"code": "MC_EQ", "config": "02_MC_EQ.json"},
            {"code": "HM", "config": "03_HM.json"},
            {"code": "RM", "config": "04_RM.json"},
        ]
    )
    (config_dir / "manifest.json").write_text(
        json.dumps(
            {
                "schema": "pipeline-phase-lag-study-v1",
                "tier": "screen",
                "label": "baseline",
                "cases": manifest_cases,
            }
        ),
        encoding="utf-8",
    )


class FakeRun:
    def __init__(
        self,
        root: Path,
        *,
        initially_complete: set[str] | None = None,
        statuses: dict[str, int] | None = None,
        phase_valid: bool = False,
    ) -> None:
        self.root = root
        self.complete = set(initially_complete or ())
        self.statuses = dict(statuses or {})
        self.phase_valid = phase_valid
        self.events: list[tuple[Any, ...]] = []
        self.next_pid = 1000
        self.active = 0
        self.maximum_active = 0

    @staticmethod
    def key(path: Path) -> str:
        return KEY_BY_FILENAME[path.name]

    def formal_validate(self, path: Path) -> None:
        self.events.append(("formal", self.key(path)))

    def config_complete(self, path: Path) -> bool:
        key = self.key(path)
        result = key in self.complete
        self.events.append(("complete", key, result))
        return result

    def exact_grid_verify(self, path: Path) -> None:
        key = self.key(path)
        if key not in self.complete:
            raise AssertionError(f"verified before completion: {key}")
        self.events.append(("verify", key))

    def phase_complete(self) -> bool:
        self.events.append(("phase-check", self.phase_valid))
        return self.phase_valid

    def phase_runner(self, environment: dict[str, str]) -> int:
        self.events.append(
            (
                "phase-run",
                environment["PIPELINE_LOCAL"],
                environment["MPM_THREADS"],
            )
        )
        self.phase_valid = True
        return 0

    def launcher(
        self, spec: local_screen.CaseSpec, environment: dict[str, str]
    ) -> "FakeProcess":
        if spec.key in {"MC_EQ", "HM", "RM"}:
            raise AssertionError(f"MC case was invoked: {spec.key}")
        self.next_pid += 1
        self.active += 1
        self.maximum_active = max(self.maximum_active, self.active)
        self.events.append(
            ("start", spec.key, environment["MPM_THREADS"], self.next_pid)
        )
        return FakeProcess(self, spec.key, self.next_pid, self.statuses.get(spec.key, 0))


class FakeProcess:
    def __init__(self, run: FakeRun, key: str, pid: int, status: int) -> None:
        self.run = run
        self.key = key
        self.pid = pid
        self.status = status

    def wait(self) -> int:
        self.run.events.append(("wait", self.key, self.pid, self.status))
        self.run.active -= 1
        if self.status == 0:
            self.run.complete.add(self.key)
        return self.status


def controller(root: Path, fake: FakeRun, **environment: str):
    values = {"MPM_THREADS": "4", "PHASE_SCREEN_MAX_TOTAL_THREADS": "8"}
    values.update(environment)
    return local_screen.LocalPhaseScreenController(
        case_root=root,
        environment=values,
        formal_validate=fake.formal_validate,
        config_complete=fake.config_complete,
        exact_grid_verify=fake.exact_grid_verify,
        phase_complete=fake.phase_complete,
        case_launcher=fake.launcher,
        phase_runner=fake.phase_runner,
    )


class LocalPhaseScreenControllerTest(unittest.TestCase):
    def test_fake_run_obeys_seven_case_dag_two_process_limit_and_never_runs_mc(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            write_registered_fixture(root)
            fake = FakeRun(root)
            controller(root, fake).run()

        starts = [event[1] for event in fake.events if event[0] == "start"]
        self.assertEqual(
            starts, ["EQ_LS", "EQ_HS", "LS", "HS", "HD", "RL", "RE"]
        )
        self.assertLessEqual(fake.maximum_active, 2)
        self.assertEqual(fake.maximum_active, 2)
        self.assertFalse({"MC_EQ", "HM", "RM"} & set(starts))
        self.assertEqual(
            {event[1] for event in fake.events if event[0] == "verify"},
            set(KEY_BY_FILENAME.values()),
        )

        first_start = next(i for i, event in enumerate(fake.events) if event[0] == "start")
        prefix = fake.events[:first_start]
        self.assertEqual(sum(event[0] == "formal" for event in prefix), 7)
        self.assertEqual(sum(event[0] == "complete" for event in prefix), 7)
        phase_run = next(i for i, event in enumerate(fake.events) if event[0] == "phase-run")
        hd_wait = next(
            i
            for i, event in enumerate(fake.events)
            if event[0] == "wait" and event[1] == "HD"
        )
        rl_start = next(
            i
            for i, event in enumerate(fake.events)
            if event[0] == "start" and event[1] == "RL"
        )
        self.assertLess(hd_wait, phase_run)
        self.assertLess(phase_run, rl_start)

    def test_valid_case_completions_and_phase_transform_are_all_skipped(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            write_registered_fixture(root)
            fake = FakeRun(
                root,
                initially_complete=set(KEY_BY_FILENAME.values()),
                phase_valid=True,
            )
            controller(root, fake).run()

        self.assertFalse(any(event[0] == "start" for event in fake.events))
        self.assertFalse(any(event[0] == "phase-run" for event in fake.events))
        self.assertEqual(sum(event[0] == "verify" for event in fake.events), 7)
        # Preflight and post-HD both independently verify an existing transform.
        self.assertGreaterEqual(sum(event[0] == "phase-check" for event in fake.events), 2)

    def test_invalid_completion_with_nonempty_result_refuses_all_launches(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            write_registered_fixture(root)
            partial = root / "results" / "screen" / "PLP_SCREEN_LS_EQ"
            partial.mkdir(parents=True)
            (partial / "particle40000.vtp").write_text("partial", encoding="utf-8")
            fake = FakeRun(root)
            with self.assertRaisesRegex(local_screen.SafetyError, "manual audit"):
                controller(root, fake).run()
            self.assertFalse(any(event[0] == "start" for event in fake.events))

    def test_partial_phase_transform_refuses_all_launches(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            write_registered_fixture(root)
            partial = root / "pressure_databases" / "screen" / "phase_erased"
            partial.mkdir(parents=True)
            (partial / "phase_erased_values.bin").write_bytes(b"partial")
            fake = FakeRun(root, phase_valid=False)
            with self.assertRaisesRegex(local_screen.SafetyError, "manual audit"):
                controller(root, fake).run()
            self.assertFalse(any(event[0] == "start" for event in fake.events))
            self.assertFalse(any(event[0] == "phase-run" for event in fake.events))

    def test_first_pid_failure_is_not_hidden_by_second_pid_success(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            write_registered_fixture(root)
            fake = FakeRun(root, statuses={"EQ_LS": 17, "EQ_HS": 0})
            with self.assertRaisesRegex(local_screen.PairRunError, "EQ_LS.*17"):
                controller(root, fake).run()

        waits = [event[1] for event in fake.events if event[0] == "wait"]
        self.assertEqual(waits, ["EQ_LS", "EQ_HS"])
        self.assertIn(("verify", "EQ_HS"), fake.events)
        self.assertFalse(
            any(
                event[0] == "start" and event[1] in {"LS", "HS", "HD", "RL", "RE"}
                for event in fake.events
            )
        )

    def test_thread_override_never_exceeds_total_limit(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            write_registered_fixture(root)
            fake = FakeRun(root)
            controller(
                root,
                fake,
                MPM_THREADS="6",
                PHASE_SCREEN_MAX_TOTAL_THREADS="8",
            ).run()

        self.assertEqual(fake.maximum_active, 1)
        self.assertTrue(
            all(event[2] == "6" for event in fake.events if event[0] == "start")
        )


if __name__ == "__main__":
    unittest.main()
