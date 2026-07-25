#!/usr/bin/env python3
"""Run generated MPM cases locally or on a Linux server."""

from __future__ import annotations

import argparse
import csv
import os
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path


def run_one(
    row: dict[str, str],
    *,
    mpm_bin: Path,
    input_dir: Path,
    log_dir: Path,
    partitions: int,
    rerun: bool,
) -> tuple[str, str]:
    result_dir = input_dir / "results" / row["uuid"]
    if result_dir.exists() and not rerun:
        return row["case_id"], "SKIPPED (result directory exists)"

    case_path = Path(row["case_file"]).resolve()
    case_argument = case_path.relative_to(input_dir).as_posix()
    command = [
        str(mpm_bin),
        "-f",
        str(input_dir) + os.sep,
        "-i",
        case_argument,
        "-p",
        str(partitions),
    ]

    log_dir.mkdir(parents=True, exist_ok=True)
    log_path = log_dir / f"{row['case_id']}.log"
    with log_path.open("w", encoding="utf-8") as stream:
        stream.write("COMMAND: " + " ".join(command) + "\n\n")
        stream.flush()
        completed = subprocess.run(
            command,
            cwd=input_dir,
            stdout=stream,
            stderr=subprocess.STDOUT,
            check=False,
        )
    if completed.returncode != 0:
        return row["case_id"], f"FAILED ({completed.returncode}); see {log_path}"
    return row["case_id"], "COMPLETED"


def main() -> None:
    script_dir = Path(__file__).resolve().parent
    input_dir = script_dir.parent

    parser = argparse.ArgumentParser()
    parser.add_argument("--mpm-bin", type=Path, required=True)
    parser.add_argument("--jobs", type=int, default=1)
    parser.add_argument("--partitions", type=int, default=24)
    parser.add_argument(
        "--family",
        choices=("all", "permeability", "gas_boundary"),
        default="all",
    )
    parser.add_argument(
        "--rerun",
        action="store_true",
        help="Run even when the result directory already exists",
    )
    parser.add_argument(
        "--skip-legacy-reference",
        action="store_true",
        help=(
            "Do not rerun the existing k/k0=1 fixed-gas and "
            "co-pressurized variable-gas reference cases"
        ),
    )
    args = parser.parse_args()

    if args.jobs < 1:
        raise ValueError("--jobs must be at least 1")
    mpm_bin = args.mpm_bin.resolve()
    if not mpm_bin.is_file():
        raise FileNotFoundError(mpm_bin)

    manifest_path = script_dir / "manifest.csv"
    with manifest_path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    if args.family != "all":
        rows = [
            row for row in rows if args.family in row["families"].split("+")
        ]
    if args.skip_legacy_reference:
        rows = [
            row
            for row in rows
            if not (
                float(row["permeability_factor"]) == 1.0
                and (
                    row["gas_model"] == "fixed"
                    or (
                        row["gas_model"] == "variable"
                        and float(row["surface_gas_pressure_ratio"]) == 1.0
                    )
                )
            )
        ]

    print(f"Running {len(rows)} cases with {args.jobs} concurrent process(es)")
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        futures = [
            pool.submit(
                run_one,
                row,
                mpm_bin=mpm_bin,
                input_dir=input_dir,
                log_dir=script_dir / "logs",
                partitions=args.partitions,
                rerun=args.rerun,
            )
            for row in rows
        ]
        failed = 0
        for future in as_completed(futures):
            case_id, status = future.result()
            print(f"{case_id}: {status}", flush=True)
            if status.startswith("FAILED"):
                failed += 1
    if failed:
        raise SystemExit(f"{failed} case(s) failed")


if __name__ == "__main__":
    main()
