"""Submit the open and closed production segment chains through Slurm."""

from __future__ import annotations

import argparse
import csv
import shutil
import subprocess
import sys
from pathlib import Path


CASE_DIR = Path(__file__).resolve().parent


def load_manifest(path: Path):
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--cases",
        nargs="+",
        choices=("open", "closed"),
        default=("open", "closed"),
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="print the dependency chains without calling sbatch",
    )
    parser.add_argument(
        "--resume",
        action="store_true",
        help="skip a contiguous prefix of segments with completion markers",
    )
    args = parser.parse_args()

    subprocess.run(
        [sys.executable, str(CASE_DIR / "prepare_production_segments.py")],
        cwd=CASE_DIR,
        check=True,
    )
    manifest_path = CASE_DIR / "production_inputs" / "manifest.csv"
    rows = [row for row in load_manifest(manifest_path) if row["case"] in args.cases]
    if not rows:
        parser.error("no selected production rows were generated")

    existing = [
        row["result_dir"]
        for row in rows
        if (CASE_DIR / row["result_dir"]).exists()
    ]
    if existing and not args.resume:
        parser.error(
            "refusing to mix a new chain with existing result directories: "
            + ", ".join(existing[:3])
            + "; inspect them, then use --resume only for completed prefixes"
        )
    if not args.dry_run and shutil.which("sbatch") is None:
        parser.error("sbatch is not available; run this script on an HPC4 login node")

    submitted = []
    for case_name in args.cases:
        case_rows = [row for row in rows if row["case"] == case_name]
        pending_rows = []
        found_gap = False
        for row in case_rows:
            result_dir = CASE_DIR / row["result_dir"]
            if result_dir.exists():
                if found_gap:
                    parser.error(
                        f"{case_name} has an existing segment after a missing one: "
                        f"{result_dir}"
                    )
                required = [
                    result_dir / "SEGMENT_COMPLETED.txt",
                    CASE_DIR / row["final_particle"],
                    CASE_DIR / row["final_checkpoint"],
                ]
                missing = [str(path) for path in required if not path.is_file()]
                if missing:
                    parser.error(
                        f"{case_name} segment {row['segment_index']} exists but is "
                        "incomplete; inspect or move it before resubmission: "
                        + ", ".join(missing)
                    )
                if not args.resume:
                    parser.error(f"existing result directory: {result_dir}")
                print(
                    f"Skipping completed {case_name} segment "
                    f"{row['segment_index']}"
                )
            else:
                found_gap = True
                pending_rows.append(row)

        if not pending_rows:
            print(f"All {case_name} production segments are already complete")
            continue

        previous_job = None
        for row in pending_rows:
            command = ["sbatch", "--parsable"]
            if previous_job is not None:
                command.append(f"--dependency=afterok:{previous_job}")
            command.extend(
                [
                    f"--export=ALL,PRODUCTION_ROW={row['row_index']}",
                    "run_production_segment_hpc4.sh",
                ]
            )
            if args.dry_run:
                job_id = f"dry-{case_name}-{int(row['segment_index']):03d}"
                print("DRY RUN:", " ".join(command))
            else:
                result = subprocess.run(
                    command,
                    cwd=CASE_DIR,
                    check=True,
                    text=True,
                    capture_output=True,
                )
                job_id = result.stdout.strip().split(";", 1)[0]
                if not job_id.isdigit():
                    raise RuntimeError(f"Unexpected sbatch response: {result.stdout!r}")
                print(
                    f"Submitted {case_name} segment {row['segment_index']}: "
                    f"job {job_id}"
                )
            submitted.append(
                {
                    "case": case_name,
                    "segment_index": row["segment_index"],
                    "row_index": row["row_index"],
                    "job_id": job_id,
                    "dependency": previous_job or "none",
                }
            )
            previous_job = job_id

    if not args.dry_run:
        if not submitted:
            print("Nothing new was submitted.")
            return
        jobs_path = CASE_DIR / "production_inputs" / "submitted_jobs.csv"
        with jobs_path.open("w", encoding="utf-8", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(submitted[0]))
            writer.writeheader()
            writer.writerows(submitted)
        print(f"Wrote {jobs_path}")
    print("Open and closed chains are independent; each chain is sequential.")


if __name__ == "__main__":
    main()
