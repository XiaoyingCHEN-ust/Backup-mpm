"""Generate restartable 50 s production segments for both Siemens columns."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path


CASE_DIR = Path(__file__).resolve().parent
DEFAULT_SEGMENT_DURATION_S = 50.0


def clone(source):
    return json.loads(json.dumps(source))


def output_name(prefix: str, step: int, max_steps: int, suffix: str) -> str:
    digits = len(str(max_steps))
    return f"{prefix}{step:0{digits}d}{suffix}"


def write_json(path: Path, config) -> None:
    with path.open("w", encoding="utf-8") as stream:
        json.dump(config, stream, indent=2)
        stream.write("\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--segment-duration",
        type=float,
        default=DEFAULT_SEGMENT_DURATION_S,
        help="nominal physical duration of each production job in seconds",
    )
    args = parser.parse_args()
    if args.segment_duration <= 0.0:
        parser.error("--segment-duration must be positive")

    output_dir = CASE_DIR / "production_inputs"
    output_dir.mkdir(exist_ok=True)
    manifest_path = output_dir / "manifest.csv"
    manifest_path.unlink(missing_ok=True)

    sources = {}
    for case_name in ("open", "closed"):
        path = CASE_DIR / f"mpm_{case_name}.json"
        if not path.is_file():
            parser.error(f"missing production source input: {path}")
        with path.open(encoding="utf-8") as stream:
            sources[case_name] = json.load(stream)

    rows = []
    for case_name, source in sources.items():
        analysis = source["analysis"]
        post = source["post_processing"]
        dt = float(analysis["dt"])
        total_steps = int(analysis["nsteps"])
        output_steps = int(post["output_steps"])
        segment_steps = round(args.segment_duration / dt)
        if not math.isclose(
            segment_steps * dt, args.segment_duration, abs_tol=1.0e-12
        ):
            parser.error(
                f"segment duration is not an integer multiple of {case_name} dt"
            )
        if total_steps % segment_steps:
            parser.error(
                f"{case_name} total steps {total_steps} are not divisible by "
                f"segment steps {segment_steps}"
            )
        if segment_steps % output_steps:
            parser.error(
                f"{case_name} segment boundary is not an output/checkpoint step"
            )

        number_of_segments = total_steps // segment_steps
        previous_uuid = "none"
        previous_end = 0
        for segment_index in range(number_of_segments):
            end_step = (segment_index + 1) * segment_steps
            uuid = f"siemens2013-production-{case_name}-seg{segment_index:03d}"
            config = clone(source)
            config["title"] = (
                f"Siemens 2013 {case_name} production segment "
                f"{segment_index + 1}/{number_of_segments}"
            )
            config["analysis"]["uuid"] = uuid
            config["analysis"]["nsteps"] = end_step
            config["post_processing"]["path"] = "results/"
            config["post_processing"]["write_hdf5"] = True
            config["post_processing"]["write_vtk"] = True
            config["post_processing"]["log_output_steps"] = False

            if segment_index == 0:
                config["analysis"]["resume"].update(
                    {
                        "resume": False,
                        "uuid": uuid,
                        "step": 0,
                        "nsteps": 0,
                        "start_from_this_step": False,
                        "this_step": 0,
                        "current_time": 0.0,
                    }
                )
                checkpoint = "none"
                resume_current_time = 0.0
            else:
                resume_current_time = (previous_end + 1) * dt
                config["analysis"]["resume"].update(
                    {
                        "resume": True,
                        "uuid": previous_uuid,
                        "step": previous_end,
                        "nsteps": previous_end,
                        "start_from_this_step": True,
                        "this_step": previous_end,
                        "current_time": resume_current_time,
                    }
                )
                checkpoint = (
                    f"results/{previous_uuid}/"
                    + output_name("particles", previous_end, previous_end, ".h5")
                )

            input_path = output_dir / f"mpm_{case_name}_seg{segment_index:03d}.json"
            write_json(input_path, config)
            result_dir = f"results/{uuid}"
            rows.append(
                {
                    "row_index": len(rows),
                    "case": case_name,
                    "segment_index": segment_index,
                    "segment_count": number_of_segments,
                    "start_step": 0 if segment_index == 0 else previous_end + 1,
                    "end_step": end_step,
                    "total_steps": total_steps,
                    "dt_s": dt,
                    "nominal_start_time_s": (
                        f"{segment_index * args.segment_duration:.9g}"
                    ),
                    "nominal_end_time_s": (
                        f"{(segment_index + 1) * args.segment_duration:.9g}"
                    ),
                    "resume_current_time_s": resume_current_time,
                    "uuid": uuid,
                    "input": input_path.relative_to(CASE_DIR).as_posix(),
                    "result_dir": result_dir,
                    "checkpoint": checkpoint,
                    "final_particle": (
                        f"{result_dir}/"
                        + output_name("particle", end_step, end_step, ".vtp")
                    ),
                    "final_checkpoint": (
                        f"{result_dir}/"
                        + output_name("particles", end_step, end_step, ".h5")
                    ),
                }
            )
            previous_uuid = uuid
            previous_end = end_step

    with manifest_path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(
        f"Wrote {len(rows)} production segments "
        f"({sum(row['case'] == 'open' for row in rows)} open, "
        f"{sum(row['case'] == 'closed' for row in rows)} closed) and "
        f"{manifest_path}"
    )


if __name__ == "__main__":
    main()
