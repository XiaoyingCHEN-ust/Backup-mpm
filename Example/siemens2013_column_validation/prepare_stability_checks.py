"""Generate short physical-time runs for selecting a stable explicit step."""

from __future__ import annotations

import argparse
import csv
import json
import re
from pathlib import Path


CASE_DIR = Path(__file__).resolve().parent
DEFAULT_TIME_STEPS = (5.0e-5, 5.0e-6, 1.0e-6)
DEFAULT_DURATION_S = 0.01
DEFAULT_REVISION = "r3-smoke"


def step_tag(dt: float) -> str:
    return f"{dt:.0e}".replace("+", "")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--duration", type=float, default=DEFAULT_DURATION_S)
    parser.add_argument(
        "--time-steps", type=float, nargs="+", default=DEFAULT_TIME_STEPS
    )
    parser.add_argument("--revision", default=DEFAULT_REVISION)
    args = parser.parse_args()

    if args.duration <= 0.0:
        parser.error("--duration must be positive")
    if any(dt <= 0.0 or dt > args.duration for dt in args.time_steps):
        parser.error("every time step must be positive and no larger than duration")
    if not re.fullmatch(r"[A-Za-z0-9_-]+", args.revision):
        parser.error("--revision may contain only letters, numbers, '_' and '-'")

    input_dir = CASE_DIR / "stability_inputs"
    input_dir.mkdir(exist_ok=True)
    rows = []

    for case_name in ("open", "closed"):
        source_path = CASE_DIR / f"mpm_{case_name}.json"
        with source_path.open(encoding="utf-8") as stream:
            source = json.load(stream)

        for dt in args.time_steps:
            tag = step_tag(dt)
            label = f"{case_name}_{tag}"
            uuid = f"siemens2013-stability-{args.revision}-{label}"
            nsteps = round(args.duration / dt)
            output_steps = max(1, nsteps // 20)

            config = json.loads(json.dumps(source))
            config["title"] = (
                f"Siemens 2013 {case_name} stability check, "
                f"dt={dt:.0e} s, duration={args.duration:g} s"
            )
            config["analysis"]["dt"] = dt
            config["analysis"]["nsteps"] = nsteps
            config["analysis"]["uuid"] = uuid
            config["analysis"]["resume"].update(
                {"resume": False, "uuid": uuid, "step": 0, "nsteps": 0}
            )
            config["post_processing"]["path"] = "stability_results/"
            config["post_processing"]["output_steps"] = output_steps

            output_path = input_dir / f"mpm_{args.revision}_{label}.json"
            with output_path.open("w", encoding="utf-8") as stream:
                json.dump(config, stream, indent=2)
                stream.write("\n")

            rows.append(
                {
                    "array_index": len(rows),
                    "label": label,
                    "case": case_name,
                    "dt_s": dt,
                    "duration_s": args.duration,
                    "nsteps": nsteps,
                    "output_steps": output_steps,
                    "uuid": uuid,
                    "input": output_path.relative_to(CASE_DIR).as_posix(),
                }
            )

    manifest = input_dir / "manifest.csv"
    with manifest.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(
        f"Wrote {len(rows)} stability inputs for revision {args.revision} "
        f"and {manifest}"
    )


if __name__ == "__main__":
    main()
