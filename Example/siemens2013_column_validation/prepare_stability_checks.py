"""Generate short physical-time runs for selecting a stable explicit step."""

from __future__ import annotations

import csv
import json
from pathlib import Path


CASE_DIR = Path(__file__).resolve().parent
TIME_STEPS = (5.0e-4, 5.0e-5, 5.0e-6, 1.0e-6)
CHECK_DURATION_S = 3.0


def step_tag(dt: float) -> str:
    return f"{dt:.0e}".replace("+", "")


def main():
    input_dir = CASE_DIR / "stability_inputs"
    input_dir.mkdir(exist_ok=True)
    rows = []

    for case_name in ("open", "closed"):
        source_path = CASE_DIR / f"mpm_{case_name}.json"
        with source_path.open(encoding="utf-8") as stream:
            source = json.load(stream)

        for dt in TIME_STEPS:
            tag = step_tag(dt)
            label = f"{case_name}_{tag}"
            uuid = f"siemens2013-stability-r2-{label}"
            nsteps = round(CHECK_DURATION_S / dt)
            output_steps = max(1, nsteps // 20)

            config = json.loads(json.dumps(source))
            config["title"] = (
                f"Siemens 2013 {case_name} stability check, dt={dt:.0e} s"
            )
            config["analysis"]["dt"] = dt
            config["analysis"]["nsteps"] = nsteps
            config["analysis"]["uuid"] = uuid
            config["analysis"]["resume"].update(
                {"resume": False, "uuid": uuid, "step": 0, "nsteps": 0}
            )
            config["post_processing"]["path"] = "stability_results/"
            config["post_processing"]["output_steps"] = output_steps

            output_path = input_dir / f"mpm_{label}.json"
            with output_path.open("w", encoding="utf-8") as stream:
                json.dump(config, stream, indent=2)
                stream.write("\n")

            rows.append(
                {
                    "array_index": len(rows),
                    "label": label,
                    "case": case_name,
                    "dt_s": dt,
                    "duration_s": CHECK_DURATION_S,
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
    print(f"Wrote {len(rows)} stability inputs and {manifest}")


if __name__ == "__main__":
    main()
