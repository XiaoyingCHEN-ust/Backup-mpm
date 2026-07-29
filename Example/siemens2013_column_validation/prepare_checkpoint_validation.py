"""Generate a two-segment restart test against the continuous 3 s open run."""

from __future__ import annotations

import json
import math
from pathlib import Path


CASE_DIR = Path(__file__).resolve().parent
DT = 2.5e-6
SPLIT_TIME_S = 1.5
TOTAL_TIME_S = 3.0
SEGMENT_A_UUID = "siemens2013-r20-checkpoint-open-segment-a"
SEGMENT_B_UUID = "siemens2013-r20-checkpoint-open-segment-b"


def clone(source):
    return json.loads(json.dumps(source))


def configure_common(config, uuid: str) -> None:
    analysis = config["analysis"]
    analysis["PIC"] = 0.0
    analysis["PIC_T"] = 0.0
    analysis["pressure_smoothing"] = False
    analysis["dt"] = DT
    analysis["uuid"] = uuid
    config["post_processing"]["path"] = "stability_results/"
    config["post_processing"]["write_hdf5"] = True
    config["post_processing"]["write_vtk"] = True
    config["post_processing"]["log_output_steps"] = False


def write_json(path: Path, config) -> None:
    with path.open("w", encoding="utf-8") as stream:
        json.dump(config, stream, indent=2)
        stream.write("\n")


def main() -> None:
    source_path = CASE_DIR / "mpm_open.json"
    if not source_path.is_file():
        raise FileNotFoundError(f"Missing open-column source input: {source_path}")
    with source_path.open(encoding="utf-8") as stream:
        source = json.load(stream)

    input_dir = CASE_DIR / "stability_inputs"
    input_dir.mkdir(exist_ok=True)

    split_step = round(SPLIT_TIME_S / DT)
    total_step = round(TOTAL_TIME_S / DT)
    if (
        not math.isclose(split_step * DT, SPLIT_TIME_S, abs_tol=1.0e-12)
        or not math.isclose(total_step * DT, TOTAL_TIME_S, abs_tol=1.0e-12)
    ):
        raise RuntimeError("Checkpoint times must be exact multiples of dt")

    # The solver updates once at step 0. The HDF5 state labelled split_step
    # has therefore completed split_step + 1 updates, exactly matching the
    # continuous reference at the same global step.
    checkpoint_time = (split_step + 1) * DT
    resume_start_step = split_step + 1
    continuous_updates = total_step + 1
    split_updates = (split_step + 1) + (total_step - resume_start_step + 1)
    if split_updates != continuous_updates:
        raise RuntimeError("Split and continuous runs have different update counts")

    segment_a = clone(source)
    configure_common(segment_a, SEGMENT_A_UUID)
    segment_a["title"] = "Siemens r20 open checkpoint segment A (0-1.5 s)"
    segment_a["analysis"]["nsteps"] = split_step
    segment_a["analysis"]["resume"].update(
        {
            "resume": False,
            "uuid": SEGMENT_A_UUID,
            "step": 0,
            "nsteps": 0,
            "start_from_this_step": False,
            "this_step": 0,
            "current_time": 0.0,
        }
    )
    segment_a["post_processing"]["output_steps"] = split_step

    segment_b = clone(source)
    configure_common(segment_b, SEGMENT_B_UUID)
    segment_b["title"] = "Siemens r20 open checkpoint segment B (1.5-3.0 s)"
    segment_b["analysis"]["nsteps"] = total_step
    segment_b["analysis"]["resume"].update(
        {
            "resume": True,
            "uuid": SEGMENT_A_UUID,
            "step": split_step,
            "nsteps": split_step,
            "start_from_this_step": True,
            "this_step": split_step,
            "current_time": checkpoint_time,
        }
    )
    segment_b["post_processing"]["output_steps"] = total_step // 20

    segment_a_path = input_dir / "mpm_r20_checkpoint_open_segment_a.json"
    segment_b_path = input_dir / "mpm_r20_checkpoint_open_segment_b.json"
    write_json(segment_a_path, segment_a)
    write_json(segment_b_path, segment_b)

    manifest = {
        "dt_s": DT,
        "split_step": split_step,
        "checkpoint_time_s": checkpoint_time,
        "resume_start_step": resume_start_step,
        "total_step": total_step,
        "total_updates": continuous_updates,
        "continuous_reference_uuid": (
            "siemens2013-stability-r19b-open-3s-open_2p5e-06"
        ),
        "segments": [
            {
                "name": "a",
                "uuid": SEGMENT_A_UUID,
                "input": segment_a_path.relative_to(CASE_DIR).as_posix(),
            },
            {
                "name": "b",
                "uuid": SEGMENT_B_UUID,
                "input": segment_b_path.relative_to(CASE_DIR).as_posix(),
                "resume_uuid": SEGMENT_A_UUID,
                "resume_step": split_step,
            },
        ],
    }
    manifest_path = input_dir / "r20_checkpoint_manifest.json"
    write_json(manifest_path, manifest)
    print(f"Wrote r20 checkpoint inputs and {manifest_path}")


if __name__ == "__main__":
    main()
