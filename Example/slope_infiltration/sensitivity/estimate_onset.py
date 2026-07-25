#!/usr/bin/env python3
"""Estimate instability time from the p99 displacement history."""

from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path


DISPLACEMENT_THRESHOLDS_M = (0.15, 0.20, 0.25)


def interpolate_crossing(
    samples: list[dict[str, str]], threshold_increment: float
) -> float | None:
    samples = sorted(samples, key=lambda row: float(row["time_s"]))
    initial = float(samples[0]["p99_displacement_m"])
    target = initial + threshold_increment

    previous_time = float(samples[0]["time_s"])
    previous_value = float(samples[0]["p99_displacement_m"])
    if previous_value >= target:
        return previous_time

    for row in samples[1:]:
        time = float(row["time_s"])
        value = float(row["p99_displacement_m"])
        if value >= target:
            if value == previous_value:
                return time
            fraction = (target - previous_value) / (value - previous_value)
            return previous_time + fraction * (time - previous_time)
        previous_time = time
        previous_value = value
    return None


def main() -> None:
    script_dir = Path(__file__).resolve().parent
    with (script_dir / "metrics.csv").open(
        newline="", encoding="utf-8"
    ) as stream:
        metrics = list(csv.DictReader(stream))
    with (script_dir / "manifest.csv").open(
        newline="", encoding="utf-8"
    ) as stream:
        manifest = {
            row["case_id"]: row for row in csv.DictReader(stream)
        }

    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in metrics:
        grouped[row["case_id"]].append(row)

    rows: list[dict[str, str | float]] = []
    for case_id in sorted(grouped):
        estimates = {
            threshold: interpolate_crossing(grouped[case_id], threshold)
            for threshold in DISPLACEMENT_THRESHOLDS_M
        }
        case = manifest[case_id]
        rows.append(
            {
                **case,
                "onset_s": (
                    "" if estimates[0.20] is None else f"{estimates[0.20]:.6g}"
                ),
                "onset_at_0p15m_s": (
                    "" if estimates[0.15] is None else f"{estimates[0.15]:.6g}"
                ),
                "onset_at_0p25m_s": (
                    "" if estimates[0.25] is None else f"{estimates[0.25]:.6g}"
                ),
            }
        )

    output_path = script_dir / "onset_times.csv"
    with output_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {len(rows)} onset estimates to {output_path}")


if __name__ == "__main__":
    main()
