#!/usr/bin/env python3
"""Cross-check a particle HDF5 checkpoint against its same-step VTP."""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
from pathlib import Path
from typing import Any

import numpy as np

from analyze_study import read_vtp
import qa_contract as qa


SCRIPT_DIR = Path(__file__).resolve().parent
REPOSITORY_ROOT = SCRIPT_DIR.parents[1]
DEFAULT_AUDITOR = (
    REPOSITORY_ROOT
    / "mpm_hpc_source"
    / "build-pipeline"
    / "hdf5_checkpoint_audit"
)
SCHEMA = qa.HDF5_VTP_AUDIT_SCHEMA
COMPARED_FIELDS = qa.HDF5_VTP_COMPARED_FIELDS
# HDF5 stores doubles while VTK may serialize Float32. These absolute
# tolerances are intentionally tight enough to detect a wrong or stale file.
TOLERANCES = qa.HDF5_VTP_ABSOLUTE_TOLERANCES


def id_sha256(identifiers: np.ndarray) -> str:
    values = np.asarray(identifiers, dtype="<u8").reshape(-1)
    return hashlib.sha256(values.tobytes(order="C")).hexdigest()


def _as_matrix(value: Any, width: int, name: str) -> np.ndarray:
    matrix = np.asarray(value, dtype=np.float64)
    if matrix.shape != (width,) or not np.all(np.isfinite(matrix)):
        raise ValueError(f"HDF5 JSONL particle has malformed {name}")
    return matrix


def read_hdf5_jsonl(
    hdf5: Path, auditor: Path
) -> tuple[dict[str, Any], dict[str, np.ndarray]]:
    completed = subprocess.run(
        [str(auditor), "--hdf5", str(hdf5), "--jsonl"],
        capture_output=True,
        text=True,
        check=False,
        timeout=120,
    )
    if completed.returncode != 0:
        raise ValueError(
            f"HDF5 structural audit failed with code {completed.returncode}: "
            f"{completed.stderr.strip()}"
        )
    lines = completed.stdout.splitlines()
    if not lines:
        raise ValueError("HDF5 auditor produced no output")
    summary = json.loads(lines[0])
    if summary.get("schema") != "mpm-hdf5-checkpoint-audit-v1" or not summary.get(
        "passed"
    ):
        raise ValueError("HDF5 auditor summary did not pass")
    count = int(summary["particle_count"])
    records = [json.loads(line) for line in lines[1:]]
    if len(records) != count:
        raise ValueError("HDF5 JSONL record count differs from summary")
    ids = np.asarray([record.get("id") for record in records], dtype=np.int64)
    if not np.array_equal(ids, np.arange(count, dtype=np.int64)):
        raise ValueError("HDF5 JSONL particle IDs are incomplete or reordered")
    return summary, {
        "ids": ids,
        "coordinates": np.vstack(
            [_as_matrix(record.get("coordinates"), 3, "coordinates") for record in records]
        ),
        "velocities": np.vstack(
            [_as_matrix(record.get("velocities"), 3, "velocities") for record in records]
        ),
        "displacements": np.vstack(
            [
                _as_matrix(record.get("displacements"), 3, "displacements")
                for record in records
            ]
        ),
        "stresses": np.vstack(
            [_as_matrix(record.get("stresses"), 6, "stresses") for record in records]
        ),
        "porosities": np.asarray(
            [record.get("porosity") for record in records], dtype=np.float64
        ).reshape(count, 1),
        "volumes": np.asarray(
            [record.get("volume") for record in records], dtype=np.float64
        ).reshape(count, 1),
    }


def _vtp_ids(arrays: dict[str, np.ndarray], count: int) -> np.ndarray:
    if "ids" not in arrays:
        raise ValueError("VTP lacks particle IDs")
    raw = np.asarray(arrays["ids"]).reshape(-1)
    if raw.size != count or not np.all(np.isfinite(raw)):
        raise ValueError("VTP particle IDs are incomplete or non-finite")
    identifiers = raw.astype(np.int64)
    if not np.array_equal(raw, identifiers):
        raise ValueError("VTP particle IDs are not integral")
    if not np.array_equal(np.sort(identifiers), np.arange(count, dtype=np.int64)):
        raise ValueError("VTP particle IDs are missing or duplicated")
    return identifiers


def crosscheck(
    hdf5: Path, vtp: Path, auditor: Path = DEFAULT_AUDITOR
) -> dict[str, Any]:
    hdf5 = hdf5.resolve()
    vtp = vtp.resolve()
    auditor = auditor.resolve()
    for description, path in (
        ("HDF5 checkpoint", hdf5),
        ("VTP checkpoint", vtp),
        ("HDF5 auditor", auditor),
    ):
        if not path.is_file():
            raise FileNotFoundError(f"{description} is missing: {path}")

    hdf5_summary, hdf5_fields = read_hdf5_jsonl(hdf5, auditor)
    points, arrays = read_vtp(vtp)
    count_vtp = len(points)
    ids_vtp = _vtp_ids(arrays, count_vtp)
    count_hdf5 = len(hdf5_fields["ids"])
    if count_hdf5 != count_vtp:
        raise ValueError(
            f"particle counts differ: HDF5={count_hdf5}, VTP={count_vtp}"
        )

    order = np.argsort(ids_vtp)
    vtp_fields: dict[str, np.ndarray] = {"coordinates": np.asarray(points)[order]}
    for name in ("velocities", "displacements", "stresses"):
        if name not in arrays:
            raise ValueError(f"VTP lacks required field {name}")
        vtp_fields[name] = np.asarray(arrays[name], dtype=np.float64)[order]
    for name in ("porosities", "volumes"):
        if name not in arrays:
            raise ValueError(f"VTP lacks required field {name}")
        vtp_fields[name] = np.asarray(arrays[name], dtype=np.float64).reshape(
            count_vtp, -1
        )[order]

    maximum_differences: dict[str, float] = {}
    for name in COMPARED_FIELDS:
        observed = vtp_fields[name]
        reference = hdf5_fields[name]
        if observed.shape != reference.shape:
            raise ValueError(
                f"HDF5/VTP {name} shapes differ: {reference.shape} vs {observed.shape}"
            )
        if not np.all(np.isfinite(observed)):
            raise ValueError(f"VTP {name} contains non-finite values")
        difference = float(np.max(np.abs(observed - reference)))
        maximum_differences[name] = difference
        if difference > TOLERANCES[name]:
            raise ValueError(
                f"HDF5/VTP {name} difference {difference:.6e} exceeds "
                f"{TOLERANCES[name]:.6e}"
            )

    digest_hdf5 = id_sha256(hdf5_fields["ids"])
    digest_vtp = id_sha256(ids_vtp[order])
    if digest_hdf5 != digest_vtp:
        raise ValueError("HDF5/VTP particle ID hashes differ")
    result = {
        "schema": SCHEMA,
        "passed": True,
        "hdf5": str(hdf5),
        "vtp": str(vtp),
        "particle_count_hdf5": count_hdf5,
        "particle_count_vtp": count_vtp,
        "ids_sha256_hdf5": digest_hdf5,
        "ids_sha256_vtp": digest_vtp,
        "compared_fields": list(COMPARED_FIELDS),
        "maximum_absolute_differences": maximum_differences,
        "tolerances": dict(TOLERANCES),
        "hdf5_structure": hdf5_summary,
    }
    return qa.validate_hdf5_vtp_audit(
        result, expected_hdf5=str(hdf5), expected_vtp=str(vtp)
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--hdf5", type=Path, required=True)
    parser.add_argument("--vtp", type=Path, required=True)
    parser.add_argument("--auditor", type=Path, default=DEFAULT_AUDITOR)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = crosscheck(args.hdf5, args.vtp, args.auditor)
    text = json.dumps(result, sort_keys=True)
    print(text)
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
