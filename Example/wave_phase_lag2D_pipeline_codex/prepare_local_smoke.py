#!/usr/bin/env python3
"""Generate exact, non-production local smoke configurations.

The smoke profile is deliberately separate from the registered screen and
production study.  Every generated JSON carries a hash of its canonical
payload and writes only below ``results/local_smoke/<label>/``.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
from pathlib import Path
from typing import Any

import prepare_study as study


CASE_ROOT = Path(__file__).resolve().parent
SMOKE_SCHEMA = "pipeline-local-smoke-v1"
SMOKE_STEPS = 5_000
SMOKE_OUTPUT_STEPS = 500
LABEL_PATTERN = re.compile(r"[a-z0-9][a-z0-9_-]{0,31}")

ROLE_FILENAMES = {
    "EQ_LS": "01_EQ_LS.json",
    "EQ_HS": "01_EQ_HS.json",
    "SANI_LS": "02_SANI_LS.json",
    "SANI_HS": "02_SANI_HS.json",
    "MC_EQ": "02_MC_EQ.json",
    "MC_DYNAMIC": "03_MC_DYNAMIC.json",
}
ROLE_DEPENDENCIES = {
    "EQ_LS": [],
    "EQ_HS": [],
    "SANI_LS": ["EQ_LS"],
    "SANI_HS": ["EQ_HS"],
    "MC_EQ": ["EQ_HS"],
    "MC_DYNAMIC": ["MC_EQ"],
}


def validate_label(label: str) -> str:
    """Return a filesystem-safe smoke label without silently normalising it."""

    if LABEL_PATTERN.fullmatch(label) is None:
        raise ValueError(
            "Smoke label must match [a-z0-9][a-z0-9_-]{0,31}; "
            f"got {label!r}"
        )
    return label


def uuid_prefix(label: str) -> str:
    return f"PLP_LOCAL_SMOKE_{validate_label(label).upper()}"


def result_path(label: str) -> str:
    return f"results/local_smoke/{validate_label(label)}/"


def _shorten_relaxation(config: dict[str, Any]) -> None:
    analysis = config["analysis"]
    analysis["nsteps"] = SMOKE_STEPS
    analysis["resume"]["step"] = SMOKE_STEPS
    analysis["resume"]["nsteps"] = SMOKE_STEPS
    analysis["rigid_pipeline"]["release_time"] = (SMOKE_STEPS + 1) * study.DT
    config["post_processing"]["output_steps"] = SMOKE_STEPS


def _payload(label: str, role: str) -> dict[str, Any]:
    """Build the exact unhashed payload for one allowed smoke role."""

    label = validate_label(label)
    if role not in ROLE_FILENAMES:
        raise ValueError(f"Unknown local smoke role: {role}")
    code_prefix = f"LOCAL_SMOKE_{label.upper()}"
    output = result_path(label)
    common = {
        "mesh_directory": ".",
        "cell_size": 0.02,
        "particle_spacing": 0.01,
        "wave_height": 0.12,
    }

    if role in {"EQ_LS", "EQ_HS"}:
        state = role.removeprefix("EQ_")
        saturation = (
            study.LOW_LAG_SATURATION
            if state == "LS"
            else study.HIGH_LAG_SATURATION
        )
        config = study.equilibrium_config(
            f"{code_prefix}_{state}",
            saturation,
            output,
            common["mesh_directory"],
            common["cell_size"],
            common["particle_spacing"],
            common["wave_height"],
        )
        _shorten_relaxation(config)
        return config

    if role == "MC_EQ":
        config = study.mc_handoff_config(
            code=f"{code_prefix}_HS",
            saturation=study.HIGH_LAG_SATURATION,
            equilibrium_uuid=f"{uuid_prefix(label)}_HS_EQ",
            result_path=output,
            **common,
        )
        _shorten_relaxation(config)
        return config

    if role in {"SANI_LS", "SANI_HS"}:
        state = role.removeprefix("SANI_")
        saturation = (
            study.LOW_LAG_SATURATION
            if state == "LS"
            else study.HIGH_LAG_SATURATION
        )
        config = study.dynamic_config(
            code=f"{code_prefix}_{state}_SANI",
            title=f"Local smoke {state} SANISAND dynamic",
            saturation=saturation,
            material_type="SANISAND2D",
            equilibrium_uuid=f"{uuid_prefix(label)}_{state}_EQ",
            result_path=output,
            nsteps=SMOKE_STEPS,
            physical_wave=True,
            equilibrium_steps=SMOKE_STEPS,
            **common,
        )
        config["post_processing"]["output_steps"] = SMOKE_OUTPUT_STEPS
        return config

    config = study.dynamic_config(
        code=f"{code_prefix}_HS_MC",
        title="Local smoke HS Mohr-Coulomb dynamic",
        saturation=study.HIGH_LAG_SATURATION,
        material_type="MohrCoulomb2D",
        equilibrium_uuid=f"{uuid_prefix(label)}_HS_MC_RELAX",
        result_path=output,
        nsteps=SMOKE_STEPS,
        physical_wave=True,
        equilibrium_steps=SMOKE_STEPS,
        **common,
    )
    config["post_processing"]["output_steps"] = SMOKE_OUTPUT_STEPS
    return config


def payload_sha256(config: dict[str, Any]) -> str:
    encoded = json.dumps(
        config, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def canonical_config(label: str, role: str) -> dict[str, Any]:
    """Return one hash-bound canonical smoke configuration."""

    config = _payload(label, role)
    config["local_smoke"] = {
        "schema": SMOKE_SCHEMA,
        "label": validate_label(label),
        "role": role,
        "baseline_sha256": payload_sha256(config),
        "result_eligibility": "local-smoke-only; excluded from manuscript evidence",
    }
    return config


def _atomic_write_json(path: Path, value: Any) -> None:
    if path.is_symlink():
        raise ValueError(f"Refusing to overwrite symlink: {path}")
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    try:
        with temporary.open("w", encoding="utf-8", newline="\n") as stream:
            json.dump(value, stream, indent=2, sort_keys=False)
            stream.write("\n")
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
    finally:
        if temporary.exists():
            temporary.unlink()


def generate(label: str, output_root: Path | None = None, *, force: bool = False) -> Path:
    """Write all smoke configs and return their manifest path."""

    label = validate_label(label)
    requested_root = (
        output_root
        if output_root is not None
        else CASE_ROOT / "configs" / "local_smoke" / label
    )
    root = requested_root.resolve()
    if output_root is None:
        try:
            root.relative_to(CASE_ROOT.resolve())
        except ValueError as error:
            raise ValueError(
                f"Smoke config path escapes the pipeline case directory: {requested_root}"
            ) from error
    if root.exists() and any(root.iterdir()) and not force:
        raise FileExistsError(
            f"Smoke config directory is not empty: {root}; choose a new label"
        )
    root.mkdir(parents=True, exist_ok=True)
    cases: list[dict[str, Any]] = []
    for role, filename in ROLE_FILENAMES.items():
        config = canonical_config(label, role)
        _atomic_write_json(root / filename, config)
        cases.append(
            {
                "role": role,
                "config": filename,
                "uuid": config["analysis"]["uuid"],
                "dependencies": ROLE_DEPENDENCIES[role],
                "baseline_sha256": config["local_smoke"]["baseline_sha256"],
            }
        )
    manifest = {
        "schema": "pipeline-local-smoke-manifest-v1",
        "label": label,
        "steps_per_stage": SMOKE_STEPS,
        "dt_s": study.DT,
        "result_path": result_path(label),
        "manuscript_evidence_eligible": False,
        "cases": cases,
    }
    manifest_path = root / "manifest.json"
    _atomic_write_json(manifest_path, manifest)
    return manifest_path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--label", required=True, help="unique lowercase run label")
    parser.add_argument(
        "--force",
        action="store_true",
        help="overwrite only the canonical config files for this exact label",
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()
    manifest = generate(args.label, force=args.force)
    print(f"Generated exact local smoke manifest: {manifest}")
    print(f"Results are isolated below: {result_path(args.label)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
