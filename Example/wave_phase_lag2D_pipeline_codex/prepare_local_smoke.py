#!/usr/bin/env python3
"""Generate exact, non-production local smoke configurations.

The smoke profile is deliberately separate from the registered screen and
production study.  Its default inventory is the four required EQ/SANISAND
stages; MC_EQ is an explicit diagnostic and no MC dynamic smoke exists. Every
generated JSON carries a hash of its canonical payload and writes only below
``results/local_smoke/<label>/``.
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
MANIFEST_SCHEMA = "pipeline-local-smoke-manifest-v2"
SMOKE_STEPS = 5_000
SMOKE_OUTPUT_STEPS = 500
LABEL_PATTERN = re.compile(r"[a-z0-9][a-z0-9_-]{0,31}")

REQUIRED_ROLES = ("EQ_LS", "EQ_HS", "SANI_LS", "SANI_HS")
DIAGNOSTIC_ROLES = ("MC_EQ",)
ROLE_FILENAMES: dict[str, str] = {
    "EQ_LS": "01_EQ_LS.json",
    "EQ_HS": "01_EQ_HS.json",
    "SANI_LS": "02_SANI_LS.json",
    "SANI_HS": "02_SANI_HS.json",
    "MC_EQ": "02_MC_EQ.json",
}
ROLE_DEPENDENCIES = {
    "EQ_LS": [],
    "EQ_HS": [],
    "SANI_LS": ["EQ_LS"],
    "SANI_HS": ["EQ_HS"],
    "MC_EQ": ["EQ_HS"],
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

    raise AssertionError(f"Unhandled local smoke role: {role}")


def payload_sha256(config: dict[str, Any]) -> str:
    encoded = json.dumps(
        config, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def manifest_payload_sha256(manifest: dict[str, Any]) -> str:
    payload = dict(manifest)
    payload.pop("manifest_sha256", None)
    return payload_sha256(payload)


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


def validate_manifest(manifest_path: Path) -> dict[str, Any]:
    """Require an exact hash-bound required/diagnostic smoke manifest."""

    manifest_path = manifest_path.resolve()
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if set(manifest) != {
        "schema",
        "label",
        "steps_per_stage",
        "dt_s",
        "result_path",
        "manuscript_evidence_eligible",
        "required_roles",
        "diagnostic_roles",
        "diagnostic_policy",
        "cases",
        "manifest_sha256",
    }:
        raise ValueError("Local smoke manifest fields are missing or unexpected")
    if manifest.get("schema") != MANIFEST_SCHEMA:
        raise ValueError("Local smoke manifest schema is stale")
    label = validate_label(str(manifest.get("label", "")))
    if manifest.get("manifest_sha256") != manifest_payload_sha256(manifest):
        raise ValueError("Local smoke manifest hash is stale")
    if (
        manifest.get("steps_per_stage") != SMOKE_STEPS
        or manifest.get("dt_s") != study.DT
        or manifest.get("result_path") != result_path(label)
        or manifest.get("manuscript_evidence_eligible") is not False
        or manifest.get("required_roles") != list(REQUIRED_ROLES)
        or manifest.get("diagnostic_roles") not in ([], list(DIAGNOSTIC_ROLES))
        or manifest.get("diagnostic_policy")
        != "explicit-opt-in; fail-fast; excluded from manuscript evidence"
    ):
        raise ValueError("Local smoke manifest contract is stale")

    selected_roles = list(REQUIRED_ROLES) + list(manifest["diagnostic_roles"])
    cases = manifest.get("cases")
    if not isinstance(cases, list) or len(cases) != len(selected_roles):
        raise ValueError("Local smoke manifest case inventory is incomplete")
    for entry, role in zip(cases, selected_roles):
        if not isinstance(entry, dict) or set(entry) != {
            "role",
            "classification",
            "required_for_smoke_pass",
            "config",
            "uuid",
            "dependencies",
            "baseline_sha256",
            "config_sha256",
        }:
            raise ValueError("Local smoke manifest case fields are malformed")
        required = role in REQUIRED_ROLES
        expected_classification = "required" if required else "diagnostic"
        config_path = manifest_path.parent / ROLE_FILENAMES[role]
        if (
            entry.get("role") != role
            or entry.get("classification") != expected_classification
            or entry.get("required_for_smoke_pass") is not required
            or entry.get("config") != ROLE_FILENAMES[role]
            or entry.get("dependencies") != ROLE_DEPENDENCIES[role]
        ):
            raise ValueError(f"Local smoke manifest role {role} is stale")
        if not config_path.is_file():
            raise FileNotFoundError(f"Local smoke config is missing: {config_path}")
        config = json.loads(config_path.read_text(encoding="utf-8"))
        canonical = canonical_config(label, role)
        if config != canonical:
            raise ValueError(f"Local smoke config {role} is not canonical")
        if (
            entry.get("uuid") != canonical["analysis"]["uuid"]
            or entry.get("baseline_sha256")
            != canonical["local_smoke"]["baseline_sha256"]
            or entry.get("config_sha256") != file_sha256(config_path)
        ):
            raise ValueError(f"Local smoke manifest hash for {role} is stale")
    return manifest


def generate(
    label: str,
    output_root: Path | None = None,
    *,
    force: bool = False,
    include_mc_diagnostic: bool = False,
) -> Path:
    """Write required smoke configs and optional MC_EQ diagnostic config."""

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
    selected_roles = list(REQUIRED_ROLES)
    if include_mc_diagnostic:
        selected_roles.extend(DIAGNOSTIC_ROLES)
    expected_names = {ROLE_FILENAMES[role] for role in selected_roles} | {
        "manifest.json"
    }
    if root.exists() and any(root.iterdir()):
        if not force:
            raise FileExistsError(
                f"Smoke config directory is not empty: {root}; choose a new label"
            )
        unexpected = sorted(
            path.name for path in root.iterdir() if path.name not in expected_names
        )
        if unexpected:
            raise ValueError(
                "Refusing to reuse a smoke label containing unexpected files: "
                + ", ".join(unexpected)
            )
    root.mkdir(parents=True, exist_ok=True)
    cases: list[dict[str, Any]] = []
    for role in selected_roles:
        filename = ROLE_FILENAMES[role]
        config = canonical_config(label, role)
        config_path = root / filename
        _atomic_write_json(config_path, config)
        required = role in REQUIRED_ROLES
        cases.append(
            {
                "role": role,
                "classification": "required" if required else "diagnostic",
                "required_for_smoke_pass": required,
                "config": filename,
                "uuid": config["analysis"]["uuid"],
                "dependencies": ROLE_DEPENDENCIES[role],
                "baseline_sha256": config["local_smoke"]["baseline_sha256"],
                "config_sha256": file_sha256(config_path),
            }
        )
    manifest = {
        "schema": MANIFEST_SCHEMA,
        "label": label,
        "steps_per_stage": SMOKE_STEPS,
        "dt_s": study.DT,
        "result_path": result_path(label),
        "manuscript_evidence_eligible": False,
        "required_roles": list(REQUIRED_ROLES),
        "diagnostic_roles": (
            list(DIAGNOSTIC_ROLES) if include_mc_diagnostic else []
        ),
        "diagnostic_policy": (
            "explicit-opt-in; fail-fast; excluded from manuscript evidence"
        ),
        "cases": cases,
    }
    manifest["manifest_sha256"] = manifest_payload_sha256(manifest)
    manifest_path = root / "manifest.json"
    _atomic_write_json(manifest_path, manifest)
    validate_manifest(manifest_path)
    return manifest_path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    action = parser.add_mutually_exclusive_group(required=True)
    action.add_argument("--label", help="unique lowercase run label")
    action.add_argument(
        "--validate-manifest",
        type=Path,
        help="revalidate an existing manifest and all hash-bound configs",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="overwrite only the canonical config files for this exact label",
    )
    parser.add_argument(
        "--include-mc-diagnostic",
        action="store_true",
        help=(
            "also generate the fail-fast MC_EQ diagnostic; this never creates "
            "an MC dynamic case and is never manuscript evidence"
        ),
    )
    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()
    if args.validate_manifest is not None:
        if args.force or args.include_mc_diagnostic:
            parser.error(
                "--force/--include-mc-diagnostic apply only when generating"
            )
        validate_manifest(args.validate_manifest)
        print(f"Validated exact local smoke manifest: {args.validate_manifest}")
        return 0
    manifest = generate(
        str(args.label),
        force=args.force,
        include_mc_diagnostic=args.include_mc_diagnostic,
    )
    print(f"Generated exact local smoke manifest: {manifest}")
    print(f"Results are isolated below: {result_path(args.label)}")
    if args.include_mc_diagnostic:
        print(
            "MC_EQ diagnostic enabled: any unchanged stability-gate failure "
            "terminates this command nonzero"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
