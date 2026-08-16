#!/usr/bin/env python3
"""Generate an isolated low-permeability/low-saturation phase-lag probe.

This generator deliberately does not alter or impersonate the registered
screen/production study.  It clones the audited HS equilibrium and physical
wave configurations into an exploratory result namespace, retains the strict
LinearElastic2D/PIC equilibrium protocol, and keeps the pipeline fixed during
the SANISAND wave probe so the pressure response can be interpreted without
pipeline-motion feedback.
"""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
import math
import re
from pathlib import Path
from typing import Any

import numpy as np


SCHEMA = "pipeline-phase-lag-exploration-v1"
DEFAULT_LABEL = "k3e-12_sw090"
DEFAULT_PERMEABILITY_M2 = 3.0e-12
DEFAULT_SATURATION = 0.90
DEFAULT_EQ_STEPS = 5000
DEFAULT_WAVE_STEPS = 39000
OUTPUT_INTERVAL = 650
PRESSURE_INTERVAL = 130
DT = 1.0e-4
GRAVITY = 9.81
K0_EFFECTIVE = 0.72


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def save_ascii_table(path: Path, values: np.ndarray) -> None:
    with path.open("w", encoding="utf-8") as stream:
        stream.write(f"{values.shape[0]}\n")
        np.savetxt(stream, values, fmt="%.10f", delimiter="\t")


def initial_effective_stresses(
    particles: np.ndarray, equilibrium: dict[str, Any], saturation: float
) -> tuple[np.ndarray, float]:
    soil, fluid = equilibrium["materials"]
    porosity = float(soil["porosity"])
    gas_density = (
        float(fluid["gas_molar_mass"])
        * float(soil["p_ref"])
        / float(fluid["gas_constant"])
        / (273.15 + float(soil["initial_temperature"]))
    )
    mixture_density = (1.0 - porosity) * float(soil["density"]) + porosity * (
        saturation * float(fluid["density"]) + (1.0 - saturation) * gas_density
    )
    effective_unit_weight = (mixture_density - float(fluid["density"])) * GRAVITY
    if not math.isfinite(effective_unit_weight) or effective_unit_weight <= 0.0:
        raise ValueError("Exploratory effective unit weight must be positive")
    seabed = float(fluid["sea_level"]) - float(fluid["depth_left"])
    burial_depth = np.maximum(seabed - particles[:, 1], 0.0)
    sigma_yy = -effective_unit_weight * burial_depth
    stresses = np.zeros((particles.shape[0], 6), dtype=np.float64)
    stresses[:, 0] = K0_EFFECTIVE * sigma_yy
    stresses[:, 1] = sigma_yy
    stresses[:, 2] = K0_EFFECTIVE * sigma_yy
    return stresses, effective_unit_weight


def set_hydraulic_parameters(
    config: dict[str, Any], permeability_m2: float, saturation: float
) -> None:
    soil, fluid = config["materials"]
    soil["intrinsic_permeability"] = permeability_m2
    fluid["liquid_saturation"] = saturation
    fluid["gas_saturation"] = 1.0 - saturation


def generate(
    root: Path,
    *,
    label: str,
    permeability_m2: float,
    saturation: float,
    equilibrium_steps: int,
    wave_steps: int,
    pressure_smoothing: bool,
    paired_replay: bool,
) -> dict[str, Any]:
    if re.fullmatch(r"[a-z0-9][a-z0-9_-]*", label) is None:
        raise ValueError("label must match [a-z0-9][a-z0-9_-]*")
    if not math.isfinite(permeability_m2) or permeability_m2 <= 0.0:
        raise ValueError("permeability must be finite and positive")
    if not math.isfinite(saturation) or not 0.001 < saturation < 0.999:
        raise ValueError("saturation must lie strictly between 0.001 and 0.999")
    if equilibrium_steps <= 0 or wave_steps < DEFAULT_WAVE_STEPS:
        raise ValueError("equilibrium steps must be positive and wave steps >= 39000")
    if wave_steps % OUTPUT_INTERVAL:
        raise ValueError("wave steps must be divisible by the 650-step output interval")

    baseline_eq_path = root / "configs/screen/01_EQ_HS.json"
    baseline_wave_path = root / "configs/screen/02_HS.json"
    equilibrium = copy.deepcopy(read_json(baseline_eq_path))
    wave = copy.deepcopy(read_json(baseline_wave_path))
    particles_path = root / equilibrium["particles"][0]["generator"]["location"]
    particles = np.loadtxt(particles_path, skiprows=1, dtype=np.float64)
    if particles.ndim != 2 or particles.shape[1] != 2:
        raise ValueError(f"Invalid particle coordinate table: {particles_path}")

    config_dir = root / "configs/phase_lag_exploratory" / label
    config_dir.mkdir(parents=True, exist_ok=True)
    result_path = f"results/phase_lag_exploratory/{label}/"
    label_token = label.upper().replace("-", "_")
    eq_uuid = f"PLP_EXP_{label_token}_EQ"
    wave_uuid = f"PLP_EXP_{label_token}_WAVE"
    stress_path = config_dir / "initial_effective_stresses.txt"

    set_hydraulic_parameters(equilibrium, permeability_m2, saturation)
    set_hydraulic_parameters(wave, permeability_m2, saturation)
    for config in (equilibrium, wave):
        config["analysis"]["pressure_smoothing"] = pressure_smoothing
        # No state write-back is allowed in this exploratory family.  When
        # pressure_smoothing is false the solver also bypasses the VTK-time
        # particle pressure-smoothing path and uses the unsmoothed PIC field.
        config["analysis"]["pressure_smoothing_in_loop"] = False
    stresses, unit_weight = initial_effective_stresses(
        particles, equilibrium, saturation
    )
    save_ascii_table(stress_path, stresses)

    equilibrium["title"] = (
        f"Exploratory phase-lag equilibrium: k={permeability_m2:.6g}, "
        f"Sw={saturation:.6g}"
    )
    equilibrium["mesh"]["particles_stresses"] = str(stress_path.relative_to(root))
    equilibrium["analysis"]["uuid"] = eq_uuid
    equilibrium["analysis"]["nsteps"] = equilibrium_steps
    equilibrium["analysis"]["rigid_pipeline"]["fixed"] = True
    equilibrium["analysis"]["rigid_pipeline"]["release_time"] = (
        equilibrium_steps + 1
    ) * DT
    equilibrium["analysis"]["resume"].update(
        {
            "resume": False,
            "uuid": eq_uuid,
            "step": equilibrium_steps,
            "nsteps": equilibrium_steps,
            "start_from_this_step": False,
            "this_step": 0,
            "current_time": 0.0,
        }
    )
    equilibrium["post_processing"]["path"] = result_path
    equilibrium["post_processing"]["output_steps"] = equilibrium_steps

    wave["title"] = (
        f"Exploratory fixed-pipeline phase-lag wave: k={permeability_m2:.6g}, "
        f"Sw={saturation:.6g}"
    )
    wave["analysis"]["uuid"] = wave_uuid
    wave["analysis"]["nsteps"] = wave_steps
    wave["analysis"]["rigid_pipeline"]["fixed"] = True
    wave["analysis"]["rigid_pipeline"]["release_time"] = (wave_steps + 1) * DT
    wave["analysis"]["resume"].update(
        {
            "resume": True,
            "uuid": eq_uuid,
            "step": equilibrium_steps,
            "nsteps": equilibrium_steps,
            "start_from_this_step": False,
            "this_step": 0,
            "current_time": 0.0,
        }
    )
    wave["post_processing"]["path"] = result_path
    wave["post_processing"]["output_steps"] = OUTPUT_INTERVAL

    eq_path = config_dir / "01_EQ.json"
    wave_path = config_dir / "02_WAVE.json"
    eq_path.write_text(json.dumps(equilibrium, indent=2) + "\n", encoding="utf-8")
    wave_path.write_text(json.dumps(wave, indent=2) + "\n", encoding="utf-8")
    generated: dict[str, Any] = {
        "equilibrium_config": str(eq_path.relative_to(root)),
        "equilibrium_sha256": sha256(eq_path),
        "wave_config": str(wave_path.relative_to(root)),
        "wave_sha256": sha256(wave_path),
        "initial_effective_stresses": str(stress_path.relative_to(root)),
        "initial_effective_stresses_sha256": sha256(stress_path),
        "equilibrium_uuid": eq_uuid,
        "wave_uuid": wave_uuid,
        "result_path": result_path,
    }

    if paired_replay:
        pressure_root = f"pressure_databases/phase_lag_exploratory/{label}"

        def replay_case(
            role: str,
            *,
            wave_pressure: bool,
            pressure_mode: str,
            pressure_leaf: str,
            pressure_prefix: str,
        ) -> tuple[dict[str, Any], Path]:
            config = copy.deepcopy(wave)
            config["title"] = (
                f"Exploratory {role} phase-lag replay: "
                f"k={permeability_m2:.6g}, Sw={saturation:.6g}, "
                f"smoothing={pressure_smoothing}"
            )
            config["analysis"]["uuid"] = f"PLP_EXP_{label_token}_{role}"
            config["analysis"]["rigid_pipeline"]["fixed"] = True
            config["analysis"]["rigid_pipeline"]["release_time"] = (
                wave_steps + 1
            ) * DT
            config["materials"][1]["wave_pressure"] = wave_pressure
            config["analysis"]["prescribed_phase_pressures"] = {
                "enable": pressure_mode == "read",
                "write": pressure_mode == "write",
                "path": f"{pressure_root}/{pressure_leaf}",
                "file_prefix": pressure_prefix,
                "source_dt": DT,
                "step_interval": PRESSURE_INTERVAL,
                "max_step": wave_steps,
                "mapping": "id",
            }
            filename = {
                "HD": "03_HD.json",
                "RL": "04_RL.json",
                "RE": "05_RE.json",
            }[role]
            return config, config_dir / filename

        replay_specs = {
            "HD": replay_case(
                "HD",
                wave_pressure=True,
                pressure_mode="write",
                pressure_leaf="lagged",
                pressure_prefix="pressure",
            ),
            "RL": replay_case(
                "RL",
                wave_pressure=False,
                pressure_mode="read",
                pressure_leaf="lagged",
                pressure_prefix="pressure",
            ),
            "RE": replay_case(
                "RE",
                wave_pressure=False,
                pressure_mode="read",
                pressure_leaf="phase_erased",
                pressure_prefix="phase_erased",
            ),
        }
        for role, (config, path) in replay_specs.items():
            path.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
            generated[f"{role.lower()}_config"] = str(path.relative_to(root))
            generated[f"{role.lower()}_sha256"] = sha256(path)
            generated[f"{role.lower()}_uuid"] = config["analysis"]["uuid"]
        generated["lagged_pressure_path"] = f"{pressure_root}/lagged"
        generated["phase_erased_pressure_path"] = (
            f"{pressure_root}/phase_erased"
        )

    manifest = {
        "schema": SCHEMA,
        "status": "exploratory-only-not-registered-manuscript-evidence",
        "label": label,
        "parameters": {
            "intrinsic_permeability_m2": permeability_m2,
            "liquid_saturation": saturation,
            "gas_saturation": 1.0 - saturation,
            "equilibrium_steps": equilibrium_steps,
            "wave_steps": wave_steps,
            "dt_s": DT,
            "wave_output_interval_steps": OUTPUT_INTERVAL,
            "pressure_smoothing": pressure_smoothing,
            "paired_replay": paired_replay,
            "pressure_database_interval_steps": (
                PRESSURE_INTERVAL if paired_replay else None
            ),
        },
        "guardrails": {
            "equilibrium_material": "LinearElastic2D",
            "equilibrium_pic": 1.0,
            "equilibrium_apic": False,
            "equilibrium_damping_s-1": 5.0,
            "maximum_velocity_m_s": 1.0e-3,
            "wave_material": "SANISAND2D",
            "wave_pipeline_fixed": True,
            "wave_damping_s-1": 0.0,
            "pressure_smoothing": pressure_smoothing,
            "pressure_smoothing_in_loop": False,
        },
        "effective_unit_weight_n_m3": unit_weight,
        "baseline": {
            "equilibrium_config": str(baseline_eq_path.relative_to(root)),
            "equilibrium_sha256": sha256(baseline_eq_path),
            "wave_config": str(baseline_wave_path.relative_to(root)),
            "wave_sha256": sha256(baseline_wave_path),
        },
        "generated": generated,
    }
    manifest_path = config_dir / "manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--label", default=DEFAULT_LABEL)
    parser.add_argument("--permeability", type=float, default=DEFAULT_PERMEABILITY_M2)
    parser.add_argument("--saturation", type=float, default=DEFAULT_SATURATION)
    parser.add_argument("--equilibrium-steps", type=int, default=DEFAULT_EQ_STEPS)
    parser.add_argument("--wave-steps", type=int, default=DEFAULT_WAVE_STEPS)
    parser.add_argument(
        "--no-pressure-smoothing",
        action="store_true",
        help="Disable both state and VTK-time pressure smoothing",
    )
    parser.add_argument(
        "--paired-replay",
        action="store_true",
        help=(
            "Also generate a fixed-pipeline pressure writer and matched "
            "lagged/phase-erased replay configs"
        ),
    )
    args = parser.parse_args()
    root = Path(__file__).resolve().parent
    manifest = generate(
        root,
        label=args.label,
        permeability_m2=args.permeability,
        saturation=args.saturation,
        equilibrium_steps=args.equilibrium_steps,
        wave_steps=args.wave_steps,
        pressure_smoothing=not args.no_pressure_smoothing,
        paired_replay=args.paired_replay,
    )
    print(json.dumps(manifest, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
