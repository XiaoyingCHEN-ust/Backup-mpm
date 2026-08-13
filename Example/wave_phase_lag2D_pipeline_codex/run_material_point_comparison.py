#!/usr/bin/env python3
"""Run, audit, and plot the constitutive material-point diagnostic.

This workflow intentionally remains separate from the 2-D field simulations.
It compares constitutive mechanisms under an identical single-point strain
history; it is not evidence of field-scale predictive superiority.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
import subprocess
from collections import defaultdict
from pathlib import Path
from typing import Iterable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


SCRIPT_DIR = Path(__file__).resolve().parent
REPOSITORY_ROOT = SCRIPT_DIR.parents[1]
DEFAULT_DRIVER = (
    REPOSITORY_ROOT
    / "mpm_hpc_source"
    / "build-pipeline"
    / "cyclic_constitutive_comparison_driver"
)
DEFAULT_OUTPUT = SCRIPT_DIR / "analysis" / "material_point"
DRIVER_SOURCE = (
    REPOSITORY_ROOT
    / "mpm_hpc_source"
    / "tests"
    / "materials"
    / "cyclic_constitutive_comparison_driver.cc"
)
EXPECTED_PATH_STEPS = 2400
EXPECTED_DENSE_STEPS = 200
MC_YIELD_TOLERANCE_PA = 1.1e-6


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def numeric(row: dict[str, str], key: str) -> float:
    try:
        return float(row[key])
    except (KeyError, ValueError) as exception:
        raise ValueError(f"invalid numeric field {key!r} in CSV row") from exception


def parse_records(csv_text: str) -> list[dict[str, str]]:
    reader = csv.DictReader(io.StringIO(csv_text))
    required = {
        "record_type",
        "material",
        "path_id",
        "step",
        "cycle_time",
        "cycle_index",
        "cycle_coordinate",
        "reversal_count",
        "gamma_xy",
        "dgamma_xy",
        "volumetric_strain",
        "deviatoric_strain",
        "p_effective_pa",
        "signed_q_pa",
        "plastic_strain",
        "void_ratio",
        "alpha_norm",
        "alpha_initial_norm",
        "fabric_norm",
        "mc_tension_residual_pa",
        "mc_shear_residual_pa",
        "reference_pressure_pa",
        "measured_tangent_shear_pa",
    }
    fields = set(reader.fieldnames or [])
    missing = sorted(required - fields)
    if missing:
        raise ValueError(f"material-point CSV missing columns: {missing}")
    rows = list(reader)
    if not rows:
        raise ValueError("material-point CSV has no data rows")
    return rows


def finite_values(rows: Iterable[dict[str, str]], key: str) -> list[float]:
    values = [numeric(row, key) for row in rows]
    if not all(math.isfinite(value) for value in values):
        raise ValueError(f"non-finite {key} in required material-point path")
    return values


def validate_records(rows: list[dict[str, str]]) -> dict[str, object]:
    grouped: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        grouped[(row["record_type"], row["material"])].append(row)

    paths: dict[str, list[dict[str, str]]] = {}
    for material in ("SANISAND", "MohrCoulomb"):
        path = grouped.get(("path", material), [])
        if len(path) != EXPECTED_PATH_STEPS + 1:
            raise ValueError(
                f"{material} path has {len(path)} rows; "
                f"expected {EXPECTED_PATH_STEPS + 1}"
            )
        steps = [int(numeric(row, "step")) for row in path]
        if steps != list(range(EXPECTED_PATH_STEPS + 1)):
            raise ValueError(f"{material} path steps are incomplete or reordered")
        if {row["path_id"] for row in path} != {"loose_cyclic_simple_shear"}:
            raise ValueError(f"{material} has an unexpected path identifier")
        finite_values(path, "p_effective_pa")
        finite_values(path, "signed_q_pa")
        finite_values(path, "plastic_strain")
        paths[material] = path

    # Exact prescribed-history equality is checked as serialized by the
    # driver, independently of constitutive response.
    history_fields = (
        "step",
        "cycle_index",
        "cycle_coordinate",
        "reversal_count",
        "loading_direction",
        "gamma_xy",
        "dgamma_xy",
        "volumetric_strain",
    )
    for left, right in zip(paths["SANISAND"], paths["MohrCoulomb"]):
        if any(left[field] != right[field] for field in history_fields):
            raise ValueError("constitutive models did not receive identical strain histories")

    sani = paths["SANISAND"]
    mc = paths["MohrCoulomb"]
    for material_rows in (sani, mc):
        initial = material_rows[0]
        if abs(numeric(initial, "p_effective_pa") - 3000.0) > 1e-9:
            raise ValueError("common initial pressure is not 3 kPa")
        if abs(numeric(initial, "gamma_xy")) > 1e-15:
            raise ValueError("cyclic path does not start at zero shear strain")
        final = material_rows[-1]
        if abs(numeric(final, "gamma_xy")) > 1e-15:
            raise ValueError("cyclic path does not close at zero shear strain")
        if int(numeric(final, "reversal_count")) != 12:
            raise ValueError("cyclic path does not contain twelve reversals")
        if int(numeric(final, "cycle_index")) != 5 or abs(
            numeric(final, "cycle_coordinate") - 1.0
        ) > 1e-15:
            raise ValueError("six-cycle endpoint is misregistered")
        if max(abs(value) for value in finite_values(material_rows, "volumetric_strain")) > 1e-15:
            raise ValueError("primary comparison is not constant volume")

    sani_pressure = finite_values(sani, "p_effective_pa")
    mc_pressure = finite_values(mc, "p_effective_pa")
    if min(sani_pressure) <= 100.0:
        raise ValueError("SANISAND path crossed the registered pressure floor")
    if min(sani_pressure) >= 0.95 * 3000.0:
        raise ValueError("loose SANISAND path did not show contractive pressure loss")
    if max(mc_pressure) - min(mc_pressure) > 1e-6:
        raise ValueError("zero-dilation Mohr-Coulomb changed mean pressure")
    if numeric(sani[-1], "plastic_strain") <= 1e-5 or numeric(
        mc[-1], "plastic_strain"
    ) <= 1e-5:
        raise ValueError("both constitutive models were not plastically mobilised")
    if max(finite_values(sani, "alpha_norm")) <= 1e-4:
        raise ValueError("SANISAND back-stress did not evolve")

    mc_positive_residual = max(
        0.0,
        *(numeric(row, key) for row in mc for key in (
            "mc_tension_residual_pa",
            "mc_shear_residual_pa",
        )),
    )
    if mc_positive_residual > MC_YIELD_TOLERANCE_PA:
        raise ValueError(
            f"Mohr-Coulomb active-surface residual {mc_positive_residual:g} Pa "
            f"exceeds {MC_YIELD_TOLERANCE_PA:g} Pa"
        )

    dense = grouped.get(("state_probe", "SANISAND"), [])
    if len(dense) != EXPECTED_DENSE_STEPS + 1:
        raise ValueError("dense SANISAND state probe is incomplete")
    if [int(numeric(row, "step")) for row in dense] != list(
        range(EXPECTED_DENSE_STEPS + 1)
    ):
        raise ValueError("dense SANISAND state-probe steps are reordered")
    if {row["path_id"] for row in dense} != {"dense_monotonic_constant_volume"}:
        raise ValueError("dense SANISAND probe has an unexpected path identifier")
    dense_pressure = finite_values(dense, "p_effective_pa")
    if min(dense_pressure[1:]) >= 3000.0 or dense_pressure[-1] <= 3600.0:
        raise ValueError("dense SANISAND probe did not cross into dilation")
    dense_fabric = finite_values(dense, "fabric_norm")
    if max(dense_fabric) <= 0.1:
        raise ValueError("dense SANISAND probe did not evolve fabric")
    dense_void = finite_values(dense, "void_ratio")
    if max(dense_void) - min(dense_void) > 1e-10:
        raise ValueError("dense state probe is not constant volume")

    tangent_rows = [row for row in rows if row["record_type"] == "tangent_audit"]
    tangent: dict[str, dict[float, float]] = defaultdict(dict)
    for row in tangent_rows:
        tangent[row["material"]][numeric(row, "reference_pressure_pa")] = numeric(
            row, "measured_tangent_shear_pa"
        )
    expected_pressures = {1500.0, 3000.0, 6000.0}
    if any(set(tangent[material]) != expected_pressures for material in paths):
        raise ValueError("pressure-tangent audit grid is incomplete")
    sani_tangent = tangent["SANISAND"]
    mc_tangent = tangent["MohrCoulomb"]
    if not sani_tangent[1500.0] < sani_tangent[3000.0] < sani_tangent[6000.0]:
        raise ValueError("SANISAND tangent is not pressure dependent")
    if abs(mc_tangent[6000.0] - mc_tangent[1500.0]) / mc_tangent[3000.0] > 1e-12:
        raise ValueError("Mohr-Coulomb tangent unexpectedly depends on pressure")
    tangent_mismatch = abs(sani_tangent[3000.0] - mc_tangent[3000.0]) / mc_tangent[
        3000.0
    ]
    if tangent_mismatch > 1e-3:
        raise ValueError("3 kPa initial tangent match exceeds 0.1%")

    return {
        "schema": "pipeline-material-point-audit-v1",
        "scope": "material_point_mechanism_diagnostic_not_field_validation",
        "status": "pass",
        "registered_path": {
            "initial_mean_effective_pressure_pa": 3000.0,
            "porosity": numeric(sani[0], "void_ratio")
            / (1.0 + numeric(sani[0], "void_ratio")),
            "gamma_amplitude": max(abs(value) for value in finite_values(sani, "gamma_xy")),
            "cycles": 6,
            "reversals": int(numeric(sani[-1], "reversal_count")),
            "steps": EXPECTED_PATH_STEPS,
            "volumetric_strain": 0.0,
        },
        "sanisand": {
            "minimum_pressure_pa": min(sani_pressure),
            "maximum_pressure_pa": max(sani_pressure),
            "final_plastic_strain": numeric(sani[-1], "plastic_strain"),
            "maximum_alpha_norm": max(finite_values(sani, "alpha_norm")),
            "maximum_alpha_initial_norm": max(
                finite_values(sani, "alpha_initial_norm")
            ),
            "maximum_fabric_norm_loose_path": max(
                finite_values(sani, "fabric_norm")
            ),
        },
        "mohr_coulomb": {
            "minimum_pressure_pa": min(mc_pressure),
            "maximum_pressure_pa": max(mc_pressure),
            "final_plastic_strain": numeric(mc[-1], "plastic_strain"),
            "maximum_positive_yield_residual_pa": mc_positive_residual,
        },
        "dense_sanisand_probe": {
            "porosity": numeric(dense[0], "void_ratio")
            / (1.0 + numeric(dense[0], "void_ratio")),
            "minimum_pressure_pa": min(dense_pressure),
            "final_pressure_pa": dense_pressure[-1],
            "maximum_fabric_norm": max(dense_fabric),
        },
        "tangent_audit": {
            "pressures_pa": sorted(expected_pressures),
            "sanisand_shear_modulus_pa": [
                sani_tangent[value] for value in sorted(expected_pressures)
            ],
            "mohr_coulomb_shear_modulus_pa": [
                mc_tangent[value] for value in sorted(expected_pressures)
            ],
            "relative_mismatch_at_3kpa": tangent_mismatch,
        },
        "interpretation_guardrail": {
            "supports": [
                "pressure-dependent elastic tangent",
                "contractive and dilative state response",
                "back-stress and fabric history evolution",
            ],
            "does_not_support": [
                "field-scale validation",
                "field-scale predictive superiority of one constitutive model",
                "replacement of the registered 2-D phase-lag comparisons",
            ],
        },
    }


def plot_records(rows: list[dict[str, str]], output_base: Path) -> None:
    paths = {
        material: [
            row
            for row in rows
            if row["record_type"] == "path" and row["material"] == material
        ]
        for material in ("SANISAND", "MohrCoulomb")
    }
    dense = [row for row in rows if row["record_type"] == "state_probe"]
    tangent = [row for row in rows if row["record_type"] == "tangent_audit"]

    colors = {"SANISAND": "#176D9C", "MohrCoulomb": "#D17A22"}
    figure, axes = plt.subplots(2, 3, figsize=(13.0, 8.2), constrained_layout=True)

    ax = axes[0, 0]
    for material, path in paths.items():
        ax.plot(
            [numeric(row, "p_effective_pa") / 1000.0 for row in path],
            [numeric(row, "signed_q_pa") / 1000.0 for row in path],
            color=colors[material],
            linewidth=1.15,
            label=material,
        )
    ax.set(xlabel="$p'$ (kPa)", ylabel="signed $q$ (kPa)", title="(a) Cyclic stress path")
    ax.legend(frameon=False)

    ax = axes[0, 1]
    for material, path in paths.items():
        ax.plot(
            [numeric(row, "cycle_time") for row in path],
            [100.0 * numeric(row, "plastic_strain") for row in path],
            color=colors[material],
            linewidth=1.2,
            label=material,
        )
    ax.set(
        xlabel="prescribed cycles",
        ylabel="accumulated plastic indicator (%)",
        title="(b) Cyclic accumulation",
    )

    ax = axes[0, 2]
    for material, path in paths.items():
        ax.plot(
            [numeric(row, "cycle_time") for row in path],
            [numeric(row, "p_effective_pa") / 3000.0 for row in path],
            color=colors[material],
            linewidth=1.2,
            label=material,
        )
    ax.axhline(1.0, color="0.65", linewidth=0.8, linestyle=":")
    ax.set(xlabel="prescribed cycles", ylabel="$p'/p'_0$", title="(c) Constant-volume response")

    ax = axes[1, 0]
    for material in paths:
        selected = sorted(
            (row for row in tangent if row["material"] == material),
            key=lambda row: numeric(row, "reference_pressure_pa"),
        )
        ax.plot(
            [numeric(row, "reference_pressure_pa") / 1000.0 for row in selected],
            [numeric(row, "measured_tangent_shear_pa") / 1.0e6 for row in selected],
            "o-",
            color=colors[material],
            linewidth=1.2,
            markersize=4,
            label=material,
        )
    ax.set(xlabel="reference $p'$ (kPa)", ylabel="initial tangent $G$ (MPa)", title="(d) Pressure-dependent stiffness")

    ax = axes[1, 1]
    sani = paths["SANISAND"]
    cycles = [numeric(row, "cycle_time") for row in sani]
    ax.plot(cycles, [numeric(row, "alpha_norm") for row in sani], color="#176D9C", label="$|\\alpha|$")
    ax.plot(
        cycles,
        [numeric(row, "alpha_initial_norm") for row in sani],
        color="#5B3C88",
        linestyle="--",
        label="$|\\alpha_{in}|$",
    )
    ax.set(xlabel="prescribed cycles", ylabel="dimensionless state norm", title="(e) SANISAND reversal memory")
    ax.legend(frameon=False)

    ax = axes[1, 2]
    dense_strain = [100.0 * numeric(row, "deviatoric_strain") for row in dense]
    pressure_line = ax.plot(
        dense_strain,
        [numeric(row, "p_effective_pa") / 3000.0 for row in dense],
        color="#2B8C5B",
        label="$p'/p'_0$",
    )[0]
    ax.set(xlabel="dense probe deviatoric strain (%)", ylabel="$p'/p'_0$", title="(f) Dilatancy and fabric probe")
    twin = ax.twinx()
    fabric_line = twin.plot(
        dense_strain,
        [numeric(row, "fabric_norm") for row in dense],
        color="#9C3D54",
        linestyle="--",
        label="$|Z|$",
    )[0]
    twin.set_ylabel("fabric norm $|Z|$")
    ax.legend([pressure_line, fabric_line], ["$p'/p'_0$", "$|Z|$"], frameon=False)

    for axis in axes.flat:
        axis.grid(True, color="0.9", linewidth=0.6)
        axis.tick_params(direction="in")
    figure.suptitle(
        "Material-point mechanism diagnostic — not field validation",
        fontsize=13,
        fontweight="semibold",
    )
    figure.savefig(output_base.with_suffix(".png"), dpi=300)
    figure.savefig(output_base.with_suffix(".pdf"))
    plt.close(figure)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--driver", type=Path, default=DEFAULT_DRIVER)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--skip-plot", action="store_true")
    args = parser.parse_args()

    driver = args.driver.resolve()
    if not driver.is_file():
        raise SystemExit(
            f"material-point driver not found: {driver}\n"
            "Build target cyclic_constitutive_comparison_driver first."
        )
    completed = subprocess.run(
        [str(driver)],
        check=False,
        capture_output=True,
        text=True,
        timeout=120,
    )
    if completed.returncode != 0:
        raise SystemExit(
            f"material-point driver failed with code {completed.returncode}:\n"
            f"{completed.stderr.strip()}"
        )

    rows = parse_records(completed.stdout)
    audit = validate_records(rows)
    audit["provenance"] = {
        "driver": str(driver),
        "driver_sha256": sha256_file(driver),
        "driver_source": str(DRIVER_SOURCE),
        "driver_source_sha256": sha256_file(DRIVER_SOURCE),
        "driver_stderr": completed.stderr.strip(),
    }

    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    csv_path = output / "cyclic_constitutive_comparison.csv"
    audit_path = output / "cyclic_constitutive_comparison.audit.json"
    log_path = output / "cyclic_constitutive_comparison.log"
    csv_path.write_text(completed.stdout, encoding="utf-8")
    audit_path.write_text(
        json.dumps(audit, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    log_path.write_text(completed.stderr, encoding="utf-8")
    if not args.skip_plot:
        plot_records(rows, output / "cyclic_constitutive_comparison")

    print(f"PASS: {audit_path}")
    print(f"CSV:  {csv_path}")
    if not args.skip_plot:
        print(f"FIG:  {output / 'cyclic_constitutive_comparison.png'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
