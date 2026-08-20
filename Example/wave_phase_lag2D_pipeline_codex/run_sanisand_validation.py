#!/usr/bin/env python3
"""Independent, publication-audited validation of the repository SANISAND.

The in-repository C++ material-point driver is compared with a separately
compiled, commit-pinned SANISAND-MS UDSM from Federico Pisano's repository.
No upstream GPL source is copied into this repository or linked into MPM.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import json
import math
import os
import shutil
import subprocess
import sys
import tempfile
from collections import defaultdict
from pathlib import Path
from typing import Any, Iterable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


SCHEMA = "pipeline-sanisand-independent-validation-v1"
EXPECTED_REFERENCE_COMMIT = "205c13b0a8fe5ffdcc1404b6fe63a59a67e613b9"
REFERENCE_URL = "https://github.com/FedericoPisano/SANISAND-MS-UDSM"
LIU_2019_DOI = "https://doi.org/10.1680/jgeot.17.P.307"
DM04_DOI = (
    "https://doi.org/10.1061/(ASCE)0733-9399(2004)130:6(622)"
)
REFERENCE_FILES = (
    "USRADDDF.f90",
    "USRLIB.for",
    "TensorOperation.f90",
    "StatedepSub.f90",
    "Initialize.f90",
    "Integration.f90",
    "SANISANDMS.f90",
)
PATH_MONOTONIC = "toyoura_monotonic_constant_volume"
PATH_CYCLIC = "toyoura_cyclic_constant_volume"
PATH_LITERATURE = "toyoura_literature_q_controlled"
IMPLEMENTATION_CURRENT = "mpm_sanisand04_style"
IMPLEMENTATION_REFERENCE = "pisano_sanisand_ms_udsm"
INITIAL_PRESSURE_PA = 294000.0
Q_AMPLITUDE_PA = 114200.0
TIGHT_TOLERANCE = 1.0e-7
LOOSE_TOLERANCE = 1.0e-5
EXPECTED_FIELDS = (
    "implementation",
    "path_id",
    "tolerance",
    "step",
    "cycle_time",
    "axial_strain",
    "p_effective_pa",
    "q_signed_pa",
    "ru",
    "void_ratio",
    "plastic_strain",
    "alpha_norm",
    "alpha_in_norm",
    "alpha_memory_norm",
    "memory_radius",
    "fabric_norm",
    "controller_target_q_pa",
    "controller_residual_pa",
    "status",
)
CORE_NUMERIC = (
    "tolerance",
    "cycle_time",
    "axial_strain",
    "p_effective_pa",
    "q_signed_pa",
    "ru",
    "void_ratio",
    "alpha_norm",
    "alpha_in_norm",
)
OPTIONAL_NUMERIC = (
    "plastic_strain",
    "alpha_memory_norm",
    "memory_radius",
    "fabric_norm",
    "controller_target_q_pa",
    "controller_residual_pa",
)


class ValidationError(RuntimeError):
    pass


def triangle(cycle_time: float) -> float:
    fraction = cycle_time - math.floor(cycle_time)
    if fraction < 0.25:
        return 4.0 * fraction
    if fraction < 0.75:
        return 2.0 - 4.0 * fraction
    return -4.0 + 4.0 * fraction


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def artifact(path: Path) -> dict[str, Any]:
    resolved = path.resolve()
    if not resolved.is_file():
        raise ValidationError(f"artifact is missing: {resolved}")
    return {
        "path": str(resolved),
        "size_bytes": resolved.stat().st_size,
        "sha256": sha256(resolved),
    }


def published_artifact(staged_path: Path, final_path: Path) -> dict[str, Any]:
    """Hash a staged file while recording its post-transaction path."""
    if not staged_path.is_file():
        raise ValidationError(f"staged artifact is missing: {staged_path}")
    return {
        "path": str(final_path.resolve()),
        "size_bytes": staged_path.stat().st_size,
        "sha256": sha256(staged_path),
    }


def run(command: list[str], *, cwd: Path | None = None) -> subprocess.CompletedProcess[str]:
    completed = subprocess.run(
        command,
        cwd=cwd,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    if completed.returncode != 0:
        rendered = " ".join(command)
        raise ValidationError(
            f"command failed ({completed.returncode}): {rendered}\n"
            f"stdout:\n{completed.stdout[-4000:]}\n"
            f"stderr:\n{completed.stderr[-4000:]}"
        )
    return completed


def parse_float(value: str, field: str, *, optional: bool = False) -> float:
    try:
        parsed = float(value)
    except ValueError as exc:
        raise ValidationError(f"invalid {field}: {value!r}") from exc
    if not optional and not math.isfinite(parsed):
        raise ValidationError(f"non-finite required field {field}")
    if optional and not (math.isfinite(parsed) or math.isnan(parsed)):
        raise ValidationError(f"invalid optional numeric field {field}")
    return parsed


def load_rows(path: Path, expected_implementation: str) -> list[dict[str, Any]]:
    with path.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        if tuple(reader.fieldnames or ()) != EXPECTED_FIELDS:
            raise ValidationError(
                f"unexpected CSV schema for {path}: {reader.fieldnames}"
            )
        raw_rows = list(reader)
    if not raw_rows:
        raise ValidationError(f"empty validation CSV: {path}")
    rows: list[dict[str, Any]] = []
    valid_paths = {PATH_MONOTONIC, PATH_CYCLIC, PATH_LITERATURE}
    for raw in raw_rows:
        if raw["implementation"] != expected_implementation:
            raise ValidationError("CSV implementation identity drift")
        if raw["path_id"] not in valid_paths:
            raise ValidationError(f"unknown validation path: {raw['path_id']}")
        if raw["status"] not in {"ok", "controller_limit"}:
            raise ValidationError(f"unknown row status: {raw['status']}")
        try:
            step = int(raw["step"])
        except ValueError as exc:
            raise ValidationError(f"invalid step: {raw['step']!r}") from exc
        if step < 0:
            raise ValidationError("negative material-point step")
        row: dict[str, Any] = dict(raw)
        row["step"] = step
        for field in CORE_NUMERIC:
            row[field] = parse_float(raw[field], field)
        for field in OPTIONAL_NUMERIC:
            row[field] = parse_float(raw[field], field, optional=True)
        if row["p_effective_pa"] <= 0.0:
            raise ValidationError("non-compressive mean effective stress")
        if abs(row["void_ratio"] - 0.808) > 1.0e-9:
            raise ValidationError("constant-volume path changed the void ratio")
        if row["status"] == "controller_limit" and row["path_id"] != PATH_LITERATURE:
            raise ValidationError("controller-limit status on a strain path")
        if row["path_id"] == PATH_LITERATURE:
            target = row["controller_target_q_pa"]
            residual = row["controller_residual_pa"]
            if not math.isfinite(target) or not math.isfinite(residual):
                raise ValidationError("q-controlled row lacks controller evidence")
            expected_target = Q_AMPLITUDE_PA * triangle(row["cycle_time"])
            if not math.isclose(target, expected_target, rel_tol=0.0, abs_tol=1e-6):
                raise ValidationError("q-controller target history drift")
            limit = 0.00501 * Q_AMPLITUDE_PA
            if row["status"] == "ok" and abs(residual) > limit:
                raise ValidationError("accepted q-controller residual exceeds limit")
            if row["status"] == "controller_limit" and abs(residual) <= limit:
                raise ValidationError("controller terminal does not exceed its limit")
        elif not (
            math.isnan(row["controller_target_q_pa"])
            and math.isnan(row["controller_residual_pa"])
        ):
            raise ValidationError("strain-controlled row contains controller evidence")
        rows.append(row)
    validate_row_grids(rows)
    return rows


def select(
    rows: Iterable[dict[str, Any]],
    implementation: str,
    path_id: str,
    tolerance: float,
) -> list[dict[str, Any]]:
    selected = [
        row
        for row in rows
        if row["implementation"] == implementation
        and row["path_id"] == path_id
        and math.isclose(row["tolerance"], tolerance, rel_tol=0.0, abs_tol=1e-12)
    ]
    return sorted(selected, key=lambda row: row["step"])


def validate_row_grids(rows: list[dict[str, Any]]) -> None:
    groups: dict[tuple[str, str, float], list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        groups[(row["implementation"], row["path_id"], row["tolerance"])].append(row)
    for key, group in groups.items():
        ordered = sorted(group, key=lambda row: row["step"])
        steps = [row["step"] for row in ordered]
        if len(steps) != len(set(steps)) or steps[0] != 0:
            raise ValidationError(f"duplicate or missing initial step in {key}")
        limits = [row for row in ordered if row["status"] == "controller_limit"]
        if len(limits) > 1 or (limits and limits[0] is not ordered[-1]):
            raise ValidationError(f"invalid controller terminal record in {key}")
        if key[1] == PATH_MONOTONIC:
            expected = list(range(0, 1001, 5))
        elif key[1] == PATH_CYCLIC:
            expected = list(range(0, 3201, 5))
        else:
            expected = list(range(0, steps[-1] + 1))
        if steps != expected:
            raise ValidationError(f"unexpected step grid in {key}")


def require_common_prescribed_strain_grid(
    current: list[dict[str, Any]], reference: list[dict[str, Any]], path_id: str
) -> None:
    for tolerance in (LOOSE_TOLERANCE, TIGHT_TOLERANCE):
        left = select(current, IMPLEMENTATION_CURRENT, path_id, tolerance)
        right = select(reference, IMPLEMENTATION_REFERENCE, path_id, tolerance)
        if [row["step"] for row in left] != [row["step"] for row in right]:
            raise ValidationError(f"current/reference step mismatch for {path_id}")
        if not np.allclose(
            [row["axial_strain"] for row in left],
            [row["axial_strain"] for row in right],
            rtol=0.0,
            atol=2e-15,
        ):
            raise ValidationError(f"prescribed strain mismatch for {path_id}")


def normalized_error(
    left: list[dict[str, Any]],
    right: list[dict[str, Any]],
    field: str,
    scale: float,
) -> dict[str, float]:
    right_by_step = {row["step"]: row for row in right if row["status"] == "ok"}
    pairs = [
        (row, right_by_step[row["step"]])
        for row in left
        if row["status"] == "ok" and row["step"] in right_by_step
    ]
    if len(pairs) < 2:
        raise ValidationError(f"insufficient paired rows for {field}")
    differences = np.asarray(
        [float(a[field]) - float(b[field]) for a, b in pairs], dtype=float
    )
    return {
        "paired_samples": len(pairs),
        "rmse_normalized": float(np.sqrt(np.mean(differences**2)) / scale),
        "maximum_absolute_normalized": float(np.max(np.abs(differences)) / scale),
    }


def convergence(
    rows: list[dict[str, Any]], implementation: str, path_id: str
) -> dict[str, Any]:
    loose = select(rows, implementation, path_id, LOOSE_TOLERANCE)
    tight = select(rows, implementation, path_id, TIGHT_TOLERANCE)
    pressure = normalized_error(loose, tight, "p_effective_pa", INITIAL_PRESSURE_PA)
    deviator = normalized_error(loose, tight, "q_signed_pa", Q_AMPLITUDE_PA)
    maximum = max(
        pressure["maximum_absolute_normalized"],
        deviator["maximum_absolute_normalized"],
    )
    return {
        "pressure": pressure,
        "deviator_stress": deviator,
        "maximum_normalized_difference": maximum,
        "limit": 0.01,
        "passed": maximum <= 0.01,
    }


def threshold_cycle(rows: list[dict[str, Any]], threshold: float) -> float | None:
    for row in rows:
        if row["ru"] >= threshold:
            return float(row["cycle_time"])
    return None


def terminal_record(rows: list[dict[str, Any]]) -> dict[str, Any]:
    terminals = [row for row in rows if row["status"] == "controller_limit"]
    if len(terminals) != 1:
        raise ValidationError("literature path must have one explicit terminal record")
    terminal = terminals[0]
    ok = [row for row in rows if row["status"] == "ok"]
    finite_residuals = [
        abs(float(row["controller_residual_pa"]))
        for row in ok
        if math.isfinite(float(row["controller_residual_pa"]))
    ]
    if not finite_residuals:
        raise ValidationError("q-controlled path has no finite controller residuals")
    maximum_residual = max(finite_residuals)
    return {
        "step": terminal["step"],
        "cycle_time": terminal["cycle_time"],
        "p_effective_pa": terminal["p_effective_pa"],
        "ru": terminal["ru"],
        "axial_strain": terminal["axial_strain"],
        "terminal_residual_pa": terminal["controller_residual_pa"],
        "maximum_accepted_controller_residual_pa": maximum_residual,
        "controller_residual_limit_pa": 0.00501 * Q_AMPLITUDE_PA,
    }


def cycle_rows(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    result: list[dict[str, Any]] = []
    groups: dict[tuple[str, str, float, int], list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        if row["step"] == 0:
            continue
        cycle = int(math.ceil(row["cycle_time"] - 1e-12))
        groups[(row["implementation"], row["path_id"], row["tolerance"], cycle)].append(row)
    for (implementation, path_id, tolerance, cycle), group in sorted(groups.items()):
        result.append(
            {
                "implementation": implementation,
                "path_id": path_id,
                "tolerance": tolerance,
                "cycle": cycle,
                "sample_count": len(group),
                "minimum_p_effective_pa": min(row["p_effective_pa"] for row in group),
                "maximum_ru": max(row["ru"] for row in group),
                "minimum_q_signed_pa": min(row["q_signed_pa"] for row in group),
                "maximum_q_signed_pa": max(row["q_signed_pa"] for row in group),
                "terminal_in_cycle": any(row["status"] != "ok" for row in group),
            }
        )
    return result


def compute_metrics(
    current: list[dict[str, Any]], reference: list[dict[str, Any]]
) -> dict[str, Any]:
    require_common_prescribed_strain_grid(current, reference, PATH_MONOTONIC)
    require_common_prescribed_strain_grid(current, reference, PATH_CYCLIC)
    all_rows = current + reference
    current_monotonic = select(
        current, IMPLEMENTATION_CURRENT, PATH_MONOTONIC, TIGHT_TOLERANCE
    )
    reference_monotonic = select(
        reference, IMPLEMENTATION_REFERENCE, PATH_MONOTONIC, TIGHT_TOLERANCE
    )
    current_cyclic = select(
        current, IMPLEMENTATION_CURRENT, PATH_CYCLIC, TIGHT_TOLERANCE
    )
    reference_cyclic = select(
        reference, IMPLEMENTATION_REFERENCE, PATH_CYCLIC, TIGHT_TOLERANCE
    )
    current_q = select(
        current, IMPLEMENTATION_CURRENT, PATH_LITERATURE, TIGHT_TOLERANCE
    )
    reference_q = select(
        reference, IMPLEMENTATION_REFERENCE, PATH_LITERATURE, TIGHT_TOLERANCE
    )
    convergence_records = {
        implementation: {
            path: convergence(
                current if implementation == IMPLEMENTATION_CURRENT else reference,
                implementation,
                path,
            )
            for path in (PATH_MONOTONIC, PATH_CYCLIC, PATH_LITERATURE)
        }
        for implementation in (IMPLEMENTATION_CURRENT, IMPLEMENTATION_REFERENCE)
    }
    convergence_passed = all(
        record["passed"]
        for implementation in convergence_records.values()
        for record in implementation.values()
    )
    monotonic = {
        "pressure": normalized_error(
            current_monotonic,
            reference_monotonic,
            "p_effective_pa",
            INITIAL_PRESSURE_PA,
        ),
        "deviator_stress": normalized_error(
            current_monotonic,
            reference_monotonic,
            "q_signed_pa",
            Q_AMPLITUDE_PA,
        ),
    }
    monotonic_limit = 0.10
    monotonic_passed = max(
        monotonic["pressure"]["maximum_absolute_normalized"],
        monotonic["deviator_stress"]["maximum_absolute_normalized"],
    ) <= monotonic_limit
    cyclic = {
        "pressure": normalized_error(
            current_cyclic, reference_cyclic, "p_effective_pa", INITIAL_PRESSURE_PA
        ),
        "deviator_stress": normalized_error(
            current_cyclic, reference_cyclic, "q_signed_pa", Q_AMPLITUDE_PA
        ),
        "current_ru_95_cycle": threshold_cycle(current_cyclic, 0.95),
        "reference_ru_95_cycle": threshold_cycle(reference_cyclic, 0.95),
        "current_ru_99_cycle": threshold_cycle(current_cyclic, 0.99),
        "reference_ru_99_cycle": threshold_cycle(reference_cyclic, 0.99),
    }
    current_terminal = terminal_record(current_q)
    reference_terminal = terminal_record(reference_q)
    memory_surface_effect_resolved = (
        reference_terminal["cycle_time"] - current_terminal["cycle_time"] >= 0.5
    )
    return {
        "initial_state": {
            "p_effective_pa": INITIAL_PRESSURE_PA,
            "void_ratio": 0.808,
            "literature_q_amplitude_pa": Q_AMPLITUDE_PA,
        },
        "integration_convergence": convergence_records,
        "integration_convergence_passed": convergence_passed,
        "monotonic_code_to_code": monotonic,
        "monotonic_maximum_difference_limit": monotonic_limit,
        "monotonic_virgin_loading_passed": monotonic_passed,
        "cyclic_same_strain": cyclic,
        "literature_q_controlled": {
            IMPLEMENTATION_CURRENT: current_terminal,
            IMPLEMENTATION_REFERENCE: reference_terminal,
            "memory_surface_extension_cycles": reference_terminal["cycle_time"]
            - current_terminal["cycle_time"],
            "memory_surface_effect_resolved": memory_surface_effect_resolved,
        },
        "verification_passed": convergence_passed and monotonic_passed,
        "sanisand_ms_equivalence_passed": False,
        "scope_decision": (
            "The repository implementation is numerically verified and agrees "
            "with the independent UDSM on virgin monotonic loading within the "
            "predeclared 10% normalized maximum-difference limit. Its cyclic "
            "response is not equivalent to SANISAND-MS because it has no "
            "memory-surface state variables."
        ),
        "row_count": len(all_rows),
    }


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ValidationError(f"cannot write empty CSV: {path}")
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _style() -> None:
    plt.rcParams.update(
        {
            "font.size": 9,
            "axes.labelsize": 9,
            "axes.titlesize": 9,
            "legend.fontsize": 8,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "axes.linewidth": 0.8,
            "lines.linewidth": 1.6,
            "savefig.dpi": 300,
            "figure.dpi": 120,
            "font.family": "DejaVu Sans",
        }
    )


def _plot_line(ax: Any, rows: list[dict[str, Any]], x: str, y: str, *,
               color: str, label: str, linestyle: str = "-",
               x_scale: float = 1.0, y_scale: float = 1.0) -> None:
    ok = [row for row in rows if row["status"] == "ok"]
    ax.plot(
        [row[x] / x_scale for row in ok],
        [row[y] / y_scale for row in ok],
        color=color,
        linestyle=linestyle,
        label=label,
        zorder=2,
    )
    terminal = [row for row in rows if row["status"] == "controller_limit"]
    if terminal:
        ax.scatter(
            [terminal[0][x] / x_scale], [terminal[0][y] / y_scale],
            marker="x", s=40,
            color=color, linewidths=1.8, zorder=4,
        )


def plot_validation_paths(
    current: list[dict[str, Any]], reference: list[dict[str, Any]], output: Path
) -> list[Path]:
    _style()
    blue = "#0072B2"
    orange = "#D55E00"
    cur_m = select(current, IMPLEMENTATION_CURRENT, PATH_MONOTONIC, TIGHT_TOLERANCE)
    ref_m = select(reference, IMPLEMENTATION_REFERENCE, PATH_MONOTONIC, TIGHT_TOLERANCE)
    cur_c = select(current, IMPLEMENTATION_CURRENT, PATH_CYCLIC, TIGHT_TOLERANCE)
    ref_c = select(reference, IMPLEMENTATION_REFERENCE, PATH_CYCLIC, TIGHT_TOLERANCE)
    fig, axes = plt.subplots(2, 2, figsize=(7.15, 5.3), constrained_layout=True)
    ax = axes[0, 0]
    _plot_line(ax, cur_m, "axial_strain", "q_signed_pa", color=blue,
               label="MPM SANISAND04-style", y_scale=1000.0)
    _plot_line(ax, ref_m, "axial_strain", "q_signed_pa", color=orange,
               label="Pisano SANISAND-MS UDSM", linestyle="--", y_scale=1000.0)
    ax.set(xlabel=r"Axial strain, $\varepsilon_a$", ylabel=r"$q$ (kPa)",
           title="(a) Virgin monotonic response")
    ax.legend(frameon=False)
    ax = axes[0, 1]
    _plot_line(ax, cur_m, "axial_strain", "p_effective_pa", color=blue,
               label="MPM SANISAND04-style", y_scale=1000.0)
    _plot_line(ax, ref_m, "axial_strain", "p_effective_pa", color=orange,
               label="Pisano SANISAND-MS UDSM", linestyle="--", y_scale=1000.0)
    ax.set(xlabel=r"Axial strain, $\varepsilon_a$", ylabel=r"$p'$ (kPa)",
           title="(b) Mean effective stress")
    ax = axes[1, 0]
    _plot_line(ax, cur_c, "axial_strain", "q_signed_pa", color=blue,
               label="MPM SANISAND04-style", y_scale=1000.0)
    _plot_line(ax, ref_c, "axial_strain", "q_signed_pa", color=orange,
               label="Pisano SANISAND-MS UDSM", linestyle="--", y_scale=1000.0)
    ax.set(xlabel=r"Axial strain, $\varepsilon_a$", ylabel=r"$q$ (kPa)",
           title="(c) Identical constant-volume cycles")
    ax = axes[1, 1]
    _plot_line(ax, cur_c, "cycle_time", "ru", color=blue,
               label="MPM SANISAND04-style")
    _plot_line(ax, ref_c, "cycle_time", "ru", color=orange,
               label="Pisano SANISAND-MS UDSM", linestyle="--")
    ax.axhline(0.95, color="0.55", linestyle=":", linewidth=1.0)
    ax.set(xlabel="Cycle number", ylabel=r"$r_u=1-p'/p'_0$",
           title="(d) Pore-pressure accumulation", ylim=(-0.03, 1.04))
    for axis in axes.flat:
        axis.grid(True, color="0.9", linewidth=0.6)
    png = output / "figure_sanisand_validation_paths.png"
    pdf = output / "figure_sanisand_validation_paths.pdf"
    fig.savefig(png, dpi=300)
    fig.savefig(pdf)
    plt.close(fig)
    return [png, pdf]


def plot_literature_path(
    current: list[dict[str, Any]], reference: list[dict[str, Any]], output: Path
) -> list[Path]:
    _style()
    blue = "#0072B2"
    orange = "#D55E00"
    cur = select(current, IMPLEMENTATION_CURRENT, PATH_LITERATURE, TIGHT_TOLERANCE)
    ref = select(reference, IMPLEMENTATION_REFERENCE, PATH_LITERATURE, TIGHT_TOLERANCE)
    fig, axes = plt.subplots(2, 2, figsize=(7.15, 5.3), constrained_layout=True)
    ax = axes[0, 0]
    _plot_line(ax, cur, "p_effective_pa", "q_signed_pa", color=blue,
               label="MPM SANISAND04-style", x_scale=1000.0, y_scale=1000.0)
    _plot_line(ax, ref, "p_effective_pa", "q_signed_pa", color=orange,
               label="Pisano SANISAND-MS UDSM", linestyle="--",
               x_scale=1000.0, y_scale=1000.0)
    ax.set(xlabel=r"$p'$ (kPa)", ylabel=r"$q$ (kPa)",
           title="(a) Liu et al. (2019) Appendix-II path")
    ax.legend(frameon=False)
    ax = axes[0, 1]
    _plot_line(ax, cur, "cycle_time", "ru", color=blue,
               label="MPM SANISAND04-style")
    _plot_line(ax, ref, "cycle_time", "ru", color=orange,
               label="Pisano SANISAND-MS UDSM", linestyle="--")
    ax.axhline(0.95, color="0.55", linestyle=":", linewidth=1.0)
    ax.set(xlabel="Cycle number", ylabel=r"$r_u=1-p'/p'_0$",
           title="(b) Effective-stress loss", ylim=(-0.03, 1.04))
    ax = axes[1, 0]
    _plot_line(ax, cur, "cycle_time", "axial_strain", color=blue,
               label="MPM SANISAND04-style")
    _plot_line(ax, ref, "cycle_time", "axial_strain", color=orange,
               label="Pisano SANISAND-MS UDSM", linestyle="--")
    ax.set(xlabel="Cycle number", ylabel=r"$\varepsilon_a$",
           title="(c) Strain demand under prescribed q")
    ax = axes[1, 1]
    for data, color, label, linestyle in (
        (cur, blue, "MPM SANISAND04-style", "-"),
        (ref, orange, "Pisano SANISAND-MS UDSM", "--"),
    ):
        grouped: dict[int, list[float]] = defaultdict(list)
        for row in data:
            value = row["controller_residual_pa"]
            if row["status"] == "ok" and math.isfinite(value) and row["step"] > 0:
                cycle = int(math.ceil(row["cycle_time"] - 1e-12))
                grouped[cycle].append(abs(value) / 1000.0)
        x_values = sorted(grouped)
        ax.plot(
            x_values,
            [max(grouped[index]) for index in x_values],
            color=color,
            linestyle=linestyle,
            marker="o",
            markersize=4,
            label=label,
        )
        terminal = [row for row in data if row["status"] == "controller_limit"]
        if terminal:
            ax.scatter(
                [terminal[0]["cycle_time"]],
                [abs(terminal[0]["controller_residual_pa"]) / 1000.0],
                marker="x", s=40, color=color, linewidths=1.8, zorder=4,
            )
    limit = 0.00501 * Q_AMPLITUDE_PA / 1000.0
    ax.axhline(limit, color="0.55", linestyle=":", linewidth=1.0)
    ax.set(xlabel="Cycle number", ylabel="Max. |q-control residual| (kPa)",
           title="(d) Per-cycle controller audit (× = terminal)")
    for axis in axes.flat:
        axis.grid(True, color="0.9", linewidth=0.6)
    png = output / "figure_sanisand_literature_validation.png"
    pdf = output / "figure_sanisand_literature_validation.pdf"
    fig.savefig(png, dpi=300)
    fig.savefig(pdf)
    plt.close(fig)
    return [png, pdf]


def compile_reference(reference_root: Path, adapter: Path, build: Path) -> dict[str, Any]:
    commit = run(["git", "-C", str(reference_root), "rev-parse", "HEAD"]).stdout.strip()
    if commit != EXPECTED_REFERENCE_COMMIT:
        raise ValidationError(
            f"reference commit must be {EXPECTED_REFERENCE_COMMIT}, got {commit}"
        )
    source_records: list[dict[str, Any]] = []
    for name in REFERENCE_FILES:
        blob = subprocess.run(
            ["git", "-C", str(reference_root), "show", f"{commit}:{name}"],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
        )
        if blob.returncode != 0:
            raise ValidationError(f"cannot read committed reference source {name}")
        data = blob.stdout
        original_sha = sha256_bytes(data)
        replacements = 0
        if name == "StatedepSub.f90":
            needle = b"FORMAT(2I, 10F15.10)"
            replacements = data.count(needle)
            if replacements != 2:
                raise ValidationError(
                    "expected exactly two compiler-format compatibility replacements"
                )
            data = data.replace(needle, b"FORMAT(2I12, 10F15.10)")
        target = build / name
        target.write_bytes(data)
        source_records.append(
            {
                "path": name,
                "committed_sha256": original_sha,
                "compiled_sha256": sha256(target),
                "compatibility_replacements": replacements,
            }
        )
    adapter_copy = build / adapter.name
    shutil.copy2(adapter, adapter_copy)
    module_dir = build / "modules"
    module_dir.mkdir()
    executable = build / "pisano_udsm_material_point_driver"
    command = [
        "gfortran",
        "-O2",
        "-ffree-line-length-none",
        "-ffixed-line-length-none",
        "-std=legacy",
        f"-J{module_dir}",
        *[str(build / name) for name in REFERENCE_FILES],
        str(adapter_copy),
        "-o",
        str(executable),
    ]
    run(command, cwd=build)
    compiler = run(["gfortran", "--version"]).stdout.splitlines()[0]
    return {
        "repository_url": REFERENCE_URL,
        "commit": commit,
        "compiler": compiler,
        "compile_command": command,
        "sources": source_records,
        "adapter": artifact(adapter),
        "executable": artifact(executable),
    }


def build_current(repo_root: Path, binary: Path, skip_build: bool) -> None:
    if not skip_build:
        run(
            [
                "cmake",
                "--build",
                str(repo_root / "mpm_hpc_source" / "build-pipeline"),
                "--target",
                "sanisand_validation_driver",
                "--parallel",
                "2",
            ]
        )
    if not binary.is_file() or not os.access(binary, os.X_OK):
        raise ValidationError(f"current validation driver is unavailable: {binary}")


def markdown_summary(metrics: dict[str, Any]) -> str:
    q = metrics["literature_q_controlled"]
    current = q[IMPLEMENTATION_CURRENT]
    reference = q[IMPLEMENTATION_REFERENCE]
    mono = metrics["monotonic_code_to_code"]
    return f"""# SANISAND independent validation result

The numerical integration checks passed at both tolerances. Under virgin,
constant-volume monotonic loading, the maximum normalized current/reference
difference was {100*max(mono['pressure']['maximum_absolute_normalized'], mono['deviator_stress']['maximum_absolute_normalized']):.2f}%, below the predeclared 10% mapping limit.

Under the Liu et al. (2019) Appendix-II condition ($p'_0=294$ kPa,
$e_0=0.808$, $q_{{ampl}}=114.2$ kPa), the current implementation reached its
audited q-controller terminal at {current['cycle_time']:.3f} cycles and
$r_u={current['ru']:.4f}$.  The independent Pisano SANISAND-MS UDSM reached
the corresponding terminal at {reference['cycle_time']:.3f} cycles and
$r_u={reference['ru']:.4f}$, an extension of {q['memory_surface_extension_cycles']:.3f}
cycles.  The cyclic curves are therefore not equivalent.

The validated claim is limited to a SANISAND04-style, fabric-enhanced model;
the current 20-state implementation must not be labelled SANISAND-MS because
it lacks the UDSM's memory-surface centre and radius variables.
"""


def validate_artifacts(records: list[dict[str, Any]]) -> None:
    seen: set[str] = set()
    for record in records:
        if set(record) != {"path", "size_bytes", "sha256"}:
            raise ValidationError("artifact record schema drift")
        path = Path(record["path"])
        if str(path) in seen:
            raise ValidationError("duplicate artifact path")
        seen.add(str(path))
        if artifact(path) != record:
            raise ValidationError(f"artifact drift: {path}")


def execute(args: argparse.Namespace) -> Path:
    script = Path(__file__).resolve()
    script_dir = script.parent
    repo_root = script_dir.parent.parent
    reference_root = args.reference_root.resolve()
    output = args.output_dir.resolve()
    binary = args.current_driver.resolve()
    adapter = (
        script_dir
        / "sanisand_validation"
        / "pisano_udsm_material_point_driver.f90"
    )
    if output.exists():
        raise ValidationError(
            f"refusing to overwrite existing validation output: {output}"
        )
    output.parent.mkdir(parents=True, exist_ok=True)
    stage = Path(tempfile.mkdtemp(prefix=f".{output.name}-", dir=output.parent))
    try:
        build_current(repo_root, binary, args.skip_build)
        with tempfile.TemporaryDirectory(prefix="pisano-udsm-build-") as temp:
            reference_build = Path(temp)
            reference_record = compile_reference(reference_root, adapter, reference_build)
            current_csv = stage / "mpm_sanisand_validation.csv"
            reference_csv = stage / "pisano_udsm_validation.csv"
            current_result = run([str(binary)])
            current_csv.write_text(current_result.stdout, encoding="utf-8")
            reference_executable = Path(reference_record["executable"]["path"])
            reference_result = run([str(reference_executable)], cwd=reference_build)
            reference_csv.write_text(reference_result.stdout, encoding="utf-8")
            current = load_rows(current_csv, IMPLEMENTATION_CURRENT)
            reference = load_rows(reference_csv, IMPLEMENTATION_REFERENCE)
            metrics = compute_metrics(current, reference)
            if not metrics["verification_passed"]:
                raise ValidationError("predeclared numerical verification gates failed")
            cycles = cycle_rows(current + reference)
            metrics_csv = stage / "validation_metrics.csv"
            flat_metrics = [
                {
                    "metric": "monotonic_pressure_max_normalized",
                    "value": metrics["monotonic_code_to_code"]["pressure"]["maximum_absolute_normalized"],
                },
                {
                    "metric": "monotonic_q_max_normalized",
                    "value": metrics["monotonic_code_to_code"]["deviator_stress"]["maximum_absolute_normalized"],
                },
                {
                    "metric": "current_q_control_terminal_cycle",
                    "value": metrics["literature_q_controlled"][IMPLEMENTATION_CURRENT]["cycle_time"],
                },
                {
                    "metric": "reference_q_control_terminal_cycle",
                    "value": metrics["literature_q_controlled"][IMPLEMENTATION_REFERENCE]["cycle_time"],
                },
                {
                    "metric": "memory_surface_extension_cycles",
                    "value": metrics["literature_q_controlled"]["memory_surface_extension_cycles"],
                },
            ]
            write_csv(metrics_csv, flat_metrics)
            cycle_csv = stage / "cycle_metrics.csv"
            write_csv(cycle_csv, cycles)
            figures = plot_validation_paths(current, reference, stage)
            figures += plot_literature_path(current, reference, stage)
            summary_md = stage / "VALIDATION_RESULT.md"
            summary_md.write_text(markdown_summary(metrics), encoding="utf-8")
            # The temporary reference executable disappears after this block;
            # preserve its digest but do not claim it as a persistent artifact.
            reference_record["executable_sha256"] = reference_record.pop("executable")["sha256"]
            staged_persistent = [
                current_csv,
                reference_csv,
                metrics_csv,
                cycle_csv,
                summary_md,
                *figures,
            ]
            persistent_artifacts = [
                published_artifact(path, output / path.name)
                for path in staged_persistent
            ]
            source_artifacts = [
                artifact(script),
                artifact(
                    repo_root
                    / "mpm_hpc_source"
                    / "tests"
                    / "materials"
                    / "sanisand_validation_driver.cc"
                ),
                artifact(adapter),
                artifact(binary),
            ]
            audit = {
                "schema": SCHEMA,
                "generated_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
                "status": "verified_with_scope_limitation",
                "model_identity": {
                    "current": {
                        "name": "MPM SANISAND04-style fabric model",
                        "constitutive_state_count": 20,
                        "has_memory_surface": False,
                    },
                    "reference": {
                        "name": "Pisano SANISAND-MS UDSM",
                        "constitutive_state_count": 31,
                        "has_memory_surface": True,
                    },
                },
                "loading_protocol": {
                    "initial_pressure_pa": INITIAL_PRESSURE_PA,
                    "initial_void_ratio": 0.808,
                    "literature_q_amplitude_pa": Q_AMPLITUDE_PA,
                    "strain_controlled_cycles": 8,
                    "literature_q_controlled_cycles_requested": 4,
                    "steps_per_cycle": 400,
                    "display_smoothing_applied": False,
                },
                "metrics": metrics,
                "reference_implementation": reference_record,
                "primary_sources": [
                    {"title": "SANISAND-MS UDSM repository", "url": REFERENCE_URL},
                    {"title": "Liu et al. (2019)", "url": LIU_2019_DOI},
                    {"title": "Dafalias and Manzari (2004)", "url": DM04_DOI},
                ],
                "source_artifacts": source_artifacts,
                "artifacts": persistent_artifacts,
            }
            validate_artifacts(source_artifacts)
            audit_path = stage / "validation_audit.json"
            audit_path.write_text(
                json.dumps(audit, indent=2, sort_keys=True, allow_nan=False) + "\n",
                encoding="utf-8",
            )
        os.replace(stage, output)
        validate_artifacts(persistent_artifacts)
        return output / "validation_audit.json"
    except BaseException:
        shutil.rmtree(stage, ignore_errors=True)
        raise


def parser() -> argparse.ArgumentParser:
    script_dir = Path(__file__).resolve().parent
    repo_root = script_dir.parent.parent
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument(
        "--reference-root",
        type=Path,
        required=True,
        help="Local clone of FedericoPisano/SANISAND-MS-UDSM at the pinned commit",
    )
    result.add_argument(
        "--output-dir",
        type=Path,
        default=script_dir / "analysis" / "sanisand_validation" / "run_205c13b0",
    )
    result.add_argument(
        "--current-driver",
        type=Path,
        default=repo_root
        / "mpm_hpc_source"
        / "build-pipeline"
        / "sanisand_validation_driver",
    )
    result.add_argument("--skip-build", action="store_true")
    return result


def main() -> int:
    try:
        audit = execute(parser().parse_args())
        print(json.dumps({"status": "ok", "audit": str(audit)}, indent=2))
        return 0
    except ValidationError as exc:
        print(f"SANISAND validation failed: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
