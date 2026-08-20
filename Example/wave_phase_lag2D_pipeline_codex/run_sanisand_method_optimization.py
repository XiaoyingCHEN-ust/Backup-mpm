#!/usr/bin/env python3
"""Audited 2x2 SANISAND parent-parameterization selection study.

The study separates the critical-state-line and elastic bulk-modulus
parameterizations.  A variant is promoted only when it improves virgin
loading and does not regress the two cyclic holdouts.  The upstream UDSM is
consumed only through the independently audited parent validation package.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import importlib.util
import json
import math
import os
import shutil
import sys
import tempfile
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


SCHEMA = "pipeline-sanisand-method-optimization-v1"
PARENT_SCHEMA = "pipeline-sanisand-independent-validation-v2"
MONOTONIC_IMPROVEMENT_FACTOR = 0.99
HOLDOUT_REGRESSION_LIMIT = 1.10

VARIANTS: tuple[tuple[str, str, tuple[str, ...]], ...] = (
    ("mapped_legacy", "mpm_sanisand04_style", ()),
    ("direct_csl", "mpm_sanisand04_direct_csl", ("--direct-csl",)),
    (
        "direct_elasticity",
        "mpm_sanisand04_direct_elasticity",
        ("--direct-elasticity",),
    ),
    (
        "direct_csl_and_elasticity",
        "mpm_sanisand04_liu2019_parent",
        ("--liu2019-parent",),
    ),
)


class MethodStudyError(RuntimeError):
    pass


def load_base(script_dir: Path) -> Any:
    path = script_dir / "run_sanisand_validation.py"
    spec = importlib.util.spec_from_file_location("sanisand_validation_base", path)
    if spec is None or spec.loader is None:
        raise MethodStudyError("cannot load parent validation module")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def load_parent(audit_path: Path, base: Any) -> tuple[dict[str, Any], Path]:
    try:
        audit = json.loads(audit_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise MethodStudyError(f"cannot read parent validation audit: {exc}") from exc
    if audit.get("schema") != PARENT_SCHEMA:
        raise MethodStudyError("parent validation schema is not current")
    if audit.get("status") != "verified_with_scope_limitation":
        raise MethodStudyError("parent validation is not verified")
    try:
        base.validate_artifacts(audit["artifacts"])
        base.validate_artifacts(audit["source_artifacts"])
    except (KeyError, base.ValidationError) as exc:
        raise MethodStudyError(f"parent validation artifact drift: {exc}") from exc
    references = [
        Path(record["path"])
        for record in audit["artifacts"]
        if Path(record["path"]).name == "pisano_udsm_validation.csv"
    ]
    if len(references) != 1:
        raise MethodStudyError("parent audit must bind one reference CSV")
    return audit, references[0]


def require_common_grid(
    current: list[dict[str, Any]],
    implementation: str,
    reference: list[dict[str, Any]],
    base: Any,
    path_id: str,
) -> None:
    for tolerance in (base.LOOSE_TOLERANCE, base.TIGHT_TOLERANCE):
        left = base.select(current, implementation, path_id, tolerance)
        right = base.select(
            reference, base.IMPLEMENTATION_REFERENCE, path_id, tolerance
        )
        if [row["step"] for row in left] != [row["step"] for row in right]:
            raise MethodStudyError(f"step grid drift for {implementation}/{path_id}")
        if not np.allclose(
            [row["axial_strain"] for row in left],
            [row["axial_strain"] for row in right],
            rtol=0.0,
            atol=2.0e-15,
        ):
            raise MethodStudyError(
                f"prescribed strain drift for {implementation}/{path_id}"
            )


def path_errors(
    current: list[dict[str, Any]],
    implementation: str,
    reference: list[dict[str, Any]],
    base: Any,
    path_id: str,
) -> dict[str, Any]:
    left = base.select(current, implementation, path_id, base.TIGHT_TOLERANCE)
    right = base.select(
        reference, base.IMPLEMENTATION_REFERENCE, path_id, base.TIGHT_TOLERANCE
    )
    return {
        "pressure": base.normalized_error(
            left, right, "p_effective_pa", base.INITIAL_PRESSURE_PA
        ),
        "deviator_stress": base.normalized_error(
            left, right, "q_signed_pa", base.Q_AMPLITUDE_PA
        ),
    }


def maximum_error(errors: dict[str, Any]) -> float:
    return max(
        errors["pressure"]["maximum_absolute_normalized"],
        errors["deviator_stress"]["maximum_absolute_normalized"],
    )


def safe_ratio(value: float, denominator: float) -> float:
    if not (math.isfinite(value) and math.isfinite(denominator)) or (
        denominator <= 0
    ):
        raise MethodStudyError("invalid positive denominator in selection metric")
    return value / denominator


def select_parameterization(records: dict[str, Any]) -> str:
    """Apply the frozen calibration/holdout contract to complete records."""
    baseline = records["mapped_legacy"]
    baseline_p_rmse = baseline["cyclic_holdout"]["pressure"]["rmse_normalized"]
    baseline_q_rmse = baseline["cyclic_holdout"]["deviator_stress"][
        "rmse_normalized"
    ]
    baseline_cycle_gap = baseline["q_control_terminal_cycle_gap"]
    eligible: list[str] = []
    for variant, _, _ in VARIANTS:
        record = records[variant]
        if variant == "mapped_legacy":
            record["selection"] = {
                "monotonic_improved": False,
                "cyclic_pressure_rmse_ratio_to_baseline": 1.0,
                "cyclic_q_rmse_ratio_to_baseline": 1.0,
                "terminal_cycle_gap_ratio_to_baseline": 1.0,
                "holdout_passed": True,
                "eligible_to_replace_baseline": False,
            }
            continue
        p_ratio = safe_ratio(
            record["cyclic_holdout"]["pressure"]["rmse_normalized"],
            baseline_p_rmse,
        )
        q_ratio = safe_ratio(
            record["cyclic_holdout"]["deviator_stress"]["rmse_normalized"],
            baseline_q_rmse,
        )
        terminal_ratio = safe_ratio(
            record["q_control_terminal_cycle_gap"], baseline_cycle_gap
        )
        monotonic_improved = (
            record["monotonic_maximum_normalized"]
            <= MONOTONIC_IMPROVEMENT_FACTOR
            * baseline["monotonic_maximum_normalized"]
        )
        holdout_passed = max(p_ratio, q_ratio, terminal_ratio) <= (
            HOLDOUT_REGRESSION_LIMIT + 1.0e-12
        )
        candidate = monotonic_improved and holdout_passed
        record["selection"] = {
            "monotonic_improved": monotonic_improved,
            "cyclic_pressure_rmse_ratio_to_baseline": p_ratio,
            "cyclic_q_rmse_ratio_to_baseline": q_ratio,
            "terminal_cycle_gap_ratio_to_baseline": terminal_ratio,
            "holdout_passed": holdout_passed,
            "eligible_to_replace_baseline": candidate,
        }
        if candidate:
            eligible.append(variant)
    return min(
        eligible,
        key=lambda name: records[name]["monotonic_maximum_normalized"],
        default="mapped_legacy",
    )


def compute_metrics(
    rows_by_variant: dict[str, list[dict[str, Any]]],
    reference: list[dict[str, Any]],
    base: Any,
) -> dict[str, Any]:
    reference_q = base.select(
        reference,
        base.IMPLEMENTATION_REFERENCE,
        base.PATH_LITERATURE,
        base.TIGHT_TOLERANCE,
    )
    reference_terminal = base.terminal_record(reference_q)
    records: dict[str, Any] = {}
    for variant, implementation, _ in VARIANTS:
        rows = rows_by_variant[variant]
        for path in (base.PATH_MONOTONIC, base.PATH_CYCLIC):
            require_common_grid(rows, implementation, reference, base, path)
        convergence = {
            path: base.convergence(rows, implementation, path)
            for path in (
                base.PATH_MONOTONIC,
                base.PATH_CYCLIC,
                base.PATH_LITERATURE,
            )
        }
        if not all(item["passed"] for item in convergence.values()):
            raise MethodStudyError(f"integration convergence failed for {variant}")
        monotonic = path_errors(
            rows, implementation, reference, base, base.PATH_MONOTONIC
        )
        cyclic = path_errors(rows, implementation, reference, base, base.PATH_CYCLIC)
        terminal = base.terminal_record(
            base.select(
                rows, implementation, base.PATH_LITERATURE, base.TIGHT_TOLERANCE
            )
        )
        records[variant] = {
            "implementation": implementation,
            "integration_convergence": convergence,
            "monotonic": monotonic,
            "monotonic_maximum_normalized": maximum_error(monotonic),
            "cyclic_holdout": cyclic,
            "q_control_terminal": terminal,
            "q_control_terminal_cycle_gap": abs(
                reference_terminal["cycle_time"] - terminal["cycle_time"]
            ),
        }

    selected = select_parameterization(records)
    return {
        "selection_contract": {
            "calibration_path": "virgin monotonic constant-volume triaxial",
            "holdout_paths": [
                "eight-cycle common-strain path",
                "Liu Appendix-II q-controlled path",
            ],
            "minimum_monotonic_improvement_fraction": 1.0
            - MONOTONIC_IMPROVEMENT_FACTOR,
            "maximum_holdout_regression_fraction": HOLDOUT_REGRESSION_LIMIT - 1.0,
            "display_smoothing_applied": False,
        },
        "reference_q_control_terminal": reference_terminal,
        "variants": records,
        "selected_parameterization": selected,
        "production_change_recommended": selected != "mapped_legacy",
    }


def metric_table(metrics: dict[str, Any]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for variant, _, _ in VARIANTS:
        record = metrics["variants"][variant]
        rows.append(
            {
                "variant": variant,
                "monotonic_p_rmse_percent": 100
                * record["monotonic"]["pressure"]["rmse_normalized"],
                "monotonic_p_max_percent": 100
                * record["monotonic"]["pressure"][
                    "maximum_absolute_normalized"
                ],
                "monotonic_q_rmse_percent": 100
                * record["monotonic"]["deviator_stress"]["rmse_normalized"],
                "monotonic_q_max_percent": 100
                * record["monotonic"]["deviator_stress"][
                    "maximum_absolute_normalized"
                ],
                "cyclic_p_rmse_percent": 100
                * record["cyclic_holdout"]["pressure"]["rmse_normalized"],
                "cyclic_q_rmse_percent": 100
                * record["cyclic_holdout"]["deviator_stress"][
                    "rmse_normalized"
                ],
                "q_control_terminal_cycle": record["q_control_terminal"][
                    "cycle_time"
                ],
                "q_control_terminal_ru": record["q_control_terminal"]["ru"],
                "holdout_passed": record["selection"]["holdout_passed"],
                "eligible_to_replace_baseline": record["selection"][
                    "eligible_to_replace_baseline"
                ],
                "selected": variant == metrics["selected_parameterization"],
            }
        )
    return rows


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def plot_metrics(metrics: dict[str, Any], output: Path) -> list[Path]:
    labels = [name for name, _, _ in VARIANTS]
    display = ["Mapped\nlegacy", "Direct\nCSL", "Direct\nelasticity", "Direct\nboth"]
    records = [metrics["variants"][name] for name in labels]
    colors = ["#1F4E79", "#5B8DB8", "#9CBAD3", "#D2E0EB"]
    plt.rcParams.update(
        {
            "font.size": 8.5,
            "axes.titlesize": 9,
            "axes.labelsize": 8.5,
            "legend.fontsize": 7.5,
            "xtick.labelsize": 7.5,
            "ytick.labelsize": 7.5,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "savefig.dpi": 300,
            "pdf.fonttype": 42,
        }
    )
    fig, axes = plt.subplots(2, 2, figsize=(7.15, 5.4), constrained_layout=True)
    x = np.arange(len(labels), dtype=float)
    width = 0.36

    ax = axes[0, 0]
    p = [
        100 * r["monotonic"]["pressure"]["maximum_absolute_normalized"]
        for r in records
    ]
    q = [
        100
        * r["monotonic"]["deviator_stress"]["maximum_absolute_normalized"]
        for r in records
    ]
    ax.bar(
        x - width / 2,
        p,
        width,
        color="#AFC9DE",
        edgecolor="#1F4E79",
        label="$p'$ max",
    )
    ax.bar(
        x + width / 2,
        q,
        width,
        color="#E69F00",
        edgecolor="#8A5C00",
        hatch="//",
        label="$q$ max",
    )
    ax.set_title("(a) Virgin-loading maximum error")
    ax.set_ylabel("Normalised error (%)")
    ax.set_xticks(x, display)
    ax.legend(frameon=False, ncol=2)
    ax.grid(axis="y", color="#D9D9D9", linewidth=0.5)

    ax = axes[0, 1]
    p = [
        100 * r["cyclic_holdout"]["pressure"]["rmse_normalized"]
        for r in records
    ]
    q = [
        100 * r["cyclic_holdout"]["deviator_stress"]["rmse_normalized"]
        for r in records
    ]
    ax.bar(
        x - width / 2,
        p,
        width,
        color="#AFC9DE",
        edgecolor="#1F4E79",
        label="$p'$ RMSE",
    )
    ax.bar(
        x + width / 2,
        q,
        width,
        color="#E69F00",
        edgecolor="#8A5C00",
        hatch="//",
        label="$q$ RMSE",
    )
    ax.set_title("(b) Eight-cycle holdout error")
    ax.set_ylabel("Normalised RMSE (%)")
    ax.set_xticks(x, display)
    ax.legend(frameon=False, ncol=2)
    ax.grid(axis="y", color="#D9D9D9", linewidth=0.5)

    ax = axes[1, 0]
    terminal = [r["q_control_terminal"]["cycle_time"] for r in records]
    ax.bar(x, terminal, color=colors, edgecolor="#1F2933")
    reference = metrics["reference_q_control_terminal"]["cycle_time"]
    ax.axhline(
        reference,
        color="#D97706",
        linestyle="--",
        linewidth=1.4,
        label="SANISAND-MS reference",
    )
    ax.set_title("(c) Prescribed-$q$ controller terminal")
    ax.set_ylabel("Cycle number")
    ax.set_xticks(x, display)
    ax.set_ylim(0.0, 1.12 * reference)
    ax.legend(frameon=False, loc="upper right", bbox_to_anchor=(1.0, 0.93))
    ax.grid(axis="y", color="#D9D9D9", linewidth=0.5)

    ax = axes[1, 1]
    p_ratio = [
        r["selection"]["cyclic_pressure_rmse_ratio_to_baseline"]
        for r in records
    ]
    terminal_ratio = [
        r["selection"]["terminal_cycle_gap_ratio_to_baseline"]
        for r in records
    ]
    ax.bar(
        x - width / 2,
        p_ratio,
        width,
        color="#AFC9DE",
        edgecolor="#1F4E79",
        label="Cyclic $p'$ RMSE",
    )
    ax.bar(
        x + width / 2,
        terminal_ratio,
        width,
        color="#E69F00",
        edgecolor="#8A5C00",
        hatch="//",
        label="Terminal-gap error",
    )
    ax.axhline(
        HOLDOUT_REGRESSION_LIMIT,
        color="#4B5563",
        linestyle=":",
        linewidth=1.2,
        label="10% regression limit",
    )
    ax.set_title("(d) Holdout regression ratios")
    ax.set_ylabel("Ratio to mapped baseline")
    ax.set_xticks(x, display)
    ax.set_ylim(
        0.0,
        1.12
        * max(
            HOLDOUT_REGRESSION_LIMIT,
            max(p_ratio),
            max(terminal_ratio),
        ),
    )
    ax.legend(frameon=False, ncol=2)
    ax.grid(axis="y", color="#D9D9D9", linewidth=0.5)

    png = output / "figure_sanisand_method_selection.png"
    pdf = output / "figure_sanisand_method_selection.pdf"
    fig.savefig(png, dpi=300)
    fig.savefig(pdf)
    plt.close(fig)
    return [png, pdf]


def result_markdown(metrics: dict[str, Any]) -> str:
    baseline = metrics["variants"]["mapped_legacy"]
    direct = metrics["variants"]["direct_elasticity"]
    return f"""# SANISAND parent-parameterization method study

The four variants were tested as a 2x2 factorial separation of the critical-
state-line and elastic bulk-modulus laws.  The mapped legacy formulation was
retained.  Direct elasticity reduced the virgin-loading maximum error from
{100*baseline['monotonic_maximum_normalized']:.3f}% to
{100*direct['monotonic_maximum_normalized']:.3f}%, but its cyclic-pressure
RMSE ratio was {direct['selection']['cyclic_pressure_rmse_ratio_to_baseline']:.3f}
and its prescribed-q terminal-gap ratio was
{direct['selection']['terminal_cycle_gap_ratio_to_baseline']:.3f}; these exceed
the predeclared 1.10 holdout limit.  No candidate both improved the calibration
path and passed both cyclic holdouts.  Therefore no production parameterization
change is recommended.

This is a model-selection diagnostic, not evidence that the mapped law is
SANISAND-MS.  The memory-surface limitation of the parent validation remains.
"""


def execute(args: argparse.Namespace) -> Path:
    script = Path(__file__).resolve()
    script_dir = script.parent
    repo_root = script_dir.parent.parent
    base = load_base(script_dir)
    audit_path = args.validation_audit.resolve()
    parent_audit, reference_csv = load_parent(audit_path, base)
    output = args.output_dir.resolve()
    if output.exists():
        raise MethodStudyError(f"refusing to overwrite method study: {output}")
    output.parent.mkdir(parents=True, exist_ok=True)
    stage = Path(tempfile.mkdtemp(prefix=f".{output.name}-", dir=output.parent))
    binary = args.current_driver.resolve()
    try:
        base.build_current(repo_root, binary, args.skip_build)
        reference = base.load_rows(reference_csv, base.IMPLEMENTATION_REFERENCE)
        rows_by_variant: dict[str, list[dict[str, Any]]] = {}
        raw_paths: list[Path] = []
        for variant, implementation, cli in VARIANTS:
            result = base.run([str(binary), *cli])
            raw = stage / f"{variant}.csv"
            raw.write_text(result.stdout, encoding="utf-8")
            rows_by_variant[variant] = base.load_rows(raw, implementation)
            raw_paths.append(raw)
        metrics = compute_metrics(rows_by_variant, reference, base)
        metric_rows = metric_table(metrics)
        metrics_csv = stage / "parameterization_metrics.csv"
        write_csv(metrics_csv, metric_rows)
        figures = plot_metrics(metrics, stage)
        result_md = stage / "METHOD_SELECTION_RESULT.md"
        result_md.write_text(result_markdown(metrics), encoding="utf-8")

        persistent = [*raw_paths, metrics_csv, result_md, *figures]
        artifacts = [
            base.published_artifact(path, output / path.name)
            for path in persistent
        ]
        sources = [
            base.artifact(script),
            base.artifact(script_dir / "run_sanisand_validation.py"),
            base.artifact(
                repo_root
                / "mpm_hpc_source"
                / "tests"
                / "materials"
                / "sanisand_validation_driver.cc"
            ),
            base.artifact(
                repo_root
                / "mpm_hpc_source"
                / "include"
                / "material"
                / "sanisand.h"
            ),
            base.artifact(
                repo_root
                / "mpm_hpc_source"
                / "include"
                / "material"
                / "sanisand.tcc"
            ),
            base.artifact(binary),
            base.artifact(audit_path),
            base.artifact(reference_csv),
        ]
        audit = {
            "schema": SCHEMA,
            "generated_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
            "status": "verified_no_production_change",
            "parent_validation": {
                "path": str(audit_path),
                "sha256": base.sha256(audit_path),
                "schema": parent_audit["schema"],
            },
            "metrics": metrics,
            "source_artifacts": sources,
            "artifacts": artifacts,
        }
        base.validate_artifacts(sources)
        audit_file = stage / "method_optimization_audit.json"
        audit_file.write_text(
            json.dumps(audit, indent=2, sort_keys=True, allow_nan=False) + "\n",
            encoding="utf-8",
        )
        os.replace(stage, output)
        base.validate_artifacts(artifacts)
        return output / audit_file.name
    except BaseException:
        shutil.rmtree(stage, ignore_errors=True)
        raise


def parser() -> argparse.ArgumentParser:
    script_dir = Path(__file__).resolve().parent
    repo_root = script_dir.parent.parent
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument(
        "--validation-audit",
        type=Path,
        default=script_dir
        / "analysis"
        / "sanisand_validation"
        / "run_205c13b0_v2"
        / "validation_audit.json",
    )
    result.add_argument(
        "--output-dir",
        type=Path,
        default=script_dir
        / "analysis"
        / "sanisand_validation"
        / "method_optimization_v1",
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
    except Exception as exc:
        print(f"SANISAND method study failed: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
