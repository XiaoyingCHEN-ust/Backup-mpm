#!/usr/bin/env python3
"""Synthesize pre-registered pipeline comparisons into manuscript evidence.

The script deliberately keeps three questions separate:

* RL--RE: does removing only the fitted fundamental phase lag reduce the
  hydraulic trigger and/or realised skeleton-stress loss?
* RL--RM: under the identical lagged pressure database, does SANISAND's cyclic
  state evolution produce a materially different stress/engineering path from
  Mohr--Coulomb?
* LS--HS: what is the fully coupled physical context when saturation changes
  phase, amplitude and drainage together?

It does not declare SANISAND "better" because of a chosen displacement order.
The defensible advantage is representation of cyclic state/fabric evolution;
whether that additional mechanism matters in this case is decided from the
observed paths and responses.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Any


EXPECTED_CODES = ("LS", "HS", "HM", "RL", "RM", "RE")
PRESSURE_RESPONSE_RELATIVE_TOLERANCE = 0.02
PHASE_REDUCTION_MINIMUM_DEG = 5.0
RESPONSE_MATERIALITY_RELATIVE = 0.05
CONSTITUTIVE_INITIAL_STATE_RELATIVE_TOLERANCE = 0.05
STATE_EVOLUTION_ABSOLUTE_TOLERANCE = 1.0e-8
TRANSFORM_INVARIANT_ABSOLUTE_TOLERANCE_PA = 1.0e-6

IF_PREFIX = "support_ROI_upward_seepage_IF_ge_1"
STRESS_PREFIX = "support_ROI_stress_loss_Rsigma_le_0p05"


def finite_float(value: Any, label: str) -> float:
    number = float(value)
    if not math.isfinite(number):
        raise ValueError(f"{label} is not finite")
    return number


def case_code(path: Path) -> str:
    stem = path.name.removesuffix("_summary.json")
    if "_" not in stem:
        raise ValueError(f"Cannot identify case code from {path.name}")
    return stem.split("_", 1)[1]


def load_summaries(directory: Path, *, require_full: bool) -> dict[str, dict[str, Any]]:
    summaries: dict[str, dict[str, Any]] = {}
    for path in sorted(directory.glob("*_summary.json")):
        code = case_code(path)
        if code in summaries:
            raise ValueError(f"Duplicate summary for case {code}")
        summaries[code] = json.loads(path.read_text(encoding="utf-8"))
    required = EXPECTED_CODES if require_full else ("LS", "HS", "HM")
    missing = [code for code in required if code not in summaries]
    if missing:
        raise FileNotFoundError(
            f"Missing manuscript analysis summaries in {directory}: {missing}"
        )
    return summaries


def particle_metric(summary: dict[str, Any], name: str) -> float:
    particles = summary["particles"]
    if name not in particles:
        raise KeyError(f"{summary['case']} has no particle metric {name}")
    return finite_float(particles[name], f"{summary['case']}:{name}")


def threshold_metric(summary: dict[str, Any], prefix: str, metric: str) -> float:
    return particle_metric(summary, f"{prefix}_{metric}")


def materially_higher(left: float, right: float) -> bool:
    scale = max(abs(left), abs(right), 1.0e-15)
    return left - right > RESPONSE_MATERIALITY_RELATIVE * scale


def materially_different(left: float, right: float) -> bool:
    scale = max(abs(left), abs(right), 1.0e-15)
    return abs(left - right) > RESPONSE_MATERIALITY_RELATIVE * scale


def direction(left: float, right: float) -> str:
    if materially_higher(left, right):
        return "left_higher"
    if materially_higher(right, left):
        return "right_higher"
    return "not_materially_different"


def phase_control_qa(
    lagged: dict[str, Any], erased: dict[str, Any], metadata_path: Path
) -> dict[str, Any]:
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    if metadata.get("schema") != "phase-erased-pressure-control-v1":
        raise ValueError("Unexpected phase-control metadata schema")
    qa = metadata["qa"]
    transform_mean_change = finite_float(
        qa["max_abs_mean_change_pa"], "phase transform mean change"
    )
    transform_amplitude_change = finite_float(
        qa["max_abs_fundamental_amplitude_change_pa"],
        "phase transform amplitude change",
    )
    point_hash_match = (
        metadata["source"]["points_sha256"]
        == metadata["output"]["points_sha256"]
    )

    lagged_phase = lagged["phase"]
    erased_phase = erased["phase"]
    if not lagged_phase.get("available") or not erased_phase.get("available"):
        raise ValueError("RL/RE phase fits are unavailable")
    common = sorted(
        set(lagged_phase["probes"]) & set(erased_phase["probes"])
    )
    if not common:
        raise ValueError("RL and RE have no common pressure probes")
    surface_scale = max(
        finite_float(lagged_phase["surface_amplitude_pa"], "RL surface amplitude"),
        finite_float(erased_phase["surface_amplitude_pa"], "RE surface amplitude"),
        1.0,
    )
    maximum_amplitude_relative_difference = 0.0
    maximum_mean_normalized_difference = 0.0
    probe_differences: dict[str, Any] = {}
    for probe in common:
        left = lagged_phase["probes"][probe]
        right = erased_phase["probes"][probe]
        left_amplitude = finite_float(left["amplitude_pa"], f"RL {probe} amplitude")
        right_amplitude = finite_float(right["amplitude_pa"], f"RE {probe} amplitude")
        amplitude_scale = max(abs(left_amplitude), abs(right_amplitude), 1.0)
        amplitude_difference = abs(left_amplitude - right_amplitude) / amplitude_scale
        mean_difference = abs(
            finite_float(left["mean_pressure_pa"], f"RL {probe} mean")
            - finite_float(right["mean_pressure_pa"], f"RE {probe} mean")
        ) / surface_scale
        maximum_amplitude_relative_difference = max(
            maximum_amplitude_relative_difference, amplitude_difference
        )
        maximum_mean_normalized_difference = max(
            maximum_mean_normalized_difference, mean_difference
        )
        probe_differences[probe] = {
            "amplitude_relative_difference": amplitude_difference,
            "mean_difference_normalized_by_surface_amplitude": mean_difference,
        }
    passed = bool(
        point_hash_match
        and transform_mean_change <= TRANSFORM_INVARIANT_ABSOLUTE_TOLERANCE_PA
        and transform_amplitude_change
        <= TRANSFORM_INVARIANT_ABSOLUTE_TOLERANCE_PA
        and maximum_amplitude_relative_difference
        <= PRESSURE_RESPONSE_RELATIVE_TOLERANCE
        and maximum_mean_normalized_difference
        <= PRESSURE_RESPONSE_RELATIVE_TOLERANCE
    )
    return {
        "passed": passed,
        "classification": metadata.get("classification"),
        "point_hash_match": point_hash_match,
        "transform_max_abs_mean_change_pa": transform_mean_change,
        "transform_max_abs_fundamental_amplitude_change_pa": (
            transform_amplitude_change
        ),
        "response_max_amplitude_relative_difference": (
            maximum_amplitude_relative_difference
        ),
        "response_max_mean_difference_normalized_by_surface_amplitude": (
            maximum_mean_normalized_difference
        ),
        "probe_differences": probe_differences,
    }


def phase_lag_contrast(lagged: dict[str, Any], erased: dict[str, Any]) -> dict[str, Any]:
    probes = ("crown", "shoulder", "invert")
    rows: dict[str, Any] = {}
    reductions: list[float] = []
    for probe in probes:
        left = lagged["phase"]["probes"].get(probe)
        right = erased["phase"]["probes"].get(probe)
        if left is None or right is None:
            continue
        lagged_abs = abs(finite_float(left["phase_lag_deg"], f"RL {probe} lag"))
        erased_abs = abs(finite_float(right["phase_lag_deg"], f"RE {probe} lag"))
        reduction = lagged_abs - erased_abs
        reductions.append(reduction)
        rows[probe] = {
            "lagged_abs_phase_lag_deg": lagged_abs,
            "phase_erased_abs_phase_lag_deg": erased_abs,
            "reduction_deg": reduction,
        }
    if not reductions:
        raise ValueError("No crown/shoulder/invert phase contrast is available")
    mean_reduction = sum(reductions) / len(reductions)
    return {
        "passed": mean_reduction >= PHASE_REDUCTION_MINIMUM_DEG,
        "mean_abs_phase_lag_reduction_deg": mean_reduction,
        "probes": rows,
    }


def phase_effect(lagged: dict[str, Any], erased: dict[str, Any]) -> dict[str, Any]:
    if_area_metric = "time_integrated_area_m2_s"
    stress_area_metric = "time_integrated_area_m2_s"
    lagged_if_area = threshold_metric(lagged, IF_PREFIX, if_area_metric)
    erased_if_area = threshold_metric(erased, IF_PREFIX, if_area_metric)
    lagged_stress_area = threshold_metric(lagged, STRESS_PREFIX, stress_area_metric)
    erased_stress_area = threshold_metric(erased, STRESS_PREFIX, stress_area_metric)
    lagged_max_if = particle_metric(lagged, "max_upward_seepage_IF")
    erased_max_if = particle_metric(erased, "max_upward_seepage_IF")
    lagged_min_rsigma = particle_metric(
        lagged, "minimum_vertical_stress_remaining_ratio"
    )
    erased_min_rsigma = particle_metric(
        erased, "minimum_vertical_stress_remaining_ratio"
    )

    hydraulic_direction = direction(lagged_if_area, erased_if_area)
    if hydraulic_direction == "left_higher" and lagged_max_if >= 1.0:
        hydraulic_status = "supported"
    elif hydraulic_direction == "right_higher":
        hydraulic_status = "opposite"
    else:
        hydraulic_status = "not_materially_resolved"

    stress_direction = direction(lagged_stress_area, erased_stress_area)
    if stress_direction == "left_higher" and lagged_min_rsigma <= 0.05:
        stress_status = "supported"
    elif stress_direction == "right_higher":
        stress_status = "opposite"
    else:
        stress_status = "not_materially_resolved"

    return {
        "hydraulic_trigger": {
            "status": hydraulic_status,
            "metric": f"{IF_PREFIX}_{if_area_metric}",
            "lagged": lagged_if_area,
            "phase_erased": erased_if_area,
            "direction": hydraulic_direction,
            "lagged_max_IF": lagged_max_if,
            "phase_erased_max_IF": erased_max_if,
        },
        "realised_skeleton_stress_loss": {
            "status": stress_status,
            "metric": f"{STRESS_PREFIX}_{stress_area_metric}",
            "lagged": lagged_stress_area,
            "phase_erased": erased_stress_area,
            "direction": stress_direction,
            "lagged_min_Rsigma": lagged_min_rsigma,
            "phase_erased_min_Rsigma": erased_min_rsigma,
        },
        "engineering_response": response_comparison(lagged, erased),
    }


def response_comparison(left: dict[str, Any], right: dict[str, Any]) -> dict[str, Any]:
    names = (
        "maximum_abs_vertical_displacement_over_D",
        "post_release_no_contact_fraction",
        "maximum_abs_rotation_rad",
    )
    metrics: dict[str, Any] = {}
    materially_separated = False
    for name in names:
        left_value = finite_float(left["pipeline"][name], f"{left['case']}:{name}")
        right_value = finite_float(right["pipeline"][name], f"{right['case']}:{name}")
        separated = materially_different(left_value, right_value)
        materially_separated = materially_separated or separated
        metrics[name] = {
            "left": left_value,
            "right": right_value,
            "direction": direction(left_value, right_value),
            "materially_different": separated,
        }
    return {"materially_separated": materially_separated, "metrics": metrics}


def read_probe_history(directory: Path, code: str) -> list[dict[str, float]]:
    paths = sorted(directory.glob(f"*_{code}_probe_history.csv"))
    if len(paths) != 1:
        raise FileNotFoundError(
            f"Expected one {code} probe history in {directory}, found {len(paths)}"
        )
    with paths[0].open(newline="", encoding="utf-8") as stream:
        rows = [
            {key: float(value) for key, value in row.items() if value != ""}
            for row in csv.DictReader(stream)
        ]
    if len(rows) < 2:
        raise ValueError(f"{code} probe history has fewer than two rows")
    return rows


def probe_path_metrics(rows: list[dict[str, float]], state_name: str) -> dict[str, Any]:
    probes = ("crown", "shoulder", "invert")
    output: dict[str, Any] = {}
    for probe in probes:
        p_key = f"{probe}_p_effective_pa"
        q_key = f"{probe}_q_pa"
        state_key = f"{probe}_{state_name}"
        selected = [row for row in rows if p_key in row and q_key in row]
        if len(selected) < 2:
            continue
        p = [row[p_key] for row in selected]
        q = [row[q_key] for row in selected]
        path_length = sum(
            math.hypot(p[i] - p[i - 1], q[i] - q[i - 1])
            for i in range(1, len(selected))
        )
        net_distance = math.hypot(p[-1] - p[0], q[-1] - q[0])
        state_values = [row[state_key] for row in selected if state_key in row]
        output[probe] = {
            "q_p_path_length_pa": path_length,
            "q_p_net_distance_pa": net_distance,
            "q_p_path_to_net_ratio": path_length / max(net_distance, 1.0e-15),
            "p_effective_min_pa": min(p),
            "p_effective_max_pa": max(p),
            "q_min_pa": min(q),
            "q_max_pa": max(q),
            "state_variable": state_name,
            "state_range": (
                max(state_values) - min(state_values) if state_values else None
            ),
            "state_max_abs": (
                max(abs(value) for value in state_values) if state_values else None
            ),
        }
    if not output:
        raise ValueError(f"No q-p' probe paths are available for {state_name}")
    return output


def constitutive_initial_state_qa(
    sanisand_rows: list[dict[str, float]],
    mohr_coulomb_rows: list[dict[str, float]],
) -> dict[str, Any]:
    probes = ("crown", "shoulder", "invert")
    left = sanisand_rows[0]
    right = mohr_coulomb_rows[0]
    left_time = finite_float(left.get("time_s"), "RL first-frame time")
    right_time = finite_float(right.get("time_s"), "RM first-frame time")
    time_match = math.isclose(left_time, right_time, rel_tol=0.0, abs_tol=1.0e-12)
    comparisons: dict[str, Any] = {}
    maximum_relative_difference = 0.0
    for probe in probes:
        for variable in ("p_effective_pa", "q_pa"):
            key = f"{probe}_{variable}"
            if key not in left or key not in right:
                raise KeyError(f"RL/RM first-frame initial-state QA is missing {key}")
            left_value = finite_float(left[key], f"RL first-frame {key}")
            right_value = finite_float(right[key], f"RM first-frame {key}")
            # The 100 Pa floor matches the minimum confinement used by the
            # R_sigma analysis and avoids magnifying harmless near-zero q noise.
            scale = max(abs(left_value), abs(right_value), 100.0)
            relative_difference = abs(left_value - right_value) / scale
            maximum_relative_difference = max(
                maximum_relative_difference, relative_difference
            )
            comparisons[key] = {
                "SANISAND": left_value,
                "Mohr_Coulomb": right_value,
                "relative_difference": relative_difference,
            }
    return {
        "passed": (
            time_match
            and maximum_relative_difference
            <= CONSTITUTIVE_INITIAL_STATE_RELATIVE_TOLERANCE
        ),
        "first_frame_time_match": time_match,
        "SANISAND_first_frame_time_s": left_time,
        "Mohr_Coulomb_first_frame_time_s": right_time,
        "maximum_relative_difference": maximum_relative_difference,
        "comparisons": comparisons,
        "note": (
            "RM passes through MC_EQ while RL resumes elastic EQ directly. A "
            "failed audit means the comparison includes handoff-induced initial-"
            "state differences and cannot isolate constitutive evolution cleanly."
        ),
    }


def constitutive_effect(
    directory: Path, sanisand: dict[str, Any], mohr_coulomb: dict[str, Any]
) -> dict[str, Any]:
    sani_rows = read_probe_history(directory, "RL")
    mc_rows = read_probe_history(directory, "RM")
    initial_state_qa = constitutive_initial_state_qa(sani_rows, mc_rows)
    sani_paths = probe_path_metrics(sani_rows, "eps_p_q")
    mc_paths = probe_path_metrics(mc_rows, "pdstrain")
    state_active = any(
        values["state_range"] is not None
        and values["state_range"] > STATE_EVOLUTION_ABSOLUTE_TOLERANCE
        for values in sani_paths.values()
    )

    stress_metrics = {
        "stress_loss_area_time": {
            "SANISAND": threshold_metric(
                sanisand, STRESS_PREFIX, "time_integrated_area_m2_s"
            ),
            "Mohr_Coulomb": threshold_metric(
                mohr_coulomb, STRESS_PREFIX, "time_integrated_area_m2_s"
            ),
        },
        "minimum_Rsigma": {
            "SANISAND": particle_metric(
                sanisand, "minimum_vertical_stress_remaining_ratio"
            ),
            "Mohr_Coulomb": particle_metric(
                mohr_coulomb, "minimum_vertical_stress_remaining_ratio"
            ),
        },
    }
    response = response_comparison(sanisand, mohr_coulomb)
    stress_separated = any(
        materially_different(values["SANISAND"], values["Mohr_Coulomb"])
        for values in stress_metrics.values()
    )
    consequence_observed = stress_separated or response["materially_separated"]
    if initial_state_qa["passed"] and state_active and consequence_observed:
        status = "mechanistic_advantage_observed"
    elif not initial_state_qa["passed"]:
        status = "initial_state_mismatch_confounds_constitutive_ablation"
    elif state_active:
        status = "cyclic_state_active_but_response_not_materially_separated"
    else:
        status = "field_evidence_insufficient"
    return {
        "status": status,
        "interpretation": (
            "SANISAND's advantage is the resolved cyclic state-dependent path, "
            "not a pre-assumed ordering of final displacement. Predictive "
            "superiority would require independent field or laboratory response data."
        ),
        "sanisand_cyclic_state_evolution_active": state_active,
        "material_response_consequence_observed": consequence_observed,
        "initial_state_QA": initial_state_qa,
        "stress_metrics": stress_metrics,
        "engineering_response": response,
        "SANISAND_probe_paths": sani_paths,
        "Mohr_Coulomb_probe_paths": mc_paths,
    }


def physical_context(low: dict[str, Any], high: dict[str, Any]) -> dict[str, Any]:
    return {
        "interpretation": (
            "LS--HS is a physical saturation contrast; phase, amplitude, fluid "
            "storage and drainage change together and it is not phase-only."
        ),
        "phase": {
            probe: {
                "LS_phase_lag_deg": low["phase"]["probes"][probe]["phase_lag_deg"],
                "HS_phase_lag_deg": high["phase"]["probes"][probe]["phase_lag_deg"],
                "LS_amplitude_ratio": low["phase"]["probes"][probe][
                    "amplitude_ratio_to_surface"
                ],
                "HS_amplitude_ratio": high["phase"]["probes"][probe][
                    "amplitude_ratio_to_surface"
                ],
            }
            for probe in ("crown", "shoulder", "invert")
            if probe in low["phase"]["probes"] and probe in high["phase"]["probes"]
        },
        "engineering_response": response_comparison(high, low),
    }


def synthesize(
    directory: Path, metadata_path: Path, *, require_full: bool = True
) -> dict[str, Any]:
    summaries = load_summaries(directory, require_full=require_full)
    constants = {
        "pressure_response_relative_tolerance": (
            PRESSURE_RESPONSE_RELATIVE_TOLERANCE
        ),
        "minimum_phase_lag_reduction_deg": PHASE_REDUCTION_MINIMUM_DEG,
        "response_materiality_relative": RESPONSE_MATERIALITY_RELATIVE,
        "constitutive_initial_state_relative_tolerance": (
            CONSTITUTIVE_INITIAL_STATE_RELATIVE_TOLERANCE
        ),
        "state_evolution_absolute_tolerance": STATE_EVOLUTION_ABSOLUTE_TOLERANCE,
        "transform_invariant_absolute_tolerance_pa": (
            TRANSFORM_INVARIANT_ABSOLUTE_TOLERANCE_PA
        ),
    }
    evidence: dict[str, Any] = {
        "schema": "pipeline-manuscript-evidence-v1",
        "pre_registered_thresholds": constants,
        "physical_LS_HS_context": physical_context(summaries["LS"], summaries["HS"]),
    }
    if require_full:
        qa = phase_control_qa(summaries["RL"], summaries["RE"], metadata_path)
        lag_contrast = phase_lag_contrast(summaries["RL"], summaries["RE"])
        phase = phase_effect(summaries["RL"], summaries["RE"])
        constitutive = constitutive_effect(
            directory, summaries["RL"], summaries["RM"]
        )
        evidence.update(
            {
                "phase_only_RL_RE": {
                    "classification": "one-way numerical counterfactual",
                    "pressure_control_QA": qa,
                    "phase_lag_contrast": lag_contrast,
                    **phase,
                },
                "constitutive_RL_RM": constitutive,
                "claim_gate": {
                    "phase_lag_hydraulic_trigger_supported": bool(
                        qa["passed"]
                        and lag_contrast["passed"]
                        and phase["hydraulic_trigger"]["status"] == "supported"
                    ),
                    "phase_lag_realised_liquefaction_supported": bool(
                        qa["passed"]
                        and lag_contrast["passed"]
                        and phase["realised_skeleton_stress_loss"]["status"]
                        == "supported"
                    ),
                    "sanisand_mechanistic_advantage_supported": (
                        constitutive["status"] == "mechanistic_advantage_observed"
                    ),
                },
            }
        )
    return evidence


def markdown_report(evidence: dict[str, Any]) -> str:
    lines = [
        "# Pipeline manuscript evidence audit",
        "",
        "This file is generated from the completed analysis outputs. It is an "
        "evidence gate, not a substitute for inspecting the histories and fields.",
        "",
    ]
    gate = evidence.get("claim_gate")
    if gate is None:
        lines.extend(
            [
                "Only the fully coupled physical cases were analysed; phase-only "
                "and matched-pressure constitutive claims remain unavailable.",
                "",
            ]
        )
        return "\n".join(lines) + "\n"
    phase = evidence["phase_only_RL_RE"]
    hydraulic = phase["hydraulic_trigger"]
    stress = phase["realised_skeleton_stress_loss"]
    constitutive = evidence["constitutive_RL_RM"]
    lines.extend(
        [
            "## Claim gates",
            "",
            f"- Phase-control pressure QA passed: `{phase['pressure_control_QA']['passed']}`.",
            f"- Fundamental phase lag was materially reduced: `{phase['phase_lag_contrast']['passed']}` "
            f"(mean reduction `{phase['phase_lag_contrast']['mean_abs_phase_lag_reduction_deg']:.4g} deg`).",
            f"- Lagged pressure increased the hydraulic trigger: "
            f"`{gate['phase_lag_hydraulic_trigger_supported']}` "
            f"(RL/RE support-zone IF area-time = `{hydraulic['lagged']:.6g}` / "
            f"`{hydraulic['phase_erased']:.6g} m2 s`).",
            f"- Lagged pressure increased realised skeleton-stress loss: "
            f"`{gate['phase_lag_realised_liquefaction_supported']}` "
            f"(RL/RE support-zone Rsigma area-time = `{stress['lagged']:.6g}` / "
            f"`{stress['phase_erased']:.6g} m2 s`).",
            f"- SANISAND mechanistic advantage was expressed in the response: "
            f"`{gate['sanisand_mechanistic_advantage_supported']}` "
            f"(status `{constitutive['status']}`; initial-state QA "
            f"`{constitutive['initial_state_QA']['passed']}`).",
            "",
            "## Permitted wording",
            "",
        ]
    )
    if gate["phase_lag_hydraulic_trigger_supported"]:
        lines.append(
            "The phase-erased counterfactual confirms that retaining the measured "
            "hydraulic phase lag increases the duration/extent of the critical "
            "upward seepage-force condition, despite matched fundamental pressure amplitude."
        )
    else:
        lines.append(
            "The phase-only comparison does not support claiming that phase lag "
            "increased the hydraulic trigger in this case. Report the null or opposite result."
        )
    lines.append("")
    if gate["phase_lag_realised_liquefaction_supported"]:
        lines.append(
            "The larger hydraulic trigger is accompanied by more extensive realised "
            "vertical skeleton-stress loss (Rsigma <= 0.05), linking phase lag to "
            "liquefaction rather than pressure timing alone."
        )
    else:
        lines.append(
            "A stronger hydraulic trigger was not shown to produce materially larger "
            "realised skeleton-stress loss; restrict the conclusion to hydraulic potential."
        )
    lines.append("")
    if gate["sanisand_mechanistic_advantage_supported"]:
        lines.append(
            "Under the identical lagged pressure history, SANISAND resolves cyclic "
            "state evolution and a materially different stress/engineering path from "
            "Mohr-Coulomb. This demonstrates a mechanistic modelling advantage, not "
            "universal superiority or a pre-assumed displacement ordering."
        )
    else:
        lines.append(
            "The matched-pressure comparison does not yet demonstrate a material field-"
            "scale consequence of SANISAND's additional cyclic state variables. Do not "
            "claim predictive superiority from model complexity alone."
        )
    lines.extend(
        [
            "",
            "RL--RE must always be labelled a one-way numerical counterfactual. "
            "Use IF >= 1 and Rsigma <= 0.05 as the primary criteria; ru is diagnostic only.",
            "",
        ]
    )
    return "\n".join(lines)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--analysis-dir", type=Path, required=True)
    parser.add_argument("--phase-metadata", type=Path, required=True)
    parser.add_argument("--output-json", type=Path, required=True)
    parser.add_argument("--output-markdown", type=Path, required=True)
    parser.add_argument(
        "--physical-only",
        action="store_true",
        help="write only LS--HS context without requiring replay cases",
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()
    evidence = synthesize(
        args.analysis_dir,
        args.phase_metadata,
        require_full=not args.physical_only,
    )
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_markdown.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(
        json.dumps(evidence, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    args.output_markdown.write_text(markdown_report(evidence), encoding="utf-8")
    print(f"Wrote manuscript evidence JSON: {args.output_json.resolve()}")
    print(f"Wrote manuscript evidence report: {args.output_markdown.resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
