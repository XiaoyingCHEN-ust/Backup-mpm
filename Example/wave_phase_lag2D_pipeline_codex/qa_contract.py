#!/usr/bin/env python3
"""Versioned completion and static-stability QA contract.

This module is intentionally light-weight.  The case validator performs the
numerical checks, while the submitter and manuscript reducer use the routines
here to reject stale or structurally incomplete audit records.
"""

from __future__ import annotations

import math
import re
from typing import Any


COMPLETION_SCHEMA = "pipeline-case-completion-v3"
STABILITY_QA_SCHEMA = "pipeline-stability-qa-v1"
HDF5_VTP_AUDIT_SCHEMA = "mpm-hdf5-vtp-checkpoint-audit-v1"
HDF5_STRUCTURE_AUDIT_SCHEMA = "mpm-hdf5-checkpoint-audit-v1"

HDF5_VTP_COMPARED_FIELDS = (
    "coordinates",
    "velocities",
    "displacements",
    "stresses",
    "porosities",
    "volumes",
)
HDF5_VTP_ABSOLUTE_TOLERANCES = {
    "coordinates": 1.0e-7,
    "velocities": 1.0e-7,
    "displacements": 1.0e-7,
    "stresses": 1.0e-3,
    "porosities": 1.0e-7,
    "volumes": 1.0e-10,
}
HDF5_STRUCTURE_COMPARED_FIELDS = (*HDF5_VTP_COMPARED_FIELDS, "status")
HDF5_FORMAL_TABLE_FIELDS = 173

LINEAR_EQUILIBRIUM_MODE = "linear-elastic-equilibrium"
MC_HANDOFF_MODE = "mohr-coulomb-handoff"
STABILITY_MODES = {LINEAR_EQUILIBRIUM_MODE, MC_HANDOFF_MODE}

MAXIMUM_VELOCITY_M_S = 1.0e-3
MAXIMUM_DISPLACEMENT_M = 2.0e-2
POROSITY_MIN_EXCLUSIVE = 0.0
POROSITY_MAX_EXCLUSIVE = 1.0

NORMAL_STRESS_QUIET_SIDE_BUFFER_CELLS = 2.0
NORMAL_STRESS_QUIET_PIPE_BUFFER_CELLS = 2.0
NORMAL_STRESS_MEDIAN_DRIFT_CELL_WEIGHTS = 1.0
NORMAL_STRESS_MAXIMUM_TENSILE_FRACTION = 0.05
NORMAL_STRESS_MINIMUM_PARTICLES_PER_ROW = 10

MC_MAXIMUM_POSITIVE_YIELD_RESIDUAL_PA = 1.0e-6

_SHA256 = re.compile(r"[0-9a-f]{64}")


def stability_contract(mode: str) -> dict[str, Any]:
    """Return the exact config/manifest contract for one static stage."""

    if mode not in STABILITY_MODES:
        raise ValueError(f"Unknown stability QA mode: {mode}")
    limits: dict[str, Any] = {
        "maximum_velocity_m_s": MAXIMUM_VELOCITY_M_S,
        "maximum_displacement_m": MAXIMUM_DISPLACEMENT_M,
        "porosity_min_exclusive": POROSITY_MIN_EXCLUSIVE,
        "porosity_max_exclusive": POROSITY_MAX_EXCLUSIVE,
    }
    if mode == LINEAR_EQUILIBRIUM_MODE:
        limits.update(
            {
                "quiet_side_buffer_cells": NORMAL_STRESS_QUIET_SIDE_BUFFER_CELLS,
                "quiet_pipe_buffer_cells": NORMAL_STRESS_QUIET_PIPE_BUFFER_CELLS,
                "median_drift_cell_weights": (
                    NORMAL_STRESS_MEDIAN_DRIFT_CELL_WEIGHTS
                ),
                "maximum_tensile_fraction": (
                    NORMAL_STRESS_MAXIMUM_TENSILE_FRACTION
                ),
            }
        )
    else:
        limits["maximum_positive_yield_residual_pa"] = (
            MC_MAXIMUM_POSITIVE_YIELD_RESIDUAL_PA
        )
    return {"schema": STABILITY_QA_SCHEMA, "mode": mode, "limits": limits}


def manifest_constants() -> dict[str, Any]:
    """Return the QA constants recorded in every registered tier manifest."""

    return {
        "completion_schema": COMPLETION_SCHEMA,
        "stability_qa_schema": STABILITY_QA_SCHEMA,
        "linear_equilibrium": stability_contract(LINEAR_EQUILIBRIUM_MODE),
        "mohr_coulomb_handoff": stability_contract(MC_HANDOFF_MODE),
        "normal_stress_tolerance_definition": (
            "gamma_effective_times_background_cell_size"
        ),
        "normal_stress_selection_basis": "initial_particle_ids_and_coordinates",
    }


def _finite_number(value: Any, description: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{description} is not numeric")
    number = float(value)
    if not math.isfinite(number):
        raise ValueError(f"{description} is not finite")
    return number


def validate_artifact_record(value: Any, description: str) -> dict[str, Any]:
    """Validate only the shape of a hash-bound artifact record."""

    if not isinstance(value, dict) or set(value) != {
        "path",
        "size_bytes",
        "sha256",
    }:
        raise ValueError(f"{description} artifact audit is malformed")
    if not isinstance(value["path"], str) or not value["path"]:
        raise ValueError(f"{description} artifact path is malformed")
    size = value["size_bytes"]
    if isinstance(size, bool) or not isinstance(size, int) or size <= 0:
        raise ValueError(f"{description} artifact size is malformed")
    digest = value["sha256"]
    if not isinstance(digest, str) or _SHA256.fullmatch(digest) is None:
        raise ValueError(f"{description} artifact hash is malformed")
    return value


def validate_hdf5_vtp_audit(
    value: Any,
    *,
    expected_hdf5: str | None = None,
    expected_vtp: str | None = None,
) -> dict[str, Any]:
    """Require a complete, passing formal HDF5-to-VTP cross-format audit."""

    expected_keys = {
        "schema",
        "passed",
        "hdf5",
        "vtp",
        "particle_count_hdf5",
        "particle_count_vtp",
        "ids_sha256_hdf5",
        "ids_sha256_vtp",
        "compared_fields",
        "maximum_absolute_differences",
        "tolerances",
        "hdf5_structure",
    }
    if not isinstance(value, dict) or set(value) != expected_keys:
        raise ValueError("HDF5/VTP cross-format audit is malformed")
    if value.get("schema") != HDF5_VTP_AUDIT_SCHEMA or value.get("passed") is not True:
        raise ValueError("HDF5/VTP cross-format audit did not formally pass")
    hdf5 = value.get("hdf5")
    vtp = value.get("vtp")
    if not isinstance(hdf5, str) or not hdf5 or not isinstance(vtp, str) or not vtp:
        raise ValueError("HDF5/VTP cross-format audit paths are malformed")
    if expected_hdf5 is not None and hdf5 != expected_hdf5:
        raise ValueError("HDF5/VTP audit references another HDF5 checkpoint")
    if expected_vtp is not None and vtp != expected_vtp:
        raise ValueError("HDF5/VTP audit references another VTP checkpoint")

    count_hdf5 = value.get("particle_count_hdf5")
    count_vtp = value.get("particle_count_vtp")
    if (
        isinstance(count_hdf5, bool)
        or not isinstance(count_hdf5, int)
        or count_hdf5 <= 0
        or isinstance(count_vtp, bool)
        or not isinstance(count_vtp, int)
        or count_vtp != count_hdf5
    ):
        raise ValueError("HDF5/VTP particle counts are invalid or unequal")
    ids_hdf5 = value.get("ids_sha256_hdf5")
    ids_vtp = value.get("ids_sha256_vtp")
    if (
        not isinstance(ids_hdf5, str)
        or _SHA256.fullmatch(ids_hdf5) is None
        or not isinstance(ids_vtp, str)
        or ids_vtp != ids_hdf5
    ):
        raise ValueError("HDF5/VTP particle ID hashes are invalid or unequal")
    fields = list(HDF5_VTP_COMPARED_FIELDS)
    if value.get("compared_fields") != fields:
        raise ValueError("HDF5/VTP compared-field contract is stale")
    differences = value.get("maximum_absolute_differences")
    tolerances = value.get("tolerances")
    if not isinstance(differences, dict) or set(differences) != set(fields):
        raise ValueError("HDF5/VTP difference audit is malformed")
    if tolerances != HDF5_VTP_ABSOLUTE_TOLERANCES:
        raise ValueError("HDF5/VTP tolerance contract is stale")
    for field in fields:
        difference = _finite_number(
            differences[field], f"HDF5/VTP {field} maximum difference"
        )
        if difference < 0.0 or difference > HDF5_VTP_ABSOLUTE_TOLERANCES[field]:
            raise ValueError(f"HDF5/VTP {field} difference exceeds its tolerance")

    structure = value.get("hdf5_structure")
    structure_keys = {
        "schema",
        "passed",
        "hdf5",
        "table",
        "table_fields",
        "particle_count",
        "minimum_id",
        "maximum_id",
        "ids_contiguous_unique",
        "formal_schema",
        "legacy_schema",
        "compared_fields",
    }
    if not isinstance(structure, dict) or set(structure) != structure_keys:
        raise ValueError("HDF5 structural audit is malformed")
    if (
        structure.get("schema") != HDF5_STRUCTURE_AUDIT_SCHEMA
        or structure.get("passed") is not True
        or structure.get("hdf5") != hdf5
        or structure.get("table") != "table"
        or structure.get("table_fields") != HDF5_FORMAL_TABLE_FIELDS
        or structure.get("particle_count") != count_hdf5
        or structure.get("minimum_id") != 0
        or structure.get("maximum_id") != count_hdf5 - 1
        or structure.get("ids_contiguous_unique") is not True
        or structure.get("formal_schema") is not True
        or structure.get("legacy_schema") is not False
        or structure.get("compared_fields")
        != list(HDF5_STRUCTURE_COMPARED_FIELDS)
    ):
        raise ValueError("HDF5 structural audit does not meet the formal contract")
    return value


def validate_stability_qa(
    value: Any, *, expected_mode: str | None = None
) -> dict[str, Any]:
    """Require a complete, passing, versioned static-stability QA record."""

    if not isinstance(value, dict):
        raise ValueError("Stability QA record is not an object")
    mode = value.get("mode")
    if mode not in STABILITY_MODES:
        raise ValueError("Stability QA mode is invalid")
    if expected_mode is not None and mode != expected_mode:
        raise ValueError(
            f"Stability QA mode {mode!r} differs from {expected_mode!r}"
        )
    detail_name = (
        "normal_stress"
        if mode == LINEAR_EQUILIBRIUM_MODE
        else "mc_feasibility"
    )
    if set(value) != {
        "schema",
        "mode",
        "limits",
        "observed",
        detail_name,
        "summary",
    }:
        raise ValueError("Stability QA fields are missing or unexpected")
    if value.get("schema") != STABILITY_QA_SCHEMA:
        raise ValueError("Stability QA schema is stale")
    expected_limits = stability_contract(mode)["limits"]
    if value.get("limits") != expected_limits:
        raise ValueError("Stability QA limits are stale")
    if not isinstance(value.get("summary"), str) or not value["summary"].strip():
        raise ValueError("Stability QA summary is absent")

    observed = value.get("observed")
    if not isinstance(observed, dict) or set(observed) != {
        "maximum_velocity_m_s",
        "maximum_displacement_m",
        "porosity_min",
        "porosity_max",
        "particle_count",
    }:
        raise ValueError("Stability QA observations are malformed")
    maximum_velocity = _finite_number(
        observed["maximum_velocity_m_s"], "maximum velocity"
    )
    maximum_displacement = _finite_number(
        observed["maximum_displacement_m"], "maximum displacement"
    )
    porosity_min = _finite_number(observed["porosity_min"], "minimum porosity")
    porosity_max = _finite_number(observed["porosity_max"], "maximum porosity")
    particle_count = observed["particle_count"]
    if (
        maximum_velocity < 0.0
        or maximum_velocity > MAXIMUM_VELOCITY_M_S
        or maximum_displacement < 0.0
        or maximum_displacement > MAXIMUM_DISPLACEMENT_M
        or not POROSITY_MIN_EXCLUSIVE < porosity_min <= porosity_max
        or not porosity_max < POROSITY_MAX_EXCLUSIVE
        or isinstance(particle_count, bool)
        or not isinstance(particle_count, int)
        or particle_count <= 0
    ):
        raise ValueError("Stability QA observations do not pass registered limits")

    if mode == LINEAR_EQUILIBRIUM_MODE:
        detail = value["normal_stress"]
        if not isinstance(detail, dict) or set(detail) != {
            "inputs",
            "particle_id_basis",
            "background_cell_size_m",
            "effective_unit_weight_n_m3",
            "median_drift_tolerance_pa",
            "selection",
            "rows",
        }:
            raise ValueError("Normal-stress QA fields are malformed")
        inputs = detail["inputs"]
        if not isinstance(inputs, dict) or set(inputs) != {
            "mesh",
            "particles",
            "initial_stresses",
            "entity_sets",
        }:
            raise ValueError("Normal-stress QA input audits are malformed")
        for name, record in inputs.items():
            validate_artifact_record(record, f"normal-stress {name}")
        if detail["particle_id_basis"] != "initial-particle-file":
            raise ValueError("Normal-stress QA does not use initial particle IDs")
        cell_size = _finite_number(
            detail["background_cell_size_m"], "background cell size"
        )
        unit_weight = _finite_number(
            detail["effective_unit_weight_n_m3"], "effective unit weight"
        )
        tolerance = _finite_number(
            detail["median_drift_tolerance_pa"], "normal-stress tolerance"
        )
        if cell_size <= 0.0 or unit_weight <= 0.0 or not math.isclose(
            tolerance,
            unit_weight * cell_size * NORMAL_STRESS_MEDIAN_DRIFT_CELL_WEIGHTS,
            rel_tol=1.0e-12,
            abs_tol=1.0e-12,
        ):
            raise ValueError("Normal-stress QA tolerance is inconsistent")
        selection = detail["selection"]
        if not isinstance(selection, dict) or set(selection) != {
            "initial_surface_pset_id",
            "initial_top_row_pset_id",
            "side_buffer_m",
            "pipe_buffer_m",
            "domain_x_m",
            "pipe_center_x_m",
            "pipe_outer_radius_m",
        }:
            raise ValueError("Normal-stress quiet selection is malformed")
        if selection["initial_surface_pset_id"] != 0:
            raise ValueError("Normal-stress QA surface band is stale")
        if selection["initial_top_row_pset_id"] != 4:
            raise ValueError("Normal-stress QA top row is stale")
        side_buffer = _finite_number(selection["side_buffer_m"], "side buffer")
        pipe_buffer = _finite_number(selection["pipe_buffer_m"], "pipe buffer")
        if not math.isclose(
            side_buffer,
            NORMAL_STRESS_QUIET_SIDE_BUFFER_CELLS * cell_size,
            rel_tol=0.0,
            abs_tol=1.0e-12,
        ) or not math.isclose(
            pipe_buffer,
            NORMAL_STRESS_QUIET_PIPE_BUFFER_CELLS * cell_size,
            rel_tol=0.0,
            abs_tol=1.0e-12,
        ):
            raise ValueError("Normal-stress quiet buffers are stale")
        domain = selection["domain_x_m"]
        if (
            not isinstance(domain, list)
            or len(domain) != 2
            or _finite_number(domain[0], "domain lower x")
            >= _finite_number(domain[1], "domain upper x")
        ):
            raise ValueError("Normal-stress domain is malformed")
        _finite_number(selection["pipe_center_x_m"], "pipe centre x")
        if _finite_number(selection["pipe_outer_radius_m"], "pipe radius") <= 0.0:
            raise ValueError("Normal-stress pipe radius is invalid")
        rows = detail["rows"]
        if not isinstance(rows, list) or len(rows) != 2:
            raise ValueError("Normal-stress QA must contain the highest two rows")
        previous_y = -math.inf
        for row in rows:
            if not isinstance(row, dict) or set(row) != {
                "initial_y_m",
                "quiet_particle_count",
                "initial_sigma_yy_median_pa",
                "final_sigma_yy_median_pa",
                "median_drift_pa",
                "tensile_fraction",
            }:
                raise ValueError("Normal-stress row audit is malformed")
            y = _finite_number(row["initial_y_m"], "initial row elevation")
            if y <= previous_y:
                raise ValueError("Normal-stress rows are not strictly ordered")
            previous_y = y
            count = row["quiet_particle_count"]
            if (
                isinstance(count, bool)
                or not isinstance(count, int)
                or count < NORMAL_STRESS_MINIMUM_PARTICLES_PER_ROW
            ):
                raise ValueError("Normal-stress quiet row is undersampled")
            _finite_number(
                row["initial_sigma_yy_median_pa"], "initial median stress"
            )
            _finite_number(row["final_sigma_yy_median_pa"], "final median stress")
            drift = _finite_number(row["median_drift_pa"], "median stress drift")
            tensile = _finite_number(row["tensile_fraction"], "tensile fraction")
            if abs(drift) > tolerance or not (
                0.0 <= tensile <= NORMAL_STRESS_MAXIMUM_TENSILE_FRACTION
            ):
                raise ValueError("Normal-stress row does not pass registered limits")
    else:
        detail = value["mc_feasibility"]
        if not isinstance(detail, dict) or set(detail) != {
            "particle_count",
            "maximum_tension_residual_pa",
            "maximum_shear_residual_pa",
            "maximum_positive_residual_pa",
            "violating_particle_count",
        }:
            raise ValueError("Mohr-Coulomb feasibility QA is malformed")
        count = detail["particle_count"]
        violations = detail["violating_particle_count"]
        if (
            isinstance(count, bool)
            or not isinstance(count, int)
            or count != particle_count
            or isinstance(violations, bool)
            or not isinstance(violations, int)
            or violations != 0
        ):
            raise ValueError("Mohr-Coulomb feasibility counts are invalid")
        tension = _finite_number(
            detail["maximum_tension_residual_pa"], "maximum tension residual"
        )
        shear = _finite_number(
            detail["maximum_shear_residual_pa"], "maximum shear residual"
        )
        positive = _finite_number(
            detail["maximum_positive_residual_pa"], "maximum positive residual"
        )
        if not math.isclose(
            positive, max(0.0, tension, shear), rel_tol=1.0e-12, abs_tol=1.0e-12
        ) or positive > MC_MAXIMUM_POSITIVE_YIELD_RESIDUAL_PA:
            raise ValueError("Mohr-Coulomb yield residual exceeds the QA limit")
    return value
