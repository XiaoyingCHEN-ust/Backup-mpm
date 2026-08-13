from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any, Iterable

import prepare_study as study


PIPELINE_HISTORY_HEADER = (
    "step,time,fixed,center_x,center_y,displacement_x,displacement_y,"
    "velocity_x,velocity_y,angle,angular_velocity,contact_force_x,"
    "contact_force_y,contact_moment,total_force_x,total_force_y,total_moment,"
    "contacts,max_penetration\n"
)


def linear_equilibrium_config(
    *,
    uuid: str = "UNIT_LS_EQ",
    nsteps: int = 10,
    result_path: str = "results/test/",
    saturation: float = study.LOW_LAG_SATURATION,
    cell_size: float = 0.02,
    particle_spacing: float = 0.01,
) -> dict[str, Any]:
    """Return a compact-step config with the full registered QA metadata."""

    config = study.equilibrium_config(
        "UNIT",
        saturation,
        result_path,
        ".",
        cell_size,
        particle_spacing,
        0.12,
    )
    config["analysis"]["uuid"] = uuid
    config["analysis"]["nsteps"] = nsteps
    config["analysis"]["resume"].update(
        {"uuid": uuid, "step": nsteps, "nsteps": nsteps}
    )
    config["post_processing"]["output_steps"] = nsteps
    return config


def _effective_unit_weight(config: dict[str, Any]) -> float:
    soil, fluid = config["materials"]
    gas_density = (
        float(fluid["gas_molar_mass"])
        * float(soil["p_ref"])
        / float(fluid["gas_constant"])
        / (273.15 + float(soil["initial_temperature"]))
    )
    mixture_density = (1.0 - float(soil["porosity"])) * float(
        soil["density"]
    ) + float(soil["porosity"]) * (
        float(fluid["liquid_saturation"]) * float(fluid["density"])
        + float(fluid["gas_saturation"]) * gas_density
    )
    return (mixture_density - float(fluid["density"])) * study.GRAVITY


def _counted_rows(rows: Iterable[Iterable[float]]) -> str:
    materialized = [tuple(row) for row in rows]
    body = "\n".join(" ".join(f"{value:.17g}" for value in row) for row in materialized)
    return f"{len(materialized)}\n{body}\n"


def _resolve(root: Path, value: str | Path) -> Path:
    path = Path(value)
    return path if path.is_absolute() else root / path


def write_pipeline_history_fixture(
    config: dict[str, Any], result_directory: Path
) -> Path:
    """Write the exact solver-named, nonempty rigid-pipeline history file."""

    digits = len(str(int(config["analysis"]["nsteps"])))
    path = result_directory / f"pipeline-history{0:0{digits}d}.csv"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        PIPELINE_HISTORY_HEADER
        + "0,0,1,0.75,0.5,0,0,0,0,0,0,0,0,0,0,0,0,0,0\n",
        encoding="utf-8",
    )
    return path


def write_fake_hdf5_auditor(root: Path) -> Path:
    """Install an executable JSONL auditor that mirrors the adjacent test VTP."""

    path = root / "fake_hdf5_checkpoint_audit.py"
    path.write_text(
        r'''#!/usr/bin/env python3
import json
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

arguments = sys.argv[1:]
if "--hdf5" not in arguments or "--jsonl" not in arguments:
    print("fake auditor requires --hdf5 PATH --jsonl", file=sys.stderr)
    raise SystemExit(2)
hdf5 = Path(arguments[arguments.index("--hdf5") + 1]).resolve()
vtp = hdf5.with_name(hdf5.name.replace("particles", "particle")).with_suffix(".vtp")
document = ET.parse(vtp).getroot()
piece = document.find(".//Piece")
point_data = {
    element.attrib.get("Name"): element
    for element in piece.findall("./PointData/DataArray")
}

def values(element):
    return [float(item) for item in (element.text or "").split()]

identifiers = [int(item) for item in values(point_data["ids"])]
count = len(identifiers)
points = values(piece.find("./Points/DataArray"))

def rows(name, width):
    flat = values(point_data[name])
    return [flat[index * width:(index + 1) * width] for index in range(count)]

coordinates = [points[index * 3:(index + 1) * 3] for index in range(count)]
velocities = rows("velocities", 3)
displacements = rows("displacements", 3)
stresses = rows("stresses", 6)
porosities = rows("porosities", 1)
volumes = rows("volumes", 1)
records = []
for row, identifier in enumerate(identifiers):
    records.append({
        "record_type": "particle",
        "id": identifier,
        "coordinates": coordinates[row],
        "velocities": velocities[row],
        "displacements": displacements[row],
        "stresses": stresses[row],
        "porosity": porosities[row][0],
        "volume": volumes[row][0],
        "status": True,
    })
records.sort(key=lambda record: record["id"])
summary = {
    "schema": "mpm-hdf5-checkpoint-audit-v1",
    "passed": True,
    "hdf5": str(hdf5),
    "table": "table",
    "table_fields": 173,
    "particle_count": count,
    "minimum_id": 0,
    "maximum_id": count - 1,
    "ids_contiguous_unique": True,
    "formal_schema": True,
    "legacy_schema": False,
    "compared_fields": [
        "coordinates", "velocities", "displacements", "stresses",
        "porosities", "volumes", "status"
    ],
}
print(json.dumps(summary, sort_keys=True))
for record in records:
    print(json.dumps(record, sort_keys=True))
''',
        encoding="utf-8",
    )
    path.chmod(0o755)
    return path


def _write_vtp(
    path: Path,
    coordinates: list[tuple[float, float]],
    identifiers: list[int],
    stresses_by_id: list[list[float]],
    *,
    maximum_velocity: float = 0.0,
    include_mc_fields: bool = False,
    phi: float = math.radians(31.0),
    cohesion: float = 0.0,
) -> None:
    count = len(identifiers)
    ordered_points = [coordinates[identifier] for identifier in identifiers]
    ordered_stresses = [stresses_by_id[identifier] for identifier in identifiers]
    velocities = [(maximum_velocity if index == 0 else 0.0, 0.0, 0.0) for index in range(count)]

    def flat(rows: Iterable[Iterable[float]]) -> str:
        return " ".join(f"{value:.17g}" for row in rows for value in row)

    point_data = [
        '<DataArray type="Int64" Name="ids" format="ascii">'
        + " ".join(str(identifier) for identifier in identifiers)
        + "</DataArray>",
        '<DataArray type="Float64" Name="porosities" format="ascii">'
        + " ".join("0.485" for _ in identifiers)
        + "</DataArray>",
        '<DataArray type="Float64" Name="volumes" format="ascii">'
        + " ".join("0.0001" for _ in identifiers)
        + "</DataArray>",
        '<DataArray type="Float64" Name="displacements" NumberOfComponents="3" format="ascii">'
        + flat((0.0, 0.0, 0.0) for _ in identifiers)
        + "</DataArray>",
        '<DataArray type="Float64" Name="velocities" NumberOfComponents="3" format="ascii">'
        + flat(velocities)
        + "</DataArray>",
        '<DataArray type="Float64" Name="stresses" NumberOfComponents="6" format="ascii">'
        + flat(ordered_stresses)
        + "</DataArray>",
    ]
    if include_mc_fields:
        point_data.extend(
            [
                '<DataArray type="Float64" Name="phi" format="ascii">'
                + " ".join(f"{phi:.17g}" for _ in identifiers)
                + "</DataArray>",
                '<DataArray type="Float64" Name="cohesion" format="ascii">'
                + " ".join(f"{cohesion:.17g}" for _ in identifiers)
                + "</DataArray>",
            ]
        )
    document = f'''<?xml version="1.0"?>
<VTKFile type="PolyData" version="0.1" byte_order="LittleEndian">
  <PolyData>
    <Piece NumberOfPoints="{count}" NumberOfVerts="0" NumberOfLines="0" NumberOfStrips="0" NumberOfPolys="0">
      <PointData>{''.join(point_data)}</PointData>
      <Points><DataArray type="Float64" NumberOfComponents="3" format="ascii">{flat((x, y, 0.0) for x, y in ordered_points)}</DataArray></Points>
    </Piece>
  </PolyData>
</VTKFile>
'''
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(document, encoding="utf-8")


def write_linear_qa_fixture(
    root: Path,
    config: dict[str, Any],
    checkpoint: Path,
    *,
    row_spacing: float = 0.01,
    x_spacing: float = 0.01,
    stress_offset_pa: float = 0.0,
    tensile_quiet_per_row: int = 0,
    maximum_velocity: float = 0.0,
    shuffled_ids: bool = False,
) -> dict[str, Any]:
    """Write a complete two-row linear-EQ QA fixture and its VTP/HDF5 pair."""

    mesh_path = _resolve(root, config["mesh"]["mesh"])
    particle_path = _resolve(
        root, config["particles"][0]["generator"]["location"]
    )
    stress_name = config["mesh"].get("particles_stresses")
    if stress_name is None:
        state = (
            "LS"
            if math.isclose(
                float(config["materials"][1]["liquid_saturation"]),
                study.LOW_LAG_SATURATION,
                rel_tol=0.0,
                abs_tol=1.0e-12,
            )
            else "HS"
        )
        stress_path = particle_path.parent / f"initial_effective_stresses_{state}.txt"
    else:
        stress_path = _resolve(root, stress_name)
    entity_sets_path = _resolve(root, config["mesh"]["entity_sets"])
    for path in (mesh_path, particle_path, stress_path, entity_sets_path):
        path.parent.mkdir(parents=True, exist_ok=True)

    mesh_path.write_text(
        "4 1\n0 0\n1.5 0\n1.5 0.5\n0 0.5\n0 1 2 3\n",
        encoding="utf-8",
    )
    seabed = float(config["materials"][1]["sea_level"]) - float(
        config["materials"][1]["depth_left"]
    )
    elevations = [seabed - 1.5 * row_spacing, seabed - 0.5 * row_spacing]
    x_count = int(round(1.5 / x_spacing))
    x_values = [(index + 0.5) * x_spacing for index in range(x_count)]
    coordinates = [(x, y) for y in elevations for x in x_values]
    particle_path.write_text(_counted_rows(coordinates), encoding="utf-8")

    lower_ids = list(range(x_count))
    top_ids = list(range(x_count, 2 * x_count))
    entity_sets_path.write_text(
        json.dumps(
            {
                "node_sets": [{"id": 0, "set": [0]}],
                "particle_sets": [
                    {"id": 0, "set": lower_ids + top_ids},
                    {"id": 4, "set": top_ids},
                ],
            }
        ),
        encoding="utf-8",
    )

    unit_weight = _effective_unit_weight(config)
    initial_stresses: list[list[float]] = []
    final_stresses: list[list[float]] = []
    for _, y in coordinates:
        sigma_yy = -unit_weight * max(seabed - y, 0.0)
        sigma_xx = study.INITIAL_EFFECTIVE_K0 * sigma_yy
        initial = [sigma_xx, sigma_yy, sigma_xx, 0.0, 0.0, 0.0]
        initial_stresses.append(initial)
        final = initial.copy()
        final[1] += stress_offset_pa
        final_stresses.append(final)
    stress_path.write_text(_counted_rows(initial_stresses), encoding="utf-8")

    cell_size = float(config["mesh"]["cellsize_min"])
    pipe = config["analysis"]["rigid_pipeline"]
    side_buffer = 2.0 * cell_size
    exclusion = float(pipe["outer_radius"]) + 2.0 * cell_size
    quiet_by_row: list[list[int]] = []
    for row_ids in (lower_ids, top_ids):
        quiet_by_row.append(
            [
                identifier
                for identifier in row_ids
                if coordinates[identifier][0] >= side_buffer - 1.0e-12
                and coordinates[identifier][0] <= 1.5 - side_buffer + 1.0e-12
                and abs(coordinates[identifier][0] - float(pipe["center"][0]))
                >= exclusion - 1.0e-12
            ]
        )
    for quiet_ids in quiet_by_row:
        for identifier in quiet_ids[:tensile_quiet_per_row]:
            final_stresses[identifier][1] = 1.0

    identifiers = list(range(len(coordinates)))
    if shuffled_ids:
        identifiers = identifiers[::2] + identifiers[1::2]
        identifiers.reverse()
    vtp = checkpoint.with_name(checkpoint.name.replace("particles", "particle")).with_suffix(
        ".vtp"
    )
    _write_vtp(
        vtp,
        coordinates,
        identifiers,
        final_stresses,
        maximum_velocity=maximum_velocity,
    )
    checkpoint.parent.mkdir(parents=True, exist_ok=True)
    checkpoint.write_bytes(b"\x89HDF\r\n\x1a\nqa-test-checkpoint")
    write_pipeline_history_fixture(config, checkpoint.parent)
    return {
        "particle_count": len(coordinates),
        "quiet_ids_by_row": quiet_by_row,
        "vtp": vtp,
        "checkpoint": checkpoint,
    }


def write_mc_qa_fixture(
    root: Path,
    config: dict[str, Any],
    checkpoint: Path,
    *,
    tensile_stress_pa: float = 0.0,
    maximum_velocity: float = 0.0,
) -> dict[str, Any]:
    """Write a minimal full-field Mohr-Coulomb handoff QA fixture."""

    # Reuse the linear input writer for a consistent particle/mesh inventory,
    # then replace the VTP with MC fields and the requested trial stress.
    fixture = write_linear_qa_fixture(root, config, checkpoint)
    particle_path = _resolve(
        root, config["particles"][0]["generator"]["location"]
    )
    rows = particle_path.read_text(encoding="utf-8").splitlines()[1:]
    coordinates = [tuple(float(value) for value in row.split()) for row in rows]
    count = len(coordinates)
    stresses = [
        [tensile_stress_pa, tensile_stress_pa, tensile_stress_pa, 0.0, 0.0, 0.0]
        for _ in range(count)
    ]
    _write_vtp(
        fixture["vtp"],
        coordinates,
        list(range(count)),
        stresses,
        maximum_velocity=maximum_velocity,
        include_mc_fields=True,
    )
    return fixture
