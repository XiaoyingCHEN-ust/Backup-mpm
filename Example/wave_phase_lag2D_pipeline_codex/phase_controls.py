#!/usr/bin/env python3
"""Create an auditable phase-erased pressure replay database.

The solver pressure database stores liquid and gas pressure for every material
point at fixed output times.  This tool rotates only the fitted wave-frequency
component at each point to the phase of the local seabed-surface point.  The
point-wise mean, fundamental amplitude, residual/higher harmonics, and the
progressive-wave phase along x are retained.

The transformed database is a one-way numerical counterfactual.  It is not a
second fully coupled physical solution and must be labelled that way in the
paper.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import shutil
import struct
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import BinaryIO, Iterator

import numpy as np


MAGIC_V1 = b"MPM_PRESSURE_V1\0"
MAGIC_V2 = b"MPM_PRESSURE_V2\0"
HEADER = struct.Struct("<16sIIQQQd")
FRAME_STEP = struct.Struct("<Q")
PHASE_ALIGNMENT_RESIDUAL_TOLERANCE_DEG = 1.0e-10


@dataclass(frozen=True)
class PressureHeader:
    dimension: int
    values_per_particle: int
    particle_count: int
    step_interval: int
    max_step: int
    source_dt: float
    format_version: str = "V2"

    @property
    def frame_count(self) -> int:
        return self.max_step // self.step_interval + 1

    @property
    def frame_bytes(self) -> int:
        return FRAME_STEP.size + 2 * self.particle_count * 8

    @property
    def sample_dt(self) -> float:
        return self.source_dt * self.step_interval

    def physical_time(self, step: int) -> float:
        """Return the physical sample time, including the legacy V1 offset."""
        legacy_offset = 1 if self.format_version == "V1" else 0
        return (step + legacy_offset) * self.source_dt


@dataclass(frozen=True)
class PressurePoints:
    particle_ids: np.ndarray
    coordinates: np.ndarray


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _temporary_file(path: Path) -> tuple[int, Path]:
    """Create an exclusive temporary file beside its eventual target."""

    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, name = tempfile.mkstemp(
        dir=path.parent, prefix=f".{path.name}.", suffix=".tmp"
    )
    return descriptor, Path(name)


def _remove_temporary(path: Path | None) -> None:
    if path is None:
        return
    try:
        path.unlink()
    except FileNotFoundError:
        pass


def _stage_json(path: Path, value: dict[str, object]) -> Path:
    descriptor, temporary = _temporary_file(path)
    try:
        mode = path.stat().st_mode & 0o7777 if path.exists() else 0o644
        os.fchmod(descriptor, mode)
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="\n") as stream:
            json.dump(value, stream, indent=2, sort_keys=True, allow_nan=False)
            stream.write("\n")
            stream.flush()
            os.fsync(stream.fileno())
    except BaseException:
        try:
            os.close(descriptor)
        except OSError:
            pass
        _remove_temporary(temporary)
        raise
    return temporary


def _stage_file_copy(source: Path, target: Path) -> Path:
    """Copy and fsync a file into a same-directory publication candidate."""

    descriptor, temporary = _temporary_file(target)
    try:
        os.fchmod(descriptor, source.stat().st_mode & 0o7777)
        with source.open("rb") as input_stream, os.fdopen(
            descriptor, "wb"
        ) as output_stream:
            shutil.copyfileobj(input_stream, output_stream)
            output_stream.flush()
            os.fsync(output_stream.fileno())
    except BaseException:
        try:
            os.close(descriptor)
        except OSError:
            pass
        _remove_temporary(temporary)
        raise
    return temporary


def atomic_write_json(path: Path, value: dict[str, object]) -> None:
    """Atomically replace a finite JSON document."""

    temporary = _stage_json(path, value)
    try:
        os.replace(temporary, path)
    finally:
        _remove_temporary(temporary)


def _paths_alias(first: Path, second: Path) -> bool:
    """Return whether two path spellings identify the same filesystem object."""

    if first.resolve() == second.resolve():
        return True
    try:
        return os.path.samefile(first, second)
    except (FileNotFoundError, NotADirectoryError):
        return False


def _validate_distinct_output_targets(
    source_paths: tuple[Path, ...], output_paths: tuple[Path, ...]
) -> None:
    for output_path in output_paths:
        for source_path in source_paths:
            if _paths_alias(output_path, source_path):
                raise ValueError(
                    "Refusing to overwrite the source pressure database or "
                    "any source alias"
                )
    for index, output_path in enumerate(output_paths):
        for other in output_paths[index + 1 :]:
            if _paths_alias(output_path, other):
                raise ValueError("Pressure-database output targets alias each other")


def _publish_transform_files(
    *,
    staged_points: Path,
    staged_values: Path,
    output_points: Path,
    output_values: Path,
    metadata_path: Path,
    metadata: dict[str, object],
) -> None:
    """Publish values, points, then metadata while preserving an old trio.

    Both data candidates must already be complete and fsynced.  Metadata is
    serialized and fsynced before any target changes, but is replaced last so
    readers never see new metadata naming only partly published data files.
    """

    metadata_temporary: Path | None = None
    backups: dict[Path, Path | None] = {}
    attempted: list[Path] = []
    try:
        pairs = (
            (staged_values, output_values),
            (staged_points, output_points),
        )
        for staged, target in pairs:
            if staged.parent.resolve() != target.parent.resolve():
                raise ValueError("Data publication candidate is not beside its target")
        for target in (output_values, output_points, metadata_path):
            if target.is_symlink():
                raise ValueError("Refusing to replace a symbolic-link output target")
            if target.exists() and not target.is_file():
                raise ValueError("Pressure-database output target is not a file")

        metadata_temporary = _stage_json(metadata_path, metadata)
        replacements = (
            (staged_values, output_values),
            (staged_points, output_points),
            (metadata_temporary, metadata_path),
        )
        for _, target in replacements:
            backups[target] = (
                _stage_file_copy(target, target) if target.exists() else None
            )

        try:
            for staged, target in replacements:
                attempted.append(target)
                os.replace(staged, target)
        except BaseException as publish_error:
            rollback_errors: list[BaseException] = []
            for target in reversed(attempted):
                backup = backups[target]
                try:
                    if backup is None:
                        if target.exists() or target.is_symlink():
                            target.unlink()
                    else:
                        os.replace(backup, target)
                        backups[target] = None
                except BaseException as rollback_error:
                    rollback_errors.append(rollback_error)
            if rollback_errors:
                raise RuntimeError(
                    "Pressure-database publication failed and rollback was "
                    f"incomplete: {rollback_errors[0]}"
                ) from publish_error
            raise
    finally:
        _remove_temporary(staged_values)
        _remove_temporary(staged_points)
        _remove_temporary(metadata_temporary)
        for backup in backups.values():
            _remove_temporary(backup)


def particle_ids_sha256(particle_ids: np.ndarray) -> str:
    """Return a platform-independent digest of an ordered particle-ID list."""

    values = np.asarray(particle_ids, dtype="<u8").reshape(-1)
    return hashlib.sha256(values.tobytes()).hexdigest()


def database_paths(directory: Path, prefix: str) -> tuple[Path, Path]:
    return (
        directory / f"{prefix}_points.txt",
        directory / f"{prefix}_values.bin",
    )


def read_header(stream: BinaryIO) -> PressureHeader:
    raw = stream.read(HEADER.size)
    if len(raw) != HEADER.size:
        raise ValueError("Pressure database header is truncated")
    magic, dimension, values_per_particle, count, interval, max_step, dt = (
        HEADER.unpack(raw)
    )
    if magic not in (MAGIC_V1, MAGIC_V2):
        raise ValueError(f"Unexpected pressure database magic {magic!r}")
    if dimension not in (2, 3):
        raise ValueError(f"Unsupported pressure database dimension {dimension}")
    if values_per_particle != 2:
        raise ValueError("Pressure database must contain liquid and gas values")
    if count == 0 or interval == 0 or dt <= 0.0:
        raise ValueError("Pressure database header contains non-positive values")
    return PressureHeader(
        dimension=dimension,
        values_per_particle=values_per_particle,
        particle_count=count,
        step_interval=interval,
        max_step=max_step,
        source_dt=dt,
        format_version="V1" if magic == MAGIC_V1 else "V2",
    )


def write_header(stream: BinaryIO, header: PressureHeader) -> None:
    stream.write(
        HEADER.pack(
            MAGIC_V1 if header.format_version == "V1" else MAGIC_V2,
            header.dimension,
            header.values_per_particle,
            header.particle_count,
            header.step_interval,
            header.max_step,
            header.source_dt,
        )
    )


def read_points(path: Path, dimension: int) -> PressurePoints:
    rows: list[list[float]] = []
    with path.open(encoding="utf-8") as stream:
        for line in stream:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            values = [float(value) for value in stripped.replace(",", " ").split()]
            if len(values) < dimension + 1:
                raise ValueError(f"Invalid pressure point row: {line.rstrip()}")
            rows.append(values[: dimension + 1])
    if not rows:
        raise ValueError(f"No pressure points found in {path}")
    array = np.asarray(rows, dtype=np.float64)
    particle_ids = array[:, 0].astype(np.uint64)
    if len(np.unique(particle_ids)) != len(particle_ids):
        raise ValueError("Pressure point IDs are not unique")
    return PressurePoints(particle_ids=particle_ids, coordinates=array[:, 1:])


def validate_database(path: Path, header: PressureHeader) -> None:
    expected = HEADER.size + header.frame_count * header.frame_bytes
    actual = path.stat().st_size
    if actual != expected:
        raise ValueError(
            f"Pressure database size is {actual} bytes; expected {expected}. "
            "The source run may be incomplete."
        )


def iter_frames(
    stream: BinaryIO, header: PressureHeader
) -> Iterator[tuple[int, np.ndarray, np.ndarray]]:
    stream.seek(HEADER.size)
    value_bytes = header.particle_count * 8
    for frame_index in range(header.frame_count):
        raw_step = stream.read(FRAME_STEP.size)
        liquid_raw = stream.read(value_bytes)
        gas_raw = stream.read(value_bytes)
        if (
            len(raw_step) != FRAME_STEP.size
            or len(liquid_raw) != value_bytes
            or len(gas_raw) != value_bytes
        ):
            raise ValueError(f"Pressure frame {frame_index} is truncated")
        step = FRAME_STEP.unpack(raw_step)[0]
        expected_step = frame_index * header.step_interval
        if step != expected_step:
            raise ValueError(
                f"Frame {frame_index} stores step {step}; expected {expected_step}"
            )
        liquid = np.frombuffer(liquid_raw, dtype="<f8")
        gas = np.frombuffer(gas_raw, dtype="<f8")
        if not np.all(np.isfinite(liquid)) or not np.all(np.isfinite(gas)):
            raise ValueError(f"Pressure frame {frame_index} contains NaN or Inf")
        yield step, liquid, gas


def fit_harmonic(
    values_path: Path,
    header: PressureHeader,
    period: float,
    fit_start: float,
    fit_end: float,
) -> np.ndarray:
    """Return [mean, cosine coefficient, sine coefficient] for both phases."""
    if not math.isfinite(period) or period <= 0.0:
        raise ValueError("Wave period must be positive")
    if (
        not math.isfinite(fit_start)
        or not math.isfinite(fit_end)
        or fit_start < 0.0
        or fit_end <= fit_start
    ):
        raise ValueError("The harmonic fit window is invalid")
    omega = 2.0 * math.pi / period
    gram = np.zeros((3, 3), dtype=np.float64)
    rhs = np.zeros((3, header.particle_count, 2), dtype=np.float64)
    samples = 0
    with values_path.open("rb") as stream:
        read_header(stream)
        for step, liquid, gas in iter_frames(stream, header):
            time = header.physical_time(step)
            if time < fit_start - 1.0e-12 or time > fit_end + 1.0e-12:
                continue
            basis = np.asarray(
                [1.0, math.cos(omega * time), math.sin(omega * time)]
            )
            gram += np.outer(basis, basis)
            rhs[:, :, 0] += basis[:, None] * liquid[None, :]
            rhs[:, :, 1] += basis[:, None] * gas[None, :]
            samples += 1
    if samples < 6:
        raise ValueError("At least six pressure frames are required in fit window")
    if np.linalg.cond(gram) > 1.0e10:
        raise ValueError("Harmonic fit is ill-conditioned; enlarge the fit window")
    coefficients = np.linalg.solve(gram, rhs.reshape(3, -1)).reshape(rhs.shape)
    if not np.all(np.isfinite(coefficients)):
        raise ValueError("Harmonic fit produced NaN or Inf")
    return coefficients


def read_counted_particle_ids(path: Path) -> np.ndarray:
    """Read the standard MPM counted particle-set format."""

    with path.open(encoding="utf-8") as stream:
        rows = [
            line.strip()
            for line in stream
            if line.strip() and not line.lstrip().startswith(("#", "!"))
        ]
    if not rows:
        raise ValueError(f"Surface-reference particle set is empty: {path}")
    declared = int(rows[0])
    particle_ids = np.asarray([int(value) for value in rows[1:]], dtype=np.uint64)
    if particle_ids.size != declared:
        raise ValueError(
            f"{path} declares {declared} particle IDs but contains "
            f"{particle_ids.size}"
        )
    if np.unique(particle_ids).size != particle_ids.size:
        raise ValueError(f"Surface-reference particle IDs are not unique: {path}")
    return particle_ids


def physical_surface_sample_indices(
    points: PressurePoints,
    surface_band: float,
    reference_ids_path: Path | None = None,
) -> tuple[np.ndarray, dict[str, object]]:
    """Return the one physical top-row sample associated with every point.

    Production transforms use the immutable top-surface particle set generated
    with the mesh.  The coordinate fallback is intentionally limited to the
    highest point at each exact x-column and exists for small synthetic/test
    databases only; it never treats an entire two-row free-surface band as the
    phase reference.
    """

    if surface_band <= 0.0:
        raise ValueError("surface_band must be positive")
    coordinates = np.asarray(points.coordinates, dtype=np.float64)
    if coordinates.ndim != 2 or coordinates.shape[1] < 2:
        raise ValueError("Pressure points must contain at least x and y coordinates")

    if reference_ids_path is not None:
        reference_ids_path = reference_ids_path.resolve()
        requested_ids = read_counted_particle_ids(reference_ids_path)
        id_to_index = {
            int(particle_id): index
            for index, particle_id in enumerate(points.particle_ids)
        }
        missing = [
            int(particle_id)
            for particle_id in requested_ids
            if int(particle_id) not in id_to_index
        ]
        if missing:
            raise ValueError(
                "Surface-reference particle set is not contained in the pressure "
                f"database; first missing ID={missing[0]}"
            )
        candidates = np.asarray(
            [id_to_index[int(particle_id)] for particle_id in requested_ids],
            dtype=np.int64,
        )
        method = "registered_top_surface_particle_set"
        source: dict[str, object] = {
            "path": str(reference_ids_path),
            "sha256": sha256(reference_ids_path),
        }
    else:
        vertical = coordinates[:, 1]
        band = np.flatnonzero(vertical >= np.max(vertical) - surface_band)
        if band.size < 2:
            raise ValueError("Too few surface points; increase --surface-band")
        # Synthetic/test databases use exact x-columns.  Retaining only the
        # maximum y for each column prevents the old two-row ambiguity.
        by_x: dict[float, int] = {}
        for index in band:
            x = float(coordinates[index, 0])
            current = by_x.get(x)
            if current is None or coordinates[index, 1] > coordinates[current, 1]:
                by_x[x] = int(index)
        candidates = np.asarray(list(by_x.values()), dtype=np.int64)
        method = "highest_point_per_exact_x_column_fallback"
        source = {"path": None, "sha256": None}

    if candidates.size < 2:
        raise ValueError("Too few physical top-row reference points")
    order = np.argsort(coordinates[candidates, 0], kind="stable")
    candidates = candidates[order]
    candidate_x = coordinates[candidates, 0]
    if np.any(np.diff(candidate_x) <= 0.0):
        raise ValueError("Physical top-row reference points must have unique x values")
    target_x = coordinates[:, 0]
    right = np.searchsorted(candidate_x, target_x, side="left")
    right = np.clip(right, 0, len(candidates) - 1)
    left = np.clip(right - 1, 0, len(candidates) - 1)
    choose_right = np.abs(candidate_x[right] - target_x) < np.abs(
        candidate_x[left] - target_x
    )
    mapped = np.where(choose_right, candidates[right], candidates[left])
    reference_ids = points.particle_ids[candidates]
    audit: dict[str, object] = {
        "method": method,
        "source": source,
        "reference_particle_count": int(reference_ids.size),
        "reference_particle_ids": [int(value) for value in reference_ids],
        "reference_particle_ids_sha256": particle_ids_sha256(reference_ids),
        "reference_y_min_m": float(np.min(coordinates[candidates, 1])),
        "reference_y_max_m": float(np.max(coordinates[candidates, 1])),
    }
    return mapped, audit


def local_surface_indices(coordinates: np.ndarray, surface_band: float) -> np.ndarray:
    """Compatibility wrapper using the highest point in each exact x-column."""

    coordinates = np.asarray(coordinates, dtype=np.float64)
    points = PressurePoints(
        particle_ids=np.arange(len(coordinates), dtype=np.uint64),
        coordinates=coordinates,
    )
    indices, _ = physical_surface_sample_indices(points, surface_band)
    return indices


def rotated_coefficients(
    coefficients: np.ndarray,
    surface_indices: np.ndarray,
    minimum_reference_amplitude: float,
) -> tuple[np.ndarray, dict[str, object]]:
    if not np.all(np.isfinite(coefficients)):
        raise ValueError("Harmonic coefficients contain NaN or Inf")
    if (
        not math.isfinite(minimum_reference_amplitude)
        or minimum_reference_amplitude <= 0.0
    ):
        raise ValueError("Minimum reference amplitude must be finite and positive")
    rotated = coefficients.copy()
    original_amplitude = np.hypot(coefficients[1], coefficients[2])
    phase_alignment: dict[str, dict[str, float | int]] = {}
    for phase in range(2):
        surface_cos = coefficients[1, surface_indices, phase]
        surface_sin = coefficients[2, surface_indices, phase]
        surface_amplitude = np.hypot(surface_cos, surface_sin)
        valid = surface_amplitude > minimum_reference_amplitude
        rotated[1, valid, phase] = (
            original_amplitude[valid, phase]
            * surface_cos[valid]
            / surface_amplitude[valid]
        )
        rotated[2, valid, phase] = (
            original_amplitude[valid, phase]
            * surface_sin[valid]
            / surface_amplitude[valid]
        )
        source_phase = np.arctan2(coefficients[2, :, phase], coefficients[1, :, phase])
        surface_phase = np.arctan2(surface_sin, surface_cos)
        rotated_phase = np.arctan2(rotated[2, :, phase], rotated[1, :, phase])
        source_difference = (
            source_phase - surface_phase + math.pi
        ) % (2.0 * math.pi) - math.pi
        residual_difference = (
            rotated_phase - surface_phase + math.pi
        ) % (2.0 * math.pi) - math.pi
        phase_defined = valid & (original_amplitude[:, phase] > minimum_reference_amplitude)
        changed = phase_defined & (np.abs(source_difference) > 1.0e-10)
        label = "liquid" if phase == 0 else "gas"
        phase_alignment[label] = {
            "surface_reference_eligible_particle_count": int(np.count_nonzero(valid)),
            "phase_defined_particle_count": int(np.count_nonzero(phase_defined)),
            "changed_particle_count": int(np.count_nonzero(changed)),
            "max_abs_source_to_surface_phase_difference_deg": (
                float(np.max(np.abs(np.degrees(source_difference[phase_defined]))))
                if np.any(phase_defined)
                else 0.0
            ),
            "max_abs_residual_to_surface_phase_difference_deg": (
                float(np.max(np.abs(np.degrees(residual_difference[phase_defined]))))
                if np.any(phase_defined)
                else 0.0
            ),
        }
    new_amplitude = np.hypot(rotated[1], rotated[2])
    qa = {
        "max_abs_mean_change_pa": float(
            np.max(np.abs(rotated[0] - coefficients[0]))
        ),
        "max_abs_fundamental_amplitude_change_pa": float(
            np.max(np.abs(new_amplitude - original_amplitude))
        ),
        "phase_alignment": phase_alignment,
    }
    return rotated, qa


def _validated_phase_alignment_record(
    alignment: dict[str, object], phase_name: str
) -> dict[str, float | int]:
    """Validate one liquid/gas phase-alignment QA record."""

    label = phase_name.capitalize()
    value = alignment.get(phase_name)
    if not isinstance(value, dict):
        raise ValueError(f"Phase metadata is missing {phase_name} phase-alignment QA")

    count_names = (
        "surface_reference_eligible_particle_count",
        "phase_defined_particle_count",
        "changed_particle_count",
    )
    counts: dict[str, int] = {}
    for name in count_names:
        count = value.get(name)
        if isinstance(count, bool) or not isinstance(count, int):
            raise ValueError(f"{label} phase-alignment QA is missing integer {name}")
        counts[name] = count
    eligible = counts["surface_reference_eligible_particle_count"]
    defined = counts["phase_defined_particle_count"]
    changed = counts["changed_particle_count"]
    if eligible <= 0 or defined <= 0 or defined > eligible:
        raise ValueError(f"{label} phase-alignment QA has no valid defined phase")
    if changed <= 0 or changed > defined:
        raise ValueError(f"{label} phase-alignment QA records no changed particles")

    numeric_names = (
        "max_abs_source_to_surface_phase_difference_deg",
        "max_abs_residual_to_surface_phase_difference_deg",
    )
    numbers: dict[str, float] = {}
    for name in numeric_names:
        number_value = value.get(name)
        if isinstance(number_value, bool):
            raise ValueError(f"{label} phase-alignment QA is missing finite {name}")
        try:
            number = float(number_value)
        except (TypeError, ValueError) as error:
            raise ValueError(
                f"{label} phase-alignment QA is missing finite {name}"
            ) from error
        if not math.isfinite(number) or number < 0.0:
            raise ValueError(f"{label} phase-alignment QA has invalid {name}")
        numbers[name] = number
    if numbers["max_abs_source_to_surface_phase_difference_deg"] <= 0.0:
        raise ValueError(f"{label} source-to-surface phase difference is not positive")
    residual = numbers["max_abs_residual_to_surface_phase_difference_deg"]
    if residual > PHASE_ALIGNMENT_RESIDUAL_TOLERANCE_DEG:
        raise ValueError(
            f"{label} phase-alignment residual exceeds "
            f"{PHASE_ALIGNMENT_RESIDUAL_TOLERANCE_DEG:.1e} deg: {residual:.6g}"
        )
    return {**counts, **numbers}


def validate_phase_alignment_qa(
    metadata: dict[str, object],
) -> dict[str, float | int]:
    """Fail closed unless both stored pressures were nontrivially phase-aligned.

    The residual is recorded in degrees.  Floating-point roundoff at or below
    ``PHASE_ALIGNMENT_RESIDUAL_TOLERANCE_DEG`` is accepted; a larger residual
    means the phase-erasure transform did not align the retained fundamental
    with its registered surface reference.
    """

    qa = metadata.get("qa")
    if not isinstance(qa, dict):
        raise ValueError("Phase metadata is missing qa")
    alignment = qa.get("phase_alignment")
    if not isinstance(alignment, dict):
        raise ValueError("Phase metadata is missing qa.phase_alignment")
    liquid = _validated_phase_alignment_record(alignment, "liquid")
    _validated_phase_alignment_record(alignment, "gas")
    return liquid


def ramp_factor(time: float, ramp_time: float) -> float:
    if ramp_time <= 0.0:
        return 1.0
    return min(1.0, max(0.0, time / ramp_time))


def validate_transform_parameters(
    *,
    period: float,
    fit_start: float,
    fit_end: float,
    ramp_time: float,
    surface_band: float,
    minimum_reference_amplitude: float,
) -> None:
    """Reject non-finite or physically invalid transform controls up front."""

    values = {
        "period": period,
        "fit_start": fit_start,
        "fit_end": fit_end,
        "ramp_time": ramp_time,
        "surface_band": surface_band,
        "minimum_reference_amplitude": minimum_reference_amplitude,
    }
    if any(
        isinstance(value, bool) or not math.isfinite(value)
        for value in values.values()
    ):
        raise ValueError("Phase-transform parameters must be finite numbers")
    if period <= 0.0:
        raise ValueError("Wave period must be positive")
    if fit_start < 0.0 or fit_end <= fit_start:
        raise ValueError("The harmonic fit window is invalid")
    if ramp_time < 0.0:
        raise ValueError("Ramp time must be non-negative")
    if surface_band <= 0.0:
        raise ValueError("Surface band must be positive")
    if minimum_reference_amplitude <= 0.0:
        raise ValueError("Minimum reference amplitude must be positive")


def transform_database(
    source_directory: Path,
    source_prefix: str,
    output_directory: Path,
    output_prefix: str,
    *,
    period: float,
    fit_start: float,
    fit_end: float,
    ramp_time: float,
    surface_band: float,
    minimum_reference_amplitude: float,
    surface_reference_ids: Path | None = None,
) -> Path:
    validate_transform_parameters(
        period=period,
        fit_start=fit_start,
        fit_end=fit_end,
        ramp_time=ramp_time,
        surface_band=surface_band,
        minimum_reference_amplitude=minimum_reference_amplitude,
    )
    source_points, source_values = database_paths(source_directory, source_prefix)
    output_points, output_values = database_paths(output_directory, output_prefix)
    metadata_path = output_directory / f"{output_prefix}_metadata.json"
    _validate_distinct_output_targets(
        (source_points, source_values),
        (output_points, output_values, metadata_path),
    )
    output_directory.mkdir(parents=True, exist_ok=True)

    with source_values.open("rb") as stream:
        header = read_header(stream)
    validate_database(source_values, header)
    points = read_points(source_points, header.dimension)
    if len(points.particle_ids) != header.particle_count:
        raise ValueError("Point and binary particle counts differ")

    coefficients = fit_harmonic(
        source_values, header, period, fit_start, fit_end
    )
    surface_indices, surface_reference = physical_surface_sample_indices(
        points, surface_band, surface_reference_ids
    )
    rotated, qa = rotated_coefficients(
        coefficients, surface_indices, minimum_reference_amplitude
    )
    validate_phase_alignment_qa({"qa": qa})

    omega = 2.0 * math.pi / period
    max_pressure_delta = 0.0
    staged_values: Path | None = None
    staged_points: Path | None = None
    try:
        descriptor, staged_values = _temporary_file(output_values)
        try:
            os.fchmod(descriptor, source_values.stat().st_mode & 0o7777)
            with source_values.open("rb") as source, os.fdopen(
                descriptor, "wb"
            ) as output:
                read_header(source)
                write_header(output, header)
                for step, liquid, gas in iter_frames(source, header):
                    time = header.physical_time(step)
                    harmonic_basis = np.asarray(
                        [math.cos(omega * time), math.sin(omega * time)]
                    )
                    original_fundamental = np.einsum(
                        "cnp,c->np", coefficients[1:], harmonic_basis
                    )
                    rotated_fundamental = np.einsum(
                        "cnp,c->np", rotated[1:], harmonic_basis
                    )
                    correction = ramp_factor(time, ramp_time) * (
                        rotated_fundamental - original_fundamental
                    )
                    transformed = np.column_stack((liquid, gas)) + correction
                    if not np.all(np.isfinite(correction)) or not np.all(
                        np.isfinite(transformed)
                    ):
                        raise ValueError("Phase-erased transform produced NaN or Inf")
                    max_pressure_delta = max(
                        max_pressure_delta, float(np.max(np.abs(correction)))
                    )
                    output.write(FRAME_STEP.pack(step))
                    output.write(
                        np.asarray(transformed[:, 0], dtype="<f8").tobytes()
                    )
                    output.write(
                        np.asarray(transformed[:, 1], dtype="<f8").tobytes()
                    )
                output.flush()
                os.fsync(output.fileno())
        except BaseException:
            try:
                os.close(descriptor)
            except OSError:
                pass
            raise

        if not math.isfinite(max_pressure_delta) or max_pressure_delta <= 0.0:
            raise ValueError("Phase-erased transform made no finite pressure change")

        staged_points = _stage_file_copy(source_points, output_points)
        validate_database(staged_values, header)

        metadata = {
            "schema": "phase-erased-pressure-control-v1",
            "classification": "one-way numerical counterfactual; not fully coupled",
            "producer": {
                "path": str(Path(__file__).resolve()),
                "sha256": sha256(Path(__file__).resolve()),
            },
            "source": {
                "points": str(source_points),
                "values": str(source_values),
                "points_sha256": sha256(source_points),
                "values_sha256": sha256(source_values),
            },
            "output": {
                "points": str(output_points),
                "values": str(output_values),
                "points_sha256": sha256(staged_points),
                "values_sha256": sha256(staged_values),
            },
            "header": {
                "format_version": header.format_version,
                "dimension": header.dimension,
                "particle_count": header.particle_count,
                "step_interval": header.step_interval,
                "max_step": header.max_step,
                "source_dt_s": header.source_dt,
                "sample_dt_s": header.sample_dt,
            },
            "transform": {
                "period_s": period,
                "fit_start_s": fit_start,
                "fit_end_s": fit_end,
                "ramp_time_s": ramp_time,
                "surface_band_m": surface_band,
                "minimum_reference_amplitude_pa": minimum_reference_amplitude,
                "surface_reference": surface_reference,
                "description": (
                    "rotate each point's fundamental liquid/gas pressure component "
                    "to its local surface phase; retain mean, amplitude, residuals, "
                    "higher harmonics, and along-wave progressive phase"
                ),
            },
            "qa": {
                **qa,
                "max_instantaneous_pressure_change_pa": max_pressure_delta,
            },
        }
        _publish_transform_files(
            staged_points=staged_points,
            staged_values=staged_values,
            output_points=output_points,
            output_values=output_values,
            metadata_path=metadata_path,
            metadata=metadata,
        )
    finally:
        _remove_temporary(staged_values)
        _remove_temporary(staged_points)
    return metadata_path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-dir", type=Path, required=True)
    parser.add_argument("--source-prefix", default="pressure")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--output-prefix", default="phase_erased")
    parser.add_argument("--period", type=float, default=1.3)
    parser.add_argument("--fit-start", type=float, default=2.6)
    parser.add_argument("--fit-end", type=float, default=7.8)
    parser.add_argument("--ramp-time", type=float, default=1.3)
    parser.add_argument("--surface-band", type=float, default=0.015)
    parser.add_argument(
        "--surface-reference-ids",
        type=Path,
        default=Path(__file__).resolve().parent
        / "top_surface_traction_particle_id.txt",
        help=(
            "counted immutable particle set defining the one physical seabed "
            "top row"
        ),
    )
    parser.add_argument("--minimum-reference-amplitude", type=float, default=1.0)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    metadata = transform_database(
        args.source_dir,
        args.source_prefix,
        args.output_dir,
        args.output_prefix,
        period=args.period,
        fit_start=args.fit_start,
        fit_end=args.fit_end,
        ramp_time=args.ramp_time,
        surface_band=args.surface_band,
        minimum_reference_amplitude=args.minimum_reference_amplitude,
        surface_reference_ids=args.surface_reference_ids,
    )
    print(f"Wrote phase-erased database and audit metadata: {metadata}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
