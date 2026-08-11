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
import shutil
import struct
from dataclasses import dataclass
from pathlib import Path
from typing import BinaryIO, Iterator

import numpy as np


MAGIC_V1 = b"MPM_PRESSURE_V1\0"
MAGIC_V2 = b"MPM_PRESSURE_V2\0"
HEADER = struct.Struct("<16sIIQQQd")
FRAME_STEP = struct.Struct("<Q")


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
    if period <= 0.0:
        raise ValueError("Wave period must be positive")
    if fit_start < 0.0 or fit_end <= fit_start:
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
    return coefficients


def local_surface_indices(coordinates: np.ndarray, surface_band: float) -> np.ndarray:
    if surface_band <= 0.0:
        raise ValueError("surface_band must be positive")
    vertical = coordinates[:, 1]
    candidates = np.flatnonzero(vertical >= np.max(vertical) - surface_band)
    if candidates.size < 2:
        raise ValueError("Too few surface points; increase --surface-band")
    order = np.argsort(coordinates[candidates, 0])
    candidates = candidates[order]
    candidate_x = coordinates[candidates, 0]
    target_x = coordinates[:, 0]
    right = np.searchsorted(candidate_x, target_x, side="left")
    right = np.clip(right, 0, len(candidates) - 1)
    left = np.clip(right - 1, 0, len(candidates) - 1)
    choose_right = np.abs(candidate_x[right] - target_x) < np.abs(
        candidate_x[left] - target_x
    )
    return np.where(choose_right, candidates[right], candidates[left])


def rotated_coefficients(
    coefficients: np.ndarray,
    surface_indices: np.ndarray,
    minimum_reference_amplitude: float,
) -> tuple[np.ndarray, dict[str, float]]:
    rotated = coefficients.copy()
    original_amplitude = np.hypot(coefficients[1], coefficients[2])
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
    new_amplitude = np.hypot(rotated[1], rotated[2])
    qa = {
        "max_abs_mean_change_pa": float(
            np.max(np.abs(rotated[0] - coefficients[0]))
        ),
        "max_abs_fundamental_amplitude_change_pa": float(
            np.max(np.abs(new_amplitude - original_amplitude))
        ),
    }
    return rotated, qa


def ramp_factor(time: float, ramp_time: float) -> float:
    if ramp_time <= 0.0:
        return 1.0
    return min(1.0, max(0.0, time / ramp_time))


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
) -> Path:
    source_points, source_values = database_paths(source_directory, source_prefix)
    output_points, output_values = database_paths(output_directory, output_prefix)
    if source_values.resolve() == output_values.resolve():
        raise ValueError("Refusing to overwrite the source pressure database")
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
    surface_indices = local_surface_indices(points.coordinates, surface_band)
    rotated, qa = rotated_coefficients(
        coefficients, surface_indices, minimum_reference_amplitude
    )

    omega = 2.0 * math.pi / period
    max_pressure_delta = 0.0
    with source_values.open("rb") as source, output_values.open("wb") as output:
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
            max_pressure_delta = max(
                max_pressure_delta, float(np.max(np.abs(correction)))
            )
            output.write(FRAME_STEP.pack(step))
            output.write(np.asarray(transformed[:, 0], dtype="<f8").tobytes())
            output.write(np.asarray(transformed[:, 1], dtype="<f8").tobytes())

    shutil.copyfile(source_points, output_points)
    validate_database(output_values, header)

    metadata = {
        "schema": "phase-erased-pressure-control-v1",
        "classification": "one-way numerical counterfactual; not fully coupled",
        "source": {
            "points": str(source_points),
            "values": str(source_values),
            "points_sha256": sha256(source_points),
            "values_sha256": sha256(source_values),
        },
        "output": {
            "points": str(output_points),
            "values": str(output_values),
            "points_sha256": sha256(output_points),
            "values_sha256": sha256(output_values),
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
            "description": (
                "rotate each point's fundamental liquid/gas pressure component "
                "to its local surface phase; retain mean, amplitude, residuals, "
                "higher harmonics, and along-wave progressive phase"
            ),
        },
        "qa": {**qa, "max_instantaneous_pressure_change_pa": max_pressure_delta},
    }
    metadata_path = output_directory / f"{output_prefix}_metadata.json"
    metadata_path.write_text(
        json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
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
    )
    print(f"Wrote phase-erased database and audit metadata: {metadata}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
