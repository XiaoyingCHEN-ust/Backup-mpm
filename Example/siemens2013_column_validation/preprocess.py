"""Generate the 2-D MPM mesh and particle sets for Siemens et al. (2013).

The model represents the initially unsaturated interval between the upper
pond and the lower wet boundary in the coarse transparent-sand column. Since
the validation targets are one-dimensional, one square cell is used across a
representative-width slice. Two particles per horizontal layer retain the
standard 2 x 2 particles per cell, while 188 particle layers resolve the
approximately 1.075 m experimental interval.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np


CASE_DIR = Path(__file__).resolve().parent
CELL_SIZE = 0.0114
NX = 1
NY = 94
WIDTH = NX * CELL_SIZE
HEIGHT = NY * CELL_SIZE
PARTICLE_SPACING = CELL_SIZE / 2.0


def gimp_mesh2d(width: float, height: float, cell_size: float):
    """Create the 16-node GIMP connectivity used by ED2Q16G."""
    x_range = np.arange(-cell_size, width + cell_size + 1.0e-14, cell_size)
    y_range = np.arange(-cell_size, height + cell_size + 1.0e-14, cell_size)
    nx = x_range.size
    ny = y_range.size

    x, y = np.meshgrid(x_range, y_range, indexing="ij")
    nodes = np.array([x, y]).T.reshape(-1, 2)

    v0 = np.add.outer(np.arange(1, nx - 2), nx * np.arange(1, ny - 2)).T.ravel()
    v1 = np.add.outer(np.arange(0, nx - 3), nx * np.arange(0, ny - 3)).T.ravel()
    group0 = np.repeat([v0], 4, axis=0)
    group1 = np.repeat([v1], 12, axis=0)
    elements = np.concatenate((group0, group1), axis=0).T

    elements[:, 0] += 0
    elements[:, 1] += 1
    elements[:, 2] += nx + 1
    elements[:, 3] += nx

    elements[:, 4] += 0
    elements[:, 5] += 1
    elements[:, 6] += 2
    elements[:, 7] += 3
    elements[:, 8] += nx + 3
    elements[:, 9] += 2 * nx + 3
    elements[:, 10] += 3 * nx + 3
    elements[:, 11] += 3 * nx + 2
    elements[:, 12] += 3 * nx + 1
    elements[:, 13] += 3 * nx
    elements[:, 14] += 2 * nx
    elements[:, 15] += nx
    return nodes, elements


def particles2d(width: float, height: float, spacing: float):
    x_range = np.arange(0.5 * spacing, width, spacing)
    y_range = np.arange(0.5 * spacing, height, spacing)
    x, y = np.meshgrid(x_range, y_range, indexing="ij")
    return np.array([x, y]).T.reshape(-1, 2)


def write_mesh(path: Path, nodes: np.ndarray, elements: np.ndarray):
    with path.open("w", encoding="utf-8", newline="\n") as stream:
        stream.write("#! elementShape quadrilateral\n")
        stream.write("#! elementNumPoints 16\n")
        stream.write(f"{nodes.shape[0]} {elements.shape[0]}\n")
        np.savetxt(stream, nodes, fmt="%.8f", delimiter="\t")
        np.savetxt(stream, elements, fmt="%i", delimiter="\t")


def write_particles(path: Path, particles: np.ndarray):
    with path.open("w", encoding="utf-8", newline="\n") as stream:
        stream.write(f"{particles.shape[0]}\n")
        np.savetxt(stream, particles, fmt="%.8f", delimiter="\t")


def write_particle_scalar(path: Path, values: np.ndarray):
    data = np.zeros((values.size, 2), dtype=float)
    data[:, 0] = values
    with path.open("w", encoding="utf-8", newline="\n") as stream:
        stream.write(f"{values.size}\n")
        np.savetxt(stream, data, fmt="%.8e", delimiter="\t")


def as_int_list(values):
    return np.asarray(values, dtype=int).ravel().tolist()


def main():
    nodes, elements = gimp_mesh2d(WIDTH, HEIGHT, CELL_SIZE)
    particles = particles2d(WIDTH, HEIGHT, PARTICLE_SPACING)

    tolerance = 1.0e-10
    node_sets = [
        {"id": 0, "set": as_int_list(np.arange(nodes.shape[0]))},
        {"id": 1, "set": as_int_list(np.where(nodes[:, 1] <= tolerance))},
        {"id": 2, "set": as_int_list(np.where(nodes[:, 0] <= tolerance))},
        {"id": 3, "set": as_int_list(np.where(nodes[:, 0] >= WIDTH - tolerance))},
    ]

    # Apply the imposed surface pressures to all four material points in the
    # top cell, rather than to only the uppermost particle row.
    top = np.where(particles[:, 1] >= HEIGHT - CELL_SIZE + tolerance)[0]
    bottom = np.where(particles[:, 1] <= 0.51 * PARTICLE_SPACING)[0]
    nonfree = np.setdiff1d(np.arange(particles.shape[0]), top)
    particle_sets = [
        {"id": 0, "set": as_int_list(nonfree)},
        {"id": 1, "set": as_int_list(top)},
        {"id": 2, "set": as_int_list(bottom)},
    ]

    entity_sets = {"node_sets": node_sets, "particle_sets": particle_sets}

    write_mesh(CASE_DIR / "gimp_mesh2d.txt", nodes, elements)
    write_particles(CASE_DIR / "particles.txt", particles)
    write_particle_scalar(
        CASE_DIR / "initial_temperature.txt",
        np.full(particles.shape[0], 25.0),
    )
    write_particle_scalar(
        CASE_DIR / "initial_volume.txt",
        np.full(particles.shape[0], PARTICLE_SPACING**2),
    )
    with (CASE_DIR / "entity_sets.json").open(
        "w", encoding="utf-8", newline="\n"
    ) as stream:
        json.dump(entity_sets, stream, indent=2)
        stream.write("\n")

    summary = {
        "width_m": WIDTH,
        "height_m": HEIGHT,
        "cell_size_m": CELL_SIZE,
        "cells": int(elements.shape[0]),
        "particles": int(particles.shape[0]),
        "top_boundary_particles": int(top.size),
        "bottom_boundary_particles": int(bottom.size),
    }
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
