"""Generate low-cost explicit/semi-implicit three-phase comparison runs."""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
from pathlib import Path


CASE_DIR = Path(__file__).resolve().parent

MESH_VARIANTS = {
    "coarse": ("gimp_mesh2d.txt", "entity_sets.json", 0.0114),
    "one_particle_per_cell": (
        "gimp_mesh2d_one_particle_per_cell.txt",
        "entity_sets_one_particle_per_cell.json",
        0.0057,
    ),
    "one_layer_per_cell": (
        "gimp_mesh2d_one_layer_per_cell.txt",
        "entity_sets_one_layer_per_cell.json",
        0.0057,
    ),
}


def step_tag(value: float) -> str:
    mantissa, exponent = f"{value:.12e}".split("e")
    mantissa = mantissa.rstrip("0").rstrip(".").replace(".", "p")
    return f"{mantissa}e{int(exponent):+03d}".replace("+", "")


def checked_steps(duration: float, dt: float) -> int:
    nsteps = round(duration / dt)
    tolerance = max(1.0e-12, duration * 1.0e-10)
    if nsteps < 1 or abs(nsteps * dt - duration) > tolerance:
        raise ValueError(
            f"duration {duration:g} s is not an integer multiple of dt {dt:g} s"
        )
    return nsteps


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--duration", type=float, default=0.05)
    parser.add_argument("--explicit-dt", type=float, default=2.5e-6)
    parser.add_argument(
        "--semi-implicit-dts",
        type=float,
        nargs="+",
        default=(2.5e-5, 1.0e-4, 5.0e-4),
    )
    parser.add_argument("--revision", default="r29-lumped-pressure-smoke")
    parser.add_argument("--boundary-penalty", type=float, default=1.0e6)
    parser.add_argument(
        "--mesh-variant",
        choices=tuple(MESH_VARIANTS),
        default="coarse",
        help="background mesh; particle positions and volumes are unchanged",
    )
    parser.add_argument(
        "--bounded-pressure-transfer",
        action="store_true",
        help=(
            "enable the diagnostic r39 nodal pressure limiter; it is disabled "
            "by default because it over-predicted the early wetting-front depth"
        ),
    )
    parser.add_argument(
        "--reconstruct-pressure-gradient",
        action="store_true",
        help=(
            "recompute particle pressure gradients from the current nodal "
            "solution without projecting particle pressure values"
        ),
    )
    parser.add_argument(
        "--reconstruct-pressure-force",
        action="store_true",
        help=(
            "use the current nodal pressure interpolation only for phase "
            "momentum internal forces; particle pressure values remain FLIP"
        ),
    )
    parser.add_argument(
        "--reconstruct-darcy-velocity",
        action="store_true",
        help=(
            "reconstruct phase velocities from the Darcy relation used by "
            "the semi-implicit pressure equation instead of fluid FLIP"
        ),
    )
    parser.add_argument(
        "--gravity-scale",
        type=float,
        default=1.0,
        help=(
            "multiply the source-input gravity vector by this value; use 0 "
            "only for a diagnostic while retaining the prescribed surface pressure"
        ),
    )
    parser.add_argument(
        "--boundary-penalties",
        type=float,
        nargs="+",
        default=None,
        help="optional semi-implicit boundary-penalty sweep",
    )
    args = parser.parse_args()

    if args.duration <= 0.0:
        parser.error("--duration must be positive")
    if args.explicit_dt <= 0.0 or any(dt <= 0.0 for dt in args.semi_implicit_dts):
        parser.error("all time steps must be positive")
    boundary_penalties = (
        args.boundary_penalties
        if args.boundary_penalties is not None
        else [args.boundary_penalty]
    )
    if any(value <= 0.0 for value in boundary_penalties):
        parser.error("all boundary penalties must be positive")
    if not math.isfinite(args.gravity_scale):
        parser.error("--gravity-scale must be finite")
    if not re.fullmatch(r"[A-Za-z0-9_-]+", args.revision):
        parser.error("--revision may contain only letters, numbers, '_' and '-'")
    reconstruct_pressure_gradient = (
        args.reconstruct_pressure_gradient
        or args.reconstruct_darcy_velocity
        or args.bounded_pressure_transfer
    )

    output_dir = CASE_DIR / "semi_implicit_inputs"
    output_dir.mkdir(exist_ok=True)
    manifest = output_dir / "manifest.csv"
    manifest.unlink(missing_ok=True)

    rows: list[dict[str, object]] = []
    methods = [("explicit", args.explicit_dt, boundary_penalties[0])] + [
        ("semi_implicit", dt, penalty)
        for dt in args.semi_implicit_dts
        for penalty in boundary_penalties
    ]
    for case_name in ("open", "closed"):
        source_path = CASE_DIR / f"mpm_{case_name}.json"
        if not source_path.is_file():
            parser.error(f"missing source input: {source_path}")
        source = json.loads(source_path.read_text(encoding="utf-8"))

        for method, dt, boundary_penalty in methods:
            try:
                nsteps = checked_steps(args.duration, dt)
            except ValueError as error:
                parser.error(str(error))
            label = f"{case_name}_{method}_{step_tag(dt)}"
            if method == "semi_implicit" and len(boundary_penalties) > 1:
                label += f"_penalty_{step_tag(boundary_penalty)}"
            uuid = f"siemens2013-{args.revision}-{label}"
            config = json.loads(json.dumps(source))
            mesh_file, entity_sets_file, cell_size = MESH_VARIANTS[
                args.mesh_variant
            ]
            for asset in (mesh_file, entity_sets_file):
                if not (CASE_DIR / asset).is_file():
                    parser.error(
                        f"missing mesh asset {asset}; run python preprocess.py"
                    )
            config["mesh"].update(
                {
                    "mesh": mesh_file,
                    "entity_sets": entity_sets_file,
                    "cellsize_min": cell_size,
                }
            )
            source_gravity = config["external_loading_conditions"]["gravity"]
            config["external_loading_conditions"]["gravity"] = [
                0.0
                if args.gravity_scale == 0.0
                else args.gravity_scale * float(component)
                for component in source_gravity
            ]
            config["title"] = (
                f"Siemens 2013 {case_name} {method} pressure smoke test, "
                f"dt={dt:g} s, duration={args.duration:g} s, "
                f"gravity scale={args.gravity_scale:g}, "
                f"mesh={args.mesh_variant}"
            )
            analysis = config["analysis"]
            analysis.update(
                {
                    "PIC": 0.0,
                    "PIC_T": 0.0,
                    "pressure_smoothing": False,
                    "pressure_integration": method,
                    "dt": dt,
                    "nsteps": nsteps,
                    "uuid": uuid,
                }
            )
            analysis["semi_implicit_pressure"] = {
                "boundary_penalty": boundary_penalty,
                "log_solver": True,
                "bounded_transfer": args.bounded_pressure_transfer,
                "reconstruct_gradient": reconstruct_pressure_gradient,
                "reconstruct_force": args.reconstruct_pressure_force,
                "reconstruct_darcy_velocity": (
                    args.reconstruct_darcy_velocity
                ),
            }
            analysis["resume"].update(
                {
                    "resume": False,
                    "uuid": uuid,
                    "step": 0,
                    "nsteps": 0,
                    "start_from_this_step": False,
                    "this_step": 0,
                    "current_time": 0.0,
                }
            )
            output_steps = max(1, nsteps // 10)
            config["post_processing"].update(
                {
                    "path": "stability_results/",
                    "write_hdf5": False,
                    "write_vtk": True,
                    "log_output_steps": False,
                    "output_steps": output_steps,
                }
            )
            if args.reconstruct_pressure_force:
                liquid_vtk = config["post_processing"].setdefault(
                    "liquid_vtk", []
                )
                for field in (
                    "force_liquid_pressures",
                    "force_gas_pressures",
                ):
                    if field not in liquid_vtk:
                        liquid_vtk.append(field)
            if reconstruct_pressure_gradient:
                liquid_vtk = config["post_processing"].setdefault(
                    "liquid_vtk", []
                )
                for field in (
                    "liquid_pressure_gradients",
                    "gas_pressure_gradients",
                ):
                    if field not in liquid_vtk:
                        liquid_vtk.append(field)

            input_path = output_dir / f"mpm_{args.revision}_{label}.json"
            input_path.write_text(
                json.dumps(config, indent=2) + "\n", encoding="utf-8"
            )
            rows.append(
                {
                    "array_index": len(rows),
                    "label": label,
                    "case": case_name,
                    "method": method,
                    "dt_s": dt,
                    "duration_s": args.duration,
                    "nsteps": nsteps,
                    "output_steps": output_steps,
                    "uuid": uuid,
                    "input": input_path.relative_to(CASE_DIR).as_posix(),
                    "boundary_penalty": boundary_penalty,
                    "gravity_scale": args.gravity_scale,
                    "mesh_variant": args.mesh_variant,
                    "bounded_pressure_transfer": (
                        args.bounded_pressure_transfer
                    ),
                    "reconstruct_pressure_gradient": (
                        reconstruct_pressure_gradient
                    ),
                    "reconstruct_pressure_force": (
                        args.reconstruct_pressure_force
                    ),
                    "reconstruct_darcy_velocity": (
                        args.reconstruct_darcy_velocity
                    ),
                }
            )

    with manifest.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {len(rows)} inputs and {manifest}")
    open_last = len(methods) - 1
    closed_first = len(methods)
    closed_last = 2 * len(methods) - 1
    print(
        "First submit only the open checks: "
        f"sbatch --array=0-{open_last} run_semi_implicit_hpc4.sh"
    )
    print(
        "After they pass, submit closed checks: "
        f"sbatch --array={closed_first}-{closed_last} "
        "run_semi_implicit_hpc4.sh"
    )


if __name__ == "__main__":
    main()
