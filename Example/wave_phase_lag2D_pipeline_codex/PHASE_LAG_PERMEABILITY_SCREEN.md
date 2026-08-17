# Unsmoothed phase-lag permeability screen

## Purpose and scope

This is a fixed four-point exploratory screen for the pipeline engineering
example. It asks whether an auditable subsurface liquid-pressure phase lag is
converted into a larger upward pressure-gradient driver near the fixed pipe.
It does not by itself establish additional liquefaction or pipe uplift. Those
claims require independently integrated lagged/phase-erased SANISAND cases and,
for the engineering response, a released-pipeline comparison.

The four intrinsic-permeability points and effect gates were fixed before
completion of the matrix. The points are
`1e-13, 3e-13, 1e-12, 3e-12 m2`. All use `Sw=0.94`. The `3e-13` calculation is
run regardless of the `3e-12` result so that the matrix is not selected after
looking at the outcome.

## Invariants

- Equilibrium: `LinearElastic2D`, pure PIC (`PIC=1`, `APIC=false`), fixed pipe,
  no wave, `dt=1e-4 s`, Cundall damping `5 s-1`.
- Exploratory equilibrium length: 5000 steps. This is a screening checkpoint,
  not the 40000-step registered production equilibrium.
- The original equilibrium gate remains `max|v|<=1e-3 m/s`; displacement,
  porosity, normal-stress, particle-count, grid-bound and HDF5/VTP checks also
  remain mandatory.
- Dynamic driver: SANISAND, APIC, fixed pipe, zero damping, 39000 steps, 60
  material-point frames, 301 pressure-database frames.
- Both solver pressure-smoothing switches are exactly JSON `false`.
- The phase-erased control rotates only the fitted fundamental to the local
  registered top-surface phase. It retains each point's mean, amplitude,
  residual/higher harmonics and progressive along-wave phase.
- Any raster smoothing is display-only (`sigma=1.25 pixels` by default). Raw
  particle pressures are differentiated and integrated before plotting.

## Locked qualification and effect rules

Pressure phase is eligible only when the harmonic fit has `R2>=0.8` and local
amplitude is at least 5% of the same-column surface amplitude. The crown must
have at least 18 degrees of resolvable same-column phase difference.

The fixed-state pressure-only screen uses the two complete post-ramp cycles
`1.3--2.6 s` and `2.6--3.9 s`. It passes only when all of the following hold:

1. lagged-minus-erased signed and net-uplift impulses are positive in both
   cycles;
2. the combined net-uplift impulse increases by at least 5%;
3. combined `IF>=1` area-time has the positive direction;
4. the shared-HD co-location of `IF>=1` and `Rsigma<=0.05` is nonzero and
   increases by at least 5%; and
5. a lagged co-location advantage of at least two reference particle areas is
   resolved for two consecutive saved frames.

Local positive-part activity and a selected instantaneous 2-D maximum are
reported as mechanism diagnostics but cannot rescue a failed net-uplift gate.
The shared-HD co-location is not causal liquefaction evidence.

## Fail-closed commands

For a fresh label, run sequentially from this directory:

```bash
python3 run_phase_lag_exploration.py prepare \
  --label LABEL --permeability K_M2 --saturation 0.94 \
  --equilibrium-steps 5000
python3 run_phase_lag_exploration.py eq --label LABEL --threads 3
python3 run_phase_lag_exploration.py hd --label LABEL --threads 3
```

The runner refuses a reused/partial namespace and records config, binary,
result and pressure-database hashes. After `hd_complete`, generate a new,
previously absent phase-erased directory, then run the transform and analyses:

```bash
python3 phase_controls.py \
  --source-dir pressure_databases/phase_lag_exploratory/LABEL/lagged \
  --source-prefix pressure \
  --output-dir pressure_databases/phase_lag_exploratory/LABEL/phase_erased \
  --output-prefix phase_erased \
  --period 1.3 --fit-start 2.6 --fit-end 3.9 --ramp-time 1.3 \
  --surface-reference-ids top_surface_traction_particle_id.txt
python3 plot_liquid_pressure_phase_lag.py \
  --output-dir analysis/phase_lag_exploratory/LABEL/phase_v3 \
  --config configs/phase_lag_exploratory/LABEL/03_HD.json \
  --result-dir results/phase_lag_exploratory/LABEL/PLP_EXP_UUID_TOKEN_HD \
  --checkpoint results/phase_lag_exploratory/LABEL/PLP_EXP_UUID_TOKEN_EQ/particle5000.vtp \
  --particles particles.txt --case-label CASE_LABEL --require-unsmoothed
python3 analyze_phase_lag_driver_screen.py \
  --label LABEL --display-smoothing-sigma 1.25
```

`CASE_LABEL` is an identity field, not a free-form caption: it must equal
`LABEL.upper()`. `UUID_TOKEN` must equal `CASE_LABEL.replace("-", "_")`, and
the runner UUIDs must be `PLP_EXP_UUID_TOKEN_EQ` and
`PLP_EXP_UUID_TOKEN_HD`. For example, `LABEL=k3e-13_sw094_nosmooth` requires
`CASE_LABEL=K3E-13_SW094_NOSMOOTH` and
`UUID_TOKEN=K3E_13_SW094_NOSMOOTH`. Do not point the phase plot at `02_WAVE` or
an `_WAVE` result when producing the current HD/pressure-control evidence.

## Current audited status (2026-08-17)

| Intrinsic permeability | Status | Crown lag / amplitude ratio | Pressure-only engineering screen |
|---:|---|---:|---|
| `1e-13 m2` | EQ/HD/phase/2-D and raw-gradient audit complete | `+49.90 deg / 0.132` | failed; net uplift `-0.0239%`, joint zero |
| `3e-13 m2` | EQ/HD/phase/2-D and raw-gradient audit complete | `+45.92 deg / 0.157` | failed; net uplift `-0.0131%`, joint advantage zero |
| `1e-12 m2` | EQ/HD/phase/2-D and raw-gradient audit complete | `+24.00 deg / 0.193` | failed; net uplift `-0.1399%`, joint advantage zero |
| `3e-12 m2` | EQ/HD/phase/2-D and raw-gradient audit complete | `-5.55 deg / 0.163` | failed; crown phase unresolved, net uplift `-0.9038%`, IF area-time `-0.8552%` |

An earlier `1e-13` run remains useful only as a reproducibility reference. Its
original v2 audit reported `+50.20 deg` from the `02_WAVE` / `_WAVE` solver
output, whereas the later phase-v3 audit reported `+50.02 deg` from the
distinct `03_HD` / `_HD` pressure-recording output used by the current driver
analysis. These values are from different solver cases, not alternative
rounding of one result; `+50.02 deg` is the applicable historical HD value.
That run failed the same net-uplift interpretation, and its historical EQ
configuration identity was not frozen by the current runner schema. It is
therefore excluded from the primary four-point summary.
No point passes. The result is therefore phase present without resolved
engineering-driver amplification for this matrix; thresholds are not relaxed
to obtain a positive manuscript claim.

## Phase-conditioned short-time result

The complete-cycle conclusion above is retained as a cancellation guardrail.
It does not answer the narrower question of short upward drag while the
same-column seabed pressure is still negative and recovering from its trough.
That comparison is generated from the same audited raw data with:

```bash
python3 analyze_phase_conditioned_uplift.py
```

The two registered recovery branches are `2.34--2.60 s` and
`3.64--3.90 s`. They use the five contiguous saved frames from the surface
pressure minimum to the cycle boundary while the excess pressure remains below
10% of its fitted amplitude in suction. No derivative or metric is smoothed.

| Intrinsic permeability | Local activity recovery 1 / 2 | Signed force recovery 1 / 2 | Short-time gate |
|---:|---:|---:|:---:|
| `1e-13 m2` | `-0.195% / -0.195%` | `-1.128% / -1.012%` | FAIL |
| `3e-13 m2` | `+0.070% / +0.061%` | `-0.467% / -0.416%` | FAIL |
| `1e-12 m2` | `+1.475% / +1.455%` | `-0.365% / -0.326%` | FAIL |
| `3e-12 m2` | `+6.771% / +6.668%` | `+1.274% / +1.156%` | PASS |

Thus retained phase structure produces a repeatable, resolved short-time
upward hydraulic-demand increase at `3e-12 m2`, even though no case passes the
two-cycle net-uplift gate. The lowest permeability has the largest crown delay
but insufficient transmission, so the effect is non-monotonic. This PASS is a
fixed-state hydraulic-risk mechanism, not realised liquefaction or pipe motion;
independent lagged/phase-erased SANISAND replays remain required for that claim.
