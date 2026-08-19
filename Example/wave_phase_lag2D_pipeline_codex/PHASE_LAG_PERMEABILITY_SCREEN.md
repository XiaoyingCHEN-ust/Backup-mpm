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

### Pre-registered upper-bound sensitivity (2026-08-18)

The four-point matrix showed a repeatable short-time upward-force advantage at
its upper point during both negative-pressure recovery branches, while the
complete-cycle result remained null.  Before running any additional solver
case, two upper-bound sensitivity points were fixed at `5e-12` and `7e-12 m2`.
Both retain `Sw=0.94`, the `0.12 m` wave height, `1.3 s` period, geometry,
constitutive parameters, time grid, fixed pipe and both solver smoothing flags
exactly as above.  Both points are run regardless of the first outcome.

The extension is successful only if the *raw* nearest-neighbour calculation
keeps the local upward-activity and signed-support differences positive in
both registered recovery branches.  A one-lattice-layer local affine gradient
fit is reported only as a secondary spatial-reconstruction sensitivity; it
must reproduce an affine pressure field exactly and agree in effect direction.
A directional decomposition separately integrates `|f_x|` and `|f_y|` over
the same fixed support cohort.  The hypothesised redistribution is supported
only when the raw calculation shows smaller horizontal activity and larger
vertical activity for the lagged field, with the one-layer reconstruction
preserving both signs.  This direction check is reported independently of the
short-time gate and cannot rescue a failed raw upward-force comparison.
A two-layer fit is prohibited for inference because the pre-run diagnostic
changed the signed-effect direction, demonstrating excessive spatial
averaging.  No solver pressure smoothing, time smoothing, window change or
threshold relaxation may be used to turn a failed raw result into a pass.

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
- Any raster smoothing is display-only (`sigma=2.0 pixels` by default). Raw
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
  --label LABEL --display-smoothing-sigma 2.0
```

`CASE_LABEL` is an identity field, not a free-form caption: it must equal
`LABEL.upper()`. `UUID_TOKEN` must equal `CASE_LABEL.replace("-", "_")`, and
the runner UUIDs must be `PLP_EXP_UUID_TOKEN_EQ` and
`PLP_EXP_UUID_TOKEN_HD`. For example, `LABEL=k3e-13_sw094_nosmooth` requires
`CASE_LABEL=K3E-13_SW094_NOSMOOTH` and
`UUID_TOKEN=K3E_13_SW094_NOSMOOTH`. Do not point the phase plot at `02_WAVE` or
an `_WAVE` result when producing the current HD/pressure-control evidence.

## Current audited status (2026-08-18)

| Intrinsic permeability | Status | Crown lag / amplitude ratio | Pressure-only engineering screen |
|---:|---|---:|---|
| `1e-13 m2` | EQ/HD/phase/2-D and raw-gradient audit complete | `+49.90 deg / 0.132` | failed; net uplift `-0.0239%`, joint zero |
| `3e-13 m2` | EQ/HD/phase/2-D and raw-gradient audit complete | `+45.92 deg / 0.157` | failed; net uplift `-0.0131%`, joint advantage zero |
| `1e-12 m2` | EQ/HD/phase/2-D and raw-gradient audit complete | `+24.00 deg / 0.193` | failed; net uplift `-0.1399%`, joint advantage zero |
| `3e-12 m2` | EQ/HD/phase/2-D and raw-gradient audit complete | `-5.55 deg / 0.163` | failed; crown phase unresolved, net uplift `-0.9038%`, IF area-time `-0.8552%` |
| `5e-12 m2` | upper-bound EQ/HD/phase/2-D and raw-gradient audit complete | `-15.74 deg / 0.146` | failed; net uplift `+0.1226%` is below 5% |
| `7e-12 m2` | upper-bound EQ/HD/phase/2-D and raw-gradient audit complete | `-22.28 deg / 0.140` | failed; net uplift `+1.9125%` is below 5% |

An earlier `1e-13` run remains useful only as a reproducibility reference. Its
original v2 audit reported `+50.20 deg` from the `02_WAVE` / `_WAVE` solver
output, whereas the later phase-v3 audit reported `+50.02 deg` from the
distinct `03_HD` / `_HD` pressure-recording output used by the current driver
analysis. These values are from different solver cases, not alternative
rounding of one result; `+50.02 deg` is the applicable historical HD value.
That run failed the same net-uplift interpretation, and its historical EQ
configuration identity was not frozen by the current runner schema. It is
therefore excluded from the primary four-point summary.
No point passes the complete-cycle engineering gate. The result is therefore
phase present without a resolved sustained net-uplift amplification; thresholds
are not relaxed to obtain a positive claim. The distinct phase-conditioned
short-time result below addresses recovery-stage vertical forcing and does not
override this guardrail.

## Phase-conditioned short-time result

The complete-cycle conclusion above is retained as a cancellation guardrail.
It does not answer the narrower question of short upward drag while the
same-column seabed pressure is still negative and recovering from its trough.
That comparison is generated from the same audited raw data with:

```bash
python3 analyze_phase_conditioned_uplift.py \
  --labels k1e-13_sw094_nosmooth_fresh k3e-13_sw094_nosmooth \
  k1e-12_sw094_nosmooth k3e-12_sw094_nosmooth \
  k5e-12_sw094_nosmooth k7e-12_sw094_nosmooth \
  --output-dir analysis/phase_conditioned_uplift_upper_sensitivity_sigma2 \
  --display-smoothing-sigma 2.0
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
| `5e-12 m2` | `+7.040% / +6.960%` | `+5.231% / +4.799%` | PASS |
| `7e-12 m2` | `+8.466% / +8.390%` | `+9.344% / +8.617%` | PASS |

The upper-bound extension also resolves the direction of the effect. In the
second recovery, raw `|f_x|` activity changes by `-3.126%` and `-3.770%` at
`5e-12` and `7e-12 m2`, whereas raw `|f_y|` activity changes by `+7.331%` and
`+8.352%`. The one-layer local-affine sensitivity retains all four signs. The
result is therefore a shift from horizontal toward vertical forcing during the
negative-pressure recovery, not a global increase of all gradient components.

Thus retained phase structure produces a repeatable, resolved short-time
upward hydraulic-demand increase at `3e-12 m2` and at both pre-registered upper
sensitivity points, with the strongest directional result at `7e-12 m2`, even
though no case passes the two-cycle net-uplift gate. The lowest permeability
has the largest crown delay but insufficient transmission, so the effect is
non-monotonic. These PASS results identify a fixed-state hydraulic-risk
mechanism, not realised liquefaction or pipe motion; independent lagged and
phase-erased SANISAND replays remain required for that claim.

## Cross-parameter physical reference

The same-parameter phase-erased comparison above is the phase-component
counterfactual. A second, deliberately different question compares two fully
physical parameter combinations: lower permeability/lower saturation versus
higher permeability/higher saturation. This brackets an engineering condition
but cannot attribute the complete difference to phase lag because both `k` and
`Sw` change.

The main pair is `k=7e-12 m2, Sw=0.94` versus
`k=9.79e-12 m2, Sw=0.993`; the `k=1e-13 m2, Sw=0.94` case is an extreme
low-transmission sensitivity. All remaining HD configuration fields agree
after removing only those two parameters and case identity paths. Both solver
smoothing switches remain false. The common external wave is unchanged and the
measured surface histories differ by only `0.0196%` for the main pair and
`0.0171%` for the extreme pair.

| Lower-k/lower-Sw case | Signed impulse change, recovery 1 / 2 | Local positive-part change, recovery 1 / 2 | Raw vertical/horizontal ratio change, recovery 1 / 2 |
|---:|---:|---:|---:|
| `7e-12 m2, 0.94` | `+221.5% / +204.8%` | `-71.5% / -72.6%` | `+3.41% / +1.21%` |
| `1e-13 m2, 0.94` | `+438.0% / +428.3%` | `-36.8% / -40.0%` | `+56.0% / +52.1%` |

The signs reconcile the earlier apparently opposite interpretations. The
higher-`k`/higher-`Sw` reference transmits a much larger absolute pressure-
gradient field, so its integral of local upward positive parts is larger. It
also contains strong alternating upward/downward and lateral layers around the
pipe. Spatial cancellation leaves signed upward impulses of only `8.361` and
`9.501 N s/m`, compared with `26.879/28.955 N s/m` for `7e-12` and
`44.980/50.190 N s/m` for `1e-13`. The one-layer affine sensitivity strengthens
the directional result: the vertical/horizontal ratio rises by
`17.0%/13.4%` for `7e-12` and `66.6%/59.9%` for `1e-13`.

The reference is not labelled no-lag after inspecting the field. Its same-
column shoulder difference is `-3.72 deg`, but crown and invert differences
are `-78.42` and `-101.84 deg`; the lower-`k` cases also do not have a larger
phase magnitude at every probe. The defensible conclusion is therefore that
the lower-`k`/lower-`Sw` combinations redistribute the short recovery-stage
gradient toward a more coherent upward resultant, which increases potential
hydraulic risk. It is not proof that phase lag alone caused the difference or
that liquefaction occurred.

```bash
python3 analyze_phase_lag_parameter_contrast.py
python3 analyze_phase_lag_parameter_contrast.py \
  --lower-label k1e-13_sw094_nosmooth_fresh \
  --higher-label k9p79e-12_sw0993_nosmooth \
  --output-dir analysis/phase_lag_parameter_contrast_extreme_low_k
```

Raw fields determine every number. A one-particle-layer affine reconstruction
and Gaussian `sigma=3 pixels` are used only for the two-dimensional display;
the alternating bands remain visible rather than being smoothed out of the
evidence. The display contains six panels at two recoveries one period apart.
Colour gives signed `fy/gamma'`, while arrows give the complete
`(fx,fy)/gamma'` vector. A single audited `q97.5` arrow-length cap is shared by
all panels, preserves direction, and is not used in any reported metric.
