# Draft manuscript section 7.5

## 7.5 Phase-lag ablation and negative-pressure uplift around a shallow-buried pipeline

Sections 7.1--7.4 showed that attenuation and phase delay are coupled outcomes
of storage and seepage. A shallow empty pipeline is used here to test the more
specific engineering hypothesis that retaining the subsurface liquid-pressure
phase structure increases the upward hydraulic driver near the pipe, especially
during the short recovery from a negative seabed-pressure trough. The pipe
diameter is 0.12 m and its cover-to-diameter ratio is 0.25. The comparison is
therefore placed in the near-surface zone where the preceding calculations
showed both pressure transmission and phase variation.

The fixed exploratory permeability screen uses four intrinsic permeabilities,
`1e-13`, `3e-13`, `1e-12` and `3e-12 m2`, all at `Sw=0.94`. The four values,
qualification criteria and effect gates were fixed before completion of the
four-point matrix. Both solver pressure-smoothing switches are false. Gaussian smoothing
with `sigma=2.0 pixels` is permitted only when rasterising a figure; no
smoothing enters a derivative, threshold, integral, extremum, frame selection
or gate decision.

After the four-point matrix had been closed, `5e-12` and `7e-12 m2` were
pre-registered as an upper-bound sensitivity extension. Both points retain
`Sw=0.94`, the 0.12-m wave height, 1.3-s period, geometry, constitutive
parameters, time grid and fixed pipe, and both are run regardless of the first
extension result. This separates a declared sensitivity test from post-hoc
parameter selection.

### 7.5.1 Initial state and numerical controls

Every current-runner case starts from a fixed-pipeline, zero-wave equilibrium.
The exploratory checkpoints were equilibrated for 0.5 s using
`LinearElastic2D`, pure PIC, `dt=1e-4 s` and Cundall damping of `5 s-1`. This
short checkpoint is a parameter-screening state rather than the registered
40,000-step production equilibrium. It is accepted only if all 7,376 particles
remain in the grid, the HDF5 and VTP states agree by particle ID, all audited
fields are finite, porosity is in `(0,1)`, the near-surface normal-stress audit
passes, and `max|v|<=1e-3 m/s`. The velocity gate is unchanged throughout the
study.

The dynamic driver restores the accepted checkpoint into SANISAND, switches
to APIC, removes damping, and keeps the pipe fixed. It runs for 3.9 s: the first
1.3-s wave period is a ramp and the following two periods are the fixed
effect window. Each accepted run contains exactly 60 material-point frames and
301 raw prescribed-pressure frames. Hash-bound runner records verify every
configuration, output frame, pressure database, solver binary and dependency.

The phase-erased database is a one-way numerical counterfactual. At each
particle it rotates only the fitted wave-frequency component to the phase of
the registered top-surface particle in the same column. It preserves the
pointwise mean, fundamental amplitude, residual and higher harmonics, and the
progressive along-wave phase. The transform is rejected unless it changes a
nonzero population, preserves mean and amplitude to numerical precision, and
reduces the residual local phase difference to the registered tolerance.
It does not represent a second two-way-coupled physical solution.

### 7.5.2 Phase, penetration and two-dimensional fields

A local phase estimate is eligible only when the harmonic fit has `R2>=0.8`
and its fundamental amplitude is at least 5% of the same-column surface
amplitude. The crown is qualified as phase delayed only when the absolute
same-column difference is at least 18 degrees. Amplitude transmission and
phase must be read together: a large delay with weak transmission need not
produce a large pressure gradient or resultant force.

The current-runner results are:

| Intrinsic permeability (m2) | Crown phase (deg) | Crown/surface amplitude | Phase-qualified | Two-cycle net-uplift change | Screen |
|---:|---:|---:|:---:|---:|:---:|
| `1e-13` | `+49.90` | `0.132` | yes | `-0.0239%` | FAIL |
| `3e-13` | `+45.92` | `0.157` | yes | `-0.0131%` | FAIL |
| `1e-12` | `+24.00` | `0.193` | yes | `-0.1399%` | FAIL |
| `3e-12` | `-5.55` | `0.163` | no | `-0.9038%` | FAIL |
| `5e-12` | `-15.74` | `0.146` | no | `+0.1226%` | FAIL |
| `7e-12` | `-22.28` | `0.140` | yes | `+1.9125%` | FAIL |

The two-dimensional fields show a resolved phase structure within and around
the shallow pipe cavity. At `3e-13 m2`, for example, the crown lags its local
surface reference by 45.92 degrees while transmitting 15.7% of the surface
amplitude. At `1e-12 m2` the crown lag decreases to 24.00 degrees and the
amplitude ratio rises to 0.193. At `3e-12 m2`, phase remains spatially
structured elsewhere in the bed, but the crown difference of 5.55 degrees is
below the fixed qualification threshold. The snapshots, amplitude maps
and phase maps use identical geometry and an amplitude/fit-quality mask; their
colour interpolation is a display operation only.

These fields demonstrate why lag magnitude alone is not an engineering
severity metric. Lower permeability can increase local delay while attenuating
the pressure that reaches the pipe. Across the completed points, increasing
permeability generally reduces the crown phase difference, whereas amplitude
transmission is non-monotonic. The useful gradient therefore depends jointly
on amplitude, neighbouring-point phase and spatial support, not on a single
probe angle.

### 7.5.3 Fixed pressure-only engineering screen

The raw hydrostatic-corrected upward force density is evaluated from the
prescribed liquid-pressure database at the saved material-point coordinates.
For the registered pipeline-support cohort, the spatially integrated signed
resultant is `F_signed(t)`. The net upward part is
`F_net+(t)=max[F_signed(t),0]`. This differs from the local positive-part
activity `integral max(f_y,0)dA`, which measures simultaneous upward hotspots
even when they are cancelled by downward forcing elsewhere. The latter is a
mechanism diagnostic and cannot establish a larger net lifting action.

Hydraulic triggering is reported with

`IF=max[0,(f_seep,l+rho_l g) dot (-g/|g|)/gamma_sub]`,

where `IF>=1` means that the upward wave-induced seepage force reaches the
submerged skeleton weight. A shared-state co-location diagnostic also records
where `IF>=1` overlaps `R_sigma<=0.05`, but both pressure alternatives use the
same HD stress field in this inexpensive screen. It is therefore a spatial
screen, not evidence that retained lag caused additional skeleton-stress loss.

The screen passes only if both post-ramp periods have positive signed and net
lagged-minus-erased impulse differences; the combined net-uplift impulse is at
least 5% larger; `IF>=1` area-time has the positive direction; and a nonzero,
at-least-5% shared-state co-location advantage persists above the two-particle
area resolution for two consecutive frames. The rule is evaluated over the
complete 1.3--3.9 s window. A favourable single frame or local positive-part
increment cannot rescue a failed net gate.

At the three completed non-fresh points, none passes. Their combined two-cycle
net-uplift changes are `-0.0131%`, `-0.1399%` and `-0.9038%` for increasing
permeability. The corresponding local positive-part changes are `+0.0415%`,
`+0.0674%` and `+0.9127%`; these small positive values describe redistribution
of local hotspots rather than a larger net load. The `IF>=1` area-time changes
are `+0.1504%`, `+0.2058%` and `-0.8552%`. Shared-state joint advantages are
0%, 0% and `+0.8982%`, respectively, and do not meet the magnitude or temporal
resolution gate.

For the fresh `1e-13 m2` point, the audited combined changes are `-0.0239%`
for net uplift, `+0.0239%` for local positive-part activity and `+0.3377%` for
`IF>=1` area-time. The shared-state co-location is zero for both pressure
histories, so its relative change is undefined rather than a fabricated zero
percentage. This point also fails the screen.

The resulting complete-cycle statement is deliberately narrow: across the four
fixed permeability points and two pre-registered upper-bound sensitivities,
resolvable subsurface phase structure did not
produce a qualifying increase in the two-cycle pipeline-zone net-uplift driver.
This integral is a cancellation guardrail; it does not answer whether a short
negative-pressure sub-cycle contains a stronger upward hydraulic demand. That
phase-conditioned question is evaluated separately below. Fixed-state pressure
results cannot by themselves establish additional liquefaction,
skeleton-stress loss or pipeline motion.

### 7.5.4 Negative-pressure recovery resolves a short-time hydraulic mechanism

The short-time comparison is registered on the same two post-ramp cycles and
uses no new solver output. In each cycle, the recovery branch begins at the
minimum saved liquid pressure at the top-surface particle directly above the
pipe and ends at the cycle boundary while the surface excess pressure remains
below `-0.10` times its fitted amplitude. Each branch contains five saved
frames over 0.26 s. The first (`2.34--2.60 s`) and second
(`3.64--3.90 s`) branches sample the same physical wave phase and therefore
provide an internal repeatability check.

Two hydraulic quantities are distinguished. The local upward-activity impulse
is the time integral of `integral max(f_y,0)dA` over the registered pipeline
support zone; it measures the simultaneous upward drag acting on different soil
regions. The signed support-force impulse integrates `integral f_y dA` over the
same zone and retains spatial cancellation. A short-time mechanism is accepted
only when a phase-resolved crown, shoulder or invert probe is present, both
recovery branches increase local upward activity by at least 5%, both signed
impulse differences have the positive direction, and the final recovery has at
least a 5% same-frame local-force advantage. This gate is separate from, and
cannot override, the complete-cycle net-uplift gate.

To test the observed direction of the pressure-gradient force, the same
support cohort also records the separate spatial integrals of `|f_x|` and
`|f_y|`. The primary values use the raw nearest-neighbour derivative. A
one-particle-layer, inverse-distance-weighted local affine fit is reported only
as a spatial-reconstruction sensitivity and for the smoother two-dimensional
raster. It exactly reproduces an affine pressure field on the translated or
sheared particle lattice and is accepted only when it preserves the raw effect
directions. A two-layer fit is excluded because the pre-run diagnostic reversed
the signed effect, demonstrating excessive spatial averaging.

| Intrinsic permeability (m2) | Recovery 1 local activity | Recovery 2 local activity | Recovery 2 signed force | Short-time mechanism |
|---:|---:|---:|---:|:---:|
| `1e-13` | `-0.195%` | `-0.195%` | `-1.012%` | FAIL |
| `3e-13` | `+0.070%` | `+0.061%` | `-0.416%` | FAIL |
| `1e-12` | `+1.475%` | `+1.455%` | `-0.326%` | FAIL |
| `3e-12` | `+6.771%` | `+6.668%` | `+1.156%` | PASS |
| `5e-12` | `+7.040%` | `+6.960%` | `+4.799%` | PASS |
| `7e-12` | `+8.466%` | `+8.390%` | `+8.617%` | PASS |

At `3e-12 m2`, the final recovery local upward-activity impulses are
`118.712` and `111.291 N s/m` for the lagged and phase-erased fields,
respectively. The lagged advantage is therefore `7.421 N s/m` (`+6.668%`).
The corresponding signed impulses are `33.289` and `32.909 N s/m`
(`+1.156%`). The preceding recovery gives closely repeated increases of
`+6.771%` and `+1.274%`. The largest same-frame local advantage in the final
branch is `+7.969%`. At the common two-dimensional comparison time
`t=3.835 s`, the same-column surface excess pressure is `-168 Pa`; the raw
lagged-minus-erased local upward activity is positive near the shallow-pipe
support region. The Gaussian `sigma=2.0 pixel` filter affects raster display
only, not these numbers.

The upper-bound extension strengthens and clarifies the direction of this
short-time result. At `7e-12 m2`, the second recovery reduces raw horizontal
absolute force activity by `3.770%` while increasing vertical absolute force
activity by `8.352%`; the first recovery gives `-3.598%` and `+8.326%`.
The one-layer local-affine reconstruction preserves the signs and gives
`-3.049%` and `+16.050%` in the second recovery. Hence the stronger upward
driver is a redistribution from horizontal toward vertical pressure-gradient
forcing, not an indiscriminate increase of every gradient component and not a
result of display smoothing.

The relevant phase structure at `3e-12 m2` is not represented by the crown
alone. The crown difference is below its earlier 18-degree qualification, but
the invert retains a same-column phase difference of `-54.54 degrees`, an
amplitude ratio of `0.0899` and `R2=0.982`. This vertical spatial phase
structure is sufficient to redistribute the upward gradient during recovery.
Conversely, `1e-13 m2` has a much larger crown delay but transmits too little
amplitude to strengthen the integrated short-time driver. The completed result
therefore supports a permeability window, not the monotonic statement that a
smaller permeability always produces a larger uplift hazard.

This PASS supports only a potential transient hydraulic-risk mechanism. At
`7e-12 m2`, the complete-window `IF>=1` area-time increases by `1.357%` and
the net-uplift impulse by `1.912%`, but the latter remains below the fixed 5%
engineering gate; the paired fields also share one fixed HD skeleton state.
The result does not yet demonstrate additional liquefaction, stress loss or
pipe motion. It instead identifies
`k_s=7e-12 m2`, `Sw=0.94` and the negative-pressure recovery branch as the
first condition to test with independently integrated lagged and phase-erased
SANISAND replays. The complete-cycle negative result remains a required
guardrail against presenting a short favourable interval as a sustained net
uplift.

### 7.5.5 Direct lower-permeability/lower-saturation engineering contrast

To test the physical interpretation without relying only on the numerical
phase-erasure operation, an additional fixed-pipe reference was calculated at
`k_s=9.79e-12 m2` and `Sw=0.993`. It was compared directly with the
`k_s=7e-12 m2`, `Sw=0.94` upper sensitivity and, separately, with the fresh
`k_s=1e-13 m2`, `Sw=0.94` point. The wave, geometry, SANISAND parameters,
time grid and both false pressure-smoothing flags are identical; only
permeability, liquid/gas saturation, equilibrium state and case identity are
allowed to differ. The measured surface-pressure histories differ by less
than `0.020%`, so the applied wave is effectively matched.

This comparison resolves why the local-positive and net-force interpretations
appeared to disagree. The higher-permeability/higher-saturation reference
transmits a larger absolute pressure-gradient field and therefore has much
larger local upward positive-part activity. Nevertheless, alternating upward
and downward layers around the pipe cancel strongly. Its signed upward
support-force impulses are only `8.361` and `9.501 N s/m` in the two recovery
branches. The corresponding `7e-12 m2`, `Sw=0.94` impulses are `26.879` and
`28.955 N s/m`, increases of `221.5%` and `204.8%`. The extreme
`1e-13 m2`, `Sw=0.94` values are `44.980` and `50.190 N s/m`, or `438.0%`
and `428.3%` above the same reference.

The directional ratios support the same reading. For `7e-12 m2`, the raw
vertical-to-horizontal absolute-force ratio is `3.41%` and `1.21%` larger in
the two recoveries; the one-layer affine reconstruction gives `17.0%` and
`13.4%`. For `1e-13 m2`, the raw increases are `56.0%` and `52.1%`, and the
one-layer values are `66.6%` and `59.9%`. Thus the lower-permeability/lower-
saturation combinations carry less total gradient activity but a more
vertically biased and less cancelled upward resultant during the negative-
pressure recovery. This is the short-duration lifting-risk quantity of
interest; it is distinct from the spatial integral of local positive parts.

The reference cannot be called a no-lag solution. Its shoulder phase relative
to the same-column surface is only `-3.72 degrees`, but the crown and invert
retain differences of `-78.42` and `-101.84 degrees`. Moreover, the lower-
permeability cases do not have a larger phase magnitude at every probe. Because
both permeability and saturation change, the direct contrast is an engineering
bracket rather than a phase-lag-only causal ablation. It supports a potential
short-time hydraulic-risk mechanism, not realised liquefaction. Numerical
metrics use the raw unsmoothed pressure field; a one-layer affine
reconstruction and Gaussian `sigma=3 pixels` are confined to the displayed
two-dimensional force maps. Colour shows the signed vertical component and
arrows show the complete pressure-gradient-force vector. To prevent a few
pipe-adjacent arrows from obscuring direction elsewhere, all six panels use
one `97.5th`-percentile display-length cap; directions and every numerical
metric remain unmodified.

### 7.5.6 Relation to the earlier registered replay

An earlier registered RL--RE replay at the baseline permeability produced a
larger raw local positive-part impulse and a larger same-particle hydraulic-
trigger/stress-loss overlap when lag was retained. That result remains a
case-specific motivating diagnostic, not the primary conclusion. Its solver
response-amplitude control failed, a derivative reconstructed from the
VTK-smoothed field reversed the `IF` area-time ordering, and the released pipe
motions were nearly identical. The current-runner unsmoothed permeability screen
was introduced specifically to prevent a favourable instant or derivative
choice from determining the claim. Its net-resultant criterion has not
reproduced an engineering-driver amplification at the completed points.

The earlier `1e-13 m2` exploratory calculation is likewise retained only as a
reproducibility reference. It resolved a large crown delay and suggested the
same attenuation-versus-lag trade-off, but its historical equilibrium config
identity was not frozen by the current runner schema. It is excluded from the
primary four-point trend and replaced by the fresh current-runner case above.

### 7.5.7 SANISAND capability and the Mohr--Coulomb limitation

The planned zero-cohesion Mohr--Coulomb field comparison cannot be admitted.
After restoring the elastic HS checkpoint, the MC equilibrium diagnostic
reached `max|v|=1.847e-3 m/s` at 0.5 s and failed the unchanged
`1e-3 m/s` gate. The incompatibility is concentrated near the shallow cavity
and free surface, where the elastic fixed-pipe stress field lies outside the
zero-cohesion, zero-tension MC admissible domain. A longer run, extra damping or
a relaxed velocity threshold is not used to conceal that handoff failure.

A separate material-point comparison therefore supports only a capability
statement. Under the same six-cycle constant-volume simple-shear history from
`p'=3 kPa` and `n=0.485`, SANISAND reduces `p'` to 249.62 Pa, accumulates
`eps_p_q=0.004193`, and evolves `max|alpha|=0.4698`. The matched-tangent
perfect-plastic Mohr--Coulomb point retains `p'=3 kPa` to numerical precision
and accumulates a plastic-strain measure of 0.000342. A dense supplemental
SANISAND path first contracts and then dilates to `p'=4.563 kPa`, while its
fabric variable reaches `max|Z|=0.7111`. These paths demonstrate the implemented
model's pressure-dependent stiffness, contraction--dilation, backstress/fabric
memory and cyclic mobility. They do not demonstrate universal or independently
validated field-scale superiority.

### 7.5.8 Figures and evidence sources

**Figure 7 -- two-dimensional pressure transmission and phase.** Four common-
scale quarter-cycle liquid-excess-pressure snapshots for one explicitly
identified current-runner case, followed by the same-column amplitude ratio
and phase-difference maps. The pipe geometry and mask are identical in every
panel. Phase is hidden where amplitude is below 5% of the local surface or the
harmonic `R2` is below 0.8.

**Figure 8 -- negative-pressure recovery and permeability window.** The raw
lagged-minus-erased local upward-force history in the final post-ramp period,
followed by the two repeated recovery-branch changes in local upward activity
and signed support force, and by the raw horizontal/vertical force-activity
decomposition with one-layer sensitivity markers. The shaded interval is the registered
negative-pressure recovery rather than a visually selected instant.

**Figure 9 -- two-dimensional short-time upward hydraulic forcing.** At the
common physical time `t=3.835 s` and same-column surface excess pressure
`-168 Pa`, each of the six permeabilities occupies one row with lagged, phase-erased and
difference fields on the same geometry and colour scales. The figure shows the
non-monotonic trade-off: the largest low-permeability phase angle is attenuated,
whereas the upper-bound sensitivity points produce the strongest short
recovery-phase vertical redistribution.

**Figure 10 -- direct physical parameter-combination contrast.** The summary
separates local positive-part activity, signed upward support-force impulse and
the vertical/horizontal force-activity ratio for lower-`k`/lower-`Sw` and
higher-`k`/higher-`Sw` physical cases. A companion `2 x 3` field figure shows
two matched negative-pressure recoveries at `t=2.535` and `3.835 s`, one period
apart. Colour uses a common diverging scale for signed `fy/gamma'`; arrows show
the complete `(fx,fy)/gamma'` vector and use one audited display-only `q97.5`
length cap. The repeated red/blue layers and arrow directions expose the
spatial cancellation that makes total local activity and signed upward force
order differently. The `1e-13 m2` comparison is shown as the extreme-
attenuation sensitivity, not as the primary parameter pair.

**Supplementary permeability guardrail.** Crown phase
magnitude and amplitude transmission for the four current-runner points,
followed by raw two-cycle changes in net uplift, local positive-part activity,
`IF>=1` area-time and shared-state co-location. Each point is labelled with its
net and overall gate status. The companion two-dimensional lagged,
phase-erased and difference fields are supplementary because their selected
instant is illustrative; acceptance uses the complete two-cycle histories.

Caption sentence common to both figures:

> Raster smoothing is used for display only. All threshold markers, force
> histories, integrals, extrema and gate decisions are computed from
> unsmoothed raw pressure-database values. The complete-cycle gate uses
> 1.3--3.9 s, whereas the short-time mechanism uses two pre-defined
> trough-to-cycle-boundary negative-pressure recovery branches.

The current-runner provenance is stored under
`analysis/phase_lag_exploratory/<label>/runner_audit.json`. Each phase plot audit
binds the exact configs, results, checkpoint, analyzer, plotting script and four
PNG/PDF artifacts. Each driver audit rebinds that phase audit, the exact runner
record and both raw pressure databases. The four-point summary rejects missing,
legacy or hash-drifted inputs. Values in this section must be copied from those
audits rather than estimated visually.

The six-point phase-conditioned short-time audit, exact CSV and PNG/PDF figures
are stored under
`analysis/phase_conditioned_uplift_upper_sensitivity_sigma2/`. Its audit
revalidates all six runner, phase-v3 and full-cycle driver sources before
calculating the recovery metrics.
