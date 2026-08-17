# Draft manuscript section 7.5

## 7.5 Phase-lag ablation in a shallow-buried pipeline: from local gradient focusing to engineering response

Sections 7.1--7.4 showed that pressure attenuation and phase lag are coupled
outcomes of storage and seepage. Here, a shallow empty pipeline is used to test
whether the lag itself increases hydraulic triggering and realised loss of
skeleton stress. The pipe diameter is 0.12 m and the cover-to-diameter ratio is
0.25. The invert is therefore 0.15 m below the seabed, within the zone in which
the preceding simulations showed measurable wave-pressure penetration.

The comparison matrix separates a physical saturation contrast from a
phase-only numerical ablation. LS (`Sw=0.993`) and HS (`Sw=0.94`) are fully
coupled SANISAND calculations, so phase, amplitude, storage and drainage change
together. RL and RE are one-way replays derived from the same fixed-pipeline
pressure database. At each material point, RE replaces the fitted fundamental
phase by that of the corresponding seabed-surface point while retaining the
point mean, fundamental amplitude, residual and higher-harmonic content. RE is
therefore a numerical counterfactual, not a second fully coupled prediction.

All accepted dynamic cases start from a fixed-pipeline, zero-wave equilibrium.
The LS and HS beds were equilibrated for 4 s using `LinearElastic2D`, pure PIC,
`dt=10^-4 s` and Cundall damping of 5 s^-1. A checkpoint was accepted only when
all 7,376 particles remained inside the grid, all audited fields were finite,
`max|v| <= 10^-3 m/s`, maximum displacement was below one background cell,
porosity was in `(0,1)`, and the near-surface normal-stress audit passed.
SANISAND then restored the accepted checkpoint directly. The stability limit
was not relaxed for any case.

The planned zero-cohesion Mohr--Coulomb field comparison could not be admitted.
After restoring the HS elastic checkpoint, the current MC_EQ diagnostic reached
`max|v|=1.847 x 10^-3 m/s` at 0.5 s and therefore failed the unchanged
`10^-3 m/s` gate. The failure was concentrated around the shallow pipe and
surface, where the elastic equilibrium stress field was not compatible with the
zero-cohesion, zero-tension Mohr--Coulomb surface. HM and RM are consequently
not reported as field comparisons.

### 7.5.1 Phase-control audit and hydraulic trigger

The prescribed-database transformation passed its pointwise invariants: the
maximum fitted-mean change was 0 Pa and the maximum fundamental-amplitude change
was `5.68 x 10^-14 Pa`. Across the crown, shoulder and invert probes, the mean
absolute phase-lag reduction was 45.37 degrees (120.91, 7.14 and 8.05 degrees,
respectively), exceeding the registered 5-degree minimum.

The input invariant does not, however, imply equal solver-output amplitudes.
The checkpoint-relative, saturation-weighted mixture pressure written to VTP is
the PIC-smoothed response field. Its largest RL--RE fundamental-amplitude
difference was 54.75% at the crown, and the largest mean difference normalised
by surface amplitude was 5.953%. Both exceed the pre-registered 2% response
similarity tolerance. Thus the raw forcing transformation is correct, but the
overall pressure-response control fails. RL--RE remains useful as a numerical
counterfactual, but it cannot be interpreted as a clean, response-level
phase-only experiment. This failure also means that a gradient reconstructed
from the VTK-smoothed output cannot be the primary hydraulic evidence.

The hydraulic trigger is the hydrostatic-corrected upward seepage-force index

`IF = max[0, (f_seep,l + rho_l g) dot (-g/|g|) / gamma_sub]`,

where `IF >= 1` means that the upward wave-induced seepage force reaches the
submerged skeleton weight. To connect this local index to the engineering load
path, the signed upward wave-induced pressure-gradient force density is defined
as

`f_y^ex = [-grad(p_l) + rho_l g] dot (-g/|g|)`,

and its positive resultant in the registered pipeline-support cohort is

`F_up^+(t) = integral_Omega_sup max(f_y^ex, 0) dA`.

`F_up^+` has units of N/m for the two-dimensional unit-thickness model. For the
primary gradient bridge, `grad(p_l)` was reconstructed by a direct
nearest-neighbour finite difference of the raw prescribed-pressure database,
with particle IDs defining the original lattice and current coordinates
defining the spacing. The saved `PIC_liquid_pressures` and
`liquid_seepage_forces` were deliberately excluded because VTK-time PIC
smoothing can affect those arrays even though it is not written back to the
solver state. Current particle volumes were used in all integrals. No pressure
averaging or spatial smoothing was used for this force, `IF`, threshold area,
frame selection or any acceptance gate; the optional Gaussian filter in Figure
8 is confined to the displayed raster.

At the common registered `t_IF=7.865 s`, RL produced
`F_up^+=287.75 N/m`, compared with `265.00 N/m` for RE, an increase of 8.59%.
The corresponding maximum raw `IF` was 9.459 for RL and 7.752 for RE (+22.02%),
while the instantaneous raw `IF >= 1` area was `1.04975 x 10^-2 m2` and
`9.29698 x 10^-3 m2`, respectively (+12.91%). Thus retaining the fitted
subsurface phase structure concentrated a stronger upward pressure gradient
around the pipeline at this wave phase. This is the engineering bridge between
the local phase-lag mechanism and an instantaneous upward lifting action.

The retained-lag increase persisted in the raw full-record measures. Peak
`F_up^+` was 393.81 N/m for RL and 362.77 N/m for RE (+8.56%), and its time
integral was 2604.72 and 2458.00 N s/m (+5.97%), respectively. Raw
`IF >= 1` area--times were `8.40026 x 10^-2` and
`8.27400 x 10^-2 m2 s` (+1.53%). These results support a larger local and
integrated upward pressure-gradient action when the fitted lag is retained in
this counterfactual. They do not by themselves prove realised liquefaction or
large pipeline movement.

For transparency, the previously registered VTK-output `IF` sensitivity gave
the opposite area--time ordering: `3.0942 x 10^-3 m2 s` for RL and
`4.08583 x 10^-3 m2 s` for RE. That index was reconstructed after the
output-only PIC smoothing and is retained as a numerical-sensitivity result,
not substituted for the raw database-gradient metric. The sign reversal is
why smoothing and derivative reconstruction must be declared rather than used
silently to support either conclusion.

### 7.5.2 Realised skeleton-stress loss and joint occurrence

Hydraulic triggering alone is not labelled liquefaction. Realised loss is
evaluated independently using

`R_sigma = sigma'_v(t) / sigma'_v(0)`,

only where `|sigma'_v(0)| >= 100 Pa`; `R_sigma <= 0.05` denotes at least 95%
loss of the initial vertical effective stress. The support-zone area--time below
this threshold was `1.73474 x 10^-2 m2 s` for RL and
`1.47971 x 10^-2 m2 s` for RE. Considered alone, this metric was 17.23% larger
for RL, and the minimum ratios were -5.288 and -3.559. The conventional `ru`
ratio is retained only as a probe diagnostic and is not used to define the
trigger or realised stress loss.

The separate IF and stress-loss integrals do not establish a causal chain.
Requiring the same particle at the same saved time to satisfy both raw
`IF >= 1` and `R_sigma <= 0.05` gives joint area--times of
`6.65331 x 10^-3 m2 s` for RL and `4.46321 x 10^-3 m2 s` for RE. The
retained-lag value is 49.07% larger. Together with the 17.23% larger separate
RL stress-loss integral, this is consistent with the proposed sequence from
lagged pressure-gradient concentration to local skeleton-stress loss. It is
not a unique causal proof because the output response-amplitude gate failed and
RE is one-way. The earlier joint result calculated with the VTK-smoothed `IF`
(`5.19624 x 10^-5` versus `2.07851 x 10^-4 m2 s`) is retained only as the
derivative/smoothing sensitivity; it cannot override or be merged with the raw
definition.

Pipeline response was correspondingly indistinguishable at engineering scale.
The maximum `|uy|/D` was 0.002914 for RL and 0.002912 for RE, the maximum
absolute rotations were `2.487 x 10^-4` and `2.455 x 10^-4 rad`, and both
post-release no-contact fractions were 0.001996. Neither case approached the
registered large-deformation threshold `max|uy|/D >= 0.5`. Equal small pipeline
responses do not prove that the local fields are identical, but they provide no
engineering-scale support for easier liquefaction caused by the retained lag.

The defensible result of this screen is therefore positive at the local
hydraulic and stress-loss levels, but null at the pipe-motion level. The
prescribed phase lag was removed successfully at the database level. Relative
to RE, RL had larger raw instantaneous and integrated upward-gradient action,
a 49.07% larger raw same-particle trigger/stress-loss area--time, and a 17.23%
larger separate stress-loss integral. This supports the scoped statement that
retained lag can make local liquefaction triggering more likely around the
pipeline by concentrating the upward pressure gradient. It does not support a
claim of proportionally larger engineering displacement: RL and RE produced
nearly identical pipeline motion. The one-way construction, failed
response-amplitude gate and derivative sensitivity preclude a universal causal
claim.

### 7.5.3 What the SANISAND comparison does and does not show

Because MC_EQ did not provide an admissible field checkpoint, this study cannot
claim a matched-pressure field-response advantage or greater predictive
accuracy for SANISAND over Mohr--Coulomb. A material-point diagnostic was used
only to compare the mechanisms resolved under an identical strain history. The
registered path applied six constant-volume simple-shear cycles
(`gamma=+/-4 x 10^-4`, 2,400 increments and 12 reversals) from `p'=3 kPa` and
`n=0.485`. The Mohr--Coulomb elastic tangent was matched to the SANISAND tangent
at the initial state.

Under that common path, SANISAND reduced `p'` from 3,000 Pa to a minimum of
249.62 Pa, accumulated `eps_p_q=0.004193`, and evolved its backstress to
`max|alpha|=0.4698`. Mohr--Coulomb retained `p'=3,000 Pa` to numerical precision
and accumulated a plastic-strain measure of 0.000342. A dense supplemental
SANISAND path first contracted (`p'` minimum 2,951.54 Pa) and then dilated to
4,563.24 Pa while the fabric variable reached `max|Z|=0.7111`. These paths
demonstrate pressure-dependent stiffness, contraction--dilation transition,
backstress/fabric memory and cyclic mobility in the implemented SANISAND model.
They demonstrate modelling capability, not independent predictive superiority;
that stronger claim would require a stable matched field chain and laboratory
or field validation for the same soil state.

### 7.5.4 Implications and limitations

The pipeline calculations isolate seabed-mediated support loss. The pipe is
loaded by self-weight, submerged buoyancy and soil contact, while direct wave
pressure and drag on an exposed pipe are not included. RE deliberately breaks
two-way coupling and must always be described as a one-way numerical
counterfactual. All reported threshold areas are based on two-dimensional
unit-thickness material-point areas; they are not three-dimensional volumes.

This mixed result is informative for the broader phase-lag argument. Phase
delay can alter the spatial gradient and stress path, but its sign and severity
cannot be inferred from lag magnitude alone. Pressure penetration, the inherent
particle--grid transfer, boundary conditions and the relative phase of
neighbouring material points jointly control the upward gradient. Display-only
smoothing can improve the readability of a two-dimensional raster, but it must
not enter any computed force, threshold, integral or acceptance gate. A
separate exploratory run at `k=1.0 x 10^-13 m2` and `Sw=0.94`, with solver
pressure smoothing disabled, resolved a 50.20-degree crown lag with a harmonic
`R2=0.984`; all 7,376 particles remained in the grid and the maximum velocity
over the 3.9-s fixed-pipeline replay was `1.37 x 10^-4 m/s`. However, at the
selected `t=3.705 s`, raw `F_up^+` was 847.93 N/m for the lagged replay and
849.24 N/m for the phase-erased replay (-0.15%), while the saved-record impulses
were 3086.22 and 3085.90 N s/m (+0.01%). Lowering permeability therefore left
the local phase delay clearly resolvable but did not amplify the pipeline-zone
resultant; reduced pressure penetration offset the delay. This diagnostic is not part of
the registered manuscript evidence and the pipe remained fixed.

A confirmatory study should not simply lower permeability until one selected
frame looks favourable. It should pre-register a small `Sw=0.94`, no-solver-
smoothing permeability series around `10^-13--10^-12 m2`, report the complete
phase-resolved force curves, and retain both phase lag and pressure penetration
as competing explanatory variables. Only after that fixed-pipeline screen
should an extended pressure database and released RL--RE pair be run, followed
by the same continuous-force, `IF`, `R_sigma`, joint-occurrence and
pipeline-response gates. No such additional production calculation is included
here.

## Figure/caption text to use

### Figure 7 -- physical saturation context and hydraulic bridge

LS and HS two-dimensional fields at the common registered `t_IF` frame, followed
by surface/crown/shoulder/invert pressure histories and fitted amplitude--phase
metrics. LS--HS is a fully coupled saturation contrast in which storage,
amplitude, phase and drainage vary together; it is not a phase-only comparison.
The pipe outline uses its instantaneous position and the orange region marks the
pre-registered support cohort.

### Figure 8 -- phase-lag-to-pipeline gradient-force bridge

(a,b) RL and RE upward excess seepage-force index fields at the common
registered `t_IF=7.865 s`; open markers identify raw particles satisfying
`IF >= 1`. (c) Lagrangian ID-matched `IF_RL - IF_RE`. (d) Raw support-zone
positive upward pressure-gradient resultant `F_up^+`; (e) its RL--RE
difference; and (f) pipeline vertical displacement normalised by diameter. The
field raster uses a shared colour scale and optional display-only Gaussian
smoothing, whereas all markers, force integrals, extrema, impulses and response
histories are calculated from raw pressure-database values. The figure
demonstrates stronger instantaneous and full-record raw upward-gradient action
in RL, while engineering-scale pipe motion remains nearly unchanged. RE remains
a one-way numerical counterfactual.

The multi-time `R_sigma` and same-particle joint-occurrence fields should be
retained as a supplementary figure rather than combined with this mechanism
bridge.

### Figure 9 -- SANISAND material-point mechanism and MC handoff audit

Identical cyclic simple-shear paths for SANISAND and Mohr--Coulomb, including
`p'`, shear stress, accumulated plastic measure and SANISAND state variables,
together with the documented MC_EQ field-handoff failure. The figure shows the
additional cyclic state mechanisms resolved by SANISAND and explains why a
matched field-response claim was withheld; it is not a field validation or a
claim of universal predictive superiority.

## Reproducible evidence sources

All values in this section come from
`analysis/screen/manuscript_evidence.json`, the audited RL/RE case summaries and
`analysis/material_point/cyclic_constitutive_comparison.audit.json`. The
continuous pipeline-gradient bridge is bound to the formal artifacts by
`analysis/screen/figures/figure8_pipeline_gradient_force_bridge.audit.json` and
its raw histories are in the adjacent CSV. The final claim wording is also
generated in `analysis/screen/manuscript_evidence.md`. Values must not be
replaced by visual estimates from the figures.
