# Draft manuscript section 7.5

## 7.5 Phase-lag ablation and pipeline-support response

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
phase-only experiment.

The hydraulic trigger is the hydrostatic-corrected upward seepage-force index

`IF = max[0, (f_seep,l + rho_l g) dot (-g/|g|) / gamma_sub]`,

where `IF >= 1` means that the upward wave-induced seepage force reaches the
submerged skeleton weight. The support-zone `IF >= 1` area--time was
`3.0942 x 10^-3 m2 s` for RL and `4.08583 x 10^-3 m2 s` for RE. The phase-erased
case was therefore 32.05% larger. The maximum indices were 2.072 and 3.243,
respectively. These values are opposite to the hypothesis that retaining the
subsurface phase lag increases the hydraulic trigger in this case.

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
Requiring the same particle at the same saved time to satisfy both `IF >= 1`
and `R_sigma <= 0.05` gives joint area--times of only
`5.19624 x 10^-5 m2 s` for RL and `2.07851 x 10^-4 m2 s` for RE. Joint
occurrence was active in 4 saved RL frames and 14 saved RE frames and was four
times larger in the phase-erased case. Moreover, only 2 of 51 RL support-zone
stress-loss particles and 2 of 43 RE particles experienced `IF >= 1` during
the preceding 1.3-s wave cycle. The lagged case therefore accumulated somewhat
more stress-loss area when that metric was integrated separately, but it did
not establish the registered hydraulic-trigger-to-stress-loss chain.

Pipeline response was correspondingly indistinguishable at engineering scale.
The maximum `|uy|/D` was 0.002914 for RL and 0.002912 for RE, the maximum
absolute rotations were `2.487 x 10^-4` and `2.455 x 10^-4 rad`, and both
post-release no-contact fractions were 0.001996. Neither case approached the
registered large-deformation threshold `max|uy|/D >= 0.5`. Equal small pipeline
responses do not prove that the local fields are identical, but they provide no
engineering-scale support for easier liquefaction caused by the retained lag.

The defensible result of this screen is therefore a null/opposite finding: the
prescribed phase lag was removed successfully at the database level, but the
PIC-smoothed pressure-response similarity gate failed; RE showed the larger IF
and same-particle joint occurrence; and RL and RE produced nearly identical
pipeline motion. The statement that phase lag makes liquefaction easier is not
supported by this registered case. The larger separate RL stress-loss integral
is reported as stress-path redistribution, not as proof of phase-lag-induced
liquefaction.

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

This negative result is informative for the broader phase-lag argument. Phase
delay can alter the spatial gradient and stress path, but its sign and severity
cannot be inferred from lag magnitude alone. Pressure penetration, smoothing,
boundary conditions and the relative phase of neighbouring material points
jointly control the upward gradient. A confirmatory study would therefore need
a redesigned control that preserves the response amplitude after PIC smoothing,
followed by the same IF, `R_sigma`, joint-occurrence and pipeline-response gates.
No such additional production calculation is included here.

## Figure/caption text to use

### Figure 7 -- physical saturation context and hydraulic bridge

LS and HS two-dimensional fields at the common registered `t_IF` frame, followed
by surface/crown/shoulder/invert pressure histories and fitted amplitude--phase
metrics. LS--HS is a fully coupled saturation contrast in which storage,
amplitude, phase and drainage vary together; it is not a phase-only comparison.
The pipe outline uses its instantaneous position and the orange region marks the
pre-registered support cohort.

### Figure 8 -- two-dimensional phase-ablation fields

RL, RE and their Lagrangian ID-matched difference at three common registered
times (`t_IF`, `t_joint` and `t_Rsigma`). Rows show checkpoint-relative excess
mixture pressure, hydrostatic-corrected upward seepage index `IF`, eligible
vertical-stress remaining ratio `R_sigma`, and the same-particle joint indicator
`IF >= 1 and R_sigma <= 0.05`. Coordinates and colour limits are shared within
each comparison, and each pipe outline uses the corresponding instantaneous
pose. RE is a one-way numerical counterfactual. The panels show the opposite or
mixed response and must not be captioned as evidence that lag makes
liquefaction easier.

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
`analysis/material_point/cyclic_constitutive_comparison.audit.json`. The final
claim wording is also generated in `analysis/screen/manuscript_evidence.md`.
Values must not be replaced by visual estimates from the figures.
