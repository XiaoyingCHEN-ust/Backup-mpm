# Draft manuscript section 7.5

## 7.5 From hydraulic phase lag to pipeline-support loss: phase and constitutive ablations

Sections 7.1--7.4 showed that pore-pressure attenuation and phase lag are
coupled outcomes of storage and seepage, and that a delayed signal may generate
stronger local pressure gradients even when its transmitted amplitude is
smaller.  The engineering consequence is examined here using a shallow empty
pipeline rather than another permeability--saturation sweep.  The pipe provides
direct observables of support loss--normalised displacement, rotation and
soil-contact loss--while retaining the hydraulic mechanism established above.
The pipe diameter is 0.12 m and the cover-to-diameter ratio is 0.25, placing the
invert 0.15 m below the seabed within the intermediate-permeability zone where
the preceding analysis predicts that delayed pressure can still penetrate to
engineering depth.

The comparison matrix separates three questions.  LS (`Sw=0.993`) and HS
(`Sw=0.94`) are fully coupled SANISAND simulations and show the physical effect
of changing saturation; phase, amplitude, storage and drainage therefore vary
together.  RL and RE are one-way replays driven by the same fixed-pipe pressure
database.  RE rotates only each material point's fitted fundamental pressure
component to the phase of its local seabed-surface point while preserving the
point mean, fundamental amplitude, residual and higher-harmonic content.  The
RL--RE comparison is consequently the phase-only diagnostic, but RE is a
numerical counterfactual rather than a second fully coupled solution.  Finally,
RL (SANISAND) and RM (Mohr--Coulomb) receive the identical lagged pressure
history and form the principal constitutive ablation, subject to a first-frame
`p'--q` initial-state audit because RM includes the MC_EQ handoff.  HS--HM is retained only
as fully coupled context because their pore pressures may diverge with their
solid responses.

All cases begin from a fixed-pipe, zero-wave equilibrium.  The two saturation
states are first equilibrated for 4 s with `LinearElastic2D`, pure PIC and a
Cundall damping coefficient of 5 s^-1.  The elastic tangent is the same local
reference tangent used for the Mohr--Coulomb ablation (`E=9.10081 MPa`,
`nu=-0.0079625`), while a separate 23.8 MPa modulus is used only for the
explicit time-step bound.  A checkpoint is accepted only if all particles
remain in the domain, all fields are finite, `max|v| <= 10^-3 m/s`, maximum
displacement is below one background cell and porosity is valid.  SANISAND
restores this checkpoint directly.  The zero-cohesion Mohr--Coulomb cases pass
through an additional 4 s fixed-pipe, pure-PIC, damped MC relaxation and the
same unchanged stability gate before dynamic loading.  This MC_EQ stage is a
numerical handoff and is not interpreted as a comparison result.

### 7.5.1 A phase-only test of the hydraulic trigger

Figure 7 first verifies the control.  Over the pre-release fitting window, the
RL and RE pressure means differ by no more than
`{{RL_RE_MAX_NORMALISED_MEAN_DIFFERENCE}}`, and their fitted fundamental
amplitudes differ by no more than `{{RL_RE_MAX_RELATIVE_AMPLITUDE_DIFFERENCE}}`.
The pressure-database audit gives zero change in the fitted point means and a
maximum amplitude change of `{{TRANSFORM_MAX_AMPLITUDE_CHANGE_PA}} Pa`, whereas
the mean crown--shoulder--invert phase lag decreases by
`{{RL_RE_MEAN_PHASE_LAG_REDUCTION_DEG}} degrees`.  These checks are prerequisites
for interpreting RL--RE as a timing ablation rather than an amplitude ablation.

The hydraulic trigger is evaluated using the hydrostatic-corrected upward
seepage-force index

`IF = max[0, (f_seep,l + rho_l g) dot (-g/|g|) / gamma_sub]`,

where `IF >= 1` denotes that the wave-induced upward seepage force reaches the
submerged skeleton weight.  The support-zone area--time above this threshold is
`{{RL_SUPPORT_IF_AREA_TIME_M2S}} m2 s` in the lagged replay and
`{{RE_SUPPORT_IF_AREA_TIME_M2S}} m2 s` in the phase-erased replay; the
corresponding maximum indices are `{{RL_MAX_IF}}` and `{{RE_MAX_IF}}`.

Use exactly one of the following result-dependent interpretations after the
evidence audit:

- If `phase_lag_hydraulic_trigger_supported=true`: Although the fundamental
  pressure amplitude is matched, retaining the measured subsurface phase lag
  increases the duration and/or spatial extent of `IF >= 1`.  The destabilising
  effect therefore arises from the phase-dependent pressure gradient, not from
  a larger pressure amplitude.
- If the gate is false: The phase-only comparison does not demonstrate a
  stronger hydraulic trigger for the lagged signal in this case.  The result
  must be reported as unresolved or opposite rather than inferred from the
  saturation comparison.

### 7.5.2 From hydraulic triggering to realised skeleton-stress loss

Hydraulic triggering is not, by itself, evidence that the soil skeleton has
liquefied.  Realised stress loss is therefore evaluated independently from

`R_sigma = sigma'_v(t) / sigma'_v(0)`,

using only points with `|sigma'_v(0)| >= 100 Pa`; `R_sigma <= 0.05` denotes at
least 95% loss of the initial vertical effective stress.  The support-zone
area--time satisfying this condition is
`{{RL_SUPPORT_RSIGMA_AREA_TIME_M2S}} m2 s` for RL and
`{{RE_SUPPORT_RSIGMA_AREA_TIME_M2S}} m2 s` for RE, with minimum ratios
`{{RL_MIN_RSIGMA}}` and `{{RE_MIN_RSIGMA}}`.  The conventional pore-pressure
ratio `ru` is retained only in probe histories and is not used to define the
trigger, affected area or duration.

- If `phase_lag_realised_liquefaction_supported=true`: The increased
  phase-controlled hydraulic trigger is accompanied by more extensive realised
  skeleton-stress loss.  This closes the causal chain from subsurface phase lag,
  through the excess pressure gradient, to liquefaction rather than pressure
  timing alone.
- If the gate is false: Retaining the lag increased, at most, the hydraulic
  potential; it did not produce a materially larger `R_sigma <= 0.05` response.
  The manuscript must stop at that narrower conclusion.

The engineering response completes this chain.  For RL and RE, report the
per-cycle pipeline uplift `uy/D`, rotation, contact number, breakout time and
post-release no-contact fraction beside the `IF` and `R_sigma` histories.  A
single final-displacement bar is insufficient because an equal endpoint can
hide different trigger timing, cyclic accumulation and contact loss.  Report
whether either case reaches the pre-registered large-deformation gate
`max|uy|/D >= 0.5` or breakout with documented contact-network loss.  Soil
`max|u|/h` and the fraction of material points with `|u| >= h` should be shown
separately for the full bed and the initial pipeline-support cohort.  These
quantities document material crossing of the background grid; they are not
strain measures and do not establish that a finite-element treatment is
impossible.

### 7.5.3 Why cyclic state evolution matters

The RL--RM comparison uses the identical lagged liquid/gas pressure database,
so the principal difference is the skeleton constitutive description.
Mohr--Coulomb uses `c=0`, `psi=0` and `phi=31.1474 degrees`, matching the
SANISAND critical-state slope `Mc=1.25`.  Its constant elastic tangent is
locally matched to the implemented SANISAND tangent at `p'=3 kPa` and
`n=0.485`; this removes a first-order stiffness bias but does not make the
models equivalent because SANISAND remains pressure and state dependent.
Because RM passes through MC_EQ while RL restores the elastic checkpoint
directly, the crown, shoulder and invert first-frame `p'` and `q` values must
also agree within the pre-registered 5% tolerance.  If this audit fails, the
comparison includes an initial-state/handoff-path difference and cannot be
presented as a clean constitutive isolation.

Figure 9 should compare representative crown, shoulder and invert `q-p'` paths,
not only final pipe motion.  SANISAND additionally reports accumulated
equivalent plastic strain `eps_p_q`, void-ratio/state evolution and the cyclic
stress path, whereas the Mohr--Coulomb comparison reports `pdstrain`.  The
matched-pressure results are
`{{RL_SUPPORT_RSIGMA_AREA_TIME_M2S}} m2 s` and
`{{RM_SUPPORT_RSIGMA_AREA_TIME_M2S}} m2 s` for the support-zone stress-loss
area--time, with maximum `|uy|/D` of `{{RL_MAX_UY_OVER_D}}` and
`{{RM_MAX_UY_OVER_D}}`, respectively.

- If `sanisand_mechanistic_advantage_supported=true`: Under the same external
  pressure history, the SANISAND state variables evolve cyclically and the
  stress-loss and/or engineering path differs materially from Mohr--Coulomb.
  The advantage demonstrated here is the ability to resolve contraction,
  pressure-dependent stiffness, fabric/state memory and cyclic mobility that a
  perfect-plastic Mohr--Coulomb model cannot represent.  It is not a claim that
  SANISAND must always predict a larger final displacement.
- If the SANISAND state evolves but the field response is not materially
  separated: The calculation demonstrates additional resolved mechanisms but
  not a material engineering consequence for this loading duration.  Do not
  describe model complexity alone as predictive superiority.
- If neither state evolution nor response separation is resolved: The field
  comparison is inconclusive; retain the single-material-point preflight and
  avoid a superiority claim.
- If `initial_state_QA.passed=false`: Describe RL--RM as a comparison of the
  complete SANISAND and MC handoff/model chains.  Do not attribute the response
  difference uniquely to cyclic constitutive evolution.

Strictly, this numerical ablation demonstrates a *mechanistic modelling
advantage*.  A claim of greater predictive accuracy would require independent
cyclic laboratory or field response data for the same soil state, which are not
available in this case study.

### 7.5.4 Implications and limitations

The pipeline calculation isolates seabed-mediated support loss.  The pipe is
loaded by self-weight, submerged buoyancy and soil contact, but direct
oscillatory wave pressure and drag on an exposed pipe are not included.
Post-breakout motion must therefore be interpreted as an idealised consequence
of support loss rather than complete wave--pipe fluid--structure interaction.
Similarly, RE deliberately breaks full two-way conservation and must always be
labelled a one-way numerical counterfactual.

Subject to the result gates above, the case study advances the hydraulic result
of Sections 7.1--7.4 in two steps.  First, matching amplitude while removing the
subsurface phase delay tests whether phase lag itself increases the critical
upward gradient and realised effective-stress loss.  Second, matching the
lagged pressure history across two constitutive models tests whether cyclic
state evolution changes the post-trigger path.  This design supports a more
specific conclusion than either a saturation sweep or a final-displacement
comparison alone: phase lag can change *when and where* hydraulic forcing
becomes critical, while the constitutive model controls how that trigger is
converted into stress loss, cyclic accumulation and pipeline migration.

## Figure/caption text to use

### Figure 7 -- case and hydraulic bridge

Geometry and registered comparisons for the shallow pipeline case, followed by
surface/crown/invert mixture-pressure histories and fitted amplitude/phase.
LS--HS is a fully coupled saturation contrast.  RL--RE is a one-way numerical
counterfactual in which the fitted fundamental phase lag is removed while each
point's mean, fundamental amplitude and residual/higher-harmonic signal are
retained.  MC_EQ is a fixed-pipe numerical handoff and is not a comparison
result.

### Figure 8 -- phase-controlled trigger-to-response chain

RL and RE at the same registered wave phase: hydrostatic-corrected upward
seepage-force index `IF`, vertical effective-stress remaining ratio `R_sigma`,
soil displacement/contact state and pipeline `uy/D`.  The primary thresholds
are `IF >= 1` and `R_sigma <= 0.05`; `ru` is diagnostic only.  RE is a one-way
numerical counterfactual, not a fully coupled physical prediction.

### Figure 9 -- constitutive and large-deformation response

Matched-pressure RL (SANISAND) and RM (Mohr--Coulomb) histories: representative
`q-p'` paths, `eps_p_q`/`pdstrain`, support-zone stress-loss area, cyclic pipe
uplift/rotation and contacts.  Soil `max|u|/h` and `fraction(|u| >= h)` document
background-grid crossing and are not interpreted as strain.  SANISAND's
advantage is assessed from cyclic path/state evolution and its observed
response consequence, not from an assumed ordering of final displacement.

## Placeholder source map

All placeholders above must be populated from `analysis/<group>/manuscript_evidence.json`
or the corresponding case summary; do not fill them from visual estimates.

| Placeholder family | JSON section |
|---|---|
| RL--RE pressure QA and lag reduction | `phase_only_RL_RE.pressure_control_QA`, `phase_lag_contrast` |
| `IF` area--time and maxima | `phase_only_RL_RE.hydraulic_trigger` |
| `R_sigma` area--time and minima | `phase_only_RL_RE.realised_skeleton_stress_loss` |
| SANISAND/MC path and response | `constitutive_RL_RM` |
| permitted conclusion | `claim_gate` and generated `manuscript_evidence.md` |
