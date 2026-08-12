# Manuscript section 7.5 plan

## Proposed title

**From hydraulic phase lag to post-liquefaction pipeline migration: phase and
constitutive ablations**

## Argument to make

The preceding sections establish the hydraulic mechanism: saturation and
intrinsic permeability control pressure attenuation and phase lag, which in
turn changes the timing and geometry of seepage-force instability.  Section
7.5 should demonstrate the engineering consequence without repeating the same
parameter sweep.

The shallow empty pipe is preferable to a generic trench slope because it
provides direct stability observables (`uy/D`, rotation and contact loss) and a
clear transition from pre-trigger response to post-liquefaction migration.  It
also provides an honest reason to retain MPM once material crosses cells and
contact topology changes.  Small-strain/Lagrangian FEM remains adequate before
that transition; ALE/CEL/remeshing could also model the post-failure stage, so
the paper should claim practical naturalness rather than impossibility of FEM.

## Figure sequence

1. **Case and controls.** Geometry, five probe locations, `C/D=0.25`, release
   time, and the LS/HS/HM/RL/RM/RE matrix.  Visually distinguish physical
   cases from one-way diagnostic replays.  Disclose `MC_EQ` as the fixed-pipe,
   damped Mohr-Coulomb handoff relaxation used by HM/RM, but do not plot or
   interpret it as a scientific comparison case.
2. **Hydraulic bridge.** Fixed-pipe pressure histories and fitted phase at the
   surface/crown/invert.  Show that LS is low lag, HS is high lag, and RL–RE has
   matching mean/amplitude with only the fundamental lag removed.
3. **Trigger-to-response chain.** At one common wave phase show the
   hydrostatic-corrected upward seepage-force index `IF`, actual vertical
   effective-stress remaining ratio `Rsigma`, soil displacement/contact and
   pipe `uy/D`.  Report soil `max|u|/h` and `fraction(|u|>=h)` separately as
   mesh-crossing evidence.  The two primary thresholds are `IF>=1` and
   `Rsigma<=0.05`; `ru` is auxiliary only.
4. **Engineering response.** Per-cycle `uy/D`, rotation and contacts for
   LS/HS; add selected MPM fields immediately before trigger, at breakout and
   after migration.  Report whether the pre-registered large-deformation gate
   was met and whether support-zone material crossed background cells; do not
   equate the latter with strain.
5. **Constitutive ablation.** Use RL vs RM as the main comparison because both
   receive the identical lagged pressure database; show representative `q-p'`,
   stress-loss area, cyclic displacement accumulation and
   `eps_p_q`/`pdstrain`.  HS vs HM is only fully coupled context.  Do not use a
   single final displacement bar as the constitutive argument.
6. **Phase-only diagnostic.** RL vs RE histories and matched `IF`/`Rsigma`/
   displacement fields.  Caption must state “one-way numerical
   counterfactual”.

Keep the main text to three composite figures if space is tight: combine the
case/trigger chain, the physical/phase controls, and the constitutive/MPM
response; place mesh refinement and material preflight in the supplement.

## Interpretation rules

- LS–HS shows the physical effect of changing saturation, including its
  inseparable effects on amplitude, drainage and phase.
- RL–RE isolates the timing contribution while holding the fitted pressure
  mean and fundamental amplitude fixed.
- RL–RM is the primary constitutive ablation under an identical prescribed
  pressure history.  HS–HM is retained only to show fully coupled context.  MC
  has no cyclic fabric/state memory, whereas SANISAND
  can represent contraction, pore-pressure accumulation and cyclic mobility.
  MC friction is matched to `Mc`, and its constant elastic tangent is locally
  matched to SANISAND at `p'=3 kPa`, `n=0.485`; state clearly that SANISAND
  remains pressure dependent.  Phrase conclusions in terms of observed paths,
  not an assumed winner.
- Very low permeability is not automatically more dangerous.  The manuscript
  already indicates intense but shallow localization at the lowest `kappa`;
  the chosen intermediate value allows the affected zone to reach the pipe.
- Use `kappa [m2]` for intrinsic permeability and `Kh [m/s]` for hydraulic
  conductivity throughout.
- Do not use `ru` to define liquefaction in the case study.  Use the
  hydrostatic-corrected upward force index `IF` for the hydraulic trigger and
  the signed-ratio-safe actual skeleton stress ratio `Rsigma` for realised
  stress loss.  The legacy code field named `liquefaction_potentials` is a
  vertical-stress-loss measure, not a seepage-force index.

## Manuscript consistency fixes required before submission

1. The paragraph calling `5e-14 m2` “high permeability” should read the actual
   high-permeability value (`5e-10 m2` in the plotted sweep).
2. Reconcile the statement of monotonic liquefaction depth with the plotted
   intermediate-permeability maximum.
3. Do not call the `LI=0.5` plotting contour `ru=0.5`; the stated critical
   seepage-force criterion is `LI=1`.
4. Make the Fig. 6 panel order in the text and caption identical.
5. Reconcile the low-permeability amplitude-ratio values reported for Figs. 1
   and 3.
6. State explicitly whether validation plots use water excess pressure and
   later plots use saturation-weighted mixture pressure.
7. Correct the phase-lag validation citation to the verified Li and Gao paper:
   *Ocean Engineering* 253 (2022), 111315,
   <https://doi.org/10.1016/j.oceaneng.2022.111315>.

## Stop conditions

Do not submit production results if any of the following occurs:

- equilibrium checkpoint is absent or dynamic resume is not explicitly logged;
- SANISAND material/checkpoint tests fail;
- measured LS lag is not low or HS lag is not appreciably higher;
- RL and RE differ materially in fitted means/amplitudes;
- phase database hashes/metadata are absent;
- particles leave the domain or the result contains NaN/Inf;
- the analysis uses `ru` as a liquefaction criterion, fails to remove the
  hydrostatic component from the seepage force, or presents the replay as a
  conservative fully coupled solution.

## Manuscript evidence gate

The completed full-screen analysis must generate `manuscript_evidence.json`
and `manuscript_evidence.md` before prose placeholders are populated.  A phase-
lag liquefaction claim requires, in order: RL/RE pressure mean/amplitude QA,
material reduction of the fitted subsurface lag, greater RL support-zone
`IF>=1` area--time, and (for actual rather than potential liquefaction) greater
RL support-zone `Rsigma<=0.05` area--time.  The SANISAND mechanistic claim uses
RL/RM under the identical lagged pressure history and requires cyclic state
evolution, an observable stress/engineering consequence and a passed first-
frame `p'-q` initial-state audit.  If the audit fails because MC_EQ changed the
initial path materially, report a model-chain comparison and do not attribute
the difference uniquely to constitutive cyclic memory.
