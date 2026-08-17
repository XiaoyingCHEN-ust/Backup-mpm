# Revision notes for `Manuscript-Discussion-0811.pdf`

These notes were prepared from the 17-page PDF created on 11 August 2026.  The
PDF is a user working file and is not committed.  Line references below follow
the visible line numbering in that PDF.

## Required scientific corrections

1. **Validation reference (pp. 1--5, lines 5, 34, 60).** Keep Li and Gao
   (2022), but give the verified record: Chang-Fei Li and Fu-Ping Gao,
   “Characterization of spatio-temporal distributions of wave-induced pore
   pressure in a non-cohesive seabed: Amplitude-attenuation and phase-lag,”
   *Ocean Engineering* 253 (2022), 111315,
   <https://doi.org/10.1016/j.oceaneng.2022.111315>.  Do not replace it with the
   PDF comment's Li--Jeng--Zhang/*Applied Ocean Research* entry; the stated DOI
   `10.1016/j.apor.2022.103143` has no Crossref work record and does not match
   the benchmark title.

2. **Boundary-loading reference (p. 2, line 32).** Replace `[add reference]`
   with the linear-wave-theory source used elsewhere in the methods.  Keep the
   distinction that the surface pressure is prescribed and the overlying water
   column is not discretised.

3. **Pressure definitions (pp. 1--5 versus pp. 6--17).** State explicitly:
   validation uses water-phase excess pressure `p'_w=p_w-p_w0` because the
   transducers measure water pressure; the parametric and pipeline phase metric
   uses the saturation-weighted mixture pressure produced by the three-phase
   model.  These quantities must not share an undefined `p'` symbol in captions.

4. **Permeability notation and units.** The table reports hydraulic
   conductivity in m/s, whereas Sections 7.1--7.4 plot intrinsic permeability
   in m2.  Use `K_h [m/s]` for hydraulic conductivity and `kappa [m2]` for
   intrinsic permeability throughout.  Do not label intrinsic permeability as
   `k_s` if `k_s` denotes conductivity in Table X.

5. **Fig. 1 versus Fig. 3 amplitude values (pp. 6 and 9).** Lines 82--83 give
   the low-permeability amplitude ratio as about 0.07, while lines 120--121 give
   about 0.35.  Recompute both statements from the same harmonic-fit dataset and
   quote one value.  Do not resolve the discrepancy by reading pixels from the
   figure.

6. **Liquefaction terminology in Section 7.2 (pp. 11--12).** The plotted field
   is a seepage-force index and should be called hydraulic liquefaction
   *potential* or trigger.  `IF >= 1` is the critical force balance.  The red
   `IF >= 0.5` contour may be retained only as a geometry visualisation contour;
   it must not be called the critical/liquefied zone.  Actual liquefaction in
   Section 7.5 is evaluated independently from `R_sigma <= 0.05`.

7. **High-permeability typo (p. 14, line 204).** Change `5 x 10^-14 m2` to
   `5 x 10^-10 m2`.  The plotted high-permeability row in Fig. 6 is
   `5 x 10^-10 m2`.

8. **Non-monotonic depth (pp. 13--17).** Lines 180--191 and 208--216 correctly
   state that maximum penetration occurs at intermediate permeability.  Delete
   or rewrite lines 233--241, which claim a systematic monotonic increase with
   permeability.  Fig. 6b itself shows the intermediate-permeability maximum
   (up to about 15.36 cm), not a monotonic trend.

9. **Wrong criterion name (p. 16, line 233).** Replace `r_u=0.5 boundary` with
   `IF=0.5 visualisation contour`.  Better: explicitly distinguish the geometry
   contour (`IF >= 0.5`) from the critical hydraulic area (`IF >= 1`).  Never use
   `ru` to label this field.

10. **Fig. 6 panel order (pp. 16--17).** The caption order is (a) width,
    (b) depth, (c) angle, (d) maximum index.  The prose currently starts with
    “Figure 6a shows maximum liquefaction depth,” which conflicts with the
    figure.  Rewrite the prose in caption order or relabel the panels.

11. **Over-broad phase-lag claim (lines 141--147 and 225--231).** Retain the
    important qualification already present in Section 7.4: larger lag tends to
    localise gradients, but severity also requires pressure penetration.  The
    defensible statement is not “more lag always causes deeper liquefaction.”
    It is “at matched/transmissible loading, lag can intensify and redirect the
    critical gradient; very low permeability produces intense but shallow
    localisation.”

12. **Section 7.5 title and placeholder (p. 17, lines 250--251).** Replace the
    trench-slope title and `xxxx. e from 0.3 to 0.46` placeholder with the
    shallow-pipeline case in `MANUSCRIPT_SECTION_7_5_DRAFT.md`.  The pipeline
    gives registered stability observables and directly tests the engineering
    consequence of the phase mechanism developed in Sections 7.1--7.4.

## Suggested replacement for the contradictory Fig. 6 prose

> Figure 6 synthesises the competition between pressure penetration and
> gradient localisation. At high intrinsic permeability, the response remains
> nearly synchronous and the maximum upward seepage-force index is generally
> small. Reducing permeability into the intermediate range preserves sufficient
> pressure penetration while increasing phase delay and spatial gradients; the
> `IF >= 0.5` visualisation contour consequently reaches its greatest depth in
> this range. At the lowest permeability, the maximum local index remains high,
> but the affected layer becomes shallow and nearly horizontal because the wave
> pressure no longer penetrates to depth. Saturation shifts these regimes by
> changing fluid storage and mobility, but it does not remove their
> non-monotonic permeability dependence.

## Wording hierarchy for the completed paper

- `IF >= 1`: critical upward hydraulic trigger.
- `R_sigma <= 0.05`: realised near-total vertical skeleton-stress loss.
- `ru`: probe-level diagnostic only.
- `IF >= 0.5`: optional visual contour for comparing geometry; not the critical
  threshold and not actual liquefaction.
- RL--RE: phase-only one-way numerical counterfactual.
- RL--RM: matched-pressure constitutive ablation.
- RL--RM must also pass the first-dynamic-frame `p'--q` initial-state audit;
  otherwise it is a comparison of complete handoff/model chains rather than a
  clean constitutive isolation.
- SANISAND advantage: ability to represent cyclic state/fabric evolution and
  mobility, demonstrated only if those mechanisms produce a resolved path and
  observable response consequence.  Do not claim universal or independently
  validated predictive superiority.

## Completed Section 7.5 result (screen analysis, 2026-08-16)

The registered seven-case screen and its independent artifact audit are
complete. Replace the Section 7.5 placeholder with
`MANUSCRIPT_SECTION_7_5_DRAFT.md` and use the generated common-scale 2-D figures
under `analysis/screen/figures/`.

The result supports a scoped local statement, not the broad statement that
phase lag always produces larger engineering motion. The database
transformation preserved fitted means
(maximum change 0 Pa) and amplitudes (maximum change
`5.68 x 10^-14 Pa`) and reduced the mean crown/shoulder/invert lag by 45.37
degrees. The PIC-smoothed response-amplitude difference reached 54.75%, so the
registered response-level control failed and the VTK-smoothed seepage-force
field must not be used as the primary gradient metric. A direct
nearest-neighbour derivative of the raw prescribed-pressure database instead
gives RL a 5.97% larger positive-force impulse, a 1.53% larger raw `IF>=1`
area--time and a 49.07% larger same-particle raw `IF>=1`/`Rsigma<=0.05`
area--time. RL also had a 17.23% larger separate stress-loss integral, while
pipeline motion remained materially indistinguishable. Report this as stronger
local hydraulic triggering and stress-loss coincidence without claiming a
large pipe-motion consequence or universal causality.

A continuous engineering bridge is now available in
`figure8_pipeline_gradient_force_bridge.{png,pdf,csv,audit.json}`. At the
registered `t_IF=7.865 s`, the retained-lag RL case produced a raw support-zone
positive upward pressure-gradient force of 287.75 N/m, versus 265.00 N/m for RE
(+8.59%); its maximum raw `IF` was 22.02% larger and its instantaneous raw
`IF>=1` area was 12.91% larger. RL also had an 8.56% larger full-record peak
force and a 5.97% larger positive-force impulse. Use this to support the scoped
statement that phase lag can concentrate a larger local lifting force and
increase local trigger/stress-loss coincidence, while retaining the null result
for pipeline motion. All metrics use raw database pressures and current
particle volumes. Gaussian smoothing is permitted only for the displayed
two-dimensional raster and must be identified as display-only in the caption.
The opposite VTK-smoothed `IF` ordering must be reported explicitly as a
derivative/smoothing sensitivity.

The completed exploratory `k=1.0 x 10^-13 m2`, `Sw=0.94` run should be used as
a parameter-selection caution, not promoted to the main engineering result. It
used no solver pressure smoothing and resolved a 50.20-degree crown lag, but
the fixed-pipeline lagged/phase-erased positive upward resultants at 3.705 s
were 847.93/849.24 N/m and their impulses were 3086.22/3085.90 N s/m. The low
permeability retained a clearly resolved delay while attenuating useful
penetration, so it did not increase the integrated pipeline-zone uplift action.
If a follow-up is run, pre-register a `Sw=0.94` permeability series around
`10^-13--10^-12 m2`, keep solver smoothing off, and use smoothing only in the
rendered raster.

The SANISAND material-point comparison supports only a model-capability claim:
under the identical cyclic shear history it resolved pressure loss,
contraction--dilation, backstress/fabric evolution and cyclic mobility that the
perfect-plastic Mohr--Coulomb baseline did not. The zero-cohesion MC_EQ field
handoff failed the unchanged `max|v|<=10^-3 m/s` gate, so the paper must not
claim matched field-scale or predictive superiority.
