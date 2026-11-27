# RANS counter-rotating vortex pair, wall-normal coolant transport and plume width: check against the stored literature

Status: search complete. Items a, b and e are described in the stored papers for the direction "RANS plume too narrow, centreline too high". Item a (CRVP too strong) is described by four papers for eddy-viscosity models, and contradicted by two for the k-epsilon model. Item c (mean-flow versus turbulent split) is described pointwise in two papers, with no RANS-versus-reference comparison. Item d is described as a model-sensitivity correlation, not as compensation.

## Scope and method

- Store: `/home/sk999/Documents/4dns_filmcooled/write_up/literature/papers/` (the path in the request, `write-up`, does not exist; the folder is `write_up`). 128 PDFs: background, methods, conference_versions, arxiv, reference.
- Extraction: `pdftotext` (default and `-layout`) for 127 PDFs. Keyword hits were listed per paper and each hit was read in a window of 250 to 500 characters, then the surrounding page where the claim mattered. Keywords: counter-rotating, CRVP, kidney, vortex pair, lift-off, penetration, lateral/spanwise spreading, diffusion, mixing, centerline effectiveness, turbulent heat/scalar flux, advective/advection, mean convection, convective flux, budget, control volume, entrainment, 7-7-7.
- `Acharya et al. 2001.pdf` is image-only. Pages 1-2, 4-7 and 11-13 were read as page images. Its statements are cited from those pages.
- Page numbers are PDF page indices, as in the earlier note. Where the printed page differs, it is given in brackets.
- Labels: **Stated** is in the paper. **Inference** is mine. Quotations are verbatim, at most two sentences.
- `Raiola & Kriegseis 2026` matched "advective" only as a wavepacket term and is not relevant.

## a. RANS over-predicting CRVP strength or coherence, or jet penetration, against LES/DNS/experiment

### Stated: CRVP or secondary vortex too strong in eddy-viscosity RANS

1. Laroche et al. 2023, "Towards an Improved Description of Film-Cooling Heat Transfer Through Anisotropic RANS Modeling" (`background/Laroche et al. 2023.pdf`), p. 11, conclusions. Penn State 7-7-7 hole, k-omega SST against LES, EBRSM and DRSM.
   > "This improvement is due to a better representation of the kidney-shaped vortices typical of a jet in crossflow interaction, and whose magnitude is usually overestimated by eddy viscosity models."
   The sentence is general. No circulation or vorticity number is given for SST. On p. 10 the paper adds: "Paysant [12] showed that the intensity of the vortices is reduced using such a model, which favors better mixing of the jet with the mainstream."
2. Azzi and Lakehal 2002, "Perspectives in Modeling Film Cooling of Turbine Blades by Transcending Conventional Two-Equation Turbulence Models", J. Turbomach. 124 (`background/Azzi & Lakehal 2002.pdf`), p. 12, conclusions. Isotropic two-layer k-epsilon, flat plate and channel.
   > "The isotropic two-layer k-epsilon model was found to underpredict the lateral spreading of the temperature field, together with an exaggerated size of the kidney vortices."
   The same page reports that non-linear algebraic models "exaggerate the size of the jet and its lifting-off from the wall". The size, not the circulation, is exaggerated.
3. Lakehal 2002, "Near-Wall Modeling of Turbulent Convective Heat Transport in Film Cooling of Turbine Blades With the Aid of Direct Numerical Simulation Data", J. Turbomach. 124 (`background/Lakehal 2002.pdf`), p. 1 (abstract) and p. 7.
   > "only with this anisotropy eddy-viscosity/diffusivity model can the spanwise spreading of the temperature field be well predicted and the strength of the secondary vortices reduced."
   p. 11: "Applying the anisotropic eddy-viscosity model tends to reduce the strength of these symmetrical vortices by flattening them towards the surface". Both Lakehal papers show the reduction in contour plots of secondary velocity. Neither gives a circulation.
4. Stratton and Shih 2019, "Evaluating Modeling Weaknesses in Predicting Film Cooling via Large-Eddy Simulations", J. Propul. Power (`background/Stratton & Shih 2019.pdf`), p. 6; conference twin 2018 (`conference_versions/Stratton & Shih 2018.pdf`), p. 6. LES against realizable k-epsilon and SST, cylindrical holes, BR 0.5 and 1.0.
   > "SST overestimates the strength, whereas k-epsilon tends to underestimate the strength."
   The strength is judged from contours of streamwise vorticity. No circulation ratio is given.

### Stated: the opposite sign

5. Coletti et al. 2013, "Turbulent transport in an inclined jet in crossflow", Int. J. Heat Fluid Flow (`background/Coletti et al. 2013.pdf`), p. 1 (abstract). MRV experiment against standard k-epsilon with Sc_t 0.85.
   > "The computations underestimate the strength of the counter-rotating vortex pair, due to an overestimated eddy viscosity."
   Ling et al. 2017 (`background/Ling et al. 2017.pdf`), p. 1, repeats this: RANS "underpredicted the strength of the counter-rotating vortex pair and significantly overpredicted the centerline film cooling effectiveness."
6. Stratton and Shih 2019, p. 6 and p. 11: k-epsilon "tends to underpredict the CRVP strength".

So the stored literature gives both signs. The sign follows the model: SST over, k-epsilon under (Stratton and Shih), standard k-epsilon under (Coletti). Inference: the 1.2 to 1.8 ratio is consistent with the SST-type result, and with Laroche's general sentence. It is not consistent with the k-epsilon results. The model used for the present RANS decides which precedent applies.

### Stated: jet penetration

7. Acharya et al. 2001, "Flow and Heat Transfer Predictions for Film Cooling", Ann. N.Y. Acad. Sci. 934 (`background/Acharya et al. 2001.pdf`), PDF p. 6 [printed p. 115], from the page image.
   > "The lateral spread and mixing of the film cooling jet is under-predicted by the two-equation models ... while the vertical penetration is over-predicted." The text adds that "the lower lateral spreading has been linked to under-estimation of the eddy viscosity in the lateral direction."
   The statement is a review of earlier work, mostly vertical or square jets. Repeated at second hand by Tyagi and Acharya 2003 (`background/Tyagi & Acharya 2003.pdf`), p. 1; Johnson and Kapat 2013 (`methods/Johnson & Kapat 2013.pdf`), p. 2; Mahesh 2013, Annu. Rev. Fluid Mech. (`background/Mahesh 2013.pdf`), p. 15 [printed p. 393]; Ellis and Xia 2024 (`background/Ellis & Xia 2024.pdf`), p. 2.
8. Contradicting: Foroutan and Yavuzkurt 2015, J. Heat Transfer 137 (`background/Foroutan & Yavuzkurt 2015.pdf`), p. 8, M = 1.5, k-epsilon (realizable) against DES and experiment.
   > "The RANS simulations not only cannot predict enough mixing between the jet and crossflow in spanwise direction, they also fail to capture the jet lift-off phenomenon in this case."
   Here RANS penetrates less, not more. Ellis and Xia 2024, p. 2, report the same for Walters and Leylek: "at higher blowing ratios, these models failed to reproduce the experimentally observed jet lift-off phenomenon." At M = 1.5 the CRVP is "considerably stronger" in DES than at low M (p. 6), and no RANS CRVP is compared.

### Not found

- No stored paper gives a RANS-to-reference circulation ratio. The 1.2 to 1.8 figure has no precedent number.
- No stored paper states that RANS over-predicts CRVP coherence (persistence downstream) for a shaped hole. Coletti (p. 11) states the reverse for persistence: a lower eddy viscosity "to capture the persistence of the CVP".

## b. Under-predicted lateral spreading, over-predicted centreline effectiveness, and the stated cause

### Stated: the effect

- Foroutan and Yavuzkurt 2015, p. 1: "RANS simulations underpredict the lateral spreading and mixing of the coolant jet". p. 6: "The RANS simulations overpredict the centerline effectiveness in the near-field region, which is directly related to the underprediction of mixing between the jet and crossflow in this region." p. 6 also records that RANS and DES converge for x/D > 15.
- Laroche et al. 2023, p. 3: "the lateral spreading of the film was shown to be largely underestimated by RANS models, whereas the centerline effectiveness was usually overestimated." p. 10 (7-7-7): "the centerline effectiveness overestimates the experimental value by approximately 25%, after X/D = 10", read there as "a containment of the vortices in the central region after injection".
- Haydt et al. 2017, J. Turbomach. 139 (`background/Haydt et al. 2017.pdf`), p. 8, realizable k-epsilon, 7-7-7 with offsets: "the centerline effectiveness is overpredicted but lateral spreading is underpredicted".
- Harrison and Bogard 2008 (`background/Harrison & Bogard 2008a.pdf`), p. 2: "simulations with standard two equation turbulence models generally result in less lateral spreading than is seen in experiments". p. 3 notes that over-predicted centreline values "are generally compensated for to some degree by lower values in between the holes".
- Azzi and Lakehal 2002, p. 10: when centreline effectiveness is over-predicted and the laterally averaged value is under-predicted, "this corresponds to a lateral spreading that is predicted as being too weak".
- Ling et al. 2016, J. Turbomach. 138 (`background/Ling et al. 2016a.pdf`), p. 6: with frozen LES mean fields, "all models overpredicted the centerline adiabatic effectiveness downstream of injection".

### Stated: the cause

| Cause given | Paper, page | Quotation |
|---|---|---|
| Isotropic eddy viscosity, anisotropic turbulence | Foroutan and Yavuzkurt 2015, p. 1 | "they fail to correctly describe the anisotropic turbulence associated with the jet-crossflow interactions" |
| Eddy-viscosity assumption | Laroche 2023, p. 3 | "The explanation for such an incorrect behavior relies in the eddy viscosity assumption mostly used for RANS models." |
| Under-estimated lateral eddy viscosity (u'w') | Acharya 2001, PDF p. 6 and p. 12 [printed p. 121] | "the lateral eddy viscosity (xz component) near the surface is significantly greater than the vertical eddy viscosity (xy component)" (LES, page image) |
| Under-estimated spanwise turbulent heat diffusion | Laroche 2023, p. 8 | "low values coincide mainly with an underestimated turbulent heat diffusion in the spanwise direction" (a priori GGDH test on 7-7-7 LES) |
| Near-wall eddy viscosity near zero; one scalar nu_t cannot give u'w' and u'v' together | Stratton and Shih 2019, p. 11 | "RANS also underpredicts the near-wall temperature gradient, this leads to the noticeable lack of lateral spreading of the coolant." |
| Constant Pr_t too high | Stratton and Shih 2019, p. 10 | "In the body of the jet, Pr_t increases to ~0.5, but remains well below the value of 0.85 imposed by RANS." (near-wall edge of jet about 0.1) |
| Constant Sc_t inadequate | Coletti 2013, p. 11 | "Clearly the constant ScT model is of limited value, even in a relatively simple configuration as the present one." |
| Counter-gradient diffusion near the wall | Stratton and Shih 2019, p. 9; Oliver et al. 2017, p. 11-12 | GDH "found to be incorrect" for the 7-7-7 hole (Oliver 2017, p. 12) |
| Pr_t variation moves centreline and span together | Qaragoez and Mazaheri 2025 (`background/Qaragoez & Mazaheri 2025.pdf`), p. 2, review of an earlier study | "A reduction in Pr_t increased effectiveness across the span for high blowing ratios, while for low blowing ratios, it decreased centerline effectiveness but improved lateral cooling." |

### Stated: corrections that restore spreading

Anisotropic eddy viscosity (Bergeles correction; Lakehal 2002; Azzi and Lakehal 2002; Yu and Yavuzkurt 2018, `background/Yu & Yavuzkurt 2018.pdf`, p. 1, a factor 4 in z on DES); EBRSM or DRSM with GGDH (Laroche 2023, p. 9-10); HOGGDH (Ling 2016, p. 6; Ellis and Xia 2024, p. 8).

### Stated: counter-evidence on cause

- Harrison and Bogard 2008, abstract: "Use of the anisotropic RSM model was not found to predict more realistic coolant spreading than the other models." Fawcett 2011 (`background/Fawcett 2011.pdf`), p. 31, reads this as contradicting the anisotropy hypothesis, and on p. 47 points to unresolved coherent unsteadiness.
- Acharya 2001, PDF p. 7 [printed p. 116]: "despite the better representation of turbulence anisotropy, the RSTM predictions did not provide the expected improvements".
- Stratton and Shih 2019, p. 6: "Interestingly, the lower CRVP strength does not translate to higher spreading."
- Milani et al. 2021, J. Fluid Mech. 906 (`background/Milani et al. 2021.pdf`), p. 20, frozen LES velocity, inclined jets at r = 1: "The isotropic model ... broadly overestimates mixing which transports scalar towards the wall immediately downstream of injection". The jet extent was larger in the span. This is the opposite sign, at fixed mean velocity.

### Inference

The papers that vary the closure at frozen mean velocity (Ling 2016, Milani 2021) find mixed signs. The papers that run full RANS (Laroche, Lakehal, Azzi and Lakehal) change the vortex pair and the spreading together. That fits the present finding that the CRVP and the plume width differ together. It is not stated in any of them as a link.

## c. Budget or control-volume split into mean-flow (advective, vortex-induced) and turbulent transport

### Stated: closest matches

1. Oliver et al. 2017, "Implicit LES for Shaped-Hole Film Cooling Flow", ASME GT2017-63314 (`background/Oliver et al. 2017.pdf`), p. 10-11, section "Downstream mean energy equation fluxes". 7-7-7 hole, BR 2, DR 1.6, iLES. The mean convective and Reynolds heat fluxes are compared on y-z planes at x/d = 1, 3 and 5, in both the wall-normal and the spanwise direction.
   > "The magnitude of this flux is significantly less than the mean convective flux, but still substantial, with the maximum Reynolds heat flux generally being approximately 15% of the maximum mean convective flux in the wall-normal direction and approximately 10% in the spanwise direction."
   Also p. 10: "for a RANS simulation to accurately capture the mean temperature field, it is necessary to get the strength of this CRVP correct." The paper identifies the convective flux pattern with the CRVP. No RANS run is compared with the split. Domain extends about 17.5d past the hole. Planes stop at x/d = 5.
2. Coletti et al. 2013, p. 4 and p. 8, p. 1. Inclined 30 degree cylindrical jet. A streamtube control volume is impermeable to the mean flux, so the loss of concentration flux through it measures turbulent diffusion.
   > "Turbulent diffusion acts to decouple mean velocity and concentration fields, as demonstrated by the drop in concentration flux within the streamtube issued from the hole." (p. 1)
   > "the concentration iso-surfaces are almost oval, indicating that the turbulent transport is dominating over the advection." (p. 4, after the double-lobed CRVP footprint "up to about three hole diameters downstream")
   The split is along the jet streamtube, not wall-normal or lateral planes. It gives a volume-averaged diffusivity, not a flux ratio.
3. Milani et al. 2021, Appendix A, p. 26-27, Fig. 19, inclined jets, r = 1, LES. The turbulent scalar flux is shown divided by the mean advection component, for the streamwise and wall-normal components.
   > "streamwise mean advection is usually thought as balancing turbulent mixing in the transverse direction"; "v'c' is relatively more important than u'c'; its magnitude is particularly high in the shear layer all around the jet and in some places close to the wall and right after injection."
   Milani et al. 2020, ASME GT2019-90498 (`background/Milani et al. 2020.pdf`), p. 3, weights the least-squares diffusivity by the inverse of the mean advection in each direction.
4. Schroeder and Thole 2017a, J. Turbomach. 139 (`background/Schroeder & Thole 2017a.pdf`), p. 5; twin 2016b (`conference_versions/Schroeder & Thole 2016b.pdf`), p. 5. Measured 7-7-7 thermal field and PIV, M = 3.
   > "This agreement for the top edge of the jet again was consistent with how coolant was transported primarily by mean-flow convection."
   p. 4 and p. 8: the jet lifts "due to action of the counterrotating vortex pair" far downstream, at low freestream turbulence. Statement on the top edge of the jet only. At high freestream turbulence (13.2%), dilution and lateral dispersion are attributed to turbulent shear stress (p. 7-8).
5. Foroutan and Yavuzkurt 2015, p. 8, M = 1.5. RANS turbulent heat flux is reconstructed from the gradient model and compared with DES.
   > "The discrepancy between the RANS and DES predictions is relatively small for this component, confirming the fact that RANS simulations are known to be adequate in predicting the mixing in the wall-normal direction."
   The streamwise and spanwise components are higher in DES. This is a statement on turbulent flux only, at the centreline plane, not a mean-versus-turbulent split.
6. Oliver et al. 2019, Int. J. Heat Mass Transfer 140 (`background/Oliver et al. 2019.pdf`), p. 8: the CRVP "is responsible for substantial wall-normal and spanwise transport of heat". p. 13: the turbulent heat flux normal to the surface sits "primarily at the top of the coolant jet", the lateral flux at the sides.

### Inference

- Oliver 2017 gives a ratio of maxima of about 0.15 (wall-normal) for the reference at x/d <= 5. The present finding places the reference turbulent flux above the mean flux downstream of x/D = 5, on a different measure (flux through the top of a near-wall layer). The two are not in conflict at x/d <= 5, but Oliver 2017 does not reach the region where the present claim begins. Coletti's statement (turbulent transport dominating after about 3D) is in the same sense as the present claim, for a cylindrical jet.
- Oliver's maxima are pointwise maxima of the flux magnitude, not integrals over a layer top. A layer-integrated ratio can differ from a ratio of maxima.

### Not found

- No stored paper splits the wall-normal flux through a plane at the top of a near-wall layer into mean and turbulent parts, for LES or RANS.
- No stored paper reports the mean-flow share of the wall-normal flux for a RANS run, or its ratio to the reference.
- No stored paper reports a RANS turbulent wall-normal flux of about half the reference.

## d. RANS compensating missing turbulent mixing by stronger mean-flow transport; stronger CRVP from too little turbulent diffusion of vorticity

### Stated: model sensitivity correlation

- Stratton and Shih 2019, p. 8; 2018, p. 8: "Earlier, it was noted that k-epsilon tends to underpredict the strength of the counter-rotating vortex pair, whereas SST tends to overpredict the strength of the counter-rotating vortex pair. Thus, these parameters are directly correlated." The parameters are the sign of the eddy-viscosity error (k-epsilon over, SST under). p. 11: "Matching the eddy viscosity predicted by LES could yield more accurate predictions of the CRVP strength in RANS simulations."
- Coletti et al. 2013, p. 11: "to capture the persistence of the CVP, a lower (or streamwise decreasing) eddy viscosity would be more appropriate", while entrainment needs "a higher (or streamwise increasing) eddy diffusivity". The abstract attributes the weak CRVP to "an overestimated eddy viscosity".
- Lakehal 2002, p. 7 and Azzi and Lakehal 2002, p. 9: raising the lateral turbulent stress reduces the vortex strength and spreads the jet. This is the direction "more turbulent diffusion, weaker vortex pair".
- Laroche 2023, p. 10-11: the EBRSM "reduces" the vortex intensity, "which favors better mixing of the jet with the mainstream", and cuts the centreline error by nearly half.
- Oliver 2017, p. 10: a RANS must get both the CRVP strength and the Reynolds heat flux right.

### Stated: a vorticity-equation argument, other flow

- Sabnis et al. 2021, "Analysis and extension of the quadratic constitutive relation for RANS methods", Aeronaut. J. (`background/Sabnis et al. 2021.pdf`), p. 14. Corner flow, not a jet. The mean streamwise vorticity equation contains viscous diffusion of vorticity and terms "related to turbulent stress anisotropy". Linear eddy-viscosity models "do not capture these vortices" (p. 3), which is the reverse of an over-strong pair.

### Stated: compensation, different sense

- Azzi and Lakehal 2002, p. 10: the laterally averaged effectiveness "may, in fact, be misleading and must therefore be seen as the consequence of a compensation effect of eta being predicted too high near the axis, but not spreading sufficiently in the lateral direction". Lakehal 2002, p. 10: "the fast decay of eta in the jet axis can be compensated by the enhanced lateral spreading." Harrison and Bogard 2008, p. 3, the same for centreline against between-hole values.
  Both concern compensation inside the effectiveness field. Neither concerns compensation between mean-flow and turbulent transport.

### Inference

The stored papers support a model-dependent link between low eddy viscosity and a strong CRVP (Stratton and Shih, Coletti), and between more lateral turbulent transport and a weaker CRVP (Lakehal, Azzi and Lakehal, Laroche). They do not state a vorticity-budget mechanism, and none uses the word compensation for a trade between mean and turbulent transport. The wall-normal reading in the present finding (stronger upwash carries coolant up while turbulent flux is low) combines the CRVP result with Oliver's observation that the convective flux is set by the CRVP. That combination is mine.

## e. 7-7-7 shaped hole

### Stated

- Laroche et al. 2023: the only stored paper that runs RANS (SST, EBRSM, DRSM) and LES on the 7-7-7 hole. Quotations in a and b: SST centreline about 25% too high after X/D = 10; vortex magnitude "usually overestimated by eddy viscosity models" (general); spanwise turbulent heat diffusion under-estimated (p. 8). One 7-7-7 result runs against a CRVP-too-strong reading: at X/D = 20, p. 9, "the reference LES simulation exhibits an outer vortex which is intense, contrary to the SST simulation". p. 5: the inner pair weakens relative to the secondary vortices on the jet sides between X/D = 5 and 10.
- Haydt et al. 2017, p. 8: realizable k-epsilon on the 7-7-7 hole, centreline over-predicted, lateral spreading under-predicted. Haydt et al. 2018, J. Turbomach. 140 (`background/Haydt et al. 2018.pdf`), p. 10: "the CFD is underpredicting the lateral diffusion of the coolant jet" and "underpredicting the mixing that occurs between the coolant jet and the mainstream flow". "These mispredicitons are well-known, and consistent with previous comparisons made for this facility and for the 7-7-7 hole in particular." Haydt and Lynch 2018 (`conference_versions/Haydt & Lynch 2018.pdf`), p. 12-13: steady RANS gets trends, not levels.
- Oliver et al. 2017 and 2019: the 7-7-7 mean energy budget in c. Both are LES, not RANS. Oliver 2017, p. 2: "common engineering RANS turbulence models, like the k-epsilon model, often dramatically overpredict the film cooling adiabatic effectiveness on the hole centerline".
- Schroeder and Thole 2017a (7-7-7 experiment): CRVP lifts the M = 3 jet far downstream; delayed lift-off to about x/D = 14 at M = 3; "primarily by mean-flow convection" at the jet top (c).
- Haydt and Lynch 2018, p. 3 (literature statement): a shaped diffuser "decreases jet penetration into the mainstream, increases lateral spread of coolant on the endwall, and decreases the strength of the CRVP".
- Bin and Yang 2024 (`background/Bin & Yang 2024.pdf`), p. 2: DNS of the 7-7-7 hole over compound angles; it records that earlier CFD of the hole "are limited to large-eddy simulations (LES) and Reynolds-averaged Navier Stokes (RANS)". No RANS comparison. Nunno et al. 2022 (`background/Nunno et al. 2022.pdf`): wall-resolved LES, 7-7-7, hole angle and BR.
- Fan-shaped: Laroche 2023, p. 3, cites Rouina (DRSM better on adiabatic effectiveness maps, BR 1 and 2). Zamiri et al. 2020 (`background/Zamiri et al. 2020.pdf`), p. 2: for a laidback fan-shaped hole "the lateral data exhibited differences as compared to the experimental results" with k-epsilon.

### Not found

- No stored paper reports CRVP circulation, or mean versus turbulent wall-normal flux, from a RANS of the 7-7-7 hole.
- No stored paper compares a 7-7-7 RANS with a 7-7-7 DNS. The only DNS (Bin and Yang 2024) has no RANS comparison.

## What has not been found described

1. A RANS-to-reference circulation ratio for the CRVP, for any hole. The ratio 1.2 to 1.8 has no precedent number in the store. Precedents are contour-based and of mixed sign.
2. A mean-flow versus turbulent split of the wall-normal coolant flux through the top of a near-wall layer, for RANS or for LES. The nearest are Oliver 2017 (ratio of maxima, x/d <= 5, 7-7-7 LES only) and Coletti 2013 (streamtube, cylindrical jet).
3. A RANS turbulent wall-normal flux of about half the reference. Foroutan and Yavuzkurt 2015 report a small RANS-DES difference in that component.
4. A statement that the downstream turbulent flux dominates beyond x/D about 5 in the reference. Coletti's "turbulent transport is dominating over the advection" after about 3D is the nearest, for a different configuration.
5. A statement that a RANS replaces missing turbulent transport by stronger secondary-flow transport. The stored papers describe compensation only inside the effectiveness field (centreline against span, Azzi and Lakehal; Harrison and Bogard).
6. A vorticity-budget argument (too little turbulent diffusion of vorticity gives a stronger CRVP) for a film-cooling jet. Sabnis 2021 gives the equation for corner flow with the opposite outcome.
7. RANS plume width quantified as a ratio (0.75 at x/D 20). Papers give the direction and, for Foroutan and Yavuzkurt 2015, footprint lengths in diameters.
8. Any 7-7-7 RANS result stated as a CRVP strength. The 7-7-7 results concern effectiveness only.
9. Not searchable: Acharya et al. 2001 was read as page images for the cited pages only (1-2, 4-7, 11-13). Schmitt 2007 is not stored.
