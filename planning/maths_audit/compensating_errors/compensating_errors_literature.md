# Compensating errors in eddy viscosity: check against the stored literature

Status: search complete. No stored paper describes the specific mechanism (RANS length scale too short, offset by assumed perfect stress-strain alignment, with equal k). Closest precedents are listed below.

## Scope and method

- Searched: 128 PDFs in `/home/sk999/Documents/4dns_filmcooled/write_up/literature/papers/` (background 68, methods 45, conference_versions 6, arxiv 7, reference 2; plus `Closure-error literature review ... .md` and `papers/index.md`).
- Extraction: `pdftotext -layout` (pypdf is not visible under `python3 -I`; it sits in the user site-packages). Text in `/tmp/claude-1001/-home-sk999-Documents-ASME-2027/f4c92ceb-6df8-4254-96c3-6552d15d0ecf/scratchpad/lit_text/`. 127 of 128 PDFs gave text; `Acharya et al. 2001.pdf` is image-only (16 bytes of text) and was not searched.
- Terms grepped: compensat, cancel, fortuit, wrong reasons, coincidental, offsetting, misalign, alignment, principal axes, Schmitt, validity, optimal / effective eddy viscosity, negative eddy viscosity, ill-conditioned, small strain, length scale, k over/under-prediction, lateral/wall-normal mixing. Every hit listed below was read in context. Page numbers are PDF page indices.
- Already cited in `references.bib`: stratton2018, stratton2019, sandberg2018, sandberg2018b, lav2019, laroche2023, ling2016 (= Ling et al. 2017, J. Turbomach. 139(2)), coletti2013, milani*, schmitt2007 (Schmitt 2007 is cited but its PDF is not stored).

## Result by category

### (a) Compensating or cancelling errors

No passage states that a RANS eddy viscosity, diffusivity or effectiveness agrees through cancelling errors of its factors. The related passages are:

| Paper (file) | Page | Quote | How it relates |
|---|---|---|---|
| Foroutan and Yavuzkurt 2015, J. Turbomach. (`background/Foroutan & Yavuzkurt 2015.pdf`) | 11 | "RANS simulations underpredict jet and crossflow mixing both in vertical and spanwise directions. This results in overprediction of centerline effectiveness in the near-field region." | Both mixing directions are too low and act in the same sense. No cancellation between directions is described. |
| Laroche et al. 2023, ASME Turbo Expo (`background/Laroche et al. 2023.pdf`) | 3 | "the lateral spreading of the film was shown to be largely underestimated by RANS models, whereas the centerline effectiveness was usually overestimated. The explanation for such an incorrect behavior relies in the eddy viscosity assumption" | Same: one-sided error, attributed to the eddy-viscosity assumption. |
| Yu and Yavuzkurt 2018 (`background/Yu & Yavuzkurt 2018.pdf`) | 1 | "To compensate for this shortcoming of DES, the eddy viscosity in the spanwise direction is increased to enhance spanwise-diffusion" | A deliberate correction (anisotropic factor added to nu_t), not a fortuitous cancellation. |
| Sandberg et al. 2018, J. Turbomach. 140(10) (`methods/Sandberg et al. 2018a.pdf`; conference twin `conference_versions/Sandberg et al. 2018b.pdf`) | 11 | "the heat flux model developed for the thick-lip case is overcompensating for the lack of mixing due to the deterministic unsteadiness" | A trained model corrects one error and over-corrects a missing mechanism in another case. Related to compensation, but between cases, not between factors. |
| Ovchinnikov et al. 2008, J. Fluid Mech. 613 (`background/Ovchinnikov et al. 2008.pdf`) | 14 | "this agreement is probably coincidental" | Bypass transition DNS; agreement attributed to a different transition stage. Same wording class, different physics. Not relevant to the finding. |
| Menter et al. 2006 (`background/Menter et al. 2006.pdf`) | 1 | "predict transition seems coincidental" | Transition modelling. Not relevant. |

### (b) Stress-strain misalignment and Boussinesq validity

| Paper (file) | Page | Quote | How it relates |
|---|---|---|---|
| Stratton and Shih 2019, J. Propul. Power (`background/Stratton & Shih 2019.pdf`) | 1, 8-9 | "Although regions of counter-gradient diffusion and stress-strain misalignment were found in the near-wall region, further above the wall, the jet behaved according to the hypotheses." Also: "By using the least-squares approach, it is possible to solve for a scalar eddy viscosity based on the Boussinesq hypothesis" (Eq. 19) | Closest precedent. Defines the least-squares optimal eddy viscosity from LES, identical in form to mu_eff = -a:S/(2\|S\|^2), and finds misalignment near the wall. Does not factor nu_t into k, L and cos(theta). |
| Stratton and Shih 2018, AIAA 2018-4734 (`conference_versions/Stratton & Shih 2018.pdf`) | 13 | "There is some evidence of misalignment of the stress and strain at higher VR ... Using the method of least-squares an 'ideal' scalar eddy-viscosity can be solved for using the LES data." | Same, conference version. |
| Sandberg et al. 2018 (`methods/Sandberg et al. 2018a.pdf`) | 5 | "the overprediction of the adiabatic wall effectiveness can be attributed to the inability of the linear RANS model to reproduce the high level of near wall anisotropy ... quantified by considering the alignment of the strain rate tensor and the anisotropy tensor" (Eq. 15, c in [-1, 1]); "a region of poor alignment occurs close to the wall, while strong alignment occurs in the wake region" | Defines a normalised stress-strain alignment from LES that is the cosine of the paper's theta (c = a_ij S_ij / (\|a\| \|S\|)). Reports poor alignment near the wall. Does not relate it to eddy-viscosity magnitude. |
| Lav et al. 2019, J. Comput. Phys. (`methods/Lav et al. 2019.pdf`) | 4-5 | "the validity of the Boussinesq hypothesis was, for example, tested by Schmitt [28] through computing the alignment of the modelled anisotropy tensor with the reference one ... A completely valid model would exhibit perfect alignment with a value of 1.0" | Definition of rho_RS and the Schmitt (2007) reference. Schmitt 2007 itself is not stored. |
| Lav et al. 2019 | 16 | "the alignment values of 0.7-0.8 in the regions x = 5.0 - 15.0 can affect the solution" | Quantifies alignment of the Boussinesq approximation in an LES-referenced case (not film cooling). |
| Laroche et al. 2023 | 3 | "Stratton and Shih [14] showed in particular that the turbulent stress and strain fields are aligned far from the wall, but their misalignment cannot be neglected near the wall, leading to large discrepancies." | Secondary statement of the Stratton and Shih result. |
| Gorle and Iaccarino 2013 (`methods/Gorle & Iaccarino 2013.pdf`) | 7 | "The specific assumption is that there is an alignment between the Reynolds stress tensor and the anisotropy tensor. For simple shear flows this is correct" | Statement that the linear model assumes alignment. |
| Iaccarino et al. 2017 (`methods/Iaccarino et al. 2017.pdf`) | 2 | "The inconsistency in the eigenvector alignments between model predictions and experimental studies has been identified as a major source of discrepancy in the RANS modeling framework" | General statement, no film cooling. |
| Milani, Ling and Eaton 2021, J. Fluid Mech. (`background/Milani et al. 2021.pdf`) | 11 | "regions of type 2 transport ... coincide with severe misalignment of -u_i c and dc/dx_i, and in such areas the angle theta is frequently higher than 90 deg" | Misalignment for the scalar flux vector (not the stress), angle theta defined; same geometry of argument. |
| Ling et al. 2016 / 2017, J. Turbomach. 139(2) (`background/Ling et al. 2017.pdf`) | 3-4 | "there are significant regions of this flow with negative eddy viscosities" and "S11 ... is greater than zero on the lee side of the jet, leading to the calculation of a negative eddy viscosity" | Anisotropy of normal stresses as the cause of the inferred value. Also see (d). |

### (c) Length-scale and k errors of two-equation models

| Paper (file) | Page | Quote | How it relates |
|---|---|---|---|
| Stratton and Shih 2018 (`conference_versions/Stratton & Shih 2018.pdf`) | 7, 13 | "SST fails to capture the magnitude, and tends to predict very low k in the middle of the jet near the wall. RANS can capture the scaling of k with VR, but underpredicts k near the wall." / "SST under predicts k throughout the entire jet region" | Direct statement on k for SST in film cooling: too low, not equal or too high. Differs from the paper's finding (k equal in the band, too high above it). Different case (round hole, 35 deg, M 0.1). |
| Stratton and Shih 2019 | 6, 9 | "SST tends to underpredict the eddy viscosity" / "both models predict very low values of eddy viscosity near the wall" | nu_t error of SST against the LES least-squares value. Not decomposed into k and L. |
| Stratton and Shih 2019 | 8 | "the eddy viscosity from LES approaches zero too far from the wall" (2018 version, p. 13); "eddy viscosity approaching zero further away from the wall as seen in Fig. 10" | Near-wall error of the LES-inferred value; related to the paper's near-wall L ratio of 0.2-0.4. |
| Azzi and Lakehal 2002, J. Heat Transfer (`background/Azzi & Lakehal 2002.pdf`) | 12 | "The isotropic two-layer k-epsilon model was found to underpredict the lateral spreading of the temperature field, together with an exaggerated size of the kidney vortices." | Isotropic eddy viscosity and lateral spreading. No length scale analysis. |

No stored paper compares the k-omega SST length scale (k^(1/2)/omega or the dissipation length) with an LES or DNS length scale in film cooling or a jet in crossflow.

### (d) Undefined, negative or ill-conditioned inferred eddy viscosity

| Paper (file) | Page | Quote | How it relates |
|---|---|---|---|
| Stratton and Shih 2019 | 8 | "Note that the contours are clipped by theta = 0.05, as Eq. (19) is not applicable in the freestream where gradients are very low." (2018 version p. 8: "the contours are clipped, as Eq. 8 is not applicable in the freestream where gradients are very low") | Mask on the inferred eddy viscosity (a temperature threshold, not a strain threshold). Direct precedent for masking. |
| Ling et al. 2017 (J. Turbomach. 139(2)) | 3 | "in real flows, it is possible for momentum to be transported counter to a velocity gradient, leading to negative eddy viscosities. An eddy viscosity field can be extracted from the LES as" (Eq. 2, a least-squares form) | The inferred value can be negative. Negative regions outlined in Fig. 5: near wall downstream of injection and in the shear layer. |
| Ling and Templeton 2015, Phys. Fluids 27 (`methods/Ling & Templeton 2015.pdf`) | 13 | "In the Boussinesq hypothesis, the existence of a non-negative eddy viscosity is presumed. However, in regions of rapidly changing velocity gradients, this assumption may not be valid." | Negative nu_t used as an error marker for the Boussinesq hypothesis. |
| Duraisamy et al. 2019, Annu. Rev. Fluid Mech. (`methods/Duraisamy et al. 2019.pdf`) | 10 | "regions where negative eddy viscosity is possibly emerging in several flows, likely indicating a failure of linear eddy viscosity models" | Same, review. |
| Wu et al. 2019 (J. Fluid Mech. 869) and Amarloo et al. 2023 | - | "ill-conditioned" | Concerns propagation of an explicit Reynolds stress into the mean-flow equations, not inference of nu_t. Not relevant to the masking question. |

No stored paper gives a strain-magnitude threshold for the inferred eddy viscosity.

## Closest precedents, in order

1. Stratton and Shih 2019 (and 2018): least-squares optimal nu_t from LES; SST underpredicts it; stress-strain misalignment found near the wall and not away from it; masking in low-gradient regions. No factorisation, no compensation.
2. Sandberg et al. 2018: LES strain-anisotropy alignment c (equal to cos theta) mapped for a film-cooling slot; poor alignment near the wall; links the effectiveness overprediction to it.
3. Ling et al. 2017: LES-inferred nu_t (Eq. 2), negative regions, cause traced to normal-stress anisotropy.
4. Lav et al. 2019: rho_RS alignment after Schmitt 2007.
5. Foroutan and Yavuzkurt 2015 and Laroche et al. 2023: one-sided RANS deficits (too little vertical and lateral mixing, centreline effectiveness too high). Context for the "too little lateral spreading" part of the question.

## Verdict

The stored literature does not describe an eddy viscosity matching the reference through a short length scale offset by an assumed alignment of cos = 1. The decomposition mu = rho k^(1/2) L cos(theta) with a per-factor error ratio is not found in any stored paper. The alignment term is found only as a diagnostic (Sandberg 2018, Lav 2019, Stratton and Shih), never as a factor that multiplies the RANS value. No stored paper says "compensating errors" for RANS film cooling.

## Neutral statement for the paper

"Stress-strain misalignment near the wall and an eddy viscosity from the SST model that departs from the least-squares value inferred from scale-resolving data have been reported for film cooling [stratton2018, stratton2019, sandberg2018]. In the present case, the ratio of the RANS to the reference eddy viscosity is close to unity in one band, and the factorisation shows that this ratio results from a shorter length scale, a perfect-alignment assumption and an equal turbulence kinetic energy acting together."

## Not stored

- Schmitt 2007, C. R. Mecanique 335, 617-627 (cited by Lav et al. 2019 for the alignment test; key `schmitt2007` is in `references.bib`; PDF not stored; existence not re-verified on the web, as the task allowed this only for a work cited for the mechanism, and Lav cites it for the alignment definition only).
