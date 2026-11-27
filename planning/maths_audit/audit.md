# Maths audit of the effective-eddy-viscosity split (E9, E11-E13), 2026-11-26

Status: audit complete, six questions answered. Read-only on 4dns_filmcooled, Tank and Elements. Scripts and logs are in this folder.
Format note: this file is a text deliverable as requested for the audit; no PDF was built.

Notation: a = anisotropic Reynolds stress (positive sign, a = rho u''u'' - 2/3 rho k I, Favre); S* = deviatoric sym(grad U~); |.| = Frobenius norm;
cos = -a:S*/(|a||S*|); L_eff = k^1/2 |b|/|S*|, |b| = |a|/(2 rho k); L_al := mu_eff/(rho k^1/2).
Region numbers below marked "file" are from `analysis/meff_rans50_c193b10/results_pp/regions_se.txt`; numbers marked "mine" are my reproduction
on the same point set without the delete-one condition (`q234_regions.log`). The two agree within 0.02 in both film regions; the jet differs by 0.12
(see Q4).

## Summary

| Q | Verdict | Key number |
|---|---|---|
| 1 | Holds | max relative residual 2.0e-7 on 5.6 M cached 193M points (float32 storage); 9e-16 in double precision on random tensors |
| 2 | Holds (algebra); the "length carries the film error" reading is partly a property of the definition | length + align = ln(L_RANS/L_al,193): film 3-10 -0.157, film 10-25 -0.355, jet +0.517; shared ln(1/cos) = the alignment term itself (0.806, 0.573, 0.643) |
| 3 | Partly: ε of the 193M is not cached; done on the 50M 3DNS (exact) and on the 193M with a proxy | L_eps split, 50M: film 3-10 length -0.888 / structure +0.695; film 10-25 -1.013 / +0.585; jet -0.131 / +0.684 |
| 4 | Partly holds | film exclusion of cos <= 0 is small (1.5 %, 0.4 % of volume); the tail cos < 0.3 holds 65 % / 42 % / 61 % of the alignment term (film 3-10 / 10-25 / jet). Jet total flips sign with the estimator: mean of ln +0.68, LS-weighted -0.51 |
| 5 | Holds | contraction QCR : S* = 0 to 3e-15; mu_eff(QCR) = mu_t to 2e-15; simple shear |a_perp|/|a_par| = 2 c_cr1 = 0.61, cos = 0.852 |
| 6 | Defects found, none changes a film sign | far-film region is compared over x/D 10-20.4 only (50M buffer), approach layer over -10.7 to -3; film-region cos median 0.62-0.72 against 0.83 on the centre plane; E13 and E10 use different estimators; SE covers sampling only |

## Q1. Is mu_eff = rho k^1/2 L_eff cos exact as coded?

Verdict: holds.

Derivation. With b = a/(2 rho k) and |b| = |a|/(2 rho k):

    rho k^1/2 L_eff cos = rho k^1/2 * k^1/2 |b|/|S*| * cos = rho k |b| cos / |S*| = |a| cos / (2|S*|) = -a:S*/(2|S*|^2) = mu_eff.

The factor 2 is in |b| and in E9; they cancel. S* is deviatoric in every file (dns_fields.py 73-77: trace/3 removed) so a:S* = a:sym(grad U) because a is traceless; the
deviatoric projection changes nothing in the numerator and is required in |S*|.

Favre against Reynolds. stat3 stores <rho u_i u_j> and <rho u_i>. The code forms a_ij = <rho u_i u_j> - rho U_i U_j with U = <rho u>/<rho> (dns_fields.py 56-59).
Since <rho u_i u_j> = rho~ U_i U_j + rho~ u''_i u''_j, this is rho u''_i u''_j exactly (Favre). Checked on random moments: 2.5e-14.
Sign: a carries the + sign of the Reynolds stress, so a Boussinesq stress is a = -2 mu_t S* and cos = +1 for an aligned stress. Checked: |cos - 1| < 7e-16, L = mu/(rho k^1/2) to 8e-16.
Bound: |b| <= sqrt(2/3) = 0.8165 (realizability). Largest cached |b| = L|S*|/k^1/2 is 0.794. No violation.

Numbers.
- Random double-precision tensors, 2e5 samples: max relative residual 8.7e-16.
- Cached 193M b010 points (`meff_planes/part_b010_*.npz`, 5 629 068 points, mu_eff stored from the independent expression -aS/(2 S2), `meff_planes_common.py` 149):
  max relative residual 2.04e-7, median 3.4e-8, 99.99 percentile 1.7e-7. These are float32 rounding.
- The plane arrays (`ref193_b010.npz`, bilinear interpolation of each factor separately) do not satisfy the product to better than the interpolation error
  (median 4e-4, tails to order 1 at sign changes of mu_eff). This is an interpolation property, not an identity failure; the identity holds point by point.

Edge case: at points with a = 0 and k = 0 (walls) cos and L are set to 0 by nan_to_num (dns_fields.py 87-88). They fail the mask (cos > 0) but are interpolated as data in at_points_ref (see Q6, item 4).

Files: `analysis/meff_c50b_c100b/dns_fields.py` 56-59, 71-77, 82-88; `analysis/meff_rans50_c50b/meff_fields.py` 83-105, 119; `figures/code/meff_planes_common.py` 121-150; `q1_identity.py`, `q1_identity.log`.

## Q2. Coupling of length and alignment

Verdict: holds. The length and alignment terms share ln(1/cos_ref) with opposite sign.

Derivation.

    L_eff = |a| / (2 rho k^1/2 |S*|),   L_al = mu_eff/(rho k^1/2) = |a| cos / (2 rho k^1/2 |S*|)   =>   L_eff = L_al / cos.

Equivalently |a|^2 = |a_par|^2 + |a_perp|^2 with |a_par| = |a| cos = 2 mu_eff |S*|. For the RANS: cos = 1 and L_R = mu_t/(rho k^1/2) = L_al,R (meff_fields.py 119). Hence

    length = <ln L_R - ln L_eff,193> = <ln(L_R/L_al,193)> + <ln cos_193>
    align  = <ln 1 - ln cos_193>                       = -<ln cos_193>
    length + align = <ln(L_R/L_al,193)>   (no cos left)
    shared part of each term: |<ln cos_193>| = the alignment term itself.

Numbers (file values; "mine" agrees to 0.02):

| region | length | align | length + align | energy | density | total | shared / length | geo-mean cos_193 |
|---|---|---|---|---|---|---|---|---|
| film 3-10 | -0.963 | +0.806 | -0.157 | -0.351 | +0.007 | -0.501 | 0.84 | 0.45 (exp(-0.806)) |
| film 10-25 | -0.928 | +0.573 | -0.355 | -0.351 | +0.004 | -0.703 | 0.62 | 0.56 |
| jet 0-3 | -0.126 | +0.643 | +0.517 | +0.051 | -0.002 | +0.566 | 5.1 | 0.53 |

Reading the table: the alignment term offsets 84 % (film 3-10) and 62 % (film 10-25) of the length term, and in the jet it is 5 times the length term. This is the
"60-85 %" quoted in the plan, and it is an identity, not a finding. The independent combination is ln(L_RANS/L_al,193): -0.157 and -0.355 in the film, +0.517 in the jet.
Energy (-0.351) is then 70 % of the total in film 3-10 and 50 % in film 10-25. The two-factor split, density + energy + ln(L_R/L_al), is exact and has no alignment term.

Is "the length scale carries most of the film error" a property of the flow or of the definition? Partly each.
- Fact: with L_eff the length term is -0.96 and -0.93, 1.9 and 1.3 times the total. This follows from assigning cos_RANS = 1 and L_eff = L_al/cos: the 193M stress magnitude
  |a|/(2|S*|) = rho k^1/2 L_eff exceeds its aligned part by 1/cos, and the RANS is compared with the magnitude. The net length effect after alignment is -0.16 and -0.36.
- Fact (Q3): an independent dissipation length gives a length term of the same size (-0.89, -1.01), so the existence of a large negative length term with a large positive offset is not an artefact of L_eff alone.
- Interpretation, confidence medium: the offset is the stress that no eddy viscosity can represent (|a_perp|). Naming it "alignment" or "structure" is a labelling choice. The statement "the model length scale is too small
  by 0.4" is a statement about stress magnitude, not about the mixing efficiency of the aligned stress.
- Recommendation: quote the combined term (L_al) and the energy term as the two independent causes, or state the shared part explicitly where length and alignment are shown.

Files: `errorbars.py` 138-143 (terms), 247 (mu from factors); `meff_fields.py` 119; `q234_regions.py`, `q234_regions.log`.

## Q3. Independent dissipation length

Verdict: partly. The 193M dissipation is not in any cache (budgets exist for C50b b006 only; `/mnt/Elements/DNS/analysis` holds only asme2027_figs and nearwall_modes). A new 193M read was not started.
Two estimates were made.

Definitions. L_eps = C_mu k^3/2/eps_m with C_mu = 0.09 and eps_m = eps/rho, eps = <tau'_ij du''_i/dx_j> (tke_budget.py 24, 89), the resolved viscous dissipation. SST gives mu_t/(rho k^1/2) = k^1/2/omega = C_mu k^3/2/eps_R when the a1 limiter is inactive,
so L_eps,R = L_R. omega is not in the RANS export (rans_on_grid.py 40-45), so this equality is assumed.

    ln(mu_t/mu_eff) = ln(rho ratio) + 1/2 ln(k ratio) + ln(L_eps ratio) + ln[ (mu/(rho k^1/2 L_eps))_R / (mu/(rho k^1/2 L_eps))_D ]
    with L_eps,R = L_R:  structure = -<ln(L_al,D/L_eps,D)>;   L_eps term + structure = <ln(L_R/L_al,D)> (the Q2 combination, independent of eps).

So eps only partitions the combination of Q2; it cannot change it.

(a) Exact: RANS 50M against the 50M 3DNS b006, same mesh, no interpolation, same point set and mask (|S*| > 0.05 U/D, cos > 0, k, eps > 0):

| region | total | energy | current length | current align | L_eps term | structure | <ln L_eff/L_eps> | <ln L_al/L_eps> | median L_eff/L_eps | f_V |
|---|---|---|---|---|---|---|---|---|---|---|
| jet 0-3 | +0.576 | +0.027 | -0.121 | +0.674 | -0.131 | +0.684 | -0.010 | -0.684 | 1.005 | 0.56 |
| film 3-10 | -0.563 | -0.377 | -0.968 | +0.775 | -0.888 | +0.695 | +0.079 | -0.695 | 1.007 | 0.54 |
| film 10-25 | -0.760 | -0.336 | -1.029 | +0.601 | -1.013 | +0.585 | +0.016 | -0.585 | 1.017 | 0.49 |
| hole exit | -0.529 | -0.287 | -0.652 | +0.383 | -0.695 | +0.425 | -0.042 | -0.425 | 0.969 | 0.99 |
| hole, diffuser | -0.131 | -0.165 | -0.775 | +0.804 | -0.990 | +1.019 | -0.215 | -1.019 | 0.748 | 0.88 |
| hole, metering | +0.969 | +0.649 | -0.963 | +1.292 | +0.158 | +0.171 | +1.121 | -0.171 | 4.299 | 0.84 |
| approach BL | -0.145 | -0.763 | -1.547 | +2.165 | -1.883 | +2.502 | -0.336 | -2.502 | 0.625 | 0.24 |

The four terms add to the total within 4e-10.

(b) 193M with eps_193 := eps_50 at the same points (proxy; the two are different meshes and batches, eps_50 is a resolved viscous value with budget residual/eps 0.16 in film 3-10, 0.51 in film 10-25, 0.70 in the diffuser, 1.06 in the metering section; C7):

| region | L_eps term | structure | <ln L_eff/L_eps> |
|---|---|---|---|
| jet 0-3 | -0.028 | +0.652 | +0.069 |
| film 3-10 | -0.807 | +0.664 | +0.156 |
| film 10-25 | -1.040 | +0.686 | -0.112 |

Findings.
- Fact (50M, resolved eps): L_eff = C_mu k^3/2/eps holds in the jet and film, as the plan says (E14 text): mean ln ratio -0.010, +0.079, +0.016, medians 1.005, 1.007, 1.017. It holds for L_eff, which contains 1/cos, and not for L_al (median 0.62-0.65, mean ln -0.58 to -0.70). The claim is unverified on the 193M.
- Fact: the L_eps split gives a length term (-0.89, -1.01, -0.13) and an offset (+0.70, +0.59, +0.68) of the same size as the current length and alignment terms (differences below 0.1). Because L_eff is about L_eps, structure = -ln cos - ln(L_eff/L_eps) is about the alignment term.
- Interpretation, confidence medium: in the film the DNS stress magnitude follows the dissipation-based scale, |a| of about 2 C_mu rho k^2 |S*|/eps, and the RANS (k, omega) dissipation length is about 0.4 of it. That is a flow-and-model property. The offset between "magnitude" and "aligned" is a property of the stress tensor (geo-mean cos 0.45). A new independent length does not remove the offset; it renames it.
- Low confidence: the RANS value L_eps,R = L_R is assumed; if the limiter is active L_eps,R > L_R and the L_eps term would be less negative.
- Metering section: eps is not trustworthy there (residual/eps 1.06); ignore L_eps for the hole interior.

Files: `meff_ransp0_c50b/fields_*.npz` (k_d, cos_d, L_d, meff_d), `budgets_c50b/b006/budget_*.npz` (eps), `q234_regions.py` (Q3 block), `q234_regions.log`.

## Q4. Averaging, exclusions, tail, estimator

Verdict: partly holds. In the film the exclusion is small and its bias is small; the alignment term is a tail average; and in the jet the sign of the total depends on the estimator.

Mask (errorbars.py 247-253): mu_t > 0, mu_eff,193 > 0, cos_193 > 0, |S*_193| > 0.05 U/D = 557 s^-1, k_RANS > 0, and k, L, cos > 0 in all ten delete-one samples. Volume weights = stride^3 dV (errorbars.py 65), stride 2.
Counter-gradient points (cos <= 0 gives mu_eff <= 0) are removed because ln(mu_t/mu_eff) is undefined. L_eff itself is defined there.

Volume fractions of each region (mine; nested classes; "used" without the delete-one condition; file f_V in brackets):

| region | non-finite (50M buffer or outside 193M, Q6.1) | |S*| <= Smin | cos <= 0 | used | cos <= 0 share of strained volume |
|---|---|---|---|---|---|
| film 3-10 | 0.000 | 0.446 | 0.015 | 0.539 (0.54) | 0.027 |
| film 10-25 | 0.301 | 0.215 | 0.004 | 0.480 (0.48) | 0.008 |
| jet 0-3 | 0.000 | 0.311 | 0.140 | 0.549 (0.53) | 0.203 |
| hole, metering | 0.016 | 0.023 | 0.114 | 0.846 (0.82) | 0.119 |
| hole, diffuser | 0.008 | 0.000 | 0.104 | 0.887 (0.87) | 0.105 |
| hole exit | 0.000 | 0.000 | 0.016 | 0.984 (0.98) | 0.016 |
| approach BL | 0.353 | 0.393 | 0.027 | 0.226 (0.22) | 0.107 |

Bias of the exclusion, sign. A point with cos <= 0 would have ln(1/cos) -> +infinity, so excluding it lowers the alignment term; it raises the length term by the same amount (Q2).
Re-including cos <= 0 points with cos floored at c0 gives the alignment term (mask A -> c0 = 0.05 / 0.01):
film 3-10 +0.82 -> +0.84 / +0.89; film 10-25 +0.57 -> +0.57 / +0.58; jet +0.72 -> +1.16 / +1.49; metering +1.25 -> +1.41 / +1.62; diffuser +0.84 -> +0.98 / +1.18.
In the film the bias is at most +0.07, which is within the noise of the term. In the jet, metering section and diffuser it is 0.15-0.8 and it has the sign of understating the alignment term.
The delete-one condition acts like an additional exclusion: my jet alignment is +0.72 and the file value is +0.643, total +0.68 against +0.566. Film values agree within 0.02. The jet number is therefore sensitive to the mask at the 0.1 level.

Tail share of the volume-mean alignment term <-ln cos_193> (mine, mask A):

| region | <-ln cos> | volume with cos < 0.3 | share of term from cos < 0.3 | volume cos < 0.5 | share from cos < 0.5 | volume cos > 0.9 | median cos (volume) |
|---|---|---|---|---|---|---|---|
| film 3-10 | 0.820 | 0.266 | 0.646 | 0.417 | 0.816 | 0.007 | 0.621 |
| film 10-25 | 0.574 | 0.123 | 0.421 | 0.241 | 0.610 | 0.005 | 0.723 |
| jet 0-3 | 0.722 | 0.203 | 0.608 | 0.341 | 0.784 | 0.019 | 0.683 |
| hole, metering | 1.253 | 0.434 | 0.739 | 0.643 | 0.897 | 0.002 | 0.352 |

Why the centre-plane median 0.83 and the volume term +0.57 to +0.81 do not conflict (`q4b_where.log`):
- The cached centre plane (x/D 3-10 and 10-25, z/D 0-2, |S*| > Smin) has median cos 0.827 and 0.826, mean -ln cos 0.29 and 0.33, and cos <= 0 fraction 0.000. The "18 %" is not a film number: it is 22.8 % over all cached plane points including the hole and free stream.
- Film volume, median cos by |y|/D (film 3-10): 0-0.5 0.83, 0.5-1 0.82, 1-1.5 0.66, 1.5-2 0.41, 2-3 0.19; by z/D: 0-0.2 0.17, 0.2-0.5 0.47, 0.5-1 0.76, 1-2 0.86. Fraction of strained film volume at |y|/D > 1.5: 36 %; z/D < 0.5: 45 %.
- The centre plane is the best-aligned part of the film. The volume term is dominated by the near-wall layer (z/D < 0.2: median cos 0.17, <-ln cos> 1.9) and the lateral flanks (|y|/D > 2: <-ln cos> 1.7). Cross-plane x/D = 10 (|y|/D < 3, z/D <= 2): median 0.645, <-ln cos> 0.71.
- Fact: the geometric-mean cos of the film volume is 0.44 (3-10) and 0.56 (10-25).

Estimators for the total (mine, same points; ln of ratio):

| region | mean of ln (pipeline) | ln of volume-mean mu ratio | LS gradient-weighted, mask A | LS, all points | LS, |S*| > Smin | older strain-weighted ln ratio (C50b b006, boxes, `meff_ransp0_c50b/regions.txt`) |
|---|---|---|---|---|---|---|
| film 3-10 | -0.485 | -0.152 | -0.448 | -0.431 | -0.437 | -0.49 |
| film 10-25 | -0.702 | -0.206 | -0.429 | -0.349 | -0.355 | -0.36 |
| jet 0-3 | +0.682 | +0.479 | -0.512 | -0.544 | -0.564 | -0.48 |
| hole, metering | +0.782 | +0.145 | +0.928 | +0.957 | +0.957 | +0.98 |
| hole, diffuser | +0.012 | -0.163 | -0.292 | -0.286 | -0.286 | -0.30 |
| hole exit | -0.419 | -0.147 | -0.517 | -0.613 | -0.613 | -0.49 |
| approach BL | -0.460 | +1.973 | +0.084 | +0.167 | +0.117 | +0.48 |

LS = sum(dV |S*|^2 mu_eff)/sum(dV |S*|^2) with each rung's own |S*| (RANS: its own strain), which equals sum(-a:S*)/sum(2|S*|^2) and is the production-weighted viscosity. It needs no threshold and includes counter-gradient points.
- Fact: the film sign (RANS lower) is robust to the estimator; the magnitude is not: far film -0.70 (mean of ln), -0.43 (LS), -0.21 (ln of volume mean). Film 3-10: -0.49, -0.45, -0.15.
- Fact: the jet sign flips. Mean of ln: RANS higher by 0.68 (file 0.57). LS: RANS lower by 0.51, and the older independent strain-weighted C50b analysis gives 0.62 (ln -0.48). The jet carries 23.5 % of its strain weight in counter-gradient points and 31 % of its volume below the strain threshold.
- Fact: approach BL changes from -0.46 to +0.08 (LS) and +1.97 (ln of mean).
- Interpretation, confidence medium: the mean of ln weights every strained, aligned point equally and drops points where the stress works against the strain; the LS value weights by production, where the mu_eff of the 193M is reduced by counter-gradient points. Statements about "RANS mu_t too low/high" in the jet and approach layer should name the estimator. The film statement survives.
- Mismatch inside the plan: E10 states the regional Pr_t as gradient-weighted least squares "so no thresholds are needed"; E13 uses a thresholded volume mean of ln. They are different estimators of the same sort of quantity.

Files: `errorbars.py` 158-181, 247-253; `meff_fields.py` 135-180 (the earlier strain-weighted version, with the covariance term ln Gamma that the E13 mean of ln does not have); `q234_regions.py`, `q4b_where.py` and logs.

## Q5. QCR and the mu_eff split

Verdict: holds. The QCR correction has zero contraction with S*.

Derivation. O_ik = 2 Omega_ik/|grad U| is antisymmetric. With tau symmetric,

    c_ij = O_ik tau_jk + O_jk tau_ik = (O tau)_ij - (tau O)_ij = [O, tau]_ij.

Boussinesq stress tau = 2 mu_t S* - 2/3 rho k I (or its negative; a sign does not matter). The identity commutes with I, so [O, tau] = 2 mu_t [O, S*]. Then

    [O,S*] : S* = tr(O S* S*) - tr(S* O S*) = tr(O S*^2) - tr(O S*^2) = 0    (cyclic trace; S* symmetric).

The correction is also symmetric (([O,S*])^T = [O,S*]) and traceless. Therefore -a_QCR:S* = -a_lin:S*, so mu_eff(QCR) = mu_t exactly, and the whole correction lies in a_perp = a + 2 mu_eff S*, with a_perp = -(correction) up to sign. It enters |a| and cos, and nothing in mu_eff.
Required condition: tau of the linear part must be proportional to S* (plus I). The statement is a property of the QCR form on a linear base; for a stress with a_perp already present, the contraction of a second correction would not vanish.

Numbers (`q5_qcr.log`):
- 1e5 random velocity-gradient tensors: |c:S*|/(|c||S*|) <= 3e-15; |mu_eff(QCR)/mu_t - 1| <= 1.6e-15; the correction is symmetric and traceless to 2e-16.
- Simple shear (U = U(z)): [O,S*] = diag(gamma, 0, -gamma), so a_perp = c_cr1-scaled diag(+1, 0, -1) mu_t gamma (computed: diag(0.614, 0, -0.614) mu_t gamma for c_cr1 = 0.307, x streamwise, z wall-normal), |a_perp|/|a_par| = 2 c_cr1 = 0.600 (c_cr1 = 0.3) or 0.614 (0.307, the measured STAR value, FINDINGS C20), cos = 0.8575 or 0.8522, ln cos = -0.154 or -0.160, mu_eff/mu_t = 1.000000.

Implications.
- Fact: in the split, a QCR RANS has cos < 1 (about 0.85 in simple shear). At fixed mu_t the alignment term would shrink by 0.16 and the length term would grow by 0.16 (Q2 coupling); the total does not change.
- Fact: the mu_eff split cannot see what QCR changes. QCR recovered 72.1 % of the x/D 7-20 eta deficit against the 193M (69.9 % against the 50M; FINDINGS LP6), and a mu_t-only account (energy, length) cannot explain that gain directly, because the local mu_eff of a QCR solution is its mu_t.
- Interpretation, confidence medium: the QCR gain acts through a_perp, i.e. the normal-stress anisotropy (uu - ww) in the momentum equation (lateral spreading, jet-wall interaction), plus the indirect change of the converged mean field (mu_t and k differ between the SST and QCR solutions). The heat flux is unchanged by QCR (gradient diffusion with mu_t/Pr_t, apriori.py 32).
- The 193M has geo-mean cos 0.45 in film 3-10, i.e. |a_perp|/|a_par| = tan(theta) about 2; the simple-shear QCR provides 0.6. The a priori test (FINDINGS C20) gives QCR recovery of the non-aligned anisotropy of 0.05/0.28 (film 3-10, outer/wall) and 0.49/0.32 (film 10-25).
- A direct test, not done: run meff_fields.py on the QCR RANS export (no QCR fields on the 50M grid were found among the caches searched; cells_qcr.csv exists per FINDINGS LS5) and show mu_t, k and cos_QCR against the SST fields.

Files: `analysis/apriori_stress/apriori.py` 8-16, 92; FINDINGS LP6, C20; `q5_qcr.py`, `q5_qcr.log`.

## Q6. Other defects

Each item is labelled fact (F) or interpretation (I) with confidence.

1. Effective extent of the far-film and approach regions (F). Non-finite values come from two sources: the 50M buffer (the RANS fields are set to NaN where the C50b mean tau = 0, meff_fields.py 130-131) and the 193M outside its own domain/buffer. Far film "x/D 10-25": RANS finite only for x/D 10.0-20.4, 193M finite for x/D 10.0-22.6; both finite for 70 % of the region volume (RANS NaN 30.1 %, 193M NaN 16.0 %). Approach layer "x/D -15 to -3": RANS finite for x/D -10.7 to -3, 193M for -12.5 to -3; both finite for 65 % of the volume. The comparison therefore covers x/D 10-20.4 and -10.7 to -3, and f_V (0.48, 0.22) is a fraction of the nominal region. Fix: state the effective range in the table and caption, or define the regions as 10-20 and -10 to -3.
2. E13 versus the regional viscosity (F). E13 is exact for the volume mean of ln (no covariance term), but that is not the mu_R of the methods appendix (strain-weighted; meff_fields.py 135-180 has a covariance factor ln Gamma of -0.24 to -0.64 per region). Both appear in the plan under the same symbol family; see Q4.
3. Units (F). All inputs are SI (S in s^-1, L in m, k in m^2/s^2, mu in Pa s); threshold 0.05 U/D = 557 s^-1 is consistent with S* in s^-1 (errorbars.py 42-43). Terms are dimensionless. No unit defect found.
4. Wall nodes in interpolation (F, I low). dns_fields.py 87-88 sets cos = L = 0 at points where a = k = 0 (wall nodes), and at_points_ref interpolates the factors trilinearly (errorbars.py 86-119) without excluding them. A 50M point whose 193M cell touches a wall node receives a biased-low cos, L and k. The affected volume is the first 193M cell layer only; its effect on the ln means was not quantified (needs a 193M read). Interpretation: small for the film regions, which are dominated by the interior; not tested.
5. Interpolation of nonlinear fields (F, not quantified). cos, L, k, rho and |S*| of the 193M are interpolated separately, and mu_eff for the LS estimator is rebuilt as the product of interpolated factors, so the product identity of Q1 is exact only at the 193M nodes. The ln terms are averaged after interpolation (Jensen: the mean of ln of an interpolated field differs from the interpolation of ln). The 3.14 M-point set at stride 2 samples the 193M at 50M spacing; the size of the effect against a native-mesh average was not measured (the native 193M volume fields are not cached as ln terms).
6. RANS side (F). The RANS values are nearest-vertex copies on the same mesh (rans_on_grid.py docstring), so no RANS interpolation exists. cos = 1 and L_R = mu_t/(rho k^1/2) hold only for the linear SST without QCR (p0dns run). The RANS is compared with a_perp = 0; it is a statement about this RANS only.
7. Jackknife (F). jackknife_fields.py 7-10 forms exact delete-one means of the cumulative sums for N = 10 flow-throughs and errorbars.py 177 uses SE^2 = (N-1)/N sum (t_i - t_bar)^2. The formula is correct. Limits: (i) it measures sampling noise of a 10 FT record only, with flow-throughs treated as independent; (ii) no bias correction for the nonlinear ln, cos and k statistics; (iii) the mask requires positive values in all ten samples, so the result is conditional on that mask (jet shift of 0.12, Q4); (iv) the ladder differences against b006 and the 50M (C25: 0.07-0.29) are larger than the quoted SE of 0.01-0.05, because those include a different reference batch, mesh and point set. Interpretation, confidence medium: the SE of 0.01-0.04 should be read as sampling error only and not as the uncertainty of the regional term. A "resolved when |value| > 2 SE" test therefore passes more easily than the data justify.
8. Weights (F). Volume weights are stride^3 dV of each retained point (errorbars.py 65), index-strided on a stretched grid; each point carries the volume of its own cell cubed up by 8. Unbiased only if dV varies smoothly across 2 cells; not checked against stride 1.
9. Plan statements to adjust (F). (a) "The typical film point on the centre plane is close to aligned (median cos 0.83)" is correct for the centre plane (0.827) and not for the film volume (median 0.62 and 0.72). (b) "L_eff equals C_mu k^3/2/eps within ln 0.1 in the jet and film of the reference": verified on the 50M, not on the 193M. (c) "No covariance term" is true for the mean of ln only.
10. Not a defect, recorded (F). The density term is at most 0.028 and the sum of the four terms reproduces the total to 0.001 (regions_se.txt; fig_meff_split.py residual print).

## Interpretation register

| Statement | Confidence |
|---|---|
| The length-alignment pair is an identity; the net length effect is ln(L_R/L_al) (-0.16 film 3-10, -0.36 film 10-25) | high |
| The dissipation-length scale confirms a RANS length about 0.4 of the DNS in the film (50M), so the large negative length term is not only a definition artefact | medium |
| Alignment and structure are the same offset under two labels | medium |
| The jet sign of the total depends on the estimator; the film sign does not | high (numbers), medium (cause) |
| The film volume term is carried by the near-wall layer and lateral flanks, not by the centre plane | high |
| The QCR gain acts through a_perp and the changed mean field, not through mu_eff | medium |
| Sampling SE understates the uncertainty of the regional terms | medium |

## Open items (needs more than cached data)

1. 193M eps at b010 for the 3DNS (a budget run on C193b; a stat3 read of the film regions only would be enough).
2. omega or the limiter factor from the RANS (add to the export) to confirm L_eps,R = L_R.
3. The QCR RANS fields on the 50M grid, to run meff_fields.py on them.
4. The wall-node interpolation test and the stride-1 weight check (193M or 50M point read).

## Files

Scripts: `audit_extract.py` (point set to scratchpad npz; imports errorbars.py read-only), `q1_identity.py`, `q234_regions.py`, `q4b_where.py`, `q5_qcr.py`.
Logs: `q1_identity.log`, `q234_regions.log`, `q4b_where.log`, `q5_qcr.log`.
Scratch data: `/tmp/claude-1001/-home-sk999-Documents-ASME-2027/f4c92ceb-6df8-4254-96c3-6552d15d0ecf/scratchpad/pts.npz`.
