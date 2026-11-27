# Definitions of the near-wall coolant layer height: literature precedent

Prepared 2026-11-27. The folder did not exist at the start, so the work began from scratch.
Verification key: VERIFIED = definition read in a fetched page, PDF or abstract (the text was extracted and inspected).
PARTIAL = existence and abstract read, definition itself not read. UNVERIFIED = from search snippet or memory only.
Heights in hole diameters D for a shaped hole at x/D 5-20 are not available from any source read, except the one
datum in section 4. All other heights must be computed from the fields.

## (a) Table of definitions

| # | Definition | Formula / rule | Source | Verified? | Same surface for both runs possible? | Fields needed | Remarks |
|---|---|---|---|---|---|---|---|
| 1 | Thermal boundary-layer thickness, 99 % | delta_T = y_99 such that T(x,y_99) = 0.99 T_inf (textbook form; a wall-temperature-referenced form is (T-Tw)/(Tinf-Tw) = 0.99) | Wikipedia, "Thermal boundary layer thickness and shape": "customarily defined as the point in the boundary layer, y99, where the temperature T reaches 99% of the free stream value" | VERIFIED (encyclopedic, secondary) | Yes if applied to the approach layer (independent of the run); for the film itself theta never reaches 0.99 inside the coolant | Approach T or U profile | Defined for an attached boundary layer without coolant. Not defined for a layer whose wall value is set by the coolant. |
| 2 | Momentum boundary-layer thickness delta_99 of approach flow | y where U = 0.99 U_inf | Common usage; used as reference length in Schroeder and Thole 2016 (they report delta_2/D = 0.14, H = 1.45, Re_delta2 = 670 at Tu 0.5 %; delta_99 not stated) | VERIFIED for the reported quantities; delta_99 itself UNVERIFIED | Yes, independent of the simulation | Approach profile at a fixed upstream station | delta_99/D for that rig is not given. A turbulent estimate from delta_2 and H would be an interpretation, not a source value. |
| 3 | Inner layer | y/delta < 0.1 | Pope, Turbulent Flows (2000), Fig. 7.13 (Cornell figure file shows lines y/delta = 0.1 and 0.3); University of Iowa lecture notes on Pope ch. 7: "Inner layer usually defined as y/delta < 0.1" | VERIFIED (Pope figure file and secondary lecture notes; the printed book text was not read) | Yes (delta from approach flow) | Approach delta | Valid for equilibrium boundary layers without film injection. |
| 4 | Viscous wall region / buffer / log-law limits | viscous sublayer y+ < 5; buffer 5 < y+ < 30; "viscous wall region y+ < 50"; log law y+ > 30 | Wikipedia "Law of the wall" (sublayer, buffer); Iowa notes on Pope: "y+ < 50 viscous wall region" | VERIFIED (secondary sources) | Only if y+ is fixed from one run. y+ varies with local wall shear, which differs between RANS and the reference and is altered by the film | u_tau (reference or approach), nu | A fixed y+ gives a height that varies in x. Wall shear under a film differs from that of a plain boundary layer. |
| 5 | Wall-jet inner layer / jet half-width | y_m = height of maximum velocity; outer layer extends from y_m to the outer edge; y_1/2 = height above y_m where U = U_max/2 | ERCOFTAC KBwiki UFR 3-10, "The plane wall jet": "inner layer (from the wall to ym, the position of maximum velocity)"; wall jet defined from Launder and Rodi (1981, 1983) | VERIFIED (ERCOFTAC page); Launder and Rodi UNVERIFIED | Yes if y_m is taken from the reference velocity field | Mean streamwise velocity | Requires a velocity maximum above the wall; a film with a mainstream velocity exceeding the coolant velocity may have none. |
| 6 | Mixing-region half-thickness (slot film, NASA TN D-7636) | half-width at which 0.99 of the stream-jet velocity difference is reached | Brown, NASA TN D-7636 (1974), p. on mixing region: "choosing for the half-width the distance from the slip line at which 0.99 times the velocity difference between stream and jet is reached" | VERIFIED (scanned OCR text) | Yes from the reference | Velocity | Slot geometry and 1/7-power boundary-layer assumptions. |
| 7 | Wieghardt integral thermal thickness (slot film) | delta_T = integral from 0 to infinity of (T - T_inf)/(T_ad,w - T_inf) dy = integral of theta dy | Voegele, U. Maryland PhD 2013, Eq. 9, describing Wieghardt: "a form similar to the more traditional enthalpy thickness"; similarity profile theta ~ exp(-0.768 (y/(delta_T + ...))^(13/6)) (exponent form garbled in the extraction, not relied on) | VERIFIED (thesis text); Wieghardt original UNVERIFIED | Yes: computed from the reference profile at each x, then used as a fixed surface in both runs | Mean T (or theta), adiabatic-wall T | Integral thickness is not a height at which theta has fallen to a fixed fraction. Ratio to the 0.1 or 0.5 heights depends on the profile shape. |
| 8 | Enthalpy thickness | Energy-balance form, integral of U (T - T_inf)/[U_inf (T_w - T_inf)] dy (curved-wall variant also given) | Furuhama, Moffat, Johnston, Kays, NASA CR-174964 (1985), section 4.2.2: "Enthalpy thickness can be calculated by two independent methods" | PARTIAL: section and method read; the printed formula is garbled in OCR and was not reconstructed with certainty | Yes from the reference | Mean U and T | Convective-flux weighted, so conceptually matches a coolant-flux balance. Used there as an energy-closure check. |
| 9 | Potential-core / effectiveness threshold | length where effectiveness = 0.99 | Voegele 2013: "approximated by the location where the effectiveness is 0.99" | VERIFIED | Streamwise length, not a height | Mean theta | Shows a 0.99 threshold used in film cooling, but for streamwise length. |
| 10 | Jet extent by concentration isocontours | contours of mean concentration 0.2, 0.4, 0.6, 0.8 plotted; "region downstream of injection where the concentration is greater than 5%" used for an uncertainty analysis | arXiv 1908.03540 (LES-enriched MRI data of inclined jets in crossflow), read via ar5iv | VERIFIED. Authors not read from the page. | Yes if the contour is taken from the reference | Mean concentration | No formal film-thickness definition in that paper; thresholds are used descriptively. |
| 11 | Scalar threshold as fraction of local maximum (jet in crossflow) | jet fluid = region where C > 1 % of the local scalar maximum on trajectory-normal planes; "A change in this cutoff value changes the curves quantitatively, but the trends remain the same" | Muppidi and Mahesh, JFM 598, 335-360 (2008), doi 10.1017/S0022112007000055, Sec. 3.3 | VERIFIED (full text) | Yes, if the fraction is applied to the reference | Mean scalar, trajectory | Only fraction of local maximum; appendix compares other cutoffs. |
| 12 | Control volume with a streamtube from the hole | "volume-averaged turbulent diffusivity is calculated using a mass-flux balance across the streamtube emanating from the jet hole" | Coletti, Benson, Ling, Elkins, Eaton, Int. J. Heat Fluid Flow 43, 149-160 (2013) | PARTIAL: abstract read at experts.umn.edu; full text not accessed, so the way the streamtube is closed is not known | Yes, if the stream surface is taken from the reference mean flow | Mean velocity, concentration | Stream surface carries zero mean-flow flux by construction, so only turbulent flux crosses it. |
| 13 | Bounding streamline of the partially mixed film | boundary = streamline separating coolant from mainstream; "the position of the hypothetical bounding streamline that preserves the condition m2/m1 = 1" | Ornano and Povey, J. Turbomachinery 142(4), 041002 (2020), doi 10.1115/1.4046277, Appendix B, Fig. 32 | VERIFIED (full text read) | Yes from the reference | Mean mass-flux distribution through a station | Used for a mixing model with constant pressure and no heat transfer across the boundary. |
| 14 | Control-volume basis of slot-film correlations | correlations derived from control-volume analyses, then tuned | Kodzwa and Eaton, Stanford TF-93 (2005): "Goldstein (1971) documented a set of successful correlations for slot injection film cooling. Such correlations were based on control volume analyses" | VERIFIED as a statement about Goldstein; Goldstein 1971 itself UNVERIFIED (not accessible) | n/a | n/a | Goldstein's film-thickness definition could not be read. |

## (b) Options

Each option lists the rule, citations, how both runs use the same surface, and weaknesses.

**Option 1. Fraction of the approach boundary-layer thickness.**
Rule: h = k delta_99 (approach boundary layer at a fixed upstream station), with k = 1 or 0.1.
Citations: row 1 and 2 (99 % thickness), row 3 (inner layer y/delta < 0.1, Pope 2000 via Cornell figure and Iowa notes).
Identical evaluation: delta_99 comes from the approach profile, independent of the simulation; the height is constant in x and z.
Weaknesses: unrelated to the coolant layer; a shaped-hole film can lie below or above it; delta_99 of the two runs must match, otherwise the approach profile of one run is used for both; delta_99 for the Schroeder-Thole rig is not given in the paper read.

**Option 2. Fixed wall-unit height.**
Rule: h = y+ value (30 or 50) based on u_tau of the approach flow or of the reference.
Citations: rows 3 and 4.
Identical evaluation: u_tau taken once, from the reference approach flow, converted to a height in D.
Weaknesses: the film changes local wall shear; a fixed y+ gives a very small height compared with the coolant layer; it is a layer of the wall flow, not of the coolant. Height in D must be computed.

**Option 3. Isosurface of mean scalar defined on the reference.**
Rule: h(x,z) is the height at which time-mean theta falls to a fraction f (0.01 or 0.1 of the local maximum, or 0.5 of the wall value) in the reference field.
Citations: row 11 (1 % of local maximum, Muppidi and Mahesh 2008), row 10 (5 % concentration threshold and contours at 0.2-0.8), row 7 and 9 (0.99 threshold on effectiveness).
Identical evaluation: surface extracted once from the reference and used as a fixed geometric surface (interpolated onto the RANS mesh) for both runs.
Weaknesses: depends on f (Muppidi and Mahesh report quantitative change with cutoff); the surface is wavy in z and has a slope, so the flux through it needs the surface normal; the RANS field can lie partly outside the surface, so the flux through the surface is not equal to the coolant leaving the layer in RANS terms.

**Option 4. Integral thickness from the reference.**
Rule: h(x) = delta_T = integral of theta dy (Wieghardt form), or the enthalpy-flux form, of the spanwise-averaged reference profile; optionally a multiple of it.
Citations: rows 7 and 8 (Voegele 2013 Eq. 9; Furuhama et al. 1985).
Identical evaluation: computed from the reference spanwise-averaged profile at each x, used for both.
Weaknesses: integral thickness is not a bounding height; its relation to the height where theta = 0.1 depends on profile shape; Wieghardt similarity applies to slot flow far downstream, not directly to a shaped-hole film at x/D 5-20 (rows 7 and 9 note boundary-layer assumptions fail near injection).

**Option 5. Stream surface carrying the coolant mass flow (streamtube or bounding streamline).**
Rule: boundary = mean-flow stream surface from the hole exit, or the streamline above which a stated fraction of the coolant mass flow lies.
Citations: row 12 (Coletti et al. 2013 streamtube mass-flux balance), row 13 (Ornano and Povey 2020 bounding streamline).
Identical evaluation: surface constructed in the reference mean flow and used for both runs. Mean-flow flux through it is zero by construction, so only the turbulent part remains.
Weaknesses: it makes the mean-flow part vanish and so removes the split asked for; stream surface of a time-mean 3D flow is not closed by entrainment; Coletti et al. closure details were not read. In RANS the same surface has a non-zero mean flux because the flow differs.

## (c) One numerical datum on height

Schroeder and Thole, J. Turbomachinery 138(9), 091001 (2016) (shaped hole, DR 1.5): at x/D = 4, velocity and turbulence contours show the jet extending "up to a height of y/D = 1.2" at M = 1.5 and 1.7 at M = 3 (laterally z/D within about +-1.2). VERIFIED (text extracted from the PSU-hosted PDF). That is a velocity-based extent, not a theta threshold, and at x/D 4 not 5-20. No source read gives theta-based heights at x/D 5-20. Those must be computed.

## Interpretation (own ranking, labelled as such)

Confidence is low to moderate because the film-cooling definitions (Goldstein 1971, Stollery and El-Ehwany, Kutateladze and Leont'ev, Haydt and Lynch, Eberly and Thole, Kohli and Bogard) could not be read.
1. Option 3 (isosurface from the reference) has the closest precedent for a jet or film scalar layer (Muppidi and Mahesh 2008, with the cutoff sensitivity reported there). Sensitivity to f should be shown.
2. Option 4 (integral thickness) is standard in the slot-film literature and is a single well-defined number per x, though less common for discrete holes.
3. Option 1 (approach delta_99) is the only definition independent of both runs; it has literature precedent as a wall-layer measure but not as a coolant layer.
4. Option 5 has precedent in control-volume analyses (Coletti 2013, Ornano and Povey 2020), but it removes the mean-flow part of the split.
5. Option 2 gives layers too thin for the coolant at x/D 5-20 (interpretation; to be checked by computing y+ = 30-50 in D).

## Not verified / not found

- Goldstein 1971 (Advances in Heat Transfer, "Film cooling"): no accessible full text; film-thickness definition not read.
- Stollery and El-Ehwany; Kutateladze and Leont'ev: only secondhand mentions (Voegele 2013 says Stollery and El-Ehwany used "a simple energy and mass balance").
- Haydt and Lynch 2018-2019, Eberly and Thole 2014 (the Eberly and Thole paper found is J. Turbomachinery 136(6) 061003, 2013, "Time-resolved film-cooling flows at high and low density ratios"; no trajectory definition found): UNVERIFIED.
- Kohli and Bogard 2005: no listing found. Kohli and Bogard measured turbulent transport at M = 0.4 (via secondary citation): UNVERIFIED.
- Yuan, Street, Ferziger 1999 (JFM 379, 71-104) and Yuan and Street 1998 (Phys. Fluids 10, 2323): existence confirmed from Muppidi and Mahesh reference text and search; entrainment definition not read. UNVERIFIED. Muppidi and Mahesh state their entrainment is computed on trajectory-normal planes with a 1 % scalar threshold, and that Yuan and Street's results "correspond to the value of scalar concentration".
- Ryan, Milani, Ling, Eaton: Milani, Ling, Eaton, JFM 906, A27 (2021), arXiv 2001.04600, was downloaded and searched; it contains no film-thickness or control-volume definition (grep on thickness / isosurface / threshold / streamtube returned only inflow boundary-layer remarks).
- Bogard and Thole 2006 review (J. Propulsion Power 22(2), 249-270, doi 10.2514/1.18034): citation confirmed by search only; text not read.

## (d) Sources (URLs)

- Muppidi and Mahesh 2008, JFM 598: https://dept.aem.umn.edu/~kmahesh/globalassets/publpdf/journal/suman_jfm_scalar.pdf ; https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/direct-numerical-simulation-of-passive-scalar-transport-in-transverse-jets/252BF2CEC3B7EBBDAA98B8C88F095039
- Coletti et al. 2013 (abstract): https://experts.umn.edu/en/publications/turbulent-transport-in-an-inclined-jet-in-crossflow/
- Ornano and Povey 2020 (author copy): https://ora.ox.ac.uk/objects/uuid:e6e3f712-91ab-4dd9-86e4-a86cbe5dbfa8/files/rnv935308n
- Voegele 2013 thesis: https://drum.lib.umd.edu/bitstreams/fa1e3223-6a27-4a7c-9de7-da771273101a/download
- Furuhama et al. 1985, NASA CR-174964: https://ntrs.nasa.gov/api/citations/19860000994/downloads/19860000994.pdf
- Brown 1974, NASA TN D-7636: https://ntrs.nasa.gov/api/citations/19740012483/downloads/19740012483.pdf
- Kodzwa and Eaton 2005, TF-93: https://web.stanford.edu/group/tfsa/TF_reports/TF-093_Kodzwa.pdf
- Schroeder and Thole 2016 (PSU-hosted): https://bpb-us-e1.wpmucdn.com/sites.psu.edu/dist/c/92759/files/2018/08/2016-Schroeder-JT-1-27gruaw.pdf
- Schroeder and Thole 2017 thermal field (abstract only): https://pure.psu.edu/en/publications/thermal-field-measurements-for-a-shaped-hole-at-low-and-high-free/
- Milani, Ling, Eaton: https://arxiv.org/abs/2001.04600
- LES-enriched MRI data: https://ar5iv.arxiv.org/html/1908.03540
- ERCOFTAC UFR 3-10: https://kbwiki.ercoftac.org/w/index.php/UFR_3-10
- Wikipedia, thermal boundary layer: https://en.wikipedia.org/wiki/Thermal_boundary_layer_thickness_and_shape ; law of the wall: https://en.wikipedia.org/wiki/Law_of_the_wall
- Pope Fig. 7.13 (Cornell): https://pope.mae.cornell.edu/TurbulentFlows/popefigures/pdf/chapter7Figures.pdf ; Iowa lecture notes on Pope ch. 7-8: https://stern.lab.uiowa.edu/sites/stern.lab.uiowa.edu/files/2025-11/Chapter-8-Pope-Part-1-final_cm.pdf
- Bogard and Thole 2006 (citation only): https://pure.psu.edu/en/publications/gas-turbine-film-cooling/
