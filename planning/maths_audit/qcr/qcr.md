# QCR analysis of the RANS 50M stress, 2026-11-27

Status: analysis complete for base, qcr and qcr_prt04. The QCR form was checked against STAR's exported Reynolds stresses: c_cr1 = 0.2999 (section 3).

## 1. Method

- Inputs (read-only): `/mnt/Elements/DNS/analysis/rans50_exports/rans50_{base,qcr,qcr_prt04}_cells_grad.csv`, each read once by `cache_csv.py` (pyarrow, 53,120,567 rows) and cached as float32 `.npy` per column in `/mnt/Elements/DNS/analysis/asme2027_figs/qcr/{base,qcr,qcr_prt04}/`.
- Row order: `perm_<tag>.log` reports identity order (fraction 1.0) for all three; no reordering applied. Row count of each cache: 53,120,567 (`DONE` marker).
- Region logic (`common.py`, from `dns/python/make_tke_inputs.py` REGIONS as used by `regions_pp/build_regions.py`), on cell centres, D = 7.75 mm (`make_tke_inputs.py` line 43, `errorbars.py`): jet 0 <= x/D <= 3 and |y/D| <= 1.5; film 3-10: 3 < x/D <= 10; film 10-20: 10 < x/D <= 20 (domain ends at x/D 20, memory note). All: 0.3 <= z/D <= 2. Differences from the point-wise set: the exit-column exclusion (hole-block cells at 0.3 <= z/D <= 0.5 above the exit) is not applied because hole blocks are not identifiable in the CSV; cells, not the 193M-grid points, are the samples; no stride.
- Wall nodes: cells with k <= 0 or non-finite values are excluded. Removed volume fraction: 0.000 in all regions (cell centres at z/D >= 0.3 carry k > 0).
- Strain: S* = deviatoric part of the symmetric velocity gradient (exported gradients of the same run). Check: exported Smag = sqrt(2 S_ij S_ij) reproduced (full S) to float32 precision.
- Weights: V (cell volume) and V|S*|^2 (S2).
- Linear part a_lin = -2 mu_t S*. QCR part a_qcr = -tau_NL = c_cr1 (O tau - tau O), tau = 2 mu_t S*, O = 2W/|grad v|, c_cr1 = 0.3 (standard QCR2000, equal to the index form O_ik tau_jk + O_jk tau_ik).
- STAR's form (README of the qcr case, User Guide Eqn. 1162): T_NL = -2 mu_t 0.04645 (O S - S O), O = W/|grad v|. Taken literally this is a factor 0.077 below the standard form with c_cr1 = 0.3 (0.04645 against 2 c_cr1 = 0.6). The earlier a-priori work (`analysis/apriori_stress/README.md`) measured (uu - ww)/(4|uw|) = 0.307 in the approach layer of the QCR run, which matches c_cr1 = 0.3. This analysis uses 0.3. Section 3 confirms c_cr1 = 0.2999 against the exported STAR stresses.
- Because a_qcr has zero contraction with S*, a_par = a_lin and a_perp = a_qcr. cos(theta) = |a_lin| / sqrt(|a_lin|^2 + |a_qcr|^2). The base and qcr_prt04 rows of the a_perp columns apply the same formula to those runs' own fields; for base this is the QCR term that would be added to the base field, not a stress the run has.
- Limiter: F2 computed per the README (Menter 2003, beta* = 0.09, a1 = 0.31): arg2 = max(2 sqrt(k)/(beta* omega d), 500 nu/(d^2 omega)), F2 = tanh(arg2^2). L_eps = sqrt(k)/omega (C_mu = beta*), ratio r = mu_t/(rho sqrt(k)) / L_eps = omega mu_t/(rho k). Limiter active: r < 0.999.

## 2. Checks

| check | result |
|---|---|
| a_qcr:S* / (|a_qcr||S*|), max over cells, all regions and runs | 2e-16 to 3e-16 |
| trace of a_qcr / |a_qcr|, max | 3e-16 |
| a_qcr symmetric | yes (to 3e-16) |
| mu_R from the full stress (-a:S*/(2|S*|^2), LS) against mu_R = sum V mu_t |S*|^2 / sum V |S*|^2 | equal; QCR contribution to mu_R: 0 (below 1e-15 relative) |
| F2 recomputed limiter ratio min(1, a1 omega/(S F2)) against the actual omega mu_t/(rho k) | active/inactive agrees in 90-96 % of cells (jet 0.96, film 3-10 0.91, film 10-20 0.93; centre band 0.9995); max difference 0.37 (jet), p99 0.25. Interpretation, confidence medium: the strain measure or the F2 variant used by STAR differs from the one assumed for the plate-near cells (z/D < 1); in the centre band, away from the wall, the match is 0.0007. The limiter numbers below carry this uncertainty. |
| sanity | rho, T and k ranges per region in `results_*.json` field `sanity`; T_vol_mean 285-289 K; k_vol_mean 40-43 J/kg |

## 3. Validation against STAR's modelled stress (2026-11-27)

`validate_star_stress.py` reads 12 byte-offset segments (12 MB each) of `cells_qcr.csv` (Tank, read-only), placed at film-dense row ranges found from a strided sample of the cached X, Z. Rows are matched to the cache by X and Y; alignment is then tested per cell on k (relative, 1e-5) and X (1e-5 D). Two segments aligned (0 misaligned rows) and contained film cells (x/D 0-20, 0.3 <= z/D <= 2, k > 0); the other ten were either misaligned or held no film cells and are excluded (the first run of the script, without the per-cell alignment test, pooled misaligned rows and is superseded). STAR stress: a = rho (u'u' - 2/3 k I) from the exported Reynolds stresses (uu+vv+ww = 2k to 1e-10). Fit of a_STAR - a_lin = c (O tau - tau O), O = 2W/|grad v|, tau = 2 mu_t S*:

- pooled c_cr1 = 0.2999 over 27,034 cells; per segment 0.2999 and 0.2999; per-cell c p10/median/p90 = 0.2997/0.2997-0.2998/0.2999-0.3000.
- relative residual |a_STAR - a_lin - c comm| / |a_STAR - a_lin| = 0.001.
- STAR's a_perp has zero contraction with S* (1e-10, float32-limited).

Fact: the exported STAR stress equals the standard QCR2000 form with c_cr1 = 0.30 (to 0.1 %), in these 27,034 cells. The documented coefficient 0.04645 with O = W/|grad v| does not reproduce it. Limit: two segments, both from the film region; the row ranges are not a random sample. No rescaling of the a_perp rows is needed; section 4 stands as computed with c_cr1 = 0.3.

## 4. QCR stress, QCR run (c_cr1 = 0.3)

|a_perp|/|a_par| is the aggregate sqrt(sum w |a_perp|^2 / sum w |a_par|^2); cos(theta) is the weighted mean and, in brackets, the weighted median. Components are weighted means over rho k of the same weights.

| region | weights | |a_perp|/|a_par| | cos(theta) mean (median) | a_perp xx / rho k | yy | zz | yz (rms) | yy - zz |
|---|---|---|---|---|---|---|---|---|
| jet 0-3 | V | 0.477 | 0.971 (0.997) | 0.062 | -0.007 | -0.055 | 0.038 | 0.048 |
| jet 0-3 | S2 | 0.530 | 0.887 (0.872) | 0.142 | -0.020 | -0.122 | 0.052 | 0.102 |
| film 3-10 | V | 0.583 | 0.902 (0.877) | 0.107 | -0.027 | -0.080 | 0.067 | 0.053 |
| film 3-10 | S2 | 0.589 | 0.860 (0.858) | 0.180 | -0.038 | -0.142 | 0.088 | 0.104 |
| film 10-20 | V | 0.582 | 0.872 (0.864) | 0.126 | -0.050 | -0.076 | 0.064 | 0.026 |
| film 10-20 | S2 | 0.586 | 0.859 (0.857) | 0.185 | -0.050 | -0.136 | 0.084 | 0.086 |
| centre band 3-20, |y/D| <= 0.25 | V | 0.587 | 0.878 (0.863) | 0.162 | -0.022 | -0.140 | 0.042 | 0.118 |

Facts: the aggregate ratio is 0.58-0.59 in the film, close to 2 c_cr1 = 0.6 (the simple-shear value, cos = 0.857). The yz mean is zero by symmetry; its rms is given. The mean a_perp is mainly xx > 0 and zz < 0 (wall-normal stress reduced), with yy slightly negative; yy - zz is positive in every region. The 193M volume medians (0.62 near film, 0.72 far film, 0.83 centre plane) are below the QCR-run medians (0.877, 0.864, 0.863 in the same order of regions), i.e. the QCR run has a smaller a_perp share than the 193M reference. Interpretation, confidence medium: the constant is the run's (section 3).

## 5. Mean-field changes, QCR against base (ratio qcr/base; qcr_prt04/base in brackets)

| region | mu_R [Pa s] base | mu_R ratio | k (vol. mean) ratio | omega_x rms ratio | v rms ratio | sum V (dT/dy)^2 ratio |
|---|---|---|---|---|---|---|
| jet 0-3 | 2.361e-3 | 0.81 (0.81) | 0.94 (0.95) | 0.81 (0.80) | 0.90 (0.91) | 0.86 (0.79) |
| film 3-10 | 2.032e-3 | 0.82 (0.78) | 0.94 (0.95) | 0.93 (0.88) | 0.86 (0.88) | 0.92 (0.56) |
| film 10-20 | 1.971e-3 | 0.74 (0.73) | 0.96 (0.97) | 1.50 (1.46) | 1.19 (1.15) | 1.01 (0.41) |

omega_x rms (base, qcr, qcr_prt04), units U/D: jet 0.0995, 0.0810, 0.0799; film 3-10 0.0693, 0.0643, 0.0610; film 10-20 0.0346, 0.0519, 0.0504. Lateral velocity v rms [m/s] (base, qcr, qcr_prt04): jet 3.78, 3.42, 3.43; film 3-10 1.04, 0.90, 0.91; film 10-20 0.487, 0.581, 0.562.

Facts: mu_R, from the run's own mu_t and |S*|^2 weights, is lower in the QCR run by 18-26 %. The mean field differs between the runs, so this is a change of the field, not of the closure at fixed field (at fixed field the QCR part changes mu_R by 0, section 2). k falls by 3-6 %. The streamwise vorticity and lateral velocity rise only in film 10-20 (+50 % and +19 %); in the jet and film 3-10 they fall by 7-19 % and 10-14 %. The lateral temperature gradient integral falls in the jet and film 3-10 and is unchanged in film 10-20 for QCR; Pr_t 0.4 reduces it by 44 % (film 3-10) and 59 % (film 10-20) against base.

## 6. Base run: dissipation length and limiter

| region | L_eps/D (S2-weighted mean) base / qcr / prt04 | r = mu_t/(rho sqrt(k) L_eps), LS-weighted base (qcr, prt04) | r, S2-weighted mean base | limiter-active volume fraction base (qcr) | strain-weight fraction base (qcr) |
|---|---|---|---|---|---|
| jet 0-3 | 0.038 / 0.034 / 0.034 | 0.919 (0.907, 0.906) | 0.883 | 0.369 (0.332) | 0.549 (0.608) |
| film 3-10 | 0.031 / 0.027 / 0.027 | 0.930 (0.932, 0.937) | 0.931 | 0.143 (0.143) | 0.752 (0.781) |
| film 10-20 | 0.036 / 0.029 / 0.029 | 0.964 (0.941, 0.943) | 0.954 | 0.117 (0.158) | 0.667 (0.789) |

Facts: r = 1 in the volume median of every region; r is 0.88-0.95 when weighted by strain. The limiter is active in 12-37 % of the volume and in 55-75 % of the strain weight. Interpretation, confidence medium (limiter state carries the F2 uncertainty of section 2): the limiter lowers the regional RANS viscosity by 4-8 % relative to the dissipation-length value, and the effect sits in the high-strain cells.

## 7. Assessment

1. Fact: the QCR part is orthogonal to S* (2e-16), so at a given mean field mu_eff is unchanged and only a_perp is added; a_perp/a_par is 0.58-0.59 in the film and 0.53 in the jet (S2 weights).
2. Fact: the QCR run differs from base in mu_R (-18 to -26 %), k (-3 to -6 %), and in the vorticity and lateral velocity (film 10-20 up, jet and film 3-10 down); mu_R is therefore not unchanged between the runs. Interpretation, confidence low to medium: the a_perp stress (yy - zz > 0 in every region) is consistent with a secondary-flow contribution in the far film, but the jet and near-film omega_x fall, so the data do not show a single mechanism in all three regions; the lower mu_R in QCR may come from the changed mean field (k, strain) and not from the closure.
3. Fact: STAR's exported stress is reproduced with c_cr1 = 0.2999 (section 3), so the a_perp values use the constant of the run. Limit: 27,034 cells from two segments.
