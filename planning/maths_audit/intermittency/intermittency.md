# Intermittency detector: methods note

Status: detector replaced in every figure script that computed gamma from the 3DNS / LES samplers; threshold policy 'band'
(owner decision, 27 Nov). Numbers are from `intermittency_policy.txt` (script `intermittency_policy.py`, same folder) and the
logs of the figure scripts.

## 1. Definition

d = (v - <v>_T)^2 / U_inf^2 at every sampler point, with

- v: spanwise (y) velocity of the sampler records, U_inf = 86.2812 m/s, D = 7.75 mm;
- <v>_T: centred running mean in time over T = 2 D/U (box, edge value repeated at the ends of the record);
- d then smoothed by a centred running mean of 1 D/U;
- laminar / turbulent by Otsu's method on log10 d (pooled histogram, 1200 bins from -12 to 0);
- gamma = fraction of time with d above the threshold (over the batch window; over span only where stated).

No normalisation by gamma_max. Filters are defined in D/U. The 193M sampler records are unevenly spaced (0.08 to 0.38 D/U, one
crash repeat dropped); both runs are resampled linearly in time to the uniform step 0.1 D/U (median record spacing of both
runs) before filtering. Samplers: z010 (plane z/D 0.1, x/D -3 to 20.5 read, y/D -3 to 3) and x10 (plane x/D 10).
Windows: 193M C193b batch b010 (3992 records, 3999 resampled), 20M LES C20bW batch b006 (2399 records).
Code: `figures/code/_gamma_extract_v.py`. Cache: `/mnt/Elements/DNS/analysis/asme2027_figs/gamma_v/` (the 5.7 GB staging copy of v
is local, in the session scratch area, because the NTFS volume was slow under load).
The first and last 1 D/U of a record use the repeated edge value.

## 2. Difference from Dellacasagrande et al., J. Turbomach. 145(4) 041015 (2023), p. 3

Same: indicator built from the spanwise velocity fluctuation squared, classification by Otsu's method on log10 of the indicator.
Different: the mean removed from v is a running mean in time (T = 2 D/U) at each point; their mean is taken over span, which is
not applicable because this flow is not homogeneous in y. A 1 D/U running mean smooths d. No normalisation by gamma_max.

## 3. Threshold policy 'band'

One threshold for the whole paper: Otsu on log10 d pooled over the intermittent band of the z/D 0.1 plane, |y/D| 2 to 3 (y/D 2..3
and -3..-2, one band across the periodic boundary), x/D 0 to 20, T = 2 D/U, from the 193M; the same value applied to both runs.

| run, band |y/D| 2 to 3 | log10 threshold | threshold | eta (Otsu) | valley / low peak | Sarle BC | skewness | excess kurtosis |
|---|---|---|---|---|---|---|---|---|
| 193M (used) | -4.510 | 3.09e-5 | 0.745 | 0.892 | 0.520 | -0.07 | -1.07 |
| 20M LES, own value | -4.570 | 2.69e-5 | 0.736 | 0.970 | 0.510 | -0.15 | -1.00 |

Histogram peaks (log10 d): 193M -5.62 and -3.42; LES -5.62 and -3.67. The LES own threshold is 0.87 times the 193M threshold.
Sarle BC is below 0.555 for both runs (0.555 is the value for a uniform distribution); the valley / low-peak ratio is 0.89 and 0.97.
The two-mode structure is weak; the classification rests on the Otsu criterion (eta 0.74).

Other pooling domains, same procedure (T = 2 D/U), for reference: z010 plane x/D 0..20 all y: 193M -4.08, LES -4.06 (Sarle BC 0.61, 0.54);
mid-pitch line y/D 3, x/D 0..20: 193M -4.79, LES -4.84 (Sarle BC 0.56, 0.53); band |y/D| 1.5 to 3: 193M -4.35, LES -4.40;
band |y/D| 2.5 to 3: 193M -4.68, LES -4.73.

Required checks, policy 'band', T = 2 D/U:

| run | approach x/D -3 to -1.5, all y: mean / p99 / max | core |y/D| < 1, x/D 5 to 20: mean / min / fraction > 0.95 |
|---|---|---|
| 193M | 0.0001 / 0.0024 / 0.0058 | 0.993 / 0.954 / 1.000 |
| 20M LES | 0 / 0 / 0 | 0.977 / 0.879 / 0.857 |

Facts: the approach-layer check (gamma < 0.02) holds for both runs. The core check (gamma > 0.95) holds for the 193M (all points
above 0.95); for the LES the mean is 0.977 and 14 % of the core points lie below 0.95 (minimum 0.879, at the centreline near x/D 20: the
y/D 0 value is 0.980, 0.991, 0.962, 0.898 at x/D 5, 10, 15, 20). Earlier policies for comparison (same checks): plane-pooled threshold
(193M value for both runs): core mean 0.969 / 0.926 (193M / LES), min 0.881 / 0.741; mid-pitch-line-pooled threshold: core mean 0.997 / 0.991,
min 0.974 / 0.931, approach max 0.050 / 0.018.

## 4. Sensitivity (mid-pitch line y/D 3, T = 2 D/U unless stated), gamma at x/D 5 / 10 / 15 / 20

Threshold factor (both runs):

| factor | 193M | 20M LES |
|---|---|---|
| x0.5 (log10 -4.81) | 0.007 / 0.293 / 0.833 / 0.970 | 0.002 / 0.337 / 0.878 / 0.981 |
| x1 (log10 -4.51) | 0.002 / 0.184 / 0.744 / 0.937 | 0.000 / 0.173 / 0.770 / 0.942 |
| x2 (log10 -4.21) | 0.000 / 0.102 / 0.626 / 0.868 | 0.000 / 0.087 / 0.623 / 0.886 |

Band limits (threshold recomputed from the 193M band, applied to both runs):

| band | log10 th | 193M | 20M LES |
|---|---|---|---|
| \|y/D\| 1.5 to 3 | -4.35 | 0.001 / 0.141 / 0.679 / 0.905 | 0.000 / 0.124 / 0.705 / 0.913 |
| \|y/D\| 2 to 3 | -4.51 | 0.002 / 0.184 / 0.744 / 0.937 | 0.000 / 0.173 / 0.770 / 0.942 |
| \|y/D\| 2.5 to 3 | -4.68 | 0.005 / 0.251 / 0.796 / 0.957 | 0.001 / 0.263 / 0.838 / 0.969 |

Running-mean length T (band threshold recomputed for each T): T = 1 D/U (log10 th -5.06): 193M 0.000 / 0.108 / 0.639 / 0.892, LES
0.000 / 0.094 / 0.629 / 0.890; T = 4 D/U (-3.99): 193M 0.013 / 0.275 / 0.778 / 0.935, LES 0.002 / 0.326 / 0.811 / 0.921.

Statement: gamma at x/D 10 changes from 0.29 to 0.10 (193M) and from 0.34 to 0.09 (LES) between threshold x0.5 and x2; the x/D at which gamma
first exceeds 0.5 is 12.8 (193M) and 11.5 (LES) at the policy threshold. The 193M / LES difference at the policy threshold
(0.01 at x/D 10, 0.03 at x/D 15) is smaller than the effect of a factor 2 in threshold. The figures show the x0.5 to x2 range as a shaded band
around gamma(x) (fig_spots panel e, fig_eta_x panel c).

## 5. Scripts switched (figures/code)

| script | gamma | action |
|---|---|---|
| fig_wall_state.py | 3DNS / LES detector gamma (panels e, f); RANS gamma (panel d) is the transported model variable | switched (e, f); (d) unchanged |
| fig_bl_profiles.py | 3DNS / LES gamma from the x10 plane; RANS gamma model variable | switched; RANS unchanged |
| fig_eta_x.py | gamma_mid of the 193M (z010 band \|y/D\| 2.81 to 3), shaded x0.5 to x2 | switched |
| fig_spots.py | detector d and gamma; RANS 50M model gamma line in panel (e) | rebuilt |
| fig_spots_dudt_stats.py, fig_spots_dudt_old.py | old \|du/dt\| detector, record of the previous fig_spots | kept unchanged |
| fig_spots_stats_v.py | spot celerities, new detector | audit script, not called by any figure |
| fig_spots_z010.py, fig_spots_wall.py, _spots_z010_extract.py | skin-friction estimate; no gamma | none |
| fig_budget_lines.py, fig_budget_maps.py | RANS transported gamma (SST gamma-Re_theta); only named in text | none |
| fig_bubble_motion.py | gamma^2 = magnitude-squared coherence | none |
| _gamma_extract.py | old detector helper | kept unchanged |

New helper: `_gamma_extract_v.py`. The earlier policies remain selectable through `POLICY` ('i', 'ii', 'iii', 'band', 'band1.5', 'band2.5').

## 6. RANS curve in fig_spots panel (e)

The green dash-dot line is the RANS 50M transported intermittency, the model variable of the SST gamma-Re_theta transition model, taken
from the z/D 0.1 layer of the cell export (same loader and cache as fig_wall_state.py panel d, line y/D 3). It is a different quantity from the
detector gamma and carries no threshold band. Values at x/D 5, 10, 15, 20: 0.020, 0.020, 0.021, 0.021. Height dependence: the cache holds the
prism layers 0.083, 0.090, 0.098, 0.105, 0.113 only; at mid-pitch, x/D 15 the value is 0.021 in each. z/D 0.05 and 0.2 are not cached.

## 7. Spot celerities: not reported

Celerities from edge slopes of connected regions are not reported. With the policy threshold the selection rule (born narrow, below 5 D wide,
at least 6 D/U) leaves N = 1 to 2 regions per run, and the values change with the threshold (policy ii: N 6 and 13; line policy: N 1 and 2). The plan
views show that the edges of the wedge between the jets meander, so the mid-pitch (x, t) regions are not isolated patches. `fig_spots_stats_v.py`
is kept as an audit script; the figures do not call it. The patch marked in fig_spots (triangle) is the connected region of the (x, t) field with the
largest downstream travel in the batch; its marker follows the centre of its x-extent at each snapshot time.

## 8. Old versus new values

Mid-pitch gamma, line y/D 3, x/D 5 / 10 / 15 / 20 (old |du/dt| -> line-pooled policy -> band policy):

| run | old | line policy | band policy |
|---|---|---|---|
| 193M | 0.008 / 0.205 / 0.670 / 0.892 | 0.006 / 0.286 / 0.828 / 0.968 | 0.002 / 0.184 / 0.744 / 0.937 |
| LES | 0.005 / 0.265 / 0.736 / 0.928 | 0.002 / 0.324 / 0.870 / 0.979 | 0.000 / 0.173 / 0.770 / 0.942 |

First x/D with mid-pitch gamma above 0.8: 193M 17.3 -> 14.3 -> 15.7; LES 16.8 -> 12.9 -> 15.3.
fig_eta_x band |y/D| 2.81 to 3 (193M), x/D 5 / 10 / 15 / 20: old 0.010 / 0.224 / 0.681 / 0.886, band policy 0.003 / 0.202 / 0.757 / 0.938;
first x/D above 0.8: old 17.2, line policy 14.3, band policy 15.5.
Plate mean (x/D -3 to 20, |y/D| < 3, fig_wall_state): 193M 0.702 -> 0.752 -> 0.724; LES 0.690 -> 0.759 -> 0.726; RANS model gamma 0.167.
x10 plane, gamma at z/D 0.1 / 0.5 (old -> line -> band): mid-pitch 193M 0.202 -> 0.285 -> 0.183 / 0.606 -> 0.881 -> 0.793;
mid-pitch LES 0.264 -> 0.322 -> 0.172 / 0.579 -> 0.857 -> 0.726; centreline z/D 0.1: 193M 0.961 -> 0.999 -> 0.996, LES 0.964 -> 0.998 -> 0.992.
