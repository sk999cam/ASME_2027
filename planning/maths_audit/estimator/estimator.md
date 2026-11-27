# Regional estimator for the effective eddy viscosity split, 2026-11-26

Status: updated with the owner decisions (see the update section); sections 1-7 below describe the previous run (mode A). Recommended estimator: gradient-weighted least squares (LS) for the regional value, with a Shapley split into density, energy, L_al and strain-weight terms. Read-only on 4dns_filmcooled and Tank. No new 193M read.
Format note: text deliverable for agents and the audit folder; no PDF built.

Notation. a = anisotropic Reynolds stress (Favre, positive sign), S* = deviatoric sym(grad U~), |.| = Frobenius norm, V = cell volume, p = -a:S*.
Per point: mu_eff = p / (2|S*|^2) = rho k^1/2 L_eff cos, L_al = mu_eff/(rho k^1/2) = L_eff cos. RANS: cos = 1, mu_t = rho k^1/2 L_R.
Regional LS value: mu_R = sum_R V p / sum_R 2V|S*|^2. Each side uses its own |S*| (RANS: its own strain).

## Update 2026-11-26 (owner decisions): wall nodes excluded, z/D >= 0.3 on plate regions

Note (owner rule 2026-11-26): the paper uses far film x/D 10-20, not 10-25 as written below. `estimator.py pts.npz C xyz.npz 20` restricts the region to x/D <= 20 (log `estimator_modeC_x20.log`; far-film total -0.855 against -0.860). The numbers below are the 10-25 run.

Decisions applied: (1) 4-player Shapley (rho, energy, L_al, strain weights) is the headline; L_eff and align are secondary (5 players, with order ranges). (2) Wall nodes excluded: k_193 = 0 (cos set to 0) or k_RANS <= 0 (0 points). (3) Plate regions (jet 0-3, film 3-10, film 10-25, approach layer) use z/D >= 0.3, as the regional Pr_t (E10). Hole regions: wall nodes excluded, no wall-distance cut yet.
The sections below this one are the previous run (mode A, no exclusion); their numbers stand for mode A only. Modes: `estimator.py pts.npz A|B|C xyz.npz`; logs `estimator_modeA.log` (previous), `estimator_modeB.log` (wall nodes excluded, no z cut), `estimator_modeC.log` (wall nodes excluded and z cut, the decided estimator).
z/D is the third coordinate of the 50M grid point (coordinates re-read from the 50M grid by `extract_xyz.py`, order verified against pts.npz; no 193M read). Hole-region z/D is negative inside the hole and is not used.

### a. Per-region table, decided estimator (mode C)

Terms: 4-player Shapley, +- delete-one jackknife over spatial blocks (not time batches). Columns B give the total with wall nodes excluded and no z cut, for comparison.

| region | total (C) | energy | rho | L_al | S-wt | total (B, uncut) | total (A, previous) | order range energy | order range L_al |
|---|---|---|---|---|---|---|---|---|---|
| jet 0-3 | -0.328+-0.028 | -0.011+-0.060 | +0.019+-0.019 | -0.367+-0.029 | +0.030+-0.076 | -0.455+-0.180 | -0.544 | [-0.06, +0.04] | [-0.48, -0.26] |
| film 3-10 | -0.390+-0.389 | -0.132+-0.246 | +0.045+-0.046 | -0.285+-0.010 | -0.018+-0.101 | -0.441+-0.458 | -0.431 | [-0.21, -0.09] | [-0.44, -0.16] |
| film 10-25 | -0.860+-0.524 | -0.220+-0.213 | +0.018+-0.025 | -0.423+-0.097 | -0.235+-0.191 | -0.424+-0.454 | -0.349 | [-0.37, -0.16] | [-0.75, -0.20] |
| approach BL | +1.109+-0.068 | +0.180+-0.028 | +0.000+-0.000 | +1.235+-0.063 | -0.306+-0.027 | +0.161+-0.090 | +0.167 | [+0.06, +0.31] | [+0.97, +1.51] |
| hole exit (wall nodes only) | -0.449+-0.071 | -0.138+-0.096 | +0.040+-0.008 | -0.251+-0.054 | -0.100+-0.139 | -0.449 | -0.613 | [-0.18, -0.11] | [-0.37, -0.15] |
| hole, diffuser (wall nodes only) | -0.284+-0.090 | -0.148+-0.010 | +0.006+-0.001 | -0.002+-0.047 | -0.139+-0.055 | -0.284 | -0.286 | [-0.16, -0.14] | [-0.06, +0.06] |
| hole, metering (wall nodes only) | +0.958+-0.044 | +0.184+-0.017 | -0.010+-0.000 | +1.040+-0.058 | -0.257+-0.009 | +0.958 | +0.957 | [+0.18, +0.19] | [+1.03, +1.05] |

Secondary 5-player split of L_al (mode C), with the minimum and maximum marginal change over 120 orders:

| region | L_eff | range | align | range | L_eff + align | 4-player L_al |
|---|---|---|---|---|---|---|
| jet 0-3 | -0.557 | [-0.67, -0.43] | +0.199 | [+0.15, +0.23] | -0.358 | -0.367 |
| film 3-10 | -0.632 | [-0.94, -0.46] | +0.356 | [+0.24, +0.50] | -0.276 | -0.285 |
| film 10-25 | -0.801 | [-1.24, -0.54] | +0.372 | [+0.27, +0.50] | -0.429 | -0.423 |
| approach BL | -0.514 | [-0.79, -0.26] | +1.734 | [+1.62, +1.81] | +1.220 | +1.235 |
| hole exit | -0.567 | [-0.73, -0.45] | +0.322 | [+0.29, +0.37] | -0.246 | -0.251 |
| hole, diffuser | -0.630 | [-0.80, -0.49] | +0.627 | [+0.54, +0.75] | -0.003 | -0.002 |
| hole, metering | -0.854 | [-1.02, -0.71] | +1.888 | [+1.74, +2.06] | +1.034 | +1.040 |

LS values after the cut (mu_R,193 / mu_R,RANS, Pa s): jet 3.339e-3 / 2.405e-3; film 3-10 2.752e-3 / 1.864e-3; film 10-25 4.229e-3 / 1.789e-3; approach 8.485e-5 / 2.573e-4. Volume removed by the z cut (beyond wall nodes): 0.077, 0.140, 0.141, 0.196. Non-finite volume: 0, 0, 0.301, 0.353 (buffer, as before). Wall-node volume: 0.003-0.005 per region.
Residual of the exact additivity: below 3e-15 in all modes.

Facts and interpretation.
- Fact (`zcut_weight.log`): z/D < 0.3 holds 0.91 (jet), 0.91 (film 3-10), 0.96 (film 10-25) and 0.94 (approach) of the 193M strain weight, on 8-20 % of the volume (RANS side: 0.91, 0.88, 0.90, 0.94). The cut therefore removes the viscous layer from the LS weights. After the cut the estimator is carried by the remaining 4-10 % of the strain weight.
- Fact: signs after the cut: jet, film 3-10, film 10-25 and hole regions as in mode B (negative), metering positive. The approach layer flips: +0.161 (B) to +1.109 (C). The approach-layer total depends on the cut; its signed counter-gradient share of sum(w mu) is -0.134 after the cut (before: -0.035), and 44 % of its volume has cos <= 0. Interpretation, confidence medium: the approach-layer sign should not be stated without naming the cut.
- Fact: film 10-25 changes from -0.424 to -0.860 and its energy term from -0.50 to -0.22 (the strain-weight term from +0.65 to -0.24). The wall-layer weights that produced the large positive S-wt term in the previous run are gone. The energy term of the jet (-0.157 to -0.011) and film 3-10 (-0.320 to -0.132) shrinks too; the film L_al (-0.285 to -0.423) is now the larger term. Interpretation, confidence medium: the earlier "energy and L_al explain the film" reading was partly a wall-layer result; after the cut L_al is the leading term in jet and film 3-10 and film 10-25 (energy -0.22 is half of it).
- Fact: the order dependence of the 4-player terms narrows after the cut (film 3-10 L_al range 0.28, previously 0.41).

### b. Hole regions

- STAR export with wall distance d: `/mnt/Elements/DNS/analysis/rans50_exports/` does not exist yet (checked 2026-11-26). Hole values above are wall-nodes-excluded only, without a d cut. Finish when it lands.
- Proposal: d/D >= 0.3, d = distance to the nearest solid wall (plate or hole wall). Reason (two lines): z/D = 0.3 in E10 is a distance from the wall in units of D, and the cut removes the viscous layer, whose thickness depends on the wall shear and not on the wall orientation; the same distance on the hole wall keeps the cut consistent between regions. Check once d lands that d/D = 0.3 corresponds to the same y+ as z/D = 0.3 on the plate (wall shear from the RANS export), and scale if it does not.
- Fact (`hole_proxy.py`, indicative only: cylinder of radius 0.5 D about the PCA axis of the metering points, axis at 30.0 degrees; this is not the STAR wall distance): d >= 0.3 D keeps 0.160 of the metering volume and 2.2e-5 of its 193M strain weight; d >= 0.2 D keeps 0.363 and 5.7e-5; d >= 0.1 D keeps 0.650 and 3.8e-4. The metering total of +0.958 is a wall-layer quantity (more than 99.9 % of the strain weight lies within 0.1 D of the wall). A d cut will change the metering LS value completely and leave a core with very small strain. Interpretation, confidence medium: for the metering section either state that the regional LS value is a wall-layer value, or use a smaller cut (d >= 0.1 D) and report the weight kept. The diffuser and exit were not tested with the proxy (non-cylindrical).

### c. Block dominance after the cut (mode C)

| region | blocks | largest block share of 193M strain weight | total | delete-one range | total without the largest block |
|---|---|---|---|---|---|
| jet 0-3 | 9 | 0.23 (block 18) | -0.328 | [-0.356, -0.319] | -0.319 |
| film 3-10 | 8 | 0.32 (block 7) | -0.390 | [-0.809, -0.313] | -0.313 |
| film 10-25 | 6 | 0.40 (block 10) | -0.860 | [-1.385, -0.652] | -0.730 |
| approach BL | 6 | 0.32 (block 15) | +1.109 | [+1.047, +1.133] | +1.124 |
| hole exit | 5 | 0.26 (block 4) | -0.449 | [-0.509, -0.406] | -0.420 |
| hole, diffuser | 5 | 0.33 (block 4) | -0.284 | [-0.315, -0.193] | -0.315 |
| hole, metering | 5 | 0.40 (block 4) | +0.958 | [+0.919, +0.973] | +0.938 |

Fact: the jet is the least block-dependent region after the cut (SE 0.028). Film 3-10 and film 10-25 still depend on single blocks (delete-one range 0.5 and 0.7 wide). All totals keep their sign under every block deletion. The block-based SE is a spatial sensitivity only.

### d. Comparison with the E10 regional Pr_t (`figures/code/fig_prt_regions.py`, `fig_prt_classes.py` `prep`/`sums`/`cls`/`ratio`, `prt_corr_common.py`)

Match:
- Estimator form: regression through the origin with squared-gradient weights, no |gradient| threshold. E10 "tensor" form `pt` = sum(mu_eff S^2)/sum(S^2) is the same expression as mu_R. E10 uses the (dU/dz)^2 and (dT/dz)^2 weights for its classic Pr_t (headline) and also computes the tensor and lateral forms.
- Wall cut: `m = (x >= x0) & (x < x1) & (z >= 0.3)` with z in D (`fig_prt_regions.py` line 42, 44). Same z/D >= 0.3 as here.
- Favre means and 193M b010 (and b006).

Differences:
1. Domain: E10 sums over points of the centre plane y = 0 (plus the cross-plane x/D = 10) interpolated from the 193M; this estimator sums over the 3-D volume of the 50M points (|y|/D up to 3 in the film). Plate z/D < 0.3 is cut in both, but the lateral flanks are only in the volume estimator.
2. Weights: E10 sums plane points as stored (no cell-volume weight; bins used only for the bootstrap); this estimator uses V|S*|^2 with V the 3-D cell volume (stride 2). On a stretched grid the two weightings differ.
3. Gradient: E10 classic uses the single wall-normal component (-rho u''w'' dU/dz over (dU/dz)^2); this estimator uses the full deviatoric tensor, -a:S* over 2|S*|^2. The E10 tensor form `pt` is closer.
4. Regions: E10 far film is x/D 10-20 (cross-plane x/D = 10 separate); this estimator far film is 10-25 with effective RANS range 10-20.4 (30 % non-finite). E10 has no approach layer or hole regions.
5. Wall nodes: E10 masks only non-finite values (`ok`); it does not exclude k = 0 nodes explicitly, but z/D >= 0.3 removes them on the plate. Here they are excluded explicitly and for the hole.
6. E10's mu_eff and the plane S are bilinear interpolations of each factor from the 193M (Q6.5 in `audit.md`); this estimator interpolates the 193M factors at the 50M points likewise.
Confidence high for the items read; not run numerically.

### e. Time-batch jackknife hook

`estimator.py` accepts `DELETE_ONE_FILES="f1.npz f2.npz ..."` (one file per delete-one sample: `c193_rho`, `c193_k`, `c193_L`, `c193_cos`, `c193_S`, shape (N,) at the points of pts.npz, same order). It recomputes the total and the 4-player terms for each sample with the full-record point mask fixed, and prints SE = sqrt((N-1)/N sum (t_i - mean)^2) in section G. Tested on three synthetic scaled copies only (runs; the SEs are meaningless). The per-batch data format from fc is not known: if the fields come per batch (not delete-one), form the delete-one record as (N record - batch) / (N - 1) first; if the cached cos and L are per-batch statistics they cannot be averaged directly, and delete-one must be formed from the batch sums of the stress and strain before cos and L are built.

## Verdicts

| Item | Verdict | Basis |
|---|---|---|
| 1a mu_R minimises sum V\|a + 2 mu S*\|^2 | Holds | derivation; numerical check below |
| 1b production identity | Holds, for the deviatoric production only | derivation; residual 0 to 7e-16 |
| 1c mu_R = strain-squared-weighted arithmetic mean of mu_eff | Holds | derivation; 1e-16 |
| 1d same principle as E10 | Holds in form; weights differ between the two LS values of Pr_t | derivation; E10 code not read |
| 2(i) weighted geometric means | Rejected: gap to ln(mu_R ratio) is +1.47 and +2.36 in the film | numbers below |
| 2(ii) Shapley split of the ratio of sums | Recommended | exact (residual < 3e-15), no mask |
| 2(iii) hierarchical (reference weights, then weight term) | Usable alternative; differs from (ii) by up to 0.15 in the film and 0.27 for the weight term in film 10-25 | numbers below |
| 3 literature | Pointwise LS "optimal eddy viscosity" verified (Wu, Xiao, Paterson 2018; McConkey et al. 2022). Regional (volume-integrated) LS: no verified source; stands on the derivation | section 3 |
| 4 table | Done; jackknife is over 22 spatial blocks, not over time batches | section 4 |

Fact (new, corrects `audit.md` Q4): the "23.5 % of jet strain weight in counter-gradient points" is wall nodes. cos <= 0 weight share in the jet is 0.235, of which 0.230 is k_193 = 0 nodes (cos set to 0 by `nan_to_num`, dns_fields.py 87-88). Strict cos < 0 holds 0.005 of the jet strain weight (`cos_neg_weight.log`). Only hole exit (0.107), approach layer (0.042) and diffuser (0.029) carry more than 1 %. Counter-gradient points change sum(w mu) by at most -0.035 (approach layer).

## 1. Maths

### 1a. Frobenius least squares

For a region R and one scalar mu,

    J(mu) = sum_R V |a + 2 mu S*|^2 = sum V|a|^2 + 4 mu sum V a:S* + 4 mu^2 sum V|S*|^2
          = sum V|a|^2 - 4 mu sum V p + 4 mu^2 sum V|S*|^2,
    dJ/dmu = 0  =>  mu_R = sum V p / (2 sum V|S*|^2)   (d2J/dmu2 = 8 sum V|S*|^2 > 0, a unique minimum).

Minimum value J_min = sum V|a|^2 - (sum V p)^2 / sum V|S*|^2. The fraction of sum V|a|^2 that one regional mu can represent is (sum V p)^2 / (sum V|S*|^2 sum V|a|^2).
The residual a_perp,R = a + 2 mu_R S* is orthogonal to S* in the weighted sense: sum V a_perp:S* = 0.

Numerical (`estimator.py` section A):
- Random traceless symmetric tensors, 2e4 points: J(mu_R(1 +- e)) - J(mu_R) = 3.39e-7, 3.39e-5, 3.39e-3 for e = 1e-3, 1e-2, 1e-1, symmetric on both sides (parabola, minimum at mu_R).
- 193M points, scalar form (|a| = 2 rho k^1/2 L_eff |S*|, p = 2 mu_eff |S*|^2): a 9-point scan over mu_R(1 +- 0.2) has its minimum at mu_R in all seven regions (scan argument / mu_R = 1.000).
- Interpretation, confidence high: the explained fraction is low, 0.036 (film 3-10), 0.024 (film 10-25), 0.053 (jet), 0.227 (hole exit). It is dominated by the normal-stress anisotropy, which no scalar viscosity represents. mu_R is the best single value; it is not a good model of the tensor. This is the same statement as the geometric-mean cos of 0.45 in `audit.md` Q2.

### 1b. Production of k

Compressible Favre flow:

    P_k = -rho u''_i u''_j dU~_i/dx_j = -(a_ij + 2/3 rho k delta_ij) dU~_i/dx_j = -a:S* - 2/3 rho k div U~.

The first term is the deviatoric production, which mu_R reproduces exactly:

    sum V(-a:S*) = sum V 2 mu_R |S*|^2   (residual 0 to 7e-16 in all regions).

- The isotropic part -2/3 rho k div U~ is not represented by mu_eff. It vanishes for incompressible flow. At the jet Mach number its size was not computed (div U~ is not in pts.npz). Confidence: unverified; state "deviatoric production" in the text.
- rho enters only through a = rho u''u'' (already in a). No extra rho weight appears in sum V p. A specific (per-mass) fit of b = a/(rho) would weight by rho^2 |S*|^2; it is a different estimator and was not computed.
- Sign: p < 0 points lower mu_R exactly as they lower production.

### 1c. Strain-squared-weighted mean

With w_i = V_i|S*_i|^2:

    mu_R = sum w_i mu_eff,i / sum w_i,   mu_eff,i = p_i / (2|S*_i|^2).

Checked on random tensors to 1e-16 and on all regions (column `sum(w mu)/sum(w)` equals mu_R to 4 digits).
Consequences:
- Points with |S*| = 0 have w = 0; no division and no |S*| threshold. The weight share below the current threshold 0.05 U/D is 0.000-0.001 in every region, so the threshold of the present mask removes volume but almost no production.
- Points with cos <= 0 stay in the sum with mu_eff <= 0.
- Weight is concentrated. 0.144 (film 3-10), 0.155 (film 10-25) and 0.230 (jet) of the 193M strain weight sit on k_193 = 0 nodes (0.3-0.5 % of the volume, |S*| 1e5 s^-1 against 7e3 s^-1 elsewhere in the jet). Both models have mu close to 0 there (RANS mu_t about 6e-9 Pa s in the jet nodes), so those nodes enter the numerators as 0 and the denominators with 14-23 % of the weight. Interpretation, confidence medium: this is a property of the estimator; it is the production-weighted regional value including the viscous-wall layer. Effect on the total of removing k_193 = 0 nodes: jet -0.544 to -0.455, film 3-10 -0.431 to -0.441, film 10-25 -0.349 to -0.424, hole exit -0.613 to -0.449, others < 0.01. The sign does not change.

Comparison with the thresholded mean of ln (current, `errorbars.py` 247-253): the mean of ln weights each strained, aligned point by volume; it needs mu_eff > 0, cos > 0, |S*| > 0.05 U/D, and all ten delete-one samples positive. LS needs no mask. Total ln ratio, same points:

| region | LS (this) | mean of ln (current) | LS - mean of ln |
|---|---|---|---|
| jet 0-3 | -0.544 | +0.682 | -1.226 |
| film 3-10 | -0.431 | -0.485 | +0.054 |
| film 10-25 | -0.349 | -0.702 | +0.353 |
| hole exit | -0.613 | -0.419 | -0.194 |
| hole, diffuser | -0.286 | +0.012 | -0.298 |
| hole, metering | +0.957 | +0.782 | +0.175 |
| approach BL | +0.167 | -0.460 | +0.627 |

Fact: the film sign is the same in both. The jet sign flips (+0.68 against -0.54); it stays negative under every delete-one block (range -0.628 to -0.437) and without wall nodes (-0.455).

### 1d. Consistency with E10

The Pr_t fit has the same structure. For scalar flux f = rho u''theta'' and gradient g = grad T~, minimising sum V|f + alpha g|^2 gives alpha_R = -sum V f.g / sum V|g|^2, Pr_t,R = mu_R/alpha_R. Both factors are LS values. They have different weights (|S*|^2 and |g|^2), so Pr_t,R is a ratio of two regressions through the origin, not one joint fit. Interpretation, confidence medium: E10 and E13 now share the principle; the E10 script was not read in this task (stated in the plan, `audit.md` Q4).

## 2. An additive split compatible with LS

### Why the ratio does not split in ln

mu_R = N/D with N = sum V|S*|^2 rho k^1/2 L_eff cos and D = sum V|S*|^2. N is a sum of products and the factors vary together inside the region, so ln N has no term-by-term form. Two exact routes exist.

### (i) Weighted geometric means: rejected

Define ln mu_G = sum w ln mu_eff / sum w with w = V|S*_193|^2 on both sides. Per point ln(mu_t/mu_eff) = ln(rho ratio) + 1/2 ln(k ratio) + ln(L_R/L_eff) - ln cos is exact, so the four terms add to ln(mu_G,RANS/mu_G,193) exactly. It needs mu_eff > 0 (cos > 0, k > 0) and the geometric mean is not the LS value.

| region | ln(mu_R ratio) | sum of the terms | gap (Jensen + mask) | share of strain weight kept |
|---|---|---|---|---|
| jet 0-3 | -0.544 | -0.352 | -0.192 | 0.765 |
| film 3-10 | -0.431 | -1.900 | +1.469 | 0.853 |
| film 10-25 | -0.349 | -2.709 | +2.360 | 0.844 |
| hole exit | -0.613 | +0.353 | -0.966 | 0.825 |
| hole, diffuser | -0.286 | +1.037 | -1.322 | 0.951 |
| hole, metering | +0.957 | +1.449 | -0.492 | 0.885 |
| approach BL | +0.167 | -2.846 | +3.013 | 0.915 |

Fact: the gap exceeds the total in six of seven regions. The geometric mean is dominated by low-mu points that the production weighting does not see. Verdict: (i) cannot be the split of an LS total.

### (ii) Shapley split of the ratio of sums: recommended

State of a region: choose for each of the players rho, k, L_eff, cos, strain weights (|S*|) the RANS or the 193M value, and compute f = ln(N/D) of that mixed state (cos = 1 and mu_t = rho k^1/2 L_R hold exactly for the RANS state, relative error 1.8e-7 float32). The full change is f(all RANS) - f(all 193M) = ln(mu_R,RANS/mu_R,193). The Shapley value of a player is its marginal change in f averaged over all orders of substitution. By the efficiency property (Shapley 1953) the terms add exactly to the total; the numerical residual is below 3e-15. It has no mask, keeps cos <= 0 points, and treats the players symmetrically, so no order has to be chosen.
(The Shapley attribution rule is standard cooperative game theory; no citation was searched or verified for it.)

Order dependence is quantified by the minimum and maximum marginal change over all orders (120 orders for five players, 24 for four).

Headline table, 4 players (rho, energy = 1/2 ln k, L_al = L_eff cos, strain weights). Terms +- delete-one block jackknife SE (section 4):

| region | total | energy | rho | L_al | S-wt | order range energy | order range L_al |
|---|---|---|---|---|---|---|---|
| jet 0-3 | -0.544+-0.176 | -0.157+-0.089 | +0.012+-0.020 | -0.286+-0.026 | -0.113+-0.095 | [-0.21, -0.10] | [-0.42, -0.16] |
| film 3-10 | -0.431+-0.436 | -0.320+-0.361 | +0.030+-0.042 | -0.334+-0.132 | +0.192+-0.110 | [-0.46, -0.22] | [-0.56, -0.15] |
| film 10-25 | -0.349+-0.438 | -0.500+-0.265 | +0.005+-0.023 | -0.582+-0.211 | +0.728+-0.073 | [-0.74, -0.35] | [-1.00, -0.26] |
| hole exit | -0.613+-0.099 | -0.138+-0.096 | +0.040+-0.008 | -0.251+-0.054 | -0.264+-0.100 | [-0.18, -0.11] | [-0.37, -0.15] |
| hole, diffuser | -0.286+-0.091 | -0.148+-0.010 | +0.006+-0.001 | -0.002+-0.047 | -0.141+-0.058 | [-0.16, -0.14] | [-0.06, +0.06] |
| hole, metering | +0.957+-0.041 | +0.184+-0.017 | -0.010+-0.000 | +1.040+-0.058 | -0.257+-0.010 | [+0.18, +0.19] | [+1.03, +1.05] |
| approach BL | +0.167+-0.090 | -0.468+-0.017 | -0.000+-0.000 | +0.675+-0.060 | -0.040+-0.023 | [-0.79, -0.03] | [+0.27, +1.19] |

5 players (L_eff and cos separate). L_al is reported as the sum L_eff + align:

| region | energy | rho | L_eff | align | L_eff + align | S-wt | order range L_eff | order range align |
|---|---|---|---|---|---|---|---|---|
| jet 0-3 | -0.171 | +0.012 | -0.749 | +0.478 | -0.271 | -0.114 | [-0.98, -0.57] | [+0.35, +0.64] |
| film 3-10 | -0.351 | +0.029 | -0.965 | +0.661 | -0.304 | +0.194 | [-1.51, -0.59] | [+0.37, +1.02] |
| film 10-25 | -0.523 | +0.004 | -1.301 | +0.759 | -0.543 | +0.713 | [-1.98, -0.73] | [+0.44, +1.18] |
| hole exit | -0.144 | +0.040 | -0.567 | +0.322 | -0.246 | -0.264 | [-0.73, -0.45] | [+0.29, +0.37] |
| hole, diffuser | -0.152 | +0.005 | -0.630 | +0.627 | -0.003 | -0.135 | [-0.80, -0.49] | [+0.54, +0.75] |
| hole, metering | +0.184 | -0.010 | -0.854 | +1.888 | +1.034 | -0.251 | [-1.02, -0.71] | [+1.74, +2.06] |
| approach BL | -0.487 | -0.000 | -1.369 | +2.056 | +0.687 | -0.033 | [-2.09, -0.63] | [+1.75, +2.39] |

Facts:
- Grouping L_eff and cos into L_al narrows the order range of the other terms (film 3-10 energy [-0.74, -0.13] with five players, [-0.46, -0.22] with four). L_eff and align alone carry order ranges of 0.6-1.2 in the film, which is larger than the terms' own difference from the sum L_al. They are large and opposite because L_eff = L_al/cos (`audit.md` Q2); only their sum is determined to about +-0.15.
- The 5-player sum L_eff + align differs from the 4-player L_al by 0.03 (film 3-10), 0.04 (film 10-25), 0.015 (jet), 0.003 (hole exit); 0.006-0.012 in the hole.
- The strain-weight term is not small. Film 10-25 +0.713 (5 players) and +0.728 (4 players). It states that the RANS strain is distributed differently inside the region from the 193M strain (e.g. in far film the RANS strain weight sits where mu_t is higher). Interpretation, confidence low-medium: part of it follows from the effective RANS range x/D 10-20.4 (30 % buffer volume, `audit.md` Q6.1).

### (iii) Hierarchical split

Level 1: ln(mu~_R,RANS/mu_R,193), with mu~ the RANS value on the 193M weights V|S*_193|^2, split by a 4-player Shapley over rho, k, L_eff, cos (24 orders). Level 2: the weight term ln(mu_R,RANS/mu~_R,RANS), the remainder. Exact by construction.

| region | total | energy | rho | L_eff | align | L_al | S-wt |
|---|---|---|---|---|---|---|---|
| jet 0-3 | -0.544 | -0.187 | +0.016 | -0.697 | +0.514 | -0.184 | -0.190 |
| film 3-10 | -0.431 | -0.328 | +0.029 | -0.873 | +0.692 | -0.181 | +0.049 |
| film 10-25 | -0.349 | -0.490 | +0.004 | -1.104 | +0.785 | -0.319 | +0.456 |
| hole exit | -0.613 | -0.139 | +0.041 | -0.495 | +0.328 | -0.167 | -0.348 |
| hole, diffuser | -0.286 | -0.154 | +0.005 | -0.567 | +0.616 | +0.050 | -0.186 |
| hole, metering | +0.957 | +0.187 | -0.010 | -0.816 | +1.850 | +1.035 | -0.255 |
| approach BL | +0.167 | -0.413 | -0.000 | -1.216 | +2.080 | +0.864 | -0.284 |

Differences from (ii): up to 0.12 in the film for energy, L_al and align (film 3-10 L_al -0.181 against -0.334), and 0.26 for the weight term in film 10-25 (+0.456 against +0.713). The hierarchical choice of the 193M weights as base is a convention; Shapley has no such choice.

Not usable (derived, not run): an LMDI (logarithmic-mean Divisia) split of the sums needs every term positive, which fails for cos <= 0 points (mu_eff <= 0). Citation not searched.

### Recommendation

Use the 4-player Shapley split of ln(mu_R,RANS/mu_R,193) into density, energy, L_al and strain-weight terms, and show L_eff and align as a further 5-player split of L_al where the alignment is discussed. It is exact without any mask or threshold, it keeps counter-gradient points like the production does, and it needs no choice of substitution order; the residual order dependence (listed above) is the quantified uncertainty of the split itself.

Interpretation, confidence medium, for the paper text: the regional total in the film is -0.43 and -0.35; energy and L_al (-0.32 and -0.33 in film 3-10; -0.50 and -0.58 in film 10-25) explain it, the density term is below 0.03, and the weight term has opposite sign (+0.19, +0.73). The current statement "the length scale carries the film error" is not supported by the LS split: L_eff alone is -0.97 and -1.30, but its offset by align is +0.66 and +0.76, and only the combination L_al (-0.30, -0.54 with five players) is determined.

## 3. Literature

Verification method: DOI resolved through Crossref (title, authors, journal, volume, pages match), and for the LS content the arXiv full text was read (`pdftotext`).

| Reference | DOI | Status | Content relevant to LS |
|---|---|---|---|
| Wu, J.-L., Xiao, H., Paterson, E., Physical Review Fluids 3, 074602 (2018), "Physics-informed machine learning approach for augmenting turbulence models: A comprehensive framework" | 10.1103/PhysRevFluids.3.074602 | Verified (arXiv 1801.02762 full text read; DOI and journal via search) | Sec. 2.1, Eqs. (4)-(5), arXiv page 7: "we introduce an optimal eddy viscosity that minimizes the discrepancy between the anisotropy Reynolds stress tensor and its linear part", nu_t^L = arg min \|\|b - nu_t S\|\| with the Frobenius norm, "computed by projecting the anisotropy stress tensor on the strain rate tensor". Sec. 3.2: nu_t^L from DNS stress and strain. Pointwise, not regional. Their b and S carry their own scaling (factor 2 in Eq. 5); not identical in notation. |
| Wu, J.-L., Xiao, H., Sun, R., Wang, Q., J. Fluid Mech. 869, 553-586 (2019), "RANS equations with explicit data-driven Reynolds stress closure can be ill-conditioned" | 10.1017/jfm.2019.205 | Verified (arXiv 1803.05581 text contains "optimal eddy viscosity", lines read by grep only; context not read in full) | Uses the optimal eddy viscosity as the linear part. Already in `references.bib` (key listed below). |
| McConkey, R., Yee, E., Lien, F. S., Physics of Fluids 34, 035110 (2022), "Deep structured neural networks for turbulence closure modeling" | 10.1063/5.0083074 | Verified (Crossref; arXiv 2201.01710 read) | Sec. II.A, Eq. (5): "The optimal eddy viscosity nu_t* is obtained as the solution of the following least-squares approximation problem", nu_t* = arg min \|\|a - (-2 nu_t S)\|\|; Eq. (6), closed form nu_t* = -(1/2)(a:S)/(S:S). Same convention as this study (a with positive stress sign, factor 2). Adds a non-negative variant (Eq. 7). Attributes the approach to Wu et al. 2018. |
| Schmitt, F. G., C. R. Mécanique 335(9-10), 617-627 (2007), "About Boussinesq's turbulent viscosity hypothesis: historical remarks and a direct evaluation of its validity" | 10.1016/j.crme.2007.08.004 | DOI and abstract verified (Crossref, publisher page) | Abstract: the hypothesis "corresponds to an alignment between the Reynolds stress and mean strain tensors"; an indicator "between 0 and 1" (French abstract) is introduced; "this hypothesis is almost never verified". The indicator's equation was not read (article PDF not retrieved): unverified that it equals cos. |
| Ling, J., Kurzawski, A., Templeton, J., J. Fluid Mech. 807, 155-166 (2016) | 10.1017/jfm.2016.615 | DOI verified | Anisotropy from a neural network with embedded invariance; no LS eddy viscosity found in what was read (title/metadata only). Not a source for the LS estimator. |
| Milani, P. M., Ling, J., Saez-Mischlich, G., Bodart, J., Eaton, J. K., J. Turbomach. 140(2), 021006 (2018) | 10.1115/1.4038275 | DOI and abstract verified | Machine-learned turbulent diffusivity assuming gradient diffusion; no LS fit was read. Unverified for LS. Milani 2019 (10.1115/1.4045389) and 2020 (10.1017/jfm.2020.820) were not read. |
| Weatheritt, J., Sandberg, R., J. Comput. Phys. 325, 22-37 (2016), "A novel evolutionary algorithm applied to algebraic modifications of the RANS stress-strain relationship" | 10.1016/j.jcp.2016.08.015 | DOI verified (Crossref) | Gene-expression programming of anisotropy; LS content not read. Unverified for LS. |
| Duraisamy, K., Iaccarino, G., Xiao, H., Annu. Rev. Fluid Mech. 51, 357-377 (2019) | 10.1146/annurev-fluid-010518-040547 | DOI verified (Crossref) | Review. Statement on optimal eddy viscosity not located. Unverified for LS. |
| Oberoi, Arias Ramirez, Larsson, APS DFD 2019 abstract "Inferring an effective eddy viscosity from High Fidelity Data" | none (APS meeting abstract) | Exists (abstract partly read); text cut off | Optimisation-based inference of an effective eddy viscosity from resolved data. Content beyond the first sentence unverified. |

Fact: pointwise LS projection of the DNS stress on the strain (our mu_eff) is verified in two papers. Not found: a published regional, strain-squared-weighted LS value (mu_R), or any film-cooling paper using it. Brener et al. 2021 and Liu et al. 2021 are cited by McConkey as applying the approach; not read.

### references.bib (read only; not edited)

Already present (keys): `ling2016b` (Ling, Kurzawski, Templeton, 10.1017/jfm.2016.615), `milani2017` (10.1115/1.4038275), `milani2019`, `milani2020`, `duraisamy2019` (10.1146/annurev-fluid-010518-040547), a Wu, Xiao, Sun, Wang 2019 entry (10.1017/jfm.2019.205; key at line about 450 of the file), `brener2021`, two Sandberg et al. 2018 entries (J. Turbomach. and GT2018-75444).
Not present, suggested:
1. Wu, Xiao, Paterson 2018, Phys. Rev. Fluids 3, 074602, doi 10.1103/PhysRevFluids.3.074602 (the primary citation for the optimal eddy viscosity).
2. McConkey, Yee, Lien 2022, Phys. Fluids 34, 035110, doi 10.1063/5.0083074 (explicit LS closed form in the same sign convention).
3. Schmitt 2007, C. R. Mécanique 335, 617-627, doi 10.1016/j.crme.2007.08.004 (alignment test; supports the cos term).
4. Weatheritt and Sandberg 2016, J. Comput. Phys. 325, 22-37, doi 10.1016/j.jcp.2016.08.015 (only if cited for anisotropy modelling).
Check by the author before use: whether the `.bib` Wu 2019 key and the two Sandberg 2018 entries are the intended ones.

## 4. Table with the recommended estimator

Excluded volume = non-finite on either side (50M buffer or outside the 193M) or k_193 < 0. Points with k_193 = 0 stay in (mu = 0 and weight kept), so they do not appear in the excluded volume; their volume share is listed.
cos <= 0 volume is kept in the LS and listed for information.

| region | mu_R,193 [Pa s] | mu_R,RANS [Pa s] | ln ratio | energy | rho | L_al | S-wt | L_eff (5 pl.) | align (5 pl.) | excluded vol. | cos <= 0 vol. | k=0 vol. | blocks |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| jet 0-3 | 3.497e-4 | 2.030e-4 | -0.544+-0.176 | -0.157+-0.089 | +0.012+-0.020 | -0.286+-0.026 | -0.113+-0.095 | -0.749+-0.111 | +0.478+-0.115 | 0.000 | 0.238 | 0.0027 | 9 |
| film 3-10 | 3.804e-4 | 2.472e-4 | -0.431+-0.436 | -0.320+-0.361 | +0.030+-0.042 | -0.334+-0.132 | +0.192+-0.110 | -0.965+-0.501 | +0.661+-0.359 | 0.000 | 0.221 | 0.0050 | 8 |
| film 10-25 | 2.998e-4 | 2.114e-4 | -0.349+-0.438 | -0.500+-0.265 | +0.005+-0.023 | -0.582+-0.211 | +0.728+-0.073 | -1.301+-0.331 | +0.759+-0.131 | 0.301 | 0.137 | 0.0050 | 6 |
| hole exit | 2.187e-3 | 1.185e-3 | -0.613+-0.099 | -0.138+-0.096 | +0.040+-0.008 | -0.251+-0.054 | -0.264+-0.100 | -0.567+-0.160 | +0.322+-0.116 | 0.000 | 0.016 | 0.0001 | 5 |
| hole, diffuser | 1.211e-3 | 9.100e-4 | -0.286+-0.091 | -0.148+-0.010 | +0.006+-0.001 | -0.002+-0.047 | -0.141+-0.058 | -0.630+-0.091 | +0.627+-0.136 | 0.008 | 0.105 | 0.0029 | 5 |
| hole, metering | 2.353e-5 | 6.127e-5 | +0.957+-0.041 | +0.184+-0.017 | -0.010+-0.000 | +1.040+-0.058 | -0.257+-0.010 | -0.854+-0.012 | +1.888+-0.048 | 0.016 | 0.120 | 0.0055 | 5 |
| approach BL | 1.753e-5 | 2.072e-5 | +0.167+-0.090 | -0.468+-0.017 | -0.000+-0.000 | +0.675+-0.060 | -0.040+-0.023 | -1.369+-0.034 | +2.056+-0.059 | 0.353 | 0.367 | 0.0069 | 6 |

Notes.
- L_eff and align in this table are the 5-player Shapley terms, so L_eff + align differs slightly from the 4-player L_al (see section 2).
- Excluded volume 0.301 (far film) and 0.353 (approach layer) is the 50M buffer / 193M-domain limit, as in `audit.md` Q6.1. Effective ranges x/D 10-20.4 and -10.7 to -3.
- Jackknife: delete-one over the spatial blocks of the point set (n = 5 to 9 blocks per region, unequal sizes). Delete-one flow-through fields (time-batch jackknife of `errorbars.py`) are not in the cache and were not read. SE = sqrt((n-1)/n sum (t_i - mean)^2). Treat as sensitivity to spatial sub-sampling, not as sampling error of the flow statistics.
- Delete-one range of the total (section E of the log): jet [-0.628, -0.437], film 3-10 [-0.907, -0.362], film 10-25 [-0.793, -0.209], hole exit [-0.676, -0.553], diffuser [-0.319, -0.194], metering [+0.920, +0.972], approach [+0.106, +0.219]. Every jet, film, exit and diffuser total stays negative under every block deletion; the metering total stays positive. The approach total is positive in every deletion but is small and the region is 35 % buffer.
- Fact: in film 3-10 the largest block holds 0.33 of the 193M strain weight, and removing it moves the total from -0.431 to -0.907. The film totals depend on a few high-strain blocks; the SE (0.44) is as large as the total.
- Reproduction: LS totals agree with `audit.md` Q4 "LS, all points" (-0.431, -0.349, -0.544, +0.957, -0.286, -0.613, +0.167).
- Confidence: totals high (arithmetic on the cached points); film decomposition medium (order ranges 0.2-0.8); strain-weight term low-medium (depends on the RANS buffer and wall nodes).

## 5. Interpretation register

| Statement | Confidence |
|---|---|
| mu_R is the Frobenius-LS scalar and the exact deviatoric production viscosity of the region | high |
| LS, mean of ln and ln of the volume mean differ strongly in the jet and approach layer; the film sign is the same in all | high |
| The 4-player Shapley split is the least ambiguous exact split; (i) fails because its gap exceeds the total | high (numbers), medium (preference) |
| The L_eff/align separation is poorly determined (order range 0.6-1.2 in the film); L_al is determined to about 0.15 | high |
| The wall-node weight (14-23 % of strain weight) is a property of the estimator and does not change any sign | medium |
| The film total depends on a few high-strain blocks | high |
| Regional LS has no verified literature precedent; the pointwise LS does | medium (search limited to what was fetched) |

## 6. Files

Inspected: `planning/maths_audit/audit.md`, `audit_extract.py`, `q234_regions.py`, `manuscript/references.bib` (read only), scratchpad `pts.npz` (3 141 253 points, 22 blocks, 7 regions), arXiv 1801.02762, 2201.01710 (text in the scratchpad `lit/`), publisher/Crossref pages for the DOIs above. Not read: E10 script, `errorbars.py` beyond the cited lines of the audit, the 193M data.
Scripts (this folder): `estimator.py` (usage `python3 -B estimator.py <scratchpad>/pts.npz A|B|C <scratchpad>/xyz.npz > estimator_mode<M>.log`, about 100 s each; xyz.npz from `extract_xyz.py`), `cos_neg_weight.py`, `extract_xyz.py`, `zcut_weight.py`, `hole_proxy.py`; logs `estimator_modeA.log` (= earlier `estimator.log`), `estimator_modeB.log`, `estimator_modeC.log`, `cos_neg_weight.log`, `zcut_weight.log`, `hole_proxy.log`.
Scratch: `/tmp/claude-1001/-home-sk999-Documents-ASME-2027/f4c92ceb-6df8-4254-96c3-6552d15d0ecf/scratchpad/pts.npz`.

## 7. Next actions

1. Decide between the 4-player split with L_al as the headline (recommended) and the hierarchical split; the difference is up to 0.26 for the weight term.
2. Decide how to treat k_193 = 0 wall nodes in the paper text: keep (estimator as defined) or state the sensitivity (changes the jet total by 0.09, hole exit by 0.16).
3. Correct the plan statement on counter-gradient weight in the jet (wall nodes, not cos < 0) and the jet sign statement (name the estimator).
4. Run the time-batch jackknife with the delete-one fields if sampling SE is needed (37 GB, not read here).
