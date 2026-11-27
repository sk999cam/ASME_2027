"""F20 fig_prt_data: does the reference match the literature, and what is new: spanwise (lateral) diffusion.

(a) Directional ratios of the eddy diffusivities of the 193M 3DNS (C193b), regional least squares, no threshold:
      nu_y/nu_z    = [sum(-rho u''v'' dU/dy) / sum(dU/dy)^2] / [sum(-rho u''w'' dU/dz) / sum(dU/dz)^2]   (momentum)
      alpha_y/alpha_z = [sum(-rho v''T'' dT/dy) / sum(dT/dy)^2] / [sum(-rho w''T'' dT/dz) / sum(dT/dz)^2]  (heat)
    on the cross-planes x/D = 10 and 20 (x/D = 4 needs a 20 GB read, not done), for the classes: all points with
    z/D >= 0.3; plume centre |y|/D < 1 and z/D >= 0.3; flanks |y|/D >= 1.5 and z/D >= 0.3; wall layer z/D < 0.3.
    Markers: batch b010 (10 FT), bars = 95 % block-bootstrap interval (0.5 D wide y bins, 1000 draws).
    Grey band: Kaszeta and Simon (2000) lateral/normal diffusivity ratio 1.5-4 (band as given in the
    brief); line at 1: isotropic eddy diffusivity of the RANS.
(b) Lateral heat flux over the film edges, RANS 50M (p0dns, iteration 122 686) against the 193M, on the same plane
    points: edge = |y|/D >= 1 and z/D >= 0.3. Measure P_y = sum(-F_y dT/dy), F_y = lateral turbulent heat flux / c_p,
    i.e. alpha_y,C times the sum of dT/dy squared; RANS F_y = -(mu_t/0.9) dT/dy. Exact identity:
      ln(P_R/P_D) = ln(mu_t,R/mu_z,D) - ln(nu_y/nu_z) + ln(Pr_y,D/0.9) + ln(sum g_R^2 / sum g_D^2),
    with mu_t,R the dT/dy-squared weighted mean of the RANS mu_t, mu_z,D = sum(-rho u'w' dU/dz)/sum(dU/dz)^2 of the
    193M (momentum, wall-normal), nu_y/nu_z the lateral to wall-normal momentum ratio, Pr_y,D = (nu_y)/(alpha_y) the
    193M lateral classic Pr_t, g = dT/dy. The four terms add to the total by construction; the script checks it.
    The RANS dT/dy is the finite difference of the plane temperature (spacing 0.01 D); the 193M dT/dy is the
    gradient of the mean field. Stack: positive terms up, negative down; marker = total (193M at b010).
    The mu_t term compares the RANS mu_t with the classic wall-normal momentum viscosity of the 193M. Triangle: the
    total without the gradient term, ln(alpha_y,R/alpha_y,D), the ratio of the lateral heat diffusivities.
Data: /mnt/Elements/DNS/film_cooling/C193b_BR15, stat3_b010 through meff_planes_common.py and
prt_corr_common.py (cache /mnt/Elements/DNS/analysis/asme2027_figs/prt_classes/); RANS fields from
meff_planes_common.py (rans50_p0dns.npz, read-only). Complete for x/D = 10 and 20 when both caches exist.
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from cycler import cycler
from matplotlib.lines import Line2D

sys.path.insert(0, '/home/sk999/Documents/3DNS_FP')
sys.path.insert(0, '/home/sk999/Documents/ASME_2027/figures/code')
from functions import set_mpl_params
import meff_planes_common as mc
import prt_corr_common as cc
import fig_prt_classes as fc

OUT = '/home/sk999/Documents/ASME_2027/figures/fig_prt_data.pdf'
OKABE_ITO = ['#000000', '#E69F00', '#56B4E9', '#009E73', '#F0E442', '#0072B2', '#D55E00', '#CC79A7']
D = mc.D
NBOOT = 1000
# row indices in fc.KEYS
CNM, CDM, CNA, CDA, LNM, LDM, LNA, LDA = 5, 6, 7, 8, 9, 10, 11, 12
NUR = (LNM, LDM, CNM, CDM)          # nu_y/nu_z
ALR = (LNA, LDA, CNA, CDA)          # alpha_y/alpha_z
rng = np.random.default_rng(5)


def ratio(S, spec):
    a, b, c, d = spec
    with np.errstate(divide='ignore', invalid='ignore'):
        return (S[a].sum() / S[b].sum()) / (S[c].sum() / S[d].sum())


def boot(S, spec):
    have = np.nonzero(S[CDM] > 0)[0]
    if len(have) < 3:
        return np.nan, np.nan
    Sh = S[:, have]
    w = rng.multinomial(len(have), np.ones(len(have)) / len(have), size=NBOOT).astype(float)
    a, b, c, d = spec
    with np.errstate(divide='ignore', invalid='ignore'):
        v = ((w @ Sh[a]) / (w @ Sh[b])) / ((w @ Sh[c]) / (w @ Sh[d]))
    v = v[np.isfinite(v)]
    return tuple(np.percentile(v, [2.5, 97.5])) if len(v) > NBOOT // 2 else (np.nan, np.nan)


CLASSES = [('all, $z/D\\geq0.3$', lambda d: d['z'] >= 0.3),
           ('centre, $|y|/D<1$', lambda d: (d['z'] >= 0.3) & (np.abs(d['x']) < 1)),
           ('flanks, $|y|/D\\geq1.5$', lambda d: (d['z'] >= 0.3) & (np.abs(d['x']) >= 1.5)),
           ('wall, $z/D<0.3$', lambda d: d['z'] < 0.3)]
EDGE = lambda d: (d['z'] >= 0.3) & (np.abs(d['x']) >= 1.0)

log = fc.log
R = {b: dict(np.load(f'{mc.CACHE}/ref193_b{b:03d}.npz')) for b in (10,)}      # 193M at b010 only (owner rule 2026-11-26)
RANS = dict(np.load(f'{mc.CACHE}/rans50_p0dns.npz'))
planes = [p for p in ('x10', 'x20') if all(os.path.exists(os.path.join(cc.CORR_DIR, f'corr_{p}_b{b:03d}.npz')) for b in (10,))]
log(f'planes available: {planes}; x4 not read (about 20 GB per batch)')

# ---- (a) directional ratios
rows = []          # (plane, class label, dict)
for p in planes:
    P = {b: fc.prep(R[b], fc.load_corr(p, b, R[b][f'{p}_A'].shape), p) for b in R}
    for lab, f in CLASSES:
        S10, n = fc.sums(P[10], f(P[10]))
        r = dict(n=n, nu10=ratio(S10, NUR), nuci=boot(S10, NUR),
                 al10=ratio(S10, ALR), alci=boot(S10, ALR),
                 pry=ratio(S10, (LNM, LDM, LNA, LDA)), prz=ratio(S10, (CNM, CDM, CNA, CDA)))
        rows.append((p, lab, r))
        log(f'{p} {lab.replace(chr(36), ""):28s} n={n:6d} nu_y/nu_z b010 {r["nu10"]:7.2f} [{r["nuci"][0]:6.2f},{r["nuci"][1]:6.2f}] | '
            f'alpha_y/alpha_z b010 {r["al10"]:7.2f} [{r["alci"][0]:6.2f},{r["alci"][1]:6.2f}] | Pr_y {r["pry"]:.2f} Pr_z {r["prz"]:.2f}')

# ---- (b) lateral flux split, RANS against 193M
split = {}
for p in planes:
    Pd = {b: fc.prep(R[b], fc.load_corr(p, b, R[b][f'{p}_A'].shape), p) for b in R}
    T = RANS[f'{p}_T']
    dy = (RANS[f'{p}_A'][0, 1] - RANS[f'{p}_A'][0, 0])
    gR = (np.gradient(T, dy, axis=1)).ravel()
    mut = RANS[f'{p}_meff'].ravel()
    for b in (10,):
        d = Pd[b]
        m = EDGE(d) & d['ok']
        S, n = fc.sums(d, EDGE(d))
        Pdat = S[LNA].sum()                           # sum(-rho v'T' dT/dy)
        gD2 = S[LDA].sum()
        PR = np.sum((mut[m] / 0.9) * gR[m] ** 2)
        gR2 = np.sum(gR[m] ** 2)
        mu_tR = np.sum(mut[m] * gR[m] ** 2) / gR2
        mu_zD = S[CNM].sum() / S[CDM].sum()
        Rnu = ratio(S, NUR)
        pry = ratio(S, (LNM, LDM, LNA, LDA))
        terms = dict(mut=np.log(mu_tR / mu_zD), aniso=-np.log(Rnu), pr=np.log(pry / 0.9), grad=np.log(gR2 / gD2))
        tot = np.log(PR / Pdat)
        split[(p, b)] = dict(terms=terms, tot=tot)
        log(f'split {p} b{b:03d}: ln(P_R/P_D) {tot:+.3f} = mu_t deficit {terms["mut"]:+.3f} + anisotropy {terms["aniso"]:+.3f} + '
            f'Pr_t setting {terms["pr"]:+.3f} + gradient {terms["grad"]:+.3f}; sum {sum(terms.values()):+.3f}; '
            f'mu_t,R {mu_tR:.3e}, mu_z,D {mu_zD:.3e}, nu_y/nu_z {Rnu:.2f}, Pr_y,D {pry:.2f}, P_R {PR:.3e}, P_D {Pdat:.3e}')
        assert abs(tot - sum(terms.values())) < 1e-9
with open(os.path.join(cc.CORR_DIR, 'prt_data_tables.txt'), 'w') as fh:
    fh.write('\n'.join(fc.LINES) + '\n')

# ---- figure
# Drawn at the printed size for full text width (7.5 in; figure* with width=\textwidth, scale 1); body text 10 pt.
set_mpl_params(10 / 22, grid_on=False)
plt.rcParams['axes.prop_cycle'] = cycler(color=OKABE_ITO)
fig, (axa, axb) = plt.subplots(1, 2, figsize=(7.5, 4.3), gridspec_kw=dict(width_ratios=[2.6, 1]))
axa.axhspan(1.5, 4.0, color='0.85', zorder=0)
axa.axhline(1.0, color='#000000', ls=':', zorder=1)
x = 0
ticks, tl = [], []
for p in planes:
    for lab, _ in CLASSES:
        r = [q for pp, ll, q in rows if pp == p and ll == lab][0]
        for dx, k, c, mk in ((-0.12, 'nu', '#D55E00', 'o'), (0.12, 'al', '#0072B2', 's')):
            lo, hi = r[f'{k}ci']
            axa.errorbar(x + dx, r[f'{k}10'], yerr=[[max(r[f'{k}10'] - lo, 0)], [max(hi - r[f'{k}10'], 0)]], color=c, marker=mk,
                         capsize=2, zorder=3)
        ticks.append(x)
        tl.append(lab.replace(', ', '\n'))
        x += 1
    axa.text(x - 2.5, -0.19, f'$x/D={p[1:]}$', transform=axa.get_xaxis_transform(), ha='center', va='top')
    x += 0.5
axa.set_yscale('log')
axa.set_ylim(0.3, 40)
axa.set_xticks(ticks)
axa.set_xticklabels(tl, fontsize=plt.rcParams['font.size'] - 1)
axa.set_ylabel(r'lateral / wall-normal')
axa.text(0.02, 0.95, '(a)', transform=axa.transAxes, va='top')
# (b)
cols = [('mut', r'$\ln(\mu_{t,\mathrm{R}}/\mu_{z,\mathrm{D}})$', '#0072B2', '/'), ('aniso', r'anisotropy', '#D55E00', '\\'),
        ('pr', r'$\mathrm{Pr}_t$ setting', '#009E73', 'x'), ('grad', r'gradient', '#CC79A7', '.')]
for i, p in enumerate(planes):
    pos = neg = 0.0
    for k, lab, c, h in cols:
        v = split[(p, 10)]['terms'][k]
        if v >= 0:
            axb.bar(i, v, bottom=pos, color=c, hatch=h, edgecolor='white', width=0.55, zorder=2)
            pos += v
        else:
            axb.bar(i, v, bottom=neg, color=c, hatch=h, edgecolor='white', width=0.55, zorder=2)
            neg += v
    axb.plot(i, split[(p, 10)]['tot'], color='k', marker='D', ls='none', zorder=4)
    axb.plot(i - 0.2, split[(p, 10)]['tot'] - split[(p, 10)]['terms']['grad'], color='k', marker='v', ls='none', zorder=4)
axb.axhline(0, color='k', lw=0.6)
axb.set_xticks(range(len(planes)))
axb.set_xticklabels([f'$x/D={p[1:]}$' for p in planes])
axb.set_xlim(-0.6, len(planes) - 0.4)
axb.set_ylabel(r'$\ln(P_{y,\mathrm{RANS}}/P_{y,\mathrm{193M}})$')
axb.text(0.04, 0.95, '(b)', transform=axb.transAxes, va='top')
h = [Line2D([], [], color='#D55E00', marker='o', ls='none', label=r'$\nu_y/\nu_z$ (momentum)'),
     Line2D([], [], color='#0072B2', marker='s', ls='none', label=r'$\alpha_y/\alpha_z$ (heat)'),
     Line2D([], [], color='k', marker='o', ls='none', label='193M, 95 % interval'),
     Line2D([], [], color='k', ls=':', label='RANS, isotropic 1'),
     Line2D([], [], color='0.85', lw=6, label='Kaszeta and Simon 1.5-4')]
h += [plt.Rectangle((0, 0), 1, 1, fc=c, hatch=hh, ec='white', label=lab) for k, lab, c, hh in cols]
h += [Line2D([], [], color='k', marker='D', ls='none', label='total'),
      Line2D([], [], color='k', marker='v', ls='none', label=r'total without gradient ($\ln\alpha_{y,R}/\alpha_{y,D}$)')]
fig.legend(handles=h, loc='lower center', ncol=4, frameon=False, bbox_to_anchor=(0.5, -0.01), fontsize=plt.rcParams['font.size'] - 1,
           handlelength=1.4, columnspacing=0.7)
fig.tight_layout(pad=0.3, rect=(0, 0.19, 1, 1))
fig.savefig(OUT)
log(f'saved {OUT}')
