"""F16 fig_prt_regions: classic regional turbulent Prandtl number of the 193M reference.

Pr_t,R = [sum(-rho u''w'' dU/dz) / sum(dU/dz)^2] / [sum(-rho w''T'' dT/dz) / sum(dT/dz)^2] (Favre means, wall-normal
classic form of Kohli and Bogard 2005, least squares over the points of a region, z/D >= 0.3, no other threshold) of the
193M 3DNS (C193b), batch b010 (10 FT), for jet (x/D 0-3), near film (3-10) and far film (10-20).
Filled markers: VOLUME values over all points of each region (paper values), from the full-record sums of
/mnt/Elements/DNS/analysis/asme2027_figs/n8_lsq/prt_classic_c193b_b010_jk.npz.
Error bars: +-2 SE of the delete-one jackknife over the 10 flow-through records, file figures/code/data/prt_volume_se.txt
(columns: region, lower, upper; regions near_jet_0_3, wake_3_10, wake_10_25 = x/D 10-20).
Open markers: centre-plane (y = 0) b010 values with 95 % block-bootstrap bars (1 D wide x bins, 1000 draws); the
far-film centre-plane value covers x/D 10-20 only (stated in the caption).
Horizontal reference lines, labelled in the plot: 0.9 RANS (solver default); 0.85 log layer (Kays 1994); 0.5 free
shear layers; Kohli and Bogard (2005) about 1 at high and about 2 at low freestream turbulence.
Centre-plane data: plane caches of meff_planes_common.py and prt_corr_common.py (read only, /mnt/Elements).
Complete, not interim.
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
import fig_prt_classes as fc

OUT = '/home/sk999/Documents/ASME_2027/figures/fig_prt_regions.pdf'
VOL = '/home/sk999/Documents/4dns_filmcooled/analysis/eps_prt_c193b/results/prt_classic.txt'
SE = '/home/sk999/Documents/ASME_2027/figures/code/data/prt_volume_se.txt'
OKABE_ITO = ['#000000', '#E69F00', '#56B4E9', '#009E73', '#F0E442', '#0072B2', '#D55E00', '#CC79A7']

# volume values, batch 10 only
vol = {}
for ln in open(VOL):
    p = ln.split()
    if not p or p[0].startswith('#') or p[0] == 'region':
        continue
    if int(p[1]) == 10:
        vol[p[0]] = float(p[4])
# paper values: full-record regional sums of the jackknife file (far film x/D 10-20); replaces the text-file values
_J = np.load('/mnt/Elements/DNS/analysis/asme2027_figs/n8_lsq/prt_classic_c193b_b010_jk.npz', allow_pickle=True)
_s = _J['sums'][0]
for _k, _v in zip(('near_jet_0_3', 'wake_3_10', 'wake_10_25'), (_s[:, 0] / _s[:, 1]) / (_s[:, 2] / _s[:, 3])):
    vol[_k] = float(_v)
vol_se = {}
if os.path.exists(SE):
    for ln in open(SE):
        p = ln.replace(',', ' ').split()
        if len(p) >= 3 and not p[0].startswith('#'):
            vol_se[p[0]] = (float(p[1]), float(p[2]))

# centre-plane b010 values with 95 % intervals
rng = np.random.default_rng(3)
R = {b: dict(np.load(f'{mc.CACHE}/ref193_b{b:03d}.npz')) for b in (6, 10)}
Y = {b: fc.prep(R[b], fc.load_corr('full', b, R[b]['y0_A'].shape), 'y0') for b in R}

classes = []
for name, key, x0, x1, c in (('jet\n$x/D$ 0-3', 'near_jet_0_3', 0.0, 3.0, '#D55E00'),
                             ('near film\n$x/D$ 3-10', 'wake_3_10', 3.0, 10.0, '#0072B2'),
                             ('far film\n$x/D$ 10-20', 'wake_10_25', 10.0, 20.0, '#009E73')):
    m = lambda d, x0=x0, x1=x1: (d['x'] >= x0) & (d['x'] < x1) & (d['z'] >= 0.3)
    classes.append((name, key, fc.cls(Y[6], Y[10], m(Y[6]), m(Y[10]), rng), c))
for name, key, r, _ in classes:
    print(name.replace('\n', ' '), f'volume b010 {vol[key]:.3f}  centre-plane b010 {r["pc10"]:.3f}  '
          f'95% [{r["ci"][0]:.3f}, {r["ci"][1]:.3f}]')

# Drawn at the printed size for a full column (3.56 in; include width=\columnwidth, scale 1); body text 10 pt.
# The reference values are named in the legend below the axes (the margin labels of the old version cannot hold 10 pt text).
set_mpl_params(10 / 22, grid_on=False)
plt.rcParams['axes.prop_cycle'] = cycler(color=OKABE_ITO)
fig, ax = plt.subplots(figsize=(3.56, 4.6))
XR = 3.0
refs = [(0.9, ':', '#000000', 'RANS, 0.9', 1.12), (0.85, '--', '#CC79A7', 'log layer, 0.85 (Kays 1994)', 0.62),
        (0.5, '-.', '#56B4E9', 'free shear, 0.5', 0.25), (1.0, '--', '0.5', r'K&B high FST, $\approx1$', 1.6),
        (2.0, '-.', '0.5', r'K&B low FST, $\approx2$', 2.0)]
for v, ls, c, lab, ty in refs:
    ax.axhline(v, ls=ls, color=c, zorder=1, label=lab)
for i, (name, key, r, c) in enumerate(classes):
    kw = {}
    if key in vol_se:
        lo, hi = vol_se[key]
        kw = dict(yerr=[[vol[key] - lo], [hi - vol[key]]], capsize=3)
    ax.errorbar(i - 0.1, vol[key], color=c, marker='o', ls='none', zorder=3, **kw)
    ax.errorbar(i + 0.1, r['pc10'], yerr=[[max(r['pc10'] - r['ci'][0], 0)], [max(r['ci'][1] - r['pc10'], 0)]], color=c,
                marker='o', mfc='none', ls='none', capsize=3, zorder=3)
ax.set_xticks(range(len(classes)))
ax.set_xticklabels([n for n, _, _, _ in classes])
ax.set_xlim(-0.5, XR - 0.5)
ax.set_ylim(0, 2.2)
ax.set_ylabel(r'$\mathrm{Pr}_{t,R}$ ($z/D\geq0.3$)')
h = [Line2D([], [], color='k', ls='none', marker='o', label=r'volume, $\pm2$ SE'),
     Line2D([], [], color='k', ls='none', marker='o', mfc='none', label='centre plane, 95 %')]
h = h + [Line2D([], [], color=c, ls=ls, label=lab) for v, ls, c, lab, ty in refs]
ax.legend(handles=h, loc='upper center', bbox_to_anchor=(0.5, -0.2), frameon=False, ncol=1,
          fontsize=plt.rcParams['font.size'], handlelength=2.2)
fig.tight_layout(pad=0.3)
fig.savefig(OUT)
print('saved', OUT)
