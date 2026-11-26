"""F16 fig_prt_regions: classic regional turbulent Prandtl number of the 193M reference, against the model values.

Shows Pr_t,R = [sum(-rho u''w'' dU/dz) / sum(dU/dz)^2] / [sum(-rho w''T'' dT/dz) / sum(dT/dz)^2] (Favre means, wall-normal
classic form of Kohli and Bogard 2005, least squares over the points of a region, no threshold) of the 193M 3DNS
(C193b) for z/D >= 0.3: jet (x/D 0-3), near film (3-10) and far film (10-20) on the centre plane y = 0, and the
cross-plane x/D = 10. Stat3 batch b010 (10 FT cumulative mean): filled markers, bars = 95 % block-bootstrap interval
(1 D wide x bins on the centre plane, 0.5 D wide y bins on the cross-plane, 1000 draws). Batch b006 (6 FT): hollow.
The weighting by the squared gradients means a point with a near-zero gradient cannot produce a near-zero
denominator; z/D < 0.3 is excluded because the viscous layer dominates the sums there.
Horizontal lines: 0.9 RANS (solver default); 0.85 log layer (Kays 1994); 0.5 free shear layers (Schlichting and
Gersten 2017; Rochhausen et al. 2015); Kohli and Bogard (2005, J. Heat Transfer 127:513) about 1 at high and about 2 at
low freestream turbulence (the 193M has Tu = 3.6 % at x/D -5, between their 0.6 % and 20 %).
Data: /mnt/Elements/DNS/film_cooling/C193b_BR15, stat3_b006 and stat3_b010, through the plane caches of
meff_planes_common.py and prt_corr_common.py (cache /mnt/Elements/DNS/analysis/asme2027_figs/prt_classes/).
Complete, not interim.
"""
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
OKABE_ITO = ['#000000', '#E69F00', '#56B4E9', '#009E73', '#F0E442', '#0072B2', '#D55E00', '#CC79A7']

rng = np.random.default_rng(3)
R = {b: dict(np.load(f'{mc.CACHE}/ref193_b{b:03d}.npz')) for b in (6, 10)}
Y = {b: fc.prep(R[b], fc.load_corr('full', b, R[b]['y0_A'].shape), 'y0') for b in R}
X = {b: fc.prep(R[b], fc.load_corr('x10', b, R[b]['x10_A'].shape), 'x10') for b in R}

classes = []
for name, x0, x1 in (('jet\n$x/D$ 0-3', 0.0, 3.0), ('near film\n$x/D$ 3-10', 3.0, 10.0), ('far film\n$x/D$ 10-20', 10.0, 20.0)):
    m = lambda d, x0=x0, x1=x1: (d['x'] >= x0) & (d['x'] < x1) & (d['z'] >= 0.3)
    classes.append((name, fc.cls(Y[6], Y[10], m(Y[6]), m(Y[10]), rng), '#D55E00' if x0 == 0 else '#0072B2' if x0 == 3 else '#009E73'))
mx = lambda d: d['z'] >= 0.3
classes.append(('cross-plane\n$x/D=10$', fc.cls(X[6], X[10], mx(X[6]), mx(X[10]), rng), '#000000'))
for name, r, _ in classes:
    print(name.replace('\n', ' '), f'b006 {r["pc6"]:.3f}  b010 {r["pc10"]:.3f}  95% [{r["ci"][0]:.3f}, {r["ci"][1]:.3f}]  n {r["n"]}')

set_mpl_params(0.5, grid_on=False)
plt.rcParams['axes.prop_cycle'] = cycler(color=OKABE_ITO)
fig, ax = plt.subplots(figsize=(3.6, 3.4))
refs = [(0.9, ':', '#000000', 'RANS 0.9 (solver default)'), (0.85, '--', '#CC79A7', 'log layer 0.85 (Kays 1994)'),
        (0.5, '-.', '#56B4E9', 'free shear 0.5'), (1.0, '--', '0.5', r'K&B high FST, $\approx1$'),
        (2.0, '-.', '0.5', r'K&B low FST, $\approx2$')]
for v, ls, c, lab in refs:
    ax.axhline(v, ls=ls, color=c, zorder=1)
for i, (name, r, c) in enumerate(classes):
    ax.errorbar(i, r['pc10'], yerr=[[max(r['pc10'] - r['ci'][0], 0)], [max(r['ci'][1] - r['pc10'], 0)]], color=c, marker='o',
                capsize=3, zorder=3)
    ax.plot(i + 0.12, r['pc6'], color=c, marker='o', mfc='none', ls='none', zorder=3)
ax.set_xticks(range(len(classes)))
ax.set_xticklabels([n for n, _, _ in classes])
ax.set_xlim(-0.5, len(classes) - 0.5)
ax.set_ylim(0, 2.4)
ax.set_ylabel(r'$\mathrm{Pr}_{t,R}$ ($z/D\geq0.3$)')
h = [Line2D([], [], color=c, ls=ls, label=lab) for v, ls, c, lab in refs]
h += [Line2D([], [], color='k', ls='none', marker='o', label='193M b010 (10 FT), 95 %'),
      Line2D([], [], color='k', ls='none', marker='o', mfc='none', label='193M b006 (6 FT)')]
fig.legend(handles=h, loc='lower center', ncol=1, frameon=False, bbox_to_anchor=(0.5, -0.01))
fig.tight_layout(rect=(0, 0.33, 1, 1))
fig.savefig(OUT)
print('saved', OUT)
