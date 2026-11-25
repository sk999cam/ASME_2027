"""F18 fig_centre_factors: where the length scale and the alignment differ, on the centre plane.

Centre plane y = 0, x/D -2..20, z/D 0..3.
 (a) ln(L_eff,RANS / L_eff,ref): RANS 50M (L = mu_t/(rho k^1/2)) over the 193M 3DNS (L_eff = k^1/2 |b|/|S*|).
 (b) cos(theta) of the 193M 3DNS reference, cos(theta) = -(a:S*)/(|a||S*|) (the RANS value is 1 by construction).
bwr with limits symmetric about zero: (a) +-2, (b) [-1, 1]. No threshold mask: every plane point is drawn, including
the free stream with near-zero strain. Panel (a) colour bar has extend='both' (values beyond +-2 take the end
colours). Blank (white) only where the value is undefined: k <= 0 or L <= 0 in either run (a), k <= 0 or S* = 0 in
the reference (b); the counts are printed.
Data: 193M 3DNS C193b stat3 batch b010 (/mnt/Elements/DNS/film_cooling/C193b_BR15); RANS 50M
rans_50M_Lm25b_BR15_FST_ZPG_p0dns (iteration 122 686, plenum p0 81 110 Pa, BR 1.592 against 1.53 of the 3DNS)
from /mnt/Tank/DNS/analysis/meff_ransp0_c50b, on the C50b mesh; both sampled on the same plane points by trilinear
interpolation (meff_planes_common.py). Complete.
Prints area medians of ln ratio and cos theta in the jet (x/D 0-3), near film (3-10) and far film (10-20).
"""
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from cycler import cycler

sys.path.insert(0, '/home/sk999/Documents/3DNS_FP')
sys.path.insert(0, '/home/sk999/Documents/ASME_2027/figures/code')
from functions import set_mpl_params
import meff_planes_common as mc

OUT = '/home/sk999/Documents/ASME_2027/figures/fig_centre_factors.png'
OKABE_ITO = ['#000000', '#E69F00', '#56B4E9', '#009E73', '#F0E442', '#0072B2', '#D55E00', '#CC79A7']
LMAX = 2.0

set_mpl_params(0.5, grid_on=False)
plt.rcParams['axes.prop_cycle'] = cycler(color=OKABE_ITO)

ref, rans = mc.load('ref'), mc.load('rans')
X, Z = ref['y0_A'] / mc.D, ref['y0_B'] / mc.D
ok = (ref['y0_k'] > 0) & (rans['y0_k'] > 0) & (ref['y0_L'] > 0) & (rans['y0_L'] > 0)
okc = (ref['y0_k'] > 0) & (ref['y0_S'] > 0)
with np.errstate(divide='ignore', invalid='ignore'):
    lr = np.where(ok, np.log(rans['y0_L'].astype(float) / ref['y0_L'].astype(float)), np.nan)
cs = np.where(okc, ref['y0_cos'].astype(float), np.nan)
print(f'undefined points (k <= 0, L <= 0 or S = 0): panel (a) {np.isnan(lr).sum()}/{lr.size}, '
      f'panel (b) {np.isnan(cs).sum()}/{cs.size}')
print(f'out of range: (a) |ln ratio| > 2: {np.mean(np.abs(lr[ok]) > 2):.3f}; (b) cos < 0: {np.mean(cs[okc] < 0):.3f}')

W, H = 7.7, 3.1
pw, ph = 6.0, 6.0 * 3 / 22
fig = plt.figure(figsize=(W, H))
a1 = fig.add_axes([0.7 / W, 1.95 / H, pw / W, ph / H])
a2 = fig.add_axes([0.7 / W, 0.62 / H, pw / W, ph / H])
c1 = fig.add_axes([(0.7 + pw + 0.15) / W, 1.95 / H, 0.14 / W, ph / H])
c2 = fig.add_axes([(0.7 + pw + 0.15) / W, 0.62 / H, 0.14 / W, ph / H])
m1 = a1.pcolormesh(X, Z, np.ma.masked_invalid(lr), cmap='bwr', vmin=-LMAX, vmax=LMAX, shading='auto',
                   rasterized=True)
m2 = a2.pcolormesh(X, Z, np.ma.masked_invalid(cs), cmap='bwr', vmin=-1, vmax=1, shading='auto', rasterized=True)
for ax, t in ((a1, r'(a) $\ln(L_\mathrm{eff,RANS}/L_\mathrm{eff,ref})$'), (a2, r'(b) $\cos\theta$, reference')):
    ax.set_facecolor('white')
    ax.set_aspect('equal')
    ax.set_xlim(-2, 20)
    ax.set_ylim(0, 3)
    ax.set_yticks([0, 1, 2, 3])
    ax.set_ylabel(r'$z/D$')
    ax.text(0.012, 0.92, t, transform=ax.transAxes, va='top')
a1.set_xticklabels([])
a2.set_xlabel(r'$x/D$')
cb1 = fig.colorbar(m1, cax=c1, extend='both')
cb1.set_label(r'$\ln(L_\mathrm{RANS}/L_\mathrm{ref})$', fontsize=9)
cb2 = fig.colorbar(m2, cax=c2)
cb2.set_label(r'$\cos\theta$', fontsize=9)
fig.savefig(OUT, dpi=300)
print('saved', OUT)

for lab, lo, hi in [('jet x/D 0-3', 0, 3), ('near film 3-10', 3, 10), ('far film 10-20', 10, 20)]:
    s = (X >= lo) & (X < hi)
    a, b = lr[s], cs[s]
    a, b = a[np.isfinite(a)], b[np.isfinite(b)]
    print('[all points] ', end='')
    print(f'{lab}: ln ratio median {np.median(a):+.2f} (IQR {np.percentile(a, 25):+.2f}..{np.percentile(a, 75):+.2f}); '
          f'cos median {np.median(b):+.2f} (IQR {np.percentile(b, 25):+.2f}..{np.percentile(b, 75):+.2f}); '
          f'valid {a.size}/{s.sum()}')
