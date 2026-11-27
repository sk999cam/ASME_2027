"""F14 fig_prt_reference: effective turbulent Prandtl number of the reference on two planes.

Shows Pr_t,eff = mu_eff / (alpha_eff/c_p) of the 193M 3DNS (C193b, stat3 batch b010, 10 FT) on the centre plane
(y = 0, x/D -2..20, z/D 0..3) and on the cross-plane x/D = 10 (y/D -3..3, z/D 0..3).
    mu_eff = -(a:S*)/(2|S*|^2),  alpha_eff = -q.grad T/|grad T|^2  (q = rho c_p u''T'', c_p cancels).
No threshold mask: every plane point is drawn. Blank (white) only where the value is undefined (alpha_eff = 0 or
NaN, i.e. no stress/flux data at the wall row); the count is printed.
Colour bar: inferno, 0 to 2, extend='both': negative values (counter-gradient flux or counter-strain stress) take
the lowest colour and values above 2 (near-zero gradient) the highest; the RANS value 0.9 is marked.
Medians are printed over all points and over the points with |S*| > 0.05 U_inf/D and |grad T| > 0.05 dT0/D
(dT0 = 100.98 K).
Fields are the 10 FT cumulative means of stat3_b010 (mean_time 0.035929 s = 400 D/U), not a snapshot.
Data: /mnt/Elements/DNS/film_cooling/C193b_BR15 stat3_<n>_b010 (via meff_planes_common.py, which holds the
definitions and the cache in /mnt/Elements/DNS/analysis/asme2027_figs/meff_planes/). Complete, not interim.
Prints the Pr_t,eff median in the jet (x/D 0-3), near film (3-10) and far film (10-25 -> to x/D 20 here)
on the centre plane (z/D 0-3) and on the cross-plane.
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

OUT = '/home/sk999/Documents/ASME_2027/figures/fig_prt_reference.png'
OKABE_ITO = ['#000000', '#E69F00', '#56B4E9', '#009E73', '#F0E442', '#0072B2', '#D55E00', '#CC79A7']
VMAX = 2.0

set_mpl_params(0.5, grid_on=False)
plt.rcParams['axes.prop_cycle'] = cycler(color=OKABE_ITO)

R = mc.load('ref')


def prt(p):
    with np.errstate(divide='ignore', invalid='ignore'):
        v = R[f'{p}_meff'] / R[f'{p}_aeff']
    return np.where(np.isfinite(v), v, np.nan)


def above(p):
    return (R[f'{p}_S'] > 0.05 * mc.U_INF / mc.D) & (R[f'{p}_G'] > 0.05 * mc.DT0 / mc.D)


P1, P2 = prt('y0'), prt('x10')
K1, K2 = above('y0'), above('x10')
print(f'undefined points: centre plane {np.isnan(P1).sum()}/{P1.size}, cross-plane {np.isnan(P2).sum()}/{P2.size}')
print(f'out of range 0-2: centre <0 {np.mean(P1[np.isfinite(P1)] < 0):.3f}, >2 {np.mean(P1[np.isfinite(P1)] > 2):.3f}; '
      f'cross <0 {np.mean(P2[np.isfinite(P2)] < 0):.3f}, >2 {np.mean(P2[np.isfinite(P2)] > 2):.3f}')
X1, Z1 = R['y0_A'] / mc.D, R['y0_B'] / mc.D
Y2, Z2 = R['x10_A'] / mc.D, R['x10_B'] / mc.D

fig = plt.figure(figsize=(7.2, 3.45))
W = 7.2
ax1 = fig.add_axes([0.7 / W, 2.45 / 3.45, 6.3 / W, 6.3 * 3 / 22 / 3.45])
ax2 = fig.add_axes([0.7 / W, 0.55 / 3.45, 3.0 / W, 1.5 / 3.45])
cax = fig.add_axes([4.6 / W, 0.75 / 3.45, 2.4 / W, 0.14 / 3.45])
for ax in (ax1, ax2):
    ax.set_facecolor('white')
kw = dict(cmap='inferno', vmin=0.0, vmax=VMAX, shading='auto', rasterized=True)
m = ax1.pcolormesh(X1, Z1, np.ma.masked_invalid(P1), **kw)
ax2.pcolormesh(Y2, Z2, np.ma.masked_invalid(P2), **kw)
# threshold line: inside it |S*| > 0.05 U/D and |grad T| > 0.05 dT0/D (where the ratio is well defined); nothing hidden
ckw = dict(levels=[0.5], colors='#56B4E9', linewidths=0.9, linestyles='--')
ax1.contour(X1, Z1, K1.astype(float), **ckw)
ax2.contour(Y2, Z2, K2.astype(float), **ckw)
ax1.set_aspect('equal')
ax2.set_aspect('equal')
ax1.set_xlim(-2, 20)
ax1.set_ylim(0, 3)
ax2.set_xlim(3, -3)                      # y/D increases to the left: plane seen looking downstream
ax2.set_xticks([3, 0, -3])
ax2.set_ylim(0, 3)
ax1.set_xlabel(r'$x/D$')
ax1.set_ylabel(r'$z/D$')
ax2.set_xlabel(r'$y/D$')
ax2.set_ylabel(r'$z/D$')
ax1.set_yticks([0, 1, 2, 3])
ax2.set_yticks([0, 1, 2, 3])
LBL = dict(color='black', va='top', bbox=dict(fc='white', ec='none', pad=1.5))
ax1.text(0.012, 0.88, r'(a) $y/D=0$', transform=ax1.transAxes, **LBL)
ax2.text(0.02, 0.95, r'(b) $x/D=10$', transform=ax2.transAxes, **LBL)
cb = fig.colorbar(m, cax=cax, orientation='horizontal', extend='both')
cb.set_ticks([0, 0.5, 0.9, 1.5, 2.0])
cb.set_ticklabels(['0', '0.5', '0.9', '1.5', '2'])
cb.ax.axvline(0.9, color='white', lw=1.4)
cb.set_label(r'$\mathrm{Pr}_{t,\mathrm{eff}}$')
fig.savefig(OUT, dpi=300)
print('saved', OUT)

# representative values (z/D 0-3)
def report(lab, P, K, s):
    for nm, sel in (('all points', s), ('above thresholds', s & K)):
        v = P[sel]
        v = v[np.isfinite(v)]
        print(f'{lab} [{nm}]: median {np.median(v):.2f}, IQR {np.percentile(v, 25):.2f}-{np.percentile(v, 75):.2f}, '
              f'n {v.size}/{s.sum()}')


for lab, lo, hi in [('jet x/D 0-3', 0, 3), ('near film 3-10', 3, 10), ('far film 10-20', 10, 20)]:
    report(f'centre plane {lab}', P1, K1, (X1 >= lo) & (X1 < hi))
report('cross-plane x/D=10', P2, K2, np.ones(P2.shape, bool))
