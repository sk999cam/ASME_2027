"""F17 fig_cross_planes: turbulence kinetic energy and effective length scale on three cross-planes.

y-z planes at x/D = 4, 10, 20 (y/D -3..3 over one pitch, z/D 0..3). Rows from the top: k/U_inf^2 of the RANS 50M and
of the 193M 3DNS, then L_eff/D of the RANS 50M and of the 193M 3DNS. White lines: contours of the mean
temperature scalar eta = (T0_fs - T)/(T0_fs - T0_cool) = 0.1, 0.2, 0.3, 0.4 (from the mean static temperature of
the same run). inferno; one colour scale per quantity.
    3DNS: L_eff = k^1/2 |b|/|S*|,  |b| = |a|/(2 rho k).  RANS: L = mu_t/(rho k^1/2).
No masking: every plane point is drawn. k: linear scale, limit the 99.5th percentile over the points with
|S*| > 0.05 U_inf/D. L: log scale, L_eff/D 0.01-1, so the large L of the weakly strained freestream is shown
without saturating (values outside 0.01-1 take the end colours).
Light grey only where k <= 0 or L <= 0 (no value defined).
Data: 193M 3DNS C193b stat3 batch b010 (/mnt/Elements/DNS/film_cooling/C193b_BR15); RANS 50M
rans_50M_Lm25b_BR15_FST_ZPG_p0dns (iteration 122 686) from /mnt/Tank/DNS/analysis/meff_ransp0_c50b on the C50b
mesh. Both are sampled on the same plane points by trilinear interpolation (meff_planes_common.py). Complete.
Prints the planewise median and 99th percentile of k and L in the region z/D 0-3.
"""
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize
from cycler import cycler

sys.path.insert(0, '/home/sk999/Documents/3DNS_FP')
sys.path.insert(0, '/home/sk999/Documents/ASME_2027/figures/code')
from functions import set_mpl_params
import meff_planes_common as mc

OUT = '/home/sk999/Documents/ASME_2027/figures/fig_cross_planes.png'
OKABE_ITO = ['#000000', '#E69F00', '#56B4E9', '#009E73', '#F0E442', '#0072B2', '#D55E00', '#CC79A7']
XS = [4, 10, 20]
ETA_LEVELS = [0.1, 0.2, 0.3, 0.4]
L_LOG = (0.01, 1.0)   # log colour limits of L_eff/D (193M 1st-99.5th percentile about 0.006-0.97)

set_mpl_params(0.5, grid_on=False)
plt.rcParams['axes.prop_cycle'] = cycler(color=OKABE_ITO)

DAT = {'RANS 50M': mc.load('rans'), '193M 3DNS': mc.load('ref')}
rows = [('k', 'RANS 50M'), ('k', '193M 3DNS'), ('L', 'RANS 50M'), ('L', '193M 3DNS')]


def field(q, rung, p, strained_only=False):
    d = DAT[rung]
    v = d[f'{p}_{q}'].astype(float)
    if q == 'k':
        return v / mc.U_INF ** 2
    ok = (d[f'{p}_k'] > 0) & (v > 0)
    if strained_only:
        ok &= d[f'{p}_S'] > 0.05 * mc.U_INF / mc.D
    return np.where(ok, v / mc.D, np.nan)


VMAX = {}
for q in ('k', 'L'):
    # colour limit from the strained points only, so that the large L of the weakly strained freestream saturates
    # instead of setting the scale; every point is drawn
    allv = np.concatenate([field(q, rg, f'x{x}', strained_only=True).ravel() for rg in DAT for x in XS])
    allv = allv[np.isfinite(allv)]
    VMAX[q] = float(np.percentile(allv, 99.5))
print('colour limits', VMAX)

W, H = 7.9, 5.2
pw, ph = 1.95, 0.98                      # axes size (in)
x0, gx = 0.62, 0.22
y0, gy, gm = 0.5, 0.12, 0.28
fig = plt.figure(figsize=(W, H))
tags = 'abcdefghijkl'
ims = {}
ytop = [y0 + (3 - r) * (ph + gy) + (gm if r < 2 else 0) for r in range(4)]
for r, (q, rung) in enumerate(rows):
    for c, x in enumerate(XS):
        p = f'x{x}'
        ax = fig.add_axes([(x0 + c * (pw + gx)) / W, ytop[r] / H, pw / W, ph / H])
        ax.set_facecolor('0.85')
        d = DAT[rung]
        Y, Z = d[f'{p}_A'] / mc.D, d[f'{p}_B'] / mc.D
        v = field(q, rung, p)
        norm = (LogNorm(vmin=L_LOG[0], vmax=L_LOG[1]) if q == 'L'
                else Normalize(vmin=0, vmax=VMAX[q]))
        ims[q] = ax.pcolormesh(Y, Z, np.ma.masked_invalid(v), cmap='inferno', norm=norm,
                               shading='auto', rasterized=True)
        e = mc.eta(d[f'{p}_T'].astype(float))
        ax.contour(Y, Z, np.ma.masked_invalid(e), levels=ETA_LEVELS, colors='white', linewidths=0.6)
        ax.set_aspect('equal')
        ax.set_xlim(-3, 3)
        ax.set_ylim(0, 3)
        ax.set_yticks([0, 1, 2, 3])
        ax.text(0.02, 0.94, f'({tags[r * 3 + c]}) $x/D={x}$', transform=ax.transAxes, color=('white' if q == 'k' else 'black'), va='top',
                fontsize=9)
        if c == 0:
            ax.set_ylabel(f'{rung}\n$z/D$', fontsize=9)
        else:
            ax.set_yticklabels([])
        if r == 3:
            ax.set_xlabel(r'$y/D$')
        else:
            ax.set_xticklabels([])
        med = np.nanmedian(v)
        print(f'{q} {rung} x/D={x}: median {med:.4g}, p99 {np.nanpercentile(v, 99):.4g}, max {np.nanmax(v):.4g}')
cx = (x0 + 3 * pw + 2 * gx + 0.15) / W
for q, (ra, rb), lab in (('k', (0, 1), r'$k/U_\infty^2$'), ('L', (2, 3), r'$L_\mathrm{eff}/D$')):
    cax = fig.add_axes([cx, ytop[rb] / H, 0.14 / W, (ytop[ra] + ph - ytop[rb]) / H])
    cb = fig.colorbar(ims[q], cax=cax, extend=('both' if q == 'L' else 'max'))
    cb.set_label(lab)
fig.savefig(OUT, dpi=300)
print('saved', OUT)
