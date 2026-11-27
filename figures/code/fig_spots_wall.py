#!/usr/bin/env python3
"""F9 fig_spots: turbulent spots between the jets seen in the wall skin friction, 193M 3DNS (b010) and 20M LES (b006).

(a), (b) Space-time diagrams of C_f (x 10^3) along the mid-pitch line y/D = 3, x/D -3 to 14.9 against t U_inf/D, 193M and 20M
         LES, same axes and colour bar. C_f = 2 mu_w u1/(z1 rho U^2) from the first grid point above the wall, z1 = 0.00516 D.
         Field shown: running mean over 1.2 D/U (the detection field of the spot statistics). Window of 40 D/U in both cases.
         Dotted band: hole-exit extent in x. Horizontal ticks t1..t4: times of the plan views below the same column.
(c), (d) Plan views of C_f at the times t1..t4, y/D folded so that the mid-pitch line y/D = 3 is the centre of the panel;
         hole outline (white fill, black border) at both ends (periodic copy); white outline: the one spot tracked in time.
Data: /mnt/Elements/DNS/analysis/asme2027_figs/n2_spots_cf/ (summary.txt). usage: python3 fig_spots.py
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Polygon  # noqa: E402

sys.path.insert(0, '/home/sk999/Documents/3DNS_FP')
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from functions import set_mpl_params  # noqa: E402
from hole_exit_outline import load  # noqa: E402

OUT = '/home/sk999/Documents/ASME_2027/figures/code/fig_spots_wall.png'   # kept copy of the wall-only version (not a paper figure)
DDIR = '/mnt/Elements/DNS/analysis/asme2027_figs/n2_spots_cf/'
OKABE_ITO = ['#000000', '#E69F00', '#56B4E9', '#009E73', '#F0E442', '#0072B2', '#D55E00', '#CC79A7']
VMIN, VMAX = 1.0, 7.0          # space-time, C_f x 1e3 (1st-99th percentile of the plotted windows)
PMIN, PMAX = 0.0, 14.0         # plan views, C_f x 1e3 (0-99.5th percentile of the snapshots, linear)
XL = (-3.0, 14.88)
XH = (-1.0, 2.35)
WIN = {'c193b': (150.0, 190.0), 'c20bw': (272.8, 312.8)}
PICK = [1, 2, 3, 5]            # plan-view records shown (of the six saved)

set_mpl_params(0.5, grid_on=False)
HX, HY = load()
CM = 'inferno'


def fold(y):
    return np.where(y < 0, y + 6.0, y)


fig = plt.figure(figsize=(7.2, 7.6))
LM, RM = 0.075, 0.915
GAP = 0.07
W = (RM - LM - GAP) / 2
cols = [LM, LM + W + GAP]
H_TOP, B_TOP = 0.255, 0.70
PH, PG = 0.105, 0.03
cases = [('c193b', '193M 3DNS'), ('c20bw', '20M LES')]
lab = ['(a)', '(b)', '(c)', '(d)']
im = None
for ci, (tag, name) in enumerate(cases):
    st = np.load(DDIR + tag + '_spacetime.npz')
    pv = np.load(DDIR + tag + '_planviews.npz')
    t0, t1 = WIN[tag]
    ax = fig.add_axes([cols[ci], B_TOP + 0.02, W, H_TOP])
    m = (st['t'] >= t0) & (st['t'] <= t1)
    im = ax.imshow(st['Cf_smooth_1p1'][m] * 1e3, origin='lower', aspect='auto', cmap=CM, vmin=VMIN, vmax=VMAX,
                   extent=[st['x'][0] - 0.02, st['x'][-1] + 0.02, st['t'][m][0], st['t'][m][-1]], interpolation='nearest')
    ax.set_xlim(*XL)
    ax.set_ylim(t0, t1)
    ax.axvspan(*XH, facecolor='none', edgecolor='w', ls=':', lw=1.0, hatch=None)
    ax.text(0.5 * sum(XH), t1 - 0.03 * (t1 - t0), 'hole exit', color='w', ha='center', va='top', fontsize=8)
    ax.set_xticklabels([])
    ax.set_ylabel(r'$tU_\infty/D$')
    ax.text(0.015, 0.97, lab[ci], transform=ax.transAxes, color='w', ha='left', va='top', fontsize=10)
    ax.text(0.985, 0.97, name, transform=ax.transAxes, color='w', ha='right', va='top', fontsize=9)
    ts = pv['times'][PICK]
    for k, tk in enumerate(ts):
        ax.plot([XL[1] - 0.45, XL[1]], [tk, tk], color='w', lw=1.0)
        ax.text(XL[1] - 0.6, tk, r'$t_%d$' % (k + 1), color='w', ha='right', va='center', fontsize=8, bbox=dict(facecolor='k', edgecolor='none', pad=0.8, alpha=0.6))
    # plan views
    for k, ip in enumerate(PICK):
        top = B_TOP - 0.0 - k * (PH + PG)
        ap = fig.add_axes([cols[ci], top - PH, W, PH])
        c = pv['Cf_smooth_1p1'][ip]                     # (x, y)
        y = pv['y']
        yf = np.concatenate([y[y >= 0], y[y < 0] + 6.0])
        cf = np.concatenate([c[:, y >= 0], c[:, y < 0]], axis=1)
        o = np.argsort(yf)
        yf, cf = yf[o], cf[:, o]
        imp = ap.pcolormesh(pv['x'], yf, cf.T * 1e3, cmap=CM, vmin=PMIN, vmax=PMAX, shading='nearest', rasterized=True)
        for dy in (0.0, 6.0):
            ap.add_patch(Polygon(np.c_[HX, HY + dy], closed=True, facecolor='w', edgecolor='k', lw=0.8, zorder=3))
        so = pv['spot_outlines'][k + 0 if False else ip]
        so = so[np.isfinite(so).all(1)]
        if len(so):
            ap.plot(so[:, 0], fold(so[:, 1]), color='w', lw=1.0, ls='-', zorder=4)
        ap.axhline(3.0, color=OKABE_ITO[2], ls='--', lw=0.9, zorder=5)
        ap.set_xlim(*XL)
        ap.set_ylim(0, 6)
        ap.set_yticks([0, 3, 6])
        ap.text(0.012, 0.88, r'$t_%d$ = %.1f' % (k + 1, pv['times'][ip]), transform=ap.transAxes, color='w',
                ha='left', va='top', fontsize=8, zorder=6, bbox=dict(facecolor='k', edgecolor='none', pad=1.5, alpha=0.75))
        if k == 0:
            ap.text(0.012, 0.12, lab[2 + ci], transform=ap.transAxes, color='w', ha='left', va='bottom', fontsize=10, zorder=6, bbox=dict(facecolor='k', edgecolor='none', pad=1.5, alpha=0.75))
        if ci == 0:
            ap.set_ylabel(r'$y/D$')
        else:
            ap.set_yticklabels([])
        if k < len(PICK) - 1:
            ap.set_xticklabels([])
        else:
            ap.set_xlabel(r'$x/D$')
        if ci == 0 and k == 0:
            ap.annotate(r'$y/D=3$ (mid-pitch)', xy=(11.5, 3.0), xytext=(8.0, 4.6), color='w', fontsize=8, ha='center',
                        va='center', bbox=dict(facecolor='k', edgecolor='none', pad=1.5, alpha=0.75), arrowprops=dict(arrowstyle='->', color='w', lw=0.8), zorder=7)
        if ci == 0 and k == 1:
            ap.annotate('tracked spot', xy=(so[:, 0].max(), fold(so[:, 1]).mean()),
                        xytext=(11.5, 4.9), color='w', fontsize=8, ha='center', va='center', bbox=dict(facecolor='k', edgecolor='none', pad=1.5, alpha=0.75),
                        arrowprops=dict(arrowstyle='->', color='w', lw=0.8), zorder=7)
cax = fig.add_axes([RM + 0.012, B_TOP + 0.02, 0.014, H_TOP])
cb = fig.colorbar(im, cax=cax)
cb.set_label(r'$C_f\times10^{3}$')
pb = B_TOP - 3 * (PH + PG) - PH
cax2 = fig.add_axes([RM + 0.012, pb, 0.014, B_TOP - pb])
cb2 = fig.colorbar(imp, cax=cax2)
cb2.set_label(r'$C_f\times10^{3}$')
fig.savefig(OUT, dpi=300, bbox_inches='tight')
print('wrote', OUT)
