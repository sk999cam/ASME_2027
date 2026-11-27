#!/usr/bin/env python3
"""F9 fig_spots: turbulent spots between the jets, 193M 3DNS (b010, left column) and 20M LES (b006, right column).

(a), (b) Space-time diagrams of the ESTIMATED wall skin friction C_f,est x 10^3 along the mid-pitch line y/D = 3, x/D -3 to 20,
         t U_inf/D 150 to 190 (193M) and 272.8 to 312.8 (LES); one inferno scale (1st to 99th percentile of the shown data).
         C_f,est = 2 mu u(z/D = 0.1) / (0.1 D) / (rho_inf U_inf^2) from the z/D 0.1 plane sampler (one-sided gradient, local
         Sutherland mu from the sampler temperature); line y/D = +3 (identical to y/D = -3, periodic). Cached by
         _spots_z010_extract.py; constants and the validation against the true wall C_f are in fig_spots_z010.py.
         Dashed lines: hole-exit extent, dotted line: end of the wall data (x/D 14.9). Ticks t1..t4 at the right edge: times of the
         plan views below.
(c), (d) Plan views of the true wall C_f x 10^3 (first grid point above the wall, smoothed over 1.2 D/U) at t1..t4, y/D folded so
         the mid-pitch line y/D = 3 is the centre; hole outline (white, black border) with its periodic copy; white outline: the
         one spot tracked in time (n2_spots_cf). Wall data end at x/D 14.9, the axes run to x/D 20 so that the columns align.
Data: /mnt/Elements/DNS/analysis/asme2027_figs/n2_spots_cf/ (plan views) and .../spots_z010/ (z/D 0.1 lines).
The wall-only version of this figure is kept in fig_spots_wall.py. usage: python3 fig_spots.py
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Polygon  # noqa: E402
from cycler import cycler  # noqa: E402

sys.path.insert(0, '/home/sk999/Documents/3DNS_FP')
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from functions import set_mpl_params  # noqa: E402
from hole_exit_outline import load as load_hole  # noqa: E402
import fig_spots_z010 as Z  # noqa: E402

OUT = '/home/sk999/Documents/ASME_2027/figures/fig_spots.png'
DDIR = '/mnt/Elements/DNS/analysis/asme2027_figs/n2_spots_cf/'
OKABE_ITO = ['#000000', '#E69F00', '#56B4E9', '#009E73', '#F0E442', '#0072B2', '#D55E00', '#CC79A7']
PMIN, PMAX = 0.0, 14.0         # plan views, C_f x 1e3
XL = (-3.0, 20.0)
XH = (-1.0, 2.35)
XWALL = 14.9
WIN = Z.WIN
PICK = [1, 2, 3, 5]            # plan-view records shown (of the six saved)
CASES = [('c193b', '193M 3DNS'), ('c20bw', '20M LES')]
FW, FH = 7.2, 6.8             # figure size, in
BG = dict(facecolor='k', edgecolor='none', pad=1.5, alpha=0.75)


def fold(y):
    return np.where(y < 0, y + 6.0, y)


def main():
    set_mpl_params(0.5, grid_on=False)
    plt.rcParams['axes.prop_cycle'] = cycler(color=OKABE_ITO)
    HX, HY = load_hole()
    V = {t: Z.load(t) for t, _ in CASES}
    # one colour scale for the top row, from the data shown
    shown = []
    for tag, _ in CASES:
        v = V[tag]
        m = (v['t'] >= WIN[tag][0]) & (v['t'] <= WIN[tag][1])
        shown.append((v['cfe'][m][:, v['x'] <= XL[1] + 1e-9] * 1e3).ravel())
    vmin, vmax = np.percentile(np.concatenate(shown), [1, 99])
    print(f'top-row colour range C_f,est x 1e3: {vmin:.2f} to {vmax:.2f}; plan views {PMIN} to {PMAX}')

    fig = plt.figure(figsize=(FW, FH))
    LM, RM, GAP = 0.075, 0.895, 0.06
    W = (RM - LM - GAP) / 2
    cols = [LM, LM + W + GAP]
    inch = lambda a: a / FH
    H_TOP, PH, PG = inch(1.75), inch(0.66), inch(0.20)
    top0 = 1 - inch(0.08)
    b_top = top0 - H_TOP
    p_top = b_top - inch(0.2)
    im = imp = None
    for ci, (tag, name) in enumerate(CASES):
        v = V[tag]
        pv = np.load(DDIR + tag + '_planviews.npz')
        t0, t1 = WIN[tag]
        ax = fig.add_axes([cols[ci], b_top, W, H_TOP])
        m = (v['t'] >= t0) & (v['t'] <= t1)
        ix = v['x'] <= XL[1] + 1e-9
        tt = v['t'][m]
        dtm = np.median(np.diff(tt))
        im = ax.imshow(v['cfe'][m][:, ix] * 1e3, origin='lower', aspect='auto', cmap='inferno', vmin=vmin, vmax=vmax,
                       extent=[v['x'][0] - 0.01, v['x'][ix][-1] + 0.01, tt[0] - dtm / 2, tt[-1] + dtm / 2], interpolation='nearest')
        ax.set_xlim(*XL)
        ax.set_ylim(t0, t1)
        for xh in XH:
            ax.axvline(xh, color='w', ls='--', lw=0.9)
        ax.axvline(XWALL, color='w', ls=':', lw=1.2)
        ax.text(0.5 * sum(XH), t1 - 0.03 * (t1 - t0), 'hole\nexit', color='w', ha='center', va='top', fontsize=8)
        ax.text(XWALL - 0.3, t1 - 0.03 * (t1 - t0), 'end of wall data', color='w', ha='right', va='top', fontsize=8)
        ax.set_xticklabels([])
        if ci == 0:
            ax.set_ylabel(r'$tU_\infty/D$')
        ax.text(0.012, 0.03, '(%s) %s' % ('ab'[ci], name), transform=ax.transAxes, color='w', ha='left', va='bottom',
                fontsize=10, bbox=BG)
        for k, tk in enumerate(pv['times'][PICK]):
            ax.plot([XL[1] - 0.5, XL[1]], [tk, tk], color='w', lw=1.0)
            ax.text(XL[1] - 0.7, tk, r'$t_%d$' % (k + 1), color='w', ha='right', va='center', fontsize=8, bbox=BG)
        # plan views
        for k, ip in enumerate(PICK):
            top = p_top - k * (PH + PG)
            ap = fig.add_axes([cols[ci], top - PH, W, PH])
            c = pv['Cf_smooth_1p1'][ip]
            y = pv['y']
            yf = np.concatenate([y[y >= 0], y[y < 0] + 6.0])
            cf = np.concatenate([c[:, y >= 0], c[:, y < 0]], axis=1)
            o = np.argsort(yf)
            yf, cf = yf[o], cf[:, o]
            imp = ap.pcolormesh(pv['x'], yf, cf.T * 1e3, cmap='inferno', vmin=PMIN, vmax=PMAX, shading='nearest', rasterized=True)
            for dy in (0.0, 6.0):
                ap.add_patch(Polygon(np.c_[HX, HY + dy], closed=True, facecolor='w', edgecolor='k', lw=0.8, zorder=3))
            so = pv['spot_outlines'][ip]
            so = so[np.isfinite(so).all(1)]
            if len(so):
                ap.plot(so[:, 0], fold(so[:, 1]), color='w', lw=1.0, zorder=4)
            ap.plot([XL[0], pv['x'][-1]], [3.0, 3.0], color=OKABE_ITO[2], ls='--', lw=0.9, zorder=5)
            ap.axvline(XWALL, color='0.5', ls=':', lw=0.8, zorder=2)
            ap.set_xlim(*XL)
            ap.set_ylim(0, 6)
            ap.set_yticks([0, 3, 6])
            ap.text(0.012, 0.88, r'$t_%d$ = %.1f' % (k + 1, pv['times'][ip]), transform=ap.transAxes, color='w',
                    ha='left', va='top', fontsize=8, zorder=6, bbox=BG)
            if k == 0:
                ap.text(0.012, 0.12, '(%s)' % 'cd'[ci], transform=ap.transAxes, color='w', ha='left', va='bottom',
                        fontsize=10, zorder=6, bbox=BG)
            if ci == 0:
                ap.set_ylabel(r'$y/D$')
            else:
                ap.set_yticklabels([])
            if k < len(PICK) - 1:
                ap.set_xticklabels([])
            else:
                ap.set_xlabel(r'$x/D$')
            if ci == 0 and k == 0:
                ap.annotate(r'$y/D=3$ (mid-pitch)', xy=(11.5, 3.0), xytext=(8.0, 4.7), color='w', fontsize=8, ha='center',
                            va='center', bbox=BG, arrowprops=dict(arrowstyle='->', color='w', lw=0.8), zorder=7)
                ap.text(0.5 * (XWALL + XL[1]), 4.6, 'no wall\ndata', color='k', fontsize=8, ha='center', va='center', zorder=7)
            if ci == 0 and k == 1:
                ap.annotate('tracked spot', xy=(so[:, 0].max(), fold(so[:, 1]).mean()), xytext=(11.5, 4.9), color='w',
                            fontsize=8, ha='center', va='center', bbox=BG,
                            arrowprops=dict(arrowstyle='->', color='w', lw=0.8), zorder=7)
    pb = p_top - 3 * (PH + PG) - PH
    cax = fig.add_axes([RM + 0.012, b_top, 0.014, H_TOP])
    fig.colorbar(im, cax=cax).set_label(r'$C_{f,\mathrm{est}}\times10^{3}$ ($z/D=0.1$)')
    cax2 = fig.add_axes([RM + 0.012, pb, 0.014, p_top - pb])
    fig.colorbar(imp, cax=cax2).set_label(r'$C_f\times10^{3}$ (wall)')
    fig.savefig(OUT, dpi=300, bbox_inches='tight')
    print('wrote', OUT)


if __name__ == '__main__':
    main()
