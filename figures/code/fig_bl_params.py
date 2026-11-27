#!/usr/bin/env python3
"""F10 fig_bl_params: shape factor H and skin friction C_f along the plate, centreline (y/D 0) and mid-pitch (|y/D| = 3).

Panels: (a) H centreline, (b) H mid-pitch, (c) C_f centreline, (d) C_f mid-pitch; x/D -5 to 20.
Rungs: RANS 5M, RANS 50M, 20M LES, 50M 3DNS, 100M 3DNS, 193M 3DNS (all six exist for H and C_f).
Grey lines: Blasius H = 2.59 and C_f = 0.664 / sqrt(Re_x); turbulent H = 1.4 and C_f = 0.0576 Re_x^-0.2.
  Re_x = Re_D (x/D - x0/D) with Re_D = U_inf D / nu from the RANS 50M approach layer (bl_summary.csv, x/D -5) and the virtual
  origin x0 fixed so that Blasius Re_theta = 0.664 sqrt(Re_x) equals the median Re_theta of the 50M 3DNS approach layer at x/D -5
  (the 3DNS has no leading edge; the plate boundary layer starts from the inflow box). The turbulent line uses the same Re_x.
Data
  3DNS / LES: stat3 means, plate blocks: batch b006 (C20bW, C50b, C100b) and b010 for the 193M (C193b, pictures_job18/stat3_b010), helper _bl_extract.py (edge rule of
              analysis/layer_state/layer_mean.py: |dU/dz| >= 3e-2 of the wall value, + 0.25 mm; Favre U; compressible
              delta*, theta; C_f = tau_w / (0.5 rho_e U_e^2)). Centreline: the wall column nearest y/D 0 at each x (|y/D| < 0.1).
              Mid-pitch: the columns at y/D = +3 and -3 (the same line, periodic) averaged at each x.
  RANS:       starccm/<case>/post_run1/bl_summary.csv (bl_post.py; probes bl_c at y/D 0 and bl_m at y/D 2.91, integer x/D,
              the hole gap x/D -1 to 2 has no probe). RANS edge rule is 1e-3 of the wall gradient, not 3e-2: the H values are
              of the same kind but not identical in the rule (layer_state README, method 2).
The centreline of the 3DNS / LES is blank inside the hole exit region (x/D -1.3 to 2.7, no plate wall). All lines end
at x/D 20 (x/D 20-25 is a buffer region and is not shown).
usage: python3 fig_bl_params.py
"""
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
from cycler import cycler  # noqa: E402

sys.path.insert(0, '/home/sk999/Documents/3DNS_FP')
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from functions import set_mpl_params  # noqa: E402
import _bl_extract as BX  # noqa: E402

OUT = '/home/sk999/Documents/ASME_2027/figures/fig_bl_params.pdf'
ST = '/home/sk999/Documents/4dns_filmcooled/starccm'
D, UINF = 0.00775, 86.2812
OKABE_ITO = ['#000000', '#E69F00', '#56B4E9', '#009E73', '#F0E442', '#0072B2', '#D55E00', '#CC79A7']
XMAX = 20.0
XLINE = 20.0       # all lines end at x/D 20 (owner rule 2026-11-26; x/D 20-25 is a buffer region, never shown)
# label, colour, style  (fixed per rung, figure brief)
RUNG = {'193M 3DNS': ('#000000', '-'), '100M 3DNS': ('#0072B2', '--'), '50M 3DNS': ('#D55E00', ':'),
        '20M LES': ('#E69F00', '-.'), 'RANS 50M': ('#009E73', '-'), 'RANS 5M': ('#CC79A7', '--')}
TAG = {'193M 3DNS': 'c193b', '100M 3DNS': 'c100b', '50M 3DNS': 'c50b', '20M LES': 'c20bw'}
RANS = {'RANS 50M': 'rans_50M_Lm25b_BR15_FST_ZPG_p0dns', 'RANS 5M': 'rans_5M_Lm25_BR15_FST_ZPG_p0dns'}


def line_c(w):
    """centreline: per unique x the column with the smallest |y| (< 0.1 D)."""
    x, y = w['x'], w['y']
    m = (np.abs(y) < 0.1) & (x > -6) & (x <= XLINE) & ~((x > -1.3) & (x < 2.7))   # hole exit region blanked
    idx = np.nonzero(m)[0]
    key = np.round(x[idx], 3)
    o = np.lexsort((np.abs(y[idx]), key))
    idx, key = idx[o], key[o]
    first = np.r_[True, np.diff(key) > 0]
    return idx[first]


def line_m(w):
    """mid-pitch: columns at y/D = +-3 (|y| > 2.995), averaged over the pair at each x."""
    x, y = w['x'], w['y']
    m = (np.abs(y) > 2.995) & (x > -6) & (x <= XLINE)
    idx = np.nonzero(m)[0]
    key = np.round(x[idx], 3)
    ux, inv = np.unique(key, return_inverse=True)
    n = np.bincount(inv)
    out = {q: np.bincount(inv, w[q][idx]) / n for q in ('H', 'cf', 'x', 'ret')}
    return out


def gapped(x, v, gap):
    o = np.argsort(x)
    x, v = x[o], v[o]
    br = np.nonzero(np.diff(x) > gap)[0]
    return np.insert(x, br + 1, np.nan), np.insert(v, br + 1, np.nan)


def rans_line(case, kind):
    d = pd.read_csv(f'{ST}/{case}/post_run1/bl_summary.csv')
    d = d[d['name'].str.startswith('bl_c' if kind == 'c' else 'bl_m')]
    d = d[(d['xD'] >= -5.5) & (d['xD'] <= 20.0)].sort_values('xD')
    return d


GREY = '0.35'
AKW = dict(fontsize=8.5, color=GREY, ha='center', va='center',
           arrowprops=dict(arrowstyle='->', color=GREY, lw=0.7, shrinkA=1, shrinkB=1))


def annotate(ax, xr, rex):
    """Text and arrow labels for every item that is not a case (no legend entries)."""
    for a in ax.ravel():                      # hole exit region, shaded
        a.axvspan(-1.3, 2.7, color='0.85', alpha=0.5, lw=0, zorder=0)
    ax[0, 0].annotate('hole exit region', xy=(0.7, 2.85), xytext=(0.7, 3.1), **AKW)
    ax[0, 1].annotate('hole exit region', xy=(0.7, 2.85), xytext=(0.7, 3.1), **AKW)
    ax[0, 0].annotate(r'laminar, $H=2.59$', xy=(13, 2.59), xytext=(13, 2.95), **AKW)
    ax[0, 0].annotate(r'turbulent, $H=1.4$', xy=(13, 1.4), xytext=(13, 1.85), **AKW)
    ax[0, 1].annotate(r'laminar, $H=2.59$', xy=(13, 2.59), xytext=(13, 2.95), **AKW)
    ax[0, 1].annotate(r'turbulent, $H=1.4$', xy=(13, 1.4), xytext=(13, 1.1), **AKW)
    for j in (0, 1):
        cl = lambda xq: 0.664 / np.sqrt(ReD_EX[0] * (xq - ReD_EX[1])) * 1e3
        ct = lambda xq: 0.0576 * (ReD_EX[0] * (xq - ReD_EX[1])) ** -0.2 * 1e3
        xl, yl = (0.7, 4.3) if j == 0 else (5.5, 4.4)
        xa = 0.7 if j == 0 else 8.0
        ax[1, j].annotate(r'laminar, $0.664\,Re_x^{-1/2}$', xy=(xa, cl(xa)), xytext=(xl, yl), **AKW)
        ax[1, j].annotate(r'turbulent, $0.0576\,Re_x^{-1/5}$', xy=(12, ct(12)), xytext=(12, 7.0), **AKW)
        ax[1, j].annotate('hole exit region', xy=(0.7, 6.6), xytext=(0.7, 7.2), **AKW)


ReD_EX = [None, None]


def main():
    set_mpl_params(0.48, grid_on=True)
    plt.rcParams['axes.prop_cycle'] = cycler(color=OKABE_ITO)
    W = {k: BX.load(t) for k, t in TAG.items()}
    # reference: Re_D from the RANS 50M approach layer, virtual origin from the 50M 3DNS Re_theta at x/D -5
    r50 = pd.read_csv(f'{ST}/{RANS["RANS 50M"]}/post_run1/bl_summary.csv')
    r0 = r50[r50['name'] == 'bl_c_xD-5'].iloc[0]
    nu = r0['ue'] * r0['theta'] / r0['Re_theta']
    ReD = UINF * D / nu
    w50 = W['50M 3DNS']
    sel = (np.abs(w50['x'] + 5.0) < 0.15) & (np.abs(w50['y']) < 3.01)
    ret0 = np.median(w50['ret'][sel])
    rex0 = (ret0 / 0.664) ** 2
    x0 = -5.0 - rex0 / ReD
    xr = np.linspace(-5, XMAX, 300)
    rex = ReD * (xr - x0)
    ReD_EX[0], ReD_EX[1] = ReD, x0

    fig, ax = plt.subplots(2, 2, figsize=(7.2, 7.4), sharex=True)
    lab = {(0, 0): '(a) centreline', (0, 1): '(b) mid-pitch', (1, 0): '(c) centreline', (1, 1): '(d) mid-pitch'}
    for j, kind in enumerate(('c', 'm')):
        for k, (col, ls) in RUNG.items():
            if k in W:
                w = W[k]
                if kind == 'c':
                    i = line_c(w)
                    xx, hh, cc = w['x'][i], w['H'][i], w['cf'][i]
                else:
                    o = line_m(w)
                    xx, hh, cc = o['x'], o['H'], o['cf']
                xh, hv = gapped(xx, hh, 0.6)
                xc, cv = gapped(xx, cc, 0.6)
                ax[0, j].plot(xh, hv, color=col, ls=ls, label=k)
                ax[1, j].plot(xc, cv * 1e3, color=col, ls=ls)
            else:
                d = rans_line(RANS[k], kind)
                xh, hv = gapped(d['xD'].values, d['H'].values, 1.5)
                xc, cv = gapped(d['xD'].values, d['cf'].values, 1.5)
                ax[0, j].plot(xh, hv, color=col, ls=ls, marker='o', ms=3, label=k)
                ax[1, j].plot(xc, cv * 1e3, color=col, ls=ls, marker='o', ms=3)
        # references
        for a in (ax[0, j],):
            a.axhline(2.59, color='0.55', ls='-', lw=1.0)
            a.axhline(1.4, color='0.55', ls=':', lw=1.0)
        ax[1, j].plot(xr, 0.664 / np.sqrt(rex) * 1e3, color='0.55', ls='-', lw=1.0)
        ax[1, j].plot(xr, 0.0576 * rex ** -0.2 * 1e3, color='0.55', ls=':', lw=1.0)
    for (i, j), t in lab.items():
        ax[i, j].text(0.02, 0.92, t, transform=ax[i, j].transAxes, ha='left', va='center')
        ax[i, j].set_xlim(-5, XMAX)
    ax[0, 0].set_ylim(0.8, 3.6)
    ax[0, 1].set_ylim(0.8, 3.6)
    for j in (0, 1):
        ax[1, j].set_ylim(-1.5, 9)
        ax[1, j].set_xlabel(r'$x/D$')
    ax[0, 0].set_ylabel(r'$H$')
    ax[1, 0].set_ylabel(r'$C_f\times10^{3}$')
    h, l = ax[0, 0].get_legend_handles_labels()
    fig.legend(h, l, loc='lower center', ncol=3, bbox_to_anchor=(0.5, -0.07))
    annotate(ax, xr, rex)
    fig.tight_layout()
    fig.savefig(OUT, bbox_inches='tight')
    print(f'Re_D {ReD:.0f}, nu {nu:.3e}, Re_theta(-5) {ret0:.1f}, x0/D {x0:.2f}')
    # numbers
    for kind in ('c', 'm'):
        for k in W:
            w = W[k]
            if kind == 'c':
                i = line_c(w); xx, hh, cc = w['x'][i], w['H'][i], w['cf'][i]
            else:
                o = line_m(w); xx, hh, cc = o['x'], o['H'], o['cf']
            for xq in (-3, 5, 10, 15, 20):
                j = np.argmin(np.abs(xx - xq))
                print(kind, k, xq, f'H {hh[j]:.2f} cf {cc[j]*1e3:.2f}')
        for k in RANS:
            d = rans_line(RANS[k], kind)
            for xq in (-3, 5, 10, 15, 20):
                j = np.argmin(np.abs(d['xD'].values - xq))
                print(kind, k, xq, f'H {d["H"].values[j]:.2f} cf {d["cf"].values[j]*1e3:.2f}')


if __name__ == '__main__':
    main()
