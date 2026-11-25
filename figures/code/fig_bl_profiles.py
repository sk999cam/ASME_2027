#!/usr/bin/env python3
"""F10 fig_bl_profiles: wall-normal profiles at x/D = 10 of U/U_inf, k/U_inf^2 and gamma against z/D.

Rows: centreline (y/D 0), y/D 1.5 and mid-pitch. Columns: (U/U_inf), (k/U_inf^2), (gamma). z/D 0 to 1.5.
Rungs: RANS 50M, 20M LES, 193M 3DNS (b010). Complete: every rung has its own gamma, no stand-in.
Data
  RANS 50M  centreline and mid-pitch: starccm/rans_50M_Lm25b_BR15_FST_ZPG_p0dns/post_run1/bl_c_xD10.csv and bl_m_xD10.csv
            (mid-pitch probe at y/D 2.91; U, modelled k, gamma). y/D 1.5 (mean of +1.5 and -1.5): the RC-S1 cell export binned in
            0.25 D boxes (/mnt/Tank/DNS/analysis/layer_state/rans_bins_sst.npz; volume-weighted mean of the two x boxes
            9.875 and 10.125); box means, so coarser in z/D than the probes.
  20M LES (stat3 b006), 193M 3DNS (stat3 b010)  Favre mean U and resolved k = 1/2 (R_xx + R_yy + R_zz), helper _bl_extract.py; column nearest
            x/D 10 (within 0.15 D) of each y/D (+-1.5 and +-3 averaged; mid-pitch is the periodic line y/D = +-3).
  20M LES, 193M 3DNS gamma  plane sampler x10 (x/D 10, y/D -3..3, z/D 0..3), records of the b006 (20M) and b010 (193M) windows
            (3 992 records for the 193M, see _gamma_extract.py); detector |du/dt| smoothed over 1.1 D/U, one global threshold
            0.0385 U^2/D (Otsu of the 50M 3DNS); y/D 0, +-1.5 and +-3 (averaged, mid-pitch is the periodic line), blanked below
            z/D 0.02 (cache /mnt/Elements/DNS/analysis/asme2027_figs/gamma/).
RANS k is the modelled k; 3DNS / LES k is resolved.
usage: python3 fig_bl_profiles.py
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
import _gamma_extract as GX  # noqa: E402

OUT = '/home/sk999/Documents/ASME_2027/figures/fig_bl_profiles.pdf'
ST = '/home/sk999/Documents/4dns_filmcooled/starccm/rans_50M_Lm25b_BR15_FST_ZPG_p0dns/post_run1'
BINS = '/mnt/Tank/DNS/analysis/layer_state/rans_bins_sst.npz'
D, UINF = 0.00775, 86.2812
OKABE_ITO = ['#000000', '#E69F00', '#56B4E9', '#009E73', '#F0E442', '#0072B2', '#D55E00', '#CC79A7']
STY = {'RANS 50M': ('#009E73', '-'), '20M LES': ('#E69F00', '-.'), '193M 3DNS': ('#000000', '-')}
ZMAX = 1.5
LINES = [(0.0, 'centreline'), (1.5, r'$y/D=1.5$'), (3.0, 'mid-pitch')]


def rans_probe(kind):
    d = pd.read_csv(f'{ST}/bl_{kind}_xD10.csv').sort_values('Z (m)')
    return dict(z=d['Z (m)'].values / D, U=d['Velocity[i] (m/s)'].values / UINF,
                k=d['Turbulent Kinetic Energy (J/kg)'].values / UINF ** 2, g=d['Intermittency'].values)


def rans_boxes(y):
    B = np.load(BINS)
    ix = [43, 44]
    iys = [12] if y == 0 else [int(round((s * y + 3.125 - 0.125) / 0.25)) for s in (1, -1)]
    s = {f: sum(B[f][i, j] for i in ix for j in iys) for f in ('V', 'Vz', 'Vu', 'Vg', 'Vk')}
    occ = s['V'] > 0
    V = s['V'][occ]
    return dict(z=s['Vz'][occ] / V, U=s['Vu'][occ] / V / UINF, k=s['Vk'][occ] / V / UINF ** 2, g=s['Vg'][occ] / V)


def dns_prof(tag, y):
    p = BX.load(tag)['prof']
    ks = [0.0] if y == 0 else [y, -y]
    z = p[ks[0]]['z']
    U = np.mean([np.interp(z, p[k]['z'], p[k]['U']) for k in ks], 0)
    k = np.mean([np.interp(z, p[k]['z'], p[k]['k']) for k in ks], 0)
    return dict(z=z, U=U, k=k, dist=max(p[k]['dist'] for k in ks))


def gamma_dns(tag, y):
    P = GX.load_x10(tag)
    G = P['gamma'][0]                                   # (y, z) at threshold 0.0385
    ys = [0.0] if y == 0 else [y, -y]
    g = np.mean([G[int(np.argmin(np.abs(P['y'] - s)))] for s in ys], 0)
    m = P['z'] >= 0.02
    return P['z'][m], g[m]


def main():
    set_mpl_params(0.7, grid_on=True)
    plt.rcParams['axes.prop_cycle'] = cycler(color=OKABE_ITO)
    fig, ax = plt.subplots(3, 3, figsize=(12.5, 12), sharey=True)
    pa = 0
    for i, (y, ylab) in enumerate(LINES):
        rp = rans_probe('c') if y == 0 else (rans_probe('m') if y == 3.0 else rans_boxes(y))
        data = {'RANS 50M': rp, '20M LES': dns_prof('c20bw', y), '193M 3DNS': dns_prof('c193b', y)}
        for nm in ('RANS 50M', '20M LES', '193M 3DNS'):
            c, ls = STY[nm]
            d = data[nm]
            ax[i, 0].plot(d['U'], d['z'], color=c, ls=ls, label=nm)
            ax[i, 1].plot(d['k'], d['z'], color=c, ls=ls)
            if nm == 'RANS 50M':
                ax[i, 2].plot(d['g'], d['z'], color=c, ls=ls)
            else:
                zg, gg = gamma_dns('c20bw' if nm == '20M LES' else 'c193b', y)
                ax[i, 2].plot(gg, zg, color=c, ls=ls)
        for j in range(3):
            a = ax[i, j]
            a.text(0.97 if j == 1 else 0.03, 0.06 if (i == 0 and j == 1) else 0.93, f'({chr(97 + pa)}) {ylab}', transform=a.transAxes,
                   ha='right' if j == 1 else 'left', va='center')
            pa += 1
            a.set_ylim(0, ZMAX)
        ax[i, 0].set_xlim(0, 1.1)
        ax[i, 2].set_xlim(-0.03, 1.05)
        ax[i, 0].set_ylabel(r'$z/D$')
    for j, xl in enumerate((r'$\tilde{U}/U_\infty$', r'$k/U_\infty^2$', r'$\gamma$')):
        ax[2, j].set_xlabel(xl)
    h, l = ax[0, 0].get_legend_handles_labels()
    fig.legend(h, l, loc='lower center', ncol=3, bbox_to_anchor=(0.5, -0.03))
    fig.tight_layout()
    fig.savefig(OUT, bbox_inches='tight')
    for i, (y, ylab) in enumerate(LINES):
        for tag in ('c20bw', 'c193b'):
            d = dns_prof(tag, y)
            j = np.argmin(np.abs(d['z'] - 0.1))
            print(y, tag, 'dist', round(d['dist'], 3), 'k(z=0.1)', d['k'][j], 'U(z=0.1)', d['U'][j], 'kmax', d['k'].max())
        for tag in ('c20bw', 'c193b'):
            zg, gg = gamma_dns(tag, y)
            print(y, tag, 'gamma z~0.1', gg[np.argmin(np.abs(zg - 0.1))], 'z~0.5', gg[np.argmin(np.abs(zg - 0.5))])


if __name__ == '__main__':
    main()
