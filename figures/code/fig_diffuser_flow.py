"""Mean axial velocity on cross-sections of the diffuser (paper Fig. F6, Results 4.1).

Usage: python3 fig_diffuser_flow.py   -> /home/sk999/Documents/ASME_2027/figures/fig_diffuser_flow.png

u_s = mean velocity component along the 30 deg hole axis (u cos30 + w sin30; Favre mean for the 3DNS/LES), normalised
by the freestream velocity U_inf = 86.28 m/s. Colour map bwr, symmetric about zero (limits +-0.05; forward flow
above 0.05 saturates, colour bar extends both ends); black line u_s = 0.
Rows: sections normal to the hole axis at s/D = -2.3 (peak of the reversed area fraction, F5 b), -1.1 and -0.5 (near
the exit). Columns: RANS 50M, 20M LES (stat3 batch b006), 193M 3DNS (reference, stat3 batch b010).
Axes: y/D, n/D (n positive towards the upstream wall), looking along the hole axis, equal aspect in every panel.
Method (diffuser_wall.py): linear interpolation from the cell centres within +-one cell-column spacing of the
section (RANS 50M 0.08 D, 20M LES 0.11 D, 193M 0.045 D), coloured only inside the hull of those cell centres (the
fluid); cached in /mnt/Elements/DNS/analysis/asme2027_figs/hole_b010/sec_<case>_<s0>.npz.
The reversed layer on the centre plane (0.03-0.07 D thick) is given as numbers in the text, not drawn.
Column labels are ax.text, not titles. Prints the reversed area fraction of each panel (pixels with u_s < 0 over
pixels in the hull) against sections_<case>.txt (F5 b). All rungs complete.
"""
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from cycler import cycler
from matplotlib.colors import Normalize

sys.path.insert(0, '/home/sk999/Documents/3DNS_FP')
sys.path.insert(0, '/home/sk999/Documents/ASME_2027/figures/code')
from functions import set_mpl_params  # noqa: E402
import diffuser_wall as W  # noqa: E402

OKABE_ITO = ['#000000', '#E69F00', '#56B4E9', '#009E73', '#F0E442', '#0072B2', '#D55E00', '#CC79A7']
OUT = '/home/sk999/Documents/ASME_2027/figures/fig_diffuser_flow.png'
RES = '/home/sk999/Documents/4dns_filmcooled/analysis/hole_discharge/results'
SEC193 = W.CACHE + '/sections_c193b_b010.txt'
U_INF = 86.28
LIM = 0.05
COLS = [('rans50m', 'RANS 50M'), ('c20bw', '20M LES'), ('c193b', '193M 3DNS')]
SLIST = (-2.3, -1.1, -0.5)

set_mpl_params(0.5, grid_on=False)
plt.rcParams['axes.prop_cycle'] = cycler(color=OKABE_ITO)

V = {k: W.sections_for(k, SLIST) for k, _ in COLS}
norm = Normalize(-LIM, LIM)
fig = plt.figure(figsize=(7.2, 5.2), constrained_layout=True)
axs = fig.subplots(3, 3, sharex=True, sharey=True)
for r, s0 in enumerate(SLIST):
    for j, (k, lab) in enumerate(COLS):
        v = V[k][s0] / U_INF
        ax = axs[r, j]
        im = ax.pcolormesh(W.GY, W.GN, v.T, cmap='bwr', norm=norm, shading='auto', rasterized=True)
        ax.contour(W.GY, W.GN, np.ma.masked_invalid(v.T), levels=[0], colors='k', linewidths=0.8)
        ax.set_aspect('equal')
        ax.set_xlim(-1.1, 1.1)
        ax.set_ylim(-1.0, 1.1)
        ax.text(0.04, 0.95, f'({"abcdefghi"[3 * r + j]}) {lab}\n$s/D={s0:.1f}$'.replace('-', '−'),
                transform=ax.transAxes, ha='left', va='top')
        d = np.genfromtxt(SEC193 if k == 'c193b' else f'{W.P.RES}/sections_{k}.txt', delimiter=',', names=True)
        f5 = float(d['rev'][np.argmin(np.abs(d['s'] - s0))])
        ok = np.isfinite(v)
        print(f'{lab:10s} s/D {s0:+.1f}: reversed area fraction {np.sum(v[ok] < 0) / ok.sum():.3f}  (F5 b {f5:.3f})')
        if r == 2:
            ax.set_xlabel(r'$y/D$')
        if j == 0:
            ax.set_ylabel(r'$n/D$')
fig.colorbar(im, ax=axs, orientation='vertical', fraction=0.04, extend='both',
             label=r'$\langle u_s\rangle/U_\infty$')
fig.savefig(OUT, dpi=300, bbox_inches='tight')
print('wrote', OUT)
