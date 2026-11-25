"""Pressure and reversed flow along the hole axis, six rungs (paper Fig. F5, Results 4.1).

Usage: python3 fig_hole_axis.py      -> /home/sk999/Documents/ASME_2027/figures/fig_hole_axis.pdf

(a) Area-mean static pressure of each cross-section normal to the 30 deg hole axis, as the coefficient
        C_p,c = (p - p_inf) / (p0c - p_inf)
    with p0c = 81 110 Pa abs (plenum total pressure, same for every rung) and p_inf = 76 596.2 Pa (reference
    pressure of the simulations, freestream static pressure). p0c - p_inf = 4513.8 Pa. C_p,c = 1 would be the
    plenum total pressure.
(b) Fraction of each section's area with reversed mean axial flow (u_s < 0).
s/D along the hole axis, 0 at the exit centre; metering section s/D < -3.5, diffuser -3.5 to the exit.
Data: analysis/hole_discharge/results/sections_<case>.txt (extract.py, every 0.2 D) for RANS 5M (rans5m),
RANS 50M (rans50m), 20M LES (c20bw), 50M / 100M 3DNS (c50b, c100b), stat3 batch b006. The 193M curve (the
reference) is batch b010 (10 FT): sections_c193b_b010.txt, the same extraction (hole_b010.py, with the functions of
extract.py / hole_loss.py), cached in /mnt/Elements/DNS/analysis/asme2027_figs/hole_b010/. All rungs complete.
Redraw of hole_discharge plots.py F3 to the CLAUDE.md figure rules.
"""
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from cycler import cycler

sys.path.insert(0, '/home/sk999/Documents/3DNS_FP')
from functions import set_mpl_params

OKABE_ITO = ['#000000', '#E69F00', '#56B4E9', '#009E73', '#F0E442', '#0072B2', '#D55E00', '#CC79A7']
RES = '/home/sk999/Documents/4dns_filmcooled/analysis/hole_discharge/results'
SEC193 = '/mnt/Elements/DNS/analysis/asme2027_figs/hole_b010/sections_c193b_b010.txt'
OUT = '/home/sk999/Documents/ASME_2027/figures/fig_hole_axis.pdf'
P0C, PINF = 81110.0, 76596.2
RUNGS = [('c193b', '193M 3DNS', '#000000', '-'), ('c100b', '100M 3DNS', '#0072B2', '--'),
         ('c50b', '50M 3DNS', '#D55E00', ':'), ('c20bw', '20M LES', '#E69F00', '-.'),
         ('rans50m', 'RANS 50M', '#009E73', '-'), ('rans5m', 'RANS 5M', '#CC79A7', '--')]

set_mpl_params(0.5, grid_on=True)
plt.rcParams['axes.prop_cycle'] = cycler(color=OKABE_ITO)

fig, axs = plt.subplots(2, 1, figsize=(5.0, 5.2), sharex=True)
for k, lab, c, ls in RUNGS:
    d = np.genfromtxt(SEC193 if k == 'c193b' else f'{RES}/sections_{k}.txt', delimiter=',', names=True)
    m = (d['s'] > -6.0) & (d['s'] < 0.0)
    axs[0].plot(d['s'][m], (d['p'][m] - PINF) / (P0C - PINF), color=c, ls=ls, label=lab)
    axs[1].plot(d['s'][m], d['rev'][m], color=c, ls=ls)
for a in axs:
    a.axvspan(-6, -3.5, color='0.5', alpha=0.12, lw=0)
    a.set_xlim(-6, 0)
axs[0].text(-4.75, 0.97, 'metering', ha='center', va='top', transform=axs[0].get_xaxis_transform())
axs[0].text(-1.75, 0.97, 'diffuser', ha='center', va='top', transform=axs[0].get_xaxis_transform())
axs[0].axvline(-3.5, color='0.4', lw=0.8)
axs[1].axvline(-3.5, color='0.4', lw=0.8)
axs[0].set_ylabel(r'$(p-p_\infty)/(p_{0c}-p_\infty)$')
axs[1].set_ylabel('reversed-flow area fraction')
axs[1].set_xlabel(r'$s/D$')
axs[0].text(0.03, 0.06, '(a)', transform=axs[0].transAxes)
axs[1].text(0.97, 0.88, '(b)', transform=axs[1].transAxes, ha='right')
axs[1].legend(*axs[0].get_legend_handles_labels(), loc='upper left', bbox_to_anchor=(0.0, 0.92))
fig.tight_layout()
fig.savefig(OUT, bbox_inches='tight')
print('wrote', OUT)
