"""Volume distribution of Delta/eta per region for the scale-resolving rungs (paper Fig. F3, interim).

Usage: python3 fig_resolution.py [OUT.pdf]
       default OUT: /home/sk999/Documents/ASME_2027/figures/fig_resolution.pdf
       copy of this script: ASME_2027/figures/code/fig_resolution.py
Input: the length_hist.npz written by dns/python/length_scales.py (Delta = cell volume^(1/3),
eta = (nu^3/eps)^(1/4) from the resolved viscous dissipation of the stat3 batch named below).

Interim: each run uses its own eta. The planned version (write-up plan, 2.3) uses Delta of each mesh over
eta of the 193M 3DNS batch b010, and adds the 5M mesh; it needs the mapping of each mesh onto the 193M points.
Style: set_mpl_params, Okabe-Ito colours with a line style per run, no titles (CLAUDE.md figure rules).
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
OUT = sys.argv[1] if len(sys.argv) > 1 else '/home/sk999/Documents/ASME_2027/figures/fig_resolution.pdf'

# (label, stat3 batch directory, colour, line style)
RUNS = [('20M', '/mnt/Tank/DNS/film_cooling/C20bW_BR15/pictures_job3/stat3_b006', OKABE_ITO[1], '-.'),
        ('50M', '/mnt/Tank/DNS/film_cooling/C50b_BR15/pictures_job9/stat3_b006', OKABE_ITO[6], ':'),
        ('100M', '/mnt/Elements/DNS/film_cooling/C100b_BR15/pictures_job9/stat3_b006', OKABE_ITO[5], '--'),
        ('193M', '/mnt/Elements/DNS/film_cooling/C193b_BR15/pictures_job18/stat3_b010', OKABE_ITO[0], '-')]
# diffuser, exit, near and far film, and the whole domain (metering left out, user 2026-11-24)
REGIONS = [('hole_diffuser', 'diffuser'), ('exit', 'exit'),
           ('wake_3_10', 'near film\n' + r'$x/D$ 3–10'), ('wake_10_25', 'far film\n' + r'$x/D$ 10–25'),
           ('domain', 'whole domain')]

set_mpl_params(0.5, grid_on=True)
plt.rcParams['axes.prop_cycle'] = cycler(color=OKABE_ITO)

fig, axs = plt.subplots(1, 5, figsize=(10.0, 2.6), sharex=True, sharey=True)
for n, (ax, (reg, name)) in enumerate(zip(axs, REGIONS)):
    for lab, path, c, ls in RUNS:
        d = np.load(path + '/kolmogorov/length_hist.npz')
        h, e = d[reg + '_Delta_eta'], d['edges_Delta_eta']
        above = 1 - np.concatenate([[0], np.cumsum(h)]) / h.sum()   # volume fraction above the edge
        ax.plot(e, above, color=c, ls=ls, label=lab)
    ax.axvline(2, color='0.5', lw=0.8, ls='--')
    ax.text(0.97, 0.95, f'({"abcde"[n]}) {name}', transform=ax.transAxes, ha='right', va='top')
    ax.set_xlim(0, 4)
    ax.set_ylim(0, 1)
    ax.set_xlabel(r'$\Delta_{mesh}/\eta_{K,193M}$')
axs[0].set_ylabel(r'volume fraction above $\Delta_{mesh}/\eta_{K,193M}$')
fig.tight_layout()
fig.legend(*axs[0].get_legend_handles_labels(), loc='lower center', ncol=4, bbox_to_anchor=(0.5, 0.98))
fig.savefig(OUT, bbox_inches='tight')
print('wrote', OUT)
