"""F14 replacement, candidate 1: per-region distributions of the point-wise effective eddy viscosity.

Usage: python3 fig_meff_pointwise.py     (table: meff_pointwise_table.py -> data/meff_pointwise_table.txt)
Output: /home/sk999/Documents/ASME_2027/figures/fig_meff_pointwise.pdf

Estimator: mu_eff = -a:S*/(2|S*|^2) (Wu, Xiao, Paterson 2018), point-wise; mask and weights in _meff_pointwise_data.py
(mode C, x/D <= 20, volume weights, wall nodes excluded, z/D >= 0.3).
(a) mu/(rho_inf U_inf D) for RANS (mu_t) and 193M b010 (mu_eff, signed): box q25-q75, line median, whiskers q10-q90.
    Linear axis; the volume fraction with mu_eff <= 0 is written under each 193M box.
(b) ln r, r = mu_t,RANS/mu_eff,193M at the same point, over the cells with mu_eff,193M > 0 (kept volume fraction written
    under each box); same box convention.
Style: set_mpl_params, Okabe-Ito, no title.
"""
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from cycler import cycler
from matplotlib.patches import Patch

sys.path.insert(0, '/home/sk999/Documents/3DNS_FP')
from functions import set_mpl_params
import _meff_pointwise_data as M

S0 = float(sys.argv[1]) if len(sys.argv) > 1 else 0.05     # strain filter |S*| D/U_inf > S0 (193M); 0 = no filter
OUT = sys.argv[2] if len(sys.argv) > 2 else '/home/sk999/Documents/ASME_2027/figures/fig_meff_pointwise.pdf'
OKABE_ITO = ['#000000', '#E69F00', '#56B4E9', '#009E73', '#F0E442', '#0072B2', '#D55E00', '#CC79A7']
C_R, C_D, C_RATIO = '#E69F00', '#0072B2', '#CC79A7'
O = M.load()
KEPT = {}
for _k, _, _ in M.REGIONS:
    _m = M.keep(O[_k], S0)
    KEPT[_k] = O[_k]['V'][_m].sum() / O[_k]['V'].sum()
    O[_k] = {a: (b[_m] if isinstance(b, np.ndarray) else b) for a, b in O[_k].items()}
SC = 1e3


def stat(a, w, label=''):
    q = M.wq(a, w)
    return dict(whislo=q[0], q1=q[1], med=q[2], q3=q[3], whishi=q[4], fliers=[], label=label)


def draw(ax, stats, pos, color, hatch, width):
    bp = ax.bxp(stats, positions=pos, widths=width, showfliers=False, patch_artist=True, manage_ticks=False,
                boxprops=dict(facecolor=color, edgecolor='k', linewidth=0.6, alpha=0.9),
                medianprops=dict(color='k', linewidth=1.4), whiskerprops=dict(linewidth=0.7),
                capprops=dict(linewidth=0.7), capwidths=width * 0.5)
    for b in bp['boxes']:
        b.set_hatch(hatch)


set_mpl_params(0.4, grid_on=True)
plt.rcParams['axes.prop_cycle'] = cycler(color=OKABE_ITO)
fig, (axa, axb) = plt.subplots(1, 2, figsize=(7.2, 3.2), gridspec_kw=dict(wspace=0.28))
pos = np.arange(3) * 1.0
w = 0.3
sR = [stat(O[k]['muR'] * SC, O[k]['V']) for k, _, _ in M.REGIONS]
sD = [stat(O[k]['muD'] * SC, O[k]['V']) for k, _, _ in M.REGIONS]
draw(axa, sR, pos - 0.19, C_R, '', w)
draw(axa, sD, pos + 0.19, C_D, '////', w)
ymin = min(s['whislo'] for s in sD)
for i, (k, _, _) in enumerate(M.REGIONS):
    o = O[k]
    f = o['V'][o['muD'] <= 0].sum() / o['V'].sum()
    axa.text(pos[i], min(-0.8, ymin - 0.4), f'kept {100 * KEPT[k]:.0f} %\n' r'$\mu_\mathrm{eff}\leq0$: ' f'{100 * f:.0f} %',
             ha='center', va='top', fontsize=5.5)
axa.axhline(0, color='k', lw=0.6)
axa.set_ylim(-5.6, 20)
axa.set_ylabel(r'$10^{3}\,\mu/(\rho_\infty U_\infty D)$')
axa.legend(handles=[Patch(facecolor=C_R, edgecolor='k', linewidth=0.6, label=r'RANS $\mu_t$'),
                    Patch(facecolor=C_D, edgecolor='k', linewidth=0.6, hatch='////', label=r'193M $\mu_\mathrm{eff}$')],
           loc='upper right', ncol=2, fontsize=6, columnspacing=0.8, handlelength=1.4)
# (b)
sr, kept = [], []
for k, _, _ in M.REGIONS:
    o = O[k]
    m = o['muD'] > 0
    sr.append(stat(np.log(o['muR'][m] / o['muD'][m]), o['V'][m]))
    kept.append(o['V'][m].sum() / o['V'].sum())
draw(axb, sr, pos, C_RATIO, '\\\\\\\\', 0.42)
axb.axhline(0, color='k', lw=0.6)
axb.set_ylim(-2.6, 4.6)
for i in range(3):
    axb.text(pos[i], sr[i]['whislo'] - 0.15, f'{100 * kept[i]:.0f} % of kept cells', ha='center', va='top', fontsize=5.5)
axb.set_ylabel(r'$\ln(\mu_{t,\mathrm{RANS}}/\mu_{\mathrm{eff},193})$')
axb.annotate('RANS larger', xy=(-0.66, 0.1), xytext=(-0.66, 1.5), ha='center', fontsize=5.5,
             arrowprops=dict(arrowstyle='->', lw=0.6))
axb.annotate('RANS smaller', xy=(-0.66, -0.1), xytext=(-0.66, -1.5), ha='center', fontsize=5.5,
             arrowprops=dict(arrowstyle='->', lw=0.6))
for ax, letter in ((axa, '(a)'), (axb, '(b)')):
    ax.set_xticks(pos)
    ax.set_xticklabels(['jet\n' r'$x/D$ 0–3', 'near film\n' r'$x/D$ 3–10', 'far film\n' r'$x/D$ 10–20'], fontsize=6.5)
    ax.set_xlim(-0.95, 2.6)
    ax.text(0.015, 0.985, letter, transform=ax.transAxes, ha='left', va='top', fontsize=8)
fig.savefig(OUT, bbox_inches='tight')
print('saved', OUT)
