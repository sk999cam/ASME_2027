"""Bubble motion and its coherence with the coolant BR, 193M 3DNS (C193b) (paper Fig. F6, Results 4.1).

Usage: python3 fig_bubble_motion.py [SIGNAL]   SIGNAL in V (default), Vp, Vm, Vd
Output: /home/sk999/Documents/ASME_2027/figures/fig_bubble_motion.pdf

Data (read only): /mnt/Elements/DNS/analysis/hole_discharge/bubble_br_long_c193b.npz, summary .txt beside it.
Longest gap-free stretch, 519 records, 103.8 D/U at 0.2 D/U. Welch, Hann, 128 records per segment, 50 % overlap,
7 segments (n_eq = 6.68). 95 % level of zero coherence = level95_eq.
(a) BR (left axis) and bubble signal (right axis, D^3), raw values, against t U/D from the start of the stretch.
(b) Magnitude-squared coherence of BR with the signal, St <= 1 (St = 0 excluded), levels marked where above the 95 % level.
"""
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from cycler import cycler

sys.path.insert(0, '/home/sk999/Documents/3DNS_FP')
from functions import set_mpl_params

SIGNAL = sys.argv[1] if len(sys.argv) > 1 else "V"
SRC = '/mnt/Elements/DNS/analysis/hole_discharge/bubble_br_long_c193b.npz'
OUT = '/home/sk999/Documents/ASME_2027/figures/fig_bubble_motion.pdf'
OKABE_ITO = ['#000000', '#E69F00', '#56B4E9', '#009E73', '#F0E442', '#0072B2', '#D55E00', '#CC79A7']
LABEL = {'V': r'$V$', 'Vp': r'$V^{+}$', 'Vm': r'$V^{-}$', 'Vd': r'$V^{+}-V^{-}$'}

Z = np.load(SRC)
rec, BR, S = Z['rec'], Z['BR'], Z[SIGNAL]
St, coh, lvl = Z['St'], Z['coh_' + SIGNAL], float(Z['level95_eq'])
tu = (rec - rec[0]) * 0.2                       # t U/D, 0.2 D/U per record

m = (St > 0) & (St <= 1.0)
Sm, Cm = St[m], coh[m]
above = Cm > lvl
ip = int(np.argmax(Cm))
print(f'signal {SIGNAL}: BR mean {BR.mean():.4f} std {BR.std():.4f}; signal mean {S.mean():.4f} std {S.std():.4f} D^3')
print(f'peak coherence {Cm[ip]:.3f} at St {Sm[ip]:.3f}; level {lvl:.3f}; bins below St 1: {m.sum()}')
print('bins above level:', ', '.join(f'{s:.3f} ({c:.3f})' for s, c in zip(Sm[above], Cm[above])))
for k in ('V', 'Vp', 'Vm', 'Vd'):
    c = Z['coh_' + k][m]; j = int(np.argmax(c))
    print(f'  {k}: peak {c[j]:.3f} at St {St[m][j]:.3f}')

sf = 0.6
set_mpl_params(sf, grid_on=True)
plt.rcParams['axes.prop_cycle'] = cycler(color=OKABE_ITO)
c1, c2 = '#0072B2', '#D55E00'

fig, (ax, bx) = plt.subplots(2, 1, figsize=(7.0, 7.0), gridspec_kw=dict(height_ratios=[1, 1.1], hspace=0.32))

# (a)
ax.plot(tu, BR, color=c1, ls='-')
ax.set_xlim(0, tu[-1])
ax.set_xlabel(r'$t\,U_\infty/D$')
ax.set_ylabel(r'$BR$', color=c1)
ax.tick_params(axis='y', colors=c1)
ax2 = ax.twinx()
ax2.grid(False)
ax2.plot(tu, S, color=c2, ls='--')
ax2.set_ylabel(LABEL[SIGNAL] + r' [$D^3$]', color=c2)
ax2.tick_params(axis='y', colors=c2)
ax.text(0.01, 0.97, '(a)', transform=ax.transAxes, va='top', ha='left')

# (b)
bx.plot(Sm, Cm, color='k', ls='-')
bx.plot(Sm[above], Cm[above], color='k', ls='none', marker='o')
bx.axhline(lvl, color=c2, ls='--', lw=1.2)
bx.set_xscale('log')
bx.set_xlim(Sm[0] * 0.9, 1.0)
bx.set_ylim(0, 1.0)
bx.set_xlabel(r'$St=fD/U_\infty$')
bx.set_ylabel(r'$\gamma^2_{BR,' + LABEL[SIGNAL].strip('$') + r'}$')
bx.annotate('95 % level', xy=(0.06, lvl), xytext=(0.06, lvl + 0.22), color=c2, ha='center',
            arrowprops=dict(arrowstyle='->', color=c2))
bx.annotate(f'$St$ = {Sm[ip]:.3f}', xy=(Sm[ip], Cm[ip]), xytext=(0.5, 0.97), ha='left', va='top',
            arrowprops=dict(arrowstyle='->', color='k'))
bx.text(0.01, 0.97, '(b)', transform=bx.transAxes, va='top', ha='left')

fig.savefig(OUT, bbox_inches='tight')
print('saved', OUT)
