"""Motion of the diffuser separation bubble, 50M 3DNS only (paper Fig. F7, Results 4.1). INTERIM.

Usage: python3 fig_bubble_motion.py   -> /home/sk999/Documents/ASME_2027/figures/fig_bubble_motion.pdf

Only the 50M 3DNS (C50b) has the snapshot data (vol_nearwall records every 0.2 D/U, 1199 records up to the end of
stat3 batch b006, 4 gap-free stretches); the other rungs are not shown. Replace by the 193M 3DNS when snapshots exist.
Data: /mnt/Tank/DNS/analysis/hole_discharge/bubble3d_c50b.npz (keys t, V, Vp, Vm; written by
analysis/hole_discharge/bubble3d.py, 3-D diffuser field, reversed = u_s < 0, node-volume weighted) and the coolant
monitor /mnt/Tank/DNS/film_cooling/C50b_BR15/monitor_dns.dat (BR = mdot/3.664 g/s, averaged over each snapshot
interval; bubble.coolant_br). Spectra and coherence use bubble.welch_pooled_n as bubble3d.analyse (Welch segments of
128 records = 25.6 D/U, 50 % overlap, Hann, pooled over the stretches; dSt = 0.039; 12 independent-ish segments pooled to
nseg below). Cached in /mnt/Elements/DNS/analysis/asme2027_figs/fig_bubble_motion/series.npz.
(a) Time trace of the separated volume V (reversed-flow volume of the diffuser) and of the two corner lobes
    V(y>0), V(y<0), in D^3, against t U_inf/D (t from the first record). Gaps between stretches are not bridged.
(b) Premultiplied spectra St*P(St), each normalised by its variance, of V and of BR, against St = f D / U_inf.
(c) Magnitude-squared coherence of V and BR; dashed line = 95 % significance level 1 - 0.05^(1/(nseg-1)).
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from cycler import cycler

sys.path.insert(0, '/home/sk999/Documents/3DNS_FP')
from functions import set_mpl_params
sys.path.insert(0, '/home/sk999/Documents/4dns_filmcooled/analysis/hole_discharge')
import bubble as B  # noqa: E402

OKABE_ITO = ['#000000', '#E69F00', '#56B4E9', '#009E73', '#F0E442', '#0072B2', '#D55E00', '#CC79A7']
SRC = '/mnt/Tank/DNS/analysis/hole_discharge/bubble3d_c50b.npz'
RUN = '/mnt/Tank/DNS/film_cooling/C50b_BR15'
CACHE_DIR = '/mnt/Elements/DNS/analysis/asme2027_figs/fig_bubble_motion'
CACHE = CACHE_DIR + '/series.npz'
OUT = '/home/sk999/Documents/ASME_2027/figures/fig_bubble_motion.pdf'
NPER = 128
TU = B.TU

if not os.path.exists(CACHE):
    os.makedirs(CACHE_DIR, exist_ok=True)
    Z = np.load(SRC)
    t, V, Vp, Vm = Z['t'], Z['V'], Z['Vp'], Z['Vm']
    st, dt = B.stretches(t)
    br, _, _ = B.coolant_br(RUN, t, dt)
    f, PVV, PBB, PVB, nseg = B.welch_pooled_n(V, br, st, dt, NPER)
    np.savez(CACHE, t=t, V=V, Vp=Vp, Vm=Vm, br=br, st=np.array(st), dt=dt, f=f, PVV=PVV, PBB=PBB, PVB=PVB, nseg=nseg)
C = np.load(CACHE)
t, V, Vp, Vm, br, st = C['t'], C['V'], C['Vp'], C['Vm'], C['br'], C['st']
f, PVV, PBB, PVB, nseg = C['f'], C['PVV'], C['PBB'], C['PVB'], int(C['nseg'])
coh = np.abs(PVB) ** 2 / (PVV * PBB)
lvl = 1 - 0.05 ** (1 / (nseg - 1))
tD = (t - t[0]) / TU
print(f'V mean {V.mean():.4f} std {V.std():.4f}; Vp {Vp.mean():.4f}, Vm {Vm.mean():.4f}; BR {br.mean():.4f} +- {br.std():.4f}')
k = 1 + np.argmax((f * PVV)[1:])
print(f'spectrum peak St {f[k]:.3f}; coh at peak {coh[k]:.3f}; max coh {coh[1:].max():.3f} at St {f[1 + np.argmax(coh[1:])]:.3f};'
      f' 95% level {lvl:.3f}; nseg {nseg}; corr lobes {np.corrcoef(Vp, Vm)[0, 1]:.3f}')

set_mpl_params(0.5, grid_on=True)
plt.rcParams['axes.prop_cycle'] = cycler(color=OKABE_ITO)

fig, axs = plt.subplots(3, 1, figsize=(5.0, 7.0))
ax = axs[0]
for a, b in st:                      # one line per gap-free stretch: no line across a gap
    ax.plot(tD[a:b], V[a:b], color='#000000', ls='-', label='$V$' if a == st[0][0] else None)
    ax.plot(tD[a:b], Vp[a:b], color='#D55E00', ls='--', label=r'$V(y>0)$' if a == st[0][0] else None)
    ax.plot(tD[a:b], Vm[a:b], color='#0072B2', ls='-.', label=r'$V(y<0)$' if a == st[0][0] else None)
ax.set_xlabel(r'$tU_\infty/D$')
ax.set_ylabel(r'separated volume ($D^3$)')
ax.set_xlim(0, tD[-1])
ax.legend(loc='upper center', ncol=3, bbox_to_anchor=(0.5, 1.2), columnspacing=1.0, handlelength=1.6)
ax.text(0.015, 0.06, '(a)', transform=ax.transAxes)

ax = axs[1]
m = f > 0
ax.semilogx(f[m], f[m] * PVV[m] / PVV[m].sum(), color='#D55E00', ls='-', label='$V$')
ax.semilogx(f[m], f[m] * PBB[m] / PBB[m].sum(), color='#0072B2', ls='--', label='BR')
ax.set_xlabel(r'$St=fD/U_\infty$')
ax.set_ylabel(r'$St\,P(St)\,/\sum P$')
ax.set_xlim(f[1], f[-1])
ax.legend(loc='upper right')
ax.text(0.015, 0.86, '(b)', transform=ax.transAxes)

ax = axs[2]
ax.semilogx(f[m], coh[m], color='#000000', ls='-')
ax.axhline(lvl, color='0.4', ls='--', lw=0.8)
ax.text(0.40, lvl + 0.02, '95 % level', transform=ax.get_yaxis_transform(), ha='left', va='bottom')
ax.set_xlabel(r'$St=fD/U_\infty$')
ax.set_ylabel(r'coherence of $V$ and BR')
ax.set_xlim(f[1], f[-1])
ax.set_ylim(0, 1)
ax.text(0.015, 0.86, '(c)', transform=ax.transAxes)
fig.tight_layout()
fig.savefig(OUT, bbox_inches='tight')
print('wrote', OUT)
