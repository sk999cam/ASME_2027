"""Motion of the diffuser separation bubble, 193M 3DNS (C193b) at batch b010 (paper Fig. F7, Results 4.1). COMPLETE.

Usage: python3 fig_bubble_motion.py   -> /home/sk999/Documents/ASME_2027/figures/fig_bubble_motion.pdf

Data: /mnt/Tank/DNS/analysis/hole_discharge/bubble3d_c193b.npz (keys t, V, Vp, Vm; written by
analysis/hole_discharge/bubble3d.py, 3-D diffuser field, reversed = u_s < 0, node-volume weighted). 1995 snapshots
every 0.2 D/U over the b010 window 4.413-40.342 ms (99.97 % covered), in 13 gap-free stretches (job restarts).
Coolant monitor /mnt/Elements/DNS/film_cooling/C193b_BR15/monitor_dns.dat (BR = mdot/3.664 g/s, averaged over each
snapshot interval; bubble.coolant_br). Spectra and coherence use bubble.welch_pooled_n as bubble3d.analyse: Welch
segments of 128 records = 25.6 D/U, 50 % overlap, Hann, constant detrend, taken inside each gap-free stretch only
(no segment crosses a gap; stretches shorter than 128 records are skipped) and pooled; dSt = 0.039; nseg segments
in total (printed by the script). Only the 193M is shown. Cached in
/mnt/Elements/DNS/analysis/asme2027_figs/fig_bubble_motion/series_c193b.npz.
(a) Time trace of the separated volume V (reversed-flow volume of the diffuser) and of the two corner lobes
    V(y>0), V(y<0), in D^3, against t U_inf/D (t from the first record). Gaps between stretches are not bridged.
The St axis of (b) and (c) is cut at St = 1: the premultiplied V spectrum holds a negligible share of its variance
    above it (0.01 in St 0.4-1, 0.00 above; b0 summary).
    Second estimate (dotted / dash-dot): Welch segments of 256 records = 51.2 D/U (dSt = 0.0195), 50 % overlap, inside
    the gap-free stretches only. Only n = 4 segments fit (95 % level 1-0.05^(1/3) = 0.63): a coarse check.
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
SRC = '/mnt/Tank/DNS/analysis/hole_discharge/bubble3d_c193b.npz'
RUN = '/mnt/Elements/DNS/film_cooling/C193b_BR15'
CACHE_DIR = '/mnt/Elements/DNS/analysis/asme2027_figs/fig_bubble_motion'
CACHE = CACHE_DIR + '/series_c193b.npz'
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
NL = 256                             # second estimate: 256 records = 51.2 D/U (4 segments fit; >= 4 required)
f2, PVV2, PBB2, PVB2, nseg2 = B.welch_pooled_n(V, br, [tuple(x) for x in st], B.stretches(t)[1], NL)
coh2 = np.abs(PVB2) ** 2 / (PVV2 * PBB2)
lvl2 = 1 - 0.05 ** (1 / (nseg2 - 1))
k2 = 1 + np.argmax((f2 * PVV2)[1:])
i27 = np.argmin(abs(f2 - 0.27))
print(f'long (256): nseg {nseg2}, dSt {f2[1]:.4f}, level {lvl2:.3f}; V peak St {f2[k2]:.3f}; coh at V peak {coh2[k2]:.3f}; '
      f'coh at St {f2[i27]:.3f} {coh2[i27]:.3f}; max coh {coh2[1:].max():.3f} at St {f2[1 + np.argmax(coh2[1:])]:.3f}')
print('long V premult (St, fP/sum):', [(round(a, 3), round(b, 4)) for a, b in zip(f2[1:12], (f2 * PVV2 / PVV2[1:].sum())[1:12])])
tD = (t - t[0]) / TU
print(f'V mean {V.mean():.4f} std {V.std():.4f}; Vp {Vp.mean():.4f}, Vm {Vm.mean():.4f}; BR {br.mean():.4f} +- {br.std():.4f}')
print(f'{len(st)} stretches, {len(t)} records; BR spectrum peak St {f[1 + np.argmax((f * PBB)[1:])]:.3f}')
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
ax.set_ylim(0.1, 1.1)
ax.text(0.015, 0.86, '(a)', transform=ax.transAxes)

ax = axs[1]
m = f > 0
m2 = f2 > 0
r = len(f2) / len(f)                 # same area per unit log St is not needed; scale so both sum to 1 per bin count
ax.semilogx(f[m], f[m] * PVV[m] / PVV[m].sum(), color='#D55E00', ls='-', label='$V$, 25.6 D/U')
ax.semilogx(f[m], f[m] * PBB[m] / PBB[m].sum(), color='#0072B2', ls='--', label='BR, 25.6 D/U')
ax.semilogx(f2[m2], f2[m2] * PVV2[m2] / PVV2[m2].sum() * r, color='#D55E00', ls=':', label='$V$, 51.2 D/U')
ax.semilogx(f2[m2], f2[m2] * PBB2[m2] / PBB2[m2].sum() * r, color='#0072B2', ls='-.', label='BR, 51.2 D/U')
ax.set_xlabel(r'$St=fD/U_\infty$')
ax.set_ylabel(r'$St\,P(St)\,/\sum P$')
ax.set_xlim(f[1], 1.0)
ax.set_ylim(0, None)
ax.legend(loc='upper left', bbox_to_anchor=(0.1, 1.0), fontsize='small')
ax.text(0.015, 0.86, '(b)', transform=ax.transAxes)

ax = axs[2]
ax.semilogx(f[m], coh[m], color='#000000', ls='-', label='25.6 D/U')
ax.semilogx(f2[m2], coh2[m2], color='#000000', ls=':', label='51.2 D/U')
ax.axhline(lvl, color='0.4', ls='-', lw=0.8, label=f'95 % level, $n={nseg}$')
ax.axhline(lvl2, color='0.4', ls=':', lw=0.8, label=f'95 % level, $n={nseg2}$')
ax.legend(loc='lower center', bbox_to_anchor=(0.5, 1.0), ncol=2, fontsize='small', columnspacing=1.0)
ax.set_xlabel(r'$St=fD/U_\infty$')
ax.set_ylabel(r'coherence of $V$ and BR')
ax.set_xlim(f[1], 1.0)
ax.set_ylim(0, 1)
ax.text(0.015, 0.86, '(c)', transform=ax.transAxes)
fig.tight_layout()
fig.savefig(OUT, bbox_inches='tight')
print('wrote', OUT)
