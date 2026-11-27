"""Table of point-wise mu_eff distributions with strain filtering (F14 replacement), and the s0 scan figure.
Usage: python3 meff_pointwise_table.py
Writes data/meff_pointwise_table.txt and the scratchpad figure meff_strain_filter_scan.pdf. Data and mask: _meff_pointwise_data.py."""
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from cycler import cycler

sys.path.insert(0, '/home/sk999/Documents/3DNS_FP')
from functions import set_mpl_params
import _meff_pointwise_data as M

OUT = '/home/sk999/Documents/ASME_2027/figures/code/data/meff_pointwise_table.txt'
FIG = M.SCR + 'meff_strain_filter_scan.pdf'
OLD = {'near_jet_0_3': -0.328, 'wake_3_10': -0.390, 'wake_10_25': -0.855}
S0S = [0.0, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0]
O = M.load()
L = []
P = L.append
P('# Point-wise effective eddy viscosity, 193M (batch b010) against RANS 50M SST; replacement data of F14')
P('# Definitions: mu_eff = -a:S*/(2|S*|^2) (Wu, Xiao, Paterson, PRF 3, 074602, 2018), a = deviatoric Reynolds stress, S* = deviatoric mean strain.')
P('#   193M: mu_eff = rho k^1/2 L_eff cos (signed). RANS: mu_t = rho k^1/2 L_R (cos = 1 by construction, checked: min = max = 1). mu normalised by rho_inf U_inf D')
P('#   (rho_inf = 0.88918 kg/m3 = p_exit/(R T_inf), U_inf = 86.2812 m/s, D = 7.75 mm). L = mu/(rho k^1/2)/D (193M: L_al = mu_eff/(rho k^1/2), signed). k/U_inf^2. |S*| D/U_inf.')
P('# Mask: mode C with x/D <= 20 (estimator.py): finite values, k_193 > 0 and k_RANS > 0 (wall nodes excluded), z/D >= 0.3, jet x/D 0-3, near film 3-10, far film 10-20.')
P('# Weights: cell volume V of the 50M grid (3.1 million point subset pts.npz). Plain volume-weighted quantiles (midpoint rule), negative values kept.')
P('# Strain filter: only cells with |S*| D/U_inf > s0, |S*| of the 193M (variant "193"), or of both rungs (variant "both"). s0 = 0 is the unfiltered set (mask above).')
P('# kept = kept volume fraction of the masked region; neg = volume fraction with mu_eff,193M <= 0 among kept cells; r = mu_t,RANS/mu_eff,193M over kept cells with mu_eff,193M > 0')
P('#   (rkept = their volume fraction of the kept cells); q25/q50/q75 are volume-weighted; mu in rho_inf U_inf D units; kR,kD = median k/U_inf^2; LR,LD = median L/D; cosD = median 193M cos;')
P('#   ratio_med = median(mu_t)/median(mu_eff) over all kept cells (negatives included).')
P('# Old regional strain-squared-weighted LS ratio mu_RANS/mu_193 (estimator_modeC_x20.log): jet 0.720, near film 0.677, far film 0.425 (ln -0.328, -0.390, -0.855).')
P('')
hd = ('variant s0 region kept neg | R_q25 R_q50 R_q75 | D_q25 D_q50 D_q75 | rkept r_q25 r_q50 r_q75 | lnr_q25 lnr_q50 lnr_q75 | ratio_med | kR kD | LR LD | cosD')
P('# columns: ' + hd)
SCAN = {}
for var in ('193', 'both'):
    for s0 in S0S:
        for key, _, lab in M.REGIONS:
            o = O[key]
            m = M.keep(o, s0, var == 'both')
            tot = o['V'].sum()
            if m.sum() < 50:
                P(f'{var} {s0:.2f} {lab.replace(" ", "_")} no cells')
                continue
            V = o['V'][m]
            kept = V.sum() / tot
            muR, muD = o['muR'][m], o['muD'][m]
            neg = V[muD <= 0].sum() / V.sum()
            qR, qD = M.wq(muR, V, (.25, .5, .75)), M.wq(muD, V, (.25, .5, .75))
            p = muD > 0
            r = muR[p] / muD[p]
            qr, ql = M.wq(r, V[p], (.25, .5, .75)), M.wq(np.log(r), V[p], (.25, .5, .75))
            md = lambda a: M.wq(a[m], V, (.5,))[0]  # noqa: E731
            rm = qR[1] / qD[1] if qD[1] > 0 else float('nan')
            SCAN[(var, s0, key)] = dict(kept=kept, medr=qr[1], rm=rm, neg=neg)
            P(f'{var} {s0:.2f} {lab.replace(" ", "_")} {kept:.4f} {neg:.4f} | {qR[0]:.4e} {qR[1]:.4e} {qR[2]:.4e} | {qD[0]:+.4e} {qD[1]:+.4e} {qD[2]:+.4e} | '
              f'{V[p].sum() / V.sum():.4f} {qr[0]:.4f} {qr[1]:.4f} {qr[2]:.4f} | {ql[0]:+.4f} {ql[1]:+.4f} {ql[2]:+.4f} | {rm:.4f} | '
              f'{md(o["kR"]):.4e} {md(o["kD"]):.4e} | {md(o["LR"]):.4e} {md(o["LD"]):+.4e} | {md(o["cD"]):+.4f}')
P('')
P('## Regional 193M |S*| D/U_inf quantiles (q50 q90 q99 max) of the masked cells, for the threshold options')
for key, _, lab in M.REGIONS:
    o = O[key]
    q = M.wq(o['SD'], o['V'], (.5, .9, .99))
    P(f'{lab:20s} {q[0]:.4f} {q[1]:.4f} {q[2]:.4f} {o["SD"].max():.4f}')
P('')
P('## s0 at which the median r crosses 1 (variant 193, linear interpolation on the scan; "none" if no crossing)')
for key, _, lab in M.REGIONS:
    xs = [s for s in S0S if ('193', s, key) in SCAN]
    ys = [SCAN[('193', s, key)]['medr'] for s in xs]
    cl = []
    for i in range(len(xs) - 1):
        if (ys[i] - 1) * (ys[i + 1] - 1) < 0:
            cl.append(f'{xs[i] + (1 - ys[i]) / (ys[i + 1] - ys[i]) * (xs[i + 1] - xs[i]):.3f} (between {xs[i]} and {xs[i + 1]}: {ys[i]:.3f}, {ys[i + 1]:.3f})')
    c = '; '.join(cl) if cl else 'none'
    P(f'{lab:20s} median r: {c}   values: ' + ' '.join(f'{y:.3f}' for y in ys))
open(OUT, 'w', encoding='utf-8').write('\n'.join(L) + '\n')
print('\n'.join(L))

# ---------------------------------------------------------------- scan figure
set_mpl_params(0.4, grid_on=True)
plt.rcParams['axes.prop_cycle'] = cycler(color=['#000000', '#E69F00', '#56B4E9', '#009E73', '#F0E442', '#0072B2', '#D55E00', '#CC79A7'])
sc = [s for s in S0S if s > 0]
xs = np.array([0.01] + sc)    # s0 = 0 drawn at 0.01 on the log axis
fig, axs = plt.subplots(1, 3, figsize=(7.2, 2.7), sharey=True, gridspec_kw=dict(wspace=0.08))
for ax, (key, _, lab), tl, old in zip(axs, M.REGIONS, ['jet $x/D$ 0–3', 'near film $x/D$ 3–10', 'far film $x/D$ 10–20'],
                                      [np.exp(OLD[k]) for k in OLD]):
    ss = [s for s in S0S if ('193', s, key) in SCAN]
    x = np.array([0.01 if s == 0 else s for s in ss])
    ax.plot(x, [SCAN[('193', s, key)]['medr'] for s in ss], '-o', color='#0072B2', ms=3, label=r'median $r$')
    ax.plot(x, [SCAN[('193', s, key)]['rm'] for s in ss], '--s', color='#D55E00', ms=3, label='ratio of medians')
    ax.plot([x[0], x[-1]], [1, 1], ':', color='#000000', lw=0.7)
    ax.axhline(old, color='#009E73', ls='-.', lw=0.9, label='regional weighted value')
    for xi, s in zip(x, ss):
        ax.text(xi, 0.03, f'{100 * SCAN[("193", s, key)]["kept"]:.0f}', ha='center', va='bottom', fontsize=5.5, transform=ax.get_xaxis_transform())
    ax.set_xscale('log')
    ax.text(0.04, 0.96, tl, transform=ax.transAxes, ha='left', va='top', fontsize=6.5)
    ax.set_xlabel(r'$s_0$ (filter $|S^*|D/U_\infty>s_0$)')
axs[0].set_ylabel(r'$\mu_{t,\mathrm{RANS}}/\mu_{\mathrm{eff},193}$')
axs[0].set_ylim(0, 11)
axs[0].text(0.03, 0.10, 'kept volume (%) under points;\nleftmost point: $s_0=0$', transform=axs[0].transAxes, fontsize=5.5, va='bottom')
axs[2].legend(loc='upper right', fontsize=5.5)
fig.savefig(FIG, bbox_inches='tight')
print('saved', FIG)
