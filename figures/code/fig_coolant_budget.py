"""fig_coolant_budget: circulation and coolant budget of the film, RANS 50M against 193M 3DNS (b010), x/D 4-20, four panels (row or 2 x 2, see LAYOUT).

Usage: python3 fig_coolant_budget.py   (writes figures/fig_coolant_budget.pdf and planning/review/44_fig_coolant_budget.png)

Merges panels of three scripts; the data handling is imported from them and they are not changed:
(a) circulation of the counter-rotating vortex pair (solid) and of the wall vortex (dash-dot): fig_mean_features.run() (fine cross-planes,
    mean of the two sides, band = between the two sides for the 193M).
(b) vertical coolant flux through z = z_x(x), mean-flow part (solid) and turbulent part (dashed): fig_flux_eta.SURF (keys Jm_w, Jt_w), band = +-1 SE of the 193M.
(c) swap test of the mean-flow flux J[w_A, c_B] through z_x: fig_swap_lateral.NPZ. Colour = source of the velocity, solid = velocity and coolant
    from the same run, dashed = swapped. Band = +-1 SE of the 193M-containing combinations. Linear axis cut at 35 g s^-1 m^-1.
(d) lateral coolant flux through |y| = 0.5 D below z_x, mean-flow (solid) and turbulent (dashed), + away from the centreline; band = +-1 SE of the 193M.
Colour convention: run by colour (193M black, RANS 50M blue), quantity by line style.
"""
import os
import sys
import warnings

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from cycler import cycler
from matplotlib.lines import Line2D

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, '/home/sk999/Documents/3DNS_FP')
from functions import set_mpl_params  # noqa: E402
import fig_mean_features as MF  # noqa: E402
import fig_flux_eta as FE  # noqa: E402
import fig_swap_lateral as SL  # noqa: E402

LAYOUT = '2x2'   # '2x2' (paper, owner 2026-11-27): 2 x 2 grid, 7.5 in wide (figure*, \textwidth); 'row': four near-square panels, 7.5 in wide. The chosen one goes to OUT.
OUT = '/home/sk999/Documents/ASME_2027/figures/fig_coolant_budget.pdf'
PNG = {'row': '/home/sk999/Documents/ASME_2027/planning/review/46_fig_coolant_budget_row.png',
       '2x2': '/home/sk999/Documents/ASME_2027/planning/review/46_fig_coolant_budget_2x2.png'}
SIZE = {'row': (7.5, 2.0), '2x2': (7.5, 6.5)}
OKABE_ITO = ['#000000', '#E69F00', '#56B4E9', '#009E73', '#F0E442', '#0072B2', '#D55E00', '#CC79A7']
C193, CRANS = OKABE_ITO[0], OKABE_ITO[5]
XMIN = 4.0           # x/D 4-20 on all panels (owner 2026-11-27: z_x confirmed from x/D 4)
YLIM_C = (0, 50)     # linear; values at x/D 4 are 22-24 (193M) and 45 (RANS mean flux)


def circulation(out, ST):
    """CRVP and wall-vortex circulation per run: (x, mean, lo, hi, open-marker flag) from the fine planes."""
    W = out['wvp']
    res = {}
    for nm in MF.NAMES:
        g = np.array([[W[nm][xd][sd]['crvp']['g'] for xd in ST] for sd in 'LR'])
        wx, wm, wlo, whi, wop = [], [], [], [], []
        for xd in ST:
            w = W[nm][xd]
            gs = [w[sd]['g'] for sd in 'LR' if not np.isnan(w[sd]['g'])]
            if not gs:
                continue
            wx.append(xd)
            wm.append(np.mean(gs))
            wlo.append(min(gs))
            whi.append(max(gs))
            wop.append(any(w[sd]['below3'] or not w[sd]['sw'] for sd in 'LR' if not np.isnan(w[sd]['g'])))
        res[nm] = dict(crvp=(ST, g.mean(0), g.min(0), g.max(0)), wv=tuple(map(np.array, (wx, wm, wlo, whi, wop))))
    return res


def panel_a(ax, out):
    ST = np.array([s for s in MF.STATIONS if s >= XMIN])
    R = circulation(out, ST)
    for nm, col in (('dns', C193), ('rans', CRANS)):
        x, m, lo, hi = R[nm]['crvp']
        ax.plot(x, m, '-', color=col, marker='o', ms=2.5)
        if nm == 'dns':
            ax.fill_between(x, lo, hi, color=col, alpha=0.2, lw=0)
        x, m, lo, hi, op = R[nm]['wv']
        ax.plot(x, m, '-.', color=col, lw=0.8)
        for xv, gv, o_ in zip(x, m, op):
            ax.plot([xv], [gv], 's', color=col, mfc='none' if o_ else col, ms=3)
        if nm == 'dns':
            ax.fill_between(x, lo, hi, color=col, alpha=0.2, lw=0)
    ax.set_ylabel(r'$\Gamma/(U_\infty D)$')
    ax.set_ylim(0, 0.2)
    ax.text(0.97, 0.62, 'CRVP', transform=ax.transAxes, ha='right', va='center')
    ax.text(9.5, 0.047, 'wall vortex', ha='center', va='center')


def panel_b(ax):
    Z0 = np.load(FE.SURF)
    m = Z0['x_D'] >= XMIN
    Z = {k: Z0[k][m] for k in Z0.files if Z0[k].ndim == 1}
    x = Z['x_D']
    for key, ls in ((FE.KEY_M, '-'), (FE.KEY_T, '--')):
        j, se, jr = Z['193_' + key], Z['193_' + key + '_se'], Z['rans_' + key]
        ax.plot(x, j, ls, color=C193)
        ax.plot(x, jr, ls, color=CRANS)
        ax.fill_between(x, np.clip(j - se, 1e-3, None), np.clip(j + se, 1e-3, None), color=C193, alpha=0.25, lw=0)
    ax.set_ylim(*YLIM_C)     # linear, same limits as panel (c)
    ax.set_ylabel(r'$J_w$ / g s$^{-1}$m$^{-1}$')


def panel_c(ax, d):
    x = d['x_D']
    combos = [('193', '193', C193, '-', r'$w$: 193, $c$: 193'),
              ('RANS', '193', CRANS, '--', r'$w$: R, $c$: 193'),
              ('193', 'RANS', C193, '--', r'$w$: 193, $c$: R'),
              ('RANS', 'RANS', CRANS, '-', r'$w$: R, $c$: R')]
    for a, b, col, ls, lab in combos:
        j = d['swap_w%s_c%s_J' % (a, b)]
        se = d['swap_w%s_c%s_J_se' % (a, b)]
        ax.plot(x, j, ls, color=col, label=lab)
        if np.any(se > 0):
            ax.fill_between(x, np.clip(j - se, 1e-3, None), j + se, color=col, alpha=0.25, lw=0)
    ax.set_ylim(*YLIM_C)
    ax.set_ylabel(r'$J_w$ / g s$^{-1}$m$^{-1}$')
    ax.legend(loc='upper right', fontsize=plt.rcParams['font.size'] - 1, handlelength=1.2, labelspacing=0.25, borderpad=0.2,
              handletextpad=0.5, borderaxespad=0.2)


def panel_d(ax, d):
    x = d['x_D']
    ax.axhline(0.0, color='gray', lw=0.7)
    for r, col in (('193', C193), ('rans', CRANS)):
        for k, ls in (('mean', '-'), ('turb', '--')):
            v = d['lat_zx_%s_%s_ye0.5' % (r, k)]
            ax.plot(x, v, ls, color=col)
            if r == '193':
                se = d['lat_zx_193_%s_ye0.5_se' % k]
                ax.fill_between(x, v - se, v + se, color=col, alpha=0.25, lw=0)
    ax.set_ylabel(r'$J_y$ / g s$^{-1}$m$^{-1}$')
    fs = plt.rcParams['font.size'] - 1
    ax.set_ylim(-14.5, 8)
    ax.text(19.8, 4.2, r'$\uparrow$ outward', ha='right', va='center', fontsize=fs, color='0.35')
    ax.text(19.8, -5.0, r'$\downarrow$ inward', ha='right', va='center', fontsize=fs, color='0.35')


def draw(layout, out, d, pdf=None):
    # Drawn at the printed size (row: 7.5 in, figure* width=\textwidth; 2x2: 5.6 in, figure* width=5.6in or 0.75\textwidth), scale 1; body text 10 pt.
    set_mpl_params(10 / 22, grid_on=True)
    plt.rcParams['axes.prop_cycle'] = cycler(color=OKABE_ITO)
    size = SIZE[layout]
    if layout == 'row':
        fig, axs = plt.subplots(1, 4, figsize=size)
    else:
        fig, axs = plt.subplots(2, 2, figsize=size)
        axs = axs.ravel()
    panel_a(axs[0], out)
    panel_b(axs[1])
    panel_c(axs[2], d)
    panel_d(axs[3], d)
    for i, (a, t) in enumerate(zip(axs, 'abcd')):
        a.set_xlim(XMIN, 20)
        a.set_xticks([5, 10, 15, 20])
        if layout == 'row' or i >= 2:
            a.set_xlabel(r'$x/D$')
        a.text(0.0, 1.02, '(%s)' % t, transform=a.transAxes, ha='left', va='bottom')
    fs = plt.rcParams['font.size'] - 1
    fig.tight_layout(pad=0.3, w_pad=0.5, h_pad=0.6, rect=(0, 0.1 if layout == 'row' else 0.06, 1, 1))
    handles = [Line2D([], [], color=C193, ls='-', label='193M 3DNS'), Line2D([], [], color=CRANS, ls='-', label='RANS 50M'),
               Line2D([], [], color='gray', ls='-', label='mean-flow (b, d)'), Line2D([], [], color='gray', ls='--', label='turbulent (b, d)')]
    fig.legend(handles=handles, loc='lower center', ncol=4, frameon=False, fontsize=fs, handlelength=1.8, columnspacing=1.5,
               bbox_to_anchor=(0.5, 0.0))
    print(layout, 'figsize', size, 'axes box w x h in:', [(round(a.get_position().width * size[0], 2), round(a.get_position().height * size[1], 2)) for a in axs])
    if pdf:
        fig.savefig(pdf)
    fig.savefig(PNG[layout], dpi=200)
    plt.close(fig)


def main():
    warnings.filterwarnings('ignore')
    out = MF.run()
    d = np.load(SL.NPZ)
    draw(LAYOUT, out, d, OUT)
    draw('2x2' if LAYOUT == 'row' else 'row', out, d)     # preview of the other layout only
    print('wrote', OUT)


if __name__ == '__main__':
    main()
