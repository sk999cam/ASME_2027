"""fig_swap_lateral: swap test of the vertical coolant flux through z_x and lateral coolant flux below z_x, RANS 50M against 193M 3DNS (b010), x/D 3-20.

Usage: python3 fig_swap_lateral.py   (writes figures/fig_swap_lateral.pdf and figures/code/data/swap_lateral_table.txt)

Source: /mnt/Elements/DNS/analysis/asme2027_figs/wall_balance/n7_swap_lateral.npz (read only; see README_n7_swap_lateral.txt there).
(a) J[w_A, c_B] = pitch-integrated <rho w>_A c~_B through z = z_x(x), g s^-1 m^-1, log axis. Colour = velocity source (193M black, RANS blue),
    solid = velocity and coolant from the same run, dashed = swapped. Band = +-1 jackknife SE of the 193M-containing combinations (10 batches).
(b) Lateral flux through the planes |y| = 0.5 D from the wall to z_x, + away from the centreline. Colour = run, solid = mean-flow part,
    dashed = turbulent part. Band = +-1 SE of the 193M parts.
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from cycler import cycler
from matplotlib.lines import Line2D

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, '/home/sk999/Documents/3DNS_FP')
from functions import set_mpl_params  # noqa: E402

NPZ = '/mnt/Elements/DNS/analysis/asme2027_figs/wall_balance/n7_swap_lateral.npz'
OUT = '/home/sk999/Documents/ASME_2027/figures/fig_swap_lateral.pdf'
TAB = os.path.join(HERE, 'data', 'swap_lateral_table.txt')
OKABE_ITO = ['#000000', '#E69F00', '#56B4E9', '#009E73', '#F0E442', '#0072B2', '#D55E00', '#CC79A7']
C193, CRANS = OKABE_ITO[0], OKABE_ITO[5]
YSCALE_A = 'linear'   # 'linear' or 'log': y axis of panel (a); the PDF in figures/ is drawn with this setting
YLIM_LIN = (0, 35)    # capped linear axis, lines run off the top near x/D 3-4
SCR = '/tmp/claude-1001/-home-sk999-Documents-ASME-2027/f4c92ceb-6df8-4254-96c3-6552d15d0ecf/scratchpad/'
REV = '/home/sk999/Documents/ASME_2027/planning/review/'


def table(d):
    x = d['x_D']
    J = lambda a, b: d['swap_w%s_c%s_J' % (a, b)]
    j0, jw, jc, jr = J('193', '193'), J('RANS', '193'), J('193', 'RANS'), J('RANS', 'RANS')
    L = lambda r, k, y: d['lat_zx_%s_%s_ye%s' % (r, k, y)]
    st = [3, 5, 7, 10, 15, 20]
    ix = [int(np.argmin(abs(x - s))) for s in st]
    o = ['# fig_swap_lateral: plotted values and derived ratios (n7_swap_lateral.npz, g s^-1 m^-1, per metre of x)',
         '# SE = 193M jackknife standard error (10 batches)', '',
         '## Plotted values',
         '#  x/D   J[193,193]  (SE)   J[R,193]  (SE)   J[193,R]  (SE)   J[R,R]']
    for i in ix:
        o.append('%6.2f  %9.2f (%.2f) %9.2f (%.2f) %9.2f (%.2f) %9.2f' % (
            x[i], j0[i], d['swap_w193_c193_J_se'][i], jw[i], d['swap_wRANS_c193_J_se'][i], jc[i], d['swap_w193_cRANS_J_se'][i], jr[i]))
    o += ['', '#  x/D   lateral at 0.5 D: mean193 (SE)  meanRANS   turb193 (SE)  turbRANS  |  1.0 D: mean193  meanRANS  turb193  turbRANS']
    for i in ix:
        o.append('%6.2f  %8.2f (%.2f) %8.2f  %8.2f (%.2f) %8.2f  |  %8.2f %8.2f %8.2f %8.2f' % (
            x[i], L('193', 'mean', '0.5')[i], d['lat_zx_193_mean_ye0.5_se'][i], L('rans', 'mean', '0.5')[i],
            L('193', 'turb', '0.5')[i], d['lat_zx_193_turb_ye0.5_se'][i], L('rans', 'turb', '0.5')[i],
            L('193', 'mean', '1')[i], L('rans', 'mean', '1')[i], L('193', 'turb', '1')[i], L('rans', 'turb', '1')[i]))
    rat = [('J[w_RANS,c_193]/J[193] (velocity effect)', jw / j0),
           ('J[w_193,c_RANS]/J[193] (coolant-distribution effect)', jc / j0),
           ('J[RANS]/J[193]', jr / j0)]
    for y in ('0.5', '1'):
        rat.append(('RANS/193M inward mean lateral flux, |y| = %s D' % ('1.0' if y == '1' else y), L('rans', 'mean', y) / L('193', 'mean', y)))
        rat.append(('RANS/193M outward turbulent lateral flux, |y| = %s D' % ('1.0' if y == '1' else y), L('rans', 'turb', y) / L('193', 'turb', y)))
    o += ['', '## Derived ratios over x/D 3-20 (min, max, and values at the stations)',
          '# ' + 'x/D'.ljust(60) + '   min    max   ' + ' '.join('%7.0f' % s for s in st)]
    for nm, r in rat:
        o.append('  %-60s %6.2f %6.2f   ' % (nm, np.nanmin(r), np.nanmax(r)) + ' '.join('%7.2f' % r[i] for i in ix))
    o += ['', '# Ratios of the lateral flux at x/D 3 are affected by sign changes of the RANS mean flux (+0.12 at 0.5 D, x/D 3).']
    with open(TAB, 'w') as f:
        f.write('\n'.join(o) + '\n')
    print('\n'.join(o))


def main(yscale=YSCALE_A, out=OUT, png=None):
    d = np.load(NPZ)
    x = d['x_D']
    # Drawn at the printed size for a full column (3.56 in; include width=\columnwidth, scale 1); body text 10 pt -> sf = 10/22.
    set_mpl_params(10 / 22, grid_on=True)
    plt.rcParams['axes.prop_cycle'] = cycler(color=OKABE_ITO)
    FS = plt.rcParams['font.size']
    fig, (A, B) = plt.subplots(2, 1, figsize=(3.56, 5.2), sharex=True)

    # (a) swap test
    combos = [('193', '193', C193, '-', r'$w$: 193M, $c$: 193M'),
              ('RANS', '193', CRANS, '--', r'$w$: RANS, $c$: 193M'),
              ('193', 'RANS', C193, '--', r'$w$: 193M, $c$: RANS'),
              ('RANS', 'RANS', CRANS, '-', r'$w$: RANS, $c$: RANS')]
    for a, b, col, ls, lab in combos:
        j = d['swap_w%s_c%s_J' % (a, b)]
        se = d['swap_w%s_c%s_J_se' % (a, b)]
        A.plot(x, j, ls, color=col, label=lab)
        if np.any(se > 0):
            A.fill_between(x, np.clip(j - se, 1e-3, None), j + se, color=col, alpha=0.25, lw=0)
    A.set_yscale(yscale)
    if yscale == 'linear':
        A.set_ylim(*YLIM_LIN)
    A.set_ylabel(r'$J$ through $z_x$ / g s$^{-1}$m$^{-1}$')
    A.legend(loc='upper right', fontsize=FS - 1, handlelength=1.6)

    # (b) lateral flux through |y| = 0.5 D
    B.axhline(0.0, color='gray', lw=0.7)
    for r, col, nm in (('193', C193, '193M'), ('rans', CRANS, 'RANS')):
        for k, ls, kn in (('mean', '-', 'mean'), ('turb', '--', 'turbulent')):
            v = d['lat_zx_%s_%s_ye0.5' % (r, k)]
            B.plot(x, v, ls, color=col)
            if r == '193':
                se = d['lat_zx_193_%s_ye0.5_se' % k]
                B.fill_between(x, v - se, v + se, color=col, alpha=0.25, lw=0)
    B.set_ylabel(r'lateral flux at $|y|=0.5D$' + '\n' + r'/ g s$^{-1}$m$^{-1}$')
    B.legend(handles=[Line2D([], [], color=C193, ls='-', label='193M'), Line2D([], [], color=CRANS, ls='-', label='RANS'),
                      Line2D([], [], color='gray', ls='-', label='mean-flow'), Line2D([], [], color='gray', ls='--', label='turbulent')],
             loc='lower right', fontsize=FS - 1, ncol=2, handlelength=1.4, columnspacing=0.8)
    B.text(19.8, 0.12, 'outward', ha='right', va='bottom', fontsize=FS - 1, color='0.35')
    B.text(19.8, -0.12, 'inward', ha='right', va='top', fontsize=FS - 1, color='0.35')
    B.set_xlabel(r'$x/D$')
    B.set_xlim(3, 20)
    B.set_xticks([3, 5, 10, 15, 20])
    for a, t in ((A, '(a)'), (B, '(b)')):
        a.text(0.0, 1.02, t, transform=a.transAxes, ha='left', va='bottom')
    fig.tight_layout()
    fig.savefig(out)
    if png:
        fig.savefig(png, dpi=200)
    plt.close(fig)
    return d


if __name__ == '__main__':
    other = 'log' if YSCALE_A == 'linear' else 'linear'
    d = main(YSCALE_A, OUT, REV + '34_swap_lateral_%s.png' % YSCALE_A)
    main(other, SCR + 'swap_lateral_%s.pdf' % other, REV + '34_swap_lateral_%s.png' % other)
    table(d)
