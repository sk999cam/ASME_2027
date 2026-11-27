"""F8c fig_flux_eta: effect of the mean-flow differences on the coolant flux and on the effectiveness, RANS 50M against 193M 3DNS (b010), x/D 3-20.

Usage: python3 fig_flux_eta.py      (writes /home/sk999/Documents/ASME_2027/figures/fig_flux_eta.pdf)

(a) Coolant flux through the inclined near-wall surface z = z_x(x) (final point-file z_x table, x/D 3-20), one pitch, g s^-1 m^-1, + away from the wall:
    VERTICAL parts only, mean (upwash) Jm_w and turbulent Jt_w; the slope term (streamwise flux times dz_x/dx) is left out. 193M solid with a
    +-2 SE band (jackknife of the _w parts, 10 batches), RANS dashed; negative stretches would be drawn as |J| dotted (none here, so the axis is log).
    Source: wall_balance/n7_surface_c193b_b010_rans50_zxpt.npz (constant SURF, keys KEY_M, KEY_T); its stored surface equals
    data/near_wall_boundary_zx.txt (maximum difference 0.000 D on x/D 3-20).
    The fixed-height version (z = 0.5 D, N7 series) stays available: `python3 fig_flux_eta.py variants` (scratchpad outputs).
(b) Ratios RANS/193M against x/D, linear axis, reference line at 1:
    - lateral half-width of the coolant plume: y where the z-integrated mean theta falls to half its centreline value, mean of the two sides, from the
      smoothed 193M and RANS volumes (fig_mean_features.py, SIG smoothing, same grid for both); ratio = RANS / 193M at each x.
    - laterally averaged effectiveness eta_bar over one pitch and centreline effectiveness eta_c (|y/D| < 0.1, bins of 0.1 D), exactly as panels (a) and (b)
      of fig_eta_x.py (same functions lateral(), lines(), plate_points(), rans_points(); RANS 50M plate faces of wall_plate.csv, 193M stat3 b010 plate
      vertices); each curve is interpolated on a common x grid (step 0.25 D, x/D 3 to 20) and the ratio formed as RANS / 193M.
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
import fig_eta_x as EX  # noqa: E402

OUT = '/home/sk999/Documents/ASME_2027/figures/fig_flux_eta.pdf'
OKABE_ITO = ['#000000', '#E69F00', '#56B4E9', '#009E73', '#F0E442', '#0072B2', '#D55E00', '#CC79A7']
XG = np.arange(3.0, 20.0001, 0.25)


def on_grid(x, v):
    ok = np.isfinite(v) & np.isfinite(x)
    return np.interp(XG, x[ok], v[ok])


def ratios():
    """RANS/193M ratios on XG: half-width, eta_bar, eta_c."""
    hw = {}
    eb = {}
    ec = {}
    for nm, rung in (('rans', 'RANS 50M'), ('dns', '193M')):
        F, ix, R = MF.analyse(nm)
        hw[nm] = on_grid(MF.XS, MF.cross_to(F, ix, R)['halfw'])
        x, y, e = EX.rans_points(rung) if nm == 'rans' else EX.plate_points(rung)
        xc, cl, _ = EX.lines(x, y, e)
        xb, lat = EX.lateral(rung)
        eb[nm] = on_grid(xb, lat)
        ec[nm] = on_grid(xc, cl)
    return hw['rans'] / hw['dns'], eb['rans'] / eb['dns'], ec['rans'] / ec['dns']


def note(a, txt, xy, xyt, col, ha='center'):
    a.annotate(txt, xy=xy, xytext=xyt, color=col, ha=ha, va='center', arrowprops=dict(arrowstyle='->', color=col, lw=0.7))


def flux_panel(A, S, h, ylabel=True):
    """Mean and turbulent flux through z = h (log axis, 193M +-2 SE band, RANS dashed, negative stretches as |J| dotted)."""
    C = OKABE_ITO
    xs = S['x_D']
    m3 = xs >= 3.0
    p = 'h%s_' % h
    for key, col in (('Jm', C[1]), ('Jt', C[3])):
        j, se, jr = S[p + '193_' + key], S[p + '193_' + key + '_se'], S[p + 'rans_' + key]
        for v, ls in ((j, '-'), (jr, '--')):
            A.plot(xs[m3], np.where(v > 0, v, np.nan)[m3], ls, color=col)
            A.plot(xs[m3], np.where(v < 0, -v, np.nan)[m3], ':', color=col, lw=0.8)
        A.fill_between(xs[m3], np.clip(j - 2 * se, 1e-3, None)[m3], np.clip(j + 2 * se, 1e-3, None)[m3], color=col, alpha=0.25, lw=0)
    A.set_yscale('log')
    A.set_ylim(0.05, 200)
    if ylabel:
        A.set_ylabel(r'flux / g s$^{-1}$m$^{-1}$')
    A.text(0.97, 0.95, r'$z=%sD$' % h, transform=A.transAxes, ha='right', va='top')
    i8 = int(np.argmin(np.abs(xs - 8.0)))
    note(A, 'mean flux', (8.0, S[p + '193_Jm'][i8]), (6.0, 0.35), C[1])
    i12 = int(np.argmin(np.abs(xs - 12.0)))
    note(A, 'turbulent flux', (12.0, S[p + '193_Jt'][i12]), (15.0, 20.0 if h == 0.5 else 40.0), C[3])
    ng = np.where(S[p + 'rans_Jt'][m3] < 0)[0]
    if len(ng):
        xn = xs[m3][ng[len(ng) // 2]]
        note(A, 'RANS $J_t<0$,\nshown as |J|', (xn, -S[p + 'rans_Jt'][m3][ng[len(ng) // 2]]), (xn + 1.0, 0.12), C[3])


SURF = '/mnt/Elements/DNS/analysis/asme2027_figs/wall_balance/n7_surface_c193b_b010_rans50_zxpt.npz'   # surface z = z_x(x) of the final point-file table (x/D 3-20)
KEY_M, KEY_T = 'Jm_w', 'Jt_w'


def surface_panel(A):
    """Vertical parts of the mean (upwash) and turbulent flux through the inclined surface z = z_x(x), log axis (all values stay positive).
    The slope term (streamwise flux times dz_x/dx) is not included. SE: jackknife of the vertical parts (_w_se)."""
    C = OKABE_ITO
    Z0 = np.load(SURF)
    m35 = Z0['x_D'] >= 3.0
    Z = {k: Z0[k][m35] for k in Z0.files if Z0[k].ndim == 1}                 # x/D 3-20
    x = Z['x_D']
    for key, col in ((KEY_M, C[1]), (KEY_T, C[3])):
        j, se, jr = Z['193_' + key], Z['193_' + key + '_se'], Z['rans_' + key]
        A.plot(x, np.where(j > 0, j, np.nan), '-', color=col)
        A.plot(x, np.where(jr > 0, jr, np.nan), '--', color=col)
        A.plot(x, np.where(j < 0, -j, np.nan), ':', color=col, lw=0.8)
        A.plot(x, np.where(jr < 0, -jr, np.nan), ':', color=col, lw=0.8)
        A.fill_between(x, np.clip(j - 2 * se, 1e-3, None), np.clip(j + 2 * se, 1e-3, None), color=col, alpha=0.25, lw=0)
    A.set_yscale('log')
    A.set_ylim(0.12, 150)
    A.set_ylabel(r'vertical flux / g s$^{-1}$m$^{-1}$')
    A.text(0.97, 0.95, 'vertical part,\n$z=z_x(x)$', transform=A.transAxes, ha='right', va='top')
    i8 = int(np.argmin(np.abs(x - 8.0)))
    note(A, 'mean flux', (13.0, Z['193_' + KEY_M][int(np.argmin(np.abs(x - 13.0)))]), (15.5, 0.27), C[1])
    i14 = int(np.argmin(np.abs(x - 14.0)))
    note(A, 'turbulent flux', (14.0, Z['193_' + KEY_T][i14]), (16.0, 14.0), C[3])


def main(heights=None, out=OUT):
    warnings.filterwarnings('ignore')
    # Drawn at the printed size (0.7 text width = 5.25 in, include scale 1); body text 10 pt -> sf = 10/22.
    set_mpl_params(10 / 22, grid_on=True)
    plt.rcParams['axes.prop_cycle'] = cycler(color=OKABE_ITO)
    C = OKABE_ITO
    n = 1 if heights is None else len(heights)
    fig, axs = plt.subplots(1, n + 1, figsize=(5.25, 2.9), sharey=False)
    B = axs[-1]
    S = np.load(MF.SERIES)
    if heights is None:
        surface_panel(axs[0])
    else:
        for i, h in enumerate(heights):
            flux_panel(axs[i], S, h, ylabel=(i == 0))
            if i > 0:
                axs[i].sharey(axs[0])
                axs[i].tick_params(labelleft=False)
    axs[0].legend(handles=[Line2D([], [], color='k', ls='-', label='193M 3DNS'), Line2D([], [], color='k', ls='--', label='RANS 50M')], loc='lower left',
                  fontsize=plt.rcParams['font.size'], handlelength=1.6)
    A = axs[0]

    rh, rb, rc = ratios()
    B.axhline(1.0, color='gray', lw=0.7)
    B.plot(XG, rh, '-', color=C[6])
    B.plot(XG, rb, '-', color=C[5])
    B.plot(XG, rc, '-', color=C[3])
    B.set_ylabel('RANS / 193M')
    k = lambda q: int(np.argmin(np.abs(XG - q)))
    B.set_ylim(0.55, 1.4)
    note(B, 'plume half-width', (15.5, rh[k(15.5)]), (15.5, 0.9), C[6])
    note(B, r'$\bar{\eta}$ (lateral mean)', (14.0, rb[k(14)]), (14.0, 0.6), C[5])
    note(B, r'$\eta_c$ (centreline)', (16.0, rc[k(16)]), (16.0, 1.35), C[3])
    for a in axs:
        a.set_xlim(3, 20)
        a.set_xticks([5, 10, 15, 20])
        a.set_xlabel(r'$x/D$')
    for a, t in zip(axs, 'abc'):
        a.text(0.0, 1.02, '(%s)' % t, transform=a.transAxes, ha='left', va='bottom')
    fig.tight_layout()
    fig.savefig(out)
    fig.savefig('/tmp/claude-1001/-home-sk999-Documents-ASME-2027/f4c92ceb-6df8-4254-96c3-6552d15d0ecf/scratchpad/flux_eta_view.png', dpi=200)
    for q in (5, 10, 20):
        print('x/D %d: half-width %.3f  eta_bar %.3f  eta_c %.3f' % (q, rh[k(q)], rb[k(q)], rc[k(q)]))


if __name__ == '__main__':
    # python3 fig_flux_eta.py            -> figures/fig_flux_eta.pdf, vertical flux through z = z_x(x) (default)
    # python3 fig_flux_eta.py variants   -> three review variants in the scratchpad (h 0.5; h 1.0; h 0.5 and 1.0 side by side)
    if len(sys.argv) > 1 and sys.argv[1] == 'variants':
        SCR = '/tmp/claude-1001/-home-sk999-Documents-ASME-2027/f4c92ceb-6df8-4254-96c3-6552d15d0ecf/scratchpad/'
        main((0.5,), SCR + 'flux_eta_h050.pdf')
        main((1.0,), SCR + 'flux_eta_h100.pdf')
        main((0.5, 1.0), SCR + 'flux_eta_h050_h100.pdf')
    else:
        main()
