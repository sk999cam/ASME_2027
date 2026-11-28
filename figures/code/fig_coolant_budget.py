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
       '2x2': '/home/sk999/Documents/ASME_2027/planning/review/59_fig_coolant_budget_les.png'}
SIZE = {'row': (7.5, 2.0), '2x2': (7.5, 6.5)}
OKABE_ITO = ['#000000', '#E69F00', '#56B4E9', '#009E73', '#F0E442', '#0072B2', '#D55E00', '#CC79A7']
C193, CRANS = OKABE_ITO[0], OKABE_ITO[5]
CLES = '#E69F00'     # 20M LES, paper-wide colour (owner 2026-11-28)
WB = '/mnt/Elements/DNS/analysis/asme2027_figs/wall_balance/'
SURF_LES, SURF_C50 = WB + 'n7_surface_c193b_b010_les20_zxpt.npz', WB + 'n7_surface_c193b_b010_c50b_zxpt.npz'
SWAP_LES, SWAP_C50 = WB + 'n7_swap_lateral_les20.npz', WB + 'n7_swap_lateral_c50b.npz'


def les_analysis():
    """20M LES by the script's own definitions (fig_mean_features), applied to les20_mean_grid.npz (same n4 grid as the 193M volume and the
    RANS cache) and les20_planes_b.npz (fine planes at x/D 3, 5, 6, 8, 12, 15, 20; the planes at x/D 2, 4, 10 do not exist for the LES).
    Returns (volume analysis dict like out['rans'], fine-plane dict {x/D: {'L'/'R': {crvp, wv ...}}})."""
    from scipy import ndimage as ndi
    MF.THR_S['les'] = MF.THR_S['rans']
    orig = MF.load

    def load(name, sig=None):
        if name != 'les':
            return orig(name, sig)
        sig = MF.SIG if sig is None else sig
        d = np.load(MF.CACHE + '/les20_mean_grid.npz')
        om, (u, v, w), th = d['omx_n'], (d[k] / MF.UINF for k in ('U', 'V', 'W')), d['theta']
        f = {'x': d['x'], 'y': d['y'], 'z': d['z']}
        for k, a in (('om', om), ('u', u), ('v', v), ('w', w), ('th', th)):
            f[k] = ndi.gaussian_filter(np.nan_to_num(np.asarray(a, dtype=np.float64)), sig, mode='nearest')
        return f
    MF.load = load
    try:
        F, ix, R = MF.analyse('les')
    finally:
        MF.load = orig
    o = MF.cross_to(F, ix, R)
    vol = dict(F=F, ix=ix, R=R, o=o, q=MF.quantities(R, o))
    pl = np.load(MF.CACHE + '/les20_planes_b.npz')
    y, z = pl['y'], pl['z']
    res = {}
    for xd in sorted(int(k[4:]) for k in pl.files if k.startswith('om_x')):      # every station the file holds (picks up added planes)
        om, V, W = pl['om_x%d' % xd], pl['V_x%d' % xd] / MF.UINF, pl['W_x%d' % xd] / MF.UINF
        oms = ndi.gaussian_filter(np.nan_to_num(np.asarray(om, dtype=np.float64)), MF.SIG_P, mode='nearest')
        Vs = ndi.gaussian_filter(np.nan_to_num(np.asarray(V, dtype=np.float64)), 3, mode='nearest')
        Ws = ndi.gaussian_filter(np.nan_to_num(np.asarray(W, dtype=np.float64)), 3, mode='nearest')
        far = (np.abs(y)[:, None] > 2.0) & (z[None, :] > 1.5)
        r = {'noise': float(np.sqrt(np.mean(oms[far] ** 2)))}
        k = int(np.argmin(np.abs(MF.XS - xd)))
        for sd, s in (('L', 1.0), ('R', -1.0)):
            yv, zv = R[sd]['yc'][k], R[sd]['zc'][k]
            wv = MF.pick(-s * oms, Vs, Ws, y, z, yv, r['noise'])
            ov = MF.pick(s * oms, Vs, Ws, y, z, yv, r['noise'], same_sign_after=abs(wv['y']) if not np.isnan(wv['y']) else abs(yv) + 1.0)
            r[sd] = dict(wv, ov=ov, crvp=MF.crvp_fine(om, y, z, yv, zv, s))
        res[xd] = r
    return vol, res
XMIN = 4.0           # x/D 4-20 on all panels (owner 2026-11-27: z_x confirmed from x/D 4)
YLIM_NEG = (-1.5, 50)   # panels (b), (c): just below the lowest plotted value (LES mean flux -0.90, lower SE edge about -1.1)
YLIM_C = (0, 50)    # linear; values at x/D 4 are 22-24 (193M) and 45 (RANS mean flux)


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


def circulation_les(planes):
    """Same quantities as circulation() for the LES from its fine planes (stations without a LES plane are absent): CRVP and WV, mean of the sides."""
    st = [s for s in sorted(planes) if s >= XMIN]
    g = np.array([[planes[xd][sd]['crvp']['g'] for xd in st] for sd in 'LR'])
    wx, wm, wop = [], [], []
    for xd in st:
        w = planes[xd]
        gs = [w[sd]['g'] for sd in 'LR' if not np.isnan(w[sd]['g'])]
        if gs:
            wx.append(xd)
            wm.append(np.mean(gs))
            wop.append(any(w[sd]['below3'] or not w[sd]['sw'] for sd in 'LR' if not np.isnan(w[sd]['g'])))
    return dict(crvp=(np.array(st), g.mean(0)), wv=(np.array(wx), np.array(wm), np.array(wop)))


def panel_a(ax, out, les=None):
    ST = np.array([s for s in MF.STATIONS if s >= XMIN])
    R = circulation(out, ST)
    if les is not None:
        L = circulation_les(les)
        x, m = L['crvp']
        ax.plot(x, m, '-', color=CLES, marker='o', ms=2.5)
        x, m, op = L['wv']
        ax.plot(x, m, '-.', color=CLES, lw=0.8)
        for xv, gv, o_ in zip(x, m, op):
            ax.plot([xv], [gv], 's', color=CLES, mfc='none' if o_ else CLES, ms=3)
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
    if os.path.exists(SURF_LES):                          # 20M LES: mean and turbulent (resolved + mean SGS) vertical flux, band = jackknife SE (6 batches)
        L0 = np.load(SURF_LES)
        ml = L0['x_D'] >= XMIN
        xl = L0['x_D'][ml]
        for key, tag, ls in (('Jm_w', 'les20_Jm_w', '-'), ('Jt_w', 'les20_Jt_w_tot', '--')):
            j, se = L0[tag][ml], L0['les20_' + key + '_se'][ml]       # SE of the resolved part (the SGS part has no SE)
            ax.plot(xl, j, ls, color=CLES)
            ax.fill_between(xl, j - se, j + se, color=CLES, alpha=0.25, lw=0)      # negative LES mean flux is drawn, not clipped
    ax.axhline(0.0, color='gray', lw=0.7)
    ax.set_ylim(*YLIM_NEG)   # same limits as panel (c)
    ax.set_ylabel(r'$J_w$ / g s$^{-1}$m$^{-1}$')


def panel_c(ax, d, dl=None):
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
    if dl is not None:                                    # 20M LES counterparts: J[w_L, c_193], J[w_193, c_L] (dashed, swapped), J[w_L, c_L] (solid)
        xl = dl['x_D']
        for a, b, ls, lab in (('les20', 'les20', '-', r'$w$: L, $c$: L'),):      # LES total only; the two swapped LES values are in coolant_budget_les_table.txt
            j, se = dl['swap_w%s_c%s_J' % (a, b)], dl['swap_w%s_c%s_J_se' % (a, b)]
            ax.plot(xl, j, ls, color=CLES, label=lab)
            if np.any(se > 0):
                ax.fill_between(xl, j - se, j + se, color=CLES, alpha=0.25, lw=0)
    ax.axhline(0.0, color='gray', lw=0.7)
    ax.set_ylim(*YLIM_NEG)
    ax.set_ylabel(r'$J_w$ / g s$^{-1}$m$^{-1}$')
    ax.legend(loc='upper right', fontsize=plt.rcParams['font.size'] - 1, handlelength=1.2, labelspacing=0.25, borderpad=0.2,
              handletextpad=0.5, borderaxespad=0.2)


def panel_d(ax, d, dl=None):
    x = d['x_D']
    ax.axhline(0.0, color='gray', lw=0.7)
    if dl is not None:                                    # 20M LES: mean, and turbulent = resolved + mean SGS
        for k, ls in (('mean', '-'), ('turb', '--')):
            v = dl['lat_zx_les20_%s_ye0.5' % k] + (dl['lat_zx_les20_sgs_ye0.5'] if k == 'turb' else 0.0)
            se = dl['lat_zx_les20_%s_ye0.5_se' % k]
            ax.plot(dl['x_D'], v, ls, color=CLES)
            ax.fill_between(dl['x_D'], v - se, v + se, color=CLES, alpha=0.25, lw=0)
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


def draw(layout, out, d, pdf=None, les=None, dl=None):
    # Drawn at the printed size (row: 7.5 in, figure* width=\textwidth; 2x2: 5.6 in, figure* width=5.6in or 0.75\textwidth), scale 1; body text 10 pt.
    set_mpl_params(10 / 22, grid_on=True)
    plt.rcParams['axes.prop_cycle'] = cycler(color=OKABE_ITO)
    size = SIZE[layout]
    if layout == 'row':
        fig, axs = plt.subplots(1, 4, figsize=size)
    else:
        fig, axs = plt.subplots(2, 2, figsize=size)
        axs = axs.ravel()
    panel_a(axs[0], out, les)
    panel_b(axs[1])
    panel_c(axs[2], d, dl)
    panel_d(axs[3], d, dl)
    for i, (a, t) in enumerate(zip(axs, 'abcd')):
        a.set_xlim(XMIN, 20)
        a.set_xticks([5, 10, 15, 20])
        if layout == 'row' or i >= 2:
            a.set_xlabel(r'$x/D$')
        a.text(0.0, 1.02, '(%s)' % t, transform=a.transAxes, ha='left', va='bottom')
    fs = plt.rcParams['font.size'] - 1
    fig.tight_layout(pad=0.3, w_pad=0.5, h_pad=0.6, rect=(0, 0.1 if layout == 'row' else 0.06, 1, 1))
    handles = [Line2D([], [], color=C193, ls='-', label='193M 3DNS'), Line2D([], [], color=CRANS, ls='-', label='RANS 50M')]
    if les is not None:
        handles.append(Line2D([], [], color=CLES, ls='-', label='20M LES'))
    handles += [Line2D([], [], color='gray', ls='-', label='mean-flow (b, d)'), Line2D([], [], color='gray', ls='--', label='turbulent (b, d)')]
    fig.legend(handles=handles, loc='lower center', ncol=len(handles), frameon=False, fontsize=fs, handlelength=1.8, columnspacing=1.5,
               bbox_to_anchor=(0.5, 0.0))
    print(layout, 'figsize', size, 'axes box w x h in:', [(round(a.get_position().width * size[0], 2), round(a.get_position().height * size[1], 2)) for a in axs])
    if pdf:
        fig.savefig(pdf)
    fig.savefig(PNG[layout], dpi=200)
    plt.close(fig)


def write_table(out, lvol, lpl, d, dl):
    """data/coolant_budget_les_table.txt: vertical flux through z_x, swap test, lateral flux, circulation; 193M, RANS 50M, 20M LES, 50M 3DNS."""
    sr, sl, sc = np.load(FE.SURF), np.load(SURF_LES), np.load(SURF_C50)
    dc = np.load(SWAP_C50)
    X = [4, 5, 7, 10, 15, 20]
    ix = lambda a, q: int(np.argmin(np.abs(a['x_D'] - q)))      # noqa: E731
    L = ['# 20M LES in the coolant budget (fig_coolant_budget.py). Fluxes g s^-1 m^-1, per metre of x, one pitch; z = z_x(x) (near_wall_boundary_zx.txt).',
         '# Vertical flux through z_x: mean-flow (Jm_w) and turbulent (Jt_w); LES turbulent = resolved + mean SGS (les20_Jt_w_tot; the SGS part is <= 0.05).',
         '# SE = jackknife over the case batches (193M 10, LES 6; LES SE is that of the resolved part, the SGS part has no SE). 50M 3DNS = c50b (not plotted).',
         '# Each ratio uses the 193M of its own file (identical to the RANS-file 193M at x/D >= 3.5 within 0.1 %; differences at x/D 3: 1.1 in Jm_w).', '',
         'VERTICAL FLUX', f'{"x/D":>5s} {"case":>9s} {"Jm_w":>8s} {"SE":>6s} {"Jt_w":>8s} {"SE":>6s} {"Jm/Jm193":>9s} {"Jt/Jt193":>9s}']
    cases = [('193M', sr, '193_Jm_w', '193_Jm_w_se', '193_Jt_w', '193_Jt_w_se', sr, '193_Jm_w', '193_Jt_w'),
             ('RANS 50M', sr, 'rans_Jm_w', None, 'rans_Jt_w', None, sr, '193_Jm_w', '193_Jt_w'),
             ('20M LES', sl, 'les20_Jm_w', 'les20_Jm_w_se', 'les20_Jt_w_tot', 'les20_Jt_w_se', sl, '193_Jm_w', '193_Jt_w'),
             ('50M 3DNS', sc, 'c50b_Jm_w', 'c50b_Jm_w_se', 'c50b_Jt_w', 'c50b_Jt_w_se', sc, '193_Jm_w', '193_Jt_w')]
    for q in X:
        for nm, f, km, kms, kt, kts, fr, rm, rt in cases:
            i = ix(f, q)
            sm_ = f'{f[kms][i]:6.2f}' if kms else f'{"-":>6s}'
            st_ = f'{f[kts][i]:6.2f}' if kts else f'{"-":>6s}'
            L.append(f'{q:5.1f} {nm:>9s} {f[km][i]:8.2f} {sm_} {f[kt][i]:8.2f} {st_} {f[km][i] / fr[rm][i]:9.2f} {f[kt][i] / fr[rt][i]:9.2f}')
    L += ['', 'x/D (from 3.0, steps of 0.25) where the turbulent flux first exceeds the mean-flow flux (Jt_w > Jm_w):']
    for nm, f, km, _, kt, *_ in cases:
        k = np.where(f[kt] > f[km])[0]
        L.append(f'  {nm:9s} {f["x_D"][k[0]]:.2f}' if len(k) else f'  {nm:9s} not before x/D 20')
    L += ['', 'SWAP TEST: mean-flow flux through z_x with the velocity w from one run and the coolant c from the other, J[w, c]',
          f'{"x/D":>5s} {"J[193,193]":>11s} | {"J[R,193]":>9s} {"J[193,R]":>9s} {"J[R,R]":>8s} | {"J[L,193]":>9s} {"J[193,L]":>9s} {"J[L,L]":>8s} | {"J[50,193]":>9s} {"J[193,50]":>9s} {"J[50,50]":>8s}']
    for q in X:
        i = ix(d, q)
        L.append(f'{q:5.1f} {d["swap_w193_c193_J"][i]:11.2f} | {d["swap_wRANS_c193_J"][i]:9.2f} {d["swap_w193_cRANS_J"][i]:9.2f} {d["swap_wRANS_cRANS_J"][i]:8.2f} | '
                 f'{dl["swap_wles20_c193_J"][i]:9.2f} {dl["swap_w193_cles20_J"][i]:9.2f} {dl["swap_wles20_cles20_J"][i]:8.2f} | '
                 f'{dc["swap_wc50b_c193_J"][i]:9.2f} {dc["swap_w193_cc50b_J"][i]:9.2f} {dc["swap_wc50b_cc50b_J"][i]:8.2f}')
    L += ['', 'LATERAL FLUX through |y| = 0.5 D and 1.0 D below z_x, mean-flow (m) and turbulent (t; LES: resolved + mean SGS), + away from the centreline',
          f'{"x/D":>5s} {"|y|/D":>6s} {"193 m":>7s} {"193 t":>7s} {"RANS m":>7s} {"RANS t":>7s} {"LES m":>7s} {"LES t":>7s} {"50M m":>7s} {"50M t":>7s}']
    for q in X:
        for ye in ('0.5', '1'):
            i = ix(d, q)
            L.append(f'{q:5.1f} {ye:>6s} {d["lat_zx_193_mean_ye" + ye][i]:7.2f} {d["lat_zx_193_turb_ye" + ye][i]:7.2f} {d["lat_zx_rans_mean_ye" + ye][i]:7.2f} '
                     f'{d["lat_zx_rans_turb_ye" + ye][i]:7.2f} {dl["lat_zx_les20_mean_ye" + ye][i]:7.2f} '
                     f'{dl["lat_zx_les20_turb_ye" + ye][i] + dl["lat_zx_les20_sgs_ye" + ye][i]:7.2f} {dc["lat_zx_c50b_mean_ye" + ye][i]:7.2f} {dc["lat_zx_c50b_turb_ye" + ye][i]:7.2f}')
    # circulation
    gq = lambda o, xq: float(o['q']['Gamma CRVP 0.2'][int(np.argmin(np.abs(MF.XS - xq)))])       # noqa: E731
    L += ['', 'CIRCULATION, CRVP and wall vortex (WV), Gamma/(U_inf D), mean of the two sides.',
          'A. Volume analysis (fig_mean_features.analyse; omega_x smoothed SIG, contour 0.2 of the local peak), n4 grid, all three runs, same definition:',
          f'{"x/D":>5s} {"193M":>8s} {"RANS":>8s} {"LES":>8s} {"LES/193M":>9s} {"RANS/193M":>10s}']
    for q in (4, 10, 20):
        a, b, c = gq(out['dns'], q), gq(out['rans'], q), gq(lvol, q)
        L.append(f'{q:5d} {a:8.4f} {b:8.4f} {c:8.4f} {c / a:9.3f} {b / a:10.3f}')
    L += ['B. Fine planes (the quantities of panel a), stations x/D 4-20 of the panel. LES planes read from les20_planes_b.npz: x/D ' +
          ', '.join(str(s) for s in sorted(lpl)) + '; "-" = no plane of that run at the station (RANS WV: none beyond x/D 4).',
          f'{"x/D":>5s} {"CRVP 193M":>10s} {"RANS":>8s} {"LES":>8s} {"LES/193M":>9s} | {"WV 193M":>8s} {"RANS":>8s} {"LES":>8s} {"LES/193M":>9s}']
    W = out['wvp']
    for q in sorted(s for s in MF.STATIONS if s >= XMIN):
        def v(nm, kind):
            if nm == 'les':
                if q not in lpl:
                    return np.nan
                w = lpl[q]
            else:
                w = W[nm][q]
            return float(np.nanmean([w[s]['crvp']['g'] if kind == 'crvp' else w[s]['g'] for s in 'LR']))
        r = [v(n_, k_) for k_ in ('crvp', 'wv') for n_ in ('dns', 'rans', 'les')]
        f_ = lambda z: f'{z:8.4f}' if np.isfinite(z) else f'{"-":>8s}'      # noqa: E731
        rt = lambda a, b: f'{a / b:9.3f}' if np.isfinite(a) and np.isfinite(b) and b else f'{"-":>9s}'      # noqa: E731
        L.append(f'{q:5d} {r[0]:10.4f} {f_(r[1])} {f_(r[2])} {rt(r[2], r[0])} | {f_(r[3])} {f_(r[4])} {f_(r[5])} {rt(r[5], r[3])}')
    L += ['', 'CHECK against the coordinator values (x/D 5 / 10 / 20, g s^-1 m^-1), this script:']
    g_ = lambda f, k: ' / '.join(f'{f[k][ix(f, q)]:.2f}' for q in (5, 10, 20))      # noqa: E731
    lat_t = ' / '.join(f'{dl["lat_zx_les20_turb_ye0.5"][ix(dl, q)] + dl["lat_zx_les20_sgs_ye0.5"][ix(dl, q)]:.2f}' for q in (5, 10, 20))
    L += [f'  LES mean flux      {g_(sl, "les20_Jm_w")}   (quoted 12.1 / 0.9 / -0.9); 193M {g_(sr, "193_Jm_w")} (13.7 / 2.4 / 0.85)',
          f'  LES turbulent flux {g_(sl, "les20_Jt_w_tot")}   (quoted 17.5 / 7.5 / 5.8); 193M {g_(sr, "193_Jt_w")} (15.8 / 6.4 / 4.9)',
          f'  swap J[L,193]      {g_(dl, "swap_wles20_c193_J")}   (quoted 13.1 / 1.3 / -0.7);  J[193,L] {g_(dl, "swap_w193_cles20_J")} (12.7 / 1.9 / 0.6)',
          f'  lateral 0.5 D mean {g_(dl, "lat_zx_les20_mean_ye0.5")}   (quoted -6.3 / -4.1 / -1.3);  turbulent {lat_t} (quoted 2.6 / 1.8 / 1.1)']
    TAB = os.path.join(HERE, 'data', 'coolant_budget_les_table.txt')
    with open(TAB, 'w') as fh:
        fh.write('\n'.join(L) + '\n')
    print('\n'.join(L))


def main():
    warnings.filterwarnings('ignore')
    out = MF.run()
    d = np.load(SL.NPZ)
    dl = np.load(SWAP_LES)
    lvol, lpl = les_analysis()
    draw(LAYOUT, out, d, OUT, les=lpl, dl=dl)
    write_table(out, lvol, lpl, d, dl)
    print('wrote', OUT)


if __name__ == '__main__':
    main()
