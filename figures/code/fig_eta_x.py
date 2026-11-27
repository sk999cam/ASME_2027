"""Adiabatic effectiveness against x/D for the six rungs and the 7-7-7 data (paper Fig. F12, Results 4.1).

Usage: python3 fig_eta_x.py          (no arguments; writes /home/sk999/Documents/ASME_2027/figures/fig_eta_x.pdf)

Panels
  (a) laterally averaged eta-bar over one pitch.
  (b) centreline eta, |y/D| < 0.1 (both sides folded), bins of 0.1 D in x.
  (c) mid-pitch eta, |y/D| = 2.91 +- 0.1 (both sides folded; as analysis/far_film/far_film.py), bins of 0.1 D in x.
eta = (T0,fs - T_w) / (T0,fs - T0,c), T0,fs 300 K, T0,c 199.02 K (adiabatic wall temperature of the means).

Data
  RANS 5M, RANS 50M: starccm/rans_{5M_Lm25,50M_Lm25b}_BR15_FST_ZPG_p0dns/post_run1/wall_plate.csv (plate faces; T).
  20M LES (C20bW), 50M (C50b), 100M (C100b) 3DNS: stat3 batch b006; 193M (C193b): stat3 batch b010
    (pictures_job18/stat3_b010). Mean T at the plate vertices (k = ktip of the tip blocks). Panel (a): eta_lateral.txt of
    dns/python/plot_stat3.py (pictures_job3 / job9 / job9 / job18); panels (b), (c): the plate vertices, cached by plate_points().
  gamma_mid(x) of the 193M (panel (c), secondary axis, dashed grey): intermittency of the plane sampler z010 (z/D 0.1) of the
    b010 window, threshold 0.0385 U^2/D, averaged over the mid-pitch band |y/D| 2.81 to 3.0 (the band of the eta line,
    |y/D| = 2.91 +- 0.1, as far as the sampler reaches), helper _gamma_extract.py, 0 to 20 D.
  Experiment: ETA_CAzeroStudies, eta_SmoothHoles_lowFST_DR1p5_M1p5__id727 (7-7-7, DR 1.5, BR 1.5, Tu 0.5 %).
    (a) mean over the map rows (as validation_50M/eta.py); (b), (c) the three pitches folded onto one
    (holes at z/D 0, +-6), as far_film.py. Markers are every 8th point of the 0.1 D map.
Cache: /mnt/Elements/DNS/analysis/asme2027_figs/fig_eta_x/plate_<rung>.npz (x/D, |y|/D, eta of the plate vertices).

Status: complete. 193M on batch b010, the other rungs on b006. Centreline and mid-pitch of 3DNS rungs are vertex values binned in x;
the RANS values are face centres binned in the same way.
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from cycler import cycler

sys.path.insert(0, '/home/sk999/Documents/3DNS_FP')
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from functions import set_mpl_params
import _gamma_extract as GX  # noqa: E402
REPO = '/home/sk999/Documents/4dns_filmcooled'
sys.path.insert(0, REPO + '/starccm/python')
sys.path.insert(0, REPO + '/dns/python')
from bl_post import D, T0_FS, T0_C, read_csv, lateral_mean  # noqa: E402

OKABE_ITO = ['#000000', '#E69F00', '#56B4E9', '#009E73', '#F0E442', '#0072B2', '#D55E00', '#CC79A7']
OUT = '/home/sk999/Documents/ASME_2027/figures/fig_eta_x.pdf'
CACHE = '/mnt/Elements/DNS/analysis/asme2027_figs/fig_eta_x'
BATCH = 6
BATCHES = {'193M': 10}
EXP = '/mnt/Tank/777_Film_Cooling_Database/ETA_CAzeroStudies/'
EXP_FILE = 'eta_SmoothHoles_lowFST_DR1p5_M1p5__id727.txt'

# label, colour, style
STYLE = {'193M': ('#000000', '-'), '100M': ('#0072B2', '--'), '50M': ('#D55E00', ':'),
         '20M LES': ('#E69F00', '-.'), 'RANS 50M': ('#009E73', '-'), 'RANS 5M': ('#CC79A7', '--')}
RANS = {'RANS 5M': REPO + '/starccm/rans_5M_Lm25_BR15_FST_ZPG_p0dns',
        'RANS 50M': REPO + '/starccm/rans_50M_Lm25b_BR15_FST_ZPG_p0dns'}
DNS = {'20M LES': ('/mnt/Tank/DNS/film_cooling/C20bW_BR15', 'pictures_job3'),
       '50M': ('/mnt/Tank/DNS/film_cooling/C50b_BR15', 'pictures_job9'),
       '100M': ('/mnt/Elements/DNS/film_cooling/C100b_BR15', 'pictures_job9'),
       '193M': ('/mnt/Elements/DNS/film_cooling/C193b_BR15', 'pictures_job18')}
ORDER = ['RANS 5M', 'RANS 50M', '20M LES', '50M', '100M', '193M']
XMAX = 20.0           # owner rule 2026-11-26: nothing beyond x/D 20 is shown or reported
XE = np.arange(-2.0, XMAX + 0.05, 0.1)
WINDOWS = ((2.5, 7), (7, 15), (15, 20))


def eta_of(T):
    return (T0_FS - T) / (T0_FS - T0_C)


def plate_points(rung):
    """x/D, |y|/D and mean eta of the plate vertices (tip blocks, k = ktip) of the stat3 batch b006 (cached)."""
    os.makedirs(CACHE, exist_ok=True)
    BATCH = BATCHES.get(rung, 6)
    f = os.path.join(CACHE, 'plate_' + rung.replace(' ', '_') + (f'_b{BATCH:03d}' if rung in BATCHES else '') + '.npz')
    if os.path.exists(f):
        Z = np.load(f)
        return Z['x'], Z['y'], Z['e']
    import json
    from tke_io import IDX
    run = DNS[rung][0]
    info = json.load(open(os.path.join(run, 'case_info.json')))
    kp = info['ktip'] - 1
    mt = [float(l.split()[1]) for l in open(os.path.join(run, 'stat3_batches.txt'))
          if l.split() and l.split()[0] == str(BATCH)][0]
    xs, ys, es = [], [], []
    for b in info['blocks']:
        if b['hole']:
            continue
        npt = b['nk'] * b['nj'] * b['ni']
        sf = os.path.join(run, f'stat3_{b["n"]}_b{BATCH:03d}')
        ntk = os.path.getsize(sf) // 8 // npt
        S = np.memmap(sf, '<f8', 'r', shape=(b['nk'], b['nj'], b['ni'], ntk))
        T = np.array(S[kp, :, :, IDX['T']]) / mt
        g = np.memmap(os.path.join(run, f'grid_{b["n"]}'), '<f8', 'r', shape=(b['nk'], b['nj'], b['ni'], 3))
        xs.append(np.array(g[kp, :, :, 0]).ravel() / D)
        ys.append(np.array(g[kp, :, :, 1]).ravel() / D)
        es.append(eta_of(T).ravel())
        print(rung, 'block', b['n'], T.size, flush=True)
    x, y, e = np.concatenate(xs), np.concatenate(ys), np.concatenate(es)
    np.savez(f, x=x, y=y, e=e, mean_time=mt)
    return x, y, e


def rans_points(rung):
    w = read_csv(RANS[rung] + '/post_run1/wall_plate.csv')
    return w['X'] / D, w['Y'] / D, eta_of(w['Temperature'])


def exp_points():
    x = np.loadtxt(EXP + 'x_D.txt').ravel()
    z = np.loadtxt(EXP + 'z_D.txt').ravel()
    e = np.loadtxt(EXP + EXP_FILE).ravel()
    ok = np.isfinite(e)
    return x[ok], np.abs(np.mod(z + 3.0, 6.0) - 3.0)[ok], e[ok]


def binned(c, v, edges, sel):
    i = np.digitize(c[sel], edges) - 1
    n = len(edges) - 1
    ok = (i >= 0) & (i < n)
    s = np.bincount(i[ok], v[sel][ok], n)
    k = np.bincount(i[ok], None, n)
    return 0.5 * (edges[1:] + edges[:-1]), np.where(k > 0, s / np.maximum(k, 1), np.nan)


def lines(x, y, e):
    xc, cl = binned(x, e, XE, np.abs(y) < 0.1)
    _, mp = binned(x, e, XE, np.abs(np.abs(y) - 2.91) < 0.1)
    return xc, cl, mp


def lateral(rung):
    """x/D and eta-bar."""
    if rung in RANS:
        w = read_csv(RANS[rung] + '/post_run1/wall_plate.csv')
        xc, lat = lateral_mean(w['X'], w['Y'], eta_of(w['Temperature']), np.arange(-2 * D, 25.5 * D, 0.1 * D))
        return xc / D, lat
    run, pic = DNS[rung]
    b = BATCHES.get(rung, 6)
    d = np.loadtxt(os.path.join(run, pic, f'stat3_b{b:03d}', 'eta_lateral.txt'))
    return d[:, 0], d[:, 1]


def window_means(x, e):
    return [np.nanmean(e[(x >= a) & (x <= b)]) for a, b in WINDOWS]


def at(x, e, q):
    """Mean of e about x/D q: +-0.26 D, except at the domain end (q = XMAX), where the window is 19.5-20.0 as in the STAR tables."""
    m = (x >= XMAX - 0.5) & (x <= XMAX) if q >= XMAX else np.abs(x - q) < 0.26
    return np.nanmean(e[m])


def main():
    set_mpl_params(0.5, grid_on=True)
    plt.rcParams['axes.prop_cycle'] = cycler(color=OKABE_ITO)
    fig, axs = plt.subplots(3, 1, figsize=(5.2, 7.2), sharex=True)
    xe = np.loadtxt(EXP + 'x_D.txt')[0]
    ee = np.nanmean(np.loadtxt(EXP + EXP_FILE), axis=0)
    ex, ey, ev = exp_points()
    xl, cl, mp = lines(ex, ey, ev)
    exp_lines = (xl, cl, mp)
    mid = {}
    out = ['window means of eta-bar (x/D 2.5-7, 7-15, 15-20):']
    for r in ORDER:
        c, ls = STYLE[r]
        x, y, e = rans_points(r) if r in RANS else plate_points(r)
        xc, cl_, mp_ = lines(x, y, e)
        xb, lat = lateral(r)
        s = np.isfinite(lat) & (xb > 2.36) & (xb <= XMAX)
        axs[0].plot(xb[s], lat[s], color=c, ls=ls, label=r)
        m = np.isfinite(cl_) & (xc > 2.36) & (xc <= XMAX)
        axs[1].plot(xc[m], cl_[m], color=c, ls=ls)
        m = np.isfinite(mp_) & (xc > 2.36) & (xc <= XMAX)
        axs[2].plot(xc[m], mp_[m], color=c, ls=ls)
        mid[r] = (xc, mp_)
        out.append(f'{r:10s} ' + ' '.join(f'{v:.3f}' for v in window_means(xb, lat))
                   + '   centreline x/D 10/20: ' + ' '.join(f'{at(xc, cl_, q):.3f}' for q in (10, 20))
                   + '   mid-pitch x/D 10/20: ' + ' '.join(f'{at(xc, mp_, q):.3f}' for q in (10, 20)))
    out.append('exp id727  ' + ' '.join(f'{v:.3f}' for v in window_means(xe, ee))
               + '   centreline x/D 10/20: ' + ' '.join(f'{at(exp_lines[0], exp_lines[1], q):.3f}' for q in (10, 20))
               + '   mid-pitch x/D 10/20: ' + ' '.join(f'{at(exp_lines[0], exp_lines[2], q):.3f}' for q in (10, 20)))
    print('\n'.join(out))
    kw = dict(marker='o', ls='none', color='k', mfc='none', ms=3.5, label='experiment')
    s = (xe >= 2.2) & (xe <= XMAX) & np.isfinite(ee)
    axs[0].plot(xe[s][::8], ee[s][::8], **kw)
    for ax, i in ((axs[1], 1), (axs[2], 2)):
        xx, vv = exp_lines[0], exp_lines[i]
        m = np.isfinite(vv) & (xx >= 2.2) & (xx <= XMAX)
        ax.plot(xx[m][::4], vv[m][::4], **kw)
    # gamma_mid(x) of the 193M on a secondary axis of panel (c)
    G = GX.load_z010('c193b')['gamma'][0]                      # (x, y) at 0.0385 U^2/D
    jm = np.where(np.abs(GX.YG) >= 2.81 - 1e-9)[0]
    gmid = G[:, jm].mean(1)
    xg = GX.XG
    ok = xg <= XMAX + 1e-9
    ax2 = axs[2].twinx()
    ax2.plot(xg[ok], gmid[ok], color='0.5', ls='--', lw=axs[2].lines[0].get_linewidth(), label=r'$\gamma_{mid}$ (193M)')
    ax2.set_ylim(0, 1)
    ax2.set_ylabel(r'$\gamma_{mid}$ (193M)', color='0.35')
    ax2.tick_params(axis='y', colors='0.35')
    ax2.grid(False)
    # numbers for the report
    xr, er = mid['193M']
    print('193M mid-pitch eta and gamma_mid (window +-0.26 D):')
    for q in (5, 10, 15, 20):
        mm = np.abs(xr - q) < 0.26
        print(f'  x/D {q}: eta_mid {np.nanmean(er[mm]):.4f}  gamma_mid {gmid[np.abs(xg - q) < 0.26].mean():.3f}  exp {np.nanmean(exp_lines[2][np.abs(exp_lines[0] - q) < 0.26]):.4f}')
    xi = np.interp(xr, xg[ok], gmid[ok])
    sel = (xi > 0.8) & (xr <= XMAX) & np.isfinite(er)
    ex_i = np.interp(xr, exp_lines[0][np.isfinite(exp_lines[2])], exp_lines[2][np.isfinite(exp_lines[2])])
    if sel.any():
        diff = er[sel] - ex_i[sel]
        print(f'gamma_mid > 0.8 from x/D {xr[sel].min():.1f} to {xr[sel].max():.1f}: 3DNS minus experiment mid-pitch eta mean {diff.mean():.4f} (min {diff.min():.4f}, max {diff.max():.4f}); 3DNS {er[sel].mean():.4f}, exp {ex_i[sel].mean():.4f}')
    else:
        print('gamma_mid never above 0.8 at x/D <= 20')
    ls_ = np.where(gmid[ok] > 0.8)[0]
    print('gamma_mid first > 0.8 at x/D', xg[ok][ls_[0]] if len(ls_) else None)
    axs[0].set_ylim(0, 0.45)
    axs[1].set_ylim(0.2, 0.9)
    axs[2].set_ylim(0, 0.06)
    axs[2].set_xlim(0, XMAX)
    for ax, lab, yl in zip(axs, '(a) (b) (c)'.split(), (r'$\bar{\eta}$', r'$\eta$ (centreline)', r'$\eta$ (mid-pitch)')):
        ax.set_ylabel(yl)
        ax.text(0.02, 0.92, lab, transform=ax.transAxes, va='top', ha='left')
    axs[2].set_xlabel(r'$x/D$')
    axs[1].legend(*axs[0].get_legend_handles_labels(), loc='lower left', frameon=False, ncol=2)
    fig.tight_layout()
    fig.savefig(OUT)


if __name__ == '__main__':
    main()
