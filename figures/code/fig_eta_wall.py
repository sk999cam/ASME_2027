"""Adiabatic effectiveness eta on the plate, one pitch, x/D 0-20 (paper Fig. F12, Results 4.1).

Usage: python3 fig_eta_wall.py       (no arguments; writes /home/sk999/Documents/ASME_2027/figures/fig_eta_wall.png)

Three stacked panels, shared colour bar (eta 0-0.8, colour map jet): (a) RANS 50M, (b) 20M LES, (c) 193M 3DNS.
eta = (T0,fs - T_w) / (T0,fs - T0,c), T0,fs 300 K, T0,c 199.02 K, from the mean temperature at the plate wall.

Data
  RANS 50M: starccm/rans_50M_Lm25b_BR15_FST_ZPG_p0dns/post_run1/wall_plate.csv (plate faces, T).
  20M LES (C20bW): stat3 batch b006, mean T at the plate vertices (k = ktip of the tip blocks), as
    dns/python/plot_stat3.py (M01_plate_eta); read and cached by plate_points() of fig_eta_x.py in
    /mnt/Elements/DNS/analysis/asme2027_figs/fig_eta_x/.
  193M 3DNS (C193b, reference): stat3 batch b010 (10 FT; C193b_BR15/stat3_<n>_b010, post-processing in
    pictures_job18/stat3_b010), same extraction, read by plate_points_b010() below and cached in
    /mnt/Elements/DNS/analysis/asme2027_figs/fig_eta_wall/plate_193M_b010.npz.
Rendering: linear interpolation of the scattered wall points (scipy LinearNDInterpolator, Delaunay) on a 0.05 D grid. Pixels farther than
  0.25 D from the nearest wall point (the hole exit, where the plate export has no faces) are left white.

Hole exit: the white polygon with a black border is the breakout outline on z = 0 (hole_exit_outline.py, from the RANS 50M
  hole-wall rim, data/hole_exit_outline.txt), same in all panels; it hides the jagged edge of the data mask.
Status: complete.
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.interpolate import LinearNDInterpolator
from scipy.spatial import cKDTree
from matplotlib.patches import Polygon

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, '/home/sk999/Documents/3DNS_FP')
import fig_eta_x as F  # noqa: E402
from functions import set_mpl_params  # noqa: E402
from hole_exit_outline import load as load_outline  # noqa: E402

OUT = '/home/sk999/Documents/ASME_2027/figures/fig_eta_wall.png'
CACHE = '/mnt/Elements/DNS/analysis/asme2027_figs/fig_eta_wall'
VMAX = 0.8
XG = np.arange(-2.0, 20.0001, 0.05)      # owner rule 2026-11-26: nothing beyond x/D 20 is shown
YG = np.arange(-3.0, 3.0001, 0.05)


def plate_points_b010():
    """x/D, y/D and mean eta of the plate vertices (tip blocks, k = ktip) of C193b stat3 batch b010 (cached)."""
    import json
    f = os.path.join(CACHE, 'plate_193M_b010.npz')
    if os.path.exists(f):
        Z = np.load(f)
        return Z['x'], Z['y'], Z['e']
    os.makedirs(CACHE, exist_ok=True)
    from tke_io import IDX
    run = F.DNS['193M'][0]
    batch = 10
    info = json.load(open(os.path.join(run, 'case_info.json')))
    kp = info['ktip'] - 1
    mt = [float(l.split()[1]) for l in open(os.path.join(run, 'stat3_batches.txt'))
          if l.split() and l.split()[0] == str(batch)][0]
    xs, ys, es = [], [], []
    for b in info['blocks']:
        if b['hole']:
            continue
        npt = b['nk'] * b['nj'] * b['ni']
        sf = os.path.join(run, f'stat3_{b["n"]}_b{batch:03d}')
        ntk = os.path.getsize(sf) // 8 // npt
        S = np.memmap(sf, '<f8', 'r', shape=(b['nk'], b['nj'], b['ni'], ntk))
        T = np.array(S[kp, :, :, IDX['T']]) / mt
        g = np.memmap(os.path.join(run, f'grid_{b["n"]}'), '<f8', 'r', shape=(b['nk'], b['nj'], b['ni'], 3))
        xs.append(np.array(g[kp, :, :, 0]).ravel() / F.D)
        ys.append(np.array(g[kp, :, :, 1]).ravel() / F.D)
        es.append(F.eta_of(T).ravel())
        print('193M b010 block', b['n'], T.size, flush=True)
    x, y, e = np.concatenate(xs), np.concatenate(ys), np.concatenate(es)
    np.savez(f, x=x, y=y, e=e, mean_time=mt)
    return x, y, e


def grid_eta(x, y, e):
    """Fold the periodic edges (copy points within 0.3 D of y = +-3 to the other side), interpolate, mask the gaps."""
    m = (x > -2.5) & (x < 20.5)
    x, y, e = x[m], y[m], e[m]
    lo, hi = y < -2.7, y > 2.7
    x = np.r_[x, x[lo], x[hi]]
    y = np.r_[y, y[lo] + 6.0, y[hi] - 6.0]
    e = np.r_[e, e[lo], e[hi]]
    if x.size > 400000:                                    # 193M: thin to ~ the pixel size
        keep = np.random.default_rng(0).random(x.size) < 400000 / x.size
        x, y, e = x[keep], y[keep], e[keep]
    X, Y = np.meshgrid(XG, YG)
    Z = LinearNDInterpolator(np.c_[x, y], e)(X, Y)
    d, _ = cKDTree(np.c_[x, y]).query(np.c_[X.ravel(), Y.ravel()])
    Z = np.ma.masked_where(d.reshape(X.shape) > 0.25, Z)
    return Z


def main():
    set_mpl_params(0.5, grid_on=False)
    HX, HY = load_outline()
    srcs = [('RANS 50M', F.rans_points('RANS 50M')), ('20M LES', F.plate_points('20M LES')),
            ('193M 3DNS', plate_points_b010())]
    fig, axs = plt.subplots(3, 1, figsize=(6.0, 5.6), sharex=True)
    for ax, lab, (n, (x, y, e)) in zip(axs, '(a) (b) (c)'.split(), srcs):
        Z = grid_eta(x, y, e)
        im = ax.pcolormesh(XG, YG, Z, cmap='jet', vmin=0, vmax=VMAX, shading='auto', rasterized=True)
        ax.add_patch(Polygon(np.c_[HX, HY], closed=True, fc='w', ec='k', lw=0.5, zorder=3))
        ax.set_xlim(-2, 20)
        ax.set_ylim(-3, 3)
        ax.set_aspect('equal')
        ax.set_ylabel(r'$y/D$')
        ax.set_yticks([-3, 0, 3])
        ax.text(0.99, 0.93, f'{lab} {n}', transform=ax.transAxes, va='top', ha='right', color='w',
                bbox=dict(fc='k', ec='none', alpha=0.55, pad=1.5))
        print(n, 'eta max (x/D>2.4):', float(np.nanmax(Z[:, XG > 2.4])))
    axs[-1].set_xlabel(r'$x/D$')
    fig.tight_layout()
    fig.colorbar(im, ax=axs, label=r'$\eta$', pad=0.02, fraction=0.04, aspect=30)
    fig.savefig(OUT, dpi=300, bbox_inches='tight')


if __name__ == '__main__':
    main()
