"""F8 fig_flow_features: mean streamwise vorticity on cross-planes and the instantaneous shear-layer roll-up (Results 3).

Usage: python3 fig_flow_features.py      (writes /home/sk999/Documents/ASME_2027/figures/fig_flow_features.png)

(a)-(c) RANS 50M and (d)-(f) 193M 3DNS (b010): mean omega_x D/U_inf (bwr, limits +-OMAX, same in all six) with in-plane mean
streamlines (V, W) on the planes x/D = 2, 4, 10; y/D -3..3 (one pitch), z/D 0..ZMAX.
(g) instantaneous 193M centre plane y = 0, x/D -2..8: temperature (T0,fs - T)/(T0,fs - T0,c) (jet) with the hole outline
and the vorticity-roll-up of the jet shear layer (Kelvin-Helmholtz).

Data
  193M: /mnt/Elements/DNS/analysis/flow_features_c193b/n3_crossplanes_c193b.npz (stat3 b010 means; U, V, W in m/s, omx in 1/s)
        and kh_frame_c193b.npz (instantaneous frame, producer resample_n3n4.py of 4dns_filmcooled).
  RANS 50M: rans50_exports/rans50_base_cells_grad.csv (STAR cell export, velocity gradients). Columns needed only:
        V, W, dVdz, dWdy, X, Y, Z. Cells within +-SLAB/D of each plane are interpolated (linear, Delaunay) on the 193M y-z grid;
        omega_x = dW/dy - dV/dz. Cached in /mnt/Elements/DNS/analysis/asme2027_figs/flow_features/rans50_planes.npz.
        Run `python3 fig_flow_features.py extract` once (long read) to build the cache.
Status: complete when the plot is checked.
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from cycler import cycler

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, '/home/sk999/Documents/3DNS_FP')
from functions import set_mpl_params  # noqa: E402

OUT = '/home/sk999/Documents/ASME_2027/figures/fig_flow_features.png'
FF = '/mnt/Elements/DNS/analysis/flow_features_c193b'
CACHE = '/mnt/Elements/DNS/analysis/asme2027_figs/flow_features'
RANS_CSV = '/mnt/Elements/DNS/analysis/rans50_exports/rans50_base_cells_grad.csv'
OKABE_ITO = ['#000000', '#E69F00', '#56B4E9', '#009E73', '#F0E442', '#0072B2', '#D55E00', '#CC79A7']
D = 7.75e-3
UINF = 86.28
PLANES = (2, 4, 10)
SLAB = 0.06          # half thickness of the RANS slab, in D
OMAX = 1.0
ZMAX = 1.6
XKH = (-2.0, 8.0)


def extract_rans():
    """Read the needed columns of the 30 GB export in chunks and keep the cells of the three slabs."""
    import pandas as pd
    cols = ['Velocity[j] (m/s)', 'Velocity[k] (m/s)', 'dVdz', 'dWdy', 'X (m)', 'Y (m)', 'Z (m)']
    keep = {p: [] for p in PLANES}
    n = 0
    for ch in pd.read_csv(RANS_CSV, usecols=cols, chunksize=2_000_000, dtype=np.float64):
        x = ch['X (m)'].to_numpy() / D
        for p in PLANES:
            m = np.abs(x - p) < SLAB
            if m.any():
                keep[p].append(ch.loc[m, cols].to_numpy())
        n += len(ch)
        print(n, flush=True)
    os.makedirs(CACHE, exist_ok=True)
    out = {}
    for p in PLANES:
        a = np.vstack(keep[p])
        out['cells_x%d' % p] = a
        print(p, a.shape, flush=True)
    np.savez(os.path.join(CACHE, 'rans50_cells_slabs.npz'), **out)


def rans_planes(y, z):
    """Interpolate the slab cells on the (y, z) grid of the 193M planes (cached)."""
    from scipy.interpolate import LinearNDInterpolator
    f = os.path.join(CACHE, 'rans50_planes.npz')
    if os.path.exists(f):
        Z = np.load(f)
        return {k: Z[k] for k in Z.files}
    S = np.load(os.path.join(CACHE, 'rans50_cells_slabs.npz'))
    Y, Zg = np.meshgrid(y, z, indexing='ij')
    out = {}
    for p in PLANES:
        a = S['cells_x%d' % p]
        v, w, dvdz, dwdy, yy, zz = a[:, 0], a[:, 1], a[:, 2], a[:, 3], a[:, 5] / D, a[:, 6] / D
        pts = np.column_stack([yy, zz])
        sel = (np.abs(yy) < 3.3) & (zz < 2.2)
        pts = pts[sel]
        for name, q in (('om', (dwdy - dvdz) * D / UINF), ('V', v), ('W', w)):
            I = LinearNDInterpolator(pts, q[sel])
            out['%s_x%d' % (name, p)] = I(Y, Zg).astype(np.float32)
    np.savez(f, **out)
    return out


def peak(om, y, z, sgn):
    """Position and value of the extremum of sign sgn (+1 / -1) of omega_x in 0.2 < z/D < ZMAX, 0 < y/D < 3 or -3 < y/D < 0."""
    Y, Zg = np.meshgrid(y, z, indexing='ij')
    m = (Zg > 0.2) & (Zg < 1.2) & (np.abs(Y) < 1.4)
    a = np.where(m & np.isfinite(om), om * sgn, -np.inf)
    i = np.unravel_index(np.argmax(a), a.shape)
    return om[i], y[i[0]], z[i[1]]


def main():
    N = np.load(os.path.join(FF, 'n3_crossplanes_c193b.npz'), allow_pickle=True)
    y, z = N['y'], N['z']
    R = rans_planes(y, z)
    set_mpl_params(0.5, grid_on=False)
    plt.rcParams['axes.prop_cycle'] = cycler(color=OKABE_ITO)
    fig = plt.figure(figsize=(11.0, 6.4))
    gs = fig.add_gridspec(2, 3, hspace=0.08, wspace=0.12, left=0.07, right=0.88, top=1.0, bottom=0.64)
    gsg = fig.add_gridspec(1, 1, left=0.07, right=0.88, top=0.54, bottom=0.08)
    zi = z <= ZMAX + 0.01
    zc = z[zi]
    letters = 'abcdef'
    res = []
    for r, rung in enumerate(('RANS 50M', '193M 3DNS')):
        for c, p in enumerate(PLANES):
            ax = fig.add_subplot(gs[r, c])
            ax.set_anchor('S' if r == 0 else 'N')
            if r == 0:
                om, V, W = R['om_x%d' % p], R['V_x%d' % p], R['W_x%d' % p]
            else:
                om = N['omx_x%d' % p] * D / UINF
                V, W = N['U_x%d' % p][..., 1], N['U_x%d' % p][..., 2]
            om, V, W = om[:, zi], V[:, zi], W[:, zi]
            im = ax.pcolormesh(y, zc, om.T, cmap='bwr', vmin=-OMAX, vmax=OMAX, shading='nearest', rasterized=True)
            Vm, Wm = np.nan_to_num(V.T), np.nan_to_num(W.T)
            ax.streamplot(y, zc, Vm, Wm, color='0.3', linewidth=0.35, density=0.8, arrowsize=0.4)
            ax.set_xlim(-3, 3)
            ax.set_ylim(0, ZMAX)
            ax.set_aspect('equal')
            ax.text(0.02, 0.95, '(%s) %s, $x/D=%d$' % (letters[3 * r + c], rung, p), transform=ax.transAxes, va='top', fontsize=8,
                    bbox=dict(fc='w', ec='none', alpha=0.85, pad=1.0), zorder=5)
            if c == 0:
                ax.set_ylabel(r'$z/D$')
            else:
                ax.set_yticklabels([])
            if r == 1:
                ax.set_xlabel(r'$y/D$')
            else:
                ax.set_xticklabels([])
            for sgn in (+1, -1):
                res.append((rung, p, sgn) + peak(om, y, zc, sgn))
            if r == 1 and c == 1:
                bb = dict(fc='w', ec='none', alpha=0.85, pad=1.0)
                ax.annotate('counter-rotating vortex pair', xy=(-0.8, 0.38), xytext=(0.0, 1.0), fontsize=7, ha='center',
                            bbox=bb, zorder=6, arrowprops=dict(arrowstyle='->', lw=0.8, color='k'))
                ax.annotate('', xy=(0.65, 0.38), xytext=(0.0, 0.95), zorder=6,
                            arrowprops=dict(arrowstyle='->', lw=0.8, color='k'))
    cax = fig.add_axes([0.895, 0.68, 0.012, 0.30])
    cb = fig.colorbar(im, cax=cax, extend='both')
    cb.set_label(r'$\overline{\omega}_x D/U_\infty$')
    # (g) instantaneous frame
    K = np.load(os.path.join(FF, 'kh_frame_c193b.npz'), allow_pickle=True)
    x, zk, T = K['x'], K['z'], K['T']
    eta = (300.0 - T) / (300.0 - 199.02)
    ax = fig.add_subplot(gsg[0, 0])
    im2 = ax.pcolormesh(x, zk, eta.T, cmap='jet', vmin=0, vmax=1, shading='nearest', rasterized=True)
    ol = K['outline']
    ax.set_facecolor('0.85')
    ax.plot(ol[:, 0], ol[:, 1], 'k-', lw=0.8, zorder=3)
    ax.set_xlim(*XKH)
    ax.set_ylim(-1.0, 2.5)
    ax.set_aspect('equal')
    ax.set_xlabel(r'$x/D$')
    ax.set_ylabel(r'$z/D$')
    ax.text(0.01, 0.95, '(g) 193M 3DNS, instantaneous, $y=0$', transform=ax.transAxes, va='top', fontsize=8, color='w')
    cax2 = fig.add_axes([0.895, 0.10, 0.012, 0.36])
    cb2 = fig.colorbar(im2, cax=cax2)
    cb2.set_label(r'$(T_{0,fs}-T)/(T_{0,fs}-T_{0,c})$')
    fig.savefig(OUT, dpi=300)
    print('t =', float(K['t']), 's, tU/D =', float(K['tU_D']), 'record', int(K['record']), 'step', int(K['step']))
    for r_ in res:
        print('%-10s x/D=%2d sign %+d  omega=%+.2f  y/D=%.2f  z/D=%.2f' % r_)


if __name__ == '__main__':
    if len(sys.argv) > 1 and sys.argv[1] == 'extract':
        extract_rans()
    else:
        main()
