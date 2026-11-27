#!/usr/bin/env python3
"""F8 fig_wall_state: plan view of the wall layer at z/D = 0.1, x/D -3 to 20, one pitch (y/D -3 to 3).

Top row, near-wall k/U_inf^2 at z/D = 0.1: (a) RANS 50M, (b) 20M LES, (c) 193M 3DNS.
Bottom row, intermittency gamma at z/D = 0.1: (d) RANS 50M (model variable of SST + gamma-Re_theta), (e) 20M LES,
(f) 193M 3DNS. Complete: all six panels are computed for their own rung.

Data
  RANS 50M  starccm rans_50M_Lm25b_BR15_FST_ZPG_p0dns, RC-S1 cell export (cells_rans_p0dns_full.csv, 53.1 M cells):
            the cell layer at z/D 0.098 (0.094 < z/D < 0.102; the layer nearest 0.1), x/D -3.6 to 21, linearly interpolated in plan;
            modelled k and gamma.
  20M LES   C20bW stat3 b006 (k) and plane sampler z010, records of the b006 window (gamma).
  193M 3DNS C193b stat3 b010 (k, pictures_job18/stat3_b010 means; helper _bl_extract.py) and plane sampler z010, records of
            the b010 window (gamma).
  k: Favre k = 1/2 (R_xx + R_yy + R_zz) of the plate blocks at z/D 0.1 (helper _bl_extract.py, cache
  /mnt/Elements/DNS/analysis/asme2027_figs/bl_extract/).
  gamma (3DNS / LES): high-pass spanwise stress d = (v - <v>_T)^2 / U_inf^2 (T = 2 D/U running mean in time, then 1 D/U
  running mean), time fraction above the Otsu threshold of log10 d pooled over the 193M intermittent band |y/D| 2 to 3,
  x/D 0 to 20 (one value for both runs); helper
  _gamma_extract_v.py, cache /mnt/Elements/DNS/analysis/asme2027_figs/gamma_v/ (definition: planning/maths_audit/intermittency).
  The old |du/dt| detector (helper _gamma_extract.py) is only read for the printed old-versus-new comparison.
The 3DNS / LES stat3 plate blocks exclude the hole, so the hole exit region (x/D about -1 to 2.5, |y/D| < 1.5) is blank
in the k maps of the 20M LES and 193M. The map stops at x/D 20 (x/D 20-25 is a buffer region, not shown). RANS k is modelled, 3DNS / LES k is resolved.
usage: python3 fig_wall_state.py
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import LogNorm, Normalize  # noqa: E402
from cycler import cycler  # noqa: E402
from scipy.interpolate import LinearNDInterpolator  # noqa: E402
from scipy.spatial import cKDTree  # noqa: E402

sys.path.insert(0, '/home/sk999/Documents/3DNS_FP')
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from functions import set_mpl_params  # noqa: E402
import _bl_extract as BX  # noqa: E402
import _gamma_extract as GX  # noqa: E402      # old detector, printed comparison only
import _gamma_extract_v as GV  # noqa: E402

OUT = '/home/sk999/Documents/ASME_2027/figures/fig_wall_state.png'
SRC = '/mnt/Tank/STARCCM/film_cooling/rans_50M_Lm25b_BR15_FST_ZPG_p0dns/cells_rans_p0dns_full.csv'
D, UINF = 0.00775, 86.2812
XLIM = (-3.0, 20.0)       # owner rule 2026-11-26: nothing beyond x/D 20 is shown
CACHE = BX.CACHE
OKABE_ITO = ['#000000', '#E69F00', '#56B4E9', '#009E73', '#F0E442', '#0072B2', '#D55E00', '#CC79A7']


def rans_slab():
    f = os.path.join(CACHE, 'rans50_z010_cells.npz')
    if os.path.exists(f):
        d = np.load(f)
        return {k: d[k] for k in d.files}
    import pandas as pd
    use = ['Turbulent Kinetic Energy (J/kg)', 'Intermittency', 'X (m)', 'Y (m)', 'Z (m)']
    acc = []
    for ch in pd.read_csv(SRC, usecols=use, chunksize=2_000_000, dtype=np.float64):
        x, y, z = ch['X (m)'].values / D, ch['Y (m)'].values / D, ch['Z (m)'].values / D
        m = (z > 0.08) & (z < 0.12) & (x > -3.6) & (x < 21.0) & (np.abs(y) < 3.2)
        if m.any():
            acc.append(np.column_stack([x[m], y[m], z[m], ch['Turbulent Kinetic Energy (J/kg)'].values[m] / UINF ** 2,
                                        ch['Intermittency'].values[m]]))
    a = np.concatenate(acc)
    d = dict(x=a[:, 0], y=a[:, 1], z=a[:, 2], k=a[:, 3], g=a[:, 4])
    os.makedirs(CACHE, exist_ok=True)
    np.savez(f, **d)
    return d


def regrid(x, y, v, gx, gy, mask_d=0.12):
    """linear interpolation on the plan grid, blank where no data point lies within mask_d (D)."""
    ok = np.isfinite(v)
    pts = np.column_stack([x[ok], y[ok]])
    # periodic copies so the pitch edges are filled
    for s in (-6.0, 6.0):
        pts = np.vstack([pts, np.column_stack([x[ok], y[ok] + s])])
    vv = np.tile(v[ok], 3)
    sel = (np.abs(pts[:, 1]) < 3.3)
    pts, vv = pts[sel], vv[sel]
    f = LinearNDInterpolator(pts, vv)
    GX, GY = np.meshgrid(gx, gy, indexing='ij')
    out = f(GX, GY)
    dist, _ = cKDTree(pts).query(np.column_stack([GX.ravel(), GY.ravel()]))
    out[dist.reshape(out.shape) > mask_d] = np.nan
    return out


def main():
    set_mpl_params(0.7, grid_on=False)
    plt.rcParams['axes.prop_cycle'] = cycler(color=OKABE_ITO)
    gx = np.arange(XLIM[0], XLIM[1] + 1e-9, 0.04)
    gy = np.arange(-3.0, 3.0 + 1e-9, 0.04)

    # --- k maps
    kmap = {}
    r = rans_slab()
    sel = (r['z'] > 0.094) & (r['z'] < 0.102)          # the prism layer at z/D 0.098 (layers: 0.083, 0.090, 0.098, 0.105, 0.113)
    r = {k: v[sel] for k, v in r.items()}
    kmap['RANS 50M'] = regrid(r['x'], r['y'], r['k'], gx, gy, 0.15)
    gmap_rans = regrid(r['x'], r['y'], r['g'], gx, gy, 0.15)
    for tag, lab in (('c20bw', '20M LES'), ('c193b', '193M 3DNS')):
        w = BX.load(tag)
        kmap[lab] = regrid(w['x'], w['y'], w['k01'], gx, gy, 0.15)
    gmap, gold, glin = {}, {}, {}
    for tag, lab in (('c20bw', '20M LES'), ('c193b', '193M 3DNS')):
        G = GV.load_z010(tag)['gamma']                 # policy threshold; (1176, 301), x -3..20.5, y -3..3 by 0.02
        ix = np.clip(np.round((gx + 3.0) / 0.02).astype(int), 0, GV.NXR - 1)
        iy = np.clip(np.round((gy + 3.0) / 0.02).astype(int), 0, GV.NY - 1)
        gmap[lab] = G[np.ix_(ix, iy)]
        Go = GX.load_z010(tag)['gamma'][0][np.ix_(ix, iy)]
        gold[lab] = Go
        glin[lab] = GV.gamma_at(GV.load(tag, 'z010')['frac'][1], GV.policy_logth(tag, 2.0, 'iii'))[np.ix_(ix, iy)]

    fig, ax = plt.subplots(2, 3, figsize=(15.5, 5.6), sharex=True, sharey=True)
    nk = LogNorm(1e-4, 3e-2)
    names = ['RANS 50M', '20M LES', '193M 3DNS']
    for j, nm in enumerate(names):
        im = ax[0, j].pcolormesh(gx, gy, kmap[nm].T, cmap='inferno', norm=nk, shading='auto', rasterized=True)
    cb = fig.colorbar(im, ax=ax[0, :], pad=0.01, fraction=0.025)
    cb.set_label(r'$k/U_\infty^2$')
    gm = {'RANS 50M': gmap_rans, '20M LES': gmap['20M LES'], '193M 3DNS': gmap['193M 3DNS']}
    for j, nm in enumerate(names):
        g2 = ax[1, j].pcolormesh(gx, gy, gm[nm].T, cmap='inferno', norm=Normalize(0, 1), shading='auto', rasterized=True)
    cb2 = fig.colorbar(g2, ax=ax[1, :], pad=0.01, fraction=0.025)
    cb2.set_label(r'$\gamma$')
    labs = {(0, 0): '(a) RANS 50M', (0, 1): '(b) 20M LES', (0, 2): '(c) 193M 3DNS',
            (1, 0): '(d) RANS 50M', (1, 1): '(e) 20M LES', (1, 2): '(f) 193M 3DNS'}
    for (i, j), t in labs.items():
        a = ax[i, j]
        a.text(0.012, 0.88, t, transform=a.transAxes, color='w', ha='left', va='center')
        a.set_xlim(*XLIM)
        a.set_ylim(-3, 3)
        a.set_yticks([-3, 0, 3])
        a.set_xticks([0, 5, 10, 15, 20])
        a.set_aspect('auto')
    for a in ax[1, :]:
        a.set_xlabel(r'$x/D$')
    for i in (0, 1):
        ax[i, 0].set_ylabel(r'$y/D$')
    fig.savefig(OUT, dpi=300, bbox_inches='tight')
    for nm, m in kmap.items():
        print(nm, 'k/U^2 median', np.nanmedian(m), 'p5', np.nanpercentile(m, 5), 'p95', np.nanpercentile(m, 95))
    for nm, m in gm.items():
        print(nm, 'plate-mean gamma (x/D -3..20, |y/D| < 3): band policy', np.nanmean(m),
              ('| line policy %.4f | old |du/dt| %.4f' % (np.nanmean(glin[nm]), np.nanmean(gold[nm]))) if nm in gold else '')


if __name__ == '__main__':
    main()
