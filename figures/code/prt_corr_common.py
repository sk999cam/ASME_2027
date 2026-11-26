#!/usr/bin/env python3
"""Targeted reads of C193b stat3_b006 and b010 for two centre-plane boxes: film x/D 5-12, z/D 0.02-1.5 and jet x/D 0-3,
z/D 0-3 (about 4 GB each).

Gives, per plane point (trilinear 8-corner interpolation as meff_planes_common), Favre fluxes u''T'', v''T'',
w''T'', the gradient of T~, the rms values u'', w'', T'' (Favre) and the stress u''w'' with dU~/dz, to obtain the
correlation coefficients R_wT, R_uT and the share of q.gradT carried by each direction (fig_prt_classes.py and
the near-wall sampling analysis). Cache: /mnt/Elements/DNS/analysis/asme2027_figs/prt_classes/corr_<film|jet>_b<batch>.npz
"""
import os
import sys
import time

import numpy as np

sys.path.insert(0, '/home/sk999/Documents/ASME_2027/figures/code')
import meff_planes_common as mc

CORR_DIR = '/mnt/Elements/DNS/analysis/asme2027_figs/prt_classes'
D = mc.D
# tag: (plane kind, plane arguments, n points, box in plane coordinates (A/D0, A/D1, B/D0, B/D1))
#   y0: A = x, B = z;  x10: A = y, B = z.  'full' covers the whole centre plane, x/D 0-20.
SPECS = {'film': ('y0', (-2 * D, 20 * D, 0.0, 3 * D), 1100, (5.0, 12.0, 0.02, 1.5)),
         'jet': ('y0', (-2 * D, 20 * D, 0.0, 3 * D), 1100, (0.0, 3.0, 0.0, 3.0)),
         'full': ('y0', (-2 * D, 20 * D, 0.0, 3 * D), 1100, (-2.0, 20.0, 0.0, 3.0)),
         'x10': ('x', (10 * D,), 600, (-3.0, 3.0, 0.0, 3.0)),
         'x20': ('x', (20 * D,), 600, (-3.0, 3.0, 0.0, 3.0)),
         'x4': ('x', (4 * D,), 600, (-3.0, 3.0, 0.0, 3.0)),
         'tu': ('y0', (-5.5 * D, -4.5 * D, 2.0 * D, 3.0 * D), 50, (-6.0, -4.0, 1.9, 3.1))}
COLS = ['rho', 'rhou', 'rhov', 'rhow', 'rhouu', 'rhovv', 'rhoww', 'rhouv', 'rhouw', 'rhoT', 'rhoTT',
        'rhouT', 'rhovT', 'rhowT', 'txx', 'tyy', 'tzz', 'txy', 'txz', 'tyz']
FIELDS = ['rho', 'uT', 'vT', 'wT', 'dTdx', 'dTdy', 'dTdz', 'urms', 'vrms', 'wrms', 'Trms', 'uw', 'dUdz', 'uv', 'dUdy',
          'k']


def work_corr(run, b, suffix, k0, k1, sp_n):
    from fd_ops_f import metrics, grad
    n = b['n']
    npl = b['nj'] * b['ni']
    C, meta = mc.read_slab(run, b, suffix, COLS, k0, k1)
    m = {c: C[c] / meta['mean_time'] for c in COLS}
    del C
    shp = (k1 - k0, b['nj'], b['ni'])
    with open(os.path.join(run, f'grid_{n}'), 'rb') as fh:
        fh.seek(k0 * npl * 24)
        X = np.fromfile(fh, '<f8', (k1 - k0) * npl * 3).reshape(*shp, 3)
    M = metrics(X)
    buf = np.all(np.stack([m[c] for c in ('txx', 'tyy', 'tzz', 'txy', 'txz', 'tyz')], -1) == 0.0, -1)
    r = m['rho']
    U = np.stack([m['rhou'], m['rhov'], m['rhow']], -1) / r[..., None]
    T = m['rhoT'] / r
    uu = m['rhouu'] / r - U[..., 0] ** 2
    vv = m['rhovv'] / r - U[..., 1] ** 2
    ww = m['rhoww'] / r - U[..., 2] ** 2
    TT = m['rhoTT'] / r - T ** 2
    uw = m['rhouw'] / r - U[..., 0] * U[..., 2]
    uv = m['rhouv'] / r - U[..., 0] * U[..., 1]
    fl = np.stack([m['rhouT'], m['rhovT'], m['rhowT']], -1) / r[..., None] - U * T[..., None]
    del m
    gT = grad(np.ascontiguousarray(T), M)
    G = grad(U, M)
    f = {'rho': r, 'uT': fl[..., 0], 'vT': fl[..., 1], 'wT': fl[..., 2], 'dTdx': gT[..., 0], 'dTdy': gT[..., 1],
         'dTdz': gT[..., 2], 'urms': np.sqrt(np.maximum(uu, 0)), 'wrms': np.sqrt(np.maximum(ww, 0)),
         'Trms': np.sqrt(np.maximum(TT, 0)), 'uw': uw, 'dUdz': G[..., 0, 2], 'vrms': np.sqrt(np.maximum(vv, 0)),
         'uv': uv, 'dUdy': G[..., 0, 1], 'k': 0.5 * (uu + vv + ww)}
    for v in f.values():
        v[buf] = np.nan
    f = {q: np.ascontiguousarray(v).astype(np.float32) for q, v in f.items()}
    pos, loc = sp_n
    return meta['mean_time'], (pos, {q: v.ravel()[loc - k0 * npl] for q, v in f.items()})


def build_corr(batch, tag, log=print):
    import block_pool
    import plane_locate as pl
    from tke_io_f import load_case
    path = os.path.join(CORR_DIR, f'corr_{tag}_b{batch:03d}.npz')
    if os.path.exists(path):
        return
    t0 = time.time()
    _, bl = load_case(mc.CASE193)
    blocks = [bl[n] for n in sorted(bl)]
    kind, args, npts, BOX = SPECS[tag]
    if kind == 'y0':
        P = pl.plane('y0', mc.CASE193, blocks, *args, n=npts, log=log)
    else:
        P = pl.plane('x', mc.CASE193, blocks, args[0], z1=3 * mc.D, n=npts, log=log)
    pl.release()
    x, z = P['X'].ravel() / mc.D, P['Y'].ravel() / mc.D
    keep = (x >= BOX[0]) & (x <= BOX[1]) & (z >= BOX[2]) & (z <= BOX[3])
    sel = np.nonzero(keep)[0]
    c = P['cor'].ravel()
    byn = {b['n']: b for b in blocks}
    jobs, npt = [], {}
    for i, n in enumerate(P['order']):
        lo, hi = P['off'][i], P['off'][i + 1]
        pos = np.nonzero((c >= lo) & (c < hi))[0]
        pos = pos[keep[pos // 8]]
        if not len(pos):
            continue
        b = byn[n]
        npl = b['nj'] * b['ni']
        loc = c[pos] - lo
        ks = loc // npl
        k0, k1 = max(int(ks.min()) - mc.PAD, 0), min(int(ks.max()) + mc.PAD + 1, b['nk'])
        npt[n] = (k1 - k0) * npl
        jobs.append((n, (mc.CASE193, b, f'_b{batch:03d}', k0, k1, (pos, loc))))
    log(f'batch {batch}: {len(jobs)} blocks, {sum(npt.values()) / 1e6:.1f} M points')
    g = {q: np.zeros(c.size, np.float32) for q in FIELDS}
    mts = set()
    for n, (mt, (pos, vals)) in block_pool.run(work_corr, jobs, lambda n: 1.0 * block_pool.GB + 1500 * npt[n], nw=1):
        mts.add(mt)
        for q in FIELDS:
            g[q][pos] = vals[q]
        log(f'block {n} done [{time.time() - t0:.0f}s]')
    assert len(mts) == 1, mts
    out = {'idx': sel, 'mean_time': mts.pop()}
    for q in FIELDS:
        v = (g[q].reshape(P['cor'].shape).astype(np.float64) * P['w']).sum(axis=1)
        out[q] = v[sel].astype(np.float32)
    np.savez_compressed(path, **out)
    log(f'cache written [{time.time() - t0:.0f}s]')


if __name__ == '__main__':
    jobs = [a.split(':') for a in sys.argv[1:]] or [('tu', '10'), ('full', '10'), ('full', '6'), ('x10', '10'), ('x10', '6')]
    for tag, bt in jobs:
        build_corr(int(bt), tag)
