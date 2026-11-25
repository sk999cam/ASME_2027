#!/usr/bin/env python3
"""Helpers for fig_prt_classes: class sums of Pr_t,eff and the Kohli and Bogard two-component form.

Not a figure script. Two parts.

1. Targeted read of the 193M stat3 batches (C193b, b006 and b010) for the centre plane y = 0, z/D <= ZMAX_KB, to
   obtain the quantities that the plane cache of meff_planes_common.py does not hold:
       uw = u''w'' (Favre), wT = w''T'' (Favre), dUdz = d(U~)/dz, dTdz = d(T~)/dz  (z = wall normal).
   Same operators, same interpolation (8 corners, trilinear weights) and same slab padding as
   meff_planes_common.work_ref. The meff and aeff of the same read are compared with the plane cache as a check.
   Only the k-slabs that hold the points with z/D <= ZMAX_KB are read (about 3.5 GB per batch).
   Cache: /mnt/Elements/DNS/analysis/asme2027_figs/prt_classes/kb_b<batch>.npz

2. Class sums (no thresholds): for a set of points C,
       mu_C    = sum(-a:S*) / sum(2|S*|^2)      = sum(meff S*^2) / sum(S*^2)
       alpha_C = sum(-rho q.gradT) / sum|gradT|^2 = sum(aeff G^2) / sum(G^2)
       Pr_C = mu_C / alpha_C.
   Two-component form: mu_C = sum(-uw dUdz)/sum(dUdz^2), alpha_C = sum(-wT dTdz)/sum(dTdz^2) (rho weights
   both; c_p cancels).
"""
import os
import sys

import numpy as np

sys.path.insert(0, '/home/sk999/Documents/ASME_2027/figures/code')
import meff_planes_common as mc

KB_DIR = '/mnt/Elements/DNS/analysis/asme2027_figs/prt_classes'
ZMAX_KB = 0.3
KB_FIELDS = ['rho', 'uw', 'wT', 'dUdz', 'dTdz', 'meff', 'aeff']


def work_kb(run, b, suffix, k0, k1, sp_n):
    from fd_ops_f import metrics, grad
    from meff_fields import heat_terms
    n = b['n']
    npl = b['nj'] * b['ni']
    C, meta = mc.read_slab(run, b, suffix, mc.COLS, k0, k1)
    m = {c: C[c] / meta['mean_time'] for c in mc.COLS}
    del C
    shp = (k1 - k0, b['nj'], b['ni'])
    with open(os.path.join(run, f'grid_{n}'), 'rb') as fh:
        fh.seek(k0 * npl * 24)
        X = np.fromfile(fh, '<f8', (k1 - k0) * npl * 3).reshape(*shp, 3)
    M = metrics(X)
    buf = np.all(np.stack([m[c] for c in ('txx', 'tyy', 'tzz', 'txy', 'txz', 'tyz')], -1) == 0.0, -1)
    r = m['rho']
    U = np.stack([m['rhou'], m['rhov'], m['rhow']], -1) / r[..., None]
    a = np.empty(r.shape + (3, 3))
    for (i, j), c in mc.PAIRS.items():
        a[..., i, j] = a[..., j, i] = m[c] - r * U[..., i] * U[..., j]
    k = 0.5 * np.trace(a, axis1=-2, axis2=-1) / r
    uw = a[..., 0, 2] / r
    for i in range(3):
        a[..., i, i] -= 2.0 / 3.0 * r * k
    T = m['rhoT'] / r
    fl = np.stack([m['rhouT'], m['rhovT'], m['rhowT']], -1) / r[..., None] - U * T[..., None]
    del m
    gT = grad(np.ascontiguousarray(T), M)
    aeff, cphi, LT, Gt = heat_terms(fl, gT, r, k)
    wT = fl[..., 2].copy()
    dTdz = gT[..., 2].copy()
    del fl, gT
    G = grad(U, M)
    del M
    dUdz = G[..., 0, 2].copy()
    Sm = 0.5 * (G + np.swapaxes(G, -1, -2))
    del G
    tr = np.trace(Sm, axis1=-2, axis2=-1)
    for i in range(3):
        Sm[..., i, i] -= tr / 3
    S2 = np.einsum('...ij,...ij', Sm, Sm)
    aS = np.einsum('...ij,...ij', a, Sm)
    del a, Sm
    with np.errstate(divide='ignore', invalid='ignore'):
        meff = -aS / (2 * S2)
    f = {'rho': r, 'uw': uw, 'wT': wT, 'dUdz': dUdz, 'dTdz': dTdz, 'meff': meff, 'aeff': aeff}
    for v in f.values():
        v[buf] = np.nan
    f = {q: v.astype(np.float32) for q, v in f.items()}
    pos, loc = sp_n
    return meta['mean_time'], (pos, {q: np.asarray(v).ravel()[loc - k0 * npl] for q, v in f.items()})


def build_kb(batch, log=print):
    """Cache {kb_b<batch>.npz}: flat plane indices of the points with z/D <= ZMAX_KB and KB_FIELDS there."""
    import time
    import block_pool
    import plane_locate as pl
    from tke_io_f import load_case
    path = os.path.join(KB_DIR, f'kb_b{batch:03d}.npz')
    if os.path.exists(path):
        return
    os.makedirs(KB_DIR, exist_ok=True)
    t0 = time.time()
    _, bl = load_case(mc.CASE193)
    blocks = [bl[n] for n in sorted(bl)]
    P = [pl.plane('y0', mc.CASE193, blocks, -2 * mc.D, 20 * mc.D, 0.0, 3 * mc.D, n=1100, log=log)]
    pl.release()
    P = P[0]
    zp = P['Y'].ravel() / mc.D
    sel = np.nonzero(zp <= ZMAX_KB)[0]
    keep = np.zeros(zp.size, bool)
    keep[sel] = True
    cor = P['cor'].ravel()
    c = cor
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
    g = {q: np.zeros(cor.size, np.float32) for q in KB_FIELDS}
    mts = set()
    for n, (mt, (pos, vals)) in block_pool.run(work_kb, jobs, lambda n: 1.0 * block_pool.GB + 1500 * npt[n], nw=1):
        mts.add(mt)
        for q in KB_FIELDS:
            g[q][pos] = vals[q]
        log(f'block {n} done [{time.time() - t0:.0f}s]')
    assert len(mts) == 1, mts
    w = P['w']
    out = {'idx': sel, 'mean_time': mts.pop()}
    for q in KB_FIELDS:
        v = (g[q].reshape(P['cor'].shape).astype(np.float64) * w).sum(axis=1)
        out[q] = v[sel].astype(np.float32)
    np.savez_compressed(path, **out)
    log(f'cache written [{time.time() - t0:.0f}s]')


def pr_ls(num_w, num, den_w, den):
    """Least-squares ratio sum(num)/sum(den) over finite points; returns value and the finite mask."""
    m = np.isfinite(num) & np.isfinite(den)
    return num[m].sum() / den[m].sum(), m


if __name__ == '__main__':
    for bt in (10, 6):
        build_kb(bt)
