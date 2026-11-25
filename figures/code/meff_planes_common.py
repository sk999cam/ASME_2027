#!/usr/bin/env python3
"""Shared plane extraction for F14 fig_prt_reference, F17 fig_cross_planes and F18 fig_centre_factors.

Not a figure script. The three figure scripts import it; the first one that runs builds the cache
(/mnt/Elements/DNS/analysis/asme2027_figs/meff_planes/), the others reuse it.

Planes (same physical points on both meshes): centre plane y = 0, x/D -2..20, z/D 0..3 (1100 x 150 points) and the
cross-planes x/D = 4, 10, 20 (y/D -3..3, z/D 0..3, 600 x 300 points).

Reference: 193M 3DNS (C193b), stat3 batch b010 (cumulative mean, first 10 FT), run
/mnt/Elements/DNS/film_cooling/C193b_BR15. Per point, from the Favre stat3 means (methods Sec. 2.4, as
analysis/meff_c50b_c100b/dns_fields.py, whose operators and column names are reused unchanged):
    a = rho R - 2/3 rho k I,  S* = deviatoric sym(grad U~),  |S*| = sqrt(S*:S*),
    mu_eff = -(a:S*)/(2|S*|^2) = rho k^1/2 L_eff cos,  cos = -(a:S*)/(|a||S*|),  L_eff = k^1/2 |b|/|S*|,
    |b| = |a|/(2 rho k),  q = u''T'' (Favre, K m/s; no sub-grid part: a 3DNS),
    alpha_eff/c_p = -rho q.grad T/|grad T|^2,  Pr_t,eff = mu_eff/(alpha_eff/c_p)  (c_p cancels).
Only the k-slab of each block that holds the plane corners (+10 planes each side for the 7-point stencil and the
metrics) is read: the stat3 file is k-major, so the slab is contiguous. Derivatives in the interior of the slab equal
those of the whole block.

RANS 50M: starccm/rans_50M_Lm25b_BR15_FST_ZPG_p0dns (p0dns, iteration 122 686), on the C50b mesh; fields from
/mnt/Tank/DNS/analysis/meff_ransp0_c50b/fields_<n>.npz (analysis/meff_ransp0_c50b): k, mu_t, own |S*|, T;
L = mu_t/(rho k^1/2), cos = 1, alpha_eff/c_p = mu_t/0.9, Pr_t = 0.9.
"""
import os
import sys
import time

import numpy as np

REPO = '/home/sk999/Documents/4dns_filmcooled'
sys.path.insert(0, os.path.join(REPO, 'dns', 'python'))
sys.path.insert(0, os.path.join(REPO, 'analysis', 'meff_rans50_c50b'))
sys.path.insert(0, os.path.join(REPO, 'analysis', 'meff_c50b_c100b'))

D = 0.00775
PITCH = 6 * D
U_INF = 86.2812
DT0 = 300.0 - 199.02            # T0 freestream - T0 coolant (K)
T0_FS, T0_COOL = 300.0, 199.02
CASE193 = '/mnt/Elements/DNS/film_cooling/C193b_BR15'
CASE50 = '/mnt/Tank/DNS/film_cooling/C50b_BR15'
RANS_DIR = '/mnt/Tank/DNS/analysis/meff_ransp0_c50b'
CACHE = '/mnt/Elements/DNS/analysis/asme2027_figs/meff_planes'
BATCH = 10                      # stat3 batch of the 193M reference (b010, 10 FT; brief of 2026-11-25)
CACHE_REF = os.path.join(CACHE, f'ref193_b{BATCH:03d}.npz')
PAD = 10                        # planes of k kept beyond the corners (stencil radius 3, metrics 2 x 3)
XPL = [4, 10, 20]
NAMES = ['y0'] + [f'x{x}' for x in XPL]
REF_FIELDS = ['rho', 'k', 'L', 'cos', 'S', 'T', 'aeff', 'G', 'meff']
RANS_FIELDS = ['rho', 'k', 'L', 'cos', 'S', 'T', 'aeff', 'G', 'meff']

COLS = ['rho', 'rhou', 'rhov', 'rhow', 'rhouu', 'rhovv', 'rhoww', 'rhouv', 'rhouw', 'rhovw', 'rhoT',
        'txx', 'tyy', 'tzz', 'txy', 'txz', 'tyz', 'rhouT', 'rhovT', 'rhowT']
PAIRS = {(0, 0): 'rhouu', (1, 1): 'rhovv', (2, 2): 'rhoww', (0, 1): 'rhouv', (0, 2): 'rhouw', (1, 2): 'rhovw'}


def define_planes(case, blocks, log=print):
    import plane_locate as pl
    P = [pl.plane('y0', case, blocks, -2 * D, 20 * D, 0.0, 3 * D, n=1100, log=log)]
    for x in XPL:
        P.append(pl.plane('x', case, blocks, x * D, z1=3 * D, n=600, log=log))
    return P


def slab_ranges(P, sp, blocks):
    """{block n: (k0, k1)} planes of k to read (None if no corner in the block)."""
    out = {}
    for b in blocks:
        n = b['n']
        npl = b['nj'] * b['ni']
        ks = [sp_q[n][1] // npl for sp_q in sp if len(sp_q[n][1])]
        if not ks:
            continue
        lo = min(int(k.min()) for k in ks)
        hi = max(int(k.max()) for k in ks)
        out[n] = (max(lo - PAD, 0), min(hi + PAD + 1, b['nk']))
    return out


def read_slab(run, b, suffix, names, k0, k1):
    from tke_io_f import stat3_ncol, stat3_meta, CHUNK
    from tke_io import IDX
    n = b['n']
    npl = b['nj'] * b['ni']
    ncol = stat3_ncol(run, b, suffix)
    npt = (k1 - k0) * npl
    cols = [IDX[c] for c in names]
    out = {c: np.empty(npt) for c in names}
    with open(os.path.join(run, f'stat3_{n}{suffix}'), 'rb') as fh:
        fh.seek(k0 * npl * ncol * 8)
        for p0 in range(0, npt, CHUNK):
            m = min(CHUNK, npt - p0)
            a = np.fromfile(fh, '<f8', m * ncol)
            if a.size != m * ncol:
                raise IOError(f'stat3_{n}{suffix}: short read at point {p0}')
            a = a.reshape(m, ncol)
            for c, i in zip(names, cols):
                out[c][p0:p0 + m] = a[:, i]
            del a
    shp = (k1 - k0, b['nj'], b['ni'])
    return {c: v.reshape(shp) for c, v in out.items()}, stat3_meta(run, n, suffix)


def work_ref(run, b, suffix, k0, k1, sp_n, outdir):
    """Reference fields on the slab k0..k1 of block b; returns the corner values of every plane."""
    import plane_locate as pl
    from fd_ops_f import metrics, grad
    from meff_fields import heat_terms
    n = b['n']
    C, meta = read_slab(run, b, suffix, COLS, k0, k1)
    m = {c: C[c] / meta['mean_time'] for c in COLS}
    del C
    shp = (k1 - k0, b['nj'], b['ni'])
    npl = b['nj'] * b['ni']
    with open(os.path.join(run, f'grid_{n}'), 'rb') as fh:
        fh.seek(k0 * npl * 24)
        X = np.fromfile(fh, '<f8', (k1 - k0) * npl * 3).reshape(*shp, 3)
    M = metrics(X)
    buf = np.all(np.stack([m[c] for c in ('txx', 'tyy', 'tzz', 'txy', 'txz', 'tyz')], -1) == 0.0, -1)
    r = m['rho']
    U = np.stack([m['rhou'], m['rhov'], m['rhow']], -1) / r[..., None]
    a = np.empty(r.shape + (3, 3))
    for (i, j), c in PAIRS.items():
        a[..., i, j] = a[..., j, i] = m[c] - r * U[..., i] * U[..., j]
    k = 0.5 * np.trace(a, axis1=-2, axis2=-1) / r
    for i in range(3):
        a[..., i, i] -= 2.0 / 3.0 * r * k
    T = m['rhoT'] / r
    fl = np.stack([m['rhouT'], m['rhovT'], m['rhowT']], -1) / r[..., None] - U * T[..., None]
    del m
    gT = grad(np.ascontiguousarray(T), M)
    aeff, cphi, LT, Gt = heat_terms(fl, gT, r, k)
    del fl, gT
    G = grad(U, M)
    del M
    Sm = 0.5 * (G + np.swapaxes(G, -1, -2))
    del G
    tr = np.trace(Sm, axis1=-2, axis2=-1)
    for i in range(3):
        Sm[..., i, i] -= tr / 3
    S2 = np.einsum('...ij,...ij', Sm, Sm)
    aS = np.einsum('...ij,...ij', a, Sm)
    an = np.sqrt(np.einsum('...ij,...ij', a, a))
    del a, Sm
    with np.errstate(divide='ignore', invalid='ignore'):
        cos = np.nan_to_num(-aS / (an * np.sqrt(S2)), nan=0.0, posinf=0.0, neginf=0.0)
        L = np.nan_to_num(np.sqrt(k) * (an / (2 * r * k)) / np.sqrt(S2), nan=0.0, posinf=0.0, neginf=0.0)
        meff = -aS / (2 * S2)
    f = {'rho': r, 'k': k, 'L': L, 'cos': cos, 'S': np.sqrt(S2), 'T': T, 'aeff': aeff, 'G': Gt, 'meff': meff}
    for v in f.values():
        v[buf] = np.nan
    f = {q: v.astype(np.float32) for q, v in f.items()}
    res = []
    for s in sp_n:
        pos, loc = s
        res.append((pos, {q: np.asarray(v).ravel()[loc - k0 * npl] for q, v in f.items()}))
    return meta['mean_time'], res


def build_ref(log=print, nw=3):
    import block_pool
    import plane_locate as pl
    from tke_io_f import load_case
    path = CACHE_REF
    if os.path.exists(path):
        return
    os.makedirs(CACHE, exist_ok=True)
    t0 = time.time()
    _, bl = load_case(CASE193)
    blocks = [bl[n] for n in sorted(bl)]
    P = define_planes(CASE193, blocks, log)
    pl.release()
    sp = [pl.split(p) for p in P]
    rng = slab_ranges(P, sp, blocks)
    tot = sum((k1 - k0) * blocks[[b['n'] for b in blocks].index(n)]['nj'] * blocks[[b['n'] for b in blocks].index(n)]['ni']
              for n, (k0, k1) in rng.items())
    log(f'blocks read: {sorted(rng)}; {tot / 1e6:.1f} M points, {tot * 74 * 8 / 1e9:.1f} GB [{time.time() - t0:.0f}s]')
    byn = {b['n']: b for b in blocks}
    parts = [[] for _ in P]
    mts = set()
    todo = []
    for n, (k0, k1) in rng.items():
        pf_ = os.path.join(CACHE, f'part_b{BATCH:03d}_{n}.npz')
        if os.path.exists(pf_):
            z = np.load(pf_)
            mts.add(float(z['mt']))
            for q in range(len(P)):
                parts[q].append((z[f'pos{q}'], {f: z[f'{f}{q}'] for f in REF_FIELDS}))
            continue
        todo.append(n)
    jobs = [(n, (CASE193, byn[n], f'_b{BATCH:03d}', rng[n][0], rng[n][1], [s[n] for s in sp], CACHE)) for n in todo]
    npt = {n: (rng[n][1] - rng[n][0]) * byn[n]['nj'] * byn[n]['ni'] for n in rng}
    for n, (mt, res) in block_pool.run(work_ref, jobs, lambda n: 1.0 * block_pool.GB + 1500 * npt[n], nw=nw):
        mts.add(mt)
        d = {'mt': mt}
        for q, (pos, vals) in enumerate(res):
            parts[q].append((pos, vals))
            d[f'pos{q}'] = pos
            for f, v in vals.items():
                d[f'{f}{q}'] = v
        np.savez(os.path.join(CACHE, f'part_b{BATCH:03d}_{n}.npz'), **d)
        log(f'block {n} done [{time.time() - t0:.0f}s]')
    assert len(mts) == 1, mts
    out = {}
    for q, p in enumerate(P):
        F = pl.assemble(p, parts[q], REF_FIELDS)
        out.update({f'{NAMES[q]}_{k}': v for k, v in F.items()})
        out[f'{NAMES[q]}_A'], out[f'{NAMES[q]}_B'] = p['X'], p['Y']
    np.savez_compressed(path, **out, mean_time=mts.pop(), batch=BATCH)
    log(f'reference cached [{time.time() - t0:.0f}s]')


def build_rans(log=print):
    import plane_locate as pl
    from tke_io_f import load_case
    path = os.path.join(CACHE, 'rans50_p0dns.npz')
    if os.path.exists(path):
        return
    os.makedirs(CACHE, exist_ok=True)
    _, bl = load_case(CASE50)
    blocks = [bl[n] for n in sorted(bl)]
    P = define_planes(CASE50, blocks, log)
    pl.release()
    sp = [pl.split(p) for p in P]
    parts = [[] for _ in P]
    for b in blocks:
        n = b['n']
        if not any(len(s[n][0]) for s in sp):
            continue
        z = np.load(os.path.join(RANS_DIR, f'fields_{n}.npz'))
        f = {'rho': z['rho_r'], 'k': z['k_r'], 'L': z['L_r'], 'S': z['S_r'], 'T': z['T_r'], 'aeff': z['aeff_r'],
             'G': z['G_r'], 'meff': z['meff_r']}
        f['cos'] = np.ones_like(f['k'])
        for q, s in enumerate(sp):
            parts[q].append(pl.take(s[n], f))
        log(f'RANS block {n}')
    out = {}
    for q, p in enumerate(P):
        F = pl.assemble(p, parts[q], RANS_FIELDS)
        out.update({f'{NAMES[q]}_{k}': v for k, v in F.items()})
        out[f'{NAMES[q]}_A'], out[f'{NAMES[q]}_B'] = p['X'], p['Y']
    np.savez_compressed(path, **out)


def load(which):
    """which: 'ref' (193M b010) or 'rans' (RANS 50M p0dns); builds the cache when missing."""
    if which == 'ref':
        build_ref()
        return dict(np.load(CACHE_REF))
    build_rans()
    return dict(np.load(os.path.join(CACHE, 'rans50_p0dns.npz')))


def eta(T):
    return (T0_FS - T) / (T0_FS - T0_COOL)


if __name__ == '__main__':
    build_rans()
    build_ref()
