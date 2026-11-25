#!/usr/bin/env python3
"""Shared extractor for fig_wall_state, fig_bl_params and fig_bl_profiles (helper module, not a figure script).

One pass over the plate blocks (6-23, the non-hole blocks) of stat3 batch b006 (193M: b010) of a 3DNS / LES case. Per wall column:
  H, C_f, Re_theta, edge height, U_e   boundary-layer rules of analysis/layer_state/layer_mean.py (EDGE = 3e-2 of the wall
                                       |dU/dz|, + 0.25 mm; compressible delta*, theta; C_f = tau_w / (0.5 rho_e U_e^2),
                                       tau_w from stat3 txz at the wall node; Favre U)
  k01                                  Favre k = 1/2 (R_xx + R_yy + R_zz) at z/D = 0.1 (linear in z between grid levels)
  profiles at x/D = 10                 z/D, U/U_inf, k/U_inf^2 of the column nearest (x/D 10, y/D target), targets
                                       y/D 0, +-1.5, +-3, accepted if the column lies within 0.15 D of the target.
Cache: /mnt/Elements/DNS/analysis/asme2027_figs/bl_extract/<tag>.npz (reused if present; the 193M file is c193b_b010.npz, c193b.npz is the old b006).
usage: python3 _bl_extract.py c20bw|c50b|c100b|c193b      # runs the extraction (stat3 files are read once)
"""
import json
import os
import sys
import time

import numpy as np

REPO = '/home/sk999/Documents/4dns_filmcooled'
sys.path.insert(0, os.path.join(REPO, 'dns', 'python'))
CACHE = '/mnt/Elements/DNS/analysis/asme2027_figs/bl_extract'
CASES = {'c20bw': '/mnt/Tank/DNS/film_cooling/C20bW_BR15',
         'c50b': '/mnt/Tank/DNS/film_cooling/C50b_BR15',
         'c100b': '/mnt/Elements/DNS/film_cooling/C100b_BR15',
         'c193b': '/mnt/Elements/DNS/film_cooling/C193b_BR15'}
BATCH = 6
BATCHES = {'c193b': 10}      # reference rung: batch b010 (10 FT), 2026-11-25; the other rungs keep b006
D, UINF = 0.00775, 86.2812
EDGE = 3e-2
COLS = ['rho', 'rhou', 'rhov', 'rhow', 'rhouu', 'rhovv', 'rhoww', 'rhoT', 'txz']
PT = (0.0, 1.5, -1.5, 3.0, -3.0)
GAS = '/mnt/Elements/DNS/analysis/nearwall_modes/c50b/layout.json'


def mu(T, g):
    return g['mu_c1'] * T ** 1.5 / (T + g['mu_s'])


def block(run, b, g, res, batch=BATCH):
    from tke_io_f import read_stat3_cols
    n = b['n']
    sh = (b['nk'], b['nj'], b['ni'])
    X = np.fromfile(os.path.join(run, f'grid_{n}'), '<f8').reshape(sh + (3,))
    C, meta = read_stat3_cols(run, b, f'_b{batch:03d}', COLS)
    mt = meta['mean_time']
    m = {c: C[c][::-1] / mt for c in COLS}       # index 0 = wall
    del C
    r = m['rho']
    u, v, w = m['rhou'] / r, m['rhov'] / r, m['rhow'] / r
    T = m['rhoT'] / r
    k = 0.5 * (np.maximum(m['rhouu'] / r - u ** 2, 0) + np.maximum(m['rhovv'] / r - v ** 2, 0)
               + np.maximum(m['rhoww'] / r - w ** 2, 0))
    z = X[::-1, 0, 0, 2]
    assert abs(z[0]) < 1e-12 and np.all(np.diff(z) > 0)
    du = np.abs(np.gradient(u, z, axis=0))
    thr = EDGE * du[:3].max(0)
    span = max(1, int(round(2.5e-4 / (z[1] - z[0]))))
    above = du >= thr[None]
    kmax = len(z) - 1 - np.argmax(above[::-1], axis=0)
    ie = np.minimum(len(z) - 1, kmax + span)
    jj, ii = np.meshgrid(np.arange(sh[1]), np.arange(sh[2]), indexing='ij')
    ue, rhoe, Te = u[ie, jj, ii], r[ie, jj, ii], T[ie, jj, ii]
    inside = np.arange(len(z))[:, None, None] <= ie[None]
    f = r * u / (rhoe * ue)[None]
    wt = np.where(inside, 1.0, 0.0)
    dz = np.diff(z)[:, None, None]
    trap = lambda q: np.sum(0.5 * (q[1:] * wt[1:] + q[:-1] * wt[:-1]) * dz * wt[1:], 0)   # noqa: E731
    dstar = trap(1 - f)
    theta = trap(f * (1 - u / ue[None]))
    cf = m['txz'][0] / (0.5 * rhoe * ue ** 2)
    ret = rhoe * ue * theta / mu(Te, g)
    xw, yw = X[-1, ..., 0] / D, X[-1, ..., 1] / D
    # k at z/D = 0.1
    zq = 0.1 * D
    i1 = int(np.searchsorted(z, zq))
    fr = (zq - z[i1 - 1]) / (z[i1] - z[i1 - 1])
    k01 = (1 - fr) * k[i1 - 1] + fr * k[i1]
    for key, a in (('x', xw), ('y', yw), ('cf', cf), ('H', dstar / theta), ('ret', ret), ('ze', z[ie] / D),
                   ('ue', ue / UINF), ('k01', k01 / UINF ** 2)):
        res.setdefault(key, []).append(a.ravel())
    for yt in PT:
        d = np.hypot(xw - 10.0, yw - yt)
        j, i = np.unravel_index(np.argmin(d), d.shape)
        if d[j, i] > 0.15:
            continue
        old = res.setdefault('prof', {}).get(yt)
        if old is None or d[j, i] < old['dist']:
            res['prof'][yt] = dict(dist=float(d[j, i]), x=float(xw[j, i]), y=float(yw[j, i]), z=z / D,
                                   U=u[:, j, i] / UINF, k=k[:, j, i] / UINF ** 2)


def cache_file(tag):
    return os.path.join(CACHE, tag + (f'_b{BATCHES[tag]:03d}' if tag in BATCHES else '') + '.npz')


def extract(tag):
    out = cache_file(tag)
    os.makedirs(CACHE, exist_ok=True)
    run = CASES[tag]
    g = json.load(open(GAS))['gas']
    bi = {b['n']: b for b in json.load(open(os.path.join(run, 'case_info.json')))['blocks']}
    res = {}
    t0 = time.time()
    for n in sorted(bi):
        if bi[n].get('hole'):
            continue
        block(run, bi[n], g, res, BATCHES.get(tag, BATCH))
        print(f'{tag} block {n} [{time.time() - t0:.0f} s]', flush=True)
    prof = res.pop('prof')
    flat = {k: np.concatenate(v) for k, v in res.items()}
    pj = {f'{yt:g}': {k: (v.tolist() if hasattr(v, 'tolist') else v) for k, v in p.items()} for yt, p in prof.items()}
    np.savez(out, prof=json.dumps(pj), **flat)
    return load(tag)


def load(tag):
    """dict of wall arrays plus 'prof' {y target: dict}; extracts if the cache is missing."""
    out = cache_file(tag)
    if not os.path.exists(out):
        return extract(tag)
    d = np.load(out)
    res = {k: d[k] for k in d.files if k != 'prof'}
    res['prof'] = {float(k): {kk: (np.array(vv) if isinstance(vv, list) else vv) for kk, vv in p.items()}
                   for k, p in json.loads(str(d['prof'])).items()}
    return res


if __name__ == '__main__':
    extract(sys.argv[1])
