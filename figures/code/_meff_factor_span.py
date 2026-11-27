"""Span means of the point-wise factors of mu = rho k^(1/2) L cos(theta), both rungs, on the grid of span_fields.npz
(x/D -2..20 step 0.05, z/D 0..2 step 0.01). Same interpolation as span_compute.py / _meff_strain_span.py (Delaunay in (x, y) on
the common z levels, linear in z, points further than 0.25 D from a data point left empty); each field is averaged over the
points of the span line that are finite for that field (plain mean, no weighting, no mask).
Lengths: RANS L_R = mu_t/(rho_R k_R^(1/2)); 193M L_eff = L of the point file (mu_eff = rho k^(1/2) L_eff cos).
Fields (SI, divide by U_inf / D outside): rhoR, k12R, LR, rhoD, k12D, LD, cD, LalD (= L_eff cos).
Cached in the scratchpad (SCR); the meff_span files carry only mu, S_mean, S_frac_gt005."""
import os

import numpy as np
from scipy.interpolate import LinearNDInterpolator
from scipy.spatial import Delaunay, cKDTree

import _meff_pointwise_data as M

CACHE = M.SCR + 'span_factors.npz'
NAMES = ['rhoR', 'k12R', 'LR', 'rhoD', 'k12D', 'LD', 'cD', 'LalD']


def get():
    if os.path.exists(CACHE):
        c = np.load(CACHE)
        return c['x'], c['z'], {k: c[k] for k in NAMES}, c['n']
    d = np.load(M.SCR + 'pts.npz')
    X = np.load(M.SCR + 'xyz.npz')['xyz']
    g = lambda n: d[n].astype(np.float64)  # noqa: E731
    with np.errstate(all='ignore'):
        F = dict(rhoR=g('rans_rho'), k12R=np.sqrt(g('rans_k')), LR=g('rans_mut') / (g('rans_rho') * np.sqrt(g('rans_k'))),
                 rhoD=g('c193_rho'), k12D=np.sqrt(g('c193_k')), LD=g('c193_L'), cD=g('c193_cos'), LalD=g('c193_L') * g('c193_cos'))
    zl = np.unique(np.round(X[d['blk'] == 7, 2], 5))
    xg = np.arange(-2, 20 + 1e-9, 0.05)
    yg = np.arange(-3, 3 + 1e-9, 0.05)
    zf = np.arange(0, 2.0 + 1e-9, 0.01)
    XX, YY = np.meshgrid(xg, yg, indexing='ij')
    tgt = np.c_[XX.ravel(), YY.ravel()]
    lev = np.abs(X[:, 2][:, None] - zl[None, :]).argmin(1)
    okz = (X[:, 2] >= -1e-6) & (np.abs(X[:, 2] - zl[lev]) < 2e-4) & (X[:, 0] > -3.5) & (X[:, 0] < 21)
    L = {k: np.full((len(zl), len(xg), len(yg)), np.nan) for k in NAMES}
    for iz in range(len(zl)):
        sel = okz & (lev == iz)
        P = X[sel][:, :2]
        for key in NAMES:
            vv = F[key][sel]
            f = np.isfinite(vv)
            pp, vf = P[f], vv[f]
            _, u = np.unique(np.round(pp, 6), axis=0, return_index=True)
            pp, vf = pp[u], vf[u]
            val = LinearNDInterpolator(Delaunay(pp), vf)(tgt)
            val[cKDTree(pp).query(tgt)[0] > 0.25] = np.nan
            L[key][iz] = val.reshape(len(xg), len(yg))
    i1 = np.clip(np.searchsorted(zl, zf), 1, len(zl) - 1)
    i0 = i1 - 1
    w = ((zf - zl[i0]) / (zl[i1] - zl[i0]))[:, None, None]
    out, n = {}, None
    for key in NAMES:
        a = (1 - w) * L[key][i0] + w * L[key][i1]
        if n is None:
            n = np.isfinite(a).sum(2)
        with np.errstate(all='ignore'):
            out[key] = np.nanmean(a, axis=2)
    np.savez(CACHE, x=xg, z=zf, n=n, **out)
    return xg, zf, out, n
