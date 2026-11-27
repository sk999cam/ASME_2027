"""Span mean of the point-wise 193M strain magnitude |S*| D/U_inf (and the span fraction of points with |S*| D/U_inf > S0) on the
grid of span_fields.npz (x/D -2..20 step 0.05, z/D 0..2 step 0.01). Same interpolation as span_compute.py (Delaunay in (x, y) on the
common z levels, linear in z, points further than 0.25 D from a data point left empty). Cached in the scratchpad (SCR).
Used for the low-strain mask of fig_meff_maps.py. Reads the point file pts.npz only."""
import os

import numpy as np
from scipy.interpolate import LinearNDInterpolator
from scipy.spatial import Delaunay, cKDTree

import _meff_pointwise_data as M

CACHE = M.SCR + 'span_strain.npz'


def get(s0=0.05):
    if os.path.exists(CACHE):
        c = np.load(CACHE)
        if float(c['s0']) == s0:
            return c['x'], c['z'], c['Smean'], c['frac']
    d = np.load(M.SCR + 'pts.npz')
    X = np.load(M.SCR + 'xyz.npz')['xyz']
    with np.errstate(all='ignore'):
        S = d['c193_S'].astype(np.float64) * M.D / M.UINF
    zl = np.unique(np.round(X[d['blk'] == 7, 2], 5))
    xg = np.arange(-2, 20 + 1e-9, 0.05)
    yg = np.arange(-3, 3 + 1e-9, 0.05)
    zf = np.arange(0, 2.0 + 1e-9, 0.01)
    XX, YY = np.meshgrid(xg, yg, indexing='ij')
    tgt = np.c_[XX.ravel(), YY.ravel()]
    lev = np.abs(X[:, 2][:, None] - zl[None, :]).argmin(1)
    okz = (X[:, 2] >= -1e-6) & (np.abs(X[:, 2] - zl[lev]) < 2e-4) & (X[:, 0] > -3.5) & (X[:, 0] < 21)
    L = {k: np.full((len(zl), len(xg), len(yg)), np.nan) for k in ('S', 'I')}
    for iz in range(len(zl)):
        sel = okz & (lev == iz)
        P = X[sel][:, :2]
        for key, v in (('S', S[sel]), ('I', (S[sel] > s0).astype(float))):
            f = np.isfinite(S[sel])
            pp, vf = P[f], v[f]
            _, u = np.unique(np.round(pp, 6), axis=0, return_index=True)
            pp, vf = pp[u], vf[u]
            val = LinearNDInterpolator(Delaunay(pp), vf)(tgt)
            val[cKDTree(pp).query(tgt)[0] > 0.25] = np.nan
            L[key][iz] = val.reshape(len(xg), len(yg))
    out = {}
    i1 = np.clip(np.searchsorted(zl, zf), 1, len(zl) - 1)
    i0 = i1 - 1
    w = ((zf - zl[i0]) / (zl[i1] - zl[i0]))[:, None, None]
    for key in L:
        out[key] = (1 - w) * L[key][i0] + w * L[key][i1]
    with np.errstate(all='ignore'):
        Sm, Fm = np.nanmean(out['S'], axis=2), np.nanmean(out['I'], axis=2)
    np.savez(CACHE, x=xg, z=zf, Smean=Sm, frac=Fm, s0=s0)
    return xg, zf, Sm, Fm
