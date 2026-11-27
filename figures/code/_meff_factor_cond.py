"""Conditional span means of mu and of the factors of mu = rho k^(1/2) L cos(theta), both rungs, strain mask applied in 3D BEFORE the span mean.

Method (grid of span_fields.npz: x/D -2..20 step 0.05, y/D -3..3 step 0.05, z/D 0..2 step 0.01; point file pts.npz):
 1. every field (mu_t, mu_eff, rho, k^(1/2), L, cos, and the 193M strain |S*| D/U_inf) is interpolated onto the regular grid exactly as
    span_compute.py does: linear in (x, y) (Delaunay) on each common z level of the point file, then linear in z between levels; grid
    points further than 0.25 D from a data point stay empty. Points that are non-finite for a field are left out of that field's
    interpolation (fields with the same finite pattern share one triangulation).
 2. the validity mask is computed on the grid from the INTERPOLATED point-wise strain: valid = |S*| D/U_inf > S0 (threshold applied to the
    interpolated strain; no 0/1 flag is interpolated). The same mask is used for both rungs.
 3. at each (x, z): conditional span mean of each field over the valid points that are finite for that field; f = valid / finite-strain
    points of the span line (valid fraction); n = finite-strain points of the span line (complete span: >= 119 of 121).
Output (cached in SCR): x, z, f, n, nvalid, and cond[key] for key in muR, muD, rhoR, k12R, LR, rhoD, k12D, LD, cD.
Lengths: RANS L_R = mu_t/(rho k^(1/2)); 193M L_eff = L of the point file (mu_eff = rho k^(1/2) L_eff cos).
"""
import os

import numpy as np
from scipy.interpolate import LinearNDInterpolator
from scipy.spatial import Delaunay, cKDTree

import _meff_pointwise_data as M

KEYS = ['muR', 'muD', 'rhoR', 'k12R', 'LR', 'rhoD', 'k12D', 'LD', 'cD']


def get(s0=0.05):
    cache = M.SCR + 'span_cond_%03d.npz' % int(round(s0 * 1000))
    if os.path.exists(cache):
        c = np.load(cache)
        return c['x'], c['z'], c['f'], c['n'], c['nvalid'], {k: c['c_' + k] for k in KEYS}
    d = np.load(M.SCR + 'pts.npz')
    X = np.load(M.SCR + 'xyz.npz')['xyz']
    g = lambda n: d[n].astype(np.float64)  # noqa: E731
    with np.errstate(all='ignore'):
        F = dict(muR=g('rans_mut') / M.MUN, muD=g('c193_rho') * np.sqrt(g('c193_k')) * g('c193_L') * g('c193_cos') / M.MUN,
                 rhoR=g('rans_rho'), k12R=np.sqrt(g('rans_k')), LR=g('rans_mut') / (g('rans_rho') * np.sqrt(g('rans_k'))),
                 rhoD=g('c193_rho'), k12D=np.sqrt(g('c193_k')), LD=g('c193_L'), cD=g('c193_cos'),
                 S=g('c193_S') * M.D / M.UINF)
    names = list(F)
    zl = np.unique(np.round(X[d['blk'] == 7, 2], 5))
    xg = np.arange(-2, 20 + 1e-9, 0.05)
    yg = np.arange(-3, 3 + 1e-9, 0.05)
    zf = np.arange(0, 2.0 + 1e-9, 0.01)
    XX, YY = np.meshgrid(xg, yg, indexing='ij')
    tgt = np.c_[XX.ravel(), YY.ravel()]
    lev = np.abs(X[:, 2][:, None] - zl[None, :]).argmin(1)
    okz = (X[:, 2] >= -1e-6) & (np.abs(X[:, 2] - zl[lev]) < 2e-4) & (X[:, 0] > -3.5) & (X[:, 0] < 21)
    L = {k: np.full((len(zl), len(xg), len(yg)), np.nan) for k in names}
    for iz in range(len(zl)):
        sel = okz & (lev == iz)
        P = X[sel][:, :2]
        fin = {k: np.isfinite(F[k][sel]) for k in names}
        groups = {}
        for k in names:
            groups.setdefault(fin[k].tobytes(), []).append(k)
        for ks in groups.values():
            f = fin[ks[0]]
            pp = P[f]
            _, u = np.unique(np.round(pp, 6), axis=0, return_index=True)
            pp = pp[u]
            V = np.stack([F[k][sel][f][u] for k in ks], axis=1)
            val = LinearNDInterpolator(Delaunay(pp), V)(tgt)
            val[cKDTree(pp).query(tgt)[0] > 0.25] = np.nan
            for j, k in enumerate(ks):
                L[k][iz] = val[:, j].reshape(len(xg), len(yg))
    i1 = np.clip(np.searchsorted(zl, zf), 1, len(zl) - 1)
    i0 = i1 - 1
    w = ((zf - zl[i0]) / (zl[i1] - zl[i0]))[:, None, None]
    grid = lambda k: (1 - w) * L[k][i0] + w * L[k][i1]  # noqa: E731
    S3 = grid('S')                                       # [nz, nx, ny] interpolated strain
    fin_s = np.isfinite(S3)
    valid = fin_s & (S3 > s0)                            # threshold on the interpolated strain
    n = fin_s.sum(2)
    nvalid = valid.sum(2)
    with np.errstate(all='ignore'):
        f_ = nvalid / np.where(n > 0, n, np.nan)
    out = {}
    for k in KEYS:
        a = grid(k)
        m = valid & np.isfinite(a)
        cnt = m.sum(2)
        with np.errstate(all='ignore'):
            out[k] = np.where(cnt > 0, np.where(m, a, 0.0).sum(2) / cnt, np.nan)
    np.savez(cache, x=xg, z=zf, f=f_, n=n, nvalid=nvalid, **{'c_' + k: v for k, v in out.items()})
    return xg, zf, f_, n, nvalid, out
