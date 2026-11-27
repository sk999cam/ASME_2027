"""Full-pitch conditional span means of the meff_span files on the 193M grid, and the near-wall split computed from them (helper module).

Source: /mnt/Elements/DNS/analysis/asme2027_figs/meff_span/span_RANS50.npz and span_C193b.npz (c_* arrays: plain mean over the points of
each (x, z) bin of the full pitch with the 193M |S*| D/U_inf > 0.05 mask applied point by point before the average; c_frac_valid = share
of the span points kept).
Grids: 193M 743 x bins (0.02 D, mesh x lines downstream of x/D 10.25) by 163 z levels; RANS 386 x bins (0.04 D) by 110 z levels.
The RANS conditional means are interpolated onto the 193M bin centres (x_D, z_D): bilinear in (z, x) over the RANS bin centres
(scipy RegularGridInterpolator), query points clamped to the RANS centre range (at most 0.003 D beyond the last RANS x centre), a NaN
in one of the four surrounding RANS bins gives NaN. The 193M arrays are used as they are. Valid fraction f and the mask are the 193M ones
for both runs (the file stores the 193M strain mask for the RANS points too).

Split (compute_curve): in a layer of bins and an x slab, with weights w = dx dz f (bin area times the valid fraction, 193M grid):
  total = ln( sum w mu_RANS / sum w mu_193M ), mu = c_mu_mean;
  four-player Shapley split of ln( sum w prod_RANS / sum w prod_193M ) over rho, k^(1/2), L_eff, cos, prod = rho k^(1/2) L_eff cos built
  from the c_<q>_mean of the bin (RANS cos = 1).
The mean of a product is not the product of the means, so the four terms sum to ln(prod ratio), not exactly to the total; the gap
(sum of the terms - total) is returned as 'res'.
"""
import itertools

import numpy as np
from scipy.interpolate import RegularGridInterpolator as RGI

SPAN_DIR = '/mnt/Elements/DNS/analysis/asme2027_figs/meff_span/'
KEYS = dict(mu='c_mu_mean', rho='c_rho_mean', k='c_sqrtk_mean', L='c_Leff_mean', c='c_cos_mean')


NKEY = dict(mu='c_mu_N', rho='c_rho_N', k='c_sqrtk_N', L='c_Leff_N', c='c_cos_N')


def _rebin_matrix(xe_src, xe_tgt):
    """W[t, s] = overlap length of source bin s with target bin t, divided by the source bin width (fraction of the source bin inside t)."""
    lo = np.maximum(xe_tgt[:-1, None], xe_src[None, :-1])
    hi = np.minimum(xe_tgt[1:, None], xe_src[None, 1:])
    return np.clip(hi - lo, 0, None) / np.diff(xe_src)[None, :]


def _rebin(f, dx, xe_tgt, fix=False):
    """sample-count-weighted rebinning in x of the conditional means of one run onto the bins xe_tgt: each source bin contributes with the
    fraction of its points that lie in the target bin (overlap fraction) times its own count of points (c_<q>_N); the valid fraction f is
    rebinned as sum(c_N_valid) / sum(finite-strain points), finite-strain points = c_N_valid / c_frac_valid (N where c_frac_valid = 0)."""
    W = _rebin_matrix(f['x_edges'], xe_tgt)
    out = {}
    xs = f['x_D']
    jb = int(np.abs(xs - 10.15).argmin()) if fix else None
    if fix:
        # 193M only: the merged bin at the block interface (x/D 10.15, width 0.04 D, 544 points per level instead of about 470, interface nodes
        # counted in both blocks) is a different sample; its means, counts and valid fraction are replaced by the linear interpolation in x
        # between the neighbouring bins before the rebinning
        assert abs(xs[jb] - 10.15) < 0.01 and np.diff(f['x_edges'])[jb] > 0.035
        t_ = (xs[jb] - xs[jb - 1]) / (xs[jb + 1] - xs[jb - 1])

        def fixc(a):
            a = a.astype(np.float64).copy()
            a[jb] = (1 - t_) * a[jb - 1] + t_ * a[jb + 1]
            return a
    else:
        def fixc(a):
            return a.astype(np.float64)
    for k, nm in KEYS.items():
        q, n = fixc(f[nm]), fixc(f[NKEY[k]])                                 # [nx, nz]
        n = np.where(np.isfinite(q), n, 0.0)
        num, den = W @ (np.where(n > 0, q, 0.0) * n), W @ n
        with np.errstate(all='ignore'):
            out[k] = np.where(den > 0, num / den, np.nan).T
    if 'c_frac_valid' in f.files:
        nv = fixc(f['c_N_valid'])
        fr = fixc(f['c_frac_valid'])
        with np.errstate(all='ignore'):
            nfin = np.where(fr > 0, nv / fr, fixc(f['N']))
        nfin = np.where(np.isfinite(nfin), nfin, 0.0)
        with np.errstate(all='ignore'):
            out['f'] = np.where((W @ nfin) > 0, (W @ np.where(np.isfinite(nv), nv, 0.0)) / (W @ nfin), np.nan).T
    return out


DX = float(__import__('os').environ.get('F14_DX', 0.1))     # x bin width of the rebinned means (D)


def load(span_dir=SPAN_DIR, dx=DX, case='RANS50'):
    """Both runs on the same x bins (width dx, centres at multiples of dx, x/D -2 to 20; dx = None keeps the native 193M x bins) and the 193M z
    levels. With dx the conditional means of each run are rebinned in x with sample-count weights before anything else (see _rebin); the RANS
    z levels are then interpolated linearly onto the 193M z levels (same x bins, so the interpolation is in z only)."""
    r, d = np.load(span_dir + f'span_{case}.npz'), np.load(span_dir + 'span_C193b.npz')       # case: 'RANS50' (default) or 'C20bW' (20M LES), the run compared with the 193M
    zc = d['z_D']
    if dx is None:
        xc, xe = d['x_D'], d['x_edges']
        q = np.stack(np.meshgrid(np.clip(zc, r['z_D'][0], r['z_D'][-1]), np.clip(xc, r['x_D'][0], r['x_D'][-1]), indexing='ij'), -1)
        R, D = {}, {}
        for k, nm in KEYS.items():
            D[k] = d[nm].T.astype(np.float64)                                  # [nz, nx]
            R[k] = RGI((r['z_D'], r['x_D']), r[nm].T.astype(np.float64), bounds_error=False, fill_value=None)(q)
        f = d['c_frac_valid'].T.astype(np.float64)
    else:
        xe = np.arange(-2.0 - dx / 2, 20.0 + dx, dx)
        xc = 0.5 * (xe[1:] + xe[:-1])
        Rb, Db = _rebin(r, dx, xe), _rebin(d, dx, xe, fix=True)
        f = Db.pop('f')
        D = Db
        R = {}
        for k in KEYS:
            ip = RGI((r['z_D'], xc), Rb[k], bounds_error=False, fill_value=None)
            R[k] = ip(np.stack(np.meshgrid(np.clip(zc, r['z_D'][0], r['z_D'][-1]), xc, indexing='ij'), -1))
        # RANS bins that were empty in x (NaN) stay NaN
    return dict(xc=xc, zc=zc, xe=xe, ze=d['z_edges'], R=R, D=D, f=f, mu_ref=float(d['mu_ref']), dx=dx,
                area=np.outer(np.diff(d['z_edges']), np.diff(xe)))


def _shapley4(fv):
    phi = np.zeros(4)
    ords = list(itertools.permutations(range(4)))
    for o in ords:
        st = [False] * 4
        prev = fv[tuple(st)]
        for i in o:
            st[i] = True
            cur = fv[tuple(st)]
            phi[i] += (cur - prev) / len(ords)
            prev = cur
    return phi


def split_cells(S, m):
    """split of the bins selected by the boolean mask m [nz, nx]."""
    w = (S['area'] * np.where(np.isfinite(S['f']), S['f'], 0.0))[m]
    R = {k: v[m] for k, v in S['R'].items()}
    D = {k: v[m] for k, v in S['D'].items()}
    keys = ['rho', 'k', 'L', 'c']
    fv = {}
    for s in itertools.product([False, True], repeat=4):
        prod = np.ones(m.sum())
        for i, kk in enumerate(keys):
            prod = prod * (D[kk] if s[i] else R[kk])
        v = (w * prod).sum()
        fv[s] = np.log(v) if v > 0 else np.nan
    # fv[s]: s[i] True = 193M value of factor i; the split is ln(RANS/193M): all False is RANS, all True is 193M
    tot = np.log((w * R['mu']).sum() / (w * D['mu']).sum())
    p4 = -_shapley4(fv)                      # Shapley of going from RANS to 193M, sign reversed: contribution to ln(RANS/193M)
    vm = {kk: ((w * R[kk]).sum() / w.sum(), (w * D[kk]).sum() / w.sum()) for kk in keys}
    vm['c'] = (1.0, vm['c'][1])
    chk = fv[(False,) * 4] - fv[(True,) * 4]
    return dict(tot=tot, p4=p4, res=p4.sum() - tot, chk=chk, vmean=vm, fneg=0.0, n=int(m.sum()), wsum=w.sum())


def compute_curve(S, zx_tab, centres, width=1.0, fmin=0.42, ztop=2.0, layer='near'):
    """split_cells in slabs of x (|x - centre| < width/2, bin centres) for the bins below z_x(x) ('near': z > 0, centre z < z_x) or at and above
    it ('above': up to ztop). Bins with f < fmin or a non-finite or non-positive mean are left out."""
    XC, ZC = np.meshgrid(S['xc'], S['zc'])
    zx = np.interp(XC, *zx_tab)
    ok = np.isfinite(S['f']) & (S['f'] >= fmin) & (ZC > 0) & (ZC <= ztop + 1e-9)
    for k in KEYS:
        ok &= np.isfinite(S['R'][k]) & np.isfinite(S['D'][k])
    ok &= (S['R']['mu'] > 0) & (S['D']['mu'] > 0) & (S['R']['rho'] > 0) & (S['D']['rho'] > 0) & (S['R']['k'] > 0) & (S['D']['k'] > 0)
    ok &= (S['R']['L'] > 0) & (S['D']['L'] > 0) & (S['D']['c'] > 0)
    lay = (ZC < zx) if layer == 'near' else (ZC >= zx)
    out = []
    for xc_ in np.atleast_1d(centres):
        slab = lay & (np.abs(XC - xc_) < width / 2) & (ZC > 0) & (ZC <= ztop + 1e-9) & np.isfinite(S['f'])
        m = slab & ok
        r = split_cells(S, m)
        r.update(x=float(xc_), kept=float((S['area'] * S['f'])[m].sum() / S['area'][slab].sum()))
        out.append(r)
    return out
