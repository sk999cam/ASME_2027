"""F8b fig_mean_features: vortex motion and strength, RANS 50M (dashed) against 193M 3DNS (solid), figure x/D 3-20 (x/D 2-3 core values are
contaminated by wall-layer vorticity). The figure has two panels: (a) CRVP core |y_c|/D, z_c/D; (b) circulation of the CRVP and of the wall vortex.
Sources: CRVP core and circulation from the fine cross-planes at the stations x/D 3, 4, 5, 6, 8, 10, 12, 15, 20 for both rungs (one source, joined by
lines; the 0.05 D volume under-reads the RANS circulation at x/D 3-8); 193M band = between the two sides. Wall vortex from the fine planes (sign rule,
closed 0.5 contour, swirl check; open marker = peak < 3 x plane noise or no swirl). The rest of the analysis (volume tracks, upwash, plume, fluxes)
is kept in the table figures/code/data/mean_features_table.txt, which this script writes; the effect on flux and eta is in fig_flux_eta.py.
Older description of the analysis follows.

Mean flow features of the jet, RANS 50M against 193M 3DNS (b010), x/D 2-20.

Usage: python3 fig_mean_features.py build   (once: RANS cell export -> regular grid cache, long read of the 30 GB csv)
       python3 fig_mean_features.py         (analysis, table figures/code/data/mean_features_table.txt, figure figures/fig_mean_features.pdf)

Data
  193M: flow_features_c193b/n4_volume_c193b.npz (x, y, z grid, 0.05 D; U, theta, omx) and n3_crossplanes_c193b.npz (0.02 D planes at x/D 2, 4, 10,
        used only to check the swirl centre).
        theta = (T0_inf - T0)/(T0_inf - T0_cool) (total temperature, as in the file note); omx in 1/s.
  RANS 50M: cache built here on the same (x, y, z) grid as the 193M volume, asme2027_figs/mean_features/rans50_mean_grid.npz.
        U, V, W and theta from the reference cell export on Tank (cells_rans_p0dns_full.csv; theta from T0 = T + |U|^2/(2 cp), cp fitted
        from the export, theta = (T0_inf - T0)/(T0_inf - T0_cool)); inverse-distance interpolation of the 8 nearest cells (as the omega_x
        cache). omega_x D/U_inf is copied from overview_3d/rans50_q_omx_grid.npz (same grid, same interpolation, from the gradient export).
  Fluxes through z = 0.5 D: wall_balance/n7_summary.txt (193M mean Jm and turbulent Jt, RANS/193M ratio), stations x/D 2, 5, 10, 20.

Processing (identical for both rungs)
  omega_x D/U_inf, U/U_inf, theta are smoothed by a Gaussian, SIG = (sx, sy, sz) grid cells (0.05 D each), applied to all of x, y, z.
  Per half plane (y<0: main sign omega_x > 0; y>0: main sign < 0), |y/D| < YMAX_C, z/D >= ZMIN_C. Core = vorticity-weighted centroid (weights
  |omega_x| of the main sign) of the connected region around the |omega_x| peak above FRAC_C of the peak. Tracking in x from x/D 2: the peak is
  searched inside a window of +-TRK_W around the previous core.
  Circulation = area integral of the main-sign omega_x in the connected region above the level (0.2 of the local peak; 0.1 and 0.4 for
  sensitivity), U_inf D units.
  Secondary vortices by the sign rule: outboard of the main core (|y| larger), the first opposite-sign extremum (WV), then the next main-sign
  extremum further out (OV).
"""
import os
import sys

import numpy as np
from scipy import ndimage as ndi

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, '/home/sk999/Documents/3DNS_FP')

ROOT = '/home/sk999/Documents/ASME_2027/figures'
OUT = ROOT + '/fig_mean_features.pdf'
TABLE = ROOT + '/code/data/mean_features_table.txt'
FF = '/mnt/Elements/DNS/analysis/flow_features_c193b/'
VOL = FF + 'n4_volume_c193b.npz'
N3 = FF + 'n3_crossplanes_c193b.npz'
CACHE = '/mnt/Elements/DNS/analysis/asme2027_figs/mean_features'
RCACHE = CACHE + '/rans50_mean_grid.npz'
RCSV = '/mnt/Tank/STARCCM/film_cooling/rans_50M_Lm25b_BR15_FST_ZPG_p0dns/cells_rans_p0dns_full.csv'   # read only, fast disk
OMX_CACHE = '/mnt/Elements/DNS/analysis/asme2027_figs/overview_3d/rans50_q_omx_grid.npz'
D = 7.75e-3
UINF = 86.28
T0INF = 300.0
T0COOL = 199.02
CP = 1005.0
OKABE_ITO = ['#000000', '#E69F00', '#56B4E9', '#009E73', '#F0E442', '#0072B2', '#D55E00', '#CC79A7']

SIG = (2.0, 3.0, 3.0)
YMAX_C = 1.5
ZMIN_C = 0.05
FRAC_C = 0.5
TRK_W = 0.3
LEVELS = (0.2, 0.1, 0.4)
XS = np.arange(2.0, 20.0001, 0.25)
TAB_X = (2, 4, 6, 10, 15, 20)


def build_rans():
    import pandas as pd
    from scipy.spatial import cKDTree
    if os.path.exists(RCACHE):
        return np.load(RCACHE)
    os.makedirs(CACHE, exist_ok=True)
    V = np.load(VOL)
    gx, gy, gz = V['x'], V['y'], V['z']
    cols = ['Velocity[i] (m/s)', 'Velocity[j] (m/s)', 'Velocity[k] (m/s)', 'Temperature (K)', 'Total Enthalpy (J/kg)', 'X (m)', 'Y (m)', 'Z (m)']
    keep = []
    n = 0
    for ch in pd.read_csv(RCSV, usecols=cols, chunksize=2_000_000, dtype=np.float64):
        n += len(ch)
        x, y, z = ch['X (m)'].to_numpy() / D, ch['Y (m)'].to_numpy() / D, ch['Z (m)'].to_numpy() / D
        m = (x > -2.3) & (x < 20.3) & (np.abs(y) < 3.3) & (z > -0.15) & (z < 3.3)
        keep.append(ch.loc[m, cols].to_numpy())
        print('rows', n, 'kept', sum(len(k) for k in keep), flush=True)
    a = np.concatenate(keep)
    del keep
    P = a[:, 5:8] / D
    um2 = (a[:, :3] ** 2).sum(1)
    cp = np.polyfit(a[::50, 3], a[::50, 4] - 0.5 * um2[::50], 1)[0]       # cp from H - |U|^2/2 = cp T + const
    print('cp from the export', cp, flush=True)
    t0 = a[:, 3] + um2 / (2 * cp)
    om = np.load(OMX_CACHE)['omx_n']                                      # same IDW grid, from cell dW/dy - dV/dz (read only)
    q = np.c_[a[:, 0], a[:, 1], a[:, 2], (T0INF - t0) / (T0INF - T0COOL)]
    tr = cKDTree(P)
    X, Y, Z = np.meshgrid(gx, gy, gz, indexing='ij')
    G = np.c_[X.ravel(), Y.ravel(), Z.ravel()]
    d, i = tr.query(G, k=8, workers=16)
    w = 1.0 / np.maximum(d, 1e-6)
    w /= w.sum(1, keepdims=True)
    far = (d[:, 0] > 0.2).reshape(X.shape)
    out = {}
    for j, name in enumerate(('U', 'V', 'W', 'theta')):
        g = (q[i, j] * w).sum(1).reshape(X.shape)
        g[far] = np.nan
        out[name] = g.astype(np.float32)
    out['omx_n'] = om
    np.savez(RCACHE, x=gx, y=gy, z=gz, n_rows=n, cp=cp, note='RANS 50M base, IDW 8 cells onto the n4 grid; U, V, W m/s; omx_n = omega_x D/U_inf '
             '(copied from overview_3d/rans50_q_omx_grid.npz); theta = (T0_inf - T0)/(T0_inf - T0_cool), T0 = T + |U|^2/(2 cp)', **out)
    return np.load(RCACHE)


N3B = FF + 'n3b_crossplanes_c193b.npz'
PLANES_B = (3, 5, 6, 8, 12, 15, 20)
SLAB = 0.06                       # half thickness of the RANS slab in D (as fig_flow_features.py)
GRAD_CSV = '/mnt/Elements/DNS/analysis/rans50_exports/rans50_base_cells_grad.csv'
RPLB = CACHE + '/rans50_planes_b.npz'
SERIES = '/mnt/Elements/DNS/analysis/asme2027_figs/wall_balance/n7_series_c193b_b010_rans50.npz'
# flux through the inclined surface z = z_x(x); swap the file when the final z_x table is written (keys: <rung>_Jm_w, _Jt_w, _Jm_u, _Jm, 193_*_se)
SURF = '/mnt/Elements/DNS/analysis/asme2027_figs/wall_balance/n7_surface_c193b_b010_rans50_zxprov.npz'


def build_planes_b():
    """RANS fine planes at x/D PLANES_B on the (y, z) grid of the 193M planes, same method as fig_flow_features.py: cells within +-SLAB of the
    plane, linear (Delaunay) interpolation, omega_x = dW/dy - dV/dz normalised by U_inf/D. Cached in CACHE (read-only source csv)."""
    import pandas as pd
    from scipy.interpolate import LinearNDInterpolator
    if os.path.exists(RPLB):
        return np.load(RPLB)
    cols = ['Velocity[j] (m/s)', 'Velocity[k] (m/s)', 'dVdz', 'dWdy', 'X (m)', 'Y (m)', 'Z (m)']
    keep = {p: [] for p in PLANES_B}
    n = 0
    for ch in pd.read_csv(GRAD_CSV, usecols=cols, chunksize=2_000_000, dtype=np.float64):
        x = ch['X (m)'].to_numpy() / D
        for p in PLANES_B:
            m = np.abs(x - p) < SLAB
            if m.any():
                keep[p].append(ch.loc[m, cols].to_numpy())
        n += len(ch)
        print(n, flush=True)
    g = np.load(N3B)
    Y, Zg = np.meshgrid(g['y'], g['z'], indexing='ij')
    out = {'y': g['y'], 'z': g['z']}
    for p in PLANES_B:
        a = np.vstack(keep[p])
        v, w, dvdz, dwdy, yy, zz = a[:, 0], a[:, 1], a[:, 2], a[:, 3], a[:, 4 + 1] / D, a[:, 6] / D
        sel = (np.abs(yy) < 3.3) & (zz < 2.2)
        pts = np.column_stack([yy[sel], zz[sel]])
        for name, q in (('om', (dwdy - dvdz) * D / UINF), ('V', v), ('W', w)):
            out['%s_x%d' % (name, p)] = LinearNDInterpolator(pts, q[sel])(Y, Zg).astype(np.float32)
        print(p, a.shape, flush=True)
    np.savez(RPLB, **out)
    return np.load(RPLB)


if __name__ == '__main__' and len(sys.argv) > 1 and sys.argv[1] == 'build':
    build_rans()
    sys.exit(0)
if __name__ == '__main__' and len(sys.argv) > 1 and sys.argv[1] == 'planes':
    build_planes_b()
    sys.exit(0)

DY = 0.05
BOX_C = 0.6      # half size of the box that clips the contour for the core centroid
BOX_G = 0.8      # half size of the box that clips the contour for the circulation
SIG_S = (1.0, 2.0, 2.0)   # lighter smoothing for the secondary cores (small, near wall); CRVP and everything else use SIG
THR_S = {'rans': 0.06, 'dns': 0.06}     # detection threshold of secondary cores, omega_x D/U_inf; about the 193M noise floor at SIG_S (0.05), same for both
ZMAX_S = 1.0     # secondary cores searched below this height only (above it the cores found are noise)
X_INIT = 4.0
NAMES = ('rans', 'dns')


def load(name, sig=None):
    """Smoothed fields on the common grid: om (omega_x D/U), u, v, w (/U_inf), th (theta), x, y, z."""
    sig = SIG if sig is None else sig
    if name == 'dns':
        d = np.load(VOL)
        om = d['omx'] * D / UINF
        u, v, w = (d['U'][..., k] / UINF for k in range(3))
        th = d['theta']
    else:
        d = build_rans()
        om = d['omx_n']
        u, v, w = d['U'] / UINF, d['V'] / UINF, d['W'] / UINF
        th = d['theta']
    f = {'x': d['x'], 'y': d['y'], 'z': d['z']}
    for k, a in (('om', om), ('u', u), ('v', v), ('w', w), ('th', th)):
        f[k] = ndi.gaussian_filter(np.nan_to_num(np.asarray(a, dtype=np.float64)), sig, mode='nearest')
    return f


def box(y, z, yc, zc, h):
    return (np.abs(y[:, None] - yc) < h) & (np.abs(z[None, :] - zc) < h)


def region(f, y, z, pk, yp, zp, frac, h):
    """Connected region around (yp, zp) where f >= frac * pk, clipped to a box of half size h and z >= ZMIN_C."""
    m = (f >= frac * pk) & box(y, z, yp, zp, h) & (z[None, :] >= ZMIN_C)
    lab, _ = ndi.label(m)
    iy, iz = np.argmin(np.abs(y - yp)), np.argmin(np.abs(z - zp))
    if lab[iy, iz] == 0:
        return np.zeros_like(m)
    return lab == lab[iy, iz]


def find_core(f, y, z, yc0, zc0, win, side_mask):
    """Peak of f (main sign positive) in the window around (yc0, zc0); returns centroid (y, z), peak, peak position."""
    m = box(y, z, yc0, zc0, win) & side_mask[:, None] & (z[None, :] >= ZMIN_C) & (np.abs(y)[:, None] < YMAX_C)
    g = np.where(m, f, -np.inf)
    iy, iz = np.unravel_index(np.argmax(g), g.shape)
    pk = f[iy, iz]
    r = region(f, y, z, pk, y[iy], z[iz], FRAC_C, BOX_C)
    wgt = np.where(r, f, 0.0)
    yc = (wgt.sum(1) * y).sum() / wgt.sum()
    zc = (wgt.sum(0) * z).sum() / wgt.sum()
    return yc, zc, pk, (y[iy], z[iz])


def circ(f, y, z, pk, yc, zc, frac):
    r = region(f, y, z, pk, yc, zc, frac, BOX_G)
    return float(f[r].sum() * DY * DY)


def secondary(om, y, z, s, yc, thr):
    """Sign rule on one half plane. om: (ny, nz); s: main sign. Returns dict for WV, OV (peak, y, z, circulation)."""
    out = {}
    mx = lambda a: ndi.maximum_filter(a, size=(5, 5), mode='nearest')
    fm = s * om
    ok = (z[None, :] >= ZMIN_C) & (z[None, :] <= ZMAX_S)
    prev = np.abs(yc)
    for key, a in (('wv', -fm), ('ov', fm)):
        loc = (a == mx(a)) & (a > thr) & ok & ((np.sign(yc) * y)[:, None] > prev + 0.0) & (np.abs(y)[:, None] < 3.0)
        idx = np.argwhere(loc)
        if len(idx) == 0:
            out[key] = None
            break
        k = np.argmin(np.abs(y[idx[:, 0]]))
        iy, iz = idx[k]
        pk = a[iy, iz]
        out[key] = dict(pk=pk, y=y[iy], z=z[iz], g=circ(a, y, z, pk, y[iy], z[iz], 0.2))
        prev = np.abs(y[iy]) + 0.1
    out.setdefault('ov', None)
    return out


def analyse(name):
    F = load(name)
    F['oms'] = load(name, SIG_S)['om']
    x, y, z = F['x'], F['y'], F['z']
    ix = np.rint((XS - x[0]) / (x[1] - x[0])).astype(int)
    sides = {'L': (-1, y < 0, 1.0), 'R': (1, y > 0, -1.0)}      # y sign, mask, main sign of omega_x
    res = {k: {} for k in sides}
    i0 = int(np.argmin(np.abs(XS - X_INIT)))
    for sd, (sy, mask, s) in sides.items():
        n = len(XS)
        R = {q: np.full(n, np.nan) for q in ('yc', 'zc', 'pk', 'g20', 'g10', 'g40')}
        for sec in ('wv', 'ov'):
            for q in ('pk', 'y', 'z', 'g'):
                R[sec + '_' + q] = np.full(n, np.nan)
        prev = None
        for order in (range(i0, n), range(i0 - 1, -1, -1)):
            prev = None
            for k in order:
                om = F['om'][ix[k]]
                f = s * om
                if prev is None:
                    yc, zc, pk, _ = find_core(f, y, z, sy * 0.8, 0.5, 0.8, mask)
                else:
                    yc, zc, pk, _ = find_core(f, y, z, prev[0], prev[1], TRK_W, mask)
                prev = (yc, zc)
                R['yc'][k], R['zc'][k], R['pk'][k] = yc, zc, pk
                for lev, key in zip(LEVELS, ('g20', 'g10', 'g40')):
                    R[key][k] = circ(f, y, z, pk, yc, zc, lev)
                sc = secondary(F['oms'][ix[k]], y, z, s, yc, THR_S[name])
                for sec in ('wv', 'ov'):
                    if sc.get(sec):
                        for q in ('pk', 'y', 'z', 'g'):
                            R[sec + '_' + q][k] = sc[sec][q]
        res[sd] = R
    return F, ix, res


def interp_z(a, z, zq):
    return np.interp(zq, z, a)


def cross_to(F, ix, res):
    """Velocity, theta and noise quantities per x (arrays over XS)."""
    y, z = F['y'], F['z']
    n = len(ix)
    o = {q: np.full(n, np.nan) for q in ('ycm', 'zcm', 'Wc', 'Wmax', 'zWmax', 'Wflank', 'yflank', 'Vin', 'zcen', 'zthmax', 'th_wall',
                                          'halfw', 'thick', 'noise', 'floor', 'noise_s')}
    iy0 = int(np.argmin(np.abs(y)))
    for k in range(n):
        yc = 0.5 * (np.abs(res['L']['yc'][k]) + np.abs(res['R']['yc'][k]))
        zc = 0.5 * (res['L']['zc'][k] + res['R']['zc'][k])
        o['ycm'][k], o['zcm'][k] = yc, zc
        w, v, th = F['w'][ix[k]], F['v'][ix[k]], F['th'][ix[k]]
        o['Wc'][k] = interp_z(w[iy0], z, zc)
        j = int(np.argmax(w[iy0]))
        o['Wmax'][k], o['zWmax'][k] = w[iy0, j], z[j]
        wz = np.array([interp_z(w[i], z, zc) for i in range(len(y))])
        m = np.abs(y) < 2.0
        jm = np.argmin(np.where(m, wz, np.inf))
        o['Wflank'][k], o['yflank'][k] = wz[jm], y[jm]
        iz = int(np.argmin(np.abs(z - 0.1)))
        o['Vin'][k] = 0.5 * (v[(y < 0) & (np.abs(y) < YMAX_C), iz].max() + (-v[(y > 0) & (np.abs(y) < YMAX_C), iz]).max())
        t = np.clip(th, 0, None)
        dz = z[1] - z[0]
        o['zcen'][k] = (t * z[None, :]).sum() / t.sum()
        c = th[iy0]
        jt = int(np.argmax(c))
        o['zthmax'][k], o['th_wall'][k] = z[jt], c[0]
        below = np.where(c[jt:] <= 0.5 * c[jt])[0]
        if len(below):
            b = jt + below[0]
            o['thick'][k] = np.interp(0.5 * c[jt], [c[b], c[b - 1]], [z[b], z[b - 1]]) if b > jt else z[jt]
        gy = t.sum(1) * dz
        hw = []
        for sgn in (1, -1):
            yy, gg = sgn * y[iy0:] if sgn == 1 else -y[:iy0 + 1][::-1], (gy[iy0:] if sgn == 1 else gy[:iy0 + 1][::-1])
            below = np.where(gg <= 0.5 * gg[0])[0]
            if len(below):
                b = below[0]
                hw.append(np.interp(0.5 * gg[0], [gg[b], gg[b - 1]], [yy[b], yy[b - 1]]))
        o['halfw'][k] = np.mean(hw) if len(hw) == 2 else np.nan
        far = (np.abs(y)[:, None] > 2.0) & (z[None, :] > 1.5)
        o['noise'][k] = np.sqrt(np.mean(F['om'][ix[k]][far] ** 2))
        reg = (np.abs(y)[:, None] > 1.0) & (np.abs(y)[:, None] < 2.5) & (z[None, :] > 1.5) & (z[None, :] < 2.45)
        o['floor'][k] = np.abs(F['oms'][ix[k]][reg]).max()
        o['noise_s'][k] = np.sqrt(np.mean(F['oms'][ix[k]][far] ** 2))
    return o


def swirl_centre(F, xd, res, k):
    """Zero of the smoothed in-plane velocity (V, W) in a window of 0.5 D around the vorticity core, per side."""
    y, z = F['y'], F['z']
    i = int(np.argmin(np.abs(F['x'] - xd)))
    out = []
    for sd in 'LR':
        yc, zc = res[sd]['yc'][k], res[sd]['zc'][k]
        sp = np.hypot(F['v'][i], F['w'][i])
        m = box(y, z, yc, zc, 0.5) & (z[None, :] >= ZMIN_C)
        g = np.where(m, sp, np.inf)
        iy, iz = np.unravel_index(np.argmin(g), g.shape)
        out.append((y[iy], z[iz], np.hypot(y[iy] - yc, z[iz] - zc)))
    return out


def fluxes():
    """Mean (Jm) and turbulent (Jt) flux through z = 0.5 D, g s^-1 m^-1, 193M and RANS (from the RANS/193M ratios), stations x/D 2, 5, 10, 20."""
    import re
    txt = open('/mnt/Elements/DNS/analysis/asme2027_figs/wall_balance/n7_summary.txt').read()
    blk = txt.split('## h/D = 0.5')[1].split('## h/D')[0]
    rows = []
    for ln in blk.splitlines():
        t = re.findall(r'[+-]?\d+\.\d+', ln)
        if re.match(r'\s*\d+\.\d\s', ln) and len(t) >= 11:
            rows.append([float(a) for a in t[:11]])
    a = np.array(rows)            # x, Jm, SE, Jt, SE, Jt/Jm, SE, RJm, SE, RJt, SE
    return dict(x=a[:, 0], Jm=a[:, 1], Jt=a[:, 3], Jm_r=a[:, 1] * a[:, 7], Jt_r=a[:, 3] * a[:, 9], rJm=a[:, 7], rJt=a[:, 9])


def first_below(x, pk, thr):
    """First x after the maximum of pk where pk < thr (NaN counts as below); NaN if it stays above."""
    pk = np.where(np.isnan(pk), 0.0, pk)
    i0 = int(np.argmax(pk))
    for i in range(i0, len(x)):
        if pk[i] < thr:
            return x[i]
    return np.nan


def sm(R, q):
    return np.nanmean([R['L'][q], R['R'][q]], axis=0)


def quantities(R, o):
    """Table quantities per x (arrays over XS) for one rung."""
    q = {}
    q['|yc|/D'], q['zc/D'] = o['ycm'], o['zcm']
    q['peak CRVP'] = sm(R, 'pk')
    q['Gamma CRVP 0.2'], q['Gamma CRVP 0.1'], q['Gamma CRVP 0.4'] = sm(R, 'g20'), sm(R, 'g10'), sm(R, 'g40')
    q['W/U centre at zc'], q['W/U max centre'], q['z of W max'] = o['Wc'], o['Wmax'], o['zWmax']
    q['W/U flank min at zc'], q['V/U inward z=0.1'] = o['Wflank'], o['Vin']
    q['z centroid'], q['z theta max (y=0)'], q['theta wall (y=0)'] = o['zcen'], o['zthmax'], o['th_wall']
    q['half width'], q['thickness (y=0, 0.5 theta max)'] = o['halfw'], o['thick']
    q['noise |omx| (SIG)'], q['noise |omx| (SIG_S)'] = o['noise'], o['noise_s']
    q['noise-floor extreme (SIG_S)'] = o['floor']
    return q


SIG_P = 1.0       # smoothing of the fine planes, cells of 0.02 D, both rungs
PLANE_X = (2, 4, 10)
BOX_P = 0.6
ZMIN_P = 0.10     # the WV core lies above the wall-layer sheet (RANS x/D 4: z/D 0.18)
ZMAX_W = 0.6
DY_W = 0.3        # outboard of the main core by at least this; strongest opposite-sign maximum in that window
THR_P = 0.02


def wv_planes_old(out):
    """Superseded by wv_planes (kept for reference; not called)."""
    n3 = np.load(N3)
    rp = np.load('/mnt/Elements/DNS/analysis/asme2027_figs/flow_features/rans50_planes.npz')
    y, z = n3['y'], n3['z']
    h = y[1] - y[0]
    res = {}
    for nm in NAMES:
        res[nm] = {}
        for xd in PLANE_X:
            om = (n3['omx_x%d' % xd] * D / UINF) if nm == 'dns' else rp['om_x%d' % xd]
            om = ndi.gaussian_filter(np.nan_to_num(np.asarray(om, dtype=np.float64)), SIG_P, mode='nearest')
            far = (np.abs(y)[:, None] > 2.0) & (z[None, :] > 1.5)
            r = {'noise': float(np.sqrt(np.mean(om[far] ** 2)))}
            k = int(np.argmin(np.abs(XS - xd)))
            for sd, s in (('L', 1.0), ('R', -1.0)):
                yc = out[nm]['R'][sd]['yc'][k]
                a = -s * om
                loc = (a == ndi.maximum_filter(a, size=(9, 9), mode='nearest')) & (a > THR_P) & (z[None, :] >= ZMIN_P) & (z[None, :] <= ZMAX_W) \
                    & ((np.sign(yc) * y)[:, None] > abs(yc) + DY_W) & (np.abs(y)[:, None] < 3.0)
                idx = np.argwhere(loc)
                if len(idx) == 0:
                    r[sd] = dict(pk=np.nan, y=np.nan, z=np.nan, g=np.nan)
                    continue
                iy, iz = idx[np.argmax(a[idx[:, 0], idx[:, 1]])]
                pk = a[iy, iz]
                m = (a >= 0.2 * pk) & (np.abs(y[:, None] - y[iy]) < BOX_P) & (np.abs(z[None, :] - z[iz]) < BOX_P) & (z[None, :] >= ZMIN_P)
                lab, _ = ndi.label(m)
                g = a[lab == lab[iy, iz]].sum() * h * h
                r[sd] = dict(pk=float(pk), y=float(y[iy]), z=float(z[iz]), g=float(g))
            res[nm][xd] = r
    return res


STATIONS = (2, 3, 4, 5, 6, 8, 10, 12, 15, 20)
ZWALL = 0.06      # a secondary core needs a closed 0.5-of-peak contour that stays above this height
SIG_C = 7.5       # fine-plane CRVP smoothing in cells of 0.02 D (= 0.15 D, the y, z smoothing of the volume analysis)
DW_S = 0.2        # half distance of the swirl check


def plane_fields(nm, xd):
    """Raw omega_x D/U_inf, V/U_inf, W/U_inf on the (y, z) grid of the 193M fine planes."""
    if nm == 'dns':
        d = np.load(N3 if xd in PLANE_X else N3B)
        U = d['U_x%d' % xd]
        return d['y'], d['z'], d['omx_x%d' % xd] * D / UINF, U[..., 1] / UINF, U[..., 2] / UINF
    if xd in PLANE_X:
        r = np.load('/mnt/Elements/DNS/analysis/asme2027_figs/flow_features/rans50_planes.npz')
        g = np.load(N3)
    else:
        r = g = build_planes_b()
    return g['y'], g['z'], r['om_x%d' % xd], r['V_x%d' % xd] / UINF, r['W_x%d' % xd] / UINF


def closed_ok(a, y, z, iy, iz, frac):
    """Connected region a >= frac * peak around (iy, iz); closed = does not touch z < ZWALL nor the array edges."""
    m = a >= frac * a[iy, iz]
    lab, _ = ndi.label(m)
    r = lab == lab[iy, iz]
    zs = np.where(r.any(0))[0]
    ys = np.where(r.any(1))[0]
    ok = z[zs.min()] >= ZWALL and zs.max() < len(z) - 1 and ys.min() > 0 and ys.max() < len(y) - 1
    return ok


def swirl_ok(V, W, y, z, y0, z0):
    """Sign change of W across the core in y and of V across the core in z (in-plane swirl); returns (w_change, v_change)."""
    iy, iz = int(np.argmin(np.abs(y - y0))), int(np.argmin(np.abs(z - z0)))
    n = int(round(DW_S / (y[1] - y[0])))
    iz0 = max(iz - n, 0)
    return bool(W[iy - n, iz] * W[iy + n, iz] < 0), bool(V[iy, iz0] * V[iy, min(iz + n, len(z) - 1)] < 0)


def pick(a, V, W, y, z, yc, noise, same_sign_after=None, ymax_off=1.6):
    """Strongest closed maximum of a (positive = the sign sought) outboard of the main core |yc| (WV), or the nearest outboard of
    same_sign_after (OV). Returns dict with pk, y, z, g, swirl flags, closed, below3."""
    thr = max(THR_P, 1.5 * noise)
    h = y[1] - y[0]
    lo = abs(yc) + DY_W if same_sign_after is None else same_sign_after + 0.2
    hi = abs(yc) + ymax_off if same_sign_after is None else 3.0
    loc = (a == ndi.maximum_filter(a, size=(9, 9), mode='nearest')) & (a > thr) & (z[None, :] >= ZWALL) & (z[None, :] <= ZMAX_W) \
        & ((np.sign(yc) * y)[:, None] > lo) & ((np.sign(yc) * y)[:, None] < hi)
    cands = [(a[i, j], i, j) for i, j in np.argwhere(loc) if closed_ok(a, y, z, i, j, 0.5)]
    nan = dict(pk=np.nan, y=np.nan, z=np.nan, g=np.nan, sw=False, below3=True)
    if not cands:
        return nan
    cands.sort(reverse=same_sign_after is None) if same_sign_after is None else cands.sort(key=lambda c: abs(y[c[1]]))
    pk, iy, iz = cands[0]
    m = (a >= 0.2 * pk) & (np.abs(y[:, None] - y[iy]) < BOX_P) & (np.abs(z[None, :] - z[iz]) < BOX_P) & (z[None, :] >= ZWALL)
    lab, _ = ndi.label(m)
    g = a[lab == lab[iy, iz]].sum() * h * h
    sw = swirl_ok(V, W, y, z, y[iy], z[iz])
    return dict(pk=float(pk), y=float(y[iy]), z=float(z[iz]), g=float(g), sw=bool(sw[0]), sw_w=sw[0], sw_v=sw[1], below3=bool(pk < 3 * noise))


def crvp_fine(om, y, z, yv, zv, s):
    """Main core from the fine plane inside a window of +-TRK_W around the volume core (yv, zv); same centroid / 0.2 circulation definitions."""
    f = s * ndi.gaussian_filter(np.nan_to_num(om), SIG_C, mode='nearest')
    h = y[1] - y[0]
    side = (y < 0) if s > 0 else (y > 0)
    yc, zc, pk, (yp, zp) = find_core(f, y, z, yv, zv, TRK_W, side)
    r = region(f, y, z, pk, yc, zc, 0.2, BOX_G)
    return dict(y=yc, z=zc, pk=pk, g=float(f[r].sum() * h * h))


def wv_planes(out):
    """Per rung, station, side: WV, OV (fine planes, SIG_P smoothing) and the CRVP core (SIG_C smoothing)."""
    res = {}
    for nm in NAMES:
        res[nm] = {}
        for xd in STATIONS:
            y, z, om, V, W = plane_fields(nm, xd)
            oms = ndi.gaussian_filter(np.nan_to_num(np.asarray(om, dtype=np.float64)), SIG_P, mode='nearest')
            Vs = ndi.gaussian_filter(np.nan_to_num(np.asarray(V, dtype=np.float64)), 3, mode='nearest')
            Ws = ndi.gaussian_filter(np.nan_to_num(np.asarray(W, dtype=np.float64)), 3, mode='nearest')
            far = (np.abs(y)[:, None] > 2.0) & (z[None, :] > 1.5)
            r = {'noise': float(np.sqrt(np.mean(oms[far] ** 2)))}
            k = int(np.argmin(np.abs(XS - xd)))
            for sd, s in (('L', 1.0), ('R', -1.0)):
                yv, zv = out[nm]['R'][sd]['yc'][k], out[nm]['R'][sd]['zc'][k]
                a = -s * oms
                wv = pick(a, Vs, Ws, y, z, yv, r['noise'])
                ov = pick(s * oms, Vs, Ws, y, z, yv, r['noise'], same_sign_after=abs(wv['y']) if not np.isnan(wv['y']) else abs(yv) + 1.0)
                r[sd] = dict(wv, ov=ov, crvp=crvp_fine(om, y, z, yv, zv, s))
            res[nm][xd] = r
    return res


def run():
    out = {}
    for nm in NAMES:
        F, ix, R = analyse(nm)
        o = cross_to(F, ix, R)
        out[nm] = dict(F=F, ix=ix, R=R, o=o, q=quantities(R, o))
    out['wvp'] = wv_planes(out)
    return out


def write_table(out, fl):
    q = {nm: out[nm]['q'] for nm in NAMES}
    lines = ['# F8b mean-flow features, RANS 50M against 193M 3DNS (b010). Table of fig_mean_features.py.',
             '# Smoothing: Gaussian, SIG = %s grid cells (0.05 D) in (x, y, z), both rungs. Cores: centroid inside the %.1f contour of the peak, |y|<%.1f, z>=%.2f.' % (SIG, FRAC_C, YMAX_C, ZMIN_C),
             '# Values are means of the two sides (y<0, y>0) except where stated. Lengths in D, W, V in U_inf, omega_x in U_inf/D, Gamma in U_inf D.',
             '# Gamma = circulation inside the connected region above 0.2 (0.1, 0.4) of the local peak, clipped to a box of +-%.1f D around the core.' % BOX_G,
             '# Wall vortex: see the section below the ratios (fine cross-planes at x/D 2, 4, 10); the volume WV/OV rows were removed (features thinner than the grid).',
             '# x/D 2-3 core values are contaminated by wall-layer vorticity; the figure starts at x/D 3.',
             '# Ratio = RANS/193M.', '']
    hdr = '%-34s' % 'quantity' + ''.join('x/D %-5d' % x + '      ' for x in TAB_X)
    for title, sel in (('RANS 50M', 'rans'), ('193M 3DNS', 'dns')):
        lines += [title, hdr]
        for k in q['rans']:
            lines.append('%-34s' % k + ''.join('%11.4f ' % q[sel][k][np.argmin(np.abs(XS - x))] for x in TAB_X))
        lines.append('')
    lines += ['Ratio RANS/193M', hdr]
    for k in q['rans']:
        r = [q['rans'][k][np.argmin(np.abs(XS - x))] / q['dns'][k][np.argmin(np.abs(XS - x))] if q['dns'][k][np.argmin(np.abs(XS - x))] != 0 else np.nan
             for x in TAB_X]
        lines.append('%-34s' % k + ''.join('%11.3f ' % v for v in r))
    lines += ['']
    lines += ['Wall vortex (WV) and outer vortex (OV) from the fine cross-planes (0.02 D; 193M n3 + n3b, RANS rans50_planes.npz + mean_features/rans50_planes_b.npz), sign rule,',
              'Gaussian sigma = %.0f cell (%.2f D) in y, z for both rungs; replaces the smoothed-volume WV/OV rows. WV = strongest opposite-sign maximum %.1f-%.1f D outboard of the main core,'
              % (SIG_P, 0.02 * SIG_P, DY_W, 1.6),
              'z/D %.2f-%.1f, with a closed 0.5-of-peak contour that stays above z/D %.2f; OV = nearest main-sign maximum outboard of the WV (same conditions). Gamma inside the 0.2-of-peak' % (ZWALL, ZMAX_W, ZWALL),
              'contour (box +-%.1f D). sw = swirl check passed (W changes sign across the core in y, +-%.1f D; V across it in z is not required, a wall-side core has no data below it; in-plane velocity smoothed 3 cells). low = peak below 3 x plane noise' % (BOX_P, DW_S),
              '(rms of |omega_x| D/U_inf at z>1.5, |y|>2). Peak in omega_x D/U_inf, positions in D, Gamma in U_inf D; nan = no valid pick.']
    for nm in NAMES:
        for xd in STATIONS:
            nz = out['wvp'][nm][xd]['noise']
            for sd in 'LR':
                for key in ('wv', 'ov'):
                    w = out['wvp'][nm][xd][sd] if key == 'wv' else out['wvp'][nm][xd][sd]['ov']
                    lines.append('  %-5s x/D %2d %s %s: peak %.3f  |y| %.2f  z %.2f  Gamma %.4f  sw %-5s low %-5s (noise %.4f)' % (
                        nm, xd, sd, key.upper(), w['pk'], abs(w['y']), w['z'], w['g'], w['sw'], w['below3'], nz))
    lines += ['', 'CRVP core and circulation from the fine planes (sigma %.1f cells = %.2f D, same definitions) against the volume analysis (mean of the two sides)' % (SIG_C, 0.02 * SIG_C),
              '  rung  x/D   |yc| fine  vol    zc fine  vol    Gamma fine  vol']
    for nm in NAMES:
        for xd in STATIONS:
            k = int(np.argmin(np.abs(XS - xd)))
            cf = [out['wvp'][nm][xd][sd]['crvp'] for sd in 'LR']
            lines.append('  %-5s %3d   %.3f %.3f   %.3f %.3f   %.4f %.4f' % (
                nm, xd, np.mean([abs(c['y']) for c in cf]), out[nm]['q']['|yc|/D'][k], np.mean([c['z'] for c in cf]), out[nm]['q']['zc/D'][k],
                np.mean([c['g'] for c in cf]), out[nm]['q']['Gamma CRVP 0.2'][k]))
    lines += ['', 'Mean-field noise (rms of the smoothed |omega_x| D/U_inf at z/D>1.5, |y/D|>2), mean over x/D 2-20, volume smoothing SIG: RANS %.4f, 193M %.4f' % (
        np.nanmean(out['rans']['o']['noise']), np.nanmean(out['dns']['o']['noise']))]
    S = np.load(SERIES)
    xs = S['x_D']
    lines += ['', 'Coolant flux through z = h, g s^-1 m^-1, one pitch, + away from the wall (N7 series, wall_balance/n7_series_c193b_b010_rans50.npz; slab x +-0.25 D; SE = jackknife, 10 batches)',
              'Ratio RANS/193M only where the 193M value exceeds 2 SE, otherwise "193M not different from zero". Turbulent share = Jt/(Jm+Jt).',
              '%-5s %-5s %9s %7s %9s %7s %9s %9s %-30s %-30s %7s %7s' % ('h/D', 'x/D', 'Jm 193M', 'SE', 'Jt 193M', 'SE', 'Jm RANS', 'Jt RANS', 'Jm RANS/193M', 'Jt RANS/193M', 'sh193', 'shRANS')]
    for h in (0.25, 0.5, 1.0):
        for xd in (4, 6, 10, 15, 20):
            i = int(np.argmin(np.abs(xs - xd)))
            jm, jt, sm_, st = S['h%s_193_Jm' % h][i], S['h%s_193_Jt' % h][i], S['h%s_193_Jm_se' % h][i], S['h%s_193_Jt_se' % h][i]
            rm, rt = S['h%s_rans_Jm' % h][i], S['h%s_rans_Jt' % h][i]
            qm = '%.3f' % (rm / jm) if abs(jm) > 2 * sm_ else '193M not different from zero'
            qt = '%.3f' % (rt / jt) if abs(jt) > 2 * st else '193M not different from zero'
            lines.append('%-5s %-5d %9.3f %7.3f %9.3f %7.3f %9.3f %9.3f %-30s %-30s %7.3f %7.3f' % (h, xd, jm, sm_, jt, st, rm, rt, qm, qt, jt / (jm + jt), rt / (rm + rt)))
    lines += ['', 'Stations with 193M closure residual > 5 % (x/D >= 3): ' + '; '.join(
        'h %s: %s' % (h, ', '.join('%.2f' % x for x in xs[(S['h%s_193_pct' % h] > 5.0) & (xs >= 3)])) for h in (0.25, 0.5, 1.0))]
    Z = np.load(SURF)
    xz = Z['x_D']
    lines += ['', 'Flux through the inclined surface z = z_x(x) (%s; z_x/D %.2f-%.2f; surface table %s), g s^-1 m^-1, one pitch, + away from the wall' % (
        os.path.basename(SURF), Z['zx_D'].min(), Z['zx_D'].max(), ', '.join('%g:%.2f' % tuple(r) for r in Z['surface_table'])),
        'Jm_w, Jt_w = vertical parts (F_z); Jm_u = streamwise-flux slope term (-F_x dz_x/dx) of the mean flux; Jm = full mean flux (Jm_w + Jm_u). The figure shows the vertical parts only.',
        '%-5s %-5s %9s %9s %9s %9s %9s %9s %9s %9s' % ('x/D', 'zx', 'Jm_w193', 'SE', 'Jt_w193', 'SE', 'Jm_u193', 'Jm193', 'Jm_wRANS', 'Jt_wRANS') + '   %9s %9s %-28s %-28s' % (
            'Jm_uRANS', 'JmRANS', 'Jm_w RANS/193M', 'Jt_w RANS/193M')]
    for xd in (4, 6, 10, 15, 20):
        i = int(np.argmin(np.abs(xz - xd)))
        g = lambda k: Z[k][i]
        qm = '%.3f' % (g('rans_Jm_w') / g('193_Jm_w')) if abs(g('193_Jm_w')) > 2 * g('193_Jm_w_se') else '193M not different from zero'
        qt = '%.3f' % (g('rans_Jt_w') / g('193_Jt_w')) if abs(g('193_Jt_w')) > 2 * g('193_Jt_w_se') else '193M not different from zero'
        lines.append('%-5d %-5.2f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f   %9.3f %9.3f %-28s %-28s' % (
            xd, g('zx_D'), g('193_Jm_w'), g('193_Jm_w_se'), g('193_Jt_w'), g('193_Jt_w_se'), g('193_Jm_u'), g('193_Jm'), g('rans_Jm_w'), g('rans_Jt_w'),
            g('rans_Jm_u'), g('rans_Jm'), qm, qt))
    for cs, nm in (('193', '193M'), ('rans', 'RANS')):
        j = np.where(Z[cs + '_Jt_w'] > Z[cs + '_Jm_w'])[0]
        lines.append('%s: turbulent vertical flux first exceeds the mean upwash (Jm_w) at x/D %s' % (nm, '%.2f' % xz[j[0]] if len(j) else 'never (x/D 3-20)'))
    os.makedirs(os.path.dirname(TABLE), exist_ok=True)
    open(TABLE, 'w').write('\n'.join(lines) + '\n')


def checks(out):
    """Swirl-centre check at x/D 4 and 10 (zero of the smoothed in-plane velocity within 0.5 D of the core); returns text lines."""
    L = ['', 'Check of the vorticity cores against the swirl centre (zero of in-plane |(V, W)| within 0.5 D of the core), sides L (y<0), R (y>0): core (y, z); swirl (y, z, distance)']
    for nm in NAMES:
        for xd in (4, 10):
            k = int(np.argmin(np.abs(XS - xd)))
            sw = swirl_centre(out[nm]['F'], xd, out[nm]['R'], k)
            L.append('  %-5s x/D %2d core %s swirl %s' % (nm, xd, [(round(out[nm]['R'][s]['yc'][k], 2), round(out[nm]['R'][s]['zc'][k], 2)) for s in 'LR'],
                                                        [tuple(round(float(a), 2) for a in t) for t in sw]))
    return L


def plot(out, fl):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from cycler import cycler
    from matplotlib.lines import Line2D
    from functions import set_mpl_params
    # Drawn at the printed size (0.7 text width = 5.25 in, include scale 1); body text 10 pt -> sf = 10/22.
    set_mpl_params(10 / 22, grid_on=True)
    plt.rcParams['axes.prop_cycle'] = cycler(color=OKABE_ITO)
    C = OKABE_ITO
    LS = {'dns': '-', 'rans': '--'}
    fig, ax = plt.subplots(1, 2, figsize=(5.25, 2.9))
    A, B = ax
    W = out['wvp']
    ST = np.array([s for s in STATIONS if s >= 3])

    def note(a, txt, xy, xyt, col, ha='left'):
        a.annotate(txt, xy=xy, xytext=xyt, color=col, ha=ha, va='center', arrowprops=dict(arrowstyle='->', color=col, lw=0.7))

    def fine(nm, key):
        """Fine-plane CRVP value per side at the stations ST, shape (2, n)."""
        return np.array([[W[nm][xd][sd]['crvp'][key] for xd in ST] for sd in 'LR'])

    # (a) core trajectory from the fine planes (one source for both rungs, ten stations joined); 193M: mean line with a band between the sides
    for nm in NAMES:
        for key, col, fn in (('y', C[5], np.abs), ('z', C[6], lambda v: v)):
            v = fn(fine(nm, key))
            A.plot(ST, v.mean(0), LS[nm], color=col, marker='o' if nm == 'dns' else '^', ms=2.5)
            if nm == 'dns':
                A.fill_between(ST, v.min(0), v.max(0), color=col, alpha=0.2, lw=0)
    A.set_ylabel(r'core position / $D$')
    A.set_ylim(0, 1.0)
    note(A, r'$|y_c|/D$', (15.0, abs(fine('dns', 'y')[:, list(ST).index(15)]).mean()), (14.0, 0.3), C[5], ha='center')
    note(A, r'$z_c/D$', (8.0, fine('dns', 'z')[:, list(ST).index(8)].mean()), (6.0, 0.9), C[6], ha='center')
    A.legend(handles=[Line2D([], [], color='k', ls='-', label='193M 3DNS'), Line2D([], [], color='k', ls='--', label='RANS 50M')],
             loc='lower right', fontsize=plt.rcParams['font.size'], handlelength=1.6)

    # (b) circulation: CRVP (fine planes, mean of the two sides) and wall vortex (fine planes)
    for nm in NAMES:
        B.plot(ST, fine(nm, 'g').mean(0), LS[nm], color=C[0], marker='o' if nm == 'dns' else '^', ms=2.5)
    # wall vortex markers from the fine planes: squares 193M (open if peak < 3 x plane noise), triangles RANS
    for nm, mk in (('dns', 's'), ('rans', '^')):
        xs, gm, lo, hi, op = [], [], [], [], []
        for xd in ST:
            w = out['wvp'][nm][xd]
            gs = [w[sd]['g'] for sd in 'LR' if not np.isnan(w[sd]['g'])]
            if not gs:
                continue
            xs.append(xd)
            gm.append(np.mean(gs))
            lo.append(min(gs))
            hi.append(max(gs))
            op.append(any(w[sd]['below3'] or not w[sd]['sw'] for sd in 'LR' if not np.isnan(w[sd]['g'])))
        xs, gm, lo, hi = map(np.array, (xs, gm, lo, hi))
        B.plot(xs, gm, LS[nm], color=C[3], lw=0.6)
        for xv, gv, o_ in zip(xs, gm, op):
            B.plot([xv], [gv], mk, color=C[3], mfc='none' if o_ else C[3], ms=5 if nm == 'rans' else 3, zorder=2 if nm == 'rans' else 3)
    B.set_ylabel(r'$\Gamma/(U_\infty D)$')
    B.set_ylim(0, 0.2)
    note(B, 'CRVP', (15.0, fine('dns', 'g')[:, list(ST).index(15)].mean()), (16.0, 0.115), C[0], ha='center')
    note(B, 'wall vortex', (15.0, np.mean([W['dns'][15][sd]['g'] for sd in 'LR']) + 0.002), (19.5, 0.036), C[3], ha='right')
    note(B, 'no RANS WV\nbeyond $x/D=4$', (4.3, W['rans'][4]['L']['g'] + 0.002), (8.0, 0.042), C[3], ha='center')
    B.text(0.97, 0.97, 'open: peak < 3 x noise\nor no swirl', transform=B.transAxes, ha='right', va='top', color=C[3])

    for a in (A, B):
        a.set_xlim(3, 20)
        a.set_xticks([5, 10, 15, 20])
    for a, t in zip((A, B), 'ab'):
        a.text(0.0, 1.02, '(%s)' % t, transform=a.transAxes, ha='left', va='bottom')
    for a in (A, B):
        a.set_xlabel(r'$x/D$')
    fig.tight_layout()
    fig.savefig(OUT)
    return fig


if __name__ == '__main__':
    import warnings
    warnings.filterwarnings('ignore')
    res = run()
    fl = fluxes()
    write_table(res, fl)
    chk = checks(res)
    print('\n'.join(chk))
    with open(TABLE, 'a') as fh:
        fh.write('\n'.join(chk) + '\n')
    plot(res, fl)
