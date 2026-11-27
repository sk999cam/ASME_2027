"""Overview of the mean flow in 3D, RANS 50M and 193M 3DNS (paper Fig. F4).

Usage: python3 fig_overview_3d.py            (writes /home/sk999/Documents/ASME_2027/figures/fig_overview_3d.png, 300 dpi)
       python3 fig_overview_3d.py --build    (only builds the RANS cache, see below)

Two panels, same camera, same isosurface level and cut: (a) RANS 50M, (b) 193M 3DNS (stat3 batch b010).
Shown: plate wall coloured by eta (jet, 0-0.8; white inside the hole exit); part-transparent isosurface of the normalised mean-flow
Q criterion, Q (D/U_inf)^2 = QLEV, coloured by the mean streamwise vorticity omega_x D/U_inf (bwr, +-OMAX), cells below z/D 0.1 removed
and cut at x/D XCUT; cross-planes x/D 4 and 10 carrying omega_x D/U_inf (bwr, same limits, opacity rising with |omega_x|) with in-plane
streamlines of the mean (V, W); empty frame at x/D 20 (no vorticity plane exists there); frame at the centre plane y = 0; hole exit
outline on z = 0. The hole walls lie under the plate and are not drawn. Domain x/D -2..20, y/D -3..3, z/D 0..3, true aspect.

Isosurface colour, panel (b) (stopgap, flagged in the panel): the 193M n4 volume has no omega_x (fields Q, l2, theta, u only), and
interpolation from the three n3 planes is not used. The 193M surface is therefore coloured by the side of the plume, omega_x D/U_inf
= -sign(y) 0.8 OMAX, the sign taken from the peaks of the 193M x/D 4 plane (omega_x < 0 for y > 0, > 0 for y < 0; same in RANS).
Panel (a) uses the RANS omega_x from the gradient export.

Data
  193M: flow_features_c193b/n4_volume_c193b.npz (Q_n), n4_walls_c193b.npz (plate eta), n3_crossplanes_c193b.npz (omx_x*, U_x*; V, W).
  RANS 50M: cell velocity gradients of rans50_base_cells_grad.csv (A_ij = dU_i/dx_j): Q = 0.5(|W|^2 - |S|^2) and omega_x = dW/dy - dV/dz,
    normalised as the 193M file ((U_inf/D)^2 and U_inf/D; U_inf 86.28 m/s, D 7.75 mm), inverse-distance interpolation (8 nearest cells) onto
    the 193M x-y-z grid, cached in /mnt/Elements/DNS/analysis/asme2027_figs/overview_3d/rans50_q_omx_grid.npz.
    Planes: asme2027_figs/flow_features/rans50_planes.npz (om_x*, V_x*, W_x*; made by fig_flow_features.py). Plate eta: wall_plate.csv.
Renderer: pyvista (off-screen) for the 3D scene; colour bars and labels by matplotlib (set_mpl_params).
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, '/home/sk999/Documents/3DNS_FP')

OUT = '/home/sk999/Documents/ASME_2027/figures/fig_overview_3d.png'
CACHE = '/mnt/Elements/DNS/analysis/asme2027_figs/overview_3d'
RCSV = '/mnt/Elements/DNS/analysis/rans50_exports/rans50_base_cells_grad.csv'
FF = '/mnt/Elements/DNS/analysis/flow_features_c193b/'
VOL = FF + 'n4_volume_c193b.npz'
WAL = FF + 'n4_walls_c193b.npz'
N3 = FF + 'n3_crossplanes_c193b.npz'
RPL = '/mnt/Elements/DNS/analysis/asme2027_figs/flow_features/rans50_planes.npz'
HOLE_CSV = '/home/sk999/Documents/4dns_filmcooled/starccm/rans_50M_Lm25b_BR15_FST_ZPG_p0dns/post_run1/wall_hole.csv'
SF = 0.709
SF_WIDE = 10.0 / 22.0                           # font 10 pt at include scale 1 (as fig_spots.py)
LAYOUT = os.environ.get('OVERVIEW_LAYOUT', 'two')    # 'two' (default): 193M above RANS, full text width; 'wide': one 7.5 in view, 193M only; 'tall': old one-column figure (include scale 0.64)
SF_TWO = 0.336                                       # font 7.4 pt in the 5.6 in design; printed at 1.35x = 10 pt (set from the measured width)
D = 7.75e-3
UINF = 86.28
NORM = (UINF / D) ** 2
ZCUT = 0.1
VMAX = 0.8
ZT = 3.0
HS = (-1.05, 0.0, 0.1)                          # horseshoe vortex anchor
CVP = [(10.0, -0.7, 0.6), (10.0, 0.6, 0.6)]     # vortex pair anchors (cores on the x/D 10 plane)

# final choice (see candidates in the scratchpad f4_candidates)
CFG = dict(cam=[(-10, -15, 11), (8, 0, 0.3), (0, 0, 1)], zoom=1.0, qlev=0.005, xcut=20.0, omax=0.15, yclip=None, minfrac=0.03, hole_ps=6,
           zcut=0.05, light=0.6, sil=3.0, color_by='omx', iso_lim=0.15, clim=(0.3, 3.0), cb_label=r'$|\bar{\omega}| D/U_\infty$',   # color_by: 'omx' | 'wmag' | 'umag' | 'z'
           sil_angle=45, shadow=0.4, shadow_sigma=3.0, iso_grey=True,
           op_iso=0.8, ontop_planes=(), ov0=0.5, ov1=0.3, sm_q=1.0, a0=0.25, a1=0.4, a_max=0.72, lab_dz=0.3, zp=2.0, zt=2.0, op_wall=1.0, planes=(2, 4, 10), win=(2400, 1300), dpi=300, out=OUT,
           sm=3.0, sm_uv=3.0, sm_iso=1.0,                     # smoothing (grid cells) of plotted 193M omega_x (plane, iso colour) and (V, W)
           core_thr=0.08, core_sm=2.5, rings=(0.12, 0.25), nring=6, gy=0.5, gz=0.4, cell=0.22, smin=0.4, lw=1.4, arrow=0.1,
           s1=(-1.6, 0.2), inset_cell=0.07,
           hcb=0.95, clear=0.1, cap_dy=0.9, cap_xmax=1.4, cap_clear=0.22, ovl=0.3, inset_pad=0.34, inset_fs=7.0, inset_w=1.4, xs_in=(0.62, 2.65, 4.1), hr=0.25, cell2=2.2, lead2='near', tag_xy=(0.01, 0.99), tag_ha='left', lead_dx=0.1, lead_dy=0.35, zmax_in=2.0, x_first=0.34, clr=0.12, lab_below=2.0, row=0.2, row0=0.06, up10=0.2, gap=0.05, tick_sp=0.78, lab_pad=0.3, lab_clip=2.6, row_dz=0.0, label_thr=0.095, img_w=5.2, w_img=5.35, in_w=1.95, in_right=0.04, in_ys=(2.95, 1.75, 0.55), wv_z=1.4, h_fig=3.9, cb_box=(3.75, 0.22, 1.5, 0.1, 0.5), dpi_wide=300, ss_wide=2, ss_two=3, dpi_two=810, top=1, w_two=6.5,
           pad_top=40, pad_bot=30, pad_hs=100, cvp_dx=-60, cvp_dy=-250)


def build_rans():
    """Read the CSV once in chunks, keep the 193M box, compute Q, lambda_2, omega_x, interpolate onto the 193M grid."""
    import pandas as pd
    from scipy.spatial import cKDTree
    f = os.path.join(CACHE, 'rans50_q_omx_grid.npz')
    if os.path.exists(f):
        return np.load(f)
    os.makedirs(CACHE, exist_ok=True)
    V = np.load(VOL)
    gx, gy, gz = V['x'], V['y'], V['z']
    cols = ['dUdx', 'dUdy', 'dUdz', 'dVdx', 'dVdy', 'dVdz', 'dWdx', 'dWdy', 'dWdz', 'X (m)', 'Y (m)', 'Z (m)']
    keep = []
    n = 0
    for ch in pd.read_csv(RCSV, usecols=cols, chunksize=2_000_000):
        n += len(ch)
        x, y, z = ch['X (m)'].to_numpy() / D, ch['Y (m)'].to_numpy() / D, ch['Z (m)'].to_numpy() / D
        m = (x > -2.3) & (x < 20.3) & (np.abs(y) < 3.3) & (z > -0.15) & (z < 3.3)
        keep.append(ch.loc[m, cols].to_numpy())
        print('rows', n, 'kept', sum(len(k) for k in keep), flush=True)
    a = np.concatenate(keep)
    del keep
    A = a[:, :9].reshape(-1, 3, 3)                      # A[i, j] = dU_i/dx_j
    P = a[:, 9:] / D
    S = 0.5 * (A + A.transpose(0, 2, 1))
    W = 0.5 * (A - A.transpose(0, 2, 1))
    Q = 0.5 * ((W ** 2).sum((1, 2)) - (S ** 2).sum((1, 2)))
    l2 = np.linalg.eigvalsh(S @ S + W @ W)[:, 1]
    om = (A[:, 2, 1] - A[:, 1, 2]) * D / UINF           # omega_x D / U_inf = (dW/dy - dV/dz) D / U_inf
    tr = cKDTree(P)
    X, Y, Z = np.meshgrid(gx, gy, gz, indexing='ij')
    G = np.c_[X.ravel(), Y.ravel(), Z.ravel()]
    d, i = tr.query(G, k=8, workers=8)
    w = 1.0 / np.maximum(d, 1e-6)
    w /= w.sum(1, keepdims=True)
    shp = X.shape
    Qg = (Q[i] * w).sum(1).reshape(shp) / NORM
    Lg = (l2[i] * w).sum(1).reshape(shp) / NORM
    Og = (om[i] * w).sum(1).reshape(shp)
    far = (d[:, 0] > 0.2).reshape(shp)                  # no cell within 0.2 D: outside the mesh
    for q in (Qg, Lg, Og):
        q[far] = np.nan
    np.savez(f, x=gx, y=gy, z=gz, Q_n=Qg.astype(np.float32), l2_n=Lg.astype(np.float32), omx_n=Og.astype(np.float32),
             n_rows=n, note='RANS 50M base; Q, lambda2 in (U_inf/D)^2, omega_x in U_inf/D; IDW 8 cells onto the n4 grid')
    return np.load(f)


def build_rans_extra():
    """RANS 50M |omega| D/U_inf and |U|/U_inf on the 193M grid (same cells and interpolation as build_rans); cached."""
    import pandas as pd
    from scipy.spatial import cKDTree
    f = os.path.join(CACHE, 'rans50_wmag_umag_grid.npz')
    if os.path.exists(f):
        return np.load(f)
    os.makedirs(CACHE, exist_ok=True)
    V = np.load(VOL)
    gx, gy, gz = V['x'], V['y'], V['z']
    cols = ['Velocity[i] (m/s)', 'Velocity[j] (m/s)', 'Velocity[k] (m/s)', 'dUdx', 'dUdy', 'dUdz', 'dVdx', 'dVdy', 'dVdz', 'dWdx', 'dWdy',
            'dWdz', 'X (m)', 'Y (m)', 'Z (m)']
    keep = []
    for ch in pd.read_csv(RCSV, usecols=cols, chunksize=2_000_000):
        x, y, z = ch['X (m)'].to_numpy() / D, ch['Y (m)'].to_numpy() / D, ch['Z (m)'].to_numpy() / D
        m = (x > -2.3) & (x < 20.3) & (np.abs(y) < 3.3) & (z > -0.15) & (z < 3.3)
        keep.append(ch.loc[m, cols].to_numpy())
    a = np.concatenate(keep)
    del keep
    A = a[:, 3:12].reshape(-1, 3, 3)                    # A[i, j] = dU_i/dx_j
    P = a[:, 12:] / D
    wx = A[:, 2, 1] - A[:, 1, 2]
    wy = A[:, 0, 2] - A[:, 2, 0]
    wz = A[:, 1, 0] - A[:, 0, 1]
    wm = np.sqrt(wx ** 2 + wy ** 2 + wz ** 2) * D / UINF
    um = np.sqrt((a[:, :3] ** 2).sum(1)) / UINF
    tr = cKDTree(P)
    X, Y, Z = np.meshgrid(gx, gy, gz, indexing='ij')
    G = np.c_[X.ravel(), Y.ravel(), Z.ravel()]
    d, i = tr.query(G, k=8, workers=8)
    w = 1.0 / np.maximum(d, 1e-6)
    w /= w.sum(1, keepdims=True)
    far = (d[:, 0] > 0.2).reshape(X.shape)
    out = {}
    for name, q in (('wmag_n', wm), ('umag_n', um)):
        g = (q[i] * w).sum(1).reshape(X.shape)
        g[far] = np.nan
        out[name] = g.astype(np.float32)
    np.savez(f, x=gx, y=gy, z=gz, note='RANS 50M base; |omega| D/U_inf and |U|/U_inf; IDW 8 cells onto the n4 grid', **out)
    return np.load(f)


def vort_193(V):
    """193M |omega| D/U_inf and |U|/U_inf from the stored mean velocity (n4 volume, 0.05 D grid, central differences, NaN = solid -> NaN)."""
    U = np.array(V['U'], float) / UINF
    dx = float(V['x'][1] - V['x'][0])
    g = [np.gradient(U[..., k], dx, axis=ax) for k in range(3) for ax in range(3)]    # g[3k+ax] = dU_k/dx_ax (x, y, z in D)
    wx, wy, wz = g[3 * 2 + 1] - g[3 * 1 + 2], g[3 * 0 + 2] - g[3 * 2 + 0], g[3 * 1 + 0] - g[3 * 0 + 1]
    return np.sqrt(wx ** 2 + wy ** 2 + wz ** 2), np.sqrt((U ** 2).sum(-1)), wx


def project(p, pts):
    """World points (D) -> pixel coordinates of the off-screen image, origin top-left."""
    import vtk
    c = vtk.vtkCoordinate()
    c.SetCoordinateSystemToWorld()
    h = p.window_size[1]
    out = []
    for q in pts:
        c.SetValue(*q)
        u, v = c.GetComputedDisplayValue(p.renderer)
        out.append((u, h - v))
    return np.array(out)


_DATA = {}


def data():
    """Everything that does not depend on the camera or opacities (loaded once per process)."""
    if _DATA:
        return _DATA
    from matplotlib.path import Path
    from scipy.interpolate import RegularGridInterpolator
    import fig_eta_x as F
    import fig_eta_wall as FW
    from hole_exit_outline import load as load_outline
    V = np.load(VOL)
    R = build_rans()
    Wd = np.load(WAL)
    n3 = np.load(N3)
    rp = np.load(RPL)
    hx, hy = load_outline()
    XG, YG = FW.XG, FW.YG
    Xw, Yw = np.meshgrid(XG, YG)
    ex, ey, ee = F.rans_points('RANS 50M')
    eta_r = np.ma.filled(FW.grid_eta(np.asarray(ex), np.asarray(ey), np.asarray(ee)), np.nan)
    eta_d = RegularGridInterpolator((Wd['y'], Wd['x']), np.array(Wd['eta']).T, bounds_error=False)((Yw, Xw))
    inh = Path(np.c_[hx, hy]).contains_points(np.c_[Xw.ravel(), Yw.ravel()]).reshape(Xw.shape)
    pl = {}
    for xp in (2, 4, 10):
        pl[xp] = dict(
            R=(rp['om_x%d' % xp], rp['V_x%d' % xp], rp['W_x%d' % xp], n3['y'], n3['z']),
            D=(n3['omx_x%d' % xp] * D / UINF, n3['U_x%d' % xp][..., 1], n3['U_x%d' % xp][..., 2], n3['y'], n3['z']))
    # x/D 20: omega_x only (the RANS cache holds omega_x on the n4 grid but not V, W), no streamlines
    i20 = int(np.argmin(np.abs(V['x'] - 20.0)))
    pl[20] = dict(R=(np.array(R['omx_n'][i20]), None, None, V['y'], V['z']),
                  D=(np.array(V['omx'][i20]) * D / UINF, None, None, V['y'], V['z']))
    import pandas as pd
    hw = pd.read_csv(HOLE_CSV, usecols=['X (m)', 'Y (m)', 'Z (m)']).to_numpy() / D
    hw = hw[(hw[:, 0] >= -2.0) & (hw[:, 2] > -1.5)]                         # inside the plotted domain
    RX = build_rans_extra()
    w193, u193, wx193 = vort_193(V)
    chk = np.isfinite(wx193) & np.isfinite(V['omx']) & (np.abs(V['omx']) > 0)
    print('check 193M omega_x from U vs stored omx: corr %.3f' % np.corrcoef(wx193[chk] * UINF / D, np.array(V['omx'])[chk])[0, 1])
    _DATA.update(wm=(np.array(RX['wmag_n']), w193), um=(np.array(RX['umag_n']), u193), V=V, R=R, hx=hx, hy=hy, Xw=Xw, Yw=Yw, eta=(eta_r, eta_d), inh=inh, planes=pl, hole=hw)
    return _DATA


def find_cores(o, y, z, thr, sigma, zp):
    """Local extrema of omega_x (smoothed) with |omega_x| > thr, z/D >= 0.08, |y/D| <= 2.8; returns [(y, z, omega_x)] sorted by |omega_x|."""
    from scipy.ndimage import gaussian_filter, maximum_filter, minimum_filter
    s = gaussian_filter(np.nan_to_num(o), sigma) if sigma > 0 else np.nan_to_num(o)
    ok = (z[None, :] >= 0.08) & (z[None, :] <= zp) & (np.abs(y[:, None]) <= 2.8)
    out = []
    for f in (maximum_filter, minimum_filter):
        m = (f(s, size=21) == s) & (np.abs(s) > thr) & ok
        out += [(y[i], z[j], s[i, j]) for i, j in zip(*np.where(m))]
    return sorted(out, key=lambda c: -abs(c[2]))


def trace_lines(y, z, v, w, seeds, cell, zp, smin, ds=0.03, lmax=4.0, lmin=0.35, ylim=3.0):
    """Evenly spaced streamlines of the in-plane velocity (v, w), unit-speed RK2 in both directions from the seeds (in priority order).
    A line stops when it enters a cell (size `cell`, in D) already used by an earlier line, leaves |y| <= ylim, 0.02 <= z <= zp,
    or where the in-plane speed is below smin (m/s). Returns a list of (n, 2) paths in (y, z), ordered along the flow."""
    from scipy.interpolate import RegularGridInterpolator
    fv = RegularGridInterpolator((y, z), np.nan_to_num(v), bounds_error=False, fill_value=0.0)
    fw = RegularGridInterpolator((y, z), np.nan_to_num(w), bounds_error=False, fill_value=0.0)

    def vel(p):
        a = np.array([fv(p)[0], fw(p)[0]])
        return a, np.hypot(*a)

    def inside(p):
        return abs(p[0]) <= ylim and 0.02 <= p[1] <= zp

    def cid(p):
        return int(np.floor((p[0] + 3) / cell)), int(np.floor(p[1] / cell))

    occ = set()
    paths = []
    for s in seeds:
        s = np.array(s, float)
        if not inside(s):
            continue
        if cid(s) in occ:
            continue
        own = {cid(s)}
        halves = []
        for sgn in (1.0, -1.0):
            p = s.copy()
            pts = []
            for _ in range(int(lmax / ds)):
                a, sp = vel(p)
                if sp < smin:
                    break
                pm = p + 0.5 * ds * sgn * a / sp
                a2, sp2 = vel(pm)
                if sp2 < smin:
                    break
                p = p + ds * sgn * a2 / sp2
                if not inside(p):
                    break
                c = cid(p)
                if c in occ:
                    break
                own.add(c)
                pts.append(p.copy())
            halves.append(pts)
        path = np.array(halves[1][::-1] + [s] + halves[0])
        if len(path) * ds < lmin:
            continue
        occ |= own
        paths.append(path)
    return paths


def plane_lines(paths, xp, arrow=0.09):
    """pyvista objects for streamline paths on the plane x = xp: one polyline set and one set of flat arrowheads at mid-length."""
    import pyvista as pv
    pts, cells = [], []
    tri, tf = [], []
    for p in paths:
        n0 = len(pts)
        pts += [(xp, a, b) for a, b in p]
        cells += [[len(p)] + list(range(n0, n0 + len(p)))]
        k = len(p) // 2
        t = p[min(k + 1, len(p) - 1)] - p[max(k - 1, 0)]
        t = t / (np.hypot(*t) + 1e-12)
        nrm = np.array([-t[1], t[0]])
        tip = p[k] + 0.5 * arrow * t
        b1 = p[k] - 0.5 * arrow * t + 0.35 * arrow * nrm
        b2 = p[k] - 0.5 * arrow * t - 0.35 * arrow * nrm
        base = len(tri)
        tri += [(xp, *tip), (xp, *b1), (xp, *b2)]
        tf += [[3, base, base + 1, base + 2]]
    lines = pv.PolyData(np.array(pts, float), lines=np.hstack(cells))
    heads = pv.PolyData(np.array(tri, float), faces=np.hstack(tf))
    return lines, heads


def scene(which, cfg):
    """Render one panel (which = 0 RANS, 1 193M); returns (image, projected anchors)."""
    import matplotlib
    import pyvista as pv
    from matplotlib.colors import Normalize
    d = data()
    V = d['V']
    z = V['z']
    omax = cfg['omax']
    Qn = np.array(d['R']['Q_n'] if which == 0 else V['Q_n'])
    e = np.array(d['eta'][which], dtype=float)
    e[d['inh']] = np.nan
    p = pv.Plotter(off_screen=True, window_size=cfg['win'])
    p.set_background('white')
    p.enable_depth_peeling(number_of_peels=12, occlusion_ratio=0.0)
    Xw, Yw = d['Xw'], d['Yw']
    pl = pv.StructuredGrid(Xw, Yw, np.zeros_like(Xw))
    # wall colour: jet(eta), darkened softly under the isosurface footprint (vertical drop shadow) so that the surface separates from the wall
    from scipy.ndimage import gaussian_filter as _gfw
    from scipy.interpolate import RegularGridInterpolator as _RGI
    qq = np.nan_to_num(Qn, nan=0.0)
    if which == 1 and cfg['sm_q'] > 0:
        qq = _gfw(qq, cfg['sm_q'])
    qq = np.where(((z >= cfg['zcut']) & True)[None, None, :] & (V['x'] <= cfg['xcut'])[:, None, None], qq, 0.0)
    fp = _gfw((qq > cfg['qlev']).any(2).astype(float), cfg['shadow_sigma'])
    shade = _RGI((V['x'], V['y']), fp, bounds_error=False, fill_value=0.0)(np.c_[Xw.ravel(), Yw.ravel()]).reshape(Xw.shape)
    wc = matplotlib.colormaps['jet'](Normalize(0, VMAX)(np.nan_to_num(e)))
    wc[..., :3] *= (1.0 - cfg['shadow'] * np.clip(shade, 0, 1))[..., None]
    wc[..., 3] = np.where(np.isfinite(e), cfg['op_wall'], 0.0)
    pl.point_data['rgba'] = (np.stack([wc[..., k].ravel(order='F') for k in range(4)], 1) * 255).astype(np.uint8)
    p.add_mesh(pl, scalars='rgba', rgba=True, lighting=False)
    if cfg.get('holewall', True):                       # hole wall of the RANS CAD (below the plate, seen through the exit)
        from matplotlib.path import Path
        hp = d['hole']
        cp = np.array(cfg['cam'][0], float)
        t = cp[2] / (cp[2] - hp[:, 2])                  # ray camera -> point hits the plate z = 0 at (hit_x, hit_y)
        hit = cp[None, :2] + t[:, None] * (hp[:, :2] - cp[None, :2])
        vis = Path(np.c_[d['hx'], d['hy']]).contains_points(hit)    # visible only through the hole exit
        gv = np.clip(0.78 + 0.22 * hp[vis][:, 2] / 1.2, 0.35, 0.78)           # grey shade darkening with depth
        p.add_points(hp[vis], scalars=(np.c_[gv, gv, gv] * 255).astype(np.uint8), rgb=True, point_size=cfg.get('hole_ps', 4) * cfg.get('ss', 1),
                     render_points_as_spheres=False, lighting=False)
    q = np.nan_to_num(Qn, nan=0.0)
    if which == 1 and cfg['sm_q'] > 0:                  # light smoothing of the 193M Q field before contouring
        from scipy.ndimage import gaussian_filter as _gf
        q = _gf(q, cfg['sm_q'])
    q[:, :, z < cfg['zcut']] = 0.0
    q[V['x'] > cfg['xcut'], :, :] = 0.0
    g = pv.ImageData(dimensions=q.shape, spacing=(0.05,) * 3, origin=(V['x'][0], V['y'][0], V['z'][0]))
    g.point_data['Q'] = q.ravel(order='F')
    if which == 0:
        om = np.nan_to_num(np.array(d['R']['omx_n']), nan=0.0)
    else:
        om = np.nan_to_num(np.array(V['omx']), nan=0.0) * D / UINF
        if cfg['sm_iso'] > 0:                           # light smoothing of the colour field only (193M)
            from scipy.ndimage import gaussian_filter
            om = gaussian_filter(om, cfg['sm_iso'])
    g.point_data['om'] = om.ravel(order='F')
    cb = cfg['color_by']
    if cb != 'omx':                                     # alternative isosurface colouring quantity (one-sided, inferno)
        if cb == 'wmag':
            sc = d['wm'][which]
        elif cb == 'umag':
            sc = d['um'][which]
        elif cb == 'sign':                              # sign of omega_x (193M: of the omega_x smoothed with sigma sm_iso)
            sc = np.sign(om)
        else:
            sc = np.broadcast_to(z[None, None, :], q.shape)
        sc = np.nan_to_num(np.array(sc, float), nan=0.0)
        if which == 1 and cfg['sm_iso'] > 0 and cb not in ('z', 'sign'):
            from scipy.ndimage import gaussian_filter
            sc = gaussian_filter(sc, cfg['sm_iso'])
        g.point_data['sc'] = sc.ravel(order='F')
    iso = g.contour([cfg['qlev']], scalars='Q')
    if cfg.get('minfrac', 0) > 0:                       # drop detached noise fragments (< minfrac of the largest piece, by cell count)
        cc = iso.connectivity('all')
        rid = np.asarray(cc.cell_data['RegionId'])
        cnt = np.bincount(rid)
        keep = cnt >= cfg['minfrac'] * cnt.max()
        iso = cc.extract_cells(np.where(keep[rid])[0]).extract_surface()
    print('isosurface', ('RANS', '193M')[which], 'Q =', cfg['qlev'], 'x/D range %.2f to %.2f' % iso.bounds[:2],
          '| y/D range %.2f %.2f' % iso.bounds[2:4])
    if cfg.get('yclip') is not None:                  # cut-away: only the far half (y >= yclip) of the isosurface is drawn
        iso = iso.clip(normal=(0, -1, 0), origin=(0, cfg['yclip'], 0), invert=True)
    import pyvista as _pv
    p.add_light(_pv.Light(position=(-4, -14, 22), focal_point=(8, 0, 0.3), color='white', intensity=cfg['light'], light_type='scene light'))
    if cfg.get('iso_grey'):                             # alternative: uniform light grey, sign carried by the planes only
        p.add_mesh(iso, color='#8c8c8c', smooth_shading=True, ambient=0.25, diffuse=0.75, specular=0.4, specular_power=30,
                   opacity=cfg['op_iso'], show_scalar_bar=False)
    elif cb == 'sign':
        p.add_mesh(iso, scalars='sc', cmap='bwr', clim=(-1, 1), n_colors=2, smooth_shading=True, ambient=0.3, diffuse=0.8, specular=0.5,
                   specular_power=30, opacity=cfg['op_iso'], show_scalar_bar=False)
    elif cb == 'omx':
        p.add_mesh(iso, scalars='om', cmap='bwr', clim=(-cfg['iso_lim'], cfg['iso_lim']), smooth_shading=True, ambient=0.3, diffuse=0.8, specular=0.5,
                   specular_power=30, opacity=cfg['op_iso'], show_scalar_bar=False)
    else:
        v_ = np.asarray(iso.point_data['sc'])
        print('iso colour', cb, ('RANS', '193M')[which], 'percentiles 2/50/98: %.3f %.3f %.3f' % tuple(np.percentile(v_, [2, 50, 98])))
        p.add_mesh(iso, scalars='sc', cmap='inferno', clim=cfg['clim'], smooth_shading=True, ambient=0.3, diffuse=0.8, specular=0.5,
                   specular_power=30, opacity=cfg['op_iso'], show_scalar_bar=False)
    if cfg['sil'] > 0:                                  # black silhouette and sharp-edge lines on the isosurface
        p.add_silhouette(iso, color='black', line_width=cfg['sil'] * cfg.get('ss', 1), feature_angle=cfg['sil_angle'])
    p.add_mesh(pv.lines_from_points(np.c_[d['hx'], d['hy'], np.full_like(d['hx'], 0.005)]), color='k', line_width=4 * cfg.get('ss', 1))
    zt = cfg['zt']
    fr = [(-2, 0, 0), (20, 0, 0), (20, 0, zt), (-2, 0, zt), (-2, 0, 0)]
    p.add_mesh(pv.lines_from_points(np.array(fr, float)), color='k', line_width=3 * cfg.get('ss', 1))
    # cross-planes: omega_x image (alpha rising with |omega_x|) and streamlines of (V, W)
    from scipy.ndimage import gaussian_filter
    cmap = matplotlib.colormaps['bwr']
    cores = {}
    top = []                                            # actors drawn in an overlay layer (plane colour and streamlines)
    for xp in (2, 4, 10, 20):
        if xp in cfg['planes']:
            o, v, w, y1, z1 = d['planes'][xp]['R' if which == 0 else 'D']
            zi = z1 <= cfg['zp'] + 1e-9                       # plane image up to z/D zp
            z1, o = z1[zi], np.array(o[:, zi], float)
            if which == 1 and cfg['sm'] > 0:                  # light smoothing of the plotted 193M vorticity (and velocity for lines)
                m = np.isnan(o)
                o = gaussian_filter(np.nan_to_num(o), cfg['sm'])
                o[m] = np.nan
            Yp, Zp = np.meshgrid(y1, z1, indexing='ij')
            on = np.nan_to_num(o)
            rgba = cmap(Normalize(-omax, omax)(on))
            a = np.clip((np.abs(on) - cfg['a0'] * omax) / (cfg['a1'] * omax), 0, 1) * cfg['a_max']   # alpha rises with |omega_x|
            rgba[..., 3] = np.where(np.isfinite(o), a, 0.0)
            sg = pv.StructuredGrid(np.full_like(Yp, xp), Yp, Zp)
            # StructuredGrid point order is x fastest: build the colour array in Fortran order
            sg.point_data['rgba'] = (np.stack([rgba[..., k].ravel(order='F') for k in range(4)], 1) * 255).astype(np.uint8)
            p.add_mesh(sg, scalars='rgba', rgba=True, lighting=False)      # true depth
            if xp in cfg['ontop_planes']:                  # overlay: only the high-vorticity part is drawn over the isosurface
                ov = rgba.copy()
                ov[..., 3] = np.where(np.isfinite(o), np.clip((np.abs(on) - cfg['ov0'] * omax) / (cfg['ov1'] * omax), 0, 1) * cfg['a_max'], 0.0)
                sg2 = pv.StructuredGrid(np.full_like(Yp, xp + 0.01), Yp, Zp)
                sg2.point_data['rgba'] = (np.stack([ov[..., k].ravel(order='F') for k in range(4)], 1) * 255).astype(np.uint8)
                top.append(p.add_mesh(sg2, scalars='rgba', rgba=True, lighting=False))
            if v is not None:
                v, w = np.array(v[:, zi], float), np.array(w[:, zi], float)
                if which == 1 and cfg['sm_uv'] > 0:
                    v, w = gaussian_filter(np.nan_to_num(v), cfg['sm_uv']), gaussian_filter(np.nan_to_num(w), cfg['sm_uv'])
                cs = find_cores(o, y1, z1, cfg['core_thr'], cfg['core_sm'] if which == 1 else 0, cfg['zp'])
                cores[xp] = cs
                seeds = []
                for yc, zc, _ in cs:                          # rings around every vortex core, then a coarse grid
                    for r in cfg['rings']:
                        for t in np.linspace(0, 2 * np.pi, cfg['nring'], endpoint=False):
                            seeds.append((yc + r * np.cos(t), zc + r * np.sin(t)))
                seeds += [(a_, b_) for a_ in np.arange(-2.9, 3.0, cfg['gy']) for b_ in np.arange(0.15, cfg['zp'], cfg['gz'])]
                paths = trace_lines(y1, z1, v, w, seeds, cfg['cell'], cfg['zp'], cfg['smin'])
                lines, heads = plane_lines(paths, xp - 0.03, cfg['arrow'])
                tl = top if xp in cfg['ontop_planes'] else []
                tl.append(p.add_mesh(lines, color='#202020', line_width=cfg['lw'] * cfg.get('ss', 1), opacity=0.95, lighting=False))
                if cfg['arrow'] > 0:
                    tl.append(p.add_mesh(heads, color='#202020', lighting=False))
        fr = [(xp, -3, 0), (xp, 3, 0), (xp, 3, zt), (xp, -3, zt), (xp, -3, 0)]
        p.add_mesh(pv.lines_from_points(np.array(fr, float)), color='#303030', line_width=3 * cfg.get('ss', 1))
    p.camera_position = cfg['cam']
    p.camera.zoom(cfg['zoom'])
    if cfg.get('ontop', True):                          # planes drawn over the isosurface wherever their alpha is > 0
        import vtk
        r2 = vtk.vtkRenderer()
        r2.SetLayer(1)
        r2.SetActiveCamera(p.camera)
        p.render_window.SetNumberOfLayers(2)
        p.render_window.AddRenderer(r2)
        for a in top:
            p.renderer.RemoveActor(a)
            r2.AddActor(a)
    p.render()
    img = p.screenshot(return_img=True)
    s1 = cfg['s1']
    anch = project(p, [HS] + CVP + [(xp, 3, zt) for xp in (2, 4, 10, 20)] + [(-2, 0, zt)] + [(2, -3, 0.0)] + [(4, -3, 0.0)] + [(10, -3, 0.0)])
    p.close()
    return img, anch, cores


def main(cfg=None):
    import matplotlib
    matplotlib.use('Agg')
    from functions import set_mpl_params
    import pyvista as pv
    pv.OFF_SCREEN = True
    c = dict(CFG)
    c.update(cfg or {})
    if LAYOUT == 'wide':
        K = c['ss_wide']
        c.update(ss=K, win=(2400 * K, 1300 * K), dpi_wide=300 * K)
        set_mpl_params(SF_WIDE, grid_on=False)      # 10 pt at print: drawn at 7.5 in, included at scale 1
        compose_wide(scene(1, c), c)
        return
    if LAYOUT == 'two':                             # two views stacked (193M above, RANS below), drawn 5.6 in wide, written as 4500 px = 7.5 in at 600 dpi
        K = c['ss_two']
        wpr = float(os.environ.get('OVERVIEW_W', c['w_two']))         # printed width (in); text stays 10 pt
        c.update(ss=K, win=(2400 * K, 1300 * K), dpi=c['dpi_two'], final_w=int(round(600 * wpr)), two_bars=True, lab_below=0.6, tick_sp=0.3, clear=0.08, cap_clear=0.15)     # WV labels close under the inset axis (leader about a third of the old length)
        set_mpl_params(SF_TWO * 7.5 / wpr, grid_on=False)               # 10 pt at print (design width 5.56 in, enlargement wpr / 5.56)
    else:
        c['top'] = 0                                # old order: RANS above
        set_mpl_params(SF, grid_on=False)           # 'tall': SF font size = 10 pt / include scale (figure is drawn 5.6 in wide, included at 3.56 in)
    res = [scene(0, c), scene(1, c)]
    compose(res, c)


_INS = {}


def compose_wide(res, c):
    """Full-text-width figure (7.5 in, scale 1): 193M 3D view on the left, the three cross-plane insets in a column at the right,
    two slim colour bars in the free lower-right corner of the 3D view."""
    import matplotlib.pyplot as plt
    from matplotlib import cm
    from matplotlib.colors import Normalize
    from PIL import Image
    im, an, _ = res
    K = c.get('ss', 1)                                  # render supersampling: image pixels per base pixel; all offsets below are in base pixels
    im = np.pad(im, ((300 * K, 0), (0, 0), (0, 0)), constant_values=255)     # room above the frame for the x/D = 20 tag
    an = (an + np.array([0, 300 * K])) / K
    W = 7.5
    cols = np.where((im < 250).any((0, 2)))[0]
    c0, c1 = cols.min() - 10 * K, cols.max() + 10 * K
    rr = np.where((im < 250).any((1, 2)))[0]
    r0, r1 = rr.min() - 60 * K, max(rr.max() + 20 * K, an[0][1] * K + 60 * K)
    wi = c['w_img']
    hp = wi * (r1 - r0) / (c1 - c0)
    H = max(hp, c['h_fig'])
    fig = plt.figure(figsize=(W, H))
    ax = fig.add_axes([0.0, (H - hp) / H, wi / W, hp / H])
    c0, c1, r0, r1 = c0 / K, c1 / K, r0 / K, r1 / K
    ax.imshow(im[int(r0 * K):int(r1 * K), int(c0 * K):int(c1 * K)], extent=(c0, c1, r1, r0))
    ax.set_xlim(c0, c1)
    ax.set_ylim(r1, r0)
    ax.axis('off')
    hs, c1p, c2p, t2, t4, t10, t20, tcl, s2p, s1p, s10p = an
    ax.text(0.0, 0.995, '193M 3DNS', transform=ax.transAxes, va='top', ha='left')
    ax.annotate('HV', xy=hs, xytext=(hs[0] - 150, hs[1] + 90), va='center', ha='center', arrowprops=dict(arrowstyle='->', lw=1.0, color='k'))
    mid = 0.5 * (c1p + c2p)
    lx, ly = mid[0] + c['cvp_dx'], mid[1] + c['cvp_dy']
    ax.text(lx, ly, 'CRVP', ha='center', va='bottom')
    for cp_ in (c1p, c2p):
        ax.annotate('', xy=cp_, xytext=(lx, ly + 12), arrowprops=dict(arrowstyle='->', lw=1.0, color='k', shrinkA=0, shrinkB=1))
    ax.text(t2[0] - 14, t2[1] + 6, r'$x/D=2$', ha='right', va='top')
    for t, s_ in ((t4, r'$x/D=4$'), (t20, r'$x/D=20$')):
        ax.text(t[0], t[1] - 15, s_, ha='center', va='bottom')
    ax.text(t10[0] - 12, t10[1] - 15, r'$x/D=10$', ha='right', va='bottom')
    ax.text(tcl[0] - 10, tcl[1], r'$y=0$', ha='right', va='center')
    fs = plt.rcParams['font.size']
    # inset column
    iw = c['in_w']
    ih = iw * c['zmax_in'] / 6.0
    x0 = W - iw - c['in_right']
    ys = c['in_ys']                                     # bottoms of the three insets (in), top to bottom: x/D 2, 4, 10
    for i, xp in enumerate((2, 4, 10)):
        o_, paths_, cs_, y_, z_ = inset_plane(1, c, xp, c['zmax_in'])
        ia = fig.add_axes([x0 / W, ys[i] / H, iw / W, ih / H])
        ia.pcolormesh(y_, z_, np.nan_to_num(o_).T, cmap='bwr', norm=Normalize(-c['omax'], c['omax']), shading='auto', rasterized=True)
        for pth in paths_:
            ia.plot(pth[:, 0], pth[:, 1], color='#202020', lw=0.6)
            m_ = len(pth) // 2
            dv = pth[min(m_ + 1, len(pth) - 1)] - pth[max(m_ - 1, 0)]
            dv = dv / (np.hypot(*dv) + 1e-12)
            ia.annotate('', xy=pth[m_] + 0.025 * dv, xytext=pth[m_] - 0.025 * dv,
                        arrowprops=dict(arrowstyle='-|>', lw=0, color='#202020', mutation_scale=5))
        ia.set_xlim(3, -3)
        ia.set_ylim(0, c['zmax_in'])
        ia.set_aspect('equal')
        ia.set_yticks([0, 2])
        ia.set_xticks([-3, 0, 3])
        ia.tick_params(labelsize=fs, length=2, pad=1.5)
        if i < 2:
            ia.set_xticklabels([])
        else:
            for lab_t, ha_t in zip(ia.get_xticklabels(), ('right', 'center', 'left')):
                lab_t.set_ha(ha_t)
            ia.set_xlabel(r'$y/D$', labelpad=1)
        ia.set_ylabel(r'$z/D$', labelpad=1)
        ia.text(0.0, 1.04, r'$x/D=%d$' % xp, transform=ia.transAxes, ha='left', va='bottom')
        for kind, q, side_ in label_items(1, xp, c):
            ia.annotate('WV', xy=(q[0], q[1] + 0.03), xytext=(q[0], c['wv_z']), ha='center', va='center', zorder=6,
                        bbox=dict(facecolor='w', edgecolor='none', alpha=0.8, pad=0.8),
                        arrowprops=dict(arrowstyle='-|>', lw=0.6, color='k', mutation_scale=6, shrinkA=0, shrinkB=0))
    cb = c['cb_box']                                    # (x, y, w, h, gap) in inches of the first bar; second bar above it
    for k, (nm, cmp_, lab) in enumerate((('eta', 'jet', r'$\eta$'), ('om', 'bwr', r'$\bar{\omega}_x D/U_\infty$'))):
        cax = fig.add_axes([cb[0] / W, (cb[1] + k * cb[4]) / H, cb[2] / W, cb[3] / H])
        nrm = Normalize(0, VMAX) if nm == 'eta' else Normalize(-c['omax'], c['omax'])
        cbar = fig.colorbar(cm.ScalarMappable(norm=nrm, cmap=cmp_), cax=cax, orientation='horizontal')
        if nm == 'om':
            cbar.set_ticks([-0.1, 0, 0.1])
        cax.tick_params(length=2, pad=1.5)
        fig.text((cb[0] - 0.08) / W, (cb[1] + k * cb[4] + cb[3] / 2) / H, lab, ha='right', va='center')
    fig.savefig(c['out'], dpi=c['dpi_wide'])
    plt.close(fig)
    print('saved', c['out'], 'figsize in: %.2f x %.2f' % (W, H))


def label_items(which, xp, c, zmax=None):
    """WV and OV labels of one inset (the CRVP is labelled in the 3D view only): [(kind, (y, z, omega_x), side)]."""
    mside = 1 if xp == 10 else -1
    zmax = c['zmax_in'] if zmax is None else zmax
    cs = [q for q in inset_plane(which, c, xp, zmax)[2] if q[1] <= 1.0]       # near-wall cores only (noise aloft is not labelled)
    return [it for it in classify_cores(cs) if it[0] == 'wall' and abs(it[1][2]) >= c['label_thr']]      # WV only (OV is not labelled)


def inset_plane(which, c, xp=4, zmax=0.8):
    """Data of the x/D xp plane for the 2D inset: smoothed omega_x (193M), sparse streamlines for z/D <= zmax, vortex cores."""
    key = (which, xp, zmax)
    if key not in _INS:
        _INS[key] = _inset_plane(which, c, xp, zmax)
    return _INS[key]


def _inset_plane(which, c, xp=4, zmax=0.8):
    from scipy.ndimage import gaussian_filter
    d = data()
    o, v, w, y1, z1 = d['planes'][xp]['R' if which == 0 else 'D']
    zi = z1 <= zmax + 1e-9
    z1, o, v, w = z1[zi], np.array(o[:, zi], float), np.array(v[:, zi], float), np.array(w[:, zi], float)
    cs = find_cores(o, y1, z1, c['core_thr'], c['core_sm'] if which == 1 else 0, zmax)       # cores from the unsmoothed field
    osm = gaussian_filter(np.nan_to_num(o), c['core_sm']) if which == 1 else np.nan_to_num(o)
    for side in (-1, 1):                                # fallback: main core = strongest correctly signed value in |y| < 1.3, z >= 0.15
        if any(q[0] * side > 0 and abs(q[0]) < 1.3 and q[2] * (-side) > 0 for q in cs):
            continue
        iy, iz = np.where((y1[:, None] * side > 0) & (np.abs(y1[:, None]) < 1.3) & (z1[None, :] >= 0.15))
        if len(iy):
            k_ = np.argmax((-side) * osm[iy, iz])
            if (-side) * osm[iy[k_], iz[k_]] >= c['core_thr']:
                cs.append((y1[iy[k_]], z1[iz[k_]], osm[iy[k_], iz[k_]]))
    if which == 1:
        m = np.isnan(o)
        o = gaussian_filter(np.nan_to_num(o), c['sm'])
        o[m] = np.nan
        v, w = gaussian_filter(np.nan_to_num(v), c['sm_uv']), gaussian_filter(np.nan_to_num(w), c['sm_uv'])
    seeds = []
    for yc, zc, _ in cs:
        for r in (0.08, 0.16):
            for t in np.linspace(0, 2 * np.pi, 6, endpoint=False):
                seeds.append((yc + r * np.cos(t), zc + 0.5 * r * np.sin(t)))
    seeds += [(a_, b_) for a_ in np.arange(-2.9, 3.0, 0.4) for b_ in np.arange(0.08, zmax, 0.16)]
    paths = trace_lines(y1, z1, v, w, seeds, c['inset_cell'] * (c['cell2'] if xp == 2 else 1.0), zmax, c['smin'], ds=0.02, lmax=2.5, lmin=0.2)
    return o, paths, cs, y1, z1


def classify_cores(cs):
    """Label rule, identical for both models and planes. Main vortex of a side = strongest core inside |y/D| < 1.3 whose sign is that of
    the main vortex of the side (omega_x < 0 for y > 0, > 0 for y < 0). Outboard cores (|y| > |y_main| + 0.3), taken in lateral order:
    the first with the opposite sign is the wall vortex; the first same-sign core farther out than the wall vortex is the outer vortex.
    Returns [(kind, (y, z, omega_x), side)]."""
    out = []
    for side in (-1, 1):
        sm_ = -side
        sc = [q for q in cs if q[0] * side > 0]
        mains = [q for q in sc if abs(q[0]) < 1.3 and q[2] * sm_ > 0]
        if not mains:
            continue
        m = max(mains, key=lambda q: abs(q[2]))
        out.append(('main', m, side))
        outb = sorted([q for q in sc if abs(q[0]) > abs(m[0]) + 0.3], key=lambda q: abs(q[0]))
        wall = next((q for q in outb if q[2] * sm_ < 0), None)
        if wall is not None:
            out.append(('wall', wall, side))
            outer = next((q for q in outb if abs(q[0]) > abs(wall[0]) + 0.15 and q[2] * sm_ > 0), None)
            if outer is not None:
                out.append(('outer', outer, side))
    return out


def compose(res, c):
    import matplotlib.pyplot as plt
    from matplotlib import cm
    from matplotlib.colors import Normalize
    W = 5.6
    K = c.get('ss', 1)                                  # render supersampling: all pixel quantities below are in base pixels (image px / K)
    res = [(im, an / K, cr) for im, an, cr in res]
    cols = [np.where((im < 250).any((0, 2)))[0] / K for im, _, _ in res]
    c0, c1 = min(q.min() for q in cols) - 10, max(q.max() for q in cols) + 10
    crops = []                                          # tight vertical crop of each 3D render (px): content, plus the horseshoe label
    for im, an, _ in res:
        rr = np.where((im < 250).any((1, 2)))[0] / K
        crops.append((rr.min() - c['pad_top'], max(rr.max() + c['pad_bot'], an[0][1] + c['pad_hs'])))
    hpx = max(b - a for a, b in crops)
    wi_ = c['img_w']                                    # width of the 3D image (in); the figure is W wide
    hp = wi_ * hpx / (c1 - c0)                          # panel height (in) so that the image fills the axes
    # inset layout: each inset sits just below-right of the lower near corner of its own plane, pushed down where needed so that the
    # three (each with a label row above) do not overlap; positions in inches relative to the bottom-left of the panel's 3D axes
    wi_in = c['inset_w']
    hi_in = wi_in * c['zmax_in'] / 6.0
    lay = []
    for k_, (im_, an_, _) in enumerate(res):
        r1o = crops[k_][0] + hpx
        A = [((a[0] - c0) / (c1 - c0) * wi_, (r1o - a[1]) / hpx * hp) for a in (an_[8], an_[9], an_[10])]
        # three insets side by side across the width, each raised along the lower (near) edge of the plate: the edge is the straight
        # line through the lower near corners of the three planes; an inset must lie below the edge at its left end (inset 1: below the
        # lower-left corner of the plate, i.e. below the 3D axes)
        lb, la = np.polyfit([a_[0] for a_ in A], [a_[1] for a_ in A], 1)
        gx = (5.47 - c['x_first'] - 3 * wi_in) / 2.0
        row = []
        for i_ in range(3):
            xl = c['xs_in'][i_]
            top = -c['clr'] if i_ == 0 else la + lb * xl - c['clr']
            row.append((xl, top - hi_in))
        lay.append(row)
        # check: no leader may pass over another inset (its axes, tick labels and WV/OV labels): boxes in inches, leader sampled densely
        boxes = [(xl - 0.16, yb - 0.62, xl + wi_in + 0.1, yb + hi_in + 0.1) for xl, yb in row]
        nbad = 0
        for i_ in range(3):
            xl, yb = row[i_]
            cx_ = min(max(A[i_][0], xl + 0.05), xl + wi_in - 0.05)
            for j_ in range(3):
                if j_ == i_:
                    continue
                for t_ in np.linspace(0, 1, 400):
                    px_, py_ = A[i_][0] + (cx_ - A[i_][0]) * t_, A[i_][1] + (yb + hi_in - A[i_][1]) * t_
                    if boxes[j_][0] <= px_ <= boxes[j_][2] and boxes[j_][1] <= py_ <= boxes[j_][3]:
                        nbad += 1
                        print('LEADER of inset %d crosses inset %d (panel %d)' % (i_, j_, k_))
                        break
        print('panel %d leader/inset check: %d crossings' % (k_, nbad), 'anchors', [(round(a_[0], 2), round(a_[1], 2)) for a_ in A],
              'insets', [(round(p_[0], 2), round(p_[1], 2)) for p_ in row])
    hr = max(0.05, c['tick_sp'] - min(p_[1] for l_ in lay for p_ in l_))      # room below the 3D axes for the lowest inset
    print('inset bottoms', [[round(p_[1], 2) for p_ in l_] for l_ in lay])
    mt = 0.15                                           # top margin so that frame labels at the very top are not clipped
    Hp = hp + hr + mt                                   # each panel (3D view + its insets) is drawn in its own figure and nested afterwards
    figs = [plt.figure(figsize=(W, Hp)) for _ in range(2)]
    axs = [f_.add_axes([0.0, hr / Hp, wi_ / W, hp / Hp]) for f_ in figs]
    for k, (ax, (im, an, _cores), lab) in enumerate(zip(axs, res, ('(b) RANS 50M', '(a) 193M 3DNS') if c['top'] == 1 else ('(a) RANS 50M', '(b) 193M 3DNS'))):
        pd_ = 400
        im = np.pad(im, ((pd_ * K, pd_ * K), (0, 0), (0, 0)), constant_values=255)
        an = an + np.array([0, pd_])
        r0 = crops[k][0] + pd_
        r1 = r0 + hpx
        ax.imshow(im[int(round(r0 * K)):int(round(r1 * K)), int(round(c0 * K)):int(round(c1 * K))], extent=(c0, c1, r1, r0))
        ax.set_xlim(c0, c1)
        ax.set_ylim(r1, r0)
        ax.axis('off')
        if k == c['top']:                                # upper panel: caption at its top left
            ax.text(c['tag_xy'][0], c['tag_xy'][1], lab, transform=ax.transAxes, va='top', ha=c['tag_ha'])
        else:
            cap_b = lab                                  # the lower caption is placed after nesting, at the first free spot top-left
        hs, c1p, c2p, t2, t4, t10, t20, tcl, s2p, s1p, s10p = an
        kw = dict(arrowprops=dict(arrowstyle='->', lw=1.0, color='k'), va='center', ha='center',
                  )
        ax.annotate('HV', xy=hs, xytext=(hs[0] - 150, hs[1] + 120), **kw)
        fig, H = figs[k], Hp
        pbot = hr                                                    # bottom of this 3D axes (in)
        fs = plt.rcParams['font.size']                  # same size as all other text of the figure (rcParams, set by sf)
        for i_in, (xp, zmax, anchor, side, mside) in enumerate(
                ((2, c['zmax_in'], s2p, 'left', -1), (4, c['zmax_in'], s1p, 'left', -1), (10, c['zmax_in'], s10p, 'left', 1))):
            wi = c['inset_w']
            hi = wi * zmax / 6.0
            x0in, yb_ = lay[k][i_in]
            ibot = pbot + yb_
            o_, paths_, cs_, y_, z_ = inset_plane(k, c, xp, zmax)
            ia = fig.add_axes([x0in / W, ibot / H, wi / W, hi / H])
            ia.pcolormesh(y_, z_, np.nan_to_num(o_).T, cmap='bwr', norm=Normalize(-c['omax'], c['omax']), shading='auto', rasterized=True)
            for pth in paths_:
                ia.plot(pth[:, 0], pth[:, 1], color='#202020', lw=0.6)
                m_ = len(pth) // 2
                dv = pth[min(m_ + 1, len(pth) - 1)] - pth[max(m_ - 1, 0)]
                dv = dv / (np.hypot(*dv) + 1e-12)
                ia.annotate('', xy=pth[m_] + 0.025 * dv, xytext=pth[m_] - 0.025 * dv,
                            arrowprops=dict(arrowstyle='-|>', lw=0, color='#202020', mutation_scale=5))
            ia.set_xlim(3, -3)                           # y/D increases to the left: the plane seen looking downstream, as in the 3D view
            ia.set_ylim(0, zmax)
            ia.set_aspect('equal')
            ia.set_xticks([-3, 0, 3])
            ia.set_yticks([0, 2])
            ia.tick_params(labelsize=fs, length=2, pad=1.5)
            for lab_t, ha_t in zip(ia.get_xticklabels(), ('right', 'center', 'left')):      # end labels inside, clear of the z tick labels
                lab_t.set_ha(ha_t)
            if i_in == 0:                                # y/D and z/D titles once, on the first (lowest-left) inset
                ia.set_xlabel(r'$y/D$', fontsize=fs, labelpad=1)
                ia.set_ylabel(r'$z/D$', fontsize=fs, labelpad=1)
            xa = 0.0                                     # x tick label at y/D = 0: labels keep clear of it
            # leader: from the plane's lower near corner straight down to the top edge of its inset (nearest point)
            ax_in = (anchor[0] - c0) / (c1 - c0) * wi_
            cx = min(max(ax_in, x0in + 0.05), x0in + wi - 0.05)
            ax.annotate('', xy=anchor, xycoords='data', xytext=(cx / W, (ibot + hi) / H), textcoords='figure fraction',
                        arrowprops=dict(arrowstyle='-', lw=0.6, color='0.25', shrinkA=0, shrinkB=0))
            kk = dict(fontsize=fs, ha='center', va='top', arrowprops=dict(arrowstyle='-|>', lw=0.6, color='k', mutation_scale=6))
            # labels follow the sign relative to the main vortex of that side and the lateral order (classify_cores)
            items = label_items(k, xp, c)
            xt = np.array([it[1][0] for it in items] + [xa], float)          # last entry: leader position (fixed)
            wd = np.array([{'main': 1.7, 'wall': 0.9, 'outer': 0.9}[it[0]] for it in items] + [0.3])   # text widths (D)
            for _ in range(300):                        # relax so that neighbouring texts (and the leader) are >= lab_gap D apart
                order = np.argsort(xt)
                for a_, b_ in zip(order[:-1], order[1:]):
                    dd = xt[b_] - xt[a_]
                    fa, fb = a_ == len(xt) - 1, b_ == len(xt) - 1
                    gp = 0.5 * (wd[a_] + wd[b_]) + c['lab_pad']
                    if dd < gp:
                        sh = gp - dd + 1e-3
                        if fa:
                            xt[b_] += sh
                        elif fb:
                            xt[a_] -= sh
                        else:
                            xt[a_] -= 0.5 * sh
                            xt[b_] += 0.5 * sh
                xt[:-1] = np.clip(xt[:-1], -c['lab_clip'], c['lab_clip'])
            tgt = np.array([it[1][0] for it in items], float)       # arrow tips: pushed apart where two cores are closer than 0.7 D
            for _ in range(50):
                o_ = np.argsort(tgt)
                for a_, b_ in zip(o_[:-1], o_[1:]):
                    if tgt[b_] - tgt[a_] < 0.7:
                        sh = 0.7 - (tgt[b_] - tgt[a_]) + 1e-3
                        tgt[a_] -= 0.5 * sh
                        tgt[b_] += 0.5 * sh
            for i, (kind, q, side_) in enumerate(items):
                x_ = xt[i]
                print('label %-5s %s  y/D x/D=%d %s: sign %+d  y/D %+.2f  z/D %.2f  omega_x D/U %+.3f' % (
                    kind, ('RANS', '193M')[k], xp, 'y<0' if side_ < 0 else 'y>0', np.sign(q[2]), q[0], q[1], q[2]))
                txt = {'main': 'CRVP', 'wall': 'WV', 'outer': 'OV'}[kind]
                ia.annotate(txt, xy=(tgt[i], q[1]), xytext=(x_, -c['lab_below']), textcoords='data', xycoords='data',
                            annotation_clip=False, **kk)
        mid = 0.5 * (c1p + c2p)
        lx, ly = mid[0] + c['cvp_dx'], mid[1] + c['cvp_dy']                # label above the plate (sky), no box
        ax.text(lx, ly, 'CRVP', ha='center', va='bottom')
        for cp_ in (c1p, c2p):
            ax.annotate('', xy=cp_, xytext=(lx, ly + 12), arrowprops=dict(arrowstyle='->', lw=1.0, color='k', shrinkA=0, shrinkB=1))
        ax.text(t2[0] - 14, t2[1] + 6, r'$x/D=2$', ha='right', va='top')        # left of its frame, clear of the inset above
        for t, s in ((t4, r'$x/D=4$'), (t20, r'$x/D=20$')):
            ax.text(t[0], t[1] - 15, s, ha='center', va='bottom')
        ax.text(t10[0] - 12, t10[1] - 15, r'$x/D=10$', ha='right', va='bottom')      # left of its frame, clear of the CRVP label
        ax.text(tcl[0] - 10, tcl[1], r'$y=0$', ha='right', va='center')
    H = c['hcb']
    fig = plt.figure(figsize=(W, H))                    # colour bars in their own strip
    three = c['color_by'] != 'omx'                      # eta, omega_x (planes, insets) and the isosurface quantity
    xs_ = [0.04, 0.37, 0.70] if three else [0.0, 0.0]
    # the two bars are stacked at the right, so the strip can rise into the white wedge below the last inset of the lower panel
    cax = fig.add_axes([0.10, 0.62 / H, 0.36, 0.12 / H])                  # side by side, spanning the width
    fig.colorbar(cm.ScalarMappable(norm=Normalize(0, VMAX), cmap='jet'), cax=cax, orientation='horizontal', label=r'$\eta$')
    cax2 = fig.add_axes([0.54, 0.62 / H, 0.36, 0.12 / H])
    fig.colorbar(cm.ScalarMappable(norm=Normalize(-c['omax'], c['omax']), cmap='bwr'), cax=cax2, orientation='horizontal',
                 label=r'$\bar{\omega}_x D/U_\infty$')
    wb = 0.26
    if three and c['color_by'] == 'sign':               # two-patch key for the isosurface colour
        from matplotlib.patches import Rectangle
        cax3 = fig.add_axes([xs_[2], 0.40 / H, wb, 0.12 / H])
        cax3.add_patch(Rectangle((0, 0), 0.5, 1, fc='#0000ff', ec='k', lw=0.8))
        cax3.add_patch(Rectangle((0.5, 0), 0.5, 1, fc='#ff0000', ec='k', lw=0.8))
        cax3.set_xlim(0, 1)
        cax3.set_ylim(0, 1)
        cax3.set_xticks([0.25, 0.75])
        cax3.set_xticklabels([r'$<0$', r'$>0$'])
        cax3.set_yticks([])
        cax3.set_xlabel(r'sign of $\bar{\omega}_x$ (surface)')
    elif three:
        cax3 = fig.add_axes([xs_[2], 0.40 / H, wb, 0.12 / H])
        fig.colorbar(cm.ScalarMappable(norm=Normalize(*c['clim']), cmap='inferno'), cax=cax3, orientation='horizontal',
                     label=c['cb_label'])
    strips = [fig]
    if c.get('two_bars'):                               # 'two': each bar in its own strip, so that each rises into the nearest free space
        plt.close(fig)
        Hs = 0.62
        strips = []
        for x_, nrm_, cmp_, lab_ in ((0.71, Normalize(-c['omax'], c['omax']), 'bwr', r'$\bar{\omega}_x D/U_\infty$'),
                                     (0.44, Normalize(0, VMAX), 'jet', r'$\eta$')):
            fs_ = plt.figure(figsize=(W, Hs))
            cx_ = fs_.add_axes([x_, 0.40 / Hs, 0.22, 0.12 / Hs])           # side by side in the free lower-right corner
            fs_.colorbar(cm.ScalarMappable(norm=nrm_, cmap=cmp_), cax=cx_, orientation='horizontal', label=lab_)
            strips.append(fs_)
    nest_save([figs[c['top']], figs[1 - c['top']]] + strips, c, caption=cap_b)


def nest_save(figs, c, caption=None):
    """Rasterise the panel figures and the colour-bar strip, then stack them as tightly as their content allows: each lower item is moved
    up until its non-white pixels (grown by a clearance margin) would touch those of the items above. The white margins are cropped."""
    import io
    import matplotlib.pyplot as plt
    from PIL import Image
    from scipy.ndimage import binary_dilation
    dpi = c['dpi']
    arrs = []
    for f_ in figs:
        b = io.BytesIO()
        f_.savefig(b, dpi=dpi, format='png', facecolor='white')
        plt.close(f_)
        arrs.append(np.array(Image.open(b).convert('RGB')))
    s = 4                                               # mask resolution reduction
    m_px = int(c['clear'] * dpi / s)                    # clearance (in) in reduced pixels
    def mask(a, clear=None):
        m = (a < 250).any(2)
        m = m[:m.shape[0] // s * s, :m.shape[1] // s * s].reshape(m.shape[0] // s, s, m.shape[1] // s, s).any((1, 3))
        return binary_dilation(m, iterations=m_px if clear is None else int(clear * dpi / s))
    def stack(canvas, a, off):
        h = max(canvas.shape[0], off + a.shape[0])
        new = np.full((h, canvas.shape[1], 3), 255, np.uint8)
        new[:canvas.shape[0]] = canvas
        new[off:off + a.shape[0]] = np.minimum(new[off:off + a.shape[0]], a)
        return new

    def highest(cm_, ma):                               # highest row (reduced) where `ma` touches nothing of `cm_`
        top = cm_.shape[0]
        while top > 0:
            n = min(ma.shape[0], cm_.shape[0] - (top - 1))
            if n > 0 and (cm_[top - 1: top - 1 + n] & ma[:n]).any():
                break
            top -= 1
        return top

    canvas = arrs[0]
    cm_ = mask(canvas)
    cap = mc = None
    if caption:                                          # caption of the lower panel: placed separately, with a larger gap
        fc = plt.figure(figsize=(1.8, 0.34))
        fc.text(0.0, 0.92, caption, va='top', ha='left')
        bc = io.BytesIO()
        fc.savefig(bc, dpi=dpi, format='png', facecolor='white')
        plt.close(fc)
        cap = np.array(Image.open(bc).convert('RGB'))
        mc = mask(cap, c['cap_clear'])
    for i_a, a in enumerate(arrs[1:]):
        ma = mask(a)
        top = highest(cm_, ma)
        if i_a == 0 and cap is not None:                # lower the panel in small steps until its caption fits top-left, near the panel top
            while True:
                tryc = stack(canvas, a, top * s)
                tm = mask(tryc)
                row0 = (top * s + int(0.15 * dpi)) // s
                spot = None
                for xr in range(0, tm.shape[1] - mc.shape[1]):
                    for dy in range(0, int(c['cap_dy'] * dpi / s)):
                        r_ = row0 + dy
                        reg = tm[r_: r_ + mc.shape[0], xr: xr + mc.shape[1]]
                        if not (reg & mc[:reg.shape[0], :reg.shape[1]]).any():
                            spot = (r_, xr)
                            break
                    if spot or xr > int(c['cap_xmax'] * dpi / s):
                        break
                if spot:
                    break
                top += 2
            canvas = tryc.copy()
            r_c, x_c = spot[0] * s, spot[1] * s
            h_c = min(cap.shape[0], canvas.shape[0] - r_c)
            canvas[r_c: r_c + h_c, x_c: x_c + cap.shape[1]] = np.minimum(canvas[r_c: r_c + h_c, x_c: x_c + cap.shape[1]], cap[:h_c])
        else:
            canvas = stack(canvas, a, top * s)
        cm_ = mask(canvas)
    nz = np.where((canvas < 250).any(2))
    pad = int(0.05 * dpi)
    r0, r1 = max(nz[0].min() - pad, 0), nz[0].max() + pad
    c0_, c1_ = max(nz[1].min() - pad, 0), nz[1].max() + pad
    out = canvas[r0:r1, c0_:c1_]
    if c.get('final_w'):                                 # 'two': written as final_w px wide = 7.5 in at 600 dpi (design width x dpi is within 1 % of it)
        fw = c['final_w']
        print('design width %.3f in; print enlargement %.3f; font %.2f pt' % (out.shape[1] / dpi, fw / 600 / (out.shape[1] / dpi),
                                                                           plt.rcParams['font.size'] * fw / 600 / (out.shape[1] / dpi)))
        out = np.array(Image.fromarray(out).resize((fw, int(round(out.shape[0] * fw / out.shape[1]))), Image.LANCZOS))
        Image.fromarray(out).save(c['out'], dpi=(600, 600))
        print('saved', c['out'], 'size in: %.2f x %.2f' % (out.shape[1] / 600, out.shape[0] / 600))
        return
    Image.fromarray(out).save(c['out'], dpi=(dpi, dpi))
    print('saved', c['out'], 'size in: %.2f x %.2f' % (out.shape[1] / dpi, out.shape[0] / dpi), 'blank fraction: %.3f' % (out >= 250).all(2).mean())


if __name__ == '__main__':
    if '--build' in sys.argv:
        build_rans()
    else:
        main()
