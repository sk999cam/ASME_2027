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
D = 7.75e-3
UINF = 86.28
NORM = (UINF / D) ** 2
ZCUT = 0.1
VMAX = 0.8
ZT = 3.0
HS = (-1.05, 0.0, 0.1)                          # horseshoe vortex anchor
CVP = [(5.0, -0.6, 0.4), (5.0, 0.6, 0.4)]       # vortex pair anchors

# final choice (see candidates in the scratchpad f4_candidates)
CFG = dict(cam=[(-15, -8, 10), (9, 0, 0.4), (0, 0, 1)], zoom=1.0, qlev=0.02, xcut=14.0, omax=0.3,
           op_iso=0.6, plane_mode='white', floor=0.1, zp=1.5, op_wall=1.0, planes=(4, 10, 20), win=(2400, 1300), dpi=300, out=OUT)


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
    _DATA.update(V=V, R=R, hx=hx, hy=hy, Xw=Xw, Yw=Yw, eta=(eta_r, eta_d), inh=inh, planes=pl)
    return _DATA


def streamlines(y, z, v, w, density=1.3, ylim=1.5, zlim=1.5):
    """2D streamlines of (v, w) inside |y| <= ylim, z <= zlim; returns the segments (n, 2, 2) in (y, z)."""
    import matplotlib.pyplot as plt
    iy = np.abs(y) <= ylim + 1e-9
    iz = z <= zlim + 1e-9
    y, z, v, w = y[iy], z[iz], v[np.ix_(iy, iz)], w[np.ix_(iy, iz)]
    fig, ax = plt.subplots()
    vm = np.ma.masked_invalid(v.T)
    wm = np.ma.masked_invalid(w.T)
    sp = ax.streamplot(y, z, vm, wm, density=density, integration_direction='both', broken_streamlines=False)
    seg = np.array(sp.lines.get_segments())
    plt.close(fig)
    return seg


def polyline_from_segments(seg, xp):
    import pyvista as pv
    n = len(seg)
    pts = np.zeros((2 * n, 3))
    pts[:, 0] = xp
    pts[:, 1:] = seg.reshape(-1, 2)
    cells = np.c_[np.full(n, 2), np.arange(0, 2 * n, 2), np.arange(1, 2 * n, 2)].ravel()
    return pv.PolyData(pts, lines=cells)


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
    pl.point_data['eta'] = e.ravel(order='F')
    p.add_mesh(pl, scalars='eta', cmap='jet', clim=(0, VMAX), nan_color='white', show_scalar_bar=False, lighting=False,
               opacity=cfg['op_wall'])
    q = np.nan_to_num(Qn, nan=0.0)
    q[:, :, z < ZCUT] = 0.0
    q[V['x'] > cfg['xcut'], :, :] = 0.0
    g = pv.ImageData(dimensions=q.shape, spacing=(0.05,) * 3, origin=(V['x'][0], V['y'][0], V['z'][0]))
    g.point_data['Q'] = q.ravel(order='F')
    if which == 0:
        om = np.nan_to_num(np.array(d['R']['omx_n']), nan=0.0)
    else:
        om = np.nan_to_num(np.array(V['omx']), nan=0.0) * D / UINF
    g.point_data['om'] = om.ravel(order='F')
    iso = g.contour([cfg['qlev']], scalars='Q')
    p.add_mesh(iso, scalars='om', cmap='bwr', clim=(-omax, omax), smooth_shading=True, specular=0.2, opacity=cfg['op_iso'],
               show_scalar_bar=False)
    p.add_mesh(pv.lines_from_points(np.c_[d['hx'], d['hy'], np.full_like(d['hx'], 0.005)]), color='k', line_width=4)
    fr = [(-2, 0, 0), (20, 0, 0), (20, 0, ZT), (-2, 0, ZT), (-2, 0, 0)]
    p.add_mesh(pv.lines_from_points(np.array(fr, float)), color='k', line_width=3)
    # cross-planes: omega_x image (alpha rising with |omega_x|) and streamlines of (V, W)
    cmap = matplotlib.colormaps['bwr']
    for xp in (4, 10, 20):
        if xp in cfg['planes']:
            o, v, w, y1, z1 = d['planes'][xp]['R' if which == 0 else 'D']
            zi = z1 <= cfg['zp'] + 1e-9                       # plane image up to z/D zp
            y1, z1, o = y1, z1[zi], o[:, zi]
            Yp, Zp = np.meshgrid(y1, z1, indexing='ij')
            rgba = cmap(Normalize(-omax, omax)(np.nan_to_num(o)))
            if cfg['plane_mode'] == 'white':                  # opaque, white background
                rgba[..., 3] = 1.0
            else:                                             # opaque where |omega_x| > floor, clear elsewhere
                rgba[..., 3] = np.where(np.abs(np.nan_to_num(o)) > cfg['floor'] * omax, 1.0, 0.0)
            rgba[..., 3] = np.where(np.isfinite(o), rgba[..., 3], 0.0)
            sg = pv.StructuredGrid(np.full_like(Yp, xp), Yp, Zp)
            # StructuredGrid point order is x fastest: build the colour array in Fortran order
            sg.point_data['rgba'] = (np.stack([rgba[..., k].ravel(order='F') for k in range(4)], 1) * 255).astype(np.uint8)
            p.add_mesh(sg, scalars='rgba', rgba=True, lighting=False)
            if v is not None:
                seg = streamlines(y1, d['planes'][xp]['R' if which == 0 else 'D'][4], v, w)
                p.add_mesh(polyline_from_segments(seg, xp - 0.03), color='#202020', line_width=1.0, opacity=0.9)
        fr = [(xp, -3, 0), (xp, 3, 0), (xp, 3, ZT), (xp, -3, ZT), (xp, -3, 0)]
        p.add_mesh(pv.lines_from_points(np.array(fr, float)), color='#303030', line_width=3)
    p.camera_position = cfg['cam']
    p.camera.zoom(cfg['zoom'])
    p.render()
    img = p.screenshot(return_img=True)
    anch = project(p, [HS] + CVP + [(xp, 3, ZT) for xp in (4, 10, 20)] + [(-2, 0, ZT)])
    p.close()
    return img, anch


def main(cfg=None):
    import matplotlib
    matplotlib.use('Agg')
    from functions import set_mpl_params
    import pyvista as pv
    pv.OFF_SCREEN = True
    set_mpl_params(0.5, grid_on=False)
    c = dict(CFG)
    c.update(cfg or {})
    res = [scene(0, c), scene(1, c)]
    compose(res, c)


def compose(res, c):
    import matplotlib.pyplot as plt
    from matplotlib import cm
    from matplotlib.colors import Normalize
    fig, axs = plt.subplots(2, 1, figsize=(5.6, 7.4))
    cols = []
    for im, _ in res:
        cols.append(np.where((im < 250).any((0, 2)))[0])
    c0, c1 = min(c.min() for c in cols) - 10, max(c.max() for c in cols) + 10
    for k, (ax, (im, an), lab) in enumerate(zip(axs, res, ('(a) RANS 50M', '(b) 193M 3DNS'))):
        im = np.pad(im, ((70, 260), (0, 0), (0, 0)), constant_values=255)       # room for the labels
        an = an + np.array([0, 70])
        r, cc = np.where((im < 250).any(2))
        r0, r1 = r.min() - 70, r.max() + 200
        ax.imshow(im[r0:r1, c0:c1], extent=(c0, c1, r1, r0))
        ax.set_xlim(c0, c1)
        ax.set_ylim(r1, r0)
        ax.axis('off')
        ax.text(0.01, 0.97, lab, transform=ax.transAxes, va='top', ha='left')
        hs, c1p, c2p, t4, t10, t20, tcl = an
        kw = dict(arrowprops=dict(arrowstyle='->', lw=1.0, color='k'), va='center', ha='center',
                  bbox=dict(fc='w', ec='none', alpha=1.0, pad=1.5))
        ax.annotate('horseshoe\nvortex', xy=hs, xytext=(hs[0] + 60, hs[1] + 190), **kw)
        mid = 0.5 * (c1p + c2p)
        ax.annotate('counter-rotating\nvortex pair', xy=c2p, xytext=(mid[0] + 470, mid[1] + 250), **kw)
        ax.annotate('', xy=c1p, xytext=(mid[0] + 470, mid[1] + 250 - 45), arrowprops=kw['arrowprops'])
        for t, s in ((t4, r'$x/D=4$'), (t10, r'$x/D=10$'), (t20, r'$x/D=20$')):
            ax.text(t[0], t[1] - 15, s, ha='center', va='bottom', fontsize=7)
        ax.text(tcl[0] - 10, tcl[1], r'$y=0$', ha='right', va='center', fontsize=7)
    fig.subplots_adjust(left=0.01, right=0.84, top=0.995, bottom=0.005, hspace=0.02)
    cax = fig.add_axes([0.89, 0.55, 0.03, 0.4])
    fig.colorbar(cm.ScalarMappable(norm=Normalize(0, VMAX), cmap='jet'), cax=cax, label=r'$\eta$')
    cax2 = fig.add_axes([0.89, 0.06, 0.03, 0.4])
    fig.colorbar(cm.ScalarMappable(norm=Normalize(-c['omax'], c['omax']), cmap='bwr'), cax=cax2,
                 label=r'$\bar{\omega}_x D/U_\infty$')
    fig.savefig(c['out'], dpi=c['dpi'], bbox_inches='tight', pad_inches=0.05)
    plt.close(fig)


if __name__ == '__main__':
    if '--build' in sys.argv:
        build_rans()
    else:
        main()
