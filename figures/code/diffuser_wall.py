"""Wall-aligned and section rasters of the mean hole-axis velocity u_s for fig_diffuser_flow.py (cached).

Cells: RANS 50M and 20M LES (batch b006) from /mnt/Tank/DNS/analysis/hole_discharge/cells_<case>.npz (read only),
193M 3DNS (batch b010) from /mnt/Elements/DNS/analysis/asme2027_figs/hole_b010/cells_c193b_b010.npz.
The mesh near the wall is made of tilted cell columns spaced COLSP apart along s (RANS 50M 0.08 D, 20M LES 0.11 D,
193M 0.045 D, read from the cell-centre plots).
Section rasters at s/D = -2.3, -1.1, -0.5: linear interpolation in (y, n) of the cell centres with |s - s0| < COLSP
(one column spacing each side), NaN outside the convex hull of those centres, so only the fluid is coloured, to mesh
resolution. (The centre-plane wall raster of earlier versions is no longer used.)
Cache: /mnt/Elements/DNS/analysis/asme2027_figs/hole_b010/sec_<case>_<s0>.npz.
"""
import os
import sys

import numpy as np
import pandas as pd
from scipy import ndimage
from scipy.interpolate import LinearNDInterpolator
from scipy.spatial import cKDTree

REPO = '/home/sk999/Documents/4dns_filmcooled'
sys.path.insert(0, REPO + '/analysis/hole_discharge')
import plots as P  # noqa: E402

CACHE = '/mnt/Elements/DNS/analysis/asme2027_figs/hole_b010'
SRC = {'rans50m': '/mnt/Tank/DNS/analysis/hole_discharge/cells_rans50m.npz',
       'c20bw': '/mnt/Tank/DNS/analysis/hole_discharge/cells_c20bw.npz',
       'c193b': CACHE + '/cells_c193b_b010.npz'}
D, C30, S30 = P.D, P.C30, P.S30
SLAB = {'c20bw': 0.02}      # half-thickness (D) of the centre-plane cell layer; 20M LES: cell layers at |y| 0.004-0.02 D
GS = np.arange(-3.3, -0.0999, 0.01)
GD = np.arange(0.0, 0.1001, 0.001)
COLSP = {'rans50m': 0.08, 'c20bw': 0.11, 'c193b': 0.045}
GN, GY = np.arange(-1.6, 1.6, 0.01), np.arange(-2.2, 2.21, 0.01)


def lower_wall():
    """Downstream (n < 0) wall of the hole on y = 0: (s, n) arrays, as plots.hole_mask."""
    w = pd.read_csv(P.WALL)
    wx, wy, wz = (w[k + ' (m)'].to_numpy() / D for k in 'XYZ')
    m = np.abs(wy) < 0.03
    ws, wn = wx[m] * C30 + wz[m] * S30, -wx[m] * S30 + wz[m] * C30
    k = (ws > -3.6) & (ws < -0.3) & (wn < 0)
    b = np.round(ws[k] / 0.05)
    pts = np.array([(ws[k][b == u][np.argmin(wn[k][b == u])], wn[k][b == u].min()) for u in np.unique(b)])
    c = np.polyfit(pts[:, 0], pts[:, 1], 1)       # the wall is straight: max residual 0.008 D (bins of 0.1 D gave a staircase)
    sg = np.array([-5.0, 0.0])
    return sg, np.polyval(c, sg)


def sections_for(case, s_list):
    """{s0: v2 (GY x GN, NaN outside the fluid hull)} of u_s on the section s/D = s0 (cached per case and s0)."""
    fs = {s0: f'{CACHE}/sec_{case}_{s0:+.1f}.npz' for s0 in s_list}
    todo = [s0 for s0 in s_list if not os.path.exists(fs[s0])]
    if todo:
        c = dict(np.load(SRC[case]))
        us = c['u'] * C30 + c['w'] * S30
        sc = (c['x'] * C30 + c['z'] * S30) / D
        N, Y = np.meshgrid(GN, GY)
        for s0 in todo:
            sec = np.abs(sc - s0) < COLSP[case]
            pn = np.column_stack([c['y'][sec] / D, (-c['x'][sec] * S30 + c['z'][sec] * C30) / D])
            v2 = LinearNDInterpolator(pn, us[sec])(np.column_stack([Y.ravel(), N.ravel()])).reshape(N.shape)
            np.savez(fs[s0], v2=v2)
    return {s0: np.load(fs[s0])['v2'] for s0 in s_list}


if __name__ == '__main__':
    for k in SRC:
        sections_for(k, (-2.3, -1.1, -0.5))
        print(k, 'done', flush=True)
