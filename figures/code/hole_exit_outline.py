"""Hole-exit (breakout) outline on the plate z = 0 of the 7-7-7 laidback fan-shaped hole, in D (paper plan-view figures).

Usage: python3 hole_exit_outline.py     (writes figures/code/data/hole_exit_outline.txt, columns x/D, y/D, closed polygon)
       from hole_exit_outline import load; x, y = load()

Source: hole-wall surface of the RANS 50M case (starccm/rans_50M_Lm25b_BR15_FST_ZPG_p0dns/post_run1/wall_hole.csv, read only),
the CAD geometry as meshed. The vertices with z/D > -0.01 are the rim of the hole at the plate. The rim is binned in 1 degree
about its centroid (radius = median), smoothed with a periodic Gaussian (sigma 1.5 bins) and written at 360 points.
D = 7.75 mm.
"""
import os

import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter1d

D = 7.75e-3
SRC = ('/home/sk999/Documents/4dns_filmcooled/starccm/rans_50M_Lm25b_BR15_FST_ZPG_p0dns/post_run1/wall_hole.csv')
OUTF = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data', 'hole_exit_outline.txt')


def build():
    w = pd.read_csv(SRC)
    x, y, z = (w[k + ' (m)'].to_numpy() / D for k in 'XYZ')
    m = z > -0.01
    x, y = x[m], y[m]
    xc, yc = 0.5 * (x.min() + x.max()), 0.0
    th = np.degrees(np.arctan2(y - yc, x - xc))
    r = np.hypot(x - xc, y - yc)
    b = np.clip(((th + 180) // 1).astype(int), 0, 359)
    rb = np.array([np.median(r[b == i]) if (b == i).any() else np.nan for i in range(360)])
    idx = np.arange(360)
    ok = ~np.isnan(rb)
    rb = np.interp(idx, idx[ok], rb[ok], period=360)
    rs = gaussian_filter1d(rb, 1.5, mode='wrap')
    t = np.radians(idx + 0.5 - 180)
    px, py = xc + rs * np.cos(t), yc + rs * np.sin(t)
    px, py = np.r_[px, px[0]], np.r_[py, py[0]]
    os.makedirs(os.path.dirname(OUTF), exist_ok=True)
    np.savetxt(OUTF, np.c_[px, py], fmt='%.5f', header='hole exit outline on z = 0, x/D y/D (closed); source ' + SRC)
    return px, py


def load():
    if not os.path.exists(OUTF):
        return build()
    a = np.loadtxt(OUTF)
    return a[:, 0], a[:, 1]


if __name__ == '__main__':
    px, py = build()
    print('x/D', px.min(), px.max(), 'y/D', py.min(), py.max())
