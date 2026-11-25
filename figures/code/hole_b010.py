"""Hole-section table and centre-plane / section rasters of the 193M 3DNS (C193b) at stat3 batch b010.

Usage: python3 hole_b010.py     (no arguments; run once, later calls reuse the cache)

Method of analysis/hole_discharge/extract.py and plots.py (4dns_filmcooled, imported, not edited):
hole blocks only (case_info.json 'hole'), node means -> hex means, Favre velocity; sections every 0.2 D normal to
the 30 deg hole axis; rasters on the 0.01 D pixel grid of plots.py F4 (centre plane) and F5 (section s/D = -1.1).
Run folder /mnt/Elements/DNS/film_cooling/C193b_BR15, files stat3_<n>_b010 (pictures_job18/stat3_b010 holds the
post-processing of the same batch).
Cache (Elements): /mnt/Elements/DNS/analysis/asme2027_figs/hole_b010/
  cells_c193b_b010.npz, sections_c193b_b010.txt, redraw_cache_c193b_b010.npz (v1, v2, ok2 as plots.py).
"""
import os
import sys

import numpy as np
import pandas as pd

REPO = '/home/sk999/Documents/4dns_filmcooled'
sys.path.insert(0, REPO + '/analysis/hole_discharge')
sys.path.insert(0, REPO + '/starccm/python')
sys.path.insert(0, REPO + '/dns/python')
from hole_loss import dns_cells, sections  # noqa: E402
import plots as P  # noqa: E402

RUN = '/mnt/Elements/DNS/film_cooling/C193b_BR15'
BATCH = 10
CACHE = '/mnt/Elements/DNS/analysis/asme2027_figs/hole_b010'
SECTIONS = CACHE + '/sections_c193b_b010.txt'
RASTER = CACHE + '/redraw_cache_c193b_b010.npz'
CELLS = CACHE + '/cells_c193b_b010.npz'


def build():
    os.makedirs(CACHE, exist_ok=True)
    if os.path.exists(CELLS):
        c = dict(np.load(CELLS))
    else:
        d = dns_cells(RUN, BATCH)
        c = {k: d[k].to_numpy() for k in d.columns}
        np.savez(CELLS, **c)
        del d
    if not os.path.exists(SECTIONS):
        S = sections(pd.DataFrame(c), 0.2)
        S.to_csv(SECTIONS, index=False, float_format='%.6g')
    if not os.path.exists(RASTER):
        D = P.D
        gx, gz = np.arange(-4.0, 3.0, 0.01), np.arange(-3.6, 0.01, 0.01)
        X, Z = np.meshgrid(gx, gz)
        P1 = np.column_stack([X.ravel() * D, np.zeros(X.size), Z.ravel() * D])
        gn, gy = np.arange(-1.6, 1.6, 0.01), np.arange(-2.2, 2.21, 0.01)
        N, Y = np.meshgrid(gn, gy)
        s0 = -1.1
        P2 = np.column_stack([(s0 * P.C30 - N.ravel() * P.S30) * D, Y.ravel() * D,
                              (s0 * P.S30 + N.ravel() * P.C30) * D])
        inside = P.hole_mask(X, Z)
        v1 = P.raster(c, P1, tol=0).reshape(X.shape)
        v1[~inside] = np.nan
        v2, ok = P.raster_ok(c, P2)
        np.savez(RASTER, v1=v1, v2=v2.reshape(N.shape), ok2=ok.reshape(N.shape))


def summary():
    d = pd.read_csv(SECTIONS)
    at = lambda s0, col: float(d.loc[(d.s - s0).abs().idxmin(), col])  # noqa: E731
    diff = d[(d.s > -3.0) & (d.s < -0.8)]
    pm, pe = at(-4.9, 'p'), at(-0.9, 'p')
    return {'BR': diff.mdot.median() / 3.664e-3, 'rev_peak': d.rev.max(), 's_rev_peak': float(d.s[d.rev.idxmax()]),
            'rise_Pa': pe - pm, 'p_m': pm, 'p_exit': pe, 'rev_exit': at(-0.9, 'rev')}


if __name__ == '__main__':
    build()
    print(summary())
