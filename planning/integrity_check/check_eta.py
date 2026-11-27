"""Read-only check of the eta-bar statements: RANS deficit range beyond x/D 7 and scale-resolving agreement.

Uses the same sources as figures/code/fig_eta_x.py (RANS wall_plate.csv, 3DNS eta_lateral.txt); writes nothing.
Run with python3 -B.
"""
import os
import sys

import numpy as np

REPO = '/home/sk999/Documents/4dns_filmcooled'
sys.path.insert(0, REPO + '/starccm/python')
from bl_post import D, T0_FS, T0_C, read_csv, lateral_mean  # noqa: E402

RANS = {'RANS 5M': REPO + '/starccm/rans_5M_Lm25_BR15_FST_ZPG_p0dns',
        'RANS 50M': REPO + '/starccm/rans_50M_Lm25b_BR15_FST_ZPG_p0dns'}
DNS = {'20M LES': ('/mnt/Tank/DNS/film_cooling/C20bW_BR15', 'pictures_job3', 6),
       '50M': ('/mnt/Tank/DNS/film_cooling/C50b_BR15', 'pictures_job9', 6),
       '100M': ('/mnt/Elements/DNS/film_cooling/C100b_BR15', 'pictures_job9', 6),
       '193M': ('/mnt/Elements/DNS/film_cooling/C193b_BR15', 'pictures_job18', 10)}


def lat(r):
    if r in RANS:
        w = read_csv(RANS[r] + '/post_run1/wall_plate.csv')
        T = w['Temperature']
        xc, la = lateral_mean(w['X'], w['Y'], (T0_FS - T) / (T0_FS - T0_C), np.arange(-2 * D, 25.5 * D, 0.1 * D))
        return xc / D, la
    run, pic, b = DNS[r]
    d = np.loadtxt(os.path.join(run, pic, f'stat3_b{b:03d}', 'eta_lateral.txt'))
    return d[:, 0], d[:, 1]


xr, er = lat('193M')
for r in ('RANS 5M', 'RANS 50M', '20M LES', '50M', '100M'):
    x, e = lat(r)
    m = (x >= 2.5) & (x <= 20)
    ei = np.interp(x[m], xr, er)
    d = e[m] - ei
    m7 = (x[m] >= 7)
    low = 1 - e[m][m7] / ei[m7]
    print(f'{r:9s} max|diff| x/D 2.5-20 {np.nanmax(np.abs(d)):.4f} at {x[m][np.nanargmax(np.abs(d))]:.2f};'
          f' deficit x/D 7-20 {100 * np.nanmin(low):.1f}-{100 * np.nanmax(low):.1f} %')
