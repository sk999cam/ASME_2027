#!/usr/bin/env python3
"""Q4/Q6: why the film volume alignment term (+0.57 to +0.81) sits with a centre-plane median cos of 0.83, and where
the points without 193M data lie. Uses pts.npz (scratchpad) with the coordinates of the point set (errorbars.points,
read-only) and the cached 193M planes (meff_planes/ref193_b010.npz).
usage: python3 -B q4b_where.py pts.npz
"""
import os
import sys

import numpy as np

os.environ['MEFF_REGIONS'] = 'pp'
sys.dont_write_bytecode = True
sys.path.insert(0, '/home/sk999/Documents/4dns_filmcooled/analysis/meff_c50b_c100b')
import errorbars as eb  # noqa: E402

D, U = 0.00775, 86.2812
SMIN = 0.05 * U / D
d = np.load(sys.argv[1])
names, P = eb.points(eb.CASES['C50b'], 2)
assert [str(s) for s in d['regnames']] == names and len(P['n']) == len(d['reg'])
assert np.array_equal(P['n'], d['blk'])
x, y, z = (P['xyz'][:, i] / D for i in range(3))
w = P['w'].astype(np.float64)
reg = P['reg']
c, S, k, L, rho = [d['c193_' + q].astype(np.float64) for q in ('cos', 'S', 'k', 'L', 'rho')]
fin = np.isfinite(c) & np.isfinite(S)


def wmed(a, ww):
    i = np.argsort(a)
    cw = np.cumsum(ww[i])
    return float(a[i][np.searchsorted(cw, 0.5 * cw[-1])])


print('--- extent of the points with no finite 193M value (interpolation not found or buffer NaN) ---')
for rn in ('wake_10_25', 'upstream_bl', 'wake_3_10', 'near_jet_0_3'):
    m = reg == names.index(rn)
    nf = m & ~fin
    print(f'{rn:12s} points {m.sum():8d}, non-finite {nf.sum():8d} ({w[nf].sum() / w[m].sum():.3f} of volume); '
          f'x/D of region [{x[m].min():.1f}, {x[m].max():.1f}], of finite [{x[m & fin].min():.1f}, {x[m & fin].max():.1f}], of non-finite '
          + (f'[{x[nf].min():.1f}, {x[nf].max():.1f}]' if nf.any() else '-'))
m = (reg == names.index('wake_10_25'))
for lo, hi in ((10, 15), (15, 20), (20, 22), (22, 24), (24, 30)):
    s = m & (x >= lo) & (x < hi)
    if s.any():
        print(f'   film 10-25, x/D {lo}-{hi}: volume share {w[s].sum() / w[m].sum():.3f}, finite fraction {w[s & fin].sum() / w[s].sum():.3f}')

print('\n--- cos distribution in the film volume by lateral position (mask: finite, S > Smin, cos > 0, k > 0) ---')
for rn in ('wake_3_10', 'wake_10_25'):
    m0 = (reg == names.index(rn)) & fin & (S > SMIN) & (c > 0) & (k > 0)
    print(rn)
    for lo, hi in ((0, 0.5), (0.5, 1), (1, 1.5), (1.5, 2), (2, 3.01)):
        s = m0 & (np.abs(y) >= lo) & (np.abs(y) < hi)
        if s.any():
            print(f'  |y|/D {lo}-{hi}: volume share {w[s].sum() / w[m0].sum():.3f}, median cos {wmed(c[s], w[s]):.3f}, <-ln cos> {np.average(-np.log(c[s]), weights=w[s]):.3f}')
    for lo, hi in ((0, 0.2), (0.2, 0.5), (0.5, 1), (1, 2.01)):
        s = m0 & (z >= lo) & (z < hi)
        if s.any():
            print(f'  z/D   {lo}-{hi}: volume share {w[s].sum() / w[m0].sum():.3f}, median cos {wmed(c[s], w[s]):.3f}, <-ln cos> {np.average(-np.log(c[s]), weights=w[s]):.3f}')
    s = m0 & (np.abs(y) < 0.3)
    print(f'  centre band |y|/D < 0.3 (all z): median cos {wmed(c[s], w[s]):.3f}; also z/D>0.2: ', end='')
    s2 = s & (z > 0.2)
    print(f'median {wmed(c[s2], w[s2]):.3f}, <-ln cos> {np.average(-np.log(c[s2]), weights=w[s2]):.3f}')
    print(f'  whole region: median {wmed(c[m0], w[m0]):.3f}')

print('\n--- cached 193M b010 planes: median cos for S > Smin (finite), z/D 0-2 ---')
Z = np.load('/mnt/Elements/DNS/analysis/asme2027_figs/meff_planes/ref193_b010.npz')
Dm = 1.0
for p in ('y0', 'x4', 'x10', 'x20'):
    cc, SS, kk = (Z[f'{p}_{q}'].astype(np.float64) for q in ('cos', 'S', 'k'))
    A, B = Z[f'{p}_A'] / D, Z[f'{p}_B'] / D
    if p == 'y0':
        xs, zs = A, B
        sel = (zs >= 0) & (zs <= 2)
        for lo, hi in ((3, 10), (10, 25), (3, 25)):
            s = sel & (xs >= lo) & (xs < hi) & np.isfinite(cc) & (SS > SMIN) & (kk > 0)
            pos = s & (cc > 0)
            print(f'  centre plane x/D {lo}-{hi}: n {s.sum()}, cos<=0 fraction {1 - pos.sum() / s.sum():.3f}, median cos (all S>Smin) {np.median(cc[s]):.3f}, '
                  f'median cos (cos>0) {np.median(cc[pos]):.3f}, <-ln cos>(cos>0) {np.mean(-np.log(cc[pos])):.3f}')
    else:
        s = (B <= 2) & np.isfinite(cc) & (SS > SMIN) & (kk > 0)
        pos = s & (cc > 0)
        print(f'  cross-plane {p} (|y|<3, z/D<=2): n {s.sum()}, cos<=0 {1 - pos.sum() / s.sum():.3f}, median cos {np.median(cc[pos]):.3f}, <-ln cos> {np.mean(-np.log(cc[pos])):.3f}; '
              f'centre strip |y|<0.3: median {np.median(cc[pos & (np.abs(A) < 0.3)]):.3f}')
