"""Read-only recomputation of the near-wall eddy-viscosity split ranges quoted in main_claude.tex.

Imports the figure module's compute_curve (no files written). Prints the ranges over the slab centres 4.5-19.5.
"""
import os
import sys

import numpy as np

CODE = '/home/sk999/Documents/ASME_2027/figures/code'
sys.path.insert(0, CODE)
import _meff_split_nearfar as SP  # noqa: E402

zt = np.loadtxt(os.path.join(CODE, 'data', 'near_wall_boundary_zx.txt'))
ZXT = (zt[:, 0], zt[:, 1])
CEN = np.arange(4.5, 19.5 + 1e-9, 0.25)
cur = SP.compute_curve(ZXT, CEN, 1.0)
mu = np.exp([r['tot'] for r in cur])
x = np.array([r['x'] for r in cur])
k = np.array([r['vmean']['k'][0] / r['vmean']['k'][1] for r in cur])
L = np.array([r['vmean']['L'][0] / r['vmean']['L'][1] for r in cur])
ic = np.array([1 / r['vmean']['c'][1] for r in cur])
print('near: mu ratio min %.3f (x %.2f) max %.3f (x %.2f)' % (mu.min(), x[mu.argmin()], mu.max(), x[mu.argmax()]))
print('near: L ratio min %.3f max %.3f' % (L.min(), L.max()))
print('near: k^1/2 ratio min %.3f max %.3f' % (k.min(), k.max()))
print('near: 1/cos193 first %.3f (x %.2f) last %.3f (x %.2f) min %.3f max %.3f' % (ic[0], x[0], ic[-1], x[-1], ic.min(), ic.max()))
fm = SP.compute_curve(ZXT, [12.0], 16.0)[0]
print('near: layer mean 4-20 mu ratio %.3f' % np.exp(fm['tot']))
cura = SP.compute_curve(ZXT, CEN, 1.0, layer='above')
mua = np.exp([r['tot'] for r in cura])
La = np.array([r['vmean']['L'][0] / r['vmean']['L'][1] for r in cura])
ka = np.array([r['vmean']['k'][0] / r['vmean']['k'][1] for r in cura])
ica = np.array([1 / r['vmean']['c'][1] for r in cura])
kept = np.array([r['kept'] for r in cura])
print('above: mu %.3f-%.3f  L %.3f-%.3f  k %.3f-%.3f  1/cos %.3f-%.3f  kept %.3f-%.3f' % (
    mua.min(), mua.max(), La.min(), La.max(), ka.min(), ka.max(), ica.min(), ica.max(), kept.min(), kept.max()))
fa = SP.compute_curve(ZXT, [12.0], 16.0, layer='above')[0]
va = fa['vmean']
print('above: layer mean mu %.3f k %.3f L %.3f 1/cos %.3f kept %.3f' % (
    np.exp(fa['tot']), va['k'][0] / va['k'][1], va['L'][0] / va['L'][1], 1 / va['c'][1], fa['kept']))
