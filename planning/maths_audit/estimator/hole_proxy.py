#!/usr/bin/env python3
"""Proxy for a wall-distance cut in the metering section: cylinder of radius 0.5 D about the PCA axis of the region's points (x-z plane).
Indicative only, NOT the STAR wall distance. usage: hole_proxy.py pts.npz xyz.npz"""
import sys
import numpy as np
d = np.load(sys.argv[1]); X = np.load(sys.argv[2])['xyz']; names = [str(x) for x in d['regnames']]
V = d['w'] / 8
g = lambda q: d[q].astype(float)  # noqa: E731
S, k = g('c193_S'), g('c193_k')
ok = np.isfinite(S) & np.isfinite(k) & (k > 0) & np.isfinite(g('rans_mut')) & np.isfinite(g('S_r')) & (g('rans_k') > 0) & np.isfinite(g('c193_L')) & np.isfinite(g('c193_cos'))
m = ok & (d['reg'] == names.index('hole_metering'))
P = X[m][:, [0, 2]]; c = np.average(P, axis=0, weights=V[m]); u, s, vt = np.linalg.svd((P - c) * np.sqrt(V[m])[:, None], full_matrices=False)
ax = vt[0]; nrm = np.array([-ax[1], ax[0]])
print('axis angle from x [deg]', np.degrees(np.arctan2(ax[1], ax[0])))
t = (X[:, [0, 2]] - c) @ nrm; y = X[:, 1]; r = np.hypot(t, y)
print('r percentiles (volume) 50/90/99/100:', [round(float(np.percentile(r[m], q)), 3) for q in (50, 90, 99, 100)])
w = V * S ** 2
for dc in (0.1, 0.2, 0.3):
    keep = m & (0.5 - r >= dc)
    print(f'd >= {dc} D (proxy): volume kept {V[keep].sum() / V[m].sum():.3f}, 193M strain-weight kept {w[keep].sum() / w[m].sum():.2e}')
