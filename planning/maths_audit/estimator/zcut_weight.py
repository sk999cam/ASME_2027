#!/usr/bin/env python3
"""Share of strain weight removed by z/D < 0.3 on plate regions (wall nodes excluded), and counter-gradient share after the cut.
usage: zcut_weight.py pts.npz xyz.npz > zcut_weight.log"""
import sys
import numpy as np
d = np.load(sys.argv[1]); z = np.load(sys.argv[2])['xyz'][:, 2]; n = [str(x) for x in d['regnames']]
V = d['w'] / 8; S = d['c193_S'].astype(float); k = d['c193_k']; c = d['c193_cos'].astype(float); SR = d['S_r'].astype(float)
ok = np.isfinite(S) & (k > 0) & np.isfinite(d['rans_mut']) & (d['rans_k'] > 0) & np.isfinite(SR) & np.isfinite(d['c193_L']) & np.isfinite(c) & np.isfinite(d['c193_rho'])
mu = d['c193_rho'] * np.sqrt(np.maximum(k, 0)) * d['c193_L'] * c
PL = ['near_jet_0_3', 'wake_3_10', 'wake_10_25', 'upstream_bl']
print('plate regions, wall nodes excluded: share removed by z/D<0.3 (volume, 193M strain weight, RANS strain weight)')
for r in PL:
    m = ok & (d['reg'] == n.index(r)); cc = m & (z < 0.3)
    print(f'{r:13s} vol {V[cc].sum() / V[m].sum():.3f}  w193 {(V * S ** 2)[cc].sum() / (V * S ** 2)[m].sum():.3f}  wRANS {(V * SR ** 2)[cc].sum() / (V * SR ** 2)[m].sum():.3f}')
print('after wall exclusion + z cut: strict cos<0 share of 193M strain weight, and signed share of sum(w mu)')
for r in PL:
    m = ok & (d['reg'] == n.index(r)) & (z >= 0.3); w = V * S ** 2
    print(f'{r:13s} w(cos<0) {w[m & (c < 0)].sum() / w[m].sum():.3f}  signed share of sum(w mu) {(w * mu)[m & (c < 0)].sum() / (w * mu)[m].sum():+.3f}')
