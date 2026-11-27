#!/usr/bin/env python3
"""Strain-weight share of strict counter-gradient points (cos < 0) against wall nodes (k_193 = 0, cos set to 0). usage: cos_neg_weight.py pts.npz"""
import sys
import numpy as np
d = np.load(sys.argv[1]); n = [str(x) for x in d['regnames']]
V = d['w'] / 8; S = d['c193_S'].astype(float); c = d['c193_cos'].astype(float); k = d['c193_k']
keys = ['c193_rho', 'c193_L', 'rans_mut', 'S_r', 'rans_L', 'rans_k', 'rans_rho']
fin = np.isfinite(S) & np.isfinite(c) & (k >= 0) & (d['rans_k'] > 0)
for q in keys:
    fin &= np.isfinite(d[q])
for r in ['near_jet_0_3', 'wake_3_10', 'wake_10_25', 'exit', 'hole_diffuser', 'hole_metering', 'upstream_bl']:
    m = fin & (d['reg'] == n.index(r)); w = V * S ** 2
    mu = d['c193_rho'] * np.sqrt(np.maximum(k, 0)) * d['c193_L'] * c
    print(f'{r:14s} cos<0 weight share {w[m & (c < 0)].sum() / w[m].sum():.3f} (volume {V[m & (c < 0)].sum() / V[m].sum():.3f}) | k=0 weight share {w[m & (k == 0)].sum() / w[m].sum():.3f} | '
          f'signed share of sum(w mu) from cos<0: {(w * mu)[m & (c < 0)].sum() / (w * mu)[m].sum():+.3f}')
